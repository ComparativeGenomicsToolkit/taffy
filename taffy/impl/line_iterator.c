#include "line_iterator.h"
#include "sonLib.h"

#ifdef USE_HTSLIB
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#endif

// Process-wide bgzf thread count, set once early in each command's main()
// via LI_set_bgzf_threads().  Default of 1 means "no extra threads" and
// preserves pre-existing behaviour.  Only consulted at LI/LW open time;
// changes after a handle is opened do not affect that handle.
static int bgzf_threads = 1;

void LI_set_bgzf_threads(int n) {
    bgzf_threads = (n > 0) ? n : 1;
}

#ifdef USE_HTSLIB
// Enable bgzf_mt on `bgzf` if threads > 1 AND the stream is actual BGZF
// (compression == 2).  Plain-gzip and uncompressed streams have no block
// structure that bgzf_mt can parallelise, so we leave them alone.  Call
// AFTER bgzf_index_build_init: the htslib docs note the index must be set
// up before workers start touching the stream.
static void maybe_enable_bgzf_threads(BGZF *bgzf, const char *what) {
    if (bgzf_threads > 1 && bgzf_compression(bgzf) == 2) {
        // n_sub_blks is unused per current htslib; pass 0.  Return code <0
        // means htslib could not spawn the pool; we soldier on single-
        // threaded rather than fail the whole command.
        int rc = bgzf_mt(bgzf, bgzf_threads, 0);
        st_logInfo("bgzf_mt(%s, threads=%d) -> %d\n", what, bgzf_threads, rc);
    }
}
#endif

LI *LI_construct(FILE *fh) {
    LI *li = st_calloc(1, sizeof(LI));
#ifdef USE_HTSLIB
    li->bgzf = bgzf_dopen(fileno(fh), "r");
    assert(li->bgzf != NULL);
    if (bgzf_compression(li->bgzf) == 2) {
        if (bgzf_index_build_init(li->bgzf) != 0) {
            assert(false);
        }
    }
    maybe_enable_bgzf_threads(li->bgzf, "read");
    kstring_t ks = KS_INITIALIZE;
    li->prev_pos = bgzf_tell(li->bgzf);
    li->pos = li->prev_pos;
    bgzf_getline(li->bgzf, '\n', &ks);
    li->line = ks_release(&ks);
#else
    li->fh = fh;
    li->prev_pos = ftell(li->fh);
    li->pos = li->prev_pos;
    li->line = stFile_getLineFromFile(fh);
#endif
    return li;
}

#ifdef USE_HTSLIB
LI *LI_construct_from_path(const char *path) {
    BGZF *bgzf = bgzf_open(path, "r");
    if (bgzf == NULL) {
        fprintf(stderr, "Unable to open input %s (htslib bgzf_open failed)\n", path);
        if (strstr(path, "://") != NULL) {
            fprintf(stderr, "  URL inputs require htslib built with libcurl support.\n"
                            "  If you built htslib yourself, rerun ./configure --enable-libcurl and rebuild.\n");
        }
        return NULL;
    }
    LI *li = st_calloc(1, sizeof(LI));
    li->bgzf = bgzf;
    if (bgzf_compression(li->bgzf) == 2) {
        if (bgzf_index_build_init(li->bgzf) != 0) {
            assert(false);
        }
    }
    kstring_t ks = KS_INITIALIZE;
    li->prev_pos = bgzf_tell(li->bgzf);
    li->pos = li->prev_pos;
    bgzf_getline(li->bgzf, '\n', &ks);
    li->line = ks_release(&ks);
    return li;
}
#else
LI *LI_construct_from_path(const char *path) {
    FILE *fh = fopen(path, "r");
    if (fh == NULL) {
        fprintf(stderr, "Unable to open input %s\n", path);
        return NULL;
    }
    return LI_construct(fh);
}
#endif

void LI_destruct(LI *li) {
#ifdef USE_HTSLIB
    bgzf_close(li->bgzf);
#endif
    free(li);
}

bool LI_indexable(LI *li) {
    assert(li != NULL);
#ifdef USE_HTSLIB
    int bc = bgzf_compression(li->bgzf);
    return bc == 0 || bc == 2;
#else
    return true;
#endif
}

char *LI_get_next_line(LI *li) {
    char *l = li->line;
    li->prev_pos = li->pos;
#ifdef USE_HTSLIB
    kstring_t ks = KS_INITIALIZE;
    li->pos = bgzf_tell(li->bgzf);
    bgzf_getline(li->bgzf, '\n', &ks);
    li->line = ks_release(&ks);
#else
    li->pos = ftell(li->fh);
    li->line = stFile_getLineFromFile(li->fh);
#endif
    return l;
}

char *LI_peek_at_next_line(LI *li) {
    return li->line;
}

void LI_seek(LI *li, int64_t position) {
    li->prev_pos = position;
    li->pos = position;
#ifdef USE_HTSLIB
    int ret = bgzf_seek(li->bgzf, position, SEEK_SET);
#else
    int ret = fseek(li->fh, position, SEEK_SET);
#endif
    assert(ret == 0);
}

int64_t LI_tell(LI *li) {
    return li->prev_pos;
}

int64_t LI_get_position(LI *li) {
    return li->pos;
}

// Initial / target capacity for the LW_put* coalescing buffer.  64 KB
// matches the BGZF block size, so each flush hands the bgzf writer roughly
// one block at a time -- amortises the per-bgzf_write locking cost across
// many emitted fields.
#define LW_BUF_TARGET (64u * 1024u)

LW *LW_construct(FILE *fh, bool use_compression) {
    LW *lw = st_calloc(1, sizeof(LW));
    lw->fh = fh;
    lw->buf_cap = LW_BUF_TARGET;
    lw->buf = st_malloc(lw->buf_cap);
    lw->buf_pos = 0;
    if(use_compression) {
#ifdef USE_HTSLIB
        lw->bgzf = bgzf_dopen(fileno(fh), "w");
        assert(lw->bgzf);
        assert(bgzf_compression(lw->bgzf) == 2);
        if (bgzf_index_build_init(lw->bgzf) != 0) {
            assert(false);
        }
        maybe_enable_bgzf_threads(lw->bgzf, "write");
#endif
    }
    return lw;
}

// Send the coalesced buffer to the underlying stream, then reset.  Callers
// that emit via the legacy LW_write OR that read back via fileno() must
// invoke this first so on-disk ordering matches call order.
void LW_flush(LW *lw) {
    if (lw->buf_pos == 0) return;
#ifdef USE_HTSLIB
    if (lw->bgzf) {
        ssize_t n = bgzf_write(lw->bgzf, lw->buf, lw->buf_pos);
        assert(n == (ssize_t)lw->buf_pos);
    } else
#endif
    {
        size_t n = fwrite(lw->buf, 1, lw->buf_pos, lw->fh);
        assert(n == lw->buf_pos);
    }
    lw->buf_pos = 0;
}

// Reserve `need` more bytes in the buffer, flushing if doing so would not
// fit and growing the buffer if a single field doesn't fit on its own.
// The growable path only triggers for absurdly long fields (eg. a single
// sequence name longer than 64 KB) -- the steady state is "flush at full".
static void lw_reserve(LW *lw, size_t need) {
    if (lw->buf_pos + need <= lw->buf_cap) return;
    LW_flush(lw);
    if (need > lw->buf_cap) {
        while (lw->buf_cap < need) lw->buf_cap *= 2;
        lw->buf = st_realloc(lw->buf, lw->buf_cap);
    }
}

// Out-of-line slow path for LW_putc.  The inline fast path in the header
// only branches here when the buffer is full; this flushes (resetting
// buf_pos to 0) and stores the byte.
void LW_putc_slow(LW *lw, char c) {
    LW_flush(lw);
    lw->buf[lw->buf_pos++] = c;
}

void LW_putn(LW *lw, const char *s, size_t n) {
    lw_reserve(lw, n);
    memcpy(lw->buf + lw->buf_pos, s, n);
    lw->buf_pos += n;
}

void LW_puts(LW *lw, const char *s) {
    LW_putn(lw, s, strlen(s));
}

void LW_putrep(LW *lw, char c, size_t n) {
    while (n > 0) {
        if (lw->buf_pos == lw->buf_cap) LW_flush(lw);
        size_t room = lw->buf_cap - lw->buf_pos;
        size_t k = (n < room) ? n : room;
        memset(lw->buf + lw->buf_pos, (unsigned char)c, k);
        lw->buf_pos += k;
        n -= k;
    }
}

// Two-digit-per-iter itoa.  ~5x faster than vsnprintf("%" PRIi64) because
// it avoids re-parsing the format string and only does one division per
// pair of digits.  Worst case (LLONG_MIN) needs 20 chars + sign.
static const char lw_itoa_pairs[200] =
    "00010203040506070809"
    "10111213141516171819"
    "20212223242526272829"
    "30313233343536373839"
    "40414243444546474849"
    "50515253545556575859"
    "60616263646566676869"
    "70717273747576777879"
    "80818283848586878889"
    "90919293949596979899";

void LW_puti64(LW *lw, int64_t v) {
    char tmp[24];
    int p = (int)sizeof(tmp);
    uint64_t u;
    bool neg = false;
    if (v < 0) { neg = true; u = (uint64_t)(-(v + 1)) + 1; }  // safe LLONG_MIN
    else u = (uint64_t)v;
    while (u >= 100) {
        uint64_t q = u / 100;
        unsigned r = (unsigned)(u - q * 100) * 2;
        tmp[--p] = lw_itoa_pairs[r + 1];
        tmp[--p] = lw_itoa_pairs[r];
        u = q;
    }
    if (u >= 10) {
        unsigned r = (unsigned)u * 2;
        tmp[--p] = lw_itoa_pairs[r + 1];
        tmp[--p] = lw_itoa_pairs[r];
    } else {
        tmp[--p] = (char)('0' + u);
    }
    if (neg) tmp[--p] = '-';
    LW_putn(lw, tmp + p, sizeof(tmp) - p);
}

void LW_destruct(LW *lw, bool clean_up_file_handle) {
    LW_flush(lw);
    free(lw->buf);
    lw->buf = NULL;
#ifdef USE_HTSLIB
    if(lw->bgzf) {
        if(bgzf_flush(lw->bgzf)) {
            assert(0); // Flush failed
        }
        if(bgzf_close(lw->bgzf)) {
            assert(0); // Close failed
        }
    }
#endif
    if(clean_up_file_handle) {
        fclose(lw->fh);
    }
    free(lw);
}

int LW_write(LW *lw, const char *string, ...) {
    LW_flush(lw); // preserve ordering against any pending LW_put* bytes
#ifdef USE_HTSLIB
    int i;
    va_list ap;
    if(lw->bgzf) { // Use bgzf compression
        // Figure out how long the string is to build the buffer
        va_start(ap, string);
        int j = vsnprintf(NULL, 0, string, ap)+1;
        assert(j >= 0);
        va_end(ap);
        // Now write the string into the buffer
        va_start(ap, string);
        char ret[j];
        i = vsnprintf(ret, j, string, ap);
        assert(ret[i] == '\0');
        // Finally, write the buffer to the bgzf stream
        i = bgzf_write(lw->bgzf, ret, i);
        assert(i+1 == j);
        va_end(ap);
    }
    else { // No compression, just fprintf to the stream
        va_start(ap, string);
        i = vfprintf(lw->fh, string, ap);
        va_end(ap);
    }
    return i;
#else
    va_list ap;
    va_start(ap, string);
    int i = vfprintf(lw->fh, string, ap);
    va_end(ap);
    return i;
#endif
}


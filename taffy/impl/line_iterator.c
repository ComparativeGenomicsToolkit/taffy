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

LW *LW_construct(FILE *fh, bool use_compression) {
    LW *lw = st_calloc(1, sizeof(LW));
    lw->fh = fh;
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

void LW_destruct(LW *lw, bool clean_up_file_handle) {
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


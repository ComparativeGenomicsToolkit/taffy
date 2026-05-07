#include "line_iterator.h"
#include "sonLib.h"

#ifdef USE_HTSLIB
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
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
#ifdef USE_HTSLIB
    lw->buf_cap = 65536;
    lw->buf = st_malloc(lw->buf_cap);
    if(use_compression) {
        lw->bgzf = bgzf_dopen(fileno(fh), "w");
        assert(lw->bgzf);
        assert(bgzf_compression(lw->bgzf) == 2);
        if (bgzf_index_build_init(lw->bgzf) != 0) {
            assert(false);
        }
    }
#endif
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
    free(lw->buf);
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
        int n;
        while (1) {
            va_start(ap, string);
            n = vsnprintf(lw->buf, lw->buf_cap, string, ap);
            va_end(ap);
            if (n >= 0 && (size_t)n < lw->buf_cap) break;
            // Buffer too small; grow and retry
            lw->buf_cap = (n >= 0) ? (size_t)(n + 1) : lw->buf_cap * 2;
            lw->buf = st_realloc(lw->buf, lw->buf_cap);
        }
        ssize_t written = bgzf_write(lw->bgzf, lw->buf, n);
        assert(written == (ssize_t)n);
        i = (int)written;
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

int LW_write_bytes(LW *lw, const char *data, int len) {
#ifdef USE_HTSLIB
    if(lw->bgzf) {
        ssize_t written = bgzf_write(lw->bgzf, data, (size_t)len);
        assert(written == (ssize_t)len);
        return (int)written;
    }
#endif
    int written = (int)fwrite(data, 1, (size_t)len, lw->fh);
    assert(written == len);
    return written;
}


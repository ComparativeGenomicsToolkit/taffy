#include "line_iterator.h"
#include "sonLib.h"

#ifdef USE_HTSLIB
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#endif

struct _LI {
#ifdef USE_HTSLIB
    BGZF *bgzf;
#else
    FILE *fh;
#endif
    char *line;
    int64_t prev_pos; // position before reading the current buffer
    int64_t pos;      // position after reading the curent buffer    
};

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

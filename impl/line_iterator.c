#include "line_iterator.h"
#include "sonLib.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"

LI *LI_construct(FILE *fh) {
    LI *li = st_calloc(1, sizeof(LI));
    li->bgzf = bgzf_dopen(fileno(fh), "r");
    assert(li->bgzf != NULL);
    if (li->bgzf->is_compressed) {
        if (bgzf_compression(li->bgzf) != 2) {
            fprintf(stderr, "Input file must be compressed with bgzip, not gzip\n");
            exit(1);
        }
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

void LI_destruct(LI *li) {
    free(li->line);
    // todo: review, as the file is currently getting closed by client
    bgzf_close(li->bgzf);
}

char *LI_get_next_line(LI *li) {
    char *l = li->line;
    kstring_t ks = KS_INITIALIZE;
    li->prev_pos = li->pos;
    li->pos = bgzf_tell(li->bgzf);
    bgzf_getline(li->bgzf, '\n', &ks);
    li->line = ks_release(&ks);
    return l;
}

char *LI_peek_at_next_line(LI *li) {
    return li->line;
}

void LI_seek(LI *li, int64_t position) {
    li->prev_pos = position;
    li->pos = position;
    int ret = bgzf_seek(li->bgzf, position, SEEK_SET);
    assert(ret == 0);
}

int64_t LI_tell(LI *li) {
    return li->prev_pos;
}

#include "line_iterator.h"
#include "sonLib.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"

LI *LI_construct(FILE *fh) {
    LI *li = st_calloc(1, sizeof(LI));
    li->bgzf = bgzf_dopen(fileno(fh), "r");
    assert(li->bgzf != NULL);
    kstring_t ks = KS_INITIALIZE;
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
    bgzf_getline(li->bgzf, '\n', &ks);
    li->line = ks_release(&ks);
    return l;
}

char *LI_peek_at_next_line(LI *li) {
    return li->line;
}


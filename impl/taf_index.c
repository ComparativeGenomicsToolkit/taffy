#include "taf.h"
#include "taf_index.h"
#include "htslib/bgzf.h"

// todo: it would be nice to be able to use something from the taf.h API for this!!
// but it's currently inside write_coordinates()
static void write_taf_tags(Alignment *alignment, FILE* fh) {
    int64_t i = 0;
    for (Alignment_Row *row = alignment->row; row; row = row->n_row, ++i) {
        fprintf(fh, " s %" PRIi64 " %s %" PRIi64 " %c %" PRIi64 "",
                i, row->sequence_name, row->start, row->strand ? '+' : '-', row->sequence_length);
    }
}

int index_taf(LI *li, FILE* idx_fp, int64_t index_block_size){
    Alignment *alignment, *p_alignment = NULL;
    char *prev_ref = NULL;
    int64_t prev_position = 0;
    int64_t prev_file_offset = bgzf_utell(li->bgzf);
    
    // read the TAF header
    stList *tags = taf_read_header(li);    
    assert(stList_length(tags) % 2 == 0);

    while((alignment = taf_read_block(p_alignment, false, li)) != NULL) {

        char *cur_ref = alignment->row->sequence_name;
        int64_t cur_offset = alignment->row->start;

        // shouldn't need to handle negative strand on reference, right?
        assert(alignment->row->strand != 0);

        // we need to update our index if we're on a new reference contig
        // or we're on the same contig but >= index_block_size bases away
        if (!prev_ref || strcmp(cur_ref, prev_ref) != 0 ||
            cur_offset - prev_position >= index_block_size) {

            fprintf(idx_fp, "%s\t%ld\t%ld\t", cur_ref, cur_offset, prev_file_offset);
            write_taf_tags(alignment, idx_fp);
            fprintf(idx_fp, "\n");
        }

        prev_ref = alignment->row->sequence_name;
        prev_position = alignment->row->start;
        prev_file_offset = bgzf_utell(li->bgzf);
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment);
        }
        p_alignment = alignment;

    }
    if(p_alignment != NULL) {
        alignment_destruct(p_alignment);
    }

    stList_destruct(tags);
    
    return 0;
}



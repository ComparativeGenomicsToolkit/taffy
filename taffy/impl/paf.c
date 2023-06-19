#include "taf.h"
#include "sonLib.h"
#include "line_iterator.h"

// write a single PAF row from the pairwise alignment between the given two MAF block rows
static void paf_write_row(Alignment_Row *q_row, Alignment_Row *t_row, int64_t num_col, LW *lw) {
    char relative_strand = q_row->strand == t_row->strand ? '+' : '-';

    // query position
    char *query_name = q_row->sequence_name;
    int64_t query_length = q_row->sequence_length;
    int64_t query_start = q_row->start;
    int64_t query_end = q_row->start + q_row->length;
    // flip it around if negative (paf coordinates are always forward)
    if (q_row->strand == false) {
        query_start = query_length - query_end;
        query_end = query_start + q_row->length;
    }
    // target position
    char *target_name = t_row->sequence_name;
    int64_t target_length = t_row->sequence_length;
    int64_t target_start = t_row->start;
    int64_t target_end = t_row->start + q_row->length;
    // flip it around if negative (paf coordinates are always forward)
    if (t_row->strand == false) {
        target_start = target_length - target_end;
        target_end = target_start + t_row->length;
    }

    int64_t num_matches = 0;
    int64_t block_length = 0;

    char *cigar_string = (char*)st_malloc(32 * num_col);
    cigar_string[0] = '\0';
    char buffer[256];
    char current_event[2] = {'.', '\0'};
    int64_t current_start = -1;
    
    for (int64_t i = 0; i < num_col; ++i) {
        char event = 'M';
        if (q_row->bases[i] == '-' && t_row->bases[i] != '-') {
            event = 'D';
        } else if (q_row->bases[i] != '-' && t_row->bases[i] == '-') {
            event = 'I';
        } else if (q_row->bases[i] == '-' && t_row->bases[i] == '-') {
            event = 'B'; // B is I and D at same time
            ++block_length;
        } else {
            ++num_matches;
        }
        ++block_length;

        if (event != current_event[0]) {
            if (current_start >= 0) {
                int64_t current_length = i - current_start;
                sprintf(buffer, "%" PRIi64 "", current_length);
                strcat(cigar_string, buffer);
                if (current_event[0] != 'B') {
                    strcat(cigar_string, current_event);
                } else {
                    block_length += current_length;
                    strcat(cigar_string, "I");
                    strcat(cigar_string, buffer);
                    strcat(cigar_string, "D");
                }
            }
            current_event[0] = event;
            current_start = i;
        }        
    }

    // finalize trailing event
    if (current_start < num_col && current_start >= 0) {
        int64_t i = num_col;
        int64_t current_length = i - current_start;
        sprintf(buffer, "%" PRIi64 "", current_length);
        strcat(cigar_string, buffer);
        if (current_event[0] != 'B') {
            strcat(cigar_string, current_event);
        } else {
            strcat(cigar_string, "I");
            strcat(cigar_string, buffer);
            strcat(cigar_string, "D");
        }        
    }

    // write it out
    LW_write(lw, "%s\t%ld\t%ld\t%ld\t%c", query_name, query_length, query_start, query_end, relative_strand);
    LW_write(lw, "\t%s\t%ld\t%ld\t%ld", target_name, target_length, target_start, target_end);
    LW_write(lw, "\t%ld\t%ld\t%ld\tcg:Z:%s\n", num_matches, block_length, 255, cigar_string);
    
    free(cigar_string);

}

void paf_write_block(Alignment *alignment, LW *lw, bool all_to_all) {
    for (Alignment_Row *t_row = alignment->row; t_row != NULL; t_row = t_row->n_row) {
        for (Alignment_Row *q_row = t_row->n_row; q_row != NULL; q_row = q_row->n_row) {            
            paf_write_row(q_row, t_row, alignment->column_number, lw);
        }        
        if (!all_to_all) {
            break;
        }
    }
}

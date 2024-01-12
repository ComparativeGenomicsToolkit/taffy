#include "taf.h"
#include "sonLib.h"
#include "line_iterator.h"

// write a single PAF row from the pairwise alignment between the given two MAF block rows
static void paf_write_row(Alignment_Row *q_row, Alignment_Row *t_row, int64_t num_col, bool cs_cigar, LW *lw) {
    char relative_strand = q_row->strand == t_row->strand ? '+' : '-';
    bool flip_cigar = t_row->strand == false;
    
    // query position
    char *query_name = q_row->sequence_name;
    int64_t query_length = q_row->sequence_length;
    int64_t query_start = q_row->start;
    int64_t query_end = q_row->start + q_row->length;
    // flip it around if negative (paf coordinates are always forward)
    if (q_row->strand == 0) {
        query_start = query_length - query_end;
        query_end = query_start + q_row->length;
    }
    // target position
    char *target_name = t_row->sequence_name;
    int64_t target_length = t_row->sequence_length;
    int64_t target_start = t_row->start;
    int64_t target_end = t_row->start + t_row->length;
    // flip it around if negative (paf coordinates are always forward)
    if (t_row->strand == 0) {
        target_start = target_length - target_end;
        target_end = target_start + t_row->length;
    }

    int64_t num_matches = 0;
    int64_t block_length = 0;

    char *cigar_string = (char*)st_malloc(32 * num_col);
    cigar_string[0] = '\0';
    char *buffer = (char*)st_malloc(32 * num_col);
    char *query_buffer = NULL;
    int64_t query_buffer_length = 0;
    char *target_buffer = NULL;
    int64_t target_buffer_length = 0;
    if (cs_cigar) {
        query_buffer = (char*)st_malloc(32 * num_col);
        target_buffer = (char*)st_malloc(32 * num_col);
    }
    char current_event[2] = {'.', '\0'};
    int64_t current_start = -1;
    int64_t current_length = 0;
    
    for (int64_t i = 0; i < num_col; ++i) {
        // scan the next column
        int64_t pos = flip_cigar ? num_col - 1 - i : i;
        char event = '.';
        if (t_row->bases[pos] != '-' && q_row->bases[pos] != '-') {
            if (t_row->bases[pos] != q_row->bases[pos] && cs_cigar) {
                event = '*';
            } else {
                event = 'M';
            }
            ++num_matches;
        } else if (t_row->bases[pos] == '-' && q_row->bases[pos] != '-') {
            event = 'I';
        } else if (t_row->bases[pos] != '-' && q_row->bases[pos] == '-') {
            event = 'D';
        } else {
            assert(t_row->bases[pos] == '-' &&  q_row->bases[pos] == '-');
        }

        // print the previous event if we need to start a new event
        if (event != '.' && (event != current_event[0] || event == '*')) {
            if (current_start >= 0) {
                if (cs_cigar) {
                    query_buffer[query_buffer_length] = '\0';
                    target_buffer[target_buffer_length] = '\0';
                    if (current_event[0] == 'M') {
                        assert(current_length == query_buffer_length && current_length == target_buffer_length);
                        sprintf(buffer, "=%s", query_buffer);
                    } else if (current_event[0] == '*') {
                        assert(current_length == 1 && query_buffer_length == 1 && target_buffer_length == 1);
                        sprintf(buffer, "*%s%s", target_buffer, query_buffer);
                    } else if (current_event[0] == 'I') {
                        assert(target_buffer_length == 0 && query_buffer_length > 0 && current_length == query_buffer_length);
                        sprintf(buffer, "+%s", query_buffer);
                    } else {
                        assert(current_event[0] == 'D');
                        assert(query_buffer_length == 0 && target_buffer_length > 0 && current_length == target_buffer_length);
                        sprintf(buffer, "-%s", target_buffer);
                    }
                    strcat(cigar_string, buffer);
                } else {
                    sprintf(buffer, "%" PRIi64 "", current_length);
                    strcat(cigar_string, buffer);
                    strcat(cigar_string, current_event);
                }
            }
            current_event[0] = event;
            current_start = i;
            current_length = 0;
            query_buffer_length = 0;
            target_buffer_length = 0;
        }

        // update the current event
        if (event != '.') {
            ++block_length;
            ++current_length;
            if (cs_cigar) {
                if (t_row->bases[pos] != '-') {
                    target_buffer[target_buffer_length++] = t_row->bases[pos];
                }
                if (q_row->bases[pos] != '-') {
                    query_buffer[query_buffer_length++] = q_row->bases[pos];
                }
            }
        }

    }

    // finalize trailing event
    if (current_start < num_col && current_start >= 0 && current_event[0] != '.') {
        int64_t i = num_col;
        // note, this code block is (must be) identical to above. todo: factor out
        if (cs_cigar) {
            query_buffer[query_buffer_length] = '\0';
            target_buffer[target_buffer_length] = '\0';
            if (current_event[0] == 'M') {
                assert(current_length == query_buffer_length && current_length == target_buffer_length);
                sprintf(buffer, "=%s", query_buffer);
            } else if (current_event[0] == '*') {
                assert(current_length == 1 && query_buffer_length == 1 && target_buffer_length == 1);
                sprintf(buffer, "*%s%s", target_buffer, query_buffer);
            } else if (current_event[0] == 'I') {
                assert(target_buffer_length == 0 && query_buffer_length > 0 && current_length == query_buffer_length);
                sprintf(buffer, "+%s", query_buffer);
            } else {
                assert(current_event[0] == 'D');
                assert(query_buffer_length == 0 && target_buffer_length > 0 && current_length == target_buffer_length);
                sprintf(buffer, "-%s", target_buffer);
            }
            strcat(cigar_string, buffer);
        } else {
            sprintf(buffer, "%" PRIi64 "", current_length);
            strcat(cigar_string, buffer);
            strcat(cigar_string, current_event);
        }
    }

    // write it out
    LW_write(lw, "%s\t%ld\t%ld\t%ld\t%c", query_name, query_length, query_start, query_end, relative_strand);
    LW_write(lw, "\t%s\t%ld\t%ld\t%ld", target_name, target_length, target_start, target_end);
    LW_write(lw, "\t%ld\t%ld\t%ld\tc%c:Z:%s\n", num_matches, block_length, 255, cs_cigar ? 's' : 'g', cigar_string);
    
    free(cigar_string);
    free(buffer);
    free(query_buffer);
    free(target_buffer);

}

void paf_write_block(Alignment *alignment, LW *lw, bool all_to_all, bool cs_cigar) {
    for (Alignment_Row *t_row = alignment->row; t_row != NULL; t_row = t_row->n_row) {
        for (Alignment_Row *q_row = t_row->n_row; q_row != NULL; q_row = q_row->n_row) {            
            paf_write_row(q_row, t_row, alignment->column_number, cs_cigar, lw);
        }        
        if (!all_to_all) {
            break;
        }
    }
}

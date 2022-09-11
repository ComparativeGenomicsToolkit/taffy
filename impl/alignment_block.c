#include "taf.h"
#include "ond.h"

void Alignment_Row_destruct(Alignment_Row *row) {
    if(row->bases != NULL) {
        free(row->bases);
    }
    if(row->sequence_name != NULL) {
        free(row->sequence_name);
    }
    if(row->left_gap_sequence != NULL) {
        free(row->left_gap_sequence);
    }
    free(row);
}

void alignment_destruct(Alignment *alignment) {
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        Alignment_Row *r = row;
        row = row->n_row;
        Alignment_Row_destruct(r);
    }
    free(alignment);
}

static stList *get_rows_in_a_list(Alignment_Row *row) {
    stList *l = stList_construct();
    while(row != NULL) {
        stList_append(l, row);
        row = row->n_row;
    }
    return l;
}

bool alignment_row_is_predecessor(Alignment_Row *left_row, Alignment_Row *right_row) {
    // Do the rows match
    return strcmp(left_row->sequence_name, right_row->sequence_name) == 0 && left_row->strand == right_row->strand &&
            left_row->start + left_row->length <= right_row->start;
}

void alignment_link_adjacent(Alignment *left_alignment, Alignment *right_alignment, bool allow_row_substitutions) {
    stList *left_rows = get_rows_in_a_list(left_alignment->row);
    stList *right_rows = get_rows_in_a_list(right_alignment->row);
    // get the alignment of the rows
    WFA *wfa = WFA_construct(left_rows, right_rows, (bool (*)(void *, void *))alignment_row_is_predecessor, 1,
                             allow_row_substitutions ? 1 : 100000000); // Use unit gap and mismatch costs for the diff
                             // unless we disallow substitutions, in which case use an arbitrarily large mismatch cost
    stList *aligned_rows = WFA_get_alignment(wfa);
    // Remove any previous links
    Alignment_Row *row = left_alignment->row;
    while(row != NULL) {
        row->r_row = NULL;
        row = row->n_row;
    }
    row = right_alignment->row;
    while(row != NULL) {
        row->l_row = NULL;
        row = row->n_row;
    }
    // connect up the rows
    assert(stList_length(aligned_rows) % 2 == 0); // must be even length
    for(int64_t i=0; i<stList_length(aligned_rows); i+=2) {
        Alignment_Row *left_row = stList_get(aligned_rows, i);
        Alignment_Row *right_row = stList_get(aligned_rows, i+1);
        left_row->r_row = right_row;
        right_row->l_row = left_row;
        if(!allow_row_substitutions) {
            assert(alignment_row_is_predecessor(left_row, right_row));
        }
    }
    // clean up
    stList_destruct(left_rows);
    stList_destruct(right_rows);
    stList_destruct(aligned_rows);
    WFA_destruct(wfa);
}

int64_t alignment_length(Alignment *alignment) {
    return alignment->row == NULL ? 0 : strlen(alignment->row->bases);
}

int64_t alignment_total_gap_length(Alignment *left_alignment, bool align_gap_sequences) {
    Alignment_Row *l_row = left_alignment->row;
    int64_t total_interstitial_gap_length = 0;
    while(l_row != NULL) {
        if(l_row->r_row != NULL && alignment_row_is_predecessor(l_row, l_row->r_row)) {
            int64_t i = l_row->r_row->start - (l_row->start + l_row->length);
            if(align_gap_sequences) {
                if (i > total_interstitial_gap_length) {
                    total_interstitial_gap_length = i;
                }
            }
            else {
                total_interstitial_gap_length += i;
            }
        }
        l_row = l_row->n_row; // Move to the next left alignment row
    }
    return total_interstitial_gap_length;
}

stList *parse_header(stList *tokens, char *header_prefix, char *delimiter) {
    stList *tags = stList_construct3(0, free);
    if(stList_length(tokens) == 0 || strcmp((char *)stList_get(tokens, 0), header_prefix) != 0) {
        st_errAbort("Header line does not start with %s\n", header_prefix);
    }
    for(int64_t i=1; i<stList_length(tokens); i++) {
        char *tag = stList_get(tokens, i);
        stList *tag_tokens = stString_splitByString(tag, delimiter);
        if(stList_length(tag_tokens) != 2) {
            st_errAbort("Header line tags not separated by ':' character: %s\n", tag);
        }
        stList_append(tags, stString_copy(stList_get(tag_tokens, 0)));
        stList_append(tags, stString_copy(stList_get(tag_tokens, 1)));
    }
    return tags;
}

int64_t alignment_number_of_common_rows(Alignment *left_alignment, Alignment *right_alignment) {
    // First un-link any rows that are substitutions as these can't be merged
    Alignment_Row *r_row = right_alignment->row;
    int64_t shared_rows = 0;
    while (r_row != NULL) {
        if (r_row->l_row != NULL && alignment_row_is_predecessor(r_row->l_row, r_row)) {
            shared_rows++;
        }
        r_row = r_row->n_row;
    }
    return shared_rows;
}

static char *make_gap(int64_t length) {
    char gap_alignment[length+1];
    for(int64_t i=0; i<length; i++) {
        gap_alignment[i] = '-';
    }
    gap_alignment[length] = '\0';
    return stString_copy(gap_alignment);
}

Alignment *alignment_merge_adjacent(Alignment *left_alignment, Alignment *right_alignment, bool align_gap_sequences) {
    // First un-link any rows that are substitutions as these can't be merged
    Alignment_Row *r_row = right_alignment->row;
    while(r_row != NULL) {
        if(r_row->l_row != NULL && !alignment_row_is_predecessor(r_row->l_row, r_row)) {
            assert(r_row->l_row->r_row == r_row);
            r_row->l_row->r_row = NULL; // unlink 1
            r_row->l_row = NULL; // unlink 2
        }
        r_row = r_row->n_row;
    }

    // Get the length of the left and right alignments
    int64_t left_alignment_length = alignment_length(left_alignment);
    int64_t right_alignment_length = alignment_length(right_alignment);

    // Add the new rows in the right alignment to the left alignment
    r_row = right_alignment->row; Alignment_Row **p_l_row = &(left_alignment->row);
    while(r_row != NULL) {
        if(r_row->l_row == NULL) { // Is an insertion
            // Make a new l_row
            Alignment_Row *l_row = st_calloc(1, sizeof(Alignment_Row));

            // Set coordinates
            l_row->sequence_name = stString_copy(r_row->sequence_name);
            l_row->start = r_row->start;
            l_row->length = 0; // is an empty alignment
            l_row->sequence_length = r_row->sequence_length;
            l_row->strand = r_row->strand;
            l_row->bases = make_gap(left_alignment_length);

            // Connect the left and right rows
            l_row->r_row = r_row;
            r_row->l_row = l_row;

            // Now connect the left row into the left alignment at the correct place
            l_row->n_row = *p_l_row;
            *p_l_row = l_row;

            // Update p_l_row to point at l_row's n_row pointer
            p_l_row = &(l_row->n_row);

            // Increase the row number
            left_alignment->row_number++;
        }
        else {
            // Update p_l_row to point at r_row->l_row's n_row pointer
            p_l_row = &(r_row->l_row->n_row);
        }
        r_row = r_row->n_row; // Move to the next right alignment row
    }

    // Calculate the sum of the length of any interstitial inserts
    int64_t total_interstitial_gap_length = alignment_total_gap_length(left_alignment, align_gap_sequences);

    // Now finally extend the left alignment rows to include the right alignment rows
    int64_t i=0; // An index into the interstitial alignment coordinate
    Alignment_Row *l_row = left_alignment->row;
    char *right_gap = make_gap(right_alignment_length + total_interstitial_gap_length); // any trailing bases needed
    while(l_row != NULL) {
        if(l_row->r_row == NULL) {
            // Is a deletion, so add in trailing gaps equal in length to the right alignment length plus any interstitial
            // gap
            char *bases = stString_print("%s%s", l_row->bases, right_gap);
            free(l_row->bases);
            l_row->bases = bases;
        }
        else {
            // Check the rows agree coordinate wise
            assert(strcmp(l_row->sequence_name, l_row->r_row->sequence_name) == 0);
            assert(l_row->strand == l_row->r_row->strand);
            assert(l_row->start + l_row->length <= l_row->r_row->start);
            if(l_row->r_row->left_gap_sequence != NULL) { // if the left gap sequence is defined, it's length must bridge the gap
                assert(l_row->r_row->start - (l_row->start + l_row->length) == strlen(l_row->r_row->left_gap_sequence));
            }

            // Make any interstitial gap
            char *gap_string = make_gap(total_interstitial_gap_length);
            int64_t interstitial_bases = l_row->r_row->start - (l_row->start + l_row->length);
            for(int64_t j=0; j<interstitial_bases; j++) {
                if(align_gap_sequences) {
                    gap_string[j] = l_row->r_row->left_gap_sequence != NULL ? l_row->r_row->left_gap_sequence[j] : 'N';
                }
                else {
                    assert(j + i < total_interstitial_gap_length);
                    gap_string[j + i] = l_row->r_row->left_gap_sequence != NULL ? l_row->r_row->left_gap_sequence[j] : 'N';
                }
            }
            i += interstitial_bases; // update i

            // Is not a deletion, so merge together two adjacent rows
            char *bases = stString_print("%s%s%s", l_row->bases, gap_string, l_row->r_row->bases);
            free(l_row->bases); // clean up
            free(gap_string);
            l_row->bases = bases;

            // Update the left row's length coordinate
            l_row->length += interstitial_bases + l_row->r_row->length;
            // Update the l_row's r_row pointer...
            if(l_row->r_row->r_row != NULL) { // Check
                assert(l_row->r_row->r_row->l_row == l_row->r_row);
            }
            l_row->r_row = l_row->r_row->r_row;
            if(l_row->r_row != NULL) {
                l_row->r_row->l_row = l_row;
            }
        }

        l_row = l_row->n_row; // Move to the next left alignment row
    }

    // Clean up
    alignment_destruct(right_alignment);  // Delete the right alignment
    free(right_gap);

    return left_alignment;
}


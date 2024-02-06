#include "taf.h"
#include "ond.h"
#include "sonLib.h"

/*
 * Method to merge together two adjacent alignments.
 */

static char *make_run(int64_t length, char c) {
    char gap_alignment[length+1];
    for(int64_t i=0; i<length; i++) {
        gap_alignment[i] = c;
    }
    gap_alignment[length] = '\0';
    return stString_copy(gap_alignment);
}

int64_t make_msa(int64_t string_no, int64_t column_no, int64_t max_alignment_length,
                 int64_t **msa, char **strings, int64_t *string_lengths,
                 char **msa_strings) {
    /*
     * Convert the list of aligned columns in msa into a 2D alignment of the strings.
     */
    // Build the MSA
    int64_t offsets[string_no];
    for (int64_t i=0; i<string_no; i++) { // Initialize the offsets
        offsets[i] = -1;
    }
    int64_t alignment_offset=0;
    for(int64_t j=0; j<column_no; j++) {
        // Work out the max indel length before position j
        int64_t max_indel=0;
        for(int64_t i=0; i<string_no; i++) {
            int64_t k=msa[i][j]; // Position in sequence i aligned to column j of longest sequence
            if(k != -1) { // Is not a gap, work out length
                max_indel = k - offsets[i] - 1 > max_indel ? k - offsets[i] - 1 : max_indel;
            }
        }

        // Now fill in the indels and column j
        for(int64_t i=0; i<string_no; i++) {
            int64_t k=msa[i][j]; // Position in sequence i aligned to column j of longest sequence
            if(k != -1) { // Is not a gap
                // Fill in the unaligned sequence in sequence i
                int64_t l=0;
                while(++offsets[i] < k) {
                    msa_strings[i][alignment_offset+l++] = strings[i][offsets[i]];
                }

                // Add trailing gaps
                while(l<max_indel) {
                    msa_strings[i][alignment_offset+l++] = '-';
                }

                // Add aligned position
                msa_strings[i][alignment_offset+max_indel] = strings[i][offsets[i]];
            }
            else { // Add gaps to alignment to account for max indel and column j
                for(int64_t l=0; l<=max_indel; l++) {
                    msa_strings[i][alignment_offset+l] = '-';
                }
            }
        }
        // Update the length of the alignment built so far
        alignment_offset += max_indel + 1; // indel length and the ref position
    }
    // Now add in any suffix gaps
    int64_t max_indel=0;
    for(int64_t i=0; i<string_no; i++) {
        int64_t j = string_lengths[i] - offsets[i] - 1;
        max_indel = j > max_indel ? j : max_indel;
    }
    for(int64_t i=0; i<string_no; i++) {
        int64_t j=0;
        while(++offsets[i] < string_lengths[i]) {
            msa_strings[i][alignment_offset+j++] = strings[i][offsets[i]];
        }
        while(j < max_indel) {
            msa_strings[i][alignment_offset+j++] = '-';
        }
    }
    alignment_offset += max_indel;
    return alignment_offset;
}

bool cmp_chars(void *a, void *b) {
    return ((char *)a)[0] == ((char *)b)[0];
}

int64_t align_interstitial_gaps(Alignment *alignment) {
    /*
     * Align the sequences that lie within the gaps between two adjacent blocks.
     * Return the length of the interstitial alignment.
     */

    // Add any missing gap strings in
    Alignment_Row *row = alignment->row;
    while (row != NULL) {
        if(row->l_row != NULL && alignment_row_is_predecessor(row->l_row, row) && row->left_gap_sequence == NULL) {
            row->left_gap_sequence = make_run(row->start - (row->l_row->start + row->l_row->length), 'N');
        }
        row = row->n_row;
    }

    // Find the longest gap sequence and number of sequences to align
    //TODO: Consider not allowing picking the longest sequence if it is all Ns
    row = alignment->row;
    char *longest_string = NULL;
    int64_t string_no=0, longest_string_length=0, total_string_length=0;
    while (row != NULL) {
        if(row->left_gap_sequence != NULL) {
            string_no++; // Increase the number of strings we are aligning
            int64_t i=strlen(row->left_gap_sequence);
            total_string_length += i;
            if(i > longest_string_length) {
                longest_string = row->left_gap_sequence;
                longest_string_length = strlen(longest_string);
            }
        }
        row = row->n_row;
    }

    // Align each sequence to the longest sequence
    
    // A longest sequence x sequence number sized integer matrix,
    int64_t **msa = st_calloc(string_no, sizeof(int64_t*));
    for (int64_t i = 0; i < string_no; ++i) {
        msa[i] = st_calloc(longest_string_length, sizeof(int64_t));
    }
    // each entry is the index of the position aligned at that node (or -1 if unaligned)
    char **row_strings = st_calloc(string_no, sizeof(char*)); // The gap strings
    int64_t *row_string_lengths = st_calloc(string_no, sizeof(int64_t));
    row = alignment->row;
    int64_t i=0;
    while (row != NULL) {
        if(row->left_gap_sequence != NULL) {
            row_strings[i] = row->left_gap_sequence;
            row_string_lengths[i] = strlen(row_strings[i]);
            // TODO: Consider making WFA have affine gaps
            WFA *wfa = WFA_construct(longest_string, row_strings[i], longest_string_length, row_string_lengths[i],
                                     sizeof(char), (bool (*)(void *, void *))cmp_chars, 1, 1);
            WFA_get_alignment(wfa, msa[i]); i++;
            WFA_destruct(wfa);
        }
        row = row->n_row;
    }

    // Now convert to a traditional MSA
    int64_t max_alignment_length = total_string_length < (longest_string_length+1)*longest_string_length ? total_string_length : (longest_string_length+1)*longest_string_length;
    // can not be longer than the minimum of the total length of the strings (where all bases
    char **msa_strings = st_calloc(string_no, sizeof(char*));
    for (int64_t i = 0; i < string_no; ++i) {
        msa_strings[i] = st_calloc(max_alignment_length, sizeof(char*));
    }
    // are unique gaps) and the square of the longest sequence length
    int64_t msa_length = make_msa(string_no, longest_string_length, max_alignment_length, msa, row_strings,
                                  row_string_lengths, msa_strings);
    assert(msa_length <= max_alignment_length);

    // Now copy the alignment strings back to the sequences
    row = alignment->row; i=0;
    while (row != NULL) {
        if(row->left_gap_sequence != NULL) {
            // Do a quick check that the sequence is as expected
            int64_t j=0, l=0;
            while(l < msa_length) {
                assert(j <= row_string_lengths[i]);
                if(msa_strings[i][l] != '-') {
                    assert(row->left_gap_sequence[j] == msa_strings[i][l]);
                    j++;
                }
                l++;
            }
            assert(j == row_string_lengths[i]);
            //
            free(row->left_gap_sequence);
            row->left_gap_sequence = stString_getSubString(msa_strings[i++], 0, msa_length);
        }
        row = row->n_row;
    }
    for (int64_t i = 0; i < string_no; ++i) {
        free(msa[i]);
    }
    free(msa);    
    free(row_strings);
    free(row_string_lengths);
    for (int64_t i = 0; i < string_no; ++i) {
        free(msa_strings[i]);
    }    
    free(msa_strings);
    return msa_length;
}

Alignment *alignment_merge_adjacent(Alignment *left_alignment, Alignment *right_alignment) {
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
            l_row->bases = make_run(left_alignment->column_number, '-');

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

    // Align the interstitial insert sequences, padding the left_gap_sequence strings with gaps to represent the alignment
    int64_t interstitial_alignment_length = align_interstitial_gaps(right_alignment);

    // Now finally extend the left alignment rows to include the right alignment rows
    Alignment_Row *l_row = left_alignment->row;
    char *right_gap = make_run(right_alignment->column_number + interstitial_alignment_length, '-'); // any trailing bases needed
    while(l_row != NULL) {
        if(l_row->r_row == NULL) {
            // Is a deletion, so add in trailing gaps equal in length to the right alignment length plus any interstitial
            // gap
            char *bases = stString_print("%s%s", l_row->bases, right_gap);
            free(l_row->bases);
            l_row->bases = bases;
        }
        else {
            Alignment_Row *r_row = l_row->r_row;

            // Check the rows agree coordinate wise
            assert(strcmp(l_row->sequence_name, r_row->sequence_name) == 0);
            assert(l_row->strand == r_row->strand);
            assert(l_row->start + l_row->length <= r_row->start);

            // Is not a deletion, so merge together two adjacent rows
            assert(r_row->left_gap_sequence != NULL);
            assert(strlen(r_row->left_gap_sequence) == interstitial_alignment_length);
            char *bases = stString_print("%s%s%s", l_row->bases, r_row->left_gap_sequence, r_row->bases);
            free(l_row->bases); // clean up
            l_row->bases = bases;

            // Update the left row's length coordinate
            int64_t interstitial_bases = r_row->start - (l_row->start + l_row->length);
            l_row->length += interstitial_bases + r_row->length;
            // Update the l_row's r_row pointer...
            if(r_row->r_row != NULL) { // Check pointers are correct
                assert(r_row->r_row->l_row == r_row);
            }
            l_row->r_row = r_row->r_row;
            if(l_row->r_row != NULL) {
                l_row->r_row->l_row = l_row;
            }
            // Null r_row's left and right pointers
            r_row->l_row = NULL;
            r_row->r_row = NULL;
        }

        l_row = l_row->n_row; // Move to the next left alignment row
    }

    // Calculate the number of columns in the merged alignment
    int64_t total_column_number = left_alignment->column_number + right_alignment->column_number + interstitial_alignment_length;

    // Fix the tags
    if(left_alignment->column_tags != NULL) {
        assert(right_alignment->column_tags != NULL);
        Tag **combined_column_tags = st_malloc(sizeof(Tag *) * total_column_number); // Allocate an expanded set of columns
        int64_t j=0;
        for(int64_t i=0; i<left_alignment->column_number; i++) { // Add the left alignment's column's tags
            combined_column_tags[j++] = left_alignment->column_tags[i];
        }
        for(int64_t i=0; i<interstitial_alignment_length; i++) { // Add empty tag lists for new columns
            combined_column_tags[j++] = NULL;
        }
        for(int64_t i=0; i<right_alignment->column_number; i++) { // Add the right alignment's column's tags
            combined_column_tags[j++] = right_alignment->column_tags[i];
        }
        free(left_alignment->column_tags); // Cleanup, but not the tag strings which we copied
        left_alignment->column_tags = combined_column_tags;
    }

    // Fix column number
    left_alignment->column_number = total_column_number;

    // Clean up
    alignment_destruct(right_alignment, 1);  // Delete the right alignment
    free(right_gap);

    return left_alignment;
}

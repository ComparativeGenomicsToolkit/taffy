#include "taf.h"
#include "ond.h"
#include "sonLib.h"

void tag_destruct(Tag *tag) {
    while(tag != NULL) {
        Tag *p_tag = tag;
        tag = tag->n_tag;
        free(p_tag);
    }
}

Tag *tag_find(Tag *tag, char *key) {
    while(tag != NULL) {
        if(strcmp(tag->key, key) == 0) {
            return tag;
        }
        tag = tag->n_tag;
    }
    return NULL;
}

Tag *tag_remove(Tag *first_tag, char *key) {
    if(strcmp(first_tag->key, key) == 0) { // If first tag is one to remove
        return first_tag->n_tag;
    }
    // Tag to remove is not the first link
    Tag *tag = first_tag;
    while(tag->n_tag != NULL) {
        if(strcmp(tag->n_tag->key, key) == 0) {
            Tag *t = tag->n_tag;
            tag->n_tag = tag->n_tag->n_tag; // Remove the link
            free(t);
            break;
        }
        tag = tag->n_tag;
    }
    return first_tag;
}

Tag *tag_construct(char *key, char *value, Tag *n_tag) {
    Tag *tag = st_calloc(1, sizeof(Tag));
    tag->key = stString_copy(key);
    tag->value = stString_copy(value);
    tag->n_tag = n_tag;
    return tag;
}

Tag *tag_parse(char *tag_string, char *delimiter, Tag *p_tag) {
    stList *tag_tokens = stString_splitByString(tag_string, delimiter);
    if (stList_length(tag_tokens) != 2) {
        st_errAbort("Tag not separated by '%s' character: %s\n", delimiter, tag_string);
    }
    Tag *tag = st_calloc(1, sizeof(Tag));
    if(p_tag != NULL) { // If the p_tag is not null, set its n_tag to point at tag
        p_tag->n_tag = tag;
    }
    tag->key = stList_get(tag_tokens, 0);
    tag->value = stList_get(tag_tokens, 1);
    stList_setDestructor(tag_tokens, NULL);
    stList_destruct(tag_tokens);
    return tag;
}

void alignment_row_destruct(Alignment_Row *row) {
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

void alignment_destruct(Alignment *alignment, bool cleanup_rows) {
    Alignment_Row *row = alignment->row;
    if(cleanup_rows) {
        while (row != NULL) {
            Alignment_Row *r = row;
            row = row->n_row;
            alignment_row_destruct(r);
        }
    }
    assert(alignment->column_tags != NULL);  // Clean up column tags
    for(int64_t i=0; i<alignment->column_number; i++) {
        tag_destruct(alignment->column_tags[i]);
    }
    free(alignment->column_tags);
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

bool alignment_row_is_predecessor_2(Alignment_Row **left_row, Alignment_Row **right_row) {
    // Do the rows match - this one is needed to work with the OND aligner which compares pointers to the objects being
    // compared
    return alignment_row_is_predecessor(left_row[0], right_row[0]);
}

void alignment_link_adjacent(Alignment *left_alignment, Alignment *right_alignment, bool allow_row_substitutions) {
    stList *left_rows = get_rows_in_a_list(left_alignment->row);
    stList *right_rows = get_rows_in_a_list(right_alignment->row);
    // get the alignment of the rows
    WFA *wfa = WFA_construct(stList_getBackingArray(left_rows), stList_getBackingArray(right_rows),
                             stList_length(left_rows), stList_length(right_rows),
                             sizeof(void *), (bool (*)(void *, void *))alignment_row_is_predecessor_2, 1,
                             allow_row_substitutions ? 1 : 100000000); // Use unit gap and mismatch costs for the diff
                             // unless we disallow substitutions, in which case use an arbitrarily large mismatch cost
    int64_t aligned_rows[stList_length(left_rows)];
    WFA_get_alignment(wfa, aligned_rows);
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
    // connect up the rows according to the alignment
    for(int64_t i=0; i<stList_length(left_rows); i++) {
        if(aligned_rows[i] != -1) {
            Alignment_Row *left_row = stList_get(left_rows, i);
            Alignment_Row *right_row = stList_get(right_rows, aligned_rows[i]);
            left_row->r_row = right_row;
            right_row->l_row = left_row;
            if(!allow_row_substitutions) {
                assert(alignment_row_is_predecessor(left_row, right_row));
            }
        }
    }
    // clean up
    stList_destruct(left_rows);
    stList_destruct(right_rows);
    WFA_destruct(wfa);
}

int64_t alignment_length(Alignment *alignment) {
    return alignment->column_number;
}

int64_t alignment_total_gap_length(Alignment *left_alignment) {
    Alignment_Row *l_row = left_alignment->row;
    int64_t total_interstitial_gap_length = 0;
    while(l_row != NULL) {
        if(l_row->r_row != NULL && alignment_row_is_predecessor(l_row, l_row->r_row)) {
            int64_t i = l_row->r_row->start - (l_row->start + l_row->length);
            if (i > total_interstitial_gap_length) {
                total_interstitial_gap_length = i;
            }
        }
        l_row = l_row->n_row; // Move to the next left alignment row
    }
    return total_interstitial_gap_length;
}

/*
 * Returns a sequence of tags from the tokens, starting at starting_token
 */
Tag *parse_tags(stList *tokens, int64_t starting_token, char *delimiter) {
    Tag *first_tag = NULL, *tag = NULL;
    if(starting_token < stList_length(tokens)) { // parse first tag pair
        tag = tag_parse(stList_get(tokens, starting_token), delimiter, NULL);
        first_tag = tag;
    }
    for(int64_t i=starting_token+1; i<stList_length(tokens); i++) {
        tag = tag_parse(stList_get(tokens, i), delimiter, tag);
    }
    return first_tag;
}

Tag *parse_header(stList *tokens, char *header_prefix, char *delimiter) {
    if(stList_length(tokens) == 0 || strcmp((char *)stList_get(tokens, 0), header_prefix) != 0) {
        st_errAbort("Header line does not start with %s\n", header_prefix);
    }
    return parse_tags(tokens, 1, delimiter);
}

void write_header(Tag *tag, FILE *fh, char *header_prefix, char *delimiter, char *end) {
    fprintf(fh, "%s", header_prefix);
    while(tag != NULL) {
        fprintf(fh, " %s%s%s", tag->key, delimiter, tag->value);
        tag = tag->n_tag;
    }
    fprintf(fh, "%s", end);
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

static char *make_run(int64_t length, char c) {
    char gap_alignment[length+1];
    for(int64_t i=0; i<length; i++) {
        gap_alignment[i] = c;
    }
    gap_alignment[length] = '\0';
    return stString_copy(gap_alignment);
}

int64_t make_msa(int64_t string_no, int64_t column_no, int64_t max_alignment_length,
                 int64_t msa[string_no][column_no], char *strings[string_no], int64_t string_lengths[string_no],
                 char msa_strings[string_no][max_alignment_length]) {
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
    int64_t msa[string_no][longest_string_length]; // A longest sequence x sequence number sized integer matrix,
    // each entry is the index of the position aligned at that node (or -1 if unaligned)
    char *row_strings[string_no]; // The gap strings
    int64_t row_string_lengths[string_no];
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
    char msa_strings[string_no][max_alignment_length]; // can not be longer than the minimum of the total length of the strings (where all bases
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
            // Check the rows agree coordinate wise
            assert(strcmp(l_row->sequence_name, l_row->r_row->sequence_name) == 0);
            assert(l_row->strand == l_row->r_row->strand);
            assert(l_row->start + l_row->length <= l_row->r_row->start);

            // Is not a deletion, so merge together two adjacent rows
            assert(l_row->r_row->left_gap_sequence != NULL);
            assert(strlen(l_row->r_row->left_gap_sequence) == interstitial_alignment_length);
            char *bases = stString_print("%s%s%s", l_row->bases, l_row->r_row->left_gap_sequence, l_row->r_row->bases);
            free(l_row->bases); // clean up
            l_row->bases = bases;

            // Update the left row's length coordinate
            int64_t interstitial_bases = l_row->r_row->start - (l_row->start + l_row->length);
            l_row->length += interstitial_bases + l_row->r_row->length;
            // Update the l_row's r_row pointer...
            if(l_row->r_row->r_row != NULL) { // Check pointers are correct
                assert(l_row->r_row->r_row->l_row == l_row->r_row);
            }
            l_row->r_row = l_row->r_row->r_row;
            if(l_row->r_row != NULL) {
                l_row->r_row->l_row = l_row;
            }
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


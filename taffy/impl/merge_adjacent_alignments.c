#include "taf.h"
#include "ond.h"
#include "sonLib.h"
#include "abpoa.h"
#include <stdlib.h>

/*
 * Method to merge together two adjacent alignments.
 */

static char *make_run(int64_t length, char c) {
    char *run = st_malloc(sizeof(char) * (length + 1));
    memset(run, c, length);
    run[length] = '\0';
    return run;
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

// char <--> uint8_t conversion copied over from abPOA example
// AaCcGgTtNn ==> 0,1,2,3,4
static unsigned char nst_nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
static const char nst_nt256_table[256] = {
       'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
};

static inline char msa_to_base(uint8_t n) {
    return (char)nst_nt256_table[n];
}

static inline uint8_t msa_to_byte(char c) {
    return nst_nt4_table[(int)c];
}

// from cactus/bar/impl/poaBarAligner.c
static abpoa_para_t* construct_abpoa_params(void) {
    abpoa_para_t *abpt = abpoa_init_para();

    // output options
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 0; // generate consensus sequence, set 0 to disable

    // alignment mode. 0:global alignment, 1:local, 2:extension
    // only global works
    abpt->align_mode = ABPOA_GLOBAL_MODE;

    // banding parameters
    abpt->wb = 300;
    abpt->wf = 0.05;

    // gap scoring model
    abpt->gap_open1 = 400;
    abpt->gap_ext1 = 30;
    abpt->gap_open2 = 1200;
    abpt->gap_ext2 = 1;
    
    // seeding paramters
    abpt->disable_seeding = 1;
    abpt->k = 19;
    abpt->w = 10;
    abpt->min_w = 500;

    // progressive toggle
    abpt->progressive_poa = 1;

    // generate the substitution matrix
    abpt->use_score_matrix = 0;
    abpoa_post_set_para(abpt);

    // optionally override the substitution matrix
    char *submat_string = stString_copy("91 -114 -61 -123 -100 -114 100 -125 -61 -100 -61 -125 100 -114 -100 -123 -61 -114 91 -100 -100 -100 -100 -100 100");
    if (submat_string && strlen(submat_string) > 0) {
        // Note, this will be used to explicitly override abpoa's subsitution matrix just before aligning
        abpt->use_score_matrix = 1;
        assert(abpt->m == 5);
        int count = 0;
        for (char* val = strtok(submat_string, " "); val != NULL; val = strtok(NULL, " ")) {
            abpt->mat[count++] = atoi(val);
        }
        assert(count == 25);
        int i; abpt->min_mis = 0, abpt->max_mat = 0;
        for (i = 0; i < abpt->m * abpt->m; ++i) {
            if (abpt->mat[i] > abpt->max_mat)
                abpt->max_mat = abpt->mat[i];
            if (-abpt->mat[i] > abpt->min_mis) 
                abpt->min_mis = -abpt->mat[i];
        }
    }
    free(submat_string);
    return abpt;
}

// input are row ranks, comparison is reverse-sort on length
static int len_rev_cmp(const void* i1, const void* i2, void* length_list) {
    int64_t len1 = (int64_t)stList_get((stList*)length_list, (int64_t)i1);
    int64_t len2 = (int64_t)stList_get((stList*)length_list, (int64_t)i2);    
    return (int64_t)len2 < (int64_t)len1 ? -1 : ((int64_t)len2 > (int64_t)len1 ? 1 : 0);
}

/*
 * abPOA is set up once and reused for every gap alignment. Building the parameters is expensive out
 * of proportion to the work: abpoa_post_set_para refills a 65536 entry global lookup table each
 * time, and doing that once per merge dominated the run time of taffy norm on alignments made of
 * many small blocks. The parameters are constant, and abpoa_msa resets the handle as its first act,
 * so one instance serves every call. This is not thread safe, and neither is abPOA's use of the
 * global tables that make it worth doing.
 */
static abpoa_t *shared_abpoa = NULL;
static abpoa_para_t *shared_abpoa_params = NULL;

static void free_shared_abpoa(void) {
    if(shared_abpoa != NULL) {
        abpoa_free(shared_abpoa);
        shared_abpoa = NULL;
    }
    if(shared_abpoa_params != NULL) {
        abpoa_free_para(shared_abpoa_params);
        shared_abpoa_params = NULL;
    }
}

static void get_shared_abpoa(abpoa_t **ab, abpoa_para_t **abpt) {
    if(shared_abpoa == NULL) {
        shared_abpoa = abpoa_init();
        shared_abpoa_params = construct_abpoa_params();
        atexit(free_shared_abpoa); // so the singleton is not reported as a leak
    }
    *ab = shared_abpoa;
    *abpt = shared_abpoa_params;
}

int64_t align_interstitial_gaps_abpoa(Alignment *alignment) {
    /*
     * Align the sequences that lie within the gaps between two adjacent blocks using abPOA.
     * Return the length of the interstitial alignment.
     */

    // sorted input really helps abpoa, so we make some arrays for random access
    stList *row_list = stList_construct2(alignment->row_number);
    stList *length_list = stList_construct2(alignment->row_number);
    stList *order_list = stList_construct2(alignment->row_number);
    
    // Add any missing gap strings in
    int64_t row_idx = 0;
    int seq_no = 0;
    for (Alignment_Row *row = alignment->row; row != NULL; row = row->n_row) {
        if(row->l_row != NULL && alignment_row_is_predecessor(row->l_row, row) && row->left_gap_sequence == NULL) {
            row->left_gap_sequence = make_run(row->start - (row->l_row->start + row->l_row->length), 'N');
        }
        int64_t row_length = 0;
        if (row->left_gap_sequence != NULL && row->left_gap_sequence[0] != '\0') {
            row_length = strlen(row->left_gap_sequence);
            ++seq_no;
        }
        stList_set(row_list, row_idx, (void*)row);
        stList_set(length_list, row_idx, (void*)row_length);
        stList_set(order_list, row_idx, (void*)row_idx);                            
        ++row_idx;
    }

    if (seq_no == 0) {
        stList_destruct(row_list);
        stList_destruct(length_list);
        stList_destruct(order_list);
        return 0;
    }

    // sort by length decreasing
    stList_sort2(order_list, len_rev_cmp, (void*)length_list);

    abpoa_t *ab;
    abpoa_para_t *abpoa_params;
    get_shared_abpoa(&ab, &abpoa_params);

    // convert into abpoa input matrix    
    int *seq_lens = (int*)st_calloc(seq_no, sizeof(int));
    uint8_t **bseqs = (uint8_t**)st_calloc(seq_no, sizeof(uint8_t*));
    for (int64_t i = 0; i < seq_no; ++i) {
        row_idx = (int64_t)stList_get(order_list, i);
        Alignment_Row *row = stList_get(row_list, row_idx);
        assert(row->left_gap_sequence != NULL && row->left_gap_sequence[0] != '\0');
        seq_lens[i] = strlen(row->left_gap_sequence);
        bseqs[i] = (uint8_t*)st_calloc(seq_lens[i], sizeof(uint8_t));
        for (int64_t col = 0; col < seq_lens[i]; ++col) {
            bseqs[i][col] = msa_to_byte(row->left_gap_sequence[col]);
        }
    }

    // run abpoa: todo try sorting on length
    // abpoa_msa resets the shared handle as its first act, but it returns before doing so when
    // handed no sequences, which would leave the previous call's alignment for us to read back. The
    // seq_no == 0 return above means that cannot happen; this is here so it stays that way.
    assert(seq_no > 0);
    abpoa_msa(ab, abpoa_params, seq_no, NULL, seq_lens, bseqs, NULL, NULL);

    // copy the results from the abpoa matrix back into the left_gap_sequences
    int64_t msa_length = ab->abc->msa_len;
    for (int64_t i = 0; i < alignment->row_number; ++i) {
        row_idx = (int64_t)stList_get(order_list, i);
        Alignment_Row *row = stList_get(row_list, row_idx);
        free(row->left_gap_sequence);
        row->left_gap_sequence = (char*)st_calloc(msa_length + 1, sizeof(char));
        if (i < seq_no) {
            for (int64_t col = 0; col < msa_length; ++col) {
                row->left_gap_sequence[col] = msa_to_base(ab->abc->msa_base[i][col]);
            }
        } else {
            for (int64_t col = 0; col < msa_length; ++col) {
                row->left_gap_sequence[col] = '-';
            }            
        }
        row->left_gap_sequence[msa_length] = '\0';
    }

    stList_destruct(row_list);
    stList_destruct(length_list);
    stList_destruct(order_list);

    // The handle and parameters are shared, so they are left for get_shared_abpoa to clean up
    free(seq_lens);
    for (int64_t i = 0; i < seq_no; ++i) {
        free(bseqs[i]);
    }
    free(bseqs);

    return msa_length;
}

/*
 * Join up to three strings whose lengths are already known. The row strings here are all exactly as
 * long as their block's column count, so going through stString_print only to have vsnprintf walk
 * them again with strlen is wasted work, and this is the innermost loop of a merge.
 */
static char *concat_known(const char *a, int64_t a_length, const char *b, int64_t b_length,
                          const char *c, int64_t c_length) {
    char *joined = st_malloc(sizeof(char) * (a_length + b_length + c_length + 1));
    memcpy(joined, a, a_length);
    memcpy(joined + a_length, b, b_length);
    memcpy(joined + a_length + b_length, c, c_length);
    joined[a_length + b_length + c_length] = '\0';
    return joined;
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
    int64_t interstitial_alignment_length = align_interstitial_gaps_abpoa(right_alignment);

    // Now finally extend the left alignment rows to include the right alignment rows
    Alignment_Row *l_row = left_alignment->row;
    int64_t left_column_number = left_alignment->column_number; // every left row is this long
    int64_t right_gap_length = right_alignment->column_number + interstitial_alignment_length;
    char *right_gap = make_run(right_gap_length, '-'); // any trailing bases needed
    while(l_row != NULL) {
        if(l_row->r_row == NULL) {
            // Is a deletion, so add in trailing gaps equal in length to the right alignment length plus any interstitial
            // gap
            char *bases = concat_known(l_row->bases, left_column_number, right_gap, right_gap_length, "", 0);
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
            char *bases = concat_known(l_row->bases, left_column_number,
                                       r_row->left_gap_sequence, interstitial_alignment_length,
                                       r_row->bases, right_alignment->column_number);
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
            right_alignment->column_tags[i] = NULL; // The merged alignment owns them now, so the
            // alignment_destruct below must not free them out from under it
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



#ifndef STTAF_H_
#define STTAF_H_

#include "sonLib.h"
#include "line_iterator.h"

/*
 * Structures to represent blocks of an alignment
 */
typedef struct _row Alignment_Row;

typedef struct _alignment {
    int64_t row_number; // Convenient counter of number rows in the alignment
    Alignment_Row *row; // An alignment is just a sequence of rows
} Alignment;

struct _row { // Each row encodes the information about an aligned sequence
    char *sequence_name; // name of sequence
    int64_t start, length, sequence_length; // zero based, half open coordinates
    bool strand; // nonzero is "+" else "-"
    char *bases; // [A-Za-z*+]* string of length "length"
    char *left_gap_sequence; // Optional interstitial gap sequence, which is the unaligned substring between this
    // sequence and the end of the previous block - may be NULL if not specified or if zero length
    Alignment_Row *l_row;  // connection to a left row (may be NULL)
    Alignment_Row *r_row;  // connection to a right row (may be NULL)
    Alignment_Row *n_row;  // the next row in the alignment
};

/*
 * Clean up the memory for an alignment
 */
void alignment_destruct(Alignment *alignment);

/*
 * Cleanup a row
 */
void Alignment_Row_destruct(Alignment_Row *row);

/*
 * Returns non-zero if left_row represents a substring on the same contig and strand as right_row, but
 * immediately before
 */
bool alignment_row_is_predecessor(Alignment_Row *left_row, Alignment_Row *right_row);

/*
 * Parse a header line that mist start with the header_prefix and then be composed of a series of key value tags,
 * each delimited by the delimiter word.
 */
stList *parse_header(stList *tokens, char *header_prefix, char *delimiter);

/*
 * Use the O(ND) alignment to diff the rows between two alignments and connect together their rows
 * so that we can determine which rows in the right_alignment are a continuation of rows in the
 * left_alignment. We use this for efficiently outputting TAF.
 */
void alignment_link_adjacent(Alignment *left_alignment, Alignment *right_alignment, bool allow_row_substitutions);

/*
 * Gets the number of columns in the alignment
 */
int64_t alignment_length(Alignment *alignment);

/*
 * Gets the sum of the interstitial gaps between this block and the next one. If align_gap_sequences
 * is non-zero the length of longest individual sequence is returned instead.
 */
int64_t alignment_total_gap_length(Alignment *left_alignment, bool align_gap_sequences);

/*
 * Number of shared rows between two alignments
 */
int64_t alignment_number_of_common_rows(Alignment *left_alignment, Alignment *right_alignment);

/*
 * Merge together adjacent blocks into one alignment. Requires that the alignment
 * rows are linked together (e.g. with alignment_link_adjacent). Destroys the input
 * alignments in the process and returns a merged alignment.
 *
 * If align_gap_sequences is non-zero then sequences are aligned.
 */
Alignment *alignment_merge_adjacent(Alignment *left_alignment, Alignment *right_alignment,
                                    bool align_gap_sequences);

/*
 * Read a maf header line
 */
stList *maf_read_header(FILE *fh);

/*
 * Read a maf alignment block from the file stream. Return NULL if none available
 */
Alignment *maf_read_block(FILE *fh);

/*
 * Write a maf header line
 */
void maf_write_header(stList *tags, FILE *fh);

/*
 * Write a maf block
 */
void maf_write_block(Alignment *alignment, FILE *fh);

/*
 * Read a taf header line
 */
stList *taf_read_header(LI *li);

/*
 * Read a taf block - that is a column with column coordinates and all subsequent coordinate-less columns that
 * are considered to be part of the block.
 */
Alignment *taf_read_block(Alignment *p_block, bool run_length_encode_bases, LI *li);

/*
 * Write a taf header line
 */
void taf_write_header(stList *tags, FILE *fh);

/*
 * Write a taf block.
 */
void taf_write_block(Alignment *p_alignment, Alignment *alignment, bool run_length_encode_bases, FILE *fh);

#endif /* STTAF_H_ */


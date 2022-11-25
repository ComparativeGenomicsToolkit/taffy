#ifndef STTAF_H_
#define STTAF_H_

#include <inttypes.h>
#include <stdio.h>
#include <stdbool.h>

#include "line_iterator.h"

/*
 * Structures to represent blocks of an alignment
 */

typedef struct _tag Tag;

struct _tag {
    char *key;
    char *value;
    Tag *n_tag; // Next tag in list
};

typedef struct _row Alignment_Row;

typedef struct _alignment {
    int64_t row_number; // Convenient counter of number rows in the alignment
    int64_t column_number; // Convenient counter of number of columns in this alignment
    Alignment_Row *row; // An alignment is just a sequence of rows
    Tag **column_tags; // The tags for each column, each stored as a sequence of tags
} Alignment;

struct _row { // Each row encodes the information about an aligned sequence
    char *sequence_name; // name of sequence
    int64_t start, length, sequence_length; // zero based, half open coordinates
    bool strand; // nonzero is "+" else "-"
    char *bases; // [A-Za-z*+]* string of bases and gaps representing the alignment of the row
    char *left_gap_sequence; // Optional interstitial gap sequence, which is the unaligned substring between this
    // sequence and the end of the previous block - may be NULL if not specified or if zero length
    Alignment_Row *l_row;  // connection to a left row (may be NULL)
    Alignment_Row *r_row;  // connection to a right row (may be NULL)
    Alignment_Row *n_row;  // the next row in the alignment
    int64_t bases_since_coordinates_reported; // this number is used by taf write coordinates to
    // indicate how many bases ago were the row's coordinates printed
};

/*
 * Make a tag
 */
Tag *tag_construct(char *key, char *value, Tag *n_tag);

/*
 * Cleanup a sequence of tags
 */
void tag_destruct(Tag *tag);

/*
 * Find a tag with a given key
 */
Tag *tag_find(Tag *tag, char *key);

/*
 * Remove a tag, cleaning it up. Returns the modified sequence.
 */
Tag *tag_remove(Tag *first_tag, char *key);

/*
 * Parse a tag from a string.
 *
 * If p_tag is not null will set p_tag->n_tag = tag, where tag is the parsed tag.
 */
Tag *tag_parse(char *tag_string, char *delimiter, Tag *p_tag);

/*
 * Clean up the memory for an alignment
 */
void alignment_destruct(Alignment *alignment, bool cleanup_rows);

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
 * Gets the max length of an interstitial gap sequence between this block and the next one.
 */
int64_t alignment_total_gap_length(Alignment *left_alignment);

/*
 * Number of shared rows between two alignments
 */
int64_t alignment_number_of_common_rows(Alignment *left_alignment, Alignment *right_alignment);

/*
 * Merge together adjacent blocks into one alignment. Requires that the alignment
 * rows are linked together (e.g. with alignment_link_adjacent). Destroys the input
 * alignments in the process and returns a merged alignment. If there are interstitial
 * sequences between the blocks, aligns these sequences together.
 */
Alignment *alignment_merge_adjacent(Alignment *left_alignment, Alignment *right_alignment);

/*
 * Cleanup a row
 */
void alignment_row_destruct(Alignment_Row *row);

/*
 * Returns non-zero if left_row represents a substring on the same contig and strand as right_row, but
 * immediately before
 */
bool alignment_row_is_predecessor(Alignment_Row *left_row, Alignment_Row *right_row);

/*
 * Read a maf header line
 */
Tag *maf_read_header(FILE *fh);

/*
 * Read a maf alignment block from the file stream. Return NULL if none available
 */
Alignment *maf_read_block(FILE *fh);

/*
 * Write a maf header line
 */
void maf_write_header(Tag *tag, FILE *fh);

/*
 * Write a maf block
 */
void maf_write_block(Alignment *alignment, FILE *fh);


/*
 * Read a taf header line
 */
Tag *taf_read_header(LI *li);

/*
 * Read a taf block - that is a column with column coordinates and all subsequent coordinate-less columns that
 * are considered to be part of the block.
 */
Alignment *taf_read_block(Alignment *p_block, bool run_length_encode_bases, LI *li);

/*
 * Write a taf header line
 */
void taf_write_header(Tag *tag, FILE *fh);

/*
 * Write a taf block.
 */
void taf_write_block(Alignment *p_alignment, Alignment *alignment, bool run_length_encode_bases,
                     int64_t repeat_coordinates_every_n_columns, FILE *fh);

#endif /* STTAF_H_ */


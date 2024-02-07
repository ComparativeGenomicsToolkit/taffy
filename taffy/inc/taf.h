#ifndef STTAF_H_
#define STTAF_H_

#include <inttypes.h>
#include <stdio.h>
#include <stdbool.h>
#include "sonLib.h"
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
    int64_t start, length, sequence_length; // zero based, half open coordinates, length is the number of non gap bases in the row
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
 * Add nucleotide coloring to a character for pretty printing
 */
char *color_base_char(char base);

/*
 * Convert a nucleotide string into a colored string suitable for pretty printing.
 */
char *color_base_string(char *bases, int64_t length);

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
 * Get the rows of the alignment in a list.
 */
stList *alignment_get_rows_in_a_list(Alignment_Row *row);

/*
 * Set the rows in the alignment given a list of rows
 */
void alignment_set_rows(Alignment *alignment, stList *rows);

/*
 * Read a column of the alignment into the buffer. The buffer must be initialized and be at least
 * of length alignment->row_number.
 */
void alignment_get_column_in_buffer(Alignment *alignment, int64_t column_index, char *buffer);

/*
 * Read a column of the alignment and return as a string
 */
char *alignment_get_column(Alignment *alignment, int64_t column_index);

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
 * Returns a pretty-printed string representing the row. Useful for debugging.
 */
char *alignment_row_to_string(Alignment_Row *row);

/*
 * Returns a pretty-printed string representing the alignment. Useful for debugging.
 */
char *alignment_to_string(Alignment *alignment);

/*
 * Replace bases that match the reference with a mask character.
 */
void alignment_mask_reference_bases(Alignment *alignment, char mask_char);

/*
 * Replace bases that match their ancestral lineage with a mask character
 */
void alignment_show_only_lineage_differences(Alignment *alignment, char mask_char, stList *sequence_prefixes, stList *tree_nodes);

/*
 * Read a maf header line
 */
Tag *maf_read_header(LI *li);

/*
 * Read a maf alignment block from the file stream. Return NULL if none available
 */
Alignment *maf_read_block(LI *li);

/*
 * Write a maf header line
 */
void maf_write_header(Tag *tag, LW *lw);

/*
 * Write a maf block
 */
void maf_write_block(Alignment *alignment, LW *lw);

/*
 * As maf write block, but with option to output pretty colored bases.
 */
void maf_write_block2(Alignment *alignment, LW *lw, bool color_bases);

/*
 * Write a block as PAF. Each PAF row reflects a pairwise alignment in the block.  The all_to_all flag
 * toggles whether we write every possible pairwise alignment, or just each non-ref to ref alignment
 * where ref is the first row in the block
 */
void paf_write_block(Alignment *alignment, LW *lw, bool all_to_all, bool cs_cigar);

/*
 * Read a taf header line
 */
Tag *taf_read_header(LI *li);

/*
 * Read a taf header line, check if the run_length_encode_bases flag is set.
 */
Tag *taf_read_header_2(LI *li, bool *run_length_encode_bases);

/*
 * Read a taf block - that is a column with column coordinates and all subsequent coordinate-less columns that
 * are considered to be part of the block.
 */
Alignment *taf_read_block(Alignment *p_block, bool run_length_encode_bases, LI *li);

/*
 * Write a taf header line
 */
void taf_write_header(Tag *tag, LW *lw);

/*
 * Write a taf block.
 */
void taf_write_block(Alignment *p_alignment, Alignment *alignment, bool run_length_encode_bases,
                     int64_t repeat_coordinates_every_n_columns, LW *lw);

/*
 * As taf write block, but with option to pretty print the output
 */
void taf_write_block2(Alignment *p_alignment, Alignment *alignment, bool run_length_encode_bases,
                      int64_t repeat_coordinates_every_n_columns, LW *lw, bool color_bases);


// the following are low-level functions used in indexing.  they could
// potentially be better put in an "internal" header

/**
 * Check if a tokenized TAF line has coordinates (ie search for ;)
 * If coordinates found, the position of the simicolon in the list is set in j
 */
bool has_coordinates(stList *tokens, int64_t *j);

/**
 * Parse the coordinates coming after an "i" or "s" field on a TAF line
 * Returns the sequence name (newly allocated) and sets the strand and
 * position and length in the parameters. *j is the position of the "i" 
 * or "s" in the tokenized line, and will be incremented for each field read
 */
char *parse_coordinates(int64_t *j, stList *tokens, int64_t *start, bool *strand,
                        int64_t *sequence_length);

/**
 * Sniff file format from header line.  returns:
 *  0: taf
 *  1: maf
 *  2: unknown
 */
int check_input_format(const char *header_line);

/**
 * it turns out just scanning for a "." doesn't work, even with the output of hal2maf.  The reason is
 * that spcecies names can contain "." characters -- at least they do in the VGP alignment.
 * so this function uses knowloedge of the genome names to greedily scan for a prefix ending in "."
 * that corresponds to a genome name in the hal. the extracted genome name is returned if found
 * (it must be freed)
 *
 * (note: this function was originally in taf_add_gap_bases.c where it took a set of species. It's
 * now used by the renaming function int taf_view.c where it uses the hash.  Only one can
 * be specified. If the hal set is used, it fails with an error if it can't find a key. If
 * the hash is used, then it just returns NULL). 
*/
char *extract_genome_name(const char *sequence_name, stSet *hal_species, stHash *genome_name_map);

/**
 * Load a two-column genome name mapping file and return
 * a hash table mapping COLUMN1 -> COLUMN2
 */
stHash *load_genome_name_mapping(char *name_mapping_path);

/**
 * Apply the name mapping to a string. If the string is found in the map,
 * the mapped name is return.  Otherwise NULL is returned.
 * If the input name has a "." in it, only the part before the "."
 * is considered (following ucsc genome.contig convention).
 *
 * So if our mapping contaings hg19 -> Homo_sapiens, then
 *
 * apply_genome_name_mapping("hg19") would return "Homo_sapiens"
 * and
 * apply_genome_name_mapping("hg19.ch10") would return "Homo_sapiens.chr10"
 *
 * If there's more than one dot, it will try splitting at each one so
 *
 * apply_genome_name_mapping("hg19.1.chr10") would return "Homo_sapiens.1.chr10"
 *
 * If a string is returned, it's up to the client to free it
 */
char *apply_genome_name_mapping(stHash *genome_name_map, char *input_name);

/**
 * Apply the name mapping to an alignment block.
 */
void apply_genome_name_mapping_to_alignment(stHash *genome_name_map, Alignment *alignment);

/*
 * Structure to represent a sequence prefix. A sequence of sequence prefixes
 * are used to order the rows in each alignment block..
 */
typedef struct _Sequence_Prefix {
    char *prefix; // The prefix string
    int64_t prefix_length; // Length of the prefix string
    int64_t index; // The index that a sequence matching the prefix should appear in an alignment block
} Sequence_Prefix;

Sequence_Prefix *sequence_prefix_construct(char *prefix, int64_t index);

void sequence_prefix_destruct(Sequence_Prefix *sequence_prefix);

/*
 * Compare two sequence prefixes by their prefix strings
 */
int sequence_prefix_cmp_fn(Sequence_Prefix *p1, Sequence_Prefix *p2);

/*
 * Loads a list of sequence prefixes from a given file handle.
 */
stList *sequence_prefix_load(FILE *sort_fh);

/*
 * Gets the index in the list of the sequence prefix of the given row's sequence name.
 */
int64_t alignment_row_get_closest_sequence_prefix(Alignment_Row *row, stList *prefixes_to_sort_by);

/*
 * Sorts the rows of an alignment according to the given sequence prefixes. Reconnects the rows
 * with the previous alignment in the process.
 */
void alignment_sort_the_rows(Alignment *p_alignment, Alignment *alignment, stList *prefixes_to_sort_by);

/*
 * Removes any rows from the alignment whose sequence name prefix matches a string in the prefixes_to_filtet_by list
 */
void alignment_filter_the_rows(Alignment *alignment, stList *prefixes_to_filter_by);

/*
 * Load sequences in fasta files into a hash from sequence names to sequences
 */
stHash *load_sequences_from_fasta_files(char **seq_file_names, int64_t seq_file_number);

/*
 * Load sequences in hal file into memory.
 */
stSet *load_sequences_from_hal_file(char *hal_file, int *hal_handle);

/*
 * Add any gap strings between representing unaligned sequences between rows of alignment and p_alignment.
 */
void alignment_add_gap_strings(Alignment *p_alignment, Alignment *alignment, stHash *fastas, int hal_handle, stSet *hal_species,
                               int64_t maximum_gap_string_length);

#endif /* STTAF_H_ */


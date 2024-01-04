from cffi import FFI
ffibuilder = FFI()

# cdef() expects a single string declaring the C types, functions and
# globals needed to use the shared object. It must be in valid C syntax.
ffibuilder.cdef("""
    typedef struct BGZF BGZF;

    typedef struct _LW {
        FILE *fh;
        BGZF *bgzf;
    } LW;
    
    LW *LW_construct(FILE *fh, bool use_compression);
    
    void LW_destruct(LW *lw, bool clean_up_file_handle);
    
    int LW_write(LW *lw, const char *string, ...);

    void free(void *ptr);

    FILE *fopen(const char *filename, const char *mode);
    
    int fclose(FILE *stream);

    typedef struct _LI {
        BGZF *bgzf;
        char *line;
        int64_t prev_pos; // position before reading the current buffer
        int64_t pos;      // position after reading the curent buffer    
    } LI;
    
    LI *LI_construct(FILE *fh);

    void LI_destruct(LI *li);
    
    char *LI_peek_at_next_line(LI *li);

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
     * Returns a pretty-printed string representing the row. Useful for debugging.
     */
    char *alignment_row_to_string(Alignment_Row *row);
    
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
     * Returns a pretty-printed string representing the alignment. Useful for debugging.
    */
    char *alignment_to_string(Alignment *alignment);

    /**
     * Sniff file format from header line.  returns:
     *  0: taf
     *  1: maf
     *  2: unknown
     */
    int check_input_format(const char *header_line);

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
    void taf_write_header(Tag *tag, LW *lw);
    
    /*
     * Write a taf block.
     */
    void taf_write_block(Alignment *p_alignment, Alignment *alignment, bool run_length_encode_bases,
                         int64_t repeat_coordinates_every_n_columns, LW *lw);
                         

   typedef struct _Tai Tai;
    
   typedef struct _TaiIt TaiIt;
    
    /*
     * Make an index of a TAF in "idx_fh" 
     * The index is made on the "reference" first contig of each block
     * For each such contig, the index will have one line for each index_block_size
     * region of it found in the TAF. 
     */
    int tai_create(LI *li, FILE* idx_fh, int64_t index_block_size);
    
    /*
     * Load the index from disk
     */
    Tai *tai_load(FILE* idx_fh, bool maf);
    
    /*
     * Free the index
     */
    void tai_destruct(Tai* idx);
    
    /*
     * Query the taf index
     * start is 0-based
     * length can be -1 for everything
     */
    TaiIt *tai_iterator(Tai* idx, LI *li, bool run_length_encode_bases, const char *contig, int64_t start, int64_t length);
    
    /*
     * Iterate through a region as obtained via the iterator
     */
    Alignment *tai_next(TaiIt *tai_it, LI *li);
                        
    /*
     * Free a tai iterator
     */
    void tai_iterator_destruct(TaiIt *tai_it);
""")

# set_source() gives the name of the python extension module to
# produce, and some C source code as a string.  This C code needs
# to make the declarated functions, types and globals available,
# so it is often just the "#include".
ffibuilder.set_source("taffy._taffy_cffi",
                      """
                           #include <stdio.h>
                           #include <stdlib.h>
                           #include "htslib/bgzf.h"
                           #include "htslib/kstring.h"
                           #include "taf.h" 
                           #include "line_iterator.h" 
                           #include "tai.h"
                      """,
                      include_dirs=["taffy/submodules/sonLib/externalTools/cutest",
                                    "taffy/submodules/sonLib/C/inc",
                                    "taffy/inc"],
                      sources=["taffy/submodules/sonLib/C/impl/stSafeC.c",
                               "taffy/submodules/sonLib/C/impl/sonLibCommon.c",
                               "taffy/submodules/sonLib/C/impl/sonLibRandom.c",
                               "taffy/submodules/sonLib/C/impl/sonLibExcept.c",
                               "taffy/submodules/sonLib/C/impl/sonLibString.c",
                               "taffy/submodules/sonLib/C/impl/hashTableC_itr.c",
                               "taffy/submodules/sonLib/C/impl/hashTableC.c",
                               "taffy/submodules/sonLib/C/impl/sonLibHash.c",
                               "taffy/submodules/sonLib/C/impl/avl.c",
                               "taffy/submodules/sonLib/C/impl/sonLibSortedSet.c",
                               "taffy/submodules/sonLib/C/impl/sonLibSet.c",
                               "taffy/submodules/sonLib/C/impl/sonLibList.c",
                               "taffy/submodules/sonLib/C/impl/sonLibFile.c",
                               "taffy/impl/line_iterator.c",
                               "taffy/impl/alignment_block.c",
                               "taffy/impl/merge_adjacent_alignments.c",
                               "taffy/impl/maf.c",
                               "taffy/impl/ond.c",
                               "taffy/impl/taf.c",
                               "taffy/impl/tai.c",
                               ],
                      extra_compile_args=["-DUSE_HTSLIB"],
                      libraries=["hts"],
                      )

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

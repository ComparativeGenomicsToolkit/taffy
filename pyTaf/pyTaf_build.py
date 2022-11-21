from cffi import FFI
import pathlib
ffibuilder = FFI()

# cdef() expects a single string declaring the C types, functions and
# globals needed to use the shared object. It must be in valid C syntax.
ffibuilder.cdef("""
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

    Alignment *maf_read_block(FILE *fh);
""")

# set_source() gives the name of the python extension module to
# produce, and some C source code as a string.  This C code needs
# to make the declarated functions, types and globals available,
# so it is often just the "#include".
ffibuilder.set_source("_pyTaf_cffi",
                      """
                           #include "inc/taf.h" // the C header of the library
                      """,
                      libraries=["sonLib", "stTaf"],   # library name, for the linker
                      library_dirs=[(pathlib.Path().absolute() / "lib").as_posix()])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
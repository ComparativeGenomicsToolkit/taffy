/*
 * taf sort: Sort the rows of a TAF file using a given ordering
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "tai.h"
#include "sonLib.h"
#include <getopt.h>
#include <time.h>

/*
 * Structure to represent a sequence prefix. A sequence of sequence prefixes
 * are used to order the rows in each alignment block..
 */
typedef struct _Sequence_Prefix {
    char *prefix; // The prefix string
    int64_t prefix_length; // Length of the prefix string
    int64_t index; // The index that a sequence matching the prefix should appear in an alignment block
} Sequence_Prefix;

static void sequence_prefix_destruct(Sequence_Prefix *sequence_prefix) {
    free(sequence_prefix->prefix);
    free(sequence_prefix);
}

static int sequence_prefix_cmp_fn(Sequence_Prefix *p1, Sequence_Prefix *p2) {
    return strcmp(p1->prefix, p2->prefix);
}

static stList *load_sort_file(FILE *sort_fh) {
    stList *prefixes_to_sort_by = stList_construct3(0, (void (*)(void *))sequence_prefix_destruct);
    int64_t index = 0;
    char *line;
    while((line = stFile_getLineFromFile(sort_fh)) != NULL) {
        stList *tokens = stString_split(line);
        if(stList_length(tokens) != 1) {
            st_errAbort("Expected exactly one string in sort file on line: %s", line);
        }
        Sequence_Prefix *sequence_prefix = st_calloc(1, sizeof(Sequence_Prefix));
        sequence_prefix->prefix = stList_pop(tokens);
        sequence_prefix->prefix_length = strlen(sequence_prefix->prefix);
        if(sequence_prefix->prefix_length == 0) {
            st_errAbort("Found an empty sequence prefix: %s", line);
        }
        sequence_prefix->index = index++;
        stList_append(prefixes_to_sort_by, sequence_prefix);
        // Clean up
        free(line);
        stList_destruct(tokens);
    }
    stList_sort(prefixes_to_sort_by, (int (*)(const void *, const void *))sequence_prefix_cmp_fn);
    return prefixes_to_sort_by;
}

static int get_closest_prefix_cmp_fn(char *sequence_name, Sequence_Prefix *sp) {
    int64_t i = strcmp(sequence_name, sp->prefix);
    if(i > 0) { // If sequence_name is lexicographically larger than sequence_prefix could
        // be a prefix (can not be a prefix is i < 0)
        for(int64_t j=0; j<sp->prefix_length; j++) {
            if(sequence_name[j] != sp->prefix[j]) {
                return 1;
            }
        }
        return 0;
    }
    return i;
}

static int64_t get_closest_prefix(Alignment_Row *row, stList *prefixes_to_sort_by) {
    // Binary search the sequence name
    Sequence_Prefix *sp = stList_binarySearch(prefixes_to_sort_by, row,
                                              (int (*)(const void *a, const void *b))get_closest_prefix_cmp_fn);
    if(sp == NULL) {
        st_logDebug("Did not find a valid prefix to match: %s\n", row->sequence_name);
    }
    return sp != NULL ? sp->index : -1; // Sequences that don't have a match will appear first in the sort
}

static int alignment_sequence_prefix_cmp_fn(Alignment_Row *a1, Alignment_Row *a2,
                                            stList *prefixes_to_sort_by) {
    int i = get_closest_prefix(a1, prefixes_to_sort_by);
    int j = get_closest_prefix(a2, prefixes_to_sort_by);
    return i < j ? -1 : (i > j ? 1 : strcmp(a1->sequence_name, a2->sequence_name));
}

/*
 * Sorts the rows of the alignment according to the given prefixes
 */
static void sort_the_rows(Alignment *p_alignment, Alignment *alignment, stList *prefixes_to_sort_by) {
    // Get the rows
    stList *rows = alignment_get_rows_in_a_list(alignment->row);
    assert(stList_length(rows) == alignment->row_number); // Quick sanity check

    // Sort the rows by the prefix ordering
    stList_sort2(rows, (int (*)(const void *, const void *, void *))alignment_sequence_prefix_cmp_fn, rows);

    // Re-connect the rows
    alignment_set_rows(alignment, rows);

    // Reset the alignment of the rows with the prior row
    alignment_link_adjacent(p_alignment, alignment, 1);
}

static void usage(void) {
    fprintf(stderr, "taffy sort [options]\n");
    fprintf(stderr, "Sort the rows of the TAF alignment file in a specified order\n");
    fprintf(stderr, "-i --inputFile : Input TAF or MAF file. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output file. If not specified outputs to stdout\n");
    fprintf(stderr, "-n --sortFile : File in which each line is a prefix of a sequence name. Rows are sorted accordingly, \n"
                    "with any ties broken by lexicographic sort of the suffixes. Can not be None\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int taf_sort_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *input_file = NULL;
    char *output_file = NULL;
    char *sort_file = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "sortFile", required_argument, 0, 'n' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:n:h", long_options, &option_index);
        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'i':
                input_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'n':
                sort_file = optarg;
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Input file string : %s\n", input_file);
    st_logInfo("Output file string : %s\n", output_file);
    st_logInfo("Sort file string : %s\n", sort_file);

    //////////////////////////////////////////////
    // Read in the taf/maf blocks and sort order file
    //////////////////////////////////////////////

    // Input taf
    FILE *input = input_file == NULL ? stdin : fopen(input_file, "r");
    if (input == NULL) {
        fprintf(stderr, "Unable to open input file: %s\n", input_file);
        return 1;
    }
    LI *li = LI_construct(input);

    // Output taf
    FILE *output_fh = output_file == NULL ? stdout : fopen(output_file, "w");
    if (output_fh == NULL) {
        fprintf(stderr, "Unable to open output file: %s\n", output_file);
        return 1;
    }
    LW *output = LW_construct(output_fh, 0);

    // Sort file
    if (sort_file == NULL) {
        fprintf(stderr, "No sort file specified!\n");
        return 1;
    }
    FILE *sort_fh = fopen(sort_file, "r");
    if (sort_fh == NULL) {
        fprintf(stderr, "Unable to open sort file: %s\n", sort_file);
        return 1;
    }
    stList *prefixes_to_sort_by = load_sort_file(sort_fh);

    // Parse the header
    Tag *tag = taf_read_header(li);
    Tag *t = tag_find(tag, "run_length_encode_bases");
    bool run_length_encode_bases = t != NULL && strcmp(t->value, "1") == 0;

    // Write the header
    taf_write_header(tag, output);
    tag_destruct(tag);

    // Write the alignment blocks
    Alignment *alignment, *p_alignment = NULL;
    while((alignment = taf_read_block(p_alignment, li, run_length_encode_bases)) != NULL) {
        // Sort the alignment block rows
        sort_the_rows(p_alignment, alignment, prefixes_to_sort_by);
        // Write the block
        taf_write_block(p_alignment, alignment, run_length_encode_bases, -1, output); // Write the block
        p_alignment = alignment;
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    LI_destruct(li);
    if(input_file != NULL) {
        fclose(input);
    }
    LW_destruct(output, output_file != NULL);

    st_logInfo("taffy sort is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}


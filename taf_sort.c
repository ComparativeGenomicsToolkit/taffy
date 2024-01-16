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
    stList *prefixes_to_sort_by = sequence_prefix_load(sort_fh);
    st_logInfo("Loaded the sort file, got %i rows\n", (int)stList_length(prefixes_to_sort_by));

    // Parse the header
    Tag *tag = taf_read_header(li);
    Tag *t = tag_find(tag, "run_length_encode_bases");
    bool run_length_encode_bases = t != NULL && strcmp(t->value, "1") == 0;

    // Write the header
    taf_write_header(tag, output);
    tag_destruct(tag);

    // Write the alignment blocks
    Alignment *alignment, *p_alignment = NULL, *pp_alignment = NULL;
    while((alignment = taf_read_block(p_alignment, run_length_encode_bases, li)) != NULL) {
        if(p_alignment) {
            // Sort the alignment block rows
            alignment_sort_the_rows(pp_alignment, p_alignment, prefixes_to_sort_by);
            // Write the block
            taf_write_block(pp_alignment, p_alignment, run_length_encode_bases, -1, output); // Write the block
        }
        pp_alignment = p_alignment;
        p_alignment = alignment;
    }
    if(p_alignment) { // Write the final block
        alignment_sort_the_rows(pp_alignment, p_alignment, prefixes_to_sort_by);
        taf_write_block(pp_alignment, p_alignment, run_length_encode_bases, -1, output); // Write the block
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


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
                    "with any ties broken by lexicographic sort of the suffixes.\n");
    fprintf(stderr, "-f --filterFile : Remove any rows with sequences matching a prefix in this file\n");
    fprintf(stderr, "-p --padFile : Add a padding row for any sequence in this file that is not a prefix of an existing row\n");
    fprintf(stderr, "-d --dupFilterFile : Remove duplicate sequences matching any prefix in this file\n");
    fprintf(stderr, "-r --ignoreFirstRow : Do not consider the first row of each maf block - useful if wanting to "
                    "preserve a reference sequence\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

stList *load_sort_file(char *sort_file) {
    if (sort_file == NULL) {
        return NULL;
    }
    FILE *sort_fh = fopen(sort_file, "r");
    if (sort_fh == NULL) {
        fprintf(stderr, "Unable to open sort/filter file: %s\n", sort_file);
        exit(1);
    }
    stList *prefixes_to_sort_by = sequence_prefix_load(sort_fh);
    st_logInfo("Loaded the sort/filter file, got %i rows\n", (int) stList_length(prefixes_to_sort_by));
    return prefixes_to_sort_by;
}



void process_alignment_block(Alignment *pp_alignment, Alignment *p_alignment, stList *prefixes_to_filter_by,
                             stList * prefixes_to_pad, stList *prefixes_to_sort_by, stList *prefixes_to_dup_filter,
                             bool run_length_encode_bases, bool ignore_first_row, LW *output) {
    if(p_alignment) {
        if(prefixes_to_filter_by) { //Remove rows matching a prefix
            alignment_filter_the_rows(p_alignment, prefixes_to_filter_by, ignore_first_row);
        }
        if(prefixes_to_pad) {
            alignment_pad_the_rows(prefixes_to_sort_by ? NULL : pp_alignment,p_alignment,
                                   prefixes_to_pad); // Note we only reconnect to the prior alignment
                                   // if we are not sorting - this is for efficiency (avoid runing O(ND) twice)
        }
        if(prefixes_to_sort_by) { // Sort the alignment block rows
            alignment_sort_the_rows(pp_alignment, p_alignment,
                                    prefixes_to_sort_by, ignore_first_row);
        }
        if(prefixes_to_dup_filter) { // Remove duplicate rows
            alignment_filter_duplicate_rows(p_alignment, prefixes_to_dup_filter, ignore_first_row);
        }
        // Write the block
        taf_write_block(pp_alignment, p_alignment,
                        run_length_encode_bases, -1, output); // Write the block
    }
    if(pp_alignment != NULL) {
        alignment_destruct(pp_alignment, 1); // Delete the left most block
    }
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
    char *filter_file = NULL;
    char *pad_file = NULL;
    char *dup_filter_file = NULL;
    bool ignore_first_row = 0;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = {{"logLevel",   required_argument, 0, 'l'},
                                               {"inputFile",  required_argument, 0, 'i'},
                                               {"outputFile", required_argument, 0, 'o'},
                                               {"sortFile",   required_argument, 0, 'n'},
                                               {"filterFile", required_argument, 0, 'f'},
                                               {"padFile", required_argument, 0, 'p'},
                                               {"dupFilterFile", required_argument, 0, 'd'},
                                               {"ignoreFirstRow", no_argument, 0, 'r'},
                                               {"help",       no_argument,       0, 'h'},
                                               {0, 0,                            0, 0}};

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:n:hrf:p:d:", long_options, &option_index);
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
            case 'f':
                filter_file = optarg;
                break;
            case 'p':
                pad_file = optarg;
                break;
            case 'd':
                dup_filter_file = optarg;
                break;
            case 'r':
                ignore_first_row = 1;
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
    st_logInfo("Filter file string : %s\n", filter_file);
    st_logInfo("Pad file string : %s\n", pad_file);
    st_logInfo("Dup filter file string : %s\n", dup_filter_file);
    st_logInfo("Ignore first row : %s\n", ignore_first_row ? "True" : "False");

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

    // Sort/filter/pad files
    stList *prefixes_to_filter_by = load_sort_file(filter_file);
    stList *prefixes_to_pad = load_sort_file(pad_file);
    stList *prefixes_to_sort_by = load_sort_file(sort_file);
    stList *prefixes_to_dup_filter = load_sort_file(dup_filter_file);

    // Parse the header
    bool run_length_encode_bases;
    Tag *tag = taf_read_header_2(li, &run_length_encode_bases);

    // Write the header
    taf_write_header(tag, output);
    tag_destruct(tag);

    // Write the alignment blocks
    Alignment *alignment, *p_alignment = NULL, *pp_alignment = NULL;
    while((alignment = taf_read_block(p_alignment, run_length_encode_bases, li)) != NULL) {
        process_alignment_block(pp_alignment, p_alignment, prefixes_to_filter_by, prefixes_to_pad,
                                prefixes_to_sort_by, prefixes_to_dup_filter, run_length_encode_bases, ignore_first_row, output);
        pp_alignment = p_alignment;
        p_alignment = alignment;
    }
    if(p_alignment) { // Write the final block
        process_alignment_block(pp_alignment, p_alignment, prefixes_to_filter_by, prefixes_to_pad,
                                prefixes_to_sort_by, prefixes_to_dup_filter, run_length_encode_bases, ignore_first_row, output);
        alignment_destruct(p_alignment, 1);
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


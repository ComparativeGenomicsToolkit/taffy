/*
 * maf_norm: Normalize a taf or maf alignment to remove unnecessary blocks
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include <getopt.h>
#include <time.h>

int64_t maximum_block_length_to_merge = 10;
int64_t maximum_gap_length = 50;

void usage() {
    fprintf(stderr, "maf_norm [options]\n");
    fprintf(stderr, "Normalize a maf format alignment to remove small blocks using the -m and -n options to determine what to merge \n");
    fprintf(stderr, "-i --inputFile : Input maf file to normalize. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output maf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-m --maximumBlockLengthToMerge : Only merge together any two adjacent blocks if one or both is less than this many bases long, by default: %" PRIi64 "\n", maximum_block_length_to_merge);
    fprintf(stderr, "-n --maximumGapLength : Only merge together two adjacent blocks if the total number of unaligned bases between the blocks is less than this many bases, by default: %" PRIi64 "\n", maximum_gap_length);
    fprintf(stderr, "-h --help : Print this help message\n");
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;
    bool run_length_encode_bases = 0;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "help", no_argument, 0, 'h' },
                                                { "maximumBlockLengthToMerge", required_argument, 0, 'm' },
                                                { "maximumGapLength", required_argument, 0, 'n' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:hm:n:", long_options, &option_index);
        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'i':
                inputFile = optarg;
                break;
            case 'o':
                outputFile = optarg;
                break;
            case 'h':
                usage();
                return 0;
            case 'm':
                maximum_block_length_to_merge = atol(optarg);
                break;
            case 'n':
                maximum_gap_length = atol(optarg);
                break;
            default:
                usage();
                return 1;
        }
    }

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Input file string : %s\n", inputFile);
    st_logInfo("Output file string : %s\n", outputFile);
    st_logInfo("Maximum block length to merge : %" PRIi64 "\n", maximum_block_length_to_merge);
    st_logInfo("Maximum gap length : %" PRIi64 "\n", maximum_gap_length);

    //////////////////////////////////////////////
    // Read in the taf blocks and merge blocks that are sufficiently small
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");

    // Pass the header line to determine parameters and write the updated taf header
    stList *tags = maf_read_header(input);
    assert(stList_length(tags) % 2 == 0);
    maf_write_header(tags, output);

    Alignment *alignment, *p_alignment = NULL;
    while((alignment = maf_read_block(input)) != NULL) {
        if(p_alignment != NULL) {
            // Try merging the blocks
            alignment_link_adjacent(p_alignment, alignment, 0); // Link to calculate the gap length
            if ((alignment_length(p_alignment) <= maximum_block_length_to_merge ||
                 alignment_length(alignment) <= maximum_block_length_to_merge) &&
                    alignment_total_gap_length(p_alignment) <= maximum_gap_length) {
                p_alignment = alignment_merge_adjacent(p_alignment, alignment);
            } else {
                maf_write_block(p_alignment, output); // Write the maf block
                alignment_destruct(p_alignment); // clean up the left-most block
                p_alignment = alignment;
            }
        }
        else {
            p_alignment = alignment;
        }
    }
    if(p_alignment != NULL) {
        maf_write_block(p_alignment, output); // Write the last maf/taf block
        alignment_destruct(p_alignment);
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    if(inputFile != NULL) {
        fclose(input);
    }
    if(outputFile != NULL) {
        fclose(output);
    }

    st_logInfo("taf_normalize is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}


/*
 * maf_norm: Normalize a taf or maf alignment to remove unnecessary blocks
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "sonLib.h"
#include <getopt.h>
#include <time.h>

int64_t maximum_block_length_to_merge = 200;
int64_t maximum_gap_length = 30;
float fraction_shared_rows = 0.6;
int64_t repeat_coordinates_every_n_columns = 1000;

static void usage() {
    fprintf(stderr, "taf_norm [options]\n");
    fprintf(stderr, "Normalize a taf format alignment to remove small blocks using the -m and -n options to determine what to merge \n");
    fprintf(stderr, "-i --inputFile : Input taf file to normalize. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output taf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-k --maf : Print maf output instead of taf\n");
    fprintf(stderr, "-m --maximumBlockLengthToMerge : Only merge together any two adjacent blocks if one or both is less than this many bases long, by default: %" PRIi64 "\n", maximum_block_length_to_merge);
    fprintf(stderr, "-n --maximumGapLength : Only merge together two adjacent blocks if the total number of unaligned bases between the blocks is less than this many bases, by default: %" PRIi64 "\n", maximum_gap_length);
    fprintf(stderr, "-q --fractionSharedRows : The fraction of rows between two blocks that need to be shared for a merge, default: %f\n", fraction_shared_rows);
    fprintf(stderr, "-s --repeatCoordinatesEveryNColumns : Repeat coordinates of each sequence at least every n columns. By default: %" PRIi64 "\n", repeat_coordinates_every_n_columns);
    fprintf(stderr, "-c --useCompression : Write the output using bgzip compression.\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

static Alignment *get_next_taf_block(LI *li, bool run_length_encode_bases) {
    static Alignment *alignments[3];
    static int64_t alignment_index=0;
    assert(alignment_index >= 0);
    while(alignment_index < 3) {
        alignments[alignment_index] = taf_read_block(alignment_index == 0 ? NULL : alignments[alignment_index-1],
                                                     run_length_encode_bases, li); // Read a block
        if(alignments[alignment_index] == NULL) { // The read block is empty
            break;
        }
        alignment_index++;
    }
    assert(alignment_index >= 0);
    if(alignment_index == 0) {
        return NULL;
    }
    Alignment *block = alignments[0];
    // Shift down the remaining alignments
    alignments[0] = alignments[1];
    alignments[1] = alignments[2];
    alignments[2] = NULL;
    // Reduce the alignment index
    alignment_index--;
    return block;
}

int taf_norm_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;
    bool run_length_encode_bases = 0;
    bool output_maf = 0;
    bool use_compression = 0;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "maf", no_argument, 0, 'k' },
                                                { "help", no_argument, 0, 'h' },
                                                { "maximumBlockLengthToMerge", required_argument, 0, 'm' },
                                                { "maximumGapLength", required_argument, 0, 'n' },
                                                { "fractionSharedRows", required_argument, 0, 'q' },
                                                { "repeatCoordinatesEveryNColumns", required_argument, 0, 's' },
                                                { "useCompression", no_argument, 0, 'c' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:hcm:n:kq:s:", long_options, &option_index);
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
            case 'k':
                output_maf = 1;
                break;
            case 'm':
                maximum_block_length_to_merge = atol(optarg);
                break;
            case 'n':
                maximum_gap_length = atol(optarg);
                break;
            case 'q':
                fraction_shared_rows = atof(optarg);
                break;
            case 'c':
                use_compression = 1;
                break;
            case 's':
                repeat_coordinates_every_n_columns = atol(optarg);
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
    st_logInfo("Output maf : %s\n", output_maf ? "true" : "false");
    st_logInfo("Repeat coordinates every n bases : %" PRIi64 "\n", repeat_coordinates_every_n_columns);
    st_logInfo("Fraction shared rows to merge adjacent blocks : %f\n", fraction_shared_rows);
    st_logInfo("Write compressed output : %s\n", use_compression ? "true" : "false");

    //////////////////////////////////////////////
    // Read in the taf blocks and merge blocks that are sufficiently small
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    LW *output = LW_construct(outputFile == NULL ? stdout : fopen(outputFile, "w"), use_compression);
    LI *li = LI_construct(input);

    // Pass the header line to determine parameters and write the updated taf header
    Tag *tag = taf_read_header(li);
    Tag *t = tag_find(tag, "run_length_encode_bases");
    if(t != NULL && strcmp(t->value, "1") == 0) {
        run_length_encode_bases = 1;
        if(output_maf) { // Remove this tag from the maf output as not relevant
            tag = tag_remove(tag, "run_length_encode_bases");
        }
    }
    output_maf ? maf_write_header(tag, output) : taf_write_header(tag, output);
    tag_destruct(tag);

    Alignment *alignment, *p_alignment = NULL, *p_p_alignment = NULL;
    while((alignment = get_next_taf_block(li, run_length_encode_bases)) != NULL) {
        if(p_alignment != NULL) {
            int64_t common_rows = alignment_number_of_common_rows(p_alignment, alignment);
            int64_t total_rows = alignment->row_number + p_alignment->row_number - common_rows;
            if (common_rows >= total_rows * fraction_shared_rows &&
                (alignment_length(p_alignment) <= maximum_block_length_to_merge ||
                 alignment_length(alignment) <= maximum_block_length_to_merge) &&
                 alignment_total_gap_length(p_alignment) <= maximum_gap_length) {
                p_alignment = alignment_merge_adjacent(p_alignment, alignment);
            } else {
                output_maf ? maf_write_block(p_alignment, output) : taf_write_block(p_p_alignment, p_alignment, run_length_encode_bases, repeat_coordinates_every_n_columns, output); // Write the maf block
                if(p_p_alignment != NULL) {
                    alignment_destruct(p_p_alignment, 1); // Clean up the left-most block
                }
                p_p_alignment = p_alignment;
                p_alignment = alignment;
            }
        }
        else {
            p_alignment = alignment;
        }
    }
    if(p_alignment != NULL) {
        output_maf ? maf_write_block(p_alignment, output) : taf_write_block(p_p_alignment, p_alignment, run_length_encode_bases, -1, output); // Write the last taf block
        alignment_destruct(p_alignment, 1);
        if(p_p_alignment != NULL) {
            alignment_destruct(p_p_alignment, 1);
        }
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    if(inputFile != NULL) {
        fclose(input);
    }
    LW_destruct(output, outputFile != NULL);

    st_logInfo("taf_norm is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}


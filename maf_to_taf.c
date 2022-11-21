/*
 * maf_to_taf: Convert a maf alignment into a taf alignment
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "sonLib.h"
#include <getopt.h>
#include <time.h>

int64_t repeat_coordinates_every_n_columns = 1000;

void usage() {
    fprintf(stderr, "maf_to_taf [options]\n");
    fprintf(stderr, "Convert a maf format alignment to taf format\n");
    fprintf(stderr, "-i --inputFile : Input maf file to invert. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output taf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-r --runLengthEncodeBases : Run length encode bases\n");
    fprintf(stderr, "-s --repeatCoordinatesEveryNColumns : Repeat coordinates of each sequence at least every n columns. By default: %" PRIi64 "\n", repeat_coordinates_every_n_columns);
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
    bool run_length_encode_bases=0;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "runLengthEncodeBases", no_argument, 0, 'r' },
                                                { "repeatCoordinatesEveryNColumns", required_argument, 0, 's' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:hrs:", long_options, &option_index);
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
            case 'r':
                run_length_encode_bases = 1;
                break;
            case 's':
                repeat_coordinates_every_n_columns = atol(optarg);
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
    st_logInfo("Input file string : %s\n", inputFile);
    st_logInfo("Output file string : %s\n", outputFile);
    st_logInfo("Run length encode bases : %s\n", run_length_encode_bases ? "True" : "False");
    st_logInfo("Repeat coordinates every n bases : %" PRIi64 "\n", repeat_coordinates_every_n_columns);

    //////////////////////////////////////////////
    // Read in the maf blocks and convert to sequence of taf columns
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");

    // Read the maf header and write the taf header
    Tag *tag = maf_read_header(input);
    if(run_length_encode_bases) {
        tag = tag_construct("run_length_encode_bases", "1", tag);
    }
    taf_write_header(tag, output);
    tag_destruct(tag);

    // Now read in the maf blocks and write the alignment
    Alignment *alignment, *p_alignment = NULL;
    while((alignment = maf_read_block(input)) != NULL) {
        if(p_alignment != NULL) {
            alignment_link_adjacent(p_alignment, alignment, 1);
        }
        taf_write_block(p_alignment, alignment, run_length_encode_bases, repeat_coordinates_every_n_columns, output);
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment);
        }
        p_alignment = alignment;
    }
    if(p_alignment != NULL) {
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

    st_logInfo("maf_to_taf is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}


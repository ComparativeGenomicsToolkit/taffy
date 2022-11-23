/*
 * taf_to_maf: Convert a taf alignment into a maf alignment
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "sonLib.h"
#include <getopt.h>
#include <time.h>

void usage() {
    fprintf(stderr, "taf_to_maf [options]\n");
    fprintf(stderr, "Convert a taf format alignment to maf format\n");
    fprintf(stderr, "-i --inputFile : Input taf file to invert. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output maf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
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
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:h", long_options, &option_index);
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

    //////////////////////////////////////////////
    // Read in the taf blocks and convert to sequence of maf blocks
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");
    LI *li = LI_construct(input);

    // Pass the header line to determine parameters and write the maf header
    Tag *tag = taf_read_header(li);
    Tag *t = tag_find(tag, "run_length_encode_bases");
    if(t != NULL && strcmp(t->value, "1") == 0) {
        run_length_encode_bases = 1;
        tag = tag_remove(tag, "run_length_encode_bases");  // Remove this tag from the maf output as not relevant
    }
    maf_write_header(tag, output);
    tag_destruct(tag);

    Alignment *alignment, *p_alignment = NULL;
    while((alignment = taf_read_block(p_alignment, run_length_encode_bases, li)) != NULL) {
        maf_write_block(alignment, output);
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment, 1);
        }
        p_alignment = alignment;
    }
    if(p_alignment != NULL) {
        alignment_destruct(p_alignment, 1);
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

    st_logInfo("taf_to_maf is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}


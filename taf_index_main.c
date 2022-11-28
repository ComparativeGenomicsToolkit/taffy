/*
 * taf_to_maf: Convert a taf alignment into a maf alignment
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "taf_index.h"
#include <getopt.h>
#include <time.h>

void usage() {
    fprintf(stderr, "taf_index [options]\n");
    fprintf(stderr, "Index a TAF file, output goes in <file>.tai\n");
    fprintf(stderr, "-i --inputFile : Input taf file to invert [REQUIRED]\n");
    fprintf(stderr, "-b --blockSize : Write an index line for intervals of this many bp [default:1000000]\n");
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
    int64_t blockSize = 1000000;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "blockSize", required_argument, 0, 'b' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:b:h", long_options, &option_index);
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
            case 'b':
                blockSize = atoi(optarg);
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

    //////////////////////////////////////////////
    // Read in the taf blocks and convert to sequence of maf blocks
    //////////////////////////////////////////////

    if (inputFile == NULL) {
        fprintf(stderr, "Input file must be specified with -i\n");
        return 1;        
    }
    FILE *input = fopen(inputFile, "r");
    char *indexFile = (char*)st_calloc(strlen(inputFile) + 5, sizeof(char));
    sprintf(indexFile, "%s.tai", inputFile);
    FILE *output = fopen(indexFile, "w");    
    LI *li = LI_construct(input);

    index_taf(li, output, blockSize);

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    if(inputFile != NULL) {
        fclose(input);
    }
    if(output != NULL) {
        fclose(output);
    }
    free(indexFile);
    
    st_logInfo("taf_index is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}



/*
 * taf view: MAF/TAF conversion and subregion extraction
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "tai.h"
#include "sonLib.h"
#include <getopt.h>
#include <time.h>

static void usage() {
    fprintf(stderr, "taffy stats [options]\n");
    fprintf(stderr, "Print statistics from a TAF or MAF file\n");
    fprintf(stderr, "-i --inputFile : Input TAF or MAF file. If not specified reads from stdin\n");
    fprintf(stderr, "-s --sequenceLengths : Print length of each *reference* sequence in the (indexed) alignment\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int taf_stats_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *taf_fn = NULL;
    bool seq_lengths = false;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "sequenceLengths", no_argument, 0, 's' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:sh", long_options, &option_index);
        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'i':
                taf_fn = optarg;
                break;
            case 's':
                seq_lengths = 1;
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
    st_logInfo("Input file string : %s\n", taf_fn);

    //////////////////////////////////////////////
    // Do the stats
    //////////////////////////////////////////////

    if (seq_lengths == false) {
        fprintf(stderr, "Please pick a stats option from { -s }\n");
        return 1;
    }
    
    // load the input
    FILE *taf_fh = taf_fn == NULL ? stdin : fopen(taf_fn, "r");
    if (taf_fh == NULL) {
        fprintf(stderr, "Unable to open input TAF file: %s\n", taf_fn);
        return 1;
    }
    LI *li = LI_construct(taf_fh);
    
    // sniff the format
    int input_format = check_input_format(LI_peek_at_next_line(li));
    if (input_format == 2) {
        fprintf(stderr, "Input not supported: unable to detect ##maf or #taf header\n");
        return 1;
    }

    // load the index if it's required by the given options
    bool index_required = seq_lengths;
    char *tai_fn = NULL;
    FILE *tai_fh = NULL;
    Tai *tai = NULL;
    if (index_required) {
        tai_fn = tai_path(taf_fn);
        tai_fh = fopen(tai_fn, "r");
        if (tai_fh == NULL) {
            fprintf(stderr, "Required index %s not found. Please run taffy index first\n", tai_fn);
            return 1;
        }
        tai = tai_load(tai_fh, input_format == 1);
    }

    // do the stats
    if (seq_lengths) {
        stHash *seq_to_len = tai_sequence_lengths(tai, li);
        stList *seq_names = stHash_getKeys(seq_to_len);
        for (int64_t i = 0; i < stList_length(seq_names); ++i) {
            void *hash_val = stHash_search(seq_to_len, stList_get(seq_names, i));
            fprintf(stdout, "%s\t%" PRIi64 "\n", (char*)stList_get(seq_names, i), (int64_t)hash_val);
        }
        stHash_destruct(seq_to_len);
        stList_destruct(seq_names);
    }

        
    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    if (index_required) {
        free(tai_fn);
        tai_destruct(tai);
    }
    
    LI_destruct(li);
    if(taf_fn != NULL) {
        fclose(taf_fh);
    }

    st_logInfo("taffy stats is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}


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
    fprintf(stderr, "taf_find [options]\n");
    fprintf(stderr, "Query a TAF file, output is TAF\n");
    fprintf(stderr, "-i --inputFile : Input taf file to invert. Must be indexed with taf index [REQUIRED]\n");
    fprintf(stderr, "-r --region  : Region to query in CONTIG:START-END format (coordinates are 1-based half-open like samtools)\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *taf_fn = NULL;
    char *region = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "region", required_argument, 0, 'r' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:r:h", long_options, &option_index);
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
            case 'r':
                region = optarg;
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
    // Read in the taf blocks and convert to sequence of maf blocks
    //////////////////////////////////////////////

    if (taf_fn == NULL) {
        fprintf(stderr, "Input file must be specified with -i\n");
        return 1;        
    }
    FILE *taf_fh = fopen(taf_fn, "r");
    LI *li = LI_construct(taf_fh);    
    char *tai_fn = tai_path(taf_fn);
    FILE *tai_fh = fopen(tai_fn, "r");
    if (tai_fh == NULL) {
        fprintf(stderr, "Index %s not found. Please run taf index first\n", tai_fn);
        return 1;
    }

    Tai* tai = tai_load(tai_fh);

    TaiIt *tai_it = tai_iterator(tai, li, region);
    Alignment *alignment = NULL;
    Alignment *p_alignment = NULL;

    while ((alignment = tai_next(tai_it, li)) != NULL) {
        fprintf(stderr, "WRITE BLOCK\n");
        if (p_alignment) {
            for (Alignment_Row *row = p_alignment->row; row; row = row->n_row) {
                fprintf(stderr, "p_row %s %ld %ld %ld\n", row->sequence_name, row->start, row->length, (int64_t)row->r_row);
            }
        }
        if (alignment) {
            for (Alignment_Row *row = alignment->row; row; row = row->n_row) {
                fprintf(stderr, "row %s %ld %ld %ld\n", row->sequence_name, row->start, row->length, (int64_t)row->l_row);
            }
        }            
        
        taf_write_block(p_alignment, alignment, false, 1000, stdout);

        if (p_alignment) {
            alignment_destruct(p_alignment);
        }
        p_alignment = alignment;
    }
    if (p_alignment) {
        alignment_destruct(p_alignment);
    }
        

    if (tai_it) fprintf(stderr, "fam");
    
    tai_destruct(tai);

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    if(taf_fh != NULL) {
        fclose(taf_fh);
    }
    if(tai_fh != NULL) {
        fclose(tai_fh);
    }
    free(tai_fn);

    LI_destruct(li);
    
    st_logInfo("taf_find is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}



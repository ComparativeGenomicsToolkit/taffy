/*
 * taf_find: Pull a subrange out of a taf.  Can only query "first row" sequences. 
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "tai.h"
#include <getopt.h>
#include <time.h>

static int64_t repeat_coordinates_every_n_columns = 1000;

void usage() {
    fprintf(stderr, "taf_find [options]\n");
    fprintf(stderr, "Query a TAF file, output is TAF\n");
    fprintf(stderr, "-i --inputFile : Input taf file to invert. Must be indexed with taf index [REQUIRED]\n");
    fprintf(stderr, "-o --outputFile : Output taf file. If not specified outputs to stdout\n");    
    fprintf(stderr, "-r --region  : Print only CONTIG:START-END, where CONTIG is a row-0 sequence name, and START-END are 0-based open-ended like BED\n");
    fprintf(stderr, "-s --repeatCoordinatesEveryNColumns : Repeat coordinates of each sequence at least every n columns. By default: %" PRIi64 "\n", repeat_coordinates_every_n_columns);    
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
    char *out_fn = NULL;
    char *region = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },                                                
                                                { "region", required_argument, 0, 'r' },
                                                { "repeatCoordinatesEveryNColumns", required_argument, 0, 's' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:r:s:h", long_options, &option_index);
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
            case 'o':
                out_fn = optarg;
                break;                
            case 'r':
                region = optarg;
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
    st_logInfo("Input file string : %s\n", taf_fn);
    if (out_fn) {
        st_logInfo("Output file string : %s\n", out_fn);
    }
    st_logInfo("Repeat coordinates every n bases : %" PRIi64 "\n", repeat_coordinates_every_n_columns);

    //////////////////////////////////////////////
    //Print the TAF blocks, using the .tai if a region is specified
    //////////////////////////////////////////////

    if (taf_fn == NULL) {
        fprintf(stderr, "Input file must be specified with -i\n");
        return 1;        
    }
    FILE *taf_fh = fopen(taf_fn, "r");
    if (taf_fh == NULL) {
        fprintf(stderr, "Unable to open input TAF file: %s\n", taf_fn);
        return 1;
    }
    LI *li = LI_construct(taf_fh);
    if (!LI_indexable(li)) {
        fprintf(stderr, "Input TAF file must be either uncompressed or bgzipped: gzip not supported: %s\n", taf_fn);
        return 1;
    }

    FILE *out_fh = out_fn == NULL ? stdout : fopen(out_fn, "w");
    if (out_fh == NULL) {
        fprintf(stderr, "Unable to ouput output TAF file: %s\n", out_fn);
        return 1;
    }

    // Pass the header line to determine parameters and write the updated taf header
    bool run_length_encode_bases = 0;
    Tag *tag = taf_read_header(li);
    Tag *t = tag_find(tag, "run_length_encode_bases");
    if(t != NULL && strcmp(t->value, "1") == 0) {
        run_length_encode_bases = 1;
    }
    taf_write_header(tag, out_fh);
    tag_destruct(tag);
    
    if (region) {
        int64_t region_start;
        int64_t region_length;
        char *region_seq = tai_parse_region(region, &region_start, &region_length);
        if (region_seq == NULL) {
            fprintf(stderr, "Invalid region: %s\n", region);
            return 1;
        }
        st_logInfo("Region: contig=%s start=%" PRIi64 " length=%" PRIi64 "\n", region_seq, region_start, region_length);
        
        char *tai_fn = tai_path(taf_fn);        
        FILE *tai_fh = fopen(tai_fn, "r");
        
        if (tai_fh == NULL) {
            fprintf(stderr, "Index %s not found. Please run taf index first\n", tai_fn);
            return 1;
        }

        Tai* tai = tai_load(tai_fh);

        TaiIt *tai_it = tai_iterator(tai, li, run_length_encode_bases, region_seq, region_start, region_length);
        if (tai_it == NULL) {
            fprintf(stderr, "Region %s not found in taf index\n", region);
            return 1;
        }
        Alignment *alignment = NULL;
        Alignment *p_alignment = NULL;

        while ((alignment = tai_next(tai_it, li)) != NULL) {
            taf_write_block(p_alignment, alignment, run_length_encode_bases, repeat_coordinates_every_n_columns, out_fh);
            
            if (p_alignment) {
                alignment_destruct(p_alignment, true);
            }
            p_alignment = alignment;
        }
        if (p_alignment) {
            alignment_destruct(p_alignment, true);
        }

        tai_destruct(tai);

        if(tai_fh != NULL) {
            fclose(tai_fh);
        }
        free(tai_fn);        
    } else {
        Alignment *alignment = NULL;
        Alignment *p_alignment = NULL;
        
        while((alignment = taf_read_block(p_alignment, run_length_encode_bases, li)) != NULL) {
            taf_write_block(p_alignment, alignment, run_length_encode_bases, repeat_coordinates_every_n_columns, out_fh);
            if (p_alignment) {
                alignment_destruct(p_alignment, true);
            }
            p_alignment = alignment;
        }
        if (p_alignment) {
            alignment_destruct(p_alignment, true);
        }            
    }
        
    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    if(taf_fh != NULL) {
        fclose(taf_fh);
    }

    LI_destruct(li);
    
    st_logInfo("taf_find is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}



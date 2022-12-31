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

static int64_t repeat_coordinates_every_n_columns = 10000;

static void usage() {
    fprintf(stderr, "taf view [options]\n");
    fprintf(stderr, "Convert between TAF and MAF formats\n");
    fprintf(stderr, "-i --inputFile : Input TAF or MAF file to convert. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output file. If not specified outputs to stdout\n");
    fprintf(stderr, "-m --maf : Output in MAF format [default=TAF format]\n");
    fprintf(stderr, "-r --region  : Print only SEQ:START-END, where SEQ is a row-0 sequence name, and START-END are 0-based open-ended like BED\n");
    fprintf(stderr, "-s --repeatCoordinatesEveryNColumns : Repeat TAF coordinates of each sequence at least every n columns. By default: %" PRIi64 "\n", repeat_coordinates_every_n_columns);
    fprintf(stderr, "-u --runLengthEncodeBases : Run length encode bases in TAF\n");
    fprintf(stderr, "-c --useCompression : Write the output using bgzip compression.\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

// sniff format
// 0: taf
// 1: maf
// 2: unknown
static int check_input_format(const char *line) {
    int ret = 2;
    assert(line != NULL);
    stList *tokens = stString_split(line);
    if (stList_length(tokens) > 0) {
        if (strcmp(stList_get(tokens, 0), "#taf") == 0) {
            ret = 0;
        } else if (strcmp(stList_get(tokens, 0), "##maf") == 0) {
            ret = 1;
        }
    }            
    stList_destruct(tokens);
    return ret;
}

int taf_view_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;
    bool run_length_encode_bases = 0;
    bool maf_output = false;
    char *region = NULL;
    bool use_compression = 0;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "maf", no_argument, 0, 'm' },                                                
                                                { "runLengthEncodeBases", no_argument, 0, 'u' },
                                                { "repeatCoordinatesEveryNColumns", required_argument, 0, 's' },
                                                { "region", required_argument, 0, 'r' },
                                                { "useCompression", no_argument, 0, 'c' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:mucs:r:h", long_options, &option_index);
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
            case 'm':
                maf_output = 1;
                break;
            case 'u':
                run_length_encode_bases = 1;
                break;
            case 's':
                repeat_coordinates_every_n_columns = atol(optarg);
                break;
            case 'r':
                region = optarg;
                break;
            case 'c':
                use_compression = 1;
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
    st_logInfo("Write compressed output : %s\n", use_compression ? "true" : "false");

    //////////////////////////////////////////////
    // Read in the taf/maf blocks and convert to sequence of taf/maf blocks
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    LW *output = LW_construct(outputFile == NULL ? stdout : fopen(outputFile, "w"), use_compression);
    LI *li = LI_construct(input);

    // sniff the format
    int input_format = check_input_format(LI_peek_at_next_line(li));
    if (input_format == 2) {
        fprintf(stderr, "Input not supported: unable to detect ##maf or #taf header\n");
        return 1;
    }
    bool maf_input = input_format == 1;
    bool taf_input = !maf_input;
    bool taf_output = !maf_output;

    if (maf_input == true && region != NULL) {
        fprintf(stderr, "-r only supported on (indexed) TAF input\n");
        return 1;
    }

    // Parse the header
    Tag *tag = maf_input ? maf_read_header(li) : taf_read_header(li);
    if (maf_input && !maf_output) {
        if(run_length_encode_bases) {
            tag = tag_construct("run_length_encode_bases", "1", tag);
        }
    }
    else if (taf_input) {
        Tag *t = tag_find(tag, "run_length_encode_bases");
        if(t != NULL && strcmp(t->value, "1") == 0) {
            run_length_encode_bases = 1;
            if (maf_output) {
                tag = tag_remove(tag, "run_length_encode_bases");  // Remove this tag from the maf output as not relevant
            }
        }
    }
    if (maf_output) {
        maf_write_header(tag, output);
    } else {
        taf_write_header(tag, output);
    }
    tag_destruct(tag);

    if (taf_input == true) {
        if (region) {
            int64_t region_start;
            int64_t region_length;
            char *region_seq = tai_parse_region(region, &region_start, &region_length);
            if (region_seq == NULL) {
                fprintf(stderr, "Invalid region: %s\n", region);
                return 1;
            }
            st_logInfo("Region: contig=%s start=%" PRIi64 " length=%" PRIi64 "\n", region_seq, region_start, region_length);
        
            char *tai_fn = tai_path(inputFile);        
            FILE *tai_fh = fopen(tai_fn, "r");
        
            if (tai_fh == NULL) {
                fprintf(stderr, "Index %s not found. Please run taffy index first\n", tai_fn);
                return 1;
            }

            Tai* tai = tai_load(tai_fh);

            TaiIt *tai_it = tai_iterator(tai, li, run_length_encode_bases, region_seq, region_start, region_length);
            if (tai_it == NULL) {
                fprintf(stderr, "Region %s not found in taffy index\n", region);
                return 1;
            }
            Alignment *alignment = NULL;
            Alignment *p_alignment = NULL;

            while ((alignment = tai_next(tai_it, li)) != NULL) {
                if (taf_output) {
                    taf_write_block(p_alignment, alignment, run_length_encode_bases, repeat_coordinates_every_n_columns, output);
                } else {
                    maf_write_block(alignment, output);
                }            
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
                if (taf_output) {
                    taf_write_block(p_alignment, alignment, run_length_encode_bases, repeat_coordinates_every_n_columns, output);
                } else {
                    maf_write_block(alignment, output);
                }
                if (p_alignment) {
                    alignment_destruct(p_alignment, true);
                }
                p_alignment = alignment;
            }
            if (p_alignment) {
                alignment_destruct(p_alignment, true);
            }            
        }
    } else {
        assert(maf_input == true);

        // Now read in the maf blocks and write the alignment
        Alignment *alignment, *p_alignment = NULL;
        while((alignment = maf_read_block(li)) != NULL) {
            if(p_alignment != NULL) {
                alignment_link_adjacent(p_alignment, alignment, 1);
            }
            taf_write_block(p_alignment, alignment, run_length_encode_bases, repeat_coordinates_every_n_columns, output);
            if(p_alignment != NULL) {
                alignment_destruct(p_alignment, 1);
            }
            p_alignment = alignment;
        }
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment, 1);
        }
    }
    
    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    LI_destruct(li);
    if(inputFile != NULL) {
        fclose(input);
    }
    LW_destruct(output, outputFile != NULL);

    st_logInfo("taf view is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}


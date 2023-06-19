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
    fprintf(stderr, "taffy view [options]\n");
    fprintf(stderr, "Convert between TAF and MAF formats\n");
    fprintf(stderr, "-i --inputFile : Input TAF or MAF file to convert. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output file. If not specified outputs to stdout\n");
    fprintf(stderr, "-m --maf : Output in MAF format [default=TAF format]\n");
    fprintf(stderr, "-p --paf : Output in all-to-one PAF format [default=TAF format]\n");
    fprintf(stderr, "-a --all-paf : Output in all-to-all PAF format [default=TAF format]\n");
    fprintf(stderr, "-r --region  : Print only SEQ:START-END, where SEQ is a row-0 sequence name, and START-END are 0-based open-ended like BED\n");
    fprintf(stderr, "-s --repeatCoordinatesEveryNColumns : Repeat TAF coordinates of each sequence at least every n columns. By default: %" PRIi64 "\n", repeat_coordinates_every_n_columns);
    fprintf(stderr, "-u --runLengthEncodeBases : Run length encode bases in TAF\n");
    fprintf(stderr, "-c --useCompression : Write the output using bgzip compression.\n");
    fprintf(stderr, "-n --nameMapFile : Apply the given two-column tab-separated name mapping to all assembly names in alignment\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
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
    bool paf_output = false;
    bool all_to_all_paf = false;
    char *region = NULL;
    bool use_compression = 0;
    char *nameMapFile = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "maf", no_argument, 0, 'm' },
                                                { "paf", no_argument, 0, 'p' },
                                                { "all-paf", no_argument, 0, 'a' },
                                                { "runLengthEncodeBases", no_argument, 0, 'u' },
                                                { "repeatCoordinatesEveryNColumns", required_argument, 0, 's' },
                                                { "region", required_argument, 0, 'r' },
                                                { "useCompression", no_argument, 0, 'c' },
                                                { "nameMapFile", required_argument, 0, 'n' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:mpaucs:r:n:h", long_options, &option_index);
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
            case 'p':
                paf_output = 1;
                break;
            case 'a':
                paf_output = 1;
                all_to_all_paf = 1;
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
            case 'n':
                nameMapFile = optarg;
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
    if (nameMapFile) {
        st_logInfo("Name map file string : %s\n", nameMapFile);
    }

    //////////////////////////////////////////////
    // Read in the taf/maf blocks and convert to sequence of taf/maf blocks
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    if (input == NULL) {
        fprintf(stderr, "Unable to open input file: %s\n", inputFile);
        return 1;
    }

    FILE *output_fh = outputFile == NULL ? stdout : fopen(outputFile, "w");
    if (output_fh == NULL) {
        fprintf(stderr, "Unable to open output file: %s\n", outputFile);
        return 1;
    }

    stHash *genome_name_map = NULL;
    if (nameMapFile != NULL) {
        genome_name_map = load_genome_name_mapping(nameMapFile);
    }
    
    LW *output = LW_construct(output_fh, use_compression);
    LI *li = LI_construct(input);

    // sniff the format
    int input_format = check_input_format(LI_peek_at_next_line(li));
    if (input_format == 2) {
        fprintf(stderr, "Input not supported: unable to detect ##maf or #taf header\n");
        return 1;
    }
    bool maf_input = input_format == 1;
    bool taf_input = !maf_input;
    bool taf_output = !maf_output && !paf_output;

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
    } else if (taf_output) {
        taf_write_header(tag, output);
    }
    tag_destruct(tag);

    // three cases below:
    // 1) generic maf/taf index lookup if (region)
    // 2) scan whole taf
    // 3) scan whole maf    
    if (region) {
        int64_t region_start;
        int64_t region_length;
        char *region_seq = tai_parse_region(region, &region_start, &region_length);
        if (region_seq == NULL) {
            fprintf(stderr, "Invalid region: %s\n", region);
            return 1;
        }
        // apply the name mapping to the region
        char *mapped_region_seq = genome_name_map != NULL ? apply_genome_name_mapping(genome_name_map, region_seq) : NULL;
        if (mapped_region_seq != NULL) {
            free(region_seq);
            region_seq = mapped_region_seq;
        }
        
        st_logInfo("Region: contig=%s start=%" PRIi64 " length=%" PRIi64 "\n", region_seq, region_start, region_length);
        
        char *tai_fn = tai_path(inputFile);        
        FILE *tai_fh = fopen(tai_fn, "r");
        
        if (tai_fh == NULL) {
            fprintf(stderr, "Index %s not found. Please run taffy index first\n", tai_fn);
            return 1;
        }

        Tai* tai = tai_load(tai_fh, !taf_input);

        TaiIt *tai_it = tai_iterator(tai, li, run_length_encode_bases, region_seq, region_start, region_length);
        if (tai_it == NULL) {
            fprintf(stderr, "Region %s:%" PRIi64 "-%" PRIi64 " not found in taffy index\n", region_seq, region_start, region_length);
            return 1;
        }
        Alignment *alignment = NULL;
        Alignment *p_alignment = NULL;

        while ((alignment = tai_next(tai_it, li)) != NULL) {
            // apply the name mapping to the alignment block
            if (genome_name_map) {
                apply_genome_name_mapping_to_alignment(genome_name_map, alignment);
            }
            if (taf_output) {
                taf_write_block(p_alignment, alignment, run_length_encode_bases, repeat_coordinates_every_n_columns, output);
            } else if (maf_output) {
                maf_write_block(alignment, output);
            } else {
                assert(paf_output == true);
                paf_write_block(alignment, output, all_to_all_paf);
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
    } else if (taf_input) {
        Alignment *alignment = NULL;
        Alignment *p_alignment = NULL;        
        while((alignment = taf_read_block(p_alignment, run_length_encode_bases, li)) != NULL) {
            // apply the name mapping to the alignment block
            if (genome_name_map) {
                apply_genome_name_mapping_to_alignment(genome_name_map, alignment);
            }
            if (taf_output) {
                taf_write_block(p_alignment, alignment, run_length_encode_bases, repeat_coordinates_every_n_columns, output);
            } else if (maf_output) {
                maf_write_block(alignment, output);
            } else {
                assert(paf_output == true);
                paf_write_block(alignment, output, all_to_all_paf);
            }                    
            if (p_alignment) {
                alignment_destruct(p_alignment, true);
            }
            p_alignment = alignment;
        }
        if (p_alignment) {
            alignment_destruct(p_alignment, true);
        }            
    } else {
        assert(maf_input == true);
        Alignment *alignment, *p_alignment = NULL;
        while((alignment = maf_read_block(li)) != NULL) {
            // apply the name mapping to the alignment block
            if (genome_name_map) {
                apply_genome_name_mapping_to_alignment(genome_name_map, alignment);
            }            
            if(p_alignment != NULL) {
                alignment_link_adjacent(p_alignment, alignment, 1);
            }
            if (taf_output) {
                taf_write_block(p_alignment, alignment, run_length_encode_bases, repeat_coordinates_every_n_columns, output);
            } else if (maf_output) {
                maf_write_block(alignment, output);
            } else {
                assert(paf_output == true);
                paf_write_block(alignment, output, all_to_all_paf);
            }
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

    if (genome_name_map != NULL) {
        stHash_destruct(genome_name_map);
    }

    st_logInfo("taffy view is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}


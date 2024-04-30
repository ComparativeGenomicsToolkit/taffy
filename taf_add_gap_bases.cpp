/*
 * maf_add_gap_bases: Add in the gap substrings between taf blocks
 *
 *  Released under the MIT license, see LICENSE.txt
*/
extern "C" {
#include "taf.h"
#include "sonLib.h"
}
#include "bioioC.h"
#include <getopt.h>
#include <time.h>
#ifdef USE_HAL
#include "halBlockViz.h"
#endif

static int64_t repeat_coordinates_every_n_columns = 10000;
static int64_t maximum_gap_string_length = 50;

static void usage() {
    fprintf(stderr, "taffy add_gap_bases SEQ_FILExN [options]\n");    
    fprintf(stderr, "Add interstitial gap strings to taf file\n");
    fprintf(stderr, "-i --inputFile : Input taf file to normalize. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output taf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-a --halFile : HAL file for extracting gap sequence (MAF must be created with hal2maf *without* --onlySequenceNames)\n");
    fprintf(stderr, "-m --maximumGapStringLength : The maximum size of a gap string to add, be default: %" PRIi64 "\n",
            maximum_gap_string_length);
    fprintf(stderr, "-s --repeatCoordinatesEveryNColumns : Repeat TAF coordinates of each sequence at least every n columns. By default: %" PRIi64 "\n", repeat_coordinates_every_n_columns);
    fprintf(stderr, "-c --useCompression : Write the output using bgzip compression.\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int taf_add_gap_bases_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;
    char *hal_file = NULL;
    bool run_length_encode_bases = 0;
    bool use_compression = 0;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "halFile", required_argument, 0, 'a' },
                                                { "repeatCoordinatesEveryNColumns", required_argument, 0, 's' },
                                                { "useCompression", no_argument, 0, 'c' },                                                
                                                { "help", no_argument, 0, 'h' },
                                                { "maximumGapStringLength", required_argument, 0, 'm' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:a:s:chm:", long_options, &option_index);
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
            case 'a':
                hal_file = optarg;
                break;
            case 's':
                repeat_coordinates_every_n_columns = atol(optarg);
                break;
            case 'c':
                use_compression = 1;
                break;                
            case 'h':
                usage();
                return 0;
            case 'm':
                maximum_gap_string_length = atol(optarg);
                break;
            default:
                usage();
                return 1;
        }
    }

    if ((hal_file == NULL) == (optind >= argc)) {
        fprintf(stderr, "[taf] Sequences must be specified either via fasta arguments OR the -a option (but not both)\n");
        return 1;
    }
#ifndef USE_HAL
    if (hal_file) {
        fprintf(stderr, "[taf] taf was not built with HAL support. Set HALDIR and recompile in order to use -a\n");
        return 1;
    }
#endif
    
    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Input file string : %s\n", inputFile);
    st_logInfo("Output file string : %s\n", outputFile);
    if (hal_file) {
        st_logInfo("HAL file string : %s\n", hal_file);
    } else {
        st_logInfo("Number of input FASTA files : %ld\n", argc - optind);
    }            
    st_logInfo("Maximum maximum gap string length : %" PRIi64 "\n", maximum_gap_string_length);

    //////////////////////////////////////////////
    // Read in the sequence files
    //////////////////////////////////////////////

    stHash *fastas = NULL;
    stSet *hal_species = NULL;
    int hal_handle = -1;
    if (optind < argc) {
        fastas = load_sequences_from_fasta_files(&(argv[optind]), argc - optind);
    }
    else {
        hal_species = load_sequences_from_hal_file(hal_file, &hal_handle);
    }

    //////////////////////////////////////////////
    // Read in the taf blocks, add the gap strings and output the taf blocks with added gap strings
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    LW *output = LW_construct(outputFile == NULL ? stdout : fopen(outputFile, "w"), use_compression);
    LI *li = LI_construct(input);

    // Pass the header line to determine parameters and write the updated taf header
    Tag *tag = taf_read_header_2(li, &run_length_encode_bases);
    taf_write_header(tag, output);
    tag_destruct(tag);

    Alignment *alignment, *p_alignment = NULL;
    while((alignment = taf_read_block(p_alignment, run_length_encode_bases, li)) != NULL) {
        // Add in the gap strings if there is a previous block
        if(p_alignment != NULL) {
            alignment_add_gap_strings(p_alignment, alignment, fastas, hal_handle, hal_species, maximum_gap_string_length);
        }

        // Write the block
        taf_write_block(p_alignment, alignment, run_length_encode_bases, repeat_coordinates_every_n_columns, output);

        // Clean up the previous alignment
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment, 1);
        }
        p_alignment = alignment; // Update the previous alignment
    }
    if(p_alignment != NULL) { // Clean up the final alignment
        alignment_destruct(p_alignment, 1);
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    LI_destruct(li);
    if(inputFile != NULL) {
        fclose(input);
    }
    LW_destruct(output, outputFile != NULL);

    if (fastas) {
        stHash_destruct(fastas);
    }
    if (hal_species) {
        stSet_destruct(hal_species);
    }

    st_logInfo("taffy add-gap-bases is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}


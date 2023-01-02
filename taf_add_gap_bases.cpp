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

static int64_t maximum_gap_string_length = 50;
static int64_t repeat_coordinates_every_n_columns = 1000;

static void usage() {
    fprintf(stderr, "taf_add_gap_bases SEQ_FILExN [options]\n");    
    fprintf(stderr, "Add interstitial gap strings to taf file\n");
    fprintf(stderr, "-i --inputFile : Input taf file to normalize. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output taf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-a --halFile : HAL file for extracting gap sequence (MAF must be created with hal2maf *without* --onlySequenceNames)\n");
    fprintf(stderr, "-m --maximumGapStringLength : The maximum size of a gap string to add, be default: %" PRIi64 "\n",
            maximum_gap_string_length);
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-s --repeatCoordinatesEveryNColumns : Repeat coordinates of each sequence at least every n columns. By default: %" PRIi64 "\n", repeat_coordinates_every_n_columns);
    fprintf(stderr, "-c --useCompression : Write the output using bgzip compression.\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

static void add_to_hash(void *fastas, const char *fasta_header, const char *sequence, int64_t length) {
    if(stHash_search((stHash *)fastas, (void *)fasta_header) != NULL) {
        // c++ gives an angry warning if we try to send our string literal directly to st_errAbort, so we do this
        char msg[8192];
        snprintf(msg, 8192, "Found duplicate sequence header: %s\n", fasta_header);
        st_errAbort(msg);
    }
    stHash_insert((stHash *)fastas, stString_copy(fasta_header), stString_copy(sequence));
}

// it turns out just scanning for a "." doesn't work, even with the output of hal2maf.  The reason is
// that spcecies names can contain "." characters -- at least they do in the VGP alignment.
// so this function uses knowloedge of the genome names to greedily scan for a prefix ending in "."
// that corresponds to a genome name in the hal. the extracted genome name is returned if found
// (it must be freed)
static char *extract_genome_name(const char *sequence_name, stSet *hal_species) {
    const char *dot = NULL;
    int64_t offset = 0;
    const char *last = sequence_name + strlen(sequence_name) - 1;

    do {
        dot = strchr(sequence_name + offset, '.');
        if (dot != NULL && dot != last && dot != sequence_name) {
            char *species_name = stString_getSubString(sequence_name, 0, dot-sequence_name);
            if (stSet_search(hal_species, species_name) != NULL) {
                return species_name;
            } else if (dot != last) {
                free(species_name);
                offset += (dot-sequence_name) + 1;
            }
        }
    } while (dot != NULL);

    // c++ gives an angry warning if we try to send our string literal directly to st_errAbort, so we do this
    char msg[8192];
    snprintf(msg, 8192, "[taf] Error: Unable to find a . that splits %s so that the left side is a genome in the HAL\n", sequence_name);
    st_errAbort(msg);
    return NULL;
}

// get a dna interval either from the fastas hash file or from the hal_handle
// note the string returned needs to be freed
static char *get_sequence_fragment(const char* sequence_name, int64_t start, int64_t length, stHash *fastas, int hal_handle, stSet *hal_species) {
    char *fragment = NULL;
    if (fastas) {
        assert(hal_handle == -1);
        char *seq = (char*)stHash_search(fastas, (void*)sequence_name);
        if(seq != NULL) {
            fragment = stString_getSubString(seq, start, length);
        }
    } else {
#ifdef USE_HAL
        assert(fastas == NULL);
        char* species_name = extract_genome_name(sequence_name, hal_species);
        char* chrom_name = (char*)sequence_name + strlen(species_name) + 1;
        fragment = halGetDna(hal_handle, species_name, chrom_name, start, start + length, NULL);
        free(species_name);
#else
        assert(false);
#endif
    }

    return fragment;
}


static void add_gap_strings(Alignment *p_alignment, Alignment *alignment, stHash *fastas, int hal_handle, stSet *hal_species) {
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        if(row->l_row != NULL && alignment_row_is_predecessor(row->l_row, row)) {
            int64_t gap_length = row->start - (row->l_row->start + row->l_row->length);
            if(gap_length <= maximum_gap_string_length && row->left_gap_sequence == NULL) {
                char* seq_interval = NULL;
                int64_t i = row->l_row->start + row->l_row->length;
                assert(i >= 0 && i < row->sequence_length);
                if(row->strand) {
                    seq_interval = get_sequence_fragment(row->sequence_name, i, gap_length, fastas, hal_handle, hal_species);
                }
                else { // Case sequence is on the negative strand
                    assert(row->sequence_length - i - gap_length >= 0);
                    char *s = get_sequence_fragment(row->sequence_name,
                                                    row->sequence_length - i - gap_length, gap_length, fastas, hal_handle, hal_species);
                    if(s != NULL) {
                        seq_interval = stString_reverseComplementString(s);
                        free(s);
                    }
                }
                if(seq_interval == NULL) {
                    st_logDebug("[taf] Missing sequence for gap, seq name: %s, skipping!\n", row->sequence_name);
                }
                else {
                    row->left_gap_sequence = seq_interval;
                    assert(strlen(row->left_gap_sequence) == (size_t)gap_length);
                }
            }
        }
        row = row->n_row;
    }
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
                                                { "help", no_argument, 0, 'h' },
                                                { "maximumGapStringLength", required_argument, 0, 'm' },
                                                { "repeatCoordinatesEveryNColumns", required_argument, 0, 's' },
                                                { "useCompression", no_argument, 0, 'c' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:a:hm:s:c", long_options, &option_index);
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
            case 'h':
                usage();
                return 0;
            case 's':
                repeat_coordinates_every_n_columns = atol(optarg);
                break;
            case 'm':
                maximum_gap_string_length = atol(optarg);
                break;
            case 'c':
                use_compression = 1;
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
    st_logInfo("Repeat coordinates every n bases : %" PRIi64 "\n", repeat_coordinates_every_n_columns);
    st_logInfo("Write compressed output : %s\n", use_compression ? "true" : "false");

    //////////////////////////////////////////////
    // Read in the sequence files
    //////////////////////////////////////////////

    stHash *fastas = NULL;
    stSet *hal_species = NULL;
    int hal_handle = -1;
    if (optind < argc) {
        fastas = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
        while(optind < argc) {
            char *seq_file = argv[optind++];
            st_logInfo("Parsing sequence file : %s\n", seq_file);
            FILE *fh = fopen(seq_file, "r");
            fastaReadToFunction(fh, fastas, add_to_hash);
            fclose(fh);
        }
        st_logInfo("Finished parsing sequence files\n");
    } else {
#ifdef USE_HAL
        hal_handle = halOpen(hal_file, NULL);
        // make a table of every species in the hal. this will be used to try to parse
        // contig names with lots of dots in them where we don't know which dot separates genome from chrom...
        hal_species = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
        struct hal_species_t *species_list = halGetSpecies(hal_handle, NULL);
        for (struct hal_species_t *species = species_list; species != NULL; species = species->next) {
            stSet_insert(hal_species, stString_copy(species->name));
        }
        halFreeSpeciesList(species_list);
#endif
    }

    //////////////////////////////////////////////
    // Read in the taf blocks, add the gap strings and output the taf blocks with added gap strings
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    LW *output = LW_construct(outputFile == NULL ? stdout : fopen(outputFile, "w"), use_compression);
    LI *li = LI_construct(input);

    // Pass the header line to determine parameters and write the updated taf header
    Tag *tag = taf_read_header(li);
    Tag *t = tag_find(tag, (char*)"run_length_encode_bases");
    if(t != NULL && strcmp(t->value, "1") == 0) {
        run_length_encode_bases = 1;
    }
    taf_write_header(tag, output);
    tag_destruct(tag);

    Alignment *alignment, *p_alignment = NULL;
    while((alignment = taf_read_block(p_alignment, run_length_encode_bases, li)) != NULL) {
        // Add in the gap strings if there is a previous block
        if(p_alignment != NULL) {
            add_gap_strings(p_alignment, alignment, fastas, hal_handle, hal_species);
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


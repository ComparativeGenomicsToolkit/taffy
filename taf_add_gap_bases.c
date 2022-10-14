/*
 * maf_add_gap_bases: Add in the gap substrings between taf blocks
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "bioioC.h"
#include <getopt.h>
#include <time.h>
#ifdef USE_HAL
#include "halBlockViz.h"
#endif

int64_t maximum_gap_string_length = 50;

void usage() {
    fprintf(stderr, "taf_add_gap_bases [options]\n");
    fprintf(stderr, "Add interstitial gap strings to taf file\n");
    fprintf(stderr, "-i --inputFile : Input taf file to normalize. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output taf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-f --fastaFile : Fasta file for extracting gap sequence (multipel allowed)\n");
    fprintf(stderr, "-a --halFile : HAL file for extracting gap sequence (MAF must be created with hal2maf *without* --onlySequenceNames)\n");
    fprintf(stderr, "-m --maximumGapStringLength : The maximum size of a gap string to add, be default: %" PRIi64 "\n",
            maximum_gap_string_length);
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

void add_to_hash(void *fastas, const char *fasta_header, const char *sequence, int64_t length) {
    if(stHash_search((stHash *)fastas, (void *)fasta_header) != NULL) {
        st_errAbort("Found duplicate sequence header: %s\n", fasta_header);
    }
    stHash_insert((stHash *)fastas, stString_copy(fasta_header), stString_copy(sequence));
}

// get a dna interval either from the fastas hash file or from the hal_handle
// note the string returned needs to be freed
char *get_sequence_fragment(const char* sequence_name, int64_t start, int64_t length, stHash *fastas, int hal_handle) {
    char *fragment = NULL;
    if (fastas) {
        assert(hal_handle == -1);
        char *seq = stHash_search(fastas, (void*)sequence_name);
        if(seq != NULL) {
            fragment = stString_getSubString(seq, start, length);
        }
    } else {
#ifdef USE_HAL
        assert(fastas == NULL);
        char *dot = strchr(sequence_name, '.');
        if (dot == NULL || *(dot+1) == '\0' || dot == sequence_name) {
            st_errAbort("[taf] Error: could not parse MAF sequence %s into GENOME.CONTIG\n", sequence_name);
        }
        char* species_name = stString_getSubString(sequence_name, 0, dot-sequence_name);
        char* chrom_name = dot + 1;
        fragment = halGetDna(hal_handle, species_name, chrom_name, start, start + length, NULL);
        free(species_name);
#else
        assert(false);
#endif
    }

    return fragment;
}


void add_gap_strings(Alignment *p_alignment, Alignment *alignment, stHash *fastas, int hal_handle) {
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        if(row->l_row != NULL && alignment_row_is_predecessor(row->l_row, row)) {
            int64_t gap_length = row->start - (row->l_row->start + row->l_row->length);
            if(gap_length <= maximum_gap_string_length && row->left_gap_sequence == NULL) {
                char* seq_interval = get_sequence_fragment(row->sequence_name, row->l_row->start + row->l_row->length, gap_length, fastas, hal_handle);
                if(seq_interval == NULL) {
                    st_logDebug("[taf] Missing sequence for gap, seq name: %s, skipping!\n", row->sequence_name);
                }
                else {
                    row->left_gap_sequence = seq_interval;
                    assert(strlen(row->left_gap_sequence) == gap_length);
                }
            }
        }
        row = row->n_row;
    }
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;
    stList* fasta_files = stList_construct();
    char *hal_file = NULL;
    bool run_length_encode_bases = 0;
    bool output_maf = 0;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "fastaFile", required_argument, 0, 'f' },
                                                { "halFile", required_argument, 0, 'a' },
                                                { "maf", no_argument, 0, 'k' },
                                                { "help", no_argument, 0, 'h' },
                                                { "maximumGapStringLength", required_argument, 0, 'm' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:f:a:hm:k", long_options, &option_index);
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
            case 'f':
                stList_append(fasta_files, optarg);
                break;
            case 'a':
                hal_file = optarg;
                break;
            case 'h':
                usage();
                return 0;
            case 'k':
                output_maf = 1;
                break;
            case 'm':
                maximum_gap_string_length = atol(optarg);
                break;
            default:
                usage();
                return 1;
        }
    }

    if ((hal_file == NULL) == (stList_length(fasta_files) == 0)) {
        fprintf(stderr, "[taf] Either -s or -a (but not both) must be used to specifiy the input sequence\n");
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
        st_logInfo("Number of input FASTA files : %ld\n", stList_length(fasta_files));
    }            
    st_logInfo("Maximum maximum gap string length : %" PRIi64 "\n", maximum_gap_string_length);

    //////////////////////////////////////////////
    // Read in the sequence files
    //////////////////////////////////////////////
    stHash *fastas = NULL;
    int hal_handle = -1;
    if (stList_length(fasta_files) > 0) {
        fastas = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
        for (int64_t i = 0; i < stList_length(fasta_files); ++i) {
            char *seq_file = stList_get(fasta_files, i);
            st_logInfo("Parsing sequence file : %s\n", seq_file);
            FILE *fh = fopen(seq_file, "r");
            fastaReadToFunction(fh, fastas, add_to_hash);
            fclose(fh);            
        }
        st_logInfo("Finished parsing sequence files\n");
    } else {
#ifdef USE_HAL
        hal_handle = halOpen(hal_file, NULL);
#endif
    }

    //////////////////////////////////////////////
    // Read in the taf blocks, add the gap strings and output the taf blocks with added gap strings
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");
    LI *li = LI_construct(input);

    // Pass the header line to determine parameters and write the updated taf header
    stList *tags = taf_read_header(li);
    assert(stList_length(tags) % 2 == 0);
    for(int64_t i=0; i<stList_length(tags); i+=2) {
        char *key = stList_get(tags, i);
        char *value = stList_get(tags, i+1);
        if(strcmp(key, "run_length_encode_bases") == 0 && strcmp(value, "1") == 0) {
            run_length_encode_bases = 1;
        }
    }
    taf_write_header(tags, output);

    Alignment *alignment, *p_alignment = NULL;
    while((alignment = taf_read_block(p_alignment, run_length_encode_bases, li)) != NULL) {
        // Add in the gap strings if there is a previous block
        if(p_alignment != NULL) {
            add_gap_strings(p_alignment, alignment, fastas, hal_handle);
        }

        // Write the block
        taf_write_block(p_alignment, alignment, run_length_encode_bases, output);

        // Clean up the previous alignment
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment);
        }
        p_alignment = alignment; // Update the previous alignment
    }
    if(p_alignment != NULL) { // Clean up the final alignment
        alignment_destruct(p_alignment);
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    LI_destruct(li);
    if(inputFile != NULL) {
        fclose(input);
    }
    if(outputFile != NULL) {
        fclose(output);
    }

    stList_destruct(fasta_files);
    if (fastas) {
        stHash_destruct(fastas);
    }

    st_logInfo("taf_add_gap_bases is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}


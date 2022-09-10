/*
 * maf_add_gap_bases: Add in the gap substrings between taf blocks
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "bioioC.h"
#include <getopt.h>
#include <time.h>

int64_t maximum_gap_string_length = 50;

void usage() {
    fprintf(stderr, "taf_add_gap_bases SEQ_FILExN [options]\n");
    fprintf(stderr, "Add interstitial gap strings to taf file\n");
    fprintf(stderr, "-i --inputFile : Input taf file to normalize. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output taf file. If not specified outputs to stdout\n");
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

void add_gap_strings(Alignment *p_alignment, Alignment *alignment, stHash *fastas) {
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        if(row->l_row != NULL && alignment_row_is_predecessor(row->l_row, row)) {
            int64_t gap_length = row->start - (row->l_row->start + row->l_row->length);
            if(gap_length <= maximum_gap_string_length && row->left_gap_sequence == NULL) {
                char *seq = stHash_search(fastas, row->sequence_name);
                if(seq == NULL) {
                    st_logDebug("Missing sequence for gap, seq name: %s, skipping!\n", row->sequence_name);
                }
                else {
                    row->left_gap_sequence = stString_getSubString(seq, row->l_row->start + row->l_row->length, gap_length);
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
    bool run_length_encode_bases = 0;
    bool output_maf = 0;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "maf", no_argument, 0, 'k' },
                                                { "help", no_argument, 0, 'h' },
                                                { "maximumGapStringLength", required_argument, 0, 'm' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:hm:k", long_options, &option_index);
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

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Input file string : %s\n", inputFile);
    st_logInfo("Output file string : %s\n", outputFile);
    st_logInfo("Maximum maximum gap string length : %" PRIi64 "\n", maximum_gap_string_length);

    //////////////////////////////////////////////
    // Read in the sequence files
    //////////////////////////////////////////////

    stHash *fastas = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    while(optind < argc) {
        char *seq_file = argv[optind++];
        st_logInfo("Parsing sequence file : %s\n", seq_file);
        FILE *fh = fopen(seq_file, "r");
        fastaReadToFunction(fh, fastas, add_to_hash);
        fclose(fh);
    }
    st_logInfo("Finished parsing sequence files\n");

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
            add_gap_strings(p_alignment, alignment, fastas);
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

    st_logInfo("taf_add_gap_bases is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}


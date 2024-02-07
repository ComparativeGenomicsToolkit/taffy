/*
 * maf_norm: Normalize a taf or maf alignment to remove unnecessary blocks
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "sonLib.h"
#include <getopt.h>
#include <time.h>

int64_t maximum_block_length_to_merge = 200;
int64_t maximum_gap_length = 30;
int64_t minimum_shared_rows = 1;
float fraction_shared_rows = 0.0;
int64_t repeat_coordinates_every_n_columns = 1000;

static void usage(void) {
    fprintf(stderr, "taffy norm [options]\n");
    fprintf(stderr, "Normalize a taf format alignment to remove small blocks using the -m and -n options to determine what to merge \n");
    fprintf(stderr, "-i --inputFile : Input taf file to normalize. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output taf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-k --maf : Print maf output instead of taf\n");
    fprintf(stderr, "-m --maximumBlockLengthToMerge : Only merge together any two adjacent blocks if one or both is less than this many bases long, by default: %" PRIi64 "\n", maximum_block_length_to_merge);
    fprintf(stderr, "-n --maximumGapLength : Only merge together two adjacent blocks if the total number of unaligned bases between the blocks is less than this many bases, by default: %" PRIi64 "\n", maximum_gap_length);
    fprintf(stderr, "-Q --minimumSharedRows : The minimum number of rows between two blocks that need to be shared for a merge, default: %" PRIi64 "\n", minimum_shared_rows);
    fprintf(stderr, "-q --fractionSharedRows : The fraction of rows between two blocks that need to be shared for a merge, default: %f\n", fraction_shared_rows);
    fprintf(stderr, "-d --filterGapCausingDupes : Reduce the number of MAF blocks by filtering out rows that induce gaps > maximumGapLength. Rows are only filtered out if they are duplications (contig of same name appears elsewhere in block, or contig with same prefix up to \".\" appears in the same block).\n");
    fprintf(stderr, "-s --repeatCoordinatesEveryNColumns : Repeat coordinates of each sequence at least every n columns. By default: %" PRIi64 "\n", repeat_coordinates_every_n_columns);
    fprintf(stderr, "-c --useCompression : Write the output using bgzip compression.\n");
    fprintf(stderr, "-a --halFile : HAL file for extracting gap sequence (MAF must be created with hal2maf *without* --onlySequenceNames)\n");
    fprintf(stderr, "-b --seqFiles : Fasta files for extracting gap sequence. Do not specify both this option and --halFile\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

static Alignment *get_next_taf_block(LI *li, bool run_length_encode_bases) {
    static Alignment *alignments[3];
    static int64_t alignment_index=0;
    assert(alignment_index >= 0);
    while(alignment_index < 3) {
        alignments[alignment_index] = taf_read_block(alignment_index == 0 ? NULL : alignments[alignment_index-1],
                                                     run_length_encode_bases, li); // Read a block
        if(alignments[alignment_index] == NULL) { // The read block is empty
            break;
        }
        alignment_index++;
    }
    assert(alignment_index >= 0);
    if(alignment_index == 0) {
        return NULL;
    }
    Alignment *block = alignments[0];
    // Shift down the remaining alignments
    alignments[0] = alignments[1];
    alignments[1] = alignments[2];
    alignments[2] = NULL;
    // Reduce the alignment index
    alignment_index--;
    return block;
}

// prototype logic to try to reduce the gap delta by greedily filtering out
// dupes with the biggest gaps.  if it doesn't find enough dupes to remove
// to cover gap_delta, it returns false and does nothing.  otherwise it returns
// true and removes the rows. 
static bool greedy_prune_by_gap(Alignment *alignment, int64_t maximum_gap_length) {

    // map row ptr to sample name, using everything up to first "." of sequence name
    // (which is quite hacky but not sure there's a choice)
    stHash *row_to_sample_name = stHash_construct2(NULL, free);
    // hash sample name to number of rows with the sample
    stHash *sample_to_count = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, free);
    int64_t i = 0;
    for (Alignment_Row *row = alignment->row; row != NULL; row = row->n_row, ++i) {
        char *dot = strchr(row->sequence_name, '.');
        char *sample_name = NULL;
        if (dot != NULL) {
            sample_name = stString_getSubString(row->sequence_name, 0, dot - row->sequence_name);
        } else {
            sample_name = stString_copy(row->sequence_name);
        }
        stHash_insert(row_to_sample_name, row, sample_name);
        
        int64_t *count = stHash_search(sample_to_count, sample_name);
        if (count == NULL) {
            count = (int64_t*)malloc(sizeof(int64_t));
            *count = 1;
            stHash_insert(sample_to_count, sample_name, count);
        } else {
            ++(*count);
        }
    }

    // check that all rows with gaps > max_gap are dupes, ie we can merge a pruned alignment
    i = 0;
    bool can_prune = true;
    // all rows that are dupes exceeding gap length are stored here
    stList *to_prune = stList_construct();
    // we only filter dupes if there's at least one row that would not get filtered
    // remember when that happens here
    stSet *samples_passing_gap_once = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);    
    for (Alignment_Row *row = alignment->row; row != NULL && can_prune; row = row->n_row, ++i) {
        int64_t gap = 0;
        if (row->l_row != NULL && alignment_row_is_predecessor(row->l_row, row)) {
            gap = row->start - (row->l_row->start + row->l_row->length);
        }
        char *sample_name = (char*)stHash_search(row_to_sample_name, row); 
        if (gap > maximum_gap_length) {            
            int64_t *count = stHash_search(sample_to_count, sample_name);            
            if (*count > 1 && row != alignment->row) {
                stList_append(to_prune, row);
            } else {
                can_prune = false;
            }
        } else {
            stSet_insert(samples_passing_gap_once, sample_name);    
        }
    }
    // if there isn't at least one copy that we can merge without dropping, then we
    // don't bother as we don't want to drop coverage below 1 for any sample    
    int64_t no_to_prune = stList_length(to_prune);
    for (i = 0; can_prune && i < no_to_prune; ++i) {
        Alignment_Row *row_to_prune = stList_get(to_prune, i);
        char *sample_name = stHash_search(row_to_sample_name, row_to_prune);
        assert(sample_name);        
        if (stSet_search(samples_passing_gap_once, sample_name) == NULL) {
            can_prune = false;
        }
    }

    bool pruned = false;
    if (can_prune) {
        // remove the rows
        assert(stList_length(to_prune) > 0);
        Alignment_Row *p_row = NULL;
        int64_t to_prune_idx = 0;
        i = 0;
        Alignment_Row *row_to_prune = stList_get(to_prune, to_prune_idx);
        Alignment_Row *row = alignment->row;
        while (row) {
            Alignment_Row *n_row = row->n_row;
            if (row == row_to_prune) {
                assert(p_row != NULL);
                p_row->n_row = row->n_row;
                alignment_row_destruct(row);
                ++to_prune_idx;
                row_to_prune = to_prune_idx < stList_length(to_prune) ? stList_get(to_prune, to_prune_idx) : NULL;
                --alignment->row_number;
                pruned = true;
            } else {
                p_row = row;
            }
            row = n_row;
            ++i;
        }
    }

    stHash_destruct(row_to_sample_name);
    stHash_destruct(sample_to_count);
    stList_destruct(to_prune);
    stSet_destruct(samples_passing_gap_once);
    
    return pruned;
}


int taf_norm_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;
    bool run_length_encode_bases = 0;
    bool output_maf = 0;
    bool use_compression = 0;
    bool filter_gap_causing_dupes = 0;
    stList *fasta_files = stList_construct();
    char *hal_file = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "maf", no_argument, 0, 'k' },
                                                { "help", no_argument, 0, 'h' },
                                                { "maximumBlockLengthToMerge", required_argument, 0, 'm' },
                                                { "maximumGapLength", required_argument, 0, 'n' },
                                                { "fractionSharedRows", required_argument, 0, 'q' },
                                                { "minimumSharedRows", required_argument, 0, 'Q' },
                                                { "filterGapCausingDupes", no_argument, 0, 'd' },
                                                { "repeatCoordinatesEveryNColumns", required_argument, 0, 's' },
                                                { "useCompression", no_argument, 0, 'c' },
                                                { "halFile", required_argument, 0, 'a' },
                                                { "seqFiles", required_argument, 0, 'b' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:hcm:n:dkQ:q:s:a:b:", long_options, &option_index);
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
                maximum_block_length_to_merge = atol(optarg);
                break;
            case 'n':
                maximum_gap_length = atol(optarg);
                break;
            case 'd':
                filter_gap_causing_dupes = 1;
                break;
            case 'Q':
                minimum_shared_rows = atol(optarg);
                break;
            case 'q':
                fraction_shared_rows = atof(optarg);
                break;
            case 'c':
                use_compression = 1;
                break;
            case 's':
                repeat_coordinates_every_n_columns = atol(optarg);
                break;
            case 'a':
                hal_file = optarg;
                break;
            case 'b':
                // Parse the set of sequence files (this is a bit fragile - files can not start with a '-' character)
                optind--;
                for( ;optind < argc && *argv[optind] != '-'; optind++){
                    stList_append(fasta_files, argv[optind]);
                }
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
    st_logInfo("Maximum block length to merge : %" PRIi64 "\n", maximum_block_length_to_merge);
    st_logInfo("Maximum gap length : %" PRIi64 "\n", maximum_gap_length);
    st_logInfo("Filter gap-causing dupes : %d\n", (int)filter_gap_causing_dupes);
    st_logInfo("Output maf : %s\n", output_maf ? "true" : "false");
    st_logInfo("Repeat coordinates every n bases : %" PRIi64 "\n", repeat_coordinates_every_n_columns);
    st_logInfo("Fraction shared rows to merge adjacent blocks : %f\n", fraction_shared_rows);
    st_logInfo("Write compressed output : %s\n", use_compression ? "true" : "false");
    if (hal_file) {
        st_logInfo("HAL file string : %s\n", hal_file);
    } else {
        st_logInfo("Number of input FASTA files : %ld\n", argc - optind);
    }

    //////////////////////////////////////////////
    // Read in the sequences if joining over unaligned gaps
    //////////////////////////////////////////////

    stHash *fastas_map = NULL;
    stSet *hal_species = NULL;
    int hal_handle = -1;
    if (hal_file) {
        hal_species = load_sequences_from_hal_file(hal_file, &hal_handle);
    }
    else if(stList_length(fasta_files) > 0) {
        fastas_map = load_sequences_from_fasta_files(stList_getBackingArray(fasta_files), stList_length(fasta_files));
        stList_destruct(fasta_files);
    }

    //////////////////////////////////////////////
    // Read in the taf blocks and merge blocks that are sufficiently small
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    LW *output = LW_construct(outputFile == NULL ? stdout : fopen(outputFile, "w"), use_compression);
    LI *li = LI_construct(input);

    // Pass the header line to determine parameters and write the updated taf header
    Tag *tag = taf_read_header_2(li, &run_length_encode_bases);
    if(output_maf && run_length_encode_bases) { // Remove this tag from the maf output as not relevant
        tag = tag_remove(tag, "run_length_encode_bases");
    }
    output_maf ? maf_write_header(tag, output) : taf_write_header(tag, output);
    tag_destruct(tag);

    Alignment *alignment, *p_alignment = NULL, *p_p_alignment = NULL;
    while((alignment = get_next_taf_block(li, run_length_encode_bases)) != NULL) {
        if(p_alignment != NULL) {
            // First realign the rows in case we in the process of merging prior blocks we have
            // identified rows that can be merged
            alignment_link_adjacent(p_alignment, alignment, 1);

            bool merged = false;
            int64_t common_rows = alignment_number_of_common_rows(p_alignment, alignment);
            int64_t total_rows = alignment->row_number + p_alignment->row_number - common_rows;
            if (common_rows >= minimum_shared_rows &&
                common_rows >= total_rows * fraction_shared_rows &&
                (alignment_length(p_alignment) <= maximum_block_length_to_merge ||
                 alignment_length(alignment) <= maximum_block_length_to_merge)) {
                int64_t total_gap = alignment_total_gap_length(p_alignment);
                if (total_gap > maximum_gap_length && filter_gap_causing_dupes) {
                    // try to greedily filter dupes in order to get the gap length down
                    bool was_pruned = greedy_prune_by_gap(alignment, maximum_gap_length);
                    total_gap = alignment_total_gap_length(p_alignment);
                    assert(was_pruned == (total_gap <= maximum_gap_length));
                }
                if (total_gap <= maximum_gap_length) {
                    if(hal_species || fastas_map) { // Now add in any gap bases if sequences are provided
                        alignment_add_gap_strings(p_alignment, alignment, fastas_map, hal_handle, hal_species, -1);
                    }
                    p_alignment = alignment_merge_adjacent(p_alignment, alignment);
                    merged = true;
                }
            }
            if (!merged) {
                output_maf ? maf_write_block(p_alignment, output) : taf_write_block(p_p_alignment, p_alignment, run_length_encode_bases, repeat_coordinates_every_n_columns, output); // Write the maf block
                if(p_p_alignment != NULL) {
                    alignment_destruct(p_p_alignment, 1); // Clean up the left-most block
                }
                p_p_alignment = p_alignment;
                p_alignment = alignment;
            }
        }
        else {
            p_alignment = alignment;
        }
    }
    if(p_alignment != NULL) {
        output_maf ? maf_write_block(p_alignment, output) : taf_write_block(p_p_alignment, p_alignment, run_length_encode_bases, -1, output); // Write the last taf block
        alignment_destruct(p_alignment, 1);
        if(p_p_alignment != NULL) {
            alignment_destruct(p_p_alignment, 1);
        }
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    if(inputFile != NULL) {
        fclose(input);
    }
    LW_destruct(output, outputFile != NULL);

    if (fastas_map) {
        stHash_destruct(fastas_map);
    }
    if (hal_species) {
        stSet_destruct(hal_species);
    }

    st_logInfo("taffy norm is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}


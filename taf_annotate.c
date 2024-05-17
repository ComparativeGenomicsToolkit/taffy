/*
 * taf annotate: TAF annotation
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "tai.h"
#include "sonLib.h"
#include <getopt.h>
#include <time.h>

static int64_t repeat_coordinates_every_n_columns = 10000;

static void usage(void) {
    fprintf(stderr, "taffy annotate [options]\n");
    fprintf(stderr, "Annotate the columns of a taf file using wiggle file\n");
    fprintf(stderr, "-i --inputFile : Input TAF file. If not specified reads from stdin\n");
    fprintf(stderr, "-w --wiggle [FILE_NAME] : REQUIRED The input wiggle file\n");
    fprintf(stderr, "-t --tagName [STRING] : REQUIRED: The name of the tag to annotate for the given wiggle\n");
    fprintf(stderr, "-s --repeatCoordinatesEveryNColumns : Repeat coordinates of each sequence at least every n columns. By default: %" PRIi64 "\n", repeat_coordinates_every_n_columns);
    fprintf(stderr, "-c --useCompression : Write the output using bgzip compression.\n");
    fprintf(stderr, "-r --refPrefix : Prefix to prepend to chrom names in annotation file to form the sequence name.\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

/*
 * Annotate an alignment block
 */
static void label_alignment(Alignment *alignment, stHash *labels, char *tag_name) {
    if(alignment->row_number > 0) {
        Alignment_Row *ref_row = alignment->row;
        assert(ref_row->strand); // Reference row is by convention on the positive strand
        stHash *seq_labels = stHash_search(labels, ref_row->sequence_name);
        if (seq_labels != NULL) {
            int64_t start = ref_row->start;
            for (int64_t i = 0; i < alignment->column_number; i++) { // For each column
                if (ref_row->bases[i] != '-') { // If a reference coordinate (not a gap)
                    double *label = stHash_search(seq_labels, (void *) start++);
                    if (label) {
                        assert(tag_find(alignment->column_tags[i], tag_name) ==
                               NULL); // No existing tag with this label
                        // Add the tag
                        alignment->column_tags[i] = tag_construct(tag_name, stString_print("%f", label[0]),
                                                                  alignment->column_tags[i]);
                    }
                }
            }
        }
    }
}

int taf_annotate_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *tag_name = NULL;
    char *taf_file = NULL;
    char *wig_file = NULL;
    char *output_file = NULL;
    bool use_compression = 0;
    char *ref_prefix = "";

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "wiggle", required_argument, 0, 'w' },
                                                { "tagName", required_argument, 0, 't' },
                                                { "repeatCoordinatesEveryNColumns", required_argument, 0, 's' },
                                                { "useCompression", no_argument, 0, 'c' },
                                                { "refPrefix", required_argument, 0, 'r' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:w:s:cr:ht:", long_options, &option_index);
        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'i':
                taf_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'w':
                wig_file = optarg;
                break;
            case 't':
                tag_name = optarg;
                break;
            case 's':
                repeat_coordinates_every_n_columns = atol(optarg);
                break;
            case 'c':
                use_compression = 1;
                break;
            case 'r':
                ref_prefix = optarg;
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

    if(!tag_name) {
        st_errAbort("No tag name given\n");
    }
    if(!wig_file) {
        st_errAbort("No wiggle file name given\n");
    }

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Input file string : %s\n", taf_file);
    st_logInfo("Output file string : %s\n", output_file);
    st_logInfo("Wig file string : %s\n", wig_file);
    st_logInfo("Tag name string : %s\n", tag_name);
    st_logInfo("Ref prefix string : %s\n", ref_prefix);

    //////////////////////////////////////////////
    // Open the inputs and outputs and parse the labels to annotate
    //////////////////////////////////////////////

    // load the input
    FILE *taf_fh = taf_file == NULL ? stdin : fopen(taf_file, "r");
    if (taf_fh == NULL) {
        fprintf(stderr, "Unable to open input TAF file: %s\n", taf_file);
        return 1;
    }
    LI *li = LI_construct(taf_fh);

    // Check is a taf file
    if (check_input_format(LI_peek_at_next_line(li)) != 0) {
        fprintf(stderr, "Input not supported: requires #taf header\n");
        return 1;
    }

    // Check if run_length_encode_bases is set and read header
    bool run_length_encode_bases = 0;
    Tag *tag = taf_read_header_2(li, &run_length_encode_bases);

    // Load the wiggle file, making coordinates 0 based
    stHash *labels = wig_parse(wig_file, ref_prefix, 1);

    // Open the output file for writing
    FILE *output_fh = output_file == NULL ? stdout : fopen(output_file, "w");
    if (output_fh == NULL) {
        fprintf(stderr, "Unable to open output file: %s\n", output_file);
        return 1;
    }
    LW *output = LW_construct(output_fh, use_compression);

    //////////////////////////////////////////////
    // Write the labelled taf
    //////////////////////////////////////////////

    // Write the taf header
    taf_write_header(tag, output);

    // Now write the body of the taf file
    Alignment *alignment = NULL;
    Alignment *p_alignment = NULL;
    // Keep reading blocks while available
    while((alignment = taf_read_block(p_alignment, run_length_encode_bases, li)) != NULL) {
        label_alignment(alignment, labels, tag_name); // Make any changes to the alignment for output

        // Write back the labelled taf
        taf_write_block2(p_alignment, alignment, run_length_encode_bases,
                         repeat_coordinates_every_n_columns, output, 0, 0);

        // Clean up
        if (p_alignment) {
            alignment_destruct(p_alignment, true);
        }
        p_alignment = alignment;
    }
    if (p_alignment) {
        alignment_destruct(p_alignment, true);
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    LI_destruct(li);
    if(taf_file != NULL) {
        fclose(taf_fh);
    }
    LW_destruct(output, output_file != NULL);
    tag_destruct(tag);

    st_logInfo("taffy annotate is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}


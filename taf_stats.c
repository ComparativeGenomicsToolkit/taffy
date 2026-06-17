/*
 * taf view: MAF/TAF conversion and subregion extraction
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "tai.h"
#include "tui.h"
#include "block_reader.h"
#include "sonLib.h"
#include "sonLibTree.h"
#include <getopt.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

// Walk an stTree, collecting every internal-node label into `out` (as
// fresh stString_copy'd char*).  Used by -s --anchorsOnly to build the
// set of valid row-0 anchor genome labels from the `# hal` tree.
static void collect_internal_labels(stTree *node, stSet *out) {
    int64_t nc = stTree_getChildNumber(node);
    if (nc > 0) {
        const char *lbl = stTree_getLabel(node);
        if (lbl != NULL && *lbl != '\0') {
            stSet_insert(out, stString_copy(lbl));
        }
    }
    for (int64_t i = 0; i < nc; i++) {
        collect_internal_labels(stTree_getChild(node, i), out);
    }
}

// Try each '.' in `seq_name` as a split point and return true if any
// prefix matches a label in `internal_labels`.  Mirrors the same
// genome-from-sequence resolution as gerp_tree_resolve_genome, but
// targeted at the internal-label set (= true row-0 anchors).
static bool seq_is_anchor(const char *seq_name, stSet *internal_labels) {
    size_t n = strlen(seq_name);
    char  stack_buf[1024];
    char *buf = (n + 1 <= sizeof(stack_buf)) ? stack_buf : st_malloc(n + 1);
    bool match = false;
    size_t off = 0;
    while (off < n) {
        const char *dot = strchr(seq_name + off, '.');
        if (dot == NULL || dot == seq_name + n - 1) break;
        if (dot == seq_name) { off++; continue; }
        size_t pre_n = (size_t)(dot - seq_name);
        memcpy(buf, seq_name, pre_n);
        buf[pre_n] = '\0';
        if (stSet_search(internal_labels, buf) != NULL) { match = true; break; }
        off = pre_n + 1;
    }
    if (buf != stack_buf) free(buf);
    return match;
}

static void usage(void) {
    fprintf(stderr, "taffy stats [options]\n");
    fprintf(stderr, "Print statistics from a TAF or MAF file\n");
    fprintf(stderr, "-i --inputFile : Input TAF or MAF file. If not specified reads from stdin\n");
    fprintf(stderr, "-s --sequenceLengths : Print length of each *reference* sequence in the (indexed) alignment.  Reads .tai if present, else falls back to the .tui's directory (which covers every sequence in a universal MAF -- both row-0 anchor chroms and leaf-genome chroms).\n");
    fprintf(stderr, "-a --alignmentStats : Print stats about block number, aligned bases, etc.\n");
    fprintf(stderr, "-b --sequenceIntervals : Print the BED intervals of each *reference* sequence covered by the alignment\n");
    fprintf(stderr, "-u --universalColumns : Print the total universal column count T (one integer line).  Requires a .tui.\n");
    fprintf(stderr, "   --anchorsOnly : Modifier on -s.  Filter the output to sequences whose genome label is an internal-node label in the `# hal` tree -- i.e. row-0 anchor chroms in a universal MAF.  Requires a .tui and a `# hal` header.\n");
    fprintf(stderr, "-T --threads N : Use N threads for bgzf I/O (default 1, only effective on bgzipped streams)\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int taf_stats_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *taf_fn = NULL;
    bool seq_lengths = false;
    bool seq_intervals = false;
    int stat_option_count = 0;
    bool alignment_stats = false;
    bool universal_columns = false;
    bool anchors_only = false;
    int bgzf_threads = 1;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    enum { OPT_ANCHORS_ONLY = 256 };
    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "sequenceLengths", no_argument, 0, 's' },
                                                { "alignmentStats", no_argument, 0, 'a' },
                                                { "sequenceIntervals", no_argument, 0, 'b' },
                                                { "universalColumns", no_argument, 0, 'u' },
                                                { "anchorsOnly", no_argument, 0, OPT_ANCHORS_ONLY },
                                                { "threads", required_argument, 0, 'T' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:sbauhT:", long_options, &option_index);
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
            case 's':
                seq_lengths = 1;
                ++stat_option_count;
                break;
            case 'a':
                alignment_stats = 1;
                ++stat_option_count;
                break;
            case 'b':
                seq_intervals = 1;
                ++stat_option_count;
                break;
            case 'u':
                universal_columns = 1;
                ++stat_option_count;
                break;
            case OPT_ANCHORS_ONLY:
                anchors_only = 1;
                break;
            case 'T':
                bgzf_threads = atoi(optarg);
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
    LI_set_bgzf_threads(bgzf_threads);
    st_logInfo("Input file string : %s\n", taf_fn);

    //////////////////////////////////////////////
    // Do the stats
    //////////////////////////////////////////////

    if (stat_option_count != 1) {
        fprintf(stderr, "Please pick a stats option from { -s, -b, -a, -u }\n");
        return 1;
    }
    if (anchors_only && !seq_lengths) {
        fprintf(stderr, "--anchorsOnly is a modifier on -s/--sequenceLengths\n");
        return 1;
    }

    // -u --universalColumns: straight read of T from the .tui.  Skips the
    // input/header dance below since we only need the index.
    if (universal_columns) {
        if (taf_fn == NULL) {
            fprintf(stderr, "-u requires -i <file> (cannot index stdin)\n");
            return 1;
        }
        char *tui_fn = tui_path(taf_fn);
        if (access(tui_fn, F_OK) != 0) {
            fprintf(stderr, "-u requires a .tui index: %s\n", tui_fn);
            free(tui_fn);
            return 1;
        }
        Tui *tui = tui_load(tui_fn);
        if (tui == NULL) {
            fprintf(stderr, "Failed to load .tui: %s\n", tui_fn);
            free(tui_fn);
            return 1;
        }
        fprintf(stdout, "%" PRIi64 "\n", tui_total_columns(tui));
        tui_destruct(tui);
        free(tui_fn);
        st_logInfo("taffy stats is done, %" PRIi64 " seconds have elapsed\n",
                   time(NULL) - startTime);
        return 0;
    }

    // load the input
    FILE *taf_fh = taf_fn == NULL ? stdin : fopen(taf_fn, "r");
    if (taf_fh == NULL) {
        fprintf(stderr, "Unable to open input TAF/MAF file: %s\n", taf_fn);
        return 1;
    }
    LI *li = LI_construct(taf_fh);

    // open a format-agnostic reader (handles MAF or TAF)
    BlockReader *reader = NULL;
    bool input_is_maf = false;
    if (!seq_lengths) {
        // seq_lengths goes through tai_sequence_lengths which manages its own header reading;
        // for the other paths we need to consume the header up-front via BlockReader
        reader = block_reader_open(li);
        if (reader == NULL) {
            LI_destruct(li);
            if (taf_fn != NULL) fclose(taf_fh);
            return 1;
        }
        input_is_maf = block_reader_is_maf(reader);
        tag_destruct(block_reader_take_header(reader));
    } else {
        // for -s we still need to know the format to load the index correctly
        int input_format = check_input_format(LI_peek_at_next_line(li));
        if (input_format != 0 && input_format != 1) {
            fprintf(stderr, "Input not supported: unable to detect ##maf or #taf header\n");
            LI_destruct(li);
            if (taf_fn != NULL) fclose(taf_fh);
            return 1;
        }
        input_is_maf = (input_format == 1);
    }

    // load the index if it's required by the given options.  For seq_lengths
    // we accept either a .tai (regular MAF) or a .tui (universal MAF) -- the
    // .tui's directory has every sequence's length so the same -s output
    // shape works on it.  Other stats options keep their .tai-only shape.
    //
    // --anchorsOnly forces the .tui path (the concept doesn't apply to
    // .tai's single-reference index) and additionally requires the input
    // header to carry a `# hal` tree, since we need its internal-node label
    // set to identify row-0 anchor chroms.
    bool index_required = seq_lengths;
    char *tai_fn = NULL;
    char *tui_fn = NULL;
    FILE *tai_fh = NULL;
    Tai *tai = NULL;
    bool used_tui = false;
    if (index_required) {
        if (taf_fn == NULL) {
            fprintf(stderr, "An index is needed to compute the requested stats so an input filename must be specified with -i\n");
            return 1;
        }
        if (!anchors_only) {
            tai_fn = tai_path(taf_fn);
            // Reject a .tai that's older than the .maf/.taf -- cactus
            // pipelines regenerate the alignment but sometimes forget to
            // re-run `taffy index`, and stale .tai sequence lengths are
            // silently wrong.  Log + fall through to the .tui path; the
            // .tui builder is always re-run with the alignment.
            struct stat tai_st, taf_st;
            if (stat(tai_fn, &tai_st) == 0 && stat(taf_fn, &taf_st) == 0 &&
                tai_st.st_mtime < taf_st.st_mtime) {
                st_logCritical("taffy stats: .tai %s is older than %s -- "
                               "treating as stale and falling back to .tui "
                               "(re-run `taffy index` if you want -s via .tai)\n",
                               tai_fn, taf_fn);
            } else {
                tai_fh = fopen(tai_fn, "r");
            }
        }
        if (tai_fh != NULL) {
            tai = tai_load(tai_fh, input_is_maf);
        } else {
            // Fall back (or force) to .tui.
            tui_fn = tui_path(taf_fn);
            if (access(tui_fn, F_OK) == 0) {
                used_tui = true;
            } else {
                fprintf(stderr, "Required index %s not found (and no .tui sibling either). Please run taffy index first\n",
                        tai_fn ? tai_fn : tui_fn);
                return 1;
            }
        }
    }

    // For --anchorsOnly: pull the `# hal` tree from the header now (the LI
    // has just had the format detected via peek; the header is still
    // unread).  Parsing produces the internal-label set used to filter the
    // -s output below.
    stSet *internal_labels = NULL;
    stTree *tree_for_labels = NULL;
    if (anchors_only) {
        bool rle = false;
        Tag *header_tag = (input_is_maf
                           ? maf_read_header(li)
                           : taf_read_header_2(li, &rle));
        Tag *hal = (header_tag != NULL)
                   ? tag_find(header_tag, (char *)TAF_HAL_TREE_KEY)
                   : NULL;
        if (hal == NULL) {
            fprintf(stderr, "--anchorsOnly requires a `# hal` tree in the input header.\n");
            tag_destruct(header_tag);
            return 1;
        }
        // stTree_parseNewickString doesn't mutate but takes char* per
        // sonLib's signature.
        tree_for_labels = stTree_parseNewickString(hal->value);
        if (tree_for_labels == NULL) {
            fprintf(stderr, "--anchorsOnly: failed to parse `# hal` Newick tree.\n");
            tag_destruct(header_tag);
            return 1;
        }
        internal_labels = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
        collect_internal_labels(tree_for_labels, internal_labels);
        tag_destruct(header_tag);
        st_logInfo("--anchorsOnly: %" PRIi64 " internal-node labels collected\n",
                   stSet_size(internal_labels));
    }

    // do the stats
    if (seq_lengths) {
        stHash *seq_to_len = NULL;
        if (used_tui) {
            // tui_sequence_lengths now reuses a loaded handle's cursor; load
            // one just for this.  Its output copies the keys and stores int
            // lengths, so the handle can be freed right after.
            Tui *stui = tui_load(tui_fn);
            if (stui != NULL) {
                seq_to_len = tui_sequence_lengths(stui);
                tui_destruct(stui);
            }
        } else {
            seq_to_len = tai_sequence_lengths(tai, li);
        }
        if (seq_to_len == NULL) {
            fprintf(stderr, "Failed to read sequence directory from %s\n",
                    used_tui ? tui_fn : tai_fn);
            free(tai_fn); free(tui_fn);
            return 1;
        }
        stList *seq_names = stHash_getKeys(seq_to_len);
        int64_t n_emitted = 0, n_filtered = 0;
        for (int64_t i = 0; i < stList_length(seq_names); ++i) {
            const char *nm = (const char *)stList_get(seq_names, i);
            if (anchors_only && !seq_is_anchor(nm, internal_labels)) {
                n_filtered++;
                continue;
            }
            void *hash_val = stHash_search(seq_to_len, (void *)nm);
            fprintf(stdout, "%s\t%" PRIi64 "\n", nm, (int64_t)hash_val);
            n_emitted++;
        }
        if (anchors_only) {
            st_logInfo("--anchorsOnly: emitted %" PRIi64 " / %" PRIi64 " sequences "
                       "(filtered %" PRIi64 " non-anchor)\n",
                       n_emitted, n_emitted + n_filtered, n_filtered);
        }
        stHash_destruct(seq_to_len);
        stList_destruct(seq_names);
    } else if (seq_intervals) {
        Alignment *alignment = NULL;
        Alignment *p_alignment = NULL;
        char *cur_seq = NULL;
        int64_t cur_start = -1;
        int64_t cur_end = 0;
        while ((alignment = block_reader_next(reader, p_alignment)) != NULL) {
            if (alignment->row_number > 0) {
                if (!cur_seq || strcmp(cur_seq, alignment->row->sequence_name) != 0 || alignment->row->start != cur_end) {
                    if (cur_seq) {
                        fprintf(stdout, "%s\t%" PRIi64 "\t%" PRIi64 "\n", cur_seq, cur_start, cur_end);
                        free(cur_seq);
                    }
                    cur_seq = stString_copy(alignment->row->sequence_name);
                    cur_start = alignment->row->start;
                    cur_end = cur_start + alignment->row->length;
                } else {
                    cur_end += alignment->row->length;
                }
            }
            if (p_alignment) {
                alignment_destruct(p_alignment, true);
            }
            p_alignment = alignment;
        }
        if (p_alignment) {
            alignment_destruct(p_alignment, true);
        }
        if (cur_seq) {
            fprintf(stdout, "%s\t%" PRIi64 "\t%" PRIi64 "\n", cur_seq, cur_start, cur_end);
            free(cur_seq);
        }
    }

    // If want column depth stats (does not currently work with any subregion)
    if(alignment_stats) {
        int64_t total_blocks = 0, total_columns = 0, total_aligned_bases = 0, total_gaps = 0, total_column_depth = 0;
        Alignment *alignment, *p_alignment = NULL;
        while ((alignment = block_reader_next(reader, p_alignment)) != NULL) {
            total_blocks++;
            total_columns += alignment_length(alignment);
            total_column_depth += alignment_length(alignment) * alignment->row_number;
            Alignment_Row *row = alignment->row;
            while (row != NULL) {
                for (int64_t i = 0; i < alignment_length(alignment); i++) {
                    if (row->bases[i] == '-') {
                        total_gaps++;
                    } else {
                        total_aligned_bases++;
                    }
                }
                row = row->n_row;
            }
            if(p_alignment != NULL) {
                alignment_destruct(p_alignment, 1);
            }
            p_alignment = alignment;
        }
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment, 1);
        }
        fprintf(stdout, "Total blocks:\t%" PRIi64 "\n", total_blocks);
        fprintf(stdout, "Total columns:\t%" PRIi64 "\n", total_columns);
        fprintf(stdout, "Avg. columns/block:\t%f\n", (float)total_columns/total_blocks);
        fprintf(stdout, "Total bases:\t%" PRIi64 "\n", total_aligned_bases);
        fprintf(stdout, "Total gaps:\t%" PRIi64 "\n", total_gaps);
        fprintf(stdout, "Avg. column depth:\t%f\n", (float)total_column_depth/total_columns);
        fprintf(stdout, "Avg. bases/column:\t%f\n", (float)total_aligned_bases/total_columns);
        fprintf(stdout, "Avg. gaps/column:\t%f\n", (float)total_gaps/total_columns);
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    if (index_required) {
        free(tai_fn);
        free(tui_fn);
        if (tai != NULL) tai_destruct(tai);
        if (tai_fh != NULL) fclose(tai_fh);
    }
    if (internal_labels != NULL) stSet_destruct(internal_labels);
    if (tree_for_labels != NULL) stTree_destruct(tree_for_labels);

    block_reader_destruct(reader);
    LI_destruct(li);
    if(taf_fn != NULL) {
        fclose(taf_fh);
    }

    st_logInfo("taffy stats is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}


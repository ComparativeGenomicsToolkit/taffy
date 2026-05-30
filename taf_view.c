/*
 * taf view: MAF/TAF conversion and subregion extraction
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "taf.h"
#include "tai.h"
#include "tui.h"
#include "remote_io.h"
#include "sonLib.h"
#include <getopt.h>
#include <time.h>
#include <unistd.h>

static int64_t repeat_coordinates_every_n_columns = 10000;
static bool show_only_reference_differences = false;
static bool show_only_lineage_differences = false;
static stTree *phylogeny = NULL;
static stList *sequence_prefixes = NULL;
static stList *tree_nodes = NULL;

static void usage(void) {
    fprintf(stderr, "taffy view [options]\n");
    fprintf(stderr, "Convert between TAF and MAF formats\n");
    fprintf(stderr, "-i --inputFile : Input TAF or MAF file to convert. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output file. If not specified outputs to stdout\n");
    fprintf(stderr, "-m --maf : Output in MAF format [default=TAF format]\n");
    fprintf(stderr, "-p --paf : Output in all-to-one (referenced on first row) PAF format [default=TAF format]\n");
    fprintf(stderr, "-P --paf-all : Output in all-to-all PAF format [default=TAF format]\n");
    fprintf(stderr, "-C --cs : Output PAF cigars in cs instead of cg format\n");
    fprintf(stderr, "-r --region  : Region SEQ:START-END (0-based half-open).  SEQ is a row-0 seq name (.tai path), or ANY genome.seq when a <inputFile>.tui is present (universal lift), or the sentinel `tcol` for a universal-column range\n");
    fprintf(stderr, "-U --universal MODE : Output mode for the universal lift {ancestor|tcol|query}.  Universal mode is AUTO-ENGAGED when <inputFile>.tui is present; -U just picks the mode.  ancestor = pass-through (default).  tcol = prepend a `tcol` sentinel row carrying the universal column at this block.  query = reorient blocks onto the queried genome (row-0, '+'); incompatible with `-r tcol:..`.  Blocks are always emitted in universal-column (file) order, which on a '-'-strand queried leaf descends in the queried-forward coordinate.\n");
    fprintf(stderr, "-s --repeatCoordinatesEveryNColumns : Repeat TAF coordinates of each sequence at least every n columns. By default: %" PRIi64 "\n", repeat_coordinates_every_n_columns);
    fprintf(stderr, "-u --runLengthEncodeBases : Run length encode output bases in TAF\n");
    fprintf(stderr, "-a --showOnlyReferenceDifferences : Replace matches with the reference (first row) with a * character\n");
    fprintf(stderr, "-b --showOnlyLineageDifferences : Show only lineage changes, replacing identity with a * character.\n "
                    "Requires the phylogeny to be specified and ancestors to be present\n");
    fprintf(stderr, "-x --colorBases : Color bases in the output. THIS IS FOR VISUALIZATION ONLY - DOES NOT PRODUCE A VALID MAF/TAF\n");
    fprintf(stderr, "-t --phylogeny [newick tree file] : Specify a file containing the phylogeny for the alignment, where species names must be prefixes of sequence names\n");
    fprintf(stderr, "-d --omitCoordinates : When printing TAF, just print the columns omitting the coordinates. THIS IS FOR VISUALIZATION ONLY - DOES NOT PRODUCE A VALID TAF \n");
    fprintf(stderr, "-c --useCompression : Write the output using bgzip compression.\n");
    fprintf(stderr, "-n --nameMapFile : Apply the given two-column tab-separated name mapping to all assembly names in alignment\n");
    fprintf(stderr, "-T --threads N : Use N threads for bgzf I/O (default 1, only effective on bgzipped streams)\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

void getTreeNodesInList(stTree *root, stList *tree_nodes) {
    stList_append(tree_nodes, root);
    for (int i = 0; i < stTree_getChildNumber(root); i++) {
        getTreeNodesInList(stTree_getChild(root, i), tree_nodes);
    }
}

void get_sequence_prefixes_for_tree_nodes(void) {
    // First get list of sequence prefixes to efficiently locate which sequence goes with which tree node
    sequence_prefixes = stList_construct3(0, (void (*)(void *))sequence_prefix_destruct);
    tree_nodes = stList_construct();
    getTreeNodesInList(phylogeny, tree_nodes);
    // For each tree node
    for(int64_t i=0; i<stList_length(tree_nodes); i++) {
        stTree *node = stList_get(tree_nodes, i);
        // Make a sequence prefix object
        stList_append(sequence_prefixes, sequence_prefix_construct(stString_copy(stTree_getLabel(node)), i));
    }
    // Sort the sequence prefixes
    stList_sort(sequence_prefixes, (int (*)(const void *, const void *))sequence_prefix_cmp_fn);
}

static void modify_alignment(Alignment *alignment) {
    if(show_only_reference_differences) {
        alignment_mask_reference_bases(alignment, '*');
    }
    else if(show_only_lineage_differences) {
        alignment_show_only_lineage_differences(alignment, '*', sequence_prefixes, tree_nodes);
    }
}

int taf_view_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;
    bool run_length_encode_output_bases = false;
    bool maf_output = false;
    bool paf_output = false;
    bool all_to_all_paf = false;
    bool paf_cs = false;
    char *region = NULL;
    // -U/--universal MODE: ancestor (default; pass-through), tcol (prepend
    // `tcol` sentinel row carrying universal-column coord), query (reorient
    // onto the queried genome).  Universal mode auto-engages when <input>.tui
    // is present; -U just picks the mode.  -U with no .tui is an error.
    enum { U_MODE_ANCESTOR, U_MODE_TCOL, U_MODE_QUERY };
    int universal_mode = U_MODE_ANCESTOR;
    bool universal_mode_set = false;
    bool use_compression = false;
    char *nameMapFile = NULL;
    char *phylogeny_file = NULL;
    static bool color_bases = false;
    bool omit_coordinates = false;
    int bgzf_threads = 1;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "maf", no_argument, 0, 'm' },
                                                { "paf", no_argument, 0, 'p' },
                                                { "paf-all", no_argument, 0, 'P' },
                                                { "cs", no_argument, 0, 'C'},
                                                { "runLengthEncodeBases", no_argument, 0, 'u' },
                                                { "showOnlyReferenceDifferences", no_argument, 0, 'a' },
                                                { "showOnlyLineageDifferences", no_argument, 0, 'b' },
                                                { "colorBases", no_argument, 0, 'x' },
                                                { "phylogeny", required_argument, 0, 't' },
                                                { "omitCoordinates", required_argument, 0, 'd' },
                                                { "repeatCoordinatesEveryNColumns", required_argument, 0, 's' },
                                                { "region", required_argument, 0, 'r' },
                                                { "universal", required_argument, 0, 'U' },
                                                { "useCompression", no_argument, 0, 'c' },
                                                { "nameMapFile", required_argument, 0, 'n' },
                                                { "threads", required_argument, 0, 'T' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:mPpCaucs:r:n:habxt:dU:T:", long_options, &option_index);
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
            case 'P':
                paf_output = 1;
                all_to_all_paf = 1;
                break;
            case 'C':
                paf_cs = true;
                break;
            case 'u':
                run_length_encode_output_bases = 1;
                break;
            case 't':
                phylogeny_file = optarg;
                break;
            case 'a':
                show_only_reference_differences = 1;
                break;
            case 'b':
                show_only_lineage_differences = 1;
                break;
            case 'x':
                color_bases = 1;
                break;
            case 'd':
                omit_coordinates = 1;
                break;
            case 's':
                repeat_coordinates_every_n_columns = atol(optarg);
                break;
            case 'U':
                if (strcmp(optarg, "ancestor") == 0)      universal_mode = U_MODE_ANCESTOR;
                else if (strcmp(optarg, "tcol") == 0)     universal_mode = U_MODE_TCOL;
                else if (strcmp(optarg, "query") == 0)    universal_mode = U_MODE_QUERY;
                else { fprintf(stderr,
                    "Invalid -U mode '%s'; expected ancestor|tcol|query\n", optarg);
                    return 1; }
                universal_mode_set = true;
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
    st_logInfo("Input file string : %s\n", inputFile);
    st_logInfo("Output file string : %s\n", outputFile);
    st_logInfo("Write compressed output : %s\n", use_compression ? "true" : "false");
    if (nameMapFile) {
        st_logInfo("Name map file string : %s\n", nameMapFile);
    }
    if(phylogeny_file) {
        st_logInfo("Phylogeny file : %s\n", phylogeny_file);
    }
    st_logInfo("Show only reference differences : %s\n", show_only_reference_differences ? "true" : "false");
    st_logInfo("Show only lineage differences : %s\n", show_only_lineage_differences ? "true" : "false");
    st_logInfo("Color bases : %s\n", color_bases ? "true" : "false");

    //////////////////////////////////////////////
    // Read in the taf/maf blocks and convert to sequence of taf/maf blocks
    //////////////////////////////////////////////

    if (paf_cs && !paf_output) {
        fprintf(stderr, "-C/--cs cannot be used without either -p or -P\n");
        return 1;
    }

    /* For URL inputs (only meaningful with -r region queries) we'll skip the
     * local fopen; LI is constructed via bgzf_open which goes through htslib's
     * URL-aware backend. */
    bool input_is_url = inputFile != NULL && is_url(inputFile);
    FILE *input = NULL;
    if (!input_is_url) {
        input = inputFile == NULL ? stdin : fopen(inputFile, "r");
        if (input == NULL) {
            fprintf(stderr, "Unable to open input file: %s\n", inputFile);
            return 1;
        }
    } else if (region == NULL) {
        fprintf(stderr, "URL inputs are only supported with -r region queries\n");
        return 1;
    }

    FILE *output_fh = outputFile == NULL ? stdout : fopen(outputFile, "w");
    if (output_fh == NULL) {
        fprintf(stderr, "Unable to open output file: %s\n", outputFile);
        return 1;
    }

    // Parse the tree if present
    if(phylogeny_file) {
        FILE *f = fopen(phylogeny_file, "r");
        char *phylogeny_string = stFile_getLineFromFile(f);
        fclose(f);
        phylogeny = stTree_parseNewickString(phylogeny_string);
        free(phylogeny_string);
        get_sequence_prefixes_for_tree_nodes();
    }

    stHash *genome_name_map = NULL;
    if (nameMapFile != NULL) {
        genome_name_map = load_genome_name_mapping(nameMapFile);
    }
    
    LW *output = LW_construct(output_fh, use_compression);
    LI *li = input_is_url ? LI_construct_from_path(inputFile) : LI_construct(input);
    if (li == NULL) {
        return 1;
    }

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
    bool run_length_encode_input_bases = false;
    Tag *tag = maf_input ? maf_read_header(li) : taf_read_header(li);
    if (maf_input && !maf_output) {
        if(run_length_encode_output_bases) {
            tag = tag_construct("run_length_encode_bases", "1", tag);
        }
    }
    else if (taf_input) {
        Tag *t = tag_find(tag, "run_length_encode_bases");
        if(t != NULL && strcmp(t->value, "1") == 0) {
            run_length_encode_input_bases = 1;
            if (maf_output || !run_length_encode_output_bases) { // If the output is a maf or we don't want run length
                // encoding in the output
                tag = tag_remove(tag, "run_length_encode_bases");  // Remove this tag from the maf
                // output as not relevant
            }
        }
        else if(!maf_output && run_length_encode_output_bases) { // If is taf input, taf output and we're turning on RLE
            tag = tag_construct("run_length_encode_bases", "1", tag);
        }
    }

    if (region && inputFile == NULL) {
        fprintf(stderr, "An input file must be specified with -i in order to perform region queries.\n");
        return 1;
    }

    if (universal_mode_set && region == NULL) {
        fprintf(stderr, "-U/--universal requires -r SEQ:START-END\n");
        return 1;
    }
    if (universal_mode_set && inputFile != NULL) {
        // fail loudly BEFORE header write, so a user-error `-U` invocation
        // doesn't leave a partial bgzip output (-c) with just the header.
        char *tp = tui_path(inputFile);
        if (access(tp, F_OK) != 0) {
            fprintf(stderr, "-U/--universal requires %s (run `taffy index -u` first)\n", tp);
            free(tp);
            return 1;
        }
        free(tp);
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

      // Auto-detect <inputFile>.tui -> engage universal mode.  Both indexes
      // present: .tui wins.  -U <mode> picks the OUTPUT shape (default
      // ancestor).  -U passed but no .tui -> error (don't silently downgrade).
      // (The -U-but-no-.tui case is already caught pre-header above.)
      // URL inputs skip the .tui path: there is no remote .tui resolver
      // yet, and access(F_OK) on a URL would fail spuriously anyway.
      char *tui_fn = input_is_url ? NULL : tui_path(inputFile);
      bool tui_present = (tui_fn != NULL && access(tui_fn, F_OK) == 0);
      if (tui_present) {
        bool tcol_input = (strcmp(region_seq, "tcol") == 0);
        if (tcol_input && universal_mode == U_MODE_QUERY) {
            fprintf(stderr, "-U query is incompatible with `-r tcol:..` "
                    "(no queried genome to reorient on)\n");
            free(tui_fn);
            return 1;
        }
        Tui *tui = tui_load(tui_fn);
        if (tui == NULL) {
            fprintf(stderr, "Could not open universal index %s\n", tui_fn);
            free(tui_fn);
            return 1;
        }
        free(tui_fn);
        int64_t T = tui_total_columns(tui);

        // Build the universal-column intervals for the extractor.
        // tcol:x-y bypasses Index B; everything else uses tui_query.
        TuiInterval *uiv = NULL;       // malloc'd iff !tcol_input
        // tcol_iv stays uninitialized on the n_uiv==0 path -- tui_extract_
        // iterator early-returns before reading any iv[i] when n_uiv==0, so
        // the garbage value is unreachable today.  Init defensively so a
        // future change to the extractor that peeks before the n-check
        // doesn't quietly read uninitialized stack.
        TuiInterval tcol_iv = {0, 0};
        int64_t n_uiv = 0;
        if (tcol_input) {
            int64_t a = region_start, b = region_start + region_length;
            if (a < 0) a = 0;
            if (b > T) b = T;
            if (a < b) { tcol_iv.start = a; tcol_iv.end = b; n_uiv = 1; }
        } else {
            uiv = tui_query(tui, region_seq, region_start,
                            region_start + region_length, &n_uiv);
        }
        const TuiInterval *iv_in = tcol_input ? &tcol_iv : uiv;
        TuiExtractIt *xit = tui_extract_iterator(tui, li, !taf_input,
                                                 run_length_encode_input_bases,
                                                 iv_in, n_uiv);
        if (!tui_extract_has_next(xit)) {
            fprintf(stderr, "Region %s:%" PRIi64 "-%" PRIi64 " maps to no "
                    "universal columns; emitting header-only output\n",
                    region_seq, region_start, region_start + region_length);
        }
        Alignment *alignment = NULL;
        Alignment *prev_alignment = NULL;  // kept alive across iterations for TAF delta
        while ((alignment = tui_extract_next(xit, li)) != NULL) {
            // Take ownership of THIS yield from the iterator so the next
            // tui_extract_next() doesn't free it -- we need it as
            // p_alignment in the NEXT iter's TAF delta write.  Without
            // this, taf_write_block2(NULL, ...) was the only safe choice
            // and the TAF emit would degenerate: each block written as
            // "all rows inserted" even when consecutive blocks share
            // contiguous rows (very common in -U query, where row-0 is
            // the queried genome and advances monotonically).  Pre-fix
            // the parser asserted on round-trip because new blocks
            // looked like "carry 9 + insert 9" with only 9 bases.
            tui_extract_take_ownership(xit);

            int64_t col_start = tui_extract_col_start(xit);
            // Apply the optional rename FIRST so the reorient match below
            // sees the post-name-map row sequence_name (same convention as
            // the rest of `taffy view`).
            if (genome_name_map) {
                apply_genome_name_mapping_to_alignment(genome_name_map, alignment);
            }
            // Mode transforms (before modify_alignment, so -a/-b see the
            // intended row-0): ancestor = no-op; tcol = prepend a `tcol`
            // sentinel row carrying this block's universal column; query =
            // reorient onto the queried genome (RC if '-', move to row-0).
            if (universal_mode == U_MODE_TCOL) {
                Alignment_Row *nr = st_calloc(1, sizeof(Alignment_Row));
                nr->sequence_name = stString_copy("tcol");
                nr->start = col_start;
                nr->length = alignment->column_number;     // gap-free by construction
                nr->strand = 1;
                nr->sequence_length = T;
                nr->bases = stString_copy(alignment->row->bases);  // lossless: per-col content
                nr->n_row = alignment->row;
                alignment->row = nr;
                alignment->row_number++;
            } else if (universal_mode == U_MODE_QUERY) {
                alignment_reorient_to_row(alignment, region_seq);
            }
            modify_alignment(alignment);
            if (taf_output) {
                // Match consecutive blocks' rows via WFA so the writer
                // can emit minimal i/d/s/g/G ops (rather than i-for-
                // every-row, which produces structurally invalid TAF).
                // alignment_link_adjacent sets r_row/l_row pointers that
                // write_coordinates reads to compute the deltas.  Same
                // helper the MAF reader's block_reader uses for the
                // streaming-read case.
                if (prev_alignment != NULL) {
                    alignment_link_adjacent(prev_alignment, alignment, 1);
                }
                taf_write_block2(prev_alignment, alignment, run_length_encode_output_bases,
                                 repeat_coordinates_every_n_columns, output, color_bases, omit_coordinates);
            } else if (maf_output) {
                // MAF has no inter-block row carry, so the link/p_alignment
                // dance is unnecessary -- preserve the pre-fix
                // optimisation.
                maf_write_block2(alignment, output, color_bases);
            } else {
                assert(paf_output == true);
                paf_write_block(alignment, output, all_to_all_paf, paf_cs);
            }

            // Rotate: drop the previously-held block, current becomes prev.
            if (prev_alignment != NULL) alignment_destruct(prev_alignment, 1);
            prev_alignment = alignment;
        }
        if (prev_alignment != NULL) alignment_destruct(prev_alignment, 1);
        tui_extract_iterator_destruct(xit);   // safe: most-recent yield was disowned
        free(uiv);                            // NULL for tcol input
        tui_destruct(tui);
      } else {
        free(tui_fn);

        char *tai_fn = input_is_url ? tai_path_for(inputFile) : tai_path(inputFile);
        FILE *tai_fh = open_tai_for_reading(tai_fn);
        if (tai_fh == NULL) {
            fprintf(stderr, "Index %s not found. Please run taffy index first\n", tai_fn);
            free(tai_fn);
            return 1;
        }
        Tai* tai = tai_load(tai_fh, !taf_input);

        TaiIt *tai_it = tai_iterator(tai, li, run_length_encode_input_bases, region_seq, region_start, region_length);
        if (!tai_has_next(tai_it)) {
            // No overlapping blocks were found. This can happen when the contig is missing from the
            // index, when the index has a stub entry for a contig with no blocks, or when the queried
            // sub-range simply doesn't intersect any indexed block. In all cases we exit cleanly so
            // batch callers don't have to special-case empty regions; the header was already written
            // above and the loop below is a no-op since tai_next returns NULL. We must NOT early-return
            // here -- the bgzip writer is only flushed/closed by LW_destruct in the cleanup at the end
            // of the function, so an early return would leave a 0-byte file when -c is set.
            fprintf(stderr, "Region %s:%" PRIi64 "-%" PRIi64 " not found in taffy index; emitting header-only output\n",
                    region_seq, region_start, region_start + region_length);
        }
        Alignment *alignment = NULL;
        Alignment *p_alignment = NULL;

        while ((alignment = tai_next(tai_it, li)) != NULL) {
            modify_alignment(alignment); // Make any changes to the alignment for output

            // apply the name mapping to the alignment block
            if (genome_name_map) {
                apply_genome_name_mapping_to_alignment(genome_name_map, alignment);
            }
            if (taf_output) {
                taf_write_block2(p_alignment, alignment, run_length_encode_output_bases,
                                 repeat_coordinates_every_n_columns, output, color_bases, omit_coordinates);
            } else if (maf_output) {
                maf_write_block2(alignment, output, color_bases);
            } else {
                assert(paf_output == true);
                paf_write_block(alignment, output, all_to_all_paf, paf_cs);
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
      }  // end non-universal single-region path (.tai)
    } else if (taf_input) {
        Alignment *alignment = NULL;
        Alignment *p_alignment = NULL;        
        while((alignment = taf_read_block(p_alignment, run_length_encode_input_bases, li)) != NULL) {
            modify_alignment(alignment); // Make any changes to the alignment for output

            // apply the name mapping to the alignment block
            if (genome_name_map) {
                apply_genome_name_mapping_to_alignment(genome_name_map, alignment);
            }
            if (taf_output) {
                taf_write_block2(p_alignment, alignment, run_length_encode_output_bases,
                                 repeat_coordinates_every_n_columns, output, color_bases, omit_coordinates);
            } else if (maf_output) {
                maf_write_block2(alignment, output, color_bases);
            } else {
                assert(paf_output == true);
                paf_write_block(alignment, output, all_to_all_paf, paf_cs);
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
            modify_alignment(alignment); // Make any changes to the alignment for output

            // apply the name mapping to the alignment block
            if (genome_name_map) {
                apply_genome_name_mapping_to_alignment(genome_name_map, alignment);
            }            
            if(p_alignment != NULL) {
                alignment_link_adjacent(p_alignment, alignment, 1);
            }
            if (taf_output) {
                taf_write_block2(p_alignment, alignment, run_length_encode_output_bases,
                                 repeat_coordinates_every_n_columns, output, color_bases, omit_coordinates);
            } else if (maf_output) {
                maf_write_block2(alignment, output, color_bases);
            } else {
                assert(paf_output == true);
                paf_write_block(alignment, output, all_to_all_paf, paf_cs);
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
    if (inputFile != NULL && !input_is_url) {
        fclose(input);
    }
    LW_destruct(output, outputFile != NULL);

    if (genome_name_map != NULL) {
        stHash_destruct(genome_name_map);
    }

    if(phylogeny_file) {
        stList_destruct(sequence_prefixes);
        stList_destruct(tree_nodes);
        stTree_destruct(phylogeny);
    }

    st_logInfo("taffy view is done, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}


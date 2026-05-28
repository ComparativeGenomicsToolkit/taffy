/*
 * taffy gerp: per-column GERP RS conservation scoring -> wig output.
 *
 * The reference (row 0) of each MAF/TAF block keys the wig record
 * (sequence_name + advancing position over non-gap row-0 bases).  In a
 * normal hg38-anchored MAF that's hg38; in a universal MAF row 0 varies
 * per block (one of many ancestor chroms) and the wig ends up with
 * multiple chroms naturally.
 *
 * See taffy/inc/gerp.h for the algorithm contract.
 */

#include "taf.h"
#include "tai.h"
#include "gerp.h"
#include "sonLib.h"
#include <ctype.h>
#include <getopt.h>
#include <time.h>

static void usage(void) {
    fprintf(stderr, "taffy gerp [options]\n");
    fprintf(stderr, "Compute per-column GERP RS conservation scores -> wig.\n");
    fprintf(stderr, "-i --inputFile : Input MAF or TAF (auto-detected; bgzipped ok).  Reads stdin if omitted.\n");
    fprintf(stderr, "-o --outputFile : Output wig.  Writes stdout if omitted.\n");
    fprintf(stderr, "-c --useCompression : Bgzip the output (matches taffy view).\n");
    fprintf(stderr, "-t --tree : Newick tree override.  Default: the `# hal` tree comment in the input header.\n");
    fprintf(stderr, "-D --depthFile : Optional second wig with the per-base count of A/C/G/T leaves at each column.\n");
    fprintf(stderr, "--branchScale : Global multiplier on branch lengths (default 1.0).\n");
    fprintf(stderr, "--minLeaves : Minimum surviving non-gap leaves to score a column (default 2).\n");
    fprintf(stderr, "--skipParalogs : Skip blocks with duplicate species (default ON).\n");
    fprintf(stderr, "--keepParalogs : Inverse of --skipParalogs (score using the first-seen row per species).\n");
    fprintf(stderr, "-T --threads N : Use N threads for bgzf I/O on bgzipped streams (default 1).\n");
    fprintf(stderr, "-l --logLevel : Set the log level.\n");
    fprintf(stderr, "-h --help : Print this help message.\n");
}

// Long-only options (no short flag).  Codes start at 256 so they don't
// collide with single-character short flags.
enum {
    OPT_BRANCH_SCALE = 256,
    OPT_MIN_LEAVES,
    OPT_SKIP_PARALOGS,
    OPT_KEEP_PARALOGS,
};

// Emit one column's wig record: "<pos> <value>\n".  Caller is responsible
// for having emitted the `variableStep chrom=...` header for this block.
static inline void emit_wig_value(LW *out, int64_t pos, double value) {
    char buf[64];
    // wig is 1-based; row->start is 0-based half-open, so position is
    // (row->start + col_offset) + 1.  Caller already passes the
    // 1-based wig position so we don't re-add here.
    LW_puti64(out, pos);
    LW_putc(out, ' ');
    int n = snprintf(buf, sizeof(buf), "%.4f", value);
    if (n < 0) n = 0;
    LW_putn(out, buf, (size_t)n);
    LW_putc(out, '\n');
}

static inline void emit_wig_header(LW *out, const char *chrom) {
    LW_putn(out, "variableStep chrom=", 19);
    LW_puts(out, chrom);
    LW_putc(out, '\n');
}

int taf_gerp_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    char *logLevelString = NULL;
    char *inputFile      = NULL;
    char *outputFile     = NULL;
    char *depthFile      = NULL;
    char *treeFile       = NULL;
    bool use_compression = false;
    double branch_scale  = 1.0;
    int64_t min_leaves   = 2;
    bool skip_paralogs   = true;
    int bgzf_threads     = 1;

    while (1) {
        static struct option long_options[] = {
            { "logLevel",       required_argument, 0, 'l' },
            { "inputFile",      required_argument, 0, 'i' },
            { "outputFile",     required_argument, 0, 'o' },
            { "useCompression", no_argument,       0, 'c' },
            { "tree",           required_argument, 0, 't' },
            { "depthFile",      required_argument, 0, 'D' },
            { "branchScale",    required_argument, 0, OPT_BRANCH_SCALE },
            { "minLeaves",      required_argument, 0, OPT_MIN_LEAVES },
            { "skipParalogs",   no_argument,       0, OPT_SKIP_PARALOGS },
            { "keepParalogs",   no_argument,       0, OPT_KEEP_PARALOGS },
            { "threads",        required_argument, 0, 'T' },
            { "help",           no_argument,       0, 'h' },
            { 0, 0, 0, 0 }
        };
        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:ct:D:T:h", long_options, &option_index);
        if (key == -1) break;
        switch (key) {
            case 'l': logLevelString = optarg; break;
            case 'i': inputFile      = optarg; break;
            case 'o': outputFile     = optarg; break;
            case 'c': use_compression = true;  break;
            case 't': treeFile       = optarg; break;
            case 'D': depthFile      = optarg; break;
            case 'T': bgzf_threads   = atoi(optarg); break;
            case OPT_BRANCH_SCALE:  branch_scale  = atof(optarg); break;
            case OPT_MIN_LEAVES:    min_leaves    = atol(optarg); break;
            case OPT_SKIP_PARALOGS: skip_paralogs = true;  break;
            case OPT_KEEP_PARALOGS: skip_paralogs = false; break;
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }

    st_setLogLevelFromString(logLevelString);
    LI_set_bgzf_threads(bgzf_threads);

    // Input.
    FILE *in_fh = (inputFile == NULL) ? stdin : fopen(inputFile, "r");
    if (in_fh == NULL) {
        fprintf(stderr, "taffy gerp: cannot open input file: %s\n", inputFile);
        return 1;
    }
    LI *li = LI_construct(in_fh);
    int input_format = check_input_format(LI_peek_at_next_line(li));
    if (input_format != 0 && input_format != 1) {
        fprintf(stderr, "taffy gerp: input must be MAF or TAF\n");
        return 1;
    }
    bool rle = false;
    Tag *header_tag = (input_format == 0)
                      ? taf_read_header_2(li, &rle)
                      : maf_read_header(li);

    // Tree: -t file overrides the `# hal` header tree.
    char *newick = NULL;
    if (treeFile != NULL) {
        FILE *tf = fopen(treeFile, "r");
        if (tf == NULL) {
            fprintf(stderr, "taffy gerp: cannot open tree file: %s\n", treeFile);
            return 1;
        }
        // Slurp the file.
        fseek(tf, 0, SEEK_END);
        long n = ftell(tf);
        fseek(tf, 0, SEEK_SET);
        newick = st_malloc((size_t)n + 1);
        size_t got = fread(newick, 1, (size_t)n, tf);
        newick[got] = '\0';
        fclose(tf);
    } else {
        Tag *hal = tag_find(header_tag, (char *) TAF_HAL_TREE_KEY);
        if (hal == NULL) {
            fprintf(stderr, "taffy gerp: no -t tree given and input header has no `# hal` tree.\n"
                            "  Either pass -t <tree.nwk> or supply input from cactus-hal2maf which preserves\n"
                            "  the tree as a `# hal` comment.\n");
            return 1;
        }
        newick = stString_copy(hal->value);
    }
    GerpTree *gt = gerp_tree_construct(newick);
    free(newick);
    if (gt == NULL) {
        fprintf(stderr, "taffy gerp: failed to parse Newick tree\n");
        return 1;
    }
    st_logInfo("taffy gerp: tree has %" PRIi64 " leaves\n", gerp_tree_n_leaves(gt));

    // Output(s).
    FILE *out_fh = (outputFile == NULL) ? stdout : fopen(outputFile, "w");
    if (out_fh == NULL) {
        fprintf(stderr, "taffy gerp: cannot open output file: %s\n", outputFile);
        return 1;
    }
    LW *out  = LW_construct(out_fh, use_compression);
    LW *dout = NULL;
    FILE *dout_fh = NULL;
    if (depthFile != NULL) {
        dout_fh = fopen(depthFile, "w");
        if (dout_fh == NULL) {
            fprintf(stderr, "taffy gerp: cannot open depth file: %s\n", depthFile);
            return 1;
        }
        dout = LW_construct(dout_fh, use_compression);
    }

    // Per-block scratch.
    GerpScratch *sc = gerp_scratch_construct(gt);
    int64_t n_leaves = gerp_tree_n_leaves(gt);
    Alignment_Row **row_by_leaf = st_malloc((size_t)n_leaves * sizeof(Alignment_Row *));
    char *leaf_bases = st_malloc((size_t)n_leaves);

    int64_t n_blocks = 0, n_skipped_paralog = 0, n_scored_cols = 0;

    Alignment *aln = NULL, *p_aln = NULL;
    while (1) {
        if (input_format == 0) aln = taf_read_block(p_aln, rle, li);
        else                   aln = maf_read_block(li);
        if (aln == NULL) break;
        n_blocks++;

        const char *unknown = NULL;
        int rc = gerp_block_resolve_rows(gt, aln, skip_paralogs, row_by_leaf, &unknown);
        if (rc == GERP_BLOCK_UNKNOWN_SPECIES) {
            fprintf(stderr, "taffy gerp: row in block has species not in tree: %s\n"
                            "  (pass -t with a tree that covers all leaves in the alignment, or\n"
                            "   drop the offending rows upstream)\n",
                    unknown ? unknown : "(unknown)");
            return 1;
        }
        if (rc == GERP_BLOCK_SKIP) {
            n_skipped_paralog++;
            goto next_block;
        }

        // Emit per-block variableStep header keyed on row-0's seq name.
        Alignment_Row *ref = aln->row;
        if (ref == NULL || ref->bases == NULL) goto next_block;
        // Reverse-strand row-0 means our anchor++ arithmetic would walk the
        // wrong direction and emit wig positions in the wrong coordinate
        // frame.  Universal MAF anchors are positive-strand by construction
        // (ancestor coords) and hg38-MAFs likewise; if a future input puts
        // a '-' on row-0 we'd need alignment_reorient_to_row first.  Fail
        // loud rather than silently produce garbage.
        if (!ref->strand) {
            fprintf(stderr, "taffy gerp: row-0 is on the reverse strand in a block "
                            "(%s).  Re-orient with taffy view -U query (or upstream) "
                            "before scoring.\n", ref->sequence_name);
            return 1;
        }
        emit_wig_header(out, ref->sequence_name);
        if (dout) emit_wig_header(dout, ref->sequence_name);

        int64_t anchor = ref->start;  // 0-based half-open; wig is 1-based, +1 below
        int64_t col_n  = aln->column_number;
        for (int64_t col = 0; col < col_n; col++) {
            char rb = ref->bases[col];
            if (rb == '-') continue;  // gap in row-0 -> no anchor coord
            // Build leaf_bases for this column.
            for (int64_t i = 0; i < n_leaves; i++) {
                Alignment_Row *r = row_by_leaf[i];
                leaf_bases[i] = (r != NULL) ? r->bases[col] : 0;
            }
            double rs = 0;
            int64_t depth = 0;
            bool scored = gerp_score_column(gt, sc, leaf_bases, min_leaves,
                                            branch_scale, &rs, &depth);
            int64_t wig_pos = anchor + 1;  // wig is 1-based
            if (scored) {
                emit_wig_value(out, wig_pos, rs);
                if (dout) emit_wig_value(dout, wig_pos, (double)depth);
                n_scored_cols++;
            } else if (dout) {
                // Depth track still records zero-depth (no score) columns
                // when the user asked for depth -- useful for masking.
                emit_wig_value(dout, wig_pos, (double)depth);
            }
            anchor++;
        }

      next_block:
        if (p_aln != NULL) alignment_destruct(p_aln, 1);
        p_aln = aln;
    }
    if (p_aln != NULL) alignment_destruct(p_aln, 1);

    st_logInfo("taffy gerp: %" PRIi64 " blocks read, %" PRIi64 " skipped (paralogs), "
               "%" PRIi64 " columns scored in %" PRIi64 " seconds\n",
               n_blocks, n_skipped_paralog, n_scored_cols,
               (int64_t)(time(NULL) - startTime));

    free(leaf_bases);
    free(row_by_leaf);
    gerp_scratch_destruct(sc);
    gerp_tree_destruct(gt);
    tag_destruct(header_tag);
    LI_destruct(li);
    if (inputFile != NULL) fclose(in_fh);
    LW_destruct(out, false);
    if (outputFile != NULL) fclose(out_fh);
    if (dout != NULL) {
        LW_destruct(dout, false);
        if (depthFile != NULL) fclose(dout_fh);
    }
    return 0;
}

/*
 * taffy gerp: per-column GERP RS conservation scoring -> wig output.
 *
 * The reference (row 0) of each MAF/TAF block keys the wig record
 * (sequence_name + advancing position over non-gap row-0 bases).  In a
 * normal hg38-anchored MAF that's hg38; in a universal MAF row 0 varies
 * per block (one of many ancestor chroms) and the wig ends up with
 * multiple chroms naturally.
 *
 * Parallelism: batched parallel-for over blocks.  A serial reader fills a
 * batch of N blocks, an OpenMP `parallel for` scores them with per-thread
 * scratch + per-block output buffers, and a serial writer drains the
 * buffers into the LW(s) in batch order.  Wig line order doesn't matter
 * for wigToBigWig, but we still emit in block-read order to keep
 * deterministic output across -T values.
 *
 * See taffy/inc/gerp.h for the algorithm contract.
 */

#include "taf.h"
#include "tai.h"
#include "tui.h"
#include "gerp.h"
#include "sonLib.h"
#include <unistd.h>
#include <ctype.h>
#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////////////////////////////////////////////////////////
// GerpBuf: a minimal growable byte buffer.  One per in-flight block so the
// worker can format wig output without touching the shared LW.  The writer
// then memcpy's the bytes into LW in batch order.
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    char  *buf;
    size_t len;
    size_t cap;
} GerpBuf;

static void gerpbuf_init(GerpBuf *b, size_t initial_cap) {
    b->cap = initial_cap;
    b->buf = st_malloc(b->cap);
    b->len = 0;
}
static void gerpbuf_destroy(GerpBuf *b) {
    free(b->buf);
    b->buf = NULL;
    b->len = b->cap = 0;
}
static void gerpbuf_reset(GerpBuf *b) { b->len = 0; }

static inline void gerpbuf_reserve(GerpBuf *b, size_t n) {
    if (b->len + n <= b->cap) return;
    while (b->len + n > b->cap) b->cap *= 2;
    b->buf = st_realloc(b->buf, b->cap);
}
static inline void gerpbuf_putn(GerpBuf *b, const char *s, size_t n) {
    gerpbuf_reserve(b, n);
    memcpy(b->buf + b->len, s, n);
    b->len += n;
}
static inline void gerpbuf_puts(GerpBuf *b, const char *s) {
    gerpbuf_putn(b, s, strlen(s));
}
static inline void gerpbuf_putc(GerpBuf *b, char c) {
    gerpbuf_reserve(b, 1);
    b->buf[b->len++] = c;
}
// Position + score line: "<pos> <value>\n".  snprintf is fine here since
// the line count per block is small; hot-path inlining would buy < 1%.
static void gerpbuf_put_score(GerpBuf *b, int64_t pos, double value, int prec) {
    char tmp[64];
    int n = snprintf(tmp, sizeof(tmp), "%" PRIi64 " %.*f\n", pos, prec, value);
    if (n > 0) gerpbuf_putn(b, tmp, (size_t)n);
}

/////////////////////////////////////////////////////////////////////////////
// Per-block scoring (worker-side).  All state passed in; no globals.  This
// function is invoked from inside an OpenMP parallel region.
/////////////////////////////////////////////////////////////////////////////

// Per-block result -- one slot per batch position.  Bufs are reused
// across batches (reset, not realloc'd, between batches).
typedef struct {
    GerpBuf rs;
    GerpBuf depth;
    int64_t cols_scored;
    int     status;              // GERP_BLOCK_OK / SKIP / UNKNOWN_SPECIES
    const char *unknown_seq;     // borrowed from aln when status == UNKNOWN
    bool    bad_strand;
    const char *bad_strand_seq;
} GerpBlockResult;

// Per-thread scratch -- one set per OpenMP worker thread.
typedef struct {
    GerpScratch    *sc;
    Alignment_Row **row_by_leaf;
    char           *leaf_bases;
} GerpThreadState;

static void score_one_block(const GerpTree *gt, GerpThreadState *ts,
                            const Alignment *aln, GerpBlockResult *res,
                            bool skip_paralogs, int64_t min_leaves,
                            double branch_scale, bool want_depth) {
    gerpbuf_reset(&res->rs);
    if (want_depth) gerpbuf_reset(&res->depth);
    res->cols_scored = 0;
    res->status = GERP_BLOCK_OK;
    res->unknown_seq = NULL;
    res->bad_strand = false;
    res->bad_strand_seq = NULL;

    int rc = gerp_block_resolve_rows(gt, aln, skip_paralogs,
                                     ts->row_by_leaf, &res->unknown_seq);
    if (rc == GERP_BLOCK_UNKNOWN_SPECIES) {
        res->status = GERP_BLOCK_UNKNOWN_SPECIES;
        return;
    }
    if (rc == GERP_BLOCK_SKIP) {
        res->status = GERP_BLOCK_SKIP;
        return;
    }

    Alignment_Row *ref = aln->row;
    if (ref == NULL || ref->bases == NULL) return;
    if (!ref->strand) {
        res->bad_strand = true;
        res->bad_strand_seq = ref->sequence_name;
        return;
    }

    // Emit per-block variableStep header keyed on row-0's seq name.
    gerpbuf_putn(&res->rs, "variableStep chrom=", 19);
    gerpbuf_puts(&res->rs, ref->sequence_name);
    gerpbuf_putc(&res->rs, '\n');
    if (want_depth) {
        gerpbuf_putn(&res->depth, "variableStep chrom=", 19);
        gerpbuf_puts(&res->depth, ref->sequence_name);
        gerpbuf_putc(&res->depth, '\n');
    }

    int64_t n_leaves = gerp_tree_n_leaves(gt);
    int64_t anchor   = ref->start;  // 0-based; wig is 1-based -> +1 below
    int64_t col_n    = aln->column_number;
    for (int64_t col = 0; col < col_n; col++) {
        char rb = ref->bases[col];
        if (rb == '-') continue;  // gap in row-0 -> no anchor coord
        for (int64_t i = 0; i < n_leaves; i++) {
            Alignment_Row *r = ts->row_by_leaf[i];
            ts->leaf_bases[i] = (r != NULL) ? r->bases[col] : 0;
        }
        double  rs    = 0;
        int64_t depth = 0;
        bool    scored = gerp_score_column(gt, ts->sc, ts->leaf_bases,
                                           min_leaves, branch_scale,
                                           &rs, &depth);
        int64_t wig_pos = anchor + 1;
        if (scored) {
            gerpbuf_put_score(&res->rs, wig_pos, rs, 4);
            if (want_depth) gerpbuf_put_score(&res->depth, wig_pos,
                                              (double)depth, 0);
            res->cols_scored++;
        } else if (want_depth) {
            gerpbuf_put_score(&res->depth, wig_pos, (double)depth, 0);
        }
        anchor++;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Driver
/////////////////////////////////////////////////////////////////////////////

static void usage(void) {
    fprintf(stderr, "taffy gerp [options]\n");
    fprintf(stderr, "Compute per-column GERP RS conservation scores -> wig.\n");
    fprintf(stderr, "-i --inputFile : Input MAF or TAF (auto-detected; bgzipped ok).  Reads stdin if omitted.\n");
    fprintf(stderr, "-o --outputFile : Output wig.  Writes stdout if omitted.\n");
    fprintf(stderr, "-c --useCompression : Bgzip the output (matches taffy view).\n");
    fprintf(stderr, "-t --tree : Newick tree override.  Default: the `# hal` tree comment in the input header.\n");
    fprintf(stderr, "-D --depthFile : Optional second wig with the per-base count of A/C/G/T leaves at each column.\n");
    fprintf(stderr, "-r --region SEQ:START-END : Restrict to one anchor chrom range.  Requires <inputFile>.tui (universal MAF, auto-detected) OR <inputFile>.tai (regular MAF).  Half-open, 0-based.\n");
    fprintf(stderr, "--branchScale : Global multiplier on branch lengths (default 1.0).\n");
    fprintf(stderr, "--minLeaves : Minimum surviving non-gap leaves to score a column (default 2).\n");
    fprintf(stderr, "--skipParalogs : Skip blocks with duplicate species (default ON).\n");
    fprintf(stderr, "--keepParalogs : Inverse of --skipParalogs (score using the first-seen row per species).\n");
    fprintf(stderr, "-T --threads N : Parallel block scoring + bgzf I/O on bgzipped streams (default 1).\n");
    fprintf(stderr, "-l --logLevel : Set the log level.\n");
    fprintf(stderr, "-h --help : Print this help message.\n");
}

enum {
    OPT_BRANCH_SCALE = 256,
    OPT_MIN_LEAVES,
    OPT_SKIP_PARALOGS,
    OPT_KEEP_PARALOGS,
};

int taf_gerp_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    char *logLevelString = NULL;
    char *inputFile      = NULL;
    char *outputFile     = NULL;
    char *depthFile      = NULL;
    char *treeFile       = NULL;
    char *region         = NULL;
    bool use_compression = false;
    double branch_scale  = 1.0;
    int64_t min_leaves   = 2;
    bool skip_paralogs   = true;
    int n_threads        = 1;

    while (1) {
        static struct option long_options[] = {
            { "logLevel",       required_argument, 0, 'l' },
            { "inputFile",      required_argument, 0, 'i' },
            { "outputFile",     required_argument, 0, 'o' },
            { "useCompression", no_argument,       0, 'c' },
            { "tree",           required_argument, 0, 't' },
            { "depthFile",      required_argument, 0, 'D' },
            { "region",         required_argument, 0, 'r' },
            { "branchScale",    required_argument, 0, OPT_BRANCH_SCALE },
            { "minLeaves",      required_argument, 0, OPT_MIN_LEAVES },
            { "skipParalogs",   no_argument,       0, OPT_SKIP_PARALOGS },
            { "keepParalogs",   no_argument,       0, OPT_KEEP_PARALOGS },
            { "threads",        required_argument, 0, 'T' },
            { "help",           no_argument,       0, 'h' },
            { 0, 0, 0, 0 }
        };
        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:ct:D:r:T:h", long_options, &option_index);
        if (key == -1) break;
        switch (key) {
            case 'l': logLevelString = optarg; break;
            case 'i': inputFile      = optarg; break;
            case 'o': outputFile     = optarg; break;
            case 'c': use_compression = true;  break;
            case 't': treeFile       = optarg; break;
            case 'D': depthFile      = optarg; break;
            case 'r': region         = optarg; break;
            case 'T': n_threads      = atoi(optarg); break;
            case OPT_BRANCH_SCALE:  branch_scale  = atof(optarg); break;
            case OPT_MIN_LEAVES:    min_leaves    = atol(optarg); break;
            case OPT_SKIP_PARALOGS: skip_paralogs = true;  break;
            case OPT_KEEP_PARALOGS: skip_paralogs = false; break;
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }
    if (n_threads < 1) n_threads = 1;

    st_setLogLevelFromString(logLevelString);
    LI_set_bgzf_threads(n_threads);

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
        fseek(tf, 0, SEEK_END);
        long n = ftell(tf);
        if (n < 0) {
            fprintf(stderr, "taffy gerp: ftell failed on tree file: %s\n", treeFile);
            fclose(tf);
            return 1;
        }
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
    bool want_depth = (depthFile != NULL);
    if (want_depth) {
        dout_fh = fopen(depthFile, "w");
        if (dout_fh == NULL) {
            fprintf(stderr, "taffy gerp: cannot open depth file: %s\n", depthFile);
            return 1;
        }
        dout = LW_construct(dout_fh, use_compression);
    }

    // Per-thread state.  One GerpScratch + row_by_leaf + leaf_bases per
    // worker.  Allocated once before the batched loop; reused across all
    // blocks the thread processes.
    int64_t n_leaves = gerp_tree_n_leaves(gt);
    GerpThreadState *ts = st_calloc((size_t)n_threads, sizeof(GerpThreadState));
    for (int t = 0; t < n_threads; t++) {
        ts[t].sc          = gerp_scratch_construct(gt);
        ts[t].row_by_leaf = st_malloc((size_t)n_leaves * sizeof(Alignment_Row *));
        ts[t].leaf_bases  = st_malloc((size_t)n_leaves);
    }

    // Per-batch slot.  4x n_threads keeps workers fed when block work is
    // uneven (heavy block ties up one worker; others can grind through
    // their next slot).
    int batch_cap = 4 * n_threads;
    if (batch_cap < 4) batch_cap = 4;
    Alignment       **batch_aln = st_calloc((size_t)batch_cap, sizeof(Alignment *));
    GerpBlockResult  *results   = st_calloc((size_t)batch_cap, sizeof(GerpBlockResult));
    for (int i = 0; i < batch_cap; i++) {
        gerpbuf_init(&results[i].rs, 4096);
        if (want_depth) gerpbuf_init(&results[i].depth, 4096);
    }

    int64_t n_blocks = 0, n_skipped_paralog = 0, n_scored_cols = 0;
    int fatal = 0;

    // Region mode: open the .tui (universal MAF) or .tai (regular MAF)
    // index and build a block iterator.  When set, the read loop pulls
    // blocks from the iterator instead of streaming the whole input -- so
    // the carry_aln chain is unused.
    Tui          *tui    = NULL;
    TuiExtractIt *tui_it = NULL;
    Tai          *tai    = NULL;
    FILE         *tai_fh = NULL;
    TaiIt        *tai_it = NULL;
    TuiInterval  *uiv    = NULL;
    char         *region_seq = NULL;
    int64_t       region_start = 0, region_length = 0;
    if (region != NULL) {
        if (inputFile == NULL) {
            fprintf(stderr, "taffy gerp: -r requires -i <file> (cannot index stdin)\n");
            return 1;
        }
        region_seq = tai_parse_region(region, &region_start, &region_length);
        if (region_seq == NULL) {
            fprintf(stderr, "taffy gerp: invalid region: %s\n", region);
            return 1;
        }
        char *tui_fn = tui_path(inputFile);
        bool tui_present = (access(tui_fn, F_OK) == 0);
        if (tui_present) {
            tui = tui_load(tui_fn);
            if (tui == NULL) {
                fprintf(stderr, "taffy gerp: cannot open .tui: %s\n", tui_fn);
                free(tui_fn); free(region_seq);
                return 1;
            }
            free(tui_fn);
            int64_t n_uiv = 0;
            uiv = tui_query(tui, region_seq, region_start,
                            region_start + region_length, &n_uiv);
            tui_it = tui_extract_iterator(tui, li, input_format == 1,
                                           rle, uiv, n_uiv);
            st_logInfo("taffy gerp: region %s:%" PRIi64 "-%" PRIi64
                       " resolved via .tui to %" PRIi64 " universal intervals\n",
                       region_seq, region_start, region_start + region_length, n_uiv);
        } else {
            free(tui_fn);
            char *tai_fn = tai_path(inputFile);
            tai_fh = fopen(tai_fn, "r");
            if (tai_fh == NULL) {
                fprintf(stderr, "taffy gerp: cannot open index (.tui or .tai must exist for %s): %s\n",
                        inputFile, tai_fn);
                free(tai_fn); free(region_seq);
                return 1;
            }
            tai = tai_load(tai_fh, input_format == 1);
            tai_it = tai_iterator(tai, li, rle, region_seq, region_start, region_length);
            st_logInfo("taffy gerp: region via .tai\n");
            free(tai_fn);
        }
    }

    // TAF read needs the previous block (p_aln) for delta-coord decoding.
    // We carry it across batches so the first block of batch N+1 has
    // batch N's last block as its predecessor.  MAF reads ignore p_aln.
    // Unused in region mode -- the iterator manages its own state.
    Alignment *carry_aln = NULL;

    while (!fatal) {
        // Phase A: serial read of up to batch_cap blocks.  TAF chains
        // through p_aln; once we've passed a block to taf_read_block as
        // p_aln, we can free the prior carry (its successor is now in
        // hand).  Region mode forces batch_cap=1: tui_extract_next's
        // returned alignment is invalidated by the NEXT call to it, so we
        // can only hold one at a time.  bgzf_mt still parallelises
        // decompress; cross-region parallelism (multiple shards) covers
        // any lost intra-region OMP throughput.
        int eff_batch_cap = (tui_it != NULL || tai_it != NULL) ? 1 : batch_cap;
        int n_read = 0;
        Alignment *p_for_read = carry_aln;
        while (n_read < eff_batch_cap) {
            Alignment *a = NULL;
            if (tui_it != NULL) {
                a = tui_extract_next(tui_it, li);
            } else if (tai_it != NULL) {
                a = tai_next(tai_it, li);
            } else if (input_format == 0) {
                a = taf_read_block(p_for_read, rle, li);
            } else {
                a = maf_read_block(li);
            }
            if (a == NULL) break;
            // The previous batch's carry is now consumed by taf_read_block
            // (it was used to decode this new block); safe to free.  Region
            // mode doesn't use the carry chain.
            if (tui_it == NULL && tai_it == NULL &&
                p_for_read == carry_aln && carry_aln != NULL) {
                alignment_destruct(carry_aln, 1);
                carry_aln = NULL;
            }
            batch_aln[n_read++] = a;
            p_for_read = a;
        }
        if (n_read == 0) {
            // No more blocks.  Drop any leftover carry.
            if (carry_aln != NULL) { alignment_destruct(carry_aln, 1); carry_aln = NULL; }
            break;
        }
        n_blocks += n_read;

        // Phase B: parallel score.  Each worker uses its own
        // GerpThreadState; results[i] is written by exactly one worker.
        #pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
        for (int i = 0; i < n_read; i++) {
#ifdef _OPENMP
            int t = omp_get_thread_num();
#else
            int t = 0;
#endif
            score_one_block(gt, &ts[t], batch_aln[i], &results[i],
                            skip_paralogs, min_leaves, branch_scale,
                            want_depth);
        }

        // Phase C: serial emit + accounting in batch order.
        for (int i = 0; i < n_read; i++) {
            GerpBlockResult *r = &results[i];
            if (r->status == GERP_BLOCK_UNKNOWN_SPECIES) {
                fprintf(stderr, "taffy gerp: row in block has species not in tree: %s\n"
                                "  (pass -t with a tree that covers all leaves in the alignment, or\n"
                                "   drop the offending rows upstream)\n",
                        r->unknown_seq ? r->unknown_seq : "(unknown)");
                fatal = 1;
                break;
            }
            if (r->bad_strand) {
                fprintf(stderr, "taffy gerp: row-0 is on the reverse strand in a block "
                                "(%s).  Re-orient with taffy view -U query (or upstream) "
                                "before scoring.\n",
                        r->bad_strand_seq ? r->bad_strand_seq : "(unknown)");
                fatal = 1;
                break;
            }
            if (r->status == GERP_BLOCK_SKIP) {
                n_skipped_paralog++;
                continue;
            }
            n_scored_cols += r->cols_scored;
            if (r->rs.len > 0)            LW_putn(out,  r->rs.buf,    r->rs.len);
            if (want_depth && r->depth.len > 0) LW_putn(dout, r->depth.buf, r->depth.len);
        }

        // Phase D: free this batch's alignments, retaining the last one
        // as carry for next batch's TAF chain.  On a fatal error we still
        // free everything cleanly.  Region mode: tui_extract_iterator
        // owns the yielded alignment (it frees on next call) so we
        // null-out without destructing; tai_iterator yields ownership so
        // we destruct normally.
        if (tui_it != NULL) {
            for (int i = 0; i < n_read; i++) batch_aln[i] = NULL;
        } else {
            int free_to = (fatal || tai_it != NULL) ? n_read : n_read - 1;
            for (int i = 0; i < free_to; i++) {
                if (batch_aln[i] != NULL) {
                    alignment_destruct(batch_aln[i], 1);
                    batch_aln[i] = NULL;
                }
            }
            if (!fatal && tai_it == NULL && n_read > 0) {
                carry_aln = batch_aln[n_read - 1];
                batch_aln[n_read - 1] = NULL;
            }
        }
    }
    if (carry_aln != NULL) alignment_destruct(carry_aln, 1);

    if (tui_it != NULL) tui_extract_iterator_destruct(tui_it);
    if (tui    != NULL) tui_destruct(tui);
    if (tai_it != NULL) tai_iterator_destruct(tai_it);
    if (tai    != NULL) tai_destruct(tai);
    if (tai_fh != NULL) fclose(tai_fh);
    free(uiv);
    free(region_seq);

    st_logInfo("taffy gerp: %" PRIi64 " blocks read, %" PRIi64 " skipped (paralogs), "
               "%" PRIi64 " columns scored in %" PRIi64 " seconds\n",
               n_blocks, n_skipped_paralog, n_scored_cols,
               (int64_t)(time(NULL) - startTime));

    for (int i = 0; i < batch_cap; i++) {
        gerpbuf_destroy(&results[i].rs);
        if (want_depth) gerpbuf_destroy(&results[i].depth);
    }
    free(results);
    free(batch_aln);
    for (int t = 0; t < n_threads; t++) {
        gerp_scratch_destruct(ts[t].sc);
        free(ts[t].row_by_leaf);
        free(ts[t].leaf_bases);
    }
    free(ts);

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
    return fatal;
}

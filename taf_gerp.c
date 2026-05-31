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
    bool    had_paralog;         // any leaf had > 1 row in this block
    const char *unknown_seq;     // borrowed from aln when status == UNKNOWN
    bool    bad_strand;
    const char *bad_strand_seq;
} GerpBlockResult;

// Per-thread scratch -- one set per OpenMP worker thread.  entries[] is
// grown on demand to fit the largest block this thread sees; leaf_csets
// is fixed-size at n_leaves bytes.
typedef struct {
    GerpScratch  *sc;
    GerpRowEntry *entries;
    int64_t       entries_cap;
    uint8_t      *leaf_csets;
} GerpThreadState;

static void score_one_block(const GerpTree *gt, GerpThreadState *ts,
                            const Alignment *aln, GerpBlockResult *res,
                            GerpParalogPolicy policy, int64_t min_leaves,
                            double branch_scale, bool want_depth) {
    gerpbuf_reset(&res->rs);
    if (want_depth) gerpbuf_reset(&res->depth);
    res->cols_scored = 0;
    res->status = GERP_BLOCK_OK;
    res->had_paralog = false;
    res->unknown_seq = NULL;
    res->bad_strand = false;
    res->bad_strand_seq = NULL;

    // entries[] needs room for every leaf row in the block (UNION mode
    // keeps them all).  Grow geometrically and remember the new cap.
    if (aln->row_number > ts->entries_cap) {
        int64_t new_cap = ts->entries_cap > 0 ? ts->entries_cap : 16;
        while (new_cap < aln->row_number) new_cap *= 2;
        ts->entries = st_realloc(ts->entries,
                                 (size_t)new_cap * sizeof(GerpRowEntry));
        ts->entries_cap = new_cap;
    }

    int64_t n_active = 0, n_paralog_dups = 0;
    int rc = gerp_block_resolve_rows(gt, aln, policy,
                                     ts->entries, ts->entries_cap,
                                     &n_active, &n_paralog_dups,
                                     &res->unknown_seq);
    if (rc == GERP_BLOCK_UNKNOWN_SPECIES) {
        res->status = GERP_BLOCK_UNKNOWN_SPECIES;
        return;
    }
    if (rc == GERP_BLOCK_SKIP) {
        res->status      = GERP_BLOCK_SKIP;
        res->had_paralog = true;  // SKIP only triggers on paralog
        return;
    }
    res->had_paralog = (n_paralog_dups > 0);

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
        // Build per-leaf character set by OR-ing each active row's base
        // bit at this column.  In UNION mode paralog rows of the same
        // species accumulate into a multi-bit cset; in FIRST mode there's
        // at most one entry per leaf so the cset is single-bit.
        memset(ts->leaf_csets, 0, (size_t)n_leaves);
        for (int64_t k = 0; k < n_active; k++) {
            const GerpRowEntry *e = &ts->entries[k];
            uint8_t bit = gerp_base_to_bit(e->row->bases[col]);
            if (bit) ts->leaf_csets[e->leaf_id] |= bit;
        }
        double  rs    = 0;
        int64_t depth = 0;
        bool    scored = gerp_score_column_csets(gt, ts->sc, ts->leaf_csets,
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
    fprintf(stderr, "-R --regionFile FILE : Like -r but a file of regions (one SEQ:START-END per line; '#' comments + blanks ignored).  Amortises one .tui/.tai load across all regions -- intended for SLURM-style sharding with many small regions per shard.\n");
    fprintf(stderr, "   --columnRange LO-HI : Restrict to a universal-column range (half-open, 0-based; HI <= T from `taffy stats -u`).  Requires .tui.  Each universal column belongs to exactly one shard, so this is the natural unit for SLURM sharding -- no chrom-by-chrom double-counting.  Mutually exclusive with -r/-R.\n");
    fprintf(stderr, "--branchScale : Global multiplier on branch lengths (default 1.0).\n");
    fprintf(stderr, "--minLeaves : Minimum surviving non-gap leaves to score a column (default 2).\n");
    fprintf(stderr, "--paralog MODE : How to treat duplicate-species rows in one block (default union).\n");
    fprintf(stderr, "                 union -- OR each species's paralog bases into a multi-state leaf cset (Hartigan).\n");
    fprintf(stderr, "                 skip  -- drop the entire block (strict GERP++ semantics).\n");
    fprintf(stderr, "                 first -- score using only the first-seen row per species.\n");
    fprintf(stderr, "--skipParalogs : Alias for --paralog skip.\n");
    fprintf(stderr, "--keepParalogs : Alias for --paralog first (kept for back-compat; prefer --paralog union).\n");
    fprintf(stderr, "-T --threads N : Parallel block scoring + bgzf I/O on bgzipped streams (default 1).\n");
    fprintf(stderr, "-l --logLevel : Set the log level.\n");
    fprintf(stderr, "-h --help : Print this help message.\n");
}

enum {
    OPT_BRANCH_SCALE = 256,
    OPT_MIN_LEAVES,
    OPT_PARALOG,
    OPT_SKIP_PARALOGS,
    OPT_KEEP_PARALOGS,
    OPT_COLUMN_RANGE,
};

static int parse_paralog_policy(const char *s, GerpParalogPolicy *out) {
    if (s == NULL) return -1;
    if (strcmp(s, "union") == 0) { *out = GERP_PARALOG_UNION; return 0; }
    if (strcmp(s, "skip")  == 0) { *out = GERP_PARALOG_SKIP;  return 0; }
    if (strcmp(s, "first") == 0) { *out = GERP_PARALOG_FIRST; return 0; }
    return -1;
}

// One parsed region for the outer loop in region mode.  `seq` is owned
// (allocated by tai_parse_region).
typedef struct {
    char   *seq;
    int64_t start;
    int64_t length;
} GerpRegion;

// Read a regions file -- one "SEQ:START-END" per line; '#' and blank
// lines are ignored.  Returns a freshly-allocated array of size *n_out,
// or NULL on error (after printing the offending line/path to stderr).
// Caller frees each entry's seq and the array.
static GerpRegion *read_regions_file(const char *path, int64_t *n_out) {
    FILE *f = fopen(path, "r");
    if (f == NULL) {
        fprintf(stderr, "taffy gerp: cannot open --regionFile: %s\n", path);
        return NULL;
    }
    char line[4096];
    int64_t n_lines = 0;
    int64_t n_regions = 0, out_cap = 0;
    GerpRegion *out = NULL;
    while (fgets(line, sizeof(line), f) != NULL) {
        n_lines++;
        size_t len = strlen(line);
        // Reject lines that didn't fit -- otherwise we'd silently parse
        // a truncated region.
        if (len == sizeof(line) - 1 && line[len-1] != '\n') {
            fprintf(stderr, "taffy gerp: --regionFile %s line %" PRIi64 ": "
                            "line exceeds %zu bytes\n",
                    path, n_lines, sizeof(line) - 1);
            for (int64_t i = 0; i < n_regions; i++) free(out[i].seq);
            free(out); fclose(f);
            return NULL;
        }
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' ||
                           line[len-1] == ' '  || line[len-1] == '\t')) {
            line[--len] = '\0';
        }
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0' || *p == '#') continue;

        int64_t start = 0, length = 0;
        char *seq = tai_parse_region(p, &start, &length);
        if (seq == NULL) {
            fprintf(stderr, "taffy gerp: --regionFile %s line %" PRIi64 ": "
                            "invalid region: %s\n", path, n_lines, p);
            for (int64_t i = 0; i < n_regions; i++) free(out[i].seq);
            free(out); fclose(f);
            return NULL;
        }
        if (n_regions >= out_cap) {
            int64_t new_cap = out_cap > 0 ? out_cap * 2 : 16;
            out = st_realloc(out, (size_t)new_cap * sizeof(GerpRegion));
            out_cap = new_cap;
        }
        out[n_regions].seq    = seq;
        out[n_regions].start  = start;
        out[n_regions].length = length;
        n_regions++;
    }
    fclose(f);
    *n_out = n_regions;
    return out;
}

int taf_gerp_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    char *logLevelString = NULL;
    char *inputFile      = NULL;
    char *outputFile     = NULL;
    char *depthFile      = NULL;
    char *treeFile       = NULL;
    char *region         = NULL;
    char *regionFile     = NULL;
    char *columnRangeArg = NULL;
    bool use_compression = false;
    double branch_scale  = 1.0;
    int64_t min_leaves   = 2;
    GerpParalogPolicy paralog_policy = GERP_PARALOG_UNION;
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
            { "regionFile",     required_argument, 0, 'R' },
            { "columnRange",    required_argument, 0, OPT_COLUMN_RANGE },
            { "branchScale",    required_argument, 0, OPT_BRANCH_SCALE },
            { "minLeaves",      required_argument, 0, OPT_MIN_LEAVES },
            { "paralog",        required_argument, 0, OPT_PARALOG },
            { "skipParalogs",   no_argument,       0, OPT_SKIP_PARALOGS },
            { "keepParalogs",   no_argument,       0, OPT_KEEP_PARALOGS },
            { "threads",        required_argument, 0, 'T' },
            { "help",           no_argument,       0, 'h' },
            { 0, 0, 0, 0 }
        };
        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:ct:D:r:R:T:h", long_options, &option_index);
        if (key == -1) break;
        switch (key) {
            case 'l': logLevelString = optarg; break;
            case 'i': inputFile      = optarg; break;
            case 'o': outputFile     = optarg; break;
            case 'c': use_compression = true;  break;
            case 't': treeFile       = optarg; break;
            case 'D': depthFile      = optarg; break;
            case 'r': region         = optarg; break;
            case 'R': regionFile     = optarg; break;
            case 'T': n_threads      = atoi(optarg); break;
            case OPT_BRANCH_SCALE:  branch_scale  = atof(optarg); break;
            case OPT_MIN_LEAVES:    min_leaves    = atol(optarg); break;
            case OPT_PARALOG:
                if (parse_paralog_policy(optarg, &paralog_policy) != 0) {
                    fprintf(stderr, "taffy gerp: --paralog must be one of union|skip|first (got %s)\n",
                            optarg);
                    return 1;
                }
                break;
            case OPT_SKIP_PARALOGS: paralog_policy = GERP_PARALOG_SKIP;  break;
            case OPT_KEEP_PARALOGS: paralog_policy = GERP_PARALOG_FIRST; break;
            case OPT_COLUMN_RANGE:  columnRangeArg = optarg; break;
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

    // CLI mutual-exclusion check goes BEFORE the output file opens so a
    // bad combination doesn't leave a stub output file behind.  -r, -R,
    // and --columnRange are pairwise exclusive (full-fledged validation
    // of the individual values happens later, after the index is loaded).
    if ((region != NULL) + (regionFile != NULL) + (columnRangeArg != NULL) > 1) {
        fprintf(stderr, "taffy gerp: -r, -R, and --columnRange are mutually exclusive\n");
        return 1;
    }

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

    // Per-thread state.  One GerpScratch + entries[] + leaf_csets per
    // worker.  entries[] grows on demand as blocks come in; leaf_csets is
    // n_leaves bytes (one 4-bit cset per leaf, per column).
    int64_t n_leaves = gerp_tree_n_leaves(gt);
    GerpThreadState *ts = st_calloc((size_t)n_threads, sizeof(GerpThreadState));
    for (int t = 0; t < n_threads; t++) {
        ts[t].sc          = gerp_scratch_construct(gt);
        ts[t].entries     = NULL;
        ts[t].entries_cap = 0;
        ts[t].leaf_csets  = st_malloc((size_t)n_leaves);
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

    int64_t n_blocks = 0, n_skipped_paralog = 0, n_paralog_blocks = 0;
    int64_t n_scored_cols = 0;
    int fatal = 0;

    // Build the regions list:
    //   no -r/-R  -> regions == NULL  (stream the whole input)
    //   -r        -> singleton list of size 1
    //   -R FILE   -> list parsed from file
    // Mutual exclusion already checked above the output opens.
    GerpRegion *regions = NULL;
    int64_t n_regions = 0;
    if (regionFile != NULL) {
        regions = read_regions_file(regionFile, &n_regions);
        if (regions == NULL) return 1;
        st_logInfo("taffy gerp: %" PRIi64 " regions read from %s\n",
                   n_regions, regionFile);
    } else if (region != NULL) {
        regions = st_calloc(1, sizeof(GerpRegion));
        regions[0].seq = tai_parse_region(region, &regions[0].start,
                                          &regions[0].length);
        if (regions[0].seq == NULL) {
            fprintf(stderr, "taffy gerp: invalid region: %s\n", region);
            free(regions);
            return 1;
        }
        n_regions = 1;
    }

    // Parse --columnRange LO-HI directly into a TuiInterval (the iterator
    // takes column intervals natively, no tui_query needed).  Half-open
    // [LO, HI), validated against T below once the .tui is loaded.
    TuiInterval col_range_iv = { 0, 0 };
    bool have_col_range = false;
    bool empty_col_range = false;
    if (columnRangeArg != NULL) {
        long long lo = 0, hi = 0;
        char extra = 0;
        if (sscanf(columnRangeArg, "%lld-%lld%c", &lo, &hi, &extra) != 2 ||
            lo < 0 || hi < lo) {
            fprintf(stderr, "taffy gerp: invalid --columnRange %s "
                            "(expected LO-HI with HI >= LO >= 0)\n", columnRangeArg);
            return 1;
        }
        col_range_iv.start = (int64_t)lo;
        col_range_iv.end   = (int64_t)hi;
        have_col_range = true;
        // lo == hi is a legal empty range: a SLURM shard whose slice of T
        // happens to be zero (T < N or rounding).  Skip all work cleanly
        // so the runner doesn't have to special-case it.
        if (lo == hi) empty_col_range = true;
    }

    // Region / column-range mode: load the .tui (universal MAF) or .tai
    // (regular MAF) index ONCE up front; the outer loop below reuses it
    // across every region.  --columnRange forces the .tui path (column
    // coords are a .tui concept).
    Tui  *tui    = NULL;
    Tai  *tai    = NULL;
    FILE *tai_fh = NULL;
    if (regions != NULL || have_col_range) {
        if (inputFile == NULL) {
            fprintf(stderr, "taffy gerp: -r/-R/--columnRange requires -i <file> (cannot index stdin)\n");
            return 1;
        }
        char *tui_fn = tui_path(inputFile);
        bool tui_present = (access(tui_fn, F_OK) == 0);
        if (have_col_range && !tui_present) {
            fprintf(stderr, "taffy gerp: --columnRange requires a .tui index: %s\n", tui_fn);
            free(tui_fn);
            return 1;
        }
        if (tui_present) {
            tui = tui_load(tui_fn);
            if (tui == NULL) {
                fprintf(stderr, "taffy gerp: cannot open .tui: %s\n", tui_fn);
                free(tui_fn);
                return 1;
            }
            free(tui_fn);
            st_logInfo("taffy gerp: loaded .tui index (universal MAF mode)\n");
        } else {
            free(tui_fn);
            char *tai_fn = tai_path(inputFile);
            tai_fh = fopen(tai_fn, "r");
            if (tai_fh == NULL) {
                fprintf(stderr, "taffy gerp: cannot open index (.tui or .tai must exist for %s): %s\n",
                        inputFile, tai_fn);
                free(tai_fn);
                return 1;
            }
            tai = tai_load(tai_fh, input_format == 1);
            free(tai_fn);
            st_logInfo("taffy gerp: loaded .tai index\n");
        }
        if (have_col_range) {
            int64_t T = tui_total_columns(tui);
            if (col_range_iv.end > T) {
                fprintf(stderr, "taffy gerp: --columnRange %" PRIi64 "-%" PRIi64
                                " exceeds T=%" PRIi64 " (total universal columns).\n",
                        col_range_iv.start, col_range_iv.end, T);
                return 1;
            }
            st_logInfo("taffy gerp: column-range mode, [%" PRIi64 ", %" PRIi64 ") of T=%" PRIi64 "\n",
                       col_range_iv.start, col_range_iv.end, T);
        }
    }

    // Outer loop: 1 iteration in streaming or column-range mode;
    // n_regions iterations in region mode (each builds its own iterator
    // over the shared tui/tai, then runs the inner read/score/emit batched
    // loop to exhaustion).  Empty column-range -> 0 iterations (success
    // no-op for SLURM shards that ended up with a zero slice).
    int64_t n_iter = (regions != NULL) ? n_regions
                   : empty_col_range   ? 0
                   :                     1;
    for (int64_t reg_idx = 0; reg_idx < n_iter && !fatal; reg_idx++) {
        TuiExtractIt *tui_it = NULL;
        TaiIt        *tai_it = NULL;
        TuiInterval  *uiv    = NULL;

        if (regions != NULL) {
            GerpRegion *r = &regions[reg_idx];
            if (tui != NULL) {
                int64_t n_uiv = 0;
                uiv = tui_query(tui, r->seq, r->start, r->start + r->length,
                                &n_uiv);
                tui_it = tui_extract_iterator(tui, li, input_format == 1,
                                              rle, uiv, n_uiv);
                st_logDebug("taffy gerp: region %s:%" PRIi64 "-%" PRIi64
                            " resolved via .tui to %" PRIi64 " universal intervals\n",
                            r->seq, r->start, r->start + r->length, n_uiv);
            } else {
                tai_it = tai_iterator(tai, li, rle, r->seq, r->start, r->length);
                st_logDebug("taffy gerp: region %s:%" PRIi64 "-%" PRIi64 " via .tai\n",
                            r->seq, r->start, r->start + r->length);
            }
        } else if (have_col_range) {
            // Direct column-range iterator: hand TuiInterval{LO,HI} straight
            // to tui_extract_iterator (no tui_query / chrom resolution).
            // `uiv` stays NULL -- col_range_iv lives on the stack across the
            // single iteration, no allocation to free in the per-iter cleanup.
            tui_it = tui_extract_iterator(tui, li, input_format == 1, rle,
                                          &col_range_iv, 1);
        }

        // TAF read needs the previous block (p_aln) for delta-coord
        // decoding.  We carry it across batches so the first block of
        // batch N+1 has batch N's last block as its predecessor.  MAF
        // reads ignore p_aln.  Unused in region mode -- the iterator
        // manages its own state.  Carry is per-region (one streaming
        // run inside each iterator).
        Alignment *carry_aln = NULL;

        while (!fatal) {
            // Phase A: serial read of up to batch_cap blocks.  TAF chains
            // through p_aln; once we've passed a block to taf_read_block as
            // p_aln, we can free the prior carry (its successor is now in
            // hand).  Region mode forces batch_cap=1: tui_extract_next's
            // returned alignment is invalidated by the NEXT call to it, so
            // we can only hold one at a time.  bgzf_mt still parallelises
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
                // Column-range exactly-once semantics are guaranteed by the
                // iterator itself: tui_extract_iterator clips each physical
                // block to the requested column range, so a boundary block
                // becomes two sub-blocks -- one in each shard, with the
                // universal columns of each sub-block disjoint from the
                // other.  No additional driver-side filter needed; verified
                // by split-vs-single equivalence test on apes (300K cols,
                // data lines byte-identical after sort).
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
                // No more blocks in this region (or in the stream).
                if (carry_aln != NULL) {
                    alignment_destruct(carry_aln, 1);
                    carry_aln = NULL;
                }
                break;
            }
            n_blocks += n_read;

            // Phase B: parallel score.
            #pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
            for (int i = 0; i < n_read; i++) {
#ifdef _OPENMP
                int t = omp_get_thread_num();
#else
                int t = 0;
#endif
                score_one_block(gt, &ts[t], batch_aln[i], &results[i],
                                paralog_policy, min_leaves, branch_scale,
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
                if (r->had_paralog) n_paralog_blocks++;
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
        if (tai_it != NULL) tai_iterator_destruct(tai_it);
        free(uiv);
    }

    if (tui    != NULL) tui_destruct(tui);
    if (tai    != NULL) tai_destruct(tai);
    if (tai_fh != NULL) fclose(tai_fh);
    if (regions != NULL) {
        for (int64_t i = 0; i < n_regions; i++) free(regions[i].seq);
        free(regions);
    }

    const char *policy_name = (paralog_policy == GERP_PARALOG_UNION) ? "union"
                            : (paralog_policy == GERP_PARALOG_SKIP)  ? "skip"
                            :                                          "first";
    st_logInfo("taffy gerp: %" PRIi64 " blocks read, %" PRIi64 " with paralogs (policy=%s; "
               "%" PRIi64 " block-skips), %" PRIi64 " columns scored in %" PRIi64 " seconds\n",
               n_blocks, n_paralog_blocks, policy_name, n_skipped_paralog,
               n_scored_cols, (int64_t)(time(NULL) - startTime));

    for (int i = 0; i < batch_cap; i++) {
        gerpbuf_destroy(&results[i].rs);
        if (want_depth) gerpbuf_destroy(&results[i].depth);
    }
    free(results);
    free(batch_aln);
    for (int t = 0; t < n_threads; t++) {
        gerp_scratch_destruct(ts[t].sc);
        free(ts[t].entries);
        free(ts[t].leaf_csets);
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

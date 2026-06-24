/*
 * taffy depth: per-column leaf depth and/or GERP RS conservation -> wig output.
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
    // --bin: this block's columns pre-aggregated per universal-column bin.
    // The block's columns are contiguous, so the touched bins are the
    // contiguous range [bin_first, bin_first+bin_n); bin_sum/bin_cnt hold the
    // depth sum and column count this block contributes to each.  Phase C
    // merges these across blocks (the boundary bin is shared with neighbours).
    int64_t  bin_first;
    int64_t *bin_sum;
    int64_t *bin_cnt;
    int64_t  bin_cap;
    int64_t  bin_n;
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

// The wig is keyed on the raw universal column [0,T), 2e9-CHUNKED into chroms
// uni<c/TUI_UNI_CHUNK> at position c%TUI_UNI_CHUNK (0-based for the --bin
// bedGraph, 1-based for the per-column wig; see TUI_UNI_CHUNK in tui.h; the
// --bin bedGraph matches the `taffy lift --bigwig` reader).  2e9 < 2^31 keeps
// every chunk a buildable signed-32-bit wigToBigWig chrom, so there is no T
// limit.  Requires the .tui (for T / column coords).

// Emit one `variableStep chrom=<chrom>` header to the rs and/or depth buffer.
static void gerp_emit_wig_header(GerpBlockResult *res, const char *chrom,
                                 bool to_rs, bool to_depth) {
    if (to_rs) {
        gerpbuf_putn(&res->rs, "variableStep chrom=", 19);
        gerpbuf_puts(&res->rs, chrom);
        gerpbuf_putc(&res->rs, '\n');
    }
    if (to_depth) {
        gerpbuf_putn(&res->depth, "variableStep chrom=", 19);
        gerpbuf_puts(&res->depth, chrom);
        gerpbuf_putc(&res->depth, '\n');
    }
}

// --bin integer coords: write one binned bedGraph line for the mean depth over
// the bin's columns.  The bin spans global universal columns [bin*N, bin*N+cnt);
// it is emitted 2e9-CHUNKED into chrom uni<chunk>, 0-based local position
// [start, start+cnt), where chunk = (bin*N)/TUI_UNI_CHUNK and start = (bin*N) %
// TUI_UNI_CHUNK.  Because TUI_UNI_CHUNK is a multiple of N (enforced at parse) no
// bin straddles a chunk boundary, so the whole bin lives in one chunk and
// start+cnt <= TUI_UNI_CHUNK < 2^31.  This matches the `taffy lift --bigwig`
// reader, which does chunk=col/TUI_UNI_CHUNK, base=chunk*TUI_UNI_CHUNK and treats
// the stored coords as 0-based.  The column axis is monotone, so this runs off
// the Phase-C running binner (already sorted, no external sort).  Drops
// all-unscored bins (awk parity).
static void gerp_flush_bin(LW *dout, int64_t bin, int64_t sum, int64_t cnt,
                           int64_t bin_size) {
    if (sum <= 0 || cnt <= 0) return;
    int64_t start  = bin * bin_size;                 // global universal column
    int64_t chunk  = start / TUI_UNI_CHUNK;          // 2e9-chunk index
    int64_t lstart = start - chunk * TUI_UNI_CHUNK;   // 0-based local position
    char line[96];
    int n = snprintf(line, sizeof line, "uni%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t%.4g\n",
                     chunk, lstart, lstart + cnt, (double)sum / (double)cnt);
    LW_putn(dout, line, (size_t)n);
}

static void score_one_block(const GerpTree *gt, GerpThreadState *ts,
                            const Alignment *aln, GerpBlockResult *res,
                            GerpParalogPolicy policy, int64_t min_leaves,
                            double branch_scale, bool want_depth,
                            bool depth_only, int64_t block_start_col,
                            int64_t bin_size) {
    gerpbuf_reset(&res->rs);
    if (want_depth) gerpbuf_reset(&res->depth);
    res->cols_scored = 0;
    res->bin_n = 0;
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

    // The `variableStep chrom=...` header is emitted LAZILY -- only when a
    // column is actually written -- so a block with no scored column (all
    // < min_leaves, common for ancestor-heavy universal-MAF blocks) writes
    // NOTHING rather than a dangling empty header (which wigToBigWig rejects).
    // The per-base integer wig keys on the chunked "uni<chunk>" chrom; one
    // block lies wholly within one chunk (blocks are << 2e9 columns), so a
    // single lazy header per block is correct.  rs and depth trigger on the
    // same column, so their (chrom,pos) structure stays byte-identical
    // (gerp-stats invariant).  No --rs leaves the rs buffer empty.
    bool hdr_emitted = false;

    int64_t n_leaves = gerp_tree_n_leaves(gt);
    int64_t col_n    = aln->column_number;
    for (int64_t col = 0; col < col_n; col++) {
        char rb = ref->bases[col];
        if (rb == '-') continue;  // gap in row-0 -> no universal column
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
        bool    scored;
        if (depth_only) {
            // Depth = number of leaves present (non-gap) at this column -- the
            // same count gerp_score_column_csets returns -- but WITHOUT the
            // O(n_nodes) Hartigan/RS tree walk, which we don't need.
            for (int64_t i = 0; i < n_leaves; i++)
                if (ts->leaf_csets[i]) depth++;
            scored = (depth >= min_leaves);
        } else {
            scored = gerp_score_column_csets(gt, ts->sc, ts->leaf_csets,
                                             min_leaves, branch_scale,
                                             &rs, &depth);
        }
        if (bin_size > 0) {
            // --bin: aggregate into this block's per-bin partials.  bin_sum
            // sums SCORED columns only (depth >= min_leaves); bin_cnt counts
            // every column the bin covers here.  Row-0 is gap-free in a
            // universal MAF, and --bin shards are N-aligned (checked at the
            // --columnRange parse), so each bin is processed from its start and
            // bin_cnt is its true universal-column count: N for a full bin,
            // fewer only for the genome's last (partial) bin.  flush emits
            // mean = bin_sum/bin_cnt (== the awk binner's sum/clamped-width)
            // and drops a bin with no scored column.  The block's columns are
            // contiguous, so the touched bins form one contiguous range
            // [bin_first, bin_first+bin_n).  bin is the GLOBAL universal-column
            // bin index over [0,T); the 2e9 chunk split happens only at emit
            // time in gerp_flush_bin.
            int64_t bin = (block_start_col + col) / bin_size;
            if (res->bin_n == 0) {
                res->bin_first  = bin;
            }
            int64_t idx = bin - res->bin_first;
            if (idx >= res->bin_n) {                   // start of the next bin
                if (idx >= res->bin_cap) {
                    int64_t nc = res->bin_cap ? res->bin_cap * 2 : 64;
                    while (nc <= idx) nc *= 2;
                    res->bin_sum = st_realloc(res->bin_sum, (size_t)nc * sizeof(int64_t));
                    res->bin_cnt = st_realloc(res->bin_cnt, (size_t)nc * sizeof(int64_t));
                    res->bin_cap = nc;
                }
                res->bin_sum[idx] = 0;
                res->bin_cnt[idx] = 0;
                res->bin_n        = idx + 1;
            }
            if (scored) { res->bin_sum[idx] += depth; res->cols_scored++; }
            res->bin_cnt[idx]++;   // every column -> bin-mean denominator
        } else if (scored) {
            // Integer universal-column wig, 2e9-CHUNKED: chrom uni<chunk>,
            // 1-based local position (gcol % TUI_UNI_CHUNK)+1.  A block is far
            // smaller than 2e9 columns, so it lies wholly within one chunk and
            // the single lazy header per block names the right chunk.
            int64_t gcol  = block_start_col + col;
            int64_t chunk = gcol / TUI_UNI_CHUNK;
            int64_t wig_pos = (gcol - chunk * TUI_UNI_CHUNK) + 1;  // 1-based local
            if (!hdr_emitted) {
                char uchrom[24];
                snprintf(uchrom, sizeof uchrom, "uni%" PRId64, chunk);
                gerp_emit_wig_header(res, uchrom, !depth_only, want_depth);
                hdr_emitted = true;
            }
            if (!depth_only) gerpbuf_put_score(&res->rs, wig_pos, rs, 4);
            if (want_depth)  gerpbuf_put_score(&res->depth, wig_pos,
                                               (double)depth, 0);
            res->cols_scored++;
        }
        // No depth-only emit when !scored: that asymmetry caused gerp-stats
        // to desync after whole-block depth-only runs (e.g. ancestor blocks
        // with < min_leaves surviving non-gap leaves -- RS emitted nothing
        // but depth emitted every column).  Per-column depth at unscored
        // sites is also not useful to the gerp-stats z-score normalisation,
        // which keys on (chrom, pos) tuples that exist in BOTH wigs.
    }
}

/////////////////////////////////////////////////////////////////////////////
// Driver
/////////////////////////////////////////////////////////////////////////////

static void usage(void) {
    fprintf(stderr, "taffy depth [options]\n");
    fprintf(stderr, "Per-column leaf DEPTH and/or GERP RS conservation from a MAF/TAF -> wig.\n");
    fprintf(stderr, "-i --inputFile : Input MAF or TAF (auto-detected; bgzipped ok).  Reads stdin if omitted.\n");
    fprintf(stderr, "Outputs -- at least one required, each opt-in:\n");
    fprintf(stderr, "  --depth FILE : Per-column (or --bin binned) count of non-gap leaves at each column.\n");
    fprintf(stderr, "  --rs FILE    : Per-column GERP RS conservation score (Hartigan tree walk).  Needs a tree.\n");
    fprintf(stderr, "Output is on the UNIVERSAL-COLUMN axis: the raw column [0,T) 2e9-chunked into chroms\n");
    fprintf(stderr, "  uni0,uni1,... at 0-based (bedGraph) / 1-based (per-column wig) position c%%2e9.  Monotone\n");
    fprintf(stderr, "  (already sorted; no downstream sort/merge) and never exceeds a 32-bit wig coord.  Requires .tui.\n");
    fprintf(stderr, "--bin N : Emit --depth as a ready-to-use binned bedGraph (mean depth per N-bp bin) instead of\n");
    fprintf(stderr, "          per-column.  N must divide 2e9 (so no bin straddles a chunk).  Only bins --depth (not --rs).\n");
    fprintf(stderr, "--sizes FILE : Also write the chrom-sizes file for the chunked uni axis (uni0..uniK from T),\n");
    fprintf(stderr, "          ready for `wigToBigWig --depth.bg --sizes out.bw`.  Requires .tui.\n");
    fprintf(stderr, "-c --useCompression : Bgzip the output(s).\n");
    fprintf(stderr, "-t --tree : Newick tree override.  Default: the `# hal` tree comment in the input header.\n");
    fprintf(stderr, "-r --region SEQ:START-END : Restrict to one anchor chrom range.  Requires <inputFile>.tui.\n");
    fprintf(stderr, "     Half-open, 0-based.\n");
    fprintf(stderr, "     Incompatible with the universal-column output (a leaf region maps to scattered columns).\n");
    fprintf(stderr, "-R --regionFile FILE : Like -r but a file of regions (one SEQ:START-END per line; '#'\n");
    fprintf(stderr, "     comments + blanks ignored).  Amortises one index load -- for SLURM region sharding.\n");
    fprintf(stderr, "   --columnRange LO-HI : Restrict to a universal-column range (half-open, 0-based; HI <= T\n");
    fprintf(stderr, "     from `taffy stats -u`).  Requires .tui.  The natural unit for SLURM sharding.  With\n");
    fprintf(stderr, "     --bin, LO must be a multiple of N (HI too, unless HI == T).  Mutually exclusive with -r/-R.\n");
    fprintf(stderr, "--minLeaves : Minimum surviving non-gap leaves to score a column (default 2).\n");
    fprintf(stderr, "--branchScale : Global multiplier on branch lengths (default 1.0).\n");
    fprintf(stderr, "--paralog MODE : Duplicate-species rows in one block (default union):\n");
    fprintf(stderr, "                 union -- OR each species's paralog bases into a multi-state leaf cset (Hartigan).\n");
    fprintf(stderr, "                 skip  -- drop the entire block (strict GERP++ semantics).\n");
    fprintf(stderr, "                 first -- score using only the first-seen row per species.\n");
    fprintf(stderr, "--skipParalogs : Alias for --paralog skip.\n");
    fprintf(stderr, "--keepParalogs : Alias for --paralog first.\n");
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
    OPT_RS,
    OPT_DEPTH,
    OPT_SIZES,
    OPT_BIN,
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
        fprintf(stderr, "taffy depth: cannot open --regionFile: %s\n", path);
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
            fprintf(stderr, "taffy depth: --regionFile %s line %" PRIi64 ": "
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
            fprintf(stderr, "taffy depth: --regionFile %s line %" PRIi64 ": "
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

int taf_depth_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    char *logLevelString = NULL;
    char *inputFile      = NULL;
    char *rsFile         = NULL;   // --rs    : GERP RS conservation wig (opt-in)
    char *depthFile      = NULL;   // --depth : leaf-depth wig (opt-in)
    char *treeFile       = NULL;
    char *region         = NULL;
    char *regionFile     = NULL;
    char *columnRangeArg = NULL;
    char *sizesFile      = NULL;   // --sizes : chrom-sizes for the chunked uni axis
    bool use_compression = false;
    double branch_scale  = 1.0;
    int64_t min_leaves   = 2;
    GerpParalogPolicy paralog_policy = GERP_PARALOG_UNION;
    int n_threads        = 1;
    int64_t bin_size     = 0;

    while (1) {
        static struct option long_options[] = {
            { "logLevel",       required_argument, 0, 'l' },
            { "inputFile",      required_argument, 0, 'i' },
            { "rs",             required_argument, 0, OPT_RS },
            { "depth",          required_argument, 0, OPT_DEPTH },
            { "useCompression", no_argument,       0, 'c' },
            { "tree",           required_argument, 0, 't' },
            { "region",         required_argument, 0, 'r' },
            { "regionFile",     required_argument, 0, 'R' },
            { "columnRange",    required_argument, 0, OPT_COLUMN_RANGE },
            { "sizes",          required_argument, 0, OPT_SIZES },
            { "branchScale",    required_argument, 0, OPT_BRANCH_SCALE },
            { "minLeaves",      required_argument, 0, OPT_MIN_LEAVES },
            { "paralog",        required_argument, 0, OPT_PARALOG },
            { "skipParalogs",   no_argument,       0, OPT_SKIP_PARALOGS },
            { "keepParalogs",   no_argument,       0, OPT_KEEP_PARALOGS },
            { "threads",        required_argument, 0, 'T' },
            { "bin",            required_argument, 0, OPT_BIN },
            { "help",           no_argument,       0, 'h' },
            { 0, 0, 0, 0 }
        };
        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:ct:r:R:T:h", long_options, &option_index);
        if (key == -1) break;
        switch (key) {
            case 'l': logLevelString = optarg; break;
            case 'i': inputFile      = optarg; break;
            case OPT_RS:    rsFile    = optarg; break;
            case OPT_DEPTH: depthFile = optarg; break;
            case 'c': use_compression = true;  break;
            case 't': treeFile       = optarg; break;
            case 'r': region         = optarg; break;
            case 'R': regionFile     = optarg; break;
            case 'T': n_threads      = atoi(optarg); break;
            case OPT_SIZES:         sizesFile     = optarg; break;
            case OPT_BRANCH_SCALE:  branch_scale  = atof(optarg); break;
            case OPT_MIN_LEAVES:    min_leaves    = atol(optarg); break;
            case OPT_PARALOG:
                if (parse_paralog_policy(optarg, &paralog_policy) != 0) {
                    fprintf(stderr, "taffy depth: --paralog must be one of union|skip|first (got %s)\n",
                            optarg);
                    return 1;
                }
                break;
            case OPT_SKIP_PARALOGS: paralog_policy = GERP_PARALOG_SKIP;  break;
            case OPT_KEEP_PARALOGS: paralog_policy = GERP_PARALOG_FIRST; break;
            case OPT_COLUMN_RANGE:  columnRangeArg = optarg; break;
            case OPT_BIN:           bin_size       = atoll(optarg); break;
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }
    if (n_threads < 1) n_threads = 1;

    // Output toggles: --rs and --depth are each opt-in (>=1 required, checked
    // below).  No --rs -> skip the RS tree walk entirely (the old --depthOnly).
    bool want_rs    = (rsFile != NULL);
    bool depth_only = !want_rs;

    st_setLogLevelFromString(logLevelString);
    LI_set_bgzf_threads(n_threads);

    // Input.
    FILE *in_fh = (inputFile == NULL) ? stdin : fopen(inputFile, "r");
    if (in_fh == NULL) {
        fprintf(stderr, "taffy depth: cannot open input file: %s\n", inputFile);
        return 1;
    }
    LI *li = LI_construct(in_fh);
    int input_format = check_input_format(LI_peek_at_next_line(li));
    if (input_format != 0 && input_format != 1) {
        fprintf(stderr, "taffy depth: input must be MAF or TAF\n");
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
            fprintf(stderr, "taffy depth: cannot open tree file: %s\n", treeFile);
            return 1;
        }
        fseek(tf, 0, SEEK_END);
        long n = ftell(tf);
        if (n < 0) {
            fprintf(stderr, "taffy depth: ftell failed on tree file: %s\n", treeFile);
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
            fprintf(stderr, "taffy depth: no -t tree given and input header has no `# hal` tree.\n"
                            "  Either pass -t <tree.nwk> or supply input from cactus-hal2maf which preserves\n"
                            "  the tree as a `# hal` comment.\n");
            return 1;
        }
        newick = stString_copy(hal->value);
    }
    GerpTree *gt = gerp_tree_construct(newick);
    free(newick);
    if (gt == NULL) {
        fprintf(stderr, "taffy depth: failed to parse Newick tree\n");
        return 1;
    }
    st_logInfo("taffy depth: tree has %" PRIi64 " leaves\n", gerp_tree_n_leaves(gt));

    // CLI mutual-exclusion check goes BEFORE the output file opens so a
    // bad combination doesn't leave a stub output file behind.  -r, -R,
    // and --columnRange are pairwise exclusive (full-fledged validation
    // of the individual values happens later, after the index is loaded).
    if ((region != NULL) + (regionFile != NULL) + (columnRangeArg != NULL) > 1) {
        fprintf(stderr, "taffy depth: -r, -R, and --columnRange are mutually exclusive\n");
        return 1;
    }
    // At least one output must be requested -- each is opt-in.
    if (rsFile == NULL && depthFile == NULL) {
        fprintf(stderr, "taffy depth: need at least one output -- pass --rs FILE and/or --depth FILE\n");
        return 1;
    }
    // Output is the raw universal column [0,T): a contiguous-column concept, so
    // -r/-R (a leaf region -> scattered, non-monotone universal columns) is out.
    // Use --columnRange (a contiguous column slice) or whole-file streaming.
    if (region != NULL || regionFile != NULL) {
        fprintf(stderr, "taffy depth: -r/-R are incompatible with the universal-column output -- a leaf\n"
                        "  region maps to non-contiguous universal columns.  Use --columnRange or whole-file\n"
                        "  streaming.\n");
        return 1;
    }
    if (bin_size < 0) {
        fprintf(stderr, "taffy depth: --bin N must be > 0\n");
        return 1;
    }
    if (bin_size > 0) {
        // --bin bins the leaf-depth output (mean depth per bin).  It does not
        // bin RS, so require --depth and forbid --rs.
        if (depthFile == NULL) {
            fprintf(stderr, "taffy depth: --bin requires --depth FILE (it bins the depth output)\n");
            return 1;
        }
        if (rsFile != NULL) {
            fprintf(stderr, "taffy depth: --bin only bins --depth; drop --rs (or run RS separately)\n");
            return 1;
        }
        // TUI_UNI_CHUNK (2e9) must be a multiple of the bin width so no bin
        // straddles a chunk boundary -- the per-line chunk selection in
        // gerp_flush_bin is correct only then (true for N=1000).
        if (TUI_UNI_CHUNK % bin_size != 0) {
            fprintf(stderr, "taffy depth: --bin N must divide %lld (the uni-axis chunk size), so no bin\n"
                            "  straddles a chunk boundary (got N=%lld)\n",
                    (long long)TUI_UNI_CHUNK, (long long)bin_size);
            return 1;
        }
    }
    // --sizes writes the chunked uni-axis chrom sizes, derived purely from T; it
    // needs the .tui (which the universal-column output requires anyway).

    // Outputs -- each opt-in (validated above: at least one is set).  There is
    // no stdout default; a wig is written only to its named file.
    LW   *out    = NULL;   // --rs
    FILE *out_fh = NULL;
    if (rsFile != NULL) {
        out_fh = fopen(rsFile, "w");
        if (out_fh == NULL) {
            fprintf(stderr, "taffy depth: cannot open --rs file: %s\n", rsFile);
            return 1;
        }
        out = LW_construct(out_fh, use_compression);
    }
    LW   *dout    = NULL;  // --depth
    FILE *dout_fh = NULL;
    bool want_depth = (depthFile != NULL);
    if (want_depth) {
        dout_fh = fopen(depthFile, "w");
        if (dout_fh == NULL) {
            fprintf(stderr, "taffy depth: cannot open --depth file: %s\n", depthFile);
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
    int64_t          *batch_col = st_calloc((size_t)batch_cap, sizeof(int64_t));  // per-block start universal column (--universal)
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
        st_logInfo("taffy depth: %" PRIi64 " regions read from %s\n",
                   n_regions, regionFile);
    } else if (region != NULL) {
        regions = st_calloc(1, sizeof(GerpRegion));
        regions[0].seq = tai_parse_region(region, &regions[0].start,
                                          &regions[0].length);
        if (regions[0].seq == NULL) {
            fprintf(stderr, "taffy depth: invalid region: %s\n", region);
            free(regions);
            return 1;
        }
        n_regions = 1;
    }

    // Parse --columnRange LO-HI directly into a TuiInterval (the iterator
    // takes column intervals natively, no tui_query needed).  Half-open
    // [LO, HI), validated against T below once the .tui is loaded.
    // Direct column range (--columnRange): no source genome involved, so
    // t_start=0, rev=0.  The extract iterator uses only (start, end).
    TuiInterval col_range_iv = { 0, 0, 0, 0 };
    bool have_col_range = false;
    bool empty_col_range = false;
    if (columnRangeArg != NULL) {
        long long lo = 0, hi = 0;
        char extra = 0;
        if (sscanf(columnRangeArg, "%lld-%lld%c", &lo, &hi, &extra) != 2 ||
            lo < 0 || hi < lo) {
            fprintf(stderr, "taffy depth: invalid --columnRange %s "
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
        // The --bin path runs a column-monotone running binner, so a shard must
        // begin/end ON bin boundaries -- else its first/last bin is processed
        // mid-bin and two shards emit overlapping records for the boundary bin
        // (wigToBigWig rejects).  (Combined with TUI_UNI_CHUNK % bin == 0, this
        // also keeps shard edges on chunk-safe positions.)  HI == T is the lone
        // exception (the axis end, no shard above it); it is allowed below once
        // T is loaded, since T is rarely bin-aligned (e.g. 577 T % 1000 == 721).
        if (bin_size > 0 && lo % bin_size != 0) {
            fprintf(stderr, "taffy depth: with --bin %lld, --columnRange LO must be a multiple\n"
                            "  of %lld (a shard must start on a bin boundary)\n",
                    (long long)bin_size, (long long)bin_size);
            return 1;
        }
    }

    // The output is the universal-column axis (chunked uni0..uniK), which needs
    // T and the per-block column coords from the .tui -- so the .tui is ALWAYS
    // required.  Load it once up front (the outer loop reuses it).  The .tai
    // (regular-MAF leaf-coord index) is no longer used by `taffy depth`.
    Tui  *tui    = NULL;
    Tai  *tai    = NULL;   // unused; kept NULL so the cleanup below is a no-op
    FILE *tai_fh = NULL;
    {
        if (inputFile == NULL) {
            fprintf(stderr, "taffy depth: requires -i <file> with a .tui index (cannot index stdin)\n");
            return 1;
        }
        char *tui_fn = tui_path(inputFile);
        if (access(tui_fn, F_OK) != 0) {
            fprintf(stderr, "taffy depth: requires a .tui index (universal MAF -- build with\n"
                            "  `taffy index -u`): %s\n", tui_fn);
            free(tui_fn);
            return 1;
        }
        tui = tui_load(tui_fn);
        if (tui == NULL) {
            fprintf(stderr, "taffy depth: cannot open .tui: %s\n", tui_fn);
            free(tui_fn);
            return 1;
        }
        free(tui_fn);
        st_logInfo("taffy depth: loaded .tui index (universal MAF mode)\n");

        int64_t T = tui_total_columns(tui);
        if (have_col_range) {
            if (col_range_iv.end > T) {
                fprintf(stderr, "taffy depth: --columnRange %" PRIi64 "-%" PRIi64
                                " exceeds T=%" PRIi64 " (total universal columns).\n",
                        col_range_iv.start, col_range_iv.end, T);
                return 1;
            }
            // HI must also land on a bin boundary (else the final shard and the
            // next shard's first bin overlap) -- EXCEPT HI == T, the axis end:
            // no shard exists above T and the binner emits a correct partial
            // final bin.  T is rarely bin-aligned, so the final shard clamps to T.
            if (bin_size > 0 && col_range_iv.end != T && col_range_iv.end % bin_size != 0) {
                fprintf(stderr, "taffy depth: with --bin %lld, --columnRange HI must be a multiple\n"
                                "  of %lld unless HI == T=%lld (the axis end)\n",
                        (long long)bin_size, (long long)bin_size, (long long)T);
                return 1;
            }
            st_logInfo("taffy depth: column-range mode, [%" PRIi64 ", %" PRIi64 ") of T=%" PRIi64 "\n",
                       col_range_iv.start, col_range_iv.end, T);
        }
        st_logInfo("taffy depth: universal-column output, chunked uni0..uni%" PRIi64
                   " (CHUNK=%lld), T=%" PRIi64 "\n",
                   T ? (T - 1) / TUI_UNI_CHUNK : 0, (long long)TUI_UNI_CHUNK, T);

        // --sizes: emit the chunked uni-axis chrom-sizes file (a pure function of
        // T and TUI_UNI_CHUNK; no block scan).  Full chunks = 2e9, last = T mod
        // 2e9; if T is an exact multiple, the last chunk is uni{K-1} (NO trailing
        // zero-length chunk -- wigToBigWig rejects a 0-size chrom).  T == 0 -> no
        // lines.  Ready for `wigToBigWig out.bg out.sizes out.bw`.
        if (sizesFile != NULL) {
            FILE *sf = fopen(sizesFile, "w");
            if (sf == NULL) {
                fprintf(stderr, "taffy depth: cannot open --sizes file: %s\n", sizesFile);
                return 1;
            }
            if (T > 0) {
                int64_t nchunks = (T + TUI_UNI_CHUNK - 1) / TUI_UNI_CHUNK;  // ceil
                for (int64_t k = 0; k < nchunks; k++) {
                    int64_t size_k = (k < nchunks - 1)
                                     ? TUI_UNI_CHUNK
                                     : (T - (nchunks - 1) * TUI_UNI_CHUNK);
                    fprintf(sf, "uni%" PRIi64 "\t%" PRIi64 "\n", k, size_k);
                }
            }
            fclose(sf);
            st_logInfo("taffy depth: wrote uni-axis sizes to %s\n", sizesFile);
        }
    }

    // Outer loop: 1 iteration in streaming or column-range mode;
    // n_regions iterations in region mode (each builds its own iterator
    // over the shared tui/tai, then runs the inner read/score/emit batched
    // loop to exhaustion).  Empty column-range -> 0 iterations (success
    // no-op for SLURM shards that ended up with a zero slice).
    // --bin running binner: depth accumulates across blocks (in column order)
    // into the current universal-column bin; emitted (chunked) when the bin
    // advances.  The bin index is GLOBAL over [0,T); the 2e9-chunk split happens
    // only at emit time in gerp_flush_bin.
    int64_t cur_bin = -1, cur_sum = 0, cur_cnt = 0;

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
                st_logDebug("taffy depth: region %s:%" PRIi64 "-%" PRIi64
                            " resolved via .tui to %" PRIi64 " universal intervals\n",
                            r->seq, r->start, r->start + r->length, n_uiv);
            } else {
                tai_it = tai_iterator(tai, li, rle, r->seq, r->start, r->length);
                st_logDebug("taffy depth: region %s:%" PRIi64 "-%" PRIi64 " via .tai\n",
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
        // Running universal column for --universal: streaming starts at 0,
        // column-range at LO; += column_number per block in read (Phase A,
        // file == column order).  Unused in region mode, which --universal
        // forbids.
        int64_t uni_col = have_col_range ? col_range_iv.start : 0;

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
                batch_col[n_read] = uni_col;       // this block's first universal column
                uni_col += a->column_number;
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
                                want_depth, depth_only,
                                batch_col[i], bin_size);
            }

            // Phase C: serial emit + accounting in batch order.
            for (int i = 0; i < n_read; i++) {
                GerpBlockResult *r = &results[i];
                if (r->status == GERP_BLOCK_UNKNOWN_SPECIES) {
                    fprintf(stderr, "taffy depth: row in block has species not in tree: %s\n"
                                    "  (pass -t with a tree that covers all leaves in the alignment, or\n"
                                    "   drop the offending rows upstream)\n",
                            r->unknown_seq ? r->unknown_seq : "(unknown)");
                    fatal = 1;
                    break;
                }
                if (r->bad_strand) {
                    fprintf(stderr, "taffy depth: row-0 is on the reverse strand in a block "
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
                if (bin_size > 0) {
                    // --bin: columns are monotone, so merge this block's per-bin
                    // partials into the running binner (the boundary bin is shared
                    // with the previous block) and flush (chunked) when the GLOBAL
                    // bin index advances.  Output is already sorted.
                    for (int64_t k = 0; k < r->bin_n; k++) {
                        int64_t bin = r->bin_first + k;
                        if (bin != cur_bin) {
                            gerp_flush_bin(dout, cur_bin, cur_sum, cur_cnt, bin_size);
                            cur_bin = bin; cur_sum = 0; cur_cnt = 0;
                        }
                        cur_sum += r->bin_sum[k];
                        cur_cnt += r->bin_cnt[k];
                    }
                } else {
                    if (out != NULL && r->rs.len > 0)   LW_putn(out,  r->rs.buf,    r->rs.len);
                    if (want_depth && r->depth.len > 0) LW_putn(dout, r->depth.buf, r->depth.len);
                }
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

    // Flush the last --bin bin (the running binner has no successor block to
    // trigger it).
    if (bin_size > 0)
        gerp_flush_bin(dout, cur_bin, cur_sum, cur_cnt, bin_size);

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
    st_logInfo("taffy depth: %" PRIi64 " blocks read, %" PRIi64 " with paralogs (policy=%s; "
               "%" PRIi64 " block-skips), %" PRIi64 " columns scored in %" PRIi64 " seconds\n",
               n_blocks, n_paralog_blocks, policy_name, n_skipped_paralog,
               n_scored_cols, (int64_t)(time(NULL) - startTime));

    for (int i = 0; i < batch_cap; i++) {
        gerpbuf_destroy(&results[i].rs);
        if (want_depth) gerpbuf_destroy(&results[i].depth);
        free(results[i].bin_sum);
        free(results[i].bin_cnt);
    }
    free(results);
    free(batch_aln);
    free(batch_col);
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
    if (out != NULL) {                 // --rs may be absent
        LW_destruct(out, false);
        fclose(out_fh);
    }
    if (dout != NULL) {
        LW_destruct(dout, false);
        if (depthFile != NULL) fclose(dout_fh);
    }
    return fatal;
}

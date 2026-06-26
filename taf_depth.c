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
#include "tui.h"
#include "gerp.h"
#include "sonLib.h"
#include "bigWig.h"   // in-tree 64-bit universal-depth bigWig writer (--bigwig)
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
    int64_t *bin_sum_vec;   // --per-species: [bin_cap * n_leaves] per-bin per-leaf covered-column counts (NULL otherwise)
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

// The wig is keyed on the raw universal column [0,T) on a single 64-bit axis,
// emitted as chrom uni0 (0-based for the --bin bedGraph, 1-based for the
// per-column wig; the --bin bedGraph matches the `taffy lift --bigwig` reader).
// The 64-bit libBigWig fork stores absolute columns directly, so there is no T
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

// ---- In-tree 64-bit universal-depth bigWig writer (single chrom uni0) -------
// Writes the binned depth DIRECTLY as a 64-bit bigWig on the single uni0 axis
// [0,T), replacing the external 32-bit `wigToBigWig` step (which cannot hold a
// >2^31 universal axis and is anyway misparsed by the 64-bit lift reader).
// Intervals arrive from the Phase-C running binner already sorted and
// non-overlapping; they are batched and streamed via the fork's
// bwAddIntervals (first batch) / bwAppendIntervals (rest).
#define UNIBW_BATCH 65536
typedef struct {
    bigWigFile_t *bw;
    const char  **chroms;   // UNIBW_BATCH pointers, all "uni0"
    uint64_t     *starts;
    uint64_t     *ends;
    float        *values;   // UNIBW_BATCH * N floats (N=1 scalar, N=#leaves --per-species)
    uint32_t      N;        // vector width: 1 = scalar bigWig, >1 = per-species vector bigWig
    int           n;        // entries buffered in the current batch
    bool          wrote;    // first batch -> bwAddIntervals, later -> bwAppend
    bool          failed;   // sticky: a flush errored -> stop writing, surface at close
} UniBW;

// N = vector width: 1 writes a scalar 64-bit bigWig (mean depth); N>1 writes the
// per-species vector bigWig (N covered-counts/interval) via the fork's vector API.
static UniBW *unibw_open(const char *path, int64_t T, uint32_t N) {
    if (N < 1) N = 1;
    if (bwInit(1 << 17) != 0) {
        fprintf(stderr, "taffy depth: bwInit failed (--bigwig)\n");
        return NULL;
    }
    bigWigFile_t *bw = bwOpen((char *) path, NULL, "w");
    if (bw == NULL) {
        fprintf(stderr, "taffy depth: cannot open --bigwig file: %s\n", path);
        bwCleanup();
        return NULL;
    }
    int hdrRv = (N > 1) ? bwCreateHdrVec(bw, 10, N) : bwCreateHdr(bw, 10);  // 10 zoom levels
    if (hdrRv) {
        fprintf(stderr, "taffy depth: bwCreateHdr%s failed (--bigwig)\n", N > 1 ? "Vec" : "");
        bwClose(bw); bwCleanup(); return NULL;
    }
    const char *cn[1] = { "uni0" };
    uint64_t    cl[1] = { (uint64_t) T };       // single axis, length T (= --sizes)
    bw->cl = bwCreateChromList(cn, cl, 1);
    if (bw->cl == NULL || bwWriteHdr(bw)) {
        fprintf(stderr, "taffy depth: bwWriteHdr failed (--bigwig)\n");
        bwClose(bw); bwCleanup(); return NULL;
    }
    UniBW *u  = st_calloc(1, sizeof(UniBW));
    u->bw     = bw;
    u->N      = N;
    u->chroms = st_malloc(UNIBW_BATCH * sizeof(char *));
    for (int i = 0; i < UNIBW_BATCH; i++) u->chroms[i] = "uni0";
    u->starts = st_malloc(UNIBW_BATCH * sizeof(uint64_t));
    u->ends   = st_malloc(UNIBW_BATCH * sizeof(uint64_t));
    u->values = st_malloc((size_t) UNIBW_BATCH * N * sizeof(float));
    u->n = 0; u->wrote = false;
    return u;
}

static int unibw_flush(UniBW *u) {
    if (u->failed) return 1;            // sticky: never write to a broken stream again
    if (u->n == 0) return 0;
    int rv;
    if (u->N > 1) {                     // per-species vector bigWig
        rv = u->wrote
            ? bwAppendIntervalsVec(u->bw, u->starts, u->ends, u->values, (uint32_t) u->n)
            : bwAddIntervalsVec(u->bw, u->chroms, u->starts, u->ends, u->values, (uint32_t) u->n);
    } else {                            // scalar bigWig
        rv = u->wrote
            ? bwAppendIntervals(u->bw, u->starts, u->ends, u->values, (uint32_t) u->n)
            : bwAddIntervals(u->bw, u->chroms, u->starts, u->ends, u->values, (uint32_t) u->n);
    }
    u->wrote = true;
    u->n = 0;
    if (rv) {
        u->failed = true;
        fprintf(stderr, "taffy depth: bwAdd/AppendIntervals failed (--bigwig)\n");
    }
    return rv;
}

// vals points to u->N floats (1 for scalar mean, N for the per-species vector).
static int unibw_add(UniBW *u, int64_t start, int64_t end, const float *vals) {
    if (u->failed) return 1;
    if (u->n == UNIBW_BATCH && unibw_flush(u)) return 1;
    u->starts[u->n] = (uint64_t) start;
    u->ends[u->n]   = (uint64_t) end;
    memcpy(u->values + (size_t) u->n * u->N, vals, (size_t) u->N * sizeof(float));
    u->n++;
    return 0;
}

static int unibw_close(UniBW *u) {
    if (u == NULL) return 0;
    int rv = unibw_flush(u);
    bwClose(u->bw);
    bwCleanup();
    free(u->chroms); free(u->starts); free(u->ends); free(u->values);
    free(u);
    return rv;
}

// --bin integer coords: emit one binned record for the mean depth over the bin's
// columns to the bedGraph text (--depth, %.4g) and/or the 64-bit bigWig
// (--bigwig, float mean) -- either output may be NULL.  The bin spans absolute
// universal columns [bin*N, bin*N+cnt) on chrom uni0 at 0-based [start, start+cnt),
// matching the `taffy lift --bigwig` reader.  The column axis is monotone, so this
// runs off the Phase-C running binner (already sorted, no external sort).  Drops
// all-unscored bins (awk parity).
static int gerp_flush_bin(LW *dout, UniBW *ubw, int64_t bin, int64_t sum,
                          int64_t cnt, int64_t bin_size) {
    if (sum <= 0 || cnt <= 0) return 0;
    int64_t start = bin * bin_size;                  // absolute universal column [0,T)
    int64_t end   = start + cnt;
    double  mean  = (double) sum / (double) cnt;
    if (dout != NULL) {
        char line[96];
        int n = snprintf(line, sizeof line, "uni0\t%" PRId64 "\t%" PRId64 "\t%.4g\n",
                         start, end, mean);
        LW_putn(dout, line, (size_t) n);
    }
    if (ubw != NULL) { float mf = (float) mean; return unibw_add(ubw, start, end, &mf); }
    return 0;
}

// --per-species: flush one bin's per-leaf covered-counts (sumvec[N]) as an
// N-vector record on uni0 [bin*bin_size, +cnt).  Drop a bin where no leaf was
// present (vsum==0; same intent as the scalar sum<=0 drop).  fbuf is N-float
// scratch.  dout (--depth alongside --per-species) gets the UNGATED scalar
// total-depth mean (= sum of the per-leaf counts / cnt).  This matches the
// STANDALONE scalar --depth track only at --minLeaves 1 -- standalone --depth is
// --minLeaves-gated, so the two diverge at --minLeaves >= 2.  (The per-species
// bigWig itself is never --minLeaves-gated, by design.)
static int gerp_flush_bin_vec(LW *dout, UniBW *ubw, int64_t bin,
                              const int64_t *sumvec, int64_t cnt,
                              int64_t bin_size, int64_t N, float *fbuf) {
    if (cnt <= 0) return 0;
    int64_t vsum = 0;
    for (int64_t c = 0; c < N; c++) vsum += sumvec[c];
    if (vsum <= 0) return 0;
    int64_t start = bin * bin_size;
    int64_t end   = start + cnt;
    if (dout != NULL) {
        char line[96];
        int n = snprintf(line, sizeof line, "uni0\t%" PRId64 "\t%" PRId64 "\t%.4g\n",
                         start, end, (double) vsum / (double) cnt);
        LW_putn(dout, line, (size_t) n);
    }
    if (ubw != NULL) {
        for (int64_t c = 0; c < N; c++) fbuf[c] = (float) sumvec[c];
        return unibw_add(ubw, start, end, fbuf);
    }
    return 0;
}

static void score_one_block(const GerpTree *gt, GerpThreadState *ts,
                            const Alignment *aln, GerpBlockResult *res,
                            GerpParalogPolicy policy, int64_t min_leaves,
                            double branch_scale, bool want_depth,
                            bool depth_only, int64_t block_start_col,
                            int64_t bin_size, bool per_species) {
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
    // The per-base integer wig keys on the single 64-bit chrom "uni0", so a
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
            // bin index over [0,T); gerp_flush_bin emits it on the single uni0
            // axis at absolute columns.
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
                    if (per_species)
                        res->bin_sum_vec = st_realloc(res->bin_sum_vec,
                                                      (size_t)nc * n_leaves * sizeof(int64_t));
                    res->bin_cap = nc;
                }
                res->bin_sum[idx] = 0;
                res->bin_cnt[idx] = 0;
                if (per_species)
                    memset(res->bin_sum_vec + idx * n_leaves, 0,
                           (size_t)n_leaves * sizeof(int64_t));
                res->bin_n        = idx + 1;
            }
            if (scored) { res->bin_sum[idx] += depth; res->cols_scored++; }
            res->bin_cnt[idx]++;   // every column -> bin-mean denominator
            if (per_species) {
                // per-species coverage: count EVERY column where a leaf is present,
                // regardless of total depth / min_leaves (a species is covered
                // independent of others).  O(n_leaves) per column.
                int64_t *bv = res->bin_sum_vec + idx * n_leaves;
                for (int64_t i = 0; i < n_leaves; i++)
                    if (ts->leaf_csets[i]) bv[i]++;
            }
        } else if (scored) {
            // Integer universal-column wig on the single 64-bit axis: chrom uni0,
            // 1-based absolute position gcol+1.
            int64_t gcol  = block_start_col + col;
            int64_t wig_pos = gcol + 1;  // 1-based absolute column
            if (!hdr_emitted) {
                gerp_emit_wig_header(res, "uni0", !depth_only, want_depth);
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
    fprintf(stderr, "Output is on the UNIVERSAL-COLUMN axis: the raw column [0,T) on a single 64-bit chrom\n");
    fprintf(stderr, "  uni0 at 0-based (bedGraph) / 1-based (per-column wig) absolute position c.  Monotone\n");
    fprintf(stderr, "  (already sorted; no downstream sort/merge).  Requires the 64-bit libBigWig fork + .tui.\n");
    fprintf(stderr, "--bin N : Emit the binned depth (mean depth per N-bp bin) instead of per-column.  Bins\n");
    fprintf(stderr, "          --depth and/or --bigwig (not --rs).\n");
    fprintf(stderr, "--bigwig FILE : Write the binned depth DIRECTLY as a 64-bit bigWig on the single uni0\n");
    fprintf(stderr, "          axis (in-tree, no external wigToBigWig -- which is 32-bit and cannot hold a\n");
    fprintf(stderr, "          >2^31 axis).  This is the file `taffy lift --bigwig` reads.  Requires --bin.\n");
    fprintf(stderr, "--sizes FILE : Also write the chrom-sizes file for the single uni0 axis (uni0 = T).\n");
    fprintf(stderr, "-c --useCompression : Bgzip the output(s).\n");
    fprintf(stderr, "-t --tree : Newick tree override.  Default: the `# hal` tree comment in the input header.\n");
    fprintf(stderr, "   --columnRange LO-HI : Restrict to a universal-column range (half-open, 0-based; HI <= T\n");
    fprintf(stderr, "     from `taffy stats -u`).  Requires .tui.  The natural unit for SLURM sharding.  With\n");
    fprintf(stderr, "     --bin, LO must be a multiple of N (HI too, unless HI == T).\n");
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
    OPT_BIGWIG,
    OPT_PER_SPECIES,
};

static int parse_paralog_policy(const char *s, GerpParalogPolicy *out) {
    if (s == NULL) return -1;
    if (strcmp(s, "union") == 0) { *out = GERP_PARALOG_UNION; return 0; }
    if (strcmp(s, "skip")  == 0) { *out = GERP_PARALOG_SKIP;  return 0; }
    if (strcmp(s, "first") == 0) { *out = GERP_PARALOG_FIRST; return 0; }
    return -1;
}

int taf_depth_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    char *logLevelString = NULL;
    char *inputFile      = NULL;
    char *rsFile         = NULL;   // --rs    : GERP RS conservation wig (opt-in)
    char *depthFile      = NULL;   // --depth : leaf-depth wig (opt-in)
    char *treeFile       = NULL;
    char *columnRangeArg = NULL;
    char *sizesFile      = NULL;   // --sizes : chrom-size for the single uni0 axis
    char *bigwigFile     = NULL;   // --bigwig: in-tree 64-bit depth bigWig (single uni0 axis)
    bool use_compression = false;
    double branch_scale  = 1.0;
    int64_t min_leaves   = 2;
    GerpParalogPolicy paralog_policy = GERP_PARALOG_UNION;
    int n_threads        = 1;
    int64_t bin_size     = 0;
    bool per_species     = false;  // --perSpecies: write a per-leaf vector bigWig (N components)

    while (1) {
        static struct option long_options[] = {
            { "logLevel",       required_argument, 0, 'l' },
            { "inputFile",      required_argument, 0, 'i' },
            { "rs",             required_argument, 0, OPT_RS },
            { "depth",          required_argument, 0, OPT_DEPTH },
            { "useCompression", no_argument,       0, 'c' },
            { "tree",           required_argument, 0, 't' },
            { "columnRange",    required_argument, 0, OPT_COLUMN_RANGE },
            { "sizes",          required_argument, 0, OPT_SIZES },
            { "branchScale",    required_argument, 0, OPT_BRANCH_SCALE },
            { "minLeaves",      required_argument, 0, OPT_MIN_LEAVES },
            { "paralog",        required_argument, 0, OPT_PARALOG },
            { "skipParalogs",   no_argument,       0, OPT_SKIP_PARALOGS },
            { "keepParalogs",   no_argument,       0, OPT_KEEP_PARALOGS },
            { "threads",        required_argument, 0, 'T' },
            { "bin",            required_argument, 0, OPT_BIN },
            { "bigwig",         required_argument, 0, OPT_BIGWIG },
            { "perSpecies",     no_argument,       0, OPT_PER_SPECIES },
            { "help",           no_argument,       0, 'h' },
            { 0, 0, 0, 0 }
        };
        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:ct:T:h", long_options, &option_index);
        if (key == -1) break;
        switch (key) {
            case 'l': logLevelString = optarg; break;
            case 'i': inputFile      = optarg; break;
            case OPT_RS:    rsFile    = optarg; break;
            case OPT_DEPTH: depthFile = optarg; break;
            case 'c': use_compression = true;  break;
            case 't': treeFile       = optarg; break;
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
            case OPT_BIGWIG:        bigwigFile     = optarg; break;
            case OPT_PER_SPECIES:   per_species    = true;   break;
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

    // At least one output must be requested -- each is opt-in.  --bigwig is a
    // binned-depth output too (the in-tree 64-bit writer), so it counts.
    if (rsFile == NULL && depthFile == NULL && bigwigFile == NULL) {
        fprintf(stderr, "taffy depth: need at least one output -- pass --rs FILE, --depth FILE and/or --bigwig FILE\n");
        return 1;
    }
    // --bigwig writes the BINNED depth directly as a 64-bit bigWig, so it requires
    // --bin (the bigWig stores per-bin mean depth on the single uni0 axis).
    if (bigwigFile != NULL && bin_size <= 0) {
        fprintf(stderr, "taffy depth: --bigwig requires --bin N (it writes the binned depth)\n");
        return 1;
    }
    if (per_species && bigwigFile == NULL) {
        fprintf(stderr, "taffy depth: --perSpecies requires --bigwig FILE (it writes the per-species vector bigWig)\n");
        return 1;
    }
    if (bin_size < 0) {
        fprintf(stderr, "taffy depth: --bin N must be > 0\n");
        return 1;
    }
    if (bin_size > 0) {
        // --bin bins the leaf-depth output (mean depth per bin).  It does not
        // bin RS, so require a binned depth sink (--depth and/or --bigwig) and
        // forbid --rs.
        if (depthFile == NULL && bigwigFile == NULL) {
            fprintf(stderr, "taffy depth: --bin requires --depth FILE and/or --bigwig FILE (it bins the depth output)\n");
            return 1;
        }
        if (rsFile != NULL) {
            fprintf(stderr, "taffy depth: --bin only bins --depth/--bigwig; drop --rs (or run RS separately)\n");
            return 1;
        }
    }
    // --sizes writes the single uni0-axis chrom size (= T), derived purely from
    // the .tui (which the universal-column output requires anyway).

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
    UniBW *ubw = NULL;     // --bigwig (opened below, once T is known)

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
        // (wigToBigWig rejects).  HI == T is the lone
        // exception (the axis end, no shard above it); it is allowed below once
        // T is loaded, since T is rarely bin-aligned (e.g. 577 T % 1000 == 721).
        if (bin_size > 0 && lo % bin_size != 0) {
            fprintf(stderr, "taffy depth: with --bin %lld, --columnRange LO must be a multiple\n"
                            "  of %lld (a shard must start on a bin boundary)\n",
                    (long long)bin_size, (long long)bin_size);
            return 1;
        }
    }

    // The output is the universal-column axis (single chrom uni0), which needs
    // T and the per-block column coords from the .tui -- so the .tui is ALWAYS
    // required.  Load it once up front.
    Tui *tui = NULL;
    int64_t T = 0;   // total universal columns = single uni0 axis length; set once below
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

        T = tui_total_columns(tui);
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
        st_logInfo("taffy depth: universal-column output, single axis uni0, T=%" PRIi64 "\n", T);

        // --sizes: emit the single uni0-axis chrom-sizes file (one line: uni0 T),
        // a pure function of T (no block scan).  T == 0 -> no lines.  Ready for
        // `wigToBigWig out.bg out.sizes out.bw`.
        if (sizesFile != NULL) {
            FILE *sf = fopen(sizesFile, "w");
            if (sf == NULL) {
                fprintf(stderr, "taffy depth: cannot open --sizes file: %s\n", sizesFile);
                return 1;
            }
            if (T > 0) fprintf(sf, "uni0\t%" PRIi64 "\n", T);
            fclose(sf);
            st_logInfo("taffy depth: wrote uni-axis sizes to %s\n", sizesFile);
        }
    }

    // --bigwig: open the in-tree 64-bit depth bigWig now that T (the single
    // uni0 axis length) is known.  Requires --bin (validated above), so the
    // Phase-C running binner drives unibw_add via gerp_flush_bin.  NB: under
    // --columnRange (SLURM sharding) each shard writes only its column slice;
    // merging per-shard bigWigs is a separate step (the whole-axis run needs no
    // merge).
    if (bigwigFile != NULL) {
        uint32_t bwN = per_species ? (uint32_t) n_leaves : 1;
        ubw = unibw_open(bigwigFile, T, bwN);
        if (ubw == NULL) {           // unibw_open already reported the reason
            remove(bigwigFile);      // drop any partial file it may have created
            return 1;
        }
        if (per_species) {
            // sidecar <bigwig>.names: one leaf name per line in gerp post-order
            // (= the vector component order), so the lift/browser can map
            // component c -> species without re-deriving the tree.
            char namesPath[4096];
            snprintf(namesPath, sizeof namesPath, "%s.names", bigwigFile);
            FILE *nf = fopen(namesPath, "w");
            if (nf == NULL) {
                fprintf(stderr, "taffy depth: cannot open per-species names sidecar: %s\n", namesPath);
                unibw_close(ubw); remove(bigwigFile); return 1;
            }
            for (int64_t i = 0; i < n_leaves; i++) {
                const char *nm = gerp_tree_leaf_name(gt, i);
                fprintf(nf, "%s\n", nm ? nm : "");
            }
            fclose(nf);
            st_logInfo("taffy depth: writing per-species vector bigWig to %s (%" PRIi64
                       " leaves, single uni0 axis, T=%" PRIi64 ")\n", bigwigFile, n_leaves, T);
        } else {
            st_logInfo("taffy depth: writing 64-bit depth bigWig to %s (single uni0 axis, T=%" PRIi64 ")\n",
                       bigwigFile, T);
        }
    }

    // One pass: column-range mode (a direct column-slice iterator) or whole-file
    // streaming (tui_it == NULL).  Empty column-range -> 0 iterations (a success
    // no-op for a SLURM shard that ended up with a zero slice).
    // --bin running binner: depth accumulates across blocks (in column order)
    // into the current universal-column bin; emitted when the bin advances.  The
    // bin index is GLOBAL over [0,T); gerp_flush_bin emits on the single uni0
    // axis at absolute columns.
    int64_t cur_bin = -1, cur_sum = 0, cur_cnt = 0;
    int64_t *cur_sum_vec = NULL;   // --perSpecies: running per-leaf counts for the open bin
    float   *fbuf = NULL;          // --perSpecies: N-float scratch for unibw_add
    if (per_species) {
        cur_sum_vec = st_calloc((size_t) n_leaves, sizeof(int64_t));
        fbuf        = st_malloc((size_t) n_leaves * sizeof(float));
    }

    int64_t n_iter = empty_col_range ? 0 : 1;
    for (int64_t iter = 0; iter < n_iter && !fatal; iter++) {
        TuiExtractIt *tui_it = NULL;

        if (have_col_range) {
            // Direct column-range iterator: hand TuiInterval{LO,HI} straight to
            // tui_extract_iterator (no tui_query / chrom resolution).  Otherwise
            // tui_it stays NULL and the inner loop streams the whole file.
            tui_it = tui_extract_iterator(tui, li, input_format == 1, rle,
                                          &col_range_iv, 1);
        }

        // TAF read needs the previous block (p_aln) for delta-coord
        // decoding.  We carry it across batches so the first block of
        // batch N+1 has batch N's last block as its predecessor.  MAF
        // reads ignore p_aln.  In column-range mode the iterator manages
        // its own state, so the carry chain is unused.
        Alignment *carry_aln = NULL;
        // Running universal column: streaming starts at 0, column-range at LO;
        // += column_number per block in read (Phase A, file order == column
        // order).
        int64_t uni_col = have_col_range ? col_range_iv.start : 0;

        while (!fatal) {
            // Phase A: serial read of up to batch_cap blocks.  TAF chains
            // through p_aln; once we've passed a block to taf_read_block as
            // p_aln, we can free the prior carry (its successor is now in
            // hand).  Column-range mode forces batch_cap=1: tui_extract_next's
            // returned alignment is invalidated by the NEXT call to it, so we
            // can only hold one at a time.  bgzf_mt still parallelises
            // decompress; cross-shard parallelism covers any lost OMP throughput.
            int eff_batch_cap = (tui_it != NULL) ? 1 : batch_cap;
            int n_read = 0;
            Alignment *p_for_read = carry_aln;
            while (n_read < eff_batch_cap) {
                Alignment *a = NULL;
                if (tui_it != NULL) {
                    a = tui_extract_next(tui_it, li);
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
                // (it was used to decode this new block); safe to free.
                // Column-range mode doesn't use the carry chain.
                if (tui_it == NULL &&
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
                // No more blocks in the stream.
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
                                batch_col[i], bin_size, per_species);
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
                    // with the previous block) and flush when the GLOBAL bin index
                    // advances.  Output is already sorted.
                    for (int64_t k = 0; k < r->bin_n; k++) {
                        int64_t bin = r->bin_first + k;
                        if (bin != cur_bin) {
                            int frv = per_species
                                ? gerp_flush_bin_vec(dout, ubw, cur_bin, cur_sum_vec, cur_cnt, bin_size, n_leaves, fbuf)
                                : gerp_flush_bin(dout, ubw, cur_bin, cur_sum, cur_cnt, bin_size);
                            if (frv) fatal = 1;   // --bigwig write failed (sticky); batch loop exits
                            cur_bin = bin; cur_sum = 0; cur_cnt = 0;
                            if (per_species) memset(cur_sum_vec, 0, (size_t) n_leaves * sizeof(int64_t));
                        }
                        if (per_species)
                            for (int64_t c = 0; c < n_leaves; c++)
                                cur_sum_vec[c] += r->bin_sum_vec[k * n_leaves + c];
                        else
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
            // free everything cleanly.  Column-range mode: tui_extract_iterator
            // owns the yielded alignment (frees it on the next call) so we
            // null-out without destructing.
            if (tui_it != NULL) {
                for (int i = 0; i < n_read; i++) batch_aln[i] = NULL;
            } else {
                int free_to = fatal ? n_read : n_read - 1;
                for (int i = 0; i < free_to; i++) {
                    if (batch_aln[i] != NULL) {
                        alignment_destruct(batch_aln[i], 1);
                        batch_aln[i] = NULL;
                    }
                }
                if (!fatal && n_read > 0) {
                    carry_aln = batch_aln[n_read - 1];
                    batch_aln[n_read - 1] = NULL;
                }
            }
        }
        if (carry_aln != NULL) alignment_destruct(carry_aln, 1);

        if (tui_it != NULL) tui_extract_iterator_destruct(tui_it);
    }

    // Flush the last --bin bin (the running binner has no successor block to
    // trigger it).
    if (bin_size > 0) {
        int frv = per_species
            ? gerp_flush_bin_vec(dout, ubw, cur_bin, cur_sum_vec, cur_cnt, bin_size, n_leaves, fbuf)
            : gerp_flush_bin(dout, ubw, cur_bin, cur_sum, cur_cnt, bin_size);
        if (frv) fatal = 1;
    }

    if (tui != NULL) tui_destruct(tui);

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
        free(results[i].bin_sum_vec);
    }
    free(results);
    free(batch_aln);
    free(batch_col);
    free(cur_sum_vec);
    free(fbuf);
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
    if (ubw != NULL) {                 // --bigwig: flush final batch + finalize index
        if (unibw_close(ubw) != 0) fatal = 1;
    }
    if (fatal && bigwigFile != NULL)   // never leave a truncated/corrupt .bw that looks complete
        remove(bigwigFile);
    return fatal;
}

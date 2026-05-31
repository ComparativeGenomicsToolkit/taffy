/*
 * taffy gerp-stats: post-process a gerp RS wig + depth wig into a
 * depth-corrected + percentile-ranked wig.
 *
 * Three streaming passes (each is one walk of the joint wig pair):
 *   Pass 1: per-depth (count, sum, sum_sq) -> mean(d) + stddev(d).
 *   Pass 2: z-score per column -> histogram of z values.
 *   Pass 3: emit (chrom, pos, percentile) wig (and/or zscore wig).
 *
 * Each pass batches N tuples and dispatches the per-line work to an
 * OpenMP worker pool.  bgzf_mt parallelises the decompress on input.
 *
 * Input wig contract: variableStep, one chrom per variableStep header,
 * data lines of "<pos> <value>".  Depth wig is a positional superset of
 * the RS wig (every non-gap row-0 base in the original universal MAF);
 * RS wig has only the scored positions.  Tool asserts the RS positions
 * are all present in the depth wig in the same chrom order.
 */

#include "taf.h"
#include "gerp.h"
#include "gerp_stats.h"
#include "sonLib.h"
#include <ctype.h>
#include <getopt.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////////////////////////////////////////////////////////
// Streaming wig reader.  Reads one wig line at a time.  variableStep
// headers update the current chrom; data lines yield (pos, value).
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    LI    *li;
    char  *cur_chrom;        // owned strdup; replaced on each variableStep
    bool   eof;
    // Pending peeked tuple (set by `_advance` when it's already buffered).
    bool    have_pending;
    int64_t pending_pos;
    double  pending_val;
    char   *pending_chrom;   // borrowed: points into cur_chrom (do not free)
} WigStream;

static WigStream *wig_stream_open(const char *path) {
    FILE *fh = fopen(path, "r");
    if (fh == NULL) {
        fprintf(stderr, "taffy gerp-stats: cannot open wig file: %s\n", path);
        return NULL;
    }
    WigStream *ws = st_calloc(1, sizeof(WigStream));
    ws->li = LI_construct(fh);
    return ws;
}

static void wig_stream_close(WigStream *ws) {
    if (ws == NULL) return;
    if (ws->li != NULL) LI_destruct(ws->li);
    free(ws->cur_chrom);
    free(ws);
}

// Parse one wig line into ws's pending slot, or set ws->eof if we ran out.
// Returns true if a data tuple is now available in `pending_*`; false on
// EOF.  Skips comment / blank / fixedStep lines and updates cur_chrom on
// variableStep.
static bool wig_stream_next(WigStream *ws) {
    if (ws->eof) return false;
    while (1) {
        char *line = LI_get_next_line(ws->li);
        if (line == NULL) { ws->eof = true; ws->have_pending = false; return false; }
        // Skip leading whitespace.
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0' || *p == '\n' || *p == '#') { free(line); continue; }
        // variableStep header.
        if (strncmp(p, "variableStep", 12) == 0) {
            // Look for "chrom=<name>".
            char *c = strstr(p, "chrom=");
            if (c == NULL) {
                fprintf(stderr, "taffy gerp-stats: variableStep line missing chrom: %s\n", p);
                free(line);
                ws->eof = true;
                return false;
            }
            c += 6;
            char *e = c;
            while (*e && *e != ' ' && *e != '\t' && *e != '\n' && *e != '\r') e++;
            *e = '\0';
            free(ws->cur_chrom);
            ws->cur_chrom = stString_copy(c);
            free(line);
            continue;
        }
        if (strncmp(p, "fixedStep", 9) == 0) {
            fprintf(stderr, "taffy gerp-stats: fixedStep wig is not supported -- "
                            "convert with wigToBedGraph or similar first.\n");
            free(line);
            ws->eof = true;
            return false;
        }
        if (strncmp(p, "track", 5) == 0 || strncmp(p, "browser", 7) == 0) {
            // UCSC header line.  Tolerant skip.
            free(line);
            continue;
        }
        // Data line: "<pos> <value>".
        char *endp = NULL;
        long long pos = strtoll(p, &endp, 10);
        if (endp == p) {
            // Garbage line; skip but warn.
            st_logDebug("taffy gerp-stats: skipping unparseable wig line: %.40s\n", p);
            free(line);
            continue;
        }
        while (*endp == ' ' || *endp == '\t') endp++;
        char *valp = endp;
        double val = strtod(valp, &endp);
        if (endp == valp) {
            st_logDebug("taffy gerp-stats: missing value on wig line: %.40s\n", p);
            free(line);
            continue;
        }
        if (ws->cur_chrom == NULL) {
            fprintf(stderr, "taffy gerp-stats: data line before any variableStep header\n");
            free(line);
            ws->eof = true;
            return false;
        }
        ws->have_pending  = true;
        ws->pending_pos   = (int64_t)pos;
        ws->pending_val   = val;
        ws->pending_chrom = ws->cur_chrom;
        free(line);
        return true;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Joint streamer: drives the RS wig; advances the depth wig in lockstep
// to find the same (chrom, pos).  Returns false on EOF of RS wig.
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    WigStream *rs;
    WigStream *depth;       // may be NULL if no depth file was supplied
    // Track when chrom changes so the caller can emit variableStep
    // headers in pass 3.  Set on every successful tuple read; toggled
    // to false after the caller acknowledges via _ack_chrom().
    //
    // last_emitted_chrom is owned (stString_copy on transition + free on
    // close): WigStream::cur_chrom is freed/recopied on every variableStep
    // line, so a borrowed pointer here would dangle across the next
    // variableStep.
    bool   chrom_changed;
    char  *last_emitted_chrom;
} JointWig;

static JointWig *joint_open(const char *rs_path, const char *depth_path) {
    WigStream *rs = wig_stream_open(rs_path);
    if (rs == NULL) return NULL;
    WigStream *de = NULL;
    if (depth_path != NULL) {
        de = wig_stream_open(depth_path);
        if (de == NULL) { wig_stream_close(rs); return NULL; }
    }
    JointWig *j = st_calloc(1, sizeof(JointWig));
    j->rs    = rs;
    j->depth = de;
    return j;
}

static void joint_close(JointWig *j) {
    if (j == NULL) return;
    wig_stream_close(j->rs);
    wig_stream_close(j->depth);
    free(j->last_emitted_chrom);
    free(j);
}

// Advance the depth stream until its current tuple matches (chrom, pos)
// from the RS stream.  Returns true on match.  Returns false (and prints
// an error) if we run off the end of depth OR pass `pos` without hitting
// it (depth was expected to be a superset).
static bool advance_depth_to(JointWig *j, const char *chrom, int64_t pos,
                             int64_t *out_depth) {
    if (j->depth == NULL) { *out_depth = 0; return true; }
    // Cap the number of cross-chrom skips so a corrupted depth wig
    // with millions of extra chroms doesn't spin silently here.  The
    // expected ratio is 1:1 same-chrom and ~K extra cross-chrom skips
    // at each chrom boundary (where K = chrom count of depth wig); a
    // 1 M cap catches truly garbage input without false-positives on
    // legitimate wigs.
    int64_t skipped = 0;
    const int64_t kSkipCap = 1000000;
    while (1) {
        if (!j->depth->have_pending && !wig_stream_next(j->depth)) {
            fprintf(stderr, "taffy gerp-stats: depth wig EOF before RS position "
                            "%s:%" PRIi64 "\n", chrom, pos);
            return false;
        }
        if (strcmp(j->depth->pending_chrom, chrom) != 0) {
            // Different chrom: discard this depth tuple, keep advancing.
            j->depth->have_pending = false;
            if (++skipped >= kSkipCap) {
                fprintf(stderr, "taffy gerp-stats: depth wig: %" PRIi64 " consecutive "
                                "cross-chrom records without finding %s:%" PRIi64 " -- "
                                "input likely corrupt or mismatched.\n",
                        skipped, chrom, pos);
                return false;
            }
            continue;
        }
        if (j->depth->pending_pos < pos) {
            j->depth->have_pending = false;
            continue;
        }
        if (j->depth->pending_pos == pos) {
            *out_depth = (int64_t)j->depth->pending_val;
            j->depth->have_pending = false;
            return true;
        }
        // depth.pos > pos -> the RS position isn't in the depth wig.
        fprintf(stderr, "taffy gerp-stats: depth wig has no entry for RS position "
                        "%s:%" PRIi64 " (next depth: %s:%" PRIi64 "). "
                        "The two wigs must come from the same gerp run.\n",
                chrom, pos,
                j->depth->pending_chrom, j->depth->pending_pos);
        return false;
    }
}

// Pull the next (chrom, pos, rs, depth) tuple.  Returns false on EOF of RS.
// On error inside depth alignment, returns false AND sets *err = 1.
static bool joint_next(JointWig *j, const char **chrom_out, int64_t *pos_out,
                       double *rs_out, int64_t *depth_out, int *err) {
    *err = 0;
    if (!j->rs->have_pending && !wig_stream_next(j->rs)) return false;
    const char *chrom = j->rs->pending_chrom;
    int64_t     pos   = j->rs->pending_pos;
    double      rs    = j->rs->pending_val;
    j->rs->have_pending = false;

    int64_t depth = 0;
    if (!advance_depth_to(j, chrom, pos, &depth)) { *err = 1; return false; }

    *chrom_out = chrom;
    *pos_out   = pos;
    *rs_out    = rs;
    *depth_out = depth;
    if (j->last_emitted_chrom == NULL ||
        strcmp(j->last_emitted_chrom, chrom) != 0) {
        j->chrom_changed = true;
        // Own the chrom copy: WigStream's cur_chrom is freed at the next
        // variableStep, so a borrowed pointer would dangle.
        free(j->last_emitted_chrom);
        j->last_emitted_chrom = stString_copy(chrom);
    } else {
        j->chrom_changed = false;
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////
// Pass 1: per-depth stats.  Stream the joint wig; tally (count, sum,
// sum_sq) per integer depth.  No depth-correction yet.
/////////////////////////////////////////////////////////////////////////////

static int pass1_collect_stats(JointWig *j, DepthStats *ds_global,
                               int n_threads, int batch_size, int64_t *n_total) {
    // Batch of (depth, rs) tuples.  Workers populate per-thread DepthStats.
    int64_t *b_depth = st_malloc((size_t)batch_size * sizeof(int64_t));
    double  *b_rs    = st_malloc((size_t)batch_size * sizeof(double));
    DepthStats **per_thread = st_malloc((size_t)n_threads * sizeof(DepthStats *));
    for (int t = 0; t < n_threads; t++) per_thread[t] = depth_stats_construct(ds_global->max_depth);

    int err = 0;
    int64_t total = 0;
    while (!err) {
        int n = 0;
        while (n < batch_size) {
            const char *chrom = NULL;
            int64_t pos = 0;
            double  rs  = 0.0;
            int64_t d   = 0;
            int     e   = 0;
            if (!joint_next(j, &chrom, &pos, &rs, &d, &e)) { err = e; break; }
            b_depth[n] = d;
            b_rs[n]    = rs;
            n++;
        }
        if (n == 0) break;
        total += n;
        #pragma omp parallel for schedule(static) num_threads(n_threads)
        for (int i = 0; i < n; i++) {
#ifdef _OPENMP
            int t = omp_get_thread_num();
#else
            int t = 0;
#endif
            depth_stats_observe(per_thread[t], b_depth[i], b_rs[i]);
        }
    }
    for (int t = 0; t < n_threads; t++) {
        depth_stats_merge(ds_global, per_thread[t]);
        depth_stats_destruct(per_thread[t]);
    }
    free(per_thread);
    free(b_depth);
    free(b_rs);
    *n_total = total;
    return err;
}

/////////////////////////////////////////////////////////////////////////////
// Pass 2: histogram of (depth-corrected) z-scores.
/////////////////////////////////////////////////////////////////////////////

static int pass2_collect_histogram(JointWig *j, const DepthStats *ds,
                                   Histogram *h_global, int n_threads,
                                   int batch_size) {
    int64_t *b_depth = st_malloc((size_t)batch_size * sizeof(int64_t));
    double  *b_rs    = st_malloc((size_t)batch_size * sizeof(double));
    Histogram **per_thread = st_malloc((size_t)n_threads * sizeof(Histogram *));
    for (int t = 0; t < n_threads; t++)
        per_thread[t] = histogram_construct(h_global->n_bins, h_global->min_val, h_global->max_val);

    int err = 0;
    while (!err) {
        int n = 0;
        while (n < batch_size) {
            const char *chrom = NULL;
            int64_t pos = 0;
            double  rs  = 0.0;
            int64_t d   = 0;
            int     e   = 0;
            if (!joint_next(j, &chrom, &pos, &rs, &d, &e)) { err = e; break; }
            b_depth[n] = d;
            b_rs[n]    = rs;
            n++;
        }
        if (n == 0) break;
        #pragma omp parallel for schedule(static) num_threads(n_threads)
        for (int i = 0; i < n; i++) {
#ifdef _OPENMP
            int t = omp_get_thread_num();
#else
            int t = 0;
#endif
            double z = depth_stats_zscore(ds, b_depth[i], b_rs[i]);
            histogram_observe(per_thread[t], z);
        }
    }
    for (int t = 0; t < n_threads; t++) {
        histogram_merge(h_global, per_thread[t]);
        histogram_destruct(per_thread[t]);
    }
    free(per_thread);
    free(b_depth);
    free(b_rs);
    return err;
}

/////////////////////////////////////////////////////////////////////////////
// Pass 3: emit output wig with the chosen (z, percentile, or both).
/////////////////////////////////////////////////////////////////////////////

typedef enum { EMIT_PERCENTILE, EMIT_ZSCORE, EMIT_BOTH } EmitMode;

// Reuse the per-block GerpBuf-style buffers from taf_gerp.c via a copy.
typedef struct {
    char  *buf;
    size_t len;
    size_t cap;
} OutBuf;

static void outbuf_init(OutBuf *b, size_t cap)  { b->cap = cap; b->buf = st_malloc(cap); b->len = 0; }
static void outbuf_destroy(OutBuf *b)           { free(b->buf); b->buf = NULL; b->len = b->cap = 0; }
static void outbuf_reset(OutBuf *b)             { b->len = 0; }
static inline void outbuf_reserve(OutBuf *b, size_t n) {
    if (b->len + n <= b->cap) return;
    while (b->len + n > b->cap) b->cap *= 2;
    b->buf = st_realloc(b->buf, b->cap);
}
static inline void outbuf_putn(OutBuf *b, const char *s, size_t n) {
    outbuf_reserve(b, n);
    memcpy(b->buf + b->len, s, n);
    b->len += n;
}

/////////////////////////////////////////////////////////////////////////////
// Clade attribution.  User supplies a set of "lineage root" labels (e.g.
// MammalsAnc0, BirdsAnc0, RayFinnedFishesAnc0); each column's anchor
// ancestor (row-0 of its block, encoded in the wig chrom name) is mapped
// to the lowest lineage root that contains it -- walking up the tree
// from the anchor until we hit a member of the roots set.  Columns whose
// anchor is above every root land in CLADE_ABOVE_ROOTS; unresolvable
// chroms (genome label not in the tree at all) land in CLADE_UNKNOWN.
//
// Cache: chrom_name -> clade_id (small int).  Populated single-threaded
// in pass 3's serial read phase, read-only in the parallel phase.
/////////////////////////////////////////////////////////////////////////////

#define CLADE_UNKNOWN     ((int64_t) -1)
#define CLADE_ABOVE_ROOTS ((int64_t) -2)

typedef struct {
    const char *name;          // borrowed from `roots_set` or sentinel string
    int64_t     depth_from_root;
} CladeBucket;

typedef struct {
    GerpTree    *gt;
    stSet       *roots_set;    // owned; lineage root labels
    CladeBucket *buckets;      // [n_buckets]; user roots first, then above + unknown
    int64_t      n_buckets;
    int64_t      idx_above;    // == n_user_roots
    int64_t      idx_unknown;  // == n_user_roots + 1
    stHash      *chrom_cache;  // chrom_str -> int64_t* clade_id (allocated)
    stHash      *label_to_idx; // root label -> int64_t* (which user-bucket)
} CladeMap;

static int cmp_bucket_by_depth(const void *a, const void *b) {
    const CladeBucket *ba = a, *bb = b;
    if (ba->depth_from_root != bb->depth_from_root) {
        return (ba->depth_from_root < bb->depth_from_root) ? -1 : 1;
    }
    return strcmp(ba->name, bb->name);
}

// Build the CladeMap.  `roots_set` is taken over (we own it from here on).
static CladeMap *clademap_construct(GerpTree *gt, stSet *roots_set) {
    CladeMap *cm = st_calloc(1, sizeof(CladeMap));
    cm->gt = gt;
    cm->roots_set = roots_set;
    int64_t n_user = (int64_t) stSet_size(roots_set);
    cm->n_buckets = n_user + 2;
    cm->idx_above   = n_user;
    cm->idx_unknown = n_user + 1;
    cm->buckets = st_calloc((size_t) cm->n_buckets, sizeof(CladeBucket));
    // Pull root labels into the user-bucket array; sort by depth_from_root.
    stList *root_labels = stSet_getList(roots_set);
    for (int64_t i = 0; i < n_user; i++) {
        const char *lbl = stList_get(root_labels, i);
        cm->buckets[i].name = lbl;  // borrowed from roots_set
        cm->buckets[i].depth_from_root = gerp_tree_depth_from_root(gt, lbl);
    }
    stList_destruct(root_labels);
    qsort(cm->buckets, (size_t) n_user, sizeof(CladeBucket), cmp_bucket_by_depth);
    cm->buckets[cm->idx_above]  .name = "__above_roots__";
    cm->buckets[cm->idx_above]  .depth_from_root = 0;
    cm->buckets[cm->idx_unknown].name = "__unknown__";
    cm->buckets[cm->idx_unknown].depth_from_root = -1;
    // Index: lineage-root label -> bucket id (after the depth sort).
    cm->label_to_idx = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, free);
    for (int64_t i = 0; i < n_user; i++) {
        int64_t *v = st_malloc(sizeof(int64_t));
        *v = i;
        stHash_insert(cm->label_to_idx, (char *) cm->buckets[i].name, v);
    }
    cm->chrom_cache = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    return cm;
}

static void clademap_destruct(CladeMap *cm) {
    if (cm == NULL) return;
    stHash_destruct(cm->chrom_cache);
    stHash_destruct(cm->label_to_idx);
    free(cm->buckets);
    stSet_destruct(cm->roots_set);
    free(cm);
}

// Resolve a chrom string to its clade bucket id.  Mutates the cache on
// first sighting; call only from single-threaded code.  Returns the
// bucket id (0..n_user_roots-1 for user lineages, idx_above for above,
// idx_unknown for unresolvable).
static int64_t clademap_lookup(CladeMap *cm, const char *chrom) {
    int64_t *cached = stHash_search(cm->chrom_cache, (void *) chrom);
    if (cached != NULL) return *cached;
    int64_t bucket;
    const char *genome = gerp_tree_resolve_genome(cm->gt, chrom);
    if (genome == NULL) {
        bucket = cm->idx_unknown;
    } else {
        const char *root = gerp_tree_walk_to_set(cm->gt, genome, cm->roots_set);
        if (root == NULL) {
            bucket = cm->idx_above;
        } else {
            int64_t *p = stHash_search(cm->label_to_idx, (void *) root);
            bucket = (p != NULL) ? *p : cm->idx_unknown;
        }
    }
    int64_t *v = st_malloc(sizeof(int64_t));
    *v = bucket;
    stHash_insert(cm->chrom_cache, stString_copy(chrom), v);
    return bucket;
}

/////////////////////////////////////////////////////////////////////////////
// Per-clade accumulators.  Streamed during pass 3.  Each thread keeps a
// per-clade array; merged into the global after each batch.
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    int64_t n_cols;
    int64_t sum_depth;
    double  sum_rs, sum_sq_rs;
    double  sum_z,  sum_sq_z;
    int64_t n_top10_global;
    int64_t n_bot10_global;
    Histogram *rs_hist;
    Histogram *z_hist;
} CladeAccum;

static CladeAccum *clade_accum_new_array(int64_t n_clades, int64_t hist_bins,
                                          double rs_min, double rs_max,
                                          double z_min,  double z_max) {
    CladeAccum *arr = st_calloc((size_t) n_clades, sizeof(CladeAccum));
    for (int64_t i = 0; i < n_clades; i++) {
        arr[i].rs_hist = histogram_construct(hist_bins, rs_min, rs_max);
        arr[i].z_hist  = histogram_construct(hist_bins, z_min,  z_max);
    }
    return arr;
}

static void clade_accum_destroy_array(CladeAccum *arr, int64_t n_clades) {
    if (arr == NULL) return;
    for (int64_t i = 0; i < n_clades; i++) {
        histogram_destruct(arr[i].rs_hist);
        histogram_destruct(arr[i].z_hist);
    }
    free(arr);
}

static void clade_accum_merge_array(CladeAccum *dst, const CladeAccum *src,
                                    int64_t n_clades) {
    for (int64_t i = 0; i < n_clades; i++) {
        dst[i].n_cols          += src[i].n_cols;
        dst[i].sum_depth       += src[i].sum_depth;
        dst[i].sum_rs          += src[i].sum_rs;
        dst[i].sum_sq_rs       += src[i].sum_sq_rs;
        dst[i].sum_z           += src[i].sum_z;
        dst[i].sum_sq_z        += src[i].sum_sq_z;
        dst[i].n_top10_global  += src[i].n_top10_global;
        dst[i].n_bot10_global  += src[i].n_bot10_global;
        histogram_merge(dst[i].rs_hist, src[i].rs_hist);
        histogram_merge(dst[i].z_hist,  src[i].z_hist);
    }
}

// Pass 3: per-clade stat collection + optional ranked-wig emit.  Caller
// passes a CladeMap (single-threaded resolver, cache built up here as
// chroms are seen) and a global CladeAccum array we accumulate INTO.
// `out` may be NULL -- then no wig is written and we only accumulate.
// `z_top10` / `z_bot10` are the global thresholds derived from the
// histogram's 0.9 / 0.1 quantiles; used to count per-clade columns in
// each global tail.
static int pass3_run(JointWig *j, const DepthStats *ds, const Histogram *h,
                     CladeMap *cm, CladeAccum *clade_global,
                     double z_top10, double z_bot10,
                     LW *out, EmitMode emit_mode,
                     int n_threads, int batch_size) {
    int64_t *b_depth = st_malloc((size_t)batch_size * sizeof(int64_t));
    double  *b_rs    = st_malloc((size_t)batch_size * sizeof(double));
    int64_t *b_pos   = st_malloc((size_t)batch_size * sizeof(int64_t));
    int64_t *b_clade = st_malloc((size_t)batch_size * sizeof(int64_t));
    char   **b_chrom = st_malloc((size_t)batch_size * sizeof(char *));
    bool    *b_first_of_chrom = st_calloc((size_t)batch_size, sizeof(bool));
    OutBuf  *b_out   = NULL;
    if (out != NULL) {
        b_out = st_calloc((size_t)batch_size, sizeof(OutBuf));
        for (int i = 0; i < batch_size; i++) outbuf_init(&b_out[i], 256);
    }
    // Per-thread per-clade accumulators.  Allocated once; merged into
    // clade_global at the end of each batch + reset.
    CladeAccum **per_thread = st_malloc((size_t)n_threads * sizeof(CladeAccum *));
    for (int t = 0; t < n_threads; t++) {
        per_thread[t] = clade_accum_new_array(cm->n_buckets,
                                              clade_global[0].rs_hist->n_bins,
                                              clade_global[0].rs_hist->min_val,
                                              clade_global[0].rs_hist->max_val,
                                              clade_global[0].z_hist->min_val,
                                              clade_global[0].z_hist->max_val);
    }

    int err = 0;
    while (!err) {
        int n = 0;
        while (n < batch_size) {
            const char *chrom = NULL;
            int64_t pos = 0;
            double  rs  = 0.0;
            int64_t d   = 0;
            int     e   = 0;
            if (!joint_next(j, &chrom, &pos, &rs, &d, &e)) { err = e; break; }
            // Own the chrom copy: `chrom` borrows into WigStream::cur_chrom,
            // which the wig reader frees at the next variableStep line --
            // and "next variableStep" can happen later in this same batch
            // on a universal MAF with many anchor chroms per 4096 cols.
            // Freed below after the parallel block emits its variableStep
            // headers (lines 619-621).
            b_chrom[n] = stString_copy(chrom);
            b_pos[n]   = pos;
            b_rs[n]    = rs;
            b_depth[n] = d;
            b_first_of_chrom[n] = j->chrom_changed;
            // SERIAL clade lookup -- mutates the cm->chrom_cache on miss.
            b_clade[n] = clademap_lookup(cm, chrom);
            n++;
        }
        if (n == 0) break;

        #pragma omp parallel for schedule(static) num_threads(n_threads)
        for (int i = 0; i < n; i++) {
#ifdef _OPENMP
            int t = omp_get_thread_num();
#else
            int t = 0;
#endif
            double z = depth_stats_zscore(ds, b_depth[i], b_rs[i]);
            double p = (h != NULL) ? histogram_percentile(h, z) : 0.0;

            // Per-clade accumulate.
            int64_t cid = b_clade[i];
            CladeAccum *ca = &per_thread[t][cid];
            ca->n_cols++;
            ca->sum_depth   += b_depth[i];
            ca->sum_rs      += b_rs[i];
            ca->sum_sq_rs   += b_rs[i] * b_rs[i];
            ca->sum_z       += z;
            ca->sum_sq_z    += z * z;
            histogram_observe(ca->rs_hist, b_rs[i]);
            histogram_observe(ca->z_hist,  z);
            if (z >= z_top10) ca->n_top10_global++;
            if (z <= z_bot10) ca->n_bot10_global++;

            // Optional wig emit.
            if (b_out != NULL) {
                outbuf_reset(&b_out[i]);
                if (b_first_of_chrom[i]) {
                    outbuf_putn(&b_out[i], "variableStep chrom=", 19);
                    outbuf_putn(&b_out[i], b_chrom[i], strlen(b_chrom[i]));
                    outbuf_putn(&b_out[i], "\n", 1);
                }
                char line[96];
                int ln = 0;
                switch (emit_mode) {
                    case EMIT_PERCENTILE:
                        ln = snprintf(line, sizeof(line), "%" PRIi64 " %.2f\n", b_pos[i], p);
                        break;
                    case EMIT_ZSCORE:
                        ln = snprintf(line, sizeof(line), "%" PRIi64 " %.4f\n", b_pos[i], z);
                        break;
                    case EMIT_BOTH:
                        ln = snprintf(line, sizeof(line), "%" PRIi64 "\t%.4f\t%.2f\n", b_pos[i], z, p);
                        break;
                }
                if (ln > 0) outbuf_putn(&b_out[i], line, (size_t)ln);
            }
        }

        if (out != NULL) {
            for (int i = 0; i < n; i++) {
                if (b_out[i].len > 0) LW_putn(out, b_out[i].buf, b_out[i].len);
            }
        }
        // Drop the owned chrom copies now that headers have been emitted.
        for (int i = 0; i < n; i++) {
            free(b_chrom[i]);
            b_chrom[i] = NULL;
        }
        // Merge per-thread per-clade tallies and reset the per-thread copies.
        for (int t = 0; t < n_threads; t++) {
            clade_accum_merge_array(clade_global, per_thread[t], cm->n_buckets);
            // Zero scalar fields; histograms get cleared per-bin.
            for (int64_t k = 0; k < cm->n_buckets; k++) {
                per_thread[t][k].n_cols = 0;
                per_thread[t][k].sum_depth = 0;
                per_thread[t][k].sum_rs = per_thread[t][k].sum_sq_rs = 0;
                per_thread[t][k].sum_z  = per_thread[t][k].sum_sq_z  = 0;
                per_thread[t][k].n_top10_global = 0;
                per_thread[t][k].n_bot10_global = 0;
                memset(per_thread[t][k].rs_hist->bin_count, 0,
                       (size_t) per_thread[t][k].rs_hist->n_bins * sizeof(int64_t));
                per_thread[t][k].rs_hist->n_total = 0;
                per_thread[t][k].rs_hist->n_clipped_lo = 0;
                per_thread[t][k].rs_hist->n_clipped_hi = 0;
                memset(per_thread[t][k].z_hist->bin_count, 0,
                       (size_t) per_thread[t][k].z_hist->n_bins * sizeof(int64_t));
                per_thread[t][k].z_hist->n_total = 0;
                per_thread[t][k].z_hist->n_clipped_lo = 0;
                per_thread[t][k].z_hist->n_clipped_hi = 0;
            }
        }
    }
    for (int t = 0; t < n_threads; t++) {
        clade_accum_destroy_array(per_thread[t], cm->n_buckets);
    }
    free(per_thread);
    if (b_out != NULL) {
        for (int i = 0; i < batch_size; i++) outbuf_destroy(&b_out[i]);
        free(b_out);
    }
    free(b_first_of_chrom);
    free(b_chrom);
    free(b_clade);
    free(b_pos);
    free(b_rs);
    free(b_depth);
    return err;
}

// Emit the per-clade TSV summary.  Columns:
//   clade  depth_from_root  n_cols  bp_relative_pct
//   mean_rs  median_rs  p5_rs  p95_rs
//   mean_z   median_z   p5_z   p95_z
//   mean_depth  frac_top10_global  frac_bot10_global
// Followed by a "__total__" row that's the sum of all buckets (for
// sanity checking the input column count).
static void emit_tsv_summary(FILE *fh, const CladeMap *cm, CladeAccum *acc,
                             int64_t n_total_columns) {
    // Finalise the per-clade histograms so quantile lookups work.
    for (int64_t i = 0; i < cm->n_buckets; i++) {
        histogram_finalize(acc[i].rs_hist);
        histogram_finalize(acc[i].z_hist);
    }
    fprintf(fh,
        "clade\tdepth_from_root\tn_cols\tbp_relative_pct\t"
        "mean_rs\tmedian_rs\tp5_rs\tp95_rs\t"
        "mean_z\tmedian_z\tp5_z\tp95_z\t"
        "mean_depth\tfrac_top10pct_global\tfrac_bot10pct_global\n");
    for (int64_t i = 0; i < cm->n_buckets; i++) {
        CladeAccum *a = &acc[i];
        if (a->n_cols == 0 && i < cm->idx_above) continue;  // empty user bucket: skip
        double n = (double) a->n_cols;
        double mean_rs = (n > 0) ? a->sum_rs / n : NAN;
        double mean_z  = (n > 0) ? a->sum_z  / n : NAN;
        double mean_d  = (n > 0) ? (double) a->sum_depth / n : NAN;
        double rel_pct = (n_total_columns > 0)
                         ? 100.0 * n / (double) n_total_columns : 0.0;
        fprintf(fh,
            "%s\t%" PRIi64 "\t%" PRIi64 "\t%.4f\t"
            "%.4f\t%.4f\t%.4f\t%.4f\t"
            "%.4f\t%.4f\t%.4f\t%.4f\t"
            "%.2f\t%.4f\t%.4f\n",
            cm->buckets[i].name, cm->buckets[i].depth_from_root, a->n_cols, rel_pct,
            mean_rs,
            histogram_quantile(a->rs_hist, 0.5),
            histogram_quantile(a->rs_hist, 0.05),
            histogram_quantile(a->rs_hist, 0.95),
            mean_z,
            histogram_quantile(a->z_hist, 0.5),
            histogram_quantile(a->z_hist, 0.05),
            histogram_quantile(a->z_hist, 0.95),
            mean_d,
            (n > 0) ? (double) a->n_top10_global / n : NAN,
            (n > 0) ? (double) a->n_bot10_global / n : NAN);
    }
}

/////////////////////////////////////////////////////////////////////////////
// Driver
/////////////////////////////////////////////////////////////////////////////

static void usage(void) {
    fprintf(stderr, "taffy gerp-stats [options]\n");
    fprintf(stderr, "Per-clade statistics over the RS wig + depth wig from `taffy gerp`.\n");
    fprintf(stderr, "Output is a TSV summary; optionally also a per-column ranked wig.\n\n");
    fprintf(stderr, "-i --inputFile      : RS wig from `taffy gerp -o`   (required)\n");
    fprintf(stderr, "-D --depthFile      : Depth wig from `taffy gerp -D` (required for --depthCorrect=zscore)\n");
    fprintf(stderr, "-t --tree           : Newick tree (required)\n");
    fprintf(stderr, "-L --lineageRoots   : Comma-separated list of internal-node labels, OR @file with one per line\n");
    fprintf(stderr, "-o --outputFile     : TSV summary (default stdout)\n");
    fprintf(stderr, "   --outputWig FILE : Optional per-column ranked wig output\n");
    fprintf(stderr, "-c --useCompression : Bgzip the --outputWig (matches taffy view convention)\n");
    fprintf(stderr, "   --depthCorrect=MODE : none | sqrt | zscore   (default zscore)\n");
    fprintf(stderr, "   --emit=MODE         : percentile | zscore | both  (default percentile; only with --outputWig)\n");
    fprintf(stderr, "   --histBins N        : Histogram bin count        (default 1000)\n");
    fprintf(stderr, "   --histRange MIN:MAX : Score range                (default -10:10)\n");
    fprintf(stderr, "   --rsHistRange MIN:MAX : Per-clade RS-histogram range (default -10:5)\n");
    fprintf(stderr, "   --minBucketN N      : Min obs per depth bucket   (default 30)\n");
    fprintf(stderr, "   --maxDepth N        : Largest expected depth     (default 4096)\n");
    fprintf(stderr, "-T --threads N      : Worker threads + bgzf I/O (default 1)\n");
    fprintf(stderr, "-l --logLevel       : Set the log level\n");
    fprintf(stderr, "-h --help\n");
}

enum {
    OPT_DEPTH_CORRECT = 256,
    OPT_EMIT,
    OPT_HIST_BINS,
    OPT_HIST_RANGE,
    OPT_RS_HIST_RANGE,
    OPT_MIN_BUCKET_N,
    OPT_MAX_DEPTH,
    OPT_OUTPUT_WIG,
};

// Load lineage roots from `spec`: either "@file" (one label per line) or a
// comma-separated CSV.  Owned by the caller via the returned stSet.
static stSet *load_lineage_roots(const char *spec) {
    stSet *s = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    if (spec == NULL || *spec == '\0') return s;
    if (spec[0] == '@') {
        FILE *fh = fopen(spec + 1, "r");
        if (fh == NULL) {
            fprintf(stderr, "taffy gerp-stats: cannot open lineage roots file: %s\n", spec + 1);
            stSet_destruct(s);
            return NULL;
        }
        char line[1024];
        while (fgets(line, sizeof(line), fh)) {
            char *p = line;
            while (*p == ' ' || *p == '\t') p++;
            size_t n = strlen(p);
            while (n > 0 && (p[n-1] == '\n' || p[n-1] == '\r' ||
                              p[n-1] == ' '  || p[n-1] == '\t')) p[--n] = '\0';
            if (n == 0 || p[0] == '#') continue;
            if (stSet_search(s, p) == NULL) stSet_insert(s, stString_copy(p));
        }
        fclose(fh);
    } else {
        char *buf = stString_copy(spec);
        char *save = NULL;
        for (char *tok = strtok_r(buf, ",", &save); tok != NULL; tok = strtok_r(NULL, ",", &save)) {
            while (*tok == ' ' || *tok == '\t') tok++;
            size_t n = strlen(tok);
            while (n > 0 && (tok[n-1] == ' ' || tok[n-1] == '\t')) tok[--n] = '\0';
            if (n == 0) continue;
            if (stSet_search(s, tok) == NULL) stSet_insert(s, stString_copy(tok));
        }
        free(buf);
    }
    return s;
}

// Read a Newick file into a heap string (caller frees).
static char *slurp_text_file(const char *path) {
    FILE *fh = fopen(path, "r");
    if (fh == NULL) return NULL;
    fseek(fh, 0, SEEK_END);
    long n = ftell(fh);
    if (n < 0) { fclose(fh); return NULL; }
    fseek(fh, 0, SEEK_SET);
    char *buf = st_malloc((size_t) n + 1);
    size_t got = fread(buf, 1, (size_t) n, fh);
    buf[got] = '\0';
    fclose(fh);
    return buf;
}

typedef enum { DC_NONE, DC_SQRT, DC_ZSCORE } DepthCorrect;

int taf_gerp_stats_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    char *logLevelString = NULL;
    char *inputFile      = NULL;
    char *depthFile      = NULL;
    char *outputFile     = NULL;   // TSV summary
    char *outputWigFile  = NULL;   // optional per-column wig
    char *treeFile       = NULL;
    char *lineageRoots   = NULL;
    bool use_compression = false;
    DepthCorrect depth_correct = DC_ZSCORE;
    EmitMode emit_mode = EMIT_PERCENTILE;
    int64_t hist_bins = 1000;
    double  hist_min = -10.0, hist_max = 10.0;
    double  rs_hist_min = -10.0, rs_hist_max = 5.0;
    int64_t min_bucket_n = 30;
    int64_t max_depth = 4096;
    int n_threads = 1;
    int batch_size = 4096;

    while (1) {
        static struct option long_options[] = {
            { "logLevel",       required_argument, 0, 'l' },
            { "inputFile",      required_argument, 0, 'i' },
            { "depthFile",      required_argument, 0, 'D' },
            { "tree",           required_argument, 0, 't' },
            { "lineageRoots",   required_argument, 0, 'L' },
            { "outputFile",     required_argument, 0, 'o' },
            { "outputWig",      required_argument, 0, OPT_OUTPUT_WIG },
            { "useCompression", no_argument,       0, 'c' },
            { "depthCorrect",   required_argument, 0, OPT_DEPTH_CORRECT },
            { "emit",           required_argument, 0, OPT_EMIT },
            { "histBins",       required_argument, 0, OPT_HIST_BINS },
            { "histRange",      required_argument, 0, OPT_HIST_RANGE },
            { "rsHistRange",    required_argument, 0, OPT_RS_HIST_RANGE },
            { "minBucketN",     required_argument, 0, OPT_MIN_BUCKET_N },
            { "maxDepth",       required_argument, 0, OPT_MAX_DEPTH },
            { "threads",        required_argument, 0, 'T' },
            { "help",           no_argument,       0, 'h' },
            { 0, 0, 0, 0 }
        };
        int idx = 0;
        int64_t key = getopt_long(argc, argv, "l:i:D:t:L:o:cT:h", long_options, &idx);
        if (key == -1) break;
        switch (key) {
            case 'l': logLevelString = optarg; break;
            case 'i': inputFile     = optarg; break;
            case 'D': depthFile     = optarg; break;
            case 't': treeFile      = optarg; break;
            case 'L': lineageRoots  = optarg; break;
            case 'o': outputFile    = optarg; break;
            case OPT_OUTPUT_WIG: outputWigFile = optarg; break;
            case 'c': use_compression = true; break;
            case 'T': n_threads     = atoi(optarg); break;
            case OPT_DEPTH_CORRECT:
                if      (strcmp(optarg, "none")    == 0) depth_correct = DC_NONE;
                else if (strcmp(optarg, "sqrt")    == 0) depth_correct = DC_SQRT;
                else if (strcmp(optarg, "zscore")  == 0) depth_correct = DC_ZSCORE;
                else { fprintf(stderr, "taffy gerp-stats: unknown --depthCorrect %s\n", optarg); return 1; }
                break;
            case OPT_EMIT:
                if      (strcmp(optarg, "percentile") == 0) emit_mode = EMIT_PERCENTILE;
                else if (strcmp(optarg, "zscore")     == 0) emit_mode = EMIT_ZSCORE;
                else if (strcmp(optarg, "both")       == 0) emit_mode = EMIT_BOTH;
                else { fprintf(stderr, "taffy gerp-stats: unknown --emit %s\n", optarg); return 1; }
                break;
            case OPT_HIST_BINS:    hist_bins    = atol(optarg); break;
            case OPT_HIST_RANGE:
                {
                    char *colon = strchr(optarg, ':');
                    if (colon == NULL) {
                        fprintf(stderr, "taffy gerp-stats: --histRange requires MIN:MAX\n");
                        return 1;
                    }
                    *colon = '\0';
                    hist_min = atof(optarg);
                    hist_max = atof(colon + 1);
                }
                break;
            case OPT_RS_HIST_RANGE:
                {
                    char *colon = strchr(optarg, ':');
                    if (colon == NULL) {
                        fprintf(stderr, "taffy gerp-stats: --rsHistRange requires MIN:MAX\n");
                        return 1;
                    }
                    *colon = '\0';
                    rs_hist_min = atof(optarg);
                    rs_hist_max = atof(colon + 1);
                }
                break;
            case OPT_MIN_BUCKET_N: min_bucket_n = atol(optarg); break;
            case OPT_MAX_DEPTH:    max_depth    = atol(optarg); break;
            case 'h': usage(); return 0;
            default:  usage(); return 1;
        }
    }
    if (n_threads < 1) n_threads = 1;
    if (hist_bins < 2 || hist_max <= hist_min || rs_hist_max <= rs_hist_min) {
        fprintf(stderr, "taffy gerp-stats: invalid histogram config\n");
        return 1;
    }
    if (inputFile == NULL || treeFile == NULL || lineageRoots == NULL) {
        fprintf(stderr, "taffy gerp-stats: -i, -t, and -L are required\n");
        usage();
        return 1;
    }
    if (depth_correct == DC_ZSCORE && depthFile == NULL) {
        fprintf(stderr, "taffy gerp-stats: --depthCorrect=zscore needs -D <depth.wig>\n");
        return 1;
    }
    if (depth_correct == DC_SQRT && depthFile == NULL) {
        fprintf(stderr, "taffy gerp-stats: --depthCorrect=sqrt needs -D <depth.wig>\n");
        return 1;
    }

    st_setLogLevelFromString(logLevelString);
    LI_set_bgzf_threads(n_threads);

    // --- Load tree + lineage roots. -----------------------------------
    char *newick = slurp_text_file(treeFile);
    if (newick == NULL) {
        fprintf(stderr, "taffy gerp-stats: cannot read tree: %s\n", treeFile);
        return 1;
    }
    GerpTree *gt = gerp_tree_construct(newick);
    free(newick);
    if (gt == NULL) {
        fprintf(stderr, "taffy gerp-stats: failed to parse tree: %s\n", treeFile);
        return 1;
    }
    stSet *roots_set = load_lineage_roots(lineageRoots);
    if (roots_set == NULL) { gerp_tree_destruct(gt); return 1; }
    if (stSet_size(roots_set) == 0) {
        fprintf(stderr, "taffy gerp-stats: --lineageRoots resolved to no labels\n");
        stSet_destruct(roots_set); gerp_tree_destruct(gt);
        return 1;
    }
    // Validate that every root label exists in the tree.
    stList *rl = stSet_getList(roots_set);
    int n_root_missing = 0;
    for (int64_t i = 0; i < stList_length(rl); i++) {
        const char *r = stList_get(rl, i);
        if (gerp_tree_depth_from_root(gt, r) < 0) {
            fprintf(stderr, "taffy gerp-stats: lineage root '%s' not found in tree\n", r);
            n_root_missing++;
        }
    }
    stList_destruct(rl);
    if (n_root_missing > 0) {
        stSet_destruct(roots_set); gerp_tree_destruct(gt);
        return 1;
    }
    CladeMap *cm = clademap_construct(gt, roots_set);  // takes over roots_set
    st_logInfo("taffy gerp-stats: %" PRIi64 " user lineage roots, tree has "
               "%" PRIi64 " leaves\n",
               cm->n_buckets - 2, gerp_tree_n_leaves(gt));

    // --- Pass 1: collect per-depth stats. -----------------------------
    JointWig *j1 = joint_open(inputFile, (depth_correct != DC_NONE) ? depthFile : NULL);
    if (j1 == NULL) return 1;
    DepthStats *ds = depth_stats_construct(max_depth);
    int64_t n_total = 0;
    st_logInfo("taffy gerp-stats: pass 1 / 3 (per-depth stats)\n");
    int rc = pass1_collect_stats(j1, ds, n_threads, batch_size, &n_total);
    joint_close(j1);
    if (rc != 0) return 1;
    depth_stats_finalize(ds, min_bucket_n);
    st_logInfo("taffy gerp-stats: pass 1 done -- %" PRIi64 " columns in %" PRIi64 "s\n",
               n_total, (int64_t)(time(NULL) - startTime));
    if (ds->n_clamped > 0) {
        // The clamped observations get z=0 in pass 3, which silently
        // depresses the global histogram + per-clade tails.  Surface to
        // stderr so the user can re-run with a larger --maxDepth.  We
        // don't know the actual max depth seen (observe drops the value
        // on the floor), so the recommendation is just "bigger."
        double pct = (n_total > 0) ? 100.0 * (double)ds->n_clamped / (double)n_total : 0.0;
        fprintf(stderr,
                "taffy gerp-stats: WARNING -- %" PRIi64 " of %" PRIi64
                " columns (%.3f%%) had depth > --maxDepth=%" PRIi64
                " and were dropped from depth-correction (z=0). "
                "Re-run with a larger --maxDepth if those columns matter.\n",
                ds->n_clamped, n_total, pct, max_depth);
    }

    if (depth_correct == DC_NONE) {
        // With no depth correction, the "z-score" is the raw RS.  Wire
        // the DepthStats so that stddev=1, mean=0 across all buckets:
        for (int64_t d = 0; d <= ds->max_depth; d++) {
            ds->mean[d]        = 0.0;
            ds->stddev[d]      = 1.0;
            ds->fallback_to[d] = d;
        }
    } else if (depth_correct == DC_SQRT) {
        // sqrt: scale RS by sqrt(depth).  Reuse the zscore plumbing by
        // setting mean=0 and stddev = 1/sqrt(d) (then z = RS*sqrt(d)).
        for (int64_t d = 0; d <= ds->max_depth; d++) {
            ds->mean[d]        = 0.0;
            ds->stddev[d]      = (d > 0) ? 1.0 / sqrt((double) d) : 1.0;
            ds->fallback_to[d] = d;
        }
    }
    // DC_ZSCORE: leave finalize()'s output as-is.

    // --- Pass 2: histogram of corrected values. ------------------------
    JointWig *j2 = joint_open(inputFile, (depth_correct != DC_NONE) ? depthFile : NULL);
    if (j2 == NULL) { depth_stats_destruct(ds); return 1; }
    Histogram *h = histogram_construct(hist_bins, hist_min, hist_max);
    st_logInfo("taffy gerp-stats: pass 2 / 3 (z-score histogram, %" PRIi64
               " bins in [%g, %g))\n", hist_bins, hist_min, hist_max);
    rc = pass2_collect_histogram(j2, ds, h, n_threads, batch_size);
    joint_close(j2);
    if (rc != 0) { depth_stats_destruct(ds); histogram_destruct(h); return 1; }
    histogram_finalize(h);
    if (h->n_clipped_lo > 0 || h->n_clipped_hi > 0) {
        st_logInfo("taffy gerp-stats: %" PRIi64 " values clipped below %g, "
                   "%" PRIi64 " clipped at/above %g (of %" PRIi64 " total)\n",
                   h->n_clipped_lo, h->min_val, h->n_clipped_hi, h->max_val, h->n_total);
    }
    st_logInfo("taffy gerp-stats: pass 2 done in %" PRIi64 "s total\n",
               (int64_t)(time(NULL) - startTime));

    // Derive global top/bot thresholds for per-clade tail counting.
    double z_top10 = histogram_quantile(h, 0.90);
    double z_bot10 = histogram_quantile(h, 0.10);
    st_logInfo("taffy gerp-stats: global thresholds top10 z=%.4f, bot10 z=%.4f\n",
               z_top10, z_bot10);

    // --- Pass 3: per-clade accumulate (+ optional wig emit). ----------
    LW   *wig_out    = NULL;
    FILE *wig_out_fh = NULL;
    if (outputWigFile != NULL) {
        wig_out_fh = fopen(outputWigFile, "w");
        if (wig_out_fh == NULL) {
            fprintf(stderr, "taffy gerp-stats: cannot open --outputWig: %s\n", outputWigFile);
            depth_stats_destruct(ds); histogram_destruct(h);
            clademap_destruct(cm); gerp_tree_destruct(gt);
            return 1;
        }
        wig_out = LW_construct(wig_out_fh, use_compression);
    }
    JointWig *j3 = joint_open(inputFile, (depth_correct != DC_NONE) ? depthFile : NULL);
    if (j3 == NULL) {
        if (wig_out != NULL) { LW_destruct(wig_out, false); fclose(wig_out_fh); }
        depth_stats_destruct(ds); histogram_destruct(h);
        clademap_destruct(cm); gerp_tree_destruct(gt);
        return 1;
    }
    CladeAccum *clade_global = clade_accum_new_array(
        cm->n_buckets, hist_bins, rs_hist_min, rs_hist_max, hist_min, hist_max);
    st_logInfo("taffy gerp-stats: pass 3 / 3 (per-clade stats%s)\n",
               wig_out ? " + wig emit" : "");
    rc = pass3_run(j3, ds, h, cm, clade_global,
                   z_top10, z_bot10,
                   wig_out, emit_mode, n_threads, batch_size);
    joint_close(j3);
    if (wig_out != NULL) {
        LW_destruct(wig_out, false);
        fclose(wig_out_fh);
    }

    // --- TSV summary out. ---------------------------------------------
    FILE *tsv_fh = (outputFile == NULL) ? stdout : fopen(outputFile, "w");
    if (tsv_fh == NULL) {
        fprintf(stderr, "taffy gerp-stats: cannot open TSV output: %s\n", outputFile);
        clade_accum_destroy_array(clade_global, cm->n_buckets);
        depth_stats_destruct(ds); histogram_destruct(h);
        clademap_destruct(cm); gerp_tree_destruct(gt);
        return 1;
    }
    emit_tsv_summary(tsv_fh, cm, clade_global, n_total);
    if (outputFile != NULL) fclose(tsv_fh);

    st_logInfo("taffy gerp-stats: done in %" PRIi64 "s\n",
               (int64_t)(time(NULL) - startTime));

    clade_accum_destroy_array(clade_global, cm->n_buckets);
    depth_stats_destruct(ds);
    histogram_destruct(h);
    clademap_destruct(cm);
    gerp_tree_destruct(gt);
    return rc;
}

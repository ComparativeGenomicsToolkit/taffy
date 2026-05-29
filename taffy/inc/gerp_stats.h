#ifndef TAFFY_GERP_RANK_H_
#define TAFFY_GERP_RANK_H_

#include <stdint.h>
#include <stdbool.h>
#include "sonLib.h"
#include "line_iterator.h"

/*
 * Post-processing for taffy gerp output: depth-correct RS scores (z-score
 * per integer depth bucket against the empirical mean/stddev at that
 * depth), then map the corrected scores to percentile ranks via a binned
 * CDF.  See gerp_rank.c for the algorithm; taf_gerp_rank.c for the CLI.
 *
 * All structures here are pure data + math.  No I/O.  Thread-safety: the
 * per-thread accumulators (DepthStats / Histogram) are designed to be
 * cheaply merged after a parallel-for; see gerp_rank.c::depth_stats_merge
 * and histogram_merge.
 */

/////////////////////////////////////////////////////////////////////////////
// DepthStats: per-integer-depth count/sum/sum_sq tallies.
//
// In universal-MAF use cases the depth ranges from 2 to (n_leaves), so
// the natural store is a dense array indexed by depth.  We cap the bucket
// count at construction so workers can pre-allocate worst-case arrays.
//
// `mean(d)` and `stddev(d)` are derived from the tallies via
// depth_stats_finalize() (one shot, after pass 1 is done across all
// threads).  Empty buckets return NaN; the caller's depth-correction step
// is responsible for the under-populated-bucket fallback policy.
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    int64_t max_depth;     // size of the arrays below
    int64_t *count;        // [max_depth + 1]
    double  *sum;          // [max_depth + 1]
    double  *sum_sq;       // [max_depth + 1]
    double  *mean;         // [max_depth + 1]; filled by _finalize
    double  *stddev;       // [max_depth + 1]; filled by _finalize
    int64_t *fallback_to;  // [max_depth + 1]; -1 = use own stats; else use this depth's
} DepthStats;

/* Construct with capacity for depths 0..max_depth inclusive. */
DepthStats *depth_stats_construct(int64_t max_depth);
void depth_stats_destruct(DepthStats *ds);

/* Tally one observation: column has integer depth d and (raw) RS value rs. */
static inline void depth_stats_observe(DepthStats *ds, int64_t d, double rs) {
    if (d < 0 || d > ds->max_depth) return;  // silently clamp out-of-range
    ds->count[d]  += 1;
    ds->sum[d]    += rs;
    ds->sum_sq[d] += rs * rs;
}

/* Merge `src` into `dst` (used to combine per-thread accumulators). */
void depth_stats_merge(DepthStats *dst, const DepthStats *src);

/*
 * Compute per-bucket mean/stddev + the fallback table.  Under-populated
 * buckets (count < min_n) get their `fallback_to[d]` set to the nearest
 * depth with count >= min_n; depths with count >= min_n use their own
 * stats (fallback_to[d] == d).  If NO bucket meets min_n the fallback
 * points at the highest-count bucket and the caller will get noisy z's --
 * we log a warning rather than refuse.
 *
 * stddev is computed with the Bessel-corrected (n-1) denominator; buckets
 * with count == 1 get stddev = NaN (only the single observation fits;
 * fallback handles the routing).
 */
void depth_stats_finalize(DepthStats *ds, int64_t min_n);

/*
 * Convenience: depth-corrected z-score for one observation.  Uses
 * fallback_to[] to route to a populated bucket if d's own bucket was
 * under-populated.  Returns 0 if even the fallback target has stddev == 0
 * (or NaN), so the caller's downstream histogram never sees NaN.
 */
static inline double depth_stats_zscore(const DepthStats *ds, int64_t d, double rs) {
    if (d < 0 || d > ds->max_depth) return 0.0;
    int64_t f = ds->fallback_to[d];
    if (f < 0) return 0.0;
    double s = ds->stddev[f];
    if (!(s > 0.0)) return 0.0;
    return (rs - ds->mean[f]) / s;
}

/////////////////////////////////////////////////////////////////////////////
// Histogram: uniform-bin score counter + CDF / percentile lookup.
//
// Pass 2 fills this from z-scored values; pass 3 reads from it.  Bins are
// uniform over [min_val, max_val); out-of-range values clamp to the edge
// bin and increment an `n_clipped` counter the CLI surfaces in the log.
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    double   min_val;
    double   max_val;
    int64_t  n_bins;        // = bin_count
    int64_t *bin_count;     // [n_bins]
    int64_t *bin_cum;       // [n_bins]; filled by _finalize (running sum)
    int64_t  n_total;       // sum over bin_count
    int64_t  n_clipped_lo;  // values below min_val
    int64_t  n_clipped_hi;  // values >= max_val
} Histogram;

Histogram *histogram_construct(int64_t n_bins, double min_val, double max_val);
void histogram_destruct(Histogram *h);

/* Drop one value into its bin.  Clamps out-of-range to the edge bin and
 * bumps the clip counter for logging. */
static inline void histogram_observe(Histogram *h, double v) {
    if (v < h->min_val) { h->bin_count[0]++; h->n_clipped_lo++; return; }
    if (v >= h->max_val) { h->bin_count[h->n_bins - 1]++; h->n_clipped_hi++; return; }
    double t = (v - h->min_val) / (h->max_val - h->min_val);
    int64_t b = (int64_t)(t * (double)h->n_bins);
    if (b < 0) b = 0;
    if (b >= h->n_bins) b = h->n_bins - 1;
    h->bin_count[b]++;
}

/* Merge `src` into `dst` (per-thread accumulators).  Both must have
 * matching n_bins / min_val / max_val. */
void histogram_merge(Histogram *dst, const Histogram *src);

/* Compute cumulative counts; must be called once pass-2 is complete. */
void histogram_finalize(Histogram *h);

/*
 * Percentile lookup.  Returns a value in [0, 100].  The reported percentile
 * is the fraction of OBSERVED values <= this value's bin, expressed as a
 * percent.  Linear interpolation within the bin would tighten the answer
 * but is overkill at 1000 bins.
 */
double histogram_percentile(const Histogram *h, double v);

#endif /* TAFFY_GERP_RANK_H_ */

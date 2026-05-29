/*
 * Pure data / math for taffy gerp-rank.  Public contract in gerp_stats.h.
 */

#include "gerp_stats.h"
#include "sonLib.h"
#include <math.h>
#include <string.h>

/////////////////////////////////////////////////////////////////////////////
// DepthStats
/////////////////////////////////////////////////////////////////////////////

DepthStats *depth_stats_construct(int64_t max_depth) {
    DepthStats *ds = st_calloc(1, sizeof(DepthStats));
    ds->max_depth   = max_depth;
    size_t n        = (size_t)(max_depth + 1);
    ds->count       = st_calloc(n, sizeof(int64_t));
    ds->sum         = st_calloc(n, sizeof(double));
    ds->sum_sq      = st_calloc(n, sizeof(double));
    ds->mean        = st_calloc(n, sizeof(double));
    ds->stddev      = st_calloc(n, sizeof(double));
    ds->fallback_to = st_malloc(n * sizeof(int64_t));
    for (size_t i = 0; i < n; i++) ds->fallback_to[i] = -1;
    return ds;
}

void depth_stats_destruct(DepthStats *ds) {
    if (ds == NULL) return;
    free(ds->count);  free(ds->sum);   free(ds->sum_sq);
    free(ds->mean);   free(ds->stddev); free(ds->fallback_to);
    free(ds);
}

void depth_stats_merge(DepthStats *dst, const DepthStats *src) {
    assert(dst->max_depth == src->max_depth);
    for (int64_t d = 0; d <= dst->max_depth; d++) {
        dst->count[d]  += src->count[d];
        dst->sum[d]    += src->sum[d];
        dst->sum_sq[d] += src->sum_sq[d];
    }
}

void depth_stats_finalize(DepthStats *ds, int64_t min_n) {
    // Compute mean + Bessel-corrected stddev for every bucket.  Empty
    // buckets get NaN sentinels; populated buckets with count==1 get
    // stddev = NaN (variance undefined).
    int64_t best_d = -1;
    int64_t best_n = -1;
    for (int64_t d = 0; d <= ds->max_depth; d++) {
        int64_t n = ds->count[d];
        if (n == 0) { ds->mean[d] = NAN; ds->stddev[d] = NAN; continue; }
        double  s = ds->sum[d];
        double ss = ds->sum_sq[d];
        ds->mean[d] = s / (double)n;
        if (n >= 2) {
            // Sample variance via the two-moment form.  Could lose
            // precision for large mean^2 -- if the user reports issues
            // at huge depths we can switch to Welford's online form.
            double var = (ss - (s * s) / (double)n) / (double)(n - 1);
            if (var < 0.0) var = 0.0;  // numerical floor
            ds->stddev[d] = sqrt(var);
        } else {
            ds->stddev[d] = NAN;
        }
        if (n > best_n) { best_n = n; best_d = d; }
    }
    if (best_d < 0) {
        // No data at all -- leave fallback_to = -1 everywhere; the
        // zscore helper short-circuits to 0 in that case.
        return;
    }
    // Wire up fallbacks.  For each d:
    //   if count[d] >= min_n -> fallback_to[d] = d (use own).
    //   else -> nearest depth with count >= min_n (walk outward).
    //   no such depth -> fall back to best_d (highest-count bucket; may
    //                     still be under-populated, with a warning).
    for (int64_t d = 0; d <= ds->max_depth; d++) {
        if (ds->count[d] >= min_n) { ds->fallback_to[d] = d; continue; }
        // Walk outward looking for the nearest qualifying bucket.
        int64_t pick = -1;
        for (int64_t step = 1; step <= ds->max_depth; step++) {
            int64_t lo = d - step, hi = d + step;
            if (lo >= 0 && ds->count[lo] >= min_n) { pick = lo; break; }
            if (hi <= ds->max_depth && ds->count[hi] >= min_n) { pick = hi; break; }
            if (lo < 0 && hi > ds->max_depth) break;
        }
        ds->fallback_to[d] = (pick >= 0) ? pick : best_d;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Histogram
/////////////////////////////////////////////////////////////////////////////

Histogram *histogram_construct(int64_t n_bins, double min_val, double max_val) {
    assert(n_bins > 0);
    assert(max_val > min_val);
    Histogram *h = st_calloc(1, sizeof(Histogram));
    h->n_bins    = n_bins;
    h->min_val   = min_val;
    h->max_val   = max_val;
    h->bin_count = st_calloc((size_t)n_bins, sizeof(int64_t));
    h->bin_cum   = st_calloc((size_t)n_bins, sizeof(int64_t));
    return h;
}

void histogram_destruct(Histogram *h) {
    if (h == NULL) return;
    free(h->bin_count);
    free(h->bin_cum);
    free(h);
}

void histogram_merge(Histogram *dst, const Histogram *src) {
    assert(dst->n_bins == src->n_bins);
    assert(dst->min_val == src->min_val && dst->max_val == src->max_val);
    for (int64_t i = 0; i < dst->n_bins; i++) dst->bin_count[i] += src->bin_count[i];
    dst->n_clipped_lo += src->n_clipped_lo;
    dst->n_clipped_hi += src->n_clipped_hi;
    dst->n_total      += src->n_total;
}

void histogram_finalize(Histogram *h) {
    int64_t cum = 0;
    for (int64_t i = 0; i < h->n_bins; i++) {
        cum += h->bin_count[i];
        h->bin_cum[i] = cum;
    }
    h->n_total = cum;  // counts include clipped (they were assigned to the edge bin)
}

double histogram_quantile(const Histogram *h, double q) {
    if (h->n_total == 0) return NAN;
    if (q <= 0.0) return h->min_val;
    if (q >= 1.0) return h->max_val;
    int64_t target = (int64_t)(q * (double) h->n_total);
    for (int64_t b = 0; b < h->n_bins; b++) {
        if (h->bin_cum[b] >= target) {
            double bin_w = (h->max_val - h->min_val) / (double) h->n_bins;
            return h->min_val + (double)b * bin_w + bin_w * 0.5;
        }
    }
    return h->max_val;
}

double histogram_percentile(const Histogram *h, double v) {
    if (h->n_total == 0) return 0.0;
    int64_t b;
    if (v < h->min_val) {
        b = 0;
    } else if (v >= h->max_val) {
        b = h->n_bins - 1;
    } else {
        double t = (v - h->min_val) / (h->max_val - h->min_val);
        b = (int64_t)(t * (double)h->n_bins);
        if (b < 0) b = 0;
        if (b >= h->n_bins) b = h->n_bins - 1;
    }
    return 100.0 * (double)h->bin_cum[b] / (double)h->n_total;
}

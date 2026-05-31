/*
 * Unit tests for taffy/impl/gerp_stats.c.  Focused on the data math:
 * DepthStats (per-depth mean/stddev + fallback table) and Histogram
 * (CDF + percentile lookup).  The streaming wig + 3-pass driver is
 * exercised end-to-end by the test harness running taffy gerp ->
 * taffy gerp-stats against evolverMammals (separate test below).
 */

#include "CuTest.h"
#include "gerp_stats.h"
#include "sonLib.h"
#include <math.h>

static void close_to(CuTest *tc, double expected, double actual, double tol) {
    if (fabs(expected - actual) > tol) {
        char msg[256];
        snprintf(msg, sizeof(msg), "expected %.6f, got %.6f (tol %g)", expected, actual, tol);
        CuFail(tc, msg);
    }
}

/////////////////////////////////////////////////////////////////////////////
// DepthStats
/////////////////////////////////////////////////////////////////////////////

static void test_depth_stats_basic_moments(CuTest *tc) {
    DepthStats *ds = depth_stats_construct(10);
    // depth=5 gets {1, 2, 3} -> mean=2, var=1, stddev=1.
    depth_stats_observe(ds, 5, 1.0);
    depth_stats_observe(ds, 5, 2.0);
    depth_stats_observe(ds, 5, 3.0);
    depth_stats_finalize(ds, 2);
    close_to(tc, 2.0, ds->mean[5], 1e-9);
    close_to(tc, 1.0, ds->stddev[5], 1e-9);
    CuAssertIntEquals(tc, 5, (int) ds->fallback_to[5]);
    // empty bucket -> NaN mean/stddev, but fallback to a populated one.
    CuAssertTrue(tc, isnan(ds->mean[3]));
    CuAssertIntEquals(tc, 5, (int) ds->fallback_to[3]);
    depth_stats_destruct(ds);
}

static void test_depth_stats_min_n_fallback(CuTest *tc) {
    DepthStats *ds = depth_stats_construct(10);
    // depth=3 gets 1 observation, depth=8 gets 30.  min_n=10 -> depth 3
    // should fall back to depth 8.
    depth_stats_observe(ds, 3, 7.0);
    for (int i = 0; i < 30; i++) depth_stats_observe(ds, 8, 1.0 * i);
    depth_stats_finalize(ds, 10);
    CuAssertIntEquals(tc, 8, (int) ds->fallback_to[3]);
    CuAssertIntEquals(tc, 8, (int) ds->fallback_to[8]);
    depth_stats_destruct(ds);
}

static void test_depth_stats_zscore_against_fallback(CuTest *tc) {
    DepthStats *ds = depth_stats_construct(10);
    // depth=5 is populated; depth=3 is sparse and falls back to 5.
    for (int i = 0; i < 40; i++) depth_stats_observe(ds, 5, (double) i);
    depth_stats_observe(ds, 3, 100.0);
    depth_stats_finalize(ds, 10);
    double z_at_5 = depth_stats_zscore(ds, 5, 19.5);  // (19.5 - mean) / stddev
    double mean5 = ds->mean[5], sd5 = ds->stddev[5];
    close_to(tc, (19.5 - mean5) / sd5, z_at_5, 1e-9);
    // depth=3 falls back to 5: (100 - mean5) / sd5.
    double z_at_3 = depth_stats_zscore(ds, 3, 100.0);
    close_to(tc, (100.0 - mean5) / sd5, z_at_3, 1e-9);
    depth_stats_destruct(ds);
}

static void test_depth_stats_merge(CuTest *tc) {
    DepthStats *a = depth_stats_construct(10);
    DepthStats *b = depth_stats_construct(10);
    depth_stats_observe(a, 4, 10.0);
    depth_stats_observe(a, 4, 20.0);
    depth_stats_observe(b, 4, 30.0);
    depth_stats_merge(a, b);
    depth_stats_finalize(a, 2);
    // {10, 20, 30}: mean=20, stddev=10.
    close_to(tc, 20.0, a->mean[4], 1e-9);
    close_to(tc, 10.0, a->stddev[4], 1e-9);
    depth_stats_destruct(a);
    depth_stats_destruct(b);
}

/////////////////////////////////////////////////////////////////////////////
// Histogram
/////////////////////////////////////////////////////////////////////////////

static void test_histogram_uniform(CuTest *tc) {
    // 10 bins over [-5, +5).  Drop values at the bin centres.
    Histogram *h = histogram_construct(10, -5.0, 5.0);
    for (int b = 0; b < 10; b++) {
        double v = -5.0 + b + 0.5;
        for (int i = 0; i < 100; i++) histogram_observe(h, v);
    }
    histogram_finalize(h);
    CuAssertIntEquals(tc, 1000, (int) h->n_total);
    // 50th percentile sits at the boundary between bins 4 and 5: 50% of
    // observations are <= bin 4 (count 500).  Value -5 + 4.999 lands in
    // bin 4 -> percentile 50.
    close_to(tc, 50.0, histogram_percentile(h, -0.001), 0.1);
    // Bin 9 covers [4, 5).  Anything inside -> 100th percentile (all
    // 1000 obs <= this bin's right edge).
    close_to(tc, 100.0, histogram_percentile(h, 4.5), 0.1);
    // Value below min -> bin 0 -> percentile = 100 * 100 / 1000 = 10.
    close_to(tc, 10.0, histogram_percentile(h, -100.0), 0.1);
    histogram_destruct(h);
}

static void test_histogram_clipping(CuTest *tc) {
    Histogram *h = histogram_construct(10, 0.0, 1.0);
    histogram_observe(h, -1.0);   // clipped low
    histogram_observe(h, 0.5);    // inside
    histogram_observe(h, 99.0);   // clipped high
    histogram_finalize(h);
    CuAssertIntEquals(tc, 1, (int) h->n_clipped_lo);
    CuAssertIntEquals(tc, 1, (int) h->n_clipped_hi);
    CuAssertIntEquals(tc, 3, (int) h->n_total);
    histogram_destruct(h);
}

static void test_histogram_merge(CuTest *tc) {
    Histogram *a = histogram_construct(4, 0.0, 4.0);
    Histogram *b = histogram_construct(4, 0.0, 4.0);
    histogram_observe(a, 0.5);
    histogram_observe(a, 1.5);
    histogram_observe(b, 2.5);
    histogram_observe(b, 3.5);
    histogram_merge(a, b);
    histogram_finalize(a);
    CuAssertIntEquals(tc, 4, (int) a->n_total);
    // 50th percentile should fall at the boundary between bins 1 and 2.
    close_to(tc, 50.0, histogram_percentile(a, 1.99), 0.1);
    close_to(tc, 100.0, histogram_percentile(a, 3.99), 0.1);
    histogram_destruct(a);
    histogram_destruct(b);
}

CuSuite* gerp_stats_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_depth_stats_basic_moments);
    SUITE_ADD_TEST(suite, test_depth_stats_min_n_fallback);
    SUITE_ADD_TEST(suite, test_depth_stats_zscore_against_fallback);
    SUITE_ADD_TEST(suite, test_depth_stats_merge);
    SUITE_ADD_TEST(suite, test_histogram_uniform);
    SUITE_ADD_TEST(suite, test_histogram_clipping);
    SUITE_ADD_TEST(suite, test_histogram_merge);
    return suite;
}

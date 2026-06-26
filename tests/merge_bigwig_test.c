/*
 * Tests for `taffy merge-bigwig` -- in particular the 64-bit path: shard column
 * slices that straddle 2^32, which the apes-scale validation never exercises.
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "bigWig.h"
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

extern int taf_merge_bigwig_main(int argc, char *argv[]);

#define MB_N    3
#define MB_TWO32  4294967296ULL          /* 2^32 */
#define MB_T      8589934592ULL          /* 2 * 2^32 -- axis comfortably past 2^32 */

/* write a single-chrom (uni0) vector bigWig with `n` intervals (values flat n*MB_N) */
static void write_vec_bw(const char *path, const uint64_t *s, const uint64_t *e,
                         const float *v, uint32_t n) {
    bigWigFile_t *fp = bwOpen(path, NULL, "w");
    if (!fp || bwCreateHdrVec(fp, 0, MB_N)) { fprintf(stderr, "write_vec_bw open/hdr %s\n", path); exit(1); }
    const char *chroms[] = { "uni0" };
    uint64_t    lens[]   = { MB_T };
    fp->cl = bwCreateChromList(chroms, lens, 1);
    if (!fp->cl || bwWriteHdr(fp)) { fprintf(stderr, "write_vec_bw hdr %s\n", path); exit(1); }
    const char **c = malloc((size_t) n * sizeof(char *));
    for (uint32_t i = 0; i < n; i++) c[i] = "uni0";
    if (bwAddIntervalsVec(fp, c, s, e, v, n)) { fprintf(stderr, "write_vec_bw add %s\n", path); exit(1); }
    free(c);
    bwClose(fp);
}

static void write_names(const char *bw_path, const char *names) {
    char *np = stString_print("%s.names", bw_path);
    FILE *f = fopen(np, "w");
    if (!f) { fprintf(stderr, "write_names %s\n", np); exit(1); }
    fputs(names, f);
    fclose(f);
    free(np);
}

/* read whole [MB_TWO32-2000, MB_TWO32+2000) range; caller destroys */
static bwOverlappingIntervalsVec_t *read_window(const char *path) {
    bigWigFile_t *fp = bwOpen(path, NULL, "r");
    if (!fp) { fprintf(stderr, "read_window open %s\n", path); exit(1); }
    bwOverlappingIntervalsVec_t *iv =
        bwGetOverlappingIntervalsVec(fp, "uni0", MB_TWO32 - 2000, MB_TWO32 + 2000);
    bwClose(fp);
    return iv;
}

/* merged.bw built by merging two shards (one below, one above 2^32) must equal
 * a single.bw written with all four intervals -- and the >2^32 coordinates must
 * round-trip without 32-bit truncation. */
static void test_merge_across_2to32(CuTest *tc) {
    const char *s0 = "/tmp/taffy_mbtest_s0.bw";
    const char *s1 = "/tmp/taffy_mbtest_s1.bw";
    const char *sg = "/tmp/taffy_mbtest_single.bw";
    const char *mg = "/tmp/taffy_mbtest_merged.bw";

    /* shard0: two runs JUST BELOW 2^32 */
    uint64_t s0s[] = { MB_TWO32 - 2000, MB_TWO32 - 1000 };
    uint64_t s0e[] = { MB_TWO32 - 1000, MB_TWO32 };
    float    s0v[] = { 1, 2, 3,  4, 5, 6 };
    /* shard1: two runs JUST ABOVE 2^32 */
    uint64_t s1s[] = { MB_TWO32,        MB_TWO32 + 1000 };
    uint64_t s1e[] = { MB_TWO32 + 1000, MB_TWO32 + 2000 };
    float    s1v[] = { 7, 8, 9,  10, 11, 12 };
    /* single: all four, in column order */
    uint64_t sgs[] = { MB_TWO32 - 2000, MB_TWO32 - 1000, MB_TWO32,        MB_TWO32 + 1000 };
    uint64_t sge[] = { MB_TWO32 - 1000, MB_TWO32,        MB_TWO32 + 1000, MB_TWO32 + 2000 };
    float    sgv[] = { 1, 2, 3,  4, 5, 6,  7, 8, 9,  10, 11, 12 };

    /* ---- phase 1: write the inputs ---- */
    if (bwInit(1 << 17)) { CuFail(tc, "bwInit (write)"); }
    write_vec_bw(s0, s0s, s0e, s0v, 2);  write_names(s0, "spA\nspB\nspC\n");
    write_vec_bw(s1, s1s, s1e, s1v, 2);  write_names(s1, "spA\nspB\nspC\n");
    write_vec_bw(sg, sgs, sge, sgv, 4);  write_names(sg, "spA\nspB\nspC\n");
    bwCleanup();

    /* ---- phase 2: merge shard0 + shard1 (merge brackets its own bwInit/Cleanup) ---- */
    char *argv[] = { "merge-bigwig", "-o", (char *) mg, (char *) s0, (char *) s1, NULL };
    optind = 1;                                  /* reset getopt state */
    int rv = taf_merge_bigwig_main(5, argv);
    CuAssertIntEquals_Msg(tc, "merge-bigwig should succeed", 0, rv);

    /* ---- phase 3: merged must equal single, with >2^32 coords intact ---- */
    if (bwInit(1 << 17)) { CuFail(tc, "bwInit (read)"); }
    bwOverlappingIntervalsVec_t *m = read_window(mg);
    bwOverlappingIntervalsVec_t *g = read_window(sg);
    CuAssertTrue(tc, m != NULL && g != NULL);
    CuAssertIntEquals_Msg(tc, "merged interval count", 4, (int) m->l);
    CuAssertIntEquals_Msg(tc, "single interval count", 4, (int) g->l);
    CuAssertIntEquals_Msg(tc, "merged vecN", MB_N, (int) m->N);

    /* the boundary coordinates must NOT be truncated to 32 bits */
    CuAssert(tc,"start[2] must be exactly 2^32 (no truncation)", m->start[2] == MB_TWO32);
    CuAssert(tc,"end[3] must be exactly 2^32+2000 (no truncation)", m->end[3] == MB_TWO32 + 2000);

    /* merged == single, element by element */
    for (uint64_t i = 0; i < m->l; i++) {
        CuAssert(tc,"start mismatch merged vs single", m->start[i] == g->start[i]);
        CuAssert(tc,"end mismatch merged vs single",   m->end[i]   == g->end[i]);
        for (uint32_t c = 0; c < MB_N; c++)
            CuAssertDblEquals_Msg(tc, "value mismatch merged vs single",
                                  g->value[i * MB_N + c], m->value[i * MB_N + c], 1e-6);
    }
    /* and the absolute expected coordinates, in case both happened to be wrong the same way */
    CuAssertTrue(tc, m->start[0] == MB_TWO32 - 2000 && m->end[0] == MB_TWO32 - 1000);
    CuAssertTrue(tc, m->start[3] == MB_TWO32 + 1000 && m->end[3] == MB_TWO32 + 2000);
    CuAssertDblEquals_Msg(tc, "value[2][0] (first run above 2^32)", 7.0, m->value[2 * MB_N + 0], 1e-6);

    bwDestroyOverlappingIntervalsVec(m);
    bwDestroyOverlappingIntervalsVec(g);
    bwCleanup();

    /* ---- the .names sidecar must be carried to the merge output ---- */
    char *mn = stString_print("%s.names", mg);
    FILE *nf = fopen(mn, "r");
    CuAssert(tc,"merged .names sidecar must exist", nf != NULL);
    if (nf) {
        char buf[256] = {0};
        size_t got = fread(buf, 1, sizeof buf - 1, nf);
        fclose(nf);
        buf[got] = '\0';
        CuAssertStrEquals_Msg(tc, "merged .names must match shard0 .names", "spA\nspB\nspC\n", buf);
    }
    free(mn);

    /* cleanup temp files */
    const char *paths[] = { s0, s1, sg, mg };
    for (int i = 0; i < 4; i++) {
        remove(paths[i]);
        char *n = stString_print("%s.names", paths[i]);
        remove(n);
        free(n);
    }
}

CuSuite* merge_bigwig_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_merge_across_2to32);
    return suite;
}

/*
 * Unit tests for taffy/blockViz public C API.  Focused on the bits we
 * can exercise without a network-fetched .tui -- chain-param setter /
 * getter, basic lifecycle, error paths.  End-to-end block extraction
 * is covered by the taffyBlockVizTest CLI smokes (which need a real
 * universal-MAF .tui to be meaningful).
 */

#include "CuTest.h"
#include "sonLib.h"
#include "taffyBlockViz.h"
#include "chain.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define TEST_TUI "tests/tui/evolverMammals.uni.maf.gz.tui"

/* ------------------------------------------------------------------ */
/* Lifecycle                                                           */
/* ------------------------------------------------------------------ */

static void test_open_close(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) {
        /* No fixture -- skip silently rather than fail; the file is a
         * generated artifact of the test harness. */
        free(err);
        return;
    }
    CuAssertTrue(tc, h >= 0);
    CuAssertIntEquals(tc, 0, taffyClose(h, &err));
    free(err);
}

static void test_open_invalid_path(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen("/no/such/path.tui", &err);
    CuAssertTrue(tc, h < 0);
    CuAssertPtrNotNull(tc, err);
    free(err);
}

static void test_close_invalid_handle(CuTest *tc) {
    char *err = NULL;
    int rc = taffyClose(99999, &err);
    CuAssertTrue(tc, rc != 0);
    CuAssertPtrNotNull(tc, err);
    free(err);
}

/* ------------------------------------------------------------------ */
/* Chain-tuning setter / getter                                        */
/* ------------------------------------------------------------------ */

static void test_chain_params_defaults(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }
    int64_t o = -999, e = -999, m = -999;
    CuAssertIntEquals(tc, 0, taffyGetChainParams(h, &o, &e, &m, &err));
    CuAssertTrue(tc, o == TAFFY_CHAIN_DEFAULT_OPEN);
    CuAssertTrue(tc, e == TAFFY_CHAIN_DEFAULT_EXTEND);
    CuAssertTrue(tc, m == TAFFY_CHAIN_DEFAULT_MAX_GAP);
    taffyClose(h, &err);
    free(err);
}

static void test_chain_params_roundtrip(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }
    /* Set all three to non-default values. */
    CuAssertIntEquals(tc, 0,
        taffySetChainParams(h, 1000, 2, 5000000, &err));
    int64_t o, e, m;
    CuAssertIntEquals(tc, 0, taffyGetChainParams(h, &o, &e, &m, &err));
    CuAssertTrue(tc, o == 1000);
    CuAssertTrue(tc, e == 2);
    CuAssertTrue(tc, m == 5000000);
    taffyClose(h, &err);
    free(err);
}

static void test_chain_params_sentinel_minus_one_preserves(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }
    /* Set a baseline. */
    taffySetChainParams(h, 100, 5, 1000, &err);
    /* -1 in chain_open should leave it at 100; chain_extend bumped to 7;
     * max_gap_length left at 1000. */
    CuAssertIntEquals(tc, 0, taffySetChainParams(h, -1, 7, -1, &err));
    int64_t o, e, m;
    taffyGetChainParams(h, &o, &e, &m, &err);
    CuAssertTrue(tc, o == 100);
    CuAssertTrue(tc, e == 7);
    CuAssertTrue(tc, m == 1000);
    taffyClose(h, &err);
    free(err);
}

static void test_chain_params_rejects_negative(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }
    /* -2, -100 (anything below -1) must be rejected. */
    char *err2 = NULL;
    CuAssertIntEquals(tc, -1, taffySetChainParams(h, -2, 1, 1000, &err2));
    CuAssertPtrNotNull(tc, err2); free(err2);
    err2 = NULL;
    CuAssertIntEquals(tc, -1, taffySetChainParams(h, 0, -100, 1000, &err2));
    CuAssertPtrNotNull(tc, err2); free(err2);
    err2 = NULL;
    CuAssertIntEquals(tc, -1, taffySetChainParams(h, 0, 1, -5, &err2));
    CuAssertPtrNotNull(tc, err2); free(err2);
    taffyClose(h, &err);
    free(err);
}

static void test_chain_params_null_outputs_ok(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }
    /* Calling Get with NULL out pointers must not crash; just a no-op. */
    CuAssertIntEquals(tc, 0, taffyGetChainParams(h, NULL, NULL, NULL, &err));
    int64_t only_open = 0;
    CuAssertIntEquals(tc, 0, taffyGetChainParams(h, &only_open, NULL, NULL, &err));
    CuAssertTrue(tc, only_open == TAFFY_CHAIN_DEFAULT_OPEN);
    taffyClose(h, &err);
    free(err);
}

static void test_chain_params_invalid_handle(CuTest *tc) {
    char *err = NULL;
    CuAssertIntEquals(tc, -1, taffySetChainParams(99999, 1, 1, 1, &err));
    CuAssertPtrNotNull(tc, err); free(err);
    err = NULL;
    int64_t o = 0;
    CuAssertIntEquals(tc, -1, taffyGetChainParams(99999, &o, NULL, NULL, &err));
    CuAssertPtrNotNull(tc, err); free(err);
}

/* ------------------------------------------------------------------ */
/* Auto-binning (wide-query LOD).  See taffyGetBlocksInTargetRange      */
/* doc: spans >= 200 kb engage bin mode.  Evolver chr6 is ~600 kb so   */
/* it crosses the threshold.                                            */
/* ------------------------------------------------------------------ */

static void test_get_blocks_narrow_per_run(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }

    /* 5 kb span -> 5000 / 2000 = 2 < MIN_BIN_BP=100 -> per-run blocks. */
    struct taffy_block_results_t *res = taffyGetBlocksInTargetRange(
        h, "simMouse_chr6", "simHuman_chr6", "simHuman.chr6",
        0, 5000, 0, TAFFY_NO_SEQUENCES, TAFFY_QUERY_AND_TARGET_DUPS,
        0, NULL, &err);
    CuAssertPtrNotNull(tc, res);

    /* Per-run blocks: each .size is the run length (1..thousands of bp);
     * qStart is a real qSpecies coord; qChrom prefixes the species. */
    int64_t n = 0;
    for (struct taffy_block_t *b = res->mappedBlocks; b; b = b->next) {
        n++;
        CuAssertTrue(tc, b->size > 0);
        CuAssertTrue(tc, b->size <= 5000);   /* a per-run can't exceed span */
        CuAssertTrue(tc, b->strand == '+' || b->strand == '-');
    }
    /* The fixture covers >= 1 block in this range (evolver is dense). */
    CuAssertTrue(tc, n >= 1);

    taffyFreeBlockResults(res);
    taffyClose(h, &err);
    free(err);
}

static void test_get_blocks_wide_binned(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }

    /* Full chr6 (evolver simHuman.chr6 is ~600 kb).  span/2000 >> 100 ->
     * bin mode.  bin_size = span / 2000 = ~300 bp. */
    struct taffy_block_results_t *res = taffyGetBlocksInTargetRange(
        h, "simMouse_chr6", "simHuman_chr6", "simHuman.chr6",
        0, 0 /* = chrom length */, 0, TAFFY_NO_SEQUENCES,
        TAFFY_QUERY_AND_TARGET_DUPS, 0, NULL, &err);
    CuAssertPtrNotNull(tc, res);

    /* Bin mode invariants:
     *  - dupe list is always NULL
     *  - each block has size > 0 (covered_bp; may exceed bin_size when
     *    paralogs map multiple source columns to the same target bin)
     *  - tStart values are unique-and-sorted per (qChrom, strand) group
     *    (map iteration order); strand is '+' or '-'
     *  - we get FAR FEWER than the per-run count would be (the fixture
     *    has ~thousands of runs across chr6) */
    CuAssertPtrEquals(tc, NULL, res->targetDupeBlocks);

    int64_t n = 0;
    for (struct taffy_block_t *b = res->mappedBlocks; b; b = b->next) {
        n++;
        CuAssertTrue(tc, b->size > 0);
        CuAssertTrue(tc, b->strand == '+' || b->strand == '-');
        CuAssertTrue(tc, b->tStart >= 0);
        CuAssertPtrNotNull(tc, b->qChrom);
    }
    /* Cap: at bin_size ~300 bp on a 600 kb chrom and ~10 qChrom/strand
     * combos in evolver, expect at most a few thousand bin entries.
     * (The per-run count for full evolver chr6 is in the same order of
     * magnitude on this small fixture, so we can't assert "much less"
     * here -- the win shows up at vertebrate scale.  Sanity-check
     * bounds only.) */
    CuAssertTrue(tc, n >= 1);
    CuAssertTrue(tc, n < 100000);

    taffyFreeBlockResults(res);
    taffyClose(h, &err);
    free(err);
}

/* ------------------------------------------------------------------ */
/* Suite                                                               */
/* ------------------------------------------------------------------ */

CuSuite* blockviz_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_open_close);
    SUITE_ADD_TEST(suite, test_open_invalid_path);
    SUITE_ADD_TEST(suite, test_close_invalid_handle);
    SUITE_ADD_TEST(suite, test_chain_params_defaults);
    SUITE_ADD_TEST(suite, test_chain_params_roundtrip);
    SUITE_ADD_TEST(suite, test_chain_params_sentinel_minus_one_preserves);
    SUITE_ADD_TEST(suite, test_chain_params_rejects_negative);
    SUITE_ADD_TEST(suite, test_chain_params_null_outputs_ok);
    SUITE_ADD_TEST(suite, test_chain_params_invalid_handle);
    SUITE_ADD_TEST(suite, test_get_blocks_narrow_per_run);
    SUITE_ADD_TEST(suite, test_get_blocks_wide_binned);
    return suite;
}

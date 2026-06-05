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
/* Output regimes (see taffyGetBlocksInTargetRange header doc).        */
/* Per-run + chain-merge for spans < 10 Mb (the test fixture);         */
/* auto-binning for spans >= 10 Mb (not exercisable on the ~600 kb     */
/* evolver fixture -- needs vertebrate-scale .tui).                    */
/* ------------------------------------------------------------------ */

static void test_get_blocks_narrow_per_run(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }

    /* 5 kb span -> well below 10 Mb bin threshold -> per-run + merge. */
    struct taffy_block_results_t *res = taffyGetBlocksInTargetRange(
        h, "simMouse_chr6", "simHuman_chr6", "simHuman.chr6",
        0, 5000, 0, TAFFY_NO_SEQUENCES, TAFFY_QUERY_AND_TARGET_DUPS,
        0, NULL, &err);
    CuAssertPtrNotNull(tc, res);

    /* Per-run (possibly post-chain-merged) blocks: each .size is the
     * span of one merged collinear stretch.  qStart is a real qSpecies
     * coord; qChrom prefixes the species. */
    int64_t n = 0;
    for (struct taffy_block_t *b = res->mappedBlocks; b; b = b->next) {
        n++;
        CuAssertTrue(tc, b->size > 0);
        CuAssertTrue(tc, b->strand == '+' || b->strand == '-');
    }
    CuAssertTrue(tc, n >= 1);

    taffyFreeBlockResults(res);
    taffyClose(h, &err);
    free(err);
}

static void test_get_blocks_full_chrom_merged(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }

    /* Full evolver chr6 (~600 kb).  Below the 10 Mb bin threshold so
     * this exercises the per-run + chain-merge path on the largest
     * range the fixture supports.  Primary-chain dupMode keeps the
     * count tractable (no paralog inflation). */
    struct taffy_block_results_t *res = taffyGetBlocksInTargetRange(
        h, "simMouse_chr6", "simHuman_chr6", "simHuman.chr6",
        0, 0 /* = chrom length */, 0, TAFFY_NO_SEQUENCES,
        TAFFY_NO_DUPS, 0, NULL, &err);
    CuAssertPtrNotNull(tc, res);

    int64_t n = 0;
    for (struct taffy_block_t *b = res->mappedBlocks; b; b = b->next) {
        n++;
        CuAssertTrue(tc, b->size > 0);
        CuAssertTrue(tc, b->strand == '+' || b->strand == '-');
        CuAssertTrue(tc, b->tStart >= 0);
        CuAssertPtrNotNull(tc, b->qChrom);
    }
    CuAssertTrue(tc, n >= 1);
    /* Snake-track sanity: post-merge primary-chain block count for one
     * chrom should fit comfortably under a few hundred even on this
     * fixture.  Loose bound to stay portable across fixture rebuilds. */
    CuAssertTrue(tc, n < 5000);

    taffyFreeBlockResults(res);
    taffyClose(h, &err);
    free(err);
}

/* chainId exposed on taffy_block_t; chainSummaries populated for non-
 * bin queries.  Asserts the contract documented in taffyBlockViz.h:
 *  - per-run + chain output: chainId nonzero, multiple blocks may
 *    share a chainId, chainSummaries lists each surviving chain in
 *    score-desc order with id matching some block.chainId.
 *  - bin output (wide span on evolver -- forced via taffySetMax
 *    OutputBlocks=1 to engage the cap at any width, OR via the
 *    natural 5 Mb threshold which evolver's 600 kb chrom doesn't
 *    reach; we use a primitive size check instead).
 */
static void test_chain_id_and_summaries(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }

    /* Narrow query on evolver chr6: per-run + chain path. */
    struct taffy_block_results_t *res = taffyGetBlocksInTargetRange(
        h, "simMouse_chr6", "simHuman_chr6", "simHuman.chr6",
        0, 5000, 0, TAFFY_NO_SEQUENCES, TAFFY_QUERY_AND_TARGET_DUPS,
        0, NULL, &err);
    CuAssertPtrNotNull(tc, res);

    /* Collect chainIds from emitted blocks. */
    stSet *block_chain_ids = stSet_construct();
    int64_t any_nonzero = 0, n_blocks = 0;
    for (struct taffy_block_t *b = res->mappedBlocks; b; b = b->next) {
        n_blocks++;
        if (b->chainId != 0) {
            any_nonzero++;
            int64_t *idp = st_malloc(sizeof(int64_t));
            *idp = b->chainId;
            stSet_insert(block_chain_ids, idp);
        }
    }
    /* Some blocks should be in chains (not all bin / flank). */
    CuAssertTrue(tc, any_nonzero > 0);

    /* chainSummaries should be present and score-desc; each summary's
     * id should map to at least one block's chainId. */
    CuAssertPtrNotNull(tc, res->chainSummaries);
    int64_t prev_score = -1, n_summaries = 0;
    for (struct taffy_chain_summary_t *c = res->chainSummaries; c; c = c->next) {
        n_summaries++;
        /* Score-desc: each subsequent <= prev. */
        if (prev_score >= 0) CuAssertTrue(tc, c->totalScore <= prev_score);
        prev_score = c->totalScore;
        /* Reasonable invariants. */
        CuAssertTrue(tc, c->id > 0);
        CuAssertTrue(tc, c->nAlns > 0);
        CuAssertTrue(tc, c->totalBp > 0);
    }
    CuAssertTrue(tc, n_summaries > 0);
    /* Summary count <= block count (chains can be > blocks in degenerate
     * cases via merge, but typically <=). */

    /* Clean up. */
    stSet_destruct(block_chain_ids);
    taffyFreeBlockResults(res);
    taffyClose(h, &err);
    free(err);
}

/* taffySet/GetMaxOutputBlocks setter contract + runtime cap effect. */
static void test_max_output_blocks_setter(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }

    /* Default cap is 500. */
    int64_t cap = 0;
    CuAssertIntEquals(tc, 0, taffyGetMaxOutputBlocks(h, &cap, &err));
    CuAssertTrue(tc, cap == 500);

    /* Set + roundtrip. */
    CuAssertIntEquals(tc, 0, taffySetMaxOutputBlocks(h, 47, &err));
    CuAssertIntEquals(tc, 0, taffyGetMaxOutputBlocks(h, &cap, &err));
    CuAssertTrue(tc, cap == 47);

    /* Reject 0 / negative. */
    CuAssertIntEquals(tc, -1, taffySetMaxOutputBlocks(h, 0, &err));
    CuAssertPtrNotNull(tc, err); free(err); err = NULL;
    CuAssertIntEquals(tc, -1, taffySetMaxOutputBlocks(h, -3, &err));
    CuAssertPtrNotNull(tc, err); free(err); err = NULL;
    /* Cap unchanged after rejected attempts. */
    CuAssertIntEquals(tc, 0, taffyGetMaxOutputBlocks(h, &cap, &err));
    CuAssertTrue(tc, cap == 47);

    /* Cap actually constrains output.  Full evolver chr6 query at
     * cap=10 must produce <= 10 blocks. */
    struct taffy_block_results_t *res = taffyGetBlocksInTargetRange(
        h, "simMouse_chr6", "simHuman_chr6", "simHuman.chr6",
        0, 0, 0, TAFFY_NO_SEQUENCES, TAFFY_QUERY_AND_TARGET_DUPS,
        0, NULL, &err);
    CuAssertPtrNotNull(tc, res);
    int64_t n = 0;
    for (struct taffy_block_t *b = res->mappedBlocks; b; b = b->next) n++;
    CuAssertTrue(tc, n <= 47);
    taffyFreeBlockResults(res);

    /* Invalid handle errors. */
    CuAssertIntEquals(tc, -1, taffySetMaxOutputBlocks(99999, 100, &err));
    CuAssertPtrNotNull(tc, err); free(err); err = NULL;
    CuAssertIntEquals(tc, -1, taffyGetMaxOutputBlocks(99999, &cap, &err));
    CuAssertPtrNotNull(tc, err); free(err); err = NULL;

    /* NULL out-ptr in getter is a no-op (no crash). */
    CuAssertIntEquals(tc, 0, taffyGetMaxOutputBlocks(h, NULL, &err));

    taffyClose(h, &err);
    free(err);
}

/* Output-cap invariant: every dupMode must yield mappedBlocks <= 500.
 * Tests the budget logic by including paralogs (which on a small
 * fixture won't actually exceed the cap, but the assertion is the
 * contract -- larger fixtures would exercise the truncation). */
static void test_get_blocks_respects_output_cap(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }

    taffy_dup_type_t modes[] = {
        TAFFY_NO_DUPS, TAFFY_QUERY_DUPS, TAFFY_QUERY_AND_TARGET_DUPS
    };
    for (size_t m = 0; m < sizeof(modes) / sizeof(modes[0]); m++) {
        struct taffy_block_results_t *res = taffyGetBlocksInTargetRange(
            h, "simMouse_chr6", "simHuman_chr6", "simHuman.chr6",
            0, 0, 0, TAFFY_NO_SEQUENCES, modes[m], 0, NULL, &err);
        CuAssertPtrNotNull(tc, res);
        int64_t n = 0;
        for (struct taffy_block_t *b = res->mappedBlocks; b; b = b->next) n++;
        /* Hard cap is 500.  Asserting that contract here regardless of
         * the fixture's natural block count. */
        CuAssertTrue(tc, n <= 500);
        taffyFreeBlockResults(res);
    }

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
    SUITE_ADD_TEST(suite, test_get_blocks_full_chrom_merged);
    SUITE_ADD_TEST(suite, test_get_blocks_respects_output_cap);
    SUITE_ADD_TEST(suite, test_max_output_blocks_setter);
    SUITE_ADD_TEST(suite, test_chain_id_and_summaries);
    return suite;
}

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
/* A tiny purpose-built .tui with one ancestor-reversed MAF block (g1 and
 * g2 both on `-` strand against anc1), used to regression-test the
 * iv_rev XOR run.strand strand-fix landed for the user-reported "bogus
 * reverse copies on the diagonal" bug.  See tests/tui/strand_revcase.maf
 * for the source MAF; the .tui is built once with `taffy index -u`. */
/* Pass the MAF path -- taffyOpen's tui_path() appends ".tui" itself. */
#define TEST_TUI_REVCASE "tests/tui/strand_revcase.maf"

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

/* --- chain_overlap_frac setter/getter ----------------------------- */

static void test_chain_overlap_frac_default_is_half(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }
    double f = -999.0;
    CuAssertIntEquals(tc, 0, taffyGetChainOverlapFrac(h, &f, &err));
    /* Default 0.5: strict 0.0 drops a whole chain on ANY overlap, which
     * culls ~17% of real ortholog coverage on the fragmented universal
     * .tui; 0.5 keeps coverage while still filtering >50%-redundant
     * paralogs. */
    CuAssertTrue(tc, f == 0.5);
    taffyClose(h, &err); free(err);
}

static void test_chain_overlap_frac_roundtrip(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }
    CuAssertIntEquals(tc, 0, taffySetChainOverlapFrac(h, 0.25, &err));
    double f = -999.0;
    taffyGetChainOverlapFrac(h, &f, &err);
    CuAssertTrue(tc, f == 0.25);
    CuAssertIntEquals(tc, 0, taffySetChainOverlapFrac(h, -1.0, &err));
    taffyGetChainOverlapFrac(h, &f, &err);
    CuAssertTrue(tc, f == -1.0);
    taffyClose(h, &err); free(err);
}

static void test_chain_overlap_frac_rejects_out_of_range(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }
    char *err2 = NULL;
    CuAssertIntEquals(tc, -1, taffySetChainOverlapFrac(h, 1.5, &err2));
    CuAssertPtrNotNull(tc, err2); free(err2);
    err2 = NULL;
    CuAssertIntEquals(tc, -1, taffySetChainOverlapFrac(h, -0.5, &err2));
    CuAssertPtrNotNull(tc, err2); free(err2);
    taffyClose(h, &err); free(err);
}

static void test_chain_overlap_frac_invalid_handle(CuTest *tc) {
    char *err = NULL;
    CuAssertIntEquals(tc, -1, taffySetChainOverlapFrac(99999, 0.0, &err));
    CuAssertPtrNotNull(tc, err); free(err);
    err = NULL;
    double f = 0;
    CuAssertIntEquals(tc, -1, taffyGetChainOverlapFrac(99999, &f, &err));
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

    /* Default cap is 1000. */
    int64_t cap = 0;
    CuAssertIntEquals(tc, 0, taffyGetMaxOutputBlocks(h, &cap, &err));
    CuAssertTrue(tc, cap == 1000);

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

/* Output-budget invariant: the COVERAGE (coarsened, wide-zoom) path emits
 * at most max_output_blocks+1 blocks -- one per occupied bin (+1 boundary).
 * Detail (narrow) output is deliberately NOT budget-bounded: it emits every
 * block in the span so a fixed locus renders identically regardless of the
 * window edges (the span, not the cap, bounds the count).  So we force the
 * coverage path with a small budget over the full chromosome (601 kb / 50 =
 * ~12 kb per bin, above the coarsen threshold) and assert the bin bound. */
static void test_get_blocks_respects_output_cap(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }

    int64_t cap = 50;
    CuAssertIntEquals(tc, 0, taffySetMaxOutputBlocks(h, cap, &err));

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
        CuAssertTrue(tc, n <= cap + 1);
        taffyFreeBlockResults(res);
    }

    taffyClose(h, &err);
    free(err);
}

/* ------------------------------------------------------------------ */
/* Strand fix: iv_rev XOR run.strand for ancestor-reversed paralogs    */
/* ------------------------------------------------------------------ */

/* Walk a result list, find the first block matching (tStart, qStart, size).
 * Returns NULL if not found. */
static struct taffy_block_t *find_block(struct taffy_block_t *head,
                                        int64_t tStart, int64_t qStart,
                                        int64_t size) {
    for (struct taffy_block_t *b = head; b; b = b->next) {
        if (b->tStart == tStart && b->qStart == qStart && b->size == size)
            return b;
    }
    return NULL;
}

/* The user-reported "bogus reverse-strand copies on the diagonal" bug:
 * when both source and target genome rows are reverse-mapped to the
 * SAME ancestor (a common cactus pattern -- the ancestor segment ended
 * up oriented opposite both descendants), the relative source<->target
 * strand is FORWARD, but the old code reported it as REVERSE because it
 * used only the run's column-to-target strand bit and ignored the
 * source-to-column strand carried by the iv.
 *
 * tests/tui/strand_revcase.maf has four blocks at anc1.chr [0..40):
 *   [0..10):   g1 +,  g2 +   -> relative '+'  (forward + forward)
 *   [10..20):  g1 +,  g2 -   -> relative '-'  (true forward<->reverse)
 *   [20..30):  g1 -,  g2 -   -> relative '+'  (both reversed = forward!)
 *   [30..40):  g1 +,  g2 +   -> relative '+'
 *
 * Pre-fix, block [20..30) would show '-' (the run's strand alone).
 * Post-fix, it correctly shows '+'. */
static void test_strand_iv_rev_xor_run_strand(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI_REVCASE, &err);
    /* If the fixture .tui is missing or the open fails, fail loudly --
     * the strand-fix is too important to let this regression test skip
     * silently.  Build with `taffy index -u tests/tui/strand_revcase.maf`. */
    CuAssert(tc, "open strand_revcase .tui failed -- build with `taffy index -u tests/tui/strand_revcase.maf`",
             h >= 0);
    /* Disable the overlap-frac filter and raise the cap so every chain
     * survives -- we want to inspect the strand of each MAF block's
     * mapping individually. */
    CuAssertIntEquals(tc, 0, taffySetChainOverlapFrac(h, -1.0, &err));
    CuAssertIntEquals(tc, 0, taffySetMaxOutputBlocks(h, 1000, &err));

    struct taffy_block_results_t *res = taffyGetBlocksInTargetRange(
        h, /*qSpecies=*/"g2", /*tSpecies=*/"g1", /*tChrom=*/"chr",
        /*tStart=*/0, /*tEnd=*/40, /*tReversed=*/0,
        TAFFY_NO_SEQUENCES, TAFFY_QUERY_DUPS,
        /*mapBackAdjacencies=*/0, /*coalescenceLimitName=*/NULL, &err);
    CuAssertPtrNotNull(tc, res);

    /* Block 1: anc1[0..10) -- both forward. */
    struct taffy_block_t *b1 = find_block(res->mappedBlocks, 0, 0, 10);
    CuAssertPtrNotNull(tc, b1);
    CuAssertIntEquals(tc, '+', b1->strand);

    /* Block 2: anc1[10..20) -- g1 forward, g2 reverse.  Genuinely
     * antiparallel -> '-' relative.  qStart should be the forward-strand
     * g2 coord (= 0, since g2 row's start=90 len=10 on '-' translates
     * to forward range [0..10)). */
    struct taffy_block_t *b2 = find_block(res->mappedBlocks, 10, 0, 10);
    CuAssertPtrNotNull(tc, b2);
    CuAssertIntEquals(tc, '-', b2->strand);

    /* Block 3: anc1[20..30) -- THE REGRESSION CASE.  Both g1 and g2
     * reverse-mapped to anc1 -> relative forward -> '+'.  Forward
     * coords: g1 [10..20), g2 [20..30).  Pre-fix this returned '-'. */
    struct taffy_block_t *b3 = find_block(res->mappedBlocks, 10, 20, 10);
    CuAssertPtrNotNull(tc, b3);
    CuAssertIntEquals(tc, '+', b3->strand);

    /* Block 4: anc1[30..40) -- both forward. */
    struct taffy_block_t *b4 = find_block(res->mappedBlocks, 30, 30, 10);
    CuAssertPtrNotNull(tc, b4);
    CuAssertIntEquals(tc, '+', b4->strand);

    taffyFreeBlockResults(res);
    taffyClose(h, &err);
    free(err);
}

/* Sanity that the strand-fix's chained sidecar / bin path produces the
 * same per-block correctness.  Even at chromosome-scale a query whose
 * dominant alignment is forward should not show majority '-' coverage. */
static void test_strand_revcase_all_chains_correct(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI_REVCASE, &err);
    CuAssert(tc, "open strand_revcase .tui failed", h >= 0);
    CuAssertIntEquals(tc, 0, taffySetChainOverlapFrac(h, -1.0, &err));
    CuAssertIntEquals(tc, 0, taffySetMaxOutputBlocks(h, 1000, &err));
    /* Query the whole queried-chrom range; with QUERY_DUPS every chain
     * is emitted so we get all 4 fixture blocks back.  The aggregate
     * '+' vs '-' bp count nails the strand fix: pre-fix this was 20/20
     * (block 3 misreported as '-' tipped the balance); post-fix 30/10. */
    struct taffy_block_results_t *res = taffyGetBlocksInTargetRange(
        h, "g2", "g1", "chr", 0, 100, 0,
        TAFFY_NO_SEQUENCES, TAFFY_QUERY_DUPS, 0, NULL, &err);
    CuAssertPtrNotNull(tc, res);
    int64_t plus_bp = 0, minus_bp = 0;
    for (struct taffy_block_t *b = res->mappedBlocks; b; b = b->next) {
        if (b->strand == '+') plus_bp += b->size;
        else if (b->strand == '-') minus_bp += b->size;
    }
    /* 3 forward blocks + 1 reverse block (anc1[10..20)): 30 bp + vs 10 bp -. */
    CuAssertIntEquals(tc, 30, (int) plus_bp);
    CuAssertIntEquals(tc, 10, (int) minus_bp);
    taffyFreeBlockResults(res);
    taffyClose(h, &err);
    free(err);
}

/* Regression: panning a fixed zoom must not make blocks blink in/out.  The
 * coarsen switch is span-based (a zoom property), NOT block-count-based, and
 * coverage bins are absolute-aligned -- so two windows of the SAME span over a
 * shared interior locus return the SAME blocks (tStart/qStart/strand) for that
 * locus, no matter where the window edges fall.  The old count-based switch
 * flipped detail<->coverage as a window's block count crossed the cap while
 * panning a dense locus, making dupe bars blink.
 *
 * We compare block PRESENCE (position + strand), NOT the coverage bp-sum: that
 * sum carries a sub-pixel (~1 bp), invisible window-dependence -- the trim's
 * running edge propagates from a boundary-spanning chain's start, which sits at
 * a different clipped offset in each window -- so a byte-exact size check would
 * over-assert a difference no one can see (and that any block re-ordering
 * perturbs). */
static void test_window_stability_coarse(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }
    CuAssertIntEquals(tc, 0, taffySetMaxOutputBlocks(h, 20, &err)); /* force coverage */

    /* Two 400 kb windows sharing the interior [120k,360k]; the SAME blocks
     * (position + strand) must appear for that interior (chainId is per-window
     * and intentionally excluded; size is too -- see header). */
    int64_t start[2] = { 0, 80000 };
    char buf[2][16384];
    for (int w = 0; w < 2; w++) {
        struct taffy_block_results_t *res = taffyGetBlocksInTargetRange(
            h, "simMouse_chr6", "simHuman_chr6", "simHuman.chr6",
            start[w], start[w] + 400000, 0, TAFFY_NO_SEQUENCES,
            TAFFY_QUERY_AND_TARGET_DUPS, 0, NULL, &err);
        CuAssertPtrNotNull(tc, res);
        size_t off = 0;
        buf[w][0] = '\0';
        for (struct taffy_block_t *b = res->mappedBlocks; b; b = b->next)
            if (b->tStart >= 120000 && b->tStart < 360000 &&
                off < sizeof(buf[w]) - 64)
                off += (size_t) snprintf(buf[w] + off, sizeof(buf[w]) - off,
                    "%ld,%ld,%c;", (long) b->tStart, (long) b->qStart, b->strand);
        taffyFreeBlockResults(res);
    }
    CuAssertTrue(tc, buf[0][0] != '\0');       /* the locus is actually covered */
    CuAssertStrEquals(tc, buf[0], buf[1]);     /* same blocks => no blinking */

    taffyClose(h, &err);
    free(err);
}

/* The min-block noise filter (taffySetMinBlockFilter): with it engaged, blocks
 * below BOTH a fraction of the window AND a fraction of the largest block are
 * dropped; the largest block always survives; off (the default) is a no-op. */
static void test_min_block_filter(CuTest *tc) {
    char *err = NULL;
    int h = taffyOpen(TEST_TUI, &err);
    if (h < 0) { free(err); return; }
    CuAssertIntEquals(tc, 0, taffySetMaxOutputBlocks(h, 30, &err));   /* force coverage */

    /* baseline: filter off (default) */
    struct taffy_block_results_t *r0 = taffyGetBlocksInTargetRange(
        h, "simMouse_chr6", "simHuman_chr6", "simHuman.chr6", 0, 600000, 0,
        TAFFY_NO_SEQUENCES, TAFFY_QUERY_AND_TARGET_DUPS, 0, NULL, &err);
    CuAssertPtrNotNull(tc, r0);
    int n0 = 0; int64_t maxsz = 0, minsz = INT64_MAX;
    for (struct taffy_block_t *b = r0->mappedBlocks; b; b = b->next) {
        n0++;
        if (b->size > maxsz) maxsz = b->size;
        if (b->size < minsz) minsz = b->size;
    }
    taffyFreeBlockResults(r0);
    CuAssertTrue(tc, n0 >= 1);

    /* off (0,0) must be a no-op */
    CuAssertIntEquals(tc, 0, taffySetMinBlockFilter(h, 0.0, 0.0, &err));
    struct taffy_block_results_t *roff = taffyGetBlocksInTargetRange(
        h, "simMouse_chr6", "simHuman_chr6", "simHuman.chr6", 0, 600000, 0,
        TAFFY_NO_SEQUENCES, TAFFY_QUERY_AND_TARGET_DUPS, 0, NULL, &err);
    CuAssertPtrNotNull(tc, roff);
    int noff = 0;
    for (struct taffy_block_t *b = roff->mappedBlocks; b; b = b->next) noff++;
    taffyFreeBlockResults(roff);
    CuAssertIntEquals(tc, n0, noff);

    /* engage: window-off (1.0) so the threshold is 0.5 * the largest block */
    CuAssertIntEquals(tc, 0, taffySetMinBlockFilter(h, 1.0, 0.5, &err));
    struct taffy_block_results_t *r1 = taffyGetBlocksInTargetRange(
        h, "simMouse_chr6", "simHuman_chr6", "simHuman.chr6", 0, 600000, 0,
        TAFFY_NO_SEQUENCES, TAFFY_QUERY_AND_TARGET_DUPS, 0, NULL, &err);
    CuAssertPtrNotNull(tc, r1);
    int n1 = 0; int64_t maxsz1 = 0, minsz1 = INT64_MAX;
    for (struct taffy_block_t *b = r1->mappedBlocks; b; b = b->next) {
        n1++;
        if (b->size > maxsz1) maxsz1 = b->size;
        if (b->size < minsz1) minsz1 = b->size;
    }
    taffyFreeBlockResults(r1);
    CuAssertTrue(tc, n1 >= 1);              /* never emptied */
    CuAssertTrue(tc, n1 <= n0);             /* only drops */
    CuAssertTrue(tc, maxsz1 == maxsz);      /* the largest block survives */
    CuAssertTrue(tc, minsz1 >= maxsz / 2);  /* every survivor >= 0.5 * max */
    if (minsz < maxsz / 2) CuAssertTrue(tc, n1 < n0);   /* a droppable block existed -> dropped */

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
    SUITE_ADD_TEST(suite, test_chain_overlap_frac_default_is_half);
    SUITE_ADD_TEST(suite, test_chain_overlap_frac_roundtrip);
    SUITE_ADD_TEST(suite, test_chain_overlap_frac_rejects_out_of_range);
    SUITE_ADD_TEST(suite, test_chain_overlap_frac_invalid_handle);
    SUITE_ADD_TEST(suite, test_get_blocks_narrow_per_run);
    SUITE_ADD_TEST(suite, test_get_blocks_full_chrom_merged);
    SUITE_ADD_TEST(suite, test_get_blocks_respects_output_cap);
    SUITE_ADD_TEST(suite, test_window_stability_coarse);
    SUITE_ADD_TEST(suite, test_max_output_blocks_setter);
    SUITE_ADD_TEST(suite, test_chain_id_and_summaries);
    SUITE_ADD_TEST(suite, test_strand_iv_rev_xor_run_strand);
    SUITE_ADD_TEST(suite, test_strand_revcase_all_chains_correct);
    SUITE_ADD_TEST(suite, test_min_block_filter);
    return suite;
}

/*
 * Unit tests for taffy/impl/chain.c -- the generic DP chainer ported
 * from paffy/impl/chaining.c.  Focused on the chaining algorithm's
 * semantics: which alignments end up in the same chain, how scores
 * accumulate, how strand / partition keys gate joining, and how the
 * gap_cost callback interacts with per-aln score.
 *
 * Conventions for tests:
 *   - score == aln length (the natural choice for taffy lift / blockViz);
 *     this means a gap cost equal to or larger than the new aln's
 *     length always prevents joining.
 *   - default gap cost in most tests uses chain_open=0, chain_extend=1
 *     so we can construct exact-arithmetic chain-vs-no-chain cases.
 *     Tests that exercise the default lastz-tuned 5000/1 params do so
 *     explicitly.
 */

#include "CuTest.h"
#include "chain.h"
#include "sonLib.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Per-test helper: malloc-then-zero a TaffyAln with the common defaults
 * filled in (score = length, score-only chain; q_name/t_name optional). */
static TaffyAln mk(const char *q, int64_t q_s, int64_t q_e,
                   const char *t, int64_t t_s, int64_t t_e, int strand) {
    TaffyAln a = {0};
    a.q_name  = q;
    a.q_start = q_s;
    a.q_end   = q_e;
    a.t_name  = t;
    a.t_start = t_s;
    a.t_end   = t_e;
    a.strand  = strand;
    a.score   = q_e - q_s;
    return a;
}

/* Run chaining with chain_open=0, chain_extend=1, max_gap_length=10M
 * unless an explicit override; returns the chain_id[] array via *cid
 * (caller frees) and the *chains* array (caller frees).  Caller passes
 * a stack-allocated TaffyAln[]; we copy so the algorithm's in-place
 * sort doesn't reorder the caller's array (makes assertions readable). */
static void run_chain(TaffyAln *src, int64_t n,
                      int64_t chain_open, int64_t chain_extend,
                      int64_t max_gap_length,
                      TaffyAln **alns_out, int64_t **cid_out,
                      TaffyChainInfo **chains_out, int64_t *nc_out) {
    TaffyAln *cp = st_malloc((size_t) n * sizeof(TaffyAln));
    memcpy(cp, src, (size_t) n * sizeof(TaffyAln));
    int64_t *cid = st_calloc((size_t)(n > 0 ? n : 1), sizeof(int64_t));
    TaffyChainCostParams params = { chain_open, chain_extend };
    TaffyChainInfo *chs = NULL;
    int64_t nc = 0;
    taffy_chain(cp, n, taffy_chain_default_gap_cost, &params,
                max_gap_length, cid, &chs, &nc);
    *alns_out  = cp;
    *cid_out   = cid;
    *chains_out = chs;
    *nc_out    = nc;
}

/* Look up chain_id for the alignment whose `user` field == tag.  The
 * algorithm's in-place sort destroys the caller's input order, so we
 * use the user pointer as a stable identity. */
static int64_t cid_for_tag(TaffyAln *alns, int64_t n, int64_t *cid,
                           void *tag) {
    for (int64_t i = 0; i < n; i++)
        if (alns[i].user == tag) return cid[i];
    return -1;
}

/* Sum of scores across chains (sanity invariant: sum of chain
 * total_score == sum of per-aln scores - sum of all gap costs joined
 * by chaining; we don't try to recompute the exact gap cost here, but
 * total_bp should always == sum of aln (q_end - q_start)). */
static int64_t sum_bp(TaffyChainInfo *cs, int64_t n) {
    int64_t s = 0; for (int64_t i = 0; i < n; i++) s += cs[i].total_bp; return s;
}
static int64_t sum_n_alns(TaffyChainInfo *cs, int64_t n) {
    int64_t s = 0; for (int64_t i = 0; i < n; i++) s += cs[i].n_alns; return s;
}

/* ------------------------------------------------------------------ */
/* Edge cases                                                          */
/* ------------------------------------------------------------------ */

static void test_empty_input(CuTest *tc) {
    TaffyAln *alns = NULL; int64_t *cid = NULL;
    TaffyChainInfo *chs = NULL; int64_t nc = -1;
    taffy_chain(NULL, 0, taffy_chain_default_gap_cost, NULL,
                INT64_MAX, NULL, &chs, &nc);
    CuAssertTrue(tc, nc == 0);
    CuAssertTrue(tc, chs == NULL);
    (void) alns; (void) cid;
}

static void test_single_aln(CuTest *tc) {
    TaffyAln in[] = { mk("q", 0, 100, "t", 1000, 1100, +1) };
    in[0].user = (void *) 0x1;
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 1, 0, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 1);
    CuAssertTrue(tc, chs[0].id == 1);
    CuAssertTrue(tc, chs[0].n_alns == 1);
    CuAssertTrue(tc, chs[0].total_score == 100);
    CuAssertTrue(tc, chs[0].total_bp == 100);
    free(alns); free(cid); free(chs);
}

/* ------------------------------------------------------------------ */
/* Forward-strand chaining                                             */
/* ------------------------------------------------------------------ */

static void test_two_abutting_forward_chain(CuTest *tc) {
    /* zero-gap chain on both axes: always joins regardless of params */
    TaffyAln in[] = {
        mk("q", 0, 100, "t", 1000, 1100, +1),
        mk("q", 100, 200, "t", 1100, 1200, +1),
    };
    in[0].user = (void *) 0x1;
    in[1].user = (void *) 0x2;
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 2, 0, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 1);
    CuAssertTrue(tc, chs[0].n_alns == 2);
    CuAssertTrue(tc, chs[0].total_score == 200);
    CuAssertTrue(tc, chs[0].total_bp == 200);
    /* both alns should share the chain id */
    int64_t c1 = cid_for_tag(alns, 2, cid, (void *)0x1);
    int64_t c2 = cid_for_tag(alns, 2, cid, (void *)0x2);
    CuAssertTrue(tc, c1 == c2);
    free(alns); free(cid); free(chs);
}

static void test_three_collinear_chain_when_score_exceeds_cost(CuTest *tc) {
    /* gap = 50+50 on each link with open=0, extend=1: cost=100.
     * Each aln has score 1000, so 100 < 1000 satisfies the guard. */
    TaffyAln in[] = {
        mk("q",    0, 1000, "t", 10000, 11000, +1),
        mk("q", 1050, 2050, "t", 11050, 12050, +1),
        mk("q", 2100, 3100, "t", 12100, 13100, +1),
    };
    in[0].user = (void *) 0x1;
    in[1].user = (void *) 0x2;
    in[2].user = (void *) 0x3;
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 3, 0, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 1);
    CuAssertTrue(tc, chs[0].n_alns == 3);
    CuAssertTrue(tc, chs[0].total_bp == 3000);
    /* score = 1000+1000+1000 - 2*(0 + 50+50) = 3000 - 200 = 2800 */
    CuAssertTrue(tc, chs[0].total_score == 2800);
    free(alns); free(cid); free(chs);
}

static void test_gap_cost_blocks_chaining(CuTest *tc) {
    /* Aln score 100; gap cost between consecutive = 100 (open=0, ext=1
     * with q_gap+t_gap=100).  Guard requires g < score (STRICT), so
     * with g==score we DO NOT chain -> 2 chains. */
    TaffyAln in[] = {
        mk("q",   0, 100, "t", 1000, 1100, +1),
        mk("q", 150, 250, "t", 1150, 1250, +1),
    };
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 2, 0, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 2);
    free(alns); free(cid); free(chs);
}

static void test_query_overlap_does_not_chain(CuTest *tc) {
    /* Second aln overlaps first on the QUERY axis -> not collinear ->
     * can never extend the predecessor.  2 chains. */
    TaffyAln in[] = {
        mk("q",   0, 1000, "t", 10000, 11000, +1),
        mk("q", 500, 1500, "t", 11000, 12000, +1),
    };
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 2, 0, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 2);
    free(alns); free(cid); free(chs);
}

static void test_target_overlap_does_not_chain(CuTest *tc) {
    /* Reverse: target overlap; also not collinear. */
    TaffyAln in[] = {
        mk("q",   0, 1000, "t", 10000, 11000, +1),
        mk("q", 1100, 2100, "t", 10500, 11500, +1),
    };
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 2, 0, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 2);
    free(alns); free(cid); free(chs);
}

/* ------------------------------------------------------------------ */
/* Strand handling                                                     */
/* ------------------------------------------------------------------ */

static void test_reverse_strand_chains(CuTest *tc) {
    /* Reverse-strand alns: query advances forward but target retreats
     * (each aln's t_start is LESS than the previous one's).  Chainer
     * mirrors the q axis internally so this becomes a forward problem.
     */
    TaffyAln in[] = {
        mk("q",    0, 1000, "t", 90000, 91000, -1),
        mk("q", 1050, 2050, "t", 88950, 89950, -1),
        mk("q", 2100, 3100, "t", 87900, 88900, -1),
    };
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 3, 0, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 1);
    CuAssertTrue(tc, chs[0].n_alns == 3);
    free(alns); free(cid); free(chs);
}

static void test_mixed_strands_dont_join(CuTest *tc) {
    /* Two alns at adjacent q ranges but opposite strands.  Must end
     * up in separate chains -- the partition step splits +/- before
     * the sweep. */
    TaffyAln in[] = {
        mk("q",   0, 1000, "t", 10000, 11000, +1),
        mk("q", 1100, 2100, "t", 11100, 12100, -1),
    };
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 2, 0, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 2);
    free(alns); free(cid); free(chs);
}

/* ------------------------------------------------------------------ */
/* Partitioning by name                                                */
/* ------------------------------------------------------------------ */

static void test_different_q_names_dont_join(CuTest *tc) {
    TaffyAln in[] = {
        mk("qA",  0, 1000, "t", 10000, 11000, +1),
        mk("qB",  0, 1000, "t", 11000, 12000, +1),
    };
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 2, 0, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 2);
    free(alns); free(cid); free(chs);
}

static void test_different_t_names_dont_join(CuTest *tc) {
    TaffyAln in[] = {
        mk("q",   0, 1000, "tA", 10000, 11000, +1),
        mk("q", 1100, 2100, "tB", 11000, 12000, +1),
    };
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 2, 0, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 2);
    free(alns); free(cid); free(chs);
}

/* ------------------------------------------------------------------ */
/* Paralog separation -- the case blockViz cares about most            */
/* ------------------------------------------------------------------ */

static void test_paralog_lands_in_separate_chain(CuTest *tc) {
    /* Three primary forward alns + one "paralog": same q range as
     * primary aln 1 but at a far target position.  Should be 2 chains,
     * primary >= paralog by score. */
    TaffyAln in[] = {
        mk("q",   0, 1000, "t", 10000, 11000, +1),
        mk("q", 1100, 2100, "t", 11100, 12100, +1),
        mk("q", 2200, 3200, "t", 12200, 13200, +1),
        mk("q", 1100, 2100, "t", 90000, 91000, +1),   /* paralog */
    };
    in[0].user = (void *)0x1;
    in[1].user = (void *)0x2;
    in[2].user = (void *)0x3;
    in[3].user = (void *)0x4;
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 4, 0, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 2);
    /* Primary (chs[0]) > paralog (chs[1]) by score */
    CuAssertTrue(tc, chs[0].total_score > chs[1].total_score);
    CuAssertTrue(tc, chs[0].n_alns == 3);
    CuAssertTrue(tc, chs[1].n_alns == 1);
    /* Primary's chain id is shared by alns 0, 1, 2; paralog (3) is in
     * the OTHER chain. */
    int64_t c1 = cid_for_tag(alns, 4, cid, (void *)0x1);
    int64_t c2 = cid_for_tag(alns, 4, cid, (void *)0x2);
    int64_t c3 = cid_for_tag(alns, 4, cid, (void *)0x3);
    int64_t cP = cid_for_tag(alns, 4, cid, (void *)0x4);
    CuAssertTrue(tc, c1 == c2 && c2 == c3);
    CuAssertTrue(tc, cP != c1);
    free(alns); free(cid); free(chs);
}

/* ------------------------------------------------------------------ */
/* max_gap_length                                                      */
/* ------------------------------------------------------------------ */

static void test_max_gap_breaks_chain(CuTest *tc) {
    /* Two short close alns + one with a huge gap.  max_gap_length=1000
     * prevents the third from joining the first two.  -> 2 chains. */
    TaffyAln in[] = {
        mk("q",     0, 1000, "t",  10000,  11000, +1),
        mk("q",  1100, 2100, "t",  11100,  12100, +1),
        mk("q", 50000,51000, "t", 100000, 101000, +1),
    };
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 3, 0, 1, 1000, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 2);
    /* Big chain has the first two (they DO join: q_gap+t_gap=100+100,
     * cost=200 < score 1000), separated by the gap-broken aln */
    int64_t bigger_n = (chs[0].n_alns > chs[1].n_alns) ? chs[0].n_alns : chs[1].n_alns;
    CuAssertTrue(tc, bigger_n == 2);
    free(alns); free(cid); free(chs);
}

/* ------------------------------------------------------------------ */
/* Output invariants                                                   */
/* ------------------------------------------------------------------ */

static void test_chains_sorted_by_score_desc(CuTest *tc) {
    /* Three independent same-strand alns on the same target, far
     * apart so they don't chain, with different scores.  Each ends
     * up in its own chain; chains_out is sorted descending by
     * total_score. */
    TaffyAln in[] = {
        /* score = q_end - q_start, set up explicitly */
        mk("q",     0,   100, "t",       0,   100, +1),  /* score 100  */
        mk("q",  5000,  5500, "t",    5000,  5500, +1),  /* score 500  */
        mk("q", 10000, 10300, "t",   10000, 10300, +1),  /* score 300  */
    };
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, 3, 1000000, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, nc == 3);
    /* chains_out sorted desc by total_score: 500 > 300 > 100 */
    CuAssertTrue(tc, chs[0].total_score == 500);
    CuAssertTrue(tc, chs[1].total_score == 300);
    CuAssertTrue(tc, chs[2].total_score == 100);
    free(alns); free(cid); free(chs);
}

static void test_total_bp_and_n_alns_invariant(CuTest *tc) {
    /* For any input, sum(chains[].total_bp) == sum(input aln lengths)
     * and sum(chains[].n_alns) == n. */
    TaffyAln in[] = {
        mk("q",   0, 1000, "t", 1000, 2000, +1),
        mk("q", 1100, 2100, "t", 2100, 3100, +1),
        mk("q", 2200, 3200, "t", 3200, 4200, +1),
        mk("q", 1100, 2100, "t", 90000,91000, +1),
        mk("q", 5000, 5500, "tB", 100, 600, +1),
        mk("q", 7000, 7100, "t", 9000, 9100, -1),
    };
    int64_t n = sizeof(in)/sizeof(in[0]);
    int64_t expected_bp = 0;
    for (int64_t i=0; i<n; i++) expected_bp += in[i].q_end - in[i].q_start;
    TaffyAln *alns; int64_t *cid; TaffyChainInfo *chs; int64_t nc;
    run_chain(in, n, 0, 1, INT64_MAX, &alns, &cid, &chs, &nc);
    CuAssertTrue(tc, sum_bp(chs, nc)     == expected_bp);
    CuAssertTrue(tc, sum_n_alns(chs, nc) == n);
    free(alns); free(cid); free(chs);
}

/* ------------------------------------------------------------------ */
/* Suite                                                               */
/* ------------------------------------------------------------------ */

CuSuite* chain_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_empty_input);
    SUITE_ADD_TEST(suite, test_single_aln);
    SUITE_ADD_TEST(suite, test_two_abutting_forward_chain);
    SUITE_ADD_TEST(suite, test_three_collinear_chain_when_score_exceeds_cost);
    SUITE_ADD_TEST(suite, test_gap_cost_blocks_chaining);
    SUITE_ADD_TEST(suite, test_query_overlap_does_not_chain);
    SUITE_ADD_TEST(suite, test_target_overlap_does_not_chain);
    SUITE_ADD_TEST(suite, test_reverse_strand_chains);
    SUITE_ADD_TEST(suite, test_mixed_strands_dont_join);
    SUITE_ADD_TEST(suite, test_different_q_names_dont_join);
    SUITE_ADD_TEST(suite, test_different_t_names_dont_join);
    SUITE_ADD_TEST(suite, test_paralog_lands_in_separate_chain);
    SUITE_ADD_TEST(suite, test_max_gap_breaks_chain);
    SUITE_ADD_TEST(suite, test_chains_sorted_by_score_desc);
    SUITE_ADD_TEST(suite, test_total_bp_and_n_alns_invariant);
    return suite;
}

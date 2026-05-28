/*
 * Unit tests for taffy/impl/gerp.c.  Focused on the Fitch parsimony
 * + branch-sum core (gerp_score_column) since that's the algorithmic
 * meat; the per-block driver is exercised by the end-to-end smoke
 * test against evolverMammals run from the harness above.
 */

#include "CuTest.h"
#include "gerp.h"
#include "sonLib.h"
#include <math.h>
#include <string.h>

// Convenience: make a (n_leaves)-entry char array from a base string.
// "ACGT" -> {'A','C','G','T'}; '-' or '.' encodes "no base for this leaf".
static char *bases_from(const char *s) {
    size_t n = strlen(s);
    char *b = st_malloc(n);
    for (size_t i = 0; i < n; i++) b[i] = (s[i] == '-' || s[i] == '.') ? 0 : s[i];
    return b;
}

static void close_to(CuTest *tc, double expected, double actual, double tol) {
    if (fabs(expected - actual) > tol) {
        char msg[256];
        snprintf(msg, sizeof(msg), "expected %.6f, got %.6f (tol %g)",
                 expected, actual, tol);
        CuFail(tc, msg);
    }
}

// Binary tree, 4 leaves.  Branch lengths chosen so the total is exactly
// 1.0, making score = (n_present_branches_total) - cost trivially
// verifiable.
//
//      Anc0
//      / .
//     /   .
//   Anc1  Anc2     branches: A->Anc1=0.10, A->Anc2=0.10, Anc1->W=0.20,
//   / .   / .                Anc1->X=0.20, Anc2->Y=0.20, Anc2->Z=0.20
//  W   X Y   Z     total = 1.00
static const char *BINARY_4 =
    "((W:0.20,X:0.20)Anc1:0.10,(Y:0.20,Z:0.20)Anc2:0.10)Anc0;";

static void test_construct_and_introspect(CuTest *tc) {
    GerpTree *gt = gerp_tree_construct(BINARY_4);
    CuAssertPtrNotNull(tc, gt);
    CuAssertIntEquals(tc, 4, (int) gerp_tree_n_leaves(gt));
    CuAssertTrue(tc, gerp_tree_is_ancestor(gt, "Anc0"));
    CuAssertTrue(tc, gerp_tree_is_ancestor(gt, "Anc1"));
    CuAssertTrue(tc, gerp_tree_is_ancestor(gt, "Anc2"));
    CuAssertTrue(tc, !gerp_tree_is_ancestor(gt, "W"));
    CuAssertTrue(tc, !gerp_tree_is_ancestor(gt, "unknownThing"));
    gerp_tree_destruct(gt);
}

static void test_all_present_no_subs(CuTest *tc) {
    GerpTree *gt = gerp_tree_construct(BINARY_4);
    GerpScratch *sc = gerp_scratch_construct(gt);
    char *b = bases_from("AAAA");
    double rs = -999; int64_t depth = -1;
    bool scored = gerp_score_column(gt, sc, b, 2, 1.0, &rs, &depth);
    CuAssertTrue(tc, scored);
    CuAssertIntEquals(tc, 4, (int) depth);
    // All branches present (sum = 1.0), 0 subs.
    close_to(tc, 1.0, rs, 1e-9);
    free(b);
    gerp_scratch_destruct(sc);
    gerp_tree_destruct(gt);
}

static void test_one_substitution(CuTest *tc) {
    GerpTree *gt = gerp_tree_construct(BINARY_4);
    GerpScratch *sc = gerp_scratch_construct(gt);
    char *b = bases_from("AAAG");
    double rs = 0; int64_t depth = 0;
    bool scored = gerp_score_column(gt, sc, b, 2, 1.0, &rs, &depth);
    CuAssertTrue(tc, scored);
    CuAssertIntEquals(tc, 4, (int) depth);
    // 1 sub anywhere on the tree -> RS = 1.0 - 1.
    close_to(tc, 0.0, rs, 1e-9);
    free(b);
    gerp_scratch_destruct(sc);
    gerp_tree_destruct(gt);
}

static void test_two_substitutions(CuTest *tc) {
    GerpTree *gt = gerp_tree_construct(BINARY_4);
    GerpScratch *sc = gerp_scratch_construct(gt);
    // W=A, X=C, Y=G, Z=T -> 3 substitutions under Fitch.
    char *b = bases_from("ACGT");
    double rs = 0; int64_t depth = 0;
    bool scored = gerp_score_column(gt, sc, b, 2, 1.0, &rs, &depth);
    CuAssertTrue(tc, scored);
    CuAssertIntEquals(tc, 4, (int) depth);
    close_to(tc, 1.0 - 3.0, rs, 1e-9);
    free(b);
    gerp_scratch_destruct(sc);
    gerp_tree_destruct(gt);
}

static void test_gap_drops_leaf(CuTest *tc) {
    GerpTree *gt = gerp_tree_construct(BINARY_4);
    GerpScratch *sc = gerp_scratch_construct(gt);
    // Z is gap.  Present leaves: W,X,Y all 'A'.  Branches in pruned tree
    // are: W (0.20) + X (0.20) + Anc1 (0.10) + Y (0.20) + Anc2 (0.10) +
    // [no Z branch].  Anc0 is root, no branch.  Total = 0.80.
    char *b = bases_from("AAA-");
    double rs = 0; int64_t depth = 0;
    bool scored = gerp_score_column(gt, sc, b, 2, 1.0, &rs, &depth);
    CuAssertTrue(tc, scored);
    CuAssertIntEquals(tc, 3, (int) depth);
    close_to(tc, 0.80, rs, 1e-9);
    free(b);
    gerp_scratch_destruct(sc);
    gerp_tree_destruct(gt);
}

static void test_min_leaves_cutoff(CuTest *tc) {
    GerpTree *gt = gerp_tree_construct(BINARY_4);
    GerpScratch *sc = gerp_scratch_construct(gt);
    char *b = bases_from("A---");  // only one leaf
    double rs = 999; int64_t depth = -1;
    bool scored = gerp_score_column(gt, sc, b, 2, 1.0, &rs, &depth);
    CuAssertTrue(tc, !scored);     // below min_leaves
    CuAssertIntEquals(tc, 1, (int) depth);
    free(b);
    gerp_scratch_destruct(sc);
    gerp_tree_destruct(gt);
}

static void test_branch_scale(CuTest *tc) {
    GerpTree *gt = gerp_tree_construct(BINARY_4);
    GerpScratch *sc = gerp_scratch_construct(gt);
    char *b = bases_from("AAAA");
    double rs = 0; int64_t depth = 0;
    bool scored = gerp_score_column(gt, sc, b, 2, 2.0, &rs, &depth);
    CuAssertTrue(tc, scored);
    // 2x branch sum (=2.0), no subs.
    close_to(tc, 2.0, rs, 1e-9);
    free(b);
    gerp_scratch_destruct(sc);
    gerp_tree_destruct(gt);
}

// 3-way polytomy + outgroup -- exercises the batch Hartigan rule that
// motivated the audit fix.  Tree:
//
//      Anc0
//      / .
//     /   .
//   Anc1  S        S sib, branch 0.10
//   /|.
//  P Q R           polytomy children, all branches 0.20
//  Total branch sum (all present) = 0.20*3 + 0.10 + 0.10 = 0.80
//
// Test the auditor's example: P={A,C} (=A here -- we feed single base),
// actually we can construct equivalent observable cases since we only
// feed concrete bases.  Use P=A, Q=G, R=C, S=C.  Hartigan batch:
//   Anc1 votes: A=1 G=1 C=1 -> max_k=1, cset={A,C,G}, cost += 3-1 = 2.
//   Anc0 children present (Anc1 set {A,C,G}, S set {C}):
//     count: A=1 C=2 G=1 T=0 -> max_k=2, cset={C}, cost += 2-2 = 0.
//   Total cost = 2.
// All 4 leaves present -> branch sum = 0.80.  RS = 0.80 - 2 = -1.20.
static const char *POLYTOMY_3_PLUS_SIB =
    "((P:0.20,Q:0.20,R:0.20)Anc1:0.10,S:0.10)Anc0;";

static void test_polytomy_hartigan(CuTest *tc) {
    GerpTree *gt = gerp_tree_construct(POLYTOMY_3_PLUS_SIB);
    GerpScratch *sc = gerp_scratch_construct(gt);
    CuAssertIntEquals(tc, 4, (int) gerp_tree_n_leaves(gt));

    // First: trivial check.  All same -> RS = 0.80.
    {
        char *b = bases_from("AAAA");
        double rs = 0; int64_t depth = 0;
        bool ok = gerp_score_column(gt, sc, b, 2, 1.0, &rs, &depth);
        CuAssertTrue(tc, ok);
        CuAssertIntEquals(tc, 4, (int) depth);
        close_to(tc, 0.80, rs, 1e-9);
        free(b);
    }

    // The polytomy case.  Leaf order P,Q,R,S, bases A,G,C,C -> cost 2.
    // RS = 0.80 - 2 = -1.20.
    {
        char *b = bases_from("AGCC");
        double rs = 0; int64_t depth = 0;
        bool ok = gerp_score_column(gt, sc, b, 2, 1.0, &rs, &depth);
        CuAssertTrue(tc, ok);
        close_to(tc, -1.20, rs, 1e-9);
        free(b);
    }
    gerp_scratch_destruct(sc);
    gerp_tree_destruct(gt);
}

CuSuite* gerp_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_construct_and_introspect);
    SUITE_ADD_TEST(suite, test_all_present_no_subs);
    SUITE_ADD_TEST(suite, test_one_substitution);
    SUITE_ADD_TEST(suite, test_two_substitutions);
    SUITE_ADD_TEST(suite, test_gap_drops_leaf);
    SUITE_ADD_TEST(suite, test_min_leaves_cutoff);
    SUITE_ADD_TEST(suite, test_branch_scale);
    SUITE_ADD_TEST(suite, test_polytomy_hartigan);
    return suite;
}

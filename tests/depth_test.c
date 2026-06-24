/*
 * Tests for the GERP scoring core (taffy/impl/gerp.c -- Fitch parsimony +
 * branch-sum, gerp_score_column) PLUS end-to-end CLI tests of the `taffy depth`
 * driver: leaf depth, the chunked universal-column axis, --bin, and the opt-in
 * output guards, run against tests/evolverMammals.maf via st_system.
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

// Paralog union: a leaf with cset {A, T} (i.e. its paralog rows showed
// both A and T) should be scored as a polymorphic leaf -- Hartigan/Fitch
// counts both A and T as "votes" at the parent.  Compare against the
// single-bit case to verify the cset path costs strictly less when an
// extra paralog base agrees with another leaf.
//
// Same BINARY_4 tree (W,X)Anc1, (Y,Z)Anc2 / Anc0.  Set up so the only
// substitution is at leaf W:
//   X = A,  Y = A,  Z = A,  W = T  (1 substitution, cost=1, RS = 1.0-1 = 0.0)
//
// Now add a "paralog" base A to W: cset_W = {A, T}.  With A available at
// W, Hartigan picks A as the consensus at Anc1, and the substitution
// vanishes (cost=0, RS = 1.0 - 0 = 1.0).
static void test_paralog_union_cset(CuTest *tc) {
    GerpTree *gt = gerp_tree_construct(BINARY_4);
    GerpScratch *sc = gerp_scratch_construct(gt);
    int64_t nl = gerp_tree_n_leaves(gt);
    CuAssertIntEquals(tc, 4, (int) nl);

    // Look up leaf order: BINARY_4 is "((W,X)Anc1,(Y,Z)Anc2)Anc0".  We
    // don't expose the leaf-id mapping directly, but the bases_from helper
    // assumes the string position == leaf id; that's how the rest of the
    // suite uses it.  So map by position: W=0, X=1, Y=2, Z=3.

    // Single-bit reference: W=T, X=A, Y=A, Z=A  ->  cost=1, RS=0.0
    {
        uint8_t cs[4];
        cs[0] = gerp_base_to_bit('T');   // W
        cs[1] = gerp_base_to_bit('A');   // X
        cs[2] = gerp_base_to_bit('A');   // Y
        cs[3] = gerp_base_to_bit('A');   // Z
        double rs = 0; int64_t depth = 0;
        bool ok = gerp_score_column_csets(gt, sc, cs, 2, 1.0, &rs, &depth);
        CuAssertTrue(tc, ok);
        CuAssertIntEquals(tc, 4, (int) depth);
        close_to(tc, 0.0, rs, 1e-9);
    }

    // Union case: W's paralog rows showed both A and T -> cset_W = {A,T}.
    // Hartigan now finds an A-consistent reconstruction at Anc1, eliminating
    // the substitution -> cost=0, RS=1.0.
    {
        uint8_t cs[4];
        cs[0] = gerp_base_to_bit('A') | gerp_base_to_bit('T');  // W = {A,T}
        cs[1] = gerp_base_to_bit('A');                          // X
        cs[2] = gerp_base_to_bit('A');                          // Y
        cs[3] = gerp_base_to_bit('A');                          // Z
        double rs = 0; int64_t depth = 0;
        bool ok = gerp_score_column_csets(gt, sc, cs, 2, 1.0, &rs, &depth);
        CuAssertTrue(tc, ok);
        // Multi-bit leaves still count as one species in depth.
        CuAssertIntEquals(tc, 4, (int) depth);
        close_to(tc, 1.0, rs, 1e-9);
    }

    // All-paralog disagreement: W = {A,T}, X = {C,G}.  Anc1 votes:
    //   A=1 (W), C=1 (X), G=1 (X), T=1 (W); max_k=1, cset_Anc1 = {A,C,G,T},
    //   cost += 2-1 = 1.
    // Anc2: Y=A, Z=A -> cset={A}, cost=0.
    // Anc0 (children Anc1 cset {A,C,G,T} and Anc2 cset {A}): counts
    //   A=2 C=1 G=1 T=1; max_k=2, cset_Anc0 = {A}, cost += 2-2 = 0.
    // Total cost = 1.  All 4 leaves present -> RS = 1.0 - 1 = 0.0.
    {
        uint8_t cs[4];
        cs[0] = gerp_base_to_bit('A') | gerp_base_to_bit('T');
        cs[1] = gerp_base_to_bit('C') | gerp_base_to_bit('G');
        cs[2] = gerp_base_to_bit('A');
        cs[3] = gerp_base_to_bit('A');
        double rs = 0; int64_t depth = 0;
        bool ok = gerp_score_column_csets(gt, sc, cs, 2, 1.0, &rs, &depth);
        CuAssertTrue(tc, ok);
        CuAssertIntEquals(tc, 4, (int) depth);
        close_to(tc, 0.0, rs, 1e-9);
    }

    gerp_scratch_destruct(sc);
    gerp_tree_destruct(gt);
}

// Char-input wrapper still works (it's how all the older tests run).
// Spot-check that a leaf_bases call equals the corresponding cset call.
static void test_char_wrapper_matches_cset(CuTest *tc) {
    GerpTree *gt = gerp_tree_construct(BINARY_4);
    GerpScratch *sc = gerp_scratch_construct(gt);

    char *b = bases_from("ACGT");  // hits the 3-substitution path
    double rs_chars = 0; int64_t d_chars = 0;
    bool ok = gerp_score_column(gt, sc, b, 2, 1.0, &rs_chars, &d_chars);
    CuAssertTrue(tc, ok);
    free(b);

    uint8_t cs[4] = {
        gerp_base_to_bit('A'),
        gerp_base_to_bit('C'),
        gerp_base_to_bit('G'),
        gerp_base_to_bit('T'),
    };
    double rs_cs = 0; int64_t d_cs = 0;
    ok = gerp_score_column_csets(gt, sc, cs, 2, 1.0, &rs_cs, &d_cs);
    CuAssertTrue(tc, ok);

    close_to(tc, rs_chars, rs_cs, 1e-12);
    CuAssertIntEquals(tc, (int) d_chars, (int) d_cs);

    gerp_scratch_destruct(sc);
    gerp_tree_destruct(gt);
}

// --- CLI / driver tests: `taffy depth` against evolverMammals.maf ---------
// 5 leaves, referenced on Anc0, carries a `# hal` tree (so --depth needs no
// -t).  st_system runs the built ./bin/taffy and returns its exit code, as in
// stats_test.c / coverage_test.c.

// At least one of --rs / --depth is required.
static void test_depth_requires_output(CuTest *tc) {
    int rc = st_system("./bin/taffy depth -i ./tests/evolverMammals.maf "
                       ">/dev/null 2>&1");
    CuAssertTrue(tc, rc != 0);
}

// Universal-column output: 2e9-chunked chroms uni0,uni1,... monotone 4-col
// bedGraph; needs a .tui.  evolverMammals (T << 2e9) all lands in chunk uni0.
// Also exercises --sizes (uni0\t<size>).
static void test_depth_integer_coords_cli(CuTest *tc) {
    const char *maf = "./tests/evolverMammals.maf";
    st_system("rm -f %s.tui ./tests/depth.int.bg ./tests/depth.int.sizes", maf);  // clean slate
    // Universal-column output without a .tui errors out.
    int rc = st_system("./bin/taffy depth -i %s --bin 1000 "
                       "--depth ./tests/depth.int.bg >/dev/null 2>&1", maf);
    CuAssertTrue(tc, rc != 0);
    // Build the universal index, then --bin emits chunked chroms uni0,uni1,...
    CuAssertIntEquals(tc, 0, st_system("./bin/taffy index -i %s -u", maf));
    rc = st_system("./bin/taffy depth -i %s --bin 1000 "
                   "--minLeaves 1 --depth ./tests/depth.int.bg "
                   "--sizes ./tests/depth.int.sizes", maf);
    CuAssertIntEquals(tc, 0, rc);
    CuAssertIntEquals(tc, 0, st_system(
        "awk '$1!~/^uni[0-9]+$/||NF!=4||$3<=$2{exit 1}' ./tests/depth.int.bg"));
    // --sizes lists at least uni0 with a positive size; evolverMammals fits one chunk.
    CuAssertIntEquals(tc, 0, st_system("test -s ./tests/depth.int.sizes"));
    CuAssertIntEquals(tc, 0, st_system(
        "awk '$1!~/^uni[0-9]+$/||NF!=2||$2<=0{exit 1}' ./tests/depth.int.sizes"));
    CuAssertIntEquals(tc, 0, st_system("grep -q '^uni0\t' ./tests/depth.int.sizes"));
    st_system("rm -f %s.tui ./tests/depth.int.bg ./tests/depth.int.sizes", maf);
}

CuSuite* depth_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_construct_and_introspect);
    SUITE_ADD_TEST(suite, test_all_present_no_subs);
    SUITE_ADD_TEST(suite, test_one_substitution);
    SUITE_ADD_TEST(suite, test_two_substitutions);
    SUITE_ADD_TEST(suite, test_gap_drops_leaf);
    SUITE_ADD_TEST(suite, test_min_leaves_cutoff);
    SUITE_ADD_TEST(suite, test_branch_scale);
    SUITE_ADD_TEST(suite, test_polytomy_hartigan);
    SUITE_ADD_TEST(suite, test_paralog_union_cset);
    SUITE_ADD_TEST(suite, test_char_wrapper_matches_cset);
    SUITE_ADD_TEST(suite, test_depth_requires_output);
    SUITE_ADD_TEST(suite, test_depth_integer_coords_cli);
    return suite;
}

/*
 * Unit tests for taffy/impl/view_chain_dup_filter.c -- the per-region
 * chain-based duplicate-row filter that drives taffy view's
 * --chainDupFilter[=N].  Focused on the per-block dup semantics
 * (genome appears 2+ times in same block triggers chain pass; lone
 * occurrences are never dropped), top_n selection, and row-0
 * pinning.
 *
 *  Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "view_chain_dup_filter.h"
#include "sonLib.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Build a single Alignment_Row.  Bases left empty (the filter only
 * looks at sequence_name + coords). */
static Alignment_Row *mk_row(const char *name, int64_t start,
                             int64_t length, int64_t seqlen, bool strand) {
    Alignment_Row *r = st_calloc(1, sizeof(Alignment_Row));
    r->sequence_name = stString_copy(name);
    r->start = start;
    r->length = length;
    r->sequence_length = seqlen;
    r->strand = strand;
    r->bases = stString_copy("");
    return r;
}

/* Build an Alignment from a NULL-terminated array of row pointers.
 * Rows are linked in order; row[0] is row-0 (reference, always pinned). */
static Alignment *mk_aln(Alignment_Row **rows, int n) {
    Alignment *a = st_calloc(1, sizeof(Alignment));
    a->row_number = n;
    a->column_number = 0;
    /* alignment_destruct asserts column_tags != NULL even when
     * column_number == 0; satisfy it with an empty allocation. */
    a->column_tags = st_calloc(1, sizeof(Tag*));
    for (int i = 0; i < n - 1; i++) rows[i]->n_row = rows[i + 1];
    a->row = rows[0];
    return a;
}

/* Count rows in an alignment's row chain. */
static int64_t count_rows(Alignment *a) {
    int64_t c = 0;
    for (Alignment_Row *r = a->row; r != NULL; r = r->n_row) c++;
    return c;
}

/* Find first row whose name first-dot prefix matches `genome`.  NULL if
 * none.  (Used to check row-0 pinning and surviving paralog identity.) */
static Alignment_Row *find_row_by_genome(Alignment *a, const char *genome) {
    size_t gl = strlen(genome);
    for (Alignment_Row *r = a->row; r != NULL; r = r->n_row) {
        if (strncmp(r->sequence_name, genome, gl) == 0 &&
            (r->sequence_name[gl] == '.' || r->sequence_name[gl] == 0)) {
            return r;
        }
    }
    return NULL;
}

/* ------------------------------------------------------------------- */

/* Per-block paralog filter: same genome appears twice in same block,
 * the chain pass picks the longer-chain paralog. */
static void test_view_chain_dup_filter_per_block_paralog(CuTest *tc) {
    /* Block A: row-0 (Anc) + dog.chr1 (paralog 1, score 100) + dog.chr2 (paralog 2, score 50) */
    Alignment_Row *aA[3] = {
        mk_row("Anc.refA", 0, 100, 10000, 1),
        mk_row("dog.chr1", 0, 100, 100000, 1),     /* primary */
        mk_row("dog.chr2", 5000, 50, 100000, 1),   /* paralog */
    };
    /* Block B: row-0 (Anc) + dog.chr1 (continuation of primary; same chain) */
    Alignment_Row *aB[2] = {
        mk_row("Anc.refA", 100, 100, 10000, 1),
        mk_row("dog.chr1", 100, 100, 100000, 1),
    };
    Alignment *blkA = mk_aln(aA, 3);
    Alignment *blkB = mk_aln(aB, 2);
    stList *blocks = stList_construct();
    stList_append(blocks, blkA);
    stList_append(blocks, blkB);

    view_chain_dup_filter(blocks, /*top_n=*/1);

    /* Block A: row-0 pinned + ONE dog row (the primary chain wins).
     * Block B unchanged (no paralog dup in it). */
    CuAssertIntEquals(tc, 2, (int) count_rows(blkA));
    CuAssertIntEquals(tc, 2, (int) count_rows(blkB));
    CuAssertPtrNotNull(tc, find_row_by_genome(blkA, "Anc"));      /* row-0 */
    Alignment_Row *survivor = find_row_by_genome(blkA, "dog");
    CuAssertPtrNotNull(tc, survivor);
    /* dog.chr1 is the primary (it abuts dog.chr1 in block B, longer chain). */
    CuAssertStrEquals(tc, "dog.chr1", survivor->sequence_name);

    alignment_destruct(blkA, 1);
    alignment_destruct(blkB, 1);
    stList_destruct(blocks);
}

/* Singleton genomes (one row per block, no per-block dup) are NEVER dropped,
 * even if their across-block alns can't be chained (tiny + gaps). */
static void test_view_chain_dup_filter_singletons_preserved(CuTest *tc) {
    /* Block A: row-0 + cat (1 row, length 100). */
    Alignment_Row *aA[2] = {
        mk_row("Anc.refA", 0, 100, 10000, 1),
        mk_row("cat.chrX", 0, 100, 50000, 1),
    };
    /* Block B: row-0 + cat (length 1, t-gap=999 from block A so chain breaks). */
    Alignment_Row *aB[2] = {
        mk_row("Anc.refA", 100, 1, 10000, 1),
        mk_row("cat.chrX", 1099, 1, 50000, 1),    /* gap of 999 bp */
    };
    Alignment *blkA = mk_aln(aA, 2);
    Alignment *blkB = mk_aln(aB, 2);
    stList *blocks = stList_construct();
    stList_append(blocks, blkA);
    stList_append(blocks, blkB);

    view_chain_dup_filter(blocks, 1);

    /* Both blocks unchanged: cat appears once per block; not a per-block
     * dup; the filter must NOT touch it regardless of chain breakage. */
    CuAssertIntEquals(tc, 2, (int) count_rows(blkA));
    CuAssertIntEquals(tc, 2, (int) count_rows(blkB));

    alignment_destruct(blkA, 1);
    alignment_destruct(blkB, 1);
    stList_destruct(blocks);
}

/* Row-0 is always preserved, even when its genome happens to appear
 * elsewhere as a paralog. */
static void test_view_chain_dup_filter_row0_pinned(CuTest *tc) {
    /* Block A: row-0 = hg.chrA + hg.chrB (sibling paralog of the same
     * genome).  The filter should keep row-0 (hg.chrA) and may drop
     * the sibling depending on chain.  Building 2 blocks so hg has 4
     * total occurrences (and 2 per block) -> triggers chain pass. */
    Alignment_Row *aA[2] = {
        mk_row("hg.chrA", 0,    100, 1000, 1),
        mk_row("hg.chrB", 0,    100, 1000, 1),
    };
    Alignment_Row *aB[2] = {
        mk_row("hg.chrA", 100,  100, 1000, 1),
        mk_row("hg.chrB", 100,  100, 1000, 1),
    };
    Alignment *blkA = mk_aln(aA, 2);
    Alignment *blkB = mk_aln(aB, 2);
    stList *blocks = stList_construct();
    stList_append(blocks, blkA);
    stList_append(blocks, blkB);

    view_chain_dup_filter(blocks, 1);

    /* Row-0 (hg.chrA in each block) must survive. */
    CuAssertStrEquals(tc, "hg.chrA", blkA->row->sequence_name);
    CuAssertStrEquals(tc, "hg.chrA", blkB->row->sequence_name);
    /* hg has 2 rows per block before, top-1 picks the primary chain of
     * the OTHER hg row (the sibling).  Both chains are colinear length-
     * 200, so chains[0] = whichever wins the tie; the loser is dropped.
     * Result: each block has 1 row (row-0 only) since the sibling lost. */
    /* (We assert just the row-0-pinned invariant; the tie-breaker for
     * which sibling chain wins is incidental.) */

    alignment_destruct(blkA, 1);
    alignment_destruct(blkB, 1);
    stList_destruct(blocks);
}

/* top_n >= n_chains: every row survives (no-op even with per-block dups). */
static void test_view_chain_dup_filter_top_n_unbounded(CuTest *tc) {
    Alignment_Row *aA[3] = {
        mk_row("Anc.refA", 0, 100, 10000, 1),
        mk_row("dog.chr1", 0, 100, 100000, 1),
        mk_row("dog.chr2", 5000, 50, 100000, 1),
    };
    Alignment *blkA = mk_aln(aA, 3);
    stList *blocks = stList_construct();
    stList_append(blocks, blkA);

    /* Both paralogs form distinct chains; top_n=100 keeps all. */
    view_chain_dup_filter(blocks, 100);

    CuAssertIntEquals(tc, 3, (int) count_rows(blkA));

    alignment_destruct(blkA, 1);
    stList_destruct(blocks);
}

/* Empty input: no crash, no work. */
static void test_view_chain_dup_filter_empty(CuTest *tc) {
    stList *blocks = stList_construct();
    view_chain_dup_filter(blocks, 1);
    CuAssertIntEquals(tc, 0, (int) stList_length(blocks));
    stList_destruct(blocks);
}

/* top_n < 1 is clamped to 1 (no crash, no underflow). */
static void test_view_chain_dup_filter_top_n_zero_clamps(CuTest *tc) {
    Alignment_Row *aA[3] = {
        mk_row("Anc.refA", 0, 100, 10000, 1),
        mk_row("dog.chr1", 0, 100, 100000, 1),
        mk_row("dog.chr2", 5000, 50, 100000, 1),
    };
    Alignment *blkA = mk_aln(aA, 3);
    stList *blocks = stList_construct();
    stList_append(blocks, blkA);

    /* top_n=0 should behave like top_n=1 (clamped). */
    view_chain_dup_filter(blocks, 0);

    /* Row-0 + one survivor = 2. */
    CuAssertIntEquals(tc, 2, (int) count_rows(blkA));

    alignment_destruct(blkA, 1);
    stList_destruct(blocks);
}

CuSuite *view_chain_dup_filter_test_suite(void) {
    CuSuite *s = CuSuiteNew();
    SUITE_ADD_TEST(s, test_view_chain_dup_filter_per_block_paralog);
    SUITE_ADD_TEST(s, test_view_chain_dup_filter_singletons_preserved);
    SUITE_ADD_TEST(s, test_view_chain_dup_filter_row0_pinned);
    SUITE_ADD_TEST(s, test_view_chain_dup_filter_top_n_unbounded);
    SUITE_ADD_TEST(s, test_view_chain_dup_filter_empty);
    SUITE_ADD_TEST(s, test_view_chain_dup_filter_top_n_zero_clamps);
    return s;
}

#include "CuTest.h"
#include "block_reader.h"
#include "sonLib.h"

// Walk a file with BlockReader and return (block_count, total_columns, first_ref_name).
static void walk_with_block_reader(const char *path, int *block_count, int64_t *total_cols, char **first_ref) {
    FILE *fh = fopen(path, "r");
    LI *li = LI_construct(fh);
    BlockReader *r = block_reader_open(li);
    Tag *tag = block_reader_take_header(r);
    tag_destruct(tag);

    *block_count = 0;
    *total_cols = 0;
    *first_ref = NULL;
    Alignment *prev = NULL;
    Alignment *a;
    while ((a = block_reader_next(r, prev)) != NULL) {
        if (*block_count == 0) {
            *first_ref = stString_copy(a->row->sequence_name);
        }
        *total_cols += a->column_number;
        (*block_count)++;
        if (prev != NULL) alignment_destruct(prev, 1);
        prev = a;
    }
    if (prev != NULL) alignment_destruct(prev, 1);
    block_reader_destruct(r);
    LI_destruct(li);
    fclose(fh);
}

static void test_block_reader_maf_input(CuTest *tc) {
    int n; int64_t cols; char *ref;
    walk_with_block_reader("./tests/evolverMammals.maf", &n, &cols, &ref);
    CuAssertTrue(tc, n > 0);
    CuAssertTrue(tc, cols > 0);
    CuAssertStrEquals(tc, "Anc0.Anc0refChr0", ref);
    free(ref);
}

static void test_block_reader_taf_input(CuTest *tc) {
    // Build a TAF on the side from the same MAF and walk it
    int rc = st_system("./bin/taffy view -i ./tests/evolverMammals.maf -o ./tests/block_reader_test.taf");
    CuAssertIntEquals(tc, 0, rc);

    int n; int64_t cols; char *ref;
    walk_with_block_reader("./tests/block_reader_test.taf", &n, &cols, &ref);
    CuAssertTrue(tc, n > 0);
    CuAssertTrue(tc, cols > 0);
    CuAssertStrEquals(tc, "Anc0.Anc0refChr0", ref);
    free(ref);

    st_system("rm -f ./tests/block_reader_test.taf");
}

static void test_block_reader_maf_and_taf_match(CuTest *tc) {
    // BlockReader walking a MAF and walking the equivalent TAF should yield the
    // same sequence of (block_count, total_columns) -- this is the contract
    // that lets downstream tools accept either format transparently.
    int rc = st_system("./bin/taffy view -i ./tests/evolverMammals.maf -o ./tests/block_reader_test.taf");
    CuAssertIntEquals(tc, 0, rc);

    int n_maf, n_taf; int64_t cols_maf, cols_taf; char *ref_maf, *ref_taf;
    walk_with_block_reader("./tests/evolverMammals.maf", &n_maf, &cols_maf, &ref_maf);
    walk_with_block_reader("./tests/block_reader_test.taf", &n_taf, &cols_taf, &ref_taf);

    CuAssertIntEquals(tc, n_maf, n_taf);
    CuAssertIntEquals(tc, cols_maf, cols_taf);
    CuAssertStrEquals(tc, ref_maf, ref_taf);
    free(ref_maf); free(ref_taf);

    st_system("rm -f ./tests/block_reader_test.taf");
}

static void test_block_reader_unrecognized_format(CuTest *tc) {
    // An input that isn't MAF or TAF should make block_reader_open return NULL.
    char *junk_path = "./tests/block_reader_junk.txt";
    FILE *junk = fopen(junk_path, "w");
    fprintf(junk, "this is not a maf or taf file\n");
    fclose(junk);

    FILE *fh = fopen(junk_path, "r");
    LI *li = LI_construct(fh);
    BlockReader *r = block_reader_open(li);
    CuAssertPtrEquals(tc, NULL, r);
    LI_destruct(li);
    fclose(fh);
    st_system("rm -f %s", junk_path);
}

static void test_block_reader_header_take(CuTest *tc) {
    // taf_read_header_2 picks up the run_length_encode_bases flag; verify the
    // reader exposes it correctly when the input is RLE-encoded TAF.
    int rc = st_system("./bin/taffy view -i ./tests/evolverMammals.maf -u -o ./tests/block_reader_test_rle.taf");
    CuAssertIntEquals(tc, 0, rc);

    FILE *fh = fopen("./tests/block_reader_test_rle.taf", "r");
    LI *li = LI_construct(fh);
    BlockReader *r = block_reader_open(li);
    CuAssertPtrNotNull(tc, r);
    CuAssertTrue(tc, !block_reader_is_maf(r));
    CuAssertTrue(tc, block_reader_run_length_encoded(r));

    Tag *tag = block_reader_take_header(r);
    CuAssertPtrNotNull(tc, tag);
    // After take, a second take returns NULL
    CuAssertPtrEquals(tc, NULL, block_reader_take_header(r));
    tag_destruct(tag);

    block_reader_destruct(r);
    LI_destruct(li);
    fclose(fh);
    st_system("rm -f ./tests/block_reader_test_rle.taf");
}

CuSuite *block_reader_test_suite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_block_reader_maf_input);
    SUITE_ADD_TEST(suite, test_block_reader_taf_input);
    SUITE_ADD_TEST(suite, test_block_reader_maf_and_taf_match);
    SUITE_ADD_TEST(suite, test_block_reader_unrecognized_format);
    SUITE_ADD_TEST(suite, test_block_reader_header_take);
    return suite;
}

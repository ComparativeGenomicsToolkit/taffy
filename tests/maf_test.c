#include "CuTest.h"
#include "taf.h"
#include "sonLib.h"

static void test_maf(CuTest *testCase, bool use_compression) {
    // Example maf file
    char *example_file = "./tests/evolverMammals.maf";
    char *temp_copy = "./tests/evolverMammals_copy.maf";

    // Read a maf and write a copy of it
    FILE *file = fopen(example_file, "r");
    LW *lw = LW_construct(fopen(temp_copy, "w"), 0);
    Alignment *alignment;
    LI *li = LI_construct(file);
    while((alignment = maf_read_block(li)) != NULL) {
        //maf_write_block(alignment, stdout);
        maf_write_block(alignment, lw);
        alignment_destruct(alignment, 1);
    }
    fclose(file);
    LI_destruct(li);
    LW_destruct(lw, 1);

    // Now check that all the maf blocks are the same
    file = fopen(example_file, "r");
    FILE *file_copy = fopen(temp_copy, "r");
    Alignment *alignment2;
    li = LI_construct(file);
    //st_uglyf("Before\n");
    LI *li_copy = LI_construct(file_copy);
    while((alignment = maf_read_block(li)) != NULL) {
        alignment2 = maf_read_block(li_copy);
        CuAssertTrue(testCase, alignment2 != NULL);
        // Check that the blocks are the same

        Alignment_Row *row = alignment->row, *row2 = alignment2->row;
        while(row != NULL) {
            CuAssertTrue(testCase, row2 != NULL);
            CuAssertStrEquals(testCase, row->sequence_name, row2->sequence_name);
            CuAssertIntEquals(testCase, row->start, row2->start);
            CuAssertIntEquals(testCase, row->length, row2->length);
            CuAssertIntEquals(testCase, row->sequence_length, row2->sequence_length);
            CuAssertIntEquals(testCase, row->strand, row2->strand);
            CuAssertStrEquals(testCase, row->bases, row2->bases);
            row = row->n_row; row2 = row2->n_row;
        }
        CuAssertTrue(testCase, row2 == NULL);

        alignment_destruct(alignment, 1);
        alignment_destruct(alignment2, 1);
    }
    CuAssertTrue(testCase, maf_read_block(li_copy) == NULL);    
    LI_destruct(li);
    LI_destruct(li_copy);

    fclose(file);
    st_system("rm %s", temp_copy);
}

static void test_maf_with_compression(CuTest *testCase) {
    test_maf(testCase, 1);
}

static void test_maf_without_compression(CuTest *testCase) {
    test_maf(testCase, 0);
}

CuSuite* maf_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_maf_with_compression);
    SUITE_ADD_TEST(suite, test_maf_without_compression);
    return suite;
}

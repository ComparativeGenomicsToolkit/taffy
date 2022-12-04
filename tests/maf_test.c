#include "CuTest.h"
#include "taf.h"
#include "sonLib.h"

static void test_maf(CuTest *testCase) {
    // Example maf file
    char *example_file = "./tests/evolverMammals.maf";
    char *temp_copy = "./tests/evolverMammals_copy.maf";

    // Read a maf and write a copy of it
    FILE *file = fopen(example_file, "r");
    FILE *out_file = fopen(temp_copy, "w");
    Alignment *alignment;
    LI *li = LI_construct(file);
    while((alignment = maf_read_block(li)) != NULL) {
        //maf_write_block(alignment, stdout);
        maf_write_block(alignment, out_file);
        alignment_destruct(alignment, 1);
    }
    fclose(file);
    fclose(out_file);
    LI_destruct(li);

    // Now check that all the maf blocks are the same
    file = fopen(example_file, "r");
    FILE *file_copy = fopen(temp_copy, "r");
    Alignment *alignment2;
    li = LI_construct(file);
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
}

CuSuite* maf_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_maf);
    return suite;
}

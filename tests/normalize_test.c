#include "CuTest.h"
#include "taf.h"

static void test_normalize(CuTest *testCase) {
    // Example maf file
    char *example_file = "./tests/evolverMammals.maf";
    // Read a maf and write a copy of it
    FILE *file = fopen(example_file, "r");
    Alignment *alignment, *p_alignment = NULL;
    while((alignment = maf_read_block(file)) != NULL) {
        if(p_alignment != NULL) {
            // Try merging the blocks
            alignment_link_adjacent(p_alignment, alignment, 0);
            if((alignment_length(p_alignment) < 50 || alignment_length(alignment) < 50) &&
                alignment_total_gap_length(p_alignment) < 50) {

                p_alignment = alignment_merge_adjacent(p_alignment, alignment);

            }
            else {
                alignment_destruct(p_alignment);
                p_alignment = alignment;
            }
        }
        else {
            p_alignment = alignment;
        }
    }
    if(p_alignment != NULL) {
        //maf_write_block(p_alignment, stderr);
        alignment_destruct(p_alignment);
    }
    fclose(file);

}

CuSuite* normalize_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_normalize);
    return suite;
}

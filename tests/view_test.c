#include "CuTest.h"
#include "taf.h"
#include "sonLib.h"
#include "bioioC.h"


static void test_paf(CuTest *testCase) {
    {
        char *example_file = "./tests/paf_test.maf";
        char *output_file = "./tests/paf_test.maf.out.paf";
        char *truth_file = "./tests/paf_test.maf.paf";
        int i = st_system("./bin/taffy view -i %s -p > %s", example_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames        
        st_system("rm -f %s", output_file);
    }

    {
        char *example_file = "./tests/paf_test_flip.maf";
        char *output_file = "./tests/paf_test_flip.maf.out.paf";
        char *truth_file = "./tests/paf_test.maf.paf";
        int i = st_system("./bin/taffy view -i %s -p > %s", example_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames                
        st_system("rm -f %s", output_file);        
    }

    {
        char *example_file = "./tests/paf_test_rev.maf";
        char *output_file = "./tests/paf_test_rev.maf.out.paf";
        char *truth_file = "./tests/paf_test_rev.maf.paf";
        int i = st_system("./bin/taffy view -i %s -p > %s", example_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames                
        st_system("rm -f %s", output_file);
    }

    {
        char *example_file = "./tests/paf_test_rev_flip.maf";
        char *output_file = "./tests/paf_test_rev_flip.maf.out.paf";
        char *truth_file = "./tests/paf_test_rev.maf.paf";
        int i = st_system("./bin/taffy view -i %s -p > %s", example_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames                
        st_system("rm -f %s", output_file);
    }

}

static void test_lineage_diffs(CuTest *testCase) {
    {
        char *example_file = "./tests/evolverMammals.maf.mini";
        char *output_file = "./tests/lineage_diffs_test.maf.out";
        char *truth_file = "./tests/evolverMammals.maf.mini.lineage_diffs";
        int i = st_system("./bin/taffy view -i %s -t '((simHuman_chr6:0.144018,(simMouse_chr6:0.084509,simRat_chr6:0.091589)mr:0.271974)Anc1:0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303)Anc2:0.032898)Anc0;' -b -m > %s",
                          example_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames
        st_system("rm -f %s", output_file);
    }
}

static void test_ref_diffs(CuTest *testCase) {
    {
        char *example_file = "./tests/evolverMammals.maf.mini";
        char *output_file = "./tests/lineage_diffs_test.maf.out";
        char *truth_file = "./tests/evolverMammals.maf.mini.ref_diffs";
        int i = st_system("./bin/taffy view -i %s -a -m > %s",
                          example_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames
        st_system("rm -f %s", output_file);
    }
}

CuSuite* view_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_paf);
    SUITE_ADD_TEST(suite, test_lineage_diffs);
    SUITE_ADD_TEST(suite, test_ref_diffs);
    return suite;
}

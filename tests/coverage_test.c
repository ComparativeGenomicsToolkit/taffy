#include "CuTest.h"
#include "taf.h"
#include "sonLib.h"
#include "bioioC.h"


static void test_coverage(CuTest *testCase) {
    {
        char *example_file = "./tests/coverage_test.maf";
        char *output_file = "./tests/coverage_test.maf.coverage.tsv";
        char *truth_file = "./tests/coverage_test.truth.tsv";
        int i = st_system("./bin/taffy view -i %s | ./bin/taffy coverage -g cat.a > %s", example_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames        
        st_system("rm -f %s", output_file);
    }
    {
        char *example_file = "./tests/coverage_test.maf";
        char *output_file = "./tests/coverage_test.a0.maf.coverage.tsv";
        char *truth_file = "./tests/coverage_test.a0.truth.tsv";
        int i = st_system("./bin/taffy view -i %s | ./bin/taffy coverage -a 0 > %s", example_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames        
        st_system("rm -f %s", output_file);
    }
    {
        char *example_file = "./tests/coverage_test.maf";
        char *output_file = "./tests/coverage_test.subset.maf.coverage.tsv";
        char *truth_file = "./tests/coverage_test.subset.truth.tsv";
        int i = st_system("./bin/taffy view -i %s | ./bin/taffy coverage --sexChr dog.chr1 > %s", example_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames        
        st_system("rm -f %s", output_file);      
    }
}

// Verifies that taffy coverage -i x.maf produces identical output to
// taffy view -i x.maf | taffy coverage. This is the per-tool dual-input
// regression test for the BlockReader migration.
static void test_coverage_maf_input(CuTest *testCase) {
    char *example_file = "./tests/coverage_test.maf";
    char *piped_out = "./tests/coverage_test.piped.tsv";
    char *direct_out = "./tests/coverage_test.direct.tsv";
    int i = st_system("./bin/taffy view -i %s | ./bin/taffy coverage -g cat.a > %s", example_file, piped_out);
    CuAssertIntEquals(testCase, 0, i);
    int j = st_system("./bin/taffy coverage -i %s -g cat.a > %s", example_file, direct_out);
    CuAssertIntEquals(testCase, 0, j);
    int diff_ret = st_system("diff %s %s", piped_out, direct_out);
    CuAssertIntEquals(testCase, 0, diff_ret);
    st_system("rm -f %s %s", piped_out, direct_out);
}

CuSuite* coverage_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_coverage);
    SUITE_ADD_TEST(suite, test_coverage_maf_input);
    return suite;
}

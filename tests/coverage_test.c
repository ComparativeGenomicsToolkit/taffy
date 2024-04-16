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

CuSuite* coverage_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_coverage);
    return suite;
}

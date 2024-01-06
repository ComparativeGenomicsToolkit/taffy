#include "CuTest.h"
#include "taf.h"
#include "sonLib.h"
#include "bioioC.h"


static void test_sort(CuTest *testCase) {
    {
        char *example_file = "./tests/evolverMammals.maf";
        char *output_file = "./tests/sort_test.taf.out";
        char *sort_file = "./tests/sort_file.txt";
        int i = st_system("./bin/taffy view -i %s | ./bin/taffy sort -n %s > %s", example_file, sort_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames        
        st_system("rm -f %s", output_file);
    }
}

CuSuite* sort_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_sort);
    return suite;
}

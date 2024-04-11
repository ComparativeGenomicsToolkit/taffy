#include "CuTest.h"
#include "taf.h"
#include "sonLib.h"
#include "bioioC.h"

static void test_sort(CuTest *testCase) {
    {
        char *example_file = "./tests/evolverMammals.maf.mini";
        char *output_file = "./tests/sort_test.maf.out";
        char *truth_file = "./tests/evolverMammals.maf.mini.sorted";
        char *sort_file = "./tests/sort_file.txt";
        int i = st_system("./bin/taffy view -i %s | ./bin/taffy sort -n %s | ./bin/taffy view -m > %s",
                          example_file, sort_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames        
        st_system("rm -f %s", output_file);
    }
}

static void test_sort_ignore_first_row(CuTest *testCase) {
    {
        char *example_file = "./tests/evolverMammals.maf.mini";
        char *output_file = "./tests/sort_test.maf.out";
        char *truth_file = "./tests/evolverMammals.maf.mini.sorted.ref.ignored";
        char *sort_file = "./tests/sort_file.txt";
        int i = st_system("./bin/taffy view -i %s | ./bin/taffy sort -n %s -r | ./bin/taffy view -m > %s",
                          example_file, sort_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames
        st_system("rm -f %s", output_file);
    }
}

static void test_filter(CuTest *testCase) {
    {
        char *example_file = "./tests/evolverMammals.maf.mini";
        char *output_file = "./tests/filter_test.maf.out";
        char *truth_file = "./tests/evolverMammals.maf.mini.filtered";
        char *filter_file = "./tests/filter_file.txt";
        int i = st_system("./bin/taffy view -i %s | ./bin/taffy sort -f %s | ./bin/taffy view -m > %s",
                          example_file, filter_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames
        st_system("rm -f %s", output_file);
    }
}

static void test_filter_ignore_first_row(CuTest *testCase) {
    {
        char *example_file = "./tests/evolverMammals.maf.mini";
        char *output_file = "./tests/filter_test.maf.out";
        char *truth_file = "./tests/evolverMammals.maf.mini.filtered.ref.ignored";
        char *filter_file = "./tests/filter_file.txt";
        int i = st_system("./bin/taffy view -i %s | ./bin/taffy sort -f %s -r | ./bin/taffy view -m > %s",
                          example_file, filter_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames
        st_system("rm -f %s", output_file);
    }
}

static void test_sort_filter_pad_and_dup_filter(CuTest *testCase) {
    {
        char *example_file = "./tests/evolverMammals.maf.mini";
        char *output_file = "./tests/sort_test.maf.out";
        char *truth_file = "./tests/evolverMammals.maf.mini.sorted.padded.dup_filtered";
        char *sort_file = "./tests/sort_file.txt";
        char *filter_file = "./tests/filter_file_2.txt";
        int i = st_system("./bin/taffy view -i %s | ./bin/taffy sort -n %s -f %s -p %s -d %s | ./bin/taffy view -m > %s",
                          example_file, sort_file, filter_file, sort_file, sort_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames
        st_system("rm -f %s", output_file);
    }
}

CuSuite* sort_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_sort);
    SUITE_ADD_TEST(suite, test_sort_ignore_first_row);
    SUITE_ADD_TEST(suite, test_filter);
    SUITE_ADD_TEST(suite, test_filter_ignore_first_row);
    SUITE_ADD_TEST(suite, test_sort_filter_pad_and_dup_filter);
    return suite;
}

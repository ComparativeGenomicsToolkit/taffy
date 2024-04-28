#include "CuTest.h"
#include "taf.h"
#include "sonLib.h"


static void test_wiggle(CuTest *testCase) {
    // Example wiggle parsing
    char *example_wig = "./tests/evolverMammals.wig.mini";
    stHash *wig = wig_parse(example_wig, "Anc0.");
    CuAssertTrue(testCase, wig_get_value(wig, "Anc0.Anc0refChr0", 0, -1) == 5.1);
    CuAssertTrue(testCase, wig_get_value(wig, "Anc0.Anc0refChr0", 1, -1) == 10);
    CuAssertTrue(testCase, wig_get_value(wig, "Anc0.Anc0refChr0", 2, -1) == 67);
    CuAssertTrue(testCase, wig_get_value(wig, "Anc0.Anc0refChr0", 3, -1) == -1);
    CuAssertTrue(testCase, wig_get_value(wig, "Anc0.Anc0refChr0", 49, -1) == -1);
    CuAssertTrue(testCase, wig_get_value(wig, "Anc0.Anc0refChr0", 50, -1) == 7);
    CuAssertTrue(testCase, wig_get_value(wig, "Anc0.Anc0refChr0", 51, -1) == 7);
    CuAssertTrue(testCase, wig_get_value(wig, "Anc0.Anc0refChr0", 52, -1) == -1);
    CuAssertTrue(testCase, wig_get_value(wig, "Anc0.Anc0refChr0", 54, -1) == -1);
    CuAssertTrue(testCase, wig_get_value(wig, "Anc0.Anc0refChr0", 55, -1) == 12);
    CuAssertTrue(testCase, wig_get_value(wig, "Anc0.Anc0refChr0", 56, -1) == 12);
    CuAssertTrue(testCase, wig_get_value(wig, "Anc0.Anc0refChr0", 57, -1) == -1);
    stHash_destruct(wig);
}

static void test_large_wiggle(CuTest *testCase) {
    // Example wiggle parsing
    char *example_wig = "./tests/447-way/447-mammalian-2022v1_hg38_chr22_22000000_22100000.phyloP.wig";
    stHash *wig = wig_parse(example_wig, "hg38.");
    CuAssertTrue(testCase, wig_get_value(wig, "hg38.chr22", 22000001, -1) == -1.301);
    CuAssertTrue(testCase, wig_get_value(wig, "hg38.chr22", 22000002, -1) == 0.663);
    CuAssertTrue(testCase, wig_get_value(wig, "hg38.chr22", 22000003, -1) == 0.663);
    CuAssertTrue(testCase, wig_get_value(wig, "hg38.chr22", 22000004, -1) == -0.1);
    stHash_destruct(wig);
}

static void test_annotate(CuTest *testCase) {
    {
        char *example_file = "./tests/evolverMammals.maf.mini";
        char *example_wig = "./tests/evolverMammals.wig.mini";
        char *output_file = "./tests/annotate_test.taf.out";
        char *truth_file = "./tests/evolverMammals.taf.mini.annotated";
        int i = st_system("./bin/taffy view -i %s | ./bin/taffy annotate -w %s -t test_label -r 'Anc0.' > %s",
                          example_file, example_wig, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames
        st_system("rm -f %s", output_file);
    }
}

CuSuite* wiggle_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_wiggle);
    SUITE_ADD_TEST(suite, test_large_wiggle);
    SUITE_ADD_TEST(suite, test_annotate);
    return suite;
}

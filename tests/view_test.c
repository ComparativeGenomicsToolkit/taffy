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
    {
        char *example_file = "./tests/paf_test_gapcol.maf";
        char *output_file = "./tests/paf_test_gapcol.maf.out.paf";
        char *truth_file = "./tests/paf_test_gapcol.maf.paf";
        int i = st_system("./bin/taffy view -i %s -p > %s", example_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames                
        st_system("rm -f %s", output_file);        
    }
    {
        // note this was verified via paftools.js view -f maf paf_test_gapcol.maf.cs.paf
        char *example_file = "./tests/paf_test_gapcol.maf";
        char *output_file = "./tests/paf_test_gapcol.maf.out.cs.paf";
        char *truth_file = "./tests/paf_test_gapcol.maf.cs.paf";
        int i = st_system("./bin/taffy view -i %s -p -C > %s", example_file, output_file);
        CuAssertIntEquals(testCase, 0, i); // return value should be zero
        int diff_ret = st_system("diff %s %s", output_file, truth_file);
        CuAssertIntEquals(testCase, 0, diff_ret); // return value should be zero if files sames                
        st_system("rm -f %s", output_file);        
    }

}

CuSuite* view_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_paf);
    return suite;
}

#include "CuTest.h"
#include "taf.h"
#include "sonLib.h"

// Verifies that taffy stats accepts MAF input directly and produces output identical to
// the equivalent `taffy view -i x.maf | taffy stats` invocation. This is the per-tool
// dual-input regression test for the BlockReader migration. Covers the three stats modes:
//  -a (alignment stats)
//  -b (sequence intervals) -- this previously errored out on MAF input
static void test_stats_maf_input(CuTest *tc) {
    char *maf = "./tests/evolverMammals.maf";

    // -a alignment stats
    int rc;
    rc = st_system("./bin/taffy stats -a -i %s > ./tests/stats.direct.a.txt", maf);
    CuAssertIntEquals(tc, 0, rc);
    rc = st_system("./bin/taffy view -i %s | ./bin/taffy stats -a > ./tests/stats.piped.a.txt", maf);
    CuAssertIntEquals(tc, 0, rc);
    CuAssertIntEquals(tc, 0, st_system("diff ./tests/stats.direct.a.txt ./tests/stats.piped.a.txt"));
    st_system("rm -f ./tests/stats.direct.a.txt ./tests/stats.piped.a.txt");

    // -b sequence intervals (previously rejected on MAF input)
    rc = st_system("./bin/taffy stats -b -i %s > ./tests/stats.direct.b.txt", maf);
    CuAssertIntEquals(tc, 0, rc);
    rc = st_system("./bin/taffy view -i %s | ./bin/taffy stats -b > ./tests/stats.piped.b.txt", maf);
    CuAssertIntEquals(tc, 0, rc);
    CuAssertIntEquals(tc, 0, st_system("diff ./tests/stats.direct.b.txt ./tests/stats.piped.b.txt"));
    st_system("rm -f ./tests/stats.direct.b.txt ./tests/stats.piped.b.txt");
}

CuSuite *stats_test_suite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stats_maf_input);
    return suite;
}

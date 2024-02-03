/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"

CuSuite* maf_test_suite(void);
CuSuite* taf_test_suite(void);
CuSuite* ond_test_suite(void);
CuSuite* normalize_test_suite(void);
CuSuite* view_test_suite(void);
CuSuite* sort_test_suite(void);

static int allTests(void) {
    CuString *output = CuStringNew();
    CuSuite* suite = CuSuiteNew();
    CuSuiteAddSuite(suite, maf_test_suite());
    CuSuiteAddSuite(suite, taf_test_suite());
    CuSuiteAddSuite(suite, ond_test_suite());
    CuSuiteAddSuite(suite, normalize_test_suite());
    CuSuiteAddSuite(suite, view_test_suite());
    CuSuiteAddSuite(suite, sort_test_suite());

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
    bool failed_tests = suite->failCount > 0;
    CuSuiteDelete(suite);
    CuStringDelete(output);
    return failed_tests;
}

int main(int argc, char *argv[]) {
    if(argc == 2) {
        st_setLogLevelFromString(argv[1]);
    }
    int i = allTests();
    //while(1);
    return i;
}

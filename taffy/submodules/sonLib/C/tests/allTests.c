/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLibGlobalsTest.h"

CuSuite* sonLib_ETreeTestSuite(void);
CuSuite* sonLib_stStringTestSuite(void);
CuSuite* sonLib_stHashTestSuite(void);
CuSuite* sonLib_stSetTestSuite(void);
CuSuite* sonLib_stSortedSetTestSuite(void);
CuSuite* sonLib_stListTestSuite(void);
CuSuite* sonLib_stCommonTestSuite(void);
CuSuite* sonLib_stIntTuplesTestSuite(void);
CuSuite* sonLib_stDoubleTuplesTestSuite(void);
CuSuite* sonLib_stExceptTestSuite(void);
CuSuite* sonLib_stRandomTestSuite(void);
CuSuite* sonLib_stCompressionTestSuite(void);
CuSuite* sonLibFileTestSuite(void);
CuSuite* stCacheSuite(void);
CuSuite* stPosetAlignmentTestSuite(void);
CuSuite* sonLibGraphTestSuite(void);
CuSuite* sonLib_stConnectivityTestSuite(void);
CuSuite* sonLib_stTreapTestSuite(void);
CuSuite* sonLib_stEulerTestSuite(void);
CuSuite* sonLib_stEdgeContainerTestSuite(void);
CuSuite* sonLib_stMatrixTestSuite(void);
CuSuite* sonLib_stPhylogenyTestSuite(void);
CuSuite* sonLib_stThreadPoolTestSuite(void);
CuSuite* sonLib_stUnionFindTestSuite(void);

int sonLibRunAllTests(void) {
    CuString *output = CuStringNew();
    CuSuite* suite = CuSuiteNew();
    CuSuiteAddSuite(suite, sonLib_stTreapTestSuite());
    CuSuiteAddSuite(suite, sonLib_stEdgeContainerTestSuite());
    CuSuiteAddSuite(suite, sonLib_stEulerTestSuite());
    CuSuiteAddSuite(suite, sonLib_stConnectivityTestSuite());
    CuSuiteAddSuite(suite, sonLib_stThreadPoolTestSuite());
    CuSuiteAddSuite(suite, sonLib_stPhylogenyTestSuite());
    CuSuiteAddSuite(suite, sonLib_stMatrixTestSuite());
    CuSuiteAddSuite(suite, sonLibGraphTestSuite());
    CuSuiteAddSuite(suite, stPosetAlignmentTestSuite());
    CuSuiteAddSuite(suite, sonLib_stCommonTestSuite());
    CuSuiteAddSuite(suite, sonLib_ETreeTestSuite());
    CuSuiteAddSuite(suite, sonLib_stStringTestSuite());
    CuSuiteAddSuite(suite, sonLib_stIntTuplesTestSuite());
    CuSuiteAddSuite(suite, sonLib_stDoubleTuplesTestSuite());
    CuSuiteAddSuite(suite, sonLib_stHashTestSuite());
    CuSuiteAddSuite(suite, sonLib_stSetTestSuite());
    CuSuiteAddSuite(suite, sonLib_stListTestSuite());
    CuSuiteAddSuite(suite, sonLib_stSortedSetTestSuite());
    CuSuiteAddSuite(suite, sonLib_stExceptTestSuite());
    CuSuiteAddSuite(suite, sonLib_stRandomTestSuite());
    CuSuiteAddSuite(suite, sonLib_stCompressionTestSuite());
    CuSuiteAddSuite(suite, sonLibFileTestSuite());
    CuSuiteAddSuite(suite, stCacheSuite());
    CuSuiteAddSuite(suite, sonLib_stUnionFindTestSuite());
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
    CuStringDelete(output);
    int status = suite->failCount > 0;
    CuSuiteDelete(suite);
    return status;
}
int main(int argc, char *argv[]) {
    if(argc == 2) {
        st_setLogLevelFromString(argv[1]);
    }
    return sonLibRunAllTests();
}

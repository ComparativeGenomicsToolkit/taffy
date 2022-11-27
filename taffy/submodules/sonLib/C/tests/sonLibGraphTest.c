/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLibGlobalsTest.h"

static const int64_t vertexNo = 7;
static stGraph *g = NULL; //Graph has vertices {0, 1, 2, 3, 4, 5, 6}

static double e12W = 0.5;
static double e13W = 0.3;
static double e14W = 0.2;
static double e24W = 0.1;
static double e34W = 0.2;
static double e35W = 0.4;
static double e05W = 0.2;

static void teardown() {
    if (g != NULL) {
        stGraph_destruct(g);
        g = NULL;
    }
}

static void setup() {
    teardown();
    g = stGraph_construct(vertexNo);
    stGraph_addEdge(g, 1, 2, e12W);
    stGraph_addEdge(g, 1, 3, e13W);
    stGraph_addEdge(g, 1, 4, e14W);
    stGraph_addEdge(g, 2, 4, e24W);
    stGraph_addEdge(g, 3, 4, e34W);
    stGraph_addEdge(g, 3, 5, e35W);
    stGraph_addEdge(g, 0, 5, e05W);
}

static void test_stGraph(CuTest *testCase) {
    setup();
    CuAssertIntEquals(testCase, stGraph_cardinality(g), vertexNo);
    stEdge *e = stGraph_getEdges(g, 0);
    CuAssertIntEquals(testCase, 5, stEdge_to(e));
    CuAssertDblEquals(testCase, e05W, stEdge_weight(e), 0.0);
    CuAssertPtrEquals(testCase, NULL, stEdge_nextEdge(e));
    e = stGraph_getEdges(g, 5);
    CuAssertIntEquals(testCase, 0, stEdge_to(e));
    CuAssertDblEquals(testCase, e05W, stEdge_weight(e), 0.0);
    e = stEdge_nextEdge(e);
    CuAssertIntEquals(testCase, 3, stEdge_to(e));
    CuAssertDblEquals(testCase, e35W, stEdge_weight(e), 0.0);
    CuAssertPtrEquals(testCase, NULL, stEdge_nextEdge(e));
    CuAssertPtrEquals(testCase, NULL, stGraph_getEdges(g, 6));
    teardown();
}

static void test_stGraph_shortestPaths(CuTest *testCase) {
    setup();
    double *dist = stGraph_shortestPaths(g, 1);
    double exDist[] = { 0.9, 0.0, 0.3, 0.3, 0.2, 0.7, INT64_MAX };
    for(int64_t v=0; v<vertexNo; v++) {
        CuAssertDblEquals(testCase, exDist[v], dist[v], 0.00001);
    }
    teardown();
    free(dist);
}

CuSuite* sonLibGraphTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stGraph);
    SUITE_ADD_TEST(suite, test_stGraph_shortestPaths);

    return suite;
}

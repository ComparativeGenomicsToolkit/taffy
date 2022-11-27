#include "CuTest.h"
#include "sonLib.h"

// Simple static test.
static void stUnionFind_staticTest(CuTest *testCase) {
    stUnionFind *unionFind = stUnionFind_construct();
    stUnionFind_add(unionFind, (void *) 1);
    stUnionFind_add(unionFind, (void *) 2);
    stUnionFind_add(unionFind, (void *) 3);
    stUnionFind_add(unionFind, (void *) 4);
    stUnionFind_add(unionFind, (void *) 5);
    stUnionFind_add(unionFind, (void *) 6);

    CuAssertTrue(testCase, 
                 stUnionFind_find(unionFind, (void *) 1) != stUnionFind_find(unionFind, (void *) 3));
    stUnionFind_union(unionFind, (void *) 1, (void *) 3);
    CuAssertTrue(testCase, 
                 stUnionFind_find(unionFind, (void *) 1) == stUnionFind_find(unionFind, (void *) 3));

    CuAssertTrue(testCase, 
                 stUnionFind_find(unionFind, (void *) 2) != stUnionFind_find(unionFind, (void *) 5));
    stUnionFind_union(unionFind, (void *) 2, (void *) 5);
    CuAssertTrue(testCase, 
                 stUnionFind_find(unionFind, (void *) 2) == stUnionFind_find(unionFind, (void *) 5));
    stUnionFind_union(unionFind, (void *) 4, (void *) 1);
    CuAssertTrue(testCase, 
                 stUnionFind_find(unionFind, (void *) 3) == stUnionFind_find(unionFind, (void *) 4));
    stUnionFind_union(unionFind, (void *) 1, (void *) 3);
    stUnionFind_union(unionFind, (void *) 1, (void *) 5);
    stUnionFind_union(unionFind, (void *) 3, (void *) 5);

    stUnionFindIt *it = stUnionFind_getIterator(unionFind);
    stSet *setOne = stUnionFindIt_getNext(it);
    stSet *setTwo = stUnionFindIt_getNext(it);
    CuAssertTrue(testCase, stUnionFindIt_getNext(it) == NULL);
    CuAssertTrue(testCase, stSet_size(setOne) == 5 || stSet_size(setTwo) == 5);
    CuAssertTrue(testCase, stSet_size(setOne) == 1 || stSet_size(setTwo) == 1);
    stUnionFind_destructIterator(it);
    stUnionFind_destruct(unionFind);
}

// Compare stUnionFind and the already battle-tested stConnectivity to
// ensure correctness.
static void stUnionFind_randomTest(CuTest *testCase) {
    for (int64_t iteration = 0; iteration < 100; iteration++) {
        stUnionFind *unionFind = stUnionFind_construct();
        stConnectivity *connectivity = stConnectivity_construct();
        int64_t numNodes = 2000;
        for (int64_t i = 1; i <= numNodes; i++) {
            stUnionFind_add(unionFind, (void *) i);
            stConnectivity_addNode(connectivity, (void *) i);
        }

        int64_t numTests = numNodes * 2;
        while (numTests > 0) {
            int64_t i = st_randomInt64(1, numNodes + 1);
            int64_t j = st_randomInt64(1, numNodes + 1);
            bool connectedInUnionFind = stUnionFind_find(unionFind, (void *) i) == stUnionFind_find(unionFind, (void *) j);
            bool connectedInConnectivity = stConnectivity_connected(connectivity, (void *) i, (void *) j);
            CuAssertIntEquals(testCase, connectedInConnectivity, connectedInUnionFind);

            stUnionFind_union(unionFind, (void *) i, (void *) j);
            if (i != j) {
                // stConnectivity requires i != j. Not adding the
                // self-loop to stConnectivity is OK when i == j,
                // because obviously self-loops cannot change the
                // connected components.
                stConnectivity_addEdge(connectivity, (void *) i, (void *) j);
            }

            CuAssertTrue(testCase, stUnionFind_find(unionFind, (void *) i) == stUnionFind_find(unionFind, (void *) j));
            numTests--;
        }

        // Check that the total number of components between the two are equal.
        int64_t numConnectivityComponents = 0;
        stConnectedComponentIterator *connIt = stConnectivity_getConnectedComponentIterator(connectivity);
        stConnectedComponent *connComponent;
        while ((connComponent = stConnectedComponentIterator_getNext(connIt)) != NULL) {
            numConnectivityComponents++;
        }
        stConnectedComponentIterator_destruct(connIt);

        int64_t numUnionFindComponents = 0;
        stUnionFindIt *unionFindIt = stUnionFind_getIterator(unionFind);
        stSet *unionFindComponent;
        while ((unionFindComponent = stUnionFindIt_getNext(unionFindIt)) != NULL) {
            numUnionFindComponents++;
        }
        stUnionFind_destructIterator(unionFindIt);
        CuAssertIntEquals(testCase, numConnectivityComponents, numUnionFindComponents);

        stUnionFind_destruct(unionFind);
        stConnectivity_destruct(connectivity);
    }
}

CuSuite *sonLib_stUnionFindTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, stUnionFind_staticTest);
    SUITE_ADD_TEST(suite, stUnionFind_randomTest);
    return suite;
}

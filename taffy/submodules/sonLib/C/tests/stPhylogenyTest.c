#include <stdlib.h>
#include <math.h>
#include "CuTest.h"
#include "sonLib.h"
#include "stPhylogeny.h"

static stIndexedTreeInfo *getIndex(stTree *node) {
    assert(node != NULL);
    stPhylogenyInfo *info = stTree_getClientData(node);
    assert(info != NULL);
    stIndexedTreeInfo *index = info->index;
    assert(index != NULL);
    return index;
}

static void testSimpleNeighborJoin(CuTest *testCase) {
    int64_t i, j;
    // Simple static test, with index 1 very far away from everything,
    // 0 and 3 very close, and 2 closer to (0,3) than to 1.
    double distances[4][4] = {{0.0, 9.0, 3.0, 0.1},
                              {9.0, 0.0, 6.0, 8.9},
                              {3.0, 6.0, 0.0, 3.0},
                              {0.1, 8.9, 3.0, 0.0}};
    stMatrix * distanceStMatrix = stMatrix_construct(4, 4);
    for(i = 0; i < 4; i++) {
        for(j = 0; j < 4; j++) {
            *stMatrix_getCell(distanceStMatrix, i, j) = distances[i][j];
        }
    }
    stTree *tree = stPhylogeny_neighborJoin(distanceStMatrix, NULL);
    stIndexedTreeInfo *info = getIndex(tree);

    // Check that leavesBelow is correct for the root (every leaf is
    // below the root)
    for(i = 0; i < 4; i++) {
        CuAssertIntEquals(testCase, info->leavesBelow[i], 1);
    }

    // We don't want to check the topology, since the root is
    // arbitrary (though currently halfway along the longest branch.)
    // So instead we check that the distances make sense. 0,3 should be the
    // closest together, etc.
    CuAssertTrue(testCase, stPhylogeny_distanceBetweenLeaves(tree, 0, 3) < stPhylogeny_distanceBetweenLeaves(tree, 0, 2));
    CuAssertTrue(testCase, stPhylogeny_distanceBetweenLeaves(tree, 0, 3) < stPhylogeny_distanceBetweenLeaves(tree, 0, 1));
    CuAssertTrue(testCase, stPhylogeny_distanceBetweenLeaves(tree, 0, 3) < stPhylogeny_distanceBetweenLeaves(tree, 3, 2));
    CuAssertTrue(testCase, stPhylogeny_distanceBetweenLeaves(tree, 0, 3) < stPhylogeny_distanceBetweenLeaves(tree, 3, 1));
    CuAssertTrue(testCase, stPhylogeny_distanceBetweenLeaves(tree, 0, 2) < stPhylogeny_distanceBetweenLeaves(tree, 0, 1));
    CuAssertTrue(testCase, stPhylogeny_distanceBetweenLeaves(tree, 0, 2) < stPhylogeny_distanceBetweenLeaves(tree, 2, 1));

    stMatrix_destruct(distanceStMatrix);
    stPhylogenyInfo_destructOnTree(tree);
    stTree_destruct(tree);
}

// Score a tree against a few simple "bootstrapped" trees to check
// that identical partitions are counted correctly.
void testSimpleBootstrapPartitionScoring(CuTest *testCase)
{
    stList *bootstraps = stList_construct();
    stTree *tree = stTree_parseNewickString("(((0,1),(2,3)),4);");
    stTree *bootstrap, *result;
    stIndexedTreeInfo *index;
    stPhylogeny_addStIndexedTreeInfo(tree);

    // An identical tree should give support to all partitions.
    bootstrap = stTree_parseNewickString("(4,((1,0),(3,2)));");
    stPhylogeny_addStIndexedTreeInfo(bootstrap);
    result = stPhylogeny_scoreFromBootstrap(tree, bootstrap);
    index = getIndex(stPhylogeny_getMRCA(result, 4, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    CuAssertTrue(testCase, index->bootstrapSupport == 1.0);
    index = getIndex(stPhylogeny_getMRCA(result, 3, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    index = getIndex(stPhylogeny_getMRCA(result, 1, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    index = getIndex(stPhylogeny_getMRCA(result, 3, 2));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    stPhylogenyInfo_destructOnTree(result);
    stTree_destruct(result);
    stList_append(bootstraps, bootstrap);

    // Swapping 0 for 4
    bootstrap = stTree_parseNewickString("(0,((1,4),(3,2)));");
    stPhylogeny_addStIndexedTreeInfo(bootstrap);
    result = stPhylogeny_scoreFromBootstrap(tree, bootstrap);
    index = getIndex(stPhylogeny_getMRCA(result, 4, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    CuAssertTrue(testCase, index->bootstrapSupport == 1.0);
    index = getIndex(stPhylogeny_getMRCA(result, 3, 0));
    CuAssertTrue(testCase, index->numBootstraps == 0);
    index = getIndex(stPhylogeny_getMRCA(result, 1, 0));
    CuAssertTrue(testCase, index->numBootstraps == 0);
    index = getIndex(stPhylogeny_getMRCA(result, 3, 2));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    stPhylogenyInfo_destructOnTree(result);
    stTree_destruct(result);
    stList_append(bootstraps, bootstrap);

    // Swapping 0 for 3
    bootstrap = stTree_parseNewickString("(4,((1,3),(0,2)));");
    stPhylogeny_addStIndexedTreeInfo(bootstrap);
    result = stPhylogeny_scoreFromBootstrap(tree, bootstrap);
    index = getIndex(stPhylogeny_getMRCA(result, 4, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    CuAssertTrue(testCase, index->bootstrapSupport == 1.0);
    index = getIndex(stPhylogeny_getMRCA(result, 3, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    index = getIndex(stPhylogeny_getMRCA(result, 1, 0));
    CuAssertTrue(testCase, index->numBootstraps == 0);
    index = getIndex(stPhylogeny_getMRCA(result, 3, 2));
    CuAssertTrue(testCase, index->numBootstraps == 0);
    stPhylogenyInfo_destructOnTree(result);
    stTree_destruct(result);
    stList_append(bootstraps, bootstrap);

    // Test all 3 together
    result = stPhylogeny_scoreFromBootstraps(tree, bootstraps);
    index = getIndex(stPhylogeny_getMRCA(result, 4, 0));
    CuAssertTrue(testCase, index->numBootstraps == 3);
    index = getIndex(stPhylogeny_getMRCA(result, 3, 0));
    CuAssertTrue(testCase, index->numBootstraps == 2);
    index = getIndex(stPhylogeny_getMRCA(result, 1, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    index = getIndex(stPhylogeny_getMRCA(result, 3, 2));
    CuAssertTrue(testCase, index->numBootstraps == 2);

    // Clean up
    stPhylogenyInfo_destructOnTree(tree);
    stTree_destruct(tree);
    stPhylogenyInfo_destructOnTree(result);
    stTree_destruct(result);
    for(int64_t i = 0; i < stList_length(bootstraps); i++) {
        stPhylogenyInfo_destructOnTree(stList_get(bootstraps, i));
        stTree_destruct(stList_get(bootstraps, i));
    }
    stList_destruct(bootstraps);
}

// Set recon info for a tree which already has its stPhylogenyInfo
// set.
static void setReconInfoManually(stTree *tree, const char *nodeLabel,
                                 stTree *species, stReconciliationEvent event) {
    stTree *node = stTree_findChild(tree, nodeLabel);
    assert(node != NULL);
    stPhylogenyInfo *info = stTree_getClientData(node);
    assert(info != NULL);
    stReconciliationInfo *recon = info->recon;
    if (recon == NULL) {
        recon = st_malloc(sizeof(stReconciliationInfo));
        info->recon = recon;
    }
    recon->species = species;
    recon->event = event;
}

// Score a tree against a few simple "bootstrapped" trees to check
// that identical partition with identical reconciliations are counted
// correctly.
void testSimpleBootstrapReconciliationScoring(CuTest *testCase)
{
    stList *bootstraps = stList_construct();
    stTree *tree = stTree_parseNewickString("(((0,1)a,(2,3)b)c,4)d;");
    stTree *bootstrap, *result;
    stIndexedTreeInfo *index;
    stPhylogeny_addStIndexedTreeInfo(tree);
    // Add dummy recon info

    // Dummy species that don't actually point anywhere. A bit nasty,
    // but it would be a bug to dereference these anyway.
    stTree *s0 = (stTree *) 0;
    stTree *s1 = (stTree *) 1;
    stTree *s2 = (stTree *) 2;
    stTree *s3 = (stTree *) 3;
    stTree *s4 = (stTree *) 4;
    stTree *s5 = (stTree *) 5;
    setReconInfoManually(tree, "0", s0, LEAF);
    setReconInfoManually(tree, "1", s1, LEAF);
    setReconInfoManually(tree, "2", s0, LEAF);
    setReconInfoManually(tree, "3", s1, LEAF);
    setReconInfoManually(tree, "4", s4, LEAF);
    setReconInfoManually(tree, "a", s2, SPECIATION);
    setReconInfoManually(tree, "b", s2, SPECIATION);
    setReconInfoManually(tree, "c", s2, DUPLICATION);
    setReconInfoManually(tree, "d", s3, SPECIATION);

    // An identical tree with identical reconciliation should give
    // support to all partitions.
    bootstrap = stTree_parseNewickString("(4,((1,0)a,(3,2)b)c)d;");
    stPhylogeny_addStIndexedTreeInfo(bootstrap);
    // Add dummy recon info identical to above
    setReconInfoManually(bootstrap, "0", s0, LEAF);
    setReconInfoManually(bootstrap, "1", s1, LEAF);
    setReconInfoManually(bootstrap, "2", s0, LEAF);
    setReconInfoManually(bootstrap, "3", s1, LEAF);
    setReconInfoManually(bootstrap, "4", s4, LEAF);
    setReconInfoManually(bootstrap, "a", s2, SPECIATION);
    setReconInfoManually(bootstrap, "b", s2, SPECIATION);
    setReconInfoManually(bootstrap, "c", s2, DUPLICATION);
    setReconInfoManually(bootstrap, "d", s3, SPECIATION);
    result = stPhylogeny_scoreReconciliationFromBootstrap(tree, bootstrap);
    // d
    index = getIndex(stPhylogeny_getMRCA(result, 4, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    CuAssertTrue(testCase, index->bootstrapSupport == 1.0);
    // c
    index = getIndex(stPhylogeny_getMRCA(result, 3, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    // a
    index = getIndex(stPhylogeny_getMRCA(result, 1, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    // b
    index = getIndex(stPhylogeny_getMRCA(result, 3, 2));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    stPhylogenyInfo_destructOnTree(result);
    stTree_destruct(result);
    stList_append(bootstraps, bootstrap);

    // Swapping 0 for 4 -- the reconciliation support level should
    // decrease, even though the reconciliation is identical, because
    // it relies on the partitions being the same
    bootstrap = stTree_parseNewickString("(0,((1,4)a,(3,2)b)c)d;");
    stPhylogeny_addStIndexedTreeInfo(bootstrap);
    // Add dummy recon info identical to above
    setReconInfoManually(bootstrap, "0", s0, LEAF);
    setReconInfoManually(bootstrap, "1", s1, LEAF);
    setReconInfoManually(bootstrap, "2", s0, LEAF);
    setReconInfoManually(bootstrap, "3", s1, LEAF);
    setReconInfoManually(bootstrap, "4", s4, LEAF);
    setReconInfoManually(bootstrap, "a", s2, SPECIATION);
    setReconInfoManually(bootstrap, "b", s2, SPECIATION);
    setReconInfoManually(bootstrap, "c", s2, DUPLICATION);
    setReconInfoManually(bootstrap, "d", s3, SPECIATION);
    result = stPhylogeny_scoreReconciliationFromBootstrap(tree, bootstrap);
    // d
    index = getIndex(stPhylogeny_getMRCA(result, 4, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    CuAssertTrue(testCase, index->bootstrapSupport == 1.0);
    // c
    index = getIndex(stPhylogeny_getMRCA(result, 3, 0));
    CuAssertTrue(testCase, index->numBootstraps == 0);
    // a
    index = getIndex(stPhylogeny_getMRCA(result, 1, 0));
    CuAssertTrue(testCase, index->numBootstraps == 0);
    // b
    index = getIndex(stPhylogeny_getMRCA(result, 3, 2));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    // NB: the reconciliation bootstrap support level for leaves is
    // not always 1!
    // 0
    index = getIndex(stPhylogeny_getLeafByIndex(result, 0));
    CuAssertTrue(testCase, index->numBootstraps == 0);
    // 4
    index = getIndex(stPhylogeny_getLeafByIndex(result, 4));
    CuAssertTrue(testCase, index->numBootstraps == 0);
    // 2
    index = getIndex(stPhylogeny_getLeafByIndex(result, 2));
    CuAssertTrue(testCase, index->numBootstraps == 1);

    stPhylogenyInfo_destructOnTree(result);
    stTree_destruct(result);
    stList_append(bootstraps, bootstrap);

    // Different reconciliation but same topology. Test that different
    // event (duplication, speciation) labels affect the support value
    // as well.
    bootstrap = stTree_parseNewickString("(4,((1,0)a,(3,2)b)c)d;");
    stPhylogeny_addStIndexedTreeInfo(bootstrap);
    // Add dummy recon info, but this time, the event for a is a
    // duplication (impossible in MP reconciliation, but we should
    // support it anyway).
    setReconInfoManually(bootstrap, "0", s0, LEAF);
    setReconInfoManually(bootstrap, "1", s1, LEAF);
    setReconInfoManually(bootstrap, "2", s0, LEAF);
    setReconInfoManually(bootstrap, "3", s1, LEAF);
    setReconInfoManually(bootstrap, "4", s4, LEAF);
    setReconInfoManually(bootstrap, "a", s2, DUPLICATION);
    setReconInfoManually(bootstrap, "b", s2, SPECIATION);
    setReconInfoManually(bootstrap, "c", s2, DUPLICATION);
    setReconInfoManually(bootstrap, "d", s5, SPECIATION);
    result = stPhylogeny_scoreReconciliationFromBootstrap(tree, bootstrap);
    // d
    index = getIndex(stPhylogeny_getMRCA(result, 4, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    CuAssertTrue(testCase, index->bootstrapSupport == 1.0);
    // c
    index = getIndex(stPhylogeny_getMRCA(result, 3, 0));
    CuAssertTrue(testCase, index->numBootstraps == 0);
    // a
    index = getIndex(stPhylogeny_getMRCA(result, 1, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    // b
    index = getIndex(stPhylogeny_getMRCA(result, 3, 2));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    // 1
    index = getIndex(stPhylogeny_getLeafByIndex(result, 1));
    CuAssertTrue(testCase, index->numBootstraps == 0);
    // 2
    index = getIndex(stPhylogeny_getLeafByIndex(result, 2));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    stPhylogenyInfo_destructOnTree(result);
    stTree_destruct(result);
    stList_append(bootstraps, bootstrap);

    // Test all 3 together
    result = stPhylogeny_scoreReconciliationFromBootstraps(tree, bootstraps);
    // d
    index = getIndex(stPhylogeny_getMRCA(result, 4, 0));
    CuAssertTrue(testCase, index->numBootstraps == 3);
    CuAssertTrue(testCase, index->bootstrapSupport == 1.0);
    // c
    index = getIndex(stPhylogeny_getMRCA(result, 3, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    // a
    index = getIndex(stPhylogeny_getMRCA(result, 1, 0));
    CuAssertTrue(testCase, index->numBootstraps == 2);
    // b
    index = getIndex(stPhylogeny_getMRCA(result, 3, 2));
    CuAssertTrue(testCase, index->numBootstraps == 3);
    // 4
    index = getIndex(stPhylogeny_getLeafByIndex(result, 4));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    // 3
    index = getIndex(stPhylogeny_getLeafByIndex(result, 3));
    CuAssertTrue(testCase, index->numBootstraps == 3);
    // 2
    index = getIndex(stPhylogeny_getLeafByIndex(result, 2));
    CuAssertTrue(testCase, index->numBootstraps == 3);
    // 1
    index = getIndex(stPhylogeny_getLeafByIndex(result, 1));
    CuAssertTrue(testCase, index->numBootstraps == 2);
    // 0
    index = getIndex(stPhylogeny_getLeafByIndex(result, 0));
    CuAssertTrue(testCase, index->numBootstraps == 1);
    stPhylogenyInfo_destructOnTree(result);
    stTree_destruct(result);

    // Clean up
    stPhylogenyInfo_destructOnTree(tree);
    stTree_destruct(tree);
    for(int64_t i = 0; i < stList_length(bootstraps); i++) {
        stPhylogenyInfo_destructOnTree(stList_get(bootstraps, i));
        stTree_destruct(stList_get(bootstraps, i));
    }
    stList_destruct(bootstraps);
}

static stMatrix *getRandomDistanceMatrix(int64_t size) {
    int64_t i, j;
    stMatrix *ret = stMatrix_construct(size, size);
    for (i = 0; i < size; i++) {
        for (j = 0; j <= i; j++) {
            double *cell = stMatrix_getCell(ret, i, j);
            double val = st_random();
            *cell = val;
            // Make sure the other half of the matrix is identical
            cell = stMatrix_getCell(ret, size - i - 1, size - j - 1);
            *cell = val;
        }
    }
    return ret;
}

// Test an assert on all nodes of a tree.
static void testOnTree(CuTest *testCase, stTree *tree,
                       void (*assertFn)(stTree *, CuTest *)) {
    int64_t i;
    assertFn(tree, testCase);
    for (i = 0; i < stTree_getChildNumber(tree); i++) {
        testOnTree(testCase, stTree_getChild(tree, i), assertFn);
    }
}

static void checkLeavesBelow(stTree *tree, CuTest *testCase) {
    int64_t i;
    stIndexedTreeInfo *index = getIndex(tree);
    for (i = 0; i < index->totalNumLeaves; i++) {
        char *label = stString_print("%" PRIi64, i);
        if (index->leavesBelow[i]) {
            // leavesBelow says it has leaf i under it -- check
            // that it actually does
            CuAssertTrue(testCase, stTree_findChild(tree, label) != NULL ||
                         index->matrixIndex == i);
        } else {
            // leavesBelow says leaf i is not under this node, check
            // that it's not
            CuAssertTrue(testCase, stTree_findChild(tree, label) == NULL);
        }
        free(label);
    }
}

// check the stPhylogeny_getLeafByIndex function on a tree.
static void checkGetLeafByIndex(stTree *tree, CuTest *testCase) {
    int64_t i;
    stIndexedTreeInfo *index = getIndex(tree);
    for (i = 0; i < index->totalNumLeaves; i++) {
        if (index->leavesBelow[i]) {
            stTree *leaf = stPhylogeny_getLeafByIndex(tree, i);
            stIndexedTreeInfo *leafIndex = getIndex(leaf);
            CuAssertTrue(testCase, leafIndex->matrixIndex == i);
        } else {
            CuAssertTrue(testCase, stPhylogeny_getLeafByIndex(tree, i) == NULL);
        }
    }
}

// Sanity-check the assorted distance functions
static void checkDistanceFunctions(stTree *tree, CuTest *testCase) {
    int64_t i, j;
    stIndexedTreeInfo *index = getIndex(tree);
    for (i = 0; i < index->totalNumLeaves; i++) {
        for (j = 0; j < i; j++) {
            if (index->leavesBelow[i] && index->leavesBelow[j]) {
                // Check that the distance from leaf i to leaf j is the same as
                // the distance from j to i
                CuAssertDblEquals(testCase, stPhylogeny_distanceBetweenLeaves(tree, i, j), stPhylogeny_distanceBetweenLeaves(tree, j, i), 0.1);
                // Check that the distance using distanceBetweenNodes
                // is the same as using distanceBetweenLeaves
                stTree *iNode = stPhylogeny_getLeafByIndex(tree, i);
                stTree *jNode = stPhylogeny_getLeafByIndex(tree, j);
                CuAssertDblEquals(testCase, stPhylogeny_distanceBetweenLeaves(tree, i, j), stPhylogeny_distanceBetweenNodes(iNode, jNode), 0.1);
            }
        }
    }

    CuAssertDblEquals(testCase, stPhylogeny_distanceBetweenNodes(tree, tree), 0.0, 0.0);
}

// Test neighbor-joining on random distance matrices, and test several
// stPhylogeny functions on the resulting trees.
static void testRandomNeighborJoin(CuTest *testCase) {
    int64_t testNum;
    for (testNum = 0; testNum < 20; testNum++) {
        // Run neighbor-joining
        int64_t numLeaves = st_randomInt64(3, 300);
        stMatrix *matrix = getRandomDistanceMatrix(numLeaves);
        bool rootOnOutgroup = st_random() > 0.5;
        stList *outgroups = NULL;

        if(rootOnOutgroup) {
            // Create outgroup list
            outgroups = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
            int64_t numOutgroups = st_randomInt64(1, numLeaves);
            for(int64_t i = 0; i < numOutgroups; i++) {
                int64_t outgroup = st_randomInt64(0, numLeaves);
                stList_append(outgroups, stIntTuple_construct1(outgroup));
            }
        }

        stTree *tree = stPhylogeny_neighborJoin(matrix, outgroups);

        if(rootOnOutgroup) {
            // Check that the root is placed on one of the outgroup branches.
            bool outgroupBranchIsRoot = false;
            for(int64_t i = 0; i < stList_length(outgroups); i++) {
                int64_t outgroup = stIntTuple_get(stList_get(outgroups, i), 0);
                if(stTree_getParent(stTree_getParent(stPhylogeny_getLeafByIndex(tree, outgroup))) == NULL) {
                    outgroupBranchIsRoot = true;
                }
            }
            CuAssertTrue(testCase, outgroupBranchIsRoot);
            stList_destruct(outgroups);
        }

        // Check that leavesBelow is set correctly
        testOnTree(testCase, tree, checkLeavesBelow);

        // Check that stPhylogeny_getLeafIndex works
        testOnTree(testCase, tree, checkGetLeafByIndex);

        // Sanity-check the distance functions.
        checkDistanceFunctions(tree, testCase);

        // Clean up
        stMatrix_destruct(matrix);
        stPhylogenyInfo_destructOnTree(tree);
        stTree_destruct(tree);
    }
}

static int64_t numBootstraps; // Totally lazy, but enables
                              // checkPartitionSupport to see the
                              // number of bootstraps used

// Check that the partition-support values on a bootstrap-scored tree
// are all sensible.
static void checkPartitionSupport(stTree *tree, CuTest *testCase) {
    stIndexedTreeInfo *index = getIndex(tree);
    CuAssertTrue(testCase, index->numBootstraps <= numBootstraps);
    CuAssertTrue(testCase, index->numBootstraps >= 0);
    CuAssertTrue(testCase, index->bootstrapSupport >= 0.0);
    CuAssertTrue(testCase, index->bootstrapSupport <= 1.0);
}

// Run faux "bootstraps" on random distance matrices of the same size,
// and test several stPhylogeny functions on the resulting trees.
static void testRandomBootstraps(CuTest *testCase) {
    stList *bootstraps = stList_construct();
    int64_t matrixSize = st_randomInt64(3, 300);
    numBootstraps = st_randomInt64(1, 200);
    stMatrix *canonicalMatrix = getRandomDistanceMatrix(matrixSize);
    stTree *canonicalTree = stPhylogeny_neighborJoin(canonicalMatrix, NULL);
    for(int64_t i = 0; i < numBootstraps; i++) {
        // Not a real bootstrap matrix, but the indices are the same
        // so the bootstrap code still needs to work.
        stMatrix *bootstrapMatrix = getRandomDistanceMatrix(matrixSize);
        stTree *bootstrapTree = stPhylogeny_neighborJoin(bootstrapMatrix, NULL);
        stList_append(bootstraps, bootstrapTree);

        stMatrix_destruct(bootstrapMatrix);
    }

    // Run the code that scores partitions by bootstraps
    stTree *bootstrappedTree = stPhylogeny_scoreFromBootstraps(canonicalTree, bootstraps);

    // Check that leavesBelow is set correctly
    testOnTree(testCase, bootstrappedTree, checkLeavesBelow);

    // Check that stPhylogeny_getLeafIndex works
    testOnTree(testCase, bootstrappedTree, checkGetLeafByIndex);

    // Sanity-check the distance functions.
    checkDistanceFunctions(bootstrappedTree, testCase);

    // Root support should always be 1
    stIndexedTreeInfo *rootIndex = getIndex(bootstrappedTree);
    CuAssertTrue(testCase, rootIndex->bootstrapSupport == 1.0);
    // Check that the parition support numbers are all reasonable
    testOnTree(testCase, bootstrappedTree, checkPartitionSupport);

    // Clean up
    stMatrix_destruct(canonicalMatrix);
    stPhylogenyInfo_destructOnTree(canonicalTree);
    stTree_destruct(canonicalTree);
    stPhylogenyInfo_destructOnTree(bootstrappedTree);
    stTree_destruct(bootstrappedTree);
    for(int64_t i = 0; i < stList_length(bootstraps); i++) {
        stPhylogenyInfo_destructOnTree(stList_get(bootstraps, i));
        stTree_destruct(stList_get(bootstraps, i));
    }
    stList_destruct(bootstraps);
}

static double getJoinCost(stMatrix *matrix, stHash *speciesToIndex, stTree *tree, const char *label1, const char *label2) {
    stTree *node1;
    if (strcmp(stTree_getLabel(tree), label1) == 0) {
        node1 = tree;
    } else {
        node1 = stTree_findChild(tree, label1);
    }
    assert(node1 != NULL);
    stTree *node2;
    if (strcmp(stTree_getLabel(tree), label2) == 0) {
        node2 = tree;
    } else {
        node2 = stTree_findChild(tree, label2);
    }
    assert(node2 != NULL);
    stIntTuple *index1 = stHash_search(speciesToIndex, node1);
    assert(index1 != NULL);
    stIntTuple *index2 = stHash_search(speciesToIndex, node2);
    assert(index2 != NULL);
    int64_t i = stIntTuple_get(index1, 0);
    int64_t j = stIntTuple_get(index2, 0);
    assert(*stMatrix_getCell(matrix, i, j) == *stMatrix_getCell(matrix, i, j));
    return *stMatrix_getCell(matrix, i, j);
}

// Make sure the join costs agree with hand calculations on a very
// simple tree.
static void testSimpleJoinCosts(CuTest *testCase) {
    stTree *tree = stTree_parseNewickString("((A,B)C,E)D;");

    // check speciesToIndex is set correctly
    stHash *speciesToIndex = stHash_construct2(NULL, (void (*)(void *)) stIntTuple_destruct);
    stMatrix *matrix = stPhylogeny_computeJoinCosts(tree, speciesToIndex, 1.0, 0.0);
    CuAssertIntEquals(testCase, stTree_getNumNodes(tree), stHash_size(speciesToIndex));
    stSet *seen = stSet_construct3((uint64_t (*)(const void *)) stIntTuple_hashKey, (int (*)(const void *, const void *)) stIntTuple_equalsFn, NULL);
    // Make sure there are no duplicate indices and that they are all >= 0.
    stHashIterator *hashIt = stHash_getIterator(speciesToIndex);
    stTree *curNode;
    while ((curNode = stHash_getNext(hashIt)) != NULL) {
        stIntTuple *index = stHash_search(speciesToIndex, curNode);
        CuAssertTrue(testCase, index != NULL);
        CuAssertTrue(testCase, !stSet_search(seen, index));
        CuAssertTrue(testCase, stIntTuple_get(index, 0) >= 0);
    }
    stHash_destructIterator(hashIt);
    stMatrix_destruct(matrix);
    stHash_destruct(speciesToIndex);

    // check dups
    speciesToIndex = stHash_construct2(NULL, (void (*)(void *)) stIntTuple_destruct);
    matrix = stPhylogeny_computeJoinCosts(tree, speciesToIndex, 1.0, 0.0);
    CuAssertDblEquals(testCase, 1.0, getJoinCost(matrix, speciesToIndex, tree, "A", "A"), 0.01);
    CuAssertDblEquals(testCase, 0.0, getJoinCost(matrix, speciesToIndex, tree, "A", "B"), 0.01);
    CuAssertDblEquals(testCase, 0.0, getJoinCost(matrix, speciesToIndex, tree, "B", "A"), 0.01);
    CuAssertDblEquals(testCase, 0.0, getJoinCost(matrix, speciesToIndex, tree, "A", "E"), 0.01);
    CuAssertDblEquals(testCase, 1.0, getJoinCost(matrix, speciesToIndex, tree, "A", "D"), 0.01);
    CuAssertDblEquals(testCase, 1.0, getJoinCost(matrix, speciesToIndex, tree, "A", "C"), 0.01);
    CuAssertDblEquals(testCase, 1.0, getJoinCost(matrix, speciesToIndex, tree, "C", "A"), 0.01);
    CuAssertDblEquals(testCase, 1.0, getJoinCost(matrix, speciesToIndex, tree, "D", "D"), 0.01);
    stMatrix_destruct(matrix);
    stHash_destruct(speciesToIndex);

    // check losses
    speciesToIndex = stHash_construct2(NULL, (void (*)(void *)) stIntTuple_destruct);
    matrix = stPhylogeny_computeJoinCosts(tree, speciesToIndex, 0.0, 1.0);
    CuAssertDblEquals(testCase, 0.0, getJoinCost(matrix, speciesToIndex, tree, "A", "A"), 0.01);
    CuAssertDblEquals(testCase, 0.0, getJoinCost(matrix, speciesToIndex, tree, "A", "B"), 0.01);
    CuAssertDblEquals(testCase, 0.0, getJoinCost(matrix, speciesToIndex, tree, "B", "A"), 0.01);
    CuAssertDblEquals(testCase, 1.0, getJoinCost(matrix, speciesToIndex, tree, "A", "E"), 0.01);
    CuAssertDblEquals(testCase, 2.0, getJoinCost(matrix, speciesToIndex, tree, "A", "D"), 0.01);
    CuAssertDblEquals(testCase, 1.0, getJoinCost(matrix, speciesToIndex, tree, "A", "C"), 0.01);
    CuAssertDblEquals(testCase, 1.0, getJoinCost(matrix, speciesToIndex, tree, "C", "A"), 0.01);
    CuAssertDblEquals(testCase, 0.0, getJoinCost(matrix, speciesToIndex, tree, "D", "D"), 0.01);

    stMatrix_destruct(matrix);
    stHash_destruct(speciesToIndex);
    stSet_destruct(seen);
    stTree_destruct(tree);
}

// get a random binary tree. numLeaves will get set to the eventual
// number of leaves.
static stTree *getRandomBinaryTree(int64_t maxDepth, int64_t *numLeaves) {
    stTree *ret = stTree_construct();
    int64_t numChildren = st_random() > 0.5 ? 2 : 0;
    if (maxDepth != 0) {
        for (int64_t i = 0; i < numChildren; i++) {
            stTree_setParent(getRandomBinaryTree(maxDepth - 1, numLeaves), ret);
        }
    }
    if (stTree_getChildNumber(ret) == 0) {
        char *name = stString_print("%" PRIi64, (*numLeaves)++);
        stTree_setLabel(ret, name);
        free(name);
    }
    return ret;
}

// Test that the join costs are equal to what the reconciliation cost
// functions would suggest.
static void testJoinCosts_random(CuTest *testCase) {
    for (int64_t testNum = 0; testNum < 100; testNum++) {
        int64_t numSpeciesLeaves = 0;
        stTree *speciesTree = getRandomBinaryTree(st_randomInt64(2, 4),
                                                  &numSpeciesLeaves);
        int64_t numSpecies = stTree_getNumNodes(speciesTree);
        stHash *speciesToIndex = stHash_construct2(NULL, (void (*)(void *)) stIntTuple_destruct);
        double dupCost = st_random();
        double lossCost = st_random();
        stMatrix *joinCosts = stPhylogeny_computeJoinCosts(speciesTree,
                                                           speciesToIndex,
                                                           dupCost, lossCost);
        stHash *indexToSpecies = stHash_invert(speciesToIndex, (uint64_t (*)(const void *)) stIntTuple_hashKey, (int (*)(const void *, const void *)) stIntTuple_equalsFn, NULL, NULL);
        for (int64_t i = 0; i < numSpecies; i++) {
            for (int64_t j = 0; j < numSpecies; j++) {
                // Check that the tree joining these two nodes really
                // would get the same join cost according to the
                // reconciliation cost functions.
                stTree *geneTree = stTree_parseNewickString("(I,J)MRCA;");
                stTree *iNode = stTree_findChild(geneTree, "I");
                assert(iNode != NULL);
                stTree *jNode = stTree_findChild(geneTree, "J");
                assert(jNode != NULL);
                // Construct the initial reconciliation for the leaves
                // of the gene tree.
                stHash *leafToSpecies = stHash_construct();
                stIntTuple *query = stIntTuple_construct1(i);
                stTree *iSpecies = stHash_search(indexToSpecies, query);
                assert(iSpecies != NULL);
                stIntTuple_destruct(query);
                stHash_insert(leafToSpecies, iNode, iSpecies);
                query = stIntTuple_construct1(j);
                stTree *jSpecies = stHash_search(indexToSpecies, query);
                assert(jSpecies != NULL);
                stIntTuple_destruct(query);
                stHash_insert(leafToSpecies, jNode, jSpecies);
                // Reconcile the tree and calculate the join cost we expect.
                stPhylogeny_reconcileAtMostBinary(geneTree, leafToSpecies, false);
                int64_t dups = 0, losses = 0;
                stPhylogeny_reconciliationCostAtMostBinary(geneTree, &dups,
                                                           &losses);
                double expectedCost = dupCost * dups + lossCost * losses;
                double actualCost = *stMatrix_getCell(joinCosts, i, j);
                CuAssertDblEquals(testCase, *stMatrix_getCell(joinCosts, i, j),
                                  *stMatrix_getCell(joinCosts, i, j), 0.01);
                CuAssertDblEquals(testCase, expectedCost, actualCost, 0.01);
                stHash_destruct(leafToSpecies);
                stPhylogenyInfo_destructOnTree(geneTree);
                stTree_destruct(geneTree);
            }
        }
        stHash_destruct(indexToSpecies);
        stHash_destruct(speciesToIndex);
        stTree_destruct(speciesTree);
        stMatrix_destruct(joinCosts);
    }
}

static stMatrix *getRandomSimilarityMatrix(int64_t size, int64_t maxNumSimilarities, int64_t maxNumDifferences) {
    stMatrix *ret = stMatrix_construct(size, size);
    for (int64_t i = 0; i < size; i++) {
        for (int64_t j = i + 1; j < size; j++) {
            double similarities = st_random() * maxNumSimilarities;
            double differences = st_random() * maxNumDifferences;
            *stMatrix_getCell(ret, i, j) = similarities;
            *stMatrix_getCell(ret, j, i) = differences;
        }
    }
    return ret;
}

static stMatrix *getDistanceMatrixFromSimilarityMatrix(stMatrix *similarityMatrix) {
    assert(stMatrix_n(similarityMatrix) == stMatrix_m(similarityMatrix));
    stMatrix *ret = stMatrix_construct(stMatrix_n(similarityMatrix), stMatrix_m(similarityMatrix));
    for (int64_t i = 0; i < stMatrix_n(ret); i++) {
        for (int64_t j = i + 1; j < stMatrix_m(ret); j++) {
            double similarities = *stMatrix_getCell(similarityMatrix, i, j);
            double differences = *stMatrix_getCell(similarityMatrix, j, i);
            double count = similarities + differences;
            *stMatrix_getCell(ret, i, j) = (count != 0.0) ? differences / count : INT64_MAX;
            *stMatrix_getCell(ret, j, i) = (count != 0.0) ? differences / count : INT64_MAX;
        }
    }
    return ret;
}

static bool isTopologyEqual(stTree *tree1, stTree *tree2) {
    if (stTree_getChildNumber(tree1) != stTree_getChildNumber(tree2)) {
        return false;
    }

    for (int64_t i = 0; i < stTree_getChildNumber(tree1); i++) {
        stTree *child1 = stTree_getChild(tree1, i);
        // Find matching child in tree2.
        stIndexedTreeInfo *index1 = getIndex(child1);
        stTree *matchingChild = NULL;
        for (int64_t j = 0; j < stTree_getChildNumber(tree2); j++) {
            stTree *child2 = stTree_getChild(tree2, j);
            stIndexedTreeInfo *index2 = getIndex(child2);
            if (index1->totalNumLeaves != index2->totalNumLeaves) {
                return false;
            }
            if (!memcmp(index1->leavesBelow, index2->leavesBelow, index1->totalNumLeaves)) {
                // Found correct child
                matchingChild = child2;
                break;
            }
        }
        if (matchingChild == NULL) {
            return false;
        }
        if (!isTopologyEqual(child1, matchingChild)) {
            return false;
        }
    }
    return true;
}

static void testGuidedNeighborJoiningReducesToNeighborJoining(CuTest *testCase) {
    for (int64_t testNum = 0; testNum < 100; testNum++) {
        // Get a random species tree.
        int64_t numSpecies = 0;
        stTree *speciesTree = getRandomBinaryTree(st_randomInt64(2, 7), &numSpecies);
        if (numSpecies < 3) {
            stTree_destruct(speciesTree);
            continue;
        }

        // Join cost stuff
        stHash *speciesToIndex = stHash_construct2(NULL, (void (*)(void *)) stIntTuple_destruct);
        stMatrix *joinCosts = stPhylogeny_computeJoinCosts(speciesTree, speciesToIndex, 0.0, 0.0);

        stMatrix *similarityMatrix = getRandomSimilarityMatrix(numSpecies, 50, 50);
        stMatrix *distanceMatrix = getDistanceMatrixFromSimilarityMatrix(similarityMatrix);
        stTree *neighborJoiningTree = stPhylogeny_neighborJoin(distanceMatrix, NULL);

        stHash *matrixIndexToJoinCostIndex = stHash_construct3((uint64_t (*)(const void *)) stIntTuple_hashKey, (int (*)(const void *, const void *)) stIntTuple_equalsFn, (void (*)(void *)) stIntTuple_destruct, (void (*)(void *)) stIntTuple_destruct);
        // assign the matrix indices to be equal to the species indices.
        for (int64_t i = 0; i < numSpecies; i++) {
            stIntTuple *iTuple = stIntTuple_construct1(i);
            char *speciesName = stString_print("%" PRIi64, i);
            stIntTuple *joinCostIndex = stHash_search(speciesToIndex, stTree_findChild(speciesTree, speciesName));
            stHash_insert(matrixIndexToJoinCostIndex, iTuple, stIntTuple_construct1(stIntTuple_get(joinCostIndex, 0)));
            free(speciesName);
        }
        // Get the MRCA matrix
        int64_t **speciesMRCAMatrix = stPhylogeny_getMRCAMatrix(speciesTree, speciesToIndex);

        stTree *guidedNeighborJoiningTree = stPhylogeny_guidedNeighborJoining(distanceMatrix, similarityMatrix, joinCosts, matrixIndexToJoinCostIndex, speciesToIndex, speciesMRCAMatrix, speciesTree);

        for(int i = 0; i < stTree_getNumNodes(speciesTree); i++) {
          free(speciesMRCAMatrix[i]);
        }
        free(speciesMRCAMatrix);

        // Root both the trees the same way (above node 0)
        stPhylogenyInfo_destructOnTree(neighborJoiningTree);
        stPhylogenyInfo_destructOnTree(guidedNeighborJoiningTree);
        stTree *tmp = stTree_findChild(neighborJoiningTree, "0");
        stTree *rootedNeighborJoiningTree = stTree_reRoot(tmp, stTree_getBranchLength(tmp)/2);
        tmp = stTree_findChild(guidedNeighborJoiningTree, "0");
        stTree *rootedGuidedNeighborJoiningTree = stTree_reRoot(tmp, stTree_getBranchLength(tmp)/2);
        stPhylogeny_addStIndexedTreeInfo(rootedNeighborJoiningTree);
        stPhylogeny_addStIndexedTreeInfo(rootedGuidedNeighborJoiningTree);

        // Check that the topologies are the same in neighbor-joining
        // vs guided neighbor-joining (when they are rooted the same way)
        CuAssertTrue(testCase, isTopologyEqual(rootedNeighborJoiningTree, rootedGuidedNeighborJoiningTree));

        // Clean up.
        stPhylogenyInfo_destructOnTree(rootedGuidedNeighborJoiningTree);
        stTree_destruct(rootedGuidedNeighborJoiningTree);
        stPhylogenyInfo_destructOnTree(rootedNeighborJoiningTree);
        stTree_destruct(rootedNeighborJoiningTree);
        stTree_destruct(guidedNeighborJoiningTree);
        stTree_destruct(neighborJoiningTree);
        stTree_destruct(speciesTree);
        stMatrix_destruct(similarityMatrix);
        stMatrix_destruct(distanceMatrix);
        stMatrix_destruct(joinCosts);
        stHash_destruct(matrixIndexToJoinCostIndex);
        stHash_destruct(speciesToIndex);
    }
}

// Check that when join costs are ratcheted up to insane levels, the
// tree produced has minimal reconciliation cost.
static void testGuidedNeighborJoiningLowersReconCost(CuTest *testCase)
{
    for (int64_t testNum = 0; testNum < 100; testNum++) {
        // Get a random species tree.
        int64_t numSpecies = 0;
        stTree *speciesTree = getRandomBinaryTree(st_randomInt64(2, 7), &numSpecies);
        if (numSpecies < 3) {
            stTree_destruct(speciesTree);
            continue;
        }

        // Join cost stuff
        stHash *speciesToIndex = stHash_construct2(NULL, (void (*)(void *)) stIntTuple_destruct);
        stMatrix *joinCosts = stPhylogeny_computeJoinCosts(speciesTree, speciesToIndex, 100000.0, 100000.0);

        // Create a set of genes of size n*(num species), where each leaf
        // species is mapped to n genes. We'd like the guided neighbor
        // joining process to create a tree that can be reconciled so that
        // there are only n dups.
        int64_t numGenesPerSpecies = st_randomInt64(2, 5);
        stMatrix *similarityMatrix = getRandomSimilarityMatrix(numSpecies*numGenesPerSpecies, 50, 50);
        // assign the matrix indices i, 2i, 3i, ... to map to species index i.
        stHash *matrixIndexToJoinCostIndex = stHash_construct3((uint64_t (*)(const void *)) stIntTuple_hashKey, (int (*)(const void *, const void *)) stIntTuple_equalsFn, (void (*)(void *)) stIntTuple_destruct, (void (*)(void *)) stIntTuple_destruct);
        for (int64_t i = 0; i < numSpecies; i++) {
            char *speciesName = stString_print("%" PRIi64, i);
            stIntTuple *joinCostIndex = stHash_search(speciesToIndex, stTree_findChild(speciesTree, speciesName));
            for (int64_t j = i; j < numGenesPerSpecies * numSpecies; j += numSpecies) {
                stIntTuple *matrixIndex = stIntTuple_construct1(j);
                stHash_insert(matrixIndexToJoinCostIndex, matrixIndex, stIntTuple_construct1(stIntTuple_get(joinCostIndex, 0)));
            }
            free(speciesName);
        }

        // get MRCA matrix
        int64_t **speciesMRCAMatrix = stPhylogeny_getMRCAMatrix(speciesTree, speciesToIndex);

        // get distance matrix
        stMatrix *distanceMatrix = getDistanceMatrixFromSimilarityMatrix(similarityMatrix);

        stTree *guidedNeighborJoiningTree = stPhylogeny_guidedNeighborJoining(distanceMatrix, similarityMatrix, joinCosts, matrixIndexToJoinCostIndex, speciesToIndex, speciesMRCAMatrix, speciesTree);

        stMatrix_destruct(distanceMatrix);
        for(int i = 0; i < stTree_getNumNodes(speciesTree); i++) {
          free(speciesMRCAMatrix[i]);
        }
        free(speciesMRCAMatrix);

        // Build the gene->species mapping for the reconciliation algorithm.
        stHash *leafToSpecies = stHash_construct();
        for (int64_t i = 0; i < numSpecies; i++) {
            char *speciesName = stString_print("%" PRIi64, i);
            stTree *species = stTree_findChild(speciesTree, speciesName);
            for (int64_t j = i; j < numGenesPerSpecies * numSpecies; j += numSpecies) {
                char *geneName = stString_print("%" PRIi64, j);
                stTree *gene = stTree_findChild(guidedNeighborJoiningTree, geneName);
                stHash_insert(leafToSpecies, gene, species);
                free(geneName);
            }
            free(speciesName);
        }

        // Now check that the gene tree we get out of the algorithm has
        // minimal reconcilation cost due to the very high join costs (in
        // this case, the number of dups should be equal to the number of
        // "extra" genes per species).
        int64_t dups = 0, losses = 0;
        stPhylogeny_reconcileAtMostBinary(guidedNeighborJoiningTree,
                                          leafToSpecies, false);
        stPhylogeny_reconciliationCostAtMostBinary(guidedNeighborJoiningTree,
                                                   &dups, &losses);
        CuAssertTrue(testCase, losses == 0);
        CuAssertTrue(testCase, dups == numGenesPerSpecies - 1);

        // Clean up.
        stTree_destruct(speciesTree);
        stHash_destruct(matrixIndexToJoinCostIndex);
        stHash_destruct(speciesToIndex);
        stHash_destruct(leafToSpecies);
        stPhylogenyInfo_destructOnTree(guidedNeighborJoiningTree);
        stTree_destruct(guidedNeighborJoiningTree);
        stMatrix_destruct(similarityMatrix);
        stMatrix_destruct(joinCosts);
    }
}

// Globals for checkMinimalReconScore
stTree *globalSpeciesTree;
stHash *globalLeafToSpecies;
int64_t bestDups;
int64_t bestLosses;

// Check that the tree rooted above this node doesn't have a lower cost
// than the "best" cost.
static void checkMinimalReconScore(stTree *tree, CuTest *testCase) {
    stTree *newRootedTree = stTree_reRoot(tree, 0.0);
    stPhylogeny_addStIndexedTreeInfo(newRootedTree);
    // This is pretty stupid, but we have to map from the
    // leafToSpecies on the old gene tree to this rerooted
    // one. TODO: Probably the leafToSpecies concept needs to be rethought.
    // Probably should use matrix index -> species node instead.
    stHash *myLeafToSpecies = stHash_construct();
    stHashIterator *hashIt = stHash_getIterator(globalLeafToSpecies);
    stTree *curGene; // Current gene in the old leafToSpecies.
    while((curGene = stHash_getNext(hashIt)) != NULL) {
        stIndexedTreeInfo *index = getIndex(curGene);
        stTree *species = stHash_search(globalLeafToSpecies, curGene);
        stTree *newGene = stPhylogeny_getLeafByIndex(newRootedTree, index->matrixIndex);
        stHash_insert(myLeafToSpecies, newGene, species);
    }

    // Check the cost.
    stPhylogeny_reconcileAtMostBinary(newRootedTree, myLeafToSpecies, false);
    int64_t dups = 0, losses = 0;
    stPhylogeny_reconciliationCostAtMostBinary(newRootedTree, &dups, &losses);
    CuAssertTrue(testCase, dups >= bestDups);
    if (dups == bestDups) {
        CuAssertTrue(testCase, losses >= bestLosses);
    }
    CuAssertTrue(testCase, dups + losses >= 0);
    stHash_destruct(myLeafToSpecies);
    stPhylogenyInfo_destructOnTree(newRootedTree);
    stTree_destruct(newRootedTree);
    stHash_destructIterator(hashIt);
}

static void testStPhylogeny_rootByReconciliationAtMostBinary_simpleTests(CuTest *testCase) {
    // First off -- sanity check that reconciling a gene tree equal to
    // a species tree is a no-op.
    int64_t numLeaves = st_randomInt64(3, 300);
    stMatrix *matrix = getRandomDistanceMatrix(numLeaves);
    stTree *speciesTree = stPhylogeny_neighborJoin(matrix, NULL);
    stTree *geneTree = stTree_clone(speciesTree);
    stHash *leafToSpecies = stHash_construct();
    for(int64_t i = 0; i < numLeaves; i++) {
        stTree *gene = stPhylogeny_getLeafByIndex(geneTree, i);
        stTree *species = stPhylogeny_getLeafByIndex(speciesTree, i);
        stHash_insert(leafToSpecies, gene, species);
    }
    stTree *rooted = stPhylogeny_rootByReconciliationAtMostBinary(geneTree,
                                                                  leafToSpecies);
    CuAssertTrue(testCase, stTree_equals(rooted, geneTree));
    // Check that the cost is 0.
    int64_t dups = 0, losses = 0;
    stPhylogeny_reconcileAtMostBinary(geneTree, leafToSpecies, false);
    stPhylogeny_reconciliationCostAtMostBinary(geneTree, &dups, &losses);
    CuAssertTrue(testCase, dups == 0);
    CuAssertTrue(testCase, losses == 0);

    stPhylogenyInfo_destructOnTree(speciesTree);
    stMatrix_destruct(matrix);
    stTree_destruct(speciesTree);
    stTree_destruct(geneTree);
    stTree_destruct(rooted);
    stHash_destruct(leafToSpecies);

    // Test an annoying example that didn't get rooted right.
    stTree *foo = stTree_parseNewickString("(((0:1,1:1):1,2:1):1,3:1);");
    speciesTree = stTree_parseNewickString("((mouse,rat)mr,human)mrh;");
    leafToSpecies = stHash_construct();
    stHash_insert(leafToSpecies, stTree_findChild(foo, "0"), stTree_findChild(speciesTree, "rat"));
    stHash_insert(leafToSpecies, stTree_findChild(foo, "1"), stTree_findChild(speciesTree, "mouse"));
    stHash_insert(leafToSpecies, stTree_findChild(foo, "2"), stTree_findChild(speciesTree, "human"));
    stHash_insert(leafToSpecies, stTree_findChild(foo, "3"), stTree_findChild(speciesTree, "human"));
    rooted = stPhylogeny_rootByReconciliationAtMostBinary(foo, leafToSpecies);
    char *s = stTree_getNewickTreeString(rooted);
    CuAssertStrEquals(testCase, "((0:1,1:1):0.5,(2:1,3:2):0.5);", s);
    free(s);
    stPhylogenyInfo_destructOnTree(foo);
    stTree_destruct(foo);
    stTree_destruct(speciesTree);
    stTree_destruct(rooted);
    stHash_destruct(leafToSpecies);
    
}

// Make sure that the tree given by rootByReconciliationAtMostBinary
// is a tree with the lowest possible dup cost, with loss cost as
// tiebreaker.
static void testStPhylogeny_rootByReconciliationAtMostBinary_random(CuTest *testCase) {
    for(int64_t testNum = 0; testNum < 1; testNum++) {
        int64_t numSpecies = st_randomInt64(3, 50);
        stMatrix *matrix = getRandomDistanceMatrix(numSpecies);
        globalSpeciesTree = stPhylogeny_neighborJoin(matrix, NULL);
        int64_t numGenes = st_randomInt64(3, 100);
        stMatrix_destruct(matrix);
        matrix = getRandomDistanceMatrix(numGenes);
        stTree *geneTree = stPhylogeny_neighborJoin(matrix, NULL);
        stMatrix_destruct(matrix);

        // Assign genes to random species.
        globalLeafToSpecies = stHash_construct();
        for(int64_t i = 0; i < numGenes; i++) {
            stTree *gene = stPhylogeny_getLeafByIndex(geneTree, i);
            stTree *species = stPhylogeny_getLeafByIndex(globalSpeciesTree, st_randomInt64(0, numSpecies));
            stHash_insert(globalLeafToSpecies, gene, species);
        }

        // Find the best rooting using either the naive or optimized version.
        stTree *rooted;
        if (st_random() > 0.5) {
            rooted = stPhylogeny_rootByReconciliationAtMostBinary(geneTree, globalLeafToSpecies);
        } else {
            rooted = stPhylogeny_rootByReconciliationNaive(geneTree, globalLeafToSpecies);
        }
        stPhylogeny_addStIndexedTreeInfo(rooted);
        // This is pretty stupid, but we have to map from the
        // leafToSpecies on the old gene tree to this rerooted
        // one.
        // Probably should use matrix index -> species node instead.
        stHash *myLeafToSpecies = stHash_construct();
        stHashIterator *hashIt = stHash_getIterator(globalLeafToSpecies);
        stTree *curGene; // Current gene in the old leafToSpecies.
        while((curGene = stHash_getNext(hashIt)) != NULL) {
            stIndexedTreeInfo *index = getIndex(curGene);
            stTree *species = stHash_search(globalLeafToSpecies, curGene);
            stTree *newGene = stPhylogeny_getLeafByIndex(rooted, index->matrixIndex);
            stHash_insert(myLeafToSpecies, newGene, species);
        }
        stHash_destructIterator(hashIt);
        stPhylogeny_reconcileAtMostBinary(rooted, myLeafToSpecies, false);
        int64_t dups = 0, losses = 0;
        stPhylogeny_reconciliationCostAtMostBinary(rooted, &dups, &losses);
        bestDups = dups;
        bestLosses = losses;
        CuAssertTrue(testCase, bestDups >= 0);
        CuAssertTrue(testCase, losses >= 0);

        // Now check all possible roots and confirm that there isn't a
        // better one.
        testOnTree(testCase, geneTree, checkMinimalReconScore);
        stPhylogenyInfo_destructOnTree(rooted);
        stTree_destruct(rooted);
        stPhylogenyInfo_destructOnTree(geneTree);
        stTree_destruct(geneTree);
        stPhylogenyInfo_destructOnTree(globalSpeciesTree);
        stTree_destruct(globalSpeciesTree);
        stHash_destruct(myLeafToSpecies);
    }
}

/*
 * Test that the atMostBinary reconciliation functions work correctly
 * with trees that have nodes that have in-degree 1, out-degree 1.
 */
static void testStPhylogeny_reconcileAtMostBinary_degree2Nodes(CuTest *testCase) {
    stTree *speciesTree = stTree_parseNewickString("(((A,B)F)G,((C,D)H,(E)I)J)K;");
    // First test a gene tree that similar to the species tree, but A
    // and B are children of G now, and E is a child of J. It should
    // have no reconciliation cost and the best root should be the
    // current one.
    stTree *geneTree = stTree_parseNewickString("((A,B)G,((C,D)H,E)J)K;");
    // Create the leaf->species hash (mapping A->A, B->B, etc.)
    stHash *leafToSpecies = stHash_construct();
    for (char label = 'A'; label <= 'E'; label++) {
        char labelStr[2];
        labelStr[0] = label;
        labelStr[1] = '\0';
        stTree *leaf = stTree_findChild(geneTree, labelStr);
        stTree *species = stTree_findChild(speciesTree, labelStr);
        stHash_insert(leafToSpecies, leaf, species);
    }
    stPhylogeny_reconcileAtMostBinary(geneTree, leafToSpecies, true);
    int64_t dups = 0, losses = 0;
    stPhylogeny_reconciliationCostAtMostBinary(geneTree, &dups, &losses);
    CuAssertIntEquals(testCase, 0, dups);
    CuAssertIntEquals(testCase, 0, losses);

    stHash_destruct(leafToSpecies);
    stPhylogenyInfo_destructOnTree(geneTree);
    stTree_destruct(geneTree);

    // Check that if E is gone we still get a (single) loss.
    geneTree = stTree_parseNewickString("((A,B)F,(C,D)H)K;");
    // Create the leaf->species hash (mapping A->A, B->B, etc.)
    leafToSpecies = stHash_construct();
    for (char label = 'A'; label <= 'D'; label++) {
        char labelStr[2];
        labelStr[0] = label;
        labelStr[1] = '\0';
        stTree *leaf = stTree_findChild(geneTree, labelStr);
        stTree *species = stTree_findChild(speciesTree, labelStr);
        stHash_insert(leafToSpecies, leaf, species);
    }
    stPhylogeny_reconcileAtMostBinary(geneTree, leafToSpecies, true);
    dups = 0;
    losses = 0;
    stPhylogeny_reconciliationCostAtMostBinary(geneTree, &dups, &losses);
    CuAssertIntEquals(testCase, 0, dups);
    CuAssertIntEquals(testCase, 1, losses);

    stHash_destruct(leafToSpecies);
    stPhylogenyInfo_destructOnTree(geneTree);
    stTree_destruct(geneTree);
    stTree_destruct(speciesTree);
}

// Build a leaf to species hash using the labels of the leaves,
// assuming that the label format is 'species-gene'.
static stHash *buildLeafToSpeciesUsingDashSeparator(stTree *geneTree, stTree *speciesTree) {
    stHash *leafToSpecies = stHash_construct();
    stList *dfStack = stList_construct();
    stList_append(dfStack, geneTree);
    while (stList_length(dfStack) != 0) {
        stTree *node = stList_pop(dfStack);
        for (int64_t i = 0; i < stTree_getChildNumber(node); i++) {
            stList_append(dfStack, stTree_getChild(node, i));
        }
        if (stTree_getChildNumber(node) == 0) {
            // Assign the leaves to their species by using the first
            // '-'-separated token.
            stList *tokens = stString_splitByString(stTree_getLabel(node), "-");
            char *speciesName = stList_get(tokens, 0);
            stTree *species = stTree_findChild(speciesTree, speciesName);
            assert(species != NULL);
            stHash_insert(leafToSpecies, node, species);
            stList_destruct(tokens);
        }
    }
    stList_destruct(dfStack);
    return leafToSpecies;
}

static void testStPhylogeny_reconcileNonBinary(CuTest *testCase) {
    // TEST 1: "exercise5" tree provided with Notung.
    stTree *speciesTree = stTree_parseNewickString("((tasmanian_devil,opossum,bandicoot,kangaroo)Metatheria,((mouse,human)Euarchontoglires,cow)Boreoeutheria)Theria;");
    stTree *geneTree = stTree_parseNewickString("((((opossum-gene1,tasmanian_devil-gene1),kangaroo-gene1),((mouse-gene,human-gene),cow-gene)),((opossum-gene2,kangaroo-gene2),(bandicoot-gene3,kangaroo-gene3)));");
    // Assign the initial reconciliation
    stHash *leafToSpecies = buildLeafToSpeciesUsingDashSeparator(geneTree, speciesTree);

    // Do the reconciliation and check that the mapping matches.
    stPhylogeny_reconcileNonBinary(geneTree, leafToSpecies, true);
    char *testNewick = stTree_getNewickTreeString(geneTree);
    CuAssertStrEquals(testCase, "((((opossum-gene1,tasmanian_devil-gene1)Metatheria,kangaroo-gene1)Metatheria,((mouse-gene,human-gene)Euarchontoglires,cow-gene)Boreoeutheria)Theria,((opossum-gene2,kangaroo-gene2)Metatheria,(bandicoot-gene3,kangaroo-gene3)Metatheria)Metatheria)Theria;", testNewick);
    free(testNewick);

    // Check that the ILSed node in the tree is not marked as a
    // duplication (since it's only a conditional duplication).
    stTree *node = stTree_getMRCA(stTree_findChild(geneTree, "opossum-gene1"), stTree_findChild(geneTree, "kangaroo-gene1"));
    stPhylogenyInfo *info = stTree_getClientData(node);
    CuAssertTrue(testCase, info->recon->event == SPECIATION);

    // Check that the duplication along lineage "2" is correctly identified.
    node = stTree_getMRCA(stTree_findChild(geneTree, "opossum-gene2"), stTree_findChild(geneTree, "kangaroo-gene3"));
    info = stTree_getClientData(node);
    CuAssertTrue(testCase, info->recon->event == DUPLICATION);

    // Check that the duplication in the root is correctly identified.
    node = geneTree;
    info = stTree_getClientData(node);
    CuAssertTrue(testCase, info->recon->event == DUPLICATION);

    // Clean up.
    stPhylogenyInfo_destructOnTree(geneTree);
    stTree_destruct(geneTree);
    stTree_destruct(speciesTree);
    stHash_destruct(leafToSpecies);

    // TEST 2: Fig 3 in Notung paper.
    speciesTree = stTree_parseNewickString("(A,B,(C,D)Beta)Alpha;");
    geneTree = stTree_parseNewickString("(A-g1,((B-g3,C-g3),D-g2));");

    leafToSpecies = buildLeafToSpeciesUsingDashSeparator(geneTree, speciesTree);

    // Reconcile and check the reconciliation gene->species
    // assignments are correct
    stPhylogeny_reconcileNonBinary(geneTree, leafToSpecies, true);
    testNewick = stTree_getNewickTreeString(geneTree);
    CuAssertStrEquals(testCase, "(A-g1,((B-g3,C-g3)Alpha,D-g2)Alpha)Alpha;", testNewick);
    free(testNewick);

    // Even though there are 3 nodes reconciled to Alpha, there should
    // be only one duplication node, and it should be where the "g3"
    // and "g2" lineages coalesce.

    node = stTree_getMRCA(stTree_findChild(geneTree, "B-g3"), stTree_findChild(geneTree, "C-g3"));
    info = stTree_getClientData(node);
    CuAssertTrue(testCase, info->recon->event == SPECIATION);

    node = stTree_getMRCA(stTree_findChild(geneTree, "B-g3"), stTree_findChild(geneTree, "D-g2"));
    info = stTree_getClientData(node);
    CuAssertTrue(testCase, info->recon->event == DUPLICATION);

    node = stTree_getMRCA(stTree_findChild(geneTree, "B-g3"), stTree_findChild(geneTree, "A-g1"));
    info = stTree_getClientData(node);
    CuAssertTrue(testCase, info->recon->event == SPECIATION);

    // Clean up.
    stPhylogenyInfo_destructOnTree(geneTree);
    stTree_destruct(geneTree);
    stTree_destruct(speciesTree);
    stHash_destruct(leafToSpecies);
}

static void testStPhylogeny_nni(CuTest *testCase) {
    stTree *geneTree = stTree_parseNewickString("(NZOHlLtJ.chr9_87081620_87115187|0|16656-16694.30:0.000120563,((C57B6J.9_87174839_87208406|16656-16694.30:0.0001,C57B6NJ.chr9_90403979_90439438|16656-16694.27:0.0001)27:0.00132877,(C57B6NJ.chr9_90403979_90439438|18548-18586.29:0.00113208,(C57B6NJ.chr1_3987551_4018377|20930-20968.30:0.0001,C57B6NJ.chr1_3987551_4018377|15210-15248.30:0.0001)branch:0.0554717)27:0.00291265)30:0.000120563)30;");
    stTree *tree1, *tree2;
    stPhylogeny_nni(stTree_findChild(geneTree, "branch"), &tree1, &tree2);
    CuAssertTrue(testCase, !stTree_equals(geneTree, tree1));
    CuAssertTrue(testCase, !stTree_equals(geneTree, tree2));
    // Truth set. OK, these are rooted and have branch lengths, so
    // it's concievable this test could fail without the method being
    // wrong. But right now we have no way to test if the unrooted
    // topology is equal.
    stTree *neighbor1 = stTree_parseNewickString("(NZOHlLtJ.chr9_87081620_87115187|0|16656-16694.30:0.000120563,((C57B6NJ.chr9_90403979_90439438|18548-18586.29:0.00113208,(C57B6NJ.chr1_3987551_4018377|20930-20968.30:0.0001,(C57B6J.9_87174839_87208406|16656-16694.30:0.0001,C57B6NJ.chr9_90403979_90439438|16656-16694.27:0.0001)27:0.00132877)branch:0.0554717)27:0.00291265,C57B6NJ.chr1_3987551_4018377|15210-15248.30:0.0001)30:0.000120563)30;");
    stTree *neighbor2 = stTree_parseNewickString("(NZOHlLtJ.chr9_87081620_87115187|0|16656-16694.30:0.000120563,((C57B6J.9_87174839_87208406|16656-16694.30:0.0001,C57B6NJ.chr9_90403979_90439438|16656-16694.27:0.0001)27:0.00132877,((C57B6NJ.chr1_3987551_4018377|20930-20968.30:0.0001,C57B6NJ.chr9_90403979_90439438|18548-18586.29:0.00113208)branch:0.0554717,C57B6NJ.chr1_3987551_4018377|15210-15248.30:0.0001)27:0.00291265)30:0.000120563)30;");
    CuAssertTrue(testCase, stTree_equals(tree1, neighbor1) || stTree_equals(tree1, neighbor2));
    CuAssertTrue(testCase, stTree_equals(tree2, neighbor1) || stTree_equals(tree2, neighbor2));
    CuAssertTrue(testCase, !stTree_equals(tree1, tree2));

    stTree_destruct(geneTree);
    stTree_destruct(tree1);
    stTree_destruct(tree2);
    stTree_destruct(neighbor1);
    stTree_destruct(neighbor2);
}

static void testStPhylogeny_getSplits(CuTest *testCase) {
    // Check against the example given in Bandelt and Dress, 1992.
    double distanceArray[7][7] = {
        { 0,  4,  5,  7, 13,  8,  6},
        { 4,  0,  1,  3,  9, 12, 10},
        { 5,  1,  0,  2,  8, 13, 11},
        { 7,  3,  2,  0,  6, 11, 13},
        {13,  9,  8,  6,  0,  5,  7},
        { 8, 12, 13, 11,  5,  0,  2},
        { 6, 10, 11, 13,  7,  2,  0}
    };
    stMatrix *distanceMatrix = stMatrix_construct(7, 7);
    for (int64_t i = 0; i < 7; i++) {
        for (int64_t j = 0; j < 7; j++) {
            *stMatrix_getCell(distanceMatrix, i, j) = distanceArray[i][j];
        }
    }
    stList *splits = stPhylogeny_getSplits(distanceMatrix, true);
    CuAssertIntEquals(testCase, 4, stList_length(splits));

    // Check the individual splits. We don't check the actual indices
    // within the splits, because that would be very long, and it
    // would be very hard to be wrong about the content of the
    // splits without the isolation index / split sizes being wrong as
    // well.
    stSplit *split1 = stList_get(splits, 0);
    CuAssertDblEquals(testCase, 6.0, split1->isolationIndex, 0.01);
    CuAssertTrue(testCase, stList_length(split1->leftSplit) == 3 || stList_length(split1->rightSplit) == 3);
    CuAssertTrue(testCase, stList_length(split1->leftSplit) == 4 || stList_length(split1->rightSplit) == 4);

    stSplit *split2 = stList_get(splits, 1);
    CuAssertDblEquals(testCase, 4.0, split2->isolationIndex, 0.01);
    CuAssertTrue(testCase, stList_length(split2->leftSplit) == 3 || stList_length(split2->rightSplit) == 3);
    CuAssertTrue(testCase, stList_length(split2->leftSplit) == 4 || stList_length(split2->rightSplit) == 4);

    stSplit *split3 = stList_get(splits, 2);
    CuAssertDblEquals(testCase, 2.0, split3->isolationIndex, 0.01);
    CuAssertTrue(testCase, stList_length(split3->leftSplit) == 3 || stList_length(split3->rightSplit) == 3);
    CuAssertTrue(testCase, stList_length(split3->leftSplit) == 4 || stList_length(split3->rightSplit) == 4);

    stSplit *split4 = stList_get(splits, 3);
    CuAssertDblEquals(testCase, 1.0, split4->isolationIndex, 0.01);
    CuAssertTrue(testCase, stList_length(split4->leftSplit) == 3 || stList_length(split4->rightSplit) == 3);
    CuAssertTrue(testCase, stList_length(split4->leftSplit) == 4 || stList_length(split4->rightSplit) == 4);

    stMatrix_destruct(distanceMatrix);
    stList_destruct(splits);
}

static void testStPhylogeny_greedySplitDecomposition(CuTest *testCase) {
    /*
  _______________ MOUSE.2
 |
 |                ________________ MOUSE.0
 |               |
 |               |                 ________________ RHESUS.0
 |               |________________|
 |_______________|                |                 ________________ CHIMP.0
 |               |                |________________|
 |               |                                 |________________ HUMAN.0
 |               |
_|               |                 ________________ CAT.0
 |               |________________|
 |                                |________________ DOG.0
 |
 |                ________________ CAT.2
 |_______________|
 |               |________________ DOG.2
 |
 |                ________________ CHIMP.2
 |               |
 |_______________|________________ HUMAN.2
                 |
                 |________________ RHESUS.2
    */
    // array ordering: ['CHIMP.0', 'CAT.2', 'CHIMP.2', 'CAT.0', 'HUMAN.2', 'HUMAN.0', 'RHESUS.0', 'RHESUS.2', 'MOUSE.0', 'DOG.0', 'MOUSE.2', 'DOG.2']
    double distanceArray[12][12] = {
        { 0.0,    0.225,  0.22,   0.06,   0.215,  0.0,    0.015,  0.21,   0.075,  0.06,   0.235,  0.255,},
        { 0.225,  0.0,    0.05,   0.22,   0.045,  0.225,  0.23,   0.04,   0.235,  0.225,   0.09,   0.03, },
        { 0.22,   0.05,   0.0,    0.215,  0.015,  0.22,   0.225,  0.01,   0.23,   0.22,   0.08,   0.08, },
        { 0.06,   0.22,   0.215,  0.0,    0.21,   0.06,   0.07,   0.205,  0.07,   0.02,   0.24,   0.25, },
        { 0.215,  0.045,  0.015,  0.21,   0.0,    0.215,  0.22,   0.005,  0.225,  0.215,   0.075,  0.075,},
        { 0.0,    0.225,  0.22,   0.06,   0.215,  0.0,    0.015,  0.21,   0.075,  0.06,   0.235,  0.255,},
        { 0.015,  0.23,   0.225,  0.07,   0.22,   0.015,  0.0,    0.215,  0.085,  0.07,   0.24,   0.26, },
        { 0.21,   0.04,   0.01,   0.205,  0.005,  0.21,   0.215,  0.0,    0.22,   0.21,   0.07,   0.07, },
        { 0.075,  0.235,  0.23,   0.07,   0.225,  0.075,  0.085,  0.22,   0.0,    0.075,   0.255,  0.265,},
        { 0.06,   0.225,  0.22,   0.02,   0.215,  0.06,   0.07,   0.21,   0.075,  0.0,    0.245,   0.255,},
        { 0.235,  0.09,   0.08,   0.24,   0.075,  0.235,  0.24,   0.07,   0.255,  0.245,  0.0,  0.115,},
        { 0.255,  0.03,   0.08,   0.25,   0.075,  0.255,  0.26,   0.07,   0.265,  0.255,   0.115,  0.0  }
    };
    stMatrix *distanceMatrix = stMatrix_construct(12, 12);
    for (int64_t i = 0; i < 12; i++) {
        for (int64_t j = 0; j < 12; j++) {
            *stMatrix_getCell(distanceMatrix, i, j) = distanceArray[i][j];
        }
    }

    stTree *tree = stPhylogeny_greedySplitDecomposition(distanceMatrix, false);
    stTree *mouse2 = stTree_findChild(tree, "10");
    stTree *rooted = stTree_reRoot(mouse2, 0.0);
    char *newick = stTree_getNewickTreeString(rooted);
    CuAssertStrEquals(testCase, "(10:0,((1:1,11:1):1,(2:1,4:1,7:1):1,(8:1,(6:1,(0:1,5:1):1):1,(3:1,9:1):1):1):1);", newick);
    free(newick);
    stMatrix_destruct(distanceMatrix);
    stTree_destruct(rooted);
    stPhylogenyInfo_destructOnTree(tree);
    stTree_destruct(tree);
}

static void testStPhylogeny_getLinkedSpeciesTree(CuTest *testCase) {
    // Hackily including the definition here, so I don't have to
    // create a whole stPhylogeny_private.h and worry about where it
    // gets included. Bad practice.
    stTree *getLinkedSpeciesTree(stTree *speciesTree, stTree *polytomy, stHash **speciesToNumGenes);
    stTree *speciesTree = stTree_parseNewickString("(((a1,(a2.1,a2.2)a2)a,b)e,((c1, c2)c,(d1, d2)d)f)g;");
    stTree *polytomy = stTree_parseNewickString("(a-1,a-2,a-3,a-4,b-1,b-2,c-1,e-1)G;");
    stHash *leafToSpecies = buildLeafToSpeciesUsingDashSeparator(polytomy, speciesTree);
    stPhylogeny_reconcileAtMostBinary(polytomy, leafToSpecies, false);
    stHash *speciesToNumGenes;
    stTree *linkedSpeciesTree = getLinkedSpeciesTree(speciesTree, polytomy, &speciesToNumGenes);
    char *newick = stTree_getNewickTreeString(linkedSpeciesTree);
    CuAssertStrEquals(testCase, "((a,b)e,(c,d)f)g;", newick);
    CuAssertIntEquals(testCase, 4, *((int64_t *) stHash_search(speciesToNumGenes, stTree_findChild(linkedSpeciesTree, "a"))));
    CuAssertIntEquals(testCase, 2, *((int64_t *) stHash_search(speciesToNumGenes, stTree_findChild(linkedSpeciesTree, "b"))));
    CuAssertIntEquals(testCase, 1, *((int64_t *) stHash_search(speciesToNumGenes, stTree_findChild(linkedSpeciesTree, "c"))));
    CuAssertIntEquals(testCase, 1, *((int64_t *) stHash_search(speciesToNumGenes, stTree_findChild(linkedSpeciesTree, "e"))));
    free(newick);
    stTree_destruct(linkedSpeciesTree);
    stHash_destruct(speciesToNumGenes);
    stHash_destruct(leafToSpecies);
    stPhylogenyInfo_destructOnTree(polytomy);
    stTree_destruct(polytomy);
    stTree_destruct(speciesTree);
}

static void testStPhylogeny_reconciliationCostAtMostBinary_polytomies(CuTest *testCase) {
    stTree *speciesTree = stTree_parseNewickString("(((a1,(a2.1,a2.2)a2)a,b)e,((c1, c2)c,(d1, d2)d)f)g;");
    stTree *polytomy = stTree_parseNewickString("(a-1,a-2,a-3,a-4,b-1,b-2,c-1,e-1)G;");
    stHash *leafToSpecies = buildLeafToSpeciesUsingDashSeparator(polytomy, speciesTree);
    stPhylogeny_reconcileAtMostBinary(polytomy, leafToSpecies, false);
    int64_t dups = 0, losses = 0;
    stPhylogeny_reconciliationCostAtMostBinary(polytomy, &dups, &losses);
    CuAssertIntEquals(testCase, 1, losses);
    CuAssertIntEquals(testCase, 4, dups);

    stHash_destruct(leafToSpecies);
    stPhylogenyInfo_destructOnTree(polytomy);
    stTree_destruct(polytomy);
    stTree_destruct(speciesTree);
}

static void testStPhylogeny_applyJukesCantorCorrection(CuTest *testCase) {
    double distanceArray[6][6] = {
        { 0.0,    0.225,  0.22,   0.745,   0.76,  0.0 },
        { 0.225,  0.0,    0.05,   0.22,   0.045,  0.225 },
        { 0.22,   0.05,   0.0,    0.215,  0.015,  0.22 },
        { 0.745,   0.22,   0.215,  0.0,    0.21,   0.06 },
        { 0.76,  0.045,  0.015,  0.21,   0.0,    0.215 },
        { 0.0,    0.225,  0.22,   0.06,   0.215,  0.0 },
    };
    stMatrix *distanceMatrix = stMatrix_construct(6, 6);
    for (int64_t i = 0; i < 6; i++) {
        for (int64_t j = 0; j < 6; j++) {
            *stMatrix_getCell(distanceMatrix, i, j) = distanceArray[i][j];
        }
    }

    stPhylogeny_applyJukesCantorCorrection(distanceMatrix);
    CuAssertDblEquals(testCase, 3.7579, *stMatrix_getCell(distanceMatrix, 3, 0), 0.001);
    CuAssertDblEquals(testCase, 0.0625, *stMatrix_getCell(distanceMatrix, 3, 5), 0.001);
    // Check that values over 0.75 don't return NaN or infinity.
    CuAssertTrue(testCase, isnormal(*stMatrix_getCell(distanceMatrix, 4, 0)));

    stMatrix_destruct(distanceMatrix);
}

// Large tree that takes forever if the reconciliation algorithm isn't
// linear time. This is sort of a hacky way to test that the
// reconciliation algorithm stays linear. The test won't fail if it is
// O(n^2) or exponential, but it will take 30+ mins.
static void testStPhylogeny_reconcileLargeTree_shouldBeFast(CuTest *testCase) {
    stTree *geneTree = stTree_parseNewickString("(((0:0.0338606,(24:0,(((((25:0,28:0.00688078):0,26:0):0,27:0):0,29:0):0.0103492,48:0.0131733):0):0.00135617):0,(49:0.00344003,51:0.0607045):0):0.0036751,((((((((((((((((((((((1:0,((15:0,(((((((((((((((((30:0,31:0):0,32:0):0,33:0):0,34:0):0,35:0):0,36:0):0,37:0):0,38:0):0,39:0):0,40:0):0,41:0):0,42:0):0,43:0):0,44:0):0,46:0):0,47:0):0.00485278,45:0.00495128):0.0185071):0.00778827,16:0):0):0,2:0):0,3:0):0,4:0):0,5:0):0,6:0):0,7:0):0,8:0):0,9:0):0,10:0):0,11:0):0,12:0):0,13:0):0,14:0):0,17:0):0,18:0):0,19:0):0,20:0):0,21:0):0,22:0):0,23:0):0.00330033,50:0):0.0036751);");
    stPhylogeny_addStIndexedTreeInfo(geneTree);
    int64_t numSpecies = 0;
    stTree *speciesTree = getRandomBinaryTree(st_randomInt64(2, 7), &numSpecies);
    stPhylogeny_addStIndexedTreeInfo(speciesTree);
    // Assign genes to random species.
    stHash *leafToSpecies = stHash_construct();
    for(int64_t i = 0; i < 52; i++) {
        stTree *gene = stPhylogeny_getLeafByIndex(geneTree, i);
        stTree *species = stPhylogeny_getLeafByIndex(speciesTree, st_randomInt64(0, numSpecies));
        stHash_insert(leafToSpecies, gene, species);
    }

    stPhylogeny_reconcileAtMostBinary(geneTree, leafToSpecies, false);
    stPhylogenyInfo_destructOnTree(geneTree);
    stTree_destruct(geneTree);
    stPhylogenyInfo_destructOnTree(speciesTree);
    stTree_destruct(speciesTree);
    stHash_destruct(leafToSpecies);
}

CuSuite* sonLib_stPhylogenyTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testStPhylogeny_reconciliationCostAtMostBinary_polytomies);
    SUITE_ADD_TEST(suite, testStPhylogeny_getLinkedSpeciesTree);
    SUITE_ADD_TEST(suite, testStPhylogeny_greedySplitDecomposition);
    SUITE_ADD_TEST(suite, testStPhylogeny_getSplits);
    SUITE_ADD_TEST(suite, testStPhylogeny_nni);
    SUITE_ADD_TEST(suite, testJoinCosts_random);
    SUITE_ADD_TEST(suite, testStPhylogeny_reconcileAtMostBinary_degree2Nodes);
    SUITE_ADD_TEST(suite, testSimpleNeighborJoin);
    SUITE_ADD_TEST(suite, testSimpleBootstrapPartitionScoring);
    SUITE_ADD_TEST(suite, testSimpleBootstrapReconciliationScoring);
    SUITE_ADD_TEST(suite, testRandomNeighborJoin);
    SUITE_ADD_TEST(suite, testRandomBootstraps);
    SUITE_ADD_TEST(suite, testSimpleJoinCosts);
    SUITE_ADD_TEST(suite, testGuidedNeighborJoiningReducesToNeighborJoining);
    SUITE_ADD_TEST(suite, testGuidedNeighborJoiningLowersReconCost);
    SUITE_ADD_TEST(suite, testStPhylogeny_rootByReconciliationAtMostBinary_simpleTests);
    SUITE_ADD_TEST(suite, testStPhylogeny_rootByReconciliationAtMostBinary_random);
    SUITE_ADD_TEST(suite, testStPhylogeny_reconcileNonBinary);
    SUITE_ADD_TEST(suite, testStPhylogeny_applyJukesCantorCorrection);
    SUITE_ADD_TEST(suite, testStPhylogeny_reconcileLargeTree_shouldBeFast);

    (void) testStPhylogeny_reconciliationCostAtMostBinary_polytomies;
    (void) testStPhylogeny_getLinkedSpeciesTree;
    (void) testStPhylogeny_greedySplitDecomposition;
    (void) testStPhylogeny_getSplits;
    (void) testStPhylogeny_nni;
    (void) testJoinCosts_random;
    (void) testStPhylogeny_reconcileAtMostBinary_degree2Nodes;
    (void) testSimpleNeighborJoin;
    (void) testSimpleBootstrapPartitionScoring;
    (void) testSimpleBootstrapReconciliationScoring;
    (void) testRandomNeighborJoin;
    (void) testRandomBootstraps;
    (void) testSimpleJoinCosts;
    (void) testGuidedNeighborJoiningReducesToNeighborJoining;
    (void) testGuidedNeighborJoiningLowersReconCost;
    (void) testStPhylogeny_rootByReconciliationAtMostBinary_simpleTests;
    (void) testStPhylogeny_rootByReconciliationAtMostBinary_random;
    (void) testStPhylogeny_reconcileNonBinary;
    return suite;
}

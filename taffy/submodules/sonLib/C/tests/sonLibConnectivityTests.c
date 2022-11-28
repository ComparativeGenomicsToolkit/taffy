#include "sonLibGlobalsTest.h"
#include <inttypes.h>
#include <time.h>

static stConnectivity *connectivity;

static void teardown(void);

static void setup(void) {
	teardown();
	connectivity = stConnectivity_construct();
	/*
	 *  1--2
	 *  |  |
	 *  3--4
	 *
	 *  5--6--7
	 */
	stConnectivity_addNode(connectivity, (void *) 1);
	stConnectivity_addNode(connectivity, (void *) 2);
	stConnectivity_addNode(connectivity, (void *) 3);
	stConnectivity_addNode(connectivity, (void *) 4);
	stConnectivity_addNode(connectivity, (void *) 5);
	stConnectivity_addNode(connectivity, (void *) 6);
	stConnectivity_addNode(connectivity, (void *) 7);

	stConnectivity_addEdge(connectivity, (void *) 1, (void *) 2);
	stConnectivity_addEdge(connectivity, (void *) 1, (void *) 3);
	stConnectivity_addEdge(connectivity, (void *) 2, (void *) 4);
	stConnectivity_addEdge(connectivity, (void *) 3, (void *) 4);

	stConnectivity_addEdge(connectivity, (void *) 5, (void *) 6);
	stConnectivity_addEdge(connectivity, (void *) 6, (void *) 7);
}

static void teardown(void) {
	if (connectivity != NULL) {
		stConnectivity_destruct(connectivity);
		connectivity = NULL;
	}
}

static bool setsEqual(stSet *set1, stSet *set2) {
	stSet *diff = stSet_getDifference(set1, set2);
	size_t size = stSet_size(diff);
	stSet_destruct(diff);
	return size == 0;
}

static void test_stConnectivity_newNodeShouldGoInANewComponent(CuTest *testCase) {
	setup();
	stConnectivity_addNode(connectivity, (void *) 0xDEADBEEF);
	stConnectedComponent *component = stConnectivity_getConnectedComponent(connectivity, (void *) 0xDEADBEEF);
	stConnectedComponentNodeIterator *nodeIt = stConnectedComponent_getNodeIterator(component);
	int64_t i = 0;
	void *node;
	while ((node = stConnectedComponentNodeIterator_getNext(nodeIt)) != NULL) {
		i++;
		CuAssertTrue(testCase, node == (void *) 0xDEADBEEF);
	}
	CuAssertIntEquals(testCase, 1, i);

	stConnectedComponentNodeIterator_destruct(nodeIt);
	teardown();
}

static void test_stConnectivity_connectedComponents(CuTest *testCase) {
	setup();

	stConnectedComponentIterator *componentIt = stConnectivity_getConnectedComponentIterator(connectivity);
	stList *components = stList_construct();
	stConnectedComponent *component;
	while ((component = stConnectedComponentIterator_getNext(componentIt)) != NULL) {
		stList_append(components, component);
	}
	stConnectedComponentIterator_destruct(componentIt);

	CuAssertIntEquals(testCase, 2, stList_length(components));

	stSet *truth1 = stSet_construct();
	stSet_insert(truth1, (void *) 1);
	stSet_insert(truth1, (void *) 2);
	stSet_insert(truth1, (void *) 3);
	stSet_insert(truth1, (void *) 4);
	stSet *test1 = stSet_construct();
	stConnectedComponent *component1 = stConnectivity_getConnectedComponent(connectivity, (void *) 4);
	stConnectedComponentNodeIterator *nodeIt = stConnectedComponent_getNodeIterator(component1);
	void *node;
	while ((node = stConnectedComponentNodeIterator_getNext(nodeIt)) != NULL) {
		CuAssertTrue(testCase, stSet_search(test1, node) == NULL);
		stSet_insert(test1, node);
	}
	stConnectedComponentNodeIterator_destruct(nodeIt);
	CuAssertTrue(testCase, setsEqual(test1, truth1));
	stSet_destruct(truth1);
	stSet_destruct(test1);

	stSet *truth2 = stSet_construct();
	stSet_insert(truth2, (void *) 5);
	stSet_insert(truth2, (void *) 6);
	stSet_insert(truth2, (void *) 7);
	stSet *test2 = stSet_construct();
	stConnectedComponent *component2 = stConnectivity_getConnectedComponent(connectivity, (void *) 5);
	nodeIt = stConnectedComponent_getNodeIterator(component2);
	while ((node = stConnectedComponentNodeIterator_getNext(nodeIt)) != NULL) {
		CuAssertTrue(testCase, stSet_search(test2, node) == NULL);
		stSet_insert(test2, node);
	}
	stConnectedComponentNodeIterator_destruct(nodeIt);
	CuAssertTrue(testCase, setsEqual(test2, truth2));
	stSet_destruct(truth2);
	stSet_destruct(test2);

	stList_destruct(components);
	teardown();
}

static void test_stConnectivity_removeNodesAndEdges(CuTest *testCase) {
	setup();

	stConnectivity_removeNode(connectivity, (void *) 6);
	stConnectivity_removeEdge(connectivity, (void *) 1, (void *) 2);
	stConnectivity_removeEdge(connectivity, (void *) 3, (void *) 1);

	stConnectedComponentIterator *componentIt = stConnectivity_getConnectedComponentIterator(connectivity);
	stList *components = stList_construct();
	stConnectedComponent *component;
	while ((component = stConnectedComponentIterator_getNext(componentIt)) != NULL) {
		stList_append(components, component);
	}
	stConnectedComponentIterator_destruct(componentIt);

	CuAssertIntEquals(testCase, 4, stList_length(components));

	stList_destruct(components);

	teardown();
}
static void test_stConnectivity_connected(CuTest *testCase) {
	setup();
	/*
	 *  1--2
	 *  |  |
	 *  3--4
	 *
	 *  5--6--7
	 */

	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 1, (void*) 3));
	CuAssertTrue(testCase, !stConnectivity_connected(connectivity, (void*) 3, (void*) 5));
	CuAssertTrue(testCase, !stConnectivity_connected(connectivity, (void*) 4, (void*) 7));
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 2, (void*) 4));
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 2, (void*) 3));
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 1, (void*) 2));
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 3, (void*) 4));

	stConnectivity_removeEdge(connectivity, (void*) 5, (void*) 6);
	CuAssertTrue(testCase, !stConnectivity_connected(connectivity, (void*) 5, (void*) 6));
	stConnectivity_removeEdge(connectivity, (void*) 1, (void*) 2);
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 1, (void*) 3));
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 3, (void*) 4));
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 4, (void*) 2));
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 1, (void*) 2));

	stConnectivity_removeEdge(connectivity, (void*) 2, (void*) 4);
	CuAssertTrue(testCase, !stConnectivity_connected(connectivity, (void*) 1, (void*) 2));
	stConnectivity_addEdge(connectivity, (void*) 2, (void *) 4);

	stConnectivity_removeEdge(connectivity, (void*) 3, (void*)4);
	CuAssertTrue(testCase, !stConnectivity_connected(connectivity, (void*) 1, (void*) 2));
	stConnectivity_addEdge(connectivity, (void*) 3, (void*) 4);
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 3, (void*) 4));
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 1, (void*) 2));


	teardown();
	connectivity = stConnectivity_construct();
	stConnectivity_addNode(connectivity, (void*) 1);
	stConnectivity_addNode(connectivity, (void*) 2);
	stConnectivity_addNode(connectivity, (void*) 3);
	stConnectivity_addNode(connectivity, (void*) 4);
	stConnectivity_addEdge(connectivity, (void*) 1, (void*) 2);
	stConnectivity_addEdge(connectivity, (void*) 1, (void*) 3);
	stConnectivity_addEdge(connectivity, (void*) 1, (void*) 4);

	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*)1, (void*)3));
	stConnectivity_removeEdge(connectivity, (void*) 1, (void*) 4);
	stConnectivity_removeEdge(connectivity, (void*) 1, (void*) 2);
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*)1, (void*) 3));
	teardown();
}
static void test_stConnectivity_nodeIterator(CuTest *testCase) {
	setup();
	stConnectedComponent *comp = stConnectivity_getConnectedComponent(connectivity, (void*)4);
	stConnectedComponentNodeIterator *it = stConnectedComponent_getNodeIterator(comp);
	void *node;
	stSet *nodesInComponent = stSet_construct();
	stSet *truth = stSet_construct();
	stSet_insert(truth, (void*)1);
	stSet_insert(truth, (void*)2);
	stSet_insert(truth, (void*)3);
	stSet_insert(truth, (void*)4);
	while((node = stConnectedComponentNodeIterator_getNext(it))) {
		stSet_insert(nodesInComponent, node);
	}
	CuAssertTrue(testCase, setsEqual(nodesInComponent, truth));
	stSet_destruct(nodesInComponent);
	stSet_destruct(truth);
	stConnectedComponentNodeIterator_destruct(it);

	teardown();
}
static void test_stConnectivity_compareWithNaive(CuTest *testCase) {
	srand(time(NULL));
	int nNodes = 500;
	int nEdgesToAdd = 800;
	int nEdgesToRemove = 500;
	int nQueries = 10000;
	stList *nodes = stList_construct3(0, free);
	stNaiveConnectivity *naive = stNaiveConnectivity_construct();
	connectivity = stConnectivity_construct();
	//add nodes
	for (int i = 1; i <= nNodes; i++) {
		int *node = malloc(sizeof(int));
		*node = i;
		stNaiveConnectivity_addNode(naive, node);
		stConnectivity_addNode(connectivity, node);
		stList_append(nodes, node);
	}
	//add edges
	while(nEdgesToAdd > 0) {
		int *node1 = stList_get(nodes, rand() % nNodes);
		int *node2 = stList_get(nodes, rand() % nNodes);
		if(node1 >= node2) continue;
		CuAssertTrue(testCase, stConnectivity_hasEdge(connectivity, node1, node2) == stNaiveConnectivity_hasEdge(naive, node1, node2));
		if(stConnectivity_hasEdge(connectivity, node1, node2)) continue;
		stConnectivity_addEdge(connectivity, node1, node2);
		stNaiveConnectivity_addEdge(naive, node1, node2);
		nEdgesToAdd--;
	}
	//remove edges
	while(nEdgesToRemove > 0) {
		int *node1 = stList_get(nodes, rand() % nNodes);
		int *node2 = stList_get(nodes, rand() % nNodes);
		if(node1 >= node2) continue;
		CuAssertTrue(testCase, stConnectivity_hasEdge(connectivity, node1, node2) == stNaiveConnectivity_hasEdge(naive, node1, node2));
		if(!stConnectivity_hasEdge(connectivity, node1, node2)) continue;
		stConnectivity_removeEdge(connectivity, node1, node2);
		stNaiveConnectivity_removeEdge(naive, node1, node2);
		nEdgesToRemove--;
	}
	//check connectivity queries
	while(nQueries > 0) {
		int *node1 = stList_get(nodes, rand() % nNodes);
		int *node2 = stList_get(nodes, rand() % nNodes);
		if(node1 >= node2) continue;
		stNaiveConnectedComponent *comp1 = stNaiveConnectivity_getConnectedComponent(naive, node1);
		stNaiveConnectedComponent *comp2 = stNaiveConnectivity_getConnectedComponent(naive, node2);
		bool naiveConnected = comp1 == comp2;
		if(!naiveConnected == stConnectivity_connected(connectivity, node1, node2)) {
			printf("%d and %d should be connected: %d\n", *node1, *node2, naiveConnected);
		}
		CuAssertTrue(testCase, naiveConnected == stConnectivity_connected(connectivity, node1, node2));
		nQueries--;
	}


	//check number of connected components
	stList *components = stList_construct();
	stConnectedComponentIterator *it = stConnectivity_getConnectedComponentIterator(connectivity);
	stConnectedComponent *comp;
	while((comp = stConnectedComponentIterator_getNext(it))) {
		stList_append(components, comp);
	}
	stConnectedComponentIterator_destruct(it);

	stList *trueComponents = stList_construct();
	stNaiveConnectedComponentIterator *itNaive = stNaiveConnectivity_getConnectedComponentIterator(naive);
	stNaiveConnectedComponent *naiveComp;
	while((naiveComp = stNaiveConnectedComponentIterator_getNext(itNaive))) {
		stList_append(trueComponents, naiveComp);
	}
	stNaiveConnectedComponentIterator_destruct(itNaive);
	CuAssertTrue(testCase, stList_length(components) == stList_length(trueComponents));
	//check the nodes in each component
	for (int i = 0; i < stList_length(components); i++) {
		stNaiveConnectedComponent *truecomponent_i = stList_get(trueComponents, i);
		stSet *trueNodesInComponent = stNaiveConnectedComponent_getNodes(truecomponent_i);
		stList *trueNodesInComponentList = stSet_getList(trueNodesInComponent);
		void *nodeInComponent = stList_get(trueNodesInComponentList, 0);
		stList_destruct(trueNodesInComponentList);

		stSet *nodesInComponent = stSet_construct();
		stConnectedComponent *comp_i = stConnectivity_getConnectedComponent(connectivity, nodeInComponent);
		stConnectedComponentNodeIterator *nodeIt = 
			stConnectedComponent_getNodeIterator(comp_i);
		void *node;
		while((node = stConnectedComponentNodeIterator_getNext(nodeIt))) {
			stSet_insert(nodesInComponent, node);
		}
		stConnectedComponentNodeIterator_destruct(nodeIt);

		CuAssertTrue(testCase, setsEqual(nodesInComponent, trueNodesInComponent));
		stSet_destruct(nodesInComponent);
	}
	stList_destruct(nodes);
	stList_destruct(trueComponents);
	stList_destruct(components);
	stNaiveConnectivity_destruct(naive);

	teardown();
}

// Very simple test that multigraphs work properly.
static void test_stConnectivity_multigraphs(CuTest *testCase) {
	setup();
	stConnectivity_addEdge(connectivity, (void*) 5, (void*) 6);

	stConnectivity_removeEdge(connectivity, (void*) 5, (void*) 6);
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 5, (void*) 6));
	stConnectivity_removeEdge(connectivity, (void*) 5, (void*) 6);
	CuAssertTrue(testCase, !stConnectivity_connected(connectivity, (void*) 5, (void*) 6));

	stConnectivity_addEdge(connectivity, (void*) 5, (void*) 6);
	stConnectivity_addEdge(connectivity, (void*) 5, (void*) 6);
	CuAssertTrue(testCase, stConnectivity_connected(connectivity, (void*) 5, (void*) 7));
	stConnectivity_removeNode(connectivity, (void*) 6);
	CuAssertTrue(testCase, !stConnectivity_connected(connectivity, (void*) 5, (void*) 7));

	teardown();
}

// Test that components always get the same pointer unless those specific
// components have been modified. (Previously the API spec dictated that
// the components must be invalidated after any edge insertion or
// deletion).
static void test_stConnectivity_constantComponentPointers(CuTest *testCase) {
	setup();

	stConnectedComponent *component1 = stConnectivity_getConnectedComponent(connectivity, (void *) 2);
	stConnectedComponent *component2 = stConnectivity_getConnectedComponent(connectivity, (void *) 5);
	stConnectivity_addEdge(connectivity, (void *) 2, (void *) 3); // Doesn't affect components

	CuAssertTrue(testCase, component1 == stConnectivity_getConnectedComponent(connectivity, (void *) 2));
	CuAssertTrue(testCase, component2 == stConnectivity_getConnectedComponent(connectivity, (void *) 5));

	stConnectivity_removeEdge(connectivity, (void *) 2, (void *) 3); // Doesn't affect components

	CuAssertTrue(testCase, component1 == stConnectivity_getConnectedComponent(connectivity, (void *) 2));
	CuAssertTrue(testCase, component2 == stConnectivity_getConnectedComponent(connectivity, (void *) 5));

	teardown();
}

typedef enum {
	CREATION,
	DELETION,
	MERGE,
	CLEAVE
} componentCallback;

componentCallback callbackResults[256];

int64_t curCallbackResultsIdx = 0;

static void createCallback(void *ignored, stConnectedComponent *adjComponent) {
	callbackResults[curCallbackResultsIdx++] = CREATION;
}

static void deleteCallback(void *ignored, stConnectedComponent *adjComponent) {
	callbackResults[curCallbackResultsIdx++] = DELETION;
}

static void mergeCallback(void *ignored, stConnectedComponent *adjComponent1, stConnectedComponent *adjComponent2) {
	callbackResults[curCallbackResultsIdx++] = MERGE;
}

static void cleaveCallback(void *ignored, stConnectedComponent *adjComponent1, stConnectedComponent *adjComponent2, stSet *vertices) {
	callbackResults[curCallbackResultsIdx++] = CLEAVE;
}

static void test_stConnectivity_callbacks(CuTest *testCase) {
	connectivity = stConnectivity_construct();

	stConnectivity_setCreationCallback(connectivity, createCallback, NULL);
	stConnectivity_setDeletionCallback(connectivity, deleteCallback, NULL);
	stConnectivity_setMergeCallback(connectivity, mergeCallback, NULL);
	stConnectivity_setCleaveCallback(connectivity, cleaveCallback, NULL);

	// Create 2 nodes
	stConnectivity_addNode(connectivity, (void *) 1);
	stConnectivity_addNode(connectivity, (void *) 2);

	// Merge
	stConnectivity_addEdge(connectivity, (void *) 1, (void *) 2);

	// Cleave
	stConnectivity_removeEdge(connectivity, (void *) 1, (void *) 2);

	// Delete
	stConnectivity_removeNode(connectivity, (void *) 1);

	CuAssertIntEquals(testCase, 5, curCallbackResultsIdx);

	CuAssertTrue(testCase, callbackResults[0] == CREATION);
	CuAssertTrue(testCase, callbackResults[1] == CREATION);
	CuAssertTrue(testCase, callbackResults[2] == MERGE);
	CuAssertTrue(testCase, callbackResults[3] == CLEAVE);
	CuAssertTrue(testCase, callbackResults[4] == DELETION);

	teardown();
}

CuSuite *sonLib_stConnectivityTestSuite(void) {
	CuSuite *suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, test_stConnectivity_newNodeShouldGoInANewComponent);
	SUITE_ADD_TEST(suite, test_stConnectivity_connectedComponents);
	SUITE_ADD_TEST(suite, test_stConnectivity_removeNodesAndEdges);
	SUITE_ADD_TEST(suite, test_stConnectivity_connected);
	SUITE_ADD_TEST(suite, test_stConnectivity_nodeIterator);
	SUITE_ADD_TEST(suite, test_stConnectivity_compareWithNaive);
	SUITE_ADD_TEST(suite, test_stConnectivity_multigraphs);
	SUITE_ADD_TEST(suite, test_stConnectivity_constantComponentPointers);
	SUITE_ADD_TEST(suite, test_stConnectivity_callbacks);
	return suite;
}


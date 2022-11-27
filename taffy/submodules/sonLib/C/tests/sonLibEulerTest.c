#include "sonLibGlobalsTest.h"

stEulerTour *et;

static void setup(void) {
	et = stEulerTour_construct();

	stEulerTour_createVertex(et, (void*)1);
	stEulerTour_createVertex(et, (void*)2);
	stEulerTour_createVertex(et, (void*)3);
	stEulerTour_createVertex(et, (void*)4);
	stEulerTour_createVertex(et, (void*)5);

}
static void teardown(void) {
	stEulerTour_destruct(et);
}
static void test_stEulerTour_makeRoot(CuTest *testCase) {
	setup();
	stEulerVertex *a = stEulerTour_getVertex(et, (void*)1);
	stEulerVertex *b = stEulerTour_getVertex(et, (void*)2);


	stEulerTour_link(et, (void*)1, (void*)2);
	stEulerTour_makeRoot(et, a);
	CuAssertTrue(testCase, stEulerTour_findRootNode(et, (void*)1) == (void*)1);
	CuAssertTrue(testCase, stEulerTour_findRootNode(et, (void*)2) == (void*)1);
	stEulerTour_makeRoot(et, b);
	CuAssertTrue(testCase, stEulerTour_findRootNode(et, (void*)1) == (void*)2);
	teardown();
}

static void test_stEulerTour_circle(CuTest *testCase) {
	setup();
	stEulerVertex *c = stEulerTour_getVertex(et, (void*)3);
	stEulerVertex *d = stEulerTour_getVertex(et, (void*)4);
	stEulerTour_link(et, (void*)1, (void*)2);
	stEulerTour_link(et, (void*)1, (void*)3);
	stEulerTour_link(et, (void*)2, (void*)4);
	stEulerTour_cut(et, (void*)1, (void*)2);
	stEulerTour_makeRoot(et, c);
	stEulerTour_makeRoot(et, d);
	stEulerTour_link(et, (void*)3, (void*)4);
	teardown();
}

static void test_stEulerTour_link(CuTest *testCase) {
	setup();
	//stEulerVertex *a = stEulerTour_getVertex(et, (void*)1);
	//stEulerVertex *b = stEulerTour_getVertex(et, (void*)2);


	stEulerTour_link(et, (void*)1, (void*)2);
	CuAssertTrue(testCase, stEulerTour_connected(et, (void*)1, (void*)2));

	CuAssertTrue(testCase, !stEulerTour_connected(et, (void*)2, (void*)3));
	CuAssertTrue(testCase, !stEulerTour_connected(et, (void*)1, (void*)3));

	/*
	   CuAssertTrue(testCase, a->leftOut->to == b);
	   CuAssertTrue(testCase, a->leftOut->from == a);
	   CuAssertTrue(testCase, b->rightIn->from == b);
	   CuAssertTrue(testCase, a->rightIn->to == a);
	   CuAssertTrue(testCase, stTreap_compare(a->leftOut->node, a->rightIn->node) < 0);
	   CuAssertTrue(testCase, stTreap_compare(a->leftOut->node, b->leftOut->node) == 0);

	   CuAssertTrue(testCase, stTreap_next(a->leftOut->node) == b->rightIn->node);
	   */

	CuAssertTrue(testCase, stEulerTour_findRootNode(et, (void*)1) == (void*)1);
	CuAssertTrue(testCase, (char*)stEulerTour_findRootNode(et, (void*)2) == (void*)1);


	/*
	 *  d--a--b
	 *      \
	 *       c
	 */

	stEulerTour_link(et, (void*)1, (void*)3);
	stEulerTour_link(et, (void*)2, (void*)4);
	//CuAssertStrEquals(testCase, stEulerVertex_getTour(a), "dacaba");
	CuAssertTrue(testCase, stEulerTour_connected(et, (void*)4, (void*)2));

	teardown();
}
static bool setsEqual(stSet *set1, stSet *set2) {
	stSet *diff = stSet_getDifference(set1, set2);
	size_t size = stSet_size(diff);
	stSet_destruct(diff);
	return size == 0;
}

static void test_stEulerTour_multipleIncidentEdges(CuTest *testCase) {
	setup();
	stSet *trueNodes = stSet_construct();
	stEulerTour_link(et, (void*)1, (void*)2);
	stEulerTour_link(et, (void*)1, (void*)3);
	stEulerTour_link(et, (void*)1, (void*)4);
	stSet_insert(trueNodes, (void*)1);
	stSet_insert(trueNodes, (void*)2);
	stSet_insert(trueNodes, (void*)3);
	stSet_insert(trueNodes, (void*)4);

	stEulerTour_cut(et, (void*)1, (void*)2);
	stSet_remove(trueNodes, (void*)2);

	CuAssertTrue(testCase, stEulerTour_connected(et, (void*)1, (void*)3));
	stSet *nodes = stEulerTour_getNodesInComponent(et, (void*)1);
	CuAssertTrue(testCase, setsEqual(nodes, trueNodes));
	stEulerTour_cut(et, (void*)1, (void*)4);

	CuAssertTrue(testCase, stEulerTour_connected(et, (void*)1, (void*)3));
	stSet_destruct(nodes);
	stSet_destruct(trueNodes);
	teardown();
}
static void test_stEulerTour_cut(CuTest *testCase) {
	setup();
	stEulerTour_link(et, (void*)1, (void*)2);
	stEulerTour_link(et, (void*)1, (void*)3);
	stEulerTour_link(et, (void*)1, (void*)4);
	stEulerTour_link(et, (void*)2, (void*)5);
	stEulerTour_cut(et, (void*)1, (void*)2);
	//CuAssertTrue(testCase, et->nComponents == 2);
	CuAssertTrue(testCase, !stEulerTour_connected(et, (void*)1, (void*)2));
	stEulerTour_link(et, (void*)1, (void*)2);
	CuAssertTrue(testCase, stEulerTour_connected(et, (void*)1, (void*)2));
	teardown();
}
static void test_stEulerTour_next(CuTest *testCase) {
	setup();
	/*
	   stEulerTour_link(et, (void*)1, (void*)2);
	   stEulerTour_link(et, (void*)1, (void*)3);
	   stEulerTour_link(et, (void*)1, (void*)4);
	//stEulerTour_link(et, (void*)2, (void*)5);
	stEulerTourIterator *it = stEulerTour_getIterator(et, (void*) 3);
	void *node;
	while((node = stEulerTourIterator_getNext(it))) {
	printf("%p \n", node);
	}
	stEulerTourIterator_destruct(it);
	stEulerTourEdgeIterator *edgeIt = stEulerTour_getEdgeIterator(et, (void*) 3);
	stEulerHalfEdge *edge;
	while((edge = stEulerTourEdgeIterator_getNext(edgeIt))) {
	printf("%p->%p ", stEulerHalfEdge_getFrom(edge), stEulerHalfEdge_getTo(edge));
	}
	stEulerTourEdgeIterator_destruct(edgeIt);
	*/
	teardown();
}


CuSuite *sonLib_stEulerTestSuite(void) {
	CuSuite *suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, test_stEulerTour_link);
	SUITE_ADD_TEST(suite, test_stEulerTour_cut);
	SUITE_ADD_TEST(suite, test_stEulerTour_next);
	SUITE_ADD_TEST(suite, test_stEulerTour_multipleIncidentEdges);
	SUITE_ADD_TEST(suite, test_stEulerTour_circle);
	SUITE_ADD_TEST(suite, test_stEulerTour_makeRoot);
	return suite;
}

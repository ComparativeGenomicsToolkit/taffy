#include "sonLibGlobalsTest.h"
stEdgeContainer *container;
void destruct(char *v) {
	free(v);
}
static void setup(void) {
	container = stEdgeContainer_construct((void(*)(void*))destruct);
	char *a = st_malloc(2);
	strcpy(a, "a");
	char *b = st_malloc(2);
	strcpy(b, "b");
	void *c = st_malloc(2);
	strcpy(c, "c");
	
	stEdgeContainer_addEdge(container, (void*)1, (void*)2, a);
	stEdgeContainer_addEdge(container, (void*)1, (void*)3, b);
	stEdgeContainer_addEdge(container, (void*)3, (void*)4, c);
}
static void teardown(void) {
	stEdgeContainer_destruct(container);
}

static void test_stEdgeContainer_getEdge(CuTest *testCase) {
	setup();
	char *a = (char*)stEdgeContainer_getEdge(container, (void*)1, (void*)2);
	CuAssertStrEquals(testCase, a, "a");
	char *b = (char*)stEdgeContainer_getEdge(container, (void*)1, (void*)3);
	CuAssertStrEquals(testCase, b, "b");
	char *c = (char*)stEdgeContainer_getEdge(container, (void*)3, (void*)4);
	CuAssertStrEquals(testCase, c, "c");


	teardown();
}
static void test_stEdgeContainer_deletionOfThirdEdge(CuTest *testCase) {
	container = stEdgeContainer_construct((void(*)(void*))destruct);
	char *a = st_malloc(2);
	strcpy(a, "a");
	char *b = st_malloc(2);
	strcpy(b, "b");
	void *c = st_malloc(2);
	strcpy(c, "c");
	stEdgeContainer_addEdge(container, (void*)1, (void*)2, a);
	stEdgeContainer_addEdge(container, (void*)1, (void*)3, b);
	stEdgeContainer_addEdge(container, (void*)1, (void*)4, c);
	stEdgeContainer_deleteEdge(container, (void*)1, (void*)3);
	stEdgeContainer_deleteEdge(container, (void*)1, (void*)4);
	CuAssertTrue(testCase, stEdgeContainer_getEdge(container, (void*)1, (void*)2) != NULL);
	teardown();
}


static void test_stEdgeContainer_deleteEdge(CuTest *testCase) {
	setup();
	stEdgeContainer_deleteEdge(container, (void*)1, (void*)2);
	CuAssertTrue(testCase, !stEdgeContainer_getEdge(container, (void*)1, (void*)2));
	char *b = (char*)stEdgeContainer_getEdge(container, (void*)1, (void*)3);
	CuAssertStrEquals(testCase, b, "b");
	teardown();
	setup();
	stEdgeContainer_deleteEdge(container, (void*)1, (void*)3);
	CuAssertTrue(testCase, !stEdgeContainer_getEdge(container, (void*)1, (void*)3));
	char *a = (char*)stEdgeContainer_getEdge(container, (void*)1, (void*)2);
	CuAssertStrEquals(testCase, a, "a");
	teardown();
}
static void test_stEdgeContainer_iterator(CuTest *testCase) {
	setup();
	stEdgeContainerIterator *it = stEdgeContainer_getIterator(container);
	void *node1, *node2;
	while(stEdgeContainer_getNext(it, &node1, &node2)) {
	}
	stEdgeContainer_destructIterator(it);
	teardown();
}

CuSuite *sonLib_stEdgeContainerTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, test_stEdgeContainer_getEdge);
	SUITE_ADD_TEST(suite, test_stEdgeContainer_deleteEdge);
	SUITE_ADD_TEST(suite, test_stEdgeContainer_deletionOfThirdEdge);
	SUITE_ADD_TEST(suite, test_stEdgeContainer_iterator);
    return suite;
}

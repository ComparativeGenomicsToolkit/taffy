#include "sonLibGlobalsTest.h"


stTreap *t;

static void teardown(void);

static void setup(void) {
	t = stTreap_construct((void*)"t");
	stTreap *a = stTreap_construct((void*)"a");
	stTreap *b = stTreap_construct((void*)"b");
	stTreap *c = stTreap_construct((void*)"c");
	stTreap *d = stTreap_construct((void*)"d");
	stTreap *e = stTreap_construct((void*)"e");
	stTreap *f = stTreap_construct((void*)"f");
	stTreap_concat(t, a);
	stTreap_concat(t, b);
	stTreap_concat(t, c);
	stTreap_concat(t, d);
	stTreap_concat(t, e);
	stTreap_concat(t, f);

}
static void teardown(void) {
	stTreap_destruct(t);
}


static void test_stTreap_ordering(CuTest *testCase) {
	setup();

	stTreap *root = stTreap_findRoot(t);
	CuAssertTrue(testCase, root->count == 7);

	char *tour = stTreap_print(t);
	CuAssertStrEquals(testCase, tour, "tabcdef");
	free(tour);

	/*
	   char *reverseTour = stTreap_printBackwards(t);
	   CuAssertStrEquals(testCase, reverseTour, "fedcbat");
	   free(reverseTour);
	   */

	stTreap *temp = stTreap_next(t);
	CuAssertTrue(testCase, stTreap_compare(t, temp) < 0);
	CuAssertTrue(testCase, stTreap_compare(temp, t) > 0);

	stTreap *node = stTreap_findMax(stTreap_findRoot(t));
	CuAssertStrEquals(testCase, node->value, "f");
	while(stTreap_prev(node)) {
		node = stTreap_prev(node);
	}
	CuAssertTrue(testCase, node == t);

	teardown();
}
static void test_stTreap_rotations(CuTest *testCase) {
	setup();
	stTreap *r = stTreap_findRoot(t);
	stTreap_rotateLeft(r);
	stTreap *newRoot = stTreap_findRoot(t);
	CuAssertTrue(testCase, newRoot->left == r);
	CuAssertTrue(testCase, newRoot->count == 7);

	char *newTour = stTreap_print(t);
	CuAssertStrEquals(testCase, newTour, "tabcdef");
	free(newTour);

	stTreap_rotateRight(newRoot);
	CuAssertTrue(testCase, stTreap_findRoot(t) == r);
	char *rotatedTour = stTreap_print(t);
	CuAssertStrEquals(testCase, rotatedTour, "tabcdef");
	free(rotatedTour);


	teardown();
}
static void test_stTreap_split(CuTest *testCase) {
	setup();
	stTreap *splitPoint = stTreap_next(t);
	splitPoint = stTreap_next(splitPoint);
	CuAssertStrEquals(testCase, splitPoint->value, "b");


	stTreap *beforeSplit = stTreap_splitBefore(splitPoint);

	char *beforeSplitStr = stTreap_print(beforeSplit);
	char *splitPointStr = stTreap_print(splitPoint);
	CuAssertStrEquals(testCase, beforeSplitStr, "ta");
	CuAssertStrEquals(testCase, splitPointStr, "bcdef");
	free(beforeSplitStr);
	free(splitPointStr);

	stTreap_concat(beforeSplit, splitPoint);
	char *recombinedTourStr = stTreap_print(t);
	CuAssertStrEquals(testCase, recombinedTourStr, "tabcdef");
	free(recombinedTourStr);

	stTreap *endNode = stTreap_findMax(stTreap_findRoot(t));
	CuAssertStrEquals(testCase, endNode->value, "f");
	char *reverseTourStr = stTreap_printBackwards(endNode);
	CuAssertStrEquals(testCase, reverseTourStr, "fedcbat");
	free(reverseTourStr);


	stTreap *aftersplit = stTreap_next(t);
	CuAssertStrEquals(testCase, aftersplit->value, "a");
	stTreap *after = stTreap_splitAfter(aftersplit);
	CuAssertTrue(testCase, after != NULL);
	char *afterSplitStr = stTreap_print(after);
	CuAssertStrEquals(testCase, afterSplitStr, "bcdef");
	free(afterSplitStr);
	stTreap_destruct(after);


	teardown();
}
static void test_stTreap_heapProperty(CuTest *testCase) {
	setup();


	stTreap *iter = t;
	while(iter) {
		if(iter->parent) CuAssertTrue(testCase, iter->parent->priority > iter->priority);
		iter = stTreap_next(iter);
	}

	teardown();
}

CuSuite *sonLib_stTreapTestSuite(void) {
	CuSuite *suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, test_stTreap_ordering);
	SUITE_ADD_TEST(suite, test_stTreap_rotations);
	SUITE_ADD_TEST(suite, test_stTreap_split);
	SUITE_ADD_TEST(suite, test_stTreap_heapProperty);
	return suite;
}


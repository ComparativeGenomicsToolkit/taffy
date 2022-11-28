/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibContainersTest.c
 *
 *  Created on: 26-May-2010
 *      Author: benedictpaten
 */

#include "sonLibGlobalsTest.h"

static stDoubleTuple *doubleTuple1 = NULL;
static stDoubleTuple *doubleTuple2, *doubleTuple3, *doubleTuple4, *doubleTuple5;

static void teardown() {
    if(doubleTuple1 != NULL) {
        stDoubleTuple_destruct(doubleTuple1);
        stDoubleTuple_destruct(doubleTuple2);
        stDoubleTuple_destruct(doubleTuple3);
        stDoubleTuple_destruct(doubleTuple4);
        stDoubleTuple_destruct(doubleTuple5);
        doubleTuple1 = NULL;
    }
}

static void setup() {
    teardown();
    doubleTuple1 = stDoubleTuple_construct(3, 1.0, 3.0, 2.0);
    doubleTuple2 = stDoubleTuple_construct(4, 1.0, 5.0, 2.0, 7.0);
    doubleTuple3 = stDoubleTuple_construct(0);
    doubleTuple4 = stDoubleTuple_construct(3, 1.0, 3.0, 2.0);
    doubleTuple5 = stDoubleTuple_construct(4, 1.0, 3.0, 2.0, 7.0);
}

static void test_stDoubleTuple_construct(CuTest *testCase) {
    assert(testCase != NULL);
    setup(); //we don't do anything, this just checks that everything can construct and destruct okay.
    teardown();
}

static void test_stDoubleTuple_hashKey(CuTest *testCase) {
    setup();
    stHash *hash = stHash_construct3((uint64_t (*)(const void *))stDoubleTuple_hashKey, (int (*)(const void *, const void *))stDoubleTuple_equalsFn, NULL, NULL);
    stHash_insert(hash, doubleTuple1, doubleTuple1);
    stHash_insert(hash, doubleTuple2, doubleTuple2);
    CuAssertTrue(testCase, stHash_search(hash, doubleTuple1) != NULL);
    CuAssertTrue(testCase, stHash_search(hash, doubleTuple2) != NULL);
    CuAssertTrue(testCase, stHash_search(hash, doubleTuple3) == NULL);
    CuAssertTrue(testCase, stHash_search(hash, doubleTuple4) != NULL);
    stHash_destruct(hash);
    teardown();
}

static void test_stDoubleTuple_cmpFn(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, stDoubleTuple_cmpFn(doubleTuple1, doubleTuple1) == 0);
    CuAssertTrue(testCase, stDoubleTuple_cmpFn(doubleTuple1, doubleTuple2) < 0);
    CuAssertTrue(testCase, stDoubleTuple_cmpFn(doubleTuple1, doubleTuple3) > 0);
    CuAssertTrue(testCase, stDoubleTuple_cmpFn(doubleTuple1, doubleTuple4) == 0);
    CuAssertTrue(testCase, stDoubleTuple_cmpFn(doubleTuple1, doubleTuple5) < 0);

    CuAssertTrue(testCase, stDoubleTuple_cmpFn(doubleTuple2, doubleTuple1) > 0);
    CuAssertTrue(testCase, stDoubleTuple_cmpFn(doubleTuple2, doubleTuple2) == 0);
    CuAssertTrue(testCase, stDoubleTuple_cmpFn(doubleTuple2, doubleTuple3) > 0);
    CuAssertTrue(testCase, stDoubleTuple_cmpFn(doubleTuple2, doubleTuple4) > 0);
    CuAssertTrue(testCase, stDoubleTuple_cmpFn(doubleTuple2, doubleTuple5) > 0);
    CuAssertTrue(testCase, stDoubleTuple_cmpFn(doubleTuple3, doubleTuple3) == 0);
    teardown();
}

static void test_stDoubleTuple_equalsFn(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, stDoubleTuple_equalsFn(doubleTuple1, doubleTuple1));
    CuAssertTrue(testCase, !stDoubleTuple_equalsFn(doubleTuple1, doubleTuple2));
    CuAssertTrue(testCase, !stDoubleTuple_equalsFn(doubleTuple1, doubleTuple3));
    CuAssertTrue(testCase, stDoubleTuple_equalsFn(doubleTuple1, doubleTuple4));
    CuAssertTrue(testCase, !stDoubleTuple_equalsFn(doubleTuple1, doubleTuple5));
    CuAssertTrue(testCase, stDoubleTuple_cmpFn(doubleTuple3, doubleTuple3) == 0);
    teardown();
}

static void test_stDoubleTuple_length(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, stDoubleTuple_length(doubleTuple1) == 3);
    CuAssertTrue(testCase, stDoubleTuple_length(doubleTuple2) == 4);
    CuAssertTrue(testCase, stDoubleTuple_length(doubleTuple3) == 0);
    CuAssertTrue(testCase, stDoubleTuple_length(doubleTuple4) == 3);
    CuAssertTrue(testCase, stDoubleTuple_length(doubleTuple5) == 4);
    teardown();
}

static void test_stDoubleTuple_getPosition(CuTest *testCase) {
    setup();
    CuAssertDblEquals(testCase, stDoubleTuple_getPosition(doubleTuple1, 0), 0.00000, 1);
    CuAssertDblEquals(testCase, stDoubleTuple_getPosition(doubleTuple1, 1), 0.00000, 3);
    CuAssertDblEquals(testCase, stDoubleTuple_getPosition(doubleTuple1, 2), 0.00000, 2);
    stDoubleTuple *i = stDoubleTuple_construct(3, 0.3, 0.2, 1001.02);
    CuAssertDblEquals(testCase, stDoubleTuple_getPosition(i, 0), 0.00000, 0.3);
    CuAssertDblEquals(testCase, stDoubleTuple_getPosition(i, 1), 0.00000, 0.2);
    CuAssertDblEquals(testCase, stDoubleTuple_getPosition(i, 2), 0.00000, 1001.02);
    stDoubleTuple_destruct(i);
    teardown();
}

CuSuite* sonLib_stDoubleTuplesTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stDoubleTuple_construct);
    SUITE_ADD_TEST(suite, test_stDoubleTuple_hashKey);
    SUITE_ADD_TEST(suite, test_stDoubleTuple_cmpFn);
    SUITE_ADD_TEST(suite, test_stDoubleTuple_equalsFn);
    SUITE_ADD_TEST(suite, test_stDoubleTuple_length);
    SUITE_ADD_TEST(suite, test_stDoubleTuple_getPosition);
    return suite;
}

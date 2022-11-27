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

static stIntTuple *intTuple1 = NULL;
static stIntTuple *intTuple2, *intTuple3, *intTuple4, *intTuple5;

static void teardown() {
    if(intTuple1 != NULL) {
        stIntTuple_destruct(intTuple1);
        stIntTuple_destruct(intTuple2);
        stIntTuple_destruct(intTuple3);
        stIntTuple_destruct(intTuple4);
        stIntTuple_destruct(intTuple5);
        intTuple1 = NULL;
    }
}

static void setup() {
    teardown();
    intTuple1 = stIntTuple_construct3( 1, 3, 2);
    intTuple2 = stIntTuple_construct4( 1, 5, 2, 7);
    intTuple3 = stIntTuple_construct0();
    intTuple4 = stIntTuple_construct3( 1, 3, 2);
    intTuple5 = stIntTuple_construct4( 1, 3, 2, 7);
}

static void test_stIntTuple_construct(CuTest *testCase) {
    assert(testCase != NULL);
    setup(); //we don't do anything, this just checks that everything can construct and destruct okay.
    teardown();
}

static void test_stIntTuple_hashKey(CuTest *testCase) {
    setup();
    stHash *hash = stHash_construct3((uint64_t (*)(const void *))stIntTuple_hashKey, (int (*)(const void *, const void *))stIntTuple_equalsFn, NULL, NULL);
    stHash_insert(hash, intTuple1, intTuple1);
    stHash_insert(hash, intTuple2, intTuple2);
    CuAssertTrue(testCase, stHash_search(hash, intTuple1) != NULL);
    CuAssertTrue(testCase, stHash_search(hash, intTuple2) != NULL);
    CuAssertTrue(testCase, stHash_search(hash, intTuple3) == NULL);
    CuAssertTrue(testCase, stHash_search(hash, intTuple4) != NULL);
    stHash_destruct(hash);
    teardown();
}

static void test_stIntTuple_cmpFn(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, stIntTuple_cmpFn(intTuple1, intTuple1) == 0);
    CuAssertTrue(testCase, stIntTuple_cmpFn(intTuple1, intTuple2) < 0);
    CuAssertTrue(testCase, stIntTuple_cmpFn(intTuple1, intTuple3) > 0);
    CuAssertTrue(testCase, stIntTuple_cmpFn(intTuple1, intTuple4) == 0);
    CuAssertTrue(testCase, stIntTuple_cmpFn(intTuple1, intTuple5) < 0);

    CuAssertTrue(testCase, stIntTuple_cmpFn(intTuple2, intTuple1) > 0);
    CuAssertTrue(testCase, stIntTuple_cmpFn(intTuple2, intTuple2) == 0);
    CuAssertTrue(testCase, stIntTuple_cmpFn(intTuple2, intTuple3) > 0);
    CuAssertTrue(testCase, stIntTuple_cmpFn(intTuple2, intTuple4) > 0);
    CuAssertTrue(testCase, stIntTuple_cmpFn(intTuple2, intTuple5) > 0);

    CuAssertTrue(testCase, stIntTuple_cmpFn(intTuple3, intTuple3) == 0);
    teardown();
}

static void test_stIntTuple_equalsFn(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, stIntTuple_equalsFn(intTuple1, intTuple1));
    CuAssertTrue(testCase, !stIntTuple_equalsFn(intTuple1, intTuple2));
    CuAssertTrue(testCase, !stIntTuple_equalsFn(intTuple1, intTuple3));
    CuAssertTrue(testCase, stIntTuple_equalsFn(intTuple1, intTuple4));
    CuAssertTrue(testCase, !stIntTuple_equalsFn(intTuple1, intTuple5));
    CuAssertTrue(testCase, stIntTuple_cmpFn(intTuple3, intTuple3) == 0);
    teardown();
}

static void test_stIntTuple_length(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, stIntTuple_length(intTuple1) == 3);
    CuAssertTrue(testCase, stIntTuple_length(intTuple2) == 4);
    CuAssertTrue(testCase, stIntTuple_length(intTuple3) == 0);
    CuAssertTrue(testCase, stIntTuple_length(intTuple4) == 3);
    CuAssertTrue(testCase, stIntTuple_length(intTuple5) == 4);
    teardown();
}

static void test_stIntTuple_getPosition(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, stIntTuple_get(intTuple1, 0) == 1);
    CuAssertTrue(testCase, stIntTuple_get(intTuple1, 1) == 3);
    CuAssertTrue(testCase, stIntTuple_get(intTuple1, 2) == 2);
    teardown();
}

CuSuite* sonLib_stTuplesTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stIntTuple_construct);
    SUITE_ADD_TEST(suite, test_stIntTuple_hashKey);
    SUITE_ADD_TEST(suite, test_stIntTuple_cmpFn);
    SUITE_ADD_TEST(suite, test_stIntTuple_equalsFn);
    SUITE_ADD_TEST(suite, test_stIntTuple_length);
    SUITE_ADD_TEST(suite, test_stIntTuple_getPosition);
    return suite;
}

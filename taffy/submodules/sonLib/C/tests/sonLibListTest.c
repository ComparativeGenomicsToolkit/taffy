/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibTest.c
 *
 *  Created on: 24 May 2010
 *      Author: benedictpaten
 */

#include "sonLibGlobalsTest.h"

static stList *list = NULL;
static int64_t stringNumber = 5;
static char *strings[5] = { "one", "two", "three", "four", "five" };

static void teardown() {
    if(list != NULL) {
        stList_destruct(list);
        list = NULL;
    }
}

static void setup() {
    teardown();
    list = stList_construct();
    int64_t i;
    for(i=0; i<stringNumber; i++) {
        stList_append(list, strings[i]);
    }
}

void test_stList_construct(CuTest *testCase) {
    setup();
    stList *list2 = stList_construct2(stringNumber);
    CuAssertTrue(testCase, stList_length(list2) == stringNumber);
    int64_t i;
    for(i=0; i<stringNumber; i++) {
        CuAssertTrue(testCase, stList_get(list2, i) == NULL);
    }
    stList_destruct(list2);
    teardown();
}

void test_stList_length(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, stList_length(list) == stringNumber);
    teardown();
}

void test_stList_get(CuTest *testCase) {
    setup();
    int64_t i;
    for(i=0; i<stringNumber; i++) {
        CuAssertTrue(testCase, stList_get(list, i) == strings[i]);
    }
    teardown();
}

void test_stList_set(CuTest *testCase) {
    setup();
    int64_t i;
    for(i=0; i<stringNumber; i++) {
        stList_set(list, i, NULL);
        CuAssertTrue(testCase, stList_get(list, i) == NULL);
    }
    teardown();
}

void test_stList_append(CuTest *testCase) {
    setup();
    stList_append(list, NULL);
    CuAssertTrue(testCase, stList_length(list) == stringNumber+1);
    CuAssertTrue(testCase, stList_get(list, stringNumber) == NULL);
    teardown();
}

void test_stList_appendAll(CuTest *testCase) {
    setup();
    stList *list2 = stList_copy(list, NULL);
    stList_appendAll(list, list2);
    CuAssertTrue(testCase, stList_length(list) == stringNumber * 2);
    int64_t i;
    for(i=0; i<stringNumber*2; i++) {
        CuAssertTrue(testCase, stList_get(list, i) == strings[i % stringNumber]);
    }
    stList_destruct(list2);
    teardown();
}

void test_stList_peek(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, stList_peek(list) == strings[stringNumber-1]);
    teardown();
}

void test_stList_pop(CuTest *testCase) {
    setup();
    int64_t i;
    for(i=stringNumber-1; i>=0; i--) {
        CuAssertTrue(testCase, stList_pop(list) == strings[i]);
        CuAssertTrue(testCase, stList_length(list) == i);
    }
    teardown();
}

void test_stList_remove(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, stList_remove(list, 0) == strings[0]);
    CuAssertTrue(testCase, stList_length(list) == stringNumber-1);
    CuAssertTrue(testCase, stList_remove(list, 1) == strings[2]);
    CuAssertTrue(testCase, stList_length(list) == stringNumber-2);
    teardown();
}

void test_stList_removeItem(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, stList_contains(list, strings[2]));
    stList_removeItem(list, strings[2]);
    CuAssertTrue(testCase, !stList_contains(list, strings[2]));
    stList_removeItem(list, strings[2]);
    CuAssertTrue(testCase, !stList_contains(list, strings[2]));
    teardown();
}

void test_stList_removeFirst(CuTest *testCase) {
    setup();
    int64_t i;
    for(i=0; i<stringNumber; i++) {
        CuAssertTrue(testCase, stList_removeFirst(list) == strings[i]);
        CuAssertTrue(testCase, stList_length(list) == stringNumber-1-i);
    }
    teardown();
}

void test_stList_removeInterval(CuTest *testCase) {
    // Tests all possible remove interval ops on the list
    for(int64_t i=0; i<stringNumber; i++) {
        for(int64_t j=i; j<stringNumber; j++) {
            setup();
            stList_removeInterval(list, i, j-i);
            CuAssertTrue(testCase, stList_length(list) == stringNumber - (j-i));
            for(int64_t k=0; k<i; k++) {
                CuAssertTrue(testCase, stList_get(list, k) == strings[k]);
            }
            for(int64_t k=j; k<stringNumber; k++) {
                CuAssertTrue(testCase, stList_get(list, k-(j-i)) == strings[k]);
            }
            teardown();
        }
    }
}

void test_stList_contains(CuTest *testCase) {
    setup();
    int64_t i;
    for (i = 0; i < stringNumber; i++) {
        CuAssertTrue(testCase, stList_contains(list, strings[i]));
    }
    CuAssertTrue(testCase, !stList_contains(list, "something"));
    CuAssertTrue(testCase, !stList_contains(list, NULL));
    teardown();
}

void test_stList_find(CuTest *testCase) {
    setup();
    int64_t i;
    for (i = 0; i < stringNumber; i++) {
        CuAssertTrue(testCase, stList_find(list, strings[i]) == i);
    }
    CuAssertTrue(testCase, stList_find(list, "something") == -1);
    CuAssertTrue(testCase, stList_find(list, NULL) == -1);
    teardown();
}

void test_stList_copy(CuTest *testCase) {
    setup();
    stList *list2 = stList_copy(list, NULL);
    CuAssertTrue(testCase, stList_length(list) == stList_length(list2));
    int64_t i;
    for (i = 0; i < stringNumber; i++) {
        CuAssertTrue(testCase, stList_get(list2, i) == strings[i]);
    }
    stList_destruct(list2);
    teardown();
}

void test_stList_reverse(CuTest *testCase) {
    setup();
    stList_reverse(list);
    int64_t i;
    for(i=0; i<stringNumber; i++) {
        CuAssertTrue(testCase, stList_get(list, i) == strings[stringNumber-1-i]);
    }
    teardown();
}

void test_stList_iterator(CuTest *testCase) {
    setup();
    stListIterator *it = stList_getIterator(list);
    int64_t i;
    for(i=0; i<stringNumber; i++) {
        CuAssertTrue(testCase, stList_getNext(it) == strings[i]);
    }
    CuAssertTrue(testCase, stList_getNext(it) == NULL);
    CuAssertTrue(testCase, stList_getNext(it) == NULL);
    stListIterator *it2 = stList_copyIterator(it);
    for(i=0; i<stringNumber; i++) {
        CuAssertTrue(testCase, stList_getPrevious(it) == strings[stringNumber-1-i]);
        CuAssertTrue(testCase, stList_getPrevious(it2) == strings[stringNumber-1-i]);
    }
    CuAssertTrue(testCase, stList_getPrevious(it) == NULL);
    CuAssertTrue(testCase, stList_getPrevious(it) == NULL);
    CuAssertTrue(testCase, stList_getPrevious(it2) == NULL);
    stList_destructIterator(it);
    stList_destructIterator(it2);
    teardown();
}

void test_stList_nullList(CuTest *testCase) {
    stList *empty = NULL;
    CuAssertTrue(testCase, stList_length(empty) == 0);
    stListIterator *it = stList_getIterator(empty);
    CuAssertTrue(testCase, stList_getNext(it) == NULL);
    CuAssertTrue(testCase, stList_getPrevious(it) == NULL);
    stList_destructIterator(it);
    stList_destruct(empty);
}

void test_stList_sort(CuTest *testCase) {
    setup();
    stList_sort(list, (int (*)(const void *, const void *))strcmp);
    CuAssertTrue(testCase, stList_length(list) == stringNumber);
    CuAssertStrEquals(testCase, "five", stList_get(list, 0));
    CuAssertStrEquals(testCase, "four", stList_get(list, 1));
    CuAssertStrEquals(testCase, "one", stList_get(list, 2));
    CuAssertStrEquals(testCase, "three", stList_get(list, 3));
    CuAssertStrEquals(testCase, "two", stList_get(list, 4));
    teardown();
}

int test_stList_sort2P(const void *a, const void *b, void *extraArg) {
    assert(strcmp(extraArg, "hello") == 0);
    return strcmp(a, b);
}

void test_stList_sort2(CuTest *testCase) {
    setup();
    stList_sort2(list, test_stList_sort2P, "hello");
    CuAssertTrue(testCase, stList_length(list) == stringNumber);
    CuAssertStrEquals(testCase, "five", stList_get(list, 0));
    CuAssertStrEquals(testCase, "four", stList_get(list, 1));
    CuAssertStrEquals(testCase, "one", stList_get(list, 2));
    CuAssertStrEquals(testCase, "three", stList_get(list, 3));
    CuAssertStrEquals(testCase, "two", stList_get(list, 4));
    teardown();
}

void test_stList_getSortedSet(CuTest *testCase) {
    setup();
    stSortedSet *sortedSet = stList_getSortedSet(list, (int (*)(const void *, const void *))strcmp);
    CuAssertTrue(testCase, stSortedSet_size(sortedSet) == stringNumber);
    stSortedSetIterator *iterator = stSortedSet_getIterator(sortedSet);
    CuAssertStrEquals(testCase, "five", stSortedSet_getNext(iterator));
    CuAssertStrEquals(testCase, "four", stSortedSet_getNext(iterator));
    CuAssertStrEquals(testCase, "one", stSortedSet_getNext(iterator));
    CuAssertStrEquals(testCase, "three", stSortedSet_getNext(iterator));
    CuAssertStrEquals(testCase, "two", stSortedSet_getNext(iterator));
    stSortedSet_destructIterator(iterator);
    stSortedSet_destruct(sortedSet);
    teardown();
}

void test_stList_filter(CuTest *testCase) {
    setup();
    stSortedSet *set = stSortedSet_construct();
    stSortedSet_insert(set, strings[0]);
    stSortedSet_insert(set, strings[4]);
    stList *list2 = stList_filterToExclude(list, set);
    stList *list3 = stList_filterToInclude(list, set);
    CuAssertTrue(testCase,stList_length(list2) == 3);
    CuAssertTrue(testCase,stList_length(list3) == 2);
    CuAssertTrue(testCase,stList_get(list2, 0) == strings[1]);
    CuAssertTrue(testCase,stList_get(list2, 1) == strings[2]);
    CuAssertTrue(testCase,stList_get(list2, 2) == strings[3]);
    CuAssertTrue(testCase,stList_get(list3, 0) == strings[0]);
    CuAssertTrue(testCase,stList_get(list3, 1) == strings[4]);

    stSortedSet_destruct(set);
    stList_destruct(list2);
    stList_destruct(list3);
    teardown();
}

CuSuite* sonLib_stListTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stList_construct);
    SUITE_ADD_TEST(suite, test_stList_append);
    SUITE_ADD_TEST(suite, test_stList_appendAll);
    SUITE_ADD_TEST(suite, test_stList_length);
    SUITE_ADD_TEST(suite, test_stList_get);
    SUITE_ADD_TEST(suite, test_stList_set);
    SUITE_ADD_TEST(suite, test_stList_peek);
    SUITE_ADD_TEST(suite, test_stList_pop);
    SUITE_ADD_TEST(suite, test_stList_remove);
    SUITE_ADD_TEST(suite, test_stList_removeFirst);
    SUITE_ADD_TEST(suite, test_stList_removeInterval);
    SUITE_ADD_TEST(suite, test_stList_contains);
    SUITE_ADD_TEST(suite, test_stList_find);
    SUITE_ADD_TEST(suite, test_stList_copy);
    SUITE_ADD_TEST(suite, test_stList_reverse);
    SUITE_ADD_TEST(suite, test_stList_iterator);
    SUITE_ADD_TEST(suite, test_stList_nullList);
    SUITE_ADD_TEST(suite, test_stList_sort);
    SUITE_ADD_TEST(suite, test_stList_getSortedSet);
    SUITE_ADD_TEST(suite, test_stList_filter);
    return suite;
}

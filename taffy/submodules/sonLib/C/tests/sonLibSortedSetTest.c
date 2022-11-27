/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLibGlobalsTest.h"

static stSortedSet *sortedSet = NULL;
static stSortedSet *sortedSet2 = NULL;
static int64_t size = 7;
static int64_t input[] = { 1, 5, -1, 10, 12, 3, -10 };
static int64_t sortedInput[] = { -10, -1, 1, 3, 5, 10, 12 };
static int64_t sortedSize = 7;


static void sonLibSortedSetTestTeardown() {
    if(sortedSet != NULL) {
        stSortedSet_destruct(sortedSet);
        stSortedSet_destruct(sortedSet2);
        sortedSet = NULL;
    }
}

static void sonLibSortedSetTestSetup() {
    sonLibSortedSetTestTeardown();
    sortedSet = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn,
            (void (*)(void *))stIntTuple_destruct);
    sortedSet2 = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn,
                (void (*)(void *))stIntTuple_destruct);
}

static void test_stSortedSet_construct(CuTest* testCase) {
    sonLibSortedSetTestSetup();
    CuAssertTrue(testCase, sortedSet != NULL);
    sonLibSortedSetTestTeardown();
}

static void test_stSortedSet_copyConstruct(CuTest* testCase) {
    sonLibSortedSetTestSetup();
    CuAssertTrue(testCase, sortedSet != NULL);
    int64_t i;
    for(i=0; i<size; i++) {
        stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
    }
    stSortedSet *sortedSet2 = stSortedSet_copyConstruct(sortedSet, NULL);
    CuAssertTrue(testCase, stSortedSet_size(sortedSet2) == stSortedSet_size(sortedSet));
    CuAssertTrue(testCase, stSortedSet_equals(sortedSet2, sortedSet));
    stSortedSet_destruct(sortedSet2);
    sonLibSortedSetTestTeardown();
}

static void test_stSortedSet(CuTest* testCase) {
    sonLibSortedSetTestSetup();
    int64_t i;
    CuAssertIntEquals(testCase, 0, stSortedSet_size(sortedSet));
    stIntTuple *oldIntTuples[size];
    for(i=0; i<size; i++) {
        oldIntTuples[i] = stIntTuple_construct1( input[i]);
        stSortedSet_insert(sortedSet, oldIntTuples[i]);
    }

    // Add in the input one more time to verify the information isn't duplicated
    for(i=0; i<size; i++) {
        stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
        stIntTuple_destruct(oldIntTuples[i]);
    }

    CuAssertIntEquals(testCase, sortedSize, stSortedSet_size(sortedSet));
    CuAssertIntEquals(testCase, sortedInput[0], stIntTuple_get(stSortedSet_getFirst(sortedSet), 0));
    CuAssertIntEquals(testCase, sortedInput[sortedSize-1], stIntTuple_get(stSortedSet_getLast(sortedSet), 0));
    for(i=0; i<sortedSize; i++) {
        CuAssertIntEquals(testCase, sortedSize-i, stSortedSet_size(sortedSet));
        stIntTuple *tuple = stIntTuple_construct1( sortedInput[i]);
        CuAssertTrue(testCase, stIntTuple_get(stSortedSet_search(sortedSet, tuple), 0) == sortedInput[i]);
        stIntTuple_destruct(stSortedSet_remove(sortedSet, tuple));
        CuAssertTrue(testCase, stSortedSet_search(sortedSet, tuple) == NULL);
        stIntTuple_destruct(tuple);
    }
    sonLibSortedSetTestTeardown();
}

static void test_stSortedSetIterator(CuTest* testCase) {
    sonLibSortedSetTestSetup();
    int64_t i;
    for(i=0; i<size; i++) {
        stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
    }
    stSortedSetIterator *iterator = stSortedSet_getIterator(sortedSet);
    CuAssertTrue(testCase, iterator != NULL);

    for(i=0; i<sortedSize; i++) {
        CuAssertIntEquals(testCase, sortedInput[i], stIntTuple_get(stSortedSet_getNext(iterator), 0));
    }
    CuAssertTrue(testCase, stSortedSet_getNext(iterator) == NULL);
    stSortedSetIterator *iterator2 = stSortedSet_copyIterator(iterator);
    CuAssertTrue(testCase, iterator2 != NULL);
    for(i=0; i<sortedSize; i++) {
        CuAssertIntEquals(testCase, sortedInput[sortedSize - 1 - i], stIntTuple_get(stSortedSet_getPrevious(iterator), 0));
        CuAssertIntEquals(testCase, sortedInput[sortedSize - 1 - i], stIntTuple_get(stSortedSet_getPrevious(iterator2), 0));
    }
    CuAssertTrue(testCase, stSortedSet_getPrevious(iterator) == NULL);
    CuAssertTrue(testCase, stSortedSet_getPrevious(iterator2) == NULL);
    stSortedSet_destructIterator(iterator);
    stSortedSet_destructIterator(iterator2);
    sonLibSortedSetTestTeardown();
}

static void test_stSortedSetIterator_getIteratorFrom(CuTest* testCase) {
    sonLibSortedSetTestSetup();
    int64_t i;
    for(i=0; i<size; i++) {
        stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
    }
    stSortedSetIterator *iterator = stSortedSet_getIterator(sortedSet);
    CuAssertTrue(testCase, iterator != NULL);
    stSortedSet_destructIterator(iterator);

    for(i=0; i<sortedSize; i++) {
        stIntTuple *queryTuple = stIntTuple_construct1( sortedInput[i]);
        stSortedSetIterator *it = stSortedSet_getIteratorFrom(sortedSet, queryTuple);
        stIntTuple *intTuple = stSortedSet_getNext(it);
        CuAssertTrue(testCase, intTuple != NULL);
        CuAssertIntEquals(testCase, sortedInput[i], stIntTuple_get(intTuple, 0));
        stSortedSet_destructIterator(it);
        stIntTuple_destruct(queryTuple);
    }

    stIntTuple *seven = stIntTuple_construct1( 7);
    stTry {
        stSortedSet_getIteratorFrom(sortedSet, seven); //This number if not in the input.
        CuAssertTrue(testCase, 0);
    } stCatch(except) {
        CuAssertTrue(testCase, stExcept_getId(except) == SORTED_SET_EXCEPTION_ID);
    } stTryEnd
    stIntTuple_destruct(seven);

    sonLibSortedSetTestTeardown();
}

static void test_stSortedSetIterator_getReverseIterator(CuTest* testCase) {
    sonLibSortedSetTestSetup();
    int64_t i;
    for(i=0; i<size; i++) {
        stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
    }
    stSortedSetIterator *iterator = stSortedSet_getIterator(sortedSet);
    CuAssertTrue(testCase, iterator != NULL);
    stSortedSet_destructIterator(iterator);

    stSortedSetIterator *it = stSortedSet_getReverseIterator(sortedSet);
    stIntTuple *intTuple = stSortedSet_getPrevious(it);
    CuAssertTrue(testCase, intTuple != NULL);
    CuAssertIntEquals(testCase, sortedInput[sortedSize-1], stIntTuple_get(intTuple, 0));
    stSortedSet_destructIterator(it);


    sonLibSortedSetTestTeardown();
}


static void test_stSortedSetEquals(CuTest* testCase) {
    sonLibSortedSetTestSetup();
    CuAssertTrue(testCase, stSortedSet_equals(sortedSet, sortedSet));
    int64_t i;
    for(i=0; i<size; i++) {
         stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
    }
    CuAssertTrue(testCase, stSortedSet_equals(sortedSet, sortedSet));

    CuAssertTrue(testCase, !stSortedSet_equals(sortedSet, sortedSet2));
    for(i=1; i<size; i++) { //first argument is unique in input..
       stSortedSet_insert(sortedSet2, stIntTuple_construct1( input[i]));
    }
    CuAssertTrue(testCase, !stSortedSet_equals(sortedSet, sortedSet2));

    stSortedSet_insert(sortedSet2, stIntTuple_construct1( input[0]));
    CuAssertTrue(testCase, stSortedSet_equals(sortedSet, sortedSet2));

    stSortedSet *sortedSet3 = stSortedSet_construct2((void (*)(void *)) stIntTuple_destruct); //diff comparator..
    CuAssertTrue(testCase, !stSortedSet_equals(sortedSet, sortedSet3));
    for(i=0; i<size; i++) {
        stSortedSet_insert(sortedSet3, stIntTuple_construct1( input[i]));
    }
    CuAssertTrue(testCase, !stSortedSet_equals(sortedSet, sortedSet3));
    stSortedSet_destruct(sortedSet3);

    sonLibSortedSetTestTeardown();
}

static void test_stSortedSetIntersection(CuTest* testCase) {
    sonLibSortedSetTestSetup();
    //Check intersection of empty sets is okay..
    stSortedSet *sortedSet3 = stSortedSet_getIntersection(sortedSet, sortedSet2);
    CuAssertTrue(testCase, stSortedSet_size(sortedSet3) == 0);
    stSortedSet_destruct(sortedSet3);

    int64_t i;
    for(i=0; i<size; i++) {
        stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
    }

    //Check intersection of empty and non-empty set is empty.
    sortedSet3 = stSortedSet_getIntersection(sortedSet, sortedSet2);
    CuAssertTrue(testCase, stSortedSet_size(sortedSet3) == 0);
    stSortedSet_destruct(sortedSet3);

    //Check intersection of two non-empty, overlapping sets in correct.
    stSortedSet_insert(sortedSet2, stIntTuple_construct1( 0));
    stSortedSet_insert(sortedSet2, stIntTuple_construct1( 1));
    stSortedSet_insert(sortedSet2, stIntTuple_construct1( 5));

    sortedSet3 = stSortedSet_getIntersection(sortedSet, sortedSet2);
    CuAssertTrue(testCase, stSortedSet_size(sortedSet3) == 2);
    stIntTuple *intTuple = stIntTuple_construct1( 1);
    CuAssertTrue(testCase, stSortedSet_search(sortedSet3, intTuple) != NULL);
    stIntTuple_destruct(intTuple);
    intTuple = stIntTuple_construct1( 5);
    CuAssertTrue(testCase, stSortedSet_search(sortedSet3, intTuple) != NULL);
    stIntTuple_destruct(intTuple);
    stSortedSet_destruct(sortedSet3);

    //Check we get an exception with sorted sets with different comparators.
    stSortedSet *sortedSet4 = stSortedSet_construct();
    stTry {
        stSortedSet_getIntersection(sortedSet, sortedSet4);
    } stCatch(except) {
        CuAssertTrue(testCase, stExcept_getId(except) == SORTED_SET_EXCEPTION_ID);
    } stTryEnd
    stSortedSet_destruct(sortedSet4);

    sonLibSortedSetTestTeardown();
}

static void test_stSortedSetUnion(CuTest* testCase) {
    sonLibSortedSetTestSetup();
    //Check intersection of empty sets is okay..
    stSortedSet *sortedSet3 = stSortedSet_getUnion(sortedSet, sortedSet2);
    CuAssertTrue(testCase, stSortedSet_size(sortedSet3) == 0);
    stSortedSet_destruct(sortedSet3);

    int64_t i;
    for(i=0; i<size; i++) {
        stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
    }

    //Check intersection of empty and non-empty set is equal to the non-empty set.
    sortedSet3 = stSortedSet_getUnion(sortedSet, sortedSet2);
    CuAssertTrue(testCase, stSortedSet_equals(sortedSet, sortedSet3));
    stSortedSet_destruct(sortedSet3);

    //Check intersection of two non-empty, overlapping sets in correct.
    stSortedSet_insert(sortedSet2, stIntTuple_construct1( 0));
    stSortedSet_insert(sortedSet2, stIntTuple_construct1( 1));
    stSortedSet_insert(sortedSet2, stIntTuple_construct1( 5));

    sortedSet3 = stSortedSet_getUnion(sortedSet, sortedSet2);
    CuAssertTrue(testCase, stSortedSet_size(sortedSet3) == stSortedSet_size(sortedSet) + 1);
    CuAssertTrue(testCase, !stSortedSet_equals(sortedSet, sortedSet3));
    stSortedSet_insert(sortedSet, stIntTuple_construct1( 0));
    CuAssertTrue(testCase, stSortedSet_size(sortedSet3) == stSortedSet_size(sortedSet));
    CuAssertTrue(testCase, stSortedSet_equals(sortedSet, sortedSet3));
    stSortedSet_destruct(sortedSet3);

    //Check we get an exception with sorted sets with different comparators.
    stSortedSet *sortedSet4 = stSortedSet_construct();
    stTry {
        stSortedSet_getUnion(sortedSet, sortedSet4);
        CuAssertTrue(testCase, 0);
    } stCatch(except) {
        CuAssertTrue(testCase, stExcept_getId(except) == SORTED_SET_EXCEPTION_ID);
    } stTryEnd
    stSortedSet_destruct(sortedSet4);

    sonLibSortedSetTestTeardown();
}

static void test_stSortedSetDifference(CuTest* testCase) {
    sonLibSortedSetTestSetup();

    //Check difference of empty sets is okay..
    stSortedSet *sortedSet3 = stSortedSet_getDifference(sortedSet, sortedSet2);
    CuAssertTrue(testCase, stSortedSet_size(sortedSet3) == 0);
    stSortedSet_destruct(sortedSet3);

    int64_t i;
    for(i=0; i<size; i++) {
        stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
    }

    //Check difference of non-empty set / empty set is the non-empty.
    sortedSet3 = stSortedSet_getDifference(sortedSet, sortedSet2);
    CuAssertTrue(testCase, stSortedSet_equals(sortedSet, sortedSet3));
    stSortedSet_destruct(sortedSet3);

    //Check difference of two non-empty, overlapping sets in correct.
    stSortedSet_insert(sortedSet2, stIntTuple_construct1( 0));
    stSortedSet_insert(sortedSet2, stIntTuple_construct1( 1));
    stSortedSet_insert(sortedSet2, stIntTuple_construct1( 5));

    sortedSet3 = stSortedSet_getDifference(sortedSet, sortedSet2);
    CuAssertTrue(testCase, stSortedSet_size(sortedSet3) == stSortedSet_size(sortedSet) - 2);
    CuAssertTrue(testCase, !stSortedSet_equals(sortedSet, sortedSet3));
    stIntTuple *one = stIntTuple_construct1( 1);
    stSortedSet_insert(sortedSet3, one);
    stIntTuple *five = stIntTuple_construct1( 5);
    stSortedSet_insert(sortedSet3, five);
    CuAssertTrue(testCase, stSortedSet_equals(sortedSet, sortedSet3));
    stSortedSet_destruct(sortedSet3);
    stIntTuple_destruct(one);
    stIntTuple_destruct(five);

    //Check we get an exception when merging sorted sets with different comparators.
    stSortedSet *sortedSet4 = stSortedSet_construct();
    stTry {
       stSortedSet_getDifference(sortedSet, sortedSet4);
       CuAssertTrue(testCase, 0);
    } stCatch(except) {
        CuAssertTrue(testCase, stExcept_getId(except) == SORTED_SET_EXCEPTION_ID);
    } stTryEnd
    stSortedSet_destruct(sortedSet4);

    sonLibSortedSetTestTeardown();
}

static void *querySortedIntTupleSet(stSortedSet *sortedSet, int64_t n) {
    stIntTuple *nTuple = stIntTuple_construct1(n);
    void *result = stSortedSet_search(sortedSet, nTuple);
    stIntTuple_destruct(nTuple);
    return result;
}

static void *querySortedIntTupleSetLessThanOrEqual(stSortedSet *sortedSet, int64_t n) {
    stIntTuple *nTuple = stIntTuple_construct1(n);
    void *result = stSortedSet_searchLessThanOrEqual(sortedSet, nTuple);
    stIntTuple_destruct(nTuple);
    return result;
}

static void test_stSortedSet_searchLessThanOrEqual(CuTest* testCase) {
    sonLibSortedSetTestSetup();
    for(int64_t i=0; i<size; i++) {
        stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
    }
    //static int64_t sortedInput[] = { -10, -1, 1, 3, 5, 10, 12 };
    CuAssertTrue(testCase, querySortedIntTupleSetLessThanOrEqual(sortedSet, -11) == NULL);
    CuAssertTrue(testCase, querySortedIntTupleSetLessThanOrEqual(sortedSet, -10) == querySortedIntTupleSet(sortedSet, -10));
    CuAssertTrue(testCase, querySortedIntTupleSetLessThanOrEqual(sortedSet, -5) == querySortedIntTupleSet(sortedSet, -10));
    CuAssertTrue(testCase, querySortedIntTupleSetLessThanOrEqual(sortedSet, 1) == querySortedIntTupleSet(sortedSet, 1));
    CuAssertTrue(testCase, querySortedIntTupleSetLessThanOrEqual(sortedSet, 13) == querySortedIntTupleSet(sortedSet, 12));
    CuAssertTrue(testCase, querySortedIntTupleSetLessThanOrEqual(sortedSet, 7) == querySortedIntTupleSet(sortedSet, 5));

    for(int64_t i=0; i<100; i++) {
        stIntTuple *k = stIntTuple_construct1(st_randomInt(-1000, 1000));
        if (stSortedSet_search(sortedSet, k)) {
            stIntTuple_destruct(k);
        } else {
            stSortedSet_insert(sortedSet, k);
        }
    }
    stList *list = stSortedSet_getList(sortedSet);
    for(int64_t i=1; i<stList_length(list); i++) {
        stIntTuple *p = stList_get(list, i-1);
        stIntTuple *j = stList_get(list, i);
        stIntTuple *k = stIntTuple_construct1( st_randomInt(stIntTuple_get(p, 0), stIntTuple_get(j, 0)));
        CuAssertTrue(testCase, stSortedSet_searchLessThanOrEqual(sortedSet, k) == p);
        stIntTuple_destruct(k);
    }
    stList_destruct(list);

    sonLibSortedSetTestTeardown();
}

static void *querySortedIntTupleSetLessThan(stSortedSet *sortedSet, int64_t n) {
    stIntTuple *nTuple = stIntTuple_construct1(n);
    void *result = stSortedSet_searchLessThan(sortedSet, nTuple);
    stIntTuple_destruct(nTuple);
    return result;
}

static void test_stSortedSet_searchLessThan(CuTest* testCase) {
    sonLibSortedSetTestSetup();
    for(int64_t i=0; i<size; i++) {
        stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
    }
    //static int64_t sortedInput[] = { -10, -1, 1, 3, 5, 10, 12 };
    CuAssertTrue(testCase, querySortedIntTupleSetLessThan(sortedSet, -11) == NULL);
    CuAssertTrue(testCase, querySortedIntTupleSetLessThan(sortedSet, -10) == NULL);
    CuAssertTrue(testCase, querySortedIntTupleSetLessThan(sortedSet, -5) == querySortedIntTupleSet(sortedSet, -10));
    CuAssertTrue(testCase, querySortedIntTupleSetLessThan(sortedSet, 13) == querySortedIntTupleSet(sortedSet, 12));
    CuAssertTrue(testCase, querySortedIntTupleSetLessThan(sortedSet, 8) == querySortedIntTupleSet(sortedSet, 5));
    CuAssertTrue(testCase, querySortedIntTupleSetLessThan(sortedSet, 1) == querySortedIntTupleSet(sortedSet, -1));
    CuAssertTrue(testCase, querySortedIntTupleSetLessThan(sortedSet, 10) == querySortedIntTupleSet(sortedSet, 5));

    for(int64_t i=0; i<100; i++) {
        stIntTuple *k = stIntTuple_construct1(st_randomInt(-1000, 1000));
        if (stSortedSet_search(sortedSet, k)) {
            stIntTuple_destruct(k);
        } else {
            stSortedSet_insert(sortedSet, k);
        }
    }
    stList *list = stSortedSet_getList(sortedSet);
    for(int64_t i=1; i<stList_length(list); i++) {
        stIntTuple *p = stList_get(list, i-1);
        stIntTuple *j = stList_get(list, i);
        stIntTuple *k = stIntTuple_construct1( st_randomInt(stIntTuple_get(p, 0)+1, stIntTuple_get(j, 0)+1));
        CuAssertTrue(testCase, stSortedSet_searchLessThan(sortedSet, k) == p);
        stIntTuple_destruct(k);
    }
    stList_destruct(list);

    sonLibSortedSetTestTeardown();
}

static void *querySortedIntTupleSetGreaterThanOrEqual(stSortedSet *sortedSet, int64_t n) {
    stIntTuple *nTuple = stIntTuple_construct1(n);
    void *result = stSortedSet_searchGreaterThanOrEqual(sortedSet, nTuple);
    stIntTuple_destruct(nTuple);
    return result;
}

static void test_stSortedSet_searchGreaterThanOrEqual(CuTest* testCase) {
    sonLibSortedSetTestSetup();

    for(int64_t i=0; i<size; i++) {
        stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
    }
    //static int64_t sortedInput[] = { -10, -1, 1, 3, 5, 10, 12 };
    CuAssertTrue(testCase, querySortedIntTupleSetGreaterThanOrEqual(sortedSet, -11) == querySortedIntTupleSet(sortedSet, -10));
    CuAssertTrue(testCase, querySortedIntTupleSetGreaterThanOrEqual(sortedSet, -10) == querySortedIntTupleSet(sortedSet, -10));
    CuAssertTrue(testCase, querySortedIntTupleSetGreaterThanOrEqual(sortedSet, -5) == querySortedIntTupleSet(sortedSet, -1));
    CuAssertTrue(testCase, querySortedIntTupleSetGreaterThanOrEqual(sortedSet, 1) == querySortedIntTupleSet(sortedSet, 1));
    CuAssertTrue(testCase, querySortedIntTupleSetGreaterThanOrEqual(sortedSet, 13) == NULL);
    CuAssertTrue(testCase, querySortedIntTupleSetGreaterThanOrEqual(sortedSet, 5) == querySortedIntTupleSet(sortedSet, 5));

    for(int64_t i=0; i<100; i++) {
        stIntTuple *k = stIntTuple_construct1(st_randomInt(-1000, 1000));
        if (stSortedSet_search(sortedSet, k)) {
            stIntTuple_destruct(k);
        } else {
            stSortedSet_insert(sortedSet, k);
        }
    }
    stList *list = stSortedSet_getList(sortedSet);
    for(int64_t i=1; i<stList_length(list); i++) {
        stIntTuple *p = stList_get(list, i-1);
        stIntTuple *j = stList_get(list, i);
        stIntTuple *k = stIntTuple_construct1( st_randomInt(stIntTuple_get(p, 0)+1, stIntTuple_get(j, 0)+1));
        CuAssertTrue(testCase, stSortedSet_searchGreaterThanOrEqual(sortedSet, k) == j);
        stIntTuple_destruct(k);
    }
    stList_destruct(list);

    sonLibSortedSetTestTeardown();
}

static void *querySortedIntTupleSetGreaterThan(stSortedSet *sortedSet, int64_t n) {
    stIntTuple *nTuple = stIntTuple_construct1(n);
    void *result = stSortedSet_searchGreaterThan(sortedSet, nTuple);
    stIntTuple_destruct(nTuple);
    return result;
}

static void test_stSortedSet_searchGreaterThan(CuTest* testCase) {
    sonLibSortedSetTestSetup();

    for(int64_t i=0; i<size; i++) {
        stSortedSet_insert(sortedSet, stIntTuple_construct1( input[i]));
    }
    //static int64_t sortedInput[] = { -10, -1, 1, 3, 5, 10, 12 };
    CuAssertTrue(testCase, querySortedIntTupleSetGreaterThan(sortedSet, -11) == querySortedIntTupleSet(sortedSet, -10));
    CuAssertTrue(testCase, querySortedIntTupleSetGreaterThan(sortedSet, -10) == querySortedIntTupleSet(sortedSet, -1));
    CuAssertTrue(testCase, querySortedIntTupleSetGreaterThan(sortedSet, -5) == querySortedIntTupleSet(sortedSet, -1));
    CuAssertTrue(testCase, querySortedIntTupleSetGreaterThan(sortedSet, 1) == querySortedIntTupleSet(sortedSet, 3));
    CuAssertTrue(testCase, querySortedIntTupleSetGreaterThan(sortedSet, 13) == NULL);
    CuAssertTrue(testCase, querySortedIntTupleSetGreaterThan(sortedSet, 12) == NULL);

    for(int64_t i=0; i<100; i++) {
        stIntTuple *k = stIntTuple_construct1(st_randomInt(-1000, 1000));
        if (stSortedSet_search(sortedSet, k)) {
            stIntTuple_destruct(k);
        } else {
            stSortedSet_insert(sortedSet, k);
        }
    }
    stList *list = stSortedSet_getList(sortedSet);
    for(int64_t i=1; i<stList_length(list); i++) {
        stIntTuple *p = stList_get(list, i-1);
        stIntTuple *j = stList_get(list, i);
        stIntTuple *k = stIntTuple_construct1( st_randomInt(stIntTuple_get(p, 0), stIntTuple_get(j, 0)));
        CuAssertTrue(testCase, stSortedSet_searchGreaterThan(sortedSet, k) == j);
        stIntTuple_destruct(k);
    }
    stList_destruct(list);

    sonLibSortedSetTestTeardown();
}

CuSuite* sonLib_stSortedSetTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stSortedSet_construct);
    SUITE_ADD_TEST(suite, test_stSortedSet_copyConstruct);
    SUITE_ADD_TEST(suite, test_stSortedSet);
    SUITE_ADD_TEST(suite, test_stSortedSetIterator);
    SUITE_ADD_TEST(suite, test_stSortedSetIterator_getIteratorFrom);
    SUITE_ADD_TEST(suite, test_stSortedSetIterator_getReverseIterator);
    SUITE_ADD_TEST(suite, test_stSortedSetEquals);
    SUITE_ADD_TEST(suite, test_stSortedSetIntersection);
    SUITE_ADD_TEST(suite, test_stSortedSetUnion);
    SUITE_ADD_TEST(suite, test_stSortedSetDifference);
    SUITE_ADD_TEST(suite, test_stSortedSet_searchLessThanOrEqual);
    SUITE_ADD_TEST(suite, test_stSortedSet_searchLessThan);
    SUITE_ADD_TEST(suite, test_stSortedSet_searchGreaterThanOrEqual);
    SUITE_ADD_TEST(suite, test_stSortedSet_searchGreaterThan);
    return suite;
}

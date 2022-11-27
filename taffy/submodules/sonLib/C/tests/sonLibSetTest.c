/*
 * Copyright (C) 2012 by Dent Earl dentearl (a) gmail com
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLibGlobalsTest.h"

static stSet *set0, *set0Prime;
static stSet *set1;
static stIntTuple *one, *two, *three, *four, *five, *six;

static void testSetup() {
    // compare by value of memory address
    set0 = stSet_construct();
    set0Prime = stSet_construct();
    // compare by value of ints.
    set1 = stSet_construct3((uint64_t(*)(const void *)) stIntTuple_hashKey, 
                            (int(*)(const void *, const void *)) stIntTuple_equalsFn,
                            (void(*)(void *)) stIntTuple_destruct);
    one = stIntTuple_construct1( 0);
    two = stIntTuple_construct1( 1);
    three = stIntTuple_construct1( 2);
    four = stIntTuple_construct1( 3);
    five = stIntTuple_construct1( 4);
    six = stIntTuple_construct1( 5);
    stSet_insert(set0, one);
    stSet_insert(set0, two);
    stSet_insert(set0, three);
    stSet_insert(set0, four);
    stSet_insert(set0, five);
    stSet_insert(set0, six);
    stSet_insert(set0Prime, one);
    stSet_insert(set0Prime, two);
    stSet_insert(set0Prime, three);
    stSet_insert(set0Prime, four);
    stSet_insert(set0Prime, five);
    stSet_insert(set0Prime, six);
    stSet_insert(set1, one);
    stSet_insert(set1, two);
    stSet_insert(set1, three);
    stSet_insert(set1, four);
    stSet_insert(set1, five);
    stSet_insert(set1, six);
}
static void testTeardown() {
    stSet_destruct(set0);
    stSet_destruct(set1);
    stSet_destruct(set0Prime);
}
static void test_stSet_construct(CuTest* testCase) {
    assert(testCase != NULL);
    testSetup();
    /* Do nothing */
    testTeardown();
}
static void test_stSet_search(CuTest* testCase) {
    testSetup();
    stIntTuple *i = stIntTuple_construct1( 0);
    stIntTuple *j = stIntTuple_construct2(10, 0);
    stIntTuple *k = stIntTuple_construct1( 5);
    //Check search by memory address
    CuAssertTrue(testCase, stSet_search(set0, one) == one);
    CuAssertTrue(testCase, stSet_search(set0, two) == two);
    CuAssertTrue(testCase, stSet_search(set0, three) == three);
    CuAssertTrue(testCase, stSet_search(set0, four) == four);
    CuAssertTrue(testCase, stSet_search(set0, five) == five);
    CuAssertTrue(testCase, stSet_search(set0, six) == six);
    //Check not present
    CuAssertTrue(testCase, stSet_search(set0, i) == NULL);
    CuAssertTrue(testCase, stSet_search(set0, j) == NULL);
    CuAssertTrue(testCase, stSet_search(set0, k) == NULL);
    //Check search by memory address
    CuAssertTrue(testCase, stSet_search(set1, one) == one);
    CuAssertTrue(testCase, stSet_search(set1, two) == two);
    CuAssertTrue(testCase, stSet_search(set1, three) == three);
    CuAssertTrue(testCase, stSet_search(set1, four) == four);
    CuAssertTrue(testCase, stSet_search(set1, five) == five);
    CuAssertTrue(testCase, stSet_search(set1, six) == six);
    //Check not present
    CuAssertTrue(testCase, stSet_search(set1, j) == NULL);
    //Check is searching by memory
    CuAssertTrue(testCase, stSet_search(set1, i) == one);
    CuAssertTrue(testCase, stSet_search(set1, k) == six);
    stIntTuple_destruct(i);
    stIntTuple_destruct(j);
    stIntTuple_destruct(k);
    testTeardown();
}
static void test_stSet_remove(CuTest* testCase) {
    testSetup();
    CuAssertTrue(testCase, stSet_remove(set0, one) == one);
    CuAssertTrue(testCase, stSet_search(set0, one) == NULL);
    CuAssertTrue(testCase, stSet_remove(set1, one) == one);
    CuAssertTrue(testCase, stSet_search(set1, one) == NULL);
    stSet_insert(set1, one);
    CuAssertTrue(testCase, stSet_search(set1, one) == one);
    testTeardown();
}
static void test_stSet_removeAndFreeKey(CuTest* testCase) {
    stSet *set2 = stSet_construct2(free);
    stList *keys = stList_construct();
    int64_t keyNumber = 1000;
    for (int64_t i = 0; i < keyNumber; i++) {
        int64_t *key = st_malloc(sizeof(*key));
        stList_append(keys, key);
        stSet_insert(set2, key);
    }
    for (int64_t i = 0; i < keyNumber; i++) {
        int64_t *key = stList_get(keys, i);
        CuAssertPtrEquals(testCase, key, stSet_removeAndFreeKey(set2, key));
    }
    CuAssertIntEquals(testCase, 0, stSet_size(set2));
    stSet_destruct(set2);
    stList_destruct(keys);
}
static void test_stSet_insert(CuTest* testCase) {
    /*
     * Tests inserting already present keys.
     */
    testSetup();
    CuAssertTrue(testCase, stSet_search(set0, one) == one);
    stSet_insert(set0, one);
    CuAssertTrue(testCase, stSet_search(set0, one) == one);
    stSet_insert(set0, three);
    CuAssertTrue(testCase, stSet_search(set0, three) == three);
    stIntTuple *seven = stIntTuple_construct2(7, 7);
    CuAssertTrue(testCase, stSet_search(set0, seven) == NULL);
    stSet_insert(set0, seven);
    CuAssertTrue(testCase, stSet_search(set0, seven) == seven);
    stIntTuple_destruct(seven);
    testTeardown();
}
static void test_stSet_size(CuTest *testCase) {
    /*
     * Tests the size function of the hash.
     */
    testSetup();
    CuAssertTrue(testCase, stSet_size(set0) == 6);
    CuAssertTrue(testCase, stSet_size(set1) == 6);
    stSet *set2 = stSet_construct();
    CuAssertTrue(testCase, stSet_size(set2) == 0);
    stSet_destruct(set2);
    testTeardown();
}
static void test_stSet_testIterator(CuTest *testCase) {
    testSetup();
    stSetIterator *iterator = stSet_getIterator(set0);
    stSetIterator *iteratorCopy = stSet_copyIterator(iterator);
    int64_t i = 0;
    stSet *seen = stSet_construct();
    for (i = 0; i < 6; i++) {
        void *o = stSet_getNext(iterator);
        CuAssertTrue(testCase, o != NULL);
        CuAssertTrue(testCase, stSet_search(set0, o) != NULL);
        CuAssertTrue(testCase, stSet_search(seen, o) == NULL);
        CuAssertTrue(testCase, stSet_getNext(iteratorCopy) == o);
        stSet_insert(seen, o);
    }
    CuAssertTrue(testCase, stSet_getNext(iterator) == NULL);
    CuAssertTrue(testCase, stSet_getNext(iterator) == NULL);
    CuAssertTrue(testCase, stSet_getNext(iteratorCopy) == NULL);
    stSet_destruct(seen);
    stSet_destructIterator(iterator);
    stSet_destructIterator(iteratorCopy);
    testTeardown();
}
static void test_stSet_testGetKeys(CuTest *testCase) {
    testSetup();
    stList *list = stSet_getKeys(set1);
    CuAssertTrue(testCase, stList_length(list) == 6);
    CuAssertTrue(testCase, stList_contains(list, one));
    CuAssertTrue(testCase, stList_contains(list, two));
    CuAssertTrue(testCase, stList_contains(list, three));
    CuAssertTrue(testCase, stList_contains(list, four));
    CuAssertTrue(testCase, stList_contains(list, five));
    CuAssertTrue(testCase, stList_contains(list, six));
    stList_destruct(list);
    testTeardown();
}
static void test_stSet_getUnion(CuTest* testCase) {
    testSetup();
    // Check union of empty sets is empty
    stSet *set2 = stSet_construct();
    stSet *set3 = stSet_construct();
    stSet *set4 = stSet_getUnion(set2, set3);
    CuAssertTrue(testCase, stSet_size(set4) == 0);
    stSet_destruct(set2);
    stSet_destruct(set3);
    stSet_destruct(set4);
    // Check union of non empty set and empty set is non-empty
    set2 = stSet_construct();
    set3 = stSet_getUnion(set0, set2);
    CuAssertTrue(testCase, stSet_size(set3) == 6);
    stSet_destruct(set2);
    stSet_destruct(set3);
    // Check union of two non-empty overlapping sets is correct
    set2 = stSet_construct();
    set3 = stSet_construct();
    stIntTuple *uniqs[4];
    uniqs[0] = stIntTuple_construct2(9, 0);
    uniqs[1] = stIntTuple_construct2(9, 1);
    uniqs[2] = stIntTuple_construct2(9, 2);
    uniqs[3] = stIntTuple_construct2(9, 3);
    stIntTuple *common[5];
    common[0] = stIntTuple_construct2(5, 0);
    common[1] = stIntTuple_construct2(5, 1);
    common[2] = stIntTuple_construct2(5, 2);
    common[3] = stIntTuple_construct2(5, 3);
    common[4] = stIntTuple_construct2(5, 4);
    for (int i = 0; i < 5; ++i) {
        stSet_insert(set2, common[i]);
        stSet_insert(set3, common[i]);
    }
    stSet_insert(set2, uniqs[0]);
    stSet_insert(set2, uniqs[1]);
    stSet_insert(set3, uniqs[2]);
    stSet_insert(set3, uniqs[3]);
    set4 = stSet_getUnion(set2, set3);
    CuAssertTrue(testCase, stSet_size(set4) == 9);
    for (int i = 0; i < 4; ++i) {
        CuAssertTrue(testCase, stSet_search(set4, uniqs[i]) != NULL);
    }
    for (int i = 0; i < 5; ++i) {
        CuAssertTrue(testCase, stSet_search(set4, common[i]) != NULL);
    }
    stSet_destruct(set2);
    stSet_destruct(set3);
    stSet_destruct(set4);
    stIntTuple_destruct(common[0]);
    stIntTuple_destruct(common[1]);
    stIntTuple_destruct(common[2]);
    stIntTuple_destruct(common[3]);
    stIntTuple_destruct(common[4]);
    stIntTuple_destruct(uniqs[0]);
    stIntTuple_destruct(uniqs[1]);
    stIntTuple_destruct(uniqs[2]);
    stIntTuple_destruct(uniqs[3]);
    // Check we get an exception with sets with different functions.
    stTry {
        stSet_getUnion(set0, set1);
    } stCatch(except) {
        CuAssertTrue(testCase, stExcept_getId(except) == SET_EXCEPTION_ID);
    } stTryEnd
    testTeardown();
}
static void test_stSet_getIntersection(CuTest* testCase) {
    testSetup();
    // Check intersection of empty sets is empty
    stSet *set2 = stSet_construct();
    stSet *set3 = stSet_construct();
    stSet *set4 = stSet_getIntersection(set2, set3);
    CuAssertTrue(testCase, stSet_size(set4) == 0);
    stSet_destruct(set2);
    stSet_destruct(set3);
    stSet_destruct(set4);
    // Check intersection of non empty set and empty set is empy
    set2 = stSet_construct();
    set3 = stSet_getIntersection(set0, set2);
    CuAssertTrue(testCase, stSet_size(set3) == 0);
    stSet_destruct(set2);
    stSet_destruct(set3);
    // Check intersection of two non-empty overlapping sets is correct
    set2 = stSet_construct();
    set3 = stSet_construct();
    stIntTuple *uniqs[4];
    uniqs[0] = stIntTuple_construct2(9, 0);
    uniqs[1] = stIntTuple_construct2(9, 1);
    uniqs[2] = stIntTuple_construct2(9, 2);
    uniqs[3] = stIntTuple_construct2(9, 3);
    stIntTuple *common[5];
    common[0] = stIntTuple_construct2(5, 0);
    common[1] = stIntTuple_construct2(5, 1);
    common[2] = stIntTuple_construct2(5, 2);
    common[3] = stIntTuple_construct2(5, 3);
    common[4] = stIntTuple_construct2(5, 4);
    for (int i = 0; i < 5; ++i) {
        stSet_insert(set2, common[i]);
        stSet_insert(set3, common[i]);
    }
    stSet_insert(set2, uniqs[0]);
    stSet_insert(set2, uniqs[1]);
    stSet_insert(set3, uniqs[2]);
    stSet_insert(set3, uniqs[3]);
    set4 = stSet_getIntersection(set2, set3);
    CuAssertTrue(testCase, stSet_size(set4) == 5);
    stSetIterator *sit = stSet_getIterator(set4);
    stIntTuple *itup;
    while ((itup = stSet_getNext(sit)) != NULL) {
        CuAssertTrue(testCase, stSet_search(set2, itup) != NULL);
        CuAssertTrue(testCase, stSet_search(set3, itup) != NULL);
    }
    for (int i = 0; i < 4; ++i) {
        CuAssertTrue(testCase, stSet_search(set4, uniqs[i]) == NULL);
    }
    stSet_destructIterator(sit);
    stSet_destruct(set2);
    stSet_destruct(set3);
    stSet_destruct(set4);
    stIntTuple_destruct(common[0]);
    stIntTuple_destruct(common[1]);
    stIntTuple_destruct(common[2]);
    stIntTuple_destruct(common[3]);
    stIntTuple_destruct(common[4]);
    stIntTuple_destruct(uniqs[0]);
    stIntTuple_destruct(uniqs[1]);
    stIntTuple_destruct(uniqs[2]);
    stIntTuple_destruct(uniqs[3]);
    // Check we get an exception with sets with different functions.
    stTry {
        stSet_getIntersection(set0, set1);
    } stCatch(except) {
        CuAssertTrue(testCase, stExcept_getId(except) == SET_EXCEPTION_ID);
    } stTryEnd
    testTeardown();
}
static void test_stSet_getDifference(CuTest* testCase) {
    testSetup();
    // Check difference of empty sets is empty
    stSet *set2 = stSet_construct();
    stSet *set3 = stSet_construct();
    stSet *set4 = stSet_getDifference(set2, set3);
    CuAssertTrue(testCase, stSet_size(set4) == 0);
    stSet_destruct(set2);
    stSet_destruct(set3);
    stSet_destruct(set4);
    // Check difference of non empty set and empty set is non-empty
    set2 = stSet_construct();
    set3 = stSet_getDifference(set0, set2);
    CuAssertTrue(testCase, stSet_size(set3) == stSet_size(set0));
    stSet_destruct(set2);
    stSet_destruct(set3);
    set2 = stSet_construct();
    set3 = stSet_getDifference(set2, set0);
    CuAssertTrue(testCase, stSet_size(set3) == 0);
    stSet_destruct(set2);
    stSet_destruct(set3);
    // Check difference of two non-empty overlapping sets is correct
    set2 = stSet_construct();
    set3 = stSet_construct();
    stIntTuple *uniqs[4];
    uniqs[0] = stIntTuple_construct2(9, 0);
    uniqs[1] = stIntTuple_construct2(9, 1);
    uniqs[2] = stIntTuple_construct2(9, 2);
    uniqs[3] = stIntTuple_construct2(9, 3);
    stIntTuple *common[5];
    common[0] = stIntTuple_construct2(5, 0);
    common[1] = stIntTuple_construct2(5, 1);
    common[2] = stIntTuple_construct2(5, 2);
    common[3] = stIntTuple_construct2(5, 3);
    common[4] = stIntTuple_construct2(5, 4);
    for (int i = 0; i < 5; ++i) {
        stSet_insert(set2, common[i]);
        stSet_insert(set3, common[i]);
    }
    stSet_insert(set2, uniqs[0]);
    stSet_insert(set2, uniqs[1]);
    stSet_insert(set3, uniqs[2]);
    stSet_insert(set3, uniqs[3]);
    set4 = stSet_getDifference(set2, set3);
    CuAssertTrue(testCase, stSet_size(set4) == 2);
    for (int i = 0; i < 2; ++i) {
        CuAssertTrue(testCase, stSet_search(set4, uniqs[i]) != NULL);
    }
    for (int i = 2; i < 4; ++i) {
        CuAssertTrue(testCase, stSet_search(set4, uniqs[i]) == NULL);
    }
    for (int i = 2; i < 5; ++i) {
        CuAssertTrue(testCase, stSet_search(set4, common[i]) == NULL);
    }
    stSet_destruct(set2);
    stSet_destruct(set3);
    stSet_destruct(set4);
    stIntTuple_destruct(common[0]);
    stIntTuple_destruct(common[1]);
    stIntTuple_destruct(common[2]);
    stIntTuple_destruct(common[3]);
    stIntTuple_destruct(common[4]);
    stIntTuple_destruct(uniqs[0]);
    stIntTuple_destruct(uniqs[1]);
    stIntTuple_destruct(uniqs[2]);
    stIntTuple_destruct(uniqs[3]);

    // Check we get an exception with sets with different functions.
    stTry {
        stSet_getDifference(set0, set1);
    } stCatch(except) {
        CuAssertTrue(testCase, stExcept_getId(except) == SET_EXCEPTION_ID);
    } stTryEnd
    testTeardown();
}

static void test_stSet_peek(CuTest *testCase) {
    /*
     * Tests the peek function of the set.
     */
    testSetup();
    CuAssertTrue(testCase, stSet_size(set0) == 6);
    CuAssertTrue(testCase, stSet_search(set0, stSet_peek(set0)) != NULL);
    CuAssertTrue(testCase, stSet_size(set1) == 6);
    CuAssertTrue(testCase, stSet_search(set1, stSet_peek(set1)) != NULL);
    stSet *set2 = stSet_construct();
    stTry {
            stSet_peek(set2);
    } stCatch(except) {
        CuAssertTrue(testCase, stExcept_getId(except) == SET_EXCEPTION_ID);
    } stTryEnd
    stSet_destruct(set2);
    testTeardown();
}

static void test_stSet_equals(CuTest *testCase) {
    /*
     * Tests the equals function of the set.
     */
    testSetup();
    CuAssertTrue(testCase, stSet_equals(set0, set0));
    CuAssertTrue(testCase, stSet_equals(set0, set0Prime));
    stSet *set2 = stSet_construct();
    CuAssertTrue(testCase, !stSet_equals(set0, set2));
    CuAssertTrue(testCase, !stSet_equals(set0Prime, set2));
    stSet_destruct(set2);
    testTeardown();
}

static void test_stSet_isSubset(CuTest *testCase) {
    /*
     * Tests the subset function of the set.
     */
    testSetup();
    CuAssertTrue(testCase, stSet_isSubset(set0, set0));
    CuAssertTrue(testCase, stSet_isSubset(set0, set0Prime));
    stSet *set2 = stSet_construct();
    CuAssertTrue(testCase, stSet_isSubset(set0, set2));
    CuAssertTrue(testCase, !stSet_isSubset(set2, set0));
    stSet_destruct(set2);
    testTeardown();
}

CuSuite* sonLib_stSetTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stSet_search);
    SUITE_ADD_TEST(suite, test_stSet_remove);
    SUITE_ADD_TEST(suite, test_stSet_removeAndFreeKey);
    SUITE_ADD_TEST(suite, test_stSet_insert);
    SUITE_ADD_TEST(suite, test_stSet_size);
    SUITE_ADD_TEST(suite, test_stSet_testIterator);
    SUITE_ADD_TEST(suite, test_stSet_construct);
    SUITE_ADD_TEST(suite, test_stSet_testGetKeys);
    SUITE_ADD_TEST(suite, test_stSet_getUnion);
    SUITE_ADD_TEST(suite, test_stSet_getIntersection);
    SUITE_ADD_TEST(suite, test_stSet_getDifference);
    SUITE_ADD_TEST(suite, test_stSet_peek);
    SUITE_ADD_TEST(suite, test_stSet_equals);
    SUITE_ADD_TEST(suite, test_stSet_isSubset);
    return suite;
}

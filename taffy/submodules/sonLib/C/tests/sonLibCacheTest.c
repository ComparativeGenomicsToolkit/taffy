/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * kVDatabaseTest.c
 *
 */

#include "sonLibGlobalsTest.h"
#include "kvDatabaseTestCommon.h"

static stCache *cache = NULL;
static int64_t recordSize;

static void teardown() {
    if(cache != NULL) {
        stCache_destruct(cache);
        cache = NULL;
    }
}

static void setup(size_t maxSize) {
    teardown();
    cache = stCache_construct2(maxSize);
}

static void readAndUpdateRecord(CuTest *testCase) {
    setup(SIZE_MAX);

    CuAssertTrue(testCase, stCache_getRecord(cache, 1, 0, INT64_MAX, &recordSize) == NULL);
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 1, 0, INT64_MAX));

    stCache_setRecord(cache, 1, 5, 6, "hello");
    stCache_setRecord(cache, 1, 12, 6, "world");

    CuAssertTrue(testCase, stCache_getRecord(cache, 1, 0, 5, &recordSize) == NULL); //Check prefix returns false
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 1, 0, 5));

    CuAssertTrue(testCase, stCache_getRecord(cache, 1, 4, 2, &recordSize) == NULL); //Check prefix overlap returns false
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 1, 4, 2));

    char *s = stCache_getRecord(cache, 1, 5, 6, &recordSize);
    CuAssertStrEquals(testCase, "hello", s); //Check we can get the first word
    free(s);
    CuAssertTrue(testCase, recordSize == 6);
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 5, 6));

    CuAssertTrue(testCase, stCache_getRecord(cache, 1, 5, 7, &recordSize) == NULL); //Check suffix overlap returns false
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 1, 5, 7));

    CuAssertTrue(testCase, stCache_getRecord(cache, 1, 11, 1, &recordSize) == NULL); //Check in between returns false
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 1, 11, 1));

    CuAssertTrue(testCase, stCache_getRecord(cache, 1, 11, 2, &recordSize) == NULL); //Check prefix of world returns false
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 1, 11, 2));

    s = stCache_getRecord(cache, 1, 12, 6, &recordSize);
    CuAssertStrEquals(testCase, "world", s); //Get the second word
    free(s);
    CuAssertTrue(testCase, recordSize == 6);
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 12, 6));

    s = stCache_getRecord(cache, 1, 17, 1, &recordSize);
    CuAssertStrEquals(testCase, "", s); //Get part of the second word
    free(s);
    CuAssertTrue(testCase, recordSize == 1);
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 17, 1));

    CuAssertTrue(testCase, stCache_getRecord(cache, 1, 17, 2, &recordSize) == NULL); //Check suffix overlap returns false
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 1, 17, 2));

    CuAssertTrue(testCase, stCache_getRecord(cache, 1, 5, 13, &recordSize) == NULL); //Check suffix overlap returns false
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 1, 5, 13));

    stCache_setRecord(cache, 1, 10, 2, "  ");

    s = stCache_getRecord(cache, 1, 5, INT64_MAX, &recordSize);
    CuAssertStrEquals(testCase, "hello  world", s); //Check we can get the first word
    free(s);
    CuAssertTrue(testCase, recordSize == 13);
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 5, 13));

    stCache_setRecord(cache, 1, 5, 6, "see ya");

    s = stCache_getRecord(cache, 1, 5, 13, &recordSize);
    CuAssertStrEquals(testCase, "see ya world", s); //Check we can get the first word
    free(s);
    CuAssertTrue(testCase, recordSize == 13);
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 5, 13));
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 5, INT64_MAX));

    stCache_clear(cache);

    CuAssertTrue(testCase, !stCache_containsRecord(cache, 1, 12, 6)); //Check its gone

    teardown();
}

static void readAndUpdateRecords(CuTest *testCase) {
    setup(SIZE_MAX);

    stCache_setRecord(cache, 1, 0, 6, "hello");
    stCache_setRecord(cache, 1, 0, 6, "world");
    stCache_setRecord(cache, INT64_MAX-1, 0, 8, "goodbye");
    stCache_setRecord(cache, 3, 0, 6, "cruel");
    stCache_setRecord(cache, INT64_MIN, 0, 6, "earth");

    char *s = stCache_getRecord(cache, 1, 0, INT64_MAX, &recordSize);
    CuAssertStrEquals(testCase, "world", s);
    free(s);
    CuAssertTrue(testCase, recordSize == 6);
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 0, INT64_MAX));

    s = stCache_getRecord(cache, INT64_MAX-1, 0, INT64_MAX, &recordSize);
    CuAssertStrEquals(testCase, "goodbye", s);
    free(s);
    CuAssertTrue(testCase, recordSize == 8);
    CuAssertTrue(testCase, stCache_containsRecord(cache, INT64_MAX-1, 0, INT64_MAX));

    s = stCache_getRecord(cache, 3, 0, INT64_MAX, &recordSize);
    CuAssertStrEquals(testCase, "cruel", s);
    free(s);
    CuAssertTrue(testCase, recordSize == 6);
    CuAssertTrue(testCase, stCache_containsRecord(cache, 3, 0, INT64_MAX));

    s = stCache_getRecord(cache, INT64_MIN, 0, INT64_MAX, &recordSize);
    CuAssertStrEquals(testCase, "earth", s);
    free(s);
    CuAssertTrue(testCase, recordSize == 6);
    CuAssertTrue(testCase, stCache_containsRecord(cache, INT64_MIN, 0, INT64_MAX));

    teardown();
}

static void limitedSizeCache(CuTest *testCase) {
    setup(12);

    stCache_setRecord(cache, 1, 0, 6, "hello");
    stCache_setRecord(cache, 1, 8, 6, "world");
    // Make sure the cache holds 12 bytes
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 0, 6));
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 8, 6));
    CuAssertIntEquals(testCase, 12, stCache_size(cache));

    // This should make the cache overflow, and evict "hello"
    stCache_setRecord(cache, 2, 0, 2, "a");
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 1, 0, 6));
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 8, 6));
    CuAssertTrue(testCase, stCache_containsRecord(cache, 2, 0, 2));
    CuAssertIntEquals(testCase, 8, stCache_size(cache));

    // Add another entry, which doesn't overflow the cache
    stCache_setRecord(cache, 3, 0, 4, "abc");
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 8, 6));
    CuAssertTrue(testCase, stCache_containsRecord(cache, 2, 0, 2));
    CuAssertTrue(testCase, stCache_containsRecord(cache, 3, 0, 4));
    CuAssertIntEquals(testCase, 12, stCache_size(cache));

    // Set the LRU order to 1->2->3
    int64_t size;
    free(stCache_getRecord(cache, 3, 0, 2, &size));
    free(stCache_getRecord(cache, 2, 0, 1, &size));
    free(stCache_getRecord(cache, 1, 8, 6, &size));

    // Overflow by 1 byte. 3 should be first to go
    stCache_setRecord(cache, 4, 0, 1, "o");
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 8, 6));
    CuAssertTrue(testCase, stCache_containsRecord(cache, 2, 0, 2));
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 3, 0, 4));
    CuAssertTrue(testCase, stCache_containsRecord(cache, 4, 0, 1));
    CuAssertIntEquals(testCase, 9, stCache_size(cache));

    // Overflow again. Now 2 should be gone
    stCache_setRecord(cache, 5, 0, 4, "abcd");
    CuAssertTrue(testCase, stCache_containsRecord(cache, 1, 8, 6));
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 2, 0, 2));
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 3, 0, 4));
    CuAssertTrue(testCase, stCache_containsRecord(cache, 4, 0, 1));
    CuAssertTrue(testCase, stCache_containsRecord(cache, 5, 0, 4));
    CuAssertIntEquals(testCase, 11, stCache_size(cache));

    // We still cache strings that are too long for the cache, but
    // they boot everything else out.
    stCache_setRecord(cache, 6, 0, 88, "this is a really long string,"
                      " way too long for the cache. but you should"
                      " still cache me");
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 1, 8, 6));
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 2, 0, 2));
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 3, 0, 4));
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 4, 0, 1));
    CuAssertTrue(testCase, !stCache_containsRecord(cache, 5, 0, 4));
    CuAssertTrue(testCase, stCache_containsRecord(cache, 6, 0, 88));
    CuAssertIntEquals(testCase, 88, stCache_size(cache));

    teardown();
}

// Test the mergeRecords functionality thoroughly (including overlap
// cases). This is also somewhat exercised by readAndUpdateRecord,
// though that only addresses the "inserted record is adjacent to two
// records" case.
static void testMergeRecords(CuTest *testCase) {
    setup(SIZE_MAX);

    // Test merging when there is an existing entry to the right.
    stCache_setRecord(cache, 1, 12, 6, "world");
    stCache_setRecord(cache, 1, 6, 7, "hello w");

    char *s = stCache_getRecord(cache, 1, 6, INT64_MAX, &recordSize);
    CuAssertStrEquals(testCase, "hello world", s);
    free(s);
    CuAssertIntEquals(testCase, 12, stCache_size(cache));

    // Test merging when there is an existing entry to the left.
    stCache_setRecord(cache, 2, 6, 7, "hello w");
    stCache_setRecord(cache, 2, 12, 6, "world");

    s = stCache_getRecord(cache, 2, 6, INT64_MAX, &recordSize);
    CuAssertStrEquals(testCase, "hello world", s);
    free(s);
    CuAssertIntEquals(testCase, 24, stCache_size(cache));

    // Test merging when there is an existing entry on both sides.
    stCache_setRecord(cache, 3, 6, 5, "hello");
    stCache_setRecord(cache, 3, 10, 3, "o w");
    stCache_setRecord(cache, 3, 12, 6, "world");

    s = stCache_getRecord(cache, 3, 6, INT64_MAX, &recordSize);
    CuAssertStrEquals(testCase, "hello world", s);
    free(s);
    CuAssertIntEquals(testCase, 36, stCache_size(cache));

    teardown();
}

CuSuite* stCacheSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, readAndUpdateRecord);
    SUITE_ADD_TEST(suite, readAndUpdateRecords);
    SUITE_ADD_TEST(suite, limitedSizeCache);
    SUITE_ADD_TEST(suite, testMergeRecords);

    return suite;
}



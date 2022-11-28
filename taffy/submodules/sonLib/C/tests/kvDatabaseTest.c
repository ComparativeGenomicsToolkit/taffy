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

static stKVDatabaseConf *conf = NULL;
static stKVDatabase *database = NULL;

static void teardown() {
    if (database != NULL) {
        stKVDatabase_deleteFromDisk(database);
        stKVDatabase_destruct(database);
        database = NULL;
    }
}

static void setup() {
    teardown();
    database = stKVDatabase_construct(conf, true);
    teardown();
    database = stKVDatabase_construct(conf, true);
}

static void constructDestructAndDelete(CuTest *testCase) {
    setup();
    //The setup and teardown functions exercise all the three named functions.
    teardown();
}

static void readWriteAndUpdateIntRecords(CuTest *testCase) {
    setup();
    CuAssertIntEquals(testCase, 0, stKVDatabase_getNumberOfRecords(database));
    //Write some int64 records
    CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, 1));
    stKVDatabase_insertInt64(database, -10000, (int64_t)50);
    stKVDatabase_insertInt64(database, 4297559944418269136, (int64_t)100);
    stKVDatabase_insertInt64(database, 8457478625225021883, (int64_t)150);
    CuAssertTrue(testCase, stKVDatabase_getInt64(database, -10000) == 50);
    CuAssertTrue(testCase, stKVDatabase_getInt64(database, 4297559944418269136) == 100);
    CuAssertTrue(testCase, stKVDatabase_getInt64(database, 8457478625225021883) == 150);

    // test update
    stKVDatabase_insertInt64(database, 4, (int64_t)100);
    stKVDatabase_updateInt64(database, 4, (int64_t)55);
    CuAssertTrue(testCase, stKVDatabase_getInt64(database, 4) == 55);

    
    //Now try removing a record
    stKVDatabase_removeRecord(database, 4);
    CuAssertIntEquals(testCase, 3, stKVDatabase_getNumberOfRecords(database));
    CuAssertPtrEquals(testCase, NULL, stKVDatabase_getRecord(database, 0));
    CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, 0));

    teardown();
}

static void readWriteAndRemoveRecords(CuTest *testCase) {
    setup();
    CuAssertIntEquals(testCase, 0, stKVDatabase_getNumberOfRecords(database));
    //Write some records
    CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, -10000));
    stKVDatabase_insertRecord(database, -10000, "Red", sizeof(char) * 4);
    CuAssertTrue(testCase, stKVDatabase_containsRecord(database, -10000));

    CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, 4297559944418269136));
    stKVDatabase_insertRecord(database, 4297559944418269136, "Green", sizeof(char) * 6);
    CuAssertTrue(testCase, stKVDatabase_containsRecord(database, 4297559944418269136));

    CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, 0));
    stKVDatabase_insertRecord(database, 0, "Black", sizeof(char) * 6);
    CuAssertTrue(testCase, stKVDatabase_containsRecord(database, 0));

    //Now read and check the records exist
    CuAssertIntEquals(testCase, 3, stKVDatabase_getNumberOfRecords(database));
    char *record = stKVDatabase_getRecord(database, -10000);
    CuAssertStrEquals(testCase, "Red", record);
    free(record);
    record = stKVDatabase_getRecord(database, 4297559944418269136);
    CuAssertStrEquals(testCase, "Green", record);
    free(record);
    record = stKVDatabase_getRecord(database, 0);
    CuAssertStrEquals(testCase, "Black", record);
    free(record);
    CuAssertTrue(testCase, stKVDatabase_containsRecord(database, -10000));
    CuAssertTrue(testCase, stKVDatabase_containsRecord(database, 4297559944418269136));
    CuAssertTrue(testCase, stKVDatabase_containsRecord(database, 0));

    //Now check we can retrieve records partially
    record = stKVDatabase_getPartialRecord(database, -10000, 2, 2, sizeof(char) * 4);
    CuAssertStrEquals(testCase, "d", record);
    free(record);
    record = stKVDatabase_getPartialRecord(database, -10000, 1, 3, sizeof(char) * 4);
    CuAssertStrEquals(testCase, "ed", record);
    free(record);
    record = stKVDatabase_getPartialRecord(database, -10000, 0, 4, sizeof(char) * 4);
    CuAssertStrEquals(testCase, "Red", record);
    free(record);
    record = stKVDatabase_getPartialRecord(database, 0, 2, 3, sizeof(char) * 6);
    record[2] = '\0';
    CuAssertStrEquals(testCase, "ac", record);
    free(record);

    //Now try removing the records
    stKVDatabase_removeRecord(database, 0);
    CuAssertIntEquals(testCase, 2, stKVDatabase_getNumberOfRecords(database));
    CuAssertPtrEquals(testCase, NULL, stKVDatabase_getRecord(database, 0));
    CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, 0));
    //Test we get exception if we remove twice.
    stTry {
            stKVDatabase_removeRecord(database, 0);
            CuAssertTrue(testCase, false);
        }
        stCatch(except)
            {
                CuAssertTrue(testCase, stExcept_getId(except) == ST_KV_DATABASE_EXCEPTION_ID);
            }stTryEnd;

    teardown();
}

static void partialRecordRetrieval(CuTest *testCase) {
    setup();

    //Make some number of large records
    stList *records = stList_construct3(0, free);
    stList *recordSizes = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < 300; i++) {
        int64_t size = st_randomInt(0, 80);
        size = size * size * size; //Use cubic size distribution
        char *randomRecord = st_malloc(size * sizeof(char));
        for (int64_t j = 0; j < size; j++) {
            randomRecord[j] = (char) st_randomInt(0, 100);
        }
        stList_append(records, randomRecord);
        stList_append(recordSizes, stIntTuple_construct1( size));
        CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, i));
        stKVDatabase_insertRecord(database, i, randomRecord, size * sizeof(char));
        CuAssertTrue(testCase, stKVDatabase_containsRecord(database, i));
        //st_uglyf("I am creating the record %" PRIi64 " %" PRIi64 "\n", i, size);
    }

    while (st_random() > 0.001) {
        int64_t recordKey = st_randomInt(0, stList_length(records));
        CuAssertTrue(testCase, stKVDatabase_containsRecord(database, recordKey));

        char *record = stList_get(records, recordKey);
        int64_t size = stIntTuple_get(stList_get(recordSizes, recordKey), 0);

        //Get partial record
        int64_t start = size > 0 ? st_randomInt(0, size) : 0;
        int64_t partialSize = size - start > 0 ? st_randomInt(start, size) - start : 0;
        assert(start >= 0);
        assert(partialSize >= 0);
        assert(partialSize + start <= size);
        //st_uglyf("I am getting record %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 "\n", recordKey, start, partialSize, size);
        char *partialRecord = stKVDatabase_getPartialRecord(database, recordKey, start * sizeof(char),
                partialSize * sizeof(char), size * sizeof(char));

        //Check they are equivalent..
        for (int64_t i = 0; i < partialSize; i++) {
            if (record[start + i] != partialRecord[i]) {
                st_uglyf("There was a difference %" PRIi64 " %" PRIi64 " for record %" PRIi64 " %" PRIi64 "\n", record[start + i], partialRecord[i], i,
                        partialSize);
            }
            CuAssertTrue(testCase, record[start + i] == partialRecord[i]);
        }

        free(partialRecord);
        //Check we can not get out of bounds.. (start less than zero)
        stTry {
                free(stKVDatabase_getPartialRecord(database, recordKey, -1, 1, size * sizeof(char)));
            }stCatch(except)
                {
                    CuAssertTrue(testCase, stExcept_getId(except) == ST_KV_DATABASE_EXCEPTION_ID);
                }stTryEnd;

        //Check we can not get out of bounds.. (start greater than index start)
        stTry {
                free(stKVDatabase_getPartialRecord(database, recordKey, size, 1, size * sizeof(char)));
            }stCatch(except)
                {
                    CuAssertTrue(testCase, stExcept_getId(except) == ST_KV_DATABASE_EXCEPTION_ID);
                }stTryEnd;

        //Check we can not get out of bounds.. (total size if greater than record length)
        stTry {
                free(stKVDatabase_getPartialRecord(database, recordKey, 0, size + 1, size * sizeof(char)));
            }stCatch(except)
                {
                    CuAssertTrue(testCase, stExcept_getId(except) == ST_KV_DATABASE_EXCEPTION_ID);
                }stTryEnd;

        //Check we can not get non existent record
        stTry {
                free(stKVDatabase_getPartialRecord(database, 1000000, 0, size, size * sizeof(char)));
            }stCatch(except)
                {
                    CuAssertTrue(testCase, stExcept_getId(except) == ST_KV_DATABASE_EXCEPTION_ID);
                }stTryEnd;
    }

    stList_destruct(records);
    stList_destruct(recordSizes);

    teardown();
}

static void testIncrementRecord(CuTest *testCase) {
    /*
     * Tests incrementing records
     */
    setup();
    // note: record we're incrementing must be the same size (8-byte integer) as
    // the one we're adding to it
    int64_t i = 100;
    int64_t key = 7;
    stKVDatabase_insertInt64(database, key, i);
    CuAssertTrue(testCase,  stKVDatabase_getInt64(database, key) == 100);

    int64_t l = 10;
    int64_t m = 15;
    int64_t n = -25;

    CuAssertIntEquals(testCase, 110, stKVDatabase_incrementInt64(database, key, l));
    CuAssertIntEquals(testCase, 110, stKVDatabase_getInt64(database, key));
    CuAssertIntEquals(testCase, 125, stKVDatabase_incrementInt64(database, key, m));
    CuAssertIntEquals(testCase, 125, stKVDatabase_getInt64(database, key));
    CuAssertIntEquals(testCase, 100, stKVDatabase_incrementInt64(database, key, n));
    CuAssertIntEquals(testCase, 100, stKVDatabase_getInt64(database, key));

    teardown();
}

static void testSetRecord(CuTest *testCase) {
    /*
     * Tests 'setting' records (where you don't know if it's an update or an insert).
     */
    setup();
    int64_t i = 100;
    stKVDatabase_setRecord(database, 1, &i, sizeof(int64_t));
    int64_t *j = stKVDatabase_getRecord(database, 1);
    CuAssertTrue(testCase, i == *j);
    free(j);
    i = 110;
    stKVDatabase_setRecord(database, 1, &i, sizeof(int64_t));
    j = stKVDatabase_getRecord(database, 1);
    CuAssertTrue(testCase, i == *j);
    free(j);
    teardown();
}


static void testBulkSetRecords(CuTest *testCase) {
    /*
     * Tests doing a bulk update of a set of records.
     */
    setup();
    int64_t i = 100, j = 110, k = 120, l = 130;
    stKVDatabase_insertRecord(database, 1, &i, sizeof(int64_t));

    stList *requests = stList_construct3(0, (void(*)(void *)) stKVDatabaseBulkRequest_destruct);
    stList_append(requests, stKVDatabaseBulkRequest_constructInsertRequest(2, &j, sizeof(int64_t)));
    stList_append(requests, stKVDatabaseBulkRequest_constructSetRequest(3, &k, sizeof(int64_t)));
    stList_append(requests, stKVDatabaseBulkRequest_constructUpdateRequest(1, &l, sizeof(int64_t)));

    stKVDatabase_bulkSetRecords(database, requests);

    stList_destruct(requests);

    int64_t *m = stKVDatabase_getRecord(database, 1);
    CuAssertTrue(testCase, m != NULL);
    CuAssertTrue(testCase, l == *m);
    free(m);

    m = stKVDatabase_getRecord(database, 2);
    CuAssertTrue(testCase, m != NULL);
    CuAssertTrue(testCase, j == *m);
    free(m);

    m = stKVDatabase_getRecord(database, 3);
    CuAssertTrue(testCase, m != NULL);
    CuAssertTrue(testCase, k == *m);
    free(m);

    teardown();
}

static void testBulkGetRecords(CuTest* testCase) {
	/*
	 * Tests the new bulk get functions
	 */
	setup();
    int64_t i = 100, j = 110, k = 120, l = 130;
    int64_t bigRecSize = 184500800;
    void *m = st_calloc(1, bigRecSize);

    int64_t ki = 4, kj = 5, kk = 3, kl = 1, km = 2;
    stKVDatabase_insertRecord(database, 1, &i, sizeof(int64_t));

    stList *requests = stList_construct3(0, (void(*)(void *)) stKVDatabaseBulkRequest_destruct);
    stList_append(requests, stKVDatabaseBulkRequest_constructInsertRequest(ki, &i, sizeof(int64_t)));
    stList_append(requests, stKVDatabaseBulkRequest_constructInsertRequest(kj, &j, sizeof(int64_t)));
    stList_append(requests, stKVDatabaseBulkRequest_constructSetRequest(kk, &k, sizeof(int64_t)));
    stList_append(requests, stKVDatabaseBulkRequest_constructUpdateRequest(kl, &l, sizeof(int64_t)));

    stKVDatabase_bulkSetRecords(database, requests);

    stList_destruct(requests);

    stKVDatabase_setRecord(database, km, m, bigRecSize);

    stList* keys = stList_construct2(5);
    stList_set(keys, 0, &ki);
    stList_set(keys, 1, &kj);
    stList_set(keys, 2, &kk);
    stList_set(keys, 3, &kl);
    stList_set(keys, 4, &km);

    stList* results = stKVDatabase_bulkGetRecords(database, keys);
    CuAssertTrue(testCase, stList_length(results) == 5);
    stList_destruct(keys);

    void* record;
    int64_t size;

    stKVDatabaseBulkResult* res0 = (stKVDatabaseBulkResult*)stList_get(results, 0);
    record = stKVDatabaseBulkResult_getRecord(res0, &size);
    CuAssertTrue(testCase, record != NULL);
    CuAssertTrue(testCase, *(int64_t*)record == i && size == sizeof(int64_t));

    stKVDatabaseBulkResult* res1 = (stKVDatabaseBulkResult*)stList_get(results, 1);
	record = stKVDatabaseBulkResult_getRecord(res1, &size);
	CuAssertTrue(testCase, record != NULL);
	CuAssertTrue(testCase, *(int64_t*)record == j && size == sizeof(int64_t));

	stKVDatabaseBulkResult* res2 = (stKVDatabaseBulkResult*)stList_get(results, 2);
	record = stKVDatabaseBulkResult_getRecord(res2, &size);
	CuAssertTrue(testCase, record != NULL);
	CuAssertTrue(testCase, *(int64_t*)record == k && size == sizeof(int64_t));

	stKVDatabaseBulkResult* res3 = (stKVDatabaseBulkResult*)stList_get(results, 3);
	record = stKVDatabaseBulkResult_getRecord(res3, &size);
	CuAssertTrue(testCase, record != NULL);
	CuAssertTrue(testCase, *(int64_t*)record == l && size == sizeof(int64_t));

	stKVDatabaseBulkResult* res4 = (stKVDatabaseBulkResult*)stList_get(results, 4);
	record = stKVDatabaseBulkResult_getRecord(res4, &size);
	CuAssertTrue(testCase, record != NULL);
	CuAssertTrue(testCase, size == bigRecSize);

    stList_destruct(results);

    results = stKVDatabase_bulkGetRecordsRange(database, 1, 6);
    CuAssertTrue(testCase, stList_length(results) == 6);

    res0 = (stKVDatabaseBulkResult*)stList_get(results, 0);
	record = stKVDatabaseBulkResult_getRecord(res0, &size);
	CuAssertTrue(testCase, record != NULL);
	CuAssertTrue(testCase, *(int64_t*)record == l && size == sizeof(int64_t));

	res1 = (stKVDatabaseBulkResult*)stList_get(results, 1);
	record = stKVDatabaseBulkResult_getRecord(res1, &size);
	CuAssertTrue(testCase, record != NULL);
	CuAssertTrue(testCase, size == bigRecSize);

	res2 = (stKVDatabaseBulkResult*)stList_get(results, 2);
	record = stKVDatabaseBulkResult_getRecord(res2, &size);
	CuAssertTrue(testCase, record != NULL);
	CuAssertTrue(testCase, *(int64_t*)record == k && size == sizeof(int64_t));

	res3 = (stKVDatabaseBulkResult*)stList_get(results, 3);
	record = stKVDatabaseBulkResult_getRecord(res3, &size);
	CuAssertTrue(testCase, record != NULL);
	CuAssertTrue(testCase, *(int64_t*)record == i && size == sizeof(int64_t));

	res4 = (stKVDatabaseBulkResult*)stList_get(results, 4);
	record = stKVDatabaseBulkResult_getRecord(res4, &size);
	CuAssertTrue(testCase, record != NULL);
	CuAssertTrue(testCase, *(int64_t*)record == j && size == sizeof(int64_t));

	stKVDatabaseBulkResult* res5 = (stKVDatabaseBulkResult*)stList_get(results, 5);
	record = stKVDatabaseBulkResult_getRecord(res5, &size);
	CuAssertTrue(testCase, record == NULL);

    stList_destruct(results);
    free(m);
    teardown();
}

static void testBulkRemoveRecords(CuTest *testCase) {
    /*
     * Tests doing a bulk update of a set of records.
     */
    setup();
    int64_t i = 100, j = 110, k = 120, l = 130;
    stKVDatabase_insertRecord(database, 1, &i, sizeof(int64_t));
    stKVDatabase_insertRecord(database, 2, &j, sizeof(int64_t));
    stKVDatabase_insertRecord(database, 3, &k, sizeof(int64_t));
    stKVDatabase_insertRecord(database, 4, &l, sizeof(int64_t));
    stKVDatabase_insertRecord(database, 5, &i, 0); //Test null record addition
    CuAssertTrue(testCase, stKVDatabase_containsRecord(database, 1));
    CuAssertTrue(testCase, stKVDatabase_containsRecord(database, 2));
    CuAssertTrue(testCase, stKVDatabase_containsRecord(database, 3));
    CuAssertTrue(testCase, stKVDatabase_containsRecord(database, 4));
    CuAssertTrue(testCase, stKVDatabase_containsRecord(database, 5));
    CuAssertTrue(testCase, stKVDatabase_getNumberOfRecords(database) == 5);

    stList *requests = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);

    // test empty request list
    stKVDatabase_bulkRemoveRecords(database, requests);

    stList_append(requests, stIntTuple_construct1( (int64_t)1));
    stList_append(requests, stIntTuple_construct1( (int64_t)2));
    stList_append(requests, stIntTuple_construct1( (int64_t)3));
    stList_append(requests, stIntTuple_construct1( (int64_t)4));
    stList_append(requests, stIntTuple_construct1( (int64_t)5));
    stKVDatabase_bulkRemoveRecords(database, requests);

    CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, 1));
    CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, 2));
    CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, 3));
    CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, 4));
    CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, 5));
    CuAssertTrue(testCase, stKVDatabase_getNumberOfRecords(database) == 0);

    stList_destruct(requests);

    teardown();
}

/*
 * Retrieves really long records from the database.
 */
static void bigRecordRetrieval(CuTest *testCase) {
    int64_t maxRecordSize = 10000000;
    stKVDatabaseConf_setMaxKTRecordSize(conf, maxRecordSize/5);
    setup();
    for (int64_t i = 0; i < 10; i++) {
        int64_t size = st_randomInt(maxRecordSize/10, maxRecordSize);
        char *randomRecord = st_malloc(size * sizeof(char));
        for (int64_t j = 0; j < size; j++) {
            randomRecord[j] = (char) st_randomInt(0, 100);
        }

        stKVDatabase_insertRecord(database, i, randomRecord, size * sizeof(char));

        if (st_random() > 0.5) {
            // Disconnect and reconnect to see if the split records
            // (in the case of KT) propagate.
            stKVDatabase_destruct(database);
            database = stKVDatabase_construct(conf, false);
        }

        //Check they are equivalent.
        int64_t size2;
        char *randomRecord2 = stKVDatabase_getRecord2(database, i, &size2);
        CuAssertTrue(testCase, size == size2);
        for (int64_t j = 0; j < size; j++) {
            CuAssertTrue(testCase, randomRecord[j] == randomRecord2[j]);
        }
        free(randomRecord);
        free(randomRecord2);
    }
    teardown();
}

/* Check that all tuple records in a set are present and have the expect
 * value.  The expected value in the set is multiplied by valueMult to get
 * the actual expected value */
static void readWriteAndRemoveRecordsLotsCheck(CuTest *testCase, stSortedSet *set, int valueMult) {
    CuAssertIntEquals(testCase, stSortedSet_size(set), stKVDatabase_getNumberOfRecords(database));
    stSortedSetIterator *it = stSortedSet_getIterator(set);
    stIntTuple *tuple;
    while ((tuple = stSortedSet_getNext(it)) != NULL) {
        int64_t *value = (int64_t *) stKVDatabase_getRecord(database, stIntTuple_get(tuple, 0));
        CuAssertTrue(testCase, stKVDatabase_containsRecord(database, stIntTuple_get(tuple, 0)));
        CuAssertIntEquals(testCase, valueMult*stIntTuple_get(tuple, 0), *value);
        free(value);
    }
    stSortedSet_destructIterator(it);
}

static void readWriteAndRemoveRecordsLotsIteration(CuTest *testCase, int numRecords, bool reopenDatabase) {
    //Make a big old list of records..
    stSortedSet *set = stSortedSet_construct3((int(*)(const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);
    while (stSortedSet_size(set) < numRecords) {
        int64_t key = st_randomInt(0, 100 * numRecords);
        stIntTuple *tuple = stIntTuple_construct1( key);
        if (stSortedSet_search(set, tuple) == NULL) {
            CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, key));
            stSortedSet_insert(set, tuple);
            stKVDatabase_insertRecord(database, key, &key, sizeof(int64_t));
            CuAssertTrue(testCase, stKVDatabase_containsRecord(database, key));
        } else {
            CuAssertTrue(testCase, stKVDatabase_containsRecord(database, key));
            stIntTuple_destruct(tuple); // already in db
        }
    }

    readWriteAndRemoveRecordsLotsCheck(testCase, set, 1);

    //Update all records to negate values
    stSortedSetIterator *it = stSortedSet_getIterator(set);
    stIntTuple *tuple;
    while ((tuple = stSortedSet_getNext(it)) != NULL) {
        int64_t *value = (int64_t *) stKVDatabase_getRecord(database, stIntTuple_get(tuple, 0));
        *value *= -1;
        stKVDatabase_updateRecord(database, stIntTuple_get(tuple, 0), value, sizeof(int64_t));
        CuAssertTrue(testCase, stKVDatabase_containsRecord(database, stIntTuple_get(tuple, 0)));
        free(value);
    }
    stSortedSet_destructIterator(it);

    readWriteAndRemoveRecordsLotsCheck(testCase, set, -1);

    //Try optionally committing the transaction and reloading the database..
    if (reopenDatabase) {
        //stKVDatabase_commitTransaction(database);
        stKVDatabase_destruct(database);
        database = stKVDatabase_construct(conf, false);
        //stKVDatabase_startTransaction(database);
    }

    //Now remove each one..
    it = stSortedSet_getIterator(set);
    while ((tuple = stSortedSet_getNext(it)) != NULL) {
        CuAssertTrue(testCase, stKVDatabase_containsRecord(database, stIntTuple_get(tuple, 0)));
        stKVDatabase_removeRecord(database, stIntTuple_get(tuple, 0));
        CuAssertTrue(testCase, !stKVDatabase_containsRecord(database, stIntTuple_get(tuple, 0)));
        //Test we get exception if we remove twice.
        stTry {
                stKVDatabase_removeRecord(database, stIntTuple_get(tuple, 0));
                CuAssertTrue(testCase, 0);
            }
            stCatch(except)
                {
                    CuAssertTrue(testCase, stExcept_getId(except) == ST_KV_DATABASE_EXCEPTION_ID);
                }stTryEnd;
    }
    stSortedSet_destructIterator(it);
    CuAssertIntEquals(testCase, 0, stKVDatabase_getNumberOfRecords(database));

    stSortedSet_destruct(set);
}

/*
 * Tests the retrieval of lots of records.
 */
static void readWriteAndRemoveRecordsLots(CuTest *testCase) {
    setup();
    readWriteAndRemoveRecordsLotsIteration(testCase, 10, false);
    readWriteAndRemoveRecordsLotsIteration(testCase, 56, true);
    readWriteAndRemoveRecordsLotsIteration(testCase, 123, false);
    readWriteAndRemoveRecordsLotsIteration(testCase, 245, true);
    teardown();
}

static void test_stKVDatabaseConf_constructFromString_tokyoCabinet(CuTest *testCase) {
    const char *xmlTestString =
            "<st_kv_database_conf type='tokyo_cabinet'><tokyo_cabinet database_dir='foo'/></st_kv_database_conf>";
    stKVDatabaseConf *conf = stKVDatabaseConf_constructFromString(xmlTestString);
    CuAssertTrue(testCase, stKVDatabaseConf_getType(conf) == stKVDatabaseTypeTokyoCabinet);
    CuAssertStrEquals(testCase, "foo", stKVDatabaseConf_getDir(conf));
    stKVDatabaseConf_destruct(conf);
}

static void test_stKVDatabaseConf_constructFromString_mysql(CuTest *testCase) {
#ifdef HAVE_MYSQL
    const char *xmlTestString =
    "<st_kv_database_conf type='mysql'><mysql host='enormous' port='5' user='foo' password='bar' database_name='mammals' table_name='flowers'/></st_kv_database_conf>";
    stKVDatabaseConf *conf = stKVDatabaseConf_constructFromString(
            xmlTestString);
    CuAssertTrue(testCase, stKVDatabaseConf_getType(conf) == stKVDatabaseTypeMySql);
    CuAssertTrue(testCase, stKVDatabaseConf_getDir(conf) == NULL);
    CuAssertStrEquals(testCase, "enormous", stKVDatabaseConf_getHost(conf));
    CuAssertIntEquals(testCase, 5, stKVDatabaseConf_getPort(conf));
    CuAssertStrEquals(testCase, "foo", stKVDatabaseConf_getUser(conf));
    CuAssertStrEquals(testCase, "bar", stKVDatabaseConf_getPassword(conf));
    CuAssertStrEquals(testCase, "mammals", stKVDatabaseConf_getDatabaseName(conf));
    CuAssertStrEquals(testCase, "flowers", stKVDatabaseConf_getTableName(conf));
#endif
}

static CuSuite* stKVDatabaseTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, readWriteAndRemoveRecords);
    SUITE_ADD_TEST(suite, readWriteAndUpdateIntRecords);
    SUITE_ADD_TEST(suite, readWriteAndRemoveRecordsLots);
    SUITE_ADD_TEST(suite, partialRecordRetrieval);
    SUITE_ADD_TEST(suite, bigRecordRetrieval);
    SUITE_ADD_TEST(suite, testSetRecord);
    SUITE_ADD_TEST(suite, testBulkRemoveRecords);
    SUITE_ADD_TEST(suite, testBulkSetRecords);
    SUITE_ADD_TEST(suite, testBulkGetRecords);
    SUITE_ADD_TEST(suite, constructDestructAndDelete);
    SUITE_ADD_TEST(suite, test_stKVDatabaseConf_constructFromString_tokyoCabinet);
    SUITE_ADD_TEST(suite, test_stKVDatabaseConf_constructFromString_mysql);
    SUITE_ADD_TEST(suite, testIncrementRecord);
    return suite;
}

static int runTests(void) {
    CuString *output = CuStringNew();
    CuSuite* suite = CuSuiteNew();
    CuSuiteAddSuite(suite, stKVDatabaseTestSuite());

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
    CuStringDelete(output);
    int status = suite->failCount > 0;
    CuSuiteDelete(suite);
    return status;
}

int main(int argc, char * const *argv) {
    const char *desc = "kvDatabaseTests [options]\n"
        "\n"
        "Run key/value database tests.\n";
    conf = kvDatabaseTestParseOptions(argc, argv, desc, 0, 0, NULL, NULL);
    return runTests();
}


/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibKVDatabase_KyotoTycoon.cpp
 *
 *  Created on: 5-1-11
 *      Author: epaull
 * 
 * Note: all the KT methods seem to have a C and a CPP version (in terms of the arguments) , 
 * and for this implementation we're using the plain C versions as much as we can.
 *
 * Update October 24, 2011 by Glenn Hickey:
 * Added "secondaryDB" which points to instance of BigRecordFile type database
 * as a fallback for big records.  Now records of a certain size (from conf)
 * don't get added to the kyoto tycoon but go into this new database instead.
 * The whole thing should be transparent to the client, and hopefully prevent
 * the dreaded network errors in the kyoto tycoon API.
 * Note that all operations that can change a record's size must check to make
 * sure that it does not get duplicated across the two db's!
 *
 * Feb 16, 2012:  Multiple database on one server deprecated to implement
 * the binary bulk functions (which require an index).  This functionality
 * could be reintroduced if we keep a name / index mapping externally (or find
 * a way to pry it out of the api)
 */

//Database functions
#ifdef HAVE_KYOTO_TYCOON
#define __STDC_FORMAT_MACROS 1
#include <inttypes.h>
#include <unistd.h>
#include <ktremotedb.h>
#include <kclangc.h>
#include <algorithm>
#include "sonLibGlobalsInternal.h"
#include "sonLibKVDatabasePrivate.h"

using namespace std;
using namespace kyototycoon;

#define SPLIT_RECORDS_KEY "_SPLIT_RECORDS"
#define SPLIT_KEY_SIZE (2 * sizeof(int64_t))

// the default expiration time: negative means indefinite, I believe
int64_t XT = kc::INT64MAX;

typedef struct {
    RemoteDB *rdb;
    stSet *splitRecords;
    int64_t *oldSplitRecords;
    size_t oldSplitRecordsSize;
} KTDB;

static void insertSplitRecordIntoCache(KTDB *db, int64_t key) {
    stIntTuple *splitRecord = stIntTuple_construct1(key);
    stSet_insert(db->splitRecords, splitRecord);
}

static void fillSplitRecordCache(KTDB *db) {
    RemoteDB *rdb = db->rdb;
    size_t splitRecordsSize;
    int64_t *splitRecords = (int64_t *) rdb->get(SPLIT_RECORDS_KEY, sizeof(SPLIT_RECORDS_KEY) - 1,
                                                 &splitRecordsSize, NULL);
    assert(splitRecordsSize % sizeof(int64_t) == 0);
    // Iterate through the split records and add them to our set.
    if (splitRecords != NULL) {
        for (size_t i = 0; i < splitRecordsSize / sizeof(int64_t); i++) {
            insertSplitRecordIntoCache(db, *(splitRecords + i));
        }
    }
    db->oldSplitRecords = splitRecords;
    db->oldSplitRecordsSize = splitRecordsSize;
}

static void removeSplitRecordFromCache(KTDB *db, int64_t key) {
    stIntTuple *tuple = stIntTuple_construct1(key);
    stSet_remove(db->splitRecords, tuple);
    stIntTuple_destruct(tuple);
}

// Atomically sync the split record cache to the database using a
// CAS. If this function returns false, a split record has been added
// / removed by someone else and the cache now reflects those
// changes. Any changes before the sync were discarded; they must be
// made again and this function must be run again to attempt to sync
// them.  If this function returns true, the DB now contains the same
// set of split records as the cache.
static bool syncSplitRecordCache(KTDB *db) {
    // Populate the new split records entry.
    size_t newSplitRecordsSize = stSet_size(db->splitRecords) * sizeof(int64_t);
    int64_t *newSplitRecords = new int64_t[stSet_size(db->splitRecords)];
    stSetIterator *recordIt = stSet_getIterator(db->splitRecords);
    stIntTuple *recordTuple;
    int64_t i = 0;
    while ((recordTuple = (stIntTuple *) stSet_getNext(recordIt)) != NULL) {
        int64_t record = stIntTuple_get(recordTuple, 0);
        newSplitRecords[i] = record;
        i++;
    }
    stSet_destructIterator(recordIt);

    // Attempt a CAS with the existing record.
    bool success = db->rdb->cas(SPLIT_RECORDS_KEY, sizeof(SPLIT_RECORDS_KEY) - 1,
                                (char *) db->oldSplitRecords, db->oldSplitRecordsSize,
                                (char *) newSplitRecords, newSplitRecordsSize);
    delete[] db->oldSplitRecords;
    if (success) {
        // Synced successfully
        db->oldSplitRecords = newSplitRecords;
        db->oldSplitRecordsSize = newSplitRecordsSize;
    } else {
        // Someone changed the value since we last saw it. Our CAS
        // failed, so we erase our cache, fill it afresh with the
        // value now in the DB, and return.
        free(newSplitRecords);
        stSet_destruct(db->splitRecords);
        db->splitRecords = stSet_construct3((uint64_t (*)(const void *)) stIntTuple_hashKey,
                                            (int (*)(const void *, const void *)) stIntTuple_equalsFn,
                                            (void (*)(void *)) stIntTuple_destruct);
        fillSplitRecordCache(db);
    }
    return success;
}

/*
 * construct in the Kyoto Tycoon case means connect to the remote DB
 */
static KTDB *constructDB(stKVDatabaseConf *conf, bool create) {
    KTDB *db = (KTDB *) st_malloc(sizeof(KTDB));
    // we actually do need a local DB dir for Kyoto Tycoon to store the sequences file
    const char *dbDir = stKVDatabaseConf_getDir(conf);
    mkdir(dbDir, S_IRWXU); // just let open of database generate error (FIXME: would be better to make this report errors)


    const char *dbRemote_Host = stKVDatabaseConf_getHost(conf);
    unsigned dbRemote_Port = stKVDatabaseConf_getPort(conf);
    int timeout = stKVDatabaseConf_getTimeout(conf);

    // create the database object
    RemoteDB *rdb = new RemoteDB();

    // tcrdb open sets the host and port for the rdb object
    if (!rdb->open(dbRemote_Host, dbRemote_Port, timeout)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Opening connection to host: %s with error: %s", dbRemote_Host, rdb->error().name());
    }
    db->rdb = rdb;

    // Get the list of split records.
    db->splitRecords = stSet_construct3((uint64_t (*)(const void *)) stIntTuple_hashKey,
                                        (int (*)(const void *, const void *)) stIntTuple_equalsFn,
                                        (void (*)(void *)) stIntTuple_destruct);

    fillSplitRecordCache(db);

    return db;
}

/* closes the remote DB connection and deletes the rdb object, but does not destroy the 
   remote database */
static void destructDB(stKVDatabase *database) {
    KTDB *db = (KTDB *) database->dbImpl;
    if (db != NULL) {
        // close the connection: first try a graceful close, then a forced close
        if (!db->rdb->close(true)) {
            if (!db->rdb->close(false)) {
                stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Closing database error: %s",
                           db->rdb->error().name());
            }
        }
        delete db->rdb;
        stSet_destruct(db->splitRecords);
        delete[] db->oldSplitRecords;
        free(db);
        database->dbImpl = NULL;
    }
}

/* WARNING: removes all records from the remote database */
static void deleteDB(stKVDatabase *database) {
    RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
    if (rdb != NULL) {
        rdb->clear();
    }
    destructDB(database);
    // this removes all records from the remove database object
}


/* check if a record already exists in the kt database*/
static bool recordInTycoon(stKVDatabase *database, int64_t key) {
    RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
    if (rdb->check((char *)&key, (size_t)sizeof(key)) == -1) {
        return false;
    } else {
        return true;
    }
}

static void buildSplitKey(int64_t key_base, int64_t suffix, char *splitKey) {
    memcpy(splitKey, &key_base, sizeof(key_base));
    memcpy(splitKey + sizeof(key_base), &suffix, sizeof(suffix));
}

/* check if a record is split */
static bool recordIsSplit(stKVDatabase *database, int64_t key_base)
{
    KTDB *db = (KTDB *) database->dbImpl;
    stIntTuple *tuple = stIntTuple_construct1(key_base);
    bool found;
    if (stSet_search(db->splitRecords, tuple) == NULL) {
        found = false;
    } else {
        found = true;
    }
    stIntTuple_destruct(tuple);
    return found;
}

static void removeRecordFromTycoonIfPresent(stKVDatabase *database, int64_t key)
{
    if (recordInTycoon(database, key) == true) {
        RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
        if (!rdb->remove((char *)&key, (size_t)sizeof(int64_t))) {
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Removing key/value to database error: %s", rdb->error().name());
        }
    }
}

static int64_t getNumberOfSplitRecords(stKVDatabase *database, int64_t key) {
    int64_t numSplits = 0;
    bool foundEnd = false;
    while (!foundEnd) {
        char splitKey[SPLIT_KEY_SIZE];
        buildSplitKey(key, numSplits, splitKey);
        RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
        size_t sp;
        char *cA;
        if ((cA = rdb->get(splitKey, sizeof(splitKey), &sp, NULL)) == NULL) {
            foundEnd = true;
        } else {
            delete[] cA;
            numSplits++;
        }
    }
    return numSplits;
}

static void removeSplitRecordIfPresent(stKVDatabase *database, int64_t key) {
    if (recordIsSplit(database, key)) {
        if (getenv("ST_SPLIT_RECORD_DBG") != NULL) {
            fprintf(stderr, "removing split record %" PRIi64, key);
        }
        RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
        int64_t numSplits = getNumberOfSplitRecords(database, key);
        char splitKey[SPLIT_KEY_SIZE];
        for (int64_t i = 0; i < numSplits; i++) {
            buildSplitKey(key, i, splitKey);
            if (!rdb->remove(splitKey, sizeof(splitKey))) {
                stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Removing key/value to database error: %s", rdb->error().name());
            }
        }
        // Remove the split record marker from the DB. This uses a
        // CAS, so it may require several attempts.
        do {
            removeSplitRecordFromCache((KTDB *) database->dbImpl, key);
        } while (!syncSplitRecordCache((KTDB *) database->dbImpl));
    }
}

static bool containsRecord(stKVDatabase *database, int64_t key) {
    bool found = recordInTycoon(database, key);
    if (!found) {
        found = recordIsSplit(database, key);
    }
    return found;
}

static void insertRecord(stKVDatabase *database, int64_t key, const void *value, int64_t sizeOfRecord) {
    stKVDatabaseConf* conf = stKVDatabase_getConf(database);
    int64_t maxRecordSize = stKVDatabaseConf_getMaxKTRecordSize(conf);
    RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
    if (sizeOfRecord > maxRecordSize)
    {
        removeSplitRecordIfPresent(database, key);
        int64_t numSplitRecords = sizeOfRecord / maxRecordSize + 1;
        if (getenv("ST_SPLIT_RECORD_DBG") != NULL) {
            fprintf(stderr, "attempting to create %" PRIi64 " records for key %" PRIi64 "\n", numSplitRecords, key);
        }
        char splitKey[SPLIT_KEY_SIZE];
        for (int64_t i = 0; i < numSplitRecords; i++) {
            buildSplitKey(key, i, splitKey);
            if (!rdb->add(splitKey, sizeof(splitKey),
                          ((const char *)value) + i * maxRecordSize,
                          (i == numSplitRecords - 1)
                            ? sizeOfRecord - i * maxRecordSize
                            : maxRecordSize)) {
                stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Inserting key/value to database error: %s", rdb->error().name());
            }
        }
        // Add the split record marker to the DB. This uses a
        // CAS, so it may require several attempts.
        do {
            insertSplitRecordIntoCache((KTDB *) database->dbImpl, key);
        } while (!syncSplitRecordCache((KTDB *) database->dbImpl));
        assert(getNumberOfSplitRecords(database, key) == numSplitRecords);
        assert(recordIsSplit(database, key));
    }
    else
    {
        removeSplitRecordIfPresent(database, key);

        // add method: If the key already exists the record will not be modified and it'll return false
        if (!rdb->add((char *)&key, sizeof(key), (const char *)value, sizeOfRecord)) {
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Inserting key/value to database error: %s", rdb->error().name());
        }
    }
}

static void insertInt64(stKVDatabase *database, int64_t key, int64_t value) {
    RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;

    // Normalize a 64-bit number in the native order into the network byte order.
    // little endian (our x86 linux machine) to big Endian....
    int64_t KCSafeIV = kyotocabinet::hton64(value);

    if (!rdb->add((char *)&key, sizeof(int64_t), (const char *)&KCSafeIV, sizeof(int64_t))) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Inserting int64 key/value to database error: %s", rdb->error().name());
    }
}

static void updateInt64(stKVDatabase *database, int64_t key, int64_t value) {
    RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;

    // Normalize a 64-bit number in the native order into the network byte order.
    // little endian (our x86 linux machine) to big Endian....
    uint64_t KCSafeIV = kyotocabinet::hton64(value);

    if (!rdb->replace((char *)&key, sizeof(int64_t), (const char *)&KCSafeIV, sizeof(int64_t))) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Updating int64 key/value to database error: %s", rdb->error().name());
    }
}

static void updateRecord(stKVDatabase *database, int64_t key, const void *value, int64_t sizeOfRecord) {
    stKVDatabaseConf* conf = stKVDatabase_getConf(database);
    int64_t maxRecordSize = stKVDatabaseConf_getMaxKTRecordSize(conf);
    if (sizeOfRecord > maxRecordSize)
    {
        if (getenv("ST_SPLIT_RECORD_DBG") != NULL) {
            fprintf(stderr, "updateRecord %" PRIi64 "\n", key);
        }
        removeSplitRecordIfPresent(database, key);
        insertRecord(database, key, value, sizeOfRecord);
    }
    else
    {
        removeSplitRecordIfPresent(database, key);
        RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
        // replace method: If the key doesn't already exist it won't be created, and we'll get an error
        if (!rdb->replace((char *)&key, (size_t)sizeof(int64_t), (const char *)value, sizeOfRecord)) {
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Updating key/value to database error: %s", rdb->error().name());
        }
    }
}

static void setRecord(stKVDatabase *database, int64_t key, const void *value, int64_t sizeOfRecord) {
    stKVDatabaseConf* conf = stKVDatabase_getConf(database);
    int64_t maxRecordSize = stKVDatabaseConf_getMaxKTRecordSize(conf);
    if (sizeOfRecord > maxRecordSize)
    {
        removeSplitRecordIfPresent(database, key);
        insertRecord(database, key, value, sizeOfRecord);
    }
    else
    {
        removeSplitRecordIfPresent(database, key);
        RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
        if (!rdb->set((char *)&key, (size_t)sizeof(int64_t), (const char *)value, sizeOfRecord)) {
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "kyoto tycoon setting key/value failed: %s", rdb->error().name());
        }
    }
}

/* increment a record by the specified numerical value: atomic operation */
/* return the new record value */
static int64_t incrementInt64(stKVDatabase *database, int64_t key, int64_t incrementAmount) {
    RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
    int64_t returnValue = kyotocabinet::INT64MIN;

    size_t sizeOfKey = sizeof(int64_t);

    if ( (returnValue = rdb->increment((char *)&key, sizeOfKey, incrementAmount, kyotocabinet::INT64MIN, XT)) == kyotocabinet::INT64MIN ) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "kyoto tycoon incremement record failed: %s", rdb->error().name());
    }

    return returnValue;
}

// sets a bulk list of records atomically 
static void bulkSetRecords(stKVDatabase *database, stList *records) {
    stKVDatabaseConf* conf = stKVDatabase_getConf(database);
    int64_t maxRecordSize = stKVDatabaseConf_getMaxKTRecordSize(conf);
    int64_t maxBulkSetSize = stKVDatabaseConf_getMaxKTBulkSetSize(conf);
    int64_t maxBulkSetNumRecords = stKVDatabaseConf_getMaxKTBulkSetNumRecords(conf);
    RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
    vector<RemoteDB::BulkRecord> recs;
    recs.reserve(stList_length(records));
    RemoteDB::BulkRecord templateRec;
    templateRec.dbidx = 0;
    templateRec.xt = XT;
    int64_t runningSize = 0;

    // copy the records from our C data structure to the CPP map needed for the Tycoon API
    for(int32_t i=0; i<stList_length(records); i++) {
        stKVDatabaseBulkRequest *request = (stKVDatabaseBulkRequest *)stList_get(records, i);

        // current batch can't get any bigger so we write and clear it
        if ((runningSize + request->size > maxBulkSetSize ||
             (int64_t)recs.size() >= maxBulkSetNumRecords) && recs.empty() == false) {
            int64_t retVal = rdb->set_bulk_binary(recs);
            if (retVal < 1) {
                assert(rdb->error().name() != NULL);
                fprintf(stderr, "Throwing an exception with the string %s\n", rdb->error().name());
                stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "kyoto tycoon set bulk record failed: %s", rdb->error().name());
            }
            recs.clear();
            runningSize = 0;
        }
        // record too big for kt, we put in the secondary
        if (request->size > maxRecordSize) {
            removeRecordFromTycoonIfPresent(database, request->key);
            removeSplitRecordIfPresent(database, request->key);
            setRecord(database, request->key, request->value, request->size);
        }
        else
        {
            removeSplitRecordIfPresent(database, request->key);
            templateRec.key = string((const char *)&(request->key), sizeof(int64_t));
            templateRec.value = string((const char *)request->value, request->size);
            recs.push_back(templateRec);
            runningSize += request->size;
        }
    }

    // test for empty list   
    if (recs.empty() == false) {
        // set values, atomic = true
        int64_t retVal = rdb->set_bulk_binary(recs);
        if (retVal < 1) {
            assert(rdb->error().name() != NULL);
            fprintf(stderr, "Throwing an exception with the string %s\n", rdb->error().name());
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "kyoto tycoon set bulk record failed: %s", rdb->error().name());
        }
    }
}

// remove a bulk list atomically 
static void bulkRemoveRecords(stKVDatabase *database, stList *records) {
    stKVDatabaseConf* conf = stKVDatabase_getConf(database);
    RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
    vector<string> keys;
    int64_t maxBulkSetNumRecords = stKVDatabaseConf_getMaxKTBulkSetNumRecords(conf);
    int64_t num_records = stList_length(records);
    
    for(int32_t i=0; i<num_records; i++) {
        int64_t key = stIntTuple_get((stIntTuple *)stList_get(records, i), 0);
        if (recordIsSplit(database, key) == true) {
            removeSplitRecordIfPresent(database, key);
        }
        else {
            keys.push_back(string((const char *)&key, sizeof(int64_t)));
        }
        if (!keys.empty() && (i == num_records - 1 || keys.size() >= maxBulkSetNumRecords)) {
            if (rdb->remove_bulk(keys, true) < 1) {
                stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "kyoto tycoon bulk remove record failed: %s", rdb->error().name());
            }
            keys.clear();
        }
    }
}

static int64_t numberOfRecords(stKVDatabase *database) {
    RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
    int64_t count = rdb->count();
    return count;
}

static void *getRecord2(stKVDatabase *database, int64_t key, int64_t *recordSize) {
    char* record = NULL;
    stKVDatabaseConf* conf = stKVDatabase_getConf(database);
    int64_t maxRecordSize = stKVDatabaseConf_getMaxKTRecordSize(conf);
    RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
    if (recordIsSplit(database, key) == true)
    {
        int64_t numSplits = getNumberOfSplitRecords(database, key);
        if (getenv("ST_SPLIT_RECORD_DBG") != NULL) {
            fprintf(stderr, "attempting to read %" PRIi64 " split records from record %" PRIi64 "\n",
                    numSplits, key);
        }
        record = (char *) st_malloc(maxRecordSize * numSplits);
        *recordSize = 0;
        char splitKey[SPLIT_KEY_SIZE];
        for (int64_t i = 0; i < numSplits; i++) {
            buildSplitKey(key, i, splitKey);
            size_t subSize;
            char *subRecord = rdb->get(splitKey, sizeof(splitKey), &subSize, NULL);
            *recordSize += subSize;
            memcpy(record + i * maxRecordSize, subRecord, subSize);
            delete[] subRecord;
        }
    }
    else
    {
        //Return value must be freed.
        size_t i;
        char* newRecord = rdb->get((char *)&key, (size_t)sizeof(int64_t), &i, NULL);
        if (newRecord == NULL) {
            // No key.
            record = NULL;
        } else {
            *recordSize = (int64_t)i;
            record = (char*)memcpy(st_malloc(*recordSize), newRecord, *recordSize);
            delete[] newRecord;
        }
    }

    return record;
}

/* get a single non-string record */
static void *getRecord(stKVDatabase *database, int64_t key) {
    int64_t i;
    return getRecord2(database, key, &i);
}

/* get a single non-string record */
static int64_t getInt64(stKVDatabase *database, int64_t key) {
    RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;

    size_t sp;
    char *newRecord = rdb->get((char *)&key, sizeof(int64_t), &sp, NULL);
    char* record = (char*)memcpy(st_malloc( sizeof(int64_t)), newRecord,  sizeof(int64_t));
    delete[] newRecord;

    // convert from KC native big-endian back to little-endian Intel...
    int64_t ret = kyotocabinet::ntoh64(*((int64_t*)record));
    free(record);
    return ret;
}

/* get part of a string record */
static void *getPartialRecord(stKVDatabase *database, int64_t key, int64_t zeroBasedByteOffset, int64_t sizeInBytes, int64_t recordSize) {
    // NB: does not treat split records specially right now, so this
    // could potentitally be more efficient on split records. But I
    // don't think anything uses this right now, so there isn't much
    // point making that optimization.
    int64_t recordSize2;
    char *record = (char *)getRecord2(database, key, &recordSize2);
    if(recordSize2 != recordSize) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "The given record size is incorrect: %lld, should be %lld", (long long)recordSize, recordSize2);
    }
    if(record == NULL) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "The record does not exist: %lld for partial retrieval", (long long)key);
    }
    if(zeroBasedByteOffset < 0 || sizeInBytes < 0 || zeroBasedByteOffset + sizeInBytes > recordSize) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Partial record retrieval to out of bounds memory, record size: %lld, requested start: %lld, requested size: %lld", (long long)recordSize, (long long)zeroBasedByteOffset, (long long)sizeInBytes);
    }
    void *partialRecord = memcpy(st_malloc(sizeInBytes), record + zeroBasedByteOffset, sizeInBytes);
    free(record);
    return partialRecord;
}

/* do a bulk get based on a list of keys.  */
static stList *bulkGetRecords(stKVDatabase *database, stList* keys)
{
    stKVDatabaseConf* conf = stKVDatabase_getConf(database);
    int64_t maxBulkSetNumRecords = stKVDatabaseConf_getMaxKTBulkSetNumRecords(conf);
    int32_t n = stList_length(keys);
    RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
    RemoteDB::BulkRecord templateRec;
    templateRec.dbidx = 0;
    templateRec.xt = XT;
    vector<RemoteDB::BulkRecord> recs;
    recs.reserve(std::min(maxBulkSetNumRecords, (int64_t)n));
    stList* results = stList_construct3(n, (void(*)(void *))stKVDatabaseBulkResult_destruct);
    int32_t prevI = 0;
    for (int32_t i = 0; i < n; ++i)
    {
        int64_t key = *(int64_t*)stList_get(keys, i);
        if (recordIsSplit(database, key) == true)
        {
            int64_t recordSize;
            char *record = (char *)getRecord2(database, key, &recordSize);
            stKVDatabaseBulkResult* result = stKVDatabaseBulkResult_construct(record, recordSize);
            stList_set(results, i, result);
        }
        else
        {
            templateRec.key = string((char*)stList_get(keys, i), (size_t)sizeof(int64_t));
            recs.push_back(templateRec);
        }

        if (recs.size() >= (uint64_t)maxBulkSetNumRecords || (i == n - 1 && !recs.empty()))
        {
            int64_t retVal = rdb->get_bulk_binary(&recs);
            if (retVal < 0)
            {
                assert(rdb->error().name() != NULL);
                fprintf(stderr, "Throwing a KT exception with the string %s\n", rdb->error().name());
                stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "kyoto tycoon get bulk record failed: %s", rdb->error().name());
            }
            int32_t recIdx = 0;
            for (int32_t j = prevI; j <= i; ++j)
            {
                if (stList_get(results, j) == NULL)
                {
                    RemoteDB::BulkRecord& curRecord = recs.at(recIdx++);
                    int64_t recordSize = curRecord.value.length() * sizeof(char);
                    void* record = st_malloc(recordSize);
                    memcpy(record, curRecord.value.data(), recordSize);
                    stKVDatabaseBulkResult* result = stKVDatabaseBulkResult_construct(record, recordSize);
                    stList_set(results, j, result);
                }
            }
            recs.clear();
            prevI = i + 1;
        }
    }
    return results;
}

static stList *bulkGetRecordsRange(stKVDatabase *database, int64_t firstKey, int64_t numRecords) {
    vector<string> keysVec;
    keysVec.reserve(numRecords);
    stList* results = stList_construct3(numRecords, (void(*)(void *))stKVDatabaseBulkResult_destruct);
    for (int64_t i = 0; i < numRecords; ++i) {
        int64_t key = firstKey + i;
        if (recordIsSplit(database, key) == true)
        {
            int64_t recordSize;
            char *record = (char *)getRecord2(database, key, &recordSize);
            stKVDatabaseBulkResult* result = stKVDatabaseBulkResult_construct(record, recordSize);
            stList_set(results, (int32_t)i, result);
        }
        else
        {
            keysVec.push_back(string((char*)&key, (size_t)sizeof(int64_t)));
        }
    }
    if (keysVec.empty() == false)
    {
        RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
        map<string, string> recs;
        int64_t retVal = rdb->get_bulk(keysVec, &recs);
        if (retVal < 0)
        {
            assert(rdb->error().name() != NULL);
            fprintf(stderr, "Throwing a KT exception with the string %s\n", rdb->error().name());
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "kyoto tycoon get bulk record failed: %s", rdb->error().name());
        }
        for (int64_t i = 0; i < numRecords; ++i)
        {
            if (stList_get(results, (int32_t)i) == NULL)
            {
                int64_t key = firstKey + i;
                string keyString = string((char*)&key, (size_t)sizeof(int64_t));
                map<string,string>::iterator mapIt = recs.find(keyString);
                void* record = NULL;
                int64_t recordSize = 0;
                if (mapIt != recs.end())
                {
                    recordSize = mapIt->second.length() * sizeof(char);
                    record = st_malloc(recordSize);
                    memcpy(record, mapIt->second.data(), recordSize);
                }
                stKVDatabaseBulkResult* result = stKVDatabaseBulkResult_construct(record, recordSize);
                stList_set(results, (int32_t)i, result);
            }
        }
    }
    return results;
}

static void removeRecord(stKVDatabase *database, int64_t key) {
    if (recordIsSplit(database, key) == true) {
        removeSplitRecordIfPresent(database->secondaryDB, key);
    }
    else
    {
        RemoteDB *rdb = ((KTDB *) database->dbImpl)->rdb;
        if (!rdb->remove((char *)&key, (size_t)sizeof(int64_t))) {
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Removing key/value to database error: %s", rdb->error().name());
        }
    }
}

void stKVDatabase_initialise_kyotoTycoon(stKVDatabase *database, stKVDatabaseConf *conf, bool create) {
    database->dbImpl = constructDB(stKVDatabase_getConf(database), create);
    database->destruct = destructDB;
    database->deleteDatabase = deleteDB;
    database->containsRecord = containsRecord;
    database->insertRecord = insertRecord;
    database->insertInt64 = insertInt64;
    database->updateRecord = updateRecord;
    database->updateInt64 = updateInt64;
    database->setRecord = setRecord;
    database->incrementInt64 = incrementInt64;
    database->bulkSetRecords = bulkSetRecords;
    database->bulkRemoveRecords = bulkRemoveRecords;
    database->numberOfRecords = numberOfRecords;
    database->getRecord = getRecord;
    database->getInt64 = getInt64;
    database->getRecord2 = getRecord2;
    database->getPartialRecord = getPartialRecord;
    database->bulkGetRecords = bulkGetRecords;
    database->bulkGetRecordsRange = bulkGetRecordsRange;
    database->removeRecord = removeRecord;
}

#endif

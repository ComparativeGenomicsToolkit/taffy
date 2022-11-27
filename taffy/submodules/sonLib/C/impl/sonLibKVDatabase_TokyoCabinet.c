/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibKVDatabase_TokyoCabinet.c
 *
 *  Created on: 18-Aug-2010
 *      Author: benedictpaten
 */

//Database functions

#include "sonLibGlobalsInternal.h"
#include "sonLibKVDatabasePrivate.h"

#ifdef HAVE_TOKYO_CABINET
#include <tcutil.h>
#include <tcbdb.h>

static int keyCmp(const char *vA1, int size1, const char *vA2,
        int size2, void *a) {
    assert(size1 == sizeof(int64_t));
    assert(size2 == sizeof(int64_t));
    assert(a == NULL);
    int64_t i = *(int64_t *) vA1;
    int64_t j = *(int64_t *) vA2;
    return i - j > 0 ? 1 : (i < j ? -1 : 0);
}

static TCBDB *constructDB(stKVDatabaseConf *conf, bool create) {
    const char *dbDir = stKVDatabaseConf_getDir(conf);
    mkdir(dbDir, S_IRWXU); // just let open of database generate error (FIXME: would be better to make this report errors)
    char *databaseName = stString_print("%s/%s", dbDir, "data");
    TCBDB *dbImpl = tcbdbnew();
    tcbdbsetcmpfunc(dbImpl, keyCmp, NULL);
    unsigned opts = BDBOWRITER | (create ? BDBOCREAT|BDBOTRUNC : 0);
    if (!tcbdbopen(dbImpl, databaseName, opts)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Opening database: %s with error: %s", databaseName, tcbdberrmsg(tcbdbecode(dbImpl)));
    }
    free(databaseName);
    return dbImpl;
}

static void destructDB(stKVDatabase *database) {
    TCBDB *dbImpl = database->dbImpl;
    if (dbImpl != NULL) {
        if (!tcbdbclose(dbImpl)) {
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Closing database error: %s", tcbdberrmsg(tcbdbecode(dbImpl)));
        }
        tcbdbdel(dbImpl);
        database->dbImpl = NULL;
    }
}

static void deleteDB(stKVDatabase *database) {
    destructDB(database);
    const char *dbDir = stKVDatabaseConf_getDir(stKVDatabase_getConf(database));
    int64_t i = st_system("rm -rf %s", dbDir);
    if (i != 0) {
        st_errAbort("Tried to delete the temporary cactus disk: %s with exit code %d", dbDir, i);
    }
}

/* check if a record already exists */
static bool recordExists(TCBDB *dbImpl, int64_t key) {
    return tcbdbvnum(dbImpl, &key, sizeof(int64_t)) > 0;
}

static bool containsRecord(stKVDatabase *database, int64_t key) {
    return recordExists(database->dbImpl, key);
}

static void insertRecord(stKVDatabase *database, int64_t key, const void *value, int64_t sizeOfRecord) {
    TCBDB *dbImpl = database->dbImpl;
    if (recordExists(dbImpl, key)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Attempt to insert a key in the database that already exists: %lld", (long long)key);
    }
    if (!tcbdbput(dbImpl, &key, sizeof(int64_t), value, sizeOfRecord)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Inserting key/value to database error: %s", tcbdberrmsg(tcbdbecode(dbImpl)));
    }
}

static void insertInt64(stKVDatabase *database, int64_t key, int64_t value) {
    TCBDB *dbImpl = database->dbImpl;
    
    if (recordExists(dbImpl, key)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Attempt to insert a key in the database that already exists: %lld", (long long)key);
    }
    if (!tcbdbput(dbImpl, &key, sizeof(int64_t), (const void *)&value, sizeof(int64_t))) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Inserting key/value to database error: %s", tcbdberrmsg(tcbdbecode(dbImpl)));
    }
}

static void updateRecord(stKVDatabase *database, int64_t key, const void *value, int64_t sizeOfRecord) {
    TCBDB *dbImpl = database->dbImpl;
    if (!recordExists(dbImpl, key)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Attempt to update a key in the database that doesn't exists: %lld", (long long)key);
    }
    if (!tcbdbput(dbImpl, &key, sizeof(int64_t), value, sizeOfRecord)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Updating key/value to database error: %s", tcbdberrmsg(tcbdbecode(dbImpl)));
    }
}

static void updateInt64(stKVDatabase *database, int64_t key, int64_t value) {
    TCBDB *dbImpl = database->dbImpl;
    if (!recordExists(dbImpl, key)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Attempt to update a key in the database that doesn't exists: %lld", (long long)key);
    }
    if (!tcbdbput(dbImpl, &key, sizeof(int64_t), (const void *)&value, sizeof(int64_t))) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Updating key/value to database error: %s", tcbdberrmsg(tcbdbecode(dbImpl)));
    }
}

static void setRecord(stKVDatabase *database, int64_t key,
        const void *value, int64_t sizeOfRecord) {
    TCBDB *dbImpl = database->dbImpl;
    if (!tcbdbput(dbImpl, &key, sizeof(int64_t), value, sizeOfRecord)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Set key/value to database error: %s", tcbdberrmsg(tcbdbecode(dbImpl)));
    }
}

static void startTransaction(stKVDatabase *database) {
    TCBDB *dbImpl = database->dbImpl;
    if (!tcbdbtranbegin(dbImpl)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Tried to start a transaction but got error: %s", tcbdberrmsg(tcbdbecode(dbImpl)));
    }
}

static void commitTransaction(stKVDatabase *database) {
    TCBDB *dbImpl = database->dbImpl;
    //Commit the transaction..
    if (!tcbdbtrancommit(dbImpl)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Tried to commit a transaction but got error: %s", tcbdberrmsg(tcbdbecode(dbImpl)));
    }
}

static void abortTransaction(stKVDatabase *database) {
    TCBDB *dbImpl = database->dbImpl;
    if (!tcbdbtranabort(dbImpl)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Tried to abort a transaction but got error: %s", tcbdberrmsg(tcbdbecode(dbImpl)));
    }
}

static void bulkSetRecords(stKVDatabase *database, stList *records) {
    startTransaction(database);
    stTry {
        for(int64_t i=0; i<stList_length(records); i++) {
            stKVDatabaseBulkRequest *request = stList_get(records, i);
            switch(request->type) {
                case UPDATE:
                updateRecord(database, request->key, request->value, request->size);
                break;
                case INSERT:
                insertRecord(database, request->key, request->value, request->size);
                break;
                case SET:
                setRecord(database, request->key, request->value, request->size);
                break;
            }
        }
        commitTransaction(database);
    }stCatch(ex) {
        abortTransaction(database);
        stThrowNewCause(
                ex,
                ST_KV_DATABASE_EXCEPTION_ID,
                "tokyo cabinet bulk set records failed");
    }stTryEnd;
}

static int64_t numberOfRecords(stKVDatabase *database) {
    TCBDB *dbImpl = database->dbImpl;
    return tcbdbrnum(dbImpl);
}

static void *getRecord2(stKVDatabase *database, int64_t key, int64_t *recordSize) {
    TCBDB *dbImpl = database->dbImpl;
    //Return value must be freed.
    int32_t i;
    void *record = tcbdbget(dbImpl, &key, sizeof(int64_t), &i);
    *recordSize = i;
    return record;
}

static int64_t getInt64(stKVDatabase *database, int64_t key) {
    TCBDB *dbImpl = database->dbImpl;
    //Return value must be freed.
    int32_t i;
    void *record = tcbdbget(dbImpl, &key, sizeof(int64_t), &i);

    if (record == NULL) {
        return -1;
    } else {
        int64_t ret = *((int64_t*)record);
        free(record);
        return ret;
    }
}

static void *getRecord(stKVDatabase *database, int64_t key) {
    int64_t i;
    return getRecord2(database, key, &i);
}

static int64_t incrementInt64(stKVDatabase *database, int64_t key, int64_t incrementAmount) {
    startTransaction(database);
    int64_t returnValue = INT64_MIN;
    stTry {
        int64_t recordSize;
        int64_t *record = getRecord2(database, key, &recordSize);
        assert(recordSize >= sizeof(int64_t));
        record[0] += incrementAmount;
        returnValue = record[0];
        updateRecord(database, key, record, recordSize);
        free(record);
        commitTransaction(database);
    }stCatch(ex) {
        abortTransaction(database);
        stThrowNewCause(
                ex,
                ST_KV_DATABASE_EXCEPTION_ID,
                "tokyo cabinet increment record failed");
    }stTryEnd;
    return returnValue;
}

static void *getPartialRecord(stKVDatabase *database, int64_t key, int64_t zeroBasedByteOffset, int64_t sizeInBytes, int64_t recordSize) {
    int64_t recordSize2;
    char *record = getRecord2(database, key, &recordSize2);
    if(recordSize2 != recordSize) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "The given record size is incorrect: %lld, should be %lld", (long long)recordSize, recordSize2);
    }
    if(record == NULL) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "The record does not exist: %lld for partial retrieval", (long long)key);
    }
    if(zeroBasedByteOffset < 0 || sizeInBytes < 0 || zeroBasedByteOffset + sizeInBytes > recordSize) {
        free(record);
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Partial record retrieval to out of bounds memory, record size: %lld, requested start: %lld, requested size: %lld", (long long)recordSize, (long long)zeroBasedByteOffset, (long long)sizeInBytes);
    }
    void *partialRecord = memcpy(st_malloc(sizeInBytes), record + zeroBasedByteOffset, sizeInBytes);
    free(record);
    return partialRecord;
}

static stList *bulkGetRecords(stKVDatabase *database, stList* keys) {
	int64_t n = stList_length(keys);
	stList* results = stList_construct3(n, (void(*)(void *))stKVDatabaseBulkResult_destruct);
	//startTransaction(database);
	stTry {
		for (int64_t i = 0; i < n; ++i)
		{
			int64_t key = *((int64_t*)stList_get(keys, i));
			int64_t recordSize;
			void* record = getRecord2(database, key, &recordSize);
			stKVDatabaseBulkResult* result = stKVDatabaseBulkResult_construct(record, recordSize);
			stList_set(results, i, result);
		}
		//commitTransaction(database);
	}stCatch(ex) {
		//abortTransaction(database);
		stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID, "tokyo cabinet bulk get records failed");
	}stTryEnd;

	return results;
}

static stList *bulkGetRecordsRange(stKVDatabase *database, int64_t firstKey, int64_t numRecords) {
	stList* results = stList_construct3(numRecords, (void(*)(void *))stKVDatabaseBulkResult_destruct);
	//startTransaction(database);
	stTry {
		for (int64_t i = 0; i < numRecords; ++i)
		{
			int64_t key = firstKey + i;
			int64_t recordSize;
			void* record = getRecord2(database, key, &recordSize);
			stKVDatabaseBulkResult* result = stKVDatabaseBulkResult_construct(record, recordSize);
			stList_set(results, i, result);
		}
		//commitTransaction(database);
	}stCatch(ex) {
		//abortTransaction(database);
		stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID, "tokyo cabinet bulk get records failed");
	}stTryEnd;

	return results;
}

static void removeRecord(stKVDatabase *database, int64_t key) {
    TCBDB *dbImpl = database->dbImpl;
    if (!tcbdbout(dbImpl, &key, sizeof(int64_t))) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Removing key/value to database error: %s", tcbdberrmsg(tcbdbecode(dbImpl)));
    }
}

static void bulkRemoveRecords(stKVDatabase *database, stList *records) {
    startTransaction(database);
    stTry {
        for(int64_t i=0; i<stList_length(records); i++) {
            stIntTuple *j = stList_get(records, i);
            removeRecord(database, stIntTuple_get(j, 0));
        }
        commitTransaction(database);
    }stCatch(ex) {
        abortTransaction(database);
        stThrowNewCause(
                ex,
                ST_KV_DATABASE_EXCEPTION_ID,
                "tokyo cabinet bulk remove records failed");
    }stTryEnd;
}

//initialisation function

void stKVDatabase_initialise_tokyoCabinet(stKVDatabase *database, stKVDatabaseConf *conf, bool create) {
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

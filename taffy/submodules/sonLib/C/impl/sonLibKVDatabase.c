/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibKVDatabase.c
 *
 *  Created on: 18-Aug-2010
 *      Author: benedictpaten
 */

#include "sonLibGlobalsInternal.h"
#include "sonLibKVDatabasePrivate.h"
#include "sonLibString.h"

const char *ST_KV_DATABASE_EXCEPTION_ID = "ST_KV_DATABASE_EXCEPTION";
const char *ST_KV_DATABASE_RETRY_TRANSACTION_EXCEPTION_ID =
        "ST_KV_DATABASE_RETRY_TRANSACTION_EXCEPTION_ID";

static bool isRetryExcept(stExcept *except) {
    return stExcept_idEq(except, ST_KV_DATABASE_RETRY_TRANSACTION_EXCEPTION_ID);
}

stKVDatabase *stKVDatabase_construct(stKVDatabaseConf *conf, bool create) {
    stKVDatabase *database = st_calloc(1, sizeof(struct stKVDatabase));
    database->conf = stKVDatabaseConf_constructClone(conf);
    database->deleted = false;

    switch (stKVDatabaseConf_getType(conf)) {
        case stKVDatabaseTypeTokyoCabinet:
#ifdef HAVE_TOKYO_CABINET
            stKVDatabase_initialise_tokyoCabinet(database, conf, create);
#else
            stThrowNew(
                    ST_KV_DATABASE_EXCEPTION_ID,
                    "requested Tokyo Cabinet database, however sonlib is not compiled with Tokyo Cabinet support");
#endif
            break;
        case stKVDatabaseTypeKyotoTycoon:
#ifdef HAVE_KYOTO_TYCOON
            stKVDatabase_initialise_kyotoTycoon(database, conf, create);
#else
            stThrowNew(
                    ST_KV_DATABASE_EXCEPTION_ID,
                    "requested Kyoto Tycoon database, however sonlib is not compiled with Kyoto Tycoon support");
#endif
            break;
        case stKVDatabaseTypeMySql:
#ifdef HAVE_MYSQL
            stKVDatabase_initialise_MySql(database, conf, create);
#else
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                    "requested MySQL database, however sonlib is not compiled with MySql support");
#endif
            break;
        case stKVDatabaseTypeRedis:
#ifdef HAVE_REDIS
            stKVDatabase_initialise_Redis(database, conf, create);
#else
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                    "requested Redis database, however sonlib is not compiled with Redis support");
#endif
            break;
        default:
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                    "BUG: unrecognized database type");
    }
    return database;
}

void stKVDatabase_destruct(stKVDatabase *database) {
    if (!database->deleted) {
        stTry {
                database->destruct(database);
            }stCatch(ex)
                {
                    stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_destruct failed");
                }stTryEnd;
    }
    stKVDatabaseConf_destruct(database->conf);
    free(database);
}

void stKVDatabase_deleteFromDisk(stKVDatabase *database) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to delete a database that has already been deleted");
    }
    stTry {
            database->deleteDatabase(database);
        }stCatch(ex)
            {
                stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID,
                        "stKVDatabase_deleteFromDisk failed");
            }stTryEnd;
    database->deleted = true;
}

bool stKVDatabase_containsRecord(stKVDatabase *database, int64_t key) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to check if a record is in a database that has been deleted");
    }
    bool containsRecord = 0;
    stTry {
            containsRecord = database->containsRecord(database, key);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_containsRecord key %lld failed",
                            (long long) key);
                }
            }stTryEnd;
    return containsRecord;
}

void stKVDatabase_insertRecord(stKVDatabase *database, int64_t key,
        const void *value, int64_t sizeOfRecord) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to insert a record into a database that has been deleted");
    }
    if(value == NULL) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                        "Trying to insert a null record into a database");
    }
    stTry {
            database->insertRecord(database, key, value, sizeOfRecord);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(
                            ex,
                            ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_insertRecord key %lld size %lld failed",
                            (long long) key, (long long) sizeOfRecord);
                }
            }stTryEnd;
}

void stKVDatabase_insertInt64(stKVDatabase *database, int64_t key, int64_t value) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to insert a int64 record into a database that has been deleted");
    }
    stTry {
            database->insertInt64(database, key, value);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(
                            ex,
                            ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_insertRecord key %lld size %d failed",
                            (long long) key, sizeof(int64_t));
                }
            }stTryEnd;
}

void stKVDatabase_updateInt64(stKVDatabase *database, int64_t key, int64_t value) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to update a int64 record into a database that has been deleted");
    }
    stTry {
            database->updateInt64(database, key, value);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(
                            ex,
                            ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_updateInt64 key %lld size %d failed",
                            (long long) key, sizeof(int64_t));
                }
            }stTryEnd;
}

void stKVDatabase_updateRecord(stKVDatabase *database, int64_t key,
        const void *value, int64_t sizeOfRecord) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to update a record in a database that has been deleted");
    }
    if(value == NULL) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                        "Trying to insert a null record into a database");
    }
    stTry {
            database->updateRecord(database, key, value, sizeOfRecord);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(
                            ex,
                            ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_updateRecord key %lld size %lld failed",
                            (long long) key, (long long) sizeOfRecord);
                }
            }stTryEnd;
}

void stKVDatabase_setRecord(stKVDatabase *database, int64_t key,
        const void *value, int64_t sizeOfRecord) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to get set a record from a database that has been deleted");
    }
    if(value == NULL) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                        "Trying to insert a null record into a database");
    }
    stTry {
            database->setRecord(database, key, value, sizeOfRecord);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_setRecord key %lld size %lld failed",
                            (long long) key, (long long) sizeOfRecord);
                }
            }stTryEnd;
}

int64_t stKVDatabase_incrementInt64(stKVDatabase *database, int64_t key,
        int64_t incrementAmount) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to increment a numerical record from a database that has been deleted");
    }
    stTry {
            return database->incrementInt64(database, key, incrementAmount);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(
                            ex,
                            ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_incrementRecord key: %lld increment size: %lld failed",
                            (long long) key, (long long) incrementAmount);
                }
            }stTryEnd;
    assert(0);
    return 0;
}

static stKVDatabaseBulkRequest *stKVDatabaseBulkRequest_construct(int64_t key,
        const void *value, int64_t sizeOfRecord,
        enum stKVDatabaseBulkRequestType type) {
    stKVDatabaseBulkRequest *buildRequest = st_malloc(
            sizeof(stKVDatabaseBulkRequest));
    buildRequest->key = key;
    if(value == NULL) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                        "Trying to insert a null record into a database");
    }
    buildRequest->value = memcpy(st_malloc(sizeOfRecord), value, sizeOfRecord);
    buildRequest->size = sizeOfRecord;
    buildRequest->type = type;
    return buildRequest;
}

stKVDatabaseBulkRequest *stKVDatabaseBulkRequest_constructInsertRequest(
        int64_t key, const void *value, int64_t sizeOfRecord) {
    return stKVDatabaseBulkRequest_construct(key, value, sizeOfRecord, INSERT);
}

stKVDatabaseBulkRequest *stKVDatabaseBulkRequest_constructUpdateRequest(
        int64_t key, const void *value, int64_t sizeOfRecord) {
    return stKVDatabaseBulkRequest_construct(key, value, sizeOfRecord, UPDATE);
}

stKVDatabaseBulkRequest *stKVDatabaseBulkRequest_constructSetRequest(
        int64_t key, const void *value, int64_t sizeOfRecord) {
    return stKVDatabaseBulkRequest_construct(key, value, sizeOfRecord, SET);
}

void stKVDatabaseBulkRequest_destruct(stKVDatabaseBulkRequest *record) {
    free(record->value);
    free(record);
}

void stKVDatabase_bulkSetRecords(stKVDatabase *database, stList *records) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to bulk set records from a database that has been deleted");
    }
    stTry {
            database->bulkSetRecords(database, records);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(
                            ex,
                            ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_bulkSetRecords with %" PRIi64 " records to update",
                            stList_length(records));
                }
            }stTryEnd;
}

void stKVDatabase_bulkRemoveRecords(stKVDatabase *database, stList *records) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to bulk remove records from a database that has been deleted");
    }

#ifndef NDEBUG
    /*This code to check the presence of records is very expensive*/
    for (int64_t i = 0; i < stList_length(records); i++) {
        int64_t key = stIntTuple_get(stList_get(records, i), 0);
        if (!stKVDatabase_containsRecord(database, key)) {
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                    "The key is not in the database which we aim to remove: %lli", key);
        }
    }
#endif
    stTry {
            database->bulkRemoveRecords(database, records);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(
                            ex,
                            ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_bulkRemoveRecords with %" PRIi64 " records to update",
                            stList_length(records));
                }
            }stTryEnd;
}

int64_t stKVDatabase_getNumberOfRecords(stKVDatabase *database) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to get the number of records from a database that has been deleted");
    }
    int64_t numRecs = 0;
    stTry {
            numRecs = database->numberOfRecords(database);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_getNumberOfRecords failed");
                }
            }stTryEnd;
    return numRecs;
}

void *stKVDatabase_getRecord(stKVDatabase *database, int64_t key) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to get a record from a database that has already been deleted");
    }
    void *data = NULL;
    stTry {
            data = database->getRecord(database, key);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_getRecord key %lld failed",
                            (long long) key);
                }
            }stTryEnd;
    return data;
}

int64_t stKVDatabase_getInt64(stKVDatabase *database, int64_t key) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to get a record from a database that has already been deleted");
    }
    int64_t value = -1;
    stTry {
            value = database->getInt64(database, key);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_getInt64 key %lld failed",
                            (long long) key);
                }
            }stTryEnd;
    return value;
}

void *stKVDatabase_getRecord2(stKVDatabase *database, int64_t key,
        int64_t *recordSize) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to get a record from a database that has already been deleted");
    }
    void *data = NULL;
    stTry {
            data = database->getRecord2(database, key, recordSize);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_getRecord2 key %lld failed",
                            (long long) key);
                }
            }stTryEnd;
    return data;
}

void *stKVDatabase_getPartialRecord(stKVDatabase *database, int64_t key,
        int64_t zeroBasedByteOffset, int64_t sizeInBytes, int64_t recordSize) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to get a record from a database that has already been deleted");
    }
    if (zeroBasedByteOffset < 0 || sizeInBytes < 0 || zeroBasedByteOffset
            + sizeInBytes > recordSize) {
        stThrowNew(
                ST_KV_DATABASE_EXCEPTION_ID,
                "Partial record retrieval to out of bounds memory, requested start: %lld, requested size: %lld, entry size: %lld",
                (long long) zeroBasedByteOffset, (long long) sizeInBytes,
                (long long) recordSize);
    }
    void *data = NULL;
    stTry {
            data = database->getPartialRecord(database, key,
                    zeroBasedByteOffset, sizeInBytes, recordSize);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(
                            ex,
                            ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_getPartialRecord key %lld offset %lld size %lld failed",
                            (long long) key, (long long) zeroBasedByteOffset,
                            (long long) sizeInBytes);
                }
            }stTryEnd;
    return data;
}

stKVDatabaseBulkResult *stKVDatabaseBulkResult_construct(void* value, int64_t sizeOfRecord)
{
	stKVDatabaseBulkResult* result = st_malloc(sizeof(stKVDatabaseBulkResult));
	result->value = value;
	result->size = sizeOfRecord;
	return result;
}

void* stKVDatabaseBulkResult_getRecord(stKVDatabaseBulkResult* bulkResult, int64_t *sizeOfRecord)
{
	*sizeOfRecord = bulkResult->size;
	return bulkResult->value;
}

void stKVDatabaseBulkResult_destruct(stKVDatabaseBulkResult* bulkResult)
{
	free(bulkResult->value);
	free(bulkResult);
}

stList *stKVDatabase_bulkGetRecords(stKVDatabase *database, stList* keys) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to get records from a database that has already been deleted");
    }
    assert(keys != NULL);
    if(stList_length(keys) == 0) {
        return stList_construct();
    }
    stList *resultsList = NULL;
    stTry {
    	resultsList = database->bulkGetRecords(database, keys);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_BulkGetRecords with %d records failed",
                            stList_length(keys));
                }
            }stTryEnd;
    return resultsList;
}

stList *stKVDatabase_bulkGetRecordsRange(stKVDatabase *database, int64_t firstKey, int64_t numRecords){
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to get records from a database that has already been deleted");
    }
    assert(numRecords > 0);
    stList *resultsList = NULL;
    stTry {
    	resultsList = database->bulkGetRecordsRange(database, firstKey, numRecords);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_BulkGetRecords2 with %lld records failed",
                            (long long int)numRecords);
                }
            }stTryEnd;
    return resultsList;
}

void stKVDatabase_removeRecord(stKVDatabase *database, int64_t key) {
    if (database->deleted) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "Trying to remove a record from a database that has already been deleted");
    }
    if (!stKVDatabase_containsRecord(database, key)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                "The key is not in the database: %lli", key);
    }
    stTry {
            database->removeRecord(database, key);
        }stCatch(ex)
            {
                if (isRetryExcept(ex)) {
                    stThrow(ex);
                } else {
                    stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID,
                            "stKVDatabase_removeRecord key %lld failed",
                            (long long) key);
                }
            }stTryEnd;
}

stKVDatabaseConf *stKVDatabase_getConf(stKVDatabase *database) {
    return database->conf;
}

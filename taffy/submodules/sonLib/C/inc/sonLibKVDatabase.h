/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef SONLIB_KV_DATABASE_H_
#define SONLIB_KV_DATABASE_H_

#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

// General database exception id 
extern const char *ST_KV_DATABASE_EXCEPTION_ID;

// Exception where transaction should be retried
extern const char *ST_KV_DATABASE_RETRY_TRANSACTION_EXCEPTION_ID;

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Database functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a non-relational database object, using the given configuration
 * information to connect to the database.  Create the database if create is
 * true, otherwise it must exist.
 */
stKVDatabase *stKVDatabase_construct(stKVDatabaseConf *conf, bool create);

/*
 * Destructs a database. If the destruction occurs during a transaction the transaction
 * is aborted and any changes are not committed to the database.
 * Also destroys any databases associated the the database.
 */
void stKVDatabase_destruct(stKVDatabase *database);

/*
 * Removes the database from the disk. Any further operations on the database will
 * create an exception, so you should generally destruct the object after calling this function.
 */
void stKVDatabase_deleteFromDisk(stKVDatabase *database);

/*
 * Returns non-zero if the database contains a record with the given key.
 */
bool stKVDatabase_containsRecord(stKVDatabase *database, int64_t key);

/*
 * Add a new key/value record into the table. Values can not be null. Throws an exception if unsuccessful.
 */
void stKVDatabase_insertRecord(stKVDatabase *database, int64_t key, const void *value, int64_t sizeOfRecord);

/*
 * Add a new int64 key/value record into the table. Throws an exception if unsuccessful.
 */
void stKVDatabase_insertInt64(stKVDatabase *database, int64_t key, int64_t value);

/*
 * Update an existing int64 key/value record into the table. Throws an exception if unsuccessful.
 */
void stKVDatabase_updateInt64(stKVDatabase *database, int64_t key, int64_t value);

/*
 * Get a int64 key/value record from the table. Throws an exception if unsuccessful.
 */
int64_t stKVDatabase_getInt64(stKVDatabase *database, int64_t key);

/*
 * Update an existing key/value record in the table. Values can not be null. Throws an exception if unsuccessful.
 */
void stKVDatabase_updateRecord(stKVDatabase *database, int64_t key, const void *value, int64_t sizeOfRecord);

/*
 * Update an existing key/value record in the table. Values can not be null. If the record does not exist it is inserted. Throws an exception if unsuccessful.
 */
void stKVDatabase_setRecord(stKVDatabase *database, int64_t key, const void *value, int64_t sizeOfRecord);

/*
 * Takes an existing record and treats it as type int64_t, incrementing the given amount in one atomic operation.
 * The result is the resulting incremented value.
 */
int64_t stKVDatabase_incrementInt64(stKVDatabase *database, int64_t key, int64_t incrementAmount);

/*
 * Creates an insert request. The value memory is copied and stored in the record, and only destroyed when freed (see stKVDatabaseBulkRequest_destruct).
 */
stKVDatabaseBulkRequest *stKVDatabaseBulkRequest_constructInsertRequest(int64_t key, const void *value, int64_t sizeOfRecord);

/*
 * Construct an update request (like stKVDatabaseBulkRequest_constructInsertRequest).
 */
stKVDatabaseBulkRequest *stKVDatabaseBulkRequest_constructUpdateRequest(int64_t key, const void *value, int64_t sizeOfRecord);

/*
 * Construct a set request (like stKVDatabaseBulkRequest_constructInsertRequest).
 */
stKVDatabaseBulkRequest *stKVDatabaseBulkRequest_constructSetRequest(int64_t key, const void *value, int64_t sizeOfRecord);

/*
 * Cleans up a kvDatabaseRecord object.
 */
void stKVDatabaseBulkRequest_destruct(stKVDatabaseBulkRequest *record);

/*
 * Updates in bulk a set of stKVDatabaseBulkRequests, represented in the list.
 * Throws a KV_DATABASE exception if unsuccessful.
 */
void stKVDatabase_bulkSetRecords(stKVDatabase *database, stList *records);

/*
 * Destruct a set of records, where each record is specified by a key, encoded in a stIntTuple.
 * Throws a KV_DATABASE exception if unsuccessful.
 */
void stKVDatabase_bulkRemoveRecords(stKVDatabase *database, stList *records);

/*
 * Gets a record from the database, given the key. The record is in newly allocated memory, and must be freed.
 * Returns NULL if the database does not contain the given record.
 */
void *stKVDatabase_getRecord(stKVDatabase *database, int64_t key);

/*
 * Gets a record from the database, given the key. The record is in newly allocated memory, and must be freed.
 * Returns NULL if the database does not contain the given record. Puts the size of the record in the in the record size field.
 */
void *stKVDatabase_getRecord2(stKVDatabase *database, int64_t key, int64_t *recordSize);


/*
 * Construct a bulk result (like a bulk request but keys are not stored)
 */
stKVDatabaseBulkResult *stKVDatabaseBulkResult_construct(void* value, int64_t sizeOfRecord);

/*
 * Get the record and size out of the bulkResult
 */
void* stKVDatabaseBulkResult_getRecord(stKVDatabaseBulkResult* bulkResult, int64_t *sizeOfRecord);

/*
 * Destruct the bulk result
 */
void stKVDatabaseBulkResult_destruct(stKVDatabaseBulkResult* bulkResult);

/*
 * Bulk get a batch of records (stored in a list of bulk results).  The nth result
 * corresponds to the nth key from the input keys list.
 */
stList *stKVDatabase_bulkGetRecords(stKVDatabase *database, stList* keys);

/*
 * Bulk get a batch of records (stored in a list of bulk results) from range of keys
 */
stList *stKVDatabase_bulkGetRecordsRange(stKVDatabase *database, int64_t firstKey, int64_t numRecords);


/*
 * Removes a record from the database. Throws an exception if unsuccessful.
 */
void stKVDatabase_removeRecord(stKVDatabase *database, int64_t key);

/*
 * Gets a record from the database, given the key. This function allows the partial retrieval of a record, using the
 * given offset and the size of the requested retrieval. The function should thrpw an exception if the record does
 * not exist (in contrast to stKVDatabase_getRecord) and similarly should throw an exception if the requested region lies
 * outside of the bounds of the record (the total size of the record must be passed the function).
 */
void *stKVDatabase_getPartialRecord(stKVDatabase *database, int64_t key, int64_t zeroBasedByteOffset, int64_t sizeInBytes, int64_t recordSize);

/*
 * Returns number of records in database.
 */
int64_t stKVDatabase_getNumberOfRecords(stKVDatabase *database);


/*
 * get the configuration object for the database.
 */
stKVDatabaseConf *stKVDatabase_getConf(stKVDatabase *database);


#ifdef __cplusplus
}
#endif
#endif

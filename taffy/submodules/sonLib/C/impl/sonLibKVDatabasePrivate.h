/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibKVDatabasePrivate.h
 *
 *  Created on: 18 Aug 2010
 *      Author: benedictpaten
 */

#ifndef SONLIBKVDATABASEPRIVATE_H_
#define SONLIBKVDATABASEPRIVATE_H_

struct stKVDatabase {
    stKVDatabaseConf *conf;
    void *dbImpl;
    struct stKVDatabase* secondaryDB;
    bool deleted;
    void (*destruct)(stKVDatabase *);
    void (*deleteDatabase)(stKVDatabase *);
    bool (*containsRecord)(stKVDatabase *, int64_t);
    void (*insertRecord)(stKVDatabase *, int64_t, const void *, int64_t);
    void (*insertInt64)(stKVDatabase *, int64_t, int64_t);
    void (*updateRecord)(stKVDatabase *, int64_t, const void *, int64_t);
    void (*updateInt64)(stKVDatabase *, int64_t, int64_t);
    void (*setRecord)(stKVDatabase *, int64_t, const void *, int64_t);
    int64_t (*incrementInt64)(stKVDatabase *, int64_t, int64_t);
    void (*bulkSetRecords)(stKVDatabase *, stList *);
    void (*bulkRemoveRecords)(stKVDatabase *, stList *);
    int64_t (*numberOfRecords)(stKVDatabase *);
    void *(*getRecord)(stKVDatabase *, int64_t key);
    int64_t (*getInt64)(stKVDatabase *, int64_t key);
    void *(*getRecord2)(stKVDatabase *database, int64_t key, int64_t *recordSize);
    void *(*getPartialRecord)(stKVDatabase *database, int64_t key, int64_t zeroBasedByteOffset, int64_t sizeInBytes, int64_t recordSize);
    stList *(*bulkGetRecords)(stKVDatabase *database, stList* keys);
    stList *(*bulkGetRecordsRange)(stKVDatabase *database, int64_t firstKey, int64_t numRecords);
    void (*removeRecord)(stKVDatabase *, int64_t key);
};

enum stKVDatabaseBulkRequestType {
    UPDATE, INSERT, SET
};

struct stKVDatabaseBulkRequest {
    int64_t key;
    void *value;
    int64_t size;
    enum stKVDatabaseBulkRequestType type;
};

struct stKVDatabaseBulkResult {
	void *value;
	int64_t size;
};

/*
 * Function initialises the pointers of the stKVDatabase object with functions for tokyoCabinet.
 */
void stKVDatabase_initialise_tokyoCabinet(stKVDatabase *database, stKVDatabaseConf *conf, bool create);

/*
 * Function initialises the pointers of the stKVDatabase object with functions for MySql.
 */
void stKVDatabase_initialise_MySql(stKVDatabase *database, stKVDatabaseConf *conf, bool create);

#ifdef __cplusplus
extern "C" {
#endif
void stKVDatabase_initialise_kyotoTycoon(stKVDatabase *database, stKVDatabaseConf *conf, bool create);
/*
 * Function initialises the pointers of the stKVDatabase object with functions for Big Record File.
 */
void stKVDatabase_initialise_bigRecordFile(stKVDatabase *database, stKVDatabaseConf *conf, bool create);

#ifdef __cplusplus
}
#endif

void stKVDatabase_initialise_Redis(stKVDatabase *database, stKVDatabaseConf *conf, bool create);

#endif /* SONLIBKVDATABASEPRIVATE_H_ */

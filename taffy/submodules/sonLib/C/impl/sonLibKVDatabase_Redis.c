#include "sonLibGlobalsInternal.h"
#include "sonLibKVDatabasePrivate.h"

#ifdef HAVE_REDIS
#include "hiredis.h"

// Thin wrapper around the redis context, in case we want to support
// something like auth, using a different DB number, etc.
typedef struct {
    redisContext *ctxt;
} RedisDB;

static redisReply *stRedisCommand(RedisDB *db, const char *string, ...) {
    va_list ap;
    va_start(ap, string);
    redisReply *reply = redisvCommand(db->ctxt, string, ap);
    va_end(ap);

    if (reply == NULL) {
        // Failure in executing the command--no reply received
        char *command = stString_print_r(string, ap);
        stExcept *except = stExcept_new(ST_KV_DATABASE_EXCEPTION_ID,
                                        "Failed to execute the Redis command '%s': %s",
                                        command, db->ctxt->errstr);
        free(command);
        stThrow(except);
    }
    if (reply->type == REDIS_REPLY_ERROR) {
        char *command = stString_print_r(string, ap);
        stExcept *except = stExcept_new(ST_KV_DATABASE_EXCEPTION_ID,
                                        "Failed to execute the Redis command '%s': %s",
                                        command, reply->str);
        free(command);
        freeReplyObject(reply);
        stThrow(except);
    }
    return reply;
}

static redisReply *stRedisGetReply(RedisDB *db) {
    redisReply *reply;
    int err = redisGetReply(db->ctxt, (void **) &reply);
    if (err == REDIS_ERR || reply == NULL) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Failed to get reply: %s", db->ctxt->errstr);
    }
    if (reply->type == REDIS_REPLY_ERROR) {
        stExcept *except = stExcept_new(ST_KV_DATABASE_EXCEPTION_ID, "Got a failure reply: %s",
                                        reply->str);
        freeReplyObject(reply);
        stThrow(except);
    }
    return reply;
}

/* connect to a database server */
static RedisDB *connect(stKVDatabaseConf *conf) {
    RedisDB *db = st_calloc(1, sizeof(RedisDB));

    db->ctxt = redisConnect(stKVDatabaseConf_getHost(conf), stKVDatabaseConf_getPort(conf));
    if (db->ctxt != NULL && db->ctxt->err) {
        stExcept *except = stExcept_new(ST_KV_DATABASE_EXCEPTION_ID, "redis connect failed: %s", db->ctxt->errstr);
        redisFree(db->ctxt);
        free(db);
        stThrow(except);
    }
    return db;
}

static void destructDB(stKVDatabase *database) {
    RedisDB *db = database->dbImpl;
    redisFree(db->ctxt);
    free(db);
    database->dbImpl = NULL;
}

static void deleteDB(stKVDatabase *database) {
    redisReply *reply = stRedisCommand(database->dbImpl, "FLUSHDB");
    freeReplyObject(reply);
    destructDB(database);
}

static void insertInt64(stKVDatabase *database, int64_t key, int64_t value) {
    redisReply *reply = stRedisCommand(database->dbImpl, "SET %" PRIi64 " %" PRIi64, key, value);
    freeReplyObject(reply);
}

static void updateInt64(stKVDatabase *database, int64_t key, int64_t value) {
    // NB: this doesn't quite follow the "update" semantics: it will
    // insert the key if it doesn't exist. This could be fixed by
    // using a lua script.
    insertInt64(database, key, value);
}

static void insertRecord(stKVDatabase *database, int64_t key, const void *value, int64_t sizeOfRecord) {
    redisReply *reply = stRedisCommand(database->dbImpl, "SET %" PRIi64 " %b", key,
                                       value, sizeOfRecord);
    freeReplyObject(reply);
}

static void updateRecord(stKVDatabase *database, int64_t key,
                         const void *value, int64_t sizeOfRecord) {
    // See caveats about the "update" in the comment to  updateInt64 above.
    insertRecord(database, key, value, sizeOfRecord);
}

static int64_t numberOfRecords(stKVDatabase *database) {
    redisReply *reply = stRedisCommand(database->dbImpl, "DBSIZE");
    assert(reply->type == REDIS_REPLY_INTEGER);
    int64_t ret = reply->integer;
    freeReplyObject(reply);
    return ret;
}

static void *getRecord2(stKVDatabase *database, int64_t key, int64_t *sizeOfRecord) {
    redisReply *reply = stRedisCommand(database->dbImpl, "GET %" PRIi64, key);
    if (reply->type == REDIS_REPLY_NIL) {
        // No key in DB.
        freeReplyObject(reply);
        return NULL;
    } else {
        assert(reply->type == REDIS_REPLY_STRING);
        void *ret = malloc(reply->len);
        memcpy(ret, reply->str, reply->len);
        if (sizeOfRecord != NULL) {
            *sizeOfRecord = reply->len;
        }
        freeReplyObject(reply);
        return ret;
    }
}

static int64_t getInt64(stKVDatabase *database, int64_t key) {
    int64_t len = 0;
    void *record = getRecord2(database, key, &len);
    // Make a null-terminated string
    char *str = malloc(len + 1);
    memcpy(str, record, len);
    free(record);
    str[len] = '\0';
    int64_t ret;
    if (sscanf(str, "%" PRIi64, &ret) == 0) {
        free(str);
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Key %" PRIi64 " does not have an integer value", key);
    }
    free(str);
    return ret;
}

static void *getRecord(stKVDatabase *database, int64_t key) {
    return getRecord2(database, key, NULL);
}

static bool containsRecord(stKVDatabase *database, int64_t key) {
    redisReply *reply = stRedisCommand(database->dbImpl, "EXISTS %" PRIi64, key);
    assert(reply->type == REDIS_REPLY_INTEGER);
    bool exists;
    if (reply->integer) {
        exists = true;
    } else {
        exists = false;
    }
    freeReplyObject(reply);
    return exists;
}

static void *getPartialRecord(stKVDatabase *database, int64_t key, int64_t zeroBasedByteOffset, int64_t sizeInBytes, int64_t recordSize) {
    redisReply *reply = stRedisCommand(database->dbImpl, "GETRANGE %" PRIi64 " %" PRIi64 " %" PRIi64,
                                       key, zeroBasedByteOffset,
                                       zeroBasedByteOffset + sizeInBytes - 1);
    assert(reply->type == REDIS_REPLY_STRING);
    if (reply->len != sizeInBytes) {
        int64_t len = reply->len;
        freeReplyObject(reply);
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "partial read of key %lld, expected %lld bytes got %lld bytes", key, sizeInBytes, len);
    }
    char *data = malloc(reply->len * sizeof(char));
    memcpy(data, reply->str, reply->len);
    freeReplyObject(reply);
    return data;
}

static void removeRecord(stKVDatabase *database, int64_t key) {
    redisReply *reply = stRedisCommand(database->dbImpl, "DEL %" PRIi64, key);
    freeReplyObject(reply);
}

static int64_t incrementInt64(stKVDatabase *database, int64_t key, int64_t incrementAmount) {
    redisReply *reply = stRedisCommand(database->dbImpl, "INCRBY %" PRIi64 " %" PRIi64, key, incrementAmount);
    assert(reply->type == REDIS_REPLY_INTEGER);
    int64_t ret = reply->integer;
    freeReplyObject(reply);
    return ret;
}

static void bulkRemoveRecords(stKVDatabase *database, stList *records) {
    RedisDB *db = database->dbImpl;
    redisAppendCommand(db->ctxt, "MULTI");
    for(int64_t i = 0; i < stList_length(records); i++) {
        stIntTuple *box = stList_get(records, i);
        redisAppendCommand(db->ctxt, "DEL %" PRIi64, stIntTuple_get(box, 0));
    }
    redisAppendCommand(db->ctxt, "EXEC");
    redisReply *reply;
    reply = stRedisGetReply(db);
    freeReplyObject(reply);
    for(int64_t i = 0; i < stList_length(records); i++) {
        reply = stRedisGetReply(db);
        freeReplyObject(reply);
    }
    reply = stRedisGetReply(db);
    freeReplyObject(reply);
}

static void setRecord(stKVDatabase *database, int64_t key,
                      const void *value, int64_t sizeOfRecord) {
    insertRecord(database, key, value, sizeOfRecord);
}

static stList *bulkGetRecords(stKVDatabase *database, stList* keys) {
    RedisDB *db = database->dbImpl;
    for (int64_t i = 0; i < stList_length(keys); i++) {
        int64_t key = *((int64_t *) stList_get(keys, i));
        redisAppendCommand(db->ctxt, "GET %" PRIi64, key);
    }
    stList* results = stList_construct3(stList_length(keys),
                                        (void(*)(void *))stKVDatabaseBulkResult_destruct);
    for(int64_t i = 0; i < stList_length(keys); i++) {
        redisReply *reply = stRedisGetReply(db);
        void *record;
        int64_t length;
        if (reply->type == REDIS_REPLY_NIL) {
            // No record found
            record = NULL;
            length = 0;
        } else {
            // Found a record
            assert(reply->type == REDIS_REPLY_STRING);
            record = malloc(reply->len);
            memcpy(record, reply->str, reply->len);
            length = reply->len;
        }
        stKVDatabaseBulkResult* result = stKVDatabaseBulkResult_construct(record, length);
        stList_set(results, i, result);
        freeReplyObject(reply);
    }
    return results;
}

static stList *bulkGetRecordsRange(stKVDatabase *database, int64_t firstKey, int64_t numRecords) {
    // Just fake this by creating the list ahead of time.
    stList *list = stList_construct3(numRecords, free);
    for (int64_t i = 0; i < numRecords; i++) {
        int64_t *n = malloc(sizeof(int64_t));
        *n = firstKey + i;
        stList_set(list, i, n);
    }
    stList *ret = bulkGetRecords(database, list);
    stList_destruct(list);
    return ret;
}

static void bulkSetRecords(stKVDatabase *database, stList *records) {
    RedisDB *db = database->dbImpl;
    stKVDatabaseConf* conf = stKVDatabase_getConf(database);
    int64_t maxBulkSetSize = stKVDatabaseConf_getMaxRedisBulkSetSize(conf);
    int64_t runningSize = 0;
    int64_t countAddedRecords = 0;
    redisAppendCommand(db->ctxt, "MULTI");
    for(int64_t i = 0; i < stList_length(records); i++) {
        stKVDatabaseBulkRequest *request = stList_get(records, i);
	if ((runningSize + request->size) > maxBulkSetSize){
		redisAppendCommand(db->ctxt, "EXEC");
    		redisReply *reply;
    		reply = stRedisGetReply(db);
    		freeReplyObject(reply);
    		for(int64_t i = 0; i < countAddedRecords; i++) {
        		reply = stRedisGetReply(db);
        		freeReplyObject(reply);
    		}
    		reply = stRedisGetReply(db);
    		freeReplyObject(reply);
		runningSize = 0;
		countAddedRecords = 0;
		redisAppendCommand(db->ctxt, "MULTI");
	}
        redisAppendCommand(db->ctxt, "SET %" PRIi64 " %b", request->key, request->value, request->size);
	runningSize += request->size;
	countAddedRecords += 1;
    }
    redisAppendCommand(db->ctxt, "EXEC");
    redisReply *reply;
    reply = stRedisGetReply(db);
    freeReplyObject(reply);
    for(int64_t i = 0; i < countAddedRecords; i++) {
        reply = stRedisGetReply(db);
        freeReplyObject(reply);
    }
    reply = stRedisGetReply(db);
    freeReplyObject(reply);
}

//initialisation function
void stKVDatabase_initialise_Redis(stKVDatabase *database, stKVDatabaseConf *conf, bool create) {
    database->dbImpl = connect(conf);
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

#endif // HAVE_REDIS

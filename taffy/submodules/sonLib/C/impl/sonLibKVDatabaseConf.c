/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibKVDatabaseConf.c
 *
 *  Created on: 2010-08-23
 *      Author: Mark Diekhans
 */

#include "sonLibKVDatabaseConf.h"
#include "sonLibKVDatabase.h"
#include "sonLibExcept.h"
#include "sonLibTypes.h"
#include "sonLibString.h"
#include "stSafeC.h"
#include "sonLib.h"

/* object storing the database configuration */
struct stKVDatabaseConf {
    stKVDatabaseType type;
    char *databaseDir;
    char *host;
    unsigned port;
    int timeout;
    int64_t maxKTRecordSize;
    int64_t maxKTBulkSetSize;
    int64_t maxRedisBulkSetSize;
    int64_t maxKTBulkSetNumRecords;
    char *user;
    char *password;
    char *databaseName;
    char *tableName;
};

stKVDatabaseConf *stKVDatabaseConf_constructTokyoCabinet(const char *databaseDir) {
    stKVDatabaseConf *conf = stSafeCCalloc(sizeof(stKVDatabaseConf));
    conf->type = stKVDatabaseTypeTokyoCabinet;
    conf->databaseDir = stString_copy(databaseDir);
    return conf;
}

stKVDatabaseConf *stKVDatabaseConf_constructKyotoTycoon(const char *host, unsigned port, int timeout,
                                                        int64_t maxRecordSize, int64_t maxBulkSetSize,
                                                        int64_t maxBulkSetNumRecords,
                                                        const char *databaseDir, const char* databaseName) {
    stKVDatabaseConf *conf = stSafeCCalloc(sizeof(stKVDatabaseConf));
    conf->type = stKVDatabaseTypeKyotoTycoon;
    conf->databaseDir = stString_copy(databaseDir);
    conf->host = stString_copy(host);
    conf->port = port;
    conf->timeout = timeout;
    conf->maxKTRecordSize = maxRecordSize;
    conf->maxKTBulkSetSize = maxBulkSetSize;
    conf->maxKTBulkSetNumRecords = maxBulkSetNumRecords;
    conf->databaseName = stString_copy(databaseName);
    return conf;
}

stKVDatabaseConf *stKVDatabaseConf_constructRedis(const char *host, unsigned port, int64_t maxRedisBulkSetSize) {
    stKVDatabaseConf *conf = stSafeCCalloc(sizeof(stKVDatabaseConf));
    conf->type = stKVDatabaseTypeRedis;
    conf->host = stString_copy(host);
    conf->port = port;
    conf->maxRedisBulkSetSize = maxRedisBulkSetSize;
    return conf;
}

static stKVDatabaseConf *constructSql(stKVDatabaseType type, const char *host, unsigned port, const char *user, const char *password,
                                      const char *databaseName, const char *tableName) {
    stKVDatabaseConf *conf = stSafeCCalloc(sizeof(stKVDatabaseConf));
    conf->type = type;
    conf->host = stString_copy(host);
    conf->port = port;
    conf->user = stString_copy(user);
    conf->password = stString_copy(password);
    conf->databaseName = stString_copy(databaseName);
    conf->tableName = stString_copy(tableName);
    return conf;
}

stKVDatabaseConf *stKVDatabaseConf_constructMySql(const char *host, unsigned port, const char *user, const char *password,
                                                  const char *databaseName, const char *tableName) {
#ifndef HAVE_MYSQL
    stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
            "requested MySQL database, however sonlib is not compiled with MySql support");
#endif
    return constructSql(stKVDatabaseTypeMySql, host, port, user, password, databaseName, tableName);
}

char *getNextToken(char **tokenStream) {
    char *token = stString_getNextWord(tokenStream);
    if(token != NULL) {
        char *cA = stString_replace(token, "\"", "");
        free(token);
        token = stString_replace(cA, "'", "");
        free(cA);
    }
    return token;
}

static void getExpectedToken(char **tokenStream, const char *expected) {
    char *cA = getNextToken(tokenStream);
    if (!stString_eq(cA, expected)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
                   "BUG: expected the token: %s in database XML string, but I got: %s from the stream %s",
                   expected, cA, *tokenStream);
    }
    free(cA);
}

static char *getKeyValue(char **tokenStream, const char *key) {
    getExpectedToken(tokenStream, key);
    return getNextToken(tokenStream);
}

/* Parse XML string into a hash.  This parses all attributes of all tags
 * into values.  st_kv_database_conf type is stored as conf_type,
 * database tag is stores as db_tag.  This does minimal error checking
 * and is really lame.
 */
static stHash *hackParseXmlString(const char *xmlString) {
    stHash *hash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    char *toReplace[5] = { "</", "<", "/>", ">", "=" };
    char *cA = stString_replace(xmlString, toReplace[0], " "), *cA2;
    for (int64_t i = 1; i < 5; i++) {
        cA2 = stString_replace(cA, toReplace[i], " ");
        free(cA);
        cA = cA2;
    }
    getExpectedToken(&cA2, "st_kv_database_conf");
    stHash_insert(hash, stString_copy("conf_type"), getKeyValue(&cA2, "type"));
    stHash_insert(hash, stString_copy("db_tag"), getNextToken(&cA2));

    char *key;
    while (((key = getNextToken(&cA2)) != NULL) && !stString_eq(key, "st_kv_database_conf")) {
        char *value = getNextToken(&cA2);
        if (value == NULL) {
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "failed to to get value for key \"%s\"", key);
        }
        if (stHash_search(hash, key) != NULL) {
            stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "got a duplicate entry in the database conf string \"%s\"", key);
        }
        stHash_insert(hash, key, value);
    }
    if(!stString_eq(key, "st_kv_database_conf")) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "got an unexpected final entry \"%s\"", key);
    }
    free(key);
    free(cA);
    return hash;
}

static const char *getXmlValueRequired(stHash *hash, const char *key) {
    const char *value = stHash_search(hash, (char*)key);
    if (value == NULL) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "did not find a \"%s\" value in the database XML string", key);
    }
    return value;
}

static int getXmlPort(stHash *hash) {
    const char *value = stHash_search(hash, "port");
    if (value == NULL) {
        return 0;
    } else {
        return stSafeStrToUInt32(value);
    }
}

static int getXmlTimeout(stHash *hash) {
    const char *value = stHash_search(hash, "timeout");
    if (value == NULL) {
        // default to -1 -- meaning no timeout
        return -1;
    } else {
        return stSafeStrToUInt32(value);
    }
}

/* Default to 50M.  It used to be 175M but since noticing
 * problems in bulk *get* within cactus secondary dbs we 
 * crank it way down.  Unlike bulk set, in the get we don't
 * know the total size of the requested records.  To prevent
 * big data transfers (which are problematic), we reduce max 
 * size of individual records and hope for the best 
 */
static int64_t getXMLMaxKTRecordSize(stHash *hash) {
    const char *value = stHash_search(hash, "max_record_size");
    if (value == NULL) {
        return (int64_t) 10000000;
    } else {
        return stSafeStrToInt64(value);
    }
}

/* Default to 175M which seems to be about where the
 * kyoto tycoon network error danger zone starts
 */
static int64_t getXMLMaxKTBulkSetSize(stHash *hash) {
    const char *value = stHash_search(hash, "max_bulkset_size");
    if (value == NULL) {
        return (int64_t) 183500800;
    } else {
        return stSafeStrToInt64(value);
    }
}

/* Default to tried-and-true value of 10000
 */
static int64_t getXMLMaxKTBulkSetNumRecords(stHash *hash) {
    const char *value = stHash_search(hash, "max_bulkset_num_records");
    if (value == NULL) {
        return (int64_t) 10000;
    } else {
        return stSafeStrToInt64(value);
    }
}

/* Default to 1.75GB which seems to be about where the
 * redis network error danger zone starts
 */
static int64_t getXMLMaxRedisBulkSetSize(stHash *hash) {
    const char *value = stHash_search(hash, "max_bulkset_size");
    if (value == NULL) {
        return (int64_t) 175000000000;
    } else {
        return stSafeStrToInt64(value);
    }
}

static stKVDatabaseConf *constructFromString(const char *xmlString) {
    stHash *hash = hackParseXmlString(xmlString);
    stKVDatabaseConf *databaseConf = NULL;
    const char *type = getXmlValueRequired(hash, "conf_type");
    const char *dbTag = getXmlValueRequired(hash, "db_tag");
    if (!stString_eq(type, dbTag)) {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Database XML tag \"%s\" did not match st_kv_database_conf type attribute", dbTag, type);
    }
    if (stString_eq(type, "tokyo_cabinet")) {
        databaseConf = stKVDatabaseConf_constructTokyoCabinet(getXmlValueRequired(hash, "database_dir"));
    } else if (stString_eq(type, "kyoto_tycoon")) {
        databaseConf = stKVDatabaseConf_constructKyotoTycoon(getXmlValueRequired(hash, "host"), 
                                                        getXmlPort(hash), 
                                                        getXmlTimeout(hash), 
                                                        getXMLMaxKTRecordSize(hash),
                                                        getXMLMaxKTBulkSetSize(hash),
                                                        getXMLMaxKTBulkSetNumRecords(hash),
                                                        getXmlValueRequired(hash, "database_dir"),
                                                        stHash_search(hash, "database_name"));
    } else if (stString_eq(type, "mysql")) {
        databaseConf = stKVDatabaseConf_constructMySql(getXmlValueRequired(hash, "host"), getXmlPort(hash),
                                                       getXmlValueRequired(hash, "user"), getXmlValueRequired(hash, "password"),
                                                       getXmlValueRequired(hash, "database_name"), getXmlValueRequired(hash, "table_name"));
    } else if (stString_eq(type, "redis")) {
        databaseConf = stKVDatabaseConf_constructRedis(getXmlValueRequired(hash, "host"), getXmlPort(hash),
						       getXMLMaxRedisBulkSetSize(hash));
    } else {
        stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "invalid database type \"%s\"", type);
    }
    stHash_destruct(hash);
    return databaseConf;
}

stKVDatabaseConf *stKVDatabaseConf_constructFromString(const char *xmlString) {
    stKVDatabaseConf *conf = NULL;
    stTry {
        conf = constructFromString(xmlString);
    } stCatch(ex) {
        stThrowNewCause(ex, ST_KV_DATABASE_EXCEPTION_ID, "Invalid database XML specification: %s", xmlString);
    } stTryEnd;
    return conf;
}


stKVDatabaseConf *stKVDatabaseConf_constructClone(stKVDatabaseConf *srcConf) {
    stKVDatabaseConf *conf = stSafeCCalloc(sizeof(stKVDatabaseConf));
    conf->type = srcConf->type;
    conf->databaseDir = stString_copy(srcConf->databaseDir);
    conf->host = stString_copy(srcConf->host);
    conf->port = srcConf->port;
    conf->timeout = srcConf->timeout;
    conf->maxKTRecordSize = srcConf->maxKTRecordSize;
    conf->maxKTBulkSetSize = srcConf->maxKTBulkSetSize;
    conf->maxRedisBulkSetSize = srcConf->maxRedisBulkSetSize;
    conf->maxKTBulkSetNumRecords = srcConf->maxKTBulkSetNumRecords;
    conf->user = stString_copy(srcConf->user);
    conf->password = stString_copy(srcConf->password);
    conf->databaseName = stString_copy(srcConf->databaseName);
    conf->tableName = stString_copy(srcConf->tableName);
    return conf;
}

void stKVDatabaseConf_destruct(stKVDatabaseConf *conf) {
    if (conf != NULL) {
        stSafeCFree(conf->databaseDir);
        stSafeCFree(conf->host);
        stSafeCFree(conf->user);
        stSafeCFree(conf->password);
        stSafeCFree(conf->databaseName);
        stSafeCFree(conf->tableName);
        stSafeCFree(conf);
    }
}

stKVDatabaseType stKVDatabaseConf_getType(stKVDatabaseConf *conf) {
    return conf->type;
}

const char *stKVDatabaseConf_getDir(stKVDatabaseConf *conf) {
    return conf->databaseDir;
}

const char *stKVDatabaseConf_getHost(stKVDatabaseConf *conf) {
    return conf->host;
}

unsigned stKVDatabaseConf_getPort(stKVDatabaseConf *conf) {
    return conf->port;
}

int stKVDatabaseConf_getTimeout(stKVDatabaseConf *conf) {
    return conf->timeout;
}

int64_t stKVDatabaseConf_getMaxKTRecordSize(stKVDatabaseConf *conf) {
    return conf->maxKTRecordSize;
}

void stKVDatabaseConf_setMaxKTRecordSize(stKVDatabaseConf *conf,
                                         int64_t maxRecordSize) {
    conf->maxKTRecordSize = maxRecordSize;
}

int64_t stKVDatabaseConf_getMaxKTBulkSetSize(stKVDatabaseConf *conf) {
    return conf->maxKTBulkSetSize;
}

int64_t stKVDatabaseConf_getMaxKTBulkSetNumRecords(stKVDatabaseConf *conf) {
    return conf->maxKTBulkSetNumRecords;
}

int64_t stKVDatabaseConf_getMaxRedisBulkSetSize(stKVDatabaseConf *conf) {
    return conf->maxRedisBulkSetSize;
}

const char *stKVDatabaseConf_getUser(stKVDatabaseConf *conf) {
    return conf->user;
}

const char *stKVDatabaseConf_getPassword(stKVDatabaseConf *conf) {
    return conf->password;
}

const char *stKVDatabaseConf_getDatabaseName(stKVDatabaseConf *conf) {
    return conf->databaseName;
}

const char *stKVDatabaseConf_getTableName(stKVDatabaseConf *conf) {
    return conf->tableName;
}


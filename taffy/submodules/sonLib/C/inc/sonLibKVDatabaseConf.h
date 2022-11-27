/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef SONLIBKVDATABASECONF_H
#define SONLIBKVDATABASECONF_H
#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    stKVDatabaseTypeTokyoCabinet,
    stKVDatabaseTypeKyotoTycoon,
    stKVDatabaseTypeMySql,
    stKVDatabaseTypeRedis,
} stKVDatabaseType;

/* 
 * Construct a new database configuration object for a Tokyo Cabinet
 * database.
 */
stKVDatabaseConf *stKVDatabaseConf_constructTokyoCabinet(const char *databaseDir);

/* 
 * Construct a new database configuration object for a Kyoto Tycoon
 * database remote object.
 */
stKVDatabaseConf *stKVDatabaseConf_constructKyotoTycoon(const char *host, unsigned port, int timeout,
                                                        int64_t maxRecordSize, int64_t maxBulkSetSize,
                                                        int64_t maxBulkSetNumRecords,
                                                        const char *databaseDir, const char* databaseName);

/* 
 * Construct a new database configuration object for a Redis
 * database remote object.
 */
stKVDatabaseConf *stKVDatabaseConf_constructRedis(const char *host, unsigned port, int64_t maxBulkSetSize);

/* 
 * Construct a new database configuration object for a MySql database.
 * password maybe NULL for no password.
 * port maybe 0 for the default port.
 */
stKVDatabaseConf *stKVDatabaseConf_constructMySql(const char *host, unsigned port, const char *user, const char *password,
                                                  const char *databaseName, const char *tableName);

/*
 * Decodes a simple piece of XML, structured as follows:
 * <st_kv_database_conf type="TYPE">
 *      <tokyo_cabinet database_dir=""/>
 *      <mysql host="" port="" user="" password="" database_name="" table_name=""/>
 *      <kyoto_cabinet host="" port=""/>
 * </st_kv_database_conf>
 *
 * Type can be "tokyo_cabinet", "mysql", or "kyoto_cabinet". If it is of that type then
 * you need to include a nested tag with the parameters for that conf constructor.
 * The labels for the nested tag are name value pairs (no order assumed) for the conf constructor
 * (see above).  The port is optional.
 */
stKVDatabaseConf *stKVDatabaseConf_constructFromString(const char *xmlString);

/* 
 * Construct a new database configuration from an existing one.
 */
stKVDatabaseConf *stKVDatabaseConf_constructClone(stKVDatabaseConf *srcConf);

/*
 * Free the object.
 */
void stKVDatabaseConf_destruct(stKVDatabaseConf *conf);

/* 
 * get the database type.
 */
stKVDatabaseType stKVDatabaseConf_getType(stKVDatabaseConf *conf);

/* get the directory for file based databases */
const char *stKVDatabaseConf_getDir(stKVDatabaseConf *conf);

/* get the host for server based databases */
const char *stKVDatabaseConf_getHost(stKVDatabaseConf *conf);

/* get the port for server based databases */
unsigned stKVDatabaseConf_getPort(stKVDatabaseConf *conf);

/* get the remote server timeout for server based databases */
int stKVDatabaseConf_getTimeout(stKVDatabaseConf *conf);

/* get the maximum size in bytes of a kyoto tycoon record */
int64_t stKVDatabaseConf_getMaxKTRecordSize(stKVDatabaseConf *conf);

/* set the maximum size in bytes of a kyoto tycoon record */
void stKVDatabaseConf_setMaxKTRecordSize(stKVDatabaseConf *conf,
                                         int64_t maxRecordSize);

/* get the maximum size in bytes of a kyoto tycoon bulk set */
int64_t stKVDatabaseConf_getMaxKTBulkSetSize(stKVDatabaseConf *conf);

/* get the maximum number of records in  kyoto tycoon bulk set */
int64_t stKVDatabaseConf_getMaxKTBulkSetNumRecords(stKVDatabaseConf *conf);

/* get the maximum size in bytes of a redis bulk set */
int64_t stKVDatabaseConf_getMaxRedisBulkSetSize(stKVDatabaseConf *conf);

/* get the user for server based databases */
const char *stKVDatabaseConf_getUser(stKVDatabaseConf *conf);

/* get the password for server based databases */
const char *stKVDatabaseConf_getPassword(stKVDatabaseConf *conf);

/* get the SQL database for server based databases */
const char *stKVDatabaseConf_getDatabaseName(stKVDatabaseConf *conf);

/* get the table name for server based databases */
const char *stKVDatabaseConf_getTableName(stKVDatabaseConf *conf);

#ifdef __cplusplus
}
#endif
#endif

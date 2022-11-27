/*
 * Copyright (C) 2006-2012 by Glenn Hickey
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * Write records directly to binary files, in a one record per
 * file scheme.  Designed to be used in conjunction with Kyoto Tycoon
 * as a work-around for "network errors" that invariably arise when
 * trying to write large records (>200MB) to the database.
 *
 * Doesn't fully implement the sonLib database interface (and is not
 * made public) but is consistent enough that it could be if needed...
 */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <limits.h>
#include "sonLibGlobalsInternal.h"
#include "sonLibKVDatabasePrivate.h"
#include "sonLibSortedSet.h"
#include "stSafeC.h"
#include "sonLibTuples.h"

/*
 * tag used to construct files for storing big records.  these all just
 * get dumped in the database dir from the conf.  note that this string
 * is used, hackily, by cactus_progressive when deleting databases..
 */
#define RECORD_FILE_TAG "BIG__RECORD__FILE__"

/*
 * Maximum number of records stored in the base directory.  If there
 * are too many files, it's probably time for a rethink.  The value below
 * is relatively high mostly just for unit testing.
 */
#define MAX_NUMBER_ENTRIES 2048

#define MAXIMUM_PATH_LENGTH 4096
/*
 * Get the full path of the file corresponding to the record
 * with the given key.  NEEDS TO BE FREED!!
 */
static char* createRecordPath(stKVDatabaseConf* conf, int64_t key)
{
	const char *basePath = stKVDatabaseConf_getDir(conf);
	const char *name = stKVDatabaseConf_getDatabaseName(conf);
	if (name == NULL)
		name = "db";
	int64_t size = strlen(basePath) + strlen(name) +
			strlen(RECORD_FILE_TAG) + 65;
	char* recordPath = (char*)malloc(size * sizeof(char));
	sprintf(recordPath, "%s/%s.%s%lld", basePath, name, RECORD_FILE_TAG,
			(long long int)key);
	return recordPath;
}

/*
 * visit every record in the directory.  currently implemented
 * as just a flat directory where records are stored in files
 * including the RECORD_FILE_TAG.
 */
static size_t visitRecords(const char* basePath,
		void (*fileFn)(const char* recordPath, void* arg), void* argument)
{
	DIR *dp;
	struct dirent *ep;
	char filePathBuffer[MAXIMUM_PATH_LENGTH];
	dp = opendir(basePath);
	int64_t numRecords = 0;
	if (dp != NULL)
	{
		while ((ep = readdir(dp)) != NULL)
		{
			if (strstr(ep->d_name, RECORD_FILE_TAG) != NULL)
			{
				sprintf(filePathBuffer, "%s/%s", basePath, ep->d_name);
				if (fileFn != NULL)
				{
					fileFn(filePathBuffer, argument);
				}
				++numRecords;
			}
		}
		(void)closedir(dp);
	}
	return numRecords;
}

/*
 * add a path to a sorted set
 */
static void add_to_sortedSet(const char* recordPath, void* arg)
{
	stSortedSet* sortedSet = (stSortedSet*)arg;
	const char* fileName = strstr(recordPath, RECORD_FILE_TAG);
	assert (fileName != NULL);
	const char* keyString = fileName + strlen(RECORD_FILE_TAG);
	assert (keyString != NULL);
	int64_t key = atol(keyString);
	stSortedSet_insert(sortedSet, stIntTuple_construct1( key));
}

/* get around complicated casting
 */
static void remove_with_arg(const char* path, void* arg)
{
	remove(path);
}

/*
 * build a sorted set of records in the given directory
 */
static stSortedSet* constructDB(stKVDatabaseConf *conf, bool create)
{
	const char *basePath = stKVDatabaseConf_getDir(conf);
    mkdir(basePath, S_IRWXU);
    stSortedSet* sortedSet = stSortedSet_construct3(
    		(int (*)(const void *, const void *))stIntTuple_cmpFn,
    		(void (*)(void *))stIntTuple_destruct);
    if (create == true)
    {
    	visitRecords(basePath, remove_with_arg, NULL);
    }
    else
    {
    	visitRecords(basePath, add_to_sortedSet, sortedSet);
    }
    return sortedSet;
}

/*
 * database in memory is just a sorted set, so we destroy that.
 */
static void destructDB(stKVDatabase *database)
{
	stSortedSet* sortedSet  = (stSortedSet*)database->dbImpl;
    if (sortedSet != NULL)
    {
    	stSortedSet_destruct(sortedSet);
    }
    database->dbImpl = NULL;
}

/* delete every file beginning with RECORD_FILE_TAG
 * in the database directory
 */
static void deleteDB(stKVDatabase *database)
{
	/* destruct called to be consistent with other dbs */
	destructDB(database);
	stKVDatabaseConf* conf = stKVDatabase_getConf(database);
	const char *basePath = stKVDatabaseConf_getDir(conf);
	visitRecords(basePath, remove_with_arg, NULL);
}

/* check if a record already exists */
static bool containsRecord(stKVDatabase *database, int64_t key)
{
	stSortedSet* sortedSet  = (stSortedSet*)database->dbImpl;
	stIntTuple* tuple = stIntTuple_construct1( key);
	void* found = stSortedSet_search(sortedSet, tuple);
	stIntTuple_destruct(tuple);
	return found != NULL;
}

/* write the record as a file in the directory, and add the key
 * to the in-memory sorted set.
 */
static void insertRecord(stKVDatabase *database, int64_t key, const void *value,
		int64_t sizeOfRecord)
{
	stSortedSet* sortedSet  = (stSortedSet*)database->dbImpl;
	if (stSortedSet_size(sortedSet) >= MAX_NUMBER_ENTRIES)
	{
		stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
				"Database capacity reached: %lld", (int64_t)MAX_NUMBER_ENTRIES);
	}
	char* recordPath = createRecordPath(stKVDatabase_getConf(database), key);
	FILE* recHandle = fopen(recordPath, "wb");
	if (recHandle == NULL)
	{
		stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Open file: %s", recordPath);
	}
	size_t retVal = fwrite(value, sizeOfRecord, 1, recHandle);
	if (retVal != 1)
	{
		stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Write file: %s", recordPath);
	}
	stSortedSet_insert(sortedSet, stIntTuple_construct1( key));
	fclose(recHandle);
	free(recordPath);
}


static void updateRecord(stKVDatabase *database, int64_t key, const void *value,
		int64_t sizeOfRecord)
{
	insertRecord(database, key, value, sizeOfRecord);
}

static void setRecord(stKVDatabase *database, int64_t key, const void *value,
		int64_t sizeOfRecord)
{
	insertRecord(database, key, value, sizeOfRecord);
}

static int64_t numberOfRecords(stKVDatabase *database)
{
	stSortedSet* sortedSet  = (stSortedSet*)database->dbImpl;
	return (int64_t)stSortedSet_size(sortedSet);
}

/*
 * NEEDS TO BE FREED
 */
static void* getRecord2(stKVDatabase *database, int64_t key, int64_t* recordSize)
{
	if (containsRecord(database, key) == false)
	{
		return NULL;
	}
	char* recordPath = createRecordPath(stKVDatabase_getConf(database), key);
	FILE* recHandle = fopen(recordPath, "rb");
	if (recHandle == NULL)
	{
		stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Read file: %s", recordPath);
	}
	fseek(recHandle, 0, SEEK_END);
	long fileSize = ftell(recHandle);
	rewind(recHandle);
	void* buffer = stSafeCMalloc(fileSize);
	size_t retVal = fread(buffer, fileSize, 1, recHandle);
	if (retVal != 1)
	{
		stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Read file: %s", recordPath);
	}
	fclose(recHandle);
	free(recordPath);
	*recordSize = fileSize;
	return buffer;
}

/*
 * NEEDS TO BE FREED
 */
static void* getRecord(stKVDatabase *database, int64_t key)
{
	int64_t i;
	return getRecord2(database, key, &i);
}

/* get part of a string record */
static void *getPartialRecord(stKVDatabase *database, int64_t key,
		int64_t zeroBasedByteOffset, int64_t sizeInBytes, int64_t recordSize)
{
	assert (zeroBasedByteOffset + sizeInBytes <= recordSize);
	if (containsRecord(database, key) == false || sizeInBytes == 0)
	{
		return NULL;
	}
	char* recordPath = createRecordPath(stKVDatabase_getConf(database), key);
	FILE* recHandle = fopen(recordPath, "rb");
	if (recHandle == NULL)
	{
		stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Read file: %s", recordPath);
	}
	fseek(recHandle, (long)zeroBasedByteOffset, 0);
	void* buffer = stSafeCMalloc(sizeInBytes);
	size_t retVal = fread(buffer, sizeInBytes, 1, recHandle);
	if (retVal != 1)
	{
		stThrowNew(ST_KV_DATABASE_EXCEPTION_ID, "Read file: %s", recordPath);
	}
	fclose(recHandle);
	free(recordPath);
	return buffer;
}

static void removeRecord(stKVDatabase *database, int64_t key)
{
	if (containsRecord(database, key) == false)
	{
		stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
				"Removing key not found: %lld", key);
	}
	stSortedSet* sortedSet  = (stSortedSet*)database->dbImpl;
	stIntTuple* tuple = stIntTuple_construct1( key);
	stSortedSet_remove(sortedSet, tuple);
	stIntTuple_destruct(tuple);
	char* recordPath = createRecordPath(stKVDatabase_getConf(database), key);
	int retVal = remove(recordPath);
	if (retVal != 0)
	{
		stThrowNew(ST_KV_DATABASE_EXCEPTION_ID,
				"Removing file returned error %d: %lld ", retVal, key);
	}
	free(recordPath);
}

void stKVDatabase_initialise_bigRecordFile(stKVDatabase *database,
		stKVDatabaseConf *conf, bool create)
{
    database->dbImpl = constructDB(stKVDatabase_getConf(database), create);
    database->secondaryDB = NULL;
    database->destruct = destructDB;
    database->deleteDatabase = deleteDB;
    database->containsRecord = containsRecord;
    database->insertRecord = insertRecord;
    database->insertInt64 = NULL;
    database->updateRecord = updateRecord;
    database->updateInt64 = NULL;
    database->setRecord = setRecord;
    database->incrementInt64 = NULL;
    database->bulkSetRecords = NULL;
    database->bulkRemoveRecords = NULL;
    database->numberOfRecords = numberOfRecords;
    database->getRecord = getRecord;
    database->getInt64 = NULL;
    database->getRecord2 = getRecord2;
    database->getPartialRecord = getPartialRecord;
    database->bulkGetRecords = NULL;
    database->bulkGetRecordsRange = NULL;
    database->removeRecord = removeRecord;
}



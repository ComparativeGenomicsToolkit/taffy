/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef SONLIB__DATABASE_H_
#define SONLIB__DATABASE_H_

#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Cache functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 *
 * Create an empty cache.
 */
stCache *stCache_construct(void);

/*
 * Create an empty cache of limited size. Overhead is not counted
 * toward the size. The cache follows a least-recently-used eviction
 * policy.
 */
stCache *stCache_construct2(size_t maxSize);

/*
 * Destructs the cache.
 */
void stCache_destruct(stCache *cache);

/*
 * Clears the cache.
 */
void stCache_clear(stCache *cache);

/*
 * Update an existing key/value record fragment in the cache. If the record does not exist it is inserted. Throws an exception if unsuccessful.
 * The offset is the start of the record fragment.
 */
void stCache_setRecord(stCache *cache, int64_t key, int64_t zeroBasedByteOffset, int64_t sizeInBytes, const void *value);

/*
 * Returns non-zero if the cache contains all of the given record fragment. If zeroBasedByteOffset=INT64_MAX and
 * sizeInBytes=INT64_MAX then no overlap is required.
 */
bool stCache_containsRecord(stCache *cache, int64_t key, int64_t zeroBasedByteOffset, int64_t sizeInBytes);

/*
 * Gets a record from the cache, given the key. This function allows the partial retrieval of a record, using the
 * given offset and the size of the requested retrieval.
 * The function returns NULL if the desired part of the record is not in the cache exist.
 * If the sizeInBytes equals INT64_MAX then the maximal size contiguous fragment of the record is returned.
 * Record size is initialised with the size of the returned buffer.
 */
void *stCache_getRecord(stCache *cache, int64_t key, int64_t zeroBasedByteOffset, int64_t sizeInBytes, int64_t *recordSize);

/*
 * Returns non-zero iff the two buffers have the same size and are identical.
 */
bool stCache_recordsIdentical(const char *value, int64_t sizeOfRecord,
        const char *updatedValue, int64_t updatedSizeOfRecord);

/*
 * Gets the current size of the cache.
 */
size_t stCache_size(stCache *cache);

#ifdef __cplusplus
}
#endif
#endif

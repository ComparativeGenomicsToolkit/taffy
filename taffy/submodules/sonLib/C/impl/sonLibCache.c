/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

//Cache functions

#include "sonLibGlobalsInternal.h"

struct _lru;

typedef struct _cacheRecord {
    /*
     * A little object for storing records in the cache.
     */
    int64_t key, start, size;
    char *record;
    struct _lru *lru;
} stCacheRecord;

// Doubly-linked list containing the LRU entries.
typedef struct _lru {
    stCacheRecord *record;
    struct _lru *prev;
    struct _lru *next;
} lru;

struct stCache {
    stSortedSet *cache;
    // Head of the LRU list (i.e., *most* recently used record).
    lru *lruHead;
    // Tail of the LRU list (i.e., *least* recently used record).
    lru *lruTail;
    // Current size of records (not including overhead).
    size_t curSize;
    // Maximum allowed size of all records (not including overhead).
    // The cache may go over this size if a single record greater than
    // maxSize is inserted.
    size_t maxSize;
};

// Free a cache record without touching the LRU list.
static void cacheRecord_destruct(stCacheRecord *i) {
    free(i->record);
    free(i);
}

static void lru_destruct(lru *lru) {
    free(lru);
}

static lru *lru_construct(stCacheRecord *record) {
    lru *ret = st_calloc(1, sizeof(lru));
    ret->record = record;
    return ret;
}

// Empty the LRU list in preparation for the cache to be destroyed.
static void clearLruList(stCache *cache) {
    while (cache->lruHead != NULL) {
        if (cache->lruTail == cache->lruHead) {
            cache->lruTail = NULL;
        }
        lru *next = cache->lruHead->next;
        lru_destruct(cache->lruHead);
        cache->lruHead = next;
    }
    assert(cache->lruTail == NULL);
}

// Remove a single LRU entry from the LRU list, sealing the gap it
// leaves behind.
static void removeLruEntry(stCache *cache, lru *entry) {
    if (cache->lruHead == entry) {
        cache->lruHead = entry->next;
    }
    if (cache->lruTail == entry) {
        cache->lruTail = entry->prev;
    }
    lru *prev = entry->prev;
    lru *next = entry->next;
    if (prev != NULL) {
        prev->next = next;
    }
    if (next != NULL) {
        next->prev = prev;
    }
    entry->next = NULL;
    entry->prev = NULL;
}

// This entry just got used, move it to the front.
static void moveLruToFront(stCache *cache, lru *entry) {
    if (cache->lruHead == entry) {
        return;
    }
    lru *prevHead = cache->lruHead;
    removeLruEntry(cache, entry);
    entry->next = prevHead;
    if (prevHead != NULL) {
        prevHead->prev = entry;
    }
    cache->lruHead = entry;
    if (cache->lruTail == NULL) {
        cache->lruTail = entry;
    }
}

// Remove and free a cache record properly, updating the LRU and used space.
static void removeRecordFromCache(stCache *cache, stCacheRecord *i) {
    stSortedSet_remove(cache->cache, i);
    assert(cache->curSize >= i->size);
    cache->curSize -= i->size;
    lru *lru = i->lru;
    removeLruEntry(cache, lru);
    lru_destruct(lru);
    cacheRecord_destruct(i);
}

// Free enough space in the cache to fit a new record with the given size.
// Least-recently-used entries will be deleted first.
static void freeSpaceInCache(stCache *cache, size_t size) {
    while (cache->lruTail != NULL && (cache->maxSize - cache->curSize) < size) {
        assert(cache->lruTail->next == NULL);
        lru *prev = cache->lruTail->prev;
        removeRecordFromCache(cache, cache->lruTail->record);
        cache->lruTail = prev;
    }
}

static int cacheRecord_cmp(const void *a, const void *b) {
    const stCacheRecord *i = a;
    const stCacheRecord *j = b;
    if (i->key > j->key) {
        return 1;
    }
    if (i->key < j->key) {
        return -1;
    }
    if (i->start > j->start) {
        return 1;
    }
    if (i->start < j->start) {
        return -1;
    }
    return 0;
}

static stCacheRecord getTempRecord(int64_t key, int64_t start, int64_t size) {
    stCacheRecord record;
    record.key = key;
    record.start = start;
    record.size = size;
    record.record = NULL;
    record.lru = NULL;
    return record;
}

static stCacheRecord *cacheRecord_construct(stCache *cache, int64_t key, const void *value,
                                            int64_t start, int64_t size, bool copyMemory) {
    assert(value != NULL);
    assert(size >= 0);
    assert(start >= 0);
    stCacheRecord *record = st_malloc(sizeof(stCacheRecord));
    record->key = key;
    record->start = start;
    record->size = size;
    record->record = copyMemory ? memcpy(st_malloc(size), value, size) : (char *) value;
    record->lru = lru_construct(record);
    // Mark as most recently used.
    moveLruToFront(cache, record->lru);
    cache->curSize += size;
    return record;
}

static stCacheRecord *getLessThanOrEqualRecord(stCache *cache,
                                               int64_t key, int64_t start, int64_t size) {
    stCacheRecord record = getTempRecord(key, start, size);
    return stSortedSet_searchLessThanOrEqual(cache->cache, &record);
}

static stCacheRecord *getGreaterThanOrEqualRecord(stCache *cache,
                                                  int64_t key, int64_t start, int64_t size) {
    stCacheRecord record = getTempRecord(key, start, size);
    return stSortedSet_searchGreaterThanOrEqual(cache->cache, &record);
}

static bool recordContainedIn(stCacheRecord *record, int64_t key,
                              int64_t start, int64_t size) {
    /*
     * Returns non zero if the record is contained in the given range.
     */
    return record->key == key && record->start >= start && record->start
        + record->size <= start + size;
}

static bool recordOverlapsWith(stCacheRecord *record, int64_t key,
                               int64_t start, int64_t size) {
    /*
     * Returns non zero if the record overlaps with the given range.
     */
    if (record->key != key) {
        return 0;
    }
    if (record->start >= start) {
        return record->start < start + size;
    }
    return record->start + record->size > start;
}

static bool recordsAdjacent(stCacheRecord *record1, stCacheRecord *record2) {
    /*
     * Returns non-zero if the records abut with record1 immediately before record2.
     */
    assert(cacheRecord_cmp(record1, record2) <= 0);
    return record1->key == record2->key && record1->start + record1->size
        == record2->start;
}

static stCacheRecord *mergeRecords(stCache *cache,
                                   stCacheRecord *record1,
                                   stCacheRecord *record2) {
    /*
     * Merges two adjacenct records.
     */
    assert(recordsAdjacent(record1, record2));
    int64_t i = record1->size + record2->size;
    char *j = memcpy(st_malloc(i), record1->record, record1->size);
    memcpy(j + record1->size, record2->record, record2->size);
    stCacheRecord *record3 = cacheRecord_construct(cache, record1->key, j,
                                                   record1->start, i, 0);

    removeRecordFromCache(cache, record1);
    removeRecordFromCache(cache, record2);
    stSortedSet_insert(cache->cache, record3);
    return record3;
}

void deleteRecord(stCache *cache, int64_t key,
                  int64_t start, int64_t size) {
    assert(!stCache_containsRecord(cache, key, start, size)); //Will not delete a record wholly contained in.
    stCacheRecord *record = getLessThanOrEqualRecord(cache, key, start,
                                                     size);
    while (record != NULL && recordOverlapsWith(record, key, start, size)) { //could have multiple fragments in there to remove.
        if (recordContainedIn(record, key, start, size)) { //We get rid of the record because it is contained in the range
            removeRecordFromCache(cache, record);
            record = getLessThanOrEqualRecord(cache, key, start, size);
        } else { //The range overlaps with, but is not fully contained in, so we trim it..
            assert(record->start < start);
            assert(record->start + record->size > start);
            size_t newSize = start - record->start;
            cache->curSize -= record->size - newSize;
            record->size = newSize;
            assert(record->size >= 0);
            break;
        }
    }
    record = getGreaterThanOrEqualRecord(cache, key, start, size);
    while (record != NULL && recordOverlapsWith(record, key, start, size)) { //could have multiple fragments in there to remove.
        if (recordContainedIn(record, key, start, size)) { //We get rid of the record because it is contained in the range
            removeRecordFromCache(cache, record);
            record = getGreaterThanOrEqualRecord(cache, key, start, size);
        } else { //The range overlaps with, but is not fully contained in, so we trim it..
            assert(record->start < start + size);
            assert(record->start > start);
            int64_t newSize = record->size - (start + size - record->start);
            cache->curSize -= record->size - newSize;
            int64_t newStart = start + size;
            assert(newSize >= 0);
            char *newMem = memcpy(st_malloc(newSize),
                                  record->record + start + size - record->start, newSize);
            free(record->record);
            record->record = newMem;
            record->start = newStart;
            record->size = newSize;
            break; //We can break at this point as we have reached the end of the range (as the record overlapped)
        }
    }
}


/*
 * Public functions
 */

stCache *stCache_construct2(size_t maxSize) {
    stCache *cache = st_malloc(sizeof(stCache));
    cache->cache = stSortedSet_construct3(cacheRecord_cmp,
                                          (void(*)(void *)) cacheRecord_destruct);
    cache->lruHead = NULL;
    cache->lruTail = NULL;
    cache->curSize = 0;
    cache->maxSize = maxSize;
    return cache;
}

stCache *stCache_construct(void) {
    return stCache_construct2(SIZE_MAX);
}

void stCache_destruct(stCache *cache) {
    stSortedSet_destruct(cache->cache);
    clearLruList(cache);
    free(cache);
}

void stCache_clear(stCache *cache) {
    stSortedSet_destruct(cache->cache);
    clearLruList(cache);
    cache->cache = stSortedSet_construct3(cacheRecord_cmp,
                                          (void(*)(void *)) cacheRecord_destruct);
    cache->curSize = 0;
}

void stCache_setRecord(stCache *cache, int64_t key,
                       int64_t start, int64_t size, const void *value) {
    //If the record is already contained we update a portion of it.
    assert(value != NULL);
    if (stCache_containsRecord(cache, key, start, size)) {
        stCacheRecord *record = getLessThanOrEqualRecord(cache, key,
                                                         start, size);
        assert(record != NULL);
        assert(record->key == key);
        assert(record->start <= start);
        assert(record->start + record->size >= start + size);
        memcpy(record->record + start - record->start, value, size);
        // Mark this entry as the most recently used.
        moveLruToFront(cache, record->lru);
        return;
    }

    //Get rid of scattered bits that are contained in this record..
    deleteRecord(cache, key, start, size);

    // Free up enough space from the least-recently-used elements to
    // fit the new record.
    freeSpaceInCache(cache, size);

    //Now get any left and right bits
    stCacheRecord *record1 = getLessThanOrEqualRecord(cache, key, start,
                                                      size);
    stCacheRecord *record2 = cacheRecord_construct(cache, key, value, start,
                                                   size, 1);
    stCacheRecord *record3 = getGreaterThanOrEqualRecord(cache, key,
                                                         start, size);
    assert(record2 != NULL);
    stSortedSet_insert(cache->cache, record2);
    if (record1 != NULL && recordsAdjacent(record1, record2)) {
        record2 = mergeRecords(cache, record1, record2);
    }
    if (record3 != NULL && recordsAdjacent(record2, record3)) {
        record2 = mergeRecords(cache, record2, record3);
    }
}

bool stCache_containsRecord(stCache *cache, int64_t key,
                            int64_t start, int64_t size) {
    assert(start >= 0);
    assert(size >= 0);
    stCacheRecord *record = getLessThanOrEqualRecord(cache, key, start,
                                                     size);
    if (record == NULL) {
        return 0;
    }
    if (record->key != key) {
        return 0;
    }
    assert(record->start <= start);
    if (size != INT64_MAX && start + size > record->start + record->size) { //If the record has a known length check we have all that we want.
        return 0;
    }
    return 1;
}

void *stCache_getRecord(stCache *cache, int64_t key,
                        int64_t start, int64_t size, int64_t *sizeRead) {
    if (stCache_containsRecord(cache, key, start, size)) {
        stCacheRecord *record = getLessThanOrEqualRecord(cache, key,
                                                         start, size);
        assert(record != NULL);
        // Mark this entry as the most recently used.
        moveLruToFront(cache, record->lru);
        int64_t j = start - record->start;
        int64_t i = size == INT64_MAX ? (record->size - j) : size;
        assert(record->start <= start);
        assert(j >= 0);
        assert(j <= record->size);
        assert(i >= 0);
        assert(j + i <= record->size);
        *sizeRead = i;
        char *cA = st_malloc(i);
        memcpy(cA, record->record + j, i);
        return cA;
    }
    return NULL;
}

bool stCache_recordsIdentical(const char *value, int64_t sizeOfRecord,
                              const char *updatedValue, int64_t updatedSizeOfRecord) {
    if (sizeOfRecord != updatedSizeOfRecord) {
        return 0;
    }
    for (int64_t j = 0; j < sizeOfRecord; j++) {
        if (value[j] != updatedValue[j]) {
            return 0;
        }
    }
    return 1;
}

size_t stCache_size(stCache *cache) {
    return cache->curSize;
}

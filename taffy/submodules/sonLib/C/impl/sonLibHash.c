/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibHash.c
 *
 *  Created on: 4 Apr 2010
 *      Author: benedictpaten
 */
#include "sonLibGlobalsInternal.h"
#include "hashTableC.h"
#include "hashTableC_itr.h"

struct _stHash {
    struct hashtable *hash;
    bool destructKeys, destructValues;
};

uint64_t stHash_pointer(const void *k) {
    // Size doesn't matter; just promote to 64 bits
    uint64_t key = (uint64_t) k;
    
    // Use the hash from <https://stackoverflow.com/a/12996028>
    // We can't just -ull these until C++11, if we're in C++
    key = (key ^ (key >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    key = (key ^ (key >> 27)) * UINT64_C(0x94d049bb133111eb);
    return key ^ (key >> 31);
}

static int stHash_equalKey(const void *key1, const void *key2) {
    return key1 == key2;
}

stHash *stHash_construct(void) {
    return stHash_construct3(stHash_pointer, stHash_equalKey, NULL, NULL);
}

stHash *stHash_construct2(void(*destructKeys)(void *), void(*destructValues)(void *)) {
    return stHash_construct3(stHash_pointer, stHash_equalKey, destructKeys, destructValues);
}

stHash *stHash_construct3(uint64_t(*hashKey)(const void *), int(*hashEqualsKey)(const void *, const void *), void(*destructKeys)(void *), void(*destructValues)(void *)) {
    stHash *hash = st_malloc(sizeof(stHash));
    hash->hash = create_hashtable(0, hashKey, hashEqualsKey, destructKeys, destructValues);
    hash->destructKeys = destructKeys != NULL;
    hash->destructValues = destructValues != NULL;
    return hash;
}

void stHash_destruct(stHash *hash) {
    hashtable_destroy(hash->hash, hash->destructValues, hash->destructKeys);
    free(hash);
}

void stHash_setDestructKeys(stHash *hash, void(*destructor)(void *)) {
    hash->destructKeys = destructor != NULL;
    hash->hash->keyFree = destructor;
}

void stHash_setDestructValues(stHash *hash, void(*destructor)(void *)) {
    hash->destructValues = destructor != NULL;
    hash->hash->valueFree = destructor;
}

void stHash_insert(stHash *hash, void *key, void *value) {
    if (stHash_search(hash, key) != NULL) { //This will ensure we don't end up with duplicate keys..
        stHash_remove(hash, key);
    }
    hashtable_insert(hash->hash, key, value);
}

void *stHash_search(stHash *hash, void *key) {
    return hashtable_search(hash->hash, key);
}

void *stHash_remove(stHash *hash, void *key) {
    return hashtable_remove(hash->hash, key, 0);
}

void *stHash_removeAndFreeKey(stHash *hash, void *key) {
    return hashtable_remove(hash->hash, key, 1);
}

int64_t stHash_size(stHash *hash) {
    return hashtable_count(hash->hash);
}

stHashIterator *stHash_getIterator(stHash *hash) {
    return hashtable_iterator(hash->hash);
}

void *stHash_getNext(stHashIterator *iterator) {
    if (iterator->e != NULL) {
        void *o = hashtable_iterator_key(iterator);
        hashtable_iterator_advance(iterator);
        return o;
    }
    return NULL;
}

stHashIterator *stHash_copyIterator(stHashIterator *iterator) {
    stHashIterator *iterator2 = st_malloc(sizeof(stHashIterator));
    iterator2->h = iterator->h;
    iterator2->e = iterator->e;
    iterator2->parent = iterator->parent;
    iterator2->index = iterator->index;
    return iterator2;
}

void stHash_destructIterator(stHashIterator *iterator) {
    free(iterator);
}

stList *stHash_getKeys(stHash *hash) {
    stList *list = stList_construct();
    stHashIterator *iterator = stHash_getIterator(hash);
    void *item;
    while ((item = stHash_getNext(iterator)) != NULL) {
        stList_append(list, item);
    }
    stHash_destructIterator(iterator);
    return list;
}

stList *stHash_getValues(stHash *hash) {
    stList *list = stList_construct();
    stHashIterator *iterator = stHash_getIterator(hash);
    void *item;
    while ((item = stHash_getNext(iterator)) != NULL) {
        stList_append(list, stHash_search(hash, item));
    }
    stHash_destructIterator(iterator);
    return list;
}

/*
 * Useful hash keys/equals functions
 */

uint64_t stHash_stringKey(const void *k) {
    // djb2
    // This algorithm was first reported by Dan Bernstein
    // many years ago in comp.lang.c
    //
    uint64_t hash = 0; //5381;
    int c;
    char *cA;
    cA = (char *) k;
    while ((c = *cA++) != '\0') {
        hash = c + (hash << 6) + (hash << 16) - hash;
        //hash = ((hash << 5) + hash) + c; // hash*33 + c
    }
    return hash;
}

int stHash_stringEqualKey(const void *key1, const void *key2) {
    return strcmp(key1, key2) == 0;
}

stHash *stHash_invert(stHash *hash, uint64_t(*hashKey)(const void *), int(*equalsFn)(const void *, const void *),
        void(*destructKeys)(void *), void(*destructValues)(void *)) {
    /*
     * Inverts the hash.
     */
    stHash *invertedHash = stHash_construct3(hashKey, equalsFn, destructKeys, destructValues);
    stHashIterator *hashIt = stHash_getIterator(hash);
    void *key;
    while ((key = stHash_getNext(hashIt)) != NULL) {
        void *value = stHash_search(hash, key);
        assert(value != NULL);
        if (stHash_search(invertedHash, value) == NULL) {
            stHash_insert(invertedHash, value, key);
        }
    }
    stHash_destructIterator(hashIt);
    return invertedHash;
}
// interface to underlying functions
uint64_t (*stHash_getHashFunction(stHash *hash))(const void *) {
    return hash->hash->hashfn;
}
int (*stHash_getEqualityFunction(stHash *hash))(const void *, const void *) {
    return hash->hash->eqfn;
}
void (*stHash_getKeyDestructorFunction(stHash *hash))(void *) {
    return hash->hash->keyFree;
}
void (*stHash_getValueDestructorFunction(stHash *hash))(void *) {
    return hash->hash->valueFree;
}

static int unsigned_cmp(const unsigned *x, const unsigned *y) {
    if (*x < *y) {
        return -1;
    } else if (*x == *y) {
        return 0;
    } else {
        return 1;
    }
}

void stHash_printDiagnostics(stHash *hash) {
    struct hashtable *h = hash->hash;
    unsigned *bucketLoad = st_malloc(h->tablelength * sizeof(unsigned));
    unsigned *occupiedBucketLoad = st_malloc(h->entrycount * sizeof(unsigned));
    size_t occupiedBuckets = 0;
    for (size_t i = 0; i < h->tablelength; i++) {
        unsigned load = 0;
        struct entry *e;
        if ((e = h->table[i]) != NULL) {
            while (e) {
                load++;
                e = e->next;
            }
            occupiedBucketLoad[occupiedBuckets++] = load;
        }
        bucketLoad[i] = load;
    }
    printf("Load: %" PRIi64 " / %" PRIi64 " (%lf%%)\n", h->entrycount, h->tablelength, ((double)h->entrycount)/h->tablelength * 100);
    printf("# occupied buckets: %zu\n", occupiedBuckets);
    qsort(bucketLoad, h->tablelength, sizeof(unsigned), (int (*)(const void *, const void *)) unsigned_cmp);
    qsort(occupiedBucketLoad, occupiedBuckets, sizeof(unsigned), (int (*)(const void *, const void *)) unsigned_cmp);
    printf("min load: %u, median load: %u, max load: %u\n", bucketLoad[0], bucketLoad[(h->tablelength - 1) / 2],
           bucketLoad[h->tablelength - 1]);
    printf("min occupied load: %u, avg occupied load: %lf, median occupied load: %u, max occupied load: %u\n",
           occupiedBucketLoad[0], ((double)h->entrycount)/occupiedBuckets,
           occupiedBucketLoad[(occupiedBuckets - 1) / 2],
           occupiedBucketLoad[occupiedBuckets - 1]);
    free(occupiedBucketLoad);
    free(bucketLoad);
}

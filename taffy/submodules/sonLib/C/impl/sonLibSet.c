/*
 * Copyright (C) 2012 by Dent Earl dentearl (a) gmail com
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibSet.c
 *
 *  Created on: 12 July 2012
 *      Author: dentearl
 */
#include "sonLibGlobalsInternal.h"
#include "sonLibHash.h"
#include "sonLibSet.h"

const char *SET_EXCEPTION_ID = "SET_EXCEPTION";

struct _stSet {
    stHash *hash;
};
struct _stSetIterator {
    stHashIterator *hashIterator;
};

uint64_t stSet_pointer(const void *k) {
    return (uint64_t) (size_t) k; // Use the low order bits
}
static int stSet_equalKey(const void *key1, const void *key2) {
    return key1 == key2;
}
stSet *stSet_construct(void) {
    return stSet_construct3(stSet_pointer, stSet_equalKey, NULL);
}
stSet *stSet_construct2(void(*destructKeys)(void *)) {
    return stSet_construct3(stSet_pointer, stSet_equalKey, destructKeys);
}
stSet *stSet_construct3(uint64_t(*hashKey)(const void *), int(*hashEqualsKey)(const void *, const void *), void(*destructKeys)(void *)) {
    stSet *set = st_malloc(sizeof(*set));
    set->hash = stHash_construct3(hashKey, hashEqualsKey, destructKeys, NULL);
    return set;
}
void stSet_destruct(stSet *set) {
    stHash_destruct(set->hash);
    free(set);
}

void stSet_setDestructor(stSet *set, void(*destructor)(void *)) {
    stHash_setDestructKeys(set->hash, destructor);
}

void stSet_insert(stSet *set, void *key) {
    if (stSet_search(set, key) != NULL) { // This will ensure we don't end up with duplicate keys..
        stSet_remove(set, key);
    }
    stHash_insert(set->hash, key, key);
}
void stSet_insertAll(stSet *set, stSet *setToAdd) {
    stSetIterator *setIt = stSet_getIterator(setToAdd);
    void *o;
    while((o = stSet_getNext(setIt)) != NULL) {
        stSet_insert(set, o);
    }
    stSet_destructIterator(setIt);
}
void *stSet_search(stSet *set, void *key) {
    return stHash_search(set->hash, key);
}
void *stSet_remove(stSet *set, void *key) {
    return stHash_remove(set->hash, key);
}
void *stSet_removeAndFreeKey(stSet *set, void *key) {
    return stHash_removeAndFreeKey(set->hash, key);
}

void stSet_removeAll(stSet *set, stSet *subset) {
    stSetIterator *it = stSet_getIterator(subset);

    void *o;
    while((o = stSet_getNext(it)) != NULL) {
        stSet_remove(set, o);
    }
    stSet_destructIterator(it);
}

int64_t stSet_size(stSet *set) {
    return stHash_size(set->hash);
}
stSetIterator *stSet_getIterator(stSet *set) {
    stSetIterator *sit = st_malloc(sizeof(*sit));
    sit->hashIterator = stHash_getIterator(set->hash);
    return sit;
}
void *stSet_getNext(stSetIterator *iterator) {
    return stHash_getNext(iterator->hashIterator);
}
stSetIterator *stSet_copyIterator(stSetIterator *iterator) {
    stSetIterator *iterator2 = st_malloc(sizeof(*iterator2));
    iterator2->hashIterator = stHash_copyIterator(iterator->hashIterator);
    return iterator2;
}
void stSet_destructIterator(stSetIterator *iterator) {
    stHash_destructIterator(iterator->hashIterator);
    free(iterator);
}
stList *stSet_getKeys(stSet *set) {
    stList *list = stList_construct();
    stSetIterator *iterator = stSet_getIterator(set);
    void *item;
    while ((item = stSet_getNext(iterator)) != NULL) {
        stList_append(list, item);
    }
    stSet_destructIterator(iterator);
    return list;
}
stList *stSet_getList(stSet *set) {
    return stSet_getKeys(set);
}
uint64_t (*stSet_getHashFunction(stSet *set))(const void *) {
    return stHash_getHashFunction(set->hash);
}
int (*stSet_getEqualityFunction(stSet *set))(const void *, const void *) {
    return stHash_getEqualityFunction(set->hash);
}
void (*stSet_getDestructorFunction(stSet *set))(void *) {
    return stHash_getKeyDestructorFunction(set->hash);
}
static int stSet_hashersEqual(stSet *set1, stSet *set2) {
    return stSet_getHashFunction(set1) == stSet_getHashFunction(set2);
}
static int stSet_equalitiesEqual(stSet *set1, stSet *set2) {
    return stSet_getEqualityFunction(set1) == stSet_getEqualityFunction(set2);
}
void *stSet_peek(stSet *set) {
    if(stSet_size(set) == 0) {
        stThrowNew(SET_EXCEPTION_ID, "Set is empty.");
    }
    stSetIterator *setIt = stSet_getIterator(set);
    void *object = stSet_getNext(setIt);
    stSet_destructIterator(setIt);
    return object;
}
static void stSet_verifySetsHaveSameFunctions(stSet *set1, stSet *set2) {
    if (!stSet_hashersEqual(set1, set2)) {
        stThrowNew(SET_EXCEPTION_ID, "Comparator functions are not equal for "
                   "two sets.");
    }
    if (!stSet_equalitiesEqual(set1, set2)) {
        stThrowNew(SET_EXCEPTION_ID, "Hash functions are not equal for "
                   "two sets.");
    }
}
stSet *stSet_getUnion(stSet *set1, stSet *set2) {
    stSet_verifySetsHaveSameFunctions(set1, set2);
    stSet *set3 = stSet_construct3(stSet_getHashFunction(set1), 
                                   stSet_getEqualityFunction(set1),
                                   NULL);
    // Add everything
    stSetIterator *sit= stSet_getIterator(set1);
    void *o;
    while ((o = stSet_getNext(sit)) != NULL) {
        stSet_insert(set3, o);
    }
    stSet_destructIterator(sit);
    sit = stSet_getIterator(set2);
    while ((o = stSet_getNext(sit)) != NULL) {
        stSet_insert(set3, o);
    }
    stSet_destructIterator(sit);
    return set3;
}

int64_t stSet_sizeOfIntersection(stSet *set1, stSet *set2) {
    stSet *s = stSet_getIntersection(set1, set2);
    int64_t i = stSet_size(s);
    stSet_destruct(s);
    return i;
}

stSet *stSet_getIntersection(stSet *set1, stSet *set2) {
    stSet_verifySetsHaveSameFunctions(set1, set2);
    stSet *set3 = stSet_construct3(stSet_getHashFunction(set1), 
                                   stSet_getEqualityFunction(set1),
                                   NULL);
    // Add those from set1 only if they are also in set2
    stSetIterator *sit= stSet_getIterator(set1);
    void *o;
    while ((o = stSet_getNext(sit)) != NULL) {
        if (stSet_search(set2, o) != NULL) {
            stSet_insert(set3, o);
        }
    }
    stSet_destructIterator(sit);
    return set3;
}
stSet *stSet_getDifference(stSet *set1, stSet *set2) {
    stSet_verifySetsHaveSameFunctions(set1, set2);
    stSet *set3 = stSet_construct3(stSet_getHashFunction(set1), 
                                   stSet_getEqualityFunction(set1),
                                   NULL);
    // Add those from set1 only if they are not in set2
    stSetIterator *sit= stSet_getIterator(set1);
    void *o;
    while ((o = stSet_getNext(sit)) != NULL) {
        if (stSet_search(set2, o) == NULL) {
            stSet_insert(set3, o);
        }
    }
    stSet_destructIterator(sit);
    return set3;
}

bool stSet_equals(stSet *set1, stSet *set2) {
    if(stSet_size(set1) != stSet_size(set2)) {
        return 0;
    }
    stSet *setI = stSet_getIntersection(set1, set2);
    bool b = stSet_size(setI) == stSet_size(set1);
    stSet_destruct(setI);
    return b;
}

bool stSet_isSubset(stSet *parentSet, stSet *putativeSubset) {
    stSet *setI = stSet_getIntersection(parentSet, putativeSubset);
    bool b = stSet_size(setI) == stSet_size(putativeSubset);
    stSet_destruct(setI);
    return b;
}

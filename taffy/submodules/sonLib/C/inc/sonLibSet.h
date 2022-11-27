/*
 * Copyright (C) 2012 by Dent Earl dentearl (a) gmail com
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef SONLIB_SET_H_
#define SONLIB_SET_H_

/*
 * sonLibSet.h
 *
 *  Created on: 12 July 2012
 *      Author: dentearl
 */

#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

// exception string
extern const char *SET_EXCEPTION_ID;

/*
 * Function which generates hash key from pointer, should work well regardless of pointer size.
 */
uint64_t stSet_pointer(const void *k);

/*
 * Constructs set, with no destructors for keys or values.
 */
stSet *stSet_construct(void);

/*
 * Constructs a set with given destructor, if null then destructor is ignored
 */
stSet *stSet_construct2(void (*destructKeys)(void *));

/*
 * Constructs a set using the given comparison functions.
 */
stSet *stSet_construct3(uint64_t (*hashKey)(const void *), int (*hashEqualsKey)(const void *, const void *),
                        void (*destructKeys)(void *));

/*
 * Destructs a set.
 */
void stSet_destruct(stSet *set);

/*
 * Set the destructor.
 */
void stSet_setDestructor(stSet *set, void(*destructor)(void *));

/*
 * Insert element, overiding if already present.
 */
void stSet_insert(stSet *set, void *key);

/*
 * Insert all elements in setToAdd to set.
 */
void stSet_insertAll(stSet *set, stSet *setToAdd);

/*
 * Search for value, returns null if not present.
 */
void *stSet_search(stSet *set, void *key);

/*
 * Removes element, returning removed element.
 */
void *stSet_remove(stSet *set, void *key);

/*
 * Removes the given subset from the given set.
 */
void stSet_removeAll(stSet *set, stSet *subset);

/*
 * Removes element, returning removed element and freeing key (using supplied function).
 */
void *stSet_removeAndFreeKey(stSet *set, void *key);

/*
 * Returns the number of keys in the set.
 */
int64_t stSet_size(stSet *set);

/*
 * Returns an arbitrary member of the set; if empty throws a set exception.
 */
void *stSet_peek(stSet *set);

/*
 * Returns an iterator of the keys in the set.
 */
stSetIterator *stSet_getIterator(stSet *set);

/*
 * Gets the next key from the iterator.
 */
void *stSet_getNext(stSetIterator *iterator);

/*
 * Duplicates the iterator.
 */
stSetIterator *stSet_copyIterator(stSetIterator *iterator);

/*
 * Destructs the iterator.
 */
void stSet_destructIterator(stSetIterator *iterator);

/*
 * Gets the keys in the set as list.
 */
stList *stSet_getKeys(stSet *set);
stList *stSet_getList(stSet *set);

// Set Functions
stSet *stSet_getUnion(stSet *set1, stSet *set2);
stSet *stSet_getIntersection(stSet *set1, stSet *set2);
stSet *stSet_getDifference(stSet *set1, stSet *set2);
int64_t stSet_sizeOfIntersection(stSet *set1, stSet *set2);

// Comparison functions
/*
 * Returns non-zero iff set1 and set2 have equal contents.
 */
bool stSet_equals(stSet *set1, stSet *set2);

/*
 * Returns non-zero iff the putativeSubset is a subset of parentSet.
 */
bool stSet_isSubset(stSet *parentSet, stSet *putativeSubset);

// Get access to the underlying functions
uint64_t (*stSet_getHashFunction(stSet *set))(const void *);
int (*stSet_getEqualityFunction(stSet *set))(const void *, const void *);
void (*stSet_getDestructorFunction(stSet *set))(void *);

#ifdef __cplusplus
}
#endif // __cplusplus
#endif // SONLIB_SET_H_

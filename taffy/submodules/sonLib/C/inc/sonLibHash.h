/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef SONLIB_HASH_H_
#define SONLIB_HASH_H_

/*
 * sonLibHash.h
 *
 *  Created on: 4 Apr 2010
 *      Author: benedictpaten
 */

#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

// FIXME: passing key as non-const is causing unnecessary casts

/*
 * Function which generates hash key from pointer, should work well regardless of pointer size.
 */
uint64_t stHash_pointer( const void *k );

/*
 * Constructs hash, with no destructors for keys or values.
 */
stHash *stHash_construct(void);

/*
 * Constructs a hash with given destructors, if null then destructors are ignored
 */
stHash *stHash_construct2(void (*destructKeys)(void *), void (*destructValues)(void *));

/*
 * Constructs a hash using the given comparison functions.
 */
stHash *stHash_construct3(uint64_t (*hashKey)(const void *), int (*hashEqualsKey)(const void *, const void *),
                          void (*destructKeys)(void *), void (*destructValues)(void *));

/*
 * Destructs a hash.
 */
void stHash_destruct(stHash *hash);

/*
 * Set destructor for the keys
 */
void stHash_setDestructKeys(stHash *hash, void(*destructor)(void *));

/*
 * Set destructor for the values
 */
void stHash_setDestructValues(stHash *hash, void(*destructor)(void *));

/*
 * Insert element, overiding if already present.
 */
void stHash_insert(stHash *hash, void *key, void *value);

/*
 * Search for value, returns null if not present.
 */
void *stHash_search(stHash *hash, void *key);

/*
 * Removes element, returning removed element.
 */
void *stHash_remove(stHash *hash, void *key);

/*
 * Removes element, returning removed element and freeing key (using supplied function).
 */
void *stHash_removeAndFreeKey(stHash *hash, void *key);

/*
 * Returns the number of key/value pairs in the hash.
 */
int64_t stHash_size(stHash *hash);

/*
 * Returns an iterator of the keys in the hash.
 */
stHashIterator *stHash_getIterator(stHash *hash);

/*
 * Gets the next key from the iterator.
 */
void *stHash_getNext(stHashIterator *iterator);

/*
 * Duplicates the iterator.
 */
stHashIterator *stHash_copyIterator(stHashIterator *iterator);

/*
 * Destructs the iterator.
 */
void stHash_destructIterator(stHashIterator *iterator);

/*
 * Gets the keys in the hash as list.
 */
stList *stHash_getKeys(stHash *hash);

/*
 * Gets the values in the hash as a list.
 */
stList *stHash_getValues(stHash *hash);

/*
 * Useful hash keys..
 */

/*
 * A hash function for a char string.
 */
uint64_t stHash_stringKey( const void *k );

/*
 * A hash equals function for two char strings.
 */
int stHash_stringEqualKey( const void *key1, const  void *key2 );

/*
 * Invert the hash, such that the values become the keys and vice versa. Where there are multiple keys
 * mapping to the same value, only the first encountered key is stored.
 */
stHash *stHash_invert(stHash *hash, uint64_t (*hashKey)(const void *),
        int(*equalsFn)(const void *, const void *), void (*destructKeys)(void *), void (*destructValues)(void *));

// access to the underlying functions:
uint64_t (*stHash_getHashFunction(stHash *hash))(const void *);
int (*stHash_getEqualityFunction(stHash *hash))(const void *, const void *);
void (*stHash_getKeyDestructorFunction(stHash *hash))(void *);
void (*stHash_getValueDestructorFunction(stHash *hash))(void *);

/*
 * Print a set of statistics (load, occupied load, etc.) for the
 * hash. Meant to be run from gdb to check if your hash function is
 * distributing properly.
 */
void stHash_printDiagnostics(stHash *hash);

#ifdef __cplusplus
}
#endif
#endif

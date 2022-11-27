/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef SONLIB_SORTED_SET_H_
#define SONLIB_SORTED_SET_H_

#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

//The exception string
extern const char *SORTED_SET_EXCEPTION_ID;

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Sorted set functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a sorted set, using pointer based comparison function.
 */
stSortedSet *stSortedSet_construct(void);

/*
 * Constructs a sorted set, using pointer based comparison function,
 * and the given destruct element function, which will be run on each element when
 * the set is destructed.
 */
stSortedSet *stSortedSet_construct2(void (*destructElementFn)(void *));

/*
 * Constructs a sorted set, using the given comparison function and destruct element function.
 * If destruct element function is null then it is ignored.
 */
stSortedSet *stSortedSet_construct3(int (*compareFn)(const void *, const void *),
                                      void (*destructElementFn)(void *));

/*
 * Clones the given sorted set, setting the element destructor to the given function.
 */
stSortedSet *stSortedSet_copyConstruct(stSortedSet *sortedSet, void (*destructElementFn)(void *));

/*
 * Set the destructor for the set.
 */
void stSortedSet_setDestructor(stSortedSet *set, void (*destructElement)(void *));

/*
 * Destructs the sorted set. Note: NOT THREAD-SAFE.
 */
void stSortedSet_destruct(stSortedSet *sortedSet);

/*
 * Inserts the object into the sorted set.
 */
void stSortedSet_insert(stSortedSet *sortedSet, void *object);

/*
 * Finds the object in the sorted set, or returns null.
 */
void *stSortedSet_search(stSortedSet *sortedSet, void *object);

/*
 * Finds the object in the the sorted set that is less than or equal to the object, or returns NULL if not found.
 */
void *stSortedSet_searchLessThanOrEqual(stSortedSet *sortedSet, void *object);

/*
 * Finds the object in the the sorted set that is less than the object, or returns NULL if not found.
 */
void *stSortedSet_searchLessThan(stSortedSet *sortedSet, void *object);

/*
 * Finds the object in the the sorted set that is greater than or equal to the object, or returns NULL if not found.
 */
void *stSortedSet_searchGreaterThanOrEqual(stSortedSet *sortedSet, void *object);

/*
 * Finds the object in the the sorted set that is greater than the object, or returns NULL if not found.
 */
void *stSortedSet_searchGreaterThan(stSortedSet *sortedSet, void *object);

/*
 * Removes and returns the object from the sorted set. NOTE: Does
 * *not* call the destruct function on the object.
 */
void *stSortedSet_remove(stSortedSet *sortedSet, void *object);

/*
 * Gets the number of elements in the sorted set.
 */
int64_t stSortedSet_size(stSortedSet *sortedSet);

/*
 * Gets the first element (with lowest value), in the sorted set.
 */
void *stSortedSet_getFirst(stSortedSet *items);

/*
 * Gets the last element in the sorted set.
 */
void *stSortedSet_getLast(stSortedSet *items);

/*
 * Constructs an iterator for the sorted set.
 */
stSortedSetIterator *stSortedSet_getIterator(stSortedSet *items);

/*
 * Gets an iterator from the given object. Creates an error if the item is not in the set.
 * The first value returns by the iterator will be the given item.
 */
stSortedSetIterator *stSortedSet_getIteratorFrom(stSortedSet *items, void *item);

/*
 * Returns an iterator so that the first call to previous will return the last member of the set.
 */
stSortedSetIterator *stSortedSet_getReverseIterator(stSortedSet *items);

/*
 * Destructs an iterator for the sorted set.
 */
void stSortedSet_destructIterator(stSortedSetIterator *iterator);

/*
 * Gets next element in the sorted set.
 */
void *stSortedSet_getNext(stSortedSetIterator *iterator);

/*
 * Gets the previous element in the sorted set.
 */
void *stSortedSet_getPrevious(stSortedSetIterator *iterator);

/*
 * Copies the iterator.
 */
stSortedSetIterator *stSortedSet_copyIterator(stSortedSetIterator *iterator);

/*
 * Gets a stList version of the sorted set, sorted in the order of the sorted set.
 * No destructor is defined for the list, so destroying the list will not destroy
 * the elements in it or the sorted set.
 */
stList *stSortedSet_getList(stSortedSet *sortedSet);

/*
 * Returns non-zero iff the two sets contain the same set of elements, in the same order (i.e. under the same comparison function).
 */
int stSortedSet_equals(stSortedSet *sortedSet1, stSortedSet *sortedSet2);

/*
 * Get the union of two sorted sets. Creates exception if they have different comparators.
 */
stSortedSet *stSortedSet_getUnion(stSortedSet *sortedSet1, stSortedSet *sortedSet2);

/*
 * Get the intersection of two sorted sets. Creates exception if they have different comparators.
 */
stSortedSet *stSortedSet_getIntersection(stSortedSet *sortedSet1, stSortedSet *sortedSet2);

/*
 * Get the set difference of sortedSet1 \ sortedSet2. Creates exception if they have different comparators.
 */
stSortedSet *stSortedSet_getDifference(stSortedSet *sortedSet1, stSortedSet *sortedSet2);

#ifdef __cplusplus
}
#endif
#endif

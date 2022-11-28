/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibstList.h
 *
 *  Created on: 24 May 2010
 *      Author: benedictpaten
 */

#ifndef SONLIB_LIST_H_
#define SONLIB_LIST_H_

#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Construct a stList with zero length.
 * The destructor will not clean up the elements in the stList.
 */
stList *stList_construct(void);

/*
 * Construct a stList with size length.
 * The destructor will not clean up the elements in the stList.
 */
stList *stList_construct2(int64_t size);

/*
 * Construct a stList with size length.
 * The destructor will call the given destructElement function
 * for each non-null entry in the stList.
 */
stList *stList_construct3(int64_t size, void(*destructElement)(void *));

/*
 * Destructs the stList and, if a destructElement function was given to the constructor,
 * calls the destruct element function for each non-null element in the stList.
 * The list maybe NULL.
 */
void stList_destruct(stList *list);

/*
 * Returns the number of elements in the stList.  The list maybe NULL,
 * in which case zero is returned.  This allows for NULL to be used
 * as an efficient way of returning an empty list.
 */
int64_t stList_length(stList *list);

/*
 * Gets item 0 <= index < stList_length(list) from the stList.
 */
void *stList_get(stList *list, int64_t index);

/*
 * Sets the item at that position in the stList.
 */
void stList_set(stList *list, int64_t index, void *item);

/*
 * Adds the item to the end of the st_list, resizing if needed.
 */
void stList_append(stList *list, void *item);

/*
 * Adds all the elements in the second st_list to the end of the first, in order.
 */
void stList_appendAll(stList *stListToAddTo, stList *stListToAdd);

/*
 * Returns the last element in the stList. Error if stList is empty.
 */
void *stList_peek(stList *list);

/*
 * Removes the last element in the stList and returns it.
 * Error if the stList is empty.
 */
void *stList_pop(stList *list);

/*
 * Removes and returns the item at the given index, returning the given item.
 */
void *stList_remove(stList *list, int64_t index);

/*
 * Removes any the first instance of this item from the given stList.
 */
void stList_removeItem(stList *list, void *item);

/*
 * Removes the first element in the stList and returns it. Creates an error if empty.
 */
void *stList_removeFirst(stList *list);

/*
 * Removes the elements in the given interval.
 */
void stList_removeInterval(stList *list, int64_t start, int64_t length);

/*
 * Find the first index of item in the list.
 */
int64_t stList_find(stList *list, void *item);

/*
 * Returns non-zero iff the stList contain one or more copies of references to the given item.
 */
int64_t stList_contains(stList *list, void *item);

/*
 * Copies the stList. Sets the given destruct item function to the new stList.. can
 * be null if you want no destruction of the items in that stList.
 */
stList *stList_copy(stList *list, void(*destructItem)(void *));

/*
 * Reverses the stList in place.
 */
void stList_reverse(stList *list);

/*
 * Gets an iterator for the stList.  The list maybe NULL, in which case a
 * iterator that only returns NULL is created..  This allows for NULL to be
 * used as an efficient way of returning an empty list.
 */
stListIterator *stList_getIterator(stList *list);

/*
 * Destruct the stList iterator.
 */
void stList_destructIterator(stListIterator *iterator);

/*
 * Gets the next item from the iterator.
 */
void *stList_getNext(stListIterator *iterator);

/*
 * Gets the previous item from the iterator.
 */
void *stList_getPrevious(stListIterator *iterator);

/*
 * Copies the iterator.
 */
stListIterator *stList_copyIterator(stListIterator *iterator);

/*
 * Sorts the stList with the given cmpFn.
 * NOTE: cmpFn takes pointers to elements, not pointer to pointers as with
 * C sort.
 */
void stList_sort(stList *list, int (*cmpFn)(const void *a, const void *b));

/*
 * Sorts the stList with the given cmpFn, passing the extra argument to the comparison function.
 * NOTE: cmpFn takes pointers to elements, not pointer to pointers as with
 * C sort.
 */
void stList_sort2(stList *list, int (*cmpFn)(const void *a, const void *b, void *extraArg), void *extraArg);

/*
 * For a sorted list, performs binary search to find the index of a given item.
 * Item is the first argument passed to the cmpFn at each step. If item is is not in the list returns -1;
 */
int64_t stList_binarySearchIndex(stList *list, void *item, int (*cmpFn)(const void *a, const void *b));

/*
 * For a sorted list, performs binary search to find an instance of a given item.
 * Item is the first argument passed to the cmpFn at each step. If item is not in the list returns NULL.
 */
void *stList_binarySearch(stList *list, void *item, int (*cmpFn)(const void *a, const void *b));

/*
 * For a sorted list, performs binary search to find the index of the first occurrence of an item in the list.
 * Item is the first argument passed to the cmpFn at each step. If item is not in the list returns -1
 */
int64_t stList_binarySearchFirstIndex(stList *list, void *item, int (*cmpFn)(const void *a, const void *b));

/*
 * Permutes the list, by iterating over each element and swapping it randomly with a new location.
 */
void stList_shuffle(stList *list);

/*
 * Gets a sorted set representation of the stList, using the given cmpFn as backing. The sorted set
 * has no defined destruct element function, so when the sorted set is destructed the elements in it and
 * in this list will not be destructed. If the cmpFn is NULL then we use the default cmpFn.
 */
stSortedSet *stList_getSortedSet(stList *list,
        int(*cmpFn)(const void *a, const void *b));

/*
 * Returns a set of the elements in a list.
 */
stSet *stList_getSet(stList *list);

/*
 * Converts list to sorted set, destroying old list in process, but transferring the destructor to the set.
 */
stSortedSet *stList_convertToSortedSet(stList *list);

/*
 * Returns a new list with elements that return 0 for the given function removed.
 */
stList *stList_filter(stList *list, bool(*fn)(void *));

/*
 * As above, but passing along extra arg as second argument to filter fn.
 */
stList *stList_filter2(stList *list, bool(*fn)(void *, void *), void *extraArg);

/*
 * Returns a new list, identical to list, but with any elements contained in set removed.
 */
stList *stList_filterToExclude(stList *list, stSortedSet *set);

/*
 * Returns a new list, identical to list, but with any elements not contained in set removed.
 */
stList *stList_filterToInclude(stList *list, stSortedSet *set);

/*
 * Sets the destructor of the list.
 */
void stList_setDestructor(stList *list, void(*destructElement)(void *));

/*
 * Returns new list which contains elements of the list of list concatenated in one list.
 */
stList *stList_join(stList *listOfLists);

/*
 * Replaces each element in the list by applying the mapFn to the element, also passing in the extraArg argument as the second list.
 */
void stList_mapReplace(stList *l, void *(*mapFn)(void *, void *), void *extraArg);

/*
 * Returns the underlying array.
 */
void *stList_getBackingArray(stList *list);

#ifdef __cplusplus
}
#endif
#endif /* SONLIBLIST_H_ */

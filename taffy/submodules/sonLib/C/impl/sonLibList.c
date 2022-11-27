/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibList.c
 *
 *  Created on: 24 May 2010
 *      Author: benedictpaten
 */

#include "safesort.h"
#include "sonLibGlobalsInternal.h"

#define MINIMUM_ARRAY_EXPAND_SIZE 5 //The minimum amount to expand the array backing a list by when it is rescaled.

/*
 * The actual datastructures backing the list
 */

/*
 * The functions..
 */

stList *stList_construct(void) {
    return stList_construct3(0, NULL);
}

stList *stList_construct2(int64_t size) {
    return stList_construct3(size, NULL);
}

stList *stList_construct3(int64_t length, void (*destructElement)(void *)) {
    assert(length >= 0);
    stList *list = st_malloc(sizeof(stList));
    list->length = length;
    list->maxLength = length;
    list->list = st_calloc(length, sizeof(void *));
    list->destructElement = destructElement;
    return list;
}

/* free elements in list */
static void destructElements(stList *list) {
    for(int64_t i=0; i<stList_length(list); i++) { //only free up to known area of list
        if(stList_get(list, i) != NULL) {
            list->destructElement(stList_get(list, i));
        }
    }
}

void stList_destruct(stList *list) {
    if (list != NULL) {
        if (list->destructElement != NULL) {
            destructElements(list);
        }
        if(list->list != NULL) {
            free(list->list);
        }
        free(list);
    }
}

int64_t stList_length(stList *list) {
    if (list == NULL) {
        return 0;
    } else {
        return list->length;
    }
}

void *stList_get(stList *list, int64_t index) {
    assert(index >= 0);
    assert(index < stList_length(list));
    return list->list[index];
}

void stList_set(stList *list, int64_t index, void *item) {
    assert(index >= 0);
    assert(index < stList_length(list));
    list->list[index] = item;
}

void stList_append(stList *list, void *item) {
    if(stList_length(list) >= list->maxLength) {
        list->maxLength = list->maxLength*1.3 + MINIMUM_ARRAY_EXPAND_SIZE;
        list->list = st_realloc(list->list, list->maxLength * sizeof(void *));
    }
    list->list[list->length++] = item;
}

void stList_appendAll(stList *stListToAddTo, stList *stListToAdd) {
    int64_t i;
    assert(stListToAdd != stListToAddTo);
    for(i=0; i<stList_length(stListToAdd); i++) {
        stList_append(stListToAddTo, stList_get(stListToAdd, i));
    }
}

void *stList_peek(stList *list) {
    assert(stList_length(list) > 0);
    return stList_get(list, stList_length(list)-1);
}

void *stList_pop(stList *list) {
    // FIXME: this would more natural to use if it return NULL when empty
    return stList_remove(list, stList_length(list)-1);
}

void *stList_remove(stList *list, int64_t index) {
    assert(index >= 0);
    assert(index < stList_length(list));
    void *o = stList_get(list, index);
    int64_t i;
    for(i=index+1; i<stList_length(list); i++) {
        stList_set(list, i-1, stList_get(list, i));
    }
    list->length--;
    return o;
}

void stList_removeInterval(stList *list, int64_t start, int64_t length) {
    assert(start >= 0);
    assert(start + length <= stList_length(list));
    if(length > 0) {
        int64_t i = start;
        int64_t j = start + length;
        while (i < start + length || j < stList_length(list)) {
            // free removed elements
            if (list->destructElement != NULL && i < start+length) {
                list->destructElement(stList_get(list, i));
            }
            // either move j to i (interval removed from start or middle of list), or clear i (interval at end of list)
            if (j < stList_length(list)) {
                stList_set(list, i, stList_get(list, j));
            } else {
                stList_set(list, i, NULL);
            }
            ++i; ++j;
        }


//        int64_t i = start;
//        for (int64_t j = start+length; j < stList_length(list); j++) {
//            if (list->destructElement != NULL && i < start+length) {
//                void* element = stList_get(list, i);
//                list->destructElement(element);
//            }
//            stList_set(list, i++, stList_get(list, j));
//        }
        list->length -= length;
    }
}

void stList_removeItem(stList *list, void *item)  {
    int64_t i;
    for(i=0; i<stList_length(list); i++) {
        if(stList_get(list, i) == item) {
            stList_remove(list, i);
            return;
        }
    }
}

void *stList_removeFirst(stList *list) {
    return stList_remove(list, 0);
}

int64_t stList_find(stList *list, void *item) {
    int64_t i;
    for (i = 0; i < stList_length(list); i++) {
        if (stList_get(list, i) == item) {
            return i;
        }
    }
    return -1;
}

int64_t stList_contains(stList *list, void *item) {
    return stList_find(list, item) != -1;
}

stList *stList_copy(stList *list, void (*destructItem)(void *)) {
    stList *list2 = stList_construct3(0, destructItem);
    stList_appendAll(list2, list);
    return list2;
}

void stList_reverse(stList *list) {
    int64_t i, j = stList_length(list);
    for(i=0; i<j/2; i++) {
        void *o = stList_get(list, j - 1 - i);
        stList_set(list, j - 1 - i, stList_get(list, i));
        stList_set(list, i, o);
    }
}

stListIterator *stList_getIterator(stList *list) {
    stListIterator *it = st_malloc(sizeof(stListIterator));
    it->list = list;
    it->index = 0;
    return it;
}

void stList_destructIterator(stListIterator *iterator) {
    free(iterator);
}

void *stList_getNext(stListIterator *iterator) {
    if ((iterator->list == NULL) || (iterator->index >= stList_length(iterator->list))) {
        return NULL;
    } else {
        return stList_get(iterator->list, iterator->index++);
    }
}

void *stList_getPrevious(stListIterator *iterator) {
    if ((iterator->list == NULL) || (iterator->index == 0)) {
        return NULL;
    } else {
        return stList_get(iterator->list, --iterator->index);
    }
}

stListIterator *stList_copyIterator(stListIterator *iterator) {
    stListIterator *it = stList_getIterator(iterator->list);
    it->index = iterator->index;
    return it;
}

/* must warp function pointer in data, as ISO C doesn't allow conversions
 * between void* and function pointers */
struct sortFuncArgs {
    int (*cmpFn)(const void *a, const void *b);
};

/* converts pointer to pointer into pointer to element */
static int sortListCmpFn(const void *a, const void *b, void* sargs) {
    return ((struct sortFuncArgs*)sargs)->cmpFn(*(char**)a, *(char**)b);
}

void stList_sort(stList *list, int (*cmpFn)(const void *a, const void *b)) {
    struct sortFuncArgs sargs = {cmpFn};
    safesort(list->list, stList_length(list), sizeof(void *), sortListCmpFn, &sargs);
}

/* must warp function pointer in data, as ISO C doesn't allow conversions
 * between void* and function pointers */
struct sort2FuncArgs {
    int (*cmpFn)(const void *a, const void *b, void* args);
    void *args;
};

/* converts pointer to pointer into pointer to element */
static int sortList2CmpFn(const void *a, const void *b, void* args) {
    struct sort2FuncArgs* sargs = (struct sort2FuncArgs*)args;
    return sargs->cmpFn(*(char**)a, *(char**)b, sargs->args);
}

void stList_sort2(stList *list, int (*cmpFn)(const void *a, const void *b, void *extraArg), void *extraArg) {
    struct sort2FuncArgs sargs = {cmpFn, extraArg};
    safesort(list->list, stList_length(list), sizeof(void *), sortList2CmpFn, &sargs);
}

int64_t stList_binarySearchIndex(stList *list, void *item, int (*cmpFn)(const void *a, const void *b)) {
    int64_t l=0, h=stList_length(list); // interval (l, h) that item can be in, l is inclusive, h is exclusive
    while(l < h) {
        int64_t m = (l + h) / 2; // Mid point
        void *listItem = stList_get(list, m);
        int i = cmpFn(item, listItem);
        if(i < 0) { // Item must occur before m in the list
            h = m;
        }
        else if(i > 0) { // Item must occur after m in the list
            l = m+1;
        } else { // else item at index i equals i
            return m;
        }
    }
    return -1;
}

void *stList_binarySearch(stList *list, void *item, int (*cmpFn)(const void *a, const void *b)) {
    int64_t i = stList_binarySearchIndex(list, item, cmpFn);
    return i == -1 ? NULL : stList_get(list, i);
}

int64_t stList_binarySearchFirstIndex(stList *list, void *item, int (*cmpFn)(const void *a, const void *b)) {
    int64_t i = stList_binarySearchIndex(list, item, cmpFn);
    while(i > 0 && cmpFn(item, stList_get(list, i-1)) == 0) {
        i--;
    }
    return i;
}

void stList_shuffle(stList *list) {
    for(int64_t i=0; i<stList_length(list); i++) {
        int64_t j = st_randomInt(i, stList_length(list));
        void *o = stList_get(list, i);
        stList_set(list, i, stList_get(list, j));
        stList_set(list, j, o);
    }
}

stSortedSet *stList_getSortedSet(stList *list, int (*cmpFn)(const void *a, const void *b)) {
    stSortedSet *sortedSet = stSortedSet_construct3(cmpFn, NULL);
    int64_t i;
    for(i=0; i<stList_length(list); i++) {
        stSortedSet_insert(sortedSet, stList_get(list, i));
    }
    return sortedSet;
}

stSet *stList_getSet(stList *list) {
    stSet *set = stSet_construct();
    for(int64_t i=0; i<stList_length(list); i++) {
        stSet_insert(set, stList_get(list, i));
    }
    return set;
}

void stList_setDestructor(stList *list, void (*destructElement)(void *)) {
    list->destructElement = destructElement;
}

stList *stList_filter(stList *list, bool(*fn)(void *)) {
    stList *list2 = stList_construct();
    for (int64_t i = 0; i < stList_length(list); i++) {
        void *o = stList_get(list, i);
        if (fn(o)) {
            stList_append(list2, o);
        }
    }
    return list2;
}

stList *stList_filter2(stList *list, bool(*fn)(void *, void *), void *extraArg) {
    stList *list2 = stList_construct();
    for (int64_t i = 0; i < stList_length(list); i++) {
        void *o = stList_get(list, i);
        if (fn(o, extraArg)) {
            stList_append(list2, o);
        }
    }
    return list2;
}

bool filterToExcludeP(void *element, void *set) {
    return stSortedSet_search(set, element) == NULL;
}

stList *stList_filterToExclude(stList *list, stSortedSet *set) {
    return stList_filter2(list, filterToExcludeP, set);
}

bool filterToIncludeP(void *element, void *set) {
    return stSortedSet_search(set, element) != NULL;
}

stList *stList_filterToInclude(stList *list, stSortedSet *set) {
    return stList_filter2(list, filterToIncludeP, set);
}

stList *stList_join(stList *listOfLists) {
    stList *joinedList = stList_construct();
    for (int64_t i = 0; i < stList_length(listOfLists); i++) {
        stList_appendAll(joinedList, stList_get(listOfLists, i));
    }
    return joinedList;
}

stSortedSet *stList_convertToSortedSet(stList *list) {
    stSortedSet *set = stList_getSortedSet(list, NULL);
    stSortedSet_setDestructor(set, list->destructElement);
    stList_setDestructor(list, NULL);
    stList_destruct(list);
    return set;
}

void stList_mapReplace(stList *l, void *(*mapFn)(void *, void *), void *extraArg) {
    int64_t j = stList_length(l);
    for(int64_t i=0; i<j; i++) {
        stList_set(l, i, mapFn(stList_get(l, i), extraArg));
    }
}

void *stList_getBackingArray(stList *list) {
    return list->list;
}

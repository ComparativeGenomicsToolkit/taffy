/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLibGlobalsInternal.h"
#include "avl.h"

const char *SORTED_SET_EXCEPTION_ID = "SORTED_SET_EXCEPTION";

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions on a sorted set and its iterator
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct _stSortedSet {
    struct avl_table *sortedSet;
    void (*destructElementFn)(void *);
    //int numberOfLiveIterators;  // number of currently allocated iterators
};

struct _stSortedSetIterator {
    stSortedSet *sortedSet;
    struct avl_traverser traverser;
};

static int st_sortedSet_cmpFn( const void *key1, const void *key2 ) {
    return key1 > key2 ? 1 : key1 < key2 ? -1 : 0;
}

/* check if the sorted set is in a state that it can be modified */
static void checkModifiable(stSortedSet *sortedSet) {
    //if (sortedSet->numberOfLiveIterators > 0) {
        //assert(0);
        //stThrowNew(SORTED_SET_EXCEPTION_ID, "attempt to modify an stSortedSet while iterators are active");
    //}
}

stSortedSet *stSortedSet_construct(void) {
    return stSortedSet_construct3(st_sortedSet_cmpFn, NULL);
}

stSortedSet *stSortedSet_construct2(void (*destructElementFn)(void *)) {
    return stSortedSet_construct3(st_sortedSet_cmpFn, destructElementFn);
}

struct _stSortedSet_construct3Fn {
    int (*compareFn)(const void *, const void *);
};

static int st_sortedSet_construct3P(const void *a, const void *b, struct _stSortedSet_construct3Fn *c) { //int (*cmpFn)(const void *, const void *)) {
    return c->compareFn(a, b);
}

stSortedSet *stSortedSet_construct3(int (*compareFn)(const void *, const void *),
                                      void (*destructElementFn)(void *)) {
    stSortedSet *sortedSet = st_malloc(sizeof(stSortedSet));
    struct _stSortedSet_construct3Fn *i = st_malloc(sizeof(struct _stSortedSet_construct3Fn));
    i->compareFn = compareFn == NULL ? st_sortedSet_cmpFn : compareFn; //this is a total hack to make the function pass ISO C compatible.
    sortedSet->sortedSet = avl_create((int (*)(const void *, const void *, void *))st_sortedSet_construct3P, i, NULL);
    sortedSet->destructElementFn = destructElementFn;
    //sortedSet->numberOfLiveIterators = 0;
    return sortedSet;
}

void stSortedSet_setDestructor(stSortedSet *set, void (*destructElement)(void *)) {
    set->destructElementFn = destructElement;
}

static struct _stSortedSet_construct3Fn *stSortedSet_getComparator(stSortedSet *sortedSet) {
    return (struct _stSortedSet_construct3Fn *)sortedSet->sortedSet->avl_param;
}

stSortedSet *stSortedSet_copyConstruct(stSortedSet *sortedSet, void (*destructElementFn)(void *)) {
    stSortedSet *sortedSet2 = stSortedSet_construct3(stSortedSet_getComparator(sortedSet)->compareFn, destructElementFn);
    stSortedSetIterator *it = stSortedSet_getIterator(sortedSet);
    void *o;
    while((o = stSortedSet_getNext(it)) != NULL) {
        stSortedSet_insert(sortedSet2, o);
    }
    stSortedSet_destructIterator(it);
    return sortedSet2;
}

static void st_sortedSet_destructP(void *a, void *b, void (*extraArg)(void*)) {
    assert(b != NULL);
    extraArg(a);
}

void stSortedSet_destruct(stSortedSet *sortedSet) {
    void *a = sortedSet->sortedSet->avl_param;
    if(sortedSet->destructElementFn != NULL) {
        avl_destroy2(sortedSet->sortedSet, st_sortedSet_destructP, sortedSet->destructElementFn);
    }
    else {
        avl_destroy(sortedSet->sortedSet, NULL);
    }
    free(a);
    free(sortedSet);
}

void stSortedSet_insert(stSortedSet *sortedSet, void *object) {
    checkModifiable(sortedSet);
    // FIXME: two passes, modify avl code.
    if(stSortedSet_search(sortedSet, object) != NULL) {
        avl_replace(sortedSet->sortedSet, object);
    }
    else {
        avl_insert(sortedSet->sortedSet, object);
    }
}

void *stSortedSet_search(stSortedSet *sortedSet, void *object) {
    return avl_find(sortedSet->sortedSet, object);
}

void *stSortedSet_searchLessThanOrEqual(stSortedSet *sortedSet, void *object) {
    return avl_find_lessThanOrEqual(sortedSet->sortedSet, object);
}

void *stSortedSet_searchLessThan(stSortedSet *sortedSet, void *object) {
    return avl_find_lessThan(sortedSet->sortedSet, object);
}

void *stSortedSet_searchGreaterThanOrEqual(stSortedSet *sortedSet, void *object) {
    return avl_find_greaterThanOrEqual(sortedSet->sortedSet, object);
}

void *stSortedSet_searchGreaterThan(stSortedSet *sortedSet, void *object) {
    return avl_find_greaterThan(sortedSet->sortedSet, object);
}

void *stSortedSet_remove(stSortedSet *sortedSet, void *object) {
    checkModifiable(sortedSet);
    return avl_delete(sortedSet->sortedSet, object);
}

int64_t stSortedSet_size(stSortedSet *sortedSet) {
    return avl_count(sortedSet->sortedSet);
}

void *stSortedSet_getFirst(stSortedSet *items) {
    struct avl_traverser traverser;
    avl_t_init(&traverser, items->sortedSet);
    return avl_t_first(&traverser, items->sortedSet);
}

void *stSortedSet_getLast(stSortedSet *items) {
    struct avl_traverser traverser;
    avl_t_init(&traverser, items->sortedSet);
    return avl_t_last(&traverser, items->sortedSet);
}

stSortedSetIterator *stSortedSet_getIterator(stSortedSet *items) {
    stSortedSetIterator *iterator;
    iterator = st_malloc(sizeof(stSortedSetIterator));
    iterator->sortedSet = items;
    avl_t_init(&iterator->traverser, items->sortedSet);
    //items->numberOfLiveIterators++;
    return iterator;
}

stSortedSetIterator *stSortedSet_getIteratorFrom(stSortedSet *items, void *item) {
    stSortedSetIterator *iterator = stSortedSet_getIterator(items);
    if(avl_t_find(&iterator->traverser, items->sortedSet, item) == NULL) {
        stSortedSet_destructIterator(iterator);
        stThrowNew(SORTED_SET_EXCEPTION_ID, "Tried to create an iterator with an item that is not in the list of items");
    }
    stSortedSet_getPrevious(iterator);
    return iterator;
}

stSortedSetIterator *stSortedSet_getReverseIterator(stSortedSet *items) {
    stSortedSetIterator *it = stSortedSet_getIteratorFrom(items, stSortedSet_getLast(items));
    stSortedSet_getNext(it);
    stSortedSet_getNext(it);
    return it;
}

void stSortedSet_destructIterator(stSortedSetIterator *iterator) {
    //assert(iterator->sortedSet->numberOfLiveIterators > 0);
    //iterator->sortedSet->numberOfLiveIterators--;
    free(iterator);
}

void *stSortedSet_getNext(stSortedSetIterator *iterator) {
    return avl_t_next(&iterator->traverser);
}

stSortedSetIterator *stSortedSet_copyIterator(stSortedSetIterator *iterator) {
    stSortedSetIterator *copyIterator;
    copyIterator = st_malloc(sizeof(stSortedSetIterator));
    copyIterator->sortedSet = iterator->sortedSet;
    //copyIterator->sortedSet->numberOfLiveIterators++;
    avl_t_copy(&copyIterator->traverser, &iterator->traverser);
    return copyIterator;
}

void *stSortedSet_getPrevious(stSortedSetIterator *iterator) {
    return avl_t_prev(&iterator->traverser);
}

static int stSortedSet_comparatorsEqual(stSortedSet *sortedSet1, stSortedSet *sortedSet2) {
    return stSortedSet_getComparator(sortedSet1)->compareFn == stSortedSet_getComparator(sortedSet2)->compareFn;
}

int stSortedSet_equals(stSortedSet *sortedSet1, stSortedSet *sortedSet2) {
    if(stSortedSet_size(sortedSet1) != stSortedSet_size(sortedSet2)) {
        return 0;
    }
    if(!stSortedSet_comparatorsEqual(sortedSet1, sortedSet2)) {
        return 0;
    }
    int (*cmpFn)(const void *, const void *) = stSortedSet_getComparator(sortedSet1)->compareFn;

    stSortedSetIterator *it1 = stSortedSet_getIterator(sortedSet1);
    stSortedSetIterator *it2 = stSortedSet_getIterator(sortedSet2);
    void *o1 = stSortedSet_getNext(it1), *o2 = stSortedSet_getNext(it2);
    while(o1 != NULL && o2 != NULL) {
        if(cmpFn(o1, o2) != 0) {
            stSortedSet_destructIterator(it1);
            stSortedSet_destructIterator(it2);
            return 0;
        }
        o1 = stSortedSet_getNext(it1);
        o2 = stSortedSet_getNext(it2);
    }
    stSortedSet_destructIterator(it1);
    stSortedSet_destructIterator(it2);
    return 1;
}

stList *stSortedSet_getList(stSortedSet *sortedSet) {
    stList *list = stList_construct2(stSortedSet_size(sortedSet));
    stSortedSetIterator *it = stSortedSet_getIterator(sortedSet);
    void *o;
    int64_t i=0;
    while((o = stSortedSet_getNext(it)) != NULL) {
        stList_set(list, i++, o);
    }
    assert(i == stSortedSet_size(sortedSet));
    stSortedSet_destructIterator(it);
    return list;
}

stSortedSet *stSortedSet_getUnion(stSortedSet *sortedSet1, stSortedSet *sortedSet2) {
    if(!stSortedSet_comparatorsEqual(sortedSet1, sortedSet2)) {
        stThrowNew(SORTED_SET_EXCEPTION_ID, "Comparators are not equal for creating the union of two sorted sets");
    }
    stSortedSet *sortedSet3 = stSortedSet_construct3(stSortedSet_getComparator(sortedSet1)->compareFn, NULL);

    //Add those from sortedSet1
    stSortedSetIterator *it= stSortedSet_getIterator(sortedSet1);
    void *o;
    while((o = stSortedSet_getNext(it)) != NULL) {
        stSortedSet_insert(sortedSet3, o);
    }
    stSortedSet_destructIterator(it);

    //Add those from sortedSet2
    it= stSortedSet_getIterator(sortedSet2);
    while((o = stSortedSet_getNext(it)) != NULL) {
        stSortedSet_insert(sortedSet3, o);
    }
    stSortedSet_destructIterator(it);

    return sortedSet3;
}

stSortedSet *stSortedSet_getIntersection(stSortedSet *sortedSet1, stSortedSet *sortedSet2) {
    if(!stSortedSet_comparatorsEqual(sortedSet1, sortedSet2)) {
        stThrowNew(SORTED_SET_EXCEPTION_ID, "Comparators are not equal for creating an intersection of two sorted sets");
    }
    stSortedSet *sortedSet3 = stSortedSet_construct3(stSortedSet_getComparator(sortedSet1)->compareFn, NULL);

    //Add those from sortedSet1 only if they are also in sortedSet2
    stSortedSetIterator *it= stSortedSet_getIterator(sortedSet1);
    void *o;
    while((o = stSortedSet_getNext(it)) != NULL) {
        if(stSortedSet_search(sortedSet2, o) != NULL) {
            stSortedSet_insert(sortedSet3, o);
        }
    }
    stSortedSet_destructIterator(it);

    return sortedSet3;
}

stSortedSet *stSortedSet_getDifference(stSortedSet *sortedSet1, stSortedSet *sortedSet2) {
    if(!stSortedSet_comparatorsEqual(sortedSet1, sortedSet2)) {
        stThrowNew(SORTED_SET_EXCEPTION_ID, "Comparators are not equal for creating the sorted set difference");
    }
    stSortedSet *sortedSet3 = stSortedSet_construct3(stSortedSet_getComparator(sortedSet1)->compareFn, NULL);

    //Add those from sortedSet1 only if they are not in sortedSet2
    stSortedSetIterator *it= stSortedSet_getIterator(sortedSet1);
    void *o;
    while((o = stSortedSet_getNext(it)) != NULL) {
        if(stSortedSet_search(sortedSet2, o) == NULL) {
            stSortedSet_insert(sortedSet3, o);
        }
    }
    stSortedSet_destructIterator(it);

    return sortedSet3;
}

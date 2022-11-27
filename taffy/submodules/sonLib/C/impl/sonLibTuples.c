/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibContainers.c
 *
 *  Created on: 26-May-2010
 *      Author: benedictpaten
 */

#include "sonLibGlobalsInternal.h"

stIntTuple *stIntTuple_construct0() {
    stIntTuple *intTuple = st_malloc(sizeof(int64_t) * 1);
    intTuple[0] = 0;
    return intTuple;
}

stIntTuple *stIntTuple_construct1(int64_t value) {
    stIntTuple *intTuple = st_malloc(sizeof(int64_t) * 2);
    intTuple[0] = 1;
    intTuple[1] = value;
    return intTuple;
}

stIntTuple *stIntTuple_construct2(int64_t value1, int64_t value2) {
    stIntTuple *intTuple = st_malloc(sizeof(int64_t) * 3);
    intTuple[0] = 2;
    intTuple[1] = value1;
    intTuple[2] = value2;
    return intTuple;
}

stIntTuple *stIntTuple_construct3(int64_t value1, int64_t value2, int64_t value3) {
    stIntTuple *intTuple = st_malloc(sizeof(int64_t) * 4);
    intTuple[0] = 3;
    intTuple[1] = value1;
    intTuple[2] = value2;
    intTuple[3] = value3;
    return intTuple;
}

stIntTuple *stIntTuple_construct4(int64_t value1, int64_t value2, int64_t value3, int64_t value4) {
    stIntTuple *intTuple = st_malloc(sizeof(int64_t) * 5);
    intTuple[0] = 4;
    intTuple[1] = value1;
    intTuple[2] = value2;
    intTuple[3] = value3;
    intTuple[4] = value4;
    return intTuple;
}

stIntTuple *stIntTuple_construct5(int64_t value1, int64_t value2, int64_t value3, int64_t value4, int64_t value5) {
    stIntTuple *intTuple = st_malloc(sizeof(int64_t) * 6);
    intTuple[0] = 5;
    intTuple[1] = value1;
    intTuple[2] = value2;
    intTuple[3] = value3;
    intTuple[4] = value4;
    intTuple[5] = value5;
    return intTuple;
}

stIntTuple *stIntTuple_constructN(int64_t length, const int64_t iA[]) {
    assert(length >= 0);
    stIntTuple *intTuple = st_malloc(sizeof(int64_t) * (length + 1));
    intTuple[0] = length;
    for(int64_t i=0; i<length; i++) {
        intTuple[i+1] = iA[i];
    }
    return intTuple;
}

void stIntTuple_destruct(stIntTuple *intTuple) {
    free(intTuple);
}

uint64_t stIntTuple_hashKey(stIntTuple *intTuple) {
    int64_t i;
    uint64_t j = 0;

    for (i = 0; i < stIntTuple_length(intTuple); i++) {
        //j += stIntTuple_getPosition(intTuple, i);

        // djb2
        // This algorithm was first reported by Dan Bernstein
        // many years ago in comp.lang.c
        //
        j = stIntTuple_get(intTuple, i) + (j << 6) + (j << 16) - j;
    }
    return j;
}

static int intCmp(int64_t i, int64_t j) {
    return i < j ? -1 : (i > j ? 1 : 0);
}

int stIntTuple_cmpFn(stIntTuple *intTuple1, stIntTuple *intTuple2) {
    int64_t
            i,
            j =
                    stIntTuple_length(intTuple1) < stIntTuple_length(intTuple2) ? stIntTuple_length(
                            intTuple1)
                            : stIntTuple_length(intTuple2);
    for (i = 0; i < j; i++) {
        int64_t k = intCmp(stIntTuple_get(intTuple1, i)
                ,stIntTuple_get(intTuple2, i));
        if (k != 0) {
            return k;
        }
    }
    return intCmp(stIntTuple_length(intTuple1), stIntTuple_length(intTuple2));
}

int stIntTuple_equalsFn(stIntTuple *intTuple1, stIntTuple *intTuple2) {
    return stIntTuple_cmpFn(intTuple1, intTuple2) == 0;
}

int64_t stIntTuple_length(stIntTuple *intTuple) {
    return intTuple[0];
}

int64_t stIntTuple_get(stIntTuple *intTuple, int64_t index) {
    assert(index < stIntTuple_length(intTuple));
    assert(index >= 0);
    return intTuple[index + 1];
}

/*
 * The following are double variants of the above functions.
 */

stDoubleTuple *stDoubleTuple_construct(int64_t length, ...) {
    assert(length >= 0);
    stDoubleTuple *doubleTuple = st_malloc(sizeof(double) * (length + 1));
    doubleTuple[0] = length;
    int64_t i;
    va_list ap;
    va_start(ap, length);
    for (i = 0; i < length; i++) {
        doubleTuple[i + 1] = va_arg(ap, double);
    }
    va_end(ap);
    return doubleTuple;
}

void stDoubleTuple_destruct(stDoubleTuple *doubleTuple) {
    free(doubleTuple);
}

uint64_t stDoubleTuple_hashKey(stDoubleTuple *doubleTuple) {
    int64_t i;
    uint64_t j = 0;

    for (i = 0; i < stDoubleTuple_length(doubleTuple); i++) {
        //j += stIntTuple_getPosition(intTuple, i);

        // djb2
        // This algorithm was first reported by Dan Bernstein
        // many years ago in comp.lang.c
        //
        j = stDoubleTuple_getPosition(doubleTuple, i) + (j << 6) + (j << 16) - j;
    }
    return j;
}

static int doubleCmp(double i, double j) {
    return i < j ? -1 : (i > j ? 1 : 0);
}

int stDoubleTuple_cmpFn(stDoubleTuple *doubleTuple1,
        stDoubleTuple *doubleTuple2) {
    int64_t i, j = stDoubleTuple_length(doubleTuple1) < stDoubleTuple_length(
            doubleTuple2) ? stDoubleTuple_length(doubleTuple1)
            : stDoubleTuple_length(doubleTuple2);
    for (i = 0; i < j; i++) {
        int64_t k = doubleCmp(stDoubleTuple_getPosition(doubleTuple1, i),
                stDoubleTuple_getPosition(doubleTuple2, i));
        if (k != 0) {
            return k;
        }
    }
    return intCmp(stDoubleTuple_length(doubleTuple1), stDoubleTuple_length(
            doubleTuple2));
}

int stDoubleTuple_equalsFn(stDoubleTuple *doubleTuple1,
        stDoubleTuple *doubleTuple2) {
    return stDoubleTuple_cmpFn(doubleTuple1, doubleTuple2) == 0;
}

int64_t stDoubleTuple_length(stDoubleTuple *doubleTuple) {
    return (int64_t) doubleTuple[0];
}

double stDoubleTuple_getPosition(stDoubleTuple *doubleTuple, int64_t index) {
    assert(index < stDoubleTuple_length(doubleTuple));
    assert(index >= 0);
    return doubleTuple[index + 1];
}


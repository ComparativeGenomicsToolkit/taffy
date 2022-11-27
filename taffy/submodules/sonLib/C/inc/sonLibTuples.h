/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibTuples.h
 *
 *  Created on: 26-May-2010
 *      Author: benedictpaten
 */

#ifndef SONLIB_TUPLES_H_
#define SONLIB_TUPLES_H_

#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Constructs a tuple of length int64_t integers.. be very careful that length equals the
 * number of subsequent arguments..
 */

stIntTuple *stIntTuple_construct0(void);

stIntTuple *stIntTuple_construct1(int64_t value);

stIntTuple *stIntTuple_construct2(int64_t value1, int64_t value2);

stIntTuple *stIntTuple_construct3(int64_t value1, int64_t value2, int64_t value3);

stIntTuple *stIntTuple_construct4(int64_t value1, int64_t value2, int64_t value3, int64_t value4);

stIntTuple *stIntTuple_construct5(int64_t value1, int64_t value2, int64_t value3, int64_t value4, int64_t value5);

stIntTuple *stIntTuple_constructN(int64_t length, const int64_t iA[]);

/*
 * Destructs the tuple.
 */
void stIntTuple_destruct(stIntTuple *intTuple);

/*
 * Creates a hash key for the tuple.
 */
uint64_t stIntTuple_hashKey(stIntTuple *intTuple);


/*
 * Compares two int tuples, comparing the first members of the tuple, then
 * the second etc.. returning if it finds a difference,
 * until one or other runs out of length. If they are equal
 * for all common positions but one is longer than the other, then the shorted
 * is deemed less than the longer tuple.
 */
int stIntTuple_cmpFn(stIntTuple *intTuple1, stIntTuple *intTuple2);

/*
 * Returns non zero iff stIntTuple_cmpFn(intTuple1, intTuple2) == 0.
 */
int stIntTuple_equalsFn(stIntTuple *intTuple1, stIntTuple *intTuple2);

/*
 * Returns the length of the tuple.
 */
int64_t stIntTuple_length(stIntTuple *intTuple);

/*
 * Returns the value of a position in the tuple, 0 <= index < stIntTuple_length(tuple).
 */
int64_t stIntTuple_get(stIntTuple *intTuple, int64_t index);

/*
 * The following are double variants of the above functions.
 *  One must be very careful to ensure that the variable arguments are of type double!
 */

stDoubleTuple *stDoubleTuple_construct(int64_t length, ...);

void stDoubleTuple_destruct(stDoubleTuple *doubleTuple);

uint64_t stDoubleTuple_hashKey(stDoubleTuple *doubleTuple);

int stDoubleTuple_cmpFn(stDoubleTuple *doubleTuple1, stDoubleTuple *doubleTuple2);

int stDoubleTuple_equalsFn(stDoubleTuple *doubleTuple1, stDoubleTuple *doubleTuple2);

int64_t stDoubleTuple_length(stDoubleTuple *doubleTuple);

double stDoubleTuple_getPosition(stDoubleTuple *doubleTuple, int64_t index);

#ifdef __cplusplus
}
#endif
#endif /* SONLIBWRAPPERS_H_ */

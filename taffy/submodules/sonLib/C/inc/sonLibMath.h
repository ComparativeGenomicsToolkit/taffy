/*
 * Copyright (C) 2012 by Benedict Paten (a) gmail com
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef SONLIB_MATH_H_
#define SONLIB_MATH_H_

/*
 * sonLibMath.h
 */

#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#define ST_MATH_LOG_ZERO -INFINITY

/*
 * Returns a reasonable approximation to log(exp(x) + exp(y)).
 */
double stMath_logAdd(double x, double y);

/*
 * Returns a more exact calculation log(exp(x) + exp(y)).
 */
double stMath_logAddExact(double x, double y);

#ifdef __cplusplus
}
#endif // __cplusplus
#endif // SONLIB_MATH_H_

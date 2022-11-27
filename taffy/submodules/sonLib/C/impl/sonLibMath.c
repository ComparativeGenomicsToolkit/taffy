/*
 * sonListPrivate.h
 *
 *  Created on: 28 Jun 2012
 *      Author: benedictpaten
 */

#include "sonLibGlobalsInternal.h"
#include <math.h>

#define logUnderflowThreshold 7.5

static inline double lookup(double x) {
    //return log (exp (x) + 1);
    assert(x >= 0.00f);
    assert(x <= logUnderflowThreshold);
    if (x <= 1.00f)
        return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
    if (x <= 2.50f)
        return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
    if (x <= 4.50f)
        return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;
    return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;
}

double stMath_logAdd(double x, double y) {
    if (x < y)
        return (x == ST_MATH_LOG_ZERO || y - x >= logUnderflowThreshold) ? y : lookup(y - x) + x;
    return (y == ST_MATH_LOG_ZERO || x - y >= logUnderflowThreshold) ? x : lookup(x - y) + y;
}

static inline double lookupExact(double x) {
    return log (exp (x) + 1);
}

double stMath_logAddExact(double x, double y) {
    if (x < y)
        return (x == ST_MATH_LOG_ZERO || y - x >= logUnderflowThreshold) ? y : lookupExact(y - x) + x;
    return (y == ST_MATH_LOG_ZERO || x - y >= logUnderflowThreshold) ? x : lookupExact(x - y) + y;
}

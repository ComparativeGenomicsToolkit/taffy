/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*This code is taken and adapted directly from the probcons library.
 */
 
/////////////////////////////////////////////////////////////////
// FLOAT_32.h
//
// Routines for doing math operations in PROBCONS.
/////////////////////////////////////////////////////////////////
#include <math.h>
#include <float.h>   
#include <stdlib.h>
#include <assert.h> 

#include "fastCMaths.h" 
#include "sonLib.h"

//typedef FLOAT_32 FLOAT_32; 

/////////////////////////////////////////////////////////////////
// LOG()
//
// Compute the logarithm of x.
/////////////////////////////////////////////////////////////////

float LOG (float x){
  return log (x); 
}

/////////////////////////////////////////////////////////////////
// EXP()
//
// Computes exp(x).
/////////////////////////////////////////////////////////////////

float EXP (float x){
  return exp(x);
  if (x > -2){
    if (x > -0.5){
      if (x > 0) {
        return exp(x);
      }
      return (((0.03254409303190190000*x + 0.16280432765779600000)*x + 0.49929760485974900000)*x + 0.99995149601363700000)*x + 0.99999925508501600000;
    }
    if (x > -1) {
      return (((0.01973899026052090000*x + 0.13822379685007000000)*x + 0.48056651562365000000)*x + 0.99326940370383500000)*x + 0.99906756856399500000;
    }
    return (((0.00940528203591384000*x + 0.09414963667859410000)*x + 0.40825793595877300000)*x + 0.93933625499130400000)*x + 0.98369508190545300000;
  }
  if (x > -8){
    if (x > -4) {
      return (((0.00217245711583303000*x + 0.03484829428350620000)*x + 0.22118199801337800000)*x + 0.67049462206469500000)*x + 0.83556950223398500000;
    }
    return (((0.00012398771025456900*x + 0.00349155785951272000)*x + 0.03727721426017900000)*x + 0.17974997741536900000)*x + 0.33249299994217400000;
  }
  if (x > -16) {
    return (((0.00000051741713416603*x + 0.00002721456879608080)*x + 0.00053418601865636800)*x + 0.00464101989351936000)*x + 0.01507447981459420000;
  }
  return 0;
} 

/*const FLOAT_32 EXP_UNDERFLOW_THRESHOLD = -4.6;
const FLOAT_32 LOG_UNDERFLOW_THRESHOLD = 7.5;*/

/////////////////////////////////////////////////////////////////
// LOOKUP()
//
// Computes log (exp (x) + 1), for 0 <= x <= 7.5.
/////////////////////////////////////////////////////////////////

float LOOKUP (float x){
  //return log (exp (x) + 1);
  assert (x >= 0.00f);
  assert (x <= LOG_UNDERFLOW_THRESHOLD);
  if (x <= 1.00f) return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
  if (x <= 2.50f) return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
  if (x <= 4.50f) return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;
  assert (x <= LOG_UNDERFLOW_THRESHOLD);
  return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;
}

/////////////////////////////////////////////////////////////////
// LOG_PLUS_EQUALS()
//
// Add two log probabilities and store in the first argument
/////////////////////////////////////////////////////////////////

void LOG_PLUS_EQUALS (float *x, float y){
  if (*x < y)
    *x = (*x == LOG_ZERO || y - *x >= LOG_UNDERFLOW_THRESHOLD) ? y : LOOKUP(y-*x) + *x;
  else
    *x = (y == LOG_ZERO || *x - y >= LOG_UNDERFLOW_THRESHOLD) ? *x : LOOKUP(*x-y) + y;
}


/////////////////////////////////////////////////////////////////
// LOG_ADD()
//
// Add two log probabilities
/////////////////////////////////////////////////////////////////

float LOG_ADD (float x, float y){
  if (x < y) return (x == LOG_ZERO || y - x >= LOG_UNDERFLOW_THRESHOLD) ? y : LOOKUP(y-x) + x;
  return (y == LOG_ZERO || x - y >= LOG_UNDERFLOW_THRESHOLD) ? x : LOOKUP(x-y) + y;
}


/////////////////////////////////////////////////////////////////
// LOG_ADD()
//
// Add three log probabilities
/////////////////////////////////////////////////////////////////

float LOG_ADD_THREE (float x1, float x2, float x3){
  return LOG_ADD (x1, LOG_ADD (x2, x3));
}

/////////////////////////////////////////////////////////////////
// MAX_EQUALS()
//
// Chooses maximum of two arguments and stores it in the first argument
/////////////////////////////////////////////////////////////////

void MAX_PLUS_EQUALS (float *x, float y){
  if ((*x) < y)
    *x = y;
}

/////////////////////////////////////////////////////////////////
// RANDOM()
//
// Return random FLOAT_32 in range [0 - 1.0 }
/////////////////////////////////////////////////////////////////
float RANDOM(void) {
    static const float i = (RAND_MAX+1.0)*1.0f;
    return rand()/i;
}

/////////////////////////////////////////////////////////////////
// RANDOM()
//
// Return random FLOAT_32 in range [0 - 1.0 }
/////////////////////////////////////////////////////////////////
float RANDOM_LOG(void) {
    return LOG(RANDOM());
}

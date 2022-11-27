/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibRandom.h
 *
 *  Created on: 22-Jun-2010
 *      Author: benedictpaten
 */

#ifndef SONLIBRANDOM_H_
#define SONLIBRANDOM_H_

#ifdef __cplusplus
extern "C" {
#endif

//The exception string
extern const char *RANDOM_EXCEPTION_ID;

//////////////////////
//Random number functions
//////////////////////

/*
 * Seed the random number generator.
 */
void st_randomSeed(int64_t seed);

/*
 * Returns a random value in the range min (inclusive) to max (exclusive), where min < max.
 */
int64_t st_randomInt(int64_t min, int64_t max);

/*
 * Like st_randomInt, but for 64 bit integers.
 */
int64_t st_randomInt64(int64_t min, int64_t max);

/*
 * Returns a random value between 0.0 (inclusive) and 1.0 (exclusive).
 */
double st_random(void);

/*
 * Returns a random value from a list.
 */
void *st_randomChoice(stList *list);

/*
 * Gets a random upper-case A,C,G, T character.
 * If includeNs is true, may additionally return an N.
 * If useLowerCase then return the character as lower case, is useRandomCase then may return either lower or upper case.
 */
char stRandom_getRandomNucleotide(bool includeNs, bool useLowerCase, bool useRandomCase);

/*
 * Get a random string of nucleotide characters using st_randomGetRandomNucleotide(), with same arguments.
 */
char *stRandom_getRandomDNAString(int64_t length, bool includeNs, bool useLowerCase, bool useRandomCase);

#ifdef __cplusplus
}
#endif
#endif /* SONLIBRANDOM_H_ */

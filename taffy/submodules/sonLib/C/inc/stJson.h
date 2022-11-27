/*
 * Copyright (C) 2018 by Benedict Paten (a) gmail com
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef SONLIB_JSON_H_
#define SONLIB_JSON_H_

/*
 * Json parsing helper functions
 */

#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/*
 * Convert json token into a string.
 */
char *stJson_token_tostr(char *js, jsmntok_t *t);

int64_t stJson_parseInt(char *js, jsmntok_t *tokens, int64_t tokenIndex);

double stJson_parseFloat(char *js, jsmntok_t *tokens, int64_t tokenIndex);

int64_t stJson_parseFloatArray(double *toArray, int64_t toArraySize, char *js, jsmntok_t *tokens, int64_t tokenIndex);

bool stJson_parseBool(char *js, jsmntok_t *tokens, int64_t tokenIndex);

/*
 * Returns the number of nested tokens for a given key. 
 */
int64_t stJson_getNestedTokenCount(jsmntok_t *tokens, int64_t tokenIndex);

/*
 * Takes a char buffer (buf) of size r and allocates a set of tokens and mutable string copy js 
 * of the buffer.
 */
size_t stJson_setupParser(char *buf, size_t r, jsmntok_t **tokens, char **js);

#ifdef __cplusplus
}
#endif // __cplusplus
#endif // SONLIB_JSON_H_

#include "sonLibGlobalsInternal.h"

/*
 * Json parsing helper functions
 */

/*
 * Convert json token into a string.
 */
char *stJson_token_tostr(char *js, jsmntok_t *t) {
    js[t->end] = '\0';
    return js + t->start;
}

int64_t stJson_parseInt(char *js, jsmntok_t *tokens, int64_t tokenIndex) {
	jsmntok_t tok = tokens[tokenIndex];
	char *tokStr = stJson_token_tostr(js, &tok);
	return atoi(tokStr);
}

double stJson_parseFloat(char *js, jsmntok_t *tokens, int64_t tokenIndex) {
	jsmntok_t tok = tokens[tokenIndex];
	char *tokStr = stJson_token_tostr(js, &tok);
	return atof(tokStr);
}

int64_t stJson_parseFloatArray(double *toArray, int64_t toArraySize, char *js, jsmntok_t *tokens, int64_t tokenIndex) {
	jsmntok_t subTok = tokens[tokenIndex];
	if(toArraySize != subTok.size) {
		st_errAbort("Reading float array from json but got an unexpected number of tokens. Expected %i but got %i \n",
				(int)toArraySize, (int)subTok.size);
	}
	for (int j = 0; j < subTok.size; j++) {
		toArray[j] = stJson_parseFloat(js, tokens, tokenIndex+j+1);
	}
	return tokenIndex+subTok.size;
}

bool stJson_parseBool(char *js, jsmntok_t *tokens, int64_t tokenIndex) {
	jsmntok_t tok = tokens[tokenIndex];
	char *tokStr = stJson_token_tostr(js, &tok);
	assert(strcmp(tokStr, "true") || strcmp(tokStr, "false"));
	return strcmp(tokStr, "true") == 0;
}

int64_t stJson_getNestedTokenCount(jsmntok_t *tokens, int64_t tokenIndex) {
	/*
	 * Gets the number of tokens for the given element and all the nested tokens it contains.
	 */
	jsmntok_t tok = tokens[tokenIndex];
	int64_t nestedTokens = 1;
	for(int64_t i=0; i<tok.size; i++) {
		nestedTokens += stJson_getNestedTokenCount(tokens, tokenIndex + nestedTokens);
	}
	return nestedTokens;
}

size_t stJson_setupParser(char *buf, size_t r, jsmntok_t **tokens, char **js) {
	// Initialise for parsing
	jsmn_parser parser;
	jsmn_init(&parser);

	// Copy buffer to a mutable array
	*js = st_malloc(r);
	strncpy(*js, buf, r);

	// Now tokenize, first calculating the token number
	int64_t tokenNumber = jsmn_parse(&parser, *js, r, NULL, 0);
	*tokens = st_calloc(tokenNumber, sizeof(jsmntok_t));
	jsmn_init(&parser);
	int64_t i = jsmn_parse(&parser, *js, r, *tokens, tokenNumber);

	if (i == JSMN_ERROR_NOMEM) {
		st_errAbort("Error when parsing json: not enough tokens allocated. Is the JSON file too big? %d\n", i);
	}

	return tokenNumber;
}

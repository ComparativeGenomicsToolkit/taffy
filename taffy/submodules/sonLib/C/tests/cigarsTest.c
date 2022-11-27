/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <string.h>

#include "fastCMaths.h"
#include "avl.h"
#include "commonC.h"
#include "bioioC.h"
#include "pairwiseAlignment.h"

int main(int argc, char *argv[]) {
    int64_t i;
    struct List *pAs;
    FILE *fileHandle;
    struct PairwiseAlignment *pA;
    int64_t keepProbs;

    assert(argc == 3);
    if(strcmp(argv[2], "True") == 0) {
        keepProbs = TRUE;
    }
    else {
        assert(strcmp(argv[2], "False") == 0);
        keepProbs = FALSE;
    }

    pAs = constructEmptyList(0, (void (*)(void *))destructPairwiseAlignment);
    fileHandle = st_fopen(argv[1], "r");
    pA = cigarRead(fileHandle);
    while(pA != NULL) {
        listAppend(pAs, pA);
        pA = cigarRead(fileHandle);
    }
    fclose(fileHandle);

    fileHandle = st_fopen(argv[1], "w");
    for(i=0; i<pAs->length; i++) {
        cigarWrite(fileHandle, pAs->list[i], keepProbs);
    }
    fclose(fileHandle);
    destructList(pAs);

    return 0;
}

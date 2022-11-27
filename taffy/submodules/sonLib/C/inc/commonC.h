/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef COMMONC_H_
#define COMMONC_H_

#include <stdio.h>
#include "fastCMaths.h"
#include "hashTableC.h"
#include "sonLib.h"

#ifdef __cplusplus
extern "C" {
#endif

//utils functions

//logging and debugging
#define DEBUG TRUE

void exitOnFailure(int64_t exitValue, const char *failureMessage, ...);

//memory
struct Chunks {
    struct List *chunkList;
    char * chunk;
    //void * chunk;
    int64_t remaining;
    int64_t chunkSize;
    int64_t elementSize;
};

struct Chunks *constructChunks(int64_t chunkSize, int64_t elementSize);

void destructChunks(struct Chunks *);

void *mallocChunk(struct Chunks *chunk);


//general data structures you always need
//lists
struct List {
    int64_t length;
    int64_t maxLength;
    void **list;
    void (*destructElement)(void *);
};

void listAppend(struct List *list, void *i);

void *listRemoveFirst(struct List *list);

void *arrayResize(void *current, int64_t *currentSize, int64_t newSize, int64_t base);

void listIntersection(struct List *list, struct List *list2, struct List *list3);

void listResize(struct List *list, int64_t newMaxSize);

int64_t listGetInt(struct List *list, int64_t index);

float listGetFloat(struct List *list, int64_t index);

void listReverse(struct List *list);

int64_t listContains(struct List *list, void *k);

void listRemove(struct List *list, void *k);

void listRemoveDuplicates(struct List *list);

int64_t listContainsDuplicates(struct List *list);

void *arrayCopyResize(void *current, int64_t *currentSize, int64_t newSize, int64_t base);

void *arrayPrepareAppend(void *current, int64_t *maxLength, int64_t currentLength, int64_t base);

void arrayShuffle(void **array, int64_t n);

void listCopyResize(struct List *list, int64_t newMaxSize);

struct List *listCopy(struct List *list);

void swapListFields(struct List *list1, struct List *list2);

struct List *cloneList(struct List *source);

void copyList(struct List *from, struct List *to);

struct hashtable *intListToHash(struct List *list, int64_t *(*getKey)(void *));

//list functions
struct List *copyConstructList(void **list, int64_t length, void (*destructElement)(void *));

struct List *constructZeroLengthList(int64_t length, void (*destructElement)(void *));

struct List *constructEmptyList(int64_t length, void (*destructElement)(void *));

void destructList(struct List *list);

void listAppendArray(struct List *list, void **array, int64_t length);

//int lists
struct IntList {
    int64_t length;
    int64_t maxLength;
    int64_t *list;
};

struct IntList *constructEmptyIntList(int64_t length);

void destructIntList(struct IntList *intList);

struct IntList *intListCopy(struct IntList *intList);

void intListAppend(struct IntList *intList, int64_t);

//ints
int64_t *constructInt(int64_t i);

void destructInt(int64_t *i);

int64_t *constructChunkInt(int64_t intValue, struct Chunks *chunks);

int64_t *constructChunkLong(int64_t longValue, struct Chunks *chunks);
//ints
int64_t *constructLong(int64_t i);

void destructLong(int64_t *i);

float *constructFloat(float i);

void destructFloat(float *i);

uint64_t hashtable_stringHashKey( const void *k );

int hashtable_stringEqualKey( const void *key1, const void *key2 );

uint64_t hashtable_intHashKey( const void *k );

int hashtable_intEqualKey( const void *key1, const void *key2 );

uint64_t hashtable_key( const void *k );

int hashtable_equalKey( const void *key1, const void *key2 );

uint64_t hashtable_intPairHashKey( const void *k );

int hashtable_intPairEqualKey( const void *key1, const void *key2 );

uint64_t hashtable_orderedIntPairHashKey( const void *k );

int hashtable_orderedIntPairEqualKey(const  void *key1, const void *key2 );

int64_t *constructIntPair(int64_t i, int64_t j);

void destructIntPair(int64_t *i);

uint64_t hashtable_longHashKey( const void *k );

int hashtable_longEqualKey( const void *key1, const void *key2 );

int64_t intComparator(int64_t *i, int64_t *j);

int64_t longComparator(int64_t *i, int64_t *j);

int intComparator_Int(int64_t *i, int64_t *j);

int longComparator_Int(int64_t *i, int64_t *j);

int64_t intsComparator(int64_t *ints1, int64_t *ints2, int64_t length);

int floatComparator(float **f, float **f2);

struct TraversalID {
    //tree traversal numbers, used as nodeIDs for identifying
    //orders in the tree
    //pre == pre order traversal
    //preEnd == max pre index + 1 of node in subtree
    //mid == mid order (in-order) traversal number
    //def __init__(self, pre, preEnd, mid):
    int64_t midStart;
    int64_t mid;
    int64_t midEnd;
    int64_t leafNo;
};

struct TraversalID *constructTraversalID(int64_t midStart, int64_t mid, int64_t midEnd, int64_t leafNo);

void destructTraversalID(struct TraversalID *traversalID);

struct BinaryTree {
    float distance;
    int64_t internal;
    char *label;
    struct TraversalID *traversalID;
    struct BinaryTree *left;
    struct BinaryTree *right;
};

int64_t leftMostLeafNo(struct TraversalID *traversalID);

int64_t rightMostLeafNo(struct TraversalID *traversalID);

int64_t leafNoInSubtree(struct TraversalID *traversalID);

struct BinaryTree *constructBinaryTree(float distance, int64_t internal,
                                        const char *label,
                                              struct BinaryTree *left,
                                              struct BinaryTree *right);

void destructBinaryTree(struct BinaryTree *binaryTree);

/*
 * Gets all the leaf sequences, in depth first, left to right traversal order.
 */
struct List *binaryTree_getOrderedLeafStrings(struct BinaryTree *binaryTree);

void binaryTree_depthFirstNumbers(struct BinaryTree *binaryTree);

void printBinaryTree(FILE *file, struct BinaryTree *binaryTree);

void annotateTree(struct BinaryTree *bT, void *(*fn)(struct BinaryTree *i), struct List *list);

void getBinaryTreeNodesInMidOrder(struct BinaryTree *binaryTree, struct BinaryTree **labels);

float linOriginRegression(struct List *pointsX, struct List *pointsY);

/*
 * Joins two paths together, somewhat intelligently, to give one concatenated path.
 */
char *pathJoin(const char *pathPrefix, const char *pathSuffix);

/*
 * Returns non zero if the difference of the two values is within the given precision,
 * else returns zero.
 */
int64_t floatValuesClose(double valueOne, double valueTwo, double precision);

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//temp files
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

int64_t constructRandomDir(const char *tempFilePath, char **tempDir);

int64_t destructRandomDir(char *tempDir);

char *getTempFile(void);

void removeTempFile(char *tempFile);

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//graph viz functions
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

void graphViz_addNodeToGraph(const char *nodeName, FILE *graphFileHandle, const char *label,
        double width, double height, const char *shape, const char *colour,
        int64_t fontsize);

void graphViz_addEdgeToGraph(const char *parentNodeName, const char *childNodeName, FILE *graphFileHandle,
        const char *label, const char *colour, double length, double weight, const char *direction);

void graphViz_setupGraphFile(FILE *graphFileHandle);

void graphViz_finishGraphFile(FILE *graphFileHandle);

const char *graphViz_getColour(void);

#ifdef __cplusplus
}
#endif
#endif /*COMMONC_H_*/

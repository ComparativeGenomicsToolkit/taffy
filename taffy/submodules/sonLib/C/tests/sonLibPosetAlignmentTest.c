/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * stPosetAlignmentTest.c
 *
 *  Created on: 1 June 2010
 *      Author: benedictpaten
 */

#include "sonLibGlobalsTest.h"

/*
 * The following test code builds a random alignment and checking as each
 * extra pair of aligned positions is added to the alignment that the alignment
 * remains partially ordered.
 */
static const int64_t MAX_SEQUENCE_NUMBER = 10;
static const int64_t MAX_SEQUENCE_SIZE = 100;
static const int64_t MAX_ALIGNMENTS = 100;
//static const int64_t MAX_ALIGNED_PAIRS = 10000;


static int64_t sequenceNumber;
static stPosetAlignment *posetAlignment = NULL;
static void teardown() {
    if(posetAlignment != NULL) {
        stPosetAlignment_destruct(posetAlignment);
        posetAlignment = NULL;
    }
}

static void setup() {
    teardown();
    sequenceNumber = st_randomInt(0, MAX_SEQUENCE_NUMBER);
    posetAlignment = stPosetAlignment_construct(sequenceNumber);
}


/*
 * Tests the constructor and the basic setup.
 */
static void test_stPosetAlignment_getSequenceNumber(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, sequenceNumber == stPosetAlignment_getSequenceNumber(posetAlignment));
    teardown();
}

/*
 * This builds an adjacency list structure for the the sequences. Every sequence-position
 * has a column in the hash with which it can be aligned with.
 */
static stHash *buildAdjacencyList(stList *pairs, int64_t sequenceNumber) {
    stHash *hash = stHash_construct3((uint64_t (*)(const void *))stIntTuple_hashKey,
            (int (*)(const void *, const void *))stIntTuple_equalsFn,
            (void (*)(void *))stIntTuple_destruct, NULL);
    for(int64_t seq=0; seq<sequenceNumber; seq++) {
        for(int64_t position=0; position<MAX_SEQUENCE_SIZE; position++) {
            stIntTuple *seqPos = stIntTuple_construct2( seq, position);
            stSortedSet *column = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, NULL);
            stSortedSet_insert(column, seqPos);
            stHash_insert(hash, seqPos, column);
        }
    }
    stListIterator *it = stList_getIterator(pairs);
    stIntTuple *pair;
    while((pair = stList_getNext(it)) != NULL) {
       stIntTuple *seqPos1 = stIntTuple_construct2( stIntTuple_get(pair, 0), stIntTuple_get(pair, 1));
       stIntTuple *seqPos2 = stIntTuple_construct2( stIntTuple_get(pair, 2), stIntTuple_get(pair, 3));
       stSortedSet *column1 = stHash_search(hash, seqPos1);
       assert(column1 != NULL);
       stSortedSet *column2 = stHash_search(hash, seqPos2);
       assert(column2 != NULL);
       if(column1 != column2) { //Merge the columns
           stSortedSetIterator *it2 = stSortedSet_getIterator(column2);
           stIntTuple *seqPos3;
           while((seqPos3 = stSortedSet_getNext(it2)) != NULL) {
               assert(stSortedSet_search(column1, seqPos3) == NULL);
               stSortedSet_insert(column1, seqPos3);
               assert(stHash_search(hash, seqPos3) == column2);
               stHash_insert(hash, seqPos3, column1);
               assert(stHash_search(hash, seqPos3) == column1);
           }
           stSortedSet_destructIterator(it2);
           stSortedSet_destruct(column2);
       }
       //Cleanup loop.
       stIntTuple_destruct(seqPos1);
       stIntTuple_destruct(seqPos2);
    }
    stList_destructIterator(it);
    return hash;
}

/*
 * Function does the actual depth first search to detect if the thing has an acyclic ordering.
 */
static int64_t dfs(stHash *adjacencyList, stIntTuple *seqPos,
                               stSortedSet *started, stSortedSet *done) {
    if(stSortedSet_search(started, seqPos) != NULL) {
        if(stSortedSet_search(done, seqPos) == NULL) {
            //We have detected a cycle
            //st_logInfo("I have cycle %" PRIi64 " %" PRIi64 "\n", stIntTuple_getPosition(seqPos, 0), stIntTuple_getPosition(seqPos, 1));
            return 1;
        }
        //We have already explored this area, but no cycle.
        return 0;
    }
    stSortedSet_insert(started, seqPos);

    int64_t cycle =0;

    stIntTuple *nextSeqPos = stIntTuple_construct2( stIntTuple_get(seqPos, 0), stIntTuple_get(seqPos, 1) + 1);
    stSortedSet *column = stHash_search(adjacencyList, nextSeqPos);
    if(column != NULL) { //It is in the adjacency list, so we can do the recursion
        assert(stSortedSet_search(column, nextSeqPos) != NULL);
        stSortedSetIterator *it = stSortedSet_getIterator(column);
        stIntTuple *seqPos2;
        while((seqPos2 = stSortedSet_getNext(it)) != NULL) {
            cycle = cycle || dfs(adjacencyList, seqPos2, started, done);
        }
        stSortedSet_destructIterator(it);
    }
    stIntTuple_destruct(nextSeqPos);
    stSortedSet_insert(done, seqPos);
    return cycle;
}

/*
 * Uses the functions above to build an adjacency list, then by DFS attempts to create
 * a valid topological sort, returning non-zero if the graph contains a cycle.
 */
static int64_t containsACycle(stList *pairs, int64_t sequenceNumber) {
    //Build an adjacency list structure..
    stHash *adjacencyList = buildAdjacencyList(pairs, sequenceNumber);

    //Do a topological sort of the adjacency list
    stSortedSet *started = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, NULL);
    stSortedSet *done = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, NULL);
    int64_t cyclic = 0;
    for(int64_t seq=0; seq<sequenceNumber; seq++) {
        stIntTuple *seqPos = stIntTuple_construct2( seq, 0); //The following hacks avoid memory cleanup..
        stSortedSet *column = stHash_search(adjacencyList, seqPos);
        assert(column != NULL);
        stIntTuple *seqPos2 = stSortedSet_search(column, seqPos);
        assert(seqPos2 != NULL);
        cyclic = cyclic || dfs(adjacencyList, seqPos2, started, done);
        stIntTuple_destruct(seqPos);
    }

    //cleanup
    stHashIterator *it = stHash_getIterator(adjacencyList);
    stIntTuple *seqPos;
    stSortedSet *columns = stSortedSet_construct2((void (*)(void *))stSortedSet_destruct);
    while((seqPos = stHash_getNext(it)) != NULL) {
        stSortedSet *column = stHash_search(adjacencyList, seqPos);
        assert(column != NULL);
        stSortedSet_insert(columns, column);
    }
    stHash_destructIterator(it);
    stHash_destruct(adjacencyList);
    stSortedSet_destruct(columns);
    stSortedSet_destruct(started);
    stSortedSet_destruct(done);

    return cyclic;
}


static void test_stPosetAlignment_addAndIsPossible(CuTest *testCase) {
    for(int64_t trial=0; trial<100; trial++) {
        setup();

        //Make random number of sequences.
        stList *sequenceLengths = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
        for(int64_t i=0; i<sequenceNumber; i++) {
            stList_append(sequenceLengths, stIntTuple_construct1( st_randomInt(0, MAX_SEQUENCE_SIZE)));
        }

        //Propose random alignment pairs...
        stList *pairs = stList_construct3(0, (void(*)(void *))stIntTuple_destruct);
        int64_t maxAlignedPairs = st_randomInt(0, MAX_ALIGNMENTS);
        if(sequenceNumber > 0) {
            for(int64_t i=0; i<maxAlignedPairs; i++) {
                int64_t seq1 = st_randomInt(0, sequenceNumber);
                int64_t seqLength1 = stIntTuple_get(stList_get(sequenceLengths, seq1), 0);
                if(seqLength1 == 0) {
                    continue;
                }
                int64_t position1 = st_randomInt(0, seqLength1);
                int64_t seq2 = st_randomInt(0, sequenceNumber);
                int64_t seqLength2 = stIntTuple_get(stList_get(sequenceLengths, seq1), 0);
                if(seqLength2 == 0) {
                    continue;
                }
                int64_t position2 = st_randomInt(0, seqLength2);
                if(seq1 != seq2) {
                    stList_append(pairs, stIntTuple_construct4( seq1, position1, seq2, position2));
                    if(stPosetAlignment_isPossible(posetAlignment, seq1, position1, seq2, position2)) {
                        st_logInfo("In %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " \n", seq1, position1, seq2, position2);
                        //For each accepted pair check it doesn't create a cycle.
                        CuAssertTrue(testCase, !containsACycle(pairs, sequenceNumber));
                        CuAssertTrue(testCase, stPosetAlignment_add(posetAlignment, seq1, position1, seq2, position2));
                    }
                    else {
                        st_logInfo("Out %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " \n", seq1, position1, seq2, position2);
                        //For each rejected pair check it creates a cycle..
                        CuAssertTrue(testCase, containsACycle(pairs, sequenceNumber));
                        CuAssertTrue(testCase, !stPosetAlignment_isPossible(posetAlignment, seq1, position1, seq2, position2));
                        stIntTuple_destruct(stList_pop(pairs)); //remove the pair which created the cycle.
                        CuAssertTrue(testCase, !containsACycle(pairs, sequenceNumber)); //Check we're back to being okay..
                    }
                }
            }
        }

        //Cleanup
        stList_destruct(sequenceLengths);
        stList_destruct(pairs);
        teardown();
        st_logInfo("Passed a random ordering test with %" PRIi64 " sequences and %" PRIi64 " aligned pairs\n", sequenceNumber, maxAlignedPairs);
    }
}

CuSuite* stPosetAlignmentTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stPosetAlignment_addAndIsPossible);
    SUITE_ADD_TEST(suite, test_stPosetAlignment_getSequenceNumber);
    return suite;
}

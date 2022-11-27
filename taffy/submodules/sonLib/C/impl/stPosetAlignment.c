/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLibGlobalsInternal.h"

struct _stPosetAlignment {
    int64_t sequenceNumber;
    stSortedSet **constraintLists;
};

static int cmpFn(int64_t i, int64_t j) {
    return i > j ? 1 : (i < j ? -1 : 0);
}

static int comparePositions(stIntTuple *position1, stIntTuple *position2) {
    if(stIntTuple_get(position1, 0) == INT64_MAX || stIntTuple_get(position2, 0) == INT64_MAX) { //Indicates we should ignore the first position and compare the second.
        assert(stIntTuple_get(position1, 1) != INT64_MAX);
        assert(stIntTuple_get(position2, 1) != INT64_MAX);
        return cmpFn(stIntTuple_get(position1, 1), stIntTuple_get(position2, 1));
    }
    return cmpFn(stIntTuple_get(position1, 0), stIntTuple_get(position2, 0));
}

stPosetAlignment *stPosetAlignment_construct(int64_t sequenceNumber) {
    stPosetAlignment *posetAlignment = st_malloc(sizeof(stPosetAlignment));
    posetAlignment->sequenceNumber = sequenceNumber;
    posetAlignment->constraintLists = st_malloc(sizeof(stSortedSet *) * sequenceNumber * sequenceNumber);
    for(int64_t i=0; i<posetAlignment->sequenceNumber; i++) {
        for(int64_t j=0; j<posetAlignment->sequenceNumber; j++) {
            if(i != j) {
                posetAlignment->constraintLists[i*posetAlignment->sequenceNumber + j] =
                    stSortedSet_construct3((int (*)(const void *, const void *))comparePositions,
                            (void (*)(void *))stIntTuple_destruct);
            }
        }
    }
    return posetAlignment;
}

void stPosetAlignment_destruct(stPosetAlignment *posetAlignment) {
    for(int64_t i=0; i<posetAlignment->sequenceNumber; i++) {
        for(int64_t j=0; j<posetAlignment->sequenceNumber; j++) {
            if(i != j) {
                stSortedSet_destruct(posetAlignment->constraintLists[i*posetAlignment->sequenceNumber + j]);
            }
        }
    }
    free(posetAlignment->constraintLists);
    free(posetAlignment);
}

int64_t stPosetAlignment_getSequenceNumber(stPosetAlignment *posetAlignment) {
    return posetAlignment->sequenceNumber;
}

static stSortedSet *getConstraintList(stPosetAlignment *posetAlignment, int64_t sequence1, int64_t sequence2) {
    assert(sequence1 >= 0 && sequence1 < posetAlignment->sequenceNumber);
    assert(sequence2 >= 0 && sequence2 < posetAlignment->sequenceNumber);
    assert(sequence1 != sequence2);
    return posetAlignment->constraintLists[sequence1 * posetAlignment->sequenceNumber + sequence2];
}

/*
 * Gets the position in sequence2 that the position in sequence2 must be less or equal to in the alignment
 */
static stIntTuple *getConstraint_lessThan(stPosetAlignment *posetAlignment, int64_t sequence1, int64_t position1, int64_t sequence2) {
    stIntTuple *pos = stIntTuple_construct2(position1, INT64_MAX);
    //Get less than or equal
    stIntTuple *constraint = stSortedSet_searchGreaterThanOrEqual(getConstraintList(posetAlignment, sequence1, sequence2), pos);
    stIntTuple_destruct(pos);
    assert(constraint == NULL || position1 <= stIntTuple_get(constraint, 0));
    return constraint;
}

/*
 * Gets the position in sequence2 that the position in sequence1 must be greater than or equal to in the alignment
 */
static stIntTuple *getConstraint_greaterThan(stPosetAlignment *posetAlignment, int64_t sequence1, int64_t position1, int64_t sequence2) {
    stIntTuple *pos = stIntTuple_construct2(INT64_MAX, position1);
    //Get less than or equal
    stIntTuple *constraint = stSortedSet_searchLessThanOrEqual(getConstraintList(posetAlignment, sequence2, sequence1), pos);
    stIntTuple_destruct(pos);
    assert(constraint == NULL || position1 >= stIntTuple_get(constraint, 1));
    return constraint;
}

/*
 * Returns non-zero iff the constraint is prime. The lessThanOrEquals argument, if non-zero specifies the constraint is less than equals.
 */
static bool lessThanConstraintIsPrime(stPosetAlignment *posetAlignment, int64_t sequence1, int64_t position1, int64_t sequence2, int64_t position2, int64_t lessThanOrEquals) {
    stIntTuple *constraint = getConstraint_lessThan(posetAlignment, sequence1, position1, sequence2);
    if(constraint == NULL) {
        return 1;
    }
    if(position2 < stIntTuple_get(constraint, 1)) { //new constraint is tighter
        return 1;
    }
    if(position2 > stIntTuple_get(constraint, 1)) { //new constraint is looser
        return 0;
    }
    if(position1 == stIntTuple_get(constraint, 0) && stIntTuple_get(constraint, 2) && !lessThanOrEquals) { //converts a less than or equals constraint to a less than constraint
        return 1;
    }
    return 0;
}

/*
 * Adds a prime less than (or equals) constraint to the list of prime constraints, removing any redundant constraints in the process.
 * The or equals is specified by making the lessThanOrEquals argument non-zero.
 */
void addConstraint_lessThan(stPosetAlignment *posetAlignment, int64_t sequence1, int64_t position1, int64_t sequence2, int64_t position2, int64_t lessThanOrEquals) {
    stSortedSet *constraintList = getConstraintList(posetAlignment, sequence1, sequence2);
    assert(position1 != INT64_MAX);
    assert(position2 != INT64_MAX);
    stIntTuple *constraint1 = stIntTuple_construct3( position1, position2, lessThanOrEquals);
    stIntTuple *constraint2;
    while((constraint2 = stSortedSet_searchLessThanOrEqual(constraintList, constraint1)) != NULL) {
        assert(stIntTuple_get(constraint2, 0) <= position1);
        if(stIntTuple_get(constraint2, 1) >= position2) {
            if(stIntTuple_get(constraint2, 1) == position2) { //Check we are not removing an equivalent or more severe constraint.
                assert((!lessThanOrEquals && stIntTuple_get(constraint2, 2)) || stIntTuple_get(constraint2, 0) < position1);
            }
            stSortedSet_remove(constraintList, constraint2);
            stIntTuple_destruct(constraint2);
        }
        else {
            assert(stIntTuple_get(constraint2, 0) < position1); //Check the constraint does not overshadow our proposed constraint.
            break;
        }
    }
    stSortedSet_insert(constraintList, constraint1);
}

bool stPosetAlignment_isPossibleP(stPosetAlignment *posetAlignment, int64_t sequence1, int64_t position1, int64_t sequence2, int64_t position2) {
    stIntTuple *constraint = getConstraint_lessThan(posetAlignment, sequence1, position1, sequence2);
    if(constraint == NULL) {
        return 1;
    }
    if(stIntTuple_get(constraint, 2) && stIntTuple_get(constraint, 0) == position1) { //less than or equals
        return position2 <= stIntTuple_get(constraint, 1);
    }
    else {
        return position2 < stIntTuple_get(constraint, 1);
    }
}

bool stPosetAlignment_isPossible(stPosetAlignment *posetAlignment, int64_t sequence1, int64_t position1, int64_t sequence2, int64_t position2) {
    return stPosetAlignment_isPossibleP(posetAlignment, sequence1, position1, sequence2, position2) &&
    stPosetAlignment_isPossibleP(posetAlignment, sequence2, position2, sequence1, position1);
}

static void stPosetAlignment_addP2(stPosetAlignment *posetAlignment, int64_t sequence1, int64_t sequence3, int64_t position3, int64_t sequence2, int64_t position2, int64_t lessThanOrEqual) {
    for(int64_t sequence4=0; sequence4<posetAlignment->sequenceNumber; sequence4++) {
        if(sequence4 != sequence1 && sequence4 != sequence2 && sequence4 != sequence3) {
            stIntTuple *constraint = getConstraint_lessThan(posetAlignment, sequence2, position2, sequence4);
            if(constraint != NULL) {
                int64_t position4 = stIntTuple_get(constraint, 1);
                int64_t transLessThanOrEqual = lessThanOrEqual && stIntTuple_get(constraint, 2) && stIntTuple_get(constraint, 0) == position2; //stuff which maintains the less than or equals
                if(lessThanConstraintIsPrime(posetAlignment, sequence3, position3, sequence4, position4, transLessThanOrEqual)) {//We have a new transitive constraint..
                    addConstraint_lessThan(posetAlignment, sequence3, position3, sequence4, position4, transLessThanOrEqual);
                }
            }
        }
    }
}

static void stPosetAlignment_addP(stPosetAlignment *posetAlignment, int64_t sequence1, int64_t position1, int64_t sequence2, int64_t position2) {
    //for all pairs do check..
    if(lessThanConstraintIsPrime(posetAlignment, sequence1, position1, sequence2, position2, 1)) {
        addConstraint_lessThan(posetAlignment, sequence1, position1, sequence2, position2, 1);
        for(int64_t sequence3=0; sequence3<posetAlignment->sequenceNumber; sequence3++) {
            if(sequence3 != sequence2) {
                if(sequence3 != sequence1) {
                    stIntTuple *constraint = getConstraint_greaterThan(posetAlignment, sequence1, position1, sequence3);
                    if(constraint != NULL) {
                        int64_t position3 = stIntTuple_get(constraint, 0); //its reversed
                        int64_t lessThanOrEqual = stIntTuple_get(constraint, 2) && stIntTuple_get(constraint, 1) == position1;
                        if(lessThanConstraintIsPrime(posetAlignment, sequence3, position3, sequence2, position2, lessThanOrEqual)) { //new constraint found, so add it to the set..
                            addConstraint_lessThan(posetAlignment, sequence3, position3, sequence2, position2, lessThanOrEqual);
                            stPosetAlignment_addP2(posetAlignment, sequence1, sequence3, position3, sequence2, position2, lessThanOrEqual);
                        }
                    }
                }
                else {
                    stPosetAlignment_addP2(posetAlignment, INT64_MAX, sequence1, position1, sequence2, position2, 1);
                }
            }
        }
    }
}

bool stPosetAlignment_add(stPosetAlignment *posetAlignment, int64_t sequence1, int64_t position1, int64_t sequence2, int64_t position2) {
    if(stPosetAlignment_isPossible(posetAlignment, sequence1, position1, sequence2, position2)) {
        stPosetAlignment_addP(posetAlignment, sequence1, position1, sequence2, position2);
        stPosetAlignment_addP(posetAlignment, sequence2, position2, sequence1, position1);
        return 1;
    }
    return 0;
}

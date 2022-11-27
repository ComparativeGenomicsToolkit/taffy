/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * basic tests of stExcept.
 */
#include "sonLibExcept.h"
#include "sonLibGlobalsTest.h"
#include "sonLibCommon.h"
#include <stdlib.h>
#include <stdio.h>
#include "stSafeC.h"

/* test throwing through two levels */
static const char *const ERR1 = "err1";
static const char *const ERR2 = "err2";
static void thrower2() {
    stThrowNew(ERR2, "error in %s", "thrower2");
}

static void thrower1(void) {
    stTry {
        thrower2();
    } stCatch(except) {
        stThrow(stExcept_newCause(except, ERR1, "error in %s", "thrower1"));
    } stTryEnd;
}

static void testThrow(CuTest *testCase) {
    volatile int gotEx = 0;
    volatile int pastThrow = 0;
    stTry {
        thrower1();
        pastThrow++;
    } stCatch(except) {
        CuAssertTrue(testCase, stExcept_getId(except) == ERR1);
        CuAssertTrue(testCase, stExcept_idEq(except, ERR1));
        CuAssertStrEquals(testCase, stExcept_getMsg(except), "error in thrower1");

        const stExcept *cause = stExcept_getCause(except); 
        CuAssertTrue(testCase, cause != NULL);
        if (cause != NULL) {
            CuAssertTrue(testCase, stExcept_getId(cause) == ERR2);
            CuAssertStrEquals(testCase, stExcept_getMsg(cause), "error in thrower2");
            CuAssertTrue(testCase, stExcept_getCause(cause) == NULL);
        }
        gotEx++;
    } stTryEnd;
    CuAssertTrue(testCase, pastThrow == 0);
    CuAssertTrue(testCase,  gotEx != 0);
    CuAssertTrue(testCase, _cexceptTOS == NULL);
}

/* test execution with no exception */
static void noop(void) {
}

static void testOk(CuTest *testCase) {
    volatile int gotEx = 0;
    volatile int pastOk = 0;
    stTry {
        noop();
        pastOk++;
    } stCatch(except) {
        stExcept_free(except);
        gotEx++;
    } stTryEnd;
    CuAssertTrue(testCase, pastOk != 0);
    CuAssertTrue(testCase,  gotEx == 0);
    CuAssertTrue(testCase,  _cexceptTOS == NULL);
}

/* test ceTryReturn */
static int returnFromTry(void) {
    stTry {
        stTryReturn(10);
    } stCatch(except) {
        stExcept_free(except);
        return 11;
    } stTryEnd;
    return 12;
}

static int returnFromCatch(void) {
    stTry {
        stThrowNew(ERR1, "throw from catch");
    } stCatch(except) {
        stExcept_free(except);
        return 11;
    } stTryEnd;
    return 12;
}

static int returnAtEnd(void) {
    stTry {
    } stCatch(except) {
        stExcept_free(except);
        return 11;
    } stTryEnd;
    return 12;
}

static void testTryReturn(CuTest *testCase) {
    int val = returnFromTry();
    CuAssertTrue(testCase, val == 10);
    val = returnFromCatch();
    CuAssertTrue(testCase, val == 11);
    val = returnAtEnd();
    CuAssertTrue(testCase, val == 12);
}

#if 0
// FIXME: finish this once there are some functions to read in all of a file

/* Child process that throws an exception with no catch, redirecting stdout/stderr to
 * the specified file. */
static void noCatchChild(const char *errFile) {
    fflush(stdout);
    fflush(stderr);
    int fn = open(errFile, O_WRONLY|O_CREAT|O_TRUNC, 0666);
    if (fn < 0) {
        st_errnoAbort("can't open: %s", errFile);
    }
    if (dup2(STDOUT_FILENO, fn) < 0) {
        st_errnoAbort("dup stdout failed");
    }
    if (dup2(STDERR_FILENO, fn) < 0) {
        st_errnoAbort("dup stderr failed");
    }
    thrower1();
    fprintf(stderr, "shouldn't make it past throw without catch\n");
    exit(1);
}

/* test handling of a throw without a catch, which is must be run in a
 * separate process */
static void testNoCatch(CuTest *testCase) {
    // FIXME: need function to get tmp file honoring TMPENV
    char *errFile = "sonLibExceptTest.tmp";
    noCatchChild
    
    unlink(errFile);  // ignore failures
}
#endif

CuSuite* sonLib_stExceptTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testThrow);
    SUITE_ADD_TEST(suite, testOk);
    SUITE_ADD_TEST(suite, testTryReturn);
    return suite;
}


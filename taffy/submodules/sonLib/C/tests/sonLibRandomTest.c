/*
 * Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibRandomTest.c
 *
 *  Created on: 22-Jun-2010
 *      Author: benedictpaten
 */

#include "sonLibGlobalsTest.h"

static int cmp32(const void *a, const void *b) {
    const int64_t *ua = (const int64_t *)a;
    const int64_t *ub = (const int64_t *)b;
    if (*ua < *ub) {
        return -1;
    } else {
        return 1;
    }
    return 0;
}
static int cmp64(const void *a, const void *b) {
    const int64_t *ua = (const int64_t *)a;
    const int64_t *ub = (const int64_t *)b;
    if (*ua < *ub) {
        return -1;
    } else {
        return 1;
    }
    return 0;
}
static double median32(int64_t *array, uint64_t n) {
    qsort(array, n, sizeof(int64_t), cmp32);
    if (n % 2) {
        return (double) array[(n + 1) / 2 - 1];
    } else {
        return (double) ((array[n / 2 - 1] + array[(n / 2)]) / 2.0);
    }
}
static double median64(int64_t *array, uint64_t n) {
    qsort(array, n, sizeof(int64_t), cmp64);
    if (n % 2) {
        return (double) array[(n + 1) / 2 - 1];
    } else {
        return (double) ((array[n / 2 - 1] + array[(n / 2)]) / 2.0);
    }
}
static double runningMean32(double m, int64_t x, uint64_t i) {
    m += (x - m) / (i + 1);
    return m;
}
static double runningMean64(double m, int64_t x, uint64_t i) {
    m += (x - m) / (i + 1);
    return m;
}

static double sv32(int64_t *array, uint64_t n, double mu) {
    // sample variance.
    double sv = 0.0;
    assert(n > 1);
    for (uint64_t i = 0; i < n; ++i) {
        sv += (array[i] - mu) * (array[i] - mu);
    }
    sv /= (n - 1);
    return sv;
}
static double sv64(int64_t *array, uint64_t n, double mu) {
    // sample variance.
    double sv = 0.0;
    assert(n > 1);
    for (uint64_t i = 0; i < n; ++i) {
        sv += (array[i] - mu) * (array[i] - mu);
    }
    sv /= (n - 1);
    return sv;
}
static void test_st_randomInt(CuTest *testCase) {
    /*
     * Exercises the random int function.
     */
	int64_t min = INT64_MIN;
	int64_t max = INT64_MAX;
    int64_t v;
    for (int64_t i = -9; i < 10; ++i) {
        CuAssertTrue(testCase, st_randomInt(i, i + 1) == i);
    }
    for(int64_t i = 0; i < 100000; i++) {
        v = st_randomInt(min, max);
        CuAssertTrue(testCase, v >= min);
        CuAssertTrue(testCase, v < max);
    }
    min = -1;
    max = 10;
    for(int64_t i = 0; i < 100000; i++) {
        v = st_randomInt(min, max);
        CuAssertTrue(testCase, v >= min);
        CuAssertTrue(testCase, v < max);
    }
    stTry {
        st_randomInt(1, 1);
    } stCatch(except) {
        CuAssertTrue(testCase, stExcept_getId(except) == RANDOM_EXCEPTION_ID);
    } stTryEnd
}
static void test_st_randomInt_distribution_0(CuTest *testCase) {
    /*
     * check the distribution of the randomInt function
     */
    double med, var, mu = 0.0;
    uint64_t n = 1000000;
    uint64_t reps = 50;
    int64_t min = 0;
    int64_t max = 101;
    int64_t *array = NULL;
    for (uint64_t j = 0; j < reps; ++j) {
        array = (int64_t *) st_malloc(sizeof(*array) * n);
        for (int64_t i = 0; i < n; i++) {
            array[i] = st_randomInt(min, max);
            mu = runningMean32(mu, array[i], i);
        }
        var = sv32(array, n, mu);
        med = median32(array, n); // median() will sort array
        CuAssertTrue(testCase, array[0] <= array[n - 1]);
        CuAssertTrue(testCase, array[0] < max);
        CuAssertTrue(testCase, array[n - 1] < max);
        // median of discrete uniform is
        // median = (a + b) / 2
        CuAssertTrue(testCase, med == 50.0);
        // mean of discrete uniform is
        // mu = (a + b) / 2
        CuAssertTrue(testCase, mu > 0.9 * (max - 1) / 2.0);
        CuAssertTrue(testCase, mu < 1.1 * (max - 1) / 2.0);
        // variance of discrete uniform is
        // var = ((b - a + 1)^2 - 1) / 12
        CuAssertTrue(testCase, var < 1.1 * (((max - min + 1) * (max - min + 1) - 1) / 12.0));
        CuAssertTrue(testCase, var > 0.9 * (((max - min + 1) * (max - min + 1) - 1) / 12.0));
        free(array);
    }
}
static void test_st_randomInt64_range_0(CuTest *testCase) {
    /*
     * Excercies the randomInt64 function.
     */
	int64_t min = INT64_MIN;
	int64_t max = INT64_MAX;
    int64_t v;
    for (int64_t i = -9; i < 10; ++i) {
        CuAssertTrue(testCase, st_randomInt64(i, i + 1) == i);
    }
    for (int64_t i = 0; i < 100000; i++) {
        v = st_randomInt64(min, max);
        CuAssertTrue(testCase, v >= min);
        CuAssertTrue(testCase, v < max);
    }
    stTry {
        st_randomInt64(1, 1);
    } stCatch(except) {
        CuAssertTrue(testCase, stExcept_getId(except) == RANDOM_EXCEPTION_ID);
    } stTryEnd
}
static void test_st_randomInt64_distribution_0(CuTest *testCase) {
    /*
     * check the distribution of the randomInt64 function
     */
    double med, var, mu = 0.0;
    uint64_t n = 1000000;
    uint64_t reps = 50;
    int64_t min = 0;
    int64_t max = 101;
    int64_t *array = NULL;
    for (uint64_t j = 0; j < reps; ++j) {
        array = (int64_t *) st_malloc(sizeof(*array) * n);
        for (int64_t i = 0; i < n; i++) {
            array[i] = st_randomInt64(min, max);
            mu = runningMean64(mu, array[i], i);
        }
        var = sv64(array, n, mu);
        med = median64(array, n); // median() will sort array
        CuAssertTrue(testCase, array[0] <= array[n - 1]);
        CuAssertTrue(testCase, array[0] < max);
        CuAssertTrue(testCase, array[n - 1] < max);
        // median of discrete uniform is
        // median = (a + b) / 2
        CuAssertTrue(testCase, med == 50.0);
        // mean of discrete uniform is
        // mu = (a + b) / 2
        CuAssertTrue(testCase, mu > 0.9 * (max - 1) / 2.0);
        CuAssertTrue(testCase, mu < 1.1 * (max - 1) / 2.0);
        // variance of discrete uniform is
        // var = ((b - a + 1)^2 - 1) / 12
        CuAssertTrue(testCase, var < 1.1 * (((max - min + 1) * (max - min + 1) - 1) / 12.0));
        CuAssertTrue(testCase, var > 0.9 * (((max - min + 1) * (max - min + 1) - 1) / 12.0));
        free(array);
    }
}
static void test_st_random(CuTest *testCase) {
    /*
     * Excercies the random int function.
     */
    for(int64_t i = 0; i < 10000; i++) {
        CuAssertTrue(testCase, st_random() >= 0);
        CuAssertTrue(testCase, st_random() < 1.0);
    }
}
static void test_st_randomChoice(CuTest *testCase) {
    /*
     * Excercies the random int function.
     */
    stList *list = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);

    stTry {
        st_randomChoice(list);
    } stCatch(except) {
        CuAssertTrue(testCase, stExcept_getId(except) == RANDOM_EXCEPTION_ID);
    } stTryEnd

    for(int64_t i = 0; i < 10; i++) {
        stList_append(list, stIntTuple_construct1( i));
    }
    for(int64_t i = 0; i < 100; i++) {
        CuAssertTrue(testCase, stList_contains(list, st_randomChoice(list)));
    }
    stList_destruct(list);
}
CuSuite* sonLib_stRandomTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_st_randomInt);
    SUITE_ADD_TEST(suite, test_st_randomInt_distribution_0);
    SUITE_ADD_TEST(suite, test_st_randomInt64_range_0);
    SUITE_ADD_TEST(suite, test_st_randomInt64_distribution_0);
    SUITE_ADD_TEST(suite, test_st_random);
    SUITE_ADD_TEST(suite, test_st_randomChoice);
    return suite;
}

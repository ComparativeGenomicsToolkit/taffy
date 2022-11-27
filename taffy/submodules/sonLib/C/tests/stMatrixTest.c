/*
 * Copyright (C) 2006-2014 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibMatrixTest.c
 *
 *  Created on: 21-May-2014
 *      Author: benedictpaten
 */

#include "sonLibGlobalsTest.h"

void test_stMatrixBasics(CuTest *testCase) {
    for (int64_t test = 0; test < 100; test++) {
        int64_t n = st_randomInt64(0, 10);
        int64_t m = st_randomInt64(0, 10);
        stMatrix *matrix = stMatrix_construct(n, m);
        CuAssertIntEquals(testCase, n, stMatrix_n(matrix));
        CuAssertIntEquals(testCase, m, stMatrix_m(matrix));
        for (int64_t i = 0; i < n; i++) {
            for (int64_t j = 0; j < m; j++) {
                CuAssertDblEquals(testCase, 0.0, *stMatrix_getCell(matrix, i, j), 0.0);
                *stMatrix_getCell(matrix, i, j) = i * j;
                CuAssertDblEquals(testCase, i * j, *stMatrix_getCell(matrix, i, j), 0.0);
            }
        }
        stMatrix_destruct(matrix);
    }
}

static stMatrix *getRandomMatrix(int64_t n, int64_t m) {
    stMatrix *matrix = stMatrix_construct(n, m);
    for (int64_t i = 0; i < n; i++) {
        for (int64_t j = 0; j < m; j++) {
            *stMatrix_getCell(matrix, i, j) = st_random();
        }
    }
    return matrix;
}

void test_stMatrixAdd(CuTest *testCase) {
    for (int64_t test = 0; test < 100; test++) {
        int64_t n = st_randomInt64(0, 10);
        int64_t m = st_randomInt64(0, 10);
        stMatrix *matrix1 = getRandomMatrix(n, m);
        stMatrix *matrix2 = getRandomMatrix(n, m);
        stMatrix *matrix3 = stMatrix_add(matrix1, matrix2);
        for (int64_t i = 0; i < n; i++) {
            for (int64_t j = 0; j < m; j++) {
                CuAssertDblEquals(testCase, *stMatrix_getCell(matrix3, i, j),
                        *stMatrix_getCell(matrix1, i, j) + *stMatrix_getCell(matrix2, i, j), 0.0);
            }
        }
        stMatrix_destruct(matrix1);
        stMatrix_destruct(matrix2);
        stMatrix_destruct(matrix3);
    }
}

void test_stMatrixScale(CuTest *testCase) {
    for (int64_t test = 0; test < 100; test++) {
        int64_t n = st_randomInt64(0, 10);
        int64_t m = st_randomInt64(0, 10);
        stMatrix *matrix = getRandomMatrix(n, m);
        stMatrix *matrix2 = stMatrix_clone(matrix);
        stMatrix_scale(matrix, 100, 50);
        for (int64_t i = 0; i < n; i++) {
            for (int64_t j = 0; j < m; j++) {
                CuAssertDblEquals(testCase, (*stMatrix_getCell(matrix2, i, j)) * 100 + 50, *stMatrix_getCell(matrix, i, j),
                        0.0);
            }
        }
        stMatrix_destruct(matrix);
        stMatrix_destruct(matrix2);
    }
}

void test_stMatrixEqual(CuTest *testCase) {
    for (int64_t test = 0; test < 100; test++) {
        int64_t n = st_randomInt64(1, 10);
        int64_t m = st_randomInt64(1, 10);
        stMatrix *matrix1 = getRandomMatrix(n, m);
        stMatrix *matrix2 = stMatrix_clone(matrix1);
        CuAssertTrue(testCase, stMatrix_equal(matrix1, matrix2, 0.0));
        stMatrix_scale(matrix1, 1.0, 1.0);
        CuAssertTrue(testCase, !stMatrix_equal(matrix1, matrix2, 0.0));
        CuAssertTrue(testCase, stMatrix_equal(matrix1, matrix2, 1.0));
        CuAssertTrue(testCase, !stMatrix_equal(matrix1, matrix2, 0.99));
        stMatrix_destruct(matrix1);
        stMatrix_destruct(matrix2);
    }
}

void test_stMatrixMultiply(CuTest *testCase) {
    stMatrix *m1 = stMatrix_construct(2, 3);
    stMatrix *m2 = stMatrix_construct(3, 1);
    *stMatrix_getCell(m1, 0, 0) = 1.0;
    *stMatrix_getCell(m1, 0, 1) = 2.0;
    *stMatrix_getCell(m1, 0, 2) = 3.0;
    *stMatrix_getCell(m1, 1, 0) = 4.0;
    *stMatrix_getCell(m1, 1, 1) = 5.0;
    *stMatrix_getCell(m1, 1, 2) = 6.0;
    *stMatrix_getCell(m2, 0, 0) = 7.0;
    *stMatrix_getCell(m2, 1, 0) = 8.0;
    *stMatrix_getCell(m2, 2, 0) = 9.0;
    stMatrix *m3 = stMatrix_multiply(m1, m2);
    CuAssertTrue(testCase, stMatrix_n(m3) == 2);
    CuAssertTrue(testCase, stMatrix_m(m3) == 1);
    CuAssertDblEquals(testCase, *stMatrix_getCell(m3, 0, 0), 50.0, 0.00001);
    CuAssertDblEquals(testCase, *stMatrix_getCell(m3, 1, 0), 122.0, 0.00001);
    stMatrix_destruct(m1);
    stMatrix_destruct(m2);
    stMatrix_destruct(m3);
}

void test_stMatrixMultiplyVector(CuTest *testCase) {
    stMatrix *m = stMatrix_construct(2, 2);
    double *v = st_calloc(2, sizeof(double));
    *stMatrix_getCell(m, 0, 0) = 2.0;
    *stMatrix_getCell(m, 0, 1) = 1.0;
    *stMatrix_getCell(m, 1, 0) = 3.0;
    *stMatrix_getCell(m, 1, 1) = 5.0;
    v[0] = 2.0;
    v[1] = 3.0;
    double *v2 = stMatrix_multiplySquareMatrixAndColumnVector(m, v);
    CuAssertDblEquals(testCase, v2[0], 7.0, 0.00001);
    CuAssertDblEquals(testCase, v2[1], 21.0, 0.00001);
    stMatrix_destruct(m);
    free(v);
    free(v2);
}

void test_stMatrixJukesCantor(CuTest *testCase) {
    stMatrix *jukesCantorMatrix = stMatrix_jukesCantor(0.5, 2);
    CuAssertTrue(testCase, stMatrix_n(jukesCantorMatrix) == 2);
    CuAssertTrue(testCase, stMatrix_m(jukesCantorMatrix) == 2);
    CuAssertDblEquals(testCase, *stMatrix_getCell(jukesCantorMatrix, 0, 0), *stMatrix_getCell(jukesCantorMatrix, 1, 1), 0.0);
    CuAssertDblEquals(testCase, *stMatrix_getCell(jukesCantorMatrix, 0, 1), *stMatrix_getCell(jukesCantorMatrix, 1, 0), 0.0);

    CuAssertDblEquals(testCase, *stMatrix_getCell(jukesCantorMatrix, 0, 0), 0.5 + 0.5*exp(-1), 0.0);
    CuAssertDblEquals(testCase, *stMatrix_getCell(jukesCantorMatrix, 0, 1), 0.5 - 0.5*exp(-1), 0.0);
    stMatrix_destruct(jukesCantorMatrix);
}

CuSuite* sonLib_stMatrixTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stMatrixBasics);
    SUITE_ADD_TEST(suite, test_stMatrixAdd);
    SUITE_ADD_TEST(suite, test_stMatrixScale);
    SUITE_ADD_TEST(suite, test_stMatrixEqual);
    SUITE_ADD_TEST(suite, test_stMatrixMultiply);
    SUITE_ADD_TEST(suite, test_stMatrixMultiplyVector);
    SUITE_ADD_TEST(suite, test_stMatrixJukesCantor);

    return suite;
}

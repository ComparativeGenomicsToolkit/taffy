/*
 * Copyright (C) 2006-2014 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * stMatrix.c
 *
 *  Created on: 21-May-2014
 *      Author: benedictpaten
 */

#include "sonLibGlobalsInternal.h"

/*
 * Represents a matrix.
 */
struct _stMatrix {
    int64_t n; //Matrix is n x m.
    int64_t m; //Matrix is n x m.
    double *M; //Matrix.
};

stMatrix *stMatrix_construct(int64_t n, int64_t m) {
    stMatrix *matrix = st_malloc(sizeof(stMatrix));
    matrix->n = n;
    matrix->m = m;
    matrix->M = st_calloc(n * m, sizeof(double));
    return matrix;
}

int64_t stMatrix_n(stMatrix *matrix) {
    return matrix->n;
}

int64_t stMatrix_m(stMatrix *matrix) {
    return matrix->m;
}

void stMatrix_destruct(stMatrix *matrix) {
    free(matrix->M);
    free(matrix);
}

double *stMatrix_getCell(stMatrix *matrix, int64_t indexN, int64_t indexM) {
    assert(indexN >= 0 && indexN < matrix->n);
    assert(indexM >= 0 && indexM < matrix->m);
    return &(matrix->M[indexN * matrix->m + indexM]);
}

stMatrix *stMatrix_multiply(stMatrix *matrix1, stMatrix *matrix2) {
    if(stMatrix_m(matrix1) != stMatrix_n(matrix2)) {
        stThrow(stExcept_new("MATRIX_EXCEPTION", "Matrices do not have equal length dimensions (%" PRIi64  "%" PRIi64 ") to multiply", stMatrix_m(matrix1), stMatrix_n(matrix2)));
    }
    stMatrix *matrix3 = stMatrix_construct(stMatrix_n(matrix1), stMatrix_m(matrix2));
    for(int64_t i=0; i<stMatrix_n(matrix1); i++) {
        for(int64_t j=0; j<stMatrix_m(matrix2); j++) {
            double *cell = stMatrix_getCell(matrix3, i, j);
            for(int64_t k=0; k<stMatrix_m(matrix1); k++) {
                *cell += *stMatrix_getCell(matrix1, i, k) * *stMatrix_getCell(matrix2, k, j);
            }
        }
    }
    return matrix3;
}

void stMatrix_multiplySquareMatrixAndColumnVector2(stMatrix *matrix, double *vector, double *output) {
    int64_t m = stMatrix_m(matrix), n = stMatrix_n(matrix);
    double *M = matrix->M;
    if(m != n) {
        stThrow(stExcept_new("MATRIX_EXCEPTION", "Matrix is not a square matrix (%" PRIi64  "%" PRIi64 ") to multiply", stMatrix_m(matrix), stMatrix_n(matrix)));
    }
    for(int64_t i=0; i<n; i++) {
        output[i] = 0.0;
    }
    for(int64_t i=0; i<n; i++) {
        for(int64_t j=0; j<m; j++) {
            output[i] += M[i * m + j] * vector[j];
        }
    }
}

double *stMatrix_multiplySquareMatrixAndColumnVector(stMatrix *matrix, double *vector) {
    double *output = st_malloc(stMatrix_n(matrix) * sizeof(double));
    stMatrix_multiplySquareMatrixAndColumnVector2(matrix, vector, output);
    return output;
}

stMatrix *stMatrix_add(stMatrix *matrix1, stMatrix *matrix2) {
    assert(matrix1->n == matrix2->n);
    assert(matrix1->m == matrix2->m);
    stMatrix *mergedMatrix = stMatrix_construct(matrix1->n, matrix1->m);
    for (int64_t i = 0; i < matrix1->n * matrix1->m; i++) {
        mergedMatrix->M[i] = matrix1->M[i] + matrix2->M[i];
    }
    return mergedMatrix;
}

stMatrix *stMatrix_clone(stMatrix *matrix) {
    stMatrix *matrix2 = stMatrix_construct(matrix->n, matrix->m);
    for (int64_t i = 0; i < matrix->n * matrix->m; i++) {
        matrix2->M[i] = matrix->M[i];
    }
    return matrix2;
}

void stMatrix_scale(stMatrix *matrix, double scaleFactor, double addFactor) {
    for (int64_t i = 0; i < matrix->n * matrix->m; i++) {
        matrix->M[i] = matrix->M[i] * scaleFactor + addFactor;
    }
}

bool stMatrix_equal(stMatrix *matrix1, stMatrix *matrix2, double close) {
    if(matrix1->n != matrix2->n || matrix1->m != matrix2->m) {
        return 0;
    }
    for(int64_t i=0; i<matrix1->n * matrix1->m; i++) {
        if(matrix1->M[i] + close < matrix2->M[i] || matrix1->M[i] - close > matrix2->M[i]) {
            return 0;
        }
    }
    return 1;
}

stMatrix *stMatrix_jukesCantor(double distance, int64_t n) {
    stMatrix *jkMatrix = stMatrix_construct(n, n);
    for(int64_t i=0; i<n; i++) {
        for(int64_t j=0; j<n; j++) {
            *stMatrix_getCell(jkMatrix, i, j) = i == j ? 1.0/n + ((n-1.0)/n) * exp(-n*distance/(n-1.0)) : 1.0/n - (1.0/n) * exp(-n*distance/(n-1.0)); //Jukes-Cantor calculation
        }
    }
    return jkMatrix;
}


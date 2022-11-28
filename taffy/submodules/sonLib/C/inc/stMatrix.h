/*
 * Copyright (C) 2009-2014 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * stMatrix.h Really toy matrix code.
 *
 *  Created on: 21-May-2014
 *      Author: benedictpaten
 */

#ifndef STMATRIX_H_
#define STMATRIX_H_

#include "sonLibTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Returns an n x n double precision matrix.
 */
stMatrix *stMatrix_construct(int64_t n, int64_t m);

void stMatrix_destruct(stMatrix *matrix);

int64_t stMatrix_n(stMatrix *matrix);

int64_t stMatrix_m(stMatrix *matrix);

/*
 * Gets the cell for the given indices.
 */
double *stMatrix_getCell(stMatrix *matrix, int64_t indexN, int64_t indexM);

/*
 * Adds together two matrices.
 */
stMatrix *stMatrix_add(stMatrix *matrix1, stMatrix *matrix2);

/*
 * Multiples two matrices, one of which is i x j and the other is j x k producing and i x k matrix.
 */
stMatrix *stMatrix_multiply(stMatrix *matrix1, stMatrix *matrix2);

/*
 *  Multiples a an n x n square matrix with a n length column vector to produce a n length output vector.
 */
double *stMatrix_multiplySquareMatrixAndColumnVector(stMatrix *matrix1, double *vector);

/*
 *  As stMatrix_multiplySquareMatrixAndColumnVector but stores the result in output.
 */
void stMatrix_multiplySquareMatrixAndColumnVector2(stMatrix *matrix, double *vector, double *output);

/*
 * Clone the matrix.
 */
stMatrix *stMatrix_clone(stMatrix *matrix);

/*
 * Scales the values in the given matrix by the given coefficients.
 */
void stMatrix_scale(stMatrix *matrix, double scaleFactor, double addFactor);

/*
 * Returns true iff the matrices have the same dimensions and are all values are within close+/- of each other.
 */
bool stMatrix_equal(stMatrix *matrix1, stMatrix *matrix2, double close);

/*
 * Return a n*n matrix representing a Jukes-Cantor substitution matrix for n characters.
 */
stMatrix *stMatrix_jukesCantor(double distance, int64_t n);

#ifdef __cplusplus
}
#endif
#endif /* STMATRIX_H_ */

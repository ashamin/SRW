#ifndef SOLVERS_H
#define SOLVERS_H

#include "Headers/srwvector.h"
#include "Headers/srwmatrix.h"
#include "Headers/eigen3vector.h"

/**
 * @brief TDMA_d
 *          Tridiagonal Matrix Algorithm
 * @param a
 *          lower diagonal
 * @param b
 *          main diagonal
 * @param c
 *          upper diagonal
 * @param d
 *          right part equation vector
 * @param x
 *          result vector x from Ax = d
 * @return
 *          vector x from Ax = d
 */
SRWVector& TDMA_d(SRWVector& a, SRWVector& b,
                  SRWVector& c, SRWVector& d, SRWVector& x);

SRWVector& TDMA(SRWMatrix& A, SRWVector&x, SRWVector& b);

SRWVector& MINCORR(SRWMatrix& A, SRWVector& b, SRWMatrix& P,
                   SRWVector& x, double epsilon, int maxit);

#endif // SOLVERS_H

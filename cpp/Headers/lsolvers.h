#ifndef LSOLVERS_H
#define LSOLVERS_H

/**
 * Solvers of linear systems of type Ax = b
 * where A is square matrix, and x,b are vectors.
 *
 */

#include "Headers/srwvector.h"
#include "Headers/srwmatrix.h"
#include "Headers/eigen3vector.h"
#include "Headers/Preconditioning/preconditioner.h"
#include "Headers/Preconditioning/par2DPreconditioner.h"

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

SRWVector& MINCORR(SRWMatrix& A, SRWVector& f, Preconditioner& P,
                   SRWVector& x, double epsilon, int &maxit);

SRWVector& seq_par_MINCORR(SRWMatrix& A, SRWVector& f, par2DPreconditioner& P,
                   SRWVector& x, double epsilon, int &maxit);

SRWVector& MINCORR_omp(SRWMatrix& A, SRWVector& f, par2DPreconditioner& P,
                       SRWVector& x, double epsilon, int &maxit);


#endif // LSOLVERS_H

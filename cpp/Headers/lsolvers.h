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

#include <vector>

#include "omp.h"

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

void TDMA_opt(const double* a, const double* b, const double* c,
              double* x, const double* d, int n);

SRWVector& TDMA(SRWMatrix& A, SRWVector&x, SRWVector& b);

SRWVector& MINCORR(SRWMatrix& A, SRWVector& f, Preconditioner& P,
                   SRWVector& x, double epsilon, int &maxit);

SRWVector& seq_par_MINCORR(SRWMatrix& A, SRWVector& f, par2DPreconditioner& P,
                   SRWVector& x, double epsilon, int &maxit);

SRWVector& MINCORR_omp_slow(SRWMatrix& A, SRWVector& f, par2DPreconditioner& P,
                       SRWVector& x, double epsilon, int &maxit);

void MINCORR_opt(double* ap, double* an, double* as, double* ae, double* aw,
                       double* f, double* x, double** iP, double* r,
                       double* corr, double* Aw,
                       double epsilon, int& maxit, int ixs, int m);

void MINCORR_omp(double* ap, double* an, double* as, double* ae, double* aw,
                       double* f, double* x, double** iP, double* r,
                       double* corr, double* tmp_v, double* Aw,
                       double* dx_d, double* dx_l, double* dx_u,
                       double* dy_d, double* dy_l, double* dy_u,
                       double epsilon, int& maxit, int ixs, int m);

void MINRES(double* ap, double* an, double* as, double* ae, double* aw,
            double* f, double* x, double** iP, double* r,
            double* corr, double* Aw,
            double epsilon, int& maxit, int ixs, int m);

void MINRES_omp(double* ap, double* an, double* as, double* ae, double* aw,
                       double* f, double* x, double** iP, double* r,
                       double* corr, double* tmp_v, double* Aw,
                       double* dx_d, double* dx_l, double* dx_u,
                       double* dy_d, double* dy_l, double* dy_u,
                       double epsilon, int& maxit, int ixs, int m);

#endif // LSOLVERS_H

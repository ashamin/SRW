#ifndef SOLVERS_H
#define SOLVERS_H

#include "srwmatrix.h"
#include "tests/Poisson2DSquareArea/Test2DPoissonSquareArea.h"

SRWMatrix& form_A_matrix(int n, double h1, double h2,
                         int I);

SRWVector& form_b_vector(Test2DPoissonSquareArea& test, int n,
                         int I, int J, double h1, double h2,
                         SRWVector& x, SRWVector& y);

SRWMatrix& solve_poiss_2D_square(Test2DPoissonSquareArea& test,
                                 int I, int J, int &maxit);


#endif // SOLVERS_H

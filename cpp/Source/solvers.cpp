
#include "Headers/solvers.h"
#include "Headers/lsolvers.h"
#include "Headers/Preconditioning/preconditioners.h"
#include "Headers/eigen3matrix.h"

SRWMatrix& form_A_matrix(int n, double h1, double h2,
                         int I){
    SRWMatrix& A = *(new Eigen3Matrix(n, n));
    A.setZero();

    SRWVector& tmp_d = A.diag(0);
    tmp_d.fill(-(2/(h1*h1) + 2/(h2*h2)));
    A.setDiag(0, tmp_d);

    tmp_d = A.diag(1);
    tmp_d.fill(1/(h1*h1));
    for (int i = I-3; i<tmp_d.length(); i+=I-2)
        tmp_d(i) = 0;
    A.setDiag(1, tmp_d);
    A.setDiag(-1, tmp_d);

    tmp_d = A.diag(I-2);
    tmp_d.fill(1/(h2*h2));
    A.setDiag(I-2, tmp_d);
    A.setDiag(-(I-2), tmp_d);

    return A;
}


SRWVector& form_b_vector(Test2DPoissonSquareArea& test, int n,
                         double h1, double h2, int I, int J,
                         SRWVector& x, SRWVector& y){
    SRWVector& b = *(new Eigen3Vector(n));
    b.fill(0);
    int k = 0;

    for (int j = 1; j<J-1; j++)
        for (int i = 1; i<I-1; i++){
            k = (I-2)*(j-1) + (i-1);
            b(k) = -test.f(x(i), y(j));
            if (j == 1) b(k) = b(k) - test.g3(x(i))/(h2*h2);
            if (j == (J-2)) b(k) = b(k) - test.g4(x(i))/(h2*h2);
            if (i == 1) b(k) = b(k) - test.g1(y(j))/(h1*h1);
            if (i == (I-2)) b(k) = b(k) - test.g2(y(j))/(h1*h1);
        }

    return b;
}

SRWMatrix& solve_poiss_2D_square(Test2DPoissonSquareArea& test,
                                 int I, int J, int &maxit){
    double h1 = test.a / (I-1), h2 = test.b / (J-1);
    SRWVector& x = *(new Eigen3Vector(I));
    SRWVector& y = *(new Eigen3Vector(J));

    for (int i = 0; i<I; i++)
        x(i) = 0 + i*h1;
    for (int j = 0; j<J; j++)
        y(j) = 0 + j*h2;

    int n = (I-2)*(J-2);

    SRWMatrix& A = form_A_matrix(n, h1, h2, I);

    SRWVector& b = form_b_vector(test, n, h1, h2, I, J, x, y);

    SRWVector& solve = *(new Eigen3Vector(A.cols()));

    /*solve = MINCORR(A, b,
                    *(new Preconditioner(Preconditioners::SSOR_precond(A, 1))),
                    solve, 1e-5, maxit);*/

    /*solve = MINCORR(A, b,
                    *(new Preconditioner(A, 1, "SSOR")),
                    solve, 1e-5, maxit);*/

    /*solve = seq_par_MINCORR(A, b,
                    *(new par2DPreconditioner(.6, n, h1, "par.SSOR")),
                    solve, 1e-5, maxit);*/

    solve = MINCORR_omp_slow(A, b,
                    *(new par2DPreconditioner(.6, n, h1, "par.SSOR")),
                    solve, 1e-5, maxit);




    //HARDCORE CALL





    //HARDCORE CALL END





    A = A.split(solve, I-2, false);

    delete &x, &y, &A, &b;

    return A;
}



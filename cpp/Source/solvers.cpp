
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

int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
        n += 1;
    return n;
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

    /*solve = MINCORR_omp_slow(A, b,
                    *(new par2DPreconditioner(.6, n, h1, "par.SSOR")),
                    solve, 1e-5, maxit);*/




    //HARDCORE CALL

    par2DPreconditioner& precond = *(new par2DPreconditioner(.6, n, h1, "par.SSOR"));

    int ixs = I-2;

    double* ap = new double[n];
    double* as = new double[n];
    double* an = new double[n];
    double* ae = new double[n];
    double* aw = new double[n];
    double* f = new double[n];
    double* x_slv = new double[n];
    double **iP = new double*[n];
    for (int i = 0; i<n; i++)
            iP[i] = new double[n];
    double* r = new double[n];
    double* corr = new double[n];
    double* tmp_v = new double[n];
    double* Aw = new double[n];
    double* dx_d = new double[n];
    double* dx_l = new double[n-1];
    double* dx_u = new double[n-1];
    double* dy_d = new double[n];
    double* dy_l = new double[n-1];
    double* dy_u = new double[n-1];

    //fill parameters of MINCORR_omp()
    SRWVector& tv = A.diag(0);
    for (int i = 0; i<n; i++)
        ap[i] = tv(i);

    tv = A.diag(1);
    for (int i = 0; i<(n-1); i++)
        an[i] = tv(i);
    an[n-1] = 0;

    tv = A.diag(-1);
    as[0] = 0;
    for (int i = 1; i<n; i++)
        as[i] = tv(i-1);

    tv = A.diag(ixs);
    for (int i = 0; i<(n-ixs); i++)
        ae[i] = tv(i);
    for (int i = n-ixs; i<n; i++)
        ae[i] = 0;

    tv = A.diag(-ixs);
    for (int i = 0; i<ixs; i++)
        aw[i] = 0;
    for (int i = ixs; i<n; i++)
        aw[i] = tv(i-ixs);

    for (int i = 0; i<n; i++)
        f[i] = b(i);
    for (int i = 0; i<n; i++)
        x_slv[i] = solve(i);

    for (int i = 0; i<n; i++)
        for (int j = 0; j<n; j++)
            iP[i][j] = precond.iP()(i, j);

    tv = precond.Dx().diag(0);
    for (int i = 0; i<n; i++)
        dx_d[i] = tv[i];
    tv = precond.Dx().diag(-1);
    for (int i = 0; i<(n-1); i++)
        dx_l[i] = tv[i];
    tv = precond.Dx().diag(1);
    for (int i = 0; i<(n-1); i++)
        dx_u[i] = tv[i];


    tv = precond.Dy().diag(0);
    for (int i = 0; i<n; i++)
        dy_d[i] = tv[i];
    tv = precond.Dy().diag(-1);
    for (int i = 0; i<(n-1); i++)
        dy_l[i] = tv[i];
    tv = precond.Dy().diag(1);
    for (int i = 0; i<(n-1); i++)
        dy_u[i] = tv[i];

    /*print(ap, n);
    print(an, n);
    print(as, n);
    print(ae, n);
    print(aw, n);*/


    /*std::cout << omp_thread_count() << std::endl;

    omp_set_dynamic(0);
    omp_set_num_threads(omp_get_max_threads());

    double time = omp_get_wtime();
    MINCORR_omp(ap, an, as, ae, aw, f, x_slv, iP, r, corr, tmp_v,
                Aw, dx_d, dx_l, dx_u, dy_d, dy_l, dy_u, 1e-5, maxit, ixs, n);
    time = omp_get_wtime() - time;

    std::cout << "OMP_TIME =" << time << std::endl << std::endl;

    std::cout << omp_thread_count() << std::endl;*/




    Preconditioner& precond1 = *(new Preconditioner(A, 1, "SSOR"));
    for (int i = 0; i<n; i++)
        for (int j = 0; j<n; j++)
            iP[i][j] = precond1.iP()(i, j);

    clock_t ctime;
    ctime = clock();

    MINCORR_opt(ap, an, as, ae, aw, f, x_slv, iP, r, corr, Aw, 1e-5, maxit, ixs, n);

    ctime = clock() - ctime;
    std::cout << "ITERATIVE_SSOR = " <<  (double)ctime/CLOCKS_PER_SEC << std::endl<< std::endl;




    for(int i =0; i<n; i++)
        solve(i) = x_slv[i];


    //HARDCORE CALL END



    A = A.split(solve, I-2, false);

    delete &x, &y, &A, &b;

    return A;
}




#include "Headers/lsolvers.h"


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
                  SRWVector& c, SRWVector& d, SRWVector& x){

    int n = d.length();
    x.resize(n);
    double local_d[n];
    double tmp = 0;

    c(0) = c(0)/b(0);
    local_d[0] = d(0)/b(0);

    for (int i = 1; i<(n-1); i++){
        tmp = b(i) - c(i-1)*a(i-1);
        c(i) = c(i) / tmp;
        local_d[i] = (d(i) - local_d[i-1]*a(i-1)) / tmp;
    }

    local_d[n-1] = (d(n-1) - local_d[n-2]*a(n-2)) / (b(n-1) - c(n-2)*a(n-2));

    x(n-1) = local_d[n-1];
    for (int i = n-2; i>=0; i--)
        x(i) = local_d[i] - c(i)*x(i+1);

    return x;
}


SRWVector& TDMA(SRWMatrix& A, SRWVector&x, SRWVector& b){
    return TDMA_d(A.diag(-1), A.diag(0), A.diag(1), b, x);
}


SRWVector& MINCORR(SRWMatrix& A, SRWVector& f, Preconditioner& P,
                   SRWVector& x, double epsilon, int &maxit){
    double max_it_local = maxit;
    maxit = 0;
    double tau = 1;
    SRWMatrix& iP = P.iP();

    SRWVector& r = (f - A*x);
    SRWVector& w = iP*r;
    SRWVector& Aw = A*w;

    clock_t ctime;
    ctime = clock();
    while(maxit++ < max_it_local){

        tau = Aw.dot(w) / static_cast<SRWVector&>(iP*Aw).dot(Aw);

        x = x + w*tau;

        r = (f - A*x);
        w = iP*r;
        Aw = A*w;

        if (r.norm("m") < epsilon) break;
        if (maxit == max_it_local)
            std::cout << "Iteration process obviously won't converge. \\n Try to increase \" maxit \" value" << std::endl;
    }
    ctime = clock() - ctime;
    std::cout << "ITERATIVE_SSOR = " <<  (double)ctime/CLOCKS_PER_SEC << std::endl<< std::endl;

    return x;
}

SRWVector& seq_par_MINCORR(SRWMatrix& A, SRWVector& f, par2DPreconditioner& P,
                   SRWVector& x, double epsilon, int &maxit){

    double max_it_local = maxit;
    maxit = 0;
    double tau = 1;
    SRWMatrix& iP = P.iP();

    int n = sqrt(iP.rows());
    int k = 0;

    SRWVector& tmp_v = *(new Eigen3Vector(P.Dx().cols()));
    SRWVector& tmp_solve = *(new Eigen3Vector(P.Dx().cols()));
    SRWVector& corr = *(new Eigen3Vector(P.Dy().cols()));
    SRWVector *r, *Aw;


    clock_t ctime;
    ctime = clock();
    while(maxit++ < max_it_local){

        r = &(f - A*x);

        tmp_v.resize(0);

        for (int i = 0; i<n; i++){
            k = (i*n);
            tmp_v = tmp_v.glue(tmp_v, TDMA(P.Dx().subMatrix(k,k,n,n), tmp_solve, r->segment(k, n)));
        }

        corr.resize(0);

        for (int i = 0; i<n; i++){
            k = (i*n);
            corr = corr.glue(corr, TDMA(P.Dy().subMatrix(k,k,n,n), tmp_solve, tmp_v.segment(k, n)));
        }

        Aw = &(A*corr);

        tau = Aw->dot(corr) / static_cast<SRWVector&>(iP**Aw).dot(*Aw);

        x = x + corr*tau;

        if (r->norm("m") < epsilon) break;
        if (maxit == max_it_local)
            std::cout << "Iteration process obviously won't converge. \\n Try to increase \" maxit \" value" << std::endl;
    }
    ctime = clock() - ctime;
    std::cout << "SEQUENT_PAR_SSOR = " <<  (double)ctime/CLOCKS_PER_SEC << std::endl<< std::endl;

    return x;
}

SRWVector& MINCORR_omp(SRWMatrix& A, SRWVector& f, par2DPreconditioner& P,
                       SRWVector& x, double epsilon, int &maxit){

    double max_it_local = maxit;
    maxit = 0;
    double tau = 1;
    SRWMatrix& iP = P.iP();

    int n = sqrt(iP.rows());
    int i;
    int k = 0;

    /*SRWVector& tmp_v = *(new Eigen3Vector(P.Dx().cols()));
    SRWVector& tmp_solve = *(new Eigen3Vector(P.Dx().cols()));*/

    omp_set_dynamic(0);
    omp_set_num_threads(10);


    // here we allocate memory for all objects to decrease memory accessing
    std::vector <SRWVector*> tmp_v;
    for (int i = 0; i<n; i++)
        tmp_v.push_back(new Eigen3Vector(0));

    std::vector <SRWVector*> tmp_solve;
    for (int i = 0; i<n; i++)
        tmp_solve.push_back(new Eigen3Vector(0));

    std::vector <SRWVector*> tmp_corr;
    for (int i = 0; i<n; i++)
        tmp_corr.push_back(new Eigen3Vector(0));

    std::vector <SRWMatrix*> sub_mtrs_Dx;
    for (int i = 0; i<n; i++){
        k = (i*n);
        sub_mtrs_Dx.push_back(&P.Dx().subMatrix(k,k,n,n));
    }

    std::vector <SRWMatrix*> sub_mtrs_Dy;
    for (int i = 0; i<n; i++){
        k = (i*n);
        sub_mtrs_Dy.push_back(&P.Dy().subMatrix(k,k,n,n));
    }

    /*std::vector <SRWVector*> seg_r;
    for (inr i = 0; i<n; i++)
        seg_r.push_back(new Eigen3Vector(0));*/



    SRWVector& corr = *(new Eigen3Vector(P.Dy().cols()));
    SRWVector *r, *Aw;


    double time = omp_get_wtime();
    while(maxit++ < max_it_local){

        r = &(f - A*x);


        #pragma omp parallel for shared(sub_mtrs_Dx, r, tmp_v, tmp_solve) \
                                    firstprivate(n) private(i, k) \
                                    schedule(static)
        for (i = 0; i<n; i++){
            k = (i*n);
            *tmp_v.at(i) = TDMA(*sub_mtrs_Dx.at(i), *tmp_solve.at(i), r->segment(k, n));
        }

        #pragma omp parallel for shared(sub_mtrs_Dy, tmp_v, tmp_corr, tmp_solve) \
                                    firstprivate(n) private(i, k) \
                                    schedule(static)
        for (i = 0; i<n; i++){
            k = (i*n);
            *tmp_corr.at(i) = TDMA(*sub_mtrs_Dy.at(i), *tmp_solve.at(i), *tmp_v.at(i));
        }

        corr.resize(0);
        for (int i = 0; i<n; i++)
            corr = corr.glue(corr, *tmp_corr.at(i));


        Aw = &(A*corr);

        tau = Aw->dot(corr) / static_cast<SRWVector&>(iP**Aw).dot(*Aw);

        x = x + corr*tau;

        if (r->norm("m") < epsilon) break;
        if (maxit == max_it_local)
            std::cout << "Iteration process obviously won't converge. \\n Try to increase \" maxit \" value" << std::endl;
    }
    time = omp_get_wtime() - time;
    std::cout << "OMP_TIME = " << time << std::endl<< std::endl;

    return x;
}





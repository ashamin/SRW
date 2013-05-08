
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

    SRWVector& r_part = *(new Eigen3Vector(d.length()));
    r_part = d;
    int n = d.length();
    x.resize(n);
    double *tmp = new double[1];
    *tmp = 0;
    SRWVector& null_v = *(new Eigen3Vector(1));
    null_v(0) = 0;

    a = a.glue(null_v, a);
    c = c.glue(c, null_v);

    /*a = a.glue(tmp, a);
    c = a.glue(c, tmp);*/

    c(0) = c(0)/b(0);
    d(0) = d(0)/b(0);

    for (int i = 1; i<n; i++){
        *tmp = b(i) - c(i-1)*a(i);
        c(i) = c(i) / *tmp;
        d(i) = (d(i) - d(i-1)*a(i)) / *tmp;
    }

    x(n-1) = d(n-1);
    for (int i = n-2; i>=0; i--)
        x(i) = d(i) - c(i)*x(i+1);

    d = r_part;

    delete tmp, &null_v, &r_part;

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

    return x;
}

SRWVector& MINCORR_omp(SRWMatrix& A, SRWVector& f, par2DPreconditioner& P,
                       SRWVector& x, double epsilon, int &maxit){
    return *(new Eigen3Vector(0));
}





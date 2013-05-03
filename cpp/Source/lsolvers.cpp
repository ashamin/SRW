
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


SRWVector& MINCORR(SRWMatrix& A, SRWVector& f, SRWMatrix& P,
                   SRWVector& x, double epsilon, int maxit){
    double tau = 1;
    SRWMatrix& iP = P.inverse();

    SRWVector& r = (f - A*x);
    SRWVector& w = iP*r;
    SRWVector& Aw = A*w;

    while(maxit-- > 0){

        tau = Aw.dot(w) / static_cast<SRWVector&>(iP*Aw).dot(Aw);

        x = x + w*tau;

        r = (f - A*x);
        w = iP*r;
        Aw = A*w;

        if (r.absv().maxv() < epsilon) break;
        if (maxit == 0)
            std::cout << "Iteration process obviously won't converge. \\n Try to increase \" maxit \" value" << std::endl;
    }

    return x;
}





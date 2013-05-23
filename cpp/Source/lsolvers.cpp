
#include "Headers/lsolvers.h"

/**
 * @brief dvacm
 *          Direct vector addition with vector constant multiplication
 *              dest = dest + source*double_const;
 * @param dest
 * @param source
 * @param c
 */
void dvacm(SRWVector& dest, SRWVector& source, double c){
    int n = dest.length();
    //if (!(n - source.length())){
        for (int i = 0; i< n; i++)
            dest(i) += source(i)*c;
    //}
}

void dvsmvm(){

}


/**
 * @brief TDMA_d
 *          Tridiagonal Matrix Algorithm
 *              you should know x.length() before call
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
    //x.resize(n);
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

void TDMA_opt(const double *a, const double *b, const double *c,
              double *x, const double *d, int n){
    double loc_c[n];
    double loc_d[n];
    double tmp = 0;
    int i = 0;

    loc_c[0] = c[0] / b[0];
    loc_d[0] = d[0] / b[0];

    for (i = 1; i<(n-1); i++){
        tmp = b[0] - loc_c[i-1]*a[i-1];
        loc_c[i] = c[i] / tmp;
        loc_d[i] = (d[i] - loc_d[i-1]*a[i-1]) / tmp;
    }

    loc_d[n-1] = (d[n-1] - loc_d[n-2]*a[n-2]) / (b[n-1] - loc_c[n-2]*a[n-2]);

    x[n-1] = loc_d[n-1];
    for (i = n-2; i>=0; i--)
        x[i] = loc_d[i] - loc_c[i]*x[i+1];

    return;
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

SRWVector& MINCORR_omp_slow(SRWMatrix& A, SRWVector& f, par2DPreconditioner& P,
                       SRWVector& x, double epsilon, int &maxit){

    double max_it_local = maxit;
    maxit = 0;
    double tau = 1;
    SRWMatrix& iP = P.iP();

    int m = iP.rows();

    int n = sqrt(m);

    m /= n;

    int i;
    int k = 0;

    // here we allocate memory for all objects to decrease memory accessing
    std::vector <SRWVector*> tmp_v;
    for (int i = 0; i<n; i++)
        tmp_v.push_back(new Eigen3Vector(m));

    std::vector <SRWVector*> tmp_solve;
    for (int i = 0; i<n; i++)
        tmp_solve.push_back(new Eigen3Vector(m));

    std::vector <SRWVector*> tmp_corr;
    for (int i = 0; i<n; i++)
        tmp_corr.push_back(new Eigen3Vector(m));

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
            //TDMA(*sub_mtrs_Dx.at(i), *tmp_v.at(i), r->segment(k, n));
            TDMA(*sub_mtrs_Dx[i], *tmp_v[i], r->segment(k, n));
        }

        #pragma omp parallel for shared(sub_mtrs_Dy, tmp_v, tmp_corr, tmp_solve) \
                                    firstprivate(n) private(i, k) \
                                    schedule(static)
        for (i = 0; i<n; i++){
            k = (i*n);
            //TDMA(*sub_mtrs_Dy.at(i), *tmp_corr.at(i), *tmp_v.at(i));
            TDMA(*sub_mtrs_Dy[i], *tmp_corr[i], *tmp_v[i]);
        }

        corr.resize(0);
        for (int i = 0; i<n; i++)
            corr = corr.glue(corr, *tmp_corr.at(i));


        Aw = &(A*corr);

        tau = Aw->dot(corr) / static_cast<SRWVector&>(iP**Aw).dot(*Aw);

        //x = x + corr*tau;
        dvacm(x, corr, tau);

        if (r->norm("m") < epsilon) break;
        if (maxit == max_it_local)
            std::cout << "Iteration process obviously won't converge. \\n Try to increase \" maxit \" value" << std::endl;
    }
    time = omp_get_wtime() - time;
    std::cout << "OMP_TIME = " << time << std::endl<< std::endl;

    return x;
}

inline double maxd(double& v1, double& v2, double& tmp){
    tmp = fabs(v2);
    return (v1 > tmp) ? v1 : tmp;
}

// ixs = int_x_splits - means internal x splits. it's value computed by substracting
//  2 from x_splits. so it's value that shows number of x_splits without splits
//  near borders

// BEFORE CALLING THIS METHOD FILL NULLS IN DX_... and DY_... WHERE
void MINCORR_omp(double* ap, double* an, double* as, double* ae, double* aw,
                       double* f, double* x, SRWMatrix& iP, double* r,
                       double* corr, double* tmp_v, double* Aw,
                       double* dx_d, double* dx_l, double* dx_u,
                       double* dy_d, double* dy_l, double* dy_u,
                       double epsilon, int& maxit, int ixs){

    double max_it_local = maxit;
    maxit = 0;
    // tau parameter
    double tau = 1;
    // norm of residual vector r
    double rnorm = 0;
    double tmp = 0;
    // (Aw, corr) and (iP*Aw, Aw) dot products
    double dp_Aw_corr = 0, dp_iPmAw_Aw = 0;

    // matrix always is square so we can split this way
    int m = iP.rows();
    int n = sqrt(m), i, k = 0;

    while(maxit++ < max_it_local){

        //computing r and norm(r)
        //r = f - A*x
        rnorm = 0;
        r[0] = f[0] - ap[0]*x[0] - an[0]*x[1] - ae[0]*x[ixs];
        rnorm = maxd(rnorm, r[0], tmp);
        for (k = 1; k<ixs; k++){
            r[k] = f[k] - as[k-1]*x[k-1] - ap[k]*x[k] - an[k]*x[k+1] - ae[k]*x[k+ixs];
            rnorm = maxd(rnorm, r[k], tmp);
        }

        for (k = ixs; k<m-ixs; k++){
            r[k] = f[k] - as[k-1]*x[k-1] - ap[k]*x[k] - an[k]*x[k+1] - ae[k]*x[k+ixs]
                    - aw[k-ixs]*x[k-ixs];
            rnorm = maxd(rnorm, r[k], tmp);
        }

        for (k = m-ixs; k<m-1; i++){
            r[k] = f[k] - as[k-1]*x[k-1] - ap[k]*x[k] - an[k]*x[k+1] - aw[k-ixs]*x[k-ixs];
            rnorm = maxd(rnorm, r[k], tmp);
        }

        k = m-1;
        r[k] = f[k] - as[k-1]*x[k-1] - ap[k]*x[k] - aw[k-ixs]*x[k-ixs];
        rnorm = maxd(rnorm, r[k], tmp);


    //parallel computing corr - correction on each step using tridiagonal matrix method

        #pragma omp parallel for shared(dx_d, dx_l, dx_u, r, tmp_v) \
                                    firstprivate(n) private(i, k) \
                                    schedule(static)
        for (i = 0; i<n; i++){
            k = i*n;
            TDMA_opt(&dx_l[k], &dx_d[k], &dx_u[k], &tmp_v[k], &r[k], n);
        }

        #pragma omp parallel for shared(dy_d, dy_l, dy_u, tmp_v, corr) \
                                firstprivate(n) private(i, k) \
                                schedule(static)
        for (i = 0; i<n; i++){
            k = i*n;
            TDMA_opt(&dy_l[k], &dy_d[k], &dy_u[k], &corr[k], &tmp_v[k], n);
        }

    //computing Aw = A*corr and dot product (Aw, corr)
        dp_Aw_corr = 0;
        Aw[0] = ap[0]*corr[0] + an[0]*corr[1] + ae[0]*corr[ixs];
        dp_Aw_corr += Aw[0]*corr[0];
        for (k = 1; k<ixs; k++){
            Aw[k] = as[k-1]*corr[k-1] + ap[k]*corr[k] + an[k]*corr[k+1] + ae[k]*corr[k+ixs];
            dp_Aw_corr += Aw[k]*corr[k];
        }

        for (k = ixs; k<m-ixs; k++){
            Aw[k] = as[k-1]*corr[k-1] + ap[k]*corr[k] + an[k]*corr[k+1] + ae[k]*corr[k+ixs]
                    + aw[k-ixs]*corr[k-ixs];
            dp_Aw_corr += Aw[k]*corr[k];
        }

        for (k = m-ixs; k<m-1; i++){
            Aw[k] = as[k-1]*corr[k-1] + ap[k]*corr[k] + an[k]*corr[k+1] + aw[k-ixs]*corr[k-ixs];
            dp_Aw_corr += Aw[k]*corr[k];
        }

        k = m-1;
        Aw[k] = as[k-1]*corr[k-1] + ap[k]*corr[k] + aw[k-ixs]*corr[k-ixs];
        dp_Aw_corr += Aw[k]*corr[k];


    //computing dot product (iP*Aw; Aw) . maybe in same time than (Aw, corr) dot product


        //here U must compute dot product!!!


    //computing tau
        tau = dp_Aw_corr / dp_iPmAw_Aw;

    //computes new approximation of x
        for (k = 0; k<m; k++)
            x[k] += corr[k]*tau;

        if (rnorm < epsilon) break;
        if (maxit == max_it_local)
            std::cout << "Iteration process obviously won't converge. \\n Try to increase \" maxit \" value" << std::endl;
    }

    return;
}





#include "minres5dOmpSSOR.h"

#include "stdlib.h"

minres5dOmpSSOR::minres5dOmpSSOR(MathArea2d* const area, const SSORpar* const precond, const double epsilon, const int maxit)
{
    ap = area->getAp();
    an = area->getAn();
    as = area->getAs();
    ae = area->getAe();
    aw = area->getAw();
    f = area->getF();
    
    dx_d = precond->dx_d;
    dx_l = precond->dx_l;
    dx_u = precond->dx_u;
    dy_d = precond->dy_d;
    dy_l = precond->dy_l;
    dy_u = precond->dy_u;
    
    this->epsilon = epsilon;
    this->maxit = maxit;
    
    m = area->getN();
    
    x = new double[area->getI() * area->getJ()];
    corr = new double[area->getI() * area->getJ()];

    for (int i = 0; i<area->getI(); i++)
        x[i] = corr[i] = 0;
    for (int i = area->getI(); i<area->getI()+m; i++){
        x[i] = (double)i/m;
    }
    for (int i = area->getI()+m; i<area->getI() * area->getJ(); i++)
        x[i] = corr[i] = 0;

    
    r = new double[m];
    Aw = new double[m];
    
    ixs = area->getI() - 2;

	loc_c = new double[m];
	loc_d = new double[m];

	this->THREAD_NUM = omp_get_max_threads();
}

minres5dOmpSSOR::minres5dOmpSSOR(MathArea2d* const area, const SSORpar* const precond, const double epsilon, const int maxit, const int thread_num)
{
	minres5dOmpSSOR(area, precond, epsilon, maxit);
	this->THREAD_NUM = thread_num;
}

minres5dOmpSSOR::~minres5dOmpSSOR(){
    delete [] corr;
    delete [] r;
    delete [] Aw;
}

double minres5dOmpSSOR::exec_time()
{
    return time;
}

int minres5dOmpSSOR::it_num()
{
    return maxit;
}

double* minres5dOmpSSOR::solve(const int thread_num){
	int threadn_tmp = THREAD_NUM;
	THREAD_NUM = thread_num;
	double * ret = solve();
	THREAD_NUM = threadn_tmp;
	return ret;
}

double* minres5dOmpSSOR::solve()
{  
    // due to OpenMP does not allow class members in OpenMP clauses
    double* ap = this->ap;
    double* as = this->as;
    double* an = this->an;
    double* ae= this->ae;
    double* aw = this->aw;
    
    double* dx_d = this->dx_d;
    double* dx_l = this->dx_l;
    double* dx_u = this->dx_u;
    double* dy_d = this->dy_d;
    double* dy_l = this->dy_l;
    double* dy_u = this->dy_u;
    
    double* x = this->x;
    double* corr = this->corr;
    double* r = this->r;
    double* f = this->f;
    double* Aw = this->Aw;

	double* loc_c = this->loc_c;
	double* loc_d = this->loc_d;
    
    double epsilon = this->epsilon;
    int ixs = this->ixs;
    int m = this->m; 
    
    double* tmp_v = new double[m];
    
    double max_it_local = maxit;
    maxit = 0;
    // tau parameter
    double tau = 1;
    // norm of residual vector r
    double rnorm = 0;
    double tmp = 0;
    // (Aw, r) and (Aw, Aw) dot products
    double dp_Aw_r = 0, dp_Aw_Aw = 0;

    // matrix always is square so we can split this way
    int n = sqrt(m), i, j, k = 0;
    int I = ixs + 2;

    omp_set_dynamic(0);
    omp_set_num_threads(THREAD_NUM);

    double *norms = (double*)malloc(THREAD_NUM*sizeof(double));
    
    time = omp_get_wtime();
    
    //while(maxit++ < max_it_local){
    while (1){
    maxit++;
    
    //computing r and norm(r)
    //r = f - A*x
    rnorm = 0;
    for (i = 0; i<THREAD_NUM; i++)
        norms[i] = 0;

#pragma omp parallel for shared(r, f, as, ap, an, ae, aw, x, norms) \
firstprivate(m, ixs) private(k, tmp) \
schedule(static)
    for (k = 0; k<m; k++){
        r[k] = f[k] - as[k]*x[I + k-1] - ap[k]*x[I + k] - an[k]*x[I + k+1] - ae[k]*x[I + k+ixs]
            - aw[k]*x[I + k-ixs];

        tmp = fabs(r[k]);
        if (*norms < tmp) *norms = tmp;
    }

    for (i = 0; i<THREAD_NUM; i++)
        if (rnorm < norms[i]) rnorm = norms[i];

    //parallel computing corr - correction on each step using tridiagonal matrix method

#pragma omp parallel for shared(dx_d, dx_l, dx_u, dy_d, dy_l, dy_u, r, tmp_v) \
firstprivate(n) private(i, k) \
schedule(static)
    for (i = 0; i<n; i++) {
        k = i*n;
        TDMA(&dx_l[k], &dx_d[k], &dx_u[k], &tmp_v[k], &r[k], n, 1, &loc_c[k], &loc_d[k]);
    }

#pragma omp parallel for shared(dx_d, dx_l, dx_u, dy_d, dy_l, dy_u, r, tmp_v, corr) \
firstprivate(n) private(i, k) \
schedule(static)
    for (k = 0; k<n; k++) {
        TDMA(&dy_l[k], &dy_d[k], &dy_u[k], &corr[I + k], &tmp_v[k], n, n, &loc_c[k], &loc_d[k]);
    }

    //computing Aw = A*r and dot products (Aw, r) and (Aw, Aw)
    dp_Aw_r = 0;
    dp_Aw_Aw = 0;

#pragma omp parallel for shared(Aw, as, ap, an, ae, aw, corr) \
firstprivate(ixs, m) private(k) \
schedule(static) \
reduction(+:dp_Aw_r, dp_Aw_Aw)
    
    for (k = 0; k<m; k++){
        Aw[k] = as[k]*corr[I + k-1] + ap[k]*corr[I + k] + an[k]*corr[I + k+1] + ae[k]*corr[I + k+ixs]
        + aw[k]*corr[I + k-ixs];

        dp_Aw_r += Aw[k]*r[k];
        dp_Aw_Aw += Aw[k]*Aw[k];
    }

    //computing tau
    tau = dp_Aw_r / dp_Aw_Aw;

    //computes new approximation of x
#pragma omp parallel for shared(corr, x) \
firstprivate(m, tau) private(k) \
schedule(static)
    
    for (k = 0; k<m; k++)
        x[I + k] += corr[I + k]*tau;

    if (rnorm < epsilon) break;
    //if (maxit == max_it_local)
    //  std::cout << "Iteration process obviously won't converge. \\n Try to increase \" maxit \" value" << std::endl;
    }
    
    time = omp_get_wtime() - time;
    
    delete [] tmp_v;
	free(norms);

    return x;

}

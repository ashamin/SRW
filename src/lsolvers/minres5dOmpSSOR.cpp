#include "minres5dOmpSSOR.h"

minres5dOmpSSOR::minres5dOmpSSOR()
{

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
  
  double epsilon = this->epsilon;
  int maxit = this->maxit;
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

  while(maxit++ < max_it_local){

    //computing r and norm(r)
    //r = f - A*x
    rnorm = 0;

#pragma omp parallel for shared(r, f, as, ap, an, ae, aw, x) \
firstprivate(m, ixs) private(k, tmp) \
schedule(static)
    for (k = 0; k<m; k++){
      r[k] = f[k] - as[k]*x[k-1] - ap[k]*x[k] - an[k]*x[k+1] - ae[k]*x[k+ixs]
	      - aw[k]*x[k-ixs];

      tmp = fabs(r[k]);
      if (rnorm < tmp)
#pragma omp critical
      if (rnorm < tmp)
	rnorm = tmp;
    }

    //parallel computing corr - correction on each step using tridiagonal matrix method

#pragma omp parallel for shared(dx_d, dx_l, dx_u, dy_d, dy_l, dy_u, r, tmp_v, corr) \
firstprivate(n) private(i, k) \
schedule(static)

    for (i = 0; i<n; i++){
      k = i*n;
      TDMA(&dx_l[k], &dx_d[k], &dx_u[k], &tmp_v[k], &r[k], n);
      TDMA(&dy_l[k], &dy_d[k], &dy_u[k], &corr[k], &tmp_v[k], n);
    }

    //computing Aw = A*r and dot products (Aw, r) and (Aw, Aw)
    dp_Aw_r = 0;
    dp_Aw_Aw = 0;

#pragma omp parallel for shared(Aw, as, ap, an, ae, aw, corr) \
firstprivate(ixs, m) private(k) \
schedule(static) \
reduction(+:dp_Aw_r, dp_Aw_Aw)

    for (k = 0; k<m; k++){
      Aw[k] = as[k]*corr[k-1] + ap[k]*corr[k] + an[k]*corr[k+1] + ae[k]*corr[k+ixs]
	    + aw[k]*corr[k-ixs];

      dp_Aw_r += Aw[k]*r[k];
      dp_Aw_Aw += Aw[k]*Aw[k];
    }

    //computing tau
    tau = dp_Aw_r / dp_Aw_Aw;

    using namespace std;

    cout << tau << endl;

    //computes new approximation of x
#pragma omp parallel for shared(corr, x) \
firstprivate(m, tau) private(k) \
schedule(static)

    for (k = 0; k<m; k++)
      x[k] += corr[k]*tau;

    if (rnorm < epsilon) break;
    if (maxit == max_it_local)
      std::cout << "Iteration process obviously won't converge. \\n Try to increase \" maxit \" value" << std::endl;
  }
  
  delete tmp_v;

  return x;

}
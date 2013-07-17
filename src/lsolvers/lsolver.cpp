#include "lsolver.h"

void TDMA(const double *a, const double *b, const double *c,
              double *x, const double *d, int n)
{
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

double* lSolver::solve()
{
}

double* lSolver5d::solve()
{
}



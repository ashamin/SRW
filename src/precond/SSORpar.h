#ifndef SSOR_PAR_H
#define SSOR_PAR_H

class SSORpar{
public:
  SSORpar(double w, int n, double h);
  
  /** Diags of tridiagonal x-related matrix */
  double* dx_d;
  double* dx_l;
  double* dx_u;
  /** Diags of tridiagonal y-related matrix */
  double* dy_d;
  double* dy_l;
  double* dy_u;
};

#endif
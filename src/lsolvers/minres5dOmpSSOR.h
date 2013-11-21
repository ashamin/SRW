#ifndef MINRES_5D_OMP_SSOR_H
#define MINRES_5D_OMP_SSOR_H

#include "minres5d.h"
#include "MathArea2d.h"
#include "SSORpar.h"
#include "omp.h"

#include <iostream>
#include "math.h"

class minres5dOmpSSOR : public minres5d
{
public:
  minres5dOmpSSOR(MathArea2d* const area, const SSORpar* const precond, const double epsilon, const int maxit);
  minres5dOmpSSOR(MathArea2d* const area, const SSORpar* const precond, const double epsilon, const int maxit, const int thread_num);
  virtual ~minres5dOmpSSOR();

  double* solve();
  double* solve(const int thread_num);
  
  virtual double exec_time();
  int it_num();

  int THREAD_NUM;

private:
  /** Main diagonal of x related part of parallel tridiagonal SSOR preconditioner */
  double* dx_d;
  /** Lower diagonal of x */
  double* dx_l;
  /** Upper diagonal of x */
  double* dx_u;
  /** Main diagonal of y */
  double* dy_d;
  /** Lower diagonal of y */
  double* dy_l;
  /** Upper diagonal of y */
  double* dy_u;

  double* loc_c;
  double* loc_d;
  
  
  /** Execution time */
  double time; 
  
};

#endif
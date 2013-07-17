#ifndef MINRES_5D_OMP_SSOR_H
#define MINRES_5D_OMP_SSOR_H

#include "minres5d.h"
#include "omp.h"

#include <iostream>
#include "math.h"

class minres5dOmpSSOR : public minres5d
{
public:
  minres5dOmpSSOR();
  double* solve();

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
};

#endif
#ifndef MINRES_5D_H
#define MINRES_5D_H

#include "lsolver.h"

class minres5d : public lSolver5d{
public:
  minres5d();
  double* solve();
  
protected:
  /** Right part of matrix equation */
  double* f;
  /** Step residual */
  double* r;
  /** Step correction */
  double* corr;
  /** 5 diagonal matrix multiplied by correction */
  double* Aw;
  /** Precision */
  double epsilon;
  /** Max allowed iterations */
  int maxit;
  /** internal x splits. it's value computed by substracting 
   * 2 from x_splits. so it's value that shows number of x_splits
   * without splits near borders 
   */
  int ixs;
  /** Size of A matrix */
  int m;
};

#endif
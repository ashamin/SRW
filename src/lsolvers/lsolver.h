#ifndef LSOLVER_H
#define LSOLVER_H


/**
 * Base class for solvers of sparse linear systems.
 * \author ashamin
 */
class lSolver{
public:
  double* solve();
};

/**
 * Tridiagonal matrix algorithm. Solves tridiagonal matrix
 * with linear time.
 * \param a lower diagonal
 * \param b main diagonal
 * \param c upper diagonal 
 * \param x solve of system
 * \param d right part of matrix equation
 * \param n dimention of diagonal
 */
void TDMA_opt(const double* a, const double* b, const double* c,
              double* x, const double* d, int n);

#endif

#ifndef LSOLVER_5D_H
#define LSOLVER_5D_H

/**
 * Base class for linear solvers of sparse five-diagonal matrices.
 * \author ashamin
 * \TODO all stuff related with memory deleting
 */
class lSolver5d : public lSolver{
public:
  double* solve();
  
protected:
  /** Center of cross schema */
  double* ap;
  /** South point of cross schema */
  double* as;
  /** North point of cross schema */
  double* an;
  /** East point of cross schema */
  double* ae;
  /** West point of cross schema */
  double* aw;
};

#endif

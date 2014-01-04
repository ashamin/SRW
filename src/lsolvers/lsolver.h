#ifndef LSOLVER_H
#define LSOLVER_H


/**
 * Base class for solvers of sparse linear systems.
 * \author ashamin
 */
class lSolver{
public:
  double* solve();
  double exec_time();
};

/**
 * Прогонка
 * @param a нижняя диагональ
 * @param b главная диагональ
 * @param c верхняя диагональ
 * @param x вектор решения
 * @param b вектор правой части
 * @param n размерность вектора x
 * @param step шаг (используется для того, чтобы прогонку можно было гнать по X и по Y)
 * @param loc_c временный вектор
 * @param loc_d временный вектор
 */
void TDMA(const double* a, const double* b, const double* c,
              double* x, const double* d, int n, int step, double *loc_c, double *loc_d);

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
  virtual ~lSolver5d();
  
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
  /** Solve of matrix equation */
  double* x;
};

#endif

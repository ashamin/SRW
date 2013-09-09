#ifndef MATHAREA_2D_H
#define MATHAREA_2D_H

#include "math.h"

/** It's strongly recommended to add double a and b parameters
 * to class which will inherits this abstract class.
 * 
 * \param	a	it's left border of boundary value problem
 * \param	b	right border of boundary value problem
 * 
 * Also it's good choice to storage your matrix A and b vector in your class.
 */
class MathArea2d{
public:
  
  virtual double* getAp() = 0;
  virtual double* getAn() = 0;
  virtual double* getAs() = 0;
  virtual double* getAe() = 0;
  virtual double* getAw() = 0;
  
  virtual int getN() = 0;
  
  virtual int getI() = 0;
  virtual int getJ() = 0;
  
  virtual void formA() = 0;
  virtual void formB() = 0;
  
  virtual double g1(double y) = 0;
  virtual double g2(double y) = 0;
  virtual double g3(double x) = 0;
  virtual double g4(double x) = 0;

  virtual double f(double x, double y) = 0;

  virtual double right_answer(double x, double y) = 0;
};

#endif
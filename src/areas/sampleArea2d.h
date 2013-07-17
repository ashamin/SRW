#ifndef SAMPLE_AREA_2D_H
#define SAPMLE_AREA_2D_H

#include "MathArea2d.h"

class sampleArea2d : public MathArea2d{
public:
  virtual double** formA(int n, double h1, double h2, int I);
  virtual double* formB(int n, double h1, double h2, int I, int J, double* x, double* y);
  
  virtual double g1(double y);
  virtual double g2(double y);
  virtual double g3(double x);
  virtual double g4(double x);

  virtual double f(double x, double y);

  virtual double right_answer(double x, double y);  
  
  double a, b;
  
private:
  double** A;
  double* b;
};

#endif
#ifndef SAMPLE_AREA_2D_H
#define SAPMLE_AREA_2D_H

#include "MathArea2d.h"

class sampleArea2d : public MathArea2d{
public:
  
  sampleArea2d(int n, int ixs, int I, int J);
  
  virtual void formA(int n, double h1, double h2, int I);
  virtual void formB(int n, double h1, double h2, int I, int J, double* x, double* y);
  
  virtual double g1(double y);
  virtual double g2(double y);
  virtual double g3(double x);
  virtual double g4(double x);

  virtual double f(double x, double y);

  virtual double right_answer(double x, double y);  
  
  double ab, bb;
 
  // diagonals of A matrix
  double* ap;
  double* as;
  double* an;
  double* ae;
  double* aw;
  
  // right part of matirx equation
  double* b;
  
  // x and y coodinates
  double* x; 
  double* y;
  
  // step size x and y axes
  double h1;
  double h2;
  
  // ixs = int_x_splits - means internal x splits. it's value computed by substracting
  //  2 from x_splits. so it's value that shows number of x_splits without splits
  //  near borders 
  int ixs;
  
private:
  int p;
};

#endif
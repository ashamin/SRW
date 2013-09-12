#ifndef SAMPLE_AREA_2D_H
#define SAPMLE_AREA_2D_H

#include "MathArea2d.h"

class sampleArea2d : public MathArea2d{
public:
  
  sampleArea2d(int I, int J);

  virtual ~sampleArea2d();
    
  virtual double* getAp();
  virtual double* getAn();
  virtual double* getAs();
  virtual double* getAe();
  virtual double* getAw();
  virtual double* getF();
  
  virtual int getN();
  
  virtual int getI();
  virtual int getJ();
  
  virtual void formA();
  virtual void formB();
  
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
  
  // size
  int n;
  
  int I, J;
  
private:
  int p;
};

#endif
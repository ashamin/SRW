#include "sampleArea2d.h"

sampleArea2d::sampleArea2d(int I, int J)
{
  this->n = (I - 2)*(J - 2);
  this->ab = 1;
  this->bb = 1;
  // here some integer parameter from [16; 20]
  this->p = 19;
  
  this->I = I;
  this->J = J;

  ixs = I - 2;
  
  x = new double[n];
  y = new double[n];
  
  ap = new double[n];
  an = new double[n];
  as = new double[n];
  ae = new double[n];
  aw = new double[n];
  
  b = new double[n];
  
  h1 = ab / (I - 1);
  h2 = bb / (J - 1);
  
  for (int i = 0; i<I; i++)
    x[i] = i*h1;
  for (int j = 0; j<J; j++)
    y[j] = j*h2;
  
  this->formA();
  this->formB();
}

sampleArea2d::~sampleArea2d(){
  delete [] ap;
  delete [] as;
  delete [] an;
  delete [] ae;
  delete [] aw;

  delete [] b;
  delete [] x;
  delete [] y;
}

double* sampleArea2d::getAp()
{
  return ap;
}

double* sampleArea2d::getAn()
{
  return an;
}

double* sampleArea2d::getAs()
{
  return as;
}

double* sampleArea2d::getAe()
{
  return ae;
}

double* sampleArea2d::getAw()
{
  return aw;
}

double* sampleArea2d::getF(){
  return b;
}

int sampleArea2d::getN()
{
  return n;
}

int sampleArea2d::getI()
{
  return I;
}

int sampleArea2d::getJ()
{
  return J;
}

void sampleArea2d::formA()
{ 
  for (int i = 0; i<n; i++)
    ap[i] = -(2.0/(h1*h1) + 2.0/(h2*h2));
  
  for (int i = 0; i<n; i++)
    an[i] = as[i] = 1/(h1*h1);
  as[0] = 0;
  an[n-1] = 0;
  
  for (int i = 0; i<n-ixs; i++)
    ae[i] = 1/(h2*h2);
  for (int i = n-ixs; i<n; i++)
    ae[i] = 0;
  
  for (int i = 0; i<ixs; i++)
    aw[i] = 0;
  for (int i = ixs; i<n; i++)
    aw[i] = 1/(h2*h2);
}

void sampleArea2d::formB()
{
  int k = 0;

  for (int j = 1; j<J-1; j++){
    for (int i = 1; i<I-1; i++){
      k = (I-2)*(j-1) + (i-1);
      b[k] = -f(x[i], y[j]);
      if (j == 1) b[k] = b[k] - g3(x[i])/(h2*h2);
      if (j == (J-2)) b[k] = b[k] - g4(x[i])/(h2*h2);
      if (i == 1) b[k] = b[k] - g1(y[j])/(h1*h1);
      if (i == (I-2)) b[k] = b[k] - g2(y[j])/(h1*h1);
    }
  }
}

double sampleArea2d::g1(double y)
{
  return (-2)*y - 4*y*y;
}

double sampleArea2d::g2(double y)
{
  return 4 - 12*y - 4*y*y;
}

double sampleArea2d::g3(double x)
{
  return x + 3*x*x;
}

double sampleArea2d::g4(double x)
{
  return 3*x*x - 9*x - 6;
}

double sampleArea2d::f(double x, double y)
{
  return 2 + 2*M_PI*M_PI*(p-15)*sin(M_PI*x)*sin(M_PI*y);
}

double sampleArea2d::right_answer(double x, double y)
{
  return x + 3*x*x - 10*x*y - 2*y - 4*y*y + (p-15)*sin(M_PI*x)*sin(M_PI*y);
}

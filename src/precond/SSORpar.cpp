#include "SSORpar.h"

SSORpar::SSORpar(double w, int n, double h)
{

  double tmp_diag[n];
  dx_d = new double[n];
  dx_l = new double[n];
  dx_u = new double[n];
  dy_d = new double[n];
  dy_l = new double[n];
  dy_u = new double[n];
  
  for (int i = 0; i< n; i++)
    dx_d[i] = dy_d[i] = -4/(h*h)/w;
  
  for (int i = 0; i<n-1; i++)
    dx_u[i] = dy_u[i] = 2/(h*h);
  dx_u[n-1] = dy_u[n-1] = 0;
  
  dx_l[0] = dx_l[0] = 0;
  for (int i = 1; i<n; i++)
    dx_l[i] = dy_l[i] = 2/(h*h);
    
  for (int i = 0; i<n; i++)
    tmp_diag[i] = 1 / (dy_d[i] / w);
  
  for (int i = 0; i<n; i++){
    dy_d[i] = tmp_diag[i] * dy_d[i];
    dy_u[i] = tmp_diag[i] * dy_u[i];
    dy_l[i] = tmp_diag[i] * dy_l[i];
  }
  
}

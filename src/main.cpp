#include <iostream>

#include "lsolvers/lsolver.h"
#include "defs.h"

// temporary main file of srw project

int main(int argc, char **argv) {
    int n = 5;
    double* l = new double[n];
    double* u = new double[n];
    double* m = new double[n];
    double* d = new double[n];
    double* x = new double[n];
    
    l[0] = 0;
    for (int i = 0; i<n-1; i++){
      l[i+1] = i+1;
      u[i] = i;
    }
    u[n-1] = 0;
    for (int i = 0; i<n; i++)
      m[i] = i*2+3;
    
    for (int i = 0; i<n; i++)
      d[i] = i+1.5;

    for (int i = 0; i<n; i++)
      x[i] = 0;
    
    TDMA(l, m, u, x, d, n);
    
    using namespace std;
    
    for (int i = 0; i<n; i++)
      cout << x[i] << endl;
  
    std::cout << "Hello, srw!" << std::endl;
    return 0;
}

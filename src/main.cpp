#include <iostream>

#include "lsolvers/minres5dOmpSSOR.h"
#include "defs.h"

#include "precond/SSORpar.h"
#include "areas/sampleArea2d.h";

// temporary main file of srw project

int main(int argc, char **argv) { 
  
    using namespace std;
  
    int I = 20, J= 20;
    int n = (I- 2)*(J -2);
    
  
    sampleArea2d* area = new sampleArea2d(I, J);
    SSORpar* precond = new SSORpar(1, n, area->h1);
    minres5dOmpSSOR* solver = new minres5dOmpSSOR(area, precond, 1e-5, 10000);
    solver->solve();
    
    cout << "OMP_TIME=" << solver->exec_time() << endl;
    cout << "ITER_NUMBER=" << solver->it_num() << endl;
    
    cout << omp_thread_count() << endl;
    
  
    std::cout << "Hello, srw!" << std::endl;
    return 0;
}

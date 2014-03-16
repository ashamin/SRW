#include <iostream>

#include "minres5dOmpSSOR.h"
#include "defs.h"

#include "SSORpar.h"
#include "sampleArea2d.h"

#include "BArea.h"
#include "implexplBouss.h"

// temporary main file of srw project

int main(int argc, char **argv) { 
  
    using namespace std;
  
    int I = 100, J = 100, T = 2;
    int n = (I- 2)*(J -2);

    using namespace Boussinesq;

    //комментарий на русском. тест
    BArea* area           = new BArea(1, 1, 1, 6, 6, 50);
    implexplBouss* solver = new implexplBouss(area, 10e-5, 5);
    
    solver->solve();
  
//    sampleArea2d* area = new sampleArea2d(I, J);
//    SSORpar* precond = new SSORpar(.4, n, area->h1);
//    minres5dOmpSSOR* solver = new minres5dOmpSSOR(area, precond, 1e-5, 100000);
//    solver->solve(4);
    
//    cout << "OMP_TIME=" << solver->exec_time() << endl;
//    cout << "ITER_NUMBER=" << solver->it_num() << endl;
    
//    cout << omp_thread_count() << endl;

    return 0;
}

#include "lsolver.h"

/** TODO
 * 4) minres constructor
 */

void TDMA(const double *a, const double *b, const double *c,
        double *x, const double *d, int n, int step, double *loc_c, double *loc_d)
{
    double tmp = 0;
    int i = 0;

    loc_c[0] = c[0] / b[0];  // c[0] - всегда 1
    loc_d[0] = d[0] / b[0]; // d[0] - всегда 0

    for (i = step; i<n*step; i+=step){
        tmp = b[i] - loc_c[i-step]*a[i];
        loc_c[i] = c[i] / tmp;
        loc_d[i] = (d[i] - loc_d[i-step]*a[i]) / tmp;
    }

    x[n*step-step] = loc_d[n*step-step];
    for (i = (n-2)*step; i>=0; i-=step)
    x[i] = loc_d[i] - loc_c[i]*x[i+step];

    return;
}

lSolver5d::~lSolver5d(){
    delete [] ap;
    delete [] an;
    delete [] as;
    delete [] ae;
    delete [] aw;
    delete [] x;
}

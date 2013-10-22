#ifndef IMPLEXPL_BOUSS
#define IMPLEXPL_BOUSS

#include "minres5d.h"
#include "MathArea2d.h"
#include "SSORpar.h"
#include "omp.h"

#include <iostream>
#include "math.h"

class implexplBouss : public minres5d
{
public:
	implexplBouss(MathArea2d* const area, const SSORpar* const precond, const double epsilon, const int maxit);
	virtual ~implexplBouss();

	double* solve();
	
	virtual double exec_time();
	int it_num();

private:

	/** Center of cross schema */
	double* ap_y;
	/** South point of cross schema */
	double* as_y;
	/** North point of cross schema */
	double* an_y;
	/** East point of cross schema */
	double* ae_y;
	/** West point of cross schema */
	double* aw_y;

	double* xs;


	/** Main diagonal of x related part of parallel tridiagonal SSOR preconditioner */
	double* dx_d;
	/** Lower diagonal of x */
	double* dx_l;
	/** Upper diagonal of x */
	double* dx_u;
	/** Main diagonal of y */
	double* dy_d;
	/** Lower diagonal of y */
	double* dy_l;
	/** Upper diagonal of y */
	double* dy_u;
	
	
	/** Execution time */
	double time; 
	
};

#endif
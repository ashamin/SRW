#ifndef IMPLEXPL_BOUSS
#define IMPLEXPL_BOUSS

#include "BArea.h"
#include "lsolver.h"
#include "omp.h"

#include <iostream>
#include "math.h"

namespace Boussinesq{

// z ceiling and z floor
double zc = 1, zf = 0;
double mu1 = 1, mu2 = 2;
double kx = 1, ky = 1;

inline double get_mu(double H){
	if (H >= zc) return mu1;
	else return mu2;
}

inline double Tx(double H){
	if (H >= zc) return kx*(zc - zf);
	else if (H < zf) return 0;
	else return kx*(H - zf);
}

inline double Ty(double H){
	if (H >= zc) return ky*(zc - zf);
	else if (H < zf) return 0;
	else return ky*(H - zf);
}

class implexplBouss
{
public:
	implexplBouss(BArea* area, const double epsilon, const int maxit);
	virtual ~implexplBouss();

	double* solve();

	double borderValue(double x, double y, double t);
	void recomputeBorderValues();

	// whole H (with borders)
	double** Hw;
	// H only borders
	double* Hb;
	// H computation area
	double* H;
	// approximation of H
	double* Ha;

	double* b;
	
	double epsilon;
	int maxit;
	// Hw vector legth
	int n;
	// H vector length
	int m;
	double* V;
	
	double exec_time();
	int it_num();

	BArea* area;

	double* x;
	double* y;

	double t;
	double dt;

private:
	void prepareIteration();

	// diags of X differential operator
	double* dx_d;
	double* dx_l;
	double* dx_u;

	// diags of Y differentioal operator
	double* dy_d;
	double* dy_l;
	double* dy_u;

	//coefficient near time differential
	double* mu;

	// corresponds to I-2 and J-2
	// x - splits and y - splits
	int xs, ys;

	int I, J;

};

}

#endif
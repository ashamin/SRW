#ifndef IMPLEXPL_BOUSS
#define IMPLEXPL_BOUSS

#include "BArea.h"
#include "lsolver.h"
#include "omp.h"

#include <iostream>
#include "math.h"

namespace Boussinesq{

// z ceiling and z floor
const double zc = 1, zf = 0;
const double mu1 = 1, mu2 = 1;
const double kx = 1, ky = 1;

#define __get_mu(mu, H) (mu) = ((H) >= zc)?mu1:mu2;
#define __Tx(ret, H) (ret) = ((H) >= zc)?kx*(zc - zf):(((H) < zf)?0:kx*((H) - zf))
#define __Ty(ret, H) (ret) = ((H) >= zc)?ky*(zc - zf):(((H) < zf)?0:ky*((H) - zf))

inline void log(char *name, double* var, int size)
{
	std::cout << "RESULT__" << name << ':' << std::endl;
	for (int i = 0; i<size; i++)
		std::cout << var[i] << ' ';
	std::cout << std::endl;
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

	double* loc_c;
	double* loc_d;

	//coefficient near time differential
	double* mu;

	// corresponds to I-2 and J-2
	// x - splits and y - splits
	int xs, ys;

	int I, J;

};

}

#endif
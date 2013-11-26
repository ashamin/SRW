#include "implexplBouss.h"

#include "stdlib.h"

using namespace std;

namespace Boussinesq
{

implexplBouss::implexplBouss(BArea* area, const double epsilon, const int maxit)
{
    this->area = area;
    this->epsilon = epsilon;
    this->maxit = maxit;
    t = 0;
    dt = area->dt;

    I = area->I;
    J = area->J;
    n = I*J;

    x = new double[n];
    y = new double[n];

    x[0] = y[0] = 0;
    for (int i = 1; i<n; i++){
        x[i] = x[i-1] + area->hx;
        y[i] = y[i-1] + area->hy;
    }

    H	= new double[n];
	Ha  = new double[n];

    for (int j = 0; j<J; j++)
        for (int i = 0; i<I; i++)
            H[i + I*j] = area->answer(x[i], y[j], 0);



    b	  = new double[n];
    V	  = new double[n];
    dx_d  = new double[n];
    dx_l  = new double[n];
    dx_u  = new double[n];
    dy_d  = new double[n];
    dy_l  = new double[n];
    dy_u  = new double[n];
    mu	  = new double[n];

	loc_c = new double[n];
	loc_d = new double[n];
}
    
implexplBouss::~implexplBouss()
{

}

double implexplBouss::borderValue(double x, double y, double t)
{
    return area->answer(x, y, t);
}

void implexplBouss::recomputeBorderValues()
{

}

void implexplBouss::prepareIteration()
{

    //recomputeBorderValues();
    // V function internal area values
    for (int j = 0; j<J; j++)
        for (int i = 0; i<I; i++)
            V[i + I*j] = area->V(x[i], y[i], t);

	// TODO: initialize lower and upper squares. see diagram in notebook
    // here we form Tx, Ty and mu vaues on each step (parallel)

	// what we need here I or not I?
	for (int i = 1; i<I-1; i++) {
		int k = i*I;
		dx_l[k] = dx_d[k] = dy_l[k] = dy_d[k] = 0;
		dx_u[k] = dy_u[k] = 1;

		for (k = i*I+1; k<i*I+I-1; k++) {
            __Tx(dx_l[k], 
                        (H[k-1] + H[k])
                        /2
            );
            __Tx(dx_u[k],
                        (H[k+1] + H[k])
                        /2
            );
            dx_d[k] = dx_l[k] + dx_u[k];

            __Ty(dy_l[k],
                        (H[k-I] + H[k])
                        /2
            );
            __Ty(dy_u[k], 
                        (H[k+I] + H[k])
                        /2
            );
			dy_d[k] = dy_l[k] + dy_u[k];
		}

		k = i * I + I;
		dx_u[k] = dx_d[k] = dy_u[k] = dy_d[k] = 0;
		dx_l[k] = dy_l[k] = 1;

	}

    for (int i = 0; i<n; i++)
		__get_mu(mu[i], H[i]);

}

double* implexplBouss::solve()
{
  
    double* tmp_v = new double[n];
    int s		  = sqrt(n);


    while (true){

        prepareIteration();
		
        // explicit step
        int k = 0;
        for (k = I; k<=n-I; k++) {
                Ha[k] = (
                         dx_l[k]*H[k-1] + dx_d[k]*H[k] + dx_u[k]*H[k+1] +
                         dy_l[k]*H[k-I] + dy_d[k]*H[k] + dy_u[k]*H[k+I] + 
                         V[k]
                        ) 
                        / mu[k];
		}


        for (int i = I; i<n-I; i++) {
            double tmp = mu[i] / dt;
            dx_d[i] -= tmp;
            b[i] = - Ha[i] * tmp;
        }

        for (int i = 1; i<s-1; i++){
            int k = i*s;
            TDMA(&dx_l[k], &dx_d[k], &dx_u[k], &tmp_v[k], &b[k], s, &loc_c[k], &loc_d[k]);
        }

        for (int i = I; i<n-I; i++) {
			double tmp = mu[i] / dt;
			dy_d[i] -= tmp;
            tmp_v[i] = -tmp_v[i] * mu[i] / dt;
		}

        for (int i = 1; i<s-1; i++){
            int k = i*s;
            TDMA(&dy_l[k], &dy_d[k], &dy_u[k], &Ha[k], &tmp_v[k], s, &loc_c[k], &loc_d[k]);
        }


        for (int i = I; i<n-I; i++)
            H[i] = H[i] + dt*Ha[i];

		log("H", H, n);

        t += dt;

        break;

    }

    return H;
}

}
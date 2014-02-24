#include "implexplBouss.h"

#include "stdlib.h"
#include <cstring>

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

    I = area->I + 2;
    J = area->J + 2;
    n = I*J;

    x = new double[area->I];
    y = new double[area->J];

    x[0] = y[0] = 0;
    for (int i = 1; i<area->I; i++)
        x[i] = x[i-1] + area->hx;
    for (int j = 1; j<area->J; j++)
        y[j] = y[j-1] + area->hy;


    H   = new double[n];
    Ha  = new double[n];

    memset(H, 0, n*sizeof(double));
    memset(Ha, 0, n*sizeof(double));

    for (int j = 1; j<J-1; j++)
        for (int i = 1; i<I-1; i++)
            H[i + I*j] = area->answer(x[i-1], y[j-1], 0);

    b     = new double[n];
    V     = new double[n];
    dx_d  = new double[n];
    dx_l  = new double[n];
    dx_u  = new double[n];
    dy_d  = new double[n];
    dy_l  = new double[n];
    dy_u  = new double[n];
    mu    = new double[n];

    loc_c = new double[n];
    loc_d = new double[n];

    memset(b, 0, n*sizeof(double));
    memset(V, 0, n*sizeof(double));
    memset(dx_d, 0, n*sizeof(double));
    memset(dx_l, 0, n*sizeof(double));
    memset(dx_u, 0, n*sizeof(double));
    memset(dy_d, 0, n*sizeof(double));
    memset(dy_l, 0, n*sizeof(double));
    memset(dy_u, 0, n*sizeof(double));
    memset(mu, 0, n*sizeof(double));
    memset(loc_c, 0, n*sizeof(double));
    memset(loc_d, 0, n*sizeof(double));
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
    // пересчитываем функцию V а каждом шаге
    for (int j = 1; j<J-1; j++)
        for (int i = 1; i<I-1; i++)
            V[i + I*j] = area->V(x[i-1], y[i-1], t);

    // here we form Tx, Ty and mu vaues on each step (parallel)
    for (int j = 1; j<J-1; j++) {
        int k = j*I + 1;
        dx_l[k] = dy_l[k] = 0;
        dx_d[k] = dy_d[k] = 1;
        dx_u[k] = dy_u[k] = -1;

        for (k = j*I+2; k<j*I+I-2; k++) {
            __Tx(dx_l[k], 
                        (H[k-1] + H[k])
                        /2
            );
            __Tx(dx_u[k],
                        (H[k+1] + H[k])
                        /2
            );
            dx_d[k] = -dx_l[k] - dx_u[k];

            __Ty(dy_l[k],
                        (H[k-I] + H[k])
                        /2
            );
            __Ty(dy_u[k], 
                        (H[k+I] + H[k])
                        /2
            );
            dy_d[k] = -dy_l[k] - dy_u[k];
        }

        // непонятно назначение данного участка кода.
        //
        // обнуляет влияние граничного значения на приграничное. и пересчитывает его
        // с учетом только внутреннегзначения.
//        k = i*I+1;
//        dx_l[k] = dy_l[k] = 0;
//        dx_d[k] = -dx_l[k] - dx_u[k];
//        dy_d[k] = -dy_l[k] - dy_u[k];
//        k = i*I+I-2;
//        dx_u[k] = dy_u[k] = 0;
//        dx_d[k] = -dx_l[k] - dx_u[k];
//        dy_d[k] = -dy_l[k] - dy_u[k];

        k = j * I + I - 2;
        dx_l[k] = dy_l[k] = -1;
        dx_d[k] = dy_d[k] = 1;
        dx_u[k] = dy_u[k] = 0;

    }
    

    for (int i = 0; i<n; i++)
        __get_mu(mu[i], H[i]);

}

double* implexplBouss::solve()
{
  
    double* tmp_v = new double[n];
    int s         = (int)sqrt(n);


    while (true){


        prepareIteration();

        log_diags_as_3dmatrix("DX", dx_l, dx_d, dx_u, n);
        log_vector("MU", mu, n);
        log_matrix("V", V, n);

        int k = 0;
        for (k = I+1; k<=n-I-1; k++) {
                Ha[k] = (
                         dx_l[k]*H[k-1] + dx_d[k]*H[k] + dx_u[k]*H[k+1] +
                         dy_l[k]*H[k-I] + dy_d[k]*H[k] + dy_u[k]*H[k+I] + 
                         V[k]
                        ) 
                        / mu[k];
        }

        log_matrix("HA", Ha, n);

        // неявная прогонка по X
        for (int i = I; i<n-I; i++) {
            double tmp = mu[i] / dt;
            dx_d[i] -= tmp;
            b[i] = - Ha[i] * tmp;
        }

        log_diags_as_3dmatrix("DX", dx_l, dx_d, dx_u, n);
        log_vector("B", b, n);

        for (int i = 1; i<s-1; i++){
            int k = i*s+1;
            TDMA(&dx_l[k], &dx_d[k], &dx_u[k], &tmp_v[k], &b[k], s-2, 1, &loc_c[k], &loc_d[k]);
        }

        log_matrix("TMP", tmp_v, n);

        // неявная прогонка по Y
        for (int i = I; i<n-I; i++) {
            double tmp = mu[i] / dt;
            dy_d[i] -= tmp;
            tmp_v[i] = -tmp_v[i] * tmp;
        }

        for (int k = I+1; k<I+s-1; k++){
            TDMA(&dy_l[k], &dy_d[k], &dy_u[k], &Ha[k], &tmp_v[k], s-2, s, &loc_c[k], &loc_d[k]);
        }

        log_matrix("HA_TDMA", Ha, n);

        for (int i = I; i<n-I; i++)
            H[i] = H[i] + dt*Ha[i];

        log_matrix("H", H, n);

        t += dt;

        break;

    }

    return H;
}

}

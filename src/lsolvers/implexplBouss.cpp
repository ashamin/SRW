#include "implexplBouss.h"

#include "stdlib.h"

namespace Boussinesq{

implexplBouss::implexplBouss(BArea* area, const double epsilon, const int maxit){
    this->area = area;
    this->epsilon = epsilon;
    this->maxit = maxit;
    t = 0;
    dt = area->dt;

    I = area->I;
    J = area->J;
    xs = I-2;
    ys = J-2;
    n = I*J;
    m = xs*ys;

    x = new double[n];
    y = new double[n];
    x[0] = y[0] = 0;
    for (int i = 1; i<n; i++){
        x[i] = x[i-1] + area->hx;
        y[i] = y[i-1] + area->hy;
    }

    Hw = new double*[n];
    Hb = new double[n-m];
    H = new double[m];
    Ha = new double[m];

    for (int i = 0; i<I; i++)
        Hw[i] = &Hb[i];
    int iB = I;
    for (int j = 1; j<J-1; j++){
        Hw[I*j] = &Hb[iB++];
        Hw[I-1 + I*j] = &Hb[iB++];
    }
    for (int i = 0; i<I; i++)
        Hw[i + I*(J-1)] = &Hb[iB++];
    int iH = 0;
    for (int j = 1; j<J-1; j++)
        for (int i = 1; i<I-1; i++)
            Hw[i + I*j] = &H[iH++];

    for (int j = 0; j<J; j++)
        for (int i = 0; i<I; i++)
            *Hw[i + I*j] = area->answer(x[i], y[j], 0);

    b = new double[m];
    V = new double[m];
    dx_d = new double[m];
    dx_l = new double[m];
    dx_u = new double[m];
    dy_d = new double[m];
    dy_l = new double[m];
    dy_u = new double[m];
    mu = new double[m];
    prepareIteration();
}

implexplBouss::~implexplBouss(){

}

void implexplBouss::prepareIteration(){

    // V function internal area values
    for (int j = 1; j<J-1; j++)
        for (int i = 1; i<I-1; i++)
            V[i + I*j] = area->V(x[i], y[i], t);

    // here we form Tx, Ty and mu vaues on each step (parallel)
    int k = 0;
    for (int j = I+1; j<=n-I-xs; j+=I)
        for (int i = j; i<j+xs; i++){
            dx_l[k] = Tx((*Hw[i-1] + *Hw[i])/2);
            dx_u[k] = Tx((*Hw[i+1] + *Hw[i])/2);
            dx_d[k] = dx_l[k]+dx_u[k];
            dy_l[k] = Ty((*Hw[i-I] + *Hw[i])/2);
            dy_u[k] = Ty((*Hw[i+J] + *Hw[i])/2);
            dy_d[k] = dy_l[k]+dy_u[k];
            k++;
        }

    for (int i = 0; i<m; i++)
        mu[i] = get_mu(H[i]);
}

double* implexplBouss::solve(){

    double* tmp_v = new double[m];
    int s = sqrt(m);


    prepareIteration();

    // explicit step
    int k = 0;
    for (int j = I+1; j<=n-I-xs; j+=I)
        for (int i = j; i<j+xs; i++){
            Ha[k] = (dx_l[k]**Hw[k-1] + dx_d[k]**Hw[k] + dx_u[k]**Hw[k+1] +
                dy_l[k]**Hw[k-I] + dy_d[k]**Hw[k] + dy_u[k]**Hw[k+I] + V[k])/mu[k];
        }
            

    for (int i = 0; i<m; i++){
        double tmp = mu[i] / dt;
        dx_d[i] -= tmp;
        dy_d[i] -= tmp;
        b[i] = Ha[i]*tmp;
    }

    for (int i = 0; i<s; i++){
        int k = i*s;
        TDMA(&dx_l[k], &dx_d[k], &dx_u[k], &tmp_v[k], &b[k], s);
    }

    for (int i = 0; i<m; i++)
        tmp_v[i] = -tmp_v[i] * mu[i] / dt;

    for (int i = 0; i<s; i++){
        int k = i*s;
        TDMA(&dy_l[k], &dy_d[k], &dy_u[k], &Ha[k], &tmp_v[k], s);
    }

    for (int i = 0; i<m; i++)
        H[i] = H[i] + dt*Ha[i];


    return H;
}

}
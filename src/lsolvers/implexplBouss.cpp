#include "implexplBouss.h"

#include "stdlib.h"

namespace Boussinesq{

implexplBouss::implexplBouss(BArea* area, const double epsilon, const int maxit){
    this->epsilon = epsilon;
    this->maxit = maxit;
    this->n = (area->I-2)*(area->J-2);
    this->H = new double[n];
    for (int i = 0; i<area->I; i++)
        for (int j = 0; j<area->J; j++)
            H[j + area->J*i] = area->answer(area->hx*i, area->hy*j, 0);
    this->dx_d = new double[n];
    this->dx_l = new double[n];
    this->dx_u = new double[n];
    this->dy_d = new double[n];
    this->dy_l = new double[n];
    this->dy_u = new double[n];
    this->mu = new double[n];
    formDiffOperators();
}

implexplBouss::~implexplBouss(){

}

void implexplBouss::formDiffOperators(){
    // here we form Tx, Ty and mu vaues on each step (parallel)
    // here we also should insert zeros to lower and upper diags

    // what to do with H values on borders - add new array

    // add all computations for dx and dy diags for i = [0, n-1] values
    //      first and last rows of tridiagonal matrices

    dx_l[0] = dy_l[0] = dx_u[n-1] = dy_u[n-1] = 0;

    for (int i = 1; i<n-1; i++){
        dx_l[i] = Tx((H[i-1] + H[i])/2);
        dx_u[i] = Tx((H[i+1] + H[i])/2);
        dx_d[i] = dx_l[i]+dx_u[i];
    }

    for (int i = 1; i<n-1; i++){
        dy_l[i] = Ty((H[i-1] + H[i])/2);
        dy_u[i] = Ty((H[i+1] + H[i])/2);
        dy_d[i] = dy_l[i]+dy_u[i];
    }

    // here must fix border values of H

    for (int i = 0; i<n; i++)
        mu[i] = get_mu(H[i]);
}

double* implexplBouss::solve(){
    return NULL;
}

}
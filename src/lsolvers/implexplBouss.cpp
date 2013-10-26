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
    formDiffOperators();
}

implexplBouss::~implexplBouss(){

}

void implexplBouss::formDiffOperators(){
    // here we form Tx, Ty and mu vaues on each step (parallel)
}

double* implexplBouss::solve(){
    return NULL;
}

}
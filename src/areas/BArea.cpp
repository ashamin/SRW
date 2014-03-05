#include "BArea.h"

BArea::BArea(double sizeX, double sizeY, double destTime, int I, int J, double T){
    this->I = I;
    this->J = J;
    this->T = T;
    hx = sizeX / (I - 1);
    hy = sizeY / (J - 1);
    dt = destTime / (T - 1);
}

double BArea::answer(double x, double y, double t){
    return 1 + 3*t*t;
}

double BArea::V(double x, double y, double t){
    return 0;
}

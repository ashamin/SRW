#include "BArea.h"
#include <cmath>

BArea::BArea(double sizeX, double sizeY, double destTime, int I, int J, double T){
    this->I = I;
    this->J = J;
    this->T = T;
    hx = sizeX / (I - 1);
    hy = sizeY / (J - 1);
    dt = destTime / (T - 1);
}

double BArea::answer(double x, double y, double t){
//    return 1 + 3*t*t;
//    return cos(x*3.141592x65359)*cos(y*3.14159265359);
    return cos(3.14159265359 * x);
//    return 1 + 3*t*t;
//    return 1 + 3*x*x;
}

double BArea::V(double x, double y, double t){
//    return 1;
//    return -2 * cos(x*3.14159265359)*cos(y*3.14159265359) * (3.14159265359*3.14159265359);
    return 3.14159265359 * 3.14159265359 * cos(3.14159265359 * x);
//    return 6*t;
//    return -6;
}

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
    return cos(x*3.14159265359)*cos(y*3.14159265359);
}

double BArea::V(double x, double y, double t){
    return -2 * cos(x*3.14159265359)*cos(y*3.14159265359) * (3.14159265359*3.14159265359);
}

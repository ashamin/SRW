#include "BArea.h"
#include <cmath>

#include <iostream>
#include <fstream>

BArea::BArea(double sizeX, double sizeY, double destTime, int I, int J, double T){
    this->I = I;
    this->J = J;
    this->T = T;
    hx = sizeX / (I - 1);
    hy = sizeY / (J - 1);
    dt = destTime / (T - 1);
}

BArea::BArea(std::string file)
{
//    std::ifstream area;
//    area.open(file);
//    if(!myFile)
//    {
//        std::cerr << "Unable to open file " << file;
//        exit(1);
//    }

//    double sizeX = area.get();
//    hx = area.get();
//    double sizeY = area.get();
//    hy = area.get();
//    double sizeT = area.get();
//    dt = area.get();
//    for (int i = 0; i<sizeX; i+=hx)
//        for (int j = 0; j<sizeY; j+=hy){

//        }

//    area.close();
}

double BArea::answer(double x, double y, double t){
    return cos(x*3.14159265359)*cos(y*3.14159265359);
//    return cos(3.14159265359 * x);
//    return 1 + 3*y*y*y;
}

double BArea::V(double x, double y, double t){
    return 2*cos(x*3.14159265359)*cos(y*3.14159265359) * (3.14159265359*3.14159265359);
//    return 3.14159265359 * 3.14159265359 * cos(3.14159265359 * x);
//    return -18*y;
}

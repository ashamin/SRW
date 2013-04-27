#include "Headers/eigen3vector.h"

Eigen3Vector::Eigen3Vector(int length)
{
    this->size = length;
    this->vector = VectorXd(this->size);
}

double Eigen3Vector::get(int index){
    return vector(index);
}

int Eigen3Vector::length(){
    return this->size;
}

Eigen3Vector::~Eigen3Vector(){

}

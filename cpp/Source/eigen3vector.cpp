#include "Headers/eigen3vector.h"

Eigen3Vector::Eigen3Vector(int length)
{
    this->vector = VectorXd(length);
}

double Eigen3Vector::get(int index){
    return this->vector(index);
}

bool Eigen3Vector::set(int index, double value){
    if (index < this->length()){
        this->vector(index) = value;
        return true;
    }
    return false;
}

int Eigen3Vector::length(){
    return this->vector.rows();
}

SRWVector* Eigen3Vector::segment(int index, int length){
    Eigen3Vector* v = new Eigen3Vector(length);
    v->vector = this->vector.segment(index, length);
    return v;
}

Eigen3Vector& Eigen3Vector::operator= (const Eigen3Vector& v){
    if (this == &v)
        return *this;
    else{
        this->vector = v.vector;
        return *this;
    }
}

Eigen3Vector::~Eigen3Vector(){

}

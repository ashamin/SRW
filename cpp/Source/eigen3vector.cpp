#include "Headers/eigen3vector.h"

Eigen3Vector::Eigen3Vector(int length)
{
    this->vector = VectorXd::Random(length);
}

int Eigen3Vector::length(){
    return this->vector.rows();
}

SRWVector& Eigen3Vector::segment(int index, int length){
    Eigen3Vector& v = *(new Eigen3Vector(length));
    v.vector = this->vector.segment(index, length);
    return v;
}

SRWVector& Eigen3Vector::operator= (SRWVector& v){
    if (this == &v)
        return *this;
    else{
        this->vector =  dynamic_cast<Eigen3Vector&>(v).vector;
        return *this;
    }
    return *this;
}

/*SRWVector& Eigen3Vector::operator* (SRWVector& v){
    Eigen3Vector& ret_v = *(new Eigen3Vector(this->length()));
    ret_v.vector = this->vector * dynamic_cast<Eigen3Vector&>(v).vector;
    return ret_v;
}*/

SRWVector& Eigen3Vector::operator+ (SRWVector& v){
    Eigen3Vector& ret_v = *(new Eigen3Vector(this->length()));
    ret_v.vector = this->vector + dynamic_cast<Eigen3Vector&>(v).vector;
    return ret_v;
}

SRWVector& Eigen3Vector::operator- (SRWVector& v){
    Eigen3Vector& ret_v = *(new Eigen3Vector(this->length()));
    ret_v.vector = this->vector - dynamic_cast<Eigen3Vector&>(v).vector;
    return ret_v;
}

double& Eigen3Vector::operator [](int index){
    return this->vector(index);
}

double& Eigen3Vector::operator ()(int index){
    return this->vector(index);
}


Eigen3Vector::~Eigen3Vector(){

}

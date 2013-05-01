#include "Headers/eigen3vector.h"

Eigen3Vector::Eigen3Vector(int length)
{
    this->vector = VectorXd::Random(length);
}

int Eigen3Vector::length(){
    return this->vector.rows();
}

void Eigen3Vector::resize(int length){
    this->vector = VectorXd(length);
}

void Eigen3Vector::print(){
    std::cout << this->vector << std::endl << std::endl;
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

/*SRWVector& Eigen3Vector::glue(double* x, SRWVector& v){
    int x_length = sizeof(x) / 2;
    Eigen3Vector& ret_v = *(new Eigen3Vector(x_length + v.length()));
    for (int i = 0; i<x_length; i++)
        ret_v(i) = x[i];
    for (int i = x_length; i<x_length + v.length(); i++)
        ret_v(i) = v(i - x_length);
    return ret_v;
}

SRWVector& Eigen3Vector::glue(SRWVector& v, double* x){
    int x_length = sizeof(x) / 2;
    Eigen3Vector& ret_v = *(new Eigen3Vector(x_length + v.length()));
    for (int i = 0; i<v.length(); i++)
        ret_v(i) = v(i);
    for (int i = v.length(); i<v.length() + x_length; i++)
        ret_v(i) = x[i - v.length()];
    return ret_v;
}*/

SRWVector& Eigen3Vector::glue(SRWVector& v1, SRWVector& v2){
    Eigen3Vector& ret_v = *(new Eigen3Vector(v1.length() + v2.length()));
    for (int i = 0; i<v1.length(); i++)
        ret_v(i) = v1(i);
    for (int i = v1.length(); i< v1.length()+v2.length(); i++)
        ret_v(i) = v2(i - v1.length());
    return ret_v;
}

Eigen3Vector::~Eigen3Vector(){

}

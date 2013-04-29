#include "Headers/eigen3vector.h"

Eigen3Vector::Eigen3Vector(int length)
{
    this->vector = VectorXd::Random(length);
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

SRWVector& Eigen3Vector::operator= (SRWVector& v){
    if (this == &v)
        return *this;
    else{
        this->vector =  dynamic_cast<Eigen3Vector&>(v).vector;
        return *this;
    }
    return *this;
}

Eigen3Vector* Eigen3Vector::operator* (const Eigen3Vector* v2){
    Eigen3Vector* v = new Eigen3Vector(1);
    v->vector = this->vector * v2->vector;
    return v;
}

Eigen3Vector* Eigen3Vector::operator+ (const Eigen3Vector* v2){
    Eigen3Vector* v = new Eigen3Vector(1);
    v->vector = this->vector + v2->vector;
    return v;
}

Eigen3Vector* Eigen3Vector::operator- (const Eigen3Vector* v2){
    Eigen3Vector* v = new Eigen3Vector(1);
    v->vector = this->vector - v2->vector;
    return v;
}

Eigen3Vector::~Eigen3Vector(){

}

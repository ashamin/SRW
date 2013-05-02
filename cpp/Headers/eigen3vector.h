#ifndef EIGEN3VECTOR_H
#define EIGEN3VECTOR_H

#include "Headers/srwvector.h"
#include "Headers/eigen3matrix.h"
#include "Eigen/Dense"

#include <iostream>
#include "stdio.h"

using namespace Eigen;

class Eigen3Vector : public SRWVector
{
public:
    Eigen3Vector(int length);
    ~Eigen3Vector();

    virtual int length();

    virtual void resize(int length);

    virtual void print();

    virtual void fill(double value);

    virtual SRWVector& segment(int index, int length);

    virtual double dot(SRWVector& v);

    virtual SRWVector& absv();
    virtual double maxv();

    virtual SRWVector& operator= (SRWVector& v);

    //virtual SRWVector& operator* (SRWVector& v);
    virtual SRWVector& operator+ (SRWVector& v);
    virtual SRWVector& operator- (SRWVector& v);

    virtual SRWVector& operator *(double c);
    virtual SRWVector& operator /(double c);

    virtual double& operator [](int index);
    virtual double& operator ()(int index);

    /*virtual SRWVector& glue(double* x, SRWVector& v);
    virtual SRWVector& glue(SRWVector& v, double* x);*/
    virtual SRWVector& glue(SRWVector& v1, SRWVector& v2);




    VectorXd getVector(){
        return this->vector;
    }

    friend class Eigen3Matrix;

private:
    VectorXd vector;
};

#endif // EIGEN3VECTOR_H

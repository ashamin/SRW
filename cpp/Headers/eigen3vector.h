#ifndef EIGEN3VECTOR_H
#define EIGEN3VECTOR_H

#include "Headers/srwvector.h"
#include "Headers/eigen3matrix.h"
#include "Eigen/Dense"

using namespace Eigen;

class Eigen3Vector : public SRWVector
{
public:
    Eigen3Vector(int length);
    ~Eigen3Vector();

    virtual double get(int index);
    virtual bool set(int index, double value);
    virtual int length();

    virtual SRWVector& segment(int index, int length);

    virtual SRWVector& operator= (SRWVector& v);

    //virtual SRWVector& operator* (SRWVector& v);
    virtual SRWVector& operator+ (SRWVector& v);
    virtual SRWVector& operator- (SRWVector& v);

    VectorXd getVector(){
        return this->vector;
    }

    friend class Eigen3Matrix;

private:
    VectorXd vector;
};

#endif // EIGEN3VECTOR_H

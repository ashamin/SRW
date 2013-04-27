#ifndef EIGEN3VECTOR_H
#define EIGEN3VECTOR_H

#include "Headers/srwvector.h"
#include "Eigen/Dense"

using namespace Eigen;

class Eigen3Vector : public SRWVector
{
public:
    Eigen3Vector(int length);
    ~Eigen3Vector();

    virtual double get(int index);
    virtual int length();

private:
    VectorXd vector;
    int size;
};

#endif // EIGEN3VECTOR_H

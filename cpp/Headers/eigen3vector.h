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
    virtual bool set(int index, double value);
    virtual int length();

    virtual SRWVector* segment(int index, int length);

    virtual SRWVector& operator= (SRWVector& v);

    virtual Eigen3Vector* operator* (const Eigen3Vector* v2);
    virtual Eigen3Vector* operator+ (const Eigen3Vector* v2);
    virtual Eigen3Vector* operator- (const Eigen3Vector* v2);

    VectorXd getVector(){
        return this->vector;
    }

private:
    VectorXd vector;
};

#endif // EIGEN3VECTOR_H

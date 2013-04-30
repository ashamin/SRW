#ifndef EIGEN3MATRIX_H
#define EIGEN3MATRIX_H

#include "Headers/srwmatrix.h"
#include "Headers/srwvector.h"
#include "Headers/eigen3vector.h"
#include "Eigen/Dense"

using namespace Eigen;

class Eigen3Matrix : public SRWMatrix
{
public:
    Eigen3Matrix(int rows, int columns);

    virtual int rows();
    virtual int cols();

    virtual SRWVector& diag(int diagnum);
    virtual bool setDiag(int diagnum, SRWVector& v);

    virtual SRWMatrix& transpose();
    virtual SRWMatrix& inverse();

    virtual SRWMatrix& operator= (SRWMatrix& m);
    virtual SRWMatrix& operator *(SRWMatrix& m);
    virtual SRWVector& operator *(SRWVector& v);
    virtual SRWMatrix& operator +(SRWMatrix& m);
    virtual SRWMatrix& operator -(SRWMatrix& m);

    //virtual SRWVector& operator [](int index);
    virtual double& operator ()(int row, int column);

    MatrixXd getMatrix(){
        return this->matrix;
    }

private:
    MatrixXd matrix;
};

#endif // EIGEN3MATRIX_H

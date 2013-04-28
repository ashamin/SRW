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

    virtual double get(int row, int column);
    virtual bool set(int row, int column, double value);

    virtual SRWVector* diag(int diagnum);
    virtual bool setDiag(int diagnum, SRWVector* v);

    virtual SRWMatrix* transpose();
    virtual SRWMatrix* inverse();

    MatrixXd getMatrix(){
        return this->matrix;
    }

private:
    MatrixXd matrix;
};

#endif // EIGEN3MATRIX_H

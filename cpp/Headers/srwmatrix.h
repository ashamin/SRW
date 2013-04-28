#ifndef SRWMATRIX_H
#define SRWMATRIX_H

#include "Headers/srwvector.h"

class SRWMatrix
{
public:
    virtual int rows() = 0;
    virtual int cols() = 0;

    virtual double get(int row, int column) = 0;
    virtual bool set(int row, int column, double value) = 0;

    virtual SRWVector* diag(int diagnum) = 0;
    virtual bool setDiag(int diagnum, SRWVector* v) = 0;

    virtual SRWMatrix* transpose() = 0;
    virtual SRWMatrix* inverse() = 0;
};

#endif // SRWMATRIX_H

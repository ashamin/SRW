#ifndef SRWMATRIX_H
#define SRWMATRIX_H

#include "Headers/srwvector.h"

class SRWMatrix
{
public:

    virtual int rowsNum() = 0;
    virtual int columnsNum() = 0;

    virtual double get(int row, int column) = 0;
    virtual double set(int row, int column, int value) = 0;

    virtual SRWVector* diag() = 0;
    virtual SRWVector* diag(int diagnum) = 0;
    virtual void setDiag(int diagnum, SRWVector* v) = 0;

    virtual SRWMatrix* transpose(SRWMatrix* matrix) = 0;
};

#endif // SRWMATRIX_H

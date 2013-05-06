#ifndef SRWMATRIX_H
#define SRWMATRIX_H

#include "Headers/srwvector.h"

class SRWMatrix
{
public:
    virtual int rows() = 0;
    virtual int cols() = 0;

    virtual void print() = 0;

    virtual void setZero() = 0;
    virtual void setOne() = 0;

    virtual SRWVector& diag(int diagnum) = 0;
    virtual bool setDiag(int diagnum, SRWVector& v) = 0;

    virtual SRWMatrix& split(SRWVector& v, int length, bool byrow) = 0;
    virtual SRWMatrix& subMatrix(int upperleft_corner_row, int upperleft_corner_col, int rows, int cols) = 0;

    virtual double norm(std::string n_type) = 0;

    virtual SRWMatrix& transpose() = 0;
    virtual SRWMatrix& inverse() = 0;

    virtual SRWMatrix& diag_m() = 0;
    virtual SRWMatrix& lower_tri(bool strictly) = 0;
    virtual SRWMatrix& upper_tri(bool strictly) = 0;

    virtual SRWMatrix& operator= (SRWMatrix& m) = 0;

    virtual bool operator ==(SRWMatrix& m) = 0;

    virtual SRWMatrix& operator* (SRWMatrix& m) = 0;
    virtual SRWVector& operator* (SRWVector& v) = 0;
    virtual SRWMatrix& operator+ (SRWMatrix& m) = 0;
    virtual SRWMatrix& operator- (SRWMatrix& m) = 0;

    virtual SRWMatrix& operator *(double c) = 0;
    virtual SRWMatrix& operator /(double c) = 0;

    //virtual SRWVector& operator [](int index) = 0;
    virtual double& operator ()(int row, int column) = 0;
};

#endif // SRWMATRIX_H

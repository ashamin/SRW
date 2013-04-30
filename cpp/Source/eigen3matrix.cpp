#include "Headers/eigen3matrix.h"

Eigen3Matrix::Eigen3Matrix(int rows, int columns)
{
    this->matrix = MatrixXd::Random(rows, columns);
}

int Eigen3Matrix::rows(){
    return this->matrix.rows();
}

int Eigen3Matrix::cols(){
    return this->matrix.cols();
}

SRWVector& Eigen3Matrix::diag(int diagnum){
    Eigen3Vector& ret_v = *(new Eigen3Vector(1));
    ret_v.vector = this->matrix.diagonal(diagnum);
    return ret_v;
}

bool Eigen3Matrix::setDiag(int diagnum, SRWVector &v){
    if (this->diag(diagnum).length() == v.length()){
        this->matrix.diagonal(diagnum) = dynamic_cast<Eigen3Vector&>(v).vector;
        return true;
    }
    return false;
}

SRWMatrix& Eigen3Matrix::transpose(){
    Eigen3Matrix& m = *(new Eigen3Matrix(this->cols(), this->rows()));
    m.matrix = this->matrix.transpose();
    return m;
}

SRWMatrix& Eigen3Matrix::inverse(){
    Eigen3Matrix& m = *(new Eigen3Matrix(this->cols(), this->rows()));
    m.matrix = this->matrix.inverse();
    return m;
}

SRWMatrix& Eigen3Matrix::operator =(SRWMatrix& m){
    if (this == &m){
        return *this;
    }
    else{
        this->matrix = dynamic_cast<Eigen3Matrix&>(m).matrix;
        return *this;
    }
}

SRWMatrix& Eigen3Matrix::operator *(SRWMatrix& m){
    Eigen3Matrix& ret_m= *(new Eigen3Matrix(this->rows(), m.cols()));
    ret_m.matrix = this->matrix * dynamic_cast<Eigen3Matrix&>(m).matrix;
    return ret_m;
}

SRWVector& Eigen3Matrix::operator *(SRWVector& v){
    Eigen3Vector& ret_v = *(new Eigen3Vector(this->rows()));
    ret_v.vector = this->matrix * dynamic_cast<Eigen3Vector&>(v).vector;
    return ret_v;
}

SRWMatrix& Eigen3Matrix::operator +(SRWMatrix& m){
    Eigen3Matrix& ret_m= *(new Eigen3Matrix(this->rows(), m.cols()));
    ret_m.matrix = this->matrix + dynamic_cast<Eigen3Matrix&>(m).matrix;
    return ret_m;
}

SRWMatrix& Eigen3Matrix::operator -(SRWMatrix& m){
    Eigen3Matrix& ret_m= *(new Eigen3Matrix(this->rows(), m.cols()));
    ret_m.matrix = this->matrix * dynamic_cast<Eigen3Matrix&>(m).matrix;
    return ret_m;
}

double& Eigen3Matrix::operator ()(int row, int column){
    return this->matrix(row, column);
}

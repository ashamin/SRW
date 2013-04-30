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

double Eigen3Matrix::get(int row, int column){
    return this->matrix(row, column);
}

bool Eigen3Matrix::set(int row, int column, double value){
    if (row < this->matrix.rows() && column < this->matrix.cols()){
        this->matrix(row, column) = value;
        return true;
    }
    return false;
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

SRWMatrix* Eigen3Matrix::transpose(){
    Eigen3Matrix* m = new Eigen3Matrix(this->cols(), this->rows());
    m->matrix = this->matrix.transpose();
    return m;
}

SRWMatrix* Eigen3Matrix::inverse(){
    Eigen3Matrix* m = new Eigen3Matrix(this->cols(), this->rows());
    m->matrix = this->matrix.inverse();
    return m;
}

Eigen3Matrix& Eigen3Matrix::operator= (const Eigen3Matrix& m){
    if (this == &m){
        return *this;
    }
    else{
        this->matrix = m.matrix;
        return *this;
    }
}

Eigen3Matrix* Eigen3Matrix::operator* (const Eigen3Matrix* m2){
    Eigen3Matrix* m = new Eigen3Matrix(1, 21);
    m->matrix = this->matrix * m2->matrix;
    return m;
}

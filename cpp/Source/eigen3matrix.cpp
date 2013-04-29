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

SRWVector* Eigen3Matrix::diag(int diagnum){
    VectorXd d = this->matrix.diagonal(diagnum);
    Eigen3Vector* v = new Eigen3Vector(d.rows());
    for (int i = 0; i<v->length(); i++)
        v->set(i, d(i));
    return v;
}

bool Eigen3Matrix::setDiag(int diagnum, SRWVector *v){
    if (((diagnum < 0 && diagnum < this->rows()) ||
            (diagnum > 0 && diagnum < this->cols())) &&
            this->diag(diagnum)->length() == v->length()){

        VectorXd d(v->length());
        for (int i = 0; i<v->length(); i++)
            d(i) = v->get(i);
        this->matrix.diagonal(diagnum) = d;
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

#ifndef TEST2DPOISSONSQUAREAREA_H
#define TEST2DPOISSONSQUAREAREA_H

#include "Headers/srwmatrix.h"
#include "Headers/eigen3matrix.h"
#include "Headers/srwvector.h"
#include "Headers/eigen3vector.h"

class Test2DPoissonSquareArea
{
public:
    virtual double g1(double y) = 0;
    virtual double g2(double y) = 0;
    virtual double g3(double x) = 0;
    virtual double g4(double x) = 0;

    virtual double f(double x, double y) = 0;

    virtual double right_answer(double x, double y) = 0;

    double a, b;

    SRWMatrix& get_right_answer_as_matrix(int I, int J){

        double h1 = this->a / (I-1), h2 = this->b / (J-1);
        SRWVector& x = *(new Eigen3Vector(I));
        SRWVector& y = *(new Eigen3Vector(J));

        for (int i = 0; i<I; i++)
            x(i) = 0 + i*h1;
        for (int j = 0; j<J; j++)
            y(j) = 0 + j*h2;

        //int n = (I-2)*(J-2);

        SRWMatrix& r_ans = *(new Eigen3Matrix(I-2, J-2));

        for (int i = 1; i<I-1; i++)
            for (int j = 1; j<J-1; j++)
                r_ans(i-1, j-1) = this->right_answer(x(i), y(j));

        delete &x, &y;

        return r_ans;
    }
};

#endif // TEST2DPOISSONSQUAREAREA_H

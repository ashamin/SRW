#ifndef TEST2DPOISSONSQUAREAREAN5_H
#define TEST2DPOISSONSQUAREAREAN5_H

#include "Test2DPoissonSquareArea.h"
#include "math.h"

class Test2DPoissonSquareAreaN5 : public Test2DPoissonSquareArea
{
public:

    Test2DPoissonSquareAreaN5(){
        this->a = 1;
        this->b = 1;
        //p parameter is integer values from [21; 25]
        this->p = 23;
    }

    virtual double g1(double y){
        return exp(-y*y);
    }

    virtual double g2(double y){
        return y/p + exp(1-y*y);
    }

    virtual double g3(double x){
        return exp(x*x);
    }

    virtual double g4(double x){
        return x/p + exp(x*x - 1);
    }

    virtual double f(double x, double y){
        return -4*x*x*exp(x*x - y*y) - 4*y*y*exp(x*x - y*y);
    }

    virtual double right_answer(double x, double y){
        return x*y/p + exp(x*x - y*y);
    }


private:
    //p parameter is integer values from [21; 25]
    int p;
};

#endif // TEST2DPOISSONSQUAREAREAN5_H

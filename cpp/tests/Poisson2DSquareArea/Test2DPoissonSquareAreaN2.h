#ifndef TEST2DPOISSONSQUAREAREAN2_H
#define TEST2DPOISSONSQUAREAREAN2_H

#include "Test2DPoissonSquareArea.h"
#include "math.h"

class Test2DPoissonSquareAreaN2 : public Test2DPoissonSquareArea
{
public:

    Test2DPoissonSquareAreaN2(){
        this->a = 1;
        this->b = M_PI_2;
        //p parameter is integer values from [1; 5]
        this->p = 3;
    }

    virtual double g1(double y){
        return sin(p*y);
    }

    virtual double g2(double y){
        return exp(p) * sin(p*y);
    }

    virtual double g3(double x){
        return 0;
    }

    virtual double g4(double x){
        return exp(p*x)*sin(p*M_PI_2);
    }

    virtual double f(double x, double y){
        return 0;
    }

    virtual double right_answer(double x, double y){
        return exp(p*x)*sin(p*y);
    }


private:
    //p parameter is integer values from [1; 5]
    int p;
};

#endif // TEST2DPOISSONSQUAREAREAN2_H

#ifndef TEST2DPOISSONSQUAREAREAN3_H
#define TEST2DPOISSONSQUAREAREAN3_H

#include "Test2DPoissonSquareArea.h"
#include "math.h"

class Test2DPoissonSquareAreaN3 : public Test2DPoissonSquareArea
{
public:

    Test2DPoissonSquareAreaN3(){
        this->a = 2;
        this->b = 1;
        //p parameter is integer values from [6; 10]
        this->p = 8;
    }

    virtual double g1(double y){
        return 0;
    }

    virtual double g2(double y){
        return p*y/4 + 4*cos(y);
    }

    virtual double g3(double x){
        return 2*x;
    }

    virtual double g4(double x){
        return p*x/8 + 2*x*cos(1);
    }

    virtual double f(double x, double y){
        return 2*x*cos(y);
    }

    virtual double right_answer(double x, double y){
        return p*x*y/8 + 2*x*cos(y);
    }


private:
    //p parameter is integer values from [6; 10]
    int p;
};

#endif // TEST2DPOISSONSQUAREAREAN3_H

#ifndef TEST2DPOISSONSQUAREAREAN1_H
#define TEST2DPOISSONSQUAREAREAN1_H

#include "Test2DPoissonSquareArea.h"
#include "math.h"

class Test2DPoissonSquareAreaN1 : public Test2DPoissonSquareArea
{
public:

    Test2DPoissonSquareAreaN1(){
        this->a = 1;
        this->b = 1;
        //p parameter is integer values from [16; 20]
        this->p = 19;
    }

    virtual double g1(double y){
        return (-2)*y - 4*y*y;
    }

    virtual double g2(double y){
        return 4 - 12*y - 4*y*y;
    }

    virtual double g3(double x){
        return x + 4*x*x;
    }

    virtual double g4(double x){
       return 3*x*x - 9*x - 6;
    }

    virtual double f(double x, double y){
        return 2 + 2*M_PI*M_PI*(p-15)*sin(M_PI*x)*sin(M_PI*y);
    }

    virtual double right_answer(double x, double y){
        return x + 3*x*x - 10*x*y - 2*y - 4*y*y + (p-15)*sin(M_PI*x)*sin(M_PI*y);
    }


private:
    //p parameter is integer values from [16; 20]
    int p;
};

#endif // TEST2DPOISSONSQUAREAREAN1_H

#ifndef TEST2DPOISSONSQUAREAREAN4_H
#define TEST2DPOISSONSQUAREAREAN4_H

#include "Test2DPoissonSquareArea.h"
#include "math.h"

class Test2DPoissonSquareAreaN4 : public Test2DPoissonSquareArea
{
public:

    Test2DPoissonSquareAreaN4(){
        this->a = 1;
        this->b = 1;
        //p parameter is integer values from [11; 15]
        this->p = 12;
    }

    virtual double g1(double y){
        return -sqrt(p)*y*y;
    }

    virtual double g2(double y){
        return (sqrt(p) + y)*(1 - y*y);
    }

    virtual double g3(double x){
        return sqrt(p)*x*x;
    }

    virtual double g4(double x){
        return (sqrt(p) + x)*(x*x - 1);
    }

    virtual double f(double x, double y){
        return 0;
    }

    virtual double right_answer(double x, double y){
        return (sqrt(p) + x*y) * (x*x - y*y);
    }


private:
    //p parameter is integer values from [11; 15]
    int p;
};

#endif // TEST2DPOISSONSQUAREAREAN4_H

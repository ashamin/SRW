#ifndef TEST2DPOISSONSQUAREAREA_H
#define TEST2DPOISSONSQUAREAREA_H

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
};

#endif // TEST2DPOISSONSQUAREAREA_H

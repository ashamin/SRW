#ifndef SRWVECTOR_H
#define SRWVECTOR_H

class SRWVector
{
public:
    virtual int length() = 0;

    virtual void resize(int length) = 0;

    virtual void print() = 0;

    /**
     * @brief segment
     * @param index
     *          index of first element of segment
     * @param length
     *          length of segment
     * @return
     */
    virtual SRWVector& segment(int index, int length) = 0;

    virtual SRWVector& operator= (SRWVector& v) = 0;
    //still didn't found how transpose vector.
    //virtual SRWVector& operator* (SRWVector& v) = 0;
    virtual SRWVector& operator+ (SRWVector& v) = 0;
    virtual SRWVector& operator- (SRWVector& v) = 0;

    virtual double& operator [](int index) = 0;
    virtual double& operator ()(int index) = 0;

    /*virtual SRWVector& glue(double* x, SRWVector& v) = 0;
    virtual SRWVector& glue(SRWVector& v, double* x) = 0;*/
    virtual SRWVector& glue(SRWVector& v1, SRWVector& v2) = 0;

};

#endif // SRWVECTOR_H

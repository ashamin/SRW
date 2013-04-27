#ifndef SRWVECTOR_H
#define SRWVECTOR_H

class SRWVector
{
public:
    virtual double get(int index) = 0;
    virtual int length() = 0;
};

#endif // SRWVECTOR_H

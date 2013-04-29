#ifndef SRWVECTOR_H
#define SRWVECTOR_H

class SRWVector
{
public:
    virtual double get(int index) = 0;
    virtual bool set(int index, double value) = 0;
    virtual int length() = 0;

    /**
     * @brief segment
     * @param index
     *          index of first element of segment
     * @param length
     *          length of segment
     * @return
     */
    virtual SRWVector* segment(int index, int length) = 0;

    //virtual SRWVector& operator= (const SRWVector& v) = 0;
};

#endif // SRWVECTOR_H

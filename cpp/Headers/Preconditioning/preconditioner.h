#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

#include "Headers/srwmatrix.h"
#include "Headers/eigen3matrix.h"

class Preconditioner{
public:
    Preconditioner(){
        this->Pr = new Eigen3Matrix(2, 2);
        this->iPr = &Pr->inverse();
   }

    Preconditioner(SRWMatrix& m){
        Pr = &m;
        iPr = &Pr->inverse();
    }

    Preconditioner(SRWMatrix&A, double w,
                   std::string p_name){
        if (p_name == "SSOR"){
            SRWMatrix& D = *(new Eigen3Matrix(A.rows(), A.cols()));
            D = A.diag_m();
            this->Pr = &((D/w + A.lower_tri(true))
                    * (D.inverse()*w)
                    * (D/w + A.upper_tri(true)));
        }
    }

    SRWMatrix& P(){
        return *Pr;
    }

    SRWMatrix& iP(){
        return *iPr;
    }

protected:
    SRWMatrix* Pr;
    SRWMatrix* iPr;
};

#endif // PRECONDITIONER_H

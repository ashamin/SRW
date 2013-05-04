#ifndef PRECONDITIONERS_H
#define PRECONDITIONERS_H

#include "Headers/srwmatrix.h"
#include "Headers/eigen3matrix.h"

class Preconditioners
{
public:
    /**
     * @brief SSOR_precond
     *      Successive Symmetric Over Ralaxation
     * @param A
     * @param w
     * @return
     */
    static SRWMatrix& SSOR_precond(SRWMatrix& A, double w){
        SRWMatrix& D = *(new Eigen3Matrix(A.rows(), A.cols()));
        D = A.diag_m();
        return (D/w + A.lower_tri(true))
                * (D.inverse()*w)
                * (D/w + A.upper_tri(true));
    }
};

#endif // PRECONDITIONERS_H

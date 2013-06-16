#ifndef PAR2DPRECONDITIONER_H
#define PAR2DPRECONDITIONER_H

#include "preconditioner.h"

class par2DPreconditioner: public Preconditioner
{
public:

    par2DPreconditioner(double w, int n, double h,
                        std::string p_type): Preconditioner(){
        if (p_type == "par.SSOR"){
            Dx_p = new Eigen3Matrix(n, n);
            Dx_p->setZero();

            SRWVector& d = Dx_p->diag(0);
            d.fill((-4/(h*h)/w));
            Dx_p->setDiag(0, d);

            d = Dx_p->diag(1);
            d.fill(2/(h*h));
            Dx_p->setDiag(1, d);
            Dx_p->setDiag(-1, d);

            Dy_p = new Eigen3Matrix(n, n);
            Dy_p->setZero();

            d = Dy_p->diag(0);
            d.fill((-4/(h*h)/w));
            Dy_p->setDiag(0, d);

            d = Dy_p->diag(1);
            d.fill(2/(h*h));
            Dy_p->setDiag(1, d);
            Dy_p->setDiag(-1, d);

            SRWMatrix& tmp_m = *(new Eigen3Matrix(n, n));
            tmp_m.setZero();
            d = Dy_p->diag(0)/w;
            tmp_m.setDiag(0, d);

            Dy_p = &(tmp_m.inverse() * *Dy_p);

            Pr = &(*Dx_p * *Dy_p);
            iPr = &Pr->inverse();

        }
    }


    par2DPreconditioner(double w, int n, double h,
                        std::string p_type, bool inv): Preconditioner(){
        if (p_type == "par.SSOR"){
            Dx_p = new Eigen3Matrix(n, n);
            Dx_p->setZero();

            SRWVector& d = Dx_p->diag(0);
            d.fill((-4/(h*h)/w));
            Dx_p->setDiag(0, d);

            d = Dx_p->diag(1);
            d.fill(2/(h*h));
            Dx_p->setDiag(1, d);
            Dx_p->setDiag(-1, d);

            Dy_p = new Eigen3Matrix(n, n);
            Dy_p->setZero();

            d = Dy_p->diag(0);
            d.fill((-4/(h*h)/w));
            Dy_p->setDiag(0, d);

            d = Dy_p->diag(1);
            d.fill(2/(h*h));
            Dy_p->setDiag(1, d);
            Dy_p->setDiag(-1, d);

            SRWMatrix& tmp_m = *(new Eigen3Matrix(n, n));
            tmp_m.setZero();
            d = Dy_p->diag(0)/w;
            tmp_m.setDiag(0, d);

            Dy_p = &(tmp_m.inverse() * *Dy_p);

            Pr = &(*Dx_p * *Dy_p);
            if (inv) iPr = &Pr->inverse();

        }
    }

    SRWMatrix& Dx(){
        return *Dx_p;
    }

    SRWMatrix& Dy(){
        return *Dy_p;
    }

private:
    SRWMatrix* Dx_p;
    SRWMatrix* Dy_p;
};

#endif // PAR2DPRECONDITIONER_H


#include "Eigen/Dense"
//#include "Eigen/Eigen"
//#include "Eigen/src/Core/SolveTriangular.h"

#include <iostream>
#include "stdio.h"
#include "Headers/eigen3vector.h"
#include "Headers/srwmatrix.h"
#include "Headers/eigen3matrix.h"

#include "Headers/solvers.h"
#include "tests/Poisson2DSquareArea/Poisson2DSquareAreaTests"
#include "Headers/preconditioners.h"

//using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{

    /*MatrixXd m1 = MatrixXd::Random(5, 5);
    VectorXd v1 = VectorXd::Random(5);

    double x = m1(1, 1);

    VectorXd v2 = m1.ldlt().solve(v1);*/

    /*MatrixXd m2(2,2);

    m1(1, 1) = 1.8767;
    m2(1, 1) = 2.2384763;

    m1 =  m1*m2;*/


   /*SRWMatrix& m = *(new Eigen3Matrix(5, 5));

    SRWVector& v = m.diag(-1);

    v = *(new Eigen3Vector(5));

    bool olol = m.setDiag(0, v);

    if (olol) cout << "OK" << endl;

    cout << dynamic_cast<Eigen3Matrix&>(m).getMatrix() << endl << endl
         << dynamic_cast<Eigen3Vector&>(v).getVector() << endl << endl;

    m = m.transpose();

    cout << dynamic_cast<Eigen3Matrix&>(m).getMatrix() << endl << endl;

    SRWMatrix& im = m.inverse();

    cout << dynamic_cast<Eigen3Matrix&>(im * m).getMatrix() << endl << endl;*/

    /*int n = 1000;

    SRWVector& x = *(new Eigen3Vector(n));
    SRWMatrix& A = *(new Eigen3Matrix(n, n));
    for (int i = 0; i<n; i++)
        for (int j = 0; j<n; j++)
            A(i, j) = 0;
    SRWVector& b = *(new Eigen3Vector(n));

    A.setDiag(0, *(new Eigen3Vector(n)));
    A.setDiag(1, *(new Eigen3Vector(n)));
    A.setDiag(-1, *(new Eigen3Vector(n)));



    x = TDMA(A, x, b);

    (A*x-b).print();*/

    SRWMatrix& m = *(new Eigen3Matrix(5, 5));
    m.print();
    Preconditioners::SSOR_precond(m, 1).print();
    m.print();


    return 0;
}

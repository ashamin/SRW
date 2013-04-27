#include <QCoreApplication>
#include "Eigen/Dense"
//#include "Eigen/Eigen"
//#include "Eigen/src/Core/SolveTriangular.h"

#include <iostream>
#include "stdio.h"
#include "Headers/eigen3vector.h"

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{

    /*MatrixXd m1 = MatrixXd::Random(5, 5);
    VectorXd v1 = VectorXd::Random(5);

    double x = m1(1, 1);

    VectorXd v2 = m1.ldlt().solve(v1);

    /*MatrixXd m2(2,2);

    m1(1, 1) = 1.8767;
    m2(1, 1) = 2.2384763;

    m1 =  m1*m2;*/

    //cout << m1 << endl << endl << m1.transpose() << endl<< endl<< v2  << endl;

    cout << "olol" << endl;

    return 0;
}


#include "Eigen/Dense"
//#include "Eigen/Eigen"
//#include "Eigen/src/Core/SolveTriangular.h"

#include <iostream>
#include "stdio.h"
#include "Headers/eigen3vector.h"
#include "Headers/srwmatrix.h"
#include "Headers/eigen3matrix.h"

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

    //cout << m1 << endl << endl << m1.transpose() << endl<< endl<< v2  << endl;

    cout << "olol" << endl;

    /*Eigen3Vector* v1 = new Eigen3Vector(5);
    Eigen3Vector* v2 = new Eigen3Vector(1);

    for (int i =0; i<5; i++)
        v1->set(i, i*5);
    v2->set(0, 20);

    cout << v1->getVector() << endl << v1 << endl << endl;
    cout << v2->getVector() << endl << v2 << endl << endl;

    *v1 = *v2;

    v1->set(0, 10000);

    cout << v1->getVector() << endl << v1 << endl << endl;
    cout << v2->getVector() << endl << v2 << endl << endl;*/

    //SRWVector* v1 = new Eigen3Vector(10);

    //cout << ((Eigen3Vector*)v1)->getVector() << endl;

    /*Eigen3Matrix* m = new Eigen3Matrix(5, 5);
    Eigen3Matrix* m1 = new Eigen3Matrix(5, 5);
    Eigen3Matrix* m2 = new Eigen3Matrix(5, 5);

    cout << m->getMatrix() << endl;

    m = *m1 * m2;

    MatrixXd x = m1->getMatrix() * m2->getMatrix();

    if (m->getMatrix() == x) cout << "OLOLOLOLOLO" << endl;*/


    /*SRWMatrix& m = Eigen3Matrix(5, 5);

    static_cast<Eigen3Matrix>(m)(5, 5);*/


    //cout << m1 << endl;
   // cout << m.getMatrix() << endl;
   /* cout << m1->getMatrix() << endl;
    cout << m->getMatrix() << endl;*/

    SRWMatrix &x = *(new Eigen3Matrix(3, 3));

    cout << dynamic_cast<Eigen3Matrix&>(x).getMatrix();

    return 0;
}

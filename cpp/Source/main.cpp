
#include "Eigen/Dense"
//#include "Eigen/Eigen"
//#include "Eigen/src/Core/SolveTriangular.h"

#include <iostream>
#include "stdio.h"
#include "Headers/eigen3vector.h"
#include "Headers/srwmatrix.h"
#include "Headers/eigen3matrix.h"

#include "Headers/lsolvers.h"
#include "tests/Poisson2DSquareArea/Poisson2DSquareAreaTests"
#include "Headers/Preconditioning/preconditioners.h"

#include "Headers/solvers.h"

#include "Headers/Preconditioning/preconditioner.h"
#include "Headers/Preconditioning/par2DPreconditioner.h"


#include <fstream>
#include <iostream>
#include "time.h"

#include <vector>

//using namespace Eigen;
using namespace std;


void print_perf_measurement(int I, double norm, int itnum, double time){
    fstream out;
    out.open("test3_perf.csv", ios::in | ios::out);


    string str;

    out.seekg(0, ios::beg);
    while (out >> str){}

    out.clear();
    out.seekp(0, ios::end);

    out << I << "," << norm << "," << itnum << "," << time << endl;

    out.close();
}

void measure_poiss(int parts){
    int I , J;
    clock_t ctime;
    I = J = parts;
    Test2DPoissonSquareArea& test = *(new Test2DPoissonSquareAreaN1);

    int itnum = 5000;


    ctime = clock();
    SRWMatrix& sl = solve_poiss_2D_square(test, I, J, itnum);
    ctime = clock() - ctime;

    SRWMatrix& ra = test.get_right_answer_as_matrix(I, J);

    double time = (double)ctime/CLOCKS_PER_SEC;

    //print_perf_measurement(I, (ra - sl).norm("F")/sl.norm("F"), itnum, time);


    /*cout << "Frobenius norm of error" << "\t\t" << "iteration number" << endl;*/
    cout << time<< endl;
    cout << I << "\t\t\t" <<  (ra - sl).norm("F")/sl.norm("F") << "\t\t\t" << itnum << endl << endl;

}

#include "omp.h"
#define N 100000

int main(int argc, char *argv[])
{

    measure_poiss(30);

    /*SRWVector& v = Eigen3Vector(5);
    v.print();*/

    /*VectorXd a(10);

    a = a.transpose();*/


    return 0;
}

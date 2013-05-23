
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

    /*double* x = new double[10];

    for (int i = 0; i<10; i++)
        x[i] = i;

    for (int i = 0; i<5; i++){
        cout << (&x[0])[i] << " "  << (&x[5])[i] << endl;
    }*/


    /*int m = 9, n = 3;
    int k = 0, i = 0;
    double* dx_l = new double[m-1];
    double* dx_d = new double[m];
    double* dx_u = new double[m-1];
    double* r = new double[m];
    double* x1= new double[m];
    double* x = new double[m];

    for (int i = 0; i<(m-1); i++)
        dx_l[i] = dx_u[i] = i;
    dx_l[2] = dx_u[2] = 0;
    dx_l[5] = dx_u[5] = 0;
    for (int i = 0; i<m; i++){
        dx_d[i] = 20;
        r[i] = 7*(m-i+2);
    }

    TDMA_opt(dx_l, dx_d, dx_u, x1, r, m);

    #pragma omp parallel for shared(dx_d, dx_l, dx_u, r, x) \
                                firstprivate(n) private(i, k) \
                                schedule(static)
    for (i = 0; i<n; i++){
        k = i*n;
        TDMA_opt(&dx_l[k], &dx_d[k], &dx_u[k], &x[k], &r[k], n);
    }

    for (int i = 0; i<m; i++)
        cout << x1[i] << " " << x[i] << endl;*/


    measure_poiss(10);

    /*SRWVector& v = Eigen3Vector(5);
    v.print();*/

    /*VectorXd a(10);

    a = a.transpose();*/


    return 0;
}

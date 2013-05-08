
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
    Test2DPoissonSquareArea& test = *(new Test2DPoissonSquareAreaN3);

    int itnum = 5000;


    ctime = clock();
    SRWMatrix& sl = solve_poiss_2D_square(test, I, J, itnum);
    ctime = clock() - ctime;

    SRWMatrix& ra = test.get_right_answer_as_matrix(I, J);

    double time = (double)ctime/CLOCKS_PER_SEC;

    print_perf_measurement(I, (ra - sl).norm("F")/sl.norm("F"), itnum, time);


    /*cout << "Frobenius norm of error" << "\t\t" << "iteration number" << endl;*/
    cout << time<< endl;
    cout << I << "\t\t\t" <<  (ra - sl).norm("F")/sl.norm("F") << "\t\t\t" << itnum << endl << endl;

}

#include "omp.h"
#define N 100000

int main(int argc, char *argv[])
{

    //measure_poiss(20);

    /*par2DPreconditioner& p = *(new par2DPreconditioner(.8, 9, .11, "par.SSOR"));

    p.P().print();*/


      double a[N], b[N], c[N];
      int i;
      omp_set_dynamic(0);      // запретить библиотеке openmp менять число потоков во время исполнения
      omp_set_num_threads(10); // установить число потоков в 10

      // инициализируем массивы
      for (i = 0; i < N; i++)
      {
          a[i] = i * 1.0;
          b[i] = i * 2.0;
      }

    //clock_t ot1 = clock();
    double t1 = omp_get_wtime();

      // вычисляем сумму массивов
    #pragma omp parallel shared(a, b, c) private(i)
      {
    #pragma omp for
        for (i = 0; i < N; i++)
            for (int j = 0; j<1000; j++)
          c[i] = a[i] + b[i];
      }

   // clock_t ot2 = clock();
    double t2 = omp_get_wtime();

    //cout << (double)(ot2-ot1)/CLOCKS_PER_SEC << endl << endl;

    cout << t2-t1 << endl << endl;

      printf ("%f\n", c[10]);
      return 0;

    return 0;
}

#include <Eigen/Core>
#include "MatrixExponential.hpp"
#include <fstream>
#include <iostream>
#include <stdlib.h> /* atoi */
#include <sys/time.h>


#ifdef EIGEN_RUNTIME_NO_MALLOC
#define EIGEN_MALLOC_ALLOWED Eigen::internal::set_is_malloc_allowed(true);
#define EIGEN_MALLOC_NOT_ALLOWED Eigen::internal::set_is_malloc_allowed(false);
#else
#define EIGEN_MALLOC_ALLOWED
#define EIGEN_MALLOC_NOT_ALLOWED
#endif

using namespace Eigen;
using namespace std;
using namespace expokit;

#define N 4
#define N_TESTS 1000
#define N_RUNS 1

int main()
{
    cout << "Start test Matrix Exponential\n";

    ofstream myfile;
    myfile.open("/home/olli/Desktop/MAT-test both interfaces use return val PT4.txt");
    myfile << "testMatrixExpo- N_TESTS: " << N_TESTS << " N_RUNS: " << N_RUNS << " size N: " << N << "\n";

    Matrix<double, N, N> A;
    A << 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8;
    A *= 0.1;
    Matrix<double, N, 1> xInit;
    xInit << 1, 2, 3, 4;

    // MatrixExponential<double, N> expUtil;
    MatrixExponential<double, Dynamic> expUtil(N);

    Matrix<double, N, 1> res1;
    Matrix<double, N, N> res2;

    struct timeval stop, start;

    //EIGEN_MALLOC_NOT_ALLOWED

    cout << "Compute run number 7: ";
    for (int k = 0; k < N_RUNS; k++) {
        gettimeofday(&start, NULL);
        for (int i = 0; i < N_TESTS; i++) {
            res2 = expUtil.compute(A);
        }
        gettimeofday(&stop, NULL);
        printf("%i ", k);
        fflush(stdout);
        myfile << "compute(A, res2) took " << ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec) << " us\n";
    }

    cout << "\nComputeExpTimesVector run number: ";
    for (int k = 0; k < N_RUNS; k++) {
        gettimeofday(&start, NULL);
        for (int i = 0; i < N_TESTS; i++) {
            res1 = expUtil.computeExpTimesVector(A, xInit);
        }
        gettimeofday(&stop, NULL);
        printf("%i ", k);
        fflush(stdout);
        myfile << "computeExpTimesVector(A, xInit, res1) took " << ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec) << " us\n";
    }

    myfile.close();

    cout << "ExpTimesVector:---->\n"
         << res1 << "\n\n";
    cout << "Compute:---->\n"
         << res2 * xInit << "\n\n";
}

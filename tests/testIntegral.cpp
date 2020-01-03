#include "Eigen/Core"
#include <fstream>
#include <iostream>
#include <sys/time.h>

#include "LDSUtility.hpp"

using namespace expokit;
using namespace std;

#define N 4
#define N_TESTS 1000000
#define N_RUNS 10

int main()
{
    cout << "Start test integral\n";

    LDSUtility<double, N> test;

    ofstream myfile;
    myfile.open("/home/olli/Desktop/INT1b5b5f56e - test return params.txt");
    myfile << "testIntegral - N_TESTS: " << N_TESTS << " N_RUNS: " << N_RUNS << "\n";

    Matrix<double, N, N> A;
    A << 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8;
    A *= 0.1;
    Matrix<double, N, 1> xInit;
    xInit << 1, 2, 3, 4;
    Matrix<double, N, 1> b;
    b << 4, 3, 2, 1;

    Matrix<double, N, 1> res;

    struct timeval stop, start;

    for (int k = 0; k < N_RUNS; k++) {
        gettimeofday(&start, NULL);
        for (int i = 0; i < N_TESTS; i++) {
            test.ComputeXt(A, b, xInit, 1, res);
        }
        gettimeofday(&stop, NULL);
        myfile << "ComputeXt(A, b, xInit, 1, res) took " << ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec) << " us\n";
    }
    cout << "x(T):\n"
         << res;

    for (int k = 0; k < N_RUNS; k++) {
        gettimeofday(&start, NULL);
        for (int i = 0; i < N_TESTS; i++) {
            test.ComputeIntegralXt(A, b, xInit, 1, res);
        }
        gettimeofday(&stop, NULL);
        myfile << "ComputeIntegralXt(A, b, xInit, 1, res) took " << ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec) << " us\n";
    }
    cout << "\n\nxint(T):\n"
         << res;

    for (int k = 0; k < N_RUNS; k++) {
        gettimeofday(&start, NULL);
        for (int i = 0; i < N_TESTS; i++) {
            test.ComputeDoubleIntegralXt(A, b, xInit, 1, res);
        }
        gettimeofday(&stop, NULL);
        myfile << "ComputeDoubleIntegralXt(A, b, xInit, 1, res) took " << ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec) << " us\n";
    }
    cout << "\n\nxintint(T):\n"
         << res;

    cout << "\n";
    return 0;
}
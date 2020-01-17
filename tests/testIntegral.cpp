#include "Eigen/Core"
#include "utils/stop-watch.h"
#include <fstream>
#include <iostream>
#include <sys/time.h>

#include "LDSUtility.hpp"

#ifdef EIGEN_RUNTIME_NO_MALLOC
#define EIGEN_MALLOC_ALLOWED Eigen::internal::set_is_malloc_allowed(true);
#define EIGEN_MALLOC_NOT_ALLOWED Eigen::internal::set_is_malloc_allowed(false);
#else
#define EIGEN_MALLOC_ALLOWED
#define EIGEN_MALLOC_NOT_ALLOWED
#endif

using namespace expokit;
using namespace std;

#define N 4
#define M N * 3 * 2
#define N_TESTS 1000
#define N_RUNS 1

int main()
{
    cout << "Start test integral\n";

    LDSUtility<double, M> test; // Static
    //LDSUtility<double, Dynamic> test(M); // Dynamic

    ofstream myfile;
    myfile.open("/home/olli/Desktop/INT-test with bigger mat PT4.txt");
    myfile << "testIntegral - N_TESTS: " << N_TESTS << " N_RUNS: " << N_RUNS << " size N: " << N << "\n";

    int m2 = int(M / 2);
    double stiffness = 1e5;
    double damping = 1e2;
    MatrixXd U = MatrixXd::Random(m2, m2);
    MatrixXd Upsilon = U * U.transpose();
    MatrixXd K = MatrixXd::Identity(m2, m2) * stiffness;
    MatrixXd B = MatrixXd::Identity(m2, m2) * damping;
    MatrixXd A = MatrixXd::Zero(M, M);
    A.topRightCorner(m2, m2) = MatrixXd::Identity(m2, m2);
    A.bottomLeftCorner(m2, m2) = -Upsilon * K;
    A.bottomRightCorner(m2, m2) = -Upsilon * B;
    MatrixXd xInit = MatrixXd::Random(M, 1);
    MatrixXd b = MatrixXd::Random(M, 1);

    Matrix<double, M, 1> res;

    struct timeval stop, start;

    EIGEN_MALLOC_NOT_ALLOWED

    cout << "ComputeXt run number: ";
    for (int k = 0; k < N_RUNS; k++) {
        gettimeofday(&start, NULL);
        for (int i = 0; i < N_TESTS; i++) {
            test.ComputeXt(A, b, xInit, 1, res);
        }
        gettimeofday(&stop, NULL);
        printf("%i ", k);
        fflush(stdout);
        myfile << "ComputeXt(A, b, xInit, 1, res) took " << ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec) << " us\n";
    }
    cout << "\nx(T):\n"
         << res;

    cout << "\n\nComputeIntegralXt run number: ";
    for (int k = 0; k < N_RUNS; k++) {
        gettimeofday(&start, NULL);
        for (int i = 0; i < N_TESTS; i++) {
            test.ComputeIntegralXt(A, b, xInit, 1, res);
        }
        gettimeofday(&stop, NULL);
        printf("%i ", k);
        fflush(stdout);
        myfile << "ComputeIntegralXt(A, b, xInit, 1, res) took " << ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec) << " us\n";
    }
    cout << "\nxint(T):\n"
         << res;

    cout << "\n\nComputeDoubleIntegralXt run number: ";
    for (int k = 0; k < N_RUNS; k++) {
        gettimeofday(&start, NULL);
        for (int i = 0; i < N_TESTS; i++) {
            test.ComputeDoubleIntegralXt(A, b, xInit, 1, res);
        }
        gettimeofday(&stop, NULL);
        printf("%i ", k);
        fflush(stdout);
        myfile << "ComputeDoubleIntegralXt(A, b, xInit, 1, res) took " << ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec) << " us\n";
    }
    cout << "\nxintint(T):\n"
         << res << "\n";

    getProfiler().report_all(3);

    return 0;
}
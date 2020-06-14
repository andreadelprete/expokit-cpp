#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "BalancingMethods.hpp"
#include "MatrixExponential.hpp"

#include "utils/errorComputing.hpp"
#include "utils/readTSV.hpp"
#include "utils/statistics.hpp"
#include "utils/stop-watch.h"

using namespace std;
using namespace expokit;
using namespace Eigen;

#ifdef EIGEN_RUNTIME_NO_MALLOC
#define EIGEN_MALLOC_ALLOWED Eigen::internal::set_is_malloc_allowed(true);
#define EIGEN_MALLOC_NOT_ALLOWED Eigen::internal::set_is_malloc_allowed(false);
#else
#define EIGEN_MALLOC_ALLOWED
#define EIGEN_MALLOC_NOT_ALLOWED
#endif

#define N 64
#define TESTS 10000

Matrix<double, N, N> james2014Generator();

int main()
{
    cout << "Running new balance" << endl;

    MatrixXd A(N, N);
    MatrixXd B(N, N);
    MatrixXd out(N, N);
    MatrixXd D(N, N);
    MatrixXd Dinv(N, N);

    BalancingMethods<double, Dynamic> util(N);
    MatrixExponential<double, Dynamic> expUtil(N);

    for (int i = 0; i < TESTS; ++i) {
        A = james2014Generator();

        START_PROFILER("EXPnotBalanced");
        expUtil.compute(A, out);
        STOP_PROFILER("EXPnotBalanced");

        START_PROFILER("NB");
        util.balanceNew(A, B, D, Dinv);
        STOP_PROFILER("NB");

        START_PROFILER("EXPNBBalanced");
        expUtil.compute(B, out);
        STOP_PROFILER("EXPNBBalanced");

        START_PROFILER("Rodney");
        util.balanceRodney(A, B, D, Dinv);
        STOP_PROFILER("Rodney");

        START_PROFILER("EXPRodneyBalanced");
        expUtil.compute(B, out);
        STOP_PROFILER("EXPRodneyBalanced");
    }

    // Print out results
    getProfiler().report_all(3);
}

// Generate test matrices
Matrix<double, N, N> james2014Generator()
{
    Matrix<double, N, N> A = Matrix<double, N, N>::Random();
    Matrix<double, N, N> D = Matrix<double, N, N>::Zero();
    Matrix<double, N, N> Dinv = Matrix<double, N, N>::Zero();

    for (int i = 0; i < N; i++) {
        unsigned int e = rand() % 20;
        unsigned int twoPow = 1U << e;

        D(i, i) = twoPow;
        Dinv(i, i) = 1.0 / twoPow;
    }

    return Dinv * A * D;
}

/*
    IOFormat CleanFmt(FullPrecision, 0, ", ", "\n", "[", "]");

    std::cout << "A:" << std::endl
              << A << std::endl
              << "it:" << it << std::endl
              << "B:" << std::endl
              << B << std::endl
              << "Norm of B:" << std::endl
              << B.lpNorm<1>() << std::endl
              << "D:" << std::endl
              << D << std::endl
              << "Dinv:" << std::endl
              << Dinv << std::endl
              << "Check:" << std::endl
              << (A - D * B * Dinv).eval() << std::endl;
*/
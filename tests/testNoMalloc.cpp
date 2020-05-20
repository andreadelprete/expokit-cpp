#include "LDSUtility.hpp"
#include "MatrixExponential.hpp"
#include "BalancingMethods.hpp"

#include <Eigen/Core>
#include <iostream>

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

int main()
{
    cout << "Start test No Malloc\n";

    Matrix<double, N, N> A;
    A << 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8;
    A *= 0.1;
    Matrix<double, N, 1> xInit;
    xInit << 1, 2, 3, 4;
    Matrix<double, N, 1> b;
    b << 4, 3, 2, 1;

    MatrixXd B(N, N);
    MatrixXd D(N, N);
    MatrixXd Dinv(N, N);

    MatrixExponential<double, Dynamic> expUtil(N);
    LDSUtility<double, Dynamic> lds(N);

    // MatrixExponential
    Matrix<double, N, 1> res1;
    Matrix<double, N, N> res2;
    BalancingMethods<double, Dynamic> balUtil(N);


    EIGEN_MALLOC_NOT_ALLOWED
    expUtil.computeExpTimesVector(A, xInit, res1);
    expUtil.compute(A, res2);

    // Balance methods
    balUtil.balanceNew(A, B, D, Dinv);
    balUtil.balanceRodney(A, B, D, Dinv);

    // LDSUtility
    lds.ComputeXt(A, b, xInit, 1, res1);
    lds.ComputeIntegralXt(A, b, xInit, 1, res1);
    lds.ComputeDoubleIntegralXt(A, b, xInit, 1, res1);
}
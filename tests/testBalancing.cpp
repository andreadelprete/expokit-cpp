#include <Eigen/Core>
#include <iostream>

#include "BalancingMethods.hpp"

using namespace expokit;
using namespace Eigen;

#define N 4

int main()
{
    MatrixXd A(N, N);
    // A << 1, 2, 300, 4, 5, 6, 0, 0, 0;
    A << -3.80466072e-01, 1.19867645e+06, 8.63170116e-02,
        -2.60920215e+05,
        1.68822169e+00, 1.56173278e+05, 4.61874965e+05,
        -9.96831380e+05,
        1.20856994e+00, -8.61700547e-01, -6.34023305e-01,
        -9.35935672e-01,
        -2.87160379e+05, 4.70987189e+05, -7.37804710e-01,
        6.75828111e-01;
    MatrixXd B(N, N);
    VectorXd D(N);
    VectorXd Dinv(N);

    BalancingMethods<double, Dynamic> util(N);

    util.balanceNew(A, B, D, Dinv);
    if ((A - D.asDiagonal() * B * Dinv.asDiagonal()).eval().norm() != 0)
        return -1;

    util.balanceRodney(A, B, D, Dinv);
    if ((A - D.asDiagonal() * B * Dinv.asDiagonal()).eval().norm() != 0)
        return -1;

    // TODO test for actual NB correctness, get results from python

    return 0;
}

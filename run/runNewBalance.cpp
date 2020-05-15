#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "MatrixExponential.hpp"

#include "utils/errorComputing.hpp"
#include "utils/readTSV.hpp"
#include "utils/statistics.hpp"
#include "utils/stop-watch.h"

using namespace std;
using namespace expokit;
using namespace Eigen;

#define N 4

int main()
{
    cout << "Running speed ups benchmark" << endl;

    IOFormat CleanFmt(FullPrecision, 0, ", ", "\n", "[", "]");

    MatrixXd A(N, N);
    // A << 1, 2, 300, 4, 5, 6, 0, 0, 0;
    A <<-3.80466072e-01,  1.19867645e+06,  8.63170116e-02,
        -2.60920215e+05,
        1.68822169e+00,  1.56173278e+05,  4.61874965e+05,
        -9.96831380e+05,
        1.20856994e+00, -8.61700547e-01, -6.34023305e-01,
        -9.35935672e-01,
       -2.87160379e+05,  4.70987189e+05, -7.37804710e-01,
         6.75828111e-01;
    MatrixXd B(N, N);
    MatrixXd D(N, N);
    MatrixXd Dinv(N, N);

    MatrixExponential<double, Dynamic> util(N);

    int it = util.newBalancing(A, B, D, Dinv);

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
}
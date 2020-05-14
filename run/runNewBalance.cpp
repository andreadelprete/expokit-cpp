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

#define N 3
#define M N * 3 * 2

int main()
{
    cout << "Running speed ups benchmark" << endl;

    MatrixXd A(N, N);
    A << 1, 2, 300, 4, 5, 6, 0, 0, 0;
    MatrixXd B(N, N);
    MatrixXd D(N, N);
    MatrixXd Dinv(N, N);

    MatrixExponential<double, Dynamic> util(3);

    int it = util.newBalancing(A, B, D, Dinv);

    std::cout << "A:" << std::endl
              << A << std::endl
              << "it:" << it << std::endl
              << "B:" << std::endl
              << B << std::endl
              << "D:" << std::endl
              << D << std::endl
              << "Dinv:" << std::endl
              << Dinv << std::endl
              << "Check:" << std::endl
              << (A - D * B * Dinv).eval() << std::endl;
}
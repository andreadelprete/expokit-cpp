#include <Eigen/Core>
#include "MatrixExponential.hpp"
#include <iostream>

#include <stdlib.h>     /* atoi */

using namespace Eigen;
using namespace std;
using namespace expokit;

#define N 3

int main()
{
     cout << "Start test Matrix Exponential\n";

     Matrix<double, N, N> A;
     A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
     Matrix<double, N, 1> xInit;
     xInit << 1, 2, 3;

     MatrixExponential<double, N> expUtil;

     Matrix<double, N, 1> res1;
     Matrix<double, N, N> res2;

     expUtil.computeExpTimesVector(A, xInit, res1);
     expUtil.compute(A, res2);

     cout << "ExpTimesVector:---->\n"
          << res1 << "\n\n";
     cout << "Compute:---->\n"
          << res2 * xInit << "\n\n";
}

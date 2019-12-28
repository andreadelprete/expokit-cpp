#include <Eigen/Core>
#include "MatrixExponential.h"

#include <Eigen/LU>
#include <iostream>

#include <stdlib.h>     /* atoi */

using namespace Eigen;
using namespace std;
// using namespace expokit;

#define N 3

int main()
{
    // Matrix<double, N, N> A = Matrix<double, N, N>::Random();
    // Matrix<double, N, 1> xInit = Matrix<double, N, 1>::Random();

    Matrix<double, N, N> A;
    A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    Matrix<double, N, 1> xInit;
    xInit << 1, 2, 3;

    MatrixExponential<Matrix<double, N, N>, Matrix<double, N, 1> > expUtil(N);

    Matrix<double, N, 1> res1;
    Matrix<double, N, N> res2;

    expUtil.computeExpTimesVector(A, xInit, res1);
    expUtil.compute(A, res2);

    cout << "ExpTimesVector:---->\n"
         << res1 << "\n\n";
    cout << "compute:---->\n"
         << res2 * xInit << "\n\n";
}

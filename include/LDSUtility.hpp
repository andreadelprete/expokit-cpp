#ifndef INTEGRAL_UTILITY
#define INTEGRAL_UTILITY

#include <Eigen/Core>
#include "MatrixExponential.h"
#include <iostream>
#include <Eigen/LU>

using namespace Eigen;
using namespace std;

namespace expokit {

/*
    Utility class that allows integrating a system with the form
    Ax + b using a matrix exponential method 
*/
template <typename type, int N>
class LDSUtility {
private:
    MatrixExponential<Matrix<type, N + 1, N + 1>, Matrix<type, N + 1, 1>> expUtil;

public:
    LDSUtility();
    ~LDSUtility();

    typedef const Eigen::Ref<const Matrix<type, N, 1>> RefVector;
    typedef const Eigen::Ref<const Matrix<type, N, N>> RefMatrix;

    /**
     * Compute the value of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    Matrix<type, N, 1> ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, type T);

    /**
     * Compute the value of the integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    Matrix<type, N, 1> ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, type T);

    /**
     * Compute the value of the double integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    Matrix<type, N, 1> ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, type T);
};

template <typename type, int N>
LDSUtility<type, N>::LDSUtility()
{
    expUtil = MatrixExponential<Matrix<type, N + 1, N + 1>, Matrix<type, N + 1, 1>> (N + 1);
}

template <typename type, int N>
LDSUtility<type, N>::~LDSUtility()
{
}

template <typename type, int N>
Matrix<type, N, 1> LDSUtility<type, N>::ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, type T)
{
    // Building aumented matrix A0
    Matrix<type, N + 1, N + 1> A0 = Matrix<type, N + 1, N + 1>::Zero();
    A0.template block<N, N>(0, 0) = A;
    A0.template block<N, 1>(0, N) = b;
    A0 *= T;

    // Building augmented state x0
    Matrix<type, N + 1, 1> x0;
    x0 << xInit, 1;

    // Building filtering matrix z0
    // TODO - couldn't we just extract a block from the result?
    Matrix<type, N, N + 1> z0;
    z0 << Matrix<type, N, N>::Identity(), Matrix<type, N, 1>::Zero();

    // Temp matrix
    Matrix<type, N + 1, 1> xTemp;

    // Matrix exponential
    expUtil.computeExpTimesVector(A0, x0, xTemp);

    // Extracting the interesting result
    return xTemp.template block<N, 1>(0, 0);
}
}

#endif
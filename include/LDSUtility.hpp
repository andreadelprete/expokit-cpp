#ifndef INTEGRAL_UTILITY
#define INTEGRAL_UTILITY

#include <Eigen/Core>
#include "MatrixExponential.hpp"
#include <iostream>

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
    MatrixExponential<type, N + 1> expUtil1;
    MatrixExponential<type, N + 2> expUtil2;
    MatrixExponential<type, N + 3> expUtil3;

public:
    LDSUtility();
    ~LDSUtility();

    typedef const Eigen::Ref<const Matrix<type, N, 1> > RefVector;
    typedef const Eigen::Ref<const Matrix<type, N, N> > RefMatrix;

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

    // Temp matrix
    Matrix<type, N + 1, 1> xTemp;

    // Matrix exponential
    expUtil1.computeExpTimesVector(A0, x0, xTemp);

    // Extracting the interesting result
    return xTemp.template block<N, 1>(0, 0);
}

template <typename type, int N>
Matrix<type, N, 1> LDSUtility<type, N>::ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, type T) {
    // Building augmented state x0
    Matrix<type, N + 1, 1> x1;
    x1 << xInit, 1;

    // Building aumented matrix A1
    Matrix<type, N + 2, N + 2> A1 = Matrix<type, N + 2, N + 2>::Zero();
    A1.template block<N, N>(0, 0) = A;
    A1.template block<N, 1>(0, N) = b;
    A1.template block<N + 1, 1>(0, N + 1) = x1;
    A1 *= T;

    // Vector z
    Matrix<type, N + 2, 1> z = Matrix<type, N + 2, 1>::Zero();
    z(N + 1, 0) = 1;

    // Temp matrix
    Matrix<type, N + 2, 1> xTemp;

    // Matrix exponential
    expUtil2.computeExpTimesVector(A1, z, xTemp);

    // Extracting the interesting result
    return xTemp.template block<N, 1>(0, 0);
}

template <typename type, int N>
Matrix<type, N, 1> LDSUtility<type, N>::ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, type T) {
    // Building aumented matrix A2
    Matrix<type, N + 3, N + 3> A2 = Matrix<type, N + 3, N + 3>::Zero();
    A2.template block<N, N>(0, 0) = A;
    A2.template block<N, 1>(0, N) = b;
    A2.template block<N, 1>(0, N + 1) = xInit;
    A2.template block<2, 2>(N, N + 1) = Matrix<type, 2, 2>::Identity();
    A2 *= T;

    // Vector z
    Matrix<type, N + 3, 1> z = Matrix<type, N + 3, 1>::Zero();
    z(N + 2, 0) = 1;

    // Temp matrix
    Matrix<type, N + 3, 1> xTemp;

    // Matrix exponential
    expUtil3.computeExpTimesVector(A2, z, xTemp);

    // Extracting the interesting result
    return xTemp.template block<N, 1>(0, 0);
}

}

#endif
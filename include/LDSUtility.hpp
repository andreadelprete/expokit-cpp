#ifndef INTEGRAL_UTILITY
#define INTEGRAL_UTILITY

#include <Eigen/Core>
#include <iostream>
#include "MatrixExponential.hpp"

using namespace Eigen;
using namespace std;

namespace expokit {

/*
    Utility class that allows integrating a system with the form
    Ax + b using a matrix exponential method 
*/
template <typename T, int N>
class LDSUtility {
private:
    MatrixExponential<T, N + 1> expUtil1;
    MatrixExponential<T, N + 2> expUtil2;
    MatrixExponential<T, N + 3> expUtil3;

    Matrix<T, N + 1, 1> res1;
    Matrix<T, N + 2, 1> res2;
    Matrix<T, N + 3, 1> res3;

public:
    LDSUtility() {}

    typedef const Ref<const Matrix<T, N, 1>> RefVector;
    typedef const Ref<const Matrix<T, N, N>> RefMatrix;

    typedef Ref<Matrix<T, N, 1>> RefOutVector;

    /**
     * Compute the value of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the double integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);
};

template <typename T, int N>
void LDSUtility<T, N>::ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    // Building aumented matrix A0
    Matrix<T, N + 1, N + 1> A0 = Matrix<T, N + 1, N + 1>::Zero();
    A0.template block<N, N>(0, 0) = A;
    A0.template block<N, 1>(0, N) = b;
    A0 *= step;

    // Building augmented state x0
    Matrix<T, N + 1, 1> x0;
    x0 << xInit, 1;

    // Matrix exponential Extracting the interesting result
    expUtil1.computeExpTimesVector(A0, x0, res1);
    out = res1.template block<N, 1>(0, 0);
}

template <typename T, int N>
void LDSUtility<T, N>::ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    // Building augmented state x0
    Matrix<T, N + 1, 1> x1;
    x1 << xInit, 1;

    // Building aumented matrix A1
    Matrix<T, N + 2, N + 2> A1 = Matrix<T, N + 2, N + 2>::Zero();
    A1.template block<N, N>(0, 0) = A;
    A1.template block<N, 1>(0, N) = b;
    A1.template block<N + 1, 1>(0, N + 1) = x1;
    A1 *= step;

    // Vector z
    Matrix<T, N + 2, 1> z = Matrix<T, N + 2, 1>::Zero();
    z(N + 1, 0) = 1;

    // Matrix exponential and extracting the interesting result
    expUtil2.computeExpTimesVector(A1, z, res2);
    out = res2.template block<N, 1>(0, 0);
}

template <typename T, int N>
void LDSUtility<T, N>::ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    // Building aumented matrix A2
    Matrix<T, N + 3, N + 3> A2 = Matrix<T, N + 3, N + 3>::Zero();
    A2.template block<N, N>(0, 0) = A;
    A2.template block<N, 1>(0, N) = b;
    A2.template block<N, 1>(0, N + 1) = xInit;
    A2.template block<2, 2>(N, N + 1) = Matrix<T, 2, 2>::Identity();
    A2 *= step;

    // Vector z
    Matrix<T, N + 3, 1> z = Matrix<T, N + 3, 1>::Zero();
    z(N + 2, 0) = 1;

    // Matrix exponential and extracting the interesting result
    expUtil3.computeExpTimesVector(A2, z, res3);
    out =  res3.template block<N, 1>(0, 0);
}

/*
        _____                              _      
        |  __ \                            (_)     
        | |  | |_   _ _ __   __ _ _ __ ___  _  ___ 
        | |  | | | | | '_ \ / _` | '_ ` _ \| |/ __|
        | |__| | |_| | | | | (_| | | | | | | | (__ 
        |_____/ \__, |_| |_|\__,_|_| |_| |_|_|\___|
                __/ |                             
                |___/                              
*/

template <typename T>
class LDSUtility<T, Dynamic> {
private:
    MatrixExponential<T, Dynamic> expUtil1;
    MatrixExponential<T, Dynamic> expUtil2;
    MatrixExponential<T, Dynamic> expUtil3;

    int n;

    typedef Matrix<T, Dynamic, 1> DynVector;
    typedef Matrix<T, Dynamic, Dynamic> DynMatrix;

    // Preallocating useful stuff
    DynVector res1, res2, res3, x0, z1, z2;
    DynMatrix A0, A1, A2;
    

public:
    // Would like to forbid creation without specifying a size, but it would prevent use as a class field
    LDSUtility()
        : LDSUtility(2)
    {
    } // Creating with dim 2
    LDSUtility(int n);

    void resize(int n);

    typedef const Ref<const DynVector> RefVector;
    typedef const Ref<const DynMatrix> RefMatrix;

    typedef Ref<DynVector> RefOutVector;

    /**
     * Compute the value of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the double integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);
};

template <typename T>
LDSUtility<T, Dynamic>::LDSUtility(int n)
{
    resize(n);
}

template <typename T>
void LDSUtility<T, Dynamic>::resize(int n)
{
    this->n = n;
    expUtil1.resize(n + 1);
    expUtil2.resize(n + 2);
    expUtil3.resize(n + 3);

    // Stuff for ComputeXt
    res1.resize(n + 1, 1);
    A0 = DynMatrix::Zero(n + 1, n + 1);
    x0.resize(n + 1, 1);

    // Stuff for ComputeIntegralXt
    res2.resize(n + 2, 1);
    A1 = DynMatrix::Zero(n + 2, n + 2);
    z1 = DynVector::Zero(n + 2, 1);
    z1(n + 1, 0) = 1; // This is constant

    // Stuff for ComputeDoubleIntegralXt
    res3.resize(n + 3, 1);
    A2 = DynMatrix::Zero(n + 3, n + 3);
    z2 = DynVector::Zero(n + 3, 1);
    z2(n + 2, 0) = 1;
}

template <typename T>
void LDSUtility<T, Dynamic>::ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    // Building aumented matrix A0
    A0.block(0, 0, n, n) = A; 
    A0.block(0, n, n, 1) = b;
    A0 *= step;

    // Building augmented state x0
    x0 << xInit, 1;

    // Matrix exponential Extracting the interesting result
    expUtil1.computeExpTimesVector(A0, x0, res1);
    out = res1.block(0, 0, n, 1);
}

template <typename T>
void LDSUtility<T, Dynamic>::ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    // Building augmented state x0
    x0 << xInit, 1;

    // Building aumented matrix A1
    A1.block(0, 0, n, n) = A;
    A1.block(0, n, n, 1) = b;
    A1.block(0, n + 1, n + 1, 1) = x0;
    A1 *= step;

    //Matrix exponential and extracting the interesting result
    expUtil2.computeExpTimesVector(A1, z1, res2);
    out = res2.block(0, 0, n, 1);
}

template <typename T>
void LDSUtility<T, Dynamic>::ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    // Building aumented matrix A2
    A2.block(0, 0, n, n) = A;
    A2.block(0, n, n, 1) = b;
    A2.block(0, n + 1, n, 1) = xInit;
    A2.block(n, n + 1, 2, 2) = Matrix<T, 2, 2>::Identity();
    A2 *= step;

    // Matrix exponential and extracting the interesting result
    expUtil3.computeExpTimesVector(A2, z2, res3);
    out = res3.block(0, 0, n, 1);
}
}

#endif
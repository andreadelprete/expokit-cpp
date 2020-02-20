#ifndef INTEGRAL_UTILITY2
#define INTEGRAL_UTILITY2

#include <Eigen/Core>
#include <iostream>
#include "LDSUtility.hpp"


using namespace Eigen;
using namespace std;

namespace expokit {

/*
    Wrapper around LDSUtility used to test timescaling optimization
    Eventually it may be extended to a generic second order system intergrator

    TODO:
    - Dynamic class skeleton
    - Test for dynamic allocation
    - Static implementation
    - Generic structure for 2 order
*/

template <typename T, int N>
class LDS2OrderUtility {
private:
    MatrixExponential<T, N + 1> expUtil1;
    MatrixExponential<T, N + 2> expUtil2;
    MatrixExponential<T, N + 3> expUtil3;

    Matrix<T, N + 1, 1> res1;
    Matrix<T, N + 2, 1> res2;
    Matrix<T, N + 3, 1> res3;

    Matrix<T, N, N> As;

    Matrix<T, N + 1, N + 1> A0;
    Matrix<T, N + 1, 1> x0;

    Matrix<T, N + 2, N + 2> A1;

    Matrix<T, N + 2, 1> z1;

    Matrix<T, N + 3, N + 3> A2;

    Matrix<T, N + 3, 1> z2;

public:
    LDS2OrderUtility();

    typedef const Ref<const Matrix<T, N, 1>> RefVector;
    typedef const Ref<const Matrix<T, N, N>> RefMatrix;
    typedef const Ref<const Matrix<T, N / 2, N / 2>> RefSubMatrix;

    typedef Ref<Matrix<T, N, 1>> RefOutVector;

    /**
     * Compute the value of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);
    // void ComputeXt(RefSubMatrix& Kbar, RefSubMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);
    // void ComputeIntegralXt(RefSubMatrix& Kbar, RefSubMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the double integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);
    // void ComputeDoubleIntegralXt(RefSubMatrix& Kbar, RefSubMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out);
};

template <typename T, int N>
LDS2OrderUtility<T, N>::LDS2OrderUtility()
{
    A0 = Matrix<T, N + 1, N + 1>::Zero();
    A1 = Matrix<T, N + 2, N + 2>::Zero();
    z1 = Matrix<T, N + 2, 1>::Zero();
    A2 = Matrix<T, N + 3, N + 3>::Zero();
    z2 = Matrix<T, N + 3, 1>::Zero();
}

template <typename T, int N>
void LDS2OrderUtility<T, N>::ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    // Building aumented matrix A0
    A0.template block<N, N>(0, 0) = A;
    A0.template block<N, 1>(0, N) = b;
    A0 *= step;

    // Building augmented state x0
    x0 << xInit, 1;

    // Matrix exponential Extracting the interesting result
    expUtil1.computeExpTimesVector(A0, x0, res1);
    out = res1.template block<N, 1>(0, 0);
}

template <typename T, int N>
void LDS2OrderUtility<T, N>::ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    // Building augmented state x0
    x0 << xInit, 1;

    // Building aumented matrix A1
    A1.template block<N, N>(0, 0) = A;
    A1.template block<N, 1>(0, N) = b;
    A1.template block<N + 1, 1>(0, N + 1) = x0;
    A1 *= step;

    // Vector z
    z1(N + 1, 0) = 1;

    // Matrix exponential and extracting the interesting result
    expUtil2.computeExpTimesVector(A1, z1, res2);
    out = res2.template block<N, 1>(0, 0);
}

template <typename T, int N>
void LDS2OrderUtility<T, N>::ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    // Building aumented matrix A2
    A2.template block<N, N>(0, 0) = A;
    A2.template block<N, 1>(0, N) = b;
    A2.template block<N, 1>(0, N + 1) = xInit;
    A2.template block<2, 2>(N, N + 1) = Matrix<T, 2, 2>::Identity();
    A2 *= step;

    // Vector z
    z2(N + 2, 0) = 1;

    // Matrix exponential and extracting the interesting result
    expUtil3.computeExpTimesVector(A2, z2, res3);
    out = res3.template block<N, 1>(0, 0);
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
class LDS2OrderUtility<T, Dynamic> {
private:
    // MatrixExponential<T, Dynamic> expUtil1;
    // MatrixExponential<T, Dynamic> expUtil2;
    // MatrixExponential<T, Dynamic> expUtil3;
    LDSUtility<T, Dynamic> firstOrder;

    int n, nHalf;
    T scalingT; // Just for debug purposes - should be local variable
    T scalingSqrd;

    typedef Matrix<T, Dynamic, 1> DynVector;
    typedef Matrix<T, Dynamic, Dynamic> DynMatrix;

    // Preallocating useful stuff
    DynVector res1, res2, res3, x0, z1, z2, mulXInit;
    DynMatrix As, A0, A1, A2;

public:
    // Would like to forbid creation without specifying a size, but it would prevent use as a class field
    LDS2OrderUtility()
        : LDS2OrderUtility(2)
    {
    } // Creating with dim 2
    LDS2OrderUtility(int n);

    void resize(int n);

    typedef const Ref<const DynVector> RefVector;
    typedef const Ref<const DynMatrix> RefMatrix;

    typedef Ref<DynVector> RefOutVector;

    /**
     * Compute the value of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeXt(RefMatrix& Kbar, RefMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeIntegralXt(RefMatrix& Kbar, RefMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the double integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeDoubleIntegralXt(RefMatrix& Kbar, RefMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out);

private:
    void setScaling(RefMatrix& Kbar);
    void composeA(RefMatrix& Kbar, RefMatrix& Bbar);
};

template <typename T>
LDS2OrderUtility<T, Dynamic>::LDS2OrderUtility(int n)
{
    resize(n);
}

template <typename T>
void LDS2OrderUtility<T, Dynamic>::resize(int n)
{
    this->n = n;
    this->nHalf = n / 2;
    firstOrder.resize(n);

    // Common
    As = DynMatrix::Zero(n, n);
    As.block(0, nHalf, nHalf, nHalf) = DynMatrix::Identity(nHalf, nHalf);

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

// Just to avoid repetition
template <typename T>
inline void LDS2OrderUtility<T, Dynamic>::setScaling(RefMatrix& Kbar)
{
    const double l1norm = Kbar.cwiseAbs().colwise().sum().maxCoeff();
    scalingT = sqrt(1 / l1norm);
    scalingSqrd = scalingT * scalingT;

    mulXInit = DynVector::Ones(n, 1);
    mulXInit.block(nHalf, 0, nHalf, 1) = scalingT * mulXInit.block(nHalf, 0, nHalf, 1);
}

template <typename T>
inline void LDS2OrderUtility<T, Dynamic>::composeA(RefMatrix& Kbar, RefMatrix& Bbar)
{
    As.block(nHalf, nHalf, nHalf, nHalf) = Bbar * -scalingT;
    As.block(nHalf, 0, nHalf, nHalf) = Kbar * -scalingSqrd;
}

// ComputeXt
template <typename T>
void LDS2OrderUtility<T, Dynamic>::ComputeXt(RefMatrix& Kbar, RefMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    setScaling(Kbar);

    composeA(Kbar, Bbar);

    firstOrder.ComputeXt(As, b * scalingSqrd, xInit.cwiseProduct(mulXInit), step / scalingT, out);
}

// ComputeIntegralXt
template <typename T>
void LDS2OrderUtility<T, Dynamic>::ComputeIntegralXt(RefMatrix& Kbar, RefMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    setScaling(Kbar);

    composeA(Kbar, Bbar);

    firstOrder.ComputeIntegralXt(As, b * scalingSqrd, xInit.cwiseProduct(mulXInit), step / scalingT, out);
}

// ComputeDoubleIntegralXt
template <typename T>
void LDS2OrderUtility<T, Dynamic>::ComputeDoubleIntegralXt(RefMatrix& Kbar, RefMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    setScaling(Kbar);

    composeA(Kbar, Bbar);

    firstOrder.ComputeDoubleIntegralXt(As, b * scalingSqrd, xInit.cwiseProduct(mulXInit), step / scalingT, out);
}
}

#endif
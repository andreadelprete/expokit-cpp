#ifndef INTEGRAL_UTILITY2
#define INTEGRAL_UTILITY2

#include "LDSUtility.hpp"
#include <Eigen/Core>
#include <iostream>

using namespace Eigen;
using namespace std;

#define NHalf N / 2

namespace expokit {

/*
    Wrapper around LDSUtility used to test timescaling optimization
    Eventually it may be extended to a generic second order system intergrator

    TODO:
    - Test for dynamic allocation
    - Rewrite outScaler to avoid matrix multiplication
*/

/*
             _____ _        _   _      
            /  ___| |      | | (_)     
            \ `--.| |_ __ _| |_ _  ___ 
             `--. \ __/ _` | __| |/ __|
            /\__/ / || (_| | |_| | (__ 
            \____/ \__\__,_|\__|_|\___|
*/
template <typename T, int N>
class LDS2OrderUtility {
private:
    LDSUtility<T, N> firstOrder;

    Matrix<T, N, 1> mulXInit, tempOut;
    Matrix<T, N, N> As, outScaler;

    T scalingT;
    T scalingSqrd;

    typedef const Ref<const Matrix<T, N, 1>> RefVector;
    typedef const Ref<const Matrix<T, N, N>> RefMatrix;
    typedef const Ref<const Matrix<T, NHalf, NHalf>> RefSubMatrix;

    typedef Ref<Matrix<T, N, 1>> RefOutVector;

public:
    LDS2OrderUtility();

    /**
     * Compute the value of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeXt(RefSubMatrix& Kbar, RefSubMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out);
    void ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeIntegralXt(RefSubMatrix& Kbar, RefSubMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out);
    void ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the double integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeDoubleIntegralXt(RefSubMatrix& Kbar, RefSubMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out);
    void ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);

private:
    void setScaling(RefSubMatrix& Kbar);
    void composeA(RefSubMatrix& Kbar, RefSubMatrix& Bbar);
};

template <typename T, int N>
LDS2OrderUtility<T, N>::LDS2OrderUtility()
{
    As = Matrix<T, N, N>::Zero();
    As.template block<NHalf, NHalf>(0, NHalf) = Matrix<T, NHalf, NHalf>::Identity();

    outScaler = Matrix<T, N, N>::Zero();
}

// Just to avoid repetition
template <typename T, int N>
inline void LDS2OrderUtility<T, N>::setScaling(RefSubMatrix& Kbar)
{
    // // Defining the scaling dependent on K. More complex stuff is possible.
    // const double l1norm = Kbar.cwiseAbs().colwise().sum().maxCoeff();
    // scalingT = sqrt(1 / l1norm);
    // scalingSqrd = scalingT * scalingT;

    // // Structure to scale only the latter half (speed) of xInit
    // mulXInit = DynVector::Ones(n, 1);
    // mulXInit.block(nHalf, 0, nHalf, 1) = scalingT * mulXInit.block(nHalf, 0, nHalf, 1);

    // // Matrix used to scale back the result - additional descaling is dependent on the operation
    // outScaler.block(0, 0, nHalf, nHalf) = DynMatrix::Identity(nHalf, nHalf);
    // outScaler.block(nHalf, nHalf, nHalf, nHalf) = (1 / scalingT) * DynMatrix::Identity(nHalf, nHalf);
}

template <typename T, int N>
inline void LDS2OrderUtility<T, N>::composeA(RefSubMatrix& Kbar, RefSubMatrix& Bbar)
{
    // As.block(nHalf, nHalf, nHalf, nHalf) = Bbar * scalingT;
    // As.block(nHalf, 0, nHalf, nHalf) = Kbar * scalingSqrd;
}

template <typename T, int N>
void LDS2OrderUtility<T, N>::ComputeXt(RefSubMatrix& Kbar, RefSubMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
}

template <typename T, int N>
void LDS2OrderUtility<T, N>::ComputeIntegralXt(RefSubMatrix& Kbar, RefSubMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
}

template <typename T, int N>
void LDS2OrderUtility<T, N>::ComputeDoubleIntegralXt(RefSubMatrix& Kbar, RefSubMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    // // Building aumented matrix A2
    // A2.template block<N, N>(0, 0) = A;
    // A2.template block<N, 1>(0, N) = b;
    // A2.template block<N, 1>(0, N + 1) = xInit;
    // A2.template block<2, 2>(N, N + 1) = Matrix<T, 2, 2>::Identity();
    // A2 *= step;

    // // Vector z
    // z2(N + 2, 0) = 1;

    // // Matrix exponential and extracting the interesting result
    // expUtil3.computeExpTimesVector(A2, z2, res3);
    // out = res3.template block<N, 1>(0, 0);
}

/*
    INTERFACE UNIFORMITY METHODS
    Be aware that the A matrix must still respect the form [0, eye; Kbar, Bbar].
    Meaning that the upper portion is ignored.
*/
// ComputeXt
template <typename T, int N>
void LDS2OrderUtility<T, N>::ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    ComputeXt(A.template block<NHalf, NHalf>(NHalf, 0) , A.template block<NHalf, NHalf>(NHalf, NHalf), b, xInit, step, out);
}
// ComputeIntegralXt
template <typename T, int N>
void LDS2OrderUtility<T, N>::ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    ComputeIntegralXt(A.template block<NHalf, NHalf>(NHalf, 0) , A.template block<NHalf, NHalf>(NHalf, NHalf), b, xInit, step, out);
}
// ComputeDoubleIntegralXt
template <typename T, int N>
void LDS2OrderUtility<T, N>::ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    ComputeDoubleIntegralXt(A.template block<NHalf, NHalf>(NHalf, 0) , A.template block<NHalf, NHalf>(NHalf, NHalf), b, xInit, step, out);
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
    LDSUtility<T, Dynamic> firstOrder;

    int n, nHalf;
    T scalingT;
    T scalingSqrd;

    typedef Matrix<T, Dynamic, 1> DynVector;
    typedef Matrix<T, Dynamic, Dynamic> DynMatrix;

    typedef const Ref<const DynVector> RefVector;
    typedef const Ref<const DynMatrix> RefMatrix;

    typedef Ref<DynVector> RefOutVector;

    // Preallocating useful stuff
    DynVector mulXInit, tempOut;
    DynMatrix As, outScaler;

public:
    // Would like to forbid creation without specifying a size, but it would prevent use as a class field
    LDS2OrderUtility()
        : LDS2OrderUtility(2)
    {
    } // Creating with dim 2
    LDS2OrderUtility(int n);

    void resize(int n);

    /**
     * Compute the value of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeXt(RefMatrix& Kbar, RefMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out);
    void ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeIntegralXt(RefMatrix& Kbar, RefMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out);
    void ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the double integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeDoubleIntegralXt(RefMatrix& Kbar, RefMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out);
    void ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);

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
    outScaler = DynMatrix::Zero(n, n);
    tempOut.resize(n, 1);
}

// Just to avoid repetition
template <typename T>
inline void LDS2OrderUtility<T, Dynamic>::setScaling(RefMatrix& Kbar)
{
    // Defining the scaling dependent on K. More complex stuff is possible.
    const double l1norm = Kbar.cwiseAbs().colwise().sum().maxCoeff();
    scalingT = sqrt(1 / l1norm);
    scalingSqrd = scalingT * scalingT;

    // Structure to scale only the latter half (speed) of xInit
    mulXInit = DynVector::Ones(n, 1);
    mulXInit.block(nHalf, 0, nHalf, 1) = scalingT * mulXInit.block(nHalf, 0, nHalf, 1);

    // Matrix used to scale back the result - additional descaling is dependent on the operation
    outScaler.block(0, 0, nHalf, nHalf) = DynMatrix::Identity(nHalf, nHalf);
    outScaler.block(nHalf, nHalf, nHalf, nHalf) = (1 / scalingT) * DynMatrix::Identity(nHalf, nHalf);
}

template <typename T>
inline void LDS2OrderUtility<T, Dynamic>::composeA(RefMatrix& Kbar, RefMatrix& Bbar)
{
    As.block(nHalf, nHalf, nHalf, nHalf) = Bbar * scalingT;
    As.block(nHalf, 0, nHalf, nHalf) = Kbar * scalingSqrd;
}

// ComputeXt
template <typename T>
void LDS2OrderUtility<T, Dynamic>::ComputeXt(RefMatrix& Kbar, RefMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    setScaling(Kbar);
    composeA(Kbar, Bbar);

    firstOrder.ComputeXt(As, b * scalingSqrd, xInit.cwiseProduct(mulXInit), step / scalingT, tempOut);

    out = outScaler * tempOut;
}

// ComputeIntegralXt
template <typename T>
void LDS2OrderUtility<T, Dynamic>::ComputeIntegralXt(RefMatrix& Kbar, RefMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    setScaling(Kbar);

    composeA(Kbar, Bbar);

    firstOrder.ComputeIntegralXt(As, b * scalingSqrd, xInit.cwiseProduct(mulXInit), step / scalingT, tempOut);

    outScaler *= scalingT;
    out = outScaler * tempOut;
}

// ComputeDoubleIntegralXt
template <typename T>
void LDS2OrderUtility<T, Dynamic>::ComputeDoubleIntegralXt(RefMatrix& Kbar, RefMatrix& Bbar, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    setScaling(Kbar);

    composeA(Kbar, Bbar);

    firstOrder.ComputeDoubleIntegralXt(As, b * scalingSqrd, xInit.cwiseProduct(mulXInit), step / scalingT, tempOut);

    outScaler *= scalingSqrd;
    out = outScaler * tempOut;
}

/*
    INTERFACE UNIFORMITY METHODS
    Be aware that the A matrix must still respect the form [0, eye; Kbar, Bbar].
    Meaning that the upper portion is ignored.
*/
// ComputeXt
template <typename T>
void LDS2OrderUtility<T, Dynamic>::ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    ComputeXt(A.block(nHalf, 0, nHalf, nHalf), A.block(nHalf, nHalf, nHalf, nHalf), b, xInit, step, out);
}
// ComputeIntegralXt
template <typename T>
void LDS2OrderUtility<T, Dynamic>::ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    ComputeIntegralXt(A.block(nHalf, 0, nHalf, nHalf), A.block(nHalf, nHalf, nHalf, nHalf), b, xInit, step, out);
}
// ComputeDoubleIntegralXt
template <typename T>
void LDS2OrderUtility<T, Dynamic>::ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    ComputeDoubleIntegralXt(A.block(nHalf, 0, nHalf, nHalf), A.block(nHalf, nHalf, nHalf, nHalf), b, xInit, step, out);
}
}

#endif
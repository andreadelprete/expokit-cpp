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

/*
             _____ _        _   _      
            /  ___| |      | | (_)     
            \ `--.| |_ __ _| |_ _  ___ 
             `--. \ __/ _` | __| |/ __|
            /\__/ / || (_| | |_| | (__ 
            \____/ \__\__,_|\__|_|\___|
*/
template <typename T, int N>
class LDSUtility {
private:
    MatrixExponential<T, N + 1> expUtil1;
    MatrixExponential<T, N + 2> expUtil2;
    MatrixExponential<T, N + 3> expUtil3;

    Matrix<T, N + 1, 1> res1;
    Matrix<T, N + 1, N + 1> res1all;

    Matrix<T, N + 2, 1> res2;
    Matrix<T, N + 2, N + 2> res2all;

    Matrix<T, N + 3, 1> res3;
    Matrix<T, N + 3, N + 3> res3all;

    Matrix<T, N + 1, N + 1> A0;
    Matrix<T, N + 1, 1> x0;

    Matrix<T, N + 2, N + 2> A1;

    Matrix<T, N + 2, 1> z1;

    Matrix<T, N + 3, N + 3> A2;

    Matrix<T, N + 3, 1> z2;

    typedef Matrix<T, N, 1> StaVector;
    typedef Matrix<T, N, N> StaMatrix;

    typedef const Ref<const StaVector> RefVector;
    typedef const Ref<const StaMatrix> RefMatrix;

    typedef Ref<StaVector> RefOutVector;

    bool timesVector;
    int TVSquarings, squaringsUsed;

public:
    LDSUtility();

    void useDelta(bool yesOrNo);
    bool isDeltaUsed();
    void setMinSquarings(int minSquarings);
    int getMinSquarings();
    void useTV(bool yesOrNo) { timesVector = yesOrNo; }
    bool isTVUsed() { return timesVector; }
    void setTVSquarings(int TVSquarings) { this->TVSquarings = TVSquarings; }
    int getTVSquarings() { return TVSquarings; }
    int getSquarings() { return squaringsUsed; }

    /**
     * Compute the value of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);
    void BalancedComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);

    /**
     * Compute the value of the double integral of x(T) given x(0)=xInit and the linear dynamics dx = Ax+b
     */
    void ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out);
};

template <typename T, int N>
LDSUtility<T, N>::LDSUtility()
{
    A0 = Matrix<T, N + 1, N + 1>::Zero();
    A1 = Matrix<T, N + 2, N + 2>::Zero();
    z1 = Matrix<T, N + 2, 1>::Zero();
    A2 = Matrix<T, N + 3, N + 3>::Zero();
    z2 = Matrix<T, N + 3, 1>::Zero();

    timesVector = false;
    TVSquarings = -1;
}

// Properties about Delta
template <typename T, int N>
void LDSUtility<T, N>::useDelta(bool yesOrNo)
{
    expUtil1.useDelta(yesOrNo);
    expUtil2.useDelta(yesOrNo);
    expUtil3.useDelta(yesOrNo);
}

template <typename T, int N>
bool LDSUtility<T, N>::isDeltaUsed()
{
    return expUtil1.isDeltaUsed();
}

// Properties about Time Scaling
template <typename T, int N>
void LDSUtility<T, N>::setMinSquarings(int minSquarings)
{
    expUtil1.setMinSquarings(minSquarings);
    expUtil2.setMinSquarings(minSquarings);
    expUtil3.setMinSquarings(minSquarings);
}

template <typename T, int N>
int LDSUtility<T, N>::getMinSquarings()
{
    return expUtil1.getMinSquarings();
}

template <typename T, int N>
void LDSUtility<T, N>::ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    // Building aumented matrix A0
    A0.template block<N, N>(0, 0) = A;
    A0.template block<N, 1>(0, N) = b;
    A0 *= step;

    // Building augmented state x0
    x0 << xInit, 1;

    // Matrix exponential Extracting the interesting result
    if (timesVector) {
        expUtil1.computeExpTimesVector(A0, x0, res1, TVSquarings);
    } else {
        expUtil1.compute(A0, res1all);
        res1.noalias() = res1all * x0;
    }

    out = res1.template block<N, 1>(0, 0);
    squaringsUsed = expUtil1.getSquarings();
}

template <typename T, int N>
void LDSUtility<T, N>::ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
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
    if (timesVector) {
        expUtil2.computeExpTimesVector(A1, z1, res2, TVSquarings);
    } else {
        expUtil2.compute(A1, res2all);
        res2.noalias() = res2all * z1;
    }

    out = res2.template block<N, 1>(0, 0);
    squaringsUsed = expUtil2.getSquarings();
}

template <typename T, int N>
void LDSUtility<T, N>::BalancedComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
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
    if (timesVector) {
        expUtil2.balanceComputeExpTimesVector(A1, z1, res2, TVSquarings);
    } else {
        expUtil2.balanceCompute(A1, res2all);
        // cout << res2all << endl;
        res2.noalias() = res2all * z1;
    }

    out = res2.template block<N, 1>(0, 0);
    squaringsUsed = expUtil2.getSquarings();
}

template <typename T, int N>
void LDSUtility<T, N>::ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
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
    if (timesVector) {
        expUtil3.computeExpTimesVector(A2, z2, res3, TVSquarings);
    } else {
        expUtil3.compute(A2, res3all);
        res3.noalias() = res3all * z2;
    }

    out = res3.template block<N, 1>(0, 0);
    squaringsUsed = expUtil3.getSquarings();
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

    typedef const Ref<const DynVector> RefVector;
    typedef const Ref<const DynMatrix> RefMatrix;

    typedef Ref<DynVector> RefOutVector;

    // Preallocating useful stuff
    DynVector res1, res2, res3, x0, z1, z2;
    DynMatrix res1all, res2all, res3all, A0, A1, A2;

    bool timesVector;
    int TVSquarings, squaringsUsed;

public:
    // Would like to forbid creation without specifying a size, but it would prevent use as a class field
    LDSUtility()
        : LDSUtility(2)
    {
    } // Creating with dim 2
    LDSUtility(int n);

    void useDelta(bool yesOrNo);
    bool isDeltaUsed();
    void setMinSquarings(int minSquarings);
    int getMinSquarings();
    void useTV(bool yesOrNo) { timesVector = yesOrNo; }
    bool isTVUsed() { return timesVector; }
    void setTVSquarings(int TVSquarings) { this->TVSquarings = TVSquarings; }
    int getTVSquarings() { return TVSquarings; }
    int getSquarings() { return squaringsUsed; }

    void resize(int n);

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

// Properties about Delta
template <typename T>
void LDSUtility<T, Dynamic>::useDelta(bool yesOrNo)
{
    expUtil1.useDelta(yesOrNo);
    expUtil2.useDelta(yesOrNo);
    expUtil3.useDelta(yesOrNo);
}

template <typename T>
bool LDSUtility<T, Dynamic>::isDeltaUsed()
{
    return expUtil1.isDeltaUsed();
}

// Properties about Time Scaling
template <typename T>
void LDSUtility<T, Dynamic>::setMinSquarings(int minSquarings)
{
    expUtil1.setMinSquarings(minSquarings);
    expUtil2.setMinSquarings(minSquarings);
    expUtil3.setMinSquarings(minSquarings);
}

template <typename T>
int LDSUtility<T, Dynamic>::getMinSquarings()
{
    return expUtil1.getMinSquarings();
}

template <typename T>
void LDSUtility<T, Dynamic>::resize(int n)
{
    this->n = n;
    timesVector = false;
    TVSquarings = -1;
    squaringsUsed = 0;
    expUtil1.resize(n + 1);
    expUtil2.resize(n + 2);
    expUtil3.resize(n + 3);

    // Stuff for ComputeXt
    res1.resize(n + 1, 1);
    res1all.resize(n + 1, n + 1);
    A0 = DynMatrix::Zero(n + 1, n + 1);
    x0.resize(n + 1, 1);

    // Stuff for ComputeIntegralXt
    res2.resize(n + 2, 1);
    res2all.resize(n + 2, n + 2);
    A1 = DynMatrix::Zero(n + 2, n + 2);
    z1 = DynVector::Zero(n + 2, 1);
    z1(n + 1, 0) = 1; // This is constant

    // Stuff for ComputeDoubleIntegralXt
    res3.resize(n + 3, 1);
    res3all.resize(n + 3, n + 3);
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
    if (timesVector) {
        expUtil1.computeExpTimesVector(A0, x0, res1, TVSquarings);
    } else {
        expUtil1.compute(A0, res1all);
        res1.noalias() = res1all * x0;
    }

    out = res1.block(0, 0, n, 1);
    squaringsUsed = expUtil1.getSquarings();
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
    if (timesVector) {
        expUtil2.computeExpTimesVector(A1, z1, res2, TVSquarings);
    } else {
        expUtil2.compute(A1, res2all);
        res2.noalias() = res2all * z1;
    }

    out = res2.block(0, 0, n, 1);
    squaringsUsed = expUtil2.getSquarings();
}

template <typename T>
void LDSUtility<T, Dynamic>::ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, T step, RefOutVector out)
{
    // Building aumented matrix A2
    A2.block(0, 0, n, n) = A;
    A2.block(0, n, n, 1) = b;
    A2.block(0, n + 1, n, 1) = xInit;
    //A2.block(n, n + 1, 2, 2) = Matrix<T, 2, 2>::Identity();
    A2(n, n + 1) = 1;
    A2(n + 1, n + 2) = 1;
    A2 *= step;

    // Matrix exponential and extracting the interesting result
    if (timesVector) {
        expUtil3.computeExpTimesVector(A2, z2, res3, TVSquarings);
    } else {
        expUtil3.compute(A2, res3all);
        res3.noalias() = res3all * z2;
    }

    out = res3.block(0, 0, n, 1);
    squaringsUsed = expUtil3.getSquarings();
}
}

#endif
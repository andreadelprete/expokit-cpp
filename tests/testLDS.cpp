#include <Eigen/Core>
#include <iostream>

#include "LDS2OrderUtility.hpp"
#include "LDSUtility.hpp"
#include "MatrixExponential.hpp"

using namespace std;
using namespace expokit;

#define INT_STEP 0.01
#define PRECISION_INT 0.1

// Prompt facilitator
bool basicAssertion(const bool test, const char message[]);

// Baseline precision implemented with vanilla Euler integration
Vector4d ComputeXt(Matrix4d& A, Vector4d& b, Vector4d& xInit, double tFinal, double dt);
Vector4d ComputeIntegralXt(Matrix4d& A, Vector4d& b, Vector4d& xInit, double tFinal, double dt);
Vector4d ComputeDoubleIntegralXt(Matrix4d& A, Vector4d& b, Vector4d& xInit, double tFinal, double dt);

// Substitute with an assert when we will be able to make tests
// Comparing with the eigen built in routine
int main()
{
    bool flag = true;

    double stiffness = 10;
    double damping = 1;
    Matrix2d U = Matrix2d::Random();
    Matrix2d Upsilon = U * U.transpose();
    Matrix2d K = Matrix2d::Identity() * stiffness;
    Matrix2d B = Matrix2d::Identity() * damping;
    Matrix4d A = Matrix4d::Zero();
    A.topRightCorner(2, 2) = Matrix2d::Identity();
    A.bottomLeftCorner(2, 2) = -Upsilon * K;
    A.bottomRightCorner(2, 2) = -Upsilon * B;
    
    A << 0,0,1,0,0,0,0,1,10,0,2,0,0,10,0,2;
    
    Vector4d v, ref, res;
    v << 5, 6, 7, 8;

    // Integration utility
    LDSUtility<double, 4> ldsStatic;
    LDSUtility<double, Dynamic> ldsDynamic(4);

    LDS2OrderUtility<double, 4> lds2Static;
    LDS2OrderUtility<double, Dynamic> lds2Dynamic(4);
    Vector4d b, xInit;
    b << 1, 2, 3, 4;
    b *= 0.1;
    xInit << 4, 3, 2, 1;
    xInit *= 0.1;

    // ComputeXt
    ref = ComputeXt(A, b, xInit, 1, INT_STEP);

    ldsStatic.ComputeXt(A, b, xInit, 1, res);
    flag &= basicAssertion(ref.isApprox(res, PRECISION_INT), "Sta ComputeXt");

    ldsDynamic.ComputeXt(A, b, xInit, 1, res);
    flag &= basicAssertion(ref.isApprox(res, PRECISION_INT), "Dyn ComputeXt");

    // lds2Static.ComputeXt(A, b, xInit, 1, res);
    // flag &= basicAssertion(ref.isApprox(res, PRECISION_INT), "Sta ComputeXt TS");

    lds2Dynamic.ComputeXt(A.block(2, 0, 2, 2), A.block(2, 2, 2, 2), b, xInit, 1, res);
    flag &= basicAssertion(ref.isApprox(res, PRECISION_INT), "Dyn ComputeXt TS");

    // ComputeIntegralXt
    ref = ComputeIntegralXt(A, b, xInit, 3, INT_STEP);

    ldsStatic.ComputeIntegralXt(A, b, xInit, 3, res);
    flag &= basicAssertion(ref.isApprox(res, PRECISION_INT), "Sta ComputeIntegralXt");
    ldsDynamic.ComputeIntegralXt(A, b, xInit, 3, res);
    flag &= basicAssertion(ref.isApprox(res, PRECISION_INT), "Dyn ComputeIntegralXt");

    // lds2Static.ComputeIntegralXt(A, b, xInit, 3, res);
    // flag &= basicAssertion(ref.isApprox(res, PRECISION_INT), "Sta ComputeIntegralXt TS");

    lds2Dynamic.ComputeIntegralXt(A.block(2, 0, 2, 2), A.block(2, 2, 2, 2), b, xInit, 3, res);
    flag &= basicAssertion(ref.isApprox(res, PRECISION_INT), "Dyn ComputeIntegralXt TS");

    // ComputeDoubleIntegralXt
    ref = ComputeDoubleIntegralXt(A, b, xInit, 2, INT_STEP);

    ldsStatic.ComputeDoubleIntegralXt(A, b, xInit, 2, res);
    flag &= basicAssertion(ref.isApprox(res, PRECISION_INT), "Sta ComputeDoubleIntegralXt");

    ldsDynamic.ComputeDoubleIntegralXt(A, b, xInit, 2, res);
    flag &= basicAssertion(ref.isApprox(res, PRECISION_INT), "Dyn ComputeDoubleIntegralXt");

    // lds2Static.ComputeDoubleIntegralXt(A, b, xInit, 2, res);
    // flag &= basicAssertion(ref.isApprox(res, PRECISION_INT), "Sta ComputeDoubleIntegralXt TS");

    lds2Dynamic.ComputeDoubleIntegralXt(A.block(2, 0, 2, 2), A.block(2, 2, 2, 2), b, xInit, 2, res);
    flag &= basicAssertion(ref.isApprox(res, PRECISION_INT), "Dyn ComputeDoubleIntegralXt TS");
    // cout << ref - res << endl;

    return !flag;
}

bool basicAssertion(const bool test, const char message[])
{
    if (test)
        cout << message << "\t\t\tOK" << endl;
    else
        cout << message << "\t\t\tFAIL" << endl;

    return test;
}

Vector4d ComputeXt(Matrix4d& A, Vector4d& b, Vector4d& xInit, double tFinal, double dt)
{
    int steps = (int)(tFinal / dt);
    Vector4d result = xInit;
    for (int s = 0; s < steps; s++) {
        result += (A * result + b) * dt;
    }

    return result;
}

Vector4d ComputeIntegralXt(Matrix4d& A, Vector4d& b, Vector4d& xInit, double tFinal, double dt)
{
    int steps = (int)(tFinal / dt);
    // Vector4d result = Vector4d::Zero();
    // for (int s = 0; s < steps; s++) {
    //     result += ComputeXt(A, b, xInit, dt * s, dt) * dt;
    // }

    Vector4d result = Vector4d::Zero();
    Vector4d dx = xInit;
    for (int s = 0; s < steps; s++) {
        dx = ComputeXt(A, b, dx, dt, dt * INT_STEP);
        result += dx * dt;
    }

    return result;
}

Vector4d ComputeDoubleIntegralXt(Matrix4d& A, Vector4d& b, Vector4d& xInit, double tFinal, double dt)
{
    int steps = (int)(tFinal / dt);
    Vector4d result = Vector4d::Zero();
    for (int s = 0; s < steps; s++) {
        result += ComputeIntegralXt(A, b, xInit, dt * s, dt) * dt;
    }

    return result;
}

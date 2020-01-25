#include <Eigen/Core>
#include <iostream>

#include "LDSUtility.hpp"
#include "MatrixExponential.hpp"

using namespace std;
using namespace expokit;

#define INT_STEP 0.01
#define PRECISION_EXP 0.001
#define PRECISION_INT 0.01

// Prompt facilitator
void basicAssertion(const bool test, const char message[]);

// Baseline precision implemented with vanilla Euler integration
Vector3d ComputeXt(Matrix3d& A, Vector3d& b, Vector3d& xInit, double tFinal, double dt);
Vector3d ComputeIntegralXt(Matrix3d& A, Vector3d& b, Vector3d& xInit, double tFinal, double dt);
Vector3d ComputeDoubleIntegralXt(Matrix3d& A, Vector3d& b, Vector3d& xInit, double tFinal, double dt);

// Substitute with an assert when we will be able to make tests
// Comparing with the eigen built in routine
int main()
{

    // Matrix exponential correctness
    Matrix3d A, result;
    A << -8, 1, -6, 3, -5, 7, -4, -9, 2;
    Vector3d v, res1, res2;
    v << 5, 6, 7;
    MatrixExponential<double, 3> util;
    Matrix3d Aref = A.exp();

    util.compute(A, result);
    basicAssertion(Aref.isApprox(result), "Matrix Compiute");

    util.computeExpTimesVector(A, v, res1);
    res2 = Aref * v;
    basicAssertion(res1.isApprox(res2, PRECISION_EXP), "Matrix TimesVector");

    // Integration utility
    LDSUtility<double, 3> lds;
    Vector3d b, xInit;
    b << 1, 2, 3;
    xInit << 3, 2, 1;

    lds.ComputeXt(A, b, xInit, 5, res1);
    res2 = ComputeXt(A, b, xInit, 5, INT_STEP);
    // cout << res1 - res2 << endl;
    basicAssertion(res1.isApprox(res2, PRECISION_INT), "ComputeXt");

    lds.ComputeIntegralXt(A, b, xInit, 3, res1);
    res2 = ComputeIntegralXt(A, b, xInit, 3, INT_STEP);
    // cout << res1 - res2 << endl;
    basicAssertion(res1.isApprox(res2, PRECISION_INT), "ComputeIntegralXt");

    lds.ComputeDoubleIntegralXt(A, b, xInit, 2, res1);
    res2 = ComputeDoubleIntegralXt(A, b, xInit, 2, INT_STEP);
    // cout << res1 - res2 << endl;
    basicAssertion(res1.isApprox(res2, PRECISION_INT), "ComputeDoubleIntegralXt");

    return 0;
}

void basicAssertion(const bool test, const char message[])
{
    if (test)
        cout << message << "\t\t\tOK" << endl;
    else
        cout << message << "\t\t\tFAIL" << endl;
}

Vector3d ComputeXt(Matrix3d& A, Vector3d& b, Vector3d& xInit, double tFinal, double dt)
{
    int steps = (int)(tFinal / dt);
    Vector3d result = xInit;
    for (int s = 0; s < steps; s++) {
        result += (A * result + b) * dt;
    }

    return result;
}

Vector3d ComputeIntegralXt(Matrix3d& A, Vector3d& b, Vector3d& xInit, double tFinal, double dt)
{
    int steps = (int)(tFinal / dt);
    // Vector3d result = Vector3d::Zero();
    // for (int s = 0; s < steps; s++) {
    //     result += ComputeXt(A, b, xInit, dt * s, dt) * dt;
    // }

    Vector3d result = Vector3d::Zero();
    Vector3d dx = xInit;
    for (int s = 0; s < steps; s++) {
        dx = ComputeXt(A, b, dx, dt, dt * INT_STEP);
        result += dx * dt;
    }

    return result;
}

Vector3d ComputeDoubleIntegralXt(Matrix3d& A, Vector3d& b, Vector3d& xInit, double tFinal, double dt)
{
    int steps = (int)(tFinal / dt);
    Vector3d result = Vector3d::Zero();
    for (int s = 0; s < steps; s++) {
        result += ComputeIntegralXt(A, b, xInit, dt * s, dt) * dt;
    }

    return result;
}

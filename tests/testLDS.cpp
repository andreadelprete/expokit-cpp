#include <Eigen/Core>
#include <iostream>

#include "LDS2OrderUtility.hpp"
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
Vector4d ComputeXt(Matrix4d& A, Vector4d& b, Vector4d& xInit, double tFinal, double dt);
Vector4d ComputeIntegralXt(Matrix4d& A, Vector4d& b, Vector4d& xInit, double tFinal, double dt);
Vector4d ComputeDoubleIntegralXt(Matrix4d& A, Vector4d& b, Vector4d& xInit, double tFinal, double dt);

// Substitute with an assert when we will be able to make tests
// Comparing with the eigen built in routine
/*
    TODO 
    - THIS CHECKS ONLY STATIC CODE, EXTEND TO DYNAMIC
*/
int main()
{
    Matrix4d A;
    A << 0.263864065047873, 0.472731157423219, 0.122861379094837, 0.672122191028925,
        0.00754540740531273, 0.438417287905947, 0.0645734168742309, 0.168649768409165,
        0.313604584288928, 0.619345999198629, 0.0768998362137855, 0.539774473871385,
        0.0104497380027225, 0.197138806124144, 0.242873531490832, 0.113272925826047;
    Vector4d v, res1, res2, res3;
    v << 5, 6, 7, 8;

    // Integration utility
    LDSUtility<double, 4> lds;
    LDS2OrderUtility<double, Dynamic> lds2(4);
    Vector4d b, xInit;
    b << 1, 2, 3, 4;
    xInit << 4, 3, 2, 1;

    lds.ComputeXt(A, b, xInit, 5, res1);
    lds2.ComputeXt(A, b, xInit, 5, res3);
    res2 = ComputeXt(A, b, xInit, 5, INT_STEP);
    cout << res1 - res2 << endl;
    basicAssertion(res1.isApprox(res2, PRECISION_INT), "ComputeXt");
    basicAssertion(res1.isApprox(res3, PRECISION_INT), "ComputeXtTimeScaling");

    lds.ComputeIntegralXt(A, b, xInit, 3, res1);
    res2 = ComputeIntegralXt(A, b, xInit, 3, INT_STEP);
    cout << res1 - res2 << endl;
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

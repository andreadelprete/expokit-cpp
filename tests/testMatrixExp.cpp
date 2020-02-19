#include <Eigen/Core>
#include <iostream>

#include "MatrixExponential.hpp"

using namespace std;
using namespace expokit;


#define PRECISION_EXP 0.001

bool basicAssertion(const bool test, const char message[]);

int main()
{
    bool flag = true;

    // Matrix exponential correctness
    Matrix4d A, result;
    A << 0.263864065047873, 0.472731157423219, 0.122861379094837, 0.672122191028925,
        0.00754540740531273, 0.438417287905947, 0.0645734168742309, 0.168649768409165,
        0.313604584288928, 0.619345999198629, 0.0768998362137855, 0.539774473871385,
        0.0104497380027225, 0.197138806124144, 0.242873531490832, 0.113272925826047;
    Vector4d v, res1, res2;
    v << 5, 6, 7, 8;
    Matrix4d Aref = A.exp();

    MatrixExponential<double, 4> utilS; // Static
    MatrixExponential<double, Dynamic> utilD(4); // Dynamic


    utilS.compute(A, result);
    flag &= basicAssertion(Aref.isApprox(result), "Static    Matrix Compute");
    utilD.compute(A, result);
    flag &= basicAssertion(Aref.isApprox(result), "Dynamic   Matrix Compute");

    res2 = Aref * v;
    utilS.computeExpTimesVector(A, v, res1);
    flag &= basicAssertion(res1.isApprox(res2, PRECISION_EXP), "Static  MatrixTimesVector");
    utilD.computeExpTimesVector(A, v, res1);
    flag &= basicAssertion(res1.isApprox(res2, PRECISION_EXP), "Dynamic MatrixTimesVector");

    // cout << Aref - result << endl;
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
#include "LDSUtility.hpp"
#include "utils/readTSV.hpp"

#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <iostream>



using namespace std;
using namespace expokit;
using namespace Eigen;

#define N 4
#define M N * 3 * 2
#define TIMESTEP 0.005

bool basicAssertion(const bool test, const char message[]);

int main()
{
    cout << "Running speed ups benchmark" << endl;

    // Load simulation data
    vector<MatrixXd> vecA = readTSV("exampleData/logStuffA", M, M);
    vector<MatrixXd> vecb = readTSV("exampleData/logStuffb", M, 1);
    vector<MatrixXd> vecxInit = readTSV("exampleData/logStuffxInit", M, 1);
    Matrix<double, M, 1> refInt;
    Matrix<double, M, 1> refDoubleInt;
    Matrix<double, M, 1> testInt;
    Matrix<double, M, 1> testDoubleInt;
    Vector2d e;

    // DeltaOFF_TSOFF_TVOFF
    // LDSUtility<double, M> util;
    LDSUtility<double, Dynamic> util(M);

    bool flag = true;
    for (unsigned int i = 0; i < vecA.size(); ++i) {
        util.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, refInt);
        util.ComputeDoubleIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, refDoubleInt);

        util.ComputeIntegrals(vecA[i], vecb[i], vecxInit[i], TIMESTEP, testInt, testDoubleInt);

        flag &= basicAssertion((refInt - testInt).norm() > 0, "Integral errore above bounds");
        
        flag &= basicAssertion((refDoubleInt - testDoubleInt).norm() > 0, "Double integral errore above bounds");

        if(!flag)
            return -1;
    }

    return 0;
}

bool basicAssertion(const bool test, const char message[])
{
    if (test)
        cout << message << "\t\t\tOK" << endl;
    else
        cout << message << "\t\t\tFAIL" << endl;

    return test;
}
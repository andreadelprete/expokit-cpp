#include "LDSUtility.hpp"
#include "utils/readTSV.hpp"

#define PROFILERYES
#include "utils/stop-watch.h"
// #include "consim/utils/stop-watch.hpp"

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
    cout << "Test computation of both integrals (first and double) with a single method" << endl;

    // Load simulation data
    // string path = "/home/student/devel/src/expokit/cpp/";
    string path = "";
    vector<MatrixXd> vecA = readTSV(path+"exampleData/logStuffA", M, M);
    vector<MatrixXd> vecb = readTSV(path+"exampleData/logStuffb", M, 1);
    vector<MatrixXd> vecxInit = readTSV(path+"exampleData/logStuffxInit", M, 1);
    Matrix<double, M, 1> refInt;
    Matrix<double, M, 1> refDoubleInt;
    Matrix<double, M, 1> testInt;
    Matrix<double, M, 1> testDoubleInt;
    Vector2d e;

    // DeltaOFF_TSOFF_TVOFF
    // LDSUtility<double, M> util;
    LDSUtility<double, Dynamic> util(M);
    
    cout<<"Matrix size "<<M<<endl;

    bool flag = true;
    for (unsigned int i = 0; i < vecA.size(); ++i) {
        START_PROFILER("computeIntegralXt + computeDoubleIntegralXt");
        util.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, refInt);
        util.ComputeDoubleIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, refDoubleInt);
        STOP_PROFILER("computeIntegralXt + computeDoubleIntegralXt");

        util.useTV(false);
        START_PROFILER("computeBothIntegrals");
        util.ComputeIntegrals(vecA[i], vecb[i], vecxInit[i], TIMESTEP, testInt, testDoubleInt);
        STOP_PROFILER("computeBothIntegrals");

        util.useTV(true);
        START_PROFILER("computeBothIntegrals with TV");
        util.ComputeIntegrals(vecA[i], vecb[i], vecxInit[i], TIMESTEP, testInt, testDoubleInt);
        STOP_PROFILER("computeBothIntegrals with TV");

        util.useBalancing(true);  
        START_PROFILER("computeBothIntegrals with TV+Balancing");
        util.ComputeIntegrals(vecA[i], vecb[i], vecxInit[i], TIMESTEP, testInt, testDoubleInt);
        STOP_PROFILER("computeBothIntegrals with TV+Balancing");
        util.useBalancing(false);

        flag &= basicAssertion((refInt - testInt).norm() > 0, "Integral errore above bounds");
        
        flag &= basicAssertion((refDoubleInt - testDoubleInt).norm() > 0, "Double integral errore above bounds");

        if(!flag)
            return -1;
    }

    getProfiler().report_all(3);

    return 0;
}

bool basicAssertion(const bool test, const char message[])
{
    if (!test)
        cout << message << "\t\t\tFAIL" << endl;
    // else
        // cout << message << "\t\t\tOK" << endl;

    return test;
}
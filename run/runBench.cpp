#include <Eigen/Core>
#include <cstdlib>
#include <iostream>

#include "LDS2OrderUtility.hpp"

#include "utils/readTSV.hpp"
#include "utils/stop-watch.h"

using namespace std;
using namespace expokit;

#define N 4
#define M N * 3 * 2

/*
    argv[0] - executable name
    argv[1] - minSquarings
*/
int main(int argc, char* argv[])
{
    cout << "Running speed ups benchmark" << endl;

    // Load simulation data
    vector<MatrixXd> vecA = readTSV("exampleData/logStuffA", M, M);
    vector<MatrixXd> vecb = readTSV("exampleData/logStuffA", M, M);
    vector<MatrixXd> vecxInit = readTSV("exampleData/logStuffA", M, M);

    Matrix<double, M, 1> out;

    // LDS2OrderUtility<double, Dynamic> util2(M);
    LDS2OrderUtility<double, M> deltaOn2;
    deltaOn2.useDelta(true);
    LDS2OrderUtility<double, M> deltaOff2;

    LDSUtility<double, M> deltaOn;
    deltaOn.useDelta(true);
    LDSUtility<double, M> deltaOff;

    // if (argc > 1) {
    //     deltaOn.setMinSquarings(atoi(argv[1]));
    // }
    for (int i = 0; i < 100; ++i) {
        for (unsigned int i = 0; i < vecA.size(); ++i) {
            START_PROFILER("DeltaOFF_TSOFF");
            deltaOff.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], 0.005, out);
            STOP_PROFILER("DeltaOFF_TSOFF");

            START_PROFILER("DeltaON_TSOFF");
            deltaOn.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], 0.005, out);
            STOP_PROFILER("DeltaON_TSOFF");

            START_PROFILER("DeltaOFF_TSON");
            deltaOff2.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], 0.005, out);
            STOP_PROFILER("DeltaOFF_TSON");

            START_PROFILER("DeltaON_TSON");
            deltaOff2.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], 0.005, out);
            STOP_PROFILER("DeltaON_TSON");
        }
    }

    getProfiler().report_all(3);
}
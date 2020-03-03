#include <Eigen/Core>
#include <cstdlib>
#include <iostream>

#include "LDS2OrderUtility.hpp"

#include "utils/readTSV.hpp"
#include "utils/stop-watch.h"
#include "utils/statistics.hpp"

using namespace std;
using namespace expokit;

#define N 4
#define M N * 3 * 2

/*
    argv[0] - executable name
    argv[1] - minSquarings
*/
int main()
{
    cout << "Running speed ups benchmark" << endl;

    Statistics stats;

    // Load simulation data
    vector<MatrixXd> vecA = readTSV("exampleData/logStuffA", M, M);
    vector<MatrixXd> vecb = readTSV("exampleData/logStuffb", M, 1);
    vector<MatrixXd> vecxInit = readTSV("exampleData/logStuffxInit", M, 1);
    Matrix<double, M, 1> out;

    // DeltaOFF_TSOFF_TVOFF
    LDSUtility<double, M> d0t0v0;

    // DeltaON_TSOFF_TVOFF
    LDSUtility<double, M> d1t0v0;
    d1t0v0.useDelta(true);

    // DeltaX_TSOFF_TVON
    LDSUtility<double, M> dxt0v1;
    dxt0v1.useTV(true);

    // DeltaOFF_TSON_TVOFF
    LDS2OrderUtility<double, M> d0t1v0;
    
    // DeltaON_TSON_TVOFF
    LDS2OrderUtility<double, M> d1t1v0;
    d1t1v0.useDelta(true);

    // DeltaX_TSON_TVON
    LDS2OrderUtility<double, M> dxt1v1;
    d1t1v0.useTV(true);



    for (unsigned int i = 0; i < vecA.size(); ++i) {
        stats.store("cojone", (double) i);
        START_PROFILER("DeltaOFF_TSOFF_TVOFF");
        d0t0v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], 0.005, out);
        STOP_PROFILER("DeltaOFF_TSOFF_TVOFF");

        START_PROFILER("DeltaON__TSOFF_TVOFF");
        d1t0v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], 0.005, out);
        STOP_PROFILER("DeltaON__TSOFF_TVOFF");

        START_PROFILER("DeltaX___TSOFF_TVON");
        dxt0v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], 0.005, out);
        STOP_PROFILER("DeltaX___TSOFF_TVON");

        START_PROFILER("DeltaOFF_TSON__TVOFF");
        d0t1v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], 0.005, out);
        STOP_PROFILER("DeltaOFF_TSON__TVOFF");

        START_PROFILER("DeltaON__TSON__TVOFF");
        d1t1v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], 0.005, out);
        STOP_PROFILER("DeltaON__TSON__TVOFF");

        START_PROFILER("DeltaX___TSON__TVON");
        dxt1v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], 0.005, out);
        STOP_PROFILER("DeltaX___TSON__TVON");                
    }


    // Finding optimal value of vec_squarings for Times Vector (TV)
    // for (int k = -1; k < 15; k++)
    // {
    //     string nopad = to_string(k);
    //     std::string padded = std::string(2 - nopad.length(), '0') + nopad;
    //     string namet0 = "T0VECHO" + padded;
    //     string namet1 = "T1VECHO" + padded;

    //     dxt0v1.setTVSquarings(k);
    //     for (unsigned int i = 0; i < vecA.size(); ++i) {
    //         START_PROFILER(namet0);
    //         dxt0v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], 0.005, out);
    //         STOP_PROFILER(namet0);
    //     }

    //     dxt1v1.setTVSquarings(k);
    //     for (unsigned int i = 0; i < vecA.size(); ++i) {
    //         START_PROFILER(namet1);
    //         dxt1v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], 0.005, out);
    //         STOP_PROFILER(namet1);
    //     }
    // }
    // The magic number seems to be 5 for dxt0v1
    // Completely irrelevant for dxt1v1
    


    getProfiler().report_all(3);
    // stats.report_all();
}
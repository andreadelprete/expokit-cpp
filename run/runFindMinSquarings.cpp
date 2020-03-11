#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "LDS2OrderUtility.hpp"

#include "utils/readTSV.hpp"
#include "utils/statistics.hpp"
#include "utils/stop-watch.h"
#include "utils/errorComputing.hpp"


using namespace std;
using namespace expokit;
using namespace Eigen;

#define N 4
#define M N * 3 * 2
#define TIMESTEP 0.005


int main()
{
    cout << "Running minimum squarings investigation" << endl;

    Statistics stats;

    // Load simulation data
    vector<MatrixXd> vecA = readTSV("exampleData/logStuffA", M, M);
    vector<MatrixXd> vecb = readTSV("exampleData/logStuffb", M, 1);
    vector<MatrixXd> vecxInit = readTSV("exampleData/logStuffxInit", M, 1);
    Matrix<double, M, 1> out;
    Matrix<double, M, 1> ref;
    Vector2d e;

    // DeltaOFF_TSOFF_TVOFF
    LDSUtility<double, M> d0t0v0;

    // DeltaON_TSOFF_TVOFF
    LDSUtility<double, M> d1t0v0;
    d1t0v0.useDelta(true);

    // DeltaON_TSON_TVOFF
    LDS2OrderUtility<double, M> d1t1v0;
    d1t1v0.useDelta(true);

    // Measuring the effect that minSquarings has on both the error and the performance
    for (int k = 5; k < 16; k++) {

        d1t0v0.setMinSquarings(k);
        d1t1v0.setMinSquarings(k);

        string nopad = to_string(k);
        std::string padded = std::string(2 - nopad.length(), '0') + nopad;
        string name = "_SQUARINGS_" + padded;

        for (unsigned int i = 0; i < vecA.size(); ++i) {
            // Using not optimazed method for error baseline
            d0t0v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, ref);

            // Delta without Time Scaling
            START_PROFILER("DeltaON__TSOFF_TVOFF" + name);
            d1t0v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
            STOP_PROFILER("DeltaON__TSOFF_TVOFF" + name);

            stats.store("DeltaON__TSOFF_TVOFF" + name, d1t0v0.getSquarings());
            e = computeErrorIntegral24(out, ref, TIMESTEP * i);
            stats.store("DeltaON__TSOFF_TVOFF" + name + "_ERROR2", e(0));
            stats.store("DeltaON__TSOFF_TVOFF" + name + "_ERRORINF", e(1));

            // Delta with Time Scaling
            START_PROFILER("DeltaON__TSON__TVOFF" + name);
            d1t1v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
            STOP_PROFILER("DeltaON__TSON__TVOFF" + name);

            stats.store("DeltaON__TSON__TVOFF" + name, d1t1v0.getSquarings());
            e = computeErrorIntegral24(out, ref, TIMESTEP * i);
            stats.store("DeltaON__TSON__TVOFF" + name + "_ERROR2", e(0));
            stats.store("DeltaON__TSON__TVOFF" + name + "_ERRORINF", e(1));
        }
    }

    getProfiler().report_all(3);
    stats.report_all();
}
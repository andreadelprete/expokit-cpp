#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "LDS2OrderUtility.hpp"

#include "utils/errorComputing.hpp"
#include "utils/readTSV.hpp"
#include "utils/statistics.hpp"
#include "utils/stop-watch.h"

using namespace std;
using namespace expokit;
using namespace Eigen;

#define N 4
#define M N * 3 * 2
#define TIMESTEP 0.005

int main()
{
    cout << "Running speed ups benchmark" << endl;

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

    // DeltaX_TSOFF_TVON
    LDSUtility<double, M> dxt0v1;
    dxt0v1.useTV(true);
    dxt0v1.setTVSquarings(5);

    for (unsigned int i = 0; i < vecA.size(); ++i) {
        START_PROFILER("ComputeIntegralXt");
        d0t0v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, ref);
        STOP_PROFILER("ComputeIntegralXt");
        // cout << ref << endl;
        stats.store("ComputeIntegralXt", d0t0v0.getSquarings());

        START_PROFILER("BalancedComputeIntegralXt");
        d0t0v0.BalancedComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
        STOP_PROFILER("BalancedComputeIntegralXt");
        // cout << out << endl;
        stats.store("BalancedComputeIntegralXt", d0t0v0.getSquarings());
        e = computeErrorIntegral24(out, ref, TIMESTEP * i);
        stats.store("ERROR2___BalancedComputeIntegralXt", e(0));
        stats.store("ERRORINF_BalancedComputeIntegralXt", e(1));

        START_PROFILER("BalancedComputeIntegralXtTV");
        dxt0v1.BalancedComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
        STOP_PROFILER("BalancedComputeIntegralXtTV");
        // cout << out << endl;
        stats.store("BalancedComputeIntegralXtTV", dxt0v1.getSquarings());
        e = computeErrorIntegral24(out, ref, TIMESTEP * i);
        stats.store("ERROR2___BalancedComputeIntegralXtTV", e(0));
        stats.store("ERRORINF_BalancedComputeIntegralXtTV", e(1));
    }

    getProfiler().report_all(3);
    stats.report_all();
}
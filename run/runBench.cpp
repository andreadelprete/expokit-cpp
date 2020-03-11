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

    // DeltaON_TSOFF_TVOFF
    LDSUtility<double, M> d1t0v0;
    d1t0v0.useDelta(true);

    // DeltaX_TSOFF_TVON
    LDSUtility<double, M> dxt0v1;
    dxt0v1.useTV(true);
    dxt0v1.setTVSquarings(5);

    // DeltaOFF_TSON_TVOFF
    LDS2OrderUtility<double, M> d0t1v0;

    // DeltaON_TSON_TVOFF
    LDS2OrderUtility<double, M> d1t1v0;
    d1t1v0.useDelta(true);

    // DeltaX_TSON_TVON
    LDS2OrderUtility<double, M> dxt1v1;
    dxt1v1.useTV(true);
    dxt1v1.setTVSquarings(5);

    for (unsigned int i = 0; i < vecA.size(); ++i) {
        START_PROFILER("DeltaOFF_TSOFF_TVOFF");
        d0t0v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, ref);
        STOP_PROFILER("DeltaOFF_TSOFF_TVOFF");
        stats.store("SQUARINGS_DeltaOFF_TSOFF_TVOFF", d0t0v0.getSquarings());

        START_PROFILER("DeltaON__TSOFF_TVOFF");
        d1t0v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
        STOP_PROFILER("DeltaON__TSOFF_TVOFF");
        stats.store("SQUARINGS_DeltaON__TSOFF_TVOFF", d1t0v0.getSquarings());
        e = computeErrorIntegral24(out, ref, TIMESTEP * i);
        stats.store("ERROR2___DeltaON__TSOFF_TVOFF", e(0));
        stats.store("ERRORINF_DeltaON__TSOFF_TVOFF", e(1));

        START_PROFILER("DeltaX___TSOFF_TVON");
        dxt0v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
        STOP_PROFILER("DeltaX___TSOFF_TVON");
        stats.store("SQUARINGS_DeltaX___TSOFF_TVON", dxt0v1.getSquarings());
        e = computeErrorIntegral24(out, ref, TIMESTEP * i);
        stats.store("ERROR2___DeltaX___TSOFF_TVON", e(0));
        stats.store("ERRORINF_DeltaX___TSOFF_TVON", e(1));

        START_PROFILER("DeltaOFF_TSON__TVOFF");
        d0t1v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
        STOP_PROFILER("DeltaOFF_TSON__TVOFF");
        stats.store("SQUARINGS_DeltaOFF_TSON__TVOFF", d0t1v0.getSquarings());
        e = computeErrorIntegral24(out, ref, TIMESTEP * i);
        stats.store("ERROR2___DeltaOFF_TSON__TVOFF", e(0));
        stats.store("ERRORINF_DeltaOFF_TSON__TVOFF", e(1));

        START_PROFILER("DeltaON__TSON__TVOFF");
        d1t1v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
        STOP_PROFILER("DeltaON__TSON__TVOFF");
        stats.store("SQUARINGS_DeltaON__TSON__TVOFF", d1t1v0.getSquarings());
        e = computeErrorIntegral24(out, ref, TIMESTEP * i);
        stats.store("ERROR2___DeltaON__TSON__TVOFF", e(0));
        stats.store("ERRORINF_DeltaON__TSON__TVOFF", e(1));

        START_PROFILER("DeltaX___TSON__TVON");
        dxt1v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
        STOP_PROFILER("DeltaX___TSON__TVON");
        stats.store("SQUARINGS_DeltaX___TSON__TVON", dxt1v1.getSquarings());
        e = computeErrorIntegral24(out, ref, TIMESTEP * i);
        stats.store("ERROR2___DeltaX___TSON__TVON", e(0));
        stats.store("ERRORINF_DeltaX___TSON__TVON", e(1));
    }

    // Print out results
    getProfiler().report_all(3);
    stats.report_all();
}


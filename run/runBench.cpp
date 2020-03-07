#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "LDS2OrderUtility.hpp"

#include "utils/readTSV.hpp"
#include "utils/statistics.hpp"
#include "utils/stop-watch.h"

using namespace std;
using namespace expokit;
using namespace Eigen;

#define N 4
#define M N * 3 * 2
#define TIMESTEP 0.005

Vector2d computeErrorIntegral(Matrix<double, M, 1>& x, Matrix<double, M, 1>& xRef, double dt);

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
        e = computeErrorIntegral(out, ref, TIMESTEP * i);
        stats.store("ERROR2___DeltaON__TSOFF_TVOFF", e(0));
        stats.store("ERRORINF_DeltaON__TSOFF_TVOFF", e(1));

        START_PROFILER("DeltaX___TSOFF_TVON");
        dxt0v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
        STOP_PROFILER("DeltaX___TSOFF_TVON");
        stats.store("SQUARINGS_DeltaX___TSOFF_TVON", dxt0v1.getSquarings());
        e = computeErrorIntegral(out, ref, TIMESTEP * i);
        stats.store("ERROR2___DeltaX___TSOFF_TVON", e(0));
        stats.store("ERRORINF_DeltaX___TSOFF_TVON", e(1));

        START_PROFILER("DeltaOFF_TSON__TVOFF");
        d0t1v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
        STOP_PROFILER("DeltaOFF_TSON__TVOFF");
        stats.store("SQUARINGS_DeltaOFF_TSON__TVOFF", d0t1v0.getSquarings());
        e = computeErrorIntegral(out, ref, TIMESTEP * i);
        stats.store("ERROR2___DeltaOFF_TSON__TVOFF", e(0));
        stats.store("ERRORINF_DeltaOFF_TSON__TVOFF", e(1));

        START_PROFILER("DeltaON__TSON__TVOFF");
        d1t1v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
        STOP_PROFILER("DeltaON__TSON__TVOFF");
        stats.store("SQUARINGS_DeltaON__TSON__TVOFF", d1t1v0.getSquarings());
        e = computeErrorIntegral(out, ref, TIMESTEP * i);
        stats.store("ERROR2___DeltaON__TSON__TVOFF", e(0));
        stats.store("ERRORINF_DeltaON__TSON__TVOFF", e(1));

        START_PROFILER("DeltaX___TSON__TVON");
        dxt1v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
        STOP_PROFILER("DeltaX___TSON__TVON");
        stats.store("SQUARINGS_DeltaX___TSON__TVON", dxt1v1.getSquarings());
        e = computeErrorIntegral(out, ref, TIMESTEP * i);
        stats.store("ERROR2___DeltaX___TSON__TVON", e(0));
        stats.store("ERRORINF_DeltaX___TSON__TVON", e(1));
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
    //         dxt0v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
    //         STOP_PROFILER(namet0);
    //     }

    //     dxt1v1.setTVSquarings(k);
    //     for (unsigned int i = 0; i < vecA.size(); ++i) {
    //         START_PROFILER(namet1);
    //         dxt1v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
    //         STOP_PROFILER(namet1);
    //     }
    // }
    // The magic number seems to be 5 for dxt0v1
    // Completely irrelevant for dxt1v1


    // Measuring the effect that minSquarings has on both the error and the performance
    // for (int k = 5; k < 16; k++) {

    //     d1t0v0.setMinSquarings(k);
    //     d1t1v0.setMinSquarings(k);

    //     string nopad = to_string(k);
    //     std::string padded = std::string(2 - nopad.length(), '0') + nopad;
    //     string name = "_SQUARINGS_" + padded;

    //     for (unsigned int i = 0; i < vecA.size(); ++i) {
    //         // Using not optimazed method for error baseline
    //         d0t0v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, ref);

    //         // Delta without Time Scaling
    //         START_PROFILER("DeltaON__TSOFF_TVOFF" + name);
    //         d1t0v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
    //         STOP_PROFILER("DeltaON__TSOFF_TVOFF" + name);

    //         stats.store("DeltaON__TSOFF_TVOFF" + name, d1t0v0.getSquarings());
    //         e = computeErrorIntegral(out, ref, TIMESTEP * i);
    //         stats.store("DeltaON__TSOFF_TVOFF" + name + "_ERROR2", e(0));
    //         stats.store("DeltaON__TSOFF_TVOFF" + name + "_ERRORINF", e(1));

    //         // Delta with Time Scaling
    //         START_PROFILER("DeltaON__TSON__TVOFF" + name);
    //         d1t1v0.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
    //         STOP_PROFILER("DeltaON__TSON__TVOFF" + name);

    //         stats.store("DeltaON__TSON__TVOFF" + name, d1t1v0.getSquarings());
    //         e = computeErrorIntegral(out, ref, TIMESTEP * i);
    //         stats.store("DeltaON__TSON__TVOFF" + name + "_ERROR2", e(0));
    //         stats.store("DeltaON__TSON__TVOFF" + name + "_ERRORINF", e(1));
    //     }
    // }

    getProfiler().report_all(3);
    stats.report_all();
}

// This function gets the result of the integration
// and computes the error in terms of Newton
Vector2d computeErrorIntegral(Matrix<double, M, 1>& x, Matrix<double, M, 1>& xRef, double dt)
{
    const double K = 1e5;
    const double B = 1e2;
    VectorXd p0(12);
    p0 << -0.190000000000000002, 0.150050000000000017, 0.000053853008907090, -0.190000000000000002, -0.150050000000000017, 0.000053853008907090, 0.190000000000000002, 0.150050000000000017, 0.000053853008907090, 0.190000000000000002, -0.150050000000000017, 0.000053853008907090;

    VectorXd fInt = K * (p0 * dt - x.block(0, 0, 12, 1) - B * x.block(12, 0, 12, 1));
    VectorXd fIntRef = K * (p0 * dt - xRef.block(0, 0, 12, 1) - B * xRef.block(12, 0, 12, 1));

    VectorXd diff = fInt - fIntRef;

    Vector2d result;
    result(0) = diff.norm();
    result(1) = diff.lpNorm<Infinity>();

    return result;
}
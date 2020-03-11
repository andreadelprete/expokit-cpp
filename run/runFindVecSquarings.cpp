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

    // DeltaX_TSOFF_TVON
    LDSUtility<double, M> dxt0v1;
    dxt0v1.useTV(true);

    // DeltaX_TSON_TVON
    LDS2OrderUtility<double, M> dxt1v1;
    dxt1v1.useTV(true);

    // Finding optimal value of vec_squarings for Times Vector (TV)
    for (int k = -1; k < 15; k++)
    {
        string nopad = to_string(k);
        std::string padded = std::string(2 - nopad.length(), '0') + nopad;
        string namet0 = "T0VECHO" + padded;
        string namet1 = "T1VECHO" + padded;

        dxt0v1.setTVSquarings(k);
        for (unsigned int i = 0; i < vecA.size(); ++i) {
            START_PROFILER(namet0);
            dxt0v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
            STOP_PROFILER(namet0);
        }

        dxt1v1.setTVSquarings(k);
        for (unsigned int i = 0; i < vecA.size(); ++i) {
            START_PROFILER(namet1);
            dxt1v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
            STOP_PROFILER(namet1);
        }
    }
    // The magic number seems to be 5 for dxt0v1
    // Completely irrelevant for dxt1v1
}
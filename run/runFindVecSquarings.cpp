#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "LDSUtility.hpp"

#include "utils/errorComputing.hpp"
#include "utils/readTSV.hpp"
#include "utils/statistics.hpp"
#include "utils/stop-watch.h"

using namespace std;
using namespace expokit;
using namespace Eigen;

#define TESTS 1000
#define N 8
/*
#define M N * 3 * 2
#define TIMESTEP 0.005
*/

int main(int argc, const char* argv[])
{

    cout << "Running finder of best vec_squarings" << endl;
    /*
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

    // Finding optimal value of vec_squarings for Times Vector (TV)
    // The old method was mapped at k == -2, the new one at k == -1
    for (int k = -2; k < 15; k++)
    {
        string nopad = to_string(k);
        std::string padded = std::string(2 - nopad.length(), '0') + nopad;
        string namet0 = "VS=" + padded;

        dxt0v1.setTVSquarings(k);
        for (unsigned int i = 0; i < vecA.size(); ++i) {
            START_PROFILER(namet0);
            dxt0v1.ComputeIntegralXt(vecA[i], vecb[i], vecxInit[i], TIMESTEP, out);
            STOP_PROFILER(namet0);
        }
    }
    // The magic number seems to be 5 for dxt0v1

    // Print out results
    getProfiler().report_all(3);
    stats.report_all();
*/

    int scale = 100;
    if (argc > 1) {
        scale = std::atoi(argv[1]);
    }

    for (int n = 1; n <= N; ++n) {
        int Ns = n * 10;

        MatrixXd A(Ns, Ns);
        VectorXd v(Ns);
        MatrixExponential<double, Dynamic> expUtil(Ns);
        VectorXd res(Ns);

        for (int i = -1; i < 20; ++i) {
            string nopad = to_string(i);
            std::string padded = std::string(2 - nopad.length(), '0') + nopad;
            string namet0 = "Ns=" + to_string(Ns) + " VS=" + padded;

            for (int k = 0; k < TESTS; ++k) {
                A = MatrixXd::Random(Ns, Ns) * scale;
                v = VectorXd::Random(Ns) * scale;

                START_PROFILER(namet0);
                expUtil.computeExpTimesVector(A, v, res, i);
                STOP_PROFILER(namet0);
            }
        }
    }

    getProfiler().report_all(3);
}
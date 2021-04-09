#include <Eigen/Core>
#include <cstdlib>
#include <iostream>

#include "LDSUtility.hpp"

#include "utils/readTSV.hpp"
#include "utils/stop-watch.h"

using namespace std;
using namespace expokit;

/*
    Testing delta optimization 
    - General correctness - test against official implementation
    - Speed against standard and official methods
*/

#define MANY 100

#define N 4
#define M N * 3 * 2

/*
    argv[0] - executable name
    argv[1] - minSquarings
*/
int main(int argc, char* argv[])
{
    cout << "Start test delta update" << endl;

    // Used with test on simulation data
    vector<MatrixXd> vecA = readTSV("exampleData/logStuffA", M, M);

    Matrix<double, M, M> res0, res1, res2;

    int m2 = int(M / 2);
    double stiffness = 1e5;
    double damping = 1e2;
    MatrixXd U = MatrixXd::Random(m2, m2);
    MatrixXd Upsilon = U * U.transpose();
    MatrixXd K = MatrixXd::Identity(m2, m2) * stiffness;
    MatrixXd B = MatrixXd::Identity(m2, m2) * damping;
    MatrixXd A = MatrixXd::Zero(M, M);
    A.topRightCorner(m2, m2) = MatrixXd::Identity(m2, m2);
    A.bottomLeftCorner(m2, m2) = -Upsilon * K;
    A.bottomRightCorner(m2, m2) = -Upsilon * B;

    // array<MatrixXd, MANY> matrices{};

    MatrixExponential<double, M> deltaOn;
    deltaOn.useDelta(true);
    if (argc > 1)
        deltaOn.setMinSquarings(atoi(argv[1]));
    cout << "Min number of squarings: " << deltaOn.getMinSquarings() << endl;

    MatrixExponential<double, M> deltaOff;

    // Checking correctness
    srand((unsigned int)time(NULL));

    bool whichData = true;
    for (int i = 0; i < MANY; ++i) {

        if (whichData) {
            A = vecA[i];
        } else {
            double smallChange = (rand() % 100) / 100.0;
            int index = rand() % 9;
            int sign = rand();
            if (sign % 2)
                A(index) += smallChange;
            else
                A(index) -= smallChange;
        }

        // Computing using all three methods
        START_PROFILER("testDelta::official");
        res0 = A.exp();
        STOP_PROFILER("testDelta::official");

        START_PROFILER("testDelta::deltaOff");
        deltaOff.compute(A, res1);
        STOP_PROFILER("testDelta::deltaOff");

        START_PROFILER("testDelta::deltaOn");
        deltaOn.compute(A, res2);
        STOP_PROFILER("testDelta::deltaOn");

        if ((res0 - res1).eval().cwiseAbs().sum() > 0) {
            cout << res0 - res1 << endl
                 << i << endl
                 << "ERROR IN STANDARD METHOD" << endl;
            return -1;
        }

        if ((res0 - res2).eval().cwiseAbs().sum() > 5) { // Arbitrary error threshold
            cout << res0 - res1 << endl
                 << i << endl
                 << "ERROR IN DELTA METHOD" << endl;
            return -2;
        }

        // Comment out to see detailed error output
        // cout << round((res0 - res2).eval().cwiseAbs().sum() * 1000.0) / 1000.0 << '\t'
        //      << round((res0 - res2).eval().maxCoeff() * 1000.0) / 1000.0 << '\t'
        //      << deltaOn.wasDeltaUsed() << '\t'
        //      << i << endl;
    }

    getProfiler().report_all(3);
}
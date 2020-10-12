#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "BalancingMethods.hpp"
#include "MatrixExponential.hpp"

#include "utils/errorComputing.hpp"
#include "utils/readTSV.hpp"
#include "utils/statistics.hpp"
#include "utils/stop-watch.h"

using namespace std;
using namespace expokit;
using namespace Eigen;

#ifdef EIGEN_RUNTIME_NO_MALLOC
#define EIGEN_MALLOC_ALLOWED Eigen::internal::set_is_malloc_allowed(true);
#define EIGEN_MALLOC_NOT_ALLOWED Eigen::internal::set_is_malloc_allowed(false);
#else
#define EIGEN_MALLOC_ALLOWED
#define EIGEN_MALLOC_NOT_ALLOWED
#endif

#define N 64
#define TESTS 10000

MatrixXd james2014Generator(int n);

int main(int argc, const char* argv[])
{
    int dim = 10;
    if (argc > 1) {
        dim = std::atoi(argv[1]);
    }

    cout << "Running new balance" << endl;

    MatrixXd A(N, N);
    MatrixXd B(N, N);
    MatrixXd out(N, N);
    VectorXd D(N, N);
    VectorXd Dinv(N, N);

    BalancingMethods<double, Dynamic> util(dim);
    MatrixExponential<double, Dynamic> expUtil(dim);

    for (int i = 0; i < TESTS; ++i) {
        A = james2014Generator(dim);

        // START_PROFILER("EXPnotBalanced");
        // expUtil.compute(A, out);
        // STOP_PROFILER("EXPnotBalanced");

        // START_PROFILER("BALNB");
        util.balanceNew(A, B, D, Dinv);
        double NBnorm = B.cwiseAbs().colwise().sum().maxCoeff();
        cout << NBnorm << "\t";
        // STOP_PROFILER("BALNB");

        // START_PROFILER("EXPNBBalanced");
        // expUtil.compute(B, out);
        // STOP_PROFILER("EXPNBBalanced");

        // START_PROFILER("BALRodney");
        util.balanceRodney(A, B, D, Dinv);
        double RodNorm = B.cwiseAbs().colwise().sum().maxCoeff();
        cout << RodNorm << "\t";
        // STOP_PROFILER("BALRodney");

        // START_PROFILER("EXPRodneyBalanced");
        // expUtil.compute(B, out);
        // STOP_PROFILER("EXPRodneyBalanced");

        // START_PROFILER("BALCombined");
        util.balanceCombined(A, B, D, Dinv);
        double CombNorm = B.cwiseAbs().colwise().sum().maxCoeff();
        cout << CombNorm << "\t";
        // STOP_PROFILER("BALCombined");

        // START_PROFILER("EXPCombinedBalanced");
        // expUtil.compute(B, out);
        // STOP_PROFILER("EXPCombinedBalanced");

        if (RodNorm < NBnorm || RodNorm < CombNorm)
            cout << "Rodney BETTER" << endl;
        else if (NBnorm < RodNorm || NBnorm < RodNorm)
            cout << "NB BETTER" << endl;
        else if (CombNorm < NBnorm || CombNorm < RodNorm)
            cout << "Comb BETTER" << endl;
        else if (CombNorm == NBnorm && CombNorm == RodNorm && RodNorm == NBnorm)
            cout << "ALL EQUAL" << endl;
    }

    // Print out results
    getProfiler().report_all(3);
}

// Generate test matrices
MatrixXd james2014Generator(int n)
{
    Matrix<double, N, N> A = Matrix<double, N, N>::Random();
    Matrix<double, N, 1> D = Matrix<double, N, 1>::Zero();
    Matrix<double, N, 1> Dinv = Matrix<double, N, 1>::Zero();

    for (int i = 0; i < n; i++) {
        unsigned int e = rand() % 20;
        unsigned int twoPow = 1U << e;

        D(i) = twoPow;
        Dinv(i) = 1.0 / twoPow;
    }

    return Dinv.asDiagonal() * A * D.asDiagonal();
}

/*
    IOFormat CleanFmt(FullPrecision, 0, ", ", "\n", "[", "]");

    std::cout << "A:" << std::endl
              << A << std::endl
              << "it:" << it << std::endl
              << "B:" << std::endl
              << B << std::endl
              << "Norm of B:" << std::endl
              << B.lpNorm<1>() << std::endl
              << "D:" << std::endl
              << D << std::endl
              << "Dinv:" << std::endl
              << Dinv << std::endl
              << "Check:" << std::endl
              << (A - D * B * Dinv).eval() << std::endl;
*/
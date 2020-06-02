#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "MatrixExponential.hpp"

#include "utils/readTSV.hpp"

using namespace std;
using namespace expokit;

#define N 4
#define M N * 3 * 2
#define TIMESTEP 0.005
#define MM 25
#define INDEX 1562

int main(int argc, const char* argv[])
{
    //     int scale = 10;
    //     if (argc > 1) {
    //         scale = std::atoi(argv[1]);
    //     }

    cout << "Running degradation on max multiplications" << endl;

    //     MatrixXd A = MatrixXd::Random(N, N) * scale;
    //     cout << "A:" << endl
    //          << A << endl;
    //     MatrixXd Atran = A.transpose();
    //     cout << "Atran:" << endl
    //          << Atran << endl;
    //     MatrixXd APos = -A * Atran;
    //     cout << "APos: " << endl
    //          << APos << endl;
    //     cout << "eig: " << endl
    //          << APos.eigenvalues() << endl;

    vector<MatrixXd> vecA = readTSV("exampleData/logStuffA", M, M);
    vector<MatrixXd> vecb = readTSV("exampleData/logStuffb", M, 1);
    vector<MatrixXd> vecxInit = readTSV("exampleData/logStuffxInit", M, 1);
    Matrix<double, M, 1> out;
    Matrix<double, M, 1> ref;

    // Building augmented state x0
    Matrix<double, M + 1, 1> x0;
    Matrix<double, M + 2, M + 2> A1 = Matrix<double, M + 2, M + 2>::Zero();
    x0 << vecxInit[INDEX], 1;

    // Building aumented matrix A1
    A1.template block<M, M>(0, 0) = vecA[INDEX];
    A1.template block<M, 1>(0, M) = vecb[INDEX];
    A1.template block<M + 1, 1>(0, M + 1) = x0;
    A1 *= TIMESTEP;

    //     cout << "A:" << endl
    //          << A1 << endl;

    MatrixXd Aref = A1.exp();
    //     cout << "Aref: " << endl
    //          << Aref << endl;
    cout << "Correct result norm: " << Aref.norm() << endl;

    MatrixExponential<double, Dynamic> expUtil(M + 2);
    MatrixXd res(M + 2, M + 2);

    for (int i = -1; i < MM; ++i) {
        expUtil.setMaxMultiplications(i);
        expUtil.compute(A1, res);

        double errNorm = (Aref - res).norm();
        cout << "Max Multiplications: " << i << " error: " << errNorm << endl;
        //    cout << res << endl;
    }

    return 0;
}
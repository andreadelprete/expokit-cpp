#include <Eigen/Core>
#include <Eigen/Eigenvalues> 


#include "MatrixExponential.hpp"

using namespace std;
using namespace expokit;

#define N 10
#define MM 20

int main (int argc, const char* argv[]) {
    int scale = 10;
    if (argc > 1) {
        scale = std::atoi(argv[1]);
    }

    cout << "Running degradation on max multiplications" << endl;

    MatrixXd A = MatrixXd::Random(N, N) * scale;
    cout << "A:" << endl
         << A << endl;
    MatrixXd Atran = A.transpose();
    cout << "Atran:" << endl
         << Atran << endl;
    MatrixXd APos = -A * Atran;
    cout << "APos: " << endl
         << APos << endl;
    cout << "eig: " << endl
         << APos.eigenvalues() << endl;

    MatrixXd Aref = APos.exp();
    cout << "Aref: " << endl
         << Aref << endl;
    cout<< "Correct result norm: " << Aref.norm() << endl;

    MatrixExponential<double, Dynamic> expUtil(N);
    MatrixXd res(N, N);


    for (int i = -1; i < MM; ++i)
    {
        expUtil.setMaxMultiplications(i);
        expUtil.compute(APos, res);
        
        double errNorm = (Aref - res).norm();
        cout << "Max Multiplications: " << i <<" error: " << errNorm << endl;
     //    cout << res << endl;
    }
    
    return 0;
}
#include <Eigen/Core>


#include "MatrixExponential.hpp"

using namespace std;
using namespace expokit;

#define N 10
#define MM 30

int main (int argc, const char* argv[]) {
    int scale = 3;
    if (argc > 1) {
        scale = std::atoi(argv[1]);
    }

    cout << "Running degradation on max multiplications" << endl;

    MatrixXd A = MatrixXd::Random(N, N) * scale;

    MatrixXd Aref = A.exp();
    cout<< "Correct result norm: " << Aref.norm() << endl;

    MatrixExponential<double, Dynamic> expUtil(N);
    MatrixXd res(N, N);


    for (int i = -1; i < MM; ++i)
    {
        expUtil.setMaxMultiplications(i);
        expUtil.compute(A, res);
        
        double errNorm = (Aref - res).norm();
        cout << "Max Multiplications: " << i <<" error: " << errNorm << endl;
    }
    
    return 0;
}
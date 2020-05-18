#include <Eigen/Core>

#include "MatExpIntegral.hpp"

#include <iostream>

#ifdef EIGEN_RUNTIME_NO_MALLOC
#define EIGEN_MALLOC_ALLOWED Eigen::internal::set_is_malloc_allowed(true);
#define EIGEN_MALLOC_NOT_ALLOWED Eigen::internal::set_is_malloc_allowed(false);
#else
#define EIGEN_MALLOC_ALLOWED
#define EIGEN_MALLOC_NOT_ALLOWED
#endif

using namespace Eigen;
using namespace std;
using namespace expokit;


int main()
{
    cout << "Start example integral of exp" << endl;

    MatrixXd A(2, 2);
    A << 1, 2, 3, 4;

    MatExpIntegral<double> util(2);

    MatrixXd out(2, 2);

    EIGEN_MALLOC_NOT_ALLOWED
    // util.computeExpIntegral(A, out, 5, 0.1);
    util.computeExpIntegral(A, out, 5);


    cout << out << endl;
}
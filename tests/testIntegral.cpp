#include "Eigen/Core"
#include <iostream>


#include "LDSUtility.hpp"

using namespace expokit;
using namespace std;

int main()
{
    cout << "Start test integral\n";
    
    LDSUtility<3> test;

    MatrixXd A = MatrixXd::Identity(3, 3);
    MatrixXd b(3, 1);
    b << 2, 2, 2;
    MatrixXd x(3, 1);
    x << 1, 1, 1;

    cout << "x(T):\n";
    cout << test.ComputeXt(A, b, x, 1);

    cout << "\n";
    return 0;
}
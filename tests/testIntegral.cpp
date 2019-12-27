#include "Eigen/Core"
#include <iostream>


#include "LDSUtility.hpp"

using namespace expokit;
using namespace std;

int main()
{
    cout << "cojone";
    
    LDSUtility<3> test;

    MatrixXd a = MatrixXd::Identity(3, 3);
    MatrixXd b(3, 1);
    b << 2, 2, 2;
    MatrixXd x(3, 1);
    x << 1, 1, 1;

    cout << "bomber\n";
    cout << test.ComputeXt(a, b, x, 5);

    cout << "\n";
    return 0;
}
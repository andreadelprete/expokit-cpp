#include <Eigen/Core>

#include "MatrixExponential.hpp"
#include "LDSUtility.hpp"

#include <iostream>


using namespace Eigen;
using namespace std;
using namespace expokit;


int main()
{
    cout << "Start test DYN Matrix Exponential\n";
    MatrixXd mRes(5, 5);
    VectorXd vRes(5);

    MatrixXd A = MatrixXd::Random(5, 5);
    VectorXd v = VectorXd::Random(5);
    
    // Dynamical usage
    MatrixExponential<double, Dynamic> test1; // Default init to dim 2
    //cout << test1.compute(A) << endl << endl; // Error
    test1.resize(5);
    test1.compute(A, mRes);
    cout << mRes << endl << endl;

    MatrixExponential<double, Dynamic> test2(5);
    test2.compute(A, mRes);
    cout << mRes << endl << endl;


    // Static usage
    MatrixExponential<double, 5> test3;
    //test3.resize(5); // Does nothing, error thrown if != 5
    test3.compute(A, mRes);
    cout << mRes << endl << endl;


    cout << "Start test DYN LDSUtility\n";
    MatrixXd b = MatrixXd::Random(5, 1);
    MatrixXd x = MatrixXd::Random(5, 1);

    // Dynamical usage
    LDSUtility<double, Dynamic> util1; // Defuailt init to dim 2
    //util1.ComputeXt(A, b, x, 1, vRes); // Error
    util1.resize(5);
    util1.ComputeXt(A, b, x, 1, vRes);
    cout << vRes << endl << endl;
    util1.ComputeIntegralXt(A, b, x, 1, vRes);
    cout << vRes << endl << endl;
    util1.ComputeDoubleIntegralXt(A, b, x, 1, vRes);
    cout << vRes << endl << endl;

    LDSUtility<double, Dynamic> util2(5);
    util1.ComputeXt(A, b, x, 1, vRes);
    cout << vRes << endl << endl;
    util1.ComputeIntegralXt(A, b, x, 1, vRes);
    cout << vRes << endl << endl;
    util1.ComputeDoubleIntegralXt(A, b, x, 1, vRes);
    cout << vRes << endl << endl;

    // Static usage
    LDSUtility<double, 5> util3;
    util1.ComputeXt(A, b, x, 1, vRes);
    cout << vRes << endl << endl;
    util1.ComputeIntegralXt(A, b, x, 1, vRes);
    cout << vRes << endl << endl;
    util1.ComputeDoubleIntegralXt(A, b, x, 1, vRes);
    cout << vRes << endl << endl;
    
}
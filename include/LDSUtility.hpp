#ifndef INTEGRAL_UTILITY
#define INTEGRAL_UTILITY

#include <Eigen/Core>
#include "MatrixExponential.h"
#include <iostream>
// #include <Eigen/Eigenvalues>
#include <Eigen/LU>
// #include "unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h"

using namespace Eigen;
using namespace std;

namespace expokit {

/*
    The idea behing templating is to instanciate a number of these
    classes, from 1 to N, where N is the number of active contacts 
*/
template <int N>
class LDSUtility {
private:
public:
    LDSUtility();
    ~LDSUtility();

    typedef const Eigen::Ref<const Matrix<double, N, 1> > RefVector;
    typedef const Eigen::Ref<const Matrix<double, N, N> > RefMatrix;

    // resize instead of new object

    /**
     * Computes the value of the Ax+b system at time T with initial conditions xInit
     */
    //VectorXd Xt(MatrixXd A, VectorXd b, VectorXd xInit, double T, double dt = 0, bool AInvertible = false);
    //RefVector ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, double T);
    Matrix<double, N, 1> ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, double T);


    /**
     * Computes the integral of the Ax+b system at time T with initial conditions xInit
     */
    //VectorXd IntegralXt(MatrixXd A, VectorXd b, VectorXd xInit, double T, double dt = 0, bool AInvertible = false);
    // RefVector& ComputeIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, double T);

    /**
     * Computes double the integral of the Ax+b system at time T with initial conditions xInit
     */
    //VectorXd DoubleIntegralXt(MatrixXd A, VectorXd b, VectorXd xInit, double T, double dt = 0, bool AInvertible = false);
    // RefVector& ComputeDoubleIntegralXt(RefMatrix& A, RefVector& b, RefVector& xInit, double T);

    /**
     * Stub proposal for computing in parallel the quantities above
     */
    //MatrixXd ComputeParallel(Vector3d what, MatrixXd A, VectorXd b, VectorXd xInit, double T, double dt = 0, bool AInvertible = false);
    //MatrixXd ComputeParallel(Vector3d what, RefMatrix &A, RefVector &b, RefVector &xInit, double T);
};

template <int N>
LDSUtility<N>::LDSUtility()
{
}

template <int N>
LDSUtility<N>::~LDSUtility()
{
}

template <int N>
//typename LDSUtility<N>::RefVector LDSUtility<N>::ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, double T)
Matrix<double, N, 1> LDSUtility<N>::ComputeXt(RefMatrix& A, RefVector& b, RefVector& xInit, double T)
{
    // Building aumented matrix A0
    Matrix<double, N, N + 1> A0temp;
    A0temp << A, b;
    Matrix<double, N + 1, N + 1> A0 = Matrix<double, N + 1, N + 1>::Zero();
    A0 << A0temp, Matrix<double, 1, N + 1>::Zero();
    A0 *= T;
    cout << "A0:---->\n" << A0 << "\n\n";
    
    // Building augmented state x0
    Matrix<double, N + 1, 1> x0;
    x0 << xInit, 1;

    // Building filtering matrix z0
    // TODO - couldn't we just extract a block from the result?
    Matrix<double, N, N + 1> z0;
    z0 << Matrix<double, N, N>::Identity(), Matrix<double, N, 1>::Zero();

    // Temp matrix
    Matrix<double, N + 1, 1> xTemp;

    // Matrix exponential - vectoring for now fixed to an arbitrary number
    MatrixExponential<Matrix<double, N + 1, N + 1>, Matrix<double, N + 1, 1> > expUtil(N+1);    
    expUtil.computeExpTimesVector(A0, x0, xTemp, 4);
    cout << "xTemp:---->\n" << xTemp << "\n\n";

    // Test whole matrix
    Matrix<double, N + 1, N + 1> xTemp2;
    expUtil.compute(A0, xTemp2);
    cout << "xTemp2:---->\n" << xTemp2 << "\n\n";

    // Result matrix 
    Matrix<double, N, 1> x = z0 * xTemp;

    return x;
}

}

#endif
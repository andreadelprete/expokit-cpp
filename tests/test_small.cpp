/* test_dense_general.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

// include Eigen before f2c to avoid compilation errors!
//#define EIGEN_RUNTIME_NO_MALLOC
//#define EIGEN_NO_MALLOC
#include <Eigen/Core>
//#include <unsupported/Eigen/MatrixFunctions>

#include <cfloat>
#include <list>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
//#include "unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h"
#include "MatrixExponential.h"

#ifdef EIGEN_RUNTIME_NO_MALLOC
  #define EIGEN_MALLOC_ALLOWED Eigen::internal::set_is_malloc_allowed(true);
  #define EIGEN_MALLOC_NOT_ALLOWED Eigen::internal::set_is_malloc_allowed(false);
#else
  #define EIGEN_MALLOC_ALLOWED
  #define EIGEN_MALLOC_NOT_ALLOWED
#endif

#ifdef __cplusplus
extern "C" {

#include "dgpadm.h"
#include "dgexpv.h"
#include "dgchbv.h"
#include "clock.h"

// you need to define MAIN__ function because f2c complains if it doesn't find one
int  MAIN__( ) {  return 0; }

}
#endif

#include <stdio.h>
#include <iostream>


using namespace Eigen;
using namespace std;

/** List of available parameters of IOFormat constructor:
precision       number of digits for floating point values, or one of the special constants StreamPrecision and FullPrecision.
flags           either 0, or DontAlignCols, which allows to disable the alignment of columns, resulting in faster code.
coeffSeparator  string printed between two coefficients of the same row
rowSeparator    string printed between two rows
rowPrefix       string printed at the beginning of each row
rowSuffix       string printed at the end of each row
matPrefix       string printed at the beginning of the matrix
matSuffix       string printed at the end of the matrix */
static const Eigen::IOFormat CleanFmt(5, 0, ", ", "\n", "[", "]");

#define MAX_PRINT_N 10
#define PRINT_VECTOR(a) if(a.cols()<=MAX_PRINT_N) std::cout<<#a<<"("<<a.rows()<<"x"<<a.cols()<<"): "<<a.transpose().format(CleanFmt)<<std::endl
#define PRINT_MATRIX(a) if(a.cols()<=MAX_PRINT_N) std::cout<<#a<<"("<<a.rows()<<"x"<<a.cols()<<"):\n"<<a.format(CleanFmt)<<std::endl

/* Computes y = A*x */
int mat_vec_mul(double *x, double *y)
{
    cout<<"mat vec mul\n";
    for(int i=0; i<5; i++)
        y[i] = x[i];
    return 0;
}


int main( int argc, const char* argv[] )
{
    printf( "\nStart test_small\n" );

    bool TEST_DGEXPV = false;
    int n_tests = 1000;
    int n_contacts = 1;
    int m = n_contacts*3*2;      // matrix size
    double t = .005;     // time step

    printf("Number of contact points %d\n\n", n_contacts);
    printf("Matrix size %d\n\n", m);

    int m2 = int(m/2);
    double stiffness = 1e5;
    double damping = 1e2;
    MatrixXd U  = MatrixXd::Random(m2, m2);
    MatrixXd Upsilon = U*U.transpose();
    MatrixXd K = MatrixXd::Identity(m2,m2)*stiffness;
    MatrixXd B = MatrixXd::Identity(m2,m2)*damping;
    MatrixXd A = MatrixXd::Zero(m, m);
    A.topRightCorner(m2,m2) = MatrixXd::Identity(m2,m2);
    A.bottomLeftCorner(m2,m2) = -Upsilon*K;
    A.bottomRightCorner(m2,m2) = -Upsilon*B;
//    A *= t;
//    MatrixXd A = matlib.block([[matlib.zeros((n2,n2)), matlib.eye(n2)],
//                      [             -Upsilon*K,      -Upsilon*B]])

    PRINT_MATRIX(A);

/* ---  Krylov-based method (for sparse matrices) ... */
    if(TEST_DGEXPV)
    {
        int mk = (m<31) ? m-2 : 30; //maximum size for the Krylov basis
        double tol = 0.;
        VectorXd y = VectorXd::Zero(m);
        y[0] = 1.;
        VectorXd res = VectorXd::Zero(m);
        double Anorm = m*0.5; // an approximation of some norm of A
        int lwsp = m*(mk+2) + 5*(mk+2)*(mk+2) + 7; // length of workspace
        VectorXd wsp(lwsp);
        int liwsp = mk+2;
        VectorXi iwsp(liwsp);
        int itrace = 1; // running mode. 0=silent, 1=print step-by-step info */
        int iflag;
        dgexpv_(&m, &mk, &t, y.data(), res.data(), &tol, &Anorm, 
                wsp.data(), &lwsp, iwsp.data(), &liwsp, (U_fp)mat_vec_mul, &itrace, &iflag);
        cout<<"end dgexpv_\n";
        cout<<"exp(t*A)e_1="<<y.transpose().format(CleanFmt)<<endl;
    }

/* ---  Pade ... */
    int ideg = 6;
    int lwsp = 4*m*m+ideg+1;
    VectorXd wsp(lwsp);
    VectorXi iwsp(m);
    int ns, iexp, iflag;
    double tic = clock_();
    for(int i=0;i<n_tests;i++)
    {
        //A = MatrixXd::Random(m, m);
        dgpadm_(&ideg, &m, &t, A.data(), &m, wsp.data(), &lwsp, iwsp.data(), &iexp, &ns, &iflag);
    }
    double tac = clock_();

    if(iflag==0)
        printf("\nWith DGPADM everything went fine:\n");
    else
        printf("\nWith DGPADM there was a problem, iflag=%d\n", iflag);
    
    if(m<=MAX_PRINT_N)
    {
        printf("exp(t*A) =\n");
        for (int i = 1; i <= m; ++i) {
	        for (int j = 1; j <= m; ++j) {
	            printf("%f ", wsp[iexp + (j - 1) * m + i - 2]);
    	    }
            printf("\n");
        }
    }
    printf("run time = %.3f ms\n", 1e3*(tac-tic)/n_tests);
    printf("Number of squaring: %d\n", ns);


/* ---  Chebyshev */
    VectorXd y = VectorXd::Zero(m);
    doublecomplex cwsp[10007];
    tic = clock_();
    for(int i=0;i<n_tests;i++)
    {
        y.tail(m-1).setZero();
        y[0] = 1.0;
        dgchbv_(&m, &t, A.data(), &m, y.data(), cwsp, iwsp.data(), &iflag);
    }
    tac = clock_();
    if(iflag==0)
        printf("\nWith DGCHBV everything went fine:\n");
    else
        printf("\nWith DGCHBV there was a problem, iflag=%d\n", iflag);
    cout<< "exp(t*A)e_1 =\n";
    PRINT_VECTOR(y);
    printf("run time = %.3f ms\n", 1e3*(tac-tic)/n_tests);


    /* Eigen Pade with scaling and squaring */
    MatrixXd expA(m,m);
    A *= t;
//    EIGEN_MALLOC_NOT_ALLOWED

    tic = clock_();
    for(int i=0;i<n_tests;i++)
    {
        expA = A.exp();
//        matrix_exp_compute(A, expA);
    }
    tac = clock_();

    printf("\nWith Eigen everything went fine:\n");
    printf("run time = %.3f ms\n", 1e3*(tac-tic)/n_tests);
    if(m<=MAX_PRINT_N)
        PRINT_MATRIX(expA);

    printf( "End test_small\n");
}


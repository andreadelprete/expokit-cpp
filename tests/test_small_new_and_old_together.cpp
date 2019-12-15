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
#include <Eigen/Core>

#ifdef __cplusplus
extern "C" {

#include "dgpadm.h"
#include "dgexpv.h"
#include "dgchbv.h"
#include "clock.c"

// you need to define MAIN__ function because f2c complains if it doesn't find one
int  MAIN__( ) {  return 0; }

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__6 = 6;
static integer c__50 = 50;
static integer c__10007 = 10007;

int test_small_old(double *a, int m, double t)
{
    /* Format strings */
    static char fmt_9000[] = "(/,a,/,a)";
    static char fmt_9001[] = "(5(1x,d11.4))";
    static char fmt_9002[] = "(a,e10.5)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static doublereal h__[2500];
    static integer i__, j, k;
    static doublereal y[50], tic, tac;
    static integer ns;
    static doublereal wsp[10007];
    static integer iexp, iwsp[50], iflag, iseed[4];
    extern doublereal clock_(void);
    extern /* Subroutine */ int dgpadm_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *), dgchbv_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    extern doublereal dlaran_(integer *);
    extern /* Subroutine */ int dspadm_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *), dschbv_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static integer mprint;

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 6, 0, 0, 0 };
    static cilist io___9 = { 0, 6, 0, 0, 0 };
    static cilist io___10 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___18 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___19 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___20 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___21 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___23 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___26 = { 0, 6, 0, 0, 0 };
    static cilist io___runtime = { 0, 6, 0, fmt_9002, 0 };

    t = .5;
    m = 5;
//    a = A;
    
/* ---  maximum number of rows/columns to be printed */
    mprint = min(5,m);
    s_wsle(&io___8);
    do_lio(&c__9, &c__1, "t =", (ftnlen)3);
    e_wsle();
    s_wsle(&io___9);
    do_lio(&c__5, &c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    e_wsle();
    s_wsfe(&io___10);
    do_fio(&c__1, "REAL SYMMETRIC CASE", (ftnlen)19);
    do_fio(&c__1, "*******************************", (ftnlen)31);
    e_wsfe();
    s_wsle(&io___11);
    do_lio(&c__9, &c__1, "A = ", (ftnlen)4);
    e_wsle();
    s_wsfe(&io___12);
    i__1 = mprint;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = mprint;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&a[i__ + j * m - m-1], (ftnlen)sizeof(
		    doublereal));
	}
    }
    e_wsfe();

/* ---  Some compliers (e.g., g77) generate 'Unsupported FORMAT specifier' */
/*     with the specification above. In this case, simply use this form: */
/* 9001 format( 5(1X,D11.4) ) */
/* ---  Pade ... */
    tic = clock_();
    dgpadm_(&c__6, &m, &t, a, &m, wsp, &c__10007, iwsp, &iexp, &ns, &iflag);
    tac = clock_();
    s_wsfe(&io___18);
    do_fio(&c__1, "With DGPADM:", (ftnlen)12);
    do_fio(&c__1, "exp(t*A) =", (ftnlen)10);
    e_wsfe();
    s_wsfe(&io___19);
    i__2 = mprint;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = mprint;
	for (j = 1; j <= i__1; ++j) {
	    do_fio(&c__1, (char *)&wsp[iexp + (j - 1) * m + i__ - 2], (ftnlen)sizeof(doublereal));
	}
    }
    e_wsfe();
    s_wsfe(&io___runtime);
    do_fio(&c__1, "runtime   = ", (ftnlen)12);
    d__1 = tac - tic;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();

/* ---  Chebyshev */
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	y[i__ - 1] = 0.;
    }
    y[0] = 1.;
    tic = clock_();
    dgchbv_(&m, &t, a, &m, y, wsp, iwsp, &iflag);
    tac = clock_();
    s_wsfe(&io___23);
    do_fio(&c__1, "With DGCHBV:", (ftnlen)12);
    do_fio(&c__1, "exp(t*A)e_1 =", (ftnlen)13);
    e_wsfe();
    i__2 = mprint;
    for (i__ = 1; i__ <= i__2; ++i__) {
	s_wsle(&io___24);
	do_lio(&c__5, &c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(doublereal));
	e_wsle();
    }
    s_wsfe(&io___runtime);
    do_fio(&c__1, "runtime   = ", (ftnlen)12);
    d__1 = tac - tic;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();

    return 0;
} /* MAIN__ */

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
    int n_tests = 100;
    int m = 5;	        // matrix size
    double t = .5;    // time step

    printf("Matrix size %d\n\n", m);

    MatrixXd A = MatrixXd::Random(m, m);
    A = -A*A.transpose();
    //A = MatrixXd::Identity(m,m);
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

    printf( "End test_small\n");
    
    test_small_old(A.data(), m, t);
}


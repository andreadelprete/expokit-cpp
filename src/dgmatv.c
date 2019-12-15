/* dgmatv.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal a[600000];
    integer ia[600000], ja[600000], nz, n;
} rmat_;

#define rmat_1 rmat_

/* -------------------------------NOTE-----------------------------------* */
/*     This is an accessory to Expokit and it is not intended to be     * */
/*     complete. It is supplied primarily to ensure an unconstrained    * */
/*     distribution and portability of the package. The matrix-vector   * */
/*     multiplication routines supplied here fit the non symmetric      * */
/*     storage and for a symmetric matrix, the entire (not half) matrix * */
/*     is required.  If the sparsity pattern is known a priori, it is   * */
/*     recommended to use the most advantageous format and to devise    * */
/*     the most advantageous matrix-vector multiplication routine.      * */
/* ----------------------------------------------------------------------* */
/* ----------------------------------------------------------------------* */
/* Subroutine */ int dgcoov_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;


/* ---  Computes y = A*x. A is passed via a fortran `common statement'. */
/* ---  A is assumed here to be under the COOrdinates storage format. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = rmat_1.n;
    for (j = 1; j <= i__1; ++j) {
	y[j] = 0.;
    }
    i__1 = rmat_1.nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[rmat_1.ia[i__ - 1]] += rmat_1.a[i__ - 1] * x[rmat_1.ja[i__ - 1]];
    }
    return 0;
} /* dgcoov_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* Subroutine */ int dgcrsv_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;


/* ---  Computes y = A*x. A is passed via a fortran `common statement'. */
/* ---  A is assumed to be under the Compress Row Storage (CRS) format. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = 0.;
	i__2 = rmat_1.ia[i__] - 1;
	for (j = rmat_1.ia[i__ - 1]; j <= i__2; ++j) {
	    y[i__] += rmat_1.a[j - 1] * x[rmat_1.ja[j - 1]];
	}
    }
    return 0;
} /* dgcrsv_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* Subroutine */ int dgccsv_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;


/* ---  Computes y = A*x. A is passed via a fortran `common statement'. */
/* ---  A is assumed to be under the Compress Column Storage (CCS) format. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = 0.;
    }
    i__1 = rmat_1.n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = rmat_1.ja[j] - 1;
	for (i__ = rmat_1.ja[j - 1]; i__ <= i__2; ++i__) {
	    y[rmat_1.ia[i__ - 1]] += rmat_1.a[i__ - 1] * x[j];
	}
    }
    return 0;
} /* dgccsv_ */


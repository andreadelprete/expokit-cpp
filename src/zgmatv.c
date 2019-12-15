/* zgmatv.f -- translated by f2c (version 20100827).
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
    doublecomplex a[50000];
    integer ia[50000], ja[50000], nz, n;
} cmat_;

#define cmat_1 cmat_

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
/* Subroutine */ int zgcoov_(doublecomplex *x, doublecomplex *y)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, j;


/* ---  Computes y = A*x. A is passed via a fortran `common statement'. */
/* ---  A is assumed here to be under the COOrdinates storage format. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = cmat_1.n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	y[i__2].r = 0., y[i__2].i = 0.;
    }
    i__1 = cmat_1.nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = cmat_1.ia[i__ - 1];
	i__3 = cmat_1.ia[i__ - 1];
	i__4 = i__ - 1;
	i__5 = cmat_1.ja[i__ - 1];
	z__2.r = cmat_1.a[i__4].r * x[i__5].r - cmat_1.a[i__4].i * x[i__5].i, 
		z__2.i = cmat_1.a[i__4].r * x[i__5].i + cmat_1.a[i__4].i * x[
		i__5].r;
	z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
	y[i__2].r = z__1.r, y[i__2].i = z__1.i;
    }
    return 0;
} /* zgcoov_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* Subroutine */ int zgcrsv_(doublecomplex *x, doublecomplex *y)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, j;


/* ---  Computes y = A*x. A is passed via a fortran `common statement'. */
/* ---  A is assumed to be under the Compress Row Storage (CRS) format. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = cmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	y[i__2].r = 0., y[i__2].i = 0.;
	i__2 = cmat_1.ia[i__] - 1;
	for (j = cmat_1.ia[i__ - 1]; j <= i__2; ++j) {
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = j - 1;
	    i__6 = cmat_1.ja[j - 1];
	    z__2.r = cmat_1.a[i__5].r * x[i__6].r - cmat_1.a[i__5].i * x[i__6]
		    .i, z__2.i = cmat_1.a[i__5].r * x[i__6].i + cmat_1.a[i__5]
		    .i * x[i__6].r;
	    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
	    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
	}
    }
    return 0;
} /* zgcrsv_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* Subroutine */ int zgccsv_(doublecomplex *x, doublecomplex *y)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, j;


/* ---  Computes y = A*x. A is passed via a fortran `common statement'. */
/* ---  A is assumed to be under the Compress Column Storage (CCS) format. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = cmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	y[i__2].r = 0., y[i__2].i = 0.;
    }
    i__1 = cmat_1.n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = cmat_1.ja[j] - 1;
	for (i__ = cmat_1.ja[j - 1]; i__ <= i__2; ++i__) {
	    i__3 = cmat_1.ia[i__ - 1];
	    i__4 = cmat_1.ia[i__ - 1];
	    i__5 = i__ - 1;
	    i__6 = j;
	    z__2.r = cmat_1.a[i__5].r * x[i__6].r - cmat_1.a[i__5].i * x[i__6]
		    .i, z__2.i = cmat_1.a[i__5].r * x[i__6].i + cmat_1.a[i__5]
		    .i * x[i__6].r;
	    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
	    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
	}
    }
    return 0;
} /* zgccsv_ */


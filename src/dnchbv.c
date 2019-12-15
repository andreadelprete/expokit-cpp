/* dnchbv.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static integer c__1 = 1;

/* ----------------------------------------------------------------------| */
/* Subroutine */ int dnchbv_(integer *m, doublereal *t, doublereal *h__, 
	integer *ldh, doublereal *y, doublecomplex *wsp)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, ih, ip, iy, iz;
    static doublecomplex tmpc, alpha[7], theta[7];
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublereal alpha0;
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);

/* -----Purpose----------------------------------------------------------| */

/* ---  DNCHBV computes y = exp(t*H)*y using the partial fraction */
/*     expansion of the uniform rational Chebyshev approximation */
/*     to exp(-x) of type (14,14). H is assumed to be upper-Hessenberg. */
/*     About 14-digit accuracy is expected if the matrix H is negative */
/*     definite. The algorithm may behave poorly otherwise. */

/* -----Arguments--------------------------------------------------------| */

/*     m       : (input) order of the Hessenberg matrix H */

/*     t       : (input) time-scaling factor (can be < 0). */

/*     H(ldh,m): (input) upper Hessenberg matrix. */

/*     y(m)    : (input/output) on input the operand vector, */
/*               on output the resulting vector exp(t*H)*y. */

/*     wsp     : (workspace). Observe that a double precision vector of */
/*               length 2*m*(m+2) can be used as well when calling this */
/*               routine (thus avoiding an idle complex array elsewhere) */

/* ----------------------------------------------------------------------| */
/*     Roger B. Sidje (rbs@maths.uq.edu.au) */
/*     EXPOKIT: Software Package for Computing Matrix Exponentials. */
/*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998 */
/* ----------------------------------------------------------------------| */

/* ---  Pointers ... */
    /* Parameter adjustments */
    --wsp;
    --y;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;

    /* Function Body */
    ih = 1;
    iy = ih + *m * *m;
    iz = iy + *m;
/* ---  Coefficients and poles of the partial fraction expansion... */
    alpha0 = 1.83216998528140087e-12;
    alpha[0].r = 55.7503973136501826, alpha[0].i = -204.295038779771857;
    alpha[1].r = -93.8666838877006739, alpha[1].i = 91.2874896775456363;
    alpha[2].r = 46.9965415550370835, alpha[2].i = -11.6167609985818103;
    alpha[3].r = -9.61424200626061065, alpha[3].i = -2.64195613880262669;
    alpha[4].r = .752722063978321642, alpha[4].i = .67036736556637777;
    alpha[5].r = -.0188781253158648576, alpha[5].i = -.0343696176445802414;
    alpha[6].r = 1.43086431411801849e-4, alpha[6].i = 2.87221133228814096e-4;
    theta[0].r = -5.62314417475317895, theta[0].i = 1.1940692161124744;
    theta[1].r = -5.0893467972821611, theta[1].i = 3.58882439228376881;
    theta[2].r = -3.99337136365302569, theta[2].i = 6.00483209099604664;
    theta[3].r = -2.26978543095856366, theta[3].i = 8.46173881758693369;
    theta[4].r = .208756929753827868, theta[4].i = 10.9912615662209418;
    theta[5].r = 3.70327340957595652, theta[5].i = 13.6563731924991884;
    theta[6].r = 8.89777151877331107, theta[6].i = 16.6309842834712071;

/* ---  Accumulation of the contribution of each pole ... */

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = iz + j - 1;
	i__3 = j;
	wsp[i__2].r = y[i__3], wsp[i__2].i = 0.;
	y[j] *= alpha0;
    }
    for (ip = 1; ip <= 7; ++ip) {
/* ---     Solve each fraction using Gaussian elimination with pivoting... */
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = iy + j - 1;
	    i__3 = iz + j - 1;
	    wsp[i__2].r = wsp[i__3].r, wsp[i__2].i = wsp[i__3].i;
/* Computing MIN */
	    i__3 = j + 1;
	    i__2 = min(i__3,*m);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = ih + (j - 1) * *m + i__ - 1;
		d__1 = -(*t) * h__[i__ + j * h_dim1];
		wsp[i__3].r = d__1, wsp[i__3].i = 0.;
	    }
	    i__2 = ih + (j - 1) * *m + j - 1;
	    i__3 = ih + (j - 1) * *m + j - 1;
	    i__4 = ip - 1;
	    z__1.r = wsp[i__3].r - theta[i__4].r, z__1.i = wsp[i__3].i - 
		    theta[i__4].i;
	    wsp[i__2].r = z__1.r, wsp[i__2].i = z__1.i;
	    i__2 = *m;
	    for (k = i__; k <= i__2; ++k) {
		i__3 = ih + (j - 1) * *m + k - 1;
		wsp[i__3].r = 0., wsp[i__3].i = 0.;
	    }
	}
	i__1 = *m - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* ---        Get pivot and exchange rows ... */
	    if (z_abs(&wsp[ih + (i__ - 1) * *m + i__ - 1]) < z_abs(&wsp[ih + (
		    i__ - 1) * *m + i__])) {
		i__2 = *m - i__ + 1;
		zswap_(&i__2, &wsp[ih + (i__ - 1) * *m + i__ - 1], m, &wsp[ih 
			+ (i__ - 1) * *m + i__], m);
		zswap_(&c__1, &wsp[iy + i__ - 1], &c__1, &wsp[iy + i__], &
			c__1);
	    }
/* ---        Forward eliminiation ... */
	    z_div(&z__1, &wsp[ih + (i__ - 1) * *m + i__], &wsp[ih + (i__ - 1) 
		    * *m + i__ - 1]);
	    tmpc.r = z__1.r, tmpc.i = z__1.i;
	    i__2 = *m - i__;
	    z__1.r = -tmpc.r, z__1.i = -tmpc.i;
	    zaxpy_(&i__2, &z__1, &wsp[ih + i__ * *m + i__ - 1], m, &wsp[ih + 
		    i__ * *m + i__], m);
	    i__2 = iy + i__;
	    i__3 = iy + i__;
	    i__4 = iy + i__ - 1;
	    z__2.r = tmpc.r * wsp[i__4].r - tmpc.i * wsp[i__4].i, z__2.i = 
		    tmpc.r * wsp[i__4].i + tmpc.i * wsp[i__4].r;
	    z__1.r = wsp[i__3].r - z__2.r, z__1.i = wsp[i__3].i - z__2.i;
	    wsp[i__2].r = z__1.r, wsp[i__2].i = z__1.i;
	}
/* ---     Backward substitution ... */
	for (i__ = *m; i__ >= 1; --i__) {
	    i__1 = iy + i__ - 1;
	    tmpc.r = wsp[i__1].r, tmpc.i = wsp[i__1].i;
	    i__1 = *m;
	    for (j = i__ + 1; j <= i__1; ++j) {
		i__2 = ih + (j - 1) * *m + i__ - 1;
		i__3 = iy + j - 1;
		z__2.r = wsp[i__2].r * wsp[i__3].r - wsp[i__2].i * wsp[i__3]
			.i, z__2.i = wsp[i__2].r * wsp[i__3].i + wsp[i__2].i *
			 wsp[i__3].r;
		z__1.r = tmpc.r - z__2.r, z__1.i = tmpc.i - z__2.i;
		tmpc.r = z__1.r, tmpc.i = z__1.i;
	    }
	    i__1 = iy + i__ - 1;
	    z_div(&z__1, &tmpc, &wsp[ih + (i__ - 1) * *m + i__ - 1]);
	    wsp[i__1].r = z__1.r, wsp[i__1].i = z__1.i;
	}
/* ---     Accumulate the partial result in y ... */
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = ip - 1;
	    i__3 = iy + j - 1;
	    z__1.r = alpha[i__2].r * wsp[i__3].r - alpha[i__2].i * wsp[i__3]
		    .i, z__1.i = alpha[i__2].r * wsp[i__3].i + alpha[i__2].i *
		     wsp[i__3].r;
	    y[j] += z__1.r;
	}
    }
    return 0;
} /* dnchbv_ */


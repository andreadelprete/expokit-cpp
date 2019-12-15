/* dgchbv.f -- translated by f2c (version 20100827).
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

/* ---  DGCHBV computes y = exp(t*H)*y using the partial fraction */
/*     expansion of the uniform rational Chebyshev approximation */
/*     to exp(-x) of type (14,14). H is a General matrix. */
/*     About 14-digit accuracy is expected if the matrix H is negative */
/*     definite. The algorithm may behave poorly otherwise. */

/* -----Arguments--------------------------------------------------------| */

/*     m       : (input) order of the matrix H */

/*     t       : (input) time-scaling factor (can be < 0). */

/*     H(ldh,m): (input) argument matrix. */

/*     y(m)    : (input/output) on input the operand vector, */
/*               on output the resulting vector exp(t*H)*y. */

/*     iwsp(m) : (workspace) */

/*     wsp     : (workspace). Observe that a double precision vector of */
/*               length 2*m*(m+2) can be used as well when calling this */
/*               routine (thus avoiding an idle complex array elsewhere) */

/* ----------------------------------------------------------------------| */
/*     Roger B. Sidje (rbs@maths.uq.edu.au) */
/*     EXPOKIT: Software Package for Computing Matrix Exponentials. */
/*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998 */
/* ----------------------------------------------------------------------| */

int dgchbv_(integer *m, doublereal *t, doublereal *h__, 
	integer *ldh, doublereal *y, doublecomplex *wsp, integer *iwsp, 
	integer *iflag);


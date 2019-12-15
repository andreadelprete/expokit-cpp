/* dgpadm.f -- translated by f2c (version 20100827).
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

/*     Computes exp(t*H), the matrix exponential of a general matrix in */
/*     full, using the irreducible rational Pade approximation to the */
/*     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ), */
/*     combined with scaling-and-squaring. */

/* -----Arguments--------------------------------------------------------| */

/*     ideg      : (input) the degre of the diagonal Pade to be used. */
/*                 a value of 6 is generally satisfactory. */

/*     m         : (input) order of H. */

/*     H(ldh,m)  : (input) argument matrix. */

/*     t         : (input) time-scale (can be < 0). */

/*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1. */

/*     ipiv(m)   : (workspace) */

/* >>>> iexph     : (output) number such that wsp(iexph) points to exp(tH) */
/*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1) */
/*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
/*                 NOTE: if the routine was called with wsp(iptr), */
/*                       then exp(tH) will start at wsp(iptr+iexph-1). */

/*     ns        : (output) number of scaling-squaring used. */

/*     iflag     : (output) exit flag. */
/*                      0 - no problem */
/*                     <0 - problem */

/* ----------------------------------------------------------------------| */
/*     Roger B. Sidje (rbs@maths.uq.edu.au) */
/*     EXPOKIT: Software Package for Computing Matrix Exponentials. */
/*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998 */
/* ----------------------------------------------------------------------| */
int dgpadm_(integer *ideg, integer *m, doublereal *t, 
	doublereal *h__, integer *ldh, doublereal *wsp, integer *lwsp, 
	integer *ipiv, integer *iexph, integer *ns, integer *iflag);


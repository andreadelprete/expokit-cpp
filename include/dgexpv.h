/* dgexpv.f -- translated by f2c (version 20100827).
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

/* ----------------------------------------------------------------------| */
/* Subroutine */ int dgexpv_(integer *n, integer *m, doublereal *t, 
	doublereal *v, doublereal *w, doublereal *tol, doublereal *anorm, 
	doublereal *wsp, integer *lwsp, integer *iwsp, integer *liwsp, S_fp 
	matvec, integer *itrace, integer *iflag);

/* ---  DGEXPV computes w = exp(t*A)*v - for a General matrix A. */

/*     It does not compute the matrix exponential in isolation but */
/*     instead, it computes directly the action of the exponential */
/*     operator on the operand vector. This way of doing so allows */
/*     for addressing large sparse problems. */

/*     The method used is based on Krylov subspace projection */
/*     techniques and the matrix under consideration interacts only */
/*     via the external routine `matvec' performing the matrix-vector */
/*     product (matrix-free method). */

/* -----Arguments--------------------------------------------------------| */

/*     n      : (input) order of the principal matrix A. */

/*     m      : (input) maximum size for the Krylov basis. */

/*     t      : (input) time at wich the solution is needed (can be < 0). */

/*     v(n)   : (input) given operand vector. */

/*     w(n)   : (output) computed approximation of exp(t*A)*v. */

/*     tol    : (input/output) the requested accuracy tolerance on w. */
/*              If on input tol=0.0d0 or tol is too small (tol.le.eps) */
/*              the internal value sqrt(eps) is used, and tol is set to */
/*              sqrt(eps) on output (`eps' denotes the machine epsilon). */
/*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol) */

/*     anorm  : (input) an approximation of some norm of A. */

/*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1 */
/*                                   +---------+-------+---------------+ */
/*              (actually, ideg=6)        V        H      wsp for PADE */

/* iwsp(liwsp): (workspace) liwsp .ge. m+2 */

/*     matvec : external subroutine for matrix-vector multiplication. */
/*              synopsis: matvec( x, y ) */
/*                        double precision x(*), y(*) */
/*              computes: y(1:n) <- A*x(1:n) */
/*                        where A is the principal matrix. */

/*     itrace : (input) running mode. 0=silent, 1=print step-by-step info */

/*     iflag  : (output) exit flag. */
/*              <0 - bad input arguments */
/*               0 - no problem */
/*               1 - maximum number of steps reached without convergence */
/*               2 - requested tolerance was too high */

/* -----Accounts on the computation--------------------------------------| */
/*     Upon exit, an interested user may retrieve accounts on the */
/*     computations. They are located in wsp and iwsp as indicated below: */

/*     location  mnemonic                 description */
/*     -----------------------------------------------------------------| */
/*     iwsp(1) = nmult, number of matrix-vector multiplications used */
/*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated */
/*     iwsp(3) = nscale, number of repeated squaring involved in Pade */
/*     iwsp(4) = nstep, number of integration steps used up to completion */
/*     iwsp(5) = nreject, number of rejected step-sizes */
/*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise */
/*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured */
/*     -----------------------------------------------------------------| */
/*     wsp(1)  = step_min, minimum step-size used during integration */
/*     wsp(2)  = step_max, maximum step-size used during integration */
/*     wsp(3)  = dummy */
/*     wsp(4)  = dummy */
/*     wsp(5)  = x_error, maximum among all local truncation errors */
/*     wsp(6)  = s_error, global sum of local truncation errors */
/*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured */
/*     wsp(8)  = t_now, integration domain successfully covered */
/*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0) */
/*     wsp(10) = ||w||/||v||, scaled norm of the solution w. */
/*     -----------------------------------------------------------------| */
/*     The `hump' is a measure of the conditioning of the problem. The */
/*     matrix exponential is well-conditioned if hump = 1, whereas it is */
/*     poorly-conditioned if hump >> 1. However the solution can still be */
/*     relatively fairly accurate even when the hump is large (the hump */
/*     is an upper bound), especially when the hump and the scaled norm */
/*     of w [this is also computed and returned in wsp(10)] are of the */
/*     same order of magnitude (further details in reference below). */

/* ----------------------------------------------------------------------| */
/* -----The following parameters may also be adjusted herein-------------| */

/*     mxstep  : maximum allowable number of integration steps. */
/*               The value 0 means an infinite number of steps. */

/*     mxreject: maximum allowable number of rejections at each step. */
/*               The value 0 means an infinite number of rejections. */

/*     ideg    : the Pade approximation of type (ideg,ideg) is used as */
/*               an approximation to exp(H). The value 0 switches to the */
/*               uniform rational Chebyshev approximation of type (14,14) */

/*     delta   : local truncation error `safety factor' */

/*     gamma   : stepsize `shrinking factor' */

/* ----------------------------------------------------------------------| */
/*     Roger B. Sidje (rbs@maths.uq.edu.au) */
/*     EXPOKIT: Software Package for Computing Matrix Exponentials. */
/*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998 */
/* ----------------------------------------------------------------------| */


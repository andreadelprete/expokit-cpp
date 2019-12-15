/* sample_p.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static integer c__1 = 1;
static integer c_b10 = 279052;
static integer c__5000 = 5000;
static integer c__5 = 5;
static doublereal c_b184 = 1.;
static doublereal c_b187 = -1.;
static integer c__9 = 9;


/* ---  sample program illustrating the use of DGPHIV ... */
/*     Nonhomogeneous problem (Example 6.5 in the Expokit report) ... */

/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_9001[] = "(a)";
    static char fmt_9002[] = "(a,e8.2)";
    static char fmt_9003[] = "(a,i9)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer i__, m;
    static doublereal t, u[5000], v[5000], w[5000], tac, tic, tol, tmp, wsp[
	    279052];
    static integer iwsp[5000];
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static integer iflag;
    extern doublereal clock_(void);
    static doublereal anorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal wsave[5000];
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), loadhb_(char *, char *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, ftnlen, ftnlen);
    static integer itrace;
    extern /* Subroutine */ int dgccsv_(doublereal *, doublereal *), dgphiv_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, S_fp, integer *, integer *), dgcnvr_(char *
	    , char *, char *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 6, 0, "(A,E8.2)", 0 };
    static cilist io___16 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___17 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___18 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___19 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___23 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___24 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___26 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___27 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___35 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___36 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___37 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___40 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___41 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___42 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___43 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___44 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___46 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___47 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___49 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___50 = { 0, 6, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___52 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___53 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___54 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___55 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___56 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___57 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___58 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___59 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___60 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___61 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___62 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___63 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___64 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___65 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___66 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___67 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___68 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___69 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___70 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___71 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___72 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___73 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___74 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___76 = { 0, 6, 0, 0, 0 };
    static cilist io___77 = { 0, 6, 0, 0, 0 };
    static cilist io___78 = { 0, 6, 0, 0, 0 };
    static cilist io___79 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___80 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___81 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___82 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___83 = { 0, 6, 0, 0, 0 };
    static cilist io___84 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___85 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___86 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___87 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___88 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___89 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___90 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___91 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___92 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___93 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___94 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___95 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___96 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___97 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___98 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___99 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___100 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___101 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___102 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___103 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___104 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___105 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___106 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___107 = { 0, 6, 0, fmt_9002, 0 };


/* ---  matrix data ... */
/* ---  BEWARE: these values must match those in dgmatv.f */
/* ---  arguments variables ... */

/* ---  Executable statements ... */
/* ---  load a Harwell-Boeing matrix ... */
    rmat_1.n = 5000;
    rmat_1.nz = 600000;
    loadhb_("../../../data/orani678$", "coo", &rmat_1.n, &rmat_1.nz, rmat_1.ia, 
	    rmat_1.ja, rmat_1.a, iwsp, (ftnlen)17, (ftnlen)3);
/* ---  compute the infinite norm of A ... */
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wsp[i__ - 1] = 0.;
    }
    i__1 = rmat_1.nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wsp[rmat_1.ia[i__ - 1] - 1] += (d__1 = rmat_1.a[i__ - 1], abs(d__1));
    }
    anorm = wsp[0];
    i__1 = rmat_1.n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (anorm < wsp[i__ - 1]) {
	    anorm = wsp[i__ - 1];
	}
    }
    s_wsfe(&io___5);
    do_fio(&c__1, "||A||_inf= ", (ftnlen)11);
    do_fio(&c__1, (char *)&anorm, (ftnlen)sizeof(doublereal));
    e_wsfe();
/* ---  back to CCS format ... */
    dgcnvr_("coo", "ccs", "n", &rmat_1.n, &rmat_1.n, &rmat_1.nz, rmat_1.ia, 
	    rmat_1.ja, rmat_1.a, iwsp, (ftnlen)3, (ftnlen)3, (ftnlen)1);
/* ---  set other input arguments ... */
    t = 10.;
    tol = 0.;
    m = 30;
    itrace = 0;
/* ---- First Run: ------------------------------------------------------ */
/* ********************************************************************** */
/* ---  the operand vector u is set to zero ... */
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	u[i__ - 1] = 0.;
    }
/* ---  the operand vector v is set to (1, ..., 1)' ... */
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[i__ - 1] = 1.;
    }
/* ---  compute w = exp(t*A)*v with DGPHIV ... */
    tic = clock_();
    dgphiv_(&rmat_1.n, &m, &t, u, v, w, &tol, &anorm, wsp, &c_b10, iwsp, &
	    c__5000, (S_fp)dgccsv_, &itrace, &iflag);
    tac = clock_();
/* --- */
    s_wsfe(&io___16);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___17);
    do_fio(&c__1, "DGPHIV (CCS) has completed:", (ftnlen)27);
    e_wsfe();
    s_wsfe(&io___18);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___19);
    do_fio(&c__1, "w(1:10) =", (ftnlen)9);
    e_wsfe();
    for (i__ = 1; i__ <= 10; ++i__) {
	s_wsle(&io___20);
	do_lio(&c__5, &c__1, (char *)&w[i__ - 1], (ftnlen)sizeof(doublereal));
	e_wsle();
    }
    dcopy_(&rmat_1.n, w, &c__1, wsave, &c__1);
/* ---  display some statistics if desired ... */
    s_wsfe(&io___22);
    do_fio(&c__1, "final report----------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___23);
    do_fio(&c__1, "runtime   = ", (ftnlen)12);
    d__1 = tac - tic;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___24);
    do_fio(&c__1, "||A||_inf = ", (ftnlen)12);
    do_fio(&c__1, (char *)&anorm, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___25);
    do_fio(&c__1, "nz        =", (ftnlen)11);
    do_fio(&c__1, (char *)&rmat_1.nz, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___26);
    do_fio(&c__1, "n         =", (ftnlen)11);
    do_fio(&c__1, (char *)&rmat_1.n, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___27);
    do_fio(&c__1, "m         =", (ftnlen)11);
    do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___28);
    do_fio(&c__1, "itrace    =", (ftnlen)11);
    do_fio(&c__1, (char *)&itrace, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___29);
    do_fio(&c__1, "iflag     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iflag, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___30);
    do_fio(&c__1, "ibrkflag  =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[5], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___31);
    do_fio(&c__1, "mbrkdwn   =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[6], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___32);
    do_fio(&c__1, "nstep     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[3], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___33);
    do_fio(&c__1, "nreject   =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[4], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___34);
    do_fio(&c__1, "nmult     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[0], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___35);
    do_fio(&c__1, "nexph     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[1], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___36);
    do_fio(&c__1, "nscale    =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[2], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___37);
    do_fio(&c__1, "tol       = ", (ftnlen)12);
    do_fio(&c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___38);
    do_fio(&c__1, "t         = ", (ftnlen)12);
    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___39);
    do_fio(&c__1, "tbrkdwn   = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[6], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___40);
    do_fio(&c__1, "step_min  = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[0], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___41);
    do_fio(&c__1, "step_max  = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[1], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___42);
    do_fio(&c__1, "max_round = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[2], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___43);
    do_fio(&c__1, "sum_round = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[3], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___44);
    do_fio(&c__1, "max_error = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[4], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___45);
    do_fio(&c__1, "sum_error = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[5], (ftnlen)sizeof(doublereal));
    e_wsfe();
/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* ---- Second Run ------------------------------------------------------ */
/* ********************************************************************** */
/* ---  the operand vector u is set to (1, ..., 1)' ... */
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	u[i__ - 1] = 1.;
    }
/* ---  the operand vector v is set to zero ... */
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[i__ - 1] = 0.;
    }
/* ---  compute w = t*phi(t*A)*u with DGPHIV ... */
    tic = clock_();
    dgphiv_(&rmat_1.n, &m, &t, u, v, w, &tol, &anorm, wsp, &c_b10, iwsp, &
	    c__5000, (S_fp)dgccsv_, &itrace, &iflag);
    tac = clock_();
/* --- */
    s_wsfe(&io___46);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___47);
    do_fio(&c__1, "DGPHIV (CCS) has completed:", (ftnlen)27);
    e_wsfe();
    s_wsfe(&io___48);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___49);
    do_fio(&c__1, "w(1:10) =", (ftnlen)9);
    e_wsfe();
    for (i__ = 1; i__ <= 10; ++i__) {
	s_wsle(&io___50);
	do_lio(&c__5, &c__1, (char *)&w[i__ - 1], (ftnlen)sizeof(doublereal));
	e_wsle();
    }
/* ---  display some statistics if desired ... */
    s_wsfe(&io___51);
    do_fio(&c__1, "final report----------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___52);
    do_fio(&c__1, "runtime   = ", (ftnlen)12);
    d__1 = tac - tic;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___53);
    do_fio(&c__1, "||A||_inf = ", (ftnlen)12);
    do_fio(&c__1, (char *)&anorm, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___54);
    do_fio(&c__1, "nz        =", (ftnlen)11);
    do_fio(&c__1, (char *)&rmat_1.nz, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___55);
    do_fio(&c__1, "n         =", (ftnlen)11);
    do_fio(&c__1, (char *)&rmat_1.n, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___56);
    do_fio(&c__1, "m         =", (ftnlen)11);
    do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___57);
    do_fio(&c__1, "itrace    =", (ftnlen)11);
    do_fio(&c__1, (char *)&itrace, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___58);
    do_fio(&c__1, "iflag     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iflag, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___59);
    do_fio(&c__1, "ibrkflag  =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[5], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___60);
    do_fio(&c__1, "mbrkdwn   =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[6], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___61);
    do_fio(&c__1, "nstep     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[3], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___62);
    do_fio(&c__1, "nreject   =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[4], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___63);
    do_fio(&c__1, "nmult     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[0], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___64);
    do_fio(&c__1, "nexph     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[1], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___65);
    do_fio(&c__1, "nscale    =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[2], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___66);
    do_fio(&c__1, "tol       = ", (ftnlen)12);
    do_fio(&c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___67);
    do_fio(&c__1, "t         = ", (ftnlen)12);
    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___68);
    do_fio(&c__1, "tbrkdwn   = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[6], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___69);
    do_fio(&c__1, "step_min  = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[0], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___70);
    do_fio(&c__1, "step_max  = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[1], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___71);
    do_fio(&c__1, "max_round = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[2], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___72);
    do_fio(&c__1, "sum_round = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[3], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___73);
    do_fio(&c__1, "max_error = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[4], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___74);
    do_fio(&c__1, "sum_error = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[5], (ftnlen)sizeof(doublereal));
    e_wsfe();
/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* ---  now check ||(u + A * w) - wsave||/||wsave||. Due to the specific */
/*     previous settings, the answer should agree with tol. */

    dgccsv_(w, v);
    daxpy_(&rmat_1.n, &c_b184, u, &c__1, v, &c__1);
    daxpy_(&rmat_1.n, &c_b187, wsave, &c__1, v, &c__1);
    tmp = dnrm2_(&rmat_1.n, v, &c__1) / dnrm2_(&rmat_1.n, wsave, &c__1);
    s_wsle(&io___76);
    e_wsle();
    s_wsle(&io___77);
    do_lio(&c__9, &c__1, "relative difference (phi vs. exp) =", (ftnlen)35);
    do_lio(&c__5, &c__1, (char *)&tmp, (ftnlen)sizeof(doublereal));
    e_wsle();
    s_wsle(&io___78);
    e_wsle();
/* ---- Third Run ------------------------------------------------------- */
/* ********************************************************************** */
/* ---  the operand vector v is set to (1, ..., 1)' ... */
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[i__ - 1] = 1.;
    }
/* ---  the operand vector u is set to (1, ..., 1)' ... */
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[i__ - 1] = 1.;
    }
/* ---  compute w = exp(t*A)*v + t*phi(t*A)*u with DGPHIV ... */
    tic = clock_();
    dgphiv_(&rmat_1.n, &m, &t, u, v, w, &tol, &anorm, wsp, &c_b10, iwsp, &
	    c__5000, (S_fp)dgccsv_, &itrace, &iflag);
    tac = clock_();
/* --- */
    s_wsfe(&io___79);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___80);
    do_fio(&c__1, "DGPHIV (CCS) has completed:", (ftnlen)27);
    e_wsfe();
    s_wsfe(&io___81);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___82);
    do_fio(&c__1, "w(1:10) =", (ftnlen)9);
    e_wsfe();
    for (i__ = 1; i__ <= 10; ++i__) {
	s_wsle(&io___83);
	do_lio(&c__5, &c__1, (char *)&w[i__ - 1], (ftnlen)sizeof(doublereal));
	e_wsle();
    }
/* ---  display some statistics if desired ... */
    s_wsfe(&io___84);
    do_fio(&c__1, "final report----------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___85);
    do_fio(&c__1, "runtime   = ", (ftnlen)12);
    d__1 = tac - tic;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___86);
    do_fio(&c__1, "||A||_inf = ", (ftnlen)12);
    do_fio(&c__1, (char *)&anorm, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___87);
    do_fio(&c__1, "nz        =", (ftnlen)11);
    do_fio(&c__1, (char *)&rmat_1.nz, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___88);
    do_fio(&c__1, "n         =", (ftnlen)11);
    do_fio(&c__1, (char *)&rmat_1.n, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___89);
    do_fio(&c__1, "m         =", (ftnlen)11);
    do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___90);
    do_fio(&c__1, "itrace    =", (ftnlen)11);
    do_fio(&c__1, (char *)&itrace, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___91);
    do_fio(&c__1, "iflag     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iflag, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___92);
    do_fio(&c__1, "ibrkflag  =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[5], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___93);
    do_fio(&c__1, "mbrkdwn   =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[6], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___94);
    do_fio(&c__1, "nstep     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[3], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___95);
    do_fio(&c__1, "nreject   =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[4], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___96);
    do_fio(&c__1, "nmult     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[0], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___97);
    do_fio(&c__1, "nexph     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[1], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___98);
    do_fio(&c__1, "nscale    =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[2], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___99);
    do_fio(&c__1, "tol       = ", (ftnlen)12);
    do_fio(&c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___100);
    do_fio(&c__1, "t         = ", (ftnlen)12);
    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___101);
    do_fio(&c__1, "tbrkdwn   = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[6], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___102);
    do_fio(&c__1, "step_min  = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[0], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___103);
    do_fio(&c__1, "step_max  = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[1], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___104);
    do_fio(&c__1, "max_round = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[2], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___105);
    do_fio(&c__1, "sum_round = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[3], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___106);
    do_fio(&c__1, "max_error = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[4], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___107);
    do_fio(&c__1, "sum_error = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[5], (ftnlen)sizeof(doublereal));
    e_wsfe();

    return 0;
} /* MAIN__ */


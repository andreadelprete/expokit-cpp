/* sample_b.f -- translated by f2c (version 20100827).
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
static integer c_b7 = 273527;
static integer c__5000 = 5000;
static integer c__5 = 5;


/* ---  sample program illustrating the use of DSEXPV ... */
/*     Forward-backward problem (Example 6.4 in the Expokit report) ... */

/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_9001[] = "(a)";
    static char fmt_9002[] = "(a,e8.2)";
    static char fmt_9003[] = "(a,i9)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer i__, j, m;
    static doublereal t, v[5000], w[5000], tac, tic, tol, tmp, wsp[273527];
    static integer iwsp[5000], iflag;
    extern doublereal clock_(void);
    static doublereal anorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), loadhb_(char *, char *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen);
    static integer itrace;
    extern /* Subroutine */ int dgcrsv_();
    extern /* Subroutine */ int dsexpv_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, U_fp, integer *, 
	    integer *);

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 6, 0, "(A,E8.2)", 0 };
    static cilist io___17 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___18 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___19 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___20 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
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
    static cilist io___46 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___47 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___49 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___50 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };
    static cilist io___53 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___54 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___55 = { 0, 6, 0, fmt_9002, 0 };
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
    static cilist io___66 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___67 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___68 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___69 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___70 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___71 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___72 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___73 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___74 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___75 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___76 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___77 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___78 = { 0, 6, 0, fmt_9002, 0 };


/* ---  matrix data ... */
/* ---- BEWARE: these values must match those in dgmatv.f */
/* ---  arguments variables ... */

/* ---  Executable statements ... */
/* ---  load the matrix ... */
    rmat_1.n = 5000;
    rmat_1.nz = 600000;
    loadhb_("../data/gr3030$", "crs", &rmat_1.n, &rmat_1.nz, rmat_1.ia, 
	    rmat_1.ja, rmat_1.a, iwsp, (ftnlen)15, (ftnlen)3);
/* ---  compute the infinite norm of A ... */
    anorm = 0.;
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tmp = 0.;
	i__2 = rmat_1.ia[i__] - 1;
	for (j = rmat_1.ia[i__ - 1]; j <= i__2; ++j) {
	    tmp += (d__1 = rmat_1.a[j - 1], abs(d__1));
	}
	if (anorm < tmp) {
	    anorm = tmp;
	}
    }
    s_wsfe(&io___6);
    do_fio(&c__1, "||A||_inf= ", (ftnlen)11);
    do_fio(&c__1, (char *)&anorm, (ftnlen)sizeof(doublereal));
    e_wsfe();
/* ---  the operand vector v is set to (1, ..., 1)' ... */
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[i__ - 1] = 1.;
    }
/* ---  set other input arguments ... */
    t = 1.;
    tol = 1e-10;
    m = 30;
    itrace = 0;
/* ---  compute exp(t*A)v with CRS format ... */
    tic = clock_();
    dsexpv_(&rmat_1.n, &m, &t, v, w, &tol, &anorm, wsp, &c_b7, iwsp, &c__5000,
	     (U_fp)dgcrsv_, &itrace, &iflag);
    tac = clock_();
/* --- */
    s_wsfe(&io___17);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___18);
    do_fio(&c__1, "DSEXPV (Forward) has completed:", (ftnlen)31);
    e_wsfe();
    s_wsfe(&io___19);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___20);
    do_fio(&c__1, "w(1:10) =", (ftnlen)9);
    e_wsfe();
    for (i__ = 1; i__ <= 10; ++i__) {
	s_wsle(&io___21);
	do_lio(&c__5, &c__1, (char *)&w[i__ - 1], (ftnlen)sizeof(doublereal));
	e_wsle();
    }
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
    s_wsfe(&io___46);
    do_fio(&c__1, "hump      = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[8], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___47);
    do_fio(&c__1, "scale-norm= ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[9], (ftnlen)sizeof(doublereal));
    e_wsfe();
/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* ---  compute exp(-t*A)w with CRS format ... */
    t = -t;
    dcopy_(&rmat_1.n, w, &c__1, v, &c__1);
    tic = clock_();
    dsexpv_(&rmat_1.n, &m, &t, v, w, &tol, &anorm, wsp, &c_b7, iwsp, &c__5000,
	     (U_fp)dgcrsv_, &itrace, &iflag);
    tac = clock_();
/* --- */
    s_wsfe(&io___48);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___49);
    do_fio(&c__1, "DSEXPV (Backward) has completed:", (ftnlen)32);
    e_wsfe();
    s_wsfe(&io___50);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___51);
    do_fio(&c__1, "w(1:10) =", (ftnlen)9);
    e_wsfe();
    for (i__ = 1; i__ <= 10; ++i__) {
	s_wsle(&io___52);
	do_lio(&c__5, &c__1, (char *)&w[i__ - 1], (ftnlen)sizeof(doublereal));
	e_wsle();
    }
/* ---  display some statistics if desired ... */
    s_wsfe(&io___53);
    do_fio(&c__1, "final report----------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___54);
    do_fio(&c__1, "runtime   = ", (ftnlen)12);
    d__1 = tac - tic;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___55);
    do_fio(&c__1, "||A||_inf = ", (ftnlen)12);
    do_fio(&c__1, (char *)&anorm, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___56);
    do_fio(&c__1, "nz        =", (ftnlen)11);
    do_fio(&c__1, (char *)&rmat_1.nz, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___57);
    do_fio(&c__1, "n         =", (ftnlen)11);
    do_fio(&c__1, (char *)&rmat_1.n, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___58);
    do_fio(&c__1, "m         =", (ftnlen)11);
    do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___59);
    do_fio(&c__1, "itrace    =", (ftnlen)11);
    do_fio(&c__1, (char *)&itrace, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___60);
    do_fio(&c__1, "iflag     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iflag, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___61);
    do_fio(&c__1, "ibrkflag  =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[5], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___62);
    do_fio(&c__1, "mbrkdwn   =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[6], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___63);
    do_fio(&c__1, "nstep     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[3], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___64);
    do_fio(&c__1, "nreject   =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[4], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___65);
    do_fio(&c__1, "nmult     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[0], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___66);
    do_fio(&c__1, "nexph     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[1], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___67);
    do_fio(&c__1, "nscale    =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[2], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___68);
    do_fio(&c__1, "tol       = ", (ftnlen)12);
    do_fio(&c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___69);
    do_fio(&c__1, "t         = ", (ftnlen)12);
    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___70);
    do_fio(&c__1, "tbrkdwn   = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[6], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___71);
    do_fio(&c__1, "step_min  = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[0], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___72);
    do_fio(&c__1, "step_max  = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[1], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___73);
    do_fio(&c__1, "max_round = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[2], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___74);
    do_fio(&c__1, "sum_round = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[3], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___75);
    do_fio(&c__1, "max_error = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[4], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___76);
    do_fio(&c__1, "sum_error = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[5], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___77);
    do_fio(&c__1, "hump      = ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[8], (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___78);
    do_fio(&c__1, "scale-norm= ", (ftnlen)12);
    do_fio(&c__1, (char *)&wsp[9], (ftnlen)sizeof(doublereal));
    e_wsfe();
    return 0;
} /* MAIN__ */


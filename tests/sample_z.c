/* sample_z.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static integer c__1 = 1;
static integer c_b8 = 299527;
static integer c__52 = 52;
static integer c__7 = 7;
static integer c__9 = 9;
static integer c__3 = 3;


/* ---  sample program illustrating the use of ZGEXPV and ZHEXPV ... */
/*     Hermitian problem (Example 6.2 in the Expokit report) ... */

/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_9001[] = "(a)";
    static char fmt_9002[] = "(a,e8.2)";
    static char fmt_9003[] = "(a,i9)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer i__, j, m;
    static doublereal t;
    static doublecomplex v[5500], w[5500];
    static doublereal s1, s2, tac, tic, tol;
    static integer nnz;
    static doublecomplex wsp[299527];
    static integer iwsp[52], iflag, iseed[4];
    extern doublereal clock_(void);
    static doublereal anorm;
    extern doublereal dlaran_(integer *);
    static integer itrace;
    extern /* Subroutine */ int getpat_(char *, integer *, integer *, integer 
	    *, integer *, ftnlen);
    extern /* Subroutine */ int zgcoov_();
    extern /* Subroutine */ int zgexpv_(integer *, integer *, doublereal *, 
	    doublecomplex *, doublecomplex *, doublereal *, doublereal *, 
	    doublecomplex *, integer *, integer *, integer *, U_fp, integer *,
	     integer *), zhexpv_(integer *, integer *, doublereal *, 
	    doublecomplex *, doublecomplex *, doublereal *, doublereal *, 
	    doublecomplex *, integer *, integer *, integer *, U_fp, integer *,
	     integer *);

    /* Fortran I/O blocks */
    static cilist io___9 = { 0, 6, 0, "(A,E8.2)", 0 };
    static cilist io___20 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___21 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___22 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___23 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___26 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___27 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___35 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___36 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___37 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___40 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___41 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___42 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___43 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___44 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___46 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___47 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___49 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___50 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___52 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___53 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___54 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___55 = { 0, 6, 0, 0, 0 };
    static cilist io___56 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___57 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___58 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___59 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___60 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___61 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___62 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___63 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___64 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___65 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___66 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___67 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___68 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___69 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___70 = { 0, 6, 0, fmt_9003, 0 };
    static cilist io___71 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___72 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___73 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___74 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___75 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___76 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___77 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___78 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___79 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___80 = { 0, 6, 0, fmt_9002, 0 };
    static cilist io___81 = { 0, 6, 0, fmt_9002, 0 };


/* ---  matrix data ... */
/* ---  BEWARE: these values must match those in zgmatv.f */
/* ---  arguments variables ... */

/* ---  Executable statements ... */
/* ---  load a symmetric pattern ... */
    cmat_1.n = 5500;
    cmat_1.nz = 25000;
    getpat_("../data/bcspwr10$", &cmat_1.n, &cmat_1.nz, cmat_1.ia, cmat_1.ja, 
	    (ftnlen)17);
/* ---  for the purpose of the experiments, expand to COOrdinates ... */
    nnz = cmat_1.nz;
    for (j = cmat_1.n; j >= 1; --j) {
	i__1 = cmat_1.ja[j] - cmat_1.ja[j - 1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    cmat_1.ja[nnz - 1] = j;
	    --nnz;
	}
    }
/* ---  fill-in an Hermitian matrix -- the conjugate part is included */
    iseed[0] = 0;
    iseed[1] = 0;
    iseed[2] = 0;
    iseed[3] = 5;
    nnz = cmat_1.nz;
    i__1 = cmat_1.nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (cmat_1.ia[i__ - 1] != cmat_1.ja[i__ - 1]) {
	    s1 = dlaran_(iseed) * 10. - 5.;
	    s2 = dlaran_(iseed) * 10. - 5.;
	    i__2 = i__ - 1;
	    z__1.r = s1, z__1.i = s2;
	    cmat_1.a[i__2].r = z__1.r, cmat_1.a[i__2].i = z__1.i;
	    ++nnz;
	    i__2 = nnz - 1;
	    d_cnjg(&z__1, &cmat_1.a[i__ - 1]);
	    cmat_1.a[i__2].r = z__1.r, cmat_1.a[i__2].i = z__1.i;
	    cmat_1.ia[nnz - 1] = cmat_1.ja[i__ - 1];
	    cmat_1.ja[nnz - 1] = cmat_1.ia[i__ - 1];
	} else {
	    s1 = dlaran_(iseed) * 10. - 5.;
	    i__2 = i__ - 1;
	    z__1.r = s1, z__1.i = 0.;
	    cmat_1.a[i__2].r = z__1.r, cmat_1.a[i__2].i = z__1.i;
	}
    }
    cmat_1.nz = nnz;
/* ---  compute the infinite norm of A ... */
    i__1 = cmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	wsp[i__2].r = 0., wsp[i__2].i = 0.;
    }
    i__1 = cmat_1.nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = cmat_1.ia[i__ - 1] - 1;
	i__3 = cmat_1.ia[i__ - 1] - 1;
	d__1 = z_abs(&cmat_1.a[i__ - 1]);
	z__1.r = wsp[i__3].r + d__1, z__1.i = wsp[i__3].i;
	wsp[i__2].r = z__1.r, wsp[i__2].i = z__1.i;
    }
    anorm = wsp[0].r;
    i__1 = cmat_1.n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	if (anorm < wsp[i__2].r) {
	    i__3 = i__ - 1;
	    anorm = wsp[i__3].r;
	}
    }
/* *--   convert from COO to CRS ... */
/*      call zgcnvr( 'coo','crs','n', n,n, nz, ia, ja, a, iwsp ) */

/* *---  compute the infinite norm of A ... */
/*      anorm = 0.0d0 */
/*      do i = 1,n */
/*         s1 = 0.0d0 */
/*         do j = ia(i),ia(i+1)-1 */
/*            s1 = s1 + ABS( a(j) ) */
/*         enddo */
/*         if ( anorm.lt.tmp ) anorm = s1 */
/*      enddo */
    s_wsfe(&io___9);
    do_fio(&c__1, "||A||_inf= ", (ftnlen)11);
    do_fio(&c__1, (char *)&anorm, (ftnlen)sizeof(doublereal));
    e_wsfe();
/* ---  the operand vector v is set to e_1 + e_n ... */
    i__1 = cmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	v[i__2].r = 0., v[i__2].i = 0.;
    }
    v[0].r = 1., v[0].i = 0.;
    i__1 = cmat_1.n - 1;
    v[i__1].r = 1., v[i__1].i = 0.;
/* ---  set other input arguments ... */
    t = 1.;
    tol = 1e-5;
    m = 30;
    itrace = 0;
/* ---  compute w = exp(t*A)v with ZGEXPV ... */
    tic = clock_();
    zgexpv_(&cmat_1.n, &m, &t, v, w, &tol, &anorm, wsp, &c_b8, iwsp, &c__52, (
	    U_fp)zgcoov_, &itrace, &iflag);
    tac = clock_();
/* --- */
    s_wsfe(&io___20);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___21);
    do_fio(&c__1, "ZGEXPV has completed:", (ftnlen)21);
    e_wsfe();
    s_wsfe(&io___22);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___23);
    do_fio(&c__1, "w(1:10) =", (ftnlen)9);
    e_wsfe();
    for (i__ = 1; i__ <= 10; ++i__) {
	s_wsle(&io___24);
	do_lio(&c__7, &c__1, (char *)&w[i__ - 1], (ftnlen)sizeof(
		doublecomplex));
	e_wsle();
    }
/* ---  display some statistics if desired ... */
    s_wsfe(&io___25);
    do_fio(&c__1, "final report----------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___26);
    do_fio(&c__1, "runtime   = ", (ftnlen)12);
    d__1 = tac - tic;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___27);
    do_fio(&c__1, "||A||_inf = ", (ftnlen)12);
    do_fio(&c__1, (char *)&anorm, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___28);
    do_fio(&c__1, "nz        =", (ftnlen)11);
    do_fio(&c__1, (char *)&cmat_1.nz, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___29);
    do_fio(&c__1, "n         =", (ftnlen)11);
    do_fio(&c__1, (char *)&cmat_1.n, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___30);
    do_fio(&c__1, "m         =", (ftnlen)11);
    do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___31);
    do_fio(&c__1, "itrace    =", (ftnlen)11);
    do_fio(&c__1, (char *)&itrace, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___32);
    do_fio(&c__1, "iflag     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iflag, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___33);
    do_fio(&c__1, "ibrkflag  =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[5], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___34);
    do_fio(&c__1, "mbrkdwn   =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[6], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___35);
    do_fio(&c__1, "nstep     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[3], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___36);
    do_fio(&c__1, "nreject   =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[4], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___37);
    do_fio(&c__1, "nmult     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[0], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___38);
    do_fio(&c__1, "nexph     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[1], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___39);
    do_fio(&c__1, "nscale    =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[2], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___40);
    do_fio(&c__1, "tol       = ", (ftnlen)12);
    do_fio(&c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___41);
    do_fio(&c__1, "t         = ", (ftnlen)12);
    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___42);
    do_fio(&c__1, "tbrkdwn   = ", (ftnlen)12);
    d__1 = wsp[6].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___43);
    do_fio(&c__1, "step_min  = ", (ftnlen)12);
    d__1 = wsp[0].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___44);
    do_fio(&c__1, "step_max  = ", (ftnlen)12);
    d__1 = wsp[1].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___45);
    do_fio(&c__1, "max_round = ", (ftnlen)12);
    d__1 = wsp[2].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___46);
    do_fio(&c__1, "sum_round = ", (ftnlen)12);
    d__1 = wsp[3].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___47);
    do_fio(&c__1, "max_error = ", (ftnlen)12);
    d__1 = wsp[4].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___48);
    do_fio(&c__1, "sum_error = ", (ftnlen)12);
    d__1 = wsp[5].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___49);
    do_fio(&c__1, "hump      = ", (ftnlen)12);
    d__1 = wsp[8].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___50);
    do_fio(&c__1, "scale-norm= ", (ftnlen)12);
    d__1 = wsp[9].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* ---  compute w = exp(t*A)v with ZHEXPV ... */
    tic = clock_();
    zhexpv_(&cmat_1.n, &m, &t, v, w, &tol, &anorm, wsp, &c_b8, iwsp, &c__52, (
	    U_fp)zgcoov_, &itrace, &iflag);
    tac = clock_();
/* --- */
    s_wsfe(&io___51);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___52);
    do_fio(&c__1, "ZHEXPV has completed:", (ftnlen)21);
    e_wsfe();
    s_wsfe(&io___53);
    do_fio(&c__1, "----------------------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___54);
    do_fio(&c__1, "w(1:10) =", (ftnlen)9);
    e_wsfe();
    for (i__ = 1; i__ <= 10; ++i__) {
	s_wsle(&io___55);
	do_lio(&c__7, &c__1, (char *)&w[i__ - 1], (ftnlen)sizeof(
		doublecomplex));
	e_wsle();
    }
/* ---  display some statistics if desired ... */
    s_wsfe(&io___56);
    do_fio(&c__1, "final report----------------------------------------", (
	    ftnlen)52);
    e_wsfe();
    s_wsfe(&io___57);
    do_fio(&c__1, "runtime   = ", (ftnlen)12);
    d__1 = tac - tic;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___58);
    do_fio(&c__1, "||A||_inf = ", (ftnlen)12);
    do_fio(&c__1, (char *)&anorm, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___59);
    do_fio(&c__1, "nz        =", (ftnlen)11);
    do_fio(&c__1, (char *)&cmat_1.nz, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___60);
    do_fio(&c__1, "n         =", (ftnlen)11);
    do_fio(&c__1, (char *)&cmat_1.n, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___61);
    do_fio(&c__1, "m         =", (ftnlen)11);
    do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___62);
    do_fio(&c__1, "itrace    =", (ftnlen)11);
    do_fio(&c__1, (char *)&itrace, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___63);
    do_fio(&c__1, "iflag     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iflag, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___64);
    do_fio(&c__1, "ibrkflag  =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[5], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___65);
    do_fio(&c__1, "mbrkdwn   =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[6], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___66);
    do_fio(&c__1, "nstep     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[3], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___67);
    do_fio(&c__1, "nreject   =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[4], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___68);
    do_fio(&c__1, "nmult     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[0], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___69);
    do_fio(&c__1, "nexph     =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[1], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___70);
    do_fio(&c__1, "nscale    =", (ftnlen)11);
    do_fio(&c__1, (char *)&iwsp[2], (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___71);
    do_fio(&c__1, "tol       = ", (ftnlen)12);
    do_fio(&c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___72);
    do_fio(&c__1, "t         = ", (ftnlen)12);
    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___73);
    do_fio(&c__1, "tbrkdwn   = ", (ftnlen)12);
    d__1 = wsp[6].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___74);
    do_fio(&c__1, "step_min  = ", (ftnlen)12);
    d__1 = wsp[0].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___75);
    do_fio(&c__1, "step_max  = ", (ftnlen)12);
    d__1 = wsp[1].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___76);
    do_fio(&c__1, "max_round = ", (ftnlen)12);
    d__1 = wsp[2].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___77);
    do_fio(&c__1, "sum_round = ", (ftnlen)12);
    d__1 = wsp[3].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___78);
    do_fio(&c__1, "max_error = ", (ftnlen)12);
    d__1 = wsp[4].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___79);
    do_fio(&c__1, "sum_error = ", (ftnlen)12);
    d__1 = wsp[5].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___80);
    do_fio(&c__1, "hump      = ", (ftnlen)12);
    d__1 = wsp[8].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___81);
    do_fio(&c__1, "scale-norm= ", (ftnlen)12);
    d__1 = wsp[9].r;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();
    return 0;
} /* MAIN__ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* Subroutine */ int getpat_(char *filename, integer *n, integer *nz, integer 
	*ia, integer *ja, ftnlen filename_len)
{
    /* Format strings */
    static char fmt_10[] = "(a72,a8/5i14/a3,11x,4i14/2a16,2a20)";
    static char fmt_11[] = "(a1,13x,2i14)";

    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer i_indx(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer f_open(olist *), s_rsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_rsfe(void), s_wsle(cilist *), do_lio(integer *, 
	    integer *, char *, ftnlen), e_wsle(void), f_clos(cllist *);

    /* Local variables */
    static integer i__, io, nn;
    static char key[8];
    static integer nnz, nrhs;
    static char type__[3];
    static integer nrow;
    static char title[72];
    static integer indcrd, valcrd;
    static char indfmt[16];
    static integer rhscrd;
    static char valfmt[20];
    static integer ptrcrd, totcrd;
    static char rhsfmt[20];
    static integer nrhsix;
    static char ptrfmt[16], rhstyp[1];

    /* Fortran I/O blocks */
    static cilist io___84 = { 0, 7, 0, fmt_10, 0 };
    static cilist io___101 = { 0, 6, 0, 0, 0 };
    static cilist io___102 = { 0, 6, 0, 0, 0 };
    static cilist io___103 = { 0, 7, 0, fmt_11, 0 };
    static cilist io___106 = { 0, 6, 0, 0, 0 };
    static cilist io___107 = { 0, 7, 0, ptrfmt, 0 };
    static cilist io___108 = { 0, 7, 0, indfmt, 0 };
    static cilist io___109 = { 0, 6, 0, 0, 0 };


/* ---  load a Harwell-Boeing pattern ... */
/* ---  arguments are fully described in LOADHB */
/* --- */
/* --- */
    /* Parameter adjustments */
    --ja;
    --ia;

    /* Function Body */
    i__ = i_indx(filename, "$", (ftnlen)80, (ftnlen)1) - 1;
    if (i__ <= 0) {
	s_stop("in GETPAT. Bad filename", (ftnlen)23);
    }
    o__1.oerr = 1;
    o__1.ounit = 7;
    o__1.ofnmlen = i__;
    o__1.ofnm = filename;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    io = f_open(&o__1);
    if (io != 0) {
	s_stop("Could not access Harwell-Boeing matrix", (ftnlen)38);
    }
    s_rsfe(&io___84);
    do_fio(&c__1, title, (ftnlen)72);
    do_fio(&c__1, key, (ftnlen)8);
    do_fio(&c__1, (char *)&totcrd, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ptrcrd, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&indcrd, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&valcrd, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&rhscrd, (ftnlen)sizeof(integer));
    do_fio(&c__1, type__, (ftnlen)3);
    do_fio(&c__1, (char *)&nrow, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nn, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nnz, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nrhs, (ftnlen)sizeof(integer));
    do_fio(&c__1, ptrfmt, (ftnlen)16);
    do_fio(&c__1, indfmt, (ftnlen)16);
    do_fio(&c__1, valfmt, (ftnlen)20);
    do_fio(&c__1, rhsfmt, (ftnlen)20);
    e_rsfe();
    s_wsle(&io___101);
    do_lio(&c__9, &c__1, title, (ftnlen)72);
    do_lio(&c__9, &c__1, "type :", (ftnlen)6);
    do_lio(&c__9, &c__1, type__, (ftnlen)3);
    do_lio(&c__9, &c__1, " size :", (ftnlen)7);
    do_lio(&c__3, &c__1, (char *)&nrow, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nn, (ftnlen)sizeof(integer));
    e_wsle();
    s_wsle(&io___102);
    do_lio(&c__9, &c__1, "order :", (ftnlen)7);
    do_lio(&c__3, &c__1, (char *)&nn, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " number of nonzero :", (ftnlen)20);
    do_lio(&c__3, &c__1, (char *)&nnz, (ftnlen)sizeof(integer));
    e_wsle();
    if (nn > *n) {
	s_stop("Please increase nmax.", (ftnlen)21);
    }
    if (nnz > *nz) {
	s_stop("Please increase nzmax.", (ftnlen)22);
    }
    *n = nn;
    *nz = nnz;
    if (rhscrd > 0) {
	s_rsfe(&io___103);
	do_fio(&c__1, rhstyp, (ftnlen)1);
	do_fio(&c__1, (char *)&nrhs, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nrhsix, (ftnlen)sizeof(integer));
	e_rsfe();
	s_wsle(&io___106);
	do_lio(&c__9, &c__1, "There is a second hand", (ftnlen)22);
	e_wsle();
    }
/* ---  read data... */
    s_rsfe(&io___107);
    i__1 = *n + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ja[i__], (ftnlen)sizeof(integer));
    }
    e_rsfe();
    s_rsfe(&io___108);
    i__1 = *nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ia[i__], (ftnlen)sizeof(integer));
    }
    e_rsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 7;
    cl__1.csta = 0;
    f_clos(&cl__1);
    s_wsle(&io___109);
    do_lio(&c__9, &c__1, "Harwell-Boeing pattern loaded", (ftnlen)29);
    e_wsle();
    return 0;
} /* getpat_ */

/* ----------------------------------------------------------------------| */
doublereal dlaran_(integer *iseed)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer it1, it2, it3, it4;


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLARAN returns a random real number from a uniform (0,1) */
/*  distribution. */

/*  Arguments */
/*  ========= */

/*  ISEED   (input/output) INTEGER array, dimension (4) */
/*          On entry, the seed of the random number generator; the array */
/*          elements must be between 0 and 4095, and ISEED(4) must be */
/*          odd. */
/*          On exit, the seed is updated. */

/*  Further Details */
/*  =============== */

/*  This routine uses a multiplicative congruential method with modulus */
/*  2**48 and multiplier 33952834046453 (see G.S.Fishman, */
/*  'Multiplicative congruential random number generators with modulus */
/*  2**b: an exhaustive analysis for b = 32 and a partial analysis for */
/*  b = 48', Math. Comp. 189, pp 331-344, 1990). */

/*  48-bit integers are stored in 4 integer array elements with 12 bits */
/*  per element. Hence the routine is portable across machines with */
/*  integers of 32 bits or more. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     multiply the seed by the multiplier modulo 2**48 */

    /* Parameter adjustments */
    --iseed;

    /* Function Body */
    it4 = iseed[4] * 2549;
    it3 = it4 / 4096;
    it4 -= it3 << 12;
    it3 = it3 + iseed[3] * 2549 + iseed[4] * 2508;
    it2 = it3 / 4096;
    it3 -= it2 << 12;
    it2 = it2 + iseed[2] * 2549 + iseed[3] * 2508 + iseed[4] * 322;
    it1 = it2 / 4096;
    it2 -= it1 << 12;
    it1 = it1 + iseed[1] * 2549 + iseed[2] * 2508 + iseed[3] * 322 + iseed[4] 
	    * 494;
    it1 %= 4096;

/*     return updated seed */

    iseed[1] = it1;
    iseed[2] = it2;
    iseed[3] = it3;
    iseed[4] = it4;

/*     convert 48-bit integer to a real number in the interval (0,1) */

    ret_val = ((doublereal) it1 + ((doublereal) it2 + ((doublereal) it3 + (
	    doublereal) it4 * 2.44140625e-4) * 2.44140625e-4) * 2.44140625e-4)
	     * 2.44140625e-4;
    return ret_val;

/*     End of DLARAN */

} /* dlaran_ */


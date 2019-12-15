/* sample_d.f -- translated by f2c (version 20100827).
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

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__6 = 6;
static integer c__50 = 50;
static integer c__10007 = 10007;
static integer c__7 = 7;


/* ---  sample program illustrating the computation of small matrix */
/*     exponentials in full with Expokit. Refer to the Expokit */
/*     documentation for more details about the methods, and */
/*     especially the domain of applicability of the Chebyshev scheme. */

/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_9000[] = "(/,a,/,a)";
    static char fmt_9001[] = "(5(1x,d11.4))";
//    static char fmt_9001[] = "(<mprint>(1x,d11.4))";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *);

    /* Local variables */
    static doublereal a[2500]	/* was [50][50] */, h__[2500]	/* was [50][
	    50] */;
    static integer i__, j, k, m;
    static doublereal t, y[50], s1, s2;
    static doublecomplex hc[2500]	/* was [50][50] */, yc[50];
    static integer ns;
    static doublereal wsp[10007];
    static integer iexp;
    static doublecomplex wspc[10007];
    static integer iwsp[50], iflag, iseed[4];
    extern /* Subroutine */ int dgpadm_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *), dgchbv_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    extern doublereal dlaran_(integer *);
    extern /* Subroutine */ int dnchbv_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *), dspadm_(
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *), dschbv_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), zgpadm_(
	    integer *, integer *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, integer *, 
	    integer *), zhpadm_(integer *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     integer *, integer *, integer *), zgchbv_(integer *, doublereal *
	    , doublecomplex *, integer *, doublecomplex *, doublereal *, 
	    integer *, integer *);
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
    static cilist io___29 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___30 = { 0, 6, 0, 0, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___35 = { 0, 6, 0, 0, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___40 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___41 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___42 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___43 = { 0, 6, 0, fmt_9001, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };
    static cilist io___50 = { 0, 6, 0, fmt_9000, 0 };
    static cilist io___51 = { 0, 6, 0, 0, 0 };




/* ------------------------- */
/*     REAL CASE */
/* ------------------------- */
/* ---  set A = random symmetric negative define matrix ... */
    t = 2.;
    m = 5;
    iseed[0] = 3;
    iseed[1] = 7;
    iseed[2] = 3;
    iseed[3] = 7;
    i__1 = m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = m;
	for (i__ = j; i__ <= i__2; ++i__) {
	    a[i__ + j * 50 - 51] = dlaran_(iseed);
	    a[j + i__ * 50 - 51] = a[i__ + j * 50 - 51];
	}
	a[j + j * 50 - 51] += -2.5;
    }
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
	    do_fio(&c__1, (char *)&a[i__ + j * 50 - 51], (ftnlen)sizeof(
		    doublereal));
	}
    }
    e_wsfe();

/* ---  Some compliers (e.g., g77) generate 'Unsupported FORMAT specifier' */
/*     with the specification above. In this case, simply use this form: */
/* 9001 format( 5(1X,D11.4) ) */
/* ---  Pade ... */
    dgpadm_(&c__6, &m, &t, a, &c__50, wsp, &c__10007, iwsp, &iexp, &ns, &
	    iflag);
    s_wsfe(&io___18);
    do_fio(&c__1, "With DGPADM:", (ftnlen)12);
    do_fio(&c__1, "exp(t*A) =", (ftnlen)10);
    e_wsfe();
    s_wsfe(&io___19);
    i__2 = mprint;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = mprint;
	for (j = 1; j <= i__1; ++j) {
	    do_fio(&c__1, (char *)&wsp[iexp + (j - 1) * m + i__ - 2], (ftnlen)
		    sizeof(doublereal));
	}
    }
/*
    e_wsfe();
    dspadm_(&c__6, &m, &t, a, &c__50, wsp, &c__10007, iwsp, &iexp, &ns, &
	    iflag);
    s_wsfe(&io___20);
    do_fio(&c__1, "With DSPADM:", (ftnlen)12);
    do_fio(&c__1, "exp(t*A) =", (ftnlen)10);
    e_wsfe();
    s_wsfe(&io___21);
    i__1 = mprint;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = mprint;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&wsp[iexp + (j - 1) * m + i__ - 2], (ftnlen)
		    sizeof(doublereal));
	}
    }
    e_wsfe();
*/
/* ---  Chebyshev */
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	y[i__ - 1] = 0.;
    }
    y[0] = 1.;
    dgchbv_(&m, &t, a, &c__50, y, wsp, iwsp, &iflag);
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
/*
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	y[i__ - 1] = 0.;
    }
    y[0] = 1.;
    dschbv_(&m, &t, a, &c__50, y, wsp, iwsp, &iflag);
    s_wsfe(&io___25);
    do_fio(&c__1, "With DSCHBV:", (ftnlen)12);
    do_fio(&c__1, "exp(t*A)e_1 =", (ftnlen)13);
    e_wsfe();
    i__2 = mprint;
    for (i__ = 1; i__ <= i__2; ++i__) {
	s_wsle(&io___26);
	do_lio(&c__5, &c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(doublereal));
	e_wsle();
    }*/
/* ---  set H = upper Hessenberg part of A ... */
    i__2 = m;
    for (j = 1; j <= i__2; ++j) {
/* Computing MIN */
	i__3 = j + 1;
	i__1 = min(i__3,m);
	for (i__ = 1; i__ <= i__1; ++i__) {
	    h__[i__ + j * 50 - 51] = a[i__ + j * 50 - 51];
	}
	i__1 = m;
	for (k = i__; k <= i__1; ++k) {
	    h__[k + j * 50 - 51] = 0.;
	}
    }
    s_wsfe(&io___29);
    do_fio(&c__1, "REAL UPPER HESSENBERG CASE", (ftnlen)26);
    do_fio(&c__1, "************************", (ftnlen)24);
    e_wsfe();
    s_wsle(&io___30);
    do_lio(&c__9, &c__1, "H =", (ftnlen)3);
    e_wsle();
    s_wsfe(&io___31);
    i__2 = mprint;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = mprint;
	for (j = 1; j <= i__1; ++j) {
	    do_fio(&c__1, (char *)&h__[i__ + j * 50 - 51], (ftnlen)sizeof(
		    doublereal));
	}
    }
    e_wsfe();
/* ---  Pade ... */
    dgpadm_(&c__6, &m, &t, h__, &c__50, wsp, &c__10007, iwsp, &iexp, &ns, &
	    iflag);
    s_wsfe(&io___32);
    do_fio(&c__1, "With DGPADM:", (ftnlen)12);
    do_fio(&c__1, "exp(t*H) =", (ftnlen)10);
    e_wsfe();
    s_wsfe(&io___33);
    i__1 = mprint;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = mprint;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&wsp[iexp + (j - 1) * m + i__ - 2], (ftnlen)
		    sizeof(doublereal));
	}
    }
    e_wsfe();
/* ---  Chebyshev */
    /*i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	y[i__ - 1] = 0.;
    }
    y[0] = 1.;
    dnchbv_(&m, &t, a, &c__50, y, wsp, &iflag);
    s_wsfe(&io___34);
    do_fio(&c__1, "With DNCHBV:", (ftnlen)12);
    do_fio(&c__1, "exp(t*A)e_1 =", (ftnlen)13);
    e_wsfe();
    i__2 = mprint;
    for (i__ = 1; i__ <= i__2; ++i__) {
	s_wsle(&io___35);
	do_lio(&c__5, &c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(doublereal));
	e_wsle();
    }*/

    return 0;
} /* MAIN__ */

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


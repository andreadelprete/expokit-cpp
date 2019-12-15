/* loadhb.f -- translated by f2c (version 20100827).
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
static integer c__9 = 9;
static integer c__3 = 3;

/* ----------------------------------------------------------------------| */
/* Subroutine */ int loadhb_(char *filename, char *spformat, integer *n, 
	integer *nz, integer *ia, integer *ja, doublereal *a, integer *iwsp, 
	ftnlen filename_len, ftnlen spformat_len)
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
	    integer *, char *, ftnlen), e_wsle(void), f_clos(cllist *), s_cmp(
	    char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, io, nn;
    static char key[8];
    static integer nnz, nrhs;
    static char type__[3];
    static integer nrow;
    static char title[72];
    static integer indcrd, valcrd;
    static char indfmt[16];
    extern /* Subroutine */ int dgcnvr_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, ftnlen, ftnlen, ftnlen);
    static integer rhscrd;
    static char valfmt[20];
    static integer ptrcrd, totcrd;
    static char rhsfmt[20];
    static integer nrhsix;
    static char ptrfmt[16], rhstyp[1];

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 7, 0, fmt_10, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 7, 0, fmt_11, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };
    static cilist io___26 = { 0, 7, 0, ptrfmt, 0 };
    static cilist io___27 = { 0, 7, 0, indfmt, 0 };
    static cilist io___28 = { 0, 7, 0, valfmt, 0 };
    static cilist io___31 = { 0, 6, 0, 0, 0 };


/* ---  Purpose ---------------------------------------------------------| */

/* ---  LOADHB loads a matrix stored under the Harwell-Boeing format */
/*     and renders it into the sparse format specified by spformat. */

/* ---  Arguments -------------------------------------------------------| */

/*     filename (input) */
/*           name of the file containing the matrix. */
/*           must end with a '$', i.e., filename is in the form: '...$' */

/*     spformat (input) */
/*           sparse format in which the matrix is forced to be on output */
/*           accepted values are: */
/*              'COO' : COOrdinates */
/*              'CRS' : Compressed Row Storage */
/*              'CCS' : Compressed Column Storage (default H-B format) */

/*     n (input/output) */
/*           On input,  the maximum allowable order */
/*           On output, the actual order of the matrix loaded */

/*     nz (input/output) */
/*           On input,  the maximum allowable number of non zero entries */
/*           On output, the actual number of non zero entries */

/*     ia,ja,a (output) */
/*           sparse matrix data stored in the format given in spformat */
/*           sufficient room is needed to achieve this: each component */
/*           must be of length >= nz. If the matrix is symmetric, both */
/*           lower and upper parts are included explicitly */

/*     iwsp (workspace) of length >= n */

/* ----------------------------------------------------------------------| */

/* --- */
    /* Parameter adjustments */
    --iwsp;
    --a;
    --ja;
    --ia;

    /* Function Body */
    i__ = i_indx(filename, "$", (ftnlen)80, (ftnlen)1) - 1;
    if (i__ <= 0) {
	s_stop("in LOADHB. Bad filename", (ftnlen)23);
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
    s_rsfe(&io___3);
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
    s_wsle(&io___20);
    do_lio(&c__9, &c__1, title, (ftnlen)72);
    do_lio(&c__9, &c__1, "type :", (ftnlen)6);
    do_lio(&c__9, &c__1, type__, (ftnlen)3);
    do_lio(&c__9, &c__1, " size :", (ftnlen)7);
    do_lio(&c__3, &c__1, (char *)&nrow, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nn, (ftnlen)sizeof(integer));
    e_wsle();
    s_wsle(&io___21);
    do_lio(&c__9, &c__1, "order :", (ftnlen)7);
    do_lio(&c__3, &c__1, (char *)&nn, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " number of nonzero :", (ftnlen)20);
    do_lio(&c__3, &c__1, (char *)&nnz, (ftnlen)sizeof(integer));
    e_wsle();
    if (nn > *n) {
	s_stop("in LOADHB. Please increase n", (ftnlen)28);
    }
    if (nnz > *nz) {
	s_stop("in LOADHB. Please increase nz", (ftnlen)29);
    }
/* ---  leave if there is no values ... */
    if (valcrd <= 0) {
	s_stop("Empty Harwell-Boeing matrix", (ftnlen)27);
    }
    if (rhscrd > 0) {
	s_rsfe(&io___22);
	do_fio(&c__1, rhstyp, (ftnlen)1);
	do_fio(&c__1, (char *)&nrhs, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nrhsix, (ftnlen)sizeof(integer));
	e_rsfe();
	s_wsle(&io___25);
	do_lio(&c__9, &c__1, "There is a second hand", (ftnlen)22);
	e_wsle();
    }
/* ---  read data... */
    s_rsfe(&io___26);
    i__1 = nn + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ja[i__], (ftnlen)sizeof(integer));
    }
    e_rsfe();
    s_rsfe(&io___27);
    i__1 = nnz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ia[i__], (ftnlen)sizeof(integer));
    }
    e_rsfe();
    s_rsfe(&io___28);
    i__1 = nnz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&a[i__], (ftnlen)sizeof(doublereal));
    }
    e_rsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 7;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* ---  for the sake of experiments, store both parts if symmetric matrix */
    if (s_cmp(type__, "RSA", (ftnlen)3, (ftnlen)3) == 0) {
/* ---     expand ja indices ... */
	k = nnz;
	for (j = nn; j >= 1; --j) {
	    i__1 = ja[j + 1] - ja[j];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ja[k] = j;
		--k;
	    }
	}
/* ---     insert the other half ... */
	k = nnz;
	i__1 = k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (ia[i__] != ja[i__]) {
		++nnz;
		if (nnz > *nz) {
		    s_stop("in LOADHB. Please increase nz", (ftnlen)29);
		}
		ia[nnz] = ja[i__];
		ja[nnz] = ia[i__];
		a[nnz] = a[i__];
	    }
	}
	s_copy(type__, "COO", (ftnlen)3, (ftnlen)3);
    } else {
	s_copy(type__, "CCS", (ftnlen)3, (ftnlen)3);
    }
    dgcnvr_(type__, spformat, "n", &nn, &nn, &nnz, &ia[1], &ja[1], &a[1], &
	    iwsp[1], (ftnlen)3, (ftnlen)3, (ftnlen)1);
    *n = nn;
    *nz = nnz;
    s_wsle(&io___31);
    do_lio(&c__9, &c__1, "Harwell-Boeing matrix loaded", (ftnlen)28);
    e_wsle();
    return 0;
} /* loadhb_ */


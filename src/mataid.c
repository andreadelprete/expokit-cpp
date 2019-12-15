/* mataid.f -- translated by f2c (version 20100827).
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

struct {
    doublecomplex a[50000];
    integer ia[50000], ja[50000], nz, n;
} cmat_;

#define cmat_1 cmat_

/* Table of constant values */

static integer c__1 = 1;
static integer c__9 = 9;
static integer c__3 = 3;

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
/* Subroutine */ int dgcoov_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;


/* ---  Computes y = A*x. A is passed via a fortran `common statement'. */
/* ---  A is assumed here to be under the COOrdinates storage format. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = rmat_1.n;
    for (j = 1; j <= i__1; ++j) {
	y[j] = 0.;
    }
    i__1 = rmat_1.nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[rmat_1.ia[i__ - 1]] += rmat_1.a[i__ - 1] * x[rmat_1.ja[i__ - 1]];
    }
    return 0;
} /* dgcoov_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* Subroutine */ int dgcrsv_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;


/* ---  Computes y = A*x. A is passed via a fortran `common statement'. */
/* ---  A is assumed to be under the Compress Row Storage (CRS) format. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = 0.;
	i__2 = rmat_1.ia[i__] - 1;
	for (j = rmat_1.ia[i__ - 1]; j <= i__2; ++j) {
	    y[i__] += rmat_1.a[j - 1] * x[rmat_1.ja[j - 1]];
	}
    }
    return 0;
} /* dgcrsv_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* Subroutine */ int dgccsv_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;


/* ---  Computes y = A*x. A is passed via a fortran `common statement'. */
/* ---  A is assumed to be under the Compress Column Storage (CCS) format. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = rmat_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = 0.;
    }
    i__1 = rmat_1.n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = rmat_1.ja[j] - 1;
	for (i__ = rmat_1.ja[j - 1]; i__ <= i__2; ++i__) {
	    y[rmat_1.ia[i__ - 1]] += rmat_1.a[i__ - 1] * x[j];
	}
    }
    return 0;
} /* dgccsv_ */

/* ----------------------------------------------------------------------| */
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

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* Subroutine */ int dgcnvr_(char *from, char *to, char *diag, integer *nrow, 
	integer *ncol, integer *nz, integer *ia, integer *ja, doublereal *a, 
	integer *iwsp, ftnlen from_len, ftnlen to_len, ftnlen diag_len)
{
    /* System generated locals */
    integer i__1;
    char ch__1[3], ch__2[3];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    logical l_le(char *, char *, ftnlen, ftnlen), l_ge(char *, char *, ftnlen,
	     ftnlen);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static char cpy_diag__[1], cpy_from__[3], c__[1];
    static integer i__, j, k, nn, iflag;
    extern /* Subroutine */ int dcmpac_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *);
    static char cpy_to__[3];

/* -----Purpose----------------------------------------------------------| */

/* ---  DGCNVR converts a sparse storage format into another sparse */
/*     storage format. The matrix can be rectangular. */

/* -----Arguments--------------------------------------------------------| */

/*     from   : (input, character*3) the storage format holding the */
/*              matrix on input. Accepted values are: */
/*              'COO' : COOrdinates */
/*              'CRS' : Compressed Row Storage */
/*              'CCS' : Compressed Column Storage. */

/*     to     : (input, character*3) the storage format holding the */
/*              matrix on output. Same accepted values as above. */

/*     diag   : (input, character*1) specifies whether or not the */
/*              entire diagonal should be lodged, i.e., including */
/*              null diagonal entries. (This may be required by */
/*              some applications that make use of the locations */
/*              of diagonal elements.) */
/*              diag = 'N', no specific attention is given to diagonal */
/*                          elements. */
/*              diag = 'D', the entire diagonal is lodged, including */
/*                          null elements on the diagonal. */
/*              if from=to & diag='D', null diagonal entries are */
/*              explicitly inserted in the actual matrix. */

/*     nrow   : (input) number of rows in the matrix. */

/*     ncol   : (input) number of columns in the matrix. */

/*     nz     : (input/output) number of non-zeros elements. */
/*              If diag='D' and null diagonal entries are inserted, then */
/*              nz is updated on exit and contains the effective number */
/*              of entries stored. In what follows, nz' (read nz prime) */
/*              denotes the updated value of nz, nz <= nz' <= nz+n. */
/*              If diag='N' then nz'=nz. */

/*     ia(*) : (input/output) of declared length .ge. nz'. */
/*             On input, */
/*                if from='CRS', ia(1:nrow+1) contains pointers for the */
/*                beginning of each row. */
/*                if from='COO', or 'CCS', ia(1:nz) contains row indices. */
/*             On output, */
/*                if to='CRS', ia(1:nrow+1) contains pointers for the */
/*                beginning of each row. */
/*                if to='COO' or 'CCS', ia(1:nz') contains row indices */
/*                in increasing order in each column. */

/*     ja(*) : (input/output) of declared length .ge. nz'. */
/*             On input, */
/*                if from='CRS', ja(1:ncol+1) contains pointers for the */
/*                beginning of each column. */
/*                if from='COO' or 'CCS', ja(1:nz) contains col. indices. */
/*             On output, */
/*                if to='CRS', ja(1:ncol+1) contains pointers for the */
/*                beginning of each column. */
/*                if to='COO' or 'CCS', ja(1:nz') contains column indices */
/*                in increasing order in each row. */

/*     a(*)  : (input/output) On input, a(1:nz) fits the input format */
/*             and on output, a(1:nz') fits the output format. */

/*   iwsp(*) : (workspace) of declared length .ge. max(nrow,ncol). */

/* ----------------------------------------------------------------------| */
/*     Roger B. Sidje (rbs@maths.uq.edu.au) */
/*     Department of Mathematics, University of Queensland. */
/*     Brisbane QLD 4072, Australia. 1996. */
/* ----------------------------------------------------------------------| */

/* ---  upper case strings ... */

    /* Parameter adjustments */
    --iwsp;
    --a;
    --ja;
    --ia;

    /* Function Body */
    for (k = 1; k <= 3; ++k) {
	*(unsigned char *)c__ = *(unsigned char *)&from[k - 1];
	if (*(unsigned char *)c__ >= 'a' && *(unsigned char *)c__ <= 'z') {
	    *(unsigned char *)c__ = (char) (*(unsigned char *)c__ - 32);
	}
	*(unsigned char *)&cpy_from__[k - 1] = *(unsigned char *)c__;
	*(unsigned char *)c__ = *(unsigned char *)&to[k - 1];
	if (*(unsigned char *)c__ >= 'a' && *(unsigned char *)c__ <= 'z') {
	    *(unsigned char *)c__ = (char) (*(unsigned char *)c__ - 32);
	}
	*(unsigned char *)&cpy_to__[k - 1] = *(unsigned char *)c__;
    }
    *(unsigned char *)c__ = *(unsigned char *)diag;
    if (*(unsigned char *)c__ >= 'a' && *(unsigned char *)c__ <= 'z') {
	*(unsigned char *)c__ = (char) (*(unsigned char *)c__ - 32);
    }
    *(unsigned char *)cpy_diag__ = *(unsigned char *)c__;
/* ---  quick return if possible ... */
    s_copy(ch__1, cpy_from__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, cpy_to__, (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3) && *(unsigned char *)cpy_diag__ == 'N') {
	return 0;
    }
    iflag = 1;
    s_copy(ch__1, cpy_from__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "COO", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	iflag = 0;
    }
    s_copy(ch__1, cpy_from__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CRS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	iflag = 0;
    }
    s_copy(ch__1, cpy_from__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CCS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	iflag = 0;
    }
    if (iflag != 0) {
	s_stop("unexpected i/o formats (in DGCNVR)", (ftnlen)34);
    }
    iflag = 1;
    s_copy(ch__1, cpy_to__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "COO", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	iflag = 0;
    }
    s_copy(ch__1, cpy_to__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CRS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	iflag = 0;
    }
    s_copy(ch__1, cpy_to__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CCS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	iflag = 0;
    }
    if (iflag != 0) {
	s_stop("unexpected i/o formats (in DGCNVR)", (ftnlen)34);
    }

/* ---  transit via COOrdinates format if input is not in COO ... */

    s_copy(ch__1, cpy_from__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CRS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
/* ---     expand ia indices ... */
	nn = *nz;
	for (i__ = *nrow; i__ >= 1; --i__) {
	    i__1 = ia[i__ + 1] - ia[i__];
	    for (k = 1; k <= i__1; ++k) {
		ia[nn] = i__;
		--nn;
	    }
	}
    }
    s_copy(ch__1, cpy_from__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CCS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
/* ---     expand ja indices ... */
	nn = *nz;
	for (j = *ncol; j >= 1; --j) {
	    i__1 = ja[j + 1] - ja[j];
	    for (k = 1; k <= i__1; ++k) {
		ja[nn] = j;
		--nn;
	    }
	}
    }

/* --   if requested, insert diagonal elements even if they are zero... */

    if (*(unsigned char *)cpy_diag__ == 'D') {
	nn = min(*nrow,*ncol);
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iwsp[i__] = 0;
	}
	i__1 = *nz;
	for (k = 1; k <= i__1; ++k) {
	    if (ia[k] == ja[k]) {
		iwsp[ia[k]] = 1;
	    }
	}
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (iwsp[i__] == 0) {
		++(*nz);
		ia[*nz] = i__;
		ja[*nz] = i__;
		a[*nz] = 0.;
	    }
	}
    }
/* ---  COO convertion ... */
    s_copy(ch__1, cpy_to__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "COO", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	return 0;
    }
/* ---  CRS convertion ... */
    s_copy(ch__1, cpy_to__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CRS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	dcmpac_(nrow, nz, &ia[1], &ja[1], &a[1], &iwsp[1]);
    }
/* ---  CCS convertion ... */
    s_copy(ch__1, cpy_to__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CCS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	dcmpac_(ncol, nz, &ja[1], &ia[1], &a[1], &iwsp[1]);
    }
    return 0;
} /* dgcnvr_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */

/* Subroutine */ int dcmpac_(integer *n, integer *nx, integer *ix, integer *
	ixx, doublereal *xx, integer *iwsp)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int idsrt1_(integer *, integer *, doublereal *), 
	    idsrt2_(integer *, integer *, integer *, doublereal *);

/* --   DCMPAC compacts the array ix and sorts ixx and xx */
/* --   (This is a gateway routine for DGCNVR) ... */
/* ----------------------------------------------------------------------| */

/* ---  sort ix and carry ixx and xx along ... */

    /* Parameter adjustments */
    --iwsp;
    --xx;
    --ixx;
    --ix;

    /* Function Body */
    idsrt2_(nx, &ix[1], &ixx[1], &xx[1]);

/* ---  adjust pointers ... */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	iwsp[k] = 0;
    }
    i__1 = *nx;
    for (k = 1; k <= i__1; ++k) {
	++iwsp[ix[k]];
    }
    ix[*n + 1] = *nx + 1;
    for (k = *n; k >= 1; --k) {
	ix[k] = ix[k + 1] - iwsp[k];
    }

/* ---  sort ixx in increasing order and carry xx along ... */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	idsrt1_(&iwsp[k], &ixx[ix[k]], &xx[ix[k]]);
    }
    return 0;
} /* dcmpac_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */

/* Subroutine */ int idsrt1_(integer *nx, integer *ix, doublereal *xx)
{
    static integer i__, j, k, l, m;
    static doublereal r__;
    static integer ij, il[21], it, iu[21];
    static doublereal tx;
    static integer iit;
    static doublereal ttx;

/* ---  IDSRT1: indirect sort -- sort ix and carry xx along */
/* ---  adapted from a SLAP (Sparse Linear Algebra Package) code. */
/* ----------------------------------------------------------------------| */
    /* Parameter adjustments */
    --xx;
    --ix;

    /* Function Body */
    if (*nx <= 1) {
	return 0;
    }
/* ---  And now...Just a little black magic... */
    m = 1;
    i__ = 1;
    j = *nx;
    r__ = .375f;
L210:
    if (r__ <= .5898437f) {
	r__ += .0390625f;
    } else {
	r__ += -.21875f;
    }
L225:
    k = i__;

/* ---  Select a central element of the array and save it in location */
/* ---  IT, TX. */

    ij = i__ + (integer) ((doublereal) (j - i__) * r__);
    it = ix[ij];
    tx = xx[ij];

/* ---  If first element of array is greater than IT, interchange with IT. */

    if (ix[i__] > it) {
	ix[ij] = ix[i__];
	ix[i__] = it;
	it = ix[ij];
	xx[ij] = xx[i__];
	xx[i__] = tx;
	tx = xx[ij];
    }
    l = j;

/* ---  If last element of array is less than IT, swap with IT. */

    if (ix[j] < it) {
	ix[ij] = ix[j];
	ix[j] = it;
	it = ix[ij];
	xx[ij] = xx[j];
	xx[j] = tx;
	tx = xx[ij];

/* ---  If first element of array is greater than IT, swap with IT. */

	if (ix[i__] > it) {
	    ix[ij] = ix[i__];
	    ix[i__] = it;
	    it = ix[ij];
	    xx[ij] = xx[i__];
	    xx[i__] = tx;
	    tx = xx[ij];
	}
    }

/* ---  Find an element in the second half of the array which is */
/* ---  smaller than IT. */

L240:
    --l;
    if (ix[l] > it) {
	goto L240;
    }

/* ---  Find an element in the first half of the array which is */
/* ---  greater than IT. */

L245:
    ++k;
    if (ix[k] < it) {
	goto L245;
    }

/* ---  Interchange these elements. */

    if (k <= l) {
	iit = ix[l];
	ix[l] = ix[k];
	ix[k] = iit;
	ttx = xx[l];
	xx[l] = xx[k];
	xx[k] = ttx;
	goto L240;
    }

/* ---  Save upper and lower subscripts of the array yet to be sorted. */

    if (l - i__ > j - k) {
	il[m - 1] = i__;
	iu[m - 1] = l;
	i__ = k;
	++m;
    } else {
	il[m - 1] = k;
	iu[m - 1] = j;
	j = l;
	++m;
    }
    goto L260;

/* ---  Begin again on another portion of the unsorted array. */

L255:
    --m;
    if (m == 0) {
	goto L300;
    }
    i__ = il[m - 1];
    j = iu[m - 1];
L260:
    if (j - i__ >= 1) {
	goto L225;
    }
    if (i__ == j) {
	goto L255;
    }
    if (i__ == 1) {
	goto L210;
    }
    --i__;
L265:
    ++i__;
    if (i__ == j) {
	goto L255;
    }
    it = ix[i__ + 1];
    tx = xx[i__ + 1];
    if (ix[i__] <= it) {
	goto L265;
    }
    k = i__;
L270:
    ix[k + 1] = ix[k];
    xx[k + 1] = xx[k];
    --k;
    if (it < ix[k]) {
	goto L270;
    }
    ix[k + 1] = it;
    xx[k + 1] = tx;
    goto L265;
L300:
    return 0;
} /* idsrt1_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* Subroutine */ int idsrt2_(integer *nx, integer *ix, integer *ixx, 
	doublereal *xx)
{
    static integer i__, j, k, l, m;
    static doublereal r__;
    static integer ij, il[21], it, iu[21], jt;
    static doublereal tx;
    static integer iit, jjt;
    static doublereal ttx;

/* ---  IDSRT2: indirect sort: sort ix and carry ixx and xx along */
/* ---  adapted from a SLAP (Sparse Linear Algebra Package) code. */
/* ----------------------------------------------------------------------| */
    /* Parameter adjustments */
    --xx;
    --ixx;
    --ix;

    /* Function Body */
    if (*nx <= 1) {
	return 0;
    }
/* ---  And now...Just a little black magic... */
    m = 1;
    i__ = 1;
    j = *nx;
    r__ = .375f;
L210:
    if (r__ <= .5898437f) {
	r__ += .0390625f;
    } else {
	r__ += -.21875f;
    }
L225:
    k = i__;

/* ---  Select a central element of the array and save it in location */
/* ---  IT, JT, TX. */

    ij = i__ + (integer) ((doublereal) (j - i__) * r__);
    it = ix[ij];
    jt = ixx[ij];
    tx = xx[ij];

/* ---  If first element of array is greater than IT, interchange with IT. */

    if (ix[i__] > it) {
	ix[ij] = ix[i__];
	ix[i__] = it;
	it = ix[ij];
	ixx[ij] = ixx[i__];
	ixx[i__] = jt;
	jt = ixx[ij];
	xx[ij] = xx[i__];
	xx[i__] = tx;
	tx = xx[ij];
    }
    l = j;

/* ---  If last element of array is less than IT, swap with IT. */

    if (ix[j] < it) {
	ix[ij] = ix[j];
	ix[j] = it;
	it = ix[ij];
	ixx[ij] = ixx[j];
	ixx[j] = jt;
	jt = ixx[ij];
	xx[ij] = xx[j];
	xx[j] = tx;
	tx = xx[ij];

/* ---  If first element of array is greater than IT, swap with IT. */

	if (ix[i__] > it) {
	    ix[ij] = ix[i__];
	    ix[i__] = it;
	    it = ix[ij];
	    ixx[ij] = ixx[i__];
	    ixx[i__] = jt;
	    jt = ixx[ij];
	    xx[ij] = xx[i__];
	    xx[i__] = tx;
	    tx = xx[ij];
	}
    }

/* ---  Find an element in the second half of the array which is */
/* ---  smaller than IT. */

L240:
    --l;
    if (ix[l] > it) {
	goto L240;
    }

/* ---  Find an element in the first half of the array which is */
/* ---  greater than IT. */

L245:
    ++k;
    if (ix[k] < it) {
	goto L245;
    }

/* ---  Interchange these elements. */

    if (k <= l) {
	iit = ix[l];
	ix[l] = ix[k];
	ix[k] = iit;
	jjt = ixx[l];
	ixx[l] = ixx[k];
	ixx[k] = jjt;
	ttx = xx[l];
	xx[l] = xx[k];
	xx[k] = ttx;
	goto L240;
    }

/* ---  Save upper and lower subscripts of the array yet to be sorted. */

    if (l - i__ > j - k) {
	il[m - 1] = i__;
	iu[m - 1] = l;
	i__ = k;
	++m;
    } else {
	il[m - 1] = k;
	iu[m - 1] = j;
	j = l;
	++m;
    }
    goto L260;

/* ---  Begin again on another portion of the unsorted array. */

L255:
    --m;
    if (m == 0) {
	goto L300;
    }
    i__ = il[m - 1];
    j = iu[m - 1];
L260:
    if (j - i__ >= 1) {
	goto L225;
    }
    if (i__ == j) {
	goto L255;
    }
    if (i__ == 1) {
	goto L210;
    }
    --i__;
L265:
    ++i__;
    if (i__ == j) {
	goto L255;
    }
    it = ix[i__ + 1];
    jt = ixx[i__ + 1];
    tx = xx[i__ + 1];
    if (ix[i__] <= it) {
	goto L265;
    }
    k = i__;
L270:
    ix[k + 1] = ix[k];
    ixx[k + 1] = ixx[k];
    xx[k + 1] = xx[k];
    --k;
    if (it < ix[k]) {
	goto L270;
    }
    ix[k + 1] = it;
    ixx[k + 1] = jt;
    xx[k + 1] = tx;
    goto L265;
L300:
    return 0;
} /* idsrt2_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* Subroutine */ int zgcnvr_(char *from, char *to, char *diag, integer *nrow, 
	integer *ncol, integer *nz, integer *ia, integer *ja, doublecomplex *
	a, integer *iwsp, ftnlen from_len, ftnlen to_len, ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2;
    char ch__1[3], ch__2[3];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    logical l_le(char *, char *, ftnlen, ftnlen), l_ge(char *, char *, ftnlen,
	     ftnlen);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static char cpy_diag__[1], cpy_from__[3], c__[1];
    static integer i__, j, k, nn, iflag;
    extern /* Subroutine */ int zcmpac_(integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *);
    static char cpy_to__[3];

/* -----Purpose----------------------------------------------------------| */

/* ---  ZGCNVR transforms a sparse storage format into another sparse */
/*     storage format. The matrix can be rectangular. */
/*     This is the complex counterpart of DGCNVR. */

/* -----Arguments--------------------------------------------------------| */

/*     from   : (input, character*3) the storage format holding the */
/*              matrix on input. Accepted values are: */
/*              'COO' : COOrdinates */
/*              'CRS' : Compressed Row Storage */
/*              'CCS' : Compressed Column Storage. */

/*     to     : (input, character*3) the storage format holding the */
/*              matrix on output. Same accepted values as above. */

/*     diag   : (input, character*1) specifies whether or not the */
/*              entire diagonal should be lodged, i.e., including */
/*              null diagonal entries. (This may be required by */
/*              some applications that make use of the locations */
/*              of diagonal elements.) */
/*              diag = 'N', no specific attention is given to diagonal */
/*                          elements. */
/*              diag = 'D', the entire diagonal is lodged, including */
/*                          null elements on the diagonal. */
/*              if from=to & diag='D', null diagonal entries are */
/*              explicitly inserted in the actual matrix. */

/*     nrow   : (input) number of rows in the matrix. */

/*     ncol   : (input) number of columns in the matrix. */

/*     nz     : (input/output) number of non-zeros elements. */
/*              If diag='D' and null diagonal entries are inserted, then */
/*              nz is updated on exit and contains the effective number */
/*              of entries stored. In what follows, nz' (read nz prime) */
/*              denotes the updated value of nz, nz <= nz' <= nz+n. */
/*              If diag='N' then nz'=nz. */

/*     ia(*) : (input/output) of declared length .ge. nz'. */
/*             On input, */
/*                if from='CRS', ia(1:nrow+1) contains pointers for the */
/*                beginning of each row. */
/*                if from='COO', or 'CCS', ia(1:nz) contains row indices. */
/*             On output, */
/*                if to='CRS', ia(1:nrow+1) contains pointers for the */
/*                beginning of each row. */
/*                if to='COO' or 'CCS', ia(1:nz') contains row indices */
/*                in increasing order in each column. */

/*     ja(*) : (input/output) of declared length .ge. nz'. */
/*             On input, */
/*                if from='CRS', ja(1:ncol+1) contains pointers for the */
/*                beginning of each column. */
/*                if from='COO' or 'CCS', ja(1:nz) contains col. indices. */
/*             On output, */
/*                if to='CRS', ja(1:ncol+1) contains pointers for the */
/*                beginning of each column. */
/*                if to='COO' or 'CCS', ja(1:nz') contains column indices */
/*                in increasing order in each row. */

/*     a(*)  : (input/output) On input, a(1:nz) fits the input format */
/*             and on output, a(1:nz') fits the output format. */

/*   iwsp(*) : (workspace) of declared length .ge. max(nrow,ncol). */

/* ----------------------------------------------------------------------| */
/*     Roger B. Sidje (rbs@maths.uq.edu.au) */
/*     Department of Mathematics, University of Queensland. */
/*     Brisbane QLD 4072, Australia. 1996. */
/* ----------------------------------------------------------------------| */

/* ---  upper case strings ... */

    /* Parameter adjustments */
    --iwsp;
    --a;
    --ja;
    --ia;

    /* Function Body */
    for (k = 1; k <= 3; ++k) {
	*(unsigned char *)c__ = *(unsigned char *)&from[k - 1];
	if (*(unsigned char *)c__ >= 'a' && *(unsigned char *)c__ <= 'z') {
	    *(unsigned char *)c__ = (char) (*(unsigned char *)c__ - 32);
	}
	*(unsigned char *)&cpy_from__[k - 1] = *(unsigned char *)c__;
	*(unsigned char *)c__ = *(unsigned char *)&to[k - 1];
	if (*(unsigned char *)c__ >= 'a' && *(unsigned char *)c__ <= 'z') {
	    *(unsigned char *)c__ = (char) (*(unsigned char *)c__ - 32);
	}
	*(unsigned char *)&cpy_to__[k - 1] = *(unsigned char *)c__;
    }
    *(unsigned char *)c__ = *(unsigned char *)diag;
    if (*(unsigned char *)c__ >= 'a' && *(unsigned char *)c__ <= 'z') {
	*(unsigned char *)c__ = (char) (*(unsigned char *)c__ - 32);
    }
    *(unsigned char *)cpy_diag__ = *(unsigned char *)c__;
/* ---  quick return if possible ... */
    s_copy(ch__1, cpy_from__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, cpy_to__, (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3) && *(unsigned char *)cpy_diag__ == 'N') {
	return 0;
    }
    iflag = 1;
    s_copy(ch__1, cpy_from__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "COO", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	iflag = 0;
    }
    s_copy(ch__1, cpy_from__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CRS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	iflag = 0;
    }
    s_copy(ch__1, cpy_from__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CCS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	iflag = 0;
    }
    if (iflag != 0) {
	s_stop("unexpected i/o formats (in ZGCNVR)", (ftnlen)34);
    }
    iflag = 1;
    s_copy(ch__1, cpy_to__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "COO", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	iflag = 0;
    }
    s_copy(ch__1, cpy_to__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CRS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	iflag = 0;
    }
    s_copy(ch__1, cpy_to__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CCS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	iflag = 0;
    }
    if (iflag != 0) {
	s_stop("unexpected i/o formats (in ZGCNVR)", (ftnlen)34);
    }

/* ---  transit via COOrdinates format if input is not in COO ... */

    s_copy(ch__1, cpy_from__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CRS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
/* ---     expand ia indices ... */
	nn = *nz;
	for (i__ = *nrow; i__ >= 1; --i__) {
	    i__1 = ia[i__ + 1] - ia[i__];
	    for (k = 1; k <= i__1; ++k) {
		ia[nn] = i__;
		--nn;
	    }
	}
    }
    s_copy(ch__1, cpy_from__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CCS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
/* ---     expand ja indices ... */
	nn = *nz;
	for (j = *ncol; j >= 1; --j) {
	    i__1 = ja[j + 1] - ja[j];
	    for (k = 1; k <= i__1; ++k) {
		ja[nn] = j;
		--nn;
	    }
	}
    }

/* --   if requested, insert diagonal elements even if they are zero... */

    if (*(unsigned char *)cpy_diag__ == 'D') {
	nn = min(*nrow,*ncol);
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iwsp[i__] = 0;
	}
	i__1 = *nz;
	for (k = 1; k <= i__1; ++k) {
	    if (ia[k] == ja[k]) {
		iwsp[ia[k]] = 1;
	    }
	}
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (iwsp[i__] == 0) {
		++(*nz);
		ia[*nz] = i__;
		ja[*nz] = i__;
		i__2 = *nz;
		a[i__2].r = 0., a[i__2].i = 0.;
	    }
	}
    }
/* ---  COO convertion ... */
    s_copy(ch__1, cpy_to__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "COO", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	return 0;
    }
/* ---  CRS convertion ... */
    s_copy(ch__1, cpy_to__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CRS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	zcmpac_(nrow, nz, &ia[1], &ja[1], &a[1], &iwsp[1]);
    }
/* ---  CCS convertion ... */
    s_copy(ch__1, cpy_to__, (ftnlen)3, (ftnlen)3);
    s_copy(ch__2, "CCS", (ftnlen)3, (ftnlen)3);
    if (l_le(ch__1, ch__2, (ftnlen)3, (ftnlen)3) && l_ge(ch__1, ch__2, (
	    ftnlen)3, (ftnlen)3)) {
	zcmpac_(ncol, nz, &ja[1], &ia[1], &a[1], &iwsp[1]);
    }
    return 0;
} /* zgcnvr_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */

/* Subroutine */ int zcmpac_(integer *n, integer *nx, integer *ix, integer *
	ixx, doublecomplex *xx, integer *iwsp)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int izsrt1_(integer *, integer *, doublecomplex *)
	    , izsrt2_(integer *, integer *, integer *, doublecomplex *);

/* --   ZCMPAC compacts the array ix and sorts ixx and xx */
/* --   (This is a gateway routine for ZGCNVR) ... */
/* ----------------------------------------------------------------------| */

/* ---  sort ix and carry ixx and xx along ... */

    /* Parameter adjustments */
    --iwsp;
    --xx;
    --ixx;
    --ix;

    /* Function Body */
    izsrt2_(nx, &ix[1], &ixx[1], &xx[1]);

/* ---  adjust pointers ... */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	iwsp[k] = 0;
    }
    i__1 = *nx;
    for (k = 1; k <= i__1; ++k) {
	++iwsp[ix[k]];
    }
    ix[*n + 1] = *nx + 1;
    for (k = *n; k >= 1; --k) {
	ix[k] = ix[k + 1] - iwsp[k];
    }

/* ---  sort ixx in increasing order and carry xx along ... */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	izsrt1_(&iwsp[k], &ixx[ix[k]], &xx[ix[k]]);
    }
    return 0;
} /* zcmpac_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */

/* Subroutine */ int izsrt1_(integer *nx, integer *ix, doublecomplex *xx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal r__;
    static integer ij, il[21], it, iu[21];
    static doublecomplex tx;
    static integer iit;
    static doublecomplex ttx;

/* ---  IZSRT1: indirect sort -- sort ix and carry xx along */
/* ---  adapted from a SLAP (Sparse Linear Algebra Package) code. */
/* ----------------------------------------------------------------------| */
    /* Parameter adjustments */
    --xx;
    --ix;

    /* Function Body */
    if (*nx <= 1) {
	return 0;
    }
/* ---  And now...Just a little black magic... */
    m = 1;
    i__ = 1;
    j = *nx;
    r__ = .375f;
L210:
    if (r__ <= .5898437f) {
	r__ += .0390625f;
    } else {
	r__ += -.21875f;
    }
L225:
    k = i__;

/* ---  Select a central element of the array and save it in location */
/* ---  IT, TX. */

    ij = i__ + (integer) ((doublereal) (j - i__) * r__);
    it = ix[ij];
    i__1 = ij;
    tx.r = xx[i__1].r, tx.i = xx[i__1].i;

/* ---  If first element of array is greater than IT, interchange with IT. */

    if (ix[i__] > it) {
	ix[ij] = ix[i__];
	ix[i__] = it;
	it = ix[ij];
	i__1 = ij;
	i__2 = i__;
	xx[i__1].r = xx[i__2].r, xx[i__1].i = xx[i__2].i;
	i__1 = i__;
	xx[i__1].r = tx.r, xx[i__1].i = tx.i;
	i__1 = ij;
	tx.r = xx[i__1].r, tx.i = xx[i__1].i;
    }
    l = j;

/* ---  If last element of array is less than IT, swap with IT. */

    if (ix[j] < it) {
	ix[ij] = ix[j];
	ix[j] = it;
	it = ix[ij];
	i__1 = ij;
	i__2 = j;
	xx[i__1].r = xx[i__2].r, xx[i__1].i = xx[i__2].i;
	i__1 = j;
	xx[i__1].r = tx.r, xx[i__1].i = tx.i;
	i__1 = ij;
	tx.r = xx[i__1].r, tx.i = xx[i__1].i;

/* ---  If first element of array is greater than IT, swap with IT. */

	if (ix[i__] > it) {
	    ix[ij] = ix[i__];
	    ix[i__] = it;
	    it = ix[ij];
	    i__1 = ij;
	    i__2 = i__;
	    xx[i__1].r = xx[i__2].r, xx[i__1].i = xx[i__2].i;
	    i__1 = i__;
	    xx[i__1].r = tx.r, xx[i__1].i = tx.i;
	    i__1 = ij;
	    tx.r = xx[i__1].r, tx.i = xx[i__1].i;
	}
    }

/* ---  Find an element in the second half of the array which is */
/* ---  smaller than IT. */

L240:
    --l;
    if (ix[l] > it) {
	goto L240;
    }

/* ---  Find an element in the first half of the array which is */
/* ---  greater than IT. */

L245:
    ++k;
    if (ix[k] < it) {
	goto L245;
    }

/* ---  Interchange these elements. */

    if (k <= l) {
	iit = ix[l];
	ix[l] = ix[k];
	ix[k] = iit;
	i__1 = l;
	ttx.r = xx[i__1].r, ttx.i = xx[i__1].i;
	i__1 = l;
	i__2 = k;
	xx[i__1].r = xx[i__2].r, xx[i__1].i = xx[i__2].i;
	i__1 = k;
	xx[i__1].r = ttx.r, xx[i__1].i = ttx.i;
	goto L240;
    }

/* ---  Save upper and lower subscripts of the array yet to be sorted. */

    if (l - i__ > j - k) {
	il[m - 1] = i__;
	iu[m - 1] = l;
	i__ = k;
	++m;
    } else {
	il[m - 1] = k;
	iu[m - 1] = j;
	j = l;
	++m;
    }
    goto L260;

/* ---  Begin again on another portion of the unsorted array. */

L255:
    --m;
    if (m == 0) {
	goto L300;
    }
    i__ = il[m - 1];
    j = iu[m - 1];
L260:
    if (j - i__ >= 1) {
	goto L225;
    }
    if (i__ == j) {
	goto L255;
    }
    if (i__ == 1) {
	goto L210;
    }
    --i__;
L265:
    ++i__;
    if (i__ == j) {
	goto L255;
    }
    it = ix[i__ + 1];
    i__1 = i__ + 1;
    tx.r = xx[i__1].r, tx.i = xx[i__1].i;
    if (ix[i__] <= it) {
	goto L265;
    }
    k = i__;
L270:
    ix[k + 1] = ix[k];
    i__1 = k + 1;
    i__2 = k;
    xx[i__1].r = xx[i__2].r, xx[i__1].i = xx[i__2].i;
    --k;
    if (it < ix[k]) {
	goto L270;
    }
    ix[k + 1] = it;
    i__1 = k + 1;
    xx[i__1].r = tx.r, xx[i__1].i = tx.i;
    goto L265;
L300:
    return 0;
} /* izsrt1_ */

/* ----------------------------------------------------------------------| */
/* ----------------------------------------------------------------------| */
/* Subroutine */ int izsrt2_(integer *nx, integer *ix, integer *ixx, 
	doublecomplex *xx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal r__;
    static integer ij, il[21], it, iu[21], jt;
    static doublecomplex tx;
    static integer iit, jjt;
    static doublecomplex ttx;

/* ---  IZSRT2: indirect sort: sort ix and carry ixx and xx along */
/* ---  adapted from a SLAP (Sparse Linear Algebra Package) code. */
/* ----------------------------------------------------------------------| */
    /* Parameter adjustments */
    --xx;
    --ixx;
    --ix;

    /* Function Body */
    if (*nx <= 1) {
	return 0;
    }
/* ---  And now...Just a little black magic... */
    m = 1;
    i__ = 1;
    j = *nx;
    r__ = .375f;
L210:
    if (r__ <= .5898437f) {
	r__ += .0390625f;
    } else {
	r__ += -.21875f;
    }
L225:
    k = i__;

/* ---  Select a central element of the array and save it in location */
/* ---  IT, JT, TX. */

    ij = i__ + (integer) ((doublereal) (j - i__) * r__);
    it = ix[ij];
    jt = ixx[ij];
    i__1 = ij;
    tx.r = xx[i__1].r, tx.i = xx[i__1].i;

/* ---  If first element of array is greater than IT, interchange with IT. */

    if (ix[i__] > it) {
	ix[ij] = ix[i__];
	ix[i__] = it;
	it = ix[ij];
	ixx[ij] = ixx[i__];
	ixx[i__] = jt;
	jt = ixx[ij];
	i__1 = ij;
	i__2 = i__;
	xx[i__1].r = xx[i__2].r, xx[i__1].i = xx[i__2].i;
	i__1 = i__;
	xx[i__1].r = tx.r, xx[i__1].i = tx.i;
	i__1 = ij;
	tx.r = xx[i__1].r, tx.i = xx[i__1].i;
    }
    l = j;

/* ---  If last element of array is less than IT, swap with IT. */

    if (ix[j] < it) {
	ix[ij] = ix[j];
	ix[j] = it;
	it = ix[ij];
	ixx[ij] = ixx[j];
	ixx[j] = jt;
	jt = ixx[ij];
	i__1 = ij;
	i__2 = j;
	xx[i__1].r = xx[i__2].r, xx[i__1].i = xx[i__2].i;
	i__1 = j;
	xx[i__1].r = tx.r, xx[i__1].i = tx.i;
	i__1 = ij;
	tx.r = xx[i__1].r, tx.i = xx[i__1].i;

/* ---  If first element of array is greater than IT, swap with IT. */

	if (ix[i__] > it) {
	    ix[ij] = ix[i__];
	    ix[i__] = it;
	    it = ix[ij];
	    ixx[ij] = ixx[i__];
	    ixx[i__] = jt;
	    jt = ixx[ij];
	    i__1 = ij;
	    i__2 = i__;
	    xx[i__1].r = xx[i__2].r, xx[i__1].i = xx[i__2].i;
	    i__1 = i__;
	    xx[i__1].r = tx.r, xx[i__1].i = tx.i;
	    i__1 = ij;
	    tx.r = xx[i__1].r, tx.i = xx[i__1].i;
	}
    }

/* ---  Find an element in the second half of the array which is */
/* ---  smaller than IT. */

L240:
    --l;
    if (ix[l] > it) {
	goto L240;
    }

/* ---  Find an element in the first half of the array which is */
/* ---  greater than IT. */

L245:
    ++k;
    if (ix[k] < it) {
	goto L245;
    }

/* ---  Interchange these elements. */

    if (k <= l) {
	iit = ix[l];
	ix[l] = ix[k];
	ix[k] = iit;
	jjt = ixx[l];
	ixx[l] = ixx[k];
	ixx[k] = jjt;
	i__1 = l;
	ttx.r = xx[i__1].r, ttx.i = xx[i__1].i;
	i__1 = l;
	i__2 = k;
	xx[i__1].r = xx[i__2].r, xx[i__1].i = xx[i__2].i;
	i__1 = k;
	xx[i__1].r = ttx.r, xx[i__1].i = ttx.i;
	goto L240;
    }

/* ---  Save upper and lower subscripts of the array yet to be sorted. */

    if (l - i__ > j - k) {
	il[m - 1] = i__;
	iu[m - 1] = l;
	i__ = k;
	++m;
    } else {
	il[m - 1] = k;
	iu[m - 1] = j;
	j = l;
	++m;
    }
    goto L260;

/* ---  Begin again on another portion of the unsorted array. */

L255:
    --m;
    if (m == 0) {
	goto L300;
    }
    i__ = il[m - 1];
    j = iu[m - 1];
L260:
    if (j - i__ >= 1) {
	goto L225;
    }
    if (i__ == j) {
	goto L255;
    }
    if (i__ == 1) {
	goto L210;
    }
    --i__;
L265:
    ++i__;
    if (i__ == j) {
	goto L255;
    }
    it = ix[i__ + 1];
    jt = ixx[i__ + 1];
    i__1 = i__ + 1;
    tx.r = xx[i__1].r, tx.i = xx[i__1].i;
    if (ix[i__] <= it) {
	goto L265;
    }
    k = i__;
L270:
    ix[k + 1] = ix[k];
    ixx[k + 1] = ixx[k];
    i__1 = k + 1;
    i__2 = k;
    xx[i__1].r = xx[i__2].r, xx[i__1].i = xx[i__2].i;
    --k;
    if (it < ix[k]) {
	goto L270;
    }
    ix[k + 1] = it;
    ixx[k + 1] = jt;
    i__1 = k + 1;
    xx[i__1].r = tx.r, xx[i__1].i = tx.i;
    goto L265;
L300:
    return 0;
} /* izsrt2_ */

/* ----------------------------------------------------------------------| */
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
    static cilist io___91 = { 0, 7, 0, fmt_10, 0 };
    static cilist io___108 = { 0, 6, 0, 0, 0 };
    static cilist io___109 = { 0, 6, 0, 0, 0 };
    static cilist io___110 = { 0, 7, 0, fmt_11, 0 };
    static cilist io___113 = { 0, 6, 0, 0, 0 };
    static cilist io___114 = { 0, 7, 0, ptrfmt, 0 };
    static cilist io___115 = { 0, 7, 0, indfmt, 0 };
    static cilist io___116 = { 0, 7, 0, valfmt, 0 };
    static cilist io___119 = { 0, 6, 0, 0, 0 };


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
    s_rsfe(&io___91);
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
    s_wsle(&io___108);
    do_lio(&c__9, &c__1, title, (ftnlen)72);
    do_lio(&c__9, &c__1, "type :", (ftnlen)6);
    do_lio(&c__9, &c__1, type__, (ftnlen)3);
    do_lio(&c__9, &c__1, " size :", (ftnlen)7);
    do_lio(&c__3, &c__1, (char *)&nrow, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nn, (ftnlen)sizeof(integer));
    e_wsle();
    s_wsle(&io___109);
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
	s_rsfe(&io___110);
	do_fio(&c__1, rhstyp, (ftnlen)1);
	do_fio(&c__1, (char *)&nrhs, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nrhsix, (ftnlen)sizeof(integer));
	e_rsfe();
	s_wsle(&io___113);
	do_lio(&c__9, &c__1, "There is a second hand", (ftnlen)22);
	e_wsle();
    }
/* ---  read data... */
    s_rsfe(&io___114);
    i__1 = nn + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ja[i__], (ftnlen)sizeof(integer));
    }
    e_rsfe();
    s_rsfe(&io___115);
    i__1 = nnz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ia[i__], (ftnlen)sizeof(integer));
    }
    e_rsfe();
    s_rsfe(&io___116);
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
    s_wsle(&io___119);
    do_lio(&c__9, &c__1, "Harwell-Boeing matrix loaded", (ftnlen)28);
    e_wsle();
    return 0;
} /* loadhb_ */


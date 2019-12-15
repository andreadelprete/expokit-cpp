/* zgcnvr.f -- translated by f2c (version 20100827).
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


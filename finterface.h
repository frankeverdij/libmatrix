/*
 * finterface.h
 * Version: 1.2
 *
 * Copyright (C) 2014 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#ifndef F_INTERFACE_H
#define F_INTERFACE_H

#include "matrixtype.h"

/* Fortran routines in BLAS, LAPACK  and PROPACK *****************************/

#define F77_FCN(f) f##_

#ifdef ATLAS
    #include "cblas.h"
    #include "clapack.h"
#else
    int F77_FCN(dgetrf)(const M_INT * m, const M_INT * n, double * a, const M_INT * lda,
            M_INT * ipiv, M_INT * info);
    int F77_FCN(dgetri)(const M_INT *, double *, const M_INT *, M_INT *, double *, 
		    const M_INT *, M_INT *);
    int F77_FCN(dscal)(const M_INT *, const double *, double *, const M_INT *);
    int F77_FCN(daxpy)(const M_INT *, const double *, double *, const M_INT *, double *, const M_INT *);
    int F77_FCN(dcopy)(const M_INT *, double *, const M_INT *, double *, const M_INT *);
    double F77_FCN(ddot)(const M_INT *, const double *, const M_INT *, const double *, const M_INT *);
    /* void F77_FCN(drot)(const M_INT *N, double *DX, const M_INT *INCX, double *DY,
            const M_INT *INCY, double *C, double *S); */
    int F77_FCN(dgemv)(const char * trans, const M_INT *m, const M_INT *n,
            const double *alpha, const double *a, const M_INT *lda,
            const double *x, const M_INT *incx, const double *beta,
            double *y, const M_INT *incy);
    int F77_FCN(dgemm)(const char *transa, const char *transb, const M_INT *m,
            const M_INT *n, const M_INT *k, const double *alpha, const double *a,
            const M_INT *lda, const double *b, const M_INT *ldb, const double *beta,
            double *c, const M_INT *ldc);
    int F77_FCN(dgesv)(const M_INT *n, const M_INT *nrhs, double *a, const M_INT *lda,
            M_INT *ipiv, double *b, const M_INT *ldb, M_INT *info);
    /* int F77_FCN(dgetrs)(const char *, const M_INT *, const M_INT *, double *, 
		    const M_INT *, M_INT *, double *, const M_INT *, M_INT *); */
    /*int F77_FCN(dgeev)(const char *jobvl, const char *jobvr, const M_INT *n, 
		    double *A, const M_INT *lda, double *wr, double *wi, 
		    double *vl, const M_INT *ldvl, double *vr, const M_INT *ldvr, 
		    double *work, const M_INT *lwork, M_INT *info);*/
    /* int F77_FCN(dpotrf)(const char *uplo, const M_INT *n, double *A,
		    const M_INT *lda, M_INT *info); */
    /* int F77_FCN(dsvdc)(double *x, const M_INT *ldx, const M_INT *n, const M_INT *p,
            double *s, double *e, double *u, const M_INT *ldu, double *v,
            const M_INT *ldv, double *work, const M_INT *job, M_INT *info); */
    /* int F77_FCN(dlansvd)(const char *jobu, const char *jobv, const M_INT *m,
            const M_INT *n, const M_INT *k, const M_INT *kmax,void (*aprod)(const char *, M_INT *,
            M_INT *, double *, double *, double *, M_INT *), double *U, const M_INT *ldu,
            double *Sigma, double *bnd, double *V, const M_INT *ldv,const double *tolin,
            double *work, const M_INT *lwork, M_INT *iwork, const M_INT *liwork,
            const double *doption, const M_INT *ioption, M_INT *info, double *dparm, M_INT *iparm); */
#endif /* ATLAS */

int F77_FCN(dgesdd)(const char *jobz, const M_INT *m, const M_INT *n, double *A,
        const M_INT *lda, double *S, double *U, const M_INT *ldu, double *VT,
        const M_INT * ldvt, double *work, const M_INT *lwork, M_INT *iwork,
        M_INT *info);
int F77_FCN(dgesvd)(const char *jobu, const char *jobvt, const M_INT *m,
        const M_INT *n, double *A, const M_INT *lda, double *S, double *U,
        const M_INT *ldu, double *VT, const M_INT *ldvt, double *work,
        const M_INT *lwork, M_INT *info);
int F77_FCN(dgelss)(const M_INT *m, const M_INT *n, const M_INT *nrhs, double *a,
        const M_INT *lda, double *b, const M_INT *ldb, double *s,
        const double *rcond, M_INT * rank, double *work, const M_INT *lwork,
        M_INT *info);
int F77_FCN(dgels)(const char *trans, const M_INT *m, const M_INT *n,
        const M_INT *nrhs, double *a, const M_INT *lda, double *b,
        const M_INT *ldb, double *work, const M_INT *lwork, M_INT *info);

#endif /* F_INTERFACE_H */

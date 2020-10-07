/*
 * sparsenull.c
 * Version: 1.1
 *
 * Copyright (C) 2012 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */

#include "sparseutils.h"
#include "SuiteSparseQR_C.h"
#include "SuiteSparse_config.h"

Mat *Sparsenull(const Sparse *p, const double tolerance)
{
    int i,j,r ;
    UMFINT k,*E, *Ai, *Ap, *Qi, *Qp;
    double *Qx;
    cholmod_sparse *A,*Q,*R;
    cholmod_common Common, *cc;
    Mat *z;
    Sparse *q;
    
    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->nz>0);
    assert(p->t==DOUBLE);
    assert(p->cc);
    assert(tolerance>0.0);

    q=Sparsetranspose(p);

    /* start CHOLMOD */
    cc = &Common;
    cholmod_l_start(cc);

    A = (cholmod_sparse*)malloc(sizeof(cholmod_sparse));

    A->nrow=q->m;
    A->ncol=q->n;
    A->nzmax=q->nz;

    A->p=(void *)malloc((A->ncol+1)*sizeof(UMFINT));
    Ap=(UMFINT *)A->p;
    for (i=0;i<=A->ncol;i++) Ap[i]=q->cc[i];

    A->i=(void *)malloc(A->nzmax*sizeof(UMFINT));
    Ai=(UMFINT *)A->i;
    for (i=0;i<A->nzmax;i++) Ai[i]=q->r[i]-1;

    A->nz=NULL;
    A->x=(double *)q->s.d;
    A->z=NULL;

    A->stype=0;
    A->itype=CHOLMOD_LONG;
    A->xtype=CHOLMOD_REAL;
    A->dtype=CHOLMOD_DOUBLE;
    A->sorted=1;
    A->packed=1;

    /*A=cholmod_l_allocate_sparse(p->m,p->n,p->nz,1,1,0,CHOLMOD_REAL,cc);*/

    /*printf("A s %d  i %d x %d  d %d\n",A->stype,A->itype,A->xtype,A->dtype);*/

    k=MIN(q->m,q->n);
    
    r=SuiteSparseQR_C_QR(SPQR_ORDERING_DEFAULT,SPQR_DEFAULT_TOL,k,A,&Q,&R,&E,cc);
    /*printf("cc rank %ld\n",cc->SPQR_istat[4]);*/
    /*printf("Q s %d  i %d x %d  d %d\n",Q->stype,Q->itype,Q->xtype,Q->dtype);*/

    z=Matalloc(Q->nrow,Q->ncol-r,DOUBLE,1);
    Qp=(UMFINT *)Q->p;
    Qi=(UMFINT *)Q->i;
    Qx=(double *)Q->x;
    for (i=r;i<Q->ncol;i++) {
        for (j=Qp[i];j<Qp[i+1];j++) {
            z->m.d[Qi[j]][i-r]=(double)Qx[j];
        }
    }
    free(A->p);
    free(A->i);
    free(A);
    Sparsefree(q);

    cholmod_l_free_sparse(&Q, cc);
    cholmod_l_free_sparse(&R, cc);
    cholmod_l_free(p->m,sizeof(UMFINT),E, cc);

    /* end CHOLMOD */
    cholmod_l_finish(cc);

    return z;
}

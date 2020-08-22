/*
 * sparseutils.c
 * Version: 1.1
 *
 * Copyright (C) 2012 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#include "sparseutils.h"

Sparse *Sparsealloc(const Mat *a, const int type)
{
    int err;
    Sparse *p;

    p=malloc(sizeof(Sparse));
    assert(p);
    if (a==NULL) p->t = type;
        else p->t = a->t;
    p->m = 0;
    p->n = 0;
    p->nz = 0;
    p->r = NULL;
    p->c = NULL;
    p->cc = NULL;
    p->s.v = NULL;
    if (a!=NULL) err=Sparseassign(a,p);
    return p;
}

int Sparseassign(const Mat *a, Sparse *p)
{
    int err,i,j,k,m,n,nz,nz_guess;
    
    assert(a);    
    m=a->x;
    assert(m>=0);
    n=a->y;
    assert(n>=0);
    assert(p);    
    assert(a->t==p->t);
    
    Sparsenullifypointers(p);
    
    if (a->m.v==NULL) {
        p->m = m;
        p->n = n;
        p->nz = 0;
        p->r = NULL;
        p->c = NULL;
        p->cc = NULL;
        p->s.v = NULL;
    } else {
        assert(m>0);
        p->m = m;        
        assert(n>0);
        p->n = n;
        nz_guess=20*n;
        p->r = ialloc(nz_guess,0);
        p->c = ialloc(nz_guess,0);
        if (p->t==DOUBLE) p->s.d = dalloc(nz_guess,0);
            else p->s.i = ialloc(nz_guess,0);
        k=0;
        nz=0;
        for (i=0; i < m; i++) { 
            for (j=0; j < n; j++) {
                if (p->t==DOUBLE) {
                    if (fabs(a->m.d[i][j])>0.0) {
                        p->s.d[nz] = a->m.d[i][j];
                        p->r[nz] = i+1; /* use fortran convention for row/column indices! */
                        p->c[nz] = j+1; /* use fortran convention for row/column indices! */
                        /* resize value and column array */
                        if (++nz==nz_guess) nz_guess=Sparsetripletrealloc(p,2*nz_guess);
                    }
                } else {
                    if (abs(a->m.i[i][j])>0) {
                        p->s.i[nz] = a->m.i[i][j];
                        p->r[nz] = i+1; /* use fortran convention for row/column indices! */
                        p->c[nz] = j+1; /* use fortran convention for row/column indices! */
                        /* resize value and column array */
                        if (++nz==nz_guess) nz_guess=Sparsetripletrealloc(p,2*nz_guess);
                    }
                }
            }
        }
        /* resize value and column array to the correct size nz */
        p->nz = nz;
        if (p->nz) {
            err=Sparsetripletrealloc(p,p->nz);
            err=mendSparse(p,1);
        } else 
            Sparsenullifypointers(p);
    }
    
    return EXIT_SUCCESS;
}

double Sparseelement(const Sparse *p, const int r, const int c)
{
    int i;
    double s=0.0;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->nz>0);
    assert(p->t==DOUBLE);
    assert(p->r);
    assert(p->c);
    assert(p->s.v);
    if ((r>=1)&&(r<=p->m)&&(c>=1)&&(c<=p->n)) {
        for (i=0; i<p->nz; i++) {
            if ((p->r[i]==r)&&(p->c[i]==c)) s=p->s.d[i];
        }
    } else {
        fprintf(stderr,"Error in Sparseelement: r=%d c=%d p->m=%d p->n=%d\n",r,c,p->m,p->n);
        exit (EXIT_FAILURE);
    }
        
    return s;
}

int Sparseclearrow(Sparse *p, const int r)
{
    int i;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->nz>0);
    assert((p->t>=0)&&(p->t<=1));
    assert(p->r);
    assert(p->s.v);

    if ((r>=1)&&(r<=p->m)) {
        if (p->t==DOUBLE) {
            for (i=0; i<p->nz; i++) {
                if (p->r[i]==r) p->s.d[i]=0.0;
            }
        } else {
            for (i=0; i<p->nz; i++) {
                if (p->r[i]==r) p->s.i[i]=0;
            }
        }
    } else {
        fprintf(stderr,"Error in Sparseclearrow: r=%d p->m=%d\n",r,p->m);
        exit (EXIT_FAILURE);
    }
        
    return EXIT_SUCCESS;
}

int Sparsesetcolumn(Sparse *p, const int c)
{
    int i,err,off;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->nz>0);
    assert((p->t>=0)&&(p->t<=1));
    assert(p->r);
    assert(p->c);
    assert(p->s.v);

    /* make sure the sparse array is cleaned up and its compressed column
       array is set */
    if (p->cc==NULL) err=mendSparse(p,1);

    if ((c>=1)&&(c<=p->n)) {
        off=p->cc[c-1];
    
        if (p->cc[c]-off==0) {
            p->r=realloc(p->r,(p->nz+1)*sizeof(int));
            memmove(p->r+off+1,p->r+off,(p->nz-off)*sizeof(int));
            p->c=realloc(p->c,(p->nz+1)*sizeof(int));
            memmove(p->c+off+1,p->c+off,(p->nz-off)*sizeof(int));
    
            if (p->t==DOUBLE) {
                p->s.d=realloc(p->s.d,(p->nz+1)*sizeof(double));
                memmove(p->s.d+off+1,p->s.d+off,(p->nz-off)*sizeof(double));
            } else {
                p->s.i=realloc(p->s.i,(p->nz+1)*sizeof(int));
                memmove(p->s.i+off+1,p->s.i+off,(p->nz-off)*sizeof(int));    
            }
            for (i=c;i<=p->n;i++) p->cc[i]++;
            p->nz++;
        }

        p->r[off]=c;
        p->c[off]=c;
        if (p->t==DOUBLE) {
            p->s.d[off]=1.0;
            for (i=p->cc[c-1]+1;i<p->cc[c];i++) p->s.d[i]=0.0;
        } else {
            p->s.i[off]=1.0;
            for (i=p->cc[c-1]+1;i<p->cc[c];i++) p->s.i[i]=0;
        }

    } else {
        fprintf(stderr,"Error in Sparsesetcolumn: c=%d p->n=%d\n",c,p->n);
        exit (EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}

Mat *Sparsevector(const Sparse *p, const int x, const char rc)
{
    int i,j;
    Mat *a;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->nz>0);
    assert(p->t==DOUBLE);
    assert(p->r);
    assert(p->c);
    assert(p->s.v);

    if (rc=='r') j=p->m; else j=p->n;

    if ((x<1)||(x>j)) {
        fprintf(stderr,"Error in Sparsevector: %c x=%d max=%d\n",rc,x,j);
        exit (EXIT_FAILURE);
    }

    if (rc=='r') {
        a=Matalloc(1,p->n,DOUBLE,1);
        for (i=0; i<p->nz; i++) {
            if (p->r[i]==x) a->v.d[p->c[i]-1]+=p->s.d[i];
        }
    } else {
        a=Matalloc(p->m,1,DOUBLE,1);
        if (p->cc!=NULL) {
            for (i=p->cc[x-1];i<p->cc[x];i++) a->v.d[p->r[i]-1]+=p->s.d[i];              
        } else {
            for (i=0; i<p->nz; i++) {
                if (p->c[i]==x) a->v.d[p->r[i]-1]+=p->s.d[i];
            }
        }
    }
        
    return a;
}

int Sparseinsert(Sparse *p, const int r, const int c, const double s)
{
    int err;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->t==DOUBLE);

    if ((r>=1)&&(r<=p->m)&&(c>=1)&&(c<=p->n)) {
    
        err=Sparsetripletrealloc(p,p->nz+1);
        p->r[p->nz]=r;
        p->c[p->nz]=c;
        p->s.d[p->nz]=s;
        p->nz=err;
    } else {
        fprintf(stderr,"Error in Sparseinsert: r=%d c=%d p->m=%d p->n=%d\n",r,c,p->m,p->n);
        exit (EXIT_FAILURE);
    }


    /* invalidate CCS array */
    if (p->cc!=NULL) free(p->cc);
    p->cc=NULL;

    return EXIT_SUCCESS;
}

/* Collapses Sparse array row (and column) p[dof[i],dof[i]]=1.0 */
int SparsesetDirichlet(Sparse *p, int *dof, int ndof, const short c)
{
    int i,j,k,err;
    double *r;
    
    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->t==DOUBLE);

    if ((p->nz>0)&&(ndof>0)) {
        assert(p->r);
        assert(p->c);
        assert(p->s.v);

        /* make sure the sparse array is cleaned up and its compressed column
           array is set */
        if (p->cc==NULL) err=mendSparse(p,1);
    
        /* sort the dof's in the array */
        err=countsortint(dof,ndof);
    
        /* check if there are no duplicates. If there are, shift and adjust ndof */
        j=0;
        for (i=1;i<ndof;i++) {
            if (dof[i]!=dof[j]) j++;
            if (j<i) dof[j]=dof[i];
        }
        ndof=j+1;

        /* good, now let's first clear rows */
        /* initialize the reverse lookup array for the rows to be cleared */
        r=dalloc(p->m+1,0);
        for (i=0;i<=p->m;i++) r[i]=1.0;
        /* loop over all dofs we need to clear */
        for (i=0;i<ndof;i++) r[dof[i]] = 0.0;
        /* clear rows */
        for (i=0;i<p->nz;i++) p->s.d[i]*=r[p->r[i]];
        free(r);

        /* check if we need to clear columns */
        if (c==1) {
            /* allocate new space at the end to fit (dof,dof)=1.0 at the back
               of the sparse matrix. The mendSparse routine will clean this up */
            err=Sparsetripletrealloc(p,p->nz+ndof);
            /* loop over all dofs we need to clear */
            for (i=0;i<ndof;i++) {
                k=dof[i];
                /* clear all column entries for each dof CCS stylee */
                for (j=p->cc[k-1];j<p->cc[k];j++) p->s.d[j]=0.0;
                p->r[p->nz+i]=k;
                p->c[p->nz+i]=k;
                p->s.d[p->nz+i]=1.0;
            }
            /* don't forget to adjust to the new size */
            p->nz+=ndof;            
        }
    }    
    
    /* invalidate CCS array, since we made it dirty */
    if (p->cc!=NULL) free(p->cc);
    p->cc=NULL;

    return EXIT_SUCCESS;
}

int clearSparserow(Sparse *p, const int r)
{
    int i;
    
    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert((p->t>=0)&&(p->t<=1));

    if ((p->nz>0)&&(r>=1)&&(r<=p->m)) {
        assert(p->r);
        assert(p->c);
        assert(p->s.v);
        for (i=0;i<p->nz;i++) {
            if (p->r[i]==r) {
                if (p->t==DOUBLE) p->s.d[i]=0.0;
                    else p->s.i[i]=0;
            }
        }
    } else {
        fprintf(stderr,"Error in clearSparserow: r=%d p->m=%d\n",r,p->m);
        exit (EXIT_FAILURE);
    }

    /* invalidate CCS array */
    if (p->cc!=NULL) free(p->cc);
    p->cc=NULL;

    return EXIT_SUCCESS;
}

int clearSparsecolumn(Sparse *p, const int c)
{
    int i;
    
    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert((p->t>=0)&&(p->t<=1));

    if ((p->nz>0)&&(c>=1)&&(c<=p->n)) {
        assert(p->r);
        assert(p->c);
        assert(p->s.v);
        for (i=0;i<p->nz;i++) {
            if (p->c[i]==c) {
                if (p->t==DOUBLE) p->s.d[i]=0.0;
                    else p->s.i[i]=0;
            }
        }
    } else {
        fprintf(stderr,"Error in clearSparsecolumn: c=%d p->n=%d\n",c,p->n);
        exit (EXIT_FAILURE);
    }

    /* invalidate CCS array */
    if (p->cc!=NULL) free(p->cc);
    p->cc=NULL;

    return EXIT_SUCCESS;
}

Sparse *Sparsecopy(const Sparse *p)
{   
    Sparse *q;

    assert(p);
    assert((p->t>=0)&&(p->t<=1));

    q=Sparsealloc(NULL,p->t);
    q->m=p->m;
    q->n=p->n;
    q->nz=p->nz;

    if (q->nz) {
        q->r=ialloc(q->nz,0);
        memcpy(q->r,p->r,q->nz*sizeof(int));
        q->c=ialloc(q->nz,0);
        memcpy(q->c,p->c,q->nz*sizeof(int));
        if (q->t==DOUBLE) {
            q->s.d=dalloc(q->nz,0);
            memcpy(q->s.d,p->s.d,q->nz*sizeof(double));
        } else {
            q->s.i=ialloc(q->nz,0);
            memcpy(q->s.i,p->s.i,q->nz*sizeof(int));
        }
        if (p->cc!=NULL){
            q->cc=ialloc(q->n+1,0);
            memcpy(q->cc,p->cc,(q->n+1)*sizeof(int));
        }    
    }
    return q;
}

int *columntoCCS(int *c, const int n, const int nz)
{
    int i,*cc,check,dummy;

    cc=ialloc(n+1,1);
    check=0;
    for (i=0; i<nz; i++) {
        if (c[i]>=check) {
            cc[c[i]-1]++;
            check=c[i];
        } else {
            fprintf(stderr,"columntoCCS: Error! column array not sorted.\n");
            exit(0);
        }
    }
    /* construct the CCS array*/
    check=0;
    for (i=0;i<=n;i++) {
        dummy=cc[i];
        cc[i]=check;
        check+=dummy;
    }

    return cc;
}

int mendSparse(Sparse *p, const short addition)
{
    int i,err,tally;
 
    /* note, use only fortran style filled vectors since the routine
       cannot distinguish between them */
 
    assert(p);
    assert((p->t>=0)&&(p->t<=1));

    if (p->nz==0) return EXIT_SUCCESS;

    assert(p->m>0);
    assert(p->n>0);
    assert(p->r);
    assert(p->c);
    assert(p->s.v);

    /* clear CCS pointer */
    if (p->cc!=NULL) free(p->cc);

    /* first, remove all zeroes and check if row and column arrays have values
       beyond the current Sparse Matrix boundaries n and m */
    tally=-1;
    if (p->t==DOUBLE) {
        for (i=0; i<p->nz; i++) {
            if (fabs(p->s.d[i])>0.0) {
                assert(p->r[i]>=1);
                assert(p->c[i]>=1);
                p->m=MAX(p->m,p->r[i]);
                p->n=MAX(p->n,p->c[i]);
                if (++tally<i) {
                    p->r[tally]=p->r[i];
                    p->c[tally]=p->c[i];
                    p->s.d[tally]=p->s.d[i];
                }
            }           
        }
    } else {
        for (i=0; i<p->nz; i++) {
            if (abs(p->s.i[i])>0) {
                assert(p->r[i]>=1);
                assert(p->c[i]>=1);
                p->m=MAX(p->m,p->r[i]);
                p->n=MAX(p->n,p->c[i]);
                if (++tally<i) {
                    p->r[tally]=p->r[i];
                    p->c[tally]=p->c[i];
                    p->s.i[tally]=p->s.i[i];
                }
            }           
        }
    }
    p->nz=Sparsetripletrealloc(p,tally+1);

    /* then sort the triplets: column first, then row */
    err=countsortsparse(p->c,p->r,p->s.v,p->t,p->nz);

    /* next, find duplicate entries and add them together */
    tally=0;
    if (p->t==DOUBLE) {
        for (i=1; i<p->nz; i++) {
            /* duplicate entry? */
            if ((p->r[i]==p->r[tally]) && (p->c[i]==p->c[tally])) {
                /* yes, so add the values */
                p->s.d[tally]+=p->s.d[i];
           } else {
                /* nope, so shift the vectors, thereby deleting the duplicates
                   but only shift if the underlying entry is non-zero */
                if (fabs(p->s.d[tally])>0.0) tally++;
                if (tally<i) {
                    p->r[tally]=p->r[i];
                    p->c[tally]=p->c[i];
                    p->s.d[tally]=p->s.d[i];
                }
            }
        }
        /* correction for the last entry */
        if (fabs(p->s.d[tally])>0.0) tally++;
    } else {
        for (i=1; i<p->nz; i++) {
            if ((p->r[i]==p->r[tally]) && (p->c[i]==p->c[tally])) {
                p->s.i[tally]+=p->s.i[i];
            } else {
                if (abs(p->s.i[tally])>0) tally++;
                if (tally<i) {
                    p->r[tally]=p->r[i];
                    p->c[tally]=p->c[i];
                    p->s.i[tally]=p->s.i[i];
                }
            }
        }
        /* correction for the last entry */
        if (abs(p->s.i[tally])>0) tally++;    
    }

    p->nz=tally;

    /* are there any non-zero entries? if not, no point in continuing */
    if (p->nz) {        
        /* re-adjust size of p->r, p->c and p->s.v */        
        err=Sparsetripletrealloc(p,p->nz);

        /* construct the CCS array*/
        p->cc=columntoCCS(p->c, p->n, p->nz);
    } else
        Sparsenullifypointers(p);

    return EXIT_SUCCESS;
}

int SparseMatadd(Sparse *p, const Mat *rc, const Mat *s)
{
    int i;
    
    assert(rc);
    assert(rc->x==2);
    assert(rc->y>0);
    assert(rc->t==INTEGER);
    assert(rc->v.i);

    assert(s);
    assert(s->x>0);
    assert(s->y>0);
    assert(s->t==DOUBLE);
    assert(s->v.d);
    assert(rc->y==s->x*s->y);

    assert(p);
    assert(p->t==s->t);
   
    if (p->nz>0) {
        assert(p->m>0);
        assert(p->n>0);
        assert(p->r);
        assert(p->c);
        assert(p->s.v);

        p->r=realloc(p->r,(p->nz+rc->y)*sizeof(int));
        memcpy(&p->r[p->nz],rc->m.i[0],rc->y*sizeof(int));
        
        p->c=realloc(p->c,(p->nz+rc->y)*sizeof(int));
        memcpy(&p->c[p->nz],rc->m.i[1],rc->y*sizeof(int));

        if (p->t==DOUBLE) {
            p->s.d=realloc(p->s.d,(p->nz+rc->y)*sizeof(double));
            memcpy(&p->s.d[p->nz],s->v.d,s->x*sizeof(double));
        } else {
            p->s.i=realloc(p->s.i,(p->nz+rc->y)*sizeof(int));
            memcpy(&p->s.i[p->nz],s->v.i,s->x*sizeof(int));
        }

        p->nz+=rc->y;
    } else {
        assert(p->m==0);
        assert(p->n==0);
        assert(p->r==NULL);
        assert(p->c==NULL);
        assert(p->s.v==NULL);

        p->nz=rc->y;

        /* create row vector */
        p->r=ialloc(p->nz,0);
        p->m=rc->m.i[0][0];
        /* create column vector */
        p->c=ialloc(p->nz,0);
        p->n=rc->m.i[1][0];
        /* copy both row and column vector */
        for (i=0; i<p->nz; i++) {
            p->r[i]=rc->m.i[0][i]; /* fortran style */
            p->c[i]=rc->m.i[1][i]; /* fortran style */
            if (rc->m.i[0][i]>p->m) p->m=rc->m.i[0][i];
            if (rc->m.i[1][i]>p->n) p->n=rc->m.i[1][i];
        }
        /* create and copy the data vector */
        if (p->t==DOUBLE) {
            p->s.d=dalloc(p->nz,0);
            memcpy(p->s.d,s->v.d,p->nz*sizeof(double));
        } else {
            p->s.i=ialloc(p->nz,0);
            memcpy(p->s.i,s->v.i,p->nz*sizeof(int));
        }
    }

    /* invalidate CCS array */
    if (p->cc!=NULL) free(p->cc);
    p->cc=NULL;

    return EXIT_SUCCESS;
}

int ipSparsetrim(Sparse *p, const Mat *a, const char rc)
{
    int err,i,j,k,l;
    int *rev;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->t==DOUBLE);
    assert((p->t>=0)&&(p->t<=1));

    assert(a);
    assert(a->x==1);
    assert(a->y>0);
    assert(a->t==INTEGER);
    assert(a->m.v);

    l = (rc=='c') ? p->n : p->m;

    j=0;
    k=0;
    rev=ialloc(l,1);
    for (i=0;i<l;i++) {
        if (i==a->v.i[j]-1) {
            if (j<a->y-1) j++;
            k++;
            rev[i]=0;
        } else
            rev[i]=i+1-k;
    }
    
    if (rc=='c') p->n-=k; else p->m-=k;

    if (p->nz) {
        assert(p->r);
        assert(p->c);
        assert(p->s.v);
        if (rc=='c') {
            for (i=0;i<p->nz;i++) {
                j=rev[p->c[i]-1];
                if (j==0) {
                    p->s.d[i]=0.0;
                    p->c[i]=1;
                } else
                    p->c[i]=j;
            }
        } else { 
            for (i=0;i<p->nz;i++) {
                j=rev[p->r[i]-1];
                if (j==0) {
                    p->s.d[i]=0.0;
                    p->r[i]=1;
                } else
                    p->r[i]=j;
            }
        }
 
        if (p->cc!=NULL) free(p->cc);
        p->cc=NULL;
    
        err=mendSparse(p,1);
    }

    free(rev);

    return EXIT_SUCCESS;
}

int ipSparsesplice(Sparse *p, const Mat *a, const char rc)
{
    int err,i,j,k,l;
    int *rev;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert((p->t>=0)&&(p->t<=1));

    assert(a);
    assert(a->x==1);
    assert(a->y>0);
    assert(a->t==INTEGER);
    assert(a->m.v);

    if (p->nz) {
        assert(p->r);
        assert(p->c);
        assert(p->s.v);
    
        l = (rc=='c') ? p->n : p->m;
        rev=ialloc(l,1);
    
        i=0;
        j=0;
        k=0;
        for (i=0;i<l;i++) {
            while (i+k==a->v.i[j]-1) {
                if (j<a->y-1) j++;
                k++;
            }
            rev[i]=i+k+1;
        }
    
        if (rc=='c') {
            p->n+=a->y;
            for (i=0;i<p->nz;i++) {
                j=rev[p->c[i]-1];
                p->c[i]=j;
            }
        } else {
            p->m+=a->y;
            for (i=0;i<p->nz;i++) {
                j=rev[p->r[i]-1];
                p->r[i]=j;
            }
        }
        free(rev);
    
        if (p->cc!=NULL) free(p->cc);
        p->cc=NULL;

        err=mendSparse(p,1);
    }

    return EXIT_SUCCESS;
}

/* B = P * A */
Mat *SparseMatmul(const Sparse *p, const Mat *a)
{
    int i,j,k,r;
    double ad;
    Mat *b;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->t==DOUBLE);

    assert(a);
    assert(p->n==a->x);
    assert(a->y>0);
    assert(a->t==p->t);
    assert(a->m.v);

    b=Matalloc(p->m,a->y,DOUBLE,1);

    if (p->nz) {
        assert(p->r);
        assert(p->c);
        assert(p->s.v);
        /* Check if sparse matrix p is proper ccs format */
        assert(p->cc);
        /* special case for dense vectors, big speedup */
        if (a->y==1) {
            for (i=0;i<p->n;i++) {
                ad=a->v.d[i];
                for (j=p->cc[i];j<p->cc[i+1];j++) {
                    b->v.d[p->r[j]-1]+=p->s.d[j]*ad;
                }
            }
        } else {
            for (i=0;i<p->n;i++) {
                for (j=p->cc[i];j<p->cc[i+1];j++) {
                    r=p->r[j]-1;
                    for (k=0;k<a->y;k++) b->m.d[r][k]+=p->s.d[j]*a->m.d[i][k];
                }
            }
        }
    }

    return b;
}

/* B = P^T * A */
Mat *SparsetransposeMatmul(const Sparse *p, const Mat *a)
{
    int i,j,k,r;
    double bd;
    Mat *b;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->t==DOUBLE);

    assert(a);
    assert(a->x==p->m);
    assert(a->y>0);
    assert(a->t==p->t);
    assert(a->m.v);

    b=Matalloc(p->n,a->y,DOUBLE,1);

    if (p->nz) {
        assert(p->r);
        assert(p->c);
        assert(p->s.v);
        /* Check if sparse matrix p is proper ccs format */
        assert(p->cc);
        /* special case for dense vectors, big speedup */
        if (a->y==1) {
            for (i=0;i<p->n;i++) {
                bd=0.0;
                for (j=p->cc[i];j<p->cc[i+1];j++) {
                    bd+=p->s.d[j]*a->v.d[p->r[j]-1];
                }
                b->v.d[i]=bd;
            }
        } else {
            for (i=0;i<p->n;i++) {
                for (j=p->cc[i];j<p->cc[i+1];j++) {
                    r=p->r[j]-1;
                    for (k=0;k<a->y;k++) b->m.d[i][k]+=p->s.d[j]*a->m.d[r][k];
                }
            }
        }
    }

    return b;
}

Sparse *Sparsemul(const Sparse *p, const Sparse *q)
{
    Sparse *x;
    int err,i,j,k,l,nz_guess,tally;
    double qs;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->t==DOUBLE);

    assert(q);
    assert(q->m>0);
    assert(q->n>0);
    assert(q->t==p->t);

    if (p->n!=q->m) {
        fprintf(stderr,"Error in Sparsemul: p->n:%d q->m:%d\n",p->n,q->m);
        exit (EXIT_FAILURE);
    }

    x=Sparsealloc(NULL,p->t);
    x->m=p->m;
    x->n=q->n;
    
    if ((p->nz==0)||(q->nz==0)) return x;

    assert(p->r);
    assert(p->c);
    assert(p->s.v);
    /* Check if sparse matrix p is proper ccs format */
    assert(p->cc);

    assert(q->r);
    assert(q->c);
    assert(q->s.v);
    /* Check if sparse matrix q is proper ccs format */
    assert(q->cc);

    nz_guess=p->nz;
    tally=0;
    x->r=ialloc(nz_guess,0);
    x->c=ialloc(nz_guess,0);
    x->s.d=dalloc(nz_guess,0);
    x->cc=NULL;

    for (i=0;i<q->n;i++) {
        for (j=q->cc[i];j<q->cc[i+1];j++) {
            k=q->r[j];
            qs=q->s.d[j];
            for (l=p->cc[k-1];l<p->cc[k];l++) {
                x->s.d[tally]=p->s.d[l]*qs;
                x->r[tally]=p->r[l];
                x->c[tally]=i+1;
                if (++tally==nz_guess) nz_guess=Sparsetripletrealloc(x,1.5*nz_guess);
            }
        }
    }
    x->nz=Sparsetripletrealloc(x,tally);
    err=mendSparse(x,1);
    
    return x;
}

void Sparsefree(Sparse *p)
{
    if (p != NULL) {
        Sparsenullifypointers(p);
        free(p);
    }
    p = NULL;

    return;
}

Sparse **PSparsealloc(const int n)
{
    int i;
    Sparse **p;

    p = malloc(n*sizeof(Sparse *));
    assert(p);
    for (i=0;i<n;i++) p[i] = NULL;

    return p;
}

void PSparsefree(Sparse **p, const int n)
{
    int i;

    if (p != NULL) {
        for (i=0;i<n;i++) Sparsefree(p[i]);
        free(p);
    }
    p = NULL;

    return;
}

int Sparsedump(const Sparse *p)
{
    int i;

    assert(p);
    assert(p->m);
    assert(p->n);
    if (p->cc!=NULL) {
    printf("cc ");
        for (i=0;i<p->n+1;i++) printf("%d ",p->cc[i]);
    }
    printf("\n");

    if (p->t==DOUBLE) {
        if (p->nz) for (i=0;i<p->nz;i++) printf("%d %d %22.18e\n",p->r[i],p->c[i],p->s.d[i]);
        printf("Sparse Matrixtype:Double x:%d y:%d nonzeros:%d\n",p->m, p->n, p->nz);
    } else {
        if (p->nz) for (i=0;i<p->nz;i++) printf("%d %d %d\n",p->r[i],p->c[i],p->s.i[i]);
        printf("Sparse Matrixtype:Integer x:%d y:%d nonzeros:%d\n",p->m, p->n, p->nz);
    }

    return EXIT_SUCCESS;
}

int Sparsetripletrealloc(Sparse *p, const int n)
{
    assert(p);
    assert(p->r);
    assert(p->c);
    assert((p->t>=0)&&(p->t<=1));
    assert(p->s.v);

    if (n) {
        p->r = realloc(p->r,n*sizeof(int));
        assert(p->r);
        p->c = realloc(p->c,n*sizeof(int));
        assert(p->c);
        if (p->t==DOUBLE) p->s.d = realloc(p->s.d,n*sizeof(double));
            else p->s.i = realloc(p->s.i,n*sizeof(int));
        assert(p->s.v);
    } else Sparsenullifypointers(p);

    return n;
}

void Sparsenullifypointers(Sparse *p)
{
    if (p->r!=NULL) free(p->r);
    if (p->c!=NULL) free(p->c);
    if (p->cc!=NULL) free(p->cc);
    if (p->s.v!=NULL) free(p->s.v);
    p->r=NULL;
    p->c=NULL;
    p->cc=NULL;
    p->s.v=NULL;

    return;
}

Mat *SparseMatdiv(const Sparse *p, const Mat *a)
{
    void *Numeric=NULL;
    
    return SparseMatdivfull(p, a, &Numeric, (short) 0);
}

/* SparseMatdivfull
 * computes x=p\q where p is sparse and x and q are Mat objects
 *
 * to take advantage of numeric factorization of p for repeated solves, the
 * Numeric pointer and the fact flag produce the following behavior:
 *
 *                  Numeric
 *           ==NULL        !=NULL
 * fact 0    create&free   free
 *      1    create        don't change
 *
 */
Mat *SparseMatdivfull(const Sparse *p, const Mat *a, void **Numeric, const short fact)
{
    int i,j,err;
    UMFINT *Ai, *Acc;
    void *Symbolic;
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
    double *u,*f;
    Mat *x;

    assert(p);
    assert(p->m>0);
    assert(p->n==p->m);
    assert(p->nz>0);
    assert(p->t==DOUBLE);
    assert(p->s.v);
    /* Check if sparse matrix p is in proper ccs format */
    assert(p->cc);

    assert((fact>=0)&&(fact<=1));
    assert(a);
    assert(a->x==p->m);
    assert(a->y>0);
    assert(a->t==p->t);

    x=Matalloc(p->n,a->y,DOUBLE,1);
    u=dalloc(x->x,1);
    f=dalloc(a->x,1);

    Ai=(UMFINT *)malloc(p->nz*sizeof(UMFINT));
    for (i=0; i<p->nz; i++) Ai[i]=(UMFINT)p->r[i]-1;

    Acc=(UMFINT *)malloc((p->n+1)*sizeof(UMFINT));
    for (i=0;i<=p->n;i++) Acc[i]=(UMFINT)p->cc[i];

    UMFGLUE(defaults)(Control);
    /*Control[UMFPACK_SCALE]=UMFPACK_SCALE_SUM;*/

    if (*Numeric==NULL) {
        err=UMFGLUE(symbolic)(p->m, p->n, Acc, Ai, p->s.d, &Symbolic, Control, Info);
        if (err<0) {
            printf("Error %d in umfpack_dx_symbolic\n",err);
            UMFGLUE(report_info)(Control, Info);
    	    UMFGLUE(report_status)(Control, err);
        }
        err=UMFGLUE(numeric)(Acc, Ai, p->s.d, Symbolic, Numeric, Control, Info);
        if (err<0) {
            printf("Error %d in umfpack_dx_numeric\n",err);
            UMFGLUE(report_info)(Control, Info);
    	    UMFGLUE(report_status)(Control, err);
        }
        UMFGLUE(free_symbolic)(&Symbolic);
    }
    for (j=0;j<a->y;j++) {
        for (i=0;i<a->x;i++) f[i]=a->m.d[i][j]; /* extract column */

        err=UMFGLUE(solve)(UMFPACK_A, Acc, Ai, p->s.d, u, f, *Numeric, Control, Info);
        if (err<0) {
            printf("Error %d in umfpack_dx_solve\n",err);
            UMFGLUE(report_info)(Control, Info);
    	    UMFGLUE(report_status)(Control, err);
        }

        for (i=0;i<x->x;i++) x->m.d[i][j]=u[i];
    }
    if (fact==0) {
        UMFGLUE(free_numeric)(Numeric);
        *Numeric=NULL;
    }
    free(Ai);
    free(Acc);
    free(f);
    free(u);

    return x;
}

Sparse *Sparsediv(const Sparse *p, const Sparse *q)
{
    void *Numeric=NULL;
    
    return Sparsedivfull(p, q, &Numeric, (short) 0);
}

/* Sparsedivfull
 * computes x=p\q where x p and q are all sparse
 *
 * to take advantage of numeric factorization of p for repeated solves, the
 * Numeric pointer and the fact flag produce the following behavior:
 *
 *                  Numeric
 *           ==NULL        !=NULL
 * fact 0    create&free   free
 *      1    create        don't change
 *
 */
Sparse *Sparsedivfull(const Sparse *p, const Sparse *q, void **Numeric, const short fact)
{
    int i,j,k,err,nz_guess,tally;
    UMFINT *Ai,*Acc;
    void *Symbolic;
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO], *u, *f;
    Sparse *x;

    assert(p);
    assert(p->m>0);
    assert(p->n==p->m);
    assert(p->nz>0);
    assert(p->t==DOUBLE);
    assert(p->s.v);
    /* Check if sparse matrix p is in proper ccs format */
    assert(p->cc);

    assert(q);
    assert(q->m==p->m);
    assert(q->n>0);
    assert(q->nz>0);
    assert(q->t==p->t);
    assert(q->s.v);
    /* Check if sparse matrix q is in proper ccs format */
    assert(q->cc);

    assert((fact>=0)&&(fact<=1));

    x=Sparsealloc(NULL,DOUBLE);
    x->m=p->n;
    x->n=q->n;
    nz_guess=p->nz;
    x->r=ialloc(nz_guess,1);
    x->c=ialloc(nz_guess,1);
    x->s.d=dalloc(nz_guess,1);

    Ai=(UMFINT *)malloc(p->nz*sizeof(UMFINT));
    for (i=0; i<p->nz; i++) Ai[i]=(UMFINT)p->r[i]-1;

    Acc=(UMFINT *)malloc((p->n+1)*sizeof(UMFINT));
    for (i=0;i<=p->n;i++) Acc[i]=(UMFINT)p->cc[i];

    UMFGLUE(defaults)(Control);
    if (*Numeric==NULL) {
        err=UMFGLUE(symbolic)(p->m, p->n, Acc, Ai, p->s.d, &Symbolic, Control, Info);
        if (err<0) {
            printf("Error %d in umfpack_dx_symbolic\n",err);
            UMFGLUE(report_info)(Control, Info);
    	    UMFGLUE(report_status)(Control, err);
        }
        err=UMFGLUE(numeric)(Acc, Ai, p->s.d, Symbolic, Numeric, Control, Info);
        if (err<0) {
            printf("Error %d in umfpack_dx_numeric\n",err);
            UMFGLUE(report_info)(Control, Info);
    	    UMFGLUE(report_status)(Control, err);
        }
        UMFGLUE(free_symbolic)(&Symbolic);
    }

    u=dalloc(x->m,1);
    f=dalloc(x->m,1);
    tally=0;
    for (j=0;j<x->n;j++) {
        if (q->cc[j]!=q->cc[j+1]) { /* column of q is empty? skip it */
            for (k=q->cc[j];k<q->cc[j+1];k++)
                f[q->r[k]-1]=q->s.d[k]; /* extract column */

            err=UMFGLUE(solve)(UMFPACK_A, Acc, Ai, p->s.d, u, f, *Numeric, Control, Info);
            if (err<0) {
                printf("Error %d in umfpack_dx_solve\n",err);
                UMFGLUE(report_info)(Control, Info);
        	    UMFGLUE(report_status)(Control, err);
            }
            memset(f,0,x->m*sizeof(double));
            /* for (k=q->cc[j];k<q->cc[j+1];k++)
                f[q->r[k]-1]=0.0; */

            for (i=0;i<x->m;i++) {
                if (fabs(u[i])>0.0) {
                    x->r[tally]=i+1;
                    x->c[tally]=j+1;
                    x->s.d[tally]=u[i];
                    if (++tally==nz_guess) nz_guess=Sparsetripletrealloc(x,2*nz_guess);
                }
            }
        }
    }
    x->nz=tally;
    free(Ai);
    free(Acc);
    free(u);
    free(f);

    /* are there any non-zero entries? */
    if (x->nz) {
        /* re-adjust size of x->r, x->c and x->s.d */
        err=Sparsetripletrealloc(x,x->nz);
        x->cc=columntoCCS(x->c, x->n, x->nz);
    } else Sparsenullifypointers(x);

    if (fact==0) {
        UMFGLUE(free_numeric)(Numeric);
        *Numeric=NULL;
    }

    return x;
}

int ipSparseadd(Sparse *p, const Sparse *q, const double scal, const short domend)
{
    int i,err;

    assert(p);
    assert(p->t==DOUBLE);
    assert(q);
    assert(q->t==p->t);

    if (q->nz==0) return EXIT_SUCCESS;

    if (p->nz) {
        assert(q->m==p->m);
        assert(q->n==p->n);
        err=Sparsetripletrealloc(p,p->nz+q->nz);
    } else {
        Sparsenullifypointers(p);
        p->m=q->m;
        p->n=q->n;
        p->r=ialloc(q->nz,0);
        p->c=ialloc(q->nz,0);
        if (p->t==DOUBLE) p->s.d=dalloc(q->nz,0);
            else p->s.i=ialloc(q->nz,0);
    }

    memcpy(p->r+(p->nz),q->r,q->nz*sizeof(int));
    memcpy(p->c+(p->nz),q->c,q->nz*sizeof(int));
    for (i=0;i<q->nz;i++) p->s.d[i+p->nz]=q->s.d[i]*scal;
    p->nz+=q->nz;

     /* invalidate CCS array */
    if (p->cc!=NULL) free(p->cc);
    p->cc=NULL;
    if (domend) err=mendSparse(p,1);
  
    return EXIT_SUCCESS;
}

int ipSparsetranspose(Sparse *p)
{
    int err,*iswap;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert((p->t>=0)&&(p->t<=1));

    if (p->nz) {
        err=p->m;
        p->m=p->n;
        p->n=err;
    
        iswap=p->r;
        p->r=p->c;
        p->c=iswap;
    
        /* invalidate CCS array */
        if (p->cc!=NULL) free(p->cc);
        p->cc=NULL;
    
        err=mendSparse(p,1);
    }

    return EXIT_SUCCESS;
}

Sparse *Sparsetranspose(const Sparse *p)
{
    int err;
    Sparse *q;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert((p->t>=0)&&(p->t<=1));

    q=Sparsealloc(NULL,p->t);
    q->m=p->n;
    q->n=p->m;
    q->nz=p->nz;

    if (q->nz) {
        q->r=ialloc(q->nz,0);
        memcpy(q->r,p->c,q->nz*sizeof(int));
        q->c=ialloc(q->nz,0);
        memcpy(q->c,p->r,q->nz*sizeof(int));
        if (q->t==DOUBLE) {
            q->s.d=dalloc(q->nz,0);
            memcpy(q->s.d,p->s.d,q->nz*sizeof(double));
        } else {
            q->s.i=ialloc(q->nz,0);
            memcpy(q->s.i,p->s.i,q->nz*sizeof(int));
        }
        err=mendSparse(q,1);
    }

    return q;
}

Mat *SparsetoMat(const Sparse *p)
{
    int i,j;
    Mat *a;

    if (p==NULL) return NULL;

    if ((p->nz==0)&&((p->m==0)||(p->n==0)))
        return NULL;

    a=Matalloc(p->m,p->n,p->t,1);

    if (p->nz) {
        /* Check if sparse matrix p is proper ccs format */
        if (p->cc) {
            if (a->t==DOUBLE) {
                for (i=0;i<p->n;i++) {
                    for (j=p->cc[i];j<p->cc[i+1];j++) {
                        a->m.d[p->r[j]-1][i]=p->s.d[j];
                    }
                }           
            } else {
                for (i=0;i<p->n;i++) {
                    for (j=p->cc[i];j<p->cc[i+1];j++) {
                        a->m.i[p->r[j]-1][i]=p->s.i[j];
                    }
                }           
            }
        } else {
            for (i=0;i<p->nz;i++) {
                if (a->t==DOUBLE) a->m.d[p->r[i]-1][p->c[i]-1]+=p->s.d[i];
                    else a->m.i[p->r[i]-1][p->c[i]-1]+=p->s.i[i];
            }
        }
    }
    return a;
}

/*
 * ipSparsecat concatenates Sparse Matrix q to p
 *
 * if rc = 'r' concatenate by rows
 * if rc = 'c' concatenate by column
 * if rc = 'd' concatenate diagonally
 *
 * if domend != 0 call mendSparse to repair the Sparse matrix
 */
int ipSparsecat(Sparse *p, const Sparse *q, const char rc, const short domend)
{
    int i,err;

    assert(p);
    assert((p->t>=0)&&(p->t<=1));
    assert(q);
    assert(q->t==p->t);

    if (q->nz==0) return EXIT_SUCCESS;

    if (p->nz) {
        err=Sparsetripletrealloc(p,p->nz+q->nz);
    } else {
        p->m=0;
        p->n=0;
        Sparsenullifypointers(p);
        p->r=ialloc(q->nz,0);
        p->c=ialloc(q->nz,0);
        if (p->t==DOUBLE) p->s.d=dalloc(q->nz,0);
            else p->s.i=ialloc(q->nz,0);
    }

    if (p->t==DOUBLE) memcpy(p->s.d+(p->nz),q->s.d,q->nz*sizeof(double));
        else memcpy(p->s.i+(p->nz),q->s.i,q->nz*sizeof(int));

    if (rc=='c') {
        memcpy(p->r+(p->nz),q->r,q->nz*sizeof(int));
        p->m=MAX(p->m,q->m);
    }
    if (rc=='r') {
        memcpy(p->c+(p->nz),q->c,q->nz*sizeof(int));
        p->n=MAX(p->n,q->n);
    }

    if ((rc=='c')||(rc=='d')) {
        for (i=0;i<q->nz;i++) p->c[p->nz+i]=q->c[i]+p->n;    
        p->n+=q->n;
    }
    if ((rc=='r')||(rc=='d')) {
        for (i=0;i<q->nz;i++) p->r[p->nz+i]=q->r[i]+p->m;
        p->m+=q->m;
    }

    p->nz+=q->nz;        
    
    if (p->cc!=NULL) free(p->cc);
    p->cc=NULL;        
    if (domend) err=mendSparse(p,1);

    return EXIT_SUCCESS;
}

Sparse *Sparseread(const char *filename)
{
    int err,i,j,m,n;
    Sparse *p;
    FILE *f;

    f=fopen(filename,"r");
    err=fscanf(f,"%d %d",&m,&n);
    assert(m>0);
    assert(n>0);

    p=Sparsealloc(NULL,DOUBLE);
    p->m=m;
    p->n=n;

    p->cc=ialloc(n+1,0);
    for (i=0;i<=n;i++) err=fscanf(f,"%d",&(p->cc[i]));
    p->nz=p->cc[n];

    p->r=ialloc(p->nz,0);
    p->c=ialloc(p->nz,0);
    p->s.d=dalloc(p->nz,0);
    for (i=0;i<p->nz;i++) err=fscanf(f,"%d %lf",&(p->r[i]),&(p->s.d[i]));
    for (i=0;i<n;i++) {
        for (j=p->cc[i];j<p->cc[i+1];j++) p->c[j]=i+1;
    }
    fclose(f);

    return p;
}

int Sparsewrite(const char *filename, const Sparse *p, const int append)
{
    int i;
    FILE *f;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert((p->t>=0)&&(p->t<=1));
    assert(p->r);
    assert(p->c);
    assert(p->s.v);
    assert(p->cc);

    if (append) f=fopen(filename,"a+");        
        else f=fopen(filename,"w");

    fprintf(f, "%d %d\n",p->m,p->n);
    for (i=0; i<=p->n; i++) fprintf(f,"%d\n",p->cc[i]);
    if (p->t==DOUBLE)
        for (i=0; i<p->nz; i++) fprintf(f,"%d %22.18e\n",p->r[i],p->s.d[i]);
    else
        for (i=0; i<p->nz; i++) fprintf(f,"%d %d",p->r[i],p->s.i[i]);

    fclose(f);

    return EXIT_SUCCESS;
}

/*
 * Sparseread_mtx(filename)
 * reads a Matrix Market file into a Sparse structure
 *
 * modified from http://math.nist.gov/MatrixMarket/mmio/c/example_read.c
 */
Sparse *Sparseread_mtx(const char *filename)
{
    Sparse *p;
    int i,err,m,n,nz;
    FILE *fp;
    MM_typecode matcode;

    fp=fopen(filename,"r");

    /* read the %% banner */
    if (mm_read_banner(fp, &matcode) != 0)
    {
        fprintf(stderr,"Could not process Matrix Market banner.\n");
        exit (EXIT_FAILURE);
    }

    /* exit if the matrix is complex or dense */
    if (mm_is_complex(matcode) || mm_is_dense(matcode)) {
        fprintf(stderr,"Sparseread_mtx() does not support ");
        fprintf(stderr,"Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(EXIT_FAILURE);
    }

    /* read size and dimensions of the sparse matrix */
    if ((err = mm_read_mtx_crd_size(fp, &m, &n, &nz)) !=0) {
        fprintf(stderr,"Could not determine Sparse matrix dimensions.");
        exit(EXIT_FAILURE);
    }

    p=Sparsealloc(NULL,DOUBLE);
    p->m=m;
    p->n=n;
    p->nz=nz;
    p->r=ialloc(nz,0);
    p->c=ialloc(nz,0);

    if (mm_is_real(matcode)) {
        p->t=DOUBLE;
        p->s.d=dalloc(nz,0);
    } else if mm_is_integer(matcode) {
        p->t=INTEGER;
        p->s.i=ialloc(nz,0);
    } else {
        fprintf(stderr,"Sparseread_mtx() does not support ");
        fprintf(stderr,"Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(EXIT_FAILURE);
    }
    
    for (i=0;i<nz;i++) {
        if (p->t==DOUBLE)
            err=fscanf(fp,"%d %d %lg",&(p->r[i]),&(p->c[i]),&(p->s.d[i]));
        else
            err=fscanf(fp,"%d %d %d",&(p->r[i]),&(p->c[i]),&(p->s.i[i]));
    }
    fclose(fp);

    mendSparse(p,1);

    return p;
}

/*
 * Sparsewrite_mtx(filename,p)
 * Writes Sparse Matrix p into a Matrix Market file
 *
 * modified from http://math.nist.gov/MatrixMarket/mmio/c/example_write.c
 */
void Sparsewrite_mtx(const char *filename, const Sparse *p)
{
    int i,j;
    FILE *fp;
    MM_typecode matcode;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert((p->t>=0)&&(p->t<=1));
    assert(p->r);
    assert(p->c);
    assert(p->s.v);
    assert(p->cc);

    fp=fopen(filename,"w");

    mm_initialize_typecode(&matcode);
    mm_set_sparse(&matcode);
    mm_set_coordinate(&matcode);
    if (p->t==DOUBLE)
        mm_set_real(&matcode);
    else
        mm_set_integer(&matcode);

    mm_write_banner(fp, matcode); 
    mm_write_mtx_crd_size(fp, p->m, p->n, p->nz);

    for (i=0;i<p->n;i++) {
        if (p->t==DOUBLE) {
            for (j=p->cc[i];j<p->cc[i+1];j++) {
                fprintf(fp,"%d %d %22.18e\n",p->r[j],i+1,p->s.d[j]);
            }
        } else {
            for (j=p->cc[i];j<p->cc[i+1];j++) {
                fprintf(fp,"%d %d %d\n",p->r[j],i+1,p->s.i[j]);
            }
        }
    }
    fclose(fp);

    return;
}

Sparse *Sparseeye(const int m, const int type)
{
    int i;
    Sparse *p;

    assert(m>0);
    assert((type>=0)&&(type<=1));
    
    p=Sparsealloc(NULL,DOUBLE);
    p->m=m;
    p->n=m;
    p->nz=m;
    p->t=type;
    p->r=ialloc(m,0);
    p->c=ialloc(m,0);
    p->cc=ialloc(m+1,1);
    if (type==DOUBLE) p->s.d=dalloc(m,0);
        else p->s.i=ialloc(m,0);

    for (i=0;i<m;i++) {
        p->r[i]=i+1;
        p->c[i]=i+1;
        p->cc[i+1]=i+1;
        if (type==DOUBLE) {
            p->s.d[i]=1.0;
        } else {
            p->s.i[i]=1;
        }
    }

    return p;
}

void ipSparsediag(Sparse *p)
{
    int i;

    assert(p);
    assert(p->m>0);
    assert(p->n==p->m);
    assert(p->r);
    assert(p->c);
    assert(p->s.v);
    assert(p->cc);

    if (p->m!=p->n) {
        fprintf(stderr,"ipSparsediag: error: Sparse Matrix should be square but is %d x %d\n",p->m,p->n);
        exit(EXIT_FAILURE);
    }

    if (p->t == DOUBLE) {
        for (i=0;i<p->nz;i++) if (p->r[i] != p->c[i]) p->s.d[i]=0.0;
    } else {
        for (i=0;i<p->nz;i++) if (p->r[i] != p->c[i]) p->s.i[i]=0;
    }
    if (p->cc!=NULL) {
        free(p->cc);
        p->cc=NULL;
    }
    mendSparse(p,1);

    return;
}

Sparse *Sparsediag(const Sparse *p)
{
    Sparse *q;

    assert(p);
    assert(p->m>0);
    assert(p->n==p->m);
    assert(p->r);
    assert(p->c);
    assert(p->s.v);
    assert(p->cc);
    
    if (p->m!=p->n) {
        fprintf(stderr,"Sparsediag: error: Sparse Matrix should be square but is %d x %d\n",p->m,p->n);
        exit(EXIT_FAILURE);
    }
    
    q=Sparsecopy(p);
    ipSparsediag(q);
    
    return q;
}

int ipSparsesymgraph(Sparse *p)
{
    int err,flag,i,j,k,m,nz=0,nz_guess;
    Sparse *q;

    assert(p);
    assert(p->m>0);
    assert(p->n==p->m);
    assert(p->r);
    assert(p->c);
    assert(p->s.v);
    assert(p->cc);

    if (p->m!=p->n) {
        fprintf(stderr,"ipSparsesymgraph: error: Sparse Matrix should be square but is %d x %d\n",p->m,p->n);
        exit(EXIT_FAILURE);
    }

    if (p->nz>0) {
        q=Sparsealloc(NULL,p->t);
        q->m=p->m;
        q->n=p->n;
        nz_guess=q->m;
        q->r = ialloc(nz_guess,0);
        q->c = ialloc(nz_guess,0);
        if (q->t==DOUBLE) q->s.d = dalloc(nz_guess,0);
            else q->s.i = ialloc(nz_guess,0);

        nz=0;
        for (i=0;i<p->n;i++) {
            for (j=p->cc[i];j<p->cc[i+1];j++) {
                m=p->r[j]-1;
                if (m!=i) {
                    flag=0;
                    for (k=p->cc[m];k<p->cc[m+1];k++) {
                        if (p->r[k]==i+1) {
                            flag=1;
                            break;
                        }
                        if (p->r[k]>(i+1)) break;
                    }
                    if (flag==0) {
                        if (q->t==DOUBLE) q->s.d[nz] = 0.0;
                            else q->s.i[nz] = 0;
                        q->r[nz] = i+1;
                        q->c[nz] = m+1;
                        /* resize value and column array */
                        if (++nz==nz_guess) nz_guess=Sparsetripletrealloc(q,2*nz_guess);
                    }
                }
            }
        }
        if (nz>0) {
            /* fprintf(stderr,"ipSparsesymgraph: adding %d zeroes\n",nz); */
            err=Sparsetripletrealloc(p,p->nz+nz);
            memcpy(&p->r[p->nz],q->r,nz*sizeof(int));
            memcpy(&p->c[p->nz],q->c,nz*sizeof(int));
            if (p->t==DOUBLE) memcpy(&p->s.d[p->nz],q->s.d,nz*sizeof(double));
                else memcpy(&p->s.i[p->nz],q->s.i,nz*sizeof(int));
            p->nz+=nz;
            Sparsefree(q);
            err=countsortsparse(p->c,p->r,p->s.v,p->t,p->nz);

            /* construct the CCS array*/
            if (p->cc!=NULL) free(p->cc);
            p->cc=columntoCCS(p->c, p->n, p->nz);
            if (p->cc[p->n]!=p->nz) {
                fprintf(stderr,"ipSparsesymgraph: error: inconsistent nz! p->nz=%d, p->cc[n]=%d\n",p->nz,p->cc[p->n]);
                exit(EXIT_FAILURE);
            }
        }
    }

    return nz;
}

Sparse *Sparseindexselect(const Sparse *p, const Mat *a, const char rc)
{
    int i, j, k, err, *rhash, nz_guess, tally=0;
    Sparse *q=NULL;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->cc);
    assert(p->r);
    assert(p->c);
    assert(p->s.v);

    assert(a);
    assert(a->x==1);
    assert(a->y>0);
    assert(a->t==INTEGER);
    assert(a->m.v);

    if ((p->nz>0)&&(a->v.i[0]>0)) {
        q=Sparsealloc(NULL,p->t);
        nz_guess=p->nz/2;
        q->r=ialloc(nz_guess,0);
        q->c=ialloc(nz_guess,0);
        if (p->t==DOUBLE) q->s.d=dalloc(nz_guess,0);
            else q->s.i=ialloc(nz_guess,0);
        q->m=p->m;
        q->n=p->n;
        if ((rc=='r')&&(a->v.i[a->y-1]<=q->m)) {
            rhash=ialloc(q->m,1);
            for (i=0;i<a->y;i++) rhash[a->v.i[i]-1]=i+1;
    
            for (i=0;i<p->n;i++) {
                for (j=p->cc[i];j<p->cc[i+1];j++) {
                    k=rhash[p->r[j]-1];
                    if (k>0) {
                        q->c[tally]=i+1;
                        q->r[tally]=k;
                        if (q->t==DOUBLE) q->s.d[tally]=p->s.d[j];
                            else q->s.i[tally]=p->s.i[j];
                        if (++tally==nz_guess)
                            nz_guess=Sparsetripletrealloc(q,nz_guess*2);
                    }
                }
            }
            free(rhash);
            q->m=a->y;
        } else if ((rc=='c')&&(a->v.i[a->y-1]<=q->n)) {
            for (i=0;i<a->y;i++) {
                k=a->v.i[i];
                for (j=p->cc[k-1];j<p->cc[k];j++) {
                    q->c[tally]=i+1;
                    q->r[tally]=p->r[j];
                    if (q->t==DOUBLE) q->s.d[tally]=p->s.d[j];
                        else q->s.i[tally]=p->s.i[j];
                    if (++tally==nz_guess)
                        nz_guess=Sparsetripletrealloc(q,nz_guess*2);
                }
            }
            q->n=a->y;
        } else {
            fprintf(stderr,"Sparseindexselect: error: wrong vector char %c\n",rc);
            exit(EXIT_FAILURE);
        }

        q->nz=tally;
        if (q->nz>0) {       
            err=Sparsetripletrealloc(q,q->nz);
            q->cc=columntoCCS(q->c, q->n, q->nz);
        } else Sparsenullifypointers(q);
    }
    
    return q;
}

void ipSparseindexselect(Sparse *p, const Mat *a, const char rc)
{
    int i, j, k, err, *rhash, tally=0;

    assert(p);
    assert(p->m>0);
    assert(p->n>0);
    assert(p->cc);
    assert(p->r);
    assert(p->c);
    assert(p->s.v);

    assert(a);
    assert(a->x==1);
    assert(a->y>0);
    assert(a->t==INTEGER);
    assert(a->m.v);

    if ((p->nz>0)&&(a->v.i[0]>0)) {
        if ((rc=='r')&&(a->v.i[a->y-1]<=p->m)) {
            rhash=ialloc(p->m,1);
            for (i=0;i<a->y;i++) rhash[a->v.i[i]-1]=i+1;
    
            for (i=0;i<p->n;i++) {
                for (j=p->cc[i];j<p->cc[i+1];j++) {
                    k=rhash[p->r[j]-1];
                    if (k>0) {
                        p->c[tally]=i+1;
                        p->r[tally]=k;
                        if (p->t==DOUBLE) p->s.d[tally]=p->s.d[j];
                            else p->s.i[tally]=p->s.i[j];
                        tally++;
                    }
                }
            }
            free(rhash);
            p->m=a->y;
        } else if ((rc=='c')&&(a->v.i[a->y-1]<=p->n)) {
            for (i=0;i<a->y;i++) {
                k=a->v.i[i];
                for (j=p->cc[k-1];j<p->cc[k];j++) {
                    p->c[tally]=i+1;
                    p->r[tally]=p->r[j];
                    if (p->t==DOUBLE) p->s.d[tally]=p->s.d[j];
                        else p->s.i[tally]=p->s.i[j];
                    tally++;
                }
            }
            p->n=a->y;
        } else {
            fprintf(stderr,"ipSparseindexselect: error: wrong vector char %c\n",rc);
            exit(EXIT_FAILURE);
        }

        p->nz=tally;
        if (p->nz>0) {       
            err=Sparsetripletrealloc(p,p->nz);
            if (p->cc!=NULL) free(p->cc);
            p->cc=columntoCCS(p->c, p->n, p->nz);
        } else Sparsenullifypointers(p);
    }
    
    return;
}

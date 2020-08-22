/*
 * matrixutils.c
 * Version: 1.2
 *
 * Copyright (C) 2014 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#include <math.h>
#include "matrixutils.h"
#include "sortutils.h"

Mat3 *Mat3alloc(const int x, const int y, const int z, const int type, const short needzero)
{
    Mat3 *a;

    assert((type>=0) && (type<=1));
    a = malloc(sizeof(Mat3));
    assert(a);

    a->t = type;    
    if ((x<=0) || (y<=0) || (z<=0)) {
        a->x = 0;
        a->y = 0;
        a->z = 0;
        a->m.v = NULL;
        a->v.v = NULL;
    } else {
        a->x = x;
        a->y = y;
        a->z = z;
        if (type==DOUBLE) {
            a->m.d = darray3Dalloc(x,y,z,needzero);
            a->v.d=a->m.d[0][0];
        } else {
            a->m.i = iarray3Dalloc(x,y,z,needzero);
            a->v.i=a->m.i[0][0];
        }
    }
    
    return a;
}

void Mat3free(Mat3 *a)
{
    if (a != NULL) {
        if (a->m.v != NULL) {
            if (a->t==DOUBLE) darray3Dfree(a->m.d,a->x);
                else iarray3Dfree(a->m.i,a->x);
        }
        free(a);
    }
    return;
}

Mat3 **PMat3alloc(const int n)
{
    int i;
    Mat3 **a;

    a = malloc(n*sizeof(Mat3 *));
    assert(a);
    for (i=0;i<n;i++) a[i] = NULL;

    return a;
}

void PMat3free(Mat3 **a, const int n)
{
    int i;

    if (a != NULL) {
        for (i=0;i<n;i++) {
            if (a[i] != NULL) Mat3free(a[i]);
        }
        free(a);  
    }  
    return;
}

Mat3 *Mat3copy(const Mat3 *a)
{
    Mat3 *b;

    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->z>0);
    assert(a->m.v);
    assert((a->t>=0) && (a->t<=1));

    b=Mat3alloc(a->x,a->y,a->z,a->t,0);
   
    if (a->t==DOUBLE) memcpy(b->v.d,a->v.d,a->x*a->y*a->z*sizeof(double));
        else memcpy(b->v.i,a->v.i,a->x*a->y*a->z*sizeof(int));
    
    return b;
}

int ipMat3copy(const Mat3 *a, Mat3 *b)
{
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->z>0);
    assert(a->m.v);
    assert((a->t>=0) && (a->t<=1));
    assert(b);
    
    b->t = a->t;
    b->y = a->y;
    b->z = a->z;
    if (a->t==DOUBLE) {
        if (b->m.v != NULL) darray3Dfree(b->m.d,b->x);
        b->m.d=darray3Dalloc(a->x,a->y,a->z,0);
        b->v.d=b->m.d[0][0];
        memcpy(b->v.d,a->v.d,a->x*a->y*a->z*sizeof(double));
    } else {
        if (b->m.v != NULL) iarray3Dfree(b->m.i,b->x);
        b->m.i=iarray3Dalloc(a->x,a->y,a->z,0);
        b->v.i=b->m.i[0][0];
        memcpy(b->v.i,a->v.i,a->x*a->y*a->z*sizeof(int));
    }
    b->x = a->x;
    
    return EXIT_SUCCESS;
}

int ipMat3zero(Mat3 *a)
{
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->z>0);
    assert(a->m.v);
    assert((a->t>=0) && (a->t<=1));

    if (a->t==DOUBLE) {
        if (a->m.v != NULL) darray3Dfree(a->m.d,a->x);
        a->m.d=darray3Dalloc(a->x,a->y,a->z,1);
        a->v.d=a->m.d[0][0];
    } else {
        if (a->m.v != NULL) iarray3Dfree(a->m.i,a->x);
        a->m.i=iarray3Dalloc(a->x,a->y,a->z,1);
        a->v.i=a->m.i[0][0];
    }

    return EXIT_SUCCESS;    
}

Mat *Matalloc(const int x, const int y, const int type, const short needzero)
{
    Mat *a;

    assert((type>=0) && (type<=1));
    a = malloc(sizeof(Mat));
    assert(a);

    a->t = type;
    if ((x<=0) || (y<=0)) {
        a->x = 0;
        a->y = 0;
        a->m.v = NULL;
        a->v.v = NULL;
    } else {
        a->x = x;
        a->y = y;
        if (type==DOUBLE) {
            a->m.d = darrayalloc(x,y,needzero);
            a->v.d = a->m.d[0];
        } else {
            a->m.i = iarrayalloc(x,y,needzero);
            a->v.i = a->m.i[0];
        }
    }
    
    return a;
}

void Matfree(Mat *a)
{
    if (a != NULL) {
        if (a->m.v != NULL) {
            if (a->m.v[0] != NULL) free(a->m.v[0]);
            free(a->m.v);
        }
        free(a);
    }
    return;
}

Mat **PMatalloc(const int n)
{
    int i;
    Mat **a;

    a = malloc(n*sizeof(Mat *));
    assert(a);
    for (i=0;i<n;i++) a[i] = NULL;

    return a;
}

void PMatfree(Mat **a, const int n)
{
    int i;

    if (a != NULL) {
        for (i=0;i<n;i++) {
            if (a[i] != NULL) Matfree(a[i]);
        }
        free(a);  
    }  
    return;
}

Mat *Matcopy(const Mat *a)
{
    Mat *b;

    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert((a->t>=0) && (a->t<=1));
    
    b=Matalloc(a->x,a->y,a->t,0);

    if (a->t==DOUBLE) memcpy(b->v.d,a->v.d,a->x*a->y*sizeof(double));
        else memcpy(b->v.i,a->v.i,a->x*a->y*sizeof(int));
    
    return b;        
}

int ipMatcopy(const Mat *a, Mat *b)
{
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert((a->t>=0) && (a->t<=1));
    assert(b);
    
    b->t = a->t;
    b->x = a->x;
    b->y = a->y;
    if (a->t==DOUBLE) {
        if (b->m.v != NULL) darrayfree(b->m.d);
        b->m.d=darrayalloc(a->x,a->y,0);
        b->v.d=b->m.d[0];
        memcpy(b->v.d,a->v.d,a->x*a->y*sizeof(double));
    } else {
        if (b->m.v != NULL) iarrayfree(b->m.i);
        b->m.i=iarrayalloc(a->x,a->y,0);
        b->v.i=b->m.i[0];
        memcpy(b->v.i,a->v.i,a->x*a->y*sizeof(int));
    }
    
    return EXIT_SUCCESS;
}

int ipMatzero(Mat *a)
{
    int i;

    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert((a->t>=0) && (a->t<=1));

    if (a->t==DOUBLE) {
        free(a->m.d[0]);
        a->m.d[0]=dalloc(a->x*a->y,1);
        for (i=1; i<a->x; i++) a->m.d[i]=a->m.d[0] + (i*a->y);
        a->v.d=a->m.d[0];
    } else {
        free(a->m.i[0]);
        a->m.i[0]=ialloc(a->x*a->y,1);
        for (i=1; i<a->x; i++) a->m.i[i]=a->m.i[0] + (i*a->y);
        a->v.i=a->m.i[0];
    }

    return EXIT_SUCCESS;    
}

double Matdot(const Mat *a, const Mat *b)
{
    double d=0.0;
    M_INT n,incx=1,incy=1;
    
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==DOUBLE);
    assert(b);
    assert(b->x>0);
    assert(b->y>0);
    assert(b->m.v);
    assert(b->t==a->t);

    n=(M_INT)a->x*a->y;

    if (n != (M_INT)(b->x*b->y)) {
        fprintf(stderr,"Error in Matdot: a->x:%d a->y:%d\n",a->x,a->y);
        fprintf(stderr,"                 b->x:%d b->y:%d\n",b->x,b->y);
        exit (EXIT_FAILURE);
    } else {
        #ifndef ATLAS
            d=F77_FCN(ddot)(&n, a->v.d, &incx, b->v.d, &incy);
        #else
            d=cblas_ddot(n, a->v.d, incx, b->v.d, incy);
        #endif /*ATLAS*/
    }

    return d;
}

Mat *Matrot(const Mat *a, const Mat *b)
{
    Mat *c;
    
    assert(a);
    assert(a->x*a->y==3);
    assert(a->m.v);
    assert(a->t==DOUBLE);
    assert(b);
    assert(b->x*b->y==3);
    assert(b->m.v);
    assert(b->t==a->t);

    c=Matalloc(3,1,DOUBLE,0);

    if ((a->x*a->y!=3)||((b->x*b->y)!=3)) {
        fprintf(stderr,"Error in Matrot: a->x:%d a->y:%d\n",a->x,a->y);
        fprintf(stderr,"                 b->x:%d b->y:%d\n",b->x,b->y);
        exit (EXIT_FAILURE);
    } else {
        c->v.d[0]=a->v.d[1]*b->v.d[2]-a->v.d[2]*b->v.d[1];
        c->v.d[1]=a->v.d[2]*b->v.d[0]-a->v.d[0]*b->v.d[2];
        c->v.d[2]=a->v.d[0]*b->v.d[1]-a->v.d[1]*b->v.d[0];
    }

    return c;
}

Mat *Matadd(const Mat *a, const Mat *b, const double alpha)
{
    int i;
    Mat *c;
    
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert((a->t>=0) && (a->t<=1));
    assert(b);
    assert(b->x==a->x);
    assert(b->y==a->y);
    assert(b->m.v);
    assert(b->t==a->t);

    c=Matalloc(a->x,a->y,a->t,0);

    if (a->t==DOUBLE) for (i=0; i<a->x*a->y; i++) c->v.d[i]=a->v.d[i]+alpha*b->v.d[i];
        else for (i=0; i<a->x*a->y; i++) c->v.i[i]=a->v.i[i]+alpha*b->v.i[i];
   
    return c;
}

int ipMatadd(Mat *a, const Mat *b, const double alpha)
{
    int i,err;
    M_INT isize,incx=1,incy=1;
    
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert((a->t>=0) && (a->t<=1));
    assert(b);
    assert(b->x==a->x);
    assert(b->y==a->y);
    assert(b->m.v);
    assert(b->t==a->t);

    isize=(M_INT)a->x*a->y;
    if (a->t==DOUBLE)
        #ifndef ATLAS
            err=F77_FCN(daxpy)(&isize, &alpha, b->v.d, &incx, a->v.d, &incy);
        #else
            cblas_daxpy(isize, alpha, b->v.d, incx, a->v.d, incy);
        #endif /* ATLAS */
        else for (i=0; i<isize; i++) a->v.i[i]+=alpha*b->v.i[i];
   
    return EXIT_SUCCESS;
}

int Matrealloc(Mat *a, const int x, const int y)
{
    int i;
    assert(a);
    assert(a->x>=0);
    assert(a->y>=0);

    if ((x==a->x)&&(y==a->y)) return EXIT_SUCCESS;
    
    if ((x<=0) || (y<=0)) {
        a->x=0;
        a->y=0;
        if (a->m.v != NULL) {
            if (a->m.v[0] != NULL) free(a->m.v[0]);
            free(a->m.v);
        }
        a->m.v=NULL;
        a->v.v=NULL;

        return EXIT_SUCCESS;
    }

    if ((a->x==0)||(a->y==0)) {
        a->x=x;
        a->y=y;
        if (a->m.v != NULL) {
            if (a->m.v[0] != NULL) free(a->m.v[0]);
            free(a->m.v);
        }
        if (a->t==DOUBLE) a->m.d = darrayalloc(x,y,1);
            else a->m.i = iarrayalloc(x,y,1);

        a->v.v = a->m.v[0];

        return EXIT_SUCCESS;
    }

    if (x!=a->x) {
        if (a->t==DOUBLE) {
            a->m.d=realloc(a->m.d,x*sizeof(double *));
            assert(a->m.d);
            a->m.d[0]=realloc(a->m.d[0],x*a->y*sizeof(double));
            assert(a->m.d[0]);
            for (i=1; i<x; i++) a->m.d[i]=a->m.d[0] + (i*a->y);
            if (x>a->x) memset(a->m.d[a->x],0,(x-a->x)*a->y*sizeof(double));
        } else {
            a->m.i=realloc(a->m.i,x*sizeof(int *));
            assert(a->m.i);
            a->m.i[0]=realloc(a->m.i[0],x*a->y*sizeof(int));
            assert(a->m.i[0]);
            for (i=1; i<x; i++) a->m.i[i]=a->m.i[0] + (i*a->y);
            if (x>a->x) memset(a->m.i[a->x],0,(x-a->x)*a->y*sizeof(int));
        }
        a->v.v=a->m.v[0];
        a->x=x;
    }
    if (y!=a->y) {
        if (a->t==DOUBLE) {
            if (y>a->y) {
                a->m.d[0]=realloc(a->m.d[0],a->x*y*sizeof(double));
                assert(a->m.d[0]);
                /* if colums are expanded, start moving row last-to-first! */
                for (i=a->x-1; i>0; i--) {
                    memmove(a->m.d[0]+(i*y),a->m.d[0]+(i*a->y),a->y*sizeof(double));
                }
                for (i=0; i<a->x; i++) {
                    memset(a->m.d[0]+(i*y)+(a->y),0,(y-a->y)*sizeof(double));
                }
            } else {
                for (i=1; i<a->x; i++) {
                    memmove(a->m.d[0]+(i*y),a->m.d[0]+(i*a->y),y*sizeof(double));
                }
                a->m.d[0]=realloc(a->m.d[0],a->x*y*sizeof(double));
                assert(a->m.d[0]);
            }
            for (i=1; i<a->x; i++) a->m.d[i]=a->m.d[0] + (i*y);
        } else {
            if (y>a->y) {
                a->m.i[0]=realloc(a->m.i[0],a->x*y*sizeof(int));
                assert(a->m.i[0]);
                /* if colums are expanded, start moving row last-to-first! */
                for (i=a->x-1; i>0; i--) {
                    memmove(a->m.i[0]+(i*y),a->m.i[0]+(i*a->y),a->y*sizeof(int));
                }
                for (i=0; i<a->x; i++) {
                    memset(a->m.i[0]+(i*y)+(a->y),0,(y-a->y)*sizeof(int));
                }
            } else {
                for (i=1; i<a->x; i++) {
                    memmove(a->m.i[0]+(i*y),a->m.i[0]+(i*a->y),y*sizeof(int));
                }
                a->m.i[0]=realloc(a->m.i[0],a->x*y*sizeof(int));
                assert(a->m.i[0]);
            }
            for (i=1; i<a->x; i++) a->m.i[i]=a->m.i[0] + (i*y);
        }
        a->v.v=a->m.v[0];
        a->y=y;
    }

    return EXIT_SUCCESS;
}

double Matdet(const Mat *a)
{
    int i, err;
    M_INT *ipvt, minmn, info, tx, ty;
    double det=1.0;
    int neg=0;
    Mat *tmp;

    assert(a);
    assert(a->x>0);
    assert(a->y==a->x);
    assert(a->v.d);
    assert(a->t==DOUBLE);

    tmp=Mattranspose(a);
    minmn=(M_INT)(tmp->x<=tmp->y ? tmp->x : tmp->y);
    ipvt=(M_INT *)calloc(minmn, sizeof(M_INT));

    tx=(M_INT)tmp->x;
    ty=(M_INT)tmp->y;
  
    /* Note width and height are reversed */
    #ifndef ATLAS
        err=F77_FCN(dgetrf)(&ty, &tx, tmp->v.d, &minmn, ipvt, &info);
    #else
        info=clapack_dgetrf(CblasColMajor, ty, tx, tmp->v.d, minmn, ipvt);        
    #endif /* ATLAS */

    if(info > 0) {
        /* singular matrix */
        fprintf(stderr, "Matdet: dgetrf reports singular matrix\n");
        return 0.0;
    }

    /* Take the product of the diagonal elements */
    for (i=0; i < minmn; i++) {
        det *= tmp->m.d[i][i];
        if (ipvt[i] != (i+1)) neg = !neg;
    }
    free(ipvt);
    Matfree(tmp);

    /* Since tmp is an LU decomposition of a rowwise permutation of A,
        multiply by appropriate sign */
    return neg?-det:det;
}

Mat *MatLU(const Mat *a, Mat * l, Mat *u, const char c)
{
    int i, j, err;
    M_INT *ipvt, minmn, info, tx, ty;
    Mat *tmp,*b;
    double dummy;

    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->v.d);
    assert(a->t==DOUBLE);

    tmp=Mattranspose(a);
    minmn=(M_INT)(tmp->x<=tmp->y ? tmp->x : tmp->y);
    ipvt=(M_INT *)calloc(minmn, sizeof(M_INT));

    tx=(M_INT)tmp->x;
    ty=(M_INT)tmp->y;

    #ifndef ATLAS
        err=F77_FCN(dgetrf)(&ty, &tx, tmp->v.d, &minmn, ipvt, &info);
    #else
        info=clapack_dgetrf(CblasColMajor, ty, tx, tmp->v.d, minmn, ipvt);
    #endif /* ATLAS */

    if(info > 0) {
        /* singular matrix */
        fprintf(stderr, "MatLU: dgetrf reports singular matrix\n");
        b = NULL;
        return b;
    }

    for (i=0;i<minmn-1;i++) {
        for (j=i+1;j<minmn;j++) {
            l->m.d[j][i]=tmp->m.d[i][j];
            u->m.d[i][j]=tmp->m.d[j][i];
        }
    }    

    b=Matalloc(1,minmn,INTEGER,1);
    for (i=0;i<minmn;i++) {
        l->m.d[i][i]=1.0;
        u->m.d[i][i]=tmp->m.d[i][i];
        b->v.i[i]=i+1;
    }
    Matfree(tmp);

    for (i=0;i<minmn;i++) {
        if (i!=ipvt[i]-1) {
            dummy=b->v.i[i];
            b->v.i[i]=b->v.i[ipvt[i]-1];
            b->v.i[ipvt[i]-1]=dummy;
        }
    }
    free(ipvt);

    if (c=='v') {
        return b;
    } else if (c=='m') {
        tmp=Matalloc(minmn,minmn,DOUBLE,1);
        for (i=0;i<minmn;i++) tmp->m.d[i][b->v.i[i]-1]=1.0;
        Matfree(b);
        return tmp;
    } else {
        fprintf(stderr, "MatLU: unknown permutation/pivot switch %c\n",c);
        exit(EXIT_FAILURE);
    }
}

double ipMatinv(Mat *a)
{
    int i,err,neg=0;
    M_INT *ipvt, info, lwork, ax, ay, minmn;
    double det=1.0,*work;

    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==DOUBLE);

    ipvt = (M_INT *)calloc(a->x, sizeof(M_INT));
    ax = (M_INT)a->x;
    ay = (M_INT)a->y;
    minmn = (ax <= ay ? ax : ay);
    
    /* Note: width and height are reversed */
    err=ipMattranspose(a);
    
    #ifndef ATLAS
        err=F77_FCN(dgetrf)(&ay, &ax, a->v.d, &ax, ipvt, &info);
    #else
        info=clapack_dgetrf(CblasColMajor, ay, ax, a->v.d, ax, ipvt);
    #endif /* ATLAS */
    if(info > 0) {
        fprintf(stderr, "Matinv: dgetrf reports singular matrix\n");
        return 0.0;
    } else {
        /* Take the product of the diagonal elements */
        for (i=0; i < minmn; i++) {
            det *= a->m.d[i][i];
            if (ipvt[i] != (i+1)) neg = !neg;
        }
        #ifndef ATLAS      
            lwork = (M_INT)a->x * 3;
            work = dalloc(lwork, 1);
            err=F77_FCN(dgetri)(&ax, a->v.d, &ay, ipvt, work, &lwork, &info);
            free(work);
        #else
            info=clapack_dgetri(CblasColMajor, ax, a->v.d, ay, ipvt);
        #endif /* ATLAS */
        err=ipMattranspose(a);
        if(info > 0) {
            fprintf(stderr, "Matinv: dgetri reports singular matrix\n");
            return 0.0;
        }
    }

    free(ipvt);
    return neg?-det:det;
}

Mat *Matmul(const Mat *a, const Mat *b)
{
    return Matmulscal(a,b,1.0);
}

Mat *Matmulscal(const Mat *a, const Mat *b, const double alpha)
{
    int err;
    M_INT ax, ay, bx, by;
    double beta=0.0;
    Mat *c;
    
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==DOUBLE);
    assert(b);
    assert(b->x>0);
    assert(b->y>0);
    assert(b->m.v);
    assert(b->t==a->t);

    ax=(M_INT)a->x;
    ay=(M_INT)a->y;
    bx=(M_INT)b->x;
    by=(M_INT)b->y;

    if (a->y==b->x) {
        c=Matalloc(b->y,a->x,DOUBLE,1);
        /* Note width and height are reversed */
        #ifndef ATLAS
            err=F77_FCN(dgemm)("T", "T", &ax, &by, &ay, &alpha, a->v.d,
                &ay, b->v.d, &by, &beta, c->v.d, &ax);
        #else
            cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, ax, by, ay, alpha, a->v.d,
                ay, b->v.d, by, beta, c->v.d, ax);
        #endif /* ATLAS */
    } else if (a->x==b->y) {
        c=Matalloc(a->y,b->x,DOUBLE,1);
        /* Note width and height are reversed */
        #ifndef ATLAS
            err=F77_FCN(dgemm)("T", "T", &bx, &ay, &by, &alpha, b->v.d,
                &by, a->v.d, &ay, &beta, c->v.d, &bx);
        #else
            cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, bx, ay, by, alpha, b->v.d,
                by, a->v.d, ay, beta, c->v.d, bx);
        #endif /* ATLAS */
    } else {
        fprintf(stderr,"Error in Matmul: a->y:%d b->x:%d\n",a->y,b->x);
        fprintf(stderr,"                 a->x:%d b->y:%d\n",a->x,b->y);
        exit (EXIT_FAILURE);
    }
    err=ipMattranspose(c);

    return c;
}

Mat *Matmulbyvec(const Mat *a, const Mat *b)
{
    return Matmulscalbyvec(a,b,1.0);
}

Mat *Matmulscalbyvec(const Mat *a, const Mat *b, const double alpha)
{
    int i,err;
    M_INT ax, ay, bx, by, incx=1, incy=1;
    double beta=0.0;
    Mat *c;
    
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==DOUBLE);
    assert(b);
    assert(b->x>0);
    assert(b->y>0);
    assert(b->m.v);
    assert(b->t==a->t);

    ax=(M_INT)a->x;
    ay=(M_INT)a->y;
    bx=(M_INT)b->x;
    by=(M_INT)b->y;

    if (a->y==b->x) {
        c=Matalloc(b->y,a->x,DOUBLE,1);
        for (i=0;i<a->x;i++) {
            /* Note width and height are reversed */
            #ifndef ATLAS
                err=F77_FCN(dgemv)("T", &by, &bx, &alpha, b->v.d,
                    &by, a->m.d[i], &incx, &beta, c->m.d[i], &incy);
            #else
                cblas_dgemv(CblasColMajor, CblasTrans, by, bx, alpha, b->v.d,
                    by, a->m.d[i], incx, beta, c->m.d[i], incy);
            #endif /* ATLAS */
        }
    } else if (a->x==b->y) {
        c=Matalloc(a->y,b->x,DOUBLE,1);
        for (i=0;i<b->x;i++) {
            /* Note width and height are reversed */
            #ifndef ATLAS
                err=F77_FCN(dgemv)("T", &ay, &ax, &alpha, a->v.d,
                    &ay, b->m.d[i], &incx, &beta, c->m.d[i], &incy);
            #else
                cblas_dgemv(CblasColMajor, CblasTrans, ay, ax, alpha, a->v.d,
                    ay, b->m.d[i], incx, beta, c->m.d[i], incy);
            #endif /* ATLAS */
        }
    } else {
        fprintf(stderr,"Error in Matmulscalbyvec: a->y:%d b->x:%d\n",a->y,b->x);
        fprintf(stderr,"                          a->x:%d b->y:%d\n",a->x,b->y);
        exit (EXIT_FAILURE);
    }
    err=ipMattranspose(c);

    return c;
}

/* Matsquare: C=AT*A */
Mat *Matsquare(const Mat *a)
{
    return Matsquarescal(a,1.0);
}

/* Matsquarescal: C=AT*A*alpha */
Mat *Matsquarescal(const Mat *a, const double alpha)
{
    int err,i,j,k=0,size;
    M_INT ax, ay;
    double *b, beta=0.0;
    Mat *c;
    
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==DOUBLE);

    ax=(M_INT)a->x;
    ay=(M_INT)a->y;

    size=a->x*a->y;
    b=dalloc(size,0);
    for (i=0;i<a->x;i++)
        for (j=i;j<size;j+=a->x)
            b[j]=a->v.d[k++];

    c=Matalloc(a->y,a->y,DOUBLE,0);
    /* Note width and height are reversed */
    #ifndef ATLAS
        err=F77_FCN(dgemm)("T", "T", &ay, &ay, &ax, &alpha, b,
            &ax, a->v.d, &ay, &beta, c->v.d, &ay);
    #else
        cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, ay, ay, ax,
            alpha, b, ax, a->v.d, ay, beta, c->v.d, ay);
    #endif
    free(b);

    return c;
}

Mat *Matscal(const Mat *a, const double alpha)
{
    int err;
    M_INT isize,incx=1,incy=1;
    Mat *b;
    
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==DOUBLE);

    b=Matalloc(a->x,a->y,a->t,1);
    isize=(M_INT)a->x*a->y;
    #ifndef ATLAS
        err=F77_FCN(daxpy)(&isize, &alpha, a->v.d, &incx, b->v.d, &incy);
    #else
        cblas_daxpy(isize, alpha, a->v.d, incx, b->v.d, incy);
    #endif /* ATLAS */
    return b;
}

int ipMatscal(Mat *a, const double alpha)
{
    int err;
    M_INT isize,inc=1;
        
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==DOUBLE);
    
    isize=(M_INT)a->x*a->y;
    #ifndef ATLAS
        err=F77_FCN(dscal)(&isize, &alpha, a->v.d, &inc);
    #else
        cblas_dscal(isize, alpha, a->v.d, inc);
    #endif /* ATLAS */
        
    return EXIT_SUCCESS;
}

/* Mattransform: C=BT*A*B */
Mat *Mattransform(const Mat *a, const Mat *b)
{
    return Mattransformscal(a, b, 1.0);
}

/* Mattransformscal: C=BT*A*B*alpha */
Mat *Mattransformscal(const Mat *a, const Mat *b, const double alpha)
{
    int err;
    M_INT ax, ay, bx, by;
    double dummy=1.0, beta=0.0, *tmp;
    Mat *c, *bb;
    
    assert(a);
    assert(a->x>0);
    assert(a->y==a->x);
    assert(a->m.v);
    assert(a->t==DOUBLE);
    assert(b);
    assert(b->m.v);
    assert((b->x==a->y)||(b->y==a->x));
    assert(b->t==a->t);
        
    ax=(M_INT)a->x;
    ay=(M_INT)a->y;
    bx=(M_INT)b->x;
    by=(M_INT)b->y;

    if (a->y==b->x) {
        /* 1) tmp=A*B */
        /* 2) CT=tmpT*B because (C)T=(BT*tmp)T */
        /* 3) C=(CT)T */
        tmp=dalloc(b->y*a->x,1);
        /* Note width and height are reversed */
        c=Matalloc(b->y,b->y,DOUBLE,1);
        #ifndef ATLAS
            err=F77_FCN(dgemm)("T", "T", &ax, &by, &ay, &dummy, a->v.d,
                &ay, b->v.d, &by, &beta, tmp, &ax);
            err=F77_FCN(dgemm)("T", "T", &by, &by, &ax, &alpha, tmp,
                &ax, b->v.d, &by, &beta, c->v.d, &by);
        #else
            cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, ax, by, ay,
                dummy, a->v.d, ay, b->v.d, by, beta, tmp, ax);
            cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, by, by, ax,
                alpha, tmp, ax, b->v.d, by, beta, c->v.d, by);
        #endif /* ATLAS */
        free(tmp);
    } else if (a->x==b->y) {
        fprintf(stderr,"Mattransformscal : transposing B to avoid losing precision\n");
        bb=Mattranspose(b);
        bx=(M_INT)bb->x;
        by=(M_INT)bb->y;
        tmp=dalloc(bb->y*a->x,1);
        /* Note width and height are reversed */
        c=Matalloc(bb->y,bb->y,DOUBLE,1);
        #ifndef ATLAS
            err=F77_FCN(dgemm)("T", "T", &ax, &by, &ay, &dummy, a->v.d,
                &ay, bb->v.d, &by, &beta, tmp, &ax);
            err=F77_FCN(dgemm)("T", "T", &by, &by, &ax, &alpha, tmp,
                &ax, bb->v.d, &by, &beta, c->v.d, &by);
        #else
            cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, ax, by, ay,
                dummy, a->v.d, ay, bb->v.d, by, beta, tmp, ax);
            cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, by, by, ax,
                alpha, tmp, ax, bb->v.d, by, beta, c->v.d, by);
        #endif /* ATLAS */
        free(tmp);
        Matfree(bb);
    } else {
        fprintf(stderr,"Error in mattransformscal\n");
        exit (EXIT_FAILURE);
    }

    err=ipMattranspose(c);
    
    return c;
}

Mat *Mattranspose(const Mat *a)
{
    int i,j;
    Mat *b;
    
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    b=Matalloc(a->y,a->x,a->t,1);
    if (a->t==DOUBLE) {
        for (i=0; i<a->x;i++)
            for (j=0; j<a->y; j++)
                b->m.d[j][i]=a->m.d[i][j];
    } else {
        for (i=0; i<a->x;i++)
            for (j=0; j<a->y; j++)
                b->m.i[j][i]=a->m.i[i][j];
    }
            
    return b;
}

/* mldivide(A,B) = A\B => x = A^-1 * B */
/* mrdivide(A,B) = A/B = (B'\A')' => x = A * B^-1 */
Mat *Matdiv(const Mat *a, const Mat *b, const char rl)
{
    Mat *c,*x;
    int dummy,err,maxmn,minmn;
    M_INT *ipiv, info, lwork, cx, cy, xx, xy;
    double *work;
    /*M_INT rank;
    double rcond=-1.0;
    double *s;*/

    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==DOUBLE);

    assert(b);
    assert(b->x>0);
    assert(b->y>0);
    assert(b->m.v);
    assert(b->t==DOUBLE);

    if (rl=='l') {
        c=Mattranspose(a);
        x=Mattranspose(b);
    } else if (rl=='r') {
        c=Matcopy(b);
        x=Matcopy(a);
    } else {
        fprintf(stderr, "Matdiv: unknown division %c\n",rl);
        exit(EXIT_FAILURE);
    }

    cx=(M_INT)c->x;
    cy=(M_INT)c->y;
    xx=(M_INT)x->x;
    xy=(M_INT)x->y;

    if (c->y==x->y) {

        if (c->x==c->y) {

            /* A (and c) is square -> use dgesv */
            ipiv=(M_INT *)calloc(c->x,sizeof(M_INT));
            #ifndef ATLAS
                err=F77_FCN(dgesv)(&cx, &xx, c->v.d, &cy, ipiv, x->v.d, &xy,
                    &info);
            #else
                info=clapack_dgesv(CblasColMajor, cx, xx, c->v.d, cy, ipiv,
                    x->v.d, xy);
            #endif /* ATLAS */
            free(ipiv);
            if(info < 0) {
                fprintf(stderr, "Matdiv: dgesv reports illegal input %d\n",(int)-info);
                exit(EXIT_FAILURE);
            }
            if(info > 0) {
                fprintf(stderr, "Matdiv: dgesv reports singular factorization\n");
                exit(EXIT_FAILURE);
            }

        } else {

            /* A (and c) is m*n -> use dgels */
            minmn = MIN(c->x,c->y);
            maxmn = MAX(c->x,c->y);
            dummy = MAX(x->x,minmn);
            lwork = 2*(minmn+dummy);
            work = dalloc(lwork, 1);
            err=Matrealloc(x,x->x,maxmn);
            xx=(M_INT)x->x;
            xy=(M_INT)x->y;
            err=F77_FCN(dgels)("N", &cy, &cx, &xx, c->v.d, &cy, x->v.d,
                &xy, work, &lwork, &info);
            free(work);
            if(info < 0) {
                fprintf(stderr, "Matdiv: dgels reports illegal input %d\n",(int)-info);
                exit(EXIT_FAILURE);
            }
            if(info > 0) {
                fprintf(stderr, "Matdiv: dgels reports singular factorization in %d\n",(int)info);
                exit(EXIT_FAILURE);
            }
            err=Matrealloc(x,x->x,c->x);

            /* A (and c) is m*n -> use dgelss */
            /*minmn = MIN(c->x,c->y);
            maxmn = MAX(c->x,c->y);
            dummy = MAX(x->x,maxmn);
            s=dalloc(minmn,1);
            lwork = MAX(2*minmn,dummy);
            lwork += 3*minmn;
            lwork *= 2;
            work = dalloc(lwork, 1);
            err=Matrealloc(x,x->x,maxmn);
            xx=(M_INT)x->x;
            xy=(M_INT)x->y;
            err=F77_FCN(dgelss)(&cy, &cx, &xx, c->v.d, &cy, x->v.d,
                &xy, s, &rcond, &rank, work, &lwork, &info);
            free(s);
            free(work);
            if(info < 0) {
                fprintf(stderr, "Matdiv: dgelss reports illegal input %d\n",(int)-info);
                exit(EXIT_FAILURE);
            }
            if(info > 0) {
                fprintf(stderr, "Matdiv: dgelss reports singular factorization\n");
                exit(EXIT_FAILURE);
            }
            err=Matrealloc(x,x->x,c->x);*/

        }
    } else {
        if (rl=='l') {
            fprintf(stderr,"Error in Matdiv: %c a->x:%d b->x:%d\n",rl,c->y,x->y);
        } else {
            fprintf(stderr,"Error in Matdiv: %c a->y:%d b->y:%d\n",rl,c->y,x->y);
        }
        exit (EXIT_FAILURE);
    }
    Matfree(c);

    if (rl=='l') err=ipMattranspose(x);

    return x;
}

double Matmaxdouble(const Mat *a)
{
    int i;
    double max;
    
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==DOUBLE);
    
    max=a->v.d[0];
    for (i=1;i<a->x*a->y;i++) max=MAX(max,a->v.d[i]);

    return max;
}

double Matmindouble(const Mat *a)
{
    int i;
    double min;
    
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==DOUBLE);

    min=a->v.d[0];
    for (i=1;i<a->x*a->y;i++) min=MIN(min,a->v.d[i]);

    return min;
}

int Matmaxint(const Mat *a)
{
    int i;
    int max;
    
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==INTEGER);

    max=a->v.i[0];
    for (i=1;i<a->x*a->y;i++) max=MAX(max,a->v.i[i]);

    return max;
}

int Matminint(const Mat *a)
{
    int i;
    int min;
    
    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==INTEGER);

    min=a->v.i[0];
    for (i=1;i<a->x*a->y;i++) min=MIN(min,a->v.i[i]);

    return min;
}

int Matdump(const Mat *a)
{
    int i,j,mini,minj;
    
    assert(a);
    assert(a->x>=0);
    assert(a->y>=0);
    assert(a->m.v);
    assert((a->t>=0)&&(a->t<=1));

    if (a->t==DOUBLE) {
        printf("Matrixtype:Double x:%d y:%d\n",a->x, a->y);
        i = 5;
    } else {
        printf("Matrixtype:Integer x:%d y:%d\n",a->x, a->y);
        i = 25;
    }
    mini = (a->x < 50 ? a->x : 50);
    minj = (a->y < i ? a->y : i);
    for (i=0;i<mini;i++) {
        printf("i:%d",i);
        for (j=0;j<minj;j++) {
            if (a->t==DOUBLE)
                printf(" %+e",a->m.d[i][j]);
            else
                printf(" %+d",a->m.i[i][j]);
        }
        printf("\n");
    }

    return EXIT_SUCCESS;
}

int ipMatrowtrim(Mat *a, const Mat *b)
{
    int i,j;
    
    assert(a);
    assert(a->x>=0);
    assert(a->y>=0);
    assert(a->m.v);
    assert((a->t>=0)&&(a->t<=1));
    
    assert(b);
    assert(b->x==1);
    assert((b->y>0)&&(b->y<a->x));
    assert(b->m.v);
    assert(b->t==INTEGER);

    if (a->t==DOUBLE) {
        for (i=0;i<b->y;i++) {
            j=b->v.i[i]-i;
            if (j<=a->x) {
                memmove(a->m.d[j-1],a->m.d[j],(a->x-j)*a->y*sizeof(double));
                a->x--;
            }
        }
        a->m.d[0]=realloc(a->m.d[0],a->x*a->y*sizeof(double));
        assert(a->m.d[0]);
    	a->v.d=a->m.d[0];
        a->m.d=realloc(a->m.d,a->x*sizeof(double *));
        assert(a->m.d);
        for (i=0;i<a->x;i++) a->m.d[i]=a->m.d[0]+(i*a->y);
    } else {
        for (i=0;i<b->y;i++) {
            j=b->v.i[i]-i;
            if (j<=a->x) {
                memmove(a->m.i[j-1],a->m.i[j],(a->x-j)*a->y*sizeof(int));
                a->x--;
            }
        }
        a->m.i[0]=realloc(a->m.i[0],a->x*a->y*sizeof(int));
        assert(a->m.i[0]);
    	a->v.i=a->m.i[0];
        a->m.i=realloc(a->m.i,a->x*sizeof(int *));
        assert(a->m.i);
        for (i=0;i<a->x;i++) a->m.i[i]=a->m.i[0]+(i*a->y);
    }

    return EXIT_SUCCESS;
}

Mat *Matrowtrim(const Mat *a, const Mat *b)
{
    int i,j,k;
    Mat *c;
    
    assert(a);
    assert(a->x>=0);
    assert(a->y>=0);
    assert(a->m.v);
    assert((a->t>=0)&&(a->t<=1));
    
    assert(b);
    assert(b->x==1);
    assert((b->y>0)&&(b->y<a->x));
    assert(b->m.v);
    assert(b->t==INTEGER);

    c=Matalloc(a->x-b->y,a->y,a->t,1);

    j=0;
    k=0;
    if (a->t==DOUBLE) {
        for (i=0;i<a->x;i++) {
            if (i==b->v.i[j]-1) {
                if (j<b->y-1) j++;
                k++;
            } else {  
                memcpy(c->m.d[i-k],a->m.d[i],a->y*sizeof(double));
            }
        }
    } else {
        for (i=0;i<a->x;i++) {
            if (i==b->v.i[j]-1) {
                if (j<b->y-1) j++;
                k++;
            } else {  
                memcpy(c->m.i[i-k],a->m.i[i],a->y*sizeof(int));
            }
        }
    }

    return c;
}

int ipMatrowsplice(Mat *a, const Mat *b)
{
    int i,j;
    
    assert(a);
    assert(a->x>=0);
    assert(a->y>=0);
    assert(a->m.v);
    assert((a->t>=0)&&(a->t<=1));
    
    assert(b);
    assert(b->x==1);
    assert((b->y>0)&&(b->y<=a->x));
    assert(b->m.v);
    assert(b->t==INTEGER);

    a->x+=b->y;
    j=b->y;
    if (a->t==DOUBLE) {
        a->m.d[0]=realloc(a->m.d[0],a->x*a->y*sizeof(double));
        assert(a->m.d[0]);
    	a->v.d=a->m.d[0];
        a->m.d=realloc(a->m.d,a->x*sizeof(double *));
        assert(a->m.d);
        for (i=1;i<a->x;i++) a->m.d[i]=a->m.d[0]+(i*a->y);

        i=a->x;
        while ((i!=0)&&(j>0)) {
            i--;
            if (i==b->v.i[j-1]-1) {
                memset(a->m.d[i],0,a->y*sizeof(double));
                j--;
            } else memmove(a->m.d[i],a->m.d[i-j],a->y*sizeof(double));
        }
   } else {
        a->m.i[0]=realloc(a->m.i[0],a->x*a->y*sizeof(int));
        assert(a->m.i[0]);
    	a->v.i=a->m.i[0];
        a->m.i=realloc(a->m.i,a->x*sizeof(int *));
        assert(a->m.i);
        for (i=1;i<a->x;i++) a->m.i[i]=a->m.i[0]+(i*a->y);

        i=a->x;
        while ((i!=0)&&(j>0)) {
            i--;
            if (i==b->v.i[j-1]-1) {
                memset(a->m.i[i],0,a->y*sizeof(int));
                j--;
            } else memmove(a->m.i[i],a->m.i[i-j],a->y*sizeof(int));
        }
    }
    
    return EXIT_SUCCESS;
}

int ipMatcat(Mat *a, const Mat *b, const char rc)
{
    int err,i,j,k;
    
    assert(a);
    assert(a->x>=0);
    assert(a->y>=0);
    assert((a->t>=0)&&(a->t<=1));
    
    assert(b);
    assert(b->x>0);
    assert(b->y>0);
    assert(b->m.v);
    assert(b->t==a->t);

    if (rc=='c') {
        j=MAX(a->x,b->x);
        k=a->y;
        err=Matrealloc(a,j,a->y+b->y);
        if (a->t==DOUBLE) for (i=0;i<b->x;i++) memcpy(a->m.d[i]+(k),b->m.d[i],b->y*sizeof(double));
            else for (i=0;i<b->x;i++) memcpy(a->m.i[i]+(k),b->m.i[i],b->y*sizeof(int));

    } else {
        j=MAX(a->y,b->y);
        k=a->x;
        err=Matrealloc(a,a->x+b->x,j);
        if (a->t==DOUBLE) for (i=0;i<b->x;i++) memcpy(a->m.d[i+k],b->m.d[i],b->y*sizeof(double));
            else for (i=0;i<b->x;i++) memcpy(a->m.i[i+k],b->m.i[i],b->y*sizeof(int));
    }

    return EXIT_SUCCESS;
}

Mat *Matcut(const Mat *a, const int offset, const int count, const char rc)
{
    int i;
    Mat *b;
    
    assert(a);
    assert(a->x>=0);
    assert(a->y>=0);
    assert(a->m.v);
    assert((a->t>=0)&&(a->t<=1));
    assert(offset>=0);
    assert(count>0);

    if (rc=='c') {
        if ((offset+count)<=a->y) {
            b=Matalloc(a->x,count,a->t,0);
            if (a->t==DOUBLE) for (i=0;i<a->x;i++) memcpy(b->m.d[i],a->m.d[i]+(offset),count*sizeof(double));
                else for (i=0;i<a->x;i++) memcpy(b->m.i[i],a->m.i[i]+(offset),count*sizeof(int));
        } else {
            fprintf(stderr,"Error: in Matcut\n");
            fprintf(stderr,"offset %d and count %d too large for matrix: a->y=%d\n",offset,count,a->y);
            exit(EXIT_FAILURE);
        }
    } else {
        if ((offset+count)<=a->x) {
            b=Matalloc(count,a->y,a->t,0);
            if (a->t==DOUBLE) memcpy(b->v.d,a->m.d[offset],count*a->y*sizeof(double));
                else memcpy(b->v.i,a->m.i[offset],count*a->y*sizeof(int));
        } else {
            fprintf(stderr,"Error: in Matcut\n");
            fprintf(stderr,"offset %d and count %d too large for matrix: a->x=%d\n",offset,count,a->x);
            exit(EXIT_FAILURE);
        }
    }        

    return b;
}

int *ipMatrowsort(Mat *a)
{
    int i,j,k,l,m,p,flag,*ib=NULL,*c,*idx,*snidx,*didx,idummy;
    void *v;
    double dummy,*db=NULL;

    assert(a);
    assert(a->x>=0);
    assert(a->y>=0);
    assert((a->t==0)&&(a->t<=1));

    /* create local containers */
    if (a->t==DOUBLE) db=dalloc(a->y,0);
        else ib=ialloc(a->y,0);
    /* create index array for counting identical supernode columns */
    c=ialloc(a->y,1);
    /* create supernode index array for returning sort order */
    snidx=ialloc(a->y,0);
    for (i=0;i<a->y;i++) snidx[i]=i;
    /* first entry is always the entire length of row 0 */
    c[0]=a->y;
    
    /* loop over all rows */
    for (i=0;i<a->x;i++) {
        p=0;
        flag=0;
        /* as long as the column pointer hasn't reached the end, do loop */
        while (p<a->y) {
            l=c[p];
            /* are there more than 1 identical supernodes? */
            if (l>1) {
                flag=1;
                /* yes, sort the supernodes */
                if (a->t==DOUBLE) v=(void *)(a->m.d[i]+p);
                    else v=(void *)(a->m.i[i]+p);
                idx=heapsortindex(v,a->t,l);
                /* update supernode index */
                didx=ialloc(l,0);
                for (j=0;j<l;j++) didx[j]=snidx[p+idx[j]];
                for (j=0;j<l;j++) snidx[p+j]=didx[j];
                free(didx);
                if (i+1 < a->x) {
                    /* sort all entries below the supernodes */
                    for (j=i+1;j<a->x;j++) {
                        if (a->t==DOUBLE) {
                            for (k=0;k<l;k++) db[k]=a->m.d[j][idx[k]+p];
                            for (k=0;k<l;k++) a->m.d[j][k+p]=db[k];
                        } else {
                            for (k=0;k<l;k++) ib[k]=a->m.i[j][idx[k]+p];
                            for (k=0;k<l;k++) a->m.i[j][k+p]=ib[k];
                        }
                    }
                    k=0;
                    m=1;
                    c[p+k]=m;
                    /* detect whether we have identical entries on the row */
                    if (a->t==DOUBLE) {
                        dummy=a->m.d[i][p+k];
                        for (j=k+1;j<l;j++) {
                            /* identical? */
                            if (a->m.d[i][p+j]==dummy) {
                                /* good, bump the count */
                                m++;
                            } else {
                                /* update c[] */
                                dummy=a->m.d[i][p+j];
                                c[p+k]=m;
                                k+=m;
                                m=1;
                            }
                            /* check for last entry */
                            if (j+1==l) c[p+k]=m;
                        }
                    } else {
                        /* similar for integers */
                        idummy=a->m.i[i][p+k];
                        for (j=k+1;j<l;j++) {
                            if (a->m.i[i][p+j]==idummy) {
                                m++;
                            } else {
                                idummy=a->m.i[i][p+j];
                                c[p+k]=m;
                                k+=m;
                                m=1;
                            }
                            if (j+1==l) c[p+k]=m;
                        }
                    }
                }
                free(idx);
            }
            p+=l;
        }
        if (flag==0) break;
/*        for (j=0;j<a->y;j++) printf("%d ",c[j]);
        printf("\n");*/
    }

    if (a->t==DOUBLE) free(db);
        else free(ib);
    free(c);

    return snidx;
}

Mat *Matread(const char *filename)
{
    int err,i,m,n,t;
    Mat *a;
    FILE *f;

    f=fopen(filename,"r");
    err=fscanf(f,"%d %d %d",&m,&n,&t);
    assert(m>0);
    assert(n>0);
    assert((t>=0)&&(t<=1));

    a=Matalloc(m,n,t,0);
    assert(a);
    assert(a->v.v);

    if (a->t==DOUBLE)
        for (i=0;i<m*n;i++) err=fscanf(f,"%lf",&(a->v.d[i]));
    else
        for (i=0;i<m*n;i++) err=fscanf(f,"%d",&(a->v.i[i]));
    fclose(f);
    
    return a;
}

int Matwrite(const char *filename, const Mat *a, const int append)
{
    int i;
    FILE *f;

    assert(a);
    assert(a->x>=0);
    assert(a->y>=0);
    assert(a->v.v);
    assert((a->t>=0)&&(a->t<=1));

    if (append) f=fopen(filename,"a+");        
        else f=fopen(filename,"w");

    fprintf(f,"%d %d %d\n",a->x,a->y,a->t);
    if (a->t==DOUBLE)
        for (i=0;i<a->x*a->y;i++) fprintf(f,"%22.18e\n",a->v.d[i]);
    else
        for (i=0;i<a->x*a->y;i++) fprintf(f,"%d\n",a->v.i[i]);
    fclose(f); 
   
    return EXIT_SUCCESS;
}

/*
 * Matread_mtx(filename)
 * reads a Matrix Market file into a Mat structure
 *
 * modified from http://math.nist.gov/MatrixMarket/mmio/c/example_read.c
 */
Mat *Matread_mtx(const char *filename)
{
    Mat *p;
    int i,j,err,m,n,di;
    double dd;
    FILE *fp;
    MM_typecode matcode;

    fp=fopen(filename,"r");

    /* read the %% banner */
    if (mm_read_banner(fp, &matcode) != 0) {
        fprintf(stderr,"Could not process Matrix Market banner.\n");
        exit (EXIT_FAILURE);
    }

    /* exit if the matrix is complex or sparse */
    if (mm_is_complex(matcode) || mm_is_sparse(matcode)) {
        fprintf(stderr,"Matread_mtx() does not support ");
        fprintf(stderr,"Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(EXIT_FAILURE);
    }

    /* read size and dimensions of the dense matrix */
    if ((err = mm_read_mtx_array_size(fp, &m, &n)) !=0) {
        fprintf(stderr,"Could not determine dense matrix dimensions.");
        exit(EXIT_FAILURE);
    }

    if (mm_is_real(matcode)) p=Matalloc(m,n,DOUBLE,0);
        else if (mm_is_integer(matcode)) p=Matalloc(m,n,INTEGER,0);
        else {
        fprintf(stderr,"Matread_mtx() does not support ");
        fprintf(stderr,"Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(EXIT_FAILURE);
    }
    
    if (p->t==DOUBLE) {
        for (j=0;j<n;j++) {
            for (i=0;i<m;i++) {
                err=fscanf(fp,"%lg",&dd);
                p->m.d[i][j]=dd;
            }
        }
    } else {
        for (j=0;j<n;j++) {
            for (i=0;i<m;i++) {
                err=fscanf(fp,"%d",&di);
                p->m.i[i][j]=di;
            }
        }
    }
    fclose(fp);

    return p;
}

/*
 * Matwrite_mtx(filename,p)
 * Writes Dense Matrix p into a Matrix Market file
 *
 * modified from http://math.nist.gov/MatrixMarket/mmio/c/example_write.c
 */
void Matwrite_mtx(const char *filename, const Mat *p)
{
    int i,j;
    FILE *fp;
    MM_typecode matcode;

    assert(p);
    assert(p->x>0);
    assert(p->x>0);
    assert((p->t>=0)&&(p->t<=1));
    assert(p->m.v);
    assert(p->v.v);

    fp=fopen(filename,"w");

    mm_initialize_typecode(&matcode);
    mm_set_dense(&matcode);
    mm_set_general(&matcode);
    if (p->t==DOUBLE)
        mm_set_real(&matcode);
    else
        mm_set_integer(&matcode);

    mm_write_banner(fp, matcode); 
    mm_write_mtx_array_size(fp, p->x, p->y);

    if (p->t==DOUBLE) {
        for (j=0;j<p->y;j++) {
            for (i=0;i<p->x;i++) fprintf(fp,"%22.18e\n",p->m.d[i][j]);
        }
    } else {
        for (j=0;j<p->y;j++) {
            for (i=0;i<p->x;i++) fprintf(fp,"%d\n",p->m.i[i][j]);
        }
    }
    fclose(fp);

    return;
}

#include "matrixutils.h"

Mat *Matnull(const Mat *a, const double tolerance, const char dc)
{
    int i, j, minmn, maxmn, err, rank;
    M_INT *iwork, info, lwork, ax, ay, bx, by;
    double *s,*u,*work;

    Mat *b,*vt;

    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);
    assert(a->t==DOUBLE);

    b=Mattranspose(a);

    minmn = MIN(b->x,b->y);
    maxmn = MAX(b->x,b->y);
    s=dalloc(minmn,1);
    u=dalloc(b->y*b->y,1);
    vt=Matalloc(b->x,b->x,DOUBLE,1);
    if (dc=='d') lwork = 2*(MAX(maxmn,4*minmn*(minmn+1))+3*minmn);
        else lwork= 2*MAX(3*minmn+maxmn,5*minmn);
    work = dalloc(lwork, 1);
    assert(work);

    ax=(M_INT)a->x;
    ay=(M_INT)a->y;
    bx=(M_INT)b->x;
    by=(M_INT)b->y;

    if (dc=='d') {
        iwork = (M_INT *)calloc(8*minmn, sizeof(M_INT));
        err=F77_FCN(dgesdd)("A", &by, &bx, b->v.d, &by, s, u, &by,
            vt->v.d, &bx, work, &lwork, iwork, &info);
        free(iwork);
    } else {
        err=F77_FCN(dgesvd)("N", "A", &by, &bx, b->v.d, &by, s, u, &by,
            vt->v.d, &bx, work, &lwork, &info);
    }
    free(u);
    free(work);
    Matfree(b);
    
    rank=0;
    for (i=0;i<minmn;i++) {
        if (s[i]>maxmn*s[0]*tolerance) rank++;
    }
    free(s);

    b=Matalloc(vt->x,vt->y-rank,DOUBLE,0);
    for (j=0;j<b->x;j++)
        for (i=0;i<b->y;i++)
            b->m.d[j][i]=vt->m.d[j][i+rank];

    Matfree(vt);

    return b;
}

/*
 * fasttranspose.c
 * Version: 1.2
 *
 * Copyright (C) 2014 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#include "matrixutils.h"

typedef void (*txpose)(double *);
txpose tfparray[8*8] = {NULL};
short tfpinit = 1;

void T23(double *a)
{
    double tmp;
    tmp=a[1]; a[1]=a[3]; a[3]=a[4]; a[4]=a[2]; a[2]=tmp;
    return;
}

void T32(double *a)
{
    double tmp;
    tmp=a[1]; a[1]=a[2]; a[2]=a[4]; a[4]=a[3]; a[3]=tmp;
    return;
}

void T24(double *a)
{
    double tmp1,tmp2;
    tmp1=a[4]; a[4]=a[2]; a[2]=a[1]; a[1]=tmp1;
    tmp2=a[5]; a[5]=a[6]; a[6]=a[3]; a[3]=tmp2;
    return;
}

void T42(double *a)
{
    double tmp1,tmp2;
    tmp1=a[1]; a[1]=a[2]; a[2]=a[4]; a[4]=tmp1;
    tmp2=a[3]; a[3]=a[6]; a[6]=a[5]; a[5]=tmp2;
    return;
}

void T34(double *a)
{
    double tmp1,tmp2;
    tmp1=a[4]; a[4]= a[5];  a[5]=a[9]; a[9]=a[3]; a[3]=a[1]; a[1]=tmp1;
    tmp2=a[8]; a[8]=a[10]; a[10]=a[7]; a[7]=a[6]; a[6]=a[2]; a[2]=tmp2;
    return;
}

void T36(double *a)
{
    double tmp;
     tmp= a[1];  a[1]= a[6];  a[6]= a[2];  a[2]=a[12]; a[12]= a[4];  a[4]=a[7];
    a[7]= a[8];  a[8]=a[14]; a[14]=a[16]; a[16]=a[11]; a[11]=a[15]; a[15]=a[5];
    a[5]=a[13]; a[13]=a[10]; a[10]= a[9];  a[9]= a[3];  a[3]=tmp;
    return;
}

void T38(double *a)
{
    double tmp1,tmp2;
     tmp1= a[1]; a[1]= a[8];  a[8]=a[18]; a[18]= a[6];  a[6]= a[2];  a[2]=a[16];
    a[16]=a[13];a[13]=a[12]; a[12]= a[4];  a[4]= a[9];  a[9]= a[3];  a[3]=tmp1;
     tmp2=a[5];  a[5]=a[17]; a[17]=a[21]; a[21]= a[7];  a[7]=a[10]; a[10]=a[11];
    a[11]=a[19];a[19]=a[14]; a[14]=a[20]; a[20]=a[22]; a[22]=a[15]; a[15]=tmp2;

    return;
}

void T63(double *a)
{
    double tmp;
      tmp= a[3];  a[3]= a[9];  a[9]=a[10]; a[10]=a[13]; a[13]=a[5]; a[5]=a[15];
    a[15]=a[11]; a[11]=a[16]; a[16]=a[14]; a[14]= a[8];  a[8]=a[7]; a[7]= a[4];
     a[4]=a[12]; a[12]= a[2];  a[2]= a[6];  a[6]= a[1];  a[1]=tmp;
    return;
}

void T43(double *a)
{
    double tmp1,tmp2;
    tmp1=a[1]; a[1]=a[3]; a[3]=a[9]; a[9]= a[5];  a[5]=a[4]; a[4]=tmp1;
    tmp2=a[2]; a[2]=a[6]; a[6]=a[7]; a[7]=a[10]; a[10]=a[8]; a[8]=tmp2;
    return;
}

void T48(double *a)
{
    double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    tmp1= a[1];  a[1]= a[8];  a[8]= a[2];  a[2]=a[16]; a[16]= a[4];  a[4]=tmp1;
    tmp2= a[3];  a[3]=a[24]; a[24]= a[6];  a[6]=a[17]; a[17]=a[12]; a[12]=tmp2;
    tmp3= a[5];  a[5]= a[9];  a[9]=a[10]; a[10]=a[18]; a[18]=a[20]; a[20]=tmp3;
    tmp4= a[7];  a[7]=a[25]; a[25]=a[14]; a[14]=a[19]; a[19]=a[28]; a[28]=tmp4;
    tmp5=a[11]; a[11]=a[26]; a[26]=a[22]; a[22]=a[21]; a[21]=a[13]; a[13]=tmp5;
    tmp6=a[15]; a[15]=a[27]; a[27]=a[30]; a[30]=a[23]; a[23]=a[29]; a[29]=tmp6;
    return;
}

void T68(double *a)
{
    double tmp1,tmp2;
     tmp1= a[1];  a[1]= a[8];  a[8]=a[17]; a[17]=a[42]; a[42]= a[7];
     a[7]= a[9];  a[9]=a[25]; a[25]=a[12]; a[12]= a[2];  a[2]=a[16];
    a[16]=a[34]; a[34]=a[37]; a[37]=a[14]; a[14]=a[18]; a[18]= a[3];
     a[3]=a[24]; a[24]= a[4];  a[4]=a[32]; a[32]=a[21]; a[21]=a[27];
    a[27]=a[28]; a[28]=a[36]; a[36]= a[6];  a[6]=tmp1;
     tmp2= a[5];  a[5]=a[40]; a[40]=a[38]; a[38]=a[22]; a[22]=a[35];
    a[35]=a[45]; a[45]=a[31]; a[31]=a[13]; a[13]=a[10]; a[10]=a[33];
    a[33]=a[29]; a[29]=a[44]; a[44]=a[23]; a[23]=a[43]; a[43]=a[15];
    a[15]=a[26]; a[26]=a[20]; a[20]=a[19]; a[19]=a[11]; a[11]=a[41];
    a[41]=a[46]; a[46]=a[39]; a[39]=a[30]; a[30]=tmp2;
    return;
}

void T83(double *a)
{
    double tmp1,tmp2;
     tmp1= a[3]; a[3]= a[9];  a[9]= a[4];  a[4]=a[12]; a[12]=a[13]; a[13]=a[16];
    a[16]= a[2]; a[2]= a[6];  a[6]=a[18]; a[18]= a[8];  a[8]=a[1];   a[1]=tmp1;
     tmp2= a[5]; a[5]=a[15]; a[15]=a[22]; a[22]=a[20]; a[20]=a[14]; a[14]=a[19];
    a[19]=a[11];a[11]=a[10]; a[10]= a[7];  a[7]=a[21]; a[21]=a[17]; a[17]=tmp2;
    return;
}

void T84(double *a)
{
    double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    tmp1= a[4];  a[4]=a[16]; a[16]= a[2];  a[2]= a[8];  a[8]= a[1];  a[1]=tmp1;
    tmp2=a[12]; a[12]=a[17]; a[17]= a[6];  a[6]=a[24]; a[24]= a[3];  a[3]=tmp2;
    tmp3=a[20]; a[20]=a[18]; a[18]=a[10]; a[10]= a[9];  a[9]= a[5];  a[5]=tmp3;
    tmp4=a[28]; a[28]=a[19]; a[19]=a[14]; a[14]=a[25]; a[25]= a[7];  a[7]=tmp4;
    tmp5=a[13]; a[13]=a[21]; a[21]=a[22]; a[22]=a[26]; a[26]=a[11]; a[11]=tmp5;
    tmp6=a[29]; a[29]=a[23]; a[23]=a[30]; a[30]=a[27]; a[27]=a[15]; a[15]=tmp6;
    return;
}

void T86(double *a)
{
    double tmp1,tmp2;
     tmp1= a[6];  a[6]=a[36]; a[36]=a[28]; a[28]=a[27]; a[27]=a[21];
    a[21]=a[32]; a[32]= a[4];  a[4]=a[24]; a[24]= a[3];  a[3]=a[18];
    a[18]=a[14]; a[14]=a[37]; a[37]=a[34]; a[34]=a[16]; a[16]= a[2];
     a[2]=a[12]; a[12]=a[25]; a[25]= a[9];  a[9]= a[7];  a[7]=a[42];
    a[42]=a[17]; a[17]= a[8];  a[8]= a[1];  a[1]=tmp1;
     tmp2=a[30]; a[30]=a[39]; a[39]=a[46]; a[46]=a[41]; a[41]=a[11];
    a[11]=a[19]; a[19]=a[20]; a[20]=a[26]; a[26]=a[15]; a[15]=a[43];
    a[43]=a[23]; a[23]=a[44]; a[44]=a[29]; a[29]=a[33]; a[33]=a[10];
    a[10]=a[13]; a[13]=a[31]; a[31]=a[45]; a[45]=a[35]; a[35]=a[22];
    a[22]=a[38]; a[38]=a[40]; a[40]= a[5];  a[5]=tmp2;
    return;
}

int ipMattranspose(Mat *a)
{
    Mat *b;
    int i,j,ax,ay,idummy;
    double dummy;
    txpose tfp;

    if (tfpinit) {
        tfparray[8*(2-1)+3-1]=&T23;
        tfparray[8*(2-1)+4-1]=&T24;
        tfparray[8*(3-1)+2-1]=&T32;
        tfparray[8*(3-1)+4-1]=&T34;
        tfparray[8*(3-1)+6-1]=&T36;
        tfparray[8*(3-1)+8-1]=&T38;
        tfparray[8*(4-1)+2-1]=&T42;
        tfparray[8*(4-1)+3-1]=&T43;
        tfparray[8*(4-1)+8-1]=&T48;
        tfparray[8*(6-1)+3-1]=&T63;
        tfparray[8*(6-1)+8-1]=&T68;
        tfparray[8*(8-1)+3-1]=&T83;
        tfparray[8*(8-1)+4-1]=&T84;
        tfparray[8*(8-1)+6-1]=&T86;
        tfpinit=0;
    }

    assert(a);
    assert(a->x>0);
    assert(a->y>0);
    assert(a->m.v);

    ax=a->x;
    ay=a->y;

    /* scalar */
    if ((ax==1)&&(ay==1)) return EXIT_SUCCESS;

    /* row vector */
    if ((ax==1)&&(ay>1)) {
        if (a->t==DOUBLE) {
            a->m.d=realloc(a->m.d,ay*sizeof(double *));
            for (i=1;i<ay;i++) a->m.d[i]=a->m.d[0]+(i);
        } else {
            a->m.i=realloc(a->m.i,ay*sizeof(int *));
            for (i=1;i<ay;i++) a->m.i[i]=a->m.i[0]+(i);
        }
        a->x=ay;
        a->y=1;
        return EXIT_SUCCESS;
    /* column vector */
    } else if ((ay==1)&&(ax>1)) {
        if (a->t==DOUBLE) {
            a->m.d=realloc(a->m.d,sizeof(double *));
        } else {
            a->m.i=realloc(a->m.i,sizeof(int *));
        }
        a->y=ax;
        a->x=1;
        return EXIT_SUCCESS;
    /* square matrix */
    } else if (ax==ay) {
        if (a->t==DOUBLE) {
            for (i=0;i<ax-1;i++) {
                for (j=i+1;j<ay;j++) {
                    dummy=a->m.d[j][i];
                    a->m.d[j][i]=a->m.d[i][j];
                    a->m.d[i][j]=dummy;
                }
            }
        } else {
            for (i=0;i<ax-1;i++) {
                for (j=i+1;j<ay;j++) {
                    idummy=a->m.i[j][i];
                    a->m.i[j][i]=a->m.i[i][j];
                    a->m.i[i][j]=idummy;
                }
            }
        }
        return EXIT_SUCCESS;
    }

    if ((a->t==DOUBLE)&&(ax<=8)&&(ay<=8)) {
        if ((tfp=tfparray[8*(ax-1)+(ay-1)])!=NULL) {
            (*tfp)(a->v.d);
            a->m.d=realloc(a->m.d,ay*sizeof(double *));
            for (i=1;i<ay;i++) a->m.d[i]=a->m.d[0]+(i*ax);
            a->x=ay;a->y=ax;
            return EXIT_SUCCESS;
        }
    }

    b=Mattranspose(a);
#ifndef NDEBUG
    printf("Slow T: %d %d\n",b->x,b->y);
#endif
    a->x=b->x;
    a->y=b->y;
    free(a->m.v[0]);
    free(a->m.v);
    a->m.v=b->m.v;
    a->v.v=b->m.v[0];
    free(b);
    
    return EXIT_SUCCESS;
}

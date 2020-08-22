/*
 * matrixutils.h
 * Version: 1.2
 *
 * Copyright (C) 2014 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#ifndef MATRIXUTILS_H
#define MATRIXUTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "mat2d3d.h"
#include "matrixtype.h"
#include "arrayutils.h"
#include "finterface.h"
#include "mmio.h"

Mat3 *Mat3alloc(const int, const int, const int, const int, const short);
void Mat3free(Mat3 *);
Mat3 **PMat3alloc(const int);
void PMat3free(Mat3 **, const int);
Mat3 *Mat3copy(const Mat3 *);
int ipMat3copy(const Mat3 *, Mat3 *);
int ipMat3zero(Mat3 *);

Mat *Matalloc(const int, const int, const int, const short);
void Matfree(Mat *);
Mat **PMatalloc(const int);
void PMatfree(Mat **, const int);
Mat *Matcopy(const Mat *);
int ipMatcopy(const Mat *, Mat *);
int ipMatzero(Mat *);
Mat *Matadd(const Mat *, const Mat *, const double);
int ipMatadd(Mat *, const Mat *, const double);
int Matrealloc(Mat *, const int, const int);


double Matdot(const Mat *, const Mat *);
Mat *Matrot(const Mat *, const Mat *);
double Matdet(const Mat *);
Mat *MatLU(const Mat *, Mat *, Mat *, const char);
double ipMatinv(Mat *);
Mat *Matscal(const Mat *, const double);
int ipMatscal(Mat *, const double);
Mat *Matmul(const Mat *, const Mat *);
Mat *Matmulscal(const Mat *, const Mat *, const double);
Mat *Matmulbyvec(const Mat *, const Mat *);
Mat *Matmulscalbyvec(const Mat *, const Mat *, const double);
/*Mat *Mattransposemul(const Mat *, const Mat *);
Mat *Mattransposemulscal(const Mat *, const Mat *, const double);*/
Mat *Matsquare(const Mat *);
Mat *Matsquarescal(const Mat *, const double);
Mat *Mattransform(const Mat *, const Mat *);
Mat *Mattransformscal(const Mat *, const Mat *, const double);
Mat *Mattranspose(const Mat *);
int ipMattranspose(Mat *);
Mat *Matnull(const Mat *,const double, const char);
Mat *Matdiv(const Mat *, const Mat *, const char);

double Matmaxdouble(const Mat *);
double Matmindouble(const Mat *);
int Matmaxint(const Mat *);
int Matminint(const Mat *);

int Matdump(const Mat *);
int ipMatrowtrim(Mat *, const Mat *);
Mat *Matrowtrim(const Mat *, const Mat *);
int ipMatrowsplice(Mat *, const Mat *);
int ipMatcat(Mat *, const Mat *, const char);
int *ipMatrowsort(Mat *);
Mat *Matcut(const Mat*, const int, const int, const char);
Mat *Matread(const char *);
int Matwrite(const char *, const Mat *, const int);
Mat *Matread_mtx(const char *);
void Matwrite_mtx(const char *, const Mat *);

#endif /* MATRIXUTILS_H */

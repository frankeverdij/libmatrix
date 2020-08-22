/*
 * mat2d3d.h
 * Version: 1.1
 *
 * Copyright (C) 2012 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#ifndef MATRIX_H
#define MATRIX_H

typedef struct {
    int t;      /* type: 0=int, 1=double */
    int x;      /* number of rows in matrix */
    int y;      /* number of columns in matrix */
    union {
        int *i;
        double *d;
        void *v;
    } v;
    union {
        int **i;
        double **d;
        void **v;
    } m;
} Mat;

typedef struct {
    int t;    /* type: 0=int, 1=double */
    int x;      /* first dimension */
    int y;      /* second dimension */
    int z;      /* third dimension */
    union {
        int *i;
        double *d;
        void *v;
    } v;
    union {
        int ***i;   /* pointer to matrix values */
        double ***d;
        void ***v;
    } m;
} Mat3;

#include "matrixutils.h"

#endif /* MATRIX_H */

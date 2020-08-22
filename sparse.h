/*
 * sparse.h
 * Version: 1.1
 *
 * Copyright (C) 2012 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#ifndef SPARSE_H
#define SPARSE_H

typedef struct {
    int t;    /* type: 0=int, 1=double */
    int m;      /* number of rows in matrix */
    int n;      /* number of columns in matrix */
    int nz;     /* number of non-zero entries */
    int *r;     /* row array size nz */
    int *c;     /* column array size nz */
    int *cc;    /* compressed column array size m+1 */
    union {
        int *i;     /* matrix values size nz */
        double *d;
        void *v;
    } s;
} Sparse;

#include "sparseutils.h"

#endif /* SPARSE_H */

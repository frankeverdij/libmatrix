/*
 * arrayutils.c
 * Version: 1.1
 *
 * Copyright (C) 2012 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#include "arrayutils.h"

int *ialloc(const int x, short needzero)
{
    int *i;

    if (needzero) i = (int *)calloc((size_t)x, sizeof(int));
        else i = malloc(x*sizeof(int));
    assert(i);
    return i;
}
double *dalloc(const int x, short needzero)
{
    double *d;
    
    if (needzero) d = (double *)calloc((size_t)x, sizeof(double));
        else d = malloc(x*sizeof(double));
    assert(d);
    return d;
}

double **darrayalloc(const int r, const int c, short needzero) {
    int i;
    double *d;
    double **darray;
    
    /* we now allocate the memory for the arrays */
    d = dalloc(c*r, needzero);

    /* next we allocate room for the pointers to the rows */
    darray = malloc(r*sizeof(double *));
    assert(darray);

    /* and now we 'point' the pointers */
    for (i=0; i<r; i++) darray[i] = d + (c*i);
    
    return darray;
}

void darrayfree(double **p)
{
    if (p!=NULL) {
        free(p[0]);
        free(p);
    }
    return;
}

void arrayfree(void **p)
{
    if (p!=NULL) {
        free(p[0]);
        free(p);
    }
    return;
}

int **iarrayalloc(const int r, const int c, short needzero) {
    int i;
    int *d;
    int **iarray;
    
    /* we now allocate the memory for the arrays */
    d = ialloc(c*r, needzero);

    /* next we allocate room for the pointers to the rows */
    iarray = malloc(r*sizeof(int *));
    assert(iarray);

    /* and now we 'point' the pointers */
    for (i=0; i<r; i++) iarray[i] = d + (c*i);
    
    return iarray;
}

void iarrayfree(int **p)
{
    if (p!=NULL) {
        free(p[0]);
        free(p);
    }
    return;
}

double ***darray3Dalloc(const int x, const int y, const int z, short needzero) {
    int i, j;
    double *d;
    double ***darray;

    /* first we set aside space for the array itself */
    d = dalloc(x*y*z, needzero);

    /* next we allocate space of an array of pointers, each
       to eventually point to the first element of a
       2 dimensional array of pointers to pointers */
    darray = malloc(x*sizeof(double **));
    assert(darray);

    /* and for each of these we assign a pointer to a newly
       allocated array of pointers to a row */
    for (i = 0; i < x; i++)
    {
        darray[i] = malloc(y * sizeof(double *));
        assert(darray[i]);

        /* and for each space in this array we put a pointer to
           the first element of each row in the array space
           originally allocated */
        for (j = 0; j < y; j++) darray[i][j] = d + (i*(z * y) + j*z);
    }
    return darray;
}

void darray3Dfree(double ***p, const int x)
{
    int i;
    
    if (p!=NULL) {
        free(p[0][0]);
        if (x>0) for (i=0; i<x; i++) free(p[i]);
        free(p);
    }
    return;
}

int ***iarray3Dalloc(const int x, const int y, const int z, short needzero) {
    int i, j;
    int *d;
    int ***iarray;

    /* first we set aside space for the array itself */
    d = ialloc(x*y*z, needzero);

    /* next we allocate space of an array of pointers, each
       to eventually point to the first element of a
       2 dimensional array of pointers to pointers */
    iarray = malloc(x*sizeof(int **));
    assert(iarray);

    /* and for each of these we assign a pointer to a newly
       allocated array of pointers to a row */
    for (i = 0; i < x; i++)
    {
        iarray[i] = malloc(y * sizeof(int *));
        assert(iarray[i]);

        /* and for each space in this array we put a pointer to
           the first element of each row in the array space
           originally allocated */
        for (j = 0; j < y; j++) iarray[i][j] = d + (i*(z * y) + j*z);
    }
    return iarray;
}

void iarray3Dfree(int ***p, const int x)
{
    int i;
    
    if (p!=NULL) {
        free(p[0][0]);
        if (x>0) for (i=0; i<x; i++) free(p[i]);
        free(p);
    }
    return;
}

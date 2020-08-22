/*
 * arrayutils.h
 * Version: 1.1
 *
 * Copyright (C) 2012 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#ifndef ARRAYUTILS_H
#define ARRAYUTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int *ialloc(const int, short);
double *dalloc(const int, short);

int **iarrayalloc(const int, const int, short);
void iarrayfree(int **);

double **darrayalloc(const int, const int, short);
void darrayfree(double **);

void arrayfree(void **p);

int ***iarray3Dalloc(const int, const int, const int, short);
void iarray3Dfree(int ***, const int);

double ***darray3Dalloc(const int, const int, const int, short);
void darray3Dfree(double ***, const int);

#endif /* ARRAYUTILS_H */

/*
 * sortutils.h
 * Version: 1.1
 *
 * Copyright (C) 2012 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#ifndef SORTUTILS_H
#define SORTUTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "arrayutils.h"
#include "matrixtype.h"

#define MAXGAP 2048

int countsortsparse (int *, int *, void *, const int, const int);
void insortsparse (int *, void *, const int, const int);
void shellsortsparse(int *, void *, const int, const int);
void heapsortsparse(int *, void *, const int, const int);
int countsortint(int *, const int);
int insortint(int *, const int);
int * heapsortindex(void *, const int, const int);

#endif /* SORTUTILS_H */

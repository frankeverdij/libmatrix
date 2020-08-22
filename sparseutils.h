/*
 * sparseutils.h
 * Version: 1.2
 *
 * Copyright (C) 2016 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#ifndef SPARSEUTILS_H
#define SPARSEUTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "umfpack.h"
#include "SuiteSparse_config.h"

#include "sparse.h"
#include "mat2d3d.h"
#include "matrixtype.h"
#include "arrayutils.h"
#include "sortutils.h"
#include "mmio.h"

#ifdef UMF32
    #define UMFGLUE(f) umfpack_di_##f
    #define UMFINT int
#else
    #define UMFGLUE(f) umfpack_dl_##f
    #define UMFINT SuiteSparse_long
#endif

Sparse *Sparsealloc(const Mat *, const int);
int Sparseassign(const Mat *, Sparse *);
int Sparseinsert(Sparse *, const int, const int, const double);
int SparsesetDirichlet(Sparse *, int *, int, const short);
int Sparsesetcolumn(Sparse *, const int);
int Sparseclearrow(Sparse *, const int);
Sparse *Sparsecopy(const Sparse *);
void Sparsefree(Sparse *);
Sparse **PSparsealloc(const int);
void PSparsefree(Sparse **, const int);

double Sparseelement(const Sparse *, const int, const int);
Mat *Sparsevector(const Sparse *, const int, const char);
int *columntoCCS(int*,int,int);
int mendSparse(Sparse *, const short);
int SparseMatadd(Sparse *, const Mat *, const Mat *);
int ipSparsetrim(Sparse *, const Mat *, const char);
int ipSparsesplice(Sparse *, const Mat *, const char);
Mat *SparseMatmul(const Sparse *, const Mat *);
Mat *SparsetransposeMatmul(const Sparse *, const Mat *);
Sparse *Sparsemul(const Sparse *, const Sparse *);
Mat *Sparserowexpand(const Sparse *);

int Sparsetripletrealloc(Sparse *, const int);
int Sparsedump(const Sparse *);
void Sparsenullifypointers(Sparse *);
Mat *SparseMatdiv(const Sparse *, const Mat *);
Sparse *Sparsediv(const Sparse *, const Sparse *);
Mat *SparseMatdivfull(const Sparse *, const Mat *, void **, const short);
Sparse *Sparsedivfull(const Sparse *, const Sparse *, void **, const short);
int ipSparseadd(Sparse *, const Sparse *, const double, const short);

void aprod(const char *, int *, int *, double *, double *, double *, int *);
int ipSparsetranspose(Sparse *);
Sparse *Sparsetranspose(const Sparse *);
Mat *SparsetoMat(const Sparse *);
int ipSparsecat(Sparse *, const Sparse *, const char, const short);
Sparse *Sparseread(const char *);
int Sparsewrite(const char *, const Sparse *, const int);
Sparse *Sparseread_mtx(const char *);
void Sparsewrite_mtx(const char *, const Sparse *);
Sparse *Sparseeye(const int, const int);
Sparse *Sparsediag(const Sparse *);
void ipSparsediag(Sparse *);
int ipSparsesymgraph(Sparse *);

void ipSparseindexselect(Sparse *, const Mat *, const char);
Sparse *Sparseindexselect(const Sparse *, const Mat *, const char);

#endif /* SPARSEUTILS_H */

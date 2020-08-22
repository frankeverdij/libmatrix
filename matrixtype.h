/*
 * matrixtype.h
 * Version: 1.1
 *
 * Copyright (C) 2012 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#ifndef MATRIXTYPE_H
#define MATRIXTYPE_H

#define INTEGER 0
#define DOUBLE 1

#define MAXTYPE 1

#ifdef MX_LONG
    #define M_INT long
#else
    #define M_INT int
#endif

#endif /* MATRIXTYPE_H */

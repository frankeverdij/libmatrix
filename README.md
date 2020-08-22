Libmatrix:  an ANSI C-library which implements full and sparse matrices
 as structures of pointer arrays

Made by Frank Everdij, Delft University of Technology

Version 1.0: First release in linear FEM code.
        1.1: Added GPL License and several bugfixes and enhancements
        1.2: removed a superfluous memmove in Matrealloc for INTEGER Matrices

For use and Copyright notices, please see the COPYRIGHT file.


This library is a KISS implementation of a Matrix type useable in ANSI C.
At the time of writing these libraries, i tried to make them as fast as
possible by not relying on C++ and it's big BOOST library.
There's nothing wrong with BLAS and LAPACK for basic Matrix and Vector
manipulation, since these are thorougly tested libraries for linear algebra
problems.

I've written this library to mimic operations in Matlab(C).
In principle you can use this matrix library to port Matlab
programs into C. I've done it several times and the results are
very comparable, if not exact.

You need an implementation of the original fortran BLAS
and LAPACK version 3 libraries as well as an implementation
of UMFPACK and AMD, preferable a version 5 or above. See:
http://netlib.org/blas/ and http://netlib.org/lapack/
for releases of BLAS and LAPACK
and:
https://people.engr.tamu.edu/davis/suitesparse.html
for a release of SuiteSparse, made available by Prof. Tim Davis,
University of Texas A&M.


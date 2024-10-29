#ifndef ReflectionMethod_h
#define ReflectionMethod_h
#include <math.h>
#include "Matrix.h"

int ReflectionInverse(int n, double* matrix, double* adjoint);

int multiplyMatrixByVector(int n, double* matrix, double* b, double* result) ;

double ResidualCalc( double* a,  double* b, double* result, int n);


int LU_decomposition(int N, const double* A, double* L, double* U);

int solve_LU(int N, const double* L, const double* U, double* X, const double* B);

#endif /* ReflectionMethod_h */

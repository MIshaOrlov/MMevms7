#ifndef ReflectionMethod_h
#define ReflectionMethod_h
#include <math.h>
#include "Matrix.h"

int ReflectionInverse(int n, double* matrix, double* adjoint);

double ResidualCalc( double* a,  double* b, double* result, int n);

#endif /* ReflectionMethod_h */

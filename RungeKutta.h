//
//  RungeKutta.hpp
//  somecode8.1
//
//  Created by Михаил on 04.11.2024.
//

#ifndef RungeKutta_hpp
#define RungeKutta_hpp

#include <iostream>
#include <cmath>



double none(double x);

void RungeKutta2( double (*func)(double), double h, const double y0[], double x_vals[], double y_vals[], int N_points, int dim) ;

int ShooterMethod(double (*func)(double),int N_points,int dim,double h, double* y_vals, double* x_vals, double* u_vals,double* phi1_full,double* phi2_full,double* psi_full);

#endif /* RungeKutta_hpp */

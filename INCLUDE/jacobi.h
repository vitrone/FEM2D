/*============================================================================+/
 | Name: jacobi.h
/+============================================================================*/

#ifndef JACOBI_H
#define JACOBI_H
/************************************************************************
 *                                                                      *
 * Calculating Jacobi polynomial $J^{(a,b)}_n(x),\,a,b>-1$ and its zeros* 
 *                                                                      *
 ************************************************************************/
#include "basic.h"
#include "matlib.h"
#include "ehandler.h"

#define JCOBI_ITER_MAX (50)

void jacobi_calcjp
(
    matlib_index n, 
    double a, 
    double b, 
    double x,
    double C[n+1],
    double D[n-1],
    double E[n-1],
    double *JP                         /* Array of size 2  */ 
);


void jacobi_find_zeros
( 
    matlib_index n, 
    double a, 
    double b, 
    double *zeros, 
    double tol
);








#define CALCJP CalcJP( n, a, b, x, C1, D, E, JP)

void CalcJP
( 
    matlib_index n, 
    double a, 
    double b, 
    double x, 
    double C[n+1], 
    double D[n-1], 
    double E[n-1], 
    double *JP
);
void find_zeros( matlib_index n, double a, double b, double *zeros, double tol);
/************************************************************************
 *                                                                      *
 * Finds the $n$-th order Pade approximation of $s^{\alpha},|\alpha|<1$ *
 * in the form of partial fractions. Supply zeros of $J^{(-a,a)}_n(x)$  *
 * to compute the numerator and the denominator coefficients:           *
 * $s^{\alpha}=\beta_0+\sum_{k=0}^n\frac{\beta_k}{s+\gamma_k}$          *
 *                                                                      *
 ************************************************************************/
void diag_pade_pf
( 
    matlib_index order, 
    double frac_pow, 
    double *pfNOM_coeff, 
    double *pfDENOM_coeff, 
    double *zeros
);

#endif

/*============================================================================+/
 | Name: jacobi.c
/+============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "jacobi.h"
/*============================================================================*/
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
)
/* 
 * Calculate Jacobi polynomials: $J^{(a,b)}_{n-1}(x),J^{(a,b)}_{n}(x)$ 
 * $a,b>-1$                                                             
 *
 * */ 
{
    matlib_index i;
    double tmp;
    
    *(JP+0) = 1;
    /* 
     * C[0] = 0.5*(a-b);
     * C[1] = 0.5*(a+b+2); 
     *
     * */
     *(JP+1) = C[0]+C[1]*x;
    if(n>1) 
    {
        for( i=2; i<(n+1); i++)
        {
            /*
             * tmp = 2*i*(i+a+b)*(2*i-2+a+b);
             * C[i] = (2*i-1+a+b)*(a*a-b*b)/tmp;
             * D[i-2] = (2*i-2+a+b)*(2*i-1+a+b)*(2*i+a+b)/tmp;
             * E[i-2] = 2*(i-1+a)*(i-1+b)*(2*i+a+b)/tmp;
             * 
             * */ 
            tmp = *(JP+1);
            *(JP+1) = (C[i]+D[i-2]*x)**(JP+1)-E[i-2]**(JP+0);
            *(JP+0) = tmp;
        }
    }
}



/*============================================================================*/

void jacobi_find_zeros
( 
    matlib_index n, 
    double a, 
    double b, 
    double *zeros, 
    double tol
)
/*
 * Description: This function calculates the zeros of the Jacobi polynomial   
 * $J^{(a,b)_n(x)}$. A combination of bisection and Newton's method is used to
 * find the roots.                                                            
 *                                                                            
 * Input(s):                                                                  
 *  n   - order of the Jacobi polynomials                                     
 *  a,b - Jacobi polynomial parameters, $a,b \in[-1/2,1/2]$                   
 *  tol - numerical tolerance                                                 
 *                                                                            
 * Output(s):                                                                 
 *  zeros - roots of the Jacobi polynomial                                    
 *                                                                            
 **/

{
    double JP[2], dJP, *C, e;
    double x, xa, xb, xmid;
    double JPa, tmp;
    matlib_index i, j = 0;
    
    double theta = M_PI/(n+0.5);

    /* Perform the check on a and b parameters. */
    bool param_valid = ((fabs(a)<0.5) & (fabs(b)<0.5));
    if (!param_valid)
    {
        term_exec( "parameters outside the range (-0.5,0.5): a=%0.4f, b=%0.4f", a, b);
    }

    /* Coefficients for the recurrence relation for Jacobi polynomials */
    double C1[n+1], D[n-1], E[n-1];
    C1[0] = 0.5*(a-b);
    C1[1] = 0.5*(a+b+2);
    for( i=2; i<(n+1); i++){
        tmp = 2*i*(i+a+b)*(2*i-2+a+b);
        C1[i] = (2*i-1+a+b)*(a*a-b*b)/tmp;
        D[i-2] = (2*i-2+a+b)*(2*i-1+a+b)*(2*i+a+b)/tmp;
        E[i-2] = 2*(i-1+a)*(i-1+b)*(2*i+a+b)/tmp;
    }
    
    C = calloc(3, sizeof(double));
    tmp = 2.0*n+a+b;
    *(C+0) = n*(a-b)/tmp;
    *(C+1) = n*(2.0*n+a+b)/tmp;
    *(C+2) = 2.0*(n+a)*(n+b)/tmp;


    if(n>0)
    {
        for( i=0; i<n; i++)
        {
            /* Determine the bracketting interval */
            xa = -cos((i+0.5)*theta);
            xb = -cos((i+1)*theta);
            jacobi_calcjp( n, a, b, xa, C1, D, E, JP);
            JPa = *(JP+1);
            jacobi_calcjp( n, a, b, xb, C1, D, E, JP);
            tmp = JPa**(JP+1);
            /* exit if the bracketting interval does not contain even one root 
             * */
            if(tmp>0)
            {
                term_exec("Root does not lie in: [% 0.16f, % 0.16f]", xa, xb);
            }
            xmid = 0.5*(xa+xb);

            jacobi_calcjp( n, a, b, xmid, C1, D, E, JP);
            
            e = fabs(*(JP+1));
            if(e>tol)
            {
                tmp = JPa**(JP+1);
                if(tmp<0)
                    xb = xmid;
                else 
                {
                    xa = xmid;
                    JPa = *(JP+1);
                }
            }
            x = xmid;
            while(e>tol && j<JCOBI_ITER_MAX)
            {
                /* Newton step is only used to find a better bracketting
                 * interval. 
                 * */ 
                dJP = ((*(C+0)-x**(C+1))**(JP+1)+*(C+2)**(JP+0))/(1-x*x);
                xmid = x-*(JP+1)/dJP;
                jacobi_calcjp( n, a, b, xmid, C1, D, E, JP);
                e = fabs(*(JP+1));
                if(e>tol)
                {
                    if(xmid>xa && xmid<xb)
                    {
                        tmp = JPa**(JP+1);
                        if(tmp<0)
                            xb = xmid;
                        else 
                        {
                            xa = xmid;
                            JPa = *(JP+1);
                        }
                    } 
                    else 
                        xmid = 0.5*(xa+xb);

                    x = xmid;
                }
                j++;
            }
            if (j>JCOBI_ITER_MAX)
            {
                term_exec("Threshold iteration reached (j=%d)", j);
            }
            j = 0;
            *zeros = xmid;
            zeros++;
        }
    }
}
/*============================================================================*/

void CalcJP
(
    matlib_index n, 
    double a, 
    double b, 
    double x,
    double C[n+1],
    double D[n-1],
    double E[n-1],
    double *JP                         /* Array of size 2  */ 
)
/************************************************************************
 *                                                                      *
 * Calculate Jacobi polynomials: $J^{(a,b)}_{n-1}(x),J^{(a,b)}_{n}(x)$  *
 * $a,b>-1$                                                             *
 *                                                                      *
 ************************************************************************/
{
    matlib_index i;
    double tmp;
    
    *(JP+0) = 1;
    /* 
     * C[0] = 0.5*(a-b);
     * C[1] = 0.5*(a+b+2); 
     *
     * */
     *(JP+1) = C[0]+C[1]*x;
    if(n>1) {
        for( i=2; i<(n+1); i++){
            /*
             * tmp = 2*i*(i+a+b)*(2*i-2+a+b);
             * C[i] = (2*i-1+a+b)*(a*a-b*b)/tmp;
             * D[i-2] = (2*i-2+a+b)*(2*i-1+a+b)*(2*i+a+b)/tmp;
             * E[i-2] = 2*(i-1+a)*(i-1+b)*(2*i+a+b)/tmp;
             * 
             * */ 
            tmp = *(JP+1);
            *(JP+1) = (C[i]+D[i-2]*x)**(JP+1)-E[i-2]**(JP+0);
            *(JP+0) = tmp;
        }
    }
}

/*============================================================================*/

void find_zeros
( 
    matlib_index n, 
    double a, 
    double b, 
    double *zeros, 
    double tol
)
/******************************************************************************
 * Description: This function calculates the zeros of the Jacobi polynomial   *
 * $J^{(a,b)_n(x)}$. Bisection method is used to find the roots.              *
 *                                                                            *
 * Input(s):                                                                  *
 *  n   - order of the Jacobi polynomials                                     *
 *  a,b - Jacobi polynomial parameters, a,b>-1                                *
 *  tol - numerical tolerance                                                 *
 *                                                                            *
 * Output(s):                                                                 *
 *  zeros - roots of the Jacobi polynomial                                   *
 *                                                                            *
 ******************************************************************************/

{
    double *JP, dJP, *C, e;
    double tmp;
    matlib_index i, j;
    
    double h = 1.0/(n*n);
    double x = -1.0;
    JP = calloc( 2, sizeof(double));
   
    /* Coefficients for the recurrence relation for Jacobi polynomials */
    double C1[n+1], D[n-1], E[n-1];
    C1[0] = 0.5*(a-b);
    C1[1] = 0.5*(a+b+2);
    for( i=2; i<(n+1); i++){
        tmp = 2*i*(i+a+b)*(2*i-2+a+b);
        C1[i] = (2*i-1+a+b)*(a*a-b*b)/tmp;
        D[i-2] = (2*i-2+a+b)*(2*i-1+a+b)*(2*i+a+b)/tmp;
        E[i-2] = 2*(i-1+a)*(i-1+b)*(2*i+a+b)/tmp;
    }
    
    C = calloc(3, sizeof(double));
    tmp = 2.0*n+a+b;
    *(C+0) = n*(a-b)/tmp;
    *(C+1) = n*(2.0*n+a+b)/tmp;
    *(C+2) = 2.0*(n+a)*(n+b)/tmp;


    if(n>0){
        for( i=0; i<n; i++){
            x = x+h*0.1;
            CALCJP;
            tmp = *(JP+1);
            x = x+h;
            CALCJP;
            tmp = *(JP+1)*tmp;
            while( tmp>0 && x<1.0){
                x = x+h;
                tmp = *(JP+1);
                CALCJP;
                tmp = *(JP+1)*tmp;
            }
            if(fabs(*(JP+1))<=tol){
                *(zeros+i) = x;
            }
            else{
                x = x-h;
                CALCJP;
                if( fabs(*(JP+1))<=tol){
                    *(zeros+i) = x-h;
                }
                else{
                    x = x+h/2.0;
                    CALCJP;
                    e = fabs(*(JP+1));
                    while(e>tol){
                        dJP = ((*(C+0)-x**(C+1))**(JP+1)+*(C+2)**(JP+0))/(1-x*x);      
                        x = x-*(JP+1)/dJP;
                        CALCJP;
                        e = fabs(*(JP+1));
                    }
                    *(zeros+i) = x;
                }
            }
        }
    }
    free(C);
}
/*============================================================================*/

/*============================================================================*/

void diag_pade_pf               
/* Partial fraction resolution of diagonal Pade approximants            */ 
( 
    matlib_index n,                          /* Order of diagonal Pade approx.   */
    double frac_pow,                /* $|\alpha|<1$ in $s^{\alpha}$     */
    double *pfNUM_coeff, 
    double *pfDENOM_coeff,          /* Numerator/Denominator coeff.     * 
                                     * of the partial fractions         * 
                                     * for the diag. Pade approximants. */
    double *zeros                   /* Zeros of the Jacobi polynomial   *
                                     * $J^{(-\alpha, \alpha)}_n(x)$     */
)
{
    double *JP, tmp, tmp1;
    matlib_index i;
    
    JP = calloc( 2, sizeof(double));
    /* Coefficients for the recurrence relation for Jacobi polynomials */
    double C1[n+1], D1[n-1], E1[n-1];
    double a = frac_pow;
    double b = -a;
    C1[0] = 0.5*(a-b);
    C1[1] = 0.5*(a+b+2);
    for( i=2; i<(n+1); i++){
        tmp = 2*i*(i+a+b)*(2*i-2+a+b);
        C1[i] = (2*i-1+a+b)*(a*a-b*b)/tmp;
        D1[i-2] = (2*i-2+a+b)*(2*i-1+a+b)*(2*i+a+b)/tmp;
        E1[i-2] = 2*(i-1+a)*(i-1+b)*(2*i+a+b)/tmp;
    }
    
    CalcJP( n, frac_pow, -frac_pow, 1.0, C1, D1, E1, JP);
    tmp = *(JP+1);
    /* Coefficients for the recurrence relation for Jacobi polynomials */
    a = -frac_pow;
    b = -a;
    double C2[n+1], D2[n-1], E2[n-1];
    C2[0] = 0.5*(a-b);
    C2[1] = 0.5*(a+b+2);
    for( i=2; i<(n+1); i++){
        tmp1 = 2*i*(i+a+b)*(2*i-2+a+b);
        C2[i] = (2*i-1+a+b)*(a*a-b*b)/tmp1;
        D2[i-2] = (2*i-2+a+b)*(2*i-1+a+b)*(2*i+a+b)/tmp1;
        E2[i-2] = 2*(i-1+a)*(i-1+b)*(2*i+a+b)/tmp1;
    }
    CalcJP( n, -frac_pow, frac_pow, 1.0, C2, D2, E2, JP);
    *(pfNUM_coeff+0) = tmp/(*(JP+1));

    if(n<=0){
        printf("diagonal Pade approximation: the order must be >0!\n");
        exit(EXIT_FAILURE);
    }
    else{
        for ( i=0; i<n; i++){
            *(pfDENOM_coeff+i) = (1.0+*(zeros+i))/(1.0-*(zeros+i));
            CalcJP( n, frac_pow,-frac_pow, *(zeros+i), C1, D1, E1, JP);
            tmp = *(JP+1);
            CalcJP( n, -frac_pow, frac_pow, *(zeros+i), C2, D2, E2, JP);
            *(pfNUM_coeff+i+1) = -(2.0*n/(n*n-frac_pow*frac_pow))
                                  **(pfDENOM_coeff+i)*tmp/(*(JP+0));
        }
    }
    free(JP);
}

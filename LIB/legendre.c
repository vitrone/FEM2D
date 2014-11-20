/*============================================================================+/
 | Name   : legendre.c                                                   
 | Author : Vishal Vaibhav                                               
 |                                                                       
 | Description : This module defines all functions for implementation of 
 | Legendre transforms.                                                  
 |                                                                       
 | History : Creation 21 July 2013                                       
/+============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


#define NDEBUG
#define MATLIB_NTRACE_DATA
#include "legendre.h"
#include "assert.h"

/*============================================================================*/

matlib_real Legendre_LGL (matlib_real x, matlib_index n)
/* 
 * Description: This function calculates the value of the LGL polynomials at a
 * given point x.                                                             
 * LGL polynomial is define as $(1-x^2)L'_n(x)$.                              
 *                  $(1-x^2)L'_n(x)=nL_{n-1}(x)-nxL_n(x)$                     
 *                                                                            
 * Input(s) :                                                                 
 *  x - coordinate point                                                      
 *  n - degree of the Legendre polynomial                                     
 *                                                                            
 * Output(s): value of the LGL polynomial                                     
 *
 * */
{
    assert(n>2);
    assert(fabs(x)<1.0);

    matlib_real LP[2], tmp, r;
    matlib_index i; 
    
    LP[0] = 1.0;
    LP[1] = x;
    for(i=2; i<(n+1); i++)
    {
        tmp = LP[1];
        /* $nL_n(x)=(2n-1)xL_{n-1}(x)-(n-1)L_{n-2}(x)$                  
         * */ 
        LP[1] = ((2*i-1)*x*LP[1]-(i-1)*LP[0])/i;
        LP[0] = tmp;
    }
    r = n**(LP+0)-n*x**(LP+1);
    return r;
}
/*============================================================================*/

void CalcLP
( 
    matlib_index n,        
    matlib_real *LP,    
    matlib_real x, 
    matlib_real C[n-1],
    matlib_real D[n-1]
)
/* 
 *
 * Description: This function calculates the value of the Legendre polynomials
 * $L_n$ and $L_{n-1}$.                                                       
 *                                                                            
 * Input(s) :                                                                 
 *  n - degree of the Legendre polynomial                                     
 *  x - coordinate pomatlib_index                                               
 *  C and D - are the array of coefficients that shows up in the recurrence     
 *  relation of the Legendre polynomials (precomputed coeff. speeds up the      
 *  computation).                                                             
 *   $L_n(x)=[(2n-1)/n]xL_{n-1}(x)-[(n-1)/n]L_{n-2}(x)$                         
 *                                                                              
 * Output(s):                                                                   
 *  LP - pointer to the array [L_{n-1}, L_{n}]                                  
 *
 * */
{

    matlib_real tmp;
    matlib_index i;
    
    *(LP+0) = 1.0;
    *(LP+1) = x;
    if(n>1) 
    {
        tmp = *(LP+1);
        *(LP+1) = *C*x**(LP+1)-*D**(LP+0);
        *(LP+0) = tmp;
        for(i=3; i<(n+1); i++)
        {
            C++; D++;
            tmp = *(LP+1);
            *(LP+1) = *C*x**(LP+1)-*D**(LP+0);
            *(LP+0) = tmp;
        }
    }
}
/*============================================================================*/

void find_Gauss_points
( 
    matlib_index n, 
    matlib_real  tol,   
    matlib_real* zeros
)
/* 
 * Description: This function computes the zeros of the Legendre polnomial or 
 * Gauss points. The bracketing interval of the zeros are known to be         
 * $x_j\in[-\cos((j-0.5)\pi/(n+1/2)), -\cos(j\pi/(n+1/2))]$.                  
 *                                                                            
 * Input(s) :                                                                 
 *  n   - degree of Legendre polynomial                                       
 *  tol - tolerance of error                                                  
 *                                                                            
 * Output(s):                                                                 
 *  zeros - pointer to the array of Gauss-zeros                               
 *
 * */
{
    matlib_real LP[2], x, e;
    matlib_real tmp, LPa, xa, xb, xmid;
    matlib_index i, j = 0;
    matlib_real theta = M_PI/(n+0.5);
    matlib_real C[n-1], D[n-1];
    /* 
     * C and D are coefficients in the recurrence relation of the Legendre
     * polynomials
     * */ 
    for ( i=0; i<n-1; i++)
    {
        C[i] = (2*i+3.0)/(i+2);
        D[i] = (i+1.0)/(i+2);
    }


    if(n>1){
        for(i=1; i<n+1; i++){
            /* begin Newton iterations once the zero has been bracketed */ 
            xa = -cos((i-0.5)*theta);
            xb = -cos(i*theta);
            CalcLP( n, LP, xa, C, D);
            LPa = *(LP+1);
            CalcLP( n, LP, xb, C, D);
            tmp = LPa**(LP+1);
            if(tmp>0)
            {
                exit(EXIT_FAILURE);
            }
            xmid = 0.5*(xa+xb);
            CalcLP( n, LP, xmid, C, D);
            e = fabs(*(LP+1));
            if(e>tol)
            {
                tmp = LPa**(LP+1);
                if(tmp<0)
                {
                    xb = xmid;
                } 
                else 
                {
                    xa = xmid;
                    LPa = *(LP+1);
                }
            }
            x = xmid;
            while(e>tol && j<1e4)
            {
                /* Newton step is only used to find a better bracketting
                 * interval. 
                 * */ 
                xmid = x+(1.0-x*x)**(LP+1)/(n*x**(LP+1)-n**(LP+0));
                CalcLP( n, LP, xmid, C, D);
                e = fabs(*(LP+1));
                if(e>tol)
                {
                    if(xmid>xa && xmid<xb)
                    {
                        tmp = LPa**(LP+1);
                        if(tmp<0)
                        {
                            xb = xmid;
                        } 
                        else 
                        {
                            xa = xmid;
                            LPa = *(LP+1);
                        }
                    } 
                    else 
                    {
                        xmid = 0.5*(xa+xb);
                    }
                    x = xmid;
                }
                j++;
            }
            j = 0;
            *zeros = xmid;
            zeros++;
        }
    }
}
/*============================================================================*/

void find_LGL_points
( 
    matlib_index n, 
    matlib_real tol, 
    matlib_real* zeros, 
    matlib_real* quadW,
    matlib_real *gzeros
)
/*
 * Description: This function finds the zeros of LGL polynomial defined as    
 * $(1-x^2)L'_n$. The zeros are known as LGL-points. The Gauss zeros bracket  
 * the LGL zeros. This information is provided in order to speed up the       
 * computation. Newton's iterations is given by                                          
 * $x_{k+1} = \frac{n}{n+1}x_k+\frac{L_{n-1}(x_k)}{(n+1)L_n(x_k)}$.           
 *                                                                            
 * Input(s) :                                                                 
 *  n   - degree of Legendre polynomial                                       
 *  tol - tolerance of error                                                  
 *  gzeros - Gauss zeros                                                      
 *                                                                            
 * Output(s):                                                                 
 *  zeros - LGL-points                                                        
 *  quadW - quadrature weights for weight function $\omega(x)=1$.             
 *  The weights are given by $\omega_j = \frac{2}{n(n+1)L_n^2(x_j)}$.         
 *
 * */
{
    debug_enter("degree of Legendre polynomials: %d", n);
    matlib_real LP[2], x, e, tmp1;
    matlib_real tmp, LGLa, xa, xb, xmid;
    matlib_index i, j = 0;
    matlib_real gamma;
    matlib_real C[n-1], D[n-1];

    /* 
     * C and D are coefficients in the recurrence relation of the Legendre
     * polynomials
     * */ 
    for ( i=0; i<n-1; i++){
        C[i] = (2*i+3.0)/(i+2);
        D[i] = (i+1.0)/(i+2);
    }
    
    gamma = 2.0/(n*(n+1.0));
    *(zeros+0) = -1.0;
    *(zeros+n) = 1.0;
    
    *(quadW+0) = gamma;
    *(quadW+n) = gamma;

    zeros++;
    quadW++;

    if(n>1){
        for(i=1; i<n; i++){
            /* root is bracketed in [xa, xb]
             * */ 
            xa = *(gzeros+i-1);
            xb = *(gzeros+i);

            CalcLP( n, LP, xa, C, D);
            LGLa = n**(LP+0)-n*xa**(LP+1);
            CalcLP( n, LP, xb, C, D);
            tmp = LGLa*(n**(LP+0)-n*xb**(LP+1));
            if(tmp>0)
                exit(EXIT_FAILURE);

            xmid = 0.5*(xa+xb);
            CalcLP( n, LP, xmid, C, D);
            tmp = n**(LP+0)-n*xmid**(LP+1);
            e = fabs(tmp);
            if(e>tol){
                tmp1 = LGLa*tmp;
                if(tmp1<0){
                    xb = xmid;
                } 
                else {
                    xa = xmid;
                    LGLa = tmp;
                }
            }
            /* A combination of bisection and Newton method is used to find the
             * root. Newton's update is used to find new bracketing interval
             * instead of the midpomatlib_index as in the pure bisection method.
             * */
            x = xmid;
            while(e>tol && j<1e4){
                /* xmid here is the outcome of a Newton step
                 * */
                xmid = n*x/(n+1)+*(LP+0)/(*(LP+1)*(n+1));
                CalcLP( n, LP, xmid, C, D);
                tmp = n**(LP+0)-n*xmid**(LP+1);
                e = fabs(tmp);
                if(e>tol){
                    if(xmid>xa && xmid<xb){
                        tmp1 = LGLa*tmp;
                        if(tmp1<0){
                            xb = xmid;
                        } 
                        else {
                            xa = xmid;
                            LGLa = tmp;
                        }
                    } 
                    else {
                        xmid = 0.5*(xa+xb);
                    }
                    x = xmid;
                }
                j++;
            }
            j = 0;
            *zeros = xmid;
            zeros++;
            *quadW = gamma/(*(LP+1)**(LP+1));
            quadW++;
        }
    }
    
    debug_exit("%s", "");
}
/*============================================================================+/
 | Following subroutines use matrices stored in row major format.
/+============================================================================*/
void backward_transform_matrix
(
    const matlib_index P,            
    const matlib_index p,            
    const matlib_real* x,          
          matlib_real* pILTM
)
/* 
 * Description: This function computes the inverse transformation matrix for 
 * Legendre transforms with LGL points. Maximum degree of Legendre polynomial
 * used for the interpolation is 'p'. The values of the transformed function 
 * is computed at a set of points 'x'.                                       
 *
 * */
{
    
    debug_enter("degree of polynomial: %d, nr. sampling points: %d", p, P+1);

    matlib_index i, j;
    /* ILTM: matrix in row major format */ 
    matlib_index stILTM = p+1;
    matlib_real tmp1, tmp2;
    /* pILTM is a pointer to an array of length p+1, note that the paranthesis 
     * is necessary!
     * */ 
    
    bool s = (p>1);

    /* Filling the first two columns */ 
    for (i=0; i<P+1; i++)
    {
        *(pILTM+i*stILTM+0) = 1.0;
        *(pILTM+i*stILTM+1) = *(x+i);
    }
    if(s)
    {
        for(j=2; j<(p+1); j++)
        {
            tmp1 = (2.0*j-1.0)/j;
            tmp2 = (j-1.0)/j;
            for(i=0; i<P+1; i++)
            {
                *(pILTM+i*stILTM+j) =  tmp1**(x+i)**(pILTM+i*stILTM+(j-1))
                                      -tmp2**(pILTM+i*stILTM+(j-2));
            }
        }
    }
    debug_exit("%s", "");

}
/*============================================================================*/

void forward_transform_matrix
(                                                                 
    const matlib_index   p,            
    const matlib_real* zeros,      
          matlib_real* pFLTM
)
/* 
 * Description: This function computes the forward transformation matrix for  
 * Legendre transforms with LGL points. Maximum degree of Legendre polynomial 
 * used for the interpolation is 'p'.                                         
 *                                                                            
 * Input(s) :                                                                 
 *  p     - maximum of the degree of Legendre polynomials used                
 *  zeros - p+1 LGL points                                                    
 *                                                                            
 * Output(s):                                                                 
 *  pFLTM - is a pointer to an array pointing to the first row of FLTM matrix 
 *  which has the dimensions (p+1)-by-(p+1) stored in a row major format.     
 *
 * */
{
    debug_enter("degree of polynomial: %d", p);

    matlib_index i, j;
    
    matlib_real LP[p+1], *pLP = &LP[0];
    matlib_real C[p-1], D[p-1], *pC = &C[0], *pD = &D[0];
    matlib_real E[p+1], *pE = &E[0];
    bool s = (p>1);

    matlib_real g = 1.0/(p*(p+1));
    matlib_index stFLTM = (p+1);

    /* Fill the first column 
     * */ 
    (*(pFLTM+0)) = g; 
    *pE = 1.0;
    pE++;
    for(i=1; i<p; i++)
    {
        g   = -g;
        *pE = 2*i+1.0;
        *(pFLTM+i*stFLTM) = *pE*g;
        pE++;
    }
    g = -g;
    *(pFLTM + p*stFLTM) = p*g;
    /* Fill the last column
     * */ 
    pE = &E[0];
    g  = 1.0/(p*(p+1));
    for(i=0; i<p; i++)
    {
        (*(pFLTM+i*stFLTM+p)) = *pE*g;
        pE++;
    }
    *pE = 1.0/(p+1);
    (*(pFLTM+p*stFLTM+p)) = *pE;
    

    if(s)
    {
        matlib_real tmp_zeros = 0;
        *pLP = 1.0;
        pLP++;
        *pLP = *(zeros+1);
        tmp_zeros = *pLP;
        
        for(j=2; j<(p+1); j++)
        {
            pLP++;
            *pC = (2*j-1.0)/j;
            *pD = (j-1.0)/j;
            *pLP = *pC*tmp_zeros**(pLP-1)-*pD**(pLP-2);
            pC++;
            pD++;
        }
        g = 1.0/(p*(p+1.0)**pLP**pLP);
        /* Fill the 2nd column
         * */ 
        i = 1;
        pE  = &E[0];
        pLP = &LP[0];
        for (j=0; j<p; j++)
        {
            (*(pFLTM+j*stFLTM+i)) = *pE*g**pLP;
            pE++;
            pLP++;
        }
        (*(pFLTM+p*stFLTM+i)) = *pE/(*pLP); 

        for (i=2; i<p; i++)
        {
            pC = &C[0];
            pD = &D[0];

            pLP = &LP[0];
            pLP++;
            *pLP = *(zeros+i);
            tmp_zeros = *pLP;
            for(j=2; j<(p+1); j++)
            {
                pLP++;
                *pLP = *pC*tmp_zeros**(pLP-1)-*pD**(pLP-2);
                pC++;
                pD++;
            }
            g = 1.0/(p*(p+1.0)**pLP**pLP);
            /* Fill the ith column
             * */ 
            pE  = &E[0];
            pLP = &LP[0];
            for (j=0; j<p; j++)
            {
                (*(pFLTM+j*stFLTM+i)) = *pE*g**pLP;
                pE++;
                pLP++;
            }
            (*(pFLTM+p*stFLTM+i)) = *pE/(*pLP); 
        }
    }
    debug_exit("degree f polynomial: %d", p);
}
/*============================================================================*/

void forward_transform_matrix2
(
    const matlib_index p,                        
    const matlib_index P,                        
    const matlib_real* zeros,                  
          matlib_real* pFLTM                   
)
/* 
 * Description: This subroutine computes the transformation matrix for        
 * Legendre transformation where the function is sampled on (P+1) LGL points  
 * while the projection is onto a Legendre polynomial-space of degree p.      
 *                                                                            
 * Input(s) :                                                                 
 *  p     - maximum of the degree of Legendre polynomials used for the final  
 *  interpolation                                                             
 *  P     - P+1 LGL points are used to sample the function, note that P would 
 *  be the highest degree of the Legendre polynomial involved                 
 *  zeros - P+1 LGL points                                                    
 *                                                                            
 * Output(s):                                                                 
 *  pFLTM - is a pointer to an array pointing to the first row of the matrix  
 *  FLTM with dimensions (p+1)-by-(P+1). It is stored in a row major format.  
 *
 * */
{

    debug_enter("degree of polynomial: %d, nr. LGL sampling points: %d", p, P+1);
    matlib_index i, j;
    
    /* If you want the integer to be interpreted as a float, use decimal 
     * notation i.e. instead of 1 use 1.0. This function was gving wrong results 
     * because of this error.
     * */


    /* p<=P
     *
     * If this is not satisfied, set p=P.
     * */

    if (p>P)
    {
        term_exec("%s","Degree of Legendre polynomials is incorrect");
    }
    bool s  = (p<P);
    bool s1 = (P>1);
    
    matlib_real LP[P+1], *pLP = &LP[0];
    matlib_real C[P-1], D[P-1], *pC = &C[0], *pD = &D[0];
    matlib_real E[p+1], *pE = &E[0];
    matlib_index stFLTM = (P+1);

    matlib_real g = 1.0/(P*(P+1.0));
    
    /* Fill the first column
     * */ 
    (*(pFLTM+0)) = g;
    *pE = 1.0;
    pE++;
    for(i=1; i<p; i++)
    {
        g = -g;
        *pE = 2.0*i+1.0;
        (*(pFLTM+i*stFLTM+0)) = *pE*g;
        pE++;
    }
    
    *pE = 2.0*p+1.0;
    g = -g;
    if (s)
    {
        (*(pFLTM+p*stFLTM+0)) = *pE*g;
    }
    else
    {
        (*(pFLTM+p*stFLTM+0)) = P*g;
    }
    /* Fill the last column
     * */
    g = 1.0/(P*(P+1.0));
    pE = &E[0];
    for(i=0; i<p; i++)
    {
        (*(pFLTM+i*stFLTM+P)) = *pE*g;
        pE++;
    }
    matlib_real tmp = 1.0/(P+1.0); 
    if(s)
        (*(pFLTM+p*stFLTM+P)) = *pE*g;
    else
        (*(pFLTM+p*stFLTM+P)) = tmp;
    
    if(s1)
    {
        /* 
         * This part of the code generates the legendre polynomials.
         * All degree polynomials values are stored.
         * */
        matlib_real tmp_zeros = 0;

        *pLP = 1.0;
        pLP++;
        *pLP = *(zeros+1);
        tmp_zeros = *pLP;
        
        for(j=2; j<(P+1); j++)
        {
            pLP++;
            *pC = (2*j-1.0)/j;
            *pD = (j-1.0)/j;
            *pLP = *pC*tmp_zeros**(pLP-1)-*pD**(pLP-2);
            pC++;
            pD++;
        }
        g = 1.0/(P*(P+1.0)**pLP**pLP);
        
        /* Fill the 2nd column
         * */
        pE  = &E[0];
        pLP = &LP[0];
        i = 1;
        for (j=0; j<p; j++)
        {
            (*(pFLTM+j*stFLTM+i)) = *pE*g**pLP;
            pE++;
            pLP++;
        }
        if (s)
            (*(pFLTM+p*stFLTM+i)) = *pE*g**pLP;
        else
            (*(pFLTM+p*stFLTM+i)) = tmp/(*pLP);
            
        for (i=2; i<P; i++)
        {
            /* 
             * Generating legendre polynomials
             *
             * */
            pC  = &C[0];
            pD  = &D[0];
            pLP = &LP[0];

            *pLP = 1.0;
            pLP++;
            *pLP = *(zeros+i);
            tmp_zeros = *pLP;
            for(j=2; j<(P+1); j++)
            {
                pLP++;
                *pLP = *pC*tmp_zeros**(pLP-1)-*pD**(pLP-2);
                pC++;
                pD++;
            }
            g = 1.0/(P*(P+1.0)**pLP**pLP);
            
            /* Fill the ith column
             * */
            pE  = &E[0];
            pLP = &LP[0];
            for (j=0; j<p; j++)
            {
                (*(pFLTM+j*stFLTM+i)) = *pE*g**pLP;
                pE++;
                pLP++;
            }
            if (s)
                (*(pFLTM+p*stFLTM+i)) = *pE*g**pLP;
            else
                (*(pFLTM+p*stFLTM+i)) = tmp/(*pLP);
        }
    }
    debug_exit("%s","");
}
/*============================================================================+/
 |Column major version of the above functions for MATLAB
/+============================================================================*/
void backward_transform_matrix_colmajor
( 
    matlib_index P, 
    matlib_index p, 
    matlib_real *x, 
    matlib_real *pILTM
)
/*
 * Description: This function computes the inverse transformation matrix for  
 * Legendre transforms with LGL points. Maximum degree of Legendre polynomial 
 * used for the interpolation is 'p'. The values of the transformed function  
 * is computed at a set of points 'x'.                                        
 *                                                                            
 * Input(s) :                                                                 
 *  P - size of x minus 1                                                     
 *  p - maximum of the degree of Legendre polynomials used                    
 *  x - pointer to the array of points where the interpolation polynomial is  
 *  evaluated, length of the array is P+1                                     
 *                                                                            
 * Output(s):                                                                 
 * pILTM: points to the first element of a (P+1)-by-(p+1) matrix stored       
 * in a column major format                                                   
 *                ILTM[i][j] = *(pILTM+i+j*col)                               
 *                                                                            
 * */
{

    matlib_index i, j;
    matlib_real tmp1, tmp2;
    bool s = (p>1);
    matlib_index col = P+1;
    matlib_index col1 = col;
    matlib_index col2 = 0;

    /* Filling the first two columns */ 
    for (i=0; i<P+1; i++)
    {
        *(pILTM+i) = 1.0;
        *(pILTM+col) = *(x+i);
        col++;
    }
    /*                                                                  
     * $nL_n = (2n-1)xL_{n-1}-(n-1)L_{n-2}$                             
     *                                                                  
     * */
    if(s)
    {
        for( j=2; j<p+1; j++)
        {
            tmp1 = (2.0*j-1.0)/j;
            tmp2 = (j-1.0)/j;
            for ( i=0; i<P+1; i++)
            {
                *(pILTM+col) = tmp1**(x+i)**(pILTM+col1)
                               -tmp2**(pILTM+col2);
                col++; 
                col1++;
                col2++;
            }
        }
    }
}
/*======================================================================*/

void forward_transform_matrix_colmajor
( 
    matlib_index p, 
    matlib_real *zeros, 
    matlib_real *pFLTM
)

/*
 * Description: This function computes the forward transformation matrix for  
 * Legendre transforms with LGL points. Maximum degree of Legendre polynomial 
 * used for the interpolation is 'p'.                                         
 *                                                                            
 * Input(s) :                                                                 
 *  p     - maximum of the degree of Legendre polynomials used                    
 *  zeros - p+1 LGL points                                                    
 *                                                                            
 * Output(s):                                                                 
 * pFLTM - pointer to the first element of FLTM                               
 * FLTM is a (p+1)-by-(p+1) matrix stored in a column major format            
 *                FLTM[i][j] = *(pFLTM+i+j*col)                               
 *                                                                            
 * */
{
    matlib_index i, j;
    matlib_real LP[p+1], *pLP = &LP[0];
    matlib_real C[p-1], D[p-1], *pC = &C[0], *pD = &D[0];
    matlib_real E[p+1], *pE = &E[0];
    bool s = (p>1);

    matlib_real g = 1.0/(p*(p+1));
    
    /* Fill the first column */ 
    *(pFLTM+0) = g; 
    *pE = 1.0;
    pE++;
    for(j=1; j<p; j++)
    {
        g   = -g;
        *pE = (2*j+1.0);
        *(pFLTM+j) = *pE*g;
        pE++;
    }
    g =   -g;
    *(pFLTM+p) = p*g;
    
    /* Fill the last column */ 
    matlib_index col = p*(p+1);
    pE = &E[0];
    g  = 1.0/(p*(p+1));
    for( j=0; j<p; j++)
    {
        *(pFLTM+col) = *pE*g;
        col++;
        pE++;
    }
    *pE = 1.0/(p+1);
    *(pFLTM+col) = *pE;

    if(s)
    {
        matlib_real tmp_zeros = 0;
        col = p+1;
        *pLP = 1.0;
        pLP++;
        *pLP = *(zeros+1);
        tmp_zeros = *pLP;

        for( j=2; j<(p+1); j++)
        {
            pLP++;
            *pC = (2.0*j-1.0)/j;
            *pD = (j-1.0)/j;
            *pLP = *pC*tmp_zeros**(pLP-1)-*pD**(pLP-2);
            pC++;
            pD++;
        }
        g = 1.0/(p*(p+1.0)**pLP**pLP);
        /* Fill the 2nd column */
        pE  = &E[0];
        pLP = &LP[0];
        for (j=0; j<p; j++)
        {
            *(pFLTM+col)= *pE*g**pLP;
            col++;
            pE++;
            pLP++;
        }
        *(pFLTM+col) = *pE/(*pLP); 
        col++;
        for ( i=2; i<p; i++)
        {
            pC = &C[0];
            pD = &D[0];

            pLP = &LP[0];
            pLP++;
            *pLP = *(zeros+i);
            tmp_zeros = *pLP;

            for( j=2; j<(p+1); j++)
            {
                pLP++;
                *pLP = *pC*tmp_zeros**(pLP-1)-*pD**(pLP-2);
                pC++;
                pD++;
            }
            g = 1.0/(p*(p+1.0)**pLP**pLP);
            /* Fill the ith column */
            pE  = &E[0];
            pLP = &LP[0];
            for (j=0; j<p; j++)
            {
                *(pFLTM+col)= *pE*g**pLP;
                col++;
                pE++;
                pLP++;
            }
            *(pFLTM+col) = *pE/(*pLP); 
            col++;
        }
    }
}
/*======================================================================*/

void forward_transform_matrix2_colmajor
( 
    matlib_index p,                      
    matlib_index P,                      
    matlib_real  *zeros,                  
    matlib_real  *pFLTM                   
)
/* 
 * Description: This subroutine computes the transformation matrix for        
 * Legendre transformation where the function is sampled on (P+1) LGL points  
 * while the projection is onto a Legendre polynomial-space of degree p.      
 *                                                                            
 * Input(s) :                                                                 
 *  p     - maximum of the degree of Legendre polynomials used for the final  
 *  interpolation                                                             
 *  P     - P+1 LGL points are used to sample the function, note that P would 
 *  be the highest degree of the Legendre polynomial involved                 
 *  zeros - P+1 LGL points                                                    
 *                                                                            
 * Output(s):                                                                 
 *  pFLTM - is a pointer to the first element of the matrix FLTM with         
 *  dimensions (p+1)-by-(P+1). It is stored in a column major format.         
 *                FLTM[i][j] = *(pFLTM+i+j*col)                               
 *                                                                            
 * */
{

    matlib_index i, j;
    
    if (p>P)
    {
        term_exec("%s","Degree of Legendre polynomials is incorrect");
    }
    bool s  = (p<P);
    bool s1 = (P>1);

    matlib_real LP[P+1], *pLP = &LP[0];
    matlib_real C[P-1], D[P-1], *pC = &C[0], *pD = &D[0];
    matlib_real E[p+1], *pE = &E[0];

    matlib_real g = 1.0/(P*(P+1));
    /* Fill the first column  */ 
    *(pFLTM+0) = g; 
    *pE = 1.0;
    pE++;
    for(j=1; j<p; j++)
    {
        g = -g;
        *pE = 2*j+1.0;
        *(pFLTM+j) = *pE*g;
        pE++;
    }
    *pE = 2.0*p+1.0;
    g = -g;
    if (s)
    {
        *(pFLTM+p) = *pE*g;
    }
    else
    {
        *(pFLTM+p) = P*g;
    }

    /* Fill the last column  */ 
    matlib_index col = P*(p+1);
    g = 1.0/(P*(P+1));
    pE = &E[0];
    for( j=0; j<p; j++)
    {
        *(pFLTM+col) = *pE*g;
        col++;
        pE++;
    }
    matlib_real tmp = 1.0/(P+1);
    if (s)
        *(pFLTM+col) = *pE*g;
    else
        *(pFLTM+col) = tmp;
    
    
    if(s1)
    {
        matlib_real tmp_zeros = 0;
        col = p+1;
        *pLP = 1.0;
        pLP++;
        *pLP = *(zeros+1);
        tmp_zeros = *pLP;

        for( j=2; j<(P+1); j++)
        {
            pLP++;
            *pC = (2.0*j-1.0)/j;
            *pD = (j-1.0)/j;
            *pLP = *pC*tmp_zeros**(pLP-1)-*pD**(pLP-2);
            pC++;
            pD++;
        }
        g = 1.0/(P*(P+1.0)**pLP**pLP);

        /* Fill the 2nd column 
         * */ 
        pE  = &E[0];
        pLP = &LP[0];
        for ( j=0; j<p; j++)
        {
            *(pFLTM+col) = *pE*g**pLP;
            col++;
            pE++;
            pLP++;
        }
        if (s)
            *(pFLTM+col) = *pE*g**pLP;
        else
            *(pFLTM+col) = tmp/(*pLP);
        col++;

        for ( i=2; i<P; i++)
        {
            pC  = &C[0];
            pD  = &D[0];
            pLP = &LP[0];

            *pLP = 1.0;
            pLP++;
            *pLP = *(zeros+i);
            tmp_zeros = *pLP;
            for(j=2; j<(P+1); j++)
            {
                pLP++;
                *pLP = *pC*tmp_zeros**(pLP-1)-*pD**(pLP-2);
                pC++;
                pD++;
            }
            g = 1.0/(P*(P+1.0)**pLP**pLP);

            /* Fill the ith column */
            pE  = &E[0];
            pLP = &LP[0];
            for (j=0; j<p; j++)
            {
                *(pFLTM+col) = *pE*g**pLP;
                col++;
                pE++;
                pLP++;
            }
            if (s)
                *(pFLTM+col) = *pE*g**pLP;
            else
                *(pFLTM+col) = tmp/(*pLP);
            col++;
        }
    }

}
/*============================================================================+/
 | Utility functions with matrix/vector objects
/+============================================================================*/

void legendre_LGLdataLT1
( 
    const matlib_index p, 
    const matlib_real  tol,
          matlib_xv*   zeros,
          matlib_xv*   quadW
)
/* 
 * zeros  - array of length p+1 containing LGL-points
 * quadW  - array of length p+1 containing quadrature weights
 *
 * */ 
{
    debug_enter("degree of polynomial: %d, tolerance: %0.16g", p, tol);
    
    matlib_index i;
    matlib_real gzeros[p];
    matlib_index nr_LGL = p+1;

    matlib_create_xv(nr_LGL, zeros, MATLIB_COL_VECT);
    matlib_create_xv(nr_LGL, quadW, MATLIB_COL_VECT);


    find_Gauss_points( p, tol, gzeros);

    BEGIN_DTRACE
        for (i=0; i<p; i++)
        {
            debug_body("Gauss zeros[%d]: %0.16f ", i, *(gzeros+i));
        }
    END_DTRACE

    find_LGL_points( p, tol, zeros->elem_p, quadW->elem_p, gzeros);

    debug_exit("%s", "");
}
/*============================================================================*/

void legendre_LGLdataLT2
( 
    const matlib_index      p, 
    const matlib_real     tol,
          matlib_xv* zeros,
          matlib_xv* quadW,
          matlib_xm* FM,
          matlib_xm* IM
)
/* 
 * zeros  - array of length p+1 containing LGL-points
 * quadW  - array of length p+1 containing quadrature weights
 * IM, FM - (p+1)-by-(p+1) matrices associated with backward and forward
 *          legendre transformation
 *
 * */ 
{
    debug_enter("degree of polynomial: %d, tolerance: %0.16g", p, tol);
    
    matlib_index i;
    matlib_real gzeros[p];
    matlib_index nr_LGL = p+1;
    
    matlib_create_xm(nr_LGL, nr_LGL, FM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    matlib_create_xm(nr_LGL, nr_LGL, IM, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);

    matlib_create_xv(nr_LGL, zeros, MATLIB_COL_VECT);
    matlib_create_xv(nr_LGL, quadW, MATLIB_COL_VECT);


    find_Gauss_points( p, tol, gzeros);

    BEGIN_DTRACE
        for (i=0; i<p; i++)
        {
            debug_body("Gauss zeros[%d]: %0.16f ", i, *(gzeros+i));
        }
    END_DTRACE

    find_LGL_points( p, tol, zeros->elem_p, quadW->elem_p, gzeros);
    
    backward_transform_matrix( p, p, zeros->elem_p, IM->elem_p);
    forward_transform_matrix( p, zeros->elem_p, FM->elem_p);

    debug_exit("%s", "");
}

/*============================================================================*/

void legendre_LGLdataFM
( 
    const matlib_xv xi,
          matlib_xm FM
)
/* 
 * FM: Forward transform matrix 
 * xi: LGL-points for sampling the function
 * FM.lenr = xi.len
 * */

{
    debug_enter( "nr. of LGL-points for sampling "
                 "size of FM: %d-by-%d", xi.len,
                 FM.lenc, FM.lenr );
    
    matlib_index p = FM.lenc-1;

    debug_body( "highest degree of polynomials: %d", p);

    assert((xi.elem_p != NULL) && (FM.elem_p != NULL));
    
    if(FM.lenr == xi.len)
    {
        if(FM.order == MATLIB_ROW_MAJOR)
        {
            forward_transform_matrix2( p, xi.len-1, xi.elem_p, FM.elem_p);
        }
        else if (FM.order==MATLIB_COL_MAJOR)
        {
            forward_transform_matrix2_colmajor(p, xi.len-1, xi.elem_p, FM.elem_p);
        }
        else
        {
            term_execb( "Storage order of FM unknown (order: %d)",
                        FM.order);
        }
    }
    else
    {
        term_execb( "dimension of matrix/vector incorrect: matrix "
                    "FM: %d-by-%d, vector xi: %d", 
                    FM.lenc, FM.lenr, xi.len );
    }
    
    debug_exit("%s", "");
}
/*============================================================================*/

void legendre_LGLdataIM
( 
    const matlib_xv xi,
          matlib_xm IM
)
/* 
 * IM: inverse transform matrix
 * xi: points at which value of the function is sought
 * 
 * IM.lenc = xi.len
 *
 * */

{
    debug_enter( "nr. of LGL-points for sampling "
                 "size of IM: %d-by-%d", xi.len,
                 IM.lenc, IM.lenr );
    matlib_index p = IM.lenr-1;

    debug_body( "highest degree of polynomials: %d", p);

    assert((xi.elem_p != NULL) && (IM.elem_p != NULL));
    
    if(IM.lenc == xi.len)
    {
        if(IM.order==MATLIB_ROW_MAJOR)
        {
            backward_transform_matrix( xi.len-1, p, xi.elem_p, IM.elem_p);
        }
        else if (IM.order==MATLIB_COL_MAJOR)
        {
            backward_transform_matrix_colmajor( xi.len-1, p, xi.elem_p, IM.elem_p);
        }
        else
        {
            term_execb( "Storage order of FM unknown (order: %d)",
                        IM.order);
        }
    }
    else
    {
        term_execb( "dimension of matrix/vector incorrect: matrix "
                    "IM: %d-by-%d, vector xi: %d", 
                    IM.lenc, IM.lenr, xi.len );
    }
    
    debug_exit("%s", "");
}
/*============================================================================*/

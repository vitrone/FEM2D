/*============================================================================+/
 | mxLGLdataLT1.c
/+============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "legendre.h"

#define TOL_DEFAULT (1e-6)

void mexFunction
(
          int      nlhs, 
          mxArray* plhs[],
          int      nrhs, 
    const mxArray* prhs[]
)
{

    if(nrhs!=2) 
    {
        mexErrMsgIdAndTxt("LGLdataLT1:nrhs","Two inputs required.");
    }
    if(nlhs!=2) 
    {
        mexErrMsgIdAndTxt("LGLdataLT1:nlhs","Two outputs required.");
    } 
    
    matlib_index p;
    matlib_real tol = (matlib_real )mxGetScalar(prhs[1]);

    if(mxGetScalar(prhs[0])>1)
    {
        p = (matlib_index)floor(mxGetScalar(prhs[0]));
    }
    else
    {
        mexErrMsgIdAndTxt( "FEM1D:mxLGLdataLT1:N", "Polynomial degree must be>1.");
    } 

    if((tol<0) || (tol>1))
    {
        mexWarnMsgIdAndTxt( "LGLdataLT1:Tol",
                            "Tolerance must be positive real number <1"
                            "(proceeding with default).");
        tol = TOL_DEFAULT;
        
    }

    plhs[0] = mxCreateDoubleMatrix( p+1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix( p+1, 1, mxREAL);

    double *zeros = mxGetPr(plhs[0]);
    double *quadW = mxGetPr(plhs[1]);
    double gzeros[p];

    find_Gauss_points( p, tol, gzeros);
    find_LGL_points( p, tol, zeros, quadW, gzeros);

}

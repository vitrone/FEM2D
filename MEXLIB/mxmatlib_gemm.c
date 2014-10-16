#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "fem2d.h"

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
        mexErrMsgTxt("Two inputs required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgTxt("One output required.");
    } 
    if(mxIsComplex(prhs[0]) && mxIsComplex(prhs[1]))
    {
        mexErrMsgTxt("Real matrices expected!");
    }
    matlib_xm A = { .lenc = mxGetM(prhs[0]), 
                    .lenr = mxGetN(prhs[0]), 
                    .elem_p = mxGetPr(prhs[0]),
                    .op     = MATLIB_NO_TRANS,
                    .order  = MATLIB_COL_MAJOR };

    matlib_xm B = { .lenc = mxGetM(prhs[1]), 
                    .lenr = mxGetN(prhs[1]), 
                    .elem_p = mxGetPr(prhs[1]),
                    .op     = MATLIB_NO_TRANS,
                    .order  = MATLIB_COL_MAJOR };


    if (A.lenr!=B.lenc)
    {
        mexErrMsgTxt("Dimension mis-match!");
    } 
    plhs[0] = mxCreateDoubleMatrix( A.lenc, B.lenr, mxREAL);

    matlib_xm C = { .lenc = A.lenc, 
                    .lenr = B.lenr, 
                    .elem_p = mxGetPr(plhs[0]),
                    .op     = MATLIB_NO_TRANS,
                    .order  = MATLIB_COL_MAJOR };

    matlib_xgemm(1.0, A, B, 0.0, C);

}

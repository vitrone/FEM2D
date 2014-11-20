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
    if (nrhs < 2) 
    {
        mexErrMsgTxt("At least two inputs required.");
    }
    if (nrhs > 3) 
    {
        mexErrMsgTxt("Too many input arguments.");
    }

    matlib_real*  plist = mxGetPr(prhs[0]);
    matlib_index* tlist = (matlib_index*) mxGetPr(prhs[1]);

    if(mxGetM(prhs[0]) != 2)
    {
        mexErrMsgTxt("Size of point list is incorrect.");
    } 
    if(mxGetM(prhs[1]) != 3)
    {
        mexErrMsgTxt("Size of triangle list is incorrect.");
    } 

    matlib_index nr_nodes   = mxGetN(prhs[0]);
    matlib_index nr_domains = mxGetN(prhs[1]);

    matlib_real tol = MATLIB_TOL;
    if (nrhs == 3)
    {
        tol = mxGetScalar(prhs[2]);
    }

    if(nr_nodes < 3)
    {
        mexErrMsgTxt("At least three nodes are needed.");
    } 
    if(nr_domains < 1)
    {
        mexErrMsgTxt("At least one triangle is needed.");
    } 

    fem2d_cc nodes = { .len    = nr_nodes, 
                       .elem_p = plist };
    fem2d_ea ea;
    fem2d_err error = fem2d_create_ea( nodes, tlist, nr_domains, &ea);
    if (error == FEM2D_FAILURE)
    {
        mexErrMsgTxt("Failed to create element array.");
    }

    const char *fnames[6] = { "len", "tol", "h_max", 
                              "angle_max", "area", "aspect_ratio"};
    matlib_index nr_fields = 6;                             
    plhs[0] = mxCreateStructMatrix( 1, 1, nr_fields, fnames);
    
    mxArray* mxA[6];
    mxA[0] = mxCreateNumericMatrix( 1, 1, mxUINT64_CLASS, mxREAL);
    mxA[1] = mxCreateDoubleMatrix(          1, 1, mxREAL);
    mxA[2] = mxCreateDoubleMatrix( nr_domains, 1, mxREAL);
    mxA[3] = mxCreateDoubleMatrix( nr_domains, 1, mxREAL);
    mxA[4] = mxCreateDoubleMatrix( nr_domains, 1, mxREAL);
    mxA[5] = mxCreateDoubleMatrix( nr_domains, 1, mxREAL);

    *((matlib_index*)mxGetPr(mxA[0])) = nr_domains;
    *(mxGetPr(mxA[1])) = tol;

    fem2d_mq mq = { .len          = nr_domains, 
                    .tol          = tol, 
                    .h_max        = mxGetPr(mxA[2]),
                    .angle_max    = mxGetPr(mxA[3]),
                    .area         = mxGetPr(mxA[4]),
                    .aspect_ratio = mxGetPr(mxA[5])  };

    error = fem2d_mesh_quality( ea, &mq);
    if (error == FEM2D_FAILURE)
    {
        mexErrMsgTxt("Failed to compute mesh quality.");
    }

    mxSetField(plhs[0], 0, "len", mxA[0]);
    mxSetField(plhs[0], 0, "tol", mxA[1]);

    mxSetField(plhs[0], 0, "h_max"       , mxA[2]);
    mxSetField(plhs[0], 0, "angle_max"   , mxA[3]);
    mxSetField(plhs[0], 0, "area"        , mxA[4]);
    mxSetField(plhs[0], 0, "aspect_ratio", mxA[5]);

    fem2d_free_ea(ea);
}

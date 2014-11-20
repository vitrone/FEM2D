#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"

#include "fem2d.h"

#define POINT_ENUM_VAL(penum)         \
    (penum == FEM2D_BOUNDARY? 1: 0)
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

    matlib_real*  plist = mxGetPr(prhs[0]);
    matlib_index* tlist = (matlib_index*) mxGetPr(prhs[1]);

    if (mxGetM(prhs[0]) != 2)
    {
        mexErrMsgTxt("Size of point list is incorrect.");
    } 
    if (mxGetM(prhs[1]) != 3)
    {
        mexErrMsgTxt("Size of triangle list is incorrect.");
    } 

    matlib_index nr_nodes  = mxGetN(prhs[0]);
    matlib_index nr_domains = mxGetN(prhs[1]);
    if (nr_nodes < 3)
    {
        mexErrMsgTxt("At least three nodes are needed.");
    } 
    if (nr_domains < 1)
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
    error = fem2d_create_vp(&ea);
    if (error == FEM2D_FAILURE)
    {
        mexErrMsgTxt("Failed to create vertex patch.");
    }


    const char *fnames[6] = { "node_index", "len", "point_enum", "vert_index", 
                              "domain_index", "node_order"};
    matlib_index nr_fields = 6;                             
    plhs[0] = mxCreateStructMatrix( nr_nodes, 1, nr_fields, fnames);
    
    matlib_index i, j, fnum;
    mxArray *mxNI, *mxlen, *mxVI, *mxpoint_enum, *mxDI, *mxNO;
    matlib_index *vptr, *dptr;
    matlib_int *nptr;

    fem2d_vp* vp_ptr;
    
    const matlib_index POS_NODE_INDEX   = mxGetFieldNumber(plhs[0], "node_index"  ); 
    const matlib_index POS_LEN          = mxGetFieldNumber(plhs[0], "len"         ); 
    const matlib_index POS_POINT_ENUM   = mxGetFieldNumber(plhs[0], "point_enum"  ); 
    const matlib_index POS_VERT_INDEX   = mxGetFieldNumber(plhs[0], "vert_index"  ); 
    const matlib_index POS_DOMAIN_INDEX = mxGetFieldNumber(plhs[0], "domain_index"); 
    const matlib_index POS_NODE_ORDER   = mxGetFieldNumber(plhs[0], "node_order"  ); 

    for (i = 0; i < ea.nr_nodes; i++)
    {
        vp_ptr = (ea.vpatch_p) + i;

        mxNI = mxCreateNumericMatrix( 1, 1, mxUINT64_CLASS, mxREAL);
        *((matlib_index*)mxGetPr(mxNI)) = (vp_ptr->node_index + 1);
        mxSetFieldByNumber(plhs[0], i, POS_NODE_INDEX, mxNI);

        mxlen = mxCreateNumericMatrix( 1, 1, mxUINT64_CLASS, mxREAL);
        *((matlib_index*)mxGetPr(mxlen)) = (vp_ptr->len);
        mxSetFieldByNumber(plhs[0], i, POS_LEN, mxlen);

        mxpoint_enum = mxCreateNumericMatrix( 1, 1, mxUINT64_CLASS, mxREAL);
        *((matlib_index*)mxGetPr(mxpoint_enum)) = POINT_ENUM_VAL(vp_ptr->point_enum);
        mxSetFieldByNumber(plhs[0], i, POS_POINT_ENUM, mxpoint_enum);

        mxVI = mxCreateNumericMatrix( vp_ptr->len, 1, mxUINT64_CLASS, mxREAL);
        mxDI = mxCreateNumericMatrix( vp_ptr->len, 1, mxUINT64_CLASS, mxREAL);
        vptr = (matlib_index*)mxGetPr(mxVI);
        dptr = (matlib_index*)mxGetPr(mxDI);
        for ( j = 0; j < vp_ptr->len; j++)
        {
            vptr[j] = (vp_ptr->vert_index[j] + 1);
            dptr[j] = ((vp_ptr->domain_p[j])->domain_index + 1);
        }
        mxSetFieldByNumber(plhs[0], i, POS_VERT_INDEX, mxVI);
        mxSetFieldByNumber(plhs[0], i, POS_DOMAIN_INDEX, mxDI);

        if (vp_ptr->point_enum == FEM2D_BOUNDARY)
        {
            mxNO = mxCreateNumericMatrix( vp_ptr->len + 1, 1, mxINT64_CLASS, mxREAL);
            nptr = (matlib_int*)mxGetPr(mxNO);
            for ( j = 0; j < vp_ptr->len +1; j++)
            {
                nptr[j] = vp_ptr->node_order[j];
            }
        }
        else
        {
            mxNO = mxCreateNumericMatrix( vp_ptr->len, 1, mxINT64_CLASS, mxREAL);
            nptr = (matlib_int*)mxGetPr(mxNO);
            for ( j = 0; j < vp_ptr->len; j++)
            {
                nptr[j] = vp_ptr->node_order[j];
            }
        }
        mxSetFieldByNumber(plhs[0], i, POS_NODE_ORDER, mxNO);
    }

    fem2d_free_ea(ea);
}

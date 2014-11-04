#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"

#define NDEBUG
#define MATLIB_NTRACE_DATA

#include "matlib.h"
/*============================================================================*/

static void adapt_complex_data
(
    matlib_index len,
    matlib_real* pr, 
    matlib_real* pi,
    matlib_complex* pz
)
{
    debug_enter("length: %d", len);
    for(matlib_index i = 0; i < len; i++)
    {
        pz[i] = pr[i]+I*pi[i];
    }
    debug_exit("%s", "");
}
/*============================================================================*/

static void fill_data
(
    matlib_io_t* mp, 
    matlib_index i,
    const mxArray* prhs,
    matlib_complex** cptr,
    matlib_index* zcnt
)
{
    debug_enter( "data index: %d, z-count: %d", i, *zcnt);
    matlib_xm *xm_ptr, xM;
    matlib_xv *xv_ptr, xV;
    matlib_zm *zm_ptr, zM;
    matlib_zv *zv_ptr, zV;
    matlib_nv *nv_ptr, nV;
    matlib_index lenc = mxGetM(prhs); 
    matlib_index lenr = mxGetN(prhs);
    matlib_index len =  lenc * lenr;


    if (mxIsComplex(prhs))
    {
        matlib_real* pr = mxGetPr(prhs); 
        matlib_real* pi = mxGetPi(prhs);

        if ((lenc == 1) || (lenr == 1))
        {
            mp->data_p[i] = calloc(1, sizeof(matlib_zv));
            zv_ptr        = (matlib_zv*) mp->data_p[i];
            zv_ptr->len   = len; 
            cptr[*zcnt]   = calloc(len, sizeof(matlib_complex));
            zv_ptr->elem_p = cptr[*zcnt];
            adapt_complex_data( len, pr, pi, zv_ptr->elem_p );

            mp->format[i] = MATLIB_ZV;
            zV = *zv_ptr;
            DEBUG_PRINT_ZV(zV, "Input nr: %d, ", i);
            (*zcnt) ++;
        }
        else
        {
            mp->data_p[i] = calloc(1, sizeof(matlib_zm));
            zm_ptr        = (matlib_zm*) mp->data_p[i];
            zm_ptr->lenc  = lenc;
            zm_ptr->lenr  = lenr;
            cptr[*zcnt]   = calloc( len, sizeof(matlib_complex));
            zm_ptr->elem_p = cptr[*zcnt];
            adapt_complex_data( len, pr, pi, zm_ptr->elem_p );

            zm_ptr->op    = MATLIB_NO_TRANS;
            zm_ptr->order = MATLIB_COL_MAJOR ;
            mp->format[i] = MATLIB_DEN_ZM;
            zM = *zm_ptr;
            DEBUG_PRINT_ZM(zM, "Input nr: %d, ", i);
            (*zcnt) ++;
        }
    }
    else if (mxIsDouble(prhs))
    {
        if ((lenc == 1) || (lenr == 1))
        {
            mp->data_p[i]  = calloc(1, sizeof(matlib_xv));
            xv_ptr         = (matlib_xv*) mp->data_p[i];
            xv_ptr->len    = len; 
            xv_ptr->elem_p = mxGetPr(prhs);
            mp->format[i]  = MATLIB_XV;
            xV = *xv_ptr;
            DEBUG_PRINT_XV(xV, "Input nr: %d, ", i);
        }
        else
        {
            mp->data_p[i]  = calloc(1, sizeof(matlib_xm));
            xm_ptr         = (matlib_xm*) mp->data_p[i];
            xm_ptr->lenc   = lenc;
            xm_ptr->lenr   = lenr;
            xm_ptr->elem_p = mxGetPr(prhs);

            xm_ptr->op    = MATLIB_NO_TRANS;
            xm_ptr->order = MATLIB_COL_MAJOR ;
            mp->format[i] = MATLIB_DEN_XM;
            xM = *xm_ptr;
            DEBUG_PRINT_XM(xM, "Input nr: %d, ", i);
        }
    }
    else if (mxIsInt64(prhs))
    {
        mp->data_p[i]  = calloc(1, sizeof(matlib_nv));
        nv_ptr         = (matlib_nv*) mp->data_p[i];
        nv_ptr->len    = len;
        nv_ptr->elem_p = (matlib_int*)mxGetPr(prhs);
        mp->format[i]  = MATLIB_NV;
        nV = *nv_ptr;
        DEBUG_PRINT_NV(nV, "Input nr: %d, ", i);
    }
    else
    {
        mexErrMsgTxt("Unsupported data format.");
    }
}
/*============================================================================*/

void mexFunction
(
          int      nlhs, 
          mxArray* plhs[],
          int      nrhs, 
    const mxArray* prhs[]
)
{


    if(nrhs<2) 
    {
        mexErrMsgTxt("At least two inputs required.");
    }
    if (!mxIsChar(prhs[0]) || (mxGetM(prhs[0]) != 1 ) )  
    {
        mexErrMsgTxt("The first input argument must be a string.");
    } 
    if(nlhs != 0) 
    {
        mexErrMsgTxt("Too many output arguments.");
    }

    char* file_name;
    size_t str_len;

    str_len   = mxGetN(prhs[0])*sizeof(mxChar)+1;
    file_name = mxMalloc(str_len);
    mxGetString(prhs[0], file_name, (mwSize)str_len);   
    mexPrintf("Saving data in file:  %s\n", file_name);
   
    matlib_io_t mp;
    matlib_err error;

    if(mxIsCell(prhs[1]))
    {
        matlib_index cell_len = mxGetNumberOfElements(prhs[1]);
        mxArray* ca[cell_len];

        error = matlib_io_create( cell_len, &mp);
        if (error == MATLIB_FAILURE)
        {
            mxFree(file_name);
            mexErrMsgTxt("Failed to initialize file data.");
        }
        matlib_complex* cptr[cell_len];
        matlib_index zcnt = 0;
        matlib_index i;
        for (i = 0; i < cell_len; i++)
        {
            ca[i] = mxGetCell(prhs[1], i);
            fill_data(&mp, i, ca[i], cptr, &zcnt);
        }
        
        error = matlib_io_fwrite(&mp, file_name);
        if (error == MATLIB_FAILURE)
        {
            matlib_io_free(&mp);
            mxFree(file_name);
            mexErrMsgTxt("Failed to write the file.");
        }
        
        for (i = 0; i < zcnt; i++ )
        {
            matlib_free(cptr[i]);
        }
    } /* Input data: cell array of matrices */ 
    else
    {
        matlib_index nr_mat = (matlib_index)(nrhs - 1);
        debug_body("Number of matrices/vectors: %d", nr_mat);

        error = matlib_io_create( nr_mat, &mp);
        if (error == MATLIB_FAILURE)
        {
            mxFree(file_name);
            mexErrMsgTxt("Failed to initialize file data.");
        }
        matlib_complex* cptr[nr_mat];
        matlib_index zcnt = 0;
        matlib_index i;
        for ( i = 0; i < nr_mat; i++)
        {
            fill_data(&mp, i, prhs[i+1], cptr, &zcnt);
        }

        error = matlib_io_fwrite(&mp, file_name);
        if (error == MATLIB_FAILURE)
        {
            matlib_io_free(&mp);
            mxFree(file_name);
            mexErrMsgTxt("Failed to write the file.");
        }

        for (i = 0; i < zcnt; i++ )
        {
            matlib_free(cptr[i]);
        }
    }
    matlib_io_free(&mp);
    mxFree(file_name);
}

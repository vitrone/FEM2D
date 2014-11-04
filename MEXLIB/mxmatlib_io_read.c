#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include "mex.h"
#include "matrix.h"

#define NDEBUG
#define MATLIB_NTRACE_DATA

#include "matlib.h"
/*============================================================================*/
typedef struct
{
    matlib_index index;
    mxArray* data_p;
    /* Following fields ensure that the data-types can be adapted to that
     * supported in MATLB 
     * */ 
    matlib_complex** zdata_p;
    matlib_index zcnt;

} mxData_t;


static matlib_err mxCreate_nv
(
    matlib_index length,
    matlib_nv*   v,
    void* extdata_p 
)
{
    v->len  = length;
    mxArray* mxptr = mxCreateNumericMatrix( length, 1, mxINT64_CLASS, mxREAL);
    err_check( (mxptr == NULL), clean_up, 
               "Memory allocation for integer vector of length %d failed!", 
               length);
    matlib_index index = ((mxData_t*)extdata_p)->index;
    mxSetCell(((mxData_t*)extdata_p)->data_p, index, mxptr);

    v->elem_p = (matlib_int*)mxGetPr(mxptr);
    
    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;

}
static matlib_err mxCreate_zv
(
    matlib_index  length,
    matlib_zv*    v,
    MATLIB_VECT_T type_enum,
    void* extdata_p
)
{
    v->len  = length;
    v->type = type_enum;

    errno = 0;
    v->elem_p = calloc( length, sizeof(matlib_complex));
    err_check( (v->elem_p == NULL), clean_up, 
               "%s: Memory allocation for complex vector of length %d failed!", 
               strerror(errno), length);
    matlib_index zcnt = ((mxData_t*)extdata_p)->zcnt;
    ((mxData_t*)extdata_p)->zdata_p[zcnt] = v->elem_p;
    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;

}
static matlib_err mxCreate_xv
(
    matlib_index  length,
    matlib_xv*    v,
    MATLIB_VECT_T type_enum,
    void* extdata_p
)
{
    v->len  = length;
    v->type = type_enum;

    mxArray* mxptr = mxCreateDoubleMatrix( 1, length, mxREAL);
    err_check( (mxptr == NULL), clean_up, 
               "Memory allocation for complex vector of length %d failed!", 
               length);
    matlib_index index = ((mxData_t*)extdata_p)->index;
    mxSetCell(((mxData_t*)extdata_p)->data_p, index, mxptr);
    v->elem_p = mxGetPr(mxptr);
    
    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}
/*============================================================================*/

static matlib_err mxCreate_zm
( 
    matlib_index lenc,
    matlib_index lenr,
    matlib_zm*   M,
    MATLIB_ORDER order_enum,
    MATLIB_OP    op_enum,
    void* extdata_p
)
{
    debug_enter("size: %d-by-%d", lenc, lenr);
    M->lenc = lenc;
    M->lenr = lenr;
    M->order = order_enum;
    M->op    = op_enum;

    errno = 0;
    M->elem_p = calloc( lenc * lenr, sizeof(matlib_complex));
    err_check( (M->elem_p == NULL), clean_up, 
               "%s: Memory allocation for complex matrix of size %d-by-%d failed!", 
               strerror(errno), lenc, lenr);
    matlib_index zcnt = ((mxData_t*)extdata_p)->zcnt;
    ((mxData_t*)extdata_p)->zdata_p[zcnt] = M->elem_p;

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

static matlib_err mxCreate_xm
( 
    matlib_index lenc,
    matlib_index lenr,
    matlib_xm*   M,
    MATLIB_ORDER order_enum,
    MATLIB_OP    op_enum,
    void* extdata_p
)
{
    M->lenc  = lenc;
    M->lenr  = lenr;
    M->order = order_enum;
    M->op    = op_enum;
    
    mxArray* mxptr;

    if (order_enum == MATLIB_COL_MAJOR)
    {
        mxptr = mxCreateDoubleMatrix( lenc, lenr, mxREAL);
    }
    else if (order_enum == MATLIB_ROW_MAJOR)
    {
        mxptr = mxCreateDoubleMatrix( lenr, lenc, mxREAL);
    }

    err_check( (mxptr == NULL), clean_up, 
               "Memory allocation for real matrix of size %d-by-%d failed!", 
               lenc, lenr);
    matlib_index index = ((mxData_t*)extdata_p)->index;
    mxSetCell(((mxData_t*)extdata_p)->data_p, index, mxptr);
    M->elem_p = mxGetPr(mxptr);
    
    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
    
}
/*============================================================================*/

static void mxDestroy(void* extdata_p)
{
    matlib_index index = ((mxData_t*)extdata_p)->index;
    mxArray* mxptr     = mxGetCell(((mxData_t*)extdata_p)->data_p, index);
    mxDestroyArray(mxptr);
}

static void mxDestroyZ(void* extdata_p)
{
    matlib_index zcnt = ((mxData_t*)extdata_p)->zcnt;
    matlib_complex* cptr = ((mxData_t*)extdata_p)->zdata_p[zcnt];
    matlib_free(cptr);
}

/*============================================================================*/

static void adapt_complex_data
(
    matlib_index lenc,
    matlib_index lenr,
    MATLIB_ORDER order_enum,
    matlib_complex* pz,
    matlib_real* pr, 
    matlib_real* pi
)
{

    debug_enter("size: %d-by-%d", lenc, lenr);
    if (order_enum == MATLIB_COL_MAJOR)
    {
        matlib_index len = lenc * lenr;
        for(matlib_index i = 0; i < len; i++)
        {
            pr[i] = creal(pz[i]);
            pi[i] = cimag(pz[i]);
        }
    }
    else if (order_enum == MATLIB_ROW_MAJOR)
    {
        matlib_complex tmp;
        for (matlib_index i = 0; i < lenc; i++)
        {
            for (matlib_index j = 0; j < lenr; j++)
            {
                tmp = pz[i*lenr+j];

                pr[i+j*lenc] = creal(tmp);
                pi[i+j*lenc] = cimag(tmp);
            }
        }
    }
    debug_exit("%s", "");
}
/*============================================================================*/
static void RowMajor2ColMajor
( 
    matlib_index lenc, 
    matlib_index lenr,
    matlib_real *ptr, 
    matlib_real *tptr
)
{
    matlib_real tmp;
    for (matlib_index i = 0; i < lenc; i++)
    {
        for (matlib_index j = 0; j < lenr; j++)
        {
            tmp = ptr[i*lenr+j];

            tptr[i+j*lenc] = tmp;
        }
    }
}

/*============================================================================*/

void clean_up(matlib_complex** zptr, matlib_index zcnt)
{
    for (matlib_index i = 0; i < zcnt; i++ )
    {
        matlib_free(zptr[i]);
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
    if(nrhs != 1)
    {
        mexErrMsgTxt("At least two inputs required.");
    }
    if (!mxIsChar(prhs[0]) || (mxGetM(prhs[0]) != 1 ) )  
    {
        mexErrMsgTxt("The first input argument must be a string.");
    } 
    if(nlhs < 1) 
    {
        mexErrMsgTxt("Too few output arguments.");
    }
    if(nlhs > 1) 
    {
        mexErrMsgTxt("Too many output arguments.");
    }
    char* file_name;
    size_t str_len;

    str_len   = mxGetN(prhs[0])*sizeof(mxChar)+1;
    file_name = mxMalloc(str_len);
    mxGetString(prhs[0], file_name, (mwSize)str_len);   
    mexPrintf("Reading data from the file:  %s\n", file_name);
   
    matlib_io_t mp;
    matlib_err error;
    mxArray* mxptr;

    errno = 0;
    FILE* fp;
    fp = fopen(file_name, "rb");
    if ( fp == NULL)
    {
        mexErrMsgTxt("Failed to open the file!");
    }


    error = matlib_io_getinfo(&mp, fp);
    if (error == MATLIB_FAILURE)
    {
        mxFree(file_name);
        mexErrMsgTxt("Failed to get file info.");
    }
    matlib_index i, tmp_len;
    matlib_int j;

    i = fseek(fp, mp.size[0], SEEK_SET);
    if ( i != 0)
    {
        mxFree(file_name);
        fclose(fp);
        mexErrMsgTxt("Error occurred moving within the file!");
    }

    plhs[0] = mxCreateCellMatrix(mp.len, 1);
    matlib_complex* zptr[mp.len];
    mxData_t extdata = { .data_p  = plhs[0], 
                         .index   = 0, 
                         .zdata_p = zptr, 
                         .zcnt = 0};

    for (i = 0; i < mp.len; i++)
    {
        extdata.index = i;
        if (mp.format[i] == MATLIB_XV)
        {
            error = matlib_io_xv_extread(&mp, i, fp, mxCreate_xv, mxDestroy, &extdata);
            if (error == MATLIB_FAILURE)
            {
                mxFree(file_name);
                fclose(fp);
                clean_up(zptr, extdata.zcnt);
                mexErrMsgTxt("Failed to read real vector!");
            }
        }
        else if (mp.format[i] == MATLIB_ZV)
        {
            error = matlib_io_zv_extread(&mp, i, fp, mxCreate_zv, mxDestroy, &extdata);
            if (error == MATLIB_FAILURE)
            {
                mxFree(file_name);
                fclose(fp);
                clean_up(zptr, extdata.zcnt);
                mexErrMsgTxt("Failed to read complex vector!");
            }
            mxptr = mxCreateDoubleMatrix(((matlib_zv*)mp.data_p[i])->len, 1, mxCOMPLEX);
            adapt_complex_data( ((matlib_zv*)mp.data_p[i])->len, 
                                1,
                                MATLIB_COL_MAJOR,
                                extdata.zdata_p[extdata.zcnt], 
                                mxGetPr(mxptr), 
                                mxGetPi(mxptr));
            mxSetCell(extdata.data_p, extdata.index, mxptr);
            (extdata.zcnt) ++;

        }
        else if (mp.format[i] == MATLIB_DEN_XM)
        {
            error = matlib_io_xm_extread(&mp, i, fp, mxCreate_xm, mxDestroy, &extdata);
            if (error == MATLIB_FAILURE)
            {
                mxFree(file_name);
                fclose(fp);
                clean_up(zptr, extdata.zcnt);
                mexErrMsgTxt("Failed to read real matrix!");
            }
            
            if (((matlib_xm*)mp.data_p[i])->order == MATLIB_ROW_MAJOR)
            {
                mxptr = mxCreateDoubleMatrix( ((matlib_xm*)mp.data_p[i])->lenc, 
                                              ((matlib_xm*)mp.data_p[i])->lenr, mxREAL);
                RowMajor2ColMajor( ((matlib_xm*)mp.data_p[i])->lenc, 
                                   ((matlib_xm*)mp.data_p[i])->lenr, 
                                   ((matlib_xm*)mp.data_p[i])->elem_p,
                                   mxGetPr(mxptr));
                mxDestroyArray(mxGetCell(extdata.data_p, extdata.index));
                mxSetCell(extdata.data_p, extdata.index, mxptr);
            }
            
        }
        else if (mp.format[i] == MATLIB_DEN_ZM)
        {
            error = matlib_io_zm_extread(&mp, i, fp, mxCreate_zm, mxDestroyZ, &extdata);
            if (error == MATLIB_FAILURE)
            {
                mxFree(file_name);
                fclose(fp);
                clean_up(zptr, extdata.zcnt);
                mexErrMsgTxt("Failed to read complex matrix!");
            }
            mxptr = mxCreateDoubleMatrix( ((matlib_zm*)mp.data_p[i])->lenc, 
                                          ((matlib_zm*)mp.data_p[i])->lenr, mxCOMPLEX);

            adapt_complex_data( ((matlib_zm*)mp.data_p[i])->lenc,
                                ((matlib_zm*)mp.data_p[i])->lenr,
                                ((matlib_zm*)mp.data_p[i])->order,
                                extdata.zdata_p[extdata.zcnt], 
                                mxGetPr(mxptr), 
                                mxGetPi(mxptr));

            mxSetCell(extdata.data_p, extdata.index, mxptr);
            (extdata.zcnt) ++;

        }
        else if (mp.format[i] == MATLIB_NV)
        {
            error = matlib_io_nv_extread(&mp, i, fp, mxCreate_nv, mxDestroy, &extdata);
            if (error == MATLIB_FAILURE)
            {
                mxFree(file_name);
                fclose(fp);
                clean_up(zptr, extdata.zcnt);
                mexErrMsgTxt("Failed to integer vector!");
            }
        }
        else
        {
           mp.format[i] == MATLIB_FORMAT_UNKNOWN;
        }
        if (mp.format[i] == MATLIB_FORMAT_UNKNOWN)
        {
            mxFree(file_name);
            fclose(fp);
            clean_up(zptr, extdata.zcnt);
            mexErrMsgTxt("Unsupported format encountered!");
        }
    }

    fclose(fp);
    matlib_io_free(&mp);
    clean_up(zptr, extdata.zcnt);
}



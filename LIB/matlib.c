/*============================================================================+
 | Name   : matlib.c                                                          |
 | Author : Vishal Vaibhav                                                    |
 |                                                                            |
 |                                                                            |
 | History : Created 20 Feb 2014                                              |
 |                                                                            |
 |                                                                            |
 |                                                                            |
 |                                                                            |
 +============================================================================*/

/*============================================================================+
 |                                                                            |
 |                              Global Includes                               |
 |                                                                            |
 +============================================================================*/
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>

/*============================================================================+
 |                                                                            |
 |                         Own Module Includes                                |
 |                                                                            |
 +============================================================================*/
//#define NDEBUG

#include "matlib.h"
#include "assert.h"
/*============================================================================+
 |                                                                            |
 |                    External Module Includes                                |
 |                                                                            |
 +============================================================================*/
#include "mkl.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

/*============================================================================+/
 |Allocation of memory
/+============================================================================*/
matlib_err matlib_create_nv
(
    matlib_index length,
    matlib_nv*   v
)
{
    v->len  = length;
    errno = 0;
    v->elem_p = calloc( length, sizeof(matlib_int));
    err_check( (v->elem_p == NULL), clean_up, 
               "%s: Memory allocation for integer vector of length %d failed!", 
               strerror(errno), length);
    
    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;

}


matlib_err matlib_create_zv
(
    matlib_index  length,
    matlib_zv*    v,
    MATLIB_VECT_T type_enum
)
{
    v->len  = length;
    v->type = type_enum;

    errno = 0;
    v->elem_p = calloc( length, sizeof(matlib_complex));
    err_check( (v->elem_p == NULL), clean_up, 
               "%s: Memory allocation for complex vector of length %d failed!", 
               strerror(errno), length);
    
    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;

}

matlib_err matlib_create_xv
(
    matlib_index  length,
    matlib_xv*    v,
    MATLIB_VECT_T type_enum
)
{
    v->len  = length;
    v->type = type_enum;

    errno = 0;
    v->elem_p = calloc( length, sizeof(matlib_real));
    err_check( (v->elem_p == NULL), clean_up, 
               "%s: Memory allocation for real vector of length %d failed!", 
               strerror(errno), length);
    
    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}
/*============================================================================*/
matlib_err matlib_create_zm
( 
    matlib_index lenc,
    matlib_index lenr,
    matlib_zm*   M,
    MATLIB_ORDER order_enum,
    MATLIB_OP    op_enum
)
{
    M->lenc = lenc;
    M->lenr = lenr;
    M->order = order_enum;
    M->op    = op_enum;

    errno = 0;
    M->elem_p = calloc( lenc * lenr, sizeof(matlib_complex));
    err_check( (M->elem_p == NULL), clean_up, 
               "%s: Memory allocation for complex matrix of size %d-by-%d failed!", 
               strerror(errno), lenc, lenr);
    
    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

matlib_err matlib_create_xm
( 
    matlib_index lenc,
    matlib_index lenr,
    matlib_xm*   M,
    MATLIB_ORDER order_enum,
    MATLIB_OP    op_enum
)
{
    M->lenc = lenc;
    M->lenr = lenr;
    M->order = order_enum;
    M->op    = op_enum;

    errno  = 0;
    M->elem_p = calloc( lenc * lenr, sizeof(matlib_real));
    err_check( (M->elem_p == NULL), clean_up, 
               "%s: Memory allocation for real matrix of size %d-by-%d failed!", 
               strerror(errno), lenc, lenr);
    
    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
    
}
/*============================================================================+/
 | Utility Functions
 +============================================================================*/
/* 
 * A complex vector of length N can be considered as a matrix of size N-by-2 in
 * a row major format. This implies accessing imaginary part of a vector 
 * requires stride 2.
 * C(i,j) = *(pC+2*i+j)
 *
 * */ 

/*============================================================================+/
 |BLAS Level I Routines
/+============================================================================*/
matlib_err matlib_xcopy
(
    const matlib_xv x,
          matlib_xv y
)
{
    debug_enter( "length of vectors: %d, %d",
                 x.len, y.len);
    matlib_index incx = 1;
    matlib_index incy = 1;

    err_check( (x.len != y.len), clean_up, 
               "Dimension mis-match (x.len: %d, y.len: %d)!",
               x.len, y.len);

    err_check( ((x.elem_p == NULL) || (y.elem_p == NULL)), clean_up, 
               "%s", "Null pointers encountered!");

    cblas_dcopy(x.len, x.elem_p, incx, y.elem_p, incy);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;

}

matlib_err matlib_zcopy
(
    const matlib_zv x,
          matlib_zv y
)
/* A X P Y
 * y <- a * x + y
 *
 * */
{
    debug_enter( "length of vectors: %d, %d",
                 x.len, y.len);
    matlib_index incx = 1;
    matlib_index incy = 1;

    err_check( (x.len != y.len), clean_up, 
               "Dimension mis-match (x.len: %d, y.len: %d)!",
               x.len, y.len);

    err_check( ((x.elem_p == NULL) || (y.elem_p == NULL)), clean_up, 
               "%s", "Null pointers encountered!");
    
    cblas_zcopy(x.len, x.elem_p, incx, y.elem_p, incy);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}


matlib_real matlib_xnrm2(matlib_xv x)
{
    debug_enter( "length of real vector: %d", x.len);
 
    matlib_index incx = 1;
    err_check( (x.elem_p == NULL), clean_up, 
               "%s", "Null pointer encountered!");
    matlib_real r = cblas_dnrm2(x.len, x.elem_p, incx);

    debug_exit("Exit Status: %s", "SUCCESS");
    return r;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_NAN;
}

matlib_real matlib_znrm2(matlib_zv x)
{
    debug_enter( "length of complex vector: %d", x.len);
    matlib_index incx = 1;
    err_check( (x.elem_p == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    matlib_real r = cblas_dznrm2(x.len, x.elem_p, incx);

    debug_exit("Exit Status: %s", "SUCCESS");
    return r;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_NAN;

}
/*============================================================================*/

matlib_err matlib_xaxpy
(
    const matlib_real alpha,
    const matlib_xv   x,
          matlib_xv   y
)
/* A X P Y
 * y <- a * x + y
 *
 * */
{
    debug_enter( "alpha: % 0.16f, length of vectors: %d, %d",
                 alpha, x.len, y.len);
    matlib_index incx = 1;
    matlib_index incy = 1;

    err_check( (x.len != y.len), clean_up, 
               "Dimension mis-match (x.len: %d, y.len: %d)!",
               x.len, y.len);
    err_check( ((x.elem_p == NULL) || (y.elem_p == NULL)), clean_up, 
               "%s", "Null pointers encountered!");

    cblas_daxpy(x.len, alpha, x.elem_p, incx, y.elem_p, incy);


    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

matlib_err matlib_zaxpy
(
    const matlib_complex alpha,
    const matlib_zv      x,
          matlib_zv      y
)
/* A X P Y
 * y <- a * x + y
 *
 * */
{
    debug_enter( "alpha: % 0.16f%+0.16fi, length of vectors: %d, %d",
                 creal(alpha), cimag(alpha), x.len, y.len);
    matlib_index incx = 1;
    matlib_index incy = 1;

    err_check( (x.len != y.len), clean_up, 
               "Dimension mis-match (x.len: %d, y.len: %d)!",
               x.len, y.len);
    err_check( ((x.elem_p == NULL) || (y.elem_p == NULL)), clean_up, 
               "%s", "Null pointers encountered!");

    cblas_zaxpy(x.len, &alpha, x.elem_p, incx, y.elem_p, incy);


    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

matlib_err matlib_zaxpby
(
    const matlib_complex alpha,
    const matlib_zv      x,
    const matlib_complex beta,
          matlib_zv      y
)
/* A X P B Y
 * y <- a * x + b * y
 *
 * */
{
    debug_enter( "alpha: % 0.16f%+0.16fi, beta: % 0.16f%+0.16fi, "
                 "length of vectors: %d, %d",
                 alpha, beta, x.len, y.len);
    matlib_index incx = 1;
    matlib_index incy = 1;

    err_check( (x.len != y.len), clean_up, 
               "Dimension mis-match (x.len: %d, y.len: %d)!",
               x.len, y.len);

    err_check( ((x.elem_p == NULL) || (y.elem_p == NULL)), clean_up, 
               "%s", "Null pointers encountered!");
    
    /* apply the scaling beta if it is not equal to 1.0 */ 
    if((creal(beta)!=1) || (cimag(beta)!=0))
    {
        cblas_zscal(y.len, &beta, y.elem_p, incy);
    }
    cblas_zaxpy(x.len, &alpha, x.elem_p, incx, y.elem_p, incy);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}
/*============================================================================*/

matlib_real matlib_xdot
(
    const matlib_xv x,
    const matlib_xv y
)
/* DOT
 * y <- a * x + y
 *
 * */
{
    debug_enter( "length of vectors: %d, %d", x.len, y.len);
    matlib_index incx = 1;
    matlib_index incy = 1;

    err_check( (x.len != y.len), clean_up, 
               "Dimension mis-match (x.len: %d, y.len: %d)!",
               x.len, y.len);

    debug_exit("exit status: %d", MATLIB_SUCCESS);
    matlib_real r = cblas_ddot(x.len, x.elem_p, incx, y.elem_p, incy);

    debug_exit("Exit Status: %s", "SUCCESS");
    return r;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_NAN;
}

/*============================================================================+/
 |BLAS Level II Routines
/+============================================================================*/
matlib_err matlib_xgemv
(
    const matlib_real alpha,
    const matlib_xm   A, 
    const matlib_xv   u,
    const matlib_real beta,
          matlib_xv   v
)
/* 
 * GE M V
 * v <-- alpha A * u + beta * v
 *
 * */

{
    debug_enter("%s", "");

    /* check if the input has NULL poindexers */
    err_check(    (A.elem_p == NULL) 
               || (u.elem_p == NULL) 
               || (v.elem_p == NULL), clean_up, 
               "%s", "Null pointer(s) encountered!");
    
    matlib_index incu = 1;
    matlib_index incv = 1;
    matlib_index strideA;

    bool dim_OK   = false;
    bool order_OK = false;

    /* Check if the order is okay! */ 
    _CBLAS_ORDER order_enum;
    
    if (A.order == MATLIB_COL_MAJOR)
    {
        strideA    = A.lenc;
        order_enum = CblasColMajor;
        order_OK   = true;
    }
    else if (A.order == MATLIB_ROW_MAJOR)
    {
        strideA    = A.lenr;
        order_enum = CblasRowMajor;
        order_OK   = true;
    }
    err_check( !order_OK, clean_up, 
               "Order of the matrix unknown (order: %s)!",
               MATLIB_ORDER_ENUM2STR(A.order));
    
    /* Check if the dimension is okay! */ 
    _CBLAS_OP op_enum;

    if ((A.op == MATLIB_TRANS) || (A.op == MATLIB_CONJ_TRANS))
    {
        dim_OK  = ((A.lenc == u.len) && (A.lenr == v.len));
        op_enum = CblasTrans;
    }
    else if (A.op == MATLIB_NO_TRANS)
    {
        dim_OK  = ((A.lenr == u.len) && (A.lenc == v.len));
        op_enum = CblasNoTrans;
    }
    err_check( !dim_OK, clean_up, 
               "Dimension mis-match in computing v <- alpha op(A) * u + beta * u"
               "(A: %d-by-%d, op(A): %s, "
               "u: %d-by-1, v: %d-by-1)!",
               A.lenc, A.lenr, MATLIB_OP_ENUM2STR(A.op), 
               u.len, v.len );

    cblas_dgemv( order_enum, 
                 op_enum, 
                 A.lenc, 
                 A.lenr, 
                 alpha, 
                 A.elem_p, 
                 strideA,
                 u.elem_p, 
                 incu, 
                 beta, 
                 v.elem_p, 
                 incv);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}


matlib_err matlib_zgemv
(
    const matlib_complex alpha,
    const matlib_zm      A, 
    const matlib_zv      u,
    const matlib_complex beta,
          matlib_zv      v
)
/* 
 * GE M V
 * v <-- alpha A * u + beta * v
 *
 * */

{
    debug_enter("%s", "");

    /* check if the input has NULL poindexers */
    err_check(    (A.elem_p == NULL) 
               || (u.elem_p == NULL) 
               || (v.elem_p == NULL), clean_up, 
               "%s", "Null pointer(s) encountered!");
    
    matlib_index incu = 1;
    matlib_index incv = 1;
    matlib_index strideA;

    bool dim_OK   = false;
    bool order_OK = false;
    
    _CBLAS_ORDER order_enum;
    
    if (A.order == MATLIB_COL_MAJOR)
    {
        strideA    = A.lenc;
        order_enum = CblasColMajor;
        order_OK   = true;
    }
    else if (A.order == MATLIB_ROW_MAJOR)
    {
        strideA    = A.lenr;
        order_enum = CblasRowMajor;
        order_OK   = true;
    }
    err_check( !order_OK, clean_up, 
               "Order of the matrix unknown (order: %s)!",
               MATLIB_ORDER_ENUM2STR(A.order));

    _CBLAS_OP op_enum;

    if (A.op == MATLIB_NO_TRANS)
    {
        dim_OK  = ((A.lenr == u.len) && (A.lenc == v.len));
        op_enum = CblasNoTrans;
    }
    else
    { 
        dim_OK = ((A.lenc == u.len) && (A.lenr == v.len));
        if (A.op == MATLIB_TRANS)
        {
            op_enum = CblasTrans;
        }
        else if (A.op == MATLIB_CONJ_TRANS)
        {
            op_enum = CblasConjTrans;
        }
    }
    err_check( !dim_OK, clean_up, 
               "Dimension mis-match in computing v <- alpha op(A) * u + beta * u"
               "(A: %d-by-%d, op(A): %s, u: %d-by-1, v: %d-by-1)!",
               A.lenc, A.lenr, MATLIB_OP_ENUM2STR(A.op), 
               u.len, v.len );

    debug_body( "alpha: % 0.16f%+0.16fi "
                "beta: % 0.16f%+0.16fi",
                creal(alpha), cimag(alpha),
                creal(beta) , cimag(beta));

    cblas_zgemv( order_enum, 
                 op_enum, 
                 A.lenc, 
                 A.lenr, 
                 &alpha, 
                 A.elem_p, 
                 strideA,
                 u.elem_p, 
                 incu, 
                 &beta, 
                 v.elem_p, 
                 incv);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}


void matlib_xcsrsymv
/* Double CSR Symmetric Matrix-Vector */ 
(
    const MATLIB_UPLO      uplo_enum,
          matlib_xm_sparse A, 
    const matlib_xv        u,
          matlib_xv        v
)
/* 
 * A: sparse matrix
 * v <-- alpha A * u + beta * v
 *
 * */

{
    debug_enter("%s", "");
    /* check if the input has NULL poindexers */
    assert(((A.elem_p != NULL) && (u.elem_p != NULL)) && (v.elem_p != NULL));
    
    /* check if the dimensions of the matrices are correct */ 
    
    char uplo_char[1];
    
    if (uplo_enum == MATLIB_UPPER)
    {
        uplo_char[0] = 'U';
    }
    else if (uplo_enum == MATLIB_LOWER)
    {
        uplo_char[0] = 'L';
    }
    else
    {
        term_execb("unknown upper/lower part (uplo_enum: %d)", uplo_enum);
    }
    debug_body( "size of A: %d-by-%d, nnz: %d", 
                A.lenc, A.lenr,
                A.rowIn[A.lenc] );

    mkl_cspblas_dcsrsymv( uplo_char, 
                          &A.lenc, 
                          A.elem_p, 
                          A.rowIn,
                          A.colIn,
                          u.elem_p, 
                          v.elem_p);

    debug_exit("%s", "");
}
void matlib_zcsrsymv
/* Double CSR Symmetric Matrix-Vector */ 
(
    const MATLIB_UPLO      uplo_enum,
          matlib_zm_sparse A, 
    const matlib_zv        u,
          matlib_zv        v
)
/* 
 * Z GE M V
 * v <-- alpha A * u + beta * v
 *
 * */

{
    debug_enter("%s", "");
    /* check if the input has NULL poindexers */
    assert(((A.elem_p != NULL) && (u.elem_p != NULL)) && (v.elem_p != NULL));
    
    /* check if the dimensions of the matrices are correct */ 
    
    char uplo_char[1];
    
    if (uplo_enum == MATLIB_UPPER)
    {
        uplo_char[0] = 'U';
    }
    else if (uplo_enum == MATLIB_LOWER)
    {
        uplo_char[0] = 'L';
    }
    else
    {
        term_execb("unknown upper/lower part (uplo_enum: %d)", uplo_enum);
    }
    debug_body( "size of A: %d-by-%d, nnz: %d", 
                A.lenc, A.lenr,
                A.rowIn[A.lenc] );

    mkl_cspblas_zcsrsymv( uplo_char, 
                          &A.lenc, 
                          A.elem_p, 
                          A.rowIn,
                          A.colIn,
                          u.elem_p, 
                          v.elem_p);

    debug_exit("%s", "");
}


/*============================================================================+/
 |BLAS Level III Routines
/+============================================================================*/
matlib_err matlib_xgemm
(
    const matlib_real alpha,
    const matlib_xm   A, 
    const matlib_xm   B, 
    const matlib_real beta,
          matlib_xm   C
)
/* 
 * GE M M
 * C <-- alpha A * B + beta * C
 *
 * */
{
    debug_enter("%s", "");

    /* Check if the input has NULL pointers! */
    err_check( (((A.elem_p == NULL) || (A.elem_p == NULL)) || (C.elem_p == NULL)),
               clean_up, "%s", "Null pointers encountered!");
    
    bool dim_OK   = false;
    bool order_OK = false;

    /* Check if the order is okay! */ 
    bool is_COL_MAJOR, is_ROW_MAJOR;

    is_COL_MAJOR = ((A.order == MATLIB_COL_MAJOR) && (B.order == MATLIB_COL_MAJOR)) && 
                    (C.order == MATLIB_COL_MAJOR);

    is_ROW_MAJOR = ((A.order == MATLIB_ROW_MAJOR) && (B.order == MATLIB_ROW_MAJOR)) && 
                    (C.order == MATLIB_ROW_MAJOR);

    CBLAS_ORDER order_enum;
    matlib_index strideA, strideB, strideC;

    if (is_COL_MAJOR)
    {
        order_enum = CblasColMajor;
        strideA = A.lenc;
        strideB = B.lenc;
        strideC = C.lenc;
        
        order_OK = true;
    }
    else if(is_ROW_MAJOR)
    {
        order_enum = CblasRowMajor;
        strideA = A.lenr;
        strideB = B.lenr;
        strideC = C.lenr;

        order_OK = true;
    }

    err_check( !order_OK, clean_up, 
               "Order of the matrices incorrect (A: %s, B: %s, C: %s)!",
               MATLIB_ORDER_ENUM2STR(A.order),
               MATLIB_ORDER_ENUM2STR(B.order),
               MATLIB_ORDER_ENUM2STR(C.order));

    /* Check if the dimensions of the matrices are correct! */ 
    matlib_index lenc_opA_and_C, lenr_opB_and_C, lenr_opA_and_lenc_opB;
    _CBLAS_OP op_enum_A;
    _CBLAS_OP op_enum_B;

    if (A.op == MATLIB_NO_TRANS)
    {
        op_enum_A = CblasNoTrans;
        if(B.op == MATLIB_NO_TRANS)
        {
            op_enum_B = CblasNoTrans;
            dim_OK    = ((A.lenr == B.lenc) && (A.lenc == C.lenc)) && (B.lenr == C.lenr);
        }
        else if(B.op == MATLIB_TRANS)
        {
            op_enum_B = CblasTrans;
            dim_OK    = ((A.lenr == B.lenr) && (A.lenc == C.lenc)) && (B.lenc == C.lenr);
        }
        /* If the dimensions are okay :
         * */ 
        lenr_opA_and_lenc_opB = A.lenr;
    }
    else if((A.op == MATLIB_TRANS) || (A.op == MATLIB_CONJ_TRANS))
    {
        op_enum_A = CblasTrans;
        if(B.op == MATLIB_NO_TRANS)
        {
            op_enum_B = CblasNoTrans;
            dim_OK    = ((A.lenc == B.lenc) && (A.lenr == C.lenc)) && (B.lenr == C.lenr);
        }
        else if(B.op == MATLIB_TRANS)
        {
            op_enum_B = CblasTrans;
            dim_OK    = ((A.lenc == B.lenr) && (A.lenr == C.lenc)) && (B.lenc == C.lenr);
        }
        /* If the dimensions are okay :
         * */ 
        lenr_opA_and_lenc_opB = A.lenc;
    }
    err_check( !dim_OK, clean_up, 
               "Dimension mismatch in computing C<-alpha op(A) * op(B) + beta * C "
               "(A:%d-by-%d, op(A): %s, B:%d-by-%d, op(B): %s, C:%d-by-%d)!",
               A.lenc, A.lenr, MATLIB_OP_ENUM2STR(A.op), 
               B.lenc, B.lenr, MATLIB_OP_ENUM2STR(B.op),
               C.lenc, C.lenr);
    
    lenc_opA_and_C = C.lenc;
    lenr_opB_and_C = C.lenr;

    cblas_dgemm( order_enum,
                 op_enum_A,  /* for A   */ 
                 op_enum_B,  /* for B   */
                 lenc_opA_and_C,
                 lenr_opB_and_C,
                 lenr_opA_and_lenc_opB,
                 alpha,
                 A.elem_p, 
                 strideA,
                 B.elem_p,
                 strideB,
                 beta,
                 C.elem_p, 
                 strideC
               );

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}


matlib_err matlib_zgemm
(
    const matlib_complex alpha,
    const matlib_zm      A, 
    const matlib_zm      B, 
    const matlib_complex beta,
          matlib_zm      C
)
/* 
 * Z GE M M
 * C <-- alpha A * B + beta * C
 *
 * */
{
    debug_enter("%s", "");

    /* Check if the input has NULL pointers! */
    err_check( (((A.elem_p == NULL) || (A.elem_p == NULL)) || (C.elem_p == NULL)),
               clean_up, "%s", "Null pointers encountered!");

    bool dim_OK   = false;
    bool order_OK = false;

    /* Check if the order is okay! */ 
    bool is_COL_MAJOR, is_ROW_MAJOR;

    is_COL_MAJOR = ((A.order == MATLIB_COL_MAJOR) && (B.order == MATLIB_COL_MAJOR)) && 
                    (C.order == MATLIB_COL_MAJOR);

    is_ROW_MAJOR = ((A.order == MATLIB_ROW_MAJOR) && (B.order == MATLIB_ROW_MAJOR)) && 
                    (C.order == MATLIB_ROW_MAJOR);

    CBLAS_ORDER order_enum;
    matlib_index strideA, strideB, strideC;

    if (is_COL_MAJOR)
    {
        order_enum = CblasColMajor;
        strideA = A.lenc;
        strideB = B.lenc;
        strideC = C.lenc;
        
        order_OK = true;
    }
    else if(is_ROW_MAJOR)
    {
        order_enum = CblasRowMajor;
        strideA = A.lenr;
        strideB = B.lenr;
        strideC = C.lenr;

        order_OK = true;
    }

    err_check( !order_OK, clean_up, 
               "Order of the matrices incorrect (A: %s, B: %s, C: %s)!",
               MATLIB_ORDER_ENUM2STR(A.order),
               MATLIB_ORDER_ENUM2STR(B.order),
               MATLIB_ORDER_ENUM2STR(C.order));


    /* Check if the dimensions of the matrices are correct! */ 
    matlib_index lenc_opA_and_C, lenr_opB_and_C, lenr_opA_and_lenc_opB;

    _CBLAS_OP op_enum_A;
    _CBLAS_OP op_enum_B;

    if (A.op == MATLIB_NO_TRANS)
    {
        op_enum_A = CblasNoTrans;
        if(B.op == MATLIB_NO_TRANS)
        {
            op_enum_B = CblasNoTrans;
            dim_OK = ((A.lenr == B.lenc) && (A.lenc == C.lenc)) && (B.lenr == C.lenr);
        }
        else
        {
            dim_OK = ((A.lenr == B.lenr) && (A.lenc == C.lenc)) && (B.lenc == C.lenr);
            if(B.op == MATLIB_TRANS)
                op_enum_B = CblasTrans;
            else
                op_enum_B = CblasConjTrans;
        }
        /* If the dimensions are okay :
         * */ 
        lenr_opA_and_lenc_opB = A.lenr;
    }
    else 
    {
        if(A.op == MATLIB_TRANS)
            op_enum_A = CblasTrans;
        else
            op_enum_A = CblasConjTrans;
        if(B.op == MATLIB_NO_TRANS)
        {
            op_enum_B = CblasNoTrans;
            dim_OK = ((A.lenc == B.lenc) && (A.lenr == C.lenc)) && (B.lenr == C.lenr);
        }
        else
        {
            dim_OK = ((A.lenc == B.lenr) && (A.lenr == C.lenc)) && (B.lenc == C.lenr);
            if(B.op == MATLIB_TRANS)
                op_enum_B = CblasTrans;
            else
                op_enum_B = CblasConjTrans;
        }
        /* If the dimensions are okay :
         * */ 
        lenr_opA_and_lenc_opB = A.lenc;
    }
    err_check( !dim_OK, clean_up, 
               "Dimension mismatch in computing C<-alpha op(A) * op(B) + beta * C "
               "(A:%d-by-%d, op(A): %d, B:%d-by-%d, op(B): %d, C:%d-by-%d)!",
               A.lenc, A.lenr, MATLIB_OP_ENUM2STR(A.op), 
               B.lenc, B.lenr, MATLIB_OP_ENUM2STR(B.op),
               C.lenc, C.lenr);

    lenc_opA_and_C = C.lenc;
    lenr_opB_and_C = C.lenr;

    cblas_zgemm( order_enum,
                 op_enum_A,  /*  for A   */ 
                 op_enum_B,  /*  for B   */
                 lenc_opA_and_C,
                 lenr_opB_and_C,
                 lenr_opA_and_lenc_opB,
                 &alpha,
                 A.elem_p, 
                 strideA,
                 B.elem_p,
                 strideB,
                 &beta,
                 C.elem_p, 
                 strideC
               );


    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

#if 0

matlib_err func
(
    matlib_
)
/* 
 * 
 *
 *
 * */

{
    debug_enter("%s", "");
    matlib_err error = MATLIB_SUCCESS;

    err_check( <expr>, clean_up, "%s", "");

    error = function_call();
    err_check( (error == MATLIB_FAILURE), clean_up, "%s", "");



    debug_exit("Exit Status: %s", MATLIB_ERR_ENUM2STR(error));
    return error;

clean_up:

    debug_exit("Exit Status: %s", "FAILURE");
    return error;

}/* End of func */ 



#endif

/*============================================================================*/


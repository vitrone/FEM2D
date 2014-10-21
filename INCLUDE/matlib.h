#ifndef MATLIB_H
#define MATLIB_H

/*============================================================================+/
 | Include all the dependencies
/+============================================================================*/
#include "complex.h"
#include "basic.h"
#include "debug.h"
#include "ehandler.h"

#define MKL_Complex16 matlib_complex 
#define _CBLAS_ORDER  CBLAS_ORDER
#define _CBLAS_OP     CBLAS_TRANSPOSE


#ifdef MATLIB_NTRACE_DATA
    #define _TRACE_DATA 0
#else
    #define _TRACE_DATA 1
#endif


/* Writing data tracing blocks */
#define BEGIN_DTRACE if(_TRACE_DATA){
#define END_DTRACE }

/*============================================================================*/

#define MATLIB_64

#ifdef MATLIB_64
    /* Indexing types */ 
    typedef unsigned long long int matlib_index;
    /* Real numbers */ 
    typedef double matlib_real;
    /* Complex numbers */ 
    typedef double complex matlib_complex;
    /* Integers */
    typedef long long int matlib_int;
#else
    /* Indexing types */ 
    typedef unsigned int matlib_index;
    /* Real numbers */ 
    typedef double matlib_real;
    /* Complex numbers */ 
    typedef double complex matlib_complex;
    /* Integers */
    typedef int matlib_int;

#endif

#define MATLIB_NAN NAN

/* Error handling */ 
typedef enum
{
    MATLIB_SUCCESS = 0,
    MATLIB_FAILURE = -1

} matlib_err;

#define MATLIB_ERR_ENUM2STR(error)               \
    (MATLIB_SUCCESS == error ? "SUCCESS":   \
     (MATLIB_FAILURE == error ? "FAILURE": "UNKOWN"))


#define matlib_dv matlib_xv
#define matlib_dm matlib_xm

/* Vector type is almost never used. Vectors are always interpreted to be of
 * appropriate form.
 * */ 
typedef enum
{
    MATLIB_COL_VECT,
    MATLIB_ROW_VECT

} MATLIB_VECT_T;

typedef enum
{
    MATLIB_NO_TRANS,
    MATLIB_TRANS,
    MATLIB_CONJ_TRANS

} MATLIB_OP;

#define MATLIB_OP_ENUM2STR(op)               \
    (MATLIB_NO_TRANS == op ? "NO TRANSPOSE":   \
     (MATLIB_TRANS == op ? "TRANSPOSE":        \
      (MATLIB_CONJ_TRANS == op? "CONJUGATE TRANSPOSE": "UNKOWN")))


/* ORDER: 
 * COL_MAJOR: elements within a column are contiguous in memory
 * ROW_MAJOR: elements  with in a row are contiguous in memory
 * 
 * Matrix elements are stored in a array as 
 * A[i][j]--> COL_MAJOR: *(pA+i+stride*j), ROW_MAJOR: *(pA+i*stride+j)
 * "stride" here represents the constant spacing two row/column elements
 * For COL_MAJOR, stride>= length of columns (lenc)
 * For ROW_MAJOR, stride>= length of rows    (lenr)
 * 
 * stride corresponds to leading dimension argument of arrays in FOTRAN
 * All matrix data types declared in this file assume stride=lenc or lenr unless
 * otherwise stated.
 *
 * */ 
typedef enum
{
    MATLIB_COL_MAJOR, 
    MATLIB_ROW_MAJOR 

} MATLIB_ORDER;

#define MATLIB_ORDER_ENUM2STR(order)                     \
    (MATLIB_COL_MAJOR == order ? "COLUMN MAJOR":         \
     (MATLIB_ROW_MAJOR == order ? "ROW MAJOR":"UNKNOWN"))

/* Data types for matrix library 
 *
 * matlib_xv : Real Vector
 * matlib_xm : Real Matrix
 *
 * */ 
typedef struct
{
    matlib_index  len;
    MATLIB_VECT_T type;
    matlib_real*  elem_p;

} matlib_xv;

typedef struct
{
    matlib_index    len;
    MATLIB_VECT_T   type;
    matlib_complex* elem_p;

} matlib_zv;

/* The operation field is to tell if the meaningful entries belong to the current
 * matrix or the transposed version of it. An actual taranspose operation is
 * carried out only if strictly needed.
 * */ 
typedef struct
{
    matlib_index lenc; /* length of columns */ 
    matlib_index lenr; /* length of rows    */ 
    MATLIB_ORDER order;
    MATLIB_OP    op; /* operation */ 
    matlib_real* elem_p;

} matlib_xm;

typedef struct
{
    matlib_index lenc; /* length of columns */ 
    matlib_index lenr; /* length of rows    */ 
    MATLIB_ORDER order;
    MATLIB_OP op; /* operation */ 
    matlib_complex* elem_p;

} matlib_zm;

/*============================================================================*/
/* SPARSE FORMATS
 * 
 * CSR3: 
 * i-th row, matrix elements: elem_p[rowIn[i]]... elem_p[rowIn[i+1]-1]
 * column matlib_index : colIn[rowIn[i]] < colIn[rowIn[i]+1] <...< col[rowIn[i+1]-1]
 * 
 * nr of non-zero elements: rowIn[M.lenc]
 *
 * Wherever "CSR" is used it refers to "CSR3".
 * */

typedef enum
{
    MATLIB_CSR3,
    MATLIB_CSC3, 
    MATLIB_COO,
    MATLIB_DIA

} MATLIB_SPARSE;

typedef enum
{
    MATLIB_UPPER,
    MATLIB_LOWER

} MATLIB_UPLO;

typedef struct
{
    matlib_index  lenc; /* length of columns */ 
    matlib_index  lenr; /* length of rows    */
    matlib_index* rowIn;
    matlib_index* colIn;
    MATLIB_SPARSE format;
    matlib_real*  elem_p;

} matlib_xm_sparse;

typedef struct
{
    matlib_index    lenc; /* length of columns */ 
    matlib_index    lenr; /* length of rows    */
    matlib_index*   rowIn;
    matlib_index*   colIn;
    MATLIB_SPARSE   format;
    matlib_complex* elem_p;

} matlib_zm_sparse;

/* This data structure is meant to store several sparse matrices which have the
 * same sparsity structure. This proves useful for solving initial value
 * problems where the potential is a time dependent function. The mass matrix
 * for several time steps can be computed and stored at once.*/ 
typedef struct
{
    matlib_index  lenc; /* length of columns */ 
    matlib_index  lenr; /* length of rows    */
    matlib_index* rowIn;
    matlib_index* colIn;
    MATLIB_SPARSE format;
    matlib_index  nsparse;
    matlib_real** elem_p;

} matlib_xm_nsparse;

typedef struct
{
    matlib_index  lenc; /* length of columns */ 
    matlib_index  lenr; /* length of rows    */
    matlib_index* rowIn;
    matlib_index* colIn;
    MATLIB_SPARSE format;
    matlib_index  nsparse;
    matlib_complex** elem_p;

} matlib_zm_nsparse;


/*============================================================================*/
/* Define MACROS */ 

#define matlib_free(ptr)                                                 \
    do{ if(ptr!=NULL)                                                    \
        free((void*) ptr);                                               \
      } while (0)

#define DEBUG_PRINT_XV(NAME, fmt,...)                                    \
    do {                                                                 \
         if(_TRACE_DATA){                                                \
            matlib_index NAME ## _i;                                     \
            for((NAME ## _i)=0;                                          \
                (NAME ## _i)<(NAME.len);                                 \
                (NAME ##_i)++){                                          \
                fprintf( stderr,                                         \
                         "%s:%d:%s: -(" fmt #NAME "[%d] = % 0.16f)\n",   \
                         __FILE__,                                       \
                         __LINE__,                                       \
                         __func__,                                       \
                        __VA_ARGS__,                                     \
                        NAME ## _i,                                      \
                        *((NAME.elem_p)+(NAME ## _i)));                  \
            }                                                            \
         }                                                               \
    } while (0)           

/* For complex Vectors */ 
#define DEBUG_PRINT_ZV(NAME, fmt,...)                                    \
    do {                                                                 \
         if(_TRACE_DATA){                                                \
            matlib_index NAME ## _i;                                     \
            for((NAME ## _i)=0;                                          \
                (NAME ## _i)<(NAME.len);                                 \
                (NAME ##_i)++){                                          \
                fprintf( stderr,                                         \
                    "%s:%d:%s: -(" fmt #NAME "[%d] = % 0.16f%+0.16fi)\n",\
                         __FILE__,                                       \
                         __LINE__,                                       \
                         __func__,                                       \
                        __VA_ARGS__,                                     \
                        NAME ## _i,                                      \
                        creal(*((NAME.elem_p)+(NAME ## _i))),            \
                        cimag(*((NAME.elem_p)+(NAME ## _i))));           \
            }                                                            \
         }                                                               \
    } while (0)                                                          

/* Debuging a matlib_real matrix */ 
#define DEBUG_PRINT_XM(NAME, fmt,...)                                    \
    do {                                                                 \
         if(_TRACE_DATA){                                                \
            matlib_index NAME ## _i, NAME ## _j;                         \
            matlib_index NAME ## _col_st, NAME ## _row_st;               \
            if ((NAME.order) == MATLIB_COL_MAJOR){                       \
                (NAME ## _col_st) = 1;                                   \
                (NAME ## _row_st) = (NAME.lenc);                         \
            }                                                            \
            else{                                                        \
                (NAME ## _col_st) = (NAME.lenr);                         \
                (NAME ## _row_st) = 1;                                   \
            }                                                            \
            for((NAME ## _i)=0;                                          \
                (NAME ## _i)<(NAME.lenc);                                \
                (NAME ##_i)++){                                          \
                for((NAME ## _j)=0;                                      \
                    (NAME ## _j)<(NAME.lenr);                            \
                    (NAME ## _j)++){                                     \
                    fprintf( stderr,                                     \
                        "%s:%d:%s: -(" fmt #NAME "[%d][%d] = % 0.16f)\n",\
                             __FILE__,                                   \
                             __LINE__,                                   \
                             __func__,                                   \
                             __VA_ARGS__,                                \
                            NAME ## _i, NAME ## _j,                      \
                            *((NAME.elem_p)+                             \
                            (NAME ## _i)*(NAME ## _col_st)+              \
                            (NAME ## _j)*(NAME ## _row_st)));            \
                }                                                        \
            }                                                            \
         }                                                               \
    } while (0)                                                          


/* Debuging a complex matrix */ 
#define DEBUG_PRINT_ZM(NAME, fmt,...)                                    \
    do {                                                                 \
         if(_TRACE_DATA){                                                \
            matlib_index NAME ## _i, NAME ## _j;                         \
            matlib_index NAME ## _col_st, NAME ## _row_st;               \
            if ((NAME.order) == MATLIB_COL_MAJOR){                       \
                (NAME ## _col_st) = 1;                                   \
                (NAME ## _row_st) = (NAME.lenc);                         \
            }                                                            \
            else{                                                        \
                (NAME ## _col_st) = (NAME.lenr);                         \
                (NAME ## _row_st) = 1;                                   \
            }                                                            \
            for((NAME ## _i)=0;                                          \
                (NAME ## _i)<(NAME.lenc);                                \
                (NAME ##_i)++){                                          \
                for((NAME ## _j)=0;                                      \
                    (NAME ## _j)<(NAME.lenr);                            \
                    (NAME ## _j)++){                                     \
                    fprintf( stderr,                                     \
                             "%s:%d:%s: -(" fmt                          \
                             #NAME "[%d][%d] = % 0.16f %+0.16fi)\n",     \
                             __FILE__,                                   \
                             __LINE__,                                   \
                             __func__,                                   \
                             __VA_ARGS__,                                \
                            NAME ## _i, NAME ## _j,                      \
                            creal(*((NAME.elem_p)+                       \
                            (NAME ## _i)*(NAME ## _col_st)+              \
                            (NAME ## _j)*(NAME ## _row_st))),            \
                            cimag(*((NAME.elem_p)+                       \
                            (NAME ## _i)*(NAME ## _col_st)+              \
                            (NAME ## _j)*(NAME ## _row_st))));           \
                }                                                        \
            }                                                            \
         }                                                               \
    } while (0)                                                          
/*============================================================================*/

/* Convert a matrix into vector 
 *
 * V = MK_VM(B);
 *
 * */ 
#define MK_VM(B) { .len  = B.lenc * B.lenr,                             \
                   .type = MATLIB_COL_VECT,                             \
                   .elem_p = B.elem_p }
 
/*============================================================================+/
 |Allocation of memory
/+============================================================================*/
matlib_err matlib_create_zv
(
    matlib_index  length,
    matlib_zv*    v,
    MATLIB_VECT_T type_enum
);

matlib_err matlib_create_xv
(
    matlib_index  length,
    matlib_xv*    v,
    MATLIB_VECT_T type_enum
);

matlib_err matlib_create_zm
( 
    matlib_index lenc,
    matlib_index lenr,
    matlib_zm*   M,
    MATLIB_ORDER order_enum,
    MATLIB_OP    op_enum
);
matlib_err matlib_create_xm
( 
    matlib_index lenc,
    matlib_index lenr,
    matlib_xm*   M,
    MATLIB_ORDER order_enum,
    MATLIB_OP    op_enum
);
/*============================================================================+/
 |BLAS Level I Routines
/+============================================================================*/

matlib_err matlib_xcopy
(
    const matlib_xv x,
          matlib_xv y
);

matlib_err matlib_zcopy
(
    const matlib_zv x,
          matlib_zv y
);


matlib_real matlib_xnrm2(matlib_xv x);
matlib_real matlib_znrm2(matlib_zv x);

matlib_err matlib_xaxpy
(
    const matlib_real alpha,
    const matlib_xv   x,
          matlib_xv   y
);
matlib_err matlib_zaxpy
(
    const matlib_complex alpha,
    const matlib_zv      x,
          matlib_zv      y
);
matlib_err matlib_zaxpby
(
    const matlib_complex alpha,
    const matlib_zv      x,
    const matlib_complex beta,
          matlib_zv      y
);
matlib_real matlib_xdot
(
    const matlib_xv x,
    const matlib_xv y
);

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
);
matlib_err matlib_zgemv
(
    const matlib_complex alpha,
    const matlib_zm      A, 
    const matlib_zv      u,
    const matlib_complex beta,
          matlib_zv      v
);

void matlib_xcsrsymv
/* Double CSR Symmetric Matrix-Vector */ 
(
    const MATLIB_UPLO      uplo_enum,
          matlib_xm_sparse A, 
    const matlib_xv        u,
          matlib_xv        v
);

void matlib_zcsrsymv
/* Double CSR Symmetric Matrix-Vector */ 
(
    const MATLIB_UPLO      uplo_enum,
          matlib_zm_sparse A, 
    const matlib_zv        u,
          matlib_zv        v
);
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
);
matlib_err matlib_zgemm
(
    const matlib_complex alpha,
    const matlib_zm      A, 
    const matlib_zm      B, 
    const matlib_complex beta,
          matlib_zm      C
);

/*============================================================================+/
 | IO functions for vectors and matrices
/+============================================================================*/

void matlib_xmwrite_csv(char* file_name, matlib_xm M);
void matlib_zmwrite_csv(char* file_name, matlib_zm M);

void matlib_xvwrite_csv
(
    char*        file_name, 
    matlib_index n,
    matlib_xv    v[n]
);

void matlib_zvwrite_csv
(
    char*        file_name, 
    matlib_index n,
    matlib_zv    v[n]
);

void matlib_xzvwrite_csv
(
    char*        file_name, 
    matlib_index m,
    matlib_xv    u[m],
    matlib_index n,
    matlib_zv    v[n]
);
/*============================================================================+/
 | Linear Solver Interface
/+============================================================================*/

/*=====================[Linear solver based on PARDISO]=======================*/
#ifdef MATLIB_64
    #define _MATLIB_PARDISO PARDISO_64 
#else
    #define _MATLIB_PARDISO PARDISO 
#endif

#define PARDISO_NIPARAM (64)

#define _E_PARDISO_INIT                 (-2)
#define _E_PARDISO_FREE                 (-1) /* release all memory */ 
#define _E_PARDISO_FREE_LU              (0)
#define _E_PARDISO_ANALYSIS_AND_FACTOR  (12)
#define _E_PARDISO_SOLVE_AND_REFINE     (33)

#define _E_PARDISO_REAL_SYM_INDEF  (-2)
#define _E_PARDISO_REAL_SYM_PDEF   ( 2)
#define _E_PARDISO_COMPLEX_SYM     ( 6)

typedef enum
{
    PARDISO_RHS,
    PARDISO_LHS

} PARDISO_SOLVEC;

typedef enum
{
    PARDISO_INIT                = -2,
    PARDISO_FREE                = -1, /* release all memory */ 
    PARDISO_FREE_LU             = 0,
    PARDISO_ANALYSIS_AND_FACTOR = 12,
    PARDISO_SOLVE_AND_REFINE    = 33

} PARDISO_PHASE;

#define PARDISO_PHASE_ENUM2STR(phase)                                          \
    (PARDISO_INIT == phase ? "INITIALIZE":                                     \
     (PARDISO_FREE == phase ? "FREE MEMORY":                                   \
      (PARDISO_FREE_LU == phase ? "FREE MEMORY LU":                            \
       (PARDISO_ANALYSIS_AND_FACTOR == phase ? "ANALYSIS AND FACTORIZATION":   \
        (PARDISO_SOLVE_AND_REFINE == phase ? "SOLVE AND REFINE": "UNKNOWN")))))

typedef enum
{
    PARDISO_REAL_SYM_INDEF = -2,
    PARDISO_REAL_SYM_PDEF  =  2,
    PARDISO_COMPLEX_SYM    =  6

} PARDISO_MTYPE;

#define PARDISO_MTYPE_ENUM2STR(mtype)                                \
    (PARDISO_REAL_SYM_INDEF == ? "REAL SYMMETRIC INDEFINITE":        \
     (PARDISO_REAL_SYM_PDEF  == ? "REAL SYMMTERIC POSITIVE DEFINITE":\
      (PARDISO_COMPLEX_SYM    == ? "COMPLEX SYMMTERIC":"UNKNOWN")))

typedef struct
{
    void*          ptr[PARDISO_NIPARAM];
    matlib_int     iparam[PARDISO_NIPARAM];
    PARDISO_PHASE  phase_enum;
    PARDISO_MTYPE  mtype;
    PARDISO_SOLVEC sol_enum;
    matlib_index   nsparse; /* number of sparse matrices with same 
                               sparsity structure */ 
    matlib_index   mnum;   /* matrix number to be used for solution */
    void*          smat_p; /* Sparse matrix struct */ 
    void*          rhs_p;  /* vector struct        */ 
    void*          sol_p;  /* vector struct        */ 

} pardiso_solver_t;

void matlib_pardiso(pardiso_solver_t* data);

/*============================================================================*/


#endif


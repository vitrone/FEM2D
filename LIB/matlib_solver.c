/*============================================================================+/
 | Name: matlib_solver.c
/+============================================================================*/
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define NDEBUG
#define MATLIB_NTRACE_DATA

#include "mkl.h"
#include "mkl_pardiso.h"
#include "matlib.h"

/*============================================================================*/
static inline matlib_err matlib_xmtype
(
    PARDISO_MTYPE mtype_enum,
    matlib_int *mtype
);
static inline matlib_err matlib_zmtype
(
    PARDISO_MTYPE mtype_enum,
    matlib_int *mtype
);

static inline matlib_err matlib_pardiso_init(pardiso_solver_t* data);

static inline matlib_err matlib_xpardiso_AF( pardiso_solver_t* data);
static inline matlib_err matlib_zpardiso_AF( pardiso_solver_t* data);
static inline matlib_err matlib_xpardiso_SR( pardiso_solver_t* data);
static inline matlib_err matlib_zpardiso_SR( pardiso_solver_t* data);

static inline matlib_err matlib_xpardiso_free(pardiso_solver_t* data);
static inline matlib_err matlib_zpardiso_free(pardiso_solver_t* data);
/*============================================================================*/
static inline matlib_err matlib_xmtype
(
    PARDISO_MTYPE mtype_enum,
    matlib_int *mtype
)
{
    matlib_err error = MATLIB_SUCCESS;
    switch(mtype_enum)
    {
        case PARDISO_REAL_SYM_INDEF:
            *mtype = _E_PARDISO_REAL_SYM_INDEF;
            break;
        case PARDISO_REAL_SYM_PDEF:
            *mtype = _E_PARDISO_REAL_SYM_PDEF;
            break;
        case PARDISO_REAL_UNSYM:
            *mtype = _E_PARDISO_REAL_UNSYM;
            break;
        default:
            *mtype = 0;
            error = MATLIB_FAILURE;
    }
    return error;
}

static inline matlib_err matlib_zmtype
(
    PARDISO_MTYPE mtype_enum,
    matlib_int *mtype
)
{
    matlib_err error = MATLIB_SUCCESS;
    switch(mtype_enum)
    {
        case PARDISO_COMPLEX_SYM:
            *mtype = _E_PARDISO_COMPLEX_SYM;
            break;
        case PARDISO_COMPLEX_HERM_INDEF:
            *mtype = _E_PARDISO_COMPLEX_HERM_INDEF;
            break;
        case PARDISO_COMPLEX_HERM_PDEF:
            *mtype = _E_PARDISO_COMPLEX_HERM_PDEF;
            break;
        case PARDISO_COMPLEX_UNSYM:
            *mtype = _E_PARDISO_COMPLEX_UNSYM;
            break;
        default:
            error = MATLIB_FAILURE;
            *mtype = 0;
    }

    return error;
}

static inline matlib_err matlib_pardiso_init(pardiso_solver_t* data)
{
    err_check(data == NULL, clean_up, "%s", "Null pointer encountered!");
    err_check(    (data->iparam == NULL) 
               || (data->ptr    == NULL), clean_up, 
               "%s", "Null pointer(s) encountered!");

    matlib_index i;
    /* SETUP PARDISO CONTROL PARAMETERS */
    for (i = 0; i < PARDISO_NIPARAM; i++)
    {
        data->iparam[i] = 0;
        data->ptr[i]    = 0;
    }
    data->iparam[0] = 1; /* Don't use default values */ 
    data->iparam[1] = 2; /* Fill-in reducing odering for input matrix */ 
    data->iparam[3] = 0; /* Preconditioning */ 
    data->iparam[4] = 0; /*  */ 
    

    if (data->sol_enum == PARDISO_RHS)
    {
        data->iparam[5] = _E_PARDISO_RHS; /* Write solution into b */ 
    }
    else if (data->sol_enum == PARDISO_LHS)
    {
        data->iparam[5] = _E_PARDISO_LHS; /* Write solution into x */
    }
    else
    {
        data->sol_enum = PARDISO_SOLVEC_UNKNOWN;
    }

    err_check( data->sol_enum == PARDISO_SOLVEC_UNKNOWN, clean_up, 
               "Incorrect option for solution vector (write solution into: %d)", 
               PARDISO_SOLVEC_ENUM2STR(data->sol_enum));

    data->iparam[7]  = 2; /* Maximum number of iterative refinement steps, output 
                             reported in iparam[6] */ 
    data->iparam[9]  = 13; /* Perturbing pivot elements */ 
    data->iparam[10] = 0;  /* Disable scaling  */ 
    data->iparam[12] = 0;  /*  */ 
    data->iparam[17] = -1; /* Disable reporting of nnz  */
    data->iparam[18] =  1; /*  */ 
    data->iparam[34] =  1; /* Zero-based indexing  */ 
    
    debug_exit( "Exit Status: %s", "SuCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

static inline matlib_err matlib_xpardiso_AF(pardiso_solver_t* data)
{
    debug_body("%s", "Start testing PARDISO");
    matlib_int mtype;
    matlib_index nrhs = 1; /* Number of right hand sides  */ 

    /* Print statistical information in file */
    matlib_index msglvl = PARDISO_MSGLVL;

    matlib_int error = 0; /* Initialize error flag */
    debug_body("nr. sparse matrices: %d", data->nsparse);
    matlib_index maxfct = 1;
    matlib_index mnum   = 1;
    matlib_int phase_enum = _E_PARDISO_ANALYSIS_AND_FACTOR;

    matlib_err err = matlib_xmtype(data->mtype, &mtype);
    err_check( err == MATLIB_FAILURE, clean_up, 
               "Incorrect option for real matrix type (mtype: %d)", 
               PARDISO_MTYPE_ENUM2STR(data->mtype));
    if (data->nsparse > 1)
    {
        
        matlib_xm_nsparse* smat_p = (matlib_xm_nsparse*) data->smat_p;
        _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                         &mtype, &phase_enum,
                         &smat_p->lenc,
                         smat_p->elem_p[data->mnum-1],
                         smat_p->rowIn,
                         smat_p->colIn,
                         NULL, &nrhs,
                         data->iparam,
                         &msglvl, NULL, NULL,
                         &error);
    }
    else
    {
        matlib_xm_sparse* smat_p = (matlib_xm_sparse*) data->smat_p;
        _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                         &mtype, &phase_enum,
                         &smat_p->lenc,
                         smat_p->elem_p,
                         smat_p->rowIn,
                         smat_p->colIn,
                         NULL, &nrhs,
                         data->iparam,
                         &msglvl, NULL, NULL,
                         &error);
    }
    err_check( error != 0, clean_up, 
               "Analysis and factorization failed (error code: %d)", error);

    debug_exit( "Exit Status: %s", "SuCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}


static inline matlib_err matlib_zpardiso_AF(pardiso_solver_t* data)
{
    debug_body("%s", "Start testing PARDISO");
    matlib_int mtype;
    matlib_index nrhs = 1; /* Number of right hand sides  */ 

    /* Print statistical information in file */
    matlib_index msglvl = PARDISO_MSGLVL;

    matlib_int error = 0; /* Initialize error flag */
    debug_body("nr. sparse matrices: %d", data->nsparse);
    matlib_index maxfct = 1;
    matlib_index mnum   = 1;
    matlib_int phase_enum = _E_PARDISO_ANALYSIS_AND_FACTOR;

    matlib_err err = matlib_zmtype(data->mtype, &mtype);
    err_check( err == MATLIB_FAILURE, clean_up, 
               "Incorrect option for complex matrix type (mtype: %d)", 
               PARDISO_MTYPE_ENUM2STR(data->mtype));

    if (data->nsparse > 1)
    {
        
        matlib_zm_nsparse* smat_p = (matlib_zm_nsparse*) data->smat_p;
        _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                         &mtype, &phase_enum,
                         &smat_p->lenc,
                         smat_p->elem_p[data->mnum-1],
                         smat_p->rowIn,
                         smat_p->colIn,
                         NULL, &nrhs,
                         data->iparam,
                         &msglvl, NULL, NULL,
                         &error);
    }
    else
    {
        matlib_zm_sparse* smat_p = (matlib_zm_sparse*) data->smat_p;
        _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                         &mtype, &phase_enum,
                         &smat_p->lenc,
                         smat_p->elem_p,
                         smat_p->rowIn,
                         smat_p->colIn,
                         NULL, &nrhs,
                         data->iparam,
                         &msglvl, NULL, NULL,
                         &error);
    }
    err_check( error != 0, clean_up, 
               "Analysis and factorization failed (error code: %d)", error);

    debug_exit( "Exit Status: %s", "SuCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}



/*============================================================================*/

static inline matlib_err matlib_xpardiso_SR(pardiso_solver_t* data)
{
    matlib_int mtype;
    matlib_index nrhs = 1; /* Number of right hand sides  */ 

    /* Print statistical information in file */
    matlib_index msglvl = PARDISO_MSGLVL; 
    matlib_int   error  = 0; /* Initialize error flag */
    matlib_index maxfct = 1;
    matlib_index mnum = 1;

    matlib_int phase_enum = _E_PARDISO_SOLVE_AND_REFINE;
    debug_body("nr. sparse matrices: %d", data->nsparse);

    matlib_err err = matlib_xmtype(data->mtype, &mtype);
    err_check( err == MATLIB_FAILURE, clean_up, 
               "Incorrect option for real matrix type (mtype: %d)", 
               PARDISO_MTYPE_ENUM2STR(data->mtype));
    if (data->nsparse > 1)
    {
        matlib_xm_nsparse* smat_p = (matlib_xm_nsparse*) data->smat_p;
        matlib_xv* rhs_p = (matlib_xv*) data->rhs_p;
        matlib_xv* sol_p = (matlib_xv*) data->sol_p;


        _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                         &mtype, &phase_enum,
                         &smat_p->lenc,
                         smat_p->elem_p[data->mnum-1],
                         smat_p->rowIn,
                         smat_p->colIn,
                         NULL, &nrhs,
                         data->iparam,
                         &msglvl, 
                         rhs_p->elem_p, 
                         sol_p->elem_p,
                         &error);

    }
    else
    {
        matlib_xm_sparse* smat_p  = (matlib_xm_sparse*) data->smat_p;
        matlib_xv* rhs_p = (matlib_xv*) data->rhs_p;
        matlib_xv* sol_p = (matlib_xv*) data->sol_p;

        _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                         &mtype, &phase_enum,
                         &smat_p->lenc,
                         smat_p->elem_p,
                         smat_p->rowIn,
                         smat_p->colIn,
                         NULL, &nrhs,
                         data->iparam,
                         &msglvl, 
                         rhs_p->elem_p, 
                         sol_p->elem_p,
                         &error);
    }
    err_check( error != 0, clean_up, 
               "ERROR during solution and refinement (error code: %d)", error);

    debug_exit( "Exit Status: %s", "SuCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}


static inline matlib_err matlib_zpardiso_SR(pardiso_solver_t* data)
{
    matlib_int mtype;
    matlib_index nrhs = 1; /* Number of right hand sides  */ 

    /* Print statistical information in file */
    matlib_index msglvl = PARDISO_MSGLVL; 
    matlib_int   error  = 0; /* Initialize error flag */
    matlib_index maxfct = 1;
    matlib_index mnum   = 1;

    matlib_int phase_enum = _E_PARDISO_SOLVE_AND_REFINE;
    debug_body("nr. sparse matrices: %d", data->nsparse);

    matlib_err err = matlib_zmtype(data->mtype, &mtype);
    err_check( err == MATLIB_FAILURE, clean_up, 
               "Incorrect option for complex matrix type (mtype: %d)", 
               PARDISO_MTYPE_ENUM2STR(data->mtype));
    if (data->nsparse > 1)
    {
        matlib_zm_nsparse* smat_p = (matlib_zm_nsparse*) data->smat_p;
        matlib_zv* rhs_p = (matlib_zv*) data->rhs_p;
        matlib_zv* sol_p = (matlib_zv*) data->sol_p;

        _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                         &mtype, &phase_enum,
                         &smat_p->lenc,
                         smat_p->elem_p[data->mnum-1],
                         smat_p->rowIn,
                         smat_p->colIn,
                         NULL, &nrhs,
                         data->iparam,
                         &msglvl, 
                         rhs_p->elem_p, 
                         sol_p->elem_p,
                         &error);

    }
    else
    {
        matlib_zm_sparse* smat_p  = (matlib_zm_sparse*) data->smat_p;
        matlib_zv* rhs_p = (matlib_zv*) data->rhs_p;
        matlib_zv* sol_p = (matlib_zv*) data->sol_p;

        _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                         &mtype, &phase_enum,
                         &smat_p->lenc,
                         smat_p->elem_p,
                         smat_p->rowIn,
                         smat_p->colIn,
                         NULL, &nrhs,
                         data->iparam,
                         &msglvl, 
                         rhs_p->elem_p, 
                         sol_p->elem_p,
                         &error);
    }
    err_check( error != 0, clean_up, 
               "ERROR during solution and refinement (error code: %d)", error);

    debug_exit( "Exit Status: %s", "SuCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

static inline matlib_err matlib_xpardiso_free(pardiso_solver_t* data)
{
    matlib_int mtype;
    matlib_index msglvl = PARDISO_MSGLVL; /* Print statistical information in file */
    matlib_int   error  = 0; /* Initialize error flag */
    matlib_index maxfct = 1;
    matlib_index mnum   = 1;
        
    matlib_int phase_enum = _E_PARDISO_FREE;
    debug_body("nr. sparse matrices: %d", data->nsparse);

    matlib_err err = matlib_xmtype(data->mtype, &mtype);
    err_check( err == MATLIB_FAILURE, clean_up, 
               "Incorrect option for real matrix type (mtype: %d)", 
               PARDISO_MTYPE_ENUM2STR(data->mtype));

    if (data->nsparse > 1)
    {
        matlib_xm_nsparse* smat_p = (matlib_xm_nsparse*) data->smat_p;

        _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                         &mtype, &phase_enum,
                         &smat_p->lenc,
                         NULL,
                         smat_p->rowIn,
                         smat_p->colIn,
                         NULL, NULL,
                         data->iparam,
                         &msglvl, 
                         NULL, NULL,
                         &error);

    }
    else
    {
        matlib_xm_sparse* smat_p  = (matlib_xm_sparse*) data->smat_p;

        _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                         &mtype, &phase_enum,
                         &smat_p->lenc,
                         NULL,
                         smat_p->rowIn,
                         smat_p->colIn,
                         NULL, NULL,
                         data->iparam,
                         &msglvl, 
                         NULL, NULL,
                         &error);
    }
    debug_exit( "Exit Status: %s", "SuCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

static inline matlib_err matlib_zpardiso_free(pardiso_solver_t* data)
{
    matlib_int mtype;
    matlib_index msglvl = PARDISO_MSGLVL; /* Print statistical information in file */
    matlib_int   error  = 0; /* Initialize error flag */
    matlib_index maxfct = 1;
    matlib_index mnum   = 1;
        
    matlib_int phase_enum = _E_PARDISO_FREE;
    debug_body("nr. sparse matrices: %d", data->nsparse);

    matlib_err err = matlib_zmtype(data->mtype, &mtype);
    err_check( err == MATLIB_FAILURE, clean_up, 
               "Incorrect option for complex matrix type (mtype: %d)", 
               PARDISO_MTYPE_ENUM2STR(data->mtype));

    if (data->nsparse > 1)
    {
        matlib_zm_nsparse* smat_p = (matlib_zm_nsparse*) data->smat_p;

        _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                         &mtype, &phase_enum,
                         &smat_p->lenc,
                         NULL,
                         smat_p->rowIn,
                         smat_p->colIn,
                         NULL, NULL,
                         data->iparam,
                         &msglvl, 
                         NULL, NULL,
                         &error);

    }
    else
    {
        matlib_zm_sparse* smat_p  = (matlib_zm_sparse*) data->smat_p;

        _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                         &mtype, &phase_enum,
                         &smat_p->lenc,
                         NULL,
                         smat_p->rowIn,
                         smat_p->colIn,
                         NULL, NULL,
                         data->iparam,
                         &msglvl, 
                         NULL, NULL,
                         &error);
    }

    debug_exit( "Exit Status: %s", "SuCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

/*============================================================================*/

matlib_err matlib_xpardiso(pardiso_solver_t* data)
/* 
 * Handles complex as well as real matrices.
 *
 * */ 
{
    debug_enter("%s", "");
    matlib_err error = MATLIB_FAILURE;

    if (data->phase_enum == PARDISO_INIT)
    {
        error = matlib_pardiso_init(data);
    }
    else if (data->phase_enum == PARDISO_ANALYSIS_AND_FACTOR)
    {
        error = matlib_xpardiso_AF(data);
    }
    else if(data->phase_enum == PARDISO_SOLVE_AND_REFINE)
    {
        error = matlib_xpardiso_SR(data);
    }
    else if(data->phase_enum == PARDISO_FREE)
    {
        error = matlib_xpardiso_free(data);
    }
    return error;
}

matlib_err matlib_zpardiso(pardiso_solver_t* data)
/* 
 * Handles complex as well as real matrices.
 *
 * */ 
{
    debug_enter("%s", "");
    matlib_err error = MATLIB_FAILURE;

    if (data->phase_enum == PARDISO_INIT)
    {
        error = matlib_pardiso_init(data);
    }
    else if (data->phase_enum == PARDISO_ANALYSIS_AND_FACTOR)
    {
        error = matlib_zpardiso_AF(data);
    }
    else if(data->phase_enum == PARDISO_SOLVE_AND_REFINE)
    {
        error = matlib_zpardiso_SR(data);
    }
    else if(data->phase_enum == PARDISO_FREE)
    {
        error = matlib_zpardiso_free(data);
    }
    return error;
}



/*============================================================================*/


void matlib_pardiso(pardiso_solver_t* data)
/* 
 * Handles complex as well as real matrices.
 *
 * */ 
{
    debug_enter("%s", "");

    matlib_int mtype, phase_enum;
    if (data->phase_enum == PARDISO_INIT)
    {
        matlib_index i;
        /* SETUP PARDISO CONTROL PARAMETERS */
        for (i = 0; i < PARDISO_NIPARAM; i++)
        {
            data->iparam[i] = 0;
            data->ptr[i]    = 0;
        }
        data->iparam[0]  = 1; /* Don't use default values */ 
        data->iparam[1]  = 2; /* Fill-in reducing odering for input matrix */ 
        data->iparam[3]  = 0; /* Preconditioning */ 
        data->iparam[4]  = 0; /*  */ 

        if(data->sol_enum==PARDISO_RHS)
        {
            data->iparam[5]  = 1; /* Write solution into b */ 
        }
        else if(data->sol_enum==PARDISO_LHS)
        {
            data->iparam[5]  = 0; /* Write solution into x */
        }
        else
        {
          term_exec( "Incorrect option for solution vector (sol_enum:%d)", 
                     data->sol_enum);
        }

        data->iparam[7]  = 2; /* Maximum number of iterative refinement steps, output 
                                 reported in iparam[6] */ 
        data->iparam[9]  = 13; /* Perturbing pivot elements */ 
        data->iparam[10] = 0;  /* Disable scaling  */ 
        data->iparam[12] = 0;  /*  */ 
        data->iparam[17] = -1; /* Disable reporting of nnz  */
        data->iparam[18] =  1; /*  */ 
        data->iparam[34] =  1; /* Zero-based indexing  */ 
        debug_body("%s", "Initialized PARDISO control parameters");
    }
    else if (data->phase_enum == PARDISO_ANALYSIS_AND_FACTOR)
    {
        debug_body("%s", "Start testing PARDISO");
        matlib_index nrhs  = 1; /* Number of right hand sides  */ 

        /* Print statistical information in file */
        matlib_index msglvl = PARDISO_MSGLVL;

        matlib_int error    = 0; /* Initialize error flag */
        debug_body("nr. sparse matrices: %d", data->nsparse);
        matlib_index maxfct = 1;
        matlib_index mnum = 1;
        phase_enum = _E_PARDISO_ANALYSIS_AND_FACTOR;

        if(data->nsparse>1)
        {
            if(data->mtype == PARDISO_COMPLEX_SYM)
            {
                mtype = _E_PARDISO_COMPLEX_SYM;
                debug_body("%s", "Factoring a complex symmetric matrix");
                matlib_zm_nsparse* smat_p = (matlib_zm_nsparse*) data->smat_p;
                _MATLIB_PARDISO( data->ptr, &maxfct, &mnum,
                          &mtype, &phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p[data->mnum-1],
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, NULL, NULL,
                          &error);
            }
            else
            {
                switch(data->mtype)
                {
                    case PARDISO_REAL_SYM_INDEF:
                        mtype = _E_PARDISO_REAL_SYM_INDEF;
                        break;
                    case PARDISO_REAL_SYM_PDEF:
                        mtype = _E_PARDISO_REAL_SYM_PDEF;
                        break;
                    default:
                        mtype = _E_PARDISO_REAL_SYM_INDEF;
                }
                
                debug_body("%s", "Factoring a real symmetric matrix");
                matlib_xm_nsparse* smat_p = (matlib_xm_nsparse*) data->smat_p;
                _MATLIB_PARDISO ( data->ptr, &maxfct, &mnum,
                          &mtype, &phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p[data->mnum-1],
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, NULL, NULL,
                          &error);
            }
        }
        else
        {
            if(data->mtype == PARDISO_COMPLEX_SYM)
            {
                mtype = _E_PARDISO_COMPLEX_SYM;
                debug_body("%s", "Factoring a complex symmetric matrix");
                matlib_zm_sparse* smat_p = (matlib_zm_sparse*) data->smat_p;
                _MATLIB_PARDISO ( data->ptr, &maxfct, &mnum,
                          &mtype, &phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p,
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, NULL, NULL,
                          &error);
            }
            else
            {
                switch(data->mtype)
                {
                    case PARDISO_REAL_SYM_INDEF:
                        mtype = _E_PARDISO_REAL_SYM_INDEF;
                        break;
                    case PARDISO_REAL_SYM_PDEF:
                        mtype = _E_PARDISO_REAL_SYM_PDEF;
                        break;
                    default:
                        mtype = _E_PARDISO_REAL_SYM_INDEF;
                }
                
                debug_body("%s", "Factoring a real symmetric matrix");
                matlib_xm_sparse* smat_p = (matlib_xm_sparse*) data->smat_p;
                _MATLIB_PARDISO ( data->ptr, &maxfct, &mnum,
                          &mtype, &phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p,
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, NULL, NULL,
                          &error);
            }
        }
        if (error != 0)
        {
          term_exec("Analysis and factorization failed (error code: %d)", error);
        }
        debug_body("%s", "Analysis and factorization completed");
    }
    else if(data->phase_enum == PARDISO_SOLVE_AND_REFINE)
    {
        matlib_index nrhs  = 1; /* Number of right hand sides  */ 

        /* Print statistical information in file */
        matlib_index msglvl = PARDISO_MSGLVL; 
        matlib_int   error  = 0; /* Initialize error flag */
        matlib_index maxfct = 1;
        matlib_index mnum = 1;

        phase_enum = _E_PARDISO_SOLVE_AND_REFINE;
        debug_body("nr. sparse matrices: %d", data->nsparse);
        if(data->nsparse>1)
        {
            if(data->mtype == PARDISO_COMPLEX_SYM)
            {

                mtype = _E_PARDISO_COMPLEX_SYM;
                matlib_zm_nsparse* smat_p = (matlib_zm_nsparse*) data->smat_p;
                matlib_zv* rhs_p  = (matlib_zv*) data->rhs_p;
                matlib_zv* sol_p  = (matlib_zv*) data->sol_p;

                BEGIN_DTRACE
                    debug_print("dimension of the sparse square matrix: %d", smat_p->lenc);
                    matlib_complex* ptr = smat_p->elem_p[data->mnum-1];
                    for(matlib_index i=0; i<smat_p->lenc; i++)
                    {
                        for(matlib_index j=smat_p->rowIn[i]; j<smat_p->rowIn[i+1]; j++)
                        {
                            debug_print( "M(%d,%d): % 0.16f %+0.16fi", 
                                         i, smat_p->colIn[j], ptr[j]);
                        }
                    }
                END_DTRACE

                _MATLIB_PARDISO ( data->ptr, &maxfct, &(mnum),
                          &mtype, &phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p[data->mnum-1],
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, 
                          rhs_p->elem_p, 
                          sol_p->elem_p,
                          &error);
            }
            else
            {
                switch(data->mtype)
                {
                    case PARDISO_REAL_SYM_INDEF:
                        mtype = _E_PARDISO_REAL_SYM_INDEF;
                        break;
                    case PARDISO_REAL_SYM_PDEF:
                        mtype = _E_PARDISO_REAL_SYM_PDEF;
                        break;
                    default:
                        mtype = _E_PARDISO_REAL_SYM_INDEF;
                }
                
                matlib_xm_nsparse* smat_p = (matlib_xm_nsparse*) data->smat_p;
                matlib_xv* rhs_p  = (matlib_xv*) data->rhs_p;
                matlib_xv* sol_p  = (matlib_xv*) data->sol_p;


                _MATLIB_PARDISO ( data->ptr, &maxfct, &mnum,
                          &mtype, &phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p[data->mnum-1],
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, 
                          rhs_p->elem_p, 
                          sol_p->elem_p,
                          &error);

            }
        }
        else
        {
            if(data->mtype == PARDISO_COMPLEX_SYM)
            {
                matlib_zm_sparse* smat_p = (matlib_zm_sparse*) data->smat_p;
                matlib_zv* rhs_p  = (matlib_zv*) data->rhs_p;
                matlib_zv* sol_p  = (matlib_zv*) data->sol_p;
                mtype = _E_PARDISO_COMPLEX_SYM;

                debug_body("%s", "Solving a complex symmetric system.");
                _MATLIB_PARDISO ( data->ptr, &maxfct, &mnum,
                          &mtype, &phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p,
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, 
                          rhs_p->elem_p, 
                          sol_p->elem_p,
                          &error);
            }
            else
            {
                switch(data->mtype)
                {
                    case PARDISO_REAL_SYM_INDEF:
                        mtype = _E_PARDISO_REAL_SYM_INDEF;
                        break;
                    case PARDISO_REAL_SYM_PDEF:
                        mtype = _E_PARDISO_REAL_SYM_PDEF;
                        break;
                    default:
                        mtype = _E_PARDISO_REAL_SYM_INDEF;
                }
                
                matlib_xm_sparse* smat_p  = (matlib_xm_sparse*) data->smat_p;
                matlib_xv* rhs_p  = (matlib_xv*) data->rhs_p;
                matlib_xv* sol_p  = (matlib_xv*) data->sol_p;

                BEGIN_DTRACE
                    debug_print("dimension of the sparse square matrix: %d", smat_p->lenc);
                    matlib_real* ptr = smat_p->elem_p;
                    for(matlib_index i=0; i<smat_p->lenc; i++)
                    {
                        for(matlib_index j=smat_p->rowIn[i]; j<smat_p->rowIn[i+1]; j++)
                        {
                            debug_print( "M(%d,%d): % 0.16f", 
                                         i, smat_p->colIn[j], ptr[j]);
                        }
                    }
                END_DTRACE
                _MATLIB_PARDISO ( data->ptr, &maxfct, &mnum,
                          &mtype, &phase_enum,
                          &smat_p->lenc,
                          smat_p->elem_p,
                          smat_p->rowIn,
                          smat_p->colIn,
                          NULL, &nrhs,
                          data->iparam,
                          &msglvl, 
                          rhs_p->elem_p, 
                          sol_p->elem_p,
                          &error);
                BEGIN_DTRACE
                    for (matlib_index j=0; j<rhs_p->len; j++)
                    {
                        debug_print("rhs[%d]: %0.16f", j, rhs_p->elem_p[j]);
                    }
                    for (matlib_index j=0; j<sol_p->len; j++)
                    {
                        debug_print("sol[%d]: %0.16f", j, sol_p->elem_p[j]);
                    }
                END_DTRACE
            }
        }

        if (error != 0)
        {
            term_exec("ERROR during solution and refinement: %d", error);
        }
    
        debug_body("%s", "Solution and refinement completed");
    }
    else if(data->phase_enum == PARDISO_FREE)
    {

        matlib_index msglvl = PARDISO_MSGLVL; /* Print statistical information in file */
        matlib_int   error  = 0; /* Initialize error flag */
        matlib_index maxfct = 1;
        matlib_index mnum = 1;
        
        phase_enum = _E_PARDISO_FREE;
        if(data->mtype == PARDISO_COMPLEX_SYM)
        {
            mtype = _E_PARDISO_COMPLEX_SYM;
            matlib_zm_nsparse* smat_p = (matlib_zm_nsparse*) data->smat_p;

            _MATLIB_PARDISO ( data->ptr, &maxfct, &mnum,
                      &mtype, &phase_enum,
                      &smat_p->lenc,
                      NULL,
                      smat_p->rowIn,
                      smat_p->colIn,
                      NULL, NULL,
                      data->iparam,
                      &msglvl, 
                      NULL, 
                      NULL,
                      &error);
        }
        else
        {
            switch(data->mtype)
            {
                case PARDISO_REAL_SYM_INDEF:
                    mtype = _E_PARDISO_REAL_SYM_INDEF;
                    break;
                case PARDISO_REAL_SYM_PDEF:
                    mtype = _E_PARDISO_REAL_SYM_PDEF;
                    break;
                default:
                    mtype = _E_PARDISO_REAL_SYM_INDEF;
            }
                
            matlib_xm_nsparse* smat_p = (matlib_xm_nsparse*) data->smat_p;

            _MATLIB_PARDISO ( data->ptr, &maxfct, &mnum,
                      &mtype, &phase_enum,
                      &smat_p->lenc,
                      NULL,
                      smat_p->rowIn,
                      smat_p->colIn,
                      NULL, NULL,
                      data->iparam,
                      &msglvl, 
                      NULL, NULL,
                      &error);
        }
    }
    
    debug_exit("%s", "");
}


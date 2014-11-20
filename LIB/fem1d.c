/*============================================================================+/
 | 
 |  Name   : fem1d.c                                                      
 |  Author : Vishal Vaibhav                                               
 |                                                                        
 |  Description : This module defines all the functions for implementing  
 |  finite element method. See the documentation for mathematical details.
 |                                                                        
 |
/+============================================================================*/

#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include "mkl.h"

#define NDEBUG
#define MATLIB_NTRACE_DATA

#include "fem1d.h"
#include "assert.h"

/*============================================================================*/
/*
 *  Reference domain: [-1, 1]
 * */ 


const matlib_index FEM1D_DIM = 1;
/* Vertex indices */ 
const matlib_index FEM1D_INDEX_V1 = 0; 
const matlib_index FEM1D_INDEX_V2 = 1;   
const matlib_index FEM1D_NV = 2;   /* nr of vertices of a triangular element */ 

/* Ref. triangle vertices */
const matlib_real FEM1D_VERT[2] = {-1.0, 1.0};

const matlib_real FEM1D_MEMI[2][2] = { {2.0/3.0,  1.0/3.0}, 
                                       {1.0/3.0,  2.0/3.0}};

const matlib_real FEM1D_MESI[2][2] = { { 1.0/2.0,  -1.0/2.0}, 
                                       {-1.0/2.0,   1.0/2.0}};
/* 
 * Combination of basis functions for mass-matrix
 * (vphi1, vphi1) --> FEM1D_INDEX_V11
 * (vphi1, vphi2) --> FEM1D_INDEX_V12
 * (vphi2, vphi2) --> FEM1D_INDEX_V22
 *
 * */ 
const matlib_index FEM1D_INDEX_V11 = 0;
const matlib_index FEM1D_INDEX_V12 = 1;
const matlib_index FEM1D_INDEX_V22 = 2;
const matlib_index FEM1D_NR_COMBI  = 3;

/*============================================================================*/

fem1d_err fem1d_refbasis
(
    const matlib_xv   xi,
          matlib_xm* vphi
)
{

    debug_enter( "length of xi: %d", xi.len);

    err_check( (xi.elem_p == NULL), clean_up, 
               "%s", "Null pointer ecountered!");

    fem1d_err error; 
    error = matlib_create_xm( xi.len, FEM1D_NV, vphi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Memory allocation for basis function-value matrix failed!");


    int i;
    matlib_real *ptr1, *ptr2;
    
    /* vphi_1 = (1 - xi)/2 
     * */ 
    ptr2 = vphi->elem_p + FEM1D_INDEX_V1 * xi.len;

    for ( ptr1 = xi.elem_p; 
          ptr1 < (xi.elem_p + xi.len);
          ptr1++, ptr2++)
    {
        *ptr2 = 0.5 * ( 1.0 - *ptr1);
    }
    
    /* vphi_2 = (1 + xi)/2 
     * */ 
    ptr2 = vphi->elem_p + FEM1D_INDEX_V2 * xi.len;

    for ( ptr1 = xi.elem_p; 
          ptr1 < (xi.elem_p + xi.len);
          ptr1++, ptr2++)
    {
        *ptr2 =  0.5 * ( 1.0 + *ptr1);
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM1D_FAILURE;
}
/*============================================================================*/
fem1d_err fem1d_create_na
(
    fem1d_na *na, 
    matlib_index nr_nodes,
    matlib_real* nodes
)
{
    na->len        = nr_nodes;
    na->nr_domains = nr_nodes - 1;
    na->elem_p     = nodes;

    errno = 0;
    na->jacob = calloc( na->nr_domains, sizeof(matlib_real));
    err_check( (na->jacob == NULL), clean_up,
               "%s: Failed allocate memory for Jacobian(length: %d)", 
                strerror(errno), na->nr_domains);

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM1D_FAILURE;

}

fem1d_err fem1d_calc_jacobian(fem1d_na *na)
{

    err_check(    (na->elem_p == NULL) 
               || (na->jacob  == NULL), clean_up, 
               "%s", "Null pointers ecountered!");
    matlib_index i;
    for ( i = 0; i < na->nr_domains; i++)
    {
        na->jacob[i] = 0.5 * (na->elem_p[i+1] - na->elem_p[i]);
        err_check( na->jacob[i] < MATLIB_TOL, clean_up, 
                   "Jacobian is zero of negative (jacobian: %0.16f)!",
                   na->jacob[i]);
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM1D_FAILURE;

}



/*============================================================================*/

fem1d_err fem1d_ref2mesh
(
    fem1d_na   na,
    matlib_xm  vphi,
    matlib_xv* x
)
{
    err_check(    (na.elem_p   == NULL) 
               || (vphi.elem_p == NULL) 
               || (x           == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check( vphi.lenr != FEM1D_NV, clean_up, 
               "Size of second input incorrect (input: %d-by-%d)!",
               vphi.lenc, vphi.lenr);

    DEBUG_PRINT_XM(vphi, "%s: ", "Ref basis");

    fem1d_err error;
    error = matlib_create_xv( vphi.lenc * na.nr_domains, x, MATLIB_COL_VECT);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "Initialization of coordinate "
               "array (of length: %d) failed!", 
               vphi.lenc * na.nr_domains);

    matlib_real alpha = 1.0;
    matlib_real beta  = 0.0;

    matlib_real* vptr; /* vertex pointer */ 

    matlib_xv nodes_vec;
    matlib_create_xv( FEM1D_NV, &nodes_vec, MATLIB_COL_VECT);
 
    matlib_xv x_tmp = { .len    = vphi.lenc,  
                        .elem_p = x->elem_p };
    
    /* 
     * vphi: size xi.len - by - FEM1D_NV 
     * Convert from COL_MAJOR to ROW_MAJOR 
     * */ 

    matlib_xm vphi_tmp = { .order = MATLIB_ROW_MAJOR,
                           .op    = MATLIB_TRANS,
                           .lenr  = vphi.lenc,
                           .lenc  = FEM1D_NV,
                           .elem_p = vphi.elem_p };

    matlib_real* ptr = na.elem_p;
    for ( ; ptr < (na.elem_p + na.nr_domains); ptr++)
    {
        nodes_vec.elem_p[FEM1D_INDEX_V1] = ptr[FEM1D_INDEX_V1];
        nodes_vec.elem_p[FEM1D_INDEX_V2] = ptr[FEM1D_INDEX_V2];
        
        error = matlib_xgemv( alpha, vphi_tmp, nodes_vec, beta, x_tmp );
        err_check( error == FEM1D_FAILURE, clean_up, 
                   "%s", "Matrix-vector multiplication failed!");

        x_tmp.elem_p = x_tmp.elem_p + vphi.lenc;
    }
    
    matlib_free(nodes_vec.elem_p);
    debug_exit( "Total nr. of sampling points: %d, Exit Status: %s",
                x->len, "SUCCESS" );
    return FEM1D_SUCCESS;
    
clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

} /* End of fem1d_ref2mesh */ 




/*============================================================================*/
fem1d_err fem1d_xinterp
(
    const fem1d_na  na,
    const matlib_xv u_nodes,
    const matlib_xm vphi,
          matlib_xv u_interp
)
/* 
 * Computes the interpolation over FEM-elements.
 *
 *
 * */

{
    debug_enter( "length of u_nodes: %d", 
                 u_nodes.len);
    
    err_check(    (na.elem_p      == NULL)   
               || (u_nodes.elem_p == NULL) 
               || (vphi.elem_p    == NULL), clean_up, 
               "%s", "Null pointer(s) ecountered!");
    
    /* Length of u_nodes = nr. of nodes */ 
    err_check( (u_nodes.len != na.len ), clean_up, 
               "Dimension mismatch (u_nodes: %d, nr nodes: )",
               u_nodes.len, na.len);

    /* Length of u_interp = nr of quad. nodes x nr of domains 
     * */ 
    err_check( (u_interp.len != (vphi.lenc * na.nr_domains) ), clean_up, 
               "Dimension mismatch (u_interp: %d)", u_interp.len);
    
    matlib_real alpha = 1.0;
    matlib_real beta  = 0.0;

    matlib_index INDEX_REAL = 0;

    matlib_xv u1_tmp = { .len    = FEM1D_NV, 
                         .type   = MATLIB_COL_VECT,
                         .elem_p = u_nodes.elem_p };

    matlib_xv u2_tmp = { .len    = vphi.lenc, 
                         .type   = MATLIB_COL_VECT,
                         .elem_p = u_interp.elem_p };



    /* Convert vphi to row-major format */ 
    matlib_xm vphi_tmp = { .order = MATLIB_ROW_MAJOR,
                           .op    = MATLIB_TRANS,
                           .lenr  = vphi.lenc,
                           .lenc  = FEM1D_NV,
                           .elem_p = vphi.elem_p };


    fem1d_err error;
    for (; u1_tmp.elem_p < u_nodes.elem_p + na.nr_domains; (u1_tmp.elem_p)++)
    {
        error = matlib_xgemv( alpha, vphi_tmp, u1_tmp, beta, u2_tmp );
        err_check( error == FEM1D_FAILURE, clean_up, 
                   "%s", "Matrix-vector multiplication failed!");
  
        u2_tmp.elem_p = u2_tmp.elem_p + vphi.lenc; 

    }
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;
    
clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

}/* End of fem1d_xinterp */ 

fem1d_err fem1d_zinterp
(
    const fem1d_na  na,
    const matlib_zv u_nodes,
    const matlib_xm vphi,
          matlib_zv u_interp
)
/* 
 * Computes the interpolation over FEM-elements.
 *
 *
 * */

{
    debug_enter( "length of u_nodes: %d", 
                 u_nodes.len);
    
    err_check(    (na.elem_p      == NULL)   
               || (u_nodes.elem_p == NULL) 
               || (vphi.elem_p    == NULL), clean_up, 
               "%s", "Null pointers ecountered!");
    
    /* Length of u_nodes = nr. of nodes */ 
    err_check( (u_nodes.len != na.len ), clean_up, 
               "Dimension mismatch (u_nodes: %d, nr nodes: )",
               u_nodes.len, na.len);

    /* Length of u_interp = nr of quad. nodes x nr of domains 
     * */ 
    err_check( (u_interp.len != (vphi.lenc * na.nr_domains) ), clean_up, 
               "Dimension mismatch (u_interp: %d)", u_interp.len);
    
    matlib_real alpha = 1.0;
    matlib_real beta  = 0.0;

    matlib_index COMPLEX_DIM = 2; /* Dimension of a compelex plane */ 

    matlib_xm u1_tmp = { .lenc   = FEM1D_NV, 
                         .lenr   = COMPLEX_DIM, 
                         .order  = MATLIB_ROW_MAJOR,
                         .op     = MATLIB_NO_TRANS,
                         .elem_p = (matlib_real*)u_nodes.elem_p };

    matlib_xm u2_tmp = { .lenc   = vphi.lenc, 
                         .lenr   = COMPLEX_DIM, 
                         .order  = MATLIB_ROW_MAJOR,
                         .op     = MATLIB_NO_TRANS,
                         .elem_p = (matlib_real*)u_interp.elem_p };


    matlib_xm vphi_tmp = { .order = MATLIB_ROW_MAJOR,
                           .op    = MATLIB_TRANS,
                           .lenr  = vphi.lenc,
                           .lenc  = FEM1D_NV,
                           .elem_p = vphi.elem_p };

    fem1d_err error;
    for (matlib_index i = 0; i < na.nr_domains; i++)
    {
        error = matlib_xgemm( alpha, vphi_tmp, u1_tmp, beta, u2_tmp );
        err_check( error == FEM1D_FAILURE, clean_up, 
                   "%s", "Matrix-vector multiplication failed!");
  
        u1_tmp.elem_p = u1_tmp.elem_p + COMPLEX_DIM; 
        u2_tmp.elem_p = u2_tmp.elem_p + COMPLEX_DIM * vphi.lenc; 

    }
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;
    
clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

}/* End of fem1d_zinterp */ 

/*============================================================================*/
matlib_real fem1d_xnormL2
(
    fem1d_na  na,
    matlib_xv u_qnodes,
    matlib_xv quadW
)
{
    debug_enter( "nr of domains: %d, "
                 "nr of nodes: %d, "
                 "nr of quadrature nodes: %d",
                 na.nr_domains, na.len, quadW.len);
    
    err_check(    (na.elem_p       == NULL) 
               || (na.jacob        == NULL) 
               || (u_qnodes.elem_p == NULL) 
               || (quadW.elem_p    == NULL), clean_up,
               "%s", "Null pointers ecountered!");

    err_check( (u_qnodes.len != quadW.len * na.nr_domains), clean_up, 
               "%s", "Dimension mis-match!");

    matlib_real* ptr = u_qnodes.elem_p;
    
    matlib_real norm = 0;
    matlib_real jacob;

    for ( matlib_index i = 0; i < na.nr_domains; i++)
    {
        jacob = na.jacob[i];
        for (matlib_index j = 0; j< quadW.len; j++, ptr++)
        {
            norm += (*ptr * *ptr * quadW.elem_p[j] * jacob);
        }
    }
    debug_exit("Exit Status: %s", "SUCCESS" );
    return sqrt(norm);

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

} /* End of fem1d_xnormL2 */ 



matlib_real fem1d_znormL2
(
    fem1d_na  na,
    matlib_zv u_qnodes,
    matlib_xv quadW
)
{
    debug_enter( "nr of domains: %d, "
                 "nr of nodes: %d, "
                 "nr of quadrature nodes: %d",
                 na.nr_domains, na.len, quadW.len);
    
    err_check(    (na.elem_p       == NULL) 
               || (na.jacob        == NULL) 
               || (u_qnodes.elem_p == NULL) 
               || (quadW.elem_p    == NULL), clean_up,
               "%s", "Null pointers ecountered!");

    err_check( (u_qnodes.len != quadW.len * na.nr_domains), clean_up, 
               "%s", "Dimension mis-match!");

    matlib_complex* ptr = u_qnodes.elem_p;
    
    matlib_real norm = 0;
    matlib_real jacob;

    for ( matlib_index i = 0; i< na.nr_domains; i++)
    {
        jacob = na.jacob[i];
        for (matlib_index j=0; j< quadW.len; j++, ptr++)
        {
            norm += (*ptr * conj(*ptr) * quadW.elem_p[j] * jacob);
        }
    }
    debug_exit("Exit Status: %s", "SUCCESS" );
    return sqrt(norm);

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

}/* End of fem1d_znormL2 */ 

/*============================================================================*/
matlib_real fem1d_xiprod 
/* Inner product iprod = (u, conj(v))_{\Omega} */ 
(
    fem1d_na  na,
    matlib_xv u_qnodes,
    matlib_xv v_qnodes,
    matlib_xv quadW
)
{
    debug_enter( "nr of domains: %d, "
                 "nr of nodes: %d, "
                 "nr of quadrature nodes: %d",
                 na.nr_domains, na.len, quadW.len);
    
    err_check(    (na.elem_p       == NULL)    
               || (na.jacob        == NULL) 
               || (u_qnodes.elem_p == NULL) 
               || (v_qnodes.elem_p == NULL) 
               || (quadW.elem_p    == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check( (u_qnodes.len != v_qnodes.len),
               clean_up, "Dimension mis-match for field vectors (u: %d, v: %d)!",
               v_qnodes.len, u_qnodes.len );

    err_check( (v_qnodes.len != quadW.len * na.nr_domains), clean_up, 
               "%s", "Dimension mis-match!");

    matlib_real* uptr = u_qnodes.elem_p;
    matlib_real* vptr = v_qnodes.elem_p;
    
    matlib_real iprod = 0;
    matlib_real jacob;

    for ( matlib_index i = 0; i< na.nr_domains; i++)
    {
        jacob = na.jacob[i];
        for (matlib_index j = 0; j < quadW.len; j++, uptr++, vptr++)
        {
            iprod += (*uptr * *vptr * quadW.elem_p[j] * jacob);
        }
    }
    debug_exit("Exit Status: %s", "SUCCESS" );
    return iprod;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

} /* End of fem1d_xiprod */ 

matlib_real fem1d_ziprod 
/* Inner product iprod = (u, conj(v))_{\Omega} */ 
(
    fem1d_na  na,
    matlib_zv u_qnodes,
    matlib_zv v_qnodes,
    matlib_xv quadW
)
{
    debug_enter( "nr of domains: %d, "
                 "nr of nodes: %d, "
                 "nr of quadrature nodes: %d",
                 na.nr_domains, na.len, quadW.len);
    
    err_check(    (na.elem_p       == NULL)    
               || (na.jacob        == NULL) 
               || (u_qnodes.elem_p == NULL) 
               || (v_qnodes.elem_p == NULL) 
               || (quadW.elem_p    == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check( (u_qnodes.len != v_qnodes.len),
               clean_up, "Dimension mis-match for field vectors (u: %d, v: %d)!",
               v_qnodes.len, u_qnodes.len );

    err_check( (v_qnodes.len != quadW.len * na.nr_domains), clean_up, 
               "%s", "Dimension mis-match!");

    matlib_complex* uptr = u_qnodes.elem_p;
    matlib_complex* vptr = v_qnodes.elem_p;
    
    matlib_complex iprod = 0;
    matlib_real jacob;

    for ( matlib_index i = 0; i< na.nr_domains; i++)
    {
        jacob = na.jacob[i];
        for (matlib_index j = 0; j< quadW.len; j++, uptr++, vptr++)
        {
            iprod += (*uptr * conj(*vptr) * quadW.elem_p[j] * jacob);
        }
    }
    debug_exit("Exit Status: %s", "SUCCESS" );
    return iprod;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

} /* fem1d_ziprod */ 

/*============================================================================*/

fem1d_err fem1d_quadP
/* Quadrature matrix for computing projection on vertex functions 
 * */ 
(
    matlib_xm  vphi,
    matlib_xv  quadW,
    matlib_xm* quadP
)
{
    debug_enter( "nr of quadrature points: %d", quadW.len);

    err_check(    (vphi.elem_p  == NULL) 
               || (quadW.elem_p == NULL), clean_up, 
               "%s", "Null pointer ecountered!");
    
    err_check( vphi.lenc != quadW.len, clean_up, 
               "Dimension mismatch (vphi: %d-by-%d, quadW: %d)!",
               vphi.lenc, vphi.lenr, quadW.len);


    fem1d_err error; 
    error = matlib_create_xm( FEM1D_NV, vphi.lenc, quadP, 
                              MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Memory allocation for projection matrix failed!");

    matlib_real* ptr1 = quadW.elem_p;
    matlib_real* ptr2 = vphi.elem_p;
    matlib_real* ptr3 = quadP->elem_p;

    for ( ; ptr1 < (quadW.elem_p + quadW.len);
            ptr1++, ptr2++, ptr3++)
    {
        *ptr3 = *ptr1 * *ptr2;
    }
    
    ptr1 = quadW.elem_p;

    for ( ; ptr1 < (quadW.elem_p + quadW.len);
            ptr1++, ptr2++, ptr3++)
    {
        *ptr3 = *ptr1 * *ptr2;
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM1D_SUCCESS;

clean_up:
    matlib_free(quadP->elem_p);
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM1D_FAILURE;

}/* End of fem1d_quadP */ 

 
/*============================================================================*/
fem1d_err fem1d_xprj
(
    fem1d_na  na,
    matlib_xv u_qnodes, /* values at the quadrature nodes */ 
    matlib_xm quadP,
    matlib_xv u_prj

)
{
    debug_enter( "nr of domains: %d, "
                 "quadP: %d-by-%d, ",
                 na.nr_domains, quadP.lenc, quadP.lenr);
    
    matlib_index mcnt = 0;
    err_check(    (na.elem_p       == NULL)    
               || (u_qnodes.elem_p == NULL) 
               || (quadP.elem_p    == NULL)  
               || (u_prj.elem_p    == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_qnodes.len != quadP.lenr * na.nr_domains) 
               || (quadP.lenc   != FEM1D_NV) 
               || (u_prj.len    != na.len), clean_up, 
               "Dimension mis-match "
               "(u_qnodes: %d, quadP: %d-by-%d, u_prj: %d, nr nodes: %d)!",
               u_qnodes.len, quadP.lenc, quadP.lenr, u_prj.len, na.len);

    matlib_xv u1_tmp = { .len    = quadP.lenr, 
                         .type   = MATLIB_COL_VECT,
                         .elem_p = u_qnodes.elem_p };

    matlib_xv u2_tmp;    
    fem1d_err error = matlib_create_xv( FEM1D_NV, &u2_tmp, MATLIB_COL_VECT);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Memory allocation for projection failed!");
    mcnt++;

    matlib_real beta  = 0.0;

    matlib_real* uprj_ptr0 = u_prj.elem_p;
    matlib_real* jptr = na.jacob;
    matlib_real* uprj_ptr;

    uprj_ptr0[FEM1D_INDEX_V1] = 0;
    for ( uprj_ptr = uprj_ptr0; 
          uprj_ptr < uprj_ptr0 + na.nr_domains; uprj_ptr++, jptr++)
    {
        error = matlib_xgemv( *jptr, quadP, u1_tmp, beta, u2_tmp );
        err_check( (error == FEM1D_FAILURE), clean_up, 
                   "%s", "Matrix multiplication failed!");
        
        (*(uprj_ptr + 0)) += u2_tmp.elem_p[FEM1D_INDEX_V1];

        (*(uprj_ptr + 1)) = u2_tmp.elem_p[FEM1D_INDEX_V2];

        u1_tmp.elem_p = u1_tmp.elem_p + quadP.lenr; 
    }


    matlib_free(u2_tmp.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:

    if (mcnt == 1)
        matlib_free(u2_tmp.elem_p);

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

} /* End of fem1d_xprj */ 

fem1d_err fem1d_periodic_xprj
(
    fem1d_na  na,
    matlib_xv u_qnodes, /* values at the quadrature nodes */ 
    matlib_xm quadP,
    matlib_xv u_prj
)
{
    debug_enter( "nr of domains: %d, "
                 "quadP: %d-by-%d, ",
                 na.nr_domains, quadP.lenc, quadP.lenr);
    
    matlib_index mcnt = 0;
    err_check(    (na.elem_p       == NULL)    
               || (u_qnodes.elem_p == NULL) 
               || (quadP.elem_p    == NULL)  
               || (u_prj.elem_p    == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_qnodes.len != quadP.lenr * na.nr_domains) 
               || (quadP.lenc   != FEM1D_NV) 
               || (u_prj.len    != na.nr_domains), clean_up, 
               "Dimension mis-match "
               "(u_qnodes: %d, quadP: %d-by-%d, u_prj: %d, nr nodes: %d)!",
               u_qnodes.len, quadP.lenc, quadP.lenr, u_prj.len, na.len);

    matlib_index INDEX_REAL  = 0;

    matlib_xv u1_tmp = { .len    = quadP.lenr, 
                         .type   = MATLIB_COL_VECT,
                         .elem_p = u_qnodes.elem_p };

    matlib_xv u2_tmp;    
    fem1d_err error = matlib_create_xv( FEM1D_NV, &u2_tmp, MATLIB_COL_VECT);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Memory allocation for projection failed!");
    mcnt++;

    matlib_real beta  = 0.0;

    matlib_real* uprj_ptr0 = u_prj.elem_p;
    matlib_real* jptr = na.jacob;
    matlib_real* uprj_ptr;

    uprj_ptr0[FEM1D_INDEX_V1] = 0;
    for ( uprj_ptr = uprj_ptr0; 
          uprj_ptr < uprj_ptr0 + (na.nr_domains-1) ; uprj_ptr++, jptr++)
    {
        error = matlib_xgemv( *jptr, quadP, u1_tmp, beta, u2_tmp );
        err_check( (error == FEM1D_FAILURE), clean_up, 
                   "%s", "Matrix multiplication failed!");
        
        uprj_ptr[FEM1D_INDEX_V1] += u2_tmp.elem_p[FEM1D_INDEX_V1];

        uprj_ptr[FEM1D_INDEX_V2] = u2_tmp.elem_p[FEM1D_INDEX_V2];

        u1_tmp.elem_p = u1_tmp.elem_p + quadP.lenr; 
    }

    error = matlib_xgemv( *jptr, quadP, u1_tmp, beta, u2_tmp );
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Matrix multiplication failed!");
    
    uprj_ptr[FEM1D_INDEX_V1] += u2_tmp.elem_p[FEM1D_INDEX_V1];

    uprj_ptr0[FEM1D_INDEX_V2] += u2_tmp.elem_p[FEM1D_INDEX_V2];

    matlib_free(u2_tmp.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:

    if (mcnt == 1)
        matlib_free(u2_tmp.elem_p);

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

} /* End of fem1d_xprj */ 


fem1d_err fem1d_zprj
(
    fem1d_na  na,
    matlib_zv u_qnodes, /* values at the quadrature nodes */ 
    matlib_xm quadP,
    matlib_zv u_prj

)
{
    debug_enter( "nr of domains: %d, "
                 "quadP: %d-by-%d, ",
                 na.nr_domains, quadP.lenc, quadP.lenr);
    
    matlib_index mcnt = 0;
    err_check(    (na.elem_p       == NULL)    
               || (u_qnodes.elem_p == NULL) 
               || (quadP.elem_p    == NULL)  
               || (u_prj.elem_p    == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_qnodes.len != quadP.lenr * na.nr_domains) 
               || (quadP.lenc   != FEM1D_NV) 
               || (u_prj.len    != na.len), clean_up, 
               "Dimension mis-match "
               "(u_qnodes: %d, quadP: %d-by-%d, u_prj: %d, nr nodes: %d)!",
               u_qnodes.len, quadP.lenc, quadP.lenr, u_prj.len, na.len);

    matlib_index INDEX_REAL  = 0;
    matlib_index INDEX_IMAG  = 1;
    matlib_index COMPLEX_DIM = 2; /* Dimension of a compelex plane */ 

    matlib_xm u1_tmp = { .lenc   = quadP.lenr, 
                         .lenr   = COMPLEX_DIM, 
                         .order  = MATLIB_ROW_MAJOR,
                         .op     = MATLIB_NO_TRANS,
                         .elem_p = (matlib_real*)u_qnodes.elem_p };

    matlib_xm u2_tmp;    
    fem1d_err error = matlib_create_xm( FEM1D_NV, COMPLEX_DIM, &u2_tmp, 
                                        MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Memory allocation for projection failed!");
    mcnt++;

    matlib_real beta  = 0.0;

    matlib_real* uprj_ptr0 = (matlib_real*)u_prj.elem_p; /* base address */ 
    matlib_real* uprj_ptr;
    matlib_real* jptr = na.jacob;

    uprj_ptr0[INDEX_REAL + COMPLEX_DIM * FEM1D_INDEX_V1] = 0;
    uprj_ptr0[INDEX_IMAG + COMPLEX_DIM * FEM1D_INDEX_V1] = 0;

    for ( uprj_ptr = uprj_ptr0; 
          uprj_ptr < (uprj_ptr0 + na.nr_domains * COMPLEX_DIM); 
          uprj_ptr += COMPLEX_DIM, jptr++)
    {
        error = matlib_xgemm( *jptr, quadP, u1_tmp, beta, u2_tmp );
        err_check( (error == FEM1D_FAILURE), clean_up, 
                   "%s", "Matrix multiplication failed!");
        
        uprj_ptr[INDEX_REAL + COMPLEX_DIM * FEM1D_INDEX_V1] 
            += u2_tmp.elem_p[INDEX_REAL + FEM1D_INDEX_V1 * COMPLEX_DIM];

        uprj_ptr[INDEX_IMAG + COMPLEX_DIM * FEM1D_INDEX_V1] 
            += u2_tmp.elem_p[INDEX_IMAG + FEM1D_INDEX_V1 * COMPLEX_DIM];

        uprj_ptr[INDEX_REAL + COMPLEX_DIM * FEM1D_INDEX_V2] 
            = u2_tmp.elem_p[INDEX_REAL + FEM1D_INDEX_V2 * COMPLEX_DIM];

        uprj_ptr[INDEX_IMAG + COMPLEX_DIM * FEM1D_INDEX_V2] 
            += u2_tmp.elem_p[INDEX_IMAG + FEM1D_INDEX_V2 * COMPLEX_DIM];

        u1_tmp.elem_p = u1_tmp.elem_p + quadP.lenr; 
    }


    matlib_free(u2_tmp.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:

    if (mcnt == 1)
        matlib_free(u2_tmp.elem_p);

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

} /* End of fem1d_zprj */ 

fem1d_err fem1d_periodic_zprj
(
    fem1d_na  na,
    matlib_zv u_qnodes, /* values at the quadrature nodes */ 
    matlib_xm quadP,
    matlib_zv u_prj
)
{
    debug_enter( "nr of domains: %d, "
                 "quadP: %d-by-%d, ",
                 na.nr_domains, quadP.lenc, quadP.lenr);
    
    matlib_index mcnt = 0;
    err_check(    (na.elem_p       == NULL)    
               || (u_qnodes.elem_p == NULL) 
               || (quadP.elem_p    == NULL)  
               || (u_prj.elem_p    == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_qnodes.len != quadP.lenr * na.nr_domains) 
               || (quadP.lenc   != FEM1D_NV) 
               || (u_prj.len    != na.nr_domains), clean_up, 
               "Dimension mis-match "
               "(u_qnodes: %d, quadP: %d-by-%d, u_prj: %d, nr nodes: %d)!",
               u_qnodes.len, quadP.lenc, quadP.lenr, u_prj.len, na.len);

    matlib_index INDEX_REAL  = 0;
    matlib_index INDEX_IMAG  = 1;
    matlib_index COMPLEX_DIM = 2; /* Dimension of a compelex plane */ 

    matlib_xm u1_tmp = { .lenc   = quadP.lenr, 
                         .lenr   = COMPLEX_DIM, 
                         .order  = MATLIB_ROW_MAJOR,
                         .op     = MATLIB_NO_TRANS,
                         .elem_p = (matlib_real*)u_qnodes.elem_p };

    matlib_xm u2_tmp;    
    fem1d_err error = matlib_create_xm( FEM1D_NV, COMPLEX_DIM, &u2_tmp, 
                                        MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Memory allocation for projection failed!");
    mcnt++;

    matlib_real beta  = 0.0;

    matlib_real* uprj_ptr0 = (matlib_real*)u_prj.elem_p; /* base address */ 
    matlib_real* uprj_ptr;
    matlib_real* jptr = na.jacob;

    uprj_ptr0[INDEX_REAL + COMPLEX_DIM * FEM1D_INDEX_V1] = 0;
    uprj_ptr0[INDEX_IMAG + COMPLEX_DIM * FEM1D_INDEX_V1] = 0;

    for ( uprj_ptr = uprj_ptr0; 
          uprj_ptr < (uprj_ptr0 + (na.nr_domains - 1) * COMPLEX_DIM); 
          uprj_ptr += COMPLEX_DIM, jptr++)
    {
        error = matlib_xgemm( *jptr, quadP, u1_tmp, beta, u2_tmp );
        err_check( (error == FEM1D_FAILURE), clean_up, 
                   "%s", "Matrix multiplication failed!");
        
        uprj_ptr[INDEX_REAL + COMPLEX_DIM * FEM1D_INDEX_V1] 
            += u2_tmp.elem_p[INDEX_REAL + FEM1D_INDEX_V1 * COMPLEX_DIM];

        uprj_ptr[INDEX_IMAG + COMPLEX_DIM * FEM1D_INDEX_V1] 
            += u2_tmp.elem_p[INDEX_IMAG + FEM1D_INDEX_V1 * COMPLEX_DIM];

        uprj_ptr[INDEX_REAL + COMPLEX_DIM * FEM1D_INDEX_V2] 
            = u2_tmp.elem_p[INDEX_REAL + FEM1D_INDEX_V2 * COMPLEX_DIM];

        uprj_ptr[INDEX_IMAG + COMPLEX_DIM * FEM1D_INDEX_V2] 
            = u2_tmp.elem_p[INDEX_IMAG + FEM1D_INDEX_V2 * COMPLEX_DIM];

        u1_tmp.elem_p = u1_tmp.elem_p + quadP.lenr; 
    }

    error = matlib_xgemm( *jptr, quadP, u1_tmp, beta, u2_tmp );
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Matrix multiplication failed!");
    
    uprj_ptr[INDEX_REAL + COMPLEX_DIM * FEM1D_INDEX_V1] 
        += u2_tmp.elem_p[INDEX_REAL + FEM1D_INDEX_V1 * COMPLEX_DIM];

    uprj_ptr[INDEX_IMAG + COMPLEX_DIM * FEM1D_INDEX_V1] 
        += u2_tmp.elem_p[INDEX_IMAG + FEM1D_INDEX_V1 * COMPLEX_DIM];

    uprj_ptr0[INDEX_REAL + COMPLEX_DIM * FEM1D_INDEX_V1] 
        = u2_tmp.elem_p[INDEX_REAL + FEM1D_INDEX_V2 * COMPLEX_DIM];

    uprj_ptr0[INDEX_IMAG + COMPLEX_DIM * FEM1D_INDEX_V1] 
        = u2_tmp.elem_p[INDEX_IMAG + FEM1D_INDEX_V2 * COMPLEX_DIM];


    matlib_free(u2_tmp.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:

    if (mcnt == 1)
        matlib_free(u2_tmp.elem_p);

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

} /* End of fem1d_periodic_zprj */ 

/*============================================================================*/
fem1d_err fem1d_NB_xprj
/* Projection from nodal basis representation */ 
(
    fem1d_na  na,
    matlib_xv u_nodes, /* values at the nodes */ 
    matlib_xv u_prj
)
{
    debug_enter( "nr of domains: %d", na.nr_domains);
    
    err_check(    (na.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL)
               || (u_prj.elem_p   == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != na.len) 
               || (u_prj.len   != na.len), clean_up, 
               "Dimension mis-match (u_nodes: %d, u_prj: %d, nr nodes: %d)!",
               u_nodes.len, u_prj.len, na.len);

    matlib_real* uptr      = u_nodes.elem_p;
    matlib_real* uprj_ptr0 = u_prj.elem_p;
    matlib_real* uprj_ptr;
    matlib_real* jptr = na.jacob;
    matlib_real f1, f2;

    uprj_ptr0[FEM1D_INDEX_V1] = 0;
    
    for ( uprj_ptr = uprj_ptr0; 
          uprj_ptr < uprj_ptr0 + na.nr_domains; uprj_ptr++, jptr++, uptr++)
    {
        
        f1 = uptr[FEM1D_INDEX_V1] * *jptr;
        f2 = uptr[FEM1D_INDEX_V2] * *jptr;

        
        uprj_ptr[FEM1D_INDEX_V1] += (   FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V1] * f1
                                      + FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V2] * f2);
        
        uprj_ptr[FEM1D_INDEX_V2]  = (   FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V1] * f1
                                      + FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V2] * f2);

    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

} /* End of fem1d_NB_xprj */ 


fem1d_err fem1d_periodic_NB_xprj
/* Projection from nodal basis representation */ 
(
    fem1d_na  na,
    matlib_xv u_nodes, /* values at the nodes */ 
    matlib_xv u_prj
)
{
    debug_enter( "nr of domains: %d", na.nr_domains);
    
    err_check(    (na.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL)
               || (u_prj.elem_p   == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != na.len) 
               || (u_prj.len   != na.nr_domains), clean_up, 
               "Dimension mis-match (u_nodes: %d, u_prj: %d, nr nodes: %d)!",
               u_nodes.len, u_prj.len, na.len);

    matlib_real* uptr      = u_nodes.elem_p;
    matlib_real* uprj_ptr0 = u_prj.elem_p;
    matlib_real* uprj_ptr;
    matlib_real* jptr = na.jacob;
    matlib_real f1, f2;

    uprj_ptr0[FEM1D_INDEX_V1] = 0;

    for ( uprj_ptr = uprj_ptr0; 
          uprj_ptr < uprj_ptr0 + (na.nr_domains - 1); uprj_ptr++, uptr++, jptr++)
    {
        
        f1 = uptr[FEM1D_INDEX_V1] * *jptr;
        f2 = uptr[FEM1D_INDEX_V2] * *jptr;

        
        uprj_ptr[FEM1D_INDEX_V1] += (   FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V1] * f1
                                      + FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V2] * f2);
        
        uprj_ptr[FEM1D_INDEX_V2]  = (   FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V1] * f1
                                      + FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V2] * f2);

    }
    f1 = uptr[FEM1D_INDEX_V1] * *jptr;
    f2 = uptr[FEM1D_INDEX_V2] * *jptr;

    
    uprj_ptr[FEM1D_INDEX_V1] += (   FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V1] * f1
                                  + FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V2] * f2);
    
    uprj_ptr0[FEM1D_INDEX_V1] += (   FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V1] * f1
                                   + FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V2] * f2);


    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

} /* End of fem1d_periodic_NB_xprj */ 


fem1d_err fem1d_NB_zprj
(
    fem1d_na  na,
    matlib_zv u_nodes, /* values at the nodes */ 
    matlib_zv u_prj
)
{
    debug_enter( "nr of domains: %d", na.nr_domains);
    
    err_check(    (na.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL)
               || (u_prj.elem_p   == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != na.len) 
               || (u_prj.len   != na.len), clean_up, 
               "Dimension mis-match (u_nodes: %d, u_prj: %d, nr nodes: %d)!",
               u_nodes.len, u_prj.len, na.len);
    
    matlib_complex* uptr      = u_nodes.elem_p;
    matlib_complex* uprj_ptr0 = u_prj.elem_p;
    matlib_complex* uprj_ptr;

    matlib_complex f1, f2;
    matlib_real* jptr = na.jacob;
    
    uprj_ptr0[FEM1D_INDEX_V1] = 0;

    for ( uprj_ptr = uprj_ptr0; 
          uprj_ptr < uprj_ptr0 + na.nr_domains; uprj_ptr++, uptr++, jptr++)
    {
        f1 = uptr[FEM1D_INDEX_V1] * *jptr;
        f2 = uptr[FEM1D_INDEX_V2] * *jptr;

        uprj_ptr[FEM1D_INDEX_V1] += (   FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V1] * f1
                                      + FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V2] * f2);
        
        uprj_ptr[FEM1D_INDEX_V2]  = (   FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V1] * f1
                                      + FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V2] * f2);

    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

} /* End of fem1d_NB_zprj */ 

fem1d_err fem1d_periodic_NB_zprj
(
    fem1d_na  na,
    matlib_zv u_nodes, /* values at the nodes */ 
    matlib_zv u_prj
)
{
    debug_enter( "nr of domains: %d", na.nr_domains);
    
    err_check(    (na.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL)
               || (u_prj.elem_p   == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != na.len) 
               || (u_prj.len   != na.nr_domains), clean_up, 
               "Dimension mis-match (u_nodes: %d, u_prj: %d, nr nodes: %d)!",
               u_nodes.len, u_prj.len, na.len);
    
    matlib_complex* uptr      = u_nodes.elem_p;
    matlib_complex* uprj_ptr0 = u_prj.elem_p;
    matlib_complex* uprj_ptr;

    matlib_complex f1, f2;
    matlib_real* jptr = na.jacob;

    uprj_ptr0[FEM1D_INDEX_V1] = 0;

    for ( uprj_ptr = uprj_ptr0; 
          uprj_ptr < uprj_ptr0 + (na.nr_domains - 1); uprj_ptr++, uptr++, jptr++)
    {
        f1 = uptr[FEM1D_INDEX_V1] * *jptr;
        f2 = uptr[FEM1D_INDEX_V2] * *jptr;

        uprj_ptr[FEM1D_INDEX_V1] += (   FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V1] * f1
                                      + FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V2] * f2);
        
        uprj_ptr[FEM1D_INDEX_V2]  = (   FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V1] * f1
                                      + FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V2] * f2);

    }
    f1 = uptr[FEM1D_INDEX_V1] * *jptr;
    f2 = uptr[FEM1D_INDEX_V2] * *jptr;

    uprj_ptr[FEM1D_INDEX_V1] += (   FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V1] * f1
                                  + FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V2] * f2);
    
    uprj_ptr0[FEM1D_INDEX_V1] += (   FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V1] * f1
                                   + FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V2] * f2);


    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

} /* End of fem1d_periodic_NB_zprj */ 
/*============================================================================*/

matlib_real fem1d_NB_xnormL2
(
    fem1d_na  na,
    matlib_xv u_nodes /* values at the nodes */ 
)
{
    debug_enter( "nr of domains: %d", na.nr_domains);
    
    err_check(    (na.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != na.len), clean_up, 
               "Dimension mis-match (u_nodes: %d, nr nodes: %d)!",
               u_nodes.len, na.len);

    matlib_real f1, f2;
    matlib_real f11, f12, f22;

    matlib_real sum = 0;
    matlib_real* jptr = na.jacob;

    matlib_real* uptr = u_nodes.elem_p;
    
    for ( ; uptr < u_nodes.elem_p + na.nr_domains; uptr++, jptr++)
    {
        f1 = uptr[FEM1D_INDEX_V1];
        f2 = uptr[FEM1D_INDEX_V2];

        f11 = f1 * f1 * *jptr;
        f22 = f2 * f2 * *jptr;

        f12 = 2 * (f1 * f2) * *jptr;
        
        sum += (   FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V1] * f11
                 + FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V2] * f22 
                 + FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V2] * f12 );

    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return sqrt(sum);

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

} /* End of fem1d_NB_xnormL2 */ 


matlib_real fem1d_NB_znormL2
(
    fem1d_na  na,
    matlib_zv u_nodes /* values at the nodes */ 
)
{
    debug_enter( "nr of domains: %d", na.nr_domains);
    
    err_check(    (na.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != na.len), clean_up, 
               "Dimension mis-match (u_nodes: %d, nr nodes: %d)!",
               u_nodes.len, na.len);
    
    matlib_complex f1, f2;
    matlib_real f11, f12, f22;

    matlib_real sum = 0;
    matlib_real* jptr = na.jacob;

    matlib_complex* uptr = u_nodes.elem_p;
    
    for ( ; uptr < u_nodes.elem_p + na.nr_domains; uptr++, jptr++)
    {
        f1 = uptr[FEM1D_INDEX_V1];
        f2 = uptr[FEM1D_INDEX_V2];

        f11 = (matlib_real) (f1 * conj(f1) * *jptr);
        f22 = (matlib_real) (f2 * conj(f2) * *jptr);

        f12 = (matlib_real) (f1 * conj(f2) + f2 * conj(f1)) * *jptr;
        
        sum += (   FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V1] * f11
                 + FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V2] * f22 
                 + FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V2] * f12 );

    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return sqrt(sum);

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

} /* End of fem1d_NB_znormL2 */ 



/*============================================================================*/
matlib_real fem1d_NB_xiprod
(
    fem1d_na  na,
    matlib_xv u_nodes, /* values at the nodes */ 
    matlib_xv v_nodes
)
{
    debug_enter( "nr of domains: %d", na.nr_domains);
    
    err_check(    (na.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL)
               || (v_nodes.elem_p == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != na.len) 
               || (v_nodes.len != na.len), clean_up, 
               "Dimension mis-match (u_nodes: %d, v_nodes: %d, nr nodes: %d)!",
               u_nodes.len, v_nodes.len, na.len);

    matlib_real* uptr = u_nodes.elem_p;
    matlib_real* vptr = v_nodes.elem_p;

    matlib_real* jptr = na.jacob;
    matlib_real f1, f2;
    matlib_real g1, g2;
    matlib_real fg11, fg22, fg12;

    matlib_real sum = 0;
    for ( ; uptr < u_nodes.elem_p + na.nr_domains; 
            uptr++, vptr++, jptr++)
    {

        f1 = uptr[FEM1D_INDEX_V1];
        f2 = uptr[FEM1D_INDEX_V2];

        g1 = vptr[FEM1D_INDEX_V1];
        g2 = vptr[FEM1D_INDEX_V2];

        fg11 = f1 * g1 * *jptr;
        fg22 = f2 * g2 * *jptr;

        fg12 = (f1 * g2 + f2 * g1) * *jptr;
        
        sum += (   FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V1] * fg11
                 + FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V2] * fg22 
                 + FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V2] * fg12 );
    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return sum;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

}/* End of fem1d_NB_xiprod */ 

matlib_complex fem1d_NB_ziprod
(
    fem1d_na  na,
    matlib_zv u_nodes, /* values at the nodes */ 
    matlib_zv v_nodes
)
{
    debug_enter( "nr of domains: %d", na.nr_domains);
    
    err_check(    (na.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL)
               || (v_nodes.elem_p == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != na.len) 
               || (v_nodes.len != na.len), clean_up, 
               "Dimension mis-match (u_nodes: %d, v_nodes: %d, nr nodes: %d)!",
               u_nodes.len, v_nodes.len, na.len);
    
    matlib_complex* uptr = u_nodes.elem_p;
    matlib_complex* vptr = v_nodes.elem_p;

    matlib_real* jptr = na.jacob;
    matlib_complex f1, f2;
    matlib_complex g1, g2;
    matlib_complex fg11, fg22, fg12;

    matlib_complex sum = 0;
    for ( ; uptr < u_nodes.elem_p + na.nr_domains; 
            uptr++, vptr++, jptr++)
    {

        f1 = uptr[FEM1D_INDEX_V1];
        f2 = uptr[FEM1D_INDEX_V2];

        g1 = vptr[FEM1D_INDEX_V1];
        g2 = vptr[FEM1D_INDEX_V2];

        fg11 = f1 * conj(g1) * *jptr;
        fg22 = f2 * conj(g2) * *jptr;

        fg12 = (f1 * conj(g2) + f2 * conj(g1)) * *jptr;
        
        sum += (   FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V1] * fg11
                 + FEM1D_MEMI[FEM1D_INDEX_V2][FEM1D_INDEX_V2] * fg22 
                 + FEM1D_MEMI[FEM1D_INDEX_V1][FEM1D_INDEX_V2] * fg12 );
    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return sum;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

}/* End of fem1d_NB_ziprod */ 

/*============================================================================*/
fem1d_err fem1d_quadM
(
    matlib_xm  vphi,
    matlib_xv  quadW,
    matlib_xm* Q
)
{
    debug_enter( "nr of quadrature points: %d", quadW.len);
    
    matlib_index mcnt = 0;
    err_check(    (vphi.elem_p  == NULL)    
               || (quadW.elem_p == NULL), clean_up, 
               "%s", "Null pointer(s) ecountered!");

    err_check(    (vphi.lenc != quadW.len) 
               || (vphi.lenr != FEM1D_NV), clean_up, 
               "Dimension mis-match (vphi: %d-by-%d, quadW: %d)!",
               vphi.lenc, vphi.lenr, quadW.len);

    fem1d_err error = matlib_create_xm( FEM1D_NR_COMBI, quadW.len, Q, 
                                        MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Memory allocation for quadrature matrix failed!");
    mcnt++;
    matlib_index i;
    matlib_index stride = quadW.len;

    for (i = 0; i < quadW.len; i++)
    {
        Q->elem_p[FEM1D_INDEX_V11 * stride + i] =   vphi.elem_p[i + FEM1D_INDEX_V1 * stride] 
                                                  * vphi.elem_p[i + FEM1D_INDEX_V1 * stride]
                                                  * quadW.elem_p[i];

        Q->elem_p[FEM1D_INDEX_V12 * stride + i] =   vphi.elem_p[i + FEM1D_INDEX_V1 * stride] 
                                                  * vphi.elem_p[i + FEM1D_INDEX_V2 * stride]
                                                  * quadW.elem_p[i];
        
        Q->elem_p[FEM1D_INDEX_V22 * stride + i] =   vphi.elem_p[i + FEM1D_INDEX_V2 * stride] 
                                                  * vphi.elem_p[i + FEM1D_INDEX_V2 * stride]
                                                  * quadW.elem_p[i];

    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    
    if (mcnt == 1)
        matlib_free(Q->elem_p);

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;
}

/*============================================================================*/
matlib_index fem1d_get_nnz(fem1d_na na)
{
    matlib_index nnz = 2 * na.nr_domains + 1;
    return nnz;
}

matlib_index fem1d_periodic_get_nnz(fem1d_na na)
{
    matlib_index nnz = 2 * na.nr_domains;
    return nnz;
}

/*============================================================================*/

fem1d_err fem1d_GMMSparsity
(
    fem1d_na      na,
    matlib_index* row, /* required to be pre-allocated */                     
    matlib_index* col  /* required to be pre-allocated */                  
)
{
    err_check(    (row == NULL) 
               || (col == NULL),  clean_up,
               "%s", "Null pointer(s) encountered!");

    /* length of col: nnz
     * length of row: nr_nodes + 1
     *
     * */ 
    matlib_index s0, s = 0;
    matlib_index i;

    for (i = 0; i < na.len - 1; i++)
    {
        row[i] = s;
        col[s] = i; 
        s++;

        col[s] = i+1; 
        s++;
    }
    row[i] = s;
    col[s] = i; 
    s++;
    i++;
    debug_body( "nr non-zero entries: %d", s);
    row[i] = s;

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;
}

fem1d_err fem1d_periodic_GMMSparsity
(
    fem1d_na      na,
    matlib_index* row, /* required to be pre-allocated */                     
    matlib_index* col  /* required to be pre-allocated */                  
)
{
    err_check(    (row == NULL) 
               || (col == NULL),  clean_up,
               "%s", "Null pointer(s) encountered!");

    /* length of col: nnz
     * length of row: nr_nodes + 1
     *
     * */ 
    matlib_index s0, s = 0;
    matlib_index i;

    row[i] = s;
    col[s] = i; 
    s++;

    col[s] = i+1; 
    s++;

    col[s] = na.nr_domains; 
    s++;
    for (i = 0; i < na.nr_domains - 1; i++)
    {
        row[i] = s;
        col[s] = i; 
        s++;

        col[s] = i+1; 
        s++;
    }
    row[i] = s;
    col[s] = i; 
    s++;
    i++;
    debug_body( "nr non-zero entries: %d", s);
    row[i] = s;

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;
}

fem1d_err fem1d_XCSRGMM
(
    fem1d_na      na,
    matlib_xv     q,
    matlib_index* row,                     
    matlib_real*  ugpmm                   
)
{
    debug_enter("nr nodes: %d, q: %d", na.len, q.len);

    err_check(    (na.elem_p   == NULL) 
               || (q.elem_p    == NULL), clean_up, 
               "%s", "Null pointer(s) encountered!");

    err_check( na.nr_domains * FEM1D_NR_COMBI != q.len, clean_up,
               "%s: Dimension mistmatch (q: %d, required size: %d)!",
               q.len, na.nr_domains * FEM1D_NR_COMBI);

    matlib_real* qptr = q.elem_p;
    
    matlib_int i, j;
    ugpmm[0] = 0;
    
    for (i = 0; i < na.nr_domains; i++, qptr += FEM1D_NR_COMBI)
    {
        j = row[i];
        ugpmm[ j + 0] += qptr[FEM1D_INDEX_V11];
        ugpmm[ j + 1] = qptr[FEM1D_INDEX_V12];

        j = row[i+1];
        ugpmm[ j + 0] = qptr[FEM1D_INDEX_V22];
    }

    /* === */ 

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;
}

fem1d_err fem1d_periodic_XCSRGMM
(
    fem1d_na      na,
    matlib_xv     q,
    matlib_index* row,                     
    matlib_real*  ugpmm                   
)
{
    debug_enter("nr nodes: %d, q: %d", na.len, q.len);

    err_check(    (na.elem_p   == NULL) 
               || (q.elem_p    == NULL), clean_up, 
               "%s", "Null pointer(s) encountered!");

    err_check( na.nr_domains * FEM1D_NR_COMBI != q.len, clean_up,
               "%s: Dimension mistmatch (q: %d, required size: %d)!",
               q.len, na.nr_domains * FEM1D_NR_COMBI);

    matlib_real* qptr = q.elem_p;
    
    matlib_int i, j;
    ugpmm[0] = 0;

    for (i = 0; i < na.nr_domains - 1; i++, qptr += FEM1D_NR_COMBI)
    {
        j = row[i];
        ugpmm[j + 0] += qptr[FEM1D_INDEX_V11];
        ugpmm[j + 1] = qptr[FEM1D_INDEX_V12];

        j = row[i + 1];
        ugpmm[j + 0] = qptr[FEM1D_INDEX_V22];
    }
    j = row[i];
    ugpmm[j + 0] += qptr[FEM1D_INDEX_V11];

    ugpmm[0 + 2] = qptr[FEM1D_INDEX_V12];
    ugpmm[0 + 0] += qptr[FEM1D_INDEX_V22];

    /* === */ 

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;
}


fem1d_err fem1d_ZCSRGMM
(
    fem1d_na      na,
    matlib_zv     q,
    matlib_index* row,                     
    matlib_complex* ugpmm                   
)
{
    debug_enter("nr nodes: %d, q: %d", na.len, q.len);

    err_check(    (na.elem_p   == NULL) 
               || (q.elem_p    == NULL), clean_up, 
               "%s", "Null pointer(s) encountered!");

    err_check( na.nr_domains * FEM1D_NR_COMBI != q.len, clean_up,
               "%s: Dimension mistmatch (q: %d, required size: %d)!",
               q.len, na.nr_domains * FEM1D_NR_COMBI);

    matlib_complex* qptr = q.elem_p;
    
    matlib_int i, j;
    ugpmm[0] = 0;
    
    for (i = 0; i < na.nr_domains; i++, qptr += FEM1D_NR_COMBI)
    {
        j = row[i];
        ugpmm[ j + 0] += qptr[FEM1D_INDEX_V11];
        ugpmm[ j + 1] = qptr[FEM1D_INDEX_V12];

        j = row[i+1];
        ugpmm[ j + 0] = qptr[FEM1D_INDEX_V22];
    }

    /* === */ 

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;
}

fem1d_err fem1d_periodic_ZCSRGMM
(
    fem1d_na      na,
    matlib_zv     q,
    matlib_index* row,                     
    matlib_complex* ugpmm                   
)
{
    debug_enter("nr nodes: %d, q: %d", na.len, q.len);

    err_check(    (na.elem_p   == NULL) 
               || (q.elem_p    == NULL), clean_up, 
               "%s", "Null pointer(s) encountered!");

    err_check( na.nr_domains * FEM1D_NR_COMBI != q.len, clean_up,
               "%s: Dimension mistmatch (q: %d, required size: %d)!",
               q.len, na.nr_domains * FEM1D_NR_COMBI);

    matlib_complex* qptr = q.elem_p;
    
    matlib_int i, j;
    ugpmm[0] = 0;

    for (i = 0; i < na.nr_domains - 1; i++, qptr += FEM1D_NR_COMBI)
    {
        j = row[i];
        ugpmm[j + 0] += qptr[FEM1D_INDEX_V11];
        ugpmm[j + 1] = qptr[FEM1D_INDEX_V12];

        j = row[i + 1];
        ugpmm[j + 0] = qptr[FEM1D_INDEX_V22];
    }
    j = row[i];
    ugpmm[j + 0] += qptr[FEM1D_INDEX_V11];

    ugpmm[0 + 2] = qptr[FEM1D_INDEX_V12];
    ugpmm[0 + 0] += qptr[FEM1D_INDEX_V22];

    /* === */ 

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;
}


/*============================================================================*/

fem1d_err fem1d_xm_sparse_GMM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xm         Q,
    matlib_xv         phi,
    matlib_xm_sparse* M
)
{
    debug_enter( "size of Q: %d-by-%d), "
                 "length of phi: %d",
                 Q.lenc, Q.lenr, phi.len);
    
    matlib_index i, mcnt = 0;
    err_check( Q.lenc != FEM1D_NR_COMBI, clean_up,
               "Dimension mismatch (Q: %d-by-%d, nr combi: %d)!", 
               Q.lenc, Q.lenr, FEM1D_NR_COMBI);

    err_check( (Q.lenr * na.nr_domains != phi.len), clean_up,
               "Dimension mismatch (phi: %d, require size: %d)!", 
               phi.len, Q.lenr * na.nr_domains);

    matlib_index nnz = fem1d_get_nnz(na);
    debug_body( "nr. of non-zero elements: %d", nnz);
    err_check( nnz == 0, clean_up, "%s", 
               "Number of non-zero elements of the mass-matrix is inconsistent!");

    fem1d_err error;
    
    matlib_xv q;
    error = matlib_create_xv( Q.lenc * na.nr_domains, &q, MATLIB_COL_VECT);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the inner product array!");
    mcnt++; /* 1 */ 
    matlib_real alpha, beta = 0.0;
    phi.len = Q.lenr;

    matlib_xv q_tmp = { .len = Q.lenc, 
                        .type = MATLIB_COL_VECT,
                        .elem_p = q.elem_p};

    for ( i = 0; i < na.nr_domains; 
          i++, (phi.elem_p) += (Q.lenr), (q_tmp.elem_p) += (Q.lenc))
    {
        alpha = na.jacob[i];
        error = matlib_xgemv(alpha, Q, phi, beta, q_tmp );
        err_check( (error == FEM1D_FAILURE), clean_up, 
                   "Failed to compute inner-products (domain index: %d)!", i);
    }
    
    error = matlib_create_xm_sparse( na.len, na.len, nnz, M, MATLIB_CSR); 
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the row array!");
    mcnt++; /* 2 */ 

    error = fem1d_GMMSparsity( na, M->rowIn, M->colIn);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to get the sparsity structure of the mass-matrix!");

    error = fem1d_XCSRGMM( na, q, M->rowIn, M->elem_p);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to compute the mass-matrix!");

    
    matlib_free(q.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    if (mcnt == 2)
    {
        matlib_free(M->elem_p);
        matlib_free(M->colIn);
        matlib_free(M->rowIn);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(q.elem_p);
    }

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

}

fem1d_err fem1d_periodic_xm_sparse_GMM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xm         Q,
    matlib_xv         phi,
    matlib_xm_sparse* M
)
{
    debug_enter( "size of Q: %d-by-%d), "
                 "length of phi: %d",
                 Q.lenc, Q.lenr, phi.len);
    
    matlib_index i, mcnt = 0;
    err_check( Q.lenc != FEM1D_NR_COMBI, clean_up,
               "Dimension mismatch (Q: %d-by-%d, nr combi: %d)!", 
               Q.lenc, Q.lenr, FEM1D_NR_COMBI);

    err_check( (Q.lenr * na.nr_domains != phi.len), clean_up,
               "Dimension mismatch (phi: %d, require size: %d)!", 
               phi.len, Q.lenr * na.nr_domains);

    matlib_index nnz = fem1d_periodic_get_nnz(na);
    debug_body( "nr. of non-zero elements: %d", nnz);
    err_check( nnz == 0, clean_up, "%s", 
               "Number of non-zero elements of the mass-matrix is inconsistent!");

    fem1d_err error;
    
    matlib_xv q;
    error = matlib_create_xv( Q.lenc * na.nr_domains, &q, MATLIB_COL_VECT);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the inner product array!");
    mcnt++; /* 1 */ 
    matlib_real alpha, beta = 0.0;
    phi.len = Q.lenr;

    matlib_xv q_tmp = { .len = Q.lenc, 
                        .type = MATLIB_COL_VECT,
                        .elem_p = q.elem_p};

    for ( i = 0; i < na.nr_domains; 
          i++, (phi.elem_p) += (Q.lenr), (q_tmp.elem_p) += (Q.lenc))
    {
        alpha = na.jacob[i];
        error = matlib_xgemv(alpha, Q, phi, beta, q_tmp );
        err_check( (error == FEM1D_FAILURE), clean_up, 
                   "Failed to compute inner-products (domain index: %d)!", i);
    }
    
    error = matlib_create_xm_sparse( na.nr_domains, na.nr_domains, nnz, M, MATLIB_CSR); 
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the row array!");
    mcnt++; /* 2 */ 

    error = fem1d_periodic_GMMSparsity( na, M->rowIn, M->colIn);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to get the sparsity structure of the mass-matrix!");

    error = fem1d_periodic_XCSRGMM( na, q, M->rowIn, M->elem_p);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to compute the mass-matrix!");

    
    matlib_free(q.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    if (mcnt == 2)
    {
        matlib_free(M->elem_p);
        matlib_free(M->colIn);
        matlib_free(M->rowIn);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(q.elem_p);
    }

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

}

fem1d_err fem1d_zm_sparse_GMM
/* Complex - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xm         Q,
    matlib_zv         phi,
    matlib_zm_sparse* M
)
{
    debug_enter( "size of Q: %d-by-%d), "
                 "length of phi: %d",
                 Q.lenc, Q.lenr, phi.len);
    
    matlib_index i, mcnt = 0;
    err_check( Q.lenc != FEM1D_NR_COMBI, clean_up,
               "Dimension mismatch (Q: %d-by-%d, nr combi: %d)!", 
               Q.lenc, Q.lenr, FEM1D_NR_COMBI);

    err_check( (Q.lenr * na.nr_domains != phi.len), clean_up,
               "Dimension mismatch (phi: %d, require size: %d)!", 
               phi.len, Q.lenr * na.nr_domains);

    matlib_index nnz = fem1d_get_nnz(na);
    debug_body( "nr. of non-zero elements: %d", nnz);
    err_check( nnz == 0, clean_up, "%s", 
               "Number of non-zero elements of the mass-matrix is inconsistent!");

    fem1d_err error;
    
    matlib_zv q;
    error = matlib_create_zv( Q.lenc * na.nr_domains, &q, MATLIB_COL_VECT);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the inner product array!");
    mcnt++; /* 1 */ 
    matlib_real alpha, beta = 0.0;

    matlib_index COMPLEX_DIM = 2;
    matlib_xm q_tmp = { .lenc   = Q.lenc, 
                        .lenr   = COMPLEX_DIM,
                        .order  = MATLIB_ROW_MAJOR,
                        .op     = MATLIB_NO_TRANS,
                        .elem_p = (matlib_real*)q.elem_p};
    matlib_xm phi_tmp = { .lenc   = Q.lenr, 
                          .lenr   = COMPLEX_DIM,
                          .order  = MATLIB_ROW_MAJOR,
                          .op     = MATLIB_NO_TRANS,
                          .elem_p = (matlib_real*)phi.elem_p};

    for ( i = 0; i < na.nr_domains; 
          i++, 
          (phi_tmp.elem_p) += (COMPLEX_DIM * Q.lenr), 
          (q_tmp.elem_p)   += (COMPLEX_DIM * Q.lenc))
    {
        alpha = na.jacob[i];
        error = matlib_xgemm(alpha, Q, phi_tmp, beta, q_tmp );
        err_check( (error == FEM1D_FAILURE), clean_up, 
                   "Failed to compute inner-product (domain index: %d)!", i);
    }
    
    error = matlib_create_zm_sparse( na.len, na.len, nnz, M, MATLIB_CSR); 
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the row array!");
    mcnt++; /* 2 */ 

    error = fem1d_GMMSparsity( na, M->rowIn, M->colIn);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to get the sparsity structure of the mass-matrix!");

    error = fem1d_ZCSRGMM( na, q, M->rowIn, M->elem_p);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to compute the mass-matrix!");

    
    matlib_free(q.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    if (mcnt == 2)
    {
        matlib_free(M->elem_p);
        matlib_free(M->colIn);
        matlib_free(M->rowIn);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(q.elem_p);
    }

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

}

fem1d_err fem1d_periodic_zm_sparse_GMM
/* Complex - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xm         Q,
    matlib_zv         phi,
    matlib_zm_sparse* M
)
{
    debug_enter( "size of Q: %d-by-%d), "
                 "length of phi: %d",
                 Q.lenc, Q.lenr, phi.len);
    
    matlib_index i, mcnt = 0;
    err_check( Q.lenc != FEM1D_NR_COMBI, clean_up,
               "Dimension mismatch (Q: %d-by-%d, nr combi: %d)!", 
               Q.lenc, Q.lenr, FEM1D_NR_COMBI);

    err_check( (Q.lenr * na.nr_domains != phi.len), clean_up,
               "Dimension mismatch (phi: %d, require size: %d)!", 
               phi.len, Q.lenr * na.nr_domains);

    matlib_index nnz = fem1d_periodic_get_nnz(na);
    debug_body( "nr. of non-zero elements: %d", nnz);
    err_check( nnz == 0, clean_up, "%s", 
               "Number of non-zero elements of the mass-matrix is inconsistent!");

    fem1d_err error;
    
    matlib_zv q;
    error = matlib_create_zv( Q.lenc * na.nr_domains, &q, MATLIB_COL_VECT);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the inner product array!");
    mcnt++; /* 1 */ 
    matlib_real alpha, beta = 0.0;

    matlib_index COMPLEX_DIM = 2;
    matlib_xm q_tmp = { .lenc   = Q.lenc, 
                        .lenr   = COMPLEX_DIM,
                        .order  = MATLIB_ROW_MAJOR,
                        .op     = MATLIB_NO_TRANS,
                        .elem_p = (matlib_real*)q.elem_p};
    matlib_xm phi_tmp = { .lenc   = Q.lenr, 
                          .lenr   = COMPLEX_DIM,
                          .order  = MATLIB_ROW_MAJOR,
                          .op     = MATLIB_NO_TRANS,
                          .elem_p = (matlib_real*)phi.elem_p};

    for ( i = 0; i < na.nr_domains; 
          i++, 
          (phi_tmp.elem_p) += (COMPLEX_DIM * Q.lenr), 
          (q_tmp.elem_p)   += (COMPLEX_DIM * Q.lenc))
    {
        alpha = na.jacob[i];
        error = matlib_xgemm(alpha, Q, phi_tmp, beta, q_tmp );
        err_check( (error == FEM1D_FAILURE), clean_up, 
                   "Failed to compute inner-product (domain index: %d)!", i);
    }
    
    error = matlib_create_zm_sparse( na.nr_domains, na.nr_domains, nnz, M, MATLIB_CSR); 
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the row array!");
    mcnt++; /* 2 */ 

    error = fem1d_periodic_GMMSparsity( na, M->rowIn, M->colIn);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to get the sparsity structure of the mass-matrix!");

    error = fem1d_periodic_ZCSRGMM( na, q, M->rowIn, M->elem_p);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to compute the mass-matrix!");

    
    matlib_free(q.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    if (mcnt == 2)
    {
        matlib_free(M->elem_p);
        matlib_free(M->colIn);
        matlib_free(M->rowIn);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(q.elem_p);
    }

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

}
/*============================================================================*/
fem1d_err fem1d_xm_sparse_GSM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xv         quadW,
    matlib_xv         phi,
    matlib_xm_sparse* M
)
{
    debug_enter( "nr quadrature points: %d, "
                 "length of phi: %d",
                 quadW.len, phi.len);
    
    matlib_index i, mcnt = 0;
    err_check( (quadW.len * na.nr_domains != phi.len), clean_up,
               "Dimension mismatch (phi: %d, required size: %d)!", 
               phi.len, quadW.len * na.nr_domains);

    matlib_index nnz = fem1d_get_nnz(na);
    debug_body( "nr. of non-zero elements: %d", nnz);
    err_check( nnz == 0, clean_up, "%s", 
               "Number of non-zero elements of the mass-matrix is inconsistent!");

    fem1d_err error;
    
    matlib_xv q;
    error = matlib_create_xv( FEM1D_NR_COMBI * na.len, &q, MATLIB_COL_VECT);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the inner product array!");
    mcnt++; /* 1 */ 
    
    phi.len = quadW.len;

    matlib_real *qptr = q.elem_p;
    matlib_real *jptr = na.jacob;
    matlib_real tmp;
    for ( i = 0; i < na.nr_domains; i++,
          (phi.elem_p) += (quadW.len), 
          qptr += FEM1D_NR_COMBI, jptr++)
    {
        tmp = matlib_xdot( quadW, phi);
        err_check( isnan(tmp), clean_up, 
                   "%s", "Failed to compute inner-products!");

        tmp = tmp/(*jptr);
        qptr[FEM1D_INDEX_V11] = tmp * FEM1D_MESI[FEM1D_INDEX_V1][FEM1D_INDEX_V1];
        qptr[FEM1D_INDEX_V12] = tmp * FEM1D_MESI[FEM1D_INDEX_V1][FEM1D_INDEX_V2];
        qptr[FEM1D_INDEX_V22] = tmp * FEM1D_MESI[FEM1D_INDEX_V2][FEM1D_INDEX_V2];
    }
    
    error = matlib_create_xm_sparse( na.len, na.len, nnz, M, MATLIB_CSR); 
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the row array!");
    mcnt++; /* 2 */ 

    error = fem1d_GMMSparsity( na, M->rowIn, M->colIn);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to get the sparsity structure of the mass-matrix!");

    error = fem1d_XCSRGMM( na, q, M->rowIn, M->elem_p);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to compute the mass-matrix!");

    
    matlib_free(q.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    if (mcnt == 2)
    {
        matlib_free(M->elem_p);
        matlib_free(M->colIn);
        matlib_free(M->rowIn);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(q.elem_p);
    }

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

}

fem1d_err fem1d_periodic_xm_sparse_GSM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xv         quadW,
    matlib_xv         phi,
    matlib_xm_sparse* M
)
{
    debug_enter( "nr quadrature points: %d, "
                 "length of phi: %d",
                 quadW.len, phi.len);
    
    matlib_index i, mcnt = 0;
    err_check( (quadW.len * na.nr_domains != phi.len), clean_up,
               "Dimension mismatch (phi: %d, required size: %d)!", 
               phi.len, quadW.len * na.nr_domains);

    matlib_index nnz = fem1d_periodic_get_nnz(na);
    debug_body( "nr. of non-zero elements: %d", nnz);
    err_check( nnz == 0, clean_up, "%s", 
               "Number of non-zero elements of the mass-matrix is inconsistent!");

    fem1d_err error;
    
    matlib_xv q;
    error = matlib_create_xv( FEM1D_NR_COMBI * na.len, &q, MATLIB_COL_VECT);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the inner product array!");
    mcnt++; /* 1 */ 
    
    phi.len = quadW.len;

    matlib_real *qptr = q.elem_p;
    matlib_real *jptr = na.jacob;
    matlib_real tmp;
    for ( i = 0; i < na.nr_domains; i++, 
          (phi.elem_p) += (quadW.len), 
          qptr += FEM1D_NR_COMBI, jptr++)
    {
        tmp = matlib_xdot( quadW, phi);
        err_check( isnan(tmp), clean_up, 
                   "%s", "Failed to compute inner-products!");

        tmp = tmp/(*jptr);
        qptr[FEM1D_INDEX_V11] = tmp * FEM1D_MESI[FEM1D_INDEX_V1][FEM1D_INDEX_V1];
        qptr[FEM1D_INDEX_V12] = tmp * FEM1D_MESI[FEM1D_INDEX_V1][FEM1D_INDEX_V2];
        qptr[FEM1D_INDEX_V22] = tmp * FEM1D_MESI[FEM1D_INDEX_V2][FEM1D_INDEX_V2];
    }
    
    error = matlib_create_xm_sparse( na.nr_domains, na.nr_domains, nnz, M, MATLIB_CSR); 
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the row array!");
    mcnt++; /* 2 */ 

    error = fem1d_periodic_GMMSparsity( na, M->rowIn, M->colIn);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to get the sparsity structure of the mass-matrix!");

    error = fem1d_periodic_XCSRGMM( na, q, M->rowIn, M->elem_p);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to compute the mass-matrix!");

    
    matlib_free(q.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    if (mcnt == 2)
    {
        matlib_free(M->elem_p);
        matlib_free(M->colIn);
        matlib_free(M->rowIn);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(q.elem_p);
    }

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

}

fem1d_err fem1d_zm_sparse_GSM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xv         quadW,
    matlib_zv         phi,
    matlib_zm_sparse* M
)
{
    debug_enter( "nr quadrature points: %d, "
                 "length of phi: %d",
                 quadW.len, phi.len);
    
    matlib_index i, mcnt = 0;
    err_check( (quadW.len * na.nr_domains != phi.len), clean_up,
               "Dimension mismatch (phi: %d, require size: %d)!", 
               phi.len, quadW.len * na.nr_domains);

    matlib_index nnz = fem1d_get_nnz(na);
    debug_body( "nr. of non-zero elements: %d", nnz);
    err_check( nnz == 0, clean_up, "%s", 
               "Number of non-zero elements of the mass-matrix is inconsistent!");

    fem1d_err error;
    
    matlib_zv q;
    error = matlib_create_zv( FEM1D_NR_COMBI * na.nr_domains, &q, MATLIB_COL_VECT);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the inner product array!");
    mcnt++; /* 1 */ 
    
    phi.len = quadW.len;

    matlib_complex *qptr = q.elem_p;
    matlib_real    *jptr = na.jacob;
    matlib_complex tmp;
    for ( i = 0; i < na.nr_domains; i++, 
          (phi.elem_p) += (quadW.len), 
          qptr += FEM1D_NR_COMBI)
    {
        tmp = matlib_xzdot( quadW, phi);
        err_check( isnan(tmp), clean_up, 
                   "%s", "Failed to compute inner-products!");

        tmp = tmp/(*jptr);
        qptr[FEM1D_INDEX_V11] = tmp * FEM1D_MESI[FEM1D_INDEX_V1][FEM1D_INDEX_V1];
        qptr[FEM1D_INDEX_V12] = tmp * FEM1D_MESI[FEM1D_INDEX_V1][FEM1D_INDEX_V2];
        qptr[FEM1D_INDEX_V22] = tmp * FEM1D_MESI[FEM1D_INDEX_V2][FEM1D_INDEX_V2];
    }
    
    error = matlib_create_zm_sparse( na.len, na.len, nnz, M, MATLIB_CSR); 
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the row array!");
    mcnt++; /* 2 */ 

    error = fem1d_GMMSparsity( na, M->rowIn, M->colIn);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to get the sparsity structure of the mass-matrix!");

    error = fem1d_ZCSRGMM( na, q, M->rowIn, M->elem_p);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to compute the mass-matrix!");

    
    matlib_free(q.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    if (mcnt == 2)
    {
        matlib_free(M->elem_p);
        matlib_free(M->colIn);
        matlib_free(M->rowIn);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(q.elem_p);
    }

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

}


fem1d_err fem1d_periodic_zm_sparse_GSM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xv         quadW,
    matlib_zv         phi,
    matlib_zm_sparse* M
)
{
    debug_enter( "nr quadrature points: %d, "
                 "length of phi: %d",
                 quadW.len, phi.len);
    
    matlib_index i, mcnt = 0;
    err_check( (quadW.len * na.nr_domains != phi.len), clean_up,
               "Dimension mismatch (phi: %d, require size: %d)!", 
               phi.len, quadW.len * na.nr_domains);

    matlib_index nnz = fem1d_periodic_get_nnz(na);
    debug_body( "nr. of non-zero elements: %d", nnz);
    err_check( nnz == 0, clean_up, "%s", 
               "Number of non-zero elements of the mass-matrix is inconsistent!");

    fem1d_err error;
    
    matlib_zv q;
    error = matlib_create_zv( FEM1D_NR_COMBI * na.nr_domains, &q, MATLIB_COL_VECT);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the inner product array!");
    mcnt++; /* 1 */ 
    
    phi.len = quadW.len;

    matlib_complex *qptr = q.elem_p;
    matlib_real    *jptr = na.jacob;
    matlib_complex tmp;
    for ( i = 0; i < na.nr_domains; i++, 
          (phi.elem_p) += (quadW.len), 
          qptr += FEM1D_NR_COMBI)
    {
        tmp = matlib_xzdot( quadW, phi);
        err_check( isnan(tmp), clean_up, 
                   "%s", "Failed to compute inner-products!");

        tmp = tmp/(*jptr);
        qptr[FEM1D_INDEX_V11] = tmp * FEM1D_MESI[FEM1D_INDEX_V1][FEM1D_INDEX_V1];
        qptr[FEM1D_INDEX_V12] = tmp * FEM1D_MESI[FEM1D_INDEX_V1][FEM1D_INDEX_V2];
        qptr[FEM1D_INDEX_V22] = tmp * FEM1D_MESI[FEM1D_INDEX_V2][FEM1D_INDEX_V2];
    }
    
    error = matlib_create_zm_sparse( na.nr_domains, na.nr_domains, nnz, M, MATLIB_CSR); 
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the row array!");
    mcnt++; /* 2 */ 

    error = fem1d_periodic_GMMSparsity( na, M->rowIn, M->colIn);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to get the sparsity structure of the mass-matrix!");

    error = fem1d_periodic_ZCSRGMM( na, q, M->rowIn, M->elem_p);
    err_check( (error == FEM1D_FAILURE), clean_up, 
               "%s", "Failed to compute the mass-matrix!");

    
    matlib_free(q.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM1D_SUCCESS;

clean_up:
    if (mcnt == 2)
    {
        matlib_free(M->elem_p);
        matlib_free(M->colIn);
        matlib_free(M->rowIn);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(q.elem_p);
    }

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM1D_FAILURE;

}

/*============================================================================*/


/*============================================================================+/
 | 
 |  Name   : fem2d.c                                                      
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

//#define NDEBUG
//#define MATLIB_NTRACE_DATA

#include "fem2d.h"
#include "assert.h"

/*============================================================================*/
/*
 * Reference triangle (K_t):
 *                                    
 *                 xi_2
 *                 ^  
 *     (-1,1)      |
 *       +         +--> xi_1 
 *       |\
 *       | \
 *       |  \
 *       |   \
 *       +----+
 * (-1,-1)     (1,-1)
 *                                    
 * vertices of the reference triangle
 * */ 

/* Cartesian coordinate indices */ 
const matlib_index FEM2D_INDEX_DIM1 = 0;
const matlib_index FEM2D_INDEX_DIM2 = 1;
const matlib_index FEM2D_DIM = 2; /* cartesian dimensions: R^2 */ 

/* Jacobian indices */ 
const matlib_index FEM2D_INDEX_J11 = 0;
const matlib_index FEM2D_INDEX_J12 = 1;
const matlib_index FEM2D_INDEX_J21 = 2;
const matlib_index FEM2D_INDEX_J22 = 3;

/* Vertex indices */ 
const matlib_index FEM2D_INDEX_V1 = 0; 
const matlib_index FEM2D_INDEX_V2 = 1;   
const matlib_index FEM2D_INDEX_V3 = 2;   
const matlib_index FEM2D_NV = 3;   /* nr of vertices of a triangular element */ 

const matlib_index FEM2D_INIT_VPATCH_SIZE = 10;
/*============================================================================*/

fem2d_err fem2d_create_cc
(
    const matlib_index length,
          fem2d_cc*    nodes
)
{
    debug_enter("length: %d", length);

    assert(length>0);
    nodes->len = length;

    errno  = 0;
    nodes->elem_p = calloc( FEM2D_DIM * length, sizeof(matlib_real));
    err_check( (nodes->elem_p == NULL), clean_up,
               "%s: initialization error: nodes of length %d", 
                strerror(errno), length );

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;
}

/*============================================================================*/

fem2d_err fem2d_refbasis
(
    const fem2d_cc   xi,
          matlib_xm* vphi
)
{

    debug_enter( "length of xi: %d", xi.len);

    err_check( (xi.elem_p==NULL), clean_up, 
               "%s", "Null pointer ecountered!");

    fem2d_err error; 
    error = matlib_create_xm( xi.len, FEM2D_NV, vphi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Memory allocation for basis function-value matrix failed!");


    int i;
    matlib_real *ptr1, *ptr2;
    
    ptr2 = vphi->elem_p;

    /* vphi_1 = -(xi_1 + xi_2)/2 
     * */ 
    for ( ptr1 = xi.elem_p; 
          ptr1 < (xi.elem_p + FEM2D_DIM * xi.len);
          ptr1 += FEM2D_DIM, ptr2++)
    {
        *ptr2 = -0.5*(   *(ptr1 + FEM2D_INDEX_DIM1) 
                       + *(ptr1 + FEM2D_INDEX_DIM2));
    }
    
    /* vphi_2 = (1 + xi_1)/2 
     * */ 
    for ( ptr1 = xi.elem_p; 
          ptr1 < (xi.elem_p + FEM2D_DIM * xi.len);
          ptr1 += FEM2D_DIM, ptr2++)
    {
        *ptr2 =  0.5*( 1.0 + *(ptr1 + FEM2D_INDEX_DIM1));
    }

    /* vphi_3 = (1 + xi_2)/2 
     * */ 
    for ( ptr1 = xi.elem_p; 
          ptr1 < (xi.elem_p + FEM2D_DIM * xi.len);
          ptr1 += FEM2D_DIM, ptr2++)
    {
        *ptr2 =  0.5*( 1.0 + *(ptr1 + FEM2D_INDEX_DIM2));
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;

}

/*============================================================================*/

fem2d_err fem2d_create_ea
(
    const fem2d_cc      nodes,
    const matlib_index  *const ia,
    const matlib_index  nr_domains,
          fem2d_ea*     ea
 
)
/* 
 * ia : nr_domain-by-3 integer array defining the vertices of a triangular elements
 * */ 
{
    debug_enter( "Nr. of nodes: %d, " 
                 "nr. of domains: %d", 
                 nodes.len, nr_domains);

    matlib_index mcnt = 0; /* allocated memory count */ 

    err_check( (nodes.elem_p==NULL) || (ia == NULL), clean_up,
               "%s", "Null pointer ecountered!");

    ea->len      = nr_domains;
    ea->nbase    = nodes.elem_p;
    ea->nr_nodes = nodes.len;

    errno = 0;
    ea->elem_p = calloc( nr_domains, sizeof(fem2d_te));
    err_check( (ea->elem_p == NULL), clean_up, 
               "%s: memory allocation for array of length %d failed!", 
                strerror(errno), nr_domains);
    mcnt++;

    matlib_index i = 0, di = 0;
    fem2d_te* dptr = ea->elem_p;

    matlib_real det = 0;
    matlib_real *vert1, *vert2, *vert3;

    for ( i = 0, di = 0;
          i < FEM2D_NV * nr_domains; 
          i += FEM2D_NV, di++, dptr++)
    {
        dptr->domain_index = di;
        debug_body("domain nr: %d", dptr->domain_index);
        errno = 0;
        dptr->vert_p = calloc( FEM2D_NV, sizeof(matlib_real*));
        err_check( (ea->elem_p == NULL), clean_up, 
                   "%s: memory allocation for vertex pointers failed (domain nr: %d)!", 
                   strerror(errno), dptr->domain_index);

        errno = 0;
        dptr->jmat = calloc( FEM2D_DIM * FEM2D_DIM, sizeof(matlib_real));
        err_check( (dptr->jmat == NULL), clean_up, 
                   "%s: memory allocation for Jacobian matrix failed (domain nr: %d)!", 
                   strerror(errno), dptr->domain_index);

        errno = 0;
        dptr->ijmat = calloc( FEM2D_DIM * FEM2D_DIM, sizeof(matlib_real));
        err_check( (dptr->ijmat == NULL), clean_up, 
                   "%s: memory allocation for inverse Jacobian matrix failed (domain nr: %d)!", 
                   strerror(errno), dptr->domain_index);

        vert1 = nodes.elem_p + FEM2D_DIM * ia[i + FEM2D_INDEX_V1];
        vert2 = nodes.elem_p + FEM2D_DIM * ia[i + FEM2D_INDEX_V2];
        vert3 = nodes.elem_p + FEM2D_DIM * ia[i + FEM2D_INDEX_V3];

        (dptr->vert_p)[FEM2D_INDEX_V1] = vert1;
        (dptr->vert_p)[FEM2D_INDEX_V2] = vert2;
        (dptr->vert_p)[FEM2D_INDEX_V3] = vert3;

        dptr->jmat[FEM2D_INDEX_J11] = 0.5*(vert2[FEM2D_INDEX_DIM1]-vert1[FEM2D_INDEX_DIM1]); 
        dptr->jmat[FEM2D_INDEX_J12] = 0.5*(vert3[FEM2D_INDEX_DIM1]-vert1[FEM2D_INDEX_DIM1]);
        dptr->jmat[FEM2D_INDEX_J21] = 0.5*(vert2[FEM2D_INDEX_DIM2]-vert1[FEM2D_INDEX_DIM2]); 
        dptr->jmat[FEM2D_INDEX_J22] = 0.5*(vert3[FEM2D_INDEX_DIM2]-vert1[FEM2D_INDEX_DIM2]);
        
        det =   dptr->jmat[FEM2D_INDEX_J11] * dptr->jmat[FEM2D_INDEX_J22] 
              - dptr->jmat[FEM2D_INDEX_J12] * dptr->jmat[FEM2D_INDEX_J21];
        
        err_check( (det <= 0), clean_up, 
                   "Zero or negative determinant encountered"
                   "(domain: %d, det = %0.16f)!", 
                    dptr->domain_index, det);
        dptr->jacob = fabs(det);
        /* 
         * M = [ a, b;                   
         *       c, d ];
         * 
         * iM = (1/det(M)) * [  d, -b;
         *                     -c,  a]
         * */
        dptr->ijmat[FEM2D_INDEX_J11] =  dptr->jmat[FEM2D_INDEX_J22]/det; 
        dptr->ijmat[FEM2D_INDEX_J12] = -dptr->jmat[FEM2D_INDEX_J12]/det;
        dptr->ijmat[FEM2D_INDEX_J21] = -dptr->jmat[FEM2D_INDEX_J21]/det; 
        dptr->ijmat[FEM2D_INDEX_J22] =  dptr->jmat[FEM2D_INDEX_J11]/det;
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    if (mcnt>0)
    {
        for ( ; (i >= 0) && (dptr != NULL); 
                i -= FEM2D_NV, dptr--)
        {
            matlib_free(dptr->vert_p);
            matlib_free(dptr->jmat);
            matlib_free(dptr->ijmat);
        }
        matlib_free(ea->elem_p);
    }
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;

}
/*============================================================================*/

fem2d_err fem2d_create_vp
(
    const matlib_index *const ia,
          fem2d_ea     *ea
)
{
    debug_enter("nr of nodes: %d", ea->nr_nodes);
    err_check( (ea->elem_p == NULL) || (ia == NULL), clean_up, 
               "%s", "Null pointers encountered!");

    matlib_index i = 0, mcnt = 0;
    matlib_index init_size = FEM2D_INIT_VPATCH_SIZE; /* maximum six domains expected in a patch */

    errno = 0;
    ea->vpatch_p = calloc( ea->nr_nodes, sizeof(fem2d_vp));
    err_check( (ea->vpatch_p == NULL), clean_up, 
                "%s: memory allocation for vertex patch failed (nr nodes: %d)!", 
                strerror(errno), ea->nr_nodes);
    mcnt++;

    errno = 0;
    fem2d_vp* vp_ptr = NULL; 
    for ( vp_ptr = ea->vpatch_p; 
          vp_ptr < ea->nr_nodes + ea->vpatch_p; 
          vp_ptr++)
    {
            vp_ptr->domain_index = calloc(init_size, sizeof(matlib_index));
            err_check( (vp_ptr->domain_index == NULL), clean_up, 
                        "%s", "Memory allocation failed!", strerror(errno));

            vp_ptr->vert_index   = calloc(init_size, sizeof(matlib_index));
            err_check( (vp_ptr->vert_index == NULL), clean_up, 
                        "%s", "Memory allocation failed!", strerror(errno));

            vp_ptr->len = 0;
    }
    debug_body("%s", "Memory allocated!");
    
    matlib_index tofill, index;
    bool init_size_OK = true;
    for ( i = 0; i< ea->len; i++)
    {
        index  = FEM2D_NV * i + FEM2D_INDEX_V1;
        tofill = (ea->vpatch_p[ia[index]]).len;
        if (tofill <= init_size )
        {
            (ea->vpatch_p[ia[index]]).domain_index[tofill] = i;
            (ea->vpatch_p[ia[index]]).vert_index[tofill]   = FEM2D_INDEX_V1;
            (ea->vpatch_p[ia[index]]).len ++;
        }
        else
        {
            init_size_OK = false;
            break;
        }

        index  = FEM2D_NV * i + FEM2D_INDEX_V2;
        tofill = (ea->vpatch_p[ia[index]]).len;
        if (tofill <= init_size )
        {
            (ea->vpatch_p[ia[index]]).domain_index[tofill] = i;
            (ea->vpatch_p[ia[index]]).vert_index[tofill]   = FEM2D_INDEX_V2;
            (ea->vpatch_p[ia[index]]).len ++;
        }
        else
        {
            init_size_OK = false;
            break;
        }

        index  = FEM2D_NV * i + FEM2D_INDEX_V3;
        tofill = (ea->vpatch_p[ia[index]]).len;
        if (tofill <= init_size )
        {
            (ea->vpatch_p[ia[index]]).domain_index[tofill] = i;
            (ea->vpatch_p[ia[index]]).vert_index[tofill]   = FEM2D_INDEX_V3;
            (ea->vpatch_p[ia[index]]).len ++;
        }
        else
        {
            init_size_OK = false;
            break;
        }
    }
    err_check( !init_size_OK, clean_up, 
                "%s", "Initial size of the memory allocated is insufficient.\n"
                "Returning with an error because this means triangulation is really bad!!!");

    /* Resize the length of vectors */
    for (i = 0; i < ea->nr_nodes; i++)
    {
        ea->vpatch_p[i].domain_index = realloc( ea->vpatch_p[i].domain_index, 
                                                ea->vpatch_p[i].len * sizeof(matlib_index)); 
        err_check( (ea->vpatch_p[i].domain_index == NULL), clean_up, 
                        "%s", "Memory reallocation failed!", strerror(errno));

        ea->vpatch_p[i].vert_index = realloc( ea->vpatch_p[i].vert_index, 
                                              ea->vpatch_p[i].len * sizeof(matlib_index)); 
        err_check( (ea->vpatch_p[i].vert_index == NULL), clean_up, 
                        "%s", "Memory reallocation failed!", strerror(errno));
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    if (mcnt == 1)
    {
        for ( ; (vp_ptr < ea->vpatch_p) & (vp_ptr != NULL); 
                vp_ptr--)
        {
            matlib_free(vp_ptr->domain_index);
            matlib_free(vp_ptr->vert_index);
        }
        matlib_free(ea->vpatch_p);
    }
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;
}

void fem2d_free_ea(fem2d_ea ea)
{
    fem2d_te* dptr;

    for ( dptr = ea.elem_p; 
          (dptr < (ea.elem_p + ea.len)) && (dptr != NULL);
          dptr++)
    {
        matlib_free(dptr->vert_p);
        matlib_free(dptr->jmat);
        matlib_free(dptr->ijmat);
    }
    fem2d_vp* vp_ptr;
    for ( vp_ptr = ea.vpatch_p; 
          (vp_ptr < ea.vpatch_p + ea.nr_nodes) & (vp_ptr != NULL); 
          vp_ptr++)
    {
        matlib_free(vp_ptr->domain_index);
        matlib_free(vp_ptr->vert_index);
    }
    matlib_free(ea.vpatch_p);

    matlib_free(ea.elem_p);
}

/*============================================================================*/
fem2d_err fem2d_create_ia
(
    const fem2d_ea      ea,
          matlib_index* ia
)
{

    debug_enter("nr of domains: %d", ea.len);
    err_check( (ea.elem_p==NULL), clean_up,
               "%s", "Null pointer ecountered!");

    fem2d_te* dptr = ea.elem_p; /* domain pointer */
    matlib_real* nptr0 = ea.nbase;

    matlib_index* ptr;
    
    matlib_index offset; 
    matlib_index vsize = (matlib_index)2*sizeof(matlib_real*);

    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); 
          dptr++, ia+=3)
    {

        /* Get the first vertex */ 
        offset = (matlib_index)(((matlib_index)*(dptr->vert_p) 
                               - (matlib_index)nptr0)/vsize);
        *(ia+0) = offset;


        /* Get the second vertex */
        offset = (matlib_index)(((matlib_index)*(dptr->vert_p+1) 
                               - (matlib_index)nptr0)/vsize);
        *(ia+1) = offset;

        /* Get the third vertex */
        offset = (matlib_index)(((matlib_index)*(dptr->vert_p+2) 
                               - (matlib_index)nptr0)/vsize);
        *(ia+2) = offset;

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;

    }
}
/*============================================================================*/

fem2d_err fem2d_centroid
(
    const fem2d_ea ea,
          fem2d_cc cen
)
{
    debug_enter( "length of cen: %d, "
                 "nr. of domains: %d",
                 cen.len, ea.len);
    err_check( (cen.len != ea.len), clean_up, 
               "Dimension mismatch (cen.len = %d, ea.len = %d)",
               cen.len, ea.len );
    err_check( (ea.elem_p==NULL) || (cen.elem_p == NULL), clean_up,
               "%s", "Null pointers ecountered!");

    fem2d_te* dptr; /* domain pointer */ 
    
    matlib_real *vptr, *cptr = cen.elem_p; 

    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); 
          dptr++, cptr += FEM2D_DIM )
    {
        cptr[FEM2D_INDEX_DIM1] = 0; 
        cptr[FEM2D_INDEX_DIM2] = 0; 

        /* Get the first vertex */ 
        vptr = *(dptr->vert_p);
        cptr[FEM2D_INDEX_DIM1] += vptr[FEM2D_INDEX_DIM1]/3.0;  
        cptr[FEM2D_INDEX_DIM2] += vptr[FEM2D_INDEX_DIM2]/3.0; 

        /* Get the second vertex */
        vptr = *(dptr->vert_p+1);
        cptr[FEM2D_INDEX_DIM1] += vptr[FEM2D_INDEX_DIM1]/3.0;  
        cptr[FEM2D_INDEX_DIM2] += vptr[FEM2D_INDEX_DIM2]/3.0; 

        /* Get the third vertex */
        vptr = *(dptr->vert_p+2);
        cptr[FEM2D_INDEX_DIM1] += vptr[FEM2D_INDEX_DIM1]/3.0;  
        cptr[FEM2D_INDEX_DIM2] += vptr[FEM2D_INDEX_DIM2]/3.0; 
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;
}
/*============================================================================*/

fem2d_err fem2d_ref2mesh
(
    const fem2d_ea  ea,
    const matlib_xm vphi,
          fem2d_cc* x
)
{
    debug_enter( "nr domains: %d, "
                 "size of vphi: %d-by-%d", 
                 ea.len, vphi.lenc, vphi.lenr);

    err_check( ((ea.elem_p==NULL) || (vphi.elem_p == NULL)) || (x==NULL), 
               clean_up, "%s", "Null pointers ecountered!");

    err_check( (vphi.lenr != FEM2D_NV), clean_up, 
               "Size of second input incorrect (input: %d-by-%d)!",
               vphi.lenc, vphi.lenr);

    DEBUG_PRINT_XM(vphi, "%s: ", "Ref basis");

    fem2d_err error;
    error = fem2d_create_cc( vphi.lenc * ea.len, x);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "Initialization of coordinate "
               "array (of length: %d) failed!", 
               vphi.lenc * ea.len);

    matlib_real alpha = 1.0;
    matlib_real beta  = 0.0;

    fem2d_te* dptr; /* domain pointer */
    matlib_real* vptr; /* vertex pointer */ 

    matlib_xm nodes_mat;
    matlib_create_xm( FEM2D_NV, FEM2D_DIM, &nodes_mat, 
                      MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    
    matlib_real* ptr; 
 
    matlib_xm x_tmp = { .lenc   = vphi.lenc, 
                        .lenr   = FEM2D_DIM, 
                        .order  = MATLIB_ROW_MAJOR,
                        .op     = MATLIB_NO_TRANS,
                        .elem_p = x->elem_p };

    
    /* 
     * vphi: size xi.len - by - 3 
     * Convert from COL_MAJOR to ROW_MAJOR 
     * */ 

    matlib_xm vphi_tmp = { .order = MATLIB_ROW_MAJOR,
                           .op    = MATLIB_TRANS,
                           .lenr  = vphi.lenc,
                           .lenc  = FEM2D_NV,
                           .elem_p = vphi.elem_p };

    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); dptr++)
    {
        ptr = nodes_mat.elem_p;

        /* Get the first vertex */ 
        vptr = *(dptr->vert_p);
        *ptr = vptr[FEM2D_INDEX_DIM1]; ptr++;
        *ptr = vptr[FEM2D_INDEX_DIM2]; ptr++;

        /* Get the second vertex */
        vptr = *(dptr->vert_p+1);
        *ptr = vptr[FEM2D_INDEX_DIM1]; ptr++;
        *ptr = vptr[FEM2D_INDEX_DIM2]; ptr++;

        /* Get the third vertex */
        vptr = *(dptr->vert_p+2);
        *ptr = vptr[FEM2D_INDEX_DIM1]; ptr++;
        *ptr = vptr[FEM2D_INDEX_DIM2];
        
        matlib_xgemm( alpha, vphi_tmp, nodes_mat, beta, x_tmp );
        
        DEBUG_PRINT_XM(x_tmp, "Points in domain[%d]: ", dptr->domain_index);

        x_tmp.elem_p = x_tmp.elem_p + FEM2D_DIM * vphi.lenc;
    }
    
    matlib_free(nodes_mat.elem_p);
    debug_exit( "Total nr. of sampling points: %d, Exit Status: %s",
                x->len, "SUCCESS" );
    return FEM2D_SUCCESS;
    
clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;
}
/*============================================================================*/

fem2d_err fem2d_interp
(
    const fem2d_ea  ea,
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
    
    err_check( ((ea.elem_p == NULL)   || (u_nodes.elem_p == NULL)) 
               || (vphi.elem_p == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    /* Length of u_interp = nr of quad. nodes x nr of domains 
     * */ 
    err_check( (u_interp.len != (vphi.lenc * ea.len) ), clean_up, 
               "Dimension mismatch (u_interp: %d)", u_interp.len);
    
    matlib_real alpha = 1.0;
    matlib_real beta  = 0.0;

    matlib_index INDEX_REAL  = 0;
    matlib_index INDEX_IMAG  = 1;
    matlib_index COMPLEX_DIM = 2; /* Dimension of a compelex plane */ 

    matlib_xm u1_tmp;    
    matlib_create_xm( FEM2D_NV, COMPLEX_DIM, &u1_tmp, 
                      MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);

    matlib_xm u2_tmp = { .lenc   = vphi.lenc, 
                         .lenr   = COMPLEX_DIM, 
                         .order  = MATLIB_ROW_MAJOR,
                         .op     = MATLIB_NO_TRANS,
                         .elem_p = (matlib_real*)u_interp.elem_p };

    fem2d_te* dptr; /* domain pointer */
    matlib_real* nptr0 = ea.nbase; /* base address */ 

    matlib_real* ptr; 
    matlib_real* uptr; /* vertex pointer */ 
    matlib_real* uptr0 = (matlib_real*)u_nodes.elem_p; /* base address */ 

    matlib_xm vphi_tmp = { .order = MATLIB_ROW_MAJOR,
                           .op    = MATLIB_TRANS,
                           .lenr  = vphi.lenc,
                           .lenc  = FEM2D_NV,
                           .elem_p = vphi.elem_p };


    matlib_index offset;
    matlib_index vsize = (matlib_index)sizeof(matlib_real*);

    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); dptr++)
    {
        ptr = u1_tmp.elem_p;

        /* Get the first vertex */ 
        offset = (matlib_index)(((matlib_index)*(dptr->vert_p) 
                               - (matlib_index)nptr0)/vsize);
        uptr = uptr0 + offset; 
        *ptr = uptr[INDEX_REAL]; ptr++;
        *ptr = uptr[INDEX_IMAG]; ptr++;
        debug_body( "domain: %d, u: %0.16f%+0.16fi, offset: %d", 
                    dptr->domain_index, uptr[INDEX_REAL], uptr[INDEX_IMAG], offset );

        /* Get the second vertex */
        offset = (matlib_index)(((matlib_index)*(dptr->vert_p+1) 
                               - (matlib_index)nptr0)/vsize);
        uptr = uptr0 + offset; 
        *ptr = uptr[INDEX_REAL]; ptr++;
        *ptr = uptr[INDEX_IMAG]; ptr++;
        debug_body( "domain: %d, u: %0.16f%+0.16fi, offset: %d", 
                    dptr->domain_index, uptr[INDEX_REAL], uptr[INDEX_IMAG], offset);

        /* Get the third vertex */
        offset = (matlib_index)(((matlib_index)*(dptr->vert_p+2) 
                               - (matlib_index)nptr0)/vsize);
        uptr = uptr0 + offset;
        *ptr = uptr[INDEX_REAL]; ptr++;
        *ptr = uptr[INDEX_IMAG];
        debug_body( "domain: %d, u: %0.16f%+0.16fi, offset: %d", 
                    dptr->domain_index, uptr[INDEX_REAL], uptr[INDEX_IMAG], offset );
        
        matlib_xgemm( alpha, vphi_tmp, u1_tmp, beta, u2_tmp );
  
        u2_tmp.elem_p = u2_tmp.elem_p + COMPLEX_DIM * vphi.lenc; 

    }
    matlib_free(u1_tmp.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;
    
clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;
}

/*============================================================================*/
matlib_real fem2d_normL2
(
    fem2d_ea  ea,
    matlib_zv u_nodes,
    matlib_xm vphi,
    matlib_xv quadW
)
{
    debug_enter( "nr of domains: %d, "
                 "nr of nodes: %d, "
                 "vphi: %d-by-%d, quadW: %d",
                 ea.len, u_nodes.len,
                 vphi.lenc, vphi.lenr, quadW.len);
    
    err_check( ((ea.elem_p == NULL)    || (u_nodes.elem_p == NULL)) 
            || ((vphi.elem_p == NULL)  || (quadW.elem_p == NULL)), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check( ((vphi.lenc != quadW.len) || (vphi.lenr != FEM2D_NV)), clean_up, 
               "Dimension mis-match (vphi: %d-by-%d, quadW: %d)!",
               vphi.lenc, vphi.lenr, quadW.len);

    matlib_index mcnt = 0;
    fem2d_err error;

    matlib_zv u_interp;
    error = matlib_create_zv( ea.len * vphi.lenc, &u_interp, MATLIB_COL_VECT);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "Memory allocation for complex vector failed (length: %d)!", 
               ea.len * vphi.lenc);
    mcnt++;

    error = fem2d_interp(ea, u_nodes, vphi, u_interp);
    err_check((error == FEM2D_FAILURE), clean_up, "%s", "Interpolation failed!");

    matlib_complex* ptr = u_interp.elem_p;
    
    matlib_real norm = 0;
    matlib_real jacob;

    for ( matlib_index i=0; i< ea.len; i++)
    {
        jacob = (ea.elem_p[i]).jacob;
        for (matlib_index j=0; j< quadW.len; j++, ptr++)
        {
            norm += (*ptr * conj(*ptr) * quadW.elem_p[j] * jacob);
        }
    }
    debug_exit("Exit Status: %s", "SUCCESS" );
    return sqrt(norm);

clean_up:

    if (mcnt>0)
        matlib_free(u_interp.elem_p);

    debug_exit("Exit Status: %s", "FAILURE" );
    return FP_NAN;
}
/*============================================================================*/

matlib_real fem2d_iprod 
/* Inner product iprod = (u, conj(v))_{\Omega} */ 
(
    fem2d_ea  ea,
    matlib_zv u_nodes,
    matlib_zv v_nodes,
    matlib_xm vphi,
    matlib_xv quadW
)
{
    debug_enter( "nr of domains: %d, "
                 "nr of nodes: %d, "
                 "vphi: %d-by-%d, quadW: %d",
                 ea.len, u_nodes.len,
                 vphi.lenc, vphi.lenr, quadW.len);
    
    err_check(    (ea.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL) 
               || (v_nodes.elem_p == NULL) 
               || (vphi.elem_p    == NULL)  
               || (quadW.elem_p   == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check( ((vphi.lenc != quadW.len) || (vphi.lenr != FEM2D_NV)), clean_up, 
               "Dimension mis-match (vphi: %d-by-%d, quadW: %d)!",
               vphi.lenc, vphi.lenr, quadW.len);

    err_check( (u_nodes.len != v_nodes.len), clean_up, 
               "Dimension mis-match for field vectors (u: %d, v: %d)!",
               v_nodes.len, u_nodes.len);

    matlib_index mcnt = 0;
    fem2d_err error;

    matlib_zv u_interp, v_interp;
    error = matlib_create_zv( ea.len * vphi.lenc, &u_interp, MATLIB_COL_VECT);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "Memory allocation for complex vector failed (length: %d)!", 
               ea.len * vphi.lenc);
    mcnt++; /* 1 */ 
    error = matlib_create_zv( ea.len * vphi.lenc, &v_interp, MATLIB_COL_VECT);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "Memory allocation for complex vector failed (length: %d)!", 
               ea.len * vphi.lenc);
    mcnt++; /* 2 */ 

    error = fem2d_interp(ea, u_nodes, vphi, u_interp);
    err_check((error == FEM2D_FAILURE), clean_up, "%s", "Interpolation failed!");

    error = fem2d_interp(ea, v_nodes, vphi, v_interp);
    err_check((error == FEM2D_FAILURE), clean_up, "%s", "Interpolation failed!");

    matlib_complex* uptr = u_interp.elem_p;
    matlib_complex* vptr = v_interp.elem_p;
    
    matlib_real iprod = 0;
    matlib_real jacob;

    for ( matlib_index i=0; i< ea.len; i++)
    {
        jacob = (ea.elem_p[i]).jacob;
        for (matlib_index j=0; j< quadW.len; j++, uptr++, vptr++)
        {
            iprod += (*uptr * conj(*vptr) * quadW.elem_p[j] * jacob);
        }
    }
    debug_exit("Exit Status: %s", "SUCCESS" );
    return iprod;

clean_up:

    if (mcnt==2)
    {
        matlib_free(u_interp.elem_p);
        mcnt--;
    }
    if (mcnt==1)
        matlib_free(u_interp.elem_p);

    debug_exit("Exit Status: %s", "FAILURE" );
    return FP_NAN;
} /* fem2d_iprod */ 

/*============================================================================*/

fem2d_err fem2d_quadP
/* Qudature matrix for computing projection on vertex functions */ 
(
    fem2d_cc   xi,
    matlib_xm* vphi,
    matlib_xv  quadW
)
{

    debug_enter( "length of xi: %d", xi.len);

    err_check( (xi.elem_p==NULL), clean_up, 
               "%s", "Null pointer ecountered!");

    fem2d_err error; 
    error = matlib_create_xm( xi.len, FEM2D_NV, vphi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Memory allocation for basis function-value matrix failed!");


    int i;
    matlib_real *ptr1, *ptr2;
    
    ptr2 = vphi->elem_p;

    /* vphi_1 = -(xi_1 + xi_2)/2 
     * */ 
    for ( ptr1 = xi.elem_p; 
          ptr1 < (xi.elem_p + FEM2D_DIM * xi.len);
          ptr1 += FEM2D_DIM, ptr2++)
    {
        *ptr2 = -0.5*(   *(ptr1 + FEM2D_INDEX_DIM1) 
                       + *(ptr1 + FEM2D_INDEX_DIM2));
    }
    
    /* vphi_2 = (1 + xi_1)/2 
     * */ 
    for ( ptr1 = xi.elem_p; 
          ptr1 < (xi.elem_p + FEM2D_DIM * xi.len);
          ptr1 += FEM2D_DIM, ptr2++)
    {
        *ptr2 =  0.5*( 1.0 + *(ptr1 + FEM2D_INDEX_DIM1));
    }

    /* vphi_3 = (1 + xi_2)/2 
     * */ 
    for ( ptr1 = xi.elem_p; 
          ptr1 < (xi.elem_p + FEM2D_DIM * xi.len);
          ptr1 += FEM2D_DIM, ptr2++)
    {
        *ptr2 =  0.5*( 1.0 + *(ptr1 + FEM2D_INDEX_DIM2));
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;

}



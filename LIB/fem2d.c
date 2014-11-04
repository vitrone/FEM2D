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

#define NDEBUG
#define MATLIB_NTRACE_DATA

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

/* Ref. triangle vertices */
const matlib_real FEM2D_VERT[6] = { -1.0, -1.0,
                                     1.0, -1.0,
                                    -1.0,  1.0};



const matlib_index FEM2D_INIT_VPATCH_SIZE = 10;

const matlib_real MEMI[3][3] = { {1.0/3.0,  1.0/6.0,  1.0/6.0}, 
                                 {1.0/6.0,  1.0/3.0,  1.0/6.0},
                                 {1.0/6.0,  1.0/6.0,  1.0/3.0}};

/* 
 * Combination of basis functions for mass-matrix
 * (vphi1, vphi1) --> FEM2D_INDEX_V11
 * (vphi1, vphi2) --> FEM2D_INDEX_V12
 * (vphi1, vphi3) --> FEM2D_INDEX_V13
 * (vphi2, vphi2) --> FEM2D_INDEX_V22
 * (vphi2, vphi3) --> FEM2D_INDEX_V23
 * (vphi3, vphi3) --> FEM2D_INDEX_V33
 *
 * */ 
const matlib_index FEM2D_INDEX_V11 = 0;
const matlib_index FEM2D_INDEX_V12 = 1;
const matlib_index FEM2D_INDEX_V13 = 2;
const matlib_index FEM2D_INDEX_V22 = 3;
const matlib_index FEM2D_INDEX_V23 = 4;
const matlib_index FEM2D_INDEX_V33 = 5;
const matlib_index FEM2D_NR_COMBI  = 6;

/*============================================================================*/

fem2d_err fem2d_create_cc
(
    const matlib_index length,
          fem2d_cc*    nodes
)
{
    debug_enter("length: %d", length);
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

    err_check( (xi.elem_p == NULL), clean_up, 
               "%s", "Null pointer ecountered!");

    fem2d_err error; 
    error = matlib_create_xm( xi.len, FEM2D_NV, vphi, MATLIB_COL_MAJOR, MATLIB_NO_TRANS);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Memory allocation for basis function-value matrix failed!");


    int i;
    matlib_real *ptr1, *ptr2;
    
    /* vphi_1 = -(xi_1 + xi_2)/2 
     * */ 
    ptr2 = vphi->elem_p + FEM2D_INDEX_V1 * xi.len;

    for ( ptr1 = xi.elem_p; 
          ptr1 < (xi.elem_p + FEM2D_DIM * xi.len);
          ptr1 += FEM2D_DIM, ptr2++)
    {
        *ptr2 = -0.5*(   *(ptr1 + FEM2D_INDEX_DIM1) 
                       + *(ptr1 + FEM2D_INDEX_DIM2));
    }
    
    /* vphi_2 = (1 + xi_1)/2 
     * */ 
    ptr2 = vphi->elem_p + FEM2D_INDEX_V2 * xi.len;

    for ( ptr1 = xi.elem_p; 
          ptr1 < (xi.elem_p + FEM2D_DIM * xi.len);
          ptr1 += FEM2D_DIM, ptr2++)
    {
        *ptr2 =  0.5*( 1.0 + *(ptr1 + FEM2D_INDEX_DIM1));
    }

    /* vphi_3 = (1 + xi_2)/2 
     * */ 
    ptr2 = vphi->elem_p + FEM2D_INDEX_V3 * xi.len;

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
static void calc_jacobian
(
    matlib_real* vert1,
    matlib_real* vert2,
    matlib_real* vert3,
    fem2d_te* dptr
)
{
        dptr->jmat[FEM2D_INDEX_J11] = 0.5*(vert2[FEM2D_INDEX_DIM1]-vert1[FEM2D_INDEX_DIM1]); 
        dptr->jmat[FEM2D_INDEX_J12] = 0.5*(vert3[FEM2D_INDEX_DIM1]-vert1[FEM2D_INDEX_DIM1]);
        dptr->jmat[FEM2D_INDEX_J21] = 0.5*(vert2[FEM2D_INDEX_DIM2]-vert1[FEM2D_INDEX_DIM2]); 
        dptr->jmat[FEM2D_INDEX_J22] = 0.5*(vert3[FEM2D_INDEX_DIM2]-vert1[FEM2D_INDEX_DIM2]);
        
        matlib_real det = (   dptr->jmat[FEM2D_INDEX_J11] * dptr->jmat[FEM2D_INDEX_J22] 
                            - dptr->jmat[FEM2D_INDEX_J12] * dptr->jmat[FEM2D_INDEX_J21]);
        dptr->jacob = fabs(det);
        dptr->sign  = (det > 0)? 1: -1;
}



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
    ea->node_p   = nodes.elem_p;
    ea->nr_nodes = nodes.len;
    ea->vpatch_p = NULL;

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
    matlib_index index1, index2, index3, pos1, pos2, pos3;

    for ( i = 0, di = 0;
          i < FEM2D_NV * nr_domains; 
          i += FEM2D_NV, di++, dptr++)
    {
        dptr->vert_p = NULL;
        dptr->jmat   = NULL;
        dptr->ijmat  = NULL;
        dptr->domain_index = di;
        errno = 0;
        dptr->vert_p = calloc( FEM2D_NV, sizeof(matlib_real*));
        err_check( (dptr->vert_p == NULL), clean_up, 
                   "%s: memory allocation for vertex pointers failed (domain nr: %d)!", 
                   strerror(errno), dptr->domain_index);
        errno = 0;
        dptr->pos_vpatch = calloc( FEM2D_NV, sizeof(matlib_real));
        err_check( (dptr->pos_vpatch == NULL), clean_up, 
                   "%s: memory allocation failed (domain nr: %d)!", 
                   strerror(errno), dptr->domain_index);

        errno = 0;
        dptr->nindex_p = calloc( FEM2D_NV, sizeof(matlib_index));
        err_check( (dptr->nindex_p == NULL), clean_up, 
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

        /* Vertex indices are sorted in the increasing order
         * */ 

        index1 = ia[i + FEM2D_INDEX_V1];
        index2 = ia[i + FEM2D_INDEX_V2];
        index3 = ia[i + FEM2D_INDEX_V3];

#if 0
        index1 = ia[i + FEM2D_INDEX_V1];
        pos1   = FEM2D_INDEX_V1;
        if (index1 > ia[i + FEM2D_INDEX_V2])
        {
            index1 = ia[i + FEM2D_INDEX_V2];
            pos1 = FEM2D_INDEX_V2;
        }
        if (index1 > ia[i + FEM2D_INDEX_V3])
        {
            index1 = ia[i + FEM2D_INDEX_V3];
            pos1 = FEM2D_INDEX_V3;
        }
        pos2   = ((pos1 + 1) % FEM2D_NV);
        pos3   = ((pos1 + 2) % FEM2D_NV);
        index2 = ia[i + pos2];
        if (index2 > ia[i + pos3])
        {
            index2 = ia[i + pos3];
            index3 = ia[i + pos2];
        }
        else
        {
            index3 = ia[i + pos3];
        }
        debug_body("indices: %d, %d, %d", index1, index2, index3);
        debug_body("pos: %d, %d, %d", pos1, pos2, pos3);

        dptr->jmat[FEM2D_INDEX_J11] = 0.5*(vert2[FEM2D_INDEX_DIM1]-vert1[FEM2D_INDEX_DIM1]); 
        dptr->jmat[FEM2D_INDEX_J12] = 0.5*(vert3[FEM2D_INDEX_DIM1]-vert1[FEM2D_INDEX_DIM1]);
        dptr->jmat[FEM2D_INDEX_J21] = 0.5*(vert2[FEM2D_INDEX_DIM2]-vert1[FEM2D_INDEX_DIM2]); 
        dptr->jmat[FEM2D_INDEX_J22] = 0.5*(vert3[FEM2D_INDEX_DIM2]-vert1[FEM2D_INDEX_DIM2]);
        
        det = (   dptr->jmat[FEM2D_INDEX_J11] * dptr->jmat[FEM2D_INDEX_J22] 
                    - dptr->jmat[FEM2D_INDEX_J12] * dptr->jmat[FEM2D_INDEX_J21]);
        
        dptr->jacob = fabs(det);
        dptr->sign  = (det > 0)? 1: -1;
#endif
        vert1 = nodes.elem_p + FEM2D_DIM * index1;
        vert2 = nodes.elem_p + FEM2D_DIM * index2;
        vert3 = nodes.elem_p + FEM2D_DIM * index3;
        
        calc_jacobian(vert1, vert2, vert3, dptr);
        err_check( (dptr->jacob < MATLIB_TOL), clean_up, 
                   "Zero determinant encountered"
                   "(domain: %d, det = %0.16g)!", 
                    dptr->domain_index, dptr->jacob);
        if (dptr->sign < 0)
        {
            calc_jacobian(vert2, vert1, vert3, dptr);

            (dptr->vert_p)[FEM2D_INDEX_V1] = vert2;
            (dptr->vert_p)[FEM2D_INDEX_V2] = vert1;
            (dptr->vert_p)[FEM2D_INDEX_V3] = vert3;

            (dptr->nindex_p)[FEM2D_INDEX_V1] = index2;
            (dptr->nindex_p)[FEM2D_INDEX_V2] = index1;
            (dptr->nindex_p)[FEM2D_INDEX_V3] = index3;
        }
        else
        {
            (dptr->vert_p)[FEM2D_INDEX_V1] = vert1;
            (dptr->vert_p)[FEM2D_INDEX_V2] = vert2;
            (dptr->vert_p)[FEM2D_INDEX_V3] = vert3;

            (dptr->nindex_p)[FEM2D_INDEX_V1] = index1;
            (dptr->nindex_p)[FEM2D_INDEX_V2] = index2;
            (dptr->nindex_p)[FEM2D_INDEX_V3] = index3;
        }

        /* 
         * M = [ a, b;                   
         *       c, d ];
         * 
         * iM = (1/det(M)) * [  d, -b;
         *                     -c,  a]
         * */
        dptr->ijmat[FEM2D_INDEX_J11] =  dptr->jmat[FEM2D_INDEX_J22]/dptr->jacob; 
        dptr->ijmat[FEM2D_INDEX_J12] = -dptr->jmat[FEM2D_INDEX_J12]/dptr->jacob;
        dptr->ijmat[FEM2D_INDEX_J21] = -dptr->jmat[FEM2D_INDEX_J21]/dptr->jacob; 
        dptr->ijmat[FEM2D_INDEX_J22] =  dptr->jmat[FEM2D_INDEX_J11]/dptr->jacob;
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    if (mcnt>0)
    {
        matlib_int ii = (matlib_int) i;
        for ( ; (ii >= 0) && (dptr != NULL); 
                ii -= FEM2D_NV, dptr--)
        {
            debug_body("Attemptig to free all the allocated memory (i: %d)!", ii);
            matlib_free(dptr->vert_p);
            matlib_free(dptr->pos_vpatch);
            matlib_free(dptr->nindex_p);
            matlib_free(dptr->jmat);
            matlib_free(dptr->ijmat);
        }
        matlib_free(ea->elem_p);
        ea->elem_p = NULL;
    }
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;

}
/*============================================================================*/

fem2d_err fem2d_create_vp(fem2d_ea *ea)
{
    debug_enter("nr of nodes: %d", ea->nr_nodes);
    err_check( (ea->elem_p == NULL), clean_up, 
               "%s", "Null pointers encountered!");

    matlib_index i = 0, mcnt = 0;
    matlib_index init_size = FEM2D_INIT_VPATCH_SIZE; /* maximum six domains expected in a patch */

    errno = 0;
    ea->vpatch_p = calloc( ea->nr_nodes, sizeof(fem2d_vp));
    err_check( (ea->vpatch_p == NULL), clean_up, 
                "%s: memory allocation for vertex patch failed (nr nodes: %d)!", 
                strerror(errno), ea->nr_nodes);
    mcnt++;

    /* Initialize memory */ 
    fem2d_vp* vp_ptr = NULL; 
    for ( vp_ptr = ea->vpatch_p; 
          vp_ptr < ea->nr_nodes + ea->vpatch_p; 
          vp_ptr++)
    {

        errno = 0;
        vp_ptr->domain_p = calloc(init_size, sizeof(fem2d_te*));
        err_check( (vp_ptr->domain_p == NULL), clean_up, 
                    "%s: Memory allocation failed!", strerror(errno));
        errno = 0;
        vp_ptr->vert_index   = calloc(init_size, sizeof(matlib_index));
        err_check( (vp_ptr->vert_index == NULL), clean_up, 
                    "%s: Memory allocation failed!", strerror(errno));

        vp_ptr->len = 0;
    }
    debug_body("%s", "Memory allocated!");
    
    /* Create the vertex patch */ 
    matlib_index offset; 
    matlib_index tofill, index;
    bool init_size_OK = true;

    fem2d_vp* vp_ptr0  = ea->vpatch_p;
    fem2d_vp* vp_ptr1;

    fem2d_te* dptr = ea->elem_p;
    for ( i = 0; i< ea->len; i++, dptr++)
    {
        vp_ptr1 = vp_ptr0 + dptr->nindex_p[FEM2D_INDEX_V1];

        tofill  = vp_ptr1->len;
        if (tofill <= init_size )
        {
            vp_ptr1->domain_p[tofill] = dptr;
            vp_ptr1->vert_index[tofill] = FEM2D_INDEX_V1;
            vp_ptr1->len ++;
        }
        else
        {
            init_size_OK = false;
            break;
        }

        vp_ptr1 = vp_ptr0 + dptr->nindex_p[FEM2D_INDEX_V2];
        tofill  = vp_ptr1->len;
        if (tofill <= init_size )
        {
            vp_ptr1->domain_p[tofill]   = dptr;
            vp_ptr1->vert_index[tofill] = FEM2D_INDEX_V2;
            vp_ptr1->len ++;
        }
        else
        {
            init_size_OK = false;
            break;
        }

        vp_ptr1 = vp_ptr0 + dptr->nindex_p[FEM2D_INDEX_V3];
        tofill  = vp_ptr1->len;
        if (tofill <= init_size )
        {
            vp_ptr1->domain_p[tofill]   = dptr;
            vp_ptr1->vert_index[tofill] = FEM2D_INDEX_V3;
            vp_ptr1->len ++;
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
    for ( i = 0; 
          (i < ea->nr_nodes) && (ea->vpatch_p[i].len != init_size);
          i++)
    {
        errno = 0;
        ea->vpatch_p[i].domain_p = realloc( ea->vpatch_p[i].domain_p, 
                                            ea->vpatch_p[i].len * sizeof(fem2d_te*)); 
        err_check( (ea->vpatch_p[i].domain_p == NULL), clean_up, 
                    "%s: Memory reallocation failed!", strerror(errno));

        errno = 0;
        ea->vpatch_p[i].vert_index = realloc( ea->vpatch_p[i].vert_index, 
                                              ea->vpatch_p[i].len * sizeof(matlib_index)); 
        err_check( (ea->vpatch_p[i].vert_index == NULL), clean_up, 
                   "%s: Memory reallocation failed!", strerror(errno));
    }
    
    mcnt++;

    /* Create edge information and put the domains in right order 
     * such that the neighboring domains are next to each other in the list.
     * */

    /* Don't reuse the following pointer. 
     * It will be used for freeing the memory in case of error!
     * */ 
    fem2d_vp *vp_bptr; 
    
    for ( vp_bptr = ea->vpatch_p; 
          (vp_bptr < ea->nr_nodes + ea->vpatch_p); vp_bptr++)
    {
        errno = 0;
        vp_bptr->bvert_index = calloc(vp_bptr->len, sizeof(matlib_index));
        err_check( (vp_bptr->bvert_index == NULL), clean_up, 
                   "%s: Memory reallocation failed!", strerror(errno));
    }

    matlib_index j, k, l, m, n, tmp;
    matlib_index next_vertex, other_vertex1, other_vertex2;
    matlib_index list1_len, list2_len;
    matlib_index *vert_index_tmp, *bvert_index_tmp;
    fem2d_te **dlist1, **dlist2, **domain_p_tmp;

    fem2d_te *dptr_tmp, *del1;

    bool found;
    

    for (i = 0; i < ea->nr_nodes; i++)
    {
        ea->vpatch_p[i].node_index = i;

        other_vertex1 = (ea->vpatch_p[i].vert_index[0] + 1) % FEM2D_NV;
        other_vertex2 = (ea->vpatch_p[i].vert_index[0] + 2) % FEM2D_NV;
        
        ea->vpatch_p[i].bvert_index[0] = other_vertex1;
        next_vertex = other_vertex2; 
        debug_body( "vert_index: %d, other vertex1: %d, "
                    "other vertex2: %d, next vertex: %d", 
                    ea->vpatch_p[i].vert_index[0], other_vertex1, 
                    other_vertex2, next_vertex);
        
        debug_body("patch length: %d", ea->vpatch_p[i].len);
        for (m = 1; m < ea->vpatch_p[i].len; m++ )
        {
            dptr    = ea->vpatch_p[i].domain_p[m-1];
            vp_ptr1 = vp_ptr0 + dptr->nindex_p[next_vertex];
        
            /* Find the common domain between: 
             *   ea->vpatch_p[i].domain_index (length: ea->vpatch_p[i].len)
             *   vp_ptr->domain_index (length: vp_ptr.len)
             *
             * Let us do a linear search.
             *
             * */ 
            list1_len = ea->vpatch_p[i].len - m;
            list2_len = vp_ptr1->len;
            debug_body("comparing lists of length: %d and %d", list1_len, list2_len);
            
            dlist1 = ea->vpatch_p[i].domain_p + m;
            dlist2 = vp_ptr1->domain_p;

            found = false;
            for (j = 0; j < list1_len; j++)
            {
                del1 = dlist1[j];
                for (k = 0; k < list2_len; k++)
                {
                   if (del1 == dlist2[k])
                   {
                        found = true;
                        debug_body( "Common domain found (j: %d, k: %d): %d",
                                    j+m, k, dlist2[k]->domain_index );
                        break;
                   }
                }
                if(found)
                    break;
            }

            if(found)
            {
                ea->vpatch_p[i].bvert_index[m] = vp_ptr1->vert_index[k];
                err_check( (vp_ptr1->vert_index[k] == ea->vpatch_p[i].vert_index[j+m]),
                           clean_up, "%s", "Inconsisten vertex ordering found!");

                if(j != 0)
                {
                    dptr_tmp = ea->vpatch_p[i].domain_p[m];
                    ea->vpatch_p[i].domain_p[m] = ea->vpatch_p[i].domain_p[j+m];
                    ea->vpatch_p[i].domain_p[j+m] = dptr_tmp;
                    debug_body("dptr_tmp: %p", dptr_tmp);

                    tmp = ea->vpatch_p[i].vert_index[m];
                    ea->vpatch_p[i].vert_index[m] = ea->vpatch_p[i].vert_index[j+m];
                    ea->vpatch_p[i].vert_index[j+m] = tmp;
                    debug_body("tmp: %d", tmp);
                }

                other_vertex1 = (ea->vpatch_p[i].vert_index[m] + 1) % FEM2D_NV;
                other_vertex2 = (ea->vpatch_p[i].vert_index[m] + 2) % FEM2D_NV;

                next_vertex =   (ea->vpatch_p[i].bvert_index[m] == other_vertex1) 
                              ? other_vertex2 : other_vertex1;
                debug_body( "vertex index: %d, other vertex1: %d, "
                            "other vertex2: %d, next vertex: %d", 
                            ea->vpatch_p[i].vert_index[m], other_vertex1,
                            other_vertex2, next_vertex);
            }
            else
            {
                debug_body("%s", "No common domain found!");
                /* copy the list in new a lists in order to
                 * continue in the backward direction
                 * */ 
                domain_p_tmp    = calloc(ea->vpatch_p[i].len, sizeof(fem2d_te*));
                vert_index_tmp  = calloc(ea->vpatch_p[i].len, sizeof(matlib_index));
                bvert_index_tmp = calloc(ea->vpatch_p[i].len, sizeof(matlib_index));

                domain_p_tmp[0]    = ea->vpatch_p[i].domain_p[m-1];
                vert_index_tmp[0]  = ea->vpatch_p[i].vert_index[m-1];
                bvert_index_tmp[0] = next_vertex;

                for( n = 1; n < m; n++)
                {
                    domain_p_tmp[n]   = ea->vpatch_p[i].domain_p[m-n-1];
                    vert_index_tmp[n] = ea->vpatch_p[i].vert_index[m-n-1];
                    
                    other_vertex1 = (vert_index_tmp[n] + 1) % FEM2D_NV;
                    other_vertex2 = (vert_index_tmp[n] + 2) % FEM2D_NV;
                    next_vertex =   (ea->vpatch_p[i].bvert_index[m-n-1] == other_vertex1) 
                                  ? other_vertex2 : other_vertex1;

                    bvert_index_tmp[n] = next_vertex;
                }
                
                next_vertex =   (ea->vpatch_p[i].bvert_index[0] == other_vertex1) 
                              ? other_vertex1 : other_vertex2;
                for( n = m; n < ea->vpatch_p[i].len; n++)
                {
                    domain_p_tmp[n]   = ea->vpatch_p[i].domain_p[n];
                    vert_index_tmp[n] = ea->vpatch_p[i].vert_index[n];
                }

                matlib_free(ea->vpatch_p[i].domain_p);
                matlib_free(ea->vpatch_p[i].vert_index);
                matlib_free(ea->vpatch_p[i].bvert_index);

                ea->vpatch_p[i].domain_p     = domain_p_tmp;
                ea->vpatch_p[i].vert_index   = vert_index_tmp;
                ea->vpatch_p[i].bvert_index  = bvert_index_tmp;
                m = m-1;
            }
        }

        m = ea->vpatch_p[i].len - 1;
        other_vertex1 = (ea->vpatch_p[i].vert_index[m] + 1) % FEM2D_NV;
        other_vertex2 = (ea->vpatch_p[i].vert_index[m] + 2) % FEM2D_NV;

        next_vertex =   (ea->vpatch_p[i].bvert_index[m] == other_vertex1) 
                      ? other_vertex2 : other_vertex1;
        debug_body( "vert_index: %d, other vertex1: %d, "
                    "other vertex2: %d, next vertex: %d", 
                    ea->vpatch_p[i].vert_index[m], other_vertex1, 
                    other_vertex2, next_vertex);

        dptr = ea->vpatch_p[i].domain_p[m];
        
        vp_ptr1 = vp_ptr0 + dptr->nindex_p[next_vertex];

        list1_len = 1;
        list2_len = vp_ptr1->len;
        debug_body("comparing lists of length: %d and %d", list1_len, list2_len);
        
        dlist1 = ea->vpatch_p[i].domain_p;
        dlist2 = vp_ptr1->domain_p;

        found = false;
        for (j = 0; j < list1_len; j++)
        {
            del1 = dlist1[j];
            for (k = 0; k < list2_len; k++)
            {
               if (del1 == dlist2[k])
               {
                    found = true;
                    debug_body( "found common domain (j: %d, k: %d): %d",
                                j+m, k, dlist2[k]->domain_index);
                    break;
               }
            }
            if(found)
                break;
        }
        if(found)
        {
            ea->vpatch_p[i].point_enum = FEM2D_INTERIOR;
            debug_body("node: %d -> Interior point", i);
        }
        else
        {
            ea->vpatch_p[i].point_enum = FEM2D_BOUNDARY;
            debug_body("node: %d -> Boundary point", i);
        }
    }


    mcnt++;
    /* Create the node order array */ 
    matlib_index list_len;
    matlib_index* node_list = NULL;
    fem2d_err error;
    fem2d_vp *vp_nptr;
    for ( vp_nptr = ea->vpatch_p; 
          vp_nptr < ea->nr_nodes + ea->vpatch_p; 
          vp_nptr++)
    {
        errno = 0;
        if(vp_nptr->point_enum == FEM2D_BOUNDARY)
        {
            vp_nptr->node_order = calloc(vp_nptr->len+1, sizeof(matlib_int));
        }
        else
        {
            vp_nptr->node_order = calloc(vp_nptr->len, sizeof(matlib_int));
        }
        err_check( (vp_nptr->node_order == NULL), clean_up, 
                    "%s: Memory allocation failed!", strerror(errno));
    }

    for (i = 0; i < ea->nr_nodes; i++)
    {
        errno = 0;
        if(ea->vpatch_p[i].point_enum == FEM2D_BOUNDARY)
        {
            node_list = calloc(ea->vpatch_p[i].len + 1, sizeof(matlib_index));
        }
        else
        {
            node_list = calloc(ea->vpatch_p[i].len, sizeof(matlib_index));
        }
        err_check( (node_list == NULL), clean_up, 
                    "%s: Memory allocation failed!", strerror(errno));

        list_len = ea->vpatch_p[i].len;
        for (j = 0; j < list_len; j++)
        {
            dptr = ea->vpatch_p[i].domain_p[j];
            /* position of the domain in the vertex patch */ 
            dptr->pos_vpatch[ea->vpatch_p[i].vert_index[j]] = j; 
            node_list[j] = dptr->nindex_p[ea->vpatch_p[i].bvert_index[j]];
        }
        if(ea->vpatch_p[i].point_enum == FEM2D_BOUNDARY)
        {
            j = list_len - 1;
            debug_body("%s", "Boundary point encountered!");
            other_vertex1 = (ea->vpatch_p[i].vert_index[j] + 1) % FEM2D_NV;
            other_vertex2 = (ea->vpatch_p[i].vert_index[j] + 2) % FEM2D_NV;

            next_vertex =   (ea->vpatch_p[i].bvert_index[j] == other_vertex1) 
                          ? other_vertex2 : other_vertex1;
            node_list[j+1] = dptr->nindex_p[next_vertex];
            list_len++;
        }
        error = fem2d_sort_node_index(i, node_list, list_len, ea->vpatch_p[i].node_order);
        err_check( (error == FEM2D_FAILURE), clean_up, 
                    "%s", "Sorting of node indices failed!");
        BEGIN_DEBUG
            debug_print( "node: %d, patch length: %d, point_enum: %s",
                         i, ea->vpatch_p[i].len,
                         FEM2D_POINT_T_ENUM2STR(ea->vpatch_p[i].point_enum));
            for (j = 0; j < list_len; j++)
            {
                debug_print( "[%d]-> node index: %d, "
                             "node order: %d",
                             j, node_list[j], ea->vpatch_p[i].node_order[j]);
            }
        END_DEBUG
        matlib_free(node_list);
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    if (mcnt == 3)
    {
        matlib_free(node_list);
        for ( ; (vp_nptr < ea->vpatch_p) & (vp_nptr != NULL); 
                vp_nptr--)
        {
            matlib_free(vp_nptr->node_order);
        }
        mcnt --;
    }
    if (mcnt == 2)
    {
        for ( ; (vp_bptr < ea->vpatch_p) & (vp_bptr != NULL); 
                vp_bptr--)
        {
            matlib_free(vp_bptr->bvert_index);
        }
        mcnt --;
    }
    if (mcnt == 1)
    {
        for ( ; (vp_ptr < ea->vpatch_p) & (vp_ptr != NULL); 
                vp_ptr--)
        {
            matlib_free(vp_ptr->domain_p);
            matlib_free(vp_ptr->vert_index);
        }
        matlib_free(ea->vpatch_p);
    }
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;
}
/*============================================================================*/

fem2d_err fem2d_sort_node_index
(
    matlib_index  zeroth,
    matlib_index* ilist,
    matlib_index  ilen,
    matlib_int*   iorder
)
{
    debug_enter("Length of the index list: %d", ilen);

    err_check(    (ilist == NULL) 
               || (iorder == NULL), clean_up,
               "%s", "Null pointers encountered!");

    matlib_index s = 0, l = 0;
    matlib_index i, j, mcnt = 0;
    
    errno = 0;
    matlib_index* plist = calloc(ilen, sizeof(matlib_index));
    err_check( plist == NULL, clean_up, 
               "%s: Memory allocation failed!", strerror(errno));
    mcnt++;

    errno = 0;
    matlib_index* nlist = calloc(ilen, sizeof(matlib_index));
    err_check( nlist == NULL, clean_up, 
               "%s: Memory allocation failed!", strerror(errno));
    mcnt++;

    for (i = 0; i < ilen; i++)
    {
        if(zeroth > ilist[i])
        {
            iorder[i] = (-1);
            nlist[s] = i;
            s++;
        }
        else
        {
            iorder[i] = 1;
            plist[l]  = i;
            l++;
        }
    }
    debug_body( "positive ones: %d, negative ones: %d", l, s);
    
    errno = 0;
    plist = realloc(plist, l * sizeof(matlib_index));
    err_check( (plist == NULL) && (l > 0), clean_up, 
               "%s: Memory allocation failed!", strerror(errno));

    errno = 0;
    nlist = realloc(nlist, s * sizeof(matlib_index));
    err_check( (nlist == NULL) && (s > 0), clean_up, 
               "%s: Memory allocation failed!", strerror(errno));
    
    if(plist != NULL)
    {
        for (j = 0; j < l; j++)
        {
            for(i = 0; i < l; i++)
            {
                if ((ilist[plist[i]] < ilist[plist[j]]) && (i != j))
                {
                    iorder[plist[j]] += 1;
                }
            }
        }
    }

    if(nlist != NULL)
    {
        for (j = 0; j < s; j++)
        {
            for(i = 0; i < s; i++)
            {
                if ((i != j) && (ilist[nlist[i]] > ilist[nlist[j]]))
                {
                    iorder[nlist[j]] += (-1); 
                }
            }
        }
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    if (mcnt == 2)
    {
        matlib_free(nlist);
        mcnt --;
    }
    if (mcnt == 1)
    {
        matlib_free(plist);
    }
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;
}
/*============================================================================*/
static bool check_edge_mid_points(fem2d_vp* vp_ptr)
/* Check if the common segments have same midpoints as computed in each
 * neighbouring domain.
 * */
{
    matlib_index i, j, nr_edges = 2; 
    
    matlib_index mlength = FEM2D_DIM * nr_edges * vp_ptr->len;

    matlib_real* edge_mid_point = calloc(mlength, sizeof(matlib_real));
    err_check( edge_mid_point == NULL, clean_up, 
               "Failed to allocate memory for the edge mid-point"
               " list (length: %d)!", mlength/nr_edges);

    matlib_real* mptr = edge_mid_point;
    matlib_real *v1, *v2, *v3;
    fem2d_te* dptr;
    
    
    for ( j = 0; j < vp_ptr->len; j++)
    {
        dptr = vp_ptr->domain_p[j];
        v1 = *(dptr->vert_p + vp_ptr->vert_index[j]);
        v2 = *(dptr->vert_p + vp_ptr->bvert_index[j]);
        v3 = *(dptr->vert_p + 
              (FEM2D_NV - vp_ptr->bvert_index[j] - vp_ptr->vert_index[j]));
        
        mptr[FEM2D_INDEX_DIM1] = 0.5 * (v1[FEM2D_INDEX_DIM1] + v2[FEM2D_INDEX_DIM1]);
        mptr[FEM2D_INDEX_DIM2] = 0.5 * (v1[FEM2D_INDEX_DIM2] + v2[FEM2D_INDEX_DIM2]);

        mptr += FEM2D_DIM;
        mptr[FEM2D_INDEX_DIM1] = 0.5 * (v1[FEM2D_INDEX_DIM1] + v3[FEM2D_INDEX_DIM1]);
        mptr[FEM2D_INDEX_DIM2] = 0.5 * (v1[FEM2D_INDEX_DIM2] + v3[FEM2D_INDEX_DIM2]);
        mptr += FEM2D_DIM;
    }
    
    matlib_real err = 0, err1;
    mptr = edge_mid_point;
    
    for ( j = 0; j < vp_ptr->len-1; j++)
    {
        mptr += FEM2D_DIM;
        err1 = sqrt(   pow(fabs(mptr[FEM2D_INDEX_DIM1] - mptr[FEM2D_DIM + FEM2D_INDEX_DIM1]), 2) 
                     + pow(fabs(mptr[FEM2D_INDEX_DIM2] - mptr[FEM2D_DIM + FEM2D_INDEX_DIM2]), 2));
        mptr += FEM2D_DIM;
        err += err1; 
        warn_if( err1 > FEM2D_MESH_TOL, "Node: %d, err1: %0.16f", vp_ptr->node_index, err1);
    }
    
    matlib_real *first_edge, *last_edge; 
    if(vp_ptr->point_enum == FEM2D_INTERIOR)
    {
        first_edge = edge_mid_point;
        last_edge  = edge_mid_point + mlength - FEM2D_DIM;

        err1 = sqrt(   pow(fabs(first_edge[FEM2D_INDEX_DIM1] - last_edge[FEM2D_INDEX_DIM1]), 2) 
                     + pow(fabs(first_edge[FEM2D_INDEX_DIM2] - last_edge[FEM2D_INDEX_DIM2]), 2));
        err += err1; 
    }

    bool is_nbhd_OK = (err < FEM2D_MESH_TOL);
    warn_if(!is_nbhd_OK, "Node: %d, err: %0.16f", vp_ptr->node_index, err);

    matlib_free(edge_mid_point);
    return (is_nbhd_OK);

clean_up:
    return (false);
}

static bool check_pos_vpatch(fem2d_vp* vp_ptr)
{
    matlib_index j, pos_vpatch; 

    fem2d_te* dptr;
    bool is_pos_OK1, is_pos_OK = true;

    for ( j = 0; j < vp_ptr->len; j++)
    {
        dptr = vp_ptr->domain_p[j];

        /* check pos_vpatch */
        pos_vpatch = dptr->pos_vpatch[vp_ptr->vert_index[j]];
        is_pos_OK1 = (pos_vpatch == j);
        is_pos_OK  = (is_pos_OK1 && is_pos_OK);

        warn_if( !is_pos_OK1, "Vertex Node: %d, domain_p[%d]: %d, pos: %d", 
                 vp_ptr->node_index, j, dptr->domain_index, pos_vpatch);
    }
    return (is_pos_OK);
}

static fem2d_err list_boundary_nodes
(
    fem2d_vp*  vp_ptr, 
    matlib_nv* bd_nodes
)
{

    matlib_index j; 
    fem2d_te* dptr;
    matlib_index nnodes; /* Number of negative nodes */ 

    for ( j = 0; j < vp_ptr->len; j++)
    {
        /* Count negative nodes */
        if (vp_ptr->node_order[j] < 0)
            nnodes++;
    }
    if (vp_ptr->point_enum == FEM2D_BOUNDARY)
    {
        if (vp_ptr->node_order[j]<0)
            nnodes++;
    }

    /* Testing order of the boundary nodes of the vertex patch 
     * */ 
    matlib_index bd_nodes_len = 0;

    if (vp_ptr->point_enum == FEM2D_INTERIOR)
    {
        bd_nodes_len = vp_ptr->len + 1;
    }
    else
    {
        bd_nodes_len = vp_ptr->len + 2;
    }
    
    fem2d_err error = matlib_create_nv(bd_nodes_len, bd_nodes);
    err_check( error == FEM2D_FAILURE, clean_up,
               "Failed to allocate memory for the boundary node"
               " list (length: %d)!", bd_nodes_len);

    matlib_index j0 = nnodes; /* index of the zeroth node in the list */ 
    matlib_index index;
    bd_nodes->elem_p[j0] = vp_ptr->node_index;

    for ( j = 0; j < vp_ptr->len; j++)
    {
        dptr  = vp_ptr->domain_p[j];
        index = j0 + vp_ptr->node_order[j];
        bd_nodes->elem_p[index] = dptr->nindex_p[vp_ptr->bvert_index[j]];
    }
    if (vp_ptr->point_enum == FEM2D_BOUNDARY)
    {
        index = (FEM2D_NV - vp_ptr->bvert_index[j-1] - vp_ptr->vert_index[j-1]);
        bd_nodes->elem_p[j0 + vp_ptr->node_order[j]] = dptr->nindex_p[index];
    }
   
    return error;

clean_up:
    return error;
}

static bool check_boundary_nodes(fem2d_vp* vp_ptr)
{

    debug_enter( "Node: %d, point type: %s, patch length: %d",
                 vp_ptr->node_index,
                 FEM2D_POINT_T_ENUM2STR(vp_ptr->point_enum), 
                 vp_ptr->len);

    matlib_index j; 
    
    bool is_order_OK = true, is_order_OK1;

    matlib_nv bd_nodes;
    fem2d_err error = list_boundary_nodes(vp_ptr, &bd_nodes);
    err_check( error == FEM2D_FAILURE, clean_up, 
               "%s", "Failed to get boundary node list!");

    matlib_index last = bd_nodes.elem_p[0];

    for (j = 1; j < bd_nodes.len; j++)
    {
        is_order_OK1 = (last < bd_nodes.elem_p[j]);
        is_order_OK  = (is_order_OK && is_order_OK1);
        
        warn_if( !is_order_OK1, 
                 "Incorrect order, boundary nodes[%d]: %d",
                 j, bd_nodes.elem_p[j]);
        
        last = bd_nodes.elem_p[j];
    }
    matlib_free(bd_nodes.elem_p);
    return (is_order_OK);

clean_up:
    return (false);
}

fem2d_err fem2d_check_vp
(
    fem2d_ea ea,
    fem2d_vpinfo_t* test_info
)
{
    err_check( (ea.vpatch_p == NULL), clean_up,
               "%s", "Null pointer ecountered!");

    fem2d_vp *vp_ptr;

    test_info->nr_total = ea.nr_nodes;

    matlib_index mcnt = 0;

    errno = 0;
    matlib_index* ptr = calloc(ea.nr_nodes, sizeof(matlib_index));
    err_check( ptr == NULL, clean_up, 
               "%s: Memory allocation failed!", strerror(errno));
    mcnt++;

    errno = 0;
    FEM2D_TEST_PF* ptr1 = calloc(ea.nr_nodes, sizeof(FEM2D_TEST_PF));   
    err_check( ptr1 == NULL, clean_up, 
               "%s: Memory allocation failed!", strerror(errno));
    mcnt++;

    errno = 0;
    FEM2D_TEST_PF* ptr2 = calloc(ea.nr_nodes, sizeof(FEM2D_TEST_PF));
    err_check( ptr2 == NULL, clean_up, 
               "%s: Memory allocation failed!", strerror(errno));
    mcnt++;

    errno = 0;
    FEM2D_TEST_PF* ptr3 = calloc(ea.nr_nodes, sizeof(FEM2D_TEST_PF));
    err_check( ptr3 == NULL, clean_up, 
               "%s: Memory allocation failed!", strerror(errno));
    mcnt++;

    matlib_index* tmp_ptr = ptr;

    FEM2D_TEST_PF* tmp_ptr1 = ptr1; 
    FEM2D_TEST_PF* tmp_ptr2 = ptr2; 
    FEM2D_TEST_PF* tmp_ptr3 = ptr3;

    matlib_index nr_failed = 0;
    matlib_index nr_passed = 0;
    bool is_nbhd_OK, is_pos_OK, is_order_OK, is_vp_OK;
    
    for ( vp_ptr = ea.vpatch_p; 
          vp_ptr < ea.vpatch_p + ea.nr_nodes; vp_ptr++)
    {
        is_nbhd_OK  = check_edge_mid_points(vp_ptr);
        is_pos_OK   = check_pos_vpatch(vp_ptr);
        is_order_OK = check_boundary_nodes(vp_ptr);

        is_vp_OK = (is_nbhd_OK && is_pos_OK && is_order_OK);
        if ((!is_vp_OK) && (test_info->option == FEM2D_TEST_FAILED))
        {
            *tmp_ptr = vp_ptr->node_index;
            nr_failed ++;

            *tmp_ptr1 = BOOL2PF(is_nbhd_OK);
            *tmp_ptr2 = BOOL2PF(is_pos_OK);
            *tmp_ptr3 = BOOL2PF(is_order_OK);

            tmp_ptr++;
            tmp_ptr1++; 
            tmp_ptr2++; 
            tmp_ptr3++;

        }
        else if ((is_vp_OK) && (test_info->option == FEM2D_TEST_PASSED))
        {
            *tmp_ptr = vp_ptr->node_index;
            nr_passed ++;

            *tmp_ptr1 = BOOL2PF(is_nbhd_OK);
            *tmp_ptr2 = BOOL2PF(is_pos_OK);
            *tmp_ptr3 = BOOL2PF(is_order_OK);

            tmp_ptr++;
            tmp_ptr1++; 
            tmp_ptr2++; 
            tmp_ptr3++;

        }
        warn_if( !is_vp_OK, 
                 "Vertex-patch is inconsistent "
                 "(vertex node: %d, node type: %s, patch length: %d)!",
                 vp_ptr->node_index,
                 FEM2D_POINT_T_ENUM2STR(vp_ptr->point_enum), 
                 vp_ptr->len);
    }

    if (test_info->option == FEM2D_TEST_FAILED)
    {
        test_info->nr_failed = nr_failed;
        test_info->nr_passed = test_info->nr_total - nr_failed;
        debug_body("Nr of inconsistent vertex patches: %d", nr_failed);

        errno = 0;
        ptr = realloc(ptr, nr_failed * sizeof(matlib_index));
        err_check( (ptr == NULL) && (nr_failed > 0), clean_up, 
                   "%s: Memory reallocation failed!", strerror(errno));

        errno = 0;
        ptr1 = realloc( ptr1, nr_failed * sizeof(FEM2D_TEST_PF));   
        err_check( (ptr1 == NULL) && (nr_failed > 0), clean_up, 
                   "%s: Memory reallocation failed!", strerror(errno));

        errno = 0;
        ptr2 = realloc( ptr2, nr_failed * sizeof(FEM2D_TEST_PF));
        err_check( (ptr2 == NULL) && (nr_failed > 0), clean_up, 
                   "%s: Memory reallocation failed!", strerror(errno));
        errno = 0;
        ptr3 = realloc( ptr3, nr_failed * sizeof(FEM2D_TEST_PF));
        err_check( (ptr3 == NULL) && (nr_failed > 0), clean_up, 
                   "%s: Memory reallocation failed!", strerror(errno));
    }
    else if (test_info->option == FEM2D_TEST_PASSED)
    {
        test_info->nr_passed = nr_passed;
        test_info->nr_failed = test_info->nr_total - nr_passed;
        debug_body("Nr of inconsistent vertex patches: %d", nr_failed);

        errno = 0;
        ptr = realloc(ptr, nr_passed * sizeof(matlib_index));
        err_check( (ptr == NULL) && (nr_passed > 0), clean_up, 
                   "%s: Memory reallocation failed!", strerror(errno));

        errno = 0;
        ptr1 = realloc( ptr1, nr_passed * sizeof(FEM2D_TEST_PF));   
        err_check( (ptr1 == NULL) && (nr_passed > 0), clean_up, 
                   "%s: Memory reallocation failed!", strerror(errno));

        errno = 0;
        ptr2 = realloc( ptr2, nr_passed * sizeof(FEM2D_TEST_PF));
        err_check( (ptr2 == NULL) && (nr_passed > 0), clean_up, 
                   "%s: Memory reallocation failed!", strerror(errno));
        errno = 0;
        ptr3 = realloc( ptr3, nr_passed * sizeof(FEM2D_TEST_PF));
        err_check( (ptr3 == NULL) && (nr_passed > 0), clean_up, 
                   "%s: Memory reallocation failed!", strerror(errno));
    }

    test_info->edge_mid_points = ptr1;
    test_info->pos_vpatch      = ptr2;
    test_info->vpatch_boundary = ptr3;
    test_info->vpatch_index    = ptr;
    

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    if (mcnt == 4)
    {
        matlib_free(ptr3);
        mcnt--;
    }
    if (mcnt == 3)
    {
        matlib_free(ptr3);
        mcnt--;
    }
    if (mcnt == 2)
    {
        matlib_free(ptr3);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(ptr3);
        mcnt--;
    }
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;
}


matlib_index find_nr_digits(matlib_index num)
{
    matlib_index digits   = 1;
    matlib_index test_num = 1;
    while (test_num < num)
    {
        test_num = 10 * test_num;
        digits++;
    }
    return digits;
}

fem2d_err fem2d_write_vpinfo
(
    fem2d_vpinfo_t* test_info,
    fem2d_ea ea,
    char* file_name
)
{
    debug_enter("Writing vertex patch info into file: %s", file_name);

    matlib_index mcnt = 0;
    errno = 0;
    FILE* fp = fopen(file_name, "w+");
    err_check( fp == NULL, clean_up, 
               "%s: Memory reallocation failed!", strerror(errno));
    mcnt++;

    fprintf( fp, "Total Number of vertex Patches       : %d\n", test_info->nr_total);
    fprintf( fp, "Number of Inconsistent Vertex Patches: %d\n\n", test_info->nr_failed);

    fem2d_err error;
    matlib_nv bd_nodes;
    matlib_index i, j;
    fem2d_te* dptr;
    fem2d_vp* vp_ptr;

    matlib_index max_iter = 0;
    if (test_info->option == FEM2D_TEST_FAILED)
    {
        max_iter = test_info->nr_failed;
    }
    if (test_info->option == FEM2D_TEST_PASSED)
    {
        max_iter = test_info->nr_passed;
    }

    matlib_index fw = find_nr_digits(ea.nr_nodes);
    for (i = 0; i < max_iter; i++)
    {
        vp_ptr = ea.vpatch_p + test_info->vpatch_index[i];
        fprintf( fp, "[%2$*1$d.] Node Index: %3$*1$d | Point Type: %4$s | Patch Length: %5$d\n",
                 fw, i, vp_ptr->node_index, FEM2D_POINT_T_ENUM2STR(vp_ptr->point_enum), vp_ptr->len);
        for (j = 0; j < vp_ptr->len; j++)
        {
            dptr = vp_ptr->domain_p[j];
            VPINFO_PRINT( fp, "[%2$d] Domain Index[%9$s]: %3$*1$d | Nodes: [%4$*1$d, %5$*1$d, %6$*1$d] "
                          "| Vertex Index: %7$d | Index Boundary: %8$d", fw, j, 
                          dptr->domain_index, 
                          dptr->nindex_p[FEM2D_INDEX_V1], 
                          dptr->nindex_p[FEM2D_INDEX_V2], 
                          dptr->nindex_p[FEM2D_INDEX_V3], 
                          vp_ptr->vert_index[j], 
                          vp_ptr->bvert_index[j], (dptr->sign >0)? "+": "-");
        }
        error = list_boundary_nodes(vp_ptr, &bd_nodes);
        err_check( error == FEM2D_FAILURE, clean_up, 
                   "%s", "Failed to get boundary node list!");
        VPINFO_PRINT( fp, "Boundary Nodes (Including Vertex Node: %2$*1$d)", fw, vp_ptr->node_index);
        for (j = 0; j < bd_nodes.len; j++)
        {
            
            VPINFO_PRINT( fp, "\t[%2$d]: %3$*1$d Point Type: %4$s", 
                          fw, j, bd_nodes.elem_p[j], 
                          FEM2D_POINT_T_ENUM2STR(ea.vpatch_p[bd_nodes.elem_p[j]].point_enum));
        }
        matlib_free(bd_nodes.elem_p);

        VPINFO_PRINT( fp, "Test Edge Mid-points    : %s",
                      FEM2D_TEST_PF_ENUM2STR(test_info->edge_mid_points[i]));
        VPINFO_PRINT( fp, "Test Domain Position    : %s",
                      FEM2D_TEST_PF_ENUM2STR(test_info->pos_vpatch[i]));
        VPINFO_PRINT( fp, "Test Boundary Node Order: %s",
                      FEM2D_TEST_PF_ENUM2STR(test_info->vpatch_boundary[i]));
        fprintf(fp, "\n");
    }

    fclose(fp);
    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    if (mcnt == 1) 
        matlib_free(bd_nodes.elem_p);
    fclose(fp);
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;
}

void fem2d_free_vpinfo(fem2d_vpinfo_t* test_info)
{
    matlib_free(test_info->edge_mid_points);
    matlib_free(test_info->pos_vpatch     );
    matlib_free(test_info->vpatch_boundary);
    matlib_free(test_info->vpatch_index   );
}

/*============================================================================*/


void fem2d_free_ea(fem2d_ea ea)
{
    fem2d_te* dptr;

    for ( dptr = ea.elem_p; 
          (dptr < (ea.elem_p + ea.len)) && (dptr != NULL);
          dptr++)
    {
        matlib_free(dptr->vert_p);
        matlib_free(dptr->pos_vpatch);
        matlib_free(dptr->nindex_p);
        matlib_free(dptr->jmat);
        matlib_free(dptr->ijmat);
    }
    fem2d_vp* vp_ptr;
    for ( vp_ptr = ea.vpatch_p; 
          (vp_ptr < ea.vpatch_p + ea.nr_nodes) & (vp_ptr != NULL); 
          vp_ptr++)
    {
        matlib_free(vp_ptr->domain_p);
        matlib_free(vp_ptr->vert_index);
        matlib_free(vp_ptr->bvert_index);
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
    matlib_real* nptr0 = ea.node_p;

    matlib_index* ptr;
    
    matlib_index offset; 
    matlib_index vsize = (matlib_index)2*sizeof(matlib_real*);

    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); 
          dptr++, ia+=3)
    {

        /* Get the first vertex */ 
        offset = (matlib_index)(((matlib_index)*(dptr->vert_p + FEM2D_INDEX_V1) 
                               - (matlib_index)nptr0)/vsize);
        *(ia+0) = offset;


        /* Get the second vertex */
        offset = (matlib_index)(((matlib_index)*(dptr->vert_p + FEM2D_INDEX_V2) 
                               - (matlib_index)nptr0)/vsize);
        *(ia+1) = offset;

        /* Get the third vertex */
        offset = (matlib_index)(((matlib_index)*(dptr->vert_p + FEM2D_INDEX_V3) 
                               - (matlib_index)nptr0)/vsize);
        *(ia+2) = offset;

    }
    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;

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
    err_check(    (ea.elem_p==NULL) 
               || (cen.elem_p == NULL), clean_up,
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
        vptr = *(dptr->vert_p + FEM2D_INDEX_V1);
        cptr[FEM2D_INDEX_DIM1] += vptr[FEM2D_INDEX_DIM1]/3.0;  
        cptr[FEM2D_INDEX_DIM2] += vptr[FEM2D_INDEX_DIM2]/3.0; 

        /* Get the second vertex */
        vptr = *(dptr->vert_p + FEM2D_INDEX_V2);
        cptr[FEM2D_INDEX_DIM1] += vptr[FEM2D_INDEX_DIM1]/3.0;  
        cptr[FEM2D_INDEX_DIM2] += vptr[FEM2D_INDEX_DIM2]/3.0; 

        /* Get the third vertex */
        vptr = *(dptr->vert_p + FEM2D_INDEX_V3);
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
        vptr = *(dptr->vert_p + FEM2D_INDEX_V1);
        *ptr = vptr[FEM2D_INDEX_DIM1]; ptr++;
        *ptr = vptr[FEM2D_INDEX_DIM2]; ptr++;

        /* Get the second vertex */
        vptr = *(dptr->vert_p + FEM2D_INDEX_V2);
        *ptr = vptr[FEM2D_INDEX_DIM1]; ptr++;
        *ptr = vptr[FEM2D_INDEX_DIM2]; ptr++;

        /* Get the third vertex */
        vptr = *(dptr->vert_p + FEM2D_INDEX_V3);
        *ptr = vptr[FEM2D_INDEX_DIM1]; ptr++;
        *ptr = vptr[FEM2D_INDEX_DIM2];
        
        error = matlib_xgemm( alpha, vphi_tmp, nodes_mat, beta, x_tmp );
        err_check( (error == FEM2D_FAILURE), clean_up, 
                   "%s", "Matrix multiplication failed!");
        
        //DEBUG_PRINT_XM(x_tmp, "Points in domain[%d]: ", dptr->domain_index);

        x_tmp.elem_p = x_tmp.elem_p + FEM2D_DIM * vphi.lenc;
    }
    
    matlib_free(nodes_mat.elem_p);
    debug_exit( "Total nr. of sampling points: %d, Exit Status: %s",
                x->len, "SUCCESS" );
    return FEM2D_SUCCESS;
    
clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;

} /* End of fem2d_ref2mesh */ 
/*============================================================================*/
fem2d_err fem2d_xinterp
(
    const fem2d_ea  ea,
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
    
    err_check(    (ea.elem_p      == NULL)   
               || (u_nodes.elem_p == NULL) 
               || (vphi.elem_p    == NULL), clean_up, 
               "%s", "Null pointer(s) ecountered!");
    
    /* Length of u_nodes = nr. of nodes */ 
    err_check( (u_nodes.len != ea.nr_nodes ), clean_up, 
               "Dimension mismatch (u_nodes: %d, nr nodes: )",
               u_nodes.len, ea.nr_nodes);

    /* Length of u_interp = nr of quad. nodes x nr of domains 
     * */ 
    err_check( (u_interp.len != (vphi.lenc * ea.len) ), clean_up, 
               "Dimension mismatch (u_interp: %d)", u_interp.len);
    
    matlib_real alpha = 1.0;
    matlib_real beta  = 0.0;

    matlib_index INDEX_REAL = 0;

    matlib_index mcnt = 0;
    matlib_xv u1_tmp;    
    fem2d_err error = matlib_create_xv( FEM2D_NV, &u1_tmp, MATLIB_COL_VECT);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Memory allocation failed!");
    mcnt++;

    matlib_xv u2_tmp = { .len    = vphi.lenc, 
                         .type   = MATLIB_COL_VECT,
                         .elem_p = u_interp.elem_p };

    fem2d_te* dptr; /* domain pointer */

    matlib_real* ptr; 
    matlib_real* uptr; /* vertex pointer */ 
    matlib_real* uptr0 = u_nodes.elem_p; /* base address */ 

    /* Convert vphi to row-major format */ 
    matlib_xm vphi_tmp = { .order = MATLIB_ROW_MAJOR,
                           .op    = MATLIB_TRANS,
                           .lenr  = vphi.lenc,
                           .lenc  = FEM2D_NV,
                           .elem_p = vphi.elem_p };


    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); dptr++)
    {
        ptr = u1_tmp.elem_p;

        /* Get the first vertex */ 
        uptr = uptr0 + dptr->nindex_p[FEM2D_INDEX_V1]; 
        *ptr = uptr[INDEX_REAL]; ptr++;
        debug_body( "domain: %d, u: %0.16f%, offset: %d", 
                    dptr->domain_index, uptr[INDEX_REAL], 
                    dptr->nindex_p[FEM2D_INDEX_V1] );

        /* Get the second vertex */
        uptr = uptr0 + dptr->nindex_p[FEM2D_INDEX_V2]; 
        *ptr = uptr[INDEX_REAL]; ptr++;
        debug_body( "domain: %d, u: %0.16f, offset: %d", 
                    dptr->domain_index, uptr[INDEX_REAL],
                    dptr->nindex_p[FEM2D_INDEX_V2]);

        /* Get the third vertex */
        uptr = uptr0 + dptr->nindex_p[FEM2D_INDEX_V3]; 
        *ptr = uptr[INDEX_REAL]; ptr++;
        debug_body( "domain: %d, u: %0.16f, offset: %d", 
                    dptr->domain_index, uptr[INDEX_REAL],
                    dptr->nindex_p[FEM2D_INDEX_V3] );
        
        error = matlib_xgemv( alpha, vphi_tmp, u1_tmp, beta, u2_tmp );
        err_check( (error == FEM2D_FAILURE), clean_up, 
                   "%s", "Matrix-vector multiplication failed!");
  
        u2_tmp.elem_p = u2_tmp.elem_p + vphi.lenc; 

    }
    matlib_free(u1_tmp.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;
    
clean_up:
    if (mcnt==1)
        matlib_free(u1_tmp.elem_p);
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;

}/* End of fem2d_xinterp */ 


fem2d_err fem2d_zinterp
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
    
    err_check(    (ea.elem_p      == NULL)   
               || (u_nodes.elem_p == NULL) 
               || (vphi.elem_p    == NULL), clean_up, 
               "%s", "Null pointers ecountered!");
    
    /* Length of u_nodes = nr. of nodes */ 
    err_check( (u_nodes.len != ea.nr_nodes ), clean_up, 
               "Dimension mismatch (u_nodes: %d, nr nodes: )",
               u_nodes.len, ea.nr_nodes);

    /* Length of u_interp = nr of quad. nodes x nr of domains 
     * */ 
    err_check( (u_interp.len != (vphi.lenc * ea.len) ), clean_up, 
               "Dimension mismatch (u_interp: %d)", u_interp.len);
    
    matlib_real alpha = 1.0;
    matlib_real beta  = 0.0;

    matlib_index INDEX_REAL  = 0;
    matlib_index INDEX_IMAG  = 1;
    matlib_index COMPLEX_DIM = 2; /* Dimension of a compelex plane */ 

    matlib_index mcnt = 0;
    matlib_xm u1_tmp;    
    fem2d_err error = matlib_create_xm( FEM2D_NV, COMPLEX_DIM, &u1_tmp, 
                                        MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Memory allocation failed!");
    mcnt++;

    matlib_xm u2_tmp = { .lenc   = vphi.lenc, 
                         .lenr   = COMPLEX_DIM, 
                         .order  = MATLIB_ROW_MAJOR,
                         .op     = MATLIB_NO_TRANS,
                         .elem_p = (matlib_real*)u_interp.elem_p };

    fem2d_te* dptr; /* domain pointer */

    matlib_real* ptr; 
    matlib_real* uptr; /* vertex pointer */ 
    matlib_real* uptr0 = (matlib_real*)u_nodes.elem_p; /* base address */ 

    matlib_xm vphi_tmp = { .order = MATLIB_ROW_MAJOR,
                           .op    = MATLIB_TRANS,
                           .lenr  = vphi.lenc,
                           .lenc  = FEM2D_NV,
                           .elem_p = vphi.elem_p };

    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); dptr++)
    {
        ptr = u1_tmp.elem_p;

        /* Get the first vertex */ 
        uptr = uptr0 + COMPLEX_DIM * dptr->nindex_p[FEM2D_INDEX_V1]; 
        *ptr = uptr[INDEX_REAL]; ptr++;
        *ptr = uptr[INDEX_IMAG]; ptr++;

        /* Get the second vertex */
        uptr = uptr0 + COMPLEX_DIM * dptr->nindex_p[FEM2D_INDEX_V2]; 
        *ptr = uptr[INDEX_REAL]; ptr++;
        *ptr = uptr[INDEX_IMAG]; ptr++;

        /* Get the third vertex */
        uptr = uptr0 + COMPLEX_DIM * dptr->nindex_p[FEM2D_INDEX_V3]; 
        *ptr = uptr[INDEX_REAL]; ptr++;
        *ptr = uptr[INDEX_IMAG];
        
        error = matlib_xgemm( alpha, vphi_tmp, u1_tmp, beta, u2_tmp );
        err_check( (error == FEM2D_FAILURE), clean_up, 
                   "%s", "Matrix multiplication failed!");
  
        u2_tmp.elem_p = u2_tmp.elem_p + COMPLEX_DIM * vphi.lenc; 

    }
    matlib_free(u1_tmp.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;
    
clean_up:
    if (mcnt==1)
        matlib_free(u1_tmp.elem_p);
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;

}/* End of fem2d_zinterp */ 

/*============================================================================*/
matlib_real fem2d_xnormL2
(
    fem2d_ea  ea,
    matlib_xv u_qnodes,
    matlib_xv quadW
)
{
    debug_enter( "nr of domains: %d, "
                 "nr of nodes: %d, "
                 "nr of quadrature nodes: %d",
                 ea.len, ea.nr_nodes, quadW.len);
    
    err_check(    (ea.elem_p       == NULL) 
               || (u_qnodes.elem_p == NULL) 
               || (quadW.elem_p    == NULL), clean_up,
               "%s", "Null pointers ecountered!");

    err_check( (u_qnodes.len != quadW.len * ea.len), clean_up, 
               "%s", "Dimension mis-match!");

    matlib_real* ptr = u_qnodes.elem_p;
    
    matlib_real norm = 0;
    matlib_real jacob;

    for ( matlib_index i = 0; i < ea.len; i++)
    {
        jacob = (ea.elem_p[i]).jacob;
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

} /* End of fem2d_xnormL2 */ 

matlib_real fem2d_znormL2
(
    fem2d_ea  ea,
    matlib_zv u_qnodes,
    matlib_xv quadW
)
{
    debug_enter( "nr of domains: %d, "
                 "nr of nodes: %d, "
                 "nr of quadrature nodes: %d",
                 ea.len, ea.nr_nodes, quadW.len);
    
    err_check(    (ea.elem_p       == NULL) 
               || (u_qnodes.elem_p == NULL) 
               || (quadW.elem_p    == NULL), clean_up,
               "%s", "Null pointers ecountered!");

    err_check( (u_qnodes.len != quadW.len * ea.len), clean_up, 
               "%s", "Dimension mis-match!");

    matlib_complex* ptr = u_qnodes.elem_p;
    
    matlib_real norm = 0;
    matlib_real jacob;

    for ( matlib_index i = 0; i< ea.len; i++)
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
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

}/* End of fem2d_znormL2 */ 

/*============================================================================*/
matlib_real fem2d_xiprod 
/* Inner product iprod = (u, conj(v))_{\Omega} */ 
(
    fem2d_ea  ea,
    matlib_xv u_qnodes,
    matlib_xv v_qnodes,
    matlib_xv quadW
)
{
    debug_enter( "nr of domains: %d, "
                 "nr of nodes: %d, "
                 "nr of quadrature nodes: %d",
                 ea.len, ea.nr_nodes, quadW.len);
    
    err_check(    (ea.elem_p       == NULL)    
               || (u_qnodes.elem_p == NULL) 
               || (v_qnodes.elem_p == NULL) 
               || (quadW.elem_p    == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check( (u_qnodes.len != v_qnodes.len),
               clean_up, "Dimension mis-match for field vectors (u: %d, v: %d)!",
               v_qnodes.len, u_qnodes.len );

    err_check( (v_qnodes.len != quadW.len * ea.len), clean_up, 
               "%s", "Dimension mis-match!");

    matlib_real* uptr = u_qnodes.elem_p;
    matlib_real* vptr = v_qnodes.elem_p;
    
    matlib_real iprod = 0;
    matlib_real jacob;

    for ( matlib_index i = 0; i< ea.len; i++)
    {
        jacob = (ea.elem_p[i]).jacob;
        for (matlib_index j = 0; j< quadW.len; j++, uptr++, vptr++)
        {
            iprod += (*uptr * *vptr * quadW.elem_p[j] * jacob);
        }
    }
    debug_exit("Exit Status: %s", "SUCCESS" );
    return iprod;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

} /* End of fem2d_xiprod */ 


matlib_real fem2d_ziprod 
/* Inner product iprod = (u, conj(v))_{\Omega} */ 
(
    fem2d_ea  ea,
    matlib_zv u_qnodes,
    matlib_zv v_qnodes,
    matlib_xv quadW
)
{
    debug_enter( "nr of domains: %d, "
                 "nr of nodes: %d, "
                 "nr of quadrature nodes: %d",
                 ea.len, ea.nr_nodes, quadW.len);
    
    err_check(    (ea.elem_p       == NULL)    
               || (u_qnodes.elem_p == NULL) 
               || (v_qnodes.elem_p == NULL) 
               || (quadW.elem_p    == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check( (u_qnodes.len != v_qnodes.len),
               clean_up, "Dimension mis-match for field vectors (u: %d, v: %d)!",
               v_qnodes.len, u_qnodes.len );

    err_check( (v_qnodes.len != quadW.len * ea.len), clean_up, 
               "%s", "Dimension mis-match!");

    matlib_complex* uptr = u_qnodes.elem_p;
    matlib_complex* vptr = v_qnodes.elem_p;
    
    matlib_real iprod = 0;
    matlib_real jacob;

    for ( matlib_index i = 0; i< ea.len; i++)
    {
        jacob = (ea.elem_p[i]).jacob;
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

} /* fem2d_ziprod */ 

/*============================================================================*/

fem2d_err fem2d_quadP
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

    fem2d_err error; 
    error = matlib_create_xm( FEM2D_NV, vphi.lenc, quadP, 
                              MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    err_check( (error == FEM2D_FAILURE), clean_up, 
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

    ptr1 = quadW.elem_p;

    for ( ; ptr1 < (quadW.elem_p + quadW.len);
            ptr1++, ptr2++, ptr3++)
    {
        *ptr3 = *ptr1 * *ptr2;
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    matlib_free(quadP->elem_p);
    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;

}/* End of fem2d_quadP */ 

matlib_real fem2d_check_quadP
(
    fem2d_cc     xi,
    matlib_xv    quadW,
    matlib_index m,
    matlib_index n
)
/* 
 *  - (xi1*xi1^m*xi2^n)/2 - (xi1^m*xi2*xi2^n)/2
 *        (xi1^m*xi2^n)/2 + (xi1*xi1^m*xi2^n)/2
 *        (xi1^m*xi2^n)/2 + (xi1^m*xi2*xi2^n)/2
 *
 * */
{

    debug_enter("m: %d, n: %d, nr of quad points: %d", m, n, xi.len);

    matlib_real q1, q2, q3;

    q1 = - fem2d_poly_symint( m + 1, n + 0)/2.0
         - fem2d_poly_symint( m + 0, n + 1)/2.0;
    debug_body("q1: %0.16f", q1);

    q2 =   fem2d_poly_symint( m + 0, n + 0)/2.0
         + fem2d_poly_symint( m + 1, n + 0)/2.0;
    debug_body("q2: %0.16f", q2);

    q3 =   fem2d_poly_symint( m + 0, n + 0)/2.0
         + fem2d_poly_symint( m + 0, n + 1)/2.0;
    debug_body("q3: %0.16f", q3);

    matlib_index mcnt = 0;
    fem2d_err error;
    matlib_xv q, phi;
    error = matlib_create_xv( 3, &q, MATLIB_COL_VECT);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Memory allocation failed!");
    mcnt++; /* 1 */ 
    error = matlib_create_xv( xi.len, &phi, MATLIB_COL_VECT);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Memory allocation failed!");
    mcnt++; /* 2 */ 

    matlib_real *xiptr;
    matlib_real *phi_ptr = phi.elem_p;
    for ( xiptr = xi.elem_p; 
          xiptr < xi.elem_p + FEM2D_DIM * xi.len; 
          xiptr += FEM2D_DIM, phi_ptr++)
    {
        *phi_ptr = pow(xiptr[FEM2D_INDEX_DIM1], m) * pow(xiptr[FEM2D_INDEX_DIM2], n);
    }

    matlib_xm vphi, quadP;
    error = fem2d_refbasis(xi, &vphi);
    //error = FEM2D_FAILURE;
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Failed to compute the values of reference basis functions!");

    mcnt++; /* 3 */

    error = fem2d_quadP(vphi, quadW, &quadP);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Failed to compute the quadrature matrix for projection!");
    mcnt++; /* 4 */ 

    error = matlib_xgemv( 1.0, quadP, phi, 0.0, q);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Matrix product failed!");


    DEBUG_PRINT_XV(q, "%s:", "Gauss quadrature");
    matlib_real err = 0;
    err += ((q1 - q.elem_p[FEM2D_INDEX_V1]) * (q1 - q.elem_p[FEM2D_INDEX_V1]));
    err += ((q2 - q.elem_p[FEM2D_INDEX_V2]) * (q2 - q.elem_p[FEM2D_INDEX_V2]));
    err += ((q3 - q.elem_p[FEM2D_INDEX_V3]) * (q3 - q.elem_p[FEM2D_INDEX_V3]));

    matlib_real norm =   (q1 * q1) 
                       + (q2 * q2)
                       + (q3 * q3);

    matlib_real e_relative = sqrt(err/norm);
    debug_body("Relative error: %0.16f", e_relative);
    debug_exit("Exit Status: %s", "SUCCESS" );

    matlib_free(quadP.elem_p);
    matlib_free(vphi.elem_p);
    matlib_free(phi.elem_p);
    matlib_free(q.elem_p);
    return e_relative;

clean_up:
    
    if (mcnt == 4)
    {
        matlib_free(quadP.elem_p);
        mcnt--;
    }
    if (mcnt == 3)
    {
        matlib_free(vphi.elem_p);
        mcnt--;
    }
    if (mcnt == 2)
    {
        matlib_free(phi.elem_p);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(q.elem_p);
    }

    debug_body("Relative error: %0.16f", MATLIB_NAN);
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;
}
 
/*============================================================================*/
fem2d_err fem2d_xprj
(
    fem2d_ea  ea,
    matlib_xv u_qnodes, /* values at the quadrature nodes */ 
    matlib_xm quadP,
    matlib_xv u_prj

)
{
    debug_enter( "nr of domains: %d, "
                 "quadP: %d-by-%d, ",
                 ea.len, quadP.lenc, quadP.lenr);
    
    err_check(    (ea.elem_p       == NULL)    
               || (u_qnodes.elem_p == NULL) 
               || (quadP.elem_p    == NULL)  
               || (u_prj.elem_p    == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_qnodes.len != quadP.lenr * ea.len) 
               || (quadP.lenc   != FEM2D_NV) 
               || (u_prj.len    != ea.nr_nodes), clean_up, 
               "Dimension mis-match "
               "(u_qnodes: %d, quadP: %d-by-%d, u_prj: %d, nr nodes: %d)!",
               u_qnodes.len, quadP.lenc, quadP.lenr, u_prj.len, ea.nr_nodes);

    matlib_index INDEX_REAL  = 0;

    matlib_xv u1_tmp = { .len    = quadP.lenr, 
                         .type   = MATLIB_COL_VECT,
                         .elem_p = u_qnodes.elem_p };

    matlib_index mcnt = 0;
    matlib_xv u2_tmp;    
    fem2d_err error = matlib_create_xv( FEM2D_NV, &u2_tmp, MATLIB_COL_VECT);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Memory allocation for projection failed!");
    mcnt++;

    matlib_real alpha = 1.0;
    matlib_real beta  = 0.0;

    fem2d_te* dptr; /* domain pointer */
    matlib_real* uprj_ptr0 = u_prj.elem_p; /* base address */ 
    matlib_real* uprj_ptr;

    /* Clear the existing values */
    for ( uprj_ptr = uprj_ptr0; 
          uprj_ptr < uprj_ptr0 + u_prj.len; uprj_ptr++)
    {
        *uprj_ptr = 0.0;
    }

    for ( dptr = ea.elem_p; dptr < (ea.len + ea.elem_p); dptr++)
    {
        alpha = dptr->jacob;
        error = matlib_xgemv( alpha, quadP, u1_tmp, beta, u2_tmp );
        err_check( (error == FEM2D_FAILURE), clean_up, 
                   "%s", "Matrix multiplication failed!");
        
        /* Get the first vertex */ 
        uprj_ptr = uprj_ptr0 + dptr->nindex_p[FEM2D_INDEX_V1]; 

        uprj_ptr[INDEX_REAL] += u2_tmp.elem_p[INDEX_REAL + FEM2D_INDEX_V1];

        /* Get the second vertex */ 
        uprj_ptr = uprj_ptr0 + dptr->nindex_p[FEM2D_INDEX_V2]; 

        uprj_ptr[INDEX_REAL] += u2_tmp.elem_p[INDEX_REAL + FEM2D_INDEX_V2];

        /* Get the third vertex */ 
        uprj_ptr = uprj_ptr0 + dptr->nindex_p[FEM2D_INDEX_V3]; 

        uprj_ptr[INDEX_REAL] += u2_tmp.elem_p[INDEX_REAL + FEM2D_INDEX_V3];


        u1_tmp.elem_p = u1_tmp.elem_p + quadP.lenr; 
    }


    matlib_free(u2_tmp.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;

clean_up:

    if (mcnt == 1)
        matlib_free(u2_tmp.elem_p);

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;

} /* End of fem2d_xprj */ 


fem2d_err fem2d_zprj
(
    fem2d_ea  ea,
    matlib_zv u_qnodes, /* values at the quadrature nodes */ 
    matlib_xm quadP,
    matlib_zv u_prj

)
{
    debug_enter( "nr of domains: %d, "
                 "quadP: %d-by-%d, ",
                 ea.len, quadP.lenc, quadP.lenr);
    
    err_check(    (ea.elem_p       == NULL)    
               || (u_qnodes.elem_p == NULL) 
               || (quadP.elem_p    == NULL)  
               || (u_prj.elem_p    == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_qnodes.len != quadP.lenr * ea.len) 
               || (quadP.lenc   != FEM2D_NV)
               || (u_prj.len    != ea.nr_nodes), clean_up, 
               "Dimension mis-match (u_qnodes: %d, quadP: %d-by-%d)!",
               u_qnodes.len, quadP.lenc, quadP.lenr);

    matlib_index INDEX_REAL  = 0;
    matlib_index INDEX_IMAG  = 1;
    matlib_index COMPLEX_DIM = 2; /* Dimension of a compelex plane */ 

    matlib_xm u1_tmp = { .lenc   = quadP.lenr, 
                         .lenr   = COMPLEX_DIM, 
                         .order  = MATLIB_ROW_MAJOR,
                         .op     = MATLIB_NO_TRANS,
                         .elem_p = (matlib_real*)u_qnodes.elem_p };

    matlib_index mcnt = 0;
    matlib_xm u2_tmp;    
    fem2d_err error = matlib_create_xm( FEM2D_NV, COMPLEX_DIM, &u2_tmp, 
                                        MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Memory allocation for projection failed!");
    mcnt++;

    matlib_real alpha = 1.0;
    matlib_real beta  = 0.0;

    fem2d_te* dptr; /* domain pointer */
    matlib_real* uprj_ptr0 = (matlib_real*)u_prj.elem_p; /* base address */ 
    matlib_real* uprj_ptr;

    /* Clear the existing values */
    for ( uprj_ptr = uprj_ptr0; 
          uprj_ptr < uprj_ptr0 + u_prj.len; uprj_ptr++)
    {
        *uprj_ptr = 0.0;
    }


    for ( dptr = ea.elem_p; dptr < (ea.len + ea.elem_p); dptr++)
    {
        alpha = dptr->jacob;
        error = matlib_xgemm( alpha, quadP, u1_tmp, beta, u2_tmp );
        err_check( (error == FEM2D_FAILURE), clean_up, 
                   "%s", "Matrix multiplication failed!");
        
        /* Get the first vertex */ 
        uprj_ptr = uprj_ptr0 + COMPLEX_DIM * dptr->nindex_p[FEM2D_INDEX_V1]; 

        uprj_ptr[INDEX_REAL] += u2_tmp.elem_p[INDEX_REAL + FEM2D_INDEX_V1 * COMPLEX_DIM];
        uprj_ptr[INDEX_IMAG] += u2_tmp.elem_p[INDEX_IMAG + FEM2D_INDEX_V1 * COMPLEX_DIM];

        /* Get the second vertex */ 
        uprj_ptr = uprj_ptr0 + COMPLEX_DIM * dptr->nindex_p[FEM2D_INDEX_V2]; 

        uprj_ptr[INDEX_REAL] += u2_tmp.elem_p[INDEX_REAL + FEM2D_INDEX_V2 * COMPLEX_DIM];
        uprj_ptr[INDEX_IMAG] += u2_tmp.elem_p[INDEX_IMAG + FEM2D_INDEX_V2 * COMPLEX_DIM];

        /* Get the third vertex */ 
        uprj_ptr = uprj_ptr0 + COMPLEX_DIM * dptr->nindex_p[FEM2D_INDEX_V3]; 

        uprj_ptr[INDEX_REAL] += u2_tmp.elem_p[INDEX_REAL + FEM2D_INDEX_V3 * COMPLEX_DIM];
        uprj_ptr[INDEX_IMAG] += u2_tmp.elem_p[INDEX_IMAG + FEM2D_INDEX_V3 * COMPLEX_DIM];


        u1_tmp.elem_p = u1_tmp.elem_p + COMPLEX_DIM * quadP.lenr; 
    }


    matlib_free(u2_tmp.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;

clean_up:

    if (mcnt == 1)
        matlib_free(u2_tmp.elem_p);

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;

} /* End of fem2d_zprj */ 

/*============================================================================*/
fem2d_err fem2d_NB_xprj
/* Projection from nodal basis representation */ 
(
    fem2d_ea  ea,
    matlib_xv u_nodes, /* values at the nodes */ 
    matlib_xv u_prj
)
{
    debug_enter( "nr of domains: %d", ea.len);
    
    err_check(    (ea.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL)
               || (u_prj.elem_p   == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != ea.nr_nodes) 
               || (u_prj.len   != ea.nr_nodes), clean_up, 
               "Dimension mis-match (u_nodes: %d, u_prj: %d, nr nodes: %d)!",
               u_nodes.len, u_prj.len, ea.nr_nodes);

    fem2d_te* dptr; /* domain pointer */
    
    matlib_real* uptr0     = u_nodes.elem_p; /* base address */ 
    matlib_real* uprj_ptr0 = u_prj.elem_p; /* base address */ 
    matlib_real* uprj_ptr;
    matlib_real* uptr;
    matlib_index offset1, offset2, offset3;
    matlib_real jacob;
    matlib_real f1, f2, f3;

    /* Clear the existing values */
    for ( uprj_ptr = uprj_ptr0; 
          uprj_ptr < uprj_ptr0 + u_prj.len; uprj_ptr++)
    {
        *uprj_ptr = 0.0;
    }

    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); dptr++)
    {
        jacob = dptr->jacob;
        
        /* Get the first vertex */ 

        offset1 = dptr->nindex_p[FEM2D_INDEX_V1];
        offset2 = dptr->nindex_p[FEM2D_INDEX_V2];
        offset3 = dptr->nindex_p[FEM2D_INDEX_V3];

        uptr = uptr0 + offset1; f1 = *uptr * jacob;
        uptr = uptr0 + offset2; f2 = *uptr * jacob;
        uptr = uptr0 + offset3; f3 = *uptr * jacob;

        /* ===== */ 
        uprj_ptr = uprj_ptr0 + offset1;
        
        (*uprj_ptr) += (   MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V1] * f1
                         + MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V2] * f2
                         + MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V3] * f3);
        
        /* ===== */ 
        uprj_ptr = uprj_ptr0 + offset2;
        
        (*uprj_ptr) += (   MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V1] * f1
                         + MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V2] * f2
                         + MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V3] * f3);

        /* ===== */ 
        uprj_ptr = uprj_ptr0 + offset3;
        
        (*uprj_ptr) += (   MEMI[FEM2D_INDEX_V3][FEM2D_INDEX_V1] * f1
                         + MEMI[FEM2D_INDEX_V3][FEM2D_INDEX_V2] * f2
                         + MEMI[FEM2D_INDEX_V3][FEM2D_INDEX_V3] * f3);
    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;

} /* End of fem2d_NB_xprj */ 



fem2d_err fem2d_NB_zprj
(
    fem2d_ea  ea,
    matlib_zv u_nodes, /* values at the nodes */ 
    matlib_zv u_prj
)
{
    debug_enter( "nr of domains: %d", ea.len);
    
    err_check(    (ea.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL)
               || (u_prj.elem_p   == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check( (u_nodes.len != ea.nr_nodes), clean_up, 
               "Dimension mis-match (u_nodes: %d, nr nodes: %d)!",
               u_nodes.len, ea.nr_nodes);

    fem2d_te* dptr; /* domain pointer */
    
    matlib_complex* uptr0     = u_nodes.elem_p; /* base address */ 
    matlib_complex* uprj_ptr0 = u_prj.elem_p; /* base address */ 
    matlib_complex* uprj_ptr;
    matlib_complex* uptr;
    matlib_index offset1, offset2, offset3;
    matlib_real jacob;
    matlib_complex f1, f2, f3;

    /* Clear the existing values */
    for ( uprj_ptr = uprj_ptr0; 
          uprj_ptr < uprj_ptr0 + u_prj.len; uprj_ptr++)
    {
        *uprj_ptr = 0.0;
    }


    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); dptr++)
    {
        jacob = dptr->jacob;
        
        /* Get the first vertex */ 

        offset1 = dptr->nindex_p[FEM2D_INDEX_V1];
        offset2 = dptr->nindex_p[FEM2D_INDEX_V2];
        offset3 = dptr->nindex_p[FEM2D_INDEX_V3];

        uptr = uptr0 + offset1; f1 = *uptr * jacob;
        uptr = uptr0 + offset2; f2 = *uptr * jacob;
        uptr = uptr0 + offset3; f3 = *uptr * jacob;

        /* ===== */ 
        uprj_ptr = uprj_ptr0 + offset1;
        
        (*uprj_ptr) += (   MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V1] * f1
                         + MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V2] * f2
                         + MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V3] * f3);
        
        /* ===== */ 
        uprj_ptr = uprj_ptr0 + offset2;
        
        (*uprj_ptr) += (   MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V1] * f1
                         + MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V2] * f2
                         + MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V3] * f3);

        /* ===== */ 
        uprj_ptr = uprj_ptr0 + offset3;
        
        (*uprj_ptr) += (   MEMI[FEM2D_INDEX_V3][FEM2D_INDEX_V1] * f1
                         + MEMI[FEM2D_INDEX_V3][FEM2D_INDEX_V2] * f2
                         + MEMI[FEM2D_INDEX_V3][FEM2D_INDEX_V3] * f3);
    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;

} /* End of fem2d_NB_zprj */ 
/*============================================================================*/

matlib_real fem2d_NB_xnormL2
(
    fem2d_ea  ea,
    matlib_xv u_nodes /* values at the nodes */ 
)
{
    debug_enter( "nr of domains: %d", ea.len);
    
    err_check(    (ea.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != ea.nr_nodes), clean_up, 
               "Dimension mis-match (u_nodes: %d, v_nodes: %d, nr nodes: %d)!",
               u_nodes.len, ea.nr_nodes);

    fem2d_te* dptr; /* domain pointer */
    
    matlib_real* uptr0 = u_nodes.elem_p; /* base address */ 
    matlib_real* uptr;
    matlib_index offset1, offset2, offset3;
    matlib_real jacob;
    matlib_real f1, f2, f3;
    matlib_real f11, f12, f13, f23, f22, f33;

    matlib_real sum = 0;
    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); dptr++)
    {
        jacob = dptr->jacob;
        
        /* Get the first vertex */ 

        offset1 = dptr->nindex_p[FEM2D_INDEX_V1];
        offset2 = dptr->nindex_p[FEM2D_INDEX_V2];
        offset3 = dptr->nindex_p[FEM2D_INDEX_V3];

        uptr = uptr0 + offset1; f1 = *uptr;
        uptr = uptr0 + offset2; f2 = *uptr;
        uptr = uptr0 + offset3; f3 = *uptr;

        f11 = f1 * f1 * jacob;
        f22 = f2 * f2 * jacob;
        f33 = f3 * f3 * jacob;

        f12 = 2 * (f1 * f2) * jacob;
        f13 = 2 * (f1 * f3) * jacob;
        f23 = 2 * (f2 * f3) * jacob;
        
        sum += (   MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V1] * f11
                 + MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V2] * f22 
                 + MEMI[FEM2D_INDEX_V3][FEM2D_INDEX_V3] * f33 
                 + MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V2] * f12 
                 + MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V3] * f13 
                 + MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V3] * f23 );

    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return sqrt(sum);

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

} /* End of fem2d_NB_xnormL2 */ 


matlib_real fem2d_NB_znormL2
(
    fem2d_ea  ea,
    matlib_zv u_nodes /* values at the nodes */ 
)
{
    debug_enter( "nr of domains: %d", ea.len);
    
    err_check(    (ea.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != ea.nr_nodes), clean_up, 
               "Dimension mis-match (u_nodes: %d, v_nodes: %d, nr nodes: %d)!",
               u_nodes.len, ea.nr_nodes);

    fem2d_te* dptr; /* domain pointer */
    
    matlib_complex* uptr0 = u_nodes.elem_p; /* base address */ 
    matlib_complex* uptr;
    matlib_index offset1, offset2, offset3;
    matlib_real jacob;
    matlib_complex f1, f2, f3;
    matlib_real f11, f12, f13, f23, f22, f33;

    matlib_real sum = 0;
    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); dptr++)
    {
        jacob = dptr->jacob;
        
        /* Get the first vertex */ 

        offset1 = dptr->nindex_p[FEM2D_INDEX_V1];
        offset2 = dptr->nindex_p[FEM2D_INDEX_V2];
        offset3 = dptr->nindex_p[FEM2D_INDEX_V3];

        uptr = uptr0 + offset1; f1 = *uptr;
        uptr = uptr0 + offset2; f2 = *uptr;
        uptr = uptr0 + offset3; f3 = *uptr;

        f11 = (matlib_real) (f1 * conj(f1) * jacob);
        f22 = (matlib_real) (f2 * conj(f2) * jacob);
        f33 = (matlib_real) (f3 * conj(f3) * jacob);

        f12 = (matlib_real) (f1 * conj(f2) + f2 * conj(f1)) * jacob;
        f13 = (matlib_real) (f1 * conj(f3) + f3 * conj(f1)) * jacob;
        f23 = (matlib_real) (f2 * conj(f3) + f3 * conj(f2)) * jacob;
        
        sum += (   MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V1] * f11
                 + MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V2] * f22 
                 + MEMI[FEM2D_INDEX_V3][FEM2D_INDEX_V3] * f33 
                 + MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V2] * f12 
                 + MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V3] * f13 
                 + MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V3] * f23 );

    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return sqrt(sum);

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

} /* End of fem2d_NB_znormL2 */ 
/*============================================================================*/
matlib_real fem2d_NB_xiprod
(
    fem2d_ea  ea,
    matlib_xv u_nodes, /* values at the nodes */ 
    matlib_xv v_nodes
)
{
    debug_enter( "nr of domains: %d", ea.len);
    
    err_check(    (ea.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL)
               || (v_nodes.elem_p == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != ea.nr_nodes) 
               || (v_nodes.len != ea.nr_nodes), clean_up, 
               "Dimension mis-match (u_nodes: %d, v_nodes: %d, nr nodes: %d)!",
               u_nodes.len, v_nodes.len, ea.nr_nodes);

    fem2d_te* dptr; /* domain pointer */
    
    matlib_real* uptr0 = u_nodes.elem_p; /* base address */ 
    matlib_real* vptr0 = v_nodes.elem_p; /* base address */ 
    matlib_real* uptr;
    matlib_real* vptr;
    matlib_index offset1, offset2, offset3;
    matlib_real jacob;
    matlib_real f1, f2, f3;
    matlib_real g1, g2, g3;
    matlib_real fg11, fg22, fg33, fg12, fg13, fg23;

    matlib_real sum = 0;
    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); dptr++)
    {
        jacob = dptr->jacob;
        
        /* Get the first vertex */ 

        offset1 = dptr->nindex_p[FEM2D_INDEX_V1];
        offset2 = dptr->nindex_p[FEM2D_INDEX_V2];
        offset3 = dptr->nindex_p[FEM2D_INDEX_V3];

        uptr = uptr0 + offset1; f1 = *uptr;
        uptr = uptr0 + offset2; f2 = *uptr;
        uptr = uptr0 + offset3; f3 = *uptr;

        vptr = vptr0 + offset1; g1 = *vptr;
        vptr = vptr0 + offset2; g2 = *vptr;
        vptr = vptr0 + offset3; g3 = *vptr;

        fg11 = f1 * g1 * jacob;
        fg22 = f2 * g2 * jacob;
        fg33 = f3 * g3 * jacob;

        fg12 = (f1 * g2 + f2 * g1) * jacob;
        fg13 = (f1 * g3 + f3 * g1) * jacob;
        fg23 = (f2 * g3 + f3 * g2) * jacob;
        
        sum += (   MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V1] * fg11
                 + MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V2] * fg22 
                 + MEMI[FEM2D_INDEX_V3][FEM2D_INDEX_V3] * fg33 
                 + MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V2] * fg12 
                 + MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V3] * fg13 
                 + MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V3] * fg23 );
    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return sum;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

}/* End of fem2d_NB_xiprod */ 

matlib_complex fem2d_NB_ziprod
(
    fem2d_ea  ea,
    matlib_zv u_nodes, /* values at the nodes */ 
    matlib_zv v_nodes
)
{
    debug_enter( "nr of domains: %d", ea.len);
    
    err_check(    (ea.elem_p      == NULL)    
               || (u_nodes.elem_p == NULL)
               || (v_nodes.elem_p == NULL), clean_up, 
               "%s", "Null pointers ecountered!");

    err_check(    (u_nodes.len != ea.nr_nodes) 
               || (v_nodes.len != ea.nr_nodes), clean_up, 
               "Dimension mis-match (u_nodes: %d, v_nodes: %d, nr nodes: %d)!",
               u_nodes.len, v_nodes.len, ea.nr_nodes);

    fem2d_te* dptr; /* domain pointer */
    
    matlib_complex* uptr0 = u_nodes.elem_p; /* base address */ 
    matlib_complex* vptr0 = v_nodes.elem_p; /* base address */ 
    matlib_complex* uptr;
    matlib_complex* vptr;
    matlib_index offset1, offset2, offset3;
    matlib_real jacob;
    matlib_complex f1, f2, f3;
    matlib_complex g1, g2, g3;
    matlib_complex fg11, fg22, fg33, fg12, fg13, fg23;

    matlib_complex sum = 0;
    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); dptr++)
    {
        jacob = dptr->jacob;
        
        /* Get the first vertex */ 

        offset1 = dptr->nindex_p[FEM2D_INDEX_V1];
        offset2 = dptr->nindex_p[FEM2D_INDEX_V2];
        offset3 = dptr->nindex_p[FEM2D_INDEX_V3];

        uptr = uptr0 + offset1; f1 = *uptr;
        uptr = uptr0 + offset2; f2 = *uptr;
        uptr = uptr0 + offset3; f3 = *uptr;

        vptr = vptr0 + offset1; g1 = *vptr;
        vptr = vptr0 + offset2; g2 = *vptr;
        vptr = vptr0 + offset3; g3 = *vptr;

        fg11 = f1 * conj(g1) * jacob;
        fg22 = f2 * conj(g2) * jacob;
        fg33 = f3 * conj(g3) * jacob;

        fg12 = (f1 * conj(g2) + f2 * conj(g1)) * jacob;
        fg13 = (f1 * conj(g3) + f3 * conj(g1)) * jacob;
        fg23 = (f2 * conj(g3) + f3 * conj(g2)) * jacob;
        
        sum += (   MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V1] * fg11
                 + MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V2] * fg22 
                 + MEMI[FEM2D_INDEX_V3][FEM2D_INDEX_V3] * fg33 
                 + MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V2] * fg12 
                 + MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V3] * fg13 
                 + MEMI[FEM2D_INDEX_V2][FEM2D_INDEX_V3] * fg23 );
    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return sum;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;

}/* End of fem2d_NB_ziprod */ 

/*============================================================================*/
fem2d_err fem2d_quadM
(
    matlib_xm  vphi,
    matlib_xv  quadW,
    matlib_xm* Q
)
{
    debug_enter( "nr of quadrature points: %d", quadW.len);
    
    err_check(    (vphi.elem_p  == NULL)    
               || (quadW.elem_p == NULL), clean_up, 
               "%s", "Null pointer(s) ecountered!");

    err_check(    (vphi.lenc != quadW.len) 
               || (vphi.lenr != FEM2D_NV), clean_up, 
               "Dimension mis-match (vphi: %d-by-%d, quadW: %d)!",
               vphi.lenc, vphi.lenr, quadW.len);

    matlib_index mcnt = 0;
    fem2d_err error = matlib_create_xm( FEM2D_NR_COMBI, quadW.len, Q, 
                                        MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Memory allocation for quadrature matrix failed!");
    mcnt++;
    matlib_index i;
    matlib_index stride = quadW.len;

    for (i = 0; i < quadW.len; i++)
    {
        Q->elem_p[FEM2D_INDEX_V11 * stride + i] =   vphi.elem_p[i + FEM2D_INDEX_V1 * stride] 
                                                  * vphi.elem_p[i + FEM2D_INDEX_V1 * stride]
                                                  * quadW.elem_p[i];

        Q->elem_p[FEM2D_INDEX_V12 * stride + i] =   vphi.elem_p[i + FEM2D_INDEX_V1 * stride] 
                                                  * vphi.elem_p[i + FEM2D_INDEX_V2 * stride]
                                                  * quadW.elem_p[i];
        
        Q->elem_p[FEM2D_INDEX_V13 * stride + i] =   vphi.elem_p[i + FEM2D_INDEX_V1 * stride] 
                                                  * vphi.elem_p[i + FEM2D_INDEX_V3 * stride]
                                                  * quadW.elem_p[i];

        Q->elem_p[FEM2D_INDEX_V22 * stride + i] =   vphi.elem_p[i + FEM2D_INDEX_V2 * stride] 
                                                  * vphi.elem_p[i + FEM2D_INDEX_V2 * stride]
                                                  * quadW.elem_p[i];

        Q->elem_p[FEM2D_INDEX_V23 * stride + i] =   vphi.elem_p[i + FEM2D_INDEX_V2 * stride] 
                                                  * vphi.elem_p[i + FEM2D_INDEX_V3 * stride]
                                                  * quadW.elem_p[i];
    
        Q->elem_p[FEM2D_INDEX_V33 * stride + i] =   vphi.elem_p[i + FEM2D_INDEX_V3 * stride] 
                                                  * vphi.elem_p[i + FEM2D_INDEX_V3 * stride]
                                                  * quadW.elem_p[i];
    }

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;

clean_up:
    
    if (mcnt == 1)
        matlib_free(Q->elem_p);

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;
}
/*============================================================================*/

matlib_real fem2d_poly_symint
(
    matlib_index m,
    matlib_index n
)
/* 
 * value of: integral_{K_t} xi1^m * xi2^n dxi1 dxi2
 * */ 
{

    matlib_int minusone_raised2m = -1;
    matlib_int minusone_raised2n = -1;

    if((m % 2) == 0)
    {
        minusone_raised2m = 1;
    }

    if((n % 2) == 0)
    {
        minusone_raised2n = 1;
    }

    matlib_real integral = (   minusone_raised2m/((n + 1.0) * (m + n + 2.0)) 
                             + minusone_raised2n/((m + 1.0) * (m + n + 2.0))
                             + minusone_raised2m * minusone_raised2n/((m + 1.0) * (n + 1.0)));
    return integral;
}

matlib_real fem2d_check_quadM
(
    fem2d_cc     xi,
    matlib_xv    quadW,
    matlib_index m,
    matlib_index n
)
/* 
 * phi = xi1^m * xi2^n
 *
 * q1 ->   (xi1^m*xi1^2*xi2^n)/4 + (xi1^m*xi2^n*xi2^2)/4 + (xi1*xi1^m*xi2*xi2^n)/2 
 * q2 -> - (xi1*xi1^m*xi2^n)/4   - (xi1^m*xi2*xi2^n)/4   - (xi1^m*xi1^2*xi2^n)/4   - (xi1*xi1^m*xi2*xi2^n)/4
 * q3 -> - (xi1*xi1^m*xi2^n)/4   - (xi1^m*xi2*xi2^n)/4   - (xi1^m*xi2^n*xi2^2)/4   - (xi1*xi1^m*xi2*xi2^n)/4
 *
 * q4 ->   (xi1^m*xi2^n)/4       + (xi1*xi1^m*xi2^n)/2   + (xi1^m*xi1^2*xi2^n)/4 
 * q5 ->   (xi1^m*xi2^n)/4       + (xi1*xi1^m*xi2^n)/4   + (xi1^m*xi2*xi2^n)/4     + (xi1*xi1^m*xi2*xi2^n)/4
 *
 * q6 ->   (xi1^m*xi2^n)/4       + (xi1^m*xi2*xi2^n)/2   + (xi1^m*xi2^n*xi2^2)/4
 *
 * */ 
{
    debug_enter("m: %d, n: %d, nr of quad points: %d", m, n, xi.len);

    matlib_real q1, q2, q3, q4, q5, q6;

    q1 =   fem2d_poly_symint( m + 2, n + 0)/4.0
         + fem2d_poly_symint( m + 0, n + 2)/4.0
         + fem2d_poly_symint( m + 1, n + 1)/2.0;
    debug_body("q1: %0.16f", q1);

    q2 = - fem2d_poly_symint( m + 1, n + 0)/4.0
         - fem2d_poly_symint( m + 0, n + 1)/4.0
         - fem2d_poly_symint( m + 2, n + 0)/4.0
         - fem2d_poly_symint( m + 1, n + 1)/4.0;
    debug_body("q2: %0.16f", q2);

    q3 = - fem2d_poly_symint( m + 1, n + 0)/4.0
         - fem2d_poly_symint( m + 0, n + 1)/4.0
         - fem2d_poly_symint( m + 0, n + 2)/4.0
         - fem2d_poly_symint( m + 1, n + 1)/4.0;
    debug_body("q3: %0.16f", q3);

    q4 =   fem2d_poly_symint( m + 0, n + 0)/4.0
         + fem2d_poly_symint( m + 1, n + 0)/2.0
         + fem2d_poly_symint( m + 2, n + 0)/4.0;
    debug_body("q4: %0.16f", q4);

    q5 =   fem2d_poly_symint( m + 0, n + 0)/4.0
         + fem2d_poly_symint( m + 1, n + 0)/4.0
         + fem2d_poly_symint( m + 0, n + 1)/4.0
         + fem2d_poly_symint( m + 1, n + 1)/4.0;
    debug_body("q5: %0.16f", q5);

    q6 =   fem2d_poly_symint( m + 0, n + 0)/4.0
         + fem2d_poly_symint( m + 0, n + 1)/2.0
         + fem2d_poly_symint( m + 0, n + 2)/4.0;
    debug_body("q6: %0.16f", q6);

    matlib_index mcnt = 0;
    fem2d_err error;
    matlib_xv q, phi;
    error = matlib_create_xv( 6, &q, MATLIB_COL_VECT);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Memory allocation failed!");
    mcnt++; /* 1 */ 
    error = matlib_create_xv( xi.len, &phi, MATLIB_COL_VECT);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Memory allocation failed!");
    mcnt++; /* 2 */ 

    matlib_real *xiptr;
    matlib_real *phi_ptr = phi.elem_p;
    for ( xiptr = xi.elem_p; 
          xiptr < xi.elem_p + FEM2D_DIM * xi.len; 
          xiptr += FEM2D_DIM, phi_ptr++)
    {
        *phi_ptr = pow(xiptr[FEM2D_INDEX_DIM1], m) * pow(xiptr[FEM2D_INDEX_DIM2], n);
    }

    matlib_xm vphi, Q;
    error = fem2d_refbasis(xi, &vphi);
    //error = FEM2D_FAILURE;
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Failed to compute the values of reference basis functions!");

    mcnt++; /* 3 */

    error = fem2d_quadM(vphi, quadW, &Q);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Failed to compute the quadrature matrix!");
    mcnt++; /* 4 */ 

    error = matlib_xgemv( 1.0, Q, phi, 0.0, q);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Matrix product failed!");


    DEBUG_PRINT_XV(q, "%s:", "Gauss quadrature");
    matlib_real err = 0;
    err += ((q1 - q.elem_p[FEM2D_INDEX_V11]) * (q1 - q.elem_p[FEM2D_INDEX_V11]));
    err += ((q2 - q.elem_p[FEM2D_INDEX_V12]) * (q2 - q.elem_p[FEM2D_INDEX_V12]));
    err += ((q3 - q.elem_p[FEM2D_INDEX_V13]) * (q3 - q.elem_p[FEM2D_INDEX_V13]));
    err += ((q4 - q.elem_p[FEM2D_INDEX_V22]) * (q4 - q.elem_p[FEM2D_INDEX_V22]));
    err += ((q5 - q.elem_p[FEM2D_INDEX_V23]) * (q5 - q.elem_p[FEM2D_INDEX_V23]));
    err += ((q6 - q.elem_p[FEM2D_INDEX_V33]) * (q6 - q.elem_p[FEM2D_INDEX_V33]));

    matlib_real norm =   (q1 * q1) 
                       + (q2 * q2)
                       + (q3 * q3)
                       + (q4 * q4)
                       + (q5 * q5)
                       + (q6 * q6);

    matlib_real e_relative = sqrt(err/norm);
    debug_body("Relative error: %0.16f", e_relative);
    debug_exit("Exit Status: %s", "SUCCESS" );

    matlib_free(Q.elem_p);
    matlib_free(vphi.elem_p);
    matlib_free(phi.elem_p);
    matlib_free(q.elem_p);
    return e_relative;

clean_up:
    
    if (mcnt == 4)
    {
        matlib_free(Q.elem_p);
        mcnt--;
    }
    if (mcnt == 3)
    {
        matlib_free(vphi.elem_p);
        mcnt--;
    }
    if (mcnt == 2)
    {
        matlib_free(phi.elem_p);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(q.elem_p);
    }

    debug_body("Relative error: %0.16f", MATLIB_NAN);
    debug_exit("Exit Status: %s", "FAILURE" );
    return MATLIB_NAN;
}

/*============================================================================*/
matlib_index fem2d_get_nnz(fem2d_ea ea)
{
    debug_enter( "nr nodes: %d, nr domains: %d",
                 ea.nr_nodes, ea.len);

    err_check( ea.vpatch_p == NULL, clean_up, 
               "%s", "Null pointer encountered!");
    matlib_index i, j;
    matlib_index nnz = 0;
    for (i = 0; i < ea.nr_nodes; i++)
    {
        nnz++;
        for (j = 0; j < ea.vpatch_p[i].len; j++)
        {
            if (ea.vpatch_p[i].node_order[j]>0)
                nnz++;
        }

        if (ea.vpatch_p[i].point_enum == FEM2D_BOUNDARY)
        {
            if (ea.vpatch_p[i].node_order[j]>0)
                nnz++;
        }
    }
    debug_exit("Exit Status: %s", "SUCCESS" );
    return nnz;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return 0;
}

/*============================================================================*/

fem2d_err fem2d_GMMSparsity
(
    fem2d_ea      ea,
    matlib_index* row, /* required to be pre-allocated */                     
    matlib_index* col  /* required to be pre-allocated */                  
)
{
    debug_enter( "nr nodes: %d, nr domains: %d", ea.nr_nodes, ea.len);

    err_check(    (ea.vpatch_p == NULL) 
               || (row         == NULL) 
               || (col         == NULL),  clean_up,
               "%s", "Null pointer(s) encountered!");

    /* length of col: nnz
     * length of row: nr_nodes + 1
     *
     * */ 
    matlib_index s0, s = 0;
    matlib_index i = 0, j;
    matlib_int node_order;

    fem2d_te* dptr;
    matlib_index next_vertex, other_vertex1, other_vertex2;

    for (i = 0; i < ea.nr_nodes; i++)
    {
        row[i] = s;
        s0 = s;

        col[s0] = i; 
        s++;
        for (j = 0; j < ea.vpatch_p[i].len; j++)
        {
            node_order = ea.vpatch_p[i].node_order[j];
            if (node_order > 0)
            {
                dptr = ea.vpatch_p[i].domain_p[j];
                next_vertex = ea.vpatch_p[i].bvert_index[j];
                col[s0 + node_order] = dptr->nindex_p[next_vertex]; 
                s++;
                debug_body( "row: %d, node order: %d, vertex: %d, boundary: %d",
                            i, node_order, 
                            dptr->nindex_p[ea.vpatch_p[i].vert_index[j]], 
                            dptr->nindex_p[next_vertex]);
            }
        }
        if(ea.vpatch_p[i].point_enum == FEM2D_BOUNDARY)
        {
            node_order = ea.vpatch_p[i].node_order[j];
            if(node_order > 0)
            {
                j = ea.vpatch_p[i].len - 1;

                /* This missing statement caused a lot of trouble, The domain
                 * pointer, dptr didmnot get update when the body of the code in
                 * the loop above did not get executed because of the "if"
                 * condition. 
                 * */ 
                dptr = ea.vpatch_p[i].domain_p[j]; 

                next_vertex = FEM2D_NV - ea.vpatch_p[i].vert_index[j] 
                                       - ea.vpatch_p[i].bvert_index[j];

                col[s0 + node_order] = dptr->nindex_p[next_vertex];
                s++;
                debug_body( "row: %d, node order: %d, vertex: %d, boundary: %d",
                            i, node_order, 
                            dptr->nindex_p[ea.vpatch_p[i].vert_index[j]], 
                            dptr->nindex_p[next_vertex]);
            }
        }
    }
    debug_body( "nr non-zero entries: %d", s);
    row[i] = s;

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;
}

fem2d_err fem2d_XCSRGMM
(
    fem2d_ea      ea,
    matlib_xv     q,
    matlib_index* row,                     
    matlib_index* col,                     
    matlib_real*  ugpmm                   
)
{
    debug_enter("nr nodes: %d, q: %d", ea.len, q.len);

    err_check(    (ea.elem_p   == NULL) 
               || (ea.vpatch_p == NULL)
               || (q.elem_p    == NULL), clean_up, 
               "%s", "Null pointer(s) encountered!");

    err_check( ea.len * FEM2D_NR_COMBI != q.len, clean_up,
               "%s: Dimension mistmatch (q: %d, required size: %d)!",
               q.len, ea.len * FEM2D_NR_COMBI);

    fem2d_te* dptr;
    matlib_real* qptr = q.elem_p;
    
    matlib_index i, pos_vpatch, other_vertex, node_index;
    matlib_int j;

    for (i = 0; i < row[ea.nr_nodes]; i++ )
    {
        ugpmm[i] = 0.0;
    }
    for ( dptr = ea.elem_p; dptr < ea.elem_p + ea.len;
          dptr++, qptr += FEM2D_NR_COMBI)
    {
        debug_body("Domain: %d", dptr->domain_index);

        node_index = dptr->nindex_p[FEM2D_INDEX_V1];

        i = row[node_index];
        ugpmm[i] += qptr[FEM2D_INDEX_V11];
        
        pos_vpatch = dptr->pos_vpatch[FEM2D_INDEX_V1];

        j = ea.vpatch_p[node_index].node_order[pos_vpatch];
        other_vertex = ea.vpatch_p[node_index].bvert_index[pos_vpatch];

        if (other_vertex == FEM2D_INDEX_V2)
        {
            if (j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V12];
            }

            if (ea.vpatch_p[node_index].point_enum == FEM2D_INTERIOR)
            {
                if (pos_vpatch < (ea.vpatch_p[node_index].len - 1))
                {
                    j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
                }
                else
                {
                    j = ea.vpatch_p[node_index].node_order[0];
                }
            }
            else
            {
                j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
            }
            if(j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V13];
            }
        }
        else
        {
            if (j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V13];
            }

            if (ea.vpatch_p[node_index].point_enum == FEM2D_INTERIOR)
            {
                if (pos_vpatch < (ea.vpatch_p[node_index].len - 1))
                {
                    j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
                }
                else
                {
                    j = ea.vpatch_p[node_index].node_order[0];
                }
            }
            else
            {
                j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
            }
            if(j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V12];
            }
        }
        /* ======== */
        
        node_index = dptr->nindex_p[FEM2D_INDEX_V2];

        i = row[node_index];
        ugpmm[i] += qptr[FEM2D_INDEX_V22];
        
        pos_vpatch = dptr->pos_vpatch[FEM2D_INDEX_V2];
        debug_body( "node: %d, pos in vertex patch: %d",
                    node_index, pos_vpatch);

        j = ea.vpatch_p[node_index].node_order[pos_vpatch];
        other_vertex = ea.vpatch_p[node_index].bvert_index[pos_vpatch];

        if (other_vertex == FEM2D_INDEX_V3)
        {
            if (j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V23];
            }

            if (ea.vpatch_p[node_index].point_enum == FEM2D_INTERIOR)
            {
                if (pos_vpatch < (ea.vpatch_p[node_index].len - 1))
                {
                    j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
                }
                else
                {
                    j = ea.vpatch_p[node_index].node_order[0];
                }
            }
            else
            {
                j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
            }
            if(j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V12];
            }
        }
        else
        {
            if (j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V12];
            }

            if (ea.vpatch_p[node_index].point_enum == FEM2D_INTERIOR)
            {
                if (pos_vpatch < (ea.vpatch_p[node_index].len - 1))
                {
                    j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
                }
                else
                {
                    j = ea.vpatch_p[node_index].node_order[0];
                }
            }
            else
            {
                j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
            }
            if(j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V23];
            }
        }

        /* ======== */
        
        node_index = dptr->nindex_p[FEM2D_INDEX_V3];

        i = row[node_index];
        ugpmm[i] += qptr[FEM2D_INDEX_V33];
        
        pos_vpatch = dptr->pos_vpatch[FEM2D_INDEX_V3];

        j = ea.vpatch_p[node_index].node_order[pos_vpatch];
        other_vertex = ea.vpatch_p[node_index].bvert_index[pos_vpatch];

        if (other_vertex == FEM2D_INDEX_V1)
        {
            if (j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V13];
            }

            if (ea.vpatch_p[node_index].point_enum == FEM2D_INTERIOR)
            {
                if (pos_vpatch < (ea.vpatch_p[node_index].len - 1))
                {
                    j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
                }
                else
                {
                    j = ea.vpatch_p[node_index].node_order[0];
                }
            }
            else
            {
                j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
            }
            if(j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V23];
            }
        }
        else
        {
            if (j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V23];
            }

            if (ea.vpatch_p[node_index].point_enum == FEM2D_INTERIOR)
            {
                if (pos_vpatch < (ea.vpatch_p[node_index].len - 1))
                {
                    j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
                }
                else
                {
                    j = ea.vpatch_p[node_index].node_order[0];
                }
            }
            else
            {
                j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
            }
            if(j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V13];
            }
        }
    }

    /* === */ 

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;
}


fem2d_err fem2d_ZCSRGMM
(
    fem2d_ea        ea,
    matlib_zv       q,
    matlib_index*   row,                     
    matlib_index*   col,                     
    matlib_complex* ugpmm                   
)
{
    debug_enter("nr nodes: %d, q: %d", ea.len, q.len);

    err_check(    (ea.elem_p   == NULL) 
               || (ea.vpatch_p == NULL)
               || (q.elem_p    == NULL), clean_up, 
               "%s", "Null pointer(s) encountered!");

    err_check( ea.len * FEM2D_NR_COMBI != q.len, clean_up,
               "%s: Dimension mistmatch (q: %d, required size: %d)!",
               q.len, ea.len * FEM2D_NR_COMBI);

    fem2d_te* dptr;
    matlib_complex* qptr = q.elem_p;
    
    matlib_index i, pos_vpatch, other_vertex, node_index;
    matlib_int j;

    for (i = 0; i < row[ea.nr_nodes]; i++ )
    {
        ugpmm[i] = 0.0;
    }
    for ( dptr = ea.elem_p; dptr < ea.elem_p + ea.len;
          dptr++, qptr += FEM2D_NR_COMBI)
    {
        debug_body("Domain: %d", dptr->domain_index);

        node_index = dptr->nindex_p[FEM2D_INDEX_V1];

        i = row[node_index];
        ugpmm[i] += qptr[FEM2D_INDEX_V11];
        
        pos_vpatch = dptr->pos_vpatch[FEM2D_INDEX_V1];

        j = ea.vpatch_p[node_index].node_order[pos_vpatch];
        other_vertex = ea.vpatch_p[node_index].bvert_index[pos_vpatch];

        if (other_vertex == FEM2D_INDEX_V2)
        {
            if (j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V12];
            }

            if (ea.vpatch_p[node_index].point_enum == FEM2D_INTERIOR)
            {
                if (pos_vpatch < (ea.vpatch_p[node_index].len - 1))
                {
                    j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
                }
                else
                {
                    j = ea.vpatch_p[node_index].node_order[0];
                }
            }
            else
            {
                j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
            }
            if(j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V13];
            }
        }
        else
        {
            if (j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V13];
            }

            if (ea.vpatch_p[node_index].point_enum == FEM2D_INTERIOR)
            {
                if (pos_vpatch < (ea.vpatch_p[node_index].len - 1))
                {
                    j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
                }
                else
                {
                    j = ea.vpatch_p[node_index].node_order[0];
                }
            }
            else
            {
                j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
            }
            if(j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V12];
            }
        }
        /* ======== */
        
        node_index = dptr->nindex_p[FEM2D_INDEX_V2];

        i = row[node_index];
        ugpmm[i] += qptr[FEM2D_INDEX_V22];
        
        pos_vpatch = dptr->pos_vpatch[FEM2D_INDEX_V2];
        debug_body( "node: %d, pos in vertex patch: %d",
                    node_index, pos_vpatch);

        j = ea.vpatch_p[node_index].node_order[pos_vpatch];
        other_vertex = ea.vpatch_p[node_index].bvert_index[pos_vpatch];

        if (other_vertex == FEM2D_INDEX_V3)
        {
            if (j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V23];
            }

            if (ea.vpatch_p[node_index].point_enum == FEM2D_INTERIOR)
            {
                if (pos_vpatch < (ea.vpatch_p[node_index].len - 1))
                {
                    j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
                }
                else
                {
                    j = ea.vpatch_p[node_index].node_order[0];
                }
            }
            else
            {
                j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
            }
            if(j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V12];
            }
        }
        else
        {
            if (j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V12];
            }

            if (ea.vpatch_p[node_index].point_enum == FEM2D_INTERIOR)
            {
                if (pos_vpatch < (ea.vpatch_p[node_index].len - 1))
                {
                    j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
                }
                else
                {
                    j = ea.vpatch_p[node_index].node_order[0];
                }
            }
            else
            {
                j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
            }
            if(j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V23];
            }
        }

        /* ======== */
        
        node_index = dptr->nindex_p[FEM2D_INDEX_V3];

        i = row[node_index];
        ugpmm[i] += qptr[FEM2D_INDEX_V33];
        
        pos_vpatch = dptr->pos_vpatch[FEM2D_INDEX_V3];

        j = ea.vpatch_p[node_index].node_order[pos_vpatch];
        other_vertex = ea.vpatch_p[node_index].bvert_index[pos_vpatch];

        if (other_vertex == FEM2D_INDEX_V1)
        {
            if (j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V13];
            }

            if (ea.vpatch_p[node_index].point_enum == FEM2D_INTERIOR)
            {
                if (pos_vpatch < (ea.vpatch_p[node_index].len - 1))
                {
                    j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
                }
                else
                {
                    j = ea.vpatch_p[node_index].node_order[0];
                }
            }
            else
            {
                j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
            }
            if(j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V23];
            }
        }
        else
        {
            if (j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V23];
            }

            if (ea.vpatch_p[node_index].point_enum == FEM2D_INTERIOR)
            {
                if (pos_vpatch < (ea.vpatch_p[node_index].len - 1))
                {
                    j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
                }
                else
                {
                    j = ea.vpatch_p[node_index].node_order[0];
                }
            }
            else
            {
                j = ea.vpatch_p[node_index].node_order[pos_vpatch + 1];
            }
            if(j>0)
            {
                ugpmm[i + j] += qptr[FEM2D_INDEX_V13];
            }
        }
    }

    /* === */ 

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;
}


fem2d_err fem2d_XCSRGMM1
(
    fem2d_ea      ea,
    matlib_index* row,                     
    matlib_index* col,                     
    matlib_real*  ugpmm                   
)
{
    debug_enter("nr nodes: %d", ea.len);

    err_check(    (ea.elem_p   == NULL) 
               || (ea.vpatch_p == NULL), clean_up,
            "%s", "Null pointer(s) encountered!");

    matlib_index i, j, k;
    matlib_int node_order;
    matlib_real jacob;

    for (i = 0; i < ea.nr_nodes; i++)
    {
        k = row[i];
        for (j = k; j < row[i+1]; j++)
        {
            ugpmm[j] = 0;
        }
        for (j = 0; j < (ea.vpatch_p[i].len - 1); j++)
        {
            jacob = ea.vpatch_p[i].domain_p[j]->jacob;
            ugpmm[k] += (jacob * MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V1]);

            node_order = ea.vpatch_p[i].node_order[j];
            if (node_order>0)
            {
                ugpmm[k + node_order] += (jacob * MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V2]);
            }
            node_order = ea.vpatch_p[i].node_order[j+1];
            if (node_order>0)
            {
                ugpmm[k + node_order] += (jacob * MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V2]);
            }
        }
        jacob = ea.vpatch_p[i].domain_p[j]->jacob;
        ugpmm[k] += (jacob * MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V1]);

        node_order = ea.vpatch_p[i].node_order[j];
        if(node_order > 0)
        {
            ugpmm[k + node_order ] += (jacob * MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V2]);
        }
        if (ea.vpatch_p[i].point_enum == FEM2D_INTERIOR)
        {
            node_order = ea.vpatch_p[i].node_order[0];
            if (node_order > 0)
            {
                ugpmm[k + node_order] += (jacob * MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V2]);
            }
        }
        else
        {
            node_order = ea.vpatch_p[i].node_order[j+1];
            if (node_order > 0)
            {
                ugpmm[k + node_order] += (jacob * MEMI[FEM2D_INDEX_V1][FEM2D_INDEX_V2]);
            }
        }
    }
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;


}

/*============================================================================*/
fem2d_err fem2d_xm_sparse_GMM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem2d_ea          ea,
    matlib_xm         Q,
    matlib_xv         phi,
    matlib_xm_sparse* M
)
{
    debug_enter( "size of Q: %d-by-%d), "
                 "length of phi: %d",
                 Q.lenc, Q.lenr, phi.len);
    
    err_check( Q.lenc != FEM2D_NR_COMBI, clean_up,
               "Dimension mismatch (Q: %d-by-%d, nr combi: %d)!", 
               Q.lenc, Q.lenr, FEM2D_NR_COMBI);

    err_check( (Q.lenr * ea.len != phi.len), clean_up,
               "Dimension mismatch (phi: %d, require size: %d)!", 
               phi.len, Q.lenr * ea.len);

    matlib_index nnz = fem2d_get_nnz(ea);
    debug_body( "nr. of non-zero elements: %d", nnz);
    err_check( nnz == 0, clean_up, "%s", 
               "Number of non-zero elements of the mass-matrix is inconsistent!");

    matlib_index i, mcnt = 0;
    fem2d_err error;
    
    matlib_xv q;
    error = matlib_create_xv( Q.lenc * ea.len, &q, MATLIB_COL_VECT);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Failed to allocate memory for the inner product array!");
    mcnt++; /* 1 */ 
    matlib_real alpha, beta = 0.0;
    phi.len = Q.lenr;

    matlib_xv q_tmp = { .len = Q.lenc, 
                        .type = MATLIB_COL_VECT,
                        .elem_p = q.elem_p};

    for ( i = 0; i < ea.len; 
          i++, (phi.elem_p) += (Q.lenr), (q_tmp.elem_p) += (Q.lenc))
    {
        alpha = ea.elem_p[i].jacob;
        error = matlib_xgemv(alpha, Q, phi, beta, q_tmp );
        err_check( (error == FEM2D_FAILURE), clean_up, 
                   "%s", "Failed to compute inner-products!");
    }
    
    M->lenc   = ea.nr_nodes;
    M->lenr   = ea.nr_nodes;
    
    M->rowIn  = calloc( ea.nr_nodes + 1, sizeof(matlib_index));
    err_check( (M->rowIn == NULL), clean_up, 
               "%s", "Failed to allocate memory for the row array!");
    mcnt++; /* 2 */ 

    M->colIn  = calloc(   nnz, sizeof(matlib_index));
    err_check( (M->colIn == NULL), clean_up,
               "%s", "Failed to allocate memory for the column array!");
    mcnt++; /* 3 */ 
    
    M->elem_p = calloc(   nnz, sizeof(matlib_real));
    err_check( (M->elem_p == NULL), clean_up,
               "%s", "Failed to allocate memory for the sparse matrix entries!");
    mcnt++; /* 4 */ 


    error = fem2d_GMMSparsity( ea, M->rowIn, M->colIn);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Failed to get the sparsity structure of the mass-matrix!");

    error = fem2d_XCSRGMM( ea, q, M->rowIn, M->colIn, M->elem_p);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Failed to compute the mass-matrix!");

    
    matlib_free(q.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;

clean_up:
    if (mcnt == 4)
    {
        matlib_free(M->elem_p);
        mcnt--;
    }
    if (mcnt == 3)
    {
        matlib_free(M->colIn);
        mcnt--;
    }
    if (mcnt == 2)
    {
        matlib_free(M->rowIn);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(q.elem_p);
    }

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;

}

fem2d_err fem2d_zm_sparse_GMM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem2d_ea          ea,
    matlib_xm         Q,
    matlib_zv         phi,
    matlib_zm_sparse* M
)
{
    debug_enter( "size of Q: %d-by-%d), "
                 "length of phi: %d",
                 Q.lenc, Q.lenr, phi.len);
    
    err_check( Q.lenc != FEM2D_NR_COMBI, clean_up,
               "Dimension mismatch (Q: %d-by-%d, nr combi: %d)!", 
               Q.lenc, Q.lenr, FEM2D_NR_COMBI);

    err_check( (Q.lenr * ea.len != phi.len), clean_up,
               "Dimension mismatch (phi: %d, require size: %d)!", 
               phi.len, Q.lenr * ea.len);

    matlib_index nnz = fem2d_get_nnz(ea);
    debug_body( "nr. of non-zero elements: %d", nnz);
    err_check( nnz == 0, clean_up, "%s", 
               "Number of non-zero elements of the mass-matrix is inconsistent!");

    matlib_index i, mcnt = 0;
    fem2d_err error;
    
    matlib_zv q;
    error = matlib_create_zv( Q.lenc * ea.len, &q, MATLIB_COL_VECT);
    err_check( (error == FEM2D_FAILURE), clean_up, 
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

    for ( i = 0; i < ea.len; 
          i++, (phi_tmp.elem_p) += (COMPLEX_DIM * Q.lenr), (q_tmp.elem_p) += (COMPLEX_DIM * Q.lenc))
    {
        alpha = ea.elem_p[i].jacob;
        error = matlib_xgemm(alpha, Q, phi_tmp, beta, q_tmp );
        err_check( (error == FEM2D_FAILURE), clean_up, 
                   "%s", "Failed to compute inner-products!");
    }
    
    M->lenc   = ea.nr_nodes;
    M->lenr   = ea.nr_nodes;
    
    M->rowIn  = calloc( ea.nr_nodes + 1, sizeof(matlib_index));
    err_check( (M->rowIn == NULL), clean_up, 
               "%s", "Failed to allocate memory for the row array!");
    mcnt++; /* 2 */ 

    M->colIn  = calloc(   nnz, sizeof(matlib_index));
    err_check( (M->colIn == NULL), clean_up,
               "%s", "Failed to allocate memory for the column array!");
    mcnt++; /* 3 */ 
    
    M->elem_p = calloc(   nnz, sizeof(matlib_complex));
    err_check( (M->elem_p == NULL), clean_up,
               "%s", "Failed to allocate memory for the sparse matrix entries!");
    mcnt++; /* 4 */ 


    error = fem2d_GMMSparsity( ea, M->rowIn, M->colIn);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Failed to get the sparsity structure of the mass-matrix!");

    error = fem2d_ZCSRGMM( ea, q, M->rowIn, M->colIn, M->elem_p);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Failed to compute the mass-matrix!");

    
    matlib_free(q.elem_p);
    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;

clean_up:
    if (mcnt == 4)
    {
        matlib_free(M->elem_p);
        mcnt--;
    }
    if (mcnt == 3)
    {
        matlib_free(M->colIn);
        mcnt--;
    }
    if (mcnt == 2)
    {
        matlib_free(M->rowIn);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(q.elem_p);
    }

    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;

}

fem2d_err fem2d_xm_sparse_GMM1
/* Real - Assemble Global Mass Matrix*/ 
(
    fem2d_ea          ea,
    matlib_xm_sparse* M
)
{
    debug_enter("", "%s");

    matlib_index nnz = fem2d_get_nnz(ea);
    debug_body( "nr. of non-zero elements: %d", nnz);
    err_check( nnz == 0, clean_up, "%s", 
               "Number of non-zero elements of the mass-matrix is inconsistent!");

    matlib_index i, mcnt = 0;
    fem2d_err error;
    
    
    M->lenc   = ea.nr_nodes;
    M->lenr   = ea.nr_nodes;
    
    M->rowIn  = calloc( ea.nr_nodes + 1, sizeof(matlib_index));
    err_check( (M->rowIn == NULL), clean_up, 
               "%s", "Failed to allocate memory for the row array!");
    mcnt++; /* 1 */ 

    M->colIn  = calloc( nnz, sizeof(matlib_index));
    err_check( (M->colIn == NULL), clean_up,
               "%s", "Failed to allocate memory for the column array!");
    mcnt++; /* 2 */ 
    
    M->elem_p = calloc( nnz, sizeof(matlib_real));
    err_check( (M->elem_p == NULL), clean_up,
               "%s", "Failed to allocate memory for the sparse matrix entries!");
    mcnt++; /* 3 */ 


    error = fem2d_GMMSparsity( ea, M->rowIn, M->colIn);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Failed to get the sparsity structure of the mass-matrix!");

    error = fem2d_XCSRGMM1( ea, M->rowIn, M->colIn, M->elem_p);
    err_check( (error == FEM2D_FAILURE), clean_up, 
               "%s", "Failed to compute the mass-matrix!");

    debug_exit("Exit Status: %s", "SUCCESS" );
    return FEM2D_SUCCESS;

clean_up:
    if (mcnt == 3)
    {
        matlib_free(M->elem_p);
        mcnt--;
    }
    if (mcnt == 2)
    {
        matlib_free(M->colIn);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(M->rowIn);
        mcnt--;
    }
    debug_exit("Exit Status: %s", "FAILURE" );
    return FEM2D_FAILURE;
}
/*============================================================================*/
fem2d_err fem2d_getmesh
(
    char* file_name,
    fem2d_cc* nodes,
    matlib_nv* iv
)
{
    debug_enter("Reading mesh data from: %s", file_name);

    matlib_io_t mp;
    matlib_index mcnt;
    fem2d_err error = matlib_io_fread(&mp, file_name);
    err_check( error == FEM2D_FAILURE, clean_up, 
               "Failed to read the file (name: %s)!", file_name);
    mcnt++;

    matlib_index pos_node, pos_ia;
    bool is_data_OK = false;
    if (mp.format[0] == MATLIB_DEN_XM) 
    {
        pos_node = 0;
        pos_ia   = 1;
        is_data_OK = true;

    }
    else if (mp.format[0] == MATLIB_NV) 
    {
        pos_node = 1;
        pos_ia   = 0;
        is_data_OK = true;
    }
    err_check( !is_data_OK, clean_up, 
               "%s", "Data is not in the expected format!");

    matlib_xm* nodes_mat = (matlib_xm*)mp.data_p[pos_node];
    if (nodes_mat->order == MATLIB_ROW_MAJOR)
    {
        err_check(    (nodes_mat->lenr != FEM2D_DIM) 
                   || (nodes_mat->lenc < FEM2D_NV) , clean_up,
                   "Mesh data is inconsistent(storage order: %s)!",
                   MATLIB_ORDER_ENUM2STR(nodes_mat->order));
        nodes->len = nodes_mat->lenc;
    }
    else if (nodes_mat->order == MATLIB_COL_MAJOR)
    {
        err_check(    (nodes_mat->lenc != FEM2D_DIM) 
                   || (nodes_mat->lenr < FEM2D_NV ) , clean_up,
                   "Mesh data is inconsistent(storage order: %s)!",
                   MATLIB_ORDER_ENUM2STR(nodes_mat->order));
        nodes->len = nodes_mat->lenr;
    }
    nodes->elem_p = nodes_mat->elem_p;

    matlib_nv* tmp = (matlib_nv*)mp.data_p[pos_ia];
    iv->len    = tmp->len;
    iv->elem_p = tmp->elem_p;

    matlib_io_free(&mp);
    debug_exit("Exit Status: %s", "SUCCESS");
    return FEM2D_SUCCESS;

clean_up:
    if (mcnt == 1)
    {
        matlib_io_freeall(&mp);
    }
    nodes->len    = 0;
    nodes->elem_p = NULL;

    iv->elem_p = NULL;
    iv->len    = 0;

    debug_exit("Exit Status: %s", "FAILURE");
    return FEM2D_FAILURE;
}

/*============================================================================*/

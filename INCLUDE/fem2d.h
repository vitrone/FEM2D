#ifndef FEM1D_H
#define FEM1D_H

/*============================================================================+/
 | Include all the dependencies
/+============================================================================*/
#include "basic.h"
#include "matlib.h"
#include "debug.h"
#include "ehandler.h"

/*============================================================================+/
 |DATA STRUCTURES AND ENUMS
/+============================================================================*/
typedef matlib_err fem2d_err;

#define FEM2D_SUCCESS MATLIB_SUCCESS
#define FEM2D_FAILURE MATLIB_FAILURE
#define FEM2D_ERR_ENUM2STR MATLIB_ERR_ENUM2STR


/* Point in R^2 (x_1, x_2)
 *
 * */ 

typedef struct
{
    matlib_index len;
    matlib_real* elem_p;

} fem2d_cc; /* cartesian coordinates in R^2*/ 

typedef enum
{
    FEM2D_INTERIOR,
    FEM2D_BOUNDARY

} FEM2D_POINT_T;

#define FEM2D_POINT_T_ENUM2STR(point_enum)                          \
    (FEM2D_INTERIOR == point_enum ? "INTERIOR POINT":               \
     (FEM2D_BOUNDARY == point_enum ? "BOUNDARY POINT": "UNKNOWN"))


typedef struct
{
    matlib_real** vert_p; /* pointer to the vertices */ 
    matlib_index* nindex_p; /* pointer to the node indices forming the triangle */ 
    matlib_index* vorder_p; 
    matlib_index  domain_index;
    matlib_index* pos_vpatch;
    matlib_real   jacob;
    matlib_real*  ijmat;
    matlib_real*  jmat;
    matlib_int    sign;

} fem2d_te; /* triangular element */

/* Define the vertex patch for a given node: TBD
 *
 *
 * */ 
typedef struct
{
    matlib_index node_index;
    matlib_index len; /* nr. of domains in the vertex patch */ 
    FEM2D_POINT_T point_enum;
    
    /* array of domain pointers arranged such that neighbhouring 
     * domains are next to each other in the list 
     * */ 
    fem2d_te**    domain_p;
    matlib_index* vert_index; /* array of vertex order indices */ 
    matlib_int*   node_order;
    /* index of vertex on the boundary not shared by the next neighbhoiring
     * domain in the vertex patch
     * */ 
    matlib_index* bvert_index; 

} fem2d_vp; /* vertex patch */ 

typedef struct
{
    matlib_real* node_p;    /* node base address */ 
    matlib_index len;      /* nr of domains */ 
    matlib_index nr_nodes; /* nr of nodes */ 
    fem2d_te*    elem_p;
    fem2d_vp*    vpatch_p;

} fem2d_ea; /* element array type */ 


/* NZE : Non-Zero Elements */ 
typedef enum
{
    FEM2D_GSMM_INIT,
    FEM2D_GET_SPARSITY_ONLY,
    FEM2D_GET_NZE_ONLY,
    FEM2D_GET_SPARSITY_NZE,
    FEM2D_GSMM_FREE,

} FEM2D_OPT_GSMM; /* OPTIONS GSMM */ 



/* Cartesian coordinate indices */ 
extern const matlib_index FEM2D_INDEX_DIM1;
extern const matlib_index FEM2D_INDEX_DIM2;
extern const matlib_index FEM2D_DIM;

/* Jacobian indices */ 
extern const matlib_index FEM2D_INDEX_J11;
extern const matlib_index FEM2D_INDEX_J12;
extern const matlib_index FEM2D_INDEX_J21;
extern const matlib_index FEM2D_INDEX_J22;

/* Ref. triangle vertices */
extern const matlib_real FEM2D_VERT[6];

/* Vertex indices */ 
extern const matlib_index FEM2D_INDEX_V1;
extern const matlib_index FEM2D_INDEX_V2;
extern const matlib_index FEM2D_INDEX_V3;
extern const matlib_index FEM2D_NV;
extern const matlib_index FEM2D_INIT_VPATCH_SIZE;
extern const matlib_real MEMI[3][3];

extern const matlib_index FEM2D_INDEX_V11;
extern const matlib_index FEM2D_INDEX_V12;
extern const matlib_index FEM2D_INDEX_V13;
extern const matlib_index FEM2D_INDEX_V22;
extern const matlib_index FEM2D_INDEX_V23;
extern const matlib_index FEM2D_INDEX_V33;
extern const matlib_index FEM2D_NR_COMBI;

/*============================================================================*/
#define FEM2D_MESH_TOL (1e-9)

typedef enum
{
    FEM2D_TEST_PASSED,
    FEM2D_TEST_FAILED

} FEM2D_TEST_PF;

#define FEM2D_TEST_PF_ENUM2STR(pf_enum)                   \
    (pf_enum == FEM2D_TEST_PASSED ? "PASSED":             \
     (pf_enum == FEM2D_TEST_FAILED ? "FAILED": "UNKLNOWN"))


typedef struct
{
    matlib_index nr_total;
    matlib_index nr_failed;
    matlib_index nr_passed;

    FEM2D_TEST_PF  option;
    FEM2D_TEST_PF* edge_mid_points;
    FEM2D_TEST_PF* pos_vpatch;
    FEM2D_TEST_PF* vpatch_boundary;
    matlib_index*  vpatch_index; 

} fem2d_vpinfo_t;

#define BOOL2PF(b)                                \
    ((b) ? FEM2D_TEST_PASSED : FEM2D_TEST_FAILED)

#define VPINFO_PRINT(fp, fmt,...)   \
        do {                        \
            fprintf( fp,            \
                     "\t" fmt "\n", \
                     __VA_ARGS__);  \
           } while (0)

/*============================================================================*/

fem2d_err fem2d_create_cc
(
    const matlib_index length,
          fem2d_cc*    nodes
);

fem2d_err fem2d_refbasis
(
    const fem2d_cc   xi,
          matlib_xm* vphi
);


fem2d_err fem2d_create_ea
(
    const fem2d_cc      nodes,
    const matlib_index* ia,
    const matlib_index  nr_domains,
          fem2d_ea*     ea
);

fem2d_err fem2d_create_vp(fem2d_ea *ea);

fem2d_err fem2d_sort_node_index
(
    matlib_index  zeroth,
    matlib_index* ilist,
    matlib_index  ilen,
    matlib_int*   iorder
);

fem2d_err fem2d_check_vp
(
    fem2d_ea ea,
    fem2d_vpinfo_t* test_info
);

fem2d_err fem2d_write_vpinfo
(
    fem2d_vpinfo_t* test_info,
    fem2d_ea ea,
    char* file_name
);

void fem2d_free_vpinfo(fem2d_vpinfo_t* test_info);

void fem2d_free_ea(fem2d_ea ea);

fem2d_err fem2d_create_ia
(
    const fem2d_ea      ea,
          matlib_index* ia
);

fem2d_err fem2d_centroid
(
    const fem2d_ea ea,
          fem2d_cc cen
);

fem2d_err fem2d_ref2mesh
(
    const fem2d_ea  ea,
    const matlib_xm vphi,
          fem2d_cc* x
);

fem2d_err fem2d_xinterp
(
    const fem2d_ea  ea,
    const matlib_xv u_nodes,
    const matlib_xm vphi,
          matlib_xv u_interp
);

fem2d_err fem2d_zinterp
(
    const fem2d_ea  ea,
    const matlib_zv u_nodes,
    const matlib_xm vphi,
          matlib_zv u_interp
);

matlib_real fem2d_xnormL2
(
    fem2d_ea  ea,
    matlib_xv u_qnodes,
    matlib_xv quadW
);

matlib_real fem2d_znormL2
(
    fem2d_ea  ea,
    matlib_zv u_qnodes,
    matlib_xv quadW
);

matlib_real fem2d_xiprod 
/* Inner product iprod = (u, conj(v))_{\Omega} */ 
(
    fem2d_ea  ea,
    matlib_xv u_qnodes,
    matlib_xv v_qnodes,
    matlib_xv quadW
);

matlib_real fem2d_ziprod 
/* Inner product iprod = (u, conj(v))_{\Omega} */ 
(
    fem2d_ea  ea,
    matlib_zv u_qnodes,
    matlib_zv v_qnodes,
    matlib_xv quadW
);

fem2d_err fem2d_quadP
/* Quadrature matrix for computing projection on vertex functions 
 * */ 
(
    matlib_xm  vphi,
    matlib_xv  quadW,
    matlib_xm* quadP
);
matlib_real fem2d_check_quadP
(
    fem2d_cc     xi,
    matlib_xv    quadW,
    matlib_index m,
    matlib_index n
);

fem2d_err fem2d_xprj
(
    fem2d_ea  ea,
    matlib_xv u_qnodes,
    matlib_xm quadP,
    matlib_xv u_prj
);

fem2d_err fem2d_zprj
(
    fem2d_ea  ea,
    matlib_zv u_qnodes,
    matlib_xm quadP,
    matlib_zv u_prj
);

/* function operating on nodal basis representation 
 * */ 
fem2d_err fem2d_NB_xprj
(
    fem2d_ea  ea,
    matlib_xv u_nodes,
    matlib_xv u_prj
);

fem2d_err fem2d_NB_zprj
(
    fem2d_ea  ea,
    matlib_zv u_nodes,
    matlib_zv u_prj
);

matlib_real fem2d_NB_xnormL2
(
    fem2d_ea  ea,
    matlib_xv u_nodes 
);

matlib_real fem2d_NB_znormL2
(
    fem2d_ea  ea,
    matlib_zv u_nodes
);

matlib_real fem2d_NB_xiprod
(
    fem2d_ea  ea,
    matlib_xv u_nodes,
    matlib_xv v_nodes
);
matlib_complex fem2d_NB_ziprod
(
    fem2d_ea  ea,
    matlib_zv u_nodes,
    matlib_zv v_nodes
);
/* ======== */ 

fem2d_err fem2d_quadM
(
    matlib_xm  vphi,
    matlib_xv  quadW,
    matlib_xm* Q
);

matlib_real fem2d_poly_symint
(
    matlib_index m,
    matlib_index n
);
matlib_real fem2d_check_quadM
(
    fem2d_cc     xi,
    matlib_xv    quadW,
    matlib_index m,
    matlib_index n
);

matlib_index fem2d_get_nnz(fem2d_ea ea);

fem2d_err fem2d_GMMSparsity
(
    fem2d_ea      ea,
    matlib_index* row, /* required to be pre-allocated */                     
    matlib_index* col  /* required to be pre-allocated */                  
);


fem2d_err fem2d_XCSRGMM
(
    fem2d_ea        ea,
    matlib_xv       q,
    matlib_index*   row,                     
    matlib_index*   col,                     
    matlib_real*    ugpmm                   
);

fem2d_err fem2d_ZCSRGMM
(
    fem2d_ea        ea,
    matlib_zv       q,
    matlib_index*   row,                     
    matlib_index*   col,                     
    matlib_complex* ugpmm                   
);

fem2d_err fem2d_xm_sparse_GMM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem2d_ea          ea,
    matlib_xm         Q,
    matlib_xv         phi,
    matlib_xm_sparse* M
);

fem2d_err fem2d_zm_sparse_GMM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem2d_ea          ea,
    matlib_xm         Q,
    matlib_zv         phi,
    matlib_zm_sparse* M
);

fem2d_err fem2d_SI_coeff
/* Stiffnes Integral Coefficients */ 
(
    fem2d_ea   ea,
    matlib_xv* S
);

fem2d_err fem2d_xm_sparse_GSM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem2d_ea          ea,
    matlib_xv         S,
    matlib_xv         quadW,
    matlib_xv         phi,
    matlib_xm_sparse* M
);

fem2d_err fem2d_zm_sparse_GSM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem2d_ea          ea,
    matlib_xv         S,
    matlib_xv         quadW,
    matlib_zv         phi,
    matlib_zm_sparse* M
);

fem2d_err fem2d_getmesh
(
    char* file_name,
    fem2d_cc* nodes,
    matlib_nv* iv
);

/*============================================================================+/
 | Quadrature on triangular domain
/+============================================================================*/

fem2d_err fem2d_symq
( 
    matlib_index n, 
    fem2d_cc*    xi, 
    matlib_xv*   quadW 
);


#endif

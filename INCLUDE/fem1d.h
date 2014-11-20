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
typedef matlib_err fem1d_err;

#define FEM1D_SUCCESS MATLIB_SUCCESS
#define FEM1D_FAILURE MATLIB_FAILURE
#define FEM1D_ERR_ENUM2STR MATLIB_ERR_ENUM2STR


/* 
 * xi in [-1, 1]
 * nodal basis 1D: l0(xi) = (1-xi)/2, l1(xi) = (1+xi)/2
 *
 *
 * */ 

typedef struct
{
    matlib_real*  elem_p;
    matlib_index len; /* nr of nodes */ 
    matlib_real* jacob;
    matlib_index nr_domains;      /* nr of domains = nr_nodes - 1 */ 

} fem1d_na; /* node array */ 


/* Ref. domain vertices */
extern const matlib_real FEM1D_VERT[2];


extern const matlib_index FEM1D_DIM;
/* Vertex indices */ 
extern const matlib_index FEM1D_INDEX_V1;
extern const matlib_index FEM1D_INDEX_V2;
extern const matlib_index FEM1D_NV;
extern const matlib_real  FEM1D_MEMI[2][2];
extern const matlib_real  FEM1D_MESI[2][2];

extern const matlib_index FEM1D_INDEX_V11;
extern const matlib_index FEM1D_INDEX_V12;
extern const matlib_index FEM1D_INDEX_V22;
extern const matlib_index FEM1D_NR_COMBI;

/*============================================================================*/

fem1d_err fem1d_refbasis
(
    const matlib_xv   xi,
          matlib_xm* vphi
);

fem1d_err fem1d_create_na
(
    fem1d_na *na, 
    matlib_index nr_nodes,
    matlib_real* nodes
);

fem1d_err fem1d_calc_jacobian(fem1d_na *na);

fem1d_err fem1d_ref2mesh
(
    fem1d_na   na,
    matlib_xm  vphi,
    matlib_xv* x
);

fem1d_err fem1d_xinterp
(
    const fem1d_na  na,
    const matlib_xv u_nodes,
    const matlib_xm vphi,
          matlib_xv u_interp
);

fem1d_err fem1d_zinterp
(
    const fem1d_na  na,
    const matlib_zv u_nodes,
    const matlib_xm vphi,
          matlib_zv u_interp
);

matlib_real fem1d_xnormL2
(
    fem1d_na  na,
    matlib_xv u_qnodes,
    matlib_xv quadW
);

matlib_real fem1d_znormL2
(
    fem1d_na  na,
    matlib_zv u_qnodes,
    matlib_xv quadW
);

fem1d_err fem1d_quadP
/* Quadrature matrix for computing projection on vertex functions 
 * */ 
(
    matlib_xm  vphi,
    matlib_xv  quadW,
    matlib_xm* quadP
);

fem1d_err fem1d_xprj
(
    fem1d_na  na,
    matlib_xv u_qnodes, /* values at the quadrature nodes */ 
    matlib_xm quadP,
    matlib_xv u_prj

);

fem1d_err fem1d_periodic_xprj
(
    fem1d_na  na,
    matlib_xv u_qnodes, /* values at the quadrature nodes */ 
    matlib_xm quadP,
    matlib_xv u_prj
);
fem1d_err fem1d_periodic_zprj
(
    fem1d_na  na,
    matlib_zv u_qnodes, /* values at the quadrature nodes */ 
    matlib_xm quadP,
    matlib_zv u_prj
);

fem1d_err fem1d_NB_xprj
/* Projection from nodal basis representation */ 
(
    fem1d_na  na,
    matlib_xv u_nodes, /* values at the nodes */ 
    matlib_xv u_prj
);

fem1d_err fem1d_NB_zprj
(
    fem1d_na  na,
    matlib_zv u_nodes, /* values at the nodes */ 
    matlib_zv u_prj
);

fem1d_err fem1d_periodic_NB_zprj
(
    fem1d_na  na,
    matlib_zv u_nodes, /* values at the nodes */ 
    matlib_zv u_prj
);
matlib_real fem1d_NB_xnormL2
(
    fem1d_na  na,
    matlib_xv u_nodes /* values at the nodes */ 
);
matlib_real fem1d_NB_znormL2
(
    fem1d_na  na,
    matlib_zv u_nodes /* values at the nodes */ 
);
matlib_real fem1d_NB_xiprod
(
    fem1d_na  na,
    matlib_xv u_nodes, /* values at the nodes */ 
    matlib_xv v_nodes
);
matlib_complex fem1d_NB_ziprod
(
    fem1d_na  na,
    matlib_zv u_nodes, /* values at the nodes */ 
    matlib_zv v_nodes
);


fem1d_err fem1d_quadM
(
    matlib_xm  vphi,
    matlib_xv  quadW,
    matlib_xm* Q
);


matlib_index fem1d_get_nnz(fem1d_na na);
matlib_index fem1d_periodic_get_nnz(fem1d_na na);

fem1d_err fem1d_GMMSparsity
(
    fem1d_na      na,
    matlib_index* row, /* required to be pre-allocated */                     
    matlib_index* col  /* required to be pre-allocated */                  
);

fem1d_err fem1d_periodic_GMMSparsity
(
    fem1d_na      na,
    matlib_index* row, /* required to be pre-allocated */                     
    matlib_index* col  /* required to be pre-allocated */                  
);

fem1d_err fem1d_XCSRGMM
(
    fem1d_na      na,
    matlib_xv     q,
    matlib_index* row,                     
    matlib_real*  ugpmm                   
);

fem1d_err fem1d_periodic_XCSRGMM
(
    fem1d_na      na,
    matlib_xv     q,
    matlib_index* row,                     
    matlib_real*  ugpmm                   
);

fem1d_err fem1d_ZCSRGMM
(
    fem1d_na      na,
    matlib_zv     q,
    matlib_index* row,                     
    matlib_complex* ugpmm                   
);

fem1d_err fem1d_periodic_ZCSRGMM
(
    fem1d_na      na,
    matlib_zv     q,
    matlib_index* row,                     
    matlib_complex* ugpmm                   
);

fem1d_err fem1d_xm_sparse_GMM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xm         Q,
    matlib_xv         phi,
    matlib_xm_sparse* M
);

fem1d_err fem1d_periodic_xm_sparse_GMM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xm         Q,
    matlib_xv         phi,
    matlib_xm_sparse* M
);

fem1d_err fem1d_zm_sparse_GMM
/* Complex - Assemble Global Mass Matrix*/ 
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xm         Q,
    matlib_zv         phi,
    matlib_zm_sparse* M
);

fem1d_err fem1d_periodic_zm_sparse_GMM
/* Complex - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xm         Q,
    matlib_zv         phi,
    matlib_zm_sparse* M
);

fem1d_err fem1d_xm_sparse_GSM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xv         quadW,
    matlib_xv         phi,
    matlib_xm_sparse* M
);

fem1d_err fem1d_periodic_xm_sparse_GSM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xv         quadW,
    matlib_xv         phi,
    matlib_xm_sparse* M
);

fem1d_err fem1d_zm_sparse_GSM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xv         quadW,
    matlib_zv         phi,
    matlib_zm_sparse* M
);

fem1d_err fem1d_periodic_zm_sparse_GSM
/* Real - Assemble Global Mass Matrix*/ 
(
    fem1d_na          na,
    matlib_xv         quadW,
    matlib_zv         phi,
    matlib_zm_sparse* M
);
#endif

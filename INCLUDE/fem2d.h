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


typedef struct
{
    matlib_index len; /* nr. of domains in the vertex patch */ 
    FEM2D_POINT_T point_enum;
    
    /* array of domain indices arrange such that neighbhouring 
     * domains are next to each other in the list 
     * */ 
    matlib_index* domain_index; 
    matlib_index* vert_index; /* array of vertex order indices */ 

    /* index of vertex on the boundary not shared by the next neighbhoiring
     * domain in the vertex patch
     * */ 
    matlib_index* bvert_index; 

} fem2d_vp; /* vertex patch */ 


typedef struct
{
    matlib_real** vert_p; /* pointer to the vertices */ 
    matlib_index  domain_index;
    matlib_real   jacob;
    matlib_real*  ijmat;
    matlib_real*  jmat;

} fem2d_te; /* triangular element */


typedef struct
{
    matlib_real* nbase; /* node base address */ 
    matlib_index len; /* nr of domains */ 
    matlib_index nr_nodes; /* nr of nodes */ 
    fem2d_te*    elem_p;
    fem2d_vp*    vpatch_p;

} fem2d_ea; /* element array type */ 

/* Cartesian coordinate indices */ 
extern const matlib_index FEM2D_INDEX_DIM1;
extern const matlib_index FEM2D_INDEX_DIM2;
extern const matlib_index FEM2D_DIM;

/* Jacobian indices */ 
extern const matlib_index FEM2D_INDEX_J11;
extern const matlib_index FEM2D_INDEX_J12;
extern const matlib_index FEM2D_INDEX_J21;
extern const matlib_index FEM2D_INDEX_J22;

/* Vertex indices */ 
extern const matlib_index FEM2D_INDEX_V1;
extern const matlib_index FEM2D_INDEX_V2;
extern const matlib_index FEM2D_INDEX_V3;
extern const matlib_index FEM2D_NV;
extern const matlib_index FEM2D_INIT_VPATCH_SIZE;
extern const matlib_real MEMI[3][3];
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

fem2d_err fem2d_interp
(
    const fem2d_ea  ea,
    const matlib_zv u_nodes,
    const matlib_xm vphi,
          matlib_zv u_interp
);

matlib_real fem2d_normL2
(
    fem2d_ea  ea,
    matlib_zv u_qnodes,
    matlib_xv quadW
);

matlib_real fem2d_iprod 
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

fem2d_err fem2d_prj
(
    fem2d_ea  ea,
    matlib_zv u_qnodes, /* values at the quadrature nodes */ 
    matlib_xm quadP,
    matlib_zv u_prj

);

fem2d_err fem2d_LEprj
(
    fem2d_ea  ea,
    matlib_zv u_nodes, /* values at the nodes */ 
    matlib_zv u_prj
);

matlib_real fem2d_LEnormL2
(
    fem2d_ea  ea,
    matlib_zv u_nodes /* values at the nodes */ 
);

matlib_complex fem2d_LEiprod
(
    fem2d_ea  ea,
    matlib_zv u_nodes, /* values at the nodes */ 
    matlib_zv v_nodes
);
#endif

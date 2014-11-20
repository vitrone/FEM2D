#ifndef LEGENDRE_H
#define LEGENDRE_H

/*============================================================================+/
 | Include all the dependencies
/+============================================================================*/
#include "basic.h"
#include "matlib.h"
#include "debug.h"
#include "ehandler.h"

/*============================================================================*/

matlib_real Legendre_LGL( matlib_real x, matlib_index n);

void CalcLP
( 
    matlib_index n, 
    matlib_real* LP, 
    matlib_real  x, 
    matlib_real  C[n-1],
    matlib_real  D[n-1]
);

void find_Gauss_points
( 
    matlib_index n, 
    matlib_real tol, 
    matlib_real* zeros
);

void find_LGL_points
( 
    matlib_index n, 
    matlib_real  tol, 
    matlib_real* zeros, 
    matlib_real* quadW,
    matlib_real* gzeros
);

void backward_transform_matrix
(
    const matlib_index P,            
    const matlib_index p,            
    const matlib_real* x,          
          matlib_real* pILTM
);
void forward_transform_matrix
(                                                                 
    const matlib_index   p,            
    const matlib_real* zeros,      
          matlib_real* pFLTM
);

void forward_transform_matrix2
(
    const matlib_index p,                        
    const matlib_index P,                        
    const matlib_real* zeros,                  
          matlib_real* pFLTM                   
);

/*======================================================================*/
void backward_transform_matrix_colmajor
( 
    matlib_index P, 
    matlib_index p, 
    matlib_real* x, 
    matlib_real* pILTM
);
void forward_transform_matrix_colmajor
( 
    matlib_index p, 
    matlib_real* zeros, 
    matlib_real* pFLTM
);
void forward_transform_matrix2_colmajor
( 
    matlib_index p, 
    matlib_index P, 
    matlib_real* zeros, 
    matlib_real* pFLTM 
);
/*============================================================================*/

void legendre_LGLdataLT1
( 
    const matlib_index p, 
    const matlib_real  tol,
          matlib_xv*   zeros,
          matlib_xv*   quadW
);
void legendre_LGLdataLT2
( 
    const matlib_index p, 
    const matlib_real  tol,
          matlib_xv*   zeros,
          matlib_xv*   quadW,
          matlib_xm*   FM,
          matlib_xm*   IM
);
void legendre_LGLdataFM
( 
    const matlib_xv xi,
          matlib_xm FM
);
void legendre_LGLdataIM
( 
    const matlib_xv xi,
          matlib_xm IM
);
#endif

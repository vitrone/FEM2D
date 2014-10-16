#ifndef PFEM1D_H
#define PFEM1D_H

/*============================================================================+/
 | Include all the dependencies
/+============================================================================*/
#include "basic.h"
#include "matlib.h"
#include "pthpool.h"
#include "debug.h"
#include "ehandler.h"
/*============================================================================+/
 |DATA STRUCTURES AND ENUMS
/+============================================================================*/
void pfem1d_XFLT
(
    const matlib_index    N,
    const matlib_xm       FM,
          matlib_xv       u,
          matlib_xv       U,
          matlib_index    num_threads,
          pthpool_data_t* mp
);

void pfem1d_ZFLT
(
    const matlib_index    N,
    const matlib_xm       FM,
          matlib_zv       u,
          matlib_zv       U,
          matlib_index    num_threads,
          pthpool_data_t* mp
);

void pfem1d_XILT
(
    const matlib_index    N,
    const matlib_xm       IM,
          matlib_xv       U,
          matlib_xv       u,
          matlib_index    num_threads,
          pthpool_data_t* mp
);

void pfem1d_ZILT
(
    const matlib_index    N,
    const matlib_xm       IM,
          matlib_zv       U,
          matlib_zv       u,
          matlib_index    num_threads,
          pthpool_data_t* mp
);

void pfem1d_XFLT2
(
    const matlib_index    N,
    const matlib_xm       FM,
          matlib_xm       u,
          matlib_xm       U,
          matlib_index    num_threads,
          pthpool_data_t* mp
);
void pfem1d_ZFLT2
(
    const matlib_index    N,
    const matlib_xm       FM,
          matlib_zm       u,
          matlib_zm       U,
          matlib_index    num_threads,
          pthpool_data_t* mp
);
void pfem1d_XF2L
(
    const matlib_index p, 
    const matlib_xv    vb,
          matlib_xv    u,
          matlib_index    num_threads,
          pthpool_data_t* mp
);

void pfem1d_ZF2L
(
    const matlib_index    p, 
    const matlib_zv       vb,
          matlib_zv       u,
          matlib_index    num_threads,
          pthpool_data_t* mp
);

void pfem1d_XPrjL2F
(
    matlib_index    p,
    matlib_xv       u,
    matlib_xv       Pvb,
    matlib_index    num_threads,
    pthpool_data_t* mp
);

void pfem1d_ZPrjL2F
(
    matlib_index    p,
    matlib_zv       u,
    matlib_zv       Pvb,
    matlib_index    num_threads,
    pthpool_data_t* mp
);

matlib_real pfem1d_XNorm2
(
    matlib_index    p,
    matlib_index    N,
    matlib_xv       u,
    matlib_index    num_threads,
    pthpool_data_t* mp
);

matlib_real pfem1d_ZNorm2
(
    matlib_index    p,
    matlib_index    N,
    matlib_zv       u,
    matlib_index    num_threads,
    pthpool_data_t* mp
);


void pfem1d_xm_nsparse_GMM
/* Double - Assemble Global Mass Matrix*/ 
(
    matlib_index       p,
    matlib_index       N,
    matlib_xm          Q,
    matlib_xm*         phi,
    matlib_xm*         q,
    matlib_xm_nsparse* M,
    matlib_index       num_threads,
    pthpool_data_t*    mp
    
);

void pfem1d_zm_nsparse_GMM
/* Double - Assemble Global Mass Matrix*/ 
(
    matlib_index       p,
    matlib_index       N,
    matlib_xm          Q,
    matlib_zm*         phi,
    matlib_zm*         q,
    matlib_zm_nsparse* M,
    matlib_index       num_threads,
    pthpool_data_t*    mp
    
);







#endif

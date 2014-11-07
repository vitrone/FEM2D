#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include "mkl.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include "mkl_pardiso.h"

//#define NDEBUG
#define MATLIB_NTRACE_DATA

#include "fem2d.h"

/* CUnit modules */
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <CUnit/TestDB.h>

/*============================================================================*/

int init_suite(void)
{
      return 0;
}
int clean_suite(void)
{
      return 0;
}              
/*============================================================================*/

static const matlib_real TOL = 1e-6;

fem2d_err Gaussian_zfunc
(
    fem2d_cc nodes_tmp, 
    matlib_zv u_nodes,
    matlib_real t
)
{
    debug_enter( "nodes: %d, u_nodes: %d", 
                 nodes_tmp.len, u_nodes.len);
    err_check(( nodes_tmp.len != u_nodes.len), clean_up, 
                "Dimension mismatch (nodes: %d, u_nodes: %d)!", 
                nodes_tmp.len, u_nodes.len);
    
    matlib_real* ptr;
    matlib_complex* uptr = u_nodes.elem_p;

    matlib_real a = 5.0;
    matlib_complex factor = 1.0/(1.0 + 4.0 * a * I * t);

    for ( ptr = nodes_tmp.elem_p; 
          ptr < (nodes_tmp.elem_p + 2 * (nodes_tmp.len)); ptr+=2, uptr++)
    {

        *uptr = factor * cexp( - a * factor * (ptr[0] * ptr[0] + ptr[1] * ptr[1]));
    }
    debug_exit("exit status: %d", FEM2D_SUCCESS);
    return FEM2D_SUCCESS;
clean_up:
    debug_exit("exit status: %d", FEM2D_FAILURE);
    return FEM2D_FAILURE;
}


void pde2d_solve_FSE(void)
{

    fem2d_err error = FEM2D_SUCCESS;
    matlib_index nr_domains;
    matlib_nv iv;
    fem2d_cc my_nodes;
    char* file_name = "mesh_data.bin";
    fem2d_getmesh(file_name, &my_nodes, &iv);
    nr_domains = iv.len/FEM2D_NV;

    matlib_index mcnt = 0;

    /* ea: element array */ 
    fem2d_ea ea;
    error = fem2d_create_ea(my_nodes, iv.elem_p, nr_domains, &ea);
    error = fem2d_create_vp(&ea);

    /* vphi: matrix containing values of reference basis functions at 
     * quadrature points
     * */ 
    /* quadW: quadrature weights */ 
    matlib_index D = 4;
    matlib_xv quadW;
    fem2d_cc xi_D;
    fem2d_symq(D, &xi_D, &quadW);

    matlib_xm vphi;
    error = fem2d_refbasis(xi_D, &vphi);

    /* x_qnodes: mesh points for interpolation 
     * */ 
    fem2d_cc x_qnodes;
    error = fem2d_ref2mesh(ea, vphi, &x_qnodes);

    matlib_xm quadP;
    error = fem2d_quadP( vphi, quadW, &quadP);
    
    matlib_xv SI;
    error = fem2d_SI_coeff( ea, &SI);
    matlib_xm Q;
    error = fem2d_quadM( vphi, quadW, &Q);

    matlib_real dt   = 0.5e-3;
    matlib_real rho  = 2.0/dt;
    matlib_real irho = dt/2.0;

    /* Assemble the mass-matrix */ 
    matlib_zv phi;
    error = matlib_create_zv( x_qnodes.len, &phi , MATLIB_COL_VECT);

    matlib_index i;
    for (i = 0; i < x_qnodes.len; i++)
    {
        phi.elem_p[i] = 1.0;
    }

    matlib_zm_sparse M;
    error = fem2d_zm_sparse_GMM( ea, Q, phi, &M);

    for (i = 0; i < x_qnodes.len; i++)
    {
        phi.elem_p[i] = irho*I;
    }
    matlib_zm_sparse S;
    error = fem2d_zm_sparse_GSM( ea, SI, quadW, phi, &S);

    matlib_zv VM = {.len = M.rowIn[M.lenc], .elem_p = M.elem_p};
    matlib_zv VS = {.len = S.rowIn[S.lenc], .elem_p = S.elem_p};

    error = matlib_zaxpy(1.0, VS, VM);

    matlib_zv u_nodes, v_nodes, u_prj;
    error = matlib_create_zv( my_nodes.len, &u_nodes , MATLIB_COL_VECT);
    error = matlib_create_zv( my_nodes.len, &v_nodes , MATLIB_COL_VECT);
    error = matlib_create_zv( my_nodes.len, &u_prj   , MATLIB_COL_VECT);

    matlib_real t = 0;
    error = Gaussian_zfunc(my_nodes, u_nodes, t);
    error = fem2d_NB_zprj(ea, u_nodes, u_prj);
    
    matlib_index Nt = 100;
    /* Initialize solver data */ 
    pardiso_solver_t data = { .nsparse  = 1, 
                              .mnum     = 1, 
                              .mtype    = PARDISO_COMPLEX_SYM,
                              .sol_enum = PARDISO_LHS, 
                              .smat_p   = (void*)&M,
                              .rhs_p    = (void*)&u_prj,
                              .sol_p    = (void*)&v_nodes};

    debug_body("%s", "Solver data initialized");
    data.phase_enum = PARDISO_INIT;
    matlib_pardiso(&data);

    data.phase_enum = PARDISO_ANALYSIS_AND_FACTOR;
    matlib_pardiso(&data);
    
    matlib_real e_relative, norm_actual;

    for (i = 1; i < Nt; i++)
    {
        data.phase_enum = PARDISO_SOLVE_AND_REFINE;
        matlib_pardiso(&data);

        error = matlib_zaxpby(2.0, v_nodes, -1.0, u_nodes);
        error = fem2d_NB_zprj(ea, u_nodes, u_prj);

        t += dt;
        error = Gaussian_zfunc(my_nodes, v_nodes, t);
        norm_actual = fem2d_NB_znormL2( ea, v_nodes);
        debug_body("Actual norm: %0.16f", norm_actual);
        
        error = matlib_zaxpy(-1.0, u_nodes, v_nodes);

        e_relative = fem2d_NB_znormL2( ea, v_nodes)/norm_actual;
        debug_body("Relative error: %0.16g", e_relative);

    }

    data.phase_enum = PARDISO_FREE;
    matlib_pardiso(&data);

    matlib_free(M.elem_p);
    matlib_free(M.rowIn);
    matlib_free(M.colIn);

    matlib_free(S.elem_p);
    matlib_free(S.rowIn);
    matlib_free(S.colIn);

    matlib_free(u_prj.elem_p);
    matlib_free(u_nodes.elem_p);
    matlib_free(v_nodes.elem_p);
    matlib_free(x_qnodes.elem_p);
    matlib_free(vphi.elem_p);
    fem2d_free_ea(ea);
}










/*============================================================================+/
 |
 |
 |
/+============================================================================*/

int main()
{
    CU_pSuite pSuite = NULL;

    /* initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
    {
        return CU_get_error();
    }

    /* Create a test array */
    CU_TestInfo test_array[] = 
    {
#if 0
#endif
        { "Solve FSE", pde2d_solve_FSE},
        CU_TEST_INFO_NULL,
    };

    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "pde2d", init_suite, clean_suite, NULL, NULL, (CU_TestInfo*) test_array },
        CU_SUITE_INFO_NULL
    }; 

    /* Register test suites */ 
    CU_ErrorCode CU_error = CU_register_suites(suites); 
    if (CU_error != CUE_SUCCESS) 
    {
        debug_body("%s", CU_get_error_msg());
        CU_cleanup_registry();
        return CU_get_error();
    }


   /* Run all tests using the CUnit Basic interface */
   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   CU_cleanup_registry();
   return CU_get_error();
}


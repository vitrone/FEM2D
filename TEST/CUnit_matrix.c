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

#include "matlib.h"

/* CUnit modules */
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <CUnit/TestDB.h>

/*============================================================================*/

static const matlib_real TOL = 1e-9;

int init_suite(void)
{
      return 0;
}
int clean_suite(void)
{
      return 0;
}              
/*============================================================================*/
/* In C99 complex numbers are defined as an adjacent pair of memory locations.
 * However, there are some syntactical differences which require type-casting:
 * (1) Interpreting a complex array as N-by-2 matlib_real matrix in rwo major format.
 * (2) Interpreting a matlib_real matrix in row major format as an array of complex
 *     numbers.
 *
 * */ 

/* The following type is fully compatible with C99 complex.*/ 
typedef matlib_real matlib_c[2];

void test_complex(void)
{
    matlib_c xa[] = { {1.0, 2.0 }, 
                      {3.0, 4.5 }};

    debug_body("xa[0] = %0.16f, %0.16f", (*(xa+0))[0], (*(xa+0))[1]);
    debug_body("xa[1] = %0.16f, %0.16f", (*(xa+1))[0], (*(xa+1))[1]);
    matlib_xm xm = { .lenc   = 2,
                     .lenr  = 2,
                     .order  = MATLIB_ROW_MAJOR,
                     .elem_p = xa[0]};
    
    DEBUG_PRINT_XM(xm, "%s", "");
    matlib_zv xcm = { .len    = 2,
                      .type   = MATLIB_ROW_MAJOR,
                      .elem_p = (complex*)xa};
    
    DEBUG_PRINT_ZV(xcm, "%s", "");
    /* convert such xa to complex type */
    matlib_c ya[2];
    *((complex*)ya[0]) = csqrt(*((complex*)xa[0]));
    *((complex*)ya[1]) = csqrt(*((complex*)xa[1]));
    debug_body("ya[0]: %0.16f %+0.16fi", *((complex*)ya[0]));
    debug_body("ya[1]: %0.16f %+0.16fi", *((complex*)ya[1]));

    matlib_xm ym = { .lenc   = 2,
                     .lenr  = 2,
                     .order  = MATLIB_ROW_MAJOR,
                     .elem_p = &(*ya)[0]};
    DEBUG_PRINT_XM(ym, "%s", "");
    matlib_zv ycv = { .len    = 2,
                      .type   = MATLIB_ROW_MAJOR,
                      .elem_p = (complex*)ya};
    
    DEBUG_PRINT_ZV(ycv, "%s", "");
    
    /* Verifying if C99 complex can be converted to matlib_real matrix.*/ 

    matlib_complex za[4] = { -1.0 + I * 2.0, 
                              1.6 + I * 2.1};

    matlib_zv z = {  .len   = 2, 
                     .type  = MATLIB_COL_VECT,
                     .elem_p = za};
    
    DEBUG_PRINT_ZV(z, "%s", "");

    debug_body("za[0] = %0.16f %+0.16fi", *((matlib_real*)(za+0)), *((matlib_real*)(za+0)+1));
    debug_body("za[1] = %0.16f %+0.16fi", *((matlib_real*)(za+1)), *((matlib_real*)(za+1)+1));

    matlib_xm zm = { .lenc  = 2,
                     .lenr  = 2,
                     .order   = MATLIB_ROW_MAJOR,
                     .elem_p = ((matlib_real*)za)};

    DEBUG_PRINT_XM(zm, "%s", "");

    CU_ASSERT_TRUE(true);

}

/*============================================================================*/
void test_zgemv(void)
{
    debug_enter("%s", "");
    matlib_complex Ma[3][4] = { { 3.0 + I*1.0, 1.0+ I*1.0, 3.0+ I*1.0, 2.0+ I*1.0},
                                { 1.0 + I*1.0, 5.0+ I*1.0, 9.0+ I*1.0, 4.0+ I*1.0},
                                { 2.0 + I*1.0, 6.0+ I*1.0, 5.0+ I*1.0, 7.0+ I*1.0}
                          };

    matlib_zm M = { .lenc   = 3, 
                    .lenr   = 4, 
                    .order  = MATLIB_ROW_MAJOR,
                    .op     = MATLIB_NO_TRANS,
                    .elem_p = Ma[0]};
    DEBUG_PRINT_ZM(M, "%s", "");

    matlib_complex xa[] = {-1.0 + I * 2.0, 
                           -1.0 + I * 2.0, 
                            1.0 + I * 2.0, 
                            1.0 + I * 2.0 };

    matlib_zv x = { .len    = 4, 
                    .type   = MATLIB_COL_VECT, 
                    .elem_p = xa};
    DEBUG_PRINT_ZV(x, "%s", "");

    matlib_complex ya[] = { 0, 0, 0 };

    matlib_zv y = { .len    = 3, 
                    .type   = MATLIB_COL_VECT, 
                    .elem_p = ya};

    matlib_zgemv( 1.0, M, x, 0, y);
    matlib_complex ya_actual[] = { Ma[0][0] * xa[0] + Ma[0][1] * xa[1] + Ma[0][2] * xa[2]+ Ma[0][3] * xa[3], 
                                   Ma[1][0] * xa[0] + Ma[1][1] * xa[1] + Ma[1][2] * xa[2]+ Ma[1][3] * xa[3],
                                   Ma[2][0] * xa[0] + Ma[2][1] * xa[1] + Ma[2][2] * xa[2]+ Ma[2][3] * xa[3]};

    matlib_zv y_actual = { .len    = 3, 
                           .type   = MATLIB_COL_VECT, 
                           .elem_p = ya_actual};

    DEBUG_PRINT_ZV(y, "%s", "");
    DEBUG_PRINT_ZV(y_actual, "%s", "");
    matlib_real norm_actual = matlib_znrm2(y_actual);
    matlib_zaxpy(-1.0, y_actual, y);
    matlib_real e_relative = matlib_znrm2(y)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);


}
/*============================================================================*/

void test_xgemv(void)
{
    debug_enter("%s", "");

    
    matlib_real Ma[3][4] = { { 3.0, 1.0, 3.0, 2.0},
                        { 1.0, 5.0, 9.0, 4.0},
                        { 2.0, 6.0, 5.0, 7.0}
                  };
    matlib_xm M = { .lenc   = 3, 
                    .lenr   = 4, 
                    .order  = MATLIB_ROW_MAJOR,
                    .op     = MATLIB_NO_TRANS,
                    .elem_p = &Ma[0][0]};

    DEBUG_PRINT_XM(M, "%s", "");
    
    matlib_real xa[] = {-1.0, -1.0, 1.0, 1.0};
    matlib_xv x = { .len    = 4, 
                    .type   = MATLIB_COL_VECT, 
                    .elem_p = xa};
    
    DEBUG_PRINT_XV(x, "%s", "");
    
    matlib_real ya[] = { 0, 0, 0 };

    matlib_xv y = { .len    = 3, 
                    .type   = MATLIB_COL_VECT, 
                    .elem_p = ya};


    matlib_xgemv( 1, M, x, 0, y);
    matlib_real ya_actual[] = { Ma[0][0] * xa[0] + Ma[0][1] * xa[1] + Ma[0][2] * xa[2]+ Ma[0][3] * xa[3], 
                           Ma[1][0] * xa[0] + Ma[1][1] * xa[1] + Ma[1][2] * xa[2]+ Ma[1][3] * xa[3],
                           Ma[2][0] * xa[0] + Ma[2][1] * xa[1] + Ma[2][2] * xa[2]+ Ma[2][3] * xa[3]};

    matlib_xv y_actual = { .len    = 3, 
                           .type   = MATLIB_COL_VECT, 
                           .elem_p = ya_actual};

    DEBUG_PRINT_XV(y, "%s", "");
    DEBUG_PRINT_XV(y_actual, "%s", "");
    matlib_real norm_actual = matlib_xnrm2(y_actual);
    matlib_xaxpy(-1.0, y_actual, y);
    matlib_real e_relative = matlib_xnrm2(y)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);
    
}

void test_xgemm(void)
{
    debug_enter("%s", "");

    
    matlib_real Ma[3][4] = { { 3.0, 1.0, 3.0, 2.0},
                        { 1.0, 5.0, 9.0, 4.0},
                        { 2.0, 6.0, 5.0, 7.0}
                  };
    matlib_xm M = { .lenc   = 3, 
                    .lenr   = 4, 
                    .order  = MATLIB_ROW_MAJOR,
                    .elem_p = &Ma[0][0]};

    DEBUG_PRINT_XM(M, "%s", "");
    
    matlib_real Na[4][3] = { {  2.0,  1.0,  3.0 },
                        {  1.0,  3.0, -7.0 },
                        {  2.0, -1.0,  5.0 },
                        { -2.0,  6.0,  8.0 }
                  };

    matlib_xm N = { .lenc   = 4, 
                    .lenr   = 3, 
                    .order  = MATLIB_ROW_MAJOR,
                    .elem_p = &Na[0][0]};
    
    DEBUG_PRINT_XM(N, "%s", "");
    matlib_real Pa[3][3];
    matlib_xm P = { .lenc   = 3, 
                    .lenr   = 3, 
                    .order  = MATLIB_ROW_MAJOR,
                    .elem_p = &Pa[0][0]};
    

    matlib_xgemm( 1, M, N, 0, P);
    DEBUG_PRINT_XM(P, "%s", "");

    matlib_real Pa_actual[3][3] = { { 9.0,    15.0,    33.0},
                               {17.0,    31.0,    45.0},
                               { 6.0,    57.0,    45.0}
                             };

    matlib_xm P_actual = { .lenc   = 3, 
                           .lenr   = 3, 
                           .order  = MATLIB_ROW_MAJOR,
                           .elem_p = &Pa_actual[0][0]};


    DEBUG_PRINT_XM(P_actual, "%s", "");
    matlib_xv u = MK_VM(P);
    matlib_xv u_actual = MK_VM(P_actual);

    matlib_real norm_actual = matlib_xnrm2(u_actual);
    matlib_xaxpy(-1.0, u_actual, u);
    matlib_real e_relative = matlib_xnrm2(u)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);
    
}
void test_zgemm(void)
{
    debug_enter("%s", "");

    
    matlib_complex Ma[3][4] = { { 3.0 + I*2.0, 1.0 + I*2.0, 3.0 + I*2.0, 2.0 + I*2.0},
                                { 1.0 + I*2.0, 5.0 + I*2.0, 9.0 + I*2.0, 4.0 + I*2.0},
                                { 2.0 + I*2.0, 6.0 + I*2.0, 5.0 + I*2.0, 7.0 + I*2.0}
                  };
    matlib_zm M = { .lenc   = 3, 
                    .lenr   = 4, 
                    .order  = MATLIB_ROW_MAJOR,
                    .elem_p = Ma[0]};

    DEBUG_PRINT_ZM(M, "%s", "");
    
    matlib_complex Na[4][3] = { {  2.0 + I*3.0,  1.0 + I*3.0,  3.0 + I*3.0 },
                                {  1.0 + I*3.0,  3.0 + I*3.0, -7.0 + I*3.0 },
                                {  2.0 + I*3.0, -1.0 + I*3.0,  5.0 + I*3.0 },
                                { -2.0 + I*3.0,  6.0 + I*3.0,  8.0 + I*3.0 }
                  };

    matlib_zm N = { .lenc   = 4, 
                    .lenr   = 3, 
                    .order  = MATLIB_ROW_MAJOR,
                    .elem_p = Na[0]};
    
    DEBUG_PRINT_ZM(N, "%s", "");
    matlib_complex Pa[3][3];
    matlib_zm P = { .lenc   = 3, 
                    .lenr   = 3, 
                    .order  = MATLIB_ROW_MAJOR,
                    .elem_p = Pa[0]};
    

    matlib_zgemm( 1, M, N, 0, P);
    DEBUG_PRINT_ZM(P, "%s", "");

    matlib_complex Pa_actual[3][3] = { {-15.0000 + I*33.0000,  -9.0000 + I*45.0000,   9.0000 + I*45.0000},
                                       {-7.0000  + I*63.0000,   7.0000 + I*75.0000,  21.0000 + I*75.0000},
                                       {-18.0000 + I*66.0000,  33.0000 + I*78.0000,  21.0000 + I*78.0000}};

    matlib_zm P_actual = { .lenc   = 3, 
                           .lenr   = 3, 
                           .order  = MATLIB_ROW_MAJOR,
                           .elem_p = Pa_actual[0]};


    DEBUG_PRINT_ZM(P_actual, "%s", "");
    
    matlib_zv u = MK_VM(P);
    matlib_zv u_actual = MK_VM(P);

    matlib_real norm_actual = matlib_znrm2(u_actual);
    matlib_zaxpy(-1.0, u_actual, u);
    matlib_real e_relative = matlib_znrm2(u)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);
    
}
/*============================================================================*/

 /* A1 -> 5-by-3 */ 

static matlib_real A1_rm[15] =
{
1.295198610346667e-01, 
4.371724372446585e-01,
1.568078875443956e-01,
8.798835824178787e-01,
3.798389169076282e-01,
3.260342840664192e-01,
4.407912846769657e-02,
9.796574398465744e-01,
3.140619286734648e-01,
6.867200184808028e-01,
3.989934226116294e-01,
8.945006639042554e-01,
7.337729049335920e-01,
4.401869142693676e-01,
2.470241792640397e-01
};

static matlib_real A1_cm[15] =
{
1.295198610346667e-01, 
8.798835824178787e-01,
4.407912846769657e-02,
6.867200184808028e-01,
7.337729049335920e-01,
4.371724372446585e-01,
3.798389169076282e-01,
9.796574398465744e-01,
3.989934226116294e-01,
4.401869142693676e-01,
1.568078875443956e-01,
3.260342840664192e-01,
3.140619286734648e-01,
8.945006639042554e-01,
2.470241792640397e-01
};

    /* A2 -> 3-by-2 */ 
static matlib_real A2_rm[6] =
{
3.106790159126492e-01,
1.436375950594925e-01,
4.088689080421710e-01,
8.713220665730218e-01,
7.080108940640530e-01,
8.315589675345991e-02};
    
static matlib_real A2_cm[6] =
{
3.106790159126492e-01,
4.088689080421710e-01,
7.080108940640530e-01,
1.436375950594925e-01,
8.713220665730218e-01,
8.315589675345991e-02
};
    
    /* A3 -> 5-by-2 */ 
static matlib_real A3_rm[10] =
{
3.300070126663544e-01,
4.125614133270096e-01,
6.595015136486246e-01,
4.844580650195631e-01,
6.366051948516696e-01,
8.860446663411344e-01,
1.009801719360099e+00,
5.206735903263244e-01,
5.828423969960751e-01,
5.094834643505544e-01};

static matlib_real A3_cm[10] =
{
3.300070126663544e-01,
6.595015136486246e-01,
6.366051948516696e-01,
1.009801719360099e+00,
5.828423969960751e-01,
4.125614133270096e-01,
4.844580650195631e-01,
8.860446663411344e-01,
5.206735903263244e-01,
5.094834643505544e-01
};

void test_xgemm_1a(void)
{
    /* 
     * All matrices are considered in row major format.
     * M3 = M1 * M2;
     *
     * */ 

    matlib_xm M1 = { .lenc   = 5, 
                     .lenr   = 3,
                     .order  = MATLIB_ROW_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = A1_rm};

    DEBUG_PRINT_XM( M1, "%s:", "M1");
    matlib_xm M2 = { .lenc   = 3, 
                     .lenr   = 2,
                     .order  = MATLIB_ROW_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = A2_rm};

    DEBUG_PRINT_XM( M2, "%s:", "M2");

    matlib_real A3_tmp[5][2];
    matlib_xm M3 = { .lenc   = 5, 
                     .lenr   = 2,
                     .order  = MATLIB_ROW_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = &A3_tmp[0][0]};

    matlib_xv VM3 = { .len    = M3.lenc*M3.lenr, 
                      .type   = MATLIB_ROW_VECT, 
                      .elem_p = M3.elem_p};

    matlib_xv VA3 = { .len    = 10, 
                      .type   = MATLIB_COL_VECT, 
                      .elem_p = A3_rm};

    matlib_xgemm( 1.0, M1, M2, 0.0, M3);
    DEBUG_PRINT_XM( M3, "%s:", "M3");

    matlib_real norm_actual = matlib_xnrm2(VA3);
    matlib_xaxpy(-1.0, VA3, VM3);
    matlib_real e_relative = matlib_xnrm2(VM3)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);


}

void test_xgemm_1b(void)
{
    /* 
     *
     * All matrices are considered in row major format.
     * M3 = M1^T * M2;
     *
     * */ 

    matlib_xm M1 = { .lenc   = 3, 
                     .lenr   = 5,
                     .order  = MATLIB_ROW_MAJOR,
                     .op     = MATLIB_TRANS, 
                     .elem_p = A1_cm}; /* Try to use the A_cm this time */ 

    DEBUG_PRINT_XM( M1, "%s:", "M1");
    matlib_xm M2 = { .lenc   = 3, 
                     .lenr   = 2,
                     .order  = MATLIB_ROW_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = A2_rm};

    DEBUG_PRINT_XM( M2, "%s:", "M2");

    matlib_real A3_tmp[5][2];
    matlib_xm M3 = { .lenc   = 5, 
                     .lenr   = 2,
                     .order  = MATLIB_ROW_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = &A3_tmp[0][0]};

    matlib_xv VM3 = { .len    = M3.lenc*M3.lenr, 
                      .type   = MATLIB_ROW_VECT, 
                      .elem_p = M3.elem_p};

    matlib_xv VA3 = { .len    = 10, 
                      .type   = MATLIB_COL_VECT, 
                      .elem_p = A3_rm};

    matlib_xgemm( 1.0, M1, M2, 0.0, M3);
    DEBUG_PRINT_XM( M3, "%s:", "M3");

    matlib_real norm_actual = matlib_xnrm2(VA3);
    matlib_xaxpy(-1.0, VA3, VM3);
    matlib_real e_relative = matlib_xnrm2(VM3)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);


}


void test_xgemm_2a(void)
{
    /* 
     * All matrices are considered in column major format.
     * M3 = M1 * M2;
     *
     * */ 

    matlib_xm M1 = { .lenc   = 5, 
                     .lenr   = 3,
                     .order  = MATLIB_COL_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = A1_cm};

    DEBUG_PRINT_XM( M1, "%s:", "M1");

    matlib_xm M2 = { .lenc   = 3, 
                     .lenr   = 2,
                     .order  = MATLIB_COL_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = A2_cm};

    DEBUG_PRINT_XM( M2, "%s:", "M2");

    matlib_real A3_tmp[10]   = {0.0,};
    matlib_xm M3 = { .lenc   = 5, 
                     .lenr   = 2,
                     .order  = MATLIB_COL_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = A3_tmp};

    matlib_xv VM3 = { .len    = M3.lenc*M3.lenr, 
                      .type   = MATLIB_COL_VECT, 
                      .elem_p = M3.elem_p};

    matlib_xv VA3 = { .len    = 10, 
                      .type   = MATLIB_COL_VECT, 
                      .elem_p = A3_cm};

    matlib_xgemm( 1.0, M1, M2, 0.0, M3);
    DEBUG_PRINT_XM( M3, "%s:", "M3");

    matlib_real norm_actual = matlib_xnrm2(VA3);
    matlib_xaxpy(-1.0, VA3, VM3);
    matlib_real e_relative = matlib_xnrm2(VM3)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);
}

void test_xgemm_2b(void)
{
    /* 
     * All matrices are considered in column major format.
     * M3 = M1^T * M2;
     *
     * */ 

    matlib_xm M1 = { .lenc   = 3, 
                     .lenr   = 5,
                     .order  = MATLIB_COL_MAJOR,
                     .op     = MATLIB_TRANS, 
                     .elem_p = A1_rm}; /* Try A1_rm this time */ 

    DEBUG_PRINT_XM( M1, "%s:", "M1");

    matlib_xm M2 = { .lenc   = 3, 
                     .lenr   = 2,
                     .order  = MATLIB_COL_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = A2_cm};

    DEBUG_PRINT_XM( M2, "%s:", "M2");

    matlib_real A3_tmp[10]   = {0.0,};
    matlib_xm M3 = { .lenc   = 5, 
                     .lenr   = 2,
                     .order  = MATLIB_COL_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = A3_tmp};

    matlib_xv VM3 = { .len    = M3.lenc*M3.lenr, 
                      .type   = MATLIB_COL_VECT, 
                      .elem_p = M3.elem_p};

    matlib_xv VA3 = { .len    = 10, 
                      .type   = MATLIB_COL_VECT, 
                      .elem_p = A3_cm};

    matlib_xgemm( 1.0, M1, M2, 0.0, M3);
    DEBUG_PRINT_XM( M3, "%s:", "M3");

    matlib_real norm_actual = matlib_xnrm2(VA3);
    matlib_xaxpy(-1.0, VA3, VM3);
    matlib_real e_relative = matlib_xnrm2(VM3)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);
}

void test_NaNs(void)
{
    float x = 0.0/0.0;
    debug_print("%f", x);
    debug_print("%d", isnan(x) );
#ifdef FP_NAN
    debug_print("%f", isnan(FP_NAN) );

#endif
    matlib_err error;
    printf("\n%s\n", MATLIB_ERR_ENUM2STR(error));

    CU_ASSERT_TRUE(true);
}


/*============================================================================+/
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
        { "Complex matrices"                     , test_complex },
        { "Real matrix-vector multiplication"    , test_xgemv },
        { "Complex matrix-vector multiplication" , test_zgemv },
        { "Real matrix-matrix multiplication"    , test_xgemm },
        { "Complex matrix-matrix multiplication" , test_zgemm },
        { "gemm 1a", test_xgemm_1a },
        { "gemm 1b", test_xgemm_1b },
        { "gemm 2a", test_xgemm_2a },
        { "gemm 2b", test_xgemm_2b },
        { "NaNs", test_NaNs },
        CU_TEST_INFO_NULL,
    };

    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "Matrix Library", init_suite, clean_suite, NULL, NULL, test_array },
        CU_SUITE_INFO_NULL,
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

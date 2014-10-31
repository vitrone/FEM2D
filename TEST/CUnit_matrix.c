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

#if 0
BEGIN_FILE_HEADER
    TIME_STAMP: YYYY-MM-DD-HH-MM-SS
    NR_MATRICES: 
END_FILE_HEADER

BEGIN_MATRIX_ENTRY: <matrix index> <number of lines>
    BEGIN_COMMENT: <number of lines>
    END_COMMENT
    FORMAT     : <format code>
    BEGIN_PPTY : <number of properties listed>
    <ppty index> <value>
    END_PPTY
    NR_ARRAY   : <number>
    BEGIN_ARRAY: <length> <data type>
    END_ARRAY
END_MATRIX_ENTRY
#endif

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

void test_zgemv1(void)
{
    debug_enter("%s", "");
    matlib_complex Ma[3][4] = { { 3.0 + I*1.0, 1.0+ I*1.0, 3.0+ I*1.0, 2.0+ I*1.0},
                                { 1.0 + I*1.0, 5.0+ I*1.0, 9.0+ I*1.0, 4.0+ I*1.0},
                                { 2.0 + I*1.0, 6.0+ I*1.0, 5.0+ I*1.0, 7.0+ I*1.0}};

    matlib_zm M = { .lenc   = 4, 
                    .lenr   = 3, 
                    .order  = MATLIB_COL_MAJOR,
                    .op     = MATLIB_TRANS,
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
                             { 2.0, 6.0, 5.0, 7.0}};
    matlib_xm M = { .lenc   = 3, 
                    .lenr   = 4, 
                    .order  = MATLIB_ROW_MAJOR,
                    .op     = MATLIB_NO_TRANS,
                    .elem_p = Ma[0]};

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
void test_xgemv1(void)
{
    debug_enter("%s", "");

    
    matlib_real Ma[3][4] = { { 3.0, 1.0, 3.0, 2.0},
                             { 1.0, 5.0, 9.0, 4.0},
                             { 2.0, 6.0, 5.0, 7.0}};

    matlib_xm M = { .lenc   = 4, 
                    .lenr   = 3, 
                    .order  = MATLIB_COL_MAJOR,
                    .op     = MATLIB_TRANS,
                    .elem_p = Ma[0]};

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
/*============================================================================*/
static void full_sym_mat
(
    matlib_xm_sparse sM,
    matlib_xm*      M
)
{
    matlib_index i, j;
    
    matlib_create_xm( sM.lenc, sM.lenr, M, MATLIB_ROW_MAJOR, MATLIB_NO_TRANS);    

    for(i=0; i<sM.lenc; i++)
    {
        j = sM.rowIn[i];
        M->elem_p[i*sM.lenr+sM.colIn[j]] = sM.elem_p[j];
        if(j != i)
        {
            M->elem_p[sM.colIn[j]*sM.lenr+i] = sM.elem_p[j];
        }

        for(j=sM.rowIn[i]+1; j<sM.rowIn[i+1]; j++)
        {
            M->elem_p[i*sM.lenr+sM.colIn[j]] = sM.elem_p[j];
            M->elem_p[sM.colIn[j]*sM.lenr+i] = sM.elem_p[j];
        }
    }

}
static matlib_index M_rows[11] = {0, 7, 11, 14, 15, 18, 20, 21, 23, 24, 25};
static matlib_index M_cols[25] = { 0, 1, 3, 5, 6, 7, 9, 1, 5, 7, 9, 5, 6,
                                7, 5, 5, 6, 7, 5, 6, 6, 8, 9, 8, 9};

static matlib_real M_vals[25] = {
     7.768417982069939e-01, 
    -4.203446056274865e-01,
     9.014914350929220e-01,
     6.216319090167377e-01,
    -1.101779357090595e-01,
    -7.586271567433341e-01,
     1.478350435489602e-01,
     1.077140113434178e+00,
    -5.013771251839112e-01,
     9.275842182056043e-01,
    -1.539453896657975e-01,
     1.788703729862607e-01,
     2.167055505203805e-01,
     5.475196872294628e-01,
    -2.077895861079801e-01,
     4.123077020405314e-01,
     4.369186280281671e-01,
     3.946755315682919e-01,
    -5.643766911304262e-01,
     2.702470598769673e+00,
     0                    ,
     4.854423718912431e-03,
    -1.801630971232774e-01,
    -2.259840351530994e+00,
    -2.084964142727976e-01};

static matlib_real spA[10] = {
     9.171936638298100e-01, 
     2.858390188203735e-01,
     7.572002291107213e-01,
     7.537290942784953e-01,
     3.804458469753567e-01,
     5.678216407252211e-01,
     7.585428956306361e-02,
     5.395011866660715e-02,
     5.307975530089727e-01,
     7.791672301020112e-01};

static matlib_real spB[10] = {
     1.690722586077362e+00, 
    -4.322474418019124e-01,
     1.475432663436200e-01,
     7.088548085546227e-01,
     2.885521797451608e-01,
     4.470567750389413e-01,
     1.763780154881953e+00,
    -3.734024189205709e-03,
    -1.199255832047895e+00,
    -8.058362796316272e-02};
    

void test_xcsrsymv(void)
{
    debug_enter("%s", "");

    matlib_xm_sparse M = { .lenc = 10, 
                           .lenr = 10,
                           .rowIn = M_rows, 
                           .colIn = M_cols, 
                           .elem_p = M_vals};

    debug_body("%s", "Sparse matrix declared");
    matlib_xv A = {.len = 10, .elem_p = spA };

    matlib_real Ba[10] = {0.0,};
    matlib_xv B = {.len = 10, .elem_p = Ba };
    matlib_xv B1 = {.len = 10, .elem_p = spB };

    debug_body("%s", "multiplying with sparse matrix");
    matlib_xcsrsymv(MATLIB_UPPER, M, A, B);

    DEBUG_PRINT_XV(B , "%s", "");
    DEBUG_PRINT_XV(B1, "%s", "");

    matlib_real norm_actual = matlib_xnrm2(B);
    matlib_xaxpy(-1.0, B1, B);
    matlib_real e_relative = matlib_xnrm2(B)/norm_actual;

    debug_body("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);

    matlib_xm fM;
    full_sym_mat(M, &fM);
    
    matlib_xmwrite_csv("fM.dat", fM);


    matlib_xgemv(1.0, fM, A, 0.0, B);


    norm_actual = matlib_xnrm2(B);
    matlib_xaxpy(-1.0, B1, B);
    e_relative = matlib_xnrm2(B)/norm_actual;

    debug_body("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);
}

/*============================================================================*/
/* Matrix I/O functuons 
 *
 * matlib_io_create
 * matlib_io_< data-type >_write
 * matlib_io_fwrite
 *
 * matlib_io_< data-type >_read
 * matlib_io_getinfo
 * matlib_io_getelem
 * matlib_io_fread
 *
 * matlib_io_elemfree
 * matlib_io_free
 * matlib_io_freeall
 *
 * */ 

/* P -> 8-by-5 matrix */ 
static matlib_complex P1_cm[40] = {
8.147236863931789e-01 + 4.387443596563982e-01*I,
9.057919370756192e-01 + 3.815584570930084e-01*I,
1.269868162935061e-01 + 7.655167881490024e-01*I,
9.133758561390194e-01 + 7.951999011370632e-01*I,
6.323592462254095e-01 + 1.868726045543786e-01*I,
9.754040499940952e-02 + 4.897643957882311e-01*I,
2.784982188670484e-01 + 4.455862007108995e-01*I,
5.468815192049838e-01 + 6.463130101112646e-01*I,
9.575068354342976e-01 + 7.093648308580726e-01*I,
9.648885351992765e-01 + 7.546866819823609e-01*I,
1.576130816775483e-01 + 2.760250769985784e-01*I,
9.705927817606157e-01 + 6.797026768536748e-01*I,
9.571669482429456e-01 + 6.550980039738407e-01*I,
4.853756487228412e-01 + 1.626117351946306e-01*I,
8.002804688888001e-01 + 1.189976815583766e-01*I,
1.418863386272153e-01 + 4.983640519821430e-01*I,
4.217612826262750e-01 + 9.597439585160811e-01*I,
9.157355251890671e-01 + 3.403857266661332e-01*I,
7.922073295595544e-01 + 5.852677509797773e-01*I,
9.594924263929030e-01 + 2.238119394911370e-01*I,
6.557406991565868e-01 + 7.512670593056529e-01*I,
3.571167857418955e-02 + 2.550951154592691e-01*I,
8.491293058687771e-01 + 5.059570516651424e-01*I,
9.339932477575505e-01 + 6.990767226566860e-01*I,
6.787351548577735e-01 + 8.909032525357985e-01*I,
7.577401305783334e-01 + 9.592914252054443e-01*I,
7.431324681249162e-01 + 5.472155299638031e-01*I,
3.922270195341682e-01 + 1.386244428286791e-01*I,
6.554778901775566e-01 + 1.492940055590575e-01*I,
1.711866878115618e-01 + 2.575082541237365e-01*I,
7.060460880196088e-01 + 8.407172559836625e-01*I,
3.183284637742068e-02 + 2.542821789715310e-01*I,
2.769229849608900e-01 + 8.142848260688164e-01*I,
4.617139063115394e-02 + 2.435249687249893e-01*I,
9.713178123584754e-02 + 9.292636231872278e-01*I,
8.234578283272926e-01 + 3.499837659848087e-01*I,
6.948286229758170e-01 + 1.965952504312082e-01*I,
3.170994800608605e-01 + 2.510838579760311e-01*I,
9.502220488383549e-01 + 6.160446761466392e-01*I,
3.444608050290876e-02 + 4.732888489027293e-01*I };

static matlib_complex P1_rm[40] = {
8.147236863931789e-01 + 4.387443596563982e-01 * I,
9.575068354342976e-01 + 7.093648308580726e-01 * I,
4.217612826262750e-01 + 9.597439585160811e-01 * I,
6.787351548577735e-01 + 8.909032525357985e-01 * I,
2.769229849608900e-01 + 8.142848260688164e-01 * I,
9.057919370756192e-01 + 3.815584570930084e-01 * I,
9.648885351992765e-01 + 7.546866819823609e-01 * I,
9.157355251890671e-01 + 3.403857266661332e-01 * I,
7.577401305783334e-01 + 9.592914252054443e-01 * I,
4.617139063115394e-02 + 2.435249687249893e-01 * I,
1.269868162935061e-01 + 7.655167881490024e-01 * I,
1.576130816775483e-01 + 2.760250769985784e-01 * I,
7.922073295595544e-01 + 5.852677509797773e-01 * I,
7.431324681249162e-01 + 5.472155299638031e-01 * I,
9.713178123584754e-02 + 9.292636231872278e-01 * I,
9.133758561390194e-01 + 7.951999011370632e-01 * I,
9.705927817606157e-01 + 6.797026768536748e-01 * I,
9.594924263929030e-01 + 2.238119394911370e-01 * I,
3.922270195341682e-01 + 1.386244428286791e-01 * I,
8.234578283272926e-01 + 3.499837659848087e-01 * I,
6.323592462254095e-01 + 1.868726045543786e-01 * I,
9.571669482429456e-01 + 6.550980039738407e-01 * I,
6.557406991565868e-01 + 7.512670593056529e-01 * I,
6.554778901775566e-01 + 1.492940055590575e-01 * I,
6.948286229758170e-01 + 1.965952504312082e-01 * I,
9.754040499940952e-02 + 4.897643957882311e-01 * I,
4.853756487228412e-01 + 1.626117351946306e-01 * I,
3.571167857418955e-02 + 2.550951154592691e-01 * I,
1.711866878115618e-01 + 2.575082541237365e-01 * I,
3.170994800608605e-01 + 2.510838579760311e-01 * I,
2.784982188670484e-01 + 4.455862007108995e-01 * I,
8.002804688888001e-01 + 1.189976815583766e-01 * I,
8.491293058687771e-01 + 5.059570516651424e-01 * I,
7.060460880196088e-01 + 8.407172559836625e-01 * I,
9.502220488383549e-01 + 6.160446761466392e-01 * I,
5.468815192049838e-01 + 6.463130101112646e-01 * I,
1.418863386272153e-01 + 4.983640519821430e-01 * I,
9.339932477575505e-01 + 6.990767226566860e-01 * I,
3.183284637742068e-02 + 2.542821789715310e-01 * I,
3.444608050290876e-02 + 4.732888489027293e-01 * I};

/* P -> 5-by-7 matrix */ 
static matlib_complex P2_cm[35] = {
3.516595070629968e-01 + 8.258169774895474e-01*I,
8.308286278962909e-01 + 5.383424352600571e-01*I,
5.852640911527243e-01 + 9.961347166268855e-01*I,
5.497236082911395e-01 + 7.817552875318368e-02*I,
9.171936638298100e-01 + 4.426782697754463e-01*I,
2.858390188203735e-01 + 1.066527701805844e-01*I,
7.572002291107213e-01 + 9.618980808550537e-01*I,
7.537290942784953e-01 + 4.634224134067444e-03*I,
3.804458469753567e-01 + 7.749104647115024e-01*I,
5.678216407252211e-01 + 8.173032206534330e-01*I,
7.585428956306361e-02 + 8.686947053635097e-01*I,
5.395011866660715e-02 + 8.443584551091032e-02*I,
5.307975530089727e-01 + 3.997826490988965e-01*I,
7.791672301020112e-01 + 2.598704028506542e-01*I,
9.340106842291830e-01 + 8.000684802243075e-01*I,
1.299062084737301e-01 + 4.314138274635446e-01*I,
5.688236608721927e-01 + 9.106475944295229e-01*I,
4.693906410582058e-01 + 1.818470283028525e-01*I,
1.190206950124140e-02 + 2.638029165219901e-01*I,
3.371226443988815e-01 + 1.455389803847170e-01*I,
1.621823081932428e-01 + 1.360685587086637e-01*I,
7.942845406839070e-01 + 8.692922076400893e-01*I,
3.112150420448049e-01 + 5.797045873655702e-01*I,
5.285331355062127e-01 + 5.498602018363320e-01*I,
1.656487294997809e-01 + 1.449547982237268e-01*I,
6.019819414016365e-01 + 8.530311177218937e-01*I,
2.629712845401443e-01 + 6.220551314850660e-01*I,
6.540790984767823e-01 + 3.509523808922709e-01*I,
6.892145031400078e-01 + 5.132495398670534e-01*I,
7.481515928237095e-01 + 4.018080337519417e-01*I,
4.505415985024978e-01 + 7.596669169084191e-02*I,
8.382137799693257e-02 + 2.399161535536580e-01*I,
2.289769687168188e-01 + 1.233189348351655e-01*I,
9.133373615016696e-01 + 1.839077882824167e-01*I,
1.523780189692230e-01 + 2.399525256649028e-01*I};

static matlib_complex P2_rm[35] = {
3.516595070629968e-01 + 8.258169774895474e-01*I,
2.858390188203735e-01 + 1.066527701805844e-01*I,
7.585428956306361e-02 + 8.686947053635097e-01*I,
1.299062084737301e-01 + 4.314138274635446e-01*I,
1.621823081932428e-01 + 1.360685587086637e-01*I,
6.019819414016365e-01 + 8.530311177218937e-01*I,
4.505415985024978e-01 + 7.596669169084191e-02*I,
8.308286278962909e-01 + 5.383424352600571e-01*I,
7.572002291107213e-01 + 9.618980808550537e-01*I,
5.395011866660715e-02 + 8.443584551091032e-02*I,
5.688236608721927e-01 + 9.106475944295229e-01*I,
7.942845406839070e-01 + 8.692922076400893e-01*I,
2.629712845401443e-01 + 6.220551314850660e-01*I,
8.382137799693257e-02 + 2.399161535536580e-01*I,
5.852640911527243e-01 + 9.961347166268855e-01*I,
7.537290942784953e-01 + 4.634224134067444e-03*I,
5.307975530089727e-01 + 3.997826490988965e-01*I,
4.693906410582058e-01 + 1.818470283028525e-01*I,
3.112150420448049e-01 + 5.797045873655702e-01*I,
6.540790984767823e-01 + 3.509523808922709e-01*I,
2.289769687168188e-01 + 1.233189348351655e-01*I,
5.497236082911395e-01 + 7.817552875318368e-02*I,
3.804458469753567e-01 + 7.749104647115024e-01*I,
7.791672301020112e-01 + 2.598704028506542e-01*I,
1.190206950124140e-02 + 2.638029165219901e-01*I,
5.285331355062127e-01 + 5.498602018363320e-01*I,
6.892145031400078e-01 + 5.132495398670534e-01*I,
9.133373615016696e-01 + 1.839077882824167e-01*I,
9.171936638298100e-01 + 4.426782697754463e-01*I,
5.678216407252211e-01 + 8.173032206534330e-01*I,
9.340106842291830e-01 + 8.000684802243075e-01*I,
3.371226443988815e-01 + 1.455389803847170e-01*I,
1.656487294997809e-01 + 1.449547982237268e-01*I,
7.481515928237095e-01 + 4.018080337519417e-01*I,
1.523780189692230e-01 + 2.399525256649028e-01*I};

static matlib_int N1[10] = {1, 2, 3, -4, 5, 6, -7, 8, -9, 10};

void read_mdata(void)
{
    char file_name[] = "bindata.dat";
    matlib_xm M1 = {.lenc   = 5, 
                    .lenr   = 3,
                    .order  = MATLIB_ROW_MAJOR,
                    .op     = MATLIB_NO_TRANS, 
                    .elem_p = A1_rm};

    matlib_xm M2 = { .lenc   = 3, 
                     .lenr   = 2,
                     .order  = MATLIB_COL_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = A2_cm};

    matlib_xv V1 = { .len    = 15,
                     .type   = MATLIB_COL_VECT, 
                     .elem_p = A1_cm};

    matlib_zm M3 = { .lenc   = 8, 
                     .lenr   = 5, 
                     .order  = MATLIB_ROW_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = P1_rm};

    matlib_zm M4 = { .lenc   = 5, 
                     .lenr   = 7, 
                     .order  = MATLIB_COL_MAJOR,
                     .op     = MATLIB_NO_TRANS, 
                     .elem_p = P2_cm};

    matlib_zv V2 = { .len    = 40,
                     .type   = MATLIB_COL_VECT, 
                     .elem_p = P1_cm};

    matlib_nv V3 = { .len    = 10,
                     .elem_p = N1};
    
    matlib_io_t mp;
    matlib_io_create( 7, &mp);
    mp.data_p[0] = &M1;
    mp.data_p[1] = &M2;
    mp.data_p[2] = &V1;
    mp.data_p[3] = &M3;
    mp.data_p[4] = &M4;
    mp.data_p[5] = &V2;
    mp.data_p[6] = &V3;

    mp.format[0] = MATLIB_DEN_XM;
    mp.format[1] = MATLIB_DEN_XM;
    mp.format[2] = MATLIB_XV;
    mp.format[3] = MATLIB_DEN_ZM;
    mp.format[4] = MATLIB_DEN_ZM;
    mp.format[5] = MATLIB_ZV;
    mp.format[6] = MATLIB_NV;

    matlib_io_fwrite(&mp, file_name);
    matlib_io_free(&mp);

    matlib_xm M;
    matlib_xv V;
    matlib_zm zM;
    matlib_zv W;

    matlib_xv u1, u2;
    matlib_zv w1, w2;
    
    matlib_real norm_actual;
    matlib_real e_relative;



#if 1
    matlib_io_t mp1;
    matlib_index data_index = 0;

    matlib_io_fread(&mp1, file_name);
    
    M = *((matlib_xm*)mp1.data_p[0]);
    DEBUG_PRINT_XM(M, "%s: ", "Data read from bin file");
    MATLIB_M2V(M, u1);
    MATLIB_M2V(M1, u2);
    norm_actual = matlib_xnrm2(u2);
    matlib_xaxpy(-1.0, u2, u1);
    e_relative = matlib_xnrm2(u1)/norm_actual;
    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);

    M = *((matlib_xm*)mp1.data_p[1]);
    DEBUG_PRINT_XM(M, "%s: ", "Data read from bin file");
    MATLIB_M2V(M, u1);
    MATLIB_M2V(M2, u2);
    norm_actual = matlib_xnrm2(u2);
    matlib_xaxpy(-1.0, u2, u1);
    e_relative = matlib_xnrm2(u1)/norm_actual;
    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);

    V = *((matlib_xv*)mp1.data_p[2]);
    DEBUG_PRINT_XV(V, "%s: ", "Data read from bin file");
    norm_actual = matlib_xnrm2(V1);
    matlib_xaxpy(-1.0, V1, V);
    e_relative = matlib_xnrm2(V)/norm_actual;
    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);

    zM = *((matlib_zm*)mp1.data_p[3]);
    DEBUG_PRINT_ZM(zM, "%s: ", "Data read from bin file");
    MATLIB_M2V(zM, w1);
    MATLIB_M2V(M3, w2);
    norm_actual = matlib_znrm2(w2);
    matlib_zaxpy(-1.0, w2, w1);
    e_relative = matlib_znrm2(w1)/norm_actual;

    zM = *((matlib_zm*)mp1.data_p[4]);
    DEBUG_PRINT_ZM(zM, "%s: ", "Data read from bin file");
    MATLIB_M2V(zM, w1);
    MATLIB_M2V(M4, w2);
    norm_actual = matlib_znrm2(w2);
    matlib_zaxpy(-1.0, w2, w1);
    e_relative = matlib_znrm2(w1)/norm_actual;
    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);

    W = *((matlib_zv*)mp1.data_p[5]);
    DEBUG_PRINT_ZV(W, "%s: ", "Data read from bin file");
    norm_actual = matlib_znrm2(V2);
    matlib_zaxpy(-1.0, V2, W);
    e_relative = matlib_znrm2(W)/norm_actual;
    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);

    matlib_nv N = *((matlib_nv*)mp1.data_p[6]);
    DEBUG_PRINT_NV(N, "%s: ", "Data read from bin file");


    FILE* fp = fopen(file_name, "rb");
    matlib_io_t mp2;
    matlib_io_getinfo(&mp2, fp);


    matlib_io_getelem(&mp2, 0, fp);
    M = *((matlib_xm*)mp2.data_p[0]);
    DEBUG_PRINT_XM(M, "%s: ", "Data read from bin file");
    MATLIB_M2V(M, u1);
    MATLIB_M2V(M1, u2);
    norm_actual = matlib_xnrm2(u2);
    matlib_xaxpy(-1.0, u2, u1);
    e_relative = matlib_xnrm2(u1)/norm_actual;
    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);

    matlib_io_getelem(&mp2, 1, fp);
    M = *((matlib_xm*)mp2.data_p[1]);
    DEBUG_PRINT_XM(M, "%s: ", "Data read from bin file");
    MATLIB_M2V(M, u1);
    MATLIB_M2V(M2, u2);
    norm_actual = matlib_xnrm2(u2);
    matlib_xaxpy(-1.0, u2, u1);
    e_relative = matlib_xnrm2(u1)/norm_actual;
    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);

    matlib_io_getelem(&mp2, 2, fp);
    V = *((matlib_xv*)mp2.data_p[2]);
    DEBUG_PRINT_XV(V, "%s: ", "Data read from bin file");
    norm_actual = matlib_xnrm2(V1);
    matlib_xaxpy(-1.0, V1, V);
    e_relative = matlib_xnrm2(V)/norm_actual;
    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);

    matlib_io_getelem(&mp2, 3, fp);
    zM = *((matlib_zm*)mp2.data_p[3]);
    DEBUG_PRINT_ZM(zM, "%s: ", "Data read from bin file");
    MATLIB_M2V(zM, w1);
    MATLIB_M2V(M3, w2);
    norm_actual = matlib_znrm2(w2);
    matlib_zaxpy(-1.0, w2, w1);
    e_relative = matlib_znrm2(w1)/norm_actual;

    matlib_io_getelem(&mp2, 4, fp);
    zM = *((matlib_zm*)mp2.data_p[4]);
    DEBUG_PRINT_ZM(zM, "%s: ", "Data read from bin file");
    MATLIB_M2V(zM, w1);
    MATLIB_M2V(M4, w2);
    norm_actual = matlib_znrm2(w2);
    matlib_zaxpy(-1.0, w2, w1);
    e_relative = matlib_znrm2(w1)/norm_actual;
    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);

    matlib_io_getelem(&mp2, 5, fp);
    W = *((matlib_zv*)mp2.data_p[5]);
    DEBUG_PRINT_ZV(W, "%s: ", "Data read from bin file");
    norm_actual = matlib_znrm2(V2);
    matlib_zaxpy(-1.0, V2, W);
    e_relative = matlib_znrm2(W)/norm_actual;
    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);

    fclose(fp);
    matlib_io_freeall(&mp1);
#endif
}
void test_read(void)
{
    matlib_io_t mp1;
    matlib_io_fread(&mp1, "test1.dat");
    
    matlib_xm M = *((matlib_xm*)mp1.data_p[0]);
    DEBUG_PRINT_XM(M, "%s: ", "Data read from bin file");
    matlib_io_free(&mp1);

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
        //{ "Complex matrices"                     , test_complex},
        //{ "Real matrix-vector multiplication"    , test_xgemv  },
        //{ "Real matrix-vector multiplication1"   , test_xgemv1 },
        //{ "Complex matrix-vector multiplication" , test_zgemv  },
        //{ "Complex matrix-vector multiplication1", test_zgemv1 },
        //{ "Real matrix-matrix multiplication"    , test_xgemm  },
        //{ "Complex matrix-matrix multiplication" , test_zgemm  },
        //{ "gemm 1a", test_xgemm_1a },
        //{ "gemm 1b", test_xgemm_1b },
        //{ "gemm 2a", test_xgemm_2a },
        //{ "gemm 2b", test_xgemm_2b },
        //{ "NaNs", test_NaNs },
        //{ "symm. CSR matrix - vector product"    , test_xcsrsymv},
        {"bin data", read_mdata},
        {"bin data1", test_read},
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

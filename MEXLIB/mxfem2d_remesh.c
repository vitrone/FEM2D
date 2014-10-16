/*============================================================================+/
 * mxfem2d_remesh.c
 *
 * option fields:
 *   algorithm
 *   verbose
 *   check
 *   quiet
 *
 * Algorithm codes: 
 *   0: Divide and Conquer Algorithm (default)
 *   1: Sweepline
 *   2: Incremental
 *   
/+============================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "triangle.h"

#define NF (sizeof(fnames)/sizeof(*fnames))

void mexFunction
(
          int      nlhs, 
          mxArray* plhs[],
          int      nrhs, 
    const mxArray* prhs[]
)
{
    if(nrhs<2) 
    {
        mexErrMsgTxt("At least two inputs required.");
    }
    if(nlhs!=2) 
    {
        mexErrMsgTxt("Two output required.");
    } 

    REAL* plist = (REAL *) mxGetPr(prhs[0]);
    int* tlist  = (int  *) mxGetPr(prhs[1]);
    

    if(mxGetM(prhs[0])!=2)
    {
        mexErrMsgTxt("Size of point list is incorrect.");
    } 
    if(mxGetM(prhs[1])!=3)
    {
        mexErrMsgTxt("Size of triangle list is incorrect.");
    } 

    int nr_points    = (int) mxGetN(prhs[0]);
    int nr_triangles = (int) mxGetN(prhs[1]);
    if(nr_points<3)
    {
        mexErrMsgTxt("At least three nodes are needed.");
    } 
    if(nr_triangles<1)
    {
        mexErrMsgTxt("At least one triangle is needed.");
    } 

    int i;
    /*
    for(i=0; i<nr_triangles; i++)
    {
        printf("%d ",  tlist[3*i + 0]);
        printf("%d ",  tlist[3*i + 1]);
        printf("%d\n", tlist[3*i + 2]);
    }*/

    int numberofpointattributes = 0;
    
    int*  pointmarkerlist    = (int* ) NULL;
    REAL* pointattributelist = (REAL*) NULL;

    struct mesh m;
    struct behavior b;
    
    triangleinit(&m);
    init_behavior(&b);
    b.refine   = 1;
    b.quality  = 1;
    b.minangle = 20;
    b.convex   = 1;
    b.noexact  = 0;
    b.splitseg = 0;
    b.conformdel = 0;

    /* Pass more options */ 
    if(nrhs==3)
    {
        int mxarray_index = 2;
        const char *fnames[] = { "algorithm", 
                                 "verbose"  ,
                                 "docheck"  ,
                                 "quiet"    , 
                                 "maxarea"  , 
                                 "minangle" ,
                                 "steiner"  };
        mxArray *mxInput[NF];
        int j, fnum;

        int nfields = mxGetNumberOfFields(prhs[mxarray_index]);
        if (nfields < NF)
        {
            mexErrMsgTxt("Not enough input fields!");
        }
        for( j=0; j<NF; j++)
        {
            fnum = mxGetFieldNumber(prhs[mxarray_index], fnames[j]);
            if(fnum==-1)
            {
                char my_str[30];
                sprintf(my_str, "Missing field: %s:", fnames[j]);
                mexErrMsgTxt(my_str);
            } 
            else
            {
                mxInput[j] = mxGetFieldByNumber(prhs[mxarray_index], 0, fnum);
            }
        }

        int algorithm = (int) mxGetScalar(mxInput[0]);

        b.verbose  =    (int) mxGetScalar(mxInput[1]);
        b.docheck  =    (int) mxGetScalar(mxInput[2]);
        b.quiet    =    (int) mxGetScalar(mxInput[3]);
        b.maxarea  = (double) mxGetScalar(mxInput[4]);
        b.minangle =    (int) mxGetScalar(mxInput[5]);
        b.steiner  =    (int) mxGetScalar(mxInput[6]);

        b.fixedarea = 1;

        if (b.verbose)
        {
            printf("options:\n");
            printf("  algorithm: %d\n", algorithm);
            printf("    verbose: %d\n", b.verbose);
            printf("    docheck: %d\n", b.docheck);
            printf("      quiet: %d\n", b.quiet);
            printf("    maxarea: %0.9g\n", b.maxarea);
            printf("\n");
        }

        switch(algorithm)
        {
            case 1:
                b.sweepline   = 1;
                b.incremental = 0;
                break; 
            case 2:
                b.sweepline   = 0;
                b.incremental = 1;
                break;
            default :
                b.sweepline   = 0;
                b.incremental = 0;
        }
    }

    m.steinerleft = b.steiner;

    transfernodes( &m, &b, 
                   plist, 
                   pointattributelist, /* NULL */ 
                   pointmarkerlist,    /* NULL */ 
                   nr_points,      
                   numberofpointattributes);

    printf("Node transfer complete!\n");
    /* Triangulate the vertices. */
    // m.hullsize = delaunay(&m, &b);
    

    b.poly   = 0;
    b.usesegments = b.poly || b.refine || b.quality || b.convex;
    b.goodangle   = cos(b.minangle * PI / 180.0);
    if (b.goodangle == 1.0) 
    {
        b.offconstant = 0.0;
    } 
    else 
    {
        b.offconstant = 0.475 * sqrt((1.0 + b.goodangle) / (1.0 - b.goodangle));
    }
    b.goodangle *= b.goodangle;
    if (b.refine && b.noiterationnum) 
    {
        mexErrMsgTxt( "Error: Iterations are need when"
                      " refining a triangulation.\n");
    }

    /* Read and reconstruct a mesh. 
     * */
    int numberoftriangleattributes = 0;

    double* trianglearealist = (double*) NULL;
    int* segmentlist       = (int*) NULL; 
    int* segmentmarkerlist = (int*) NULL;
    REAL* talist           = (REAL*) NULL;
    int numberofsegments   = 0;

    m.hullsize = reconstruct( &m, &b, 
                              tlist,
                              talist, 
                              trianglearealist,
                              nr_triangles, 
                              3,
                              numberoftriangleattributes,
                              segmentlist, 
                              segmentmarkerlist,
                              numberofsegments);
        
    printf("Mesh reconstruction complete with hullsize = %d!\n", m.hullsize);
    m.edges = (3l * m.triangles.items + m.hullsize) / 2l;
    /* 
     * Ensure that no vertex can be mistaken for a triangular bounding 
     * box vertex in insertvertex().                                 
     * */
    m.infvertex1 = (vertex) NULL;
    m.infvertex2 = (vertex) NULL;
    m.infvertex3 = (vertex) NULL;
    m.holes   = 0;
    m.regions = 0;

    if (b.quality && (nr_triangles> 0)) 
    {
        /* Enforce angle and area constraints. 
         * */
        enforcequality(&m, &b);           
    }
    int nr_triangles1 = m.triangles.items;
    int nr_points1    = m.vertices.items;

    plhs[0] = mxCreateDoubleMatrix(2, nr_points1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(3, nr_triangles1, mxREAL);
    
    double* plist1 = mxGetPr(plhs[0]);
    double* tlist1 = mxGetPr(plhs[1]);


    writenodes( &m, &b, &plist1, 
                &pointattributelist,
                &pointmarkerlist);

    int* tlist_tmp = (int*) NULL;
    writeelements( &m, &b, &tlist_tmp, &talist);

    for(i=0; i<nr_triangles1; i++)
    {
        tlist1[3*i + 0] = tlist_tmp[3*i + 0] + 1;
        tlist1[3*i + 1] = tlist_tmp[3*i + 1] + 1;
        tlist1[3*i + 2] = tlist_tmp[3*i + 2] + 1;
    }

    if (!b.quiet)
    {
        statistics(&m, &b);
    }

    if (b.docheck)
    {
        m.checksegments = 1;
        checkmesh(&m, &b);         
        int horrors = checkdelaunay(&m, &b);
        if(horrors>0)
        {
            printf("Nr. of non-delaunay triangles found: %d\n", horrors);
            mexErrMsgTxt("Triangulation failed.");
        } 

    }

    triangledeinit(&m, &b);
    //free(tlist_tmp);
}








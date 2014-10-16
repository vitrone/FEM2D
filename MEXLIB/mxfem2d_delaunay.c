/*============================================================================+/
 * mxfem2d_delaunay.c
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
    if(nrhs<1) 
    {
        mexErrMsgTxt("At least one input required.");
    }
    if(nlhs!=1) 
    {
        mexErrMsgTxt("One output required.");
    } 

    REAL* plist = (REAL *) mxGetPr(prhs[0]);

    if(mxGetM(prhs[0])!=2)
    {
        mexErrMsgTxt("Size of point list incorrect.");
    } 

    int nr_points = (int) mxGetN(prhs[0]);
    if(nr_points<3)
    {
        mexErrMsgTxt("At least three nodes are needed.");
    } 

    int numberofpointattributes = 0;
    
    int*  pointmarkerlist    = (int* ) NULL;
    REAL* pointattributelist = (REAL*) NULL;

    struct mesh m;
    struct behavior b;
    
    triangleinit(&m);
    init_behavior(&b);
    /* Turn on exact arithmetic */
    b.noexact = 0;

    b.convex = 1;
    b.usesegments = b.poly || b.refine || b.quality || b.convex;

    /* Pass more options */ 
    if(nrhs==2)
    {
        int mxarray_index = 1;
        const char *fnames[] = { "algorithm", 
                                 "verbose",
                                 "docheck",
                                 "quiet"};
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

        b.verbose = (int) mxGetScalar(mxInput[1]);
        b.docheck = (int) mxGetScalar(mxInput[2]);
        b.quiet   = (int) mxGetScalar(mxInput[3]);


        if (b.verbose)
        {
            printf("options:\n");
            printf("  algorithm: %d\n", algorithm);
            printf("    verbose: %d\n", b.verbose);
            printf("      check: %d\n", b.docheck);
            printf("      quiet: %d\n", b.quiet);
            printf("\n");
        }
    }

    m.steinerleft = b.steiner;

    transfernodes( &m, &b, 
                   plist, 
                   pointattributelist, /* NULL */ 
                   pointmarkerlist,    /* NULL */ 
                   nr_points,      
                   numberofpointattributes);

    /* Triangulate the vertices. */
    m.hullsize = delaunay(&m, &b);

    /* 
     * Ensure that no vertex can be mistaken for a triangular bounding 
     * box vertex in insertvertex().                                 
     * */
    m.infvertex1 = (vertex) NULL;
    m.infvertex2 = (vertex) NULL;
    m.infvertex3 = (vertex) NULL;

    /* Without a PSLG, there can be no holes or regional attributes or 
     * area constraints. The following are set to zero to avoid an 
     * accidental free() later.                                  
     * */
    m.holes   = 0;
    m.regions = 0;
    
    /* Output */ 
    int i;

    REAL* talist = (REAL*) NULL;
    
    int nr_triangles = m.triangles.items;

    plhs[0] = mxCreateDoubleMatrix(3, nr_triangles, mxREAL);
    
    double* tlist = mxGetPr(plhs[0]);

    numbernodes( &m, &b);

    int* tlist_tmp = (int*) NULL;
    writeelements( &m, &b, &tlist_tmp, &talist);

    for(i=0; i<nr_triangles; i++)
    {
        tlist[3*i + 0] = tlist_tmp[3*i + 0] + 1;
        tlist[3*i + 1] = tlist_tmp[3*i + 1] + 1;
        tlist[3*i + 2] = tlist_tmp[3*i + 2] + 1;
    }
    
    if (!b.quiet)
    {
        statistics(&m, &b);
    }

    if (b.docheck)
    {
        if (b.usesegments) 
        {
            int* segmentlist       = (int*) NULL;
            int* segmentmarkerlist = (int*) NULL;
            int numberofsegments   = 0;

            /* Segments will be introduced next. 
             * */
            formskeleton( &m, &b, 
                          segmentlist,
                          segmentmarkerlist, 
                          numberofsegments);
        }

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

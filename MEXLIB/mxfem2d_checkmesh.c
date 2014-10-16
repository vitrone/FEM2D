/*============================================================================+/
 * mxfem2d_checkmesh.c
 *
 * option fields:
 *   algorithm
 *   verbose
 *   check
 *   quiet
 *
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
/* 
    if(nlhs!=1) 
    {
        mexErrMsgTxt("One output required.");
    } 
 * */ 

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
    b.noexact = 0;
    b.convex  = 1;
    b.usesegments = b.poly || b.refine || b.quality || b.convex;

    /* Pass more options */ 
    if(nrhs==3)
    {
        int mxarray_index = 2;
        const char *fnames[] = { "verbose"  ,
                                 "quiet"    };
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


        b.verbose = (int) mxGetScalar(mxInput[0]);
        b.quiet   = (int) mxGetScalar(mxInput[1]);


        if (b.verbose)
        {
            printf("options:\n");
            printf("    verbose: %d\n", b.verbose);
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


    /* Reconstruct a mesh. 
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
    
    statistics(&m, &b);
    m.checksegments = 1;
    checkmesh(&m, &b);         
    int horrors = checkdelaunay(&m, &b);
    if(horrors>0)
    {
        printf("Nr. of non-delaunay triangles found: %d\n", horrors);
    } 
    triangledeinit(&m, &b);
    //free(tlist_tmp);
}









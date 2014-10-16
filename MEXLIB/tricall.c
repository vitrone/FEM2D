/*****************************************************************************/
/*                                                                           */
/*  (tricall.c)                                                              */
/*                                                                           */
/*  Example program that demonstrates how to call Triangle.                  */
/*                                                                           */
/*  Accompanies Triangle Version 1.6                                         */
/*  July 19, 1996                                                            */
/*                                                                           */
/*  This file is placed in the public domain (but the file that it calls     */
/*  is still copyrighted!) by                                                */
/*  Jonathan Richard Shewchuk                                                */
/*  2360 Woolsey #H                                                          */
/*  Berkeley, California  94705-1927                                         */
/*  jrs@cs.berkeley.edu                                                      */
/*                                                                           */
/*****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "triangle.h"

/*****************************************************************************/
/*                                                                           */
/*  main()   Create and refine a mesh.                                       */
/*                                                                           */
/*****************************************************************************/

void mexFunction
(
          int      nlhs, 
          mxArray* plhs[],
          int      nrhs, 
    const mxArray* prhs[]
)
{
    struct triangulateio in, mid, out;
    int i, j;

    /* get input points. */
    if(nrhs!=1) 
    {
        mexErrMsgTxt("One input required.");
    }
    if(nlhs!=2) 
    {
        mexErrMsgTxt("Two output required.");
    } 

    in.pointlist = (REAL *) mxGetPr(prhs[0]);

    if(mxGetM(prhs[0])!=2)
    {
        mexErrMsgTxt("Size of point list incorrect.");
    } 

    in.numberofpoints = (int) mxGetN(prhs[0]);
    in.numberofpointattributes = 1;
    in.pointattributelist = (REAL *) malloc( in.numberofpoints *
                                             in.numberofpointattributes *
                                             sizeof(REAL));
    
    for(i=0; i<in.numberofpointattributes; i++)
    {
        in.pointattributelist[i] = 0.0;
    }
    in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));

    for(i=0; i<in.numberofpoints; i++)
    {
        in.pointmarkerlist[0] = 0;
    }
    in.numberofsegments = 0;
    in.numberofholes    = 0;
    in.numberofregions  = 1;
    /* Make necessary initializations so that Triangle can return a 
     * triangulation in `mid'.  
     * */
    mid.pointlist             = (REAL *) NULL; /* Not needed if -N switch used. */
    mid.pointattributelist    = (REAL *) NULL;
    mid.pointmarkerlist       = (int  *) NULL; /* Not needed if -N or -B switch used. */
    mid.trianglelist          = (int  *) NULL; /* Not needed if -E switch used. */
    mid.triangleattributelist = (REAL *) NULL;
    mid.neighborlist          = (int  *) NULL; /* Needed only if -n switch used. */
    mid.segmentlist           = (int  *) NULL;
    mid.segmentmarkerlist     = (int  *) NULL;
    mid.edgelist              = (int  *) NULL; /* Needed only if -e switch used. */
    mid.edgemarkerlist        = (int  *) NULL; /* Needed if -e used and -B not used. */

    /* Not needed if -N switch used or number of point attributes is zero: */
    /* Not needed if -E switch used or number of triangle attributes is zero: */
    /* Needed only if segments are output (-p or -c) and -P not used: */
    /* Needed only if segments are output (-p or -c) and -P and -B not used: */

    /* Triangulate the points.  Switches are chosen to read and write a  
     *  - PSLG (p), 
     *  - preserve the convex hull (c), 
     *  - number everything from zero (z), 
     *  - assign a regional attribute to each element (A), 
     *  - and produce an edge list (e), 
     *  - a Voronoi diagram (v), 
     *  - and a triangle neighbor list (n).
     * */

    triangulate("zenV", &in, &mid, (struct triangulateio *)NULL);
    
    printf("triangulation complete!\n");
        
    plhs[0] = mxCreateDoubleMatrix( 2, mid.numberofpoints, mxREAL);
    
    REAL *pointlist = mxGetPr(plhs[0]);
    for(i=0; i<mid.numberofpoints; i++)
    {
        pointlist[2*i]   = mid.pointlist[2*i];
        pointlist[2*i+1] = mid.pointlist[2*i+1];
    }
    printf("triangulation complete!\n");

    plhs[1] = mxCreateDoubleMatrix( 3, mid.numberoftriangles, mxREAL);
    
    REAL *trianglelist = mxGetPr(plhs[1]);
    for(i=0; i<mid.numberoftriangles; i++)
    {
        trianglelist[3*i + 0] = mid.trianglelist[3*i + 0];
        trianglelist[3*i + 1] = mid.trianglelist[3*i + 1];
        trianglelist[3*i + 2] = mid.trianglelist[3*i + 2];
    }
    printf("triangulation complete!\n");

    /* Free all allocated arrays, including those allocated by Triangle. */

    //free(in.pointlist);
    //free(in.pointattributelist);
    //free(in.pointmarkerlist);
    //free(in.regionlist);
    //free(mid.pointlist);
    //free(mid.pointattributelist);
    //free(mid.pointmarkerlist);
    //free(mid.trianglelist);
    //free(mid.triangleattributelist);
    //free(mid.trianglearealist);
    //free(mid.neighborlist);
    //free(mid.segmentlist);
    //free(mid.segmentmarkerlist);
    //free(mid.edgelist);
    //free(mid.edgemarkerlist);
}

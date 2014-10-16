/*============================================================================+/
 *
/+============================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mex.h"
#include "matrix.h"
#include "triangle.h"


void mxtriangulate
(
    struct triangulateio *in,
    struct triangulateio *out, 
)
{
    struct mesh m;
    struct behavior b;
    REAL *holearray;                                        /* Array of holes. */
    REAL *regionarray;   /* Array of regional attributes and area constraints. */
    
    triangleinit(&m);
    parsecommandline(1, &triswitches, &b);
    m.steinerleft = b.steiner;

    transfernodes( &m, &b, in->pointlist, 
                   in->pointattributelist,
                   in->pointmarkerlist, 
                   in->numberofpoints,
                   in->numberofpointattributes);

    #ifdef CDT_ONLY
    m.hullsize = delaunay(&m, &b);                /* Triangulate the vertices. */
    #else /* not CDT_ONLY */
    if (b.refine) 
    {
        /* Read and reconstruct a mesh. 
         * */
        m.hullsize = reconstruct( &m, &b, in->trianglelist,
                                  in->triangleattributelist, 
                                  in->trianglearealist,
                                  in->numberoftriangles, 
                                  in->numberofcorners,
                                  in->numberoftriangleattributes,
                                  in->segmentlist, 
                                  in->segmentmarkerlist,
                                  in->numberofsegments);
    } 
    else 
    {
        /* Triangulate the vertices. 
         * */
        m.hullsize = delaunay(&m, &b);
    }
    #endif /* not CDT_ONLY */

    /* 
     * Ensure that no vertex can be mistaken for a triangular bounding 
     * box vertex in insertvertex().                                 
     * */
    m.infvertex1 = (vertex) NULL;
    m.infvertex2 = (vertex) NULL;
    m.infvertex3 = (vertex) NULL;

    if (b.usesegments) 
    {
        m.checksegments = 1;                
        /* Segments will be introduced next. 
         * */
        if (!b.refine) 
        {
            /* Insert PSLG segments and/or convex hull segments. 
             * */
            formskeleton( &m, &b, in->segmentlist,
                          in->segmentmarkerlist, 
                          in->numberofsegments);
        }
    }
    
    if (b.poly && (m.triangles.items > 0)) 
    {
        holearray   = in->holelist;
        m.holes     = in->numberofholes;
        regionarray = in->regionlist;
        m.regions   = in->numberofregions;

        if (!b.refine) 
        {
            /* Carve out holes and concavities. 
             * */
            carveholes(&m, &b, holearray, m.holes, regionarray, m.regions);
        }
    } 
    else 
    {
        /* Without a PSLG, there can be no holes or regional attributes or 
         * area constraints. The following are set to zero to avoid an 
         * accidental free() later.                                  
         * */
        m.holes = 0;
        m.regions = 0;
    }

    #ifndef CDT_ONLY
    if (b.quality && (m.triangles.items > 0)) 
    {
        /* Enforce angle and area constraints. 
         * */
        enforcequality(&m, &b);           
    }
    #endif /* not CDT_ONLY */
    
    /* Calculate the number of edges. */
    m.edges = (3l * m.triangles.items + m.hullsize) / 2l;

    if (b.order > 1) 
    {
        /* Promote elements to higher polynomial order. 
         * */
        highorder(&m, &b);       
    }
    if (!b.quiet) 
    {
        printf("\n");
    }

    if (b.jettison) 
    {
        out->numberofpoints = m.vertices.items - m.undeads;
    } 
    else 
    {
        out->numberofpoints = m.vertices.items;
    }
    out->numberofpointattributes    = m.nextras;
    out->numberoftriangles          = m.triangles.items;
    out->numberofcorners            = (b.order + 1) * (b.order + 2) / 2;
    out->numberoftriangleattributes = m.eextras;
    out->numberofedges              = m.edges;
    
    if (b.usesegments) 
    {
        out->numberofsegments = m.subsegs.items;
    } 
    else 
    {
        out->numberofsegments = m.hullsize;
    }
    if (vorout != (struct triangulateio *) NULL) 
    {
        vorout->numberofpoints          = m.triangles.items;
        vorout->numberofpointattributes = m.nextras;
        vorout->numberofedges           = m.edges;
    }
    /* If not using iteration numbers, don't write a .node file if one was 
     * read, because the original one would be overwritten!              
     * */
    if (b.nonodewritten || (b.noiterationnum && m.readnodefile)) 
    {
        if (!b.quiet) 
        {
            printf("NOT writing vertices.\n");
        }
        /* We must remember to number the vertices. 
         * */
        numbernodes(&m, &b);         
    } 
    else 
    {
        /* writenodes() numbers the vertices too. 
         * */
        writenodes( &m, &b, &out->pointlist, 
                    &out->pointattributelist,
                    &out->pointmarkerlist);
    }
    if (b.noelewritten) 
    {
        if (!b.quiet) 
        {
            printf("NOT writing triangles.\n");
        }
    } 
    else 
    {
        writeelements(&m, &b, &out->trianglelist, &out->triangleattributelist);
    }
    /* The -c switch (convex switch) causes a PSLG to be written 
     * even if none was read.                                  
     * */
    if (b.poly || b.convex) 
    {
        /* If not using iteration numbers, don't overwrite the .poly file. 
         * */
        if (b.nopolywritten || b.noiterationnum) 
        {
            if (!b.quiet) 
            {
                printf("NOT writing segments.\n");
            }
        } 
        else 
        {
            writepoly(&m, &b, &out->segmentlist, &out->segmentmarkerlist);
            out->numberofholes   = m.holes;
            out->numberofregions = m.regions;
            if (b.poly) 
            {
                out->holelist   = in->holelist;
                out->regionlist = in->regionlist;
            } 
            else 
            {
                out->holelist   = (REAL *) NULL;
                out->regionlist = (REAL *) NULL;
            }
        }
    }
    if (b.edgesout) 
    {
        writeedges(&m, &b, &out->edgelist, &out->edgemarkerlist);
    }
    if (b.voronoi) 
    {
        writevoronoi( &m, &b, &vorout->pointlist, 
                      &vorout->pointattributelist,
                      &vorout->pointmarkerlist, &vorout->edgelist,
                      &vorout->edgemarkerlist, &vorout->normlist);
    }
    if (b.neighbors) 
    {
        writeneighbors(&m, &b, &out->neighborlist);
    }

    if (!b.quiet)
    {
        statistics(&m, &b);
    }

    #ifndef REDUCED
    if (b.docheck) 
    {
        checkmesh(&m, &b);
        checkdelaunay(&m, &b);
    }
    #endif /* not REDUCED */

    triangledeinit(&m, &b);
}






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

    in.pointlist      = (REAL *) mxGetPr(prhs[0]);

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


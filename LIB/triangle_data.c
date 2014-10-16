/*============================================================================+/
 *
/+============================================================================*/

#include "triangle.h"

/*****************************************************************************/
/*                                                                           */
/*  parsecommandline()   Read the command line, identify switches, and set   */
/*                       up options and file names.                          */
/*                                                                           */
/*****************************************************************************/
void init_behavior(struct behavior* b)
{
    b->poly      = 0; /* Polygon or PLSG */ 
    b->refine    = 0;
    b->quality   = 0;
    b->vararea   = 0;
    b->fixedarea = 0;
    b->usertest  = 0;

    b->regionattrib   =  0;
    b->convex         =  0;
    b->weighted       =  0;
    b->jettison       =  0;
    b->firstnumber    =  0;
    b->edgesout       =  0;
    b->voronoi        =  0;
    b->neighbors      =  0;
    b->geomview       =  0;
    b->nobound        =  0;
    b->nopolywritten  =  0;
    b->nonodewritten  =  0;
    b->noelewritten   =  0;
    b->noiterationnum =  0;
    b->noholes        =  0;
    b->noexact        =  0;
    b->incremental    =  0;
    b->sweepline      =  0;
    b->dwyer          =  1;
    b->splitseg       =  0;
    b->nobisect       =  0;
    b->conformdel     =  0;
    b->steiner        = -1;
    b->usesegments    =  0;
    b->order          =  1;

    b->minangle =  0.0;
    b->maxarea  = -1.0;
    
    b->docheck = 0;
    b->quiet   = 1;
    b->verbose = 0;

}

void parsecommandline
(
    int argc, 
    char **argv, 
    struct behavior *b
)
{
    #define STARTINDEX 0
    int i, j, k;
    char workstring[FILENAMESIZE];

    b->poly      = 0;
    b->refine    = 0;
    b->quality   = 0;
    b->vararea   = 0;
    b->fixedarea = 0;
    b->usertest  = 0;

    b->regionattrib   =  0;
    b->convex         =  0;
    b->weighted       =  0;
    b->jettison       =  0;
    b->firstnumber    =  1;
    b->edgesout       =  0;
    b->voronoi        =  0;
    b->neighbors      =  0;
    b->geomview       =  0;
    b->nobound        =  0;
    b->nopolywritten  =  0;
    b->nonodewritten  =  0;
    b->noelewritten   =  0;
    b->noiterationnum =  0;
    b->noholes        =  0;
    b->noexact        =  0;
    b->incremental    =  0;
    b->sweepline      =  0;
    b->dwyer          =  1;
    b->splitseg       =  0;
    b->docheck        =  0;
    b->nobisect       =  0;
    b->conformdel     =  0;
    b->steiner        = -1;
    b->order          =  1;

    b->minangle =  0.0;
    b->maxarea  = -1.0;
    
    b->quiet   = 0;
    b->verbose = 0;

    for (i = STARTINDEX; i < argc; i++) 
    {
        for (j = STARTINDEX; argv[i][j] != '\0'; j++) 
        {
            if (argv[i][j] == 'p') 
            {
                b->poly = 1;
            }
            #ifndef CDT_ONLY
            if (argv[i][j] == 'r') 
            {
                b->refine = 1;
            }
            if (argv[i][j] == 'q') 
            {
                b->quality = 1;
                if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) 
                        || (argv[i][j + 1] == '.')) 
                {
                    k = 0;
                    while ( ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                            (argv[i][j + 1] == '.')) 
                    {
                        j++;
                        workstring[k] = argv[i][j];
                        k++;
                    }
                    workstring[k] = '\0';
                    b->minangle = (REAL) strtod(workstring, (char **) NULL);
                } 
                else 
                {
                    b->minangle = 20.0;
                }
            }
            if (argv[i][j] == 'a') 
            {
                b->quality = 1;
                if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                    (argv[i][j + 1] == '.')) 
                {
                    b->fixedarea = 1;
                    k = 0;
                    while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                            (argv[i][j + 1] == '.')) 
                    {
                        j++;
                        workstring[k] = argv[i][j];
                        k++;
                    }
                    workstring[k] = '\0';
                    b->maxarea = (REAL) strtod(workstring, (char **) NULL);
                    if (b->maxarea <= 0.0) 
                    {
                        printf("Error:  Maximum area must be greater than zero.\n");
                        triexit(1);
	    }
                } 
                else 
                {
                    b->vararea = 1;
                }
            }
                
            if (argv[i][j] == 'u') 
            {
                b->quality = 1;
                b->usertest = 1;
            }
            #endif /* not CDT_ONLY */

            if (argv[i][j] == 'A') 
            {
                b->regionattrib = 1;
            }
            if (argv[i][j] == 'c') 
            {
                b->convex = 1;
            }
            if (argv[i][j] == 'w') 
            {
                b->weighted = 1;
            }
            if (argv[i][j] == 'W') 
            {
                b->weighted = 2;
            }
            if (argv[i][j] == 'j') 
            {
                b->jettison = 1;
            }
            if (argv[i][j] == 'z') 
            {
                b->firstnumber = 0;
            }
            if (argv[i][j] == 'e') 
            {
                b->edgesout = 1;
            }
            if (argv[i][j] == 'v') 
            {
                b->voronoi = 1;
            }
            if (argv[i][j] == 'n') 
            {
                b->neighbors = 1;
            }
            if (argv[i][j] == 'g') 
            {
                b->geomview = 1;
            }
            if (argv[i][j] == 'B') 
            {
              b->nobound = 1;
            }
            if (argv[i][j] == 'P') 
            {
              b->nopolywritten = 1;
            }
            if (argv[i][j] == 'N') 
            {
              b->nonodewritten = 1;
            }
            if (argv[i][j] == 'E') 
            {
              b->noelewritten = 1;
            }
            if (argv[i][j] == 'O')
            {
                b->noholes = 1;
            }
            if (argv[i][j] == 'X') 
            {
                b->noexact = 1;
            }
            if (argv[i][j] == 'o') 
            {
                if (argv[i][j + 1] == '2') 
                {
                    j++;
                    b->order = 2;
                }
            }
            #ifndef CDT_ONLY
            if (argv[i][j] == 'Y') 
            {
                b->nobisect++;
            }
            if (argv[i][j] == 'S') 
            {
                b->steiner = 0;
                while ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) 
                {
                    j++;
                    b->steiner = b->steiner * 10 + (int) (argv[i][j] - '0');
                }
            }
            #endif /* not CDT_ONLY */

            #ifndef REDUCED
            if (argv[i][j] == 'i') 
            {
                b->incremental = 1;
            }
            if (argv[i][j] == 'F') 
            {
                b->sweepline = 1;
            }
            #endif /* not REDUCED */

            if (argv[i][j] == 'l') 
            {
                b->dwyer = 0;
            }
            #ifndef REDUCED
            #ifndef CDT_ONLY
            if (argv[i][j] == 's') 
            {
                b->splitseg = 1;
            }
            if ((argv[i][j] == 'D') || (argv[i][j] == 'L')) 
            {
                b->quality = 1;
                b->conformdel = 1;
            }
            #endif /* not CDT_ONLY */
            
            if (argv[i][j] == 'C') 
            {
                b->docheck = 1;
            }
            #endif /* not REDUCED */

            if (argv[i][j] == 'Q') 
            {
                b->quiet = 1;
            }
            if (argv[i][j] == 'V') 
            {
                b->verbose++;
            }
        }
    }
    b->usesegments = b->poly || b->refine || b->quality || b->convex;
    b->goodangle = cos(b->minangle * PI / 180.0);
    if (b->goodangle == 1.0) 
    {
        b->offconstant = 0.0;
    } 
    else 
    {
        b->offconstant = 0.475 * sqrt((1.0 + b->goodangle) / (1.0 - b->goodangle));
    }
    b->goodangle *= b->goodangle;

    if (b->refine && b->noiterationnum) 
    {
        printf(
            "Error:  You cannot use the -I switch when refining a triangulation.\n");
        triexit(1);
    }
  /* Be careful not to allocate space for element area constraints that */
  /*   will never be assigned any value (other than the default -1.0).  */
  if (!b->refine && !b->poly) {
    b->vararea = 0;
  }
  /* Be careful not to add an extra attribute to each element unless the */
  /*   input supports it (PSLG in, but not refining a preexisting mesh). */
  if (b->refine || !b->poly) {
    b->regionattrib = 0;
  }
  /* Regular/weighted triangulations are incompatible with PSLGs */
  /*   and meshing.                                              */
  if (b->weighted && (b->poly || b->quality)) {
    b->weighted = 0;
    if (!b->quiet) {
      printf("Warning:  weighted triangulations (-w, -W) are incompatible\n");
      printf("  with PSLGs (-p) and meshing (-q, -a, -u).  Weights ignored.\n"
             );
    }
  }
  if (b->jettison && b->nonodewritten && !b->quiet) {
    printf("Warning:  -j and -N switches are somewhat incompatible.\n");
    printf("  If any vertices are jettisoned, you will need the output\n");
    printf("  .node file to reconstruct the new node indices.");
  }

}

/*****************************************************************************/
/*                                                                           */
/*  transfernodes()   Read the vertices from memory.                         */
/*                                                                           */
/*****************************************************************************/

void transfernodes
(
    struct mesh     *m, 
    struct behavior *b, 
    REAL* pointlist,
    REAL* pointattriblist, 
    int*  pointmarkerlist,
    int   numberofpoints, 
    int   numberofpointattribs
)
{
    vertex vertexloop;
    REAL x, y;
    int i, j;
    int coordindex;
    int attribindex;

    m->invertices = numberofpoints;
    m->mesh_dim   = 2;
    m->nextras    = numberofpointattribs;
    m->readnodefile = 0;
    if (m->invertices < 3) 
    {
        printf("Error:  Input must have at least three input vertices.\n");
        triexit(1);
    }
    if (m->nextras == 0) 
    {
        b->weighted = 0;
    }

    initialize_vertex_pool(m, b);

    /* Read the vertices. */
    coordindex  = 0;
    attribindex = 0;
    for (i = 0; i < m->invertices; i++) 
    {
        vertexloop = (vertex) poolalloc(&m->vertices);
        /* Read the vertex coordinates. */
        x = vertexloop[0] = pointlist[coordindex++];
        y = vertexloop[1] = pointlist[coordindex++];
        /* Read the vertex attributes. */
        for (j = 0; j < numberofpointattribs; j++) 
        {
            vertexloop[2 + j] = pointattriblist[attribindex++];
        }
        if (pointmarkerlist != (int *) NULL) 
        {
            /* Read a vertex marker. */
            setvertexmark(vertexloop, pointmarkerlist[i]);
        } 
        else 
        {
            /* If no markers are specified, they default to zero. */
            setvertexmark(vertexloop, 0);
        }
        setvertextype(vertexloop, INPUTVERTEX);
        /* Determine the smallest and largest x and y coordinates. */
        if (i == 0) 
        {
            m->xmin = m->xmax = x;
            m->ymin = m->ymax = y;
        } 
        else 
        {
            m->xmin = (x < m->xmin) ? x : m->xmin;
            m->xmax = (x > m->xmax) ? x : m->xmax;
            m->ymin = (y < m->ymin) ? y : m->ymin;
            m->ymax = (y > m->ymax) ? y : m->ymax;
        }
    }

    /* Nonexistent x value used as a flag to mark circle events in sweepline */
    /*   Delaunay algorithm.                                                 */
    m->xminextreme = 10 * m->xmin - 9 * m->xmax;
}


/*****************************************************************************/
/*                                                                           */
/*  writenodes()   Number the vertices and write them to a .node file.       */
/*                                                                           */
/*  To save memory, the vertex numbers are written over the boundary markers */
/*  after the vertices are written to a file.                                */
/*                                                                           */
/*****************************************************************************/

void writenodes
(
    struct mesh *m, 
    struct behavior *b, 
    REAL** pointlist,
    REAL** pointattriblist, 
    int**  pointmarkerlist
)
/* INPUT: m, b
 * OUTPUT: pointlist, pointattriblist, pointmarkerlist
 *
 * */
{
    REAL *plist;
    REAL *palist;
    int *pmlist;
    int coordindex;
    int attribindex;
    vertex vertexloop;
    long outvertices;
    int vertexnumber;
    int i;

    if (b->jettison) 
    {
        outvertices = m->vertices.items - m->undeads;
    } 
    else 
    {
        outvertices = m->vertices.items;
    }

    if (!b->quiet) 
    {
        printf("Writing vertices.\n");
    }
    /* Allocate memory for output vertices if necessary. */
    if (*pointlist == (REAL *) NULL) 
    {
        if (!b->quiet) 
        {
            printf("Allocating memory for vertices.\n");
        }
        *pointlist = (REAL *) trimalloc((int) (outvertices * 2 * sizeof(REAL)));
    }
    /* Allocate memory for output vertex attributes if necessary. */
    if ((m->nextras > 0) && (*pointattriblist == (REAL *) NULL)) 
    {
        if (!b->quiet) 
        {
            printf("Allocating memory for point attribute list.\n");
        }
        *pointattriblist = (REAL *) trimalloc((int) ( outvertices * m->nextras *
                                                      sizeof(REAL)));
    }
    /* Allocate memory for output vertex markers if necessary. */
    if (!b->nobound && (*pointmarkerlist == (int *) NULL)) 
    {
        *pointmarkerlist = (int *) trimalloc((int) (outvertices * sizeof(int)));
    }
    plist       = *pointlist;
    palist      = *pointattriblist;
    pmlist      = *pointmarkerlist;
    coordindex  = 0;
    attribindex = 0;

    traversalinit(&m->vertices);
    vertexnumber = b->firstnumber;
    vertexloop   = vertextraverse(m);

    while (vertexloop != (vertex) NULL) 
    {
        if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX)) 
        {
            /* X and y coordinates. */
            plist[coordindex++] = vertexloop[0];
            plist[coordindex++] = vertexloop[1];
            /* Vertex attributes. */
            for (i = 0; i < m->nextras; i++) 
            {
                palist[attribindex++] = vertexloop[2 + i];
            }
            if (!b->nobound) 
            {
                /* Copy the boundary marker. */
                pmlist[vertexnumber - b->firstnumber] = vertexmark(vertexloop);
            }

            setvertexmark(vertexloop, vertexnumber);
            vertexnumber++;
        }
        vertexloop = vertextraverse(m);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  numbernodes()   Number the vertices.                                     */
/*                                                                           */
/*  Each vertex is assigned a marker equal to its number.                    */
/*                                                                           */
/*  Used when writenodes() is not called because no .node file is written.   */
/*                                                                           */
/*****************************************************************************/

void numbernodes(struct mesh *m, struct behavior *b)
{
    vertex vertexloop;
    int vertexnumber;

    traversalinit(&m->vertices);
    vertexnumber = b->firstnumber;
    vertexloop = vertextraverse(m);
    while (vertexloop != (vertex) NULL) 
    {
        setvertexmark(vertexloop, vertexnumber);
        if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX)) 
        {
            vertexnumber++;
        }
        vertexloop = vertextraverse(m);
    }
}

/*============================================================================*/

void writeelements
(
    struct mesh *m, 
    struct behavior *b,
    int**  trianglelist, 
    REAL** triangleattriblist
)
{
    int *tlist;
    REAL *talist;
    int vertexindex;
    int attribindex;
    struct otri triangleloop;
    vertex p1, p2, p3;
    vertex mid1, mid2, mid3;
    long elementnumber;
    int i;

    if (!b->quiet) 
    {
        printf("Writing triangles.\n");
    }
    /* Allocate memory for output triangles if necessary. */
    if (*trianglelist == (int *) NULL) 
    {
        if (!b->quiet) 
        {
            printf("Allocating memory for triangles.\n");
        }
        *trianglelist = (int *) trimalloc((int)( 
            m->triangles.items*((b->order + 1)*(b->order + 2)/2)*sizeof(int)));
    }
    /* Allocate memory for output triangle attributes if necessary. */
    if ((m->eextras > 0) && (*triangleattriblist == (REAL *) NULL)) 
    {
        if (!b->quiet) 
        {
            printf("Allocating memory for triangles atttribute list.\n");
        }
        *triangleattriblist = (REAL *) trimalloc((int)(
                    m->triangles.items*m->eextras*sizeof(REAL)));
    }
    tlist  = *trianglelist;
    talist = *triangleattriblist;

    vertexindex = 0;
    attribindex = 0;

    traversalinit(&m->triangles);

    triangleloop.tri    = triangletraverse(m);
    triangleloop.orient = 0;
    elementnumber       = b->firstnumber;
    while (triangleloop.tri != (triangle *) NULL) 
    {
        org(triangleloop, p1);
        dest(triangleloop, p2);
        apex(triangleloop, p3);
        if (b->order == 1) 
        {
            tlist[vertexindex++] = vertexmark(p1);
            tlist[vertexindex++] = vertexmark(p2);
            tlist[vertexindex++] = vertexmark(p3);
        } 
        else 
        {
            mid1 = (vertex) triangleloop.tri[m->highorderindex + 1];
            mid2 = (vertex) triangleloop.tri[m->highorderindex + 2];
            mid3 = (vertex) triangleloop.tri[m->highorderindex];

            tlist[vertexindex++] = vertexmark(p1);
            tlist[vertexindex++] = vertexmark(p2);
            tlist[vertexindex++] = vertexmark(p3);
            tlist[vertexindex++] = vertexmark(mid1);
            tlist[vertexindex++] = vertexmark(mid2);
            tlist[vertexindex++] = vertexmark(mid3);
        }

        for (i = 0; i < m->eextras; i++) 
        {
            talist[attribindex++] = elemattribute(triangleloop, i);
        }

        triangleloop.tri = triangletraverse(m);
        elementnumber++;
    }
}


void writepoly
(
    struct mesh *m, 
    struct behavior *b,
    int**  segmentlist, 
    int**  segmentmarkerlist
)
{
    int *slist;
    int *smlist;
    int index;
    struct osub subsegloop;
    vertex endpoint1, endpoint2;
    long subsegnumber;

    if (!b->quiet) 
    {
        printf("Writing segments.\n");
    }
    /* Allocate memory for output segments if necessary. */
    if (*segmentlist == (int *) NULL) 
    {
        *segmentlist = (int *) trimalloc((int)(m->subsegs.items*2*sizeof(int)));
    }
    /* Allocate memory for output segment markers if necessary. */
    if (!b->nobound && (*segmentmarkerlist == (int *) NULL)) 
    {
        *segmentmarkerlist = (int *) trimalloc((int)(m->subsegs.items*sizeof(int)));
    }
    slist  = *segmentlist;
    smlist = *segmentmarkerlist;
    index  = 0;

    traversalinit(&m->subsegs);

    subsegloop.ss       = subsegtraverse(m);
    subsegloop.ssorient = 0;
    subsegnumber        = b->firstnumber;
    
    while (subsegloop.ss != (subseg *) NULL) 
    {
        sorg(subsegloop, endpoint1);
        sdest(subsegloop, endpoint2);
        /* Copy indices of the segment's two endpoints. */
        slist[index++] = vertexmark(endpoint1);
        slist[index++] = vertexmark(endpoint2);
        if (!b->nobound) 
        {
          /* Copy the boundary marker. */
          smlist[subsegnumber - b->firstnumber] = mark(subsegloop);
        }

        subsegloop.ss = subsegtraverse(m);
        subsegnumber++;
    }

}

/*****************************************************************************/
/*                                                                           */
/*  writeedges()   Write the edges to an .edge file.                         */
/*                                                                           */
/*****************************************************************************/

void writeedges
(
    struct mesh *m, 
    struct behavior *b,
    int** edgelist, 
    int** edgemarkerlist
)
{
    int *elist;
    int *emlist;
    int index;
    struct otri triangleloop, trisym;
    struct osub checkmark;
    vertex p1, p2;
    long edgenumber;
    triangle ptr;                         /* Temporary variable used by sym(). */
    subseg sptr;                      /* Temporary variable used by tspivot(). */

    if (!b->quiet) 
    {
        printf("Writing edges.\n");
    }
    /* Allocate memory for edges if necessary. */
    if (*edgelist == (int *) NULL) 
    {
        *edgelist = (int *) trimalloc((int) (m->edges * 2 * sizeof(int)));
    }
    /* Allocate memory for edge markers if necessary. */
    if (!b->nobound && (*edgemarkerlist == (int *) NULL)) 
    {
        *edgemarkerlist = (int *) trimalloc((int) (m->edges * sizeof(int)));
    }
    elist  = *edgelist;
    emlist = *edgemarkerlist;
    index  = 0;

    traversalinit(&m->triangles);
    triangleloop.tri = triangletraverse(m);
    edgenumber       = b->firstnumber;

    /* To loop over the set of edges, loop over all triangles, and look at   
     * the three edges of each triangle.  If there isn't another triangle  
     * adjacent to the edge, operate on the edge.  If there is another     
     * adjacent triangle, operate on the edge only if the current triangle 
     * has a smaller pointer than its neighbor.  This way, each edge is    
     * considered only once.                                               
     * */
    while (triangleloop.tri != (triangle *) NULL) 
    {
        for ( triangleloop.orient = 0; 
              triangleloop.orient < 3;
              triangleloop.orient++) 
        {
            sym(triangleloop, trisym);
            if ((triangleloop.tri < trisym.tri) || (trisym.tri == m->dummytri)) 
            {
                org(triangleloop, p1);
                dest(triangleloop, p2);
                elist[index++] = vertexmark(p1);
                elist[index++] = vertexmark(p2);
                if (!b->nobound)
                {
                    /* Edge number, indices of two endpoints, and a boundary marker. 
                     * If there's no subsegment, the boundary marker is zero.      */
                    if (b->usesegments) 
                    {
                        tspivot(triangleloop, checkmark);

                        if (checkmark.ss == m->dummysub) 
                        {
                            emlist[edgenumber - b->firstnumber] = 0;
                        } 
                        else 
                        {
                            emlist[edgenumber - b->firstnumber] = mark(checkmark);
                        }
                    } 
                    else 
                    {
                        emlist[edgenumber - b->firstnumber] = trisym.tri == m->dummytri;
                    }
                }
                edgenumber++;
            }
        }
        triangleloop.tri = triangletraverse(m);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  writevoronoi()   Write the Voronoi diagram to a .v.node and .v.edge      */
/*                   file.                                                   */
/*                                                                           */
/*  The Voronoi diagram is the geometric dual of the Delaunay triangulation. */
/*  Hence, the Voronoi vertices are listed by traversing the Delaunay        */
/*  triangles, and the Voronoi edges are listed by traversing the Delaunay   */
/*  edges.                                                                   */
/*                                                                           */
/*  WARNING:  In order to assign numbers to the Voronoi vertices, this       */
/*  procedure messes up the subsegments or the extra nodes of every          */
/*  element.  Hence, you should call this procedure last.                    */
/*                                                                           */
/*****************************************************************************/

void writevoronoi
(
    struct mesh *m, 
    struct behavior *b, 
    REAL** vpointlist,
    REAL** vpointattriblist, 
    int**  vpointmarkerlist,
    int**  vedgelist, 
    int**  vedgemarkerlist, 
    REAL** vnormlist
)
{
    REAL *plist;
    REAL *palist;
    int *elist;
    REAL *normlist;
    int coordindex;
    int attribindex;
    struct otri triangleloop, trisym;
    vertex torg, tdest, tapex;
    REAL circumcenter[2];
    REAL xi, eta;
    long vnodenumber, vedgenumber;
    int p1, p2;
    int i;
    /* Temporary variable used by sym(). */
    triangle ptr;                         

    if (!b->quiet) 
    {
        printf("Writing Voronoi vertices.\n");
    }
    /* Allocate memory for Voronoi vertices if necessary. */
    if (*vpointlist == (REAL *) NULL) 
    {
        *vpointlist = (REAL *) trimalloc((int)( m->triangles.items * 2 *
                                                sizeof(REAL)));
    }
    /* Allocate memory for Voronoi vertex attributes if necessary. */
    if (*vpointattriblist == (REAL *) NULL) 
    {
        *vpointattriblist = (REAL *) trimalloc((int)( m->triangles.items *
                                                      m->nextras * sizeof(REAL)));
    }
    *vpointmarkerlist = (int *) NULL;
    
    plist  = *vpointlist;
    palist = *vpointattriblist;

    coordindex  = 0;
    attribindex = 0;

    traversalinit(&m->triangles);
    
    triangleloop.tri    = triangletraverse(m);
    triangleloop.orient = 0;
    vnodenumber         = b->firstnumber;

    while (triangleloop.tri != (triangle *) NULL) 
    {
        org(triangleloop, torg);
        dest(triangleloop, tdest);
        apex(triangleloop, tapex);
        findcircumcenter(m, b, torg, tdest, tapex, circumcenter, &xi, &eta, 0);
        /* X and y coordinates. */
        plist[coordindex++] = circumcenter[0];
        plist[coordindex++] = circumcenter[1];
        for (i = 2; i < 2 + m->nextras; i++) 
        {
            /* Interpolate the vertex attributes at the circumcenter. */
            palist[attribindex++] = torg[i] + xi * (tdest[i] - torg[i])
                                       + eta * (tapex[i] - torg[i]);
        }

        * (int *) (triangleloop.tri + 6) = (int) vnodenumber;
        triangleloop.tri = triangletraverse(m);
        vnodenumber++;
    }


    if (!b->quiet) 
    {
        printf("Writing Voronoi edges.\n");
    }
    /* Allocate memory for output Voronoi edges if necessary. */
    if (*vedgelist == (int *) NULL) 
    {
        *vedgelist = (int *) trimalloc((int) (m->edges * 2 * sizeof(int)));
    }
    *vedgemarkerlist = (int *) NULL;
    /* Allocate memory for output Voronoi norms if necessary. */
    if (*vnormlist == (REAL *) NULL) 
    {
        *vnormlist = (REAL *) trimalloc((int) (m->edges * 2 * sizeof(REAL)));
    }
    elist      = *vedgelist;
    normlist   = *vnormlist;
    coordindex = 0;

    traversalinit(&m->triangles);

    triangleloop.tri = triangletraverse(m);
    vedgenumber      = b->firstnumber;

    /* To loop over the set of edges, loop over all triangles, and look at   
     * the three edges of each triangle.  If there isn't another triangle  
     * adjacent to the edge, operate on the edge.  If there is another     
     * adjacent triangle, operate on the edge only if the current triangle 
     * has a smaller pointer than its neighbor.  This way, each edge is    
     * considered only once.                                               
     * */
    while (triangleloop.tri != (triangle *) NULL) 
    {
        for ( triangleloop.orient = 0; 
              triangleloop.orient < 3;
              triangleloop.orient++) 
        {
            sym(triangleloop, trisym);
            if ((triangleloop.tri < trisym.tri) || (trisym.tri == m->dummytri)) 
            {
                /* Find the number of this triangle (and Voronoi vertex). 
                 * */
                p1 = * (int *) (triangleloop.tri + 6);
                if (trisym.tri == m->dummytri) 
                {
                    org(triangleloop, torg);
                    dest(triangleloop, tdest);
                    /* Copy an infinite ray.  Index of one endpoint, and -1. 
                     * */
                    elist[coordindex] = p1;
                    normlist[coordindex++] = tdest[1] - torg[1];
                    elist[coordindex] = -1;
                    normlist[coordindex++] = torg[0] - tdest[0];
                } 
                else 
                {
                    /* Find the number of the adjacent triangle (and Voronoi vertex). 
                     * */
                    p2 = * (int *) (trisym.tri + 6);
                    /* Finite edge.  Write indices of two endpoints. 
                     * */
                    elist[coordindex]      = p1;
                    normlist[coordindex++] = 0.0;
                    elist[coordindex]      = p2;
                    normlist[coordindex++] = 0.0;
                }
                vedgenumber++;
            }
        }
        triangleloop.tri = triangletraverse(m);
    }
}
/*============================================================================*/

void writeneighbors
(
    struct mesh *m, 
    struct behavior *b, 
    int **neighborlist
)
{
    int *nlist;
    int index;
    struct otri triangleloop, trisym;
    long elementnumber;
    int neighbor1, neighbor2, neighbor3;

    /* Temporary variable used by sym(). */
    triangle ptr;                         

    if (!b->quiet) 
    {
        printf("Writing neighbors.\n");
    }
    /* Allocate memory for neighbors if necessary. */
    if (*neighborlist == (int *) NULL) 
    {
        *neighborlist = (int *) trimalloc((int)( m->triangles.items * 3 *
                                                 sizeof(int)));
    }
    nlist = *neighborlist;
    index = 0;

    traversalinit(&m->triangles);

    triangleloop.tri    = triangletraverse(m);
    triangleloop.orient = 0;
    elementnumber       = b->firstnumber;
    
    while (triangleloop.tri != (triangle *) NULL) 
    {
        *(int *) (triangleloop.tri + 6) = (int) elementnumber;
        triangleloop.tri = triangletraverse(m);
        elementnumber++;
    }
    * (int *) (m->dummytri + 6) = -1;

    traversalinit(&m->triangles);
    triangleloop.tri = triangletraverse(m);
    elementnumber    = b->firstnumber;
    while (triangleloop.tri != (triangle *) NULL) 
    {
        triangleloop.orient = 1;
        sym(triangleloop, trisym);
        neighbor1 = * (int *) (trisym.tri + 6);
        triangleloop.orient = 2;
        sym(triangleloop, trisym);
        neighbor2 = * (int *) (trisym.tri + 6);
        triangleloop.orient = 0;
        sym(triangleloop, trisym);
        neighbor3 = * (int *) (trisym.tri + 6);

        nlist[index++] = neighbor1;
        nlist[index++] = neighbor2;
        nlist[index++] = neighbor3;

        triangleloop.tri = triangletraverse(m);
        elementnumber++;
    }

}

/*============================================================================*/

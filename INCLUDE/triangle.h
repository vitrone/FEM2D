#ifndef TRIANGLE_H
#define TRIANGLE_H

/*============================================================================+/
 |
 | DT: A Delaunay triangulation (DT) of a vertex set is a triangulation of the 
 | vertex set with the property that no vertex in the vertex set falls in the 
 | interior of the circumcircle of any triangle in the triangulation.
 |
 | PSLG: A Planar Straight Line Graph (PSLG) is a collection of vertices and 
 | segments. Segments are edges whose endpoints are vertices in the PSLG, and 
 | whose presence in any mesh generated from the PSLG is enforced.
 |
 | CDT: A conforming Delaunay triangulation (CDT) of a PSLG is a true Delaunay
 | triangulation in which each PSLG segment may have been subdivided into 
 | several edges by the insertion of additional vertices, called Steiner points. 
 | Steiner points are necessary to allow the segments to exist in the mesh while
 | maintaining the Delaunay property. Steiner points are also inserted to meet 
 | constraints on the minimum angle and maximum triangle area.
 |
 | _CDT: A constrained Delaunay triangulation of a PSLG is similar to a Delaunay 
 | triangulation, but each PSLG segment is present as a single edge in the 
 | triangulation. A constrained Delaunay triangulation is not truly a Delaunay 
 | triangulation. Some of its triangles might not be Delaunay, but they are all 
 | constrained Delaunay.
 |
 | CCDT: A constrained conforming Delaunay triangulation (CCDT) of a PSLG is a 
 | constrained Delaunay triangulation that includes Steiner points. It usually 
 | takes fewer vertices to make a good-quality CCDT than a good-quality CDT, 
 | because the triangles do not need to be Delaunay (although they still must 
 | be constrained Delaunay).
/+============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LINUX

#ifdef CPU86
#include <float.h>
#endif /* CPU86 */
#ifdef LINUX
#include <fpu_control.h>
#endif /* LINUX */
#include "triangle.h"


#define VOID void
#define REAL double
#define tri_real double

/* To insert lots of self-checks for internal errors, define the SELF_CHECK 
 * symbol.  This will slow down the program significantly.  It is best to 
 * define the symbol using the -DSELF_CHECK compiler switch, but you could 
 * write "#define SELF_CHECK" below.  If you are modifying this code, I 
 * recommend you turn self-checks on until your work is debugged.
 * */

#define SELF_CHECK

/* #define REDUCED */
/* #define CDT_ONLY */


//#define INEXACT /* Nothing */
#define INEXACT volatile

#define FILENAMESIZE 2048

/* For efficiency, a variety of data structures are allocated in bulk.  The  */
/*   following constants determine how many of each structure is allocated   */
/*   at once.                                                                */

#define TRIPERBLOCK 4092         /* Number of triangles allocated at once.   */
#define SUBSEGPERBLOCK 508       /* Number of subsegments allocated at once. */
#define VERTEXPERBLOCK 4092      /* Number of vertices allocated at once.    */
#define VIRUSPERBLOCK 1020       /* Number of virus triangles allocated at once. */

/* Number of encroached subsegments allocated at once. */
#define BADSUBSEGPERBLOCK 252
/* Number of skinny triangles allocated at once. */
#define BADTRIPERBLOCK 4092
/* Number of flipped triangles allocated at once. */
#define FLIPSTACKERPERBLOCK 252
/* Number of splay tree nodes allocated at once. */
#define SPLAYNODEPERBLOCK 508

/* The vertex types. A DEADVERTEX has been deleted entirely. 
 * An UNDEADVERTEX is not part of the mesh, but is written to the output 
 * .node file and affects the node indexing in the other output files.
 * */

#define INPUTVERTEX 0
#define SEGMENTVERTEX 1
#define FREEVERTEX 2
#define DEADVERTEX -32768
#define UNDEADVERTEX -32767


/* Two constants for algorithms based on random sampling.  Both constants    */
/*   have been chosen empirically to optimize their respective algorithms.   */

/* Used for the point location scheme of Mucke, Saias, and Zhu, to decide    */
/*   how large a random sample of triangles to inspect.                      */

#define SAMPLEFACTOR 11

/* Used in Fortune's sweepline Delaunay algorithm to determine what fraction */
/*   of boundary edges should be maintained in the splay tree for point      */
/*   location on the front.                                                  */

#define SAMPLERATE 10

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308
#define SQUAREROOTTWO 1.4142135623730950488016887242096980785696718753769480732
#define ONETHIRD 0.333333333333333333333333333333333333333333333333333333333333


/* Labels that signify the result of point location. The result of a search 
 * indicates that the point falls in the interior of a triangle, on an edge, 
 * on a vertex, or outside the mesh.                              
 * */
enum locateresult 
{
    INTRIANGLE, 
    ONEDGE, 
    ONVERTEX, 
    OUTSIDE
};

/* Labels that signify the result of vertex insertion.  The result indicates */
/*   that the vertex was inserted with complete success, was inserted but    */
/*   encroaches upon a subsegment, was not inserted because it lies on a     */
/*   segment, or was not inserted because another vertex occupies the same   */
/*   location.                                                               */

enum insertvertexresult 
{
    SUCCESSFULVERTEX, 
    ENCROACHINGVERTEX, 
    VIOLATINGVERTEX,                     
    DUPLICATEVERTEX
};

/* Labels that signify the result of direction finding.  The result          */
/*   indicates that a segment connecting the two query points falls within   */
/*   the direction triangle, along the left edge of the direction triangle,  */
/*   or along the right edge of the direction triangle.                      */

enum finddirectionresult 
{
    WITHIN, 
    LEFTCOLLINEAR, 
    RIGHTCOLLINEAR
};

/*============================================================================+/
 | The basic mesh data structures
 |
 | There are three:  
 | - vertices, 
 | - triangles, and 
 | - subsegments (abbreviated `subseg')
 |
 | These three data structures, linked by pointers, comprise the mesh.  
 | 
 | A vertex simply represents a mesh vertex and its properties.
 | 
 | A triangle is a triangle.
 |
 | A subsegment is a special data structure used to represent an impenetrable 
 | edge of the mesh (perhaps on the outer boundary, on the boundary of a 
 | hole, or part of an internal boundary separating two triangulated 
 | regions). Subsegments represent boundaries, defined by the user, that 
 | triangles may not lie across.
 |
 | A triangle consists of a list of three vertices, a list of three adjoining 
 | triangles, a list of three adjoining subsegments (when segments exist), 
 | an arbitrary number of optional user-defined floating-point attributes, 
 | and an optional area constraint. The latter is an upper bound on the 
 | permissible area of each triangle in a region, used for mesh refinement.
 |
 | For a triangle on a boundary of the mesh, some or all of the neighboring 
 | triangles may not be present. For a triangle in the interior of the mesh, 
 | often no neighboring subsegments are present. Such absent triangles and 
 | subsegments are never represented by NULL pointers; they are represented 
 | by two special records:  
 | - `dummytri', the triangle that fills "outer space", and 
 | - `dummysub', the omnipresent subsegment.
 |
 | `dummytri' and `dummysub' are used for several reasons; for instance, 
 | they can be dereferenced and their contents examined without violating 
 | protected memory.                                                        
                                                                           
  However, it is important to understand that a triangle includes other    
  information as well. The pointers to adjoining vertices, triangles, and 
  subsegments are ordered in a way that indicates their geometric relation 
  to each other.  Furthermore, each of these pointers contains orientation 
  information.  Each pointer to an adjoining triangle indicates which face 
  of that triangle is contacted.  Similarly, each pointer to an adjoining  
  subsegment indicates which side of that subsegment is contacted, and how 
  the subsegment is oriented relative to the triangle.                     
                                                                           
  The data structure representing a subsegment may be thought to be        
  abutting the edge of one or two triangle data structures:  either        
  sandwiched between two triangles, or resting against one triangle on an  
  exterior boundary or hole boundary.                                      
                                                                           
  A subsegment consists of a list of four vertices--the vertices of the    
  subsegment, and the vertices of the segment it is a part of--a list of   
  two adjoining subsegments, and a list of two adjoining triangles.  One   
  of the two adjoining triangles may not be present (though there should   
  always be one), and neighboring subsegments might not be present.        
  Subsegments also store a user-defined integer "boundary marker".         
  Typically, this integer is used to indicate what boundary conditions are 
  to be applied at that location in a finite element simulation.           
                                                                           
  Like triangles, subsegments maintain information about the relative      
  orientation of neighboring objects.                                      
                                                                           
  Vertices are relatively simple.  A vertex is a list of floating-point    
  numbers, starting with the x, and y coordinates, followed by an          
  arbitrary number of optional user-defined floating-point attributes,     
  followed by an integer boundary marker.  During the segment insertion    
  phase, there is also a pointer from each vertex to a triangle that may   
  contain it.  Each pointer is not always correct, but when one is, it     
  speeds up segment insertion.  These pointers are assigned values once    
  at the beginning of the segment insertion phase, and are not used or     
  updated except during this phase.  Edge flipping during segment          
  insertion will render some of them incorrect.  Hence, don't rely upon    
  them for anything.                                                       
                                                                           
  Other than the exception mentioned above, vertices have no information   
  about what triangles, subfacets, or subsegments they are linked to.
/+============================================================================*/

/*****************************************************************************/
/*                                                                           */
/*  Handles                                                                  */
/*                                                                           */
/*  The oriented triangle (`otri') and oriented subsegment (`osub') data     */
/*  structures defined below do not themselves store any part of the mesh.   */
/*  The mesh itself is made of `triangle's, `subseg's, and `vertex's.        */
/*                                                                           */
/*  Oriented triangles and oriented subsegments will usually be referred to  */
/*  as "handles."  A handle is essentially a pointer into the mesh; it       */
/*  allows you to "hold" one particular part of the mesh.  Handles are used  */
/*  to specify the regions in which one is traversing and modifying the mesh.*/
/*  A single `triangle' may be held by many handles, or none at all.  (The   */
/*  latter case is not a memory leak, because the triangle is still          */
/*  connected to other triangles in the mesh.)                               */
/*                                                                           */
/*  An `otri' is a handle that holds a triangle.  It holds a specific edge   */
/*  of the triangle.  An `osub' is a handle that holds a subsegment.  It     */
/*  holds either the left or right side of the subsegment.                   */
/*                                                                           */
/*  Navigation about the mesh is accomplished through a set of mesh          */
/*  manipulation primitives, further below.  Many of these primitives take   */
/*  a handle and produce a new handle that holds the mesh near the first     */
/*  handle.  Other primitives take two handles and glue the corresponding    */
/*  parts of the mesh together.  The orientation of the handles is           */
/*  important.  For instance, when two triangles are glued together by the   */
/*  bond() primitive, they are glued at the edges on which the handles lie.  */
/*                                                                           */
/*  Because vertices have no information about which triangles they are      */
/*  attached to, I commonly represent a vertex by use of a handle whose      */
/*  origin is the vertex.  A single handle can simultaneously represent a    */
/*  triangle, an edge, and a vertex.                                         */
/*                                                                           */
/*****************************************************************************/

/* The triangle data structure.  Each triangle contains three pointers to    */
/*   adjoining triangles, plus three pointers to vertices, plus three        */
/*   pointers to subsegments (declared below; these pointers are usually     */
/*   `dummysub').  It may or may not also contain user-defined attributes    */
/*   and/or a floating-point "area constraint."  It may also contain extra   */
/*   pointers for nodes, when the user asks for high-order elements.         */
/*   Because the size and structure of a `triangle' is not decided until     */
/*   runtime, I haven't simply declared the type `triangle' as a struct.     */

typedef REAL** triangle;            /* Really:  typedef triangle *triangle   */

/* An oriented triangle:  includes a pointer to a triangle and orientation.  */
/*   The orientation denotes an edge of the triangle.  Hence, there are      */
/*   three possible orientations.  By convention, each edge always points    */
/*   counterclockwise about the corresponding triangle.                      */

struct otri 
{
    triangle* tri;
    int       orient; /* Ranges from 0 to 2. */
};

/* The subsegment data structure.  Each subsegment contains two pointers to  */
/*   adjoining subsegments, plus four pointers to vertices, plus two         */
/*   pointers to adjoining triangles, plus one boundary marker, plus one     */
/*   segment number.                                                         */

typedef REAL** subseg;                  /* Really:  typedef subseg *subseg   */

/* An oriented subsegment:  includes a pointer to a subsegment and an        */
/*   orientation.  The orientation denotes a side of the edge.  Hence, there */
/*   are two possible orientations.  By convention, the edge is always       */
/*   directed so that the "side" denoted is the right side of the edge.      */

struct osub 
{
  subseg* ss;
  int     ssorient;                                   /* Ranges from 0 to 1. */
};

/* The vertex data structure.  Each vertex is actually an array of REALs.    */
/*   The number of REALs is unknown until runtime.  An integer boundary      */
/*   marker, and sometimes a pointer to a triangle, is appended after the    */
/*   REALs.                                                                  */

typedef REAL *vertex;

/* A queue used to store encroached subsegments.  Each subsegment's vertices */
/*   are stored so that we can check whether a subsegment is still the same. */

struct badsubseg {
  subseg encsubseg;                             /* An encroached subsegment. */
  vertex subsegorg, subsegdest;                         /* Its two vertices. */
};

/* A queue used to store bad triangles.  The key is the square of the cosine */
/*   of the smallest angle of the triangle.  Each triangle's vertices are    */
/*   stored so that one can check whether a triangle is still the same.      */

struct badtriang {
  triangle poortri;                       /* A skinny or too-large triangle. */
  REAL key;                             /* cos^2 of smallest (apical) angle. */
  vertex triangorg, triangdest, triangapex;           /* Its three vertices. */
  struct badtriang *nexttriang;             /* Pointer to next bad triangle. */
};

/* A stack of triangles flipped during the most recent vertex insertion.     */
/*   The stack is used to undo the vertex insertion if the vertex encroaches */
/*   upon a subsegment.                                                      */

struct flipstacker {
  triangle flippedtri;                       /* A recently flipped triangle. */
  struct flipstacker *prevflip;               /* Previous flip in the stack. */
};

/* A node in a heap used to store events for the sweepline Delaunay          */
/*   algorithm.  Nodes do not point directly to their parents or children in */
/*   the heap.  Instead, each node knows its position in the heap, and can   */
/*   look up its parent and children in a separate array.  The `eventptr'    */
/*   points either to a `vertex' or to a triangle (in encoded format, so     */
/*   that an orientation is included).  In the latter case, the origin of    */
/*   the oriented triangle is the apex of a "circle event" of the sweepline  */
/*   algorithm.  To distinguish site events from circle events, all circle   */
/*   events are given an invalid (smaller than `xmin') x-coordinate `xkey'.  */

struct event {
  REAL xkey, ykey;                              /* Coordinates of the event. */
  VOID *eventptr;      /* Can be a vertex or the location of a circle event. */
  int heapposition;              /* Marks this event's position in the heap. */
};

/* A node in the splay tree.  Each node holds an oriented ghost triangle     */
/*   that represents a boundary edge of the growing triangulation.  When a   */
/*   circle event covers two boundary edges with a triangle, so that they    */
/*   are no longer boundary edges, those edges are not immediately deleted   */
/*   from the tree; rather, they are lazily deleted when they are next       */
/*   encountered.  (Since only a random sample of boundary edges are kept    */
/*   in the tree, lazy deletion is faster.)  `keydest' is used to verify     */
/*   that a triangle is still the same as when it entered the splay tree; if */
/*   it has been rotated (due to a circle event), it no longer represents a  */
/*   boundary edge and should be deleted.                                    */

struct splaynode {
  struct otri keyedge;                     /* Lprev of an edge on the front. */
  vertex keydest;           /* Used to verify that splay node is still live. */
  struct splaynode *lchild, *rchild;              /* Children in splay tree. */
};

/* A type used to allocate memory.  firstblock is the first block of items.  */
/*   nowblock is the block from which items are currently being allocated.   */
/*   nextitem points to the next slab of free memory for an item.            */
/*   deaditemstack is the head of a linked list (stack) of deallocated items */
/*   that can be recycled.  unallocateditems is the number of items that     */
/*   remain to be allocated from nowblock.                                   */
/*                                                                           */
/* Traversal is the process of walking through the entire list of items, and */
/*   is separate from allocation.  Note that a traversal will visit items on */
/*   the "deaditemstack" stack as well as live items.  pathblock points to   */
/*   the block currently being traversed.  pathitem points to the next item  */
/*   to be traversed.  pathitemsleft is the number of items that remain to   */
/*   be traversed in pathblock.                                              */
/*                                                                           */
/* alignbytes determines how new records should be aligned in memory.        */
/*   itembytes is the length of a record in bytes (after rounding up).       */
/*   itemsperblock is the number of items allocated at once in a single      */
/*   block.  itemsfirstblock is the number of items in the first block,      */
/*   which can vary from the others.  items is the number of currently       */
/*   allocated items.  maxitems is the maximum number of items that have     */
/*   been allocated at once; it is the current number of items plus the      */
/*   number of records kept on deaditemstack.                                */

struct memorypool 
{
    VOID **firstblock, **nowblock;
    VOID *nextitem;
    VOID *deaditemstack;
    VOID **pathblock;
    VOID *pathitem;
    int alignbytes;
    int itembytes;
    int itemsperblock;
    int itemsfirstblock;
    long items, maxitems;
    int unallocateditems;
    int pathitemsleft;
};


/* Global constants.                                                         */

REAL splitter;       /* Used to split REAL factors for exact multiplication. */
REAL epsilon;                             /* Floating-point machine epsilon. */
REAL resulterrbound;
REAL ccwerrboundA, ccwerrboundB, ccwerrboundC;
REAL iccerrboundA, iccerrboundB, iccerrboundC;
REAL o3derrboundA, o3derrboundB, o3derrboundC;

/* Random number seed is not constant, but I've made it global anyway.       */

unsigned long randomseed;                     /* Current random number seed. */


/* Mesh data structure.  Triangle operates on only one mesh, but the mesh    */
/*   structure is used (instead of global variables) to allow reentrancy.    */

struct mesh 
{

    /* Variables used to allocate memory for triangles, subsegments, vertices,
     * viri (triangles being eaten), encroached segments, bad (skinny or too large) 
     * triangles, and splay tree nodes.
     * */

    struct memorypool triangles;
    struct memorypool subsegs;
    struct memorypool vertices;
    struct memorypool viri;
    struct memorypool badsubsegs;
    struct memorypool badtriangles;
    struct memorypool flipstackers;
    struct memorypool splaynodes;

    /* Variables that maintain the bad triangle queues.  The queues are 
     * ordered from 4095 (highest priority) to 0 (lowest priority).            
     * */

    struct badtriang *queuefront[4096];
    struct badtriang *queuetail[4096];
    int nextnonemptyq[4096];
    int firstnonemptyq;

    /* Variable that maintains the stack of recently flipped triangles.          
     * */

    struct flipstacker *lastflip;

    /* Other variables. */

    REAL xmin, xmax, ymin, ymax;                            /* x and y bounds. */
    REAL xminextreme;      /* Nonexistent x value used as a flag in sweepline. */
    int invertices;                               /* Number of input vertices. */
    int inelements;                              /* Number of input triangles. */
    int insegments;                               /* Number of input segments. */
    int holes;                                       /* Number of input holes. */
    int regions;                                   /* Number of input regions. */
    int undeads;    /* Number of input vertices that don't appear in the mesh. */
    long edges;                                     /* Number of output edges. */
    int mesh_dim;                                /* Dimension (ought to be 2). */
    int nextras;                           /* Number of attributes per vertex. */
    int eextras;                         /* Number of attributes per triangle. */
    long hullsize;                          /* Number of edges in convex hull. */
    int steinerleft;                 /* Number of Steiner points not yet used. */
    int vertexmarkindex;         /* Index to find boundary marker of a vertex. */
    int vertex2triindex;     /* Index to find a triangle adjacent to a vertex. */
    int highorderindex;  /* Index to find extra nodes for high-order elements. */
    int elemattribindex;            /* Index to find attributes of a triangle. */
    int areaboundindex;             /* Index to find area bound of a triangle. */
    int checksegments;         /* Are there segments in the triangulation yet? */
    int checkquality;                  /* Has quality triangulation begun yet? */
    int readnodefile;                           /* Has a .node file been read? */
    long samples;              /* Number of random samples for point location. */

    long incirclecount;                 /* Number of incircle tests performed. */
    long counterclockcount;     /* Number of counterclockwise tests performed. */
    long orient3dcount;           /* Number of 3D orientation tests performed. */
    long hyperbolacount;      /* Number of right-of-hyperbola tests performed. */
    long circumcentercount;  /* Number of circumcenter calculations performed. */
    long circletopcount;       /* Number of circle top calculations performed. */

    /* Triangular bounding box vertices.                                         
     * */
    vertex infvertex1, infvertex2, infvertex3;

    /* Pointer to the `triangle' that occupies all of "outer space."             
     * */

    triangle *dummytri;

    /* Keep base address so we can free() it later. */
    triangle *dummytribase;    

    /* Pointer to the omnipresent subsegment.  Referenced by any triangle or
     * subsegment that isn't really connected to a subsegment at that
     * location.                                                               
     * */
  
    subseg *dummysub;
    subseg *dummysubbase;      /* Keep base address so we can free() it later. */

    /* Pointer to a recently visited triangle. Improves point location if proximate 
     * vertices are inserted sequentially.                           
     * */
    struct otri recenttri;

};                                                  /* End of `struct mesh'. */


/* Data structure for command line switches and file names.  This structure  */
/*   is used (instead of global variables) to allow reentrancy.              */
/* Switches for the triangulator.                                            */
/*   poly: -p switch.  refine: -r switch.                                    */
/*   quality: -q switch.                                                     */
/*     minangle: minimum angle bound, specified after -q switch.             */
/*     goodangle: cosine squared of minangle.                                */
/*     offconstant: constant used to place off-center Steiner points.        */
/*   vararea: -a switch without number.                                      */
/*   fixedarea: -a switch with number.                                       */
/*     maxarea: maximum area bound, specified after -a switch.               */
/*   usertest: -u switch.                                                    */
/*   regionattrib: -A switch.  convex: -c switch.                            */
/*   weighted: 1 for -w switch, 2 for -W switch.  jettison: -j switch        */
/*   firstnumber: inverse of -z switch.  All items are numbered starting     */
/*     from `firstnumber'.                                                   */
/*   edgesout: -e switch.  voronoi: -v switch.                               */
/*   neighbors: -n switch.  geomview: -g switch.                             */
/*   nobound: -B switch.  nopolywritten: -P switch.                          */
/*   nonodewritten: -N switch.  noelewritten: -E switch.                     */
/*   noiterationnum: -I switch.  noholes: -O switch.                         */
/*   noexact: -X switch.                                                     */
/*   order: element order, specified after -o switch.                        */
/*   nobisect: count of how often -Y switch is selected.                     */
/*   steiner: maximum number of Steiner points, specified after -S switch.   */
/*   incremental: -i switch.  sweepline: -F switch.                          */
/*   dwyer: inverse of -l switch.                                            */
/*   splitseg: -s switch.                                                    */
/*   conformdel: -D switch.  docheck: -C switch.                             */
/*   quiet: -Q switch.  verbose: count of how often -V switch is selected.   */
/*   usesegments: -p, -r, -q, or -c switch; determines whether segments are  */
/*     used at all.                                                          */
/*                                                                           */
/* Read the instructions to find out the meaning of these switches.          */

struct behavior 
{
    int poly, refine, quality, vararea, fixedarea, usertest;
    int regionattrib, convex, weighted, jettison;
    int firstnumber;
    int edgesout, voronoi, neighbors, geomview;
    int nobound, nopolywritten, nonodewritten, noelewritten, noiterationnum;
    int noholes, noexact, conformdel;
    int incremental, sweepline, dwyer;
    int splitseg;
    int docheck;
    int quiet, verbose;
    int usesegments;
    int order;
    int nobisect;
    int steiner;
    REAL minangle, goodangle, offconstant;
    REAL maxarea;

};                                              /* End of `struct behavior'. */


/*****************************************************************************/
/*                                                                           */
/*  Mesh manipulation primitives.  Each triangle contains three pointers to  */
/*  other triangles, with orientations.  Each pointer points not to the      */
/*  first byte of a triangle, but to one of the first three bytes of a       */
/*  triangle.  It is necessary to extract both the triangle itself and the   */
/*  orientation.  To save memory, I keep both pieces of information in one   */
/*  pointer.  To make this possible, I assume that all triangles are aligned */
/*  to four-byte boundaries.  The decode() routine below decodes a pointer,  */
/*  extracting an orientation (in the range 0 to 2) and a pointer to the     */
/*  beginning of a triangle.  The encode() routine compresses a pointer to a */
/*  triangle and an orientation into a single pointer.  My assumptions that  */
/*  triangles are four-byte-aligned and that the `unsigned long' type is     */
/*  long enough to hold a pointer are two of the few kludges in this program.*/
/*                                                                           */
/*  Subsegments are manipulated similarly.  A pointer to a subsegment        */
/*  carries both an address and an orientation in the range 0 to 1.          */
/*                                                                           */
/*  The other primitives take an oriented triangle or oriented subsegment,   */
/*  and return an oriented triangle or oriented subsegment or vertex; or     */
/*  they change the connections in the data structure.                       */
/*                                                                           */
/*  Below, triangles and subsegments are denoted by their vertices.  The     */
/*  triangle abc has origin (org) a, destination (dest) b, and apex (apex)   */
/*  c.  These vertices occur in counterclockwise order about the triangle.   */
/*  The handle abc may simultaneously denote vertex a, edge ab, and triangle */
/*  abc.                                                                     */
/*                                                                           */
/*  Similarly, the subsegment ab has origin (sorg) a and destination (sdest) */
/*  b.  If ab is thought to be directed upward (with b directly above a),    */
/*  then the handle ab is thought to grasp the right side of ab, and may     */
/*  simultaneously denote vertex a and edge ab.                              */
/*                                                                           */
/*  An asterisk (*) denotes a vertex whose identity is unknown.              */
/*                                                                           */
/*  Given this notation, a partial list of mesh manipulation primitives      */
/*  follows.                                                                 */
/*                                                                           */
/*                                                                           */
/*  For triangles:                                                           */
/*                                                                           */
/*  sym:  Find the abutting triangle; same edge.                             */
/*  sym(abc) -> ba*                                                          */
/*                                                                           */
/*  lnext:  Find the next edge (counterclockwise) of a triangle.             */
/*  lnext(abc) -> bca                                                        */
/*                                                                           */
/*  lprev:  Find the previous edge (clockwise) of a triangle.                */
/*  lprev(abc) -> cab                                                        */
/*                                                                           */
/*  onext:  Find the next edge counterclockwise with the same origin.        */
/*  onext(abc) -> ac*                                                        */
/*                                                                           */
/*  oprev:  Find the next edge clockwise with the same origin.               */
/*  oprev(abc) -> a*b                                                        */
/*                                                                           */
/*  dnext:  Find the next edge counterclockwise with the same destination.   */
/*  dnext(abc) -> *ba                                                        */
/*                                                                           */
/*  dprev:  Find the next edge clockwise with the same destination.          */
/*  dprev(abc) -> cb*                                                        */
/*                                                                           */
/*  rnext:  Find the next edge (counterclockwise) of the adjacent triangle.  */
/*  rnext(abc) -> *a*                                                        */
/*                                                                           */
/*  rprev:  Find the previous edge (clockwise) of the adjacent triangle.     */
/*  rprev(abc) -> b**                                                        */
/*                                                                           */
/*  org:  Origin          dest:  Destination          apex:  Apex            */
/*  org(abc) -> a         dest(abc) -> b              apex(abc) -> c         */
/*                                                                           */
/*  bond:  Bond two triangles together at the resepective handles.           */
/*  bond(abc, bad)                                                           */
/*                                                                           */
/*                                                                           */
/*  For subsegments:                                                         */
/*                                                                           */
/*  ssym:  Reverse the orientation of a subsegment.                          */
/*  ssym(ab) -> ba                                                           */
/*                                                                           */
/*  spivot:  Find adjoining subsegment with the same origin.                 */
/*  spivot(ab) -> a*                                                         */
/*                                                                           */
/*  snext:  Find next subsegment in sequence.                                */
/*  snext(ab) -> b*                                                          */
/*                                                                           */
/*  sorg:  Origin                      sdest:  Destination                   */
/*  sorg(ab) -> a                      sdest(ab) -> b                        */
/*                                                                           */
/*  sbond:  Bond two subsegments together at the respective origins.         */
/*  sbond(ab, ac)                                                            */
/*                                                                           */
/*                                                                           */
/*  For interacting tetrahedra and subfacets:                                */
/*                                                                           */
/*  tspivot:  Find a subsegment abutting a triangle.                         */
/*  tspivot(abc) -> ba                                                       */
/*                                                                           */
/*  stpivot:  Find a triangle abutting a subsegment.                         */
/*  stpivot(ab) -> ba*                                                       */
/*                                                                           */
/*  tsbond:  Bond a triangle to a subsegment.                                */
/*  tsbond(abc, ba)                                                          */
/*                                                                           */
/*****************************************************************************/

/********* Mesh manipulation primitives begin here                   *********/
/**                                                                         **/
/**                                                                         **/

/* Fast lookup arrays to speed some of the mesh manipulation primitives.     */
extern int plus1mod3[3]; 
extern int minus1mod3[3]; 

/********* Primitives for triangles                                  *********/
/*                                                                           */
/*                                                                           */

/* decode() converts a pointer to an oriented triangle.  The orientation is  */
/*   extracted from the two least significant bits of the pointer.           */

#define decode(ptr, otri)                                                     \
  (otri).orient = (int) ((unsigned long) (ptr) & (unsigned long) 3l);         \
  (otri).tri = (triangle *)                                                   \
                  ((unsigned long) (ptr) ^ (unsigned long) (otri).orient)

/* encode() compresses an oriented triangle into a single pointer.  It       */
/*   relies on the assumption that all triangles are aligned to four-byte    */
/*   boundaries, so the two least significant bits of (otri).tri are zero.   */

#define encode(otri)                                                          \
  (triangle) ((unsigned long) (otri).tri | (unsigned long) (otri).orient)

/* The following handle manipulation primitives are all described by Guibas  */
/*   and Stolfi.  However, Guibas and Stolfi use an edge-based data          */
/*   structure, whereas I use a triangle-based data structure.               */

/* sym() finds the abutting triangle, on the same edge.  Note that the edge  */
/*   direction is necessarily reversed, because the handle specified by an   */
/*   oriented triangle is directed counterclockwise around the triangle.     */

#define sym(otri1, otri2)                                                     \
  ptr = (otri1).tri[(otri1).orient];                                          \
  decode(ptr, otri2);

#define symself(otri)                                                         \
  ptr = (otri).tri[(otri).orient];                                            \
  decode(ptr, otri);

/* lnext() finds the next edge (counterclockwise) of a triangle.             */

#define lnext(otri1, otri2)                                                   \
  (otri2).tri = (otri1).tri;                                                  \
  (otri2).orient = plus1mod3[(otri1).orient]

#define lnextself(otri)                                                       \
  (otri).orient = plus1mod3[(otri).orient]

/* lprev() finds the previous edge (clockwise) of a triangle.                */

#define lprev(otri1, otri2)                                                   \
  (otri2).tri = (otri1).tri;                                                  \
  (otri2).orient = minus1mod3[(otri1).orient]

#define lprevself(otri)                                                       \
  (otri).orient = minus1mod3[(otri).orient]

/* onext() spins counterclockwise around a vertex; that is, it finds the     */
/*   next edge with the same origin in the counterclockwise direction.  This */
/*   edge is part of a different triangle.                                   */

#define onext(otri1, otri2)                                                   \
  lprev(otri1, otri2);                                                        \
  symself(otri2);

#define onextself(otri)                                                       \
  lprevself(otri);                                                            \
  symself(otri);

/* oprev() spins clockwise around a vertex; that is, it finds the next edge  */
/*   with the same origin in the clockwise direction.  This edge is part of  */
/*   a different triangle.                                                   */

#define oprev(otri1, otri2)                                                   \
  sym(otri1, otri2);                                                          \
  lnextself(otri2);

#define oprevself(otri)                                                       \
  symself(otri);                                                              \
  lnextself(otri);

/* dnext() spins counterclockwise around a vertex; that is, it finds the     */
/*   next edge with the same destination in the counterclockwise direction.  */
/*   This edge is part of a different triangle.                              */

#define dnext(otri1, otri2)                                                   \
  sym(otri1, otri2);                                                          \
  lprevself(otri2);

#define dnextself(otri)                                                       \
  symself(otri);                                                              \
  lprevself(otri);

/* dprev() spins clockwise around a vertex; that is, it finds the next edge  */
/*   with the same destination in the clockwise direction.  This edge is     */
/*   part of a different triangle.                                           */

#define dprev(otri1, otri2)                                                   \
  lnext(otri1, otri2);                                                        \
  symself(otri2);

#define dprevself(otri)                                                       \
  lnextself(otri);                                                            \
  symself(otri);

/* rnext() moves one edge counterclockwise about the adjacent triangle.      */
/*   (It's best understood by reading Guibas and Stolfi.  It involves        */
/*   changing triangles twice.)                                              */

#define rnext(otri1, otri2)                                                   \
  sym(otri1, otri2);                                                          \
  lnextself(otri2);                                                           \
  symself(otri2);

#define rnextself(otri)                                                       \
  symself(otri);                                                              \
  lnextself(otri);                                                            \
  symself(otri);

/* rprev() moves one edge clockwise about the adjacent triangle.             */
/*   (It's best understood by reading Guibas and Stolfi.  It involves        */
/*   changing triangles twice.)                                              */

#define rprev(otri1, otri2)                                                   \
  sym(otri1, otri2);                                                          \
  lprevself(otri2);                                                           \
  symself(otri2);

#define rprevself(otri)                                                       \
  symself(otri);                                                              \
  lprevself(otri);                                                            \
  symself(otri);

/* These primitives determine or set the origin, destination, or apex of a   */
/* triangle.                                                                 */

#define org(otri, vertexptr)                                                  \
  vertexptr = (vertex) (otri).tri[plus1mod3[(otri).orient] + 3]

#define dest(otri, vertexptr)                                                 \
  vertexptr = (vertex) (otri).tri[minus1mod3[(otri).orient] + 3]

#define apex(otri, vertexptr)                                                 \
  vertexptr = (vertex) (otri).tri[(otri).orient + 3]

#define setorg(otri, vertexptr)                                               \
  (otri).tri[plus1mod3[(otri).orient] + 3] = (triangle) vertexptr

#define setdest(otri, vertexptr)                                              \
  (otri).tri[minus1mod3[(otri).orient] + 3] = (triangle) vertexptr

#define setapex(otri, vertexptr)                                              \
  (otri).tri[(otri).orient + 3] = (triangle) vertexptr

/* Bond two triangles together.                                              */

#define bond(otri1, otri2)                                                    \
  (otri1).tri[(otri1).orient] = encode(otri2);                                \
  (otri2).tri[(otri2).orient] = encode(otri1)

/* Dissolve a bond (from one side).  Note that the other triangle will still */
/*   think it's connected to this triangle.  Usually, however, the other     */
/*   triangle is being deleted entirely, or bonded to another triangle, so   */
/*   it doesn't matter.                                                      */

#define dissolve(otri)                                                        \
  (otri).tri[(otri).orient] = (triangle) m->dummytri

/* Copy an oriented triangle.                                                */

#define otricopy(otri1, otri2)                                                \
  (otri2).tri = (otri1).tri;                                                  \
  (otri2).orient = (otri1).orient

/* Test for equality of oriented triangles.                                  */

#define otriequal(otri1, otri2)                                               \
  (((otri1).tri == (otri2).tri) &&                                            \
   ((otri1).orient == (otri2).orient))

/* Primitives to infect or cure a triangle with the virus.  These rely on    */
/*   the assumption that all subsegments are aligned to four-byte boundaries.*/

#define infect(otri)                                                          \
  (otri).tri[6] = (triangle)                                                  \
                    ((unsigned long) (otri).tri[6] | (unsigned long) 2l)

#define uninfect(otri)                                                        \
  (otri).tri[6] = (triangle)                                                  \
                    ((unsigned long) (otri).tri[6] & ~ (unsigned long) 2l)

/* Test a triangle for viral infection.                                      */

#define infected(otri)                                                        \
  (((unsigned long) (otri).tri[6] & (unsigned long) 2l) != 0l)

/* Check or set a triangle's attributes.                                     */

#define elemattribute(otri, attnum)                                           \
  ((REAL *) (otri).tri)[m->elemattribindex + (attnum)]

#define setelemattribute(otri, attnum, value)                                 \
  ((REAL *) (otri).tri)[m->elemattribindex + (attnum)] = value

/* Check or set a triangle's maximum area bound.                             */

#define areabound(otri)  ((REAL *) (otri).tri)[m->areaboundindex]

#define setareabound(otri, value)                                             \
  ((REAL *) (otri).tri)[m->areaboundindex] = value

/* Check or set a triangle's deallocation.  Its second pointer is set to     */
/*   NULL to indicate that it is not allocated.  (Its first pointer is used  */
/*   for the stack of dead items.)  Its fourth pointer (its first vertex)    */
/*   is set to NULL in case a `badtriang' structure points to it.            */

#define deadtri(tria)  ((tria)[1] == (triangle) NULL)

#define killtri(tria)                                                         \
  (tria)[1] = (triangle) NULL;                                                \
  (tria)[3] = (triangle) NULL

/********* Primitives for subsegments                                *********/
/*                                                                           */
/*                                                                           */

/* sdecode() converts a pointer to an oriented subsegment.  The orientation  */
/*   is extracted from the least significant bit of the pointer.  The two    */
/*   least significant bits (one for orientation, one for viral infection)   */
/*   are masked out to produce the real pointer.                             */

#define sdecode(sptr, osub)                                                   \
  (osub).ssorient = (int) ((unsigned long) (sptr) & (unsigned long) 1l);      \
  (osub).ss = (subseg *)                                                      \
              ((unsigned long) (sptr) & ~ (unsigned long) 3l)

/* sencode() compresses an oriented subsegment into a single pointer.  It    */
/*   relies on the assumption that all subsegments are aligned to two-byte   */
/*   boundaries, so the least significant bit of (osub).ss is zero.          */

#define sencode(osub)                                                         \
  (subseg) ((unsigned long) (osub).ss | (unsigned long) (osub).ssorient)

/* ssym() toggles the orientation of a subsegment.                           */

#define ssym(osub1, osub2)                                                    \
  (osub2).ss = (osub1).ss;                                                    \
  (osub2).ssorient = 1 - (osub1).ssorient

#define ssymself(osub)                                                        \
  (osub).ssorient = 1 - (osub).ssorient

/* spivot() finds the other subsegment (from the same segment) that shares   */
/*   the same origin.                                                        */

#define spivot(osub1, osub2)                                                  \
  sptr = (osub1).ss[(osub1).ssorient];                                        \
  sdecode(sptr, osub2)

#define spivotself(osub)                                                      \
  sptr = (osub).ss[(osub).ssorient];                                          \
  sdecode(sptr, osub)

/* snext() finds the next subsegment (from the same segment) in sequence;    */
/*   one whose origin is the input subsegment's destination.                 */

#define snext(osub1, osub2)                                                   \
  sptr = (osub1).ss[1 - (osub1).ssorient];                                    \
  sdecode(sptr, osub2)

#define snextself(osub)                                                       \
  sptr = (osub).ss[1 - (osub).ssorient];                                      \
  sdecode(sptr, osub)

/* These primitives determine or set the origin or destination of a          */
/*   subsegment or the segment that includes it.                             */

#define sorg(osub, vertexptr)                                                 \
  vertexptr = (vertex) (osub).ss[2 + (osub).ssorient]

#define sdest(osub, vertexptr)                                                \
  vertexptr = (vertex) (osub).ss[3 - (osub).ssorient]

#define setsorg(osub, vertexptr)                                              \
  (osub).ss[2 + (osub).ssorient] = (subseg) vertexptr

#define setsdest(osub, vertexptr)                                             \
  (osub).ss[3 - (osub).ssorient] = (subseg) vertexptr

#define segorg(osub, vertexptr)                                               \
  vertexptr = (vertex) (osub).ss[4 + (osub).ssorient]

#define segdest(osub, vertexptr)                                              \
  vertexptr = (vertex) (osub).ss[5 - (osub).ssorient]

#define setsegorg(osub, vertexptr)                                            \
  (osub).ss[4 + (osub).ssorient] = (subseg) vertexptr

#define setsegdest(osub, vertexptr)                                           \
  (osub).ss[5 - (osub).ssorient] = (subseg) vertexptr

/* These primitives read or set a boundary marker.  Boundary markers are     */
/*   used to hold user-defined tags for setting boundary conditions in       */
/*   finite element solvers.                                                 */

#define mark(osub)  (* (int *) ((osub).ss + 8))

#define setmark(osub, value)                                                  \
  * (int *) ((osub).ss + 8) = value

/* Bond two subsegments together.                                            */

#define sbond(osub1, osub2)                                                   \
  (osub1).ss[(osub1).ssorient] = sencode(osub2);                              \
  (osub2).ss[(osub2).ssorient] = sencode(osub1)

/* Dissolve a subsegment bond (from one side).  Note that the other          */
/*   subsegment will still think it's connected to this subsegment.          */

#define sdissolve(osub)                                                       \
  (osub).ss[(osub).ssorient] = (subseg) m->dummysub

/* Copy a subsegment.                                                        */

#define subsegcopy(osub1, osub2)                                              \
  (osub2).ss = (osub1).ss;                                                    \
  (osub2).ssorient = (osub1).ssorient

/* Test for equality of subsegments.                                         */

#define subsegequal(osub1, osub2)                                             \
  (((osub1).ss == (osub2).ss) &&                                              \
   ((osub1).ssorient == (osub2).ssorient))

/* Check or set a subsegment's deallocation.  Its second pointer is set to   */
/*   NULL to indicate that it is not allocated.  (Its first pointer is used  */
/*   for the stack of dead items.)  Its third pointer (its first vertex)     */
/*   is set to NULL in case a `badsubseg' structure points to it.            */

#define deadsubseg(sub)  ((sub)[1] == (subseg) NULL)

#define killsubseg(sub)                                                       \
  (sub)[1] = (subseg) NULL;                                                   \
  (sub)[2] = (subseg) NULL

/********* Primitives for interacting triangles and subsegments      *********/
/*                                                                           */
/*                                                                           */

/* tspivot() finds a subsegment abutting a triangle.                         */

#define tspivot(otri, osub)                                                   \
  sptr = (subseg) (otri).tri[6 + (otri).orient];                              \
  sdecode(sptr, osub)

/* stpivot() finds a triangle abutting a subsegment.  It requires that the   */
/*   variable `ptr' of type `triangle' be defined.                           */

#define stpivot(osub, otri)                                                   \
  ptr = (triangle) (osub).ss[6 + (osub).ssorient];                            \
  decode(ptr, otri)

/* Bond a triangle to a subsegment.                                          */

#define tsbond(otri, osub)                                                    \
  (otri).tri[6 + (otri).orient] = (triangle) sencode(osub);                   \
  (osub).ss[6 + (osub).ssorient] = (subseg) encode(otri)

/* Dissolve a bond (from the triangle side).                                 */

#define tsdissolve(otri)                                                      \
  (otri).tri[6 + (otri).orient] = (triangle) m->dummysub

/* Dissolve a bond (from the subsegment side).                               */

#define stdissolve(osub)                                                      \
  (osub).ss[6 + (osub).ssorient] = (subseg) m->dummytri

/********* Primitives for vertices                                   *********/
/*                                                                           */
/*                                                                           */

#define vertexmark(vx)  ((int *) (vx))[m->vertexmarkindex]

#define setvertexmark(vx, value)                                              \
  ((int *) (vx))[m->vertexmarkindex] = value

#define vertextype(vx)  ((int *) (vx))[m->vertexmarkindex + 1]

#define setvertextype(vx, value)                                              \
  ((int *) (vx))[m->vertexmarkindex + 1] = value

#define vertex2tri(vx)  ((triangle *) (vx))[m->vertex2triindex]

#define setvertex2tri(vx, value)                                              \
  ((triangle *) (vx))[m->vertex2triindex] = value

/**                                                                         **/
/**                                                                         **/
/********* Mesh manipulation primitives end here                     *********/

/********* User-defined triangle evaluation routine begins here      *********/
/**                                                                         **/
/**                                                                         **/
/*****************************************************************************/
/*                                                                           */
/*  The `triangulateio' structure.                                           */
/*                                                                           */
/*  Used to pass data into and out of the triangulate() procedure.           */
/*                                                                           */
/*                                                                           */
/*  Arrays are used to store points, triangles, markers, and so forth.  In   */
/*  all cases, the first item in any array is stored starting at index [0].  */
/*  However, that item is item number `1' unless the `z' switch is used, in  */
/*  which case it is item number `0'.  Hence, you may find it easier to      */
/*  index points (and triangles in the neighbor list) if you use the `z'     */
/*  switch.  Unless, of course, you're calling Triangle from a Fortran       */
/*  program.                                                                 */
/*                                                                           */
/*  Description of fields (except the `numberof' fields, which are obvious): */
/*                                                                           */
/*  `pointlist':  An array of point coordinates.  The first point's x        */
/*    coordinate is at index [0] and its y coordinate at index [1], followed */
/*    by the coordinates of the remaining points.  Each point occupies two   */
/*    REALs.                                                                 */
/*  `pointattributelist':  An array of point attributes.  Each point's       */
/*    attributes occupy `numberofpointattributes' REALs.                     */
/*  `pointmarkerlist':  An array of point markers; one int per point.        */
/*                                                                           */
/*  `trianglelist':  An array of triangle corners.  The first triangle's     */
/*    first corner is at index [0], followed by its other two corners in     */
/*    counterclockwise order, followed by any other nodes if the triangle    */
/*    represents a nonlinear element.  Each triangle occupies                */
/*    `numberofcorners' ints.                                                */
/*  `triangleattributelist':  An array of triangle attributes.  Each         */
/*    triangle's attributes occupy `numberoftriangleattributes' REALs.       */
/*  `trianglearealist':  An array of triangle area constraints; one REAL per */
/*    triangle.  Input only.                                                 */
/*  `neighborlist':  An array of triangle neighbors; three ints per          */
/*    triangle.  Output only.                                                */
/*                                                                           */
/*  `segmentlist':  An array of segment endpoints.  The first segment's      */
/*    endpoints are at indices [0] and [1], followed by the remaining        */
/*    segments.  Two ints per segment.                                       */
/*  `segmentmarkerlist':  An array of segment markers; one int per segment.  */
/*                                                                           */
/*  `holelist':  An array of holes.  The first hole's x and y coordinates    */
/*    are at indices [0] and [1], followed by the remaining holes.  Two      */
/*    REALs per hole.  Input only, although the pointer is copied to the     */
/*    output structure for your convenience.                                 */
/*                                                                           */
/*  `regionlist':  An array of regional attributes and area constraints.     */
/*    The first constraint's x and y coordinates are at indices [0] and [1], */
/*    followed by the regional attribute at index [2], followed by the       */
/*    maximum area at index [3], followed by the remaining area constraints. */
/*    Four REALs per area constraint.  Note that each regional attribute is  */
/*    used only if you select the `A' switch, and each area constraint is    */
/*    used only if you select the `a' switch (with no number following), but */
/*    omitting one of these switches does not change the memory layout.      */
/*    Input only, although the pointer is copied to the output structure for */
/*    your convenience.                                                      */
/*                                                                           */
/*  `edgelist':  An array of edge endpoints.  The first edge's endpoints are */
/*    at indices [0] and [1], followed by the remaining edges.  Two ints per */
/*    edge.  Output only.                                                    */
/*  `edgemarkerlist':  An array of edge markers; one int per edge.  Output   */
/*    only.                                                                  */
/*  `normlist':  An array of normal vectors, used for infinite rays in       */
/*    Voronoi diagrams.  The first normal vector's x and y magnitudes are    */
/*    at indices [0] and [1], followed by the remaining vectors.  For each   */
/*    finite edge in a Voronoi diagram, the normal vector written is the     */
/*    zero vector.  Two REALs per edge.  Output only.                        */
/*                                                                           */
/*                                                                           */
/*  Any input fields that Triangle will examine must be initialized.         */
/*  Furthermore, for each output array that Triangle will write to, you      */
/*  must either provide space by setting the appropriate pointer to point    */
/*  to the space you want the data written to, or you must initialize the    */
/*  pointer to NULL, which tells Triangle to allocate space for the results. */
/*  The latter option is preferable, because Triangle always knows exactly   */
/*  how much space to allocate.  The former option is provided mainly for    */
/*  people who need to call Triangle from Fortran code, though it also makes */
/*  possible some nasty space-saving tricks, like writing the output to the  */
/*  same arrays as the input.                                                */
/*                                                                           */
/*  Triangle will not free() any input or output arrays, including those it  */
/*  allocates itself; that's up to you.  You should free arrays allocated by */
/*  Triangle by calling the trifree() procedure defined below.  (By default, */
/*  trifree() just calls the standard free() library procedure, but          */
/*  applications that call triangulate() may replace trimalloc() and         */
/*  trifree() in triangle.c to use specialized memory allocators.)           */
/*                                                                           */
/*  Here's a guide to help you decide which fields you must initialize       */
/*  before you call triangulate().                                           */
/*                                                                           */
/*  `in':                                                                    */
/*                                                                           */
/*    - `pointlist' must always point to a list of points; `numberofpoints'  */
/*      and `numberofpointattributes' must be properly set.                  */
/*      `pointmarkerlist' must either be set to NULL (in which case all      */
/*      markers default to zero), or must point to a list of markers.  If    */
/*      `numberofpointattributes' is not zero, `pointattributelist' must     */
/*      point to a list of point attributes.                                 */
/*    - If the `r' switch is used, `trianglelist' must point to a list of    */
/*      triangles, and `numberoftriangles', `numberofcorners', and           */
/*      `numberoftriangleattributes' must be properly set.  If               */
/*      `numberoftriangleattributes' is not zero, `triangleattributelist'    */
/*      must point to a list of triangle attributes.  If the `a' switch is   */
/*      used (with no number following), `trianglearealist' must point to a  */
/*      list of triangle area constraints.  `neighborlist' may be ignored.   */
/*    - If the `p' switch is used, `segmentlist' must point to a list of     */
/*      segments, `numberofsegments' must be properly set, and               */
/*      `segmentmarkerlist' must either be set to NULL (in which case all    */
/*      markers default to zero), or must point to a list of markers.        */
/*    - If the `p' switch is used without the `r' switch, then               */
/*      `numberofholes' and `numberofregions' must be properly set.  If      */
/*      `numberofholes' is not zero, `holelist' must point to a list of      */
/*      holes.  If `numberofregions' is not zero, `regionlist' must point to */
/*      a list of region constraints.                                        */
/*    - If the `p' switch is used, `holelist', `numberofholes',              */
/*      `regionlist', and `numberofregions' is copied to `out'.  (You can    */
/*      nonetheless get away with not initializing them if the `r' switch is */
/*      used.)                                                               */
/*    - `edgelist', `edgemarkerlist', `normlist', and `numberofedges' may be */
/*      ignored.                                                             */
/*                                                                           */
/*  `out':                                                                   */
/*                                                                           */
/*    - `pointlist' must be initialized (NULL or pointing to memory) unless  */
/*      the `N' switch is used.  `pointmarkerlist' must be initialized       */
/*      unless the `N' or `B' switch is used.  If `N' is not used and        */
/*      `in->numberofpointattributes' is not zero, `pointattributelist' must */
/*      be initialized.                                                      */
/*    - `trianglelist' must be initialized unless the `E' switch is used.    */
/*      `neighborlist' must be initialized if the `n' switch is used.  If    */
/*      the `E' switch is not used and (`in->numberofelementattributes' is   */
/*      not zero or the `A' switch is used), `elementattributelist' must be  */
/*      initialized.  `trianglearealist' may be ignored.                     */
/*    - `segmentlist' must be initialized if the `p' or `c' switch is used,  */
/*      and the `P' switch is not used.  `segmentmarkerlist' must also be    */
/*      initialized under these circumstances unless the `B' switch is used. */
/*    - `edgelist' must be initialized if the `e' switch is used.            */
/*      `edgemarkerlist' must be initialized if the `e' switch is used and   */
/*      the `B' switch is not.                                               */
/*    - `holelist', `regionlist', `normlist', and all scalars may be ignored.*/
/*                                                                           */
/*  `vorout' (only needed if `v' switch is used):                            */
/*                                                                           */
/*    - `pointlist' must be initialized.  If `in->numberofpointattributes'   */
/*      is not zero, `pointattributelist' must be initialized.               */
/*      `pointmarkerlist' may be ignored.                                    */
/*    - `edgelist' and `normlist' must both be initialized.                  */
/*      `edgemarkerlist' may be ignored.                                     */
/*    - Everything else may be ignored.                                      */
/*                                                                           */
/*  After a call to triangulate(), the valid fields of `out' and `vorout'    */
/*  will depend, in an obvious way, on the choice of switches used.  Note    */
/*  that when the `p' switch is used, the pointers `holelist' and            */
/*  `regionlist' are copied from `in' to `out', but no new space is          */
/*  allocated; be careful that you don't free() the same array twice.  On    */
/*  the other hand, Triangle will never copy the `pointlist' pointer (or any */
/*  others); new space is allocated for `out->pointlist', or if the `N'      */
/*  switch is used, `out->pointlist' remains uninitialized.                  */
/*                                                                           */
/*  All of the meaningful `numberof' fields will be properly set; for        */
/*  instance, `numberofedges' will represent the number of edges in the      */
/*  triangulation whether or not the edges were written.  If segments are    */
/*  not used, `numberofsegments' will indicate the number of boundary edges. */
/*                                                                           */
/*****************************************************************************/

struct triangulateio {
  REAL *pointlist;                                               /* In / out */
  REAL *pointattributelist;                                      /* In / out */
  int *pointmarkerlist;                                          /* In / out */
  int numberofpoints;                                            /* In / out */
  int numberofpointattributes;                                   /* In / out */

  int *trianglelist;                                             /* In / out */
  REAL *triangleattributelist;                                   /* In / out */
  REAL *trianglearealist;                                         /* In only */
  int *neighborlist;                                             /* Out only */
  int numberoftriangles;                                         /* In / out */
  int numberofcorners;                                           /* In / out */
  int numberoftriangleattributes;                                /* In / out */

  int *segmentlist;                                              /* In / out */
  int *segmentmarkerlist;                                        /* In / out */
  int numberofsegments;                                          /* In / out */

  REAL *holelist;                        /* In / pointer to array copied out */
  int numberofholes;                                      /* In / copied out */

  REAL *regionlist;                      /* In / pointer to array copied out */
  int numberofregions;                                    /* In / copied out */

  int *edgelist;                                                 /* Out only */
  int *edgemarkerlist;            /* Not used with Voronoi diagram; out only */
  REAL *normlist;                /* Used only with Voronoi diagram; out only */
  int numberofedges;                                             /* Out only */
};

/*============================================================================*/

/* Delaunay triangulation algorithms
 * */ 
typedef enum
{
    FEM2D_TRI_DIVIDE_AND_CONQUER, /* Divide-and-conquer */ 
    FEM2D_TRI_INCREAMENTAL,
    FEM2D_TRI_SWEEP

} FEM2D_TRI_ALG;

long delaunay
( 
    struct mesh *m, 
    struct behavior *b
);

int reconstruct
(
    struct mesh *m, 
    struct behavior *b, 
    int *trianglelist,
    REAL *triangleattriblist, 
    REAL *trianglearealist,
    int elements, 
    int corners, 
    int attribs,
    int *segmentlist,
    int *segmentmarkerlist, 
    int numberofsegments
);
void enforcequality(struct mesh *m, struct behavior *b);

void formskeleton
(
    struct mesh *m, 
    struct behavior *b, 
    int *segmentlist,
    int *segmentmarkerlist, 
    int numberofsegments
);

void triangulate
( 
    char *, 
    struct triangulateio *, 
    struct triangulateio *,
    struct triangulateio *
);

/*============================================================================*/
/********* Geometric primitives begin here                           *********/
/**                                                                         **/
/**                                                                         **/

/* The adaptive exact arithmetic geometric predicates implemented herein are */
/*   described in detail in my paper, "Adaptive Precision Floating-Point     */
/*   Arithmetic and Fast Robust Geometric Predicates."  See the header for a */
/*   full citation.                                                          */

/* Which of the following two methods of finding the absolute values is      */
/*   fastest is compiler-dependent.  A few compilers can inline and optimize */
/*   the fabs() call; but most will incur the overhead of a function call,   */
/*   which is disastrously slow.  A faster way on IEEE machines might be to  */
/*   mask the appropriate bit, but that's difficult to do in C without       */
/*   forcing the value to be stored to memory (rather than be kept in the    */
/*   register to which the optimizer assigned it).                           */

#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a))
/* #define Absolute(a)  fabs(a) */

/* Many of the operations are broken up into two pieces, a main part that    */
/*   performs an approximate operation, and a "tail" that computes the       */
/*   roundoff error of that operation.                                       */
/*                                                                           */
/* The operations Fast_Two_Sum(), Fast_Two_Diff(), Two_Sum(), Two_Diff(),    */
/*   Split(), and Two_Product() are all implemented as described in the      */
/*   reference.  Each of these macros requires certain variables to be       */
/*   defined in the calling routine.  The variables `bvirt', `c', `abig',    */
/*   `_i', `_j', `_k', `_l', `_m', and `_n' are declared `INEXACT' because   */
/*   they store the result of an operation that may incur roundoff error.    */
/*   The input parameter `x' (or the highest numbered `x_' parameter) must   */
/*   also be declared `INEXACT'.                                             */

#define Fast_Two_Sum_Tail(a, b, x, y) \
  bvirt = x - a; \
  y = b - bvirt

#define Fast_Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Fast_Two_Sum_Tail(a, b, x, y)

#define Two_Sum_Tail(a, b, x, y) \
  bvirt = (REAL) (x - a); \
  avirt = x - bvirt; \
  bround = b - bvirt; \
  around = a - avirt; \
  y = around + bround

#define Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Two_Sum_Tail(a, b, x, y)

#define Two_Diff_Tail(a, b, x, y) \
  bvirt = (REAL) (a - x); \
  avirt = x + bvirt; \
  bround = bvirt - b; \
  around = a - avirt; \
  y = around + bround

#define Two_Diff(a, b, x, y) \
  x = (REAL) (a - b); \
  Two_Diff_Tail(a, b, x, y)

#define Split(a, ahi, alo) \
  c = (REAL) (splitter * a); \
  abig = (REAL) (c - a); \
  ahi = c - abig; \
  alo = a - ahi

#define Two_Product_Tail(a, b, x, y) \
  Split(a, ahi, alo); \
  Split(b, bhi, blo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

#define Two_Product(a, b, x, y) \
  x = (REAL) (a * b); \
  Two_Product_Tail(a, b, x, y)

/* Two_Product_Presplit() is Two_Product() where one of the inputs has       */
/*   already been split.  Avoids redundant splitting.                        */

#define Two_Product_Presplit(a, b, bhi, blo, x, y) \
  x = (REAL) (a * b); \
  Split(a, ahi, alo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

/* Square() can be done more quickly than Two_Product().                     */

#define Square_Tail(a, x, y) \
  Split(a, ahi, alo); \
  err1 = x - (ahi * ahi); \
  err3 = err1 - ((ahi + ahi) * alo); \
  y = (alo * alo) - err3

#define Square(a, x, y) \
  x = (REAL) (a * a); \
  Square_Tail(a, x, y)

/* Macros for summing expansions of various fixed lengths.  These are all    */
/*   unrolled versions of Expansion_Sum().                                   */

#define Two_One_Sum(a1, a0, b, x2, x1, x0) \
  Two_Sum(a0, b , _i, x0); \
  Two_Sum(a1, _i, x2, x1)

#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b , _i, x0); \
  Two_Sum( a1, _i, x2, x1)

#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b0, _j, _0, x0); \
  Two_One_Sum(_j, _0, b1, x3, x2, x1)

#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Diff(a1, a0, b0, _j, _0, x0); \
  Two_One_Diff(_j, _0, b1, x3, x2, x1)

/* Macro for multiplying a two-component expansion by a single component.    */

#define Two_One_Product(a1, a0, b, x3, x2, x1, x0) \
  Split(b, bhi, blo); \
  Two_Product_Presplit(a0, b, bhi, blo, _i, x0); \
  Two_Product_Presplit(a1, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x1); \
  Fast_Two_Sum(_j, _k, x3, x2)

void maketriangle(struct mesh *m, struct behavior *b, struct otri *newotri);
void makesubseg(struct mesh *m, struct osub *newsubseg);
void exactinit();
int fast_expansion_sum_zeroelim(int elen, REAL *e, int flen, REAL *f, REAL *h);
int scale_expansion_zeroelim(int elen, REAL *e, REAL b, REAL *h);
REAL estimate(int elen, REAL *e);
REAL counterclockwiseadapt(vertex pa, vertex pb, vertex pc, REAL detsum);
REAL counterclockwise(struct mesh *m, struct behavior *b,
                      vertex pa, vertex pb, vertex pc);
REAL incircleadapt(vertex pa, vertex pb, vertex pc, vertex pd, REAL permanent);

REAL incircle(struct mesh *m, struct behavior *b,
              vertex pa, vertex pb, vertex pc, vertex pd);
REAL orient3dadapt(vertex pa, vertex pb, vertex pc, vertex pd,
                   REAL aheight, REAL bheight, REAL cheight, REAL dheight,
                   REAL permanent);
REAL orient3d(struct mesh *m, struct behavior *b,
              vertex pa, vertex pb, vertex pc, vertex pd,
              REAL aheight, REAL bheight, REAL cheight, REAL dheight);

REAL nonregular(struct mesh *m, struct behavior *b,
                vertex pa, vertex pb, vertex pc, vertex pd);

void findcircumcenter(struct mesh *m, struct behavior *b,
                      vertex torg, vertex tdest, vertex tapex,
                      vertex circumcenter, REAL *xi, REAL *eta, int offcenter);


/*============================================================================*/
/* Data manipulation routines                                                 */
/*============================================================================*/

void init_behavior(struct behavior* b);
void parsecommandline
(
    int argc, 
    char **argv, 
    struct behavior *b
);

void transfernodes
(
    struct mesh     *m, 
    struct behavior *b, 
    REAL* pointlist,
    REAL* pointattriblist, 
    int*  pointmarkerlist,
    int   numberofpoints, 
    int   numberofpointattribs
);

void writenodes
(
    struct mesh *m, 
    struct behavior *b, 
    REAL** pointlist,
    REAL** pointattriblist, 
    int**  pointmarkerlist
);

void numbernodes(struct mesh *m, struct behavior *b);

void writeelements
(
    struct mesh *m, 
    struct behavior *b,
    int**  trianglelist, 
    REAL** triangleattriblist
);

void writepoly
(
    struct mesh *m, 
    struct behavior *b,
    int**  segmentlist, 
    int**  segmentmarkerlist
);

void writeedges
(
    struct mesh *m, 
    struct behavior *b,
    int** edgelist, 
    int** edgemarkerlist
);

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
);

void writeneighbors
(
    struct mesh *m, 
    struct behavior *b, 
    int **neighborlist
);
/*============================================================================*/
/* Memory Management routines */ 

void triexit(int status);
VOID *trimalloc(int size);
void trifree(VOID *memptr);
void poolzero(struct memorypool *pool);
void pool_restart(struct memorypool *pool);
void pool_init
(
    struct memorypool *pool, 
    int bytecount, 
    int itemcount,
    int firstitemcount, 
    int alignment
);
void pooldeinit(struct memorypool *pool);
VOID *poolalloc(struct memorypool *pool);
void pooldealloc(struct memorypool *pool, VOID *dyingitem);
void traversalinit(struct memorypool *pool);
void dummyinit
(
    struct mesh *m, 
    struct behavior *b, 
    int trianglebytes,
    int subsegbytes
);
void initialize_vertex_pool(struct mesh *m, struct behavior *b);
void initializetrisubpools(struct mesh *m, struct behavior *b);
void triangledealloc(struct mesh *m, triangle *dyingtriangle);
void subsegdealloc(struct mesh *m, subseg *dyingsubseg);
subseg *subsegtraverse(struct mesh *m);
void vertexdealloc(struct mesh *m, vertex dyingvertex);

VOID *traverse(struct memorypool *pool);
vertex vertextraverse(struct mesh *m);
vertex getvertex(struct mesh *m, struct behavior *b, int number);
triangle *triangletraverse(struct mesh *m);

#ifndef CDT_ONLY
void badsubsegdealloc(struct mesh *m, struct badsubseg *dyingseg);
struct badsubseg *badsubsegtraverse(struct mesh *m);
#endif /* not CDT_ONLY */

void triangleinit(struct mesh *m);
void triangledeinit(struct mesh *m, struct behavior *b);

/*============================================================================*/
void internalerror(void);
void printtriangle(struct mesh *m, struct behavior *b, struct otri *t);
void printsubseg(struct mesh *m, struct behavior *b, struct osub *s);
/*============================================================================*/
/* Verification routines */ 
void checkmesh(struct mesh *m, struct behavior *b);

int checkdelaunay
(
    struct mesh *m, 
    struct behavior *b
);
void enqueuebadtriang(struct mesh *m, struct behavior *b,
                      struct badtriang *badtri);
void enqueuebadtri(struct mesh *m, struct behavior *b, struct otri *enqtri,
                   REAL minedge, vertex enqapex, vertex enqorg, vertex enqdest);
struct badtriang *dequeuebadtriang(struct mesh *m);
int checkseg4encroach(struct mesh *m, struct behavior *b,
                      struct osub *testsubseg);
void testtriangle(struct mesh *m, struct behavior *b, struct otri *testtri);
int triunsuitable(vertex triorg, vertex tridest, vertex triapex, REAL area);
void quality_statistics(struct mesh *m, struct behavior *b);
void statistics(struct mesh *m, struct behavior *b);
#endif


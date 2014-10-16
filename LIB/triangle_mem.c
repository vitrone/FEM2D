/*============================================================================+/
 *
 * Memory allocation and program exit wrappers begin here
 *
/+============================================================================*/
#include "triangle.h"

/*============================================================================*/





/*============================================================================*/



void triexit(int status)
{
  exit(status);
}

VOID *trimalloc(int size)
{
    VOID *memptr;

    memptr = (VOID *) malloc((unsigned int) size);
    if (memptr == (VOID *) NULL) 
    {
        printf("Error:  Out of memory.\n");
        triexit(1);
    }
    return(memptr);
}

void trifree(VOID *memptr)
{
    if (memptr != NULL)
    {
        free(memptr);
    }
}

/*============================================================================*/
/*****************************************************************************/
/*                                                                           */
/*  poolzero()   Set all of a pool's fields to zero.                         */
/*                                                                           */
/*  This procedure should never be called on a pool that has any memory      */
/*  allocated to it, as that memory would leak.                              */
/*                                                                           */
/*****************************************************************************/

void poolzero(struct memorypool *pool)
{
    pool->firstblock       = (VOID **) NULL;
    pool->nowblock         = (VOID **) NULL;
    pool->nextitem         = (VOID  *) NULL;
    pool->deaditemstack    = (VOID  *) NULL;
    pool->pathblock        = (VOID **) NULL;
    pool->pathitem         = (VOID  *) NULL;
    pool->alignbytes       = 0;
    pool->itembytes        = 0;
    pool->itemsperblock    = 0;
    pool->itemsfirstblock  = 0;
    pool->items            = 0;
    pool->maxitems         = 0;
    pool->unallocateditems = 0;
    pool->pathitemsleft    = 0;
}

/*****************************************************************************/
/*                                                                           */
/*  pool_restart()   Deallocate all items in a pool.                          */
/*                                                                           */
/*  The pool is returned to its starting state, except that no memory is     */
/*  freed to the operating system.  Rather, the previously allocated blocks  */
/*  are ready to be reused.                                                  */
/*                                                                           */
/*****************************************************************************/

void pool_restart(struct memorypool *pool)
{
    unsigned long alignptr;

    pool->items    = 0;
    pool->maxitems = 0;

    /* Set the currently active block. */
    pool->nowblock = pool->firstblock;
    /* Find the first item in the pool.  Increment by the size of (VOID *). 
     * */
    alignptr = (unsigned long) (pool->nowblock + 1);
    /* Align the item on an `alignbytes'-byte boundary. */
    pool->nextitem = (VOID *) (alignptr + (unsigned long) pool->alignbytes -
                     (alignptr % (unsigned long) pool->alignbytes));

    /* There are lots of unallocated items left in this block. */
    pool->unallocateditems = pool->itemsfirstblock;
    /* The stack of deallocated items is empty. */
    pool->deaditemstack = (VOID *) NULL;
}

/*****************************************************************************/
/*                                                                           */
/*  pool_init()   Initialize a pool of memory for allocation of items.        */
/*                                                                           */
/*  This routine initializes the machinery for allocating items.  A `pool'   */
/*  is created whose records have size at least `bytecount'.  Items will be  */
/*  allocated in `itemcount'-item blocks.  Each item is assumed to be a      */
/*  collection of words, and either pointers or floating-point values are    */
/*  assumed to be the "primary" word type.  (The "primary" word type is used */
/*  to determine alignment of items.)  If `alignment' isn't zero, all items  */
/*  will be `alignment'-byte aligned in memory.  `alignment' must be either  */
/*  a multiple or a factor of the primary word size; powers of two are safe. */
/*  `alignment' is normally used to create a few unused bits at the bottom   */
/*  of each item's pointer, in which information may be stored.              */
/*                                                                           */
/*  Don't change this routine unless you understand it.                      */
/*                                                                           */
/*****************************************************************************/

void pool_init
(
    struct memorypool *pool, 
    int bytecount, 
    int itemcount,
    int firstitemcount, 
    int alignment
)
{
    /* Find the proper alignment, which must be at least as large as:
     * - the parameter `alignment'
     * - sizeof(VOID *), so the stack of dead items can be maintained 
     *   without unaligned accesses                                
     * */
    if (alignment > sizeof(VOID *)) 
    {
        pool->alignbytes = alignment;
    } 
    else 
    {
        pool->alignbytes = sizeof(VOID *);
    }
    pool->itembytes = ((bytecount - 1) / pool->alignbytes + 1) * pool->alignbytes;
    
    pool->itemsperblock = itemcount;
    if (firstitemcount == 0) 
    {
        pool->itemsfirstblock = itemcount;
    } 
    else 
    {
        pool->itemsfirstblock = firstitemcount;
    }

    /* Allocate a block of items.  Space for `itemsfirstblock' items and one 
     * pointer (to point to the next block) are allocated, as well as space 
     * to ensure alignment of the items.
     * */
    pool->firstblock = (VOID **) trimalloc( pool->itemsfirstblock * pool->itembytes 
                                            + (int) sizeof(VOID *) 
                                            + pool->alignbytes);
    /* Set the next block pointer to NULL. */
    *(pool->firstblock) = (VOID *) NULL;
    pool_restart(pool);
}

/*****************************************************************************/
/*                                                                           */
/*  pooldeinit()   Free to the operating system all memory taken by a pool.  */
/*                                                                           */
/*****************************************************************************/

void pooldeinit(struct memorypool *pool)
{
    while (pool->firstblock != (VOID **) NULL) 
    {
        pool->nowblock = (VOID **) *(pool->firstblock);
        trifree((VOID *) pool->firstblock);
        pool->firstblock = pool->nowblock;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  poolalloc()   Allocate space for an item.                                */
/*                                                                           */
/*****************************************************************************/

VOID *poolalloc(struct memorypool *pool)
{
    VOID *newitem;
    VOID **newblock;
    unsigned long alignptr;

    /* First check the linked list of dead items.  If the list is not empty, 
     * allocate an item from the list rather than a fresh one. 
     * */
    if (pool->deaditemstack != (VOID *) NULL) 
    {
        /* Take first item in list. */
        newitem = pool->deaditemstack;               
        pool->deaditemstack = * (VOID **) pool->deaditemstack;
    } 
    else 
    {
        /* Check if there are any free items left in the current block. */
        if (pool->unallocateditems == 0) 
        {
            /* Check if another block must be allocated. */
            if (*(pool->nowblock) == (VOID *) NULL) 
            {
                /* Allocate a new block of items, pointed to by the previous block. 
                 * */
                newblock = (VOID **) trimalloc( pool->itemsperblock * pool->itembytes +
                                                (int) sizeof(VOID *) +
                                                pool->alignbytes);
                *(pool->nowblock) = (VOID *) newblock;
                /* The next block pointer is NULL. */
                *newblock = (VOID *) NULL;
            }

            /* Move to the new block. */
            pool->nowblock = (VOID **) *(pool->nowblock);
            /* Find the first item in the block. Increment by the size of (VOID *). 
             * */
            alignptr = (unsigned long) (pool->nowblock + 1);
            /* Align the item on an `alignbytes'-byte boundary. */
            pool->nextitem = (VOID *)(alignptr + (unsigned long) pool->alignbytes -
                             (alignptr % (unsigned long) pool->alignbytes));
            /* There are lots of unallocated items left in this block. 
             * */
            pool->unallocateditems = pool->itemsperblock;
        }

        /* Allocate a new item. */
        newitem = pool->nextitem;
        /* Advance `nextitem' pointer to next free item in block. */
        pool->nextitem = (VOID *) ((char *) pool->nextitem + pool->itembytes);
        pool->unallocateditems--;
        pool->maxitems++;
    }
    pool->items++;
    return newitem;
}

/*****************************************************************************/
/*                                                                           */
/*  pooldealloc()   Deallocate space for an item.                            */
/*                                                                           */
/*  The deallocated space is stored in a queue for later reuse.              */
/*                                                                           */
/*****************************************************************************/

void pooldealloc(struct memorypool *pool, VOID *dyingitem)
{
  /* Push freshly killed item onto stack. */
  *((VOID **) dyingitem) = pool->deaditemstack;
  pool->deaditemstack = dyingitem;
  pool->items--;
}

/*****************************************************************************/
/*                                                                           */
/*  traversalinit()   Prepare to traverse the entire list of items.          */
/*                                                                           */
/*  This routine is used in conjunction with traverse().                     */
/*                                                                           */
/*****************************************************************************/

void traversalinit(struct memorypool *pool)
{
  unsigned long alignptr;

  /* Begin the traversal in the first block. */
  pool->pathblock = pool->firstblock;
  /* Find the first item in the block.  Increment by the size of (VOID *). */
  alignptr = (unsigned long) (pool->pathblock + 1);
  /* Align with item on an `alignbytes'-byte boundary. */
  pool->pathitem = (VOID *)
    (alignptr + (unsigned long) pool->alignbytes -
     (alignptr % (unsigned long) pool->alignbytes));
  /* Set the number of items left in the current block. */
  pool->pathitemsleft = pool->itemsfirstblock;
}

/*****************************************************************************/
/*                                                                           */
/*  traverse()   Find the next item in the list.                             */
/*                                                                           */
/*  This routine is used in conjunction with traversalinit().  Be forewarned */
/*  that this routine successively returns all items in the list, including  */
/*  deallocated ones on the deaditemqueue.  It's up to you to figure out     */
/*  which ones are actually dead.  Why?  I don't want to allocate extra      */
/*  space just to demarcate dead items.  It can usually be done more         */
/*  space-efficiently by a routine that knows something about the structure  */
/*  of the item.                                                             */
/*                                                                           */
/*****************************************************************************/

VOID *traverse(struct memorypool *pool)
{
    VOID *newitem;
    unsigned long alignptr;

    /* Stop upon exhausting the list of items. */
    if (pool->pathitem == pool->nextitem) 
    {
        return (VOID *) NULL;
    }

    /* Check whether any untraversed items remain in the current block. 
     * */
    if (pool->pathitemsleft == 0) 
    {
        /* Find the next block. */
        pool->pathblock = (VOID **) *(pool->pathblock);

        /* Find the first item in the block. 
         * Increment by the size of (VOID *). 
         * */
        alignptr = (unsigned long) (pool->pathblock + 1);
        /* Align with item on an `alignbytes'-byte boundary. 
         * */
        pool->pathitem = (VOID *) ( alignptr + 
                                    (unsigned long) pool->alignbytes - 
                                (alignptr % (unsigned long) pool->alignbytes));
        /* Set the number of items left in the current block. 
         * */
        pool->pathitemsleft = pool->itemsperblock;
    }

    newitem = pool->pathitem;
    /* Find the next item in the block. 
     * */
    pool->pathitem = (VOID *) ((char *) pool->pathitem + pool->itembytes);
    pool->pathitemsleft--;
    return newitem;
}

/*****************************************************************************/
/*                                                                           */
/*  dummyinit()   Initialize the triangle that fills "outer space" and the   */
/*                omnipresent subsegment.                                    */
/*                                                                           */
/*  The triangle that fills "outer space," called `dummytri', is pointed to  */
/*  by every triangle and subsegment on a boundary (be it outer or inner) of */
/*  the triangulation.  Also, `dummytri' points to one of the triangles on   */
/*  the convex hull (until the holes and concavities are carved), making it  */
/*  possible to find a starting triangle for point location.                 */
/*                                                                           */
/*  The omnipresent subsegment, `dummysub', is pointed to by every triangle  */
/*  or subsegment that doesn't have a full complement of real subsegments    */
/*  to point to.                                                             */
/*                                                                           */
/*  `dummytri' and `dummysub' are generally required to fulfill only a few   */
/*  invariants:  their vertices must remain NULL and `dummytri' must always  */
/*  be bonded (at offset zero) to some triangle on the convex hull of the    */
/*  mesh, via a boundary edge.  Otherwise, the connections of `dummytri' and */
/*  `dummysub' may change willy-nilly.  This makes it possible to avoid      */
/*  writing a good deal of special-case code (in the edge flip, for example) */
/*  for dealing with the boundary of the mesh, places where no subsegment is */
/*  present, and so forth.  Other entities are frequently bonded to          */
/*  `dummytri' and `dummysub' as if they were real mesh entities, with no    */
/*  harm done.                                                               */
/*                                                                           */
/*****************************************************************************/

void dummyinit
(
    struct mesh *m, 
    struct behavior *b, 
    int trianglebytes,
    int subsegbytes
)
{
  unsigned long alignptr;

  /* Set up `dummytri', the `triangle' that occupies "outer space." */
  m->dummytribase = (triangle *) trimalloc(trianglebytes +
                                           m->triangles.alignbytes);
  /* Align `dummytri' on a `triangles.alignbytes'-byte boundary. */
  alignptr = (unsigned long) m->dummytribase;
  m->dummytri = (triangle *)
    (alignptr + (unsigned long) m->triangles.alignbytes -
     (alignptr % (unsigned long) m->triangles.alignbytes));
  /* Initialize the three adjoining triangles to be "outer space."  These  */
  /*   will eventually be changed by various bonding operations, but their */
  /*   values don't really matter, as long as they can legally be          */
  /*   dereferenced.                                                       */
  m->dummytri[0] = (triangle) m->dummytri;
  m->dummytri[1] = (triangle) m->dummytri;
  m->dummytri[2] = (triangle) m->dummytri;
  /* Three NULL vertices. */
  m->dummytri[3] = (triangle) NULL;
  m->dummytri[4] = (triangle) NULL;
  m->dummytri[5] = (triangle) NULL;

  if (b->usesegments) {
    /* Set up `dummysub', the omnipresent subsegment pointed to by any */
    /*   triangle side or subsegment end that isn't attached to a real */
    /*   subsegment.                                                   */
    m->dummysubbase = (subseg *) trimalloc(subsegbytes +
                                           m->subsegs.alignbytes);
    /* Align `dummysub' on a `subsegs.alignbytes'-byte boundary. */
    alignptr = (unsigned long) m->dummysubbase;
    m->dummysub = (subseg *)
      (alignptr + (unsigned long) m->subsegs.alignbytes -
       (alignptr % (unsigned long) m->subsegs.alignbytes));
    /* Initialize the two adjoining subsegments to be the omnipresent      */
    /*   subsegment.  These will eventually be changed by various bonding  */
    /*   operations, but their values don't really matter, as long as they */
    /*   can legally be dereferenced.                                      */
    m->dummysub[0] = (subseg) m->dummysub;
    m->dummysub[1] = (subseg) m->dummysub;
    /* Four NULL vertices. */
    m->dummysub[2] = (subseg) NULL;
    m->dummysub[3] = (subseg) NULL;
    m->dummysub[4] = (subseg) NULL;
    m->dummysub[5] = (subseg) NULL;
    /* Initialize the two adjoining triangles to be "outer space." */
    m->dummysub[6] = (subseg) m->dummytri;
    m->dummysub[7] = (subseg) m->dummytri;
    /* Set the boundary marker to zero. */
    * (int *) (m->dummysub + 8) = 0;

    /* Initialize the three adjoining subsegments of `dummytri' to be */
    /*   the omnipresent subsegment.                                  */
    m->dummytri[6] = (triangle) m->dummysub;
    m->dummytri[7] = (triangle) m->dummysub;
    m->dummytri[8] = (triangle) m->dummysub;
  }
}

/*****************************************************************************/
/*                                                                           */
/*  initializevertexpool()   Calculate the size of the vertex data structure */
/*                           and initialize its memory pool.                 */
/*                                                                           */
/*  This routine also computes the `vertexmarkindex' and `vertex2triindex'   */
/*  indices used to find values within each vertex.                          */
/*                                                                           */
/*****************************************************************************/

void initialize_vertex_pool(struct mesh *m, struct behavior *b)
{
    int vertexsize;

    /* The index within each vertex at which the boundary marker is found,
     * followed by the vertex type. Ensure the vertex marker is aligned to 
     * a sizeof(int)-byte address.                                         
     * */
    m->vertexmarkindex = ((m->mesh_dim + m->nextras) * sizeof(REAL) +
                           sizeof(int) - 1) / sizeof(int);
    vertexsize = (m->vertexmarkindex + 2) * sizeof(int);
    if (b->poly) 
    {
        /* The index within each vertex at which a triangle pointer is found. 
         * Ensure the pointer is aligned to a sizeof(triangle)-byte address. 
         * */
        m->vertex2triindex = (vertexsize + sizeof(triangle) - 1) / 
                             sizeof(triangle);
        vertexsize = (m->vertex2triindex + 1) * sizeof(triangle);
    }

    /* Initialize the pool of vertices. */
    pool_init( &m->vertices, 
              vertexsize, 
              VERTEXPERBLOCK,
              m->invertices > VERTEXPERBLOCK ? m->invertices : VERTEXPERBLOCK,
              sizeof(REAL));
}

/*****************************************************************************/
/*                                                                           */
/*  initializetrisubpools()   Calculate the sizes of the triangle and        */
/*                            subsegment data structures and initialize      */
/*                            their memory pools.                            */
/*                                                                           */
/*  This routine also computes the `highorderindex', `elemattribindex', and  */
/*  `areaboundindex' indices used to find values within each triangle.       */
/*                                                                           */
/*****************************************************************************/

void initializetrisubpools(struct mesh *m, struct behavior *b)
{
  int trisize;

  /* The index within each triangle at which the extra nodes (above three)  */
  /*   associated with high order elements are found.  There are three      */
  /*   pointers to other triangles, three pointers to corners, and possibly */
  /*   three pointers to subsegments before the extra nodes.                */
  m->highorderindex = 6 + (b->usesegments * 3);
  /* The number of bytes occupied by a triangle. */
  trisize = ((b->order + 1) * (b->order + 2) / 2 + (m->highorderindex - 3)) *
            sizeof(triangle);
  /* The index within each triangle at which its attributes are found, */
  /*   where the index is measured in REALs.                           */
  m->elemattribindex = (trisize + sizeof(REAL) - 1) / sizeof(REAL);
  /* The index within each triangle at which the maximum area constraint  */
  /*   is found, where the index is measured in REALs.  Note that if the  */
  /*   `regionattrib' flag is set, an additional attribute will be added. */
  m->areaboundindex = m->elemattribindex + m->eextras + b->regionattrib;
  /* If triangle attributes or an area bound are needed, increase the number */
  /*   of bytes occupied by a triangle.                                      */
  if (b->vararea) {
    trisize = (m->areaboundindex + 1) * sizeof(REAL);
  } else if (m->eextras + b->regionattrib > 0) {
    trisize = m->areaboundindex * sizeof(REAL);
  }
  /* If a Voronoi diagram or triangle neighbor graph is requested, make    */
  /*   sure there's room to store an integer index in each triangle.  This */
  /*   integer index can occupy the same space as the subsegment pointers  */
  /*   or attributes or area constraint or extra nodes.                    */
  if ((b->voronoi || b->neighbors) &&
      (trisize < 6 * sizeof(triangle) + sizeof(int))) {
    trisize = 6 * sizeof(triangle) + sizeof(int);
  }

  /* Having determined the memory size of a triangle, initialize the pool. */
  pool_init(&m->triangles, trisize, TRIPERBLOCK,
           (2 * m->invertices - 2) > TRIPERBLOCK ? (2 * m->invertices - 2) :
           TRIPERBLOCK, 4);

  if (b->usesegments) {
    /* Initialize the pool of subsegments.  Take into account all eight */
    /*   pointers and one boundary marker.                              */
    pool_init(&m->subsegs, 8 * sizeof(triangle) + sizeof(int),
             SUBSEGPERBLOCK, SUBSEGPERBLOCK, 4);

    /* Initialize the "outer space" triangle and omnipresent subsegment. */
    dummyinit(m, b, m->triangles.itembytes, m->subsegs.itembytes);
  } else {
    /* Initialize the "outer space" triangle. */
    dummyinit(m, b, m->triangles.itembytes, 0);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  triangledealloc()   Deallocate space for a triangle, marking it dead.    */
/*                                                                           */
/*****************************************************************************/

void triangledealloc(struct mesh *m, triangle *dyingtriangle)
{
  /* Mark the triangle as dead.  This makes it possible to detect dead */
  /*   triangles when traversing the list of all triangles.            */
  killtri(dyingtriangle);
  pooldealloc(&m->triangles, (VOID *) dyingtriangle);
}

/*****************************************************************************/
/*                                                                           */
/*  triangletraverse()   Traverse the triangles, skipping dead ones.         */
/*                                                                           */
/*****************************************************************************/

triangle *triangletraverse(struct mesh *m)
{
    triangle *newtriangle;

    do 
    {
        newtriangle = (triangle *) traverse(&m->triangles);

        if (newtriangle == (triangle *) NULL) 
        {
            return (triangle *) NULL;
        }
    } while (deadtri(newtriangle)); /* Skip dead ones. */

    return newtriangle;
}

/*****************************************************************************/
/*                                                                           */
/*  subsegdealloc()   Deallocate space for a subsegment, marking it dead.    */
/*                                                                           */
/*****************************************************************************/

void subsegdealloc(struct mesh *m, subseg *dyingsubseg)
{
  /* Mark the subsegment as dead.  This makes it possible to detect dead */
  /*   subsegments when traversing the list of all subsegments.          */
  killsubseg(dyingsubseg);
  pooldealloc(&m->subsegs, (VOID *) dyingsubseg);
}

/*****************************************************************************/
/*                                                                           */
/*  subsegtraverse()   Traverse the subsegments, skipping dead ones.         */
/*                                                                           */
/*****************************************************************************/

subseg *subsegtraverse(struct mesh *m)
{
  subseg *newsubseg;

  do {
    newsubseg = (subseg *) traverse(&m->subsegs);
    if (newsubseg == (subseg *) NULL) {
      return (subseg *) NULL;
    }
  } while (deadsubseg(newsubseg));                        /* Skip dead ones. */
  return newsubseg;
}

/*****************************************************************************/
/*                                                                           */
/*  vertexdealloc()   Deallocate space for a vertex, marking it dead.        */
/*                                                                           */
/*****************************************************************************/

void vertexdealloc(struct mesh *m, vertex dyingvertex)
{
  /* Mark the vertex as dead.  This makes it possible to detect dead */
  /*   vertices when traversing the list of all vertices.            */
  setvertextype(dyingvertex, DEADVERTEX);
  pooldealloc(&m->vertices, (VOID *) dyingvertex);
}

/*****************************************************************************/
/*                                                                           */
/*  vertextraverse()   Traverse the vertices, skipping dead ones.            */
/*                                                                           */
/*****************************************************************************/

vertex vertextraverse(struct mesh *m)
{
    vertex newvertex;

    do 
    {
        newvertex = (vertex) traverse(&m->vertices);

        if (newvertex == (vertex) NULL) 
        {
            return (vertex) NULL;
        }

    } while (vertextype(newvertex) == DEADVERTEX);          /* Skip dead ones. */
    return newvertex;
}

/*****************************************************************************/
/*                                                                           */
/*  badsubsegdealloc()   Deallocate space for a bad subsegment, marking it   */
/*                       dead.                                               */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY
void badsubsegdealloc(struct mesh *m, struct badsubseg *dyingseg)
{
  /* Set subsegment's origin to NULL.  This makes it possible to detect dead */
  /*   badsubsegs when traversing the list of all badsubsegs             .   */
  dyingseg->subsegorg = (vertex) NULL;
  pooldealloc(&m->badsubsegs, (VOID *) dyingseg);
}

/*****************************************************************************/
/*                                                                           */
/*  badsubsegtraverse()   Traverse the bad subsegments, skipping dead ones.  */
/*                                                                           */
/*****************************************************************************/

struct badsubseg *badsubsegtraverse(struct mesh *m)
{
  struct badsubseg *newseg;

  do {
    newseg = (struct badsubseg *) traverse(&m->badsubsegs);
    if (newseg == (struct badsubseg *) NULL) {
      return (struct badsubseg *) NULL;
    }
  } while (newseg->subsegorg == (vertex) NULL);           /* Skip dead ones. */
  return newseg;
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  getvertex()   Get a specific vertex, by number, from the list.           */
/*                                                                           */
/*  The first vertex is number 'firstnumber'.                                */
/*                                                                           */
/*  Note that this takes O(n) time (with a small constant, if VERTEXPERBLOCK */
/*  is large).  I don't care to take the trouble to make it work in constant */
/*  time.                                                                    */
/*                                                                           */
/*****************************************************************************/

vertex getvertex
(
    struct mesh *m, 
    struct behavior *b, 
    int number
)
{
    VOID **getblock;
    char *foundvertex;
    unsigned long alignptr;
    int current;

    getblock = m->vertices.firstblock;
    current  = b->firstnumber;

    /* Find the right block. */
    if (current + m->vertices.itemsfirstblock <= number) 
    {
        getblock = (VOID **) *getblock;

        current += m->vertices.itemsfirstblock;
        
        while (current + m->vertices.itemsperblock <= number) 
        {
            getblock = (VOID **) *getblock;
            current += m->vertices.itemsperblock;
        }
    }

    /* Now find the right vertex. */
    alignptr = (unsigned long) (getblock + 1);

    foundvertex = (char *) ( alignptr + 
                            (unsigned long) m->vertices.alignbytes -
                (alignptr % (unsigned long) m->vertices.alignbytes));

    return (vertex) (foundvertex + m->vertices.itembytes * (number - current));
}

/*****************************************************************************/
/*                                                                           */
/*  triangleinit()   Initialize some variables.                              */
/*                                                                           */
/*****************************************************************************/

void triangleinit(struct mesh *m)
{
    poolzero( &m->vertices    );
    poolzero( &m->triangles   );
    poolzero( &m->subsegs     );
    poolzero( &m->viri        );
    poolzero( &m->badsubsegs  );
    poolzero( &m->badtriangles);
    poolzero( &m->flipstackers);
    poolzero( &m->splaynodes  );

    /* No triangle has been visited yet. 
     * */
    m->recenttri.tri = (triangle *) NULL; 
    /* No eliminated input vertices yet. 
     * */
    m->undeads = 0;                       
    /* Point location should take at least one sample. 
     * */
    m->samples = 1;         
    /* There are no segments in the triangulation yet. 
     * */
    m->checksegments = 0;   
    /* The quality triangulation stage has not begun. 
     * */
    m->checkquality      = 0;     
    m->incirclecount     = 0;
    m->counterclockcount = 0;
    m->orient3dcount     = 0;
    m->hyperbolacount    = 0;
    m->circletopcount    = 0;
    m->circumcentercount = 0;
    randomseed = 1;

    /* Initialize exact arithmetic constants. 
     * */
    exactinit();                     
}

/*****************************************************************************/
/*                                                                           */
/*  triangledeinit()   Free all remaining allocated memory.                  */
/*                                                                           */
/*****************************************************************************/

void triangledeinit(struct mesh *m, struct behavior *b)
{
  pooldeinit(&m->triangles);
  trifree((VOID *) m->dummytribase);
  if (b->usesegments) {
    pooldeinit(&m->subsegs);
    trifree((VOID *) m->dummysubbase);
  }
  pooldeinit(&m->vertices);
#ifndef CDT_ONLY
  if (b->quality) {
    pooldeinit(&m->badsubsegs);
    if ((b->minangle > 0.0) || b->vararea || b->fixedarea || b->usertest) {
      pooldeinit(&m->badtriangles);
      pooldeinit(&m->flipstackers);
    }
  }
#endif /* not CDT_ONLY */
}

/* Memory management routines end here                                        */
/* Memory allocation and program exit wrappers end here                       */
/*============================================================================*/


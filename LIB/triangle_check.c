/*============================================================================+/
 * Mesh quality testing routines
/+============================================================================*/
#include "triangle.h"

/*****************************************************************************/
/*                                                                           */
/*  checkmesh()   Test the mesh for topological consistency.                 */
/*                                                                           */
/*****************************************************************************/
void checkmesh
(
    struct mesh *m, 
    struct behavior *b
)
{
  struct otri triangleloop;
  struct otri oppotri, oppooppotri;
  vertex triorg, tridest, triapex;
  vertex oppoorg, oppodest;
  int horrors;
  int saveexact;
  triangle ptr;                         /* Temporary variable used by sym(). */

  /* Temporarily turn on exact arithmetic if it's off. */
  saveexact = b->noexact;
  b->noexact = 0;
  if (!b->quiet) {
    printf("  Checking consistency of mesh...\n");
  }
  horrors = 0;
  /* Run through the list of triangles, checking each one. */
  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  while (triangleloop.tri != (triangle *) NULL) {
    /* Check all three edges of the triangle. */
    for (triangleloop.orient = 0; triangleloop.orient < 3;
         triangleloop.orient++) {
      org(triangleloop, triorg);
      dest(triangleloop, tridest);
      if (triangleloop.orient == 0) {       /* Only test for inversion once. */
        /* Test if the triangle is flat or inverted. */
        apex(triangleloop, triapex);
        if (counterclockwise(m, b, triorg, tridest, triapex) <= 0.0) {
          printf("  !! !! Inverted ");
          printtriangle(m, b, &triangleloop);
          horrors++;
        }
      }
      /* Find the neighboring triangle on this edge. */
      sym(triangleloop, oppotri);
      if (oppotri.tri != m->dummytri) {
        /* Check that the triangle's neighbor knows it's a neighbor. */
        sym(oppotri, oppooppotri);
        if ((triangleloop.tri != oppooppotri.tri)
            || (triangleloop.orient != oppooppotri.orient)) {
          printf("  !! !! Asymmetric triangle-triangle bond:\n");
          if (triangleloop.tri == oppooppotri.tri) {
            printf("   (Right triangle, wrong orientation)\n");
          }
          printf("    First ");
          printtriangle(m, b, &triangleloop);
          printf("    Second (nonreciprocating) ");
          printtriangle(m, b, &oppotri);
          horrors++;
        }
        /* Check that both triangles agree on the identities */
        /*   of their shared vertices.                       */
        org(oppotri, oppoorg);
        dest(oppotri, oppodest);
        if ((triorg != oppodest) || (tridest != oppoorg)) {
          printf("  !! !! Mismatched edge coordinates between two triangles:\n"
                 );
          printf("    First mismatched ");
          printtriangle(m, b, &triangleloop);
          printf("    Second mismatched ");
          printtriangle(m, b, &oppotri);
          horrors++;
        }
      }
    }
    triangleloop.tri = triangletraverse(m);
  }
  if (horrors == 0) {
    if (!b->quiet) {
      printf("  In my studied opinion, the mesh appears to be consistent.\n");
    }
  } else if (horrors == 1) {
    printf("  !! !! !! !! Precisely one festering wound discovered.\n");
  } else {
    printf("  !! !! !! !! %d abominations witnessed.\n", horrors);
  }
  /* Restore the status of exact arithmetic. */
  b->noexact = saveexact;
}

/*****************************************************************************/
/*                                                                           */
/*  checkdelaunay()   Ensure that the mesh is (constrained) Delaunay.        */
/*                                                                           */
/*****************************************************************************/
int checkdelaunay
(
    struct mesh *m, 
    struct behavior *b
)
{
    struct otri triangleloop;
    struct otri oppotri;
    struct osub opposubseg;
    vertex triorg, tridest, triapex;
    vertex oppoapex;
    int shouldbedelaunay;
    int horrors;
    int saveexact;
    triangle ptr;  /* Temporary variable used by sym().     */
    subseg sptr;   /* Temporary variable used by tspivot(). */

    /* Temporarily turn on exact arithmetic if it's off. 
     * */
    saveexact = b->noexact;
    b->noexact = 0;
    if (!b->quiet) 
    {
      printf("Checking Delaunay property of mesh...\n");
    }
    horrors = 0;
    /* Run through the list of triangles, checking each one. 
     * */
    traversalinit(&m->triangles);
    triangleloop.tri = triangletraverse(m);
    while (triangleloop.tri != (triangle *) NULL) 
    {
        /* Check all three edges of the triangle. */
        for ( triangleloop.orient = 0; 
              triangleloop.orient < 3;
              triangleloop.orient++) 
        {
            org(triangleloop, triorg);
            dest(triangleloop, tridest);
            apex(triangleloop, triapex);
            sym(triangleloop, oppotri);
            apex(oppotri, oppoapex);
            /* Only test that the edge is locally Delaunay if there is an 
             * adjoining triangle whose pointer is larger (to ensure that 
             * each pair isn't tested twice).                             
             * */
            shouldbedelaunay = (oppotri.tri != m->dummytri) &&
                  !deadtri(oppotri.tri) && (triangleloop.tri < oppotri.tri) &&
                  (triorg   != m->infvertex1) && (triorg   != m->infvertex2) &&
                  (triorg   != m->infvertex3) &&
                  (tridest  != m->infvertex1) && (tridest  != m->infvertex2) &&
                  (tridest  != m->infvertex3) &&
                  (triapex  != m->infvertex1) && (triapex  != m->infvertex2) &&
                  (triapex  != m->infvertex3) &&
                  (oppoapex != m->infvertex1) && (oppoapex != m->infvertex2) &&
                  (oppoapex != m->infvertex3);

            if (m->checksegments && shouldbedelaunay) 
            {
                /* If a subsegment separates the triangles, then the edge is 
                 * constrained, so no local Delaunay test should be done.  */
                tspivot(triangleloop, opposubseg);
                if (opposubseg.ss != m->dummysub)
                {
                    shouldbedelaunay = 0;
                }
            }
            if (shouldbedelaunay) 
            {
                if (nonregular(m, b, triorg, tridest, triapex, oppoapex) > 0.0)
                {
                    if(!b->quiet)
                    {
                        if (!b->weighted) 
                        {
                            printf("Non-Delaunay pair of triangles:\n");
                            printf("First non-Delaunay ");
                            printtriangle(m, b, &triangleloop);
                            printf("Second non-Delaunay ");
                        } 
                        else
                        {
                            printf("Non-regular pair of triangles:\n");
                            printf("First non-regular ");
                            printtriangle(m, b, &triangleloop);
                            printf("Second non-regular ");
                        }
                        printtriangle(m, b, &oppotri);
                    }
                    horrors++;
                }
            }
        }
        triangleloop.tri = triangletraverse(m);
    }
    if(!b->quiet)
    {
        if (horrors == 0) 
        {
            printf( "By virtue of my perceptive intelligence,"
                    " I declare the mesh Delaunay.\n");
        } 
        else if (horrors == 1)
        {
            printf("Precisely one terrifying transgression identified.\n");
        } 
        else
        {
            printf("%d obscenities viewed with horror.\n", horrors);
        }
    }
    /* Restore the status of exact arithmetic. */
    b->noexact = saveexact;
    return horrors;

}

/*****************************************************************************/
/*                                                                           */
/*  enqueuebadtriang()   Add a bad triangle data structure to the end of a   */
/*                       queue.                                              */
/*                                                                           */
/*  The queue is actually a set of 4096 queues.  I use multiple queues to    */
/*  give priority to smaller angles.  I originally implemented a heap, but   */
/*  the queues are faster by a larger margin than I'd suspected.             */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

void enqueuebadtriang(struct mesh *m, struct behavior *b,
                      struct badtriang *badtri)
{
  REAL length, multiplier;
  int exponent, expincrement;
  int queuenumber;
  int posexponent;
  int i;

  if (b->verbose > 2) {
    printf("  Queueing bad triangle:\n");
    printf("    (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
           badtri->triangorg[0], badtri->triangorg[1],
           badtri->triangdest[0], badtri->triangdest[1],
           badtri->triangapex[0], badtri->triangapex[1]);
  }

  /* Determine the appropriate queue to put the bad triangle into.    */
  /*   Recall that the key is the square of its shortest edge length. */
  if (badtri->key >= 1.0) {
    length = badtri->key;
    posexponent = 1;
  } else {
    /* `badtri->key' is 2.0 to a negative exponent, so we'll record that */
    /*   fact and use the reciprocal of `badtri->key', which is > 1.0.   */
    length = 1.0 / badtri->key;
    posexponent = 0;
  }
  /* `length' is approximately 2.0 to what exponent?  The following code */
  /*   determines the answer in time logarithmic in the exponent.        */
  exponent = 0;
  while (length > 2.0) {
    /* Find an approximation by repeated squaring of two. */
    expincrement = 1;
    multiplier = 0.5;
    while (length * multiplier * multiplier > 1.0) {
      expincrement *= 2;
      multiplier *= multiplier;
    }
    /* Reduce the value of `length', then iterate if necessary. */
    exponent += expincrement;
    length *= multiplier;
  }
  /* `length' is approximately squareroot(2.0) to what exponent? */
  exponent = 2.0 * exponent + (length > SQUAREROOTTWO);
  /* `exponent' is now in the range 0...2047 for IEEE double precision.   */
  /*   Choose a queue in the range 0...4095.  The shortest edges have the */
  /*   highest priority (queue 4095).                                     */
  if (posexponent) {
    queuenumber = 2047 - exponent;
  } else {
    queuenumber = 2048 + exponent;
  }

  /* Are we inserting into an empty queue? */
  if (m->queuefront[queuenumber] == (struct badtriang *) NULL) {
    /* Yes, we are inserting into an empty queue.     */
    /*   Will this become the highest-priority queue? */
    if (queuenumber > m->firstnonemptyq) {
      /* Yes, this is the highest-priority queue. */
      m->nextnonemptyq[queuenumber] = m->firstnonemptyq;
      m->firstnonemptyq = queuenumber;
    } else {
      /* No, this is not the highest-priority queue. */
      /*   Find the queue with next higher priority. */
      i = queuenumber + 1;
      while (m->queuefront[i] == (struct badtriang *) NULL) {
        i++;
      }
      /* Mark the newly nonempty queue as following a higher-priority queue. */
      m->nextnonemptyq[queuenumber] = m->nextnonemptyq[i];
      m->nextnonemptyq[i] = queuenumber;
    }
    /* Put the bad triangle at the beginning of the (empty) queue. */
    m->queuefront[queuenumber] = badtri;
  } else {
    /* Add the bad triangle to the end of an already nonempty queue. */
    m->queuetail[queuenumber]->nexttriang = badtri;
  }
  /* Maintain a pointer to the last triangle of the queue. */
  m->queuetail[queuenumber] = badtri;
  /* Newly enqueued bad triangle has no successor in the queue. */
  badtri->nexttriang = (struct badtriang *) NULL;
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  enqueuebadtri()   Add a bad triangle to the end of a queue.              */
/*                                                                           */
/*  Allocates a badtriang data structure for the triangle, then passes it to */
/*  enqueuebadtriang().                                                      */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

void enqueuebadtri(struct mesh *m, struct behavior *b, struct otri *enqtri,
                   REAL minedge, vertex enqapex, vertex enqorg, vertex enqdest)
{
  struct badtriang *newbad;

  /* Allocate space for the bad triangle. */
  newbad = (struct badtriang *) poolalloc(&m->badtriangles);
  newbad->poortri = encode(*enqtri);
  newbad->key = minedge;
  newbad->triangapex = enqapex;
  newbad->triangorg = enqorg;
  newbad->triangdest = enqdest;
  enqueuebadtriang(m, b, newbad);
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  dequeuebadtriang()   Remove a triangle from the front of the queue.      */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

struct badtriang *dequeuebadtriang(struct mesh *m)
{
  struct badtriang *result;

  /* If no queues are nonempty, return NULL. */
  if (m->firstnonemptyq < 0) {
    return (struct badtriang *) NULL;
  }
  /* Find the first triangle of the highest-priority queue. */
  result = m->queuefront[m->firstnonemptyq];
  /* Remove the triangle from the queue. */
  m->queuefront[m->firstnonemptyq] = result->nexttriang;
  /* If this queue is now empty, note the new highest-priority */
  /*   nonempty queue.                                         */
  if (result == m->queuetail[m->firstnonemptyq]) {
    m->firstnonemptyq = m->nextnonemptyq[m->firstnonemptyq];
  }
  return result;
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  checkseg4encroach()   Check a subsegment to see if it is encroached; add */
/*                        it to the list if it is.                           */
/*                                                                           */
/*  A subsegment is encroached if there is a vertex in its diametral lens.   */
/*  For Ruppert's algorithm (-D switch), the "diametral lens" is the         */
/*  diametral circle.  For Chew's algorithm (default), the diametral lens is */
/*  just big enough to enclose two isosceles triangles whose bases are the   */
/*  subsegment.  Each of the two isosceles triangles has two angles equal    */
/*  to `b->minangle'.                                                        */
/*                                                                           */
/*  Chew's algorithm does not require diametral lenses at all--but they save */
/*  time.  Any vertex inside a subsegment's diametral lens implies that the  */
/*  triangle adjoining the subsegment will be too skinny, so it's only a     */
/*  matter of time before the encroaching vertex is deleted by Chew's        */
/*  algorithm.  It's faster to simply not insert the doomed vertex in the    */
/*  first place, which is why I use diametral lenses with Chew's algorithm.  */
/*                                                                           */
/*  Returns a nonzero value if the subsegment is encroached.                 */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

int checkseg4encroach(struct mesh *m, struct behavior *b,
                      struct osub *testsubseg)
{
  struct otri neighbortri;
  struct osub testsym;
  struct badsubseg *encroachedseg;
  REAL dotproduct;
  int encroached;
  int sides;
  vertex eorg, edest, eapex;
  triangle ptr;                     /* Temporary variable used by stpivot(). */

  encroached = 0;
  sides = 0;

  sorg(*testsubseg, eorg);
  sdest(*testsubseg, edest);
  /* Check one neighbor of the subsegment. */
  stpivot(*testsubseg, neighbortri);
  /* Does the neighbor exist, or is this a boundary edge? */
  if (neighbortri.tri != m->dummytri) {
    sides++;
    /* Find a vertex opposite this subsegment. */
    apex(neighbortri, eapex);
    /* Check whether the apex is in the diametral lens of the subsegment */
    /*   (the diametral circle if `conformdel' is set).  A dot product   */
    /*   of two sides of the triangle is used to check whether the angle */
    /*   at the apex is greater than (180 - 2 `minangle') degrees (for   */
    /*   lenses; 90 degrees for diametral circles).                      */
    dotproduct = (eorg[0] - eapex[0]) * (edest[0] - eapex[0]) +
                 (eorg[1] - eapex[1]) * (edest[1] - eapex[1]);
    if (dotproduct < 0.0) {
      if (b->conformdel ||
          (dotproduct * dotproduct >=
           (2.0 * b->goodangle - 1.0) * (2.0 * b->goodangle - 1.0) *
           ((eorg[0] - eapex[0]) * (eorg[0] - eapex[0]) +
            (eorg[1] - eapex[1]) * (eorg[1] - eapex[1])) *
           ((edest[0] - eapex[0]) * (edest[0] - eapex[0]) +
            (edest[1] - eapex[1]) * (edest[1] - eapex[1])))) {
        encroached = 1;
      }
    }
  }
  /* Check the other neighbor of the subsegment. */
  ssym(*testsubseg, testsym);
  stpivot(testsym, neighbortri);
  /* Does the neighbor exist, or is this a boundary edge? */
  if (neighbortri.tri != m->dummytri) {
    sides++;
    /* Find the other vertex opposite this subsegment. */
    apex(neighbortri, eapex);
    /* Check whether the apex is in the diametral lens of the subsegment */
    /*   (or the diametral circle, if `conformdel' is set).              */
    dotproduct = (eorg[0] - eapex[0]) * (edest[0] - eapex[0]) +
                 (eorg[1] - eapex[1]) * (edest[1] - eapex[1]);
    if (dotproduct < 0.0) {
      if (b->conformdel ||
          (dotproduct * dotproduct >=
           (2.0 * b->goodangle - 1.0) * (2.0 * b->goodangle - 1.0) *
           ((eorg[0] - eapex[0]) * (eorg[0] - eapex[0]) +
            (eorg[1] - eapex[1]) * (eorg[1] - eapex[1])) *
           ((edest[0] - eapex[0]) * (edest[0] - eapex[0]) +
            (edest[1] - eapex[1]) * (edest[1] - eapex[1])))) {
        encroached += 2;
      }
    }
  }

  if (encroached && (!b->nobisect || ((b->nobisect == 1) && (sides == 2)))) {
    if (b->verbose > 2) {
      printf(
        "  Queueing encroached subsegment (%.12g, %.12g) (%.12g, %.12g).\n",
        eorg[0], eorg[1], edest[0], edest[1]);
    }
    /* Add the subsegment to the list of encroached subsegments. */
    /*   Be sure to get the orientation right.                   */
    encroachedseg = (struct badsubseg *) poolalloc(&m->badsubsegs);
    if (encroached == 1) {
      encroachedseg->encsubseg = sencode(*testsubseg);
      encroachedseg->subsegorg = eorg;
      encroachedseg->subsegdest = edest;
    } else {
      encroachedseg->encsubseg = sencode(testsym);
      encroachedseg->subsegorg = edest;
      encroachedseg->subsegdest = eorg;
    }
  }

  return encroached;
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  testtriangle()   Test a triangle for quality and size.                   */
/*                                                                           */
/*  Tests a triangle to see if it satisfies the minimum angle condition and  */
/*  the maximum area condition.  Triangles that aren't up to spec are added  */
/*  to the bad triangle queue.                                               */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

void testtriangle(struct mesh *m, struct behavior *b, struct otri *testtri)
{
  struct otri tri1, tri2;
  struct osub testsub;
  vertex torg, tdest, tapex;
  vertex base1, base2;
  vertex org1, dest1, org2, dest2;
  vertex joinvertex;
  REAL dxod, dyod, dxda, dyda, dxao, dyao;
  REAL dxod2, dyod2, dxda2, dyda2, dxao2, dyao2;
  REAL apexlen, orglen, destlen, minedge;
  REAL angle;
  REAL area;
  REAL dist1, dist2;
  subseg sptr;                      /* Temporary variable used by tspivot(). */
  triangle ptr;           /* Temporary variable used by oprev() and dnext(). */

  org(*testtri, torg);
  dest(*testtri, tdest);
  apex(*testtri, tapex);
  dxod = torg[0] - tdest[0];
  dyod = torg[1] - tdest[1];
  dxda = tdest[0] - tapex[0];
  dyda = tdest[1] - tapex[1];
  dxao = tapex[0] - torg[0];
  dyao = tapex[1] - torg[1];
  dxod2 = dxod * dxod;
  dyod2 = dyod * dyod;
  dxda2 = dxda * dxda;
  dyda2 = dyda * dyda;
  dxao2 = dxao * dxao;
  dyao2 = dyao * dyao;
  /* Find the lengths of the triangle's three edges. */
  apexlen = dxod2 + dyod2;
  orglen = dxda2 + dyda2;
  destlen = dxao2 + dyao2;

  if ((apexlen < orglen) && (apexlen < destlen)) {
    /* The edge opposite the apex is shortest. */
    minedge = apexlen;
    /* Find the square of the cosine of the angle at the apex. */
    angle = dxda * dxao + dyda * dyao;
    angle = angle * angle / (orglen * destlen);
    base1 = torg;
    base2 = tdest;
    otricopy(*testtri, tri1);
  } else if (orglen < destlen) {
    /* The edge opposite the origin is shortest. */
    minedge = orglen;
    /* Find the square of the cosine of the angle at the origin. */
    angle = dxod * dxao + dyod * dyao;
    angle = angle * angle / (apexlen * destlen);
    base1 = tdest;
    base2 = tapex;
    lnext(*testtri, tri1);
  } else {
    /* The edge opposite the destination is shortest. */
    minedge = destlen;
    /* Find the square of the cosine of the angle at the destination. */
    angle = dxod * dxda + dyod * dyda;
    angle = angle * angle / (apexlen * orglen);
    base1 = tapex;
    base2 = torg;
    lprev(*testtri, tri1);
  }

  if (b->vararea || b->fixedarea || b->usertest) {
    /* Check whether the area is larger than permitted. */
    area = 0.5 * (dxod * dyda - dyod * dxda);
    if (b->fixedarea && (area > b->maxarea)) {
      /* Add this triangle to the list of bad triangles. */
      enqueuebadtri(m, b, testtri, minedge, tapex, torg, tdest);
      return;
    }

    /* Nonpositive area constraints are treated as unconstrained. */
    if ((b->vararea) && (area > areabound(*testtri)) &&
        (areabound(*testtri) > 0.0)) {
      /* Add this triangle to the list of bad triangles. */
      enqueuebadtri(m, b, testtri, minedge, tapex, torg, tdest);
      return;
    }

    if (b->usertest) {
      /* Check whether the user thinks this triangle is too large. */
      if (triunsuitable(torg, tdest, tapex, area)) {
        enqueuebadtri(m, b, testtri, minedge, tapex, torg, tdest);
        return;
      }
    }
  }

  /* Check whether the angle is smaller than permitted. */
  if (angle > b->goodangle) {
    /* Use the rules of Miller, Pav, and Walkington to decide that certain */
    /*   triangles should not be split, even if they have bad angles.      */
    /*   A skinny triangle is not split if its shortest edge subtends a    */
    /*   small input angle, and both endpoints of the edge lie on a        */
    /*   concentric circular shell.  For convenience, I make a small       */
    /*   adjustment to that rule:  I check if the endpoints of the edge    */
    /*   both lie in segment interiors, equidistant from the apex where    */
    /*   the two segments meet.                                            */
    /* First, check if both points lie in segment interiors.               */
    if ((vertextype(base1) == SEGMENTVERTEX) &&
        (vertextype(base2) == SEGMENTVERTEX)) {
      /* Check if both points lie in a common segment.  If they do, the */
      /*   skinny triangle is enqueued to be split as usual.            */
      tspivot(tri1, testsub);
      if (testsub.ss == m->dummysub) {
        /* No common segment.  Find a subsegment that contains `torg'. */
        otricopy(tri1, tri2);
        do {
          oprevself(tri1);
          tspivot(tri1, testsub);
        } while (testsub.ss == m->dummysub);
        /* Find the endpoints of the containing segment. */
        segorg(testsub, org1);
        segdest(testsub, dest1);
        /* Find a subsegment that contains `tdest'. */
        do {
          dnextself(tri2);
          tspivot(tri2, testsub);
        } while (testsub.ss == m->dummysub);
        /* Find the endpoints of the containing segment. */
        segorg(testsub, org2);
        segdest(testsub, dest2);
        /* Check if the two containing segments have an endpoint in common. */
        joinvertex = (vertex) NULL;
        if ((dest1[0] == org2[0]) && (dest1[1] == org2[1])) {
          joinvertex = dest1;
        } else if ((org1[0] == dest2[0]) && (org1[1] == dest2[1])) {
          joinvertex = org1;
        }
        if (joinvertex != (vertex) NULL) {
          /* Compute the distance from the common endpoint (of the two  */
          /*   segments) to each of the endpoints of the shortest edge. */
          dist1 = ((base1[0] - joinvertex[0]) * (base1[0] - joinvertex[0]) +
                   (base1[1] - joinvertex[1]) * (base1[1] - joinvertex[1]));
          dist2 = ((base2[0] - joinvertex[0]) * (base2[0] - joinvertex[0]) +
                   (base2[1] - joinvertex[1]) * (base2[1] - joinvertex[1]));
          /* If the two distances are equal, don't split the triangle. */
          if ((dist1 < 1.001 * dist2) && (dist1 > 0.999 * dist2)) {
            /* Return now to avoid enqueueing the bad triangle. */
            return;
          }
        }
      }
    }

    /* Add this triangle to the list of bad triangles. */
    enqueuebadtri(m, b, testtri, minedge, tapex, torg, tdest);
  }
}

#endif /* not CDT_ONLY */
/*****************************************************************************/
/*                                                                           */
/*  triunsuitable()   Determine if a triangle is unsuitable, and thus must   */
/*                    be further refined.                                    */
/*                                                                           */
/*  You may write your own procedure that decides whether or not a selected  */
/*  triangle is too big (and needs to be refined).  There are two ways to do */
/*  this.                                                                    */
/*                                                                           */
/*  (1)  Modify the procedure `triunsuitable' below, then recompile          */
/*  Triangle.                                                                */
/*                                                                           */
/*  (2)  Define the symbol EXTERNAL_TEST (either by adding the definition    */
/*  to this file, or by using the appropriate compiler switch).  This way,   */
/*  you can compile triangle.c separately from your test.  Write your own    */
/*  `triunsuitable' procedure in a separate C file (using the same prototype */
/*  as below).  Compile it and link the object code with triangle.o.         */
/*                                                                           */
/*  This procedure returns 1 if the triangle is too large and should be      */
/*  refined; 0 otherwise.                                                    */
/*                                                                           */
/*****************************************************************************/
int triunsuitable(vertex triorg, vertex tridest, vertex triapex, REAL area)
/* The triangle's origin vertex. */
/* The triangle's destination vertex. */
/* The triangle's apex vertex. */
/* The area of the triangle. */
{
  REAL dxoa, dxda, dxod;
  REAL dyoa, dyda, dyod;
  REAL oalen, dalen, odlen;
  REAL maxlen;

  dxoa = triorg[0] - triapex[0];
  dyoa = triorg[1] - triapex[1];
  dxda = tridest[0] - triapex[0];
  dyda = tridest[1] - triapex[1];
  dxod = triorg[0] - tridest[0];
  dyod = triorg[1] - tridest[1];
  /* Find the squares of the lengths of the triangle's three edges. */
  oalen = dxoa * dxoa + dyoa * dyoa;
  dalen = dxda * dxda + dyda * dyda;
  odlen = dxod * dxod + dyod * dyod;
  /* Find the square of the length of the longest edge. */
  maxlen = (dalen > oalen) ? dalen : oalen;
  maxlen = (odlen > maxlen) ? odlen : maxlen;

  if (maxlen > 0.05 * (triorg[0] * triorg[0] + triorg[1] * triorg[1]) + 0.02) {
    return 1;
  } else {
    return 0;
  }
}


/*****************************************************************************/
/*                                                                           */
/*  quality_statistics()   Print statistics about the quality of the mesh.   */
/*                                                                           */
/*****************************************************************************/

void quality_statistics(struct mesh *m, struct behavior *b)
{
  struct otri triangleloop;
  vertex p[3];
  REAL cossquaretable[8];
  REAL ratiotable[16];
  REAL dx[3], dy[3];
  REAL edgelength[3];
  REAL dotproduct;
  REAL cossquare;
  REAL triarea;
  REAL shortest, longest;
  REAL trilongest2;
  REAL smallestarea, biggestarea;
  REAL triminaltitude2;
  REAL minaltitude;
  REAL triaspect2;
  REAL worstaspect;
  REAL smallestangle, biggestangle;
  REAL radconst, degconst;
  int angletable[18];
  int aspecttable[16];
  int aspectindex;
  int tendegree;
  int acutebiggest;
  int i, ii, j, k;

  printf("Mesh quality statistics:\n\n");
  radconst = PI / 18.0;
  degconst = 180.0 / PI;
  for (i = 0; i < 8; i++) {
    cossquaretable[i] = cos(radconst * (REAL) (i + 1));
    cossquaretable[i] = cossquaretable[i] * cossquaretable[i];
  }
  for (i = 0; i < 18; i++) {
    angletable[i] = 0;
  }

  ratiotable[0]  =      1.5;      ratiotable[1]  =     2.0;
  ratiotable[2]  =      2.5;      ratiotable[3]  =     3.0;
  ratiotable[4]  =      4.0;      ratiotable[5]  =     6.0;
  ratiotable[6]  =     10.0;      ratiotable[7]  =    15.0;
  ratiotable[8]  =     25.0;      ratiotable[9]  =    50.0;
  ratiotable[10] =    100.0;      ratiotable[11] =   300.0;
  ratiotable[12] =   1000.0;      ratiotable[13] = 10000.0;
  ratiotable[14] = 100000.0;      ratiotable[15] =     0.0;
  for (i = 0; i < 16; i++) {
    aspecttable[i] = 0;
  }

  worstaspect = 0.0;
  minaltitude = m->xmax - m->xmin + m->ymax - m->ymin;
  minaltitude = minaltitude * minaltitude;
  shortest = minaltitude;
  longest = 0.0;
  smallestarea = minaltitude;
  biggestarea = 0.0;
  worstaspect = 0.0;
  smallestangle = 0.0;
  biggestangle = 2.0;
  acutebiggest = 1;

  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  triangleloop.orient = 0;
  while (triangleloop.tri != (triangle *) NULL) {
    org(triangleloop, p[0]);
    dest(triangleloop, p[1]);
    apex(triangleloop, p[2]);
    trilongest2 = 0.0;

    for (i = 0; i < 3; i++) {
      j = plus1mod3[i];
      k = minus1mod3[i];
      dx[i] = p[j][0] - p[k][0];
      dy[i] = p[j][1] - p[k][1];
      edgelength[i] = dx[i] * dx[i] + dy[i] * dy[i];
      if (edgelength[i] > trilongest2) {
        trilongest2 = edgelength[i];
      }
      if (edgelength[i] > longest) {
        longest = edgelength[i];
      }
      if (edgelength[i] < shortest) {
        shortest = edgelength[i];
      }
    }

    triarea = counterclockwise(m, b, p[0], p[1], p[2]);
    if (triarea < smallestarea) {
      smallestarea = triarea;
    }
    if (triarea > biggestarea) {
      biggestarea = triarea;
    }
    triminaltitude2 = triarea * triarea / trilongest2;
    if (triminaltitude2 < minaltitude) {
      minaltitude = triminaltitude2;
    }
    triaspect2 = trilongest2 / triminaltitude2;
    if (triaspect2 > worstaspect) {
      worstaspect = triaspect2;
    }
    aspectindex = 0;
    while ((triaspect2 > ratiotable[aspectindex] * ratiotable[aspectindex])
           && (aspectindex < 15)) {
      aspectindex++;
    }
    aspecttable[aspectindex]++;

    for (i = 0; i < 3; i++) {
      j = plus1mod3[i];
      k = minus1mod3[i];
      dotproduct = dx[j] * dx[k] + dy[j] * dy[k];
      cossquare = dotproduct * dotproduct / (edgelength[j] * edgelength[k]);
      tendegree = 8;
      for (ii = 7; ii >= 0; ii--) {
        if (cossquare > cossquaretable[ii]) {
          tendegree = ii;
        }
      }
      if (dotproduct <= 0.0) {
        angletable[tendegree]++;
        if (cossquare > smallestangle) {
          smallestangle = cossquare;
        }
        if (acutebiggest && (cossquare < biggestangle)) {
          biggestangle = cossquare;
        }
      } else {
        angletable[17 - tendegree]++;
        if (acutebiggest || (cossquare > biggestangle)) {
          biggestangle = cossquare;
          acutebiggest = 0;
        }
      }
    }
    triangleloop.tri = triangletraverse(m);
  }

  shortest = sqrt(shortest);
  longest = sqrt(longest);
  minaltitude = sqrt(minaltitude);
  worstaspect = sqrt(worstaspect);
  smallestarea *= 0.5;
  biggestarea *= 0.5;
  if (smallestangle >= 1.0) {
    smallestangle = 0.0;
  } else {
    smallestangle = degconst * acos(sqrt(smallestangle));
  }
  if (biggestangle >= 1.0) {
    biggestangle = 180.0;
  } else {
    if (acutebiggest) {
      biggestangle = degconst * acos(sqrt(biggestangle));
    } else {
      biggestangle = 180.0 - degconst * acos(sqrt(biggestangle));
    }
  }

  printf("  Smallest area: %16.5g   |  Largest area: %16.5g\n",
         smallestarea, biggestarea);
  printf("  Shortest edge: %16.5g   |  Longest edge: %16.5g\n",
         shortest, longest);
  printf("  Shortest altitude: %12.5g   |  Largest aspect ratio: %8.5g\n\n",
         minaltitude, worstaspect);

  printf("  Triangle aspect ratio histogram:\n");
  printf("  1.1547 - %-6.6g    :  %8d    | %6.6g - %-6.6g     :  %8d\n",
         ratiotable[0], aspecttable[0], ratiotable[7], ratiotable[8],
         aspecttable[8]);
  for (i = 1; i < 7; i++) {
    printf("  %6.6g - %-6.6g    :  %8d    | %6.6g - %-6.6g     :  %8d\n",
           ratiotable[i - 1], ratiotable[i], aspecttable[i],
           ratiotable[i + 7], ratiotable[i + 8], aspecttable[i + 8]);
  }
  printf("  %6.6g - %-6.6g    :  %8d    | %6.6g -            :  %8d\n",
         ratiotable[6], ratiotable[7], aspecttable[7], ratiotable[14],
         aspecttable[15]);
  printf("  (Aspect ratio is longest edge divided by shortest altitude)\n\n");

  printf("  Smallest angle: %15.5g   |  Largest angle: %15.5g\n\n",
         smallestangle, biggestangle);

  printf("  Angle histogram:\n");
  for (i = 0; i < 9; i++) {
    printf("    %3d - %3d degrees:  %8d    |    %3d - %3d degrees:  %8d\n",
           i * 10, i * 10 + 10, angletable[i],
           i * 10 + 90, i * 10 + 100, angletable[i + 9]);
  }
  printf("\n");
}

/*****************************************************************************/
/*                                                                           */
/*  statistics()   Print all sorts of cool facts.                            */
/*                                                                           */
/*****************************************************************************/

void statistics(struct mesh *m, struct behavior *b)
{
  printf("\nStatistics:\n\n");
  printf("  Input vertices: %d\n", m->invertices);
  if (b->refine) {
    printf("  Input triangles: %d\n", m->inelements);
  }
  if (b->poly) {
    printf("  Input segments: %d\n", m->insegments);
    if (!b->refine) {
      printf("  Input holes: %d\n", m->holes);
    }
  }

  printf("\n  Mesh vertices: %ld\n", m->vertices.items - m->undeads);
  printf("  Mesh triangles: %ld\n", m->triangles.items);
  printf("  Mesh edges: %ld\n", m->edges);
  printf("  Mesh exterior boundary edges: %ld\n", m->hullsize);
  if (b->poly || b->refine) {
    printf("  Mesh interior boundary edges: %ld\n",
           m->subsegs.items - m->hullsize);
    printf("  Mesh subsegments (constrained edges): %ld\n",
           m->subsegs.items);
  }
  printf("\n");

  if (b->verbose) {
    quality_statistics(m, b);
    printf("Memory allocation statistics:\n\n");
    printf("  Maximum number of vertices: %ld\n", m->vertices.maxitems);
    printf("  Maximum number of triangles: %ld\n", m->triangles.maxitems);
    if (m->subsegs.maxitems > 0) {
      printf("  Maximum number of subsegments: %ld\n", m->subsegs.maxitems);
    }
    if (m->viri.maxitems > 0) {
      printf("  Maximum number of viri: %ld\n", m->viri.maxitems);
    }
    if (m->badsubsegs.maxitems > 0) {
      printf("  Maximum number of encroached subsegments: %ld\n",
             m->badsubsegs.maxitems);
    }
    if (m->badtriangles.maxitems > 0) {
      printf("  Maximum number of bad triangles: %ld\n",
             m->badtriangles.maxitems);
    }
    if (m->flipstackers.maxitems > 0) {
      printf("  Maximum number of stacked triangle flips: %ld\n",
             m->flipstackers.maxitems);
    }
    if (m->splaynodes.maxitems > 0) {
      printf("  Maximum number of splay tree nodes: %ld\n",
             m->splaynodes.maxitems);
    }
    printf("  Approximate heap memory use (bytes): %ld\n\n",
           m->vertices.maxitems * m->vertices.itembytes +
           m->triangles.maxitems * m->triangles.itembytes +
           m->subsegs.maxitems * m->subsegs.itembytes +
           m->viri.maxitems * m->viri.itembytes +
           m->badsubsegs.maxitems * m->badsubsegs.itembytes +
           m->badtriangles.maxitems * m->badtriangles.itembytes +
           m->flipstackers.maxitems * m->flipstackers.itembytes +
           m->splaynodes.maxitems * m->splaynodes.itembytes);

    printf("Algorithmic statistics:\n\n");
    if (!b->weighted) {
      printf("  Number of incircle tests: %ld\n", m->incirclecount);
    } else {
      printf("  Number of 3D orientation tests: %ld\n", m->orient3dcount);
    }
    printf("  Number of 2D orientation tests: %ld\n", m->counterclockcount);
    if (m->hyperbolacount > 0) {
      printf("  Number of right-of-hyperbola tests: %ld\n",
             m->hyperbolacount);
    }
    if (m->circletopcount > 0) {
      printf("  Number of circle top computations: %ld\n",
             m->circletopcount);
    }
    if (m->circumcentercount > 0) {
      printf("  Number of triangle circumcenter computations: %ld\n",
             m->circumcentercount);
    }
    printf("\n");
  }
}

/**                                                                         **/
/**                                                                         **/
/********* Mesh quality testing routines end here                    *********/

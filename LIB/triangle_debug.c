/*============================================================================+/
 * Debugging routines
/+============================================================================*/
#include "triangle.h"

/*****************************************************************************/
/*                                                                           */
/*  internalerror()   Ask the user to send me the defective product.  Exit.  */
/*                                                                           */
/*****************************************************************************/

void internalerror(void)
{
    printf(" Please report this bug to jrs@cs.berkeley.edu\n");
    printf(" Include the message above, your input data set, and the exact\n");
    printf(" command line you used to run Triangle.\n");
    triexit(1);
}

/*****************************************************************************/
/*                                                                           */
/*  printtriangle()   Print out the details of an oriented triangle.         */
/*                                                                           */
/*  I originally wrote this procedure to simplify debugging; it can be       */
/*  called directly from the debugger, and presents information about an     */
/*  oriented triangle in digestible form.  It's also used when the           */
/*  highest level of verbosity (`-VVV') is specified.                        */
/*                                                                           */
/*****************************************************************************/

void printtriangle(struct mesh *m, struct behavior *b, struct otri *t)
{
  struct otri printtri;
  struct osub printsh;
  vertex printvertex;

  printf("triangle x%lx with orientation %d:\n", (unsigned long) t->tri,
         t->orient);
  decode(t->tri[0], printtri);
  if (printtri.tri == m->dummytri) {
    printf("    [0] = Outer space\n");
  } else {
    printf("    [0] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }
  decode(t->tri[1], printtri);
  if (printtri.tri == m->dummytri) {
    printf("    [1] = Outer space\n");
  } else {
    printf("    [1] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }
  decode(t->tri[2], printtri);
  if (printtri.tri == m->dummytri) {
    printf("    [2] = Outer space\n");
  } else {
    printf("    [2] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }

  org(*t, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Origin[%d] = NULL\n", (t->orient + 1) % 3 + 3);
  else
    printf("    Origin[%d] = x%lx  (%.12g, %.12g)\n",
           (t->orient + 1) % 3 + 3, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);
  dest(*t, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Dest  [%d] = NULL\n", (t->orient + 2) % 3 + 3);
  else
    printf("    Dest  [%d] = x%lx  (%.12g, %.12g)\n",
           (t->orient + 2) % 3 + 3, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);
  apex(*t, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Apex  [%d] = NULL\n", t->orient + 3);
  else
    printf("    Apex  [%d] = x%lx  (%.12g, %.12g)\n",
           t->orient + 3, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);

  if (b->usesegments) {
    sdecode(t->tri[6], printsh);
    if (printsh.ss != m->dummysub) {
      printf("    [6] = x%lx  %d\n", (unsigned long) printsh.ss,
             printsh.ssorient);
    }
    sdecode(t->tri[7], printsh);
    if (printsh.ss != m->dummysub) {
      printf("    [7] = x%lx  %d\n", (unsigned long) printsh.ss,
             printsh.ssorient);
    }
    sdecode(t->tri[8], printsh);
    if (printsh.ss != m->dummysub) {
      printf("    [8] = x%lx  %d\n", (unsigned long) printsh.ss,
             printsh.ssorient);
    }
  }

  if (b->vararea) {
    printf("    Area constraint:  %.4g\n", areabound(*t));
  }
}

/*****************************************************************************/
/*                                                                           */
/*  printsubseg()   Print out the details of an oriented subsegment.         */
/*                                                                           */
/*  I originally wrote this procedure to simplify debugging; it can be       */
/*  called directly from the debugger, and presents information about an     */
/*  oriented subsegment in digestible form.  It's also used when the highest */
/*  level of verbosity (`-VVV') is specified.                                */
/*                                                                           */
/*****************************************************************************/

void printsubseg(struct mesh *m, struct behavior *b, struct osub *s)
{
  struct osub printsh;
  struct otri printtri;
  vertex printvertex;

  printf("subsegment x%lx with orientation %d and mark %d:\n",
         (unsigned long) s->ss, s->ssorient, mark(*s));
  sdecode(s->ss[0], printsh);
  if (printsh.ss == m->dummysub) {
    printf("    [0] = No subsegment\n");
  } else {
    printf("    [0] = x%lx  %d\n", (unsigned long) printsh.ss,
           printsh.ssorient);
  }
  sdecode(s->ss[1], printsh);
  if (printsh.ss == m->dummysub) {
    printf("    [1] = No subsegment\n");
  } else {
    printf("    [1] = x%lx  %d\n", (unsigned long) printsh.ss,
           printsh.ssorient);
  }

  sorg(*s, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Origin[%d] = NULL\n", 2 + s->ssorient);
  else
    printf("    Origin[%d] = x%lx  (%.12g, %.12g)\n",
           2 + s->ssorient, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);
  sdest(*s, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Dest  [%d] = NULL\n", 3 - s->ssorient);
  else
    printf("    Dest  [%d] = x%lx  (%.12g, %.12g)\n",
           3 - s->ssorient, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);

  decode(s->ss[6], printtri);
  if (printtri.tri == m->dummytri) {
    printf("    [6] = Outer space\n");
  } else {
    printf("    [6] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }
  decode(s->ss[7], printtri);
  if (printtri.tri == m->dummytri) {
    printf("    [7] = Outer space\n");
  } else {
    printf("    [7] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }

  segorg(*s, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Segment origin[%d] = NULL\n", 4 + s->ssorient);
  else
    printf("    Segment origin[%d] = x%lx  (%.12g, %.12g)\n",
           4 + s->ssorient, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);
  segdest(*s, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Segment dest  [%d] = NULL\n", 5 - s->ssorient);
  else
    printf("    Segment dest  [%d] = x%lx  (%.12g, %.12g)\n",
           5 - s->ssorient, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);
}
/*============================================================================*/


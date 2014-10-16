/*============================================================================+/
 | Name: matlib_io.c
 |
/+============================================================================*/

#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>

#define NDEBUG

#include "matlib.h"
#include "assert.h"

/*============================================================================*/

void matlib_xmwrite_csv(char* file_name, matlib_xm M)
{

    FILE *fp = fopen(file_name, "w+");
#if 1
    if(fp==NULL)
    {
        term_exec("failed to create the file '%s'", file_name );
    }
#endif

    matlib_index i, j, col_st, row_st;
    if(M.order == MATLIB_COL_MAJOR)
    {
        col_st = 1;
        row_st = M.lenc;
    }
    else if(M.order == MATLIB_ROW_MAJOR)
    {
        col_st = M.lenr;
        row_st = 1;
    }
    else
    {
        term_exec( "Storage order unknown (order: %d)", M.order);
    }

    for (i=0; i<M.lenc; i++)
    {
        for (j=0; j<M.lenr-1; j++)
        {
            fprintf(fp, "% 0.16f\t", M.elem_p[i*col_st+j*row_st]);
        }
        fprintf(fp, "% 0.16f\n", M.elem_p[i*col_st+j*row_st]);
    }
    fclose(fp);
}
void matlib_zmwrite_csv(char* file_name, matlib_zm M)
{

    FILE *fp = fopen(file_name, "w+");
    if(fp==NULL)
    {
        term_exec("failed to create the file '%s'", file_name );
    }

    matlib_index i, j, col_st, row_st;
    if(M.order == MATLIB_COL_MAJOR)
    {
        col_st = 1;
        row_st = M.lenc;
    }
    else if(M.order == MATLIB_ROW_MAJOR)
    {
        col_st = M.lenr;
        row_st = 1;
    }
    else
    {
        term_exec( "Storage order unknown (order: %d)", M.order);
    }


    for (i=0; i<M.lenc; i++)
    {
        for (j=0; j<M.lenr-1; j++)
        {
            fprintf(fp, "% 0.16f%+0.16fi\t", M.elem_p[i*col_st+j*row_st]);
        }
        fprintf(fp, "% 0.16f%+0.16fi\n", M.elem_p[i*col_st+j*row_st]);
    }
    fclose(fp);
}

void matlib_xvwrite_csv
(
    char*        file_name, 
    matlib_index n,
    matlib_xv    v[n]
)
{
    FILE *fp = fopen(file_name, "w+");
    if(fp==NULL)
    {
        term_exec("failed to create the file '%s'", file_name );
    }

    matlib_index i, j;

    for (i=0; i<v->len; i++)
    {
        for (j=0; j<n-1; j++)
        {
            fprintf(fp, "% 0.16f\t", v[j].elem_p[i]);
        }
        fprintf(fp, "% 0.16f\n", v[j].elem_p[i]);
    }
    fclose(fp);
}

void matlib_zvwrite_csv
(
    char*        file_name, 
    matlib_index n,
    matlib_zv    v[n]
)
{

    FILE *fp = fopen(file_name, "w+");
    if(fp==NULL)
    {
        term_exec("failed to create the file '%s'", file_name );
    }

    matlib_index i, j;

    for (i=0; i<v->len; i++)
    {
        for (j=0; j<n-1; j++)
        {
            fprintf(fp, "% 0.16f%+0.16fi\t", v[j].elem_p[i]);
        }
        fprintf(fp, "% 0.16f%+0.16fi\n", v[j].elem_p[i]);
    }
    fclose(fp);
}

void matlib_xzvwrite_csv
(
    char*        file_name, 
    matlib_index m,
    matlib_xv    u[m],
    matlib_index n,
    matlib_zv    v[n]
)
{
    FILE *fp = fopen(file_name, "w+");
    if(fp==NULL)
    {
        term_exec("failed to create the file '%s'", file_name );
    }

    matlib_index i, j, k;

    for (i=0; i<u->len; i++)
    {
        for (j=0; j<m-1; j++)
        {
            fprintf(fp, "% 0.16f\t", u[j].elem_p[i]);
        }
        if(n>0)
        {
            fprintf(fp, "% 0.16f\t", u[j].elem_p[i]);
            for (k=0; k<n-1; k++)
            {
                fprintf(fp, "% 0.16f%+0.16fi\t", v[k].elem_p[i]);
            }
            fprintf(fp, "% 0.16f%+0.16fi\n", v[k].elem_p[i]);
        }
        else
        {
            fprintf(fp, "% 0.16f\n", u[j].elem_p[i]);
        }
    }
    fclose(fp);
}


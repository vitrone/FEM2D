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

//#define NDEBUG

#include "matlib.h"
#include "assert.h"

/*============================================================================*/
matlib_err matlib_io_create(matlib_index len, matlib_io_t* mp)
{

    matlib_index mcnt = 0;
    mp->len = len;
    errno = 0;
    
    mp->data_p = calloc(len, sizeof(void*));
    err_check( (mp->data_p == NULL), clean_up, 
               "%s: Memory allocation failed!", strerror(errno));
    mcnt++;
    errno = 0;
    mp->format = calloc(len, sizeof(MATLIB_IO_FORMAT));
    err_check( (mp->format == NULL), clean_up, 
               "%s: Memory allocation failed!", strerror(errno));
    mcnt++;

    errno = 0;
    mp->size = calloc(len + 1, sizeof(size_t));
    err_check( (mp->size == NULL), clean_up, 
               "%s: Memory allocation failed!", strerror(errno));

    for(matlib_index i = 0; i < len; i++)
    {
        mp->data_p[i] = NULL;
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    if (mcnt == 2)
    {
        matlib_free(mp->format);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(mp->data_p);
    }

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}


matlib_err matlib_io_size(matlib_io_t* mp)
{
    size_t mem_size = 0;

    err_check(    (mp->data_p == NULL) 
               || (mp->format == NULL)
               || (mp->size   == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    matlib_index i;
    mp->size[0] = MATLIB_INDEX_SIZE;

    for (i = 0; i < mp->len; i++)
    {
        mem_size = 0;
        err_check( mp->format[i] == MATLIB_FORMAT_UNKNOWN, clean_up,
                   "%s", "Unknown matrix format encountered!");

        if (mp->format[i] == MATLIB_XV)
        {
            matlib_xv* V = (matlib_xv*)mp->data_p[i];
            mem_size = MATLIB_INDEX_SIZE + MATLIB_ENUM_SIZE;
            mem_size += (V->len * MATLIB_REAL_SIZE);
        }
        else if (mp->format[i] == MATLIB_ZV)
        {
            matlib_zv* V = (matlib_zv*)mp->data_p[i];
            mem_size = MATLIB_INDEX_SIZE + MATLIB_ENUM_SIZE;
            mem_size += (V->len * MATLIB_COMPLEX_SIZE);
        }
        else if (mp->format[i] == MATLIB_DEN_XM)
        {
            matlib_xm* M = (matlib_xm*)mp->data_p[i];
            mem_size = 2 * MATLIB_INDEX_SIZE + 3 * MATLIB_ENUM_SIZE;
            mem_size += (M->lenc * M->lenr * MATLIB_REAL_SIZE); 
        }
        else if (mp->format[i] == MATLIB_DEN_ZM)
        {
            matlib_zm* M = (matlib_zm*)mp->data_p[i];
            mem_size = 2 * MATLIB_INDEX_SIZE + 3 * MATLIB_ENUM_SIZE;
            mem_size += (M->lenc * M->lenr * MATLIB_COMPLEX_SIZE); 
        }
        else if (mp->format[i] == MATLIB_NV)
        {
            matlib_nv* V = (matlib_nv*)mp->data_p[i];
            mem_size = MATLIB_INDEX_SIZE + MATLIB_ENUM_SIZE;
            mem_size += (V->len * MATLIB_INT_SIZE);
        }
        else
        {
            mp->format[i] = MATLIB_FORMAT_UNKNOWN;
        }
        err_check( mp->format[i] == MATLIB_FORMAT_UNKNOWN, clean_up, 
                   "%s", "Unknown matrix format encountered!");
        mp->size[i+1] = mem_size;
    }

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;

}
/*============================================================================*/

char matlib_io_order_char(MATLIB_ORDER order_enum)
{
    debug_enter( "Order: %s", MATLIB_ORDER_ENUM2STR(order_enum));
    char tmpc;

    switch (order_enum)
    {
        case MATLIB_COL_MAJOR:
            tmpc = 'A';
            break;

        case MATLIB_ROW_MAJOR:
            tmpc = 'B';
            break;

        default:
            tmpc = 'U';
    }
    debug_exit( "Order character: %c", tmpc);
    return tmpc;
}
MATLIB_ORDER matlib_io_order_enum(char orderc)
{
    debug_enter("Order character: %c", orderc);
    MATLIB_ORDER tmp_enum;
    switch (orderc)
    {
        case 'A':
            tmp_enum = MATLIB_COL_MAJOR;
            break;

        case 'B':
            tmp_enum = MATLIB_ROW_MAJOR;
            break;

        case 'U':
            tmp_enum = MATLIB_ORDER_UNKNOWN;

        default:
            tmp_enum = MATLIB_ORDER_UNKNOWN;
    }
    debug_exit( "Order %s", MATLIB_ORDER_ENUM2STR(tmp_enum));
    return tmp_enum;
}

char matlib_io_op_char(MATLIB_OP op_enum)
{
    char tmpc;
    switch (op_enum)
    {
        case MATLIB_NO_TRANS:
            tmpc = 'A';
            break;

        case MATLIB_TRANS:
            tmpc = 'B';
            break;
        case MATLIB_CONJ_TRANS:
            tmpc = 'C';
            break;
        default:
            tmpc = 'U';
    }
    return tmpc;
}
MATLIB_OP matlib_io_op_enum(char opc)
{
    MATLIB_OP tmp_enum;
    switch (opc)
    {
        case 'A': 
            tmp_enum = MATLIB_NO_TRANS;
            break;

        case 'B': 
            tmp_enum = MATLIB_TRANS;
            break;
        case 'C': 
            tmp_enum = MATLIB_CONJ_TRANS;
            break;
        default:
            tmp_enum = MATLIB_OP_UNKNOWN;
    }
    return tmp_enum;
}
char matlib_io_format_char(MATLIB_IO_FORMAT format_enum)
{
    char tmpc;
    switch (format_enum)
    {
        case MATLIB_XV:
            tmpc = 'A';
            break;

        case MATLIB_ZV:
            tmpc = 'B';
            break;
        case MATLIB_DEN_XM:
            tmpc = 'C';
            break;
        case MATLIB_DEN_ZM:
            tmpc = 'D';
            break;
        case MATLIB_NV:
            tmpc = 'E';
            break;
        default:
            tmpc = 'U';
    }
    return tmpc;
}

MATLIB_IO_FORMAT matlib_io_format_enum(char formatc)
{
    MATLIB_IO_FORMAT tmp_enum;
    switch (formatc)
    {
        case 'A': 
            tmp_enum = MATLIB_XV;
            break;

        case 'B': 
            tmp_enum = MATLIB_ZV;
            break;
        case 'C': 
            tmp_enum = MATLIB_DEN_XM;
            break;
        case 'D': 
            tmp_enum = MATLIB_DEN_ZM;
            break;
        case 'E': 
            tmp_enum = MATLIB_NV;
            break;
        case 'U': 
            tmp_enum = MATLIB_FORMAT_UNKNOWN;
            break;
        default:
            tmp_enum = MATLIB_FORMAT_UNKNOWN;
    }
    return tmp_enum;
}

/*============================================================================*/
matlib_err matlib_io_nv_write
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    matlib_index i = data_index;
    size_t nr_elems;

    char formatc = matlib_io_format_char(mp->format[i]);
    err_check( formatc == 'U', clean_up, 
               "%s", "Unknown matrix format encountered!");
    matlib_nv* V = (matlib_nv*)mp->data_p[i];
    nr_elems = fwrite( &formatc, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( &(V->len), MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( V->elem_p, MATLIB_INT_SIZE, V->len, fp);
    err_check( nr_elems != V->len, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, V->len);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

matlib_err matlib_io_xv_write
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    matlib_index i = data_index;
    size_t nr_elems;

    char formatc = matlib_io_format_char(mp->format[i]);
    err_check( formatc == 'U', clean_up, 
               "%s", "Unknown matrix format encountered!");
    matlib_xv* V = (matlib_xv*)mp->data_p[i];
    nr_elems = fwrite( &formatc, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( &(V->len), MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( V->elem_p, MATLIB_REAL_SIZE, V->len, fp);
    err_check( nr_elems != V->len, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, V->len);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

matlib_err matlib_io_zv_write
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    matlib_index i = data_index;
    size_t nr_elems;

    char formatc = matlib_io_format_char(mp->format[i]);
    err_check( formatc == 'U', clean_up, 
               "%s", "Unknown matrix format encountered!");
    matlib_zv* V = (matlib_zv*)mp->data_p[i];
    nr_elems = fwrite( &formatc, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( &(V->len), MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( V->elem_p, MATLIB_COMPLEX_SIZE, V->len, fp);
    err_check( nr_elems != V->len, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, V->len);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

matlib_err matlib_io_xm_write
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    size_t nr_elems;
    matlib_index i = data_index;
    char formatc = matlib_io_format_char(mp->format[i]);
    err_check( formatc == 'U', clean_up, 
               "%s", "Unknown matrix format encountered!");

    matlib_xm* M = (matlib_xm*)mp->data_p[i];
    nr_elems = fwrite( &formatc, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( &(M->lenc), MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( &(M->lenr), MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    char order = matlib_io_order_char(M->order);
    char op    = matlib_io_op_char(M->op);
    nr_elems = fwrite( &(order), MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( &(op), MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( M->elem_p, MATLIB_REAL_SIZE, M->lenc * M->lenr, fp);
    err_check( nr_elems != (M->lenc * M->lenr), clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, M->lenc * M->lenr);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

matlib_err matlib_io_zm_write
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    size_t nr_elems;
    matlib_index i = data_index;
    char formatc = matlib_io_format_char(mp->format[i]);
    err_check( formatc == 'U', clean_up, 
               "%s", "Unknown matrix format encountered!");

    matlib_zm* M = (matlib_zm*)mp->data_p[i];
    nr_elems = fwrite( &formatc, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( &(M->lenc), MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( &(M->lenr), MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    char order = matlib_io_order_char(M->order);
    char op    = matlib_io_op_char(M->op);
    nr_elems = fwrite( &(order), MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( &(op), MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);
    nr_elems = fwrite( M->elem_p, MATLIB_COMPLEX_SIZE, M->lenc * M->lenr, fp);
    err_check( nr_elems != (M->lenc * M->lenr), clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, M->lenc * M->lenr);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}
/*============================================================================*/

matlib_err matlib_io_fwrite
(
    matlib_io_t* mp,
    char* file_name
)
{
    err_check(    (mp->data_p == NULL) 
               || (mp->format == NULL)
               || (mp->size   == NULL)
               || (file_name  == NULL), clean_up, 
               "%s", "Null pointer(s) encountered!");
    errno = 0;
    FILE* fp;
    fp = fopen(file_name, "wb+");
    err_check( fp == NULL, clean_up, 
               "%s: Failed to open the file!", strerror(errno));

    size_t nr_elems;
    nr_elems = fwrite(&(mp->len), MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, 1);

    matlib_err error = matlib_io_size(mp);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to computed size!");
    
    nr_elems = fwrite(&mp->size[1], sizeof(size_t), mp->len, fp);
    err_check( nr_elems != mp->len, clean_up, 
               "A writing error occurred "
               "(nr. elements written: %d, total nr. elements: %d)!",
               nr_elems, mp->len);
    matlib_index i;
    char formatc;

    for (i = 0; i < mp->len; i++)
    {
        err_check( mp->format[i] == MATLIB_FORMAT_UNKNOWN, clean_up, 
                   "%s", "Unknown matrix format encountered!");
        
        if (mp->format[i] == MATLIB_XV)
        {
            error = matlib_io_xv_write(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to write real vector (data index: %d)!", i);
        }
        else if (mp->format[i] == MATLIB_ZV)
        {
            error = matlib_io_zv_write(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to write complex vector (data index: %d)!", i);
        }
        else if (mp->format[i] == MATLIB_DEN_XM)
        {
            error = matlib_io_xm_write(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to write real matrix (data index: %d)!", i);
        }
        else if (mp->format[i] == MATLIB_DEN_ZM)
        {
            error = matlib_io_zm_write(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to write complex matrix (data index: %d)!", i);
        }
        else if (mp->format[i] == MATLIB_NV)
        {
            error = matlib_io_nv_write(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to write integer vector (data index: %d)!", i);
        }
        else
        {
            mp->format[i] = MATLIB_FORMAT_UNKNOWN;
        }
        err_check( mp->format[i] == MATLIB_FORMAT_UNKNOWN, clean_up, 
                   "Unknown matrix format encountered (data index: %d)!", i);
    }
    fclose(fp);
    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;
clean_up:
    fclose(fp);
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

/*============================================================================*/
matlib_err matlib_io_nv_read
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    size_t offset = MATLIB_ENUM_SIZE;
    matlib_index r = fseek(fp, offset, SEEK_CUR);
    err_check( r != 0, clean_up, 
               "Error occurred moving within the file (nr bytes: %zu)!",
               offset);

    matlib_index len, nr_elems;
    matlib_index i = data_index;

    matlib_index mcnt = 0;
    matlib_err error;

    nr_elems = fread( &len, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);

    errno = 0;
    mp->data_p[i] = calloc( 1, sizeof(matlib_nv));
    err_check( mp->data_p[i] == NULL, clean_up, 
               "%s: Failed to allocate memory!", strerror(errno));
    mcnt++;

    matlib_nv* nv_ptr = (matlib_nv*)mp->data_p[i];

    error = matlib_create_nv(len, nv_ptr);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to allocate memory!");
    mcnt++;

    nr_elems = fread( nv_ptr->elem_p, MATLIB_INT_SIZE, len, fp);
    err_check( nr_elems != len, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, len);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    
    if (mcnt == 2)
    {
        matlib_free(nv_ptr->elem_p);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(mp->data_p[i]);
    }

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}


matlib_err matlib_io_xv_read
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    size_t offset = MATLIB_ENUM_SIZE;
    matlib_index r = fseek(fp, offset, SEEK_CUR);
    err_check( r != 0, clean_up, 
               "Error occurred moving within the file (nr bytes: %zu)!",
               offset);

    matlib_index len, nr_elems;
    matlib_index i = data_index;

    matlib_index mcnt = 0;
    matlib_err error;

    nr_elems = fread( &len, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);

    errno = 0;
    mp->data_p[i] = calloc( 1, sizeof(matlib_xv));
    err_check( mp->data_p[i] == NULL, clean_up, 
               "%s: Failed to allocate memory!", strerror(errno));
    mcnt++;

    matlib_xv* xv_ptr = (matlib_xv*)mp->data_p[i];

    error = matlib_create_xv(len, xv_ptr, MATLIB_COL_VECT);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to allocate memory!");
    mcnt++;

    nr_elems = fread( xv_ptr->elem_p, MATLIB_REAL_SIZE, len, fp);
    err_check( nr_elems != len, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, len);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    
    if (mcnt == 2)
    {
        matlib_free(xv_ptr->elem_p);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(mp->data_p[i]);
    }

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

matlib_err matlib_io_zv_read
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    size_t offset = MATLIB_ENUM_SIZE;
    matlib_index r = fseek(fp, offset, SEEK_CUR);
    err_check( r != 0, clean_up, 
               "Error occurred moving within the file (nr bytes: %zu)!",
               offset);

    matlib_index len, nr_elems;
    matlib_index i = data_index;

    matlib_index mcnt = 0;
    matlib_err error;

    nr_elems = fread( &len, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);

    errno = 0;
    mp->data_p[i] = calloc( 1, sizeof(matlib_zv));
    err_check( mp->data_p[i] == NULL, clean_up, 
               "%s: Failed to allocate memory!", strerror(errno));
    mcnt++;

    matlib_zv* zv_ptr = (matlib_zv*)mp->data_p[i];

    error = matlib_create_zv(len, zv_ptr, MATLIB_COL_VECT);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to allocate memory!");
    mcnt++;

    nr_elems = fread( zv_ptr->elem_p, MATLIB_COMPLEX_SIZE, len, fp);
    err_check( nr_elems != len, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, len);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    
    if (mcnt == 2)
    {
        matlib_free(zv_ptr->elem_p);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(mp->data_p[i]);
    }

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}
/*============================================================================*/

matlib_err matlib_io_xm_read
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    size_t offset = MATLIB_ENUM_SIZE;
    matlib_index r = fseek(fp, offset, SEEK_CUR);
    err_check( r != 0, clean_up, 
               "Error occurred moving within the file (nr bytes: %zu)!",
               offset);

    matlib_index lenc, lenr, nr_elems;
    matlib_index i = data_index;

    matlib_index mcnt = 0;
    matlib_err error;

    nr_elems = fread( &lenc, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    
    nr_elems = fread( &lenr, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    debug_body("size: %d-by-%d", lenc, lenr);

    char order;
    nr_elems = fread( &order, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    MATLIB_ORDER order_enum = matlib_io_order_enum(order);
    debug_body( "Order: %s", MATLIB_ORDER_ENUM2STR(order_enum));

    char op;
    nr_elems = fread( &op, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    MATLIB_OP op_enum = matlib_io_op_enum(op);
    debug_body("Operation: %s", MATLIB_OP_ENUM2STR(op_enum));
    
    errno = 0;
    mp->data_p[i] = calloc(1, sizeof(matlib_xm));
    err_check( mp->data_p[i] == NULL, clean_up, 
               "%s: Failed to allocate memory!", strerror(errno));
    mcnt++;
    
    matlib_xm* xm_ptr = (matlib_xm*)mp->data_p[i];
    error = matlib_create_xm( lenc, lenr, xm_ptr, order_enum, op_enum);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to allocate memory!");
    mcnt++;

    nr_elems = fread( xm_ptr->elem_p, MATLIB_REAL_SIZE, lenc * lenr, fp);
    err_check( nr_elems != (lenc * lenr), clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, (lenc * lenr));

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    
    if (mcnt == 2)
    {
        matlib_free(xm_ptr->elem_p);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(mp->data_p[i]);
    }

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;

} /* End of matlib_io_xm_read */ 

matlib_err matlib_io_zm_read
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp
)
{
    debug_enter("data index: %d", data_index);
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    size_t offset = MATLIB_ENUM_SIZE;
    matlib_index r = fseek(fp, offset, SEEK_CUR);
    err_check( r != 0, clean_up, 
               "Error occurred moving within the file (nr bytes: %zu)!",
               offset);

    matlib_index lenc, lenr, nr_elems;
    matlib_index i = data_index;

    matlib_index mcnt = 0;
    matlib_err error;

    nr_elems = fread( &lenc, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    
    nr_elems = fread( &lenr, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    debug_body("size: %d-by-%d", lenc, lenr);

    char order;
    nr_elems = fread( &order, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    MATLIB_ORDER order_enum = matlib_io_order_enum(order);
    debug_body( "Order: %s", MATLIB_ORDER_ENUM2STR(order_enum));

    char op;
    nr_elems = fread( &op, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    MATLIB_OP op_enum = matlib_io_op_enum(op);
    debug_body("Operation: %s", MATLIB_OP_ENUM2STR(op_enum));
    
    errno = 0;
    mp->data_p[i] = calloc(1, sizeof(matlib_zm));
    err_check( mp->data_p[i] == NULL, clean_up, 
               "%s: Failed to allocate memory!", strerror(errno));
    mcnt++;
    
    matlib_zm* zm_ptr = (matlib_zm*)mp->data_p[i];
    error = matlib_create_zm( lenc, lenr, zm_ptr, order_enum, op_enum);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to allocate memory!");
    mcnt++;

    nr_elems = fread( zm_ptr->elem_p, MATLIB_COMPLEX_SIZE, lenc * lenr, fp);
    err_check( nr_elems != (lenc * lenr), clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, (lenc * lenr));

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    
    if (mcnt == 2)
    {
        matlib_free(zm_ptr->elem_p);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(mp->data_p[i]);
    }

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}
/*============================================================================*/
matlib_err matlib_io_getinfo
(
    matlib_io_t* mp,
    FILE* fp
)
{
    err_check( (mp == NULL) || (fp == NULL), clean_up,
               "%s", "Null pointer(s) encountered!");

    matlib_index mcnt = 0;
    size_t nr_elems;
    matlib_index alen;
    nr_elems = fread( &(alen), MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);

    matlib_err error = matlib_io_create(alen, mp);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to allocate memory!");
    mcnt++;

    mp->size[0] = (mp->len + 1) * MATLIB_INDEX_SIZE;
    nr_elems = fread( &(mp->size[1]), MATLIB_INDEX_SIZE, mp->len, fp);
    err_check( nr_elems != mp->len, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, mp->len);

    matlib_index i, r;
    char formatc;
    size_t offset = mp->size[0];
    for (i = 0; i < mp->len; i++)
    {
        nr_elems = fread( &formatc, MATLIB_ENUM_SIZE, 1, fp);
        err_check( nr_elems != 1, clean_up, 
                   "A reading error occurred "
                   "(nr. elements read: %d, total nr. elements: %d)!",
                   nr_elems, 1);
        debug_body("Format: %c", formatc);

        mp->format[i] = matlib_io_format_enum(formatc);
        err_check( mp->format[i] == MATLIB_FORMAT_UNKNOWN, clean_up, 
                   "Unknown matrix format encountered (data index: %d)!", i);

        offset = (mp->size[i+1] - MATLIB_ENUM_SIZE);
        debug_body("Offset: %zu", offset);
        r = fseek(fp, offset, SEEK_CUR);
        err_check( r != 0, clean_up, 
                   "Error occurred moving within the file (nr bytes: %zu)!",
                   offset);
    }
    
    /* Return to the begining of the file */ 
    r = fseek(fp, 0L, SEEK_SET);
    err_check( r != 0, clean_up, 
               "%s", "Failed to return to the begining of the file!");

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;

}/* End of matlib_io_getinfo */ 


matlib_err matlib_io_getelem
(
    matlib_io_t* mp,
    matlib_index data_index,
    FILE* fp
)
{
    err_check( (mp == NULL) || (fp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( (mp->size == NULL) || (mp->format == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    matlib_err error;
    matlib_index i;
    matlib_int j;

    size_t offset = 0;
    for (i = 0; i < (data_index + 1); i++)
    {
        offset += (mp->size[i]);
    }
    debug_body("offset: %zu", offset);
    i = fseek(fp, offset, SEEK_SET);
    err_check( i != 0, clean_up, 
               "Error occurred moving within the file (nr bytes: %zu)!",
               mp->size[0]);

    i = data_index;
    if (mp->format[i] == MATLIB_XV)
    {
        error = matlib_io_xv_read(mp, i, fp);
        err_check( error == MATLIB_FAILURE, clean_up, 
                   "Failed to read real vector (data index: %d)!", i);
    }
    else if (mp->format[i] == MATLIB_ZV)
    {
        error = matlib_io_zv_read(mp, i, fp);
        err_check( error == MATLIB_FAILURE, clean_up, 
                   "Failed to read complex vector (data index: %d)!", i);
    }
    else if (mp->format[i] == MATLIB_DEN_XM)
    {
        error = matlib_io_xm_read(mp, i, fp);
        err_check( error == MATLIB_FAILURE, clean_up, 
                   "Failed to read real matrix (data index: %d)!", i);
    }
    else if (mp->format[i] == MATLIB_DEN_ZM)
    {
        error = matlib_io_zm_read(mp, i, fp);
        err_check( error == MATLIB_FAILURE, clean_up, 
                   "Failed to read complex matrix (data index: %d)!", i);
    }
    else if (mp->format[i] == MATLIB_NV)
    {
        error = matlib_io_nv_read(mp, i, fp);
        err_check( error == MATLIB_FAILURE, clean_up, 
                   "Failed to read integer matrix (data index: %d)!", i);
    }
    else
    {
       mp->format[i] == MATLIB_FORMAT_UNKNOWN;
    }
    err_check( mp->format[i] == MATLIB_FORMAT_UNKNOWN, clean_up, 
               "Unknown matrix format encountered (data index: %d)!", i);


    debug_exit("Read data (index: %d), Exit Status: %s", i, "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:

    matlib_io_elemfree(mp, i);
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}
/*============================================================================*/
matlib_err matlib_io_fread
(
    matlib_io_t* mp,
    char* file_name
)
{
    err_check(    (mp        == NULL)
               || (file_name == NULL), clean_up, 
               "%s", "Null pointer encountered!");
    errno = 0;
    FILE* fp;
    fp = fopen(file_name, "rb");
    err_check( fp == NULL, clean_up, 
               "%s: Failed to open the file!", strerror(errno));

    matlib_err error = matlib_io_getinfo(mp, fp);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to get info from file!");

    matlib_index i;
    matlib_int j;

    i = fseek(fp, mp->size[0], SEEK_SET);
    err_check( i != 0, clean_up, 
               "Error occurred moving within the file (nr bytes: %zu)!",
               mp->size[0]);

    for (i = 0; i < mp->len; i++)
    {
        if (mp->format[i] == MATLIB_XV)
        {
            error = matlib_io_xv_read(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to read real vector (data index: %d)!", i);
        }
        else if (mp->format[i] == MATLIB_ZV)
        {
            error = matlib_io_zv_read(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to read complex vector (data index: %d)!", i);
        }
        else if (mp->format[i] == MATLIB_DEN_XM)
        {
            error = matlib_io_xm_read(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to read real matrix (data index: %d)!", i);
        }
        else if (mp->format[i] == MATLIB_DEN_ZM)
        {
            error = matlib_io_zm_read(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to read complex matrix (data index: %d)!", i);
        }
        else if (mp->format[i] == MATLIB_NV)
        {
            error = matlib_io_nv_read(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to read integer matrix (data index: %d)!", i);
        }
        else
        {
           mp->format[i] == MATLIB_FORMAT_UNKNOWN;
        }
        err_check( mp->format[i] == MATLIB_FORMAT_UNKNOWN, clean_up, 
                   "Unknown matrix format encountered (data index: %d)!", i);
    }

    fclose(fp);
    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:

    j = (matlib_int)i;
    for (; j>=0; j--)
    {
        matlib_io_elemfree(mp, j);
    }
    fclose(fp);
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}
/*============================================================================*/
matlib_err matlib_io_nv_extread
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp,
    matlib_err (*create)(matlib_index, matlib_nv *, void*),
    void (*destroy)(void*),
    void* extdata_p
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    size_t offset = MATLIB_ENUM_SIZE;
    matlib_index r = fseek(fp, offset, SEEK_CUR);
    err_check( r != 0, clean_up, 
               "Error occurred moving within the file (nr bytes: %zu)!",
               offset);

    matlib_index len, nr_elems;
    matlib_index i = data_index;

    matlib_index mcnt = 0;
    matlib_err error;

    nr_elems = fread( &len, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);

    errno = 0;
    mp->data_p[i] = calloc( 1, sizeof(matlib_nv));
    err_check( mp->data_p[i] == NULL, clean_up, 
               "%s: Failed to allocate memory!", strerror(errno));
    mcnt++;

    matlib_nv* nv_ptr = (matlib_nv*)mp->data_p[i];

    error = create(len, nv_ptr, extdata_p);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to allocate memory!");
    mcnt++;

    nr_elems = fread( nv_ptr->elem_p, MATLIB_INT_SIZE, len, fp);
    err_check( nr_elems != len, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, len);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    
    if (mcnt == 2)
    {
        destroy(extdata_p);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(mp->data_p[i]);
    }

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}


matlib_err matlib_io_xv_extread
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp,
    matlib_err (*create)( matlib_index, 
                          matlib_xv *, 
                          MATLIB_VECT_T, 
                          void*),
    void (*destroy)(void*),
    void* extdata_p
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    size_t offset = MATLIB_ENUM_SIZE;
    matlib_index r = fseek(fp, offset, SEEK_CUR);
    err_check( r != 0, clean_up, 
               "Error occurred moving within the file (nr bytes: %zu)!",
               offset);

    matlib_index len, nr_elems;
    matlib_index i = data_index;

    matlib_index mcnt = 0;
    matlib_err error;

    nr_elems = fread( &len, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);

    errno = 0;
    mp->data_p[i] = calloc( 1, sizeof(matlib_xv));
    err_check( mp->data_p[i] == NULL, clean_up, 
               "%s: Failed to allocate memory!", strerror(errno));
    mcnt++;

    matlib_xv* xv_ptr = (matlib_xv*)mp->data_p[i];

    error = create(len, xv_ptr, MATLIB_COL_VECT, extdata_p);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to allocate memory!");
    mcnt++;

    nr_elems = fread( xv_ptr->elem_p, MATLIB_REAL_SIZE, len, fp);
    err_check( nr_elems != len, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, len);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    
    if (mcnt == 2)
    {
        destroy(extdata_p);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(mp->data_p[i]);
    }

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

matlib_err matlib_io_zv_extread
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp,
    matlib_err (*create)( matlib_index, 
                          matlib_zv *, 
                          MATLIB_VECT_T, 
                          void*),
    void (*destroy)(void*),
    void* extdata_p
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    size_t offset = MATLIB_ENUM_SIZE;
    matlib_index r = fseek(fp, offset, SEEK_CUR);
    err_check( r != 0, clean_up, 
               "Error occurred moving within the file (nr bytes: %zu)!",
               offset);

    matlib_index len, nr_elems;
    matlib_index i = data_index;

    matlib_index mcnt = 0;
    matlib_err error;

    nr_elems = fread( &len, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);

    errno = 0;
    mp->data_p[i] = calloc( 1, sizeof(matlib_zv));
    err_check( mp->data_p[i] == NULL, clean_up, 
               "%s: Failed to allocate memory!", strerror(errno));
    mcnt++;

    matlib_zv* zv_ptr = (matlib_zv*)mp->data_p[i];

    error = create(len, zv_ptr, MATLIB_COL_VECT, extdata_p);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to allocate memory!");
    mcnt++;

    nr_elems = fread( zv_ptr->elem_p, MATLIB_COMPLEX_SIZE, len, fp);
    err_check( nr_elems != len, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, len);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    
    if (mcnt == 2)
    {
        destroy(extdata_p);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(mp->data_p[i]);
    }

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}
/*============================================================================*/

matlib_err matlib_io_xm_extread
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp,
    matlib_err (*create)( matlib_index,
                          matlib_index, 
                          matlib_xm *, 
                          MATLIB_ORDER, 
                          MATLIB_OP, 
                          void*),
    void (*destroy)(void*),
    void* extdata_p
)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    size_t offset = MATLIB_ENUM_SIZE;
    matlib_index r = fseek(fp, offset, SEEK_CUR);
    err_check( r != 0, clean_up, 
               "Error occurred moving within the file (nr bytes: %zu)!",
               offset);

    matlib_index lenc, lenr, nr_elems;
    matlib_index i = data_index;

    matlib_index mcnt = 0;
    matlib_err error;

    nr_elems = fread( &lenc, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    
    nr_elems = fread( &lenr, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    debug_body("size: %d-by-%d", lenc, lenr);

    char order;
    nr_elems = fread( &order, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    MATLIB_ORDER order_enum = matlib_io_order_enum(order);
    debug_body( "Order: %s", MATLIB_ORDER_ENUM2STR(order_enum));

    char op;
    nr_elems = fread( &op, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    MATLIB_OP op_enum = matlib_io_op_enum(op);
    debug_body("Operation: %s", MATLIB_OP_ENUM2STR(op_enum));
    
    errno = 0;
    mp->data_p[i] = calloc(1, sizeof(matlib_xm));
    err_check( mp->data_p[i] == NULL, clean_up, 
               "%s: Failed to allocate memory!", strerror(errno));
    mcnt++;
    
    matlib_xm* xm_ptr = (matlib_xm*)mp->data_p[i];
    error = create( lenc, lenr, xm_ptr, order_enum, op_enum, extdata_p);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to allocate memory!");
    mcnt++;

    nr_elems = fread( xm_ptr->elem_p, MATLIB_REAL_SIZE, lenc * lenr, fp);
    err_check( nr_elems != (lenc * lenr), clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, (lenc * lenr));

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    
    if (mcnt == 2)
    {
        destroy(extdata_p);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(mp->data_p[i]);
    }

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;

} /* End of matlib_io_xm_read */ 

matlib_err matlib_io_zm_extread
(
    matlib_io_t* mp, 
    matlib_index data_index, 
    FILE* fp,
    matlib_err (*create)( matlib_index,
                          matlib_index, 
                          matlib_zm *, 
                          MATLIB_ORDER, 
                          MATLIB_OP, 
                          void*),
    void (*destroy)(void*),
    void* extdata_p
)
{
    debug_enter("data index: %d", data_index);
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    err_check( data_index >= mp->len, clean_up, 
               "Error occurred in accessing the data "
               "element (sought: %d, length: %d)!",
               data_index, mp->len);

    size_t offset = MATLIB_ENUM_SIZE;
    matlib_index r = fseek(fp, offset, SEEK_CUR);
    err_check( r != 0, clean_up, 
               "Error occurred moving within the file (nr bytes: %zu)!",
               offset);

    matlib_index lenc, lenr, nr_elems;
    matlib_index i = data_index;

    matlib_index mcnt = 0;
    matlib_err error;

    nr_elems = fread( &lenc, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    
    nr_elems = fread( &lenr, MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    debug_body("size: %d-by-%d", lenc, lenr);

    char order;
    nr_elems = fread( &order, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    MATLIB_ORDER order_enum = matlib_io_order_enum(order);
    debug_body( "Order: %s", MATLIB_ORDER_ENUM2STR(order_enum));

    char op;
    nr_elems = fread( &op, MATLIB_ENUM_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);
    MATLIB_OP op_enum = matlib_io_op_enum(op);
    debug_body("Operation: %s", MATLIB_OP_ENUM2STR(op_enum));
    
    errno = 0;
    mp->data_p[i] = calloc(1, sizeof(matlib_zm));
    err_check( mp->data_p[i] == NULL, clean_up, 
               "%s: Failed to allocate memory!", strerror(errno));
    mcnt++;
    
    matlib_zm* zm_ptr = (matlib_zm*)mp->data_p[i];
    error = create( lenc, lenr, zm_ptr, order_enum, op_enum, extdata_p);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to allocate memory!");
    mcnt++;

    nr_elems = fread( zm_ptr->elem_p, MATLIB_COMPLEX_SIZE, lenc * lenr, fp);
    err_check( nr_elems != (lenc * lenr), clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, (lenc * lenr));

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    
    if (mcnt == 2)
    {
        destroy(extdata_p);
        mcnt--;
    }
    if (mcnt == 1)
    {
        matlib_free(mp->data_p[i]);
    }

    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

#if 0
matlib_err matlib_io_mem_analyze
(
    matlib_io_t* mp,
    FILE* fp
)
{
    err_check( (mp == NULL) || (fp == NULL), clean_up,
               "%s", "Null pointer(s) encountered!");

    matlib_index mcnt = 0;
    size_t nr_elems;
    matlib_index alen;
    nr_elems = fread( &(alen), MATLIB_INDEX_SIZE, 1, fp);
    err_check( nr_elems != 1, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, 1);

    matlib_err error = matlib_io_create(alen, mp);
    err_check( error == MATLIB_FAILURE, clean_up, 
               "%s", "Failed to allocate memory!");
    mcnt++;

    mp->size[0] = (mp->len + 1) * MATLIB_INDEX_SIZE;
    nr_elems = fread( &(mp->size[1]), MATLIB_INDEX_SIZE, mp->len, fp);
    err_check( nr_elems != mp->len, clean_up, 
               "A reading error occurred "
               "(nr. elements read: %d, total nr. elements: %d)!",
               nr_elems, mp->len);

    matlib_index i, r;
    char formatc;
    size_t offset = mp->size[0];
    for (i = 0; i < mp->len; i++)
    {
        nr_elems = fread( &formatc, MATLIB_ENUM_SIZE, 1, fp);
        err_check( nr_elems != 1, clean_up, 
                   "A reading error occurred "
                   "(nr. elements read: %d, total nr. elements: %d)!",
                   nr_elems, 1);
        debug_body("Format: %c", formatc);

        mp->format[i] = matlib_io_format_enum(formatc);
        err_check( mp->format[i] == MATLIB_FORMAT_UNKNOWN, clean_up, 
                   "Unknown matrix format encountered (data index: %d)!", i);

        if (mp->format[i] == MATLIB_XV)
        {
            error = matlib_io_xv_read(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to read real vector (data index: %d)!", i);
        }
        else if (mp->format[i] == MATLIB_ZV)
        {
            error = matlib_io_zv_read(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to read complex vector (data index: %d)!", i);
        }
        else if (mp->format[i] == MATLIB_DEN_XM)
        {
            error = matlib_io_xm_read(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to read real matrix (data index: %d)!", i);
        }
        else if (mp->format[i] == MATLIB_DEN_ZM)
        {
            error = matlib_io_zm_read(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to read complex matrix (data index: %d)!", i);
        }
        else if (mp->format[i] == MATLIB_NV)
        {
            error = matlib_io_nv_read(mp, i, fp);
            err_check( error == MATLIB_FAILURE, clean_up, 
                       "Failed to read integer matrix (data index: %d)!", i);
        }




        offset = (mp->size[i+1] - MATLIB_ENUM_SIZE);
        debug_body("Offset: %zu", offset);
        r = fseek(fp, offset, SEEK_CUR);
        err_check( r != 0, clean_up, 
                   "Error occurred moving within the file (nr bytes: %zu)!",
                   offset);
    }
    
    /* Return to the begining of the file */ 
    r = fseek(fp, 0L, SEEK_SET);
    err_check( r != 0, clean_up, 
               "%s", "Failed to return to the begining of the file!");

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;

    

}
#endif

/*============================================================================*/

matlib_err matlib_io_elemfree
(
    matlib_io_t* mp,
    matlib_index data_index
)
{
    debug_enter("data index: %d", data_index);
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    matlib_index i = data_index;
    err_check( (mp->data_p[i] == NULL), clean_up, 
               "%s", "Null pointer encountered!");
    if (mp->format[i] == MATLIB_XV)
    {
        matlib_xv* xv_ptr = (matlib_xv*)mp->data_p[i];
        matlib_free(xv_ptr->elem_p);
        matlib_free(xv_ptr);

    }
    else if (mp->format[i] == MATLIB_ZV)
    {
        matlib_zv* zv_ptr = (matlib_zv*)mp->data_p[i];
        matlib_free(zv_ptr->elem_p);
        matlib_free(zv_ptr);
    }
    else if (mp->format[i] == MATLIB_DEN_XM)
    {
        matlib_xm* xm_ptr = (matlib_xm*)mp->data_p[i];
        matlib_free(xm_ptr->elem_p);
        //matlib_free(xm_ptr);
    }
    else if (mp->format[i] == MATLIB_DEN_ZM)
    {
        matlib_zm* zm_ptr = (matlib_zm*)mp->data_p[i];
        matlib_free(zm_ptr->elem_p);
        matlib_free(zm_ptr);
    }
    else
    {
       mp->format[i] == MATLIB_FORMAT_UNKNOWN;
    }
    err_check( mp->format[i] == MATLIB_FORMAT_UNKNOWN, clean_up, 
               "%s", "Unknown matrix format encountered!");

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

matlib_err matlib_io_free(matlib_io_t* mp)
{
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");

    matlib_free(mp->data_p);
    matlib_free(mp->size);
    matlib_free(mp->format);

    mp->data_p = NULL;
    mp->size   = NULL;
    mp->format = NULL;

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

matlib_err matlib_io_freeall(matlib_io_t* mp)
{
    debug_enter("%s", "Free all memory!");
    err_check( (mp == NULL), clean_up, 
               "%s", "Null pointer encountered!");
    for (matlib_index i = 0; i < mp->len; i++)
    {
        debug_enter("data index: %d", i);
        matlib_io_elemfree(mp, i);
    }
    matlib_io_free(mp);

    debug_exit("Exit Status: %s", "SUCCESS");
    return MATLIB_SUCCESS;

clean_up:
    debug_exit("Exit Status: %s", "FAILURE");
    return MATLIB_FAILURE;
}

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


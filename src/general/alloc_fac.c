/*

Copyright (c) 1991 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2005 by:
Eero Häkkinen
*/

/*
  Full control of all allocation 
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "memalloc.h"
#include "printmsg.h"

#include "stdafx.h"

/***************************************************************************/
float *gomp_ReallocateFloatVector(float *vec, size_t nh)
/* add more memory to a floating pointer array */
    /* name of array to get more memory */
    /* new length of array              */
    /* Leif Laaksonen 1990              */
/***************************************************************************/
{
    float *v;

    v = realloc(vec , (nh*sizeof(*v)));
    if( v == NULL) {
        gomp_PrintMessage("?gOpenMol ERROR");
        gomp_PrintMessage("** ERROR in gomp_AllocateFloatVector ... Reallocation failure ");
        return NULL;
    }

    return v;
}
/***************************************************************************/
double  *gomp_ReallocateDoubleVector(double  *vec, size_t nh) 
/* add more memory to a double pointer array
   const double  *vec;              name of array to get more memory 
   size_t nh;                   new length of array              
   Leif Laaksonen 1990              
*/
/***************************************************************************/
{
    double  *v;

    v = realloc(vec , (nh*sizeof(*v)));
    if( v == NULL) {
        gomp_PrintMessage("?gOpenMol ERROR");
        gomp_PrintMessage("** ERROR in gomp_AllocateDoubleVector ... Reallocation failure ");
        return NULL;
    }

    return v;
}

/***************************************************************************/
double  *gomp_AllocateDoubleVector(size_t nh)
/* reserve memory for double array */
    /* expected length of array       */
    /* Leif Laaksonen 1990            */
/***************************************************************************/
{
    double  *v;

    v = malloc(nh*sizeof(*v));
    if( v == NULL) {
        gomp_PrintMessage("?gOpenMol ERROR");
        gomp_PrintMessage("** ERROR in gomp_AllocateDoubleVector ... Allocation failure ");
        return NULL;
    }

    return v;
}


/***************************************************************************/
float *gomp_AllocateFloatVector(size_t nh)
/* reserve memory for float array */
    /* expected length of array       */
    /* Leif Laaksonen 1990            */
/***************************************************************************/
{
    float *v;

    v = malloc(nh*sizeof(*v));
    if( v == NULL) {
        gomp_PrintMessage("?gOpenMol ERROR");
        gomp_PrintMessage("** ERROR in gomp_AllocateFloatVector ... Allocation failure ");
        return NULL;
    }

    return v;
}

/***************************************************************************/
int *gomp_ReallocateIntVector(int *vec,size_t nh)
/* add more memory to an integer pointer array */
    /* name of array to get more memory */
    /* new length of array              */
    /* Leif Laaksonen 1990              */
/***************************************************************************/
{
    int *v;
    
    v = realloc(vec , (nh*sizeof(*v)));
    if( v == NULL) {
        gomp_PrintMessage("?gOpenMol ERROR");
        gomp_PrintMessage("** ERROR in gomp_AllocateIntVector ... Reallocation failure ");
        return NULL;
    }

    return v;
}

/***************************************************************************/
int *gomp_AllocateIntVector(size_t nh)
/* reserve memory for integer array */
    /* expected length of array         */
    /* Leif Laaksonen 1990              */
/***************************************************************************/
{
    int *v;

    v = malloc(nh*sizeof(*v));
    if( v == NULL) {
        gomp_PrintMessage("?gOpenMol ERROR");
        gomp_PrintMessage("** ERROR in gomp_AllocateIntVector ... Allocation failure ... \n\n");
        return NULL;
    }

    return v;
}


/***************************************************************************/
char *gomp_ReallocateCharVector(char *vec, size_t nh)
/* add more space to character array */
    /* array to get the space */
    /* expected length        */
    /* Leif Laaksonen 1990    */
/***************************************************************************/
{
    char *v;

    v = realloc(vec , (nh*sizeof(*v)));
    if( v == NULL) {
        gomp_PrintMessage("?gOpenMol ERROR");
        gomp_PrintMessage("** ERROR - in gomp_AllocateCharVector ... Reallocation failure");
        return NULL;
    }

    return v;
}

/***************************************************************************/
char *gomp_AllocateCharVector(size_t nh)
/* reserve memory for character array */
    /* length of character array          */
    /* Leif Laaksonen 1990                */
/***************************************************************************/
{
    char *v;

    v = malloc(nh*sizeof(*v));
    if( v == NULL) {
        gomp_PrintMessage("?gOpenMol ERROR");
        gomp_PrintMessage("** ERROR in cvector ... Allocation failure ...");
        return NULL;
    }

    return v;
}

/***************************************************************************/
void *gomp_AllocateVoidVector(size_t nh)
/* reserve memory for a void  array */
    /* expected length of array         */
    /* Eero Häkkinen 2003               */
/***************************************************************************/
{
    void *v;

    v = malloc( nh );
    if( v == NULL) {
        gomp_PrintMessage("?gOpenMol ERROR");
        gomp_PrintMessage("** ERROR in vvector ... Allocation failure ...");
        return NULL;
    }

    return v;
}

/***************************************************************************/
void *gomp_ReallocateVoidVector(void *vec, size_t nh)
/* add more memory to an void array */
    /* name of array to get more memory */
    /* new length of array              */
    /* Eero Häkkinen 2003               */
/***************************************************************************/
{
    short *v;

    v = realloc(vec , nh);
    if( v == NULL) {
        gomp_PrintMessage("?gOpenMol ERROR");
        gomp_PrintMessage("** ERROR in vvector ... Reallocation failure ");
        return NULL;
    }

    return v;
}

/***************************************************************************/
void gomp_FreeVector(void *vec)
/***************************************************************************/
{
    free(vec);
}

#define GET_DATA(vector, item_size) \
    ((DataVectorData*)(((char*)(vector)) - \
    gomp_DataVectorGetDataSize(item_size)))
/***************************************************************************/
void gomp_DataVectorCopyInitFunc(
    void  *data,
    size_t item_size,
    const  DataVectorHandler *pHandler)
/***************************************************************************/
{
    memcpy(data, pHandler->cdata, item_size);
}
/***************************************************************************/
void gomp_DataVectorZeroInitFunc(
    void  *data,
    size_t item_size,
    const DataVectorHandler *pHandler)
/***************************************************************************/
{
    memset(data, 0, item_size);
}
/***************************************************************************/
void *gomp_DataVectorCreateI(
    const DataVectorHandler *pHandler,
    size_t item_size,
    size_t size)
/***************************************************************************/
{
    char *vector;
    /* Alloc memory. */
    DataVectorData *data = gomp_AllocateVoidVector(
        gomp_DataVectorGetDataSize(item_size) + size * item_size);
    if ( data == NULL )
        return NULL;
    vector = ((char*)data) + gomp_DataVectorGetDataSize(item_size);

    /* Store properties. */
    data->total_size = data->size = size;
    data->handler    = *pHandler;
    /* Init array. */
    if ( data->handler.InitFunc ) {
        size_t i;
        for ( i = 0 ; i < size ; ++i )
            data->handler.InitFunc(
                vector + i * item_size, item_size, &data->handler);
    }

    return vector;
}
/***************************************************************************/
void *gomp_DataVectorCreateOrAppendI(
    void  *vector,
    const DataVectorHandler *pHandler,
    size_t item_size,
    size_t relative_size)
/***************************************************************************/
{
    if ( vector ) {
        /* Vector is already created. */
        return gomp_DataVectorSetSizeI(
            vector, item_size,
            GET_DATA(vector, item_size)->size + relative_size);
    }
    else
        return gomp_DataVectorCreateI(pHandler, item_size, relative_size);
}
/***************************************************************************/
void *gomp_DataVectorCreateOrAppendExactI(
    void  *vector,
    const DataVectorHandler *pHandler,
    size_t item_size,
    size_t relative_size)
/***************************************************************************/
{
    if (vector) {
        /* Vector is already created. */
        return gomp_DataVectorSetExactSizeI(
            vector, item_size,
            GET_DATA(vector, item_size)->size + relative_size);
    }
    else
        return gomp_DataVectorCreateI(pHandler, item_size, relative_size);
}
/***************************************************************************/
static void *DataVectorSetSizeI(
    void  *vector,
    size_t item_size,
    size_t size,
    size_t total_size)
/***************************************************************************/
{
    DataVectorData *data = GET_DATA(vector, item_size);

    if ( size < data->size ) {
        /* Delete unused items. */
        if ( data->handler.DestroyFunc ) {
            char *array = vector;
            size_t i;
            for ( i = size ; i < data->size ; i++ )
                data->handler.DestroyFunc(
                    array + i * item_size, item_size, &data->handler);
        }
        data->size = size;
    }

    if ( total_size != data->total_size ) {
        /* We need a new array. */
        DataVectorData *new_data =
            gomp_ReallocateVoidVector(
                data,
                gomp_DataVectorGetDataSize(item_size) +
                total_size * item_size);
        if ( new_data == NULL )
            return NULL;
        data             = new_data;
        data->total_size = total_size;
        vector           =
            (void*)(((char*)data) + gomp_DataVectorGetDataSize(item_size));
    }

    if ( size > data->size ) {
        /* Init new items. */
        if ( data->handler.InitFunc ) {
            char *array = vector;
            size_t i;
            for ( i=data->size; i<size; ++i )
                data->handler.InitFunc(
                    array + i * item_size, item_size, &data->handler);
        }
        data->size = size;
    }

    return vector;
}
/***************************************************************************/
void *gomp_DataVectorRemoveI(
    void  *vector,
    size_t item_size,
    int    index,
    size_t count)
/***************************************************************************/
{
    DataVectorData *data = GET_DATA(vector, item_size);
    size_t orig_size     = data->size;

    if ( (size_t)index + count > orig_size ) {
        if ( (size_t)index >= orig_size )
            return vector;
        count = orig_size - (size_t)index;
    }

    /* Delete items. Set the size in a way that DataVectorSetSizeI
     * frees count items starting from index.
     */
    data->size = index + count; 
    DataVectorSetSizeI(vector, item_size, (size_t)index, data->total_size);
    /* Set the real size. DataVectorSetSizeI deleted so items. */
    data->size = orig_size - count;
    /* Move the last items to the place of item index. */
    if ( orig_size > (size_t)index + count )
        memmove(
            (char*)vector + item_size * index,
            (char*)vector + item_size * (index + count),
            item_size * (orig_size - (index + count)));

    /* Let gomp_DataVectorSetSizeI decrease the total size using its
     * own size policy. That will always succeed because items are
     * already allocated.
     */
    return gomp_DataVectorSetSizeI(vector, item_size, data->size);
}
/***************************************************************************/
void *gomp_DataVectorSetSizeI(void *vector, size_t item_size, size_t size)
/***************************************************************************/
{
    DataVectorData *data = GET_DATA(vector, item_size);
    size_t total_size;
    void *new_vector;

    if ( data->total_size / 2 <= size && size <= data->total_size )
        /* Array size is quite OK. Keep the old array. */
        total_size = data->total_size;
    else {
        /* Request a new array so that number of used items may vary
         * a little bit without need of reallocation.
         */
        total_size = size + size / 2;
        if ( total_size < 10 )
            total_size = 10;
    }

    new_vector = DataVectorSetSizeI(vector, item_size, size, total_size);
    if ( new_vector )
        /* We got a new array. */
        return new_vector;
    if ( size <= data->total_size )
        /* We can't decrease the size.                */
        /* I don't know if this is actually possible. */
        /* Delete items but keep the old array.       */
        return DataVectorSetSizeI(vector, item_size, size, data->total_size);
    /* Don't try to allocate unnecessary items. */
    return DataVectorSetSizeI(vector, item_size, size, size);
}
/***************************************************************************/
void *gomp_DataVectorSetExactSizeI(void *vector, size_t item_size, size_t size)
/***************************************************************************/
{
    void *new_vector = DataVectorSetSizeI(vector, item_size, size, size);
    if ( ! new_vector ) {
        DataVectorData *data = GET_DATA(vector, item_size);
        if ( size <= data->total_size )
            /* We can't decrease the size.             */
            /* I don't know is this actually possible. */
            /* Delete items but keep the old array.    */
            return DataVectorSetSizeI(
                vector, item_size, size, data->total_size);
    }
    return new_vector;
}
/***************************************************************************/
void *gomp_DataVectorSetTotalSizeI(
    void *vector, size_t item_size, size_t total_size)
/***************************************************************************/
{
    DataVectorData *data = GET_DATA(vector, item_size);
    size_t          size = data->size;
    void *new_vector;
    if ( total_size < size )
        /* We won't destroy elements. */
        total_size = size;
    new_vector = DataVectorSetSizeI(vector, item_size, size, total_size);
    if ( ! new_vector ) {
        if ( data->total_size < total_size )
            /* Reallocation failed. */
            return NULL;
        /* Reallocation failed. But the table size is big enough. */
        return vector;
    }
    return new_vector;
}
/***************************************************************************/
void gomp_DataVectorFreeI(void *vector, size_t item_size)
/***************************************************************************/
{
    if ( vector ) {
        /* Delete items. */
        DataVectorData *data = GET_DATA(vector, item_size);
        DataVectorSetSizeI(vector, item_size, 0, data->total_size);
        /* Free the array. */
        gomp_FreeVector(data);
    }
}

void *gomp_DataVectorTmp;

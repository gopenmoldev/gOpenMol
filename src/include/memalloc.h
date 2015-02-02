#ifndef INC_GOPENMOL_MEMALLOC
#define INC_GOPENMOL_MEMALLOC

/***** PUBLIC GOMAPI BEGIN *****/

/*

                       Copyright (c) 1994 - 2005 by:
        Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
                   Confidential unpublished property of 
                           Leif Laaksonen
                        All rights reserved

        Enhancements 2003 - 2005 by:
            Eero Häkkinen

    Allocation and deallocation of memory
*/

/***** PUBLIC GOMAPI END *****/

#include "gomdefs.h"
#ifdef ENABLE_EXTENSIONS
#  include "gomlib/gommemalloc.h"
#  include "gomext/memalloc.h"
#endif

/***** PUBLIC GOMAPI BEGIN *****/

#include <stddef.h>

/** @weakgroup gom_doc
***/

/** @defgroup gom_doc_memory Memory Management
*** @ingroup  gom_doc
***/

/** @file
*** @ingroup  gom_doc_memory
***/

/** @name Allocate Functions
***/
/* { */
/** @brief Allocate a vector.
***
*** @ingroup  gom_doc_memory
*** This function is similar to malloc().
*** @param count  Number of elements to allocate
*** @return       A newly allocated vector.
*** @retval NULL  Allocation failed.
***/
extern float  *gomp_AllocateFloatVector(
    GOM_ARG( size_t  , count   ) );
extern double *gomp_AllocateDoubleVector(
    GOM_ARG( size_t  , count   ) );
extern int    *gomp_AllocateIntVector(
    GOM_ARG( size_t  , count   ) );
extern char   *gomp_AllocateCharVector(
    GOM_ARG( size_t  , count   ) );
/** @brief Allocate a vector.
*** @ingroup  gom_doc_memory
*** This function is similar to malloc().
*** @param bytes  Number of bytes to allocate
*** @return       A newly allocated vector.
*** @retval NULL  Allocation failed.
***/
extern void   *gomp_AllocateVoidVector(
    GOM_ARG( size_t  , bytes   ) );
/* } */

/** @name Reallocate Functions
***/
/* @{ */
/** @brief Re-allocate a vector.
***
*** @ingroup  gom_doc_memory
***
*** This function is similar to realloc().
***
*** Changes the size of the vector.
*** The elements up to the minimum of the old and the new size
*** will be unchanged.
*** If allocation fails the original vector is left unchanged and must be freed.
*** @param pointer Pointer previously returned
***                by this function or
***                by Allocate function or
***                a NULL pointer
*** @param count   Total number of elements to allocate
*** @return        A newly allocated vector.
*** @retval NULL   Allocation failed.
***/
extern float  *gomp_ReallocateFloatVector(
    GOM_ARG( float  *, pointer ),
    GOM_ARG( size_t  , count   ) );
extern double *gomp_ReallocateDoubleVector(
    GOM_ARG( double *, pointer ),
    GOM_ARG( size_t  , count   ) );
extern int    *gomp_ReallocateIntVector(
    GOM_ARG( int    *, pointer ),
    GOM_ARG( size_t  , count   ) );
extern char   *gomp_ReallocateCharVector(
    GOM_ARG( char    *, pointer ),
    GOM_ARG( size_t  , count   ) );
/** @brief Re-allocate a vector.
***
*** @ingroup  gom_doc_memory
***
*** This function is similar to realloc().
***
*** Changes the size of the vector.
*** The elements up to the minimum of the old and the new size
*** will be unchanged.
*** If allocation fails the original vector is left unchanged and must be freed.
*** @param pointer Pointer previously returned
***                by this function or
***                by Allocate function or
***                a NULL pointer
*** @param bytes   Total number of bytes to allocate
*** @return        A newly allocated vector.
*** @retval NULL   Allocation failed.
***/
extern void   *gomp_ReallocateVoidVector(
    GOM_ARG( void   *, pointer ),
    GOM_ARG( size_t  , bytes   ) );
/* @} */

/** @brief Free a vector.
***
*** @ingroup  gom_doc_memory
*** This function is similar to free().
*** @param pointer Pointer previously returned
***                by gomp_Allocate* function or
***                by gomp_Reallocate* function or
***                a NULL pointer
***/
extern void gomp_FreeVector(
    GOM_ARG( void   *, pointer ) );

/***** PUBLIC GOMAPI END *****/

/**********************************************************/
/* DataVectorHandler contains information needed to       */
/* initialize and destroy vector elements.                */
/* Fields:                                                */
/*     InitFunc:                                          */
/*         Function to initialize a vector element.       */
/*         May be gomp_DataVectorZeroInitFunc or          */
/*                gomp_DataVectorCopyInitFunc             */
/*         The first parametre is address to the element. */
/*     DestroyFunc:                                       */
/*         Function to destroy a vector element.          */
/*         The first parametre is address to the element. */
/*     cdata:                                             */
/*         Source data for gomp_DataVectorCopyInitFunc.   */
/**********************************************************/
typedef struct DataVectorHandler_ {
    void (*InitFunc   )(void*, size_t, const struct DataVectorHandler_ *);
    void (*DestroyFunc)(void*, size_t, const struct DataVectorHandler_ *);
    const void *cdata;
} DataVectorHandler;

/**********************************************************/
/* Internal data structure. Please do not use this.       */
/**********************************************************/
typedef struct DataVectorData_ {
    size_t            size;
    size_t            total_size;
    DataVectorHandler handler;
} DataVectorData;

/**********************************************************/
/* Internal data variable. Please do not use this.        */
/**********************************************************/
extern void *gomp_DataVectorTmp;

/**********************************************************/
/* May be used as an initialization function in           */
/* DataVectorHandler structures. Will initialize vector   */
/* elements with zeros.                                   */
/**********************************************************/
extern void  gomp_DataVectorZeroInitFunc(
    void*, size_t, const DataVectorHandler*);

/**********************************************************/
/* May be used as an initialization function in           */
/* DataVectorHandler structures.                          */
/* Will initialize vector elements with bytes from cdata. */
/**********************************************************/
extern void  gomp_DataVectorCopyInitFunc(
    void*, size_t, const DataVectorHandler*);

/**********************************************************/
/* Allocates a new table for vector. Old table won't be   */
/* freed so vector must be uninitialized.                 */
/* Will create a table of size elements and will          */
/* initialize all of them.                                */
/* Returns the new table on success, NULL on failure.     */ 
/**********************************************************/
extern void *gomp_DataVectorCreateI(
    const DataVectorHandler * /* handler            */,
    size_t                    /* item size          */,
    size_t                    /* size               */);
#define gomp_DataVectorCreate(                \
    /* <item type>**             */ pArray,   \
    /* const DataVectorHandler*  */ pHandler, \
    /* size_t                    */ size)     \
(*(pArray) = gomp_DataVectorCreateI(pHandler,sizeof(**(pArray)),size))

/**********************************************************/
/* If vector is uninitialized does the same as            */
/* gomp_DataVectorCreate does. Otherwise calls the        */
/* gomp_DataVectorSetSize function with increased size.   */
/* Returns the new table on success, NULL on failure.     */ 
/**********************************************************/
extern void *gomp_DataVectorCreateOrAppendI(
    void *                    /* array pointer             */,
    const DataVectorHandler * /* handler                   */,
    size_t                    /* item size                 */,
    size_t                    /* number of items to append */);
#define gomp_DataVectorCreateOrAppend(        \
    /* <item type>**             */ pArray,   \
    /* const DataVectorHandler*  */ pHandler, \
    /* size_t                    */ count)    \
(gomp_DataVectorTmp =                         \
 gomp_DataVectorCreateOrAppendI(*(pArray),pHandler,sizeof(**(pArray)),count), \
 gomp_DataVectorTmp ? (*(pArray)=gomp_DataVectorTmp) : NULL)

/**********************************************************/
/* If vector is uninitialized does the same as            */
/* gomp_DataVectorCreate does. Otherwise calls the        */
/* gomp_DataVectorSetExactSize function with increased    */
/* size.                                                  */
/* Returns the new table on success, NULL on failure.     */ 
/**********************************************************/
extern void *gomp_DataVectorCreateOrAppendExactI(
    void *                    /* array pointer             */,
    const DataVectorHandler * /* handler                   */,
    size_t                    /* item size                 */,
    size_t                    /* number of items to append */);
#define gomp_DataVectorCreateOrAppendExact(   \
    /* <item type>**             */ pArray,   \
    /* const DataVectorHandler*  */ pHandler, \
    /* size_t                    */ count)    \
(gomp_DataVectorTmp =                         \
 gomp_DataVectorCreateOrAppendExactI(         \
     *(pArray),pHandler,sizeof(**(pArray)),count), \
 gomp_DataVectorTmp ? (*(pArray)=gomp_DataVectorTmp) : NULL)

/**********************************************************/
/* Will remove count item from the vector starting from   */
/* the item index. Will destroy unused elements which are */
/* initialized. Vector must be initialized.               */
/* Returns the new table. Will never fail.                */ 
/**********************************************************/
extern void *gomp_DataVectorRemoveI(
    void *                    /* array pointer     */,
    size_t                    /* item size         */,
    int                       /* index             */,
    size_t                    /* count             */);
#define gomp_DataVectorRemove(                \
    /* <item type>**             */ pArray,   \
    /* int                       */ index,    \
    /* size_t                    */ count)    \
(*(pArray) = gomp_DataVectorRemoveI(*(pArray),sizeof(**(pArray)),index,count))

/**********************************************************/
/* Will set number of initialized elements equal to size. */
/* Will reallocate the table if the table is too small.   */
/* Will initialize newly created elements. Will destroy   */
/* unused elements which are initialized. May allocate    */
/* bigger table than absolutely necessary to prevent      */
/* frequent reallocations. Vector must be initialized.    */
/* Returns the new table on success, NULL on failure.     */ 
/* Decreasing the size will never fail.                   */
/**********************************************************/
extern void *gomp_DataVectorSetSizeI(
    void *                    /* array pointer     */,
    size_t                    /* item size         */,
    size_t                    /* size              */);
#define gomp_DataVectorSetSize(  \
    /* <item type>**  */ pArray, \
    /* size_t         */ size)   \
(gomp_DataVectorTmp = gomp_DataVectorSetSizeI( \
    *(pArray),sizeof(**(pArray)),size), \
 gomp_DataVectorTmp ? (*(pArray)=gomp_DataVectorTmp) : NULL)

/**********************************************************/
/* Will create a new table without initializing the new   */
/* elements. This way gomp_DataVectorSetSize function do  */
/* not have to reallocate the table again.                */
/* If total size is smaller than number of initialized    */
/* elements in the vector, then total size is set to be   */
/* equal to that. So gomp_DataVectorSetTotalSize(vector,0)*/
/* will trim the vector. Vector must be initialized.      */
/* Returns the new table on success, NULL on failure.     */ 
/* Decreasing the size will never fail.                   */
/**********************************************************/
extern void *gomp_DataVectorSetTotalSizeI(
    void *                    /* array pointer     */,
    size_t                    /* item size         */,
    size_t                    /* total size        */);
#define gomp_DataVectorSetTotalSize( \
    /* <item type>**  */ pArray,     \
    /* size_t         */ size)       \
(gomp_DataVectorTmp = gomp_DataVectorSetTotalSizeI( \
    *(pArray),sizeof(**(pArray)),size), \
 gomp_DataVectorTmp ? (*(pArray)=gomp_DataVectorTmp) : NULL)

/**********************************************************/
/* Does the same as gomp_DataVectorSetSize except that    */
/* will allocate a new table which size is equal to size. */
/* Vector must be initialized.                            */
/* Returns the new table on success, NULL on failure.     */ 
/* Decreasing the size will never fail.                   */
/**********************************************************/
extern void *gomp_DataVectorSetExactSizeI(
    void *                    /* array pointer     */,
    size_t                    /* item size         */,
    size_t                    /* size              */);
#define gomp_DataVectorSetExactSize( \
    /* <item type>**  */ pArray,     \
    /* size_t         */ size)       \
(gomp_DataVectorTmp = gomp_DataVectorSetExactSizeI( \
    *(pArray),sizeof(**(pArray)),size), \
 gomp_DataVectorTmp ? (*(pArray)=gomp_DataVectorTmp) : NULL)

/**********************************************************/
/* Will free the vector and mark it as uninitialized.     */
/**********************************************************/
extern void        gomp_DataVectorFreeI(
    void *                    /* array pointer     */,
    size_t                    /* item size         */);
#define gomp_DataVectorFree(     \
    /* <item type>**  */ pArray) \
(gomp_DataVectorFreeI(*(pArray),sizeof(**(pArray))),*(pArray) = NULL)

/** Keep a correct allignment. */
#define gomp_DataVectorGetDataSize(item_size) \
    ((sizeof(DataVectorData)+item_size-1)/item_size*item_size)
/**********************************************************/
/* Returns the number of initialized elements in the      */
/* vector.                                                */
/* Returns 0 if vector is uninitialized.                  */
/**********************************************************/
#define gomp_DataVectorGetSize(  \
    /* <item type>**  */ pArray) \
((*(pArray)) ? ((const DataVectorData*)(((const char*)*(pArray))-gomp_DataVectorGetDataSize(sizeof(**(pArray)))))->size : 0)

#endif /* INC_GOPENMOL_MEMALLOC */

/***** PUBLIC GOMAPI BEGIN *****/
#include "gomlistener.h"
/****** PUBLIC GOMAPI END ******/

/* define the basic functions */
#define NEW     0    /* New atom structure will be loaded */
#define APPEND  1    /* Structure will be appended to old list */

#include "gomdefs.h"
#ifdef ENABLE_EXTENSIONS
#  include "gomlib/gommolecstruct.h"
#  include "gomext/molecstruct.h"
#else

/***** PUBLIC GOMAPI BEGIN *****/

/** @weakgroup gom_doc
***/

/** @file
*** @ingroup  gom_doc_struct
*** @ingroup  gom_doc_listener_struct_deleted
***/

/** @defgroup gom_doc_struct Mocule structure
*** @ingroup  gom_doc
*** @{
***/

/** @brief Number of molecule structures.
***/
extern int         gomp_GetNumMolecStructs(void);

/** @brief Molecule structure name.
***/
extern const char *gomp_GetMolecStructName(
    GOM_ARG( int, structure ) );
/** @brief Molecule structure name.
***/
extern       int    gomp_PutMolecStructName(
    GOM_ARG( int         , structure ),
    GOM_ARG( const char *, name      ) );

/** @brief Molecule structure file name.
***/
extern const char *gomp_GetMolecStructFileName(
    GOM_ARG( int, structure ) );
/** @brief Molecule structure file name.
***/
extern       int    gomp_PutMolecStructFileName(
    GOM_ARG( int         , structure ),
    GOM_ARG( const char *, fileName  ) );

/** @}
***/

/** @weakgroup gom_doc_listener
***/

/** @defgroup gom_doc_listener_struct_delete Molecule Structure Delete Listeners
*** @ingroup  gom_doc_listener
*** @ingroup  gom_doc_struct
***
*** Molecule structure delete listeners are called while
*** molecule structure is about to be deleted or merged.
***/

/** @brief Molecule structure pre delete listener handle type.
*** @ingroup gom_doc_listener_struct_delete
***/
GOM_DECLARE_HANDLE(gom_MolecStructPreDeleteListener);

/** @brief Molecule structure delete listener handle type.
*** @ingroup gom_doc_listener_struct_delete
***/
GOM_DECLARE_HANDLE(gom_MolecStructDeleteListener);

/** @brief Molecule structure post delete listener handle type.
*** @ingroup gom_doc_listener_struct_delete
***/
GOM_DECLARE_HANDLE(gom_MolecStructPostDeleteListener);

/** @brief Molecule structure pre delete listener callback type.
*** @ingroup gom_doc_listener_struct_delete
***
*** @par
*** Will be called before a deletion or merge.
*** Deletion or merge will be canceled if the listener returns
*** an error (a positive value).
***
*** @param callbackData  Custom data pointer
*** @param structure     Molecule structure index.
*** @param destination   Structure index of the destinaion structure
***                      in the case of merging or -1.
*** @sa @ref gom_ListenerFunc
***/
typedef int (*gom_MolecStructPreDeleteListenerFunc)(
    GOM_ARG( void *, callbackData ),
    GOM_ARG( int   , structure    ),
    GOM_ARG( int   , destination  ) );

/** @brief Molecule structure delete listener callback type.
*** @ingroup gom_doc_listener_struct_delete
***
*** @par
*** Will be called during a deletion or merge.
*** Original structures still exists during the call.
***
*** @param callbackData  Custom data pointer
*** @param structure     Molecule structure index.
*** @param destination   Structure index of the destinaion structure
***                      in the case of merging or -1.
*** @sa @ref gom_ListenerFunc
***/
typedef int (*gom_MolecStructDeleteListenerFunc)(
    GOM_ARG( void *, callbackData ),
    GOM_ARG( int   , structure    ),
    GOM_ARG( int   , destination  ) );

/** @brief Molecule structure pre delete listener callback type.
*** @ingroup gom_doc_listener_struct_delete
***
*** @par
*** Will be called after a deletion or merge.
***
*** @param callbackData  Custom data pointer
*** @param structure     Molecule structure index.
*** @param destination   Structure index of the destinaion structure
***                      in the case of merging or -1.
*** @sa @ref gom_ListenerFunc
***/
typedef int (*gom_MolecStructPostDeleteListenerFunc)(
    GOM_ARG( void *, callbackData ),
    GOM_ARG( int   , structure    ),
    GOM_ARG( int   , destination  ) );

/** @name Registration Functions
***/
/** @brief Register a new molecule structure deleted listener callback.
*** @ingroup gom_doc_listener_struct_delete
***
*** @param callback      Listener callback function
*** @param callbackData  Pointer to be passed to @p listener
*** @return       New listener handle
*** @retval NULL  An error occured
***/
/* @{ */
extern gom_MolecStructPreDeleteListener* gomp_AddMolecStructPreDeleteListener(
    GOM_ARG( gom_MolecStructPreDeleteListenerFunc, callback     ),
    GOM_ARG( void                               *, callbackData ) );
extern gom_MolecStructDeleteListener* gomp_AddMolecStructDeleteListener(
    GOM_ARG( gom_MolecStructDeleteListenerFunc, callback     ),
    GOM_ARG( void                            *, callbackData ) );
extern gom_MolecStructPostDeleteListener* gomp_AddMolecStructPostDeleteListener(
    GOM_ARG( gom_MolecStructPostDeleteListenerFunc, callback     ),
    GOM_ARG( void                                *, callbackData ) );
/* @} */
/** @brief Unregister a molecule structure deleted listener callback by using the handle.
*** @ingroup gom_doc_listener_struct_delete
***
*** @param listener  Listener handle
*** @return  Number of canceled listeners (always 1).
***/
/* @{ */
extern int gomp_CancelMolecStructPreDeleteListener(
    GOM_ARG( gom_MolecStructPreDeleteListener *, listener ) );
extern int gomp_CancelMolecStructDeleteListener(
    GOM_ARG( gom_MolecStructDeleteListener *, listener ) );
extern int gomp_CancelMolecStructPostDeleteListener(
    GOM_ARG( gom_MolecStructPostDeleteListener *, listener ) );
/* @} */
/** @brief Unregister a molecule structure deleted listener callback by using the callback.
*** @ingroup gom_doc_listener_struct_delete
***
*** @param callback  Listener callback
*** @return  Number of canceled listeners (always non-negative).
***/
/* @{ */
extern int gomp_CancelMolecStructPreDeleteListenersByFunc(
    GOM_ARG( gom_MolecStructPreDeleteListenerFunc, callback ) );
extern int gomp_CancelMolecStructDeleteListenersByFunc(
    GOM_ARG( gom_MolecStructDeleteListenerFunc, callback ) );
extern int gomp_CancelMolecStructPostDeleteListenersByFunc(
    GOM_ARG( gom_MolecStructPostDeleteListenerFunc, callback ) );
/* @} */

/****** PUBLIC GOMAPI END ******/

#endif /* ! ENABLE_EXTENSIONS */

extern int         gomp_HasMolecStructs(void);

extern int         gomp_CreateMolecStruct(const char *, int, int);
extern int         gomp_MergeMolecStruct(int, int);
extern int         gomp_MergeMolecStructs(void);
extern int         gomp_DeleteMolecStruct(int);
extern int         gomp_DeleteMolecStructs(void);

extern int         gomp_PostReadAtoms(int,int);

extern const char *gomp_ShowAtomNameStack(void);
extern const char *gomp_ShowResidueNameStack(void);
extern const char *gomp_ShowSegmentNameStack(void);

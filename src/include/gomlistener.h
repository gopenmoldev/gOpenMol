#ifndef INC_GOPENMOL_GOMLISTENER
#define INC_GOPENMOL_GOMLISTENER

typedef struct gom_Listener_struct gom_Listener;

#include "gomdefs.h"
#ifdef ENABLE_EXTENSIONS
#  include "gomlib/gomlistener.h"
#  include "gomext/gomlistener.h"
#else

/***** PUBLIC GOMAPI BEGIN *****/

/** @weakgroup gom_doc
***/

/** @defgroup gom_doc_listener Listeners
*** @ingroup  gom_doc
***/

/** @file
*** @ingroup  gom_doc_listener
*** @ingroup  gom_doc_listener_impl
*** @ingroup  gom_doc_listener_update_data
***/

/** @defgroup gom_doc_listener_impl Listener Implementation
*** @ingroup  gom_doc_listener
***
*** Listeners are usually implemented either
*** by @ref GOM_IMPLEMENT_LISTENER macro or
*** by @ref GOM_IMPLEMENT_SIMPLE_LISTENER macro so
*** normally there is no need to call these function directly.
***
*** @{
***/

#ifndef GOPENMOL
/** @brief Listener handle type.
***/
GOM_DECLARE_HANDLE(gom_Listener);
#endif

/** @brief Listener list type.
***/
typedef gom_Listener* gom_ListenerList;

/** @brief Generic listener callback type.
***
*** A listener callback function normally have a type
*** @code int CallbackFunction(<args>, void *userData) @endcode
*** The number and the types of the @p \<args\> are listener type specific.
***
*** @retval >0 An error occured
*** @retval 0  Listener callback succeeded
*** @retval <0 Listener callback succeeded and wants to be canceled
***/
typedef int (*gom_ListenerFunc)();

/** @brief Register a new listener callback.
***
*** @param listPointer   Pointer to a listener list
*** @param callback      Listener callback function
*** @param callbackData  Pointer to be passed to @p callback
*** @return       New listener handle
*** @retval NULL  An error occured
***/
extern gom_Listener* gomp_AddListener(
    GOM_ARG( gom_ListenerList *, listPointer  ),
    GOM_ARG( gom_ListenerFunc  , callback     ),
    GOM_ARG( void             *, callbackData ) );

/** @brief Unregister a listener callback using the handle.
***
*** @param listPointer  Pointer to a listener list
*** @param listener     Listener handle
*** @return  Number of canceled listeners (always 1).
***/
extern int gomp_CancelListener(
    GOM_ARG( gom_ListenerList *, listPointer ),
    GOM_ARG( gom_Listener     *, listener    ) );

/** @brief Unregister a listener callback using the callback.
***
*** @param listPointer  Pointer to a listener list
*** @param callback     Listener callback
*** @return  Number of canceled listeners (always non-negative).
***/
extern int gomp_CancelListenersByFunc(
    GOM_ARG( gom_ListenerList *, listPointer ),
    GOM_ARG( gom_ListenerFunc  , callback    ) );

/** @brief Listener helper function type.
***
*** Helper function should do something like this:
*** @code
*** static int ListenerCallHelper(
***     void             *helperData,
***     gom_ListenerFunc  callback,
***     void             *callbackData )
*** {
***     return ((ListenerFuncType)callback)(
***         callbackData,
***         ((HelperDataType*)helperData)->args1,
***         ... );
*** }
*** @endcode
***
*** @retval >0 An error occured
*** @retval 0  Listener callback succeeded
*** @retval <0 Listener callback succeeded and wants to be canceled
***/
typedef int (*gom_ListenerHelperFunc)(
    GOM_ARG( void           *, helperData   ),
    GOM_ARG( gom_ListenerFunc, callback     ),
    GOM_ARG( void           *, callbackData ) );

/** @brief Call listeners.
***
*** Walk through all the listeners in a listener list and pass
*** each of them at a time to a helper function
*** which should call the listener callback with correct arguments.
***
*** @param listPointer Pointer to a listener list
*** @param helper      Helper function
*** @param helperData  Pointer to be passed to @p helper
*** @retval 0  All listeners succeeded (or no listeners).
*** @retval >0 At least one listener failed.
***
*** @sa gomp_CallSimpleListeners
***/
extern int gomp_CallListeners(
    GOM_ARG( gom_ListenerList     *, listPointer ),
    GOM_ARG( gom_ListenerHelperFunc, helper      ),
    GOM_ARG( void                 *, helperData  ) );

/** @brief Simple listener callback type.
***
*** If the listener callback function has only one argument
*** (the callback data pointer),
*** @ref gomp_CallSimpleListeners can be used instead of
*** @ref gomp_CallListeners.
***
*** @sa @ref gom_ListenerFunc
*** @sa gomp_CallListeners
*** @sa gomp_CallSimpleListeners
***/
typedef int (*gom_SimpleListenerFunc)(
    GOM_ARG( void *, callbackData ) );

/** @brief Call simple listeners.
***
*** @param listPointer  Pointer to a listener list
*** @retval 0  All listeners succeeded (or no listeners).
*** @retval >0 At least one listener failed.
***
*** @sa gom_SimpleListenerFunc
*** @sa gomp_CallListeners
***/
extern int gomp_CallSimpleListeners(
    GOM_ARG( gom_ListenerList *, listPointer ) );

#ifdef __cplusplus
    /** @brief Convert a listener callback to a generic callback type.
    ***/
    /* In C++ int(*)() means function with no arguments.
    **/
#   define GOM_GET_LISTENER_CALLBACK(callback) \
        static_cast<gom_ListenerFunc>(callback)
#else
    /** @brief Convert a listener callback to a generic callback type.
    ***/
    /* In C int(*)() means function with no arguments.
    **/
#   define GOM_GET_LISTENER_CALLBACK(callback) (callback)
#endif

/** @brief Get a pointer to the type specific listener list.
***
*** Get a pointer to the type specific listener list implemented
*** by @ref GOM_IMPLEMENT_LISTENER or
*** by @ref GOM_IMPLEMENT_SIMPLE_LISTENER.
***
*** @note Listener list is declared static.
***/
#define GOM_GET_LISTENER_LIST( Type ) &Type##ListenerList

/** @brief Implement type specific listener functions.
***
*** Function
*** @code int <Prefix>Call<Type>Listeners(<args>); @endcode
*** is not implemented. There is no generic implementation for that.
***/
#define GOM_IMPLEMENT_LISTENER( Prefix, Type )   \
                                                           \
    static gom_ListenerList Type##ListenerList = NULL;     \
                                                           \
    Prefix##Type##Listener *Prefix##Add##Type##Listener(   \
        Prefix##Type##ListenerFunc callback,               \
        void                      *callbackData )          \
    {                                                      \
        return (Prefix##Type##Listener *)gomp_AddListener( \
            GOM_GET_LISTENER_LIST( Type ),                 \
            GOM_GET_LISTENER_CALLBACK( callback ),         \
            callbackData );                                \
    }                                                      \
                                                           \
    int Prefix##Cancel##Type##Listener(                    \
        Prefix##Type##Listener *listener )                 \
    {                                                      \
        return gomp_CancelListener(                        \
            GOM_GET_LISTENER_LIST( Type ),                 \
            (gom_Listener *)listener );                    \
    }                                                      \
                                                           \
    int Prefix##Cancel##Type##ListenersByFunc(             \
        Prefix##Type##ListenerFunc callback )              \
    {                                                      \
        return gomp_CancelListenersByFunc(                 \
            GOM_GET_LISTENER_LIST( Type ),                 \
            GOM_GET_LISTENER_CALLBACK( callback ) );       \
    }

/** @brief Implement type specific listener functions for simple listeners.
***/
#define GOM_IMPLEMENT_SIMPLE_LISTENER( Prefix, Type )              \
                                                                   \
    GOM_IMPLEMENT_LISTENER( Prefix, Type )                         \
                                                                   \
    int Prefix##Call##Type##Listeners(void)                        \
    {                                                              \
        return gomp_CallSimpleListeners(                           \
            GOM_GET_LISTENER_LIST( Type ) );                       \
    }

/** @}
***/

/** @defgroup gom_doc_listener_update_data Update Data Listeners
*** @ingroup  gom_doc_listener
***
*** Update data listeners are called to update changed data.
*** They are called by @ref gomp_UpdateData.
***
*** @note This is a main level listener system.
***       Do not create unnecessary update data listeners.
***       Update data listeners should always cancel themselves once
***       they have been called and re-register not until
***       it is absolutely necessary.
***       Use only to create new listener systems.
***
*** Use
*** @ref gomp_AddAtomDataChangedListener,
*** @ref gomp_AddAtomCoordinateChangedListener,
*** @ref gomp_AddMolecStructDeleteListener
*** or similar function to add normal listeners.
***
*** @{
***/

/** @brief Update data listener handle type.
***/
GOM_DECLARE_HANDLE(gom_UpdateDataListener);

/** @brief Update data listener callback type.
*** @param callbackData  Custom data pointer
***/
typedef int (*gom_UpdateDataListenerFunc)(
    GOM_ARG( void *, callbackData ) );

/** @brief Register a new update data listener callback.
***
*** @param callback      Listener callback function
*** @param callbackData  Pointer to be passed to @p listener
*** @return       New listener handle
*** @retval NULL  An error occured
***/
extern gom_UpdateDataListener* gomp_AddUpdateDataListener(
    GOM_ARG( gom_UpdateDataListenerFunc, callback    ),
    GOM_ARG( void                      *, callbackData ) );

/** @brief Unregister an update data listener callback by using the handle.
***
*** @param listener  Listener handle
*** @return  Number of canceled listeners (always 1).
***/
extern int gomp_CancelUpdateDataListener(
    GOM_ARG( gom_UpdateDataListener *, listener ) );

/** @brief Unregister an update data listener callback by using the callback.
***
*** @param callback  Listener callback
*** @return  Number of canceled listeners (always non-negative).
***/
extern int gomp_CancelUpdateDataListenersByFunc(
    GOM_ARG( gom_UpdateDataListenerFunc, callback ) );
    
/** @brief Update all out dated data.
***
*** Calls update data listeners which in turn will call more listeners.
***
*** This function will be called before the display is updated and by
*** a Tcl command @c refresh.
***/
extern int gomp_UpdateData( void );

/** @}
***/
            
/***** PUBLIC GOMAPI END *****/
 
#endif /* ! ENABLE_EXTENSIONS */

/** @brief Implement type specific gOpenMol listener functions.
*** @ingroup gom_doc_listener_impl
***/
#define GOM_IMPLEMENT_GOM_LISTENER(Type) \
    typedef gom_##Type##Listener     gomp_##Type##Listener;     \
    typedef gom_##Type##ListenerFunc gomp_##Type##ListenerFunc; \
    GOM_IMPLEMENT_LISTENER(gomp_,Type)

/** @brief Implement type specific gOpenMol listener functions for simple listeners.
*** @ingroup gom_doc_listener_impl
***/
#define GOM_IMPLEMENT_SIMPLE_GOM_LISTENER(Type) \
    typedef gom_##Type##Listener     gomp_##Type##Listener;     \
    typedef gom_##Type##ListenerFunc gomp_##Type##ListenerFunc; \
    GOM_IMPLEMENT_SIMPLE_LISTENER(gomp_,Type)
    
struct gom_Listener_struct {
    struct gom_Listener_struct *prev;
    struct gom_Listener_struct *next;
    gom_ListenerFunc callback;
    void            *callbackData;
};

#endif /* INC_GOPENMOL_GOMLISTENER */

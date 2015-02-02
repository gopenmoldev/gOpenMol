/***** PUBLIC GOMAPI BEGIN *****/
/*
                       Copyright (c) 2002 - 2003 by:
        Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
                      Confidential unpublished property of 
                              Leif Laaksonen  
                            All rights reserved

        Coded by: Eero Häkkinen
*/

#include "gomlistener.h"
/****** PUBLIC GOMAPI END ******/

typedef struct gom_Plotter_struct gom_Plotter;

#include "gomdefs.h"
#ifdef ENABLE_EXTENSIONS
#  include "gomlib/gomplot.h"
#  include "gomext/plot.h"
#else

/***** PUBLIC GOMAPI BEGIN *****/

/** @weakgroup gom_doc
***/

/** @defgroup gom_doc_plot Plot System 
*** @ingroup  gom_doc
*** @{
***/

/** @file
***/

#ifndef GOPENMOL
/** @brief Plotter handle type.
***/
GOM_DECLARE_HANDLE(gom_Plotter);
#endif

/** @brief Plot flags
***/
enum gom_PlotFlag {
    /** Draw only simple elements (such as lines).
    *** These objects are drown always and they are never
    *** compiled into display lists.
    ***/
    gom_PlotSimpleElements = 0x01,
    /** Draw complex elements in a simple way (as lines).
    *** These objects are drown then user rotates
    *** the molecules using fast drawing mode.
    ***/
    gom_PlotSimplifiedElements = 0x02,
    /** Draw complex elements (not lines).
    *** These objects are drawn then the molecules is
    *** in place and then user rotates the molecules using low drawing mode.
    *** Drawing system may compile these objects into displaylists.
    ***/
    gom_PlotComplexElements = 0x04,
    /** User is rotating the structures.
    *** Heavy fine tunings (such as sorting polygons in z-order)
    *** should be avoided.
    ***/
    gom_PlotRealTimeRotation = 0x10
};

/** @brief Plotter callback type
*** @param callbackData  Custom data pointer
*** @param structure     Structure index
*** @param flags         ORed @ref gom_Plotter flags
*** @retval >0 An error occured.
*** @retval 0  Success.
*** @retval <0 Nothing was drown.
***
*** The purpose of a negative return value is to let the plot system to
*** do some optimization. Plot system has no obligation to follow
*** the return value.
*** Callback function which has returned a negative value
*** will be called with the same parameters at latest after the next call
*** @ref gomp_InvalidatePlotter or
*** @ref gomp_InvalidatePlotters (at latest).
***/
typedef int (*gom_PlotterFunc)(
    GOM_ARG( void *, callbackData ),
    GOM_ARG( int   , structure    ),
    GOM_ARG( int   , flags        ) );

/** @brief Plot order base values.
***/
enum gom_PlotOrderBase {
    gom_PlotOrderBasePlotSequence = 200,
    gom_PlotOrderBaseNormal       = 500,
    gom_PlotOrderBaseSurface      = 500,
    gom_PlotOrderBaseBlend        = 800
};

/** @brief Registers a new plot callback.
***
*** @param callback      Plotter callback function
*** @param callbackData  Pointer to be passed to @p callback
*** @param name          Name (must be a static string)
*** @param order         Plot order
*** @return Handle
***
*** The name is used by the Tcl parser and can include $ sign to
*** mark minimal unique prefix.
*** callback functions will be called in ascending plot order.
***/
extern gom_Plotter* gomp_RegisterPlotter(
    GOM_ARG( gom_PlotterFunc, callback     ),
    GOM_ARG( void          *, callbackData ),
    GOM_ARG( const char    *, name         ),
    GOM_ARG( int            , order        ) );

/** @brief Set plotter specific display list state.
***
*** @param plotter  Plotter handle
*** @param set      Whether to set the state on (!=0) or off (0)
*** @retval 0  Success
***/
extern int gomp_SetPlotterDisplayListState(
    GOM_ARG( gom_Plotter *, plotter ),
    GOM_ARG( int          , set     ) );

/** @brief Invalidate cached plotter specific data.
***
*** The plotter callback will be called to get the new data then needed.
***
*** @param plotter  Plotter handle
*** @retval 0  Success
***
*** @sa @ref gom_PlotterFunc
***/
extern int gomp_InvalidatePlotter(
    GOM_ARG( gom_Plotter *, plotter ) );

/** @brief Invalidate all plotter data.
***
*** Plotter callbacks will be called to get new data then needed.
***
*** @retval 0  Success
***
*** @sa @ref gom_PlotterFunc
***/
extern int gomp_InvalidatePlotters( void );

/** @brief Structure for @ref gomp_InvalidatePlotterDelayed
***/
typedef struct gom_PlotterData_struct {
    gom_Plotter  *plotter;      /**< Plotter handle      */
    gom_Listener *invalidation; /**< Invalidation handle */
} gom_PlotterData;

/** @brief Delayed invalidation of cached plotter specific data.
***
*** @ref gomp_InvalidatePlotter will be called by @ref gomp_UpdateData.
***
*** It is not an error, if @p state->plotter is NULL
*** during the call of this function or 
*** during the call of @ref gomp_UpdateData.
***
*** @param plotterData  Pointer to a plotter data structure.
*** @retval 0  Success
***
*** @sa @ref gomp_InvalidatePlotter
*** @sa @ref gomp_UpdateData
***/
extern int gomp_InvalidatePlotterDelayed(
    GOM_ARG( gom_PlotterData *, plotterData ) );

#ifndef DOXYGEN
#define gom_InvalidatePlotterDelayed( plotterData ) \
    ( \
        ! (plotterData)->invalidation && \
        (gom_InvalidatePlotterDelayed)( plotterData ) )
/****** PUBLIC GOMAPI END ******/
#define gomp_InvalidatePlotterDelayed gom_InvalidatePlotterDelayed
/****** PUBLIC GOMAPI BEGIN ******/
#endif

/** @brief Unregisters a plotter.
***
*** Handle will be invalid after the call.
***
*** @param plotter  Plotter handle
*** @retval 0  Success
***/
extern int gomp_UnregisterPlotter(
    GOM_ARG( gom_Plotter *, plotter ) );

/** @brief Register or unregister a plotter callback.
***
*** @li If @p set is non-zero and @p state->plotter is non-NULL,
***     does nothing and returns 0.
*** @li If @p set is non-zero and @p state->plotter is NULL,
***     calls @ref gomp_RegisterPlotter,
***     stores the result to @p state->plotter and
***     returns 0 on success.
*** @li If @p set is zero and @p state->plotter is non-NULL,
***     calls @ref gomp_UnregisterPlotter,
***     sets @p state->plotter to NULL and
***     returns 0.
*** @li If @p set is zero and @p state->plotter is NULL,
***     does nothing and returns 0.
***
*** @param set           Whether to register (!=0) or unregister (0)
*** @param plotterData   Pointer to a plotter data structure.
*** @param callback      Plotter callback function
*** @param callbackData  Pointer to be passed to @p callback
*** @param name          Name (must be a static string)
*** @param order         Plot order
*** @retval 0  Success
***/
extern int gomp_SetPlotterRegistrationState(
    GOM_ARG( int              , set          ),
    GOM_ARG( gom_PlotterData *, plotterData  ),
    GOM_ARG( gom_PlotterFunc  , callback     ),
    GOM_ARG( void            *, callbackData ),
    GOM_ARG( const char      *, name         ),
    GOM_ARG( int              , order        ) );

/** @brief Draws all structures.
***
*** May use display lists and call plotter callback functions when necessarily.
***
*** @param realTimeRotation  User is rotating the structures
*** @param drawFast          Fast drawing mode while rotating
*** @retval 0  Success
***/
extern int gomp_CallPlotters(
    GOM_ARG( int, realTimeRotation ),
    GOM_ARG( int, drawFast         ) );


/****** PUBLIC GOMAPI END ******/

#endif /* ! ENABLE_EXTENSIONS */

#define PLOTTER_NAME_PLOT_LINE         "plot_l$ine"
#define PLOTTER_NAME_PLOT_SPHERE       "plot_s$phere"
#define PLOTTER_NAME_PLOT_CYLINDER     "plot_c$ylinder"
#define PLOTTER_NAME_PLOT_TRIANGLE     "plot_t$riangle"
#define PLOTTER_NAME_PLOT_ARROW        "plot-a$rrow"
#define PLOTTER_NAME_PLOT_PLANE        "plot_p$lane"
#define PLOTTER_NAME_CUTPLANE          "cutp$lane"
#define PLOTTER_NAME_TRACE             "trac$e"
#define PLOTTER_NAME_MOLECULE          "mole$cule"
#define PLOTTER_NAME_ATOMPICKING       "pick$ing"
#define PLOTTER_NAME_ATOMSELECTION     "sele$ction"
#define PLOTTER_NAME_CELLBOX           "cell"
#define PLOTTER_NAME_COORD_AXIS        "coor$daxis"
#define PLOTTER_NAME_MONITOR_DISTANCE  "monitor_d$istance"
#define PLOTTER_NAME_MONITOR_ANGLE     "monitor_a$ngle"
#define PLOTTER_NAME_MONITOR_TORSION   "monitor_t$orsion"
#define PLOTTER_NAME_LABEL             "labe$l"
#define PLOTTER_NAME_PLUMBER           "plum$ber"
#define PLOTTER_NAME_VECTOR            "vect$or"
#define PLOTTER_NAME_CONTOUR           "cont$our"

#define PLOTTER_ORDER_PLOT_LINE        gom_PlotOrderBasePlotSequence
#define PLOTTER_ORDER_PLOT_SPHERE      gom_PlotOrderBasePlotSequence
#define PLOTTER_ORDER_PLOT_CYLINDER    gom_PlotOrderBasePlotSequence
#define PLOTTER_ORDER_PLOT_TRIANGLE    gom_PlotOrderBasePlotSequence
#define PLOTTER_ORDER_PLOT_ARROW       gom_PlotOrderBasePlotSequence
#define PLOTTER_ORDER_PLOT_PLANE       gom_PlotOrderBasePlotSequence
#define PLOTTER_ORDER_CUTPLANE         (gom_PlotOrderBasePlotSequence+10)
#define PLOTTER_ORDER_TRACE            gom_PlotOrderBaseNormal
#define PLOTTER_ORDER_MOLECULE         gom_PlotOrderBaseNormal
#define PLOTTER_ORDER_ATOMPICKING      gom_PlotOrderBaseNormal
#define PLOTTER_ORDER_ATOMSELECTION    gom_PlotOrderBaseNormal
#define PLOTTER_ORDER_CELLBOX          gom_PlotOrderBaseNormal
#define PLOTTER_ORDER_COORD_AXIS       gom_PlotOrderBaseNormal
#define PLOTTER_ORDER_MONITOR_DISTANCE gom_PlotOrderBaseNormal
#define PLOTTER_ORDER_MONITOR_ANGLE    gom_PlotOrderBaseNormal
#define PLOTTER_ORDER_MONITOR_TORSION  gom_PlotOrderBaseNormal
#define PLOTTER_ORDER_LABEL            gom_PlotOrderBaseNormal
#define PLOTTER_ORDER_PLUMBER          gom_PlotOrderBaseSurface
#define PLOTTER_ORDER_VECTOR           gom_PlotOrderBaseSurface
#define PLOTTER_ORDER_CONTOUR          gom_PlotOrderBaseBlend

extern int gomp_ParseObjectDisplayListState(const char *,int);
extern int gomp_ParseGetObjectDisplayListTypes(void);
extern int gomp_ParseGetObjectDisplayListState(const char *);

extern int gomp_GetDisplayListState(void);
extern int gomp_SetDisplayListState(int);
extern int gomp_GetDefaultDisplayListState(void);
extern int gomp_SetDefaultDisplayListState(int);

extern gom_Plotter **gomp_GetPlotterIterator(void);
extern int gomp_FreePlotterDisplayLists(gom_Plotter*);
extern int gomp_CallPreparePlottersListeners(void);

struct gom_Plotter_struct {
    struct gom_Plotter_struct *next;
    gom_PlotterFunc            callback;
    void                      *callbackData;
    const char                *name;
    int                        order;
    int                        dispListStart;
    int                        dispListCount;
    unsigned short int         bypassSimple;
    unsigned short int         bypassComplex;
};

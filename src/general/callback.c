/*
  Copyright (c) 2002 - 2005 by:
  Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
  Confidential unpublished property of 
  Leif Laaksonen  
  All rights reserved

  Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <tcl.h>

#include "memalloc.h"
#include "molecstruct.h"
#include "plot.h"
#include "tclutils.h"

#include "stdafx.h"

static struct {
    int State;
    int DefaultState;
} DisplayListStates = { 0 , 0 };

/****************************************************************************/
static gom_ListenerList PreparePlottersListenerList;
/****************************************************************************/
int gomp_CallPreparePlottersListeners()
/****************************************************************************/
{
     return gomp_CallSimpleListeners( &PreparePlottersListenerList );
}

/****************************************************************************/
static int InvalidatePlotter( void *callbackData )
/****************************************************************************/
{
    gom_PlotterData *plotterData = callbackData;
    if ( plotterData->plotter )
        gomp_InvalidatePlotter( plotterData->plotter );
    /* Cancel this listener. */
    plotterData->invalidation = NULL;
    return(-1);
}

/****************************************************************************/
#undef gom_InvalidatePlotterDelayed
/****************************************************************************/
int gomp_InvalidatePlotterDelayed(
    gom_PlotterData *plotterData )
/****************************************************************************/
{
    if ( plotterData->invalidation )
        /* Listener is already created. */
        return(0);

    /* Create the listener. */
    plotterData->invalidation = gomp_AddListener(
        &PreparePlottersListenerList,
        InvalidatePlotter,
        plotterData );
    return(plotterData->invalidation ? 0 : 1);
}

/****************************************************************************/
gom_Plotter* gomp_RegisterPlotter(
    gom_PlotterFunc callback,
    void           *callbackData,
    const char     *name,
    int             order )
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    gom_Plotter **stackPointer;
    gom_Plotter  *newPlotter;

    /* Allocate space for the new function. */
    newPlotter = malloc(sizeof(*newPlotter));

    if ( ! newPlotter )
        return NULL;

    /* Find the right place for the new callback function. */
    stackPointer = gomp_GetPlotterIterator();
    while ( *stackPointer && (*stackPointer)->order < order )
        stackPointer = &(*stackPointer)->next;

    /* Store the function. */
    newPlotter->next          = *stackPointer;
    *stackPointer             = newPlotter;

    newPlotter->callback      = callback;
    newPlotter->callbackData  = callbackData;
    newPlotter->name          = name;
    newPlotter->order         = order;
    newPlotter->dispListStart = 0;
    newPlotter->dispListCount = -1;
    newPlotter->bypassSimple  = 0;
    newPlotter->bypassComplex = 0;
    
    gomp_SetPlotterDisplayListState(
        newPlotter, gomp_GetDefaultDisplayListState() );

    return newPlotter;
#else
    return NULL;
#endif /* ENABLE_GRAPHICS */
}

/****************************************************************************/
int gomp_SetPlotterDisplayListState(
    gom_Plotter *plotter, int set )
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    int NStruct;

    if ( ! set ) {
        /* Free display lists. */
        gomp_FreePlotterDisplayLists( plotter );
        plotter->dispListCount = -1;
    }
    else if ( plotter->dispListCount < 0 ) {
        /* We don't have display lists yet. */
        NStruct = gomp_GetNumMolecStructs();
        if ( NStruct < 1 )
            NStruct = 1;
        plotter->dispListCount = NStruct;
    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}

/****************************************************************************/
int gomp_InvalidatePlotter( gom_Plotter *plotter )
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    if ( ! plotter )
        return(0);

    /* Try to draw all structures. */
    plotter->bypassSimple  = 0;
    plotter->bypassComplex = 0;

    /* Reset display lists. */
    gomp_FreePlotterDisplayLists( plotter );
    
    if ( plotter->dispListCount >= 0 ) {
        int NStruct = gomp_GetNumMolecStructs();
        if ( NStruct < 1 )
            NStruct = 1;
        plotter->dispListCount = NStruct;
    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}

/****************************************************************************/
int gomp_InvalidatePlotters()
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    gom_Plotter *iter;

    for ( iter = *gomp_GetPlotterIterator() ; iter ; iter = iter->next )
        gomp_InvalidatePlotter( iter );
#endif /* ENABLE_GRAPHICS */

    return(0);
}

/****************************************************************************/
int gomp_UnregisterPlotter( gom_Plotter *plotter )
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    gom_Plotter *iter,**pIter;

    if ( ! plotter )
        return(0);

    pIter = gomp_GetPlotterIterator();

    for ( iter = *pIter ; iter ; pIter = &iter->next , iter = *pIter ) {
        if ( iter == plotter ) {
            gomp_FreePlotterDisplayLists( iter );
            *pIter = iter->next;
            free(iter);
            return(0);
        }
    }
#endif /* ENABLE_GRAPHICS */

    return(1);
}
/****************************************************************************/
int gomp_SetPlotterRegistrationState(
    int              set,
    gom_PlotterData *plotterData,
    gom_PlotterFunc  callback,
    void            *callbackData,
    const char      *name,
    int              order )
/****************************************************************************/
{
    if ( set ) {
        if ( plotterData->plotter )
            /* Drawing callback is already registered. */
            return(0);
        /* Register the drawing callback. */
        plotterData->plotter = gomp_RegisterPlotter(
            callback, callbackData, name, order );
        return(plotterData->plotter ? 0 : 1);
    }
    else {
        if ( ! plotterData->plotter )
            /* Drawing callback is not registered. */
            return(0);
        /* Unregister the drawing callback. */
        gomp_UnregisterPlotter( plotterData->plotter );
        plotterData->plotter = NULL;
        return(0);
    }
}
/****************************************************************************/
int gomp_ParseObjectDisplayListState( const char *type, int set )
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    gom_Plotter *iter;
    int count;

    count = 0;

    for ( iter = *gomp_GetPlotterIterator() ; iter ; iter = iter->next ) {
        if ( gomp_StringMatch(type, iter->name) ) {
            gomp_SetPlotterDisplayListState(iter, set );
            count++;
        }
    }

    if ( count )
        return(0);
#endif /* ENABLE_GRAPHICS */

    /* unrecognized type */
    return(1);
}

/****************************************************************************/
int gomp_ParseGetObjectDisplayListTypes()
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    gom_Plotter *iter,*iter2;
    char *list,*p;
    const char *q;
    size_t length = 0;

    for ( iter = *gomp_GetPlotterIterator() ; iter ; iter = iter->next ) {
        /* Count each type only once. */
        for( iter2 = iter->next ; iter2 ; iter2 = iter2->next ) {
            if ( strcmp(iter->name, iter2->name) == 0 )
                break;
        }
        if ( iter2 )
            continue;

        /* Omit $-signs. */
        for ( q = iter->name ; *q ; q++ ) {
            if ( *q == '$' )
                continue;
            ++length;
        }
        /* Space for deliminators. */
        length += 3;
    }

    if ( length == 0 ) {
#endif /* ENABLE_GRAPHICS */
        (void)gomp_SendTclReturn("{}");
        return(0);
#ifdef ENABLE_GRAPHICS
    }

    list = gomp_AllocateCharVector(length);
    p    = list;

    for ( iter = *gomp_GetPlotterIterator() ; iter; iter = iter->next ) {
        /* Copy each type only once. */
        for( iter2 = iter->next ; iter2 ; iter2 = iter2->next ) {
            if ( strcmp(iter->name, iter2->name) == 0 )
                break;
        }
        if ( iter2 )
            continue;

        /* Omit $-signs. */
        *p++ = '{';
        for ( q = iter->name ; *q ; q++ ) {
            if ( *q == '$' )
                continue;
            *p++ = *q;
        }
        /* Add deliminator. */
        *p++ = '}';
        *p++ = ' ';
    }
    *--p = '\0';
    (void)gomp_SendTclReturn(list);
    free(list);

    return(0);
#endif /* ENABLE_GRAPHICS */
}

/****************************************************************************/
int gomp_ParseGetObjectDisplayListState( const char *type )
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    gom_Plotter *iter;

    for ( iter = *gomp_GetPlotterIterator() ; iter ; iter = iter->next ) {
        if ( gomp_StringMatch(type, iter->name) ) {
            if ( iter->dispListCount >= 0 )
                (void)gomp_SendTclReturn("1");
            else
                (void)gomp_SendTclReturn("0");
            return(0);
        }
    }
#endif /* ENABLE_GRAPHICS */

    /* unrecognized type */
    return(1);
}

/****************************************************************************/
int gomp_GetDisplayListState()
/****************************************************************************/
{
    return(DisplayListStates.State);
}

/****************************************************************************/
int gomp_SetDisplayListState( int set )
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    if ( ! set )
        gomp_InvalidatePlotters();
#endif /* ENABLE_GRAPHICS */

    DisplayListStates.State = set;

    return(0);
}

/****************************************************************************/
int gomp_GetDefaultDisplayListState()
/****************************************************************************/
{
    return(DisplayListStates.DefaultState);
}

/****************************************************************************/
int gomp_SetDefaultDisplayListState( int set )
/****************************************************************************/
{
    DisplayListStates.DefaultState = set;

    return(0);
}

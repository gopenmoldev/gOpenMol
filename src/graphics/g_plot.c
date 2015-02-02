/*
  Copyright (c) 2002 - 2005 by:
  Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
  Confidential unpublished property of 
  Leif Laaksonen  
  All rights reserved

  Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#ifdef ENABLE_GRAPHICS
#include "gomstdio.h"
#include <string.h>
#include <ctype.h>
#include <math.h>

#if defined(WIN32)
#include <windows.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>

#include "molecstruct.h"
#include "plot.h"
#include "projview.h"
#endif /* ENABLE_GRAPHICS */

#include "stdafx.h"

#ifdef ENABLE_GRAPHICS

static int InitPlotterDisplayLists( gom_Plotter * );

static gom_Plotter *PlotterStart;

/****************************************************************************/
gom_Plotter **gomp_GetPlotterIterator(void)
/****************************************************************************/
{
    return &PlotterStart;
}

/****************************************************************************/
int gomp_FreePlotterDisplayLists( gom_Plotter* plotter )
/****************************************************************************/
{
    if ( plotter->dispListStart > 0 ) {
        glDeleteLists( plotter->dispListStart, plotter->dispListCount );
        plotter->dispListStart = 0;
    }

    return(0);
}

/****************************************************************************/
int InitPlotterDisplayLists( gom_Plotter *plotter )
/****************************************************************************/
{
    const float *RotMP;
    unsigned short int bypass;
    int Wstr,value;

    /* Create type specific display lists. */
    if ( plotter->dispListStart > 0 )
        /* Display lists are already created. */
        return(0);

    if ( plotter->dispListCount <= 0 )
        /* User don't want to use display list for this type or */
        /* there is no need for display lists.                  */
        return(0);
        
    /* Try to create new display lists. */
    plotter->dispListStart = glGenLists( plotter->dispListCount );
    if( plotter->dispListStart <= 0 )
        return(1);

    /* We got new display lists.                 */
    /* Let's compile complex elements into them. */
    for ( Wstr = plotter->dispListCount-1 ; Wstr >= 0 ; Wstr-- ) {
        bypass = (unsigned short int)1<<Wstr;

        /* Call back function may want to retrieve */
        /* matrix and viewport data using glGet.   */
        /* To enable that we apply current view    */
        /* structure specific matrix.              */
        RotMP = gomp_GetSavedModelViewMatrixMT( Wstr );
        glPushMatrix();
        glMultMatrixf(RotMP);
        glNewList( plotter->dispListStart + Wstr, GL_COMPILE );
        if ( plotter->bypassComplex & bypass )
            value = -1;
        else
            value = plotter->callback(
                plotter->callbackData, Wstr, gom_PlotComplexElements );
        glEndList();
        glPopMatrix();

        if ( value < 0 ) {
            plotter->bypassComplex |= bypass;
            /* Test if this is the last display list so that */
            /* we can easily delete it.                      */
            if ( Wstr + 1 == plotter->dispListCount ) {
                glDeleteLists( plotter->dispListStart + Wstr, 1 );
                --plotter->dispListCount;
            }
        }
    }
    if ( plotter->dispListCount == 0 )
        plotter->dispListStart = 0;

    return(0);
}

/****************************************************************************/
int gomp_CallPlotters( int rotating, int drawFast )
/****************************************************************************/
{
    gom_Plotter* iter;
    unsigned short int   bypass;
    int                  Wstr, NStruct,DrawingFlags;
    const float *RotMP;

    /* Invalidate invalid display lists. */
    gomp_CallPreparePlottersListeners();

    NStruct = gomp_GetNumMolecStructs();

    if ( rotating && drawFast ) {
        /* Draw only simple elements. */
        DrawingFlags =
            gom_PlotSimpleElements |
            gom_PlotSimplifiedElements |
            gom_PlotRealTimeRotation;
        for ( iter = *gomp_GetPlotterIterator() ; iter ; iter = iter->next ) {
            for ( Wstr = 0, bypass = 1
                      ; Wstr < NStruct ; Wstr++, bypass <<= 1 ) {
                if ( (iter->bypassSimple  & bypass) &&
                     (iter->bypassComplex & bypass) )
                    continue;
                RotMP = gomp_GetSavedModelViewMatrixMT( Wstr );
                glPushMatrix();
                glMultMatrixf(RotMP);
                if ( iter->callback(
                    iter->callbackData, Wstr, DrawingFlags ) < 0 )
                    iter->bypassSimple |= bypass;
                glPopMatrix();
            }
        }

        return(0);
    }

    if ( gomp_GetDisplayListState() ) {
        DrawingFlags =
            gom_PlotSimpleElements |
            ( rotating ? gom_PlotRealTimeRotation : 0 );
        for ( iter = *gomp_GetPlotterIterator() ; iter ; iter = iter->next ) {
            if ( iter->dispListCount >= 0 &&
                 ( iter->dispListStart > 0 ||
                   InitPlotterDisplayLists(iter) == 0 ) ) {
                /* We have a display list.               */
                /* Let's draw complex elements using it. */
                for ( Wstr = 0, bypass = 1 ;
                      Wstr < NStruct ; Wstr++, bypass <<= 1 ) {
                    RotMP = gomp_GetSavedModelViewMatrixMT(Wstr);
                    glPushMatrix();
                    glMultMatrixf(RotMP);
                    if ( Wstr < iter->dispListCount )
                        glCallList( iter->dispListStart + Wstr );
                    if ( iter->bypassSimple  & bypass )
                        ;
                    else if ( iter->callback(
                        iter->callbackData, Wstr, DrawingFlags ) < 0 )
                        iter->bypassSimple |= bypass;
                    glPopMatrix();
                }
            }
            else {
                /* Draw both the simple and complex elements. */
                for ( Wstr = 0, bypass = 1 ;
                      Wstr < NStruct ; Wstr++, bypass <<= 1 ) {
                    RotMP = gomp_GetSavedModelViewMatrixMT(Wstr);
                    glPushMatrix();
                    glMultMatrixf(RotMP);
                    if ( ( iter->bypassSimple  & bypass ) &&
                         ( iter->bypassComplex & bypass ) )
                        ;
                    else if ( iter->callback(
                        iter->callbackData,Wstr,
                        DrawingFlags | gom_PlotComplexElements) < 0 ) {
                        iter->bypassSimple  |= bypass;
                        iter->bypassComplex |= bypass;
                    }
                    glPopMatrix();
                }
            }
        }

        return(0);
    }
    
    /* Display lists are not in use. */
    DrawingFlags =
        gom_PlotSimpleElements |
        gom_PlotComplexElements |
        ( rotating ? gom_PlotRealTimeRotation : 0 );
    for ( iter = *gomp_GetPlotterIterator() ; iter ; iter = iter->next ) {
        /* Draw both the simple and complex elements. */
        for ( Wstr = 0, bypass = 1 ;
              Wstr < NStruct ; Wstr++, bypass <<= 1 ) {
            RotMP = gomp_GetSavedModelViewMatrixMT(Wstr);
            glPushMatrix();
            glMultMatrixf(RotMP);
            if ( ( iter->bypassSimple  & bypass ) &&
                 ( iter->bypassComplex & bypass ) )
                ;
            else if ( iter->callback(
            iter->callbackData, Wstr, DrawingFlags ) < 0 ) {
                iter->bypassSimple  |= bypass;
                iter->bypassComplex |= bypass;
            }
            glPopMatrix();
        }
    }

    return(0);
}
#else /* ENABLE_GRAPHICS */
extern int i;
#endif

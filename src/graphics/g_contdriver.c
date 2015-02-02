/*

Copyright (c) 1995 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <string.h>
#include <sys/types.h>
#include <stdlib.h>

#include "contour.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"

#include "stdafx.h"

#define RED_COLOUR    0
#define GREEN_COLOUR  1
#define BLUE_COLOUR   2

#define CONTOUR_LEVELS_MAX   10

#define   CONTDRIVER "contour_widget.hlp"

/* ............................................................ */


static int  Contour2StructureMappingEntries = 0;
static int *Contour2StructureMappingValues;

/*  structure to handle the display on/off toggle button status */


static struct {
    int   Which;
    int   Type;
    int   Method; /* == 0 default (direct), != 0 save contour data */
    unsigned long int change_mask;
    gom_PlotterData Plotter;
} DispSurface = { 0 , 0 , 0 , 0 , { NULL, NULL } };

#if 0
static struct {
    int   WidgetActive;    /* keeps track of the widget activation */
    float red;             /* current red colour   */
    float green;           /* current green colour */
    float blue;            /* current blue colour  */
    float alpha;
} ContourColourBrowser = { 0 , 0.0 , 0.0 , 0.0 , 1.0};
#endif

#if 0
static int AccessingContourLevel = 0;
static int GetAccessingContourLevel(void);
#endif

/****************************************************************************/
int gomp_GetSurfaceMethod()
/****************************************************************************/
{
    return(DispSurface.Method);
}
/****************************************************************************/
int gomp_SetSurfaceMethod(int Method)
/****************************************************************************/
{
    gomp_ContourIsChanging(-1);

    DispSurface.Method = Method;

    return(0);
}
/****************************************************************************/
int gomp_DeleteInvalidatedPolygonData()
/****************************************************************************/
{
    if ( DispSurface.change_mask ) {
        /* delete reserved spaces */
        int Contour, Level;
        unsigned long int bit = 1;
        for ( Contour = 0 ;
              Contour < gomp_GetContoursDefined() ; ++Contour, bit <<= 1 ) {
            if ( bit && ! ( DispSurface.change_mask & bit ) )
                /* This contour haven't changed. */
                continue;
            for ( Level = 0 ;
                  Level < gomp_GetContourLevels(Contour) ; ++Level )
                gomp_DataVectorFree(
                    &ContourInfo[Contour].levels[Level].polygons);
        }
        DispSurface.change_mask = 0;
    }

    return(0);
}
#if 0
/****************************************************************************/
int SurfaceControl()
/****************************************************************************/
{
    return(DispSurface.Plotter.Plotter!=NULL);
}
#endif
/****************************************************************************/
int gomp_InvalidateContourPlotter(unsigned long int change_mask)
/****************************************************************************/
{
    /* Invalidate saved polygon data. But don't delete it before */
    /* the the next drawing. That data might be needed by        */
    /* "show polygon" command.                                   */

    /* change_mask may have been overflowed. */
    DispSurface.change_mask = change_mask ? change_mask : 0x1;

    gomp_InvalidatePlotterDelayed(&DispSurface.Plotter);

    return(0);
}
/****************************************************************************/
int gomp_SetSurfaceControlON()
/****************************************************************************/
{
    gomp_SetPlotterRegistrationState(
        1, &DispSurface.Plotter, gomp_PlotIsoSurf, NULL,
        PLOTTER_NAME_CONTOUR,PLOTTER_ORDER_CONTOUR);

    return(0);
}
/****************************************************************************/
int gomp_SetSurfaceControlOFF()
/****************************************************************************/
{
    gomp_SetPlotterRegistrationState(
        0, &DispSurface.Plotter, gomp_PlotIsoSurf, NULL,
        PLOTTER_NAME_CONTOUR,PLOTTER_ORDER_CONTOUR);

    return(0);
}
/****************************************************************************/
int gomp_DeleteAllContours()
/****************************************************************************/
{
/* the order here is important (it looks uggly I know) */
    gomp_ContourIsChanging(-1);
    DispSurface.Which  = 0;
    DispSurface.Type   = 0;
    gomp_SetSurfaceControlOFF();

/* delete reserved contour memory */
    gomp_DataVectorFree(&ContourInfo);

/* delete the cut plane information as well */
    (void)gomp_DeleteCutPlaneDataX();
    (void)gomp_DeleteCutPlaneDataY();
    (void)gomp_DeleteCutPlaneDataZ();
    
/* delete contour to structure mapping */
    (void)gomp_DeleteGetContour2StructureMapping();

/* delete any possible XYZ cutplanes */
    (void)gomp_DisableCutPlanePlotStateXYZ(1);
    (void)gomp_DisableCutPlanePlotStateXYZ(2);
    (void)gomp_DisableCutPlanePlotStateXYZ(3);

    (void)gomp_DeleteGetContour2StructureMapping();

    return(0);
}
#if 0
/****************************************************************************/
int GetAccessingContourLevel()
/****************************************************************************/
{
    return(AccessingContourLevel);
}
/****************************************************************************/
int StartupClean()
/****************************************************************************/
{
    printf("Dummy in gomp_StartupClean\n");

    return(0);
}
#endif
/****************************************************************************/
int  gomp_AddContour2StructureMappingSpace()
/****************************************************************************/
{
    if(Contour2StructureMappingEntries) {

        Contour2StructureMappingValues = 
            gomp_ReallocateIntVector(Contour2StructureMappingValues,
                                 Contour2StructureMappingEntries + 1);
        Contour2StructureMappingEntries++;

    } else {

        Contour2StructureMappingValues = 
            gomp_AllocateIntVector(Contour2StructureMappingEntries + 1);
        Contour2StructureMappingEntries = 1;

    }


    return(0);
}
/****************************************************************************/
int  gomp_GetContour2StructureMapping(int ContourNr)
/****************************************************************************/
{
    if((ContourNr < 0) || (ContourNr >= Contour2StructureMappingEntries)) {
        gomp_PrintERROR("contour to structure [contour] out of allowed range");
        return(-1);
    }

    return(Contour2StructureMappingValues[ContourNr]);
}
/****************************************************************************/
int  gomp_DeleteGetContour2StructureMapping()
/****************************************************************************/
{
    if(Contour2StructureMappingEntries) {

        free(Contour2StructureMappingValues);
        Contour2StructureMappingEntries = 0;
    }

    return(0);
}
/****************************************************************************/
int  gomp_SetContour2StructureMapping(int ContourNr, int StructureNr)
/****************************************************************************/
{
    if((ContourNr < 0) || (ContourNr >= Contour2StructureMappingEntries)) {
        gomp_PrintERROR("contour to structure [contour] out of allowed range [error 1]");
        return(1);
    }
    if((StructureNr < 0) || (StructureNr >= gomp_GetNumMolecStructs())) {
        gomp_PrintERROR("contour to structure [structure] out of allowed range [error 2]");
        return(1);
    }

    gomp_ContourIsChanging(ContourNr);

    Contour2StructureMappingValues[ContourNr] = StructureNr;

    return(0);
}

/*

Copyright (c) 1990 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include "gomstdio.h"
#include <math.h>
#include <string.h>
#include <tcl.h>

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>
#endif /* ENABLE_GRAPHICS */

#if defined(MOTIF)
#include <GL/GLwMDrawA.h>
#endif

#include "colouring.h"
#include "gomstring.h"
#include "gomtext.h"
#include "label.h"
#include "molecoord.h"
#include "molecule.h"
#include "plot.h"
#include "selection.h"
#include "tclutils.h"

#include "stdafx.h"

static struct {
    gom_PlotterData Plotter;
    int             Type;
} Label = { { NULL , NULL } , 0 };

static int  PlotLabel(void*,int,int);
#define AtomLabelDataIsChanging() \
    gomp_InvalidatePlotterDelayed(&Label.Plotter)


/***********************************************************************/
int  PlotLabel(void* userData,int Wstr,int drawFlags)
/***********************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int    i;
    static const float *x;
    static const float *y;
    static const float *z;
    static const char *label;
    static char   text[BUFF_LEN];
    static char   font[BUFF_LEN];
    static const float *red;
    static const float *green;
    static const float *blue;
    static const char *Value;
    static float CRi;
    static float CGi;
    static float CBi;
/*  static int   ITemp;*/
    static int   INumLabels;
    static int   LabelSizeLimit;

    if ( ! ( drawFlags & gom_PlotComplexElements ) )
        return(-1);

    if(!gomp_GetSelectedStructure(Wstr))
        return(0);

    Value = Tcl_GetVar(gomp_GetTclInterp(),
                       "gomAtomLabelFont",TCL_GLOBAL_ONLY);

    if(Value) {
        gomp_CopyString(font,Value,BUFF_LEN);
    } else {
        gomp_CopyString(font,"BITMAP_HELVETICA_12",BUFF_LEN);
    }

    Value = Tcl_GetVar(gomp_GetTclInterp(),
                       "gomAtomLabelNumberLimit",TCL_GLOBAL_ONLY);

    if(Value) {
        LabelSizeLimit = atoi(Value);
    } else {
        LabelSizeLimit = 200;
    }

    glDisable(GL_LIGHTING);

    label     = gomp_GetAtomLabelDisplayStatePointer(Wstr);
    x         = gomp_GetAtomXCoordPointer(Wstr);
    y         = gomp_GetAtomYCoordPointer(Wstr);
    z         = gomp_GetAtomZCoordPointer(Wstr);

    red       = gomp_GetAtomColourRedPointer(Wstr);
    green     = gomp_GetAtomColourGreenPointer(Wstr);
    blue      = gomp_GetAtomColourBluePointer(Wstr);

/* total number of labels */
    INumLabels = 0;
    for(i = 0   ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
        if((int)label[i] < 1) continue;
        INumLabels++;
    }

    for(i = 0   ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if((int)label[i] < 1) continue;

        CRi = red[i];
        CGi = green[i];
        CBi = blue[i];

        if(!gomp_GetDisplayColourType())
            (void)gomp_RGB2Grayscale(&CRi , &CGi , &CBi);

        glColor3f(CRi , CGi , CBi);
        glRasterPos3f(x[i] , y[i] , z[i]);

/*        if( INumLabels <= LabelSizeLimit) {

sprintf(text,"lulPrepareAtomLabel %d %d %d", i+1 , j+1 , gomp_GetPlotLabelType());
ITemp = Tcl_GlobalEval(gomp_GetTclInterp(),text);
if(ITemp != TCL_OK) {
gomp_PrintERROR("can't execute the script 'lulPrepareAtomLabel'");
return(1);
}
} else {
*/
        switch(gomp_GetPlotLabelType()) {
        case FULL_LABEL_TYPE:
            sprintf(text,"%s:%s(%d):%s(%d)",
                    gomp_GetAtomSegName(Wstr,i),
                    gomp_GetAtomResName(Wstr,i),gomp_GetAtomResNum1(Wstr,i),
                    gomp_GetAtomAtmName(Wstr,i),(i+1));
            break;
        case ATOM_LABEL_TYPE:
            sprintf(text,"%s(%d)",gomp_GetAtomAtmName(Wstr,i),(i+1));
            break;
        case RESIDUE_LABEL_TYPE:
            sprintf(text,"%s(%d)",gomp_GetAtomResName(Wstr,i),gomp_GetAtomResNum1(Wstr,i));
            break;
        }

        /*      }*/

        gomp_PrintString(text , font);
    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}
/***********************************************************************/
int  gomp_AtomLabelDataIsChanging()
/***********************************************************************/
{
    AtomLabelDataIsChanging();

    return(0);
}
/***********************************************************************/
int  gomp_SetPlotLabelState(int State)
/***********************************************************************/
{
    AtomLabelDataIsChanging();
    
    return gomp_SetPlotterRegistrationState(
        State, &Label.Plotter, PlotLabel, NULL,
        PLOTTER_NAME_LABEL, PLOTTER_ORDER_LABEL);
}

/***********************************************************************/
int gomp_GetPlotLabelState()
/***********************************************************************/
{
    return(Label.Plotter.plotter!=NULL);
}

/***********************************************************************/
int  gomp_SetPlotLabelType(int Value)
/***********************************************************************/
{
    AtomLabelDataIsChanging();

    Label.Type = Value;

    return(0);
}

/***********************************************************************/
int gomp_GetPlotLabelType()
/***********************************************************************/
{
    return(Label.Type);
}

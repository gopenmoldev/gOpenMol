/*

Copyright (c) 1995-2004 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

in collaboration with

OpenMol Molecular Astrophysics Group
Max-Planck-Institut  fuer Astrophysik
Garching, GERMANY

Enhancements 2002 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <string.h>

#include "trace.h"

#include "stdafx.h"

#define TRACE_ON  1
#define TRACE_OFF 0

#define   TRACE_DRIVER "trace_driver_widget.hlp"

/*  structure to handle the display on/off toggle button status */
static struct {
    int MatrixDisplay;
} TraceSelectionMask = { 0 };

/****************************************************************************/
int  gomp_GetDisplayTraceAtoms()
/****************************************************************************/
{
    return(TraceSelectionMask.MatrixDisplay);
}

/****************************************************************************/
int  gomp_SetDisplayTraceAtoms(int LdpStatus)
/****************************************************************************/
{
    TraceSelectionMask.MatrixDisplay = LdpStatus;

    return(0);
}

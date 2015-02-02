/*

Copyright (c) 1995 - 2004 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>

#include "ldp.h"

#include "stdafx.h"

/*  structure to handle the display on/off toggle button status */
static struct {
    int WidgetDisplay;
    int MatrixDisplay;
} LdpSelectionMask = { 1 , 0 };

#define ATOM1 1
#define ATOM2 2

/****************************************************************************/
int  gomp_GetDisplayLDPmatrix()
/****************************************************************************/
{
    return(LdpSelectionMask.MatrixDisplay);
}

/****************************************************************************/
int  gomp_SetDisplayLDPmatrix(int LdpStatus)
/****************************************************************************/
{
    LdpSelectionMask.MatrixDisplay = LdpStatus;

    return(0);
}
/****************************************************************************/
int gomp_ApplyLdpSelectionMaskLine(const char *GetTextSeg ,
                                 const char *GetTextRes ,
                                 const char *GetTextAtm)
/****************************************************************************/
{

    if(gomp_PushAtomToLDP(ATOM1 , GetTextSeg , GetTextRes , GetTextAtm))
        return(1);

    if(gomp_PushAtomToLDP(ATOM2 , GetTextSeg , GetTextRes , GetTextAtm))
        return(1);

/*    if(gomp_BuildLDParray()) return(1); */

    return(0);
}

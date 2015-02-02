/*

Copyright (c) 1995 - 2004 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Enhancements 2002 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <string.h>
#include <stdlib.h>

#if defined(WIN32)
#include <windows.h>
#endif

#include "contour.h"
#include "plot_molec.h"

#include "stdafx.h"

/* this section contains various parameters related to the display
   of various things */

static struct {
    int MoleculeWidth;
    int ContourWidth;
} LineDrawing = { 1 , 1 };

/****************************************************************************/
int gomp_SetMoleculeLineWidth(int Width)
/****************************************************************************/
{
    LineDrawing.MoleculeWidth = Width;

    return(0);
}
/****************************************************************************/
int gomp_GetMoleculeLineWidth()
/****************************************************************************/
{
    return(LineDrawing.MoleculeWidth);
}
/****************************************************************************/
int gomp_SetContourLineWidth(int Width)
/****************************************************************************/
{
    LineDrawing.ContourWidth = Width;

    return(0);
}
/****************************************************************************/
int gomp_GetContourLineWidth()
/****************************************************************************/
{
    return(LineDrawing.ContourWidth);
}

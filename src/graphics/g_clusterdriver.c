/*

Copyright (c) 1995 - 2004 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

*/

#include "maindefs.h"

#include "gomstdio.h"

#include "cluster.h"

#include "stdafx.h"

#define CLUSTER_ON  1
#define CLUSTER_OFF 0

int  DisplayClusterSelectionMaskWidget(int);

#if 0
static int ClusterSelectionMaskWidgetState(void);
static int ClusterSelectionMaskWidgetActive = 0;
#endif

/*  structure to handle the display on/off toggle button status */
static struct {
    int WidgetDisplay;
    int MatrixDisplay;
} ClusterSelectionMask = { 1 , 0 };


#if 0
/****************************************************************************/
int ClusterSelectionMaskWidgetState(void)
/****************************************************************************/
{
    return(ClusterSelectionMaskWidgetActive);
}
#endif
/****************************************************************************/
int  gomp_GetDisplayCLUSTERmatrix()
/****************************************************************************/
{
    return(ClusterSelectionMask.MatrixDisplay);
}

/****************************************************************************/
int  gomp_SetDisplayCLUSTERmatrix(int ClusterStatus)
/****************************************************************************/
{
    ClusterSelectionMask.MatrixDisplay = ClusterStatus;

    return(0);
}


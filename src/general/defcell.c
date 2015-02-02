/*

Copyright (c) 1995 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
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
#include <stdlib.h>

#include "cell.h"
#include "plot.h"

#include "stdafx.h"

/*  cell dimensions */
static struct {
    int   display;
    float a;
    float b;
    float c;
    float alpha;
    float beta;
    float gamma;
    float Xtrans;
    float Ytrans;
    float Ztrans;
} Cell  = { 0    ,
            1.e+20f, 1.e+20f, 1.e+20f,
              90.0f,   90.0f,   90.0f,
               0.0,     0.0,     0.0};

/*  structure to handle the display on/off toggle button status */
static struct {
    gom_PlotterData Plotter; /* default display off */
    int   linewidth;        /* default linewidth 1 pixel wide */
} PlotCell = { { NULL , NULL } , 1 };

#define CellDataIsChanging() \
    gomp_InvalidatePlotterDelayed(&PlotCell.Plotter)

#define IMPLEMENT_GETSET_CELL_VALUE(title,var) \
                                               \
    float gomp_GetCell##title(void)            \
    {                                          \
        return(Cell.var);                      \
    }                                          \
                                               \
    int  gomp_SetCell##title(float Value)      \
    {                                          \
        CellDataIsChanging();                  \
        Cell.var = Value;                      \
        return(0);                             \
    }

/***********************************************************************/
IMPLEMENT_GETSET_CELL_VALUE(A,a)
IMPLEMENT_GETSET_CELL_VALUE(B,b)
IMPLEMENT_GETSET_CELL_VALUE(C,c)
/***********************************************************************/
IMPLEMENT_GETSET_CELL_VALUE(Alpha,alpha)
IMPLEMENT_GETSET_CELL_VALUE(Beta,beta)
IMPLEMENT_GETSET_CELL_VALUE(Gamma,gamma)
/***********************************************************************/
IMPLEMENT_GETSET_CELL_VALUE(Xtrans,Xtrans)
IMPLEMENT_GETSET_CELL_VALUE(Ytrans,Ytrans)
IMPLEMENT_GETSET_CELL_VALUE(Ztrans,Ztrans)
/***********************************************************************/

/***********************************************************************/
int gomp_GetPlotCell()
/***********************************************************************/
{
    return(PlotCell.Plotter.plotter!=NULL);
}
/***********************************************************************/
int gomp_SetPlotCell(int State)
/***********************************************************************/
{
    return(gomp_SetPlotterRegistrationState(
        State,  &PlotCell.Plotter,
        gomp_PlotCellBox, NULL,
        PLOTTER_NAME_CELLBOX, PLOTTER_ORDER_CELLBOX));
}
/***********************************************************************/
int gomp_SetCellLinewidth(int Linewidth)
/***********************************************************************/
{
    CellDataIsChanging();
    PlotCell.linewidth = Linewidth;
    return(0);
}
/***********************************************************************/
int gomp_GetCellLinewidth()
/***********************************************************************/
{
    return(PlotCell.linewidth);
}

/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2001 - 2005 by:
Eero Häkkinen

Juha Ruokolainen <jpr@csc.fi> Center for Scientific Computing 1995
Leif Laaksonen <Leif.Laaksonen@csc.fi> Center for Scientific Computing 1995
Eero Häkkinen Center for Scientific Computing 2001,2002

Tässä on tämmönen pieni pätkä joka piirtää semmottei liinoja ku
joskus näkee näissä proteiinikuvissa yms.
kutsu on:

gomp_DrawRibbon( np,x,y,z,Color,Width,dx,dy,dz,ddx,ddy,ddz )

dx,dy,dz on käyrän tangenttivektori
ddx,ddy,ddz on  vektori käyrän leveyssuuntaan (esim toka derva)

jos piirtelet samaa nauhaa monta kertaa,niin
tossa ohjelmassa kannattaa ehkä erottaa
toi normaalin laskenta & normalisointi erikseen
(tekee kerran ja tallettaa) ja piirto erikseen...

Siinä tube ohjelman esimerkissä muuten on pieni bugi
niiden dervojen laskussa, noi kakkosella jaot alla
on liikaa...


for( i=1; i<NP-1; i++ )
{
dx[i] = (x[i+1]-x[i-1])/2;
dy[i] = (y[i+1]-y[i-1])/2;
dz[i] = (z[i+1]-z[i-1])/2;
}

...

for( i=1; i<NP-1; i++ )
{
    ddx[i] = (dx[i+1]-dx[i-1])/2;
    ddy[i] = (dy[i+1]-dy[i-1])/2;
    ddz[i] = (dz[i+1]-dz[i-1])/2;
}

...

t. Juha

for( i=1; i<NP-1; i++ )
{
    dx[i] = (x[i+1]-x[i-1]);
    dy[i] = (y[i+1]-y[i-1]);
    dz[i] = (z[i+1]-z[i-1]);
}

dx[i] =  3*x[i] - 4*x[i-1] + x[i-2];
dx[0] = -3*x[0] +  4*x[1]  - x[2];

dy[i] =  3*y[i] - 4*y[i-1] + y[i-2];
dy[0] = -3*y[0] +  4*y[1]  - y[2];

dz[i] =  3*z[i] - 4*z[i-1] + z[i-2];
dz[0] = -3*z[0] +  4*z[1]  - z[2];

for( i=1; i<NP-1; i++ )
{
    ddx[i] = (dx[i+1]-dx[i-1]);
    ddy[i] = (dy[i+1]-dy[i-1]);
    ddz[i] = (dz[i+1]-dz[i-1]);
}
ddx[i] =  3*dx[i] - 4*dx[i-1] + dx[i-2];
ddx[0] = -3*dx[0] +  4*dx[1]  - dx[2];

ddy[i] =  3*dy[i] - 4*dy[i-1] + dy[i-2];
ddy[0] = -3*dy[0] +  4*dy[1]  - dy[2];

ddz[i] =  3*dz[i] - 4*dz[i-1] + dz[i-2];
ddz[0] = -3*dz[0] +  4*dz[1]  - dz[2];

*/

#include "maindefs.h"

#include "gomstdio.h"
#include "gommath.h"
#include <stdlib.h>
#include <string.h>

#if defined(WIN32)
#include <windows.h>
#endif

#include "gomfile.h"
#include "gomlistener.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "plumber.h"
#include "plot.h"
#include "printmsg.h"

#include "stdafx.h"

#define DEFAULT_WIDTH_SCALE 1.0
#define SPLINE_POINTS        50

static void DestroyPlumberDataType(void*,size_t,const DataVectorHandler *);

/* Zero and NULL init */
static const PlumberDataType NullPlumberDataType = { 0 };
static const DataVectorHandler PlumberDataVectorHandler = {
    gomp_DataVectorCopyInitFunc,
    DestroyPlumberDataType,
    &NullPlumberDataType
};

static struct {
    PlumberDataType *data;
    gom_PlotterData  Plotter;
    gom_AtomCoordinateChangedListener *AtomCoordinateChangedListener;
    gom_MolecStructDeleteListener     *MolecStructDeleteListener;
    double *dx;
    double *dy;
    double *dz;
} Plumber;

static int ResetMolecStructureSpecificData(void *, int, int);
static int ReserveSpaceForModelPlumber(int, int, int);
static int ReserveSpaceForPlumber(int, int, int);
static int CalculatePlumberDC(int);
static int CalculatePlumberNWC(int , int);
static int CalculatePlumberHelixNWC(int);
static int FixPlumberNWC(int);
static int UpdatePlumberData(int, int, const int *);
static int UpdateGluedPlumbers(void *, unsigned long int);

#define PlumberDataIsChanging() \
    gomp_InvalidatePlotterDelayed(&Plumber.Plotter)

/**************************************************************************/
static void DestroyPlumberDataType(
    void *data,size_t size,const DataVectorHandler *pHandler)
/**************************************************************************/
{
    PlumberDataType *item = data;

    /* free vertex tables */
    free(item->atom_list);

    free(item->x);
    free(item->y);
    free(item->z);
    
    free(item->NX);
    free(item->NY);
    free(item->NZ);
    
    free(item->WX);
    free(item->WY);
    free(item->WZ);
    
    free(item->width);
    free(item->thickness);
}
/**************************************************************************/
int ResetMolecStructureSpecificData(void *userData, int Wstr, int Dstr)
/**************************************************************************/
{
    if ( Wstr < 0 )
        /* gOpenMol is about to be reset. */
        gomp_DeletePlumbers();
    else {
        int i;
        for ( i = gomp_GetPlumberSets() - 1 ; i >= 0 ; i-- ) {
            if ( Plumber.data[i].Wstr == Wstr ) {
                PlumberDataIsChanging();

                if ( Dstr < 0 )
                    /* Structure is about to be deleted. */
                    gomp_DeletePlumber( i );
                else {
                    /* Structure Wstr is about to be merged into Dstr. */
                    int j;
                    int NAtoms = gomp_GetNumAtomsInMolecStruct(Dstr);
                    Plumber.data[i].Wstr = Dstr;
                    for ( j = Plumber.data[i].Atoms - 1 ; j >= 0 ; j-- )
                        Plumber.data[i].atom_list[j] += NAtoms;
                }
            }
            if ( Plumber.data[i].Wstr > Wstr )
                --Plumber.data[i].Wstr;
        }
    }
    return(gomp_GetPlumberSets() ? 0 : -1);
}
/**************************************************************************/
int gomp_DeletePlumber(int Set)
/**************************************************************************/
{
    if( gomp_DataVectorGetSize(&Plumber.data) <= 1 )
        return gomp_DeletePlumbers();

    gomp_SetPlumberStructure(Set, -1);
    (void)gomp_DataVectorRemove(&Plumber.data, Set, 1);
    (void)gomp_DataVectorSetTotalSize(&Plumber.data,0);

    PlumberDataIsChanging();

    return(0);
}

/**************************************************************************/
int gomp_DeletePlumbers()
/**************************************************************************/
{
    PlumberDataIsChanging();

    gomp_DataVectorFree(&Plumber.data);

    gomp_SetPlumberDisplay(0);

    gomp_CancelAtomCoordinateChangedListener(
        Plumber.AtomCoordinateChangedListener);
    Plumber.AtomCoordinateChangedListener = NULL;

    return(0);
}

/**************************************************************************/
int ReserveSpaceForModelPlumber(int Points, int Atoms, int Type)
/**************************************************************************/
{
    int i;

    if ( ! gomp_DataVectorCreateOrAppend(
        &Plumber.data, &PlumberDataVectorHandler,1) )
        return(-1);

    if ( ! Plumber.MolecStructDeleteListener )
        Plumber.MolecStructDeleteListener =
            gomp_AddMolecStructDeleteListener(
                ResetMolecStructureSpecificData, NULL);

    i = gomp_DataVectorGetSize(&Plumber.data) - 1;

    if( Type != TRACE_TYPE ) {
        if ( ! ( Plumber.data[i].width     = gomp_AllocateDoubleVector(Points) ) ||
             ! ( Plumber.data[i].thickness = gomp_AllocateDoubleVector(Points) ) ||
             ! ( Plumber.data[i].NX        = gomp_AllocateDoubleVector(Points) ) ||
             ! ( Plumber.data[i].NY        = gomp_AllocateDoubleVector(Points) ) ||
             ! ( Plumber.data[i].NZ        = gomp_AllocateDoubleVector(Points) ) ||
             ! ( Plumber.data[i].WX        = gomp_AllocateDoubleVector(Points) ) ||
             ! ( Plumber.data[i].WY        = gomp_AllocateDoubleVector(Points) ) ||
             ! ( Plumber.data[i].WZ        = gomp_AllocateDoubleVector(Points) ) ) {
            (void)gomp_DataVectorSetExactSize(&Plumber.data,i);
            return(-1);
        }
    }

    Plumber.data[i].IsUpdated = 1;
    Plumber.data[i].Atoms     = Atoms;
    if ( Atoms ) {
        Plumber.data[i].atom_list = gomp_AllocateIntVector(Atoms);
        if ( ! Plumber.AtomCoordinateChangedListener )
            Plumber.AtomCoordinateChangedListener =
                gomp_AddAtomCoordinateChangedListener(
                    UpdateGluedPlumbers, NULL);
    }

    Plumber.data[i].Points = Points;
    Plumber.data[i].Type   = Type;

    if ( ! ( Plumber.data[i].x      = gomp_AllocateDoubleVector(Points) ) ||
         ! ( Plumber.data[i].y      = gomp_AllocateDoubleVector(Points) ) ||
         ! ( Plumber.data[i].z      = gomp_AllocateDoubleVector(Points) ) ) {
        (void)gomp_DataVectorSetExactSize(&Plumber.data,i);
        return(-1);
    }

    return(0);
}

/**************************************************************************/
int ReserveSpaceForPlumber(int Atoms, int Type, int Glue)
/**************************************************************************/
{
    int NP;

    if( Type == TRACE_TYPE )
        NP = Atoms;
    else
        NP = SPLINE_POINTS * (Atoms - 1) + 1;

    return ReserveSpaceForModelPlumber(NP,Glue ? Atoms : 0,Type);
}

/**************************************************************************/
const PlumberDataType *gomp_GetPlumbers(void)
/**************************************************************************/
{
    return(Plumber.data);
}
/**************************************************************************/
int gomp_GetPlumberSets()
/**************************************************************************/
{
    return(gomp_DataVectorGetSize(&Plumber.data));
}

/**************************************************************************/
const double  *gomp_GetPlumberXp(int Position)
/**************************************************************************/
{
    return(Plumber.data[Position].x);
}
/**************************************************************************/
const double  *gomp_GetPlumberYp(int Position)
/**************************************************************************/
{
    return(Plumber.data[Position].y);
}

/**************************************************************************/
const double  *gomp_GetPlumberZp(int Position)
/**************************************************************************/
{
    return(Plumber.data[Position].z);
}

/**************************************************************************/
float gomp_GetPlumberRed(int Position)
/**************************************************************************/
{
    return(Plumber.data[Position].red);
}

/**************************************************************************/
float gomp_GetPlumberGreen(int Position)
/**************************************************************************/
{
    return(Plumber.data[Position].green);
}

/**************************************************************************/
float gomp_GetPlumberBlue(int Position)
/**************************************************************************/
{
    return(Plumber.data[Position].blue);
}

/**************************************************************************/
int  gomp_GetPlumberStructure(int Position)
/**************************************************************************/
{
    return(Plumber.data[Position].Wstr);
}

/**************************************************************************/
int  gomp_GetPlumberAtoms(int Position)
/**************************************************************************/
{
    return(Plumber.data[Position].Atoms);
}

/**************************************************************************/
const int * gomp_GetPlumberAtomList(int Position)
/**************************************************************************/
{
    return(Plumber.data[Position].atom_list);
}

/**************************************************************************/
int    CalculatePlumberDC(int Set)
/**************************************************************************/
{
    const double  *x;
    const double  *y;
    const double  *z;
    int    NP;
    int    i;

    if ( (size_t)Set >= gomp_DataVectorGetSize(&Plumber.data))
        return(1);

    x  = gomp_GetPlumberXp(Set);
    y  = gomp_GetPlumberYp(Set);
    z  = gomp_GetPlumberZp(Set);
    NP = Plumber.data[Set].Points;

    for( i=1; i<NP-1; i++ )
    {
        Plumber.dx[i] = x[i+1]-x[i-1];
        Plumber.dy[i] = y[i+1]-y[i-1];
        Plumber.dz[i] = z[i+1]-z[i-1];
    }

    Plumber.dx[i] =  3*x[i] - 4*x[i-1] + x[i-2];
    Plumber.dx[0] = -3*x[0] +  4*x[1]  - x[2];

    Plumber.dy[i] =  3*y[i] - 4*y[i-1] + y[i-2];
    Plumber.dy[0] = -3*y[0] +  4*y[1]  - y[2];

    Plumber.dz[i] =  3*z[i] - 4*z[i-1] + z[i-2];
    Plumber.dz[0] = -3*z[0] +  4*z[1]  - z[2];

    return(0);

}
/**************************************************************************/
int    CalculatePlumberNWC(int Set , int Type)
/**************************************************************************/
{
    const double  *dx;
    const double  *dy;
    const double  *dz;
    double s,*NX,*NY,*NZ,*WX,*WY,*WZ;
    double Ax,Ay,Az;
    int    NP;
    int    i;

    if( (size_t)Set >= gomp_DataVectorGetSize(&Plumber.data) ) return(1);

    dx = Plumber.dx;
    dy = Plumber.dy;
    dz = Plumber.dz;
    NP = Plumber.data[Set].Points;
    NX = Plumber.data[Set].NX;
    NY = Plumber.data[Set].NY;
    NZ = Plumber.data[Set].NZ;
    WX = Plumber.data[Set].WX;
    WY = Plumber.data[Set].WY;
    WZ = Plumber.data[Set].WZ;


    /* Vector [NX NY NZ] defines direction of the ribbon normal.    */
    /* Vector [dx dy dz] is the direction of the ribbon. So its     */
    /* derivate is parallel to the ribbon normal.                   */

    for( i=1; i<NP-1; i++ )
    {
        NX[i] = (dx[i+1]-dx[i-1]);
        NY[i] = (dy[i+1]-dy[i-1]);
        NZ[i] = (dz[i+1]-dz[i-1]);
    }

    /* We suppose that derivate of NX is approximately constant.    */
    /* From linear approximation we get:                            */
    /*     NX[-1]  = 3*NX[0] - 3*NX[1] + NX[2]                      */
    /*     NX[i+1] = 3*NX[i] - 3*NX[i-1] + NX[i-2]                  */
    NX[i] =  3.0*dx[i] - 4.0*dx[i-1] + dx[i-2];
    NX[0] = -3.0*dx[0] + 4.0*dx[1]   - dx[2];

    NY[i] =  3.0*dy[i] - 4.0*dy[i-1] + dy[i-2];
    NY[0] = -3.0*dy[0] + 4.0*dy[1]   - dy[2];

    NZ[i] =  3.0*dz[i] - 4.0*dz[i-1] + dz[i-2];
    NZ[0] = -3.0*dz[0] + 4.0*dz[1]   - dz[2];

    Ax=Ay=Az=0.0;
    for( i=0; i<NP; i++ )
    {
        /* Vector [WX WY WZ] is the second crossection vector.
           It also is a width vector for a ribbon.
           Use cross production to calculate that vector. */
        WX[i] = dy[i]*NZ[i] - dz[i]*NY[i];
        WY[i] = dz[i]*NX[i] - dx[i]*NZ[i];
        WZ[i] = dx[i]*NY[i] - dy[i]*NX[i];

        s  = WX[i]*WX[i];
        s += WY[i]*WY[i];
        s += WZ[i]*WZ[i];
        if( s > 0.0 )
            s  = 1.0 / sqrt( s );

        WX[i] *= s;
        WY[i] *= s;
        WZ[i] *= s;

        /* Calculate the average [WX WY WZ] vector. */
        Ax += WX[i];
        Ay += WY[i];
        Az += WZ[i];
    }

    s  = Ax*Ax;
    s += Ay*Ay;
    s += Az*Az;
    s  = 1.0 / sqrt( s );

    Ax *= s;
    Ay *= s;
    Az *= s;

    for( i=0; i<NP; ++i )
    {
        /* Direct each [WX WY WZ] vector towards */
        /* the average vector.                   */
        WX[i] = Ax;
        WY[i] = Ay;
        WZ[i] = Az;

        /* Recalculate [NX NY NZ] vectors.               */
        /* The order of cross production is significant. */
        /* If you change the order you have to multiply  */
        /* one of the s'es by -1.                        */
        NX[i] = WY[i]*dz[i] - WZ[i]*dy[i];
        NY[i] = WZ[i]*dx[i] - WX[i]*dz[i];
        NZ[i] = WX[i]*dy[i] - WY[i]*dx[i];

        s  = NX[i]*NX[i];
        s += NY[i]*NY[i];
        s += NZ[i]*NZ[i];
        s  = 1.0 / sqrt( s );

        NX[i] *= s;
        NY[i] *= s;
        NZ[i] *= s;
    }

/* we don't need these anymore ... */
    free(Plumber.dx);
    free(Plumber.dy);
    free(Plumber.dz);
    Plumber.dx=Plumber.dy=Plumber.dz=(double  *)NULL;

    return(0);
}
/**************************************************************************/
int    CalculatePlumberHelixNWC(int Set)
/**************************************************************************/
{
    const double  *x;
    const double  *y;
    const double  *z;
    const double  *dx;
    const double  *dy;
    const double  *dz;
    double s,*NX,*NY,*NZ,*WX,*WY,*WZ;
    int    NP;
    int    i;
    /* helix has a hydrogen bond between every 4th amino acid. */
    const int points_per_round=4*SPLINE_POINTS;

    if((size_t)Set >= gomp_DataVectorGetSize(&Plumber.data)) return(1);

    x  = gomp_GetPlumberXp(Set);
    y  = gomp_GetPlumberYp(Set);
    z  = gomp_GetPlumberZp(Set);
    dx = Plumber.dx;
    dy = Plumber.dy;
    dz = Plumber.dz;
    NP = Plumber.data[Set].Points;
    NX = Plumber.data[Set].NX;
    NY = Plumber.data[Set].NY;
    NZ = Plumber.data[Set].NZ;
    WX = Plumber.data[Set].WX;
    WY = Plumber.data[Set].WY;
    WZ = Plumber.data[Set].WZ;


    /* Vector [WX WY WZ] defines direction of the helix width.  */
    /* Define it to be parallel to the helix axis.                 */
    WX[0] = x[NP-1] + x[NP-1-points_per_round/2] - x[points_per_round/2] - x[0];
    WY[0] = y[NP-1] + y[NP-1-points_per_round/2] - y[points_per_round/2] - y[0];
    WZ[0] = z[NP-1] + z[NP-1-points_per_round/2] - z[points_per_round/2] - z[0];

    s  = WX[0]*WX[0];
    s += WY[0]*WY[0];
    s += WZ[0]*WZ[0];
    s  = 1.0 / sqrt( s );

    WX[0] = s*WX[0];
    WY[0] = s*WY[0];
    WZ[0] = s*WZ[0];

    for( i=1; i<NP; ++i)
    {
        WX[i] = WX[0];
        WY[i] = WY[0];
        WZ[i] = WZ[0];
    }

    for( i=0; i<NP; ++i ) {
        /* Vector [NX NY NZ] is the second crossection vector.
           It also is a normal vector for flat helicies.
           Use cross production to calculate that vector. */
        NX[i] = WY[i]*dz[i] - WZ[i]*dy[i];
        NY[i] = WZ[i]*dx[i] - WX[i]*dz[i];
        NZ[i] = WX[i]*dy[i] - WY[i]*dx[i];

        s  = NX[i]*NX[i];
        s += NY[i]*NY[i];
        s += NZ[i]*NZ[i];
        s  = 1.0 / sqrt( s );

        NX[i] = s*NX[i];
        NY[i] = s*NY[i];
        NZ[i] = s*NZ[i];
    }

/* we don't need them anymore ... */
    free(Plumber.dx);
    free(Plumber.dy);
    free(Plumber.dz);
    Plumber.dx=Plumber.dy=Plumber.dz=(double  *)NULL;

    return(0);
}
/**************************************************************************/
int    FixPlumberNWC(int Set)
/**************************************************************************/
{
    double  *NX,*NY,*NZ,*WX,*WY,*WZ;
    int    NP;
    int    i;

    if((size_t)Set >= gomp_DataVectorGetSize(&Plumber.data)) return(1);

    NP = Plumber.data[Set].Points;
    NX = Plumber.data[Set].NX;
    NY = Plumber.data[Set].NY;
    NZ = Plumber.data[Set].NZ;
    WX = Plumber.data[Set].WX;
    WY = Plumber.data[Set].WY;
    WZ = Plumber.data[Set].WZ;


    /* Let's smooth the ribbon. */
    for( i=1; i<NP; ++i )
    {
        if( NX[i-1]*NX[i] + NY[i-1]*NY[i] + NZ[i-1]*NZ[i] < 0.0 ) {
            NX[i]=-NX[i];
            NY[i]=-NY[i];
            NZ[i]=-NZ[i];

            WX[i]=-WX[i];
            WY[i]=-WY[i];
            WZ[i]=-WZ[i];
        }
    }

    return(0);
}
/**************************************************************************/
int    gomp_SetPlumberDisplay(int StatusValue)
/**************************************************************************/
{
    PlumberDataIsChanging();

    return(gomp_SetPlotterRegistrationState(
        StatusValue, &Plumber.Plotter, gomp_PlotPlumber,NULL,
        PLOTTER_NAME_PLUMBER, PLOTTER_ORDER_PLUMBER));
}
/**************************************************************************/
int    gomp_GetPlumberDisplay()
/**************************************************************************/
{
    return(Plumber.Plotter.plotter!=NULL);
}

/**************************************************************************/
int    gomp_SetPlumberRed(int set, float value)
/**************************************************************************/
{
    PlumberDataIsChanging();

    Plumber.data[set].red = value;

    return(0);
}

/**************************************************************************/
int    gomp_SetPlumberGreen(int set, float value)
/**************************************************************************/
{
    PlumberDataIsChanging();

    Plumber.data[set].green = value;

    return(0);
}

/**************************************************************************/
int    gomp_SetPlumberBlue(int set, float value)
/**************************************************************************/
{
    PlumberDataIsChanging();

    Plumber.data[set].blue = value;

    return(0);
}

/**************************************************************************/
int gomp_SetPlumberStructure(int set,int Wstr)
/**************************************************************************/
{
    PlumberDataIsChanging();

    if ( Plumber.data[set].Atoms ) {
        int i;
        /* Last chance to update the plumber. */
        gomp_UpdateData();
        free(Plumber.data[set].atom_list);
        Plumber.data[set].Atoms     = 0;
        Plumber.data[set].atom_list = NULL;
        /* Check if there still is need for the listener. */
        for ( i = 0 ; i < gomp_GetPlumberSets() ; i++ ) {
            if ( Plumber.data[i].Atoms )
                break;
        }
        if ( i >= gomp_GetPlumberSets() ) {
            /* No need for the listener. */
            gomp_CancelAtomCoordinateChangedListener(
                Plumber.AtomCoordinateChangedListener);
            Plumber.AtomCoordinateChangedListener = NULL;
        }
    }

    Plumber.data[set].Wstr = Wstr;

    return(0);
}

/**************************************************************************/
int UpdateGluedPlumbers(void *userData, unsigned long int change_mask)
/**************************************************************************/
{
    int i, Sets;
    unsigned long int bit;

    Sets = gomp_DataVectorGetSize(&Plumber.data);

    for ( i = 0 ; i < Sets ; i++ ) {
        bit   = 1;
        bit <<= Plumber.data[i].Wstr;
        if ( change_mask & bit || bit == 0 ) {
            UpdatePlumberData(
                i, Plumber.data[i].Atoms, Plumber.data[i].atom_list);
            PlumberDataIsChanging();
        }
    }

    return(0);
}

/**************************************************************************/
int gomp_LoadPlumberAtoms(int Wstr , int slong, const int *sel_list ,
                          float red , float green , float blue ,
                          float FRad, int Type, float Width, float Thickness ,
                          int Glue)
/**************************************************************************/
{
    int    i;
    int    Set;
    int    NP;
    double param;

    switch(Type) {
    case RIBBON_TYPE:
    case CYLINDER_TYPE:
    case TRACE_TYPE:
        if( slong < 2 ) {
            gomp_PrintERROR("At least 2 points have to be defined in the ribbon atom list");
            return(1);
        }
        break;
    default:
        if ( slong < 3 ) {
            gomp_PrintERROR("At least 3 points have to be defined in the atom list");
            return(1);
        }
    }

    if ( ReserveSpaceForPlumber(slong,Type,Glue) )
        return(1);

    Set = gomp_DataVectorGetSize(&Plumber.data) - 1;
    NP  = Plumber.data[Set].Points;

    if( Glue )
        memcpy( Plumber.data[Set].atom_list,
                sel_list, slong*sizeof(*sel_list) );        


    Plumber.data[Set].Wstr = Wstr;

    (void)gomp_SetPlumberRed(   Set , red  );
    (void)gomp_SetPlumberGreen( Set , green);
    (void)gomp_SetPlumberBlue(  Set , blue );
    (void)gomp_SetPlumberDisplayType( Set , Type);

    switch( Type )
    {
    case FLAT_HELIX_TYPE:
    case SOLID_HELIX_TYPE:
    case STRAND_TYPE:
        for(i=SPLINE_POINTS; i< NP-SPLINE_POINTS ; ++i) {
            Plumber.data[Set].width[i]     = Width;
            Plumber.data[Set].thickness[i] = Thickness;
        }
        /* Create a smoothing curve to the both ends of the plumber. */
        for(i=0; i<SPLINE_POINTS; ++i) {
            /* param increases from 0 to 1 */
            param = 0.5-0.5*cos( i*M_PI/SPLINE_POINTS );

            Plumber.data[Set].width[i]          = (1.0-param)*FRad + param*Width;
            Plumber.data[Set].width[NP-i-1]     = (1.0-param)*FRad + param*Width;
            Plumber.data[Set].thickness[i]      = (1.0-param)*FRad + param*Thickness;
            Plumber.data[Set].thickness[NP-i-1] = (1.0-param)*FRad + param*Thickness;;
        }
        break;
    case ARROW_TYPE:
        for(i=0; i< NP ; ++i) {
            Plumber.data[Set].width[i]      = Width;
            Plumber.data[Set].thickness[i]  = Thickness;
        }
        /* Create an array head to the last end of the sheet.    */
        for(i=0; i<SPLINE_POINTS; ++i)
            Plumber.data[Set].width[NP-i-1] = (2.0*i*Width)/SPLINE_POINTS;
        break;
    case TRACE_TYPE:
        break;
    default:
        for(i = 0 ; i < NP ; i++) {
            Plumber.data[Set].width[i]      = FRad;
            Plumber.data[Set].thickness[i]  = FRad;
        }
    }

    if ( UpdatePlumberData(Set, slong, sel_list) != 0 ) {
        /* out of memory */
        (void)gomp_DataVectorSetExactSize(&Plumber.data,Set);
        return(-1);
    }

    PlumberDataIsChanging();

    return(0);
}

/**************************************************************************/
int UpdatePlumberData(int Set , int slong, const int *sel_list)
/**************************************************************************/
{
    const float *x;
    const float *y;
    const float *z;
    int         *d;
    int          i;
    int         NP;

    NP = Plumber.data[Set].Points;

    if ( ! ( Plumber.dx = gomp_AllocateDoubleVector(NP) ) ||
         ! ( Plumber.dy = gomp_AllocateDoubleVector(NP) ) ||
         ! ( Plumber.dz = gomp_AllocateDoubleVector(NP) ) ||
         ! ( d          = gomp_AllocateIntVector(NP) ) )
        return(-1);

    PlumberDataIsChanging();

    x = gomp_GetAtomXCoordPointer(Plumber.data[Set].Wstr);
    y = gomp_GetAtomYCoordPointer(Plumber.data[Set].Wstr);
    z = gomp_GetAtomZCoordPointer(Plumber.data[Set].Wstr);

    switch( Plumber.data[Set].Type ) {
    case ARROW_TYPE:
    case STRAND_TYPE:
        d[0] = 0;
        Plumber.data[Set].x[d[0]] = x[sel_list[0]];
        Plumber.data[Set].y[d[0]] = y[sel_list[0]];
        Plumber.data[Set].z[d[0]] = z[sel_list[0]];

        /* Points in the middle of the ribbon can't go along atoms because */
        /* arrow must not contain many curves.                             */
        for(i=1 ; i<slong-1 ; ++i) {
            d[i] = SPLINE_POINTS * i;
            Plumber.data[Set].x[d[i]] =
                ( x[sel_list[i-1]] + 2.0*x[sel_list[i]] + x[sel_list[i+1]] )/4.0;
            Plumber.data[Set].y[d[i]] =
                ( y[sel_list[i-1]] + 2.0*y[sel_list[i]] + y[sel_list[i+1]] )/4.0;
            Plumber.data[Set].z[d[i]] =
                ( z[sel_list[i-1]] + 2.0*z[sel_list[i]] + z[sel_list[i+1]] )/4.0;
        }

        d[i] = SPLINE_POINTS * i;
        Plumber.data[Set].x[d[i]] = x[sel_list[i]];
        Plumber.data[Set].y[d[i]] = y[sel_list[i]];
        Plumber.data[Set].z[d[i]] = z[sel_list[i]];
        break;
    case TRACE_TYPE:
        for(i = 0 ; i < slong ; i++) {
            Plumber.data[Set].x[i] = x[sel_list[i]];
            Plumber.data[Set].y[i] = y[sel_list[i]];
            Plumber.data[Set].z[i] = z[sel_list[i]];
        }
        return(0);
    default:
        for(i = 0 ; i < slong ; i++) {
            d[i] = SPLINE_POINTS * i;
            Plumber.data[Set].x[d[i]] = x[sel_list[i]];
            Plumber.data[Set].y[d[i]] = y[sel_list[i]];
            Plumber.data[Set].z[d[i]] = z[sel_list[i]];
        }
    }

    if( slong < 3 ) {
        (void)gomp_LinearSpline( Plumber.data[Set].x, slong, d );
        (void)gomp_LinearSpline( Plumber.data[Set].y, slong, d );
        (void)gomp_LinearSpline( Plumber.data[Set].z, slong, d );
    }
    else {
        (void)gomp_CubicSpline( Plumber.data[Set].x, slong, d );
        (void)gomp_CubicSpline( Plumber.data[Set].y, slong, d );
        (void)gomp_CubicSpline( Plumber.data[Set].z, slong, d );
    }

    (void)CalculatePlumberDC(Set);
    switch( Plumber.data[Set].Type )
    {
    case FLAT_HELIX_TYPE:
    case SOLID_HELIX_TYPE:
        if(slong >= 4)
        {
            /* We need at least 4 points to determine helix axis. */
            (void)CalculatePlumberHelixNWC(Set);
            break;
        }
        /* Continue to the default segment. */
    default:
        (void)CalculatePlumberNWC(Set , Plumber.data[Set].Type);
    }

    free(d);
    
    return(0);
}

/**************************************************************************/
int    gomp_SetPlumberDisplayType(int Set , int Value)
/**************************************************************************/
{
    Plumber.data[Set].Type = Value;

    PlumberDataIsChanging();

    return(0);
}
/**************************************************************************/
int    gomp_GetPlumberDisplayType(int Set)
/**************************************************************************/
{
    return(Plumber.data[Set].Type);
}

/**************************************************************************/
int gomp_RetrievePlumberInfoFromModelFile(FILE *Model_f)
/**************************************************************************/
{
    int   i,j;
    int   Sets;
    int   Local_SplinePoints;
    int   PDisplay;
    int   PType;
    int   PPoints;
    int   PLocation;
    int   PAtoms;
    int   PWstr;
    float Pred;
    float Pgreen;
    float Pblue;
    char  InputText[BUFF_LEN];
    int    PIndex;
    double Pwidth;
    double Pthickness;
    double Temp1;
    double Temp2;
    double Temp3;
    double Temp4;
    double Temp5;
    double Temp6;
    double Temp7;
    double Temp8;
    double Temp9;


/* here we go ... */

    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText, "%d %d %d",&Sets,&PDisplay,&Local_SplinePoints);

    if ( Local_SplinePoints != SPLINE_POINTS )
        gomp_PrintEXIT("internal spline point error, number of points does not match");

    (void)gomp_DeletePlumbers();

    for ( i = 0 ; i < Sets ; i++ ) {
        
        /* for backwards compatibility */
        PAtoms = PWstr = 0;

        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText, "%d %d %d %f %f %f %d %d %d",
               &PType    ,
               &PPoints  ,
               &PLocation,
               &Pred     ,
               &Pgreen   ,
               &Pblue    ,
               &Local_SplinePoints,
               &PAtoms   ,
               &PWstr);

        if ( ReserveSpaceForModelPlumber(PPoints,PAtoms,PType) )
            return(1);

        Plumber.data[i].red     = Pred;
        Plumber.data[i].green   = Pgreen;
        Plumber.data[i].blue    = Pblue;
        Plumber.data[i].Wstr    = PWstr;

        if( PAtoms ) {
            /* read atom index list */
            for(j = 0 ; j < PAtoms ; j++) {
                gomp_Fgets(InputText,BUFF_LEN,Model_f);
                sscanf(InputText , "%d\n",&PIndex);
                Plumber.data[i].atom_list[j] = PIndex;
            }
            if( Plumber.data[i].Type != TRACE_TYPE ) {
                /* read width and optionally thickness */
                for(j = 0 ; j < PPoints ; j++) {
                    gomp_Fgets(InputText,BUFF_LEN,Model_f);
                    switch(Plumber.data[i].Type) {
                    case SOLID_HELIX_TYPE:
                    case ARROW_TYPE:
                    case STRAND_TYPE:
                        sscanf(InputText , "%lf %lf\n",&Pwidth,&Pthickness);
                        Plumber.data[i].width[j]     = Pwidth;
                        Plumber.data[i].thickness[j] = Pthickness;
                        break;
                    default:
                        sscanf(InputText , "%lf\n",    &Pwidth);
                        Plumber.data[i].width[j]     = Pwidth;
                        Plumber.data[i].thickness[j] = Pwidth;
                    }
                }
            }
            /* plumber coordinates have to be recalculated */
            Plumber.data[i].IsUpdated = 0;
        }
        else {
            for(j = 0 ; j < PPoints ; j++) {
                gomp_Fgets(InputText,BUFF_LEN,Model_f);
                switch(Plumber.data[i].Type) {
                case SOLID_HELIX_TYPE:
                case ARROW_TYPE:
                case STRAND_TYPE:
                    sscanf(InputText , "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                           &Pwidth,&Pthickness,
                           &Temp1,&Temp2,&Temp3,
                           &Temp4,&Temp5,&Temp6,
                           &Temp7,&Temp8,&Temp9);
                    Plumber.data[i].width[j]     = Pwidth;
                    Plumber.data[i].thickness[j] = Pthickness;
                    break;
                default:
                    sscanf(InputText , "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                           &Pwidth,
                           &Temp1,&Temp2,&Temp3,
                           &Temp4,&Temp5,&Temp6,
                           &Temp7,&Temp8,&Temp9);                
                    Plumber.data[i].width[j]     = Pwidth;
                    Plumber.data[i].thickness[j] = Pwidth;
                }
                Plumber.data[i].x[j]     = Temp1;
                Plumber.data[i].y[j]     = Temp2;
                Plumber.data[i].z[j]     = Temp3;
                Plumber.data[i].NX[j]    = Temp4;
                Plumber.data[i].NY[j]    = Temp5;
                Plumber.data[i].NZ[j]    = Temp6;
                Plumber.data[i].WX[j]    = Temp7;
                Plumber.data[i].WY[j]    = Temp8;
                Plumber.data[i].WZ[j]    = Temp9;
            }
            FixPlumberNWC( i );
            /* plumber coordinates have been read from the file */
            Plumber.data[i].IsUpdated = 1;
        }
    }

    gomp_SetPlumberDisplay( PDisplay );
    
    PlumberDataIsChanging();

    return(0);
}

/**************************************************************************/
int gomp_StorePlumberInfo2ModelFile(FILE *Model_f)
/**************************************************************************/
{
    int   i,j,Sets;
    int Display;


/* here we go ... */
    Sets = gomp_DataVectorGetSize(&Plumber.data);
    if(!Sets) return(0);

/* PLUMBER - tag */
    fprintf(Model_f , "[Plumber]\n");
    Display = Plumber.Plotter.plotter != NULL;
    fprintf(Model_f , "%d %d %d\n",Sets   ,
            Display,
            SPLINE_POINTS);

    for(i = 0 ; i < Sets ; i++) {
        fprintf(Model_f , "%d %d %d %f %f %f %d %d %d\n",
                Plumber.data[i].Type    ,
                Plumber.data[i].Points  ,
                0                      , /* Location */
                Plumber.data[i].red     ,
                Plumber.data[i].green   ,
                Plumber.data[i].blue    ,
                SPLINE_POINTS          ,
                Plumber.data[i].Atoms   ,
                Plumber.data[i].Wstr);

        if( Plumber.data[i].Atoms ) {
            for(j = 0 ; j < Plumber.data[i].Atoms ; j++)
                fprintf(Model_f , "%d\n" , Plumber.data[i].atom_list[j]);
            if( Plumber.data[i].Type != TRACE_TYPE ) {
                for(j = 0 ; j < Plumber.data[i].Points ; j++) {
                    switch(Plumber.data[i].Type) {
                    case SOLID_HELIX_TYPE:
                    case ARROW_TYPE:
                    case STRAND_TYPE:
                        fprintf(Model_f , "%.10g %.10g\n",
                                Plumber.data[i].width[j],
                                Plumber.data[i].thickness[j]);
                        break;
                    default:
                        fprintf(Model_f , "%.10g\n",
                                Plumber.data[i].width[j]);
                    }
                }
            }
        }
        else {
            for(j = 0 ; j < Plumber.data[i].Points ; j++) {
                if( Plumber.data[i].Type != TRACE_TYPE ) {
                    fprintf(Model_f , "%.10g " , Plumber.data[i].width[j]);
                    switch(Plumber.data[i].Type) {
                    case SOLID_HELIX_TYPE:
                    case ARROW_TYPE:
                    case STRAND_TYPE:
                        fprintf(Model_f , "%.10g " , Plumber.data[i].thickness[j]);
                        break;
                    }
                }
                fprintf(Model_f ,
                        "%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n",
                        Plumber.data[i].x[j]  , Plumber.data[i].y[j]  , Plumber.data[i].z[j] ,
                        Plumber.data[i].NX[j] , Plumber.data[i].NY[j] , Plumber.data[i].NZ[j],
                        Plumber.data[i].WX[j] , Plumber.data[i].WY[j] , Plumber.data[i].WZ[j]);
            }
        }
    }
    
    return(0);
}

/*
  Enhancements 2001 - 2005 by:
  Eero Häkkinen
*/
/* From: Juha Ruokolainen <jpr@csc.fi> */

#include "maindefs.h"

#ifdef ENABLE_GRAPHICS
#if defined(WIN32)
#include <windows.h>
#endif

#if defined(GLUT)
#include <GL/glut.h>
#else
#include <GL/glaux.h>
#endif

#include <GL/gl.h>
#endif /* ENABLE_GRAPHICS */

#include "gomstdio.h"
#include "gommath.h"
#include <stdlib.h>
#include "gomstdio.h"
#include <math.h>

#include "colouring.h"
#include "memalloc.h"
#include "plot_plumber.h"
#include "plumber.h"
#include "plot.h"

#include "stdafx.h"

#define MAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MIN(x,y) ( (x) > (y) ? (y) : (x) )
#ifndef ABS
#  define ABS(x) ( (x) > 0 ? (x) : (-(x) ) )
#endif

#define AEPS 1.0e-15

#define ALLOCMEM( s ) calloc( s, 1 )
#define FREEMEM( p ) free( p ) 

#define A( i, j ) a[n * ( i ) + ( j )]

#if 0
static void ludecomp( double  *, int , int * );
static void lu_mtrinv( double  *a, int );

void lu_mtrinv( double  *a, int n )
{
    int i,j,k;
    int *pivot;

    double s;
  
    pivot = ALLOCMEM(n*sizeof(int));

    /*
     *  AP = LU
     */
    ludecomp( a,n,pivot );

    for( i=0; i<n; i++ )
    {
        if ( ABS(A(i,i))<AEPS )
        {
            fprintf( stderr, "Inv: Matrix is singular.\n" );
            return;
        }
        A(i,i) = 1.0/A(i,i);
    }

    /*  
     *  INV(U)
     */
    for( i=n-2; i>=0; i-- )
        for( j=n-1; j>i; j-- )
        {
            s = -A(i,j);
            for( k=i+1; k<j; k++ ) s -= A(i,k) * A(k,j);
            A(i,j) = s;
        }

    /*
     * INV(L)
     */
    for( i=n-2; i>=0; i-- )
        for( j=n-1; j>i; j-- )
        {
            s = 0.0;
            for( k=i+1; k<=j; k++ ) s -= A(j,k) * A(k,i);
            A(j,i) = A(i,i)*s;
        }
  
    /* 
     * A  = INV(AP)
     */
    for( i=0; i<n; i++ )
        for( j=0; j<n; j++ )
        {
            s = 0.0;
            for( k=MAX(i,j); k<n; k++ )
            {
                if ( k!=i )
                    s += A(i,k)*A(k,j);
                else
                    s += A(k,j);

                A(i,j) = s;
            }
        }

    /*
     * A = INV(A) (at last)
     */
    for( i=n-1; i>=0; i-- )
        if ( pivot[i]!=i )
            for( j=0; j<n; j++ )
            {
                s = A(i,j);
                A(i,j) = A(pivot[i],j);
                A(pivot[i],j) = s;
            }

    FREEMEM(pivot);
}

/*
 * LU- decomposition by gaussian elimination. Row pivoting is used.
 * 
 * result : AP = L'U ; L' = LD; pivot[i] is the swapped column number
 * for column i.
 *
 * Result is stored in place of original matrix.
 *
 */
static void ludecomp( double  *a, int n, int *pivot )
{
    double swap;
    int i,j,k,l;

    for( i=0; i<n-1; i++ )
    {
        j = i;
        for( k=i+1; k<n; k++ ) if ( ABS(A(i,k)) > ABS(A(i,j)) ) j = k;

        if ( ABS(A(i,j)) < AEPS )
        {
            fprintf( stderr, "LUDecomp: Matrix is (at least almost) singular, %d %d.\n", i, j );
            return;
        }

        pivot[i] = j;
    
        if ( j != i )
            for( k=0; k<=i; k++ )
            {
                swap = A(k,j);
                A(k,j) = A(k,i);
                A(k,i) = swap;
            }

        for( k=i+1; k<n; k++ ) A(i,k) /= A(i,i);
    
        for( k=i+1; k<n; k++ )
        {
            if ( j != i )
            {
                swap = A(k,i);
                A(k,i) = A(k,j); 
                A(k,j) = swap;
            }

            for( l=i+1; l<n; l++ ) A(k,l) -= A(k,i) * A(i,l);
        }
    }
  
    pivot[n-1] = n-1;
    if ( ABS(A(n-1,n-1))<AEPS )
    {
        fprintf( stderr, "LUDecomp: Matrix is (at least almost) singular.\n" );
    }
}
#endif


/* #define DEBUG_TUBE */

#define DEFAULT_RAD_SCALE 1.0

typedef struct { double x,y,z; } point_t;

#if 0
static void VectorToAngles(double vx,double vy,double vz,double  *px,double  *py,double  *pz,double r)
{
    double a,b,x,y,z,RadToDeg=180.0/M_PI;

    x = (*px-vx)*r;
    y = (*py-vy)*r;
    z = (vz-*pz)*r;

    b = atan2(y,z);
    r = y*sin(b) + z*cos(b);
    a = atan2(r,x);

    *px = RadToDeg*b;
    *py = RadToDeg*a;
    *pz = 0.0;
}
#endif

#ifdef ENABLE_GRAPHICS
static void ComputeTubeCircle
( const double  *x, const double  *y, const double  *z,
  const double  *NX,const double  *NY,const double  *NZ,
  const double  *WX,const double  *WY,const double  *WZ,
  const double  *Radius,
  point_t *Ellipse,point_t *Normal,int NPoints )
{
    double a,s,ca,sa;
    int i;

    for( i=0; i<NPoints-1; ++i )
    {
        a = 2*M_PI*i / (NPoints-1.0);
        ca = cos( a );
        sa = sin( a );

        Normal[i].x = WX[0]*ca + NX[0]*sa;
        Normal[i].y = WY[0]*ca + NY[0]*sa;
        Normal[i].z = WZ[0]*ca + NZ[0]*sa;

        s  = Normal[i].x*Normal[i].x;
        s += Normal[i].y*Normal[i].y;
        s += Normal[i].z*Normal[i].z;

        s = 1.0 / sqrt( s );

        Normal[i].x *= s;
        Normal[i].y *= s;
        Normal[i].z *= s;

        Ellipse[i].x = x[0] + DEFAULT_RAD_SCALE*Radius[0]*
            ( WX[0]*ca + NX[0]*sa );
        Ellipse[i].y = y[0] + DEFAULT_RAD_SCALE*Radius[0]*
            ( WY[0]*ca + NY[0]*sa );
        Ellipse[i].z = z[0] + DEFAULT_RAD_SCALE*Radius[0]*
            ( WZ[0]*ca + NZ[0]*sa );
    }

    Normal[i] = Normal[0];
    Ellipse[i] = Ellipse[0];
}

static void ComputeTubeEllipse(
    const double  *x, const double  *y, const double  *z,
    const double  *NX,const double  *NY,const double  *NZ,
    const double  *WX,const double  *WY,const double  *WZ,
    const double  *width,const double  * thickness,
    point_t *Ellipse,point_t *Normal,int NPoints )
{
    double a,s,ca,sa;
    int i;

    for( i=0; i<NPoints-1; ++i )
    {
        a = 2*M_PI*i / (NPoints-1.0);
        ca = cos( a );
        sa = sin( a );

        Normal[i].x = thickness[0]*WX[0]*ca + width[0]*NX[0]*sa;
        Normal[i].y = thickness[0]*WY[0]*ca + width[0]*NY[0]*sa;
        Normal[i].z = thickness[0]*WZ[0]*ca + width[0]*NZ[0]*sa;

        s  = Normal[i].x*Normal[i].x;
        s += Normal[i].y*Normal[i].y;
        s += Normal[i].z*Normal[i].z;

        s = 1.0 / sqrt( s );

        Normal[i].x *= s;
        Normal[i].y *= s;
        Normal[i].z *= s;

        Ellipse[i].x = x[0] + DEFAULT_RAD_SCALE*
            ( width[0]*WX[0]*ca + thickness[0]*NX[0]*sa );
        Ellipse[i].y = y[0] + DEFAULT_RAD_SCALE*
            ( width[0]*WY[0]*ca + thickness[0]*NY[0]*sa );
        Ellipse[i].z = z[0] + DEFAULT_RAD_SCALE*
            ( width[0]*WZ[0]*ca + thickness[0]*NZ[0]*sa );
    }

    Normal[i] = Normal[0];
    Ellipse[i] = Ellipse[0];
}
/*********************************************************************************/
static int DrawRectangleTubePlane(int NCurve,
                                  const double  *x,const double  *y,const double  *z,
                                  const double  *Color,
                                  const double  *Width,const double  *Thickness,
                                  const double  *NX,const double  *NY,const double  *NZ,
                                  int Ndir,
                                  const double  *WX,const double  *WY,const double  *WZ,
                                  int Wdir)
/*********************************************************************************/
{
    int i;

    glBegin( GL_QUAD_STRIP );
    for( i=0; i<NCurve; i++ )
    {

        glNormal3d( NX[i]*Ndir , NY[i]*Ndir , NZ[i]*Ndir );
        if ( Color ) glTexCoord1d( Color[i] );

        glVertex3d(
            x[i] + Thickness[i]*NX[i]*Ndir - Width[i]*WX[i]*Wdir,
            y[i] + Thickness[i]*NY[i]*Ndir - Width[i]*WY[i]*Wdir,
            z[i] + Thickness[i]*NZ[i]*Ndir - Width[i]*WZ[i]*Wdir
            );

        glNormal3d( NX[i]*Ndir , NY[i]*Ndir , NZ[i]*Ndir );
        if ( Color ) glTexCoord1d( Color[i] );

        glVertex3d(
            x[i] + Thickness[i]*NX[i]*Ndir + Width[i]*WX[i]*Wdir,
            y[i] + Thickness[i]*NY[i]*Ndir + Width[i]*WY[i]*Wdir,
            z[i] + Thickness[i]*NZ[i]*Ndir + Width[i]*WZ[i]*Wdir
            );
    }
    glEnd();

    return 0;
}

/*
 *
 */
/**************************************************************************/
int gomp_DrawRibbon
( 
    int NCurve,          /* number of curve points */
                         /* space coordinates of the curve                */
    const double  *x,const double  *y,const double  *z,
    const double  *Color, /* color coordinate of the curve (scaled to 0-1) */
    const double  *Width, /* Width scale of the curve                      */
    const double  *NX,const double  *NY,const double  *NZ,
    const double  *WX,const double  *WY,const double  *WZ
    )
/**************************************************************************/
{
    int i;

    glBegin( GL_QUAD_STRIP );
    for( i=0; i<NCurve; i++ )
    {
        glNormal3d( NX[i],NY[i],NZ[i] );
        if ( Color ) glTexCoord1d( Color[i] );

        glVertex3d(
            x[i] - Width[i]*WX[i],
            y[i] - Width[i]*WY[i],
            z[i] - Width[i]*WZ[i]
            );

        glNormal3d( NX[i],NY[i],NZ[i] );
        if ( Color ) glTexCoord1d( Color[i] );

        glVertex3d(
            x[i] + Width[i]*WX[i],
            y[i] + Width[i]*WY[i],
            z[i] + Width[i]*WZ[i]
            );
    }
    glEnd();

    return 0;
}

/*
 *
 */
/**************************************************************************/
int gomp_DrawTube
( 
    int NCurve,int NCircle, /* number of curve points, number of circle points */
                            /* space coordinates of the curve                  */
    const double  *x,const double  *y,const double  *z,       
    const double  *Color,    /* color coordinate of the curve (scaled to 0-1)   */
    const double  *Radius,   /* Radius scale of the curve                       */
    const double  *NX,const double  *NY,const double  *NZ,
    const double  *WX,const double  *WY,const double  *WZ
    )
/**************************************************************************/
{
    point_t *Ellipses[2],*EllipseNormals[2];
    int swap=0;

    int i,j;

    Ellipses[0]  = malloc( NCircle*sizeof(point_t) );
    Ellipses[1]  = malloc( NCircle*sizeof(point_t) );

    EllipseNormals[0] = malloc( NCircle*sizeof(point_t) );
    EllipseNormals[1] = malloc( NCircle*sizeof(point_t) );

    if ( !( Ellipses[0] && Ellipses[1] && EllipseNormals[0] && EllipseNormals[0] ) )
    {
        fprintf( stderr, "gomp_DrawTube: Can't alloc memory...\n" );

        if ( Ellipses[0] ) free( Ellipses[0] );
        if ( Ellipses[1] ) free( Ellipses[1] );

        if ( EllipseNormals[0] ) free( EllipseNormals[0] );
        if ( EllipseNormals[1] ) free( EllipseNormals[1] );

        return 1;
    } 

    ComputeTubeCircle( &x[0],&y[0],&z[0], &NX[0], &NY[0], &NZ[0],
                       &WX[0], &WY[0], &WZ[0], &Radius[0],
                       Ellipses[swap],EllipseNormals[swap],NCircle);

    for( i=1; i<NCurve; ++i )
    {
        ComputeTubeCircle( &x[i],&y[i],&z[i], &NX[i], &NY[i], &NZ[i],
                           &WX[i], &WY[i], &WZ[i], &Radius[i],
                           Ellipses[1-swap],EllipseNormals[1-swap],NCircle);

        glBegin( GL_QUAD_STRIP );

        for( j=0; j<NCircle; ++j )
        {
            glNormal3dv( &EllipseNormals[swap][j].x );
            if ( Color ) glTexCoord1d( Color[i-1] );
            glVertex3dv( &Ellipses[swap][j].x );

            glNormal3dv( &EllipseNormals[1-swap][j].x );
            if ( Color ) glTexCoord1d( Color[i] );
            glVertex3dv( &Ellipses[1-swap][j].x );
        }

        glEnd();

        swap=1-swap;
    }

    free( Ellipses[0] );
    free( Ellipses[1] );

    free( EllipseNormals[0] );
    free( EllipseNormals[1] );

    return 0;
}

/*
 *
 */
/**************************************************************************/
int gomp_DrawEllipseTube
( 
    int NCurve,int NCircle, /* number of curve points, number of circle points */
                            /* space coordinates of the curve                  */
    const double  *x,const double  *y,const double  *z,       
    const double  *Color,    /* color coordinate of the curve (scaled to 0-1)   */
    const double  *Width ,   /* Width scale of the curve (major axis)           */
    const double  *Thickness,/* Thickness scale of the curve (minor axis)       */
                            /* Direction of the minor axis                     */
    const double  *NX,const double  *NY,const double  *NZ,
                            /* Direction of the major axis                     */
    const double  *WX,const double  *WY,const double  *WZ
    )
/**************************************************************************/
{
    point_t *Ellipses[2],*EllipseNormals[2];
    int swap=0;

    int i,j;

    Ellipses[0]  = malloc( NCircle*sizeof(point_t) );
    Ellipses[1]  = malloc( NCircle*sizeof(point_t) );

    EllipseNormals[0] = malloc( NCircle*sizeof(point_t) );
    EllipseNormals[1] = malloc( NCircle*sizeof(point_t) );

    if ( !( Ellipses[0] && Ellipses[1] && EllipseNormals[0] && EllipseNormals[0] ) )
    {
        fprintf( stderr, "gomp_DrawEllipseTube: Can't alloc memory...\n" );

        if ( Ellipses[0] ) free( Ellipses[0] );
        if ( Ellipses[1] ) free( Ellipses[1] );

        if ( EllipseNormals[0] ) free( EllipseNormals[0] );
        if ( EllipseNormals[1] ) free( EllipseNormals[1] );

        return 1;
    } 

    ComputeTubeEllipse( &x[0],&y[0],&z[0], &NX[0], &NY[0], &NZ[0],
                        &WX[0], &WY[0], &WZ[0], &Width[0],&Thickness[0],
                        Ellipses[swap],EllipseNormals[swap],NCircle);

    for( i=1; i<NCurve; ++i )
    {
        ComputeTubeEllipse( &x[i],&y[i],&z[i], &NX[i], &NY[i], &NZ[i],
                            &WX[i], &WY[i], &WZ[i], &Width[i],&Thickness[i],
                            Ellipses[1-swap],EllipseNormals[1-swap],NCircle);

        glBegin( GL_QUAD_STRIP );

        for( j=0; j<NCircle; ++j )
        {
            glNormal3dv( &EllipseNormals[swap][j].x );
            if ( Color ) glTexCoord1d( Color[i-1] );
            glVertex3dv( &Ellipses[swap][j].x );

            glNormal3dv( &EllipseNormals[1-swap][j].x );
            if ( Color ) glTexCoord1d( Color[i] );
            glVertex3dv( &Ellipses[1-swap][j].x );
        }

        glEnd();

        swap=1-swap;
    }

    free( Ellipses[0] );
    free( Ellipses[1] );

    free( EllipseNormals[0] );
    free( EllipseNormals[1] );

    return 0;
}

/*
 *
 */
/**************************************************************************/
int gomp_DrawRectangleTube
( 
    int NCurve,              /* number of curve points                          */
                             /* space coordinates of the curve                  */
    const double  *x,const double  *y,const double  *z,
    const double  *Color,     /* color coordinate of the curve (scaled to 0-1)   */
    const double  *Width ,    /* Width scale of the curve (major axis)           */
    const double  *Thickness, /* Thickness scale of the curve (minor axis)       */
                             /* Direction of the minor axis                     */
    const double  *NX,const double  *NY,const double  *NZ,
                             /* Direction of the major axis                     */
    const double  *WX,const double  *WY,const double  *WZ
    )
/**************************************************************************/
{
    DrawRectangleTubePlane(
        NCurve,x,y,z,
        Color,Width,Thickness,
        NX,NY,NZ,1,
        WX,WY,WZ,1);

    DrawRectangleTubePlane(
        NCurve,x,y,z,
        Color,Thickness,Width,
        WX,WY,WZ,1,
        NX,NY,NZ,-1);

    DrawRectangleTubePlane(
        NCurve,x,y,z,
        Color,Width,Thickness,
        NX,NY,NZ,-1,
        WX,WY,WZ,-1);

    DrawRectangleTubePlane(
        NCurve,x,y,z,
        Color,Thickness,Width,
        WX,WY,WZ,-1,
        NX,NY,NZ,1);

    return 0;
}

/**************************************************************************/
int gomp_DrawStructureTrace
( 
    int NCurve,         /* number of curve points                         */
                        /* space coordinates of the curve                 */
    const double  *x,const double  *y,const double  *z,
    const double  *Color /* color coordinate of the curve (scaled to 0-1)  */
    )
/**************************************************************************/
{
    int i;

    glBegin(GL_LINE_STRIP);
    for( i=0; i<NCurve; ++i )
    {
        if ( Color ) glTexCoord1d( Color[i-1] );
        glVertex3d( x[i], y[i], z[i] );
    }
    glEnd();

    return 0;
}
#endif /* ENABLE_GRAPHICS */
/**************************************************************************/
int gomp_PlotPlumber(void *userData,int Wstr,int drawFlags)
/**************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static size_t i, Sets;
    static float  CRh,CGh,CBh;
    static const  PlumberDataType *Plumbers;

    if ( ! ( drawFlags & gom_PlotComplexElements    ) &&
         ! ( drawFlags & gom_PlotSimplifiedElements ) )
        return(-1);

    Plumbers = gomp_GetPlumbers();
    Sets     = gomp_DataVectorGetSize(&Plumbers);

    for ( i = 0 ; i < Sets ; i++ ) {

        if ( Plumbers[i].Wstr != Wstr )
            continue;

        CRh = Plumbers[i].red;
        CGh = Plumbers[i].green;
        CBh = Plumbers[i].blue;

        if ( ! gomp_GetDisplayColourType() )
            (void)gomp_RGB2Grayscale(&CRh , &CGh , &CBh);

        glColor3f(CRh , CGh , CBh);

        if( drawFlags & gom_PlotComplexElements ) {
            glEnable(GL_LIGHTING);

            switch ( Plumbers[i].Type )
            {
            case RIBBON_TYPE:
            case FLAT_HELIX_TYPE:
                (void)gomp_DrawRibbon(
                    Plumbers[i].Points,
                    Plumbers[i].x,  Plumbers[i].y,  Plumbers[i].z,
                    (const double  *)NULL,
                    Plumbers[i].width,
                    Plumbers[i].NX, Plumbers[i].NY, Plumbers[i].NZ,
                    Plumbers[i].WX, Plumbers[i].WY, Plumbers[i].WZ);
                break;
            case CYLINDER_TYPE:
                (void)gomp_DrawTube(
                    Plumbers[i].Points,
                    16 , 
                    Plumbers[i].x,  Plumbers[i].y,  Plumbers[i].z,
                    (const double  *)NULL,
                    Plumbers[i].width,
                    Plumbers[i].NX, Plumbers[i].NY, Plumbers[i].NZ,
                    Plumbers[i].WX, Plumbers[i].WY, Plumbers[i].WZ);
                break;
            case SOLID_HELIX_TYPE:
            case STRAND_TYPE:
                (void)gomp_DrawEllipseTube(
                    Plumbers[i].Points,
                    16 , 
                    Plumbers[i].x,  Plumbers[i].y,  Plumbers[i].z,
                    (const double  *)NULL,
                    Plumbers[i].width,
                    Plumbers[i].thickness,
                    Plumbers[i].NX, Plumbers[i].NY, Plumbers[i].NZ,
                    Plumbers[i].WX, Plumbers[i].WY, Plumbers[i].WZ);
                break;
            case ARROW_TYPE:
                (void)gomp_DrawRectangleTube(
                    Plumbers[i].Points,
                    Plumbers[i].x,  Plumbers[i].y,  Plumbers[i].z,
                    (const double  *)NULL,
                    Plumbers[i].width,
                    Plumbers[i].thickness,
                    Plumbers[i].NX, Plumbers[i].NY, Plumbers[i].NZ,
                    Plumbers[i].WX, Plumbers[i].WY, Plumbers[i].WZ);
                break;
            case TRACE_TYPE:
                glDisable(GL_LIGHTING);
                (void)gomp_DrawStructureTrace(
                    Plumbers[i].Points,
                    Plumbers[i].x,  Plumbers[i].y,  Plumbers[i].z,
                    (const double  *)NULL);
            }
        }
        else if( drawFlags & gom_PlotSimplifiedElements ) {
            glDisable(GL_LIGHTING);
            (void)gomp_DrawStructureTrace(
                Plumbers[i].Points,
                Plumbers[i].x,  Plumbers[i].y,  Plumbers[i].z,
                (const double  *)NULL);
        }
    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}

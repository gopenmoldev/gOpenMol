/*
  Enhancements 2002 - 2005 by:
  Eero Häkkinen

  From: Juha Ruokolainen <jpr@csc.fi>

  Niin vielä tällaisen tiedon tietty tarvitset vielä, eli
  miten tuo laatikon kulmapisteiden järjestys on määritelty


  Eight node brick volume element

  Type code: 808

    7---------6
   /|        /|
  4---------5 |
  | |       | |
  | |       | |
  | 3-------|-2     w v
 w|/v       |/      |/
  0---------1       ---u
  u


  (u,v,w ovat lokaali koordinaatisto, eivät viittaa normaalin komponentteihin)

  Yksi sovellus joka ei ehkä ensimmäisenä tule mieleen on hakea mielivaltaisen
  suuntaisia (ja paikka ei ole rajoitettu määriteltyihin tasoihin) 2d leikkauksia
  geometriasta...

  t. Juha
*/

#include "maindefs.h"

#include <string.h>
#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>

#if defined(WIN32)
#include <windows.h>
#endif

#include <tcl.h>

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>
#endif /* ENABLE_GRAPHICS */

#include "colouring.h"
#include "contour.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "tclutils.h"

#include "stdafx.h"

#ifdef ENABLE_GRAPHICS
static const int ElmBrickFace[6][9] = 
{
    { 0,1,2,3, 8, 9,10,11,20 },
    { 4,5,6,7,16,17,18,19,21 },
    { 0,1,5,4, 8,13,16,12,22 },
    { 3,2,6,7,10,14,18,15,24 },
    { 0,3,7,4,11,15,19,12,25 },
    { 1,2,6,5, 9,14,17,13,23 }
};

static struct {
    GLdouble modelMatrix[16];
    GLdouble projMatrix[16];
    GLint    viewport[4];
} ViewMatrices;

#if defined(SINGLE)
typedef float real;
#else
typedef double real;
#endif

static int elm_8node_brick_isosurface(real  ,const real *,const real *,
                                      const real *,const real *,const real *,
                                      const real *,const real *,const real *,
                                      polygon_t *);
static int elm_8node_brick_isosurface1(real  ,const real *,const real *,
                                       const real *,const real *,const real *,
                                       const real *,const real *,const real *,
                                       polygon_t *);
static int elm_4node_tetra_isosurface(real  ,const real *,const real *,
                                      const real *,const real *,const real *, 
                                      const real *,const real *,const real *,
                                      polygon_t *);

static void PlotJPRMeshSurfaceTriangles(size_t, polygon_t *, float);
static void PlotJPRSolidSurfaceTriangles(size_t, polygon_t *, float);
static void PlotJPRTextureMeshSurfaceTriangles(size_t, polygon_t *, float);
static void PlotJPRRainbowMeshSurfaceTriangles(size_t, polygon_t *, float);
static void PlotJPRTextureSolidSurfaceTriangles(size_t, polygon_t *, float);
static void PlotJPRRainbowSolidSurfaceTriangles(size_t, polygon_t *, float);

static GLdouble GetViewCoordinateZ(double x,double y,double z);
static int      SortInZOrder(const void *a,const void *b);


/******************************************************************************
 *
 *     Name:        elm_8node_brick_isosurface
 *
 *     Purpose:     Extract isosurfaces for brick element.
 *
 *     Parameters:
 *
 *         Input:   (double )  K: contour threshold
 *                  (const double  *) F: contour quantity values at nodes
 *                  (const double  *) C: color quantity values at nodes
 *                  (const double  *) X,Y,Z: node coordinates
 *                  (const double  *) U,V,W: normal vector at nodes
 *
 *         Output:  (polygon_t *)Polygon: output triangles.
 *
 *   Return value:  How many triangles we've got (possible values are 0-48)...
 *
 *****************************************************************************/
/****************************************************************************/
int elm_8node_brick_isosurface(real K,const real *F,const real *C,
                               const real *X,const real *Y,const real *Z,
                               const real *U,const real *V,const real *W,
                               polygon_t *Polygon)
/****************************************************************************/
{
    real tx[4],ty[4],tz[4],tu[4],tv[4],tw[4],tf[4],tc[4];
    int i,j,n;

    int above = 0, below = 0;

    for( i=0; i<8; i++ ) above += F[i]>K;
    for( i=0; i<8; i++ ) below += F[i]<K;
    if ( below == 8 || above == 8 ) return 0;

    tx[0] = 0.125*(X[0]+X[1]+X[2]+X[3]+X[4]+X[5]+X[6]+X[7]);
    ty[0] = 0.125*(Y[0]+Y[1]+Y[2]+Y[3]+Y[4]+Y[5]+Y[6]+Y[7]);
    tz[0] = 0.125*(Z[0]+Z[1]+Z[2]+Z[3]+Z[4]+Z[5]+Z[6]+Z[7]);
    tu[0] = 0.125*(U[0]+U[1]+U[2]+U[3]+U[4]+U[5]+U[6]+U[7]);
    tv[0] = 0.125*(V[0]+V[1]+V[2]+V[3]+V[4]+V[5]+V[6]+V[7]);
    tw[0] = 0.125*(W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6]+W[7]);
    tf[0] = 0.125*(F[0]+F[1]+F[2]+F[3]+F[4]+F[5]+F[6]+F[7]);
    tc[0] = 0.125*(C[0]+C[1]+C[2]+C[3]+C[4]+C[5]+C[6]+C[7]);
    
    n = 0;
    for( i=0; i<6; i++ )
    {
        tx[1] = 0.0;
        ty[1] = 0.0;
        tz[1] = 0.0;
        tu[1] = 0.0;
        tv[1] = 0.0;
        tw[1] = 0.0;
        tf[1] = 0.0;
        tc[1] = 0.0;
        for( j=0; j<4; j++ )
        {
            tx[1] += X[ElmBrickFace[i][j]];
            ty[1] += Y[ElmBrickFace[i][j]];
            tz[1] += Z[ElmBrickFace[i][j]];
            tu[1] += U[ElmBrickFace[i][j]];
            tv[1] += V[ElmBrickFace[i][j]];
            tw[1] += W[ElmBrickFace[i][j]];
            tf[1] += F[ElmBrickFace[i][j]];
            tc[1] += C[ElmBrickFace[i][j]];
        }
        tx[1] /= 4.0;
        ty[1] /= 4.0;
        tz[1] /= 4.0;
        tu[1] /= 4.0;
        tv[1] /= 4.0;
        tw[1] /= 4.0;
        tf[1] /= 4.0;
        tc[1] /= 4.0;
        
        for( j=0; j<3; j++ )
        {
            tx[2] = X[ElmBrickFace[i][j]];
            ty[2] = Y[ElmBrickFace[i][j]];
            tz[2] = Z[ElmBrickFace[i][j]];
            tu[2] = U[ElmBrickFace[i][j]];
            tv[2] = V[ElmBrickFace[i][j]];
            tw[2] = W[ElmBrickFace[i][j]];
            tf[2] = F[ElmBrickFace[i][j]];
            tc[2] = C[ElmBrickFace[i][j]];

            tx[3] = X[ElmBrickFace[i][j+1]];
            ty[3] = Y[ElmBrickFace[i][j+1]];
            tz[3] = Z[ElmBrickFace[i][j+1]];
            tu[3] = U[ElmBrickFace[i][j+1]];
            tv[3] = V[ElmBrickFace[i][j+1]];
            tw[3] = W[ElmBrickFace[i][j+1]];
            tf[3] = F[ElmBrickFace[i][j+1]];
            tc[3] = C[ElmBrickFace[i][j+1]];
            n += elm_4node_tetra_isosurface(
                K,tf,tc,tx,ty,tz,tu,tv,tw,&Polygon[n] );
        }

        tx[2] = X[ElmBrickFace[i][3]];
        ty[2] = Y[ElmBrickFace[i][3]];
        tz[2] = Z[ElmBrickFace[i][3]];
        tu[2] = U[ElmBrickFace[i][3]];
        tv[2] = V[ElmBrickFace[i][3]];
        tw[2] = W[ElmBrickFace[i][3]];
        tf[2] = F[ElmBrickFace[i][3]];
        tc[2] = C[ElmBrickFace[i][3]];

        tx[3] = X[ElmBrickFace[i][0]];
        ty[3] = Y[ElmBrickFace[i][0]];
        tz[3] = Z[ElmBrickFace[i][0]];
        tu[3] = U[ElmBrickFace[i][0]];
        tv[3] = V[ElmBrickFace[i][0]];
        tw[3] = W[ElmBrickFace[i][0]];
        tf[3] = F[ElmBrickFace[i][0]];
        tc[3] = C[ElmBrickFace[i][0]];
        n += elm_4node_tetra_isosurface(
            K,tf,tc,tx,ty,tz,tu,tv,tw,&Polygon[n] );
    }
    return n;
}


/******************************************************************************
 *
 *     Name:        elm_8node_brick_isosurface1
 *
 *     Purpose:     Extract isosurfaces for brick element.
 *
 *     Parameters:
 *
 *         Input:   (double )  K: contour threshold
 *                  (const double  *) F: contour quantity values at nodes
 *                  (const double  *) C: color quantity values at nodes
 *                  (const double  *) X,Y,Z: node coordinates
 *                  (const double  *) U,V,W: normal vector at nodes
 *
 *         Output:  (polygon_t *)Polygon: output triangles.
 *
 *   Return value:  How many triangles we've got (possible values are 0-24)...
 *
 *****************************************************************************/
int elm_8node_brick_isosurface1(real  K,const real *F,const real *C,
                                const real *X,const real *Y,const real *Z,
                                const real *U,const real *V,const real *W,
                                polygon_t *Polygon)
{
    real tx[4],ty[4],tz[4],tu[4],tv[4],tw[4],tf[4],tc[4];
    int i,j,l,n;

    static int map[12][3] =
        {
            { 0,1,2 }, { 0,2,3 }, { 4,5,6 }, { 4,6,7 }, { 3,2,6 }, { 3,6,7 },
            { 1,5,6 }, { 1,6,2 }, { 0,4,7 }, { 0,7,3 }, { 0,1,5 }, { 0,5,4 },
        };

    int above = 0, below = 0;

    for( i=0; i<8; i++ ) above += F[i]>K;
    for( i=0; i<8; i++ ) below += F[i]<K;
    if ( below == 8 || above == 8 ) return 0;

    tx[0] = 0.125*( X[0] + X[1] + X[2] + X[3] + X[4] + X[5] + X[6] + X[7]);
    ty[0] = 0.125*( Y[0] + Y[1] + Y[2] + Y[3] + Y[4] + Y[5] + Y[6] + Y[7]);
    tz[0] = 0.125*( Z[0] + Z[1] + Z[2] + Z[3] + Z[4] + Z[5] + Z[6] + Z[7]);
    tu[0] = 0.125*( U[0] + U[1] + U[2] + U[3] + U[4] + U[5] + U[6] + U[7]);
    tv[0] = 0.125*( V[0] + V[1] + V[2] + V[3] + V[4] + V[5] + V[6] + V[7]);
    tw[0] = 0.125*( W[0] + W[1] + W[2] + W[3] + W[4] + W[5] + W[6] + W[7]);
    tf[0] = 0.125*( F[0] + F[1] + F[2] + F[3] + F[4] + F[5] + F[6] + F[7]);
    tc[0] = 0.125*( C[0] + C[1] + C[2] + C[3] + C[4] + C[5] + C[6] + C[7]);

    n = 0;
    for( i=0; i<12; i++ )
    {
        for( j=1; j<4; j++ )
        {
            l = map[i][j-1];
            tx[j] = X[l];
            ty[j] = Y[l];
            tz[j] = Z[l];
            tu[j] = U[l];
            tv[j] = V[l];
            tw[j] = W[l];
            tf[j] = F[l];
            tc[j] = C[l];
        }
        n += elm_4node_tetra_isosurface(
            K,tf,tc,tx,ty,tz,tu,tv,tw,&Polygon[n] );
    }
    return n;
}

/******************************************************************************
 *
 *     Name:        elm_4node_tetra_isosurface
 *
 *     Purpose:     Extract isosurfaces for element
 *
 *     Parameters:
 *
 *         Input:   (double )  K: contour threshold
 *                  (const double  *) F: contour quantity values at nodes 
 *                  (const double  *) C: color quantity values at nodes 
 *                  (const double  *) X,Y,Z: node coordinates
 *                  (const double  *) U,V,W: normal vector at nodes
 *
 *         Output:  (polygon_t *)Polygon, output triangles (0,1 or 2) triangles
 *
 *   Return value:  How many triangles we've got...
 *
 *****************************************************************************/
/****************************************************************************/
int elm_4node_tetra_isosurface(real K,const real *F,const real *C,
                               const real *X,const real *Y,const real *Z,
                               const real *U,const real *V,const real *W,
                               polygon_t *Polygon)
/****************************************************************************/
{
    real t,tx[4],ty[4],tz[4],tu[4],tv[4],tw[4],tf[4],tc[4];
    real ax,ay,az,bx,by,bz,nx,ny,nz;
    int S0 = F[0] > K;
    int S1 = F[1] > K;
    int S2 = F[2] > K;
    int S3 = F[3] > K;

    int S = S0+S1+S2+S3,I[4],j;

    if ( S==0 || S==4 ) return 0;

    if ( S==1 || S==3 )
    {
        if ( (S==1 && S0) || (S==3 && !S0) )
        {
            I[0] = 0;
            I[1] = 1;
            I[2] = 2;
            I[3] = 3;
        } else if ( (S==1 && S1) || (S==3 && !S1) )
        {
            I[0] = 1;
            I[1] = 0;
            I[2] = 2;
            I[3] = 3;
        } else if ( (S==1 && S2) || (S==3 && !S2) )
        {
            I[0] = 2;
            I[1] = 0;
            I[2] = 1;
            I[3] = 3;
        } else if ( (S==1 && S3) || (S==3 && !S3) )
        {
            I[0] = 3;
            I[1] = 0;
            I[2] = 1;
            I[3] = 2;
        } else { return 0; }
   
        for( j=1; j<4; j++ )
        {
            t = (K-F[I[0]]) / (F[I[j]]-F[I[0]]);
            Polygon->x[j-1] = t*(X[I[j]]-X[I[0]]) + X[I[0]];
            Polygon->y[j-1] = t*(Y[I[j]]-Y[I[0]]) + Y[I[0]];
            Polygon->z[j-1] = t*(Z[I[j]]-Z[I[0]]) + Z[I[0]];
            Polygon->u[j-1] = t*(U[I[j]]-U[I[0]]) + U[I[0]];
            Polygon->v[j-1] = t*(V[I[j]]-V[I[0]]) + V[I[0]];
            Polygon->w[j-1] = t*(W[I[j]]-W[I[0]]) + W[I[0]];
            Polygon->c[j-1] = t*(C[I[j]]-C[I[0]]) + C[I[0]];
            Polygon->f[j-1] = K;
        }

        ax = Polygon->x[1] - Polygon->x[0];
        ay = Polygon->y[1] - Polygon->y[0];
        az = Polygon->z[1] - Polygon->z[0];

        bx = Polygon->x[2] - Polygon->x[0];
        by = Polygon->y[2] - Polygon->y[0];
        bz = Polygon->z[2] - Polygon->z[0];

        nx = ay*bz - az*by;
        ny = az*bx - ax*bz;
        nz = ax*by - ay*bx;

        ax = Polygon->u[0] + Polygon->u[1] + Polygon->u[2];
        ay = Polygon->v[0] + Polygon->v[1] + Polygon->v[2];
        az = Polygon->w[0] + Polygon->w[1] + Polygon->w[2];

        if ( nx*ax + ny*ay + nz*az < 0.0 )
        {
            real s;

#define swap( x,y ) { s=x; x=y; y=s; }

            swap( Polygon->x[1], Polygon->x[2] );
            swap( Polygon->y[1], Polygon->y[2] );
            swap( Polygon->z[1], Polygon->z[2] );
            swap( Polygon->u[1], Polygon->u[2] );
            swap( Polygon->v[1], Polygon->v[2] );
            swap( Polygon->w[1], Polygon->w[2] );
            swap( Polygon->f[1], Polygon->f[2] );
            swap( Polygon->c[1], Polygon->c[2] );

#undef swap
        }

        return 1;
    } else
    {
        if ( (S0 && S1) || (!S0 && !S1) )
        {
            t = (K-F[0])/ (F[2]-F[0]);
            tx[0] = t*(X[2]-X[0]) + X[0];
            ty[0] = t*(Y[2]-Y[0]) + Y[0];
            tz[0] = t*(Z[2]-Z[0]) + Z[0];
            tu[0] = t*(U[2]-U[0]) + U[0];
            tv[0] = t*(V[2]-V[0]) + V[0];
            tw[0] = t*(W[2]-W[0]) + W[0];
            tc[0] = t*(C[2]-C[0]) + C[0];
            tf[0] = K;

            t = (K-F[1]) / (F[2]-F[1]);
            tx[1] = t*(X[2]-X[1]) + X[1];
            ty[1] = t*(Y[2]-Y[1]) + Y[1];
            tz[1] = t*(Z[2]-Z[1]) + Z[1];
            tu[1] = t*(U[2]-U[1]) + U[1];
            tv[1] = t*(V[2]-V[1]) + V[1];
            tw[1] = t*(W[2]-W[1]) + W[1];
            tc[1] = t*(C[2]-C[1]) + C[1];
            tf[1] = K;

            t = (K-F[1]) / (F[3]-F[1]);
            tx[2] = t*(X[3]-X[1]) + X[1];
            ty[2] = t*(Y[3]-Y[1]) + Y[1];
            tz[2] = t*(Z[3]-Z[1]) + Z[1];
            tu[2] = t*(U[3]-U[1]) + U[1];
            tv[2] = t*(V[3]-V[1]) + V[1];
            tw[2] = t*(W[3]-W[1]) + W[1];
            tc[2] = t*(C[3]-C[1]) + C[1];
            tf[2] = K;

            t = (K-F[0]) / (F[3]-F[0]);
            tx[3] = t*(X[3]-X[0]) + X[0];
            ty[3] = t*(Y[3]-Y[0]) + Y[0];
            tz[3] = t*(Z[3]-Z[0]) + Z[0];
            tu[3] = t*(U[3]-U[0]) + U[0];
            tv[3] = t*(V[3]-V[0]) + V[0];
            tw[3] = t*(W[3]-W[0]) + W[0];
            tc[3] = t*(C[3]-C[0]) + C[0];
            tf[3] = K;
        }
        else if ( (S0 && S2) || (!S0 && !S2) )
        {
            t = (K-F[0]) / (F[1]-F[0]);
            tx[0] = t*(X[1]-X[0]) + X[0];
            ty[0] = t*(Y[1]-Y[0]) + Y[0];
            tz[0] = t*(Z[1]-Z[0]) + Z[0];
            tu[0] = t*(U[1]-U[0]) + U[0];
            tv[0] = t*(V[1]-V[0]) + V[0];
            tw[0] = t*(W[1]-W[0]) + W[0];
            tc[0] = t*(C[1]-C[0]) + C[0];
            tf[0] = K;

            t = (K-F[2]) / (F[1]-F[2]);
            tx[1] = t*(X[1]-X[2]) + X[2];
            ty[1] = t*(Y[1]-Y[2]) + Y[2];
            tz[1] = t*(Z[1]-Z[2]) + Z[2];
            tu[1] = t*(U[1]-U[2]) + U[2];
            tv[1] = t*(V[1]-V[2]) + V[2];
            tw[1] = t*(W[1]-W[2]) + W[2];
            tc[1] = t*(C[1]-C[2]) + C[2];
            tf[1] = K;

            t = (K-F[2]) / (F[3]-F[2]);
            tx[2] = t*(X[3]-X[2]) + X[2];
            ty[2] = t*(Y[3]-Y[2]) + Y[2];
            tz[2] = t*(Z[3]-Z[2]) + Z[2];
            tu[2] = t*(U[3]-U[2]) + U[2];
            tv[2] = t*(V[3]-V[2]) + V[2];
            tw[2] = t*(W[3]-W[2]) + W[2];
            tc[2] = t*(C[3]-C[2]) + C[2];
            tf[2] = K;

            t = (K-F[0]) / (F[3]-F[0]);
            tx[3] = t*(X[3]-X[0]) + X[0];
            ty[3] = t*(Y[3]-Y[0]) + Y[0];
            tz[3] = t*(Z[3]-Z[0]) + Z[0];
            tu[3] = t*(U[3]-U[0]) + U[0];
            tv[3] = t*(V[3]-V[0]) + V[0];
            tw[3] = t*(W[3]-W[0]) + W[0];
            tc[3] = t*(C[3]-C[0]) + C[0];
            tf[3] = K;
        }
        else if ( (S0 && S3) || (!S0 && !S3) )
        {
            t = (K-F[0]) / (F[1]-F[0]);
            tx[0] = t*(X[1]-X[0]) + X[0];
            ty[0] = t*(Y[1]-Y[0]) + Y[0];
            tz[0] = t*(Z[1]-Z[0]) + Z[0];
            tu[0] = t*(U[1]-U[0]) + U[0];
            tv[0] = t*(V[1]-V[0]) + V[0];
            tw[0] = t*(W[1]-W[0]) + W[0];
            tc[0] = t*(C[1]-C[0]) + C[0];
            tf[0] = K;

            t = (K-F[3]) / (F[1]-F[3]);
            tx[1] = t*(X[1]-X[3]) + X[3];
            ty[1] = t*(Y[1]-Y[3]) + Y[3];
            tz[1] = t*(Z[1]-Z[3]) + Z[3];
            tu[1] = t*(U[1]-U[3]) + U[3];
            tv[1] = t*(V[1]-V[3]) + V[3];
            tw[1] = t*(W[1]-W[3]) + W[3];
            tc[1] = t*(C[1]-C[3]) + C[3];
            tf[1] = K;

            t = (K-F[3]) / (F[2]-F[3]);
            tx[2] = t*(X[2]-X[3]) + X[3];
            ty[2] = t*(Y[2]-Y[3]) + Y[3];
            tz[2] = t*(Z[2]-Z[3]) + Z[3];
            tu[2] = t*(U[2]-U[3]) + U[3];
            tv[2] = t*(V[2]-V[3]) + V[3];
            tw[2] = t*(W[2]-W[3]) + W[3];
            tc[2] = t*(C[2]-C[3]) + C[3];
            tf[2] = K;

            t = (K-F[0]) / (F[2]-F[0]);
            tx[3] = t*(X[2]-X[0]) + X[0];
            ty[3] = t*(Y[2]-Y[0]) + Y[0];
            tz[3] = t*(Z[2]-Z[0]) + Z[0];
            tu[3] = t*(U[2]-U[0]) + U[0];
            tv[3] = t*(V[2]-V[0]) + V[0];
            tw[3] = t*(W[2]-W[0]) + W[0];
            tc[3] = t*(C[2]-C[0]) + C[0];
            tf[3] = K;
        }

        Polygon[0].x[0] = tx[0];
        Polygon[0].y[0] = ty[0];
        Polygon[0].z[0] = tz[0];
        Polygon[0].u[0] = tu[0];
        Polygon[0].v[0] = tv[0];
        Polygon[0].w[0] = tw[0];
        Polygon[0].f[0] = tf[0];
        Polygon[0].c[0] = tc[0];

        ax = tx[1] - tx[0];
        ay = ty[1] - ty[0];
        az = tz[1] - tz[0];

        bx = tx[2] - tx[0];
        by = ty[2] - ty[0];
        bz = tz[2] - tz[0];

        nx = ay*bz - az*by;
        ny = az*bx - ax*bz;
        nz = ax*by - ay*bx;

        ax = tu[0] + tu[1] + tu[2] + tu[3];
        ay = tv[0] + tv[1] + tv[2] + tv[3];
        az = tw[0] + tw[1] + tw[2] + tw[3];

        if ( nx*ax + ny*ay + nz*az >= 0.0 )
        {
            Polygon[0].x[1] = tx[1];
            Polygon[0].y[1] = ty[1];
            Polygon[0].z[1] = tz[1];
            Polygon[0].u[1] = tu[1];
            Polygon[0].v[1] = tv[1];
            Polygon[0].w[1] = tw[1];
            Polygon[0].f[1] = tf[1];
            Polygon[0].c[1] = tc[1];

            Polygon[0].x[2] = tx[2];
            Polygon[0].y[2] = ty[2];
            Polygon[0].z[2] = tz[2];
            Polygon[0].u[2] = tu[2];
            Polygon[0].v[2] = tv[2];
            Polygon[0].w[2] = tw[2];
            Polygon[0].f[2] = tf[2];
            Polygon[0].c[2] = tc[2];
        } else 
        {
            Polygon[0].x[1] = tx[2];
            Polygon[0].y[1] = ty[2];
            Polygon[0].z[1] = tz[2];
            Polygon[0].u[1] = tu[2];
            Polygon[0].v[1] = tv[2];
            Polygon[0].w[1] = tw[2];
            Polygon[0].f[1] = tf[2];
            Polygon[0].c[1] = tc[2];

            Polygon[0].x[2] = tx[1];
            Polygon[0].y[2] = ty[1];
            Polygon[0].z[2] = tz[1];
            Polygon[0].u[2] = tu[1];
            Polygon[0].v[2] = tv[1];
            Polygon[0].w[2] = tw[1];
            Polygon[0].f[2] = tf[1];
            Polygon[0].c[2] = tc[1];
        }

        Polygon[1].x[0] = tx[0];
        Polygon[1].y[0] = ty[0];
        Polygon[1].z[0] = tz[0];
        Polygon[1].u[0] = tu[0];
        Polygon[1].v[0] = tv[0];
        Polygon[1].w[0] = tw[0];
        Polygon[1].f[0] = tf[0];
        Polygon[1].c[0] = tc[0];

        ax = tx[2] - tx[0];
        ay = ty[2] - ty[0];
        az = tz[2] - tz[0];

        bx = tx[3] - tx[0];
        by = ty[3] - ty[0];
        bz = tz[3] - tz[0];

        nx = ay*bz - az*by;
        ny = az*bx - ax*bz;
        nz = ax*by - ay*bx;

        ax = tu[0] + tu[1] + tu[2] + tu[3];
        ay = tv[0] + tv[1] + tv[2] + tv[3];
        az = tw[0] + tw[1] + tw[2] + tw[3];

        if ( nx*ax + ny*ay + nz*az >= 0.0 )
        {
            Polygon[1].x[1] = tx[2];
            Polygon[1].y[1] = ty[2];
            Polygon[1].z[1] = tz[2];
            Polygon[1].u[1] = tu[2];
            Polygon[1].v[1] = tv[2];
            Polygon[1].w[1] = tw[2];
            Polygon[1].f[1] = tf[2];
            Polygon[1].c[1] = tc[2];

            Polygon[1].x[2] = tx[3];
            Polygon[1].y[2] = ty[3];
            Polygon[1].z[2] = tz[3];
            Polygon[1].u[2] = tu[3];
            Polygon[1].v[2] = tv[3];
            Polygon[1].w[2] = tw[3];
            Polygon[1].f[2] = tf[3];
            Polygon[1].c[2] = tc[3];
        } else
        {
            Polygon[1].x[1] = tx[3];
            Polygon[1].y[1] = ty[3];
            Polygon[1].z[1] = tz[3];
            Polygon[1].u[1] = tu[3];
            Polygon[1].v[1] = tv[3];
            Polygon[1].w[1] = tw[3];
            Polygon[1].f[1] = tf[3];
            Polygon[1].c[1] = tc[3];

            Polygon[1].x[2] = tx[2];
            Polygon[1].y[2] = ty[2];
            Polygon[1].z[2] = tz[2];
            Polygon[1].u[2] = tu[2];
            Polygon[1].v[2] = tv[2];
            Polygon[1].w[2] = tw[2];
            Polygon[1].f[2] = tf[2];
            Polygon[1].c[2] = tc[2];
        }

        return 2;
    }

    return 0;
}

#define GRID_ACCEPT 1.e-4

#if defined(NORMALIZE)
static inline void NormalizeJRPTrianleUVW(polygon_t *Triangle)
{
    float Sum;
    int i;

    for ( i = 0 ; i < 3 ; i++ ) {
        Sum   = 1.0 / sqrt(
            Triangle->u[i] * Triangle->u[i] +
            Triangle->v[i] * Triangle->v[i] +
            Triangle->w[i] * Triangle->w[i]);

        Triangle->u[i] *= Sum;
        Triangle->v[i] *= Sum;
        Triangle->w[i] *= Sum;
    }
}
#define NORMALIZE_JRP_TRIANGLE_UVW(Triangle) \
    NormalizeJRPTrianleUVW(&(Triangle))
#else
#define NORMALIZE_JRP_TRIANGLE_UVW(pTriangle)
#endif

static void PlotJPRMeshSurfaceTriangles(size_t TriangleCount,
                                        polygon_t *Triangles,
                                        float AlphaBlend)
{
    size_t i;
    for (i = 0; i < TriangleCount; i++) {
        glBegin(GL_LINE_LOOP);
        glVertex3f(Triangles[i].x[0] , Triangles[i].y[0] , Triangles[i].z[0]);
        glVertex3f(Triangles[i].x[1] , Triangles[i].y[1] , Triangles[i].z[1]);
        glVertex3f(Triangles[i].x[2] , Triangles[i].y[2] , Triangles[i].z[2]);
        glEnd();
    }
}
static void PlotJPRSolidSurfaceTriangles(size_t TriangleCount,
                                         polygon_t *Triangles,
                                         float AlphaBlend)
{
    size_t i;
    glBegin(GL_TRIANGLES);
    for (i = 0; i < TriangleCount; i++) {
        NORMALIZE_JRP_TRIANGLE_UVW(Triangles[i]);

        glNormal3f(Triangles[i].u[0] , Triangles[i].v[0] , Triangles[i].w[0]);
        glVertex3f(Triangles[i].x[0] , Triangles[i].y[0] , Triangles[i].z[0]);

        glNormal3f(Triangles[i].u[1] , Triangles[i].v[1] , Triangles[i].w[1]);
        glVertex3f(Triangles[i].x[1] , Triangles[i].y[1] , Triangles[i].z[1]);

        glNormal3f(Triangles[i].u[2] , Triangles[i].v[2] , Triangles[i].w[2]);
        glVertex3f(Triangles[i].x[2] , Triangles[i].y[2] , Triangles[i].z[2]);
    }
    glEnd();
}
static void PlotJPRTextureMeshSurfaceTriangles(size_t TriangleCount,
                                               polygon_t *Triangles,
                                               float AlphaBlend)
{
    size_t i;
    for (i = 0; i < TriangleCount; i++) {
        glBegin(GL_LINE_LOOP);

        glTexCoord1f(Triangles[i].c[0]);
        glVertex3f(Triangles[i].x[0] , Triangles[i].y[0] , Triangles[i].z[0]);

        glTexCoord1f(Triangles[i].c[1]);
        glVertex3f(Triangles[i].x[1] , Triangles[i].y[1] , Triangles[i].z[1]);

        glTexCoord1f(Triangles[i].c[2]);
        glVertex3f(Triangles[i].x[2] , Triangles[i].y[2] , Triangles[i].z[2]);

        glEnd();
    }
}
static void PlotJPRRainbowMeshSurfaceTriangles(size_t TriangleCount,
                                               polygon_t *Triangles,
                                               float AlphaBlend)
{
    size_t i;
    float  Red, Green, Blue;
    for (i = 0; i < TriangleCount; i++) {
        glBegin(GL_LINE_LOOP);

        gomp_PreRainbow(Triangles[i].c[0] , &Red , &Green , &Blue);
        glColor3f(Red , Green , Blue);
        glVertex3f(Triangles[i].x[0] , Triangles[i].y[0] , Triangles[i].z[0]);

        gomp_PreRainbow(Triangles[i].c[1] , &Red , &Green , &Blue);
        glColor3f(Red , Green , Blue);
        glVertex3f(Triangles[i].x[1] , Triangles[i].y[1] , Triangles[i].z[1]);

        gomp_PreRainbow(Triangles[i].c[2] , &Red , &Green , &Blue);
        glColor3f(Red , Green , Blue);
        glVertex3f(Triangles[i].x[2] , Triangles[i].y[2] , Triangles[i].z[2]);

        glEnd();
    }
}
static void PlotJPRTextureSolidSurfaceTriangles(size_t TriangleCount,
                                                polygon_t *Triangles,
                                                float AlphaBlend)
{
    size_t i;
    glBegin(GL_TRIANGLES);
    for (i = 0; i < TriangleCount; i++) {
        NORMALIZE_JRP_TRIANGLE_UVW(Triangles[i]);

        glTexCoord1f(Triangles[i].c[0]);
        glNormal3f(Triangles[i].u[0] , Triangles[i].v[0] , Triangles[i].w[0]);
        glVertex3f(Triangles[i].x[0] , Triangles[i].y[0] , Triangles[i].z[0]);

        glTexCoord1f(Triangles[i].c[1]);
        glNormal3f(Triangles[i].u[1] , Triangles[i].v[1] , Triangles[i].w[1]);
        glVertex3f(Triangles[i].x[1] , Triangles[i].y[1] , Triangles[i].z[1]);

        glTexCoord1f(Triangles[i].c[2]);
        glNormal3f(Triangles[i].u[2] , Triangles[i].v[2] , Triangles[i].w[2]);
        glVertex3f(Triangles[i].x[2] , Triangles[i].y[2] , Triangles[i].z[2]);
    }
    glEnd();
}
static void PlotJPRRainbowSolidSurfaceTriangles(size_t TriangleCount,
                                                polygon_t *Triangles,
                                                float AlphaBlend)
{
    size_t i;
    float  Red, Green, Blue;
    glBegin(GL_TRIANGLES);
    for (i = 0; i < TriangleCount; i++) {
        NORMALIZE_JRP_TRIANGLE_UVW(Triangles[i]);

        gomp_PreRainbow(Triangles[i].c[0] , &Red , &Green , &Blue);
        glColor4f(Red , Green , Blue , AlphaBlend);
        glNormal3f(Triangles[i].u[0] , Triangles[i].v[0] , Triangles[i].w[0]);
        glVertex3f(Triangles[i].x[0] , Triangles[i].y[0] , Triangles[i].z[0]);

        gomp_PreRainbow(Triangles[i].c[1] , &Red , &Green , &Blue);
        glColor4f(Red , Green , Blue , AlphaBlend);
        glNormal3f(Triangles[i].u[1] , Triangles[i].v[1] , Triangles[i].w[1]);
        glVertex3f(Triangles[i].x[1] , Triangles[i].y[1] , Triangles[i].z[1]);

        gomp_PreRainbow(Triangles[i].c[2] , &Red , &Green , &Blue);
        glColor4f(Red , Green , Blue , AlphaBlend);
        glNormal3f(Triangles[i].u[2] , Triangles[i].v[2] , Triangles[i].w[2]);
        glVertex3f(Triangles[i].x[2] , Triangles[i].y[2] , Triangles[i].z[2]);
    }
    glEnd();
}

static GLdouble GetViewCoordinateZ(double x,double y,double z)
{
    GLdouble view_x,view_y,view_z;

    gluProject(x,y,z,
               ViewMatrices.modelMatrix,
               ViewMatrices.projMatrix,
               ViewMatrices.viewport,
               &view_x,&view_y,&view_z);

    return view_z;
}

static int SortInZOrder(const void *a,const void *b)
{
    const polygon_t *ta = a;
    const polygon_t *tb = b;
    GLdouble za,zb,x,y;

    if ( gluProject(ta->x[0], ta->y[0], ta->z[0],
                    ViewMatrices.modelMatrix,
                    ViewMatrices.projMatrix,
                    ViewMatrices.viewport,
                    &x,&y,&za) != GL_FALSE &&
         gluProject(tb->x[0], tb->y[0], tb->z[0],
                    ViewMatrices.modelMatrix,
                    ViewMatrices.projMatrix,
                    ViewMatrices.viewport,
                    &x,&y,&zb) != GL_FALSE ) {
        if ( za < zb )
            return 1;
        else if ( za > zb )
            return -1;
    }
    return 0;
}

/****************************************************************************/
int gomp_CalculateJPRisosurface1
(
    int Which1  , int Which2    , float IsoValue , float AlphaBlend, 
    int ContSmooth, int ContMesh, float FScaleMin, float FScaleMax,
    int SortPolygons, polygon_t **polygon_vector
    )
/****************************************************************************/
{
    static int     PointsX;
    static int     PointsY;
    static int     PointsZ;
    static int     PointsXY;

    static struct {
        int dir;   /* 1       or -1            */
        int first; /* first index              */
        int end;   /* first index out of range */
        int curr;
        int neg_step[2];
        int pos_step[2];
    } i,j,k;

    static int       Ascending_i,l;
    static int       Point[8];
    static real      K;
    static real      F[8];
    static real      C[8];
    static real      X[8];
    static real      Y[8];
    static real      Z[8];
    static real      U[8];
    static real      V[8];
    static real      W[8];
    static polygon_t Triangles[48];
    static int       TriangleCount;
    static float     Dx;
    static float     Dy;
    static float     Dz;
    static float     Diff;

    static int       JumpJK;
    static int       ColorMapping;
    static const float *sumxyz;
    void (*PlotFunction)(size_t,polygon_t*,float);

    if(Which1 <  0 ||
       Which1 >= gomp_GetContoursDefined()) {
        gomp_FormatERROR(
            "contour data index 1 (%d) is out of allowed range (1 , %d)",
            Which1+1,gomp_GetContoursDefined());
        return(1);
    }

    if(Which2 <  0 ||
       Which2 >= gomp_GetContoursDefined()) {
        gomp_FormatERROR(
            "contour data index 2 (%d) is out of allowed range (1 , %d)",
            Which2+1,gomp_GetContoursDefined());
        return(1);
    }

    if(IsoValue < ContourInfo[Which1].min ||
       IsoValue > ContourInfo[Which1].max) {
        gomp_FormatERROR(
            "isovalue (%f) is outside defined range (%f , %f)",
            IsoValue,ContourInfo[Which1].min , ContourInfo[Which1].max); 
        return(1);
    }

    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);
    (void)gomp_Prepare1DTexture();
    glDisable(GL_TEXTURE_1D);

    if ( ContMesh ) /* mesh type */
        glDisable(GL_LIGHTING);

    ColorMapping = gomp_GetColorMappingType();

    /* Get projection matrices. */
    glMatrixMode(GL_PROJECTION);
    glGetDoublev(GL_PROJECTION_MATRIX,ViewMatrices.projMatrix);
    glMatrixMode(GL_MODELVIEW);
    glGetDoublev(GL_MODELVIEW_MATRIX, ViewMatrices.modelMatrix);
    glGetIntegerv(GL_VIEWPORT,        ViewMatrices.viewport);

    if(Which1 != Which2) {
/* check that both grid data files share the same space */
        if(ContourInfo[Which1].xdim != ContourInfo[Which2].xdim) {
            gomp_PrintERROR(
                "Number of points in grid data files in the "
                "x-direction does not match");
            return(1);
        }
        if(ContourInfo[Which1].ydim != ContourInfo[Which2].ydim) {
            gomp_PrintERROR(
                "Number of points in grid data files in the "
                "y-direction does not match");
            return(1);
        }
        if(ContourInfo[Which1].zdim != ContourInfo[Which2].zdim) {
            gomp_PrintERROR(
                "Number of points in grid data files in the "
                "z-direction does not match");
            return(1);
        }

        if(fabs(ContourInfo[Which1].Xmin - ContourInfo[Which2].Xmin) >
           GRID_ACCEPT) {
            gomp_PrintERROR(
                "Min. value in x-direction in grid data files does not match");
            return(1);
        }
        if(fabs(ContourInfo[Which1].Ymin - ContourInfo[Which2].Ymin) >
           GRID_ACCEPT) {
            gomp_PrintERROR(
                "Min. value in y-direction in grid data files does not match");
            return(1);
        }
        if(fabs(ContourInfo[Which1].Zmin - ContourInfo[Which2].Zmin) >
           GRID_ACCEPT) {
            gomp_PrintERROR(
                "Min. value in z-direction in grid data files does not match");
            return(1);
        }
        Diff           = FScaleMax - FScaleMin;
    } else {
        C[0] = C[1] = C[2] = C[3] = C[4] = C[5] = C[6] = C[7] = 1.0f;
    }

    if(Which1 == Which2) {
        glDisable(GL_TEXTURE_1D);
        if ( ContMesh )
            PlotFunction = PlotJPRMeshSurfaceTriangles;
        else
            PlotFunction = PlotJPRSolidSurfaceTriangles;
    } else { /* Which1 != Which2 */
        if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
            glEnable(GL_TEXTURE_1D);
            if ( ContMesh )
                PlotFunction = PlotJPRTextureMeshSurfaceTriangles;
            else
                PlotFunction = PlotJPRTextureSolidSurfaceTriangles;
        }
        else {
            if ( ContMesh )
                PlotFunction = PlotJPRRainbowMeshSurfaceTriangles;
            else
                PlotFunction = PlotJPRRainbowSolidSurfaceTriangles;
        }
    }

    if ( ! polygon_vector || gomp_DataVectorGetSize(polygon_vector) <= 0 ) {

        /* Create the vector if it is uninitialized. */
        if ( polygon_vector &&
             ! gomp_DataVectorCreateOrAppend(polygon_vector,
                                             &gomp_PolygonTypeHandler,0) ) {
            gomp_PrintERROR("can't allocate 'polygon_t' space");
            return(1);
        }

        sumxyz    = gomp_GetTranslateArray();

        PointsX   = ContourInfo[Which1].xdim;
        PointsY   = ContourInfo[Which1].ydim;
        PointsZ   = ContourInfo[Which1].zdim;

        Dx   = (ContourInfo[Which1].Xmax - ContourInfo[Which1].Xmin) / 
            (float)(ContourInfo[Which1].xdim - 1);

        Dy   = (ContourInfo[Which1].Ymax - ContourInfo[Which1].Ymin) / 
            (float)(ContourInfo[Which1].ydim - 1);

        Dz   = (ContourInfo[Which1].Zmax - ContourInfo[Which1].Zmin) / 
            (float)(ContourInfo[Which1].zdim - 1);

        K        = IsoValue;
        PointsXY = PointsX * PointsY;

        /* Find the right direction. We want to draw from back to front. */
        /* These are right criteria if view is ortogonal.                */
        /* In perspective view, result is good enough.                   */
        i.dir = ( GetViewCoordinateZ(ContourInfo[Which1].Xmin,0,0) >
                  GetViewCoordinateZ(ContourInfo[Which1].Xmax,0,0) ) ? 1 : -1;
        j.dir = ( GetViewCoordinateZ(0,ContourInfo[Which1].Ymin,0) >
                  GetViewCoordinateZ(0,ContourInfo[Which1].Ymax,0) ) ? 1 : -1;
        k.dir = ( GetViewCoordinateZ(0,0,ContourInfo[Which1].Zmin) >
                  GetViewCoordinateZ(0,0,ContourInfo[Which1].Zmax) ) ? 1 : -1;

        /* i goes from 1 to PointsX-1 or             */
        /*        from PointsX-2 to 0                */
        /* That is because previous value is copied. */
        i.first = ( i.dir > 0 ) ? 1 : PointsX;
        i.end   = ( i.dir > 0 ) ? ( PointsX - 2 ) : -1;

        /* These are symmetrical.               */
        /* j belongs to range 0...PointsY-1     */
        /* k belongs to range 0...PointsZ-1     */
        /* End is the first index out of range. */
        j.first = ( j.dir > 0 ) ? 0 : ( PointsY - 2 );
        j.end   = ( j.dir > 0 ) ? ( PointsY - 1 ) : -1;

        k.first = ( k.dir > 0 ) ? 0 : ( PointsZ - 2 );
        k.end   = ( k.dir > 0 ) ? ( PointsZ - 1 ) : -1;

        /*    7---------6      (0,1,1)---------(1,1,1)                  */
        /*   /|        /|           /|        /|                        */
        /*  4---------5 |    (0,0,1)---------(1,0,1)                    */
        /*  | |       | |          | |       | |                        */
        /*  | |       | |          | |       | |       z       w        */
        /*  | 3-------|-2      (0,1|0)-------|-(1,1,0) | y     | v      */
        /*  |/        |/           |/        |/        |/      |/       */
        /*  0---------1      (0,0,0)---------(1,0,0)   +--- x  +--- u   */
        /*                                                              */
        /* edge numbers      relative coordinates      axis             */
        /*                                                              */
        /* We can always do a negative jump if relative coordinate is 1.*/
        /* We can always do a positive jump if relative coordinate is 0.*/
        i.neg_step[1] = i.pos_step[0] = 1;
        j.neg_step[1] = j.pos_step[0] = PointsX;
        k.neg_step[1] = k.pos_step[0] = PointsXY;

        Ascending_i = ( i.dir > 0 ) ? 1 : 0;

        for ( k.curr = k.first ; k.curr != k.end ; k.curr += k.dir ) {

            Z[0]   =  k.curr * Dz + ContourInfo[Which1].Zmin - sumxyz[2];
            Z[1]   =  Z[0];
            Z[2]   =  Z[0];
            Z[3]   =  Z[0];
            Z[4]   =  Z[0] + Dz;
            Z[5]   =  Z[4];
            Z[6]   =  Z[4];
            Z[7]   =  Z[4];

            if( k.curr == 0 )
                k.neg_step[0] = 0;
            else
                k.neg_step[0] = PointsXY;

            if( k.curr == (PointsZ - 2) )
                k.pos_step[1] = 0;
            else
                k.pos_step[1] = PointsXY;

            for ( j.curr = j.first ; j.curr != j.end ; j.curr += j.dir ) {
                static const int EdgesInXOrder[2][4] = {
                    {0,3,4,7}, /* edges on the left   */
                    {1,2,5,6}  /* edges on the right  */
                };
                /* These correspond the above edges.                    */
                /* Relative Y and Z coordinates are the same for both   */
                /* the left and the right edges. X coordinate differs.  */
                static const int RelativeYCoords[4] = {0,1,0,1};
                static const int RelativeZCoords[4] = {0,0,1,1};
                int RelativeIndeces[4];
                RelativeIndeces[0] = 0;
                RelativeIndeces[1] = PointsX;
                RelativeIndeces[2] = PointsXY;
                RelativeIndeces[3] = PointsXY + PointsX;

                Y[0]   =  j.curr * Dy + ContourInfo[Which1].Ymin - sumxyz[1];
                Y[1]   =  Y[0];
                Y[2]   =  Y[0] + Dy;
                Y[3]   =  Y[2];
                Y[4]   =  Y[0];
                Y[5]   =  Y[0];
                Y[6]   =  Y[2];
                Y[7]   =  Y[2];

                JumpJK =  k.curr * PointsXY + j.curr * PointsX;

                if( j.curr == 0 )
                    j.neg_step[0] = 0;
                else
                    j.neg_step[0] = PointsX;

                if( j.curr == (PointsY - 2) )
                    j.pos_step[1] = 0;
                else
                    j.pos_step[1] = PointsX;

                if ( PointsX <= 2 )
                    i.neg_step[0] = i.pos_step[0] = 0;
                else {
                    i.neg_step[0] = 1 - Ascending_i;
                    i.pos_step[1] = Ascending_i;
                }

                X[EdgesInXOrder[Ascending_i][0]] =
                    i.first * Dx + ContourInfo[Which1].Xmin - sumxyz[0];
                Point[EdgesInXOrder[Ascending_i][0]] =
                    JumpJK + i.first;
                for ( l = 0 ; l < 4 ; l++ ) {
                    int index    = EdgesInXOrder[Ascending_i][l];
                    X[index]     = X[EdgesInXOrder[Ascending_i][0]];
                    Point[index] = Point[EdgesInXOrder[Ascending_i][0]] +
                        RelativeIndeces[l];
                    /* surface normals */
                    U[index] =
                        ContourInfo[Which1].data[
                            Point[index] + i.pos_step[1] ] -
                        ContourInfo[Which1].data[
                            Point[index] - i.neg_step[1] ];
                    V[index] =
                        ContourInfo[Which1].data[
                            Point[index] + j.pos_step[RelativeYCoords[l]] ] -
                        ContourInfo[Which1].data[
                            Point[index] - j.neg_step[RelativeYCoords[l]] ];
                    W[index] =
                        ContourInfo[Which1].data[
                            Point[index] + k.pos_step[RelativeZCoords[l]] ] -
                        ContourInfo[Which1].data[
                            Point[index] - k.neg_step[RelativeZCoords[l]] ];
                    F[index] = ContourInfo[Which1].data[Point[index]];
                }

                for ( i.curr = i.first ; i.curr != i.end ; i.curr += i.dir ) {
                    for ( l = 0 ; l < 4 ; l++ ) {
                        int from = EdgesInXOrder[Ascending_i  ][l];
                        int to   = EdgesInXOrder[1-Ascending_i][l];

                        X[to]     = X[from];
                        U[to]     = U[from];
                        V[to]     = V[from];
                        W[to]     = W[from];
                        F[to]     = F[from];
                        Point[to] = Point[from];
                    }

                    if( i.curr == 0 )
                        i.neg_step[0] = 0;
                    else
                        i.neg_step[0] = 1;

                    if( i.curr == (PointsX - 2) )
                        i.pos_step[1] = 0;
                    else
                        i.pos_step[1] = 1;

                    X[EdgesInXOrder[Ascending_i][0]] =
                        ( i.curr + i.dir ) * Dx +
                        ContourInfo[Which1].Xmin - sumxyz[0];
                    Point[EdgesInXOrder[Ascending_i][0]] =
                        JumpJK + ( i.curr + i.dir );
                    for ( l = 0 ; l < 4 ; l++ ) {
                        int index    = EdgesInXOrder[Ascending_i][l];
                        X[index]     = X[EdgesInXOrder[Ascending_i][0]];
                        Point[index] = Point[EdgesInXOrder[Ascending_i][0]] +
                            RelativeIndeces[l];
                        /* surface normals */
                        U[index] =
                            ContourInfo[Which1].data[
                                Point[index] + i.pos_step[1] ] -
                            ContourInfo[Which1].data[
                                Point[index] - i.neg_step[1] ];
                        V[index] =
                            ContourInfo[Which1].data[
                                Point[index] + j.pos_step[RelativeYCoords[l]] ] -
                            ContourInfo[Which1].data[
                                Point[index] -
                                j.neg_step[RelativeYCoords[l]] ];
                        W[index] =
                            ContourInfo[Which1].data[
                                Point[index] +
                                k.pos_step[RelativeZCoords[l]] ] -
                            ContourInfo[Which1].data[
                                Point[index] -
                                k.neg_step[RelativeZCoords[l]] ];
                        F[index] = ContourInfo[Which1].data[Point[index]];
                    }

                    if(Which1 == Which2) {
                        if(ContSmooth)
                            TriangleCount = elm_8node_brick_isosurface(
                                K, F, C, X, Y, Z, U, V, W, Triangles);
                        else
                            TriangleCount = elm_8node_brick_isosurface1(
                                K, F, C, X, Y, Z, U, V, W, Triangles);
                    } else { /* Which1!=Which2 */
                        C[0] = (ContourInfo[Which2].data[Point[0]] -
                                FScaleMin)/Diff;
                        C[1] = (ContourInfo[Which2].data[Point[1]] -
                                FScaleMin)/Diff;
                        C[2] = (ContourInfo[Which2].data[Point[2]] -
                                FScaleMin)/Diff;
                        C[3] = (ContourInfo[Which2].data[Point[3]] -
                                FScaleMin)/Diff;
                        C[4] = (ContourInfo[Which2].data[Point[4]] -
                                FScaleMin)/Diff;
                        C[5] = (ContourInfo[Which2].data[Point[5]] -
                                FScaleMin)/Diff;
                        C[6] = (ContourInfo[Which2].data[Point[6]] -
                                FScaleMin)/Diff;
                        C[7] = (ContourInfo[Which2].data[Point[7]] -
                                FScaleMin)/Diff;

                        if(ContSmooth)
                            TriangleCount = elm_8node_brick_isosurface(
                                K, F, C, X, Y, Z, U, V, W, Triangles);
                        else 
                            TriangleCount = elm_8node_brick_isosurface1(
                                K, F, C, X, Y, Z, U, V, W, Triangles);
                    }

                    if (TriangleCount) {
                        if ( polygon_vector ) {
                            size_t VectorLength = gomp_DataVectorGetSize(
                                polygon_vector);
                            if ( ! gomp_DataVectorSetSize(
                                polygon_vector,
                                VectorLength + TriangleCount ) ) {
                                gomp_PrintERROR(
                                    "can't allocate 'polygon_t' space");
                                return(1);
                            }
                            memcpy(
                                &(*polygon_vector)[VectorLength],
                                Triangles,
                                TriangleCount * sizeof(polygon_t));
                        } else
                            PlotFunction(TriangleCount,Triangles,AlphaBlend);
                    }
                }
            }
        }
    }

    if ( polygon_vector ) {
        size_t     TriangleCount = gomp_DataVectorGetSize(polygon_vector);
        polygon_t *Triangles     = *polygon_vector;

        if ( SortPolygons )
            qsort(Triangles,TriangleCount,sizeof(*Triangles),SortInZOrder);

        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        PlotFunction(TriangleCount,Triangles,AlphaBlend);
    }

    glDisable( GL_TEXTURE_1D );
    glEnable(GL_LIGHTING);

    return(0);
}
#endif /* ENABLE_GRAPHICS */

/****************************************************************************/
int gomp_ShowNumberOfPolygons(int Contour,int Level)
/****************************************************************************/
{
    /* if direct method is no return 0 */
    if(!gomp_GetSurfaceMethod())
        return(0);
    else
        return(gomp_DataVectorGetSize(
            &ContourInfo[Contour].levels[Level].polygons));
}
/****************************************************************************/
int  gomp_ReturnSurfacePolygonData(int Contour, int Level, int Which)
/****************************************************************************/
{
    static const float *sumxyz;
    static polygon_t *Polygon;

    sumxyz      = gomp_GetTranslateArray();

    Polygon     = &ContourInfo[Contour].levels[Level].polygons[Which];

    gomp_SendTclObjReturn(gomp_CreateTclList(
                            /* position */  /* normal */
                            "%f %f %f       %f %f %f" /* point 1 */
                            "%f %f %f       %f %f %f" /* point 2 */
                            "%f %f %f       %f %f %f" /* point 3 */
                            /* colour */
                            "%f %f %f",
                            /* point 1 */
                            Polygon->x[0], Polygon->y[0], Polygon->z[0],
                            Polygon->u[0], Polygon->v[0], Polygon->w[0],
                            /* point 2 */
                            Polygon->x[1], Polygon->y[1], Polygon->z[1],
                            Polygon->u[1], Polygon->v[1], Polygon->w[1],
                            /* point 3 */
                            Polygon->x[2], Polygon->y[2], Polygon->z[2],
                            Polygon->u[2], Polygon->v[2], Polygon->w[2],
                            /* colour */
                            Polygon->c[0], Polygon->c[1], Polygon->c[2]));

    return(0);
}

#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int gomp_CalculateJPRisosurface2
(
    int Which1 , int Which2,
    float pA, float pB, float pC, float pD,
    float IsoValue , float AlphaBlend , 
    int ContSmooth, float ScaleMin , float ScaleMax
    )
/****************************************************************************/
{
    static int     PointsX,i;
    static int     PointsY,j;
    static int     PointsZ,k;
    static int     Point0,Point1,Point2,Point3,Point4,Point5,Point6,Point7;

    static real  K;
    static real  F[8];
    static real  C[8];
    static real  X[8];
    static real  Y[8];
    static real  Z[8];
    static real  U[8];
    static real  V[8];
    static real  W[8];
    static polygon_t  Polygon[48];
    static int     Triangles;
    static float   Dx;
    static float   Dy;
    static float   Dz;
    static float   Diff;
    static float   FScaleMin;
    static float   FScaleMax;

    static int    Loop;
    static int    JumpJK;
    static int    tri;
    static int    PXY;
    static int    ColorMapping;
    static const float *sumxyz;

    IsoValue = 0.0;

    if(Which1 <  0 ||
       Which1 >= gomp_GetContoursDefined()) {
        gomp_PrintERROR("contour data index 1 is out of range");
        return(1);
    }

    FScaleMin = ScaleMin;
    FScaleMax = ScaleMax;

    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);
    (void)gomp_Prepare1DTexture();
    glDisable(GL_TEXTURE_1D);

    if(gomp_GetContourDisplayTypeGlobal())     /* mesh type */
        glDisable(GL_LIGHTING);

    ColorMapping = gomp_GetColorMappingType();

    if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
        glEnable(GL_TEXTURE_1D);
        U[0] = pA;
        W[0] = pB;
        V[0] = pC;

        U[1] = pA;
        W[1] = pB;
        V[1] = pC;

        U[2] = pA;
        W[2] = pB;
        V[2] = pC;

        U[3] = pA;
        W[3] = pB;
        V[3] = pC;

        U[4] = pA;
        W[4] = pB;
        V[4] = pC;

        U[5] = pA;
        W[5] = pB;
        V[5] = pC;

        U[6] = pA;
        W[6] = pB;
        V[6] = pC;

        U[7] = pA;
        W[7] = pB;
        V[7] = pC;
    } else {
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);
    }

    Diff      = FScaleMax - FScaleMin;

    sumxyz    = gomp_GetTranslateArray();

    PointsX   = ContourInfo[Which1].xdim;
    PointsY   = ContourInfo[Which1].ydim;
    PointsZ   = ContourInfo[Which1].zdim;

    Dx   = (ContourInfo[Which1].Xmax - ContourInfo[Which1].Xmin) / 
        (float)(ContourInfo[Which1].xdim - 1);

    Dy   = (ContourInfo[Which1].Ymax - ContourInfo[Which1].Ymin) / 
        (float)(ContourInfo[Which1].ydim - 1);

    Dz   = (ContourInfo[Which1].Zmax - ContourInfo[Which1].Zmin) / 
        (float)(ContourInfo[Which1].zdim - 1);

    sumxyz    = gomp_GetTranslateArray();
    Loop      = 0;

    K   = IsoValue;
    PXY = PointsX * PointsY;

    for(k = 0 ; k < PointsZ - 1 ; k++)    {

        Z[1]   =  k * Dz + ContourInfo[Which1].Zmin - sumxyz[2];
        Z[0]   =  Z[1];
        Z[2]   =  Z[1] + Dz;
        Z[3]   =  Z[2];
        Z[5]   =  Z[1];
        Z[4]   =  Z[5];
        Z[6]   =  Z[5] + Dz;
        Z[7]   =  Z[6];

        for(j = 0 ; j < PointsY - 1 ; j++)    {

            Y[1]   =  j * Dy + ContourInfo[Which1].Ymin - sumxyz[1];
            Y[0]   =  Y[1];
            Y[2]   =  Y[1];
            Y[3]   =  Y[2];
            Y[5]   =  Y[1] + Dy;
            Y[4]   =  Y[5];
            Y[6]   =  Y[5];
            Y[7]   =  Y[6];

            JumpJK  =  k*PXY + j*PointsX;

            for( i = 0 ; i < PointsX - 1; i++) {

/* i,j,k */
                X[1]   =  i * Dx + ContourInfo[Which1].Xmin - sumxyz[0]; 
                Point1 =  JumpJK + i;
/* i+1,j,k */
                X[0]   =  X[1] + Dx;
                Point0 =  JumpJK + i + 1;
/* i,j,k+1 */
                X[2]   =  X[1];
                Point2 =  JumpJK + PXY + i;
/* i+1,j,k+1 */
                X[3]   =  X[0];
                Point3 =  JumpJK + PXY + (i + 1);
/* i,j+1,k */
                X[5]   =   i      * Dx + ContourInfo[Which1].Xmin - sumxyz[0];
                Point5 =  JumpJK + PointsX + i;
/* i+1,j+1,k */
                X[4]   =  X[5] + Dx;
                Point4 =  JumpJK +  PointsX + i + 1;
/* i,j+1,k+1 */
                X[6]   =  X[5];
                Point6 =  JumpJK + PXY + PointsX + i;
/* i+1,j+1,k+1 */
                X[7]   =  X[4];
                Point7 =  JumpJK + PXY + PointsX + (i + 1);

                F[0] =
                    pA * (X[0] + sumxyz[0]) +
                    pB * (Y[0] + sumxyz[1]) +
                    pC * (Z[0] + sumxyz[2]) +
                    pD;
                F[1] =
                    pA * (X[1] + sumxyz[0]) +
                    pB * (Y[1] + sumxyz[1]) +
                    pC * (Z[1] + sumxyz[2]) +
                    pD;
                F[2] =
                    pA * (X[2] + sumxyz[0]) +
                    pB * (Y[2] + sumxyz[1]) +
                    pC * (Z[2] + sumxyz[2]) +
                    pD;
                F[3] =
                    pA * (X[3] + sumxyz[0]) +
                    pB * (Y[3] + sumxyz[1]) +
                    pC * (Z[3] + sumxyz[2]) +
                    pD;
                F[4] =
                    pA * (X[4] + sumxyz[0]) +
                    pB * (Y[4] + sumxyz[1]) +
                    pC * (Z[4] + sumxyz[2]) +
                    pD;
                F[5] =
                    pA * (X[5] + sumxyz[0]) +
                    pB * (Y[5] + sumxyz[1]) +
                    pC * (Z[5] + sumxyz[2]) +
                    pD;
                F[6] =
                    pA * (X[6] + sumxyz[0]) +
                    pB * (Y[6] + sumxyz[1]) +
                    pC * (Z[6] + sumxyz[2]) +
                    pD;
                F[7] =
                    pA * (X[7] + sumxyz[0]) +
                    pB * (Y[7] + sumxyz[1]) +
                    pC * (Z[7] + sumxyz[2]) +
                    pD;

                C[0] = (ContourInfo[Which1].data[Point0] - FScaleMin)/Diff;
                C[1] = (ContourInfo[Which1].data[Point1] - FScaleMin)/Diff;
                C[2] = (ContourInfo[Which1].data[Point2] - FScaleMin)/Diff;
                C[3] = (ContourInfo[Which1].data[Point3] - FScaleMin)/Diff;
                C[4] = (ContourInfo[Which1].data[Point4] - FScaleMin)/Diff;
                C[5] = (ContourInfo[Which1].data[Point5] - FScaleMin)/Diff;
                C[6] = (ContourInfo[Which1].data[Point6] - FScaleMin)/Diff;
                C[7] = (ContourInfo[Which1].data[Point7] - FScaleMin)/Diff;

                if(ContSmooth)
                    Triangles = elm_8node_brick_isosurface(
                        K , F , C , X , Y , Z ,
                        U , W , V , Polygon);
                else 
                    Triangles = elm_8node_brick_isosurface1(
                        K , F , C , X , Y , Z , 
                        U , W , V , Polygon);

                if(gomp_GetContourDisplayTypeGlobal() && Triangles) {
                    /* mesh type */

                    Loop++;

                    for (tri = 0; tri < Triangles; tri++) {

                        if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
                            glBegin(GL_LINE_LOOP);
                            glTexCoord1f(Polygon[tri].c[0]);
                            glVertex3f(Polygon[tri].x[0] ,
                                       Polygon[tri].y[0] , Polygon[tri].z[0]);

                            glTexCoord1f(Polygon[tri].c[1]);
                            glVertex3f(Polygon[tri].x[1] ,
                                       Polygon[tri].y[1] , Polygon[tri].z[1]);

                            glTexCoord1f(Polygon[tri].c[2]);
                            glVertex3f(Polygon[tri].x[2] ,
                                       Polygon[tri].y[2] , Polygon[tri].z[2]);
                            glEnd();
                        } else {
                            float  Red, Green, Blue;
                            glBegin(GL_LINE_LOOP);
                            gomp_PreRainbow((double)Polygon[tri].c[0] ,
                                            &Red , &Green , &Blue);
                            glColor3f(Red , Green , Blue);
                            glVertex3f(Polygon[tri].x[0] ,
                                       Polygon[tri].y[0] , Polygon[tri].z[0]);

                            gomp_PreRainbow((double)Polygon[tri].c[1] ,
                                            &Red , &Green , &Blue);
                            glColor3f(Red , Green , Blue);
                            glVertex3f(Polygon[tri].x[1] ,
                                       Polygon[tri].y[1] , Polygon[tri].z[1]);

                            gomp_PreRainbow((double)Polygon[tri].c[2] ,
                                            &Red , &Green , &Blue);
                            glColor3f(Red , Green , Blue);
                            glVertex3f(Polygon[tri].x[2] ,
                                       Polygon[tri].y[2] , Polygon[tri].z[2]);
                            glEnd();
                        }
                    }
                } else if(Triangles) {

                    Loop++;

                    for (tri = 0; tri < Triangles; tri++) {
                        NORMALIZE_JRP_TRIANGLE_UVW(Polygon[tri]);

                        if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
                            glBegin(GL_POLYGON);
                            glTexCoord1f(Polygon[tri].c[0]);
                            glNormal3f(Polygon[tri].u[0] ,
                                       Polygon[tri].v[0] ,
                                       Polygon[tri].w[0]);
                            glVertex3f(Polygon[tri].x[0] ,
                                       Polygon[tri].y[0] ,
                                       Polygon[tri].z[0]);

                            glTexCoord1f(Polygon[tri].c[1]);
                            glNormal3f(Polygon[tri].u[1] ,
                                       Polygon[tri].v[1] ,
                                       Polygon[tri].w[1]);
                            glVertex3f(Polygon[tri].x[1] ,
                                       Polygon[tri].y[1] ,
                                       Polygon[tri].z[1]);


                            glTexCoord1f(Polygon[tri].c[2]);
                            glNormal3f(Polygon[tri].u[2] ,
                                       Polygon[tri].v[2] ,
                                       Polygon[tri].w[2]);
                            glVertex3f(Polygon[tri].x[2] ,
                                       Polygon[tri].y[2] ,
                                       Polygon[tri].z[2]);
                            glEnd();
                        } else {
                            float  Red, Green, Blue;
                            glBegin(GL_POLYGON);
                            gomp_PreRainbow((double)Polygon[tri].c[0] ,
                                            &Red , &Green , &Blue);
                            glColor4f(Red , Green , Blue , AlphaBlend);
                            glVertex3f(Polygon[tri].x[0] ,
                                       Polygon[tri].y[0] ,
                                       Polygon[tri].z[0]);

                            gomp_PreRainbow((double)Polygon[tri].c[1] ,
                                            &Red , &Green , &Blue);
                            glColor4f(Red , Green , Blue , AlphaBlend);
                            glVertex3f(Polygon[tri].x[1] ,
                                       Polygon[tri].y[1] ,
                                       Polygon[tri].z[1]);

                            gomp_PreRainbow((double)Polygon[tri].c[2] ,
                                            &Red , &Green , &Blue);
                            glColor4f(Red , Green , Blue , AlphaBlend);
                            glVertex3f(Polygon[tri].x[2] ,
                                       Polygon[tri].y[2] ,
                                       Polygon[tri].z[2]);
                            glEnd();
                        }
                    }
                }
            }
        }
    }
/*
  if(gomp_GetDisplayListState() && 
  (gomp_GetDisplayListRegenerateState() & DISPLAY_LIST_CUTPLANE)) {
  gomp_SetDisplayListRegenerateState(~DISPLAY_LIST_CUTPLANE);
  glEndList();
  glCallList(SURFACE2_DISPLAY_LIST);

  }
*/
    glDisable( GL_TEXTURE_1D );
    glEnable(GL_LIGHTING);

    return(0);
}
#endif /* ENABLE_GRAPHICS */

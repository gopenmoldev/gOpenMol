/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Original coding by Juha Ruokolainen.

Enhancements 2002 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <stdlib.h>

#include "plumber.h"

#include "stdafx.h"

static void SolveTriDiag( int , const double  *, const double  *, double  * );

/**************************************************************************
 *
 * Solve a tridiagonal system resulting from cubic spline interpolation.
 *
 **************************************************************************/
void SolveTriDiag( int n, const double  *y, const double  *h, double  *r )
{
    double s,*b;
    int i;

    b = malloc( n*sizeof(double) );

    for( i=1; i<n-1; i++ )
    {
        b[i] = 2 * ( h[i-1] + h[i] );
        r[i] = 3 * ( h[i]   * ( y[i]-y[i-1] ) / h[i-1] +
                     h[i-1] * ( y[i+1]-y[i] ) / h[i] );
    }

    r[1] -= h[1]*r[0];
    for( i=1; i<n-2; i++ )
    {
        s = -h[i+1] / b[i];
        r[i+1] += s*r[i];
        b[i+1] += s*h[i-1];
    }

    for( i=n-2; i>=1; i-- ) r[i] = (r[i] - h[i-1]*r[i+1]) / b[i];

    free(b);
}


/************************************************************************
 *
 * Compute cubic spline from given n knots. The interpolated values
 * are written to vector x, which also contains the knot values at
 * x[Knots[i]] on entry. The tangent directions at endpoints are
 * approximated by one sided differences.
 *
 ************************************************************************/
int gomp_CubicSpline( double  *x, int n, const int *Knots )
{
    double a,b,t,*y,*r,*h;
    int i,j,k;

    y = malloc( n*sizeof(double) );
    h = malloc( n*sizeof(double) );
    r = malloc( n*sizeof(double) );

    for( i=0; i<n; i++ )   y[i] = x[Knots[i]];
    for( i=0; i<n-1; i++ ) h[i] = Knots[i+1] - Knots[i];

    r[0]   = ( y[1] - y[0] ) / h[0];
    r[n-1] = ( y[n-1] - y[n-2] ) / h[n-2];

    SolveTriDiag( n,y,h,r );

    for( i=0,k=0; i<n-1; i++ )
    {
        a = -2 * ( y[i+1] - y[i] ) + (   r[i] + r[i+1] ) * h[i];
        b =  3 * ( y[i+1] - y[i] ) - ( 2*r[i] + r[i+1] ) * h[i];

        for( j=Knots[i]; j<Knots[i+1]; j++,k++ )
        {
            t = (double)(j - Knots[i]) / (double)(Knots[i+1] - Knots[i]);
            x[k] = ((a*t + b) * t + r[i]*h[i]) * t + y[i];
        }
    }

    free( r );
    free( h );
    free( y );

    return(0);
}

/************************************************************************
 *
 * Compute linear spline from given n knots. The interpolated values
 * are written to vector x, which also contains the knot values at
 * x[Knots[i]] on entry.
 *
 ************************************************************************/
int gomp_LinearSpline( double  *x, int n, const int *Knots )
{
    double t;
    int i,j,k;

    for( i=0,k=0; i<n-1; ++i )
    {
        for( j=Knots[i]; j<Knots[i+1]; ++j,++k )
        {
            t = (double)(j- Knots[i]) / (double)(Knots[i+1] - Knots[i]);
            x[k] = ( 1.0 - t ) * x[Knots[i]] + t * x[Knots[i+1]];
        }
    }

    return(0);
}

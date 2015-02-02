
/*
  The program to superimpose atoms of two molecules by quaternion method

  David J. Heisterberg
  The Ohio Supercomputer Center
  1224 Kinnear Rd.
  Columbus, OH  43212-1163
  (614)292-6036
  djh@osc.edu    djh@ohstpy.bitnet    ohstpy::djh

  Translated to C from fitest.f program and interfaced with Xmol program
  by Jan Labanowski,  jkl@osc.edu   jkl@ohstpy.bitnet   ohstpy::jkl

  To complile: 
  cc -o quatfit quatfit.c -lm

  Copyright: Ohio Supercomputer Center, David J. Heisterberg, 1990.
  The program can be copied and distributed freely, provided that
  this copyright in not removed. You may acknowledge the use of the
  program in published material as:
  David J. Heisterberg, 1990, unpublished results.


  Enhancements 2002 - 2004 by:
  Eero Häkkinen
*/

/* there was a typo in the formula for RMS. In the part calculating the
RMS, there was a line:
   wnorm += s;  
while it should be:
   wnorm += s*s;
Martin Lema mlema[-at-]unq.edu.ar was kind to find it. THANKS!!! jkl.
It did not affect fitting. But the RMS value was not right when weights
were different than 1.
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <ctype.h>
#include <string.h>

#include <stdlib.h>

#include <sys/types.h>

#include "gomstring.h"
#include "measure.h"
#include "memalloc.h"
#include "printmsg.h"

#include "stdafx.h"

#define MAXPOINTS     400
#define MAXLINELEN    250

#define ATOM_DISP_ON     1
#define ATOM_DISP_OFF    0

/* options
   -r refmol    reference molecule xmol file. If this option is not given, the
   information is read from standard input.

   -f fitmol    fitted molecule input xmol file If this option is not given, the
   information is read from standard input.

   -p pairs     file with the list of fitted atom pairs and weights. If this
   option is not specified, pairs and weights are taken from
   stdin. If file name "none" is used (i.e. -p none), atoms of the
   fitted molecule are fitted to atoms of the reference
   molecule with the same number, and all weights are assumed 1.0.
   If molecules do not have the same number of atoms, the
   smaller number of atoms is fitted. 

   -o outmol    output file for transformed fitted molecule. If this option
   is not given, the information is written to standard output.

   -s statfile  file with fit statistics (distances between fitted atoms, etc).
   If this option is not given, the information is written to
   standard output.

   If any files are read from stdin, the order is: refmol, fitmol, pairs.
   If any files are written to stdout, the order is: outmol, statfile.

   The file formats are:
   The refmol, fitmol and outmol files are in the XYZ format used with xmol
   program from Minnesota Supercomputer Institute. The format is:
   1st line: number of atoms
   2nd line: title
   3rd and next lines have the format depending on the kind of information:
   T  X  Y  Z                (total of 4 columns)
   T  X  Y  Z  C             (total of 5 columns)
   T  X  Y  Z  Mx My Mz      (total of 7 columns)
   T  X  Y  Z  C Mx My Mz    (total of 8 columns)
   where T is atom type (usually elemnt symbol), X, Y, Z are cartesian
   coordinates of an atom in Angstroms, C is atomic charge, and Mx, My, Mz
   are normal modes.

   The pairs file format is:
   1st line: number of pairs
   2nd and next lines:
   Ar   Af    W
   where Ar is the atom number of the reference molecule, Af is the atom
   number of fitted molecule, w is the statistical weight. Weights W are
   related to the square of expected deviation "sigma" between the reference
   and fitted molecule atoms and allow to make fit of some atom pairs more
   tight. W is proportional to 1/sigma^2. The larger the weight, the more
   tight will be the resulting fit for the given pair.

   The statfile lists results of the fit with explanation.

*/   


/*====================================================================
  CENTER
  center or translate a molecule. 
  n - number of atoms
  x - on input  - original xyz coordinates of a molecule
  on output - moved xyz coordinates (see io for modes).

  w - if io=1, weights of atoms
  if io=2 or 3, unused

  io - 1 weighted geometric center of the molecule will be at (0,0,0)
  2 molecule will be moved by a vector -o (i.e., components of a vector o
  will be subtracted from atom coordinates). 
  3 molecule will be moved by a vector +o (i.e., components of a vector o
  will be added atom coordinates). 

  o - if io=1, output, center of original coordinates
  if io=2, input, vector o will be subtracted from atomic coordinates
  if io=3, input, vector o will be added to atomic coordinates

  =====================================================================*/
static void center(int n, double  *x[], const double  *w, int io, double o[4]);
static void rotmol(int n, double  *x[4], double  *y[4], double u[4][4]);
static void jacobi(double a[4][4], double d[4], double v[4][4], int nrot);

void center(int n, double  *x[], const double  *w, int io, double o[4])
{
    double wnorm, modif;
    int i;

    if (io == 2) {
        modif = -1.0;
    }
    else if (io == 3) {
        modif = 1.0;
    }
    else {
        modif = -1.0;
        o[1] = 0.0;
        o[2] = 0.0;
        o[3] = 0.0;
        wnorm = 0.0;
        for (i = 1; i <= n; i++) {
            o[1] = o[1] + x[1][i] * sqrt(w[i]);
            o[2] = o[2] + x[2][i] * sqrt(w[i]);
            o[3] = o[3] + x[3][i] * sqrt(w[i]);
            wnorm = wnorm + sqrt(w[i]);
        }
        o[1] = o[1] / wnorm;
        o[2] = o[2] / wnorm;
        o[3] = o[3] / wnorm;
    }


    for (i = 1; i <= n; i++) {
        x[1][i] = x[1][i] + modif*o[1];
        x[2][i] = x[2][i] + modif*o[2];
        x[3][i] = x[3][i] + modif*o[3];
    }

}

/*================================
  ROTMOL
  rotate a molecule
  n - number of atoms
  x - input coordinates
  y - rotated coordinates y = u * x
  u - left rotation matrix
  ==================================*/

void rotmol (int n, double  *x[4], double  *y[4], double u[4][4])
/*
  int n;
  const double  *x[4];
  const double  *y[4];
  double u[4][4];
*/
{
    double yx, yy, yz;
    int i;

    for (i = 1; i <= n; i++) {
        yx = u[1][1] * x[1][i] + u[1][2] * x[2][i] + u[1][3] * x[3][i];
        yy = u[2][1] * x[1][i] + u[2][2] * x[2][i] + u[2][3] * x[3][i];
        yz = u[3][1] * x[1][i] + u[3][2] * x[2][i] + u[3][3] * x[3][i];

        y[1][i] = yx;
        y[2][i] = yy;
        y[3][i] = yz;
    }
}

/*=======================================================
  JACOBI
  Jacobi diagonalizer with sorted output. It is only good for 4x4 matrices.
  (was too lazy to do pointers...)
  a - input: matrix to diagonalize
  v - output: eigenvectors
  d - output: eigenvalues
  nrot - input: maximum number of sweeps
  =========================================================*/

void jacobi(double a[4][4], double d[4], double v[4][4], int nrot)
/*
  int nrot;
  double a[4][4];
  double d[4];
  double v[4][4];
*/
{
    double onorm, dnorm;
    double b, dma, q, t, c, s;
    double atemp, vtemp, dtemp;
    int i, j, k, l;

    for (j = 0; j <= 3; j++) {
        for (i = 0; i <= 3; i++) {
            v[i][j] = 0.0;
        }
        v[j][j] = 1.0;
        d[j] = a[j][j];
    }

    for (l = 1; l <= nrot; l++) {
        dnorm = 0.0;
        onorm = 0.0;
        for (j = 0; j <= 3; j++) {
            dnorm = dnorm + fabs(d[j]);
            for (i = 0; i <= j - 1; i++) {
                onorm = onorm + fabs(a[i][j]);
            }
        }
        if((onorm/dnorm) <= 1.0e-12) goto Exit_now;
        for (j = 1; j <= 3; j++) {
            for (i = 0; i <= j - 1; i++) {
                b = a[i][j];
                if(fabs(b) > 0.0) {
                    dma = d[j] - d[i];
                    if((fabs(dma) + fabs(b)) <=  fabs(dma)) {
                        t = b / dma;
                    }
                    else {
                        q = 0.5 * dma / b;
                        t = 1.0/(fabs(q) + sqrt(1.0+q*q));
                        if(q < 0.0) {
                            t = -t;
                        }
                    }
                    c = 1.0/sqrt(t * t + 1.0);
                    s = t * c;
                    a[i][j] = 0.0;
                    for (k = 0; k <= i-1; k++) {
                        atemp = c * a[k][i] - s * a[k][j];
                        a[k][j] = s * a[k][i] + c * a[k][j];
                        a[k][i] = atemp;
                    }
                    for (k = i+1; k <= j-1; k++) {
                        atemp = c * a[i][k] - s * a[k][j];
                        a[k][j] = s * a[i][k] + c * a[k][j];
                        a[i][k] = atemp;
                    }
                    for (k = j+1; k <= 3; k++) {
                        atemp = c * a[i][k] - s * a[j][k];
                        a[j][k] = s * a[i][k] + c * a[j][k];
                        a[i][k] = atemp;
                    }
                    for (k = 0; k <= 3; k++) {
                        vtemp = c * v[k][i] - s * v[k][j];
                        v[k][j] = s * v[k][i] + c * v[k][j];
                        v[k][i] = vtemp;
                    }
                    dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b;
                    d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b;
                    d[i] = dtemp;
                }  /* end if */
            } /* end for i */
        } /* end for j */
    } /* end for l */
 
  Exit_now:

    nrot = l;

    for (j = 0; j <= 2; j++) {
        k = j;
        dtemp = d[k];
        for (i = j+1; i <= 3; i++) {
            if(d[i] < dtemp) {
                k = i;
                dtemp = d[k];
            }
        }

        if(k > j) {
            d[k] = d[j];
            d[j] = dtemp;
            for (i = 0; i <= 3; i++) {
                dtemp = v[i][k];
                v[i][k] = v[i][j];
                v[i][j] = dtemp;
            }
        }
    }
}



/*==========================================
  Q2MAT
  Generate a left rotation matrix from a normalized quaternion

  INPUT
  q      - normalized quaternion

  OUTPUT
  u      - the rotation matrix
  ===========================================*/

static void q2mat (double q[4], double u[4][4])
/*
  double q[4];
  double u[4][4];
*/
{
    u[1][1] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
    u[2][1] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
    u[3][1] = 2.0 * (q[1] * q[3] + q[0] * q[2]);

    u[1][2] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
    u[2][2] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
    u[3][2] = 2.0 * (q[2] * q[3] - q[0] * q[1]);

    u[1][3] = 2.0 *(q[3] * q[1] - q[0] * q[2]);
    u[2][3] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
    u[3][3] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}


/*==========================================
  QTRFIT
  Find the quaternion, q,[and left rotation matrix, u] that minimizes

  |qTXq - Y| ^ 2  [|uX - Y| ^ 2]

  This is equivalent to maximizing Re (qTXTqY).

  This is equivalent to finding the largest eigenvalue and corresponding
  eigenvector of the matrix

  [A2   AUx  AUy  AUz ]
  [AUx  Ux2  UxUy UzUx]
  [AUy  UxUy Uy2  UyUz]
  [AUz  UzUx UyUz Uz2 ]

  where

  A2   = Xx Yx + Xy Yy + Xz Yz
  Ux2  = Xx Yx - Xy Yy - Xz Yz
  Uy2  = Xy Yy - Xz Yz - Xx Yx
  Uz2  = Xz Yz - Xx Yx - Xy Yy
  AUx  = Xz Yy - Xy Yz
  AUy  = Xx Yz - Xz Yx
  AUz  = Xy Yx - Xx Yy
  UxUy = Xx Yy + Xy Yx
  UyUz = Xy Yz + Xz Yy
  UzUx = Xz Yx + Xx Yz

  The left rotation matrix, u, is obtained from q by

  u = qT1q

  INPUT
  n      - number of points
  x      - fitted molecule coordinates
  y      - reference molecule coordinates
  w      - weights

  OUTPUT
  q      - the best-fit quaternion
  u      - the best-fit left rotation matrix
  nr     - max number of gomp_jacobi sweeps

  =====================================*/

static void qtrfit (int n, double  *x[4], double  *y[4], const double  *w, double q[4], double u[4][4], int nr)
/*
  int n;
  const double  *x[4];
  const double  *y[4];
  const double  *w;
  double q[4];
  double u[4][4];
  int nr;
*/
{
    double xxyx, xxyy, xxyz;
    double xyyx, xyyy, xyyz;
    double xzyx, xzyy, xzyz;
    double c[4][4], v[4][4];
    double d[4];
    int i, j;


/* generate the upper triangle of the quadratic form matrix */

    xxyx = 0.0;
    xxyy = 0.0;
    xxyz = 0.0;
    xyyx = 0.0;
    xyyy = 0.0;
    xyyz = 0.0;
    xzyx = 0.0;
    xzyy = 0.0;
    xzyz = 0.0;
 
    for (i = 1; i <= n; i++) {
        xxyx = xxyx + x[1][i] * y[1][i] * w[i];
        xxyy = xxyy + x[1][i] * y[2][i] * w[i];
        xxyz = xxyz + x[1][i] * y[3][i] * w[i];
        xyyx = xyyx + x[2][i] * y[1][i] * w[i];
        xyyy = xyyy + x[2][i] * y[2][i] * w[i];
        xyyz = xyyz + x[2][i] * y[3][i] * w[i];
        xzyx = xzyx + x[3][i] * y[1][i] * w[i];
        xzyy = xzyy + x[3][i] * y[2][i] * w[i];
        xzyz = xzyz + x[3][i] * y[3][i] * w[i];
    }
 
    for(i = 0; i <= 3; i++) {
        for(j = 0; j <= 3; j++) {
            c[i][j] = 0.0;
        }
    }

    c[0][0] = xxyx + xyyy + xzyz;

    c[0][1] = xzyy - xyyz;
    c[1][1] = xxyx - xyyy - xzyz;

    c[0][2] = xxyz - xzyx;
    c[1][2] = xxyy + xyyx;
    c[2][2] = xyyy - xzyz - xxyx;

    c[0][3] = xyyx - xxyy;
    c[1][3] = xzyx + xxyz;
    c[2][3] = xyyz + xzyy;
    c[3][3] = xzyz - xxyx - xyyy;

/* diagonalize c */

    jacobi (c, d, v, nr);

/* extract the desired quaternion */

    q[0] = v[0][3];
    q[1] = v[1][3];
    q[2] = v[2][3];
    q[3] = v[3][3];

/* generate the rotation matrix */

    q2mat (q, u);

}

/*=======================================================*/


/* FITEST
   rigid fit test driver
   reads in data, fits, and writes out
*/

float gomp_QuatFit(const float *XCoord1 , const float *YCoord1    , const float *ZCoord1 , int Long1 ,
                 float *XCoord2 , float *YCoord2    , float *ZCoord2 , int Long2 ,
                 const int *AList1  , int    AListLong1 , 
                 const int *AList2  , int    AListLong2 ,
                 int    DisplayON)

{
    int n_fields_r;                /* no of fields in xmol file for ref. molec */
    int n_fields_f;                /* no of fields in xmol file for fit. molec */
    int nat_r;                     /* number of all atoms in reference molecule */
    int nat_f;                     /* number of all atoms in fitted molecule */
    double  *xyz_r[4];              /* coordinates for reference molecule */
    double  *xyz_f[4];              /* coordinates for fitted molecule */
    int npairs;                    /* no of fitted atom pairs */
    const int *atoms_r;                  /* atoms of ref. molecule to be superimposed */
    const int *atoms_f;                  /* atoms of fit. molecule to be superimposed */
    double  *ref_xyz[4];            /* ref. molecule atom coordinates to fit */
    double  *fit_xyz[4];            /* fit. molecule atom coordinates to fit */
    double  *weight;                /* fitted atom pair weights */
    double ref_center[4];          /* center of ref. molecule fitted atoms */
    double fit_center[4];          /* center of ref. molecule fitted atoms */
    double q[4];                   /* quaternion */
    double u[4][4];                /* left rotation matrix for coordinates */
    int i, j;                      /* aux variables */
    double s, d, wd, rms, wnorm;   /* aux variables */
    int max_sweeps;                /* max number of iterations in gomp_jacobi */
    int read_pairs;                /* 1 - read pairs, 0 - make pairs, weights */
    int ic,jc;
    char OutText[BUFF_LEN];
    char Name1[BUFF_LEN];
    char Name2[BUFF_LEN];

    /* set defaults */
    max_sweeps = 30;
    read_pairs = 1;

    /* Now read in the ref molecule */


    nat_r = Long1;
    n_fields_r = 4;

    /* read coordinates of ref molecule */

    xyz_r[1] = gomp_AllocateDoubleVector(Long1 + 1);
    xyz_r[2] = gomp_AllocateDoubleVector(Long1 + 1);
    xyz_r[3] = gomp_AllocateDoubleVector(Long1 + 1);

    if(xyz_r[1] == NULL ||
       xyz_r[2] == NULL ||
       xyz_r[3] == NULL) {
        gomp_PrintMessage("?ERROR - can't get space for the reference molecule");
        return((float)-1.0);
    }

    for (i = 1; i <=  nat_r; i++) {
        xyz_r[1][i] = XCoord1[i-1];
        xyz_r[2][i] = YCoord1[i-1];
        xyz_r[3][i] = ZCoord1[i-1];
    }
   
    /* Now read in the fitted molecule */

    xyz_f[1] = gomp_AllocateDoubleVector(Long2 + 1);
    xyz_f[2] = gomp_AllocateDoubleVector(Long2 + 1);
    xyz_f[3] = gomp_AllocateDoubleVector(Long2 + 1);

    if(xyz_f[1] == NULL ||
       xyz_f[2] == NULL ||
       xyz_f[3] == NULL) {
        gomp_PrintMessage("?ERROR - can't get space for the reference molecule");
        return((float)-1.0);
    }

    nat_f = Long2;
    n_fields_f = 4;

    for (i = 1; i <=  nat_f; i++) {
        xyz_f[1][i] = XCoord2[i-1];
        xyz_f[2][i] = YCoord2[i-1];
        xyz_f[3][i] = ZCoord2[i-1];
    }
 
    npairs  = AListLong1;

    if(npairs < 2) {
        sprintf(OutText,
                "Error: Cannot fit a single atom. Need at least 2");
        gomp_PrintMessage(OutText);
        free(xyz_r[1]);
        free(xyz_r[2]);
        free(xyz_r[3]);
        free(xyz_f[1]);  
        free(xyz_f[2]);  
        free(xyz_f[3]);  
        return((float)-1.0);
    }

    atoms_r = AList1;
    atoms_f = AList2;

    weight = gomp_AllocateDoubleVector(AListLong1 + 1);

    if(weight == NULL) {
        gomp_PrintMessage("?ERROR - can't get space for the weights");
        free(xyz_r[1]);
        free(xyz_r[2]);
        free(xyz_r[3]);
        free(xyz_f[1]);  
        free(xyz_f[2]);  
        free(xyz_f[3]);  
        return((float)-1.0);
    }

    for(i = 1; i <= npairs; i++) {
        weight[i] = 1.0;
    }

    ref_xyz[1] = gomp_AllocateDoubleVector(AListLong1 + 1);
    ref_xyz[2] = gomp_AllocateDoubleVector(AListLong1 + 1);
    ref_xyz[3] = gomp_AllocateDoubleVector(AListLong1 + 1);

    fit_xyz[1] = gomp_AllocateDoubleVector(AListLong1 + 1);
    fit_xyz[2] = gomp_AllocateDoubleVector(AListLong1 + 1);
    fit_xyz[3] = gomp_AllocateDoubleVector(AListLong1 + 1);

    /* extract fitted atoms to tables */
    for (i = 1; i <= npairs; i++) {
        for (j = 1; j <= 3; j++) {
            ref_xyz[j][i] = xyz_r[j][atoms_r[i-1] + 1];
            fit_xyz[j][i] = xyz_f[j][atoms_f[i-1] + 1];
        }
    }

    /* ===  Atom coordinates are fit in both modes === */
    /* center ref molecule fitted atoms around (0,0,0) */
    center (npairs, ref_xyz, weight, 1, ref_center);

    /* center fitted molecule fitted atoms around (0,0,0) */
    center (npairs, fit_xyz, weight, 1, fit_center);

    /* fit specified atoms of fit_molecule to those of ref_molecule */
    qtrfit(npairs, fit_xyz, ref_xyz, weight, q, u, max_sweeps);

    /* subtract coordinates of the center of fitted atoms of the fitted molecule
       from all atom coordinates of the fitted molecule (note that weight is
       a dummy parameter) */
    center(nat_f, xyz_f, weight, 2, fit_center);

    /* rotate the fitted molecule by the rotation matrix u */
    rotmol(nat_f, xyz_f, xyz_f, u);
    /* same with set of fitted atoms of the fitted molecule */
    rotmol(npairs, fit_xyz, fit_xyz, u);

    /* translate atoms of the fitted molecule to the center
       of fitted atoms of the reference molecule */
    center(nat_f, xyz_f, weight, 3, ref_center);
    /* same with set of fitted atoms of the fitted molecule */
    center(npairs, fit_xyz, weight, 3, ref_center);
    /* translate fitted atoms of reference molecule to their orig. location */
    center(npairs, ref_xyz, weight, 3, ref_center);
  
    /* find distances between fitted and reference atoms and print them in
       out file */

    if(DisplayON) {  
        sprintf(OutText,
                "Distances and weighted distances between fitted atoms");
        gomp_PrintMessage(OutText);
        sprintf(OutText,"Ref.At.             Fit.At.              Distance  Dist*sqrt(weight)  weight");
        gomp_PrintMessage(OutText);
    }

    rms = 0.0;
    wnorm = 0.0;
    for (i = 1; i <= npairs; i++) {
        d = 0.0;
        for (j = 1; j <= 3; j++) {
            s = ref_xyz[j][i] - fit_xyz[j][i];
            d += s*s;
        }
        d = sqrt(d);
        s = sqrt(weight[i]);
        wd = s*d;
        rms += wd*wd;
        wnorm += s*s;
        ic = atoms_r[i-1];
        jc = atoms_f[i-1];

/*
  strncpy(Name1,sprint_names(ic),BUFF_LEN-1);
  strncpy(Name2,sprint_names(jc),BUFF_LEN-1);
*/
        gomp_CopyString(Name1,"Atom1",BUFF_LEN);
        gomp_CopyString(Name2,"Atom2",BUFF_LEN);

        if(DisplayON) {
            sprintf(OutText, "%s %s  %11.6f  %11.6f  %11.6f",
                    Name1, Name2, d, wd, weight[i]);
            gomp_PrintMessage(OutText);
        }
    }
  
    rms = sqrt(rms/wnorm);

    if(DisplayON) {
        gomp_PrintMessage(" ");
        gomp_PrintMessage(" ");
        sprintf(OutText, "Weighted root mean square=%10.6f", rms);
        gomp_PrintMessage(OutText);
        gomp_PrintMessage(" ");
        gomp_PrintMessage("Center of reference molecule fitted atoms");
        sprintf(OutText, "Xc = %11.6f Yc = %11.6f Zc = %11.6f",
                ref_center[1], ref_center[2], ref_center[3]);
        gomp_PrintMessage(OutText);
        gomp_PrintMessage(" ");
        gomp_PrintMessage(" ");
        sprintf(OutText, "Center of fitted molecule fitted atoms");
        gomp_PrintMessage(OutText);
        sprintf(OutText, "Xc = %11.6f Yc = %11.6f Zc = %11.6f",
                fit_center[1], fit_center[2], fit_center[3]);
        gomp_PrintMessage(OutText);

        gomp_PrintMessage(" ");
        gomp_PrintMessage(" ");
        sprintf(OutText,"Left rotation matrix");
        gomp_PrintMessage(OutText);
        for (i = 1; i <= 3; i++) {
            sprintf(OutText, " %11.6f  %11.6f  %11.6f", 
                    u[1][i], u[2][i], u[3][i]);
            gomp_PrintMessage(OutText);
        }
    }

    for (i = 1; i <=  nat_f; i++) {
        XCoord2[i-1] =  xyz_f[1][i];
        YCoord2[i-1] =  xyz_f[2][i];
        ZCoord2[i-1] =  xyz_f[3][i];
    }


    free(xyz_r[1]);  
    free(xyz_r[2]);  
    free(xyz_r[3]);  

    free(xyz_f[1]);  
    free(xyz_f[2]);  
    free(xyz_f[3]);  

    free(ref_xyz[1]);
    free(ref_xyz[2]);
    free(ref_xyz[3]);

    free(fit_xyz[1]);
    free(fit_xyz[2]);
    free(fit_xyz[3]);

    free(weight);

    return((float)rms);

}


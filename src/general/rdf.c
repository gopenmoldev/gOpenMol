/*

Copyright (c) 1990 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

/*
  calculate atomic radial distribution function for a single
  configuration. 

  Where box is the length of the cubic periodic box 
  cutoff is the limit up to which the rdf is calculated
  (the lower limit is zero)

  fmp oct 89
  modified to c by 
  Leif Laaksonen dec 89


  Na      = 6.0221 10^23/mol
  density = 1000 kg/mol
  water   = 18.0073 g/mol

  ====>  0.033442 molecules/A^3


  1999-09-27: Changed to look for connected atoms in stead of
  same residue number, which works really only for
  CHARMM crd files.
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include "gommath.h"
#include <string.h>

#include "gomclipbrd.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "printmsg.h"
#include "rdf.h"
#include "selection.h"

#include "stdafx.h"

static struct {
    int    Mean;                   /* != 0 average is calulated     */
    int    Plot;                   /* != 0 plot is on               */
    int    Set;                    /* != 0 collection is on         */
    int    Numbers;                /* number of points along x-axis */
    float *XValues;                /* X - pointer to values         */
    float *YValues;                /* Y - pointer to values         */
} RaDiFu = { 0 , 0 , 0 , 0 , NULL , NULL};

/*************************************************************************/
int  gomp_CalcRDF(const int *fromi,int fromlong,float boxl1,float boxl2,float boxl3,
                float rcut,int mbin,
                const char *Text1,
                const char *Text2,
                const char *Text3)
/*
  const int *fromi;      distribution function is calculated for atom i
  int   fromlong;   length of the fromi vector
  float boxl1;      box dimension x 
  float boxl2;      box dimension y 
  float boxl3;      box dimension z 
  float rcut;       cut off distance
  int   mbin;       number of rdf values (bins)
  const char *Text1;      Segment name (of the atom(s) around i) 
  const char *Text2;      Residue name (of the atom(s) around i)
  const char *Text3;      Atom name    (of the atom(s) around i)  
*/
/*************************************************************************/
{

    static int i,j,k,l,id,il;
    static int *ibin;                 /* pointer to number array */
    static float *g;                  /* pointer to distribution funct */
    static float *x_axis;
    static float *y_axis;
    static float delta;               /* step length             */
    static float d,dx,dy,dz,xboxl,yboxl,zboxl;
    static float fact,rl,ru,rideal;
    static float Tmp1,Tmp2;
    static int atom_list,slong;
    static int *sel_list;
    static const float *x;
    static const float *y;
    static const float *z;
    static int    Wstr;
    static const int *cnct_run1;
    static int    c_1;
    static const int *cnct_run2;
    static int    c_2;
    static int    hit1;
    static int    hit2;

    if(RaDiFu.Set && (mbin != RaDiFu.Numbers)) {
        gomp_PrintMessage("?ERROR - RDFs are collected and new nbin != old");
        gomp_PrintMessage(" Will gomp_rap old values and start collecting new");
        (void)gomp_DeleteRDF();
    }

    if(RaDiFu.Set && RaDiFu.Mean) {
        gomp_PrintMessage(" Average RDF is already calulated");
        gomp_PrintMessage(" Will start collecting new series        ");
        (void)gomp_DeleteRDF();
    }

    if(!RaDiFu.Set) {
        RaDiFu.XValues = gomp_AllocateFloatVector(mbin);
        RaDiFu.YValues = gomp_AllocateFloatVector(mbin);

        for(i = 0 ; i < mbin ; i++) RaDiFu.YValues[i] = 0.0;
    }

    atom_list = gomp_GetNumAtomsInMolecStruct(0);

    sel_list  = gomp_AllocateIntVector(atom_list);

    slong     = gomp_MakeSelectionList(0,Text1,Text2,Text3,sel_list);

    if(!slong) {
        gomp_PrintERROR("no atoms in the second selection list");
        free(RaDiFu.XValues);
        free(RaDiFu.YValues);
        free(sel_list);
        return(1);
    }

    RaDiFu.Numbers = mbin;

    delta     = rcut / (float)(mbin-1) ;

    Wstr      = 0;
    x         = gomp_GetAtomXCoordPointer(Wstr);
    y         = gomp_GetAtomYCoordPointer(Wstr);
    z         = gomp_GetAtomZCoordPointer(Wstr);

/*     ... initialise */
       
    ibin   = gomp_AllocateIntVector(mbin);       /* reserve memory for mbin elements */
    g      = gomp_AllocateFloatVector(mbin);       /* reserve memory for dist funct    */
    x_axis = gomp_AllocateFloatVector(mbin);
    y_axis = gomp_AllocateFloatVector(mbin);

    for(i = 0 ; i < mbin  ; i++)  ibin[i] = 0;

    xboxl = 1. / boxl1;
    yboxl = 1. / boxl2;
    zboxl = 1. / boxl3;
/*     ... calculate rdf      */

/* apply always to FIRST structure */

    for(il = 0 ; il < fromlong ; il++) {
        l = fromi[il];
/*         i = gomp_GetAtomResNum1(Wstr , l);*/

        for(k = 0 ; k < slong ; k++) {
            j = sel_list[k];

            if(l == j) continue;

/*              if(gomp_GetAtomResNum1(Wstr , j) == i) continue;*/
            cnct_run1 = gomp_GetAtomConnection(Wstr, l);
            c_1 = *cnct_run1++;
            hit1 = 0;
            if(c_1) {
                for(i = 0 ; i < c_1 ; i++) {
                    if(*cnct_run1++ == j) {
                        hit1 = 1;
                        break;
                    }
                }
            }

            if(hit1) continue;
              
            cnct_run2 = gomp_GetAtomConnection(Wstr, j);
            c_2 = *cnct_run2++;
            hit2 = 0;
            if(c_2) {
                for(i = 0 ; i < c_2 ; i++) {
                    if(*cnct_run2++ == l) { 
                        hit2 = 1;
                        break;
                    }
                }
            }

            if(hit2) continue;

            dx = x[l] - x[j];
            dy = y[l] - y[j];
            dz = z[l] - z[j];
  
            dx = dx - boxl1 * nearbyint(xboxl * dx);
            dy = dy - boxl2 * nearbyint(yboxl * dy);
            dz = dz - boxl3 * nearbyint(zboxl * dz);
  
            d = sqrt (dx * dx + dy * dy + dz * dz);
            id = (int)(d / delta) ;

            if (id < mbin ) ibin[id] = ibin[id] + 1; /* I use 1 here not 2!!*/

        }
    }

/*
  fact = M_PI * 4.0 * WATER_DENS / 3.0;
*/
    Tmp1 = (float)(fromlong);
    Tmp2 = (float)(slong);

    fact = M_PI * 4.0 * Tmp2 / (3.0 * boxl1 * boxl2 * boxl3);

    for(i = 0 ; i < mbin ; i++) {
        rl = i * delta;
        x_axis[i] = rl;
/*         y_axis[i] = (float)(ibin[i]); */
        ru = rl + delta;
        rideal = fact * (ru*ru*ru - rl*rl*rl);
        g[i] = y_axis[i] =  (float)(ibin[i]) / (Tmp1 * rideal);
    }

    for(i = 0 ; i < mbin ; i++) {
        RaDiFu.XValues[i]  = i * delta;
        RaDiFu.YValues[i] += g[i];
    }

    RaDiFu.Set++;

    free(g);
    free(ibin);
    free(x_axis);
    free(y_axis);
    free(sel_list);

    return(0);
}

/*************************************************************************/
int gomp_GetNumRDFObs()
/*************************************************************************/
{
    return(RaDiFu.Numbers);
}
/*************************************************************************/
const float *gomp_GetRDFVecX()
/*************************************************************************/
{
    return(RaDiFu.XValues);
}
/*************************************************************************/
const float *gomp_GetRDFVecY()
/*************************************************************************/
{
    return(RaDiFu.YValues);
}

/************************************************************************/
int gomp_DeleteRDF()
/************************************************************************/
{
    if(!RaDiFu.Set) {
        return(1);
    }
    RaDiFu.Set     = 0;
    RaDiFu.Numbers = 0;
    RaDiFu.Plot    = 0;
    RaDiFu.Mean    = 0;
    if(RaDiFu.Numbers) {
        free(RaDiFu.XValues);
        free(RaDiFu.YValues);
    }

    return(0);
}
/************************************************************************/
int gomp_CalcMeanRDF()
/************************************************************************/
{
    int i;
    float RDFsets;
    char OutText[BUFF_LEN];

    if(!RaDiFu.Set || !RaDiFu.Numbers) {
        gomp_PrintERROR("no rdf defined ");
        return(1);
    }

    if(RaDiFu.Set && RaDiFu.Mean) {
        gomp_PrintERROR("the average value is already calculated");
        return(1);
    }

    gomp_PrintMessage("Calculating the average values for the RDF");
    sprintf(OutText,"Number of sets: %d",RaDiFu.Set);
    gomp_PrintMessage(OutText);   

    RDFsets = (float)(RaDiFu.Set);

    for(i = 0 ; i < RaDiFu.Numbers ; i++) 
        RaDiFu.YValues[i] /= RDFsets;


    RaDiFu.Mean = 1;  /* it's  done now */

    return(0);
}

/************************************************************************/
int gomp_WriteRDF(const char *FileName)
/************************************************************************/
{
    FILE *write_p;
    int i;
    char OutText[BUFF_LEN];

    if(!RaDiFu.Set || !RaDiFu.Numbers) {
        gomp_PrintERROR("no RDF calculated to be written to disk");
        return(1);
    }

    if(FileName[0] == (char)NULL) {
        gomp_PrintERROR("rdf file name is missing");
        return(1);
    }

    write_p = fopen(FileName,"w");
    if(write_p == NULL) {
        sprintf(OutText,"?ERROR - unable to open file : %s ",FileName);
        gomp_PrintERROR(OutText);
        return(1);
    }

/* ready to write now ... */

    sprintf(OutText,"Writing the RDF to file '%s'",FileName);
    gomp_PrintMessage(OutText);

    for(i = 0 ; i < RaDiFu.Numbers ; i++) {
        fprintf(write_p," %f %f \n",RaDiFu.XValues[i],
                RaDiFu.YValues[i]);
    }


    fclose(write_p);

    return(0);
}

/************************************************************************/
int   gomp_CopyRDFarray2Clipboard()
/************************************************************************/
{
    static int  i;
    static char Text[BUFF_LEN];

    if(!gomp_GetNumRDFObs()) {
        gomp_PrintERROR("no RDF information is available");
        return(1);
    }

    {
        char *String = NULL;
        const float *ValueX;
        const float *ValueY;

        ValueX  = gomp_GetRDFVecX();
        ValueY  = gomp_GetRDFVecY();

        for (i = 0 ; i < gomp_GetNumRDFObs() ; i++) {
#if defined(WIN32)
            sprintf(Text,"%f  %f\r\n",ValueX[i],ValueY[i]);
#else
            sprintf(Text,"%f  %f\n",ValueX[i],ValueY[i]);
#endif

            if(!i) {
                String = gomp_AllocateCharVector(strlen(Text) + 1);
                strncpy(String,Text,strlen(Text));
                String[strlen(Text)] = (char)NULL;
            }
            else {
                String = gomp_ReallocateCharVector(String , 
                                                strlen(Text) + strlen(String) + 1);
                strncat(String,Text,strlen(Text));
                String[strlen(String)] = (char)NULL;
            }
        }
        (void)gomp_CopyText2Clipboard(String);
        free(String);
    }

    return(0);
}
/************************************************************************/
int    gomp_GetRDFstatus()
/************************************************************************/
{
    return(RaDiFu.Set);
}
/************************************************************************/
int    gomp_GetRDFobservations()
/************************************************************************/
{
    return(RaDiFu.Set);
}

/************************************************************************/
int    gomp_GetRDFaverageStatus()
/************************************************************************/
{
    return(RaDiFu.Mean);
}

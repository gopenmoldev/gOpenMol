/*

Copyright (c) 1990 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2005 by:
Eero Häkkinen
*/

/*
  This program reads a Charmm and flat file vector file

  Leif Laaksonen 1989, 1995, 2001
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <tcl.h>

#include "colouring.h"
#include "gomfile.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "objseg.h"
#include "plot.h"
#include "printmsg.h"
#include "rforce.h"
#include "tclutils.h"

#include "stdafx.h"

#define KARP_LINE_LEN   120   /* vector file line length */

#define VECTOR_ON  1
#define VECTOR_OFF 0

#define RABS(a)    ( ( a ) > 0 ? (a) : -(a))

#if 0
/* define the CHARMm structure                          */
struct CHARMm { /* Charmm structure */
    int numat;          /* number of atoms in this structure */
    int *res1;          /* residue number 1 list */
    int *res2;          /* residue number 2 list */
    char *resnam;       /* residue name list */
    char *atnam;        /* atom name list */
    char *segment;      /* segment name list */
    float *x;           /* x,y and z coordinates */
    float *y;
    float *z;
    float *bvalue;
};      /* bvalue list */

static struct  CHARMm TCHARMm;

static int GetCHARMmStructSpace(int);
static int DelCHARMmStruct(struct CHARMm);
#endif

/* end of CHARMm structure                              */

/* pointers for vector/force calculation */

static struct {
    gom_Plotter *CallbackHandle; /* NULL , no plot , != 0 a plot */
    float *fx;      /* pointers to the vectors */
    float *fy;
    float *fz;
    int      wstr;     /* structure number */
    int *sel_list; /* selection list */
    int      ent_list; /* entries in the selection list */
    int      maxfi;    /* index of max vector atom */
    float    maxfa;    /* value of the max vector  */
    int      minfi;    /* index on min vector atom */
    float    minfa;    /* value of the min vector  */ 
    float    scale;    /* scale factor */
    float    radius;   /* radius of cylinder */
    float    range_min; /* min value for display range */
    float    range_max; /* max value for display range */
} atm_vector;

/*
  Read in vectors from a file in charmm format

  Leif Laaksonen 1990, 1995, 2001
*/
/*************************************************************************/
int gomp_ReadCharmmVector(const char *inp_file)
/*************************************************************************/
{
    static char inputl[KARP_LINE_LEN];
    static int tatomn,tres1,tres2;
    static int i,j;
    static int type_warning;
    static float tf;

    static char tseg[BUFF_LEN];
    static char tres[BUFF_LEN];
    static char tatm[BUFF_LEN];

    static char OutText[BUFF_LEN];

    static FILE *chm_in;

    static int   Wstr;


    type_warning = 0;
    Wstr         = 0;

    chm_in=fopen(inp_file,"r");
    if(chm_in == NULL) {
        sprintf(OutText,"Can't open input file : %s ",inp_file);
        gomp_PrintMessage(OutText);
        return(1);
    }

/*
  Start reading file

*/

/*  1.0 A title is expected    */

    fgets(inputl,KARP_LINE_LEN,chm_in);
    gomp_PrintMessage(inputl);

    while(strncmp(inputl,"*",1) == 0) {
        fgets(inputl,KARP_LINE_LEN,chm_in);
        if(strncmp(inputl,"*",1) != 0) break;
        if(inputl[strlen(inputl) - 1] == '\n') 
            inputl[strlen(inputl) - 1] = '\0';
        sprintf(OutText,"%s",inputl);
        gomp_PrintMessage(OutText);
    }
/*   2.0 Numbers of atoms */

    sscanf(inputl,"%d",&tatomn);

/* check that the number of atoms match the current molecule */
    if(tatomn != gomp_GetNumAtomsInMolecStruct(Wstr)) {
        gomp_PrintMessage("?ERROR - number of atoms in the Vector file does not match current molecule\n");
        return(1);
    } 

/*   3.0 Read vector cards  */
    (void)gomp_DeleteVectorStructure();

    atm_vector.sel_list = gomp_AllocateIntVector(tatomn);
    atm_vector.fx       = gomp_AllocateFloatVector(tatomn);
    atm_vector.fy      = gomp_AllocateFloatVector(tatomn);
    atm_vector.fz     = gomp_AllocateFloatVector(tatomn);

    atm_vector.ent_list  = tatomn;
    atm_vector.wstr      = 0;

    for(i = 0 ; i < tatomn; i++ ) {

        atm_vector.sel_list[i] = i;

        fgets(inputl,KARP_LINE_LEN,chm_in);
        sscanf(inputl,"%d %d %s %s %f %f %f %s %d",
               &j,&tres1,tres,tatm,&atm_vector.fx[i],
               &atm_vector.fy[i],
               &atm_vector.fz[i],tseg,&tres2);

    }
/* look for the min/max */
    atm_vector.maxfi =  -1;
    atm_vector.maxfa = 0.0;
    atm_vector.minfi =  -1;
    atm_vector.minfa = 1.e+20f;

    for( i = 0 ; i < tatomn ; i++) {
        tf = sqrt(atm_vector.fx[i]*atm_vector.fx[i] +
                  atm_vector.fy[i]*atm_vector.fy[i] +
                  atm_vector.fz[i]*atm_vector.fz[i]);
        if(tf > atm_vector.maxfa) {
            atm_vector.maxfa = tf;
            atm_vector.maxfi = i ;
        }

        if(tf < atm_vector.minfa) {
            atm_vector.minfa = tf;
            atm_vector.minfi = i ;
        }
    }

    sprintf(OutText,"*** Vector statistics from file '%s' ",inp_file);
    gomp_PrintMessage(OutText);

    if(atm_vector.maxfi >= 0) {
        i = atm_vector.maxfi;
        sprintf(OutText,"Max Vector is for %s:%s(%d):%s(%d) : %f",
                gomp_GetAtomSegName(Wstr,i),gomp_GetAtomResName(Wstr,i),gomp_GetAtomResNum1(Wstr,i),
                gomp_GetAtomAtmName(Wstr,i),(i+1),atm_vector.maxfa);
        gomp_PrintMessage(OutText);
    }
    else
        gomp_PrintMessage("?ERROR - problems in assigning max Vector");

    if(atm_vector.minfi >= 0) {
        i = atm_vector.minfi;
        sprintf(OutText,"Min Vector is for %s:%s(%d):%s(%d) : %f",
                gomp_GetAtomSegName(Wstr,i),gomp_GetAtomResName(Wstr,i),gomp_GetAtomResNum1(Wstr,i),
                gomp_GetAtomAtmName(Wstr,i),(i+1),atm_vector.minfa);
        gomp_PrintMessage(OutText);
    }
    else
        gomp_PrintMessage("?ERROR - problems in assigning min Vector");

/*      atm_vector.scale = (atm_vector.maxfa - atm_vector.minfa) / 3.;*/
    atm_vector.scale  = 1.0;
    atm_vector.radius = 0.5;

    gomp_PrintMessage("**********   Done   **********\n");

    fclose(chm_in);

    return(0);
}
#if 0
/*************************************************************************/
int GetCHARMmStructSpace(int Atoms)
/*************************************************************************/
{

/* residue number 1 list */
    TCHARMm.res1   =   gomp_AllocateIntVector(Atoms);
/* residue number 2 list */
    TCHARMm.res2   =   gomp_AllocateIntVector(Atoms);
/* residue name list     */
    TCHARMm.resnam =   gomp_AllocateCharVector(Atoms * MAX_RES_NAME_LEN);
/* atom name list        */
    TCHARMm.atnam  =   gomp_AllocateCharVector(Atoms * MAX_ATM_NAME_LEN);
/* segment name list     */
    TCHARMm.segment =  gomp_AllocateCharVector(Atoms * MAX_SEG_NAME_LEN);
/* x coordinate          */
    TCHARMm.x       =  gomp_AllocateFloatVector(Atoms);
/* y coordinate          */
    TCHARMm.y       =  gomp_AllocateFloatVector(Atoms);
/* z coordinate          */
    TCHARMm.z       =  gomp_AllocateFloatVector(Atoms);
/* bvalue list           */
    TCHARMm.bvalue  =  gomp_AllocateFloatVector(Atoms);

    return(0);
}
/*************************************************************************/
int DelCHARMmStruct(struct CHARMm DelCHARMm)
/*************************************************************************/
{
    if(DelCHARMm.numat < 1) {
        gomp_PrintMessage("?ERROR - no structure defined to be deleted");
        return(1);
    }

    free(DelCHARMm.res1);
    free(DelCHARMm.res2);
    free(DelCHARMm.resnam);
    free(DelCHARMm.atnam);
    free(DelCHARMm.segment);
    free(DelCHARMm.x);
    free(DelCHARMm.y);
    free(DelCHARMm.z);
    free(DelCHARMm.bvalue);

    return(0);
}
#endif
/************************************************************************/
int gomp_ColorByVector(int Wstr, int slong, const int *sel_list ,
                     double *Tmin,double *Tmax)
/************************************************************************/
{
    float delta,Vector;
    double step;
    int  i,j;
    float  rr,gg,bb;
    char OutText[BUFF_LEN];
    float  Vmin;
    float  Vmax;

    if(!atm_vector.ent_list) {
        gomp_PrintERROR("No data available in the CHARMm vector array");
        return(1);
    }

    if(Tmax == NULL) 
        Vmax = atm_vector.maxfa;
    else
        Vmax = *Tmax;

    if(Tmin == NULL) 
        Vmin = atm_vector.minfa;
    else
        Vmin = *Tmin;

    delta   = Vmax - Vmin;

    sprintf(OutText,"Max Vector: %7.3f (red)",Vmax);
    gomp_PrintMessage(OutText);
    sprintf(OutText,"Min Vector: %7.3f (blue)",Vmin);
    gomp_PrintMessage(OutText);

    for(i = 0 ; i < slong ; i++) {
        j = sel_list[i];

        Vector =  sqrt(atm_vector.fx[j]*atm_vector.fx[j] +
                       atm_vector.fy[j]*atm_vector.fy[j] +
                       atm_vector.fz[j]*atm_vector.fz[j]);

        step = (double)((Vector - Vmin) / delta);

        gomp_PreRainbow(step,&rr,&gg,&bb);
        gomp_PutAtomColour(Wstr, rr, gg, bb, j);
    }

    return(0);
}
/*************************************************************************/
int gomp_DeleteVectorStructure()
/*************************************************************************/
{
    if( atm_vector.sel_list ) {
        free(atm_vector.fx);
        free(atm_vector.fy);
        free(atm_vector.fz);

        free(atm_vector.sel_list);

        atm_vector.sel_list = NULL;
    }

    gomp_UnregisterPlotter(atm_vector.CallbackHandle);
    atm_vector.CallbackHandle = NULL;

    (void)gomp_DelLineArrowSeg();

    return(0);
}
/*************************************************************************/
const float *gomp_GetVectorForceXp()
/*************************************************************************/
{
    return(atm_vector.fx);
}
/*************************************************************************/
const float *gomp_GetVectorForceYp()
/*************************************************************************/
{
    return(atm_vector.fy);
}

/*************************************************************************/
const float *gomp_GetVectorForceZp()
/*************************************************************************/
{
    return(atm_vector.fz);
}

/*************************************************************************/
int     gomp_SetVectorListLength(int Length)
/*************************************************************************/
{
    atm_vector.ent_list = Length;

    return(0);
}
/*************************************************************************/
int     gomp_GetVectorListLength()
/*************************************************************************/
{
    return(atm_vector.ent_list);
}
/*************************************************************************/
/* gomp_GetVectorListArray                                                 */
/*************************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(int,VectorListArray,(void),
    atm_vector.sel_list,;)
/*************************************************************************/
int     gomp_GetVectorStructureIndex()
/*************************************************************************/
{
    return(atm_vector.wstr);
}

/*************************************************************************/
int     gomp_SetPlotVectorStatus(int PStatus)
/*************************************************************************/
{
    if( PStatus ) {
        if(gomp_GetVectorListLength() || gomp_GetTotalLineArrowEntries()) {
            if( !atm_vector.CallbackHandle )
                atm_vector.CallbackHandle =
                    gomp_RegisterPlotter(
                        gomp_PlotVectorData,NULL,
                        PLOTTER_NAME_VECTOR,PLOTTER_ORDER_VECTOR);
        }
    }
    else {
        gomp_UnregisterPlotter(atm_vector.CallbackHandle);
        atm_vector.CallbackHandle = NULL;
    }

    return(0);
}
/*************************************************************************/
int gomp_GetPlotVectorStatus()
/*************************************************************************/
{
    return(atm_vector.CallbackHandle!=0);
}
/*************************************************************************/
int     gomp_SetVectorScale(float Scale)
/*************************************************************************/
{
    atm_vector.scale = Scale;

    return(0);
}
/*************************************************************************/
float   gomp_GetVectorScale()
/*************************************************************************/
{
    return(atm_vector.scale);
}

/*************************************************************************/
int     gomp_SetVectorRadius(float Radius)
/*************************************************************************/
{
    atm_vector.radius = Radius;

    return(0);
}
/*************************************************************************/
float   gomp_GetVectorRadius()
/*************************************************************************/
{
    return(atm_vector.radius);
}
/*************************************************************************/
int gomp_ReadFlatFileVector(const char *inp_file)
/*************************************************************************/
{

    static char   inputl[KARP_LINE_LEN];
    static int    i;
    static float  xc1,yc1,zc1,xc2,yc2,zc2;
    static char   OutText[BUFF_LEN];
    static FILE  *chm_in;
    static float  tf;
    static const char *value;
    static int    DefColorIsGiven;


    chm_in=fopen(inp_file,"r");
    if(chm_in == NULL) {
        sprintf(OutText,"Can't open input file : %s ",inp_file);
        gomp_PrintMessage(OutText);
        return(1);
    }

    atm_vector.scale      = 1.0;
    atm_vector.maxfa      = 0.0;
    atm_vector.minfa      = 1.e+25f;
    atm_vector.range_min  = 0.0;
    atm_vector.range_max  = 1.e+25f;

    DefColorIsGiven = 0;
    value  = Tcl_GetVar(gomp_GetTclInterp() , "gomDefaultGradienColor", TCL_GLOBAL_ONLY);

    if(value) {
        DefColorIsGiven = 1;
    } 


/*
  Start reading file

*/
/* first five lines of general information */
    gomp_PrintMessage("File header information:");
    gomp_Fgets(inputl,KARP_LINE_LEN,chm_in);
    gomp_PrintMessage(inputl);
    gomp_Fgets(inputl,KARP_LINE_LEN,chm_in);
    gomp_PrintMessage(inputl);
    gomp_Fgets(inputl,KARP_LINE_LEN,chm_in);
    gomp_PrintMessage(inputl);
    gomp_Fgets(inputl,KARP_LINE_LEN,chm_in);
    gomp_PrintMessage(inputl);
    gomp_Fgets(inputl,KARP_LINE_LEN,chm_in);
    gomp_PrintMessage(inputl);
/* done! */
    i = 0;
    while(fgets(inputl,KARP_LINE_LEN,chm_in) != NULL) {

        sscanf(inputl,"%f %f %f %f %f %f", &xc1 , &yc1 , &zc1,
               &xc2 , &yc2 , &zc2);

        if(DefColorIsGiven) {
            if(gomp_PushLineArrowStack(xc1,yc1,zc1, xc2,yc2,zc2,
                                  value,"append")) {
                gomp_PrintERROR("can't push new entries into the line arrow stack");
                return(1);
            }
        } else {
            if(gomp_PushLineArrowStack(xc1,yc1,zc1, xc2,yc2,zc2,
                                  "append","")) {
                gomp_PrintERROR("can't push new entries into the line arrow stack");
                return(1);
            }
        }

        tf = sqrt(xc2 * xc2 + yc2 * yc2 + zc2 * zc2);
        if(tf > atm_vector.maxfa) {
            atm_vector.maxfa = tf;
        }

        if(tf < atm_vector.minfa) {
            atm_vector.minfa = tf;
        }


        i++;
    }

    atm_vector.range_max = atm_vector.maxfa;
    atm_vector.range_min = atm_vector.minfa;

    gomp_PrintMessage("**********   Done   **********");
    sprintf(OutText,"Found %d entries in the input file\n",i);
    gomp_PrintMessage(OutText);

    sprintf(OutText,"*** Vector statistics from file '%s' ",inp_file);
    gomp_PrintMessage(OutText);

    sprintf(OutText,"Max vector norm: %f",atm_vector.maxfa);
    gomp_PrintMessage(OutText);

    sprintf(OutText,"Min vector norm: %f",atm_vector.minfa);
    gomp_PrintMessage(OutText);

    fclose(chm_in);

    return(0);
}
/*************************************************************************/
int     gomp_SetVectorDisplayRange(float MinValue , float MaxValue)
/*************************************************************************/
{
    atm_vector.range_min = MinValue;
    atm_vector.range_max = MaxValue;

    return(0);
}
/*************************************************************************/
int     gomp_GetVectorDisplayRange(float *MinValue, float *MaxValue)
/*************************************************************************/
{
    *MinValue = atm_vector.range_min;
    *MaxValue = atm_vector.range_max;

    return(0);
}


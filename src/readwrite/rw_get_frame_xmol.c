/*

Copyright (c) 1996 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002, 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <tcl.h>

#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"

#include "stdafx.h"

#define MOD(a,b)  (((a+1)/b)*b - (a+1))

static int  Look4XmolFrames(FILE *, int * , int *);
static int  PlaceXMOLFrame(FILE * , int , int );

static int *XmolFrameStartPoint = (int *)NULL;

/***************************************************************************/
int gomp_GetFrameXmol(int alt,FILE *File_p , int iappend)  
    /* read one frame from a xmol trajectory   
       mode of operation (=0) first time in read
       (>0) trajectory number 
       if = 0 no append , if = 1 append */
/***************************************************************************/
{
    static int    natom;
    static int    nstep;
    static int    Wstr;
    static char   title[BUFF_LEN];

/* just to be sure ... */
    rewind(File_p);
/* determine the "record length"    */
    natom = gomp_GetNumAtomsInMolecStruct(0);

/*  start reading  */

    if(!alt) {

        (void)Look4XmolFrames(File_p , &natom , &nstep);

        (void)gomp_SetNumberOfTrajectoryAtoms(natom);
        (void)gomp_SetNumberOfFrames(nstep);
        (void)gomp_SetTrajectoryTimeInfo(0 , 0);
        (void)gomp_SetNumberOfFreeAtoms(0 , 0);
        (void)gomp_SetTrajectoryDisplayParams(1 , nstep , 1);
        (void)gomp_PutDisplayFrameNumber(1);

        sprintf(title," Atoms found                : %d   ",natom);
        gomp_PrintMessage(title);
        sprintf(title," Free atoms                 : %d   ",natom);
        gomp_PrintMessage(title);
        sprintf(title," Dynamics steps             : %d   ",nstep);
        gomp_PrintMessage(title);
      
        return(0);
    }


    if(iappend == 0)
        Wstr = 0;
    else {
        sprintf(title,"Xmol frame (%d)",alt);
        Wstr = gomp_CreateMolecStruct(title , natom , APPEND);
        if ( Wstr < 0 )
            return(1);
    }

/*  get pointer to coordinate vectors */ 
    (void)PlaceXMOLFrame(File_p , alt , Wstr);

    return(0);
/*                                  */
}

#define XMOL_LINE_LEN   120   /*XMOL xyz file line length */

/*********************************************************************/
int Look4XmolFrames(FILE *chm_in , int *Atoms , int *Frames)
/*********************************************************************/
{
    static int i;
    static int XmolAtoms;

    float   TXc,TYc,TZc,TBv;
    char    TAtmN[BUFF_LEN];

    char inputl[XMOL_LINE_LEN];
    int  items;


    rewind(chm_in);
    *Atoms  = 0;
    *Frames = 0;

    free(XmolFrameStartPoint);
    XmolFrameStartPoint  = (int *)NULL;

/*
  Start reading file

*/

/* Save internal pointer in XMOL file */
    XmolFrameStartPoint     = gomp_AllocateIntVector(1);
    XmolFrameStartPoint[0]  = ftell(chm_in);

/*  Two title lines are expected    */

    while(fgets(inputl,XMOL_LINE_LEN,chm_in) != NULL) {

/*   Numbers of atoms */
        sscanf(inputl,"%d",&XmolAtoms);
        fgets(inputl,XMOL_LINE_LEN,chm_in);

        *Atoms = XmolAtoms;

/*   Read atom cards   */

        for(i = 0 ; i < XmolAtoms ; i++ ) {

            fgets(inputl,XMOL_LINE_LEN,chm_in);

            items = sscanf(inputl," %s %f %f %f %f",
                           TAtmN,&TXc,&TYc,&TZc,&TBv);

        }

        *Frames = *Frames + 1;

        XmolFrameStartPoint = 
            gomp_ReallocateIntVector(XmolFrameStartPoint , *Frames + 1);
        XmolFrameStartPoint[*Frames]  = ftell(chm_in);

    }

    rewind(chm_in);

    return(0);
}

/*********************************************************************/
int PlaceXMOLFrame(FILE *chm_in , int Alt , int Wstr)
/*********************************************************************/
{
    static int i,j,k;
    static int XmolAtoms;

    float   TXc,TYc,TZc,TBv;
    const float *sumxyz;
    char    TAtmN[BUFF_LEN];

    char inputl[XMOL_LINE_LEN];
    int  items;
    int  method;

/*
  Start reading file

*/
    rewind(chm_in);
    sumxyz   = gomp_GetTranslateArray();

    method = gomp_GetFormattedTrajectoryReader();

    if(method) {

        for ( k = 0 ; k < Alt ; k++) {

/*  Two title lines are expected    */

/*   Numbers of atoms */
            fgets(inputl,XMOL_LINE_LEN,chm_in);
            sscanf(inputl,"%d",&XmolAtoms);
            fgets(inputl,XMOL_LINE_LEN,chm_in);


/*   Read atom cards   */

            for(i = 0 ; i < XmolAtoms ; i++ ) {
                fgets(inputl,XMOL_LINE_LEN,chm_in);

                if(k == (Alt - 1)) {
                    items = sscanf(inputl," %s %f %f %f %f",
                                   TAtmN,&TXc,&TYc,&TZc,&TBv);

                    j = gomp_PutAtomXCoord(Wstr , TXc - sumxyz[0] , i);
                    j = gomp_PutAtomYCoord(Wstr , TYc - sumxyz[1] , i);
                    j = gomp_PutAtomZCoord(Wstr , TZc - sumxyz[2] , i);

                    j = gomp_PutAtomCharge(Wstr ,  TBv , i);
                }
            }
        }
    } else {

        if(fseek(chm_in , (long)XmolFrameStartPoint[Alt - 1] , SEEK_SET)) {
            gomp_PrintERROR("can't position at frame for XMOL trajectory");
            rewind(chm_in);
            return(1);
        }

        fgets(inputl,XMOL_LINE_LEN,chm_in);
        sscanf(inputl,"%d",&XmolAtoms);
        fgets(inputl,XMOL_LINE_LEN,chm_in);

/*   Read atom cards   */

        for(i = 0 ; i < XmolAtoms ; i++ ) {

            fgets(inputl,XMOL_LINE_LEN,chm_in);

            items = sscanf(inputl," %s %f %f %f %f",
                           TAtmN,&TXc,&TYc,&TZc,&TBv);

            j = gomp_PutAtomXCoord(Wstr , TXc - sumxyz[0] , i);
            j = gomp_PutAtomYCoord(Wstr , TYc - sumxyz[1] , i);
            j = gomp_PutAtomZCoord(Wstr , TZc - sumxyz[2] , i);

            j = gomp_PutAtomCharge(Wstr ,  TBv , i);
        }
    }

    rewind(chm_in);

    return(0);
}

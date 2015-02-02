/*

Copyright (c) 1998 - 2005 by:
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
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"

#include "stdafx.h"

#define MOD(a,b)  (((a+1)/b)*b - (a+1))

static int Look4AMBERAFrames(FILE *, int * , int *);
static int PlaceAMBERAFrame(FILE * , int , int );

static int *AmberaFrameStartPoint = NULL;

/***************************************************************************/
int gomp_GetFrameFAmber(int alt,FILE *File_p , int iappend)  
    /* read one frame from a AMBER formatted trajectory   
       mode of operation (=0) first time in read
       (>0) trajectory number 
       if = 0 no append , if = 1 append */
/***************************************************************************/
{
    static int    natom;
    static int    natomx;
    static int    nstep;
    static int    Wstr;
    static char   title[BUFF_LEN];

/* just to be sure ... */
    rewind(File_p);
/* determine the "record length"    */
    natomx = gomp_GetNumAtomsInMolecStruct(0);

/*  start reading  */

    if(!alt) {

        (void)Look4AMBERAFrames(File_p , &natom , &nstep);

        if(natomx != natom) {
            gomp_PrintERROR("number of atoms in trajectory file do no match current system");
            return(1);
        }
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
        sprintf(title,"AMBERA frame (%d)",alt);
        Wstr = gomp_CreateMolecStruct(title, natom, APPEND);
        if ( Wstr < 0 )
            return(1);
    }

/*  get pointer to coordinate vectors */ 
    PlaceAMBERAFrame(File_p , alt , Wstr);

    return(0);
/*                                  */
}

#define AMBER_LINE_LEN   120   /* AMBER coord file line length */

/*********************************************************************/
int Look4AMBERAFrames(FILE *chm_in , int *Atoms , int *Frames)
/*********************************************************************/
{
    static int i;
    static int AMBERAtoms;
    static int AMBERFrames;

    float   TXc;

    char inputl[AMBER_LINE_LEN];
    int  Trigger;


    rewind(chm_in);
    *Atoms      = 0;
    *Frames     = 0;
    AMBERAtoms = gomp_GetNumAtomsInMolecStruct(0);

    if(AmberaFrameStartPoint != NULL) {
        free(AmberaFrameStartPoint);
        AmberaFrameStartPoint  = NULL;
    }
/*
  Start reading file

*/

/* read title line */
    fgets(inputl,AMBER_LINE_LEN,chm_in);
    gomp_PrintMessage("AMBER file title:");
    gomp_PrintMessage(inputl);

    AMBERFrames = 0;

    while(!feof(chm_in)) {

/* save the pointers in the AMBER file */
        if(AMBERFrames < 1) {
            AmberaFrameStartPoint     = gomp_AllocateIntVector(1);
            AmberaFrameStartPoint[0]  = ftell(chm_in);
        } else {
            AmberaFrameStartPoint = 
                gomp_ReallocateIntVector(AmberaFrameStartPoint , AMBERFrames + 1);
            AmberaFrameStartPoint[AMBERFrames]  = ftell(chm_in);
        }

/* read coordinates (3 * AMBERAtoms) */
        Trigger = 0;
        for (i = 0 ; i < 3 * AMBERAtoms ; i++) {
            fscanf(chm_in,"%f",&TXc);
            if(feof(chm_in)) {
                Trigger = 1;
                break;
            }
        }

        if(Trigger) break;
/* read box size */
        fscanf(chm_in,"%f",&TXc);
        fscanf(chm_in,"%f",&TXc);
        fscanf(chm_in,"%f",&TXc);
            
        AMBERFrames++;
    }

/*   Numbers of atoms */

    *Atoms   = AMBERAtoms;
    *Frames  = AMBERFrames;

    rewind(chm_in);

    return(0);
}

/*********************************************************************/
int PlaceAMBERAFrame(FILE *chm_in , int Alt , int Wstr)
/*********************************************************************/
{
    static int i,j;
    static int AMBERAtoms;
    static int AMBERFrames;

    float   TXc,TYc,TZc;
    const float *sumxyz;

    char inputl[AMBER_LINE_LEN];
    int  items;
    int  method;

/*
  Start reading file

*/
    rewind(chm_in);
    sumxyz         = gomp_GetTranslateArray();
    AMBERFrames    = 0;
    items          = Alt - 1;
    AMBERAtoms     = gomp_GetNumAtomsInMolecStruct(0);

    method = gomp_GetFormattedTrajectoryReader();

    if(method) {

/* read title line */
        fgets(inputl,AMBER_LINE_LEN,chm_in);

/* position -1 from the right one */
        for(j = 0 ; j < items ; j++) {

/* read coordinates (3 * AMBERAtoms) */
            for (i = 0 ; i < 3 * AMBERAtoms ; i++) {
                fscanf(chm_in,"%f",&TXc);
            }
/* read box size */
            fscanf(chm_in,"%f",&TXc);
            fscanf(chm_in,"%f",&TXc);
            fscanf(chm_in,"%f",&TXc);
        }

/* grab now the right one */
/* X */
/* Y */
/* Z */
        for (i = 0 ; i < AMBERAtoms ; i++) {
            fscanf(chm_in,"%f %f %f",&TXc,&TYc,&TZc);
            j = gomp_PutAtomXCoord(Wstr , TXc - sumxyz[0] , i);
            j = gomp_PutAtomYCoord(Wstr , TYc - sumxyz[1] , i);
            j = gomp_PutAtomZCoord(Wstr , TZc - sumxyz[2] , i);
        }

    } else {
 
        if(fseek(chm_in , (long)AmberaFrameStartPoint[Alt - 1] , SEEK_SET)) {
            gomp_PrintERROR("can't position at frame for AMBER trajectory");
            rewind(chm_in);
            return(1);
        }

/* grab now the right one */
/* X */
/* Y */
/* Z */
        for (i = 0 ; i < AMBERAtoms ; i++) {
            fscanf(chm_in,"%f %f %f",&TXc,&TYc,&TZc);
            j = gomp_PutAtomXCoord(Wstr , TXc  - sumxyz[0] , i);
            j = gomp_PutAtomYCoord(Wstr , TYc  - sumxyz[1] , i);
            j = gomp_PutAtomZCoord(Wstr , TZc  - sumxyz[2] , i);
        }
    }


    rewind(chm_in);
    return(0);

}

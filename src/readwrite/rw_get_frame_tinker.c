/*

Copyright (c) 1997 - 2005 by:
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

static int  Look4TINKERFrames(FILE *, int * , int *);
static int  PlaceTINKERFrame(FILE * , int , int );

static int *TinkerFrameStartPoint = (int *)NULL;

/***************************************************************************/
int gomp_GetFrameTINKER(int alt,FILE *File_p , int iappend)  
    /* read one frame from a xmol trajectory   
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

        (void)Look4TINKERFrames(File_p , &natom , &nstep);

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
        sprintf(title,"TINKER frame (%d)",alt);
        Wstr = gomp_CreateMolecStruct(title , natom , APPEND);
        if ( Wstr < 0 )
            return(1);
    }

/*  get pointer to coordinate vectors */ 
    (void)PlaceTINKERFrame(File_p , alt , Wstr);

    return(0);
/*                                  */
}

#define TINKER_LINE_LEN   120   /* TINKER coord file line length */

/*********************************************************************/
int Look4TINKERFrames(FILE *chm_in , int *Atoms , int *Frames)
/*********************************************************************/
{
    static int i;
    static int TinkerAtoms;

    float   TXc,TYc,TZc;
    char    TAtmN[BUFF_LEN];

    char inputl[TINKER_LINE_LEN];
    int  items;


    rewind(chm_in);
    *Atoms  = 0;
    *Frames = 0;

    free(TinkerFrameStartPoint);
    TinkerFrameStartPoint  = (int *)NULL;
/*
  Start reading file

*/

/* Save internal pointer in Tinker file */
    TinkerFrameStartPoint     = gomp_AllocateIntVector(1);
    TinkerFrameStartPoint[0]  = ftell(chm_in);

    while(fgets(inputl,TINKER_LINE_LEN,chm_in) != NULL) {

/*   Numbers of atoms */
        sscanf(inputl,"%d",&TinkerAtoms);

        *Atoms = TinkerAtoms;

/*   Read atom cards   */

        for(i = 0 ; i < TinkerAtoms ; i++ ) {

            fgets(inputl,TINKER_LINE_LEN,chm_in);

            items = sscanf(inputl,"%*d %s %f %f %f",
                           TAtmN,&TXc,&TYc,&TZc);

        }

        *Frames = *Frames + 1;

        TinkerFrameStartPoint = 
            gomp_ReallocateIntVector(TinkerFrameStartPoint , *Frames + 1);
        TinkerFrameStartPoint[*Frames]  = ftell(chm_in);

    }

    rewind(chm_in);

    return(0);
}

/*********************************************************************/
int PlaceTINKERFrame(FILE *chm_in , int Alt , int Wstr)
/*********************************************************************/
{
    static int i,j,k;
    static int TinkerAtoms;

    float   TXc,TYc,TZc;
    const float *sumxyz;
    char    TAtmN[BUFF_LEN];

    char inputl[TINKER_LINE_LEN];
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
            fgets(inputl,TINKER_LINE_LEN,chm_in);
            sscanf(inputl,"%d",&TinkerAtoms);


/*   Read atom cards   */

            for(i = 0 ; i < TinkerAtoms ; i++ ) {
                fgets(inputl,TINKER_LINE_LEN,chm_in);

                if(k == (Alt - 1)) {
                    items = sscanf(inputl,"%*d %s %f %f %f",
                                   TAtmN,&TXc,&TYc,&TZc);

                    j = gomp_PutAtomXCoord(Wstr , TXc - sumxyz[0] , i);
                    j = gomp_PutAtomYCoord(Wstr , TYc - sumxyz[1] , i);
                    j = gomp_PutAtomZCoord(Wstr , TZc - sumxyz[2] , i);

                    j = gomp_PutAtomCharge(Wstr ,  0.0 , i);
                }
            }
        }
    } else {

        if(fseek(chm_in , (long)TinkerFrameStartPoint[Alt - 1] , SEEK_SET)) {
            gomp_PrintERROR("can't position at frame for TINKER trajectory");
            rewind(chm_in);
            return(1);
        }

        fgets(inputl,TINKER_LINE_LEN,chm_in);
        sscanf(inputl,"%d",&TinkerAtoms);

        for(i = 0 ; i < TinkerAtoms ; i++ ) {
            fgets(inputl,TINKER_LINE_LEN,chm_in);

            items = sscanf(inputl,"%*d %s %f %f %f",
                           TAtmN,&TXc,&TYc,&TZc);

            j = gomp_PutAtomXCoord(Wstr , TXc - sumxyz[0] , i);
            j = gomp_PutAtomYCoord(Wstr , TYc - sumxyz[1] , i);
            j = gomp_PutAtomZCoord(Wstr , TZc - sumxyz[2] , i);

            j = gomp_PutAtomCharge(Wstr ,  0.0 , i);
        }
    }

    rewind(chm_in);

    return(0);
}

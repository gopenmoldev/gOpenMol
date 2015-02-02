/*

Copyright (c) 1997 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

1999-02-11: Corrected a bug with the frame counter. Number of atoms returned was always
= 0!

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

#include "coord_file.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"

#include "stdafx.h"

#define MOD(a,b)  (((a+1)/b)*b - (a+1))

static int Look4GROMOS96AFrames(FILE *, int * , int *);
static int PlaceGROMOS96AFrame(FILE * , int , int );

static float GROMOS96CoordAmplifier = 10.0;

static int *Gromos96aFrameStartPoint = (int *)NULL;

/***************************************************************************/
int gomp_GetFrameGROMOS96A(int alt,FILE *File_p , int iappend)  
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

        (void)Look4GROMOS96AFrames(File_p , &natom , &nstep);

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
        sprintf(title,"GROMOS96 frame (%d)",alt);
        Wstr = gomp_CreateMolecStruct(title , natom , APPEND);
        if( Wstr < 0 )
            return(1);
    }
/*  get pointer to coordinate vectors */ 
    (void)PlaceGROMOS96AFrame(File_p , alt , Wstr);

    return(0);
/*                                  */
}

#define GROMOS96_LINE_LEN   120   /*GROMOS96 coord file line length */

/*********************************************************************/
int Look4GROMOS96AFrames(FILE *chm_in , int *Atoms , int *Frames)
/*********************************************************************/
{
    static int GROMOS96Atoms;
    static int GROMOS96Frames;

    char inputl[GROMOS96_LINE_LEN];

    rewind(chm_in);
    *Atoms  = 0;
    *Frames = 0;

    free(Gromos96aFrameStartPoint);
    Gromos96aFrameStartPoint  = NULL;
/*
  Start reading file

*/

    GROMOS96Frames = 0;

    while(fgets(inputl,GROMOS96_LINE_LEN,chm_in) != NULL) {

        GROMOS96Atoms  = 0;

        if(Tcl_StringMatch(inputl,"POSITION*") ||
           Tcl_StringMatch(inputl,"POSITIONRED*")) {

/* save the pointers in the GROMSO96 file */
            if(GROMOS96Frames < 1) {
                Gromos96aFrameStartPoint     = gomp_AllocateIntVector(1);
                Gromos96aFrameStartPoint[0]  = ftell(chm_in);
            } else {
                Gromos96aFrameStartPoint = 
                    gomp_ReallocateIntVector(Gromos96aFrameStartPoint , GROMOS96Frames + 1);
                Gromos96aFrameStartPoint[GROMOS96Frames]  = ftell(chm_in);
            }

            while(fgets(inputl,GROMOS96_LINE_LEN,chm_in) != NULL) {

                if(inputl[0] == '#' || inputl[1] == '#') 
                    continue;

                if(Tcl_StringMatch(inputl,"END*")) {
                    GROMOS96Frames++;
                    break;
                }

                GROMOS96Atoms++;

            }

/*   Numbers of atoms */
            *Atoms   = GROMOS96Atoms;

        }
    }

    *Frames  = GROMOS96Frames;

    rewind(chm_in);

    return(0);
}

/*********************************************************************/
int PlaceGROMOS96AFrame(FILE *chm_in , int Alt , int Wstr)
/*********************************************************************/
{
    static int j;
    static int GROMOS96Atoms;
    static int GROMOS96Frames;

    float   TXc,TYc,TZc;
    const float *sumxyz;

    char inputl[GROMOS96_LINE_LEN];
    int  items;
    int  method;

/*
  Start reading file

*/
    rewind(chm_in);
    sumxyz         = gomp_GetTranslateArray();
    GROMOS96Frames = 0;
    items          = Alt - 1;

    method = gomp_GetFormattedTrajectoryReader();

    if(method) {
/* position -1 from the right one */
        while(fgets(inputl,GROMOS96_LINE_LEN,chm_in) != NULL) {

            GROMOS96Atoms  = 0;

            if(Tcl_StringMatch(inputl,"POSITION*") ||
               Tcl_StringMatch(inputl,"POSITIONRED*")) {

                while(fgets(inputl,GROMOS96_LINE_LEN,chm_in) != NULL) {

                    if(inputl[0] == '#' || inputl[1] == '#') 
                        continue;

                    if(Tcl_StringMatch(inputl,"END*")) {
                        GROMOS96Frames++;
                        break;
                    }

                    GROMOS96Atoms++;

                }
            }
            if(items == GROMOS96Frames) break;
        }

/* grab now the right one */

        while(fgets(inputl,GROMOS96_LINE_LEN,chm_in) != NULL) {

            GROMOS96Atoms  = 0;

            if(Tcl_StringMatch(inputl,"POSITION*") ||
               Tcl_StringMatch(inputl,"POSITIONRED*")) {

                while(fgets(inputl,GROMOS96_LINE_LEN,chm_in) != NULL) {

                    if(inputl[0] == '#' || inputl[1] == '#') 
                        continue;

                    if(Tcl_StringMatch(inputl,"END*")) {
                        rewind(chm_in);
                        return(0);
                    }

                    items = sscanf(inputl,"%f %f %f",
                                   &TXc,&TYc,&TZc);

                    j = gomp_PutAtomXCoord(Wstr , GROMOS96CoordAmplifier * TXc 
                                  - sumxyz[0] , GROMOS96Atoms);
                    j = gomp_PutAtomYCoord(Wstr , GROMOS96CoordAmplifier * TYc 
                                  - sumxyz[1] , GROMOS96Atoms);
                    j = gomp_PutAtomZCoord(Wstr , GROMOS96CoordAmplifier * TZc 
                                  - sumxyz[2] , GROMOS96Atoms);

                    GROMOS96Atoms++;
                }

            }
        }
    } else {
 
        if(fseek(chm_in , (long)Gromos96aFrameStartPoint[Alt - 1] , SEEK_SET)) {
            gomp_PrintERROR("can't position at frame for GROMOS96 trajectory");
            rewind(chm_in);
            return(1);
        }

        GROMOS96Atoms  = 0;

        while(fgets(inputl,GROMOS96_LINE_LEN,chm_in) != NULL) {

            if(inputl[0] == '#' || inputl[1] == '#') 
                continue;

            if(Tcl_StringMatch(inputl,"END*")) {
                rewind(chm_in);
                return(0);
            }

            items = sscanf(inputl,"%f %f %f",
                           &TXc,&TYc,&TZc);

            j = gomp_PutAtomXCoord(Wstr , GROMOS96CoordAmplifier * TXc 
                          - sumxyz[0] , GROMOS96Atoms);
            j = gomp_PutAtomYCoord(Wstr , GROMOS96CoordAmplifier * TYc 
                          - sumxyz[1] , GROMOS96Atoms);
            j = gomp_PutAtomZCoord(Wstr , GROMOS96CoordAmplifier * TZc 
                          - sumxyz[2] , GROMOS96Atoms);

            GROMOS96Atoms++;
        }
    }


    rewind(chm_in);
    return(0);

}
/*********************************************************************/
int   gomp_SetGROMOS96CoordAmplifier(float Value)
/*********************************************************************/
{
    GROMOS96CoordAmplifier = Value;

    return(0);
}
/*********************************************************************/
float gomp_GetGROMOS96CoordAmplifier()
/*********************************************************************/
{
    return(GROMOS96CoordAmplifier);
}

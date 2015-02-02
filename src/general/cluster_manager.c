/*

Copyright (c) 1990 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "cluster.h"
#include "measure.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "selection.h"
#include "trajectory.h"

#include "stdafx.h"

#define ON  1
#define OFF 0

static struct {
    int    Display;
    int    NumSets;
    float *DataArray;
} ClusterData = { 0 , 0 , (float *)NULL};

static struct {
    float min1;
    float max1;
    float min2;
    float max2;
    float min3;
    float max3;
    float min4;
    float max4;
    float rest;
} ClusterBounds = { 0.4f , 0.8f , 0.8f , 1.2f , 1.2f , 1.6f , 1.6f ,
                                   2.0f , 2.0f};
static int GetClusterSpace(int);
static int PutClusterData(int , int , const float *);

#ifdef ENABLE_GRAPHICS
/***********************************************************************/
int gomp_PreCluster()
/***********************************************************************/
{
    char TempS[BUFF_LEN];
    static float min1,min2,min3,min4;
    static float max1,max2,max3,max4;
    static float rest;


    min1  = ClusterBounds.min1;
    max1 = ClusterBounds.max1;
    min2  = ClusterBounds.min2;
    max2 = ClusterBounds.max2;
    min3  = ClusterBounds.min3;
    max3 = ClusterBounds.max3;
    min4  = ClusterBounds.min4;
    max4 = ClusterBounds.max4;
    rest  = ClusterBounds.rest;

/*
  gomp_PrintMessage("The displayed intervals are:");
  sprintf(TempS,"          dist <  %4.2f \n",min1);
  gomp_PrintMessage(TempS);
  sprintf(TempS,"%4.2f   < dist <  %4.2f \n",min1,max1);
  gomp_PrintMessage(TempS);
  sprintf(TempS,"%4.2f   < dist <  %4.2f \n",min2,max2);
  gomp_PrintMessage(TempS);
  sprintf(TempS,"%4.2f   < dist <  %4.2f \n",min3,max3);
  gomp_PrintMessage(TempS);
  sprintf(TempS,"%4.2f   < dist <  %4.2f \n",min4,max4);
  gomp_PrintMessage(TempS);
  sprintf(TempS,"%4.2f   < dist          \n\n",rest);
  gomp_PrintMessage(TempS);
*/

    sprintf(TempS,"(1) Number of frames is: %d",ClusterData.NumSets);
    gomp_PrintMessage(TempS);

    if(ClusterData.NumSets < 1) {
        gomp_PrintMessage("?WARNING - no frames to be displayed ");
        return(1);
    }

    return(gomp_PlotClusterMatrix(
               ClusterData.NumSets , 
               ClusterData.NumSets ,
               min1,min2,min3,min4,max1,max2,max3,max4,rest));


}
#endif /* ENABLE_GRAPHICS */
#if 0
/***********************************************************************/
int gomp_ChangeClusterBounds( const char *text1 , const char *text2 , const char *text3 ,
                            const char *text4 , const char *text5 , const char *text6 ,
                            const char *text7 , const char *text8 , const char *text9)
/***********************************************************************/
{

    ClusterBounds.min1 = atof(text1);
    ClusterBounds.max1 = atof(text2);
    ClusterBounds.min2 = atof(text3);
    ClusterBounds.max2 = atof(text4);
    ClusterBounds.min3 = atof(text5);
    ClusterBounds.max3 = atof(text6);
    ClusterBounds.min4 = atof(text7);
    ClusterBounds.max4 = atof(text8);
    ClusterBounds.rest = atof(text9);

    return(0);
}

/***********************************************************************/
int gomp_ResetClusterBounds()
/***********************************************************************/
{

    ClusterBounds.min1 = 0.4f;
    ClusterBounds.max1 = 0.8f;
    ClusterBounds.min2 = 0.8f;
    ClusterBounds.max2 = 1.2f;
    ClusterBounds.min3 = 1.2f;
    ClusterBounds.max3 = 1.6f;
    ClusterBounds.min4 = 1.6f;
    ClusterBounds.max4 = 2.0f;
    ClusterBounds.rest = 2.0f;

    return(0);      
}
#endif


/************************************************************************/
int gomp_CalcCluster(const char *text1 , const char *text2 , const char *text3,
                   const char *text4 , const char *text5 , const char *text6)
/************************************************************************/
{
    int    i,j,k;
    int    iLoop,jLoop;
    float  Temp;
    char   OutText[BUFF_LEN];
    float *x;
    float *y;
    float *z;
    float *xt;
    float *yt;
    float *zt;
    FILE  *File_p;
    int    FirstFrame;
    int    LastFrame;
    int    StepFrame;
    int    Wstr;
    int *sel_list1;
    int *sel_list2;
    int    slong1,slong2,atom_max;
    float  temp1,temp2,temp3,mmass;    

    if(gomp_GetNumberOfFrames() < 1) {
        gomp_PrintERROR("Number of frames is not defined ");
        return(1);
    }

    if(gomp_GetTrajectoryFileType() == XMOL_TRAJ      ||
       gomp_GetTrajectoryFileType() == GROMOS96A_TRAJ ||
       gomp_GetTrajectoryFileType() == TINKER_TRAJ    ||
       gomp_GetTrajectoryFileType() == FAMBER_TRAJ    ||
       gomp_GetTrajectoryFileType() == FDL_POLY_TRAJ) {
        File_p = fopen(gomp_GetTrajectoryFileName(),"r");
    } else {
        File_p = fopen(gomp_GetTrajectoryFileName(),"rb");
    }

    if(File_p    == NULL) {
        sprintf(OutText,"?ERROR - can't open trajectory file : %s\n",
                gomp_GetTrajectoryFileName());
        gomp_PrintMessage(OutText);
        return(1);
    }

    if(text4[0] == '\0')
        text4 = text1;
    if(text5[0] == '\0')
        text5 = text2;
    if(text6[0] == '\0')
        text6 = text3;

    sprintf(OutText,"Selection1: '%s' '%s' '%s' . Selection2: '%s' '%s' '%s'",
            text1,text2,text3,text4,text5,text6);
    gomp_PrintMessage(OutText);

    Wstr     = 0;
    atom_max = gomp_GetNumAtomsInMolecStruct(0);

    sel_list1 = gomp_AllocateIntVector(atom_max);
    sel_list2 = gomp_AllocateIntVector(atom_max);

    temp1 = 0.0;
    temp2 = 0.0;
    temp3 = 0.0;
    mmass = 0.0;

    if(text1[0] == '\0') {

        for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
            sel_list1[i] = i;
            sel_list2[i] = i;
        }

        slong1 = gomp_GetNumAtomsInMolecStruct(0);
        slong2 = gomp_GetNumAtomsInMolecStruct(1);

        if(slong1 != slong2) {
            gomp_PrintERROR("number of atoms in list 1 and 2 must be equal");
            free(sel_list1);
            free(sel_list2);
            return(1);
        }
    }
    else {
/* select from the first structure (cheat the selection list) */
        slong1 = gomp_MakeSelectionList(0, text1,text2,text3,sel_list1);
        slong2 = gomp_MakeSelectionList(0, text4,text5,text6,sel_list2);
    }

    if(slong1 < 1) {
        gomp_PrintERROR("no atoms selected from list 1");
        free(sel_list1);
        free(sel_list2);
        return(1);
    }

    if(slong2 < 1) {
        gomp_PrintERROR("no atoms selected from list 2");
        free(sel_list1);
        free(sel_list2);
        return(1);
    }

    if(slong1 != slong2) {
        gomp_PrintERROR("number of atoms in list 1 and 2 must be equal");
        free(sel_list1);
        free(sel_list2);
        return(1);
    }

    (void)gomp_GetTrajectoryDisplayParams(&FirstFrame ,
                                        &LastFrame  ,
                                        &StepFrame);

    if(ClusterData.NumSets) /* delete old data */
        (void)gomp_DeleteClusterData();

    ClusterData.NumSets = gomp_GetNumberOfFrames();

/*    SaveCoordinates  */
    (void)gomp_SaveAtomCoords(Wstr);
/*    done .........   */

    xt  = gomp_AllocateFloatVector(gomp_GetNumAtomsInMolecStruct(Wstr));
    yt  = gomp_AllocateFloatVector(gomp_GetNumAtomsInMolecStruct(Wstr));
    zt  = gomp_AllocateFloatVector(gomp_GetNumAtomsInMolecStruct(Wstr));

    (void)GetClusterSpace(ClusterData.NumSets);

    Wstr = 0;
    x    = gomp_GetModifiableAtomXCoordPointer(Wstr);
    y    = gomp_GetModifiableAtomYCoordPointer(Wstr);
    z    = gomp_GetModifiableAtomZCoordPointer(Wstr);

    iLoop = 0;

/* first loop */
    for(i  =  FirstFrame; 
        i <=  LastFrame - 1;
        i += StepFrame) {

        if(gomp_GetOneFrame(i , File_p , TRAJ_OLD)) {
            gomp_PrintERROR("?ERROR - can't read first frame");
            fclose(File_p);
            free(xt);
            free(yt);
            free(zt);
            free(sel_list1);
            free(sel_list2);
            return(1);
        }

        rewind(File_p);

/* put array to reference array */
        for(k = 0; k < gomp_GetNumAtomsInMolecStruct(Wstr) ; k++) {
            xt[k] = x[k];
            yt[k] = y[k];
            zt[k] = z[k];
        }

        jLoop = iLoop + StepFrame;
/* second loop */
        for(j  =  i + StepFrame;
            j <=  LastFrame;
            j += StepFrame) {

            if(gomp_GetOneFrame(j , File_p , TRAJ_OLD)) {
                gomp_PrintERROR("?ERROR - can't read second frame");
                free(xt);
                free(yt);
                free(zt);
                free(sel_list1);
                free(sel_list2);
                fclose(File_p);
                return(1);
            }

            rewind(File_p);

            Temp = gomp_QuatFit(xt ,
                              yt ,
                              zt ,
                              gomp_GetNumAtomsInMolecStruct(Wstr) ,
                              x ,
                              y ,
                              z ,
                              gomp_GetNumAtomsInMolecStruct(Wstr) ,
                              sel_list1 , slong1 ,
                              sel_list2 , slong2 ,
                              OFF);

            if(PutClusterData(iLoop , jLoop , &Temp)) {
                gomp_PrintMessage("?ERROR - can't save cluster data");
                return(1);
            }

            jLoop++;
        }
        iLoop++;
    }

    fclose(File_p);
    free(xt);
    free(yt);
    free(zt);
    free(sel_list1);
    free(sel_list2);

    (void)gomp_GetSavedAtomCoords(Wstr);

    return(0);
}
/************************************************************************/
int gomp_DeleteClusterData()
/************************************************************************/
{
    if(ClusterData.NumSets) 
        free(ClusterData.DataArray);
    ClusterData.NumSets = 0;

    return(0);
}
/************************************************************************/
int GetClusterSpace(int Observ)
/************************************************************************/
{
    if(ClusterData.NumSets) {
        (void)gomp_DeleteClusterData();
    }

    ClusterData.DataArray = gomp_AllocateFloatVector(Observ * (Observ - 1) / 2);

    ClusterData.NumSets   = Observ;

    return(0);
}
/************************************************************************/
int PutClusterData(int ip , int jp , const float *Value)
/************************************************************************/
{
    int i;

    if(!ClusterData.NumSets) {
        gomp_PrintMessage("?ERROR - no space to put Cluster Data");
        return(1);
    }

    i = ip + jp * (jp - 1)/2 ;

    ClusterData.DataArray[i] = *Value;

    return(0);
}

/*************************************************************************/
int gomp_WriteClusterData(const char *FileName , const char *FileTitle)
/*************************************************************************/
{

    FILE *File_p;
    char OutText[BUFF_LEN];
    int ContrlR;

    if(!ClusterData.NumSets) {
        gomp_PrintMessage("?ERROR - no cluster data available");
        return(1);
    }

    if(FileName[0] == '\0') {
        gomp_PrintMessage("?ERROR - file name missing");
        return(1);
    }


    File_p = fopen(FileName,"wb");
    if(File_p == NULL) {
        sprintf(OutText,"?ERROR - can't open file '%s'",FileName);
        gomp_PrintMessage(OutText);
        return(1);
    }   

/* title write */
    ContrlR = fwrite(FileTitle,sizeof(char),80,File_p);
/* size record  */
    ContrlR = fwrite(&ClusterData.NumSets,sizeof(int),1,File_p);
/* data         */
    ContrlR = fwrite(ClusterData.DataArray,sizeof(float),
                     ClusterData.NumSets * (ClusterData.NumSets - 1) /2 ,File_p);

    fclose(File_p);
    return(0);
}

/*************************************************************************/
int gomp_ReadClusterData(const char *FileName)
/*************************************************************************/
{

    FILE *File_p;
    char OutText[BUFF_LEN];
    int ContrlR;
    int CSize;
    int TotalSize;
    int i;
    float MaxC;
    float MinC;

    if(FileName[0] == '\0') {
        gomp_PrintERROR("cluster file name missing");
        return(1);
    }


    File_p = fopen(FileName,"rb");
    if(File_p == NULL) {
        sprintf(OutText,"?ERROR - can't open file '%s'",FileName);
        gomp_PrintMessage(OutText);
        return(1);
    }   

    if(ClusterData.NumSets) 
        (void)gomp_DeleteClusterData();

/* title record */
    ContrlR = fread(OutText,sizeof(char),80,File_p);
    gomp_PrintMessage(OutText);
/* size record  */
    ContrlR = fread(&CSize,sizeof(int),1,File_p);
    ClusterData.NumSets = CSize;
/* data         */
    (void)GetClusterSpace(CSize);
    TotalSize = CSize * (CSize - 1) / 2;
    ContrlR = fread(ClusterData.DataArray,sizeof(float),TotalSize,File_p);
    fclose(File_p);

    MaxC = 0.0;
    MinC = 1.e+10;

    for(i = 0 ; i < TotalSize ; i++) {
        if(ClusterData.DataArray[i] < MinC) MinC = ClusterData.DataArray[i];
        if(ClusterData.DataArray[i] > MaxC) MaxC = ClusterData.DataArray[i];
    }

    sprintf(OutText,"Min cluster value: %f",MinC);
    gomp_PrintMessage(OutText);
    sprintf(OutText,"Max cluster value: %f",MaxC);
    gomp_PrintMessage(OutText);

    return(0);
}
#if 0
/****************************************************************************/
int gomp_ClusterAddNames()
/****************************************************************************/
{

    int   Wstr0 = 0;
    int   Wstr1 = 1;
     

    int   i;

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr0); i++) {

        (void)gomp_PutAtomSegName(Wstr1,gomp_GetAtomSegName(Wstr0,i),i);
        (void)gomp_PutAtomResName(Wstr1,gomp_GetAtomResName(Wstr0,i),i);
        (void)gomp_PutAtomAtmName(Wstr1,gomp_GetAtomAtmName(Wstr0,i),i);

        (void)gomp_PutAtomResNum1(Wstr1,gomp_GetAtomResNum1(Wstr0,i),i);
        (void)gomp_PutAtomResNum2(Wstr1,gomp_GetAtomResNum2(Wstr0,i),i);

    }

    return(0);
}
#endif
/***********************************************************************/
const float *gomp_GetClusterData()
/***********************************************************************/
{
    if(ClusterData.NumSets) 
        return(ClusterData.DataArray);
    else
        return((const float *)NULL);
}

/***********************************************************************/
int    gomp_GetClusterStatus()
/***********************************************************************/
{
    return(ClusterData.NumSets);
}

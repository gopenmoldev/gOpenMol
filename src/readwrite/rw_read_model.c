/*

Copyright (c) 1991 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "contour.h"
#include "gomfile.h"
#include "gomtext.h"
#include "model_file.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "objseg.h"
#include "plumber.h"
#include "printmsg.h"
#include "projview.h"
#include "selection.h"

#include "stdafx.h"

int gomp_ModelFileVersion;

/*********************************************************************/
int gomp_ReadOldModel(const char *FileName)
/*********************************************************************/
{

    int   Version = 0;  /* in fact it is 100 * version in float) */
    char  InputText[BUFF_LEN];
    char  InputText1[BUFF_LEN];
    char  OutText[BUFF_LEN];
    int   NumStruct;
    int   NumAtoms;
    int   Wstr = 0;
    int   i,j,k;

    int   ResN1;
    int   ResN2;
    char  SegNam[BUFF_LEN];
    char  ResNam[BUFF_LEN];
    char  AtmNam[BUFF_LEN];
    float XC;
    float YC;
    float ZC;
    float BV;
    float AtmCha;
    float AtmNucCha;   
    float AtmCovar;
    float XMove;
    float YMove;
    float ZMove;

    int   AtomType;
    float AtomBndRad;
    float AtomVdwRad;
    float AtomPluRad;
    char  AtomGlobal;
    float AtomEmin;
    float AtomRmin;
    float AtomPatom;
    float AtomMass;
    int   AtomCnct;
    char  AtomHbond;
    char  AtomAtype[BUFF_LEN];

    char *DispList;
    int   DL;
    char *CPKList;
    float VDWr;
    int   CL;
    char *LicoList;
    int   LL;

    int   AtmC[BUFF_LEN];
    int   AtmCN;
    int *AtmConn;

    float MSize;
    float Mnear;
    float Mfar;
    float Mdamp;

    float  Matrix[16];
    FILE  *Model_f;

    Model_f = fopen(FileName,"r");
    if(Model_f == NULL) {
        sprintf(OutText,"Can't open input file : %s",FileName);
        gomp_PrintMessage(OutText);
        return(1);
    }

/* VERSION - tag */
    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    if(!strncmp(InputText,"[Version]",strlen("[Version]"))) {
        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        Version = atoi(InputText);
        gomp_ModelFileVersion = Version;
    }
/* TAG - tag */
    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    if(!strncmp(InputText,"[Tag]",strlen("[Tag]"))) {
        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        (void)gomp_PutTagText(InputText);
    }

/* DESCRIPTION - tag */
    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    if(!strncmp(InputText,"[Description]",strlen("[Description]"))) {
        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        (void)gomp_PutDescText(InputText);
    }

/* AVAILABLE - tag */
    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    if(!strncmp(InputText,"[Available]",strlen("[Available]"))) {
        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        (void)gomp_PutAvailText(InputText);
    }


/* NUMSYST - tag   */
    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    if(!strncmp(InputText,"[NumSyst]",strlen("[NumSyst]"))) {
        if(gomp_Fgets(InputText,BUFF_LEN,Model_f) != NULL) 
            sscanf(InputText,"%d",&NumStruct);
        else {
            gomp_PrintERROR("Incomplete 'gom'-file");
            return(1);
        }
    }

/* loop over the systems */

    for(i = 0 ; i < NumStruct ; i++) {

/* from version 2 there is a Name Tag ... */
        if(Version > 1) {
            (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
            if(strncmp(InputText,"[Name]",strlen("[Atoms]"))) {
                gomp_PrintERROR("Problems in [Name] - tag");
                return(1);
            }

            (void)gomp_Fgets(InputText1,BUFF_LEN,Model_f);
        }

        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        if(strncmp(InputText,"[Atoms]",strlen("[Atoms]"))) {
            gomp_PrintERROR("Problems in [Atoms] - tag");
            return(1);
        }

        if(Version > 5) {
            (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
            sscanf(InputText,"%d %f %f %f %d %d",&NumAtoms , &XMove, &YMove, &ZMove, &j ,&k);
            (void)gomp_SetSystemTranslateState(j);
            (void)gomp_SetObjectCenterType(k);
            if(k) 
                (void)gomp_SetLocalTransformationState();
            else
                (void)gomp_SetGlobalTransformationState();

            if(i)
                Wstr = gomp_CreateMolecStruct(InputText1 , NumAtoms , APPEND); 
            else
                Wstr = gomp_CreateMolecStruct(InputText1 , NumAtoms , NEW);
            if ( Wstr < 0 )
                goto end;
       
            (void)gomp_SaveTranslateArray(XMove , YMove , ZMove);

        } else {
            (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
            sscanf(InputText,"%d %f %f %f",&NumAtoms , &XMove, &YMove, &ZMove);
            (void)gomp_SaveTranslateArray(XMove , YMove , ZMove);
        }

        if(Version < 2) {
            sprintf(InputText1,"Structure nr: %d",(i+1));
        }

        (void)gomp_ActivateSelectedStructure( i , STRUCTURE_SELECTION_ON);

        for(j = 0 ; j < NumAtoms ; j++) {
            (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);

            sscanf(InputText," %*d %d %s %s %f %f %f %s %d %f %f %f %f\n",
                   &ResN1,ResNam,AtmNam,&XC,&YC,&ZC,SegNam,&ResN2,&BV,
                   &AtmCha,&AtmNucCha,&AtmCovar);

            gomp_PutAtomResNum1( i , ResN1 , j);
            gomp_PutAtomResNum2( i , ResN2 , j);
            gomp_PutAtomSegName( i , SegNam ,j); 
            gomp_PutAtomResName( i , ResNam ,j);
            gomp_PutAtomAtmName( i , AtmNam ,j);
            gomp_PutAtomXCoord(i , (XC - XMove) , j); 
            gomp_PutAtomYCoord(i , (YC - YMove) , j);
            gomp_PutAtomZCoord(i , (ZC - ZMove) , j);
            gomp_PutAtomBValue( i , BV , j);
            gomp_PutAtomCharge( i , AtmCha ,j); 
            gomp_PutAtomNucCharge( i , AtmNucCha , j);
            gomp_PutAtomCovar( i , AtmCovar , j);
        }

        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        if(strncmp(InputText,"[BasisSetTags]",strlen("[BasisSetTags]")))
            return(1);

        for(j = 0 ; j < NumAtoms ; j++) {
            (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
            gomp_PutAtomBasisSetTag( i , InputText , j);
        }

        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        if(strncmp(InputText,"[AtomParameters]",strlen("[AtomParameters]")))
            return(1);

        for(j = 0 ; j < NumAtoms ; j++) {
            (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
            sscanf(InputText,"%d %f %f %f %c %f %f %f %f %d %c %s\n",
                   &AtomType,
                   &AtomBndRad,
                   &AtomVdwRad,
                   &AtomPluRad,
                   &AtomGlobal,
                   &AtomEmin,
                   &AtomRmin,
                   &AtomPatom,
                   &AtomMass,
                   &AtomCnct,
                   &AtomHbond,
                   AtomAtype);

            gomp_PutAtomType(i ,   AtomType   , j);
            gomp_PutAtomBndRad(i , AtomBndRad , j);
            gomp_PutAtomVdwRad(i , AtomVdwRad , j);
            gomp_PutAtomPluRad(i , AtomPluRad , j);
            gomp_PutAtomGlobal(i , AtomGlobal , j);
            gomp_PutAtomEmin(i ,   AtomEmin   , j);
            gomp_PutAtomRmin(i ,   AtomRmin   , j);
            gomp_PutAtomPatom(i ,  AtomPatom  , j);
            gomp_PutAtomMass(i ,   AtomMass   , j);
            gomp_PutAtomCnct(i ,   AtomCnct   , j);
            gomp_PutAtomHbond(i ,  AtomHbond  , j);
            gomp_PutAtomAtype(i ,  AtomAtype  , j);

        }

        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        if(strncmp(InputText,"[AtomColour]",strlen("[AtomColour]")))
            return(1);

        for(j = 0 ; j < NumAtoms ; j++) {
            float rr, gg, bb;
            (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
            sscanf(InputText,"%f %f %f",&rr,&gg,&bb);
            gomp_PutAtomColour(i, rr, gg, bb, j);
        }

        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        if(strncmp(InputText,"[DisplayLists]",strlen("[DisplayLists]")))
            return(1);

        DispList = gomp_GetModifiableAtomDisplayStatePointer( i );
        CPKList  = gomp_GetModifiableAtomCPKDisplayStatePointer( i );
        LicoList = gomp_GetModifiableAtomLicoDisplayStatePointer( i );

        for(j = 0 ; j < NumAtoms ; j++) {
            (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
            sscanf(InputText,"%d %d %f %d",&DL,&CL,&VDWr,&LL);

            DispList[j] = DL;
            CPKList[j]  = CL;
            (void)gomp_PutAtomVdwRad(i , VDWr , j);
            LicoList[j] = LL;
        }

        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        if(strncmp(InputText,"[AtomConnectivity]",strlen("[AtomConnectivity]")))
            return(1);

        for(j = 0 ; j < NumAtoms ; j++) {

            AtmConn = gomp_GetModifiableAtomConnection(i, j);

            (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
            sscanf(InputText,"%*d %d %d %d %d %d %d %d %d %d %d %d",
                   &AtmCN,
                   &AtmC[0],
                   &AtmC[1],
                   &AtmC[2],
                   &AtmC[3],
                   &AtmC[4],
                   &AtmC[5],
                   &AtmC[6],
                   &AtmC[7],
                   &AtmC[8],
                   &AtmC[9]);

            AtmConn[0] = AtmCN;
            for(k = 1 ; k <= AtmCN ; k++)
                AtmConn[k] = AtmC[k-1];
        }

        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        if(strncmp(InputText,"[Display]",strlen("[Display]")))
            return(1);

        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText,"%f %f %f %f",&MSize,&Mnear,&Mfar,&Mdamp);

        (void)gomp_SetSizeOfSystem(MSize);
        (void)gomp_SetPerspectiveNear(Mnear);
        (void)gomp_SetPerspectiveFar(Mfar);
        (void)gomp_SetTranslationDamping(Mdamp);

        if(Version > 7) {
            (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
            sscanf(InputText,"%f %f",&Mnear,&Mfar);
            (void)gomp_PutAtomLicoRadS(i , Mnear );
            (void)gomp_PutAtomLicoRadC(i , Mnear);
        }

/* from Version 2 up there is some extra information */
        if(Version > 2)
            (void)gomp_ReadDisplayAttributesGOM(Model_f);

        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText,"%f %f %f",&XMove, &YMove, &ZMove);
        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText,"%f %f %f %f",&Matrix[0],
               &Matrix[1],
               &Matrix[2],
               &Matrix[3]);
        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText,"%f %f %f %f",&Matrix[4],
               &Matrix[5],
               &Matrix[6],
               &Matrix[7]);
        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText,"%f %f %f %f",&Matrix[8],
               &Matrix[9],
               &Matrix[10],
               &Matrix[11]);
        (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText,"%f %f %f %f",&Matrix[12],
               &Matrix[13],
               &Matrix[14],
               &Matrix[15]);

        (void)gomp_SaveTranslateArrayMT(    i , XMove , YMove , ZMove);
        (void)gomp_SaveModelViewMatrixMT( i , Matrix);

    }

    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText,"%f %f %f %f",&Matrix[0],
           &Matrix[1],
           &Matrix[2],
           &Matrix[3]);
    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText,"%f %f %f %f",&Matrix[4],
           &Matrix[5],
           &Matrix[6],
           &Matrix[7]);
    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText,"%f %f %f %f",&Matrix[8],
           &Matrix[9],
           &Matrix[10],
           &Matrix[11]);
    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText,"%f %f %f %f",&Matrix[12],
           &Matrix[13],
           &Matrix[14],
           &Matrix[15]);
    (void)gomp_PutModelviewMatrix(Matrix);


    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText,"%f %f %f %f",&Matrix[0],
           &Matrix[1],
           &Matrix[2],
           &Matrix[3]);
    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText,"%f %f %f %f",&Matrix[4],
           &Matrix[5],
           &Matrix[6],
           &Matrix[7]);
    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText,"%f %f %f %f",&Matrix[8],
           &Matrix[9],
           &Matrix[10],
           &Matrix[11]);
    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText,"%f %f %f %f",&Matrix[12],
           &Matrix[13],
           &Matrix[14],
           &Matrix[15]);
    (void)gomp_PutProjectionMatrix(Matrix);

/* an infinite loop here to grab the rest of the information */

    while(gomp_Fgets(InputText,BUFF_LEN,Model_f) != NULL) {

/* contour tag */
        if(!strncmp(InputText,"[Contour]",strlen("[Contour]")))
            (void)gomp_ReadContourInfo2ModelFile(Model_f);

/* tcl tag */
        else if(!strncmp(InputText,"[Tcl]",strlen("[Tcl]")))
            (void)gomp_ReadTclInfo2FromFile(Model_f);

/* plumber tag */
        else if(!strncmp(InputText,"[Plumber]",strlen("[Plumber]")))
            (void)gomp_RetrievePlumberInfoFromModelFile(Model_f);

/* sphere tag */
        else if(!strncmp(InputText,"[Sphere]",strlen("[Sphere]")))
            (void)gomp_ReadSphereSegFromModelFile(Model_f);

/* clinder tag */
        else if(!strncmp(InputText,"[Cylinder]",strlen("[Cylinder]")))
            (void)gomp_ReadCylinderSegFromModelFile(Model_f);

/* arrow tag */
        else if(!strncmp(InputText,"[Arrow]",strlen("[Arrow]")))
            (void)gomp_ReadArrowSegFromModelFile(Model_f);

/* plane tag */
        else if(!strncmp(InputText,"[Plane]",strlen("[Plane]")))
            (void)gomp_ReadPlaneSegFromModelFile(Model_f);

/* x cutplane */
        else if(!strncmp(InputText,"[Cut plane X]",strlen("[Cut plane X]")))
            (void)gomp_ReadCutPlaneXFromModelFile(Model_f);

/* y cutplane */
        else if(!strncmp(InputText,"[Cut plane Y]",strlen("[Cut plane Y]")))
            (void)gomp_ReadCutPlaneYFromModelFile(Model_f);

/* z cutplane */
        else if(!strncmp(InputText,"[Cut plane Z]",strlen("[Cut plane Z]")))
            (void)gomp_ReadCutPlaneZFromModelFile(Model_f);

/* xyz cutplane(s) */
        else if(!strncmp(InputText,"[Cut plane XYZ]",strlen("[Cut plane XYZ]")))
            (void)gomp_ReadCutPlaneXYZFromModelFile(Model_f);

/* text stack */
        else if(!strncmp(InputText,"[Text]",strlen("[Text]")))
            (void)gomp_ReadText2ModelFile(Model_f);

/* contour clipplane */
        else if(!strncmp(InputText,"[Contour Clipplane]",strlen("[Contour Clipplane]")))
            (void)gomp_ReadContourClipplaneInfoFromModelFile(Model_f);

    }

/* now all infomation is collected ......................... */

end:
    fclose(Model_f);

    return(Wstr >= 0 ? 0 : 1);
}

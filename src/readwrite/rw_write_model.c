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
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "contour.h"
#include "gomtext.h"
#include "gomversion.h"
#include "model_file.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "objseg.h"
#include "plumber.h"
#include "printmsg.h"
#include "projview.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_WriteOldModel(const char *FileName)
/*********************************************************************/
{

    int   Version = GOPENMOL_MODEL_FILE_VERSION;
    char  OutText[BUFF_LEN];
    int   i,j,k;
    const int *AtmConn;
    const float *Red,*Green,*Blue;
    const float *MatrixP;
    const float *Move;

    FILE *Model_f;

    Model_f = fopen(FileName,"w");
    if(Model_f == NULL) {
        sprintf(OutText,"Can't open output file : '%s'",FileName);
        gomp_PrintMessage(OutText);
        return(1);
    }

/* VERSION - tag */
    fprintf(Model_f , "[Version]\n");
    fprintf(Model_f , "%d\n",Version);

/* TAG - tag */
    fprintf(Model_f , "[Tag]\n");
    fprintf(Model_f , "%s\n",gomp_GetTagText());

/* DESCRIPTION - tag */
    fprintf(Model_f , "[Description]\n");
    fprintf(Model_f , "%s\n",gomp_GetDescText());

/* AVAILABLE - tag */
    fprintf(Model_f , "[Available]\n");
    fprintf(Model_f , "%s\n",gomp_GetAvailText());

/* NUMSYST - tag   */
    fprintf(Model_f , "[NumSyst]\n");
    fprintf(Model_f , "%d\n",gomp_GetNumMolecStructs());

/* loop over the systems */
    Move = gomp_GetTranslateArray();
 
    for(i = 0 ; i < gomp_GetNumMolecStructs() ; i++) {

        fprintf(Model_f , "[Name]\n");
        fprintf(Model_f , "%s\n",gomp_GetMolecStructName(i));

        fprintf(Model_f , "[Atoms]\n");
        fprintf(Model_f , "%d %f %f %f %d %d\n",gomp_GetNumAtomsInMolecStruct(i),
                Move[0] , Move[1] , Move[2]  ,
                gomp_GetSystemTranslateState() ,
                gomp_GetObjectCenterType());

        for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++) {
            fprintf(Model_f," %d %d %.4s %.4s %f %f %f %.4s %d %f %f %f %f\n",
                    j+1,gomp_GetAtomResNum1( i , j),gomp_GetAtomResName( i , j) ,gomp_GetAtomAtmName(i , j),
                    gomp_GetAtomXCoord(i , j) + Move[0], 
                    gomp_GetAtomYCoord(i , j) + Move[1],
                    gomp_GetAtomZCoord(i , j) + Move[2],
                    gomp_GetAtomSegName( i , j), gomp_GetAtomResNum1( i , j),gomp_GetAtomBValue( i , j),
                    gomp_GetAtomCharge( i , j) , gomp_GetAtomNucCharge( i , j) , gomp_GetAtomCovar(i , j));
        }

        fprintf(Model_f , "[BasisSetTags]\n");
        for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++) 
            fprintf(Model_f,"%s\n",gomp_GetAtomBasisSetTag( i , j));

        fprintf(Model_f , "[AtomParameters]\n");
        for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++) 
            fprintf(Model_f,"%d %f %f %f %c %f %f %f %f %d %c %s\n",
                    gomp_GetAtomType(i , j),
                    gomp_GetAtomBndRad(i , j),
                    gomp_GetAtomVdwRad(i , j),
                    gomp_GetAtomPluRad(i , j),
                    gomp_GetAtomGlobal(i , j),
                    gomp_GetAtomEmin(i , j),
                    gomp_GetAtomRmin(i , j),
                    gomp_GetAtomPatom(i , j),
                    gomp_GetAtomMass(i , j),
                    gomp_GetAtomCnct(i , j),
                    gomp_GetAtomHbond(i , j),
                    gomp_GetAtomAtype(i , j));

        fprintf(Model_f , "[AtomColour]\n");
        Red   = gomp_GetAtomColourRedPointer(i);
        Green = gomp_GetAtomColourGreenPointer(i);
        Blue  = gomp_GetAtomColourBluePointer(i);

        for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++) 
            fprintf(Model_f,"%f %f %f\n",Red[j],Green[j],Blue[j]);

        fprintf(Model_f , "[DisplayLists]\n");
        for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++) 
            fprintf(Model_f,"%d %d %f %d\n",
                    gomp_GetAtomDisplayState(i , j),
                    gomp_GetAtomCPKDisplayState( i , j),
                    gomp_GetAtomVdwRad(i , j),
                    gomp_GetAtomLicoDisplayState(i , j));
  
        fprintf(Model_f , "[AtomConnectivity]\n");
        for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++) {
            AtmConn = gomp_GetAtomConnection(i, j);
            fprintf(Model_f,"%d %d",j,AtmConn[0]);
            for(k = 1 ; k <= AtmConn[0] ; k++) 
                fprintf(Model_f," %d",AtmConn[k]);
            fprintf(Model_f,"\n");
        }

        fprintf(Model_f , "[Display]\n");
        fprintf(Model_f,"%f %f %f %f\n",gomp_GetSizeOfSystem(),
                gomp_GetPerspectiveNear(),
                gomp_GetPerspectiveFar(),
                gomp_GetTranslationDamping());

        fprintf(Model_f,"%f %f\n",
                gomp_GetAtomLicoRadS(i),gomp_GetAtomLicoRadC(i));
      

        (void)gomp_WriteDisplayAttributesGOM(Model_f);

        if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {

            MatrixP = gomp_GetSavedModelViewMatrixMT(i);
            Move    = gomp_GetTranslateArrayMT(i);
            fprintf(Model_f,"%f %f %f\n",Move[0],Move[1],Move[2]);
            fprintf(Model_f,"%f %f %f %f\n",
                    MatrixP[0],MatrixP[1],MatrixP[2],MatrixP[3]);
            fprintf(Model_f,"%f %f %f %f\n",
                    MatrixP[4],MatrixP[5],MatrixP[6],MatrixP[7]);
            fprintf(Model_f,"%f %f %f %f\n",
                    MatrixP[8],MatrixP[9],MatrixP[10],MatrixP[11]);
            fprintf(Model_f,"%f %f %f %f\n",
                    MatrixP[12],MatrixP[13],MatrixP[14],MatrixP[15]);
        }
        else {
            fprintf(Model_f,"0.0 0.0 0.0\n");
            fprintf(Model_f,"1.0 0.0 0.0 0.0\n");
            fprintf(Model_f,"0.0 1.0 0.0 0.0\n");
            fprintf(Model_f,"0.0 0.0 1.0 0.0\n");
            fprintf(Model_f,"0.0 0.0 0.0 1.0\n");
        }
    }

#ifdef ENABLE_GRAPHICS
    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        float  Matrix[16];
        (void)gomp_GetModelviewMatrix(Matrix);
        fprintf(Model_f,"%f %f %f %f\n",Matrix[0],Matrix[1],Matrix[2],Matrix[3]);
        fprintf(Model_f,"%f %f %f %f\n",Matrix[4],Matrix[5],Matrix[6],Matrix[7]);
        fprintf(Model_f,"%f %f %f %f\n",Matrix[8],Matrix[9],Matrix[10],Matrix[11]);
        fprintf(Model_f,"%f %f %f %f\n",Matrix[12],Matrix[13],Matrix[14],Matrix[15]);

        (void)gomp_GetProjectionMatrix(Matrix);
        fprintf(Model_f,"%f %f %f %f\n",Matrix[0],Matrix[1],Matrix[2],Matrix[3]);
        fprintf(Model_f,"%f %f %f %f\n",Matrix[4],Matrix[5],Matrix[6],Matrix[7]);
        fprintf(Model_f,"%f %f %f %f\n",Matrix[8],Matrix[9],Matrix[10],Matrix[11]);
        fprintf(Model_f,"%f %f %f %f\n",Matrix[12],Matrix[13],Matrix[14],Matrix[15]);

    }
    else {
#endif /* ENABLE_GRAPHICS */
        fprintf(Model_f,"1.0 0.0 0.0 0.0\n");
        fprintf(Model_f,"0.0 1.0 0.0 0.0\n");
        fprintf(Model_f,"0.0 0.0 1.0 0.0\n");
        fprintf(Model_f,"0.0 0.0 0.0 1.0\n");


        fprintf(Model_f,"1.0 0.0 0.0 0.0\n");
        fprintf(Model_f,"0.0 1.0 0.0 0.0\n");
        fprintf(Model_f,"0.0 0.0 1.0 0.0\n");
        fprintf(Model_f,"0.0 0.0 0.0 1.0\n");
#ifdef ENABLE_GRAPHICS
    }
#endif /* ENABLE_GRAPHICS */

/* from here on everything has to be grabbed by an "infinite loop" */

    (void)gomp_WriteContourInfo2ModelFile(Model_f);

    (void)gomp_StorePlumberInfo2ModelFile(Model_f);

    (void)gomp_WriteSphereSeg2ModelFile(Model_f);
    (void)gomp_WriteCylinderSeg2ModelFile(Model_f);
    (void)gomp_WriteArrowSeg2ModelFile(Model_f);
    (void)gomp_WritePlaneSeg2ModelFile(Model_f);
    (void)gomp_WriteCutPlane2ModelFile(Model_f);
    (void)gomp_WriteText2ModelFile(Model_f);
    (void)gomp_WriteContourClipplaneInfo2ModelFile(Model_f);

    fclose(Model_f);

/* model file is saved, reset trigger ... */
    (void)gomp_SetFileSavingState(0);

    return(0);
}

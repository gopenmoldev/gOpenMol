/*

Copyright (c) 1990 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002, 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "coord_file.h"
#include "memalloc.h"
#include "model_file.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"

#include "stdafx.h"

#define CONV_ATOMS 57

static int PrepareBaSAtoms( char * , int);

/************************************************************************/
int gomp_WriteCoordFree(int Which , const char *file_name , int alt)
/************************************************************************/
{

    FILE *write_p;
    int i;
    int ihelp;
    char OutText[BUFF_LEN];
    const char *DispList;
    const float *Move;

/* check first structures available */
    i = gomp_GetNumMolecStructs();
    if(!i) {
        gomp_PrintERROR("no structure available");
        return(1);
    }
    if(Which < 1 || Which > i) {
        gomp_PrintERROR("structure index out of range");
        return(1);
    }

    if(file_name[0] == '?') {
        gomp_PrintMessage("?ERROR - file name can't contain '?'");
        return(1);
    }

    write_p = fopen(file_name,"w");
    if(write_p == NULL) {
        sprintf(OutText,"?ERROR - unable to open file : %s ",file_name);
        gomp_PrintMessage(OutText);
        return(1);
    }

    sprintf(OutText,"Writing current molecule/system (nr:%d)",Which);
    gomp_PrintMessage(OutText);

    Which    = Which - 1;
    Move     = gomp_GetTranslateArray();
    DispList = gomp_GetAtomDisplayStatePointer(Which);

    ihelp = 0;
    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {

        if(alt) {
            if(DispList[i] == 0) continue;
        }
        ihelp++;
    }

    fprintf(write_p,"* Output in CHARMM free format \n");
    fprintf(write_p,"*        \n");
    fprintf(write_p," %d \n",ihelp);

    ihelp = 1;
    for(i = 0  ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {

        if(alt) {
            if(DispList[i] == 0) continue;
        }
        fprintf(write_p," %d %d %.4s %.4s %f %f %f %.4s %d %f\n",
                ihelp,gomp_GetAtomResNum1( Which , i),gomp_GetAtomResName( Which , i) ,gomp_GetAtomAtmName(Which , i),
                (gomp_GetAtomXCoord(Which , i) + Move[0]), 
                (gomp_GetAtomYCoord(Which , i) + Move[1]),
                (gomp_GetAtomZCoord(Which , i) + Move[2]),
                gomp_GetAtomSegName(Which , i), gomp_GetAtomResNum1(Which , i),gomp_GetAtomBValue(Which , i));
        ihelp++;
    }

    fclose(write_p);

    return(0);
}

/************************************************************************/
int gomp_WriteCoordBaS(int Which , const char *file_name , int alt)
/************************************************************************/
{

    FILE *write_p;
    int i;
    int ihelp;
    char OutText[BUFF_LEN];
    char *element_vec;
    int ret_val;
    const float *Move;

/* check first structures available */
    i = gomp_GetNumMolecStructs();
    if(!i) {
        gomp_PrintERROR("no structure available");
        return(1);
    }
    if(Which < 1 || Which > i) {
        gomp_PrintERROR("structure index out of range");
        return(1);
    }

    if(file_name[0] == '?') {
        gomp_PrintMessage("?ERROR - file name can't contain '?'");
        return(1);
    }

    write_p = fopen(file_name,"w");
    if(write_p == NULL) {
        sprintf(OutText,"?ERROR - unable to open file : %s ",file_name);
        gomp_PrintMessage(OutText);
        return(1);
    }

/* ready to write now ... */

    sprintf(OutText,"Writing current molecule/system (nr:%d)",Which);
    gomp_PrintMessage(OutText);

    Which   = Which - 1;
    element_vec = 
        gomp_AllocateCharVector(2*(gomp_GetNumAtomsInMolecStruct(Which) + 1));

    ihelp = 0;

    Move = gomp_GetTranslateArray();

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {

        if(alt) {
            if(gomp_GetAtomDisplayState( Which , i) == 0)
                continue;
        }
        strncpy(element_vec+2*ihelp,gomp_GetAtomAtmName( Which , i),2);
        ihelp++;
    }

    ret_val = PrepareBaSAtoms(element_vec , ihelp); 

    fprintf(write_p,"Title SCARECROW output to B&S\n");

    ihelp = 0;

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {

        if(alt) {
            if(gomp_GetAtomDisplayState( Which , i) == 0)
                continue;
        }
        fprintf(write_p,"%.2s  %f  %f  %f\n",
                element_vec+2*ihelp, gomp_GetAtomXCoord( Which , i)+Move[0],
                gomp_GetAtomYCoord( Which , i)+Move[1],
                gomp_GetAtomYCoord( Which , i)+Move[2]);
        ihelp++;
    }

    fprintf(write_p,"END\n");

    fclose(write_p);

    free(element_vec);

    return(0);
}
/************************************************************************/
int gomp_WriteCoordPDB(int Which , const char *file_name , int alt)
/************************************************************************/
{

    FILE *write_p;
    int i;
    int ihelp;
    char OutText[BUFF_LEN];
    const float *Move;

/* check first structures available */
    i = gomp_GetNumMolecStructs();
    if(!i) {
        gomp_PrintERROR("no structure available");
        return(1);
    }
    if(Which < 1 || Which > i) {
        gomp_PrintERROR("structure index out of range");
        return(1);
    }

    if(file_name[0] == '?') {
        gomp_PrintMessage("?ERROR - file name can't contain '?'");
        return(1);
    }

    write_p = fopen(file_name,"w");
    if(write_p == NULL) {
        sprintf(OutText,"?ERROR - unable to open file : %s ",file_name);
        gomp_PrintMessage(OutText);
        return(1);
    }

/* ready to write now ... */
    sprintf(OutText,"Writing current molecule/system (nr:%d)",Which);
    gomp_PrintMessage(OutText);

    fprintf(write_p,"REMARK   1 PDB output from SCARECROW\n");
    fprintf(write_p,"REMARK   1 =========================\n");

    Which = Which - 1;
    Move  = gomp_GetTranslateArray();
    ihelp = 1;
    for(i = 0  ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {

        if(alt) {
            if(gomp_GetAtomDisplayState( Which , i) == 0)
                continue;
        }

        strncpy(OutText,gomp_GetAtomAtmName( Which , i),MAX_ATM_NAME_LEN);

        if(isalpha(OutText[0])) {
            fprintf(write_p,
                    "ATOM  %5d  %-3.3s %-3.3s %1.1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    ihelp,gomp_GetAtomAtmName( Which , i),gomp_GetAtomResName( Which , i),gomp_GetAtomSegName( Which , i),
                    gomp_GetAtomResNum1( Which , i),
                    gomp_GetAtomXCoord(Which , i)+Move[0],
                    gomp_GetAtomYCoord(Which , i)+Move[1],
                    gomp_GetAtomZCoord(Which , i)+Move[2],0.0,gomp_GetAtomBValue( Which , i));
        } else {
            fprintf(write_p,
                    "ATOM  %5d %-4.4s %-3.3s %1.1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    ihelp,gomp_GetAtomAtmName( Which , i),gomp_GetAtomResName( Which , i),gomp_GetAtomSegName( Which , i),
                    gomp_GetAtomResNum1( Which , i),
                    gomp_GetAtomXCoord(Which , i)+Move[0],
                    gomp_GetAtomYCoord(Which , i)+Move[1],
                    gomp_GetAtomZCoord(Which , i)+Move[2],0.0,gomp_GetAtomBValue( Which , i));
        }

        ihelp++;
    }

    fclose(write_p);

    return(0);
}


/************************************************************************/
int gomp_WriteCoordCHARMM(int Which , const char *file_name , int alt)
/************************************************************************/
{

    FILE *write_p;
    int i;
    int ihelp;
    char OutText[BUFF_LEN];
    const float *Move;

/* check first structures available */
    i = gomp_GetNumMolecStructs();
    if(!i) {
        gomp_PrintERROR("no structure available");
        return(1);
    }
    if(Which < 1 || Which > i) {
        gomp_PrintERROR("structure index out of range");
        return(1);
    }

    if(file_name[0] == '?') {
        gomp_PrintMessage("?ERROR - file name can't contain '?'");
        return(1);
    }

    write_p = fopen(file_name,"w");
    if(write_p == NULL) {
        sprintf(OutText,"?ERROR - unable to open file : %s ",file_name);
        gomp_PrintMessage(OutText);
        return(1);
    }

    sprintf(OutText,"Writing current molecule/system (nr:%d)",Which);
    gomp_PrintMessage(OutText);

    Which = Which - 1;
    Move  = gomp_GetTranslateArray();
    ihelp = 0;
    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {

        if(alt) {
            if(gomp_GetAtomDisplayState( Which , i) == 0)
                continue;
        }
        ihelp++;
    }

    fprintf(write_p,"* Output in CHARMM format \n");
    fprintf(write_p,"* ======================= \n");
    fprintf(write_p,"*        \n");
    fprintf(write_p,"%5d \n",ihelp);

    ihelp = 1;
    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {

        if(alt) {
            if(gomp_GetAtomDisplayState( Which , i) == 0)
                continue;
        }
        fprintf(write_p,
                "%5d%5d %-4.4s %-4.4s%10.5f%10.5f%10.5f %4.4s %-4d%10.5f\n",
                ihelp,gomp_GetAtomResNum1( Which , i),gomp_GetAtomResName( Which , i) ,gomp_GetAtomAtmName(Which , i),
                gomp_GetAtomXCoord(Which , i)+Move[0], 
                gomp_GetAtomYCoord(Which , i)+Move[1],
                gomp_GetAtomZCoord(Which , i)+Move[2],
                gomp_GetAtomSegName( Which , i), gomp_GetAtomResNum1( Which , i),gomp_GetAtomBValue( Which , i));
        ihelp++;
    }

    fclose(write_p);

    return(0);
}
/*
  Write coordinates out as a single xyz (XMOL) file 
*/
/************************************************************************/
int gomp_WriteCoordXYZ(int Which , const char *file_name , int alt)
/************************************************************************/
{

    FILE *write_p;
    int i;
    int ihelp;
    char OutText[BUFF_LEN];
    const char *DispList;
    const float *Move;

/* check first structures available */
    i = gomp_GetNumMolecStructs();
    if(!i) {
        gomp_PrintERROR("no structure available");
        return(1);
    }
    if(Which < 1 || Which > i) {
        gomp_PrintERROR("structure index out of range");
        return(1);
    }

    if(file_name[0] == '?') {
        gomp_PrintMessage("?ERROR - file name can't contain '?'");
        return(1);
    }

    write_p = fopen(file_name,"w");
    if(write_p == NULL) {
        sprintf(OutText,"?ERROR - unable to open file : %s ",file_name);
        gomp_PrintMessage(OutText);
        return(1);
    }

    sprintf(OutText,"Writing current molecule/system (nr:%d)",Which);
    gomp_PrintMessage(OutText);

    Which    = Which - 1;
    Move     = gomp_GetTranslateArray();
    DispList = gomp_GetAtomDisplayStatePointer(Which);

    ihelp = 0;
    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {

        if(alt) {
            if(DispList[i] == 0) continue;
        }
        ihelp++;
    }

    fprintf(write_p," %d \n",ihelp);
    fprintf(write_p,"* File written in XYZ (XMOL) format \n");

    ihelp = 1;
    for(i = 0  ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {

        if(alt) {
            if(DispList[i] == 0) continue;
        }
        fprintf(write_p,"%.4s %f %f %f \n",
                gomp_GetAtomAtmName(Which , i),
                (gomp_GetAtomXCoord(Which , i) + Move[0]), 
                (gomp_GetAtomYCoord(Which , i) + Move[1]),
                (gomp_GetAtomZCoord(Which , i) + Move[2]));
        ihelp++;
    }

    fclose(write_p);

    return(0);
}

/**********************************************************************/
int PrepareBaSAtoms(char *element_vec , int Iatoms)
/**********************************************************************/
{


    static int nsymbl= CONV_ATOMS; /* atom symbol list */
/*
 *
 */
    static const char *pt =
        "?? HHeLiBe B C N O FNeNaMgAlSi P SClArCACBCGCDCECZNANBNGNDNENZOAOBOGODOEOZOHSASBSGSDHHHGHZNHCHOCHAHBHCHDHEHNHTHOCa";
    static const char *PT =
        " HHeLiBe B C N O FNeNaMgAlSi P SClArCa";

    static int ihelpv[CONV_ATOMS] = {0,
                                     1,2,
                                     3,4,
                                     5,6,7,8,9,10,
                                     11,12,
                                     13,14,15,16,17,18,
                                     6,6,6,6,6,6,
                                     7,7,7,7,7,7,
                                     8,8,8,8,8,8,8,
                                     16,16,16,16,
                                     1,1,1,7,6,8,1,1,1,1,1,1,1,1,19};

    static int i,j,swtch,nerror;

    static char OutText[BUFF_LEN];

/* clear all numbers in the element vector */
    for(i = 0 ; i < Iatoms ; i++) {

        if(element_vec[2*i+1] == ' ' || 
           element_vec[2*i+1] == '\0') {
            element_vec[2*i+1] = element_vec[2*i];
            element_vec[2*i]   = ' ';
        }

        if(isdigit(element_vec[2*i+1])) {
            element_vec[2*i+1] = element_vec[2*i];
            element_vec[2*i] = ' ';
        }
    }

    for(i = 0 ; i < Iatoms ; i++) {

        swtch = 0;

        for(j = 0 ; j < nsymbl ; j++) { 
            if(strncmp((element_vec+2*i),(pt+2*j),2) == 0)  {
                strncpy((element_vec+2*i),PT+2*(ihelpv[j]-1),2);
                swtch = 1;
                break;
                break;
            }
        }
        if(swtch != 0) continue;

        sprintf(OutText," ?WARNING - can't translate element  >%.2s<",
                element_vec+2*i);
        gomp_PrintMessage(OutText);
        gomp_PrintMessage("Known atoms are:\n");
        for(j = 0 ; j < (signed int)strlen(PT)/ 30 ; j++) {
            sprintf(OutText,"%s",PT+j*30);
            gomp_PrintMessage(OutText);
        }
        if(strlen(PT) - 30 * ( strlen(PT) / 30)) {
            sprintf(OutText,"%s",PT+(strlen(PT) / 30)*30);
            gomp_PrintMessage(OutText);
        }
        nerror++;
    }
    return(0);
}


/************************************************************************/
int gomp_WriteCoordOPENMOL(int Which , const char *file_name , int alt)
/************************************************************************/
{

    FILE *write_p;
    int i;
    int ihelp;
    char OutText[BUFF_LEN];
    const char *DispList;
    float conv = 0.52917715f;
    const float *Move;

/* check first structures available */
    i = gomp_GetNumMolecStructs();
    if(!i) {
        gomp_PrintERROR("no structure available");
        return(1);
    }
    if(Which < 1 || Which > i) {
        gomp_PrintERROR("structure index out of range");
        return(1);
    }

    if(file_name[0] == '?') {
        gomp_PrintMessage("?ERROR - file name can't contain '?'");
        return(1);
    }

    write_p = fopen(file_name,"w");
    if(write_p == NULL) {
        sprintf(OutText,"?ERROR - unable to open file : %s ",file_name);
        gomp_PrintMessage(OutText);
        return(1);
    }

    sprintf(OutText,"Writing current molecule/system (nr:%d)",Which);
    gomp_PrintMessage(OutText);

    Which    = Which - 1;
    DispList = gomp_GetAtomDisplayStatePointer(Which);

    Move  = gomp_GetTranslateArray();
    ihelp = 0;
    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {

        if(alt) {
            if(DispList[i] == 0) continue;
        }
        ihelp++;
    }

    fprintf(write_p,": label\n");
    if(gomp_GetDescText() == '\0') 
        fprintf(write_p,"** default title from gOpenMol **\n");
    else
        fprintf(write_p,"%s\n",gomp_GetDescText());
    fprintf(write_p,": gaussian_basis_set_charge_centers\n");

    ihelp = 1;
    for(i = 0  ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {

        if(alt) {
            if(DispList[i] == 0) continue;
        }
        fprintf(write_p," %s %f %f %f %f %s \n",
                gomp_GetAtomAtmName(Which , i), gomp_GetAtomNucCharge( Which , i),
                (gomp_GetAtomXCoord(Which , i) + Move[0])/conv, 
                (gomp_GetAtomYCoord(Which , i) + Move[1])/conv ,
                (gomp_GetAtomZCoord(Which , i) + Move[2])/conv,
                gomp_GetAtomBasisSetTag( Which , i));
        ihelp++;
    }

    fclose(write_p);

    return(0);
}

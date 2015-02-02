/*

Copyright (c) 1992 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2005 by:
Eero Häkkinen
*/
 
 
/*
  This program reads a MUMOD coordinate file
 
  Leif Laaksonen 1990
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "coord_file.h"
#include "gomstring.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"

#include "stdafx.h"

#define MUMOD_LINE_LEN   120   /* MUMOD file line length */

#define Kill_Comma(InText)  {size_t Cloop; for(Cloop = 0 ; \
                                           Cloop < strlen(InText);\
                                           Cloop++) {\
                                if(InText[Cloop] == ',') InText[Cloop] = ' ';}}

#define FREE1    free(molecu); free(n); free(nsite);
#define FREE2    free(name); free(sitmas); free(charge); free(xpr); free(ypr); free(zpr);
#define FREE3    free(xm); free(ym); free(zm);

static int CalculateMUMODatoms(FILE *);

/*************************************************************************/
int gomp_ReadCoordinatesMUMOD(const char *inp_file , int Append)
/*************************************************************************/
{
 

    float TXc,TYc,TZc;
    int   Mumod_Cards = 1;
    int   Mumod_Atoms = 0;
    int   Atom_point;
    int   Resi_point;
 
    char input_text[MUMOD_LINE_LEN];
    char temp[MUMOD_LINE_LEN];
    int i,j,k,l,poc;
    char OutText[BUFF_LEN];
/* mumod variables */
    int   ntype;
    char **molecu;
    int *n;
    int *nsite;
    char **name;
    float *sitmas;
    float *charge;
    float *xpr;
    float *ypr;
    float *zpr;
    float *xm;
    float *ym;
    float *zm;
 
    FILE *mumod_in;

    int numat;
    int Wstr;
 
    mumod_in=fopen(inp_file,"r");
    if(mumod_in == NULL) {
        sprintf(OutText,"Can't open input file : %s",inp_file);
        gomp_PrintMessage(OutText);
        return(1);
    }

/* ................................................................ */
    rewind(mumod_in);
    if(!(numat = CalculateMUMODatoms(mumod_in)))
        return(1);
    rewind(mumod_in);
/* ................................................................ */
    printf("Number of atoms: %d\n",numat);

    if(Append)
        Wstr = gomp_CreateMolecStruct(inp_file , numat , APPEND); 
    else
        Wstr = gomp_CreateMolecStruct(inp_file , numat , NEW);
    if ( Wstr < 0 )
        goto end;

/* ................................................................ */

    sprintf(OutText,"********** Reading : %s **********",inp_file);
    gomp_PrintMessage(OutText);
/*
  Start reading file
 
*/
    Mumod_Cards = 0;

    if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
        gomp_PrintMessage("?ERROR - can't read first line");
        fclose(mumod_in);
        return(1);
    }

    Mumod_Cards++;
    strcpy(temp,input_text);
    gomp_String2Lower(temp);
    if( (k = gomp_Indexo(temp,"new")) == 0) {
        gomp_PrintMessage("?ERROR - on first line in mumod input file");
        gomp_PrintMessage(input_text);
        fclose(mumod_in);
        return(1);
    }
 
    while(fgets(input_text,MUMOD_LINE_LEN,mumod_in) != NULL) {
        Mumod_Cards++;
        strcpy(temp,input_text);

        gomp_String2Lower(temp);
        if(gomp_Indexo(temp,"end") == 1) break;

        i = strlen(input_text);
        if(input_text[i - 1] == '\n') input_text[i - 1] = '\0';
        gomp_PrintMessage(input_text);
        continue;
    }

    if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
        gomp_PrintMessage("?ERROR - can't read input file");
        fclose(mumod_in);
        return(1);
    }
    gomp_PrintMessage(input_text);
    if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
        gomp_PrintMessage("?ERROR - can't read input file");
        fclose(mumod_in);
        return(1);
    }
    gomp_PrintMessage(input_text);
        
    Mumod_Cards += 2;

    if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
        gomp_PrintMessage("?ERROR - can't read input file");
        fclose(mumod_in);
        return(1);
    }
    ntype = atoi(input_text);

    molecu = malloc(ntype * sizeof(const char *));

    for(i = 0 ; i < ntype ; i++) {  
        molecu[i] = malloc(BUFF_LEN);
        if(molecu[i] == NULL) {
            gomp_PrintMessage("?ERROR - can't alloc memory in 'rmumod.c'");
            fclose(mumod_in);
            return(1);
        }
    }

    n      = gomp_AllocateIntVector(ntype);
    nsite  = gomp_AllocateIntVector(ntype);

    Resi_point = 1;

    Atom_point = 0;

    for(i = 0 ; i < ntype ; i++) {      /* ntype loop */

        if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
            gomp_PrintMessage("?ERROR - can't read input file");
            fclose(mumod_in);
            FREE1;
            return(1);
        }

        sscanf(input_text,"%s %d %d %*s",molecu[i],&n[i],&nsite[i]);

        name = malloc(nsite[i] * sizeof(const char *));

        for(poc = 0 ; poc < nsite[i] ; poc++) {
            name[poc] = malloc(BUFF_LEN);
            if(name[poc] == NULL) {
                gomp_PrintMessage("?ERROR - can't alloc memory in 'rmumod.c'");
                fclose(mumod_in);
                FREE1;
                return(1);
            }
        }

        sitmas = gomp_AllocateFloatVector(nsite[i]);
        charge = gomp_AllocateFloatVector(nsite[i]);
        xpr    = gomp_AllocateFloatVector(nsite[i]);
        ypr    = gomp_AllocateFloatVector(nsite[i]);
        zpr    = gomp_AllocateFloatVector(nsite[i]);

        for(j = 0 ; j < nsite[i] ; j++) {     /* nsite(i) loop */

            if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
                gomp_PrintMessage("?ERROR - can't read input file");
                fclose(mumod_in);
                FREE1;
                FREE2;
                return(1);
            }

            sscanf(input_text,"%s %f %*d %f %f %f %f",name[j],&sitmas[j],&charge[j],
                   &xpr[j],&ypr[j],&zpr[j]);
        }

        for(j = 0 ; j < nsite[i] ; j++) {      /* nsite(i) loop */
            if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
                gomp_PrintMessage("?ERROR - can't read input file");
                fclose(mumod_in);
                FREE1;
                FREE2;
                return(1);
            }
        }


/* build the coordinate data */
        for(k = 0 ; k < n[i] ; k++) {
            for(l = 0 ; l < nsite[i] ; l++) {

                (void)gomp_PutAtomSegName(Wstr , molecu[i],Atom_point);
                (void)gomp_PutAtomResName(Wstr , molecu[i],Atom_point);
                (void)gomp_PutAtomAtmName(Wstr , name[l],Atom_point);
                (void)gomp_PutAtomResNum1(Wstr , Resi_point,Atom_point);
                (void)gomp_PutAtomResNum2(Wstr , Resi_point,Atom_point);
                (void)gomp_PutAtomBValue(Wstr  , 0.0,Atom_point);
                (void)gomp_PutAtomXCoord(Wstr  , xpr[l],Atom_point);
                (void)gomp_PutAtomYCoord(Wstr ,ypr[l],Atom_point);
                (void)gomp_PutAtomZCoord(Wstr,zpr[l],Atom_point);
                (void)gomp_PutAtomCharge(Wstr , charge[l],Atom_point);
                Atom_point++;
            }
            Resi_point++;
        }

        FREE2;
    }

    Atom_point = 0;

    for(i = 0 ; i < ntype ; i++) {

        xm = gomp_AllocateFloatVector(n[i]);
        ym = gomp_AllocateFloatVector(n[i]);
        zm = gomp_AllocateFloatVector(n[i]);

        for(j = 0 ; j < n[i] ; j++) {

            if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
                gomp_PrintMessage("?ERROR - can't read input file (3)");
                fclose(mumod_in);
                FREE1;
                FREE3;
                return(1);
            }

            Kill_Comma(input_text);

            sscanf(input_text,"%f %f %f",&xm[j],&ym[j],&zm[j]);
        }

        for(j = 0 ; j < n[i] ; j++) {
            for(k = 0 ; k < nsite[i] ; k++) {

                TXc = gomp_GetAtomXCoord(Wstr , Atom_point) + xm[j];
                TYc = gomp_GetAtomYCoord(Wstr , Atom_point) + ym[j];
                TZc = gomp_GetAtomZCoord(Wstr , Atom_point) + zm[j];
                (void)gomp_PutAtomXCoord(Wstr , TXc,Atom_point);
                (void)gomp_PutAtomYCoord(Wstr , TYc,Atom_point);
                (void)gomp_PutAtomZCoord(Wstr , TZc,Atom_point);

                Atom_point++;
            }
        }
        FREE3;
    }

    Mumod_Atoms = 0;
    for(i = 0 ; i < ntype ; i++) {
        Mumod_Atoms += n[i] * nsite[i];
    }
    
    if(numat !=  Mumod_Atoms) {
        gomp_PrintMessage("$ERROR - number of atoms does not match");
        exit(1);
    }
 
    gomp_PrintMessage("**********   Done   **********");

    FREE1;
end:
    fclose(mumod_in);
    return(Wstr >= 0 ? 0 : 1);    
}


/*************************************************************************/
int CalculateMUMODatoms(FILE *mumod_in)
/*************************************************************************/
{
 

    int   Mumod_Cards = 1;
    int   Mumod_Atoms = 0;
    int   Atom_point;
    int   Resi_point;
 
    char input_text[MUMOD_LINE_LEN];
    char temp[MUMOD_LINE_LEN];
    int i,j,k,l,poc;
    char OutText[BUFF_LEN];
/* mumod variables */
    int   ntype;
    char **molecu;
    int *n;
    int *nsite;
    char **name;
    float *sitmas;
    float *charge;
    float *xpr;
    float *ypr;
    float *zpr;
    float *xm;
    float *ym;
    float *zm;

/* ................................................................ */
    rewind(mumod_in);
 
    sprintf(OutText,"*** Calculating number of atoms ***");
    gomp_PrintMessage(OutText);
/*
  Start reading file
 
*/
    Mumod_Cards = 0;

    if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
        gomp_PrintMessage("?ERROR - can't read first line");
        fclose(mumod_in);
        return(0);
    }

    Mumod_Cards++;
    strcpy(temp,input_text);
    gomp_String2Lower(temp);
    if( (k = gomp_Indexo(temp,"new")) == 0) {
        gomp_PrintMessage("?ERROR - on first line in mumod input file");
        gomp_PrintMessage(input_text);
        fclose(mumod_in);
        return(0);
    }
 
    while(fgets(input_text,MUMOD_LINE_LEN,mumod_in) != NULL) {
        Mumod_Cards++;
        strcpy(temp,input_text);

        gomp_String2Lower(temp);
        if(gomp_Indexo(temp,"end") == 1) break;

        i = strlen(input_text);
        if(input_text[i - 1] == '\n') input_text[i - 1] = '\0';
        gomp_PrintMessage(input_text);
        continue;
    }

    if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
        gomp_PrintMessage("?ERROR - can't read input file");
        fclose(mumod_in);
        return(0);
    }
    gomp_PrintMessage(input_text);
    if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
        gomp_PrintMessage("?ERROR - can't read input file");
        fclose(mumod_in);
        return(0);
    }
    gomp_PrintMessage(input_text);
        
    Mumod_Cards += 2;

    if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
        gomp_PrintMessage("?ERROR - can't read input file");
        fclose(mumod_in);
        return(0);
    }

    ntype = atoi(input_text);

    molecu = malloc(ntype * sizeof(const char *));

    for(i = 0 ; i < ntype ; i++) {  
        molecu[i] = malloc(BUFF_LEN);
        if(molecu[i] == NULL) {
            gomp_PrintMessage("?ERROR - can't alloc memory in 'rmumod.c'");
            fclose(mumod_in);
            return(0);
        }
    }

    n      = gomp_AllocateIntVector(ntype);
    nsite  = gomp_AllocateIntVector(ntype);

    Resi_point = 1;

    Atom_point = 0;

    for(i = 0 ; i < ntype ; i++) {      /* ntype loop */

        if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
            gomp_PrintMessage("?ERROR - can't read input file");
            fclose(mumod_in);
            FREE1;
            return(0);
        }

        sscanf(input_text,"%s %d %d %*s",molecu[i],&n[i],&nsite[i]);

        name = malloc(nsite[i] * sizeof(const char *));

        for(poc = 0 ; poc < nsite[i] ; poc++) {
            name[poc] = malloc(BUFF_LEN);
            if(name[poc] == NULL) {
                gomp_PrintMessage("?ERROR - can't alloc memory in 'rmumod.c'");
                fclose(mumod_in);
                FREE1;
                return(0);
            }
        }

        sitmas = gomp_AllocateFloatVector(nsite[i]);
        charge = gomp_AllocateFloatVector(nsite[i]);
        xpr    = gomp_AllocateFloatVector(nsite[i]);
        ypr    = gomp_AllocateFloatVector(nsite[i]);
        zpr    = gomp_AllocateFloatVector(nsite[i]);

        for(j = 0 ; j < nsite[i] ; j++) {     /* nsite(i) loop */

            if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
                gomp_PrintMessage("?ERROR - can't read input file");
                fclose(mumod_in);
                FREE1;
                FREE2;
                return(0);
            }

            sscanf(input_text,"%s %f %*d %f %f %f %f",name[j],&sitmas[j],&charge[j],
                   &xpr[j],&ypr[j],&zpr[j]);
        }

        for(j = 0 ; j < nsite[i] ; j++) {      /* nsite(i) loop */
            if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
                gomp_PrintMessage("?ERROR - can't read input file");
                fclose(mumod_in);
                FREE1;
                FREE2;
                return(0);
            }
        }


/* build the coordinate data */
        for(k = 0 ; k < n[i] ; k++) {
            for(l = 0 ; l < nsite[i] ; l++) {

/*
  (void)gomp_PutAtomSegName(Wstr , molecu[i],Atom_point);
  (void)gomp_PutAtomResName(Wstr , molecu[i],Atom_point);
  (void)gomp_PutAtomAtmName(Wstr , name[l],Atom_point);
  (void)gomp_PutAtomResNum1(Wstr , Resi_point,Atom_point);
  (void)gomp_PutAtomResNum2(Wstr , Resi_point,Atom_point);
  (void)gomp_PutAtomBValue(Wstr  , 0.0,Atom_point);
  (void)gomp_PutAtomXCoord(Wstr  , xpr[l],Atom_point);
  (void)gomp_PutAtomYCoord(Wstr ,ypr[l],Atom_point);
  (void)gomp_PutAtomZCoord(Wstr,zpr[l],Atom_point);
  (void)gomp_PutAtomCharge(Wstr,charge[l],Atom_point);
*/
                Atom_point++;
            }
            Resi_point++;
        }

        FREE2;
    }

    for(i = 0 ; i < ntype ; i++) {

        xm = gomp_AllocateFloatVector(n[i]);
        ym = gomp_AllocateFloatVector(n[i]);
        zm = gomp_AllocateFloatVector(n[i]);

        for(j = 0 ; j < n[i] ; j++) {

            if(fgets(input_text,MUMOD_LINE_LEN,mumod_in) == NULL) {
                gomp_PrintMessage("?ERROR - can't read input file (3)");
                fclose(mumod_in);
                FREE1;
                FREE3;
                return(0);
            }

            Kill_Comma(input_text);

            sscanf(input_text,"%f %f %f",&xm[j],&ym[j],&zm[j]);
        }

        for(j = 0 ; j < n[i] ; j++) {
            for(k = 0 ; k < nsite[i] ; k++) {

                Atom_point++;
            }
        }
        FREE3;
    }

/* it's possible to already update now */

    Mumod_Atoms = 0;
    for(i = 0 ; i < ntype ; i++) {
        Mumod_Atoms += n[i] * nsite[i];
    }
    
    gomp_PrintMessage("**********   Done   **********");

    FREE1;
    return(Mumod_Atoms);    
}



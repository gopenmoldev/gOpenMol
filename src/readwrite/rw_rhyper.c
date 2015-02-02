/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  
Enhancements 2003, 2005 by:
Eero Häkkinen

*/
 
 
/*
  This program reads a HYPERCHEM coordinate file
 
  Leif Laaksonen 1993, 1994
*/
 
#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <string.h>

#include "coord_file.h"
#include "gomstring.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"

#include "stdafx.h"

#define HYPER_LINE_LEN   132   /* hyper file line length */
 
static int AtomsInHYPERCHEMFile(FILE *);

/*************************************************************************/
int gomp_ReadCoordinatesHYPERCHEM(const char *Text1 , int Append)
/*************************************************************************/
{
 

    int   TRs1,TRs2,NumAtomMol;
    float TXc,TYc,TZc,TBv,TCh;
    char  TResN[BUFF_LEN];
    char  TAtmN[BUFF_LEN];
    char  TSegN[BUFF_LEN];

 
    char inputl[HYPER_LINE_LEN];
    int  j,loop;
    char OutText[BUFF_LEN];
    int  MolNum = 0;
    int  HyperChemAtoms;
    int  Wstr;
 
    FILE *hyper_in;
 
    hyper_in=fopen(Text1,"r");
    if(hyper_in == NULL) {
        sprintf(OutText,"$Can't open input file : %s",Text1);
        gomp_PrintERROR(OutText);
        return(1);
    }
 
    TRs1 = 1;
    TRs2 = 1;
    TBv  = 0.0;

    sprintf(OutText,"********** Reading : %s **********",Text1);
    gomp_PrintMessage(OutText);
/*
  Start reading file
 
*/
/* it is not possible say how many atom will be read at this stage so
   reserve space for one (1) atom first and add any further space needed
   This is ugly and slow but a "quick and dirty hack */

    HyperChemAtoms = AtomsInHYPERCHEMFile(hyper_in);

    if(Append)
        Wstr        = gomp_CreateMolecStruct(Text1 , HyperChemAtoms , APPEND);
    else
        Wstr        = gomp_CreateMolecStruct(Text1 , HyperChemAtoms , NEW);
    if ( Wstr < 0 )
        goto end;

    loop = 0;

    while(fgets(inputl,HYPER_LINE_LEN,hyper_in) != NULL) { 

/* Title/comment card (;) */
        if(inputl[0] == ';') 
            gomp_PrintMessage(&inputl[1]);

/* Coordinates card (mol) */
        if(gomp_Indexo(inputl,"mol") == 1) {

            sscanf(inputl,"%*s %d %s",&TRs1,TSegN);

            MolNum = TRs1;

            if(TSegN[0] != '\0')
                gomp_CopyString(TResN,TSegN,BUFF_LEN);
            else   {
                strncpy(TSegN,"HYPE",4);
                strncpy(TResN,"HYPE",4);
            }
        }

/* Coordinates card (endmol) */
        else if(gomp_Indexo(inputl,"endmol") == 1) {

            sscanf(inputl,"%*s %d",&TRs2);

            if(TRs2 < 1) TRs2 = 1;

            if(MolNum != TRs2) {
                gomp_PrintWARNING("?mol/endmol numbers do not match");
            }
        }

/* residue card           */
        else if(gomp_Indexo(inputl,"endres") == 1) {

            sscanf(inputl,"%d",&TRs1);

            if(TRs1 != TRs2) {
                gomp_PrintWARNING("?residue number in res/endres does not mach");
            }
        }

        else if(gomp_Indexo(inputl,"res") == 1) {

            sscanf(inputl,"%*s %d %s",&TRs1,OutText);

            TRs2 = TRs1;

            if(OutText[0] == '\0')
                strncpy(TResN,"HYPE",4);
            else
                gomp_CopyString(TResN,OutText,BUFF_LEN);
        }


/* atom card              */
        else if(gomp_Indexo(inputl,"atom") == 1) {

            sscanf(inputl,"%*s %d %s %s %*s %*c %f %f %f %f", &NumAtomMol,
                   TAtmN, OutText,
                   &TCh,&TXc,&TYc,&TZc);

            TRs2 = TRs1;

            j = gomp_PutAtomResNum1(Wstr , TRs1 , loop);
            j = gomp_PutAtomResNum2(Wstr , TRs2 , loop);

            j = gomp_PutAtomResName(Wstr , TResN , loop);

            if(TAtmN[0] == '-') {
                j = gomp_PutAtomAtmName(Wstr , OutText, loop);
            }
            else {
                j = gomp_PutAtomAtmName(Wstr , TAtmN, loop);
            }

            j = gomp_PutAtomXCoord(Wstr , TXc , loop);
            j = gomp_PutAtomYCoord(Wstr , TYc , loop);
            j = gomp_PutAtomZCoord(Wstr , TZc , loop);


            j = gomp_PutAtomSegName(Wstr , TSegN , loop);

            j = gomp_PutAtomBValue(Wstr , TBv , loop);
            j = gomp_PutAtomCharge(Wstr , TCh , loop);

            loop++;
        }
    }

    if(loop != HyperChemAtoms) {
        gomp_PrintWARNING("Number of atoms read from HyperChem file does not match");
        return(1);
    }

    gomp_PrintMessage("**********   Done   **********");
    
end:
    fclose(hyper_in);

    return(Wstr >= 0 ? 0 :1);
}
/*************************************************************************/
int AtomsInHYPERCHEMFile(FILE *hyper_in)
/*************************************************************************/
{ 
    char inputl[HYPER_LINE_LEN];
    int  HyperChemAtoms;
 
/* rewind file (just for sure) */
    rewind(hyper_in);

/*
  Start reading file
 
*/

    HyperChemAtoms = 0 ; /* start counting number of atoms */

    while(fgets(inputl,HYPER_LINE_LEN,hyper_in) != NULL) { 

/* atom card              */
        if(gomp_Indexo(inputl,"atom") == 1) {
            HyperChemAtoms++;
        }
    }

/* rewind file for the real processing */
    rewind(hyper_in);

    return(HyperChemAtoms);
}

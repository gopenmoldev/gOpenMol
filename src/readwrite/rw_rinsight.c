/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  
Enhancements 2005 by:
Eero Häkkinen

*/


/*
  This program reads a Insight *.car coordinate file

  Leif Laaksonen 1990, 1994 (extended for gOpenMol)
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "coord_file.h"
#include "gomstring.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"

#include "stdafx.h"

#define INST_LINE_LEN   120 /* insight file line length is 80 char
                               but this prevents from getting less than 80 const char */

static int AtomsInINSIGHTFile(FILE *);

/*************************************************************************/
int gomp_ReadCoordinatesINSIGHT(const char *Text1 , int Append)
/*************************************************************************/
{
    static int i,j,numat,Wstr;
    static int InsightAtoms;

    char    input[INST_LINE_LEN];
    char    temp[BUFF_LEN];
    char    OutText[BUFF_LEN];

    FILE *chm_in;

    chm_in=fopen(Text1,"r");
    if(chm_in == NULL) {
        sprintf(OutText,"Can't open input file : %s",Text1);
        gomp_PrintMessage(OutText);
        return(1);
    }

    sprintf(OutText,"********** Reading : %s **********",Text1);
    gomp_PrintMessage(OutText);
    gomp_PrintMessage("           Title   : ");
/*
  Start reading file

*/

    InsightAtoms = AtomsInINSIGHTFile(chm_in);

/*  1     */

    fgets(input,INST_LINE_LEN,chm_in);
    gomp_PrintMessage(input);

/*  2     */

    fgets(input,INST_LINE_LEN,chm_in);
    gomp_PrintMessage(input);

    i = 0;
    if(gomp_Indexo(input,"on") > 0 ||
       gomp_Indexo(input,"ON") > 0) i = 1; /* PBC is on */
/*  3     */

    fgets(input,INST_LINE_LEN,chm_in);
    gomp_PrintMessage(input);

/*  4     */

    fgets(input,INST_LINE_LEN,chm_in);
    gomp_PrintMessage(input);

/*  5     */

    if(i) {
        fgets(input,INST_LINE_LEN,chm_in);
        gomp_PrintMessage(input);
    }

/*  6 Read atom cards  */

    if(Append)
        Wstr = gomp_CreateMolecStruct(Text1 , InsightAtoms , APPEND);
    else
        Wstr = gomp_CreateMolecStruct(Text1 , InsightAtoms , NEW);
    if ( Wstr < 0 )
        goto end;

    numat = 0 ; /* start counting number of atoms */

    while(fgets(input,INST_LINE_LEN,chm_in) != NULL) { 

/* check for empty lines containing just space */
        if(sscanf(input,"%s",temp) < 1) continue;
        strncpy(temp,input,INST_LINE_LEN);
        gomp_String2Lower(temp);
        if(gomp_Indexo(temp,"end") == 1) continue;

        j = gomp_PutAtomAtmName(Wstr , input , numat);   /* atom name */
        strncpy(temp,input+5,15);    /* x - coordinate */
        j = gomp_PutAtomXCoord(Wstr , atof(temp) , numat);
        strncpy(temp,input+20,15);    /* y - coordinate */
        j = gomp_PutAtomYCoord(Wstr , atof(temp) , numat); 
        strncpy(temp,input+35,15);    /* z - coordinate */
        j = gomp_PutAtomZCoord(Wstr , atof(temp) , numat);
        j = gomp_PutAtomResName(Wstr , input+51 ,  numat); /* residue name */
        strncpy(temp,input+55,5);       /* residue sequence number relative to
                                           the beginning of the current mol */
        j = gomp_PutAtomResNum1(Wstr , atoi(temp) , numat);

        strncpy(temp,input+65,10);   /* partial charge */
/*       atm_charge[i] = atof(temp); */
        strncpy(temp,input+75,5);  /* absolute atomic sequence number relative
                                      to the beginning of the entire molecule */ 
        temp[5] = '\0'; /* prevent from including the <cr> */
        j = gomp_PutAtomResNum2(Wstr , atoi(temp) , numat);

/* value added stuff      */
        j = gomp_PutAtomSegName(Wstr , "INTH" , numat); /* give it a segment name */
/* ...................... */


/*         j = gomp_PutAtomCharge( 0.0 , i); */
        numat++;

    }

    if(numat != InsightAtoms) {
        gomp_PrintMessage("?ERROR - can't read correct number of atoms from Insight file");
        return(1);
    }

    gomp_PrintMessage("**********   Done   **********");
    
end:
    fclose(chm_in);

    return(Wstr >= 0 ? 0 : 1);
}
/*************************************************************************/
int AtomsInINSIGHTFile(FILE *chm_in)
/*************************************************************************/
{
    static int InsightAtoms;

    int  i;
    char input[INST_LINE_LEN];
    char temp[INST_LINE_LEN];

/*
  Start reading file

*/

/* rewind just for sure */
    rewind(chm_in);

/*  1     */

    fgets(input,INST_LINE_LEN,chm_in);

/*  2     */

    fgets(input,INST_LINE_LEN,chm_in);

    i = 0;
    if(gomp_Indexo(input,"on") > 0 ||
       gomp_Indexo(input,"ON") > 0) i = 1; /* PBC is on */
/*  3     */

    fgets(input,INST_LINE_LEN,chm_in);

/*  4     */

    fgets(input,INST_LINE_LEN,chm_in);

/*  5     */

    if(i) {
        fgets(input,INST_LINE_LEN,chm_in);
    }

/*  6 Read atom cards  */

    InsightAtoms = 0;  /* start counting number of atoms */

    while(fgets(input,INST_LINE_LEN,chm_in) != NULL) { 

/* check for empty lines containing just space */
        if(sscanf(input,"%s",temp) < 1) continue;

        strncpy(temp,input,INST_LINE_LEN);
        gomp_String2Lower(temp);
        if(gomp_Indexo(temp,"end") == 1) continue;

        InsightAtoms++;

    }

/* rewind file for the real processing */
    rewind(chm_in);

    return(InsightAtoms);
}

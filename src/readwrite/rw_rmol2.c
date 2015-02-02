/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Enhancements 2003, 2005 by:
Eero Häkkinen
*/


/*
  This program reads a SYBYL mol2 coordinate file

  Leif Laaksonen 1995

  This code is based on the "rdsybyl.c" program from the Babel Program.
  The original code is by Jussi Eloranta.

  Supported record types are:
  (1) molecule
  (2) atom

  others are ignored.

*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "coord_file.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"

#include "stdafx.h"

#define MOL2_LINE_LEN   132   /* pdb file line length */

#define MOL  "@<TRIPOS>MOLECULE"
#define ATO  "@<TRIPOS>ATOM"

/*************************************************************************/
int gomp_ReadCoordinatesMOL2(const char *inp_file , int Append)
/*************************************************************************/
{


    static char inputl[MOL2_LINE_LEN];
    static int tatomn;
    static int i,j;
    static int type_warning;
    static int Wstr;
    static int numat;
    static int tatoms;
    static int sets;

    int     TRs1;
    float   TXc,TYc,TZc,TBv;
    char    TResN[BUFF_LEN];
    char    TAtmN[BUFF_LEN];

    char OutText[BUFF_LEN];

    FILE *chm_in;

    type_warning = 0;

    chm_in=fopen(inp_file,"r");
    if(chm_in == NULL) {
        sprintf(OutText,"Can't open input file : %s",inp_file);
        gomp_PrintMessage(OutText);
        return(1);
    }

    sprintf(OutText,"********** Reading : %s **********",inp_file);
    gomp_PrintMessage(OutText);
    gomp_PrintMessage("           Title   : ");
/*
  Start reading file

*/
/* figure out total number of atoms */
    tatoms = 0;
    sets   = 0;
    while(fgets(inputl,MOL2_LINE_LEN,chm_in)) {

        if(!*inputl || *inputl == '#') continue;

/* @<TRIPOS>MOLECULE record (6 data lines) */
        if(!strncmp(inputl,MOL,sizeof(MOL) - 1)) {
            /* molecule name */
            fgets(inputl,MOL2_LINE_LEN,chm_in);
            /* # atoms, # bonds, # substructures, # features , # sets */
            fgets(inputl,MOL2_LINE_LEN,chm_in);
            sscanf(inputl,"%d %*d %*d %*d %*d",&numat);
            sets++;
            sprintf(OutText,"Atoms: %d in set: %d",numat,sets);
            gomp_PrintMessage(OutText);
            tatoms = tatoms + numat;
        }
    }

    sprintf(OutText,"Total number of sets: %d",sets);
    gomp_PrintMessage(OutText);
    sprintf(OutText,"Total number of atoms: %d",tatoms);
    gomp_PrintMessage(OutText);

    rewind(chm_in);

    if(sets > 1) gomp_PrintWARNING("can only read first molecule in a multi molecule MOL2 file!");

    sets = 0;

    while(fgets(inputl,MOL2_LINE_LEN,chm_in) != NULL) {

        if(*inputl =='#') gomp_PrintMessage(inputl);

        if(!*inputl || *inputl == '#') continue;

/* @<TRIPOS>MOLECULE record (6 data lines) */
        if(!strncmp(inputl,MOL,sizeof(MOL) - 1)) {

            if(sets) break;

            /* molecule name */
            fgets(inputl,MOL2_LINE_LEN,chm_in);
            /* # atoms, # bonds, # substructures, # features , # sets */
            fgets(inputl,MOL2_LINE_LEN,chm_in);
            sscanf(inputl,"%d %*d %*d %*d %*d",&numat);

            if(Append)
                Wstr = gomp_CreateMolecStruct(inp_file , numat , APPEND);
            else
                Wstr = gomp_CreateMolecStruct(inp_file , numat , NEW);
            if ( Wstr < 0 )
                goto end;

            sets++;

/*
 * molecule type *
 fgets(inputl,MOL2_LINE_LEN,chm_in);
 * types of charges * 
 fgets(inputl,MOL2_LINE_LEN,chm_in);
 * internal status bit *
 fgets(inputl,MOL2_LINE_LEN,chm_in);
 * comment *
 fgets(inputl,MOL2_LINE_LEN,chm_in);
*/
            while(fgets(inputl,MOL2_LINE_LEN,chm_in) != NULL) {
                if(!strncmp(inputl,ATO,sizeof(ATO) - 1)) {

                    for(i = 0 ; i < numat; i++ ) {

                        fgets(inputl,MOL2_LINE_LEN,chm_in);

                        sscanf(inputl,"%d %s %f %f %f %*s %d %s %f",
                               &tatomn,TAtmN,&TXc,&TYc,&TZc,&TRs1,TResN,&TBv);
                        j = gomp_PutAtomSegName(Wstr , "SYB" , i);
                        j = gomp_PutAtomResName(Wstr , TResN , i);
                        j = gomp_PutAtomAtmName(Wstr , TAtmN , i);
                        j = gomp_PutAtomXCoord(Wstr , TXc , i);
                        j = gomp_PutAtomYCoord(Wstr , TYc , i);
                        j = gomp_PutAtomZCoord(Wstr , TZc , i);

                        j = gomp_PutAtomResNum1(Wstr , TRs1 , i);
                        j = gomp_PutAtomResNum2(Wstr , TRs1 , i);

                        j = gomp_PutAtomBValue(Wstr , 0.0 , i);
                        j = gomp_PutAtomCharge(Wstr , TBv , i);

                    }
                    break;
                }
            }
        }
    }
    gomp_PrintMessage("**********   Done   **********");

end:
    fclose(chm_in);

    return(Wstr >= 0 ? 0 : 1);
}

/*  

Copyright (c) 1995 - 2005 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Enhancements 2003, 2005 by:
Eero Häkkinen
*/


/*
  This program reads a formatted Gaussian checkpoint file

  Leif Laaksonen 1995
  
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "coord_file.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"

#include "stdafx.h"

#define GAUSSIAN_LINE_LEN   120   /* GAUSSIAN file line length */
#define NUMBER_OF_ATOMS               "Number of atoms"
#define ATOMIC_NUMBERS                "Atomic numbers"
#define CURRENT_CARTESIAN_COORDINATES "Current cartesian coordinates"

/*************************************************************************/
int gomp_ReadCoordinatesGAUSSIAN(const char *Text1, int Append)  
    /* GAUSSIAN format file reader */
/*************************************************************************/
{
    static int i,j;
    char inputl[GAUSSIAN_LINE_LEN];   /* input line */
    char OutText[BUFF_LEN];
    int  GaussianAtoms;

    int     Wstr;
    int     AtNum;
    float   TXc,TYc,TZc;
    char    TAtmN[BUFF_LEN];

/* AUDIT: Are these two really allocated before use? */
    int    *AtmNumber_p = NULL;
    float  *AtmCoord_p  = NULL;
    FILE   *gaussian_in;

    gaussian_in=fopen(Text1,"r");
    if(gaussian_in == NULL) {
        sprintf(OutText,"Can't open input file : %s",Text1);
        gomp_PrintMessage(OutText);
        return(1);
    }

    printf("********** Reading : %s **********\n",Text1);

/* loop through the whole file to get number of atoms ... */
    GaussianAtoms = 0;
    while(fgets(inputl,GAUSSIAN_LINE_LEN,gaussian_in) != NULL) {
        if(strncmp(inputl,NUMBER_OF_ATOMS,strlen(NUMBER_OF_ATOMS)) == 0) {
            sscanf(&inputl[strlen(NUMBER_OF_ATOMS)],"%*s %d",&GaussianAtoms);
        }
    }

    if(!GaussianAtoms) {
        gomp_PrintERROR("Can't find number of atoms in Gaussian input file");
        return(1);
    }

/* it's possible to already update now */

    if(Append)
        Wstr = gomp_CreateMolecStruct(Text1 , GaussianAtoms , APPEND);
    else
        Wstr = gomp_CreateMolecStruct(Text1 , GaussianAtoms , NEW);
    if ( Wstr < 0 )
        goto end;

    rewind(gaussian_in);

    while(fgets(inputl,GAUSSIAN_LINE_LEN,gaussian_in) != NULL) {

/* get the atomic numbers to get the atomic symbols */
        if(strncmp(inputl,ATOMIC_NUMBERS,strlen(ATOMIC_NUMBERS)) == 0) {

            AtmNumber_p = gomp_AllocateIntVector(GaussianAtoms);
            for(i = 0 ; i < GaussianAtoms ; i++) {
                fscanf(gaussian_in,"%d",&AtmNumber_p[i]);
/*printf("Atom number: %d\n",AtmNumber_p[i]); */
            }
        }

/* get the cartesian coordinates */
        if(strncmp(inputl,CURRENT_CARTESIAN_COORDINATES,
                   strlen(CURRENT_CARTESIAN_COORDINATES))  == 0) {
            AtmCoord_p = gomp_AllocateFloatVector(3 * GaussianAtoms);

            for(i = 0 ; i < 3 * GaussianAtoms ; i++) {
                fscanf(gaussian_in,"%f",&AtmCoord_p[i]);
/*printf("'%d' %f\n",i,AtmCoord_p[i]); */
            }

            AtNum = 0;
            for(i = 0 ; i < 3 * GaussianAtoms ; i += 3) {
                TXc = BOHR_RADIUS * AtmCoord_p[i];
                TYc = BOHR_RADIUS * AtmCoord_p[i + 1];
                TZc = BOHR_RADIUS * AtmCoord_p[i + 2];

                j = gomp_PutAtomXCoord(Wstr , TXc , AtNum);
                j = gomp_PutAtomYCoord(Wstr , TYc , AtNum);
                j = gomp_PutAtomZCoord(Wstr , TZc , AtNum);

                j = gomp_PutAtomBValue(Wstr , 0.0 , AtNum);
                j = gomp_PutAtomCharge(Wstr ,  0.0 , AtNum);

                j = gomp_PutAtomSegName(Wstr , "GAUS" , AtNum);
                j = gomp_PutAtomResName(Wstr , "GAUS" , AtNum);

                strncpy(TAtmN,gomp_Number2Name(AtmNumber_p[AtNum]),BUFF_LEN-1);
/*printf("'%s'%s'\n",TAtmN,gomp_Number2Name(AtmNumber_p[AtNum]));*/
                j = gomp_PutAtomAtmName(Wstr , TAtmN , AtNum);

                j = gomp_PutAtomResNum1(Wstr , 1 , AtNum);
                j = gomp_PutAtomResNum2(Wstr , 1, AtNum);
                AtNum++;
            }
        }
    }

end:
    free(AtmCoord_p);
    free(AtmNumber_p);
    fclose(gaussian_in);

    return(Wstr >= 0 ? 0 : 1);

}

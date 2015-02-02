/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  
Enhancements 2003, 2005 by:
Eero Häkkinen

*/


/*
  This program reads a AMBER coordinate files.
  Files are very close to pdb files

  Leif Laaksonen 1989,1994
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

#define AMBER_LINE_LEN   120   /* amber file line length */

/*************************************************************************/
int gomp_ReadCoordinatesAMBER(const char *Text1 , int Append)  
    /* AMBER format file reader */
/*************************************************************************/
{

    char inputl[AMBER_LINE_LEN];   /* input line */

    char OutText[BUFF_LEN];

    float TXc,TYc,TZc,TBv;
    int   TRs1;

    static int i,j,loop;
    static int type_warning;
    static int AmberAtoms;
    static int Wstr;


    FILE *amber_in;

    type_warning = 0;

    amber_in=fopen(Text1,"r");
    if(amber_in == NULL) {
        sprintf(OutText,"can't open input file : %s",Text1);
        gomp_PrintERROR(OutText);
        return(1);
    }

    sprintf(OutText,"********** Reading : %s **********",Text1);
    gomp_PrintMessage(OutText);

/*
  Start reading file
  Tags recogniced by this routine are:

  ATOM:   Atom coordinate records for "standard groups"
  HETATM: Atom coordinate records for "non-standard" groups
  END:    End-of-entry record
*/

/* text label */
    fgets(inputl,AMBER_LINE_LEN,amber_in);
    gomp_PrintMessage(inputl);

/* read number of atoms */
    fgets(inputl,AMBER_LINE_LEN,amber_in);
    sscanf(inputl,"%d",&AmberAtoms);

    if(Append)
        Wstr = gomp_CreateMolecStruct(Text1 , AmberAtoms , APPEND);
    else
        Wstr = gomp_CreateMolecStruct(Text1 , AmberAtoms , NEW);
    if ( Wstr < 0 )
        goto end;
         
    loop = 0;
    i    = 0;

    while(fgets(inputl,AMBER_LINE_LEN,amber_in) != NULL) { 

/* ATOM and HETATM */
        if((strncmp(inputl,"ATOM",4)   == 0) ||
           (strncmp(inputl,"HETATM",6) == 0) ) { /*start atom or hetatm*/

            if(AmberAtoms == loop) {
                sprintf(OutText,"number of atoms defined do not match real number\nWill only import first '%d' atoms",loop);
                gomp_PrintERROR(OutText);
                fclose(amber_in);
                return(1);
            }

            j = gomp_PutAtomAtmName(Wstr , (inputl+12) , i);
            j = gomp_PutAtomResName(Wstr , (inputl+17) , i);
            j = gomp_PutAtomSegName(Wstr , (inputl+21) , i);

            sscanf(inputl+22,"%d",&TRs1);
            j = gomp_PutAtomResNum1(Wstr , TRs1 , i);
            j = gomp_PutAtomResNum2(Wstr , TRs1 , i);
            sscanf(inputl+30,"%f %f %f",&TXc,&TYc,&TZc);
            j = gomp_PutAtomXCoord(Wstr , TXc , i);
            j = gomp_PutAtomYCoord(Wstr , TYc , i);
            j = gomp_PutAtomZCoord(Wstr , TZc , i);
            sscanf(inputl+60,"%f",&TBv);
            j = gomp_PutAtomBValue(Wstr , TBv , i);
            j = gomp_PutAtomCharge(Wstr ,  0.0 , i);

            if(inputl[21] == ' ')
                gomp_PutAtomSegName(Wstr , "S1",i);

            i++;
            loop++;

        } /* end atom or hetatm */

/* END */

        if(strncmp(inputl,"END",3) == 0) break; /*the end*/

    }

    gomp_PrintMessage("**********   Done   **********");

end:
    fclose(amber_in);

    return(Wstr >= 0 ? 0 : 1);

}

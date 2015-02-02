/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Enhancements 2003, 2005 by:
Eero Häkkinen
*/


/*
  This program reads a plain xyz coordinate file

  Leif Laaksonen 1994
  
  The format is:

  (1) number_of_atoms
  (2 ..n) Atom symbom x y z coordinates

  Example:
  3
  N  26.465  27.452  -2.490  
  CA 25.497  26.862  -1.573 
  C  26.193  26.179   -.437 

*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "coord_file.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"

#include "stdafx.h"

#define FREE_LINE_LEN   120   /* FREE file line length */

/*************************************************************************/
int gomp_ReadCoordinatesGXYZ(const char *Text1, int Append)  
    /* GXYZ format file reader */
/*************************************************************************/
{
    static int i,j;
    char inputl[FREE_LINE_LEN];   /* input line */
    char OutText[BUFF_LEN];
    int FreeAtoms;

    int     Wstr;
    float   TXc,TYc,TZc;
    char    TAtmN[BUFF_LEN];

    FILE *free_in;

    free_in=fopen(Text1,"r");
    if(free_in == NULL) {
        sprintf(OutText,"Can't open input file : %s",Text1);
        gomp_PrintMessage(OutText);
        return(1);
    }

    printf("********** Reading : %s **********\n",Text1);

    fgets(inputl,FREE_LINE_LEN,free_in);
    sscanf(inputl,"%s",OutText);
    if(isdigit(OutText[0])) {
        sscanf(inputl,"%d",&FreeAtoms);
    } else {
        FreeAtoms = 1;
        while(fgets(inputl,FREE_LINE_LEN,free_in) != NULL) {
            if(!sscanf(inputl,"%s",OutText)) continue;
            FreeAtoms++;
        }
        rewind(free_in);
    }

/* it's possible to already update now */

    if(Append)
        Wstr = gomp_CreateMolecStruct(Text1 , FreeAtoms , APPEND);
    else
        Wstr = gomp_CreateMolecStruct(Text1 , FreeAtoms , NEW);
    if ( Wstr < 0 )
        goto end;

    for(i = 0 ; i < FreeAtoms; i++ ) {
        fgets(inputl,FREE_LINE_LEN,free_in);

        sscanf(inputl,"%s %f %f %f",TAtmN,&TXc,&TYc,&TZc);

        j = gomp_PutAtomSegName(Wstr , "XYZ" , i);
        j = gomp_PutAtomResName(Wstr , "XYZ" , i);
        j = gomp_PutAtomAtmName(Wstr , TAtmN , i);

        j = gomp_PutAtomXCoord(Wstr , TXc , i);
        j = gomp_PutAtomYCoord(Wstr , TYc , i);
        j = gomp_PutAtomZCoord(Wstr , TZc , i);

        j = gomp_PutAtomResNum1(Wstr , 1 , i);
        j = gomp_PutAtomResNum2(Wstr , 1, i);

        j = gomp_PutAtomBValue(Wstr , 0.0 , i);
        j = gomp_PutAtomCharge(Wstr ,  0.0 , i);


    }

end:
    fclose(free_in);

    return(Wstr >= 0 ? 0 : 1);

}

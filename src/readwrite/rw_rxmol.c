
/*

Copyright (c) 1996 - 2005 by:
Leif Laaksonen , Center for Scientific Computing , ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved

Enhancements 2003, 2005 by:
Eero Häkkinen
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

#define XMOL_LINE_LEN   120   /*XMOL xyz file line length */

/*********************************************************************/
int gomp_ReadCoordinatesXMOL(const char *Text1 , int Append)
/*********************************************************************/
{
    static int i,j;
    static int XmolAtoms;

    int     Wstr;
    float   TXc,TYc,TZc,TBv;
    char    TAtmN[BUFF_LEN];

    char inputl[XMOL_LINE_LEN];
    char OutText[BUFF_LEN];
    int  items;

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

/*  Two title lines are expected    */

/*   Numbers of atoms */
    fgets(inputl,XMOL_LINE_LEN,chm_in);
    sscanf(inputl,"%d",&XmolAtoms);
    fgets(inputl,XMOL_LINE_LEN,chm_in);
    gomp_PrintMessage(inputl);


/* it's possible to already update now */

    if(Append)
        Wstr = gomp_CreateMolecStruct(Text1 , XmolAtoms , APPEND); 
    else
        Wstr = gomp_CreateMolecStruct(Text1 , XmolAtoms , NEW);
    if ( Wstr < 0 )
        goto end;

/*   Read atom cards   */

    for(i = 0 ; i < XmolAtoms ; i++ ) {
        fgets(inputl,XMOL_LINE_LEN,chm_in);

        items = sscanf(inputl," %s %f %f %f %f ",
                       TAtmN,&TXc,&TYc,&TZc,&TBv);

        j = gomp_PutAtomSegName(Wstr , "xmol" , i);
        j = gomp_PutAtomResName(Wstr , "xmol" , i);
        j = gomp_PutAtomAtmName(Wstr , TAtmN , i);

        j = gomp_PutAtomXCoord(Wstr , TXc , i);
        j = gomp_PutAtomYCoord(Wstr , TYc , i);
        j = gomp_PutAtomZCoord(Wstr , TZc , i);

        j = gomp_PutAtomResNum1(Wstr , 1 , i);
        j = gomp_PutAtomResNum2(Wstr , 1 , i);

        j = gomp_PutAtomBValue(Wstr , 0.0 , i);
        j = gomp_PutAtomCharge(Wstr ,  TBv , i);

    }

    gomp_PrintMessage("**********   Done   **********");

end:
    fclose(chm_in);

/*     ShowDataStructure(0 , XmolAtoms);  */

    return(Wstr >= 0 ? 0 : 1);
}


/*

Copyright (c) 1991 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing , ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved

Enhancements 2002 - 2005 by:
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

#define KARP_LINE_LEN   120   /* karp/charmm file line length */

/*********************************************************************/
int gomp_ReadCoordinatesCHARMM(const char *Text1 , int Append)
/*********************************************************************/
{
    static int i,j;
    static int CharmmAtoms;

    int     TRs1,TRs2;
    int     Wstr;
/*  int     AtNum; */
    float   TXc,TYc,TZc,TBv;
    char    TSegN[BUFF_LEN];
    char    TResN[BUFF_LEN];
    char    TAtmN[BUFF_LEN];
    char    Temp[BUFF_LEN];

    char inputl[KARP_LINE_LEN];
    char OutText[BUFF_LEN];

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

/*  1.0 A title is expected    */

    fgets(inputl,KARP_LINE_LEN,chm_in);
    gomp_PrintMessage(inputl);

    while(strncmp(inputl,"*",1) == 0) {
        fgets(inputl,KARP_LINE_LEN,chm_in);
        if(strncmp(inputl,"*",1) != 0) break;
        gomp_PrintMessage(inputl);
    }

/*   2.0 Numbers of atoms */

    sscanf(inputl,"%d",&CharmmAtoms);

/* it's possible to already update now */

    if(Append)
        Wstr = gomp_CreateMolecStruct(Text1 , CharmmAtoms , APPEND); 
    else
        Wstr = gomp_CreateMolecStruct(Text1 , CharmmAtoms , NEW);
    if ( Wstr < 0 )
        goto end;

/*   3.0 Read atom cards   */

    for(i = 0 ; i < CharmmAtoms ; i++ ) {
        fgets(inputl,KARP_LINE_LEN,chm_in);

        memset(Temp , 0 , (size_t)BUFF_LEN);
        strncpy(Temp , (inputl+5) , 5);
        sscanf(Temp,"%d",&TRs1);

        memset(Temp , 0 , (size_t)BUFF_LEN);
        strncpy(Temp , (inputl+11) , MAX_RES_NAME_LEN);
        sscanf(Temp,"%s",TResN);

        memset(Temp , 0 , (size_t)BUFF_LEN);
        strncpy(Temp , (inputl+16) , MAX_ATM_NAME_LEN);
        sscanf(Temp,"%s",TAtmN);

        memset(Temp , 0 , (size_t)BUFF_LEN);
        strncpy(Temp , (inputl+20) , 10);
        sscanf(Temp,"%f",&TXc);

        memset(Temp , 0 , (size_t)BUFF_LEN);
        strncpy(Temp , (inputl+30) , 10);
        sscanf(Temp,"%f",&TYc);

        memset(Temp , 0 , (size_t)BUFF_LEN);
        strncpy(Temp , (inputl+40) , 10);
        sscanf(Temp,"%f",&TZc);

        memset(Temp , 0 , (size_t)BUFF_LEN);
        strncpy(Temp , (inputl+51) , MAX_SEG_NAME_LEN);
        sscanf(Temp,"%s",TSegN);

        memset(Temp , 0 , (size_t)BUFF_LEN);
        strncpy(Temp , (inputl+56) , 4);
        sscanf(Temp,"%d",&TRs2);

        memset(Temp , 0 , (size_t)BUFF_LEN);
        strncpy(Temp , (inputl+60) , 10);
        sscanf(Temp,"%f",&TBv);

/*
  sscanf(inputl," %d %d %s %s %f %f %f %s %d %f ",
  &AtNum,&TRs1,TResN,TAtmN,&TXc,&TYc,&TZc,TSegN,&TRs2,&TBv);

  fprintf(write_p,
  "%5d%5d %-4.4s %-4.4s%10.5f%10.5f%10.5f %4.4s %-4d%10.5f\n",
  ihelp,gomp_GetAtomResNum1( Which , i),gomp_GetAtomResName( Which , i) ,gomp_GetAtomAtmName(Which , i),
  gomp_GetAtomXCoord(Which , i)+Move[0], 
  gomp_GetAtomYCoord(Which , i)+Move[1],
  gomp_GetAtomZCoord(Which , i)+Move[2],
  gomp_GetAtomSegName( Which , i), gomp_GetAtomResNum1( Which , i),gomp_GetAtomBValue( Which , i));
  ihelp++;
  }

*/
        j = gomp_PutAtomSegName(Wstr , TSegN , i);
        j = gomp_PutAtomResName(Wstr , TResN , i);
        j = gomp_PutAtomAtmName(Wstr , TAtmN , i);

        j = gomp_PutAtomXCoord(Wstr , TXc , i);
        j = gomp_PutAtomYCoord(Wstr , TYc , i);
        j = gomp_PutAtomZCoord(Wstr , TZc , i);

        j = gomp_PutAtomResNum1(Wstr , TRs1 , i);
        j = gomp_PutAtomResNum2(Wstr , TRs2 , i);

        j = gomp_PutAtomBValue(Wstr , TBv , i);
        j = gomp_PutAtomCharge(Wstr ,  0.0 , i);

    }

    gomp_PrintMessage("**********   Done   **********");

end:
    fclose(chm_in);

/*     ShowDataStructure(0 , CharmmAtoms);  */

    return(Wstr >= 0 ? 0 : 1);
}
/*********************************************************************/
/*static int ShowDataStructure(int Alt , int Atoms)*/
/*********************************************************************/
/*{
    int i;

    for(i = 0 ; i < Atoms ; i++) {

        printf("'%s'%s'%s'%f,%f,%f\n",gomp_GetAtomSegName(Alt , i),
               gomp_GetAtomResName(Alt , i),
               gomp_GetAtomAtmName(Alt , i),
               gomp_GetAtomXCoord(Alt , i),
               gomp_GetAtomYCoord(Alt , i),
               gomp_GetAtomZCoord(Alt , i));
    }
 
    return(0);
}
*/

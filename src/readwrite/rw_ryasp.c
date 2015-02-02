/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Enhancements 2003, 2005 by:
Eero Häkkinen
*/
 
 
/*
  This program reads a YASP coordinate file
 
  Leif Laaksonen 1990, 1994
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

#define YASP_LINE_LEN   120   /* yasp file line length */
 
/*************************************************************************/
int gomp_ReadCoordinatesYASP(const char *Text1 , int Append)  
    /* routine for reading a yasp file */
/*************************************************************************/
{
    int   TRs1,TRs2;
    float TXc,TYc,TZc,TBv;
 
    char   input[YASP_LINE_LEN];
    char   temp[YASP_LINE_LEN];
    int    i,j,k;
    char   OutText[BUFF_LEN];
    int    YASPatoms;
    int    Wstr = 0;
 
    FILE *yasp_in;
 
    yasp_in=fopen(Text1,"r");
    if(yasp_in == NULL) {
        sprintf(OutText,"$Can't open input file : %s",Text1);
        gomp_PrintERROR(OutText);
        return(1);
    }
 
    TRs1 = 1;
    TRs2 = 1;
    TBv  = 0.0;
    TXc  = 0.0;
    TYc  = 0.0;
    TZc  = 0.0;

    sprintf(OutText,"********** Reading : %s **********",Text1);
    gomp_PrintMessage(OutText);
/*
  Start reading file
 
*/
 
    while(fgets(input,YASP_LINE_LEN,yasp_in) != NULL) {
        strcpy(temp,input);
        gomp_String2Lower(temp);
 
/* Title card */
 
        if(( k = gomp_Indexo(temp,"title:")) == 1) {
            fgets(input,YASP_LINE_LEN,yasp_in);
            gomp_PrintMessage(input);
        }
 
/* Coordinates card */
 
        if(( k = gomp_Indexo(temp,"coordinates:")) == 1) {
            fgets(input,YASP_LINE_LEN,yasp_in);
            sscanf(input,"%d",&YASPatoms);

/* it's possible to already update now */

            if(Append)
                Wstr = gomp_CreateMolecStruct(Text1 , YASPatoms , APPEND);
            else
                Wstr = gomp_CreateMolecStruct(Text1 , YASPatoms , NEW);
            if ( Wstr < 0 )
                goto end;

            for(i = 0 ; i < YASPatoms ; i++ ) {
                fgets(input,YASP_LINE_LEN,yasp_in);
                sscanf(input," %*d %s %f %f %f ",
                       temp,
                       &TXc ,&TYc ,&TZc );
 
                TXc *= 10.;
                TYc *= 10.;
                TZc *= 10.;
 
                j = gomp_PutAtomXCoord(Wstr , TXc , i);
                j = gomp_PutAtomYCoord(Wstr , TYc , i);
                j = gomp_PutAtomZCoord(Wstr , TZc , i);


/* truncate and put atom name in right place */
                k = strlen(temp);
                for(j = 1 ; j < k ; j++) {
                    if(temp[j] == 39) {
                        temp[j] = '\0';
                        break;
                    }
                }
 
                if(temp[0] == 39)
                    j = gomp_PutAtomAtmName(Wstr , (temp+1) , i);
                else
                    j = gomp_PutAtomAtmName(Wstr , temp , i);
/* add dummies */
                j = gomp_PutAtomSegName(Wstr , "YASP" , i);
                j = gomp_PutAtomResName(Wstr , "YASP" , i);

                j = gomp_PutAtomResNum1(Wstr , TRs1 , i);
                j = gomp_PutAtomResNum2(Wstr , TRs2 , i);

                j = gomp_PutAtomBValue(Wstr , TBv , i);
                j = gomp_PutAtomCharge(Wstr , 0.0 , i);
            }
        }
/* Coordinates/velocities card */
 
        if(( k = gomp_Indexo(temp,"coordinates/velocities:")) == 1) {
            fgets(input,YASP_LINE_LEN,yasp_in);
            sscanf(input,"%d",&YASPatoms);
 
/* it's possible to already update now */

            Wstr = gomp_CreateMolecStruct(Text1 , YASPatoms , NEW);
            if ( Wstr < 0 )
                goto end;

            for(i = 0 ; i < YASPatoms ; i++ ) {
                fgets(input,YASP_LINE_LEN,yasp_in);
                sscanf(input,"%*d %s %f %f %f",
                       temp,&TXc ,&TYc ,&TZc );
 
                TXc *= 10.;
                TYc *= 10.;
                TZc *= 10.;
 
                j = gomp_PutAtomXCoord(Wstr , TXc , i);
                j = gomp_PutAtomYCoord(Wstr , TYc , i);
                j = gomp_PutAtomZCoord(Wstr , TZc , i);

                fgets(input,YASP_LINE_LEN,yasp_in);  /* read velocities (not used) */
 
 
/* truncate and put atom name in right place */
                k = strlen(temp);
                for(j = 1 ; j < k ; j++) {
                    if(temp[j] == 39) {
                        temp[j] = '\0';
                        break;
                    }
                }
 
                if(temp[0] == 39)
                    j = gomp_PutAtomAtmName(Wstr , (temp+1) , i);
                else
                    j = gomp_PutAtomAtmName(Wstr , temp , i);
/* add dummies */
                j = gomp_PutAtomSegName(Wstr , "YASP" , i);
                j = gomp_PutAtomResName(Wstr , "YASP" , i);

                j = gomp_PutAtomResNum1(Wstr , TRs1 , i);
                j = gomp_PutAtomResNum2(Wstr , TRs2 , i);

                j = gomp_PutAtomBValue(Wstr , TBv , i);
                j = gomp_PutAtomCharge(Wstr , 0.0 , i);
            }
 
        }
 
/* basta card */
        if(( k = gomp_Indexo(temp,"basta:")) == 1) {
            gomp_PrintMessage("**********   Done   **********");
            return(0);
        }
 
    }
 
    gomp_PrintMessage("**********   Done   **********");
    
end:
    fclose(yasp_in);

    return(Wstr >= 0 ? 0 : 1);
}

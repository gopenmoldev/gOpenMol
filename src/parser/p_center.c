/*

Copyright (c) 1999 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <sys/types.h>

#include "gomstring.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "measure.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "selection.h"
#include "tclutils.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_CenterCommand(ClientData clientdata, Tcl_Interp *interp,
                       int argc, const char **argv)
/*********************************************************************/
{
    static int    i,j;
    static char   Text[BUFF_LEN];
    static char   Text1[BUFF_LEN];
    static char   Text2[BUFF_LEN];
    static char   Text3[BUFF_LEN];
    static float  xc;
    static float  yc;
    static float  zc;
    static float *x;
    static float *y;
    static float *z;
    static const float *sumxyz;
    static const float *temp;

/* #1   center  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);

    if(gomp_StringMatch(Text , "cent$er")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

/* system */
        if(gomp_StringMatch(Text , "syst$em")) {

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);


            sumxyz    = gomp_GetTranslateArray();

            if(!gomp_CalcCoordinateCenter(Text1 , Text2 , Text3 ,
                                        &xc , &yc , &zc)) {

                for(i = 0 ; i < gomp_GetNumMolecStructs() ; i++) {

                    if(!gomp_GetSelectedStructure(i)) continue;

                    x         = gomp_GetModifiableAtomXCoordPointer(i);
                    y         = gomp_GetModifiableAtomYCoordPointer(i);
                    z         = gomp_GetModifiableAtomZCoordPointer(i);

                    for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++) {
                        x[j]  = x[j] - xc + sumxyz[0];
                        y[j]  = y[j] - yc + sumxyz[1];
                        z[j]  = z[j] - zc + sumxyz[2];
                    }
                }
                (void)gomp_SaveTranslateArray( xc , yc , zc);
                return(TCL_OK);
            }
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "rota$tion")) {

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);


            sumxyz    = gomp_GetTranslateArray();

            if(!gomp_CalcCoordinateCenter(Text1 , Text2 , Text3 ,
                                        &xc , &yc , &zc)) {

                for(i = 0 ; i < gomp_GetNumMolecStructs() ; i++) {

                    if(!gomp_GetSelectedStructure(i)) continue;

                    temp = gomp_GetTranslateArrayMT(i);
                    (void)gomp_SaveTranslateArrayMT(    i , temp[0] + xc - sumxyz[0], 
                                                   temp[1] + yc - sumxyz[1], 
                                                   temp[2] + zc - sumxyz[2]);
                }
                return(TCL_OK);
            }
            else
                return(TCL_ERROR);
        }
        else {
            gomp_PrintERROR("command 'center' not recognized");
            return(TCL_ERROR);
        }

    }
/*  E R R O R command not recognized         */
    gomp_PrintERROR("'center' command not recognized");

    return(TCL_ERROR);

}


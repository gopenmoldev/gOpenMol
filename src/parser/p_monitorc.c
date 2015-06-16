/*

Copyright (c) 1994 - 2005 by:
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

#include "colouring.h"
#include "gommonitor.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "molecstruct.h"
#include "plot.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

#define ON  1
#define OFF 0

/*********************************************************************/
int gomp_MonitorCommand(ClientData clientdata, Tcl_Interp *interp,
                      int argc, const char **argv)
/*********************************************************************/
{
    static char Text[BUFF_LEN];
    static char Text1[BUFF_LEN];
    static char Text2[BUFF_LEN];
    static char Text3[BUFF_LEN];
    static char Text4[BUFF_LEN];
    static char Text5[BUFF_LEN];
    static char Text6[BUFF_LEN];
    static char Text7[BUFF_LEN];
    static char Text8[BUFF_LEN];
    static char Text9[BUFF_LEN];
    static char Text10[BUFF_LEN];
    static char Text11[BUFF_LEN];
    static char Text12[BUFF_LEN];
    static char Text13[BUFF_LEN];
    static char Text14[BUFF_LEN];

/* #1   monitor ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "moni$tor")) {

        if(!gomp_GetNumMolecStructs()) {
            gomp_PrintERROR("no atom structures defined");
            return(TCL_ERROR);
        }

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

        if(     gomp_StringMatch(Text , "dist$ance"))  {

            strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text4,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text5,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text6,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text7,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text8,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);

/* check first if the intention is to change the colour */
            if(gomp_StringMatch(Text1 , "colo$ur") ||
               gomp_StringMatch(Text1 , "colo$r")) {
                {
                    float RedC;
                    float GreenC;
                    float BlueC;

                    if(strlen(Text2) == 0) {
                        gomp_PrintERROR("line colour value/nale is missing");
                        return(TCL_ERROR);
                    }

                    if(gomp_ColourName2RGB(Text2 , &RedC , &GreenC , &BlueC)) {
                        sprintf(Text4,"can't resolve colour name '%s'",Text2);
                        gomp_PrintERROR(Text4);
                        return(TCL_ERROR);
                    }

                    if(strlen(Text3) == 0) {
		      strcpy(Text3, "1"); 
                    }

                    if(!gomp_SetDistMonitorColor((atoi(Text3) - 1), 
                                               RedC, GreenC, BlueC))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
            }
/* check first if the intention is to change the type */
            if(gomp_StringMatch(Text1 , "type")) {
                {

		  if(strlen(Text2) == 0) {
                        gomp_PrintERROR("line type value/nale is missing");
                        return(TCL_ERROR);
                    }

		  if(strlen(Text3) == 0) {
		    strcpy(Text3, "1");
                    }

                    if(!gomp_SetDistMonitorType((atoi(Text3) - 1), atoi(Text2)))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
            }

            if(gomp_SelectDistArray(
                   Text1, 
                   Text2,
                   Text3,
                   Text4,
                   Text5,
                   Text6,
                   Text7,
                   Text8))
                return(TCL_ERROR);

            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "-dis$tance")) {
            gomp_ResetMonitorDistanceData();
            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "angl$e"))     {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text4,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text5,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text6,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text7,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text8,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text9,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text10,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text11,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);

/* check first if the intention is to change the colour */
            if(gomp_StringMatch(Text1 , "colo$ur") ||
               gomp_StringMatch(Text1 , "colo$r")) {
                {
                    float RedC;
                    float GreenC;
                    float BlueC;

                    if(strlen(Text2) == 0) {
                        gomp_PrintERROR("line colour value/nale is missing");
                        return(TCL_ERROR);
                    }

                    if(gomp_ColourName2RGB(Text2 , &RedC , &GreenC , &BlueC)) {
                        sprintf(Text4,"can't resolve colour name '%s'",Text2);
                        gomp_PrintERROR(Text4);
                        return(TCL_ERROR);
                    }

                    if(strlen(Text3) == 0) {
		      strcpy(Text3, "1");
                    }

                    if(!gomp_SetAngMonitorColor((atoi(Text3) - 1), 
                                              RedC, GreenC, BlueC))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
            }

/* check first if the intention is to change the type */
            if(gomp_StringMatch(Text1 , "type")) {
                {
		  if(strlen(Text2) == 0) {
                        gomp_PrintERROR("line type value/nale is missing");
                        return(TCL_OK);
                    }

		  if(strlen(Text3) == 0) {
		    strcpy(Text3, "1");
		  }

                    if(!gomp_SetAngMonitorType((atoi(Text3) - 1), atoi(Text2)))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
            }
            if(gomp_SelectAngArray(
                   Text1, 
                   Text2,
                   Text3,
                   Text4,
                   Text5,
                   Text6,
                   Text7,
                   Text8,
                   Text9,
                   Text10,
                   Text11))
                return(TCL_ERROR);

            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "-ang$le"))    {
            gomp_ResetMonitorAngleData();
            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "tors$ion") ||
                gomp_StringMatch(Text , "dihe$dral"))   {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text4,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text5,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text6,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text7,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text8,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text9,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text10,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text11,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text12,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text13,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            strncpy(Text14,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);

/* check first if the intention is to change the colour */
            if(gomp_StringMatch(Text1 , "colo$ur") ||
               gomp_StringMatch(Text1 , "colo$r")) {
                {
                    float RedC;
                    float GreenC;
                    float BlueC;

                    if(strlen(Text2) == 0) {
                        gomp_PrintERROR("line colour value/nale is missing");
                        return(TCL_ERROR);
                    }

                    if(gomp_ColourName2RGB(Text2 , &RedC , &GreenC , &BlueC)) {
                        sprintf(Text4,"can't resolve colour name '%s'",Text2);
                        gomp_PrintERROR(Text4);
                        return(TCL_ERROR);
                    }

                    if(strlen(Text3) == 0) {
		      strcpy(Text3, "1");
                    }

                    if(!gomp_SetTorsMonitorColor((atoi(Text3) - 1), 
                                               RedC, GreenC, BlueC))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
            }

/* check first if the intention is to change the type */
            if(gomp_StringMatch(Text1 , "type")) {
                {
		  if(strlen(Text2) == 0) {
                        gomp_PrintERROR("line type value/nale is missing");
                        return(TCL_ERROR);
                    }

		  if(strlen(Text3) == 0) {
		    strcpy(Text3, "1");
                    }

                    if(!gomp_SetTorsMonitorType((atoi(Text3) - 1), atoi(Text2)))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
            }
            if(gomp_SelectTorsArray(
                   Text1, 
                   Text2,
                   Text3,
                   Text4, 
                   Text5,
                   Text6,
                   Text7,
                   Text8,
                   Text9,
                   Text10,
                   Text11,
                   Text12,
                   Text13,
                   Text14))
                return(TCL_ERROR);

            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "-tor$sion") ||
                gomp_StringMatch(Text , "-dih$edral"))  {
            gomp_ResetMonitorTorsionData();
            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "disp$lay"))  {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text,"dist$ance")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(gomp_StringMatch(Text1,"on")) {
                    (void)gomp_SetDistMonitor(ON);
                    return(TCL_OK);
                }
                else if(gomp_StringMatch(Text1,"off")) {
                    (void)gomp_SetDistMonitor(OFF);
                    return(TCL_OK);
                }

                else {
                    gomp_PrintERROR("unknown state in 'monitor display distance' command");
                    return(TCL_ERROR);
                }
            }
/* .... */
            else if(gomp_StringMatch(Text,"angl$e")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(gomp_StringMatch(Text1,"on")) {
                    (void)gomp_SetAngMonitor(ON);
                    return(TCL_OK);
                }
                else if(gomp_StringMatch(Text1,"off")) {
                    (void)gomp_SetAngMonitor(OFF);
                    return(TCL_OK);
                }

                else {
                    gomp_PrintERROR("unknown state in 'monitor display distance' command");
                    return(TCL_ERROR);
                }
            }
/* .... */
            else if(gomp_StringMatch(Text,"tors$ion") ||
                    gomp_StringMatch(Text,"dihe$dral")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(gomp_StringMatch(Text1,"on")) {
                    (void)gomp_SetTorsMonitor(ON);
                    return(TCL_OK);
                }
                else if(gomp_StringMatch(Text1,"off")) {
                    (void)gomp_SetTorsMonitor(OFF);
                    return(TCL_OK);
                }

                else {
                    gomp_PrintERROR("unknown state in 'monitor display distance' command");
                    return(TCL_ERROR);
                }
            }
        }
        else if(gomp_StringMatch(Text , "edit"))  {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text , "dist$ance"))  {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(gomp_EditDistVector( atoi(Text1) - 1 , 
                                      atoi(Text2) - 1 , atoi(Text3) - 1 , 
                                      atoi(Text4)     , atof(Text5) , atof(Text6) , atof(Text7))) {
                    return(TCL_ERROR);
                } else {
                    return(TCL_OK);
                }
            } else if(gomp_StringMatch(Text , "angl$e"))  {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(gomp_EditAngVector( atoi(Text1) - 1 , 
                                     atoi(Text2) - 1 , atoi(Text3) - 1 , atoi(Text4) - 1 ,
                                     atoi(Text5)     , atof(Text6) , atof(Text7) , atof(Text8))) {
                    return(TCL_ERROR);
                } else {
                    return(TCL_OK);
                }
            } else if(gomp_StringMatch(Text , "tors$ion"))  {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(gomp_EditTorsVector( atoi(Text1) - 1 , 
                                      atoi(Text2) - 1 , atoi(Text3) - 1 , atoi(Text4) - 1 , atoi(Text5) - 1 ,
                                      atoi(Text6)     , atof(Text7) , atof(Text8) , atof(Text9))) {
                    return(TCL_ERROR);
                } else {
                    return(TCL_OK);
                }
            } else {
                gomp_PrintERROR("unknow 'monitor edit' parameter");
                return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "dele$te"))  {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("index number to list is missing");
                return(TCL_ERROR);
            }
            if(gomp_StringMatch(Text , "dist$ance"))  {
                if(gomp_DeleteDistanceIndex(atoi(Text1)))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);

            } else if(gomp_StringMatch(Text , "angl$e"))  {
                if(gomp_DeleteAngleIndex(atoi(Text1)))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            } else if(gomp_StringMatch(Text , "torsion"))  {
                if(gomp_DeleteTorsionIndex(atoi(Text1)))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            } else {
                gomp_PrintERROR("unknown 'monitor delete' parameter");
                return(TCL_ERROR);
            }
        }
        else {
            gomp_PrintERROR("'monitor' command not recognized");
            return(TCL_ERROR);
        }

    }

/*  E R R O R command not recognized         */
    return(TCL_ERROR);
}


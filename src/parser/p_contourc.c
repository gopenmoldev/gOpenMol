/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2004 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tcl.h>
#include <sys/types.h>
#include <ctype.h>

#include "colouring.h"
#include "contour.h"
#include "gomfile.h"
#include "gomstring.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "plot.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

#define TRIGGER 1.e-06

static int ParseContourNameList(const char * , const char * , const char *);

/*********************************************************************/
int gomp_ContourCommand(ClientData clientdata, Tcl_Interp *interp,
                      int argc, const char **argv)
/*********************************************************************/
{
    static char  Text[BUFF_LEN];
    static char  Text1[BUFF_LEN];
    static char  Text2[BUFF_LEN];
    static char  Text3[BUFF_LEN];
    static char  Text4[BUFF_LEN];
    static char  Text5[BUFF_LEN];
    static char  Text6[BUFF_LEN];
    static char  Text7[BUFF_LEN];
    static char  Text8[BUFF_LEN];
    static char  Text9[BUFF_LEN];
    static char  Text10[BUFF_LEN];
    static int   ITemp;
    static int   ITemp1;
    static int   ITemp2;
    static int   Loop;
    static int   i;
    static float FTemp;
    static float FTemp1;
    static float FTemp2;
    static float FTemp3;
    static float FTempO;
    static float FTempO1;
    static float FTempO2;
    static float FTempO3;
    static const char *TValue;
    static int   iContourType;
    static float fContourOpaque;
    static int   iContourSmooth;

/*   contour  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "conto$ur")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

/* file */
        if(gomp_StringMatch(Text , "file")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text1[0] == (char)NULL) {
                gomp_PrintERROR("file name is missing");
                return(TCL_ERROR);
            }

            if(gomp_FileNameIsURL(Text1)) {
                sprintf(Text1,"?Can't handle the given URL '%s'",Text1);
                gomp_PrintERROR(Text1);
                return(TCL_ERROR);
            }

            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text2[0] == (char)NULL) {
                sprintf(Text2,"%d",(gomp_GetContoursDefined()+1));
                sprintf(Text,"will use default name '%s' for the contour",Text2);
                gomp_PrintMessage(Text);
            }

            ITemp = gomp_ContourDriver(Text1 , Text2);

            if(!ITemp)
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-fil$e")) {
            if(!gomp_DeleteAllContours())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "dele$te")) {
            if(!gomp_DeleteAllContours())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "info")) {
            if(!gomp_ShowContourInfo())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "meth$od")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_StringMatch(Text , "dire$ct")) {
                (void)gomp_SetSurfaceMethod(0);
            }
            else if(gomp_StringMatch(Text , "save")) {
                (void)gomp_SetSurfaceMethod(1);
            }
            else {
                gomp_PrintERROR("wrong method (has to be 'direct' or 'save'");
                return(TCL_ERROR);
            }
            return(TCL_OK);
        }
/* plot */
        else if(gomp_StringMatch(Text , "plot")) {
            int OldContourLevels;

            gomp_CopyString(Text ,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(Text[0] == (char)NULL) {
                gomp_PrintERROR("contour name missing");
                return(TCL_ERROR);
            }

/* check into the contour name text to see if there are two names in the
   list. This means that the grid data of the second name will be projected
   on the first one.
*/
            if(ParseContourNameList(Text, Text1 , Text2)) {
                return(TCL_OK);
            }

            if(!(ITemp = gomp_CheckContourName(Text))) {
                sprintf(Text,"no contour name match '%s'",Text);
                gomp_PrintERROR(Text);
                return(TCL_ERROR);
            }

/* reset levels before start ...*/
/* seve the old number of contour levels ... */
            OldContourLevels = gomp_GetContourLevels(ITemp - 1);

            (void)gomp_SetContourProjection(ITemp - 1, ITemp - 1);

/* type */
            TValue = Tcl_GetVar(gomp_GetTclInterp(),"gomDefaultContourType",TCL_GLOBAL_ONLY);

            iContourType = CONTOUR_TYPE_SOLID; /* default */

            if(TValue) {
                if(isdigit(TValue[0])) {
                    if(atoi(TValue)) {
                        iContourType = CONTOUR_TYPE_MESH;
                    } else {
                        iContourType = CONTOUR_TYPE_SOLID;
                    }
                } else {
                    if(Tcl_StringCaseMatch(TValue, "solid", 1)) {
                        iContourType = CONTOUR_TYPE_SOLID;
                    } else if(Tcl_StringCaseMatch(TValue, "mesh", 1)) {
                        iContourType = CONTOUR_TYPE_MESH;
                    } else if(gomp_StringMatch(TValue , "line")) {
                        iContourType = CONTOUR_TYPE_LINE;
                    } else {
                        gomp_PrintERROR("wrong option (has to be 'solid' or 'mesh'");
                        return(TCL_ERROR);
                    }
                }
            } 

/* opaque */

            TValue = Tcl_GetVar(gomp_GetTclInterp(),"gomDefaultContourOpaque",TCL_GLOBAL_ONLY);

            fContourOpaque = 1.0; /* default */

            if(TValue) {
                fContourOpaque = atof(TValue);
            } 

/* smooth */
            TValue = Tcl_GetVar(gomp_GetTclInterp(),"gomDefaultContourSmooth",TCL_GLOBAL_ONLY);

            iContourSmooth = 0; /* default */

            if(TValue) {
                if(isdigit(TValue[0])) {
                    if(atoi(TValue)) {
                        iContourSmooth = 1;
                    } else {
                        iContourSmooth = 0;
                    }
                } else {
                    if(Tcl_StringCaseMatch(TValue, "on", 1)) {
                        iContourSmooth = 1;
                    } else if(Tcl_StringCaseMatch(TValue, "off", 1)) {
                        iContourSmooth = 0;
                    } else {
                        gomp_PrintERROR("wrong option (has to be 'on' or 'off'");
                        return(TCL_ERROR);
                    }
                }
            } 

            Loop = 0;
            while(Text1[0] != (char)NULL) {

                FTemp      = atof(Text1);
                if(gomp_ColourName2RGB(Text2 , &FTemp1 , &FTemp2 , &FTemp3)) {
                    gomp_FinalizeContourLevels(ITemp - 1,Loop);
                    sprintf(Text,"can't resolve colour name '%s'",Text2);
                    gomp_PrintERROR(Text2);
                    return(TCL_ERROR);
                }
/* check the old defined values */
                if(Loop < OldContourLevels) {

                    FTempO  = gomp_GetContourValue(ITemp - 1 , Loop);
                    FTempO1 = gomp_GetContourColourRed(ITemp - 1 , Loop);
                    FTempO2 = gomp_GetContourColourGreen(ITemp - 1 , Loop);
                    FTempO3 = gomp_GetContourColourBlue(ITemp - 1 , Loop);
                    
                    if((fabs(FTemp  - FTempO)  > TRIGGER) ||
                       (fabs(FTemp1 - FTempO1) > TRIGGER) ||
                       (fabs(FTemp2 - FTempO2) > TRIGGER) ||
                       (fabs(FTemp3 - FTempO3) > TRIGGER)) {
/* update values ... */
                        (void)gomp_ParseContourLevels(ITemp - 1,
                                                    FTemp ,
                                                    FTemp1,FTemp2,FTemp3,
                                                    Loop);
                    }
                } else {
                    (void)gomp_ParseContourLevels(ITemp - 1,
                                                FTemp ,
                                                FTemp1,FTemp2,FTemp3,
                                                Loop);

                    if(iContourSmooth)
                        (void) gomp_ContourSmoothON(ITemp - 1 , Loop);
                    else
                        (void) gomp_ContourSmoothOFF(ITemp - 1 , Loop);

                    (void) gomp_ContourDisplayON(ITemp - 1 , Loop);
                    (void) gomp_SetContourDisplayType(ITemp - 1 ,
                                                    Loop      ,
                                                    iContourType);
                    (void) gomp_SetContourAlpha(ITemp - 1 , Loop , fContourOpaque);
                    (void) gomp_SetContourCullFace(ITemp - 1 , Loop , 0);
                }
                Loop++;
                (void) gomp_SetSurfaceControlON();
                strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                        BUFF_LEN-1);
                strncpy(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                        BUFF_LEN-1);
            }
            (void)gomp_FinalizeContourLevels(ITemp - 1,Loop);

            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "alph$ablend")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(!(ITemp = gomp_CheckContourName(Text))) {
                sprintf(Text,"no contour name match '%s'",Text);
                gomp_PrintERROR(Text);
                return(TCL_ERROR);
            }
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            FTemp = atof(Text);
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] == (char)NULL) {
                if(!gomp_GetContourLevels(ITemp - 1)) {
                    gomp_PrintERROR("no contour levels defined");
                    return(TCL_ERROR);
                }
                for(i = 0 ; i < gomp_GetContourLevels(ITemp - 1) ; i++) 
                    (void)gomp_SetContourAlpha(ITemp - 1,i,FTemp);
                return(TCL_OK);
            }
            ITemp1 = atoi(Text);
            if(ITemp1 < 1 || ITemp1 > gomp_GetContourLevels(ITemp - 1)) {
                gomp_PrintERROR("given contour level is not in the allowed range");
                return(TCL_ERROR);
            }

            if(!gomp_SetContourAlpha(ITemp - 1,ITemp1 - 1,FTemp))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "smoo$th")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(!(ITemp = gomp_CheckContourName(Text))) {
                sprintf(Text,"no contour name match '%s'",Text);
                gomp_PrintERROR(Text);
                return(TCL_ERROR);
            }
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            ITemp1 = 0;
            if(gomp_StringMatch(Text , "on")) {
                ITemp1 = 1;
            }
            else if(gomp_StringMatch(Text , "off")) {
                ITemp1 = 0;
            }
            else {
                gomp_PrintERROR("wrong state (has to be 'on' or 'off'");
                return(TCL_ERROR);
            }
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] == (char)NULL) {
                if(!gomp_GetContourLevels(ITemp - 1)) {
                    gomp_PrintERROR("no contour levels defined");
                    return(TCL_ERROR);
                }
                for(i = 0 ; i < gomp_GetContourLevels(ITemp - 1) ; i++) {
                    if(ITemp1)
                        (void)gomp_ContourSmoothON( ITemp - 1 , i);
                    else
                        (void)gomp_ContourSmoothOFF(ITemp - 1 , i);
                }
                return(TCL_OK);
            }
            ITemp2 = atoi(Text);
            if(ITemp2 < 1 || ITemp2 > gomp_GetContourLevels(ITemp - 1)) {
                gomp_PrintERROR("given contour level is not in the allowed range");
                return(TCL_ERROR);
            }

            if(ITemp1)
                (void)gomp_ContourSmoothON( ITemp - 1 , ITemp2 - 1);
            else
                (void)gomp_ContourSmoothOFF(ITemp - 1 , ITemp2 - 1);
            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "cull$face")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(!(ITemp = gomp_CheckContourName(Text))) {
                sprintf(Text,"no contour name match '%s'",Text);
                gomp_PrintERROR(Text);
                return(TCL_ERROR);
            }
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            ITemp1 = 0;
            if(gomp_StringMatch(Text , "on")) {
                ITemp1 = 1;
            }
            else if(gomp_StringMatch(Text , "off")) {
                ITemp1 = 0;
            }
            else {
                gomp_PrintERROR("wrong state (has to be 'on' or 'off'");
                return(TCL_ERROR);
            }
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] == (char)NULL) {
                if(!gomp_GetContourLevels(ITemp - 1)) {
                    gomp_PrintERROR("no contour levels defined");
                    return(TCL_ERROR);
                }
                for(i = 0 ; i < gomp_GetContourLevels(ITemp - 1) ; i++) {
                    if(ITemp1)
                        (void)gomp_SetContourCullFace( ITemp - 1 , i , 1);
                    else
                        (void)gomp_SetContourCullFace( ITemp - 1 , i , 0);
                }
                return(TCL_OK);
            }
            ITemp2 = atoi(Text);
            if(ITemp2 < 1 || ITemp2 > gomp_GetContourLevels(ITemp - 1)) {
                gomp_PrintERROR("given contour level is not in the allowed range");
                return(TCL_ERROR);
            }

            if(ITemp1)
                (void)gomp_SetContourCullFace( ITemp - 1 , ITemp2 - 1 , 1);
            else
                (void)gomp_SetContourCullFace( ITemp - 1 , ITemp2 - 1 , 0);
            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "disp$lay")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(!(ITemp = gomp_CheckContourName(Text))) {
                sprintf(Text,"no contour name match '%s'",Text);
                gomp_PrintERROR(Text);
                return(TCL_ERROR);
            }
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            ITemp1 = 0;
            if(gomp_StringMatch(Text , "on")) {
                ITemp1 = 1;
            }
            else if(gomp_StringMatch(Text , "off")) {
                ITemp1 = 0;
            }
            else {
                gomp_PrintERROR("wrong state (has to be 'on' or 'off'");
                return(TCL_ERROR);
            }
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] == (char)NULL) {
                if(!gomp_GetContourLevels(ITemp - 1)) {
                    gomp_PrintERROR("no contour levels defined");
                    return(TCL_ERROR);
                }
                for(i = 0 ; i < gomp_GetContourLevels(ITemp - 1) ; i++) {
                    if(ITemp1)
                        (void)gomp_ContourDisplayON( ITemp - 1 , i);
                    else
                        (void)gomp_ContourDisplayOFF(ITemp - 1 , i);
                }
                return(TCL_OK);
            }
            ITemp2 = atoi(Text);
            if(ITemp2 < 1 || ITemp2 > gomp_GetContourLevels(ITemp - 1)) {
                gomp_PrintERROR("given contour level is not in the allowed range");
                return(TCL_ERROR);
            }
            if(ITemp1)
                (void)gomp_ContourDisplayON( ITemp - 1 , ITemp2 - 1);
            else
                (void)gomp_ContourDisplayOFF(ITemp - 1 , ITemp2 - 1);
            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "type")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(!(ITemp = gomp_CheckContourName(Text))) {
                sprintf(Text,"no contour name match '%s'",Text);
                gomp_PrintERROR(Text);
                return(TCL_ERROR);
            }
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            ITemp1 = CONTOUR_TYPE_SOLID;

            if(gomp_StringMatch(Text , "soli$d")) {
                ITemp1 = CONTOUR_TYPE_SOLID;
            }
            else if(gomp_StringMatch(Text , "mesh")) {
                ITemp1 = CONTOUR_TYPE_MESH;
            }
            else if(gomp_StringMatch(Text , "line")) {
                ITemp1 = CONTOUR_TYPE_LINE;
            }
            else {
                gomp_PrintERROR("wrong state (has to be 'mesh' or 'solid'");
                return(TCL_ERROR);
            }
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] == (char)NULL) {
                if(!gomp_GetContourLevels(ITemp - 1)) {
                    gomp_PrintERROR("no contour levels defined");
                    return(TCL_ERROR);
                }
                for(i = 0 ; i < gomp_GetContourLevels(ITemp - 1) ; i++) {
                    (void)gomp_SetContourDisplayType(ITemp - 1 , i , ITemp1);
                }
                return(TCL_OK);
            }
            ITemp2 = atoi(Text);
            if(ITemp2 < 1 || ITemp2 > gomp_GetContourLevels(ITemp - 1)) {
                gomp_PrintERROR("given contour level is not in the allowed range");
                return(TCL_ERROR);
            }
            (void)gomp_SetContourDisplayType(ITemp - 1 , ITemp2 - 1 , ITemp1);
            
            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "mapp$ing")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(Text1[0] == (char)NULL) {
                gomp_PrintERROR("contour name missing for #1");
                return(TCL_ERROR);
            }

            if(!(ITemp1 = gomp_CheckContourName(Text1))) {
                sprintf(Text,"no contour name match '%s'",Text1);
                gomp_PrintERROR(Text1);
                return(TCL_ERROR);
            }
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(Text2[0] == (char)NULL) {
                gomp_PrintERROR("contour name missing for #2 (map this on #1)");
                return(TCL_ERROR);
            }

            if(!(ITemp2 = gomp_CheckContourName(Text2))) {
                sprintf(Text,"no contour name match '%s'",Text2);
                gomp_PrintERROR(Text2);
                return(TCL_ERROR);
            }

            (void)gomp_SetContourProjection(ITemp1 - 1, ITemp2 - 1);

            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "comb$ine")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(Text1[0] == (char)NULL) {
                gomp_PrintERROR("contour name missing for #1");
                return(TCL_ERROR);
            }

            if(!(ITemp1 = gomp_CheckContourName(Text1))) {
                sprintf(Text,"no contour name match '%s'",Text1);
                gomp_PrintERROR(Text1);
                return(TCL_ERROR);
            }
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(Text2[0] == (char)NULL) {
                gomp_PrintERROR("structure number missing");
                return(TCL_ERROR);
            }

            ITemp2 = atoi(Text2);

            if(gomp_SetContour2StructureMapping(ITemp1 - 1, ITemp2 - 1)) {
                return(TCL_ERROR);
            }

            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "clip$plane")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text1[0] == (char)NULL) {
                gomp_PrintERROR("clipping plane parameter(s) missing");
                return(TCL_ERROR);
            }
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text2[0] == (char)NULL) {
                gomp_PrintERROR("clipping plane option parameters missing");
                return(TCL_ERROR);
            }
            if ( Text1[1] == '\0' ) {
                switch ( Text1[0] ) {
                case 'x':
                case 'X': ITemp2 = 1; break;
                case 'y':
                case 'Y': ITemp2 = 2; break;
                case 'z':
                case 'Z': ITemp2 = 3; break;
                default: goto error;
                }
                if(gomp_StringMatch(Text2 , "on")) {
                    (void)gomp_SetContourClippingPlaneState( ITemp2 , ON);
                    return(TCL_OK);
                } else if(gomp_StringMatch(Text2 , "off")) {
                    (void)gomp_SetContourClippingPlaneState( ITemp2 , OFF);
                    return(TCL_OK);
                }
                if(!(ITemp = gomp_CheckContourName(Text2))) {
                    gomp_FormatERROR("no contour name match '%s'", Text2);
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                ITemp1 = atoi(Text4);
                if ( gomp_SetContourLevelClippingPlaneAxis(ITemp - 1, ITemp1 - 1, "xyz"[ITemp2 - 1]) ||
                     gomp_SetContourLevelClippingPlanePosition(ITemp - 1, ITemp1 - 1, atof(Text3)) )
                    return(TCL_ERROR);
                return(TCL_OK);
/* x */
            } else if(gomp_StringMatch(Text1 , "x+") ||
                      gomp_StringMatch(Text1 , "X+")) {
                (void)gomp_SetContourClippingPlaneParameters( 1 , '+', atof(Text2));
                return(TCL_OK);
            } else if(gomp_StringMatch(Text1 , "x-") ||
                      gomp_StringMatch(Text1 , "X-")) {
                (void)gomp_SetContourClippingPlaneParameters( 1 , '-', atof(Text2));
                return(TCL_OK);
/* y */
            } else if(gomp_StringMatch(Text1 , "y+") ||
                      gomp_StringMatch(Text1 , "Y+")) {
                (void)gomp_SetContourClippingPlaneParameters( 2 , '+', atof(Text2));
                return(TCL_OK);
            } else if(gomp_StringMatch(Text1 , "y-") ||
                      gomp_StringMatch(Text1 , "Y-")) {
                (void)gomp_SetContourClippingPlaneParameters( 2 , '-', atof(Text2));
                return(TCL_OK);
/* z */
            } else if(gomp_StringMatch(Text1 , "z+") ||
                      gomp_StringMatch(Text1 , "Z+")) {
                (void)gomp_SetContourClippingPlaneParameters( 3 , '+', atof(Text2));
                return(TCL_OK);
            } else if(gomp_StringMatch(Text1 , "z-") ||
                      gomp_StringMatch(Text1 , "Z-")) {
                (void)gomp_SetContourClippingPlaneParameters( 3 , '-', atof(Text2));
                return(TCL_OK);
/* xyz 1 - 3 */
            } else if(gomp_StringMatch(Text1 , "xyz1") ||
                      gomp_StringMatch(Text1 , "XYZ1")) {
                if(gomp_StringMatch(Text2 , "on")) {
                    (void)gomp_SetContourXYZClippingPlaneState( 1 , ON);
                    return(TCL_OK);
                } else if(gomp_StringMatch(Text2 , "off")) {
                    (void)gomp_SetContourXYZClippingPlaneState( 1 , OFF);
                    return(TCL_OK);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text10,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                (void)gomp_SetContourXYZClippingPlaneParameters( 1 , atof(Text2) , atof(Text3) , atof(Text4) , 
                                                               atof(Text5) , atof(Text6) , atof(Text7) ,
                                                               atof(Text8) , atof(Text9) , atof(Text10));

                return(TCL_OK);
            } else if(gomp_StringMatch(Text1 , "xyz2") ||
                      gomp_StringMatch(Text1 , "XYZ2")) {
                if(gomp_StringMatch(Text2 , "on")) {
                    (void)gomp_SetContourXYZClippingPlaneState( 2 , ON);
                    return(TCL_OK);
                } else if(gomp_StringMatch(Text2 , "off")) {
                    (void)gomp_SetContourXYZClippingPlaneState( 2 , OFF);
                    return(TCL_OK);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text10,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                (void)gomp_SetContourXYZClippingPlaneParameters( 2 , atof(Text2) , atof(Text3) , atof(Text4) , 
                                                               atof(Text5) , atof(Text6) , atof(Text7) ,
                                                               atof(Text8) , atof(Text9) , atof(Text10));

                return(TCL_OK);
            } else if(gomp_StringMatch(Text1 , "xyz3") ||
                      gomp_StringMatch(Text1 , "XYZ3")) {
                if(gomp_StringMatch(Text2 , "on")) {
                    (void)gomp_SetContourXYZClippingPlaneState( 3 , ON);
                    return(TCL_OK);
                } else if(gomp_StringMatch(Text2 , "off")) {
                    (void)gomp_SetContourXYZClippingPlaneState( 3 , OFF);
                    return(TCL_OK);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text10,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                (void)gomp_SetContourXYZClippingPlaneParameters( 3 , atof(Text2) , atof(Text3) , atof(Text4) , 
                                                               atof(Text5) , atof(Text6) , atof(Text7) ,
                                                               atof(Text8) , atof(Text9) , atof(Text10));

                return(TCL_OK);
            } else {
                gomp_PrintERROR("parameter to turn clipping plane has to be x, y, z, x+, x-, y+, y-, z+, z-, xyz1, xyz2 or xyz3");
                return(TCL_ERROR);
            }
        }
        else {
            gomp_PrintERROR("'contour' command not recognized");
            return(TCL_ERROR);
        }
    }
/*  E R R O R command not recognized         */
error:
    gomp_PrintERROR("'contour' command not recognized");

    return(TCL_ERROR);

}

/*********************************************************************/
int ParseContourNameList(const char *TextI , const char *TextI1 , const char *TextI2)
/*********************************************************************/
{
    static int argc, code;
    const char **argv;

    static char  Text[BUFF_LEN];
    static char  Text1[BUFF_LEN];
    static char  Text2[BUFF_LEN];
    static char  Text3[BUFF_LEN];
    static int   ITemp;
    static int   OldContourLevels;
    static int   ITemp1;
    static int   Loop;
    static float FTemp;
    static float FTemp1;
    static float FTemp2;
    static float FTempO;
    static float FTempO1;
    static float FTempO2;

    code = Tcl_SplitList(gomp_GetTclInterp(), TextI, &argc, &argv);

/* only one name coming in */
    if(argc == 1) {  
        Tcl_Free((char *)CONST_CAST(char **, argv));
        return(0);
    }

    if(argc > 2) {
        gomp_PrintERROR("too many entries ( > 2) in the list. Will take first two");
    }

    if(!(ITemp = gomp_CheckContourName(argv[0]))) {
        sprintf(Text,"Entry #1: no contour name match '%s'",argv[0]);
        gomp_PrintERROR(Text);
        return(1);
    }
    if(!(ITemp1 = gomp_CheckContourName(argv[1]))) {
        sprintf(Text,"Entry #2: no contour name match '%s'",argv[1]);
        gomp_PrintERROR(Text);
        return(1);
    }

    Tcl_Free((char *)CONST_CAST(char **, argv));

/* reset levels before start ...*/
/* seve the old number of contour levels ... */
    OldContourLevels = gomp_GetContourLevels(ITemp - 1);

    (void)gomp_SetContourProjection(ITemp - 1, ITemp1 - 1);

    gomp_CopyString(Text1,TextI1,BUFF_LEN);
    gomp_CopyString(Text2,TextI2,BUFF_LEN);
    strncpy(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),
            BUFF_LEN-1);
    Loop = 0;
    while(Text1[0] != (char)NULL) {

        FTemp      = atof(Text1);

        if(*Text2 != (char)NULL) {
            FTemp1 = atof(Text2);
        } else {
            FTemp1 = gomp_GetContourMin(ITemp1 - 1);
        }
        if(*Text3 != (char)NULL) {
            FTemp2 = atof(Text3);
        } else {
            FTemp2 = gomp_GetContourMax(ITemp1 - 1);
        }

/* check the old defined values */
        if(Loop < OldContourLevels) {

            FTempO  = gomp_GetContourValue(ITemp - 1 , Loop);
            FTempO1 = gomp_GetContourProjectionMin(ITemp - 1 , Loop);
            FTempO2 = gomp_GetContourProjectionMax(ITemp - 1 , Loop);
            
            if((fabs(FTemp  - FTempO)  > TRIGGER) ||
               (fabs(FTemp1 - FTempO1) > TRIGGER) ||
               (fabs(FTemp2 - FTempO2) > TRIGGER)) {
/* update values ... */
                (void)gomp_ParseContourLevelsProjection(
                    ITemp - 1,
                    FTemp ,
                    FTemp1,FTemp2,
                    Loop);
            }
        } else {
            (void)gomp_ParseContourLevelsProjection(
                ITemp - 1,
                FTemp ,
                FTemp1,FTemp2,
                Loop);
            (void) gomp_ContourSmoothOFF(ITemp - 1 , Loop);
            (void) gomp_ContourDisplayON(ITemp - 1 , Loop);
            (void) gomp_SetContourDisplayType(ITemp - 1 ,
                                            Loop      ,
                                            CONTOUR_TYPE_SOLID);
            (void) gomp_SetContourAlpha(ITemp - 1 , Loop , 1.0);
        }
        Loop++;
        (void) gomp_SetSurfaceControlON();

        strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                BUFF_LEN-1);
        strncpy(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                BUFF_LEN-1);
        strncpy(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                BUFF_LEN-1);
    }
    (void)gomp_FinalizeContourLevels(ITemp - 1 , Loop);

    return(1);
}

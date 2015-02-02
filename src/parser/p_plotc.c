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
#include <math.h>

#include "axis.h"
#include "cell.h"
#include "cluster.h"
#include "colouring.h"
#include "contour.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "gomtext.h"
#include "ldp.h"
#include "objseg.h"
#include "plot.h"
#include "printmsg.h"
#include "rforce.h"
#include "tclutils.h"

#include "stdafx.h"

#define LDP_ON  1
#define LDP_OFF 0

#define CLUSTER_ON  1
#define CLUSTER_OFF 0

#define VECTOR_ON   1
#define VECTOR_OFF  0

/*********************************************************************/
int gomp_PlotCommand(ClientData clientdata, Tcl_Interp *interp,
                   int argc, const char **argv)
/*********************************************************************/
{
    static int   ITemp;
    static int   ITemp1;
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
    static char  Text11[BUFF_LEN];
    static char  Text12[BUFF_LEN];
    static char  Text13[BUFF_LEN];
    static char  Text14[BUFF_LEN];
    static char  Text15[BUFF_LEN];
    static char  Text16[BUFF_LEN];
    static char  Text17[BUFF_LEN];
    static char  Text18[BUFF_LEN];
    static char  Text19[BUFF_LEN];
    static char  Text20[BUFF_LEN];
    static char  Text21[BUFF_LEN];
    static char  Text22[BUFF_LEN];
    static char  Text23[BUFF_LEN];
    static float FTemp;
    static float FTemp1;
    static float FTemp2;
    static float FTemp3;

/* #1   plot  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "plot")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

/* sphere */
        if(gomp_StringMatch(Text , "sphe$re")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text10,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
      
            if(!gomp_PushSphereStack(Text1,
                                Text2,
                                Text3,
                                Text4,
                                Text5,
                                Text6,
                                Text7,
                                Text8,
                                Text9,
                                Text10))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-sph$ere")) {
            (void)gomp_DelSphereSeg();
            return(TCL_OK);
        }
/* contour cut plane */
        else if(gomp_StringMatch(Text , "cutp$lane")) {
            if(!gomp_GetContoursDefined()) {
                gomp_PrintERROR("no contour file defined defined");
                return(TCL_ERROR);
            }
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(gomp_StringMatch(Text1 , "x") ||
               gomp_StringMatch(Text1 , "X")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(!(ITemp = gomp_CheckContourName(Text2))) {
                    sprintf(Text,"no contour name match '%s'",Text);
                    gomp_PrintERROR(Text);
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text3[0] == (char)NULL) {
                    gomp_PrintERROR("x-value for plane is undefined");
                    return(TCL_ERROR);
                }
                ITemp--;
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(!gomp_PrepareCutPlaneX(ITemp,atof(Text3),Text4,Text5,Text6))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "y") ||
                    gomp_StringMatch(Text1 , "Y")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(!(ITemp = gomp_CheckContourName(Text2))) {
                    sprintf(Text,"no contour name match '%s'",Text);
                    gomp_PrintERROR(Text);
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text3[0] == (char)NULL) {
                    gomp_PrintERROR("y-value for plane is undefined");
                    return(TCL_ERROR);
                }
                ITemp--;
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(!gomp_PrepareCutPlaneY(ITemp,atof(Text3),Text4,Text5,Text6))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "z") ||
                    gomp_StringMatch(Text1 , "Z")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(!(ITemp = gomp_CheckContourName(Text2))) {
                    sprintf(Text,"no contour name match '%s'",Text);
                    gomp_PrintERROR(Text);
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text3[0] == (char)NULL) {
                    gomp_PrintERROR("z-value for plane is undefined");
                    return(TCL_ERROR);
                }
                ITemp--;
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(!gomp_PrepareCutPlaneZ(ITemp,atof(Text3),Text4,Text5,Text6))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "xyz1") ||
                    gomp_StringMatch(Text1 , "XYZ1")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(!(ITemp = gomp_CheckContourName(Text2))) {
                    sprintf(Text,"no contour name match '%s'",Text);
                    gomp_PrintERROR(Text);
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text3[0] == (char)NULL) {
                    gomp_PrintERROR("xyz-value for plane is undefined");
                    return(TCL_ERROR);
                }
                ITemp--;
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text10,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text11,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text12,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text13,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                ITemp1 = gomp_PrepareCutPlaneXYZ( 1 , ITemp,
                                                Text3,Text4,Text5,
                                                Text6,Text7,Text8,
                                                Text9,Text10,Text11,
                                                Text12,Text13);
                if(!ITemp1) {
                    (void) gomp_SetSurfaceControlON();
                    return(TCL_OK);
                } else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "xyz2") ||
                    gomp_StringMatch(Text1 , "XYZ2")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(!(ITemp = gomp_CheckContourName(Text2))) {
                    sprintf(Text,"no contour name match '%s'",Text);
                    gomp_PrintERROR(Text);
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text3[0] == (char)NULL) {
                    gomp_PrintERROR("xyz-value for plane is undefined");
                    return(TCL_ERROR);
                }
                ITemp--;
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text10,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text11,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text12,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text13,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                ITemp1 = gomp_PrepareCutPlaneXYZ( 2 , ITemp,
                                                Text3,Text4,Text5,
                                                Text6,Text7,Text8,
                                                Text9,Text10,Text11,
                                                Text12,Text13);
                if(!ITemp1) {
                    (void) gomp_SetSurfaceControlON();
                    return(TCL_OK);
                } else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "xyz3") ||
                    gomp_StringMatch(Text1 , "XYZ3")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(!(ITemp = gomp_CheckContourName(Text2))) {
                    sprintf(Text,"no contour name match '%s'",Text);
                    gomp_PrintERROR(Text);
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text3[0] == (char)NULL) {
                    gomp_PrintERROR("xyz-value for plane is undefined");
                    return(TCL_ERROR);
                }
                ITemp--;
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text10,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text11,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text12,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text13,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                ITemp1 = gomp_PrepareCutPlaneXYZ( 3 , ITemp,
                                                Text3,Text4,Text5,
                                                Text6,Text7,Text8,
                                                Text9,Text10,Text11,
                                                Text12,Text13);
                if(!ITemp1) { 
                    (void) gomp_SetSurfaceControlON();
                    return(TCL_OK);
                } else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "3dx") ||
                    gomp_StringMatch(Text1 , "3DX")) {
                (void)gomp_SetCutPlaneType_3D_X(1);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "3dy") ||
                    gomp_StringMatch(Text1 , "3DY")) {
                (void)gomp_SetCutPlaneType_3D_Y(1);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "3dz") ||
                    gomp_StringMatch(Text1 , "3DZ")) {
                (void)gomp_SetCutPlaneType_3D_Z(1);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "2dx") ||
                    gomp_StringMatch(Text1 , "2DX")) {
                (void)gomp_SetCutPlaneType_3D_X(0);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "2dy") ||
                    gomp_StringMatch(Text1 , "2DY")) {
                (void)gomp_SetCutPlaneType_3D_Y(0);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "2dz") ||
                    gomp_StringMatch(Text1 , "2DZ")) {
                (void)gomp_SetCutPlaneType_3D_Z(0);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "spec$trum") ||
                    gomp_StringMatch(Text1 , "prof$ile")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

                if(!(ITemp = gomp_CheckContourName(Text2))) {
                    sprintf(Text,"no contour name match '%s'",Text);
                    gomp_PrintERROR(Text);
                    return(TCL_ERROR);
                }

                if(Text3[0] == (char)NULL) {
                    gomp_PrintERROR("please define x-,y- or z-axis");
                    return(TCL_ERROR);
                }

                ITemp--;

                if(gomp_StringMatch(Text3 , "x") || gomp_StringMatch(Text3 , "x")) {
                    if(!gomp_PrepareCutPlaneSpectrumX(ITemp,Text4,Text5,Text6,Text7))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else if(gomp_StringMatch(Text3 , "y") || gomp_StringMatch(Text3 , "Y")) {
                    if(!gomp_PrepareCutPlaneSpectrumY(ITemp,Text4,Text5,Text6,Text7))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else if(gomp_StringMatch(Text3 , "z") || gomp_StringMatch(Text3 , "Z")) {
                    if(!gomp_PrepareCutPlaneSpectrumZ(ITemp,Text4,Text5,Text6,Text7))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else {
                    gomp_PrintERROR("please define x-,y- or z-axis");
                    return(TCL_ERROR);
                }
            }
            else {
                gomp_PrintERROR("undefinied 'plot cutplane' option");
                return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "-cutp$lane")) {
            if(!gomp_GetContoursDefined()) {
                gomp_PrintERROR("no contour file defined defined");
                return(TCL_ERROR);
            }
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(gomp_StringMatch(Text1 , "x") ||
               gomp_StringMatch(Text1 , "X")) {
                if(!gomp_DeleteCutPlaneDataX())
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            if(gomp_StringMatch(Text1 , "y") ||
               gomp_StringMatch(Text1 , "Y")) {
                if(!gomp_DeleteCutPlaneDataY())
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            if(gomp_StringMatch(Text1 , "z") ||
               gomp_StringMatch(Text1 , "Z")) {
                if(!gomp_DeleteCutPlaneDataZ())
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            if(gomp_StringMatch(Text1 , "xyz1") ||
               gomp_StringMatch(Text1 , "XYZ1")) {
                if(!gomp_DisableCutPlanePlotStateXYZ(1))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            if(gomp_StringMatch(Text1 , "xyz2") ||
               gomp_StringMatch(Text1 , "XYZ2")) {
                if(!gomp_DisableCutPlanePlotStateXYZ(2))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            if(gomp_StringMatch(Text1 , "xyz3") ||
               gomp_StringMatch(Text1 , "XYZ3")) {
                if(!gomp_DisableCutPlanePlotStateXYZ(3))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else {
                (void)gomp_DeleteCutPlaneDataX();
                (void)gomp_DeleteCutPlaneDataY();
                (void)gomp_DeleteCutPlaneDataZ();
                (void)gomp_DisableCutPlanePlotStateXYZ(1);
                (void)gomp_DisableCutPlanePlotStateXYZ(2);
                (void)gomp_DisableCutPlanePlotStateXYZ(3);
                return(TCL_OK);
            }
        }
/* cylinder */
        else if(gomp_StringMatch(Text , "cyli$nder")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text10,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(!gomp_PushCylinderStack(Text1,
                                  Text2,
                                  Text3,
                                  Text4,
                                  Text5,
                                  Text6,
                                  Text7,
                                  Text8,
                                  Text9,
                                  Text10))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-cyl$inder")) {
            (void)gomp_DelCylinderSeg();
            return(TCL_OK);
        }
/* arrow */
        else if(gomp_StringMatch(Text , "arro$w")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text10,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(!gomp_PushArrowStack(Text1,
                               Text2,
                               Text3,
                               Text4,
                               Text5,
                               Text6,
                               Text7,
                               Text8,
                               Text9,
                               Text10))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-arro$w")) {
            (void)gomp_DelArrowSeg();
            return(TCL_OK);
        }
/* line */
        else if(gomp_StringMatch(Text , "line")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(!gomp_PushLineStack(Text1,
                              Text2,
                              Text3,
                              Text4,
                              Text5,
                              Text6,
                              Text7,
                              Text8))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-lin$e")) {
            (void)gomp_DelLineSeg();
            return(TCL_OK);
        }
/* plane */
        else if(gomp_StringMatch(Text , "plan$e")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(!gomp_PushPlaneStack(Text1,
                               Text2,
                               Text3,
                               Text4,
                               Text5,
                               Text6,
                               Text7))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-pla$ne")) {
            (void)gomp_DelPlaneSeg();
            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "ldp")) {
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(gomp_StringMatch(Text4 , "atom$s")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(!gomp_ApplyLdpSelectionMaskLine(Text1,Text2,Text3))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text4 , "on")) {
                if(!gomp_SetDisplayLDPmatrix(LDP_ON))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text4 , "off")) {
                if(!gomp_SetDisplayLDPmatrix(LDP_OFF))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else {
                gomp_PrintERROR("unrecognized 'plot ldp' command");
                return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "-clu$ster")) {
            if(!gomp_SetDisplayCLUSTERmatrix(CLUSTER_OFF))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "clus$ter")) {

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_StringMatch(Text1 , "on")) {
                if(!gomp_SetDisplayCLUSTERmatrix(CLUSTER_ON))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "off")) {
                if(!gomp_SetDisplayCLUSTERmatrix(CLUSTER_OFF))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else {
                gomp_PrintERROR("unrecognized 'plot cluster' command");
                return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "vect$or")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_StringMatch(Text1 , "on")) {
                if(!gomp_SetPlotVectorStatus(VECTOR_ON))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "off")) {
                if(!gomp_SetPlotVectorStatus(VECTOR_OFF))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }

            else if(gomp_StringMatch(Text1 , "atom$s")) {
                if(!gomp_ParsePlotVectorList(Text2,Text3,Text4,Text5,Text6))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "scal$e")) {
                if(Text2[0] == (char)NULL) {
                    gomp_PrintERROR("scale value missing");
                    return(TCL_ERROR);
                } else {
                    FTemp = atof(Text2);
                }
                (void)gomp_SetVectorScale(FTemp);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "rang$e")) {
                {
                    float RangeMin;
                    float RangeMax;

                    if(Text2[0] == (char)NULL) {
                        gomp_PrintERROR("min range value missing");
                        return(TCL_ERROR);
                    } else {
                        RangeMin = atof(Text2);
                    }

                    if(Text3[0] == (char)NULL) {
                        RangeMax = 1.0e+25f;
                    } else {
                        RangeMax = atof(Text3);
                    }


                    if(!gomp_SetVectorDisplayRange(RangeMin , RangeMax))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
            }
            else if(gomp_StringMatch(Text1 , "byco$lours") || 
                    gomp_StringMatch(Text1 , "byco$lors")) {
                if(!gomp_SetGradientDisplayStyle(1))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "-byc$olours") || 
                    gomp_StringMatch(Text1 , "-byc$olors")) {
                if(!gomp_SetGradientDisplayStyle(0))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else {
                gomp_PrintERROR("unrecognized 'plot vector' command");
                return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "-vec$tor")) {
            (void)gomp_DeleteVectorStructure();
            if(!gomp_SetPlotVectorStatus(VECTOR_OFF))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* triangles */
        else if(gomp_StringMatch(Text , "tria$ngles")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text10,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text11,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text12,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text13,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text14,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text15,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text16,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text17,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text18,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text19,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text20,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text21,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text22,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text23,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(!gomp_PushTriangleStack(Text1,
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
                                  Text14,
                                  Text15,
                                  Text16,
                                  Text17,
                                  Text18,
                                  Text19,
                                  Text20,
                                  Text21,
                                  Text22,
                                  Text23))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-tri$angles")) {
            (void)gomp_DelTriangleSeg();
            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "csca$le")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            (void)gomp_SetColourScalePlotStatus(ON);

            if(Text1[0] != (char)NULL) {
                gomp_SetColourScalePlotLevels(atoi(Text1));
            }
            if(Text2[0] != (char)NULL) {
                gomp_SetColourScalePlotMin(atof(Text2));
            }
            if(Text3[0] != (char)NULL) {
                gomp_SetColourScalePlotMax(atof(Text3));
            }

            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "-csc$ale")) {
            if(!gomp_SetColourScalePlotStatus(OFF))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "cell")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(gomp_StringMatch(Text1 , "on")) {
                if(!gomp_SetPlotCell(ON))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "off")) {
                if(!gomp_SetPlotCell(OFF))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else {
                gomp_PrintERROR("plot cell command must be 'on' or 'off'");
                return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "text3")) {
/* text colour */
            strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
            if(Text1[0] == (char)NULL) {
                gomp_PrintERROR("colour name is missing");
                return(TCL_ERROR);
            }
            strncpy(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
/* x coordinate */
            if(Text2[0] == (char)NULL) {
                gomp_PrintERROR("x coordinate is missing");
                return(TCL_ERROR);
            }
            strncpy(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
/* y coordinate */
            if(Text3[0] == (char)NULL) {
                gomp_PrintERROR("y coordinate  is missing");
                return(TCL_ERROR);
            }
            strncpy(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
/* z coordinate */
            if(Text4[0] == (char)NULL) {
                gomp_PrintERROR("z coordinate  is missing");
                return(TCL_ERROR);
            }
            strncpy(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
/* text         */
            if(Text5[0] == (char)NULL) {
                gomp_PrintERROR("text is missing");
                return(TCL_ERROR);
            }

            if(gomp_ColourName2RGB(Text1,&FTemp1,&FTemp2,&FTemp3)) {
                gomp_PrintERROR("can't process the supplied colour");
                return(TCL_ERROR);
            }
            if(!gomp_PushText2AnnotateStack(3 , Text5,"*",
                                          FTemp1,FTemp2,FTemp3,
                                          (float)atof(Text2),(float)atof(Text3), (float)atof(Text4)))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* text and text2 */
        else if(gomp_StringMatch(Text , "text")  ||
                gomp_StringMatch(Text , "text2")) {
/* text colour */
            strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
            if(Text1[0] == (char)NULL) {
                gomp_PrintERROR("colour name is missing");
                return(TCL_ERROR);
            }
            strncpy(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
/* x coordinate */
            if(Text2[0] == (char)NULL) {
                gomp_PrintERROR("x coordinate (0...1) is missing");
                return(TCL_ERROR);
            }
            strncpy(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
/* y coordinate */
            if(Text3[0] == (char)NULL) {
                gomp_PrintERROR("y coordinate (0...1) is missing");
                return(TCL_ERROR);
            }
            strncpy(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
/* text         */
            if(Text4[0] == (char)NULL) {
                gomp_PrintERROR("text is missing");
                return(TCL_ERROR);
            }

            if(gomp_ColourName2RGB(Text1,&FTemp1,&FTemp2,&FTemp3)) {
                gomp_PrintERROR("can't process the supplied colour");
                return(TCL_ERROR);
            }
            if(!gomp_PushText2AnnotateStack(2 , Text4,"*",
                                          FTemp1,FTemp2,FTemp3,
                                          (float)atof(Text2),(float)atof(Text3), (float)0.0))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-tex$t")) {
            if(gomp_DeleteTextStack())
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "axis")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(!gomp_ParsePlotAxisList(Text1 , Text2 , Text3 ,
                                     Text4 , Text5 , Text6))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-axi$s")) {
            if(!gomp_DeletePlotAxisList())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else {
            gomp_PrintERROR("command 'plot' not recognized");
            return(TCL_ERROR);
        }

    }
/*  E R R O R command not recognized         */
    gomp_PrintERROR("'plot' command not recognized");

    return(TCL_ERROR);

}


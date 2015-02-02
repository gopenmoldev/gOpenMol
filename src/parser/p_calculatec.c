/*

Copyright (c) 1995 - 2005 by:
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

#include "bond.h"
#include "cluster.h"
#include "gomclipbrd.h"
#include "gommonitor.h"
#include "gomstring.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "math_oper.h"
#include "measure.h"
#include "molecstruct.h"
#include "msd.h"
#include "plot.h"
#include "plumber.h"
#include "printmsg.h"
#include "rdf.h"
#include "tclutils.h"
#include "trajectory.h"

#include "stdafx.h"

#define DIST_TYPE  0
#define ANG_TYPE   1
#define TORS_TYPE  2

/*********************************************************************/
int gomp_CalculateCommand(ClientData clientdata, Tcl_Interp *interp,
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
    static int  ITemp;
    static float FTemp;

/* calculate command */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "calc$ulate")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

        if(gomp_StringMatch(Text , "test")) {

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyBitmap2Clipboard();
            return(TCL_OK);
        }
/* connectivity */
        if(gomp_StringMatch(Text , "conn$ectivity")) {

            if(!gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("no molecular structure(s) defined");
                return(TCL_ERROR);
            }

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_ParseCalcConnList(Text1 , Text2 , Text3 , "Dummy"))
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
/* h-bonds */
        else if(gomp_StringMatch(Text , "hbon$ds")) {
            if(!gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("no molecular structure(s) defined");
                return(TCL_ERROR);
            }

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_StringMatch(Text5 , "nohy$drogens")) {
                if(gomp_ParseCalcHbondsList(Text1 , Text2 , Text3 , 0 , Text4))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
            else if (gomp_StringMatch(Text4 , "nohy$drogens")) {
                if(gomp_ParseCalcHbondsList(Text1 , Text2 , Text3 , 0 , ""))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
            else if (gomp_StringMatch(Text5 , "hydr$ogens")) {
                if(gomp_ParseCalcHbondsList(Text1 , Text2 , Text3 , 1 , Text4))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
            else if (gomp_StringMatch(Text4 , "hydr$ogens")) {
                if(gomp_ParseCalcHbondsList(Text1 , Text2 , Text3 , 1 , ""))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
            else {
                if(gomp_ParseCalcHbondsList(Text1 , Text2 , Text3 , 1 , Text4))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
        }
        else if(gomp_StringMatch(Text , "-hbo$nds")) {
            if(!gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("no molecular structure(s) defined");
                return(TCL_ERROR);
            }

            (void)gomp_DeleteAllHbonds();
 
            return(TCL_OK);
        }
/* geometrical center */
        else if(gomp_StringMatch(Text , "geomc$enter")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            {
                float xc = 1.e+20f;
                float yc = 1.e+20f;
                float zc = 1.e+20f;
                ITemp = gomp_CalcCoordinateCenter(
                    Text1,
                    Text2,
                    Text3,
                    &xc, &yc ,&zc);
                if(xc > 1.e+18 &&
                   yc > 1.e+18 &&
                   zc > 1.e+18)
                    *Text = '\0';
                else
                    sprintf(Text,"%f %f %f",xc,yc,zc);
/*        gomp_PrintMessage(Text);*/
                gomp_SendTclReturn(Text);

                if(ITemp)
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
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
            gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            {
                float A = 0.0;
                float B = 0.0;
                float C = 0.0;
                float D = 0.0;
                float p1[3];
                float p2[3];
                float p3[3];

                p1[0] = atof(Text1);
                p1[1] = atof(Text2);
                p1[2] = atof(Text3);

                p2[0] = atof(Text4);
                p2[1] = atof(Text5);
                p2[2] = atof(Text6);

                p3[0] = atof(Text7);
                p3[1] = atof(Text8);
                p3[2] = atof(Text9);

                ITemp = gomp_CalcPlaneFrom3Points(p1, p2 , p3 , 
                                                &A , &B , &C , &D);

                sprintf(Text,"%f %f %f %f",A , B , C , D);
/*        gomp_PrintMessage(Text);*/
                gomp_SendTclReturn(Text);

                if(ITemp)
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
        }

/* center of mass */
        else if(gomp_StringMatch(Text , "massc$enter")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            {
                float xc = 1.e+20f;
                float yc = 1.e+20f;
                float zc = 1.e+20f;
                ITemp = gomp_CalculateCenterofMass(
                    Text1,
                    Text2,
                    Text3,
                    &xc, &yc ,&zc);
                if(xc > 1.e+18 &&
                   yc > 1.e+18 &&
                   zc > 1.e+18) 
                    *Text = '\0';
                else
                    sprintf(Text,"%f %f %f",xc,yc,zc);
/*        gomp_PrintMessage(Text);*/
                gomp_SendTclReturn(Text);
                if(ITemp)
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
        }
/* average structure */
        else if(gomp_StringMatch(Text , "avst$ructure")) {
            if(gomp_TrajAvStructure())
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
/* cluster */
        else if(gomp_StringMatch(Text , "clus$ter")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(gomp_CalcCluster(Text1,Text2,Text3,
                              Text1,Text2,Text3))
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
/* -cluster */
        else if(gomp_StringMatch(Text , "-clu$ster")) {
            if(gomp_DeleteClusterData())
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
/* distance */
        else if(gomp_StringMatch(Text , "dist$ance")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            sprintf(Text,"%f",
                    gomp_CalculateDistance(Text1 , Text2 , Text3 ,
                                         Text4 , Text5 , Text6));
/*        gomp_PrintMessage(Text);*/
            if(gomp_SendTclReturn(Text))
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
/* angle   */
        else if(gomp_StringMatch(Text , "angl$e")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            sprintf(Text,"%f",
                    gomp_CalculateAngle(Text1 , Text2 , Text3 ,
                                      Text4 , Text5 , Text6 ,
                                      Text7 , Text8 , Text9));
/*        gomp_PrintMessage(Text);*/
            if(gomp_SendTclReturn(Text))
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
/* torsion */
        else if(gomp_StringMatch(Text , "tors$ion")) {
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

            sprintf(Text,"%f",
                    gomp_CalculateTorsionAngle(Text1 , Text2 , Text3 ,
                                             Text4 , Text5 , Text6 ,
                                             Text7 , Text8 , Text9 ,
                                             Text10, Text11, Text12));
/*        gomp_PrintMessage(Text);*/
            if(gomp_SendTclReturn(Text))
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "rmsf$luctuation")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_StringMatch(Text4 , "fit") || Text4[0] == (char)NULL) {
                if(gomp_RMS_Fluctuation( 1 , Text1,Text2,Text3))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
            else if(gomp_StringMatch(Text4 , "nofi$t")) {
                if(gomp_RMS_Fluctuation( 0 , Text1,Text2,Text3))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
            else {
                gomp_PrintERROR("unknown parameter in command");
                return(TCL_ERROR);
            }
        }
/* mean square displacement */
        else if(gomp_StringMatch(Text , "msdi$splacement")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_StringMatch(Text1 , "atom")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(gomp_MeanSquareDisplacement( 0 , Text2,Text3,Text4))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "mass$center")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(gomp_MeanSquareDisplacement( 1 , Text2,Text3,Text4))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
            else {
                gomp_PrintERROR("unrecognized 'calculate msdisplacement' command");
                return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "rdf")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(gomp_ParseRDFList(Text1, Text2, Text3, 
                               Text4, Text5, Text6,
                               Text7, Text8))
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "-rdf")) {
            if(gomp_DeleteRDF())
                return(TCL_OK);
            else
                return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "rdfm$ean")) {
            if(gomp_CalcMeanRDF())
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "bbto$rsion") ||
                gomp_StringMatch(Text , "bbdi$hedrals")) {
            if(gomp_SecondaryStructure())
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }   
        else if(gomp_StringMatch(Text , "corr$elation")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(gomp_StringMatch(Text1 , "dist$ance")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text2[0] == (char)NULL) {
                    gomp_PrintERROR("index value is missing");
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text3[0] == (char)NULL)
                    gomp_CopyString(Text3,Text2,BUFF_LEN);
                if(gomp_CalculateCorrelation(DIST_TYPE, atoi(Text2), atoi(Text3)))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "angl$e")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text2[0] == (char)NULL) {
                    gomp_PrintERROR("index value is missing");
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text3[0] == (char)NULL)
                    gomp_CopyString(Text3,Text2,BUFF_LEN);
                if(gomp_CalculateCorrelation(ANG_TYPE, atoi(Text2), atoi(Text3)))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "tors$ion")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text2[0] == (char)NULL) {
                    gomp_PrintERROR("index value is missing");
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text3[0] == (char)NULL)
                    gomp_CopyString(Text3,Text2,BUFF_LEN);
                if(gomp_CalculateCorrelation(TORS_TYPE, atoi(Text2), atoi(Text3)))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
        }
        else if(gomp_StringMatch(Text , "quat$fit")) {
            if(gomp_GetNumMolecStructs() < 2) {
                gomp_PrintERROR("number of structures must be >= 2");
                (void)gomp_SendTclReturn("ERROR: number of structures must be >= 2");
                return(TCL_ERROR);
            }
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
          
            ITemp = 1;

            if(gomp_StringMatch(Text9 , "on")) {
                ITemp = 1;
            }
            else if(gomp_StringMatch(Text9 , "off")) {
                ITemp = 0;
            }
            else if(Text9[0] == (char)NULL) {
                ITemp = 1;
            }
            else {
                gomp_PrintERROR("option has to be 'on' or 'off'");
                return(TCL_ERROR);
            }
         
            FTemp = gomp_CalculateQuatfit(Text1 , Text2 , Text3 , Text4 ,
                                        Text5 , Text6 , Text7 , Text8 , ITemp);
            sprintf(Text9,"%f",FTemp);
/*         gomp_PrintMessage(Text9);*/
            (void)gomp_SendTclReturn(Text9);
            return(TCL_OK);
        }
    }

/*  E R R O R command not recognized         */
    gomp_PrintERROR("problems in 'calculate' command"); 

    return(TCL_ERROR);
}


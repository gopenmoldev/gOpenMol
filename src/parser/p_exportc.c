/*

Copyright (c) 1995 - 2005 by:
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
#include <tcl.h>
#include <math.h>

#include "cluster.h"
#include "coord_file.h"
#include "correl.h"
#include "gommonitor.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "model_file.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "msd.h"
#include "printmsg.h"
#include "projview.h"
#include "rdf.h"
#include "tclutils.h"

#include "stdafx.h"

#define GET_COMMAND_LENGTH 5

static int WriteCoordinates(int Which ,
                            const char *Text1 ,
                            const char *Text2,
                            const char *Text3);

/*********************************************************************/
int gomp_ExportCommand(ClientData clientdata, Tcl_Interp *interp, 
                     int argc, const char **argv)
/*********************************************************************/
{
    static char  Text[BUFF_LEN];
    static char  Text1[BUFF_LEN];
    static char  Text2[BUFF_LEN];
    static char  Text3[BUFF_LEN];
    static char  Text4[BUFF_LEN];
    static int   ITemp;
    static const char *value;

    float Xmin , Xmax;
    float Ymin , Ymax;
    float Zmin , Zmax;
    int   Xpts , Ypts, Zpts;
    float ProbeRad;

    const float *sumxyz;

    int    gomp_InputView;

/* #1   export ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "expo$rt")) {

/* #1.1 coor$dinates */
        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

        if(gomp_StringMatch(Text , "coor$dinates")) {

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            ITemp = atoi(Text1);
            if(ITemp < 1 || ITemp > gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("structure number out of range");
                return(TCL_ERROR);
            }
       
            if(!WriteCoordinates(ITemp , Text2 , Text3 , Text4))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else  if(gomp_StringMatch(Text , "dict$ionary")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            if(!gomp_ImportDictionaryAndApply(Text1))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* #1.2 rmsfluctuation */
        else  if(gomp_StringMatch(Text , "rmsf$luctuation")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("file name missing");
                return(TCL_ERROR);
            }
            if(!gomp_ExportMSFset(Text))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* #1.3 rdf */
        else if(gomp_StringMatch(Text , "rdf")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("rdf file name is missing");
                return(TCL_ERROR);
            }
            if(!gomp_WriteRDF(Text1))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* means square displacement */
        else if(gomp_StringMatch(Text , "msdi$splacement")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("file name is missing");
                return(TCL_ERROR);
            }
            if(!gomp_MSDdataWrite(Text1))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* #1.4 mode$l */
        else if(gomp_StringMatch(Text , "mode$l")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("file name missing");
                return(TCL_ERROR);
            }
            if(!gomp_WriteOldModel(Text))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* cluster */
        else if(gomp_StringMatch(Text , "clus$ter")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("cluster file name missing");
                return(TCL_ERROR);
            }
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text2) == 0) {
                gomp_CopyString(Text2,"** Default cluster file label **",BUFF_LEN);
            }
            if(!gomp_WriteClusterData(Text1 , Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* correlation data */
        else if(gomp_StringMatch(Text , "corr$elation")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("correlation file name missing");
                return(TCL_ERROR);
            }
            if(!gomp_WriteCorrelationArray(Text1))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* time series */
        else if(gomp_StringMatch(Text , "dist$ance")) {
            strncpy(Text,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            if(gomp_StringMatch(Text , "list") ||
               gomp_StringMatch(Text , "seri$es")) {
                strncpy(Text,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text) == 0) {
                    gomp_PrintERROR("index to list is missing");
                    return(TCL_ERROR);
                }
                ITemp = atoi(Text);
                gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text) == 0) {
                    gomp_PrintERROR("file name is missing");
                    return(TCL_ERROR);
                }
                if(!gomp_WriteDistanceVector(ITemp,Text))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "angl$e")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text , "list") ||
               gomp_StringMatch(Text , "seri$es")) {
                gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                ITemp = atoi(Text);
                gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text) == 0) {
                    gomp_PrintERROR("file name is missing");
                    return(TCL_ERROR);
                }
                if(!gomp_WriteAngleVector(ITemp,Text))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "tors$ion")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text , "list") ||
               gomp_StringMatch(Text , "seri$es")) {
                gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                ITemp = atoi(Text);
                gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text) == 0) {
                    gomp_PrintERROR("file name is missing");
                    return(TCL_ERROR);
                }
                if(!gomp_WriteTorsionVector(ITemp,Text))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
        }
/* protein backbone torsion */
        else if(gomp_StringMatch(Text , "bbto$rsion") ||
                gomp_StringMatch(Text , "bbdi$hedrals")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("backbone torsion file name missing");
                return(TCL_ERROR);
            }
            if(!gomp_WriteBackboneDihedrals(Text1))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* input */
        else if(gomp_StringMatch(Text , "input")) {
            if(gomp_GetNumMolecStructs() < 1) {
                gomp_PrintERROR("no structure(s) available");
                return(TCL_ERROR);
            }

            value  = Tcl_GetVar(gomp_GetTclInterp() , 
                                        "lulInputView", 
                                        TCL_GLOBAL_ONLY);

            if(!value) {
                gomp_InputView = 0;
            } else {
                if(atoi(value) > 0) {
                    gomp_InputView = 1;
                } else {
                    gomp_InputView = 0;
                }
            }

            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
/* .. USER */
            if(gomp_StringMatch(Text , "user")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("structure number missing");
                    return(TCL_ERROR);
                }
                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("input file name is missing");
                    return(TCL_ERROR);
                }

                if(!gomp_InputView) {
                    if(gomp_PushAtomCoordinates())
                        return(TCL_ERROR);
                    if(gomp_DoViewingTransformationOverStructures( -1 )) {
                        gomp_PopAtomCoordinates();
                        return(TCL_ERROR);
                    }
                }

                sprintf(Text2,"lulMakeUSERinput %d {%s}",atoi(Text2),Text1);

                ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text2);
                if(ITemp != TCL_OK) {
                    gomp_PrintERROR("can't execute the make USER input script");
                    return(TCL_ERROR);
                }

                if(!gomp_InputView)
                    gomp_PopAtomCoordinates();

                return(TCL_OK);
            }

/* .. ICON8 */
            else if(gomp_StringMatch(Text , "icon$8")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("structure number missing");
                    return(TCL_ERROR);
                }
                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("input file is missing");
                    return(TCL_ERROR);
                }
                if(!gomp_InputView) {
                    if(gomp_PushAtomCoordinates())
                        return(TCL_ERROR);
                    if(gomp_DoViewingTransformationOverStructures( -1 )) {
                        gomp_PopAtomCoordinates();
                        return(TCL_ERROR);
                    }
                }

                if(!gomp_ExternalInput4ICON8(atoi(Text2) - 1 , Text1)) {
                    if(!gomp_InputView)
                        gomp_PopAtomCoordinates();
                    return(TCL_OK);
                } else {
                    if(!gomp_InputView)
                        gomp_PopAtomCoordinates();
                    return(TCL_ERROR);
                }
            }
/* .. GAMESS */ 
            else if(gomp_StringMatch(Text , "game$ss")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("structure number missing");
                    return(TCL_ERROR);
                }
                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("input file name is missing");
                    return(TCL_ERROR);
                }

                if(!gomp_InputView) {
                    if(gomp_PushAtomCoordinates())
                        return(TCL_ERROR);
                    if(gomp_DoViewingTransformationOverStructures( -1 )) {
                        gomp_PopAtomCoordinates();
                        return(TCL_ERROR);
                    }
                }

                sprintf(Text2,"lulMakeGAMESSinput %d {%s}",atoi(Text2) , Text1);

                ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text2);
                if(ITemp != TCL_OK) {
                    if(!gomp_InputView)
                        gomp_PopAtomCoordinates();
                    gomp_PrintERROR("can't execute the make GAMESS input script");
                    return(TCL_ERROR);
                }
                if(!gomp_InputView)
                    gomp_PopAtomCoordinates();
                return(TCL_OK);
            }
/* .. MOPAC */
            else if(gomp_StringMatch(Text , "mopa$c")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("structure number missing");
                    return(TCL_ERROR);
                }
                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("input file name is missing");
                    return(TCL_ERROR);
                }

                if(!gomp_InputView) {
                    if(gomp_PushAtomCoordinates())
                        return(TCL_ERROR);
                    if(gomp_DoViewingTransformationOverStructures( -1 )) {
                        gomp_PopAtomCoordinates();
                        return(TCL_ERROR);
                    }
                }

                sprintf(Text2,"lulMakeMOPACinput %d {%s}",atoi(Text2),Text1);

                ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text2);
                if(ITemp != TCL_OK) {
                    if(!gomp_InputView)
                        gomp_PopAtomCoordinates();
                    gomp_PrintERROR("can't execute the make MOPAC input script");
                    return(TCL_ERROR);
                }
                if(!gomp_InputView)
                    gomp_PopAtomCoordinates();
                return(TCL_OK);
            }
/* .. OpenMol */
            else if(gomp_StringMatch(Text , "open$mol")) {

                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("structure index missing");
                    return(TCL_ERROR);
                }
                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("input file name is missing");
                    return(TCL_ERROR);
                }

                if(!gomp_InputView) {
                    if(gomp_PushAtomCoordinates())
                        return(TCL_ERROR);
                    if(gomp_DoViewingTransformationOverStructures( -1 )) {
                        gomp_PopAtomCoordinates();
                        return(TCL_ERROR);
                    }
                }

                if(!gomp_WriteCoordOPENMOL(atoi(Text2) - 1 , Text1 , 1)) {
                    if(!gomp_InputView)
                        gomp_PopAtomCoordinates();
                    return(TCL_OK);
                } else {
                    if(!gomp_InputView)
                        gomp_PopAtomCoordinates();
                    return(TCL_ERROR);
                }
            }
/* .. PROBESURF */
            else if(gomp_StringMatch(Text , "prob$esurf")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("structure index missing");
                    return(TCL_ERROR);
                }
                ITemp = atoi(Text2) - 1;

                if(!gomp_InputView) {
                    if(gomp_PushAtomCoordinates())
                        return(TCL_ERROR);
                    if(gomp_DoViewingTransformationOverStructures( -1 )) {
                        gomp_PopAtomCoordinates();
                        return(TCL_ERROR);
                    }
                }
            
            
                if(strlen(Text1) == 0) {
                    if(!gomp_InputView)
                        gomp_PopAtomCoordinates();
                    gomp_PrintERROR("input file is missing");
                    return(TCL_ERROR);
                }

                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                sumxyz    = gomp_GetTranslateArray();
                /* Xmin .. */
                if(strlen(Text2) != 0) {
                    Xmin = atof(Text2);
                }
                else 
                    Xmin = -(gomp_GetSizeOfSystem() + sumxyz[0]);

                /* Xmax .. */
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) != 0) {
                    Xmax = atof(Text2);
                }
                else
                    Xmax =  gomp_GetSizeOfSystem() + sumxyz[0];

                /* Ymin .. */
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) != 0) {
                    Ymin = atof(Text2);
                }
                else  
                    Ymin  = -(gomp_GetSizeOfSystem() + sumxyz[1]);

                /* Ymax .. */
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) != 0) {
                    Ymax = atof(Text2);
                }
                else
                    Ymax =   gomp_GetSizeOfSystem() + sumxyz[1];

                /* Zmin .. */
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) != 0) {
                    Zmin = atof(Text2);
                }
                else 
                    Zmin =  -(gomp_GetSizeOfSystem() + sumxyz[2]);

                /* Zmax .. */
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) != 0) {
                    Zmax = atof(Text2);
                }
                else
                    Zmax =   gomp_GetSizeOfSystem() + sumxyz[2];

                /* Xpts .. */
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) != 0) {
                    Xpts = atoi(Text2);
                }
                else 
                    Xpts = 60;
                /* Ypts .. */
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) != 0) {
                    Ypts = atoi(Text2);
                }
                else 
                    Ypts = 60;
                /* Zpts .. */
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) != 0) {
                    Zpts = atoi(Text2);
                }
                else 
                    Zpts = 60;


                /* ProbeRad .. */
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) != 0) {
                    ProbeRad = atof(Text2);
                }
                else 
                    ProbeRad = 2.0;

                if(!gomp_ExternalInput4PROBESURF(
                       ITemp , 
                       Text1 , Xmin , Xmax ,
                       Ymin , Ymax ,
                       Zmin , Zmax ,
                       Xpts , Ypts , Zpts ,
                       ProbeRad)) {
                    if(!gomp_InputView)
                        gomp_PopAtomCoordinates();
                    return(TCL_OK);
                } else {
                    if(!gomp_InputView)
                        gomp_PopAtomCoordinates();
                    return(TCL_ERROR);
                }
            }
        }
        else {
            gomp_PrintERROR("unrecognized 'export' command");
            return(TCL_ERROR);
        }
    }
/* .................. */ 

/*  E R R O R command not recognized         */
    gomp_PrintERROR("problems with the 'export' command");

    return(TCL_ERROR);
}

/*********************************************************************/
int WriteCoordinates(int Which ,
                     const char *Text1 ,
                     const char *Text2,
                     const char *Text3)
/*********************************************************************/
{
    static char Text[3*BUFF_LEN];
    int type;
    int ITemp;

    if(gomp_StringMatch(Text3 , "disp$laymask"))
        type = 1;
    else
        type = 0;


/*  coordinate file command */
/* ball and stick  */
    if(gomp_StringMatch(Text1 , "bands$tick")) {

      if(strlen(Text2) == 0) {
            gomp_PrintERROR("output file name is missing");
            return(TCL_ERROR);
        }

        if(gomp_PushAtomCoordinates())
            return(TCL_ERROR);

        if(gomp_DoViewingTransformationOverStructures( -1))  {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        }

        if(gomp_WriteCoordBaS(Which , Text2, type)) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        } else {
            gomp_PopAtomCoordinates();
            return(TCL_OK);
        }
    }
/* charmm */
    else if(gomp_StringMatch(Text1 , "char$mm") ||
            gomp_StringMatch(Text1 , "karp$lus")) {

      if(strlen(Text2) == 0) {
            gomp_PrintERROR("output file name is missing");
            return(TCL_ERROR);
        }
        if(gomp_PushAtomCoordinates())
            return(TCL_ERROR);

        if(gomp_DoViewingTransformationOverStructures( -1)) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        }

        if(gomp_WriteCoordCHARMM(Which , Text2, type)) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        } else { 
            gomp_PopAtomCoordinates();
            return(TCL_OK);
        }
    }
/* free  */
    else if(gomp_StringMatch(Text1 , "free")) {

      if(strlen(Text2) == 0) {
            gomp_PrintERROR("output file name is missing");
            return(TCL_ERROR);
        }

        if(gomp_PushAtomCoordinates())
            return(TCL_ERROR);

        if(gomp_DoViewingTransformationOverStructures( -1)) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        }

        if(gomp_WriteCoordFree(Which , Text2, type)) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        } else {
            gomp_PopAtomCoordinates();
            return(TCL_OK);
        }
    }
/* pdbq  */
    else if(gomp_StringMatch(Text1 , "pdbq")) {

      if(strlen(Text2) == 0) {
            gomp_PrintERROR("output file name is missing");
            return(TCL_ERROR);
        }

        if(gomp_PushAtomCoordinates())
            return(TCL_ERROR);

        if(gomp_DoViewingTransformationOverStructures( - 1)) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        }

        sprintf(Text,"lulAutoDock::MakePDBQCoordinates %d {%s} %d",Which,Text2,type);

        ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text);
        if(ITemp != TCL_OK) {
            gomp_PopAtomCoordinates();
            gomp_PrintERROR("can't execute the make PDBQ coordinates script 'lulMakeTinkerCoordinates'");
            return(TCL_ERROR);
        }
        else {
            gomp_PopAtomCoordinates();
            return(TCL_OK);
        }
    }
/* pdb  */
    else if(gomp_StringMatch(Text1 , "pdb")) {

      if(strlen(Text2) == 0) {
            gomp_PrintERROR("output file name is missing");
            return(TCL_ERROR);
        }

        if(gomp_PushAtomCoordinates())
            return(TCL_ERROR);

        if(gomp_DoViewingTransformationOverStructures( - 1)) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        }

        if(gomp_WriteCoordPDB(Which , Text2, type)) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        } else { 
            gomp_PopAtomCoordinates();
            return(TCL_OK);
        }
    }
/* openmol  */
    else if(gomp_StringMatch(Text1 , "open$mol")) {

      if(strlen(Text2) == 0) {
            gomp_PrintERROR("output file name is missing");
            return(TCL_ERROR);
        }

        if(gomp_PushAtomCoordinates())
            return(TCL_ERROR);

        if(gomp_DoViewingTransformationOverStructures( - 1)) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        }

        if(gomp_WriteCoordOPENMOL(Which , Text2, type)) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        } else {
            gomp_PopAtomCoordinates();
            return(TCL_OK);
        }
    }
/* xyz  */
    else if(gomp_StringMatch(Text1 , "xyz")) {

      if(strlen(Text2) == 0) {
            gomp_PrintERROR("output file name is missing");
            return(TCL_ERROR);
        }

        if(gomp_PushAtomCoordinates())
            return(TCL_ERROR);

        if(gomp_DoViewingTransformationOverStructures( - 1)) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        }

        if(gomp_WriteCoordXYZ(Which , Text2, type)) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        } else { 
            gomp_PopAtomCoordinates();
            return(TCL_OK);
        }
    }
    else {
        Tcl_Obj *list;

        if ( gomp_PushAtomCoordinates() )
            return(TCL_ERROR);

        if ( gomp_DoViewingTransformationOverStructures(-1) ) {
            gomp_PopAtomCoordinates();
            return(TCL_ERROR);
        }

        list = gomp_CreateTclList(
            "%s %d %s %s %d",
            "gom::WriteAtomCoordinates",
            Which, Text1, Text2, type);
        ITemp = Tcl_GlobalEvalObj(gomp_GetTclInterp(),list);
        if( ITemp == TCL_OK ) {
            if( atoi( Tcl_GetStringResult( gomp_GetTclInterp() ) ) != 0 )
                ITemp = TCL_ERROR;
        }

        gomp_PopAtomCoordinates();

        if(ITemp != TCL_OK)
            return(TCL_ERROR);
        return(TCL_OK);
    }
}


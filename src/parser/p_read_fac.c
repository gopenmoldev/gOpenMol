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

#include "cluster.h"
#include "coord_file.h"
#include "gaussian.h"
#include "gomfile.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "model_file.h"
#include "molecstruct.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "rforce.h"
#include "selection.h"
#include "tclutils.h"
#include "trajectory.h"

#include "stdafx.h"

#define GET_COMMAND_LENGTH 5

#if 0
static int GetTypeOfCoordinateFile(void);
static int TypeOfInputCoordinateFile;
#endif
static int ReadCoordinatesUsingBuildinFilter(const char * ,
                                             const char * , const char *);
static int ReadCoordinatesFile(const char * , int , int);

/*********************************************************************/
int gomp_ImportCommand(ClientData clientdata, Tcl_Interp *interp, 
                     int argc, const char **argv)
/*********************************************************************/
{
    static char Text[4*BUFF_LEN];
    static char Text1[BUFF_LEN];
    static char Text2[BUFF_LEN];
    static char Text3[BUFF_LEN];
    static int  Level;
    static int  ITemp;

/* #1   import ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "impo$rt")) {

/* #1.1 coor$dinates */
        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

        if(gomp_StringMatch(Text , "coor$dinates")) {

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_FileNameIsURL(Text2)) {
                sprintf(Text1,"?Can't handle the given URL '%s'",Text2);
                gomp_PrintERROR(Text1);
                return(TCL_ERROR);
            }

            if( gomp_ReadCoordinates(Text1,Text2,Text3) == 0 )
                return(TCL_OK);
            return(TCL_ERROR);
        }
/* postprocessing */
        else  if(gomp_StringMatch(Text , "post$process")) {
            if(gomp_PostReadAtoms(1,gomp_GetNumMolecStructs()-1)) {
                return(TCL_ERROR);
            }

#ifdef ENABLE_GRAPHICS
            if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {/* check if there is a display ... */
                (void)gomp_PrepareDisplay(gomp_GetNumMolecStructs());
            } 
#endif /* ENABLE_GRAPHICS */

            (void)gomp_SetFileSavingState(1);

            return(TCL_OK);
        }
/* cluster */
        else  if(gomp_StringMatch(Text , "clus$ter")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
            if(Text1[0] == (char)NULL) {
                gomp_PrintERROR("cluster file name missing");
                sprintf(Text3,"1");
                (void)gomp_SendTclReturn(Text3);
                return(TCL_ERROR);
            }

            if(gomp_FileNameIsURL(Text1)) {
                sprintf(Text1,"?Can't handle the given URL '%s'",Text1);
                gomp_PrintERROR(Text1);
                return(TCL_ERROR);
            }

            ITemp = gomp_ReadClusterData(Text1); 
            sprintf(Text2,"%d",ITemp);
            gomp_SendTclReturn(Text2);
            if(!ITemp)
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else  if(gomp_StringMatch(Text , "dict$ionary")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
            if(Text1[0] == (char)NULL) {
                gomp_PrintERROR("file name is missing");
                (void)gomp_SendTclReturn("1");
                return(TCL_ERROR);
            }
            if(gomp_FileNameIsURL(Text1)) {
                sprintf(Text1,"?Can't handle the given URL '%s'",Text1);
                gomp_PrintERROR(Text1);
                return(TCL_ERROR);
            }
            ITemp = gomp_ImportDictionaryAndApply(Text1);
            sprintf(Text2,"%d",ITemp);
            (void)gomp_SendTclReturn(Text2);
            if(!ITemp)
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/*  vector */
        else  if(gomp_StringMatch(Text , "vect$or")) {
            strncpy(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
            if(Text2[0] == (char)NULL) {
                gomp_PrintERROR("file type is missing");
                (void)gomp_SendTclReturn("1");
                return(TCL_ERROR);
            }

            ITemp = 0;

            if(gomp_StringMatch(Text2 , "charm$m")) {
                ITemp = 1;
            } else if(gomp_StringMatch(Text2 , "flat$file")) {
                ITemp = 2;
            } else {
                gomp_PrintERROR("unknow vector file supplied");
                (void)gomp_SendTclReturn("1");
                return(TCL_ERROR);
            }

            strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);

            if(gomp_FileNameIsURL(Text1)) {
                sprintf(Text1,"?Can't handle the given URL '%s'",Text1);
                gomp_PrintERROR(Text1);
                return(TCL_ERROR);
            }

            if(ITemp == 1) {
                ITemp = gomp_ReadCharmmVector(Text1);
            } else if (ITemp == 2) {
                ITemp = gomp_ReadFlatFileVector(Text1);
            } else {
                gomp_PrintERROR("unknow vector file supplied");
                (void)gomp_SendTclReturn("1");
                return(TCL_ERROR);
            }
            sprintf(Text3,"%d",ITemp);
            (void)gomp_SendTclReturn(Text3);
            if(!ITemp)
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* #1.3 gbas$is */
        else if(gomp_StringMatch(Text , "gbas$is")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);

            if(Text1[0] == (char)NULL) {
                gomp_PrintERROR("file name is missing");
                (void)gomp_SendTclReturn("1");
                return(TCL_ERROR);
            }

            if(gomp_FileNameIsURL(Text1)) {
                sprintf(Text1,"?Can't handle the given URL '%s'",Text1);
                gomp_PrintERROR(Text1);
                return(TCL_ERROR);
            }

            ITemp = gomp_ReadGaussianBasisLineInput(Text1);
            sprintf(Text3,"%d",ITemp);
            (void)gomp_SendTclReturn(Text3);
            if(!ITemp)
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* model */
        else  if(gomp_StringMatch(Text , "mode$l")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);

            if(Text1[0] == (char)NULL) {
                gomp_PrintERROR("file name is missing");
                (void)gomp_SendTclReturn("1");
                return(TCL_ERROR);
            }

            if(gomp_FileNameIsURL(Text1)) {
                sprintf(Text1,"?Can't handle the given URL '%s'",Text1);
                gomp_PrintERROR(Text1);
                return(TCL_ERROR);
            }
            Level = gomp_ReadOldModel(Text1);

            sprintf(Text2,"%d",Level);
            (void)gomp_SendTclReturn(Text2);
            if(!Level)
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* charge */
        else  if(gomp_StringMatch(Text , "char$ge")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
            if(Text1[0] == (char)NULL) {
                gomp_PrintERROR("source type missing");
                return(TCL_ERROR);
            }
/* ICON8 */
            if(gomp_StringMatch(Text1 , "icon8") ||
               gomp_StringMatch(Text1 , "ICON8")) {
                strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                        BUFF_LEN-1);
                if(Text1[0] == (char)NULL) {
                    gomp_PrintERROR("ICON8 output file name missing");
                    return(TCL_ERROR);
                }

                if(gomp_FileNameIsURL(Text1)) {
                    sprintf(Text1,"?Can't handle the given URL '%s'",Text1);
                    gomp_PrintERROR(Text1);
                    return(TCL_ERROR);
                }

                sprintf(Text,"lulImportChargeFromICON8output {%s}",
                        Text1 );
                ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text);
                if(ITemp != TCL_OK) {
                    gomp_PrintERROR("can't execute the read ICON8 charge script 'lulImportChargeFromICON8output'");
                    return(TCL_ERROR);
                }

                return(TCL_OK);
            }
/* MOPAC6 */
            else if(gomp_StringMatch(Text1 , "mopa$c6") ||
                    gomp_StringMatch(Text1 , "MOPA$C6")) {
                strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                        BUFF_LEN-1);
                if(Text1[0] == (char)NULL) {
                    gomp_PrintERROR("MOPAC6 output file name missing");
                    return(TCL_ERROR);
                }
                if(gomp_FileNameIsURL(Text1)) {
                    sprintf(Text1,"?Can't handle the given URL '%s'",Text1);
                    gomp_PrintERROR(Text1);
                    return(TCL_ERROR);
                }
                sprintf(Text,"lulImportChargeFromMOPAC6output {%s}",
                        Text1 );
                ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text);
                if(ITemp != TCL_OK) {
                    gomp_PrintERROR("can't execute the read MOPAC6 charge script 'lulImportChargeFromMOPAC6output'");
                    return(TCL_ERROR);
                }

                return(TCL_OK);
            }
/* GAUSSIAN */
            else if(gomp_StringMatch(Text1 , "gaus$sian") ||
                    gomp_StringMatch(Text1 , "GAUS$SIAN")) {
                strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                        BUFF_LEN-1);
                if(Text1[0] == (char)NULL) {
                    gomp_PrintERROR("GAUSSIAN output file name missing");
                    return(TCL_ERROR);
                }
                if(gomp_FileNameIsURL(Text1)) {
                    sprintf(Text1,"?Can't handle the given URL '%s'",Text1);
                    gomp_PrintERROR(Text1);
                    return(TCL_ERROR);
                }
                sprintf(Text,"lulImportChargeFromGAUSSIANoutput {%s}",
                        Text1 );
                ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text);
                if(ITemp != TCL_OK) {
                    gomp_PrintERROR("can't execute the read GAUSSIAN charge script 'lulImportChargeFromGAUSSIANoutput'");
                    return(TCL_ERROR);
                }

                return(TCL_OK);
            }
/* USER */
            else if(gomp_StringMatch(Text1 , "user") ||
                    gomp_StringMatch(Text1 , "MOPA$C6")) {
                strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                        BUFF_LEN-1);
                if(Text1[0] == (char)NULL) {
                    gomp_PrintERROR("USER output file name missing");
                    return(TCL_ERROR);
                }
                if(gomp_FileNameIsURL(Text1)) {
                    sprintf(Text1,"?Can't handle the given URL '%s'",Text1);
                    gomp_PrintERROR(Text1);
                    return(TCL_ERROR);
                }
                sprintf(Text,"lulImportChargeFromUSERoutput {%s}",
                        Text1 );
                ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text);
                if(ITemp != TCL_OK) {
                    gomp_PrintERROR("can't execute the read USER charge script 'lulImportChargeFromUSERoutput'");
                    return(TCL_ERROR);
                }

                return(TCL_OK);
            }
            else {
                gomp_PrintERROR("'import charge' command not recognized");
                return(TCL_ERROR);
            }
        }
    }
/* .................. */ 

/*  E R R O R command not recognized         */
    gomp_PrintERROR("problems with the 'import' command");
    sprintf(Text3,"problems with the 'import' command");
    (void)gomp_SendTclReturn(Text3);
    return(TCL_ERROR);
}

/*********************************************************************/
int gomp_FileNameIsURL(char *Text)
/*********************************************************************/
{
    int   ITemp;
    char  NextText[BUFF_LEN];
    const char *value;

/* 
   check to see if this is a http link if yes download the file into the
   gOpenMol temp directory and read it in from there 
*/
    if(Tcl_StringCaseMatch(Text,"http://*",1)) {

        sprintf(NextText,"lulUtility::gomGetCoordinatesURL {%s}",Text);
        ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), NextText);
        if(ITemp != TCL_OK) {
            gomp_PrintERROR("can't execute the read download URL coordinates script 'lulUtility::gomGetCoordinatesURL'");
            return(1);
        }
          
        value  = Tcl_GetVar(gomp_GetTclInterp() , "gomURLFileName", TCL_GLOBAL_ONLY);

        if(!value || (value[0] == (char)NULL)) {
            gomp_PrintERROR("can't retrieve given URL");
            return(1);
        } 
          
        gomp_CopyString(Text,value,BUFF_LEN);
    }

    return(0);

}

/*********************************************************************/
int gomp_ReadCoordinates(
    const char *Type, const char *FileName, const char *Action)
/*********************************************************************/
{
    static int  NStruct,i;
    static int  Level;
    static int  ITemp;
    static int  FileType;
    Tcl_Obj    *list;

    NStruct = gomp_GetNumMolecStructs();

    if( gomp_StringMatch( Action , "appe$nd" ) )
        Action = "append";
    else {
        NStruct = 0;
        Action  = "new";
    }

    /* call Tcl pre read procedure */
    list = gomp_CreateTclList(
        "%s %s %s %s",
        "gom::PreReadAtomCoordinates",
        Type, FileName, Action);
    ITemp = Tcl_GlobalEvalObj(gomp_GetTclInterp(), list);
    if (ITemp != TCL_OK)
        return(1);

    Level = ReadCoordinatesUsingBuildinFilter(Type, FileName, Action);

    if ( Level < 0 ) {
        /* Not a build in type. */
        /* call Tcl read procedure */
        list = gomp_CreateTclList(
            "%s %s %s %s",
            "gom::ReadAtomCoordinates",
            Type, FileName, Action);
        ITemp = Tcl_GlobalEvalObj(gomp_GetTclInterp(), list);
        if(ITemp != TCL_OK)
            return(1);
        
        ITemp = atoi( Tcl_GetStringResult( gomp_GetTclInterp() ) );
    }
    else if ( Level > 0 )
        ITemp = 0;
    else
        ITemp = 1;

    if ( Level != 0 ) {
        /* select a structure */
        (void)gomp_ActivateSelectedStructure(
            gomp_GetNumMolecStructs() - 1 , STRUCTURE_SELECTION_ON);

        if ( Level != DYNAFRAME_COORD ) {
            for( i = NStruct ; i < gomp_GetNumMolecStructs() ; i++ ) {
                if ( gomp_PostReadAtoms( 1 , i ) ) {
                    gomp_PrintERROR("PROBLEMS in 'gomp_PostReadAtoms'\n");
                    return(1);
                }
            }
        
#ifdef ENABLE_GRAPHICS
            if ( gomp_GetTermType() == GRAPHICS_AVAILABLE ) {
                /* check if there is a display ... */
                (void)gomp_PrepareDisplay(gomp_GetNumMolecStructs());
                (void)gomp_PrepareStatusDisplay(FileName);
            }
#endif /* ENABLE_GRAPHICS */
        }
    }

    /* call Tcl post read procedure */
    list = gomp_CreateTclList(
        "%s %s %s %s %d",
        "gom::PostReadAtomCoordinates",
        Type, FileName, Action, ITemp);
    if( Tcl_GlobalEvalObj(gomp_GetTclInterp(), list) != TCL_OK)
        return(1);

/* if the input file is an XMOL file read it also into the trajectory input */
    if( ( Level == XMOL_COORD ) && (gomp_GetNumMolecStructs() == 1)) {
        FileType = gomp_ParseTrajType("xmol"); 
        (void)gomp_GetTrajectory(FileName, "\0", FileType);
/* if only one frame drop the trajectory */
        if(gomp_GetNumberOfFrames() < 2) {
            (void)gomp_SetTrajectoryStructureOff();
        }
    }

    gomp_SendTclObjReturn(Tcl_NewIntObj(ITemp));
    
    return(0);
}


/*********************************************************************/
int ReadCoordinatesUsingBuildinFilter(const char *Text1 ,
                                      const char *Text2 ,
                                      const char *Text3)
/*********************************************************************/
{
    static int   append;

    append = 0;

/*  coordinate file command */
/* amber */
    if(gomp_StringMatch(Text1 , "ambe$r")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2, AMBER_COORD , append))
            return(0);
        else 
            return(AMBER_COORD);
    }
/* charmm */
    else if(gomp_StringMatch(Text1 , "char$mm") ||
            gomp_StringMatch(Text1 , "karp$lus")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2, CHARMM_COORD , append))
            return(0);
        else 
            return(CHARMM_COORD);
    }
/*  get frame from file command */
    else if(gomp_StringMatch(Text1 , "fram$e")) {
        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2, DYNAFRAME_COORD , append))
            return(0);
        else 
            return(DYNAFRAME_COORD);
    }
/* gaussian */
    else if(gomp_StringMatch(Text1 , "gaus$sian")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2, GAUSSIAN_COORD , append))
            return(0);
        else 
            return(GAUSSIAN_COORD);
    }
/* hyper */
    else if(gomp_StringMatch(Text1 , "hype$rchem")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2,HYPER_COORD , append))
            return(0);
        else 
            return(HYPER_COORD);
    }
/* insight */
    else if(gomp_StringMatch(Text1 , "insi$ght")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2, INSIGHT_COORD , append))
            return(0);
        else 
            return(INSIGHT_COORD);
    }
/* mol2 */
    else if(gomp_StringMatch(Text1 , "mol2")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2, MOL2_COORD , append))
            return(0);
        else 
            return(MOL2_COORD);
    }
/* mopac graphics file */
    else if(gomp_StringMatch(Text1 , "mopa$c")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2, MOPAC_COORD , append))
            return(0);
        else 
            return(MOPAC_COORD);
    }
/* mumod */
    else if(gomp_StringMatch(Text1 , "mumo$d")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2, MUMOD_COORD , append))
            return(0);
        else 
            return(MUMOD_COORD);
    }
/* openmol */
    else if(gomp_StringMatch(Text1 , "open$mol")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2, OPENMOL_COORD , append))
            return(0);
        else 
            return(OPENMOL_COORD);
    }
/* pdb */
    else if(gomp_StringMatch(Text1 , "pdb")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2,PDB_COORD , append))
            return(0);
        else 
            return(PDB_COORD);
    }
/* XMOL */
    else if(gomp_StringMatch(Text1 , "xmol") || gomp_StringMatch(Text1 , "xyz")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2, XMOL_COORD , append))
            return(0);
        else 
            return(XMOL_COORD);
    }
/* GXYZ */
    else if(gomp_StringMatch(Text1 , "gxyz")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2, GXYZ_COORD , append))
            return(0);
        else 
            return(GXYZ_COORD);
    }
/* yasp */
    else if(gomp_StringMatch(Text1 , "yasp")) {

        append = gomp_StringMatch(Text3 , "appe$nd");

        if(ReadCoordinatesFile(Text2, YASP_COORD , append))
            return(0);
        else 
            return(YASP_COORD);
    }
    else
        return -1;
}

/*********************************************************************/
int ReadCoordinatesFile(const char *FileName , int Alt , int Append)
/*********************************************************************/
{
    if((Alt != OPENMOL_COORD) &&
       (Alt != DYNAFRAME_COORD)) {
        if(gomp_Check_if_file_exists(FileName)) {
            gomp_PrintERROR(" File does not exist");
            gomp_PrintMessage(FileName);
            return(1);
        }
    }
#if 0
    TypeOfInputCoordinateFile = Alt;
#endif
/*  coordinates command */
    switch(Alt) {

    case AMBER_COORD:     /* AMBER coordinates */

        if(gomp_ReadCoordinatesAMBER( FileName , Append))
            return(1);
        break;

    case CHARMM_COORD: /* CHARMM coordinates */

        if(gomp_ReadCoordinatesCHARMM( FileName , Append))
            return(1);
        break;

    case DYNAFRAME_COORD:     /* Dynamics frame coordinates */

        if(gomp_ReadCoordinatesFRAME( atoi(FileName) , Append))
            return(1);
        break;

    case GAUSSIAN_COORD:     /* Gaussian coordinates */

        if(gomp_ReadCoordinatesGAUSSIAN( FileName , Append))
            return(1);
        break;

    case HYPER_COORD:     /* HyperChem coordinates */

        if(gomp_ReadCoordinatesHYPERCHEM( FileName , Append))
            return(1);
        break;

    case INSIGHT_COORD:     /* Insight coordinates */

        if(gomp_ReadCoordinatesINSIGHT( FileName , Append))
            return(1);
        break;

    case MOL2_COORD:     /* Mol2 coordinates */

        if(gomp_ReadCoordinatesMOL2( FileName , Append))
            return(1);
        break;

    case MOPAC_COORD:     /* Mopac graphics coordinates */

        if(gomp_ReadCoordinatesMOPACgraph( FileName , Append))
            return(1);
        break;

    case MUMOD_COORD:     /* Mumod coordinates */

        if(gomp_ReadCoordinatesMUMOD( FileName , Append))
            return(1);
        break;

    case OPENMOL_COORD: /* OPENMOL coordinates */

        if(gomp_ReadCoordinatesOPENMOL( FileName , Append))
            return(1);
        break;

    case PDB_COORD:     /* PDB coordinates */

        if(gomp_ReadCoordinatesPDB( FileName , Append))
            return(1);
        break;

    case XMOL_COORD:   /* XMOL coordinates */

        if(gomp_ReadCoordinatesXMOL( FileName , Append))
            return(1);
        break;

    case GXYZ_COORD:     /* GXYZ coordinates */

        if(gomp_ReadCoordinatesGXYZ( FileName , Append))
            return(1);
        break;

    case YASP_COORD:     /* YASP coordinates */

        if(gomp_ReadCoordinatesYASP( FileName , Append))
            return(1);
        break;

    default:
        gomp_PrintERROR("Unknown coordinate file type");
        return(1);
        break;
    }

    return(0);
}
#if 0
/*********************************************************************/
int GetTypeOfCoordinateFile()
/*********************************************************************/
{
    return(TypeOfInputCoordinateFile);
}
#endif
/***************************************************************************/

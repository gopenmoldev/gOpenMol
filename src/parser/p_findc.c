/*

Copyright (c) 1996 - 2005 by:
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
#include <math.h>
#include <tcl.h>
#include <sys/types.h>

#include "bond.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "listutils.h"
#include "molecstruct.h"
#include "plot.h"
#include "plumber.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_FindCommand(ClientData clientdata, Tcl_Interp *interp,
                   int argc, const char **argv)
/*********************************************************************/
{
    static int  ITemp,ITemp2;
    static char Text[BUFF_LEN];
    static char Text1[BUFF_LEN];
    static char Text2[BUFF_LEN];
    static char Text3[BUFF_LEN];
    static char Text4[BUFF_LEN];
    static char Text5[BUFF_LEN];
    static char Text6[BUFF_LEN];
    static char Text7[BUFF_LEN];
    static char Text8[BUFF_LEN];

/* #1   find  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "find")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
/* S-S bonds */
        if(gomp_StringMatch(Text , "ssbo$nds")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_StringMatch(Text1 , "all")) {
                if(!gomp_FindSSbondsAll(atof(Text)))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            } else {
                if(!gomp_FindSSbondsI((atoi(Text1) - 1) , atof(Text)))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
        }
/* H bonds */
        else if(gomp_StringMatch(Text , "hbon$ds")) {
            if(!gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("no molecular structure(s) defined");
                return(TCL_ERROR);
            }

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_StringMatch(Text4 , "-hyd$rogens")) {
                if(gomp_ParseCalcHbondsList(Text1 , Text2 , Text3 , 0 , ""))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
            else {
                if(gomp_ParseCalcHbondsList(Text1 , Text2 , Text3 , 1 , ""))
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
        else if(gomp_StringMatch(Text , "chai$ns")) {
            if(!gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("no molecular structure(s) defined");
                return(TCL_ERROR);
            }

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(*Text3 == (char)NULL)
                ITemp = 0;
            else if(gomp_StringMatch(Text3 , "all"))
                ITemp = LIST_ALL_STRUCTURES;
            else if(gomp_StringMatch(Text3 , "sele$cted"))
                ITemp = LIST_SELECTED_STRUCTURES;
            else {
                ITemp = atoi(Text3);

                if(ITemp < 1 || ITemp > gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("structure number out of range");
                    return(TCL_ERROR);
                }

                --ITemp;
            }

            if(gomp_ParseListProteinChains(ITemp , Text1 , Text2))
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "atom$s") ||
                gomp_StringMatch(Text , "resi$due") ||
                gomp_StringMatch(Text , "segm$ent") ||
                gomp_StringMatch(Text , "sel$list")) {
            if(!gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("no molecular structure(s) defined");
                return(TCL_ERROR);
            }

            if( gomp_StringMatch(Text , "atom$s") )
                ITemp2 = ATOM_LIST;
            else if( gomp_StringMatch(Text , "resi$due") )
                ITemp2 = RESIDUE_LIST;
            else if( gomp_StringMatch(Text , "segm$ent") )
                ITemp2 = SEGMENT_LIST;
            else
                ITemp2 = SELECTION_LIST;

            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_StringMatch(Text , "around")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text7,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

                gomp_CopyString(Text8,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

                if(*Text8 == (char)NULL)
                    ITemp = 0;
                else if(gomp_StringMatch(Text8 , "all"))
                    ITemp = LIST_ALL_STRUCTURES;
                else if(gomp_StringMatch(Text8 , "sele$cted"))
                    ITemp = LIST_SELECTED_STRUCTURES;
                else {
                    ITemp = atoi(Text8);

                    if(ITemp < 1 || ITemp > gomp_GetNumMolecStructs()) {
                        gomp_PrintERROR("structure number out of range");
                        return(TCL_ERROR);
                    }

                    --ITemp;
                }

                if(gomp_ParseListAtomsAround( ITemp , Text1 , Text2 , Text3 , Text4 ,
                                            Text5 , Text6 , Text7 , ITemp2 ))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
            else if(gomp_StringMatch(Text , "join")) {
                if(*argv[argc-1] == (char)NULL)
                    ITemp = 0;
                else if(gomp_StringMatch(argv[argc-1] , "all"))
                    ITemp = LIST_ALL_STRUCTURES;
                else if(gomp_StringMatch(argv[argc-1] , "sele$cted"))
                    ITemp = LIST_SELECTED_STRUCTURES;
                else {
                    ITemp = atoi(argv[argc-1]);

                    if(ITemp < 1 || ITemp > gomp_GetNumMolecStructs()) {
                        gomp_PrintERROR("structure number out of range");
                        return(TCL_ERROR);
                    }

                    --ITemp;
                }

                if( (argc - 4) % 3 ) {
                    gomp_PrintERROR("'find *** join' command parametre count mismatch");
                    return(TCL_ERROR);
                }

                if(gomp_ParseListAtomsJoin( ITemp , (argc - 4)/3 , argv + 3 , ITemp2 ))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
        }
        else {
            gomp_PrintERROR("'find' command not recognized");
            return(TCL_ERROR);
        }

    }
/*  E R R O R command not recognized         */
    gomp_PrintERROR("'find' command not recognized");
    return(TCL_ERROR);

}


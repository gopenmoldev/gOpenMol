/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2001 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <sys/types.h>

#include "bond.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "plot.h"
#include "tclutils.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_EditCommand(ClientData clientdata, Tcl_Interp *interp,
                   int argc, const char **argv)
/*********************************************************************/
{
    static int    bond_type;
    static char   Text[BUFF_LEN];
    static char   Text1[BUFF_LEN];
    static char   Text2[BUFF_LEN];
    static char   Text3[BUFF_LEN];
    static char   Text4[BUFF_LEN];
    static char   Text5[BUFF_LEN];
    static char   Text6[BUFF_LEN];

/* #1   edit  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "edit")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

/* bond */
        if(gomp_StringMatch(Text , "bond") ||
           gomp_StringMatch(Text , "hbon$d")) {
            if(!gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("no atom structures defined");
                return(TCL_ERROR);
            }

            if(gomp_StringMatch(Text , "bond"))
                bond_type = 0;
            else
                bond_type = 1;

            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(gomp_StringMatch(Text , "brea$k")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text4[0] != (char)NULL) {
                    if(!gomp_EditBondM(bond_type , 2 , Text1, Text2, Text3,
                                     Text4, Text5, Text6))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else {
                    if(!gomp_BreakBondM(bond_type , Text1, Text2, Text3))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
            }
            else if(gomp_StringMatch(Text , "crea$te")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(Text4[0] != (char)NULL) { 
                    if(!gomp_EditBondM(bond_type , 1 , Text1, Text2, Text3, 
                                     Text4, Text5, Text6))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
            }
        }
/* hbsubset */
        else if(gomp_StringMatch(Text , "hbsu$bset")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_ParseCalcHbondSubset(Text1, Text2, Text3))
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "-hbs$ubset")) {
            if(gomp_DeleteAllHbondSubsets())
                return(TCL_ERROR);
            else
                return(TCL_OK);
        }
    }


/*  E R R O R command not recognized         */
    gomp_PrintERROR("'edit' command not recognized");
    return(TCL_ERROR);

}


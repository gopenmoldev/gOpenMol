/*

Copyright (c) 1994 - 2004 by:
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

#include "gommonitor.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

#define DIST_TYPE  0
#define ANG_TYPE   1
#define TORS_TYPE  2

static const char *ManiList[] = {
    "dave$rage"  ,
    "squa$re"    ,
    "cos"        ,
    "cos2"       ,
    "sqrt"       ,
    "dini$tial"  ,
    "copy"       ,
    "add"        ,
    "log"        ,
    "exp"        ,
    "powe$rreal" ,
    "mult$real"  ,
    "divi$dereal",
    "shif$treal" ,
    "dmin"       ,
    "abs"        ,
    "divf$irst"  ,
    "divm$aximum",
    "pspe$ctrum" ,
    "zero"};

static int  ManipulationAction(const char *);

/*********************************************************************/
int gomp_ManipulateCommand(ClientData clientdata, Tcl_Interp *interp,
                         int argc, const char **argv)
/*********************************************************************/
{
    static int    type;
    static char   Text[BUFF_LEN];
    static char   Text1[BUFF_LEN];
    static char   Text2[BUFF_LEN];
    static char   Text3[BUFF_LEN];
    static char   Text4[BUFF_LEN];
    static char   Text5[BUFF_LEN];

/* #1   manipulate  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);

    if(gomp_StringMatch(Text , "mani$pulate")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

/* time series */
        if(gomp_StringMatch(Text , "time$serie")) {

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            if(gomp_StringMatch(Text1 , "dist$ance")) {
                type = ManipulationAction(Text2);
                if(strlen(Text3) == 0) {
                    gomp_PrintERROR("index into list is missing");
                    return(TCL_ERROR);
                }
                if(type) {
                    if(!gomp_ManipulateTimeSeries(DIST_TYPE , type , 
                                                Text3 , Text4 , Text5))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else {
                    gomp_PrintERROR("unrecognized manipulation action");
                    return(TCL_ERROR);
                }
            }
            else if(gomp_StringMatch(Text1 , "angl$e")) {
                type = ManipulationAction(Text2);
                if(strlen(Text3) == 0) {
                    gomp_PrintERROR("index into list is missing");
                    return(TCL_ERROR);
                }
                if(type) {
                    if(!gomp_ManipulateTimeSeries(ANG_TYPE , type , 
                                                Text3 , Text4 , Text5))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else {
                    gomp_PrintERROR("unrecognized manipulation action");
                    return(TCL_ERROR);
                }
            }
            else if(gomp_StringMatch(Text1 , "tors$ion")) {
                type = ManipulationAction(Text2);
                if(strlen(Text3) == 0) {
                    gomp_PrintERROR("index into list is missing");
                    return(TCL_ERROR);
                }
                if(type) {
                    if(!gomp_ManipulateTimeSeries(TORS_TYPE , type , 
                                                Text3 , Text4 , Text5))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else {
                    gomp_PrintERROR("unrecognized manipulation action");
                    return(TCL_ERROR);
                }
            }
        }
        else {
            gomp_PrintERROR("command 'manipulate' not recognized");
            return(TCL_ERROR);
        }

    }
/*  e r r o r command not recognized         */
    gomp_PrintERROR("'manipulate' command not recognized");
    return(TCL_ERROR);

}

/*********************************************************************/
int  ManipulationAction(const char *action)
/*********************************************************************/
{
    static int i;
     
    for(i = 0 ; i < (int)(sizeof(ManiList)/sizeof(const char *)) ; i++) {

        if(gomp_StringMatch(action , ManiList[i])) {
            return(i+1);
        }
    }

    return(0);
}

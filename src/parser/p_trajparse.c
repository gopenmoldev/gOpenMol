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

#include "gomfile.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "printmsg.h"
#include "projview.h"
#include "tclutils.h"
#include "trajectory.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_TrajectoryCommand(ClientData clientdata, Tcl_Interp *interp,
                         int argc, const char **argv)
/*********************************************************************/
{
    static int   FileType;
    static char  Text[BUFF_LEN];
    static char  Text1[BUFF_LEN];
    static char  Text2[BUFF_LEN];
    static char  Text3[BUFF_LEN];
    static char  Text4[BUFF_LEN];

/* #1   trajectory ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "traj$ectory")) {

/* #1.1 file */

        gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

        if(gomp_StringMatch(Text1,"file"))
        {

            if(!(FileType = gomp_ParseTrajType(gomp_GetNextFromParserStack(argc,(const char **)NULL)))) 
                return(TCL_ERROR);

            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(*Text == (char)NULL) {
                gomp_PrintERROR("Incomplete 'trajectory file' command, file name missing");
                return(TCL_ERROR);
            }

            if(gomp_FileNameIsURL(Text)) {
                sprintf(Text,"?Can't handle the given URL '%s'",Text);
                gomp_PrintERROR(Text);
                return(TCL_ERROR);
            }

            (void)gomp_GetTrajectory(Text,
                                   "\0",
                                   FileType);

            return(TCL_OK);
        }
/* delete trajectory information */
        else if(gomp_StringMatch(Text1,"-fil$e")) {
                 (void)gomp_SetTrajectoryStructureOff();
        }

/* limits */
        else if(gomp_StringMatch(Text1,"limi$ts")) {
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text2[0] == (char)NULL) {
                gomp_PrintERROR("values missing");
                return(TCL_ERROR);
            }
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(!gomp_SetTrajectoryDisplayParams(
                   atoi(Text2),atoi(Text3),atoi(Text4)))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text1,"acti$on")) {
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text2[0] == (char)NULL) {
                gomp_PrintERROR("value missing");
                return(TCL_ERROR);
            }

            if(gomp_StringMatch(Text2,"slow")) {
                (void)gomp_SetFormattedTrajectoryReader(1);
                return(TCL_OK);
            } else if (gomp_StringMatch(Text2,"fast")) {
                (void)gomp_SetFormattedTrajectoryReader(0);
                return(TCL_OK);
            }

            (void)gomp_SetFormattedTrajectoryReader(atoi(Text2));

            return(TCL_OK);
        }
/* retrieve velocities/forces */
        else if(gomp_StringMatch(Text1,"retr$ieve")) {
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text2[0] == (char)NULL) {
                gomp_PrintERROR("retrieve either 'velocities of forces'");
                return(TCL_ERROR);
            }

            if(gomp_StringMatch(Text2,"velo$cities")) {
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(gomp_StringMatch(Text3 , "on")) {
                    if(!gomp_SetVelocityRetrieveState(ON)) {
                        return(TCL_OK);
                    }
                    else {
                        return(TCL_ERROR);
                    }
                }
                else if(gomp_StringMatch(Text3 , "off")) {
                    if(!gomp_SetVelocityRetrieveState(OFF)) {
                        return(TCL_OK);
                    }
                    else {
                        return(TCL_ERROR);
                    }
                }
                else {
                    gomp_PrintERROR("option has to be 'on' or 'off'");
                    return(TCL_ERROR);
                }
            } else if(gomp_StringMatch(Text2,"forc$es")) {
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
                if(gomp_StringMatch(Text3 , "on")) {
                    if(!gomp_SetForceRetrieveState(ON)) {
                        return(TCL_OK);
                    }
                    else {
                        return(TCL_ERROR);
                    }
                }
                else if(gomp_StringMatch(Text3 , "off")) {
                    if(!gomp_SetForceRetrieveState(OFF)) {
                        return(TCL_OK);
                    }
                    else {
                        return(TCL_ERROR);
                    }
                }
                else {
                    gomp_PrintERROR("option has to be 'on' or 'off'");
                    return(TCL_ERROR);
                }
            } else  {
                gomp_PrintERROR("retrieve either 'velocities or forces'");
                return(TCL_ERROR);
            }
        }
/* loop */
        else if(gomp_StringMatch(Text1,"loop")) {
#ifdef ENABLE_GRAPHICS
/* check that there is graphics */
            if(gomp_GetTermType() != GRAPHICS_AVAILABLE) return(TCL_OK);
/*     done ...                 */
            strncpy(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
            if(gomp_StringMatch(Text1,"play")) {

/* check that a trajectory file is defined ... */
                if(!gomp_GetTrajectoryStructureState()) {
                    gomp_PrintERROR("no trajectory information available");
                    return(TCL_ERROR);
                }
                (void)gomp_DisplayTrajectory(TRAJ_PLAY);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1,"stop")) {
                (void)gomp_DisplayTrajectory(TRAJ_STOP_LOOP);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1,"forw$ard")) {
                (void)gomp_DisplayTrajectory(TRAJ_FORWARD_FRAME);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1,"back$ward")) {
                (void)gomp_DisplayTrajectory(TRAJ_BACKWARD_FRAME);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1,"firs$t")) {
                (void)gomp_DisplayTrajectory(TRAJ_FIRST_FRAME);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1,"last")) {
                (void)gomp_DisplayTrajectory(TRAJ_LAST_FRAME);
                return(TCL_OK);
            }
            else {
                gomp_PrintERROR("uncomplete 'trajectory loop' command");
                gomp_SendTclReturn("uncomplete 'trajectory loop' command");
                return(TCL_ERROR);
            }
#else
            return(TCL_OK);
#endif /* ENABLE_GRAPHICS */
        } else {
            gomp_PrintERROR("unknown 'trajectory' command defined");
            return(TCL_ERROR);
        }
    }

/*  E R R O R command not recognized         */

    return(TCL_OK);
}


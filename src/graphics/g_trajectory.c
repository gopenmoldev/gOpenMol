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
#include <math.h>
#include <ctype.h>
#include <tcl.h>

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>

#if defined(GLUT)
#include <GL/glut.h>
#else
#include <GL/glaux.h>
#endif

#if !defined(WIN32)
#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <X11/cursorfont.h>
#include <X11/Intrinsic.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/glx.h>
#endif

#if defined(MOTIF)
#include <GL/GLwMDrawA.h>
#endif
#endif /* ENABLE_GRAPHICS */

#include "colouring.h"
#include "drawscene.h"
#include "gomfile.h"
#include "gommain.h"
#include "gomtext.h"
#include "gomfile.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"
#include "tclutils.h"

#include "stdafx.h"

/*
  #include "play_button.xbm"
  #include "fforward_button.xbm"
  #include "freverse_button.xbm"
  #include "stop_button.xbm"
  #include "pause_button.xbm"
  #include "last_button.xbm"
  #include "first_button.xbm"
*/

#define   MAIN_TRAJECTORY "main_trajectory_widget.hlp"

static struct {
    int    DisplayState;
    char   FontName[BUFF_LEN];
    float  xp; 
    float  yp;
    float  Red;
    float  Green;
    float  Blue;
} FrameNumberDisplay = {
    1 , "BITMAP_HELVETICA_12" , 0.75 , 0.9 , 1.0 , 0.0 , 0.0};

static int PostReadFrame(int);

/*********************************************************************/
int gomp_GetTrajectory(const char *FileName , const char *Append , int Alt)
/*********************************************************************/
{
    FILE *File_p;

    if(Alt != OPENMOL_TRAJ) {
        if(gomp_Check_if_file_exists(FileName)) {
            gomp_PrintERROR("File does not exist");
            gomp_PrintMessage(FileName);
            return(1);
        }
    }

    (void)gomp_SetTrajectoryFileType(Alt);

    (void)gomp_SetTrajectoryFileName(FileName);

    if(gomp_GetTrajectoryFileType() == XMOL_TRAJ      ||
       gomp_GetTrajectoryFileType() == GROMOS96A_TRAJ ||
       gomp_GetTrajectoryFileType() == TINKER_TRAJ    ||
       gomp_GetTrajectoryFileType() == FAMBER_TRAJ    ||
       gomp_GetTrajectoryFileType() == FDL_POLY_TRAJ) {

           if(gomp_CheckTextFileLineEnding(FileName)) {
             File_p = fopen(FileName , "rb");
             gomp_PrintMessage("\n********************************************************\nERROR - it looks as if your trajectory is not a proper text file!\nWill try to solve the problem!\n********************************************************\n");
           } else {
             File_p = fopen(FileName , "r");
           }
    } else {
        File_p = fopen(FileName , "rb");
    }

    if(File_p == NULL) return(1);

/* this reads now only the header information */
    if(gomp_GetOneFrame(0 , File_p , 0)) {
        fclose(File_p);
        return(1);
    }

    (void)gomp_SetTrajectoryStructureInUse();

/* get frame number one (first frame) */
    if(gomp_GetOneFrame(1 , File_p , 0)) {
        fclose(File_p);
        return(1);
    }

    fclose(File_p);
    return(0);
}

static struct {
    FILE *File;
    int   ID;
    int   LoopDisplayDelay;
    char  HasIdleProc;
    char  HasTimerProc;
} TrajectoryPlay = { NULL, 0, 0, 0, 0 };

static void PlayTrajectoryTimer(ClientData changeFrame);

static void PlayTrajectory(ClientData changeFrame)
{
    int UsePostRead;
    int FirstFrame, Frame, LastFrame, FrameStep;
    Tcl_Obj *var;
    int MyID = TrajectoryPlay.ID;

    TrajectoryPlay.HasIdleProc = 0;

    if ( ! TrajectoryPlay.File )
        /* Playing is stopped. */
        return;

    /* Retrieve parametres. */
    (void)gomp_GetTrajectoryDisplayParams(&FirstFrame , 
                                          &LastFrame  ,
                                          &FrameStep);

    UsePostRead = 0;
    var = Tcl_GetVar2Ex(
        gomp_GetTclInterp(), "lulPostReadFrameTrigger",
        NULL, TCL_GLOBAL_ONLY);
    if ( var )
        Tcl_GetIntFromObj(NULL, var, &UsePostRead);

    Frame = gomp_GetDisplayFrameNumber();
    if ( changeFrame ) {
        Frame += FrameStep;
        if ( Frame > LastFrame )
            Frame = FirstFrame;
    }

    /* Get the frame. */
    if ( gomp_GetOneFrame(Frame, TrajectoryPlay.File, 0) ) {
        fclose(TrajectoryPlay.File);
        TrajectoryPlay.File = NULL;
        return;
    }

    rewind(TrajectoryPlay.File);

    if ( UsePostRead ) {
        if ( PostReadFrame(Frame) ) {
            fclose(TrajectoryPlay.File);
            TrajectoryPlay.File = NULL;
            return;
        }
    }

    if ( TrajectoryPlay.LoopDisplayDelay ) {
        /* Check that PostReadFrame didn't stop playing the
           trajectory. */
        if ( TrajectoryPlay.ID != MyID )
            return;
    }

    if ( ! TrajectoryPlay.HasTimerProc ) {
        /* Create the timer handler.
         * PlayTrajectoryTimer will normally create it by self
         * but now loop delay is so small that it leaved creation to
         * us. Otherwise us would have been called (we are an idle proc).
         */
        Tcl_CreateTimerHandler(
            TrajectoryPlay.LoopDisplayDelay,
            PlayTrajectoryTimer,
            (ClientData)1); /* change frame */
        TrajectoryPlay.HasTimerProc = 1;
    }

#ifdef ENABLE_GRAPHICS
    gomp_DrawSceneCallback();
#else
    gomp_UpdateData();
#endif
}

static void PlayTrajectoryTimer(ClientData changeFrame)
{
    Tcl_Obj *var;

    TrajectoryPlay.HasTimerProc = 0;

    if ( TrajectoryPlay.HasIdleProc )
        return;

    /* Get loop delay. */
    TrajectoryPlay.LoopDisplayDelay = 0;
    var = Tcl_GetVar2Ex(
        gomp_GetTclInterp(), "gomDisplayLoopSleep",
        NULL, TCL_GLOBAL_ONLY);
    if ( var ) {
        Tcl_GetIntFromObj(NULL, var, &TrajectoryPlay.LoopDisplayDelay);
        if ( TrajectoryPlay.LoopDisplayDelay < 0 )
            TrajectoryPlay.LoopDisplayDelay = 0;
    }

    /* Create the timer handler now so that time used by */
    /* gomp_DrawSceneCallback will be counted.             */
    Tcl_CreateTimerHandler(
        TrajectoryPlay.LoopDisplayDelay,
        PlayTrajectoryTimer, (ClientData)1); /* change frame */

    TrajectoryPlay.HasTimerProc = 1;

    Tcl_DoWhenIdle(PlayTrajectory, changeFrame);

    TrajectoryPlay.HasIdleProc  = 1;
}

/*********************************************************************/
int  gomp_DisplayTrajectory(int Method)
/*********************************************************************/
{
    int   FirstFrame, Frame = 0, LastFrame, FrameStep;
    char   InputText[BUFF_LEN];
    const char *Value;
    int    UsePostRead;
    FILE  *NewTrajectoryFile;

    if ( Method == TRAJ_STOP_LOOP ) {
        /* Stop playing. */
        if ( TrajectoryPlay.File )
            fclose(TrajectoryPlay.File);
        TrajectoryPlay.File = NULL;
        TrajectoryPlay.ID++;
        return(0);
    }

    if ( strlen(gomp_GetTrajectoryFileName()) == 0 )
        return(1);

    (void)gomp_GetTrajectoryDisplayParams(&FirstFrame , 
                                          &LastFrame  ,
                                          &FrameStep);

    if ( gomp_GetTrajectoryFileType() == XMOL_TRAJ      ||
         gomp_GetTrajectoryFileType() == GROMOS96A_TRAJ ||
         gomp_GetTrajectoryFileType() == TINKER_TRAJ    ||
         gomp_GetTrajectoryFileType() == FAMBER_TRAJ    ||
         gomp_GetTrajectoryFileType() == FDL_POLY_TRAJ ) {
        NewTrajectoryFile = fopen(gomp_GetTrajectoryFileName() , "r");
    } else {
        NewTrajectoryFile = fopen(gomp_GetTrajectoryFileName() , "rb");
    }

    if ( NewTrajectoryFile == NULL )
        return(1);

    /* Stop playing. */
    if ( TrajectoryPlay.File )
        fclose(TrajectoryPlay.File);
    TrajectoryPlay.File = NULL;
    TrajectoryPlay.ID++;

/* check for the Tcl variable to see if there is something coming in */
    Value = Tcl_GetVar(gomp_GetTclInterp(),"gomCurrentFrameProperty",TCL_GLOBAL_ONLY);

    if(Value) {
        sscanf(Value,"%f %f %s %s",&FrameNumberDisplay.xp , 
               &FrameNumberDisplay.yp ,
               InputText             ,
               FrameNumberDisplay.FontName);
                                                
        if(gomp_ColourName2RGB(InputText , &FrameNumberDisplay.Red   , 
                          &FrameNumberDisplay.Green , 
                          &FrameNumberDisplay.Blue)) {
            gomp_PrintMessage("?ERROR - can't resolve the colour");
            return(1);
        }
    }
/* post read frame string enables the exe of tcl after reading a frame */
    Value = Tcl_GetVar(gomp_GetTclInterp(),"lulPostReadFrameTrigger",TCL_GLOBAL_ONLY);

    UsePostRead = 0;
    if(Value) {
        if(atoi(Value))
            UsePostRead = 1;
    }

    switch( Method ) {

    case TRAJ_PLAY:
        /* display whole frame as a loop */
        TrajectoryPlay.File = NewTrajectoryFile;
        Tcl_CreateTimerHandler(0, PlayTrajectoryTimer, NULL);
        return(0);

    case TRAJ_FORWARD_FRAME: /* forward frame */
        Frame = gomp_GetDisplayFrameNumber() + FrameStep;
        if( Frame > LastFrame)
            Frame = LastFrame;
        break;

    case TRAJ_BACKWARD_FRAME: /* backward frame */
        Frame = gomp_GetDisplayFrameNumber() - FrameStep;
        if( Frame < FirstFrame)
            Frame = FirstFrame;
        break;

    case TRAJ_FIRST_FRAME: /* first frame */
        Frame = FirstFrame;
        break;

    case TRAJ_LAST_FRAME: /* last  frame */
        Frame = LastFrame;
        break;
    }

    if ( gomp_GetOneFrame(Frame , NewTrajectoryFile , 0) ) {
        fclose(NewTrajectoryFile);
        return(1);
    }

    fclose(NewTrajectoryFile);

    if ( UsePostRead ) {
        if ( PostReadFrame(Frame) )
            return(1);
    }

#ifdef ENABLE_GRAPHICS
    gomp_DrawSceneCallback();
#else
    gomp_UpdateData();
#endif

    return(0);
}

/****************************************************************************/
int gomp_PeekRunningFrameNumberProperty()
/****************************************************************************/
{
    static const char *Value;
    static char  InputText[BUFF_LEN];

    Value = Tcl_GetVar(gomp_GetTclInterp(),"gomCurrentFrameProperty",TCL_GLOBAL_ONLY);

    if(Value) {
        sscanf(Value,"%f %f %s",&FrameNumberDisplay.xp , 
               &FrameNumberDisplay.yp ,
               InputText);
                                                
        if(gomp_ColourName2RGB(InputText , &FrameNumberDisplay.Red   , 
                          &FrameNumberDisplay.Green , 
                          &FrameNumberDisplay.Blue)) {
            gomp_PrintMessage("?ERROR - can't resolve the colour");
            return(1);
        }
    }

    return(0);
}
/****************************************************************************/
int gomp_PutRunningFrameNumber(int Frame)
/****************************************************************************/
{
    static char InputText[BUFF_LEN];
    static const char *Value;

    sprintf(InputText,"%d",Frame);
    Value = Tcl_SetVar(gomp_GetTclInterp(),"gomCurrentFrame",InputText,0);

    if(!Value) {
        gomp_PrintERROR("can't set frame number into variable 'gomCurrentFrame'");
        return(1);
    }

    return(0);
}

/****************************************************************************/
int gomp_DisplayRunningFrameNumber()
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static char InputText[BUFF_LEN];
    static int  FrameNumber;
    static int  mm;

/* look for the display state */
    if(!FrameNumberDisplay.DisplayState) return(0);

/* is there a graphics display attached? */
    if(gomp_GetTermType() != GRAPHICS_AVAILABLE) return(0);

/* check that a trajectory file is defined ... */
    if(!gomp_GetTrajectoryStructureState()) {
        return(0);
    }

    FrameNumber = gomp_GetDisplayFrameNumber();

    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);

    glDisable(GL_LIGHTING);

    glPushMatrix();
    glLoadIdentity();
    glOrtho(0.0, 1.0 , 0.0, 1.0 , 0.0 , 1.0);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(mm);

    glColor3f(FrameNumberDisplay.Red , 
              FrameNumberDisplay.Green , 
              FrameNumberDisplay.Blue);
    glRasterPos3f( FrameNumberDisplay.xp , FrameNumberDisplay.yp , 0.0);
    sprintf(InputText,"Frame #: %d", FrameNumber);
    (void)gomp_PrintString(InputText , FrameNumberDisplay.FontName);
    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(mm);

    glEnable(GL_LIGHTING);
#endif /* ENABLE_GRAPHICS */

    return(0);
}

/****************************************************************************/
int  gomp_SetDisplayRunningFrameNumberState(int State)
/****************************************************************************/
{
    FrameNumberDisplay.DisplayState = State;

    return(0);
}
/****************************************************************************/
int  gomp_GetDisplayRunningFrameNumberState()
/****************************************************************************/
{

    return(FrameNumberDisplay.DisplayState);
}

/****************************************************************************/
static int PostReadFrame(int i)
/****************************************************************************/
{
    static int  ITemp;
    static char Text[BUFF_LEN];

    sprintf(Text,"lulPostReadFrame %d",i);
    ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text);
    if(ITemp != TCL_OK) {
        gomp_PrintERROR("can't execute the post read script 'lulPosReadFrame'");
        return(1);
    }

    return(0);

}

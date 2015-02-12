/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2005 by:
Eero Häkkinen
  


  (1) Version 1.0 was released on 28.11.1997 for Digital, IBM, SGI and
      Windows 95/NT.
  (2) Development of next version started 30.11.1997.
       * Perspective and Orthogonal projection was added

  (3) Version 1.1 was released on 29.01.1998 for Windows 95/NT.

  (4) Version 1.11 was released on 30.03.1998 for Windows 95/NT,
      Alpha NT and Unix (SGI, DEC, Linux).

  (5) Version 1.2 was release on 01.10.1998 for Windows 95/NT.

  (6) Version 1.21 was released on 01.11.1998 for Windows 95/NT and Unix

  (7) Version 1.3 is under work. Release date around March/April 1999 
       * Tcl/Tk 8.1b1 in use 1999-01-02

  (8) Version 1.3a3 was released on 19.5.1999 with Tcl/Tk version 8.1
       * The WinNT/95/98 version was released to the network
       * Source was distributed to some trusted people for Linux versions

  (9) Version 1.30 was released on 23.9.1999 with Tcl/Tk version 8.2.
      The source code was also released but not the Tcl/TK/MESA/GLUT
      libraries.
       * WinNT/95/98
       * Linux

  (10) Version 1.31 was released on 30.10.1999. It corrects some drastic
       bugs in the version 1.3.
        * WinNT/95/98
        * Intel Linux

  (11) Version 1.40 is under construction.
       Relesed on 1.6.2000.

  (12) Next version (1.5 or 2.0) has started (31.05.2000)

  (13) Version 2.0 was released on 31.03.2001.

  (14) Version 2.1 is now in the pipeline.
        * Additions made by Eero Häkkinen CSC during the summer 2001
       Significant contributions in:
        - calculation of hydrogen bonds
        - display protein secondary structure using helixes and arrows

  (15) Version 2.2 is slowly growing 

  (16) Version 2.2 was released 23.9.2002

  (17) Version 2.2.1 

  (18) Version 2.30 was released 23.9.2003

  (19) Version 3.0 was released 23.9.2005

################################################################################

  Bug fix history

  Version 1.4

  1) In the routine "gomp_CreateMolecStructure" in the "data_structure.c"
     file the incoming parameter file name ("Name") had to be copied
     to a static temporary variable. If not copied the value would
     be corrupted.

     The reason is a bit unclear to me but it has to do with the fact that
     there is a complicated calling procedure:
     gOpenMol (parser) ==> tcl-script ==> gOpenMol (parser) ==> tcl-script.
     
     Corrected: 26.05.2000 
     Reported:  LUL
     Fixed:     LUL

  2) A typo error in the gopenmol_gui.tcl file in the data directory.
     The files with the extensions defined do not show in the file
     list.
     Changed line {"GROMACS trajectory" "trn trr xtc TRN TRR XTC" TEXT} to
     {"GROMACS trajectory" ".trn .trr .xtc .TRN .TRR .XTC" TEXT}

     Correct: 26.05.2000
     Reported: Berk Hess
     Fixed:    LUL

  3) In the gopenmol_gui.tcl file there is a file extension .trn defined that
     does not exist. It should read .trj.
     Changed line {"GROMACS trajectory" ".trn .trr .xtc .TRN .TRR .XTC" TEXT} to
     {"GROMACS trajectory" ".trj .trr .xtc .TRJ .TRR .XTC" TEXT}
     and
static      set TrajExt(Gromacs)    ".trn .trr .xtc .TRN .TRR .XTC" to
static      set TrajExt(Gromacs)    ".trj .trr .xtc .TRJ .TRR .XTC"



     A correction has to go also into the rw_get_frame_gromacs.c.
     Changed
     if(Tcl_StringCaseMatch(TempText, "*.trr", 1) || 
        Tcl_StringCaseMatch(TempText, "*.trn", 1)) {
     to
     if(Tcl_StringCaseMatch(TempText, "*.trr", 1) || 
        Tcl_StringCaseMatch(TempText, "*.trj", 1)) {

     Corrected: 26.05.2000
     Reported:  Berk Hess
     Fixed:     LUL

  4) The company behind Tcl/Tk has changed their name, yet an other time. 
     The new name is now Ajuba Solutions. This means that all the links to 
     Tcl/Tk in the gOpenMol documentation is wrong! The new link is:
     Tcl/Tk and documentation (http://dev.ajubasolutions.com/) 

     Corrected: 26.05.2000
     Reported:  LUL
     Fixed:     LUL

     Features added:

  1) It's possible to choose two coloring types:
       RGB coloring or
       grayscale
     This should make it easier to make pictures for publications.

     define coloring color
                     grayscale
  
*/

#include "maindefs.h"

#include <math.h>
#include <setjmp.h>
#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>

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

#ifndef WIN32
#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <X11/cursorfont.h>
#include <X11/Intrinsic.h>
#endif

#if defined(IRIX)
#include <X11/StringDefs.h>
#endif

#endif /* ENABLE_GRAPHICS */

#include <tcl.h>
#ifdef ENABLE_GRAPHICS
#include <tk.h>
#endif /* ENABLE_GRAPHICS */

#include "cluster.h"
#include "colouring.h"
#include "coord_man.h"
#include "drawscene.h"
#include "g_status.h"
#include "gomcast.h"
#include "gomclipbrd.h"
#include "gomfile.h"
#include "gommain.h"
#include "gomstring.h"
#include "gomtext.h"
#include "gomversion.h"
#include "gomwindow.h"
#include "ldp.h"
#include "light_model.h"
#include "listutils.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "objseg.h"
#include "picking.h"
#include "plot_cpk.h"
#include "plot_molec.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "selection.h"
#include "stereo.h"
#include "tclutils.h"
#include "text_stack.h"
#include "trajectory.h"

#include "stdafx.h"

#define YES 1
#define NO  0
#define RABS(a)    ( ( a ) > 0.0 ? (a) : -(a))

#define GRAPHICS_YES 1
#define GRAPHICS_NO  0

#define NEW     0       /* New atom structure will be loaded      */
#define APPEND  1       /* Structure will be appended to old list */

#define COLOR_MATCH   1.0e-3

/*  **** VERY important parameters for the display *** */

#ifdef ENABLE_GRAPHICS
static Tk_Window  gOpenMolWindow;

#if defined(WIN32)
static HWND       GraphicsWindowHandle;
#endif

#if defined(GLUT)

static struct {
    int    WindowingStyle;
    int    Windows;
    int   *WindowStack;
    int   *WindowType;
    char **WindowName;
} WindowHandle;

static int WindowID;
#if 0
static const char *GetWindowNameFromStack(int);
#endif
static int        DeleteWindowFromStack(int);

#endif /* GLUT */
#endif /* ENABLE_GRAPHICS */

static int AllowIndividualScaling = 0;

#ifdef ENABLE_GRAPHICS
static struct {
    int RealBitsPerPixel;
    int UsingBitsPerPixel;
    int MouseButtons;
} SystemInfo;
#endif /* ENABLE_GRAPHICS */

/* Buffer value == 0 is BACK BUFFER, !=0 is FRONT BUFFER */
static struct {
  int Buffer;
} DrawBuffer= {0};


/* Projection transformation */
static struct {
  int Type;
} ProjectionTransformation = {0};

static int  AbortRedrawDisplay;

static int SystemRedisplayMode = 0; /* FAST */
  
/* ... ... ... END ... ... ... */

/* in the PseudoColor mode we have to use a colour map with "n" number
   of entries */

#if 0
static struct {
    int *ColorMap_p;
    int   Entries;
    int   DefinedUp;
    int *Red;
    int *Green;
    int *Blue;
    int   Step;
    char *Name;
} Key2ColorMap = {NULL , 0 , 0 , NULL , NULL , NULL , 50 , NULL};

static int Push2ColorMapKey(int , int , int , int);
static int MakeColorMapModifications(void);
#endif
/* .................................................................... */

#ifdef ENABLE_GRAPHICS
/* define callback procedures */
static int  DrawSceneStereoPair(void);
static int  DrawSceneQuadStereo (void);
static int  ChangeCursor(int);
#endif
/* .............. */

/*  picked atom information */
static struct {
    int Active;
    int Set;
    int HitAtom;
} PickAtom = { 0 , 0 , 0};

static int  SetPickedAtom(int , int);
static const int *GetPickedAtom(void);
static int  PlotPickedAtom(int , int , double);
static int  PlotPickedAtoms(void*, int, int);


static int PickedAtomListLength = 0;

static struct {
    gom_Plotter *CallbackHandle;
    int *Set;
    int *HitAtom;
} PickedAtomList = { NULL, NULL, NULL };

static int  IsAtomPicked(int , int);
static const int *GetPickedAtomListFromStack(void);
static const int *GetPickedStructureListFromStack(void);
static int  FillPickedAtomSegmentResidueAtomList(SegmentResidueAtomList_t*,int,int);

static int UpdateIdentifyAtomWidget(int, int);

#ifdef ENABLE_GRAPHICS
static struct {
    int X;
    int Y;
} MouseLocation; 

static int  MouseButtonState;
#endif

/* define some colourse 'colors.h' */
const float gomp_BLACKv[3]   = {0.0 , 0.0 , 0.0};   /* black    */
const float gomp_REDv[3]     = {1.0 , 0.0 , 0.0};   /* red      */
const float gomp_GREENv[3]   = {0.0 , 1.0 , 0.0};   /* green    */
const float gomp_YELLOWv[3]  = {1.0 , 1.0 , 0.0};   /* yellow   */
const float gomp_BLUEv[3]    = {0.0 , 0.0 , 1.0};   /* blue     */
const float gomp_MAGENTAv[3] = {1.0 , 0.0 , 1.0};   /* magenta  */
const float gomp_CYANv[3]    = {0.0 , 1.0 , 1.0};   /* cyan     */
const float gomp_WHITEv[3]   = {1.0 , 1.0 , 1.0};   /* white    */

#ifdef ENABLE_GRAPHICS
static struct {
    int Width;
    int Height;
} WindowBeforeFullScreen = { 500 , 500 };
#endif /* ENABLE_GRAPHICS */

static char  GlobalTextString[BUFF_LEN];

/* ........................ */

/* ............................................ */

#ifdef ENABLE_GRAPHICS
static int LeftMouseState  = 0;
static int RightMouseState = 0;
static int MiddleMouseState = 0;

#if 0
static int GenericEventHandler(ClientData , XEvent *);
#endif

#if defined(WIN32)
long gOMWndProc(HWND hWnd, UINT message, DWORD wParam, LONG lParam)
{
    printf("Hello leif\n");
    return((long)0);
}
#endif

#if defined(WIN32)
#if defined(GLUT)
static void MouseMotion(int , int);
static void MouseMotionPassive(int , int);
static void HandleMouse(int , int , int , int);

static void LeftMouseTrapDown(int , int);
static void LeftMouseTrapUp(int , int);

static void MiddleMouseTrapDown(int , int);
static void MiddleMouseTrapUp(int , int);

static void RightMouseTrapDown(int , int);
static void RightMouseTrapUp(int , int);
#if 0
static void GetMouseLoc( int * , int * );
#endif

static void HandleSpecialFunc(int, int , int); 

static void TclEventQueue(void);

static void NULLFunction(void);

static struct {
    int State;
} WindowUpdateDisplayMode = { 1 };
#else
static void CALLBACK LeftMouseTrapDown(AUX_EVENTREC *event);
static void CALLBACK LeftMouseTrapUp(AUX_EVENTREC *event);

static void CALLBACK RightMouseTrapDown(AUX_EVENTREC *event);
static void CALLBACK RightMouseTrapUp(AUX_EVENTREC *event);

static void CALLBACK TclEventQueue(AUX_EVENTREC *event);
#endif
#else
#if defined(GLUT)
static void MouseMotion(int , int);
static void MouseMotionPassive(int , int);
static void HandleMouse(int , int , int , int);

static void LeftMouseTrapDown(int , int);
static void LeftMouseTrapUp(int , int);

static void MiddleMouseTrapDown(int , int);
static void MiddleMouseTrapUp(int , int);

static void RightMouseTrapDown(int , int);
static void RightMouseTrapUp(int , int);
#if 0
static void GetMouseLoc( int * , int * );
#endif

static void HandleSpecialFunc(int, int , int); 


static void TclEventQueue(void);

static void NULLFunction(void);

static struct {
    int State;
} WindowUpdateDisplayMode = { 1 };

#else
static void LeftMouseTrapDown(AUX_EVENTREC *event);
static void LeftMouseTrapUp(AUX_EVENTREC *event);

static void RightMouseTrapDown(AUX_EVENTREC *event);
static void RightMouseTrapUp(AUX_EVENTREC *event);

static void TclEventQueue(AUX_EVENTREC *event);
#endif
#endif

static jmp_buf ExitGlutMessageLoopJump;
static int     ExitGlutMessageLoopCount = -1;

/****************************************************************************/
int gomp_Mmain(int argc, const char *argv[])
/*
  int argc; 
  const char *argv[];
*/
/****************************************************************************/
{
    char Temp[BUFF_LEN];

/* some initialization */

    /* create the toplevel shell */

#if defined(USE_TK_STUBS)
    if(Tk_InitStubs(gomp_GetTclInterp() , "8.3" , 0) == NULL)
#else
    if(Tk_Init(gomp_GetTclInterp()) != TCL_OK)
#endif
    {
        Tcl_Interp  *interp;
        interp = gomp_GetTclInterp();
        sprintf(Temp,"'%s' problems creating the Tk window", Tcl_GetStringResult(interp));
        gomp_PrintERROR(Temp);
        return(1);
    }

    gOpenMolWindow = Tk_MainWindow(gomp_GetTclInterp());

    if(gOpenMolWindow == (Tk_Window)NULL) {
        gomp_PrintMessage("can't open graphics display");
        gomp_PrintMessage("will now go into the line command mode");
        gomp_SetTermType(GRAPHICS_NOT_AVAILABLE);
        return(1);
    }

/* put flag on for a graphics display attached ... */
    gomp_SetTermType(GRAPHICS_AVAILABLE);

   SystemInfo.RealBitsPerPixel = Tk_Depth(gOpenMolWindow);

#if !defined(GLUT)
#if defined(WIN32)
    if( SystemInfo.RealBitsPerPixel  < 9) {
       auxInitDisplayMode (AUX_DOUBLE | AUX_RGB | AUX_DEPTH16);
       SystemInfo.UsingBitsPerPixel =  16;
       }
    else {
       auxInitDisplayMode (AUX_DOUBLE | AUX_RGB | AUX_DEPTH24);
       SystemInfo.UsingBitsPerPixel =  24;
       }
#else
       auxInitDisplayMode (AUX_DOUBLE | AUX_RGB | AUX_DEPTH);
#endif
#endif

#if defined(GLUT)
    glutInitWindowPosition( 0 , 0);
    glutInitWindowSize( 500 , 500 );
    glutInit(&argc,CONST_CAST(char**,argv)); 
    if (gomp_GetStereoDisplayState() == STEREO_DISPLAY_ON) {
        glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STEREO);
    } else {
        glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    }
    sprintf(Temp,"(%d) : %s",gomp_GetNumDefinedWindows() + 1,GOPENMOL_INFO);
    WindowID = glutCreateWindow(Temp);
    (void)gomp_PushToWindowStack(WindowID , STRUCTURE_WINDOW , GOPENMOL_INFO);
#else
    auxInitPosition (0, 0, 500, 500);
    auxInitWindow (GOPENMOL_INFO);
#endif

#if defined(WIN32)
#if !defined(GLUT)
    GraphicsWindowHandle = auxGetHWND();
#endif
#endif

/*    auxReshapeFunc (myReshape); */
#if defined(GLUT)
    glutMouseFunc(HandleMouse);
    glutMotionFunc(MouseMotion);
    glutPassiveMotionFunc(MouseMotionPassive);
    glutSpecialFunc(HandleSpecialFunc);

#else
    auxMouseFunc(AUX_LEFTBUTTON  , AUX_MOUSEDOWN , gomp_LeftMouseTrapDown);
    auxMouseFunc(AUX_LEFTBUTTON  , AUX_MOUSEUP   , gomp_LeftMouseTrapUp);

    auxMouseFunc(AUX_RIGHTBUTTON , AUX_MOUSEDOWN , gomp_RightMouseTrapDown);
    auxMouseFunc(AUX_RIGHTBUTTON , AUX_MOUSEUP   , gomp_RightMouseTrapUp);
#endif

#if defined(GLUT)
    glutReshapeFunc(gomp_Reshape);

    if(gomp_GetUpdateDisplayMode())
        glutDisplayFunc(gomp_DrawSceneCallback);
    else
        glutDisplayFunc(NULLFunction);

    glutIdleFunc(TclEventQueue);
#else
    auxReshapeFunc(gomp_Reshape);
    auxExposeFunc(gomp_Reshape);
#endif
/*
  Tk_CreateGenericHandler(gomp_GenericEventHandler , (ClientData) 1);
*/
/* initialize the OpenGL window */
    gomp_WindowInit();

/* read the startup tcl gomp_rip */
    (void)gomp_OpenStartFile();

    if(argc > 1) 
        gomp_OpenInputFiles(argc,argv);

/*  Tk_MainLoop();*/
/* first empty the Tcl/Tk event stack 
   TclEventQueue();
*/

/* now enter the glutMainLoop ...     */
    glutMainLoop();

/* newer reached. */
    return(0);

}

/****************************************************************************/
int gomp_DrawScene()
/****************************************************************************/
{

/* mark the current window to be redisplayed (in fact it redisplays all windows) */

    if(gomp_GetUpdateDisplayMode())
        glutPostRedisplay();
    else
        gomp_DrawSceneCallback();


    return(0);
}

/****************************************************************************/
void gomp_DrawSceneCallback()
/****************************************************************************/
{
    static int    i,j;
    static float  RedC,BlueC,GreenC;
    static int    SaveWindowID;

#if defined(GLUT)
    glutSetWindow(gomp_GetWindowIDFromStack(0));
    (void)ChangeCursor(1);
#else
#if defined(WIN32)
    SetCursor(LoadCursor(NULL , IDC_WAIT));
    ShowCursor(TRUE);
#else
    (void)ChangeCursor(1);
#endif
#endif

/*  look to see if a tcl script should be run... */
/*    (void)gomp_TclRunScript();*/

/*  set the background colour */
    if(gomp_GetWindowingStyle() == MULTI_WINDOWING) {
        SaveWindowID = glutGetWindow();
        for(i = 0 ; i < gomp_GetNumDefinedWindows() ; i++) {
            glutSetWindow(gomp_GetWindowIDFromStack(i));
            (void)gomp_GetBGColor(&RedC , &GreenC , &BlueC);

            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&RedC , &GreenC , &BlueC);

            glClearColor((GLclampf)RedC , 
                         (GLclampf)GreenC , 
                         (GLclampf)BlueC,
                         0.0);

            glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
        }
        glutSetWindow(SaveWindowID);
    } else {
        (void)gomp_GetBGColor(&RedC , &GreenC , &BlueC);

        if(!gomp_GetDisplayColourType())
            (void)gomp_RGB2Grayscale(&RedC , &GreenC , &BlueC);

        glClearColor((GLclampf)RedC , 
                     (GLclampf)GreenC , 
                     (GLclampf)BlueC,
                     0.0);
        glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
    }
/* ....... done ............. */

/* no strcutures defined, do just the basics ... */
    if(!gomp_GetNumMolecStructs()) {
#if defined(GLUT)
        (void)ChangeCursor(0);
#else
#if defined(WIN32)
        ShowCursor(FALSE);
        SetCursor(LoadCursor(NULL , IDC_ARROW));
#else
        (void)ChangeCursor(0);
#endif
#endif

#if defined(GLUT)
        if(gomp_GetWindowingStyle() == MULTI_WINDOWING ) {
            for(i = 0 ; i < gomp_GetNumDefinedWindows() ; i++) {
                glutSetWindow(gomp_GetWindowIDFromStack(i));
                glutSwapBuffers();
            }
        } else {
            glutSwapBuffers();
        }
#else
        auxSwapBuffers();
#endif
        return;
    } /* no molecular structure to display */ 

/* structures are defined ... */

    glLineWidth((GLfloat)1.0);

/* ......... */
    if(gomp_GetDisplayLDPmatrix()) {
#if defined(GLUT)
        if(gomp_GetWindowingStyle()) {
            (void)gomp_UpdateScreen();
            if(!gomp_GetMouseButtonState()) {
                if(gomp_BuildLDParray())
                    return;
            }
        } else {
            (void)gomp_DisplayRunningFrameNumber();
            if(gomp_GetEntriesInTextStack())
                (void)gomp_PlotTextStack();
            if(gomp_BuildLDParray())
                return;
        }
#else
        (void)gomp_DisplayRunningFrameNumber();
        if(gomp_GetEntriesInTextStack())
            (void)gomp_PlotTextStack();
        if(gomp_BuildLDParray()) return;
#endif
    }
    else if(gomp_GetDisplayCLUSTERmatrix()) {
        if(gomp_PreCluster()) return;
    }
    else {
/* check to see if this is a stereo pair plot ... */
        if(gomp_GetStereoPlotState()) {      /* yes it is ... */
            (void)DrawSceneStereoPair();
        }
        else if(gomp_QuadStereoIsOn()) {
            DrawSceneQuadStereo();
        } else {                             /* no it is not ... */
            (void)gomp_UpdateScreen();
        }
    }

/*
  if(!gomp_GetDrawBuffer())
  glXSwapBuffers(XtDisplay(gWindow.gOpenMolWidgetStructureForm),
  XtWindow(gWindow.gOpenMolWidgetStructureForm));
*/
#if defined(GLUT)
    if(!gomp_GetDrawBuffer())
    {
        if((j = gomp_GetNumDefinedWindows())) {
            for(i = 0 ; i < j ; i++) {
                glutSetWindow(gomp_GetWindowIDFromStack(i));
                glutSwapBuffers();
            }
        } else {
            glutSwapBuffers();
        }
    }
#else
    if(!gomp_GetDrawBuffer())
        auxSwapBuffers();
#endif

#if defined(GLUT)
    (void)ChangeCursor(0);
#else
#if defined(WIN32)
    ShowCursor(FALSE);
    SetCursor(LoadCursor(NULL , IDC_ARROW));
#else
    (void)ChangeCursor(0);
#endif
#endif

/*    (void)ChangeCursor(0);*/
/*    (void)gomp_SetDisplayInterrupt(OFF);*/

    return;
}

/****************************************************************************/
int gomp_UpdateScreen()
/****************************************************************************/
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {
    } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {
        gluLookAt(0.0 , 0.0 , 
                  (gomp_GetPerspectiveNear() + 0.5 * gomp_GetPerspectiveWindow()) ,
                  0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0);
    } else {
        gomp_PrintERROR("undefined Projection Transformation defined");
        return(1);
    }

    gomp_UpdateData();

    gomp_CallPlotters(
        gomp_GetMouseButtonState(),      /* are we rotating */
        !gomp_GetSystemRedisplayMode()); /* do we want to draw fast while rotating */

    glLineWidth((GLfloat)1.0);

    (void)gomp_DisplayRunningFrameNumber();

    if(gomp_GetColourScalePlotStatus())
        (void)gomp_DisplayColourScale();

    if(gomp_GetEntriesInTextStack())
        (void)gomp_PlotTextStack();

    glLineWidth((GLfloat)1.0);

    return(0);
}
/****************************************************************************/
int  gomp_PrepareDisplay(int What)
/****************************************************************************/
{
    static float  Size;
    static float  Scale;
    static int    WSize[4];
    static float  RotMS[16];
    static int    i;
/*    static float  Msin = 0.806563;  (1/(2 * sin(45/2)) - 0.5) */
    static float  Msin = 0.306563;   /* (1/(2 * sin(45/2)) - 1.0) */

/*  nothing available to display  */
    if(!gomp_GetNumMolecStructs())
        return(0);
/*  ............................  */

/* Check if GLX is available, if not return gracefully from here */ 
/*  if(gomp_OpenGLsupport()) {*/
/* ............................................................. */
    (void)gomp_ResetPerspectiveWindowAttributes();

    Size = gomp_GetSizeOfSystem();

/* this is an "uggly" way of handling it but I make the assumption that
   the system can have an atom very far from the center of mass so the
   used "size" is twice the real size. */


    if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {

        (void)gomp_SetPerspectiveNear(-Size/2.0);

        (void)gomp_SetPerspectiveWindowAttributes(gomp_GetPerspectiveNear() , Size);
        glGetIntegerv(GL_VIEWPORT , WSize);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        glViewport((GLint)WSize[0],(GLint)WSize[1],
                   (GLsizei)WSize[2],
                   (GLsizei)WSize[3]);

        Scale = (GLfloat)(WSize[2])/
            (GLfloat)(WSize[3]);

        glOrtho(
            (GLdouble)(Scale * gomp_GetPerspectiveNear()),(GLdouble)(Scale * gomp_GetPerspectiveFar()), 
            (GLdouble)gomp_GetPerspectiveNear(),(GLdouble)gomp_GetPerspectiveFar(),
            (GLdouble)gomp_GetPerspectiveNear(),(GLdouble)gomp_GetPerspectiveFar());

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glGetFloatv(GL_MODELVIEW_MATRIX, RotMS);

        if(What) {
            gomp_SaveModelViewMatrixMT(What - 1 , RotMS);
        } else {
            for(i = 0; i < gomp_GetNumMolecStructs() ; i++) {
                gomp_SaveModelViewMatrixMT(i , RotMS);
            }
        }

    } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {

        (void)gomp_SetPerspectiveNear(2.0 * Size * Msin);

        (void)gomp_SetPerspectiveWindowAttributes(gomp_GetPerspectiveNear() , 
                                                (Size + Size));
        glGetIntegerv(GL_VIEWPORT , WSize);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        glViewport((GLint)WSize[0],(GLint)WSize[1],
                   (GLsizei)WSize[2],
                   (GLsizei)WSize[3]);

        Scale = (GLfloat)(WSize[2])/
            (GLfloat)(WSize[3]);

        gluPerspective((GLdouble)gomp_GetPerspectiveAngle(),  (GLdouble)Scale,
                       (GLdouble)gomp_GetPerspectiveNear(),
                       (GLdouble)gomp_GetPerspectiveFar());

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glGetFloatv(GL_MODELVIEW_MATRIX, RotMS);

        if(What) {
            gomp_SaveModelViewMatrixMT(What - 1 , RotMS);
        } else {
            for(i = 0; i < gomp_GetNumMolecStructs() ; i++) {
                gomp_SaveModelViewMatrixMT(i , RotMS);
            }
        }

    } else {
        gomp_PrintERROR("undefined Projection Transformation defined");
        return(1);
    }

    return(0);
}


/****************************************************************************/
int gomp_SetMouseButtonState(int State)
/****************************************************************************/
{
    MouseButtonState = State;

    return(0);
}

/****************************************************************************/
int gomp_GetMouseButtonState()
/****************************************************************************/
{
    return(MouseButtonState);
}
#endif /* ENABLE_GRAPHICS */

/****************************************************************************/
int  gomp_SetDrawBuffer(int Buffer)
/****************************************************************************/
{
    DrawBuffer.Buffer = Buffer;

#ifdef ENABLE_GRAPHICS
    if(!Buffer) 
        glDrawBuffer(GL_BACK);
    else
        glDrawBuffer(GL_FRONT);
#endif /* ENABLE_GRAPHICS */

    return(0);
}
/****************************************************************************/
int  gomp_GetDrawBuffer()
/****************************************************************************/
{
    return(DrawBuffer.Buffer);
}

/****************************************************************************/
int  gomp_SetDisplayInterrupt(int State)
/****************************************************************************/
{
    AbortRedrawDisplay = State;

    return(0);
}

/****************************************************************************/
int  gomp_GetDisplayInterrupt()
/****************************************************************************/
{
    return(AbortRedrawDisplay);
}

#ifdef ENABLE_GRAPHICS
/****************************************************************************/
#if defined(WIN32)
#if defined(GLUT)
void TclEventQueue(void)
#else
void CALLBACK TclEventQueue(AUX_EVENTREC *event)
#endif
#else
#if defined(GLUT)
void TclEventQueue()
#else
void TclEventQueue(AUX_EVENTREC *event)
#endif
#endif
/****************************************************************************/
{
    unsigned int count = 0;

    if ( ExitGlutMessageLoopCount >= 0 ) {
        /* We are peeking the message queue. */
        --ExitGlutMessageLoopCount;
        if ( ExitGlutMessageLoopCount < 0 )
            /* Quit from the inner message loop. */
            longjmp(ExitGlutMessageLoopJump,1);
    }

    /* Process at most 33 Tcl events. */
    while ( Tcl_DoOneEvent(TCL_ALL_EVENTS | TCL_DONT_WAIT) &&
            ! ( count & 0x20 ) )
        count++;

/* fool the callback to sleep 1 ms to prevent 100% cpu utilization */
    Tcl_Sleep(1);
}

/****************************************************************************/
#if defined(WIN32)
#if defined(GLUT)
void gomp_Reshape(GLsizei w, GLsizei h)
#else
void CALLBACK gomp_Reshape(GLsizei w, GLsizei h)
#endif
#else
void gomp_Reshape(GLsizei w, GLsizei h)
#endif
/****************************************************************************/
{
    static float Scale;
    static struct {
        int x;
        int y;
        int width;
        int height;
    } Wattribs; 
/*  static int WSize[4];*/

/*    glGetIntegerv(GL_VIEWPORT , WSize);*/

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    Wattribs.x = 0 /*WSize[0]*/;
    Wattribs.y = 0 /*WSize[1]*/;

    Wattribs.width  = w;
    Wattribs.height = h;

    glViewport((GLint)Wattribs.x,(GLint)Wattribs.y,
               (GLsizei)Wattribs.width ,
               (GLsizei)Wattribs.height);

    Scale = (GLfloat)(Wattribs.width)/
        (GLfloat)(Wattribs.height);

    if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {

        glOrtho(
            (GLdouble)(Scale * gomp_GetPerspectiveNear()),(GLdouble)(Scale * gomp_GetPerspectiveFar()), 
            (GLdouble)gomp_GetPerspectiveNear(),(GLdouble)gomp_GetPerspectiveFar(),
            (GLdouble)gomp_GetPerspectiveNear(),(GLdouble)gomp_GetPerspectiveFar());

    } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {

        gluPerspective((GLdouble)gomp_GetPerspectiveAngle() , (GLdouble)Scale,
                       (GLdouble)gomp_GetPerspectiveNear() , 
                       (GLdouble)gomp_GetPerspectiveFar());

    } else {
        gomp_PrintERROR("undefined Projection Transformation defined");
        return;
    }

    glMatrixMode(GL_MODELVIEW);
/*    (void)gomp_DrawScene ();*/
}

/****************************************************************************/
void gomp_WindowInit()
/****************************************************************************/
{
    static int i,j;
    static float  RedC,BlueC,GreenC;    

/*
  XSelectInput(XtDisplay(w), XtWindow(w),
  XtBuildEventMask(w) | PointerMotionMask);
*/
/*  set the background colour */
    (void)gomp_GetBGColor(&RedC , &BlueC , &GreenC);
    glClearColor((GLclampf)RedC ,
                 (GLclampf)BlueC ,
                 (GLclampf)GreenC,
                 0.0);

    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);

/* ....... done .............. */

/* Use the Widget Size to set aspect ratio */
/*----------------------------------------*/
    glEnable(GL_DEPTH_TEST);    /* use Z buffer      */  

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {

        glOrtho(
            (GLdouble)-1.0,(GLdouble)1.0, 
            (GLdouble)-1.0,(GLdouble)1.0, 
            (GLdouble)-1.0,(GLdouble)1.0);

    } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {

        glFrustum(-1.0 , 1.0 , -1.0 , 1.0 , 10.0 , 30.0);

    } else {
        gomp_PrintERROR("undefined Projection Transformation defined");
        return;
    }
    
    glMatrixMode(GL_MODELVIEW);

    glTranslatef(0.0 , 0.0 , -20.0); 

/*------------------------------------
  Do the setting of lighting
  -------------------------------------*/
    (void)gomp_SetUpLightMaterialModels();

    glEnable(GL_AUTO_NORMAL);
    glEnable(GL_NORMALIZE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LINE_SMOOTH);

    glHint(GL_LINE_SMOOTH_HINT , GL_NICEST);

/*    glEnable(GL_DITHER); */

    (void)gomp_SetSphereQuad();

    glPixelStorei(GL_PACK_ALIGNMENT , 1);

/* .............................. */

/* Check the visual, and print the results */
/*
  (void)gomp_ShowVisual(); 
*/

#if defined(GLUT)
    if((j = gomp_GetNumDefinedWindows())) {
        for(i = 0 ; i < j ; i++) {
            glutSetWindow(gomp_GetWindowIDFromStack(i));
            glutSwapBuffers();
        }
    } else {
        glutSwapBuffers();
    }
#else
    auxSwapBuffers();
#endif

    gomp_DrawScene();
}
#if 0
/****************************************************************************/
int GenericEventHandler(ClientData clientdata, XEvent *eventPtr)
/****************************************************************************/
{
    printf("Generic event handler\n");

    return(0);
}
#endif
#endif /* ENABLE_GRAPHICS */
/****************************************************************************/
int gomp_ResetView()
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    float  Size;
    float  Size2;
    float  Scale;
    int    WSize[4];
    float  RotMS[16];
    const float *RotMP;
    float  i;
    static float  Msin = 0.306563;   /* (1/(2 * sin(45/2)) - 1.0) */

/* check if there is a graphics display associated to the process */
    if(gomp_GetTermType() != GRAPHICS_AVAILABLE) {
        gomp_PrintERROR("no graphics display available");
        return(1);
    }

/* Check if GLX is available, if not return gracefully from here */ 
/*  if(gomp_OpenGLsupport()) { */
/* ............................................................. */

    (void)gomp_ResetPerspectiveWindowAttributes();

    Size = gomp_GetSizeOfSystem();

/* this is an "uggly" way of handling it but I make the assumption that
   the system can have an atom very far from the center of mass so the
   used "size" is twice the real size. */

    glGetIntegerv(GL_VIEWPORT , WSize);
    Scale = (GLfloat)(WSize[2])/
        (GLfloat)(WSize[3]);

    if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {

        Size2 = Size / 2.0;

        (void)gomp_SetPerspectiveNear(-Size2);
/*       (void)gomp_SetPerspectiveFar(Size);*/

        (void)gomp_SetPerspectiveWindowAttributes(gomp_GetPerspectiveNear() , Size);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        glOrtho(
            -(GLdouble)(Scale * Size2),(GLdouble)(Scale * Size2), 
            -(GLdouble)Size2,(GLdouble)Size2, 
            -(GLdouble)Size2,(GLdouble)Size2);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glGetFloatv(GL_MODELVIEW_MATRIX, RotMS);

        for(i = 0; i < gomp_GetNumMolecStructs() ; i++) {
            gomp_SaveModelViewMatrixMT(i , RotMS);
        }
        RotMP = gomp_CalculateGeometricalCenterMT( - 1);

    } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {

        (void)gomp_SetPerspectiveNear(2.0 * Size * Msin);

        (void)gomp_SetPerspectiveWindowAttributes(gomp_GetPerspectiveNear() ,
                                                (Size + Size));

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        gluPerspective((GLdouble)gomp_GetPerspectiveAngle(),  Scale,
                       gomp_GetPerspectiveNear(),
                       gomp_GetPerspectiveFar());
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glGetFloatv(GL_MODELVIEW_MATRIX, RotMS);
 
        for(i = 0; i < gomp_GetNumMolecStructs() ; i++) {
            gomp_SaveModelViewMatrixMT(i , RotMS);
        }
        RotMP = gomp_CalculateGeometricalCenterMT( - 1);

    } else {
        gomp_PrintERROR("undefined Projection Transformation defined");
        return(1);
    }

/*
  }

  gomp_UpdateAtomInputWidget();
*/
    (void)gomp_DrawScene ();
#endif /* ENABLE_GRAPHICS */

    return(0);
}

#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int DrawSceneStereoPair ()
/****************************************************************************/
{


/* ......... */
    (void)gomp_PushModelViewingData();
/* ......... */

    glLineWidth((GLfloat)(gomp_GetMoleculeLineWidth()+0.0));

/* Picture #1 */    
    glMatrixMode(GL_MODELVIEW);

    gomp_Rotate(       gomp_GetStereoPlotAngle()     , 0.0 , 1.0 , 0.0);
    gomp_Translate(   -gomp_GetStereoPlotTranslate() , 0.0 , 0.0);
    
    if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {
    } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {
        gluLookAt(0.0 , 0.0 , 
                  (gomp_GetPerspectiveNear() + 0.5 * gomp_GetPerspectiveWindow()) ,
                  0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0);
    } else {
        gomp_PrintERROR("undefined Projection Transformation defined");
        return(1);
    }

/* ......... */

    (void)gomp_UpdateScreen();

/* ......... */
    (void)gomp_PopModelViewingData();
/* ......... */

/* Picture #2 */    

    gomp_Rotate(       -gomp_GetStereoPlotAngle()     , 0.0 , 1.0 , 0.0 );
    gomp_Translate(     gomp_GetStereoPlotTranslate() , 0.0 , 0.0);
    
    glMatrixMode(GL_MODELVIEW);

    if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {
    } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {
        gluLookAt(0.0 , 0.0 , 
                  (gomp_GetPerspectiveNear() + 0.5 * gomp_GetPerspectiveWindow()) ,
                  0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0);
    } else {
        gomp_PrintERROR("undefined Projection Transformation defined");
        return(1);
    }

/* ......... */


    (void)gomp_UpdateScreen();

/* ......... */
    (void)gomp_PopModelViewingData();
/* ......... */

    return(0);
}
#endif /* ENABLE_GRAPHICS */

/****************************************************************************/
int SetPickedAtom(int Set, int Atom)
/****************************************************************************/
{
    PickAtom.Set     = Set;
    PickAtom.HitAtom = Atom;

    return(0);
}
/****************************************************************************/
const int *GetPickedAtom()
/****************************************************************************/
{
    static int RetVal[2];

    RetVal[0] = PickAtom.Set;
    RetVal[1] = PickAtom.HitAtom;

    return(RetVal);
}
#if 0
/****************************************************************************/
int  Push2ColorMapKey(int Index ,int Red, int Green, int Blue)
/****************************************************************************/
{
    int   i;
    int   Hit;

    if(!gomp_GetColourTableLength()) return(1);

    if(!Key2ColorMap.DefinedUp) {
        Key2ColorMap.Name   = 
            gomp_AllocateCharVector(Key2ColorMap.Step * BUFF_LEN);
        Key2ColorMap.Red    =
            gomp_AllocateIntVector(Key2ColorMap.Step);
        Key2ColorMap.Green  =
            gomp_AllocateIntVector(Key2ColorMap.Step);
        Key2ColorMap.Blue   =
            gomp_AllocateIntVector(Key2ColorMap.Step);
        Key2ColorMap.ColorMap_p =
            gomp_AllocateIntVector(Key2ColorMap.Step);

        Key2ColorMap.DefinedUp += Key2ColorMap.Step;
    }
    else {
        if(Key2ColorMap.Entries == Key2ColorMap.DefinedUp) {
            Key2ColorMap.DefinedUp += Key2ColorMap.Step;
            Key2ColorMap.Name   = gomp_ReallocateCharVector(Key2ColorMap.Name,
                                                Key2ColorMap.DefinedUp * BUFF_LEN);
            Key2ColorMap.Red    = gomp_ReallocateIntVector(Key2ColorMap.Red,
                                                Key2ColorMap.DefinedUp);
            Key2ColorMap.Green  = gomp_ReallocateIntVector(Key2ColorMap.Green,
                                                Key2ColorMap.DefinedUp);
            Key2ColorMap.Blue   = gomp_ReallocateIntVector(Key2ColorMap.Blue,
                                                Key2ColorMap.DefinedUp);
            Key2ColorMap.ColorMap_p   = gomp_ReallocateIntVector(Key2ColorMap.ColorMap_p,
                                                      Key2ColorMap.DefinedUp);
        }
    }

    Hit = 0;
    for(i = 0 ; i < Key2ColorMap.Entries ; i++) {
        if((Key2ColorMap.Red[i]   == Red)   && 
           (Key2ColorMap.Green[i] == Green) &&
           (Key2ColorMap.Blue[i]  == Blue)) {
            Hit = 1;
            break;
        }
    }

    if(!Hit) {
        Key2ColorMap.Red[Key2ColorMap.Entries]        = Red;
        Key2ColorMap.Green[Key2ColorMap.Entries]      = Green;
        Key2ColorMap.Blue[Key2ColorMap.Entries]       = Blue;
        Key2ColorMap.ColorMap_p[Key2ColorMap.Entries] = Index;

        Key2ColorMap.Entries++;
/*
  for(i = 0 ; gomp_GetColourTableLength() ; i++) {
  (void)gomp_ReturnColourFromIndex(i , TempColour , 
  &TempRed , &TempGreen , &TempBlue);
  }
*/
    }

    return(0);
}
/****************************************************************************/
int   MakeColorMapModifications()
/****************************************************************************/
{

/* modify color selection in 'coloursel.c' */
    (void)gomp_ModifyColourSel();

/* modify colour selection in 'contdriver.c' */
/*      (void)gomp_ModifyContourSel(); */


    return(0);
}
#endif
#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int gomp_ScaleDisplay(float scale_x , float scale_y , float scale_z)
/****************************************************************************/
{
    static float RotMS[16];

    if(gomp_GetTermType() != GRAPHICS_AVAILABLE) {
        gomp_PrintERROR("no graphics display available");
        return(1);
    }

    glScalef(scale_x,
             scale_y,
             scale_z);

    glPushMatrix();
    glLoadMatrixf(gomp_GetSavedModelViewMatrixMT(0));
    glScalef(scale_x,
             scale_y,
             scale_z);
    glGetFloatv(GL_MODELVIEW_MATRIX, RotMS);
    gomp_SaveModelViewMatrixMT(0 , RotMS);

    glPopMatrix();

    return(0);
}

/****************************************************************************/
int gomp_ResetProjection()
/****************************************************************************/
{
    static float Scale;
    static int    WSize[4];

    
    if(gomp_GetTermType() != GRAPHICS_AVAILABLE) return(1); /* no graphics */

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glGetIntegerv(GL_VIEWPORT , WSize);

    glViewport((GLint)WSize[0],(GLint)WSize[1],
               (GLsizei)WSize[2],
               (GLsizei)WSize[3]);

    Scale = (GLfloat)(WSize[2])/
        (GLfloat)(WSize[3]);

    if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {

        glOrtho(
            (GLdouble)(Scale * gomp_GetPerspectiveNear()),(GLdouble)(Scale * gomp_GetPerspectiveFar()), 
            (GLdouble)gomp_GetPerspectiveNear(),(GLdouble)gomp_GetPerspectiveFar(),
            (GLdouble)gomp_GetPerspectiveNear(),(GLdouble)gomp_GetPerspectiveFar());

    } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {
        gluPerspective((GLdouble)gomp_GetPerspectiveAngle() , (GLdouble)Scale,
                       (GLdouble)gomp_GetPerspectiveNear() ,
                       (GLdouble)gomp_GetPerspectiveFar());

    } else {
        gomp_PrintERROR("undefined Projection Transformation defined");
        return(1);
    }

    glMatrixMode(GL_MODELVIEW);

/*    (void)gomp_DrawScene ();*/

    return(0);
}

/****************************************************************************/
#if defined(WIN32)
#if defined(GLUT)
void LeftMouseTrapDown(int ix , int iy)
#else
    void CALLBACK LeftMouseTrapDown(AUX_EVENTREC *event)
#endif
#else
#if defined(GLUT)
    void LeftMouseTrapDown(int ix , int iy)
#else
    void LeftMouseTrapDown(AUX_EVENTREC *event)
#endif
#endif
/****************************************************************************/
{
#if defined(GLUT)
    static int   sl;
    static int   HitAtom;
    static int   ITemp; 
    static char  Text[BUFF_LEN];

    LeftMouseState = 1;
    (void)gomp_SetMouseButtonState(ON);

    if(gomp_GetIdentifyAtomActive()) {

/* loop over the structures */
        for(sl = 0 ; sl < gomp_GetNumMolecStructs() ; sl++) {

            if(!gomp_GetSelectedStructure(sl)) continue;

            HitAtom = gomp_IdentifyAtomFromCoords(sl , ix , iy);
            if(HitAtom >= 0) {

                if( !gomp_GetAtomDisplayState(sl,HitAtom) ) {
                    /* Atom is invisible */
                    if( !IsAtomPicked(sl,HitAtom) )
                        continue;
                }

                glDrawBuffer(GL_FRONT);
                (void)PlotPickedAtom(sl , HitAtom , 0.30);
                glFlush();
                glDrawBuffer(GL_BACK);

                (void)SetPickedAtom( sl , HitAtom);
                (void)UpdateIdentifyAtomWidget(sl , HitAtom);

                (void)ChangeCursor(2);

                sprintf(Text,"lulHandleClickedAtomDown %d %d",sl+1,HitAtom+1);

                ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text );
                if(ITemp != TCL_OK) {
                    gomp_PrintERROR("can't execute the script 'lulHandleClickedAtomDown <down>'");
                    return;
                }

                break;
            }
        }

        sprintf(Text,"lulHandlePickAtomEventStart");

        ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text );
        if(ITemp != TCL_OK) {
            gomp_PrintERROR("can't execute the script 'lulHandlePickAtomEventStart'");
            return;
        }

        LeftMouseState = 0;
    }


    return;
#else
#if defined(WIN32)
    static MSG   msg;
#else
    XEvent eventret;
    int x,y,x_root,y_root;
    Window root,child;
#endif

    static GLint dx,dy;
    static float delta_x,delta_y;
    static int newx1,oldx1,newy1,oldy1;
    static int   sl;
    static int   HitAtom;

    LeftMouseState = 1;
    (void)gomp_SetMouseButtonState(ON);

    return;

    oldx1 = event->data[AUX_MOUSEX];
    oldy1 = event->data[AUX_MOUSEY];

    while(LeftMouseState) {

        auxGetMouseLoc( &newx1 , &newy1 );

        if(gomp_GetIdentifyAtomActive()) {

/* loop over the structures */
            for(sl = 0 ; sl < gomp_GetNumMolecStructs() ; sl++) {

                if(!gomp_GetSelectedStructure(sl)) continue;

                HitAtom = gomp_IdentifyAtom(sl , newx1 , newy1);
                if(HitAtom >= 0) {
                    LeftMouseState = 0;
                    (void)SetPickedAtom( sl , HitAtom);
                    (void)UpdateIdentifyAtomWidget(sl , HitAtom);
                    break;
                }
            }
        }
        else {
            dx      = (newx1 - oldx1);
            dy      = (newy1 - oldy1);
            delta_x = (float)dx;
            delta_y = (float)dy;

            if(gomp_GetSelectionModeStatus()) {
                if(gomp_GetRotationState()) {
                    (void)gomp_RotateCoordinates1X(delta_x , 'y');
                    (void)gomp_RotateCoordinates1X(delta_y , 'x');
                }
                else if(gomp_GetTranslationState()) {
                    (void)gomp_TranslateCoordinatesX(
                        delta_x *
                        gomp_GetTranslationDamping() ,
                        -delta_y *
                        gomp_GetTranslationDamping(), 0.0);
                }
            }
            else {
                if(gomp_GetRotationState()) {
                    (void)gomp_Rotate(delta_x , 0.0 , 1.0 , 0.0);
                    (void)gomp_Rotate(delta_y , 1.0 , 0.0 , 0.0);
                }
                else if(gomp_GetTranslationState()) {
                    (void)gomp_Translate( delta_x *
                                        gomp_GetTranslationDamping() ,
                                        -delta_y *
                                        gomp_GetTranslationDamping(), 0.0);
                }
            }
        }
#if defined(WIN32)
        while(PeekMessage(&msg, (HWND) NULL, 0, 0 , PM_REMOVE)) { 
            TranslateMessage(&msg); 
            DispatchMessage(&msg); 
        }
#else
/*
  while(XtAppPending(gomp_GetMainWidgetContext())) {
  XtAppNextEvent(gomp_GetMainWidgetContext(),
  &eventret);
  XtDispatchEvent(&eventret);
  }
*/

        XQueryPointer(auxXDisplay(),auxXWindow(),&root,&child,
                      &x_root,&y_root,&x,&y,&LeftMouseState );

#endif

        if(dx || dy)
            gomp_DrawScene();
        
        oldx1 = newx1;
        oldy1 = newy1;
    }
#endif
}

/****************************************************************************/
#if defined(WIN32)
#if defined(GLUT)
void LeftMouseTrapUp(int ix , int iy)
#else
    void CALLBACK LeftMouseTrapUp(AUX_EVENTREC *event)
#endif
#else
#if defined(GLUT)
    void LeftMouseTrapUp(int ix , int iy)
#else
    void LeftMouseTrapUp(AUX_EVENTREC *event)
#endif
#endif
/****************************************************************************/
{
    static int   ITemp;
    static const int *Picked;
    static int   Structure,AtomIndex;
    static char  Text[BUFF_LEN];
    static char *Text2;
    static char *Structures;
    static char *Atoms;
    static char *SelList;
    static int   NAtoms;
    static SegmentResidueAtomList_t List;

    LeftMouseState = 0;

#if defined(GLUT)
    NAtoms    = gomp_GetPickedAtomListLength();
    Picked    = GetPickedAtom();
    Structure = Picked[0];
    AtomIndex = Picked[1];

    if(gomp_GetIdentifyAtomActive()) {

        if( Structure >= 0 ) {
            sprintf(Text,"lulHandleClickedAtomRelease %d %d",
                    Structure+1,AtomIndex+1);
            ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text );
            if(ITemp != TCL_OK) {
                gomp_PrintERROR("can't execute the script 'lulHandleClickedAtomRelease <up>'");
                return;
            }
            if( gomp_IdentifyAtomFromCoords( Structure , ix , iy ) == AtomIndex ) {
                /* Toggle the picking state. */
                if( IsAtomPicked( Structure , AtomIndex ) )
                    gomp_PopPickedAtomFromList( Structure, AtomIndex );
                else
                    gomp_PushPickedAtomToList( Structure, AtomIndex );

                /* End picking event. */
                strcpy(Text,"lulHandlePickAtomEventEnd {} {} {} {} {}");
                ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text );
                if(ITemp != TCL_OK) {
                    gomp_PrintERROR("can't execute the script 'lulHandlePickAtomEventEnd <up>'");
                    return;
                }
            }
            else {
                /* Do dropping. */
                sprintf(Text,"lulHandlePickAtomEventEnd %d %d {%s} %d %d",
                        Structure+1,AtomIndex+1,
                        gomp_GetAtomSegName(Structure,AtomIndex),
                        gomp_GetAtomResNum1(Structure,AtomIndex),
                        AtomIndex+1);
                ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text );
                if(ITemp != TCL_OK) {
                    gomp_PrintERROR("can't execute the script 'lulHandlePickAtomEventEnd <up>'");
                    return;
                }
            }
        }
        else {
            /* Do dropping. */
            if( gomp_GetPickedAtomListLength() > 0 ) {
                FillPickedAtomSegmentResidueAtomList(
                    &List,LIST_ALL_STRUCTURES,SELECTION_LIST);
                SelList = gomp_GetSegmentResidueAtomList(&List);
                gomp_FreeSegmentResidueAtomList(&List);
            }
            else {
                SelList = gomp_AllocateCharVector(9);
                strcpy(SelList,"{} {} {}");
            }

            Structures = gomp_MakeIndexList(
                gomp_GetPickedAtomListLength(),
                GetPickedStructureListFromStack(),1,0,' ');
            Atoms      = gomp_MakeIndexList(
                gomp_GetPickedAtomListLength(),
                GetPickedAtomListFromStack(),1,0,' ');

            strcpy(Text,"lulHandlePickAtomEventEnd");
            Text2      = gomp_AllocateCharVector( strlen(Text) +
                                    strlen(Structures) + strlen(Atoms) +
                                    strlen(SelList)    + 8 );
            sprintf(Text2,"%s {%s} {%s} %s",Text,Structures,Atoms,SelList);

            gomp_FreeVector(Structures);
            gomp_FreeVector(Atoms);
            gomp_FreeVector(SelList);

            ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), Text2 );
            gomp_FreeVector(Text2);
            if(ITemp != TCL_OK) {
                gomp_PrintERROR("can't execute the script 'lulHandlePickAtomEventEnd <up>'");
                return;
            }
        }
/* reset picked atom */
        (void)SetPickedAtom( -1, -1);
    }
    (void)ChangeCursor(0);
#endif
    (void)gomp_SetMouseButtonState(OFF);
    (void)gomp_DrawScene();
}

/****************************************************************************/
#if defined(WIN32)
#if defined(GLUT)
void RightMouseTrapDown(int ix , int iy)
#else
    void CALLBACK RightMouseTrapDown(AUX_EVENTREC *event)
#endif
#else
#if defined(GLUT)
    void RightMouseTrapDown(int ix , int iy)
#else
    void RightMouseTrapDown(AUX_EVENTREC *event)
#endif
#endif
/****************************************************************************/
{
#if defined(GLUT)
    RightMouseState = 1;
    (void)gomp_SetMouseButtonState(ON);

    return;
#else
#if defined(WIN32)
    static MSG   msg;
#else
    XEvent eventret;
    int x,y,x_root,y_root;
    Window root,child;
#endif

    static GLint dx,dy;
    static float delta_x,delta_y;
    static int   newx1,oldx1,newy1,oldy1;
    static float RotMS[16];
    static float scale_x,scale_y,scale_z;

    RightMouseState = 1;
    (void)gomp_SetMouseButtonState(ON);

    oldx1 = event->data[AUX_MOUSEX];
    oldy1 = event->data[AUX_MOUSEY];

    while(RightMouseState) {

        auxGetMouseLoc( &newx1 , &newy1 );

        dx      = (newx1 - oldx1);
        dy      = (newy1 - oldy1);
        delta_x = (float)dx;
        delta_y = (float)dy;

        scale_z = scale_y = scale_x = 
            RABS((float)oldy1/(float)newy1);

        glScalef(scale_x,
                 scale_y,
                 scale_z);

        glPushMatrix();
        glLoadMatrixf(gomp_GetSavedModelViewMatrixMT(0));
        glScalef(scale_x,
                 scale_y,
                 scale_z);
        glGetFloatv(GL_MODELVIEW_MATRIX, RotMS);
        gomp_SaveModelViewMatrixMT(0 , RotMS);

        glPopMatrix();

#if defined(WIN32)
        while(PeekMessage(&msg, (HWND) NULL, 0, 0 , PM_REMOVE)) { 
            TranslateMessage(&msg); 
            DispatchMessage(&msg); 
        }
#else
/*
  while(XtAppPending(gomp_GetMainWidgetContext())) {
  XtAppNextEvent(gomp_GetMainWidgetContext(),
  &eventret);
  XtDispatchEvent(&eventret);
  }
*/
        XQueryPointer(auxXDisplay(),auxXWindow(),&root,&child,
                      &x_root,&y_root,&x,&y,&RightMouseState );

#endif
    
        if(dx || dy)
            gomp_DrawScene();
        
        oldx1 = newx1;
        oldy1 = newy1;
    }
#endif
}

/****************************************************************************/
#if defined(WIN32)
#if defined(GLUT)
void RightMouseTrapUp(int ix , int iy)
#else
    void CALLBACK RightMouseTrapUp(AUX_EVENTREC *event)
#endif
#else
#if defined(GLUT)
    void RightMouseTrapUp(int ix , int iy)
#else
    void RightMouseTrapUp(AUX_EVENTREC *event)
#endif
#endif
/****************************************************************************/
{
    RightMouseState = 0;
    (void)gomp_SetMouseButtonState(OFF);
    (void)gomp_DrawScene();
}

/****************************************************************************/
#if defined(WIN32)
#if defined(GLUT)
void MiddleMouseTrapDown(int ix , int iy)
#else
    void CALLBACK MiddleMouseTrapDown(AUX_EVENTREC *event)
#endif
#else
#if defined(GLUT)
    void MiddleMouseTrapDown(int ix , int iy)
#else
    void MiddleMouseTrapDown(AUX_EVENTREC *event)
#endif
#endif
/****************************************************************************/
{
#if defined(GLUT)
    MiddleMouseState = 1;
    (void)gomp_SetMouseButtonState(ON);

    return;
#else
#if defined(WIN32)
    static MSG   msg;
#else
    XEvent eventret;
    int x,y,x_root,y_root;
    Window root,child;
#endif

    static GLint dx,dy;
    static float delta_x,delta_y;
    static int   newx1,oldx1,newy1,oldy1;
    static float RotMS[16];

    MiddleMouseState = 1;
    (void)gomp_SetMouseButtonState(ON);

    oldx1 = event->data[AUX_MOUSEX];
    oldy1 = event->data[AUX_MOUSEY];

    while(MiddleMouseState) {

        auxGetMouseLoc( &newx1 , &newy1 );

        dx      = (newx1 - oldx1);
        dy      = (newy1 - oldy1);
        delta_x = (float)dx;
        delta_y = (float)dy;

        (void)gomp_Rotate(delta_x , 0.0 , 0.0 , 1.0);

#if defined(WIN32)
        while(PeekMessage(&msg, (HWND) NULL, 0, 0 , PM_REMOVE)) { 
            TranslateMessage(&msg); 
            DispatchMessage(&msg); 
        }
#else
/*
  while(XtAppPending(gomp_GetMainWidgetContext())) {
  XtAppNextEvent(gomp_GetMainWidgetContext(),
  &eventret);
  XtDispatchEvent(&eventret);
  }
*/
        XQueryPointer(auxXDisplay(),auxXWindow(),&root,&child,
                      &x_root,&y_root,&x,&y,&MiddleMouseState );

#endif
    
        if(dx || dy)
            gomp_DrawScene();
        
        oldx1 = newx1;
        oldy1 = newy1;
    }

#endif
}

/****************************************************************************/
#if defined(WIN32)
#if defined(GLUT)
void MiddleMouseTrapUp(int ix , int iy)
#else
    void CALLBACK MiddleMouseTrapUp(AUX_EVENTREC *event)
#endif
#else
#if defined(GLUT)
    void MiddleMouseTrapUp(int ix , int iy)
#else
    void MiddleMouseTrapUp(AUX_EVENTREC *event)
#endif
#endif
/****************************************************************************/
{
    MiddleMouseState = 0;
    (void)gomp_SetMouseButtonState(OFF);
    (void)gomp_DrawScene();
}

/****************************************************************************/
int gomp_PrepareStatusDisplay(const char *FileName)
/****************************************************************************/
{
    int  FileSize;
    int  NumStruc;
    char OutText[BUFF_LEN];

    FileSize = gomp_GetFileSize(FileName);
    if(FileSize > 0) {
        if(FileSize < 1000)
            sprintf(OutText,"File size: %d bytes",FileSize);
        else 
            sprintf(OutText,"File size: %d Kb",FileSize/1000);
        
        (void)gomp_PutText2Info1(OutText);
    }
    sprintf(OutText,"Atoms: %d",gomp_GetTotalNumberOfAtoms());
    (void)gomp_PutText2Info3(OutText);
          
    NumStruc = gomp_GetNumMolecStructs();
    return(0);
}
#endif /* ENABLE_GRAPHICS */

/****************************************************************************/
int gomp_SetIdentifyAtomActive(int State)
/****************************************************************************/
{
    static const char *ITemp;

    PickAtom.Active = State;

    if(State) {
        ITemp = Tcl_SetVar(gomp_GetTclInterp(),"gomPickSwitch","1",TCL_GLOBAL_ONLY);
        if(!ITemp) {
            gomp_PrintERROR("can't set tcl variable 'gomPickSwitch'");
        }
    } else {
        ITemp = Tcl_SetVar(gomp_GetTclInterp(),"gomPickSwitch","0",TCL_GLOBAL_ONLY);
        if(!ITemp) {
            gomp_PrintERROR("can't set tcl variable 'gomPickSwitch'");
        }
    }

    return(0);
}
/****************************************************************************/
int gomp_GetIdentifyAtomActive()
/****************************************************************************/
{
    return(PickAtom.Active);

}
/*
  send the atom and structure IDs to the GUI
*/
/****************************************************************************/
int UpdateIdentifyAtomWidget(int Set , int AtomID)
/****************************************************************************/
{
    static char OutText[BUFF_LEN];
    static int  Code;

/* check if graphics is available */
    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        sprintf(OutText,"lulPassAtomID2GUI %d %d",Set+1,AtomID+1);
        Code = Tcl_GlobalEval(gomp_GetTclInterp(), OutText);
        if(Code != TCL_OK) {
            gomp_PrintERROR("can't pass atom information to GUI");
            return(1);
        }
    }

    return(0);
}

/****************************************************************************/
int gomp_CopyBitmap2Clipboard()
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    int     WSize[4];

#if defined(WIN32)
    HBITMAP  hBitmap;
    HWND     hwnd;
    HDC      hdc,hdcMem;
#endif 

/* check for the graphics ... */

    if(gomp_GetTermType() != GRAPHICS_AVAILABLE)
        return(0);

    glGetIntegerv(GL_VIEWPORT , WSize);

#if defined(WIN32)
/*

The BITMAP structure defines the type, width, height, color format, and bit values of a bitmap. 

typedef struct tagBITMAP {  // bm  
LONG   bmType; 
LONG   bmWidth; 
LONG   bmHeight; 
LONG   bmWidthBytes; 
WORD   bmPlanes; 
WORD   bmBitsPixel; 
LPVOID bmBits; 
} BITMAP; 
 
Members
bmType
Specifies the bitmap type. This member must be zero. 
bmWidth
Specifies the width, in pixels, of the bitmap. The width must be greater 
than zero. 
bmHeight
Specifies the height, in pixels, of the bitmap. The height must be greater 
than zero. 
bmWidthBytes
Specifies the number of bytes in each gomp_an line. This value must be 
divisible by 2, because Windows assumes that the bit values of a bitmap 
form an array that is word aligned. 
bmPlanes
Specifies the count of color planes. 
bmBitsPixel
Specifies the number of bits required to indicate the color of a pixel. 
bmBits
Points to the location of the bit values for the bitmap. The bmBits member 
must be a long pointer to an array of character (1-byte) values. 
 
Remarks
The bitmap formats currently used are monochrome and color. The monochrome 
bitmap uses a one-bit, one-plane format. Each gomp_an is a multiple of 32 bits. 
Scans are organized as follows for a monochrome bitmap of height n: 
Scan 0  
Scan 1 
. 
. 
. 
Scan n-2 
Scan n-1 
 
The pixels on a monochrome device are either black or white. If the 
corresponding bit in the bitmap is 1, the pixel is set to the foreground 
color; if the corresponding bit in the bitmap is zero, the pixel is set 
to the background color. 

*/
/*
  hwnd   = auxGetHWND();
  hdc    = GetDC(hwnd);

  hBitmap = CreateBitmapIndirect(&bitmap);

  printf("%d %d %d %d %d %d %d\n",
  bitmap.bmType, 
  bitmap.bmWidth, 
  bitmap.bmHeight, 
  bitmap.bmWidthBytes, 
  bitmap.bmPlanes, 
  bitmap.bmBitsPixel, 
  bitmap.bmBits); 

  bitmap.bmType       =             0; 
  bitmap.bmWidth      =      WSize[2]; 
  bitmap.bmHeight     =      WSize[3]; 
  bitmap.bmWidthBytes =  3 * WSize[2]; 
  bitmap.bmPlanes     =             3; 
  bitmap.bmBitsPixel  =            24; 
  bitmap.bmBits       =        pixels; 
*/
#if defined(GLUT)
    gomp_SetWindow(1);
    hwnd   = GetActiveWindow();
/*  hdc    = GetDC(hwnd);        */
    hdc    = wglGetCurrentDC();
#else
    hwnd   = auxGetHWND();
    hdc    = auxGetHDC();
#endif
    hdcMem = CreateCompatibleDC(hdc);

    hBitmap = CreateCompatibleBitmap(hdc , WSize[2] , WSize[3]);

    SelectObject(hdcMem , hBitmap);

    BitBlt(hdcMem, 0 , 0 , WSize[2] , WSize[3] , hdc , 0 , 0 , SRCCOPY);

    OpenClipboard(hwnd);
    EmptyClipboard();
    SetClipboardData(CF_BITMAP , hBitmap);
    CloseClipboard();

    DeleteDC(hdcMem);
    DeleteObject(hBitmap);

    return(0);
#else
    gomp_PrintERROR("the clipboard is not available");
    return(1);
#endif
#else /* !ENABLE_GRAPHICS */
    return(0);
#endif
}

/****************************************************************************/
int   gomp_DeleteClipboardData()
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
#if defined(WIN32)
    HWND     hwnd;
#endif
/* check for the graphics ... */

    if(gomp_GetTermType() != GRAPHICS_AVAILABLE)
        return(0);

#if defined(WIN32)
#if defined(GLUT)
    hwnd   = GetActiveWindow();
#else
    hwnd   = auxGetHWND();
#endif
    OpenClipboard(hwnd);
    EmptyClipboard();
    CloseClipboard();

    return(0);
#else
    gomp_PrintERROR("the clipboard is not available");
    return(1);
#endif
#else /* !ENABLE_GRAPHICS */
    return(0);
#endif
}

/****************************************************************************/
int gomp_CopyText2Clipboard(const char *Text)
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
#if defined(WIN32)
    int      i;
    HWND     hwnd;
    HGLOBAL  hGlobalMemory;
    char *pGlobalMemory;

#endif 
/* check for the graphics ... */

    if(gomp_GetTermType() != GRAPHICS_AVAILABLE)
        return(0);

#if defined(WIN32)

    hGlobalMemory = GlobalAlloc(GHND , (strlen(Text) + 1));
    if(hGlobalMemory == NULL) {
        gomp_PrintERROR("can't allocate clipboard memory");
        return(1);
    }

    pGlobalMemory = GlobalLock(hGlobalMemory);
        
    for(i = 0 ; i < (int)strlen(Text); i++)
        pGlobalMemory[i] = Text[i];

    GlobalUnlock(hGlobalMemory);

#if defined(GLUT)
    hwnd   = GetActiveWindow();
#else
    hwnd   = auxGetHWND();
#endif

    OpenClipboard(hwnd);
    EmptyClipboard();
    SetClipboardData(CF_TEXT , hGlobalMemory);
    CloseClipboard();
    GlobalFree(hGlobalMemory);


    return(0);
#else
    gomp_PrintERROR("the clipboard is not available");
    return(1);
#endif
#else /* !ENABLE_GRAPHICS */
    return(0);
#endif
}
#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int gomp_ResizeWindow(int Width, int Height)
/****************************************************************************/
{

#if defined(OBSOLATE)

#if defined(WIN32)
    int      WSize[4];
    HWND     hwnd;

#if defined(GLUT)
    hwnd   = GetActiveWindow();
#else /* !GLUT */
    hwnd   = auxGetHWND();
#endif /* GLUT */

    glGetIntegerv(GL_VIEWPORT , WSize);

    (void)SetWindowPos(
        hwnd,          /* handle of window            */
        HWND_TOP,      /* placement-order handle      */
        0 ,           /* horizontal position         */
        0 ,           /* vertical position           */
/* reserve extra space for borders ...                */
        Width +  2 * GetSystemMetrics(SM_CXFRAME),
        /* width                       */
/* reserve extra space for borders and caption bar... */
        Height + 2 * GetSystemMetrics(SM_CYFRAME) + GetSystemMetrics(SM_CYCAPTION) ,
        /* height                      */
        SWP_NOMOVE|SWP_DRAWFRAME|SWP_FRAMECHANGED|SWP_SHOWWINDOW
        /* window-positioning flags    */
        );

/*
  (void)MoveWindow(hwnd , WSize[0] , WSize[1] , Width , Height , TRUE);
*/
#else /* !WIN32 */
    gomp_PrintERROR("not supported using AUX library");
    return(1);
#endif /* WIN32 */

    (void)gomp_Reshape( Width , Height);
#else /* !OBSOLATE */
#ifndef GLUT
    gomp_PrintERROR("not supported using AUX library");
    return(1);
#else
    glutReshapeWindow(Width , Height);
    (void)gomp_Reshape( Width , Height);
#endif /* GLUT */
#endif /* OBSOLATE */

    return(0);
}

/****************************************************************************/
int gomp_MoveWindow(int Xc, int Yc)
/****************************************************************************/
{

    int      WSize[4];

#if defined(OBSOLATE)

#if defined(WIN32)
    HWND     hwnd;
    HDC      hdc;

#if defined(GLUT)
    hwnd   = GetActiveWindow();
/*  hdc    = GetDC(hwnd);        */
    hdc    = wglGetCurrentDC();
#else
    hwnd   = auxGetHWND();
    hdc    = GetDC(hwnd);
#endif

    glGetIntegerv(GL_VIEWPORT , WSize);


    (void)MoveWindow(hwnd , Xc, 
                     abs(GetDeviceCaps(hdc,VERTRES) - Yc), 
                     WSize[2] , WSize[3] , TRUE);

    ReleaseDC(hwnd,hdc);

#else
    gomp_PrintERROR("not supported using AUX library");
    return(1);
#endif

    (void)gomp_Reshape( WSize[2] , WSize[3]);
#else
    glutPositionWindow(Xc , Yc);
    glGetIntegerv(GL_VIEWPORT , WSize);
    (void)gomp_Reshape( WSize[2] , WSize[3]);
#endif

    return(0);
}

/****************************************************************************/
int gomp_IconifyWindow()
/****************************************************************************/
{


#if defined(GLUT)

    glutIconifyWindow();
    return(0);
#else
    gomp_PrintERROR("not supported using AUX library");
    return(1);
#endif

}
/****************************************************************************/
int gomp_DeIconifyWindow()
/****************************************************************************/
{


#if defined(GLUT)

    glutShowWindow();
    return(0);
#else
    gomp_PrintERROR("not supported using AUX library");
    return(1);
#endif

}

/****************************************************************************/
int gomp_FullScreenWindow()
/****************************************************************************/
{

#if defined(GLUT)

    int      WSize[4];

    glGetIntegerv(GL_VIEWPORT , WSize);

    WindowBeforeFullScreen.Width  = WSize[0];
    WindowBeforeFullScreen.Height = WSize[1];

    glutFullScreen();

    return(0);
#else
    gomp_PrintERROR("not supported using AUX library");
    return(1);
#endif

}

/****************************************************************************/
int gomp_PopWindow()
/****************************************************************************/
{

#if defined(GLUT)

    glutPopWindow();

    return(0);
#else
    gomp_PrintERROR("not supported using AUX library");
    return(1);
#endif

}

/****************************************************************************/
int  gomp_GetWindowInfo(int *Xpos, int *Ypos, int *Width, int *Height)
/****************************************************************************/
{
    int      WSize[4];

#if defined(WIN32)
    HWND     hwnd;
    RECT     rect;
#endif

/* check that there is a graphics available */
    if(gomp_GetTermType() != GRAPHICS_AVAILABLE) {
        *Xpos  =  0;
        *Ypos  =  0;
        *Width =  0;
        *Height = 0;
        return(1);
    }

#if defined(WIN32)
#if defined(GLUT)
    hwnd   = GetActiveWindow();
#else
    hwnd   = auxGetHWND();
#endif

    (void)GetWindowRect(hwnd, /* handle of window                            */
                        &rect /* address of structure for window coordinates */
        );

/*printf("(1) %d %d %d %d\n",rect.left,rect.right,rect.top,rect.bottom);*/
#endif

    glGetIntegerv(GL_VIEWPORT , WSize);
    *Xpos  =  WSize[0];
    *Ypos  =  WSize[1];
    *Width =  WSize[2];
    *Height = WSize[3];

    return(0);
}
/*
  Empty the message queue ...
  1998-07-09
*/
/****************************************************************************/
int gomp_PeekMessageQueue()
/****************************************************************************/
{
#if defined(WIN32)
    static MSG   msg;
#endif

    if ( setjmp(ExitGlutMessageLoopJump) == 0 ) {
        /* Exit glut message queue after 3 idle calls. */
        ExitGlutMessageLoopCount = 3;
        glutMainLoop();
    }

    return(0);
}

#if defined(WIN32)
#if defined(GLUT)
/****************************************************************************/
int ChangeCursor(int State)
/****************************************************************************/
{
    switch(State) {

    case 0: 
        /* set default cursor back */
        glutSetWindow(WindowID);
        glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
        break;
    case 1:   /* change to predefined cursor */
        glutSetWindow(WindowID);
        glutSetCursor(GLUT_CURSOR_WAIT);
        break;
    case 2:   /* change cursor to "clicked" atom */
        glutSetWindow(WindowID);
        glutSetCursor(GLUT_CURSOR_TEXT);
        break;

    }

    return(0);
}    
#endif
#else
#if !defined(GLUT)
/****************************************************************************/
Display    *gomp_GetStructureDisplay()
/****************************************************************************/
{
    return(auxXDisplay());
}
/****************************************************************************/
XtAppContext  gomp_GetMainWidgetContextObsolate()
/****************************************************************************/
{ 
    return(XtDisplayToApplicationContext(auxXDisplay()));
} 
#endif
/****************************************************************************/
int ChangeCursor(int State)
/****************************************************************************/
{
#ifndef GLUT
    Cursor c1;
#endif
    switch(State) {

    case 0: 
        /* set default cursor back */
#if defined(GLUT)
        glutSetWindow(WindowID);
        glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
#else
        XUndefineCursor(auxXDisplay(),auxXWindow());
        XFlush(auxXDisplay());
#endif
        break;
    case 1:   /* change to predefined cursor */
#if defined(GLUT)
        glutSetWindow(WindowID);
        glutSetCursor(GLUT_CURSOR_WAIT);
#else
        c1 = XCreateFontCursor(auxXDisplay() , XC_watch);
        XDefineCursor(auxXDisplay(),auxXWindow(),c1);
        XFlush(auxXDisplay());
#endif
        break;
    case 2:   /* change cursor to "clicked" atom */
        glutSetWindow(WindowID);
        glutSetCursor(GLUT_CURSOR_TEXT);
        break;
    }

    return(0);
}    

#endif
#endif /* ENABLE_GRAPHICS */

/****************************************************************************/
int gomp_SetProjectionTransformation(int Type)
/****************************************************************************/
{
    ProjectionTransformation.Type = Type;

    return(0);
}
/****************************************************************************/
int gomp_GetProjectionTransformation()
/****************************************************************************/
{

    return(ProjectionTransformation.Type);
}

#ifdef ENABLE_GRAPHICS
#if defined(GLUT)
/****************************************************************************/
#if defined(WIN32)
void HandleMouse(int button, int state, int x, int y)
#else
void HandleMouse(int button, int state, int x, int y)
#endif
/****************************************************************************/
{

    if(button == GLUT_LEFT_BUTTON) {
        if(state == GLUT_DOWN) {
            if(glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
                if(gomp_SetRotationState(OFF)) {
                    gomp_PrintERROR("can't change manipulation state");
                }
            }
            else if(glutGetModifiers() == GLUT_ACTIVE_CTRL) {
                if(gomp_SetIdentifyAtomActive(ON)) {
                    gomp_PrintERROR("can't activate atom picking");
                }
            }
            else if(glutGetModifiers() == GLUT_ACTIVE_ALT) {
                MiddleMouseTrapDown(x , y);
                return;
            }
            LeftMouseTrapDown(x , y);
        } else if (state == GLUT_UP) {
            if(glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
                if(gomp_SetRotationState(ON)) {
                    gomp_PrintERROR("can't change manipulation state");
                }
            }
            else if(glutGetModifiers() == GLUT_ACTIVE_CTRL) {
                LeftMouseTrapUp(x , y);
                if(gomp_SetIdentifyAtomActive(OFF)) {
                    gomp_PrintERROR("can't activate atom picking");
                }
                return;
            }
            else if(glutGetModifiers() == GLUT_ACTIVE_ALT) {
                MiddleMouseTrapUp(x , y);
            }

            LeftMouseTrapUp(x , y);
            MiddleMouseTrapUp(x , y);
        }
    } else if (button == GLUT_RIGHT_BUTTON) {
        if(state == GLUT_DOWN) {
            RightMouseTrapDown(x , y);
        } else if (state == GLUT_UP) {
            RightMouseTrapUp(x , y);
        }
    }  else if (button == GLUT_MIDDLE_BUTTON) {
        if(state == GLUT_DOWN) {
            MiddleMouseTrapDown(x , y);
        } else if (state == GLUT_UP) {
            MiddleMouseTrapUp(x , y);
        }
    }
}
#if 0
/****************************************************************************/
void GetMouseLoc( int *x , int *y )
/****************************************************************************/
{
    *x = MouseLocation.X;
    *y = MouseLocation.Y;
}
#endif
/****************************************************************************/
void MouseMotion(int x, int y)
/****************************************************************************/
{

    static GLint dx,dy;
    static float delta_x,delta_y;
    static float scale_x,scale_y,scale_z;
    static int   newx1,oldx1,newy1,oldy1;
    static float RotMS[16];
    static int   i;

    newx1 = x;
    newy1 = y;
    oldx1 = MouseLocation.X;
    oldy1 = MouseLocation.Y;

/* if left mouse button is pressed */
    if(LeftMouseState) {

        dx      = (newx1 - oldx1);
        dy      = (newy1 - oldy1);
        delta_x = (float)dx;
        delta_y = (float)dy;

        if(gomp_GetSelectionModeStatus()) {
            if(gomp_GetRotationState()) {
                (void)gomp_RotateCoordinates1X(delta_x , 'y');
                (void)gomp_RotateCoordinates1X(delta_y , 'x');
            }
            else if(gomp_GetTranslationState()) {
                (void)gomp_TranslateCoordinatesX(
                    delta_x *
                    gomp_GetTranslationDamping() ,
                    -delta_y *
                    gomp_GetTranslationDamping(), 0.0);
            }
        }
        else {
            if(gomp_GetRotationState()) {
                (void)gomp_Rotate(delta_x , 0.0 , 1.0 , 0.0);
                (void)gomp_Rotate(delta_y , 1.0 , 0.0 , 0.0);
            }
            else if(gomp_GetTranslationState()) {
                (void)gomp_Translate( delta_x *
                                    gomp_GetTranslationDamping() ,
                                    -delta_y *
                                    gomp_GetTranslationDamping(), 0.0);
            }
        }

        if(dx || dy)
            gomp_DrawScene();
    }

/* if right mouse button is pressed */
    if(RightMouseState) {

        newx1 = x;
        newy1 = y;
        oldx1 = MouseLocation.X;
        oldy1 = MouseLocation.Y;

        dx      = (newx1 - oldx1);
        dy      = (newy1 - oldy1);
        delta_x = (float)dx;
        delta_y = (float)dy;

        scale_z = scale_y = scale_x = 
            RABS((float)oldy1/(float)newy1);


        for(i = 0 ; i < gomp_GetNumMolecStructs(); i++) {

            if(gomp_GetSelectedStructure(i) || !gomp_GetAllowIndividualScaling()) { 
                glPushMatrix();
                glLoadMatrixf(gomp_GetSavedModelViewMatrixMT(i));
                glScalef(scale_x,
                         scale_y,
                         scale_z);
                glGetFloatv(GL_MODELVIEW_MATRIX, RotMS);
                gomp_SaveModelViewMatrixMT(i , RotMS);

                glPopMatrix();
            }
        }
       
        if(dx || dy)
            gomp_DrawScene();

    }

/* if middle mouse button is pressed */
    if(MiddleMouseState) {

        newx1 = x;
        newy1 = y;
        oldx1 = MouseLocation.X;
        oldy1 = MouseLocation.Y;

        dx      = (newx1 - oldx1);
        dy      = (newy1 - oldy1);
        delta_x = (float)dx;
        delta_y = (float)dy;

        (void)gomp_Rotate(delta_x , 0.0 , 0.0 , 1.0);

    
        if(dx || dy)
            gomp_DrawScene();   
    }

    MouseLocation.X = x;
    MouseLocation.Y = y;
}

/****************************************************************************/
#if defined(WIN32)
void MouseMotionPassive(int x, int y)
#else
void MouseMotionPassive(int x, int y)
#endif
/****************************************************************************/
{
    static int   sl;
    static int   HitAtom;


    MouseLocation.X = x;
    MouseLocation.Y = y;

    if(gomp_GetIdentifyAtomActive()) {
/* loop over the structures */
        for(sl = 0 ; sl < gomp_GetNumMolecStructs() ; sl++) {

            if(!gomp_GetSelectedStructure(sl)) continue;

            HitAtom = gomp_IdentifyAtomFromCoords(sl , x , y);
            if(HitAtom >= 0) {
                if( !gomp_GetAtomDisplayState(sl,HitAtom) ) {
                    /* Atom is invisible */
                    if( !IsAtomPicked(sl,HitAtom) )
                        continue;
                }

                glDrawBuffer(GL_FRONT);
                (void)PlotPickedAtom(sl , HitAtom , 0.15);
                glFlush();
                glDrawBuffer(GL_BACK);

                (void)UpdateIdentifyAtomWidget(sl , HitAtom);
                return;
            }
        }
    }
}
/****************************************************************************/
int gomp_GetWindowID()
/****************************************************************************/
{
    return(glutGetWindow());
}
/****************************************************************************/
int gomp_SetWindow(int Window)
/****************************************************************************/
{
    if(!WindowHandle.Windows) {
        gomp_PrintERROR("no windows so far defined!");
        return(1);
    }
    if(Window < 1 || Window > WindowHandle.Windows) {
        gomp_PrintERROR("windows entry out of allowed range");
        return(1);
    }

    glutSetWindow(gomp_GetWindowIDFromStack(Window - 1));

    return(0);
}
/****************************************************************************/
int gomp_PushToWindowStack(int WindowID, int WindowType, const char *WindowName)
/****************************************************************************/
{
    if(!WindowHandle.Windows) {

        WindowHandle.WindowStack    = gomp_AllocateIntVector(1);
        WindowHandle.WindowStack[0] = WindowID;
        WindowHandle.WindowType     = gomp_AllocateIntVector(1);
        WindowHandle.WindowType[0]  = WindowType;
        WindowHandle.WindowName     = malloc(sizeof(const char *));

        if(WindowHandle.WindowName == NULL) {
            gomp_PrintERROR("can't allocate memory in 'push to window stack'");
            return(1);
        }

        WindowHandle.WindowName[0]  = malloc(strlen(WindowName) + 1);

        if(WindowHandle.WindowName[0] == NULL) {
            gomp_PrintERROR("can't allocate memory in 'push to window stack'");
            return(1);
        }

        strncpy(WindowHandle.WindowName[0] , WindowName , strlen(WindowName));
        WindowHandle.Windows = 1;
    } else {
        WindowHandle.WindowStack = gomp_ReallocateIntVector(WindowHandle.WindowStack ,
                                                 WindowHandle.Windows + 1);
        WindowHandle.WindowStack[WindowHandle.Windows] = WindowID;
        WindowHandle.WindowType  = gomp_ReallocateIntVector(WindowHandle.WindowType ,
                                                 WindowHandle.Windows + 1);
        WindowHandle.WindowType[WindowHandle.Windows]  = WindowType;
        WindowHandle.WindowName     = realloc(
            WindowHandle.WindowName , 
            (WindowHandle.Windows + 1) * 
            sizeof(const char *));

        if(WindowHandle.WindowName == NULL) {
            gomp_PrintERROR("can't allocate memory in 'push to window stack'");
            return(1);
        }

        WindowHandle.WindowName[WindowHandle.Windows]  = 
            malloc(strlen(WindowName) + 1);

        if(WindowHandle.WindowName[WindowHandle.Windows] == NULL) {
            gomp_PrintERROR("can't allocate memory in 'push to window stack'");
            return(1);
        }

        strncpy(WindowHandle.WindowName[WindowHandle.Windows] , 
                WindowName , strlen(WindowName));
        WindowHandle.Windows++;
    }

    return(0);
}
/*
  Entries are indexed from 0 ... N-1
*/
/****************************************************************************/
int gomp_GetWindowIDFromStack(int WEntry)
/****************************************************************************/
{
    if(!WindowHandle.Windows) {
        gomp_PrintERROR("no windows so far defined!");
        return(-1);
    }
    if(WEntry < 0 || WEntry > (WindowHandle.Windows - 1)) {
        gomp_PrintERROR("windows entry out of allowed range");
        return(-1);
    }

    return(WindowHandle.WindowStack[WEntry]);
}
/*
  Entries are indexed from 0 ... N-1
*/
/****************************************************************************/
int gomp_GetWindowTypeFromStack(int WEntry)
/****************************************************************************/
{
    if(!WindowHandle.Windows) {
        gomp_PrintERROR("no windows so far defined!");
        return(-1);
    }
    if(WEntry < 0 || WEntry > (WindowHandle.Windows - 1)) {
        gomp_PrintERROR("windows entry out of allowed range!");
        return(-1);
    }

    return(WindowHandle.WindowType[WEntry]);
}
/****************************************************************************/
int gomp_GetNumDefinedWindows()
/****************************************************************************/
{
    return(WindowHandle.Windows);
}
/*
  Entries are indexed from 0 ... N-1
*/
/****************************************************************************/
int        DeleteWindowFromStack(int Which)
/****************************************************************************/
{
    int i;

    if(!WindowHandle.Windows) {
        gomp_PrintERROR("no windows so far defined!");
        return(1);
    }
    if(Which < 0 || Which > (WindowHandle.Windows - 1)) {
        gomp_PrintERROR("windows entry out of allowed range!");
        return(1);
    }

    glutDestroyWindow(WindowHandle.WindowStack[Which]);

    if(Which != (WindowHandle.Windows - 1)) {

        for(i = Which + 1 ; i < WindowHandle.Windows; i++) {

            WindowHandle.WindowStack[i - 1] = WindowHandle.WindowStack[i];
            WindowHandle.WindowType[i - 1]  = WindowHandle.WindowType[i];
            WindowHandle.WindowName[i - 1]  = 
                realloc(WindowHandle.WindowName[i - 1] , strlen(WindowHandle.WindowName[i]) + 1);

            if(WindowHandle.WindowName[i - 1] == (const char *)NULL) {
                gomp_PrintERROR("can't allocate memory in 'push to window stack'");
                return(1);
            }
               
            strncpy(WindowHandle.WindowName[i - 1] , 
                    WindowHandle.WindowName[i]     ,
                    strlen(WindowHandle.WindowName[i]));
        }
    }

    WindowHandle.WindowStack = gomp_ReallocateIntVector(WindowHandle.WindowStack ,
                                             WindowHandle.Windows - 1);
    WindowHandle.WindowType  = gomp_ReallocateIntVector(WindowHandle.WindowType ,
                                             WindowHandle.Windows - 1);
    free(WindowHandle.WindowName[WindowHandle.Windows - 1]);
    WindowHandle.WindowName = realloc(WindowHandle.WindowName ,
                                      (WindowHandle.Windows - 1) * sizeof(const char *));

    if(WindowHandle.WindowName == NULL) {
        gomp_PrintERROR("can't allocate memory in 'push to window stack'");
        return(1);
    }
        
    WindowHandle.Windows--;

    return(0);
}
/*
  Style == 0 , single windowing
  Style != 0 , multi windowing
*/
/****************************************************************************/
int gomp_SetWindowingStyle(int Style)
/****************************************************************************/
{
    int i;

    WindowHandle.WindowingStyle = Style;

    if(Style == SINGLE_WINDOWING) {
        for(i = 1 ; i < WindowHandle.Windows ; i++) {
            if(DeleteWindowFromStack(i)) {
                gomp_PrintERROR("can not delete window from stack");
                return(1);
            }
        }
    }

    return(0);
}
/****************************************************************************/
int gomp_GetWindowingStyle()
/****************************************************************************/
{
    return(WindowHandle.WindowingStyle);
}
#if 0
/****************************************************************************/
const char *GetWindowNameFromStack(int Which)
/****************************************************************************/
{
    if(!WindowHandle.Windows) {
        gomp_PrintERROR("no windows so far defined!");
        return((const char *)NULL);
    }
    if(Which < 0 || Which > (WindowHandle.Windows - 1)) {
        gomp_PrintERROR("windows entry out of allowed range!");
        return((const char *)NULL);
    }

    return(WindowHandle.WindowName[Which]);
}
#endif
/****************************************************************************/
void HandleSpecialFunc(int key, int x, int y)
/****************************************************************************/
{
/* all Fx keys ... */
    if(key == GLUT_KEY_F1) {
        if(TCL_OK != gomp_SendCommand2Parser("eval $F1")){
            gomp_PrintERROR("can't evaluate command in F1 key");
        }
    } else if(key == GLUT_KEY_F2) {
        if(TCL_OK != gomp_SendCommand2Parser("eval $F2")){
            gomp_PrintERROR("can't evaluate command in F2 key");
        }
    }  else if(key == GLUT_KEY_F3) {
        if(TCL_OK != gomp_SendCommand2Parser("eval $F3")){
            gomp_PrintERROR("can't evaluate command in F3 key");
        }
    } else if(key == GLUT_KEY_F4) {
        if(TCL_OK != gomp_SendCommand2Parser("eval $F4")){
            gomp_PrintERROR("can't evaluate command in F4 key");
        }
    } else if(key == GLUT_KEY_F5) {
        if(TCL_OK != gomp_SendCommand2Parser("eval $F5")){
            gomp_PrintERROR("can't evaluate command in F5 key");
        }
    } else if(key == GLUT_KEY_F6) {
        if(TCL_OK != gomp_SendCommand2Parser("eval $F6")){
            gomp_PrintERROR("can't evaluate command in F6 key");
        }
    } else if(key == GLUT_KEY_F7) {
        if(TCL_OK != gomp_SendCommand2Parser("eval $F7")){
            gomp_PrintERROR("can't evaluate command in F7 key");
        }
    } else if(key == GLUT_KEY_F8) {
        if(TCL_OK != gomp_SendCommand2Parser("eval $F8")){
            gomp_PrintERROR("can't evaluate command in F8 key");
        }
    } else if(key == GLUT_KEY_F9) {
        if(TCL_OK != gomp_SendCommand2Parser("eval $F9")){
            gomp_PrintERROR("can't evaluate command in F9 key");
        }
    } else if(key == GLUT_KEY_F10) {
        if(TCL_OK != gomp_SendCommand2Parser("eval $F10")){
            gomp_PrintERROR("can't evaluate command in F10 key");
        }
    } else if(key == GLUT_KEY_F11) {
        if(TCL_OK != gomp_SendCommand2Parser("eval $F11")){
            gomp_PrintERROR("can't evaluate command in F11 key");
        }
    } else if(key == GLUT_KEY_F12) {
        if(TCL_OK != gomp_SendCommand2Parser("eval $F12")){
            gomp_PrintERROR("can't evaluate command in F12 key");
        }
/* arrow up/down keys */
    } else if(key == GLUT_KEY_UP) {
        if(TCL_OK != gomp_SendCommand2Parser("trajectory loop forward")){
            gomp_PrintERROR("can't execute 'trajectory loop forward'");
        }
    } else if(key == GLUT_KEY_DOWN) {
        if(TCL_OK != gomp_SendCommand2Parser("trajectory loop backward")){
            gomp_PrintERROR("can't execute 'trajectory loop backward'");
        }
/* control main Tcl/Tk widget up/down states */
    } else if(key == GLUT_KEY_PAGE_UP) {
        if(TCL_OK != gomp_SendCommand2Parser("wm deiconify .")){
            gomp_PrintERROR("can't evaluate command 'window 0 deiconify'");
        }
    } else if(key == GLUT_KEY_PAGE_DOWN) {
        if(TCL_OK != gomp_SendCommand2Parser("wm iconify .")){
            gomp_PrintERROR("can't evaluate command 'window 0 iconify'");
        }
    } else if(key == GLUT_KEY_END) {
        (void)gomp_DeletePickedAtomStack();
        (void)gomp_DrawScene();
    }


}
/****************************************************************************/
void NULLFunction()
/****************************************************************************/
{
    return;
}
/*
  If State != 0 window update will be done automatically
  == 0 window update has to be done through "display" command
*/
/****************************************************************************/
int gomp_GetUpdateDisplayMode()
/****************************************************************************/
{
    return(WindowUpdateDisplayMode.State);
}
/****************************************************************************/
int gomp_SetUpdateDisplayMode(int State)
/****************************************************************************/
{
    WindowUpdateDisplayMode.State = State;

    if(State)
        glutDisplayFunc(gomp_DrawSceneCallback);
    else
        glutDisplayFunc(NULLFunction);

    return(0);
}
/****************************************************************************/
int gomp_GetWindowParameters(int *Win_X, int *Win_Y, int *Win_W, int *Win_H)
/****************************************************************************/
{
    *Win_X = glutGet(GLUT_WINDOW_X);
    *Win_Y = glutGet(GLUT_WINDOW_Y);
    *Win_W = glutGet(GLUT_WINDOW_WIDTH);
    *Win_H = glutGet(GLUT_WINDOW_HEIGHT);

    return(0);
}
#endif /* GLUT */
#endif /* ENABLE_GRAPHICS */

/****************************************************************************/
int gomp_SetSystemRedisplayMode(int State)
/****************************************************************************/
{
    SystemRedisplayMode = State;

    return(0);
}
/****************************************************************************/
int gomp_GetSystemRedisplayMode()
/****************************************************************************/
{

    return(SystemRedisplayMode);
}
/****************************************************************************/
int gomp_GetAllowIndividualScaling()
/****************************************************************************/
{
    return(AllowIndividualScaling);
}
/****************************************************************************/
int gomp_SetAllowIndividualScaling(int Value)
/****************************************************************************/
{
    AllowIndividualScaling = Value;

    return(0);
}

/****************************************************************************/
int  gomp_PushPickedAtomToList(int Structure , int Atom)
/****************************************************************************/
{
    if( PickedAtomListLength ) {

/* check to se that it's not already in the list */
        if( IsAtomPicked(Structure , Atom) )
            return(0);

        PickedAtomList.Set     =
            gomp_ReallocateIntVector(PickedAtomList.Set , sizeof(int) *
                             (PickedAtomListLength + 1));
        PickedAtomList.HitAtom =
            gomp_ReallocateIntVector(PickedAtomList.HitAtom , sizeof(int) *
                             (PickedAtomListLength + 1));
        PickedAtomList.Set[PickedAtomListLength]     = Structure;
        PickedAtomList.HitAtom[PickedAtomListLength] = Atom;
        PickedAtomListLength++;
    } else {
        PickedAtomList.Set     = gomp_AllocateIntVector(sizeof(int));
        PickedAtomList.HitAtom = gomp_AllocateIntVector(sizeof(int));
        PickedAtomList.Set[0]     = Structure;
        PickedAtomList.HitAtom[0] = Atom;
        PickedAtomListLength   = 1;
    }

    if( !PickedAtomList.CallbackHandle )
        PickedAtomList.CallbackHandle = gomp_RegisterPlotter(
            PlotPickedAtoms,NULL,
            PLOTTER_NAME_ATOMPICKING,PLOTTER_ORDER_ATOMPICKING);

    return(0);
}

/****************************************************************************/
int  IsAtomPicked(int Structure , int Atom)
/****************************************************************************/
{
    static int i;

    for(i = 0 ; i < PickedAtomListLength ; i++) {
        if((PickedAtomList.Set[i]     == Structure) &&
           (PickedAtomList.HitAtom[i] == Atom))
            return i+1;
    }

    return(0);
}

/****************************************************************************/
int  gomp_PopPickedAtomFromList(int Structure , int Atom)
/****************************************************************************/
{
    static int i;

    i = IsAtomPicked(Structure, Atom)-1;

    if( i < 0 ) return(1);

    /* remove atom from the list */
    if( PickedAtomListLength < 1 )
        /* no need for the list anymore */
        gomp_DeletePickedAtomStack();
    else {
        if( PickedAtomListLength-(i+1) > 0 ) {
            memmove(&PickedAtomList.Set[i],
                    &PickedAtomList.Set[i+1],
                    (PickedAtomListLength-(i+1))*
                    sizeof(PickedAtomList.Set[0]));
            memmove(&PickedAtomList.HitAtom[i],
                    &PickedAtomList.HitAtom[i+1],
                    (PickedAtomListLength-(i+1))*
                    sizeof(PickedAtomList.HitAtom[0]));
        }
        --PickedAtomListLength;
    }
    return(0);
}

/****************************************************************************/
int  gomp_GetPickedAtomListLength()
/****************************************************************************/
{
    return(PickedAtomListLength);
}
/****************************************************************************/
const int *GetPickedAtomListFromStack()
/****************************************************************************/
{
    return(PickedAtomList.HitAtom);
}
/****************************************************************************/
const int *GetPickedStructureListFromStack()
/****************************************************************************/
{
    return(PickedAtomList.Set);
}
/****************************************************************************/
int FillPickedAtomSegmentResidueAtomList(SegmentResidueAtomList_t *List,
                                            int Which, int listType)
/****************************************************************************/
{
    static int  i;
    static const int *AtomList,*StructureList;

    StructureList = GetPickedStructureListFromStack();
    AtomList      = GetPickedAtomListFromStack();

    if( gomp_InitSegmentResidueAtomList(List,Which,listType,0) )
        return(1);

    for( i=0; i<gomp_GetPickedAtomListLength(); i++ ) {

        if( StructureList[i] >= List->from &&
            StructureList[i] <= List->to   &&
            gomp_PreSegmentResidueAtomListStructure(List,StructureList[i])==0 ) {

            (void)gomp_PostSegmentResidueAtomListStructure(
                List,StructureList[i],1,&AtomList[i]);
        }
    }

    return(0);
}

/****************************************************************************/
int gomp_ListPickedAtoms(int Which, int listType)
/****************************************************************************/
{
    static SegmentResidueAtomList_t List;

    if( gomp_GetPickedAtomListLength() <= 0 ) {
        gomp_SendTclReturn("");
        return(0);
    }

    if(FillPickedAtomSegmentResidueAtomList(&List,Which,listType))
        return(1);

    return(gomp_ReturnAndFreeSegmentResidueAtomList(&List));
}
/****************************************************************************/
int  gomp_DeletePickedAtomStack()
/****************************************************************************/
{
    if(PickedAtomListLength) {
        free(PickedAtomList.Set);
        free(PickedAtomList.HitAtom);
        PickedAtomListLength = 0;
        
        gomp_UnregisterPlotter(
            PickedAtomList.CallbackHandle);
        PickedAtomList.CallbackHandle = NULL;
    }

    return(0);
}

#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int PlotPickedAtom(int Structure, int Atom, double radius)
/****************************************************************************/
{
    static int   SphereQ;
    static int   SphereQ2;
    static const float *ColRED,*ColGREEN,*ColBLUE;
    static float CRh;
    static float CGh;
    static float CBh;
    const float *RotMP;

    /* plot a sphere */
    glDisable(GL_LIGHTING);
    SphereQ  = gomp_GetSphereQuality();
    SphereQ2 = 2 * SphereQ;
 
    glPushMatrix();
    glLoadIdentity();
    if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {
    } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {
        gluLookAt(0.0 , 0.0 , 
                  (gomp_GetPerspectiveNear() + 0.5 * gomp_GetPerspectiveWindow()) ,
                  0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0);
    } else {
        glPopMatrix();
        gomp_PrintERROR("undefined Projection Transformation defined");
        return(1);
    }

    RotMP = gomp_GetSavedModelViewMatrixMT(Structure);
    glMultMatrixf(RotMP);
    glTranslatef(
        gomp_GetAtomXCoord(Structure,Atom),
        gomp_GetAtomYCoord(Structure,Atom),
        gomp_GetAtomZCoord(Structure,Atom) );

    ColRED   = gomp_GetAtomColourRedPointer(Structure);
    ColGREEN = gomp_GetAtomColourGreenPointer(Structure);
    ColBLUE  = gomp_GetAtomColourBluePointer(Structure);

    CRh = ColRED[Atom];
    CGh = ColGREEN[Atom];
    CBh = ColBLUE[Atom];

    if(!gomp_GetDisplayColourType())
        (void)gomp_RGB2Grayscale(&CRh , &CGh , &CBh);
        
    glColor4f(CRh , CGh , CBh , 1.0);

    gluSphere(gomp_SphereQuad, radius , SphereQ2 , SphereQ);
    glPopMatrix();
    glEnable(GL_LIGHTING);

    return(0);
}
#endif /* ENABLE_GRAPHICS */

/****************************************************************************/
int PlotPickedAtoms(void* userData,int Wstr,int drawFlags)
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int   i,sl;
    static int   HitAtom;
    static const float *ColRED,*ColGREEN,*ColBLUE;
    static const float *XC,*YC,*ZC;
    static int   SphereQ;
    static int   SphereQ2;
    static const int *AtomListP;
    static const int *StructureListP;
    static float CRh;
    static float CGh;
    static float CBh;

    if ( ! ( drawFlags & gom_PlotSimpleElements ) )
        return(-1);

    AtomListP       = GetPickedAtomListFromStack();
    StructureListP  = GetPickedStructureListFromStack();

    glDisable(GL_LIGHTING);

/* loop over the structures */
    for ( i = 0 ; i < gomp_GetPickedAtomListLength() ; i++ ) {

        if ( StructureListP[i] != Wstr )
            continue;

        HitAtom = AtomListP[i];

        ColRED    = gomp_GetAtomColourRedPointer(sl);
        ColGREEN = gomp_GetAtomColourGreenPointer(sl);
        ColBLUE = gomp_GetAtomColourBluePointer(sl);
        XC        = gomp_GetAtomXCoordPointer(sl);
        YC       = gomp_GetAtomYCoordPointer(sl);
        ZC      = gomp_GetAtomZCoordPointer(sl);
         
/* plot a sphere */

        SphereQ  = gomp_GetSphereQuality();
        SphereQ2 = 2 * SphereQ;

        glPushMatrix();

        glTranslatef(XC[HitAtom],YC[HitAtom],ZC[HitAtom]);

        CRh = ColRED[HitAtom];
        CGh = ColGREEN[HitAtom];
        CBh = ColBLUE[HitAtom];

        if(!gomp_GetDisplayColourType())
            (void)gomp_RGB2Grayscale(&CRh , &CGh , &CBh);
        
        glColor4f(CRh , CGh , CBh , 1.0);
        /*   gluSphere(SphereQuad,(double)PlotSphere.Radius[i], SphereQ2 , SphereQ);*/
        gluSphere(gomp_SphereQuad, 0.3 , SphereQ2 , SphereQ);
        glPopMatrix();
    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}

/****************************************************************************/
int  gomp_PlotSelectedAtoms(void* userData,int Wstr,int drawFlags)
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    const float *ColRED,*ColGREEN,*ColBLUE;
    const float *XC,*YC,*ZC;
    int   SphereQ;
    int   SphereQ2;
    float CRh;
    float CGh;
    float CBh;
    int NAtoms, i, iAtom;
    const int *sel_list;

    if ( ! ( drawFlags & gom_PlotSimpleElements ) ||
         Wstr >= gomp_GetNumMolecStructs() )
        return(-1);

    glDisable(GL_LIGHTING);

    sel_list = gomp_GetSelectionList(Wstr);

    /* get atom information */
    ColRED    = gomp_GetAtomColourRedPointer(Wstr);
    ColGREEN = gomp_GetAtomColourGreenPointer(Wstr);
    ColBLUE = gomp_GetAtomColourBluePointer(Wstr);
    XC        = gomp_GetAtomXCoordPointer(Wstr);
    YC       = gomp_GetAtomYCoordPointer(Wstr);
    ZC      = gomp_GetAtomZCoordPointer(Wstr);

    NAtoms    = gomp_GetSelectionListLength(Wstr);

    /* loop over selected atoms */
    for( i=0; i<NAtoms; i++ ) {
        /* plot an octahedral */
        iAtom = sel_list[i];
            
        SphereQ  = gomp_GetSphereQuality();
        SphereQ2 = 2 * SphereQ;
        
        glPushMatrix();
        glTranslatef(XC[iAtom],YC[iAtom],ZC[iAtom]);

        CRh = ColRED[iAtom];
        CGh = ColGREEN[iAtom];
        CBh = ColBLUE[iAtom];

        if(!gomp_GetDisplayColourType())
            (void)gomp_RGB2Grayscale(&CRh , &CGh , &CBh);
            
        glColor4f(CRh , CGh , CBh , 1.0);
        gluSphere(gomp_SphereQuad, 0.5 , 4, 2);
        glPopMatrix();
    }
#endif /* ENABLE_GRAPHICS */
    
    return(0);
}

#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int DrawSceneQuadStereo ()
/****************************************************************************/
{


/* ......... */
    (void)gomp_PushModelViewingData();
/* ......... */

    glLineWidth((GLfloat)(gomp_GetMoleculeLineWidth()+0.0));

/* Picture #1 */    
    glDrawBuffer(GL_BACK_LEFT);
    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);

    gomp_Rotate(  +gomp_GetQuadStereoHalfAngle()     , 0.0 , 1.0 , 0.0);
    
    if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {
    } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {
        gluLookAt(0.0 , 0.0 , 
                  (gomp_GetPerspectiveNear() + 0.5 * gomp_GetPerspectiveWindow()) ,
                  0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0);
    } else {
        gomp_PrintERROR("undefined Projection Transformation defined");
        return(1);
    }

/* ......... */

    (void)gomp_UpdateScreen();

/* ......... */
    (void)gomp_PopModelViewingData();
/* ......... */

/* Picture #2 */    
    glDrawBuffer(GL_BACK_RIGHT);
    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);

    gomp_Rotate(  -gomp_GetQuadStereoHalfAngle()     , 0.0 , 1.0 , 0.0 );
    
    glMatrixMode(GL_MODELVIEW);

    if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {
    } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {
        gluLookAt(0.0 , 0.0 , 
                  (gomp_GetPerspectiveNear() + 0.5 * gomp_GetPerspectiveWindow()) ,
                  0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0);
    } else {
        gomp_PrintERROR("undefined Projection Transformation defined");
        return(1);
    }

/* ......... */


    (void)gomp_UpdateScreen();

/* ......... */
    (void)gomp_PopModelViewingData();
/* ......... */

    glDrawBuffer(GL_BACK); /* this is default for double buffered views */

    return(0);

}
/****************************************************************************/
int gomp_CheckHardwareStereo()
/****************************************************************************/
{
    static GLboolean Stereo;

    glGetBooleanv(GL_STEREO , & Stereo);

    return((int)Stereo);
}
#endif /* ENABLE_GRAPHICS */

/*********************************************************************/
int   gomp_SetGlobalTextString(const char *InString)
/*********************************************************************/
{
    gomp_CopyString(GlobalTextString,InString,BUFF_LEN);

    return(0);
}
/*********************************************************************/
const char *gomp_GetGlobalTextString()
/*********************************************************************/
{
    return(GlobalTextString);
}

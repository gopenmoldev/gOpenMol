/*
  Copyright (c) 1994 - 2004 by:
  Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
  Confidential unpublished property of 
  Leif Laaksonen  
  All rights reserved
  


  This will eventually be a graphical interface into the 
  OpenMol system of software.

  The system will be built on Tcl/Tk 

  Work started summer 1992,
  continuing on summer 1993, 1994, 1995, 1996

  Leif Laaksonen
  Center for Scientific Computing / MPI Astrophysics

  Enhancements 2002 - 2003 by:
  Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#if HAVE_UNISTD_H
#include <unistd.h>
#endif
#if defined(WIN32)
#include <windows.h>
#endif

#include "coord_file.h"
#include "gomlib/gomdefs.h"
#include "gomenv.h"
#include "gomlog.h"
#include "gommain.h"
#include "gomproc.h"
#include "gomstring.h"
#include "printmsg.h"
#include "projview.h"
#include "stereo.h"
#include "tclutils.h"

#include "stdafx.h"

/* SIGISMONDO: flag -s */
static int    StereoDisplayState = STEREO_DISPLAY_OFF;
 
/* END SIGISMONDO */

/* define some global variables */
static int TermType = GRAPHICS_STATUS_UNKNOWN;
/*   0 means ==>     no graphics
     1               X-windows/OpenGL/WIN32
     2               status is not yet determined
*/

static int DebugLevel = 0;          /*   debug level
                                         = 0 , no debug on
                                         = 1 , gives some extra output
                                         > 1 , gives a massive amount of data
                                */ 
static struct {
    int Width;
    int Height;
    int Xc;
    int Yc;
} WindowInfo = { 500 , 500 , 1 , 1};

static char StartUpDir[BUFF_LEN];

static int GetInputParams(int , const char *[]);
static int GetStartUpWindow(const char *);

/* END of global variables      */

/* possibility to create also gOpenMol as a library */

#if defined(GOMLIB)
extern int gom_main(int, const char *[]);
int gom_main(int argc, const char *argv[])
#else
/*****************************************************************/
int main(int argc, const char *argv[])
/*****************************************************************/
#endif
{
/* Read parameters from the input line                    */
    if(argc > 1) 
        GetInputParams(argc,argv);

/* initialize tcl */

    if(gomp_Tcl_AppInit(argc , argv)) {
        printf("Tcl/Tk ERROR, can't continue\n");
        Tcl_Exit(1);
    }

/* handle the environment definitions */
    (void) gomp_GetEnv(argc , argv);

/*   Tcl_FindExecutable(argv[0]); this used to work for version before 8.1b2 */

/* open log file (if needed) */
    if(gomp_OpenLogFile()) {
        gomp_PrintERROR("can't open log file");
    }
    (void) gomp_WriteToLogFile("OpenMol Starting ...");

/* here we go ...                                         */
    /* (void)gomp_PrintMessage(OPENMOL_label);*/

/* Tell them who we are ...                               */
    (void)gomp_PrintBANNER();

/* Be prepared for the log file                           */
    (void)gomp_SaveModuleName(argv[0]);
    (void)gomp_Get_proc_info();

/* initialize the cputime counter */
    gomp_Get_cpu_secs(NULL, NULL);
/* .............................. */

/* Go and check the type of terminal device available ... */
    if(TermType == GRAPHICS_STATUS_UNKNOWN) {
#ifdef ENABLE_GRAPHICS
        if(gomp_Mmain(argc,argv))  /* can't open X-Window display */
            TermType = GRAPHICS_NOT_AVAILABLE;
#else
        TermType = GRAPHICS_NOT_AVAILABLE;
#endif
    }

/* read the startup tcl script */
    (void)gomp_OpenStartFile();

    if(argc > 1) 
        gomp_OpenInputFiles(argc,argv);

    (void)gomp_CommandLineInput();

    gomp_Exit(0);
    return(0); /* Will never be reached. */
}
/*****************************************************************/
int GetInputParams(int argc, const char *argv[])
/*
  int   argc;
  const char *argv[];
*/
/*****************************************************************/
{
    int i;

/* We need this then reading coordinate files. */
#ifdef WIN32
    if ( ! GetCurrentDirectory(BUFF_LEN,StartUpDir) )
#else
    if ( getcwd(StartUpDir,BUFF_LEN) == NULL )
#endif
        strcpy(StartUpDir,".");

    for(i = 1 ; i < argc ; i++) {

        if(argv[i][0] == '-')     {

            switch(argv[i][1])       {

/* -c input coordinate file */
            case 'c':          /* read a coordinate file */
                if ( argv[i][2] )
                    break;
                if ( (i+1) >= argc || !argv[i+1][0] )
                    gomp_PrintEXIT("file name is missing");
                i++;
                break;
/* -d debug level           */
            case 'd':          /* define debug level     */
                if(argv[i][2] != (char)NULL) {
                    DebugLevel = atoi(&argv[i][2]);
                } else {
                    if((i+1) >= argc) {
                        gomp_PrintERROR("debug level is missing");
                        exit(1);
                    } else {
                        i++;
                        DebugLevel = atoi(&argv[i][0]);
                    }
                }
                break;
/* -e environment file input */
            case 'e':
                if ( argv[i][2] ) {
                    (void)gomp_SetEnvironmentFile(&argv[i][2]);
                } else {
                    if ( (i+1) >= argc || !argv[i+1][0] )
                        gomp_PrintEXIT("file name is missing");
                    else
                        (void)gomp_SetEnvironmentFile(argv[++i]);
                }
                break;               
/* -h print help ...        */
            case 'h':
                gomp_PrintUSAGE();
                exit(0);
/* -l log file control      */
            case 'l':
                switch ( argv[i][2] ) {
                case 'n':
                case 'N':
                    (void)gomp_SetLogFileStatus(1);
                    break;
                case 'y':
                case 'Y':
                    (void)gomp_SetLogFileStatus(0); 
                    break;
                case '\0':
                    if ( (i+1) >= argc || !argv[i+1][0] )
                        gomp_PrintEXIT("parameter is missing");
                    else {
                        switch ( argv[++i][0] ) {
                        case 'n':
                        case 'N':
                            (void)gomp_SetLogFileStatus(1);
                            break;
                        case 'y':
                        case 'Y':
                            (void)gomp_SetLogFileStatus(0);
                            break;
                        }
                    }
                }
                break;
/* -m start with model file */
            case 'm':         /* read a model file */
                if ( argv[i][2] )
                    break;
                if( (i+1) >= argc || !argv[i+1][0] )
                    gomp_PrintEXIT("file name is missing");
                i++;
                break;
/* -p start with parameter file */
            case 'p':
                if ( argv[i][2] )
                    (void)gomp_SetAtomParametersFile(&argv[i][2]);
                else {
                    if( (i+1) >= argc || !argv[i+1][0] )
                        gomp_PrintEXIT("file name is missing");
                    else
                        (void)gomp_SetAtomParametersFile(argv[++i]);
                }
                break;
/* -x define server port number */
            case 'x':          /* define new port number */
/*
  ServerPort = atoi(&argv[i][2]);
*/
                break;
/* -t use no graphics           */
            case 't':          /* no graphics */
                TermType = 0;
                break;
/* -w startup window            */
            case 'w':          /* startup window parameters */
                if ( argv[i][2] ) {
                    (void)GetStartUpWindow(&argv[i][2]);
                } else {
                    if( (i+1) >= argc || !argv[i+1][0] )
                        gomp_PrintEXIT("file name is missing");
                    else
                        (void)GetStartUpWindow(argv[++i]);
                }
                break;
/* -s STEREO             */
            case 's':          /* stereo */
                gomp_SetStereoDisplayState(STEREO_DISPLAY_ON);
                break;
            default:
                gomp_FormatMessage("?ERROR - unknow descriptor '%s'\n",
                                   argv[i]);
            }
        }
    }

    return(0);
}
/*****************************************************************/
int gomp_OpenInputFiles(int argc, const char *argv[])
/*
  int   argc;
  const char *argv[];
*/
/*****************************************************************/
{
    int i;
    char File[BUFF_LEN];
    const char *name;
    int (*Read)(const char*);

    for(i = 1 ; i < argc ; i++) {

        if(argv[i][0] == '-') {

            switch ( argv[i][1] ) {

/* -c input coordinate file */
            case 'c':          /* read a coordinate file */
                Read = gomp_StartWithCoordinateFile;
                break;
/* -m start with model file */
            case 'm':         /* read a model file */
                Read = gomp_StartWithModelFile;
                break;
            default:
                continue;
            }
            if ( argv[i][2] != '\0' )
                name = &argv[i][2];
            else
                name = argv[++i];
        }
        else {
/* look for both input coordinate name and tcl script */
            Read = gomp_StartWithCoordinateFile;
            name = argv[i];
        }

        if ( ! strstr(name,"://") ) {
/* not an URL */
#ifdef WIN32
            if ( name[0] &&
                 strncmp(name+1,":\\",2) == 0 &&
                 strncmp(name+1,":/" ,2) == 0 ) {
                /* relative WIN32 path */
                snprintf(File,BUFF_LEN,"%s\\%s",StartUpDir,name);
                name = File;
            }
#else
            if ( name[0] != '/' ) {
                /* relative UNIX path */
                snprintf(File,BUFF_LEN,"%s/%s",StartUpDir,name);
                name = File;
            }
#endif
        }

        if ( Read(name) )
            gomp_Exit(1);
    }

    return(0);
}
/*****************************************************************/
int gomp_GetTermType()
/*****************************************************************/
{
    return(TermType);
}
/*****************************************************************/
int gomp_SetTermType(int State)
/*****************************************************************/
{
    TermType = State;
    return(0);
}

/***********************************************************************/
int GetStartUpWindow(const char *String)
/***********************************************************************/
{

    int j;
    char  IString[BUFF_LEN];
    const char *ISparsed[BUFF_LEN];


    j = gomp_SplitString(IString , "x" , ISparsed );

    if(j < 1) {
        gomp_PrintERROR("please supply all parameters (-wXCxYCxWidth[xHeight])");
        return(1);
    }

    WindowInfo.Xc = atoi(ISparsed[0]);
    WindowInfo.Yc = atoi(ISparsed[1]);

    WindowInfo.Width = atoi(ISparsed[2]);
    WindowInfo.Height = atoi(ISparsed[3]);
       
    if(WindowInfo.Height < 1) WindowInfo.Height = WindowInfo.Width;

    if(WindowInfo.Xc < 0 || WindowInfo.Yc < 0) {
        gomp_PrintERROR("in window coordinates (out of window)");
        return(1);
    }

    if(WindowInfo.Width < 1 || WindowInfo.Height < 1) {
        gomp_PrintERROR("window is too small");
        return(1);
    }

/* ok return ... */

    return(0);
}
/*****************************************************************/
int gomp_GetStereoDisplayState()
/*****************************************************************/
{
    return(StereoDisplayState);
}
/*****************************************************************/
int gomp_SetStereoDisplayState(int State)
/*****************************************************************/
{
    StereoDisplayState = State;

    return(0);
}



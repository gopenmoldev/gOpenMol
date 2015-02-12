/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2004 by:
Eero Häkkinen


This code is based on the tcl-command language
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

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
#endif /* ENABLE_GRAPHICS */

#include <tcl.h>

#include "gomfile.h"
#include "gomlog.h"
#include "gommain.h"
#include "gomproc.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "listutils.h"
#include "molecstruct.h"
#include "picking.h"
#include "plot.h"
#include "printmsg.h"
#include "parser.h"
#include "tclutils.h"

#include "stdafx.h"

#define GRAPHICS_YES 1
#define GRAPHICS_NO  0

#define PARSER 1

/* Tcl - interface */
#define      TCL_RUN_DISPLAY_FILE_VAR  "tcl_file"  
#if 0
static int         SetInputTclScriptFileName(const char *);
#endif
static const char *GetInputTclScriptFileName(void);

static Tcl_DString Tcl_Cmd;
static Tcl_Interp *Interp;

static char         gomInputScriptName[BUFF_LEN];

/* commands */
static int CatchCommand(ClientData , Tcl_Interp *, int , Tcl_Obj * const[]);
static int Quit(ClientData , Tcl_Interp *, int , const char **);
static int PauseCommand(ClientData , Tcl_Interp *, int , const char **);
static int Quit(ClientData , Tcl_Interp *, int , const char **);
static int FilterCommand(ClientData , Tcl_Interp *, int , const char **);
static int WebURLCommand(ClientData , Tcl_Interp *, int , const char **);
static int PrintInterface(ClientData , Tcl_Interp *, int , const char **);
static int ErrorInterface(ClientData , Tcl_Interp *, int , const char **);
static int PipeCommand(ClientData , Tcl_Interp *, int , const char **);
static int TestCommand(ClientData , Tcl_Interp *, int , const char **);

/*********************************************************************/
int  gomp_Tcl_AppInit(int argc , const char *argv[])
/*********************************************************************/
{

    int code;

    Interp = Tcl_CreateInterp();

    if(Interp == (Tcl_Interp *) NULL) {
        printf("ERROR: can't create interpreter\n");
        return TCL_ERROR;
    }

    Tcl_FindExecutable(argv[0]);

#if defined(USE_TCL_STUBS)
    if (Tcl_InitStubs(Interp , "8.4" , 0) == NULL)
#else
        if (Tcl_Init(Interp) != TCL_OK)
#endif
        {   
            printf("ERROR: can't initialize the Tcl shell\n");
            return TCL_ERROR;
        }

/* define the commands */
/* catch */
// Catch is already defined in Tcl, collision...
//    Tcl_CreateObjCommand(Interp,"catch",CatchCommand,(ClientData)NULL,
//                         (Tcl_CmdDeleteProc *)NULL);
/* quit */
    Tcl_CreateCommand(Interp,"quit",Quit,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);
/* get */
    Tcl_CreateCommand(Interp,"import",gomp_ImportCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);
/* show */
    gomp_CreateGomParser(Interp,"show",gomp_ShowCommand);
/* trajectory */
    Tcl_CreateCommand(Interp,"trajectory",
                      gomp_TrajectoryCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* rotate */
    Tcl_CreateCommand(Interp,"rotate",
                      gomp_RotateCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* atom */
    gomp_CreateGomParser(Interp,"atom",gomp_AtomCommand);

/* gomSyntax */
    gomp_CreateGomParser(Interp,"gomSyntax",gomp_SyntaxCommand);

/* manipulate */
    Tcl_CreateCommand(Interp,"manipulate",
                      gomp_ManipulateCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* monitor */
    Tcl_CreateCommand(Interp,"monitor",
                      gomp_MonitorCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* plot */
    Tcl_CreateCommand(Interp,"plot",
                      gomp_PlotCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* calculate */
    Tcl_CreateCommand(Interp,"calculate",
                      gomp_CalculateCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* export */
    Tcl_CreateCommand(Interp,"export",
                      gomp_ExportCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* define */
    Tcl_CreateCommand(Interp,"define",
                      gomp_DefineCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* contour */
    Tcl_CreateCommand(Interp,"contour",
                      gomp_ContourCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* scale */
    Tcl_CreateCommand(Interp,"gscale",
                      gomp_ScaleCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* select */
    Tcl_CreateCommand(Interp,"select",
                      gomp_SelectCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* translate */
    Tcl_CreateCommand(Interp,"translate",
                      gomp_TranslateCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* edit */
    Tcl_CreateCommand(Interp,"edit",
                      gomp_EditCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* reset */
    Tcl_CreateCommand(Interp,"reset",
                      gomp_ResetCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* fill */
    Tcl_CreateCommand(Interp,"fill",
                      gomp_FillCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* help */
    Tcl_CreateCommand(Interp,"help",
                      gomp_HelpCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* mopen */
    Tcl_CreateCommand(Interp,"mopen",
                      gomp_MOpenCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* msave */
    Tcl_CreateCommand(Interp,"msave",
                      gomp_MSaveCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* plumber */
    Tcl_CreateCommand(Interp,"plumber",
                      gomp_PlumberCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* trace */
    Tcl_CreateCommand(Interp,"ptrace",
                      gomp_TraceCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* center  */
    Tcl_CreateCommand(Interp,"center",
                      gomp_CenterCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* diagonalize  */
    Tcl_CreateCommand(Interp,"diagonalize",
                      gomp_DiagonalizeCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);
#ifdef ENABLE_GRAPHICS
/* display */
    Tcl_CreateCommand(Interp,"display",gomp_DisplayCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* window      */
    Tcl_CreateCommand(Interp,"window",
                      gomp_WindowCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);
/* copy  */
    Tcl_CreateCommand(Interp,"copy",
                      gomp_CopyCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* hardcopy */
    Tcl_CreateCommand(Interp,"hardcopy",
                      gomp_HardcopyCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);
#endif /* ENABLE_GRAPHICS */

/* find      */
    Tcl_CreateCommand(Interp,"find",
                      gomp_FindCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* pause      */
    Tcl_CreateCommand(Interp,"pause",
                      PauseCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* refresh    */
    gomp_CreateGomParser(Interp,"refresh",gomp_RefreshCommand);

/* filter can be used instead of 'eval' */
    Tcl_CreateCommand(Interp,"filter",
                      FilterCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* WebURL */
    Tcl_CreateCommand(Interp,"webbrowse",
                      WebURLCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* PrintInterface */
    Tcl_CreateCommand(Interp,"gomPrint",
                      PrintInterface,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* ErrorInterface */
    Tcl_CreateCommand(Interp,"gomError",
                      ErrorInterface,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* Picking Interface */
    Tcl_CreateCommand(Interp,"picking",
                      gomp_PickingCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* Pipe Interface */
    Tcl_CreateCommand(Interp,"pipe",
                      PipeCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);
/* Python Interface */
    Tcl_CreateCommand(Interp,"python",
                      gomp_PythonCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* gomTest Interface */
    Tcl_CreateCommand(Interp,"gomtest",
                      TestCommand,(ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

/* some set commands before we start */
    code = Tcl_GlobalEval(Interp,"set tcl_interactive 1");
    code = Tcl_GlobalEval(Interp,"set auto_noexec 1");

    return(TCL_OK);
}
/*********************************************************************/
int  gomp_SendCommand2Parser(const char *Command)
/*********************************************************************/
{
    static int  code, length;
    static int  Switch = 0;

#if defined(JUNK)
    if(gomp_GetPythonState()) {

        return(PyRun_SimpleString(Command));

    }
#endif

    if(Switch==0) {
        Tcl_DStringInit(&Tcl_Cmd);
        Switch = 1;
    }

    Tcl_DStringAppend(&Tcl_Cmd , Command, - 1);

    if(!Tcl_CommandComplete(Tcl_DStringValue(&Tcl_Cmd)))
        return(TCL_OK);

    code = Tcl_GlobalEval(Interp,Tcl_DStringValue(&Tcl_Cmd));
    if ( *Tcl_GetStringResult(Interp) )
        gomp_PrintMessage(Tcl_GetStringResult(Interp));

    Switch = 0;
    length = Tcl_DStringLength(&Tcl_Cmd);
    printf(
        (length && (Tcl_DStringValue(&Tcl_Cmd))[length-1] == '\n') ?
        ">>> %s" : ">>> %s\n",
        Tcl_DStringValue(&Tcl_Cmd));
    Tcl_DStringFree(&Tcl_Cmd);
    return(code);
}
/*********************************************************************/
int  gomp_SendTclObj2Parser( Tcl_Obj *obj )
/*********************************************************************/
{
    int  code;

    code = Tcl_EvalObjEx(Interp,obj,TCL_EVAL_GLOBAL);
    if ( *Tcl_GetStringResult(Interp) )
        gomp_PrintMessage(Tcl_GetStringResult(Interp));

    printf(">>> %s\n",Tcl_GetString(obj));
    return(code);
}
/*********************************************************************/
int gomp_CommandLineInput()
/*********************************************************************/
{
    static char  CommandLine[BUFF_LEN];
    static const char *InputTclScript;

/* if there is a tcl script coming execute that first */
    InputTclScript = GetInputTclScriptFileName();
    if(InputTclScript[0] != (char)NULL) {
        sprintf(CommandLine,"source {%s}",InputTclScript);
        if(TCL_OK != gomp_SendCommand2Parser(CommandLine)) {
            gomp_PrintERROR("can't evaluate supplied input tcl script");
        }
    }

    while(1) {

        printf("> ");
        if ( fgets(CommandLine,BUFF_LEN,stdin) == NULL )
            exit(1);

        (void)gomp_SendCommand2Parser(CommandLine); 
    }

    return(0);
}
/*********************************************************************/
int CatchCommand(ClientData clientdata , Tcl_Interp *interp,
                 int objc, Tcl_Obj * const objv[])
/*********************************************************************/
{
    int      code,oldSuppress;
    Tcl_Obj *result;

    if ( objc < 2 || objc > 3 ) {
        Tcl_WrongNumArgs(interp,1,objv,"command ?varName?");
        return(TCL_ERROR);
    }

    /* Evaluate. Don't show errors. */
    oldSuppress = gomp_SuppressErrors(1);
    code        = Tcl_EvalObjEx(interp,objv[1],0);
    gomp_SuppressErrors(oldSuppress);
    if ( objc == 3 ) {
        /* store the result */
        result = Tcl_GetObjResult(interp);
        if ( Tcl_ObjSetVar2(interp,objv[2],NULL,result,
                            TCL_LEAVE_ERR_MSG) == NULL )
            return(TCL_ERROR);
    }
    /* return the error code */
    result = Tcl_NewIntObj(code);
    Tcl_SetObjResult(interp,result);

    return(TCL_OK);
}
/*********************************************************************/
int Quit(ClientData clientdata , Tcl_Interp *interp,
         int argc , const char *argv[])
/*********************************************************************/
{
    static const char *value;

    if(argc != 1) {
        gomp_PrintERROR("Wrong # arguments: (should be 1)");
        return(TCL_ERROR);
    }

/* save information into a startup file at exit */
    value  = Tcl_GetVar(interp , 
                        "lulProgramShutdownString" , 
                        TCL_GLOBAL_ONLY);

    if(value != (const char *)NULL) {
        if(TCL_OK != gomp_SendCommand2Parser(value)) {
            gomp_PrintERROR("can't evaluate supplied input tcl script at exit");
        }
    }

/* shot down the various parts of the program ... */
    (void)gomp_Get_date(2);
    (void)gomp_CloseLogFile();

/* and last shot down Tcl/Tk and exit ...        */
    Tcl_Exit(0);
/* newer reached                                 */
    return(-1);
}

/*********************************************************************/
int PauseCommand(ClientData clientdata , Tcl_Interp *interp,
                 int argc , const char *argv[])
/*********************************************************************/
{
    static char  Text[BUFF_LEN];
    static int   Itime;

    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

    if(Text[0] == (char)NULL) {
        gomp_PrintERROR("pause time missing");
        return(1);
    }

    Itime = (int)(atof(Text) * 1000.0);

/* sleep for milliseconds */
    Tcl_Sleep(Itime);

    return(0);
}


/*********************************************************************/
static int ParseRefresh(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
/*********************************************************************/
{
    GOM_PARSER_VERIFY( gomp_UpdateData() == 0 );
    GOM_PARSER_SUCCEEDED;
}
/*********************************************************************/
const gom_ParserArgumentList gomp_RefreshCommand[]={
    GOM_PARSER_FINAL_ARGS(GOM_PARSER_NO_MORE_ARGS,
                          ParseRefresh,
                          GOM_PARSER_UNUSED_VALUE),
    GOM_PARSER_END_ARGUMENT_LIST
};

/*********************************************************************/
int FilterCommand(ClientData clientdata , Tcl_Interp *interp,
                  int argc , const char *argv[])
/*********************************************************************/
{
    static char  Text[BUFF_LEN];
    static char  Text1[BUFF_LEN];
    Tcl_DString  gomp_Tcl_Collect;
    static int   i;

    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

    if(Text1[0] == (char)NULL) 
        return(1);

    Tcl_DStringInit(&gomp_Tcl_Collect);

    for(i = 1 ; i < argc ; i++) {
        Tcl_DStringAppend(&gomp_Tcl_Collect , "{"      , - 1);
        Tcl_DStringAppend(&gomp_Tcl_Collect , argv[i]  , - 1);
        Tcl_DStringAppend(&gomp_Tcl_Collect , "} "     , - 1);
    }

    if(gomp_SendCommand2Parser(gomp_Tcl_Collect.string)) {
        Tcl_DStringFree(&gomp_Tcl_Collect);
        return(TCL_ERROR);
    }

    Tcl_DStringFree(&gomp_Tcl_Collect);
    return(TCL_OK);

}
/*********************************************************************/
int gomp_PickingCommand(ClientData clientdata , Tcl_Interp *interp,
                        int argc , const char *argv[])
/*********************************************************************/
{
    static char  Text[BUFF_LEN];
    static char  Text1[BUFF_LEN];
    static char  Text2[BUFF_LEN];
    static int   ITemp1;
    static int   ITemp2;

    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

    if(Text[0] == (char)NULL) {
        gomp_PrintERROR("picking parameter missing");
        return(TCL_ERROR);
    }

    if(gomp_StringMatch(Text , "atom")) {
        gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
        if(Text1[0] == (char)NULL) {
            gomp_PrintERROR("atom index value missing");
            return(TCL_ERROR);
        }
        ITemp1 = atoi(Text1) - 1;
        gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
        if(Text2[0] == (char)NULL) {
            ITemp2 = 0;
        } 
        else {
            ITemp2 = atoi(Text2) - 1;
            if(ITemp2 < 0 || ITemp2 >= gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("structure index outside allowed range");
                return(TCL_ERROR);
            }
        }

        if(gomp_PushPickedAtomToList(ITemp2 , ITemp1))
            return(TCL_ERROR);
        else
            return(TCL_OK);
    }
    if(gomp_StringMatch(Text , "-atom")) {
        gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
        if(Text1[0] == (char)NULL) {
            gomp_PrintERROR("atom index value missing");
            return(TCL_ERROR);
        }
        ITemp1 = atoi(Text1) - 1;
        gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
        if(Text2[0] == (char)NULL) {
            ITemp2 = 0;
        } 
        else {
            ITemp2 = atoi(Text2) - 1;
            if(ITemp2 < 0 || ITemp2 >= gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("structure index outside allowed range");
                return(TCL_ERROR);
            }
        }

        if(gomp_PopPickedAtomFromList(ITemp2 , ITemp1))
            return(TCL_ERROR);
        else
            return(TCL_OK);
    }
    else if(gomp_StringMatch(Text , "get")) {
        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
        
        gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
        if(Text1[0] == (char)NULL) {
            gomp_PrintERROR("atom index value missing");
            return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text1 , "all"))
            ITemp1 = LIST_ALL_STRUCTURES;
        else if(gomp_StringMatch(Text1 , "sele$cted"))
            ITemp1 = LIST_SELECTED_STRUCTURES;
        else {
            ITemp1 = atoi(Text1) - 1;
            if(ITemp1 < 0 || ITemp1 >= gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("structure index outside allowed range");
                return(TCL_ERROR);
            }
        }

        if(gomp_StringMatch(Text , "atom$s"))
            ITemp2 = ATOM_LIST;
        else if(gomp_StringMatch(Text , "resi$due"))
            ITemp2 = RESIDUE_LIST;
        else if(gomp_StringMatch(Text , "segm$ent"))
            ITemp2 = SEGMENT_LIST;
        else if(gomp_StringMatch(Text , "sel$list"))
            ITemp2 = SELECTION_LIST;
        else {
            gomp_PrintERROR("wrong picking parameter suppied, allowed are 'atoms/residue/segment'");
            return(TCL_ERROR);
        }

        if(gomp_ListPickedAtoms(ITemp1 , ITemp2))
            return(TCL_ERROR);
        else
            return(TCL_OK);
    }
    else if(gomp_StringMatch(Text , "rese$t")) {
        if(gomp_DeletePickedAtomStack()) 
            return(TCL_ERROR);
        else {
            return(TCL_OK);
        }
    }
    else if(gomp_StringMatch(Text , "stat$us")) {

        sprintf(Text,"%d",gomp_GetPickedAtomListLength());
        (void)gomp_SendTclReturn(Text);
        return(TCL_OK);
    } 
    else {

        gomp_PrintERROR("wrong picking parameter suppied, allowed are 'atom/get/reset/status'");

    }

    return(TCL_ERROR);
}

/*********************************************************************/
int PipeCommand(ClientData clientdata , Tcl_Interp *interp,
                int argc , const char *argv[])
/*********************************************************************/
{
    static char  Text[BUFF_LEN];
    static char  Text1[BUFF_LEN];
    static char  Text2[BUFF_LEN];

    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

    if(Text[0] == (char)NULL) {
        gomp_PrintERROR("action parameter missing");
        return(TCL_ERROR);
    }

    gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

    if(gomp_StringMatch(Text , "open")) {
        if(Text1[0] == (char)NULL) {
            gomp_PrintERROR("program parameter missing");
            return(TCL_ERROR);
        }
        gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
        if(Text2[0] == (char)NULL) {
            gomp_PrintERROR("type parameter missing");
            return(TCL_ERROR);
        }

        if(gomp_OpenPipe2Program(Text1 , Text2)) {
            return(TCL_ERROR);
        } else {
            return(TCL_OK);
        }
    }
    else if(gomp_StringMatch(Text , "close")) {
        if(gomp_ClosePipe2Program()) {
            return(TCL_ERROR);
        } else {
            return(TCL_OK);
        }
    }
    else  {
        if(gomp_SendPipe2Program(Text)) 
            return(TCL_ERROR);
        else 
            return(TCL_OK);
    } 

}
#if defined(JUNK)
void gomp_Junk()
{
    printf("Hello Leif\n");
}
#endif
/*********************************************************************/
int TestCommand(ClientData clientdata , Tcl_Interp *interp,
                int argc , const char *argv[])
/*********************************************************************/
{
    static char  Text[BUFF_LEN];
/*    unsigned const int *address;*/

    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

    printf("%s\n",Text);
/*    address = (unsigned const int *)atoi(Text);*/
/*    *address = &gomp_Junk;*/

    return(TCL_OK);
}

/*********************************************************************/
const char *gomp_GetNextFromParserStack(int argc, const char **argv)
/*********************************************************************/
{
    static const char **point;
    static int          index;
    static int          total;
    static char         Empty = (char)NULL;

    if(  argv    != NULL) {
        index    = 0;
        total    = argc;
        point    = argv;
        if(argc  < 1) 
            return(&Empty);
        index  = 1;
    }
    else {
        if(index == total)
            return(&Empty);
        index++;
    }
    return(point[index-1]);     
}     
/*********************************************************************/
int   gomp_SendTclObjReturn(Tcl_Obj *obj)
/*********************************************************************/
{
    if(Interp == (Tcl_Interp *)NULL) {
        /* Free obj if it isn't shared. */
        Tcl_IncrRefCount(obj);
        Tcl_DecrRefCount(obj);
        return(1);
    }

    Tcl_SetObjResult(Interp, obj);
    return(0);
}
/*********************************************************************/
int   gomp_SendTclReturn(const char *Text)
/*********************************************************************/
{
    if(Interp == (Tcl_Interp *)NULL)
        return(1);

    Tcl_SetResult(Interp, CONST_CAST(char *, Text), TCL_VOLATILE);
    return(0);
}

/*********************************************************************/
int gomp_SendFile2TclParser(const char *File)
/*********************************************************************/
{
    return(Tcl_EvalFile(Interp , File));
}

/*********************************************************************/
Tcl_Interp  *gomp_GetTclInterp()
/*********************************************************************/
{
    return(Interp);
}

/*********************************************************************/
int          gomp_TclRunScript()
/*********************************************************************/
{
    static const char *value;
    static int         code;
    static Tcl_Interp *interp;

    interp = gomp_GetTclInterp();
    value  = Tcl_GetVar(interp , 
                                      TCL_RUN_DISPLAY_FILE_VAR , 
                                      TCL_GLOBAL_ONLY);

    if(value != (const char *)NULL) {
        code = Tcl_EvalFile(interp , value);
        if(code != TCL_OK) {
            if(Tcl_GetStringResult(interp) != (char *)NULL) {
                gomp_PrintERROR(Tcl_GetStringResult(interp));
                return(1);
            }
        }
    }
    return(0);
}

/*********************************************************************/
const char *gomp_GetOS()
/*********************************************************************/
{
    static const char *value;
    
    value  = Tcl_GetVar2(gomp_GetTclInterp() , "tcl_platform","os", TCL_GLOBAL_ONLY);
    
    if(!value) {
        gomp_PrintERROR("can't retrieve variable 'tcl_platform(os)'");
        return("");
    }

    return(value);
}
/*********************************************************************/
int WebURLCommand(ClientData clientdata , Tcl_Interp *interp,
                  int argc , const char *argv[])
/*********************************************************************/
{
    static char  Text[BUFF_LEN];

#if defined(WIN32)
    HWND     hwnd;
#endif

    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

    if(Text[0] == (char)NULL) {
        gomp_PrintERROR("Web URL is missing");
        return(1);
    }

#if defined(WIN32)
#ifdef ENABLE_GRAPHICS
#if defined(GLUT)
    hwnd   = GetActiveWindow();
#else
    hwnd   = auxGetHWND();
#endif
#else
    hwnd   = NULL;
#endif /* ENABLE_GRAPHICS */
    ShellExecute(hwnd, NULL, Text, NULL, NULL, SW_SHOW);
#endif

    return(0);
}

/*********************************************************************/
int PrintInterface(ClientData clientdata , Tcl_Interp *interp,
                   int argc , const char *argv[])
/*********************************************************************/
{
    static char  Text[BUFF_LEN];

    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

    if(Text[0] == (char)NULL) {
        return(0);
    }

    (void)gomp_PrintMessage(Text);

    return(0);
}

/*********************************************************************/
int ErrorInterface(ClientData clientdata , Tcl_Interp *interp,
                      int argc , const char *argv[])
/*********************************************************************/
{
    static char  Text[BUFF_LEN];

    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

    if(Text[0] == (char)NULL) {
        return(0);
    }

    (void)gomp_PrintERROR(Text);

    return(0);
}
#if 0
/*********************************************************************/
int          SetInputTclScriptFileName(const char *Name)
/*********************************************************************/
{
    memset(gomInputScriptName, 0, BUFF_LEN);

    if( Name!=NULL && strlen(Name) ) {
        strncpy(gomInputScriptName , Name , BUFF_LEN-1);
    }

    return(0);
}
#endif
/*********************************************************************/
const char *GetInputTclScriptFileName()
/*********************************************************************/
{
    return(gomInputScriptName);
}
/***************************************************************************/
int gomp_RungOpenMolResetScript()
/***************************************************************************/
{
    static int   ITemp;

    ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), "lulUtility::gOpenMolReset");
    if(ITemp != TCL_OK) {
        gomp_PrintERROR("can't execute 'lulUtility::gOpenMolReset' script");
        return(1);
    }

    return(0);
}

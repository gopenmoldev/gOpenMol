/*

Copyright (c) 1994 - 2005 by:
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

#ifndef WIN32
#include <sys/types.h>
/* Will choke on AIX if included by <Python.h>. */
#include <unistd.h>
#endif

#include "gomstring.h"
#include "gomtcl.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

#ifdef HAVE_LIBPYTHON
/* Python tries to link pythonxx_d.lib in debugging mode in Windows. */
#undef _DEBUG
#undef _POSIX_C_SOURCE
#undef _XOPEN_SOURCE
#undef vsnprintf
#undef snprintf
#include <Python.h>
#endif

static int SetPythonState(int);
static int GetPythonState(void);
static int PythonState = 0;

static void InitPythonExtensions(void);

/*********************************************************************/
int gomp_PythonCommand(ClientData clientdata, Tcl_Interp *interp,
                     int argc, const char **argv)
/*********************************************************************/
{
#ifdef HAVE_LIBPYTHON
    static char  Text[BUFF_LEN];
    static const char *Value;

/* mopen */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "pyth$on")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

        if(strlen(Text) == 0) {
            gomp_PrintERROR("python state value missing (on/off)");
            return(TCL_ERROR);
        }

        if(gomp_StringMatch(Text , "on")) {
            if(!GetPythonState()) {
                Py_Initialize();
                InitPythonExtensions();
            }
            Value = Tcl_SetVar(gomp_GetTclInterp(),"gomParser","python",TCL_GLOBAL_ONLY);
            if(!Value) {
                gomp_PrintERROR("can't set tcl variable 'gomParser'");
                return(TCL_ERROR);
            }
            SetPythonState(ON);
        } else if(gomp_StringMatch(Text , "off")) {
            if(GetPythonState())
                Py_Finalize();
            Value = Tcl_SetVar(gomp_GetTclInterp(),"gomParser","tcl",TCL_GLOBAL_ONLY);
            if(!Value) {
                gomp_PrintERROR("can't set tcl variable 'gomParser'");
                return(TCL_ERROR);
            }
            SetPythonState(OFF);
        } else {
            if(!GetPythonState()) {
                Py_Initialize();
                InitPythonExtensions();
            }
/*          printf("'%s'\n",Py_GetPath()); */
            if(!PyRun_SimpleString(CONST_CAST(char *, argv[1])))
                return(TCL_OK);
            else 
                return(TCL_ERROR);
        }

        return(TCL_OK);
    }
/*  E R R O R command not recognized         */
    gomp_PrintERROR("'python' command not recognized");
    return(TCL_ERROR);
#else /* HAVE_LIBPYTHON */
    gomp_PrintERROR("Python support is not compiled into gOpenMol");
    return(1);
#endif
}
/*********************************************************************/
int SetPythonState(int Value)
/*********************************************************************/
{
    PythonState = Value;

    return(0);
}
/*********************************************************************/
int GetPythonState()
/*********************************************************************/
{
    return(PythonState);
}

#ifdef HAVE_LIBPYTHON
/*********************************************************************/
/* module functions */
/*********************************************************************/
static PyObject *                                 /* returns object */
tclparser_get(PyObject *self, PyObject *args)           /* self unused in modules */
/*********************************************************************/
{                                                 /* args from python call */
    const char *fromPython;
    const char *Value;

    /* convert Python -> C */
    if (! PyArg_Parse(args, CONST_CAST(char *, "(s)"), &fromPython)) 
        return NULL;                              /* null=raise exception */
    else {
        Value = Tcl_GetVar(gomp_GetTclInterp(),fromPython,TCL_GLOBAL_ONLY);
        if(!Value) {
            return NULL;      
        }
        return Py_BuildValue(CONST_CAST(char *, "s"), Value); /* convert C -> Python */
    }
}

/*********************************************************************/
static PyObject *                                 /* returns object */
tclparser_set(PyObject *self, PyObject *args)           /* self unused in modules */
/*********************************************************************/
{                                                 /* args from python call */
    const char *fromPython1,*fromPython2;
    const char *Value;

    if (! PyArg_Parse(args, CONST_CAST(char *, "(ss)"),
                      &fromPython1, &fromPython2))  
        /* convert Python -> C */
        return NULL;                              /* null=raise exception */
    else {
        Value = Tcl_SetVar(gomp_GetTclInterp(),fromPython1,fromPython2,TCL_GLOBAL_ONLY);
        if(!Value) {
            gomp_PrintERROR("can't set the tcl variable");
            return(NULL);
        }
        /* convert C -> Python */
        return Py_BuildValue(CONST_CAST(char *, "s"), fromPython2);
    }
}

/* registration table  */
static char get[] = "get";
static char set[] = "set";
static struct PyMethodDef tclparser_methods[] = {
    {get, tclparser_get, 1}, /* method name, C func ptr, always-tuple */
    {set, tclparser_set, 1}, /* method name, C func ptr, always-tuple */
    {NULL, NULL}                     /* end of table marker */
};



/* module initializer */
/* called on first import */
void InitPythonExtensions()
{
    /* name matters if loaded dynamically */
    static char name[] = "tclparser";
    (void)Py_InitModule(name, tclparser_methods);
    /* mod name, table ptr */
}

#endif /* HAVE_LIBPYTHON */

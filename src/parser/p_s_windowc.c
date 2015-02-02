/*

Copyright (c) 2003 - 2004 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#include "gomproc.h"
#include "gomwindow.h"
#include "parser.h"
#include "projview.h"
#include "text_stack.h"

#include "stdafx.h"

#ifdef ENABLE_GRAPHICS

/* show window defined */
#define SHOW_WINDOW_DEFINED_ENTRY \
GOM_PARSER_FINAL_INT_CMD("defined",&parseShowWindowDefined)
const static IntFunc parseShowWindowDefined = { gomp_GetNumDefinedWindows };

/* show window parameters */
#define  SHOW_WINDOW_PARAMETERS_ENTRY \
GOM_PARSER_FINAL_CMD("parameters",\
                     GOM_PARSER_NEED_ARGS_RANGE("?WinID?",0,1),\
                     ParseShowWindowParameters,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowWindowParameters(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    int WinID;
    int winX, winY, winW, winH;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(WinID)
        GOM_PARSER_END_ARGS;
    GOM_PARSER_VERIFY_GOM_PARSER(
        "windows id",
        GOM_PARSER_DEFAULT(WinID,1) ||
        GOM_PARSER_RETRIEVE_INT(WinID) );

    if ( gomp_SetWindow(WinID) ) {
        gomp_PrintERROR("window index out of allowed range");
        GOM_PARSER_FAILED;
    }

    gomp_GetWindowParameters(&winX, &winY, &winW, &winH);
    GOM_PARSER_RETURN_LIST(("%d %d %d %d",winX,winY,winW,winH));
    GOM_PARSER_SUCCEEDED;
}

/* show window redraw */
#define SHOW_WINDOW_REDRAW_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD("redraw",&parseShowWindowRedraw)
/* show window update */
#define SHOW_WINDOW_UPDATE_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD("update",&parseShowWindowRedraw)
static const BooleanEnumFunc parseShowWindowRedraw = {
    gomp_GetUpdateDisplayMode, { "manual", "automatic" } };

/* show window style */
#define SHOW_WINDOW_STYLE_ENTRY \
GOM_PARSER_FINAL_CMD("style",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowWindowStyle,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowWindowStyle(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    switch ( gomp_GetWindowingStyle() ) {
    case SINGLE_WINDOWING:
        GOM_PARSER_RETURN_STRING("single");
        break;
    case MULTI_WINDOWING:
        GOM_PARSER_RETURN_STRING("multi");
        break;
    default:
        gomp_PrintERROR("wrong window mode (single/multi)");
        GOM_PARSER_FAILED;
    }
    GOM_PARSER_SUCCEEDED;
}

const gom_ParserArgumentList gomp_ShowWindowCommand[] = {
    SHOW_WINDOW_DEFINED_ENTRY,
    SHOW_WINDOW_PARAMETERS_ENTRY,
    SHOW_WINDOW_REDRAW_ENTRY,
    SHOW_WINDOW_STYLE_ENTRY,
    SHOW_WINDOW_UPDATE_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};
#else /* ENABLE_GRAPHICS */
extern int i;
#endif

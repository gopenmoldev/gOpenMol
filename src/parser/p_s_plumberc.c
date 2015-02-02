/*

Copyright (c) 2003 - 2004 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#include <stdlib.h>

#include "listutils.h"
#include "parser.h"
#include "plumber.h"

#include "stdafx.h"

#define ALIAS GOM_PARSER_CMD_ALIAS

static int ParseGetPlumberIndex(int *pIplumber,GOM_PARSER_ARGLIST)
{
    int Iplumber;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Iplumber)
    GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "plumber index",
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Iplumber, 1, gomp_GetPlumberSets() ) );

    *pIplumber = Iplumber;

    GOM_PARSER_SUCCEEDED;
}

#define GET_PLUMBER_INDEX \
    int Iplumber; \
    GOM_PARSER_VERIFY_CHILD(\
        ParseGetPlumberIndex(&Iplumber,GOM_PARSER_PASS_ARGS))

/* show plumber count */
#define SHOW_PLUMBER_COUNT_ENTRY \
GOM_PARSER_FINAL_INT_CMD("count",&parseShowPlumberCount)
static const IntFunc parseShowPlumberCount = { gomp_GetPlumberSets };

/* show plumber displaystate */
#define SHOW_PLUMBER_DISPLAY_STATE_ENTRY \
GOM_PARSER_FINAL_INT_CMD("displaystate",&parseShowPlumberDisplayState)
static const IntFunc parseShowPlumberDisplayState = { gomp_GetPlumberDisplay };

/* show plumber colour Iplumber */
/* show plumber color  Iplumber */
#define SHOW_PLUMBER_COLOUR_ENTRY \
GOM_PARSER_FINAL_CMD("colour" ALIAS "color",\
                     GOM_PARSER_NEED_ARGS("Iplumber",1),\
                     ParseShowPlumberColour,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowPlumberColour(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_PLUMBER_INDEX;
    GOM_PARSER_RETURN_LIST(("%f %f %f",
                            gomp_GetPlumberRed  (Iplumber - 1),
                            gomp_GetPlumberGreen(Iplumber - 1),
                            gomp_GetPlumberBlue (Iplumber - 1)));
    GOM_PARSER_SUCCEEDED;
}

/* show plumber structure Iplumber */
#define SHOW_PLUMBER_STRUCTURE_ENTRY \
GOM_PARSER_FINAL_CMD("structure",\
                     GOM_PARSER_NEED_ARGS("Iplumber",1),\
                     ParseShowPlumberStructure,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowPlumberStructure(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_PLUMBER_INDEX;
    GOM_PARSER_RETURN_INT(gomp_GetPlumberStructure(Iplumber - 1) + 1);
    GOM_PARSER_SUCCEEDED;
}

/* show plumber atoms Iplumber */
#define SHOW_PLUMBER_ATOMS_ENTRY \
GOM_PARSER_FINAL_CMD("atoms",\
                     GOM_PARSER_NEED_ARGS("Iplumber",1),\
                     ParseShowPlumberAtoms,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowPlumberAtoms(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    char *atom_list;
    GET_PLUMBER_INDEX;
    atom_list = gomp_MakeIndexList(
        gomp_GetPlumberAtoms(Iplumber - 1),
        gomp_GetPlumberAtomList(Iplumber - 1),1,1,',');
    if ( atom_list ) {
        GOM_PARSER_RETURN_STRING(atom_list);
    }
    free(atom_list);
    GOM_PARSER_SUCCEEDED;
}

/* show plumber type Iplumber */
#define SHOW_PLUMBER_TYPE_ENTRY \
GOM_PARSER_FINAL_CMD("type",\
                     GOM_PARSER_NEED_ARGS("Iplumber",1),\
                     ParseShowPlumberType,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowPlumberType(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_PLUMBER_INDEX;
    switch ( gomp_GetPlumberDisplayType(Iplumber - 1) ) {
    case RIBBON_TYPE:      GOM_PARSER_RETURN_STRING("ribbon"); break;
    case CYLINDER_TYPE:    GOM_PARSER_RETURN_STRING("tube");   break;
    case FLAT_HELIX_TYPE:  GOM_PARSER_RETURN_STRING("fhelix"); break;
    case SOLID_HELIX_TYPE: GOM_PARSER_RETURN_STRING("shelix"); break;
    case ARROW_TYPE:       GOM_PARSER_RETURN_STRING("arrow");  break;
    case STRAND_TYPE:      GOM_PARSER_RETURN_STRING("strand"); break;
    case TRACE_TYPE:       GOM_PARSER_RETURN_STRING("trace");  break;
    }
    GOM_PARSER_SUCCEEDED;
}

const gom_ParserArgumentList gomp_ShowPlumberCommand[] = {
    SHOW_PLUMBER_COUNT_ENTRY,
    SHOW_PLUMBER_DISPLAY_STATE_ENTRY,
    SHOW_PLUMBER_COLOUR_ENTRY,
    SHOW_PLUMBER_STRUCTURE_ENTRY,
    SHOW_PLUMBER_ATOMS_ENTRY,
    SHOW_PLUMBER_TYPE_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

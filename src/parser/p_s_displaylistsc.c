/*

Copyright (c) 2003 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#include "parser.h"
#include "plot.h"

#include "stdafx.h"

/* show displaylists default */
#define SHOW_DISPLAYLISTS_DEFAULT_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_CMD("default",&parseShowDisplayListsDefault)
const static IntFunc parseShowDisplayListsDefault = {
    gomp_GetDefaultDisplayListState };

/* show displaylists types */
#define SHOW_DISPLAYLISTS_TYPES_ENTRY \
GOM_PARSER_FINAL_VERIFY_FUNC_CMD("types",&parseShowDisplayListsTypes)
const static IntFunc parseShowDisplayListsTypes = {
    gomp_ParseGetObjectDisplayListTypes };

/* show displaylists Type */
/* show displaylists      */
#define SHOW_DISPLAYLISTS_STATUS_ENTRY \
GOM_PARSER_FINAL_ARGS(GOM_PARSER_NEED_ARGS("Type",1),\
                      ParseDisplayListsStatus,\
                      GOM_PARSER_UNUSED_VALUE),\
GOM_PARSER_FINAL_ARGS(GOM_PARSER_NO_MORE_ARGS,\
                      gomp_ParserReturnBooleanValueFromFunc,\
                      GOM_PARSER_SET_POINTER_VALUE(&parseShowDisplayStatus))
static int ParseDisplayListsStatus(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    const char *Type;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Type)
    GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "object type",
        GOM_PARSER_RETRIEVE_STRING(Type));
    
    GOM_PARSER_VERIFY(gomp_ParseGetObjectDisplayListState(Type) == 0);

    GOM_PARSER_SUCCEEDED;
}
const static IntFunc parseShowDisplayStatus = {
    gomp_GetDisplayListState };

const gom_ParserArgumentList gomp_ShowDisplaylistsCommand[] = {
    SHOW_DISPLAYLISTS_DEFAULT_ENTRY,
    SHOW_DISPLAYLISTS_TYPES_ENTRY,
    SHOW_DISPLAYLISTS_STATUS_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

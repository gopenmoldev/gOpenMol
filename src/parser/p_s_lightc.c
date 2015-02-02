/*

Copyright (c) 2003 - 2004 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#include "light_model.h"
#include "parser.h"

#include "stdafx.h"

/* show light diffuse {red|green|blue} */
#define SHOW_LIGHT_DIFFUSE_ENTRY \
GOM_PARSER_CMD_PART("diffuse",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowLightDiffuse,\
                    GOM_PARSER_UNUSED_VALUE)
const static FloatFunc parseShowLightDiffuseRed   = {
    gomp_GetLightDiffuseRed };
const static FloatFunc parseShowLightDiffuseGreen = {
    gomp_GetLightDiffuseGreen };
const static FloatFunc parseShowLightDiffuseBlue  = {
    gomp_GetLightDiffuseBlue };
const static gom_ParserArgumentList parseShowLightDiffuse[] = {
    GOM_PARSER_FINAL_FLOAT_CMD("red",  &parseShowLightDiffuseRed),
    GOM_PARSER_FINAL_FLOAT_CMD("green",&parseShowLightDiffuseGreen),
    GOM_PARSER_FINAL_FLOAT_CMD("blue", &parseShowLightDiffuseBlue),
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show light position */
#define SHOW_LIGHT_POSITION_ENTRY \
GOM_PARSER_FINAL_CMD("position",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowLightPosition,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowLightPosition(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    const static char *pos[9] = {
        "c",
        "n","ne","e","se","s","sw","w","nw"
    };
    GOM_PARSER_RETURN_STRING(pos[gomp_GetLightPosition()]);
    GOM_PARSER_SUCCEEDED;
}

const gom_ParserArgumentList gomp_ShowLightCommand[] = {
    SHOW_LIGHT_DIFFUSE_ENTRY,
    SHOW_LIGHT_POSITION_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

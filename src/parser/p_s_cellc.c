/*

Copyright (c) 2003 - 2004 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#include "cell.h"
#include "parser.h"

#include "stdafx.h"

#define ALIAS GOM_PARSER_CMD_ALIAS

/* show cell colour */
/* show cell color  */
#define SHOW_CELL_COLOUR_ENTRY \
GOM_PARSER_FINAL_CMD("colour" ALIAS "color",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowCellColour,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowCellColour(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    float red, green, blue;
    (void)gomp_GetCellColour(&red, &green, &blue);
    GOM_PARSER_RETURN_LIST(("%f %f %f", red,  green, blue));
    GOM_PARSER_SUCCEEDED;
}

/* show cell linewidth */
#define SHOW_CELL_LINE_WIDTH_ENTRY \
GOM_PARSER_FINAL_INT_CMD("linewidth",&parseShowCellLineWidth)
static const IntFunc parseShowCellLineWidth = { gomp_GetCellLinewidth };

/* show cell dimensions */
#define SHOW_CELL_DIMENSION_ENTRY \
GOM_PARSER_FINAL_CMD("dimensions",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowCellDimensions,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowCellDimensions(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GOM_PARSER_RETURN_LIST(("%f %f %f",
                            gomp_GetCellA(),
                            gomp_GetCellB(),
                            gomp_GetCellC()));
    GOM_PARSER_SUCCEEDED;
}

/* show cell position    */
/* show cell translation */
#define SHOW_CELL_POSITION_ENTRY \
GOM_PARSER_FINAL_CMD("position",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowCellPosition,\
                     GOM_PARSER_UNUSED_VALUE),\
GOM_PARSER_FINAL_CMD("translation",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowCellPosition,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowCellPosition(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GOM_PARSER_RETURN_LIST(("%f %f %f",
                            gomp_GetCellXtrans(),
                            gomp_GetCellYtrans(),
                            gomp_GetCellZtrans()));
    GOM_PARSER_SUCCEEDED;
}

const gom_ParserArgumentList gomp_ShowCellCommand[] = {
    SHOW_CELL_COLOUR_ENTRY,
    SHOW_CELL_LINE_WIDTH_ENTRY,
    SHOW_CELL_DIMENSION_ENTRY,
    SHOW_CELL_POSITION_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

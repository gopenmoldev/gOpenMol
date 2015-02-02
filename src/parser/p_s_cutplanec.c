/*

Copyright (c) 2003 - 2004 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#include "contour.h"
#include "parser.h"

#include "stdafx.h"

#define ALIAS GOM_PARSER_CMD_ALIAS

/* show cutplane smooth x  */
/* show cutplane smooth yz */
#define SHOW_CUT_PLANE_SMOOTH_X_ENTRY \
GOM_PARSER_FINAL_INT_CMD("x" ALIAS "yz",&parseShowCutPlaneSmoothX)
static const IntFunc parseShowCutPlaneSmoothX = { gomp_GetCutplaneSmoothTypeX};

/* show cutplane smooth y  */
/* show cutplane smooth xz */
#define SHOW_CUT_PLANE_SMOOTH_Y_ENTRY \
GOM_PARSER_FINAL_INT_CMD("y" ALIAS "xz",&parseShowCutPlaneSmoothY)
static const IntFunc parseShowCutPlaneSmoothY = { gomp_GetCutplaneSmoothTypeY};

/* show cutplane smooth z  */
/* show cutplane smooth xy */
#define SHOW_CUT_PLANE_SMOOTH_Z_ENTRY \
GOM_PARSER_FINAL_INT_CMD("z" ALIAS "xy",&parseShowCutPlaneSmoothZ)
static const IntFunc parseShowCutPlaneSmoothZ = { gomp_GetCutplaneSmoothTypeZ};

/* show cutplane smooth ... */
#define SHOW_CUT_PLANE_SMOOTH_ENTRY \
GOM_PARSER_CMD_PART("smooth",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowCutPlaneSmooth,\
                    GOM_PARSER_UNUSED_VALUE)
static const gom_ParserArgumentList parseShowCutPlaneSmooth[] = {
    SHOW_CUT_PLANE_SMOOTH_X_ENTRY,
    SHOW_CUT_PLANE_SMOOTH_Y_ENTRY,
    SHOW_CUT_PLANE_SMOOTH_Z_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show cutplane type */
#define SHOW_CUT_PLANE_TYPE_ENTRY \
GOM_PARSER_FINAL_CMD("type",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowCutPlaneType,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowCutPlaneType(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GOM_PARSER_RETURN_LIST(("%s %s %s",
                            gomp_GetCutPlaneType_3D_X() ? "3dx" : "2dx",
                            gomp_GetCutPlaneType_3D_Y() ? "3dy" : "2dy",
                            gomp_GetCutPlaneType_3D_Z() ? "3dz" : "2dz"));
    GOM_PARSER_SUCCEEDED;
}

/* show cutplane damping */
#define SHOW_CUT_PLANE_DAMPING_ENTRY \
GOM_PARSER_FINAL_FLOAT_CMD("damping",&parseShowCutPlaneDamping)
static const FloatFunc parseShowCutPlaneDamping = { gomp_GetCutPlaneDamping };

/* show cutplane xyz1 */
/* show cutplane xyz2 */
/* show cutplane xyz3 */
#define SHOW_CUT_PLANE_XYZ_ENTRY \
GOM_PARSER_FINAL_CMD("xyz1",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowCutPlaneXYZ,\
                     GOM_PARSER_SET_VALUE(1)),\
GOM_PARSER_FINAL_CMD("xyz2",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowCutPlaneXYZ,\
                     GOM_PARSER_SET_VALUE(2)),\
GOM_PARSER_FINAL_CMD("xyz3",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowCutPlaneXYZ,\
                     GOM_PARSER_SET_VALUE(3))
static int ParseShowCutPlaneXYZ(GOM_PARSER_ARGLIST,intptr_t Which)
{
    GOM_PARSER_RETURN_INT(gomp_GetContourXYZClippingPlaneState(Which));
    GOM_PARSER_SUCCEEDED;
}

const gom_ParserArgumentList gomp_ShowCutPlaneCommand[] = {
    SHOW_CUT_PLANE_SMOOTH_ENTRY,
    SHOW_CUT_PLANE_TYPE_ENTRY,
    SHOW_CUT_PLANE_DAMPING_ENTRY,
    SHOW_CUT_PLANE_XYZ_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

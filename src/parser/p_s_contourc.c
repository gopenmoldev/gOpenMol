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

#define GET_CONTOUR_INDEX_BASE \
    GOM_PARSER_VERIFY_GOM_PARSER( \
        "contour index", \
        GOM_PARSER_DEFAULT(Scontour,"1") || \
        GOM_PARSER_RETRIEVE_STRING(Scontour) ); \
    Icontour = gomp_CheckContourName(Scontour); \
    GOM_PARSER_VERIFY( Icontour > 0 )

/* ... Icontour */
static int ParseGetContourIndex(int *pIcontour,GOM_PARSER_ARGLIST)
{
    const char *Scontour;
    int         Icontour;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Scontour)
    GOM_PARSER_END_ARGS;

    GET_CONTOUR_INDEX_BASE;

    *pIcontour = Icontour;
    
    GOM_PARSER_SUCCEEDED;
}

#define GET_CONTOUR_INDEX \
    int Icontour; \
    GOM_PARSER_VERIFY_CHILD( \
        ParseGetContourIndex(&Icontour,GOM_PARSER_PASS_ARGS))

#define GET_CONTOUR_AND_LEVEL_INDECES_BASE \
    GET_CONTOUR_INDEX_BASE; \
    GOM_PARSER_VERIFY_GOM_PARSER( \
        "contour level index", \
        ( GOM_PARSER_DEFAULT(Ilevel,1) || \
          GOM_PARSER_RETRIEVE_INT(Ilevel) ) && \
        GOM_PARSER_CHECK_INT_RANGE( \
            Ilevel, 1, gomp_GetContourLevels(Icontour - 1) ) )

/* ... Icontour Ilevel */   
static int ParseGetContourAndLevelIndeces(int *pIcontour, int *pIlevel,
                                          GOM_PARSER_ARGLIST)
{
    const char *Scontour;
    int         Icontour, Ilevel;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Scontour)
        GOM_PARSER_ARG(Ilevel)
    GOM_PARSER_END_ARGS;

    GET_CONTOUR_AND_LEVEL_INDECES_BASE;

    *pIcontour = Icontour;
    *pIlevel   = Ilevel;

    GOM_PARSER_SUCCEEDED;
}

#define GET_CONTOUR_AND_LEVEL_INDECES \
    int Icontour, Ilevel; \
    GOM_PARSER_VERIFY_CHILD( \
        ParseGetContourAndLevelIndeces(&Icontour,&Ilevel,GOM_PARSER_PASS_ARGS))

/* show contour defined */
#define SHOW_CONTOUR_DEFINED_ENTRY \
GOM_PARSER_FINAL_INT_CMD("defined",&parseShowContourDefined)
static const IntFunc parseShowContourDefined = { gomp_GetContoursDefined };

/* show contour filename Icontour */
#define SHOW_CONTOUR_FILE_NAME_ENTRY \
GOM_PARSER_FINAL_CMD("filename",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Icontour?",0,1),\
                     ParseShowContourFileName,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourFileName(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_INDEX;
    GOM_PARSER_RETURN_STRING(gomp_GetContourFileName(Icontour - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show contour structure Icontour */
#define SHOW_CONTOUR_STRUCTURE_ENTRY \
GOM_PARSER_FINAL_CMD("structure",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Icontour?",0,1),\
                     ParseShowContourStructure,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourStructure(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_INDEX;
    GOM_PARSER_RETURN_INT(gomp_GetContour2StructureMapping(Icontour - 1) + 1);
    GOM_PARSER_SUCCEEDED;
}

/* show contour minimum Icontour */
#define SHOW_CONTOUR_MINIMUM_ENTRY \
GOM_PARSER_FINAL_CMD("minimum",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Icontour?",0,1),\
                     ParseShowContourMinimum,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourMinimum(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_INDEX;
    GOM_PARSER_RETURN_DOUBLE(gomp_GetContourMin(Icontour - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show contour maximum Icontour */
#define SHOW_CONTOUR_MAXIMUM_ENTRY \
GOM_PARSER_FINAL_CMD("maximum",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Icontour?",0,1),\
                     ParseShowContourMaximum,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourMaximum(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_INDEX;
    GOM_PARSER_RETURN_DOUBLE(gomp_GetContourMax(Icontour - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show contour header Icontour */
#define SHOW_CONTOUR_HEADER_ENTRY \
GOM_PARSER_FINAL_CMD("header",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Icontour?",0,1),\
                     ParseShowContourHeader,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourHeader(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_INDEX;
    GOM_PARSER_RETURN_LIST(("%d %d %d   %f %f   %f %f   %f %f",
                            gomp_GetContourPointsX(Icontour - 1),
                            gomp_GetContourPointsY(Icontour - 1),
                            gomp_GetContourPointsZ(Icontour - 1),

                            gomp_GetContourMinX(Icontour - 1),
                            gomp_GetContourMaxX(Icontour - 1),

                            gomp_GetContourMinY(Icontour - 1),
                            gomp_GetContourMaxY(Icontour - 1),

                            gomp_GetContourMinZ(Icontour - 1),
                            gomp_GetContourMaxZ(Icontour - 1)));
    GOM_PARSER_SUCCEEDED;
}

/* show contour levels Icontour */
#define SHOW_CONTOUR_LEVELS_ENTRY \
GOM_PARSER_FINAL_CMD("levels",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Icontour?",0,1),\
                     ParseShowContourLevels,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourLevels(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_INDEX;
    GOM_PARSER_RETURN_INT(gomp_GetContourLevels(Icontour - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show contour values Icontour Ilevel */
#define SHOW_CONTOUR_VALUES_ENTRY \
GOM_PARSER_FINAL_CMD("values",\
                     GOM_PARSER_NEED_ARGS("Icontour Ilevel",2),\
                     ParseShowContourValues,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourValues(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_AND_LEVEL_INDECES;
    GOM_PARSER_RETURN_LIST(("%f %f %f %f",
                            gomp_GetContourValue(Icontour - 1, Ilevel - 1),
                            gomp_GetContourColourRed  (Icontour-1,Ilevel-1),
                            gomp_GetContourColourGreen(Icontour-1,Ilevel-1),
                            gomp_GetContourColourBlue (Icontour-1,Ilevel-1)));
    GOM_PARSER_SUCCEEDED;
}

/* show contour smooth Icontour Ilevel */
#define SHOW_CONTOUR_SMOOTH_ENTRY \
GOM_PARSER_FINAL_CMD("smooth",\
                     GOM_PARSER_NEED_ARGS("Icontour Ilevel",2),\
                     ParseShowContourSmooth,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourSmooth(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_AND_LEVEL_INDECES;
    GOM_PARSER_RETURN_BOOLEAN(gomp_GetContourSmooth(Icontour - 1, Ilevel - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show contour cullface Icontour Ilevel */
#define SHOW_CONTOUR_CULLFACE_ENTRY \
GOM_PARSER_FINAL_CMD("cullface",\
                     GOM_PARSER_NEED_ARGS("Icontour Ilevel",2),\
                     ParseShowContourCullface,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourCullface(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_AND_LEVEL_INDECES;
    GOM_PARSER_RETURN_STRING(
        gomp_GetContourCullFace(Icontour - 1, Ilevel - 1) ? "on" : "off");
    GOM_PARSER_SUCCEEDED;
}

/* show contour clipplane direction Icontour Ilevel */
#define SHOW_CONTOUR_CLIPPLANE_DIRECTION_ENTRY \
GOM_PARSER_FINAL_CMD("direction",\
                     GOM_PARSER_NEED_ARGS("Icontour Ilevel",2),\
                     ParseShowContourClipDirection,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourClipDirection(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    char dir[2] = "";
    GET_CONTOUR_AND_LEVEL_INDECES;
    dir[0] = gomp_GetContourLevelClippingPlaneAxis(Icontour - 1, Ilevel - 1);
    GOM_PARSER_RETURN_STRING(dir);
    GOM_PARSER_SUCCEEDED;
}

/* show contour clipplane position Icontour Ilevel */
#define SHOW_CONTOUR_CLIPPLANE_POSITION_ENTRY \
GOM_PARSER_FINAL_CMD("position",\
                     GOM_PARSER_NEED_ARGS("Icontour Ilevel",2),\
                     ParseShowContourClipPosition,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourClipPosition(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_AND_LEVEL_INDECES;
    GOM_PARSER_RETURN_DOUBLE(
        gomp_GetContourLevelClippingPlanePosition(Icontour - 1, Ilevel - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show contour clipplane ... */
#define SHOW_CONTOUR_CLIPPLANE_ENTRY \
GOM_PARSER_CMD_PART("clipplane",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowContourClip,\
                    GOM_PARSER_UNUSED_VALUE)
static const gom_ParserArgumentList parseShowContourClip[] = {
    SHOW_CONTOUR_CLIPPLANE_DIRECTION_ENTRY,
    SHOW_CONTOUR_CLIPPLANE_POSITION_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show contour display Icontour Ilevel */
#define SHOW_CONTOUR_DISPLAY_ENTRY \
GOM_PARSER_FINAL_CMD("display",\
                     GOM_PARSER_NEED_ARGS("Icontour Ilevel",2),\
                     ParseShowContourDisplay,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourDisplay(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_AND_LEVEL_INDECES;
    GOM_PARSER_RETURN_BOOLEAN(
        gomp_GetContourDisplayState(Icontour - 1, Ilevel - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show contour name Icontour */
#define SHOW_CONTOUR_NAME_ENTRY \
GOM_PARSER_FINAL_CMD("name",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Icontour?",0,1),\
                     ParseShowContourName,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourName(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_INDEX;
    GOM_PARSER_RETURN_STRING(gomp_GetContourName(Icontour - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show contour mapping    Icontour */
/* show contour projection Icontour */
#define SHOW_CONTOUR_MAPPING_ENTRY \
GOM_PARSER_FINAL_CMD("mapping",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Icontour?",0,1),\
                     ParseShowContourMapping,\
                     GOM_PARSER_UNUSED_VALUE),\
GOM_PARSER_FINAL_CMD("projection",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Icontour?",0,1),\
                     ParseShowContourMapping,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourMapping(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    int projIndex;
    GET_CONTOUR_INDEX;
    projIndex = gomp_GetContourProjection(Icontour - 1);
    if (  projIndex == Icontour - 1 ) {
        GOM_PARSER_RETURN_LIST(("%s %d",
                                gomp_GetContourName(projIndex),0));
        GOM_PARSER_SUCCEEDED;
    }
    else {
        gom_ParserList list;
        int i;
        if ( ! GOM_PARSER_LIST_INIT(list) )
            goto failed;
        
        /* name of a contour */
        if ( ! GOM_PARSER_LIST_APPEND_STRING(
                 list,gomp_GetContourName(projIndex)) )
            goto failed;
        /* number of levels */
        if ( ! GOM_PARSER_LIST_APPEND_INT(
                 list, gomp_GetContourLevels(Icontour - 1)) )
            goto failed;
            
        /* projection values for each level. */
        for ( i = 0 ; i < gomp_GetContourLevels(Icontour - 1) ; i++ ) {
            if ( ! GOM_PARSER_LIST_APPEND_DOUBLE(
                     list,gomp_GetContourProjectionMin(Icontour - 1, i)) ||
                 ! GOM_PARSER_LIST_APPEND_DOUBLE(
                     list,gomp_GetContourProjectionMax(Icontour - 1, i)) )
                goto failed;
        }
        GOM_PARSER_LIST_RETURN(list);
        GOM_PARSER_SUCCEEDED;
      failed:
        GOM_PARSER_LIST_FREE(list);
        GOM_PARSER_FAILED;
    }
}

/* show contour type Icontour Ilevel */
#define SHOW_CONTOUR_TYPE_ENTRY \
GOM_PARSER_FINAL_CMD("type",\
                     GOM_PARSER_NEED_ARGS("Icontour Ilevel",2),\
                     ParseShowContourType,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourType(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_AND_LEVEL_INDECES;
    GOM_PARSER_RETURN_INT(gomp_GetContourDisplayType(Icontour - 1,Ilevel - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show contour transparency Icontour Ilevel */
/* show contour alphablend   Icontour Ilevel */
#define SHOW_CONTOUR_TRANSPARENCY_ENTRY \
GOM_PARSER_FINAL_CMD("transparency",\
                     GOM_PARSER_NEED_ARGS("Icontour Ilevel",2),\
                     ParseShowContourTransparency,\
                     GOM_PARSER_UNUSED_VALUE),\
GOM_PARSER_FINAL_CMD("alphablend",\
                     GOM_PARSER_NEED_ARGS("Icontour Ilevel",2),\
                     ParseShowContourTransparency,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourTransparency(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_CONTOUR_AND_LEVEL_INDECES;
    GOM_PARSER_RETURN_DOUBLE(gomp_GetContourAlpha(Icontour - 1, Ilevel - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show contour cube Icontour */
/* show contour cell Icontour */
#define SHOW_CONTOUR_CUBE_ENTRY \
GOM_PARSER_FINAL_CMD("cube",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Icontour?",0,1),\
                     ParseShowContourCube,\
                     GOM_PARSER_UNUSED_VALUE),\
GOM_PARSER_FINAL_CMD("cell",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Icontour?",0,1),\
                     ParseShowContourCube,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourCube(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    struct {
        float min, max;
        int   dim;
    } x,y,z;
    GET_CONTOUR_INDEX;
    (void)gomp_GetContourCube(Icontour - 1,
                              &x.min, &x.max, &y.min, &y.max, &z.min, &z.max,
                              &x.dim,         &y.dim,         &z.dim);
    GOM_PARSER_RETURN_LIST(("%f %f   %f %f   %f %f   %d %d %d",
                            x.min, x.max,
                            y.min, y.max,
                            z.min, z.max,
                            x.dim, y.dim, z.dim));
    GOM_PARSER_SUCCEEDED;
}

/* show contour polygon method */
#define SHOW_CONTOUR_POLYGON_METHOD_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD(\
    "method",&parseShowContourPolygonMethod)
static const BooleanEnumFunc parseShowContourPolygonMethod = {
    gomp_GetSurfaceMethod, { "direct", "save" } };

/* show contour polygon entries Ilevel Icontour */
#define SHOW_CONTOUR_POLYGON_ENTRIES_ENTRY \
GOM_PARSER_FINAL_CMD("entries",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Ilevel Icontour?",0,2),\
                     ParseShowContourPolygonEntries,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourPolygonEntries(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    const char *Scontour;
    int         Icontour, Ilevel;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Ilevel)
        GOM_PARSER_ARG(Scontour)
    GOM_PARSER_END_ARGS;

    GET_CONTOUR_AND_LEVEL_INDECES_BASE;

    GOM_PARSER_RETURN_INT(gomp_ShowNumberOfPolygons(Icontour - 1, Ilevel - 1));

    GOM_PARSER_SUCCEEDED;
}

/* show contour polygon data Ientry Ilevel Icontour */
#define SHOW_CONTOUR_POLYGON_DATA_ENTRY \
GOM_PARSER_FINAL_CMD("data",\
                     GOM_PARSER_NEED_ARGS_RANGE(\
                         "?Ientry Ilevel Icontour?",0,3),\
                     ParseShowContourPolygonData,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowContourPolygonData(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    const char *Scontour;
    int         Icontour, Ilevel, Ientry;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Ientry)
        GOM_PARSER_ARG(Ilevel)
        GOM_PARSER_ARG(Scontour)
    GOM_PARSER_END_ARGS;

    GET_CONTOUR_AND_LEVEL_INDECES_BASE;
    GOM_PARSER_VERIFY_GOM_PARSER(
        "polygon entry index",
        ( GOM_PARSER_DEFAULT(Ientry,1) ||
          GOM_PARSER_RETRIEVE_INT(Ientry) ) &&
        GOM_PARSER_CHECK_INT_RANGE(
            Ientry, 1, gomp_ShowNumberOfPolygons(Icontour - 1, Ilevel - 1) ) );

    GOM_PARSER_VERIFY(
        gomp_ReturnSurfacePolygonData(
            Icontour - 1, Ilevel - 1, Ientry - 1) == 0 );

    GOM_PARSER_SUCCEEDED;
}

/* show contour polygon ... */
#define SHOW_CONTOUR_POLYGON_ENTRY \
GOM_PARSER_CMD_PART("polygon",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowContourPolygon,\
                    GOM_PARSER_UNUSED_VALUE)
static const gom_ParserArgumentList parseShowContourPolygon[] = {
    SHOW_CONTOUR_POLYGON_METHOD_ENTRY,
    SHOW_CONTOUR_POLYGON_ENTRIES_ENTRY,
    SHOW_CONTOUR_POLYGON_DATA_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};
    
const gom_ParserArgumentList gomp_ShowContourCommand[] = {
    SHOW_CONTOUR_DEFINED_ENTRY,
    SHOW_CONTOUR_FILE_NAME_ENTRY,
    SHOW_CONTOUR_STRUCTURE_ENTRY,
    SHOW_CONTOUR_MINIMUM_ENTRY,
    SHOW_CONTOUR_MAXIMUM_ENTRY,
    SHOW_CONTOUR_HEADER_ENTRY,
    SHOW_CONTOUR_LEVELS_ENTRY,
    SHOW_CONTOUR_VALUES_ENTRY,
    SHOW_CONTOUR_SMOOTH_ENTRY,
    SHOW_CONTOUR_CULLFACE_ENTRY,
    SHOW_CONTOUR_DISPLAY_ENTRY,
    SHOW_CONTOUR_CLIPPLANE_ENTRY,
    SHOW_CONTOUR_NAME_ENTRY,
    SHOW_CONTOUR_MAPPING_ENTRY,
    SHOW_CONTOUR_TYPE_ENTRY,
    SHOW_CONTOUR_TRANSPARENCY_ENTRY,
    SHOW_CONTOUR_CUBE_ENTRY,
    SHOW_CONTOUR_POLYGON_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

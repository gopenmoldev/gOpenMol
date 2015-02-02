/*

Copyright (c) 2003 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#include <stdlib.h>

#include "cluster.h"
#include "colouring.h"
#include "coord_file.h"
#include "gomenv.h"
#include "gomproc.h"
#include "gomtcl.h"
#include "gomversion.h"
#include "ldp.h"
#include "model_file.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "msd.h"
#include "objseg.h"
#include "parser.h"
#include "plot_molec.h"
#include "projview.h"
#include "rdf.h"
#include "rforce.h"
#include "selection.h"
#include "stereo.h"
#include "text_stack.h"
#include "trace.h"
#include "trajectory.h"

#include "stdafx.h"

#define ALIAS GOM_PARSER_CMD_ALIAS

int gomp_ParserReturnFloatValueFromFunc(GOM_PARSER_ARGLIST,ptrdiff_t Offset)
{
    FloatFunc *floatf = GOM_PARSER_GET_POINTER_VALUE(Offset);
    GOM_PARSER_RETURN_DOUBLE(floatf->Func());
    GOM_PARSER_SUCCEEDED;
}

int gomp_ParserReturnIntValueFromFunc(GOM_PARSER_ARGLIST,ptrdiff_t Offset)
{
    IntFunc *intf = GOM_PARSER_GET_POINTER_VALUE(Offset);
    GOM_PARSER_RETURN_INT(intf->Func());
    GOM_PARSER_SUCCEEDED;
}

int gomp_ParserReturnBooleanValueFromFunc(GOM_PARSER_ARGLIST,ptrdiff_t Offset)
{
    IntFunc *intf = GOM_PARSER_GET_POINTER_VALUE(Offset);
    GOM_PARSER_RETURN_BOOLEAN(intf->Func());
    GOM_PARSER_SUCCEEDED;
}

int gomp_ParserReturnBooleanEnumValueFromFunc(GOM_PARSER_ARGLIST,
                                              ptrdiff_t Offset)
{
    BooleanEnumFunc *booleanf = GOM_PARSER_GET_POINTER_VALUE(Offset);
    GOM_PARSER_RETURN_STRING(
        (booleanf->Func()) ?
        booleanf->strings.true_value :
        booleanf->strings.false_value);
    GOM_PARSER_SUCCEEDED;
}

int gomp_ParserReturnStringValueFromFunc(GOM_PARSER_ARGLIST,ptrdiff_t Offset)
{
    StringFunc *stringf = GOM_PARSER_GET_POINTER_VALUE(Offset);
    GOM_PARSER_RETURN_STRING(stringf->Func());
    GOM_PARSER_SUCCEEDED;
}

int gomp_ParserVerifyFuncResult(GOM_PARSER_ARGLIST,ptrdiff_t Offset)
{
    IntFunc *intf = GOM_PARSER_GET_POINTER_VALUE(Offset);
    GOM_PARSER_VERIFY( intf->Func() == 0 );
    GOM_PARSER_SUCCEEDED;
}

/* >>>>> atom >>>>> */
/* show atom ... */
#define SHOW_ATOM_ENTRY \
GOM_PARSER_CMD_PART("atom",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    gomp_ShowAtomCommand,\
                    GOM_PARSER_UNUSED_VALUE)

/* show atomnumber Seg Res Atm Istruct */
#define SHOW_ATOM_NUMBER_ENTRY \
GOM_PARSER_FINAL_CMD("atomnumber",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm Istruct?",0,4),\
                     ParseShowAtomNumber,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowAtomNumber(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    gom_SelectionList atom;
    int Istruct, IAtom;
    
    GOM_PARSER_START_ARGS
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_ARG(Istruct)
    GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY( gomp_HasMolecStructs() )

    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom selection list",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,atom));
    GOM_PARSER_VERIFY_GOM_PARSER(
        "structure index",
        GOM_PARSER_DEFAULT(Istruct,1) ||
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Istruct,1,gomp_GetNumMolecStructs()));

    IAtom = gomp_ShowAtomNumber(Istruct - 1, &atom);
    if ( IAtom < 0 )
        GOM_PARSER_FAILED;

    GOM_PARSER_RETURN_INT(IAtom);

    GOM_PARSER_SUCCEEDED;
}

/* show bondstyle */
#define SHOW_BOND_STYLE_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD("bondstyle",&parseShowBondStyle)
static const BooleanEnumFunc parseShowBondStyle = {
    gomp_GetBondDisplayStyle, { "smooth", "half" } };

/* show cylinderquality */
#define SHOW_CYLINDER_QUALITY_ENTRY \
GOM_PARSER_FINAL_INT_CMD("cylinderquality",&parseShowCylinderQuality)
static const IntFunc parseShowCylinderQuality = { gomp_GetCylinderQuality };

/* show spherequality */
#define SHOW_SPHERE_QUALITY_ENTRY \
GOM_PARSER_FINAL_INT_CMD("spherequality",&parseShowSphereQuality)
static const IntFunc parseShowSphereQuality = { gomp_GetSphereQuality };

typedef struct {
    float (*Func)(int);
} FloatFuncInt;

/* show licocylinder */
#define SHOW_LICO_CYLINDER_ENTRY \
GOM_PARSER_FINAL_CMD("licocylinder",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Istruct?",0,1),\
                     ParseShowLico,\
                     GOM_PARSER_SET_POINTER_VALUE(&parseShowLicoCylinder))
static const FloatFuncInt parseShowLicoCylinder = { gomp_GetAtomLicoRadC };

/* show licosphere */
#define SHOW_LICO_SPHERE_ENTRY \
GOM_PARSER_FINAL_CMD("licosphere",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Istruct?",0,1),\
                     ParseShowLico,\
                     GOM_PARSER_SET_POINTER_VALUE(&parseShowLicoSphere))
static const FloatFuncInt parseShowLicoSphere = { gomp_GetAtomLicoRadS };

static int ParseShowLico(GOM_PARSER_ARGLIST,ptrdiff_t Offset)
{
    int Istruct;
    FloatFuncInt *floatfint = GOM_PARSER_GET_POINTER_VALUE(Offset);
    
    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Istruct)
    GOM_PARSER_END_ARGS;
    
    GOM_PARSER_VERIFY( gomp_HasMolecStructs() );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "structure index",
        GOM_PARSER_DEFAULT(Istruct,1) ||
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Istruct,1,gomp_GetNumMolecStructs()));

    GOM_PARSER_RETURN_DOUBLE(floatfint->Func(Istruct - 1));

    GOM_PARSER_SUCCEEDED;
}

/* show molstructures */
#define SHOW_MOLSTRUCTURES_ENTRY \
GOM_PARSER_FINAL_INT_CMD("molstructures",&parseShowMolstructures)
static const IntFunc parseShowMolstructures = { gomp_GetNumMolecStructs };

/* show numatoms Istruct */
#define SHOW_NUM_ATOMS_ENTRY \
GOM_PARSER_FINAL_CMD("numatoms",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Istruct?",0,1),\
                     ParseShowNumAtoms,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowNumAtoms(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    int Istruct;
    
    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Istruct)
    GOM_PARSER_END_ARGS;
    
    GOM_PARSER_VERIFY( gomp_HasMolecStructs() );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "structure index",
        GOM_PARSER_DEFAULT(Istruct,1) ||
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Istruct,1,gomp_GetNumMolecStructs()));
    
    GOM_PARSER_RETURN_INT(gomp_GetNumAtomsInMolecStruct(Istruct - 1));

    GOM_PARSER_SUCCEEDED;
}

/* show residue stack */
#define SHOW_RESIDUE_STACK_ENTRY \
GOM_PARSER_FINAL_STRING_CMD("stack",&parseShowResidueStack)
static const StringFunc parseShowResidueStack = { gomp_ShowResidueNameStack };

/* show residue ... */
#define SHOW_RESIDUE_ENTRY \
GOM_PARSER_CMD_PART("residue",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowResidue,\
                    GOM_PARSER_UNUSED_VALUE)
static const gom_ParserArgumentList parseShowResidue[] = {
    SHOW_RESIDUE_STACK_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show segment stack */
#define SHOW_SEGMENT_STACK_ENTRY \
GOM_PARSER_FINAL_STRING_CMD("stack",&parseShowSegmentStack)
static const StringFunc parseShowSegmentStack = { gomp_ShowSegmentNameStack };

/* show segment ... */
#define SHOW_SEGMENT_ENTRY \
GOM_PARSER_CMD_PART("segment",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowSegment,\
                    GOM_PARSER_UNUSED_VALUE)
static const gom_ParserArgumentList parseShowSegment[] = {
    SHOW_SEGMENT_STACK_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};
/* <<<<< atom <<<<< */

/* >>>>> colour >>>>> */
/* show bgcolour */
/* show bgcolor  */
#define SHOW_BGCOLOUR_ENTRY \
GOM_PARSER_FINAL_CMD("bgcolour" ALIAS "bgcolor",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowBgColour,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowBgColour(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    float red, green, blue;
    (void)gomp_GetBGColor(&red, &green, &blue);
    
    GOM_PARSER_RETURN_LIST(("%f %f %f",red,green,blue));
    GOM_PARSER_SUCCEEDED;
}

/* show colourmapping */
/* show colormapping  */
#define SHOW_COLOUR_MAPPING_ENTRY \
GOM_PARSER_FINAL_CMD("colourmapping" ALIAS "colormapping",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowColourMapping,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowColourMapping(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    GOM_PARSER_RETURN_STRING(
        (gomp_GetColorMappingType() == COLOR_MAPPING_TYPE_TEXTURE) ?
        "texture" : "rainbow");
    GOM_PARSER_SUCCEEDED;
}

/* show colourtype */
/* show colortype  */
#define SHOW_COLOUR_TYPE_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD(\
    "colourtype" ALIAS "colortype",&parseShowColourType)
static const BooleanEnumFunc parseShowColourType = {
    gomp_GetDisplayColourType, { "gray", "color" } };
/* <<<<< colour <<<<< */


/* >>>>> directories >>>>> */
/* show bindirectory */
#define SHOW_BIN_DIRECTORY_ENTRY \
GOM_PARSER_FINAL_STRING_CMD("bindirectory",&parseShowBinDirectory)
static const StringFunc parseShowBinDirectory = { gomp_ShowBinDir };

/* show datadirectory */
#define SHOW_DATA_DIRECTORY_ENTRY \
GOM_PARSER_FINAL_STRING_CMD("datadirectory",&parseShowDataDirectory)
static const StringFunc parseShowDataDirectory = { gomp_ShowDataDir };

/* show helpdirectory */
#define SHOW_HELP_DIRECTORY_ENTRY \
GOM_PARSER_FINAL_STRING_CMD("helpdirectory",&parseShowHelpDirectory)
static const StringFunc parseShowHelpDirectory = { gomp_ShowHelpDir };

/* show homedirectory */
#define SHOW_HOME_DIRECTORY_ENTRY \
GOM_PARSER_FINAL_STRING_CMD("homedirectory",&parseShowHomeDirectory)
static const StringFunc parseShowHomeDirectory = { gomp_ShowHomeDir };
/* <<<<< directories <<<<< */

/* show cell ... */
#define SHOW_CELL_ENTRY \
GOM_PARSER_CMD_PART("cell",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    gomp_ShowCellCommand,\
                    GOM_PARSER_UNUSED_VALUE)

/* show cluster status */
#define SHOW_CLUSTER_STATUS_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_CMD("status",&parseShowClusterStatus)
static const IntFunc parseShowClusterStatus = { gomp_GetClusterStatus };

/* show cluster ... */
#define SHOW_CLUSTER_ENTRY \
GOM_PARSER_CMD_PART("cluster",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowCluster,\
                    GOM_PARSER_UNUSED_VALUE)
static const gom_ParserArgumentList parseShowCluster[] = {
    SHOW_CLUSTER_STATUS_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show cputime */
#define SHOW_CPUTIME_ENTRY \
GOM_PARSER_FINAL_CMD("cputime",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowCpuTime,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowCpuTime(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    float parent, child;
    gomp_Get_cpu_secs(&parent, &child);
    GOM_PARSER_RETURN_DOUBLE(parent);
    GOM_PARSER_SUCCEEDED;
}

/* show contour ... */
#define SHOW_CONTOUR_ENTRY \
GOM_PARSER_CMD_PART("contour",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    gomp_ShowContourCommand,\
                    GOM_PARSER_UNUSED_VALUE)

/* show cutplane */
#define SHOW_CUT_PLANE_ENTRY \
GOM_PARSER_CMD_PART("cutplane",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    gomp_ShowCutPlaneCommand,\
                    GOM_PARSER_UNUSED_VALUE)

/* show cscale */
#define SHOW_CSCALE_ENTRY \
GOM_PARSER_CMD_PART("cscale",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowColourScale,\
                    GOM_PARSER_UNUSED_VALUE)
static const IntFunc   parseShowColourScaleBins    = {
    gomp_GetColourScalePlotLevels };
static const FloatFunc parseShowColourScaleMinimum = {
    gomp_GetColourScalePlotMin };
static const FloatFunc parseShowColourScaleMaximum = {
    gomp_GetColourScalePlotMax };
static const IntFunc   parseShowColourScaleStatus  = {
    gomp_GetColourScalePlotStatus };
static const gom_ParserArgumentList parseShowColourScale[] = {
    GOM_PARSER_FINAL_INT_CMD    ("bins",   &parseShowColourScaleBins),
    GOM_PARSER_FINAL_FLOAT_CMD  ("minimum",&parseShowColourScaleMinimum),
    GOM_PARSER_FINAL_FLOAT_CMD  ("maximum",&parseShowColourScaleMaximum),
    GOM_PARSER_FINAL_BOOLEAN_CMD("state",  &parseShowColourScaleStatus),
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show displaylists */
#define SHOW_DISPLAYLISTS_ENTRY \
GOM_PARSER_CMD_PART("displaylists",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    gomp_ShowDisplaylistsCommand,\
                    GOM_PARSER_UNUSED_VALUE)

/* show drawbuffer */
#define SHOW_DRAWBUFFER_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD("drawbuffer",&parseShowDrawBuffer)
static const BooleanEnumFunc parseShowDrawBuffer = {
    gomp_GetDrawBuffer, { "back", "front" } };

/* >>>>> projection >>>>> */
/* show farplane {step|value} */
#define SHOW_FAR_PLANE_ENTRY \
GOM_PARSER_CMD_PART("farplane",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowFarPlane,\
                    GOM_PARSER_UNUSED_VALUE)
static const FloatFunc parseShowFarPlaneStep  = { gomp_GetPerspectiveFarStep };
static const FloatFunc parseShowFarPlaneValue = { gomp_GetPerspectiveFar };
static const gom_ParserArgumentList parseShowFarPlane[] = {
    GOM_PARSER_FINAL_FLOAT_CMD("step", &parseShowFarPlaneStep),
    GOM_PARSER_FINAL_FLOAT_CMD("value",&parseShowFarPlaneValue),
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show nearplane {step|value} */
#define SHOW_NEAR_PLANE_ENTRY \
GOM_PARSER_CMD_PART("nearplane",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowNearPlane,\
                    GOM_PARSER_UNUSED_VALUE)
static const FloatFunc parseShowNearPlaneStep  = {gomp_GetPerspectiveNearStep};
static const FloatFunc parseShowNearPlaneValue = { gomp_GetPerspectiveNear };
static const gom_ParserArgumentList parseShowNearPlane[] = {
    GOM_PARSER_FINAL_FLOAT_CMD("step", &parseShowNearPlaneStep),
    GOM_PARSER_FINAL_FLOAT_CMD("value",&parseShowNearPlaneValue),
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show viewangle */
#define SHOW_VIEW_ANGLE_ENTRY \
GOM_PARSER_FINAL_FLOAT_CMD("viewangle",&parseShowViewAngle)
static const FloatFunc parseShowViewAngle = { gomp_GetPerspectiveAngle };

/* show projection */
#define SHOW_PROJECTION_ENTRY \
GOM_PARSER_FINAL_INT_CMD("projection",&parseShowProjection)
static const IntFunc parseShowProjection = { gomp_GetProjectionTransformation};

/* show graphics */
#define SHOW_GRAPHICS_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_CMD("graphics",&parseShowGraphics)
static const IntFunc parseShowGraphics = { gomp_GetTermType };

/* show gbasis ... */
#define SHOW_GBASIS_ENTRY \
GOM_PARSER_CMD_PART("gbasis",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    gomp_ShowGbasisCommand,\
                    GOM_PARSER_UNUSED_VALUE)

/* show gromos96 coordamplifier */
#define SHOW_GROMOS96_COORDAMPLIFIER_ENTRY \
GOM_PARSER_FINAL_FLOAT_CMD("coordamplifier",&parseShowGromos96Coordamplifier)
static const FloatFunc parseShowGromos96Coordamplifier = {
    gomp_GetGROMOS96CoordAmplifier };

/* show gromos96 ... */
#define SHOW_GROMOS96_ENTRY \
GOM_PARSER_CMD_PART("gromos96",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowGromos96,\
                    GOM_PARSER_UNUSED_VALUE)
static const gom_ParserArgumentList parseShowGromos96[] = {
    SHOW_GROMOS96_COORDAMPLIFIER_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* tool to move text between c-code and tcl */
/* show gtext */
#define SHOW_GTEXT_ENTRY \
GOM_PARSER_FINAL_STRING_CMD("gtext",&parseShowGtext)
static const StringFunc parseShowGtext = { gomp_GetGlobalTextString };

/* show ldpstatus */
#define SHOW_LDPSTATUS_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD("ldpstatus",&parseShowLdpStatus)
static const BooleanEnumFunc parseShowLdpStatus = {
    gomp_GetDisplayLDPmatrix, { "off", "on" }
};

/* show light */
#define SHOW_LIGHT_ENTRY \
GOM_PARSER_CMD_PART("light",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    gomp_ShowLightCommand,\
                    GOM_PARSER_UNUSED_VALUE)

/* show material */
#define SHOW_MATERIAL_ENTRY \
GOM_PARSER_CMD_PART("material",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    gomp_ShowMaterialCommand,\
                    GOM_PARSER_UNUSED_VALUE)

/* show mlinewidth */
#define SHOW_MLINE_WIDTH_ENTRY \
GOM_PARSER_FINAL_INT_CMD("mlinewidth",&parseShowMlineWidth)
static const IntFunc parseShowMlineWidth = { gomp_GetMoleculeLineWidth };

/* show monitor ... */
#define SHOW_MONITOR_ENTRY \
GOM_PARSER_CMD_PART("monitor",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    gomp_ShowMonitorCommand,\
                    GOM_PARSER_UNUSED_VALUE)

/* show numframes */
#define SHOW_NUM_FRAMES_ENTRY \
GOM_PARSER_FINAL_INT_CMD("numframes",&parseShowNumFrames)
static const IntFunc parseShowNumFrames = { gomp_GetNumberOfFrames };

/* show os */
#define SHOW_OS_ENTRY \
GOM_PARSER_FINAL_STRING_CMD("os",&parseShowOS)
static const StringFunc parseShowOS = { gomp_GetOS };

/* show output */
#define SHOW_OUTPUT_ENTRY \
GOM_PARSER_FINAL_VERIFY_FUNC_CMD("output",&parseShowOutput)
static const IntFunc parseShowOutput = { gomp_SetOutput2Widget };

/* show plumber */
#define SHOW_PLUMBER_ENTRY \
GOM_PARSER_CMD_PART("plumber",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    gomp_ShowPlumberCommand,\
                    GOM_PARSER_UNUSED_VALUE)

/* show process */
#define SHOW_PROCESS_ENTRY \
GOM_PARSER_FINAL_VERIFY_FUNC_CMD("process",&parseShowProcess)
static const IntFunc parseShowProcess = { gomp_Get_proc_info };

#ifdef ENABLE_GRAPHICS
/* show quadstereo ... */
#define SHOW_QUADSTEREO_ENTRY \
GOM_PARSER_CMD_PART("quadstereo",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowQuadStereo,\
                    GOM_PARSER_UNUSED_VALUE)
static const BooleanEnumFunc parseShowQuadStereoState     = {
    gomp_QuadStereoIsOn, { "off", "on" } };
static const FloatFunc       parseShowQuadStereoAngle     = {
    gomp_GetQuadStereoHalfAngle };
static const IntFunc         parseShowQuadStereoAvailable = {
    gomp_CheckHardwareStereo };
static const gom_ParserArgumentList parseShowQuadStereo[] = {
    GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD("state",&parseShowQuadStereoState),
    GOM_PARSER_FINAL_FLOAT_CMD("angle",&parseShowQuadStereoAngle),
    GOM_PARSER_FINAL_BOOLEAN_CMD("available",&parseShowQuadStereoAvailable),
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show spair ... */
#define SHOW_STEREO_PAIR_ENTRY \
GOM_PARSER_CMD_PART("spair",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowStereoPair,\
                    GOM_PARSER_UNUSED_VALUE)
static const BooleanEnumFunc parseShowStereoPairState    = {
    gomp_GetStereoPlotState, { "off", "on" } };
static const FloatFunc       parseShowStereoPairDistance = {
    gomp_GetStereoPlotTranslate };
static const FloatFunc       parseShowStereoPairAngle    = {
    gomp_GetStereoPlotAngle };
static const gom_ParserArgumentList parseShowStereoPair[] = {
    GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD("state",&parseShowStereoPairState),
    GOM_PARSER_FINAL_FLOAT_CMD("distance",&parseShowStereoPairDistance),
    GOM_PARSER_FINAL_FLOAT_CMD("angle",&parseShowStereoPairAngle),
    GOM_PARSER_END_ARGUMENT_LIST
};  
#endif /* ENABLE_GRAPHICS */

/* show rainbow Fvalue */
#define SHOW_RAINBOW_ENTRY \
GOM_PARSER_FINAL_CMD("rainbow",\
                     GOM_PARSER_NEED_ARGS("Fvalue",1),\
                     ParseShowRainbow,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowRainbow(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    double Fvalue;
    float red, green, blue;
    
    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Fvalue)
    GOM_PARSER_END_ARGS;
    
    GOM_PARSER_VERIFY_GOM_PARSER(
        "scale factor",
        GOM_PARSER_RETRIEVE_DOUBLE_CHECK_RANGE(Fvalue, 0.0, 1.0) );
    
    gomp_PreRainbow(Fvalue, &red, &green, &blue);
    GOM_PARSER_RETURN_LIST(("%f %f %f",red,green,blue));
    GOM_PARSER_SUCCEEDED;
}

/* show rdf ... */
#define SHOW_RDF_ENTRY \
GOM_PARSER_CMD_PART("rdf",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowRDF,\
                    GOM_PARSER_UNUSED_VALUE)
static const IntFunc parseShowRDFstatus       = { gomp_GetRDFstatus };
static const IntFunc parseShowRDFobservations = { gomp_GetRDFobservations };
static const IntFunc parseShowRDFpoints       = { gomp_GetNumRDFObs };
static const IntFunc parseShowRDFaverage      = { gomp_GetRDFaverageStatus };
static const gom_ParserArgumentList parseShowRDF[] = {
    GOM_PARSER_FINAL_BOOLEAN_CMD("status",  &parseShowRDFstatus),
    GOM_PARSER_FINAL_INT_CMD("observations",&parseShowRDFobservations),
    GOM_PARSER_FINAL_INT_CMD("points",      &parseShowRDFpoints),
    GOM_PARSER_FINAL_INT_CMD("average",     &parseShowRDFaverage),
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show redisplay */
#define SHOW_REDISPLAY_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD("redisplay",&parseShowRedisplay)
static const BooleanEnumFunc parseShowRedisplay = {
    gomp_GetSystemRedisplayMode, { "fast", "slow" } };

/* show release */
#define SHOW_RELEASE_ENTRY \
GOM_PARSER_FINAL_CMD("release",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowRelease,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowRelease(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    GOM_PARSER_RETURN_STRING(GOPENMOL_RELEASE);
    GOM_PARSER_SUCCEEDED;
}

/* show rmsfluctuation length */
#define SHOW_RMS_FLUCTUATION_LENGTH_ENTRY \
GOM_PARSER_FINAL_INT_CMD("length",&parseShowRMSfluctuationLength)
static const IntFunc parseShowRMSfluctuationLength = { gomp_GetRMSlength };

/* show rmsfluctuation atomindex Irmsfatomindex */
#define SHOW_RMS_FLUCTUATION_ATOMINDEX_ENTRY \
GOM_PARSER_FINAL_CMD("atomindex",\
                     GOM_PARSER_NEED_ARGS("Irmsfatomindex",1),\
                     ParseShowRMSfluctuationAtomIndex,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowRMSfluctuationAtomIndex(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    int Irmsfatomindex, Iatom;
    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Irmsfatomindex)
    GOM_PARSER_END_ARGS;
    GOM_PARSER_VERIFY_GOM_PARSER(
        "rms fluctuation atom index",
        GOM_PARSER_RETRIEVE_INT(Irmsfatomindex) );
    GOM_PARSER_VERIFY( ( Iatom = gomp_GetRMSatomindex(Irmsfatomindex) ) > 0 );
    GOM_PARSER_RETURN_INT( Iatom );
    GOM_PARSER_SUCCEEDED;
}

/* show rmsfluctuation Irmsfatomindex */
#define SHOW_RMS_FLUCTUATION_ARGS_ENTRY \
GOM_PARSER_FINAL_ARGS(GOM_PARSER_NEED_ARGS("Irmsfatomindex",1),\
                      ParseShowRMSfluctuationArgs,\
                      GOM_PARSER_UNUSED_VALUE)
static int ParseShowRMSfluctuationArgs(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    int Irmsfatomindex;
    const float *FTemp_p;
    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Irmsfatomindex)
    GOM_PARSER_END_ARGS;
    GOM_PARSER_VERIFY_GOM_PARSER(
        "rms fluctuation atom index",
        GOM_PARSER_RETRIEVE_INT(Irmsfatomindex) );
    GOM_PARSER_VERIFY( ( FTemp_p = gomp_GetRMSfluctuation(Irmsfatomindex) ) );
    GOM_PARSER_RETURN_LIST(("%f %f %f %f",
                            FTemp_p[0], FTemp_p[1], FTemp_p[2], FTemp_p[3]));
    GOM_PARSER_SUCCEEDED;
}

/* show rmsfluctuation */
#define SHOW_RMS_FLUCTUATION_ENTRY \
GOM_PARSER_CMD_PART("rmsfluctuation",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowRmsFluctuation,\
                    GOM_PARSER_UNUSED_VALUE)
static const gom_ParserArgumentList parseShowRmsFluctuation[] = {
    SHOW_RMS_FLUCTUATION_LENGTH_ENTRY,
    SHOW_RMS_FLUCTUATION_ATOMINDEX_ENTRY,
    SHOW_RMS_FLUCTUATION_ARGS_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* >>>>> transformation >>>>> */
/* show rotation state */
#define SHOW_ROTATION_ENTRY \
GOM_PARSER_CMD_PART("rotation",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowRotation,\
                    GOM_PARSER_UNUSED_VALUE)
static const IntFunc parserShowRotationState = { gomp_GetRotationState };
static const gom_ParserArgumentList parseShowRotation[] = {
    GOM_PARSER_FINAL_BOOLEAN_CMD("state",&parserShowRotationState),
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show scaling type */
#define SHOW_SCALING_ENTRY \
GOM_PARSER_CMD_PART("scaling",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowScaling,\
                    GOM_PARSER_UNUSED_VALUE)
static const BooleanEnumFunc parseShowScalingType = {
    gomp_GetAllowIndividualScaling, { "global", "local" } };
static const gom_ParserArgumentList parseShowScaling[] = {
    GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD("type",&parseShowScalingType),
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show system translation */
#define SHOW_SYSTEM_ENTRY \
GOM_PARSER_CMD_PART("system",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowSystem,\
                    GOM_PARSER_UNUSED_VALUE)
static const BooleanEnumFunc parseShowSystemTranslation = {
    gomp_GetSystemTranslateState, { "on", "off" } };
static const gom_ParserArgumentList parseShowSystem[] = {
    GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD("translation",
                                            &parseShowSystemTranslation),
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show transformation type */
#define SHOW_TRANSFORMATION_ENTRY \
GOM_PARSER_CMD_PART("transformation",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowTransformation,\
                    GOM_PARSER_UNUSED_VALUE)
static const BooleanEnumFunc parseShowTransformationType = {
    gomp_GetObjectCenterType, { "global", "local" } };
static const gom_ParserArgumentList parseShowTransformation[] = {
    GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD("type",
                                            &parseShowTransformationType),
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show translation array  */
/* show translation vector */
#define SHOW_TRANSLATION_ARRAY_ENTRY \
GOM_PARSER_FINAL_CMD("array",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowTranslationArray,\
                     GOM_PARSER_UNUSED_VALUE),\
GOM_PARSER_FINAL_CMD("vector",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowTranslationArray,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowTranslationArray(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    const float *sumxyz = gomp_GetTranslateArray();
    GOM_PARSER_RETURN_LIST(("%f %f %f",sumxyz[0],sumxyz[1],sumxyz[2]));
    GOM_PARSER_SUCCEEDED;
}

/* show translation state */
#define SHOW_TRANSLATION_STATE_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_CMD("state",&parseShowTranslationState)
static const IntFunc parseShowTranslationState = { gomp_GetTranslationState };

/* show translation ... */
#define SHOW_TRANSLATION_ENTRY \
GOM_PARSER_CMD_PART("translation",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowTranslation,\
                    GOM_PARSER_UNUSED_VALUE)
static const gom_ParserArgumentList parseShowTranslation[] = {
    SHOW_TRANSLATION_ARRAY_ENTRY,
    SHOW_TRANSLATION_STATE_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};
/* <<<<< transformation <<<<< */

/* show savestatus */
#define SHOW_SAVE_STATUS_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_CMD("savestatus",&parseShowSaveStatus)
static const IntFunc parseShowSaveStatus = { gomp_gOpenMolNeedsSaving };

/* show selection */
#define SHOW_SELECTION_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD("selection",&parseShowSelection)
static const BooleanEnumFunc parseShowSelection = {
    gomp_GetSelectionModeStatus, { "off", "on" } };

/* show sizeofsystem */
#define SHOW_SIZE_OF_SYSTEM_ENTRY \
GOM_PARSER_FINAL_CMD("sizeofsystem",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowSizeOfSystem,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowSizeOfSystem(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    GOM_PARSER_VERIFY( gomp_HasMolecStructs() );
    GOM_PARSER_RETURN_DOUBLE( gomp_GetSizeOfSystem() );
    GOM_PARSER_SUCCEEDED;
}


/* show status */
#define SHOW_STATUS_ENTRY \
GOM_PARSER_FINAL_CMD("status",\
                     GOM_PARSER_NEED_ARGS("Istruct",1),\
                     ParseShowStatus,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowStatus(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    int Istruct;
    
    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Istruct)
        GOM_PARSER_END_ARGS;
    
    GOM_PARSER_VERIFY( gomp_HasMolecStructs() );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "molecule structure index",
        GOM_PARSER_DEFAULT(Istruct,1) ||
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Istruct,1,gomp_GetNumMolecStructs()));
    GOM_PARSER_RETURN_INT(gomp_GetSelectedStructure(Istruct - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show time */
#define SHOW_TIME_ENTRY \
GOM_PARSER_FINAL_CMD("time",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowTime,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowTime(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    GOM_PARSER_VERIFY(gomp_Get_date(3));
    GOM_PARSER_SUCCEEDED;
}

/* show traces  */
/* show ptraces */
#define SHOW_TRACES_ENTRY \
GOM_PARSER_FINAL_INT_CMD("traces",&parseShowTraces),\
GOM_PARSER_FINAL_INT_CMD("ptraces",&parseShowTraces)
static const IntFunc parseShowTraces = { gomp_GetTraceSets };

/* show trajectory display */
#define SHOW_TRAJECTORY_DISPLAY_ENTRY \
GOM_PARSER_FINAL_CMD("display",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowTrajectoryDisplay,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowTrajectoryDisplay(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    int FirstFrame, LastFrame, FrameStep;
    (void)gomp_GetTrajectoryDisplayParams(&FirstFrame,&LastFrame,&FrameStep);
    GOM_PARSER_RETURN_LIST(("%d %d %d",FirstFrame,LastFrame,FrameStep));
    GOM_PARSER_SUCCEEDED;
}

/* show trajectory ... */
#define SHOW_TRAJECTORY_ENTRY \
GOM_PARSER_CMD_PART("trajectory",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowTrajectory,\
                    GOM_PARSER_UNUSED_VALUE)
static const IntFunc    parseShowTrajectoryDefined  = {
    gomp_GetTrajectoryStructureState };
static const IntFunc    parseShowTrajectoryFrames   = {
    gomp_GetNumberOfFrames };
static const StringFunc parseShowTrajectoryFileName = {
    gomp_GetTrajectoryFileName };
static const IntFunc    parseShowTrajectoryAction   = {
    gomp_GetFormattedTrajectoryReader };
static const IntFunc    parseShowTrajectoryCurrent  = {
    gomp_GetDisplayFrameNumber };
static const IntFunc    parseShowTrajectoryFid      = {
    gomp_GetDisplayRunningFrameNumberState };
static const StringFunc parseShowTrajectoryType     = {
    gomp_GetTrajectoryFileTypeName };
static gom_ParserArgumentList parseShowTrajectory[] = {
    GOM_PARSER_FINAL_BOOLEAN_CMD("defined", &parseShowTrajectoryDefined),
    GOM_PARSER_FINAL_INT_CMD    ("frames",  &parseShowTrajectoryFrames),
    GOM_PARSER_FINAL_STRING_CMD ("filename",&parseShowTrajectoryFileName),
    GOM_PARSER_FINAL_INT_CMD    ("action",  &parseShowTrajectoryAction),
    GOM_PARSER_FINAL_INT_CMD    ("current", &parseShowTrajectoryCurrent),
    GOM_PARSER_FINAL_BOOLEAN_CMD("fid",     &parseShowTrajectoryFid),
    GOM_PARSER_FINAL_STRING_CMD ("type",    &parseShowTrajectoryType),
    SHOW_TRAJECTORY_DISPLAY_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show vector */
#define SHOW_VECTOR_ENTRY \
GOM_PARSER_FINAL_CMD("vector",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowVector,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowVector(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    float min,max;
    (void)gomp_GetVectorDisplayRange(&min,&max);
    GOM_PARSER_RETURN_LIST(("%f %f %f %f",
                            gomp_GetVectorScale(),
                            gomp_GetVectorRadius(),
                            min,max));
    GOM_PARSER_SUCCEEDED;
}

/* show version */
#define SHOW_VERSION_ENTRY \
GOM_PARSER_FINAL_CMD("version",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowVersion,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowVersion(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    GOM_PARSER_RETURN_INT(GOPENMOL_VERSION);
    GOM_PARSER_SUCCEEDED;
}

#ifdef ENABLE_GRAPHICS
/* show window */
#define SHOW_WINDOW_ENTRY \
GOM_PARSER_CMD_PART("window",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    gomp_ShowWindowCommand,\
                    GOM_PARSER_UNUSED_VALUE)
#endif

const gom_ParserArgumentList gomp_ShowCommand[] = {
    /* atom */
    SHOW_ATOM_ENTRY,
    SHOW_ATOM_NUMBER_ENTRY,
    SHOW_BOND_STYLE_ENTRY,
    SHOW_CYLINDER_QUALITY_ENTRY,
    SHOW_SPHERE_QUALITY_ENTRY,
    SHOW_LICO_CYLINDER_ENTRY,
    SHOW_LICO_SPHERE_ENTRY,
    SHOW_RESIDUE_ENTRY,
    SHOW_SEGMENT_ENTRY,

    /* structure */
    SHOW_MOLSTRUCTURES_ENTRY,
    SHOW_NUM_ATOMS_ENTRY,

    /* colour */
    SHOW_BGCOLOUR_ENTRY,
    SHOW_COLOUR_MAPPING_ENTRY,
    SHOW_COLOUR_TYPE_ENTRY,

    /* directories */
    SHOW_BIN_DIRECTORY_ENTRY,
    SHOW_DATA_DIRECTORY_ENTRY,
    SHOW_HELP_DIRECTORY_ENTRY,
    SHOW_HOME_DIRECTORY_ENTRY,

    SHOW_CELL_ENTRY,
    SHOW_CLUSTER_ENTRY,
    SHOW_CPUTIME_ENTRY,

    /* contour */
    SHOW_CONTOUR_ENTRY,
    SHOW_CUT_PLANE_ENTRY,

    SHOW_CSCALE_ENTRY,

    SHOW_DISPLAYLISTS_ENTRY,
    SHOW_DRAWBUFFER_ENTRY,

    /* projection */
    SHOW_FAR_PLANE_ENTRY,
    SHOW_NEAR_PLANE_ENTRY,
    SHOW_VIEW_ANGLE_ENTRY,
    SHOW_PROJECTION_ENTRY,

    SHOW_GBASIS_ENTRY,
    SHOW_GROMOS96_ENTRY,
    SHOW_GTEXT_ENTRY,
    SHOW_GRAPHICS_ENTRY,
    SHOW_LDPSTATUS_ENTRY,
    SHOW_LIGHT_ENTRY,
    SHOW_MATERIAL_ENTRY,
    SHOW_MLINE_WIDTH_ENTRY,
    SHOW_MONITOR_ENTRY,
    SHOW_NUM_FRAMES_ENTRY,
    SHOW_OS_ENTRY,
    SHOW_OUTPUT_ENTRY,
    SHOW_PLUMBER_ENTRY,
    SHOW_PROCESS_ENTRY,

    /* stereo */
#ifdef ENABLE_GRAPHICS
    SHOW_QUADSTEREO_ENTRY,
    SHOW_STEREO_PAIR_ENTRY,
#endif

    SHOW_RAINBOW_ENTRY,
    SHOW_RDF_ENTRY,
    SHOW_REDISPLAY_ENTRY,
    SHOW_RELEASE_ENTRY,
    SHOW_RMS_FLUCTUATION_ENTRY,

    /* transformation */
    SHOW_ROTATION_ENTRY,
    SHOW_SCALING_ENTRY,
    SHOW_SYSTEM_ENTRY,
    SHOW_TRANSFORMATION_ENTRY,
    SHOW_TRANSLATION_ENTRY,

    SHOW_SAVE_STATUS_ENTRY,

    SHOW_SELECTION_ENTRY,
    SHOW_SIZE_OF_SYSTEM_ENTRY,
    SHOW_STATUS_ENTRY,
    SHOW_TIME_ENTRY,
    SHOW_TRACES_ENTRY,
    SHOW_TRAJECTORY_ENTRY,
    SHOW_VECTOR_ENTRY,
    SHOW_VERSION_ENTRY,
#ifdef ENABLE_GRAPHICS
    SHOW_WINDOW_ENTRY,
#endif
    GOM_PARSER_END_ARGUMENT_LIST
};

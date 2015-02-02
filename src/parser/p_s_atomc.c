/*

Copyright (c) 2003 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#include "bond.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "parser.h"
#include "picking.h"
#include "plot_molec.h"
#include "projview.h"
#include "selection.h"
#include "trajectory.h"

#include "stdafx.h"

#define ALIAS GOM_PARSER_CMD_ALIAS

/* ... ?Istruct? */
static int ParseGetStructureIndex(int *pIstruct,
                                  GOM_PARSER_ARGLIST)
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

    *pIstruct = Istruct;

    GOM_PARSER_SUCCEEDED;
}

#define GET_STRUCTURE_INDEX \
    int Iatom, Istruct; \
    GOM_PARSER_VERIFY_CHILD( \
        ParseGetStructureIndex(&Istruct,GOM_PARSER_PASS_ARGS))

/* ... Iatom ?Istruct? */
static int ParseGetAtomAndStructureIndeces(int *pIatom, int *pIstruct,
                                           GOM_PARSER_ARGLIST)
{
    int Iatom, Istruct;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Iatom)
        GOM_PARSER_ARG(Istruct)
    GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY( gomp_HasMolecStructs() );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "structure index",
        GOM_PARSER_DEFAULT(Istruct,1) ||
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Istruct,1,gomp_GetNumMolecStructs()));
    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom index",
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Iatom,1,gomp_GetNumAtomsInMolecStruct(Istruct - 1)));

    *pIatom   = Iatom;
    *pIstruct = Istruct;

    GOM_PARSER_SUCCEEDED;
}

#define GET_ATOM_AND_STRUCTURE_INDECES \
    int Iatom, Istruct; \
    GOM_PARSER_VERIFY_CHILD( \
        ParseGetAtomAndStructureIndeces(&Iatom,&Istruct,GOM_PARSER_PASS_ARGS))

/* show atom coordinates Iatom Istruct */
#define SHOW_ATOM_COORDINATES_ENTRY \
GOM_PARSER_FINAL_CMD("coordinates",\
                     GOM_PARSER_NEED_ARGS_RANGE("Iatom ?Istruct?",1,2),\
                     ParseShowAtomCoordinates,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowAtomCoordinates(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    const float *sumxyz = gomp_GetTranslateArray();
    GET_ATOM_AND_STRUCTURE_INDECES;
    GOM_PARSER_RETURN_LIST(
        ("%f %f %f",
         gomp_GetAtomXCoord( Istruct - 1, Iatom - 1 ) + sumxyz[0],
         gomp_GetAtomYCoord( Istruct - 1, Iatom - 1 ) + sumxyz[1],
         gomp_GetAtomZCoord( Istruct - 1, Iatom - 1 ) + sumxyz[2] ) );
    GOM_PARSER_SUCCEEDED;
}

/* show atom colour Iatom Istruct */
/* show atom color  Iatom Istruct */
#define SHOW_ATOM_COLOUR_ENTRY \
GOM_PARSER_FINAL_CMD("colour" ALIAS "color",\
                     GOM_PARSER_NEED_ARGS_RANGE("Iatom ?Istruct?",1,2),\
                     ParseShowAtomColour,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowAtomColour(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    const float *red,*green,*blue;
    GET_ATOM_AND_STRUCTURE_INDECES;
    red   = gomp_GetAtomColourRedPointer  (Istruct - 1);
    green = gomp_GetAtomColourGreenPointer(Istruct - 1);
    blue  = gomp_GetAtomColourBluePointer (Istruct - 1);
    GOM_PARSER_RETURN_LIST(("%f %f %f",
                            red  [Iatom - 1],
                            green[Iatom - 1],
                            blue [Iatom - 1]));
    GOM_PARSER_SUCCEEDED;
}

/* show atom force Iatom */
#define SHOW_ATOM_FORCE_ENTRY \
GOM_PARSER_FINAL_CMD("force",\
                     GOM_PARSER_NEED_ARGS("Iatom",1),\
                     ParseShowAtomForce,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowAtomForce(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    int Iatom;
    const float *Xf,*Yf,*Zf;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Iatom)
    GOM_PARSER_END_ARGS;

    if ( ! gomp_GetForceRetrieveReadyState() ) {
        gomp_PrintERROR("no atom force information available");
        GOM_PARSER_FAILED;
    }
    GOM_PARSER_VERIFY( gomp_HasMolecStructs() );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom index",
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Iatom,1,gomp_GetNumAtomsInMolecStruct(0)));

    Xf = gomp_GetForceXComponentPointer();
    Yf = gomp_GetForceYComponentPointer();
    Zf = gomp_GetForceZComponentPointer();

    GOM_PARSER_RETURN_LIST(("%f %f %f",
                            Xf[Iatom - 1], Yf[Iatom - 1], Zf[Iatom - 1]));
    GOM_PARSER_SUCCEEDED;
}

/* show atom velocity Iatom */
#define SHOW_ATOM_VELOCITY_ENTRY \
GOM_PARSER_FINAL_CMD("velocity",\
                     GOM_PARSER_NEED_ARGS("Iatom",1),\
                     ParseShowAtomVelocity,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowAtomVelocity(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    int Iatom;
    const float *Xv,*Yv,*Zv;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Iatom)
    GOM_PARSER_END_ARGS;

    if ( ! gomp_GetVelocityRetrieveReadyState() ) {
        gomp_PrintERROR("no atom velocity information available");
        GOM_PARSER_FAILED;
    }
    GOM_PARSER_VERIFY( gomp_HasMolecStructs() );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom index",
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Iatom,1,gomp_GetNumAtomsInMolecStruct(0)));

    Xv = gomp_GetVelocityXComponentPointer();
    Yv = gomp_GetVelocityYComponentPointer();
    Zv = gomp_GetVelocityZComponentPointer();

    GOM_PARSER_RETURN_LIST(("%f %f %f",
                            Xv[Iatom - 1], Yv[Iatom - 1], Zv[Iatom - 1]));
    GOM_PARSER_SUCCEEDED;
}

/* show atom cross */
#define SHOW_ATOM_CROSS_ENTRY \
GOM_PARSER_FINAL_FLOAT_CMD("cross",&parseShowAtomCross)
static const FloatFunc parseShowAtomCross = { gomp_GetCrossLen };

/* show atom indentify */
/* "indentify" corrected to "identify", 14.10.2003, LUL */
#define SHOW_ATOM_IDENTIFY_ENTRY \
GOM_PARSER_FINAL_BOOLEAN_CMD("identify",&parseShowAtomIdentify)
static const IntFunc parseShowAtomIdentify = { gomp_GetIdentifyAtomActive };

/* show atom selection */
#define SHOW_ATOM_SELECTION_ENTRY \
GOM_PARSER_FINAL_CMD("selection",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowAtomSelection,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowAtomSelection(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    switch(gomp_GetAtomSelectionMode()) {            
    case RESIDUE_SELECTION:
        GOM_PARSER_RETURN_STRING("residue");
        break;
    case ATOM_SELECTION:
        GOM_PARSER_RETURN_STRING("atom");
        break;
    }

    GOM_PARSER_SUCCEEDED;
}

/* show atom window */
#define SHOW_ATOM_WINDOW_ENTRY_ENTRY \
GOM_PARSER_FINAL_INT_CMD("window",&parseShowAtomWindow)
static const IntFunc parseShowAtomWindow = { gomp_GetConnectionSearchWindow };

/* show atom reconnectivity */
#define SHOW_ATOM_RECONNECTIVITY_ENTRY \
GOM_PARSER_FINAL_INT_CMD("reconnectivity",&parseShowAtomReconnectivity)
static const IntFunc parseShowAtomReconnectivity = {
    gomp_GetBondReconnectivityState };

/* show atom hbreconnectivity */
#define SHOW_ATOM_HBRECONNECTIVITY_ENTRY \
GOM_PARSER_FINAL_INT_CMD("hbreconnectivity",&parseShowAtomHbreconnectivity)
static const IntFunc parseShowAtomHbreconnectivity = {
    gomp_GetHBondReconnectivityState };

/* show atom maxconnectivity */
#define SHOW_ATOM_MAXCONNECTIVITY_ENTRY \
GOM_PARSER_FINAL_INT_CMD("maxconnectivity",&parseShowAtomMaxconnectivity)
static const IntFunc parseShowAtomMaxconnectivity = {
    gomp_GetMaxAtomConnections };

/* show atom info Elem */
#define SHOW_ATOM_INFO_ENTRY \
GOM_PARSER_FINAL_CMD("info",\
                     GOM_PARSER_NEED_ARGS("Elem",1),\
                     ParseShowAtomInfo,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowAtomInfo(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    const char *Elem;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Elem)
    GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "element type",
        GOM_PARSER_RETRIEVE_STRING(Elem) );

    GOM_PARSER_VERIFY( gomp_ParseGetAtomInfo(Elem) == 0 );

    GOM_PARSER_SUCCEEDED;
}

/* show atom stack */
#define SHOW_ATOM_STACK_ENTRY \
GOM_PARSER_FINAL_STRING_CMD("stack",&parseShowAtomStack)
static const StringFunc parseShowAtomStack = { gomp_ShowAtomNameStack };

/* show atom structure atoms Istruct */
#define SHOW_ATOM_STRUCT_ATOMS_ENTRY \
GOM_PARSER_FINAL_CMD("atoms",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Istruct?",0,1),\
                     ParseShowAtomStructAtoms,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowAtomStructAtoms(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_STRUCTURE_INDEX;
    GOM_PARSER_RETURN_INT(gomp_GetNumAtomsInMolecStruct(Istruct - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show atom structure defined */
#define SHOW_ATOM_STRUCT_DEFINED_ENTRY \
GOM_PARSER_FINAL_INT_CMD("defined",&parseShowAtomStructDefined)
static const IntFunc parseShowAtomStructDefined = { gomp_GetNumMolecStructs };

/* show atom structure filename Istruct */
#define SHOW_ATOM_STRUCT_FILE_NAME_ENTRY \
GOM_PARSER_FINAL_CMD("filename",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Istruct?",0,1),\
                     ParseShowAtomStructFileName,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowAtomStructFileName(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_STRUCTURE_INDEX;
    GOM_PARSER_RETURN_STRING(gomp_GetMolecStructFileName(Istruct - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show atom structure name Istruct */
#define SHOW_ATOM_STRUCT_NAME_ENTRY \
GOM_PARSER_FINAL_CMD("name",\
                     GOM_PARSER_NEED_ARGS_RANGE("?Istruct?",0,1),\
                     ParseShowAtomStructName,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowAtomStructName(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    GET_STRUCTURE_INDEX;
    GOM_PARSER_RETURN_STRING(gomp_GetMolecStructName(Istruct - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show atom structure ... */
#define SHOW_ATOM_STRUCT_ENTRY \
GOM_PARSER_CMD_PART("structure" ALIAS "structures",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowAtomStruct,\
                    GOM_PARSER_UNUSED_VALUE)
static const gom_ParserArgumentList parseShowAtomStruct[] = {
    SHOW_ATOM_STRUCT_ATOMS_ENTRY,
    SHOW_ATOM_STRUCT_DEFINED_ENTRY,
    SHOW_ATOM_STRUCT_FILE_NAME_ENTRY,
    SHOW_ATOM_STRUCT_NAME_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* Show string data. */
typedef struct {
    const char* (*Func)(int,int);
} AtomStringFunc;

static int ParseShowAtomStringData(GOM_PARSER_ARGLIST,intptr_t Ptr)
{
    AtomStringFunc *asf = GOM_PARSER_GET_POINTER_VALUE(Ptr);
    GET_ATOM_AND_STRUCTURE_INDECES;
    GOM_PARSER_RETURN_STRING(asf->Func(Istruct - 1, Iatom - 1));
    GOM_PARSER_SUCCEEDED;
}

/* Show integer data. */
typedef struct {
    int (*Func)(int,int);
} AtomIntFunc;

static int ParseShowAtomIntData(GOM_PARSER_ARGLIST,intptr_t Ptr)
{
    AtomIntFunc *aif = GOM_PARSER_GET_POINTER_VALUE(Ptr);
    GET_ATOM_AND_STRUCTURE_INDECES;
    GOM_PARSER_RETURN_INT(aif->Func(Istruct - 1, Iatom - 1));
    GOM_PARSER_SUCCEEDED;
}

/* Show state data. */
typedef struct {
    char (*Func)(int,int);
} AtomStateFunc;

static int ParseShowAtomStateData(GOM_PARSER_ARGLIST,intptr_t Ptr)
{
    AtomStateFunc *asf = GOM_PARSER_GET_POINTER_VALUE(Ptr);
    GET_ATOM_AND_STRUCTURE_INDECES;
    GOM_PARSER_RETURN_BOOLEAN(asf->Func(Istruct - 1, Iatom - 1));
    GOM_PARSER_SUCCEEDED;
}

/* Show character as a string. */
typedef struct {
    char (*Func)(int,int);
} AtomCharFunc;

static int ParseShowAtomCharData(GOM_PARSER_ARGLIST,intptr_t Ptr)
{
    AtomCharFunc *acf = GOM_PARSER_GET_POINTER_VALUE(Ptr);
    char result[2] = {'\0','\0'};
    GET_ATOM_AND_STRUCTURE_INDECES;
    result[0] = acf->Func(Istruct - 1, Iatom - 1);
    GOM_PARSER_RETURN_STRING(result);
    GOM_PARSER_SUCCEEDED;
}

/* Show real data. */
typedef struct {
    float (*Func)(int,int);
} AtomFloatFunc;

static int ParseShowAtomFloatData(GOM_PARSER_ARGLIST,intptr_t Ptr)
{
    AtomFloatFunc *aff = GOM_PARSER_GET_POINTER_VALUE(Ptr);
    GET_ATOM_AND_STRUCTURE_INDECES;
    GOM_PARSER_RETURN_DOUBLE(aff->Func(Istruct - 1, Iatom - 1));
    GOM_PARSER_SUCCEEDED;
}

/* Show index list. */
typedef struct {
    const int* (*Func)(int,int);
} AtomConnFunc;

static int ParseShowAtomConnData(GOM_PARSER_ARGLIST,intptr_t Ptr)
{
    AtomConnFunc *acf = GOM_PARSER_GET_POINTER_VALUE(Ptr);
    const int *conn;
    int i;
    gom_ParserList list;
    /* Retrieve indeces. */
    GET_ATOM_AND_STRUCTURE_INDECES;

    conn = acf->Func(Istruct - 1, Iatom - 1);
    if ( ! GOM_PARSER_LIST_INIT(list) ) {
        GOM_PARSER_FAILED;
    }
    /* Append element count and all valid elements.
     * conn[0] is element count.
     * 1 ... conn[0] are valid elements.
     */
    if ( ! GOM_PARSER_LIST_APPEND_INT(list,conn[0]) ) {
        GOM_PARSER_LIST_FREE(list);
        GOM_PARSER_FAILED;
    }
    for ( i = 1 ; i <= conn[0] ; i++ ) {
        if ( ! GOM_PARSER_LIST_APPEND_INT(list,conn[i]+1) ) {
            GOM_PARSER_LIST_FREE(list);
            GOM_PARSER_FAILED;
        }
    }
    GOM_PARSER_LIST_RETURN(list);
    GOM_PARSER_SUCCEEDED;
}

static const AtomStringFunc GetAtomSegName     = { gomp_GetAtomSegName };
static const AtomStringFunc GetAtomResName     = { gomp_GetAtomResName };
static const AtomStringFunc GetAtomAtmName     = { gomp_GetAtomAtmName };
static const AtomIntFunc    GetAtomResNum      = { gomp_GetAtomResNum1 };
static const AtomFloatFunc  GetAtomCharge      = { gomp_GetAtomCharge };
static const AtomFloatFunc  GetAtomNucCharge   = { gomp_GetAtomNucCharge };
static const AtomFloatFunc  GetAtomCovar       = { gomp_GetAtomCovar };
static const AtomFloatFunc  GetAtomBValue      = { gomp_GetAtomBValue };
static const AtomStringFunc GetAtomBasisSetTag = { gomp_GetAtomBasisSetTag };
static const AtomIntFunc    GetAtomType        = { gomp_GetAtomType };
static const AtomFloatFunc  GetAtomBndRad      = { gomp_GetAtomBndRad };
static const AtomFloatFunc  GetAtomVdwRad      = { gomp_GetAtomVdwRad };
static const AtomFloatFunc  GetAtomPluRad      = { gomp_GetAtomPluRad };
static const AtomCharFunc   GetAtomGlobal      = { gomp_GetAtomGlobal };
static const AtomFloatFunc  GetAtomEmin        = { gomp_GetAtomEmin };
static const AtomFloatFunc  GetAtomRmin        = { gomp_GetAtomRmin };
static const AtomFloatFunc  GetAtomPatom       = { gomp_GetAtomPatom };
static const AtomFloatFunc  GetAtomMass        = { gomp_GetAtomMass };
static const AtomIntFunc    GetAtomCnct        = { gomp_GetAtomCnct };
static const AtomCharFunc   GetAtomHbond       = { gomp_GetAtomHbond };
static const AtomStringFunc GetAtomAtype       = { gomp_GetAtomAtype };
static const AtomStateFunc  GetAtomDisplayState     = { gomp_GetAtomDisplayState };
static const AtomStateFunc  GetAtomCPKDisplayState  = {
    gomp_GetAtomCPKDisplayState };
static const AtomFloatFunc  GetAtomCPKScale         = {
    gomp_GetAtomCPKScale };
static const AtomStateFunc  GetAtomLicoDisplayState = {
    gomp_GetAtomLicoDisplayState };
static const AtomConnFunc   GetAtomConnection   = { gomp_GetAtomConnection };
static const AtomConnFunc   GetAtomHydrogenBond = { gomp_GetAtomHydrogenBond };

#define SHOW_ATOM_DATA_ENTRY(Type,cmd,pointer) \
GOM_PARSER_FINAL_CMD(cmd,\
                     GOM_PARSER_NEED_ARGS_RANGE("Iatom ?Istruct?",1,2),\
                     ParseShowAtom##Type##Data,\
                     GOM_PARSER_SET_POINTER_VALUE(pointer))

const gom_ParserArgumentList gomp_ShowAtomCommand[] = {
    /* atom specific info */
    SHOW_ATOM_DATA_ENTRY(String, "segmentname",    &GetAtomSegName),
    SHOW_ATOM_DATA_ENTRY(String, "residuename",    &GetAtomResName),
    SHOW_ATOM_DATA_ENTRY(String, "atomname",       &GetAtomAtmName),
    SHOW_ATOM_DATA_ENTRY(Int   , "resnumber",      &GetAtomResNum),
    SHOW_ATOM_DATA_ENTRY(Float , "charge",         &GetAtomCharge),
    SHOW_ATOM_DATA_ENTRY(Float , "nuclearcharge",  &GetAtomNucCharge),
    SHOW_ATOM_DATA_ENTRY(Float , "covar",          &GetAtomCovar),
    SHOW_ATOM_DATA_ENTRY(Float , "bvalue",         &GetAtomBValue),
    SHOW_ATOM_DATA_ENTRY(String, "gbasis",         &GetAtomBasisSetTag),
    SHOW_ATOM_DATA_ENTRY(Int   , "type",           &GetAtomType),
    SHOW_ATOM_DATA_ENTRY(Float , "bnd",            &GetAtomBndRad),
    SHOW_ATOM_DATA_ENTRY(Float , "vdw",            &GetAtomVdwRad),
    SHOW_ATOM_DATA_ENTRY(Float , "plu",            &GetAtomPluRad),
    SHOW_ATOM_DATA_ENTRY(Char  , "global",         &GetAtomGlobal),
    SHOW_ATOM_DATA_ENTRY(Float , "rmin",           &GetAtomRmin),
    SHOW_ATOM_DATA_ENTRY(Float , "emin",           &GetAtomEmin),
    SHOW_ATOM_DATA_ENTRY(Float , "patom",          &GetAtomPatom),
    SHOW_ATOM_DATA_ENTRY(Float , "mass",           &GetAtomMass),
    SHOW_ATOM_DATA_ENTRY(Float , "amass",          &GetAtomMass), /* ALIAS */
    SHOW_ATOM_DATA_ENTRY(Int   , "cnct",           &GetAtomCnct),
    SHOW_ATOM_DATA_ENTRY(Char  , "hbond",          &GetAtomHbond),
    SHOW_ATOM_DATA_ENTRY(String, "atype",          &GetAtomAtype),
    SHOW_ATOM_DATA_ENTRY(State , "displaystate",   &GetAtomDisplayState),
    SHOW_ATOM_DATA_ENTRY(State , "cpkstate" ALIAS "cpkdisplay",
                                                   &GetAtomCPKDisplayState),
    SHOW_ATOM_DATA_ENTRY(Float , "cpkscale",       &GetAtomCPKScale),
    SHOW_ATOM_DATA_ENTRY(State , "licoricestate",  &GetAtomLicoDisplayState),
    SHOW_ATOM_DATA_ENTRY(Conn  , "connectivity",   &GetAtomConnection),
    SHOW_ATOM_DATA_ENTRY(Conn  , "hbconnectivity", &GetAtomHydrogenBond),
    SHOW_ATOM_COORDINATES_ENTRY,
    SHOW_ATOM_COLOUR_ENTRY,
    SHOW_ATOM_FORCE_ENTRY,
    SHOW_ATOM_VELOCITY_ENTRY,
    /* general atom data */
    SHOW_ATOM_CROSS_ENTRY,
    SHOW_ATOM_IDENTIFY_ENTRY,
    SHOW_ATOM_SELECTION_ENTRY,
    SHOW_ATOM_WINDOW_ENTRY_ENTRY,
    SHOW_ATOM_RECONNECTIVITY_ENTRY,
    SHOW_ATOM_HBRECONNECTIVITY_ENTRY,
    SHOW_ATOM_MAXCONNECTIVITY_ENTRY,
    /* element data */
    SHOW_ATOM_INFO_ENTRY,
    /* name stack */
    SHOW_ATOM_STACK_ENTRY,
    /* structure */
    SHOW_ATOM_STRUCT_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

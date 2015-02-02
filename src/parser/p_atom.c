/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include <limits.h>
#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>

#include "gomtcl.h"
#include "label.h"
#include "molecule.h"
#include "molecstruct.h"
#include "plot.h"
#include "selection.h"
#include "parser.h"

#include "stdafx.h"

int gomp_HasMolecStructs(void)
{
    if ( gomp_GetNumMolecStructs() )
        return(1);
    gomp_PrintERROR("no molecule structure defined");
    return(0);
}

#define ALIAS GOM_PARSER_CMD_ALIAS

/* atom colour bycharge Seg Res Atm Fmin Fmax */
/* atom color  bycharge Seg Res Atm Fmin Fmax */
/* atom colour charge   Seg Res Atm Fmin Fmax */
/* atom color  charge   Seg Res Atm Fmin Fmax */
#define PARSE_ATOM_COLOUR_BYCHARGE_ENTRY \
    GOM_PARSER_FINAL_CMD(\
        "bycharge",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm Fmin Fmax?",0,5),\
        ParseAtomColourByCharge,GOM_PARSER_UNUSED_VALUE),\
    GOM_PARSER_FINAL_CMD(\
        "charge",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm Fmin Fmax?",0,5),\
        ParseAtomColourByCharge,GOM_PARSER_UNUSED_VALUE)
static int ParseAtomColourByCharge(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    gom_SelectionList atoms;
    double Fmin,Fmax,*pFmin,*pFmax;

    GOM_PARSER_START_ARGS
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_ARG(Fmin) GOM_PARSER_ARG(Fmax)
        GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom selection list",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,atoms) );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "minimum atomic partial charge",
        GOM_PARSER_NULL_POINTER(Fmin) ||
        GOM_PARSER_RETRIEVE_DOUBLE(Fmin) );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "maximum atomic partial charge",
        GOM_PARSER_NULL_POINTER(Fmax) ||
        GOM_PARSER_RETRIEVE_DOUBLE(Fmax) );
    GOM_PARSER_VERIFY(
        gomp_ParseColourListByCharge(&atoms,pFmin,pFmax) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom colour byfourth Seg Res Atm */
/* atom color  fourth   Seg Res Atm */
#define PARSE_ATOM_COLOUR_BYFOURTH_ENTRY \
    GOM_PARSER_FINAL_CMD(\
        "byfourth",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm?",0,3),\
        ParseAtomColourByFourth,GOM_PARSER_UNUSED_VALUE),\
    GOM_PARSER_FINAL_CMD(\
        "fourth",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm?",0,3),\
        ParseAtomColourByFourth,GOM_PARSER_UNUSED_VALUE)
static int ParseAtomColourByFourth(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    gom_SelectionList atoms;

    GOM_PARSER_START_ARGS
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom selection list",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,atoms) );
    GOM_PARSER_VERIFY(
        gomp_ParseColourListFourth(&atoms) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom colour byresnumber Seg Res Atm */
/* atom color  resnumber   Seg Res Atm */
#define PARSE_ATOM_COLOUR_BYRESNUMBER_ENTRY \
    GOM_PARSER_FINAL_CMD(\
        "byresnumber",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm?",0,3),\
        ParseAtomColourByResidue,GOM_PARSER_UNUSED_VALUE),\
    GOM_PARSER_FINAL_CMD(\
        "resnumber",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm?",0,3),\
        ParseAtomColourByResidue,GOM_PARSER_UNUSED_VALUE)
static int ParseAtomColourByResidue(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    gom_SelectionList atoms;

    GOM_PARSER_START_ARGS
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom selection list",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,atoms) );
    GOM_PARSER_VERIFY(
        gomp_ParseColourListResidueNumber(&atoms) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom colour byvector Seg Res Atm Fmin Fmax */
/* atom color  byvector Seg Res Atm Fmin Fmax */
/* atom colour vector   Seg Res Atm Fmin Fmax */
/* atom color  vector   Seg Res Atm Fmin Fmax */
#define PARSE_ATOM_COLOUR_BYVECTOR_ENTRY \
    GOM_PARSER_FINAL_CMD(\
        "byvector",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm Fmin Fmax?",0,5),\
        ParseAtomColourByVector,GOM_PARSER_UNUSED_VALUE),\
    GOM_PARSER_FINAL_CMD(\
        "vector",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm Fmin Fmax?",0,5),\
        ParseAtomColourByVector,GOM_PARSER_UNUSED_VALUE)
static int ParseAtomColourByVector(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    gom_SelectionList atoms;
    double Fmin,Fmax,*pFmin,*pFmax;

    GOM_PARSER_START_ARGS
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_ARG(Fmin) GOM_PARSER_ARG(Fmax)
        GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom selection list",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,atoms) );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "minimum vector array value",
        GOM_PARSER_NULL_POINTER(Fmin) ||
        GOM_PARSER_RETRIEVE_DOUBLE(Fmin) );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "maximum vector array value",
        GOM_PARSER_NULL_POINTER(Fmax) ||
        GOM_PARSER_RETRIEVE_DOUBLE(Fmax) );
    GOM_PARSER_VERIFY(
        gomp_ParseColourListByVector(&atoms,pFmin,pFmax) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom colour Seg Res Atm Colour */
/* atom color  Seg Res Atm Colour */
#define PARSE_ATOM_COLOUR_ARGS_ENTRY \
    GOM_PARSER_FINAL_ARGS(\
        GOM_PARSER_NEED_ARGS("Seg Res Atm Colour",4),\
        ParseAtomColour,GOM_PARSER_UNUSED_VALUE)
static int ParseAtomColour(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    gom_SelectionList atoms;
    gom_FloatColour   Colour;

    GOM_PARSER_START_ARGS
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_ARG(Colour)
        GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom selection list",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,atoms) );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "colour",
        GOM_PARSER_RETRIEVE_COLOUR(Colour) );
    GOM_PARSER_VERIFY(
        gomp_ParseColourList(&atoms,&Colour) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom colour ... */
/* atom color  ... */
#define PARSE_ATOM_COLOUR_ENTRY \
    GOM_PARSER_CMD_PART(\
        "colour" ALIAS "color",GOM_PARSER_NO_MORE_ARGS,\
        parseAtomColour,GOM_PARSER_UNUSED_VALUE)
const static gom_ParserArgumentList parseAtomColour[] = {
    PARSE_ATOM_COLOUR_BYCHARGE_ENTRY,
    PARSE_ATOM_COLOUR_BYVECTOR_ENTRY,
    PARSE_ATOM_COLOUR_BYFOURTH_ENTRY,
    PARSE_ATOM_COLOUR_BYRESNUMBER_ENTRY,
    PARSE_ATOM_COLOUR_ARGS_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* atom cpk  Seg Res Atm */
/* atom -cpk Seg Res Atm */
#define PARSE_ATOM_CPK_ENTRY \
    GOM_PARSER_FINAL_CMD(\
        "cpk",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm?",0,3),\
        ParseAtomCpk,GOM_PARSER_SET_VALUE(ON)),\
    GOM_PARSER_FINAL_CMD(\
        "-cpk",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm?",0,3),\
        ParseAtomCpk,GOM_PARSER_SET_VALUE(OFF))
static int ParseAtomCpk(GOM_PARSER_ARGLIST,intptr_t Show)
{
    gom_SelectionList atoms;

    GOM_PARSER_START_ARGS
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom selection list",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,atoms) );
    GOM_PARSER_VERIFY(
        gomp_ParseTypeCPKList(Show,&atoms) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom display  Seg1 Res1 Atm1 around FRadius Seg2 Res2 Atm2 ?Colour? */
/* atom display  Seg1 Res1 Atm1 around FRadius Fx   Fy   Fz   ?Colour? */
#define PARSE_ATOM_NOMINUS_DISPLAY_AROUND_ATOMS_ENTRY \
    GOM_PARSER_FINAL_ARGS(\
        GOM_PARSER_NEED_ARGS_RANGE("?Seg2 Res2 Atm2 Colour?",0,4),\
        ParseAtomDisplayAroundAtoms,GOM_PARSER_INHERIT_VALUE),\
    GOM_PARSER_FINAL_ARGS(\
        GOM_PARSER_NEED_ARGS_RANGE("Fx Fy Fz ?Colour?",3,4),\
        ParseAtomDisplayAroundAtoms,GOM_PARSER_INHERIT_VALUE)
/* atom -display Seg1 Res1 Atm1 around FRadius Seg2 Res2 Atm2          */
/* atom -display Seg1 Res1 Atm1 around FRadius Fx   Fy   Fz            */
#define PARSE_ATOM_MINUS_DISPLAY_AROUND_ATOMS_ENTRY \
    GOM_PARSER_FINAL_ARGS(\
        GOM_PARSER_NEED_ARGS_RANGE("?Seg2 Res2 Atm2?",0,3),\
        ParseAtomDisplayAroundAtoms,GOM_PARSER_INHERIT_VALUE),\
    GOM_PARSER_FINAL_ARGS(\
        GOM_PARSER_NEED_ARGS("Fx Fy Fz",3),\
        ParseAtomDisplayAroundAtoms,GOM_PARSER_INHERIT_VALUE)
static int ParseAtomDisplayAroundAtoms(GOM_PARSER_ARGLIST,intptr_t Display)
{
    gom_SelectionList      centre;
    double                 FRadius;
    gom_SelectionList      atoms;
    gom_FloatColour        Colour;
    const gom_FloatColour *pColour = &Colour;

    GOM_PARSER_START_ARGS
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_ARG(FRadius)
        GOM_PARSER_SELECTIONLIST_ARG(2)
        GOM_PARSER_ARG(Colour)
        GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "centre atoms",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,centre) );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "search radius",
        GOM_PARSER_RETRIEVE_DOUBLE(FRadius) );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom selection list",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(2,atoms) );
    if ( GOM_PARSER_HAS_ARG(Colour) ) {
        GOM_PARSER_VERIFY_GOM_PARSER(
            "colour",
            GOM_PARSER_RETRIEVE_COLOUR(Colour) );
    }
    else
        pColour = NULL;
    GOM_PARSER_VERIFY(
        gomp_ControlSelectRoundAtoms(
            Display,&centre,FRadius,&atoms,pColour) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom display  Seg1 Res1 Atm1 around FRadius */
/* atom -display Seg1 Res1 Atm1 around FRadius */
#define PARSE_ATOM_DISPLAY_AROUND_ARGS_ENTRY \
    GOM_PARSER_FINAL_ARGS(\
        GOM_PARSER_NO_MORE_ARGS,\
        ParseAtomDisplayAround,GOM_PARSER_INHERIT_VALUE)
static int ParseAtomDisplayAround(GOM_PARSER_ARGLIST,intptr_t Display)
{
    gom_SelectionList centre;
    double FRadius;
    const static gom_SelectionList all_atoms = {"*","*","*"};

    GOM_PARSER_START_ARGS
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_ARG(FRadius)
        GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "centre atoms",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,centre) );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "search radius",
        GOM_PARSER_RETRIEVE_DOUBLE(FRadius) );
    GOM_PARSER_VERIFY(
        gomp_ControlSelectRoundAtoms(
            Display,&centre,FRadius,&all_atoms,NULL) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom display  Seg1 Res1 Atm1 around FRadius ... */
#define PARSE_ATOM_NOMINUS_DISPLAY_AROUND_ENTRY \
    GOM_PARSER_CMD_PART(\
        "around",GOM_PARSER_NEED_ARGS("FRadius",1),\
        parseAtomNoMinusDisplayAround,GOM_PARSER_INHERIT_VALUE)
/* atom -display Seg1 Res1 Atm1 around FRadius ... */
#define PARSE_ATOM_MINUS_DISPLAY_AROUND_ENTRY \
    GOM_PARSER_CMD_PART(\
        "around",GOM_PARSER_NEED_ARGS("FRadius",1),\
        parseAtomMinusDisplayAround,GOM_PARSER_INHERIT_VALUE)
const static gom_ParserArgumentList parseAtomNoMinusDisplayAround[] = {
    PARSE_ATOM_NOMINUS_DISPLAY_AROUND_ATOMS_ENTRY,
    PARSE_ATOM_DISPLAY_AROUND_ARGS_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};
const static gom_ParserArgumentList parseAtomMinusDisplayAround[] = {
    PARSE_ATOM_MINUS_DISPLAY_AROUND_ATOMS_ENTRY,
    PARSE_ATOM_DISPLAY_AROUND_ARGS_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* atom display  Seg Res Atm */
/* atom -display Seg Res Atm */
#define PARSE_ATOM_DISPLAY_ARGS_ENTRY \
    GOM_PARSER_FINAL_ARGS(\
        GOM_PARSER_NO_MORE_ARGS,\
        ParseAtomDisplay,GOM_PARSER_INHERIT_VALUE)
static int ParseAtomDisplay(GOM_PARSER_ARGLIST,intptr_t Display)
{
    gom_SelectionList atoms;

    GOM_PARSER_START_ARGS
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom selection list",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,atoms) );
    GOM_PARSER_VERIFY(
        gomp_ParseDisplayList(Display,&atoms) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom display  Seg1 Res1 Atm1 ... */
/* atom -display Seg1 Res1 Atm1 ... */
#define PARSE_ATOM_DISPLAY_ENTRY \
    GOM_PARSER_CMD_PART(\
        "display",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm?",0,3),\
        parseAtomNoMinusDisplay,GOM_PARSER_SET_VALUE(ON)),\
    GOM_PARSER_CMD_PART(\
        "-display",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm?",0,3),\
        parseAtomMinusDisplay,GOM_PARSER_SET_VALUE(OFF))
const static gom_ParserArgumentList parseAtomNoMinusDisplay[] = {
    PARSE_ATOM_NOMINUS_DISPLAY_AROUND_ENTRY,
    PARSE_ATOM_DISPLAY_ARGS_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};
const static gom_ParserArgumentList parseAtomMinusDisplay[] = {
    PARSE_ATOM_MINUS_DISPLAY_AROUND_ENTRY,
    PARSE_ATOM_DISPLAY_ARGS_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* atom label full    */
/* atom label atom    */
/* atom label residue */
#define PARSE_ATOM_NOMINUS_LABEL_TYPE_ENTRY \
    GOM_PARSER_FINAL_CMD(\
        "full",GOM_PARSER_NO_MORE_ARGS,\
        ParseAtomLabelType,GOM_PARSER_SET_VALUE(FULL_LABEL_TYPE)),\
    GOM_PARSER_FINAL_CMD(\
        "atom",GOM_PARSER_NO_MORE_ARGS,\
        ParseAtomLabelType,GOM_PARSER_SET_VALUE(ATOM_LABEL_TYPE)),\
    GOM_PARSER_FINAL_CMD(\
        "residue",GOM_PARSER_NO_MORE_ARGS,\
        ParseAtomLabelType,GOM_PARSER_SET_VALUE(RESIDUE_LABEL_TYPE))
static int ParseAtomLabelType(GOM_PARSER_ARGLIST,intptr_t Type)
{
    GOM_PARSER_VERIFY(gomp_SetPlotLabelType(Type) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom label  Seg Res Atm */
#define PARSE_ATOM_NOMINUS_LABEL_ARGS_ENTRY \
    GOM_PARSER_FINAL_ARGS(\
        GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm?",0,3),\
        ParseAtomLabel,GOM_PARSER_INHERIT_VALUE)
/* atom label  ...         */
/* atom -label Seg Res Atm */
#define PARSE_ATOM_LABEL_ENTRY \
    GOM_PARSER_CMD_PART(\
        "label",GOM_PARSER_NO_MORE_ARGS,\
        parseAtomNoMinusLabel,GOM_PARSER_SET_VALUE(ON)),\
    GOM_PARSER_FINAL_CMD(\
        "-label",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm?",0,3),\
        ParseAtomLabel,GOM_PARSER_SET_VALUE(OFF))
static int ParseAtomLabel(GOM_PARSER_ARGLIST,intptr_t Show)
{
    gom_SelectionList atoms;

    GOM_PARSER_START_ARGS
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom selection list",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,atoms) );
    GOM_PARSER_VERIFY(gomp_ParseAtomLabelList(Show,&atoms) == 0);
    GOM_PARSER_SUCCEEDED;
}
const static gom_ParserArgumentList parseAtomNoMinusLabel[] = {
    PARSE_ATOM_NOMINUS_LABEL_TYPE_ENTRY,
    PARSE_ATOM_NOMINUS_LABEL_ARGS_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* atom licorice  Seg Res Atm */
/* atom -licorice Seg Res Atm */
#define PARSE_ATOM_LICORICE_ENTRY \
    GOM_PARSER_FINAL_CMD(\
        "licorice",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm?",0,3),\
        ParseAtomLicorice,GOM_PARSER_SET_VALUE(ON)),\
    GOM_PARSER_FINAL_CMD(\
        "-licorice",GOM_PARSER_NEED_ARGS_RANGE("?Seg Res Atm?",0,3),\
        ParseAtomLicorice,GOM_PARSER_SET_VALUE(OFF))
static int ParseAtomLicorice(GOM_PARSER_ARGLIST,intptr_t Show)
{
    gom_SelectionList atoms;

    GOM_PARSER_START_ARGS
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom selection list",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,atoms) );
    GOM_PARSER_VERIFY(gomp_ParseTypeLicoList(Show,&atoms) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom scale cpk Fscale Seg Res Atm */
#define PARSE_ATOM_SCALE_CPK_ENTRY \
    GOM_PARSER_FINAL_CMD(\
        "cpk",GOM_PARSER_NEED_ARGS_RANGE("Fscale ?Seg Res Atm?",1,4),\
        ParseAtomScaleCpk,GOM_PARSER_UNUSED_VALUE)
static int ParseAtomScaleCpk(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    double Fscale;
    gom_SelectionList atoms;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Fscale)
        GOM_PARSER_SELECTIONLIST_ARG(1)
        GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "CPK scale",
        GOM_PARSER_RETRIEVE_DOUBLE(Fscale) );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom selection list",
        GOM_PARSER_RETRIEVE_SELECTIONLIST(1,atoms) );
    GOM_PARSER_VERIFY(gomp_ParseCPKScaleList(Fscale,&atoms) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom scale ... */
#define PARSE_ATOM_SCALE_ENTRY \
    GOM_PARSER_CMD_PART(\
        "scale",GOM_PARSER_NO_MORE_ARGS,\
        parseAtomScale,GOM_PARSER_UNUSED_VALUE)
const static gom_ParserArgumentList parseAtomScale[] = {
    PARSE_ATOM_SCALE_CPK_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* atom selection ... */
#define PARSE_ATOM_SELECTION_ENTRY \
    GOM_PARSER_CMD_PART(\
        "selection",GOM_PARSER_NO_MORE_ARGS,\
        parseAtomSelection,GOM_PARSER_UNUSED_VALUE)
static int ParseAtomSelection(GOM_PARSER_ARGLIST,intptr_t Mode)
{
    GOM_PARSER_VERIFY(gomp_SetAtomSelectionMode(Mode) == 0);
    GOM_PARSER_SUCCEEDED;
}
const static gom_ParserArgumentList parseAtomSelection[] = {
    GOM_PARSER_FINAL_CMD("atom",GOM_PARSER_NO_MORE_ARGS,
                         ParseAtomSelection,
                         GOM_PARSER_SET_VALUE(ATOM_SELECTION)),
    GOM_PARSER_FINAL_CMD("residue",GOM_PARSER_NO_MORE_ARGS,
                         ParseAtomSelection,
                         GOM_PARSER_SET_VALUE(RESIDUE_SELECTION)),
    GOM_PARSER_END_ARGUMENT_LIST
};

/* atom structure append Name Natoms */
/* atom structure new    Name Natoms */
#define PARSE_ATOM_STRUCT_CREATE_ENTRY \
    GOM_PARSER_FINAL_CMD(\
        "append",GOM_PARSER_NEED_ARGS("Name Natoms",2),\
        ParseAtomStructCreate,GOM_PARSER_SET_VALUE(APPEND)),\
    GOM_PARSER_FINAL_CMD(\
        "new",GOM_PARSER_NEED_ARGS("Name Natoms",2),\
        ParseAtomStructCreate,GOM_PARSER_SET_VALUE(NEW))
static int ParseAtomStructCreate(GOM_PARSER_ARGLIST,intptr_t Type)
{
    char const *Name;
    int         Natoms;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Name)
        GOM_PARSER_ARG(Natoms)
    GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "structure name",
        GOM_PARSER_RETRIEVE_STRING(Name));
    GOM_PARSER_VERIFY_GOM_PARSER(
        "atom count",
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(Natoms,1,INT_MAX));
    GOM_PARSER_VERIFY(gomp_CreateMolecStruct(Name,Natoms,Type) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom structure delete Istruct */
#define PARSE_ATOM_STRUCT_DELETE_ENTRY \
    GOM_PARSER_FINAL_CMD(\
        "delete",GOM_PARSER_NEED_ARGS_RANGE("?Istruct?",0,1),\
        ParseAtomStructDelete,GOM_PARSER_UNUSED_VALUE)
static int ParseAtomStructDelete(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    int Istruct;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Istruct)
    GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        "structure index",
        GOM_PARSER_DEFAULT(Istruct,-1) ||
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Istruct,1,gomp_GetNumMolecStructs()));

    GOM_PARSER_VERIFY(
        ( Istruct > 0 ?
          gomp_DeleteMolecStruct(Istruct-1) :
          gomp_DeleteMolecStructs() ) == 0 );
    GOM_PARSER_SUCCEEDED;
}

/* atom structure merge Istruct Dstruct */
#define PARSE_ATOM_STRUCT_MERGE_ENTRY \
    GOM_PARSER_FINAL_CMD(\
        "merge",GOM_PARSER_NEED_ARGS_RANGE("?Istruct Dstruct?",0,2),\
        ParseAtomStructMerge,GOM_PARSER_UNUSED_VALUE)
static int ParseAtomStructMerge(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    int Istruct, Dstruct;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Istruct)
        GOM_PARSER_ARG(Dstruct)
    GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY( gomp_HasMolecStructs() );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "structure index",
        GOM_PARSER_DEFAULT(Istruct,-1) ||
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Istruct,1,gomp_GetNumMolecStructs()));
    GOM_PARSER_VERIFY_GOM_PARSER(
        "destination structure index",
        GOM_PARSER_DEFAULT(Dstruct,1) ||
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Dstruct,1,gomp_GetNumMolecStructs()));

    GOM_PARSER_VERIFY(
        ( Istruct > 0 ?
          gomp_MergeMolecStruct(Istruct-1,Dstruct-1) :
          gomp_MergeMolecStructs() ) == 0 );
    GOM_PARSER_SUCCEEDED;
}

/* atom structure rename Name Istruct */
#define PARSE_ATOM_STRUCT_RENAME_ENTRY \
    GOM_PARSER_FINAL_CMD(\
        "rename",GOM_PARSER_NEED_ARGS_RANGE("Name ?Istruct?",1,2),\
        ParseAtomStructRename,GOM_PARSER_UNUSED_VALUE)
static int ParseAtomStructRename(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    char const *Name;
    int         Istruct;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Name)
        GOM_PARSER_ARG(Istruct)
    GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY( gomp_HasMolecStructs() );
    GOM_PARSER_VERIFY_GOM_PARSER(
        "structure name",
        GOM_PARSER_RETRIEVE_STRING(Name));
    GOM_PARSER_VERIFY_GOM_PARSER(
        "structure index",
        GOM_PARSER_DEFAULT(Istruct,1) ||
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Istruct,1,gomp_GetNumMolecStructs()));
    GOM_PARSER_VERIFY(gomp_PutMolecStructFileName(Istruct-1,Name) == 0);
    GOM_PARSER_SUCCEEDED;
}

/* atom structure ... */
#define PARSE_ATOM_STRUCT_ENTRY \
    GOM_PARSER_CMD_PART(\
        "structure" ALIAS "structures",GOM_PARSER_NO_MORE_ARGS,\
        parseAtomStruct,GOM_PARSER_UNUSED_VALUE)
static const gom_ParserArgumentList parseAtomStruct[] = {
    PARSE_ATOM_STRUCT_CREATE_ENTRY,
    PARSE_ATOM_STRUCT_DELETE_ENTRY,
    PARSE_ATOM_STRUCT_MERGE_ENTRY,
    PARSE_ATOM_STRUCT_RENAME_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

const gom_ParserArgumentList gomp_AtomCommand[] = {
    PARSE_ATOM_CPK_ENTRY,
    PARSE_ATOM_COLOUR_ENTRY,
    PARSE_ATOM_DISPLAY_ENTRY,
    PARSE_ATOM_LABEL_ENTRY,
    PARSE_ATOM_LICORICE_ENTRY,
    PARSE_ATOM_SCALE_ENTRY,
    PARSE_ATOM_SELECTION_ENTRY,
    PARSE_ATOM_STRUCT_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

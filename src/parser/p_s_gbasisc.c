/*

Copyright (c) 2003 - 2004 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#include "gaussian.h"
#include "parser.h"

#include "stdafx.h"

/* show gbasis set */
#define SHOW_GBASIS_SET_ENTRY \
GOM_PARSER_FINAL_INT_CMD("set",&parseShowGbasisSet)
static const IntFunc parseShowGbasisSet = { gomp_GetNumberOfGaussianBasisSets};

/* show gbasis tag Itag */
#define SHOW_GBASIS_TAG_ENTRY \
GOM_PARSER_FINAL_CMD("tag",\
                     GOM_PARSER_NEED_ARGS("Itag",1),\
                     ParseShowGbasisTag,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowGbasisTag(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    int Itag;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Itag)
    GOM_PARSER_END_ARGS;
    
    GOM_PARSER_VERIFY_GOM_PARSER(
        "tag index",
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Itag, 1, gomp_GetNumberOfGaussianBasisSets()) );

    GOM_PARSER_RETURN_STRING( gomp_GetGaussianBasisTag(Itag - 1) );

    GOM_PARSER_SUCCEEDED;
}

/* show gbasis entry Ientry */
#define SHOW_GBASIS_ENTRY_ENTRY \
GOM_PARSER_FINAL_CMD("entry",\
                     GOM_PARSER_NEED_ARGS("Ientry",1),\
                     ParseShowGbasisEntry,\
                     GOM_PARSER_UNUSED_VALUE)
static int ParseShowGbasisEntry(GOM_PARSER_ARGLIST,intptr_t Flags)
{
    int Ientry;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Ientry)
    GOM_PARSER_END_ARGS;
    
    GOM_PARSER_VERIFY_GOM_PARSER(
        "entry index",
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Ientry, 1, gomp_GetNumberOfGaussianBasisSets()) );

    GOM_PARSER_RETURN_STRING( gomp_GetGaussianBasisSetEntry(Ientry - 1) );

    GOM_PARSER_SUCCEEDED;
}

const gom_ParserArgumentList gomp_ShowGbasisCommand[] = {
    SHOW_GBASIS_SET_ENTRY,
    SHOW_GBASIS_TAG_ENTRY,
    SHOW_GBASIS_ENTRY_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

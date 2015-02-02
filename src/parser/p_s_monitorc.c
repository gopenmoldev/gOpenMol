/*

Copyright (c) 2003 - 2004 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#include "gommonitor.h"
#include "parser.h"

#include "stdafx.h"

#define ALIAS GOM_PARSER_CMD_ALIAS

typedef struct {
    int Knots;
    const char *argmsg;
    int (*Samples)(void);
    int (*State)(void);
    int (*Colour)(int,float*,float*,float*);
    int (*Type)(int);
    const int* (*List)(void);
    int (*Series)(void);
} MonitorFunc;

/* show monitor <type> state */
#define SHOW_MONITOR_GENERAL_STATE_ENTRY \
GOM_PARSER_FINAL_CMD("state",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowMonitorGeneralState,\
                     GOM_PARSER_INHERIT_VALUE)
static int ParseShowMonitorGeneralState(GOM_PARSER_ARGLIST,intptr_t Ptr)
{
    MonitorFunc *mf = GOM_PARSER_GET_POINTER_VALUE(Ptr);
    GOM_PARSER_RETURN_BOOLEAN(mf->State());
    GOM_PARSER_SUCCEEDED;
}

/* show monitor <type> colour Index */
/* show monitor <type> color  Index */
#define SHOW_MONITOR_GENERAL_COLOUR_ENTRY \
GOM_PARSER_FINAL_CMD("colour" ALIAS "color",\
                     GOM_PARSER_NEED_ARGS("Index",1),\
                     ParseShowMonitorGeneralColour,\
                     GOM_PARSER_INHERIT_VALUE)
static int ParseShowMonitorGeneralColour(GOM_PARSER_ARGLIST,intptr_t Ptr)
{
    MonitorFunc *mf = GOM_PARSER_GET_POINTER_VALUE(Ptr);
    int Index;
    float red, green, blue;
    
    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Index)
    GOM_PARSER_END_ARGS;

    GOM_PARSER_VERIFY_GOM_PARSER(
        mf->argmsg,
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Index, 1, mf->Samples() / mf->Knots));

    (void)mf->Colour(Index - 1, &red, &green, &blue);
    GOM_PARSER_RETURN_LIST(("%f %f %f",red,green,blue));
    GOM_PARSER_SUCCEEDED;
}

/* show monitor <type> type Index */
#define SHOW_MONITOR_GENERAL_TYPE_ENTRY \
GOM_PARSER_FINAL_CMD("type",\
                     GOM_PARSER_NEED_ARGS("Index",1),\
                     ParseShowMonitorGeneralType,\
                     GOM_PARSER_INHERIT_VALUE)
static int ParseShowMonitorGeneralType(GOM_PARSER_ARGLIST,intptr_t Ptr)
{
    MonitorFunc *mf = GOM_PARSER_GET_POINTER_VALUE(Ptr);
    int Index;
    
    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Index)
    GOM_PARSER_END_ARGS;
    
    GOM_PARSER_VERIFY_GOM_PARSER(
        mf->argmsg,
        GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(
            Index, 1, mf->Samples() / mf->Knots));

    GOM_PARSER_RETURN_INT(mf->Type(Index - 1));
    GOM_PARSER_SUCCEEDED;
}

/* show monitor <type> list  */
/* show monitor <type> arroy */
#define SHOW_MONITOR_GENERAL_LIST_ENTRY \
GOM_PARSER_FINAL_CMD("list",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowMonitorGeneralList,\
                     GOM_PARSER_INHERIT_VALUE),\
GOM_PARSER_FINAL_CMD("arroy",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowMonitorGeneralList,\
                     GOM_PARSER_INHERIT_VALUE)
static int ParseShowMonitorGeneralList(GOM_PARSER_ARGLIST,intptr_t Ptr)
{
    MonitorFunc *mf   = GOM_PARSER_GET_POINTER_VALUE(Ptr);
    gom_ParserList list;
    const int   *slist;
    int i,j;

    if ( ! GOM_PARSER_LIST_INIT(list) )
        goto list_failed;
    /* number of items */
    if ( ! GOM_PARSER_LIST_APPEND_INT(list, mf->Samples() / mf->Knots) )
        goto list_failed;
    /* atom indeces of the entry */
    slist = mf->List();
    for ( i = 0 ; i < mf->Samples() ; i += mf->Knots ) {
        gom_ParserList sublist;
        if ( ! GOM_PARSER_LIST_INIT(sublist) )
            goto sublist_failed;
        for ( j = 0 ; j < mf->Knots ; j++ ) {
            if ( ! GOM_PARSER_LIST_APPEND_INT(sublist, slist[i+j]+1) )
                goto sublist_failed;
        }
        if ( ! GOM_PARSER_LIST_APPEND_LIST(list,sublist) )
            goto sublist_failed;
        continue;
  sublist_failed:
        GOM_PARSER_LIST_FREE(sublist);
        goto list_failed;
    }

    GOM_PARSER_LIST_RETURN(list);
    GOM_PARSER_SUCCEEDED;
  list_failed:
    GOM_PARSER_LIST_FREE(list);
    GOM_PARSER_FAILED;
}

/* show monitor <type> timeseries */
#define SHOW_MONITOR_GENERAL_TIME_SERIES_ENTRY \
GOM_PARSER_FINAL_CMD("timeseries",\
                     GOM_PARSER_NO_MORE_ARGS,\
                     ParseShowMonitorGeneralTimeSeries,\
                     GOM_PARSER_INHERIT_VALUE)
static int ParseShowMonitorGeneralTimeSeries(GOM_PARSER_ARGLIST, intptr_t Ptr)
{
    MonitorFunc *mf = GOM_PARSER_GET_POINTER_VALUE(Ptr);
    GOM_PARSER_RETURN_INT(mf->Series());
    GOM_PARSER_SUCCEEDED;
}


/* show monitor <type> */
static const gom_ParserArgumentList parseShowMonitorGeneral[] = {
    SHOW_MONITOR_GENERAL_STATE_ENTRY,
    SHOW_MONITOR_GENERAL_COLOUR_ENTRY,
    SHOW_MONITOR_GENERAL_TYPE_ENTRY,
    SHOW_MONITOR_GENERAL_LIST_ENTRY,
    SHOW_MONITOR_GENERAL_TIME_SERIES_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

/* show monitor distance ... */
#define SHOW_MONITOR_DISTANCE_ENTRY \
GOM_PARSER_CMD_PART("distance",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowMonitorGeneral,\
                    GOM_PARSER_SET_POINTER_VALUE(&DistFunc))
static const MonitorFunc DistFunc = {
    2,
    "distance monitor index",
    gomp_GetDistMonitorSamples,
    gomp_GetDistMonitor,
    gomp_GetDistMonitorColor,
    gomp_GetDistMonitorType,
    gomp_GetDistMonitorSamplesList,
    gomp_GetNumberDistanceSeries
};

/* show monitor angle ... */
#define SHOW_MONITOR_ANGLE_ENTRY \
GOM_PARSER_CMD_PART("angle",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowMonitorGeneral,\
                    GOM_PARSER_SET_POINTER_VALUE(&AngleFunc))
static const MonitorFunc AngleFunc = {
    3,
    "angle monitor index",
    gomp_GetAngMonitorSamples,
    gomp_GetAngMonitor,
    gomp_GetAngMonitorColor,
    gomp_GetAngMonitorType,
    gomp_GetAngMonitorSamplesList,
    gomp_GetNumberAngleSeries
};

/* show monitor torsion ... */
#define SHOW_MONITOR_TORSION_ENTRY \
GOM_PARSER_CMD_PART("torsion",\
                    GOM_PARSER_NO_MORE_ARGS,\
                    parseShowMonitorGeneral,\
                    GOM_PARSER_SET_POINTER_VALUE(&TorsFunc))
static const MonitorFunc TorsFunc = {
    4,
    "torsion monitor index",
    gomp_GetTorsMonitorSamples,
    gomp_GetTorsMonitor,
    gomp_GetTorsMonitorColor,
    gomp_GetTorsMonitorType,
    gomp_GetTorsMonitorSamplesList,
    gomp_GetNumberTorsionSeries
};

const gom_ParserArgumentList gomp_ShowMonitorCommand[] = {
    SHOW_MONITOR_DISTANCE_ENTRY,
    SHOW_MONITOR_ANGLE_ENTRY,
    SHOW_MONITOR_TORSION_ENTRY,
    GOM_PARSER_END_ARGUMENT_LIST
};

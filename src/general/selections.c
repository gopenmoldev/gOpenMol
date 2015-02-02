/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2004 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#if !defined(WIN32)
#include <sys/types.h>
#include <sys/time.h>
/*#include <sys/resource.h>*/
#endif

#include <tcl.h>

#if defined(IRIX)
#include <sys/procfs.h>
#endif

#include "axis.h"
#include "bond.h"
#include "cell.h"
#include "colouring.h"
#include "coord_man.h"
#include "g_status.h"
#include "gomstring.h"
#include "label.h"
#include "listutils.h"
#include "measure.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "plot.h"
#include "plumber.h"
#include "printmsg.h"
#include "projview.h"
#include "rdf.h"
#include "rforce.h"
#include "selection.h"
#include "tclutils.h"

#include "stdafx.h"

#define MAX_FLOAT  1.e+35    /* This max float is just a guess */
#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define Rabs(a)    ( ( a ) > 0   ? (a) : -(a))
#define Fabs(a)    ( ( a ) > 0.0 ? (a) : -(a))
#define SCALE(a,b,c)     scaleO(a,b,c) /* use own scale function */

#define   SEGMENT      MAX_SEG_NAME_LEN
#define   RESIDUE      MAX_RES_NAME_LEN
#define   ATOM         MAX_ATM_NAME_LEN
#define   SEP_STRUCT   '^'

/* definitions                   */
#define  NEW_SELECTION_LIST    1 /* create new selection list     */
#define  OLD_SELECTION_LIST    0 /* use old selection list        */
#define  DEL_SELECTION_LIST   -1 /* free the selection list space */

/* functions */

static int gomAtomSelectionMode = 0;
static int gomSelectionModeStatus = 0;
static int ShowAtomsAround(int, int, int, const int *, int, float,
                         const gom_SelectionList *, const gom_FloatColour *);
static int ShowAtomsAround_xyz(int, float, float, float, float,
                             const gom_SelectionList *,
                             const gom_FloatColour *);

static int collect_selected_atoms;  /* help variable to collect atoms */

/* functions and variables to support a save of the atoms that
   are selected during a command */

static int *AtomHitList           = (int *)NULL;
static int  AtomHitListTextLength = 0;
static int  AtomHitListActive     = 0;

static int  InitAtomHitList(void);
static int  SaveAtomHitList(int , const int *);

/* structures */

static struct {
    int  Length;
    int  State; /* if != 0 all are selected, if == 0 at least one structure is not selected */ 
    int *List;
} Selected = { 0  , 0 , (int *)NULL};

static struct {
    int   Structures;
    float Xaxis;
    float Yaxis;
    float Zaxis;
    int *Structure;
    int *Length;
    int *List; 
    int   CoordPoints;
    float *Xc;
    float *Yc;
    float *Zc;
    gom_PlotterData Plotter;
} CoordAxisList = {
    0, 0.0, 0.0, 0.0, 
    (int *)NULL, (int *)NULL, (int *)NULL,
    0, (float *)NULL, (float *)NULL, (float *)NULL,
    { NULL, NULL }
};


#define CoordAxisDataIsChanging() \
    gomp_InvalidatePlotterDelayed(&CoordAxisList.Plotter)

/* end of declarations */

/*   Check the input names against the list and mark a '1' in the list if 
     it is found 

     This is the main engine for making selections from atom sets

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     This routine uses the Tcl string match routine:

     Tcl_StringMatch(String1,String2)

     that makes a glob-style matching

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*/

typedef struct {
    char *seg_list;
    char *res_list;
    int  (*res_ranges)[2];
    char *atm_list;
    int  (*atm_ranges)[2];
} SegResAtmListsAndRanges;

static void FreeSegResAtmListsAndRanges(SegResAtmListsAndRanges *lists)
{
    gomp_DataVectorFree(&lists->seg_list);
    gomp_DataVectorFree(&lists->res_list);
    gomp_DataVectorFree(&lists->res_ranges);
    gomp_DataVectorFree(&lists->atm_list);
    gomp_DataVectorFree(&lists->atm_ranges);
}

static int CheckStructureIndex(const char *input,int Which)
{
    int from,to,pos1,pos2;

    while ( *input == ',' || isspace((unsigned char)*input) )
        ++input;

    while ( *input != SEP_STRUCT ) {
        pos1 = pos2 = -1;
        /* Execution of a %n directive shouldn't increment      */
        /* the assignment count. Not all libraries follow       */
        /* the standard though, so don't test the return value. */ 
        sscanf(input,"%d%n-%d%n",&from,&pos1,&to,&pos2);
        if ( pos2 > 0 && (
                 input[pos2] == ',' || input[pos2] == SEP_STRUCT ||
                 isspace((unsigned char)input[pos2]) ) ) {
            /* We have a range. */
            input += pos2;
        }
        else if ( pos1 > 0 && (
                      input[pos1] == ',' || input[pos1] == SEP_STRUCT ||
                      isspace((unsigned char)input[pos1]) ) ) {
            /* We have an index. */
            input  += pos1;
            to      = from;
        }
        else {
            /* Not a number. Skip the item. */
            while ( *input != ',' && *input != SEP_STRUCT &&
                    !isspace((unsigned char)*input) )
                ++input;
        }
        if ( from - 1 <= Which && Which <= to - 1 )
            /* Structure index is in range. */
            return 0;
        /* Skip separator. */
        while ( *input == ',' || isspace((unsigned char)*input) )
            ++input;
    }

    /* Structure index is not in range. */
    return 1;
}

static int ParseListsAndRanges(const char *input,
                               char **list,int (**ranges)[2],
                               size_t list_item_size)
{
    static const DataVectorHandler vectorHandle = { NULL, NULL, NULL };
    int from,to,pos1,pos2,IsRange;

    while ( *input == ',' || isspace((unsigned char)*input) )
        ++input;

    if ( ! *input ) {
        /* Append asterix. */
        char *dst;
        /* Increase the array. */
        if ( ! gomp_DataVectorCreateOrAppend(
                 list,&vectorHandle,list_item_size) )
            return 1;
        /* Copy the name. */
        dst = *list + gomp_DataVectorGetSize(list) - list_item_size;
        strcpy(dst,"*");
        return 0;
    }

    IsRange = 0;

    while ( *input ) {
        if ( ranges ) {
            IsRange = 1;
            pos1 = pos2 = -1;
            /* Execution of a %n directive shouldn't increment      */
            /* the assignment count. Not all libraries follow       */
            /* the standard though, so don't test the return value. */ 
            sscanf(input,"%d%n-%d%n",&from,&pos1,&to,&pos2);
            if ( pos2 > 0 &&
                 ( input[pos2] == ',' || input[pos2] == '\0' ||
                   isspace((unsigned char)input[pos2]) ) ) {
                /* We have a range. */
                input += pos2;
                if ( from > to ) {
                    int temp = from;
                    from     = to;
                    to       = temp;
                }
            }
            else if ( pos1 > 0 && (
                          input[pos1] == ',' || input[pos1] == '\0' ||
                          isspace((unsigned char)input[pos1]) ) ) {
                /* We have an index. */
                input  += pos1;
                to      = from;
            }
            else
                IsRange = 0;
        }
        if ( IsRange ) {
            /* Increase the array. */
            if ( ! gomp_DataVectorCreateOrAppend(ranges,&vectorHandle,1) )
                return 1;
            /* Store. */
            (*ranges)[gomp_DataVectorGetSize(ranges)-1][0] = from;
            (*ranges)[gomp_DataVectorGetSize(ranges)-1][1] = to;
        }
        else {
            /* We have a name. */
            char *dst;
            size_t i;
            /* Increase the array. */
            if ( ! gomp_DataVectorCreateOrAppend(
                     list,&vectorHandle,list_item_size) )
                return 1;
            /* Copy the name. Truncate if the name is too long. */
            dst = *list + gomp_DataVectorGetSize(list) - list_item_size;
            for ( i = 0 ; i < list_item_size - 1 &&
                      *input != ',' && *input != '\0' &&
                      !isspace((unsigned char)*input) ; ++i )
                dst[i] = *input++;
            dst[i] = '\0';
            /* Omit the tail of the name. */
            while ( *input != ',' && *input != '\0' &&
                    !isspace((unsigned char)*input) )
                ++input;
        }
        /* Skip a separator. */
        while ( *input == ',' || isspace((unsigned char)*input) )
            ++input;
    }

    return 0;
}

/************************************************************************/
int gomp_MakeSelectionList(
    int        Which,          /* structure number    */
    const char *xsegment,      /* segment name (list) */
    const char *xresnam,       /* resideu name (list) */
    const char *xatnam,        /* atom name (list)    */
    int *sel_list)      /* index list to atoms */
/************************************************************************/
{
    static int i,j;
    static int sel_long;
    static int iserie,nseries,Count,from,to,atm_hit,res_hit,seg_hit;
/* different lists */
    static SegResAtmListsAndRanges lists;
    static const segment_name_t *SegPointer;
    static const residue_name_t *ResPointer;
    static const atom_name_t    *AtmPointer;
    static const int *Num1Pointer;
    static int Low,High;

    if( Which >= gomp_GetNumMolecStructs() ) 
        /* no structures ==> no atoms */
        return(0);

    sel_long = 0;

/*
  In the xsegment variable there can be as the first characters
  a number and the hat '^' where the number refers to the current
  structure.
*/
    if(xsegment[0] != '\0') {
        const char *hat = strchr(xsegment,SEP_STRUCT);
        if ( hat ) {
            if ( CheckStructureIndex(xsegment,Which) != 0 )
                return(0);
            xsegment = hat + 1;
        }
    }

    /* check for easy cases */
    if ( ( xsegment[0] == '\0' || strcmp(xsegment,"*") == 0 ) &&
         ( xresnam[0]  == '\0' || strcmp(xresnam,"*")  == 0 ) &&
         ( xatnam[0]   == '\0' || strcmp(xatnam,"*")   == 0 ) ) {
        High = gomp_GetNumAtomsInMolecStruct(Which);
        for(i = 0 ; i < High ; i++)
            sel_list[i] = i;
        return(High);
    };

/* go and hunt for the segment ... */
/*     tupper(wsegment);           */
/* check first for series */
    memset(&lists,0,sizeof(lists));
    if ( ParseListsAndRanges(xsegment,&lists.seg_list,NULL,
                             MAX_SEG_NAME_LEN+1) ||
         ParseListsAndRanges(xresnam,&lists.res_list,&lists.res_ranges,
                             MAX_RES_NAME_LEN+1) ||
         ParseListsAndRanges(xatnam,&lists.atm_list,&lists.atm_ranges,
                             MAX_ATM_NAME_LEN+1) ) {
        FreeSegResAtmListsAndRanges(&lists);
        return -1;
    }

/* do it now! */

    sel_long = 0;

    SegPointer  = gomp_GetAtomSegNamePointer(Which);
    ResPointer  = gomp_GetAtomResNamePointer(Which);
    AtmPointer  = gomp_GetAtomAtmNamePointer(Which);
    Num1Pointer = gomp_GetAtomResNum1Pointer(Which);

    Low  = 0;
    High = gomp_GetNumAtomsInMolecStruct(Which);

    if ( gomp_DataVectorGetSize(&lists.atm_list) == 0 ) {

        nseries = gomp_DataVectorGetSize(&lists.atm_ranges);
        for ( iserie = 0 ; iserie < nseries ; iserie++ ) {

            from = lists.atm_ranges[iserie][0] - 1;
            to   = lists.atm_ranges[iserie][1] - 1;

            if ( from < Low )
                from = Low;
            if ( to >= High )
                to = High - 1;

            for ( i = from ; i <= to ; i++ ) {

                res_hit = 0;
                seg_hit = 0;

/* now bang it together  .... */

                Count = gomp_DataVectorGetSize(&lists.seg_list);
                switch ( Count ) {

                case 0:
                    break;

                case 1:
                    if ( Tcl_StringMatch(SegPointer[i],lists.seg_list) )
                        seg_hit = 1;
                    break;

                default:
                    for ( j = 0 ; j < Count ; j++ ) {
                        if ( Tcl_StringMatch(SegPointer[i],
                                             lists.seg_list+
                                             j*(MAX_SEG_NAME_LEN+1)) ) {
                            seg_hit = 1;
                            break;
                        }
                    }
                }

                Count = gomp_DataVectorGetSize(&lists.res_list);
                switch ( Count ) {
                    
                case 0:
                    break;

                case 1:
                    if ( Tcl_StringMatch(ResPointer[i],lists.res_list) ) 
                        res_hit = 1;
                    break;

                default:
                    for ( j = 0 ; j < Count ; j++ ) {
                        if ( Tcl_StringMatch(ResPointer[i],
                                             lists.res_list+
                                             j*(MAX_RES_NAME_LEN+1)) ) {
                            res_hit = 1;
                            break;
                        }
                    }
                }
            
                Count = gomp_DataVectorGetSize(&lists.res_ranges);
                switch ( Count ) {

                case 0:
                    break;

                case 1:
                    if ( Num1Pointer[i] >= lists.res_ranges[0][0] && 
                         Num1Pointer[i] <= lists.res_ranges[0][1] )
                        res_hit = 1;
                    break;

                default:
                    for ( j = 0 ; j < Count ; j++ ) {
                        if ( Num1Pointer[i] >= lists.res_ranges[j][0] && 
                             Num1Pointer[i] <= lists.res_ranges[j][1] ) {
                            res_hit = 1;
                            break;
                        }
                    }
                }

                if ( seg_hit && res_hit ) {
                    if ( ( sel_long == 0 || sel_list[sel_long-1] != i ) &&
                         sel_long < High ) {
                        sel_list[sel_long] = i;
                        sel_long++;
                    }
                }
            }
        }
    }
    else if ( gomp_DataVectorGetSize(&lists.res_list) == 0 ) {

        nseries = gomp_DataVectorGetSize(&lists.res_ranges);
        for ( iserie = 0 ; iserie < nseries ; iserie++ ) {

            for( i = Low ; i < High ; i++ ) {

                if( Num1Pointer[i] < lists.res_ranges[iserie][0] ||
                    Num1Pointer[i] > lists.res_ranges[iserie][1] )
                    continue;

                atm_hit = 0;
                seg_hit = 0;

/* now bang it together  .... */
                Count = gomp_DataVectorGetSize(&lists.seg_list);
                switch ( Count ) {

                case 0:
                    break;

                case 1:
                    if ( Tcl_StringMatch(SegPointer[i],lists.seg_list) )
                        seg_hit = 1;
                    break;

                default:
                    for ( j = 0 ; j < Count ; j++ ) {
                        if ( Tcl_StringMatch(SegPointer[i],
                                             lists.seg_list+
                                             j*(MAX_SEG_NAME_LEN+1)) ) {
                            seg_hit = 1;
                            break;
                        }
                    }
                }

                Count = gomp_DataVectorGetSize(&lists.atm_list);
                switch ( Count ) {

                case 0:
                    break;

                case 1:
                    if ( Tcl_StringMatch(AtmPointer[i],lists.atm_list) )
                        atm_hit = 1;
                    break;

                default:
                    for ( j = 0 ; j < Count ; j++ ) {
                        if ( Tcl_StringMatch(AtmPointer[i],
                                             lists.atm_list+
                                             j*(MAX_ATM_NAME_LEN+1)) ) {
                            atm_hit = 1;
                            break;
                        }
                    }
                }

                Count = gomp_DataVectorGetSize(&lists.atm_ranges);
                switch ( Count ) {

                case 0:
                    break;

                case 1:
                    if ( i >= lists.atm_ranges[0][0] - 1 &&
                         i <= lists.atm_ranges[0][1] - 1 )
                        atm_hit = 1;
                    break;

                default:
                    for ( j = 0 ; j < Count ; j++ ) {
                        if ( i >= lists.atm_ranges[j][0] - 1 &&
                             i <= lists.atm_ranges[j][1] - 1 ) {
                            atm_hit = 1;
                            break;
                        }
                    }
                }

                if ( seg_hit && atm_hit ) {
                    if ( ( sel_long == 0 || sel_list[sel_long-1] != i ) &&
                         sel_long < High ) {
                        sel_list[sel_long] = i;
                        sel_long++;
                    }
                }
            }
        }
    }
    else {

        for ( i = Low ; i < High ; i++ ) {

            atm_hit = 0;
            res_hit = 0;
            seg_hit = 0;

/* now bang it together  .... */

            Count = gomp_DataVectorGetSize(&lists.seg_list);
            switch ( Count ) {

            case 0:
                break;

            case 1:
                if ( Tcl_StringMatch(SegPointer[i],lists.seg_list) )
                    seg_hit = 1;
                break;

            default:
                for ( j = 0 ; j < Count ; j++ ) {
                    if ( Tcl_StringMatch(SegPointer[i],
                                         lists.seg_list+
                                         j*(MAX_SEG_NAME_LEN+1)) ) {
                        seg_hit = 1;
                        break;
                    }
                }
            }

            Count = gomp_DataVectorGetSize(&lists.res_list);
            switch ( Count ) {

            case 0:
                break;

            case 1:
                if ( Tcl_StringMatch(ResPointer[i],lists.res_list) )
                    res_hit = 1;
                break;

            default:
                for ( j = 0 ; j < Count ; j++ ) {
                    if ( Tcl_StringMatch(ResPointer[i],
                                         lists.res_list+
                                         j*(MAX_RES_NAME_LEN+1)) ) {
                        res_hit = 1;
                        break;
                    }
                }
            }
            
            Count = gomp_DataVectorGetSize(&lists.res_ranges);
            switch ( Count ) {

            case 0:
                break;

            case 1:
                if ( Num1Pointer[i] >= lists.res_ranges[0][0] &&
                     Num1Pointer[i] <= lists.res_ranges[0][1] )
                    res_hit = 1;
                break;

            default:
                for ( j = 0 ; j < Count ; j++ ) {
                    if ( Num1Pointer[i] >= lists.res_ranges[j][0] &&
                         Num1Pointer[i] <= lists.res_ranges[j][1] ) {
                        res_hit = 1;
                        break;
                    }
                }
            }

            Count = gomp_DataVectorGetSize(&lists.atm_list);
            switch ( Count ) {

            case 0:
                break;

            case 1:
                if ( Tcl_StringMatch(AtmPointer[i],lists.atm_list) )
                    atm_hit = 1;
                break;

            default:
                for ( j = 0 ; j < Count ; j++ ) {
                    if ( Tcl_StringMatch(AtmPointer[i],
                                         lists.atm_list+
                                         j*(MAX_ATM_NAME_LEN+1)) ) {
                        atm_hit = 1;
                        break;
                    }
                }
            }

            Count = gomp_DataVectorGetSize(&lists.atm_ranges);
            switch ( Count ) {

            case 0:
                break;

            case 1:
                if ( i >= lists.atm_ranges[0][0] - 1 &&
                     i <= lists.atm_ranges[0][1] - 1 )
                    atm_hit = 1;
                break;

            default:
                for ( j = 0 ; j < Count ; j++ ) {
                    if ( i >= lists.atm_ranges[j][0] - 1 &&
                         i <= lists.atm_ranges[j][1] - 1 ) {
                        atm_hit = 1;
                        break;
                    }
                }
            }

            if ( seg_hit && res_hit && atm_hit ) {
                sel_list[sel_long] = i;
                sel_long++;
            }
        }
    }

    FreeSegResAtmListsAndRanges(&lists);
    return(sel_long); /* return with the length of the select list */
}       
/***********************************************************************/
int gomp_Sign_char(char *input , int HowLong) 
    /* look for significant characters */

    /* The input string can contain following characters:
       ? and * where "?" means that the position can be
       any character and "*" which means that the rest
       of the string can be any characters.
       The input string is max characters long
       depending on if it is a segment, residue or
       atom name */
/***********************************************************************/
{

    register int i;
    register int ast;
    register int qmk;
    register int ast_pos;
    register char mask;
    static   int slong;

    /* check first if there are any "*" and "?" in the string */

    ast  =   0;
    qmk  =   0;
    mask =   1;

    slong = strlen(input); /* input string is 'long' characters long */
    if(slong > HowLong) slong = HowLong;

    for(i = 0 ; i < slong ; i++) {
        if(input[i] == '*') ast++;
        if(input[i] == '?') qmk++;
    }
    if(slong < HowLong) {
        for(i = slong ; i < HowLong ; i++) input[i] = ' ';
        slong = HowLong;
    }

    if(ast == 0 && qmk == 0) return(0);

    if(ast) { /* interpret the "*" */
        for(i = 0 ; i < slong ; i++ ) if(input[i] == '*') break;
    }

    ast_pos = i; 
    input[i] = mask;
    if(slong-(i+1)) {
        for(i = (ast_pos+1) ; i < slong ; i++)  {

            if(!isalnum(input[i])) input[i] = mask;
        }

    }

    if(qmk == 0) return(0);

    if(qmk) { /* interpret the "?" */

        for(i = 0 ; i < slong ; i++) {
            if(input[i] == '?') input[i] = mask;
        }
    }

    return(0);
}

/*
  If Strings are equal     return 1 
  not equal return 0
*/

/***********************************************************************/
int gomp_CompareStrings(const char *input1, const char *input2 , int HowLong)
/***********************************************************************/
{
    register  int i;
    register  char mask;

    mask = 1;

    for(i = 0 ; i < HowLong ; i++)  {
        if(input1[i] == mask || input2[i] == mask) continue;
        if(input1[i] != input2[i]) return(0);
    }
    return(1);
}

/************************************************************************/
static int gomp_MakeSelectionList_xyz(
        int   Which,        /* structure number    */
        float xc,           /* Center for the      */
        float yc,           /* search              */
        float zc,
        float rad,          /* Max distance        */
        const char *Segment,/* segment name (list) */
        const char *Residue,/* resideu name (list) */
        const char *Atom,   /* atom name (list)    */
        int   AcceptPoint,  /* Can distance be 0   */                   
        int *sel_list)     /* index list to atoms */
/************************************************************************/
{
    static int   i,slong,*sel_target;
    static float dist2,rad2;
    static const float *x,*y,*z;
    float xb,yb,zb;

    slong  = gomp_MakeSelectionList(Which , Segment , Residue , Atom , sel_list);

    rad2     = rad*rad;

    x        = gomp_GetAtomXCoordPointer(Which);
    y        = gomp_GetAtomYCoordPointer(Which);
    z        = gomp_GetAtomZCoordPointer(Which);

    sel_target = sel_list;

    for(i = 0 ; i < slong ; i++) {
        xb = x[sel_list[i]];
        yb = y[sel_list[i]];
        zb = z[sel_list[i]];

        dist2 = (xb-xc)*(xb-xc) + (yb-yc)*(yb-yc) + (zb-zc)*(zb-zc);
        
        if( dist2 < rad2 ) {
            if( AcceptPoint || dist2 > 0 )
                *sel_target++ = sel_list[i];
        }
    }

    return(sel_target-sel_list);
}
/************************************************************************/
static int gomp_MakeSelectionListAround(
        int   Which,        /* structure number    */
        const char *SegmentBasicSet, /* center subset   */
        const char *ResidueBasicSet, /* for the search  */
        const char *AtomBasicSet,
        float rad,          /* Max distance        */
        const char *Segment,/* segment name (list) */
        const char *Residue,/* resideu name (list) */
        const char *Atom,   /* atom name (list)    */
        int   AcceptCenter, /* Whether include center */                   
        int *sel_list)     /* index list to atoms */
/************************************************************************/
{
    static int   i,j,k;
    static int   NAtoms,slong_basic,slong;
    static int *basic_sel_list,*sel_target;
    static float rad2;
    static const float *x,*y,*z;
    static float xb,yb,zb;
    static float xc,yc,zc;

    NAtoms = gomp_GetNumAtomsInMolecStruct(Which);

    basic_sel_list = gomp_AllocateIntVector(NAtoms);

    slong_basic = gomp_MakeSelectionList(Which ,
                                       SegmentBasicSet , ResidueBasicSet , AtomBasicSet , basic_sel_list);

    if( slong_basic <= 0 ) {
        gomp_PrintMessage("?ERROR - no atoms in the basic set");
        return(0);
    }

    slong  = gomp_MakeSelectionList(Which , Segment , Residue , Atom , sel_list);

    rad2     = rad*rad;

    x        = gomp_GetAtomXCoordPointer(Which);
    y        = gomp_GetAtomYCoordPointer(Which);
    z        = gomp_GetAtomZCoordPointer(Which);

    sel_target = sel_list;

    for(i = 0 ; i < slong ; i++) {
        xc = x[sel_list[i]];
        yc = y[sel_list[i]];
        zc = z[sel_list[i]];

        for(j = 0; j < slong_basic ; j++) {
            xb = x[basic_sel_list[j]];
            yb = y[basic_sel_list[j]];
            zb = z[basic_sel_list[j]];

            if(   (xb-xc)*(xb-xc)
                  + (yb-yc)*(yb-yc)
                  + (zb-zc)*(zb-zc) < rad2 ) {

                if( !AcceptCenter ) {
                    for(k = 0; k < slong_basic ; k++) {
                        if( sel_list[i] == basic_sel_list[k] )
                            break;
                    }
                    if( k < slong_basic )
                        continue;
                }

                *sel_target++ = sel_list[i];
                break;
            }
        }
    }

    free(basic_sel_list);

    return(sel_target-sel_list);
}

/*
  If alt > 0 add the atoms to the display list
  alt < 1 remove the atoms from the display list
*/
/************************************************************************/
int gomp_ParseDisplayList(int Display,
                          const gom_SelectionList *atoms)
/************************************************************************/
{
    register int i;
    static int StruL;
    static int ihelp,atom_list;
    static int *sel_list;
    static int slong;
    static char *disp_list;
    static char chelp[BUFF_LEN];
    static int  Wstr;

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);

    for(StruL = 0 ; StruL < gomp_GetNumMolecStructs() ; StruL++) {

        if(!gomp_GetSelectedStructure(StruL))
            continue;

        slong     = gomp_MakeSelectionList(StruL ,
                                         atoms->Segment, atoms->Residue, atoms->Atom, 
                                         sel_list);
        ihelp     = 1;

        if(slong < 0)
            continue;

        if(slong > 0) {      /* #1 */
            Wstr      = StruL;
            disp_list = gomp_GetModifiableAtomDisplayStatePointer(Wstr);

/* check if the idea was to turn the selected atoms on */
            if(Display) { /* yes it was ... */
                for(i = 0 ; i < slong ; i++)
                    disp_list[sel_list[i]] = 1;
            }
            else {
                for(i = 0 ; i < slong ; i++)
                    disp_list[sel_list[i]] = 0;
            }
        }                   /* #1 */
        else {
            sprintf(chelp,"structure (%d): no atoms in the selection list",
                    (StruL+1));
            gomp_PrintMessage(chelp);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    free(sel_list);
    sel_list = NULL;

    return(0);
}
/************************************************************************/
int gomp_ParseTypeCPKList(int Show,
                          const gom_SelectionList *atoms)
/************************************************************************/
{

    register int  i;
    static int    ihelp,atom_list;
    static int *sel_list;
    static int    slong;
    static char *disp_list;
    static char   chelp[BUFF_LEN];
    static int    j;

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);

    for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {
         
        if(!gomp_GetSelectedStructure(j)) continue;

        slong     = gomp_MakeSelectionList(
            j ,
            atoms->Segment,atoms->Residue,atoms->Atom,sel_list);
        ihelp     = 1;

        if(slong < 0) continue;

        if(slong > 0) {      /* #1 */

            disp_list = gomp_GetModifiableAtomCPKDisplayStatePointer(j);

/* check if the idea was to turn the selected atoms on */
            if( Show ) { /* yes it was ... */
                for(i = 0 ; i < slong ; i++)
                    disp_list[sel_list[i]] = 1;
            }
            else {
                for(i = 0 ; i < slong ; i++)
                    disp_list[sel_list[i]] = 0;     
            }
        }                   /* #1 */
        else {
            sprintf(chelp,"structure (%d) no atoms in the selection list",
                    (j + 1));
            gomp_PrintMessage(chelp);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    if(sel_list) {
        free(sel_list);
        sel_list = NULL;
    }

    return(0);
}
/************************************************************************/
int gomp_ParseCPKScaleList(float scale,
                           const gom_SelectionList *atoms)
/************************************************************************/
{
    register int  i;
    static int    ihelp,atom_list;
    static int *sel_list;
    static int    slong;
    static char   chelp[BUFF_LEN];
    static int    j;

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);

    for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {

        if(!gomp_GetSelectedStructure(j))
            continue;

        slong     = gomp_MakeSelectionList(
            j ,
            atoms->Segment, atoms->Residue, atoms->Atom, sel_list);
        ihelp     = 1;

        if(slong < 0)
            continue;

        if(slong > 0) {      /* #1 */
            for ( i = 0 ; i < slong ; i++ )
                gomp_PutAtomCPKScale(j, scale, sel_list[i]);
        }                   /* #1 */
        else {
            sprintf(chelp,"structure (%d) no atoms in the selection list",
                    (j + 1));
            gomp_PrintMessage(chelp);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    free(sel_list);

    return(0);
}
/************************************************************************/
int gomp_ParseTypeLicoList(int Show,
                           const gom_SelectionList *atoms)
/************************************************************************/
{
    register int i;
    static int ihelp,atom_list;
    static int *sel_list;
    static int slong;
    static char *disp_list;
    static char chelp[BUFF_LEN];
    static int  j;

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);

    for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {

        if(!gomp_GetSelectedStructure(j))
            continue;

        slong     = gomp_MakeSelectionList(
            j ,
            atoms->Segment, atoms->Residue, atoms->Atom, sel_list);
        ihelp     = 1;

        if(slong < 0)
            continue;

        if(slong > 0) {      /* #1 */

            disp_list = gomp_GetModifiableAtomLicoDisplayStatePointer(j);

/* check if the idea was to turn the selected atoms on */
            if(Show) { /* yes it was ... */
                for(i = 0 ; i < slong ; i++)
                    disp_list[sel_list[i]] = 1;
            }
            else {
                for(i = 0 ; i < slong ; i++)
                    disp_list[sel_list[i]] = 0;     
            }
        }                   /* #1 */
        else {
            sprintf(chelp,"structure (%d) no atoms in the selection list",
                    (j + 1));
            gomp_PrintMessage(chelp);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    free(sel_list);

    return(0);
}

/************************************************************************/
int gomp_ParseColourList(const gom_SelectionList *atoms,
                         const gom_FloatColour   *colour)
/************************************************************************/
{
    register int i;
    static int ihelp,atom_list;
    static int *sel_list;
    static int slong;
    static char chelp[BUFF_LEN];
    static int  j;

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);

    for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {

        if(!gomp_GetSelectedStructure(j))
            continue;

        slong     = gomp_MakeSelectionList(
            j ,
            atoms->Segment, atoms->Residue, atoms->Atom, sel_list);
        ihelp     = 1;

        if(slong < 0)
            continue;

        if(slong > 0) {      /* #1 */

            for(i = 0 ; i < slong ; i++)
                gomp_PutAtomColour(
                    j, colour->red, colour->green, colour->blue, sel_list[i]);
        }                   /* #1 */
        else {
            sprintf(chelp,"structure (%d) no atoms in the selection list",
                    (j + 1));
            gomp_PrintMessage(chelp);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    free(sel_list);

    return(0);
}
/************************************************************************/
int gomp_ParseSelectionList(int Which ,
                            const char *Segment , 
                            const char *Residue , 
                            const char *Atom)
/************************************************************************/
{

    static int atom_list;
    static int *sel_list;
    static int slong;
    static char chelp[BUFF_LEN];

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);

    slong     = 
        gomp_MakeSelectionList(Which , Segment , Residue, Atom, sel_list);

    if(slong > 0) {      /* #1 */

        (void)gomp_UpdateSelectionList(Which , slong , sel_list);
        (void)gomp_CalcSelectionCenter();

    }                    /* #1 */
    else {
        sprintf(chelp,"structure (%d) no atoms in the selection list",
                (Which+1));
        gomp_PrintMessage(chelp);
    }

    if(SaveAtomHitList(slong , sel_list)) {
        gomp_PrintERROR("can't save selection list");
    }

    if(sel_list) {
        free(sel_list);
        sel_list = NULL;
    }

    return(0);
}

/************************************************************************/
int gomp_CalculateCenterofMass(const char *segement,
                               const char *residue,
                               const char *atom,
                               float *xc      , float *yc     , float *zc)
    /* calculate centre of mass */
/************************************************************************/
{

    static int     i,j,si;
    static int *sel_list;
    static int     slong;
    static int     ihelp,atom_max;
    static int     strucs;
    static float   mass1[3],mass2,mass3;
    static char    OutText[BUFF_LEN];
    static const float *x,*y,*z;
    static const char *disp_list;
    static int     Wstr;
    static const float *sumxyz;
     
/* check to see if there are already coordinates coming in ... */
    if(gomp_IsStringAFloat(segement) && gomp_IsStringAFloat(residue) && gomp_IsStringAFloat(atom)) {

        *xc = atof(segement);  
        *yc = atof(residue);  
        *zc = atof(atom);  
  
        return(0);
    }

    if(InitAtomHitList())
        return(1);

    atom_max = gomp_GetTotalNumberOfAtoms();
    sel_list = gomp_AllocateIntVector(atom_max);

    *xc    = 0.0;
    *yc    = 0.0;
    *zc    = 0.0;
    strucs = 0;

    for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {

        if(!gomp_GetSelectedStructure(j)) continue;

        Wstr = j;

        sprintf(OutText,
                "Centre of mass for structure # %d (x, y and z) calculated from selection list : >%.4s<>%.4s<>%.4s<",
                (j+1),segement,residue,atom);
        gomp_PrintMessage(OutText);

        slong = gomp_MakeSelectionList(j,segement,residue,atom,sel_list);

        if(slong < 0) continue;

        if(slong > 0) {
            strucs++;
            ihelp = 0;
            mass1[0] = mass1[1] = mass1[2] = 0.0;
            mass2 = 0.0;

            sumxyz    = gomp_GetTranslateArray();
            disp_list = gomp_GetAtomDisplayStatePointer(Wstr);
            x         = gomp_GetAtomXCoordPointer(Wstr);
            y         = gomp_GetAtomYCoordPointer(Wstr);
            z         = gomp_GetAtomZCoordPointer(Wstr);

            for(i = 0 ; i < slong ; i++) {
                si = sel_list[i];

/* now bang it together  .... */
                ihelp++;
                mass3     = gomp_GetAtomMass(Wstr , si);
                mass1[0] += mass3*x[si];
                mass1[1] += mass3*y[si];
                mass1[2] += mass3*z[si];
                mass2 +=mass3;
            }

            if(mass2 < 0.01) {
                gomp_PrintMessage("?ERROR - atomic masses are not defined or set ");
                return(1);
            }

            mass1[0] /= mass2;
            mass1[1] /= mass2;
            mass1[2] /= mass2;

/*
  sprintf(OutText," %f  %f  %f ",
  mass1[0]+sumxyz[0],mass1[1]+sumxyz[1],mass1[2]+sumxyz[2]);
  gomp_PrintMessage(OutText);
  sprintf(OutText,"Total mass is: %f for %d atoms",mass2,ihelp);
  gomp_PrintMessage(OutText);
*/
            *xc = *xc + mass1[0]+sumxyz[0];
            *yc = *yc + mass1[1]+sumxyz[1];
            *zc = *zc + mass1[2]+sumxyz[2];
        }
        else {
            sprintf(OutText,"no atoms in the selection list for struct # %d",(j+1));
            gomp_PrintMessage(OutText);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    if(strucs) {
        *xc = *xc / strucs;
        *yc = *yc / strucs;
        *zc = *zc / strucs;
    }

    free(sel_list);

    return(0);
}

/************************************************************************/
int gomp_CalcCoordinateCenter(const char *Segment,
                              const char *Residue,
                              const char *Atom ,
                              float *xc     ,float *yc     ,float *zc)
    /* calculate coordinate centre */
/************************************************************************/
{
    static int     i,j,si;
    static int *sel_list;
    static int     slong;
    static int     ihelp,atom_max;
    static int     strucs;
    static float   mass1[3];
    static char    OutText[BUFF_LEN];
    static const float *x,*y,*z;
    static int     Wstr;
    static const float *sumxyz;


/* check to see if there are already coordinates coming in ... */
    if(gomp_IsStringAFloat(Segment) && gomp_IsStringAFloat(Residue) && gomp_IsStringAFloat(Atom)) {

        *xc = atof(Segment);  
        *yc = atof(Residue);  
        *zc = atof(Atom);  
  
        return(0);
    }

    if(InitAtomHitList())
        return(1);

    atom_max = gomp_GetTotalNumberOfAtoms();
    sel_list = gomp_AllocateIntVector(atom_max);

    *xc    = 0.0;
    *yc    = 0.0;
    *zc    = 0.0;
    strucs = 0;

    for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {

        if(!gomp_GetSelectedStructure(j)) continue;

        Wstr = j;

        sprintf(OutText,"Coordinate Centre for structure # %d (x, y and z) calculated from selection list : >%.4s<>%.4s<>%.4s<",
                (j+1),Segment,Residue,Atom);
        gomp_PrintMessage(OutText);

        slong = gomp_MakeSelectionList(j,Segment,Residue,Atom,sel_list);

        if(slong < 0) continue;

        if(slong > 0) {
            strucs++;
            ihelp = 0;
            mass1[0] = mass1[1] = mass1[2] = 0.0;

            sumxyz    = gomp_GetTranslateArray();
            x         = gomp_GetAtomXCoordPointer(Wstr);
            y         = gomp_GetAtomYCoordPointer(Wstr);
            z         = gomp_GetAtomZCoordPointer(Wstr);

            for(i = 0 ; i < slong ; i++) {
                si = sel_list[i];

/* now bang it together  .... */
                ihelp++;

                mass1[0] += x[si];
                mass1[1] += y[si];
                mass1[2] += z[si];
            }

            mass1[0] /= (float)slong;
            mass1[1] /= (float)slong;
            mass1[2] /= (float)slong;

/*
  sprintf(OutText," %f  %f  %f ",
  mass1[0]+sumxyz[0],mass1[1]+sumxyz[1],mass1[2]+sumxyz[2]);
  gomp_PrintMessage(OutText);
*/
            *xc = *xc + mass1[0]+sumxyz[0];
            *yc = *yc + mass1[1]+sumxyz[1];
            *zc = *zc + mass1[2]+sumxyz[2];
        }
        else {
            sprintf(OutText,"no atoms in the selection list for struct # %d",(j+1));
            gomp_PrintMessage(OutText);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }
    
    if(strucs) {
        *xc = *xc / strucs;
        *yc = *yc / strucs;
        *zc = *zc / strucs;
    } else {
        gomp_PrintERROR("no atoms at all matched");
        free(sel_list);
        return(1);
    }
    free(sel_list);

    return(0);
}

/************************************************************************/
int ShowAtomsAround(int Display,
                         int ListControl,
                         int pindex, const int *sel_list, int sel_list_long,
                         float rad,
                         const gom_SelectionList *atoms,
                         const gom_FloatColour   *colour)
/*
  int pindex;              index to the select list (sel_list)
  const int *sel_list;           selection list
  int sel_list_long;       length of selection list 
  float rad;               radius for the search 
*/
/************************************************************************/
{

    static int *disp_res,*disp_res_x,*disp_res_y,*Comp_set,Comp_set_long;
    static int     i,j;
    static int     disp_atoms = 0;
    static float   xf,yf,zf,rad2;
    static int     resij,from,atom_max;
    static const float *x,*y,*z;
    static char *disp_list;
    static const segment_name_t *segment;
    static const residue_name_t *residue;
    static int     Wstr;
    static int    AtomsFirstRound;

    if ( ListControl == DEL_SELECTION_LIST ) {
        /* clean up the space ... */

        if(SaveAtomHitList(collect_selected_atoms , disp_res_y)) {
            gomp_PrintERROR("can't save selection list");
        }

        free(disp_res);
        free(disp_res_x);
        free(disp_res_y);
        free(Comp_set);

        disp_res   = NULL;
        disp_res_x = NULL;
        disp_res_y = NULL;
        Comp_set   = NULL;

        return(0);
    }

    Wstr      = 0;
    disp_list = gomp_GetModifiableAtomDisplayStatePointer(Wstr);
    x         = gomp_GetAtomXCoordPointer(Wstr);
    y         = gomp_GetAtomYCoordPointer(Wstr);
    z         = gomp_GetAtomZCoordPointer(Wstr);
    segment   = gomp_GetAtomSegNamePointer(Wstr);
    residue   = gomp_GetAtomResNamePointer(Wstr);

    from      = sel_list[pindex];
    xf        = x[from];
    yf        = y[from];
    zf        = z[from];
    rad2      = rad*rad;

    if ( ListControl == NEW_SELECTION_LIST ) {
        /* create new list ... */
        if(disp_res)
            free(disp_res);
        atom_max    = gomp_GetNumAtomsInMolecStruct(Wstr);
        disp_res    = gomp_AllocateIntVector(atom_max);
        disp_res_x  = gomp_AllocateIntVector(atom_max);
        disp_res_y  = gomp_AllocateIntVector(atom_max);

        for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
            disp_res[i]   = OFF;   /* first all off ... */
            disp_res_x[i] = OFF;   /* first all off ... */
            disp_res_y[i] = OFF;   /* first all off ... */
        }

        free(Comp_set);
        Comp_set      = gomp_AllocateIntVector(atom_max);
        Comp_set_long = gomp_MakeSelectionList(
            Wstr,atoms->Segment,atoms->Residue,atoms->Atom,Comp_set);

        if(Comp_set_long < 1) {
            gomp_PrintMessage("?ERROR - no atoms in the comparison set");
            
            free(disp_res);
            free(disp_res_x);
            free(disp_res_y);
            free(Comp_set);
            
            disp_res   = NULL;
            disp_res_x = NULL;
            disp_res_y = NULL;
            Comp_set   = NULL;
            
            return(1);
        }
    }

    AtomsFirstRound = 0;

    for( i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
        if(disp_res[i] || (i == from))
            continue;              /* first level     */

        if(   ( (x[i]-xf)*(x[i]-xf) 
                + (y[i]-yf)*(y[i]-yf)
                + (z[i]-zf)*(z[i]-zf))  < rad2 ) {

            resij = gomp_GetAtomResNum1(Wstr , i);
/* check the selection mode */
            if(gomp_GetAtomSelectionMode( ) == RESIDUE_SELECTION) {
                for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(Wstr) ; j++) {
                    if((strcmp(segment[i],segment[j])==0) && 
                       (strcmp(residue[i],residue[j])==0) &&
                       (gomp_GetAtomResNum1(Wstr , j) == resij)) {
                        disp_res_x[AtomsFirstRound] = j;
                        disp_res[j]                 = ON;
                        AtomsFirstRound++;
                    }
                }
            }
            else {
                disp_res_x[AtomsFirstRound] = i;
                disp_res[i]                 = ON;
                AtomsFirstRound++;
            }
        }
    }

/* now I know how many atoms 'AtomsFirstRound' is inside the search distance 
   then I have to see if these atoms are in the matching atom list
*/

    disp_atoms = 0;

    for(i = 0 ; i < AtomsFirstRound ; i++) {

        for(j = 0 ; j < Comp_set_long ; j++) {

            if(Comp_set[j] == disp_res_x[i]) {
                if(!Display) {
                    disp_list[disp_res_x[i]] = OFF;
                } else {
                    disp_list[disp_res_x[i]] = ON;
                }
                disp_res_y[collect_selected_atoms] = disp_res_x[i];
                Comp_set[j] = -1;
                if(Display && colour)
                    gomp_PutAtomColour(Wstr,
                                       colour->red,
                                       colour->green,
                                       colour->blue,
                                       disp_res_x[i]);
                collect_selected_atoms++;
                break;
            }
        }
    }

    return(0);
}

/************************************************************************/
int ShowAtomsAround_xyz(int Display,
                             float xc, float yc, float zc,
                             float rad,
                             const gom_SelectionList *atoms,
                             const gom_FloatColour   *colour)
/*
  float xc,yc,zc;                 center for the search 
  float rad;                      radius for the search  
*/
/************************************************************************/
{
    static int *disp_res,*disp_res_x,*disp_res_y,*Comp_set,Comp_set_long;
    static int     i,j,resij,atom_max;
    static int     disp_atoms = 0;
    static float   rad2;
    static const float *x,*y,*z;
    static char *disp_list;
    static int     Wstr;
    static const segment_name_t *segment;
    static const residue_name_t *residue;
    static int    AtomsFirstRound;

    atom_max    = gomp_GetNumAtomsInMolecStruct(Wstr);
    disp_res    = gomp_AllocateIntVector(atom_max);
    disp_res_x  = gomp_AllocateIntVector(atom_max);
    disp_res_y  = gomp_AllocateIntVector(atom_max);

    rad2      = rad*rad;

    disp_list = gomp_GetModifiableAtomDisplayStatePointer(Wstr);
    x         = gomp_GetAtomXCoordPointer(Wstr);
    y         = gomp_GetAtomYCoordPointer(Wstr);
    z         = gomp_GetAtomZCoordPointer(Wstr);
    segment   = gomp_GetAtomSegNamePointer(Wstr);
    residue   = gomp_GetAtomResNamePointer(Wstr);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++)
        disp_res[i] = 0;

    AtomsFirstRound = 0;

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

/* don't include atoms already selected */
        if(disp_res[i]) continue;

        if(   ( (x[i]-xc)*(x[i]-xc) 
                + (y[i]-yc)*(y[i]-yc)
                + (z[i]-zc)*(z[i]-zc) )  < rad2 ) {

            resij = gomp_GetAtomResNum1(Wstr , i);
/* check the selection mode */
            if(gomp_GetAtomSelectionMode( ) == RESIDUE_SELECTION) {
                for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(Wstr) ; j++) {
                    if((strcmp(segment[i],segment[j])==0) &&
                       (strcmp(residue[i],residue[j])==0) &&
                       (gomp_GetAtomResNum1(Wstr , j) == resij)) {
                        disp_res_x[AtomsFirstRound] = j;
                        disp_res[j]                 = ON;
                        AtomsFirstRound++;
                    }
                }
            } else { 
                disp_res_x[AtomsFirstRound] = i;
                disp_res[i]                 = ON;
                AtomsFirstRound++;
            }
        }
    }

    disp_atoms = 0;

    Comp_set      = gomp_AllocateIntVector(atom_max);
    Comp_set_long = gomp_MakeSelectionList(Wstr,
                                         atoms->Segment,atoms->Residue,atoms->Atom,Comp_set);

    if(Comp_set_long < 1) {
        gomp_PrintMessage("?ERROR - no atoms in the comparison set");
        free(disp_res);
        free(Comp_set);
        return(1);
    }

    for ( i = 0 ; i < AtomsFirstRound ; i++ ) {

        for(j = 0 ; j < Comp_set_long ; j++) {

            if(Comp_set[j] == disp_res_x[i]) {
                if(!Display) {
                    disp_list[disp_res_x[i]] = OFF;
                } else {
                    disp_list[disp_res_x[i]] = ON;
                }
                disp_res_y[collect_selected_atoms] = disp_res_x[i];
                Comp_set[j] = -1;
                if(Display && colour)
                    gomp_PutAtomColour(Wstr,
                                       colour->red,
                                       colour->green,
                                       colour->blue,
                                       disp_res_x[i]);
                collect_selected_atoms++;
                break;
            }
        }
    }

    if(SaveAtomHitList(collect_selected_atoms , disp_res_y)) {
        gomp_PrintERROR("can't save selection list");
    }

    free(Comp_set);
    free(disp_res);
    free(disp_res_x);
    free(disp_res_y);

    return(0);
}
#if 0
/* 
   Calculate the centre of current displayed coordinates

   If alt != 0 Calculate new centre and make the translation

   If alt == 0 Calculate new centre but make no translation
*/
/***********************************************************************/
static const float *CalcCoordCentre(int alt)  
    /* Calculate centre of given coordinates  */
/***********************************************************************/
{
    int     i;
    float   numsum;
    int     from,to;
    char    OutText[BUFF_LEN];
    int     Wstr;
    const char *disp_list;
    float   Tsumx,Tsumy,Tsumz;
    float *x,*y,*z;
    static  float  sumxyz[3];

    Wstr = 0;

    Tsumx=0.0; Tsumy=0.0; Tsumz=0.0;
    numsum = 0.0;

    from      = 0;
    to        = gomp_GetNumAtomsInMolecStruct(Wstr);
    disp_list = gomp_GetAtomDisplayStatePointer(Wstr);
    x         = gomp_GetModifiableAtomXCoordPointer(Wstr);
    y         = gomp_GetModifiableAtomYCoordPointer(Wstr);
    z         = gomp_GetModifiableAtomZCoordPointer(Wstr);

    for(i = from ; i < to ; i++) {
        if(disp_list[i] == 0) continue;
        numsum = numsum + 1.;
        Tsumx = Tsumx + x[i]; 
        Tsumy = Tsumy + y[i]; 
        Tsumz = Tsumz + z[i]; }

    Tsumx = Tsumx/numsum; 
    Tsumy = Tsumy/numsum; 
    Tsumz = Tsumz/numsum;

    if(alt) {
        for(i = from ; i < to ; i++) {
            x[i] = x[i] - Tsumx; 
            y[i] = y[i] - Tsumy; 
            z[i] = z[i] - Tsumz;
        }

        sumxyz[0]   = Tsumy;
        sumxyz[1]  = Tsumy;
        sumxyz[2] = Tsumz;
    }

    if(gomp_GetDebugLevel()) {
        sprintf(OutText,"Translating ... x: %f , y: %f , z: %f \n",
                -Tsumx,-Tsumy,-Tsumz);
        gomp_PrintMessage(OutText);
    }

    return(sumxyz);
}
#endif

/************************************************************************/
int gomp_ControlSelectRoundAtoms(int Select,
                                 const gom_SelectionList *centre,
                                 float radius,
                                 const gom_SelectionList *atoms,
                                 const gom_FloatColour *colour)
/************************************************************************/
{
    static int    i;
    static int *sel_list;
    static int    slong;
    static char   chelp[BUFF_LEN];
    static char   chelp1[BUFF_LEN];
    static int    atom_list;
    static int    ListControl;
    static const float *TranslateXYZ;

    if(gomp_GetAtomSelectionMode( ) == RESIDUE_SELECTION) 
        gomp_PrintMessage("Selecting by whole residues ...");
    else
        gomp_PrintMessage("Selecting by atoms ...");

    if(InitAtomHitList())
        return(1);

    if(gomp_IsStringAFloat(centre->Segment)) {

/* given as x,y and z coordinates */

/* these can contain a colour */

        collect_selected_atoms = 0;

        TranslateXYZ = gomp_GetTranslateArray();

        (void)ShowAtomsAround_xyz(
            Select,
            atof(centre->Segment) - TranslateXYZ[0],
            atof(centre->Residue) - TranslateXYZ[1],
            atof(centre->Atom)    - TranslateXYZ[2],
            radius,
            atoms,
            colour);
    }
    else {
/* given as res:seg:atom */

/* these can contain a colour */

        atom_list = gomp_GetNumAtomsInMolecStruct(0);

        sel_list  = gomp_AllocateIntVector(atom_list);

        slong     = gomp_MakeSelectionList(
            0,
            centre->Segment,centre->Residue,centre->Atom,sel_list);

        if(slong > 0) {
            if(radius < 0.001) {
                gomp_PrintERROR("?ERROR - search radius too small");
                free(sel_list);
                return(1);
            }

            collect_selected_atoms = 0;

            ListControl = NEW_SELECTION_LIST;

            for(i = 0 ; i < slong ; i++ ) {

                (void)gomp_PutText2StatusLine2((10 * (i + 1))/slong);

                if ( ShowAtomsAround(
                         Select,ListControl,i,sel_list,slong,radius,
                         atoms,colour) )
                    break;
             
                ListControl = OLD_SELECTION_LIST;
            }

            (void)gomp_PutText2StatusLine2(0);

/* clean up the space ... */
            ListControl = DEL_SELECTION_LIST;
            (void)ShowAtomsAround(
                Select,ListControl,i,sel_list,slong,radius,
                atoms,colour);
            free(sel_list);
        }
        else
            gomp_PrintMessage("?ERROR - no atoms in the selection list");
    }

    sprintf(chelp1,"Selected %d atoms from the system",collect_selected_atoms);
    gomp_PrintMessage(chelp1);

    sprintf(chelp,"%d",collect_selected_atoms);
    (void)gomp_SendTclReturn(chelp);
    return(0);
}
#if 0
/************************************************************************/
int gomp_CalcAtomMassCenter(const char *seg, const char *res, const char *atm,
                            float *xcm, float *ycm, float *zcm)    
    /* calculate centre of mass
       const char *seg;      segment name 
       const char *res;      residue name 
       const char *atm;      atom name    
       const float *xcm;     x-coordinate of centre of mass
       const float *ycm;     y-           -"-               
       const float *zcm;     z-           -"-               


       Leif Laaksonen 1990
       Modified for gOpenMol 1995
    */
/************************************************************************/
{

    static int    i,si;
    static int *sel_list;
    static int    slong;
    static int    ihelp,atom_list;
    static float  mass1[3],mass2,mass3;
    static char   OutText[BUFF_LEN];
    static const float *x,*y,*z;
    static const float *sumxyz;
    static int    Wstr;

    atom_list = gomp_GetNumAtomsInMolecStruct(0);
    sel_list  = gomp_AllocateIntVector(atom_list);

    Wstr      = 0;
    x         = gomp_GetAtomXCoordPointer(Wstr);
    y         = gomp_GetAtomYCoordPointer(Wstr);
    z         = gomp_GetAtomZCoordPointer(Wstr);

    sumxyz    = gomp_GetTranslateArray();

    sprintf(OutText,"Centre of mass for structure # 1 calculated from selection list : >%.4s<>%.4s<>%.4s< ",seg,res,atm);
    gomp_PrintMessage(OutText);

    slong = gomp_MakeSelectionList(Wstr,seg,res,atm,sel_list);

    if(slong > 0) {
        ihelp = 0;
        mass1[0] = mass1[1] = mass1[2] = 0.0;
        mass2 = 0.0;

        for(i = 0 ; i < slong ; i++) {
            si = sel_list[i];

/* now bang it together  .... */
            ihelp++;
            mass3 = gomp_GetAtomMass(Wstr , si);
            mass1[0] += mass3*x[si];
            mass1[1] += mass3*y[si];
            mass1[2] += mass3*z[si];
            mass2 +=mass3;
        }

        if(mass2 < 0.01) {
            gomp_PrintMessage("?ERROR - atomic masses are not defined or set ");
            return(1);
        }

        mass1[0] /= mass2;
        mass1[1] /= mass2;
        mass1[2] /= mass2;

        sprintf(OutText,"Centre of mass is: x = %f , y = %f , z = %f",
                mass1[0]+sumxyz[0],mass1[1]+sumxyz[1],mass1[2]+sumxyz[2]);
        gomp_PrintMessage(OutText);
        sprintf(OutText,"Total mass is: %f for %d atoms",mass2,ihelp);
        gomp_PrintMessage(OutText);
    }
    else {
        gomp_PrintMessage("?ERROR - no atoms in the selection list ");
        free(sel_list);
        return(1);
    }

/* save centre of mass and return */
    *xcm = mass1[0];
    *ycm = mass1[1];
    *zcm = mass1[2];
/* .............................. */
    free(sel_list);
    return(0);
}

/************************************************************************/
int gomp_CalcAtomCoordCenter(const char *seg, const char *res, const char *atm,
                             float *xcm, float *ycm, float *zcm) 
    /* calculate coordinate centre 
       const char *seg;      segment name 
       const char *res;      residue name 
       const char *atm;      atom name    
       const float *xcm;     x-coordinate coordinate centre
       const float *ycm;     y-           -"-              
       const float *zcm;     z-           -"-               


       Leif Laaksonen 1990 
       Modified for gOpenMol 1995
    */
/************************************************************************/
{

    static int    i,si;
    static int *sel_list;
    static int    slong;
    static int    ihelp,atom_list;
    static float  mass1[3],mass2;
    static char   OutText[BUFF_LEN];
    static const float *x,*y,*z;
    static int    Wstr;
    static const float *sumxyz;

    Wstr = 0;
    atom_list = gomp_GetNumAtomsInMolecStruct(Wstr);
    sel_list  = gomp_AllocateIntVector(atom_list);

    x         = gomp_GetAtomXCoordPointer(Wstr);
    y         = gomp_GetAtomYCoordPointer(Wstr);
    z         = gomp_GetAtomZCoordPointer(Wstr);

    sumxyz    = gomp_GetTranslateArray();

    if(InitAtomHitList())
        return(1);

    sprintf(OutText,"Coordinate centre for structure # 1 calculated from selection list : >%.4s<>%.4s<>%.4s< ",seg,res,atm);
    gomp_PrintMessage(OutText);

    slong = gomp_MakeSelectionList(Wstr,seg,res,atm,sel_list);

    if(slong > 0) {
        ihelp = 0;
        mass1[0] = mass1[1] = mass1[2] = 0.0;
        mass2 = (float)slong;

        for(i = 0 ; i < slong ; i++) {
            si = sel_list[i];

/* now bang it together  .... */
            ihelp++;
            mass1[0] += x[si];
            mass1[1] += y[si];
            mass1[2] += z[si];
        }


        mass1[0] /= mass2;
        mass1[1] /= mass2;
        mass1[2] /= mass2;

        sprintf(OutText,"Coordinate center is: x = %f , y = %f , z = %f",
                mass1[0] + sumxyz[0],mass1[1]+sumxyz[2],mass1[2]+sumxyz[2]);
        gomp_PrintMessage(OutText);
/* save hitlist */
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }

    }
    else {
        gomp_PrintMessage("?ERROR - no atoms in the selection list ");
        free(sel_list);
        return(1);
    }

/* save coordinate centre and return */
    *xcm = mass1[0];
    *ycm = mass1[1];
    *zcm = mass1[2];
/* .............................. */
    if(sel_list) {
        free(sel_list);
        sel_list = NULL;
    }
    return(0);
}
#endif
/************************************************************************/
int gomp_ParsePlumberList(const char *Segment , 
                        const char *Residue , 
                        const char *Atom    ,
                        float Red     ,
                        float Green   ,
                        float Blue    ,
                        float FRad    ,
                        int   Type    ,
                        float Width   ,
                        float Thickness,
                        int   Glue)
/************************************************************************/
{

    static int   atom_list;
    static int *sel_list;
    static int   slong;
    static int   StruL;
    static char  chelp[BUFF_LEN];

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);
  
    for(StruL = 0 ; StruL < gomp_GetNumMolecStructs() ; StruL++) {

        if(!gomp_GetSelectedStructure(StruL)) continue;

        slong     = gomp_MakeSelectionList(StruL , Segment , Residue, Atom, 
                                         sel_list);

        if(slong < 0) continue;

        if(slong > 0) {      /* #1 */

            if(gomp_LoadPlumberAtoms(StruL , slong, sel_list, Red, Green, Blue,
                                   FRad, Type, Width, Thickness, Glue )) {
                free(sel_list);
                return(1);
            }

        }                    /* #1 */
        else {
            sprintf(chelp,"structure (%d): no atoms in the selection list",
                    (StruL+1));
            gomp_PrintMessage(chelp);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    if(sel_list) {
        free(sel_list);
        sel_list = NULL;
    }

    return(0);
}

/************************************************************************/
int gomp_ParseColourListByCharge(const gom_SelectionList *atoms,
                                 double *min_charge,double *max_charge)
/************************************************************************/
{
    static int   atom_list;
    static int *sel_list;
    static int   slong;
    char         chelp[BUFF_LEN];
    static int   StruL;

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);

    for(StruL = 0 ; StruL < gomp_GetNumMolecStructs() ; StruL++) {

        if(!gomp_GetSelectedStructure(StruL))
            continue;

        slong     = gomp_MakeSelectionList(StruL ,
                                         atoms->Segment,atoms->Residue,atoms->Atom,sel_list);

        if(slong < 0)
            continue;

        if(slong > 0) {      /* #1 */

            if(gomp_ColorByCharge(StruL , slong , sel_list,
                                min_charge,max_charge)) {
                free(sel_list);
                return(1);
            }
        }                    /* #1 */
        else {
            sprintf(chelp,"structure (%d) no atoms in the selection list",
                    (StruL+1));
            gomp_PrintMessage(chelp);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    free(sel_list);
    sel_list = NULL;

    return(0);
}

/************************************************************************/
int gomp_ParseColourListFourth(const gom_SelectionList *atoms)
/************************************************************************/
{
    static int   atom_list;
    static int *sel_list;
    static int   slong;
    static int   StruL;
    static char  chelp[BUFF_LEN];

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);

    for(StruL = 0 ; StruL < gomp_GetNumMolecStructs() ; StruL++) {

        if(!gomp_GetSelectedStructure(StruL))
            continue;

        slong     = gomp_MakeSelectionList(StruL ,
                                         atoms->Segment,atoms->Residue,atoms->Atom,sel_list);

        if(slong < 0)
            continue;

        if(slong > 0) {      /* #1 */

            if(gomp_ColorByFourth(StruL , slong , sel_list)) {
                free(sel_list);
                return(1);
            }
        }                   /* #1 */
        else {
            sprintf(chelp,"structure (%d): no atoms in the selection list",
                    (StruL+1));
            gomp_PrintMessage(chelp);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    free(sel_list);
    sel_list = NULL;

    return(0);
}

/************************************************************************/
int gomp_ParseColourListResidueNumber(const gom_SelectionList *atoms)
/************************************************************************/
{
    static int atom_list;
    static int *sel_list;
    static int slong;
    static int   StruL;
    static char  chelp[BUFF_LEN];

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);

    for(StruL = 0 ; StruL < gomp_GetNumMolecStructs() ; StruL++) {

        if(!gomp_GetSelectedStructure(StruL))
            continue;

        slong     = gomp_MakeSelectionList(
            StruL,
            atoms->Segment, atoms->Residue, atoms->Atom,
            sel_list);

        if(slong < 0)
            continue;

        if(slong > 0) {      /* #1 */

            if(gomp_ColorByResidueNumber(StruL , slong , sel_list)) {
                free(sel_list);
                return(1);
            }
        }                    /* #1 */
        else {
            sprintf(chelp,"structure (%d): no atoms in the selection list",
                    (StruL+1));
            gomp_PrintMessage(chelp);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    free(sel_list);

    return(0);
}


/************************************************************************/
int   gomp_ShowAtomNumber(int StruL ,
                          const gom_SelectionList *atom)
/************************************************************************/
{

    static int atom_list;
    static int *sel_list;
    static int slong;

    atom_list = gomp_GetNumAtomsInMolecStruct(StruL);

    sel_list  = gomp_AllocateIntVector(atom_list);
    slong     = gomp_MakeSelectionList(
        StruL, atom->Segment, atom->Residue, atom->Atom, sel_list);

    if(slong > 1) { 
        gomp_PrintERROR("more than one atom in the list, only one allowed");
        return(-1);
    }
    if(!slong) { 
        gomp_PrintERROR("no atom in the selection list");
        return(-1);
    }

    slong = sel_list[0];

    if(sel_list) {
        free(sel_list);
        sel_list = NULL;
    }

    return(slong);

}

/************************************************************************/
int gomp_ParseRDFList(const char *Segment1 , 
                    const char *Residue1 , 
                    const char *Atom1    ,
                    const char *Segment2 ,
                    const char *Residue2 ,
                    const char *Atom2    ,
                    const char *TRcut     ,
                    const char *TNbin)
/************************************************************************/
{

    static int    Nbin,atom_list;
    static int *sel_list;
    static int    slong;
    static float  Rcut;
    static float  CellA;
    static float  CellB;
    static float  CellC;
    static char   chelp[BUFF_LEN];
    static int    Retval;

    atom_list = gomp_GetNumAtomsInMolecStruct(0);

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    Retval    = 1;
    Rcut      = 10.0;
    Nbin      = 101;
    sel_list  = gomp_AllocateIntVector(atom_list);
    slong     = gomp_MakeSelectionList(0,Segment1 , Residue1, Atom1, sel_list);

    if(slong > 0) {      /* #1 */

        if(TRcut[0] != (char)NULL) 
            Rcut = atof(TRcut);

        if(TNbin[0] != (char)NULL)
            Nbin = atoi(TNbin);

        CellA = gomp_GetCellA();
        if(CellA > 1.0e+18) {
            gomp_PrintERROR("your dimension a for the cell is > 1.0e+18");
            return(1);
        }
        CellB = gomp_GetCellB();
        if(CellB > 1.0e+18) {
            gomp_PrintERROR("your dimension b for the cell is > 1.0e+18");
            return(1);
        }
        CellC = gomp_GetCellC();
        if(CellC > 1.0e+18) {
            gomp_PrintERROR("your dimension c for the cell is > 1.0e+18");
            return(1);
        }

        sprintf(chelp,"Structure # 1: Rcut =  %f, BoxL =  %g * %g * %g, Nbin =  %d",
                Rcut,CellA,CellB,CellC,Nbin);
        gomp_PrintMessage(chelp);

        Retval = gomp_CalcRDF(sel_list,slong , CellA , CellB , CellC ,Rcut,Nbin,
                            Segment2, Residue2, Atom2);
/* save hitlist */
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }                   /* #1 */
    else
        gomp_PrintMessage("?ERROR - no atoms in the selection list ");

    if(sel_list) {
        free(sel_list);
        sel_list = NULL;
    }

    return(Retval);
}
/************************************************************************/
int gomp_GetSelectedStructure(int Which)
/************************************************************************/
{
    return(Selected.List[Which]);
}
/************************************************************************/
int gomp_DeleteSelectedStructureList()
/************************************************************************/
{
    if(Selected.Length) {
        free(Selected.List);
        Selected.Length = 0;
    }

    return(0);
}

/************************************************************************/
int gomp_PushSelectedStructure(int State)
/************************************************************************/
{
    static int i;

    if(!Selected.Length) {
        Selected.List    = gomp_AllocateIntVector(1);
        Selected.List[0] = State;
        Selected.Length  = 1;
        Selected.State   = 1;
    }
    else {
        Selected.List = gomp_ReallocateIntVector(Selected.List , (Selected.Length + 1));
        Selected.List[Selected.Length] = State;
        Selected.Length++;

        Selected.State   = 1;

        for(i = 0 ; i < Selected.Length ; i++) {
            if(!Selected.List[i]) {
                Selected.State  = 0;
                break;
            }
        }

    }

    return(0);
}
/************************************************************************/
int gomp_ActivateSelectedStructure(int Which , int State)
/************************************************************************/
{
    static int i;

    Selected.List[Which] = State;
    Selected.State       = 1;

    for(i = 0 ; i < Selected.Length ; i++) {
        if(!Selected.List[i]) {
            Selected.State      = 0;
            break;
        }
    }

    return(0);
}    

/************************************************************************/
int gomp_GetSelectedStructureStatus()
/************************************************************************/
{
    return(Selected.State);
}
#if 0
/************************************************************************/
int gomp_GetAtomIndexFromName(const char *Input , int *Wstr , int *AtomI)
/************************************************************************/
{

    register int i;
    static int atom_list;
    static int *sel_list;
    static int slong;
    static int  j;
    static char Segment[BUFF_LEN];
    static char Residue[BUFF_LEN];
    static char Atom[BUFF_LEN];

    if(Input[0] == (char)NULL) {
        gomp_PrintERROR("can't parse an empty string");
        return(1);
    }

    i = sscanf(Input,"%d %s %s %s",Wstr,Segment,Residue,Atom);
    if(i != 4) {
        gomp_PrintERROR("wrong number of parameters in input string (4)");
        return(1);
    }
    if(*Wstr < 0 || *Wstr > gomp_GetNumMolecStructs()) {
        gomp_PrintERROR("structure number out of range");
        return(1);
    }

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list) {
        return(1);
    }

    *Wstr = *Wstr - 1;

    sel_list  = gomp_AllocateIntVector(atom_list);

    j = *Wstr;

    slong  = gomp_MakeSelectionList( j , Segment , Residue, Atom, sel_list);

    if(slong > 1) {
        gomp_PrintWARNING("selection list is > 1, taking first hit");
    }

    *AtomI = sel_list[0];

    free(sel_list);
    sel_list = NULL;

    return(0);
}
#endif
/************************************************************************/
int gomp_ParseColourListByVector(const gom_SelectionList *atoms,
                                 double *Vmin,double *Vmax)
/************************************************************************/
{
    static int atom_list;
    static int *sel_list;
    static int slong;
    char chelp[BUFF_LEN];
    static int StruL;

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);

    for(StruL = 0 ; StruL < gomp_GetNumMolecStructs() ; StruL++) {

        if(!gomp_GetSelectedStructure(StruL))
            continue;

        slong     = gomp_MakeSelectionList(StruL ,
                                         atoms->Segment, atoms->Residue, atoms->Atom, sel_list);

        if(slong < 0)
            continue;

        if(slong > 0) {      /* #1 */

            if(gomp_ColorByVector(StruL , slong , sel_list , Vmin , Vmax)) {
                free(sel_list);
                return(1);
            }
        }                    /* #1 */
        else {
            sprintf(chelp,"structure (%d) no atoms in the selection list",
                    (StruL+1));
            gomp_PrintMessage(chelp);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    free(sel_list);
    sel_list = NULL;

    return(0);
}

/************************************************************************/
int gomp_ParsePlotVectorList(  const char *Segment ,
                               const char *Residue ,
                               const char *Atom    ,
                               const char *Radius  ,
                               const char *Scale)
/************************************************************************/
{

    register int i;
    static int atom_list;
    static int *sel_list;
    static int *temp_list;
    static int slong;
    static int   StruL;
    char chelp[BUFF_LEN];

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);
    temp_list = gomp_GetModifiableVectorListArray();

    StruL     = gomp_GetVectorStructureIndex();

    slong     = gomp_MakeSelectionList(StruL , Segment , Residue, Atom, 
                                     sel_list);

    if(slong > 0) {      /* #1 */

        if(*Scale  != (char)NULL)
            (void)gomp_SetVectorScale(atof(Scale));
        if(*Radius != (char)NULL)
            (void)gomp_SetVectorRadius(atof(Radius));

        (void)gomp_SetVectorListLength(slong);

        for(i = 0 ; i < slong ; i++)
            temp_list[i] = sel_list[i];
    }                    /* #1 */
    else {
        sprintf(chelp,"structure (%d) no atoms in the selection list",
                (StruL+1));
        gomp_PrintMessage(chelp);
    }

    if(SaveAtomHitList(slong , sel_list)) {
        gomp_PrintERROR("can't save selection list");
    }

    if(sel_list) {
        free(sel_list);
        sel_list = NULL;
    }

    return(0);
}

/************************************************************************/
int gomp_SetAtomSelectionMode(int SMode)
/************************************************************************/
{
    gomAtomSelectionMode = SMode;

    return(0);
}

/************************************************************************/
int gomp_GetAtomSelectionMode()
/************************************************************************/
{
    return(gomAtomSelectionMode);
}


/************************************************************************/
int gomp_SetSelectionModeStatus(int SMode)
/************************************************************************/
{
    gomSelectionModeStatus = SMode;

    return(0);
}

/************************************************************************/
int gomp_GetSelectionModeStatus()
/************************************************************************/
{
    return(gomSelectionModeStatus);
}

/************************************************************************/
int gomp_ParseAtomLabelList(int Show,
                            const gom_SelectionList *atoms)
/************************************************************************/
{
    register int i;
    static int   atom_list;
    static int *sel_list;
    static char *label;
    static int  slong;
    static int  StruL;
    char chelp[BUFF_LEN];

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);

    for(StruL = 0 ; StruL < gomp_GetNumMolecStructs() ; StruL++) {

        if(!gomp_GetSelectedStructure(StruL))
            continue;

        label     = gomp_GetModifiableAtomLabelDisplayStatePointer(StruL);
        slong     = gomp_MakeSelectionList(
            StruL, atoms->Segment, atoms->Residue, atoms->Atom, sel_list);

        if(slong < 0)
            continue;

        if(slong > 0) {      /* #1 */
            if ( Show ) {
                for(i = 0 ; i < slong ; i++)
                    label[sel_list[i]] = 1;

                (void)gomp_SetPlotLabelState(ON);
            }
            else {
                for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(StruL) ; i++)
                    label[i] = 0;

                (void)gomp_SetPlotLabelState(OFF);
            }
        }/* #1 */
        else {
            sprintf(chelp,"structure (%d) no atoms in the selection list",
                    (StruL+1));
            gomp_PrintMessage(chelp);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    free(sel_list);

    return(0);
}
/************************************************************************/
int gomp_ResetAtmConn(int allStructs)
/************************************************************************/
{
    static int Wstr;

    for(Wstr = 0 ; Wstr < gomp_GetNumMolecStructs() ; Wstr++) {

        if(!allStructs && !gomp_GetSelectedStructure(Wstr)) continue;

        if(gomp_CalcAtomConn(Wstr))
            return(1);
    }
    return(0);
}
/************************************************************************/
int gomp_ParseCalcConnList(const char *Segment , 
                           const char *Residue , 
                           const char *Atom    ,
                           const char *Extra)
/************************************************************************/
{

    register int i;
    static int StruL;
    static int atom_list;
    static int *sel_list;
    static int slong;
    static char chelp[BUFF_LEN];
    static int  Wstr;

    if((strcmp(Segment,"") == 0 || strcmp(Segment,"*") == 0) &&
       (strcmp(Residue,"") == 0 || strcmp(Residue,"*") == 0) &&
       (strcmp(Atom,"") == 0 || strcmp(Atom,"*") == 0))
        return gomp_ResetAtmConn(0);

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list) {
        return(1);
    }

    if(InitAtomHitList())
        return(1);

    sel_list  = gomp_AllocateIntVector(atom_list);

    for(StruL = 0 ; StruL < gomp_GetNumMolecStructs() ; StruL++) {

        if(!gomp_GetSelectedStructure(StruL)) continue;

        slong     = gomp_MakeSelectionList(StruL , Segment , Residue, Atom, 
                                         sel_list);
        if(slong < 0) continue;

        if(slong > 0) {      /* #1 */

            Wstr      = StruL;

            for(i = 0 ; i < slong ; i++) {   
                if(gomp_CalcAtomConnI(Wstr , sel_list[i]))
                    return(1);
            }
        }                   /* #1 */
        else 
        {
            sprintf(chelp,"structure (%d) no atoms in the selection list",
                    (StruL+1));
            gomp_PrintMessage(chelp);
        }
        if(SaveAtomHitList(slong , sel_list)) {
            gomp_PrintERROR("can't save selection list");
        }
    }

    if(sel_list) {
        free(sel_list);
        sel_list = NULL;
    }

    return(0);
}

/************************************************************************/
int gomp_ParseCalcHbondsList(const char *Segment ,
                             const char *Residue ,
                             const char *Atom    ,
                             int   SearchHydrogens ,
                             const char *Extra)
/************************************************************************/
{

    static int   NAtoms;
    static int *sel_list;
    static int   slong;
    static int   Wstr;
    static float RedC;
    static float GreenC;
    static float BlueC;

    NAtoms = gomp_GetTotalNumberOfAtoms();

    if(!NAtoms)
        return(1);

    if(InitAtomHitList())
        return(1);

    if(*Extra != (char)NULL) {
        if(gomp_ColourName2RGB(Extra , &RedC , &GreenC , &BlueC)) {
            gomp_FormatERROR("can't resolve colour name '%s'",Extra);
            return(1);
        }
        (void)gomp_SetHBondColour(RedC , GreenC , BlueC);
    }

    sel_list = gomp_AllocateIntVector(NAtoms);

    for(Wstr = 0 ; Wstr < gomp_GetNumMolecStructs() ; Wstr++) {

        if(!gomp_GetSelectedStructure(Wstr)) continue;

        slong    = gomp_MakeSelectionList( Wstr,
                                         Segment, Residue, Atom, sel_list );

        if(slong<=0)
            gomp_FormatMessage(
                "structure (%d) no atoms in the selection list",(Wstr+1));

        /* Do always. */
        if(gomp_FindHbonds(Wstr, SearchHydrogens, slong, sel_list))
            return(1);
    }

    free(sel_list);

    return(0);
}
/************************************************************************/
int gomp_ParseCalcHbondSubset(const char *Segment ,
                              const char *Residue ,
                              const char *Atom)
/************************************************************************/
{
    static int   NAtoms;
    static int *sel_list;
    static int   slong;
    static int   Wstr;

    NAtoms = gomp_GetTotalNumberOfAtoms();

    if(!NAtoms)
        return(1);

    if(InitAtomHitList())
        return(1);

    sel_list = gomp_AllocateIntVector(NAtoms);

    for(Wstr = 0 ; Wstr < gomp_GetNumMolecStructs() ; Wstr++) {

        if(!gomp_GetSelectedStructure(Wstr)) continue;

        slong = gomp_MakeSelectionList(
            Wstr, Segment, Residue, Atom, sel_list );

        if( slong <= 0 )
            gomp_FormatMessage(
                "structure (%d) no atoms in the selection list",(Wstr+1));

        /* Do always */
        if(gomp_SetHbondSubset(Wstr, slong, sel_list))
            return(1);
    }

    return(0);
}

/************************************************************************/
char *gomp_MakeIndexList(int slong,const int *list,int add,
                         int AllowRanges,char separator)
/************************************************************************/
{
    static int    i,IsRange;
    static char   Text[BUFF_LEN];
    static char *tcl_list,*new_list;
    static size_t length,length2;

    if( slong == 0 ) {
        tcl_list  = gomp_AllocateCharVector(1);
        if(!tcl_list)
            return NULL;
        *tcl_list = '\0';
        return tcl_list;
    }

    /* add first index */
    sprintf(Text,"%d",list[0]+add);
    length = strlen(Text);
    tcl_list = gomp_AllocateCharVector(length+1);
    strcpy(tcl_list,Text);

    IsRange = 0;

    for( i=1; i<slong; i++ ) {

        if( AllowRanges ) {
            if( list[i] == list[i-1] + 1 ) {
                /* Start or continue index range */
                IsRange = 1;

                if( i < slong - 1 )
                    continue;

                sprintf(Text,"-%d",list[i]+add);
            }
            else if( list[i] == list[i-1] ) {
                if( i == slong - 1 && IsRange )
                    /* End index range */
                    sprintf(Text,"-%d",list[i]+add);
                else
                    /* Eliminate duplicates. */
                    continue;
            }
            else if( IsRange ) {
                sprintf(Text,"-%d%c%d",list[i-1]+add,separator,list[i]+add);
                IsRange = 0;
            }
            else
                sprintf(Text,"%c%d",separator,list[i]+add);
        }
        else
            sprintf(Text,"%c%d",separator,list[i]+add);

        length2  = length + strlen(Text);
        new_list = gomp_ReallocateCharVector(tcl_list,length2+1);
        if( new_list )
            tcl_list = new_list;
        else {
            free(tcl_list);
            return NULL;
        }
        strcpy(tcl_list+length,Text);
        length   = length2;
    }

    return tcl_list;
}

/************************************************************************/
int gomp_ParseListProteinChains(int Which , const char *Segment ,
                                const char *Residue)
/************************************************************************/
{
    static int   i,j,from,to;
    static int   Wstr,NAtoms,NResidue;
    static int   alpha_c,*alpha_list,*res_list;
    static char  Text[BUFF_LEN],*residues;
    static int   chains,*chain_list;
    static char *tcl_list;
    static const char *SegName;

    if( Which >= 0 ) {
        NAtoms = gomp_GetNumAtomsInMolecStruct(Which);
        from   = Which;
        to     = Which + 1;
    }
    else {
        NAtoms = gomp_GetTotalNumberOfAtoms();
        from   = 0;
        to     = gomp_GetNumMolecStructs();
    }

    alpha_list = gomp_AllocateIntVector(NAtoms);

    /* make TCL list from residue numbers of protein chains */
    tcl_list  = gomp_AllocateCharVector(1);
    *tcl_list = '\0';

    for( Wstr=from; Wstr<to; Wstr++ ) {
        
        if( Which == LIST_SELECTED_STRUCTURES &&
            !gomp_GetSelectedStructure(Wstr) )
            continue;

        res_list   = alpha_list;    /* use the same buffer although */
                                    /* they have nothing in common  */
    
        alpha_c = gomp_MakeSelectionList(Wstr , Segment , Residue, "CA", alpha_list);

        if(alpha_c <= 0)
            continue;

        chains   = gomp_ProteinChainList(Wstr , &chain_list, &alpha_list, &alpha_c);
        
        for( i=0; i<chains; i++ ) {
            /* Count the number of residues in a chain. */
            if( i < chains - 1 )
                NResidue = chain_list[i+1] - chain_list[i];
            else
                NResidue = alpha_c         - chain_list[i];

            /* Get segment name. */
            SegName = gomp_GetAtomSegName(Wstr,alpha_list[chain_list[i]]);

            /* Make a residue list. */
            for( j=0; j<NResidue; j++ )
                res_list[j] = gomp_GetAtomResNum1(Wstr,alpha_list[chain_list[i]+j]);

            residues = gomp_MakeIndexList(NResidue,res_list,0,1,',');
            sprintf(Text,"%s{{%s} {",(*tcl_list ? " " : ""),SegName);
            tcl_list = gomp_ReallocateCharVector(tcl_list,
                                     strlen(tcl_list)+strlen(Text)+strlen(residues)+3);
            strcat(tcl_list,Text);
            strcat(tcl_list,residues);
            strcat(tcl_list,"}}");
            free(residues);
        }

        free(chain_list);
    }
    
    free(alpha_list);
    
    gomp_SendTclReturn(tcl_list);

    free(tcl_list);

    return(0);
}

/************************************************************************/
int gomp_AppendSegmentList(char **pList,int Wstr,int slong,const int *sel_list)
/************************************************************************/
{
    static int    i;
    static const char *seg_name;
    static char *p,*new_list;
    static size_t seg_length,list_length;

    if( !*pList ) {
        *pList  = gomp_AllocateCharVector(1);
        if( !*pList )
            return(1);
        **pList = '\0';
    }

    list_length = strlen(*pList);

    for( i=0; i<slong; i++ ) {
        seg_name   = gomp_GetAtomSegName(Wstr,sel_list[i]);
        seg_length = strlen(seg_name);

        /* Check if it is already there. */
        for( p=*pList-1; (p=strstr(p+1,seg_name)); ) {
            if( (p==*pList || *(p-1)==',') &&
                (p[seg_length]=='\0' || p[seg_length]==',') )
                break;
        }
        if( p )
            continue;

        /* Append segment name */
        new_list = gomp_ReallocateCharVector(
            *pList,
            list_length + seg_length + (( list_length > 0 ) ? 2 : 1) );

        if( !new_list )
            return(1);
        *pList = new_list;

        if( list_length )
            (*pList)[list_length++] = ',';

        strcpy( *pList + list_length , seg_name );
        list_length += seg_length;
    }

    return(0);
}

static int CompareInts(const void *a,const void *b)
{
    return *(const int *)a-*(const int *)b;
}

/************************************************************************/
int gomp_InitSegmentResidueAtomList(SegmentResidueAtomList_t *List,
                                  int Which,
                                  int listType,
                                  int AllocSelList)
/************************************************************************/
{
    int failed;

    failed = 0;

    List->Which     = Which;
    List->listType  = listType;
    List->SortLists = 0;

    /* Check which structures have to walked throw, and */
    /* count how many atom they contain.                */
    if( List->Which >= 0 ) {
        List->NAtoms = gomp_GetNumAtomsInMolecStruct(List->Which);
        List->from   = List->Which;
        List->to     = List->Which;
    }
    else {
        List->NAtoms = gomp_GetTotalNumberOfAtoms();
        List->from   = 0;
        List->to     = gomp_GetNumMolecStructs() - 1;
    }

    List->slong_total = 0;

    List->sel_list  = NULL;
    List->atom_list = NULL;
    List->res_list  = NULL;
    List->seg_list  = NULL;

    switch( List->listType ) {
    case SELECTION_LIST:
        /* Do ATOM_LIST, RESIDUE_LIST and SEGMENT_LIST sections. */
    case ATOM_LIST:
        /* We need atom list. */
        List->atom_list = gomp_AllocateIntVector(List->NAtoms);
        if( !List->atom_list )
            failed = 1;

        if( List->listType != SELECTION_LIST)
            break;
    case RESIDUE_LIST:
        /* We need res_list. */
        List->res_list = gomp_AllocateIntVector(List->NAtoms);
        if( !List->res_list )
            failed = 1;

        if( List->listType != SELECTION_LIST)
            break;
    case SEGMENT_LIST:
        /* We don't need lists (yet). */
        break;
    }

    if( AllocSelList ) {
        /* Use existing list if we have one. */
        List->sel_list = List->atom_list ? List->atom_list : List->res_list;
        if( !List->sel_list )
        {
            /* Create new list for this purpose. */
            List->sel_list = List->atom_list = gomp_AllocateIntVector(List->NAtoms);
            if( !List->sel_list )
                failed = 1;
        }
    }

    if( failed ) {
        gomp_FreeSegmentResidueAtomList(List);
        return(1);
    }

    return(0);
}

/************************************************************************/
int gomp_PreSegmentResidueAtomListStructure(SegmentResidueAtomList_t *List,
                                          int Wstr)
/************************************************************************/
{
    if( List->Which == LIST_SELECTED_STRUCTURES &&
        !gomp_GetSelectedStructure(Wstr) )
        return(1);
    return(0);
}

/************************************************************************/
int gomp_PostSegmentResidueAtomListStructure(SegmentResidueAtomList_t *List,
                                           int Wstr,int slong,const int *sel_list)
/************************************************************************/
{
    int i,*target_list;

    if( slong <= 0 )
        return(0);

    switch( List->listType ) {
    case SELECTION_LIST:
        /* Do ATOM_LIST, RESIDUE_LIST and SEGMENT_LIST sections. */
    case ATOM_LIST:
        target_list = List->atom_list + List->slong_total;
        if( sel_list != target_list )
            /* sel_list isn't List->sel_list. */
            /* We have to copy.               */
            memcpy(target_list, sel_list, slong*sizeof(*sel_list));

        if( List->listType != SELECTION_LIST)
            break;
    case RESIDUE_LIST:
        /* Convert atom numbers to residue numbers. */
        target_list = List->res_list + List->slong_total;
        for( i=0; i<slong; i++ )
            target_list[i] = gomp_GetAtomResNum1(Wstr,sel_list[i]);

        if( List->listType != SELECTION_LIST)
            break;
    case SEGMENT_LIST:
        /* Find segment names. */
        if( gomp_AppendSegmentList(&List->seg_list,Wstr,slong,sel_list) )
            return(1);

        break;
    }

    if( List->slong_total > 0 )
        /* Multiple additions. */
        List->SortLists = 1;

    List->slong_total += slong;
    if( List->sel_list )
        List->sel_list += slong;

    return(0);
}

/************************************************************************/
char *gomp_GetSegmentResidueAtomList(SegmentResidueAtomList_t *List)
/************************************************************************/
{
    char *atom_list=NULL,*res_list=NULL,*seg_list,*tcl_list=NULL;
    int failed;

    failed = 0;

    switch( List->listType ) {
    case SELECTION_LIST:
        /* Do ATOM_LIST, RESIDUE_LIST and SEGMENT_LIST sections. */
    case ATOM_LIST:
        if( List->SortLists )
            qsort(List->atom_list, List->slong_total,
                  sizeof(*List->atom_list), CompareInts);
        /* Make index list. Increase indexes by one and allow ranges. */
        tcl_list = atom_list = gomp_MakeIndexList(
            List->slong_total,List->atom_list,1,1,',');
        if( !atom_list )
            failed = 1;

        if( List->listType != SELECTION_LIST)
            break;
    case RESIDUE_LIST:
        if( List->SortLists )
            qsort(List->res_list, List->slong_total,
                  sizeof(*List->res_list), CompareInts);
        /* Make index list. Allow ranges. */
        tcl_list = res_list = gomp_MakeIndexList(
            List->slong_total,List->res_list,0,1,',');
        if( !res_list )
            failed = 1;

        if( List->listType != SELECTION_LIST)
            break;
    case SEGMENT_LIST:
        tcl_list = seg_list = List->seg_list;
        List->seg_list      = NULL;
        if( !seg_list ) {
            tcl_list = seg_list = gomp_AllocateCharVector(1);
            if( !seg_list )
                failed = 1;
            else
                *seg_list='\0';
        }

        if( List->listType != SELECTION_LIST)
            break;

        /* Rest of the SELECTION_LIST section */
        if( !failed ) {
            tcl_list = gomp_AllocateCharVector(
                strlen(seg_list )+3+
                strlen(res_list )+3+
                strlen(atom_list)+3);

            if( tcl_list )
                sprintf(tcl_list,"{%s} {%s} {%s}",seg_list,res_list,atom_list);
            else
                failed = 1;
        }

        /* Some allocations may have succeeded. */
        free(atom_list);
        free(res_list );
        free(seg_list );
    }

    if(failed)
        return NULL;

    return tcl_list;
}

/************************************************************************/
int gomp_FreeSegmentResidueAtomList(SegmentResidueAtomList_t *List)
/************************************************************************/
{
    /* List->sel_list isn't an independed list. Don't free it. */
    free(List->atom_list);
    free(List->res_list);
    free(List->seg_list);

    List->atom_list = NULL;
    List->res_list  = NULL;
    List->seg_list  = NULL;

    return(0);
}

/************************************************************************/
int gomp_ReturnAndFreeSegmentResidueAtomList(SegmentResidueAtomList_t *List)
/************************************************************************/
{
    char *tcl_list;

    tcl_list = gomp_GetSegmentResidueAtomList(List);
    if( !tcl_list )
        return(1);

    gomp_SendTclReturn(tcl_list);
    free(tcl_list);

    return(0);
}

/************************************************************************/
int gomp_ParseListAtomsAround(int Which, const char *Segment,
                            const char *Residue,
                            const char *Atom,
                            const char *Radius,
                            const char *SegmentBasicSet,
                            const char *ResidueBasicSet,
                            const char *AtomBasicSet,
                            int   listType)
/************************************************************************/
{
    static int    Wstr;
    static int    AroundXYZ;
    static int    slong;
    static float  xc,yc,zc,rad;
    static const float *sumxyz;
    static SegmentResidueAtomList_t List;

    if( gomp_IsStringAFloat(SegmentBasicSet) &&
        gomp_IsStringAFloat(ResidueBasicSet) &&
        gomp_IsStringAFloat(AtomBasicSet) ) {

        AroundXYZ = 1;

        sumxyz = gomp_GetTranslateArray();

        xc = atof(SegmentBasicSet) - sumxyz[0];
        yc = atof(ResidueBasicSet) - sumxyz[1];
        zc = atof(AtomBasicSet)    - sumxyz[2];
    }
    else
        AroundXYZ = 0;

    rad = atof(Radius);
    if(rad < 0.001) {
        gomp_PrintERROR("?ERROR - search radius too small");
        return(1);
    }

    if( gomp_InitSegmentResidueAtomList(&List,Which,listType,1) )
        return(1);

    /* Loop over structures. */
    for( Wstr=List.from; Wstr<=List.to; Wstr++ ) {

        if( gomp_PreSegmentResidueAtomListStructure(&List,Wstr) )
            continue;

        if(AroundXYZ)
            slong = gomp_MakeSelectionList_xyz(
                Wstr , xc, yc, zc, rad,
                Segment, Residue, Atom, 1, List.sel_list);
        else {
            slong = gomp_MakeSelectionListAround(
                Wstr,
                SegmentBasicSet, ResidueBasicSet, AtomBasicSet, rad,
                Segment,         Residue,         Atom,
                0, List.sel_list);
        }

        (void)gomp_PostSegmentResidueAtomListStructure(
            &List,Wstr,slong,List.sel_list);
    }

    return(gomp_ReturnAndFreeSegmentResidueAtomList(&List));
}

/************************************************************************/
int gomp_ParseListAtomsJoin(int Which, int Sets, const char **argv,int listType)
/************************************************************************/
{
    static int   Wstr,i,slong;
    static const char *Segment,*Residue,*Atom;
    static SegmentResidueAtomList_t List;

    if( Sets <= 0 ) {
        gomp_SendTclReturn("");
        return(0);
    }

    if( gomp_InitSegmentResidueAtomList(&List,Which,listType,1) )
        return(1);

    /* Loop over sets. */
    for( i=0; i<Sets; i++ ) {
        Segment = argv[3*i+0];
        Residue = argv[3*i+1];
        Atom    = argv[3*i+2];

        /* Loop over structures. */
        for( Wstr=List.from; Wstr<=List.to; Wstr++ ) {

            if( gomp_PreSegmentResidueAtomListStructure(&List,Wstr) )
                continue;

            slong = gomp_MakeSelectionList(Wstr ,
                                         Segment , Residue, Atom, List.sel_list);

            (void)gomp_PostSegmentResidueAtomListStructure(
                &List,Wstr,slong,List.sel_list);
        }
    }

    return(gomp_ReturnAndFreeSegmentResidueAtomList(&List));
}

/************************************************************************/
int gomp_ParsePlotAxisList(const char *Segment , 
                         const char *Residue , 
                         const char *Atom    ,
                         const char *Xaxis   ,
                         const char *Yaxis   ,
                         const char *Zaxis)
/************************************************************************/
{
    register int  i;
    static int    ihelp,atom_list;
    static int *sel_list;
    static int    slong;
    static char   chelp[BUFF_LEN];
    static int    j,k;
    static int    Total;
    static const float *sumxyz;

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(!atom_list)
        /* nothing to do */
        return(0);

    if(InitAtomHitList())
        return(1);

    CoordAxisDataIsChanging();

    if(*Xaxis != (char)NULL)
        CoordAxisList.Xaxis = atof(Xaxis);
    else
        CoordAxisList.Xaxis = 2.0;
    if(*Yaxis != (char)NULL)
        CoordAxisList.Yaxis = atof(Yaxis);
    else
        CoordAxisList.Yaxis = 2.0;
    if(*Zaxis != (char)NULL)
        CoordAxisList.Zaxis = atof(Zaxis);
    else
        CoordAxisList.Zaxis = 2.0;

    sel_list  = gomp_AllocateIntVector(atom_list);

/* if there is an old list, delete it! */
/*
  if(CoordAxisList.Structures || CoordAxisList.CoordPoints) {
  (void)gomp_DeletePlotAxisList();
  }
*/
/* check to see see if in comes x,y and z coordinates */
    if(gomp_IsStringAFloat(Segment) && gomp_IsStringAFloat(Residue) && gomp_IsStringAFloat(Atom)) {

        sumxyz   = gomp_GetTranslateArray();

        if(!CoordAxisList.CoordPoints) {

            CoordAxisList.Xc          = gomp_AllocateFloatVector(1);
            CoordAxisList.Yc          = gomp_AllocateFloatVector(1);
            CoordAxisList.Zc          = gomp_AllocateFloatVector(1);

            CoordAxisList.CoordPoints = 1;
            CoordAxisList.Xc[0]       = atof(Segment) - sumxyz[0];
            CoordAxisList.Yc[0]       = atof(Residue) - sumxyz[1];
            CoordAxisList.Zc[0]       = atof(Atom)    - sumxyz[2];

        } else {
            CoordAxisList.Xc          = gomp_ReallocateFloatVector(
                CoordAxisList.Xc , 
                CoordAxisList.CoordPoints + 1);
            CoordAxisList.Yc          = gomp_ReallocateFloatVector(
                CoordAxisList.Yc ,
                CoordAxisList.CoordPoints + 1);
            CoordAxisList.Zc          = gomp_ReallocateFloatVector(
                CoordAxisList.Zc ,
                CoordAxisList.CoordPoints + 1);

            CoordAxisList.Xc[CoordAxisList.CoordPoints]  =
                atof(Segment) - sumxyz[0];
            CoordAxisList.Yc[CoordAxisList.CoordPoints]  =
                atof(Residue) - sumxyz[1];
            CoordAxisList.Zc[CoordAxisList.CoordPoints]  =
                atof(Atom)    - sumxyz[2];

            CoordAxisList.CoordPoints++;

        }

        return gomp_SetPlotterRegistrationState(
            1, &CoordAxisList.Plotter, gomp_PlotCoordAxis, NULL,
            PLOTTER_NAME_COORD_AXIS,PLOTTER_ORDER_COORD_AXIS);
    }

    Total = 0;

    for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {

        if(!gomp_GetSelectedStructure(j)) continue;

        slong     = gomp_MakeSelectionList(
            j , Segment , Residue, Atom, sel_list);
        ihelp     = 1;

        if(slong < 0) continue;

        if(slong > 0) {      /* #1 */

            if(!CoordAxisList.Structures) {
                CoordAxisList.Structure = gomp_AllocateIntVector(1);
                CoordAxisList.Structure[CoordAxisList.Structures] = j;
                CoordAxisList.Length    = gomp_AllocateIntVector(1);
                CoordAxisList.Length[CoordAxisList.Structures] = slong;
                CoordAxisList.List      = gomp_AllocateIntVector(slong);
                for(i = 0 ; i < slong ; i++)
                    CoordAxisList.List[i] = sel_list[i];
                Total += slong;
                CoordAxisList.Structures++;
            } else {

                Total = 0;
                for(k = 0 ; k < CoordAxisList.Structures ; k++) {
                    Total += CoordAxisList.Length[k]; 
                }
                CoordAxisList.Structure =
                    gomp_ReallocateIntVector(CoordAxisList.Structure , 
                                     CoordAxisList.Structures + 1);
                CoordAxisList.Structure[CoordAxisList.Structures] = j;
                CoordAxisList.Length =
                    gomp_ReallocateIntVector(CoordAxisList.Length , 
                                     CoordAxisList.Structures + 1);
                CoordAxisList.Length[CoordAxisList.Structures] = slong;
                CoordAxisList.List =
                    gomp_ReallocateIntVector(CoordAxisList.List ,Total + slong);
                for(i = 0 ; i < slong ; i++)
                    CoordAxisList.List[Total + i] = sel_list[i];
                Total += slong;
                CoordAxisList.Structures++;
            }

            gomp_SetPlotterRegistrationState(
                1, &CoordAxisList.Plotter, gomp_PlotCoordAxis, NULL,
                PLOTTER_NAME_COORD_AXIS, PLOTTER_ORDER_COORD_AXIS);

            if(SaveAtomHitList(slong , sel_list))
                gomp_PrintERROR("can't save selection list");
        }                   /* #1 */
        else {
            sprintf(chelp,"structure (%d) no atoms in the selection list",
                    (j + 1));
            gomp_PrintMessage(chelp);
        }
    }

    if(sel_list) {
        free(sel_list);
        sel_list = NULL;
    }

    return(0);
}

/************************************************************************/
int gomp_GetPlotAxisListLength()
/************************************************************************/
{
    return(CoordAxisList.Structures);
}
/************************************************************************/
int gomp_DeletePlotAxisList()
/************************************************************************/
{
    CoordAxisDataIsChanging();

    if(CoordAxisList.Structures) {
        free(CoordAxisList.Structure);
        free(CoordAxisList.Length);
        free(CoordAxisList.List);
        CoordAxisList.Structures = 0;
    }

    if(CoordAxisList.CoordPoints) {
        free(CoordAxisList.Xc);
        free(CoordAxisList.Yc);
        free(CoordAxisList.Zc);
        CoordAxisList.CoordPoints = 0;
    }

    gomp_UnregisterPlotter(
        CoordAxisList.Plotter.plotter);
    CoordAxisList.Plotter.plotter = NULL;

    return(0);
}

/************************************************************************/
const int *gomp_GetPlotAxisListP()
/************************************************************************/
{
    return(CoordAxisList.List);
}
/************************************************************************/
const int *gomp_GetPlotAxisLengthP()
/************************************************************************/
{
    return(CoordAxisList.Length);
}
/************************************************************************/
const int *gomp_GetPlotAxisStructureP()
/************************************************************************/
{
    return(CoordAxisList.Structure);
}

/************************************************************************/
int  gomp_GetPlotAxisLength(float *Xaxis , float *Yaxis , float *Zaxis)
/************************************************************************/
{
    *Xaxis = CoordAxisList.Xaxis;
    *Yaxis = CoordAxisList.Yaxis;
    *Zaxis = CoordAxisList.Zaxis;

    return(0);
}

/************************************************************************/
int   gomp_GetPlotAxisCoordPoints()
/************************************************************************/
{
    return(CoordAxisList.CoordPoints);
}
/************************************************************************/
const float *gomp_GetPlotAxisXCoord()
/************************************************************************/
{
    return(CoordAxisList.Xc);
}
/************************************************************************/
const float *gomp_GetPlotAxisYCoord()
/************************************************************************/
{
    return(CoordAxisList.Yc);
}

/************************************************************************/
const float *gomp_GetPlotAxisZCoord()
/************************************************************************/
{
    return(CoordAxisList.Zc);
}

/************************************************************************/
int  InitAtomHitList()
/************************************************************************/
{
    int   IITemp;
    const char *ITemp;

    ITemp  = Tcl_GetVar(gomp_GetTclInterp() , "gomAtomHitListActive",  TCL_GLOBAL_ONLY);

    AtomHitListActive = 1;

    if(!ITemp) {
        AtomHitListActive = 0;
    } else {
        if(!atoi(ITemp))
            AtomHitListActive = 0;
    }

    if(AtomHitListTextLength) {
        AtomHitListTextLength = 0;
        free(AtomHitList);

        IITemp = Tcl_UnsetVar(gomp_GetTclInterp(),"gomAtomHitList",TCL_GLOBAL_ONLY);
    }

    return(0);
}
/************************************************************************/
int SaveAtomHitList(int Entries , const int *EntryPointer)
/************************************************************************/
{
    int   loop;
    char *String;
    const char *ITemp;
    char  Temp[BUFF_LEN];

    if(!AtomHitListActive) {
        return(0);
    }

    if(Entries < 1) {
        sprintf(Temp,"{} ");
        String = gomp_AllocateCharVector(strlen(Temp) + 1);
        strncpy(String,Temp,strlen(Temp));
        String[strlen(Temp)] = (char)NULL;
        AtomHitListTextLength = strlen(String);
        ITemp = Tcl_SetVar(gomp_GetTclInterp(),"gomAtomHitList",String,TCL_GLOBAL_ONLY);
        if(!ITemp) {
            gomp_PrintERROR("can't set tcl variable 'gomAtomHitList'");
            return(1);
        }
        return(0);
    }

    if(Entries == 1) {
        sprintf(Temp,"{%d} ",EntryPointer[0]+1);
    } else {
        sprintf(Temp,"{%d ",EntryPointer[0]+1);
    }
    String = gomp_AllocateCharVector(strlen(Temp) + 1);
    strncpy(String,Temp,strlen(Temp));
    String[strlen(Temp)] = (char)NULL;

    for (loop = 1 ; loop < Entries ; loop++) {

        if(loop != (Entries - 1))
            sprintf(Temp,"%d ", EntryPointer[loop]+1);
        else
            sprintf(Temp,"%d} ",EntryPointer[loop]+1);

        String = gomp_ReallocateCharVector(String , 
                                        strlen(Temp) + strlen(String) + 1);
        strncat(String,Temp,strlen(Temp));
        String[strlen(String)] = (char)NULL;
    }

    if(AtomHitListTextLength) {
        ITemp = Tcl_SetVar(gomp_GetTclInterp(),"gomAtomHitList",String,TCL_APPEND_VALUE&TCL_GLOBAL_ONLY);
        if(!ITemp) {
            gomp_PrintERROR("can't append to tcl variable 'gomAtomHitList'");
            return(1);
        }
    } else {
        ITemp = Tcl_SetVar(gomp_GetTclInterp(),"gomAtomHitList",String,TCL_GLOBAL_ONLY);
        if(!ITemp) {
            gomp_PrintERROR("can't set tcl variable 'gomAtomHitList'");
            return(1);
        }
    }
    AtomHitListTextLength = strlen(String);
    free(String);

    return(0);
}

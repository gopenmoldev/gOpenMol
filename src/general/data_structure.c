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

#include <assert.h>
#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <tcl.h>

#ifndef WIN32
#include <sys/types.h>
#endif

#include "memalloc.h"
#include "axis.h"
#include "bond.h"
#include "colouring.h"
#include "contour.h"
#include "coord_man.h"
#include "gomlistener.h"
#include "gommain.h"
#include "gommonitor.h"
#include "gomstring.h"
#include "gomtext.h"
#include "label.h"
#include "model_file.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "msd.h"
#include "parser.h"
#include "picking.h"
#include "plot_molec.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "rdf.h"
#include "rforce.h"
#include "selection.h"
#include "tclutils.h"
#include "trajectory.h"

#include "stdafx.h"

#define MAX_ATM_CONN 20 /* max atom connections per atom
                           can be redined with the comand:
                           define atom maxcon Ivalue
                        */

/* define the Molecule structure                       */
typedef struct {
    char  *display_state;  /* atoms display state (ON/OFF) */
/* CPK parameters         */
    char  *CPK_state;      /* CPK display  state (ON/OFF)  */
    float *CPK_scale;      /* scale factor for scaling the cpk display */
/* licorice parameters    */
    char  *lico_state;     /* licorice display state       */
    float  lico_rads;      /* radius of liquorice sphere   */
    float  lico_radc;      /* radius of liquorice cylinder */
/* label properties       */
    char  *label_state;    /* label display state (ON/OFF) */
/* colour definitions     */
    float *red;            /* red colour                   */
    float *green;          /* green colour                 */
    float *blue;           /* blue  colour                 */
} AtomDisplayProperties;


typedef struct {
    int   *res1;         /* residue number 1 list */
    int   *res2;         /* residue number 2 list */
    atom_name_t    *x_atm; /* atom name list        */
    residue_name_t *x_res; /* residue name list     */
    segment_name_t *x_seg; /* segment name list     */
} AtomNames;

/* atom location and name info structures  */
typedef struct {
    float *x;                /* x,y and z coordinates             */
    float *y;
    float *z;
} AtomCoordinates;

typedef char GBasisSetTag_t[BUFF_LEN];

typedef struct {     /* Molec structure                   */
    float *charge;                /* atom partial charge */
    float *ncharge;               /* atom nuclear charge */
    float *covar;                 /* covalent radius     */
    float *bvalue;                /* bvalue list         */
    GBasisSetTag_t *x_GBasisSetTag; /* Basis set tag       */
    AtomData *AtomLevelInfo;      /* structure with atom type information */
} AtomParametres;

typedef struct {
    int NSize;
    int **Conn;              /* atom connection matrix            */
    int **Hbond;             /* hydrogen bond matrix              */
} AtomConnections;

typedef enum {
    AtomDataChangedListener,
    AtomCoordinateChangedListener,
    AtomDataChangedListenerCount
} AtomDataChangedListenerType;

typedef struct {
    int NAtoms;
    char                 *name;     /* structure name                 */
    char                 *fileName; /* structure file name            */
    AtomNames             names;
    AtomParametres        params;
    AtomCoordinates       coords;
    AtomConnections       connections;
    AtomDisplayProperties disp_states;
} Molecule;

static Molecule *MolecStructs = NULL;
static const Molecule NullMolecule = { 0 }; /* Zero and NULL init */

static void FreeMolecStruct(void *, size_t, const DataVectorHandler *);

static DataVectorHandler MolecStructsHandler = {
    gomp_DataVectorCopyInitFunc, /* creator function        */
    FreeMolecStruct,             /* destructor function     */
    &NullMolecule
};

static gom_PlotterData AtomPlotter = { NULL, NULL };

/* structure to save atom coodinates "push atom coordinates" */
static struct {
    int NStruct;
    struct {
        int NAtoms;
        AtomCoordinates coords;
    } *molec;
} PushedAtomCoordinates = {0,NULL};

/* General gOpenMol control structure */

static struct {
    int  Modified;              /* controls if the gOpenMol structure
                                   has changed ( == needs saving) */
    char Tag[BUFF_LEN];         /* general session TAG            */
    char Description[BUFF_LEN]; /* general session Description    */
    char Available[BUFF_LEN];   /* general session Available      */
} gOpenMol = { 1 , "" , "" , "" };

#if 0
/*  debug structure */
static struct {
    char Location[BUFF_LEN];
    int  Level;
} Debug;
#endif

static struct {
    int NAtoms;
    AtomCoordinates coords;
} SavedAtomCoords = {0, {NULL, NULL, NULL}};

/* ................ */

static int MaxAtmConnections  = MAX_ATM_CONN;

#if 0
static int SetMolecStructOrder(int , int , int);
#endif

/* end of Molec structure                           */

static int PlotMolecule(void*,int,int);

/* quality factor for the sphere */
static struct {
    int Value;
} SphereQuality = {10};

/* quality factor for the cylinder */
static struct {
    int Value;
} CylinderQuality = {10};

static void CopyWord( char *dst, const char *src, size_t len )
{
    char format[13];
    sprintf( format, "%%%zds", len - 1 );
    sscanf( src, format, dst );
}
#define CopyWord(dst,src) CopyWord(dst,src,sizeof(dst))
static int DuplicateString( char **dst, const char *src )
{
    char *tmp = gomp_AllocateCharVector( strlen(src) + 1 );
    if ( ! tmp )
        return(1);
    strcpy(tmp,src);
    gomp_FreeVector(*dst);
    *dst = tmp;
    return(0);
}

/***************************************************************************/
/* Define some macros to reduce the work needed to implement               */
/* atom parametre functions.                                               */
/***************************************************************************/

/***************************************************************************/
/* Implement following functions:                                          */
/*     const char *gomp_GetAtom<title>(int Wstr, int Point);               */
/*     int         gomp_PutAtom<title>(int Wstr, const char *Value,        */
/*                                     int Point);                         */
/*     const char *gomp_GetAtom<title>Pointer(int Wstr);                   */
/*           char *gomp_GetModifiableAtom<title>Pointer(int Wstr);         */
/* gomp_PutAtom<title> and gomp_GetModifiableAtom<title>Pointer functions  */
/* will execute the <invalidate> block in addition to the normal function. */
/***************************************************************************/
#define IMPLEMENT_GETPUT_ATOM_STRING(type,title,var,invalidate)  \
const char *gomp_GetAtom##title(int Wstr, int Point)             \
{                                                                \
    return(MolecStructs[Wstr].var[Point]);                       \
}                                                                \
int gomp_PutAtom##title(int Wstr, const char *String, int Point) \
{                                                                \
    invalidate;                                                  \
    String += strspn(String," \t\n");                            \
    CopyWord(MolecStructs[Wstr].var[Point],String);              \
    return(0);                                                   \
}                                                                \
IMPLEMENT_GET_MODIFIABLE_POINTER(                                \
    type,Atom##title##Pointer,(int Wstr),           \
    MolecStructs[Wstr].var,invalidate)

/***************************************************************************/
/* Implement following functions:                                          */
/*     type        gomp_GetAtom<title>(int Wstr, int Point);               */
/*     int         gomp_PutAtom<title>(int Wstr, type Value, int Point);   */
/*     const type *gomp_GetAtom<title>Pointer(int Wstr);                   */
/*           type *gomp_GetModifiableAtom<title>Pointer(int Wstr);         */
/* gomp_PutAtom<title> and gomp_GetModifiableAtom<title>Pointer functions  */
/* will execute the <invalidate> block in addition to the normal function. */
/***************************************************************************/
#define IMPLEMENT_GETPUT_ATOM_VALUE(type,title,var,invalidate) \
type gomp_GetAtom##title(int Wstr, int Point)                  \
{                                                              \
    return(MolecStructs[Wstr].var[Point]);                     \
}                                                              \
int  gomp_PutAtom##title(int Wstr, type Value, int Point)      \
{                                                              \
    invalidate;                                                \
    MolecStructs[Wstr].var[Point] = Value;                     \
    return(0);                                                 \
}                                                              \
IMPLEMENT_GET_MODIFIABLE_POINTER(                              \
    type,Atom##title##Pointer,(int Wstr),                      \
    MolecStructs[Wstr].var,invalidate)

/***************************************************************************/
/* Implement following functions:                                          */
/*     type        gomp_GetAtom<title>(int Wstr, int Point);               */
/*     int         gomp_PutAtom<title>(int Wstr, type Value, int Point);   */
/* gomp_PutAtom<title> function will execute the <invalidate> block        */
/* in addition to the normal function.                                     */
/***************************************************************************/
#define IMPLEMENT_GETPUT_ATOM_PARAM(type,title,var,invalidate)      \
type gomp_GetAtom##title(int Wstr, int Point)                       \
{                                                                   \
    return(MolecStructs[Wstr].params.AtomLevelInfo[Point].var);  \
}                                                                   \
int  gomp_PutAtom##title(int Wstr, type Value, int Point)           \
{                                                                   \
    invalidate;                                                     \
    MolecStructs[Wstr].params.AtomLevelInfo[Point].var = Value;  \
    return(0);                                                      \
}

/***************************************************************************/
/* Implement atom related notifiers.                                       */
/***************************************************************************/

/***************************************************************************/
/* Implement AtomDataChangedListener and AtomCoordinateChangedListener.        */
GOM_IMPLEMENT_GOM_LISTENER(AtomDataChanged)
GOM_IMPLEMENT_GOM_LISTENER(AtomCoordinateChanged)

static int CallAtomListener(
    void            *pStructMask,
    gom_ListenerFunc callback,
    void            *callbackData )
{
    return ((gom_AtomDataChangedListenerFunc)callback)(
        callbackData,*(unsigned long int*)pStructMask);
}

/***************************************************************************/
/* Define DataStateListener which will call AtomDataChangedListeners or    */
/* AtomCoordinateChangedListeners.                                         */
static struct {
    gom_UpdateDataListener* updateDataListener;
    unsigned long int change_mask;
} AtomChangeMasks[AtomDataChangedListenerCount];

static int CallTypedAtomListeners(AtomDataChangedListenerType Type,
                                  gom_ListenerList *list)
{
    gom_UpdateDataListener *listener    = AtomChangeMasks[Type].updateDataListener;
    unsigned long int       change_mask = AtomChangeMasks[Type].change_mask;
    int result;

    /**
     * Reinitialize the mask and the event. Listeners may change the
     * atom state so we can't do that after we have called the
     * notifiers, anymore.
     */
    memset(&AtomChangeMasks[Type],0,sizeof(AtomChangeMasks[Type]));
    /* Call atom notifiers. */
    result = gomp_CallListeners(list,CallAtomListener,&change_mask);
    if ( result != 0 )
        /* Error occured. */
        return result;
    gomp_CancelUpdateDataListener(listener);
    return(0);
}

static int CallAtomListeners(void *data) {
    int result = 0;
    
    /**
     * Inform that atom graphics is invalidated.
     * This could be part of AtomPropertyIsChanging macro. But to
     * reduced the work needed by AtomPropertyIsChanging it is changed
     * to here.
     */
    gomp_InvalidatePlotterDelayed(&AtomPlotter);

    switch ( (ptrdiff_t)data ) {
    case AtomDataChangedListener:
        /* Call atom property listeners. */
        result = CallTypedAtomListeners(
            AtomDataChangedListener,
            GOM_GET_LISTENER_LIST(AtomDataChanged));
        break;
    case AtomCoordinateChangedListener:
        /* Call AtomCoordinateChanged listeners and */
        /* AtomDataChanged listeners.               */
        result = CallTypedAtomListeners(
            AtomCoordinateChangedListener,
            GOM_GET_LISTENER_LIST(AtomCoordinateChanged));
        if ( result )
            /* An error occured. */
            return(result);
        result = CallTypedAtomListeners(
            AtomDataChangedListener,
            GOM_GET_LISTENER_LIST(AtomDataChanged));
        break;
    default:
        assert(0);
    }

    return(result);
}

/***************************************************************************/
/* Implement AtomStructureDestructorListener.                              */
/* The structure number Wstr passed to the notifier may be negative        */
/* meaning that all structures are deleted.                                */
GOM_IMPLEMENT_GOM_LISTENER(MolecStructPreDelete)
GOM_IMPLEMENT_GOM_LISTENER(MolecStructDelete)
GOM_IMPLEMENT_GOM_LISTENER(MolecStructPostDelete)

static int CallMolecStructDeleteListener(
    void            *VStructL,
    gom_ListenerFunc callback,
    void            *callbackData )
{
    int* StructL = VStructL;
    return ((gom_MolecStructDeleteListenerFunc)callback)(
        callbackData,StructL[0],StructL[1]);
}
static int CallMolecStructDeleteListeners3(
    int Wstr, int Dstr, gom_ListenerList *list)
{
    int StructL[2];
    StructL[0] = Wstr;
    StructL[1] = Dstr;
    return gomp_CallListeners(list,CallMolecStructDeleteListener,StructL);
}
#define CallMolecStructPreDeleteListeners(Wstr,Dstr) \
    CallMolecStructDeleteListeners3( \
        Wstr,Dstr,GOM_GET_LISTENER_LIST(MolecStructPreDelete))
#define CallMolecStructDeleteListeners(Wstr,Dstr) \
    CallMolecStructDeleteListeners3( \
        Wstr,Dstr,GOM_GET_LISTENER_LIST(MolecStructDelete))
#define CallMolecStructPostDeleteListeners(Wstr,Dstr) \
    CallMolecStructDeleteListeners3( \
        Wstr,Dstr,GOM_GET_LISTENER_LIST(MolecStructPostDelete))

/***************************************************************************/
/* Mark the atom data changed. Note that atom graphics is invalidated by   */
/* CallAtomListeners listener.                                             */
#define AtomPropertyIsChanging(Wstr,type)                              \
{                                                                      \
    AtomChangeMasks[type].change_mask |= 1<<Wstr;                      \
    if ( ! AtomChangeMasks[type].updateDataListener )                  \
        AtomChangeMasks[type].updateDataListener =                     \
            gomp_AddUpdateDataListener(CallAtomListeners,(void*)type); \
}
/***************************************************************************/


/*
  Reserve space for number of 'Atoms' atoms

  If there is already an atom structure reserved reserve a new structure
  and start filling that

  Leif Laaksonen 1994
  Eero Häkkinen  2003, 2005

*/

/***************************************************************************/
static int ReallocateMolecStruct(
    Molecule *molec, const char *FileName, int NAtoms)
/***************************************************************************/
{
    int ConnListLength, i, NTotAtoms;
    Tcl_Obj *var;
    double dvalue;
    void    *tmp;

    NTotAtoms = molec->NAtoms + NAtoms;
    
    /* structure file name */
    if ( FileName ) {
        char Name[BUFF_LEN];

        gomp_SplitFile(Name, FileName);

        if ( DuplicateString(&molec->name    ,Name    ) ||
             DuplicateString(&molec->fileName,FileName) )
            return(1);
    }
#define REALLOCATE(x) \
    ( ! ( ( tmp      = gomp_ReallocateVoidVector( \
                molec->x,NTotAtoms*sizeof(*molec->x)) ) && \
          ( molec->x = tmp ) ) )
    /* names */
    if ( REALLOCATE(names.res1) ||
         REALLOCATE(names.res2) ||
         REALLOCATE(names.x_atm ) ||
         REALLOCATE(names.x_res ) ||
         REALLOCATE(names.x_seg ) )
        return(1);
    /* properties  */
    if ( REALLOCATE(params.charge       ) ||
         REALLOCATE(params.ncharge      ) ||
         REALLOCATE(params.covar        ) ||
         REALLOCATE(params.bvalue       ) ||
         REALLOCATE(params.x_GBasisSetTag ) ||
         REALLOCATE(params.AtomLevelInfo) )
        return(1);
    /* coordinates */
    if ( REALLOCATE(coords.x) ||
         REALLOCATE(coords.y) ||
         REALLOCATE(coords.z) )
        return(1);
    /* connections */
    if ( REALLOCATE(connections.Conn ) ||
         REALLOCATE(connections.Hbond) )
        return(1);
    ConnListLength = gomp_GetMaxAtomConnections() + 1;
    while ( molec->connections.NSize < NTotAtoms ) {
        i = molec->connections.NSize;
        molec->connections.Conn[i]  = gomp_AllocateIntVector(ConnListLength);
        molec->connections.Hbond[i] = gomp_AllocateIntVector(ConnListLength);
        if ( ! molec->connections.Conn[i] ||
             ! molec->connections.Hbond[i] ) {
            gomp_FreeVector(molec->connections.Conn[i] );
            gomp_FreeVector(molec->connections.Hbond[i]);
            return(1);
        }
        ++molec->connections.NSize;
    }
    /* display properties */
    if ( REALLOCATE(disp_states.display_state) ||
         REALLOCATE(disp_states.CPK_state    ) ||
         REALLOCATE(disp_states.CPK_scale    ) ||
         REALLOCATE(disp_states.lico_state   ) ||
         REALLOCATE(disp_states.label_state  ) ||
         REALLOCATE(disp_states.red          ) ||
         REALLOCATE(disp_states.green        ) ||
         REALLOCATE(disp_states.blue         ) )
        return(1);
#undef REALLOCATE

    /* Init new atoms. */
    if ( molec->NAtoms == 0 ) {
        /* take the initial value from the tcl variable */
        var = Tcl_GetVar2Ex(
            gomp_GetTclInterp(), "gomLicoSphere", NULL, TCL_GLOBAL_ONLY);
        if ( ! var || Tcl_GetDoubleFromObj(NULL, var, &dvalue) != TCL_OK )
            dvalue = 0.3;
        molec->disp_states.lico_rads = (float)dvalue;

        var = Tcl_GetVar2Ex(
            gomp_GetTclInterp(), "gomLicoCylinder", NULL, TCL_GLOBAL_ONLY);
        if ( ! var || Tcl_GetDoubleFromObj(NULL, var, &dvalue) != TCL_OK )
            dvalue = 0.3;
        molec->disp_states.lico_radc = (float)dvalue;
    }

    for ( i = molec->NAtoms ; i < NTotAtoms ; i++ ) {
        /* Init Basis set tag. */
        memset(
            molec->params.x_GBasisSetTag[i],
            0,
            sizeof(molec->params.x_GBasisSetTag[i]) );

        /* Init connections. */
        molec->connections.Conn[i][0]  = 0;
        molec->connections.Hbond[i][0] = 0;

        /* Init display properties. */
        molec->disp_states.display_state[i] = 1;
        molec->disp_states.CPK_state[i]     = 0;
        molec->disp_states.CPK_scale[i]     = 1.0f;
        molec->disp_states.lico_state[i]    = 0;
        molec->disp_states.label_state[i]   = 0;
    }

    /* Save the new size.
     */
    molec->NAtoms = NTotAtoms;

    return(0);
}

/***************************************************************************/
void FreeMolecStruct(
    void *data, size_t size, const DataVectorHandler *pHandler)
/***************************************************************************/
{
    Molecule *molec = data;
    int i;

    /* Free arrays. */
    /* names */
    gomp_FreeVector(molec->name);
    gomp_FreeVector(molec->names.res1);
    gomp_FreeVector(molec->names.res2);
    gomp_FreeVector(molec->names.x_atm);
    gomp_FreeVector(molec->names.x_res);
    gomp_FreeVector(molec->names.x_seg);
    /* parameters  */
    gomp_FreeVector(molec->params.charge);
    gomp_FreeVector(molec->params.ncharge);
    gomp_FreeVector(molec->params.covar);
    gomp_FreeVector(molec->params.bvalue);
    gomp_FreeVector(molec->params.x_GBasisSetTag);
    gomp_FreeVector(molec->params.AtomLevelInfo);
    /* coordinates */
    gomp_FreeVector(molec->coords.x);
    gomp_FreeVector(molec->coords.y);
    gomp_FreeVector(molec->coords.z);
    /* connections */
    while ( molec->connections.NSize > 0 ) {
        i = --molec->connections.NSize;
        gomp_FreeVector(molec->connections.Conn[i] );
        gomp_FreeVector(molec->connections.Hbond[i]);
    }
    gomp_FreeVector(molec->connections.Conn);
    gomp_FreeVector(molec->connections.Hbond);
    /* display states */
    gomp_FreeVector(molec->disp_states.display_state);
    gomp_FreeVector(molec->disp_states.CPK_state);
    gomp_FreeVector(molec->disp_states.CPK_scale);
    gomp_FreeVector(molec->disp_states.lico_state);
    gomp_FreeVector(molec->disp_states.label_state);
    gomp_FreeVector(molec->disp_states.red);
    gomp_FreeVector(molec->disp_states.green);
    gomp_FreeVector(molec->disp_states.blue);
}

/***************************************************************************/
static int GetSpaceForMolecStruct(const char *Name , int NAtoms)
/***************************************************************************/
{
    if ( NAtoms < 1 ) {
        gomp_PrintERROR("number of atoms has to be > 0");
        return(1);
    }

    /* There is usually need for just one structure so there is no
     * point to allocate memory in advance. Let's use Exact
     * version.
     */
    if ( ! gomp_DataVectorCreateOrAppendExact(
             &MolecStructs,&MolecStructsHandler,1) )
        /* Allocatation failed. Error is printed. */
        return(1);

    /* Will be zero and NULL initialized by gomp_DataVectorCreateOrAppendExact
     */
    if ( ReallocateMolecStruct(
             &MolecStructs[gomp_DataVectorGetSize(&MolecStructs) - 1],
             Name, NAtoms ) ||
         gomp_AppendSpaceModelViewMatrixMT(1) ) {
        /* Allocatation failed. Error is printed. */
        gomp_DataVectorSetSize(
            &MolecStructs,
            gomp_DataVectorGetSize(&MolecStructs) - 1 );
        return(1);
    }
    
    return(0);
}
/***********************************************************************/
int gomp_CreateMolecStruct(const char *Name , int Atoms , int Type)
/***********************************************************************/
{
    switch(Type) {

    case NEW:     /* trash everything else and start new structure */
        if(gomp_GetNumMolecStructs())
            /* there are old structures */
            (void)gomp_DeleteMolecStructs();
        /* continue */

    case APPEND: /* append to old data structure */
        if(GetSpaceForMolecStruct(Name , Atoms)) {
            gomp_PrintERROR("Can't get space for atoms");
            return(-1);
        }
        (void)gomp_UpdateMolecStructList();
        break;

    default:
        gomp_PrintEXIT("ERROR type for new structure");

    }

    if( ! AtomPlotter.plotter )
        AtomPlotter.plotter = gomp_RegisterPlotter(
            PlotMolecule,NULL,PLOTTER_NAME_MOLECULE,PLOTTER_ORDER_MOLECULE);

    return(gomp_GetNumMolecStructs() - 1);
}
/***********************************************************************/
static void DeleteMolecStructRelatedData( int Wstr )
/***********************************************************************/
{
/* delete selection lists */
    gomp_DeleteSelectionList();
    gomp_DeletePickedAtomStack();
/* delete hydrogen bond search lists */
    gomp_DeleteAllHbondSearchData();
/* delete contour related information */
    gomp_DeleteAllContours();
    switch ( Wstr ) {
    case -1: /* all           */
    case  0: /* the first one */
/* delete trajectory informatiopn     */
        gomp_DeleteTrajectoryStructure();
/* delete monitor and time series */
        gomp_ResetMonitorData();
        gomp_ResetMonitorDistanceData();
        gomp_ResetMonitorAngleData();
        gomp_ResetMonitorTorsionData();
        gomp_DeleteDistanceSeries();
        gomp_DeleteAngleSeries();
        gomp_DeleteTorsionSeries();
/* delete rms fluctuation data */
        gomp_DeleteMSFset();
/* delete rdf data */
        gomp_DeleteRDF();
/* delete vector structure */
        gomp_DeleteVectorStructure();
/* delete cut plane information */
        gomp_DeleteCutPlaneDataX();
        gomp_DeleteCutPlaneDataY();
        gomp_DeleteCutPlaneDataZ();
/* reset the axis display lists        */
        gomp_DeletePlotAxisList();
        break;
    }
    if ( gomp_GetNumMolecStructs() == 0 ) {
        /* reset */
        static const float Identity[16] = {
            1.0 , 0.0 , 0.0 , 0.0,
            0.0 , 1.0 , 0.0 , 0.0,
            0.0 , 0.0 , 1.0 , 0.0,
            0.0 , 0.0 , 0.0 , 1.0 };
        gomp_UnregisterPlotter(AtomPlotter.plotter);
        AtomPlotter.plotter = NULL;
/* delete translation array */
        gomp_SaveTranslateArray(0.0 , 0.0 , 0.0);
        gomp_SaveModelViewMatrix(Identity);
        gomp_RemoveSpaceModelViewMatrixMT(0,~(size_t)0);
/* reset the need to save the gom file */
        gomp_SetFileSavingState(0);
/* delete all text stacks */
        gomp_DeleteTextStack();
/* delete the atom, residue and segment label stacks */
        gomp_ResetNameStack();
/* At the end of the reset procedure run a Tcl/Tk script 

The script is in the data/utility.tcl script and it is called:

lulUtility::gOpenMolReset

*/
        gomp_RungOpenMolResetScript();
    }
    else
        gomp_RemoveSpaceModelViewMatrixMT(Wstr,1);
/* update structure list in the widget */
    gomp_UpdateMolecStructList();
}
/***********************************************************************/
int gomp_MergeMolecStruct( int Wstr, int Dstr )
/***********************************************************************/
{
    int i, j, NWAtoms, NDAtoms, NTotAtoms;

    if ( Wstr == Dstr )
        return(0);
    
    NWAtoms   = MolecStructs[Wstr].NAtoms;
    NDAtoms   = MolecStructs[Dstr].NAtoms;
    NTotAtoms = NWAtoms + NDAtoms;
   
    /* Pre merge. */
    if ( CallMolecStructPreDeleteListeners(Wstr,Dstr) )
        return(1);
    
    if ( ReallocateMolecStruct(
            &MolecStructs[Dstr],NULL,NWAtoms) )
        return(1);
    MolecStructs[Dstr].NAtoms = NDAtoms;
        
    /* Merge. */
    CallMolecStructDeleteListeners(Wstr,Dstr);

#define MOVE(x) \
    memcpy( \
        &MolecStructs[Dstr].x[NDAtoms], \
        &MolecStructs[Wstr].x[0], \
        NWAtoms * sizeof( MolecStructs[Wstr].x[0] ) ) 
    /* names */
    MOVE(names.res1);
    MOVE(names.res2);
    MOVE(names.x_atm );
    MOVE(names.x_res );
    MOVE(names.x_seg );
    /* properties  */
    MOVE(params.charge       );
    MOVE(params.ncharge      );
    MOVE(params.covar        );
    MOVE(params.bvalue       );
    MOVE(params.x_GBasisSetTag );
    MOVE(params.AtomLevelInfo);
    /* coordinates */
    MOVE(coords.x);
    MOVE(coords.y);
    MOVE(coords.z);
    /* connections */
    MOVE(connections.Conn );
    MOVE(connections.Hbond);
    /* display properties */
    MOVE(disp_states.display_state);
    MOVE(disp_states.CPK_state    );
    MOVE(disp_states.CPK_scale    );
    MOVE(disp_states.lico_state   );
    MOVE(disp_states.label_state  );
    MOVE(disp_states.red          );
    MOVE(disp_states.green        );
    MOVE(disp_states.blue         );
#undef MOVE
    MolecStructs[Dstr].NAtoms = NTotAtoms;

    /* Convert. */
    for ( i = NDAtoms ; i < NTotAtoms ; i++ ) {
        /* Convert connections. */
        for ( j = 1 ; j <= MolecStructs[Dstr].connections.Conn[i][0]  ; j++ )
            MolecStructs[Dstr].connections.Conn[i][j]  += NDAtoms;
        for ( j = 1 ; j <= MolecStructs[Dstr].connections.Hbond[i][0] ; j++ )
            MolecStructs[Dstr].connections.Hbond[i][j] += NDAtoms;
    }

    /* Delete. */
    /* connection pointers are copied. Do not free them. */
    MolecStructs[Wstr].connections.NSize = 0;
    gomp_DataVectorRemove(&MolecStructs,Wstr,1);

    /* Post merge */
    DeleteMolecStructRelatedData(Wstr);
    CallMolecStructPostDeleteListeners(Wstr,Dstr);

    return(0);
}
/***********************************************************************/
int gomp_MergeMolecStructs()
/***********************************************************************/
{
    while ( gomp_GetNumMolecStructs() > 1 ) {
        if ( gomp_MergeMolecStruct(1,0) )
            return(1);
    }
    return(0);
}
/***************************************************************************/
int gomp_DeleteMolecStruct(int Wstr)
/***************************************************************************/
{
    /* Pre delete. */
    CallMolecStructPreDeleteListeners(Wstr,-1);

    /* Delete. */
    CallMolecStructDeleteListeners(Wstr,-1);

    if ( Wstr >= 0 )
        gomp_DataVectorRemove(&MolecStructs,Wstr,1);
    else
        gomp_DataVectorFree(&MolecStructs);

    /* Post delete. */
    DeleteMolecStructRelatedData(Wstr);
    CallMolecStructPostDeleteListeners(Wstr,-1);
    
    return(0);
}
/***************************************************************************/
int gomp_DeleteMolecStructs()
/***************************************************************************/
{
    return gomp_DeleteMolecStruct(-1);
}
/***************************************************************************/
int gomp_GetNumMolecStructs()
/***************************************************************************/
{
    return(gomp_DataVectorGetSize(&MolecStructs));
}
/***************************************************************************/
int gomp_GetNumAtomsInMolecStruct(int Structure)
/***************************************************************************/
{
    int NStruct = gomp_GetNumMolecStructs();

    if ( Structure > NStruct ) {
        gomp_PrintMessage("?Error in structure number");
        return(-1);
    }

/* if there are no structures defined there are no atoms */
    if ( NStruct < 1 )
        return(0);

    return(MolecStructs[Structure].NAtoms);
}
/*
  This routine builds the internal order of the atoms in the
  structure.

  The index means the next atom in the array. A negative index
  means that it is the last atom in the structure.

  Leif Laaksonen 1994

  WARNING!!!! no check of the index is made here
*/
#if 0
/***************************************************************************/
int SetMolecStructOrder( int Structure , int Index , int Value)
/***************************************************************************/
{
    MolecStruct.MolecAtoms[Structure].order[Index] = Value;

    return(0);
}
#endif

/***************************************************************************/
const char *gomp_GetMolecStructName(int Wstr)
/***************************************************************************/
{
    return(MolecStructs[Wstr].name);
}
/***************************************************************************/
int   gomp_PutMolecStructName(int Wstr, const char *Name)
/***************************************************************************/
{
    int Ret = DuplicateString(&MolecStructs[Wstr].name,Name);
    if ( Ret == 0 ) {
        AtomPropertyIsChanging(Wstr, AtomDataChangedListener);
        gomp_UpdateMolecStructList();
    }
    return(Ret);
}
/***************************************************************************/
const char *gomp_GetMolecStructFileName(int Wstr)
/***************************************************************************/
{
    return(MolecStructs[Wstr].fileName);
}
/***************************************************************************/
int   gomp_PutMolecStructFileName(int Wstr, const char *FileName)
/***************************************************************************/
{
    char Name[BUFF_LEN];
    if ( DuplicateString(&MolecStructs[Wstr].fileName,FileName) )
        return(1);
    gomp_SplitFile(Name, FileName);
    return(gomp_PutMolecStructName(Wstr,Name));
}
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_STRING(atom_name_t,AtmName,names.x_atm,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_STRING(residue_name_t,ResName,names.x_res,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_STRING(segment_name_t,SegName,names.x_seg,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(int,ResNum1,names.res1,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(int,ResNum2,names.res2,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/

/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(float,Charge,params.charge,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(float,NucCharge,params.ncharge,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(float,Covar,params.covar,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(float,BValue,params.bvalue,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
const GBasisSetTag_t *gomp_GetAtomBasisSetTagPointer(int);
      GBasisSetTag_t *gomp_GetModifiableAtomBasisSetTagPointer(int);
IMPLEMENT_GETPUT_ATOM_STRING(GBasisSetTag_t,BasisSetTag,params.x_GBasisSetTag,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/

/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_PARAM(int,Type,type,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_PARAM(float,BndRad,bndrad,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_PARAM(float,VdwRad,vdwrad,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_PARAM(float,PluRad,plurad,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_PARAM(char,Global,global,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_PARAM(float,Emin,emin,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_PARAM(float,Rmin,rmin,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_PARAM(float,Patom,patom,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_PARAM(float,Mass,mass,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_PARAM(int,Cnct,cnct,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_PARAM(char,Hbond,hbond,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
const char *gomp_GetAtomAtype(int Wstr, int Hit)
/***********************************************************************/
{
    return(MolecStructs[Wstr].params.AtomLevelInfo[Hit].atype);
}
/***********************************************************************/
int   gomp_PutAtomAtype(int Wstr, const char *Value , int Hit)
/***********************************************************************/
{
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener);
    CopyWord(MolecStructs[Wstr].params.AtomLevelInfo[Hit].atype,Value);
    return(0);
}
/***********************************************************************/

/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(float,XCoord,coords.x,
    AtomPropertyIsChanging(Wstr, AtomCoordinateChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(float,YCoord,coords.y,
    AtomPropertyIsChanging(Wstr, AtomCoordinateChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(float,ZCoord,coords.z,
    AtomPropertyIsChanging(Wstr, AtomCoordinateChangedListener))
/***********************************************************************/

/***********************************************************************/
/* gomp_GetAtomConnection, gomp_GetModifiableAtomConnection            */
/***********************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(
    int,AtomConnection,(int Wstr, int Point),
    MolecStructs[Wstr].connections.Conn[Point],
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
/* gomp_GetAtomHydrogenBond, gomp_GetModifiableAtomHydrogenBond        */
/***********************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(
    int,AtomHydrogenBond,(int Wstr, int Point),
    MolecStructs[Wstr].connections.Hbond[Point],
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/

/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(
    char, DisplayState,     disp_states.display_state,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(
    char, CPKDisplayState,  disp_states.CPK_state,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(
    float,CPKScale,         disp_states.CPK_scale,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(
    char, LicoDisplayState, disp_states.lico_state,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
float gomp_GetAtomLicoRadS(int Wstr)
/***********************************************************************/
{
    return(MolecStructs[Wstr].disp_states.lico_rads);
}
/***********************************************************************/
int gomp_PutAtomLicoRadS(int Wstr , float Value)
/***********************************************************************/
{
    MolecStructs[Wstr].disp_states.lico_rads = Value;
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener);
    return(0);
}
/***********************************************************************/
float gomp_GetAtomLicoRadC(int Wstr)
/***********************************************************************/
{
    return(MolecStructs[Wstr].disp_states.lico_radc);
}
/***********************************************************************/
int gomp_PutAtomLicoRadC(int Wstr , float Value)
/***********************************************************************/
{
    MolecStructs[Wstr].disp_states.lico_radc = Value;
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener);
    return(0);
}
/***********************************************************************/
IMPLEMENT_GETPUT_ATOM_VALUE(
    char, LabelDisplayState, disp_states.label_state,{
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener);
    gomp_AtomLabelDataIsChanging();
})
/***********************************************************************/
int gomp_GetAtomColour(int Wstr, float *r, float *g, float *b, int Point)
/***********************************************************************/
{
    *r = MolecStructs[Wstr].disp_states.red  [Point];
    *g = MolecStructs[Wstr].disp_states.green[Point];
    *b = MolecStructs[Wstr].disp_states.blue [Point];

    return(0);
}
/***********************************************************************/
int gomp_PutAtomColour(int Wstr, float r, float g, float b, int Point)
/***********************************************************************/
{
    MolecStructs[Wstr].disp_states.red  [Point] = r;
    MolecStructs[Wstr].disp_states.green[Point] = g;
    MolecStructs[Wstr].disp_states.blue [Point] = b;

    AtomPropertyIsChanging(Wstr, AtomDataChangedListener);

    return(0);
}
/***********************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(
    float, AtomColourRedPointer, (int Wstr),
    MolecStructs[Wstr].disp_states.red,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(
    float, AtomColourGreenPointer, (int Wstr),
    MolecStructs[Wstr].disp_states.green,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(
    float, AtomColourBluePointer, (int Wstr),
    MolecStructs[Wstr].disp_states.blue,
    AtomPropertyIsChanging(Wstr, AtomDataChangedListener))
/***********************************************************************/

/***********************************************************************/
int gomp_GetMaxAtomConnections()
/***********************************************************************/
{
    return(MaxAtmConnections);
}
/***********************************************************************/
int gomp_SetMaxAtomConnections(int AtmC)
/***********************************************************************/
{
    int Temp;

    if(AtmC < 1) {
        gomp_PrintERROR("max atom connectivity has to be > 0");
        return(0);
    }
/* must delete entire gOpenMol structure to do this */
    if(gomp_GetNumMolecStructs())
        Temp = gomp_ResetgOpenMol();
    else
        Temp = 0;

    MaxAtmConnections = AtmC;

    return(Temp);
}
/***********************************************************************/
int gomp_GetTotalNumberOfAtoms()
/***********************************************************************/
{
    int i,j;
    static int Total;

    if(!(j = gomp_GetNumMolecStructs())) 
        return(0);

    Total = 0;

    for(i = 0 ; i < j ; i++)
        Total += gomp_GetNumAtomsInMolecStruct(i);

    return(Total);
}
/***********************************************************************/
int PlotMolecule(void *userData,int Wstr,int drawFlags)
/***********************************************************************/
{
    if( Wstr >= gomp_GetNumMolecStructs() ) return(-1);

    if( drawFlags & gom_PlotSimpleElements )
        gomp_PlotMoleculeStick(Wstr);
    if( drawFlags & gom_PlotComplexElements ) {
        gomp_PlotMoleculeCPK(Wstr);
        gomp_PlotMoleculeLicorice(Wstr);
    }

    return(0);
}

/***********************************************************************/
int gomp_AssignAtomProperties(int Wstr, int Hit, int SmashHit)
/***********************************************************************/
{
    float FNCharge;

    MolecStructs[Wstr].params.AtomLevelInfo[Hit] = 
        gomp_AtomTypes.AtomParams[SmashHit];

    FNCharge = (float)gomp_AtomTypes.AtomParams[SmashHit].ncharge;

    (void)gomp_PutAtomNucCharge( Wstr , FNCharge , Hit);

    AtomPropertyIsChanging(Wstr, AtomDataChangedListener);

    return(0);
}

/***********************************************************************/
int gomp_AssignAtomColourProperties(int Wstr, int Hit, const char *Colour)
/***********************************************************************/
{
    static float r,g,b;
    static int retv;

    retv = gomp_ColourName2RGB(Colour, &r, &g , &b);

    (void)gomp_PutAtomColour(Wstr , r , g , b , Hit);

    AtomPropertyIsChanging(Wstr, AtomDataChangedListener);

    return(retv);
}
/***********************************************************************/
int gomp_AssignAtomBasisSet(int Wstr, int Hit, const char *BasisSet)
/***********************************************************************/
{
    return(gomp_PutAtomBasisSetTag(Wstr , BasisSet , Hit));
}
     



/*
  This function returns the number of atoms up to structure 'Place'
*/
#if 0
/***********************************************************************/
int gomp_AtomsUptoStructure(int Place)
/***********************************************************************/
{
    static int i;
    static int t;

    t = 0;

    if(Place > gomp_GetNumMolecStructs()) Place = gomp_GetNumMolecStructs();
    if(Place < 0) Place = 0;
    if(!Place)   return(0);
    
    for(i = 0 ; i < Place ; i++) 
        t += gomp_GetNumAtomsInMolecStruct(i);        

    return(t);
}
#endif
/* ..................... */


/****************************************************************************/
const char *gomp_GetTagText()
/****************************************************************************/
{
    return(gOpenMol.Tag);
}
/****************************************************************************/
const char *gomp_GetDescText()
/****************************************************************************/
{
    return(gOpenMol.Description);
}
/****************************************************************************/
const char *gomp_GetAvailText()
/****************************************************************************/
{
    return(gOpenMol.Available);
}

/****************************************************************************/
int  gomp_PutTagText(const char *Text)
/****************************************************************************/
{
    gomp_CopyString(gOpenMol.Tag, Text, BUFF_LEN);

    return(0);
}

/****************************************************************************/
int   gomp_PutDescText(const char *Text)
/****************************************************************************/
{
    gomp_CopyString(gOpenMol.Description, Text, BUFF_LEN);

    return(0);
}
/****************************************************************************/
int   gomp_PutAvailText(const char *Text)
/****************************************************************************/
{
    gomp_CopyString(gOpenMol.Available, Text, BUFF_LEN);

    return(0);
}
#if 0
/****************************************************************************/
const char *gomp_GetgOpenMolFileName()
/****************************************************************************/
{
    return(gOpenMol.FileName);
}
/****************************************************************************/
int   gomp_PutgOpenMolFileName(const char *FileName)
/****************************************************************************/
{
    gomp_CopyString(gOpenMol.FileName,FileName,BUFF_LEN);

    return(0);
}
#endif
/****************************************************************************/
int   gomp_gOpenMolNeedsSaving()
/****************************************************************************/
{
    return(gOpenMol.Modified);
}
/****************************************************************************/
int   gomp_SetFileSavingState(int Value)
/****************************************************************************/
{
    gOpenMol.Modified = Value;

    return(0);
}
/****************************************************************************/
int    gomp_GetMinRes1Num()
/****************************************************************************/
{
    int  i,j;
    int  MinRes;
    const int *Res;

    MinRes = 100000;
    for ( i = 0 ; i < gomp_GetNumMolecStructs() ; i++ ) {
        Res  = gomp_GetAtomResNum1Pointer(i);
        for ( j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++ ) {
            if ( Res[j] < MinRes )
                MinRes = Res[j];
        }
    }

    return(MinRes);
}
/****************************************************************************/
int gomp_GetMaxRes1Num()
/****************************************************************************/
{
    int  i,j;
    int  MaxRes;
    const int *Res;

    MaxRes = 0;
    for ( i = 0 ; i < gomp_GetNumMolecStructs() ; i++ ) {
        Res = gomp_GetAtomResNum1Pointer(i);
        for ( j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++ ) {
            if ( Res[j] > MaxRes )
                MaxRes = Res[j];
        }
    }

    return(MaxRes);
}

#if 0
/****************************************************************************/
int gomp_GetMinRes2Num()
/****************************************************************************/
{
    int  i,j;
    int  MinRes;
    const int *Res;

    MinRes = 100000;
    for ( i = 0 ; i < gomp_GetNumMolecStructs() ; i++ ) {
        Res   = gomp_GetAtomResNum2Pointer(i);
        for ( j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++ ) {
            if ( Res[j] < MinRes )
                MinRes = Res[j];
        }
    }

    return(MinRes);
}


/****************************************************************************/
int gomp_GetMaxRes2Num()
/****************************************************************************/
{
    int  i,j;
    int  MaxRes;
    const int *Res;

    MaxRes = 0;
    for ( i = 0 ; i < gomp_GetNumMolecStructs() ; i++ ) {
        Res   = gomp_GetAtomResNum2Pointer(i);
        for ( j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++ ) {
            if ( Res[j] > MaxRes )
                MaxRes = Res[j];
        }
    }

    return(MaxRes);
}
/****************************************************************************/
int gomp_SetDebugLevel(int Level)
/****************************************************************************/
{
    Debug.Level = Level;

    return(0);
}
/****************************************************************************/
int gomp_GetDebugLevel()
/****************************************************************************/
{
    return(Debug.Level);
}
#endif
/***************************************************************************/
int   gomp_SaveAtomCoords(int Wstr)
/***************************************************************************/
{
    const float *x;
    const float *y;
    const float *z;

    static int i, NAtoms;

    SavedAtomCoords.NAtoms = 0;
    free(SavedAtomCoords.coords.x);
    free(SavedAtomCoords.coords.y);
    free(SavedAtomCoords.coords.z);

    NAtoms = gomp_GetNumAtomsInMolecStruct(Wstr);
    if ( ! ( SavedAtomCoords.coords.x = gomp_AllocateFloatVector(NAtoms) ) ||
         ! ( SavedAtomCoords.coords.y = gomp_AllocateFloatVector(NAtoms) ) ||
         ! ( SavedAtomCoords.coords.z = gomp_AllocateFloatVector(NAtoms) ) )
        return(1);

    x  = gomp_GetAtomXCoordPointer(Wstr);
    y  = gomp_GetAtomYCoordPointer(Wstr);
    z  = gomp_GetAtomZCoordPointer(Wstr);

    for(i = 0; i < NAtoms ; i++) {
        SavedAtomCoords.coords.x[i] = x[i];
        SavedAtomCoords.coords.y[i] = y[i];
        SavedAtomCoords.coords.z[i] = z[i];
    }

    SavedAtomCoords.NAtoms = NAtoms;

    return(0);
}
/***************************************************************************/
int   gomp_GetSavedAtomCoords(int Wstr)
/***************************************************************************/
{
    float *x;
    float *y;
    float *z;

    static int i;

    if ( ! SavedAtomCoords.NAtoms ) 
        return(1); /* no saved set available */
    if ( SavedAtomCoords.NAtoms != gomp_GetNumAtomsInMolecStruct(Wstr) )
        return(2); /* number of saved atoms != actual number */

    x  = gomp_GetModifiableAtomXCoordPointer(Wstr);
    y  = gomp_GetModifiableAtomYCoordPointer(Wstr);
    z  = gomp_GetModifiableAtomZCoordPointer(Wstr);

    for ( i = 0; i < SavedAtomCoords.NAtoms ; i++ ) {
        x[i] = SavedAtomCoords.coords.x[i];
        y[i] = SavedAtomCoords.coords.y[i];
        z[i] = SavedAtomCoords.coords.z[i];
    }

    SavedAtomCoords.NAtoms = 0;
    free(SavedAtomCoords.coords.x);
    free(SavedAtomCoords.coords.y);
    free(SavedAtomCoords.coords.z);

    SavedAtomCoords.coords.x =
        SavedAtomCoords.coords.y =
        SavedAtomCoords.coords.z = NULL;
    
    AtomPropertyIsChanging(Wstr, AtomCoordinateChangedListener);

    return(0);
}

/***************************************************************************/
int    gomp_GetSphereQuality()
/***************************************************************************/
{
    return(SphereQuality.Value);
}

/***************************************************************************/
int    gomp_SetSphereQuality(int Value)
/***************************************************************************/
{
    SphereQuality.Value = Value;

    return(0);
} 
/***************************************************************************/
int    gomp_GetCylinderQuality()
/***************************************************************************/
{
    return(CylinderQuality.Value);
}

/***************************************************************************/
int    gomp_SetCylinderQuality(int Value)
/***************************************************************************/
{
    CylinderQuality.Value = Value;

    return(0);
}
#if 0
/***************************************************************************/
int    gomp_CopyAtomInfoInStructure(int From , int To)
/***************************************************************************/
{
    int i,j;
    int Swop;
    int ConnMax;
    int Atoms;

    if(gomp_GetNumAtomsInMolecStruct(From) != gomp_GetNumAtomsInMolecStruct(To))
        return(0);

    if(From >= To) {
        gomp_PrintERROR("from index has to < than to index");
        return(1);
    }

    if(From < 0 || To < 0) {
        gomp_PrintERROR("either index is < 1");
        return(1);
    }

    if(From >= gomp_GetNumMolecStructs() || To >= gomp_GetNumMolecStructs()) {
        gomp_PrintERROR("either index is > number of structures");
        return(1);
    }

/* do the job ... */

    ConnMax  = gomp_GetMaxAtomConnections() + 1;

    Swop     = gomp_AtomsUptoStructure(From);

    Atoms                                = MolecStruct.MolecAtoms[From].numat;
    MolecStruct.MolecAtoms[To].numat  = MolecStruct.MolecAtoms[From].numat;

    strncpy(MolecStruct.MolecAtoms[To].name ,
            MolecStruct.MolecAtoms[From].name, 
            strlen(MolecStruct.MolecAtoms[From].name));

    strncpy(MolecStruct.MolecAtoms[To].segment, 
            MolecStruct.MolecAtoms[From].segment,
            MAX_SEG_NAME_LEN * Atoms);
    strncpy(MolecStruct.MolecAtoms[To].resname,    
            MolecStruct.MolecAtoms[From].resname,
            MAX_RES_NAME_LEN * Atoms);
    strncpy(MolecStruct.MolecAtoms[To].atname,
            MolecStruct.MolecAtoms[From].atname,
            MAX_ATM_NAME_LEN * Atoms);

    strncpy(MolecStruct.MolecAtoms[To].GBasisSetTag,
            MolecStruct.MolecAtoms[From].GBasisSetTag,
            BUFF_LEN * Atoms);

    MolecStruct.DispLists[To].lico_rads = 
        MolecStruct.DispLists[From].lico_rads;
    MolecStruct.DispLists[To].lico_radc = 
        MolecStruct.DispLists[From].lico_radc;

    for( i = 0 ; i < gomp_GetNumAtomsInMolecStruct(From) ; i++) {

        for(j = 0 ; j < ConnMax ; j++) {
            MolecStruct.AtmConn[Swop + i + gomp_GetNumAtomsInMolecStruct(From)][j]      =
                MolecStruct.AtmConn[Swop + i][j];
            MolecStruct.AtmHbond[Swop + i + gomp_GetNumAtomsInMolecStruct(From)][j]     =
                MolecStruct.AtmHbond[Swop + i][j];
        }

        MolecStruct.MolecAtoms[To].order[i]          = 
            MolecStruct.MolecAtoms[From].order[i];
        MolecStruct.MolecAtoms[To].res1[i]          = 
            MolecStruct.MolecAtoms[From].res1[i];
        MolecStruct.MolecAtoms[To].res2[i]          = 
            MolecStruct.MolecAtoms[From].res2[i];
        MolecStruct.MolecAtoms[To].x[i]             = 
            MolecStruct.MolecAtoms[From].x[i];
        MolecStruct.MolecAtoms[To].y[i]             = 
            MolecStruct.MolecAtoms[From].y[i];
        MolecStruct.MolecAtoms[To].z[i]             = 
            MolecStruct.MolecAtoms[From].z[i];
        MolecStruct.MolecAtoms[To].charge[i]        = 
            MolecStruct.MolecAtoms[From].charge[i];
        MolecStruct.MolecAtoms[To].ncharge[i]       = 
            MolecStruct.MolecAtoms[From].ncharge[i];
        MolecStruct.MolecAtoms[To].covar[i]         = 
            MolecStruct.MolecAtoms[From].covar[i];
        MolecStruct.MolecAtoms[To].bvalue[i]        = 
            MolecStruct.MolecAtoms[From].bvalue[i];

        MolecStruct.MolecAtoms[To].AtomLevelInfo[i].type   =
            MolecStruct.MolecAtoms[From].AtomLevelInfo[i].type;
        MolecStruct.MolecAtoms[To].AtomLevelInfo[i].bndrad =
            MolecStruct.MolecAtoms[From].AtomLevelInfo[i].bndrad;
        MolecStruct.MolecAtoms[To].AtomLevelInfo[i].vdwrad =
            MolecStruct.MolecAtoms[From].AtomLevelInfo[i].vdwrad;
        MolecStruct.MolecAtoms[To].AtomLevelInfo[i].plurad =
            MolecStruct.MolecAtoms[From].AtomLevelInfo[i].plurad;
        MolecStruct.MolecAtoms[To].AtomLevelInfo[i].global =
            MolecStruct.MolecAtoms[From].AtomLevelInfo[i].global;
        MolecStruct.MolecAtoms[To].AtomLevelInfo[i].emin   =
            MolecStruct.MolecAtoms[From].AtomLevelInfo[i].emin;
        MolecStruct.MolecAtoms[To].AtomLevelInfo[i].rmin   =
            MolecStruct.MolecAtoms[From].AtomLevelInfo[i].rmin;
        MolecStruct.MolecAtoms[To].AtomLevelInfo[i].patom  =
            MolecStruct.MolecAtoms[From].AtomLevelInfo[i].patom;
        MolecStruct.MolecAtoms[To].AtomLevelInfo[i].mass   =
            MolecStruct.MolecAtoms[From].AtomLevelInfo[i].mass;
        MolecStruct.MolecAtoms[To].AtomLevelInfo[i].cnct   =
            MolecStruct.MolecAtoms[From].AtomLevelInfo[i].cnct;
        MolecStruct.MolecAtoms[To].AtomLevelInfo[i].hbond  =
            MolecStruct.MolecAtoms[From].AtomLevelInfo[i].hbond;
        strncpy(MolecStruct.MolecAtoms[To].AtomLevelInfo[i].atype,
                MolecStruct.MolecAtoms[From].AtomLevelInfo[i].atype,
                4);

        MolecStruct.DispLists[To].disp_list[i] = 
            MolecStruct.DispLists[From].disp_list[i];
        MolecStruct.DispLists[To].CPK_list[i]  = 
            MolecStruct.DispLists[From].CPK_list[i];
        MolecStruct.DispLists[To].cpk_scale[i] = 
            MolecStruct.DispLists[From].cpk_scale[i];
        MolecStruct.DispLists[To].lico_list[i] = 
            MolecStruct.DispLists[From].lico_list[i];

        MolecStruct.DispLists[To].red[i]       = 
            MolecStruct.DispLists[From].red[i];
        MolecStruct.DispLists[To].green[i]     = 
            MolecStruct.DispLists[From].green[i];
        MolecStruct.DispLists[To].blue[i]      = 
            MolecStruct.DispLists[From].blue[i];

        MolecStruct.DispLists[To].disp_list[i] = 
            MolecStruct.DispLists[From].disp_list[i];
        MolecStruct.DispLists[To].CPK_list[i]  = 
            MolecStruct.DispLists[From].CPK_list[i];
        MolecStruct.DispLists[To].cpk_scale[i] = 
            MolecStruct.DispLists[From].cpk_scale[i];
        MolecStruct.DispLists[To].lico_list[i] = 
            MolecStruct.DispLists[From].lico_list[i];
        MolecStruct.DispLists[To].label_toggle[i] = 
            MolecStruct.DispLists[From].label_toggle[i];
    }

    return(0);
}

/* observe that position runs from 1 ... n */
/***************************************************************************/
int gomp_DeleteMolecStructPosition(int Position)
/***************************************************************************/
{
    int i,j;
    int Atoms;
    int AtomsD;
    int Wstr;
    int Last;
    const int **TempConn;
    int MaxC;

    if(!(MolecStruct.NumMolecStruct)) { /* no previous structures */
        gomp_PrintMessage("?gOpenMol ERROR");
        gomp_PrintMessage("?No structure defined to be deleted");
        return(1);
    }

    MaxC = MAX_ATM_CONN;

/* look if the position is != last ... */
    if(Position != MolecStruct.NumMolecStruct) { /* shuffle list */

/* handle connectivity list ... */
        AtomsD = 0;
        for( i = 0 ; i < MolecStruct.NumMolecStruct ; i++) 
            AtomsD += gomp_GetNumAtomsInMolecStruct(i);

        Atoms      = AtomsD - gomp_GetNumAtomsInMolecStruct(Position - 1);
        TempConn   = 
            gomp_AllocateVoidVector(Atoms * sizeof(const int *));

        for(j = 0 ; j < Atoms ; j++) {
            TempConn[j] = gomp_AllocateIntVector(MaxC);
        }

/* handle the rest ...          */

        Wstr  = Position - 1;
        Last  = MolecStruct.NumMolecStruct - 1;

/* number of atoms in the position to delete */
        AtomsD= gomp_GetNumAtomsInMolecStruct(Wstr);
/* number of atoms in the last position      */
        Atoms = MolecStruct.MolecAtoms[MolecStruct.NumMolecStruct - 1].numat;

        MolecStruct.MolecAtoms[Wstr].numat         = Atoms;

        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].order);
        MolecStruct.MolecAtoms[Wstr].order         = 
            MolecStruct.MolecAtoms[Last].order;

        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].res1);
        MolecStruct.MolecAtoms[Wstr].res1          = 
            MolecStruct.MolecAtoms[Last].res1;

        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].res2);
        MolecStruct.MolecAtoms[Wstr].res2          = 
            MolecStruct.MolecAtoms[Last].res2;

        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].segment);
        MolecStruct.MolecAtoms[Wstr].segment       = 
            MolecStruct.MolecAtoms[Last].segment;

        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].resname);
        MolecStruct.MolecAtoms[Wstr].resname       = 
            MolecStruct.MolecAtoms[Last].resname;

        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].atname);
        MolecStruct.MolecAtoms[Wstr].atname        = 
            MolecStruct.MolecAtoms[Last].atname;

        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].x);
        MolecStruct.MolecAtoms[Wstr].x             = 
            MolecStruct.MolecAtoms[Last].x;
        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].y);
        MolecStruct.MolecAtoms[Wstr].y             = 
            MolecStruct.MolecAtoms[Last].y;
        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].z);
        MolecStruct.MolecAtoms[Wstr].z             = 
            MolecStruct.MolecAtoms[Last].z;

        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].charge);
        MolecStruct.MolecAtoms[Wstr].charge        = 
            MolecStruct.MolecAtoms[Last].charge;
        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].ncharge);
        MolecStruct.MolecAtoms[Wstr].ncharge       = 
            MolecStruct.MolecAtoms[Last].ncharge;

        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].covar);
        MolecStruct.MolecAtoms[Wstr].covar         = 
            MolecStruct.MolecAtoms[Last].covar;

        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].bvalue);
        MolecStruct.MolecAtoms[Wstr].bvalue        = 
            MolecStruct.MolecAtoms[Last].bvalue;

        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].AtomLevelInfo);
        MolecStruct.MolecAtoms[Wstr].AtomLevelInfo = 
            MolecStruct.MolecAtoms[Last].AtomLevelInfo;

        gomp_FreeVector(MolecStruct.MolecAtoms[Wstr].GBasisSetTag);
        MolecStruct.MolecAtoms[Wstr].GBasisSetTag  = 
            MolecStruct.MolecAtoms[Last].GBasisSetTag;

        gomp_FreeVector(MolecStruct.DispLists[Wstr].disp_list);
        MolecStruct.DispLists[Wstr].disp_list= 
            MolecStruct.DispLists[Last].disp_list;

        gomp_FreeVector(MolecStruct.DispLists[Wstr].CPK_list);
        MolecStruct.DispLists[Wstr].CPK_list = 
            MolecStruct.DispLists[Last].CPK_list;
        gomp_FreeVector(MolecStruct.DispLists[Wstr].cpk_scale);
        MolecStruct.DispLists[Wstr].cpk_scale= 
            MolecStruct.DispLists[Last].cpk_scale;
        gomp_FreeVector(MolecStruct.DispLists[Wstr].lico_list);
        MolecStruct.DispLists[Wstr].lico_list= 
            MolecStruct.DispLists[Last].lico_list;
        gomp_FreeVector(MolecStruct.DispLists[Wstr].label_toggle);
        MolecStruct.DispLists[Wstr].lico_list= 
            MolecStruct.DispLists[Last].label_toggle;

        gomp_FreeVector(MolecStruct.DispLists[Wstr].red);
        MolecStruct.DispLists[Wstr].red      = 
            MolecStruct.DispLists[Last].red;
        gomp_FreeVector(MolecStruct.DispLists[Wstr].green);
        MolecStruct.DispLists[Wstr].green    = 
            MolecStruct.DispLists[Last].green;
        gomp_FreeVector(MolecStruct.DispLists[Wstr].blue);
        MolecStruct.DispLists[Wstr].blue     = 
            MolecStruct.DispLists[Last].blue;

/* loop to put the data from last position into the 'position' */
        for( i = 0 ; i < Atoms ; i++) {
        }
    }
    else {  /* delete the last position */

        for(i = 0 ; i < gomp_GetTotalNumberOfAtoms() ; i++) {
            gomp_FreeVector(MolecStruct.AtmConn[i]);
            gomp_FreeVector(MolecStruct.AtmHbond[i]);
        }
        Last  = MolecStruct.NumMolecStruct - 1;

        gomp_FreeVector(MolecStruct.AtmConn);
        gomp_FreeVector(MolecStruct.AtmHbond);

        gomp_FreeVector(MolecStruct.MolecAtoms[Last].res1);
        gomp_FreeVector(MolecStruct.MolecAtoms[Last].res2);
        gomp_FreeVector(MolecStruct.MolecAtoms[Last].segment);
        gomp_FreeVector(MolecStruct.MolecAtoms[Last].resname);
        gomp_FreeVector(MolecStruct.MolecAtoms[Last].atname);
        gomp_FreeVector(MolecStruct.MolecAtoms[Last].x);
        gomp_FreeVector(MolecStruct.MolecAtoms[Last].y);
        gomp_FreeVector(MolecStruct.MolecAtoms[Last].z);
        gomp_FreeVector(MolecStruct.MolecAtoms[Last].charge);
        gomp_FreeVector(MolecStruct.MolecAtoms[Last].ncharge);
        gomp_FreeVector(MolecStruct.MolecAtoms[Last].covar);
        gomp_FreeVector(MolecStruct.MolecAtoms[Last].bvalue);

        gomp_FreeVector(MolecStruct.MolecAtoms[Last].GBasisSetTag);
        gomp_FreeVector(MolecStruct.MolecAtoms[Last].AtomLevelInfo);

        gomp_FreeVector(MolecStruct.DispLists[Last].disp_list);
        gomp_FreeVector(MolecStruct.DispLists[Last].CPK_list);
        gomp_FreeVector(MolecStruct.DispLists[Last].cpk_scale);
        gomp_FreeVector(MolecStruct.DispLists[Last].lico_list);

        gomp_FreeVector(MolecStruct.DispLists[Last].red);
        gomp_FreeVector(MolecStruct.DispLists[Last].green);
        gomp_FreeVector(MolecStruct.DispLists[Last].blue);


        gomp_FreeVector(&MolecStruct.MolecAtoms[Last]);
        gomp_FreeVector(MolecStruct.DispLists);

        MolecStruct.NumMolecStruct = 0;

        gomp_PutgOpenMolFileName("\0");
    }

    return(0);
}
#endif
/***************************************************************************/
static void FreePushedAtomCoordinates()
/***************************************************************************/
{
    int i;
    
    for ( i = 0 ; i < PushedAtomCoordinates.NStruct ; i++ ) {
        free(PushedAtomCoordinates.molec[i].coords.x);
        free(PushedAtomCoordinates.molec[i].coords.y);
        free(PushedAtomCoordinates.molec[i].coords.z);
    }
    PushedAtomCoordinates.NStruct = 0;
    free(PushedAtomCoordinates.molec);
    PushedAtomCoordinates.molec = NULL; 
}
/***************************************************************************/
int gomp_PushAtomCoordinates()
/***************************************************************************/
{
    static int i,j;
    static const float *XC,*YC,*ZC;
    static int NumStruct, NumAtoms;


/* Structure loop */
    NumStruct = gomp_GetNumMolecStructs();

    if ( NumStruct < 1 )
        return(0);

    if ( PushedAtomCoordinates.NStruct ) {
        gomp_PrintWARNING("stack is already reserved! "
                          "Old stack will be deleted");
        FreePushedAtomCoordinates();
    }

    PushedAtomCoordinates.molec = malloc(
        sizeof(*PushedAtomCoordinates.molec) * NumStruct);

    if ( PushedAtomCoordinates.molec == NULL ) {
        gomp_PrintERROR("can't reserve space for push coordinates stack");
        return(1);
    }

    PushedAtomCoordinates.NStruct = NumStruct;

    for ( i = 0 ; i < NumStruct ; i++ ) {

        NumAtoms = gomp_GetNumAtomsInMolecStruct(i);

        XC       = gomp_GetAtomXCoordPointer(i);
        YC       = gomp_GetAtomYCoordPointer(i);
        ZC       = gomp_GetAtomZCoordPointer(i);

        PushedAtomCoordinates.molec[i].NAtoms = NumAtoms;
        if ( ! ( PushedAtomCoordinates.molec[i].coords.x =
                 gomp_AllocateFloatVector(NumAtoms) ) ||
             ! ( PushedAtomCoordinates.molec[i].coords.y =
                 gomp_AllocateFloatVector(NumAtoms) ) ||
             ! ( PushedAtomCoordinates.molec[i].coords.z =
                 gomp_AllocateFloatVector(NumAtoms) ) ) {
            FreePushedAtomCoordinates();
            gomp_PrintERROR("can't reserve space for push coordinates stack");
            return(1);
        }
        
        for ( j = 0 ; j < NumAtoms ; j++ ) {
            PushedAtomCoordinates.molec[i].coords.x[j] = XC[j];
            PushedAtomCoordinates.molec[i].coords.y[j] = YC[j];
            PushedAtomCoordinates.molec[i].coords.z[j] = ZC[j];
        }
    }

    return(0);
}
/***************************************************************************/
int gomp_PopAtomCoordinates()
/***************************************************************************/
{
    static int i,j;
    static float *XC,*YC,*ZC;

    if ( ! PushedAtomCoordinates.NStruct ) {
        gomp_PrintWARNING("stack has no entries");
        return(1);
    }

    if ( PushedAtomCoordinates.NStruct != gomp_GetNumMolecStructs() ) {
        gomp_PrintWARNING("number of saved atoms != actual number");
        return(1);
    }

    for ( i = 0 ; i < PushedAtomCoordinates.NStruct ; i++ ) {

        if ( PushedAtomCoordinates.molec[i].NAtoms !=
             gomp_GetNumAtomsInMolecStruct(i) ) {
            gomp_PrintWARNING("number of saved atoms != actual number");
            continue;
        }

        XC = gomp_GetModifiableAtomXCoordPointer(i);
        YC = gomp_GetModifiableAtomYCoordPointer(i);
        ZC = gomp_GetModifiableAtomZCoordPointer(i);

        for ( j = 0 ; j < PushedAtomCoordinates.molec[i].NAtoms ; j++ ) {
            XC[j] = PushedAtomCoordinates.molec[i].coords.x[j];
            YC[j] = PushedAtomCoordinates.molec[i].coords.y[j];
            ZC[j] = PushedAtomCoordinates.molec[i].coords.z[j];
        }

        AtomPropertyIsChanging(i, AtomCoordinateChangedListener);
    }

    FreePushedAtomCoordinates();

    return(0);
}

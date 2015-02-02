#include <tcl.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#ifndef NDEBUG
#  define DebugPrint(text)
#else
#  define DebugPrint(text) printf("%s\n",text)
#endif

#include <GL/glu.h>

#include <gopenmolext.h>

#define SEPARATED_DRAWING_CALLBACKS

typedef struct Ball_t {
    /* Ball position. */
    double x,y,z;
    /* Ball size. */
    double radius;
    int   relocatable; /* 0: static position, 1: relocatable */
    struct Ball_t *next;
} Ball_t;

typedef struct {
    gom_PlotterData plotter;
    GLUquadricObj  *gluObject;
    Ball_t         *first_ball;
    gom_AtomCoordinateChangedListener *coordinate_listener;
    gom_MolecStructDeleteListener     *molecstruct_listener;
} PluginData_t;

static PluginData_t PluginData; /* zero inited */

static int AddBall(ClientData ,Tcl_Interp *,int ,Tcl_Obj *CONST[]);
static int RemoveBalls(ClientData,Tcl_Interp *,int, Tcl_Obj *CONST[]);
static int PlotBalls(void *, int, int);
static int UpdateBallCoordinates(void *, unsigned long int);
static int ResetBallInMolecStruct(void *, int, int);
static void ResetBalls(PluginData_t *);

DYNEXPORT_C int Ball_Init(Tcl_Interp *Interp)
{
    printf("%s\n",   "Creating Tcl extensions:");
    printf("\t%s\n", "AddBall Fradius ?Fx Fy Fz?");
    Tcl_CreateObjCommand(Interp,"AddBall",AddBall,NULL,NULL);
    printf("\t%s\n", "RemoveBalls");
    Tcl_CreateObjCommand(Interp,"RemoveBalls",RemoveBalls,NULL,NULL);
    printf("%s\n",   "Finished.");

    return(TCL_OK);
}

int AddBall(ClientData userData,Tcl_Interp *Interp,int objc,Tcl_Obj *CONST objv[])
{
    Ball_t *ball,**pBall;

    /* Check parameters. */
    if ( objc != 2 && objc != 5 ) {
        Tcl_WrongNumArgs(Interp,1,objv,"Fradius ?Fx Fy Fz?");
        return(TCL_ERROR);
    }

    /* Find the end of the ball list. */
    pBall = &PluginData.first_ball;
    while ( *pBall )
        pBall = &(*pBall)->next;
    /* Append a new ball. */
    ball = *pBall = malloc(sizeof(Ball_t));
    if ( ! ball )
        Tcl_Panic("Error out of memory.");

    /* Get values. */
    if ( Tcl_GetDoubleFromObj(Interp,objv[1],&ball->radius) != TCL_OK ||
         ( objc == 5 && (
           Tcl_GetDoubleFromObj(Interp,objv[2],&ball->x) != TCL_OK ||
           Tcl_GetDoubleFromObj(Interp,objv[3],&ball->y) != TCL_OK ||
           Tcl_GetDoubleFromObj(Interp,objv[4],&ball->z) != TCL_OK ) ) ) {
        free(*pBall);
        *pBall = NULL;
        return(TCL_ERROR);
    }

    if ( objc == 2 )
        /* Mark ball as relocatable. */
        ball->relocatable = 1;
    else
        ball->relocatable = 0;
    
    ball->next = NULL;

    /* Add new function to drawing stack if it isn't there already. */
    gom_SetPlotterRegistrationState(
        1 /* ON */, &PluginData.plotter,
        PlotBalls, &PluginData, "ball_plugin",
        gom_PlotOrderBaseNormal);

    if( ! PluginData.gluObject )
        PluginData.gluObject = gluNewQuadric();

    /* Create listeners. */
    if ( ball->relocatable && ! PluginData.coordinate_listener )
        /* We need to know then atom coordinates are changed. */
        PluginData.coordinate_listener =
            gom_AddAtomCoordinateChangedListener(
                UpdateBallCoordinates, &PluginData);
    if ( ! PluginData.molecstruct_listener )
        /* We need to know then gOpenMol is reseted. */
        PluginData.molecstruct_listener =
            gom_AddMolecStructDeleteListener(
                ResetBallInMolecStruct, &PluginData);

    return(TCL_OK);
}

int RemoveBalls(ClientData userData,Tcl_Interp *Interp,int objc,Tcl_Obj *CONST objv[])
{
    ResetBalls(&PluginData);
    return(TCL_OK);
}

int PlotBalls(void *userData, int Wstr, int drawFlags)
{
    PluginData_t *data = userData;
    Ball_t *ball;
    GLint   slices;

    /* We will never draw anything to structure number 1 and up. */
    if ( Wstr != 0 )
        return(-1);

    if( drawFlags & gom_PlotComplexElements ) {
        DebugPrint("Drawing handler: drawing complex balls");

        glEnable(GL_LIGHTING);

        slices = 20;
    }
    else if ( drawFlags & gom_PlotSimplifiedElements ) {
        DebugPrint("Drawing handler: drawing simplified balls");

        glDisable(GL_LIGHTING);

        slices = 2;
    }
    else
        /* We will never draw simple elements. */
        return(-1);

    /* Put to white. */
    glColor4f(1.0 , 1.0 , 1.0 , 1.0);

    /* Plot balls. Rotation matrix is already */
    /* applied by gOpenMol plot system.   */
    for( ball = data->first_ball ; ball ; ball = ball->next ) {
        glPushMatrix();
        glTranslated(ball->x, ball->y, ball->z); 
        gluSphere(data->gluObject, ball->radius, 2*slices, slices);
        glPopMatrix();
    }

    return(0);
}

int UpdateBallCoordinates(void *userData, unsigned long int StructMask)
{
    PluginData_t *data = userData;
    Ball_t *ball;
    int i,NAtoms;

    /* We will never draw anything to structure number 1 and up. */
    if ( ! ( StructMask & ( 1 << 0 ) ) )
        /* The first structure hasn't changed. There is nothing to do. */
        return(0);

    /* Relocate balls. */
    DebugPrint("Drawing handler: relocating balls\n");

    NAtoms = gom_GetNumAtomsInMolecStruct(0);

    ball = data->first_ball;
    for ( i = 0 ; i < NAtoms && ball ; i++ , ball = ball->next ) {
        if( ball->relocatable ) {
            ball->x = gom_GetAtomXCoord(0,i);
            ball->y = gom_GetAtomYCoord(0,i);
            ball->z = gom_GetAtomZCoord(0,i);
        }
    }

    gom_InvalidatePlotterDelayed(&data->plotter);

    return(0);
}

int ResetBallInMolecStruct(void *userData, int Wstr, int Dstr)
{
    /* We will never draw anything to structure number 1 and up. */
    switch ( Wstr ) {
    case -1: /* Full reset */
    case 0:  /* The first structure is about to be deleted or merged. */
        if ( Dstr != 1 ) {
            ResetBalls(userData);
            /* Cancel this listener. */
            PluginData.molecstruct_listener = NULL;
            return(-1);
        }
        break;
/*  default:
        if ( Dstr < 0 )
            / Structure Wstr is deleted. /
            ;
        else
            / Structure Wstr is merged into Dstr /
            ;
*/
    }
    return(0);
}

void ResetBalls(PluginData_t *data)
{
    Ball_t *next_ball;

    while( data->first_ball ) {
        /* Free memory. */
        next_ball = data->first_ball->next;
        free(data->first_ball);
        data->first_ball = next_ball;
    }

    /* Remove function from drawing stack. */
    gom_SetPlotterRegistrationState(
        0 /* OFF */, &PluginData.plotter, NULL, NULL, NULL, 0);

    /* Cancel coordinate listener. */
    if ( data->coordinate_listener )
        gom_CancelAtomCoordinateChangedListener(data->coordinate_listener);
    data->coordinate_listener = NULL;
}

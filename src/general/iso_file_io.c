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

#include "gomstdio.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <tcl.h>

#include <sys/types.h>

#define ISOVIS_MAIN 1
#include "colouring.h"
#include "contour.h"
#include "gomendian.h"
#include "gomfile.h"
#include "gommain.h"
#include "gomstring.h"
#include "memalloc.h"
#include "model_file.h"
#include "molecoord.h"
#include "molecule.h"
#include "printmsg.h"
#include "projview.h"
#include "tclutils.h"

#include "stdafx.h"

#include "isovis.h"

#define Rabs(a)       ( ( a ) > 0.0 ? (a) : -(a))

/* contour structures */
static void DestroyContourLevel(void *,size_t, const DataVectorHandler *);
static void DestroyContourStruct(void *,size_t, const DataVectorHandler *);
const DataVectorHandler gomp_PolygonTypeHandler = { NULL, NULL };
static const Contour_Level NullContourLevel = {
    0.0f, "",
    1.0f, 1.0f, 1.0f,
    CONTOUR_TYPE_SOLID, 1.0f, 0, 1, 0,
    'x', 0.0f,
    0.0f, 0.0f,
    NULL
};
const DataVectorHandler gomp_ContourLevelHandler  = {
    gomp_DataVectorCopyInitFunc,
    DestroyContourLevel,
    &NullContourLevel
};
static const Contour_Level NullContour = { 0 }; /* Zero and NULL init */
const DataVectorHandler gomp_ContourStructHandler = {
    gomp_DataVectorCopyInitFunc,
    DestroyContourStruct,
    &NullContour
};
Contour_Struct *gomp_ContourInfo;

static struct {
    gom_UpdateDataListener* updateDataListener;
    unsigned long int change_mask;
} ContourChangedMask;

static int ContourInternalType; /* = 0, binary, != 0, unformatted */

static void GetMaxMin(const float *data,
                      int xdim, int ydim, int zdim, float *maxp, float *minp);
static int  SCAREgetdims(const char *filename,int *rank, int shape[]);
static int  SCAREgetdata(const char *,int,const int *,float *);

/************************************************************************/
void DestroyContourLevel(
    void *data, size_t size, const DataVectorHandler *pHandler)
/************************************************************************/
{
    Contour_Level *level = data;
    gomp_DataVectorFree(&level->polygons);
}
/************************************************************************/
void DestroyContourStruct(
    void *data, size_t size, const DataVectorHandler *pHandler)
/************************************************************************/
{
    Contour_Struct *info = data;
    gomp_FreeVector(info->data);
    gomp_DataVectorFree(&info->levels);
}
/***************************************************************************/
GOM_IMPLEMENT_GOM_LISTENER(ContourDataChanged)
/***************************************************************************/
static int CallContourDataChangedListener(
    void            *pContourMask,
    gom_ListenerFunc callback,
    void            *callbackData )
{
    return ((gom_ContourDataChangedListenerFunc)callback)(
        callbackData,*(unsigned long int*)pContourMask);
}
/***************************************************************************/
static int CallContourDataChangedListeners(void *data)
/***************************************************************************/
{
    int               result      = 0;
    unsigned long int change_mask = ContourChangedMask.change_mask;
    
    /**
     * Reinitialize the mask and the event. Listeners may change the
     * atom state so we can't do that after we have called the
     * notifiers, anymore.
     */
    memset(&ContourChangedMask,0,sizeof(ContourChangedMask));

    /**
     * Inform that contour graphics is invalidated.
     * This could be part of ContourPropertyIsChanging macro. But to
     * reduced the work needed by ContourPropertyIsChanging it is changed
     * to here.
     */
    gomp_InvalidateContourPlotter(change_mask);
    /* Call contour listeners. */
    result = gomp_CallListeners(
        GOM_GET_LISTENER_LIST(ContourDataChanged),
        CallContourDataChangedListener,
        &ContourChangedMask.change_mask);
    if ( result != 0 )
        /* Error occured. */
        return result;
    /* Cancel this listener. */
    return(-1);
}
/***************************************************************************/
/* Mark the contour data changed. Note that contour graphics is            */
/* invalidated by CallContourListeners listener.                           */
/***************************************************************************/
#define ContourIsChanging(IContour)                                           \
{                                                                             \
    ContourChangedMask.change_mask |= 1 << IContour;                          \
    if ( ! ContourChangedMask.updateDataListener )                           \
        ContourChangedMask.updateDataListener =                               \
            gomp_AddUpdateDataListener(CallContourDataChangedListeners,NULL); \
}

/************************************************************************/
int  gomp_ContourIsChanging(int IContour)
/************************************************************************/
{
    int min = IContour, max = IContour;
    int i;
    if ( IContour < 0 ) {
        min = 0;
        max = gomp_GetContoursDefined() - 1;
    }
    for ( i = min ; i <= max ; i++ ) {
        ContourIsChanging(i);
    }
    return(0);
}

/************************************************************************/
int  gomp_ContourDriver(const char *ContFile, const char *ContName)
/************************************************************************/
{
    static float *Cdata;
    static int    Cxdim;
    static int    Cydim;
    static int    Czdim;
    static float  Cmax;
    static float  Cmin;
    static int    LastContour;
    static char   OutText[BUFF_LEN];

/* reserve space for contour ... */
    if ( ! gomp_DataVectorCreateOrAppendExact(
             &ContourInfo,&gomp_ContourStructHandler,1) )
        return(-1);

    LastContour = gomp_GetContoursDefined() - 1;

/* get the data from file ... */

    if ( gomp_GetContourData(ContFile,
                             &Cdata,
                             &Cxdim,
                             &Cydim,
                             &Czdim) == -1 ) {
        (void)gomp_DataVectorSetExactSize(&ContourInfo,LastContour);
        gomp_PrintMessage("Error from get_data");
        return(-1);
    }

/* save grabbed information           */
    ContourInfo[LastContour].xdim = Cxdim;
    ContourInfo[LastContour].ydim = Cydim;
    ContourInfo[LastContour].zdim = Czdim;
    ContourInfo[LastContour].data = Cdata;

/* these values can be used later ... */
    ContourInfo[LastContour].Xmin = ISO_xmin;
    ContourInfo[LastContour].Xmax = ISO_xmax;
    ContourInfo[LastContour].Ymin = ISO_ymin;
    ContourInfo[LastContour].Ymax = ISO_ymax;
    ContourInfo[LastContour].Zmin = ISO_zmin;
    ContourInfo[LastContour].Zmax = ISO_zmax;

    gomp_CopyString(ContourInfo[LastContour].ContFile,ContFile,BUFF_LEN);

    GetMaxMin(
        ContourInfo[LastContour].data,
        ContourInfo[LastContour].xdim,
        ContourInfo[LastContour].ydim,
        ContourInfo[LastContour].zdim,
        &Cmax,
        &Cmin);

    ContourInfo[LastContour].max = Cmax;
    ContourInfo[LastContour].min = Cmin;

    sprintf(OutText,"Minimum value %f maximum value %f",Cmin,Cmax);
    gomp_PrintMessage(OutText);

    ContourInfo[LastContour].ProjectionIndex = LastContour;

/* ready with the file ...  */

    if(ContName[0] == '\0') { /* no contour name defined */
        sprintf(OutText,"%d",LastContour+1);
        gomp_CopyString(ContourInfo[LastContour].Name,OutText,BUFF_LEN);
    } else {
        (void)gomp_ParseContourName(ContName);
        gomp_CopyString(ContourInfo[LastContour].Name,ContName,BUFF_LEN);
    }

    if(gomp_AddContour2StructureMappingSpace()) {
        (void)gomp_DataVectorSetExactSize(&ContourInfo,LastContour);
        gomp_PrintERROR("Error from adding space to 'contour to structure' mapping");
        return(-1);
    }

    if(gomp_SetContour2StructureMapping( LastContour, 0)) {
        (void)gomp_DataVectorSetExactSize(&ContourInfo,LastContour);
        gomp_PrintERROR("Error from 'set contour to structure' mapping");
        return(-1);
    }

    ContourIsChanging(LastContour);
    
    return(0);
}

/************************************************************************/
int  gomp_FillContourStructure(const char *ContFile, const char *ContName,
                               int   Cxdim, int   Cydim, int Czdim,
                               float Cxmin, float Cxmax,
                               float Cymin, float Cymax,
                               float Czmin, float Czmax,
                               float *Cdata)
/************************************************************************/
{
    static float  Cmax;
    static float  Cmin;
    static int    LastContour;
    static char   OutText[BUFF_LEN];

/* reserve space for contour ... */
    if ( ! gomp_DataVectorCreateOrAppend(
             &ContourInfo,&gomp_ContourStructHandler,1) )
        return(-1);

    LastContour = gomp_GetContoursDefined() - 1;

/* save structure information           */
    ContourInfo[LastContour].xdim = Cxdim;
    ContourInfo[LastContour].ydim = Cydim;
    ContourInfo[LastContour].zdim = Czdim;
    ContourInfo[LastContour].data = Cdata;
        

/* these values can be used later ... */
    ContourInfo[LastContour].Xmin = Cxmin;
    ContourInfo[LastContour].Xmax = Cxmax;
    ContourInfo[LastContour].Ymin = Cymin;
    ContourInfo[LastContour].Ymax = Cymax;
    ContourInfo[LastContour].Zmin = Czmin;
    ContourInfo[LastContour].Zmax = Czmax;

    gomp_CopyString(ContourInfo[LastContour].ContFile,ContFile,BUFF_LEN);

    GetMaxMin(
        ContourInfo[LastContour].data,
        ContourInfo[LastContour].xdim,
        ContourInfo[LastContour].ydim,
        ContourInfo[LastContour].zdim,
        &Cmax,
        &Cmin);

    ContourInfo[LastContour].max = Cmax;
    ContourInfo[LastContour].min = Cmin;

    sprintf(OutText,"Minimum value %f maximum value %f",Cmin,Cmax);
    gomp_PrintMessage(OutText);

    ContourInfo[LastContour].ProjectionIndex = LastContour;

/* ready with the file ...  */

    if(ContName[0] == '\0') { /* no contour name defined */
        sprintf(OutText,"%d",LastContour+1);
        gomp_CopyString(ContourInfo[LastContour].Name,OutText,BUFF_LEN);
    } else {
        (void)gomp_ParseContourName(ContName);
        gomp_CopyString(ContourInfo[LastContour].Name,ContName,BUFF_LEN);
    }

    if(gomp_AddContour2StructureMappingSpace()) {
        (void)gomp_DataVectorSetExactSize(&ContourInfo,LastContour);
        gomp_PrintERROR("Error from adding space to 'contour to structure' mapping");
        return(-1);
    }

    if(gomp_SetContour2StructureMapping( LastContour, 0)) {
        (void)gomp_DataVectorSetExactSize(&ContourInfo,LastContour);
        gomp_PrintERROR("Error from 'set contour to structure' mapping");
        return(-1);
    }

    ContourIsChanging(LastContour);
    
    return(0);
}
#if defined(JUNK)
/************************************************************************/
struct ContourInfo  *gomp_GetContourStructurePointer(int which)
/************************************************************************/
{
    if(which < 1) {
        gomp_PrintERROR("contour structure index < 0");
        return((struct ContourInfo *)NULL);
    }
    if(which >= ContoursDefined) {
        gomp_PrintERROR("contour structure index > contours defined");
        return((struct ContourInfo *)NULL);
    }

    return((struct ContourInfo *) &ContourInfo[which]);
}
#endif
/**************************** get_data ****************************/
/**************************** get_data ****************************/
/**************************** get_data ****************************/
/**************************** get_data ****************************/

/* This subroutine will read in the data file. */
int gomp_GetContourData(const char *filename,
                        float **data,
                        int *xdim, int *ydim, int *zdim)
{

    int rank;
    int shape[3];
    int size;
    int ret;
    char OutText[BUFF_LEN];

    ret=SCAREgetdims(filename,&rank,shape);
    if (ret != 0) {
        sprintf(OutText,"%s: error from SCAREgetdims %d for %s",
                MY_NAME,ret,filename);
        gomp_PrintMessage(OutText);
        fclose(ISO_file);
        return -1;
    }

    if (rank != 3) {
        sprintf(OutText,"%s: rank error from SCAREgetdims %d for %s",
                MY_NAME,ret,filename);
        gomp_PrintMessage(OutText);
        fclose(ISO_file);
        return -1;
    }

    *xdim = shape[0];  *ydim = shape[1];  *zdim = shape[2];
    if (VERBOSE)
        printf("%s: data set size xdim=%d ydim=%d zdim=%d\n",
               MY_NAME,*xdim,*ydim,*zdim);

    size = (*xdim) * (*ydim) * (*zdim) * sizeof(float);

    if ((*data = malloc(size)) == NULL) {
        sprintf(OutText,"%s: error, not enough memory for the data set\n",MY_NAME);
        gomp_PrintMessage(OutText);
        sprintf(OutText,"Need space for '%d * %d * %d = %d' entries",
                *xdim,*ydim,*zdim,
                (*xdim)*(*ydim)*(*zdim));
        gomp_PrintMessage(OutText);
        fclose(ISO_file);
        return -1;
    }


    ret=SCAREgetdata(filename,rank,shape,*data);
    if (ret != 0) {
        sprintf(OutText,"%s: error from SCAREgetdata %d file %s",
                MY_NAME,ret,filename);
        gomp_PrintMessage(OutText);
        return -1;
    }

    return 0;
}

#define RECORD()   { Items = \
                   fread(&record, sizeof(int) , 1L , ISO_file);\
                   if(Items < 1) {\
                     gomp_PrintMessage("?ERROR1 - in reading contour file (*)");\
                     fclose(ISO_file);\
                     return(1);}}

#define FREAD(value_p , size)    { Items = \
                                 fread(value_p, size , 1L , ISO_file);\
                                 if(Items < 1) {\
                     gomp_PrintMessage("?ERROR2 - in reading contour file (**)");\
                     fclose(ISO_file);\
                     return(2);}}

#define FREADN(value_p , num , size) { Items = \
                                fread(value_p, size , num , ISO_file);\
                   if(Items != num) {\
                     gomp_PrintMessage("?ERROR3 - in reading contour file (***)");\
                     fclose(ISO_file);\
                     return(3);}}


/**************************** SCAREgetdims ****************************/
/**************************** SCAREgetdims ****************************/
/**************************** SCAREgetdims ****************************/
/**************************** SCAREgetdims ****************************/
int SCAREgetdims(const char *filename, int *rank, int shape[])
{
    char    OutText[BUFF_LEN];
    int     Items;
    int     record;
    int     test_byte;
    int     LastContour = gomp_GetContoursDefined() - 1;

    ISO_file = fopen(filename,"rb");
    if(ISO_file == NULL) {
        sprintf(OutText,"?ERROR - can't open file contour file '%s'",
                filename);
        gomp_PrintERROR(OutText);
        return(1);
    }

    ContourIsChanging(LastContour);
    
/* first position is the rank (if this is a binary file) */
    FREAD(rank, sizeof(int));

    test_byte = *rank;
    gomp_Reverse_int( & test_byte );
    gomContourSwapBytes = 0;

/* pure binary file */
    if(*rank == 3) { /* this is a binary file */
        ContourInternalType = 0;

/* type of surface            */
        FREAD(&TypeOfSurface, sizeof(int));

        if(VERBOSE) {
            switch(TypeOfSurface) {
            case 0:
                gomp_PrintMessage("?WARNING - unknown surface type ");
                break;
            case 1:
                gomp_PrintMessage("=> Reading a VSS surface ");
                break;
            case 2:
                gomp_PrintMessage("=> Reading an orbital/density surface");
                break;
            case 3:
                gomp_PrintMessage("=> Reading a probe surface");
                break;
            }
        }

/* read shape ...             */
        FREAD(&shape[2] , sizeof(int));
        FREAD(&shape[1] , sizeof(int));
        FREAD(&shape[0] , sizeof(int));

        if(VERBOSE)
            printf("zdim: %d , ydim: %d , xdim: %d\n",shape[2],shape[1],shape[0]);
/* done with shape */

/* read min/max ... */
        FREAD(&ISO_zmin , sizeof(float));
        FREAD(&ISO_zmax , sizeof(float));
        FREAD(&ISO_ymin , sizeof(float));
        FREAD(&ISO_ymax , sizeof(float));
        FREAD(&ISO_xmin , sizeof(float));
        FREAD(&ISO_xmax , sizeof(float));
    }
/* this is a fortran unformatted file */
    else if(*rank == 2 * sizeof(int)) { /* this is a file produced with fortran */

        ContourInternalType = 1;

        FREAD(rank, sizeof(int));
        FREAD(&TypeOfSurface, sizeof(int));

        RECORD();   /* control record */

        RECORD();   /* control record */
        FREAD(&shape[2] , sizeof(int));
        FREAD(&shape[1] , sizeof(int));
        FREAD(&shape[0] , sizeof(int));
        RECORD();   /* control record */

/* read min/max ... */
        RECORD();   /* control record */
        FREAD(&ISO_zmin , sizeof(float));
        FREAD(&ISO_zmax , sizeof(float));
        FREAD(&ISO_ymin , sizeof(float));
        FREAD(&ISO_ymax , sizeof(float));
        FREAD(&ISO_xmin , sizeof(float));
        FREAD(&ISO_xmax , sizeof(float));
        RECORD();   /* control record */
    }
    else if(test_byte == 3) { /* byte swap */

        gomp_PrintMessage("Enabling automatic byte_swapping...");
        gomContourSwapBytes = 1;
        ContourInternalType = 0;
        *rank = 3;
        FREAD(&TypeOfSurface, sizeof(int));
        gomp_Reverse_int( &TypeOfSurface );
/* read shape ...             */
        FREAD(&shape[2] , sizeof(int));
        FREAD(&shape[1] , sizeof(int));
        FREAD(&shape[0] , sizeof(int));

        gomp_Reverse_int_array( shape , 3 );

        if(VERBOSE)
            printf("zdim: %d , ydim: %d , xdim: %d\n",shape[2],shape[1],shape[0]);
/* done with shape */

/* read min/max ... */
        FREAD(&ISO_zmin , sizeof(float));
        gomp_Reverse_float( &ISO_zmin );
        FREAD(&ISO_zmax , sizeof(float));
        gomp_Reverse_float( &ISO_zmax );
        FREAD(&ISO_ymin , sizeof(float));
        gomp_Reverse_float( &ISO_ymin );
        FREAD(&ISO_ymax , sizeof(float));
        gomp_Reverse_float( &ISO_ymax );
        FREAD(&ISO_xmin , sizeof(float));
        gomp_Reverse_float( &ISO_xmin );
        FREAD(&ISO_xmax , sizeof(float));
        gomp_Reverse_float( &ISO_xmax );
    }
    else if(test_byte == 2 * sizeof(int)) {

        gomp_PrintMessage("Enabling automatic byte_swapping...");
        gomContourSwapBytes = 1;
        ContourInternalType = 1;

        FREAD(rank, sizeof(int));
        gomp_Reverse_int( rank );
        FREAD(&TypeOfSurface, sizeof(int));
        gomp_Reverse_int( & TypeOfSurface );

        RECORD();   /* control record */

        RECORD();   /* control record */
        FREAD(&shape[2] , sizeof(int));
        FREAD(&shape[1] , sizeof(int));
        FREAD(&shape[0] , sizeof(int));
        RECORD();   /* control record */
        gomp_Reverse_int_array( shape , 3 );

/* read min/max ... */
        RECORD();   /* control record */
        FREAD(&ISO_zmin , sizeof(float));
        gomp_Reverse_float( &ISO_zmin );
        FREAD(&ISO_zmax , sizeof(float));
        gomp_Reverse_float( &ISO_zmax );
        FREAD(&ISO_ymin , sizeof(float));
        gomp_Reverse_float( &ISO_ymin );
        FREAD(&ISO_ymax , sizeof(float));
        gomp_Reverse_float( &ISO_ymax );
        FREAD(&ISO_xmin , sizeof(float));
        gomp_Reverse_float( &ISO_xmin );
        FREAD(&ISO_xmax , sizeof(float));
        gomp_Reverse_float( &ISO_xmax );
        RECORD();   /* control record */
    } else {
        gomp_PrintERROR("can't determine the internal structure of the plt file");
        return(1);
    }

/* make some conversion (if necessary) */
    if(TypeOfSurface == 100) {  /* OpenMol data in atomic units */
        ISO_zmin *= ATM2ANG;
        ISO_zmax *= ATM2ANG;
        ISO_ymin *= ATM2ANG;
        ISO_ymax *= ATM2ANG;
        ISO_xmin *= ATM2ANG;
        ISO_xmax *= ATM2ANG;
    }

    if(VERBOSE)
        printf("x: %f %f \ny: %f %f \nz: %f %f\n",
               ISO_xmin,
               ISO_xmax,
               ISO_ymin,
               ISO_ymax,
               ISO_zmin,
               ISO_zmax);
          
/* done with min/max */

/* the translation from the centering of the system is taken into
   account at the plotting stage (g_contour.c) */
    ISO_StepX   = (ISO_xmax - ISO_xmin) / 2.0 + ISO_xmin;
    ISO_StepY  = (ISO_ymax - ISO_ymin) / 2.0 + ISO_ymin;
    ISO_StepZ = (ISO_zmax - ISO_zmin) / 2.0 + ISO_zmin;

    if(VERBOSE)
        printf("Translate %f %f %f \n",ISO_StepX,ISO_StepY,ISO_StepZ);

    ContourInfo[LastContour].Xtrans   = ISO_StepX;
    ContourInfo[LastContour].Ytrans  = ISO_StepY;
    ContourInfo[LastContour].Ztrans = ISO_StepZ;

    ISO_StepX   = (ISO_xmax - ISO_xmin) / ((float) (shape[0] - 1));
    ISO_StepY  = (ISO_ymax - ISO_ymin) / ((float) (shape[1] - 1));
    ISO_StepZ = (ISO_zmax - ISO_zmin) / ((float) (shape[2] - 1));

    if(VERBOSE)
        printf("Scale %f %f %f \n",ISO_StepX,ISO_StepY,ISO_StepZ);

    ContourInfo[LastContour].Xscale   = ISO_StepX;
    ContourInfo[LastContour].Yscale  = ISO_StepY;
    ContourInfo[LastContour].Zscale = ISO_StepZ;

/* shift further half the step length */
    ContourInfo[LastContour].Xtrans   += ISO_StepX / 2.0;
    ContourInfo[LastContour].Ytrans  += ISO_StepY / 2.0;
    ContourInfo[LastContour].Ztrans += ISO_StepZ / 2.0;
      
    return(0);
}

/**************************** SCAREgetdata ****************************/
/**************************** SCAREgetdata ****************************/
/**************************** SCAREgetdata ****************************/
/**************************** SCAREgetdata ****************************/
int SCAREgetdata(const char *filename,int rank,const int *shape,float *data)
{   
    int Items;
    int il;
    int size;
    int record;
    float *Temp;

/* unformatted file */
    if(ContourInternalType) {
        size = shape[1] * shape[0];
        Temp = data;
        for(il = 0 ; il < shape[2] ; il++) {
            RECORD();   /* control record */
            FREADN(Temp , size , sizeof(float));
            RECORD();   /* control record */
            Temp += size;
        }
        if(gomContourSwapBytes) {
            gomp_Reverse_float_array( data , size * shape[2] );
        }
    } else {
/* pure binary file */
        size = shape[2] * shape[1] * shape[0];

        FREADN(data , size , sizeof(float));

        if(gomContourSwapBytes) {
            gomp_Reverse_float_array( data , size );
        }
    }

    fclose(ISO_file);

    return(0);
}


/************************************************************************/
int  gomp_ParseContourLevels(int Contour , float Value ,
                             float RedC  ,
                             float GreenC,
                             float BlueC , int Which)
/************************************************************************/
{
    if (gomp_GetContourLevels(Contour) <= Which) {
        /* We need more levels. */
        if ( ! gomp_DataVectorCreateOrAppend(
                 &ContourInfo[Contour].levels,
                 &gomp_ContourLevelHandler,
                 Which + 1 - gomp_GetContourLevels(Contour)) )
            /* error out of memory */
            return(-1);
    }

    ContourIsChanging(Contour);
    
    ContourInfo[Contour].levels[Which].ColVal = Value;
    ContourInfo[Contour].levels[Which].RedC   = RedC;
    ContourInfo[Contour].levels[Which].GreenC = GreenC;
    ContourInfo[Contour].levels[Which].BlueC  = BlueC;

/*      ContourInfo[Level].NumVal++; */
/*      ContourInfo[Level].Display   = 1; */
    strcpy(ContourInfo[Contour].levels[Which].ColNam, "**not avail**"); 

    return(0);
}
/************************************************************************/
int  gomp_FinalizeContourLevels(int Contour , int LevelCount)
/************************************************************************/
{
    /* Not a real change.        */
    /*ContourIsChanging(Contour);*/
    
    /* Free unused levels. */
    (void)gomp_DataVectorSetExactSize(&ContourInfo[Contour].levels,LevelCount);

    return(0);
}
/*
  WARNING! These functions contain no index check. That has to be
  done outside from here!

  LUL Sepetember 1994
*/
/************************************************************************/
float gomp_GetContourMin(int Which)
/************************************************************************/
{
    if(Which >= 0 && Which < gomp_GetContoursDefined())
        return(ContourInfo[Which].min);
    else {
        gomp_PrintERROR("index out of range in 'gomp_GetContourMin'");
        return((float)0.0);
    }
}
      
/************************************************************************/
float gomp_GetContourMax(int Which)
/************************************************************************/
{
    if(Which >= 0 && Which < gomp_GetContoursDefined())
        return(ContourInfo[Which].max);
    else {
        gomp_PrintERROR("index out of range in 'gomp_GetContourMax'");
        return((float)0.0);
    }

}
/************************************************************************/
int   gomp_ContourSmoothON(int Which , int Level)
/************************************************************************/
{
    if( ( Which >= 0 && Which < gomp_GetContoursDefined() ) &&
        Level < gomp_GetContourLevels(Which) ) {
        ContourIsChanging(Which);
        ContourInfo[Which].levels[Level].ContSmooth = 1;
        return(0);
    }
    else {
        gomp_PrintERROR("index out of range in 'gomp_ContourSmoothON'");
        return(1);
    }
}
/************************************************************************/
int   gomp_ContourSmoothOFF(int Which , int Level)
/************************************************************************/
{
    if( ( Which >= 0 && Which < gomp_GetContoursDefined() ) &&
        Level < gomp_GetContourLevels(Which) ) {
        ContourIsChanging(Which);
        ContourInfo[Which].levels[Level].ContSmooth = 0;
        return(0);
    }
    else {
        gomp_PrintERROR("index out of range in 'gomp_ContourSmoothOFF'");
        return(1);
    }
}
/************************************************************************/
int  gomp_ContourDisplayON(int Which , int Level)
/************************************************************************/
{
    if ( ( Which >= 0 && Which < gomp_GetContoursDefined() ) &&
         Level < gomp_GetContourLevels(Which) ) {
        ContourIsChanging(Which);
        ContourInfo[Which].levels[Level].Display = 1;
        return(0);
    }
    else {
        gomp_PrintERROR("index out of range in 'gomp_ContourDisplayON'");
        return(1);
    }
}
/************************************************************************/
int  gomp_ContourDisplayOFF(int Which , int Level)
/************************************************************************/
{
    if ( ( Which >= 0 && Which < gomp_GetContoursDefined() ) &&
         Level < gomp_GetContourLevels(Which) ) {
        ContourIsChanging(Which);
        ContourInfo[Which].levels[Level].Display = 0;
        return(0);
    }
    else {
        gomp_PrintERROR("index out of range in 'gomp_ContourDisplayOFF'");
        return(1);
    }
}
/************************************************************************/
int   gomp_SetContourDisplayType(int Which, int Level, int Value)
/************************************************************************/
{
    if ( ( Which >= 0 && Which < gomp_GetContoursDefined() ) &&
         Level < gomp_GetContourLevels(Which) ) {
        ContourIsChanging(Which);
        ContourInfo[Which].levels[Level].DisplayType = Value;
        return(0);
    }
    else {
        gomp_PrintERROR("index out of range in 'gomp_SetContourDisplayType'");
        return(1);
    }
}
/************************************************************************/
int   gomp_GetContourDisplayType(int Which, int Level)
/************************************************************************/
{
    if((Which >= 0 && Which < gomp_GetContoursDefined()) &&
       Level < gomp_GetContourLevels(Which))
        return ContourInfo[Which].levels[Level].DisplayType;
    else
        gomp_PrintERROR("index out of range in 'gomp_GetContourDisplayType'");

    return(0);
}

/************************************************************************/
int   gomp_SetContourDisplayTypeGlobal(int Value)
/************************************************************************/
{
    SurfaceDisplayType = Value;

    gomp_ContourIsChanging(-1);

    return(0);      
}
/************************************************************************/
int   gomp_GetContourDisplayTypeGlobal()
/************************************************************************/
{
    return(SurfaceDisplayType);
}
/************************************************************************/
int   gomp_GetContoursDefined()
/************************************************************************/
{
    return(gomp_DataVectorGetSize(&ContourInfo));
}
/************************************************************************/
const char *gomp_GetContourFileName(int Which)
/************************************************************************/
{
    if(Which >= 0 && Which < gomp_GetContoursDefined())
        return(ContourInfo[Which].ContFile);
    else
        return((const char *)NULL);
}
/************************************************************************/
int   gomp_CheckContourName(const char *Name)
/************************************************************************/
{
    int i;

    if(!gomp_GetContoursDefined()) {
        gomp_PrintERROR("no contours defined");
        return(0);
    }

/* check first as a character ...*/
    for(i = 0 ; i < gomp_GetContoursDefined() ; i++) {
        if(gomp_StringMatch(ContourInfo[i].Name , Name)) {
            return(i+1);
        }
    }
/* no hit as a character see if it is a number */
    if(isdigit(Name[0])) {
        i = atoi(Name);
        if(i < 1 || i > gomp_GetContoursDefined()) {
            gomp_PrintERROR("contour index out of range");
            return(0);
        }
        return(i);
    }
/* no match */
    gomp_FormatERROR("no contour name match '%s'",Name);
    return(0);
}
/************************************************************************/
int   gomp_ParseContourName(const char *Name)
/************************************************************************/
{
    int i;

    if(!gomp_GetContoursDefined())
        return(0);

    for(i = 0 ; i < gomp_GetContoursDefined() ; i++) {
        if(gomp_StringMatch(ContourInfo[i].Name , Name)) {
            gomp_PrintWARNING("contour name already in use");
            return(i+1);
        }
    }
/* no name match */
    return(0);
}
/************************************************************************/
int   gomp_GetContourLevels(int Which)
/************************************************************************/
{
    return(gomp_DataVectorGetSize(&ContourInfo[Which].levels));
}
/************************************************************************/
int   gomp_SetContourAlpha(int Which ,int Level ,float Value)
/************************************************************************/
{
    ContourInfo[Which].levels[Level].AlphaBlend = Value;

    ContourIsChanging(Which);

    return(0);
}
/************************************************************************/
float gomp_GetContourAlpha(int Which ,int Level)
/************************************************************************/
{
    return(ContourInfo[Which].levels[Level].AlphaBlend);
}
/************************************************************************/
int   gomp_SetContourCullFace(int Which ,int Level ,int Value)
/************************************************************************/
{
    ContourInfo[Which].levels[Level].CullFace = Value;

    ContourIsChanging(Which);

    return(0);
}
/************************************************************************/
int  gomp_GetContourCullFace(int Which ,int Level)
/************************************************************************/
{
    return(ContourInfo[Which].levels[Level].CullFace);
}
/************************************************************************/
int   gomp_WriteContourInfo2ModelFile(FILE *Model_f)
/************************************************************************/
{
    int   i,j,k,lb;

    if(!gomp_GetContoursDefined()) /* no contours defined */
        return(0);

/* CONTOUR - tag */
    fprintf(Model_f , "[Contour]\n");
    fprintf(Model_f , "%d\n",gomp_GetContoursDefined());
    fprintf(Model_f , "%d\n",gomp_GetColorMappingType());

    for(i = 0 ; i < gomp_GetContoursDefined() ; i++) {
        fprintf(Model_f,"%s\n",ContourInfo[i].ContFile);
        fprintf(Model_f,"%s\n",ContourInfo[i].Name);
        fprintf(Model_f,"%f %f %d %d %d\n",
                ContourInfo[i].min,ContourInfo[i].max,
                ContourInfo[i].xdim,ContourInfo[i].ydim,ContourInfo[i].zdim);

        fprintf(Model_f,"%f %f %f %f %f %f\n",
                ContourInfo[i].Xmin,ContourInfo[i].Xmax,
                ContourInfo[i].Ymin,ContourInfo[i].Ymax,
                ContourInfo[i].Zmin,ContourInfo[i].Zmax);

/*
  fprintf(Model_f,"%f %f %f %f %f %f\n",
  ContourInfo[i].Xtrans,ContourInfo[i].Ytrans,ContourInfo[i].Ztrans,
  ContourInfo[i].Xscale,ContourInfo[i].Yscale,ContourInfo[i].Zscale);
*/
        fprintf(Model_f,"%d\n",ContourInfo[i].ProjectionIndex);
        fprintf(Model_f,"%d\n",gomp_GetContour2StructureMapping(i));
        fprintf(Model_f,"%d\n",gomp_GetContourLevels(i));

        for(j = 0 ; j < gomp_GetContourLevels(i) ; j++) {
            fprintf(Model_f,"%s\n",ContourInfo[i].levels[j].ColNam);
            fprintf(Model_f,"%f %f %f %f %d %f %d %d\n",
                    ContourInfo[i].levels[j].ColVal      ,
                    ContourInfo[i].levels[j].RedC        ,
                    ContourInfo[i].levels[j].GreenC      ,
                    ContourInfo[i].levels[j].BlueC       ,
                    ContourInfo[i].levels[j].DisplayType ,
                    ContourInfo[i].levels[j].AlphaBlend  ,
                    ContourInfo[i].levels[j].ContSmooth  ,
                    ContourInfo[i].levels[j].Display);

            fprintf(Model_f,"%f %f %d %c %f\n",
                    ContourInfo[i].levels[j].FScaleMin , 
                    ContourInfo[i].levels[j].FScaleMax , 
                    ContourInfo[i].levels[j].CullFace  ,
                    ContourInfo[i].levels[j].ClipAxis,
                    ContourInfo[i].levels[j].ClipPos);
        }
        k  = ContourInfo[i].xdim * ContourInfo[i].ydim * ContourInfo[i].zdim;
        lb = 0;
        for(j = 0 ; j < k ; j++) {
            fprintf(Model_f,"%f ",ContourInfo[i].data[j]);
            lb++;
            if(lb == 8) {
                fprintf(Model_f,"\n");
                lb =  0;
            }
        }
        if(lb)
            fprintf(Model_f,"\n");
    }

    return(0);
}
/************************************************************************/
int   gomp_ReadContourInfo2ModelFile(FILE *Model_f)
/************************************************************************/
{
    int   i,j,k;
    char  InputText[BUFF_LEN];
    int   ITemp,ContoursDefined,ContourLevelsDefined;
    int   ColourMappingType;

    if(gomp_GetContoursDefined())
        (void)gomp_DeleteAllContours(); /* in contdriver.c */ 

    (void)gomp_DeleteGetContour2StructureMapping();

/* CONTOUR - tag */
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText, "%d",&ContoursDefined);

/* reserve space for contour ... */
    if ( ! gomp_DataVectorCreate(
             &ContourInfo,
             &gomp_ContourStructHandler,
             ContoursDefined) )
        return(-1);

    if(gomp_ModelFileVersion >= 5) {
        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText, "%d",&ColourMappingType);
        gomp_SetColorMappingType(ColourMappingType);
    }

    for(i = 0 ; i < ContoursDefined ; i++) {

        if(gomp_AddContour2StructureMappingSpace()) {
            gomp_PrintERROR("Error from adding space to 'contour to structure' mapping");
            return(-1);
        }

        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        gomp_CopyString(ContourInfo[i].ContFile,InputText,BUFF_LEN);
        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        gomp_CopyString(ContourInfo[i].Name,InputText,BUFF_LEN);
        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText,"%f %f %d %d %d",
               &ContourInfo[i].min,&ContourInfo[i].max,
               &ContourInfo[i].xdim,&ContourInfo[i].ydim,&ContourInfo[i].zdim);

        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText,"%f %f %f %f %f %f",
               &ContourInfo[i].Xmin,&ContourInfo[i].Xmax,
               &ContourInfo[i].Ymin,&ContourInfo[i].Ymax,
               &ContourInfo[i].Zmin,&ContourInfo[i].Zmax);
    
        if(gomp_ModelFileVersion < 5) {
            gomp_Fgets(InputText,BUFF_LEN,Model_f);
            sscanf(InputText,"%f %f %f %f %f %f",
                   &ContourInfo[i].Xtrans,&ContourInfo[i].Ytrans,&ContourInfo[i].Ztrans,
                   &ContourInfo[i].Xscale,&ContourInfo[i].Yscale,&ContourInfo[i].Zscale);
        }

        if(gomp_ModelFileVersion >= 5) {
            gomp_Fgets(InputText,BUFF_LEN,Model_f);
            sscanf(InputText,"%d",&ContourInfo[i].ProjectionIndex);
        }

        if(gomp_ModelFileVersion >= 7) {
            gomp_Fgets(InputText,BUFF_LEN,Model_f);
            sscanf(InputText,"%d",&ITemp);
            if(gomp_SetContour2StructureMapping(i , ITemp)) {
                gomp_PrintERROR("can't solve contour mapping");
                gomp_Exit(1);
            }
        }

        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText,"%d",&ContourLevelsDefined);

/* reserve space for contour levels ... */
        if ( ! gomp_DataVectorCreate(
                 &ContourInfo[i].levels,
                 &gomp_ContourLevelHandler,
                 ContourLevelsDefined) ) {
            (void)gomp_DeleteAllContours(); /* in contdriver.c */
            return(-1);
        }

        for(j = 0 ; j < ContourLevelsDefined ; j++) {
            gomp_Fgets(InputText,BUFF_LEN,Model_f);
            gomp_CopyString(ContourInfo[i].levels[j].ColNam,InputText,BUFF_LEN);

            gomp_Fgets(InputText,BUFF_LEN,Model_f);
            sscanf(InputText,"%f %f %f %f %d %f %d %d",
                   &ContourInfo[i].levels[j].ColVal      ,
                   &ContourInfo[i].levels[j].RedC        ,
                   &ContourInfo[i].levels[j].GreenC      ,
                   &ContourInfo[i].levels[j].BlueC       ,
                   &ContourInfo[i].levels[j].DisplayType ,
                   &ContourInfo[i].levels[j].AlphaBlend  ,
                   &ContourInfo[i].levels[j].ContSmooth  ,
                   &ContourInfo[i].levels[j].Display);

            if(gomp_ModelFileVersion >= 5) {
                gomp_Fgets(InputText,BUFF_LEN,Model_f);
                sscanf(InputText,"%f %f %d %c %f",
                       &ContourInfo[i].levels[j].FScaleMin , 
                       &ContourInfo[i].levels[j].FScaleMax , 
                       &ContourInfo[i].levels[j].CullFace,
                       &ContourInfo[i].levels[j].ClipAxis,
                       &ContourInfo[i].levels[j].ClipPos);
            }
        }

        k  = ContourInfo[i].xdim * ContourInfo[i].ydim * ContourInfo[i].zdim;

        ContourInfo[i].data = gomp_AllocateFloatVector(k);

        for(j = 0 ; j < k ; j++) {
            fscanf(Model_f,"%f",&ContourInfo[i].data[j]);
        }
        gomp_Fgets(InputText,BUFF_LEN,Model_f);
    }

    if(gomp_GetTermType() == GRAPHICS_AVAILABLE)
        (void) gomp_SetSurfaceControlON();

    return(0);
}
/************************************************************************/
int   gomp_ShowContourInfo()
/************************************************************************/
{
    int  i,j;
    char Text[BUFF_LEN];

    if(!gomp_GetContoursDefined()) {
        gomp_PrintWARNING("no contours defined");
        return(1);
    }

    for(i = 0 ; i < gomp_GetContoursDefined() ; i++) {
        sprintf(Text,"(%d): File: %s, Min: %f, Max: %f, Name: %s",
                (i+1),
                ContourInfo[i].ContFile,
                ContourInfo[i].min,ContourInfo[i].max,
                ContourInfo[i].Name);
        gomp_PrintMessage(Text);
        for(j = 0 ; j < gomp_GetContourLevels(i) ; j++) {
            sprintf(Text,"Iso Value: %f, Red: %f , Green: %f , Blue: %f",
                    ContourInfo[i].levels[j].ColVal,
                    ContourInfo[i].levels[j].RedC,
                    ContourInfo[i].levels[j].GreenC,
                    ContourInfo[i].levels[j].BlueC);
            gomp_PrintMessage(Text);
        }
    }

    return(0);
}

/************************************************************************/
int   gomp_CheckContourIndex(int Input)
/************************************************************************/
{
/*    It's outside the range ... */
    if(Input >= gomp_GetContoursDefined() || Input < 0) {
        return(1);
    }
/*    It's inside the range ...  */
    return(0);
}

/************************************************************************/
float gomp_GetContourValue(int Contour, int Index)
/************************************************************************/
{
    return(ContourInfo[Contour].levels[Index].ColVal);
}
/************************************************************************/
float gomp_GetContourColourRed(int Contour, int Index)
/************************************************************************/
{
    return(ContourInfo[Contour].levels[Index].RedC);
}
/************************************************************************/
float gomp_GetContourColourGreen(int Contour, int Index)
/************************************************************************/
{
    return(ContourInfo[Contour].levels[Index].GreenC);
}
/************************************************************************/
float gomp_GetContourColourBlue(int Contour, int Index)
/************************************************************************/
{
    return(ContourInfo[Contour].levels[Index].BlueC);
}

/************************************************************************/
const char *gomp_GetContourName(int Contour)
/************************************************************************/
{
    return(ContourInfo[Contour].Name);
}
/*
  this is a temporary solutions because this returns only the
  first one. in fact the full range could be set
*/
/************************************************************************/
int   gomp_GetContourSmooth(int Contour , int Level)
/************************************************************************/
{
    return(ContourInfo[Contour].levels[Level].ContSmooth);
}
/************************************************************************/
int   gomp_GetContourDisplayState(int Contour , int Level)
/************************************************************************/
{
    return(ContourInfo[Contour].levels[Level].Display);
}
/************************************************************************/
int   gomp_GetContourCube(int Contour,
                          float *Minx ,float *Maxx ,
                          float *Miny ,float *Maxy ,
                          float *Minz ,float *Maxz ,
                          int *Xdim ,int *Ydim , int *Zdim)
/************************************************************************/
{
    *Minx = ContourInfo[Contour].Xmin;
    *Maxx = ContourInfo[Contour].Xmax;

    *Miny = ContourInfo[Contour].Ymin;
    *Maxy = ContourInfo[Contour].Ymax;

    *Minz = ContourInfo[Contour].Zmin;
    *Maxz = ContourInfo[Contour].Zmax;

    *Xdim = ContourInfo[Contour].xdim;
    *Ydim = ContourInfo[Contour].ydim;
    *Zdim = ContourInfo[Contour].zdim;

    return(0);
}

/************************************************************************/
int   gomp_SetContourProjection(int Master, int Slave)
/************************************************************************/
{
    ContourInfo[Master].ProjectionIndex = Slave;

    ContourIsChanging(Master);

    return(0);
}
/************************************************************************/
int   gomp_GetContourProjection(int Master)
/************************************************************************/
{
    return(ContourInfo[Master].ProjectionIndex);
}

/************************************************************************/
int  gomp_ParseContourLevelsProjection(int Contour , float Value ,
                                       float FMin  ,
                                       float FMax,
                                       int Which)
/************************************************************************/
{
    if (gomp_GetContourLevels(Contour) <= Which) {
        /* We need more levels. */
        if ( ! gomp_DataVectorCreateOrAppend(
                 &ContourInfo[Contour].levels,
                 &gomp_ContourLevelHandler,
                 Which + 1 - gomp_GetContourLevels(Contour)) )
            /* error out of memory */
            return(-1);
    }

    ContourIsChanging(Contour);

    ContourInfo[Contour].levels[Which].ColVal    = Value;
    ContourInfo[Contour].levels[Which].FScaleMin = FMin;
    ContourInfo[Contour].levels[Which].FScaleMax = FMax;
/* default colour applies just to the GUI button colours */
    ContourInfo[Contour].levels[Which].RedC      = 1.0;
    ContourInfo[Contour].levels[Which].GreenC    = 1.0;
    ContourInfo[Contour].levels[Which].BlueC     = 1.0;

/*      ContourInfo[Level].NumVal++; */
/*      ContourInfo[Level].Display   = 1; */
    strcpy(ContourInfo[Contour].levels[Which].ColNam, "**not avail**"); 

    return(0);
}

/************************************************************************/
float gomp_GetContourProjectionMin(int Contour , int Level)
/************************************************************************/
{
    return(ContourInfo[Contour].levels[Level].FScaleMin);
}
/************************************************************************/
float gomp_GetContourProjectionMax(int Contour , int Level)
/************************************************************************/
{
    return(ContourInfo[Contour].levels[Level].FScaleMax);
}
/************************************************************************/
void GetMaxMin(const float *data,
               int xdim, int ydim, int zdim, float *maxp, float *minp)
/* This subroutine finds the maximum & minimum data values */
{
    int   i;
    float max, min;

    max = min = *data;

    for (i = 0; i < (xdim * ydim * zdim); i++) {
        if (data[i] > max) max = data[i];
        if (data[i] < min) min = data[i];
    }
    *maxp = max; *minp = min;
}
/************************************************************************/
const float *gomp_GetContourDataPointer(int which)
/************************************************************************/
{

    if ( which < 0 || which >= gomp_GetContoursDefined() )
        return((const float *)NULL);

    return(ContourInfo[which].data);

}
/************************************************************************/
int   gomp_GetContourPointsX(int Which)
/************************************************************************/
{
    if(Which < 0 || Which > gomp_GetContoursDefined()) {
        gomp_PrintERROR("contour index out of allowed range (points in x direction)");
        return(-1);
    }

    return(ContourInfo[Which].xdim);
}
/************************************************************************/
int   gomp_GetContourPointsY(int Which)
/************************************************************************/
{
    if ( Which < 0 || Which > gomp_GetContoursDefined() ) {
        gomp_PrintERROR("contour index out of allowed range (points in y direction)");
        return(-1);
    }

    return(ContourInfo[Which].ydim);
}
/************************************************************************/
int   gomp_GetContourPointsZ(int Which)
/************************************************************************/
{
    if ( Which < 0 || Which > gomp_GetContoursDefined() ) {
        gomp_PrintERROR("contour index out of allowed range (points in z direction)");
        return(-1);
    }

    return(ContourInfo[Which].zdim);
}
/************************************************************************/
float gomp_GetContourMinX(int Which)
/************************************************************************/
{
    if ( Which < 0 || Which > gomp_GetContoursDefined() ) {
        gomp_PrintERROR("contour index out of allowed range (min x value)");
        return(-1.0);
    }

    return(ContourInfo[Which].Xmin);
}

/************************************************************************/
float gomp_GetContourMaxX(int Which)
/************************************************************************/
{
    if ( Which < 0 || Which > gomp_GetContoursDefined() ) {
        gomp_PrintERROR("contour index out of allowed range (max x value)");
        return(-1.0);
    }

    return(ContourInfo[Which].Xmax);
}

/************************************************************************/
float gomp_GetContourMinY(int Which)
/************************************************************************/
{
    if ( Which < 0 || Which > gomp_GetContoursDefined() ) {
        gomp_PrintERROR("contour index out of allowed range (min y value)");
        return(-1.0);
    }

    return(ContourInfo[Which].Ymin);
}

/************************************************************************/
float gomp_GetContourMaxY(int Which)
/************************************************************************/
{
    if(Which < 0 || Which > gomp_GetContoursDefined()) {
        gomp_PrintERROR("contour index out of allowed range (max y value)");
        return(-1.0);
    }

    return(ContourInfo[Which].Ymax);
}

/************************************************************************/
float gomp_GetContourMinZ(int Which)
/************************************************************************/
{
    if(Which < 0 || Which > gomp_GetContoursDefined()) {
        gomp_PrintERROR("contour index out of allowed range (min z value)");
        return(-1.0);
    }

    return(ContourInfo[Which].Zmin);
}

/************************************************************************/
float gomp_GetContourMaxZ(int Which)
/************************************************************************/
{
    if(Which < 0 || Which > gomp_GetContoursDefined()) {
        gomp_PrintERROR("contour index out of allowed range (max z value)");
        return(-1.0);
    }

    return(ContourInfo[Which].Zmax);
}
/************************************************************************/
int gomp_SetContourLevelClippingPlanePosition(
    int Which, int Level, float pos)
/************************************************************************/
{
    if( Which >= 0 && Which < gomp_GetContoursDefined() &&
        Level >= 0 && Level < gomp_GetContourLevels(Which) ) {
        ContourIsChanging(Which);
        ContourInfo[Which].levels[Level].ClipPos = pos;
        return(0);
    }
    else {
        gomp_PrintERROR("index out of range in 'gomp_SetContourLevelClippingPlanePosition'");
        return(1);
    }
}
/************************************************************************/
int gomp_SetContourLevelClippingPlaneAxis(
    int Which, int Level, char axis)
/************************************************************************/
{
    switch ( axis ) {
    case 'x':
    case 'y':
    case 'z':
        break;
    default:
        gomp_PrintERROR(
            "invalid axis in 'gomp_SetContourLevelClippingPlaneAxis'");
        return(1);
    }
    if( Which >= 0 && Which < gomp_GetContoursDefined() &&
        Level >= 0 && Level < gomp_GetContourLevels(Which) ) {
        ContourIsChanging(Which);
        ContourInfo[Which].levels[Level].ClipAxis = axis;
        return(0);
    }
    else {
        gomp_PrintERROR("index out of range in 'gomp_SetContourLevelClippingPlaneAxis'");
        return(1);
    }
}
/************************************************************************/
float gomp_GetContourLevelClippingPlanePosition(int Which, int Level)
/************************************************************************/
{
    if( Which >= 0 && Which < gomp_GetContoursDefined() &&
        Level >= 0 && Level < gomp_GetContourLevels(Which) )
        return(ContourInfo[Which].levels[Level].ClipPos);
    else {
        gomp_PrintERROR("index out of range in 'gomp_GetContourLevelClippingPlanePosition'");
        return(-1.0);
    }
}
/************************************************************************/
char gomp_GetContourLevelClippingPlaneAxis(int Which, int Level)
/************************************************************************/
{
    if( Which >= 0 && Which < gomp_GetContoursDefined() &&
        Level >= 0 && Level < gomp_GetContourLevels(Which) )
        return(ContourInfo[Which].levels[Level].ClipAxis);
    else {
        gomp_PrintERROR("index out of range in 'gomp_SetContourLevelClippingPlaneAxis'");
        return(' ');
    }
}
/************************************************************************/

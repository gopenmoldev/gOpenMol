/***** PUBLIC GOMAPI BEGIN *****/
/*

                       Copyright (c) 1995 - 2005 by:
        Leif Laaksonen, Center for Scientific Computing, ESPOO, FINLAND
            Confidential unpublished property of Leif Laaksonen
                        All rights reserved


                       in collaboration with

                 OpenMol Molecular Astrophysics Group
                Max-Planck-Institut  fuer Astrophysik
                        Garching, GERMANY
*/
/****** PUBLIC GOMAPI END ******/

#include "gomstdio.h"

#include "memalloc.h"

#define CONTOUR_TYPE_SOLID    0  /* display object as a solid surface       */
#define CONTOUR_TYPE_MESH     1  /* display object as a mesh surface        */
#define CONTOUR_TYPE_LINE     2

#define   CUTPLANE_NO_ACTION      0
#define   CUTPLANE_SQRT_ACTION    1
#define   CUTPLANE_LOG10_ACTION   2

typedef struct
{
    float x[3], y[3], z[3];
    float u[3], v[3], w[3];
    float c[3], f[3];
} polygon_t;

extern const DataVectorHandler gomp_PolygonTypeHandler;

/* contour structure */

typedef struct Contour_Level_ {
    float ColVal;
    char  ColNam[BUFF_LEN];
    float RedC;
    float GreenC;
    float BlueC;
    int   DisplayType;
    float AlphaBlend;
    int   ContSmooth;
    int   Display;
    int   CullFace;
    char  ClipAxis;
    float ClipPos;
/* data needed for the special case of projection of grid data from
   one file upn an other */
    float FScaleMin;
    float FScaleMax;
    polygon_t *polygons; /* contour level polygons */
} Contour_Level;

extern const DataVectorHandler gomp_ContourLevelHandler;

typedef struct Contour_Struct_ {
    float *data;           /* contour data */
    float min;
    float max;
    int   xdim;
    int   ydim;
    int   zdim;
    float Xmin;
    float Xmax;
    float Ymin;
    float Ymax;
    float Zmin;
    float Zmax;
    float Xtrans;
    float Ytrans;
    float Ztrans;
    float Xscale;
    float Yscale;
    float Zscale;
    Contour_Level *levels; /* contour level data */
    char  ContFile[BUFF_LEN];
    char  Name[BUFF_LEN];
    int   ProjectionIndex;
} Contour_Struct;

extern Contour_Struct *gomp_ContourInfo;

#define ContourInfo gomp_ContourInfo

extern const DataVectorHandler gomp_ContourStructHandler;

/* contour functions */

#include "gomdefs.h"
#ifdef ENABLE_EXTENSIONS
#   include "gomlib/gomcontour.h"
#   include "gomext/contour.h"
#else

/***** PUBLIC GOMAPI BEGIN *****/

/** @weakgroup gom_doc_listener
***/

/** @weakgroup gom_doc
***/

/** @defgroup gom_doc_contour Contours
*** @ingroup  gom_doc
***/

/** @name Contour Functions
***/
/** @brief Get contour data.
*** @ingroup  gom_doc_contour
*** @par Example
*** @code
*** const int    PX  = gomp_GetContourPointsX( contour );
*** const int    PY  = gomp_GetContourPointsY( contour );
*** const int    PZ  = gomp_GetContourPointsZ( contour );
*** const float *Ptr = gomp_GetContourDataPointer( contour );
*** int i, j, k;
*** for ( i = 0 ; i < PX ; ++i )
***     for ( j = 0 ; j < PY ; ++j )
***         for ( k = 0 ; k < PZ ; ++k )
***             printf(
***                 "data[%d][%d][%d] = %f\n",
***                 i, j, k, Ptr[i + DX * ( j + DY * k )];
*** @endcode
***/
/* @{ */
extern const float *gomp_GetContourDataPointer( GOM_ARG( int, contour ) );
/* @} */

/** @brief Fill contour data.
*** @ingroup  gom_doc_contour
*** @param contourFile  File name for a record
*** @param contourName  Contour name
*** @param pointsX      Contour grid dimension
*** @param pointsY      Contour grid dimension
*** @param pointsZ      Contour grid dimension
*** @param minX         Contour x dimension and position
*** @param maxX         Contour x dimension and position
*** @param minY         Contour y dimension and position
*** @param maxY         Contour y dimension and position
*** @param minZ         Contour z dimension and position
*** @param maxZ         Contour z dimension and position
*** @param contourData  Contour data.
***                     Must be allocated by @ref gomp_AllocateFloatVector.
***                     Will be managed and released by gOpenMol
***                     after a succcessful call.
*** @retval 0  Success.
*** @sa @ref gomp_AllocateFloatVector
***/
/* @{ */
extern int          gomp_FillContourStructure(
    GOM_ARG( const char *, contourFile ),
    GOM_ARG( const char *, contourName ),
    GOM_ARG( int         , pointsX     ),
    GOM_ARG( int         , pointsY     ),
    GOM_ARG( int         , pointsZ     ),
    GOM_ARG( float       , minX        ),
    GOM_ARG( float       , maxX        ),
    GOM_ARG( float       , minY        ),
    GOM_ARG( float       , maxY        ),
    GOM_ARG( float       , minZ        ),
    GOM_ARG( float       , maxZ        ),
    GOM_ARG( float      *, contourData ) );
/* @} */

/** @brief Contour grid dimension.
*** @ingroup  gom_doc_contour
***/
/* @{ */
extern int          gomp_GetContourPointsX( GOM_ARG( int, contour ) );
extern int          gomp_GetContourPointsY( GOM_ARG( int, contour ) );
extern int          gomp_GetContourPointsZ( GOM_ARG( int, contour ) );
/* @} */
/** @brief Contour dimension and position.
*** @ingroup  gom_doc_contour
***/
/* @{ */
extern float        gomp_GetContourMinX( GOM_ARG( int, contour ) );
extern float        gomp_GetContourMaxX( GOM_ARG( int, contour ) );
extern float        gomp_GetContourMinY( GOM_ARG( int, contour ) );
extern float        gomp_GetContourMaxY( GOM_ARG( int, contour ) );
extern float        gomp_GetContourMinZ( GOM_ARG( int, contour ) );
extern float        gomp_GetContourMaxZ( GOM_ARG( int, contour ) );
/* @} */

/** @weakgroup gom_doc_listener
***/

/** @defgroup gom_doc_listener_contour_data Contour Data Changed Listeners
*** @ingroup gom_doc_listener
*** @ingroup gom_doc_contour
***
*** Contour data changed listeners are called after data update is requested
*** if contour data has been changed.
***
*** @sa @ref gomp_UpdateData.
***
*** @{
***/

/** @brief Contour data changed listener handle type.
***/
GOM_DECLARE_HANDLE(gom_ContourDataChangedListener);

/** @brief Contour data changed listener callback type.
*** @param callbackData  Custom data pointer
*** @param mask          Contour mask
***
*** Note that the mask may have overflowed.
*** To check if a contour may have changed, use something similar to
*** @code
*** int               contour = ...; // contour index
*** unsigned long int bit       = ((unsigned long int)1) << contour;
*** if ( ( contourMask & bit ) || ! bit ) {
***     // The contour may have changed.
*** }
*** @endcode
***/
typedef int (*gom_ContourDataChangedListenerFunc)(
    GOM_ARG( void            *, callbackData  ),
    GOM_ARG( unsigned long int, contourMask ) );

/** @brief Register a new contour data changed listener callback.
***
*** @param callback      Listener callback function
*** @param callbackData  Pointer to be passed to @p listener
*** @return       New listener handle
*** @retval NULL  An error occured
***/
extern gom_ContourDataChangedListener* gomp_AddContourDataChangedListener(
    GOM_ARG( gom_ContourDataChangedListenerFunc, callback    ),
    GOM_ARG( void                             *, callbackData ) );

/** @brief Unregister an contour data changed listener callback by using the handle.
***
*** @param listener  Listener handle
*** @return  Number of canceled listeners (always 1).
***/
extern int gomp_CancelContourDataChangedListener(
    GOM_ARG( gom_ContourDataChangedListener *, listener ) );

/** @brief Unregister an contour data changed listener callback by using the callback.
***
*** @param callback  Listener callback
*** @return  Number of canceled listeners (always non-negative).
***/
extern int gomp_CancelContourDataChangedListenersByFunc(
    GOM_ARG( gom_ContourDataChangedListenerFunc, callback ) );

/** @}
***/

/***** PUBLIC GOMAPI END *****/

#endif


extern int          gomp_GetContourData(const char *, float **,
                                        int *, int *, int *);
extern int          gomp_ContourDriver(const char *, const char *);
extern int          gomp_ParseContourLevels(int,
                                            float, float, float, float, int);
extern int          gomp_ParseContourLevelsProjection(int, float, float, float,
                                                      int);
extern int          gomp_FinalizeContourLevels(int, int);
extern int          gomp_ContourSmoothON(int, int);
extern int          gomp_ContourSmoothOFF(int, int);
extern int          gomp_ContourDisplayON(int, int);
extern int          gomp_ContourDisplayOFF(int, int);
extern int          gomp_SetContourDisplayType(int, int, int);
extern int          gomp_GetContourDisplayType(int, int);
extern int          gomp_SetContourDisplayTypeGlobal(int);
extern int          gomp_GetContourDisplayTypeGlobal(void);
extern const char *gomp_GetContourFileName(int);
extern int          gomp_GetContoursDefined(void);
extern float        gomp_GetContourMin(int);
extern float        gomp_GetContourMax(int);
extern int          gomp_CheckContourName(const char *);
extern int          gomp_ParseContourName(const char *);
extern int          gomp_GetContourLevels(int);
extern int          gomp_SetContourAlpha(int, int, float);
extern float        gomp_GetContourAlpha(int, int);
extern int          gomp_SetContourCullFace(int, int, int);
extern int          gomp_GetContourCullFace(int, int);
extern int          gomp_ShowContourInfo(void);
extern int          gomp_CheckContourIndex(int);
extern float        gomp_GetContourValue(int, int);
extern float        gomp_GetContourColourRed(int, int);
extern float        gomp_GetContourColourGreen(int, int);
extern float        gomp_GetContourColourBlue(int, int);
extern int          gomp_GetContourSmooth(int, int);
extern int          gomp_GetContourDisplayState(int, int);
extern const char *gomp_GetContourName(int);
extern int          gomp_GetContourCube(int,
                                        float *, float *,
                                        float *, float *,
                                        float *, float *,
                                        int *, int *, int *);

extern int          gomp_SetCutPlaneDamping(float);
extern float        gomp_GetCutPlaneDamping(void);

extern int          gomp_GetCutPlaneType_3D_X(void);
extern int          gomp_GetCutPlaneType_3D_Y(void);
extern int          gomp_GetCutPlaneType_3D_Z(void);

extern int          gomp_GetCutplaneSmoothTypeX(void);
extern int          gomp_GetCutplaneSmoothTypeY(void);
extern int          gomp_GetCutplaneSmoothTypeZ(void);

extern int          gomp_SetCutPlaneType_3D_X(int);
extern int          gomp_SetCutPlaneType_3D_Y(int);
extern int          gomp_SetCutPlaneType_3D_Z(int);

extern int          gomp_PrepareCutPlaneX(int, float, const char *,
                                          const char *, const char *);
extern int          gomp_PrepareCutPlaneY(int, float, const char *,
                                          const char *, const char *);
extern int          gomp_PrepareCutPlaneZ(int, float, const char *,
                                          const char *, const char *);
extern int          gomp_DeleteCutPlaneDataX(void);
extern int          gomp_DeleteCutPlaneDataY(void);
extern int          gomp_DeleteCutPlaneDataZ(void);
extern int          gomp_CheckContourName(const char *);
extern int          gomp_PrepareCutPlaneXYZ(int, int, const char *,
                                            const char *, const char *, 
                                            const char *, const char *,
                                            const char *, 
                                            const char *, const char *,
                                            const char *, 
                                            const char *, const char *);
extern int          gomp_GetCutPlaneXYZ(int, int *, 
                                        float *, float *, float *, 
                                        float *, float *, float *);
extern int          gomp_GetCutPlanePlotStateXYZ(int);
extern int          gomp_DisableCutPlanePlotStateXYZ(int);

extern int          gomp_PrepareCutPlaneSpectrumZ(int, const char *,
                                                  const char *, const char *,
                                                  const char *);
extern int          gomp_PlotCutPlaneSpectrumZ(void);

extern int          gomp_PrepareCutPlaneSpectrumY(int, const char *,
                                                  const char *, const char *,
                                                  const char *);
extern int          gomp_PlotCutPlaneSpectrumY(void);

extern int          gomp_PrepareCutPlaneSpectrumX(int,
                                                  const char *, const char *,
                                                  const char *, const char *);
extern int          gomp_PlotCutPlaneSpectrumX(void);

extern int          gomp_CalculateJPRisosurface1(int, int, float, float,
                                                 int, int, float, float,
                                                 int, polygon_t **);
extern int          gomp_CalculateJPRisosurface2(int, int, float, float, float,
                                                 float, float, float, int,
                                                 float, float);

extern int          gomp_SetContourLineWidth(int);
extern int          gomp_GetContourLineWidth(void);

extern int          gomp_SetContourProjection(int, int);
extern int          gomp_GetContourProjection(int);

extern float        gomp_GetContourProjectionMin(int, int);
extern float        gomp_GetContourProjectionMax(int, int);

extern int          gomp_ShowNumberOfPolygons(int, int);
extern int          gomp_SetSurfaceControlON(void);
extern int          gomp_SetSurfaceControlOFF(void);
extern int          gomp_GetSurfaceMethod(void);
extern int          gomp_SetSurfaceMethod(int);
extern int          gomp_DeleteInvalidatedPolygonData(void);
extern int          gomp_ReturnSurfacePolygonData(int, int, int);
extern int          gomp_SetContour2StructureMapping(int, int);
extern int          gomp_GetContour2StructureMapping(int);
extern int          gomp_DeleteGetContour2StructureMapping(void);
extern int          gomp_AddContour2StructureMappingSpace(void);

extern int          gomp_SetContourClippingPlaneState(int, int);
extern int          gomp_GetContourXClippingPlaneState(void);
extern int          gomp_GetContourYClippingPlaneState(void);
extern int          gomp_GetContourZClippingPlaneState(void);
extern int          gomp_SetContourXYZClippingPlaneState(int, int);
extern int          gomp_GetContourXYZClippingPlaneState(int);

extern int          gomp_SetContourClippingPlaneParameters( int, char, float);
extern int          gomp_SetContourXYZClippingPlaneParameters(
    int, float, float, float, float, float, float, float, float, float );
extern int          gomp_SetContourLevelClippingPlanePosition(int, int, float);
extern int          gomp_SetContourLevelClippingPlaneAxis(int, int, char);
extern float        gomp_GetContourLevelClippingPlanePosition(int, int);
extern char         gomp_GetContourLevelClippingPlaneAxis(int, int);

extern int          gomp_ReadContourInfo2ModelFile(FILE *);
extern int          gomp_WriteContourInfo2ModelFile(FILE *);
extern int          gomp_ReadContourClipplaneInfoFromModelFile(FILE *);
extern int          gomp_WriteContourClipplaneInfo2ModelFile(FILE *);

extern int          gomp_WriteCutPlane2ModelFile(FILE *);
extern int          gomp_ReadCutPlaneXFromModelFile(FILE *);
extern int          gomp_ReadCutPlaneYFromModelFile(FILE *);
extern int          gomp_ReadCutPlaneZFromModelFile(FILE *);
extern int          gomp_ReadCutPlaneXYZFromModelFile(FILE *);

extern int          gomp_DeleteAllContours(void);

extern int          gomp_PlotIsoSurf(void *,int,int);
extern int          gomp_InvalidateContourPlotter(unsigned long int);
extern int          gomp_ContourIsChanging(int);

extern int          gomp_CalculateIsoCurves(int, float, char, float);

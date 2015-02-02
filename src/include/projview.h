#ifdef ENABLE_EXTENSIONS
#  include "gomlib/gomprojview.h"
#  include "gomext/projview.h"
#endif

/* Projection transformations */
#define ORTHOGRAPHIC_VIEW     1
#define PERSPECTIVE_VIEW      0

#define GRAPHICS_NOT_AVAILABLE  0
#define GRAPHICS_AVAILABLE      1
#define GRAPHICS_STATUS_UNKNOWN 2

extern int          gomp_GetTermType(void);
extern int          gomp_SetTermType(int);

extern int          gomp_GetUpdateDisplayMode(void);
extern int          gomp_SetUpdateDisplayMode(int);

extern int          gomp_SetSystemRedisplayMode(int);
extern int          gomp_GetSystemRedisplayMode(void);

extern int          gomp_SetDrawBuffer(int);
extern int          gomp_GetDrawBuffer(void);

extern int          gomp_GetAllowIndividualScaling(void);
extern int          gomp_SetAllowIndividualScaling(int);

extern int          gomp_SetMouseButtonState(int);
extern int          gomp_GetMouseButtonState(void);

extern int          gomp_Rotate(float, float, float, float);
extern int          gomp_Translate(float, float, float);

extern int          gomp_RotateCoordinates1X(float, char);
extern int          gomp_TranslateCoordinatesX(float, float, float);

extern int          gomp_GetRotationState(void);
extern int          gomp_GetTranslationState(void);
extern int          gomp_SetRotationState(int);
extern int          gomp_SetTranslationState(int);

extern int          gomp_SaveTranslateArray(float, float, float);
extern int          gomp_SaveProjectionMatrix(const float *);
extern int          gomp_SaveModelViewMatrix(const float *);

/***** PUBLIC GOMAPI BEGIN *****/

/** @weakgroup gom_doc_plot
***/

/** @file
*** @ingroup gom_doc_plot
***/

/** @brief Get a translate array.
*** @ingroup gom_doc_plot
*** @return [x y z] translate vector.
***/
extern const float *gomp_GetTranslateArray(void);

/****** PUBLIC GOMAPI END ******/

extern const float *gomp_GetSavedProjectionMatrix(void);
extern const float *gomp_GetSavedModelViewMatrix(void);
extern void         gomp_ScaleModelView(float, float, float);

extern int          gomp_SaveTranslateArrayMT(int, float, float, float);
extern int          gomp_SaveProjectionMatrixMT(int, const float *);
extern int          gomp_SaveModelViewMatrixMT(int, const float *);
extern const float *gomp_GetTranslateArrayMT(int);
extern const float *gomp_GetSavedProjectionMatrixMT(int);
extern const float *gomp_GetSavedModelViewMatrixMT(int);    
extern void         gomp_ScaleModelViewMT(int, float, float, float);
extern int          gomp_AppendSpaceModelViewMatrixMT(size_t);
extern int          gomp_RemoveSpaceModelViewMatrixMT(int,size_t);

extern int          gomp_GetModelviewMatrix(float *);
extern int          gomp_PutModelviewMatrix(const float *);
extern int          gomp_GetProjectionMatrix(float *);
extern int          gomp_PutProjectionMatrix(const float *);

extern int          gomp_GetObjectCenterType(void);
extern int          gomp_SetObjectCenterType(int);

extern int          gomp_SetSizeOfSystem(float);
extern float        gomp_GetSizeOfSystem(void);

const char         *gomp_SetGlobalTransformationState(void);
const char         *gomp_SetLocalTransformationState(void);

extern int          gomp_GetSystemTranslateState(void);
extern int          gomp_SetSystemTranslateState(int);

extern float        gomp_GetTranslationDamping(void);
extern int          gomp_SetTranslationDamping(float);

extern int          gomp_SetPerspectiveNear(float);
extern int          gomp_SetPerspectiveNearStep(float);
extern float        gomp_GetPerspectiveNearStep(void);

extern int          gomp_SetPerspectiveFar(float);
extern int          gomp_SetPerspectiveFarStep(float);
extern float        gomp_GetPerspectiveFarStep(void);

extern int          gomp_SetPerspectiveStep(float);
extern float        gomp_GetPerspectiveNear(void);
extern float        gomp_GetPerspectiveFar(void);
extern int          gomp_SetPerspectiveAngle(float);
extern float        gomp_GetPerspectiveAngle(void);
extern float        gomp_GetPerspectiveStep(void);
extern int          gomp_SetPerspectiveWindow(float);
extern int          gomp_SetPerspectiveWindow(float);
extern int          gomp_ResetPerspectiveWindowAttributes(void);
extern float        gomp_GetPerspectiveWindow(void);
extern float        gomp_GetPerspectiveWindowStep(void);
extern int          gomp_SetPerspectiveWindowStep(float);
extern float        gomp_GetPerspectiveDistanceStep(void);
extern int          gomp_SetPerspectiveDistanceStep(float);
extern int          gomp_SetPerspectiveWindowAttributes(float , float);

extern int          gomp_SetProjectionTransformation(int);
extern int          gomp_GetProjectionTransformation(void);

extern int          gomp_DoViewingTransformationOverStructures(int);

extern int          gomp_PrepareDisplay(int);
extern int          gomp_PrepareStatusDisplay(const char *);

extern int          gomp_ResetProjection(void);

extern int          gomp_WriteDisplayAttributesGOM(FILE *);
extern int          gomp_ReadDisplayAttributesGOM(FILE *);

extern int          gomp_ScaleDisplay(float , float , float);
extern int          gomp_ResetProjection(void);

extern const float *gomp_CalculateGeometricalCenter(int);
extern const float *gomp_CalculateGeometricalCenterMT(int);

extern int          gomp_PushModelViewingData(void);
extern int          gomp_PopModelViewingData(void);

extern int          gomp_IdentifyAtomFromCoords(int , int , int);

extern int          gomp_ResetView(void);

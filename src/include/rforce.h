extern int        gomp_DeleteVectorStructure(void);
extern int        gomp_SetVectorListLength(int);
extern int        gomp_GetVectorListLength(void);
extern const int *gomp_GetVectorListArray(void);
extern int       *gomp_GetModifiableVectorListArray(void);
extern int        gomp_GetVectorStructureIndex(void);

extern int        gomp_SetVectorScale(float);
extern float      gomp_GetVectorScale(void);
extern int        gomp_SetVectorRadius(float);
extern float      gomp_GetVectorRadius(void);

extern int        gomp_SetPlotVectorStatus(int);
extern int        gomp_GetPlotVectorStatus(void);
extern int        gomp_SetVectorDisplayRange(float  , float);

extern int        gomp_ReadCharmmVector(const char *);
extern int        gomp_ReadFlatFileVector(const char *);

extern int        gomp_GetVectorDisplayRange(float *, float *);

extern const float *gomp_GetVectorForceXp(void);
extern const float *gomp_GetVectorForceYp(void);
extern const float *gomp_GetVectorForceZp(void);

extern int        gomp_PlotVectorData(void*,int,int);

extern int        gomp_ParsePlotVectorList(const char *,const char *,
                                           const char *, const char *,
                                           const char *);


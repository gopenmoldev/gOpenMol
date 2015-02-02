extern int gomp_DelLineSeg(void);

extern int gomp_DelArrowSeg(void);
extern int gomp_WriteArrowSeg2ModelFile(FILE *);
extern int gomp_ReadArrowSegFromModelFile(FILE *);

extern int gomp_DelSphereSeg(void);
extern int gomp_WriteSphereSeg2ModelFile(FILE *);
extern int gomp_ReadSphereSegFromModelFile(FILE *);

extern int gomp_DelCylinderSeg(void);
extern int gomp_WriteCylinderSeg2ModelFile(FILE *);
extern int gomp_ReadCylinderSegFromModelFile(FILE *);

extern int gomp_DelPlaneSeg(void);
extern int gomp_WritePlaneSeg2ModelFile(FILE *);
extern int gomp_ReadPlaneSegFromModelFile(FILE *);

extern int gomp_DelTriangleSeg(void);
extern int gomp_WriteTriangleSeg2ModelFile(FILE *);
extern int gomp_ReadTriangleSegFromModelFile(FILE *);

extern int gomp_GetTotalLineArrowEntries(void);
extern int gomp_DelLineArrowSeg(void);
extern int gomp_WriteLineArrowSeg2ModelFile(FILE *);
extern int gomp_ReadLineArrowSegFromModelFile(FILE *);

extern int gomp_SetGradientDisplayStyle(int);
extern int gomp_SetColourScalePlotStatus(int);

extern int   gomp_SetColourScalePlotMin(float);
extern float gomp_GetColourScalePlotMin(void);
extern int   gomp_SetColourScalePlotMax(float);
extern float gomp_GetColourScalePlotMax(void);
extern int   gomp_SetColourScalePlotStep(float);
extern float gomp_GetColourScalePlotStep(void);
extern int   gomp_SetColourScalePlotLevels(int);
extern int   gomp_GetColourScalePlotLevels(void);
extern int   gomp_SetColourScalePlotStatus(int);
extern int   gomp_GetColourScalePlotStatus(void);

extern int gomp_PushArrowStack(const char *,const char *,const char *,
                          const char *,const char *,const char *,
                          const char *,const char *,const char *,
                          const char *);
extern int gomp_PushSphereStack(const char *, const char *,const char *,
                           const char *, const char *,const char *,
                           const char *, const char *,const char *, 
                           const char *);
extern int gomp_PushCylinderStack(const char *,const char *,const char *,
                             const char *,const char *,const char *,
                             const char *,const char *,const char *,
                             const char *);
extern int gomp_PushLineStack(const char *,const char *,const char *,
                         const char *,const char *,const char *,
                         const char *,const char *);
extern int gomp_PushPlaneStack(const char *,const char *,const char *,
                          const char *,const char *,const char *,
                          const char *);
extern int gomp_PushTriangleStack(const char *, const char *,const char *,
                             const char *, const char *,const char *,
                             const char *, const char *,const char *,
                             const char *, const char *,const char *,
                             const char *, const char *,const char *,
                             const char *, const char *,const char *,
                             const char *, const char *,const char *,
                             const char *, const char *);
extern int gomp_PushLineArrowStack(float,float,float,
                              float,float,float,
                              const char *,const char *);

extern int gomp_DisplayColourScale(void);

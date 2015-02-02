extern int gomp_CalcRDF(const int *,int ,float ,float ,float ,float,int,
                      const char *,const char *,const char *);
extern int gomp_ParseRDFList(const char *, const char *, const char *,
                           const char *, const char *, const char *,
                           const char *, const char *);
extern int gomp_GetRDFobservations(void);
extern int gomp_GetNumRDFObs(void);
extern int gomp_CalcMeanRDF(void);
extern int gomp_DeleteRDF(void);
extern int gomp_GetNumRDFObs(void);
extern const float *gomp_GetRDFVecX(void);
extern const float *gomp_GetRDFVecY(void);
extern int gomp_WriteRDF(const char *);
extern int gomp_GetRDFstatus(void);
extern int gomp_GetRDFplotStatus(void);
extern int gomp_GetRDFobservations(void);
extern int gomp_GetRDFaverageStatus(void);

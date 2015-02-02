extern int gomp_PreCluster(void);
extern int gomp_PlotClusterMatrix(int , int ,
                                  float, float ,float ,float ,
                                  float, float ,float ,float ,float);
extern int gomp_DeleteClusterData(void);

extern int gomp_CalcCluster(
    const char *, const char *,
    const char *, const char *,
    const char *, const char *);

extern int gomp_ReadClusterData(const char *);
extern int gomp_WriteClusterData(const char * , const char * );
extern int gomp_GetDisplayCLUSTERmatrix(void);
extern int gomp_SetDisplayCLUSTERmatrix(int);
extern int gomp_GetClusterStatus(void);
extern const float *gomp_GetClusterData(void);

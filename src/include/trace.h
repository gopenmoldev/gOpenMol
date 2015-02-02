/* type of export input types */
#define CHARMM_INPUT      1
#define MOPAC_INPUT       2
#define GAMESS_INPUT      3
#define GAUSSIANXX_INPUT  4
#define OPENMOL_INPUT     5
#define PROBESURF_INPUT   6
#define YASP_INPUT        7
#define ICON8_INPUT       8

extern int gomp_PlotAtomTrace(void*, int, int);
extern int gomp_DeleteTrace(void);
extern int gomp_WriteAtomTrace(int , FILE *);
extern int gomp_AtomsInTrace(void);

extern int gomp_GetDisplayTraceAtoms(void);
extern int gomp_SetDisplayTraceAtoms(int);

extern int gomp_GetTraceSets(void);

extern int gomp_TraceAtoms(const char *,const char *, const char *, int);

extern int gomp_GetTraceState(void);

extern const int *gomp_GetTraceAtomsInSet(void);
extern const int *gomp_GetTraceAtomList(void);
extern const float *gomp_GetTraceAtomXCoord(void);
extern const float *gomp_GetTraceAtomYCoord(void);
extern const float *gomp_GetTraceAtomZCoord(void);

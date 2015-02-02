/* plumber types */
#define RIBBON_TYPE        0
#define CYLINDER_TYPE      1
#define FLAT_HELIX_TYPE    2
#define SOLID_HELIX_TYPE   3
#define ARROW_TYPE         4
#define STRAND_TYPE        5
#define TRACE_TYPE         6

typedef struct {
    int     Wstr;
    int     IsUpdated;
    int     Atoms;
    int    *atom_list;
    int     Type;
    int     Points;
    float   red;
    float   green;
    float   blue;
    double *width;
    double *thickness;
    double *x;
    double *y;
    double *z;
    double *NX;
    double *NY;
    double *NZ;
    double *WX;
    double *WY;
    double *WZ;
} PlumberDataType;

extern int           gomp_DeletePlumber(int);
extern int           gomp_DeletePlumbers(void);
extern int           gomp_RetrievePlumberInfoFromModelFile(FILE *);
extern int           gomp_StorePlumberInfo2ModelFile(FILE *);
extern int           gomp_GetPlumberSets(void);
extern const double *gomp_GetPlumberXp(int);
extern const double  *gomp_GetPlumberYp(int);
extern const double  *gomp_GetPlumberZp(int);
extern float         gomp_GetPlumberRed(int);
extern float         gomp_GetPlumberGreen(int);
extern float         gomp_GetPlumberBlue(int);
extern int           gomp_GetPlumberStructure(int);
extern int           gomp_GetPlumberAtoms(int);
extern const int    *gomp_GetPlumberAtomList(int);
extern int           gomp_SetPlumberRed(int , float);
extern int           gomp_SetPlumberGreen(int , float);
extern int           gomp_SetPlumberBlue(int , float);
extern int           gomp_SetPlumberStructure(int,int);
extern int           gomp_SetPlumberDisplay(int);
extern int           gomp_GetPlumberDisplay(void);
extern int           gomp_PlotPlumber(void*,int,int);
extern int           gomp_LoadPlumberAtoms(int, int, const int *, float, float, float, float, int , float , float , int );
extern int           gomp_SetPlumberDisplayType(int , int);
extern int           gomp_GetPlumberDisplayType(int);
extern int           gomp_ParsePlumberList(const char *,
                                         const char *,
                                         const char *,
                                         float, float, float,
                                         float, int,
                                         float, float, int);
extern const PlumberDataType *gomp_GetPlumbers(void);

extern int gomp_LinearSpline(double  *, int, const int *);
extern int gomp_CubicSpline( double  *, int, const int *);

extern int gomp_ProteinChainList(int, int **, int **, int *);
extern int gomp_ParseListProteinChains(int, const char * , const char *);
extern int gomp_ParseListAtomsAround(int,
                                   const char *, const char *, const char *,
                                   const char *, const char *, const char *,
                                   const char *,int);
extern int gomp_ParseListAtomsJoin(int, int, const char **,int);
extern int gomp_SecondaryStructure(void);

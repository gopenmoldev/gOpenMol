/*

  Information about the monitor and time series stuff
  Leif Laaksonen, December 1996
*/

#define DIST_LIST   0
#define ANGLE_LIST  1

#define MAX_FLOAT  1.e+35    /* This max float is just a guess */

/* manipulation list */
#define     DAVERAGE      1
#define     SQUARE        2
#define     COS           3
#define     COS2          4
#define     SQRT          5
#define     DINITIAL      6
#define     COPYNR        7
#define     ADDNR         8
#define     LOG           9
#define     EXP          10
#define     POWERREAL    11
#define     MULTREAL     12
#define     DIVIDEREAL   13
#define     SHIFTREAL    14
#define     DMIN         15
#define     ABS          16
#define     DIVFIRST     17
#define     DIVMAXIMUM   18
#define     PSPECTRUM    19
#define     ZERO         20

#define DIST_TYPE  0
#define ANG_TYPE   1
#define TORS_TYPE  2

#define   SEGMENT      MAX_SEG_NAME_LEN
#define   RESIDUE      MAX_RES_NAME_LEN
#define   ATOM         MAX_ATM_NAME_LEN
#define   SEP_STRUCT   '^'


extern    int    gomp_ExpandDistVec(int);
extern    int    gomp_SelectDistArray(
    const char *Seg1, const char *Res1, const char *Atm1,
    const char *Seg2, const char *Res2, const char *Atm2,
    const char *Type, const char *Colour);
extern    int    gomp_SetDistMonitor(int);
extern    int    gomp_GetDistMonitor(void);
extern    int    gomp_GetDistMonitorSamples(void);
extern    const int *gomp_GetDistMonitorSamplesList(void);
extern    int    gomp_GetDistMonitorType(int);
extern    int    gomp_GetDistMonitorColor(int , float *, float *, float *);
extern    int    gomp_SetDistMonitorType(int , int);
extern    int    gomp_SetDistMonitorColor(int , float, float, float);
extern    int    gomp_FillDistanceSeries(void);

/* angle monitor array */
extern    int    gomp_ExpandAngVec(int);
extern    int    gomp_SelectAngArray(
    const char *Seg1, const char *Res1, const char *Atm1,
    const char *Seg2, const char *Res2, const char *Atm2,
    const char *Seg3, const char *Res3, const char *Atm3,
    const char *Type, const char *Color);
extern    int    gomp_SetAngMonitor(int);
extern    int    gomp_GetAngMonitor(void);
extern    int    gomp_GetAngMonitorSamples(void);
extern    const int *gomp_GetAngMonitorSamplesList(void);
extern    int    gomp_GetAngMonitorType(int);
extern    int    gomp_GetAngMonitorColor(int , float *, float *, float *);
extern    int    gomp_SetAngMonitorType(int , int);
extern    int    gomp_SetAngMonitorColor(int , float, float, float);
extern    int    gomp_FillAngleSeries(void);

/* torsion monitor array */
extern    int    gomp_ExpandTorsVec(int);
extern    int    gomp_SelectTorsArray(
    const char *Seg1, const char *Res1, const char *Atm1,
    const char *Seg2, const char *Res2, const char *Atm2,
    const char *Seg3, const char *Res3, const char *Atm3,
    const char *Seg4, const char *Res4, const char *Atm4,
    const char *Type, const char *Color);
extern    int    gomp_SetTorsMonitor(int);
extern    int    gomp_GetTorsMonitor(void);
extern    int    gomp_GetTorsMonitorSamples(void);
extern    const int *gomp_GetTorsMonitorSamplesList(void);
extern    int    gomp_GetTorsMonitorType(int);
extern    int    gomp_GetTorsMonitorColor(int , float *, float *, float *);
extern    int    gomp_SetTorsMonitorType(int , int);
extern    int    gomp_SetTorsMonitorColor(int , float, float, float);
extern    int    gomp_FillTorsionSeries(void);


extern    int    gomp_ResetMonitorData(void);
extern    int    gomp_ResetMonitorDistanceData(void);
extern    int    gomp_ResetMonitorAngleData(void);
extern    int    gomp_ResetMonitorTorsionData(void);

extern    int    gomp_CheckTimeSeriesIndex(int , int);
extern    int    gomp_ManipulateTimeSeries(int , int , const char *, const char * , const char *);
extern    int    gomp_Mantime(float *, float *, int , int , float , int );

/* time series handle */
typedef struct _TimeSeries {
    int    Length;
    float *Value;
} TimeSeries;

/* distance */
extern    int         gomp_ExpandeDistanceSeries(int);
extern    int         gomp_DeleteDistanceSeries(void);
extern    int         gomp_StoreValues2DistanceSeries(int , int , const float *);
extern    int         gomp_FillDistanceSeries(void);
extern    int         gomp_WriteDistanceVector(int,const char *);
extern    int         gomp_GetNumberDistanceSeries(void);
extern    int         gomp_EditDistVector( int , int , int , int , float, float , float);
extern    int         gomp_DeleteDistanceIndex(int);


/* angle */
extern    int         gomp_ExpandeDistanceSeries(int);
extern    int         gomp_DeleteAngleSeries(void);
extern    int         gomp_FillAngleSeries(void);
extern    int         gomp_WriteAngleVector(int,const char *);
extern    int         gomp_GetNumberAngleSeries(void);
extern    int         gomp_EditAngVector( int , int , int , int , int , float, float , float);
extern    int         gomp_DeleteAngleIndex(int);

/* torsion */
extern    int         gomp_ExpandeDistanceSeries(int);
extern    int         gomp_DeleteTorsionSeries(void);
extern    int         gomp_FillTorsionSeries(void);
extern    int         gomp_WriteTorsionVector(int,const char *);
extern    int         gomp_GetNumberTorsionSeries(void);
extern    int         gomp_EditTorsVector( int , int , int , int , int , int , float, float , float);
extern    int         gomp_DeleteTorsionIndex(int);

/* correlation data structure */
typedef struct {
      int    corr_wind;   
             /* switch to indicate that window is on = 1 , off = 0 */
      float *corr_vec1;   /* pointer to time series 1 */
      float *corr_vec2;   /* pointer to time series 2 */
      int    corr_obs;    /* number of points in vectors 1 and 2 */
      int    corr_length; /* length of the correlation array */
      float *corr_val;    /* pointer to the correlation values */
} corr_info_t;

extern      int    gomp_PutCorrelationVectorAddress1(float *);
extern      int    gomp_PutCorrelationVectorAddress2(float *);
extern      const float *gomp_GetCorrelationResultVector(void);
extern      int    gomp_CalculateCorrelation(int , int , int);
extern      int    gomp_PreCorrel(void);
extern      int    gomp_GetCorrelationResultVectorLength(void);

extern      int    gomp_WriteBackboneDihedrals(const char *);

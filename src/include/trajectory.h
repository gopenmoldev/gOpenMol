/*
    Leif Laaksonen, Center for Scientific Computing 1995 - 2003
*/

#define TRAJ_OLD  0   /* Coordinate data is loaded into first (0) position  */
#define TRAJ_NEW  1   /* Coordinate data is loaded into first free position */
#define TRAJ_TMP  2   /* Coordinate data is loaded into temporary position  */

/*  trajectory functions and handles */

#define MAXplen 300      /* max number of points */
#define MAXdlen 600      /* max number of distances = MAXdlen/2 */
#define MAXalen 900      /* max number of angles    = MAXalen/3 */
#define MAXtlen 1200     /* max number of torsions  = MAXtlen/4 */
#define MAXcclen 900     /* max number of coordinate centres is 300 */
#define MAXmclen 900     /* max number of mass centres is 300       */

/* trajectory file type accepted */
#define  AMBER_TRAJ      1
#define  CERIUS2_TRAJ    2
#define  CHARMM_TRAJ     3
#define  DISCOVER_TRAJ   4
#define  FDL_POLY_TRAJ   5
#define  UDL_POLY_TRAJ   6
#define  FAMBER_TRAJ     7
#define  GROMACS_TRAJ    8
#define  GROMOS_TRAJ     9
#define  GROMOS96A_TRAJ  10
#define  GROMOS96B_TRAJ  11
#define  HYPERCHEM_TRAJ  12
#define  MUMOD_TRAJ      13
#define  OPENMOL_TRAJ    14
#define  TINKER_TRAJ     15
#define  XMOL_TRAJ       16
#define  XPLOR_TRAJ      17
#define  YASP_TRAJ       18

#define  TRAJ_PLAY            1
#define  TRAJ_FORWARD_FRAME   2
#define  TRAJ_BACKWARD_FRAME  3
#define  TRAJ_STOP_LOOP       4
#define  TRAJ_FIRST_FRAME     5
#define  TRAJ_LAST_FRAME      6

/* structure for the dynamics trajectory file */

typedef struct {
    int inuse;                                /* in use/not in use  */
    char traj_file[BUFF_LEN];                 /* file name          */
    int display;                              /* display frame number */
    int type;                                 /* trajectory type    */
    int natom;                                /* number of atoms    */
    int nstep;                                /* number of steps    */
    int time_bw_steps;                        /* time between steps */
    int time_first_frame;                     /* time of first frame */
    int first_frame;                          /* first frame to be displayed */
    int last_frame;                           /* last frame to be displayed  */
    int delta_frame;                          /* display every delta frame */
    int retrieve_type;                        /* fast/slow formatted reader */
    int retrieve_velocity;                    /* switch to control the retrieve of
                                                 the velocities. 
                                                 Remember that not all trajectories
                                                 have this information included. */
    int retrieve_velocity_ready;              /* velocities are ready */
    float *Xvelocities;                       /* velocities array */
    float *Yvelocities;                       /* velocities array */
    float *Zvelocities;                       /* velocities array */
    float  Vmin;                              /* velocity min */
    float  Vmax;                              /* velocity max */
    int retrieve_force;                       /* switch to control the retrieve of
                                                 the forces.
                                                 Remember that not all trajectories
                                                 (in fact very few) have this information
                                                 included */
    int retrieve_force_ready;                 /* forces are ready */
    float *Xforces;                           /* forces array     */
    float *Yforces;                           /* forces array     */
    float *Zforces;                           /* forces array     */
    float  Fmin;                              /* force  min */
    float  Fmax;                              /* force  max */
} TrajectoryInfo_t;

extern int         gomp_ParseTrajType(const char *);
extern int         gomp_SetTrajectoryStructureInUse(void);
extern int         gomp_GetTrajectoryStructureState(void);
extern int         gomp_SetTrajectoryStructureOff(void);
extern int         gomp_DeleteTrajectoryStructure(void);
extern int         gomp_SetTrajectoryFileName(const char *);
extern const char *gomp_GetTrajectoryFileName(void);
extern int         gomp_SetTrajectoryFileType(int);
extern int         gomp_GetTrajectoryFileType(void);
extern const char *gomp_GetTrajectoryFileTypeName(void);
extern int         gomp_SetNumberOfFrames(int);
extern int         gomp_GetNumberOfFrames(void);
extern int         gomp_SetTrajectoryTimeInfo(int  , int);
extern int         gomp_GetTrajectoryTimeInfo(int *, int *);
extern int         gomp_SetTrajectoryDisplayParams(int  , int  , int);
extern int         gomp_GetTrajectoryDisplayParams(int *, int *, int *);
extern int         gomp_GetTrajectoryFirstFrame(void);
extern int         gomp_GetTrajectoryLastFrame(void);
extern int         gomp_GetTrajectoryDisplayFrames(void);

extern int         gomp_GetTrajectoryDeltaFrame(void);
extern int         gomp_SetNumberOfTrajectoryAtoms(int);
extern int         gomp_GetNumberOfTrajectoryAtoms(void);

extern int         gomp_GetNumberOfFreeAtoms(int);
extern int         gomp_SetNumberOfFreeAtoms(int , int);
extern const int  *gomp_GetFreeAtomListPointer(int);
extern int        *gomp_GetModifiableFreeAtomListPointer(int);
extern int         gomp_SetFreeAtomListPointer(int , int);

extern int         gomp_GetDisplayFrameNumber(void);
extern int         gomp_PutDisplayFrameNumber(int);
extern int         gomp_PutRunningFrameNumber(int);
extern int         gomp_DisplayRunningFrameNumber(void);
extern int         gomp_SetDisplayRunningFrameNumberState(int);
extern int         gomp_GetDisplayRunningFrameNumberState(void);

extern int         gomp_GetOneFrame(int , FILE* , int);

extern int         gomp_SetFormattedTrajectoryReader(int);
extern int         gomp_GetFormattedTrajectoryReader(void);

extern int         gomp_PeekRunningFrameNumberProperty(void);

extern const char *gomp_GetTrajectoryTypeFileExtension(int);
extern int         gomp_PutTrajectoryTypeFileExtension(int, const char *);
extern int         gomp_GetTrajectoryTypeDefault(void);
extern int         gomp_PutTrajectoryTypeDefault(int);

extern int  gomp_GetFrameAmber(int , FILE* , int);
extern int  gomp_GetFrameCharmm(int , FILE* , int);
extern int  gomp_GetFrameCerius2(int , FILE* , int);
extern int  gomp_GetFrameDiscover(int , FILE* , int);
extern int  gomp_GetFrameDL_PolyFORMATTED(int , FILE* , int);
extern int  gomp_GetFrameDL_PolyUNFORMATTED(int , FILE* , int);
extern int  gomp_GetFrameFAmber(int , FILE* , int);
extern int  gomp_GetFrameGromos(int , FILE* , int);
extern int  gomp_GetFrameGromacs(int , FILE* , int);
extern int  gomp_GetFrameGROMOS96A(int , FILE* , int);
extern int  gomp_GetFrameGROMOS96B(int , FILE* , int);
extern int  gomp_GetFrameHyperChem(int , FILE* , int);
extern int  gomp_GetFrameMumod(int , FILE* , int);
extern int  gomp_GetFrameTINKER(int , FILE* , int);
extern int  gomp_GetFrameXmol(int , FILE* , int);
extern int  gomp_GetFrameXplor(int , FILE* , int);
extern int  gomp_GetFrameYasp(int , FILE* , int);

extern int          gomp_SetVelocityRetrieveState(int);
extern int          gomp_GetVelocityRetrieveState(void);
extern int          gomp_SetForceRetrieveState(int);
extern int          gomp_GetForceRetrieveState(void);
extern int          gomp_DeleteVelocityForceData(void);
extern int          gomp_SetVelocityRetrieveReadyState(int);
extern int          gomp_GetVelocityRetrieveReadyState(void);
extern int          gomp_GetVelocitySpace(int);
extern const float *gomp_GetVelocityXComponentPointer(void);
extern const float *gomp_GetVelocityYComponentPointer(void);
extern const float *gomp_GetVelocityZComponentPointer(void);
extern float       *gomp_GetModifiableVelocityXComponentPointer(void);
extern float       *gomp_GetModifiableVelocityYComponentPointer(void);
extern float       *gomp_GetModifiableVelocityZComponentPointer(void);
extern int          gomp_SetVelocityMinMax(float  , float);
extern int          gomp_GetVelocityMinMax(float *, float *);
extern int          gomp_SetForceRetrieveReadyState(int);
extern int          gomp_GetForceRetrieveReadyState(void);
extern int          gomp_GetForceSpace(int);
extern const float *gomp_GetForceXComponentPointer(void);
extern const float *gomp_GetForceYComponentPointer(void);
extern const float *gomp_GetForceZComponentPointer(void);
extern float       *gomp_GetModifiableForceXComponentPointer(void);
extern float       *gomp_GetModifiableForceYComponentPointer(void);
extern float       *gomp_GetModifiableForceZComponentPointer(void);
extern int          gomp_SetForceMinMax(float  , float);
extern int          gomp_GetForceMinMax(float *, float *);
extern int          gomp_CalculateVelocityMinMax(void);
extern int          gomp_CalculateForceMinMax(void);
extern int          gomp_TrajAvStructure(void);

extern int          gomp_GetTrajectory(const char *, const char *, int);
extern int          gomp_DisplayTrajectory(int);

/*  end of trajectory functions and handles */

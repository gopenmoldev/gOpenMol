/* type of coordinate files supported */

#define CHARMM_COORD     1
#define FREE_COORD       2
#define INSIGHT_COORD    3
#define AMBER_COORD      4
#define YASP_COORD       5
#define MUMOD_COORD      6
/*#define GROMOS_COORD     7*/
#define MOPAC_COORD      8
#define PDB_COORD        9
/*#define GAMESS_COORD    10*/
#define OPENMOL_COORD   11
#define HYPER_COORD     12
#define GXYZ_COORD      13
#define SOCKET_COORD    14
#define DYNAFRAME_COORD 15
#define GAUSSIAN_COORD  16
#define MOL2_COORD      17
#define XMOL_COORD      18
/*#define TXYZ_COORD      19*/
/*#define USER_COORD      20*/
#define GROMOS96_COORD  21
/*#define UHBD_COORD      22*/
#define PDBQ_COORD      23
/*#define ADF_COORD       24*/
/*#define GROMACS_COORD   25*/
/*#define DL_POLY_COORD   26*/
/*#define JAGUAR_COORD    27*/
/*#define TURBOMOLE_COORD 28*/

extern int         gomp_ResetAllDefaultFileTypes(void);

extern const char *gomp_GetAMBERcoordFileType(void);
extern int         gomp_PutAMBERcoordFileType(const char *);
extern int         gomp_GetAMBERdefault(void);
extern int         gomp_PutAMBERdefault(void);

extern const char *gomp_GetCHARMMcoordFileType(void);
extern int         gomp_PutCHARMMcoordFileType(const char *);
extern int         gomp_GetCHARMMdefault(void);
extern int         gomp_PutCHARMMdefault(void);

extern const char *gomp_GetGAMESScoordFileType(void);
extern int         gomp_PutGAMESScoordFileType(const char *);
extern int         gomp_GetGAMESSdefault(void);
extern int         gomp_PutGAMESSdefault(void);

extern const char *gomp_GetGAUSSIANcoordFileType(void);
extern int         gomp_PutGAUSSIANcoordFileType(const char *);
extern int         gomp_GetGAUSSIANdefault(void);
extern int         gomp_PutGAUSSIANdefault(void);

extern const char *gomp_GetGROMACScoordFileType(void);
extern int         gomp_PutGROMACScoordFileType(const char *);
extern int         gomp_GetGROMACSdefault(void);
extern int         gomp_PutGROMACSdefault(void);

extern const char *gomp_GetHYPERCHEMcoordFileType(void);
extern int         gomp_PutHYPERCHEMcoordFileType(const char *);
extern int         gomp_GetHYPERCHEMdefault(void);
extern int         gomp_PutHYPERCHEMdefault(void);

extern const char *gomp_GetINSIGHTcoordFileType(void);
extern int         gomp_PutINSIGHTcoordFileType(const char *);
extern int         gomp_GetINSIGHTdefault(void);
extern int         gomp_PutINSIGHTdefault(void);

extern const char *gomp_GetJAGUARcoordFileType(void);
extern int         gomp_PutJAGUARcoordFileType(const char *);
extern int         gomp_GetJAGUARdefault(void);
extern int         gomp_PutJAGUARdefault(void);

extern const char *gomp_GetMOL2coordFileType(void);
extern int         gomp_PutMOL2coordFileType(const char *);
extern int         gomp_GetMOL2default(void);
extern int         gomp_PutMOL2default(void);

extern const char *gomp_GetMOPACgraphcoordFileType(void);
extern int         gomp_PutMOPACgraphcoordFileType(const char *);
extern int         gomp_GetMOPACgraphdefault(void);
extern int         gomp_PutMOPACgraphdefault(void);

extern const char *gomp_GetMUMODcoordFileType(void);
extern int         gomp_PutMUMODcoordFileType(const char *);
extern int         gomp_GetMUMODdefault(void);
extern int         gomp_PutMUMODdefault(void);

extern const char *gomp_GetOPENMOLcoordFileType(void);
extern int         gomp_PutOPENMOLcoordFileType(const char *);
extern int         gomp_GetOPENMOLdefault(void);
extern int         gomp_PutOPENMOLdefault(void);

extern const char *gomp_GetPDBcoordFileType(void);
extern int         gomp_PutPDBcoordFileType(const char *);
extern int         gomp_GetPDBdefault(void);
extern int         gomp_PutPDBdefault(void);

extern const char *gomp_GetTINKERcoordFileType(void);
extern int         gomp_PutTINKERcoordFileType(const char *);
extern int         gomp_GetTINKERdefault(void);
extern int         gomp_PutTINKERdefault(void);

extern const char *gomp_GetUHBDcoordFileType(void);
extern int         gomp_PutUHBDcoordFileType(const char *);
extern int         gomp_GetUHBDdefault(void);
extern int         gomp_PutUHBDdefault(void);

extern const char *gomp_GetXMOLcoordFileType(void);
extern int         gomp_PutXMOLcoordFileType(const char *);
extern int         gomp_GetXMOLdefault(void);
extern int         gomp_PutXMOLdefault(void);

extern const char *gomp_GetXYZcoordFileType(void);
extern int         gomp_PutXYZcoordFileType(const char *);
extern int         gomp_GetXYZdefault(void);
extern int         gomp_PutXYZdefault(void);

extern const char *gomp_GetYASPcoordFileType(void);
extern int         gomp_PutYASPcoordFileType(const char *);
extern int         gomp_GetYASPdefault(void);
extern int         gomp_PutYASPdefault(void);

extern const char *gomp_GetOPENMOLcenterFileType(void);
extern int         gomp_PutOPENMOLcenterFileType(const char *);
extern int         gomp_GetOPENMOLcenterdefault(void);
extern int         gomp_PutOPENMOLcenterdefault(void);

extern int         gomp_SetGROMOS96CoordAmplifier(float);
extern float       gomp_GetGROMOS96CoordAmplifier(void);

extern int         gomp_ImportDictionaryAndApply(const char *);

extern int         gomp_ReadCoordinates(const char *,const char *,const char *);

extern int         gomp_ReadCoordinatesAMBER(const char * , int);
extern int         gomp_ReadCoordinatesCHARMM(const char * , int);
extern int         gomp_ReadCoordinatesFRAME(int, int);
extern int         gomp_ReadCoordinatesGAUSSIAN(const char * , int);
extern int         gomp_ReadCoordinatesHYPERCHEM(const char * , int);
extern int         gomp_ReadCoordinatesINSIGHT(const char * , int);
extern int         gomp_ReadCoordinatesMOL2(const char * , int);
extern int         gomp_ReadCoordinatesMOPACgraph(const char * , int);
extern int         gomp_ReadCoordinatesMUMOD(const char * , int);
extern int         gomp_ReadCoordinatesOPENMOL(const char * , int);
extern int         gomp_ReadCoordinatesPDB(const char * , int);
extern int         gomp_ReadCoordinatesXMOL(const char * , int);
extern int         gomp_ReadCoordinatesGXYZ(const char * , int);
extern int         gomp_ReadCoordinatesYASP(const char * , int);

extern int         gomp_ExternalInput4ICON8(int, const char *);
extern int         gomp_WriteCoordOPENMOL(int, const char *, int);
extern int         gomp_ExternalInput4PROBESURF(int, const char *,
                                                float, float,
                                                float, float,
                                                float, float,
                                                int,int,int,
                                                float);

extern int         gomp_WriteCoordBaS(int, const char *, int);
extern int         gomp_WriteCoordCHARMM(int, const char *, int);
extern int         gomp_WriteCoordFree(int, const char *, int);
extern int         gomp_WriteCoordPDB(int, const char *, int);
extern int         gomp_WriteCoordXYZ(int, const char *, int);

extern int         gomp_StartWithCoordinateFile(const char *);
extern int         gomp_StartWithModelFile(const char *);

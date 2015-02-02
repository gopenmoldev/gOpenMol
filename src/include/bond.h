extern int gomp_SetBondReconnectivityState(int);
extern int gomp_GetBondReconnectivityState(void);
extern int gomp_SetHBondReconnectivityState(int);
extern int gomp_GetHBondReconnectivityState(void);

extern int gomp_CheckBond(int, int, int, int);
extern int gomp_BreakBond(int, int, int, int);

extern int gomp_BreakBondM(int , const char *, const char *, const char *);

extern int gomp_EditBondI(int , int , int ,
                          const char *, const char *, const char *,
                          const char *, const char *,const char *);
extern int gomp_EditBondM(int , int ,
                          const char *, const char *, const char *,
                          const char *, const char *,const char *);


extern int gomp_CalcAtomConn(int);
extern int gomp_CalcAtomConnI(int, int);
extern int gomp_RecalculateHbonds(int);

extern int gomp_DeleteAllHbonds(void);
extern int gomp_DeleteAllHbondSearchData(void);
extern int gomp_DeleteAllHbondSubsets(void);

extern int gomp_SetHBondColour(float, float, float);
extern int gomp_GetHBondColour(float *, float *, float *);

extern int gomp_FindHbonds(int, int, int, const int *);
extern int gomp_SetHbondSubset(int, int, const int *);

extern int gomp_FindSSbondsAll(float);
extern int gomp_FindSSbondsI(int, float);

extern int gomp_ParseCalcConnList(
    const char *,const char *,const char *,const char *);
extern int gomp_ParseCalcHbondsList(
    const char *,const char *,const char *,int, const char *);
extern int gomp_ParseCalcHbondSubset(
    const char *,const char *,const char *);

extern int gomp_SetBondReconnectivityState(int);
extern int gomp_GetBondReconnectivityState(void);
extern int gomp_SetHBondReconnectivityState(int);
extern int gomp_GetHBondReconnectivityState(void);

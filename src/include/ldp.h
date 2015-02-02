#define LDP_ON  1
#define LDP_OFF 0

extern int gomp_CalculateLDParraydraw_ldp(int, int, float, float,
                                          float, float,
                                          float, float,
                                          float, float,
                                          float );
extern int  gomp_ApplyLdpSelectionMaskLine(
    const char *, const char *, const char *);
extern int  gomp_GetDisplayLDPmatrix(void);
extern int  gomp_SetDisplayLDPmatrix(int);
extern int  gomp_BuildLDParray(void);
extern int  gomp_GetNumLdpAtoms1(void);
extern int  gomp_GetNumLdpAtoms2(void);
extern const int *gomp_GetLdpAtomList1(void);
extern const int *gomp_GetLdpAtomList2(void);
int         gomp_PushAtomToLDP(int, const char *, const char *, const char *);

/* OS types                           */
#define WINDOWS_NT  "Windows NT"
#define WINDOWS_95  "Windows 95"

extern void        gomp_Get_cpu_secs(float * , float *);
extern const char *gomp_GetOS(void);
extern int         gomp_Get_date(int);
extern int         gomp_Get_proc_info(void);

#if HAVE_GETRUSAGE
extern struct rusage gomp_RUsage;
#elif defined(IRIX)
extern prusage_t     gomp_ProcessInfo;
#endif
extern int gomp_RunStatistics(void);

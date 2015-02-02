extern const char *gomp_ShowRootDir(void);
extern const char *gomp_ShowDataDir(void);
extern const char *gomp_ShowHomeDir(void);
extern const char *gomp_LongHostName(void);
extern const char *gomp_ShowHelpDir(void);
extern const char *gomp_ShowBinDir(void);
extern const char *gomp_ShowTempDir(void);

extern int gomp_GetEnv(int, const char *[]);

extern int gomp_SetEnvironmentFile(const char *);
extern int gomp_SetAtomParametersFile(const char *);

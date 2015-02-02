#include "gomstdio.h"

extern const char *gomp_GetDescText(void);
extern const char *gomp_GetTagText(void);
extern const char *gomp_GetDescText(void);
extern const char *gomp_GetAvailText(void);

extern int         gomp_PutTagText(const char *);
extern int         gomp_PutDescText(const char *);
extern int         gomp_PutAvailText(const char *);

extern int         gomp_SetFileSavingState(int);
extern int         gomp_gOpenMolNeedsSaving(void);

extern int         gomp_SplitFile(char *, const char *);

extern int         gomp_ReadOldModel(const char *);
extern int         gomp_WriteOldModel(const char *);

extern int         gomp_ReadTclInfo2FromFile(FILE *);

extern int gomp_ModelFileVersion;

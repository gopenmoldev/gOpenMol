#include "gomstdio.h"

extern int         gomp_ReadGaussianBasisLineInput(const char *);
extern int         gomp_ReadBasisSetRecords(FILE *);
extern int         gomp_PushBasis2Stack(const char *);
extern int         gomp_PushGBasisStartingPoint(int);
extern int         gomp_DeleteGaussianBasisLineInput(void);
extern int         gomp_GetNumberOfGaussianBasisSets(void);
extern const char *gomp_GetGaussianBasisTag(int);
extern const char *gomp_GetGaussianBasisSetEntry(int);

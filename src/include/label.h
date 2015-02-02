#include "parser_types.h"

/* label types */
#define RESIDUE_LABEL_TYPE  0
#define ATOM_LABEL_TYPE     1
#define FULL_LABEL_TYPE     2

extern int gomp_SetPlotLabelState(int);
extern int gomp_GetPlotLabelState(void);
extern int gomp_SetPlotLabelType(int);
extern int gomp_GetPlotLabelType(void);
extern int gomp_AtomLabelDataIsChanging(void);
extern int gomp_ParseAtomLabelList(int, const gom_SelectionList *);

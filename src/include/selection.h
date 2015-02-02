#define RESIDUE_SELECTION   0  /* atoms are selection by the whole residue */
#define ATOM_SELECTION      1  /* atoms are selected  by the atoms itself */

/* for selecting/unselecting structures */
#define STRUCTURE_SELECTION_OFF   0
#define STRUCTURE_SELECTION_ON    1

extern int gomp_GetAtomSelectionMode(void);
extern int gomp_SetAtomSelectionMode(int);

extern int gomp_MakeSelectionList(int,
                                const char *,
                                const char *,
                                const char *, int *);

extern int gomp_SetSelectionModeStatus(int);
extern int gomp_GetSelectionModeStatus(void);

extern int gomp_GetSelectedStructure(int);
extern int gomp_GetSelectedStructureStatus(void);
extern int gomp_ActivateSelectedStructure(int , int);
extern int gomp_SelectedStructure(int);
extern int gomp_DeleteSelectedStructureList(void);
extern int gomp_PushSelectedStructure(int);

extern int gomp_UpdateMolecStructList(void);


#define ATOM_LIST       0
#define RESIDUE_LIST    1
#define SEGMENT_LIST    2
#define SELECTION_LIST  3

#define LIST_ALL_STRUCTURES       -2
#define LIST_SELECTED_STRUCTURES  -1

typedef struct {
    int   Which,listType,SortLists;
    int   NAtoms;
    int   from,to,slong_total;
    int  *sel_list;
    int  *atom_list;
    int  *res_list;
    char *seg_list;
} SegmentResidueAtomList_t;

extern char *gomp_MakeIndexList(int,const int *,int,int,char);
extern int   gomp_AppendSegmentList(char **,int,int,const int *);
extern int   gomp_InitSegmentResidueAtomList(
    SegmentResidueAtomList_t *,int,int,int);
extern int   gomp_PreSegmentResidueAtomListStructure(
    SegmentResidueAtomList_t *,int);
extern int   gomp_PostSegmentResidueAtomListStructure(
    SegmentResidueAtomList_t *,int,int,const int *);
extern char *gomp_GetSegmentResidueAtomList(
    SegmentResidueAtomList_t *);
extern int   gomp_ReturnAndFreeSegmentResidueAtomList(
    SegmentResidueAtomList_t *);
extern int   gomp_FreeSegmentResidueAtomList(
    SegmentResidueAtomList_t *);

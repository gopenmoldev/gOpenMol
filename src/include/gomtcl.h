#include "parser.h"

extern int gomp_CalculateCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_CenterCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_ContourCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_CopyCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_DefineCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_DiagonalizeCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_DisplayCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_EditCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_ExportCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_FillCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_FindCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_HardcopyCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_HelpCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_ManipulateCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_MonitorCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_MOpenCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_PlotCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_PlumberCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_PythonCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_ImportCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_ResetCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_RotateCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_MSaveCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_ScaleCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_SelectCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_OldShowCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_TraceCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_TrajectoryCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_TranslateCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_WindowCommand(ClientData, Tcl_Interp *, int, const char *[]);
extern int gomp_PickingCommand(ClientData, Tcl_Interp *, int, const char *[]);

extern const gom_ParserArgumentList gomp_AtomCommand[];
extern const gom_ParserArgumentList gomp_RefreshCommand[];
extern const gom_ParserArgumentList gomp_ShowCommand[];
extern const gom_ParserArgumentList gomp_ShowAtomCommand[];
extern const gom_ParserArgumentList gomp_ShowCellCommand[];
extern const gom_ParserArgumentList gomp_ShowContourCommand[];
extern const gom_ParserArgumentList gomp_ShowCutPlaneCommand[];
extern const gom_ParserArgumentList gomp_ShowDisplaylistsCommand[];
extern const gom_ParserArgumentList gomp_ShowGbasisCommand[];
extern const gom_ParserArgumentList gomp_ShowLightCommand[];
extern const gom_ParserArgumentList gomp_ShowMaterialCommand[];
extern const gom_ParserArgumentList gomp_ShowMonitorCommand[];
extern const gom_ParserArgumentList gomp_ShowPlumberCommand[];
extern const gom_ParserArgumentList gomp_ShowWindowCommand[];
extern const gom_ParserArgumentList gomp_SyntaxCommand[];

#define gomp_CreateGomParser(Interp,cmd,Table) \
    Tcl_CreateObjCommand(Interp,cmd,gomp_ParserCommand, \
        CONST_CAST(void*,Table),NULL)

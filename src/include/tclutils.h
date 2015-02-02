#include <stdarg.h>

#include "gomcext.h"

#include "gomdefs.h"
#ifdef ENABLE_EXTENSIONS
#  include "gomlib/gomtclutils.h"
#  include "gomext/tclutils.h"
#endif

/***** PUBLIC GOMAPI BEGIN *****/

#include <tcl.h>

/** @weakgroup gom_doc
***/

/** @defgroup gom_doc_tcl Tcl Parser
*** @ingroup  gom_doc
*** @{
***/

/** @brief Send a Tcl script command to a Tcl parser.
***/
extern int         gomp_SendCommand2Parser(
    GOM_ARG( const char *, script ) );
/** @brief Send a Tcl object to a Tcl parser.
***/
extern int         gomp_SendTclObj2Parser(
    GOM_ARG( Tcl_Obj    *, obj ) );

/** @}
***/

/****** PUBLIC GOMAPI END ******/

extern int         gomp_SendFile2TclParser(const char *);
extern int         gomp_Tcl_AppInit(int , const char *[]);
extern int         gomp_TclRunScript(void);
extern Tcl_Interp *gomp_GetTclInterp(void);
extern int         gomp_SendTclReturn(const char *);
extern int         gomp_SendTclObjReturn(Tcl_Obj *);
extern int         gomp_StringMatch(const char *, const char *);

extern Tcl_Obj*    gomp_CreateTclList(const char *, ...)
                        CHECK_FORMAT_PRINTF;
extern Tcl_Obj*    gomp_CreateTclListv(const char *, va_list);

extern const char *gomp_GetNextFromParserStack(int, const char **);

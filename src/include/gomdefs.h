/***** PUBLIC GOMAPI BEGIN *****/
/*
                           Copyright (c) 2002 - 2003 by:
        Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
                      Confidential unpublished property of 
                              Leif Laaksonen  
                            All rights reserved

        Coded by: Eero Häkkinen
*/

/***** PUBLIC GOMAPI END *****/
#ifdef ENABLE_EXTENSIONS
#   include "gomlib/gomdefs.h"
#   include "gomext/gomdefs.h"
#else
#   ifndef INC_GOPENMOL_GOMDEFS
#   define INC_GOPENMOL_GOMDEFS
/***** PUBLIC GOMAPI BEGIN *****/

#if defined(WIN32)
#   ifdef __cplusplus
#       define DYNEXPORT_C extern "C" __declspec(dllexport)
#       define DYNIMPORT_C extern "C" __declspec(dllimport)
#   else
#       define DYNEXPORT_C extern     __declspec(dllexport)
#       define DYNIMPORT_C extern     __declspec(dllimport)
#   endif
#   define   DYNEXPORT   extern     __declspec(dllexport)
#   define   DYNIMPORT   extern     __declspec(dllimport)
#else
#   ifdef __cplusplus
#       define DYNEXPORT_C extern "C"
#       define DYNIMPORT_C extern "C"
#   else
#       define DYNEXPORT_C extern
#       define DYNIMPORT_C extern
#   endif
#   define   DYNEXPORT   extern
#   define   DYNIMPORT   extern
#endif

#ifdef DOXYGEN
    /** @mainpage
    *** See @ref gom_doc.
    ***/
#   define GOPENMOLAPI
#   define GOM_DECLARE_HANDLE(type) typedef type##_struct type
#   define GOM_ARG(type,name)       type name
#   ifdef GOPENMOL
        /** @defgroup gom_doc gOpenMol API
        ***/
#   else
        /** @defgroup gom_doc gOpenMol Extension API
        ***/
#   endif
#else
#   ifdef GOPENMOL
#       ifdef ENABLE_EXTENSIONS
#           define GOPENMOLAPI     DYNEXPORT_C
#       else
#           undef  GOPENMOLAPI
#       endif
#   else
#       define GOPENMOLAPI       DYNIMPORT_C
#   endif
#   define GOM_DECLARE_HANDLE(type) \
        typedef struct type##Handle { int unused; } type
#   define GOM_ARG(type,name) type
#endif

/***** PUBLIC GOMAPI END *****/
#   endif /* INC_GOPENMOL_GOMDEFS */
#endif /* ! ENABLE_EXTENSIONS */

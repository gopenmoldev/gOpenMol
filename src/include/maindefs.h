#if defined(HAVE_CONFIG_H)
#  if defined(_MSC_VER)
#    include <winconfig.h>

typedef __int64           intmax_t;
typedef unsigned __int64 uintmax_t;
#  else
#    include <config.h>

#    if ! SIZEOF_INTMAX_T
#      undef  SIZEOF_INTMAX_T
#      if SIZEOF_LONG_LONG_INT
#        define SIZEOF_INTMAX_T SIZEOF_LONG_LONG_INT
typedef long long int           intmax_t;
typedef unsigned long long int uintmax_t;
#      else
#        define SIZEOF_INTMAX_T SIZEOF_LONG_INT
typedef long int           intmax_t;
typedef unsigned long int uintmax_t;
#      endif
#    endif /* ! SIZEOF_INTMAX_T */

#    if ! SIZEOF_INTPTR_T
#      undef  SIZEOF_INTPTR_T
#      if SIZEOF_INT >= SIZEOF_VOID_P
#        define SIZEOF_INTPTR_T SIZEOF_INT
typedef int intptr_t;
#      elif SIZEOF_LONG_INT >= SIZEOF_VOID_P
#        define SIZEOF_INTPTR_T SIZEOF_VOID_P
typedef long int intptr_t;
#      elif SIZEOF_LONG_LONG_INT >= SIZEOF_VOID_P
#        define SIZEOF_LONG_LONG_INT SIZEOF_INTPTR_T
typedef long long int intptr_t;
#      elif SIZEOF_INTMAX_T >= SIZEOF_VOID_P
#        define SIZEOF_INTMAX_T SIZEOF_INTPTR_T
typedef intmax_t intptr_t;
#      else
#        define SIZEOF_INTPTR_T 0
#      endif
#    endif /* ! SIZEOF_INTPTR_T */

#    if ! SIZEOF_PTRDIFF_T
#      undef  SIZEOF_PTRDIFF_T
#      define SIZEOF_PTRDIFF_T SIZEOF_INT
typedef int ptrdiff_t;
#    endif /* ! SIZEOF_PTRDIFF_T */

#    if ! SIZEOF_SIZE_T
#      undef  SIZEOF_SIZE_T
#      define SIZEOF_SIZE_T SIZEOF_INT
typedef unsigned size_t;
#    endif /* ! SIZEOF_SIZE_T */
#  endif
#endif

#define GLUT
#define GOPENMOL

#define BUFF_LEN    512         /* default buffer length for text */

#define ON    1
#define OFF   0

/* predefine some files which will be needed at various stages */
#define LOGFILE            "gOPENMOL.LOG"
#define STARTUP            "gopenmolrc.tcl"
#define STREAMER_DATA_FILE "streamer.data"
#define ELEM_PARAM         "element.data"
#define ATOM_EQ            "atom_conversion.data"
#define ATOM_PARAM         "atom_param.data"
#define COLOUR_TABLE       "colour_table.data"

/* end of file names */

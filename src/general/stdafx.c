#include "maindefs.h"

#ifdef USE_PRECOMPILED_HEADERS

#ifdef _MSC_VER

/* C headers */
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <float.h>
#include <limits.h>
#include "gommath.h"
#include <setjmp.h>
#include <stdarg.h>
#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Windows headers */
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
/* Could somebody please tell me why windows.h defines this?       */
/* Some source files uses it as a squared radius so undef it.      */
#undef rad2
#undef small
#include <process.h>
#include <shellapi.h>
#include <winsock.h>

/* OpenGL headers */
#ifdef ENABLE_GRAPHICS
#if defined(GLUT)
#include <GL/glut.h>
#else
#include <GL/glaux.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>
#endif /* ENABLE_GRAPHICS */

/* Tcl/Tk headers */
#include <tcl.h>
#ifdef ENABLE_GRAPHICS
#include <tk.h>
#endif /* ENABLE_GRAPHICS */

/* gOpenMol headers */
#include "atom_param.h"
#include "axis.h"
#include "bond.h"
#include "cell.h"
#include "cluster.h"
#include "colors.h"
#include "colouring.h"
#include "contour.h"
#include "coord_file.h"
#include "coord_man.h"
#include "correl.h"
#include "drawscene.h"
#include "g_status.h"
#include "gaussian.h"
#include "gomclipbrd.h"
#include "gomdefs.h"
#include "gomendian.h"
#include "gomenv.h"
#include "gomfile.h"
#include "gomhelp.h"
#include "gomimage.h"
#include "gomstdio.h"
#include "gomlistener.h"
#include "gomlog.h"
#include "gommain.h"
#include "gommonitor.h"
#include "gomproc.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "gomtext.h"
#include "gomversion.h"
#include "gomwindow.h"
#include "hardcopy.h"
#include "label.h"
#include "ldp.h"
#include "light_model.h"
#include "listutils.h"
#include "maindefs.h"
#include "math_oper.h"
#include "measure.h"
#include "memalloc.h"
#include "model_file.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "msd.h"
#include "objseg.h"
#include "parser.h"
#include "parser_types.h"
#include "pbmplus_x.h"
#include "picking.h"
#include "plot_cpk.h"
#include "plot_molec.h"
#include "plot_plumber.h"
#include "plot_prop.h"
#include "plot.h"
#include "plumber.h"
#include "printmsg.h"
#include "projview.h"
#include "rdf.h"
#include "rforce.h"
#include "selection.h"
#include "sgi.h"
#include "stereo.h"
#include "tclutils.h"
#include "text_stack.h"
#include "tga.h"
#include "trace.h"
#include "trajectory.h"
#include "winconfig.h"
#include "x11wd.h"

#endif /* WIN32 */

#endif /* USE_PRECOMPILED_HEADERS */

#include "stdafx.h"

/* For ANSI compatibility.                  */
/* ANSI doesn't support empty source files. */
#include "memalloc.h"

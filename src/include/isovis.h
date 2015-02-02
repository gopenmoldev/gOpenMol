/*****************************************************************************
* 
* The following source code is in the public domain.
* Specifically, we give to the public domain all rights for future licensing
* of the source code, all resale rights, and all publishing rights.
* 
* We ask, but do not require, that the following message be included in all
* derived works:
* 
* Portions developed at the National Center for Supercomputing Applications at
* the University of Illinois at Urbana-Champaign.
* 
* THE UNIVERSITY OF ILLINOIS GIVES NO WARRANTY, EXPRESSED OR IMPLIED, FOR THE
* SOFTWARE AND/OR DOCUMENTATION PROVIDED, INCLUDING, WITHOUT LIMITATION,
* WARRANTY OF MERCHANTABILITY AND WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE
* 
****************************************************************************/

#if 0
static char OUTNAME[80] = "";
#endif
static char MY_NAME[80];
static int VERBOSE = 0;
#if 0
static char VSET_NAME[80] = "";
static char WFT_NAME[80] = "";
static int STORE_POLYGONS = 0;
static int DTM_OUTPUT = 0;
static char BYU_NAME[80] = "";
static int XDIM, YDIM, ZDIM, RAW_INPUT;
static int BBOX = 0, AUXLIGHT = 0;
static int SMOOTH = 0, NTSC = 0;
#endif
static FILE  *ISO_file;

static float ISO_xmin,ISO_xmax;
static float ISO_ymin,ISO_ymax;
static float ISO_zmin,ISO_zmax;
static float ISO_StepX;
static float ISO_StepY;
static float ISO_StepZ;
static int TypeOfSurface;        /* = 0 , unknown
                             = 1 , VSS surface
                             = 2 , orbital- or density surface
                             = 3 , probe surface
                           */
static int gomContourSwapBytes;

/* flag to allocate memory from 'malloc' first time through */
#if 0
static int ISO_first_alloc = 1;
#endif

#if defined(SGI_IMAGE)
static char SOUTNAME[80] = "";
#endif

#if 0
static int DISPLAY = 1;
#endif
static int SurfaceDisplayType = CONTOUR_TYPE_SOLID;  /* solid type by default */

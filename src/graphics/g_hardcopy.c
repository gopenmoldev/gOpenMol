/*

Copyright (c) 1995 - 2004 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#ifdef ENABLE_GRAPHICS
#include "gomstdio.h"
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <tcl.h>

#if defined(WIN32)
#include <windows.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>

#if defined(MOTIF)
#include <GL/GLwMDrawA.h>
#endif

#include "gomimage.h"
#include "hardcopy.h"
#include "printmsg.h"
#include "tclutils.h"
#endif /* ENABLE_GRAPHICS */

#include "stdafx.h"

#ifdef ENABLE_GRAPHICS

#if defined(HAVE_LIBJPEG)
#  undef EXTERN
#  undef HAVE_STDLIB_H
#  undef HAVE_STDDEF_H
#  ifdef WIN32
#    pragma warning (disable : 4005)
#    pragma warning (disable : 4142)
#  endif
#  include "jpeglib.h"
#endif

#define FILE_SELECTION_OK      10
#define FILE_SELECTION_CANCEL  11
#define FILE_SELECTION_HELP    12
#define FORMAT_PS             100
#define FORMAT_XWD            101
#define FORMAT_BMP            102
#define FORMAT_JPG            103
#define FORMAT_SGI            104
#define FORMAT_TGA            105
#define PS_EXTENSION          "*.ps"
#define XWD_EXTENSION         "*.xwd"
#define BMP_EXTENSION         "*.bmp"
#define SGI_EXTENSION         "*.rgb"
#define TGA_EXTENSION         "*.tga"
#define JPG_EXTENSION         "*.jpg"

#define   HARDCOPY "hardcopy_widget.hlp"

static struct {
    int Ready;
    int x1,y1;
    int x2,y2;
} ScrCoords;

static int  HardcopyIsActive = 0;

/**************************************************************************/
int gomp_HardcopyPostScript(const char *hard_file , const char *InputLine)
    /* make a PostScript hardcopy */
/**************************************************************************/
{

    int     WSize[4];
    static  unsigned char *pixels;
    char    OutText[BUFF_LEN];

    if(hard_file[0] == '\0') 
        return(1);

    glGetIntegerv(GL_VIEWPORT , WSize);

    pixels = calloc((WSize[2])  * 
                                     (WSize[3]) * 3 , 
                                     sizeof(unsigned char));
    if(pixels == (unsigned const char *)NULL)
    {
        gomp_PrintERROR("Memory allocation failed, can't save picture");
        return(1);
    }

    sprintf(OutText,"Preparing and writing PostScrip file '%s' ...",
            hard_file);
    gomp_PrintMessage(OutText);

    glPixelStorei(GL_PACK_ALIGNMENT , 1);

    glReadPixels((GLint)WSize[0],(GLint)WSize[1],
                 (GLsizei)(WSize[2]) ,
                 (GLsizei)(WSize[3]),
                 GL_RGB,GL_UNSIGNED_BYTE,(GLvoid*)pixels);

    (void) gomp_2ColorPs(hard_file , WSize[2] , WSize[3] , pixels ,
                       InputLine);

    free(pixels);

    gomp_PrintMessage("Done!");

    return(0);
}

/**************************************************************************/
int gomp_PrepareHardcopy()
/**************************************************************************/
{
    HardcopyIsActive = 0;
/*
  (void)gomp_DisplayMakeHardcopyWidget(1);
*/
    return(0);
}
/****************************************************************************/
int gomp_GetHardcopyActiveState()
/****************************************************************************/
{
    return(HardcopyIsActive);
}
/****************************************************************************/
int gomp_PutHardcopyCoordinates(int x, int y)
/****************************************************************************/
{
    if(!ScrCoords.Ready) {
        ScrCoords.x1 = x;
        ScrCoords.y1 = y;
        ScrCoords.Ready = 1;
    }
    else if(ScrCoords.Ready) {
        ScrCoords.x2 = x;
        ScrCoords.y2 = y;
        ScrCoords.Ready = 2;
    }
    return(0);
}

/**************************************************************************/
int gomp_HardcopyXWD(const char *hard_file)
    /* make a XWD hardcopy */
/**************************************************************************/
{
    int     WSize[4];
    static  unsigned char *pixels;
    char    InputLine[BUFF_LEN];
    char    OutText[BUFF_LEN];

    if(hard_file[0] == '\0') 
        return(1);

    glGetIntegerv(GL_VIEWPORT , WSize);

    pixels = calloc((WSize[2])  * 
                                     (WSize[3]) * 3 , 
                                     sizeof(unsigned char));
    if(pixels == (unsigned const char *)NULL)
    {
        gomp_PrintERROR("Memory allocation failed, can't save picture");
        return(1);
    }

    sprintf(OutText,"Preparing and writing XWD file '%s' ...",
            hard_file);
    gomp_PrintMessage(OutText);

    glPixelStorei(GL_PACK_ALIGNMENT , 1);

    glReadPixels((GLint)WSize[0],(GLint)WSize[1],
                 (GLsizei)(WSize[2]) ,
                 (GLsizei)(WSize[3]),
                 GL_RGB,GL_UNSIGNED_BYTE,(GLvoid*)pixels);

    (void) gomp_2XWD(hard_file , WSize[2] , WSize[3] , pixels ,
                   InputLine);

    free(pixels);

    gomp_PrintMessage("Done!");

    return(0);
}

/**************************************************************************/
int gomp_HardcopyBMP(const char *hard_file)
    /* make a BMP hardcopy */
/**************************************************************************/
{     
    int     WSize[4];
    static  unsigned char *pixels;
    char    OutText[BUFF_LEN];
    char    InputLine[BUFF_LEN];

    if(hard_file[0] == '\0') 
        return(1);

    glGetIntegerv(GL_VIEWPORT , WSize);

    pixels = calloc((WSize[2]) * 
                                     (WSize[3]) * 3 , 
                                     sizeof(unsigned char));
    if(pixels == (unsigned const char *)NULL)
    {
        gomp_PrintERROR("Memory allocation failed, can't save picture");
        return(1);
    }

    sprintf(OutText,"Preparing and writing BMP file '%s' ...",
            hard_file);
    gomp_PrintMessage(OutText);

    glPixelStorei(GL_PACK_ALIGNMENT , 1);

    glReadPixels((GLint)WSize[0],(GLint)WSize[1],
                 (GLsizei)(WSize[2]) ,
                 (GLsizei)(WSize[3]),
                 GL_RGB,GL_UNSIGNED_BYTE,(GLvoid*)pixels);

    (void) gomp_2BMP(hard_file , WSize[2] , WSize[3] , pixels ,
                   InputLine);

    free(pixels);

    gomp_PrintMessage("Done!");
    return(0);
}
/**************************************************************************/
int gomp_HardcopyRGB(const char *hard_file)
    /* make a XWD hardcopy */
/**************************************************************************/
{
    int   WSize[4];
    static  unsigned char *pixels;
    char    InputLine[BUFF_LEN];
    char    OutText[BUFF_LEN];

    if(hard_file[0] == '\0') 
        return(1);

    glGetIntegerv(GL_VIEWPORT , WSize);

/*
  w = gomp_GetStructureWidget();
  CDisplay = XtDisplay(w);
  CWindow  = XtWindow(w);

  XGetWindowAttributes( CDisplay , CWindow , &Wattribs);

  printf("%d %d %d %d\n",Wattribs.x,Wattribs.y,Wattribs.width,
  Wattribs.height);
*/

    pixels = calloc((WSize[2])  * 
                                     (WSize[3]) * 3 , 
                                     sizeof(unsigned char));

    if(pixels == (unsigned const char *)NULL)
    {
        gomp_PrintERROR("Memory allocation failed, can't save picture");
        return(1);
    }

    sprintf(OutText,"Preparing and writing RGB file '%s' ...",
            hard_file);
    gomp_PrintMessage(OutText);

    glPixelStorei(GL_PACK_ALIGNMENT , 1);

    glReadPixels((GLint)WSize[0],(GLint)WSize[1],
                 (GLsizei)(WSize[2]) ,
                 (GLsizei)(WSize[3]),
                 GL_RGB,GL_UNSIGNED_BYTE,(GLvoid*)pixels);

    (void) gomp_2SGI(hard_file , WSize[2] , WSize[3] , pixels ,
                   InputLine);

    free(pixels);

    gomp_PrintMessage("Done!");
    return(0);
}

/**************************************************************************/
int gomp_HardcopyTGA(const char *hard_file)
    /* make a TGA hardcopy */
/**************************************************************************/
{
    int     WSize[4];
    static  unsigned char *pixels;
    char    InputLine[BUFF_LEN];
    char    OutText[BUFF_LEN];

    if(hard_file[0] == '\0') 
        return(1);

    glGetIntegerv(GL_VIEWPORT , WSize);

    pixels = calloc((WSize[2])  * 
                                     (WSize[3]) * 3 , 
                                     sizeof(unsigned char));
    if(pixels == (unsigned const char *)NULL)
    {
        gomp_PrintERROR("Memory allocation failed, can't save picture");
        return(1);
    }

    sprintf(OutText,"Preparing and writing Targa file '%s' ...",
            hard_file);
    gomp_PrintMessage(OutText);

    glPixelStorei(GL_PACK_ALIGNMENT , 1);

    glReadPixels((GLint)WSize[0],(GLint)WSize[1],
                 (GLsizei)(WSize[2]) ,
                 (GLsizei)(WSize[3]),
                 GL_RGB,GL_UNSIGNED_BYTE,(GLvoid*)pixels);

    (void) gomp_2TGA(hard_file , WSize[2] , WSize[3] , pixels ,
                   InputLine);

    free(pixels);

    gomp_PrintMessage("Done!");

    return(0);
}

/**************************************************************************/
int gomp_HardcopyWindowsMetafile(const char *hard_file)
    /* make a Windows Metafile hardcopy */
/**************************************************************************/
{     
#if defined(WIN32)

    HDC     hdcMeta;

    if(hard_file[0] == '\0') 
        return(1);

    hdcMeta  = CreateMetaFile(hard_file);

    (void)CloseMetaFile(hdcMeta);

    gomp_PrintMessage("Done!");

#else

    gomp_PrintERROR("Windows Metafile is not implemented on this architecture");
#endif

    return(0);
}
/**************************************************************************/
int gomp_HardcopyJPEG(const char *hard_file)
    /* make a BMP hardcopy */
/**************************************************************************/
{     
#ifdef HAVE_LIBJPEG
    int     WSize[4];
    static  unsigned char *pixels;
    char    OutText[BUFF_LEN];
    FILE   *outfile;
    struct  jpeg_compress_struct cinfo;
    struct  jpeg_error_mgr jerr;
    int     JPEGquality;
    const char *value;


    if(hard_file[0] == '\0') 
        return(1);

    glGetIntegerv(GL_VIEWPORT , WSize);

    pixels = calloc((WSize[2]) * 
                                     (WSize[3]) * 3 , 
                                     sizeof(unsigned char));
    if(pixels == (unsigned const char *)NULL)
    {
        gomp_PrintERROR("Memory allocation failed, can't save picture");
        return(1);
    }

    sprintf(OutText,"Preparing and writing JPG file '%s' ...",
            hard_file);
    gomp_PrintMessage(OutText);

    glPixelStorei(GL_PACK_ALIGNMENT , 1);

    glReadPixels((GLint)WSize[0],(GLint)WSize[1],
                 (GLsizei)(WSize[2]) ,
                 (GLsizei)(WSize[3]),
                 GL_RGB,GL_UNSIGNED_BYTE,(GLvoid*)pixels);

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    if ((outfile = fopen(hard_file, "wb")) == NULL) {
        gomp_PrintERROR("can't open output file");
        return(1);
    }

    jpeg_stdio_dest(&cinfo, outfile);

    cinfo.image_width  = WSize[2];  /* image width and height, in pixels */
    cinfo.image_height = WSize[3];
    cinfo.input_components = 3;     /* # of color components per pixel */
    cinfo.in_color_space = JCS_RGB;    /* colorspace of input image */

    jpeg_set_defaults(&cinfo);

    value  = Tcl_GetVar(gomp_GetTclInterp() , "gomJPEGquality", TCL_GLOBAL_ONLY);

    if(value) {
        JPEGquality = atoi(value);
        if(JPEGquality < 1 || JPEGquality > 100) JPEGquality = 75;
        sprintf(OutText,"Saving JPG file using quality %d",JPEGquality);
        gomp_PrintMessage(OutText);
        jpeg_set_quality(&cinfo, JPEGquality, TRUE /* limit to baseline-JPEG values */);
    }

    jpeg_start_compress(&cinfo, TRUE);

    {
        JSAMPROW row_pointer[1];    /* pointer to a single row */
        int row_stride;         /* physical row width in buffer */
        int next_scanline;

        row_stride = WSize[2] * 3;  /* JSAMPLEs per row in image_buffer */

        next_scanline = cinfo.image_height;
        while (cinfo.next_scanline < cinfo.image_height) {
            next_scanline--; /* flip image */
            row_pointer[0] = & pixels[next_scanline * row_stride];
            jpeg_write_scanlines(&cinfo, row_pointer, 1);
        }
    }
    jpeg_finish_compress(&cinfo);

    jpeg_destroy_compress(&cinfo);
 
    fclose(outfile);

    free(pixels);

    gomp_PrintMessage("Done!");
    return(0);
#else /* HAVE_LIBJPEG */
    gomp_PrintERROR("JPEG support is not compiled into gOpenMol");
    return(1);
#endif
}
#else /* ENABLE_GRAPHICS */
extern int i;
#endif

/*
  Enhancements 2002, 2004 by:
  Eero Häkkinen
*/
/* gotoxwd.c - produce a color X11 window dump
**
** To be used with gOpenMol. 
** Leif Laaksonen  Center for Scientific Computing 1995
**
** Copyright (C) 1989, 1991 by Jef Poskanzer.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <string.h>

#include "gomimage.h"
#include "pbmplus_x.h"
#include "printmsg.h"
#include "x11wd.h"

#include "stdafx.h"

/***************************************************************************/
int gomp_2XWD(const char *FileName , int xsize , int ysize ,
            unsigned const char *pixels   , const char *InputLine)
/***************************************************************************/
{
    const char * dumpname;
    int argn, rows, cols, row, col;
    int pseudodepth;
    int forcedirect, direct;
    X11WDFileHeader h11;
    FILE *File_p;
    int  Loop;

    argn        = 1;
    pseudodepth = 8;
    forcedirect = 1;
    direct      = 1;
    cols        = xsize;
    rows        = ysize;
    dumpname    = "stdin";

    File_p = fopen(FileName , "wb");
    if(File_p == NULL) {
        gomp_PrintERROR("Can't open output file");
        return(1);
    }

    /* Set up the header. */
    h11.header_size      = sizeof(h11) + strlen( dumpname ) + 1;
    h11.file_version     = X11WD_FILE_VERSION;
    h11.pixmap_format    = ZPixmap;
    h11.pixmap_width     = cols;
    h11.pixmap_height    = rows;
    h11.xoffset          = 0;
    h11.byte_order       = MSBFirst;
    h11.bitmap_bit_order = MSBFirst;
    h11.window_width     = cols;
    h11.window_height    = rows;
    h11.window_x         = 0;
    h11.window_y         = 0;
    h11.window_bdrwidth  = 0;

    h11.pixmap_depth     = 24;
    h11.bitmap_unit      = 32;
    h11.bitmap_pad       = 32;
    h11.bits_per_pixel   = 32;
    h11.visual_class     = DirectColor;
    h11.colormap_entries = 256;
    h11.ncolors          = 0;
    h11.red_mask         = 0xff0000;
    h11.green_mask       = 0xff00;
    h11.blue_mask        = 0xff;
    h11.bytes_per_line   = cols * 4;
    h11.bits_per_rgb     = h11.pixmap_depth;

    /* Write out the header in big-endian order. */
    pm_writebiglong( File_p, h11.header_size );
    pm_writebiglong( File_p, h11.file_version );
    pm_writebiglong( File_p, h11.pixmap_format );
    pm_writebiglong( File_p, h11.pixmap_depth );
    pm_writebiglong( File_p, h11.pixmap_width );
    pm_writebiglong( File_p, h11.pixmap_height );
    pm_writebiglong( File_p, h11.xoffset );
    pm_writebiglong( File_p, h11.byte_order );
    pm_writebiglong( File_p, h11.bitmap_unit );
    pm_writebiglong( File_p, h11.bitmap_bit_order );
    pm_writebiglong( File_p, h11.bitmap_pad );
    pm_writebiglong( File_p, h11.bits_per_pixel );
    pm_writebiglong( File_p, h11.bytes_per_line );
    pm_writebiglong( File_p, h11.visual_class );
    pm_writebiglong( File_p, h11.red_mask );
    pm_writebiglong( File_p, h11.green_mask );
    pm_writebiglong( File_p, h11.blue_mask );
    pm_writebiglong( File_p, h11.bits_per_rgb );
    pm_writebiglong( File_p, h11.colormap_entries );
    pm_writebiglong( File_p, h11.ncolors );
    pm_writebiglong( File_p, h11.window_width );
    pm_writebiglong( File_p, h11.window_height );
    pm_writebiglong( File_p, h11.window_x );
    pm_writebiglong( File_p, h11.window_y );
    pm_writebiglong( File_p, h11.window_bdrwidth );

    /* Write out the dump name. */
    fwrite( dumpname, 1, strlen( dumpname ) + 1, File_p );

    /* Finally, write out the data. */
    for ( row = 0; row < rows; row++ )
    {
        Loop = rows - row - 1;
        for ( col = 0; col < 3 * cols; col += 3)
        {
            unsigned long ul;
            register unsigned long valr,valg,valb;

            valr = (int)pixels[    col + 3 * Loop * cols];
            valg = (int)pixels[1 + col + 3 * Loop * cols];
            valb = (int)pixels[2 + col + 3 * Loop * cols];
            ul = ( ( valr ) << 16 ) |
                ( ( valg ) <<  8 ) |
                ( ( valb ) );
            fwrite( &ul, sizeof(ul), 1, File_p );
        }
    }

    fclose(File_p);

    return( 0 );
}

/***************************************************************************/
int pm_writebiglong(FILE *out, long l )
/***************************************************************************/
{
    (void) putc( ( l >> 24 ) & 0xff, out );
    (void) putc( ( l >> 16 ) & 0xff, out );
    (void) putc( ( l >> 8 ) & 0xff, out );
    (void) putc( l & 0xff, out );
    return 0;
}

/***************************************************************************/
int pm_writelittlelong(FILE *out, long l )
/***************************************************************************/
{
    (void) putc( l & 0xff, out );
    (void) putc( ( l >> 8 ) & 0xff, out );
    (void) putc( ( l >> 16 ) & 0xff, out );
    (void) putc( ( l >> 24 ) & 0xff, out );
    return 0;
}

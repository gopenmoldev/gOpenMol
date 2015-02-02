/* ppmtotga.c - read a portable pixmap and produce a TrueVision Targa file
**
** Copyright (C) 1989, 1991 by Mark Shand and Jef Poskanzer
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

#include "gomimage.h"
#include "printmsg.h"
#include "tga.h"

#include "stdafx.h"

/* Max number of colors allowed for colormapped output. */
#define MAXCOLORS 256

/* Forward routines. */
static void writetga( struct ImageHeader* tgaP, const char * id );
static int writechar( FILE* out, unsigned char s );

static FILE* ifp_tga;

/***************************************************************************/
int gomp_2TGA(const char *FileName , int xsize , int ysize ,
            unsigned const char *pixels   , const char *InputLine)
/***************************************************************************/
{
    int    rows, cols, row, col;
    struct ImageHeader tgaHeader;
    int  maxval;

    ifp_tga = fopen(FileName , "wb");
    if(ifp_tga == NULL) {
        gomp_PrintERROR("Can't open output file");
        return(1);
    }

    cols                = xsize;
    rows                = ysize;
    maxval              = 255;

    tgaHeader.ImgType   = TGA_RGB;
    tgaHeader.PixelSize = 24;
    tgaHeader.X_org_lo  = tgaHeader.X_org_hi = 0;
    tgaHeader.Y_org_lo  = tgaHeader.Y_org_hi = 0;
    tgaHeader.Width_lo  = cols % 256;
    tgaHeader.Width_hi  = cols / 256;
    tgaHeader.Height_lo = rows % 256;
    tgaHeader.Height_hi = rows / 256;
    tgaHeader.AttBits   = 0;
    tgaHeader.Rsrvd     = 0;
    tgaHeader.IntrLve   = 0;
    tgaHeader.OrgBit    = 0;
    tgaHeader.IDLength  = 0;
    tgaHeader.Index_lo  = 0;
    tgaHeader.Index_hi  = 0;
    tgaHeader.CoMapType = 0;
    tgaHeader.Length_lo = 0;
    tgaHeader.Length_hi = 0;
    tgaHeader.CoSize    = 0;


    /* Write out the Targa header. */
    writetga( &tgaHeader, (const char *) 0 );

    /* Write out the pixels */
    for ( row = 0; row < rows; ++row )
    {
        for ( col = 0; col < cols; ++col  ) {

/* in order blue , green , red */
            (void)writechar(ifp_tga , pixels[2 + 3 * col + row * (3 * cols)]);  
            (void)writechar(ifp_tga , pixels[1 + 3 * col + row * (3 * cols)]);  
            (void)writechar(ifp_tga , pixels[    3 * col + row * (3 * cols)]);  
        }
    }

    fclose(ifp_tga);

    return(0);
}

static void writetga( struct ImageHeader* tgaP ,const char * id)
{
    unsigned char flags;

/* general */
    writechar(ifp_tga , tgaP->IDLength );
    writechar(ifp_tga , tgaP->CoMapType );
    writechar(ifp_tga , tgaP->ImgType );
/* color map */
    writechar(ifp_tga , tgaP->Index_lo );
    writechar(ifp_tga , tgaP->Index_hi );
    writechar(ifp_tga , tgaP->Length_lo );
    writechar(ifp_tga , tgaP->Length_hi );
    writechar(ifp_tga , tgaP->CoSize );
    /* image data */
    writechar(ifp_tga , tgaP->X_org_lo );
    writechar(ifp_tga , tgaP->X_org_hi );
    writechar(ifp_tga , tgaP->Y_org_lo );
    writechar(ifp_tga , tgaP->Y_org_hi );
    writechar(ifp_tga , tgaP->Width_lo );
    writechar(ifp_tga , tgaP->Width_hi );
    writechar(ifp_tga , tgaP->Height_lo );
    writechar(ifp_tga , tgaP->Height_hi );
    writechar(ifp_tga , tgaP->PixelSize );
    flags = ( tgaP->AttBits & 0xf ) | ( ( tgaP->Rsrvd & 0x1 ) << 4 ) |
        ( ( tgaP->OrgBit & 0x1 ) << 5 ) | ( ( tgaP->OrgBit & 0x3 ) << 6 );
    writechar( ifp_tga , flags );

    if ( tgaP->IDLength )
        fwrite( id, 1, (int) tgaP->IDLength, ifp_tga );
}

static int writechar( FILE* out, unsigned char s )
{
    (void) putc( s , out );
    return 0;
}

/*  

Copyright (c) 1995 - 2004 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

*/

#include "maindefs.h"

#include "gomstdio.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "gomimage.h"
#include "printmsg.h"

#include "stdafx.h"

/*

Code originating from Matti Gröhn at CSC
modified by Leif Laaksonen CSC 1996

*/

#if 0
static int bmp_writebiglong(FILE *, long );
static int bmp_writebigshort( FILE *, short );
#endif
static int bmp_writelittlelong(FILE *, long );
static int bmp_writelittleshort( FILE *, short );

/****************************************************************************/
int gomp_2BMP(const char *outfile , int xsize , int ysize , unsigned const char *pixels, const char *Line)
/****************************************************************************/
{
    FILE *fp;

    long nbytes, npixels;
    long header_incl, header_excl;
    short ncolorbits;
    long zerobytes;
    short intone;
    int i,j,loop;
    char cbzero;
    char ch1, ch2;
    int  xsize4;
    int  xsize3;

    i = xsize * 3; 
    while(i - (i / 4) * 4) i++;
    xsize4 = i;
    xsize3 = xsize * 3;

    ch1          = 'B';                /* prepare header */
    ch2          = 'M';
    npixels      = xsize * ysize * 3;
    header_incl  = 54;
    header_excl  = 40;
    nbytes       = (xsize4 * ysize) + header_incl;
    intone       = 1;
    ncolorbits   = 24;
    zerobytes    = 0;
    cbzero        = (char)0;

    fp = fopen(outfile , "wb");
    if(fp == NULL) {
        gomp_PrintERROR("Can't open output file");
        return(1);
    }

/* windows bmp file header information */

    fwrite(&ch1,sizeof(char),1,fp);
    fwrite(&ch2,sizeof(char),1,fp);

    bmp_writelittlelong(fp,  nbytes );
    bmp_writelittlelong(fp,  zerobytes );
    bmp_writelittlelong(fp,  header_incl );

    bmp_writelittlelong(fp,  header_excl );
    bmp_writelittlelong(fp,  xsize );
    bmp_writelittlelong(fp,  ysize );
    bmp_writelittleshort(fp, intone );
    bmp_writelittleshort(fp, ncolorbits );

    bmp_writelittlelong(fp,  zerobytes );
    bmp_writelittlelong(fp,  npixels );
    bmp_writelittlelong(fp,  zerobytes );
    bmp_writelittlelong(fp,  zerobytes );
    bmp_writelittlelong(fp,  zerobytes );
    bmp_writelittlelong(fp,  zerobytes );
  
    /* saving the pixel information */

    loop = 0;
    for(i = 0 ; i < ysize ; i++) {
        for(j = 0; j < xsize ; j++)
        {
            fwrite(&pixels[loop + 2],sizeof(unsigned char),1,fp);
            fwrite(&pixels[loop + 1],sizeof(unsigned char),1,fp);
            fwrite(&pixels[loop    ],sizeof(unsigned char),1,fp);
            loop += 3;
        }
        for(j = 0 ; j < (xsize4 - xsize3) ; j++)
            fwrite(&cbzero,sizeof(unsigned char),1,fp);
    }

    fclose(fp);

    return(0);
}
#if 0
/***************************************************************************/
int bmp_writebiglong(FILE *out, long l )
/***************************************************************************/
{
    (void) putc( ( l >> 24 ) & 0xff, out );
    (void) putc( ( l >> 16 ) & 0xff, out );
    (void) putc( ( l >> 8 ) & 0xff, out );
    (void) putc( l & 0xff, out );
    return 0;
}
/***************************************************************************/
int bmp_writebigshort( FILE* out, short s )
/***************************************************************************/
{
    (void) putc( ( s >> 8 ) & 0xff, out );
    (void) putc( s & 0xff, out );

    return 0;
}
#endif
/***************************************************************************/
int bmp_writelittleshort( FILE* out, short s )
/***************************************************************************/
{
    (void) putc( s & 0xff, out );
    (void) putc( ( s >> 8 ) & 0xff, out );
    return 0;
}
/***************************************************************************/
int bmp_writelittlelong(FILE* out, long l )
/***************************************************************************/
{
    (void) putc( l & 0xff, out );
    (void) putc( ( l >> 8 ) & 0xff, out );
    (void) putc( ( l >> 16 ) & 0xff, out );
    (void) putc( ( l >> 24 ) & 0xff, out );
    return 0;
}

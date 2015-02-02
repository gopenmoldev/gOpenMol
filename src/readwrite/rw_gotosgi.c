/*
  Enhancements 2002, 2004 by:
  Eero Häkkinen
*/
/* gotosgi.c -  write an SGI image
**
** Copyright (C) 1994 by Ingo Wilken (Ingo.Wilken@informatik.uni-oldenburg.de)
**
** Based on the SGI image description v0.9 by Paul Haeberli (paul@sgi.comp)
** Available via ftp from sgi.com:graphics/SGIIMAGESPEC
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
**
** 29Jan94: first version
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>

#include "gomimage.h"
#include "pbmplus_x.h"
#include "printmsg.h"
#include "sgi.h"

#include "stdafx.h"

typedef short       ScanElem;
typedef struct {
    ScanElem *  data;
    long        length;
} ScanLine;

/* prototypes */
/*static void put_big_short ARGS((short s));*/
static int writebigshort( FILE* out, short s );
static void put_big_short(short);
static void put_big_long(long l);
#define put_byte(b)     (void)(putc((unsigned char)(b), ifp_sgi))
/*static void put_short_as_byte ARGS((short s));*/
static void put_short_as_byte(short);
static void write_header(int cols, int rows, int maxval, int bpc, int dimensions, int channels, const char *imagename);
static void write_table(long *table, int tabsize);
static void write_channels(int cols, int rows, int channels, void (*put)(short) );
static long * build_channels(unsigned const char *pixels, int cols, int rows, int maxval, int format, int bpc, int channels);
static ScanElem *compress(ScanElem *temp, int row, int rows, int cols, int chan_no, long *table, int bpc);
static int rle_compress(ScanElem *inbuf, int cols);
static void * xmalloc(int bytes);
#define MALLOC(n, type)     (type *)xmalloc((n) * sizeof(type))

#define WORSTCOMPR(x)   (2*(x) + 2)


#define MAXVAL_BYTE     255
#define MAXVAL_WORD     65535

static char storage = STORAGE_RLE;
static ScanLine * channel[3];
static ScanElem * rletemp;
static FILE *ifp_sgi;

/***************************************************************************/
int gomp_2SGI(const char *FileName , int xsize , int ysize ,
            unsigned const char *pixels   , const char *InputLine)
/***************************************************************************/
{
    int cols, rows, format;
    int maxval;
    const char *imagename;
    int bpc, dimensions, channels;
    long *table = NULL;

    ifp_sgi = fopen(FileName , "wb");
    if(ifp_sgi == NULL) {
        gomp_PrintERROR("Can't open output file");
        return(1);
    }

    imagename = FileName;

    cols       = xsize;
    rows       = ysize;
    maxval     = 255;
    format     = 0; /* dummy, not used */

    dimensions = 3;
    channels   = 3;

    bpc = 1;

    table = build_channels(pixels, cols, rows, maxval, format, bpc, channels);

/* write the stuff now */
    write_header(cols, rows, maxval, bpc, dimensions, channels, imagename);
    if( table )
        write_table(table, rows * channels);
    if( bpc == 1 )
        write_channels(cols, rows, channels, put_short_as_byte);
    else
        write_channels(cols, rows, channels, put_big_short);

    fclose(ifp_sgi);

    return(0);
}


static void
write_header(int cols, int rows, int maxval, int bpc, int dimensions, int channels, const char *imagename)
{
    int i;

#ifdef DEBUG
    pm_message("writing header");
#endif

    put_big_short(SGI_MAGIC);
    put_byte(storage);
    put_byte((char)bpc);
    put_big_short(dimensions);
    put_big_short(cols);
    put_big_short(rows);
    put_big_short(channels);
    put_big_long(0);                /* PIXMIN */
    put_big_long(maxval);           /* PIXMAX */
    for( i = 0; i < 4; i++ )
        put_byte(0);
    for( i = 0; i < 79 && imagename[i] != '\0'; i++ )
        put_byte(imagename[i]);
    for(; i < 80; i++ )
        put_byte(0);
    put_big_long(CMAP_NORMAL);
    for( i = 0; i < 404; i++ )
        put_byte(0);
}


static void
write_table(long *table, int tabsize)
{
    int i;
    long offset;

#ifdef DEBUG
    pm_message("writing table");
#endif

    offset = HeaderSize + tabsize * 8;
    for( i = 0; i < tabsize; i++ ) {
        put_big_long(offset);
        offset += table[i];
    }
    for( i = 0; i < tabsize; i++ )
        put_big_long(table[i]);
}


static void
write_channels(int cols, int rows, int channels, void (*put)(short))
{
    int i, row, col;

#ifdef DEBUG
    pm_message("writing image data");
#endif

    for( i = 0; i < channels; i++ ) {
        for( row = 0; row < rows; row++ ) {
            for( col = 0; col < channel[i][row].length; col++ ) {
                (*put)(channel[i][row].data[col]);
            }
        }
    }
}

static void *xmalloc(int bytes)
{
    void *mem;

    if( bytes == 0 )
        return NULL;

    mem = malloc(bytes);
    if( mem == NULL )
        printf("out of memory allocating %d bytes", bytes);
    return mem;
}

static void
put_big_short(short s)
{
    if ( writebigshort( ifp_sgi, s ) == -1 )
        gomp_PrintERROR("SGI file: write error" );
}


static void
put_big_long(long l)
{
    if ( pm_writebiglong( ifp_sgi, l ) == -1 )
        gomp_PrintERROR("SGI file: write error" );
}


static void
put_short_as_byte(short s)
{
    put_byte((unsigned char)s);
}


static long *
build_channels(unsigned const char *pixels, int cols, int rows, int maxval, int format, int bpc, int channels)
/*
  unsigned const char *pixels;
  int cols, rows;
  int maxval;
  int format, bpc, channels;
*/
{
    int i, row, col, sgirow;
    long *table = NULL;
    ScanElem *temp;
    int Loop;
    int colx;

#ifdef DEBUG
    pm_message("building channels");
#endif

    if( storage != STORAGE_VERBATIM ) {
        table = MALLOC(channels * rows, long);
        rletemp = MALLOC(WORSTCOMPR(cols), ScanElem);
    }
    temp = MALLOC(cols, ScanElem);

    for( i = 0; i < channels; i++ )
        channel[i] = MALLOC(rows, ScanLine);

    for( row = 0, sgirow = rows-1; row < rows; row++, sgirow-- ) {

        Loop = rows - row - 1;
        for( col  = 0; col < cols; col ++ ) {
            colx = 3 * col;
            temp[col] = (ScanElem)pixels[    colx + 3 * Loop * cols];
        }
        temp = compress(temp, sgirow, rows, cols, 0, table, bpc);
        for( col  = 0; col < cols; col ++ ) {
            colx = 3 * col;
            temp[col] = (ScanElem)pixels[1 + colx + 3 * Loop * cols];
        }
        temp = compress(temp, sgirow, rows, cols, 1, table, bpc);
        for( col  = 0; col < cols; col ++ ) {
            colx = 3 * col;
            temp[col] = (ScanElem)pixels[2 + colx + 3 * Loop * cols];
        }
        temp = compress(temp, sgirow, rows, cols, 2, table, bpc);
    }

    free(temp);
    if( table )
        free(rletemp);
    return table;
}


static ScanElem *
compress(ScanElem *temp, int row, int rows, int cols, int chan_no, long *table, int bpc)
/*
  ScanElem *temp;
  int row, rows, cols, chan_no;
  long *table;
  int bpc;
*/
{
    int len, i, tabrow;
    ScanElem *p;

    switch( storage ) {
    case STORAGE_VERBATIM:
        channel[chan_no][row].length = cols;
        channel[chan_no][row].data = temp;
        temp = MALLOC(cols, ScanElem);
        break;
    case STORAGE_RLE:
        tabrow = chan_no * rows + row;
        len = rle_compress(temp, cols);     /* writes result into rletemp */
        channel[chan_no][row].length = len;
        channel[chan_no][row].data = p = MALLOC(len, ScanElem);
        for( i = 0; i < len; i++, p++ )
            *p = rletemp[i];
        table[tabrow] = len * bpc;
        break;
    default:
        gomp_PrintERROR("unknown storage type - can't happen");
    }
    return temp;
}


/*
  slightly modified RLE algorithm from ppmtoilbm.c
  written by Robert A. Knop (rknop@mop.caltech.edu)
*/
static int
rle_compress(ScanElem *inbuf, int size)
/*
  ScanElem *inbuf;
  int size;
*/
{
    int in, out, hold, count;
    ScanElem *outbuf = rletemp;

    in=out=0;
    while( in<size ) {
        if( (in<size-1) && (inbuf[in]==inbuf[in+1]) ) {     /*Begin replicate run*/
            for( count=0,hold=in; in<size && inbuf[in]==inbuf[hold] && count<127; in++,count++)
                ;
            outbuf[out++]=(ScanElem)(count);
            outbuf[out++]=inbuf[hold];
        }
        else {  /*Do a literal run*/
            hold=out; out++; count=0;
            while( ((in>=size-2)&&(in<size)) || ((in<size-2) && ((inbuf[in]!=inbuf[in+1])||(inbuf[in]!=inbuf[in+2]))) ) {
                outbuf[out++]=inbuf[in++];
                if( ++count>=127 )
                    break;
            }
            outbuf[hold]=(ScanElem)(count | 0x80);
        }
    }
    outbuf[out++] = (ScanElem)0;     /* terminator */
    return(out);
}

int writebigshort( FILE* out, short s )
{
    (void) putc( ( s >> 8 ) & 0xff, out );
    (void) putc( s & 0xff, out );
    return 0;
}

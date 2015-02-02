/*
  Enhancements 2002, 2004 by:
  Eero Häkkinen
*/
 
/*      tocolps -
 *              Convert a color  image to color PostScript.
 *      This knows how to print images with either 1, 2, 4, or 8 bits per
 *      pixel, and how to generate different screen densities and gomp_reen
 *      angles.  Postscript data is written to standard out, so use a
 *      command like:
 *
 *                      tocolps blat.rgb | lp
 *
 *                              to actually print a picture (portrait).
 *
 *                      tocolpsl blat.rgb | lp
 *
 *                              to actually print a picture (landscape).
 *
 *   compile with:
 *      cc -o tocolps -O tocolps.c -lc_s -limage -lm
 *      ln tocolps tocolpsl
 *   to create portrait and landscape versions.  Note that the landscape
 *   version can be accessed from tocolps by the -l switch.  This also swaps
 *   the values of maxxsize and maxysize.
 *
 *   This version also now uses different screen angles for each color,
 *   uses a default screen density of 50, and has an improved spot function.
 *   Also, screen angles and frequencies can be changed by tches.
 *
 *   Improvements and bug fixes will be gratefully accepted.  Send to the
 *   addresses below.  Complaints will be forwarded to /dev/null.
 *
 *
 *                                 Steven H. Izen, 4/26/90
 *
 *      tops.c by                       Paul Haeberli - 1988
 *
 *      Adapted from P. H.'s tops.c to work with color
 *      January, April,  1990.
 *        by Steven H. Izen, Dept. of Math. & Stat.,
 *           Case Western Reserve University, Cleveland, OH 44106
 *           izen@cwru.cwru.edu or steve@pitacat.math.cwru.edu or
 *           steve@izen386.math.cwru.edu
 */

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <tcl.h>

#include "gomcast.h"
#include "gomimage.h"
#include "printmsg.h"
#include "tclutils.h"
 
#include "stdafx.h"

static int   PShi[256], PSlow[256];
static int   PSlandscape_flag = 0;
static int   PSpos = 0;
static FILE *Hardcopy_p;


static int tops(    int   , int   , unsigned const char * ,
                float , float ,
                float , float ,
                float , float ,
                float , float ,
                float , float , float , int );

static int PSmaketables(void);
static int PSpsputchar(int);

/***************************************************************************/
int gomp_2ColorPs(const char *FileName , int xsize , int ysize , 
                unsigned const char *pixels   , const char *InputLine)
/***************************************************************************/
{
    static float  pixperinch, maxxsize, maxysize, temp;
    static float  cyanscreendensity, cyanscreenangle;
    static float  magentascreendensity, magentascreenangle;
    static float  yellowscreendensity, yellowscreenangle;
    static float  blackscreendensity, blackscreenangle;
    static int    bitsper, i;
    static int    listArgc;
    const char **listArgv;

    Hardcopy_p = fopen(FileName , "w");
    if(Hardcopy_p == NULL) {
        gomp_PrintERROR("Can't open output file");
        return(1);
    }

/* 
   if(argc<2) {
   fprintf(stderr,"usage: %s inimage [-b bitsperpixel]\n",argv[0]);
   fprintf(stderr,"                    [-l ] # landscape mode\n");
   fprintf(stderr,"                    [-k blackscreendensity]\n");
   fprintf(stderr,"                    [-K blackscreenangle]\n");
   fprintf(stderr,"                    [-c cyanscreendensity]\n");
   fprintf(stderr,"                    [-C cyanscreenangle]\n");
   fprintf(stderr,"                    [-g magentascreendensity]\n");
   fprintf(stderr,"                    [-G magentascreenangle]\n");
   fprintf(stderr,"                    [-y yellowscreendensity]\n");
   fprintf(stderr,"                    [-Y yellowscreenangle]\n");
   fprintf(stderr,"                    [-p pixelsperinch]\n");
   fprintf(stderr,"                    [-m maxxinches maxyinches]\n");
   exit(1);
   }
   if (strcmp(argv[0],"tocolpsl")==0)
   PSlandscape_flag = 1;

   image = iopen(argv[1],"r");
   if(!image) {
   fprintf(stderr,"%s: can't open input image file\n",argv[0]);
   exit(1);
   }
*/

    cyanscreendensity     = 50.0;
    cyanscreenangle       = 15.0;
    magentascreendensity  = 50.0;
    magentascreenangle    = 75.0;
    yellowscreendensity   = 50.0;
    yellowscreenangle     = 0.0;
    blackscreendensity    = 50.0;
    blackscreenangle      = 45.0;
    pixperinch            = -1.0;
    maxxsize              = 528.0;
    maxysize              = 700.0;
    bitsper               = 8;
    temp                  = 0.0;
    PSlandscape_flag   = 0;

/* check to see if there are some extra paramaters */

    if(Tcl_SplitList(gomp_GetTclInterp() , InputLine , 
                     &listArgc , &listArgv) != TCL_OK) {

        gomp_PrintERROR("can't parse the input to the PostScrip parser");
        return(1);
    }

    for(i = 0; i < listArgc; i++) {

        if(listArgv[i][0] == '-') {
            switch(listArgv[i][1]) {
            case 'l':
                PSlandscape_flag = 1;
                break;
            case 'k':
                blackscreendensity = atof(&listArgv[i][2]);
                break;
            case 'K':
                blackscreenangle = atof(&listArgv[i][2]);
                break;
            case 'c':
                cyanscreendensity = atof(&listArgv[i][2]);
                break;
            case 'C':
                cyanscreenangle = atof(&listArgv[i][2]);
                break;
            case 'g':
                magentascreendensity = atof(&listArgv[i][2]);
                break;
            case 'G':
                magentascreenangle = atof(&listArgv[i][2]);
                break;
            case 'y':
                yellowscreendensity = atof(&listArgv[i][2]);
                break;
            case 'Y':
                yellowscreenangle = atof(&listArgv[i][2]);
                break;
            case 'p':
                pixperinch = atof(&listArgv[i][2]);
                break;
            }
        }
    }

    Tcl_Free((char *)CONST_CAST(char **, listArgv));

    if (PSlandscape_flag)       /* swap roles of x and y */
    {
        temp = maxxsize;
        maxxsize = maxysize;
        maxysize = temp;
    }

    (void) tops(xsize , ysize , pixels ,
                cyanscreendensity,cyanscreenangle,
                magentascreendensity,magentascreenangle,
                yellowscreendensity,yellowscreenangle,
                blackscreendensity,blackscreenangle,
                pixperinch,maxxsize,maxysize,bitsper);


    return(0);
}
 
/***************************************************************************/
int tops(int xsize , int ysize , unsigned const char *pixels ,
         float cyanscreendensity,   float cyanscreenangle,
         float magentascreendensity,float magentascreenangle,
         float yellowscreendensity, float yellowscreenangle,
         float blackscreendensity,  float blackscreenangle,
         float pixperinch, float maxxsize, float maxysize, int bitsper)
/***************************************************************************/
{
    static  int   x, y, val, plane;
    static  int   picstrlen;
    static  float doscale, ppiscale;

    PSmaketables();

    PSpos       = 0;
    picstrlen = xsize * bitsper;
    picstrlen = (picstrlen+7)/8;

    if(ysize/(float)xsize < maxysize/maxxsize)
        doscale = maxxsize/xsize;
    else
        doscale = maxysize/ysize;

    if(pixperinch > 0.0) {
        ppiscale = 72.0/pixperinch;
        if(ppiscale<doscale)
            doscale = ppiscale;
        else {
            gomp_PrintERROR("Can't fit image into print area, increase print area with -m option\n");
            return(1);
        }
    }
 
/* put out the header */
    fprintf(Hardcopy_p,"%c",'%');
    fprintf(Hardcopy_p,"!");
    fprintf(Hardcopy_p,"P");
    fprintf(Hardcopy_p,"S");
    fprintf(Hardcopy_p,"\n");
    fprintf(Hardcopy_p,"%c%c Output by tocolps- Steve Izen's hack to tops.c\n",'%','%');
    fprintf(Hardcopy_p,"%c%c Adapted to gOpenMol by Leif Laaksonen CSC 1995\n",'%','%');
    fprintf(Hardcopy_p,"initgraphics\n");

/* define the half tone screen */

    fprintf(Hardcopy_p,"%f %f \n", cyanscreendensity,cyanscreenangle);
    fprintf(Hardcopy_p,"{ abs exch abs 2 copy add 1 gt\n");
    fprintf(Hardcopy_p,"{ 1 sub dup mul exch 1 sub dup mul add 1 sub }\n");
    fprintf(Hardcopy_p,"{ dup mul exch dup mul add 1 exch sub} ifelse }\n");
    fprintf(Hardcopy_p,"%f %f \n", magentascreendensity,magentascreenangle);
    fprintf(Hardcopy_p,"{ abs exch abs 2 copy add 1 gt\n");
    fprintf(Hardcopy_p,"{ 1 sub dup mul exch 1 sub dup mul add 1 sub }\n");
    fprintf(Hardcopy_p,"{ dup mul exch dup mul add 1 exch sub} ifelse }\n");
    fprintf(Hardcopy_p,"%f %f \n", yellowscreendensity,yellowscreenangle);
    fprintf(Hardcopy_p,"{ abs exch abs 2 copy add 1 gt\n");
    fprintf(Hardcopy_p,"{ 1 sub dup mul exch 1 sub dup mul add 1 sub }\n");
    fprintf(Hardcopy_p,"{ dup mul exch dup mul add 1 exch sub} ifelse }\n");
    fprintf(Hardcopy_p,"%f %f \n", blackscreendensity,blackscreenangle);
    fprintf(Hardcopy_p,"{ abs exch abs 2 copy add 1 gt\n");
    fprintf(Hardcopy_p,"{ 1 sub dup mul exch 1 sub dup mul add 1 sub }\n");
    fprintf(Hardcopy_p,"{ dup mul exch dup mul add 1 exch sub} ifelse }\n setcolorscreen\n");
 
 
/* allocate the pixel buffer */
    fprintf(Hardcopy_p,"/rpicstr %d string def\n",picstrlen);
    fprintf(Hardcopy_p,"/gpicstr %d string def\n",picstrlen);
    fprintf(Hardcopy_p,"/bpicstr %d string def\n",picstrlen);

/* rotate image if in landscape mode */
    if (PSlandscape_flag)
        fprintf(Hardcopy_p,"0 792 translate -90 rotate\n");


/* do the proper image translation and gomp_aling */
    fprintf(Hardcopy_p,"45 %f translate\n",maxysize+50.0);
    fprintf(Hardcopy_p,"%f %f scale\n",doscale*xsize, -doscale*ysize);
    fprintf(Hardcopy_p,"%d %d %d\n",xsize,ysize,bitsper);
    fprintf(Hardcopy_p,"[%d 0 0 -%d 0 %d]\n",xsize,ysize,ysize);
    fprintf(Hardcopy_p,"{ currentfile\n rpicstr readhexstring pop}");
    fprintf(Hardcopy_p,"{ currentfile\n gpicstr readhexstring pop}");
    fprintf(Hardcopy_p,"{ currentfile\n bpicstr readhexstring pop}");
    fprintf(Hardcopy_p,"true 3 colorimage\n");
 
/* send out the picure */
    for( y      = 0; y     <     ysize ; y++ ) {
        for( x     = 0; x     <         3 ; x++)  {
            for (plane = 0; plane < 3 * xsize ; plane += 3) {
                switch(bitsper) {
                case 8:                        
                    switch(x) {
                    case 0:
/* red */
                        val = (int)pixels[    plane + 3 * y * xsize];
                        if(val > 255 || val < 0)
                            printf("bad poop\n");
                        PSpsputchar(PShi[val]);
                        PSpsputchar(PSlow[val]);
                        break;

                    case 1:
/* green */
                        val = (int)pixels[1 + plane + 3 * y * xsize];
                        if(val > 255 || val < 0)
                            printf("bad poop\n");
                        PSpsputchar(PShi[val]);
                        PSpsputchar(PSlow[val]);
                        break;

                    case 2:
/* blue */
                        val = (int)pixels[2 + plane + 3 * y * xsize];
                        if(val > 255 || val < 0)
                            printf("bad poop\n");
                        PSpsputchar(PShi[val]);
                        PSpsputchar(PSlow[val]);
                        break;
                    }
                    break;
                default:
                    fprintf(stderr,"bits per pixel must be a power of 2!!\n");
                    exit(1);
                }
            }
        }
    }
    fprintf(Hardcopy_p,"\nshowpage\n");
    fclose(Hardcopy_p);
    return (0);
}
 
/***************************************************************************/
int PSmaketables()
/***************************************************************************/
{
    register int i;
    static   int Switch = 0;

    if(Switch) return (0);
 
    for(i=0; i<256; i++) {
        PShi[i] = "0123456789abcdef"[i>>4];
        PSlow[i] = "0123456789abcdef"[i&0xf];
    }

    Switch = 1;

    return(0);
}
 
 
/***************************************************************************/
int PSpsputchar(int c)
/***************************************************************************/
{
    fprintf(Hardcopy_p,"%c",c);
    if(++PSpos == 50) {
        fprintf(Hardcopy_p,"\n");
        PSpos = 0;
    }

    return(0);
}
 

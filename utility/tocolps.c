 
/*      tocolps -
 *              Convert a color  image to color PostScript.
 *      This knows how to print images with either 1, 2, 4, or 8 bits per
 *      pixel, and how to generate different screen densities and screen
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
#include <gl/image.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
 
int hi[256], low[256];
int reverse_flag = 0;
int landscape_flag = 0;
short buf[4096];
 
main(argc,argv)
int argc;
char **argv;
{
    float pixperinch, maxxsize, maxysize, temp;
    float cyanscreendensity, cyanscreenangle;
    float magentascreendensity, magentascreenangle;
    float yellowscreendensity, yellowscreenangle;
    float blackscreendensity, blackscreenangle;
    int bitsper, i;
    IMAGE *image;
 
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
      landscape_flag = 1;
    image = iopen(argv[1],"r");
    if(!image) {
	fprintf(stderr,"%s: can't open input image file\n",argv[0]);
	exit(1);
    }
    cyanscreendensity = 50.0;
    cyanscreenangle = 15.0;
    magentascreendensity = 50.0;
    magentascreenangle = 75.0;
    yellowscreendensity = 50.0;
    yellowscreenangle = 0.0;
    blackscreendensity = 50.0;
    blackscreenangle = 45.0;
    pixperinch = -1.0;
    maxxsize = 528.0;
    maxysize = 700.0;
    bitsper = 8;
    for(i=2; i<argc; i++) {
	if(argv[i][0] == '-') {
	    switch(argv[i][1]) {
		case 'b':
		    i++;
		    bitsper = atoi(argv[i]);
		    switch(bitsper) {
			case 1:
			case 2:
			case 4:
			case 8:
			    break;
		 	default:
			    fprintf(stderr,"tops: bits per pixel must be 1, 2, 4, or 8\n");
			    exit(1);
		    }
		    break;
		case 'l':
		    landscape_flag = 1;
		    break;
		case 'k':
		    i++;
		    blackscreendensity = atof(argv[i]);
		    break;
		case 'K':
		    i++;
		    blackscreenangle = atof(argv[i]);
		    break;
		case 'c':
		    i++;
		    cyanscreendensity = atof(argv[i]);
		    break;
		case 'C':
		    i++;
		    cyanscreenangle = atof(argv[i]);
		    break;
		case 'g':
		    i++;
		    magentascreendensity = atof(argv[i]);
		    break;
		case 'G':
		    i++;
		    magentascreenangle = atof(argv[i]);
		    break;
		case 'y':
		    i++;
		    yellowscreendensity = atof(argv[i]);
		    break;
		case 'Y':
		    i++;
		    yellowscreenangle = atof(argv[i]);
		    break;
		case 'p':
		    i++;
		    pixperinch = atof(argv[i]);
		    break;
		case 'm':
		    i++;
		    maxxsize = 72.0*atof(argv[i]);
		    i++;
		    maxysize = 72.0*atof(argv[i]);
		    break;
		case 'r':
		    reverse_flag = 1; /* not implemented yet */
		    break;
	    }
	}
    }
    if (landscape_flag)       /* swap roles of x and y */
      {
	temp = maxxsize;
	maxxsize = maxysize;
	maxysize = temp;
      }
    tops(image,
	 cyanscreendensity,cyanscreenangle,
	 magentascreendensity,magentascreenangle,
	 yellowscreendensity,yellowscreenangle,
	 blackscreendensity,blackscreenangle,
	 pixperinch,maxxsize,maxysize,bitsper);
}
 
tops(image,
	 cyanscreendensity,cyanscreenangle,
	 magentascreendensity,magentascreenangle,
	 yellowscreendensity,yellowscreenangle,
	 blackscreendensity,blackscreenangle,
	 pixperinch,maxxsize,maxysize,bitsper)
 
float pixperinch, maxxsize, maxysize;
float cyanscreendensity, cyanscreenangle;
float magentascreendensity, magentascreenangle;
float yellowscreendensity, yellowscreenangle;
float blackscreendensity, blackscreenangle;
int bitsper;
IMAGE *image;
 
{
    register int x, y, n, i, val, plane;
    int picstrlen, xsize, ysize;
    float doscale, ppiscale;
 
    xsize = image->xsize;
    ysize = image->ysize;
    maketables();
 
    picstrlen = xsize*bitsper;
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
        fprintf(stderr,"tops: can't fit image into print area, increase print area with -m option\n");
	    exit(1);
	}
    }
 
/* put out the header */
    putchar('%');
    putchar('!');
    putchar('P');
    putchar('S');
    putchar('\n');
    printf("%% Output by tocolps- Steve Izen's hack to tops.c\n");
    printf("initgraphics\n");
 
/* define the half tone screen */
 
    printf("%f %f \n", cyanscreendensity,cyanscreenangle);
    printf("{ abs exch abs 2 copy add 1 gt\n");
    printf("{ 1 sub dup mul exch 1 sub dup mul add 1 sub }\n");
    printf("{ dup mul exch dup mul add 1 exch sub} ifelse }\n");
    printf("%f %f \n", magentascreendensity,magentascreenangle);
    printf("{ abs exch abs 2 copy add 1 gt\n");
    printf("{ 1 sub dup mul exch 1 sub dup mul add 1 sub }\n");
    printf("{ dup mul exch dup mul add 1 exch sub} ifelse }\n");
    printf("%f %f \n", yellowscreendensity,yellowscreenangle);
    printf("{ abs exch abs 2 copy add 1 gt\n");
    printf("{ 1 sub dup mul exch 1 sub dup mul add 1 sub }\n");
    printf("{ dup mul exch dup mul add 1 exch sub} ifelse }\n");
    printf("%f %f \n", blackscreendensity,blackscreenangle);
    printf("{ abs exch abs 2 copy add 1 gt\n");
    printf("{ 1 sub dup mul exch 1 sub dup mul add 1 sub }\n");
    printf("{ dup mul exch dup mul add 1 exch sub} ifelse }\n setcolorscreen\n");
 
 
/* allocate the pixel buffer */
    printf("/rpicstr %d string def\n",picstrlen);
    printf("/gpicstr %d string def\n",picstrlen);
    printf("/bpicstr %d string def\n",picstrlen);
 
/* rotate image if in landscape mode */
    if (landscape_flag)
      printf("0 792 translate -90 rotate\n");
 
 
/* do the proper image translation and scaling */
    printf("45 %f translate\n",maxysize+50.0);
    printf("%f %f scale\n",doscale*xsize,-doscale*ysize);
    printf("%d %d %d\n",xsize,ysize,bitsper);
    printf("[%d 0 0 -%d 0 %d]\n",xsize,ysize,ysize);
    printf("{ currentfile\n rpicstr readhexstring pop}");
    printf("{ currentfile\n gpicstr readhexstring pop}");
    printf("{ currentfile\n bpicstr readhexstring pop}");
    printf("true 3 colorimage\n");
 
/* send out the picure */
    for( y=0; y<ysize; y++ ) {
      for (plane=0;plane<3;++plane) {
	getrow(image,buf,y,plane);
	switch(bitsper) {
	    case 1:
		x=0;
		for(n=2*picstrlen; n--; ) {
		    val = 0;
		    for(i=0; i<4; i++) {
			val <<= 1;
			val |= (buf[x]&0x80) >> 7;
			x++;
		    }
		    psputchar("0123456789abcdef"[val]);
		}
		break;
	    case 2:
		x=0;
		for(n=2*picstrlen; n--; ) {
		    val = 0;
		    for(i=0; i<2; i++) {
			val <<= 2;
			val |= (buf[x]&0xc0) >> 6;
			x++;
		    }
		    psputchar("0123456789abcdef"[val]);
		}
		break;
	    case 4:
		x=0;
		for(n=2*picstrlen; n--; ) {
		    val = (buf[x]&0xf0) >> 4;
		    x++;
		    psputchar("0123456789abcdef"[val]);
		}
		break;
	    case 8:
		x=0;
		for(n=2*picstrlen; n--; ) {
		    val = buf[x];
		    if(val > 255)
			fprintf(stderr,"bad poop\n");
		    x++;
		    n--;
		    psputchar(hi[val]);
		    psputchar(low[val]);
		}
		break;
	    default:
		fprintf(stderr,"bits per pixel must be a power of 2!!\n");
		exit(1);
	}
    }
    }
    printf("\nshowpage\n");
}
 
maketables()
{
    register int i;
 
    for(i=0; i<256; i++) {
	hi[i] = "0123456789abcdef"[i>>4];
	low[i] = "0123456789abcdef"[i&0xf];
    }
}
 
static int pos = 0;
 
psputchar(c)
int c;
{
    putchar(c);
    if(++pos == 50) {
	putchar('\n');
	pos = 0;
    }
}
 

/*
 *	tomac -
 *		Convert an IRIS image file to macpaint format.
 *
 *	to compile:
 *		cc tomac.c -o tomac -limage -I/usr/include/gl
 *
 *				Paul Haeberli - 1989
 *
 */
#include <gl/image.h>
 
#define HEADERSIZE	512
 
#define MAXXSIZE	576
#define MAXYSIZE	720
 
short sbuf[4096];
unsigned char ibits[MAXXSIZE/8+20];
unsigned char pbits[MAXXSIZE/8+20];
char header[HEADERSIZE];
 
main(argc,argv)
int argc;
char *argv[];
{
    if(argc<3) {
	fprintf(stderr,"usage: tomac image.bw image.mac\n");
	exit(1);
    }
    writemac(argv[1],argv[2]);
}
 
writemac(iname,oname)
char *iname, *oname;
{
    FILE *outf;
    IMAGE *iimage;
    int i, y, n;
    int xmargin, ymargin;
    int xsize, ysize;
 
    iimage = iopen(iname,"r");
    if(!iimage) {
	fprintf(stderr,"tomac: can't open input file %s\n",iname);
	exit(1);
    }
    outf = fopen(oname,"w");
    if(!outf) {
	fprintf(stderr,"tomac: can't open output file %s\n",iname);
	exit(1);
    }
    xsize = iimage->xsize;
    ysize = iimage->ysize;
    xmargin = (MAXXSIZE-xsize)/2.0;
    if(xmargin<0)
	xmargin = 0;
    ymargin = (MAXYSIZE-ysize)/2.0;
    if(ymargin<0)
	ymargin = 0;
    for (i=0; i<HEADERSIZE; i++)
	fputc(0,outf);
    setrow(sbuf,sbuf,255,MAXXSIZE);
    for(y=0; y<MAXYSIZE; y++) {
 	if(y>ymargin && y<(ymargin+ysize))
	    getrow(iimage,sbuf+xmargin,ysize-1-(y-ymargin),0);
 	else
	    setrow(sbuf,255,MAXXSIZE);
	rowtobits(sbuf,ibits,MAXXSIZE);
	n = packbits(ibits,pbits,MAXXSIZE);
  	fwrite(pbits,n,1,outf);
    }
    iclose(iimage);
    fclose(outf);
    return 1;
}
 
packbits(ibits,pbits,nbits)
unsigned char *ibits, *pbits;
int nbits;
{
    int bytes;
    unsigned char *sptr;
    unsigned char *ibitsend;
    unsigned char *optr = pbits;
    int nbytes, todo, cc, count;
 
    nbytes = ((nbits-1)/8)+1;
    ibitsend = ibits+nbytes;
    while(ibits<ibitsend) {
	sptr = ibits;
	ibits += 2;
	while((ibits<ibitsend)&&((ibits[-2]!=ibits[-1])||(ibits[-1]!=ibits[0])))
	    ibits++;
 	if(ibits != ibitsend) {
	    ibits -= 2;
	}
	count = ibits-sptr;
	while(count) {
	    todo = count>127 ? 127:count;
	    count -= todo;
	    *optr++ = todo-1;
	    while(todo--)
		*optr++ = *sptr++;
	}
	if(ibits == ibitsend)
	    break;
	sptr = ibits;
	cc = *ibits++;
	while( (ibits<ibitsend) && (*ibits == cc) )
	    ibits++;
	count = ibits-sptr;
	while(count) {
	    todo = count>128 ? 128:count;
	    count -= todo;
	    *optr++ = 257-todo;
	    *optr++ = cc;
	}
    }
    return optr - pbits;
}
 
/*
 *	row -
 *		support for operations on image rows.
 *
 */
zerorow(sptr,n)
short *sptr;
int n;
{
    bzero(sptr,n*sizeof(short));
}
 
setrow(sptr,val,n)
short *sptr;
int val, n;
{
    if(val==0)
	zerorow(sptr,n);
    else {
	while(n>=8) {
	    sptr[0] = val;
	    sptr[1] = val;
	    sptr[2] = val;
	    sptr[3] = val;
	    sptr[4] = val;
	    sptr[5] = val;
	    sptr[6] = val;
	    sptr[7] = val;
	    sptr += 8;
	    n -= 8;
	}
	while(n--)
	    *sptr++ = val;
    }
}
 
bitstorow(bits,sbuf,n)
unsigned char *bits;
short *sbuf;
int n;
{
    int i, val, nbytes;
 
    nbytes = ((n-1)/8)+1;
    for(i = 0; i<nbytes; i++ ) {
	val = *bits++;
	if(val&0x80)
	    sbuf[0] = 0;
	else
	    sbuf[0] = 255;
	if(val&0x40)
	    sbuf[1] = 0;
	else
	    sbuf[1] = 255;
	if(val&0x20)
	    sbuf[2] = 0;
	else
	    sbuf[2] = 255;
	if(val&0x10)
	    sbuf[3] = 0;
	else
	    sbuf[3] = 255;
	if(val&0x08)
	    sbuf[4] = 0;
	else
	    sbuf[4] = 255;
	if(val&0x04)
	    sbuf[5] = 0;
	else
	    sbuf[5] = 255;
	if(val&0x02)
	    sbuf[6] = 0;
	else
	    sbuf[6] = 255;
	if(val&0x01)
	    sbuf[7] = 0;
	else
	    sbuf[7] = 255;
	sbuf += 8;
    }
}
 
rowtobits(sbuf,bits,n)
short *sbuf;
unsigned char *bits;
int n;
{
    int i, val, nbytes, thresh;
 
    nbytes = ((n-1)/8)+1;
    thresh = 128;
    for(i = 0; i<nbytes; i++) {
	val = 0;
	if(sbuf[0]<thresh)
	    val |= 0x80;
	if(sbuf[1]<thresh)
	    val |= 0x40;
	if(sbuf[2]<thresh)
	    val |= 0x20;
	if(sbuf[3]<thresh)
	    val |= 0x10;
	if(sbuf[4]<thresh)
	    val |= 0x08;
	if(sbuf[5]<thresh)
	    val |= 0x04;
	if(sbuf[6]<thresh)
	    val |= 0x02;
	if(sbuf[7]<thresh)
	    val |= 0x01;
	sbuf += 8;
	*bits++ = val;
    }
}


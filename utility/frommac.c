
/*
 *	frommac -
 *		Convert macpaint image files to IRIS image files.
 *
 *				Paul Haeberli - 1989
 */
#include "image.h"
 
#define MAXXSIZE	576
#define MAXYSIZE	720
 
short sbuf[4096];
unsigned char cbuf[72+256];
 
main(argc,argv)
int argc;
char *argv[];
{
    if(argc<3) {
	fprintf(stderr,"usage: frommac image.mac image.bw\n");
	exit(1);
    }
    if(readmac(argv[1],argv[2],512))
	exit(0);
    else if(readmac(argv[1],argv[2],512+128))
	exit(0);
    else {
	fprintf(stderr,"frommac: bad macpaint file %s\n",argv[1]);
	exit(1);
    }
}
 
readmac(iname,oname,offset)
char *iname, *oname;
int offset;
{
    FILE *inf;
    IMAGE *oimage;
    int i, y;
 
    inf = fopen(iname,"r");
    if(!inf) {
	fprintf(stderr,"frommac: can't open input file %s\n",iname);
	exit(1);
    }
    oimage = iopen(oname,"w",RLE(1),2,MAXXSIZE,MAXYSIZE,1);
    for (i=0; i<offset; i++)
	getc(inf);
    for(y=0; y<MAXYSIZE; y++) {
	if(!readline(inf)) {
	    iclose(oimage);
	    fclose(inf);
	    return 0;
	    break;
	}
	bitstorow(cbuf,sbuf,MAXXSIZE);
	putrow(oimage,sbuf,MAXYSIZE-1-y,0);
    }
    iclose(oimage);
    fclose(inf);
    fprintf(stderr,"head size was %d\n",offset);
    return 1;
}
 
readline(inf)
FILE *inf;
{
    int pos, cnt, val;
 
    pos = 0;
    while(pos < 72) {
	cnt = getc(inf);
	if((cnt&0x80)==0) {
	    cnt++;
	    while(cnt--)
		cbuf[pos++] = getc(inf);
	} else {
	    cnt = 257-cnt;
	    val = getc(inf);
	    while (cnt--)
		cbuf[pos++] = val;
	}
    }
    if(pos==72)
	return 1;
    else
	return 0;
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


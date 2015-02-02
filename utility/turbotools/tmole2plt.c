/* 
 *                            tmole2plt 
 *
 *                   Written by Jonas Juselius,
 *             Department of Chemistry, University of Helsinki,
 *                           FINLAND 1998
 * 
 *   This is a small program to convert TURBOMOLE moloch grids 
 *   into a format recongized by gOpenMol. The proram is loosely
 *   based on code written originally by Leif Laaksonen at
 *   Centre for Scientific Computing , ESPOO, FINLAND.
 *  
 *   Usage: $ tmole2plt {infile} {outfile}
 *     ex.: $ tmole2plt density.dat density.plt
 *
 *   Oh, yes... This program is supplied as is, and comes without any kind 
 *   of WARRANTY what so ever! If you have problems, you can try to contact
 *   me by mail (if you can find me that iz :^). 
 *   Also, sorry for the lack of comments ;)
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define IWRITE(value_p)    { iostat = \
                                 fwrite(value_p, sizeof(int), 1L , outfd);\
                                 if(iostat < 1) {\
                     printf("?ERROR - in writing contour file (*)\n");\
                     return(1);}}
#define FWRITE(value_p)    { iostat = \
                                 fwrite(value_p, sizeof(float), 1L , outfd);\
                                 if(iostat < 1) {\
                     printf("?ERROR - in writing contour file (*)\n");\
                     return(1);}}

#define BUFF_LEN   256

int surface;              		/* = 0 , unknown
                                   = 1 , VSS surface
                                   = 2 , orbital- or density surface
                                   = 3 , probe surface
                                */

char   infile[BUFF_LEN];
char   outfile[BUFF_LEN];
FILE  *infd;
FILE  *outfd;
float  xmin, xmax;
float  ymin, ymax;
float  zmin, zmax;
int  xdim, ydim, zdim;
float  *databuf;

int iostat;

int write_plt(FILE *outfd);
int read_tmole(FILE *inffd);

int      rank;
int      surface;
float    au2a  = 0.52917715;

main(int argc , char **argv)
{
	/* sanity check... */
	if (argc != 3) {
		fprintf(stderr, "*** tmole2plt by j0nas, 1998\n");
		fprintf(stderr, "usage: tmole2plt infile outfile\n");
		exit(1);
	}

	/* input file */
	strncpy(infile, argv[1], BUFF_LEN);
	printf("Input file:     '%s'\n", infile);

	/* output file */
	strncpy(outfile, argv[2], BUFF_LEN);
	printf("Output file:    '%s'\n", outfile);

	infd = fopen(infile, "r");

	if(infd == NULL) {
		printf("?: can't open input file: '%s'\n", infile);
		exit(1);
	}

	(void) read_tmole(infd);

	printf("Xpoints: %d, YPoints: %d, Zpoints: %d\n", xdim, ydim, zdim);
	printf("Xmin: %f, Xmax: %f\nYmin: %f, Ymax: %f\nZmin: %f, Zmax: %f\n",
			xmin, xmax, ymin, ymax, zmin, zmax);

	outfd = fopen(outfile,"wb");

	if (outfd == NULL) {
		printf("?: can't open output file : '%s'\n", outfile);
		exit(1);
	}

	(void) write_plt(outfd);
	free(databuf);

	return(0);
}

int write_plt(FILE *outfd)
{
	int i;
	float max = -1.0e10;
	float min =  1.0e10;

	/* print out number of points in x- , y- and z-directions */

	IWRITE(&rank);
	IWRITE(&surface);
	IWRITE(&zdim);
	IWRITE(&ydim);
	IWRITE(&xdim);

	FWRITE(&zmin);
	FWRITE(&zmax);
	FWRITE(&ymin);
	FWRITE(&ymax);
	FWRITE(&xmin);
	FWRITE(&xmax);

	for (i=0; i < xdim*ydim*zdim; i++) {
		if (databuf[i] < min) 
			min = databuf[i];
		if (databuf[i] > max) 
			max = databuf[i];

		FWRITE(&databuf[i]);
	}      

	printf("Max value: %f \n", max);
	printf("Min value: %f \n", min);
	printf("! Done.\n");

	fclose(outfd);
	return 1;
}


int read_tmole(FILE *infd)
{
	int i;
	float delta;
	char line[256], str[256], *cp;

	rank=3;
	surface=200;

	while (1) {
		fgets(line, 256, infd);
		sscanf(line,"%s", str);

		if (!strcmp(str, "$grid1")) {
			sscanf(line, "%*s %*s %f %*s %f %*s %d", 
					&xmin, &delta, &xdim);
			xmin*=au2a;
			delta*=au2a;
			xmax=delta*(xdim-1)+xmin;
		}
		else if (!strcmp(str, "$grid2")) {
			sscanf(line, "%*s %*s %f %*s %f %*s %d",
					&ymin, &delta, &ydim);
			ymin*=au2a;
			delta*=au2a;
			ymax=delta*(ydim-1)+ymin;
		}
		else if (!strcmp(str, "$grid3")) {
			sscanf(line, "%*s %*s %f %*s %f %*s %d",
					&zmin, &delta, &zdim);
			zmin*=au2a;
			delta*=au2a;
			zmax=delta*(zdim-1)+zmin;
		}
		else if (!strcmp(str, "$plotdata")) 
			break;
	}

	databuf = (float *) malloc(xdim*ydim*zdim*sizeof(float));

	if (databuf == NULL) {
		printf("$ERROR - can't allocate memory\n");
		exit(1);
	}
	
	fscanf(infd, "%15s", str);

	if ( (cp=strchr(str, 'D')) != NULL ) {
		*cp='e';
		databuf[0] = atof(str);
		for (i=1; i < xdim*ydim*zdim; i++) {
			fscanf(infd,"%15s",str);
			cp=strchr(str, 'D');
			*cp='e';
			databuf[i] = atof(str);
		}      
	}
	else if ( (cp=strchr(str, 'd')) != NULL ) {
		*cp='e';
		databuf[0] = atof(str);
		for (i=1; i < xdim*ydim*zdim; i++) {
			fscanf(infd,"%15s",str);
			cp=strchr(str, 'd');
			*cp='e';
			databuf[i] = atof(str);
		}      
	}
	else {
		databuf[0] = atof(str);
		for (i=1; i < xdim*ydim*zdim; i++) {
			fscanf(infd,"%15s",str);
			databuf[i] = atof(str);
		}      
	}

	fclose(infd);

	return(0);
}


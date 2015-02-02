/*
                   PROGRAM VSSMOD    
C---ORIGINAL PROGRAMME WRITTEN BY C. GIESSNER-PRETTRE.                      
C
C--MODIFED BY  SAL FEB 1986
C--Hacked further by LUL 1989 (all CDL parts removed)
C
C   This version of the main block is considerably different to the 
C   original VSS programme. Code for the production of line-printer
C   plots has been removed. The remaining code has been relaid and 
C   modified to give potentials at user defined intervals along 
C   all the axes.
C
C
C-Modified by SAL JULY 1986 to allow finer grid and adjustable array
C- for the potentials
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>

#define LINE_LEN 80
#define BUFF_LEN 500
#define BUFFER_LIM 500
#define CONTOUR_FILE "vsscont.plt"
#define SEP_STRING ","

#define FWRITE(value_p , size)    { Items = \
                                 fwrite(value_p, size , 1L , VSS_p);\
                                 if(Items < 1) {\
                     printf("?ERROR - in writing contour file (*)\n");\
                     return(1);}}

#define FWRITEN(value_p , num , size) { Items = \
                                fwrite(value_p, size , num , VSS_p);\
                   if(Items < 1) {\
                     printf("?ERROR - in writing contour file (**)\n");\
                     return(1);}}
void uotm(float , float , float , float *);                                             
void uota1(float, float);
void uota2(float, float);
void ParseVSSlimits(int *, float *, float *);


      float *ar3pot; /* potential array */
      float *x,*y,*z; /* x,y and z-coordinate pointers */
      float *q;       /* gross atomic population */
      float *zet;     /* screening constant pointer */
      int    *iz;      /* number of valence electrons */
      float  vs;
      int    nat;      /* total number of atoms */
      int    nap;      /* number of heavy atoms (?) */
      int    nah;      /* number of hydrogen atoms */
      int npts[3];     /* number of points in x,y and z */
      int BuffMult;    /* number of inner loops in buffer */
      float xmin[3];  /* min value */
      float xmax[3];  /* max value */
      float rinc[3];  /* step length */
      float autang = 0.52917715;   /*             */
      int    iqmap;
      int    OldType = 0; /* if !=0 old type format in input file for
                             grossp and screen */
      char VSS_file[BUFF_LEN] = CONTOUR_FILE;
      char VSS_limits[BUFF_LEN] = "\0";

           char *announce =
"      **********************************************************\n\
      *                                                        *\n\
      *    Electrostatic Potential program                     *\n\
      *    VSS                                                 *\n\
      *                                                        *\n\
      *    Leif Laaksonen    Version 1.0                       *\n\
      *                      1991                              *\n\
      *                                                        *\n\
      **********************************************************\n\n";

main(argc,argv)
    int argc;
    char **argv;
{

      float angtau;
      int iu,nh,i;
      char ctitle[LINE_LEN];

      printf("%s",announce);

/* parse input line */
       if(argc > 1) {
               for(i = 1 ; i < argc ; i++) {
                  if(argv[i][0] == '-') {
                     switch(argv[i][1]) {
                         case 'o':    /* output file name */
                                  strncpy(VSS_file,&argv[i][2],BUFF_LEN);
                                  break;
                         case 'z':    /* old type input format */
                                  OldType = 1;
                                  break;
			 case 'l':    /* points and coordinates limits
                         numx,numy,numz,xmin,xmax,ymin,ymax,zmin,zmax */
                                  strncpy(VSS_limits,&argv[i][2],BUFF_LEN);
                                  break;
                         default:
                                  printf("?ERROR - unknown parameter '-%s'\n",
                                  argv[i][1]);
                                  break;
                      } /* switch */
                   } /* if */
                } /* for */
        } /* if */ 

      angtau=1/autang;

/* read title. */
      gets(ctitle);
      printf("Title: %s \n",ctitle);                                  
      gets(ctitle);
      printf("Title: %s \n",ctitle);                                  
      gets(ctitle);
      printf("Title: %s \n",ctitle);                                  
/* read total number of atoms, number of hydrogens and units flag. */
      gets(ctitle);
      sscanf(ctitle,"%d %d %d",&nat,&nh,&iu);
       nap = nat - nh;                                            

      printf("Number of atoms    : %d \n",nat);
      printf("Number of hydrogens: %d \n",nh);

/* read atom coordinates and number of valence electrons */
      x = (float *) malloc(nat * sizeof(float));
      y = (float *) malloc(nat * sizeof(float));
      z = (float *) malloc(nat * sizeof(float));
      iz= (int *) malloc(nat * sizeof(int));

      if(x == NULL || y == NULL || z == NULL || iz == NULL) {
        printf("?ERROR - memory allocation error \n");
         exit(1);}

      for(i = 0; i < nat ; i++ ) {
      gets(ctitle);
      sscanf(ctitle,"%f %f %f %d",&x[i],&y[i],&z[i],&iz[i]);}

      printf("\n**** Coordinates and number of electrons ****\n");

      for(i = 0; i < nat ; i++) 
      printf("Nr: %d x: %f y: %f z: %f nel: %d \n",i+1,x[i],y[i],z[i],iz[i]);  

/* if necessary convert angstom units to atomic units */
      if(iu == 1) {                                                          
        for(i = 0 ; i < nat ; i++) {                                    
        x[i]=x[i]*angtau;                                                
        y[i]=y[i]*angtau;                                                
        z[i]=z[i]*angtau;}                                                
      }

      q = (float *) malloc(nat * sizeof(float));
      zet = (float *) malloc(nat * sizeof(float));
      if(q == NULL || zet == NULL) {
        printf("?ERROR - memory allocation error \n");
        exit(1);
      }

      if(OldType) {
/* read gross atomic populations */
      for(i = 0 ; i < nat ; i++ ) {
        gets(ctitle);
        sscanf(ctitle,"%f",&q[i]);}

/* read  atomic orbital screening constants. */
      for(i = 0 ; i < nat ; i++ ) {
       gets(ctitle);
       sscanf(ctitle,"%f",&zet[i]);}}
      else {
      for(i = 0 ; i < nat ; i++ ) {
        gets(ctitle);
        sscanf(ctitle,"%f %f",&q[i],&zet[i]);}
      }

/* read type of map */
      gets(ctitle);
      sscanf(ctitle,"%d",&iqmap);  

/* read coordinate bounds and number of intervals along the axes. */
      gets(ctitle);
      sscanf(ctitle,"%f %f %f %f %f %f",
      &xmin[0],&xmax[0],&xmin[1],&xmax[1],&xmin[2],&xmax[2]);

      gets(ctitle);
      sscanf(ctitle,"%d %d %d",
      &npts[0],&npts[1],&npts[2]);

/*    check inner loop (x) */
      BuffMult = BUFFER_LIM / npts[0];
      if(BuffMult < 1) {
       printf("?ERROR - internal buffer (BUFFER_LIM) has to be deeper than number of points in x direction\n");
       exit(1);}

      if(VSS_limits[0] != '\0') ParseVSSlimits(npts,xmin,xmax);

      printf("Xmin: %f , Xmax: %f\n",xmin[0],xmax[0]);
      printf("Ymin: %f , Ymax: %f\n",xmin[1],xmax[1]);
      printf("Zmin: %f , Zmax: %f\n",xmin[2],xmax[2]);
      printf("Xpts: %d , Ypts: %d , Zpts: %d\n",npts[0],npts[1],npts[2]);

/* change to atomic units  */
      for(i=0;i<3;i++) {
       xmax[i] /= autang;
        xmin[i] /= autang;}

/* calculate step-lengths  */
        for(i = 0 ; i < 3 ; i++) {        
          if(npts[i] > 1) 
            rinc[i] = (xmax[i]-xmin[i])/(npts[i]-1);
          else
            rinc[i]=0.0;
	}
/* generate the binary map file and listing file */
      slpmap(); /* go and get it now ....*/ 

      }


slpmap()
{
      static int i,j,k,Loop,Loop1,Items;
      static float xc,yc,zc,pot;
      static float Help1,Help2;
      FILE *VSS_p;
      static float Buffer[BUFFER_LIM];
      static float MaxV = -1.0e10;
      static float MinV =  1.0e10;
      static int TypeOfSurface;

/* print out number of points in x- , y- and z-directions */

      VSS_p = fopen(VSS_file,"wb");
      if(VSS_p == NULL) {
        printf("?ERROR - can't open contour file '%s' \n",VSS_file);
        return;}

      printf("\n==> Writing contour to file '%s' \n",VSS_file);


      i = 3;
      FWRITE(&i , sizeof(int));
      TypeOfSurface = 1;
      FWRITE(&TypeOfSurface , sizeof(int));
      FWRITE(&npts[2] , sizeof(int));
       FWRITE(&npts[1] , sizeof(int));
        FWRITE(&npts[0] , sizeof(int));

      Help1 = xmin[2]*autang;
       Help2 = xmax[2]*autang;
        FWRITE(&Help1 , sizeof(float));
        FWRITE(&Help2 , sizeof(float));
      Help1 = xmin[1]*autang;
       Help2 = xmax[1]*autang;
        FWRITE(&Help1 , sizeof(float));
        FWRITE(&Help2 , sizeof(float));
      Help1 = xmin[0]*autang;
       Help2 = xmax[0]*autang;
        FWRITE(&Help1 , sizeof(float));
        FWRITE(&Help2 , sizeof(float));
/*

      fprintf(VSS_p," 3         ! Rank value\n");
      fprintf(VSS_p," 3         ! Type of surface\n");
      fprintf(VSS_p," %d        ! number of points in the z-direction\n",npts[2]);
      fprintf(VSS_p," %d        ! number of points in the y-direction\n",npts[1]);
      fprintf(VSS_p," %d        ! number of points in the x-direction\n",npts[0]);
      fprintf(VSS_p," %f  %f    ! zmin and zmax\n",xmin[2]*autang,xmax[2]*autang);
      fprintf(VSS_p," %f  %f    ! ymin and ymax\n",xmin[1]*autang,xmax[1]*autang);
      fprintf(VSS_p," %f  %f    ! xmin and xmax\n",xmin[0]*autang,xmax[0]*autang);

*/

      MaxV = -1.0e10;
      MinV =  1.0e10;

      Loop = 0;
       Loop1 = 0;

      for(i = 0 ; i < npts[2] ; i++) { /* loop over z */
            zc = xmin[2] + rinc[2] * i;
         for(j = 0 ; j < npts[1] ; j++) { /* loop over y */
            yc = xmin[1] + rinc[1] * j;
            for(k = 0; k < npts[0] ; k++) { /* loop over x */
            xc = xmin[0] + rinc[0] * k;

            uotm(xc,yc,zc,&pot);

            Help1 = pot * 627.5396;

            if(Help1 < MinV) MinV = Help1;
            if(Help1 > MaxV) MaxV = Help1;

              Buffer[Loop] = Help1;
              Loop++;

/*
            FWRITE(&Help1 , sizeof(float));

            fprintf(VSS_p,"%f \n",pot);
*/
            }   /* x */

            Loop1++;

            if(Loop1 == BuffMult){
              FWRITEN(Buffer , Loop , sizeof(float));
                Loop = 0;
                 Loop1 = 0;}

         }      /* y */
      }         /* z */

      if(Loop)
            FWRITEN(Buffer , Loop , sizeof(float));

      printf("Max value: %f \n",MaxV);
      printf("Min value: %f \n",MinV);
      printf("! Done.\n");

      fclose(VSS_p);
}

/***************************************************************************/
void uotm(float xh, float yh, float zh, float *pot)                                             
/***************************************************************************/
{
      register int i;
      static   float rep,ro;
      static   float xi,yi,zi;
      static   float xf,yf,zf;

      *pot=0.0;
       rep=0.0;
  
      for(i = 0 ; i < nat ; i++) {
       xi=x[i];
        yi=y[i]; 
         zi=z[i];

       xf = xi - xh;
        yf = yi - yh;
         zf = zi - zh;

         ro=sqrt( xf * xf + yf * yf + zf * zf);

        if(ro - 0.1 <= 0.0) {
         ro = 0.0;
         *pot = 0.0;
         return;}
        else {
           if(i-nap <= 0) {                                         
           uota2(ro,zet[i]);                                    
           *pot =  *pot - vs*q[i];}
           else {                                      
           uota1(ro,zet[i]);                                
           *pot =  *pot - vs*q[i]; }                                   
        }
      if(ro  >  0.0)  rep=rep+iz[i]/ro;                                    
     }
        *pot = *pot + rep;
}

void uota2(float r, float zcc)
{                                                    
     static float ro;
     static float dero;

      if(r == 0.0 )   vs=zcc/2.;                                     
      else {
      ro=zcc*r;                                                        
      dero=exp(-2.*ro);                                                        
      vs=(1./r)*(1.-(1.+ro*(1.5+ro*(1.+ro/3.)))*dero);}                
   }

void uota1(float r, float zcc)                                                    
{
      static float ro;
      static float dero;

      if(r == 0.0 ) vs = zcc;                                     
      else {
      ro=zcc*r;                                                     
      dero=exp(-2.*ro);                                             
      vs=(1./r)*(1.-(1.+ro)*dero); }
}

/*************************************************************************/
void ParseVSSlimits(int *GetNpts, float *GetMin, float *GetMax)
/*************************************************************************/
{
     char *find;
     int loop;

     printf(".... > Redefining input limits < ....\n");

     find = strtok(VSS_limits,SEP_STRING);

     if(find ==NULL) return;

     GetNpts[0] = atoi(find);

     for(loop = 1; loop < 9 ; loop++) {

     find = strtok(NULL,SEP_STRING);

       if(find == NULL) return;

       switch(loop) {
       case 1: /* ypts */
              GetNpts[1] = atoi(find);
              break;
       case 2: /* zpts */
              GetNpts[2] = atoi(find);
              break;
       case 3: /* xmin */
              GetMin[0] = atof(find);
              break;
       case 4: /* xmax */
              GetMin[1] = atof(find);
              break;
       case 5: /* ymin */
              GetMin[2] = atof(find);
              break;
       case 6: /* ymax */
              GetMin[3] = atof(find);
              break;
       case 7: /* zmin */
              GetMin[4] = atof(find);
              break;
       case 8: /* zmax */
              GetMin[5] = atof(find);
              break;}
   }
}







/*

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved

Enhancements 2003, 2005 by:
Eero Häkkinen

*/

/*

This FORTAN code originates from Biosym/MSI. The c-hack is a product
of my own imagination.

Leif Laaksonen, CSC (1995)

C**********************************************************************
C
C                       PROPRIETARY AND CONFIDENTIAL
C              Copyright (c) 1994, Molecular Simulations Inc.
C
C  THIS PROGRAM IS AN UNPUBLISHED WORK FULLY PROTECTED BY THE UNITED
C  STATES COPYRIGHT LAWS AND IS NOT TO BE COPIED IN ANY FORM WITHOUT
C  PERMISSION OF THE COPYRIGHT HOLDER, MOLECULAR SIMULATIONS INC.
C  IT IS CONSIDERED PROPRIETARY AND IS NOT TO BE DIVULGED OR USED BY
C  PARTIES WHO HAVE NOT RECEIVED WRITTEN AUTHORIZATION FROM THE OWNERS.
C
C       $Header: /usr/source/c2cap_c16/src/trajconv/RCS/trj2ascii.f,v 1.2 94/11/16 12:20:58 mpinches Exp $
C
C       Change record:
C
C       $Log:   trj2ascii.f,v $
C   Revision 1.2  94/11/16  12:20:58  mpinches
C   Submission Number: 3807 -- Remote Site Checkin
C   
C   Revision 1.1  94/04/08  06:33:36  jeremy
C   Initial revision
C   
C       Born from Polygraf rdtrj
C
C       Purpose:
C       Conversion utility between Polygraf trajectory file and ascii.
C
C**********************************************************************

program trj2ascii_main
      implicit none

c     Local variables
      integer       in, out, dest
      parameter     (in=20)
      parameter     (out=22)
      character*256 infile, outfile
      character*1   yesno
      integer       narg
      logical       exists
      integer       lenin, lenout

c.....Functions
      integer       iargc, leng1

c.....Code.
c
c     get the input file name.
c
      narg=iargc()
      if (narg.lt.1) then
         write(*,*) 'SYNTAX: trj2ascii infile [outfile]'
      elseif (narg.gt.2) then
         write(*,*) 'Maximum of two arguments should be supplied!'
      else

c.......Establish input file
         call getarg(1,infile)
         lenin = leng1(infile)
         open (unit=in,file=infile(1:lenin),err=7777,
 1       status='OLD',form='unformatted')

         if (narg.gt.1) then

c.........Output to file in argument 2

            call getarg(2,outfile)
            lenout = leng1(outfile)
            exists = .false.
            inquire(file=outfile(1:lenout),exist=exists)
            if (exists) then
               write(*,75) outfile(1:lenout)
               write(*,80)
               read(*,84) yesno
               write(*,*)
               if (yesno(1:1).eq.'n'.or.yesno(1:1).eq.'N') go to 9999
               open (unit=out,file=outfile(1:lenout),err=8888,
 1             status='OLD',form='formatted')
            else
               open (unit=out,file=outfile(1:lenout),err=8888,
 1             status='NEW',form='formatted')
            endif
            dest = out
         else

c.........Output to standard output

            dest = 6
         endif
         call convert(in, dest)

c.......Close files
         close (in)
         if (narg.gt.1) close(out)

      endif

      go to 9999
c
c     error opening input file.
c
 7777 continue
      write(*,7780)infile(1:lenin)
      go to 9999
c
c     error opening output file.
c
 8888 continue
      write(*,8890)outfile(1:lenout)
      close(unit = in)
      go to 9999

 9999 continue
      stop

 75   format(/,1x,' The output file ',a,' already exists.')
 80   format(1x,' Do you want to overwrite (y/n) : ',$)
 84   format(a1)
 7780 format(/,1x,' Error opening input file : ',a)
 8890 format(/,1x,' Error opening output file : ',a)
      end
c

      subroutine convert(in, dest)
      implicit none

c.....Arguments
      integer            in     ! [i] input unit
      integer            dest   ! [i] destput unit

c.....Local variables
c
      integer i
      integer j
c
      integer MAXATM
      parameter (MAXATM=20000)
      integer MXFILE
      parameter (MXFILE=24)

      character*4 hdr
      integer icntrl(20)
      integer ntrjti
      character*80 trjtic(10)
      integer neexti
      character*80 eextic(10)
      logical period
      logical molxtl
      logical lcanon
      logical defcel
      logical prtthrm
      logical lnose
      logical lnpecan
      logical ltmpdamp
      LOGICAL l1, l2, l3, l4, l5, l6, l7, l8
      integer natom
      integer nmovatm
      integer movatm1
      integer movatmn
      integer leexti
      character*80 eextit
      integer lparti
      character*80 partit
c
      character*80 yesno
      character*80 infile
      character*80 destfile
      integer      idx
      logical      exists
      integer      len,index
      real         zfrdis,zprfds
c
      real curtim
      integer itstep
      real system
      real avetem
      real dstep
      real firstt
      real finalt
      double precision e
      double precision eb
      double precision et
      double precision ep
      double precision ei
      double precision enb
      double precision eel
      double precision ehb
      double precision ec
      double precision eu
      double precision tint
      double precision tnb
      double precision ea
      double precision eba
      double precision eta
      double precision epa
      double precision eia
      double precision enba
      double precision eela
      double precision ehba
      double precision eca
      double precision eua
      double precision tinta
      double precision tnba
      double precision tote
      double precision totke
      double precision totea
      double precision tkea
      integer  iconmp
      integer imstep
      logical lvelwr
      integer iconfs
      integer icstep
c
      double precision pressura
      double precision vola
      double precision pvtota
      double precision pvkina
      double precision pvpota
      double precision radgyra
c
      double precision signose
      double precision zfrict
      real   zprfrict
      double precision snose,snoseh,ssdot
      real qcanon
      double precision sigdyn(2)
      real gamtmp
c
      double precision tcela
      real s2r(6)
      real s2rdot(6)
c
      integer natmcel
      double precision strsa(6)
      double precision extstrsa ! added after 220 but version called 210
c
      double precision eabtota
      double precision eabvala
      double precision eabelha
      double precision eabnba
      double precision eabmisa
      double precision dltfaba
      double precision expprta
c
      integer  nflusd
      integer  mvatmpfu(MXFILE)
      integer  natmpfu(MXFILE)
      integer  totmov
      integer  mvatmofst(maxatm)
c
      character  decusd(MXFILE)*8
c
      integer version
c
      real x(MAXATM),y(MAXATM),z(MAXATM)
      real velx(MAXATM),vely(MAXATM),velz(MAXATM)
      write(*,88)
 88   format(/,1x,' Working....')
c
      read (in,err=9990,end=9000) hdr,icntrl
      version=icntrl(1)
c
      write (dest,10) hdr
 10   format ('Header:',t10,a4)
      write (dest,20) icntrl
 20   format ('Control:',t10,20i5)
c
      if (version.ge.311) then
         czmc 2 tc looking at thermodyn.f for reason
c     read (in,err=9990,end=9000) ntrjti,(trjtic(i)(1:60),i=1,ntrjti)
         read (in,err=9990,end=9000) ntrjti,(trjtic(i),i=1,ntrjti)
         write (dest,25) ntrjti
         do 12 j=1,ntrjti
            write (dest,26) trjtic(j)(1:60)
 12      continue
      endif
 25   format ('Nremarks:',t10,i5)
 26   format ('REMARK :',1x,a60)
c     
      read (in,err=9990,end=9000) neexti,(eextic(j),j=1,neexti)
c
      write (dest,30) neexti,(eextic(j),j=1,neexti)
 30   format ('Neexti:',t10,i5,/,(t10,a80))
c
      if (version.le.150) then
         period=.false.
         molxtl=.false.
         lcanon=.false.
         defcel=.false.
         prtthrm=.false.
      elseif ( version.lt. 300 ) then
c<mrsp>         read (in,err=9990,end=9000) period,molxtl,lcanon,defcel,
c<mrsp>     $        prtthrm
         read (in,err=9990,end=9000) l1, l2, l3, l4, l5
         period  = (.NOT.(l1 .EQV. .FALSE.))
         molxtl  = (.NOT.(l2 .EQV. .FALSE.))
         lcanon  = (.NOT.(l3 .EQV. .FALSE.))
         defcel  = (.NOT.(l4 .EQV. .FALSE.))
         prtthrm = (.NOT.(l5 .EQV. .FALSE.))
         lnose = .false.
         lnpecan=.false.
         ltmpdamp=.false.
         write (dest,40,err=9995) period,molxtl,lcanon,defcel,prtthrm
 40      format ('Period:',t15,l2,/,
     $        'MolXtl:',t15,l2,/,
     $        'Canonical:',t15,l2,/,
     $        'DefCell:',t15,l2,/,
     $        'PertTheory:',t15,l2)
      else
c<mrsp>         read (in,err=9990,end=9000) period,molxtl,lcanon,defcel,
c<mrsp>     $        prtthrm,lnose,lnpecan,ltmpdamp
         read (in,err=9990,end=9000) l1, l2, l3, l4, l5, l6, l7, l8
         period   = (.NOT.(l1 .EQV. .FALSE.))
         molxtl   = (.NOT.(l2 .EQV. .FALSE.))
         lcanon   = (.NOT.(l3 .EQV. .FALSE.))
         defcel   = (.NOT.(l4 .EQV. .FALSE.))
         prtthrm  = (.NOT.(l5 .EQV. .FALSE.))
         lnose    = (.NOT.(l6 .EQV. .FALSE.))
         lnpecan  = (.NOT.(l7 .EQV. .FALSE.))
         ltmpdamp = (.NOT.(l8 .EQV. .FALSE.))
         write (dest,41,err=9995) period,molxtl,lcanon,defcel,
     $        prtthrm,lnose,lnpecan,ltmpdamp
 41      format ('Period:',t15,l2,/,
     $        'MolXtl:',t15,l2,/,
     $        'Canonical:',t15,l2,/,
     $        'DefCell:',t15,l2,/,
     $        'PertTheory:',t15,l2,/,
     $        'NoseorHoover:',t15,l2,/,
     $        'NpTCanon:',t15,l2,/,
     $        'TempDamping:',t15,l2)
      end if
c     
      if (version.ge.200) then
         read(in,err=9990,end=9000)nflusd
         write(dest,200,err=9995) nflusd
 200     format('Nfileused:',t15,i5)
         if (nflusd.gt.MXFILE) goto 9990
         backspace(in)
         read(in,err=9990,end=9000)nflusd,(mvatmpfu(i),i=1,nflusd),
     $        (natmpfu(i),i=1,nflusd),(decusd(i),i=1,nflusd)
c
         do 300 i = 1,nflusd
            write(dest,220,err=9995) i,mvatmpfu(i),natmpfu(i),decusd(i)
 220        format('Filnum:',3x,i5,1x,'Movatms:',2x,i5,
 1          1x,'Totatms:',2x,i5,1x,'Descriptor:',4x,a8)
 300      continue
c     
          read(in,err=9990,end=9000)totmov
          write(dest,350,err=9995)totmov
 350      format('Totmovatm:',t15,i5)

          if (totmov.gt.MAXATM) goto 9990
          backspace(in)
          read(in,err=9990,end=9000)totmov,(mvatmofst(i),i=1,totmov)
c
          write(dest,400,err=9995)(mvatmofst(i),i=1,totmov)
 400      format('Movatmofst:',t15,10i5,/,(t15,10i5) )
c
          nmovatm = totmov
c
       else
          read (in,err=9990,end=9000) natom,nmovatm,movatm1,movatmn
c
          write (dest,500,err=9995) natom,nmovatm,movatm1,movatmn
 500      format ('Natom:',t10,i5,/,
     $         'NMovatm:',t10,i5,/,
     $         'NMovatm1:',t10,i5,/,
     $         'NMovatmn:',t10,i5)

          if (nmovatm.gt.MAXATM) goto 9990
       endif
c
       read (in,err=9990,end=9000) leexti,eextit(1:leexti)
c
       write (dest,600) leexti,(eextit(i:i),i=1,leexti)
 600   format ('Leexti:',t10,i5,/,t10,80a1)
c
       read (in,err=9990,end=9000) lparti,partit(1:lparti)
c
       write (dest,700) lparti,(partit(i:i),i=1,lparti)
 700   format ('Parti:',t10,i5,/,t10,80a1)
c
c     ----- end of header information -----
c
 1000  continue
c<mrsp>         read (in,err=9990,end=9000) curtim,itstep,system,avetem,
c<mrsp>     $     dstep,firstt,
c<mrsp>     1     finalt,e,eb,et,ep,ei,enb,eel,ehb,ec,eu,tint,tnb,
c<mrsp>     3     ea,eba,eta,epa,eia,enba,eela,ehba,eca,eua,tinta,tnba,
c<mrsp>     4     tote,totke,totea,tkea,iconmp,imstep,lvelwr,iconfs,
c<mrsp>     5     icstep
       read (in,err=9990,end=9000) curtim,itstep,system,avetem,
     $      dstep,firstt,
 1     finalt,e,eb,et,ep,ei,enb,eel,ehb,ec,eu,tint,tnb,
 3     ea,eba,eta,epa,eia,enba,eela,ehba,eca,eua,tinta,tnba,
 4     tote,totke,totea,tkea,iconmp,imstep,l1,iconfs,
 5     icstep
       lvelwr = (.NOT.(l1 .EQV. .FALSE.))
c
       write (dest,170,err=9995) curtim,itstep,system,avetem,
     $      dstep,firstt,
 1     finalt,e,eb,et,ep,ei,enb,eel,ehb,ec,eu,tint,tnb,
 3     ea,eba,eta,epa,eia,enba,eela,ehba,eca,eua,tinta,tnba,
 4     tote,totke,totea,tkea,iconmp,imstep,lvelwr,iconfs,
 5     icstep
 170   format ('Time:',t10,f14.4,i10,7(/,t10,4f14.4),
     $      /,t10,5f14.4,/,2i5,l5,2i5)
c     format ('Time:',t10,f10.4,i10,4(/,t10,7f10.4),
c     $        /,t10,5f10.4,2i5,l5,2i5)
c     
       if (version.ge.155) then
          read (in,err=9990,end=9000) pressura,vola,pvtota,pvkina,
     $         pvpota,radgyra
          write (dest,800,err=9995)  pressura,vola,pvtota,pvkina,
     $         pvpota,radgyra
 800      format ('Pressure:',t10,3f20.4,/,t10,3f20.4)
       end if
c     
       if (lcanon) then
          if (version.lt.300) then
             read (in,err=9990,end=9000) signose,zfrict,zprfrict
             write (dest,900,err=9995) signose,zfrict,zprfrict
 900         format ('Canonical:',t10,3f10.4)
          else
             if ( lnose ) then
                read (in,err=9990,end=9000) snose,snoseh,ssdot
     $               ,qcanon
                write (dest,901,err=9995) snose,snoseh,ssdot
     $               ,qcanon
 901            format ('Canonical(Nose) :',t10,4f10.4)
             else
                read (in,err=9990,end=9000) signose,zfrict,zprfrict
     $               ,qcanon
                write (dest,902,err=9995) signose,zfrict,zprfrict
     $               ,qcanon
 902            format ('Canonical(Hoover):',t10,4f10.4)
             endif
          endif
       end if
       if(version.ge.220)then
          if(period)then
             read (in,end=9000,err=9990) tcela,s2r,s2rdot
             write (dest,1005,err=9995) tcela,s2r,s2rdot
 1005        format ('DefCell:',t10,7f10.4,6f12.4)
          end if
       else
          if (defcel) then
             read (in,end=9000,err=9990) tcela,s2r
             write (dest,1050,err=9995) tcela,s2r
 1050        format ('DefCell:',t10,7f10.4)
          end if
       end if
       if (period) then
          if ( (version .eq. 210) .or.
     $         (version .ge. 300) ) then
             read (in,err=9990,end=9000) natmcel,strsa,extstrsa
             write (dest,1101,err=9995) natmcel,strsa,extstrsa
 1101        format ('Period:',t10,i5,7f10.4)
          else
             read (in,err=9990,end=9000) natmcel,strsa
             write (dest,1100,err=9995) natmcel,strsa
 1100        format ('Period:',t10,i5,6f10.4)
          endif
       end if
       if ( version .ge. 300 ) then
          if (period.and. lnpecan ) then
             read (in,err=9990,end=9000) sigdyn(1),sigdyn(2),qcanon
             write (dest,1140,err=9995) sigdyn(1),sigdyn(2),qcanon
 1140        format ('NPTCanonical:',t10,3f10.4)
          endif
          if(ltmpdamp) then
             read (in,err=9990,end=9000) GAMTMP
             write (dest,1160,err=9995) GAMTMP
 1160        format ('TempDamp:',t10,f10.4)
          endif
       endif
       if (prtthrm) then
          read (in,err=9990,end=9000) eabtota,eabvala,
$           eabelha,eabnba,eabmisa,dltfaba,expprta
      write (dest,1200,err=9995) eabtota,eabvala,eabelha,eabnba,
     $eabmisa,dltfaba,expprta
 1200 format ('PertTheory:',t10,7f10.4)
      end if
c     
      read (in,err=9990,end=9000) (x(i),i=1,nmovatm)
      read (in,err=9990,end=9000) (y(i),i=1,nmovatm)
      read (in,err=9990,end=9000) (z(i),i=1,nmovatm)
      write (dest,1300,err=9995) (i,x(i),y(i),z(i),i=1,nmovatm)
 1300 format ('Cordinat:',t10,2(i5,3f10.4),/,(t10,i5,3f10.4,
     $     i5,3f10.4))
*     
*        ----- velocities if needed -----
*
      if ( lvelwr )  then
         read (in,err=9990,end=9000) (velx(i),i=1,nmovatm)
         read (in,err=9990,end=9000) (vely(i),i=1,nmovatm)
         read (in,err=9990,end=9000) (velz(i),i=1,nmovatm)
         write (dest,1400,err=9995)
     $        (i,velx(i),vely(i),velz(i),i=1,nmovatm)
 1400    format ('Velocities:',t10,2(i5,3f14.4),/,(t10,i5,3f14.4,
     $        i5,3f14.4))
c     1400         format ('Velocities:',t10,2(i5,3f10.4),/,(t10,i5,3f10.4,
c     $        i5,3f10.4))
      endif
      go to 1000
c
c     normal completion.
c
 9000 continue
      write (*,9001)
 9001 format (/,1x,' Binary to ascii conversion completed...')
      go to 20000
c
c     I/O error.
c
c      write (*,10000)
10000 format (/,1x,' I/O error ....')
c      go to 20000

 9990 continue
      write(dest,9992)
 9992 format(/,' Input error ....')
      go to 20000
c
 9995 continue
      write (dest,9998)
 9998 format (/,' Output error ....')
      go to 20000
c     
20000 continue
      return
      end

      integer function leng1(STRING)
c       sh made integer 12:19:89 5pm
c       determines the length of STRING and adds 1.
c       this function is an alternative to LENGTH when you want to
c       avoid bounds errors
      character string*(*)
      leng1=len(string)
      do while(string(leng1:leng1).eq.' ')
         leng1=leng1-1
         if(leng1.eq.0)then
            leng1=1
            return
         endif
      enddo
c      leng1=leng1+1
      return
      end

*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "gomendian.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"

#include "stdafx.h"

#define MOD(a,b)  (((a+1)/b)*b - (a+1))
#define MXFILE             24
#define CERIUS2_TEXT_FIELD 80
#define CERIUS2_REMARKS    10

/* Read a frame from the trajectory file. */

/***************************************************************************/
int gomp_GetFrameCerius2(int alt, FILE *File_p, int iappend)  
    /* read frame number 'alt' from charmm trajectory   */
    /* mode of operation:
       ( = 0) Check file and print trajectory information
       ( > 0) get frame number trajectory number        */
    /* *File_p is file pointer to the trajectory file   */
    /* append, if = 0 no append , if = 1 append         */
/***************************************************************************/
{
/*  
    Program to read Cerius2 binary files

    leif laaksonen  made for gOpenMol 1995

    Modified 1998 for the Cerius2 version 3.x trajectories.

*/

/* 1 HDR , ICONTR  */
    static char hdr[4];      /* unknown header information */
    static int  icntrl[20];
    /* contains information about the datasets held in file

    (1) Version number * 100
    */

/* 2 TITLE         */
/*  static int ntitl;         (title(i,j), i=1,10),j=1,ntitl)
                              double precision title(10,10) */

/* 3 more TITLE    */
    static int neexti;


/* 4 */
    static int period;
    static int molxtl;
    static int lcanon;
    static int defcel;
    static int prtthrm;
    static int lnose;
    static int lnpecan;
    static int ltmpdamp;
    static int l1, l2, l3, l4, l5, l6, l7, l8;
    static int natom;
    static int nmovatm;
    static int leexti;
    static char eextit[CERIUS2_TEXT_FIELD];
    static int lparti;
    static char partit[CERIUS2_TEXT_FIELD];

    static int swap_bytes = 0;

/* version > 2010 */
    static double Dcurtim;

    static double pressura;
    static double vola;
    static double pvtota;
    static double pvkina;
    static double pvpota;
    static double radgyra;
    static double avpressura;
    static double avvola;
    static double avpvtota;
    static double avpvkina;
    static double avpvpota;
    static double avradgyra;

    static double signose;
    static double zfrict;
    static double   zprfrict;
    static double snose,snoseh,ssdot;
    static double qcanon;
    static double sigdyn[2];
    static double gamtmp;

    static double tcela;
    static double avtcela;
    static double mcell;
    static double exstr[6];
    static double Ds2r[6];
    static double Ds2rdot[6];

    static double natmcel;
    static double strsa[6];
    static double extstrsa; /* added after 220 but version called 210 */
    static double avstrsa[6];
    static double avextstrsa;

    static double express;
    static int  nflusd;
    static int  mvatmpfu[MXFILE];
    static int  natmpfu[MXFILE];
    static int  totmov;
    static int *mvatmofst;
    static int  lvelwr;
    static int  lforwr;
    static int  itstep;
    static int  iconmp;
    static int  imstep;
    static int  iconfs;
    static int  icstep;

    static char  decusd[80][MXFILE];

    static int version;


    static int icount,i;
    static int record;
    static int nstep;
    static long record_len1=0;        /* the "record" length in bytes of one
                                         record containing the x,y and z
                                         coordinates plus the information
                                         in between the coordinates (first 
                                         record)*/
    static long record_len2;          /* rest of the records if icntrl[8] > 0 */

    static long ret_fseek;
    static int idx1;
    static int nfreat; /* number of free atoms */

    static double  *tax; /* temp vectors in case of fixed atoms */
    static double  *tay;
    static double  *taz;
    static float *tax1; /* temp vectors for velocity/force */
    static float *tay1;
    static float *taz1;

    static char label[BUFF_LEN];
    static int  numset;  /* number of data sets */
    static int  StartRecord;
    static int  Wstr;

/* pointers to coordinates */
    static float *Xcoord;
    static float *Ycoord;
    static float *Zcoord;
    static const float *sumxyz;

    rewind(File_p);

    sumxyz   = gomp_GetTranslateArray();

/*  start reading  */
  
    if(!alt) {

        Wstr  = 0;
        natom = gomp_GetNumAtomsInMolecStruct(Wstr);

/* #1 */
        icount = fread(&record,sizeof(int), 1 ,File_p);

/* check for right type of trajectory
   record has to be = (4 chars + 20 * sizeof(int))!
*/
        swap_bytes = 0;
        if(record != (4 + 20 * sizeof(int))) {
            gomp_Reverse_int( & record );
            if(record != (4 + 20 * sizeof(int))) {
                gomp_PrintERROR("wrong internal structure of trajectory file");
                return(1);
            } else {
                gomp_PrintMessage("Enabling automatic byte_swapping...");
                swap_bytes = 1;
            }
        }
    
        icount = fread(hdr,sizeof(char), 4 ,File_p);
        icount = fread(icntrl,sizeof(int), 20 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int_array( icntrl, 20 );
        }
        icount = fread(&record,sizeof(int), 1 ,File_p);

        version = icntrl[0]; /* version number */

/* #2 */
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&neexti,sizeof(int), 1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &neexti );
        }

        if(neexti > CERIUS2_REMARKS) {
            gomp_PrintERROR("Too many Cerius Title Records (> 10)");
            return(1);
        }

        for(i = 0 ; i < neexti ; i++) {
            icount = fread(label, CERIUS2_TEXT_FIELD ,1,File_p);
            label[CERIUS2_TEXT_FIELD - 1] = '\0';
            gomp_PrintMessage(label);
        }
        icount = fread(&record,sizeof(int), 1 ,File_p);

/* this is an unclear record that I don't know */
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

/* #3 */
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&l1,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &l1 );
        }
        icount = fread(&l2,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &l2 );
        }
        icount = fread(&l3,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &l3 );
        }
        icount = fread(&l4,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &l4 );
        }
        icount = fread(&l5,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &l5 );
        }
        icount = fread(&l6,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &l6 );
        }
        icount = fread(&l7,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &l7 );
        }
        icount = fread(&l8,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &l8 );
        }
        icount = fread(&record,sizeof(int), 1 ,File_p);
        period  = l1;
        molxtl  = l2;
        lcanon  = l3;
        defcel  = l4;
        prtthrm = l5;
        lnose   = l6;
        lnpecan = l7;
        ltmpdamp= l8;
/* #4 ...*/
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&nflusd,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &nflusd );
        }
        if (nflusd > MXFILE) {
            gomp_PrintERROR("nflusd > MXFILE");
            return(1);
        }

        icount = fread(mvatmpfu,sizeof(int), nflusd , File_p);
        icount = fread(natmpfu,sizeof(int),  nflusd , File_p);

        for(i = 0 ; i < nflusd ; i++)
            icount = fread(decusd[i],sizeof(char), 8 , File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&totmov,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &totmov );
        }
        (void)gomp_SetFreeAtomListPointer( 0 , totmov);
        mvatmofst = gomp_GetModifiableFreeAtomListPointer(0);

        icount = fread(mvatmofst,sizeof(int), totmov ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int_array( mvatmofst, totmov );
        }
        nmovatm = totmov;
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&leexti,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &leexti );
        }
        icount = fread(eextit,sizeof(char), leexti , File_p);
        if(leexti) {
            eextit[79] = '\0';
            gomp_PrintMessage(eextit);
        }
        else 
            eextit[0]  = '\0';
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&lparti,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &lparti );
        }
        icount = fread(partit,sizeof(char), lparti , File_p);
        if(lparti) {
            partit[79] = '\0';
            gomp_PrintMessage(partit);
        } else 
            partit[0] = '\0';
/* end of record #4*/
        icount = fread(&record,sizeof(int), 1 ,File_p);

        StartRecord = ftell(File_p);

/* first time read everything */

        icount = fread(&record,sizeof(int), 1 ,File_p);

        lvelwr = 0;
        lforwr = 0;

        icount = fread(&Dcurtim,sizeof(double), 1 , File_p);
        icount = fread(&itstep,sizeof(int), 1 , File_p);

        icount = fseek(File_p, 57 * sizeof(double),SEEK_CUR);

        icount = fread(&iconmp,sizeof(int), 1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &iconmp );
        }
        icount = fread(&imstep,sizeof(int), 1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &imstep );
        }
        icount = fread(&lvelwr,sizeof(int), 1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &lvelwr );
        }
        icount = fread(&lforwr,sizeof(int), 1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &lforwr );
        }
        icount = fread(&iconfs,sizeof(int), 1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &iconfs );
        }
        icount = fread(&icstep,sizeof(int), 1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( &icstep );
        }
        icount = fread(&record,sizeof(int), 1 ,File_p);
/* #end */
/* next # */

        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&pressura,sizeof(double), 1 , File_p);
        icount = fread(&vola,sizeof(double), 1 , File_p);
        icount = fread(&pvtota,sizeof(double), 1 , File_p);
        icount = fread(&pvkina,sizeof(double), 1 , File_p);
        icount = fread(&pvpota,sizeof(double), 1 , File_p);
        icount = fread(&radgyra,sizeof(double), 1 , File_p);
        icount = fread(&avpressura,sizeof(double), 1 , File_p);
        icount = fread(&avvola,sizeof(double), 1 , File_p);
        icount = fread(&avpvtota,sizeof(double), 1 , File_p);
        icount = fread(&avpvkina,sizeof(double), 1 , File_p);
        icount = fread(&avpvpota,sizeof(double), 1 , File_p);
        icount = fread(&avradgyra,sizeof(double), 1 , File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
/* next # */
        if (lcanon) {
            if ( lnose ) {
                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(&snose,sizeof(double), 1 , File_p);
                icount = fread(&snoseh,sizeof(double), 1 , File_p);
                icount = fread(&ssdot,sizeof(double), 1 , File_p);
                icount = fread(&qcanon,sizeof(double), 1 , File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);
            }
            else {
                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(&signose,sizeof(double), 1 , File_p);
                icount = fread(&zfrict,sizeof(double), 1 , File_p);
                icount = fread(&zprfrict,sizeof(double), 1 , File_p);
                icount = fread(&qcanon,sizeof(double), 1 , File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);
            }
        }
/* next # */
        if(period) {
            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(&tcela,sizeof(double), 1 , File_p);
            icount = fread(&avtcela,sizeof(double), 1 , File_p);
            icount = fread(Ds2r,sizeof(double), 6 , File_p);
            icount = fread(Ds2rdot,sizeof(double), 6 , File_p);
            icount = fread(&mcell,sizeof(double), 1 ,File_p);
            icount = fread(&express,sizeof(double), 1 ,File_p);
            icount = fread(exstr,sizeof(double), 6 , File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);
        }
/* next # */
        if (period) {
            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(&natmcel,sizeof(double), 1 , File_p);
            icount = fread(strsa,sizeof(double), 6 , File_p);
            icount = fread(&extstrsa,sizeof(double), 1 , File_p);
            icount = fread(avstrsa,sizeof(double), 6 , File_p);
            icount = fread(&avextstrsa,sizeof(double), 1 , File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);
        }
/* next # */
        if (period && lnpecan ) {
            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(&sigdyn[0],sizeof(double), 1 , File_p);
            icount = fread(&sigdyn[1],sizeof(double), 1 , File_p);
            icount = fread(&qcanon,sizeof(double), 1 , File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);
        }
/* next # */
        if(ltmpdamp) {
            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(&gamtmp,sizeof(double), 1 , File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);
        }

/* determine the "record length"    */

        record_len1 = 6 * sizeof(int) +    /* the record contribution  */
            3 * totmov * sizeof(double); /* x-,y- and z-contribution */

        if(lvelwr) record_len1 += record_len1;

        if(lforwr) record_len1 += record_len1;

        record_len2 = ftell(File_p) - StartRecord;

/*                                  */
        ret_fseek    = fseek(File_p,0L,SEEK_END);
        if(ret_fseek) {
            sprintf(label,"?ERROR - can't read trajectory file : %s ",
                    gomp_GetTrajectoryFileName());
            gomp_PrintMessage(label);
            return(1);
        }

        numset  = (ftell(File_p) - StartRecord)/(record_len1 + record_len2);

        nstep   = numset;

        sprintf(label," Info for trajectory file   : %s   ",
                gomp_GetTrajectoryFileName());
        gomp_PrintMessage(label);
        sprintf(label," CERIUS2 version             : %d   ",version);
        gomp_PrintMessage(label);
        sprintf(label," Atoms found                : %d   ",natom);
        gomp_PrintMessage(label);
        sprintf(label," Free atoms                 : %d   ",totmov);
        gomp_PrintMessage(label);
        sprintf(label," Dynamics steps             : %d   ",nstep);
        gomp_PrintMessage(label);
        if(lvelwr) {
            gomp_PrintMessage(" Velocities available       : YES");
        } else {
            gomp_PrintMessage(" Velocities available       : NO");
        }
        if(lforwr) {
            gomp_PrintMessage(" Force components available : YES");
        } else {
            gomp_PrintMessage(" Force components available : NO");
        }      
/*
  sprintf(label," Time between data sets     : %d   ",icntrl[2]);
  gomp_PrintMessage(label);
  sprintf(label," Time of the first data set : %d   ",icntrl[1]);
  gomp_PrintMessage(label);
*/

/* update trajectory info ... */
        (void)gomp_SetNumberOfTrajectoryAtoms(natom);
        (void)gomp_SetNumberOfFrames(nstep);
        (void)gomp_SetTrajectoryTimeInfo(0 , 0);
        (void)gomp_SetNumberOfFreeAtoms(0 , totmov);
        (void)gomp_SetTrajectoryDisplayParams(1 , nstep , 1);
        (void)gomp_PutDisplayFrameNumber(1);
        return(0);
    }
/*     ----- end of header information ----- */

    ret_fseek    = fseek(File_p,(StartRecord + record_len2 + 
                                 (alt - 1) * (record_len1 + record_len2)) ,SEEK_CUR);

    if(ret_fseek) {
        sprintf(label,"?ERROR - can't read trajectory file : %s ",
                gomp_GetTrajectoryFileName());
        gomp_PrintMessage(label);
        return(1);
    }

/* loop starts here */
/*
  icount = fread(&record,sizeof(int), 1 ,File_p);
  icount = fread(x,sizeof(float), nmovatm , File_p);
  icount = fread(&record,sizeof(int), 1 ,File_p);

  icount = fread(&record,sizeof(int), 1 ,File_p);
  icount = fread(y,sizeof(float), nmovatm , File_p);
  icount = fread(&record,sizeof(int), 1 ,File_p);

  icount = fread(&record,sizeof(int), 1 ,File_p);
  icount = fread(z,sizeof(float), nmovatm , File_p);
  icount = fread(&record,sizeof(int), 1 ,File_p);
*/
/*        ----- velocities if needed ----- */


    if(!iappend) {  /* start NO append */

/*  get pointer to coordinate vectors */
        Wstr     = 0;
        Xcoord   = gomp_GetModifiableAtomXCoordPointer(Wstr);
        Ycoord  = gomp_GetModifiableAtomYCoordPointer(Wstr);
        Zcoord = gomp_GetModifiableAtomZCoordPointer(Wstr);

        if((natom - nmovatm)) {

            tax   = gomp_AllocateDoubleVector(nfreat);
            tay  = gomp_AllocateDoubleVector(nfreat);
            taz = gomp_AllocateDoubleVector(nfreat);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(tax,sizeof(double),nmovatm,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(tay,sizeof(double),nmovatm,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(taz,sizeof(double),nmovatm,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            if (swap_bytes ) {
                gomp_Reverse_double_array( tax, nmovatm );
                gomp_Reverse_double_array( tay, nmovatm );
                gomp_Reverse_double_array( taz, nmovatm );
            }

            for(i = 0 ; i < natom ; i++) {
                Xcoord[i]   += sumxyz[0];
                Ycoord[i]  += sumxyz[1];
                Zcoord[i] += sumxyz[2];
            }

/* put new values in */

            for(i = 0 ; i < nfreat; i++) {
                Xcoord[mvatmofst[i] - 1]   = (float)tax[i];
                Ycoord[mvatmofst[i] - 1]  = (float)tay[i];
                Zcoord[mvatmofst[i] - 1] = (float)taz[i];
            }


/* atom velocities */
/* retrieve velocities and forces if requested */
            if(lvelwr && gomp_GetVelocityRetrieveState()) {

                if(gomp_GetVelocitySpace(natom))
                    return(1);

                tax1 = gomp_GetModifiableVelocityXComponentPointer();
                tay1 = gomp_GetModifiableVelocityYComponentPointer();
                taz1 = gomp_GetModifiableVelocityZComponentPointer();

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tax,sizeof(double),nmovatm,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tay,sizeof(double),nmovatm,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(taz,sizeof(double),nmovatm,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                if (swap_bytes ) {
                    gomp_Reverse_double_array( tax, nmovatm );
                    gomp_Reverse_double_array( tay, nmovatm );
                    gomp_Reverse_double_array( taz, nmovatm );
                }

/* put new values in */

                for(i = 0 ; i < nfreat; i++) {
                    tax1[mvatmofst[i] - 1] = (float)tax[i];
                    tay1[mvatmofst[i] - 1] = (float)tay[i];
                    taz1[mvatmofst[i] - 1] = (float)taz[i];
                }

                if(gomp_CalculateVelocityMinMax())
                    gomp_PrintERROR("can't calculate velocity min/max values");
            }
/* atom forces     */
            if(lforwr && gomp_GetForceRetrieveState()) {

                if(gomp_GetForceSpace(natom))
                    return(1);

                tax1 = gomp_GetModifiableForceXComponentPointer();
                tay1 = gomp_GetModifiableForceYComponentPointer();
                taz1 = gomp_GetModifiableForceZComponentPointer();

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tax,sizeof(double),nmovatm,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tay,sizeof(double),nmovatm,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(taz,sizeof(double),nmovatm,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                if (swap_bytes ) {
                    gomp_Reverse_double_array( tax, nmovatm );
                    gomp_Reverse_double_array( tay, nmovatm );
                    gomp_Reverse_double_array( taz, nmovatm );
                }

/* put new values in */

                for(i = 0 ; i < nfreat; i++) {
                    tax1[mvatmofst[i] - 1]   = (float)tax[i];
                    tay1[mvatmofst[i] - 1]  = (float)tay[i];
                    taz1[mvatmofst[i] - 1] = (float)taz[i];
                }

                if(gomp_CalculateForceMinMax())
                    gomp_PrintERROR("can't calculate force min/max values");
            }
        }
        else {

            tax   = gomp_AllocateDoubleVector(natom);
            tay  = gomp_AllocateDoubleVector(natom);
            taz = gomp_AllocateDoubleVector(natom);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(tax,sizeof(double),natom,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(tay,sizeof(double),natom,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(taz,sizeof(double),natom,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            if (swap_bytes ) {
                gomp_Reverse_double_array( tax, natom );
                gomp_Reverse_double_array( tay, natom );
                gomp_Reverse_double_array( taz, natom );
            }

            for(i = 0 ; i < natom ; i++) {
                Xcoord[i]   = (float)tax[i];
                Ycoord[i]  = (float)tay[i];
                Zcoord[i] = (float)taz[i];
            }

/* atom velocities */
/* retrieve velocities and forces if requested */
            if(lvelwr && gomp_GetVelocityRetrieveState()) {

                if(gomp_GetVelocitySpace(natom))
                    return(1);

                tax1 = gomp_GetModifiableVelocityXComponentPointer();
                tay1 = gomp_GetModifiableVelocityYComponentPointer();
                taz1 = gomp_GetModifiableVelocityZComponentPointer();

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tax,sizeof(double),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tay,sizeof(double),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(taz,sizeof(double),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                if (swap_bytes ) {
                    gomp_Reverse_double_array( tax, natom );
                    gomp_Reverse_double_array( tay, natom );
                    gomp_Reverse_double_array( taz, natom );
                }

/* put new values in */

                for(i = 0 ; i < natom; i++) {
                    tax1[i]   = (float)tax[i];
                    tay1[i]  = (float)tay[i];
                    taz1[i] = (float)taz[i];
                }

                if(gomp_CalculateVelocityMinMax())
                    gomp_PrintERROR("can't calculate velocity min/max values");
            }
/* atom forces     */
            if(lforwr && gomp_GetForceRetrieveState()) {

                if(gomp_GetForceSpace(natom))
                    return(1);

                tax1 = gomp_GetModifiableForceXComponentPointer();
                tay1 = gomp_GetModifiableForceYComponentPointer();
                tay1 = gomp_GetModifiableForceZComponentPointer();

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tax,sizeof(double),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tay,sizeof(double),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(taz,sizeof(double),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                if (swap_bytes ) {
                    gomp_Reverse_double_array( tax, natom );
                    gomp_Reverse_double_array( tay, natom );
                    gomp_Reverse_double_array( taz, natom );
                }

/* put new values in */

                for(i = 0 ; i < natom; i++) {
                    tax1[i]   = (float)tax[i];
                    tay1[i]  = (float)tay[i];
                    taz1[i] = (float)taz[i];
                }

                if(gomp_CalculateForceMinMax())
                    gomp_PrintERROR("can't calculate force min/max values");
            }

        }
    }/*   end of NO append */
    else {   /*   start of append  */

        idx1     = 0;
/*  get pointer to coordinate vectors */
        sprintf(label,"Cerius2 frame (%d)",alt);
        Wstr   = gomp_CreateMolecStruct(label , natom , APPEND);
        if ( Wstr < 0 )
            return(1);
        Xcoord = gomp_GetModifiableAtomXCoordPointer(Wstr);
        Ycoord = gomp_GetModifiableAtomYCoordPointer(Wstr);
        Zcoord = gomp_GetModifiableAtomZCoordPointer(Wstr);

        if((natom - nmovatm)) {

/* temp vectors */
            tax   = gomp_AllocateDoubleVector(nfreat);
            tay  = gomp_AllocateDoubleVector(nfreat);
            taz = gomp_AllocateDoubleVector(nfreat);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(tax,sizeof(double),nmovatm,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(tay,sizeof(double),nmovatm,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(taz,sizeof(double),nmovatm,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            if (swap_bytes ) {
                gomp_Reverse_double_array( tax, nmovatm );
                gomp_Reverse_double_array( tay, nmovatm );
                gomp_Reverse_double_array( taz, nmovatm );
            }

            for(i = 0 ; i < natom ; i++) {
                Xcoord[i + idx1]   = Xcoord[i] + sumxyz[0];
                Ycoord[i + idx1]  = Ycoord[i] + sumxyz[1];
                Zcoord[i + idx1] = Zcoord[i] + sumxyz[2];
            }

/* put new values in */

            for(i = 0 ; i < nfreat; i++) {
                Xcoord[mvatmofst[i] - 1]   = (float)tax[i];
                Ycoord[mvatmofst[i] - 1]  = (float)tay[i];
                Zcoord[mvatmofst[i] - 1] = (float)taz[i];
            }
/* atom velocities */
/* retrieve velocities and forces if requested */
            if(lvelwr && gomp_GetVelocityRetrieveState()) {

                if(gomp_GetVelocitySpace(natom))
                    return(1);

                tax1 = gomp_GetModifiableVelocityXComponentPointer();
                tay1 = gomp_GetModifiableVelocityYComponentPointer();
                taz1 = gomp_GetModifiableVelocityZComponentPointer();

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tax,sizeof(double),nmovatm,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tay,sizeof(double),nmovatm,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(taz,sizeof(double),nmovatm,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                if (swap_bytes ) {
                    gomp_Reverse_double_array( tax, nmovatm );
                    gomp_Reverse_double_array( tay, nmovatm );
                    gomp_Reverse_double_array( taz, nmovatm );
                }

/* put new values in */

                for(i = 0 ; i < nfreat; i++) {
                    tax1[mvatmofst[i] - 1]   = (float)tax[i];
                    tay1[mvatmofst[i] - 1]  = (float)tay[i];
                    taz1[mvatmofst[i] - 1] = (float)taz[i];
                }

                if(gomp_CalculateVelocityMinMax())
                    gomp_PrintERROR("can't calculate velocity min/max values");
            }
/* atom forces     */
            if(lforwr && gomp_GetForceRetrieveState()) {

                if(gomp_GetForceSpace(natom))
                    return(1);
        
                Xcoord = gomp_GetModifiableForceXComponentPointer();
                Ycoord = gomp_GetModifiableForceYComponentPointer();
                Zcoord = gomp_GetModifiableForceZComponentPointer();

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tax,sizeof(double),nmovatm,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tay,sizeof(double),nmovatm,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(taz,sizeof(double),nmovatm,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                if (swap_bytes ) {
                    gomp_Reverse_double_array( tax, nmovatm );
                    gomp_Reverse_double_array( tay, nmovatm );
                    gomp_Reverse_double_array( taz, nmovatm );
                }

/* put new values in */

                for(i = 0 ; i < nfreat; i++) {
                    tax1[mvatmofst[i] - 1]   = (float)tax[i];
                    tay1[mvatmofst[i] - 1]  = (float)tay[i];
                    taz1[mvatmofst[i] - 1] = (float)taz[i];
                }

                if(gomp_CalculateForceMinMax())
                    gomp_PrintERROR("can't calculate force min/max values");
            }
        }
        else {

            tax   = gomp_AllocateDoubleVector(natom);
            tay  = gomp_AllocateDoubleVector(natom);
            taz = gomp_AllocateDoubleVector(natom);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(tax,sizeof(double),natom,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(tay,sizeof(double),natom,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(taz,sizeof(double),natom,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            if (swap_bytes ) {
                gomp_Reverse_double_array( tax, natom );
                gomp_Reverse_double_array( tay, natom );
                gomp_Reverse_double_array( taz, natom );
            }

            for(i = 0 ; i < natom ; i++) {
                Xcoord[i]   = (float)tax[i];
                Ycoord[i]  = (float)tay[i];
                Zcoord[i] = (float)taz[i];
            }
/* atom velocities */
/* retrieve velocities and forces if requested */
            if(lvelwr && gomp_GetVelocityRetrieveState()) {

                if(gomp_GetVelocitySpace(natom))
                    return(1);

                tax1 = gomp_GetModifiableVelocityXComponentPointer();
                tay1 = gomp_GetModifiableVelocityYComponentPointer();
                taz1 = gomp_GetModifiableVelocityZComponentPointer();

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tax,sizeof(double),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tay,sizeof(double),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(taz,sizeof(double),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                if (swap_bytes ) {
                    gomp_Reverse_double_array( tax, natom );
                    gomp_Reverse_double_array( tay, natom );
                    gomp_Reverse_double_array( taz, natom );
                }

/* put new values in */

                for(i = 0 ; i < natom; i++) {
                    tax1[i]   = (float)tax[i];
                    tay1[i]  = (float)tay[i];
                    taz1[i] = (float)taz[i];
                }

                if(gomp_CalculateVelocityMinMax())
                    gomp_PrintERROR("can't calculate velocity min/max values");
            }
/* atom forces     */
            if(lforwr && gomp_GetForceRetrieveState()) {

                if(gomp_GetForceSpace(natom))
                    return(1);

                tax1 = gomp_GetModifiableForceXComponentPointer();
                tay1 = gomp_GetModifiableForceYComponentPointer();
                taz1 = gomp_GetModifiableForceZComponentPointer();

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tax,sizeof(double),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tay,sizeof(double),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(taz,sizeof(double),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                if (swap_bytes ) {
                    gomp_Reverse_double_array( tax, natom );
                    gomp_Reverse_double_array( tay, natom );
                    gomp_Reverse_double_array( taz, natom );
                }

/* put new values in */

                for(i = 0 ; i < natom; i++) {
                    tax1[i]   = (float)tax[i];
                    tay1[i]  = (float)tay[i];
                    taz1[i] = (float)taz[i];
                }

                if(gomp_CalculateForceMinMax())
                    gomp_PrintERROR("can't calculate force min/max values");
            }
        }
    }
    /* end of append */


/* shift back using the translation info */

    for(i = 0 ; i < natom ; i++) {
        Xcoord[i]   -= sumxyz[0];
        Ycoord[i]  -= sumxyz[1];
        Zcoord[i] -= sumxyz[2];
    }
/* end of translation                    */

    free(tax);
    free(tay);
    free(taz);

    return(0);
/*                                  */

}




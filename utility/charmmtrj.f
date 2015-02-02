C
C  this is a small program to convert a CHARMM trajectory to an formatted
C  file and back to unformatted
C
C  Copyright:
C  Leif Laaksonen 1996
C
C IF YOU USE DIGITAL VISUAL FORTRAN THERE IS NO NEED TO ACTIVATE
C THE WIN32 PART IN THIS CODE
C
C  If compiled on Windows use the WIN32 definition
C  to select the right parts of the program
C
      parameter (maxa=50000)

      character*80 title(32)
      double precision xtlabc

      dimension xc(maxa),yc(maxa),zc(maxa)
      dimension icntrl(20)
      dimension xtlabc(6)
      dimension nfreelist(maxa)

      character*80 infile
      character*80 outfile
      character*4  ahdr

      data maxatom /50000/
C
C define action type
C
      write(6,*) 'Binary to formatted (1) or formatted to binary (2)'
	read(5,*) Action
	if(Action .gt. 2) then 
	   write(6,*) 'only 1 and 2 are alloved'
	   stop
	endif
	if(Action .lt. 1) then
	   write(6,*) 'only 1 and 2 are alloved'
	   stop
	endif

C
C input file name
C

      write(6,*) 'Give input file name'
      read(5,'(A)')   infile

C
C output file name
C      
      write(6,*) 'Give output file name'
      read(5,'(A)')   outfile

C
      if(Action .eq. 1) then
	else
      endif	 

c
c binary to formatted
c
      if(Action .eq. 1) then

      write(6,*) 'Processing now BINARY ==> FORMATTED'

      open(unit=20,status="OLD",file=infile,
     1 form="UNFORMATTED",ERR = 1000)
      open(unit=22,file=outfile,form="FORMATTED",
     1 ERR = 1000)
c
c reading input file 
c

      read(20) ahdr,icntrl
      read(20) ntitl,(title(i),i=1,ntitl)
      read(20) natom

      write(6,*) 'Statistics for file input file:'
      write(6,*) (title(i),i=1,ntitl)

c too many atoms (change value and recompile)
      if((natom .gt. maxatom)) then
         write(6,*) 'Too many atoms ',natom,
     1   ' (edit program file and recompile)'
          stop
      endif

      write(22,'(A)') ahdr
      write(22,*) (icntrl(i),i=1,20)
      write(22,*) ntitl
      write(22,'(A)') (title(i),i=1,ntitl)
      write(22,*) natom

      write(6,*) 'Number of atoms:  ',natom
      write(6,*) 'Number of frames: ',icntrl(1)

c
c crystal
c
      icrystal = 0

      if(icntrl(11) .gt. 0) then
	icrystal = 1
      endif

c
c free atoms
c

      nfree    = natom
      ifixed   = 0
      if(icntrl(9) .gt. 0)  then
	nfree    = natom - icntrl(9)
        read(20) (nfreelist(i),i=1,nfree)
	 write(22,*) (nfreelist(i),i=1,nfree)
	ifixed   = 1
      endif

      do i=1,icntrl(1)

       if(icrystal .gt. 0) then
          read(20) xtlabc
	    write(22,*) xtlabc
       endif

c      first frame is always saved complete
       if( i .eq . 1) then
           iatoms = natom
       else
           iatoms = nfree
       endif

         read(20) (xc(j),j=1,iatoms)
          write(22,*) (xc(j),j=1,iatoms)
         read(20) (yc(j),j=1,iatoms)
          write(22,*) (yc(j),j=1,iatoms)
         read(20) (zc(j),j=1,iatoms)
          write(22,*) (zc(j),j=1,iatoms)

	enddo

      else
c
c formatted to binary
c
c$if defined(WIN32)
c      write(6,*) "Making a special WIN32 binary file"
c      write(6,*) "=================================="
c$else
      write(6,*) "Making a UNIX/WIN32 unformatted file"
      write(6,*) "=============================="
c$endif

      write(6,*) 'Processing now FORMATTED  ==> BINARY '

      open(unit=20,status="OLD",file=infile,
     1 form="FORMATTED",ERR = 1000)
c$if defined(WIN32) 
c      open(unit=22,file=outfile,form="BINARY",
c     1 ERR = 1000)
c$else
      open(unit=22,file=outfile,form="UNFORMATTED",
     1 ERR = 1000)
c$endif

      read(20,'(A)') ahdr
      read(20,*) (icntrl(i),i=1,20)
      read(20,*) ntitl
      read(20,'(A)') (title(i),i=1,ntitl)
      read(20,*) natom

c$if defined(WIN32)
c      irecord = 2 * 4
c	write(22) irecord
c      write(22) ahdr,icntrl
c	write(22) irecord
c
c      irecord = 1
c	write(22) irecord
c      write(22) ntitl,(title(i),i=1,ntitl)
c	write(22) irecord
c
c      irecord = 4
c	write(22) irecord
c      write(22) natom
c	write(22) irecord
c$else
      write(22) ahdr,icntrl
      write(22) ntitl,(title(i),i=1,ntitl)
      write(22) natom
c$endif
      write(6,*) 'Statistics for file input file:'
      write(6,*) (title(i),i=1,ntitl)
      write(6,*) 'Number of atoms:  ',natom
      write(6,*) 'Number of frames: ',icntrl(1)

c
c crystal
c
      icrystal = 0
      if(icntrl(11) .gt. 0) then
	icrystal = 1
      endif

c
c free atoms
c
      nfree    = natom
      ifixed   = 0
      if(icntrl(9) .gt. 0)  then
	nfree    = natom - icntrl(9)
      write(6,*) 'Number of fixed atoms: ',icntrl(9)
      read(20,*) (nfreelist(i),i=1,nfree)
c$if defined(WIN32)
c      irecord = 4 * nfree
c	write(22) irecord
c	write(22) (nfreelist(i),i=1,nfree)
c	write(22) irecord
c$else
	write(22) (nfreelist(i),i=1,nfree)
c$endif
	ifixed   = 1
	endif

      do i=1,icntrl(1)

       if(icrystal .gt. 0) then

          read(20 , *) xtlabc

c$if defined(WIN32)
c          irecord = 8 * 6
c          write(22) irecord
c	    write(22) xtlabc
c          write(22) irecord
c$else
	    write(22) xtlabc
c$endif
       endif

c      first frame is always saved complete
       if( i .eq . 1) then
           iatoms = natom
       else
           iatoms = nfree
       endif

          read(20,*) (xc(j),j=1,iatoms)
          read(20,*) (yc(j),j=1,iatoms)
          read(20,*) (zc(j),j=1,iatoms)

c$if defined(WIN32)
c           irecord = 4 * iatoms
c           write(22) irecord
c           write(22) (xc(j),j=1,iatoms)
c           write(22) irecord
c
c           write(22) irecord
c           write(22) (yc(j),j=1,iatoms)
c           write(22) irecord
c
c           write(22) irecord
c           write(22) (zc(j),j=1,iatoms)
c           write(22) irecord
c$else
           write(22) (xc(j),j=1,iatoms)
           write(22) (yc(j),j=1,iatoms)
           write(22) (zc(j),j=1,iatoms)
c$endif

      enddo

      endif

      close(unit=20)
      close(unit=22)

      STOP 'Job is DONE!'

1000  write(6,*) 'ERROR in opening one of the files'
      STOP
      end


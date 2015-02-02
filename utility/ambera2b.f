C
C   This is a small program for transforming AMBER formatted
C   trajectories to binary for gOpenMol
C
C   Leif Laaksonen 1993
C
C IF YOU USE DIGITAL VISUAL FORTRAN THERE IS NO NEED TO ACTIVATE
C THE WIN32 PART IN THIS CODE
C
C To use this program on the WIN32 platform use the
C WIN32 definition to select right parts of the code
C
      parameter (MAXVEC = 30000)

      character*80 Infile,Outfile,label
      dimension   x(MAXVEC)

      iBox = 0
      write(6,*) 'Give input file name'
      read(5,'(A)') Infile
      write(6,*) 'Give output file name'
      read(5,'(A)') Outfile
      write(6,*) 'Give number of atoms'
      read(5,*) NumAtm
      write(6,*) 'Do you have box size included? [0 = no, 1 = yes]'
      read(5,*) iBox

      open(unit=21,status="OLD",file=Infile,
     1 form="FORMATTED",ERR = 1001)

c#if defined(WIN32)
c      open(unit=22,status="NEW",file=Outfile,
c     1 form="BINARY",ERR = 1002)
c#else
      open(unit=22,status="NEW",file=Outfile,
     1 form="UNFORMATTED",ERR = 1002)
c#endif

      iframes = 0
      iDummy  = 1

      read(21,'(a)') label
      write(6,*) 'Title: ',label
c#if defined(WIN32)
c      write(22) iDummy
c      write(22) label
c      write(22) Idummy
c#else
      write(22) label
c#endif

5     continue
C      read(21,'(10F8.3)',end=2000) (x(i),i = 1,3*Numatm)
c
c     This should also now work for non-standard files as long
c     as the data is writtenn as one long array
c
      read(21,*,end=2000) (x(i),i = 1,3*Numatm)

c#if defined(WIN32)
c
c      write(22) iDummy
c      write(22) (x(i),i = 1,3*NumAtm)
c      write(22) iDummy
c
c      if(iBox .gt. 0) then
c       read(21,*,end=2000) BoxX,BoxY,BoxZ
c        write(22) iDummy
c        write(22) BoxX,BoxY,BoxZ
c        write(22) iDummy
c      endif
c
c#else
      write(22) (x(i),i = 1,3*NumAtm)
      if(iBox .gt. 0) then
       read(21,*,end=2000) BoxX,BoxY,BoxZ
       write(22) BoxX,BoxY,BoxZ
      endif
c#endif
      iframes = iframes + 1
      goto 5

1001  write(6,*) '?ERROR - in formatted file'
      stop
1002  write(6,*) '?ERROR - in unformatted file'
      stop
2000  write(6,*) 'Found ',iframes,' frames'

      end

C
C    This will not work for trajectories produced with CHARMM22
C    with constant pressure. I will try to fix it as soon as possible!!!!
C
c
c    small program to merge two charmm trajectories into one
c
c    leif laaksonen, centre for scientific computing, 1992
c
c    This program asks for two input files and one output file.
c    
c    The new output file will contain the header from input file one
c    and the data from input file one first and input file two
c    appended to it.
c
c    input1 + input2 ==> {input1input2}
c

      parameter (maxa=10000)

      character*80 title1(32)
      character*80 title2(32)

      dimension xc(maxa),yc(maxa),zc(maxa)
      dimension icntrl(20)
      dimension icntrl1(20)
      dimension icntrl2(20)

      character*80 infile1
      character*80 infile2

      character*80 outfile

      data maxatom / 10000/

      write(6,*) 'Give input file 1'
      read(5,'(A)')   infile1
      write(6,*) 'Give input file 2'
      read(5,'(A)')   infile2

      write(6,*) 'Give output file'
      read(5,'(A)')   outfile

      open(unit=20,status="OLD",file=infile1,
     1 form="UNFORMATTED",ERR = 1000)
      open(unit=21,status="OLD",file=infile2,
     1 form="UNFORMATTED",ERR = 1000)
      open(unit=22,file=outfile,form="UNFORMATTED",
     1 ERR = 1000)

c start with file 1

      read(20) ahdr,icntrl1
      read(20) ntitl,(title1(i),i=1,ntitl)

      read(20) natom1
      write(6,*) 'Statistics for file 1:'
      write(6,*) (title1(i),i=1,ntitl)
      write(6,*) 'NATOM: ',natom1

c start with file 2

      read(21) ahdr,icntrl2
      read(21) ntitl,(title2(i),i=1,ntitl)

      read(21) natom2
      write(6,*) 'Statistics for file 2:'
      write(6,*) (title2(i),i=1,ntitl)
      write(6,*) 'NATOM: ',natom2

      if(natom1 .ne. natom2) then
        write(6,*) 'Number of atoms does not match'
         stop
      endif

      if(icntrl1(3) .ne. icntrl2(3)) then
        write(6,*) 'Time step in data sets does not match'
         stop
      endif

      if((natom1 .gt. maxatom)) then
         write(6,*) 'Too many atoms ',natom1,
     1   ' (edit program file and recompile)'
          stop
      endif

c write to output
      do i=1,10
       icntrl(i) = icntrl1(i)
      enddo

      icntrl(1) = icntrl1(1) + icntrl2(1)
      write(22) ahdr,icntrl
      write(22) ntitl,(title1(i),i=1,ntitl)
      write(22) natom1

      write(6,*) 'Writing set 1...'

      do i=1,icntrl1(1)
       read(20) (xc(j),j=1,natom1)
        write(22) (xc(j),j=1,natom1)
       read(20) (yc(j),j=1,natom1)
        write(22) (yc(j),j=1,natom1)
       read(20) (zc(j),j=1,natom1)
        write(22) (zc(j),j=1,natom1)
      enddo

      write(6,*) 'Writing set 2...'

      do i=1,icntrl2(1)
       read(21) (xc(j),j=1,natom2)
        write(22) (xc(j),j=1,natom2)
       read(21) (yc(j),j=1,natom2)
        write(22) (yc(j),j=1,natom2)
       read(21) (zc(j),j=1,natom2)
        write(22) (zc(j),j=1,natom2)
      enddo

      close(unit=20)
      close(unit=21)
      close(unit=22)

      STOP

1000  write(6,*) 'ERROR in opening of file'
      STOP
      end


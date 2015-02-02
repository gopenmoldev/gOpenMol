C
C   NO TABS IN THE DATA FILE!!!!!!!
C
C   This is a small program for transforming XMOL formatted
C   trajectories into AMBER binary ones
C
C   Leif Laaksonen 1996
C
      parameter (MAXVEC = 30000)

      character*80 Infile,Outfile,label
      character*80 Aname(MAXVEC)
      character*80 Temp
      dimension   x(MAXVEC)
      dimension   BoxSize(3)

      write(6,*) 'Give input file name'
      read(5,'(A)') Infile
      write(6,*) 'Give output file name'
      read(5,'(A)') Outfile

      open(unit=21,status="OLD",file=Infile,
     1 form="FORMATTED",ERR = 1001)

      open(unit=22,status="NEW",file=Outfile,
     1 form="UNFORMATTED",ERR = 1002)

      iframes = 0

      BoxSize(1) = 1.0
      BoxSize(2) = 1.0
      BoxSize(3) = 1.0

5     continue

      read(21,'(i6)',end = 3000) NumAtm
      read(21,'(a)',end = 2000) label
      if(iframes .eq. 0) then
         write(6,*) 'Title: ',label
         write(22) label
      endif

      do i = 1,NumAtm
      jj = 3 * (i -1)

      read(21,'(a)',end = 2000) Temp
C
C I'm not that good in fortran. Most likely
C this can be done in a more clever way ...
C
      icharhit = 0
      do ii = 1,len(Temp)
      if((icharhit .eq. 1) .and. (Temp(ii:ii) .eq. " ")) then
         goto 10
      endif
      if(Temp(ii:ii) .ne. " ") then
         icharhit = 1
      endif
      enddo
 10   continue

      read(Temp(ii:len(Temp)),*) x(jj+1),x(jj+2),x(jj+3)

      enddo

      write(22) (x(k),k = 1,3*NumAtm)
      write(22) BoxSize
      iframes = iframes + 1
      goto 5

1001  write(6,*) '?ERROR - in formatted file'
      stop
1002  write(6,*) '?ERROR - in unformatted file'
      stop
2000  write(6,*) 'corrupted trajectory file, uncompleted frame'
      stop
3000  write(6,*) 'Found ',iframes,' frames'

      end

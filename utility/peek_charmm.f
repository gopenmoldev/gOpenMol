C
C     Small fortran program to scan the trajectory
C
C     Leif laaksonen 1993, 1994
C
c     The icntrl array contains important information:
c
c     if icntrl(9)  > 0 there are fixed atoms
c        icntrl(11) > 0 there is crystal information
c
      integer MAXATOMS, MAXSTEPS
      parameter (MAXATOMS = 50000)
      parameter (MAXSTEPS = 15000)
      dimension x(MAXATOMS)
      dimension y(MAXATOMS)
      dimension z(MAXATOMS)
      dimension ifree(MAXATOMS)
      integer near(MAXATOMS),count,n
      integer atomcent
      double precision xtlabc(6)
      dimension icntrl(20)
      character*256 InFile
      character hdr*4
      character*4 type(MAXATOMS)
      character*80 title(32)
c
c Input binary .crd trajectory filename and output filename
c
      write (6,1001)
      read (5,1002) InFile
      open(unit = 21 , status = "OLD" , file = InFile ,
     -     form="UNFORMATTED" , ERR = 100)
c     Header
      read (21) hdr,icntrl
      write (6,1000) hdr
      write (6,*) 'Number of data sets: ',icntrl(1)
c     Tile lines
      read (21) ntitl,(title(i),i=1,ntitl)
      write (6,1003) (title(i),i=1,ntitl)
      read(21) natom
      write(6,*) 'Number of atoms: ',natom
      write (6,1111) icntrl(1),icntrl(3),icntrl(2),icntrl(5),
     -  icntrl(8),icntrl(9) 
      nfreat = natom - icntrl(9)
      print *,'Number of free atoms: ',nfreat
C free atom array
      if (icntrl(9) .gt. 0) read (21) (ifree(i),i=1,nfreat)
      if (icntrl(11) .eq. 0) print *,'No crystal information'
      if (icntrl(11) .gt. 0) then
        print *,'Crystal information is present'
        write (6,1115)
        read *,ivers
      end if
      if (natom .gt. MAXATOMS) then
        print *,'Redimension with maxatom=',natom,
     -    ' if you want to scan through the file'
        stop
      end if
      nread=0
      do while (.true.)
c       Crystal information
        if (nread .eq. 0) then
          if (icntrl(11) .gt. 0) read (21,end=999) (xtlabc(j),j=1,6)
          read(21,end=999) (x(j),j=1,natom)
          read(21,end=998) (y(j),j=1,natom)
          read(21,end=998) (z(j),j=1,natom)
        else
          if (icntrl(11) .gt. 0 .and. ivers .ge. 27)
     -      read (21,end=999) (xtlabc(j),j=1,6)
          if (icntrl(9) .gt. 0) then
            read(21,end=999) (x(ifree(j)),j=1,nfreat)
            read(21,end=998) (y(ifree(j)),j=1,nfreat)
            read(21,end=998) (z(ifree(j)),j=1,nfreat)
          else
            read(21,end=999) (x(j),j=1,natom)
            read(21,end=998) (y(j),j=1,natom)
            read(21,end=998) (z(j),j=1,natom)
          end if
        end if
        nread=nread+1
      end do
      stop
999   if (icntrl(1) .eq. 0 .or. icntrl(1) .eq. nread) then
        write (6,1112) nread,natom
      else
        write (6,1114) nread
      end if 
      stop
998   write (6,1113) nread+1
      stop
100   print *, 'Cant open input file properly'
      stop
1000  format(' HDR is: ',a)
1001  format(' Name of the trajectory file=',$)
1002  format(a)
1003  format(' Title lines:',/,(1x,a))
1111  format(' Number of coordinate sets in file=',i8,/,
     -  ' Frequency for saving coordinates=',i8,/,
     -  ' Number of steps for creation run=',i8,/,
     -  ' Frequency for saving velocities=',i8,/,
     -  ' Number of degrees of freedom=',i8,/,
     -  ' Number of fixed atoms=',i8)
1112  format(' Succesfully read ',i8,' records of ',i8,' atoms')
1113  format(' Trajectory file is truncated at the ',i6,'-th snapshot')
1114  format(' Trajectory file is short. Last snapshot is the',i6,'-th')
1115  format(' Crystal information is present',/,
     -  ' Which Charmm version (19? ... 27?) created this file ? ',$)
      end
 

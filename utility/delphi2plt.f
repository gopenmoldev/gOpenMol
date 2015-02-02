C
C read an unformatted potential file from Delphi
C in the phi format and convert it into a plt
C file recognized by gOpenMol
C
C Leif Laaksonen/CSC
C 2001-04-19
C
C               unformatted (binary file)
      character*20 uplbl
      character*10 nxtlbl
      character*60 toplbl
      character*16 botlbl
      real*4 scale,oldmid(3)

      real*4 phimap(200 , 200 , 200)

      character*132 delphifile,pltfile

      do i1 = 1,200
         do i2 = 1,200
            do i3 = 1,200
               phimap(i3,i2,i1) = 0.0
            end do
         end do
      end do

      write(6,*) 'Delphi .phi grid file name:'
      read (5,991) delphifile
      write(6,*) 'Desired .plt output file name:'
      read (5,991) pltfile
 991  format (a132)
      write(6,*) 'Gridpoints in x, y and z directions:'
      read(5,*) intx,inty,intz

      if(intx .gt. 200) then
        write(6,*) 'Array (x dim) defined for < 200'
        call exit
      endif
      if(inty .gt. 200) then
        write(6,*) 'Array (y dim) defined for < 200'
        call exit
      endif
      if(intz .gt. 200) then
        write(6,*) 'Array (z dim) defined for < 200'
        call exit
      endif

      open (unit=14,file=delphifile,form='unformatted')

      read(14) uplbl
      read(14) nxtlbl,toplbl
      write(6,*) "Title: ",toplbl

      read(14) (((phimap(i,j,k),i=1,intx),j=1,inty),k=1,intz)

      read(14) botlbl
      read(14) scale, oldmidx, oldmidy, oldmidz

      ixmid = (intx + 1)/2
      iymid = (inty + 1)/2
      izmid = (intz + 1)/2

      xstart  = (1    - ixmid)/scale + oldmidx
      xend    = (intx - ixmid)/scale + oldmidx
      ystart  = (1    - iymid)/scale + oldmidy
      yend    = (intx - iymid)/scale + oldmidy
      zstart  = (1    - izmid)/scale + oldmidz
      zend    = (intx - izmid)/scale + oldmidz
C
      write(6,*) 'NptsX, NptsY, NptsZ:', intx,inty,intz
      write(6,*) 'min/max x, y and z:',
     1xstart,xend,ystart,yend,zstart,zend
C
      open (unit=11,file=pltfile,form='unformatted')
C this defines that it is 3D information and 204 = Delphi type file
      write (11) 3,204
      write (11) intz,inty,intx
      write (11) zstart,zend,ystart,yend,xstart,xend

      do  200  k = 1 , intz

        write(11) ((phimap(i,j,k), i=1,intx),j=1,inty)

200   continue

      close(14)
      close(11)



      end

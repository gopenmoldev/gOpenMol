C
C read an unformatted potential file from Delphi
C in the Insight format and convert it into a plt
C file recognized by gOpenMol
C
C Leif Laaksonen/CSC
C 2001-04-19
C
C This is odd because according to the manual this shoul be
C character*132 !?!?!?!?!?!?
C
      character*60  toplbl    
C!ascii header

      integer*4  ivary         
C!0  =>  x  index  varys  most rapidly

      integer*4 nbyte          
C!=4, # of bytes in data

      integer*4 inddat         
C!=0, floating point data

      real*4 xang,yang,zang    
C!=90,90,90 unit cell angles

      integer*4 intx,inty,intz 
C!=igrid-1, # of intervals/grid side
                  
      real*4 extent            
C!maximum extent of grid

      real*4 xstart,xend       
C!beginning, end of grid sides

      real*4 ystart,yend       
C!in fractional

      real*4 zstart,zend       
C!units of extent

      real*4 phimap(200 , 200 , 200)

      character*132 delphifile,pltfile

      do i1 = 1,200
         do i2 = 1,200
            do i3 = 1,200
               phimap(i3,i2,i1) = 0.0
            end do
         end do
      end do

       write(6,*) 'Delphi Insight .phi grid file name:'
      read (5,991) delphifile
      write(6,*) 'Desired .plt output file name:'
      read (5,991) pltfile
 991  format (a132)

      open (unit=14,file=delphifile,form='unformatted')

      read(14) toplbl
      write(6,*) 'Title: ',toplbl
      read(14) ivary, nbyte, intdat, extent, extent, extent,
     &xang, yang, zang, xstart, xend, ystart,  yend,  zstart,
     &zend, intx, inty, intz

      intx = intx + 1
      inty = inty + 1
      intz = intz + 1

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

      do k = 1,intz
       do j = 1,inty
              read(14) (phimap(i,j,k),i=1,intx)
       end do
      end do

      xend  = extent * xend
      xstart= extent * xstart
      yend  = extent * yend
      ystart= extent * ystart
      zend  = extent * zend
      zstart= extent * zstart
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

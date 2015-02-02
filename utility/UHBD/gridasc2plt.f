C
C Unknown origin
C
C Modified by Leif Laaksonen/CSC 1997 to do the conversion into
C the plot file format supported by gOpenMmol
C Modified by Eero Häkkinen/CSC 2003 to compile by g77.
C 

      program uhbd_to_plt
convert a 200**3 or less ascii uhbd .grd to plt unformatted form

      real phi(200,200,200),oldmid(3),scale,h,ox,oy,oz
      character*10 nxtlbl
      character*80 toplbl,uhbdfile,pltfile
      character*72 title

      integer im,jm,km,idum1,idum2,idum3,idum4,grdflg,kdum1,kdum2
      real dum2,dum3,dum4,dum5,dum6,dum7,dum8,phimax,phimin

      do i1 = 1,200
         do i2 = 1,200
            do i3 = 1,200
               phi(i3,i2,i1) = 0.
            end do
         end do
      end do

      write(6,*)
     .'I assume ascii uhbd output; must be 200x200x200 or less'
 100  continue
      write(6,*) 'UHBD .phi grid file name:'
      read (5,991) uhbdfile
      write(6,*) 'Desired .plt output file name:'
      read (5,991) pltfile
 991  format (a80)

      open (unit=10,file=uhbdfile,form='formatted')

 77   format(a72/2e12.6,5i7/3i7,4e12.6/4e12.6/2e12.6,2i7)
 78   format(3i7)
 79   format( (6(1x,e12.6)) )

      read (10,77) title, scale, dum2, grdflg, idum2,
     &                    km, idum5, km, im, jm, km,
     &                    h, ox, oy, oz,
     &                    dum3, dum4, dum5, dum6, dum7,
     &                    dum8, idum3, idum4

      write(6,*) title
      write(6,*) 'im,jm,km:', im,jm,km
      write(6,*) 'h,ox,oy,oz:', h,ox,oy,oz

      do  20  k = 1 , km
         read (10,78) k1, im, jm

         read (10,79) ( ( phi (i,j,k), i=1,im ), j=1,jm )
 20   continue

      phimin = 99999.
      phimax = -99999.

      do i1 = 1,km
         do i2 = 1,jm
            do i3 = 1,im
               phimax = max(phi(i3,i2,i1), phimax)
               phimin = min(phi(i3,i2,i1), phimin)
            endd o
         end do
      end do
c
C
      ox1 = ox + h
      oy1 = oy + h
      oz1 = oz + h
C
      write(6,*) 'min/max phi:', phimin,phimax
      write(6,*) 'NptsX, NptsY, NptsZ:', im,jm,km
      write(6,*) 'min/max x, y and z:',
     1ox1,ox+h*im,oy1,oy+h*jm,oz1,oz+h*km
C
      open (unit=11,file=pltfile,form='unformatted')
C this defines that it is 3D information and 120 = UHBD type file
      write (11) 3,120
      write (11) km,jm,im
      write (11) oz1,oz+h*km,oy1,oy+h*jm,ox1,ox+h*im

      do  200  k = 1 , km

        write(11) ((phi(i,j,k), i=1,im),j=1,jm)

200   continue

      close(10)
      close(11)
      end


C
C
C Modified by Leif Laaksonen/CSC 1999 to do the conversion into
C the plot file format supported by gOpenMmol
C 
C Information atke from: 
C www.scripps.edu/pub/olson-web/doc/autodock/ad3/AD305UG.21.html#pgfId=34041
C
C AUTODOCK Version 3
C III 4. Grid Map File
C Extension: .map
C The first six lines of each grid map hold header information which describe 
C the spatial features of the maps and the files used or created. These headers
C are checked by AutoDock to ensure that they are appropriate for the requested
C docking. The remainder of the file contains grid point energies, written as 
C floating point numbers, one per line. They are ordered according to the 
C nested loops z( y( x ) ). A sample header from a grid map is shown below:
C______________________________________________________________________________ 
C
C GRID_PARAMETER_FILE vac1.nbc.gpf
C GRID_DATA_FILE 4phv.nbc_maps.fld
C MACROMOLECULE 4phv.new.pdbq
C SPACING 0.375
C NELEMENTS 50 50 80
C CENTER -0.026 4.353 -0.038
C 125.095596
C 123.634560
C 116.724602
C 108.233879
C :
C______________________________________________________________________________


      program autodock_to_plt

      real phi(200,200,200),h,ox,oy,oz,ox1,oy1,oz1,spacing
      character*80 toplbl,autodockfile,pltfile
      character*80 title

      integer im,jm,km
      real temp1,phimax,phimin

      write(6, *) 
     x'Convert an AutoDock map file to a plt file for gOpenMol'
      write(6, *) 'Leif Laaksonen/2004'
      write(6, *) 'Give AutoDock map file name:'
      read (5,991) autodockfile
      write(6, *) 'Desired .plt output file name:'
      read (5,991) pltfile
 991  format (a80)

      open (unit=10,file=autodockfile,form='formatted')
C GRID_PARAMETER_FILE
      read(10,991) title
c      write(6,991) title
C GRID_DATA_FILE
      read(10,991) title
c      write(6,991) title
C MACROMOLECULE
      read(10,991) title
c      write(6,991) title
C SPACING
      read(10,*) title,spacing
      Write(6, *) 'Spacing: ',spacing
C NELEMENTS
      read(10,*) title,im,jm,km
      write(6, *) 'Gridpoints (x,y,z): ',(im+1),(jm+1),(km+1)

      if((im+1 .gt. 200) .or. (jm+1 .gt. 200) .or. (km+1 .gt. 200)) then
        write(6,*) 'dimension is > 200 x 200 x 200'
        stop
      endif

C CENTER
      read(10,*) title,ox,oy,oz
      write(6, *) 'Center (x,y,z): ',ox,oy,oz

      ox1 = ox - im * (spacing / 2.0)
      oy1 = oy - jm * (spacing / 2.0)
      oz1 = oz - km * (spacing / 2.0)

c      write(6,*) spacing,im,jm,km,ox,oy,oz

      im = im + 1
	jm = jm + 1
	km = km + 1

      do  20  k = 1 , km
         read (10,*) ( ( phi (i,j,k), i=1,im ), j=1,jm )
 20   continue

      phimin = 99999
      phimax = -99999

      do i1 = 1,km
         do i2 = 1,jm
            do i3 = 1,im
               phimax = max(phi(i3,i2,i1), phimax)
               phimin = min(phi(i3,i2,i1), phimin)
            end do
         end do
      end do
c
      h = spacing
      write(6, *) 'min/max phi:', phimin,phimax
      write(6, *) 'NptsX, NptsY, NptsZ:',im,jm,km
      write(6, *) 'min/max x, y and z:',
     z ox1,ox1 + h * im,
     z oy1,oy1 + h * jm,
     z oz1,oz1 + h * km
c

      open (unit=11,file=pltfile,form='unformatted')
C this defines that it is 3D information and 203 = AUTODOCK type file
      write (11) 3,203
      write (11) km,jm,im
      write (11) 
     z oz1,oz1 + h * km,
     z oy1,oy1 + h * jm,
     z ox1,ox1 + h * im

      do  200  k = 1 , km

        write(11) ((phi(i,j,k), i=1,im),j=1,jm)

200   continue

      write(6,*) 'Done!'

      close(10)
      close(11)
      end


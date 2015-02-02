c 
c Program to convert the Grid kont files to gOpenMol plt files
c
c Leif Laaksonen 1998
c

      program kont2plt
c max number of grids
      parameter (maxgrd=200)
      character*120 filein,fileout
      dimension e(maxgrd,maxgrd,maxgrd)
      integer       nx
      integer       ny
      integer       nz
      character*60  head
      character*72  headlong
      real*4        x_len, y_len, z_len
      real*4        alpha, beta,  gamma
      real*4        x_start, y_start, z_start
      real*4        x_end  , y_end  , z_end
      mxgrd=maxgrd
      infile   = 20
      ioutfile = 21

      write (6,100)
 100   format (
     +' Type the name of your unformatted gridkont file: ')
 160   read (5,170) filein
 170   format (a)
      write(6,*) 
     >'Convert to formatted (1) or unformatted (2) plt file: '
      read(5,*) Ialt

      if(Ialt .eq. 0) Ialt = 1

      if((Ialt .lt. 1) .or. (Ialt .gt. 2)) then
       write(6,*) 'Wrong option: ',Ialt,' it has to be 1 or 2'
       stop
      endif 
      write (6,190)
 190       format (/,
     +' Type the name of the output plot file: ')
      read (5,170) fileout
c
c start processing the file
c
      open (infile,file=filein,form='unformatted')

C----- Get the title and test in which format the file is in?
      read ( infile , err = 900) head

      goto 910
900   continue
      rewind(infile)
      read ( infile , err = 990) 
     +headlong,cler,emax,neta,nixnix,npla,nz1,nz2,nx,ny,nz,
     +ra,rx,ry,rz,vdwrj,effnj,alphj,qj,eminj,rminj,jd,ja

      write(6,'(A)') 'Title: ',headlong
      xstart = rx + ra
      ystart = ry + ra
      zstart = rz + ra

      xend   = rx+ra*real(nx)
      yend   = ry+ra*real(ny)
      zend   = rz+ra*real(nz)

      do 400 ii=1,nz
        read(infile) idz1,idz2,idz3
        read (infile) ((e(kk,jj,ii), kk=1,nx),jj=1,ny)
 400  continue

      goto 920
910   continue

      write(6,'(A)') 'Title: ',head
C----- Get the second line
      read ( infile ) vary, data_size, data_type,
     > x_len, y_len, z_len, alpha, beta,  gamma,
     > x_start, x_end,
     > y_start, y_end,
     > z_start, z_end,
     > nx, ny, nz

      xend  =x_len*x_end
      xstart=x_len*x_start
      yend  =y_len*y_end
      ystart=y_len*y_start
      zend  =z_len*z_end
      zstart=z_len*z_start

C----- Get the number of points in each direction
      nx = nx + 1
      ny = ny + 1
      nz = nz + 1

      do 300 ii=1,nz
       do 310 jj=1,ny
        read (infile) (e(kk,jj,ii), kk=1,nx)
 310   continue
 300  continue

920   continue

C
      write(6,*) 'NptsX, NptsY, NptsZ:', nx,ny,nz
      write(6,*) 'Min/max x, y and x:' , 
     1xstart,xend,ystart,yend,zstart,zend
C
      if(Ialt .eq. 1) then
      open (ioutfile,file=fileout,form='formatted')
C FORMATTED
C this defines that it is 3D information and 205 = GRID type file
      write (ioutfile,*) 3,205
      write (ioutfile,*) nz,ny,nx
      write (ioutfile,*) zstart,zend,ystart,yend,xstart,xend

      do  500  k = 1 , nz

        write(ioutfile,*) ((e(i,j,k), i=1,nx),j=1,ny)

 500       continue
      else 
      open (ioutfile,file=fileout,form='unformatted')
C UNFORMATTED
C this defines that it is 3D information and 190 = GRID type file
      write (ioutfile) 3,190
      write (ioutfile) nz,ny,nx
      write (ioutfile) zstart,zend,ystart,yend,xstart,xend

      do  510  k = 1 , nz

        write(ioutfile) ((e(i,j,k), i=1,nx),j=1,ny)

 510       continue

      endif

C Done!

      close (infile)
      close (ioutfile)
      stop
990   continue
      write(6,*) 'ERROR! Can not read header information'
      stop
      end






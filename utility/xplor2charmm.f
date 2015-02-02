
C Date: Tue, 2 Aug 1994 18:36:40 -0700 (PDT)
C From: INGRID@uoxray.uoregon.edu
C
C The XPLOR to CHARMM-converter is a really primitive program,
C but in case it might still be useful (it is at least easy to 
C modify..), I'll include it here. I need the intermediate conversion 
C to an ASCII file since I'm running XPLOR on a non-SGI machine with
C a different byte order. I'm not sure if the program will treat
C fixed atoms correctly, I didn't check that. Another problem
C is the timestep, it is *not* the correct one I've used in the
C XPLOR commandfile, it seems to be only vaguely related. This 
C didn't matter for me, so I didn't look into it.
C
      program xplor2charm
c
c converts xplor-binary-trajectories to charmm-trajectories 
c for input into dynal
c
      character*80 infile,outfile,title(32)
      character*4 hdr
      integer io,istart,nsavc,diff,iqfirst,natom,nfreat,ntitle
      integer freeat(10000)
      double precision x(10000),y(10000),z(10000),delta,dum(10,10)
      real temp(10000),rdelta,rx(10000),ry(10000),rz(10000)

      write(*,*) 'input file ? '
      read(5,10) infile
      write(*,*) 'output file ? '
      read(5,10) outfile
10    format(a80)
      write(*,*) 'xplor to ascii (1) or ascii to charmm ? (2)'
      read(*,*) ians

      if(ians.eq.1) then

      ntra=0
      open(1,file=infile,form='unformatted',status='old')
      open(2,file=outfile,recl=32700,status='unknown')
  
      iqfirst=1
c read header of xplor trajectory
      read(1,end=500) hdr,io,istart,nsavc,io,io,io,io,io,diff,delta,
     &        io,io,io,io,io,io,io,io,io
c      read(1,end=500) idum,idum
      read(1) ntitle,(title(i)(1:80),i=1,ntitle)
      read(1,end=500) natom
      nfreat=natom-diff 
      if(nfreat.ne.natom) then 
          read(1,end=500) (freeat(i),I=1,nfreat) 
      else
       do 300 i=1,natom
300       freeat(i)=i
      end if
      write(*,*) natom,' atoms found, ',nfreat,' free. '
      delta=delta/2.0
      write(*,*) 'timestep is ',delta,' picoseconds.'
c read start coordinates of atoms
      if(iqfirst.eq.1) then
         read(1,end=500) (temp(i),i=1,natom)
         do 400 i=1,natom
400         x(i)=temp(i)
         read(1,end=500) (temp(i),i=1,natom)
         do 401 i=1,natom
401         y(i)=temp(i)
         read(1,end=500) (temp(i),i=1,natom)
         do 402 i=1,natom
402         z(i)=temp(i)
         iqfirst=0
      end if
c write header of output ascii trajectory
      write(*,*) 'writing header...'
      write(2,1000) hdr
1000  format(a4)
      write(2,1001) io,istart,nsavc,diff,ntitle,natom,nfreat,real(delta) 
1001  format(7i10,e20.8)
      write(2,1010) title(1)
      write(2,1010) title(2)
      write(2,1010) title(3)
1010  format(a80)
      write(2,1012) natom
1012  format(i10)
      if(nfreat.ne.natom) write(2,1020) (freeat(i),i=1,nfreat)
c write start coordinates
      write(*,*) 'writing start coordinates...' 
      write(2,1020) (real(x(i)),i=1,natom)
      write(2,1020) (real(y(i)),i=1,natom)
      write(2,1020) (real(z(i)),i=1,natom)

c loop through rest of trajectory
      write(*,*) 'writing trajectory...'
100   read(1,end=600) (temp(i),i=1,nfreat) 
         do 408 i=1,nfreat
408         x(freeat(i))=temp(i)
      read(1,end=500) (temp(i),i=1,nfreat) 
         do 409 i=1,nfreat
409         y(freeat(i))=temp(i)
      read(1,end=500) (temp(i),i=1,nfreat) 
         do 410 i=1,nfreat
410         z(freeat(i))=temp(i)
c      write(*,*) (x(i),i=1,nfreat)
c      write(*,*) (y(i),i=1,nfreat)
c      write(*,*) (z(i),i=1,nfreat)
      write(2,1020) (real(x(freeat(i))),i=1,nfreat)
      write(2,1020) (real(y(freeat(i))),i=1,nfreat) 
      write(2,1020) (real(z(freeat(i))),i=1,nfreat)
1020  format(10000f7.2)
      ntra=ntra+1
      go to 100

500   write(*,*) 'unexpected end of file !!!!'
600   continue
      nstep=ntra+1
      write(*,*) 'number of frames: ',nstep
      write(*,*) 'edit the ascii-traj !!!! '
      write(*,*) 'put no.of frames as 1st number!'
      write(*,*) '(or give number of frames when'
      write(*,*) 'the program asks...)'
      close(1)
      close(2)

      else
c***** end part xplor --> ascii
c***** start part ascii --> charmm

      ntra=0
      open(1,file=infile,status='old')
      open(2,file=outfile,form='unformatted',status='unknown')
  
      iqfirst=1
c read header of ascii trajectory
      read(1,2001,end=2500) hdr
2001  format(a4)
      read(1,2005,end=2500) nstep,istart,nsavc,diff,ntitle,natom,
     &                      nfreat,rdelta
2005  format(7i10,e20.8)
      read(1,2010,end=2500) title(1)
      read(1,2010,end=2500) title(2)
      read(1,2010,end=2500) title(3)
2010  format(a80)
      read(1,2015,end=2500) natom
2015  format(i10)
      nfreat=natom-diff 
      if(nfreat.ne.natom) then 
          read(1,3020,end=2500) (freeat(i),I=1,nfreat) 
      else
       do 2300 i=1,natom
2300       freeat(i)=i
      end if
c read start coordinates of atoms
      if(iqfirst.eq.1) then
         read(1,3020,end=2500) (rx(i),i=1,natom)
         read(1,3020,end=2500) (ry(i),i=1,natom)
         read(1,3020,end=2500) (rz(i),i=1,natom)
         iqfirst=0
      end if
c write header of output charmm trajectory
      write(*,*) 'number of frames ? '
      read(*,*) nstep
      io=0
      istep=int(rdelta*1000.0)
      write(*,*) 'nstep,diff,natom,istep: ',nstep,diff,natom,istep
      write(2) hdr,nstep,0,istep,io,io,io,io,io,diff,io, 
     &              0,io,io,io,io,io,io,io,io,0
c this is the number of title lines:
      write(2) 0 
      write(2) natom
      nfreat=natom-diff
      if(nfreat.ne.natom) then
         write(2) (freeat(i),i=1,nfreat)
      end if
c write start coordinates 
      write(2) (rx(i),i=1,natom)
      write(2) (ry(i),i=1,natom)
      write(2) (rz(i),i=1,natom)

c loop through rest of trajectory
3100   read(1,3020,end=2600) (rx(freeat(i)),i=1,nfreat) 
      read(1,3020,end=2500) (ry(freeat(i)),i=1,nfreat) 
      read(1,3020,end=2500) (rz(freeat(i)),i=1,nfreat) 
c      write(*,*) (rx(i),i=1,nfreat)
c      write(*,*) (ry(i),i=1,nfreat)
c      write(*,*) (rz(i),i=1,nfreat)
      write(2) (rx(freeat(i)),i=1,nfreat)
      write(2) (ry(freeat(i)),i=1,nfreat) 
      write(2) (rz(freeat(i)),i=1,nfreat)
3020  format(10000f7.2)
      ntra=ntra+1
      go to 3100

2500   write(*,*) 'unexpected end of file !!!!'
2600   continue
      write(*,*) 'number of frames written: ',ntra+1
      close(1)
      close(2)

      end if
      end


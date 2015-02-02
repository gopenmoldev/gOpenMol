*
* I piece of code for moving data from CSD to be displayed by CHEM.
*
* This program reads the FDAT file (output from CSD) and converts it
* to CSSR format for CHEMX or PDB for almost any program.
*
* The input needed is the name of the FADT file and which segment
* in the file will be converted.
*
*
*   Should work with CSD Release 3.2 (1988)
*
*   Leif Laaksonen
*
*   Version 1.0     1988-07-04
*
*   1989-04-17 : PDB format added.
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*   Fortran code for the format of a CSSR file
*    ------------------------------------------
*      Read(Idisk_out,'(22x,A6,10x,3f8.3)')
*     1 Refcode,a,b,c
*      Read(Idisk_out,'(21x,3f8.3)')
*     1 Alpha,Beta,Gamma
*
*      Read(Idisk_out,'(2I4,1x,A60)')
*     1 Iatom,Inorm,Title1
*      Read(Idisk_out,'(I4,I2,I1,A60)')
*     1 Junk1,Iqset,Junk2,Title2
*
*      Do 53 I=1,iatom
*
*      Read(Idisk_out,'(I4,1x,A4,1x,1x,3(F9.5,1x),8I4,1x,F7.3,1x,I3)')
*     1 num(i),Name(i),xpr(i),ypr(i),zpr(i),(Icx(m),M=1,4),
*     2 Ic5,Ic6,Ic7,Ic8,
*     3 Charge(i),IGroup
*
*53    Continue
*
*---------------------------------oo-------------------------------------
*
* Cambridge structural data base format (FDAT)
*
*Card type 1: Directory record
*
*    (1H#,A8,2I1,I6,6X,11I3,22I1,I2)
*
*  Look at page 79 CSD System (3.1) User's Manual Part II
*
*Card type 2: Unit Cell Parameter Record
*
*    (6I6,6I1,6I2,2I3,I3,A8,I3,I24X)
*
*Card Type 3: Text Field
*
*Card Type 4: Symmetry Positions
*
*    (5(3I1,I2,3I1,I2,3I1,I2),5X)
*
*Card Type 5: Radius Values
*
*    (16(A2,I3))
*
*Card Type 6: Atomic Coordinate Records
*
*
*a)  (4(A,4P,3F5.0))
*
*b)  (A,4P,3F7.0,1X,A,3F7.0,1X,A,3F7.0
*
* !!!!!! INORM HAS to be 0 in the CSSR file in this case
*
*Card Type 7: Reported Bond Lengths
*
*    (8(2I3,I4))
*
*Card Type 8: Crystallographic Connectivity Records
*
*a)  If NAT+NSAT < 100 :   (40I2)
*
*b)  If NAT+NSAT > 100 :   (26I3,2X)
*
*
*
*
*  Program to convert CSD FDAT-files to CHEMX CSSR-files
*
      Character Hash*1,Refcod*8,Spg*8,Text(4)*80,Name*5,Rname*2
      Character Input_file*80,Output_file*80,Atom*6,AlternLI*1,
     1          ResidueN*3,ChainI*1,CodeIR*1
      Logical yes_no
      Dimension  Icsd(10),Icelld(6),Iprec(6),Icesd(6),Idens(2)
      Common/GetM/RS(200,3,4),Name(400),Co(400,3),ICN(1000),
     1 Rname(400),Rvalue(400),ESD(6),Icx(400,8),Icp(400)
*
*  Introduce ...
*
      Write(6,'(///A/A/A/A/A/)')
     1'  ************************************************************',
     2'  * This is the CSD FDAT file to CHEMX CSSR or PDB format    *',
     3'  * file conversion program (V1.0).                          *',
     4'  *                                                          *',
     5'  ************************************************************'

881      Write(6,'(A$)') ' Choose from CSSR (1) or PDB (2) file: '
         Read(5,*) Ifile
         if(Ifile.lt.1.or.Ifile.gt.2) GoTo 881
*
800      Write(6,'(/A$)') ' Give input file name : '
         Read(5,'(A)') Input_file
* Check it
         Inquire(file=Input_file,exist=yes_no)
         If(.not. yes_no) Then
         write(6,'(A)') ' %ERROR-Non existent file '
         GoTo 800
         End if
*
         Do 801 i=80,1,-1
          If(Input_file(I:I).ne.' ') GoTo 802
801        Continue
*
802        Long=i
         Do 803 i=Long,1,-1
          If(Input_file(I:I).eq.'.') GoTo 804
803        Continue
804      Iper=i
*       If(Iper.eq.Long+1) Then
*       Output_file=Input_file(1:Long)//'.XI'
*       Else
*       Output_file=Input_file(1:Iper)//'XI'
*       End If
*
*       Write(6,'(A,A)') ' Output file name(s) is/are : ',Output_file     
*
* Check for hash records
*
*    Read Card Type 1: Directory Record
*
      Open(unit=21,file=Input_file)
*
      Write(6,'(A)') ' Scaning input file for REFCODE fields '
      IFDAT=0
810   Read(21,'(A)',End=820) Text(1)
      If(Text(1)(1:1).eq.'#') Then
      IFDAT=IFDAT+1
      Write(6,'(1x,i3,1x,A)') IFDAT,Text(1)(2:9)
      End if
      GoTo 810
820   Continue
* No hash...
      If(IFDAT.eq.0) Then
      Write(6,'(//A//)')
     1' %ERROR-Are you sure this is a FDAT file? (No # in file) '
      Stop
      End if
*
      Write(6,'(/1x,A,A,I3)') 
     1Input_file(1:Long),' contains nr. of FDAT fields : ',IFDAT
      Write(6,'(A$)') ' Which do you want converted (0=EXIT) (nr ?) : '      
      Read(5,*) Iwhich
      If(Iwhich.lt.0.or.Iwhich.gt.IFDAT) GoTo 820
      If(Iwhich.eq.0) Call Exit
*
      Write(6,'(/A,1x,I3)') ' Scaning ... for field nr. ',Iwhich
      Rewind 21
      Loop=0
830   Read(21,'(A)',END=840) Text(1)
      If(Text(1)(1:1).eq.'#') Then
       Loop=Loop+1
        If(Loop.eq.Iwhich) Then
        Do 900 i=9,2,-1
         If(Text(1)(I:I).ne.' ') GoTo 910
900       Continue
         If(Ifile.eq.1) Then
         If(i.eq.1) Output_file='Unknown.XI'
         EndIf
         If(Ifile.eq.2) Then
         If(i.eq.1) Output_file='Unknown.pdb'
         EndIf
910        Lref=i
        If(Ifile.eq.1) Then
        Output_file=Text(1)(2:Lref)//'.XI'
        Write(6,'(A,A)') 
     1' **** Working...  output file name is : ',Output_file
        else
        Output_file=Text(1)(2:Lref)//'.PDB'
        Write(6,'(A,A)') 
     1' **** Working...  output file name is : ',Output_file
        EndIf        
        Backspace 21
        GoTo 840
        End If
       End If
       GoTo 830
840    Continue
*
      Read(21,5)
     1Hash,Refcod,Isys,ICat,IAdat,Ncards,Nrfac,Nrem,Ndis,Nerr,Nopr,
     2Nrad,Nat,Nsat,Nbnd,Ncon,ICell,Intf,IAtfor,ICent,IErr,IRpa,
     3ITd,IPd,Nu,ICbl,IAs,IPol,ICSD,IYear
5     Format(A1,A8,2I1,I6,6X,11I3,22I1,I2)
*
      Write(6,'(/A,A/A,I6/)')
     1' Refcode for this structure is : ',Refcod,
     2' Accession date is             : ',IAdat
*
      If(IATfor.lt.1) Then
      Write(6,'(/A/)')
     1' %WARNING-This structure contains no coordinates !!!'
      End If
* 
*       
*    Read Card Type 2: Unit cell Parameter
*
      Read(21,10)
     1Icelld,Iprec,Icesd,Idens,Nspg,Spg,IZ,Itol
10    Format(6I6,6I1,6I2,2I3,I3,A8,I3,I2)
*
*  Make the needed conversions...
*
      a=float(Icelld(1))/10**Iprec(1)
      b=float(Icelld(2))/10**Iprec(2)
      c=float(Icelld(3))/10**Iprec(3)

      alpha= float(Icelld(4))/10**Iprec(4)
      beta=  float(Icelld(5))/10**Iprec(5)
      gamma= float(Icelld(6))/10**Iprec(6)

*      Write(6,'(A,3(F8.3,1x),/A,3(F8.3,1x))')
*     1' a , b , c : ',a,b,c,' alpha , beta , gamma : ',
*     2alpha,beta,gamma
 
      Do 100 i=1,6
       ESD(i)=float(Icesd(i))*10**Iprec(i)
100     Continue

*           
*    Read Card Type 3: Text Fields
*
      Ntotal=Nrfac+Nrem+Ndis+Nerr
      loop=Ntotal/80+1
      If((((Ntotal/80)*80)-Ntotal).eq.0) loop=loop-1
*
*
      Read(21,'(A)') (Text(i),i=1,loop)

      Write(6,'(A)') ' Text field is:'
      write(6,'(A)') (Text(i),i=1,loop)
*
      If(Nopr.gt.200) Then
       Write(6,'(/A)') ' %ERROR-Dimension for RS is too small'
        Call Exit
         End if
*
      IF(Nopr.gt.0) Then
      READ(21,180) (((RS(N,I,J),J=1,4), I=1,3), N=1,Nopr)
180    FORMAT(5(3F1.0,F2.0,3F1.0,F2.0,3F1.0,F2.0),5X)
        End If
*
      If(Nrad.gt.400) Then 
      Write(6,'(/A)') 
     1 ' %ERROR-Dimension for Rname and Rvalue is too small '
      End if
*
C RADIUS

      IF(Nrad.gt.0) Then
       READ(21,200) (Rname(i),Rvalue(I),I=1,NRad)
200     FORMAT(16(A,F3.0))
         End if
*
C ATOMS
*
      NATD=Nat+Nsat
*                   
      If(Natd.gt.400) Then
      Write(6,'(/A)') 
     1 ' %ERROR-Dimension for Name and Co is too small '
      End if
*
C IN PLUTO THE FIRST 8 ATOMS ARE RESERVED FOR THE CELL CORNER POINTS
*      NAT=NATD+8
*      IF(NATD.EQ.0) GO TO 240
C ATOM FORMAT 1
      IF(IAtfor.EQ.1) Then 
      READ(21,220)(NAME(I),(Co(I,K),K=1,3),I=1,NATD)
*
*       Do 907 I=1,3
*        Do 907 J=1,Natd
*907      Co(j,i)=Co(j,i)/10000.
      End If
C ATOM FORMAT 2
      IF(IAtfor.EQ.2) Then
      READ(21,230)(NAME(I),(Co(I,K),K=1,3),I=1,NATD)
*
*       Do 904 I=1,3
*        Do 904 J=1,Natd
*904      Co(j,i)=Co(j,i)/100000.
      End If
*  220 FORMAT(4(A,3F5.0))
*  230 FORMAT(A,3F7.0,1X,A,3F7.0,1X,A,3F7.0)
  220 FORMAT(4(A,4P,3F5.0))
  230 FORMAT(A,5P,3F7.0,1X,A,3F7.0,1X,A,3F7.0)
C BONDS
      IF(NBND.EQ.0) GO TO 260
C SKIP THE BOND CARDS  - CALCULATE THE NUMBER OF CARDS
      IF(NATD.GT.0) NC=(NBND+7)/8
      IF(NATD.EQ.0) NC=(NBND+4)/5
      DO 250 I=1,NC
  250 READ(21,220)
C CONNECTIVITY
260   Continue
*
      If(NCON.gt.1000) Then
      Write(6,'(/A)')
     1' %ERROR-Dimension for ICN is too small'
      End if
*
      IF(NCON.gt.0) Then
       IF(NATD.LT.100) READ(21,270)(ICN(I),I=1,NCON)
        IF(NATD.GE.100) READ(21,280)(ICN(I),I=1,NCON)
  270    FORMAT(40I2)
  280    FORMAT(26I3,2X)
      End If
*
*  Build connectivity list
*
      Idisk_out=22
      Open(Unit=Idisk_out,File=Output_file,Status='NEW')
*
      Do 420 i=1,Natd
*
       Do 423 j=1,8
423     Icx(i,j)=0
*
420    Icp(i)=1
*
      Do 400 i=1,Natd
      If(Icn(i).eq.0) GoTo 400
       Icx(i,Icp(i))=Icn(i)
        Icp(i)=Icp(i)+1
         Icx(Icn(i),Icp(Icn(i)))=i
          Icp(Icn(i))=Icp(Icn(i))+1
400      Continue
      Do 410 i=Natd+1,Ncon,2
*
*      Icn(i) is bond to Icn(i+1)
*
      If(Icn(i).eq.0.or.Icn(i+1).eq.0) GoTo 410 
*
      Icx(Icn(i),Icp(Icn(i)))=Icn(i+1)
      Icp(Icn(i))=Icp(Icn(i))+1
*
      Icx(Icn(i+1),Icp(Icn(i+1)))=Icn(i)
      Icp(Icn(i+1))=Icp(Icn(i+1))+1
410   Continue
C
C CHECK FOR PADDING 0 AT END OF THIS RECORD
*      IF (ITO(NCON).EQ.0) THEN
*        IDC(11)=IDC(11)-1
*      ENDIF
C END OF FDAT ENTRY
         If(Ifile.eq.1) Then
*
*   Write out CSSR file
*
*   Fortran code for the format of a CSSR file
*    ------------------------------------------
      Iqset=0
      Inorm=0
      Write(Idisk_out,'(22x,A6,10x,3f8.3)')
     1 Refcod,a,b,c
      Write(Idisk_out,'(21x,3f8.3)')
     1 Alpha,Beta,Gamma

      Write(Idisk_out,'(2I4,1x,A60)')
     1 Natd,Inorm,Text(1)
      Write(Idisk_out,'(I4,I2,I1,A60)')
     1 Junk1,Iqset,Junk2,Text(2)

      Ic5=0
       Ic6=0
        Ic7=0
         Ic8=0
      Charge=0.0
      Igroup=1

      Do 53 I=1,Natd

      Write(Idisk_out,'(I4,1x,A4,1x,1x,3(F9.5,1x),8I4,1x,F7.3,1x,I3)')
     1 I,Name(i),Co(i,1),Co(i,2),Co(i,3),(Icx(i,m),M=1,8),
     2 Charge,IGroup

53    Continue
      Else 
* PDB format
*ATOM      9  C1          1       0.000   0.000   0.000
*ATOM      1 N    GLN     1       0.800 -14.340  51.470
      Write(Idisk_out,'(A)') 'HEADER'
      ISN=1
      IfootN=0
      Occ=0.0
      TempF=0.0
      ResidueN='MOL'
      CodeIR=' '
      Atom='ATOM'
      ChainI=' '
      AlternLI=' '
      Do 554 i=1,Natd
      Write(Idisk_out,555) 
     1Atom,I,Name(i)(1:3),AlternLI,ResidueN,ChainI,ISN,CodeIR,
     1 a*Co(i,1),
     3 b*Co(i,2),c*Co(i,3),
     2Occ,TempF,IFootN
555    Format(A,I5,2x,A,A1,A,1x,A,I4,A,3x,3F8.3,2f6.2,1x,I3)
554   Continue
      EndIf
      Close(unit=Idisk_out)
      Write(6,'(/A,A,A/)') ' ====> ',Refcod,' is converted '
      GoTo 820
     
      End

    

                   PROGRAM VSSMOD    
C---ORIGINAL PROGRAMME WRITTEN BY C. GIESSNER-PRETTRE.                      
C
C--MODIFED BY  SAL FEB 1986
C
C   This version of the main block is considerably different to the 
C   original VSS programme. Code for the production of line-printer
C   plots has been removed. The remaining code has been relaid and 
C   modified to give potentials at user defined intervals along 
C   all the axes.
C
C
C-Modified by SAL JULY 1986 to allow finer grid and adjustable array
C- for the potentials
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                       
      INCLUDE '(QCDLFL)'
      INCLUDE '(VSSDIM)'
      INCLUDE '(VSSMPD)'
      DIMENSION AR3POT(NVSSPX,NVSSPY,NVSSPZ)
      COMMON/VSSCOM/ X(NVSSTM),Y(NVSSTM),Z(NVSSTM),IZ(NVSSTM),Q(NVSSTM),
     1 VS, NAT,NAP,ZET(NVSSTM)
      COMMON/PMAP/NPTS(3),XMIN(3),XMAX(3),RINC(3),AUTANG,IQMAP
      CHARACTER CTITLE*80,CQMRUN*25                            
      CALL CDSRT2
      AUTANG=0.52917715
      ANGTAU=1/AUTANG
      RMINPT=0.0
      CBINFL=' '
      CSSRFL=' '
      CTITLE=' '
C---Read binary filename
      READ(ICDLIN,150)CBINFL(1:60)
      IF(INDEX(CBINFL(1:60),' ').NE.0)GOTO 100
      READ(ICDLIN,'(A)')CBINFL(61:120)
      IF(INDEX(CBINFL(61:120),' ').NE.0)GOTO 100
      READ(ICDLUT,'(A)')CBINFL(121:200)
      IF(INDEX(CBINFL(121:200),' ').NE.0)GOTO 100
      READ(ICLUT,'(A)')CBINFL(201:)
C---Read CSSR filename
100   READ(ICDLIN,150)CSSRFL(1:60)
      IF(INDEX(CBINFL(1:60),' ').NE.0)GOTO 200
      READ(ICDLIN,'(A)')CSSRFL(61:120)
      IF(INDEX(CSSRFL(61:120),' ').NE.0)GOTO 200
      READ(ICDLIN,'(A)')CSSRFL(121:200)
      IF(INDEX(CSSRFL(121:200),' ').NE.0)GOTO 200
      READ(ICDLIN,'(A)')CSSRFL(201:)
150   FORMAT(19X,A60)                        
C---Read title.
200   READ(ICDLIN,'(A)')CTITLE
      WRITE(ICDLUT,'(1X,A)')CTITLE                                  
C---Read total number of atoms, number of hydrogens and units flag.
      READ(ICDLIN,*)NAT,NH,IU
C---Check total number of atoms does not exceed allowed maximum.
      IF(NAT.GT.NVSSTM)THEN
        WRITE(MOFILE,250)
        GOTO 999
      ENDIF
250   FORMAT(' Number of atoms exceeds allowed maximum ')               
       NAP=NAT-NH                                            
C---Read atom coordinates and number of valence electrons
      READ(ICDLIN,*)(X(I),Y(I),Z(I),IZ(I),I=1,NAT)                              
C--- If necessary convert angstom units to atomic units 
      IF(IU.EQ.1)THEN                                                          
        DO 5 I=1,NAT                                                    
        X(I)=X(I)*ANGTAU                                                
        Y(I)=Y(I)*ANGTAU                                                
        Z(I)=Z(I)*ANGTAU                                                
    5 CONTINUE                                                                 
      ENDIF
      WRITE(ICDLUT,1020)(I,X(I),Y(I),Z(I),IZ(I),I=1,NAT)  
 1020 FORMAT(I4,3F16.6,I4)                                      
C---Read gross atomic populations
      READ(ICDLIN,*)(Q(I),I=1,NAT)                                              
C---Read  atomic orbital screening constants.
      READ(ICDLIN,*)(ZET(I),I=1,NAT)                                            
      WRITE(ICDLUT,1030)(Q(I),I=1,NAT)                    
1030  FORMAT (10X,10F10.5)
C---  Read type of map
      READ(ICDLIN,*)IQMAP  
C---Read coordinate bounds and number of intervals along the axes.
      READ(ICDLIN,*) (XMIN(I),XMAX(I),I=1,3),NPTS
      WRITE(ICDLUT,1040) (XMIN(I),XMAX(I),I=1,3), NPTS
1040  FORMAT(10X,' X BOUNDS',2F6.1,' Y BOUNDS', 2F6.1,' Z BOUNDS',            
     1 2F6.1 / ' NPTSX',I5,' NPTSY',I5,' NPTSZ',I5)
C--- If the array length is exceeded print error message.
      IF(IQMAP.EQ.1.AND.NPTS(1)*NPTS(2).GT.64000 .OR.
     1   IQMAP.EQ.2.AND.NPTS(1)*NPTS(2)*NPTS(3).GT.64000)THEN
         WRITE(MOFILE,1045)
         GOTO 999
      ENDIF 
1045  FORMAT(1X,'Potential array bounds exceeded in VSS')        
      IF(IQMAP.EQ.1)NPTS(3)=1
C---Calculate step-lengths 
        DO 6 I=1,3        
          IF(NPTS(I).GT.1)THEN
            RINC(I)=(XMAX(I)-XMIN(I))/(NPTS(I)-1)
          ELSE
            RINC(I)=0.0
          ENDIF
    6   CONTINUE
C-Generate the binary map file and listing file
      CALL SLPMAP(AR3POT,0)
      IF(IERR.EQ.1)GOTO 999
      CQMRUN='Electrostatic potential'
      CALL QCDLHS(CTITLE,'VSS',CQMRUN,' ',' ',' ')    
      CALL QCDLFN
999   STOP                                                           
      END
C                                                                       
C
      SUBROUTINE UOTM(XH,YH,ZH,POT)                                             
      IMPLICIT REAL*8(A-H,O-Z)                                                  
      INCLUDE '(VSSDIM)'
      COMMON/VSSCOM/ X(NVSSTM),Y(NVSSTM),Z(NVSSTM),IZ(NVSSTM),Q(NVSSTM),
     1 VS,NAT,NAP,ZET(NVSSTM)                  
      POT=0.                                                                    
      REP=0.                                                                    
      DO 10 I=1,NAT                                                             
      XI=X(I)                                                                   
      YI=Y(I)                                                                   
      ZI=Z(I)                                                                   
      RO=DSQRT((XH-XI)**2+(YH-YI)**2+(ZH-ZI)**2)                                
      IF(RO-.01)25,25,20                                                        
   25 RO=0.                                                                     
      POT=0.                                                                    
      GO TO 100                                                                 
   20 CONTINUE                                                                  
      IF(I-NAP)4,4,5                                                            
    4 CONTINUE                                                                  
      CALL UOTA2(RO,ZET(I))                                                     
      POT=POT-VS*Q(I)                                                           
      GO TO 8                                                                   
    5 CONTINUE                                                                  
      CALL UOTA1(RO,ZET(I))                                                     
      POT=POT-VS*Q(I)                                                           
    8 CONTINUE                                                                  
      IF(RO)10,10,9                                                             
    9 REP=REP+IZ(I)/RO                                                          
   10 CONTINUE                                                                  
      POT=POT+REP                                                               
  100 RETURN                                                                    
      END                                                                       
C
C
      SUBROUTINE UOTA2(R,Z)                                                     
      IMPLICIT REAL*8(A-H,O-Z)                                                  
      INCLUDE '(VSSDIM)'
      COMMON/VSSCOM/ X(NVSSTM),Y(NVSSTM),W(NVSSTM),IZ(NVSSTM),Q(NVSSTM),
     1 VS,NAT,NAP,ZET(NVSSTM)                  
      IF(R)5,2,5                                                                
    2 CONTINUE                                                                  
      VS=Z/2.                                                                   
      GO TO 50                                                                  
    5 CONTINUE                                                                  
      RO=Z*R                                                                    
      DERO=DEXP(-2.*RO)                                                         
      VS=(1./R)*(1.-(1.+RO*(1.5+RO*(1.+RO/3.)))*DERO)                           
   50 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C
C
      SUBROUTINE UOTA1(R,Z)                                                     
      IMPLICIT REAL*8(A-H,O-Z)                                                  
      INCLUDE '(VSSDIM)'
      COMMON/VSSCOM/ X(NVSSTM),Y(NVSSTM),W(NVSSTM),IZ(NVSSTM),Q(NVSSTM),
     1 VS,NAT,NAP,ZET(NVSSTM)                  
      IF(R)5,2,5                                                                
    2 CONTINUE                                                                  
      VS=Z                                                                      
      GO TO 50                                                                  
    5 CONTINUE                                                                  
      RO=Z*R                                                                    
      DERO=DEXP(-2.*RO)                                                         
      VS=(1./R)*(1.-(1.+RO)*DERO)                                               
   50 CONTINUE                                                                  
      RETURN                                                                    
      END

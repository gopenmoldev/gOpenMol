*
*  This version supports the input generation to the QCPE program VSS
*
*       Leif Laaksonen 1992
*
*
************************************************************************
*
C   =517    FORTICON8 (VAX VERSION)
C**********************************************************************
C                                                                     *
C    PROGRAM FORTICON0  COMPLETE FORTRAN VERSION OF ICON8             *
C                                                                     *
C    THE FOLLOWING SUBPROGRAMS, WHICH EXIST AS ASSEMBLER ROUTINES     *
C    IN ICON8, ARE TRANSLATED INTO FORTRAN0  MATRIX, ABFNS,           *
C    LOVLAP, GRMSCH, TRNFRM, DSUM, ROTATE, DOT, VECSUM, REDUCE,       *
C    FULCHM, AND REDCHM.  FORTICON INLUDES THESE AS WELL AS ALL       *
C    THE FORTRAN SUBPROGRAMS OF ICON8.                                *
C                                                                     *
C**********************************************************************
C
           IMPLICIT REAL*8(A-H,O-Z)
C
C     PROGRAM ICON FOR PERFORMING EXTENDED HUCKEL CALCULATIONS
C     WITH OR WITHOUT CHARGE ITERATION.
C     ** QCPE VERSION **
C
C     ** SAMPLE DECK **
C ....0....1....0....2....0....3....0....4....0....5....0....6....0
C
C ETHYLENE
C   4  2        2
C  0.92665        1.205          0.0
C -0.92665        1.205          0.0
C  0.92665       -1.205          0.0
C -0.92665       -1.205          0.0
C  0.0            0.67           0.0
C  0.0           -0.67           0.0
C  C C
C
C
           LOGICAL*4 L1,L2,L3,L4,L5,ONEMAT,ITERAT
           LOGICAL*4 PRT,PUN
           INTEGER*4 IOVPOP,IENRGY
           INTEGER*4 AC,SYMBOL,VELEC
           REAL*8 LAMPRI
           INTEGER*4 PRTCYC
           LOGICAL*4 PARTIT,PRINTX,ITABLE
           INTEGER*4 STAR,STAR2
           COMMON/TITLE/AB(10)
           COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,
     .     IPRINT,IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT
           COMMON/OUT/PRT(20),PUN(20),IOVPOP(24),IENRGY(24)
           COMMON/ATOM/AC(1000),SYMBOL(100),VELEC(40),NS(40),NP(40),
     .     ND(40),EXPS(40),EXPP(40),EXPD(40),EXPD2(40),C1(40),C2(40),
     .     COULS(40),COULP(40),COULD(40),X(1000),Y(1000),Z(1000)
           COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,
     .     PRTCYC,ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)
           COMMON/ABC/AS2(5),BS2(5),CS2(5),AP2(5),BP2(5),CP2(5),AD2(5),
     .     BD2(5),CD2(5),AS3(5),BS3(5),CS3(5),AP3(5),BP3(5),CP3(5),
     .     AD3(5),BD3(5),CD3(5)
           COMMON/STARS/STAR,STAR2
C
C     CALL INPUT TO READ IN AND PRINT OUT INPUT DATA.
C
1000       CALL INPUT(NATOM,NDIM,NTYPE)
           IF(IPRINT.LT.-3) GO TO 2000
C
C     CALCULATE MATRIX DIMENSIONS.
C
           NC=NDIM*(NDIM+1)/2
           NHS=NC+NC-NDIM
           NSS=NHS
C                                   
C     IF NO WAVEFUNCTIONS NEEDED, SET NHS=0. THIS HAS THE EFFECT
C     OF EQUIVALENCING THE H AND S MATRICES.
C
           ONEMAT=.FALSE.
           IF(.NOT.ITERAT.AND.IPRINT.LT.-1) ONEMAT=.TRUE.
           IF(IPUNCH.NE.0) GO TO 600
           DO 400 I=6,20
           IF(PUN(I)) GO TO 600
400        CONTINUE
           GO TO 500
600        ONEMAT=.FALSE.
500        IF(ONEMAT) NHS=0
           IF(METH.LT.3) NTYPE=1
           NMD=NTYPE*NTYPE
           NCL=NDIM
           IF(NDIM.LT.10) NCL=10
           NOCC=(NDIM+1)/2
           NHDG=1
           IF(METH.GT.2.AND.L5) NHDG=NDIM
C
C     CALL MATRIX TO ALLOCATE SPACE FOR MATRICES.
C     ORDER OF MATRICES0 H S MAD C SP PD MAXS MAXP MAXD
C                        COUL0 SORB IOCC HDG
C
200        CALL MATRIX(13,NHS,NSS,NMD,NC,NDIM,NDIM,NDIM,NDIM,2*NDIM,
     .     NCL,NDIM,NOCC,NHDG,     NDIM,NDIM,NC,NATOM,NTYPE,NHDG,NC,NC,
     .     NC,NC,NC,NC,NDIM)
2000       CONTINUE
*           GO TO 1000
           END
           BLOCK DATA
C
C     INITIALIZATION OF INTERNAL ATOMIC DATA. THERE ARE PROVISIONS
C     FOR 20 USER DEFINED ATOMS, 15 INTERNALLY DEFINED ATOMS, AND
C     SPACE FOR 5 MORE TO BE USED EITHER WAY.
C
           IMPLICIT REAL*8(A-H,O-Z)
           INTEGER*4 STAR,STAR2
           INTEGER*4 KEY,SYMBOL,VELEC,SymbH
           COMMON/ATOM/KEY(1000),SYMBOL(100),VELEC(40),NS(40),NP(40),
     .     ND(40),EXPS(40),EXPP(40),EXPD(40),EXPD2(40),C1(40),C2(40),
     .     COULS(40),COULP(40),COULD(40),X(1000),Y(1000),Z(1000)
           COMMON/STARS/STAR,STAR2
           Common/LUL/IhelpV,SymbH,XX,YY,ZZ,MaxDim
           COMMON/START/NUSER
*           DATA SYMBOL/20*'**',' C',' N',' O',' F','SI',' P',' S',
*     .     'CL','LI','BE',' B','NA','MG','AL',' H',5*'  '/
           DATA SYMBOL/20*'**',' C',' N',' O',' F','Si',' P',' S',
     .     'Cl','Li','Be',' B','Na','Mg','Al',' H',
     &           'CA','CB','CG','CD','CE','CZ',
     &           'NA','NB','NG','ND','NE','NZ',
     &           'OA','OB','OG','OD','OE','OZ','OH',
     &           'SA','SB','SG','SD',
     &           'HH','HG','HZ','HE','NH','CH','OC',
     &           'HA','HB','HD','HN','Ca',30 * '  '/
      Dimension IhelpV(100),SymbH(1000),XX(1000),YY(1000),ZZ(1000)
      Data IhelpV /20*0,
     &             21,22,23,24,25,26,27,
     &             28,29,30,31,32,33,34,35,
     &             21,21,21,21,21,21,
     &             22,22,22,22,22,22,
     &             23,23,23,23,23,23,23,
     &             27,27,27,27,
     &             35,35,35,35,22,21,23,
     &             35,35,35,35,36,30*0/
      Data SymbH/1000*'**'/
           DATA VELEC/20*0,4,5,6,7,4,5,6,7,1,2,3,1,2,3,1,2,4*0/
********
           DATA NS/20*0,4*2,4*3,3*2,3*3,1,4,4*0/
           DATA EXPS/20*0.0D0,1.625D0,1.950D0,2.275D0,2.425D0,1.383D0,
     .     1.6D0,1.817D0,2.033D0,0.65D0,0.975D0,1.3D0,0.733D0,0.95D0,
     .     1.167D0,1.3D0,1.2D0,4*0.0D0/
           DATA COULS/20*0.0D0,-21.4D0,-26.0D0,-32.3D0,-40.0D0,-17.3D0,
     .     -18.6D0,-20.0D0,-30.0D0,-5.4D0,-10.0D0,-15.2D0,-5.1D0,-9.0D0,
     .     -12.3D0,-13.6D0,-7.0D0,4*0.0D0/
           DATA NP/20*0,4*2,4*3,3*2,3*3,0,4,4*0/
           DATA EXPP/20*0.0D0,1.625D0,1.950D0,2.275D0,2.425D0,1.383D0,
     .     1.6D0,1.817D0,2.033D0,0.65D0,0.975D0,1.3D0,0.733D0,0.95D0,
     .     1.167D0,0.0,1.2D0,4*0.0D0/
           DATA COULP/20*0.0D0,-11.4D0,-13.4D0,-14.8D0,-18.1D0,-9.2D0,
     .     -14.0D0,-13.3D0,-15.0D0,-3.5D0,-6.0D0,-8.5D0,-3.0D0,-4.5D0,
     .     -6.5D0,0.0,-4.0,4*0.0D0/
           DATA ND/24*0,4*3,12*0/
           DATA EXPD/24*0.0D0,1.383D0,1.4D0,1.5D0,2.033D0,12*0.0D0/
           DATA COULD/24*0.0D0,-6.0D0,-7.0D0,-8.0D0,-9.0D0,12*0.0D0/
*
***********
*           DATA NS/20*0,4*2,4*3,3*2,3*3,1,5*0/
*           DATA EXPS/20*0.0D0,1.608D0,1.924D0,2.246D0,2.564D0,1.634D0,
*     .     1.881D0,2.122D0,2.356D0,0.65D0,0.975D0,1.3D0,0.733D0,0.95D0,
*     .     1.167D0,1.3D0,5*0.0D0/
*           DATA COULS/20*0.0D0,-19.2D0,-25.7D0,-33.9D0,-42.8D0,-14.7D0,
*     .     -18.9D0,-29.2D0,-30.0D0,-5.4D0,-10.0D0,-15.2D0,-5.1D0,-9.0D0,
*     .     -12.3D0,-13.6D0,5*0.0D0/
*           DATA NP/20*0,4*2,4*3,3*2,3*3,6*0/
*           DATA EXPP/20*0.0D0,1.568D0,1.917D0,2.227D0,2.55D0,1.428D0,
*     .     1.628D0,1.827D0,2.039D0,0.65D0,0.975D0,1.3D0,0.733D0,0.95D0,
*     .     1.167D0,6*0.0D0/
*           DATA COULP/20*0.0D0,-11.8D0,-15.4D0,-17.2D0,-19.9D0,-8.1D0,
*     .     -10.7D0,-11.9D0,-13.8D0,-3.5D0,-6.0D0,-8.5D0,-3.0D0,-4.5D0,
*     .     -6.5D0,6*0.0D0/
*           DATA ND/24*0,4*3,12*0/
*           DATA EXPD/24*0.0D0,1.383D0,1.4D0,1.5D0,2.033D0,12*0.0D0/
*           DATA COULD/24*0.0D0,-6.0D0,-7.0D0,-8.0D0,-9.0D0,12*0.0D0/
*
*********
           DATA C1/40*0.0D0/
           DATA C2/40*0.0D0/
           DATA EXPD2/40*0.0D0/
           DATA STAR,STAR2 / ' *','**'/
           DATA NUSER/21/
           Data MaxDim/1000/
           END
           SUBROUTINE INPUT(NATOM,NDIM,NTYPE)
C
C     SUBROUTINE FOR READING IN AND PRINTING OUT INPUT DATA.
C
           IMPLICIT REAL*8(A-H,O-Z)
           LOGICAL*4 PRT,PUN
           INTEGER*4 IOVPOP,IENRGY
           INTEGER*4 AC,SYMBOL,VELEC,SymbH
           LOGICAL*4 L1,L2,L3,L4,L5,ONEMAT,ITERAT
           REAL*8 LAMPRI
           INTEGER*4 PRTCYC
           LOGICAL*4 PARTIT,PRINTX,ITABLE
           INTEGER*4 CHANGE(20)
           REAL*4 MADS(20),MADP(20),MADD(20)
           INTEGER*4 STAR,STAR2
           INTEGER*4 CONTIN
           INTEGER*4 HYDROG
           COMMON/TITLE/AB(10)
           COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,
     .     IPRINT,IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT
           Common/LUL/IhelpV(100),SymbH(1000),XX(1000),YY(1000),
     .     ZZ(1000),MaxDim
           COMMON/OUT/PRT(20),PUN(20),IOVPOP(24),IENRGY(24)
           COMMON/ATOM/AC(1000),SYMBOL(100),VELEC(40),NS(40),NP(40),
     .     ND(40),EXPS(40),EXPP(40),EXPD(40),EXPD2(40),C1(40),C2(40),
     .     COULS(40),COULP(40),COULD(40),X(1000),Y(1000),Z(1000)
           COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,
     .     PRTCYC,ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)
C
C     SINCE THE ITERNAL ATOMIC PARAMETERS ( EXPS, EXPP, ETC. ) ARE
C     NOT USED WHEN DOING CHARGE ITERATION ( METH >1 ) THE SPACE
C     ALLOCATED TO THEM CAN BE USED FOR THE VSIE CHARGE ITERATION
C     PARAMETERS.
C
           DIMENSION AS1(20),BS1(20),CS1(20),AP1(20),BP1(20),CP1(20),
     .     AD1(20),BD1(20),CD1(20)
           EQUIVALENCE (AS1(1),EXPS(21)),(BS1(1),EXPP(21)),
     .     (CS1(1),EXPD(21)),(AP1(1),EXPD2(21)),(BP1(1),C1(21)),
     .     (CP1(1),C2(21)),(AD1(1),COULS(21)),(BD1(1),COULP(21)),
     .     (CD1(1),COULD(21))
           COMMON/ABC/AS2(5),BS2(5),CS2(5),AP2(5),BP2(5),CP2(5),AD2(5),
     .     BD2(5),CD2(5),AS3(5),BS3(5),CS3(5),AP3(5),BP3(5),CP3(5),
     .     AD3(5),BD3(5),CD3(5)
           EQUIVALENCE (CHANGE(1),AS2(1))
           EQUIVALENCE (MADS(1),NS(21)),(MADP(1),NP(21)),
     .     (MADD(1),ND(21))
           COMMON/STARS/STAR,STAR2
           COMMON/START/NUSER
           DIMENSION EXTRA(9)
           EQUIVALENCE (X(1),EXTRA(1))
           EQUIVALENCE (AB(10),CONTIN)
           DATA HYDROG/' H'/
C
C     READ AND WRITE TITLE.
C     IF CONTIN IS EQUAL TO STAR THEN ANOTHER TITLE CARD WILL
C     BE READ AND PRINTED. HOWEVER ONLY THE FIRST IS STORED
C     FOR PRINTING LATER ON. GET YOUR GOODIES ON THE FIRST.
C
1000       READ(5,1,END=115) AB
1          FORMAT(8A8,A6,A2)
           WRITE(6,2) AB
2          FORMAT('1',T10,8A8,A6,A2)
11         IF(CONTIN.NE.STAR) GO TO 9
           READ(5,1) EXTRA,CONTIN
           WRITE(6,12) EXTRA,CONTIN
12         FORMAT(T10,8A8,A6,A2)
           GO TO 11
9          CONTINUE
C
C     READ PARAMETER CARD.
C
           READ(5,3) NH,NA,KA,METH,IPRINT,IPUNCH,L1,L2,L3,L4,L5,CON,
     .     PEEP,COULH,(PRT(I),I=1,20),(PUN(J),J=1,20)
3          FORMAT(6I3,5L1,F5.2,2F6.3,40L1)
C
C Check if dimenions are allowed
           If((na+nh) .gt. MaxDim) Stop 'Too many atoms '
C
C     INSERT DEFAULT PARAMETERS.
C
           IF(CON.LT.1.E-05) CON=1.75D0
           IF(PEEP.LT.1.E-05) PEEP=1.3D0
           IF(COULH.GT.-1.E-05) COULH=-13.6D0
           ITERAT=METH.NE.0
           NATOM=NH+NA
           NATM=NATOM
C
C     SET IPRINT OPTION.
C
           IF(IPRINT.GT.1) GO TO 250
           PRT(6)=.TRUE.
           PRT(7)=.TRUE.
           PRT(11)=.TRUE.
           PRT(17)=.TRUE.
           PRT(19)=.TRUE.
           IF(IPRINT.GT.0) GO TO 250
           PRT(12)=.TRUE.
           PRT(14)=.TRUE.
           PRT(20)=.TRUE.
           IF(IPRINT.GT.-1) GO TO 250
           PRT(13)=.TRUE.
           PRT(15)=.TRUE.
           PRT(16)=.TRUE.
           PRT(18)=.TRUE.
           IF(IPRINT.GT.-2) GO TO 250
           PRT(10)=.TRUE.
C
C     READ COORDINATES AND HEAVY ATOM CARD.
C
250        READ(5,5) (X(I),Y(I),Z(I),I=1,NATOM)
5          FORMAT(3F15.6)
           READ(5,8) (AC(I),I=1,NA)
8          FORMAT(40A2)
*
* Save coordinates for later use
*
           Do 8000 i=1,Natom
           XX(i)=X(i)
            YY(i)=Y(i)
8000         ZZ(i)=Z(i)
C
C     READ AND DECODE ATOM DEFINITION CARDS.
C
           NDIM=NH
           NTYPE=NH
           NELEC=NH-KA
           K=NUSER
*           NUSER2=40
* NUSER2 has to be the same as IhelpV or SYMBOL is defined
*
*
*  Max atom labels in the list to be searched
*
           NUSER2 = 100
           IF(METH.GE.2) NUSER2=20
           DO 100 I=1,NA
            SymbH(i)=Ac(i)
           IF(NUSER.GT.NUSER2) GO TO 103
           DO 102 J=NUSER,NUSER2
*           JSAVE=J
           JSAVE=IhelpV(j)
           IF(AC(I) .EQ. SYMBOL(J)) GO TO 101
102        CONTINUE
C
C     PROVISION FOR USER SPECIFIED DATA.
C
           IF(AC(I).EQ.STAR) GO TO 103
           IF(AC(I).EQ.STAR2) GO TO 105
           WRITE(6,6) I,AC(I)
6          FORMAT(//,T10,'HEAVY ATOM',I3,' NOT RECOGNIZED. SYMBOL0',A2)
           IF(METH.GE.2) WRITE(6,13)
13         FORMAT(/,T10,'REMEMBER IF USING METH > 1 ALL ATOMIC',
     .     ' PARAMETERS MUST BE DEFINED BY THE USER.')
115        REWIND 7
           STOP
103        NUSER=NUSER-1
105        READ(5,7) SYMBOL(NUSER),VELEC(NUSER),NS(NUSER),EXPS(NUSER),
     .     COULS(NUSER),NP(NUSER),EXPP(NUSER),COULP(NUSER),ND(NUSER),
     .     EXPD(NUSER),COULD(NUSER),C1(NUSER),EXPD2(NUSER),C2(NUSER)
7          FORMAT(A2,I3,3(I3,2F6.3),F6.4,F6.3,F6.4)
           JSAVE=NUSER
C
C     NORMALIZE USER SPECIFIED CONTRACTED D ORBITAL.
C
           IF(C2(NUSER).EQ.0.) GO TO 101
           S=(4.D0*EXPD(NUSER)*EXPD2(NUSER)/(EXPD(NUSER)+EXPD2(NUSER)
     .     )**2)**(ND(NUSER)+.5D0)
           S=1.D0/SQRT(C1(NUSER)**2+C2(NUSER)**2+(S+S)*C1(NUSER)
     .     *C2(NUSER))
           C1(NUSER)=S*C1(NUSER)
           C2(NUSER)=S*C2(NUSER)
101        NELEC=NELEC+VELEC(JSAVE)
C
C     AC, LATER REFERENCED AS KEY, IS A POINTER TO THE PARAMETER TABLES.
C
111        AC(I)=JSAVE
           NDIM=NDIM+4
           IF(NP(JSAVE).EQ.0) NDIM=NDIM-3
           IF(ND(JSAVE).NE.0) NDIM=NDIM+5
           NTYPE=NTYPE+2
           IF(NP(JSAVE).EQ.0) NTYPE=NTYPE-1
           IF(ND(JSAVE).NE.0) NTYPE=NTYPE+1
100        CONTINUE
C
C     READ IN CHARGE ITERATION PARAMETERS. SET DEFAULT VALUES.
C
           IF(.NOT.ITERAT) GO TO 60
           IF(K.NE.21) GO TO 60
           READ(5,61) DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,
     .     PRTCYC,NCON,PARTIT
61         FORMAT(6F10.5,3I5,4X,L1)
           IF(.NOT.PARTIT.OR.METH.EQ.2) GO TO 65
           WRITE(6,66)
66         FORMAT(///,T10,'PARTIAL ITERATION ( PARTIT = TRUE ) MAY',
     .     ' ONLY BE USED IF METH = 2.')
           STOP
65         IF(DELTAC.EQ.0.0D0) DELTAC=0.0001D0
           IF(SENSE.EQ.0.0D0) SENSE=2.0D0
           IF(MAXCYC.EQ.0) MAXCYC=25
           IF(PRTCYC.EQ.0) PRTCYC=MAXCYC
           IF(NCON.EQ.0) NCON=3
           IF(DAMP1.EQ.0.0D0) DAMP1=0.1D0
           IF(METH.GE.3) GO TO 62
           IF(DAMP2.EQ.0.0D0) DAMP2=0.25D0
           IF(LAMPRI.EQ.0.0D0) LAMPRI=0.25D0
           GO TO 63
62         IF(DAMP2.EQ.0.0D0) DAMP2=0.75D0
           IF(LAMPRI.EQ.0.0D0) LAMPRI=0.75D0
63         IF(METH.LT.2) GO TO 60
           DO 32 I=1,20
32         ITABLE(I)=.FALSE.
           NUSER2=21-NUSER
C
C     READ IN SYMBOLS OF ATOMS ON WHICH CHARGE ITERATION
C     IS TO BE PERFORMED.
C
           IF(.NOT.PARTIT) GO TO 30
           READ(5,31) CHANGE
31         FORMAT(20A2)
           DO 33 I=1,NUSER2
           J=21-I
           DO 33 K=1,NUSER2
           IF(SYMBOL(J).EQ.CHANGE(K)) ITABLE(J)=.TRUE.
33         CONTINUE
           GO TO 34
30         DO 35 I=1,NUSER2
           J=21-I
35         ITABLE(J)=.TRUE.
C
C     READ IN VSIE AND MADELUNG PARAMETERS.
C
34         DO 36 I=1,NUSER2
           J=21-I
           IF(.NOT.ITABLE(J)) GO TO 36
           READ(5,37) AS1(I),BS1(I),CS1(I),MADS(I)
37         FORMAT(4F10.8)
           IF(NP(J).EQ.0) GO TO 36
           IF(ND(J).NE.0) GO TO 38
           READ(5,37) AP1(I),BP1(I),CP1(I),MADP(I)
           GO TO 36
38         IF(NCON.EQ.3) GO TO 39
           READ(5,37) AP1(I),BP1(I),CP1(I),MADP(I),AD1(I),BD1(I),
     .     CD1(I),MADD(I)
           GO TO 36
39         READ(5,40) AS2(I),BS2(I),CS2(I),AS3(I),BS3(I),CS3(I)
40         FORMAT(3F10.8)
           READ(5,37) AP1(I),BP1(I),CP1(I),MADP(I)
           READ(5,40) AP2(I),BP2(I),CP2(I),AP3(I),BP3(I),CP3(I)
           READ(5,37) AD1(I),BD1(I),CD1(I),MADD(I)
           READ(5,40) AD2(I),BD2(I),CD2(I),AD3(I),BD3(I),CD3(I)
36         CONTINUE
C
C     READ IN IOVPOP(I) AND IENRGY(I). INDIVIDULAL OVERLAP POPULATION
C     ANALYSES ARE PERFORMED FROM ORBITAL IOVPOP(N) TO ORBITAL
C     IOVPOP(N+1). INDIVIDUAL ENERGY MATRIX ANALYSES ARE PERFORMED
C     FROM ORBITAL IENRGY(N) TO ORBITAL IENRGY(N+1).
C
60         IF(L3) READ(5,67) IOVPOP
67         FORMAT(24I3)
           IF(L4) READ(5,67) IENRGY
C
C     PRINT OUT TYPE OF CALCULATION.
C
           IF(METH.EQ.0) WRITE(6,90)
           IF(METH.EQ.1) WRITE(6,91)
           IF(METH.GE.2) WRITE(6,92)
           IF(METH.GT.2) WRITE(6,93)
           IF(L5) WRITE(6,94)
90         FORMAT(///,T10,'EXTENDED HUCKEL CALCULATION.')
91         FORMAT(///,T10,'EXTENDED HUCKEL CALCULATION WITH CHARGE',
     .     ' ITERATION.',/,T10,'LINEAR CHARGE DEPENDENCE OF SENSE*CHARG'
     .     ,'E FOR H(I,I)''S.')
92         FORMAT(///,T10,'EXTENDED HUCKEL CALCULATION WITH CHARGE',
     .     ' ITERATION.')
93         FORMAT(T10,'MADELUNG CORRECTION INCLUDED.')
94         FORMAT(T10,'WEIGHTED HIJ FORMULA USED.')
C
C     PRINT OUT ATOMIC COORDINATES AND PARAMETERS.
C
           IF(PRT(1)) GO TO 80
           WRITE(6,74)
74         FORMAT(///,T5,'ATOM',T17,'X',T29,'Y',T41,'Z',T56,'S',T76,'P',
     .     T96,'D',T113,'CONTRACTED D'/T47,'N',T50,'EXP',T59,'COUL',
     .     T67,'N',T70,'EXP',T79,'COUL',T87,'N',T90,'EXPD1',T99,'COUL',
     .     T109,'C1',T118,'C2',T125,'EXPD2')
           IF(NH.EQ.0) GO TO 72
           J=1
           DO 76 I=1,NH
76         WRITE(6,53) HYDROG,I,X(I),Y(I),Z(I),J,PEEP,COULH
53         FORMAT(T4,A2,I3,3F12.5,3(I3,F8.4,F9.4),2F9.5,F8.4)
72         CONTINUE
           DO 151 I=1,NA
           KEYI=AC(I)
           INH=I+NH
           IF(NP(KEYI).NE.0) GO TO 152
*           WRITE(6,53) SYMBOL(KEYI),INH,X(INH),Y(INH),Z(INH),NS(KEYI),
*     .     EXPS(KEYI),COULS(KEYI)
           WRITE(6,53) SYMBH(I),INH,X(INH),Y(INH),Z(INH),NS(KEYI),
     .     EXPS(KEYI),COULS(KEYI)
           GO TO 151
152        IF(ND(KEYI).NE.0) GO TO 153
*           WRITE(6,53) SYMBOL(KEYI),INH,X(INH),Y(INH),Z(INH),NS(KEYI),
*     .     EXPS(KEYI),COULS(KEYI),NP(KEYI),EXPP(KEYI),COULP(KEYI)
           WRITE(6,53) SYMBH(I),INH,X(INH),Y(INH),Z(INH),NS(KEYI),
     .     EXPS(KEYI),COULS(KEYI),NP(KEYI),EXPP(KEYI),COULP(KEYI)
           GO TO 151
153        IF(C2(KEYI).NE.0.0D0) GO TO 154
*           WRITE(6,53) SYMBOL(KEYI),INH,X(INH),Y(INH),Z(INH),NS(KEYI),
*     .     EXPS(KEYI),COULS(KEYI),NP(KEYI),EXPP(KEYI),COULP(KEYI),
*     .     ND(KEYI),EXPD(KEYI),COULD(KEYI)
           WRITE(6,53) SYMBH(I),INH,X(INH),Y(INH),Z(INH),NS(KEYI),
     .     EXPS(KEYI),COULS(KEYI),NP(KEYI),EXPP(KEYI),COULP(KEYI),
     .     ND(KEYI),EXPD(KEYI),COULD(KEYI)
           GO TO 151
*154        WRITE(6,53) SYMBOL(KEYI),INH,X(INH),Y(INH),Z(INH),NS(KEYI),
*     .     EXPS(KEYI),COULS(KEYI),NP(KEYI),EXPP(KEYI),COULP(KEYI),
*     .     ND(KEYI),EXPD(KEYI),COULD(KEYI),C1(KEYI),C2(KEYI),EXPD2(KEYI)
154        WRITE(6,53) SYMBH(I),INH,X(INH),Y(INH),Z(INH),NS(KEYI),
     .     EXPS(KEYI),COULS(KEYI),NP(KEYI),EXPP(KEYI),COULP(KEYI),
     .     ND(KEYI),EXPD(KEYI),COULD(KEYI),C1(KEYI),C2(KEYI),EXPD2(KEYI)
151        CONTINUE
           WRITE(6,160) KA,IPRINT,IPUNCH,CON
160        FORMAT(///,T10,'CHARGE =',I3,8X,'IPRINT =',I3,8X,'IPUNCH =',
     .     I3,8X,'HUCKEL CONSTANT =',F7.3)
C
C     PRINT OUT ITERATION PARAMETERS.
C
           IF(.NOT.ITERAT) GO TO 80
           WRITE(6,81) DAMP1,DAMP2,DAMP3,LAMPRI,MAXCYC,PRTCYC,
     .     SENSE,DELTAC    
81         FORMAT(/,T10,'DAMP1 =',F6.3,6X,'DAMP2 =',F6.3,6X,'DAMP3 =',
     .     F6.3,6X,'LAMPRI =',F6.3,//,T10,'MAXCYC =',I3,8X,'PRTCYC =',
     .     I3,8X,'SENSE =',F6.3,6X,'DELTAC =',F10.7)
C
C     PRINT OUT VSIE PARAMETERS.
C
           IF(METH.LT.2) GO TO 80
           WRITE(6,82)
82         FORMAT(///,' VSIE PARAMETERS',//,T10,'ATOM',T26,'A',T39,'B',
     .     T52,'C')
           NUSER2=21-NUSER
           DO 83 I=1,NUSER2
           J=21-I
           IF(.NOT.ITABLE(J)) GO TO 83
           WRITE(6,84) SYMBOL(J),AS1(I),BS1(I),CS1(I)
84         FORMAT(/,T11,A2,4X,3F13.5)
           IF(NP(J).EQ.0) GO TO 83
           IF(ND(J).NE.0) GO TO 85
           WRITE(6,86) AP1(I),BP1(I),CP1(I)
86         FORMAT(T17,3F13.5)
           GO TO 83
85         IF(NCON.EQ.3) GO TO 87
           WRITE(6,86) AP1(I),BP1(I),CP1(I),AD1(I),BD1(I),CD1(I)
           GO TO 83
87         WRITE(6,86) AS2(I),BS2(I),CS2(I),AS3(I),BS3(I),CS3(I),AP1(I),
     .     BP1(I),CP1(I),AP2(I),BP2(I),CP2(I),AP3(I),BP3(I),CP3(I),
     .     AD1(I),BD1(I),CD1(I),AD2(I),BD2(I),CD2(I),AD3(I),BD3(I),
     .     CD3(I)
83         CONTINUE
80         WRITE(6,99)
99         FORMAT(///)
           RETURN
           END
           SUBROUTINE MOVLAP(H,S,MAD,C,SP,PD,MAXS,MAXP,MAXD,COUL0,SORB,
     .     IOCC,HDG,     NDIM, ND1 ,NC,NATOM,NTYPE,NHDG)
C
C     SUBROUTINE TO CALCULATE INTERATOMIC DISTANCES, OVERLAP
C     INTEGRALS, AND MADELUNG PARAMETERS.
C
           IMPLICIT DOUBLE PRECISION (A-H,O-Z)
           REAL*8 MAD
           LOGICAL*4 SP,PD
           INTEGER*4 COUL0
           INTEGER*4 SORB
           INTEGER*4 STAR,STAR2
           REAL*4 MADS(20),MADP(20),MADD(20)
           LOGICAL*4 JGO
           DIMENSION H(NDIM,NDIM),S(NDIM,NDIM),MAD(NTYPE,NTYPE),C(NC),
     .     SP(NDIM),PD(NDIM),MAXS(NDIM),MAXP(NDIM),MAXD(NDIM),COUL0(20),
     .     SORB(NATOM),IOCC(NDIM),HDG(NHDG)
C
C     COUL0 DIMENSIONED AT 20 FOR EASY READING DURING PROCESSING
C     OF DELETION INPUT. DELETIONS DONE IN SUBROUTINE DELETS.
C
           COMMON/TITLE/AB(10)
           COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,
     .     IPRINT,IPUNCH,LA,LB,L3,L4,L5,ONEMAT,ITERAT
           LOGICAL*4 LA,LB,L3,L4,L5,ONEMAT,ITERAT
           COMMON/OUT/PRT(20),PUN(20),IOVPOP(24),IENRGY(24)
           LOGICAL*4 PRT,PUN
           COMMON/ATOM/KEY(1000),SYMBOL(100),VELEC(40),NS(40),NP(40),
     .     ND(40),EXPS(40),EXPP(40),EXPD(40),EXPD2(40),C1(40),C2(40),
     .     COULS(40),COULP(40),COULD(40),X(1000),Y(1000),Z(1000)
           INTEGER*4 SYMBOL,KEY,VELEC
           COMMON/STARS/STAR,STAR2
           COMMON /LOCLAP/ SK1,SK2,R,L1,L2,M,N1,N2,MAX
           EQUIVALENCE (MADS(1),NS(21)),(MADP(1),NP(21)),
     .     (MADD(1),ND(21))
           DIMENSION PTR(9),DTR(25)
           DIMENSION A(20),B(20),A1(20),B1(20)
           EQUIVALENCE (PTR(3),CA),(PTR(8),CB)
           DATA SQRT3/1.7320508075688770/
           DATA AUI/1.889644746D0/
           DATA PTR(9)/0.D0/,DTR(12)/0.D0/,DTR(22)/0.D0/
           NH1=NH+1
C
C     HYDROGEN-HYDROGEN OVERLAPS.
C
           IF(NH.LE.1) GO TO 106
           DO 107 I=2,NH
           IM1=I-1
           DO 107 J=1,IM1
           R=SQRT((X(I)-X(J))**2+(Y(I)-Y(J))**2+(Z(I)-Z(J))**2)
C
C     STORE OVERLAPS IN UPPER RIGHT TRIANGLE OF S(I,J). PUT
C     DISTANCES IN LOWER RIGHT TRIANGLE.
C
           S(I,J)=R
           R=R*PEEP*AUI
           IF(R.GT.50) GO TO 105
           SIGMA=(1.D0+R*(1.D0+R/3.D0))/EXP(R)
           GO TO 107
105        SIGMA=0.D0
107        S(J,I)=SIGMA
C
C     HEAVY ATOM-HEAVY ATOM OVERLAPS. LOCAL COORDINATE SYSTEM
C     CENTERED ON ATOM J. FILL IN UPPER RIGHT TRIANGLE OF S(I,J).
C
106        IORB=NH1
           DO 130 I=1,NA
           INH=I+NH
           IM1=I-1
           KEYI=KEY(I)
           MAXD(I)=ND(KEYI)
           MAXP(I)=NP(KEYI)
           MAXS(I)=NS(KEYI)
           SP(I)=EXPS(KEYI) .EQ. EXPP(KEYI)
           PD(I)=EXPP(KEYI) .EQ. EXPD(KEYI)
           IF(PD(I)) MAXP(I)=MAX0(MAXP(I),MAXD(I))
           IF(SP(I)) MAXS(I)=MAX0(NS(KEYI),MAXP(I))
           SORB(I)=IORB
C
C     SORB(I) CONTAINS A POINTER TO THE S ORBITAL ON ATOM I.
C
           IORBS=IORB
           IORB=IORB+4
           IF(NP(KEYI).EQ.0) IORB=IORB-3
           IF(ND(KEYI).NE.0) IORB=IORB+5
           IF(NP(KEYI).EQ.0) GO TO 298
           JD=IORB-1
           JD1=JD-1
           DO 280 JC = IORBS,JD1
           ID=JC+1
           DO 280 IC=ID,JD
280        S(JC,IC)=0.D0
298        CONTINUE
           IF(I.EQ.1) GO TO 300
           DO 131 J=1,IM1
           KEYJ=KEY(J)
           JNH=J+NH
           JORBS=SORB(J)
           DELX=X(INH)-X(JNH)
           DELY=Y(INH)-Y(JNH)
           DELZ=Z(INH)-Z(JNH)
           RT2=DELX**2+DELY**2
           R=SQRT(RT2+DELZ**2)
           S(INH,JNH)=R
C
C     STORE DISTANCES IN LOWER LEFT TRIANGLE OF S(I,J).
C
           IF(R.GT.0.0D0) GO TO 102
           ID=IORB-1
           JD=SORB(J+1)-1
           DO 103 IC=IORBS,ID
           DO 103 JC=JORBS,JD
103        S(JC,IC)=0.0D0
           GO TO 131
102        IF(RT2.GT.1.E-10) GO TO 135
           CB=1.D0
           SB=0.D0
           SA=0.D0
           GOTO 136
135        T=SQRT(RT2)
           CB=DELX/T
           SB=DELY/T
           SA=T/R
136        CA=DELZ/R
C
C     THE TRANSFORMATION MATRICES ARE CALCULATED EXPLICITLY.
C     PTR IS THE MATRIX FOR PROJECTING THE X,Y,Z ORBITALS
C     ONTO THE LOCAL SYSTEM. THE ELEMENTS ARE ORDERED SO THAT FIRST
C     X THEN Y THEN Z IS PROJECTED ONTO THE Z' AXIS (SIGMA).
C     THEN THE 3 ARE PROJECTED ONTO THE X' AXIS AND THEN THE Y' (PI).
C     THE D ORBITALS ARE HANDLED SIMILARLY. THE ORDER OF PROJECTION
C     IS X2-Y2,Z2,XY,XZ,YZ FIRST ONTO Z2'(SIGMA)AND THEN ONTO XZ' AND
C     YZ'(PI). FINALLY THE 5 ORBITALS ARE PROJECTED ONTO X'2-Y'2 AND
C     THEN XY' (DELTA).
C
C     THOSE PTR AND DTR WHICH ARE ZERO ARE INITIALIZED IN A DATA STATE-
C     MENT.  CA AND CB HAVE BEEN EQUIVALENCED TO PTR(3) AND PTR(8)
C     RESPECTIVELY TO SAVE TIME.
C
           PTR(1)= SA*CB
           PTR(2)= SA*SB
C ...      PTR(3)= CA
           PTR(4)= CA*CB
           PTR(5)= CA*SB
           PTR(6)= -SA
           PTR(7)= -SB
C ...      PTR(8)= CB
C ...      PTR(9)= 0.D0
           IF(ND(KEYI)+ND(KEYJ).EQ.0) GO TO 180
           CA2=CA**2
           SA2=SA*SA
           CB2=CB*CB
           SB2=SB*SB
           CBSB= CB*SB
           CASA= CA*SA
           CB2SB2= CB2-SB2
           DTR(1)= SQRT3*.5D0*SA2*CB2SB2
           DTR(2)= 1.D0-1.5D0*SA2
           DTR(3)= SQRT3*CBSB*SA2
           DTR(4)= SQRT3*CASA*CB
           DTR(5)= SQRT3*CASA*SB
           DTR(6)= CASA*CB2SB2
           DTR(7)= -SQRT3*CASA
           DTR(8)= 2.D0*CASA*CBSB
           DTR(9)= CB*(CA2-SA2)
           DTR(10)= SB*(CA2-SA2)
           DTR(11)= -2.D0*SA*CBSB
C ...      DTR(12)= 0.D0
           DTR(13)= SA* CB2SB2
           DTR(14)= -PTR(5)
           DTR(15)= PTR(4)
           IF(ND(KEYI)*ND(KEYJ).EQ.0) GO TO 180
           DTR(16)=.5D0*(1.D0+CA2)*CB2SB2
           DTR(17)= .5D0*SQRT3*SA2
           DTR(18)= CBSB*(1.D0+CA2)
           DTR(19)= -CASA*CB
           DTR(20)= -CASA*SB
           DTR(21)= -2.D0*CA*CBSB
C ...      DTR(22)= 0.D0
           DTR(23)= CA*CB2SB2
           DTR(24)= PTR(2)
           DTR(25)= -PTR(1)
180        R=R*AUI
C
C     (S(I) S(J)).
C
           N2=NS(KEYJ)
           N1=NS(KEYI)
           L2=0
           L1=0
           M=0
           MAX=MAXS(I)+MAXS(J)
           SK1=EXPS(KEYI)
           SK2=EXPS(KEYJ)
           CALL ABFNS(A,B)
           CALL LOVLAP(SIGMA,A,B)
           S(JORBS,IORBS)=SIGMA
C
C     IF THE S EXPONENT OF ATOM I EQUALS THE P EXPONENT WE NEED
C     NOT CALCULATE THE A AND B FUNCTIONS AGAIN.
C
C     (P(I) S(J)).
C
           JGO=.FALSE.
           IF(KEYI.EQ.KEYJ) GO TO 126
           IF((.NOT.SP(I)).OR.(NP(KEYI).EQ.0)) GO TO 126
220        N1=NP(KEYI)
           L1=1
           CALL LOVLAP(SIGMA,A,B)
           SIGMA=-SIGMA
           DO 200 IC=1,3
200        S(JORBS,IORBS+IC)=PTR(IC)*SIGMA
           IF(PD(I)) GO TO 221
           IF(JGO) GO TO 217
           GO TO 137
C
C     (D(I) S(J)) CONDITIONALLY AT FIRST CHANCE.
C
221        N1=ND(KEYI)
           L1=2
168        CALL LOVLAP(SIGMA,A,B)
           IF(C2(KEYI).EQ.0.D0) GO TO 167
           SK1=EXPD2(KEYI)
           CALL ABFNS(A1,B1)
           CALL LOVLAP(PART2,A1,B1)
           SIGMA=C1(KEYI)*SIGMA+C2(KEYI)*PART2
           SK1=EXPD(KEYI)
167        ID=IORBS+3
           DO 201 IC=1,5
201        S(JORBS,ID+IC)=DTR(IC)*SIGMA
C
C     CALCULATE (D(I) P(J)) IF CAN USE SAME A'S AND B'S.
C
           IF(SP(J)) GO TO 222
           IF(JGO) GO TO 228
           GO TO 137
222        N2=NP(KEYJ)
           L2=1
           M=0
           CALL LOVLAP(SIGMA,A,B)
           M=1
           CALL LOVLAP(PI,A,B)
           IF(C2(KEYI).EQ.0.D0) GO TO 1169
           SK1=EXPD2(KEYI)
           CALL LOVLAP(PART2,A1,B1)
           PI=C1(KEYI)*PI+C2(KEYI)*PART2
           M=0
           CALL LOVLAP(PART2,A1,B1)
           SK1=EXPD(KEYI)
           SIGMA=C1(KEYI)*SIGMA+C2(KEYI)*PART2
1169       PI=-PI
           ID=IORBS+3
           DO 195 JC=1,3
           DO 195 IC=1,5
195        S(JORBS+JC,ID+IC)=PTR(JC)*DTR(IC)*SIGMA+(PTR(JC+3)*DTR(IC+5)
     .     +PTR(JC+6)*DTR(IC+10))*PI
           IF(JGO) GO TO 131
C
C     NOW TEST FOR DUPLICATE EXPONENTS ON ATOM J.
C     HOWEVER DO CALCULATIONS ANYHOW.
C
137        N1=NS(KEYI)
           L1=0
C
C     (S(I) P(J)).
C
126        IF(SP(J)) GO TO 138
           IF(NP(KEYJ).EQ.0) GO TO 210
           MAX=MAXS(I)+MAXP(J)
           SK2=EXPP(KEYJ)
           CALL ABFNS(A,B)
138        N2=NP(KEYJ)
           L2=1
           M=0
           CALL LOVLAP(SIGMA,A,B)
           DO 202 IC=1,3
202        S(JORBS+IC,IORBS)=PTR(IC)*SIGMA
           IF(SP(I)) GO TO 156
           JGO=.TRUE.
           IF(ND(KEYJ).NE.0) GO TO 149
C
C     BRANCH TO TEST FOR EXPP(J) .EQ. EXPD(J). CALCULATE (S D) ANYHOW.
C     RETURN WILL BE MADE TO THE NEXT STATEMENT.
C
C     (P(I) P(J))   EXPP(I) EQ,NE EXPS(I).
C
           GO TO 646
146        N2=NP(KEYJ)
           L2=1
           SK2=EXPP(KEYJ)
646        IF(NP(KEYI).EQ.0) GO TO 210
           SK1=EXPP(KEYI)
C
C     THESE STATEMENTS USED ONLY IF HAVE ALREADY CALCULATED (S(I) D(J))
C     WHICH MEANS THAT SP(I) IS FALSE.
C
           MAX=MAXP(I)+MAXP(J)
           CALL ABFNS(A,B)
156        N1=NP(KEYI)
           L1=1
148        M=0
           CALL LOVLAP(SIGMA,A,B)
           SIGMA=-SIGMA
           M=1
           CALL LOVLAP(PI,A,B)
           DO 204 JC=1,3
           DO 204 IC=JC,3
           S(JORBS+JC,IORBS+IC)=PTR(JC)*PTR(IC)*SIGMA + (PTR(JC+3)*
     .     PTR(IC+3)+PTR(JC+6)*PTR(IC+6))*PI
204        S(JORBS+IC,IORBS+JC)=S(JORBS+JC,IORBS+IC)
147        IF(ND(KEYJ).EQ.0) GO TO 210
C
C     BRANCH AROUND (S(I) D(J)) IF ALREADY DONE.
C
           IF(JGO) GO TO 160
C
C     (S(I) D(J)).
C
           N1=NS(KEYI)
           L1=0
149        N2=ND(KEYJ)
           L2=2
           IF(PD(J)) GO TO 142
           SK2=EXPD(KEYJ)
           MAX=MAXS(I)+MAXD(J)
           CALL ABFNS(A,B)
142        M=0
           CALL LOVLAP(SIGMA,A,B)
           IF(C2(KEYJ).EQ.0.D0) GO TO 151
           SK2=EXPD2(KEYJ)
           CALL ABFNS(A1,B1)
           CALL LOVLAP(PART2,A1,B1)
           SIGMA=C1(KEYJ)*SIGMA+C2(KEYJ)*PART2
           SK2=EXPD(KEYJ)
151        JD=JORBS+3
           DO 205 IC=1,5
205        S(JD+IC,IORBS)=DTR(IC)*SIGMA
150        IF(JGO) GO TO 146
C
C     SP(I) IS TRUE IF HERE SO BRANCH AS WE ALSO HAVE D ON ATOM J.
C
           GO TO 170
160        JGO=.FALSE.
C
C          (P(I) D(J)).
C
           N2=ND(KEYJ)
           L2=2
           IF(PD(J)) GO TO 178
           SK2=EXPD(KEYJ)
           MAX=MAXP(I)+MAXD(J)
           CALL ABFNS(A,B)
178        IF(C2(KEYJ).EQ.0.D0) GO TO 170
           SK2=EXPD2(KEYJ)
           CALL ABFNS(A1,B1)
           SK2=EXPD(KEYJ)
170        N1=NP(KEYI)
           L1=1
           M=0
           CALL LOVLAP(SIGMA,A,B)
           M=1
           CALL LOVLAP(PI,A,B)
           IF(C2(KEYJ).EQ.0.D0) GO TO 171
           SK2=EXPD2(KEYJ)
           CALL LOVLAP(PART2,A1,B1)
           PI=C1(KEYJ)*PI+C2(KEYJ)*PART2
           M=0
           CALL LOVLAP(PART2,A1,B1)
           SIGMA=C1(KEYJ)*SIGMA+C2(KEYJ)*PART2
           SK2=EXPD(KEYJ)
171        SIGMA=-SIGMA
           DO 206 IC=1,3
           DO 206 JC=1,5
206        S(JD+JC,IORBS+IC)=DTR(JC)*PTR(IC)*SIGMA+(DTR(JC+5)*PTR(IC+3)
     .     +DTR(JC+10)*PTR(IC+6))*PI
C
C     (D(I) D(J)).
C
           IF(ND(KEYI).EQ.0) GO TO 210
           MAX=MAXD(I)+MAXD(J)
           IF(PD(I)) GO TO 208
           SK1=EXPD(KEYI)
           CALL ABFNS(A,B)
           IF(C2(KEYJ).EQ.0.D0) GO TO 208
           SK2=EXPD2(KEYJ)
           CALL ABFNS(A1,B1)
           SK2=EXPD(KEYJ)
208        N1=ND(KEYI)
           L1=2
           M=0
           CALL LOVLAP(SIGMA,A,B)
           M=1
           CALL LOVLAP(PI,A,B)
           M=2
           CALL LOVLAP(DELTA,A,B)
           CC=C2(KEYI)
           IF(C2(KEYJ).EQ.0.D0) GO TO 173
           CC=C1(KEYJ)*CC
           SK2=EXPD2(KEYJ)
           CALL LOVLAP(PART2,A1,B1)
           DELTA=C1(KEYJ)*DELTA+C2(KEYJ)*PART2
           M=1
           CALL LOVLAP(PART3,A1,B1)
           PI=C1(KEYJ)*PI+C2(KEYJ)*PART3
           M=0
           CALL LOVLAP(PART4,A1,B1)
           SIGMA=C1(KEYJ)*SIGMA+C2(KEYJ)*PART4
           SK2=EXPD(KEYJ)
           M=2
173        IF(C2(KEYI).EQ.0.D0) GO TO 172
           IF(KEYI.EQ.KEYJ) GO TO 176
           SK1=EXPD2(KEYI)
           CALL ABFNS(A1,B1)
           CALL LOVLAP(PART2,A1,B1)
           M=1
           CALL LOVLAP(PART3,A1,B1)
           M=0
           CALL LOVLAP(PART4,A1,B1)
176        SIGMA=C1(KEYI)*SIGMA+CC*PART4
           PI =C1(KEYI)*PI+CC*PART3
           DELTA=C1(KEYI)*DELTA+CC*PART2
           IF(C2(KEYJ).EQ.0.D0) GO TO 172
           SK1=EXPD2(KEYI)
           SK2=EXPD2(KEYJ)
           CALL ABFNS(A1,B1)
           M=0
           CALL LOVLAP(PART2,A1,B1)
           CC=C2(KEYI)*C2(KEYJ)
           SIGMA=SIGMA+CC*PART2
           M=1
           CALL LOVLAP(PART2,A1,B1)
           PI=PI+CC*PART2
           M=2
           CALL LOVLAP(PART2,A1,B1)
           DELTA=DELTA+CC*PART2
172        PI=-PI
           JD=JORBS+3
           DO 211 IC=1,5
           ID=IORBS+3
           DO 211 JC=1,5
           S(JD+JC,ID+IC) = DTR(IC)*DTR(JC)*SIGMA+(DTR(IC+5)*DTR(JC+5)
     .     +DTR(IC+10)*DTR(JC+10))*PI+(DTR(IC+15)*DTR(JC+15)+DTR(IC+20)
     .     *DTR(JC+20))*DELTA
211        S(JD+IC,ID+JC)=S(JD+JC,ID+IC)
C
C     FILLING IN OTHER HALF OF OVERLAPS FOR (J I) AS NEEDED.
C
210        IF(KEYI.EQ.KEYJ) GO TO 213
           N2=NS(KEYJ)
           L2=0
           SK2=EXPS(KEYJ)
           M=0
           JGO=.TRUE.
           IF(NP(KEYI).EQ.0) GO TO 131
           IF(SP(I)) GO TO 215
           MAX=MAXP(I)+MAXS(J)
           SK1=EXPP(KEYI)
           CALL ABFNS(A,B)
           GO TO 220
215        IF(PD(I)) GO TO 227
217        IF(ND(KEYI).EQ.0) GO TO 131
           MAX=MAXD(I)+MAXS(J)
           SK1=EXPD(KEYI)
           CALL ABFNS(A,B)
           GO TO 221
227        IF(SP(J)) GO TO 131
           N1=ND(KEYI)
           L1=2
           SK1=EXPD(KEYI)
228        IF(NP(KEYJ).EQ.0) GO TO 131
           SK2=EXPP(KEYJ)
           MAX=MAXD(I)+MAXP(J)
           CALL ABFNS(A,B)
           IF(C2(KEYI).EQ.0.D0) GO TO 222
           SK1=EXPD2(KEYI)
           CALL ABFNS(A1,B1)
           SK1=EXPD(KEYI)
           GO TO 222
213        IF(NP(KEYI).EQ.0) GO TO 131
           DO 237 IC=1,3
237        S(JORBS,IORBS+IC)=-S(JORBS+IC,IORBS)
           IF(ND(KEYI).EQ.0) GO TO 131
           DO 238 IC=4,8
           S(JORBS,IORBS+IC)=S(JORBS+IC,IORBS)
           DO 238 JC=1,3
238        S(JORBS+JC,IORBS+IC)=-S(JORBS+IC,IORBS+JC)
131        CONTINUE
300        IF(NH.EQ.0) GO TO 130
           N2=1
           L2=0
           M=0
           SK2=PEEP
           DO 301 J=1,NH
           DELX=X(J)-X(INH)
           DELY=Y(J)-Y(INH)
           DELZ=Z(J)-Z(INH)
           RT2=DELX**2+DELY**2
           R=SQRT(RT2+DELZ**2)
C
C     STORE DISTANCES IN LOWER LEFT TRIANGLE OF S(I,J).
C
           S(INH,J)=R
           IF(RT2.GT.1.D-10) GO TO 303
           CB=1.D0
           SB=0.D0
           SA=0.D0
           GO TO 302
303        T=SQRT(RT2)
           CB=DELX/T
           SB=DELY/T
           SA=T/R
302        CA=DELZ/R
           R=R*AUI
C
C     H(J) S(I)).
C
           N1=NS(KEYI)
           L1=0
           MAX=1+MAXS(I)
           SK1=EXPS(KEYI)
           CALL ABFNS(A,B)
           CALL LOVLAP(SIGMA,A,B)
           S(J,IORBS)=SIGMA
           IF(NP(KEYI).EQ.0) GO TO 301
           IF(SP(I)) GO TO 304
           SK1=EXPP(KEYI)
           MAX=1+MAXP(I)
           CALL ABFNS(A,B)
304        N1=NP(KEYI)
           L1=1
           CALL LOVLAP(SIGMA,A,B)
           S(J,IORBS+3)=CA*SIGMA
           SIGMA=SIGMA*SA
           S(J,IORBS+2)=SB*SIGMA
           S(J,IORBS+1)=CB*SIGMA
           IF(ND(KEYI).EQ.0) GO TO 301
           IF(PD(I)) GO TO 305
           SK1=EXPD(KEYI)
           MAX=1+ND(KEYI)
           CALL ABFNS(A,B)
305        N1=ND(KEYI)
           L1=2
           CALL LOVLAP(SIGMA,A,B)
           IF(C2(KEYI).EQ.0.D0) GO TO 181
           SK1=EXPD2(KEYI)
           CALL ABFNS(A1,B1)
           CALL LOVLAP(PART2,A1,B1)
           SK1=EXPD(KEYI)
           SIGMA=C1(KEYI)*SIGMA+C2(KEYI)*PART2
181        CONTINUE
           S(J,IORBS+5)=(1.D0-1.5D0*SA*SA)*SIGMA
           SIGMA=SIGMA*SQRT3*SA
           S(J,IORBS+4)=.5D0*SA*(CB*CB-SB*SB)*SIGMA
           S(J,IORBS+6)=CB*SB*SA*SIGMA
           SIGMA=SIGMA*CA
           S(J,IORBS+7)=CB*SIGMA
           S(J,IORBS+8)=SB*SIGMA
301        CONTINUE
130        CONTINUE
C
C     CALL DELETS TO SET CERTAIN OVERLAP INTEGRALS = 0.
C
           IF(.NOT.LB) GO TO 835
           CALL DELETS(S,COUL0,SORB,NDIM)
           WRITE(6,2010)
2010       FORMAT(///)
C
C     CALCULATE INTERATOMIC MADELUNG PARAMETERS.
C
835        IF(METH.LT.3) GO TO 450
           IC=1
           DO 401 I=1,NA
           KEYI=KEY(I)
           RS=0.0D0
           ID=IC      
           MAD(ID,ID)=REAL(MADS(21-KEYI))
           IF(NP(KEYI).EQ.0) GO TO 402
           ID=IC+1
           MAD(ID,ID)=REAL(MADP(21-KEYI))
           IF(ND(KEYI).EQ.0) GO TO 403
           ID=IC+2
           MAD(ID,ID)=REAL(MADD(21-KEYI))
403        M=IC+1
           DO 404 K=M,ID
           K1=K-1
           DO 404 L=IC,K1
           CA=MAD(K,K)
           CB=MAD(L,L)
           SA=VALMAD(CA,CB,RS)
           MAD(K,L)=SA
404        MAD(L,K)=SA
402        IF(I.EQ.1) GO TO 401
           IM1=I-1
           JC=1
           DO 406 J=1,IM1
           KEYJ=KEY(J)
           RS=S(I,J)*AUI/27.21D0
           JD=JC
           IF(NP(KEYJ).NE.0) JD=JC+1
           IF(ND(KEYJ).NE.0) JD=JC+2
           DO 407 K=IC,ID
           CA=MAD(K,K)
           DO 407 L=JC,JD
           CB=MAD(L,L)
           SA=VALMAD(CA,CB,RS)
           MAD(K,L)=SA
407        MAD(L,K)=SA
406        JC=JD+1
401        IC=ID+1
C
C     SET UP DISTANCE MATRIX FOR PRINTING.
C     STUFF ELEMENTS OF S INTO C TO GET THEM OUT OF THE WAY.
C
450        ISUB=1
C
C     ZERO DISTANCE ALONG DIAGONAL.
C
           S(1,1)=0.D0
           DO 1010 I=2,NATOM
           S(I,I)=0.D0
           IM1=I-1
           DO 1005 J=1,IM1
           C(ISUB)=S(J,I)
           ISUB=ISUB+1
1005       S(J,I)=S(I,J)
1010       CONTINUE
           IF(PRT(3)) GO TO 2004
           WRITE(6,2000)
2000       FORMAT('DISTANCE MATRIX')
           CALL PEGLEG(S,NATOM,NDIM)
2004       IF(PUN(3)) WRITE(7,2050) ((S(I,J),I=1,NATOM),J=1,NATOM)
2050       FORMAT(8F9.6)
C
C     SET UP OVERLAP MATRIX FOR PRINTING.
C     REPLACE ELEMENTS IN OVERLAP MATRIX FROM C.
C
1015       S(1,1)=1.D0
           ISUB=1
           DO 1025 I=2,NDIM
           S(I,I)=1.D0
           IM1=I-1
           DO 1020 J=1,IM1
           IF(I.GT.NATOM) GO TO 1020
           S(J,I)=C(ISUB)
           ISUB=ISUB+1
1020       S(I,J)=S(J,I)
1025       CONTINUE
           IF(PRT(4)) GO TO 2005
           WRITE(6,2001)
2001       FORMAT('OVERLAP MATRIX')
           CALL PEGLEG(S,NDIM,NDIM)
2005       IF(PUN(4)) WRITE(7,2050) S
C
C     PRINT OUT MADELUNG PARAMETERS.
C
           IF(METH.LT.3) GO TO 460
           IF(PRT(5)) GO TO 2006
           WRITE(6,2002)
2002       FORMAT('MADELUNG PARAMETERS')
           CALL PEGLEG(MAD,NTYPE,NTYPE)
2006       IF(PUN(5)) WRITE(7,2050) MAD
460        RETURN
           END
      SUBROUTINE DELETS(S,COUL0,SORB,NDIM)
C
C  SUBROUTINE FOR SETTING CERTAIN OVERLAP INTEGRALS EQUAL TO ZERO.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S(NDIM,NDIM),COUL0(20),SORB(NDIM)
      INTEGER*4 COUL0
      INTEGER*4 SORB
      LOGICAL*4 L1,L2,L3,L4,L5,ONEMAT,ITERAT
      LOGICAL*4 PRT,PUN
      LOGICAL*4 IERR
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,
     1 IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT
      COMMON/OUT/PRT(20),PUN(20),IOVPOP(24),IENRGY(24)
      DATA ORBTL/' ORBITAL'/,ATMPR/'  ATOM  '/
      SORB(NA+1)=NDIM+1
      IERR=.FALSE.
C
C  READ IN NUMBERS INDICATING WHICH OVERLAP INTEGRALS ARE TO BE SET
C  TO ZERO. A POSITIVE NUMBER REFERS TO AN ATOM, A NEGATIVE ONE TO
C  AN ORBITAL.
C
10    READ(5,1000,END=300) COUL0
1000  FORMAT(20I3)
      DO 200 IDEL=1,19,2
      I=COUL0(IDEL)
      J=COUL0(IDEL+1)
C
C  TERMINATE ON ENCOUNTERING A ZERO.
C
      IF(I.EQ.0.OR.J.EQ.0) GO TO 400
      IF(IERR) GO TO 200
      IABSV=IABS(I)
      JABSV=IABS(J)
      IF(I.GT.NH) GO TO 20
C
C  I REFERS TO ORBITAL (NEGATIVE) OR H ATOM (POSITIVE,LE NH).
C
      ILOW=IABSV
      IHIGH=IABSV
C
C  ERROR IF I OUT OF RANGE OF ORBITAL NUMBERS.
C
      IF(IABSV.GT.NDIM) GO TO 160
      GO TO 30
C
C  ERROR IF I OUT OF RANGE OF ATOMS.
C
20    IF(I.GT.NATM) GO TO 160
      ILOW=SORB(IABSV-NH)
      IHIGH=SORB(IABSV-NH+1)-1
30    IF(J.GT.NH) GO TO 40
      JLOW=JABSV
      JHIGH=JABSV
C
C  CHECK TO SEE IF J IS IN RANGE.
C
      IF(JABSV.GT.NDIM) GO TO 160
      GO TO 50
40    IF(J.GT.NATM) GO TO 160
      JLOW=SORB(JABSV-NH)
      JHIGH=SORB(JABSV-NH+1)-1
50    X1=ATMPR
      IF(I.LT.0) X1=ORBTL
      X2=ATMPR
      IF(J.LT.0) X2=ORBTL
      IF(.NOT.PRT(2)) WRITE(6,1002) X1,IABSV,X2,JABSV
1002  FORMAT('ALL S(I,J) SET TO ZERO BETWEEN',A8,I4,' AND',A8,I4,'.')
C
C  J MUST BE LESS THAN OR EQUAL TO I SINCE WE ONLY HAVE A HALF
C  MATRIX AT THIS POINT.
C
      IF(JHIGH.LE.IHIGH) GO TO 60
      I=JHIGH
      JHIGH=IHIGH
      IHIGH=I
      I=JLOW
      JLOW=ILOW
      ILOW=I
60    DO 100 I=ILOW,IHIGH
      DO 100 J=JLOW,JHIGH
100   S(J,I)=0.D0
      GO TO 200
160   IERR=.TRUE.
      WRITE(6,1003) COUL0
1003  FORMAT('NUMBER OUT OF RANGE IN FOLLOWING DELETION CARD'/,20I5/
     *'NO FURTHER DELETIONS, BUT SCANNING FOR ZERO TO TERMINATE CARD RE
     *ADING.')
200   CONTINUE
      GO TO 10
300   WRITE(6,1004)
1004  FORMAT('END OF FILE IN DELETION CARDS, NO FURTHER DELETIONS.')
400   RETURN
      END
           DOUBLE PRECISION FUNCTION VALMAD(A,B,R)
C
C     FUNCTION ROUTINE FOR CALCULATING MADELUNG PARAMETERS.
C
           IMPLICIT REAL*8(A-H,O-Z)
           IF(A.LT.0.01D0.OR.B.LT.0.01D0) GO TO 1
           AB=(A+B)/(2.0D0*A*B)
           VALMAD=1.0D0/SQRT(R*R+AB*AB)
           GO TO 2
1          VALMAD=0.0D0
           IF(R.GT.0.001D0) VALMAD=1.0D0/R
2          RETURN
           END
      SUBROUTINE PEGLEG(A,N,NL)
C
C  SUBROUTINE TO PRINT OUT MATRICES IN READABLE FORMAT.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NL,NL)
      NROW=N
      NCOL=N
      GO TO 10
      ENTRY OUTMAT(A,NL,NR,NC)
      NROW=NR
      NCOL=NC
10    KITE=0
20    LOW=KITE+1
      KITE=KITE+14
      IF(KITE.GT.NCOL) KITE=NCOL
      WRITE(6,1000) (I,I=LOW,KITE)
1000  FORMAT(/5X,14I8,//)
      DO 30 I=1,NROW
30    WRITE(6,1001) I,(A(I,J),J=LOW,KITE)
1001  FORMAT(I5,2X,14F8.4)
      IF(KITE.LT.NCOL) GO TO 20
      RETURN
      END
      SUBROUTINE HUCKEL(H,S,MAD,C,SP,PD,MAXS,MAXP,MAXD,COUL0,SORB,IOCC,
     1 HDG,     NDIM, ND1 ,NC,NATOM,NTYPE,NHDG)
C
C  SUBROUTINE TO 1) DETERMINE ORBITAL OCCUPATION NUMBERS, AND 2) SETUP
C  HUCKEL MATRIX.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 MAD
      REAL*4 IOCC
      LOGICAL*4 L1,L2,L3,L4,L5,ONEMAT,ITERAT
      INTEGER*4 SYMBOL,KEY,VELEC
      REAL*8 LAMPRI
      INTEGER*4 PRTCYC
      LOGICAL*4 PARTIT,PRINTX,ITABLE
      INTEGER*4 STAR,STAR2
      INTEGER*4 ONE,TWO,STAR1,TRUE
      DIMENSION H(NDIM,NDIM),S(NDIM,NDIM),MAD(NTYPE,NTYPE),C(NC),
     1 SP(NDIM),PD(NDIM),MAXS(NDIM),MAXP(NDIM),MAXD(NDIM),COUL0(NATOM),
     2 SORB(NDIM),IOCC(NDIM),HDG(NHDG)
      COMMON/TITLE/AB(10)
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,
     1 IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT
      COMMON/OUT/PRT(20),PUN(20),IOVPOP(24),IENRGY(24)
      LOGICAL*4 PRT,PUN
      COMMON/ATOM/KEY(1000),SYMBOL(100),VELEC(40),NS(40),NP(40),ND(40),
     1 EXPS(40),EXPP(40),EXPD(40),EXPD2(40),C1(40),C2(40),COULS(40),
     2 COULP(40),COULD(40),X(1000),Y(1000),Z(1000)
      COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,PRTCYC,
     1 ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)
      COMMON/STARS/STAR,STAR2
      DATA ONE,TWO,STAR1/'1','2','*'/
      IF(ONEMAT.AND.IPRINT.LT.-2) RETURN
      ICYCLE=1
C
C  SETUP DEFAULT ORBITAL OCCUPATIONS.
C
      IF(L1) GO TO 3
      IH=NELEC/2
      JH=NDIM+1-IH
      DO 1 I=1,JH
1     IOCC(I)=0.0
      DO 2 I=JH,NDIM
2     IOCC(I)=2.0
      IF(IH+IH.NE.NELEC) IOCC(JH-1)=1.0
      GO TO 500
C
C  PROVISION FOR READING IN USER SPECIFIED OCCUPATIONS.
C  ALSO PROVISION FOR NON-INTEGER OCCUPATIONS.
C
3     READ(5,2000) (MAXD(I),I=1,NDIM)
2000  FORMAT(80A1)
      TEL=0.0D0
      DO 4 I=1,NDIM
      J=NDIM+1-I
      IOCC(J)=0.0
      IF(MAXD(I).EQ.ONE) IOCC(J)=1.0
      IF(MAXD(I).EQ.TWO) IOCC(J)=2.0
      IF(MAXD(I).EQ.STAR1) READ(5,900) IOCC(J)
900   FORMAT(F15.8)
4     TEL=TEL+REAL(IOCC(J))
      TEL=ABS(TEL-FLOAT(NELEC))
      IF(TEL.LT.0.001D0) GO TO 500
      WRITE(6,2001)
2001  FORMAT('*** WARNING ****    ORBITAL POPULATIONS INCONSISTENT WITH
     1 ASSUMED CHARGE ON MOLECULE',////,T10,'I',T25,'IOCC(I)',/)
      DO 99 I=1,NDIM
99    WRITE(6,2002) I,IOCC(I)
2002  FORMAT(T8,I3,T22,F12.8)
      STOP
C
C  CALL GRMSCH TO ORTHOGONALISE BASIS SET.
C
500   CALL GRMSCH(S,C,NDIM)
      CON=.5D0*CON
      IF(.NOT.ITERAT) GO TO 15
      DO 5 I=1,NDIM
5     SORB(I)=0.D0
      DO 10 I=1,NATOM
      X(I)=0.D0
      Z(I)=0.D0
10    COUL0(I)=0.D0
15    IF(.NOT.ONEMAT) GO TO 25
C
C  IN ONE MATRIX CASE, STUFF DIAGONAL ELEMENTS OF S INTO SP.
C
      DO 20 I=1,NDIM
20    SP(I)=S(I,I)
C
C  SETUP DIAGONAL ELEMENTS OF HUCKEL MATRIX IN H(I,J).
C
25    IF(NH.EQ.0) GO TO 35
      ET=COULH
      DO 30 I=1,NH
      IF(ITERAT) ET=COULH-COUL0(I)
30    C(I)=ET
35    IH=NH+1
      ET=0.D0
      DO 50 I=1,NA
      KEYI=KEY(I)
      IF(ITERAT) ET=COUL0(I+NH)
      C(IH)=COULS(KEYI)-ET
      IH=IH+1
      IF(NP(KEYI).EQ.0) GO TO 50
      HH=COULP(KEYI)-ET
      JH=IH+2
      ASSIGN 40 TO IL
      GO TO 42
40    IF(ND(KEYI).EQ.0) GO TO 50
      HH=COULD(KEYI)-ET
      JH=IH+4
      ASSIGN 50 TO IL
42    DO 45 J=IH,JH
45    C(J)=HH
      IH=JH+1
      GO TO IL,(40,50)
50    CONTINUE
      DO 55 I=1,NDIM
55    H(I,I)=C(I)
      IF(NHDG.EQ.1) GO TO 59
      DO 56 I=1,NDIM
56    HDG(I)=C(I)
C
C  SETUP OFF-DIAGONAL ELEMENTS OF HUCKEL MATRIX.
C
59    CNST=CON
      DO 58 I=2,NDIM
      IL=I-1
      DO 58 J=1,IL
      HH=C(I)+C(J)
      IF(.NOT.L5) GO TO 58
      ET=(C(I)-C(J))/HH
      ET=ET*ET
      CNST=CON+ET/2.0D0+ET*ET*(0.5D0-CON)
58    H(I,J)=CNST*HH*S(I,J)
      IF(ONEMAT) GO TO 100
      DO 60 I=2,NDIM
      IL=I-1
      DO 60 J=1,IL
60    H(J,I)=H(I,J)
C
C  PRINT OUT HUCKEL MATRIX. PRINT OUT TITLE IF METH IS NOT
C  EQUAL TO ZERO.
C
806   IF(ICYCLE.GT.1) GO TO 800
      IF(PRT(6)) GO TO 805
      IF(METH.EQ.0) GO TO 801
      GO TO 802
800   IF(ICYCLE.GE.10000) WRITE(6,701) AB
701   FORMAT('RESULTS OF CALCULATION0  ',8A8,A6,A2,//)
      IF(PRT(7).OR..NOT.PRINTX) GO TO 805
801   WRITE(6,803)
803   FORMAT('HUCKEL MATRIX')
      CALL PEGLEG(H,NDIM,NDIM)
      GO TO 805
802   WRITE(6,804)
804   FORMAT('INPUT HUCKEL MATRIX')
      CALL PEGLEG(H,NDIM,NDIM)
805   IF(ICYCLE.EQ.1.AND.PUN(6)) WRITE(7,825) H
      IF(ICYCLE.GT.1.AND.PUN(7).AND.PRINTX) WRITE(7,825) H
825   FORMAT(8F9.5)
C
C  IF CALCULATING ENERGY MATRIX, STORE H(I,I) IN X(I),Y(I),Z(I).
C
      IF(ICYCLE.LE.MAXCYC.AND.METH.NE.0) GO TO 100
      IF(NH.EQ.0) GO TO 369
      DO 370 I=1,NH
370   X(I)=H(I,I)
369   IH=NH+1
      JH=NH+1
      DO 371 I=1,NA
      KEYI=KEY(I)
      X(IH)=H(JH,JH)
      JH=JH+1
      IF(NP(KEYI).EQ.0) GO TO 371
      Y(IH)=H(JH,JH)
      JH=JH+3
      IF(ND(KEYI).EQ.0) GO TO 371
      Z(IH)=H(JH,JH)
      JH=JH+5
371   IH=IH+1
C
C  CALL TRNFRM TO TRANSFORM HUCKEL MATRIX TO ORTHOGONAL BASIS SET.
C  THEN CALL GIVENS TO PERFORM DIAGONALIZATION.
C
100   IH=1
      IF(ONEMAT) IH=2
      CALL TRNFRM(S,H,C,COUL0,NDIM,SP,IH)
      IF(ONEMAT) GO TO 110
      IH=NDIM
      GO TO 120
110   IH=-NDIM
*120   CALL GIVENS(NDIM,IH,NDIM,C,SP,COUL0,H)
120   Call PreRSP(C,H,Ndim,Coul0)
130   IF(ICYCLE.GE.10000) ITERAT=.FALSE.
C
C  PRINT OUT TITLE, ENERGY LEVELS, AND OCCUPATION NUMBERS.
C
      IF(ITERAT) GO TO 700
      IF(METH.EQ.0) WRITE(6,701) AB
      IF(PRT(8)) GO TO 710
*  Print also energies in atomic units
      enconv = 1./27.21161

      WRITE(6,702) 
     1    ((NDIM-I+1),COUL0(I),COUL0(i)*enconv,IOCC(I),i=1,ndim)
702   FORMAT(////,10X,'ENERGY LEVELS (EV)     and   (AU)',/,
     1(/10X,'E(',I4,') =',F12.5,4x,F12.5,
     18X,F6.4))
710   IF(PUN(8)) WRITE(7,825) (COUL0(I),I=1,NDIM)
700   IF(ONEMAT) GO TO 200
C
C  DIDDLE WITH C,H.
C
      DO 160 J=1,NDIM
      DO 140 K=1,NDIM
140   C(K)=H(K,IH)
      DO 155 I=1,NDIM
      ET=0.D0
      DO 150 K=I,NDIM
150   ET=ET+S(I,K)*C(K)
155   H(I,IH)=ET
160   IH=IH-1
      K=1
      DO 180 I=2,NDIM
      IL=I-1
      DO 180 J=1,IL
      C(K)=S(I,J)
180   K=K+1
200   IF(METH.GT.1.AND.ITERAT) GO TO 210
C
C  CALL OUTPUT FOR FINAL PRINT OUT OF RESULTS.
C
205   CALL OUTPUT(H,S,MAD,C,COUL0,SORB,IOCC,HDG,NDIM,NTYPE,NC,NHDG)
      IF(.NOT.ITERAT) GO TO 999
      GO TO 220
C
C  IF DOING CHARGE ITERATIVE CALCULATION ( METH >1 ), CALL ITRATE
C  TO SETUP HUCKEL MATRIX.
C
210   CALL ITRATE(H,S,MAD,C,COUL0,SORB,IOCC,HDG,NDIM,NTYPE,NC,NHDG)
      IF(.NOT.ITERAT) GO TO 205
220   IF(ICYCLE.GT.MAXCYC) ICYCLE=10000
      IF(METH.GT.1) GO TO 806
      GO TO 15
  999 RETURN
      END
      SUBROUTINE GIVENS (NX,NROOTX,NJX,A,B,ROOT,VECT)
C
C      SUBROUTINE TO CALCULATE THE EIGENVALUES AND EIGENVECTORS
C      OF A REAL SYMMETRIC MATRIX.
C
C
C      THE PARAMETERS FOR THE ROUTINE ARE0
C
C          NX     ORDER OF MATRIX.
C
C          NROOTX NUMBER OF ROOTS FOR WHICH EIGENVECTORS ARE WANTED.
C                 IF NO VECTORS ARE WANTED, MAKE NROOTX NEGATIVE.
C
C          NJX    ROW DIMENSION OF VECT ARRAY.  SEE 'VECT' BELOW.
C                 NJX MUST BE NOT LESS THAN NX.
C
C          A      MATRIX STORED BY COLUMNS IN PACKED UPPER TRIANGULAR
C                 FORM, I.E. OCCUPYING NX*(NX+1)/2 CONSECUTIVE
C                 LOCATIONS.
C
C          B      SCRATCH ARRAY USED BY GIVENS.  MUST BE AT LEAST
C                 NX*6 CELLS.
C
C          ROOT   ARRAY TO HOLD THE EIGENVALUES.  MUST BE AT LEAST
C                 NX CELLS LONG.  THE ROOTS ARE ORDERED LARGEST FIRST
C                 IN THIS ARRAY.
C
C          VECT   EIGENVECTOR ARRAY.  EACH COLUMN WILL HOLD AN
C                 EIGENVECTOR FOR THE CORRESPONDING ROOT.  MUST BE
C                 DIMENSIONED WITH 'NJX' ROWS AND AT LEAST 'NJX'
C                 COLUMNS, UNLESS NO VECTORS ARE REQUESTED (NEGATIVE
C                 NROOTX).  IN THIS LATTER CASE, THE ARGUMENT VECT
C                 IS JUST A DUMMY, AND THE STORAGE IS NOT USED.
C                 THE EIGENVECTORS ARE NORMALIZED TO UNIT LENGTH.
C
C      THE ARRAYS A AND B ARE DESTROYED BY THE COMPUTATION.  THE
C      RESULTS APPEAR IN ROOT AND VECT.
C
C      FOR PROPER FUNCTIONING OF THIS ROUTINE, THE RESULT OF A FLOATING
C      POINT UNDERFLOW SHOULD BE A ZERO.
C
C      THE ORIGINAL REFERENCE TO THE GIVENS TECHNIQUE IS IN OAK RIDGE
C      REPORT NUMBER ORNL 1574 (PHYSICS), BY WALLACE GIVENS.
C
C      THE METHOD AS PRESENTED IN THIS PROGRAM CONSISTS OF FOUR STEPS0
C
C      FIRST, THE INPUT MATRIX IS REDUCED TO TRIDIAGONAL FORM BY THE
C      HOUSEHOLDER TECHNIQUE (J. H. WILKINSON, COMP. J. 3, 23 (1960)).
C      THE EIGENVALUES OF THE TRIDIAGONAL MATRIX ARE THEN FOUND USING
C      THE QR TRANSFORM METHOD.  SEE J. H. WILKINSON, THE ALGEBRAIC
C      EIGENVALUE PROBLEM(1965) FOR A DESCRIPTION OF THIS ALGORITHM.
C      THE EIGENVECTORS OF THE TRIDIAGONAL FORM ARE THEN EVALUATED
C      (J. H. WILKINSON, COMP. J. 1, 90 (1958)), BY THE METHOD OF
C      INVERSE ITERATION, FOR NONDEGENERATE MATRICES.
C      FOR MATRICES WITH DEGENERATE OR NEAR-DEGENERATE EIGENVALUES,
C      THE EIGENVECTORS ARE EVALUATED INSTEAD BY FURTHER QR TRANSFORMS.
C      THIS METHOD GIVES ORTHOGONAL VECTORS EVEN FOR DEGENERATE ROOTS.
C      FINALLY THE TRIDIAGONAL VECTORS ARE ROTATED TO VECTORS OF THE
C      ORIGINAL ARRAY (FIRST REFERENCE).
C
C      THE INVERSE ITERATION PORTION OF THIS PROGRAM WAS ADAPTED
C      FROM THE QUANTUM CHEMISTRY PROGRAM EXCHANGE NUMBER 62.1, BY
C      FRANKLIN PROSSER.  THE EIGENVALUE SUBROUTINE (EVQR) WAS WRITTEN
C      BY WALTER NIELSEN.
C                                      ROY GORDON, SEPT. 1969
C
C      AN EXCELLENT PRESENTATION OF THE GIVENS TECHNIQUE IS FOUND IN
C      J. M. ORTEGA'S ARTICLE IN 'MATHEMATICS FOR DIGITAL COMPUTERS,'
C      VOLUMD 2, ED. BY RALSTON AND WILF, WILEY (1967), PAGE 94.
C
      IMPLICIT DOUBLE PRECISION(A-H,R-Z)
      COMMON /VECTOR/ FACT,IDIF
      DIMENSION B(NX,6),A(1),ROOT(NX),VECT(NJX,NROOTX)
C
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  *                                                                 *
C  *   USERS PLEASE NOTE0 TWO PARAMETERS, ETA AND THETA, SHOULD BE   *
C  *   ADJUSTED BY THE USER FOR HIS PARTICULAR MACHINE.              *
C  *                                                                 *
C  *   ETA IS AN INDICATION OF THE PRECISION OF THE FLOATING POINT   *
C  *   REPRESENTATION ON THE COMPUTER BEING USED (ROUGHLY 10**(-M),  *
C  *   WHERE M IS THE NUMBER OF DECIMALS OF PRECISION ).             *
C  *                                                                 *
C  *   THETA IS AN INDICATION OF THE RANGE OF NUMBERS THAT CAN BE    *
C  *   EXPRESSED IN THE FLOATING POINT REPRESENTATION (ROUGHLY THE   *
C  *   LARGEST NUMBER).                                              *
C  *                                                                 *
C  *   SOME RECOMMENDED VALUES FOLLOW.                               *
C  *                                                                 *
C  *   FOR IBM 7094, UNIVAC 1108, ETC. (27-BIT BINARY FRACTION,      *
C  *   8-BIT BINARY EXPONENT), ETA=1.E-8, THETA=1.E37.               *
C  *   FOR CONTROL DATA 3600 (36-BIT BINARY FRACTION, 11-BIT BINARY  *
C  *   EXPONENT), ETA=1.E-11, THETA=1.E307.                          *
C  *   FOR CONTROL DATA 6600 (48-BIT BINARY FRACTION, 11-BIT BINARY  *
C  *   EXPONENT), ETA=1.E-14, THETA=1.E307.                          *
C  *   FOR IBM 360/50 AND 360/65 DOUBLE PRECISION (56-BIT HEXA-      *
C  *   DECIMAL FRACTION, 7-BIT HEXADECIMAL EXPONENT), ETA=1.E-16,    *
C  *   THETA=1.E75.                                                  *
C  *                                                                 *
C  *   OTHER PARAMETERS WHICH MUST BE ADJUSTED ARE0                  *
C  *                                                                 *
C  *   DEL1 = ETA/1.D2, DELTA = ETA**2*1.D2, SMALL = ETA**2/1.D2,    *
C  *   DELBIG = THETA*DELTA/1.D3, THETA1 = 1.D3/THETA, EMAG = ETA,   *
C  *   TOLER = 1.D2*DSQRT(ETA)                                       *
C  *                                                                 *
C  *   TOLER IS A FACTOR USED TO DETERMINE IF ANY ROOTS ARE CLOSE    *
C  *   ENOUGH TOGETHER TO BE CONSIDERED DEGENERATE FOR PURPOSES OF   *
C  *   CALCULATING EIGENVECTORS.  FOR THE MATRIX NORMED TO UNITY, IF *
C  *   THE DIFFERENCE BETWEEN TWO ROOTS IS LESS THAN TOLER, THEN THE *
C  *   QR TRANSFORMATION IS USED TO FORM THE EIGENVECTORS.           *
C  *                                                                 *
C  *   EMAG IS A TOLERANCE FOR NEGLIGIBLE ELEMENTS IN THE QR         *
C  *   ITERATION FOR EIGENVECTORS FOR DEGENERATE EIGENVALUES.        *
C  *                                                                 *
C  *   IN THE FOLLOWING ROUTINE, ETA = 1.D-16 AND THETA = 1.D75.     *
C  *                                                                 *
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DATA ETA/1.D-16/,THETA/1.D38/,DEL1/1.D-18/,DELTA/1.D-30/,
     *SMALL/1.D-34/,DELBIG/1.D38/,THETA1/1.D-38/,TOLER/1.D-6/,
     *EMAG/1.D-16/
      N = NX
      FLOATN = FLOAT(N)
      NROOT = IABS(NROOTX)
      IF (NROOT.EQ.0) GO TO 1001
      IF (N-1) 1001,1003,105
 1003 ROOT(1) = A(1)
      IF (NROOTX.GT.0) VECT(1,1) = 1.D0
      GO TO 1001
  105 CONTINUE
C
C  NSIZE IS THE NUMBER OF ELEMENTS IN THE PACKED ARRAY.
C
      NSIZE = (N*(N+1))/2
      NM1 = N-1
      NM2 = N-2
      NP1 = N+1
C
C  COMPUTE TRACE.
C
      TRACE = 0.D0
      JUMP = 1
      DO 1 J=2,NP1
      TRACE = TRACE + A(JUMP)
    1 JUMP = JUMP + J
      TRACE = TRACE/FLOATN
C
C  SUBTRACT TRACE FROM DIAGONAL ELEMENTS TO GIVE A MORE RELIABLE NORM
C  WHEN THERE ARE LARGE DIAGONAL ELEMENTS.
C
      JUMP = 1
      DO 2 J=2,NP1
      A(JUMP) = A(JUMP) - TRACE
    2 JUMP = JUMP + J
C
C  SCALE MATRIX TO EUCLIDEAN NORM OF 1.  SCALE FACTOR IS ANORM.
C
      FACTOR = 0.D0
      DO 70 I=1,NSIZE
   70 FACTOR = MAX(FACTOR,ABS(A(I)))
      IF (FACTOR.NE.0.D0) GO TO 72
C
C  NULL MATRIX.  FIX UP ROOTS AND VECTORS, THEN EXIT.
C
      DO 78 I=1,NROOT
      IF (NROOTX.LT.0) GO TO 78
      DO 77 J=1,N
   77 VECT(J,I) = 0.D0
      VECT(I,I) = 1.D0
   78 ROOT(I) = 0.D0
      GO TO 1001
   72 ANORM = 0.D0
   86 SCALE = 1.D0/FACTOR
      DO 80 I=1,NSIZE
   80 ANORM = ANORM + (A(I)*SCALE)**2
      ANORM = ANORM+ANORM
C
C  SUBTRACT DIAGONAL CONTRIBUTIONS WHICH WERE COUNTED TWICE.
C
      JUMP = 1
      DO 81 J=2,NP1
      ANORM = ANORM -(A(JUMP)*SCALE)**2
   81 JUMP = JUMP + J
   83 ANORM = FACTOR*SQRT(ANORM)
      SCALE = 1.D0/ANORM
      DO 91 I=1,NSIZE
   91 A(I) = A(I)*SCALE
      ALIMIT = 1.D0
C
C  TRIDIAGONALIZATION OF SYMMETRIC MATRIX.
C
      ID = 0
      IA = 1
      IF (NM2.EQ.0) GO TO 201
      DO 200 J=1,NM2
C
C  J COUNTS ROW OF A MATRIX TO BE DIAGONALIZED. IA INDICATES START OF
C  NON-CODIAGONAL ELEMENTS IN THE ROW. ID IS THE INDEX OF CODIAGONAL
C  ELEMENT ON THE ROW BEING CODIAGONALIZED.
C
      IA = IA+J+2
      ID = ID+J+1
      JP2 = J + 2
      J1 = J + 1
C
C  FIND LIMITS FOR BAND OF SIGNIFICANT MATRIX ELEMENTS.
C
      LIMIT = J1
      II = IA
      DO 99 I=JP2,N
      B(I,5) = A(II)
      IF (ABS(B(I,5)).GT.DEL1) LIMIT = I
   99 II = II + I
      DTEMP = A(ID)
      IF (LIMIT.GT.J1) GO TO 110
C
C  NO TRANSFORMATION NECESSARY IF ALL THE NON-CODIAGONAL
C  ELEMENTS ARE TINY.
C
  120 B(J,1) = DTEMP
      A(ID) = 0.D0
      GO TO 200
C
C  SUM SQUARES OF SIGNIFICANT NON-CODIAGONAL ELEMENTS OF ROW J.
C
  110 IDIF = LIMIT -JP2
      SUM = DOT(B(JP2,5),B(JP2,5))
C
C  NOW COMPLETE THE SUM OF OFF-DIAGONAL SQUARES.
C
      SUM = SQRT(SUM + DTEMP**2)
C
C  NEW CODIAGONAL ELEMENT.
C
      B(J,1) = -SIGN(SUM,DTEMP)
C
C  FIRST NON-ZERO ELEMENT OF THIS W-VECTOR.
C
      B(J+1,2) = SQRT((1.D0 + ABS(DTEMP)/SUM)*5.D-1)
C
C  FORM REST OF THE W-VECTOR ELEMENTS.
C
      TEMP = SIGN(5.D-1/(B(J+1,2)*SUM),DTEMP)
      II = IA
      DO 130 I=JP2,LIMIT
      B(I,2) = A(II)*TEMP
  130 II = II + I
C
C  FORM P-VECTOR AND SCALAR.  P-VECTOR = A-MATRIX*W-VECTOR.
C  SCALAR = W-VECTOR*P-VECTOR.
C
      DAK = 0.D0
C
C  IC IS THE LOCATION OF THE NEXT DIAGONAL ELEMENT. I RUNS OVER THE
C  NON-ZERO P-ELEMENTS. CASES FOR I LESS THAN LIMIT.
C
      IC = ID + 1
      LIMLES = LIMIT - 1
      DO 188 I=J1,LIMLES
C
C  FORM FIRST PART OF P ELEMENT THEN MOVE IC TO TOP OF NEXT
C  A-MATRIX 'ROW'.
C
      IDIF = I - J1
      DTEMP = DOT(B(J1,2),A(IC))
      IC = IC + I
C
C  COMPLETE P ELEMENT. CHANGE INCREMENTING MODE AT DIAGONAL ELEMENT.
C
  178 IP1 = I + 1
      JJ = IC + IDIF
      DTEMP = DTEMP + DSUM(B(N,1),A(JJ),IP1,LIMIT)
C
C  BUILD UP THE K-SCALAR (AK).
C
      DAK = DAK + DTEMP*B(I,2)
  188 B(I,1) = DTEMP
C
C  CASE FOR I = LIMIT.
C
      IDIF = LIMIT - J1
      DTEMP = DOT(B(J1,2),A(IC))
      DAK = DAK + DTEMP*B(LIMIT,2)
      B(LIMIT,1) = DTEMP
      IDIF = LIMIT - J1
C
C  TEST TO SEE IF ANY I VALUES REMAIN. DO REMAINING VALUES.
C
      IF (LIMIT.EQ.N) GO TO 190
      IC = IC + LIMIT
      LIMLO = LIMIT + 1
      DO 189 I=LIMLO,N
      B(I,1) = DOT(B(J1,2),A(IC))
      B(I,2) = 0.D0
  189 IC = IC + I
C
C  FORM THE Q-VECTOR.
C
  190 FACT = -DAK
      CALL VECSUM(B(J1,1),B(J1,2))
C
C  TRANSFORM THE REST OF THE A-MATRIX. JJ INDICATES START-1 OF THE
C  REST OF THE A-MATRIX. MOVE W-VECTOR INTO THE OLD A-MATRIX LOCATIONS
C  TO SAVE SPACE. I RUNS OVER THE SIGNIFICANT ELEMENTS OF THE W-VECTOR.
C
      JJ = ID
      DO 160 I=J1,N
      A(JJ) = B(I,2)
      IF (I.GT.LIMIT) GO TO 161
      B2 = B(I,2)
      FACT = -B2 - B2
      IDIF = I - J1
      CALL VECSUM(A(JJ+1),B(J1,1))
  161 B1 = B(I,1)
      FACT = -B1 - B1
      IDIF = MIN0(I,LIMIT) - J1
      CALL VECSUM(A(JJ+1),B(J1,2))
  160 JJ = JJ + I
C
C  STORE AWAY LIMIT FOR LATER USE IN BACK TRANSFORMATION. MOVE LAST
C  CODIAGONAL ELEMENT OUT INTO ITS PROPER PLACE.
C
  200 B(J,6) = LIMIT
  201 CONTINUE
      B(NM1,1) = A(NSIZE-1)
      A(NSIZE-1) = 0.D0
C
C  USE QR TRANSFORM METHOD TO FIND EIGENVALUES OF THE TRIDIAGONAL
C  MATRIX. MOVE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX INTO
C  ROOT ARRAY. THIS IS A MORE CONVENIENT INDEXING POSITION. ALSO,
C  PUT SQUARE OF CODIAGONAL ELEMENTS IN THIRD N ELEMENTS.
C
      JUMP = 0
      DO 320 J=1,NM1
      JUMP = JUMP + J
      ROOT(J) = A(JUMP)
  320 B(J,3) = B(J,1)**2
      ROOT(N) = A(NSIZE)
      CALL EVQR(ROOT,B(1,3),N,30,SMALL)
C
C  ROOT NOW CONTAINS THE SHIFTED AND SCALED EIGENVALUES. STORE
C  EIGENVALUES FOR POSSIBLE LATER USE AS SHIFTS IN EVALUATING
C  EIGENVECTORS FOR DEGENERATE MATRICES.
C
      DO 325 J=1,N
  325 B(J,2) = ROOT(J)
C
C  SORT THE EIGENVALUES INTO DESCENDING ALGEBRAIC ORDER.
C
      DO 330 I=1,NM1
      IP1 = I + 1
      DO 330 J=IP1,N
      IF (ROOT(I).GE.ROOT(J)) GO TO 330
      TEMP = ROOT(I)
      ROOT(I) = ROOT(J)
      ROOT(J) = TEMP
  330 CONTINUE
C
C  QUIT NOW IF NO VECTORS WERE REQUESTED. OTHERWISE, TEST FOR
C  DEGENERACY OR NEAR DEGENERACY OF EIGENVALUES FOR WHICH EIGENVECTORS
C  WERE REQUESTED. IF ONLY ONE VECTOR REQUESTED, DEGENERACY DOESN'T
C  MATTER.
C
      IF (NROOTX.LT.0) GO TO 1002
      IF (NROOTX.EQ.1) GO TO 807
      NTOP = NROOT - 1
      DO 400 I=1,NTOP
      IF (ABS(ROOT(I+1)-ROOT(I)).LE.TOLER) GO TO 410
  400 CONTINUE
C
C  NEXT STATEMENT IS REACHED IF ALL EIGENVALUES FOR WHICH EIGENVECTORS
C  WERE REQUESTED ARE WELL SEPARATED.
C
      GO TO 807
C
C  THE FOLLOWING IS REACHED IF THERE ARE ANY DEGENERATE CLUSTERS
C  OF EIGENVALUES.  USE FURTHER QR TRANSFORMS TO EVALUATE
C  THE EIGENVECTORS OF THE TRIDIAGONAL MATRIX.  THIS METHOD
C  GIVES ORTHOGONAL EIGENVECTORS EVEN WHEN THE EIGENVALUES ARE
C  DEGENERATE.  HOWEVER, IT TAKES MORE ARITHMETIC THAN THE METHOD
C  OF INVERSE ITERATION (AT LEAST AT LARGE N).
C
C  PUT DIAGONAL ELEMENTS OF TRIDIAGONAL MATRIX INTO ROOT.  PUT OFF-
C  DIAGONAL ELEMENTS INTO B(I,3).
C
  410 JUMP = 0
      DO 440 J=1,NM1
      JUMP = JUMP + J
      ROOT(J) = A(JUMP)
  440 B(J,3) = B(J,1)
C
C  LAST DIAGONAL ELEMENT.
C
      ROOT(N) = A(NSIZE)
C
C  INITIALIZE VECTORS TO A UNIT MATRIX.
C
      DO 450 I=1,N
      DO 445 J=1,N
  445 VECT(J,I) = 0.D0
  450 VECT(I,I) = 1.D0
C
C  FORM EIGENVECTORS OF TRIDIAGONAL MATRIX FOR DEGENERATE
C  MATRICES AND TRANSPOSE THE VECTORS.
C
      CALL QRTN(ROOT,B(1,3),VECT,B(1,2),N,25,EMAG,NJX)
      DO 456 I=1,NM1
      IP1 = I + 1
      DO 455 J=IP1,N
      FLIP = VECT(I,J)
      VECT(I,J) = VECT(J,I)
  455 VECT(J,I) = FLIP
  456 CONTINUE
C
C  IF ROOTS WERE NOT LOCATED IN DESCENDING ORDER, INTERCHANGE ROOTS
C  AND VECTORS.
C
      ITOP = NROOT-1
      DO 480 I=1,ITOP
      IP1 = I+1
      DO 480 J=IP1,N
      IF (ROOT(I).GE.ROOT(J)) GO TO 480
      TEMP = ROOT(I)
      ROOT(I) = ROOT(J)
      ROOT(J) = TEMP
      DO 470 K=1,N
      TEMP = VECT(K,I)
      VECT(K,I) = VECT(K,J)
  470 VECT(K,J) = TEMP
  480 CONTINUE
C
C  DEGENERATE VECTORS ARE NOW COMPLETE.
C
      GO TO 940
C
C  EIGENVECTORS OF TRIDIAGONAL MATRIX FOR NONDEGENERATE MATRICES.
C  INITIALIZE VECTOR ARRAY.
C
  807 CONTINUE
      DO 705 I=1,NROOT
      DO 15 J=1,N
   15 VECT(J,I) = 1.D0
  705 CONTINUE
      DO 700 I=1,NROOT
C
C  USE INVERSE ITERATION TO FIND VECTORS.
C
  701 AROOT = ROOT(I)
      ELIM1 = A(1) - AROOT
      ELIM2 = B(1,1)
      JUMP = 1
      DO 750 J=1,NM1
      JUMP = JUMP + J + 1
C
C  GET THE CORRECT PIVOT EQUATION FOR THIS STEP.
C
      IF (ABS(ELIM1).LE.ABS(B(J,1))) GO TO 760
C
C  FIRST (ELIM1) EQUATION IS THE PIVOT THIS TIME.  CASE 1.
C
      B(J,2) = ELIM1
      B(J,3) = ELIM2
      B(J,4) = 0.D0
      TEMP = B(J,1)/ELIM1
      ELIM1 = A(JUMP) - AROOT - TEMP*ELIM2
      ELIM2 = B(J+1,1)
      GO TO 755
C
C  SECOND EQUATION IS THE PIVOT THIS TIME.  CASE 2.
C
  760 B(J,2) = B(J,1)
      B(J,3) = A(JUMP) - AROOT
      B(J,4) = B(J+1,1)
      TEMP = 1.D0
      IF (ABS(B(J,1)).GT.THETA1) TEMP = ELIM1/B(J,1)
      ELIM1 = ELIM2 - TEMP*B(J,3)
      ELIM2 = -TEMP*B(J+1,1)
C
C  SAVE FACTOR FOR THE SECOND ITERATION.
C
  755 B(J,5) = TEMP
  750 CONTINUE
      B(N,2) = ELIM1
      B(N,3) = 0.D0
      B(N,4) = 0.D0
      B(NM1,4) = 0.D0
      ITER = 1
C
C  BACK SUBSTITUTE TO GET THIS VECTOR.
C
  790 L = N + 1
      DO 780 J=1,N
      L = L - 1
  786 CONTINUE
      ELIM1 = VECT(L,I)-VECT(L+1,I)*B(L,3)-VECT(L+2,I)*B(L,4)
C
C  IF OVERFLOW IS CONCEIVABLE, SCALE THE VECTOR DOWN. THIS APPROACH
C  IS USED TO AVOID MACHINE-DEPENDENT AND SYSTEM-DEPENDENT CALLS TO
C  OVERFLOW ROUTINES.
C
      IF (ABS(ELIM1).GT.DELBIG) GO TO 782
      TEMP = B(L,2)
      IF (ABS(B(L,2)).LT.DELTA) TEMP = DELTA
      VECT(L,I) = ELIM1/TEMP
      GO TO 780
  782 DO 784 K=1,N
  784 VECT(K,I) = VECT(K,I)/DELBIG
      GO TO 786
  780 CONTINUE
      GO TO (820,900), ITER
C
C  SECOND ITERATION.
C
  820 ITER = ITER + 1
  890 ELIM1 = VECT(1,I)
      DO 830 J=1,NM1
      IF (B(J,2).EQ.B(J,1)) GO TO 840
C
C  CASE ONE.
C
  850 VECT(J,I) = ELIM1
      ELIM1 = VECT(J+1,I) - ELIM1*B(J,5)
      GO TO 830
C
C  CASE TWO.
C
  840 VECT(J,I) = VECT(J+1,I)
      ELIM1 = ELIM1 - VECT(J+1,I)*TEMP
  830 CONTINUE
      VECT(N,I) = ELIM1
      GO TO 790
C
C  NORMALIZE THE VECTOR.
C
  900 ELIM1 = 0.D0
      DO 904 J=1,N
  904 ELIM1 = MAX(ABS(VECT(J,I)),ELIM1)
      TEMP = 0.D0
      DO 910 J=1,N
      ELIM2 = VECT(J,I)/ELIM1
  910 TEMP = TEMP + ELIM2**2
      TEMP = 1.D0/(SQRT(TEMP)*ELIM1)
      DO 920 J=1,N
      VECT(J,I) = VECT(J,I)*TEMP
      IF (ABS(VECT(J,I)).LT.DEL1) VECT(J,I) = 0.D0
  920 CONTINUE
  700 CONTINUE
C
C  ROTATE THE CODIAGONAL VECTORS INTO VECTORS OF ORIGINAL ARRAY.
C  LOOP OVER ALL THE TRANSFORMATION VECTORS.
C
  940 IF (NM2.EQ.0) GO TO 1002
      JUMP = NSIZE - NP1
      IM = NM1
      DO 950 I=1,NM2
      LIMIT = IDINT(B(IM-1,6))
      J1 = JUMP
C
C  MOVE A TRANSFORMATION VECTOR OUT INTO BETTER INDEXING POSITION.
C
      DO 955 J=IM,LIMIT
      B(J,2) = A(J1)
  955 J1 = J1 + J
      IDIF = LIMIT - IM
C
C  MODIFY ALL REQUESTED VECTORS.
C
      DO 960 K=1,NROOT
C
C  FORM SCALAR PRODUCT OF TRANSFORMATION VECTOR WITH EIGENVECTOR.
C
      TMP = DOT(B(IM,2),VECT(IM,K))
      FACT = -TMP - TMP
      CALL VECSUM(VECT(IM,K),B(IM,2))
  960 CONTINUE
      JUMP = JUMP - IM
  950 IM = IM - 1
 1002 CONTINUE
C
C  RESTORE ROOTS TO THEIR PROPER SIZE AND ADD BACK TRACE.
C
      DO 95 I=1,N
   95 ROOT(I) = ROOT(I)*ANORM + TRACE
 1001 RETURN
      END
      SUBROUTINE EVQR(A,B,N,M,TOL)
C
C  SUBROUTINE FOR PERFORMING A QR TRANSFORM ON A REAL SYMMETRIC
C  TRIDIAGONAL MATRIX.
C
C  ON ENTERING, A CONTAINS THE N DIAGONAL ELEMENTS OF THE TRIDIAGONAL
C  MATRIX. B CONTAINS THE SQUARES OF THE N-1 OFF-DIAGONAL ELEMENTS.
C  ITERATION IS CONTINUED UNTIL THE SQUARES OF THE OFF-DIAGONAL
C  ELEMENTS ARE LESS THAN TOL. TYPICALLY LESS THAN TWO ITERATIONS PER
C  EIGENVALUE ARE REQUIRED. THUS THE UPPER LIMIT M TO THE NUMBER OF
C  ITERATIONS PER EIGENVALUE MAY BE SAFELY SET AT 20 OR SO.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(1),B(1)
      NX = N
      NI = 1
      SH = 0.D0
C
C  K COUNTS THE NUMBER OF ITERATIONS PER EIGENVALUE.
C
      K = 0
      IF (NX-2) 50,60,85
C
C  EACH NEW ITERATION BEGINS HERE.
C
  100 K = K + 1
      IF (K.LT.M) GO TO 101
      WRITE (6,1000) K
 1000 FORMAT(35H NO CONVERGENCE OF QR ALGORITHM IN   ,I4,11H ITERATIONS)
      CALL EXIT
C
C  SOLVE THE TWO BY TWO IN THE LOWER RIGHT CORNER AND USE SMALLER ROOT
C  AS THE NEXT SHIFT.
C
  101 AT = A(NX) + A(NX-1)
      ST = AT*5.D-1
      DISC = AT**2 - 4.D0*(A(NX)*A(NX-1)-B(NX-1))
      IF (DISC.LE.0.D0) GO TO 15
      ST = ST - SIGN(SQRT(DISC),ST)*5.D-1
C
C  INCREASE THE TOTAL SHIFT BY THE TEMPORARY SHIFT.
C
   15 SH = SH + ST
C
C  THIS LOOP SUBTRACTS THE TEMPORARY SHIFT FROM THE DIAGONAL ELEMENTS.
C
      DO 20 I=1,NX
   20 A(I) = A(I) - ST
C
C  INITIALIZE.
C
      G = A(NI)
      PS = G**2
      RS = PS + B(NI)
      SX = B(NI)/RS
      CXS = 1.D0
      CX = PS/RS
      U = SX*(G+A(NI+1))
      A(NI) = G + U
      NTOP = NX - 2
C
C  THIS LOOP COMPLETES ONE ITERATION, THAT IS ONE QR TRANSFORM.
C
      DO 10 I=NI,NTOP
C
C  G IS THE GAMMA IN THE NOTATION OF WILKINSON.
C
      G = A(I+1) - U
      IF (CX.GT.TOL) GO TO 12
      PS = B(I)*CXS
      GO TO 16
   12 PS = G**2/CX
   16 RS = PS + B(I+1)
C
C  ROTATE AN OFF-DIAGONAL ELEMENT.
C
      B(I) = SX*RS
      SX = B(I+1)/RS
      CXS = CX
      CX = PS/RS
      U = SX*(G+A(I+2))
C
C  ROTATE A DIAGONAL ELEMENT.
C
      A(I+1) = G + U
   10 CONTINUE
C
C  COMPUTE THE LAST DIAGONAL ELEMENT.
C
      A(NX) = A(NX) - U
C
C  COMPUTE THE LAST OFF-DIAGONAL ELEMENT.
C
      IF (CX.GT.TOL) GO TO 112
      PS = B(NTOP+1)*CXS
      GO TO 116
  112 PS = ((A(NX))**2)/CX
  116 B(NTOP+1) = SX*PS
C
C  END OF ONE ITERATION.
C
   85 IT = NX
C
C  CHECK UPWARD THROUGH THE OFF-DIAGONAL ELEMENTS TO FIND THOSE LESS
C  THAN TOL. IF NO OFF-DIAGONAL ELEMENTS LESS THAN TOL ARE FOUND,
C  PERFORM ANOTHER ITERATION.
C
   30 IT = IT-1
      IF (ABS(B(IT)).LE.TOL) GO TO 40
      IF (IT-NI) 100,100,30
C
C  BRANCH ACCORDING TO WHETHER THE MATRIX ISOLATED BY THE SMALL
C  OFF-DIAGONAL ELEMENT IS OF DIMENSION ONE, TWO, OR MORE.
C
   40 IF (NX-IT-2) 50,60,70
C
C  EXTRACT THE EIGENVALUE OF A ONE BY ONE MATRIX BY ADDING BACK
C  THE SHIFT.
C
   50 A(NX) = A(NX) + SH
C
C  DECREASE THE SIZE OF THE PORTION OF THE MATRIX AFFECTED BY
C  LATER ITERATIONS.
C
      NX = NX-1
C
C  RESET THE ITERATION COUNTER.
C
      K = 1
      GO TO 80
C
C  EXTRACT THE EIGENVALUES FROM A TWO BY TWO MATRIX.
C
   60 AL = B(NX-1)
      AM = 5.D-1*(A(NX-1)-A(NX))
      AMS = AM**2
      SAM = SIGN(1.0,AM)
      AN = SQRT(AL+AM**2)
      CX = (AN+ABS(AM))/(2.D0*AN)
      SX = B(NX-1)/(4.D0*AN**2*CX)
      TA = A(NX-1)
      TB = A(NX)
      TC = B(NX-1)
      I = NX
C
C  ROTATE THE DIAGONAL ELEMENTS AND THE OFF-DIAGONAL ELEMENTS.
C
      A(NX-1) = TA*CX+TB*SX+TC*SAM/AN+SH
      A(NX) = TA*SX+TB*CX-TC*SAM/AN+SH
      B(NX-1) = 4.D0*AMS*CX*SX-ABS(AM)*TC/AN+TC*(CX-SX)**2
C
C  RESET THE ITERATION COUNTER.
C
      K = 1
C
C  DECREASE THE SIZE OF THE PORTION OF THE MATRIX AFFECTED BY THE LATER
C  ITERATIONS.
C
      NX = NX-2
      GO TO 80
C
C  THE NEXT STATEMENT IS REACHED WHEN THE PORTION OF THE MATRIX
C  ISOLATED IS GREATER THAN TWO BY TWO.  IT CHANGES THE LOWER LIMIT
C  OF THE ITERATION SO THAT ONLY THIS PORTION WILL BE AFFECTED BY
C  SUBSEQUENT ROTATIONS UNTIL ALL ITS EIGENVALUES ARE FOUND.
C
   70 NI = IT + 1
C
C  TRANSFER TO BEGINNING OF ANOTHER ITERATION.
C
      GO TO 85
C
C  NEXT STATEMENT IS REACHED AFTER EITHER ONE OR TWO EIGENVALUES
C  HAVE JUST BEEN FOUND.  IT TRANSFERS IF ALL THE EIGENVALUES IN
C  THIS PORTION OF THE MATRIX HAVE BEEN FOUND.
C
   80 IF (NX.LT.NI) GO TO 90
C
C  BRANCH ACCORDING TO WHETHER ONE, TWO, OR MORE EIGENVALUES REMAIN
C  TO BE FOUND.
C
   95 IF (NX-NI-1) 50,60,85
C
C  THE NEXT STATEMENT IS REACHED WHEN ALL EIGENVALUES IN THIS PART
C  OF THE MATRIX HAVE BEEN FOUND. IT RETURNS IF THIS IS THE LAST
C  PART OF THE MATRIX.
C
   90 IF (NI.EQ.1) RETURN
C
C  ENLARGE THE PORTION OF THE MATRIX BEING TREATED TO INCLUDE THE
C  BEGINNING OF THE MATRIX.
C
      NI = 1
      GO TO 95
      END                               
      SUBROUTINE QRTN (A,B,V,EIG,N,M,TOL,NJX)
C
C  SUBROUTINE FOR PERFORMING A QR TRANSFORM ON A REAL SYMMETRIC
C  TRIDIAGONAL MATRIX.
C
C  N IS THE DIMENSION OF THE MATRIX. A CONTAINS THE N DIAGONAL
C  ELEMENTS OF THE TRIDIAGONAL MATRIX. B CONTAINS THE N-1 OFF-
C  DIAGONAL ELEMENTS. M IS THE MAXIMUM NUMBER OF ITERATIONS ( SAY
C  20 ). NJX IS THE PHYSICAL ROW DIMENSION OF THE EIGENVECTOR
C  MATRIX V.
C
C  THE EIGENVALUES ARE ASSUMED KNOWN AND PLACED IN THE FIRST N
C  ELEMENTS OF EIG.
C
C  N VECTORS V ( EACH OF LENGTH N ) ARE TRANSFORMED INTO THE
C  BASIS IN WHICH THE MATRIX IS DIAGONAL.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),V(1),EIG(1)
      NX = N
      NNM1 = NJX*(NX-1)
      NI = 1
C
C  SET INITIAL TOTAL SHIFT.
C
      SH = 0
      IF (NX-2) 50,60,1
C
C  K COUNTS THE NUMBER OF ITERATIONS PER EIGENVALUE.
C
    1 K=0
C
C  SET INITIAL TEMPORARY SHIFT.
C
   98 ST = EIG(NX)-SH
C
C  CHECK FOR SMALL OFF-DIAGONAL ELEMENTS.
C
      IT = NX
   99 IT = IT-1
      IF (ABS(B(IT)).LE.TOL) GO TO 40
      IF (IT.GT.NI) GO TO 99
C
C  NO SMALL OFF-DIAGONAL ELEMENTS FOUND.  ITERATE.  EACH NEW
C  ITERATION BEGINS HERE.
C
  100 K = K + 1
      IF (K.LT.M) GO TO 11
      WRITE (6,1000) K
 1000 FORMAT(35H NO CONVERGENCE OF QR ALGORITHM IN   ,I4,11H ITERATIONS)
      CALL EXIT
   11 IF (K.EQ.1) GO TO 15
C
C  SOLVE THE TWO BY TWO IN THE LOWER RIGHT CORNER AND USE SMALLER ROOT
C  AS THE NEXT SHIFT.
C
   12 AT = A(NX) + A(NX-1)
      ST = AT*5.D-1
      DISC = AT**2-4.D0*(A(NX)*A(NX-1)-B(NX-1)**2)
      IF (DISC.LE.0.D0) GO TO 15
      ST = ST-SIGN(SQRT(DISC),ST)*5.D-1
C
C  INCREASE THE TOTAL SHIFT BY THE TEMPORARY SHIFT.
C
   15 SH = SH + ST
C
C  THIS LOOP SUBTRACTS THE TEMPORARY SHIFT FROM THE DIAGONAL ELEMENTS.
C
      DO 20 I=1,NX
   20 A(I) = A(I) - ST
      R = SQRT(A(NI)**2+B(NI)**2)
      S = B(NI)/R
      CS = S
      C = A(NI)/R
      U = (S**2)*(A(NI)+A(NI+1))
      A(NI) = A(NI) + U
      CALL ROTATE(V(NI),C,S,NJX,NNM1)
      NTOP = NX - 2
C
C  THIS LOOP COMPLETES ONE ITERATION, THAT IS ONE QR TRANSFORM.
C
      DO 10 I=NI,NTOP
C
C  G IS THE GAMMA AND Q IS THE P IN THE NOTATION OF WILKINSON.
C
      G = A(I+1) - U
      Q = C*A(I+1) - CS*B(I)
      R = SQRT(Q**2+B(I+1)**2)
C
C  ROTATE AN OFF-DIAGONAL ELEMENT.
C
      B(I) = S*R
C
C  FIND THE NEW SINE AND COSINE FOR THE JACOBI ROTATION. THEN
C  COMPUTE A NEW U.
C
      S = B(I+1)/R
      CS = C*S
      C = Q/R
      U = (S**2)*(G+A(I+2))
C
C  ROTATE A DIAGONAL ELEMENT.
C
      A(I+1) = G + U
C
C  ROTATE THE VECTORS.
C
      CALL ROTATE (V(I+1),C,S,NJX,NNM1)
   10 CONTINUE
C
C  COMPUTE THE LAST OFF DIAGONAL ELEMENT.
C
      B(NTOP+1) = S*(C*A(NX)-CS*B(NTOP+1))
C
C  COMPUTE THE LAST DIAGONAL ELEMENT.
C
      A(NX) = A(NX) - U
C
C  END OF ONE ITERATION.
C
   85 IT = NX
C
C  CHECK UPWARD THROUGH THE OFF DIAGONAL ELEMENTS TO FIND THOSE LESS
C  THAN TOL. IF NO OFF-DIAGONAL ELEMENTS LESS THAN TOL ARE FOUND,
C  PERFORM ANOTHER ITERATION.
C
   30 IT = IT - 1
      IF (ABS(B(IT)).LE.TOL) GO TO 40
      IF (IT-NI) 100,100,30
C
C  BRANCH ACCORDING TO WHETHER THE MATRIX ISOLATED BY THE SMALL
C  OFF-DIAGONAL ELEMENT IS OF DIMENSION ONE, TWO, OR MORE.
C
   40 IF (NX-IT-2) 50,60,70
C
C  EXTRACT THE EIGENVALUE OF A ONE BY ONE MATRIX BY ADDING BACK
C  THE SHIFT.
C
   50 A(NX) = A(NX) + SH
C
C  DECREASE THE SIZE OF THE PORTION OF THE MATRIX AFFECTED BY
C  LATER ITERATIONS.
C
      NX = NX - 1
C
C  RESET THE ITERATION COUNTER.
C
      K = 0
      GO TO 80
C
C  EXTRACT THE EIGENVALUES FROM A TWO BY TWO MATRIX AND PERFORM THE
C  CORRESPONDING ROTATIONS ON THE VECTORS.
C
   60 AL = -B(NX-1)
      AM = 5.D-1*(A(NX-1)-A(NX))
      AN = SQRT(AL**2+AM**2)
      C = SQRT((AN+ABS(AM))/(2.D0*AN))
      S = SIGN(5.E-1,AM)*AL/(AN*C)
      TA = A(NX-1)
      TB = A(NX)
      TC = B(NX-1)
      CX = C**2
      SX = S**2
      CS = C*S
C
C  ROTATE THE DIAGONAL ELEMENTS, THE OFF-DIAGONAL ELEMENTS, AND
C  THE VECTORS.
C
      A(NX-1) = TA*CX+TB*SX-2.D0*TC*CS+SH
      A(NX) = TA*SX+TB*CX+2.D0*TC*CS+SH
      B(NX-1) = 2.D0*AM*CS+TC*(CX-SX)
      I = NX-1
      S = -S
      CALL ROTATE (V(I),C,S,NJX,NNM1)
C
C  RESET THE ITERATION COUNTER.
C
      K = 0
C
C  DECREASE THE SIZE OF THE PORTION OF THE MATRIX AFFECTED BY THE
C  LATER ITERATIONS.
C
      NX = NX-2
      GO TO 80
C
C  THE NEXT STATEMENT IS REACHED WHEN THE PORTION OF THE MATRIX
C  ISOLATED IS GREATER THAN TWO BY TWO.  IT CHANGES THE LOWER LIMIT
C  OF THE ITERATION SO THAT ONLY THIS PORTION WILL BE AFFECTED BY
C  SUBSEQUENT ROTATIONS UNTIL ALL ITS EIGENVALUES ARE FOUND.
C
   70 NI = IT + 1
C
C  TRANSFER TO BEGINNING OF ANOTHER ITERATION.
C
      GO TO 100
C
C  NEXT STATEMENT IS REACHED AFTER EITHER ONE OR TWO EIGENVALUES
C  HAVE JUST BEEN FOUND.  IT TRANSFERS IF ALL THE EIGENVALUES IN
C  THIS PORTION OF THE MATRIX HAVE BEEN FOUND.
C
   80 IF (NX.LT.NI) GO TO 90
C
C  BRANCH ACCORDING TO WHETHER ONE, TWO, OR MORE EIGENVALUES REMAIN
C  TO BE FOUND.
C
   95 IF (NX-NI-1) 50,60,98
C
C  THE NEXT STATEMENT IS REACHED WHEN ALL EIGENVALUES IN THIS PART
C  OF THE MATRIX HAVE BEEN FOUND. IT RETURNS IF THIS IS THE LAST
C  PART OF THE MATRIX.
C
   90 IF (NI.EQ.1) RETURN
C
C  ENLARGE THE PORTION OF THE MATRIX BEING TREATED TO INCLUDE THE
C  BEGINNING OF THE MATRIX.
C
      NI = 1
      GO TO 95
      END
      SUBROUTINE ITRATE(H,U,MAD,C,E,W,IOCC,HDG,NDIM,NTYPE,NC,NHDG)
C
C  SUBROUTINE TO SETUP HUCKEL MATRIX WHEN USING CHARGE ITERATION
C  OPTION ( METH = 2 OR 3 ).
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION H(NDIM,NDIM),U(NDIM,NDIM),MAD(NTYPE,NTYPE),C(NC),
     1 E(NDIM),W(NDIM),IOCC(NDIM),HDG(NHDG)
      REAL*8 MAD
      REAL*4 IOCC
      Common/LUL/IhelpV(100),SymbH(1000),XX(1000),YY(1000),ZZ(1000)
     1,MaxDim    
      COMMON/TITLE/AB(10)
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,
     1 IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT
      LOGICAL*4 L1,L2,L3,L4,L5,ONEMAT,ITERAT
      COMMON/OUT/PRT(20),PUN(20),IOVPOP(24),IENRGY(24)
      LOGICAL*4 PRT,PUN
      COMMON/ATOM/KEY(1000),SYMBOL(100),VELEC(40),NS(40),NP(40),ND(40),
     1 EXPS(40),EXPP(40),EXPD(40),EXPD2(40),C1(40),C2(40),COULS(40),
     2 COULP(40),COULD(40),X(1000),Y(1000),Z(1000)
      INTEGER*4 SYMBOL,KEY,VELEC,SymbH
      COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,PRTCYC,
     1 ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)
      REAL*8 LAMPRI
      INTEGER*4 PRTCYC
      LOGICAL*4 PARTIT,PRINTX,ITABLE
C
C  SINCE THE INTERNAL ATOMIC PARAMETERS ( EXPS,EXPP, ETC. ) ARE
C  NOT USED WHEN DOING CHARGE ITERATION, THE SPACE ALLOCATED TO
C  THEM HAS BEEN USED FOR THE VSIE PARAMETERS.
C
      DIMENSION AS1(20),BS1(20),CS1(20),AP1(20),BP1(20),CP1(20),
     1 AD1(20),BD1(20),CD1(20)
      EQUIVALENCE (AS1(1),EXPS(21)),(BS1(1),EXPP(21)),(CS1(1),EXPD(21)),
     1 (AP1(1),EXPD2(21)),(BP1(1),C1(21)),(CP1(1),C2(21)),
     2 (AD1(1),COULS(21)),(BD1(1),COULP(21)),(CD1(1),COULD(21))
      COMMON/ABC/AS2(5),BS2(5),CS2(5),AP2(5),BP2(5),CP2(5),AD2(5),
     1 BD2(5),CD2(5),AS3(5),BS3(5),CS3(5),AP3(5),BP3(5),CP3(5),AD3(5),
     2 BD3(5),CD3(5)
      REAL*4 ZS(20),ZP(20),ZD(20)
      EQUIVALENCE (ZS(1),NS(21)),(ZP(1),NP(21)),(ZD(1),ND(21))
      DIMENSION IL(3),JL(3)
      EQUIVALENCE (PEEP,ADJUST),(COULH,SIGNN),(SENSE,DENSE),
     1 (VELEC(21),DCHG),(PRT(1),NMIN),(PUN(1),NLN)
      DIMENSION F(6)
      EQUIVALENCE (A1,F(1)),(A2,F(3)),(A3,F(5))
      REAL*8 LAMBDA
      QON=FLOAT(KA)/FLOAT(NELEC)
      IF(ICYCLE.GT.1) GO TO 696
      WRITE(6,750)
750   FORMAT('1')
      DO 5 I=1,NDIM
      NMIN=I
      IF(IOCC(I).GT.0.0001) GO TO 697
5     CONTINUE
C
C  IF CHARGE ITERATION WITH MADELUNG CORRECTION ( METH>2 ) IS BEING
C  USED, COMPUTE THE TOTAL ORBITAL OCCUPATIONS FOR THE FREE ATOMS.
C  FILL ATOMIC ORBITALS IN ORDER0 1S 2S 2P 3S 3P 3D 4S 4P 4D 5S 5P
C  5D 6S 6P.
C
697   IF(METH.LT.3) GO TO 696
      DO 700 I=1,NA
      KEYI=KEY(I)
      J=NP(KEYI)
      IF(J.EQ.0) GO TO 701
      IL(2)=J*10+1
      JL(2)=2
      J=NS(KEYI)
      IL(1)=J*10
      JL(1)=1
      J=ND(KEYI)
      IF(J.EQ.0) GO TO 702
      IL(3)=J*10+2
      JL(3)=3
      DO 710 J=2,3
      IJ=4-J
      K=IJ+1
      DO 710 L=1,IJ
      N=K-L
      IF(IL(K).GT.IL(N)) GO TO 710
      I1=IL(K)
      J1=JL(K)
      IL(K)=IL(N)
      JL(K)=JL(N)
      IL(N)=I1
      JL(N)=J1
710   CONTINUE
      GO TO 720
701   JL(1)=1
      JL(2)=2
      JL(3)=3
      GO TO 720
702   JL(3)=3
      IF(IL(1).LT.IL(2)) GO TO 720
      JL(2)=1
      JL(1)=2
720   L=VELEC(KEYI)
      DO 700 J=1,3
      K=JL(J)
      J1=4*(K-1)+2
      IF(L.LT.J1) J1=L
      L=L-J1
      IF(K-2) 725,726,727
725   ZS(KEYI)=FLOAT(J1)
      GO TO 700
726   ZP(KEYI)=FLOAT(J1)
      GO TO 700
727   ZD(KEYI)=FLOAT(J1)
700   CONTINUE
696   IF(PRINTX) WRITE(6,751)
751   FORMAT(////)
      PRINTX=((ICYCLE/PRTCYC)*PRTCYC.EQ.ICYCLE)
      IF(ICYCLE.EQ.MAXCYC) PRINTX=.TRUE.
C
C  CALCULATE SUM OF ONE-ELECTRON ENERGIES.
C
753   SUM=0.0D0
      DO 13 I=1,NDIM
      W(I)=REAL(IOCC(I))
13    SUM=SUM+E(I)*W(I)
C
C  COMPUTE ATOMIC ORBITAL OCCUPATIONS AND STORE IN E(I).
C
3004  IJ=1
      DO  319  I=1,NDIM
      E(I)=0.0D0
      DO  319  J=1,I
      UB=0.0D0
      DO 320 K=NMIN,NDIM
320   UB=UB+H(I,K)*H(J,K)*W(K)
      UB=UB*0.50D0
      IF(I.EQ.J)  GO  TO  28
      UB=(UB+UB)*C(IJ)
      IJ=IJ+1
28    E(I)=E(I)+UB
      E(J)=E(J)+UB
319   CONTINUE
C
C  COMPUTE THE ORBITAL OCCUPATION OF A GIVEN TYPE (S,P,D) WHICH
C  VARIES MOST FROM THE LAST CYCLE.
C
3005  DENOM=0.0D0
      DENSE2=DENSE
      DCHG2=DCHG
      J=1
      DO  5000 I=1,NA
      KEYI=KEY(I)
      SDENSE=E(J)
      SIGNS=SDENSE-X(I)
      C(J)=SIGNS
      DIFFS=ABS(SIGNS)
      N=J
      IF(DENOM.GT.DIFFS) GO TO 5001
      DENSE=X(I)
      DCHG=SIGNS
      DENOM=DIFFS
      NLF=10*I
5001  IF(NP(KEYI).EQ.0) GO TO 5000
      PDENSE=E(J+1)+E(J+2)+E(J+3)
      SIGNP=PDENSE-Y(I)
      C(J+3)=SIGNP
      DIFFP=ABS(SIGNP)
      N=J+3
      IF(DENOM.GT.DIFFP) GO TO 5002
      DENSE=Y(I)
      DCHG=SIGNP
      DENOM=DIFFP
      NLF=10*I+1
5002  IF(ND(KEYI).EQ.0) GO TO 5000
      DDENSE=E(J+4)+E(J+5)+E(J+6)+E(J+7)+E(J+8)
      SIGND=DDENSE-Z(I)
      C(J+8)=SIGND
      DIFFD=ABS(SIGND)
      N=J+8
      IF(DENOM.GT.DIFFD) GO TO 5000
      DENSE=Z(I)
      DCHG=SIGND
      DENOM=DIFFD
      NLF=10*I+2
5000  J=N+1
      SIGNF=DCHG/DENOM
C
C  DETERMINE LAMBDA, THE DAMPING FACTOR.
C
      IF(ICYCLE.LE.2) SIGNN=SIGNF
      IF(ICYCLE.LE.2) NLN=NLF
      ISIGNN=SIGNN
      ISIGNF=SIGNF
      IF(ICYCLE.EQ.2) ADJUST=DAMP1*DENOM
      IF(ISIGNN.EQ.ISIGNF.OR.NLN.NE.NLF) GO TO 8001
      IF(DAMP3.NE.0.0D0) GO TO 401
      LAMBDA=(DENSE2-DENSE)/(DCHG-DCHG2)
      GO TO 402
401   ADJUST=DAMP3*ADJUST
8001  IF((ADJUST/DENOM).GE.LAMPRI) ADJUST=DAMP2*DENOM
      IF(ICYCLE.EQ.1) ADJUST=DENOM
      LAMBDA=ADJUST/DENOM
402   SIGNN=SIGNF
      NLN=NLF
C
C
C  IN THIS SECTION OF THE PROGRAM THREE CALCULATIONS ARE PERFORMED0
C
C    1. THE TOTAL ORBITAL OCCUPATIONS OF A GIVEN TYPE (S,P,D) ARE
C       DAMPED AND STORED (IN X(I),Y(I),Z(I) RESPECTIVELY).
C    2. THE NET CHARGES ARE CALCULATED AND STORED IN C(I).
C    3. -VSIE'S ARE CALCULATED AND STORED IN W(I).
C
C
      J=1
      DO 803 I=1,NA
      KEYI=KEY(I)
      SDENSE=X(I)+LAMBDA*C(J)
      X(I)=SDENSE
      UB=SDENSE
      N=J
      IF(NP(KEYI).EQ.0) GO TO 804
      PDENSE=Y(I)+LAMBDA*C(J+3)
      Y(I)=PDENSE
      UB=UB+PDENSE
      N=J+3
      IF(ND(KEYI).EQ.0) GO TO 805
      N=J+8
      DDENSE=Z(I)+LAMBDA*C(J+8)
      Z(I)=DDENSE
      UB=UB+DDENSE
      GO TO 806
804   Q=VELEC(KEYI)-UB
      IF(.NOT.PARTIT) GO TO 1111
      IF(ITABLE(KEYI)) GO TO 1111
      W(J)=COULS(KEYI)
      GO TO 807
1111  KEYI=21-KEY(I)
      W(J)=-((AS1(KEYI)*Q+BS1(KEYI))*Q+CS1(KEYI))
      GO TO 807
805   Q=VELEC(KEYI)-UB
      IF(.NOT.PARTIT) GO TO 1113
      IF(ITABLE(KEYI)) GO TO 1113
      W(J)=COULS(KEYI)
      W(J+1)=COULP(KEYI)
      W(J+2)=COULP(KEYI)
      W(J+3)=COULP(KEYI)
      GO TO 807
1113  KEYI=21-KEY(I)
      VSIES1=(AS1(KEYI)*Q+BS1(KEYI))*Q+CS1(KEYI)
      W(J)=-VSIES1
      VSIEP1=(AP1(KEYI)*Q+BP1(KEYI))*Q+CP1(KEYI)
      W(J+1)=-VSIEP1
      W(J+2)=-VSIEP1
      W(J+3)=-VSIEP1
      GO TO 807
806   Q=VELEC(KEYI)-UB
      IF(.NOT.PARTIT) GO TO 1115
      IF(ITABLE(KEYI)) GO TO 1115
      W(J)=COULS(KEYI)
      W(J+1)=COULP(KEYI)
      W(J+2)=COULP(KEYI)
      W(J+3)=COULP(KEYI)
      W(J+4)=COULD(KEYI)
      W(J+5)=COULD(KEYI)
      W(J+6)=COULD(KEYI)
      W(J+7)=COULD(KEYI)
      W(J+8)=COULD(KEYI)
      GO TO 807
1115  IF(ND(KEYI).EQ.NP(KEYI)) GO TO 1116
      IF(NCON.NE.3) GO TO 1116
      KEYI=21-KEY(I)
      VSIES1=(AS1(KEYI)*Q+BS1(KEYI))*Q+CS1(KEYI)
      VSIES2=(AS2(KEYI)*Q+BS2(KEYI))*Q+CS2(KEYI)
      VSIES3=(AS3(KEYI)*Q+BS3(KEYI))*Q+CS3(KEYI)
      W(J)=(SDENSE+PDENSE-2.0D0)*VSIES1+(1.0D0-SDENSE)*VSIES2-
     *PDENSE*VSIES3
      VSIEP1=(AP1(KEYI)*Q+BP1(KEYI))*Q+CP1(KEYI)
      VSIEP2=(AP2(KEYI)*Q+BP2(KEYI))*Q+CP2(KEYI)
      VSIEP3=(AP3(KEYI)*Q+BP3(KEYI))*Q+CP3(KEYI)
      W(J+1)=(SDENSE+PDENSE-2.0D0)*VSIEP1+(1.0D0-PDENSE)*VSIEP2-
     *SDENSE*VSIEP3
      W(J+2)=W(J+1)
      W(J+3)=W(J+1)
      VSIED1=(AD1(KEYI)*Q+BD1(KEYI))*Q+CD1(KEYI)
      VSIED2=(AD2(KEYI)*Q+BD2(KEYI))*Q+CD2(KEYI)
      VSIED3=(AD3(KEYI)*Q+BD3(KEYI))*Q+CD3(KEYI)
      W(J+4)=(SDENSE+PDENSE-1.0D0)*VSIED1-SDENSE*VSIED2-PDENSE*VSIED3
      W(J+5)=W(J+4)
      W(J+6)=W(J+4)
      W(J+7)=W(J+4)
      W(J+8)=W(J+4)
      GO TO 807
1116  KEYI=21-KEY(I)
      VSIES1=(AS1(KEYI)*Q+BS1(KEYI))*Q+CS1(KEYI)
      W(J)=-VSIES1
      VSIEP1=(AP1(KEYI)*Q+BP1(KEYI))*Q+CP1(KEYI)
      W(J+1)=-VSIEP1
      W(J+2)=W(J+1)
      W(J+3)=W(J+1)
      VSIED1=(AD1(KEYI)*Q+BD1(KEYI))*Q+CD1(KEYI)
      W(J+4)=-VSIED1
      W(J+5)=W(J+4)
      W(J+6)=W(J+4)
      W(J+7)=W(J+4)
      W(J+8)=W(J+4)
807   J=N+1
      C(I)=Q
803   CONTINUE
C
C  IF CHARGE ITERATION WITHOUT MADELUNG CORRECTION ( METH=2 ) IS
C  BEING USED, SETUP HUCKEL MATRIX. OTHERWISE SKIP THIS SECTION.
C
      IF(METH.GT.2) GO TO 999
      H(1,1)=W(1)
      CNST=CON
      DO 760 I=2,NDIM
      H(I,I)=W(I)
      J1=I-1
      DO 760 J=1,J1
      UB=W(I)+W(J)
      IF(.NOT.L5) GO TO 761
      UC=(W(I)-W(J))/UB
      UC=UC*UC
      CNST=CON+UC/2.0D0+UC*UC*(0.5D0-CON)
761   UB=CNST*U(I,J)*UB
      H(J,I)=UB
760   H(I,J)=UB
      GO TO 850
C
C  IF CHARGE ITERATION WITH MADELUNG CORRECTION ( METH>2 ) IS BEING
C  USED, SETUP HUCKEL MATRIX. OTHERWISE SKIP THIS SECTION.
C
999   DGSUM=0.0D0
      N=1
      M=1
      DO 880 I=1,NA
      KEYI=KEY(I)
      N1=N
      IF(NP(KEYI).NE.0) N1=N+1
      IF(ND(KEYI).NE.0) N1=N+2
      DO 881 J=N,N1
      DG1=0.0D0
      DG2=0.0D0
      L=1
      DO 882 K=1,NA
      KEYK=KEY(K)
      UB=(X(K)-REAL(ZS(KEYK)))*MAD(J,L)
      UC=X(K)*MAD(J,L)
      L=L+1
      IF(NP(KEYK).EQ.0) GO TO 883
      UB=UB+(Y(K)-REAL(ZP(KEYK)))*MAD(J,L)
      UC=UC+Y(K)*MAD(J,L)
      L=L+1
      IF(ND(KEYK).EQ.0) GO TO 883
      UB=UB+(Z(K)-REAL(ZD(KEYK)))*MAD(J,L)
      UC=UC+Z(K)*MAD(J,L)
      L=L+1
883   IF(K.NE.I) DG1=DG1+UB
882   DG2=DG2+UC
      J1=M+2*(J-N)
      DO 884 L=M,J1
      H(L,L)=W(L)+DG1
      IF(L5) HDG(L)=H(L,L)+QON*DG2
      W(L)=DG2
884   DGSUM=DGSUM+DG2
881   M=J1+1
880   N=N1+1
      CNST=CON
      DO 885 I=2,NDIM
      J1=I-1
      DO 885 J=1,J1
      IF(.NOT.L5) GO TO 886
      UB=(HDG(I)-HDG(J))/(HDG(I)+HDG(J))
      UB=UB*UB
      CNST=CON+UB/2.0D0+UB*UB*(0.5D0-CON)
886   UB=U(I,J)*(CNST*(H(I,I)+H(J,J))-QON*(0.5D0-CNST)*(W(I)+W(J)))
      H(I,J)=UB
885   H(J,I)=UB
      DGSUM=-(QON*DGSUM)/FLOAT(NDIM)
C
C  IF DOING LAST CYCLE CALCULATE ENERGY CORRECTIONS.
C
      IF(ICYCLE.NE.MAXCYC) GO TO 850
      N=1
      UB=0.0D0
      UC=0.0D0
      DO 887 I=1,NA
      KEYI=KEY(I)
      K=N
      A1=X(I)
      E(1)=REAL(ZS(KEYI))
      UB=UB-0.5D0*(A1*A1-A1)*MAD(N,N)
      IF(NP(KEYI).EQ.0) GO TO 888
      N=N+1
      A2=Y(I)
      E(3)=REAL(ZP(KEYI))
      UB=UB-0.5D0*(A2*A2-A2)*MAD(N,N)-A1*A2*MAD(N-1,N)
      IF(ND(KEYI).EQ.0) GO TO 888
      N=N+1
      A3=Z(I)
      E(5)=REAL(ZD(KEYI))
      UB=UB-0.5D0*(A3*A3-A3)*MAD(N,N)-A1*A3*MAD(N-2,N)-A2*A3*MAD(N-1,N)
888   M=1
      I1=I-1
      IF(I1.EQ.0) GO TO 887
      DO 889 J=1,I1
      KEYJ=KEY(J)
      L=M
      F(2)=X(J)
      E(2)=REAL(ZS(KEYJ))
      IF(NP(KEYJ).EQ.0) GO TO 890
      M=M+1
      F(4)=Y(J)
      E(4)=REAL(ZP(KEYJ))
      IF(ND(KEYJ).EQ.0) GO TO 890
      M=M+1
      F(6)=Z(J)
      E(6)=REAL(ZD(KEYJ))
890   DO 891 IJ=K,N
      DO 891 JK=L,M
      N1=2*(IJ-K)+1
      M1=2*(JK-L)+2
891   UC=UC-(F(N1)*F(M1)-E(N1)*E(M1))*MAD(IJ,JK)
889   M=M+1
887   N=N+1
      A1=UB
      A2=UC
      A3=UB+UC
C
C  SAVE MADELUNG TERMS FOR USE IN SUBROUTINE OUTPUT IF DOING
C  THE LAST CYCLE.
C
      K=1
      L=1
      DO 792 I=1,NA
      KEYI=KEY(I)
      MAD(K,K)=W(L)
      K=K+1
      L=L+1
      IF(NP(KEYI).EQ.0) GO TO 792
      MAD(K,K)=W(L)
      K=K+1
      L=L+3
      IF(ND(KEYI).EQ.0) GO TO 792
      MAD(K,K)=W(L)
      K=K+1
      L=L+5
792   CONTINUE
C
C  PRINT OUT RESULTS TO SHOW PROGRESS OF ITERATION PROCEDURE.
C
850   IF(PRTCYC.GT.0) WRITE(6,793) ICYCLE
793   FORMAT('CYCLE NO.',I3,'0')
      IF(PRTCYC.LT.0) WRITE(6,794)
794   FORMAT('CONVERGENCE REACHED - FINAL CYCLE FOLLOWS',///)
      PRTCYC=IABS(PRTCYC)
      J=NLF/10
      K=NLF-10*J
      WRITE(6,795) SUM,LAMBDA,J,ISIGNF,DENOM,ADJUST,K
795   FORMAT('+',T25,'ENERGY =',F15.8,T52,'LAMBDA =',F8.5,T78,'ATOM =',
     1 I3,T92,'SIGN =',I3,/,T25,'DENOM =',D16.8,T52,'ADJUST =',D14.7,
     2 T78,'NL   =',I3)
C
C  PRINT OUT ATOMIC CHARGES, ORBITAL OCCUPATIONS, AND CORRECTED
C  H(I,I)'S IF PRINTX IS TRUE.
C
      IF(.NOT.PRINTX) GO TO 500
      WRITE(6,600)
600   FORMAT(////,T16,'ATOM',T30,'NET CHG.-DAMPED',T60,'SUMMED ORBITAL O
     1CCUPATIONS-DAMPED'/T65,'S',T75,'P',T85,'D'/)
      DO 650 I=1,NA
      KEYI=KEY(I)
      WRITE(6,601) X(I)
601   FORMAT(T60,F10.5)
      IF(NP(KEYI).EQ.0) GO TO 625
      WRITE(6,602) Y(I)
602   FORMAT('+',T70,F10.5)
      IF(ND(KEYI).EQ.0) GO TO 625
      WRITE(6,603) Z(I)
603   FORMAT('+',T80,F10.5)
625   UB=C(I)
*650   WRITE(6,604) SYMBOL(KEYI),I,UB
650   WRITE(6,604) SYMBH(I),I,UB
604   FORMAT('+',T15,A2,I3,T30,F10.5)
      WRITE(6,809)
809   FORMAT(///,T16,'ATOM',T70,'CORRECTED H(I,I)''S',/T40,'S',T50,
     1 'X',T60,'Y',T70,'Z',T80,'X2-Y2',T90,'Z2',T100,'XY',T110,'XZ',
     2 T120,'YZ'/)
      J=1
      DO 810 I=1,NATM
      KEYI=KEY(I)
      N=J
      IF(NP(KEYI).NE.0) N=J+3
      IF(ND(KEYI).NE.0) N=J+8
815   WRITE(6,816) SYMBOL(KEYI),I,(H(K,K),K=J,N)
816   FORMAT(T15,A2,I3,T35,9F10.5)
810   J=N+1
      IF(METH.GE.3.AND.ABS(QON).GT.0.0001D0) WRITE(6,870) DGSUM
870   FORMAT(///,T15,'AVERAGE SHIFT OF MO''S DUE TO NON-ZERO TOTAL CHARG
     1E =',F12.8,' EV.')
      IF(METH.GE.3.AND.ICYCLE.EQ.MAXCYC) WRITE(6,808) A1,A2,A3
808   FORMAT(///,T15,'ENERGY CORRECTIONS0',//,T20,'ONE-CENTER',T40,
     1F16.8,' EV.',//,T20,'TWO-CENTER',T40,F16.8,' EV.',//,T20,'TOTAL',
     2T40,F16.8,' EV.')
      WRITE(6,752)
752   FORMAT(///)
500   ICYCLE=ICYCLE+1
C
C  CHECK FOR CONVERGENCE ( IE. DENOM LESS THAN DELTAC ).
C
      IF(ICYCLE.GE.MAXCYC) GO TO 433
      IF(ICYCLE.LE.2) GO TO 433
      IF(DENOM.GE.DELTAC) GO TO 433
      ICYCLE=MAXCYC
      PRTCYC=-PRTCYC
433   RETURN
      END
      SUBROUTINE LITER(NH,NA,E,W,NCYC)
C
C  SUBROUTINE FOR CALCULATING Q*SENSE WHEN USING CHARGE ITERATION
C  OPTION ( METH = 1 ).
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION E(1),W(1)
      COMMON/ATOM/KEY(1000),SYMBOL(100),VELEC(40),NS(40),NP(40),ND(40),
     1 EXPS(40),EXPP(40),EXPD(40),EXPD2(40),C1(40),C2(40),COULS(40),
     2 COULP(40),COULD(40),X(1000),Y(1000),Z(1000)
      INTEGER*4 SYMBOL,KEY,VELEC
      COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,PRTCYC,
     1 ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)
      REAL*8 LAMPRI
      INTEGER*4 PRTCYC
      LOGICAL*4 PARTIT,PRINTX,ITABLE
      NCYC=ICYCLE+1-(ICYCLE/4)*4
      ICYCLE=ICYCLE+1
      DELTA=0.D0
      NATOM=NH+NA
      INDEX=1
      DO 60 I=1,NATOM
      INCR=1
      IF(I.GT.NH) GO TO 100
      CHG=1.0D0-E(INDEX)
      GO TO 150
100   KEYI=KEY(I-NH)
      CHG=E(INDEX)
      IF(ND(KEYI).EQ.0) GO TO 120
      INCR=9
      CHG=CHG+E(INDEX+4)+E(INDEX+5)+E(INDEX+6)+E(INDEX+7)+E(INDEX+8)
      GO TO 130
120   IF(NP(KEYI).EQ.0) GO TO 140
      INCR=4
130   CHG=CHG+E(INDEX+1)+E(INDEX+2)+E(INDEX+3)
140   CHG=VELEC(KEYI)-CHG
150   INDEX=INDEX+INCR
      GO TO (10,20,30,40),NCYC
10    CHG=0.25D0*(CHG+Y(I)+Z(I)+W(I))
      X(I)=CHG
      Y(I)=CHG
      Z(I)=CHG
      W(I)=CHG
      GO TO 60
20    DELTA=DELTA+ABS(CHG-X(I))
      Y(I)=CHG
      GO TO 60
30    DELTA=DELTA+ABS(CHG-Y(I))
      Z(I)=CHG
      GO TO 60
40    DELTA=DELTA+ABS(CHG-Z(I))
      W(I)=CHG
60    E(I)=CHG*SENSE
      IF(DELTAC.LT.DELTA.OR.NCYC.EQ.1) RETURN
      ICYCLE=15000
      RETURN
      END
      SUBROUTINE MATRIX(N,NS1,NS2,NS3,NS4,NS5,NS6,NS7,NS8,NS9,NS10,
     X  NS11,NS12,NS13,ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9,ND10,
     X  ND11,ND12,ND13)
      IMPLICIT REAL*8(A-H,O-Z)

      Integer  Iaddress

C FORTRAN90 stuff ...
      double precision, allocatable::a(:)

C      DIMENSION A(4500000)
C      NSIZE=4500000
C
C    SUBROUTINE TO ALLOCATE STORAGE FOR MATRICES.
C
C      CALL ERRSET(74,0,0,0,1,0)
      I1=1
      I2=I1+NS1
      I3=I2+NS2
      I4=I3+NS3
      I5=I4+NS4
      I6=I5+NS5
      I7=I6+NS6
      I8=I7+NS7
      I9=I8+NS8
      I10=I9+NS9
      I11=I10+NS10
      I12=I11+NS11
      I13=I12+NS12
      IFF=I13+NS13
C      IF(IFF.GT.NSIZE) GO TO 200

      Write(6,*) '***  Dynamic memory allocation ***'
      Write(6,*) 'Need now memory (in kB): ',(IFF*8)/1000
      Write(6,*) ' '

C   If you ever use this program on a 64 bit machine check if the '8'
C   needs to be changed
C      CALL FMALLOC(Iaddress , IFF * 8)

      allocate(a(iff))
       
      CALL BANGIT(a,
     X  I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I13,
     X  ND1,ND2,ND3,ND4,ND5,ND6)

C      CALL BANGIT(%val(Iaddress),
C     X  I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I13,
C     X  ND1,ND2,ND3,ND4,ND5,ND6)

C      CALL MOVLAP (A(I1),A(I2),A(I3),A(I4),A(I5),A(I6),A(I7),
C     X  A(I8),A(I9),A(I10),A(I11),A(I12),A(I13),ND1,ND2,ND3,ND4,
C     X  ND5,ND6)
C      CALL HUCKEL (A(I1),A(I2),A(I3),A(I4),A(I5),A(I6),A(I7),
C     X  A(I8),A(I9),A(I10),A(I11),A(I12),A(I13),ND1,ND2,ND3,ND4,
C     X  ND5,ND6)
      RETURN
C 200  WRITE(6,1001) IFF,NSIZE
C 1001 FORMAT('*** INSUFFICIENT SPACE FOR MATRICES'/'PROGRAM REQUESTS '
C     1,I7,' DOUBLE WORDS BUT ONLY ',I7,' ARE AVAILABLE'/'RECOMPILE WITH
C     2LARGER DIMENSION FOR A AND INCREASED VALUE OF NSIZE')
C      RETURN
      END
C**   Stuff added
      SUBROUTINE BANGIT(A,I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,
     X  I11,I12,I13,ND1,ND2,ND3,ND4,ND5,ND6)

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(*)

      CALL MOVLAP (A(I1),A(I2),A(I3),A(I4),A(I5),A(I6),A(I7),
     X  A(I8),A(I9),A(I10),A(I11),A(I12),A(I13),ND1,ND2,ND3,ND4,
     X  ND5,ND6)
      CALL HUCKEL (A(I1),A(I2),A(I3),A(I4),A(I5),A(I6),A(I7),
     X  A(I8),A(I9),A(I10),A(I11),A(I12),A(I13),ND1,ND2,ND3,ND4,
     X  ND5,ND6)
      RETURN
      END
C** End of added stuff (LUL 1992-07-01)
      SUBROUTINE ABFNS(A,B)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(20),B(20)
      COMMON/LOCLAP/SK1,SK2,RR,L1,L2,M,N1,N2,MAXCAL
C
C    SUBROUTINE FOR CALCULATING A AND B FUNCTIONS FOR USE IN LOVLAP.
C
      J=MAXCAL+1
      RHO1=0.5D0*(SK1+SK2)*RR
      RHO2=0.5D0*(SK1-SK2)*RR
      IF(ABS(RHO1).GT.165.D0) GO TO 100
      IF(ABS(RHO2).GT.165.D0) GO TO 100
      C=EXP(-RHO1)
      A(1)=C/RHO1
      DO 15 I=2,J
 15   A(I)=(FLOAT(I-1)*A(I-1)+C)/RHO1
      IX=J
      IR=ABS(2.*RHO2)
      IS=MIN0(IR+1,19)
      IF(RHO2) 25,35,25
 25   D=EXP(RHO2)
      H=1.D0/D
C
C    USE THE DSINH ROUTINE INSTEAD OF SUMMING THE INFINITE SERIES.
C
      R=2.D0*SINH(RHO2)
C
C    AS MANY SUCCESSIVE B-FUNCTIONS ARE GENERATED FROM B-0 BY THE
C    RECURRENCE FORMULAE AS ACCURACY WILL PERMIT.
C
      B(1)=R/RHO2
      DO 51 I=2,IX,IS
      IF(IR.EQ.0) GO TO 40
      IL=IS-1
C
C    MODIFICATION TO AVOID EXCEEDING STORAGE LIMITS.
C    D. WALLACE 04/14/71
C
      DO 31 K=I,IX
      IF((-1)**K) 29,29,30
 29   B(K)=(R+FLOAT(K-1)*B(K-1))/RHO2
      GO TO 31
 30   B(K)=-(D+H-FLOAT(K-1)*B(K-1))/RHO2
 31   CONTINUE
 40   IN=I+IS-1
C
C    AFTER THE RECURRENCE FORMULAE HAVE BEEN APPLIED AN APPROPRIATE
C    NUMBER OF TIMES THE NEXT B-FUNCTION IS OBTAINED BY SUMMATION
C    OF THE INFINITE SERIES.
C
      IF(IN-IX) 39,39,38
 39   IF((-1)**IN) 44,44,42
 42   TR=RHO2
 105  B(IN)=-2.*TR/FLOAT(IN+1)
      DO 43 J=1,500
      TR=TR*RHO2**2/FLOAT((2*J)*(2*J+1))
      IF(ABS(TR/B(IN))-1.0D-7 ) 51,51,43
 43   B(IN)=B(IN)-2.*TR/FLOAT(IN+1+2*J)
      GO TO 51
 44   TR=1.
 107  B(IN)=2.*TR/FLOAT(IN)
      DO 46 J=1,500
      TR=TR*RHO2**2/FLOAT((2*J)*(2*J-1))
      IF(ABS(TR/B(IN))-1.0D-7 ) 51,51,46
 46   B(IN)=B(IN)+2.*TR/FLOAT(IN+2*J)
 51   CONTINUE
C
C    IF THE ARGUMENT IS ZERO A SEPARATE FORMULA MUST BE USED.
C
      GO TO 38
 35   DO 36 I=1,IX,2
      B(I)=2.D0/FLOAT(I)
 36   B(I+1)=0.D0
 38   RETURN
 100  DO 101 I=1,20
      A(I)=0.D0
 101  B(I)=0.D0
      GO TO 38
      END
      SUBROUTINE LOVLAP (STRAD,A,B)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(20),B(20)
      DIMENSION FACT(25)
      COMMON/LOCLAP/SK1,SK2,R,L1,L2,M1,N1,N2,MAX
      DIMENSION BINCOE(7,7)
      DATA BINCOE/7*1.D0,   0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     1 2*0.D0,1.D0,3.D0,6.D0,10.D0,15.D0,   3*0.D0,1.D0,4.D0,10.D0,
     2 20.D0,   4*0.D0,1.D0,5.D0,15.D0,   5*0.D0,1.D0,6.D0,   6*0.D0,
     3 1.D0/
      LOGICAL*4 JGO
      DATA JGO/.FALSE./
C
C    SUBROUTINE TO CALCULATE OVERLAP INTEGRALS IN A LOCAL
C    COORDINATE SYSTEM.
C
C    INTEGRALS ARE CALCULATED BY TRANSFORMATION TO ELLIPSOIDAL
C    COORDINATES AND THEREBY EXPRESSED IN TERMS OF C-FUNCTIONS.
C    SEE J.C.P.,24,201. ORIGINALLY WRITTEN BY R.M.STEVENS.
C
C
C    GENERATE FACTORIALS ONLY ONCE.
C
      IF(JGO) GO TO 10
      JGO=.TRUE.
      FACT(1)=1.D0
      DO 5 I=2,25
 5    FACT(I)=FACT(I-1)*FLOAT(I-1)
 10   CONTINUE
      M2=M1
      STRAD=0.D0
      RHOA=R*SK1
      RHOB=R*SK2
      TERMA=0.5D0**(L1+L2+1) * SQRT(FLOAT((L1+L1+1)*(L2+L2+1))*
     1 FACT(L1-M1+1)*FACT(L2-M1+1)/(FACT(N1+N1+1)*FACT(N2+N2+1)*
     2 FACT(L1+M1+1)*FACT(L2+M1+1))*RHOA**(N1+N1+1)*RHOB**(N2+N2+1))
      JEND=1+((L1-M1)/2)
      KEND=1+((L2-M2)/2)
      IEB=M1+1
      DO 50 J=1,JEND
      JU=J-1
      IAB=N1-L1+JU+JU+1
      ICB=L1-M1-JU-JU+1
      CON1=FACT(L1+L1-JU-JU+1)/(FACT(L1-M1-JU-JU+1)*FACT(JU+1)*
     1 FACT(L1-JU+1))
      DO 50 K=1,KEND
      KU=K-1
      CON12=CON1*FACT(L2+L2-KU-KU+1)/(FACT(L2-M2-KU-KU+1)*FACT(KU+1)*
     1 FACT(L2-KU+1))
      IEV=JU+KU+L2
      IF(2*(IEV/2).NE.IEV) CON12=-CON12
      IBB=N2-L2+KU+KU+1
      IDB=L2-M2-KU-KU+1
      VALUE=0.D0
      DO 90 I6=1,IEB
      DO 90 I5=1,IEB
      VALUE1=BINCOE(IEB,I6)*BINCOE(IEB,I5)
      IEV=I5+I6
      IF(2*(IEV/2).NE.IEV) VALUE1=-VALUE1
      DO 90 I4=1,IDB
      VALUE1=-VALUE1
      VALUE2=BINCOE(IDB,I4)*VALUE1
      DO 90 I3=1,ICB
      VALUE3=BINCOE(ICB,I3)*VALUE2
      DO 90 I2=1,IBB
      VALUE3=-VALUE3
      VALUE4=BINCOE(IBB,I2)*VALUE3
      DO 90 I1=1,IAB
      TERM=VALUE4*BINCOE(IAB,I1)
      IR=I1+I2+IEB+IEB-I6-I6-I3+IDB-I4+ICB-1
      IP=IAB-I1+IBB-I2+IEB+IEB-I5-I5+ICB-I3+IDB-I4+1
 90   VALUE=VALUE+A(IP)*B(IR)*TERM
 50   STRAD=STRAD+VALUE*CON12
      STRAD=STRAD*TERMA
      RETURN
      END
      SUBROUTINE GRMSCH(U,C,NDIM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(NDIM,NDIM),C(NDIM)
C
C    SUBROUTINE TO CARRY OUT A GRAM SCHMIDT ORTHOGONALIZATION ON A
C    SET OF VECTORS X(I). THEY ARE CONVERTED INTO A SET OF ORTHONORMAL
C    VECTORS Y(J).  THE UNITARY TRANSFORMATION IS DEVELOPED.
C
C    U   CONTAINS THE OVERLAP MATRIX IN THE UPPER RIGHT TRIANGULAR
C        PART. THE LOWER TRIANGULAR PART IS NOT USED. THE UNITARY
C        TRANSFORMATION IS RETURNED IN THE UPPER RIGHT TRIANGLE,
C        INCLUDING THE DIAGONAL.
C    C   IS A WORK VECTOR OF LENGTH NDIM.
C
C    NDIM IS THE DIMENSION OF C AND U (IE. C(NDIM), U(NDIM,NDIM)).
C
C    ACTUALLY THE FIRST COLUMN OF U MAY BE USED AS THE WORK AREA
C    PROVIDED THAT THE ELEMENT U(1,1) IS SET EQUAL TO 1.D0 AFTER
C    RETURN IS MADE TO THE CALLING PROGRAM.
C
      U(1,1)=1.D0
      DO 100 I=2,NDIM
      I1=I-1
      XNORM=1.D0
      DO 40 J=1,I1
      SUM=0.D0
      DO 30 K=1,J
      SUM=SUM+U(K,J)*U(K,I)
 30   CONTINUE
      C(J)=SUM
      XNORM=XNORM-SUM**2
 40   CONTINUE
      XNORM=SQRT(1.D0/XNORM)
      DO 70 J=1,I1
      SUM=0.D0
      DO 60 K=J,I1
      SUM=SUM+C(K)*U(J,K)
 60   CONTINUE
      U(J,I)=-SUM*XNORM
 70   CONTINUE
      U(I,I)=XNORM
 100  CONTINUE
      RETURN
      END
      SUBROUTINE TRNFRM(S,H,C,COUL0,NDIM,SP,IEXIT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S(NDIM,NDIM),H(NDIM,NDIM),C(NDIM),COUL0(NDIM),SP(NDIM)
      LOGICAL*4 ONEMAT
C
C    SUBROUTINE FOR CARRYING OUT A CHANGE IN BASIS SET ON A MATRIX H
C    BY MEANS OF A MATRIX U.
C
C                C=U'*H*U
C
C    H IS ASSUMED REAL SYMMETRIC AND ONLY THE LOWER LEFT TRIANGLE
C    IS REFERENCED. U IS ASSUMED TO BE ENTIRELY CONTAINED IN ITS
C    UPPER RIGHT TRIANGLE AND ONLY THIS PART IS REFERENCED. THE
C    RESULTING REAL SYMMETRIC MATRIX, C, IS STORED IN PACKED FORM.
C
C    INDEX1 OF LOOP1 POINTS TO U(1,J), J=1,NDIM
C    INDEX2 OF LOOP2 POINTS TO H(I,1), I=1,J
C    INDEX3 OF LOOP3 POINTS TO U(K,J), K=1,I
C    INDEX4 OF LOOP4 POINTS TO U(K,J), K=I+1,J
C    INDEX5 OF LOOP5 POINTS TO U(1,I), I=1,J
C    INDEX6 OF LOOP6 POINTS TO U(K,I), K=1,I
C
 
      ONEMAT=IEXIT.EQ.2.OR.IEXIT.EQ.3
      ISUB=1
      DO 100 J=1,NDIM
      DO 40 I=1,J
      SUM=0.D0
      ILIM=I
      IF(.NOT.ONEMAT) GO TO 10
      IF(I.LT.J) GO TO 10
      ILIM=I-1
C
C    IN THE ONE MATRIX CASE, THE S AND H MATRICES OVERLAP ALONG THE
C    DIAGONAL.  THE DIAGONAL OF S IS THEN PUT INTO SP.
C
      SUM=SP(I)*H(I,I)
      IF(I.EQ.1) GO TO 25
 10   DO 20 K=1,ILIM
 20   SUM=SUM+S(K,J)*H(I,K)
 25   IF(I.EQ.J) GO TO 40
      I1=I+1
      ILIM=J
      IF(.NOT.ONEMAT) GO TO 27
      ILIM=J-1
      SUM=SUM+SP(J)*H(J,I)
      IF(ILIM.LT.I1) GO TO 40
 27   DO 30 K=I1,ILIM
 30   SUM=SUM+S(K,J)*H(K,I)
 40   COUL0(I)=SUM
      DO 60 I=1,J
      SUM=0.D0
      ILIM=I
      IF(.NOT.ONEMAT) GO TO 45
      ILIM=I-1
      SUM=COUL0(I)*SP(I)
      IF(ILIM.EQ.0) GO TO 55
 45   DO 50 K=1,ILIM
 50   SUM=SUM+COUL0(K)*S(K,I)
 55   C(ISUB)=SUM
 60   ISUB=ISUB+1
 100  CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DSUM(B,A,IP1,LIMIT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION B(1),A(1)
C
C    FUNCTION FOR USE IN GIVENS DIAGNOLIZATION PACKAGE.
C
      JJ=1
      DSUM=0.D0
      DO 180 II=IP1,LIMIT
      DSUM=DSUM+B(II+1)*A(JJ)
 180  JJ=JJ+II
      RETURN
      END
      SUBROUTINE ROTATE(V,C,S,NJX,JTOP)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION V(1)
C
C    FUNCTION FOR USE IN GIVENS DIAGNOLIZATION PACKAGE.
C
      JLIM=JTOP+1
      DO 10 J=1,JLIM,NJX
      TA=V(J)
      TB=V(J+1)
      V(J)=TA*C+TB*S
 10   V(J+1)=TB*C-TA*S
      RETURN
      END
      DOUBLE PRECISION FUNCTION DOT(A,B)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(1),B(1)
      COMMON /VECTOR/ FACTOR,LIMIT
C
C    FUNCTION FOR USE IN GIVENS DIAGNOLIZATION PACKAGE.
C
      ITOP=LIMIT+1
      DOT=0.D0
      DO 10 I=1,ITOP
 10   DOT=DOT+A(I)*B(I)
      RETURN
      END
      SUBROUTINE VECSUM(A,B)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(1),B(1)
      COMMON /VECTOR/ FACTOR,LIMIT
C
C    FUNCTION FOR USE IN GIVENS DIAGNOLIZATION PACKAGE.
C
      ITOP=LIMIT+1
      DO 10 I=1,ITOP
 10   A(I)=A(I)+FACTOR*B(I)
      RETURN
      END
      SUBROUTINE REDUCEM(U,NDIM,IDUMMY1,IDUMMY2) 
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(NDIM,NDIM)
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,
     1 IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT
      LOGICAL*4 L1,L2,L3,L4,L5,ONEMAT,ITERAT
      COMMON/ATOM/KEY(1000),SYMBOL(100),VELEC(40),NS(40),NP(40),
     1 ND(40),
     1 EXPS(40),EXPP(40),EXPD(40),EXPD2(40),C1(40),C2(40),
     2 COULS(40),COULP(40),COULD(40),X(1000),Y(1000),Z(1000)
      INTEGER*4 SYMBOL,KEY,VELEC
C
C    SUBROUTINE TO CALCULATE REDUCED MATRIX.
C
      CALL REDCHM(U,NDIM)
      ISUB=NH+1
      IAO=ISUB
      DO 100 I=1,NA
      KEYI=KEY(I)
      IF(ND(KEYI)) 10,20,10
 10   IATOP=IAO+8
      GO TO 50
 20   IF(NP(KEYI)) 30,40,30
 30   IATOP=IAO+3
      GO TO 50
 40   IATOP=IAO
 50   DO 70 J=1,NATM
      SUM=0.D0
      DO 60 K=IAO,IATOP
      SUM=SUM+U(J,K)
 60   CONTINUE
      U(J,ISUB)=SUM
 70   CONTINUE
      ISUB=ISUB+1
 100  IAO=IATOP+1
      RETURN
      END
      SUBROUTINE FULCHM(E,U,C,H,N5)
      REAL*8 E(N5),C(N5),U(N5,N5),H(N5,N5)
      LOGICAL*4 JGO
C
C    SUBROUTINE TO CALCULATE COMPLETE CHARGE MATRIX.
C
      KJ=1
      DO 31 I=1,N5
      JGO=.FALSE.
      IJ=KJ
      DO 32 K=1,N5
      IF(JGO) GO TO 34
 33   IF(I.NE.K) GO TO 35
      JGO=.TRUE.
      E(K)=1.0D0
      GO TO 32
 35   E(K)=C(IJ)
      IJ=IJ+1
      GO TO 32
 34   IJ=IJ+K-2
      E(K)=C(IJ)
 32   CONTINUE
      KJ=KJ+I-1
      DO 31 J=1,N5
      UB=0.0D0
      DO 36 K=1,N5
 36   UB=UB+H(K,J)*E(K)
 31   U(I,J)=2.0D0*H(I,J)*UB
      RETURN
      END
      SUBROUTINE REDCHM(U,N5)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(N5,N5)
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,
     1 IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT
      LOGICAL*4 L1,L2,L3,L4,L5,ONEMAT,ITERAT
      COMMON/ATOM/KEY(1000),SYMBOL(100),VELEC(40),NS(40),NP(40),
     1 ND(40),
     1 EXPS(40),EXPP(40),EXPD(40),EXPD2(40),C1(40),C2(40),COULS(40),
     2 COULP(40),COULD(40),X(1000),Y(1000),Z(1000)
      INTEGER*4 SYMBOL,KEY,VELEC
C
C    SUBROUTINE TO CALCULATE REDUCED CHARGE MATRIX.
C
      IF(NA.EQ.0) RETURN
      NH1=NH+1
      IAO=NH1
      ISUB=NH1
      DO 128 I=1,NA
      KEYI=KEY(I)
      IF(NP(KEYI)) 122,121,122
 121  IATOP=IAO
      GO TO 125
 122  IF(ND(KEYI)) 124,123,124
 123  IATOP=IAO+3
      GO TO 125
 124  IATOP=IAO+8
 125  DO 127 J=1,N5
      SUM=0.D0
      DO 126 K=IAO,IATOP
 126  SUM=SUM+U(K,J)
 127  U(ISUB,J)=SUM
      ISUB=ISUB+1
 128  IAO=IATOP+1
      RETURN
      END
      SUBROUTINE CORECT(A)
      RETURN
      END
      SUBROUTINE Jacobi(A,R,N,MV,ROOT)
      Implicit Real*8 (A-H,O-Z)
      Parameter (Mxinv=1000)
      DIMENSION A(1),ROOT(1),R(1)
      Dimension INV(Mxinv)
C....   GENERATE IDENTITY MATRIX
*
      If(N.gt.Mxinv) Then
       Write(6,'(///A,2x,i8,A)')
     1 ' Requesting N= ',N,' in JACOBI'
       Write(6,'(A)')
     2 ' Please change MXINV to match N'
       STOP
      End if
*
      Do 2345 I=1,N+1
2345   Inv(i)=i*(i-1)/2
*
      IF(MV-1) 10,25,10
10    IQ=-N      
      DO 20 J=1,N
      IQ=IQ+N
      DO 20 I=1,N
      IJ=IQ+I
      R(IJ)=0.0
      IF(I-J) 20,15,20
15    R(IJ)=1.
20    CONTINUE
C....  COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANRMX)
25    ANORM=0.0
      DO 35 I=1,N
      DO 35 J=I,N
      IF(I-J) 30,35,30
30    IA=I+INV(J)
      ANORM=ANORM+A(IA)*A(IA)
35    CONTINUE
      IF(ANORM) 165,165,40
40    ANORM=1.414*SQRT(ANORM)
      XN=N
      ANRMX=ANORM*1.E-08/XN
C....  INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR
      IND=0
      THR=ANORM
45    THR=THR/XN
50    L=1
55    M=L+1
C....  COMPUTE SIN AND COS
60    MQ=INV(M)
      LQ=INV(L)
      LM=L+MQ
62    IF(ABS(A(LM))-THR) 130,65,65
65    IND=1
      LL=L+LQ
      MM=M+MQ
      X=0.5*(A(LL)-A(MM))
68    Y=-A(LM)/SQRT(A(LM)*A(LM)+X*X)
      IF(X) 70,75,75
70    Y=-Y
75    SINX=Y/SQRT(2.*(1.+(SQRT(ABS(1.-Y*Y)))))
      SINX2=SINX*SINX
78    COSX=SQRT(ABS(1.-SINX2))
      COSX2=COSX*COSX
      SINCS=SINX*COSX
C....   ROTATE L AND M COLUMNS
      ILQ=N*(L-1)
      IMQ=N*(M-1)
      DO 125 I=1,N
      IQ=INV(I)
      IF(I-L) 80,115,80
80    IF(I-M) 85,115,90
85    IM=I+MQ
      GO TO 95
90    IM=M+IQ
95    IF(I-L) 100,105,105
100   IL=I+LQ
      GO TO 110
105   IL=L+IQ
110   X=A(IL)*COSX-A(IM)*SINX
      A(IM)=A(IL)*SINX+A(IM)*COSX
      A(IL)=X
115   IF(MV-1) 120,125,120
120   ILR=ILQ+I
      IMR=IMQ+I
      X=R(ILR)*COSX-R(IMR)*SINX
      R(IMR)=R(ILR)*SINX+R(IMR)*COSX
      R(ILR)=X
125   CONTINUE
      X=2.*A(LM)*SINCS
      Y=A(LL)*COSX2+A(MM)*SINX2-X
      X=A(LL)*SINX2+A(MM)*COSX2+X
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
      A(LL)=Y
      A(MM)=X
C...   TEST FOR COMPLETION
C...   TEST FOR M = LAST COLUMN
130   IF(M-N) 135,140,135
135   M=M+1
      GO TO 60
C....  TEST FOR L = SECOND FROM LAST COLUMN
140   IF(L-(N-1)) 145,150,145
145   L=L+1
      GO TO 55
150   IF(IND-1) 160,155,160
155   IND=0
      GO TO 50
C....  COMPARE THRESHOLD WITH FINAL NORM
160   IF(THR-ANRMX) 165,165,45
C....  SORT EIGENVALUES AND EIGENVECTORS
165   IQ=-N
      DO 185 I=1,N
      IQ=IQ+N
      LL=INV(I+1)
      JQ=N*(I-2)
      DO 185 J=I,N
      JQ=JQ+N
      MM=INV(J+1)
      IF(A(LL)-A(MM)) 170,185,185
170   X=A(LL)
      A(LL)=A(MM)
      A(MM)=X
      IF(MV-1) 175,185,175
175   DO 180 K=1,N
      ILR=IQ+K   
      IMR=JQ+K
      X=R(ILR)
      R(ILR)=R(IMR)
180   R(IMR)=X
185   CONTINUE
      DO 190 I=1,N
      II=INV(I+1)
190   ROOT(I)=A(II)
      RETURN
      END
      Subroutine PreRSP(C,H,Ndim,Coul0)
      Implicit Real*8 (A-h,O-z)
      Parameter (NMAX=5000)
      Dimension C(Ndim),H(Ndim,Ndim),Coul0(Ndim)
      Dimension E(NMAX),EE(NMAX)
      Data Irev/1/
*
      If(Ndim.gt.Nmax) Then
       Write(6,*) ' NDIM > NMAX in PreRSP '
       Stop
      EndIf
*
      Move=0
*
*      Do 5 j=1,Ndim
*       Do 5 i=1,j
*        Move=Move+1 
*5     H(j,i)=C(Move)
*  
*      Call QL77(Ndim,H,Coul0,E,Ndim,Irev)
*
*      Do 5 i=1,Ndim
*       Do 5 j=1,i
*        Move=Move+1 
*5     H(j,i)=C(Move)
*
*      Do 10 j=1,Ndim
*       Do 10 i=1,Ndim
*10      AA=H(j,i)
*      STOP 
*  
*      Call QL77(Ndim,H,Coul0,E,Ndim,Irev)
       Call RSP(C,NDIM,IREV,Coul0,H,E,EE)
*
       Num=Ndim/2
        Do 50 i=1,Num
        Ndimi=Ndim-i+1
         Switch=Coul0(i)
         Coul0(i)=Coul0(Ndimi)
          Coul0(Ndimi)=Switch
        Do 60 j=1,Ndim
         Switch=H(j,i)
         H(j,i)=H(j,Ndimi)
          H(j,Ndimi)=Switch
60      Continue
50      Continue
*
*       IF(MND.ne.0) Then
*       Write(6,*) ' MND not .eq. to 0 '
*       End if
*
      Return
      End
C
C  SAMPLE INPUT FOR QCPE 517 (REMOVE THESE THREE COMMENT LINES)
C
*ETHYLENE
*  4  2        5
* 0.92665        1.205          0.0
*-0.92665        1.205          0.0
* 0.92665       -1.205          0.0
*-0.92665       -1.205          0.0
* 0.0            0.67           0.0
* 0.0           -0.67           0.0
* C C

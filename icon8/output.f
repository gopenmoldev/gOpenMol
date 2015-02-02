      SUBROUTINE OUTPUT(H,U,MAD,C,E,W,IOCC,HDG,NDIM,NTYPE,NC,NHDG)
C
C  SUBROUTINE TO ANALYSE AND PRINT OUT RESULTS.
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMPRI
      INTEGER*4 PRTCYC
      LOGICAL*4 PARTIT,PRINTX,ITABLE
      INTEGER*4 SYMB,HYDROG
      REAL*8 MAD
      REAL*4 IOCC
      LOGICAL*4 PRT,PUN
      INTEGER*4 IOVPOP,IENRGY
      LOGICAL*4 L1,L2,L3,L4,L5,ONEMAT,ITERAT
      INTEGER*4 SYMBOL,KEY,VELEC,SymbH
      DIMENSION H(NDIM,NDIM),U(NDIM,NDIM),MAD(NTYPE,NTYPE),C(NC),
     1 E(NDIM),W(NDIM),IOCC(NDIM),HDG(NHDG)
      Common/LULX/Ipointer(1000),GrossP(1000),Screen(1000)
      Common/LUL/IhelpV(100),SymbH(1000),XX(1000),YY(1000),ZZ(1000)
     1,MaxDim    
      COMMON/TITLE/AB(10)
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,
     1 IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT
      COMMON/OUT/PRT(20),PUN(20),IOVPOP(24),IENRGY(24)
      COMMON/ATOM/KEY(1000),SYMBOL(100),VELEC(40),NS(40),NP(40),ND(40),
     1 EXPS(40),EXPP(40),EXPD(40),EXPD2(40),C1(40),C2(40),COULS(40),
     2 COULP(40),COULD(40),X(1000),Y(1000),Z(1000)
      COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,PRTCYC,
     1 ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)
      DATA HYDROG/' H'/
      DO 5 I=1,NDIM
      NMIN=I
      IF(IOCC(I).GT.0.0001) GO TO 7
5     CONTINUE
7     IF(ITERAT) GO TO 38
C
C  CALCULATE AND PRINT OUT SUM OF ONE-ELECTRON ENERGIES.
C                                                             
10    SUM=0.0D0
      enconv = 1./27.21161
C
      open(22,FILE='energies.dat',STATUS="UNKNOWN")
      Write(22,'(i5,2x,i5)') NDIM,NATM
C
      DO 13 I=1,NDIM
      W(I)=DBLE(IOCC(I))
C
C Own addition to flush out the orbital energies and occupation
C
      Write(22,'(f16.8,2x,f16.8,2x,f10.4)') E(NDIM - i + 1),
     1     E(NDIM - i + 1)*enconv,
     1     DBLE(IOCC(NDIM - i + 1))
13    SUM=SUM+E(I)*W(I)
      IF(.NOT.PRT(9)) WRITE(6,2001) SUM,SUM*enconv
2001  FORMAT(//,5x,'SUM OF ONE-ELECTRON ENERGIES =',F16.8,' EV.',
     13x,F16.8,' AU. ',///)

      open(7,FILE="punch.dat",STATUS="UNKNOWN")
      IF(PUN(9)) WRITE(7,2002) SUM
2002  FORMAT(F20.8)
      IF(ONEMAT) GO TO 9999
C
C  PRINT OUT WAVE FUNCTIONS.
C
      IF(PRT(10)) GO TO 1003
      WRITE(6,1002)
1002  FORMAT('WAVE FUNCTIONS'/'MO''S IN COLUMNS, AO''S IN ROWS')
      CALL PEGLEG(H,NDIM,NDIM)
*
* Changed the following output format 
* LUL 1990
*
1003  IF(PUN(10)) Then
       write(7,*) '*** Wave function matrix *** '
        WRITE(7,2003) H
2003     FORMAT(8(1x,F9.6))
      Endif 
C
C  CALCULATE AND PRINT OUT DENSITY MATRIX.
C
500   IF(PRT(11).AND..NOT.PUN(11)) GO TO 38
      DO 300 I=1,NDIM
      DO 300 J=1,I
      U(I,J)=0.0D0
      DO 310 K=1,NDIM
310   U(I,J)=U(I,J)+H(I,K)*H(J,K)*W(K)
300   U(J,I)=U(I,J)
      IF(PRT(11)) GO TO 360
      WRITE(6,350)
350   FORMAT('DENSITY MATRIX')
      CALL PEGLEG(U,NDIM,NDIM)
360   IF(PUN(11)) WRITE(7,2003) U
C
C  CALCULATE ATOMIC ORBITAL OCCUPATIONS AND STORE IN E(I).
C  CALCULATE OVERLAP POPULATION MATRIX.
C
38    IJ=1
      DO 60 I=1,NDIM
      E(I)=0.0D0
      DO 60 J=1,I
      UB=0.0D0
      DO 41 K=NMIN,NDIM
41    UB=UB+H(I,K)*H(J,K)*DBLE(IOCC(K))
      UB=UB*0.5D0
      IF(I.EQ.J) GO TO 50
      UB=(UB+UB)*C(IJ)
      IJ=IJ+1
50    E(I)=E(I)+UB
      E(J)=E(J)+UB
      IF(ITERAT) GO TO 60
      UB=UB+UB
      U(I,J)=UB
      U(J,I)=UB
60    CONTINUE
C
C  IF DOING CHARGE ITERATION ( METH=1 ) CALL LITER.
C
      IF(.NOT.ITERAT) GO TO 80
      K=ICYCLE
      PRINTX=((ICYCLE/PRTCYC)*PRTCYC.EQ.ICYCLE)
      IF(ICYCLE.EQ.MAXCYC) PRINTX=.TRUE.
      CALL LITER(NH,NA,E,W,J)
      IF(ICYCLE.EQ.15000) PRINTX=.TRUE.
      IF(K.EQ.1) WRITE(6,355)
355   FORMAT('1     ATOMIC CHARGES',//)
      GO TO (81,82,83,84),J
81    WRITE(6,375) K,(X(I),I=1,NATM)
375   FORMAT(/,T3,'CYCLE NO.',I3,(T20,10F10.5))
      GO TO 9000
82    WRITE(6,375) K,(Y(I),I=1,NATM)
      GO TO 9000
83    WRITE(6,375) K,(Z(I),I=1,NATM)
      GO TO 9000
84    WRITE(6,375) K,(W(I),I=1,NATM)
9000  RETURN          
C
C  PRINT OUT OVERLAP POPULATION MATRIX.
C
80    IF(PRT(12)) GO TO 1009
      WRITE(6,1006) NELEC
1006  FORMAT('OVERLAP POPULATION MATRIX FOR',I4,' ELECTRONS')
      CALL PEGLEG(U,NDIM,NDIM)
1009  IF(PUN(12)) WRITE(7,2003) U
C
C  CALL REDUCE TO CALCULATE REDUCED OVERLAP MATRIX. PRINT OUT
C  REDUCED OVERLAP MATRIX.
C
      IF(PRT(13).AND..NOT.PUN(13)) GO TO 2005
      CALL REDUCEM(U,NDIM,NA,NH)
      DO 100 I=2,NATM
      K=I-1
      DO 100 J=1,K
100   U(J,I)=U(I,J)
      IF(PRT(13)) GO TO 1015
      WRITE(6,1007)
1007  FORMAT('REDUCED OVERLAP POPULATION MATRIX, ATOM BY ATOM')
      CALL PEGLEG(U,NATM,NDIM)
1015  IF(PUN(13)) WRITE(7,2003) ((U(I,J),I=1,NATM),J=1,NATM)
C
C  IF L3 IS TRUE, CALCULATE AND PRINT OUT OVERLAP POPULATION ANALYSIS,
C  ORBITAL BY ORBITAL, FOR EACH MOLECULAR ORBITAL SPECIFIED.
C
2005  IF(.NOT.L3) GO TO 21
      DO 600 N=1,23,2
      IF(IOVPOP(N).EQ.0) GO TO 25
      KMIN=IOVPOP(N)
      KMAX=IOVPOP(N+1)
      DO 600 K=KMIN,KMAX
      WRITE(6,2004) K,W(K)
2004  FORMAT(///'OVERLAP POPULATION MATRIX, ORBITAL BY ORBITAL, FOR MOL
     1ECULAR ORBITAL',I4,5X,'OCCUPATION IS',F7.4,' ELECTRONS')
      IF(IOCC(K).LT.0.0001) GO TO 600
      SUM=0.5D0*W(K)
      IJ=1
      DO 460 I=1,NDIM
      DO 460 J=1,I
      UB=H(I,K)*H(J,K)
      IF(I.EQ.J) GO TO 450
      UB=(UB+UB)*C(IJ)
      IJ=IJ+1
450   UB=(UB+UB)*SUM
      U(J,I)=UB
      U(I,J)=UB
460   CONTINUE
      CALL PEGLEG(U,NDIM,NDIM)
600   CONTINUE
25    WRITE(6,7003)
7003  FORMAT(///)
C
C  CALL FULCHM TO CALCULATE COMPLETE CHARGE MATRIX. PRINT OUT
C  COMPLETE CHARGE MATRIX.
C
21    L1=PRT(14).AND..NOT.PUN(14)
      L2=PRT(15).AND..NOT.PUN(15)
      IF(L1.AND.L2) GO TO 1020
      CALL FULCHM(W,U,C,H,NDIM)
      IF(PRT(14)) GO TO 2021
      WRITE(6,1008)
1008  FORMAT('COMPLETE CHARGE MATRIX FOR EACH MO, NORMALIZED TO TWO ELE
     1CTRONS REGARDLESS OF OCCUPATION')
      CALL PEGLEG(U,NDIM,NDIM)
2021  IF(PUN(14)) WRITE(7,2003) U
C
C  CALL REDCHM TO CALCULATE REDUCED CHARGE MATRIX. PRINT OUT
C  REDUCED CHARGE MATRIX.
C
      IF(L2) GO TO 1020
      CALL REDCHM(U,NDIM)
      IF(PRT(15)) GO TO 1022
      WRITE(6,1019)
1019  FORMAT('REDUCED CHARGE MATRIX, MO''S IN COLUMNS, ATOMS IN ROWS')
      CALL OUTMAT(U,NDIM,NATM,NDIM)
1022  IF(PUN(15)) WRITE(7,2003) ((U(I,J),I=1,NATM),J=1,NDIM)

C
C  PRINT OUT ATOMIC CHARGES AND ORBITAL OCCUPATIONS.
C
1020  IF(PRT(16)) GO TO 40
C
C     open file to flush out the input file to vss program
C
      open(unit=21 , file = 'vss.inp' , status = 'unknown',
     1     access = 'sequential' , err = 12345)    
      WRITE(6,1010)
1010  FORMAT(/'ATOM',T12,'NET CHG.',T35,'ATOMIC ORBITAL OCCUPATION FOR
     *GIVEN MO OCCUPATION'/T35,'S',T45,'X',T55,'Y',T65,'Z',T75,'X2-Y2',
     *T85,'Z2',T95,'XY',T105,'XZ',T115,'YZ'/)
      SYMB=HYDROG
      J=1
      DO 140 I=1,NATM
       Grossp(i) = 0.0
      IF(I.GT.NH) GO TO 120
      UB=1.0-E(I)
      N=I
      GO TO 130
120   KITE=I-NH
      KEYI=KEY(KITE)
*      SYMB=SYMBOL(KEYI)
      SYMB=SYMBH(Kite)
      UB=VELEC(KEYI)-E(J)
      N=J
      IF(NP(KEYI).EQ.0) GO TO 130
      UB=UB-E(J+1)-E(J+2)-E(J+3)
      N=J+3
      IF(ND(KEYI).EQ.0) GO TO 130
      UB=UB-E(J+4)-E(J+5)-E(J+6)-E(J+7)-E(J+8)
      N=J+8
130   WRITE(6,1011) SYMB,I,UB,(E(K),K=J,N)
1011  FORMAT(1X,A2,I3,T10,F10.5,T30,9F10.5)
* Own additions to flush out the total charges
      Write(22,1011) SYMB,I,UB
*      close(22)
*
      Do 133 k=j,n
133    Grossp(i)=Grossp(i)+E(k)      
140   J=N+1
      close(22)
*
* Prepare input to the VSS program in the CHEMX system
*
      InumH=0
      InumO=0
      IuIu=1
      Xmax=-1.e20
      Xmin=1.e20
      Ymax=-1.e20
      Ymin=1.e20
      Zmax=-1.e20
      Zmin=1.e20
      Do 5000 i=1,NATM
      Screen(i)=1.
      pot=1.
      Keyi=Key(i)
* if Key is = 0 it is assumed that the missing label is HYDROGEN 
      if(keyi.eq.0) keyi = 35
*
      If(xx(i).gt.xmax) xmax=xx(i)
      if(xx(i).lt.xmin) xmin=xx(i)
      If(yy(i).gt.ymax) ymax=yy(i)
      if(yy(i).lt.ymin) ymin=yy(i)
      If(zz(i).gt.zmax) zmax=zz(i)
      if(zz(i).lt.zmin) zmin=zz(i)
*
       If(Keyi.eq.35) Then
        InumH=InumH+1         
        Ipointer(i)=1
       Else
        InumO=InumO+1
        Ipointer(i)=2          
       EndIf
      Screen(i)=Screen(i)*Exps(Keyi)
      If(NP(keyi).gt.0) Then
       Screen(i)=Screen(i)*Expp(Keyi)
        pot=pot+1.
      End if
      If(ND(keyi).gt.0) Then
       Screen(i)=Screen(i)*Expd(keyi)
        pot=pot+1.
      End if
      Screen(i)=Screen(i)**(1./pot)
5000  Continue

      open(21,FILE="vss.inp",STATUS="UNKNOWN")
      Write(21,'(A)') ' *******NO_INFO*******'
      Write(21,'(A)') ' *******NO_INFO*******'
      Write(21,'(10A8)') AB
      Write(21,'(3(1X,I5))') NATM,InumH,IuIu
*
* First non-hydrogen
*
      Do 5010 i=1,Natm
      If(Ipointer(i).eq.2) Then
       Kay=Key(i)
        Write(21,'(3(1X,f10.4),1x,i5)') xx(i),yy(i),zz(i),velec(kay)
      End if
5010  Continue
      Do 5020 i=1,Natm
      If(Ipointer(i).eq.1) Then
       Kay=Key(i)
        Write(21,'(3(1X,f10.4),1x,i5)') xx(i),yy(i),zz(i),velec(kay)
      End if
5020  Continue

      Do 5030 i=1,Natm
       If(Ipointer(i).eq.2) Then
        write(21,'(F8.5,3x,F8.5)') Grossp(i),Screen(i)
       End if
5030  Continue
      Do 5040 i=1,Natm
       If(Ipointer(i).eq.1) Then                     
        write(21,'(F8.5,3x,F8.5)') Grossp(i),Screen(i)
       End if
5040  Continue 
      Write(21,'(A)') ' 2'
      Write(21,'(6(1x,F8.4))') xmin,xmax,ymin,ymax,zmin,zmax      
      Write(21,'(A)') ' 60  60  60'
C
      close(unit = 21)
C
C  IF CALCULATING ENERGY MATRIX, PUT DIAGONAL ELEMENTS OF HUCKEL
C  MATRIX IN W(I).
C
40    PRT(1)=PRT(17).AND..NOT.PUN(17)
      PRT(2)=PRT(18).AND..NOT.PUN(18)
      PRT(3)=PRT(19).AND..NOT.PUN(19)
      PRT(4)=PRT(20).AND..NOT.PUN(20)
      L1=PRT(1).AND.PRT(2)
      L2=PRT(3).AND.PRT(4)
      IF(L1.AND.L2.AND..NOT.L4) GO TO 9999
      IF(NH.EQ.0) GO TO 23
      DO 22 I=1,NH
22    W(I)=X(I)
23    J=NH+1
      K=NH+1
      DO 24 I=1,NA
      KEYI=KEY(I)
      W(J)=X(K)
      J=J+1
      IF(NP(KEYI).EQ.0) GO TO 24
      UB=Y(K)
      W(J)=UB
      W(J+1)=UB
      W(J+2)=UB
      J=J+3
      IF(ND(KEYI).EQ.0) GO TO 24
      UB=Z(K)
      W(J)=UB
      W(J+1)=UB
      W(J+2)=UB
      W(J+3)=UB
      W(J+4)=UB
      J=J+5
24    K=K+1
C
C  IF DOING CHARGE ITERATION WITH MADELUNG CORRECTION ON MOLECULE
C  WITH NON-ZERO CHARGE, PUT MADELUNG TERMS IN E(I).
C          
      QON= FLOAT(KA)/FLOAT(NELEC)
      IF(METH.LT.3.OR.ABS(QON).LT.0.0001D0) GO TO 180
      IF(NH.EQ.0) GO TO 181
      DO 182 I=1,NH
182   E(I)=MAD(I,I)
181   J=NH+1
      K=NH+1
      DO 183 I=1,NA
      KEYI=KEY(I)
      E(J)=MAD(K,K)
      J=J+1
      K=K+1
      IF(NP(KEYI).EQ.0) GO TO 183
      UB=MAD(K,K)
      E(J)=UB
      E(J+1)=UB
      E(J+2)=UB
      J=J+3
      K=K+1
      IF(ND(KEYI).EQ.0) GO TO 183
      UB=MAD(K,K)
      E(J)=UB
      E(J+1)=UB
      E(J+2)=UB
      E(J+3)=UB
      E(J+4)=UB
      J=J+5
      K=K+1
183   CONTINUE
C
C  CALCULATE AND PRINT OUT ENERGY MATRIX.
C
180   ONEMAT=.TRUE.
      IF(L1) GO TO 7000
170   SUM=1.0D0
      IF(ONEMAT) SUM=0.0D0
      IJ=1
      CNST=2.0D0*CON
      CN2=CNST
      DO 28 I=1,NDIM
      U(I,I)=0.0D0
      DO 28 J=1,I
      UB=0.0D0
      DO 26 K=NMIN,NDIM
26    UB=UB+H(I,K)*H(J,K)*DBLE(IOCC(K))
      IF(I.EQ.J) GO TO 28
      IF(.NOT.L5) GO TO 35
      UC=W(I)
      ET=W(J)
      IF(NHDG.EQ.1) GO TO 36
      UC=HDG(I)
      ET=HDG(J)
36    UC=(UC-ET)/(UC+ET)
      UC=UC*UC
      CNST=CN2+UC+UC*UC*(1.0D0-CN2)
35    IF(METH.LT.3.OR.ABS(QON).LT.0.0001D0) GO TO 31
      UC=UB*C(IJ)*((CNST-SUM)*(W(I)+W(J))-QON*(1.0D0-CNST)*(E(I)+E(J)))
      GO TO 32
31    UC=UB*(CNST-SUM)*C(IJ)*(W(I)+W(J))
32    UB=SUM*UB*C(IJ)
      IJ=IJ+1
      U(J,I)=UC
      U(I,J)=UC
      U(J,J)=U(J,J)+UB*W(J)
28    U(I,I)=U(I,I)+UB*W(I)
      IF(PRT(17)) GO TO 3020
      IF(ONEMAT) WRITE(6,3005)
3005  FORMAT(///,'ENERGY MATRIX')
      IF(.NOT.ONEMAT) WRITE(6,3006)
3006  FORMAT('ENERGY PARTITIONING')
      CALL PEGLEG(U,NDIM,NDIM)
3020  IF(PUN(17)) WRITE(7,2800) U
2800  FORMAT(8F9.5)
C
C  CALL REDUCE TO CALCULATE REDUCED ENERGY MATRIX. PRINT OUT
C  REDUCED ENERGY MATRIX.
C
      IF(PRT(2)) GO TO 7000
      CALL REDUCEM(U,NDIM,NA,NH)
      DO 700 I=2,NATM
      K=I-1
      DO 700 J=1,K
700   U(J,I)=U(I,J)
      IF(PRT(18)) GO TO 3021
      IF(ONEMAT) WRITE(6,3007)
3007  FORMAT('REDUCED ENERGY MATRIX, ATOM BY ATOM')
      IF(.NOT.ONEMAT) WRITE(6,3008)
3008  FORMAT('REDUCED ENERGY PARTITIONING, ATOM BY ATOM')
      KITE=0
701   LOW=KITE+1
      KITE=KITE+13
      IF(KITE.GT.NATM) KITE=NATM
      WRITE(6,702) (I,I=LOW,KITE)
702   FORMAT(/5X,13I9,//)
      DO 703 I=1,NATM
703   WRITE(6,704) I,(U(I,J),J=LOW,KITE)
704   FORMAT(I5,2X,13F9.4)
      IF(KITE.LT.NATM) GO TO 701
3021  IF(PUN(18)) WRITE(7,2111) ((U(I,J),I=1,NATM),J=1,NATM)
2111  FORMAT(7F10.5)
C
C  IF L4 IS TRUE, CALCULATE AND PRINT OUT ENERGY MATRIX ANALYSIS,
C  ORBITAL BY ORBITAL, FOR EACH MOLECULAR ORBITAL SPECIFIED.
C
7000  IF(.NOT.ONEMAT) GO TO 9999
      IF(.NOT.L4) GO TO 71
      DO 246 N=1,23,2
      IF(IENRGY(N).EQ.0) GO TO 72
      KMIN=IENRGY(N)
      KMAX=IENRGY(N+1)
      DO 246 K=KMIN,KMAX
      WRITE(6,3004) K,IOCC(K)
3004  FORMAT(///'ENERGY MATRIX, ORBITAL BY ORBITAL, FOR MOLECULAR ORBIT
     1AL',I4,5X,'OCCUPATION IS',F7.4,' ELECTRONS')
      IF(IOCC(K).LT.0.0001) GO TO 246
      IJ=1
      CNST=CON
      EX=DBLE(IOCC(K))
      DO 244 I=1,NDIM
      DO 244 J=1,I
      UB=H(I,K)*H(J,K)
      IF(I.EQ.J) GO TO 242
      IF(.NOT.L5) GO TO 236
      UC=W(I)
      ET=W(J)
      IF(NHDG.EQ.1) GO TO 237
      UC=HDG(I)
      ET=HDG(J)
237   UC=(UC-ET)/(UC+ET)
      UC=UC*UC
      CNST=CON+UC/2.0D0+UC*UC*(0.5D0-CON)
236   IF(METH.LT.3.OR.ABS(QON).LT.0.0001D0) GO TO 247
      UB=UB*C(IJ)*(CNST*(W(I)+W(J))-QON*(0.5D0-CNST)*(E(I)+E(J)))
      GO TO 248
247   UB=UB*CNST*C(IJ)*(W(I)+W(J))
248   IJ=IJ+1
      UB=2.0D0*UB*EX
      U(I,J)=UB
      U(J,I)=UB
      GO TO 244
242   U(I,I)=UB*W(I)*EX
244   CONTINUE
      CALL PEGLEG(U,NDIM,NDIM)
246   CONTINUE
72    WRITE(6,7003)

C
      close(7)
C
C  CALCULATE AND PRINT OUT ENERGY PARTITIONING AND REDUCED
C  ENERGY PARTITIONING.
C
71    IF(L2) GO TO 9999
      ONEMAT=.FALSE.
      PRT(17)=PRT(19)
      PUN(17)=PUN(19)
      PRT(18)=PRT(20)
      PUN(18)=PUN(20)
      PRT(2)=PRT(4)
      GO TO 170
9999  RETURN
12345 STOP
      END

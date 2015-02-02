 
      PROGRAM SURFAS
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C----                      COPYRIGHT (C) 1988 BY                            ----
C----           RIJKS UNIVERSITEIT GRONINGEN, THE NETHERLANDS               ----
C----                       ALL RIGHTS RESERVED.                            ----
C----                                                                       ----
C---- THIS SOFTWARE IS FURNISHED UNDER A LICENSE AND MAY BE USED AND COPIED ----
C---- ONLY IN ACCORDANCE WITH THE TERMS OF SUCH LICENSE AND WITH THE        ----
C---- INCLUSION OF THE ABOVE COPYRIGHT NOTICE.  THIS SOFTWARE OR ANY OTHER  ----
C---- COPIES THEREOF MAY NOT BE PROVIDED OR OTHERWISE MADE AVAILABLE TO ANY ----
C---- OTHER PERSON.  NO TITLE TO AND OWNERSHIP OF THE SOFTWARE IS HEREBY    ----
C---- TRANSFERRED.                                                          ----
C----                                                                       ----
C---- THE INFORMATION IN THIS SOFTWARE IS SUBJECT TO CHANGE WITHOUT NOTICE  ----
C---- AND SHOULD NOT BE CONSTRUED AS A COMMITMENT BY THE UNIVERSITY OF      ----
C---- GRONINGEN, OR BY THE AUTHORS.                                         ----
C----                                                                       ----
C---- NEITHER THE UNIVERSITY OF GRONINGEN NOR THE AUTHORS ASSUME            ----
C---- RESPONSIBILITY FOR THE USE OR RELIABILITY OF THIS SOFTWARE PRODUCT.   ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- NOTE FROM THE AUTHORS:                                                ----
C---- THE ABOVE IS THE FORMAL TERMINOLOGY APLYING TO ALL UNIVERSITY OF      ----
C---- GRONINGEN'S BIO-MOLECULAR SOFTWARE. THERE ARE HOWEVER NO OBJECTIONS   ----
C---- FOR FREE OF CHARGE REDISTRIBUTION TO NON-PROFIT ORGANISATIONS FOR     ----
C---- THIS SOFTWARE PRODUCT.                                                ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- REMARKABLE FAST VANDER WAALS OR SURFACE ACCESSIBILITY PROGRAM.        ----
C----                                                                       ----
C---- THE SURFACE CALCULATION PART HAS BEEN DESIGNED AND IMPLEMENTED BY:    ----
C---- RICHARD VOORINTHOLT.                                                  ----
C---- THE CONVERSION TO FRODO MAPS, THE DOCUMENTATION AND THE FINAL         ----
C---- ASSEMBLY HAVE BEEN DONE BY G. VRIEND.                                 ----
C----                                                                       ----
C---- THE ROUTINES WHICH HAVE THEIR NAME STARTING WITH 'GV' ARE EXTRACTED   ----
C---- FROM G.VRIEND'S TOOLS LIBRARY. THIS LIBRARY IS AVAILABLE ON REQUEST.  ----
C----                                                                       ----
C---- NO LIMITATIONS ARE PUT ON THE USAGE OF THIS PROGRAM. HOWEVER, IT      ----
C---- WOULD BE APPRECIATED IF PROPER ACKNOWLEDGEMENTS WERE MADE IN ARTICLES ----
C---- IN WHICH RESULTS OF THIS PROGRAM ARE USED.                            ----
C----                                                                       ----
C---- THE PROGRAM CAN BE OBTAINED FREE OF CHARGE FROM                       ----
C---- G. VRIEND                                                             ----
C---- DEPT. CHEMICAL PHYSICS                                                ----
C---- RIJKS UNIVERSITEIT GRONINGEN                                          ----
C---- NIJENBORGH 16                                                         ----
C---- 9747 AG GRONINGEN                                                     ----
C---- THE NETHERLANDS                                                       ----
C----                                                                       ----
C---- BITNET: VRIEND@EMBL                                                   ----
C----                                                                       ----
C---- EXAMPLE FOR USAGE IN BATCH AT A VAX:                                  ----
C----                                                                       ----
C---- THE COORDINATES ARE FOR CRAMBIN IN [PDB.FILES]1CRN.PDB                ----
C---- THE WATER MOLECULES WILL BE SKIPPED UPON READING (Y OR NO IN LINE 2)  ----
C---- ACCESSIBILITIES ARE CALCULATED FOR A MAXIMAL PROBE SIZE OF 4.5 A      ----
C---- THE RESOLUTION IS RATHER FINE (.4 ANGSTROM)                           ----
C---- THE FRODO DSN6 FILE WILL BE CALLED: 1CRNSURF.DN6                      ----
C----                                                                       ----
C---- $ SET DEF [YOUR DIRECTORY]                                            ----
C---- $ RUN CFB:SURFAS                                                      ----
C---- [PDB.FILES]1CRN.PDB                                                   ----
C---- Y                                                                     ----
C---- 4.5                                                                   ----
C---- 0.4                                                                   ----
C---- 1CRNSURF.DN6                                                          ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
C----
C---- REDUNDANT DECLARATIONS JUST FOR CLARITY:
C---- (AND TO HELP PASCAL PROGRAMMERS TO CONVERT THIS PROGRAM)
C---- (REAL MAN DON'T EAT KIRSCH, BUT DO USE FORTRAN)
C----
      REAL           X,Y,Z,VDWRAD,ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,
     +               STEPA,STEPB,STEPC,UVWORG,RESDEF,RESOLU,GVFMIN,
     +               GVFMAX,PROD,PRBRAD,BIGRAD,REX2,FU,RAD2,XR,YR,ZR,
     +               XM,YM,ZM,XG,YG,ZG,DX,DY,DZ,DX2,DY2,DZ2,CELL,
     +               XMDX2,YMDY2,ZMDZ2,XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,
     +               PROBE,FFX,F
      INTEGER        I,NUMBER,NDIVA,NDIVB,NDIVC,NU,NV,NW,ITEMP,IXMIN,
     +               IYMIN,IZMIN,IXMAX,IYMAX,IZMAX,MAXUVW,IOXYZ,IMXYZ,
     +               IGRID,IPLUS,INDEX,IUVW,MP,NX,NY,NZZ,IX,IY,IZ,
     +               GVFFFF,GVFFFC,J,K,II,JJ,KK,I1,I2,I3,ICT,NUMNAM,
     +               IAT,IIX,ICONTR
C----
C---- SERIOUS DECLARATIONS
C----
      DIMENSION      UVWORG(3)
      COMMON /COORD/ X(32000),Y(32000),Z(32000),VDWRAD(32000),NUMBER
C----
C---- READ THE COORDINATES
C----
      CALL GETMCF
C----
C---- CALCULATE THE SURFACE
C----
      CALL SURFS2 (ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,STEPA,STEPB,
     + STEPC,NDIVA,NDIVB,NDIVC,NU,NV,NW,UVWORG)
C----
C---- CREATE THE DSN6 FRODO STYLE MAP
C----
      CALL MAKDN6 (ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,STEPA,STEPB,
     + STEPC,NDIVA,NDIVB,NDIVC,NU,NV,NW,UVWORG)
 
      END
      SUBROUTINE GETMCF
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THIS SUBROUTINE READS AN MCF. AN MCF IS A GRONINGEN STYLE COORDINATE  ----
C---- FILE. THIS SUBROUTINE APPLIES A LITTLE DIRTY TRICK SO THAT IT ALSO    ----
C---- READS PDB (BROOKHAVEN) FILES.                                         ----
C----                                                                       ----
C---- VARIABLES USED:                                                       ----
C---- ATYPE   : CHARACTER*4 ATOM TYPE AS READ IN                            ----
C---- RESIDU  : CHARACTER*4 NAME OF THE RESIDUE (ALA, CYS, ETC) AS READ IN  ----
C---- LINE    : CHARACTER*80 USED TO READ IN THE MCF/PDB LINES. THIS LINE   ----
C----           IS THEN ANALYSED FOR ITS CONTENTS                           ----
C---- GVFFFC  : INTEGER FUNCTION TO READ CHARACTERS FROM THE INPUT STREAM   ----
C----           (UNIT 5). CAN BE REPLACED BY A SIMPLE READ, BUT THEN THE    ----
C----           ERROR DETECTION POSSIBILITIES OF GVFFFC ARE LOST.           ----
C---- GVFYON  : LOGICAL FUNCTION WHICH BECOMES TRUE IF THE QUESTION PASSED  ----
C----           ON TO IT IS ANSWERED WITH Y(ES) OTHERWISE FALSE             ----
C---- X, Y, Z : ARRAYS TO STORE THE X-, Y-, AND Z- COORDINATES OF THE ATOMS ----
C---- VDWRAD  : ARRAY TO STORE THE VAN DER WAALS RADIUS PER ATOM. AS THE    ----
C----           PROGRAM ONLY USES ATOMS, AND DOES NOTHING WITH RESIDUES,    ----
C----           THE ARRAYS X, Y, Z, AND VDWRAD ARE EVERYTHING THE REST OF   ----
C----           THE PROGRAM NEEDS.                                          ----
C---- NUMBER  : TOTAL NUMBER OF ATOMS READ BY THIS SUBROUTINE (MAX = 32000) ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
C----
C---- REDUNDANT DECLARATIONS JUST FOR CLARITY:
C----
      REAL           X,Y,Z,VDWRAD,ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,
     +               STEPA,STEPB,STEPC,UVWORG,RESDEF,RESOLU,GVFMIN,
     +               GVFMAX,PROD,PRBRAD,BIGRAD,REX2,FU,RAD2,XR,YR,ZR,
     +               XM,YM,ZM,XG,YG,ZG,DX,DY,DZ,DX2,DY2,DZ2,CELL,
     +               XMDX2,YMDY2,ZMDZ2,XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,
     +               PROBE,FFX,F
      INTEGER        I,NUMBER,NDIVA,NDIVB,NDIVC,NU,NV,NW,ITEMP,IXMIN,
     +               IYMIN,IZMIN,IXMAX,IYMAX,IZMAX,MAXUVW,IOXYZ,IMXYZ,
     +               IGRID,IPLUS,INDEX,IUVW,MP,NX,NY,NZZ,IX,IY,IZ,
     +               GVFFFF,GVFFFC,J,K,II,JJ,KK,I1,I2,I3,ICT,NUMNAM,
     +               IAT,IIX,ICONTR,NUMSKP
 
      CHARACTER*30   FILE
      CHARACTER*4    ATYPE,RESIDU
      CHARACTER*80   LINE
      LOGICAL        GVFISF,GVFYON,SKIP
      COMMON /COORD/ X(32000),Y(32000),Z(32000),VDWRAD(32000),NUMBER
 
      NUMBER = 0
C----
C---- OPEN THE INPUT FILE
C---- THE DOLLAR CAN BE REMOVED FROM FORMAT 1000 FOR FULL FORTRAN-77
C---- COMPATIBILITY
C----
   10 WRITE (6,1000)
C----
C---- GVFFFC READS CHARACTER STRING. THE FUNCTION VALUE RETURNED IS THE
C---- NUMBER OF CHARACTERS TYPED. FOR FULL FORTRAN-77 COMPATIBILITY THE
C---- NEXT FEW LINES SHOULD BE REPLACED BY:
C----       READ (5,123) LINE
C----   123 FORMAT (A80)
C----
C---- GVSBEL IS A VAX-VMS SPECIFIC SUBROUTINE THAT RINGS THE BELL IN THE
C---- KEYBOARD. CALLS TO THIS SUBROUTINE SHOULD BE REMOVED FOR FULL FORTRAN-77
C---- COMPATIBILITY
C----
      I=GVFFFC(79,LINE)
      IF (I.EQ.0) THEN
         CALL GVSBEL
         WRITE (6,1010)
         GOTO 10
      END IF
C----
C---- THE FOLLOWING FEW LINES ARE USED TO OPEN THE FORMATTED COORDINATE FILE.
C---- THEY CAN, IF WANTED, BE REPLACED BY:
C----
C---- OPEN (UNIT=11,FILE=LINE,FORM='FORMATTED',STATUS='OLD')
C----
C---- BUT THIS IS NOT NEEDED AS GVFISF (LOGICAL FUNCTION TESTING THE EXISTENCE
C---- OF A FILE) AND GVSOPF (OPENING A FILE) ARE SUPPOSED TO BE FULLY
C---- FORTRAN-77 COMPATIBLE.
C----
      IF (GVFISF(LINE)) THEN
         CALL GVSOPF (LINE,11,'OLD')
         REWIND 11
      ELSE
         CALL GVSBEL
         WRITE (6,1020)
         GOTO 10
      END IF
C----
C---- ASK IF THE USER WANTS TO SKIP WATERS
C----
      SKIP=GVFYON('DO YOU WANT TO SKIP WATER MOLECULES (Y/N) ')
      NUMSKP=0
C----
C---- READ NEXT DATA LINE LOOKING FOR AN ATOM/HET RECORD
C----
   20 READ (11,1030,END=40,ERR=30) LINE
      IF (LINE(1:4).EQ.'ATOM'.OR.LINE(1:6).EQ.'HETATM') THEN
         NUMBER=NUMBER+1
         READ (LINE,1040) ATYPE,RESIDU,NUMNAM,X(NUMBER),Y(NUMBER),
     +    Z(NUMBER)
C-------
C------- DIRTY FIX TO READ BOTH GRONINGEN STYLE COORDINATE FILES (MCF'S) AND
C------- PDB FILES WITH THE SAME SUBROUTINE. OTHER USERS MIGHT WANT TO CHANGE
C------- FORMAT 1040 AND REMOVE THIS NEXT STATEMENT
C------- (GVSSHL REMOVES LEADING BLANKS FROM TEXT STRINGS. IT IS FORTRAN-77)
C-------
         CALL GVSSHL(ATYPE)
C-------
C------- IN CASE WATERS SHOULD BE SKIPPED, SEE IF THIS IS A WATER
C-------
         IF (SKIP) THEN
            IF (INDEX('HOHHOH2OWATSOLV',RESIDU(1:3)).NE.0) THEN
               NUMSKP=NUMSKP+1
               NUMBER=NUMBER-1
               GOTO 20
            END IF
         END IF
C-------
C------- SET THE VAN DER WAALS RADIUS AS FUNCTION OF THE ATOM TYPE.
C------- THIS IS THE PLACE WHERE THE PROGRAM SHOULD BE CHANGED IN CASE THE USER
C------- WOULD LIKE TO TREAT WATERS SPECIALLY, ETC.
C-------
         IF (ATYPE(1:1).EQ.'C') THEN
            VDWRAD(NUMBER)=1.9
         ELSE IF (ATYPE(1:1).EQ.'N') THEN
            VDWRAD(NUMBER)=1.7
         ELSE IF (ATYPE(1:1).EQ.'O') THEN
            VDWRAD(NUMBER)=1.4
         ELSE IF (ATYPE(1:1).EQ.'S') THEN
            VDWRAD(NUMBER)=1.8
         ELSE IF (ATYPE(1:1).EQ.'P') THEN
            VDWRAD(NUMBER)=1.8
         ELSE IF (ATYPE(1:1).EQ.'I') THEN
            VDWRAD(NUMBER)=2.0
         ELSE IF (ATYPE(1:1).EQ.'Z') THEN
            VDWRAD(NUMBER)=1.8
         ELSE IF (ATYPE(1:1).EQ.'A') THEN
            VDWRAD(NUMBER)=1.8
         ELSE IF (ATYPE(1:1).EQ.'M') THEN
            VDWRAD(NUMBER)=1.7
         END IF
      END IF
      GOTO 20
C----
C---- HERE UPON ERROR
C---- GVSCLF IS A SUBROUTINE TO CLOSE A FILE
C---- IT CAN BE REPLACED BY:
C----
C---- CLOSE (UNIT=11)
C----
C---- BUT THIS IS NOT NEEDED AS GVSCLF IS FULLY FORTRAN-77 COMPATIBLE
C----
   30 CALL GVSBEL
      WRITE (6,1050)
      CALL GVSCLF (11)
      RETURN
C----
C---- HERE UPON NORMAL END
C----
   40 WRITE (6,1060) NUMBER
      IF (NUMSKP.NE.0) THEN
         WRITE (6,1070) NUMSKP
      ELSE
         WRITE (6,1080)
      END IF
      CALL GVSCLF (11)
C----
C---- FORMATS:
C----
C---- REMOVE THE DOLLAR FROM FORMAT 1000 FOR FULL FORTRAN-77 COMPATIBILITY
C---- CHANGE FORMAT 1040 IS CASE YOU WANT TO USE PDB FILES ONLY AS INPUT
C----
 1000 FORMAT ('$GIVE THE NAME OF THE COORDINATE FILE ')
 1010 FORMAT (' ERROR. YOU SHOULD HAVE GIVEN A FILE NAME')
 1020 FORMAT (' ERROR. THIS FILE DOES NOT EXIST.')
 1030 FORMAT (A80)
 1040 FORMAT (12X,A4,1X,A4,1X,I4,4X,3F8.3,2F6.2)
 1050 FORMAT (' ERROR READING COORDINATE FILE')
 1060 FORMAT (' ',I5,' ATOMS READ FROM COORDINATE FILE')
 1070 FORMAT (' ',I5,' WATER MOLECULES REJECTED')
 1080 FORMAT ('    NO WATER MOLECULES REJECTED')
 
      RETURN
      END
      SUBROUTINE MAKDN6 (ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,STEPA,
     + STEPB,STEPC,NDIVA,NDIVB,NDIVC,NU,NV,NW,UVWORG)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- MAKDN6 WRITE THE GRID OUT IN FRODO DSN6 STYLE. IN CASE YOU WANT OTHER ----
C---- (MORE GENERAL) MAP TYPES, YOU SHOULD CONVERT THR THREE DIMENSIONAL    ----
C---- ARRAY GRID FROM BYTES TO INTEGERS OR REALS, AND USE THE MAP WRITING   ----
C---- LINES THAT ARE PRESENT IN NEXT SUBROUTINE (SURFS2).                   ----
C---- VARIABLES:                                                            ----
C---- MP     : SIZE OF THE GRID ARRAY                                       ----
C---- FILE   : CHARACTER*30 TO STORE THE NAME OF THE OUTPUT FILE NAME       ----
C----          ONLY 30 CHARACTERS BECAUSE FRODO CAN NOT HANDLE MORE IN MOST ----
C----          OF THE VERSIONS FLOATING AROUND                              ----
C---- IOXYZ  :
C---- IUVW   :
C---- IGRID  :
C---- IMXYZ  :
C---- UVWORG :
C---- CELL   : CELL DIMENSIONS A, B, C, ALPHA, BETA, GAMMA, IN ANGSTROMS    ----
C----          AND DEGREES RESPECTIVELY                                     ----
C---- FIRST  : FIRST RECORD OF THE DSN6 MAP (SEE YOUR LOCAL FRODO WRITEUP)  ----
C---- GRID   : THREE DIMENSIONAL BYTE ARRAY TO HOLD THE GENERATED MAP       ----
C---- RECORD : BYTE ARRAY TO STORE BRICKS FOR THE DSN6 FILE                 ----
C---- IRECRD : INTEGER*2 ARRAY TO BE EQUIVALENCED WITH RECORD FOR EASY      ----
C----          MANIPULATION OF THE BYTES AND INTEGER*2'S USED BY FRODO'S    ----
C----          DSN6 FILES. BE AWARE THAT THIS WHOLE BYTE AND INTEGER*2      ----
C----          STUFF HAS NOTHING TO DO WITH FORTRAN-77, BUT THAT IS WHAT    ----
C----          FRODO MAPS ARE MADE OFF....                                  ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
C----
C---- REDUNDANT DECLARATIONS JUST FOR CLARITY:
C----
      REAL           X,Y,Z,VDWRAD,ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,
     +               STEPA,STEPB,STEPC,UVWORG,RESDEF,RESOLU,GVFMIN,
     +               GVFMAX,PROD,PRBRAD,BIGRAD,REX2,FU,RAD2,XR,YR,ZR,
     +               XM,YM,ZM,XG,YG,ZG,DX,DY,DZ,DX2,DY2,DZ2,CELL,
     +               XMDX2,YMDY2,ZMDZ2,XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,
     +               PROBE,FFX,F
      INTEGER        I,NUMBER,NDIVA,NDIVB,NDIVC,NU,NV,NW,ITEMP,IXMIN,
     +               IYMIN,IZMIN,IXMAX,IYMAX,IZMAX,MAXUVW,IOXYZ,IMXYZ,
     +               IGRID,IPLUS,INDEX,IUVW,MP,NX,NY,NZZ,IX,IY,IZ,
     +               GVFFFF,GVFFFC,J,K,II,JJ,KK,I1,I2,I3,ICT,NUMNAM,
     +               IAT,IIX,ICONTR
 
      PARAMETER      (MP=200)
      CHARACTER*30   FILE
      DIMENSION      IOXYZ(3),IUVW(3),IGRID(3),IMXYZ(3),UVWORG(3),
     +               CELL(6)
      INTEGER*2      FIRST(256)
      BYTE           GRID,RECORD(512)
      COMMON /BIGGY/ GRID(MP,MP,MP)
      INTEGER*2      IRECRD(256)
      EQUIVALENCE    (IRECRD(1),RECORD(1))
 
C----
C---- THIS CONVERSION FROM ONE SET OF VARIABLES TO ANOTHER IS MAILY DONE TO
C---- MAKE IS EASIER FOR THOSE WHO KNOW THE FRODO SOURCE CODE A LITTLE TO
C---- UNDERSTAND WHAT ALL VARIABLES MEAN. THE SAME NOMENCLATURE IS USED
C---- AS IS USED IN FRODO ITSELF.
C----
      IUVW(1)  = 1
      IUVW(2)  = 2
      IUVW(3)  = 3
      CELL(1)  = ACELL
      CELL(2)  = BCELL
      CELL(3)  = CCELL
      CELL(4)  = ALPHA
      CELL(5)  = BETA
      CELL(6)  = GAMMA
      IGRID(1) = NDIVA
      IGRID(2) = NDIVB
      IGRID(3) = NDIVC
      IOXYZ(1) = NINT(UVWORG(1)*FLOAT(IGRID(1)))
      IOXYZ(2) = NINT(UVWORG(2)*FLOAT(IGRID(2)))
      IOXYZ(3) = NINT(UVWORG(3)*FLOAT(IGRID(3)))
      IMXYZ(1) = NU
      IMXYZ(2) = NV
      IMXYZ(3) = NW
C----
C---- PROD AND PLUS DON'T NEED ANY VALUE OFCOURSE BECAUSE WE HAVE A SYNTHETIC
C---- MAP FOR WHICH THE EXTREMES ARE KNOWN TO BE ZERO AND HUNDRED.
C----
      PROD=1.0
      IPLUS=0
C----
C---- OPEN THE RANDOM ACCESS OUTPUT
C---- SEE THE COMMENTS ABOUT GVFFFC AND GVSBEL IN THE ROUTINE GETMCF
C----
   10 WRITE (6,1000)
      I=GVFFFC(30,FILE)
      IF (I.EQ.0) THEN
         CALL GVSBEL
         WRITE (6,1010)
         GOTO 10
      END IF
      IF (I.GT.30) THEN
         CALL GVSBEL
         WRITE (6,1020)
         GOTO 10
      END IF
C----
C---- OPEN THE DSN6 FILE
C----
      OPEN (UNIT=2,FILE=FILE,STATUS='NEW',FORM='UNFORMATTED',RECL=128,
     1  ACCESS='DIRECT',MAXREC=32767,RECORDTYPE='FIXED')
C----
C---- GET THE NUMBER OF BRICKS IN THREE DIRECTIONS
C----
      NX=IMXYZ(1)/8
      NY=IMXYZ(2)/8
      NZZ=IMXYZ(3)/8
C----
C---- WRITE OUTPUT FILE HEADER BLOCK (= THE FIRST DSN6 FILE RECORD)
C----
      DO 20 I = 1,3
         FIRST(I  ) = IOXYZ(I)
         FIRST(I+3) = IMXYZ(I)
         FIRST(I+6) = IGRID(I)
         FIRST(I+9) = 100*CELL(I)
         FIRST(I+12)= 100*CELL(I+3)
   20 CONTINUE
      FIRST(16) = 100
      FIRST(17) = 0
      FIRST(18) = 100
      FIRST(19) = 100
      DO 30 I = 20,256
         FIRST(I) = 0
   30 CONTINUE
      WRITE (2,REC=1) FIRST
      INDEX = 1
C----
C---- LOOP OVER THE SECTIONS IN GROUPS OF 8
C----
      DO 90 K=1,NW,8
C-------
C------- LOOP OVER THIS GROUP OF EIGTH PLATES
C------- THE NEXT LOOP IN LOOP IN LOOP CONSTRUCTION CAN BE MADE
C------- MUCH SHORTER AND A LITTLE BIT FASTER IF ALL THESE LOOPS
C------- ARE PLACED IN THE WRITE STATEMENT, BUT THIS WAY IT IS MUCH
C------- CLEARER WHAT IS HAPPENING, AND IT TAKES ONLY A FEW PERCENT
C------- OF THE TOTAL EXECUTION TIME.
C-------
C------- I, J, K RUN OVER THE MAP (IN SEQUENCE FASTEST, MIDDLE, SLOWEST
C------- RUNNING INDEX, IN STEPS OF 8) II, JJ, KK DO THE SAME OVER THE
C------- BRICKS OF 8*8*8 BYTES THAT EACH TIME MAKE UP A BRICK
C-------
         DO 80 J=1,NV,8
            DO 70 I=1,NU,8
               ICT=0
               DO 60 KK=0,7
                  DO 50 JJ=0,7
                     DO 40 II=0,7
                        ICT=ICT+1
                        RECORD(ICT)=GRID(I+II,J+JJ,K+KK)
   40                CONTINUE
   50             CONTINUE
   60          CONTINUE
               INDEX=INDEX+1
               WRITE (2,REC=INDEX) IRECRD
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
C----
C---- CLOSE THE DSN6 FILE (SEE THE NOTE ON GVSCLF IN ROUTINE GETMCF)
C----
      CALL GVSCLF (2)
C----
C---- FORMATS:
C---- REMOVE THE DOLLAR FROM FORMAT 1000 FOR FORTRAN-77 COMPATIBILITY
C----
 1000 FORMAT ('$ENTER NAME OF OUTPUT DSN6 FILE : ')
 1010 FORMAT (' ERROR. YOU SHOULD GIVE A FILE NAME.')
 1020 FORMAT (' ERROR. FILE NAMES SHOULD NOT HAVE OVER 30 CHARACTERS')
 
      RETURN
      END
 
      SUBROUTINE SURFS2 (ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,STEPA,
     + STEPB,STEPC,NDIVA,NDIVB,NDIVC,NU,NV,NW,UVWORG)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THIS SUBROUTINE IS THE MAP GENERATOR. IT ADMINISTRATES THE DIMENSIONS ----
C---- AND THE ORIGIN AND EXTENT OF THE MAP. THE PARTS NEEDED TO GENERATE A  ----
C---- A COMPLETELY GENERAL DENSITY MAP INSTEAD OF A FRODO DSN6 FILE ARE     ----
C---- LOCATED IN THIS SUBROUTINE, BUT COMMENTED OUT.                        ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
C----
C---- REDUNDANT DECLARATIONS JUST FOR CLARITY:
C----
      REAL           X,Y,Z,VDWRAD,ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,
     +               STEPA,STEPB,STEPC,UVWORG,RESDEF,RESOLU,GVFMIN,
     +               GVFMAX,PROD,PRBRAD,BIGRAD,REX2,FU,RAD2,XR,YR,ZR,
     +               XM,YM,ZM,XG,YG,ZG,DX,DY,DZ,DX2,DY2,DZ2,CELL,
     +               XMDX2,YMDY2,ZMDZ2,XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,
     +               FFX,F
      INTEGER        I,NUMBER,NDIVA,NDIVB,NDIVC,NU,NV,NW,ITEMP,IXMIN,
     +               IYMIN,IZMIN,IXMAX,IYMAX,IZMAX,MAXUVW,IOXYZ,IMXYZ,
     +               IGRID,IPLUS,INDEX,IUVW,MP,NX,NY,NZZ,IX,IY,IZ,
     +               GVFFFF,GVFFFC,J,K,II,JJ,KK,I1,I2,I3,ICT,NUMNAM,
     +               IAT,IIX,ICONTR
 
      PARAMETER      (MP=200)
      CHARACTER*30   FILE
      DIMENSION      IOXYZ(3),IUVW(3),IGRID(3),IMXYZ(3),UVWORG(3),
     +               CELL(6)
      CHARACTER*80   LINE
      COMMON /COORD/ X(32000),Y(32000),Z(32000),VDWRAD(32000),NUMBER
      BYTE           GRID
      COMMON /BIGGY/ GRID (MP,MP,MP)
C----
C---- DETERMINE THE EXTREMES OF THE MOLECULE
C---- THE FUNCTIONS GVFMIN AND GVFMAX ARE COMPLETELY FORTRAN-77 COMPATIBLE
C---- FUNCTIONS THAT RETURN THE MINIMAL AND MAXIMAL VALUES OF THE ARRAYS
C---- PASSED ON TO THEM RESPECTIVELY
C----
      XMIN=GVFMIN(X,NUMBER)
      YMIN=GVFMIN(Y,NUMBER)
      ZMIN=GVFMIN(Z,NUMBER)
      XMAX=GVFMAX(X,NUMBER)
      YMAX=GVFMAX(Y,NUMBER)
      ZMAX=GVFMAX(Z,NUMBER)
C----
C---- ASK FOR THE PROBE RADIUS
C---- THE FUNCTION GVFFFF READS A FREE FORMATTED REAL. IT DOES SOME TESTING.
C---- FOR FULL FORTRAN-77 COMPATIBILITY THE FEW LINES INVOLVING GVFFFF SHOULD
C---- BE REPLACED BY:
C----
C---- READ (5,*) PRBRAD
C----
      WRITE (6,1000)
      I=GVFFFF(1,PRBRAD)
      IF (I.EQ.0) PRBRAD=1.4
C----
C---- BIGRAD IS SOMEWHAT BIGGER THAN PRBRAD IN ORDER TO PREVENT PROBLEMS
C---- AT THE EDGE OF THE MAP LATER ON. THE PROGRAMMER SHOULD ONLY TAKE
C---- CARE THAT BIGRAD STAYS IN SIZE APPROXIMATELY THE SAME AS THE LARGEST
C---- POSSIBLE VAN DER WAALS RADIUS
C----
      BIGRAD=1.75+PRBRAD
C----
C---- WRITE OUT WHAT CONTOUR LEVELS SHOULD BE USED FOR DIFFERENT
C---- PROBE SIZES. SEE LATER DOWN IN THE PROGRAM FOR THE MEANING OF
C---- VARIABLES.
C----
      WRITE (6,*) '      PROBE      CONTOUR '
      WRITE (6,*) '      RADIUS:      AT:'
      DO 10 I=0,4
         REX2=(1.8+PRBRAD)**2
         FU=(1.8+PRBRAD*I/4.0)**2
         RAD2=1.8**2
         ICONTR=NINT(100*(REX2-FU)/(REX2-RAD2))
         FU=PRBRAD*I/4.0
         WRITE (6,1010) FU,ICONTR
   10 CONTINUE
C----
C---- UPDATE THE EXTREMES (TO PREVENT ARRAY DIMENSION OVERFLOWS)
C----
      XMIN=XMIN-BIGRAD
      XMAX=XMAX+BIGRAD
      YMIN=YMIN-BIGRAD
      YMAX=YMAX+BIGRAD
      ZMIN=ZMIN-BIGRAD
      ZMAX=ZMAX+BIGRAD
      IXMIN=NINT(XMIN)
      IYMIN=NINT(YMIN)
      IZMIN=NINT(ZMIN)
      IXMAX=NINT(XMAX)
      IYMAX=NINT(YMAX)
      IZMAX=NINT(ZMAX)
C----
C---- DEFINE THE `CELL DIMENSIONS`
C----
      ACELL=XMAX-XMIN
      BCELL=YMAX-YMIN
      CCELL=ZMAX-ZMIN
C----
C---- ASK FOR THE RESOLUTION
C---- RESDEF, THE DEFAULT RESOLUTION IS CHOSEN SUCH AS TO KEEP THE CPU TIME
C---- IN THE ORDER OF A FEW MINUTES ON A MICRO-VAX
C---- SEE THE NOTE ON GVFFFF EARLIER IN THIS SUBROUTINE
C----
   20 RESDEF=(ACELL+BCELL+CCELL)/192.0
      WRITE (6,1020) RESDEF
      I=GVFFFF(1,RESOLU)
      IF (I.EQ.0) RESOLU=RESDEF
C----
C---- CREATE THE CELL ANGLES
C----
      ALPHA=90.0
      BETA=90.0
      GAMMA=90.0
C----
C---- CREATE THE STEP SIZES
C----
      NDIVA=NINT(ACELL/RESOLU)
      NDIVB=NINT(BCELL/RESOLU)
      NDIVC=NINT(CCELL/RESOLU)
      STEPA=ACELL/NDIVA
      STEPB=BCELL/NDIVB
      STEPC=CCELL/NDIVC
C----
C---- CREATE NU,NV,NW
C----
      NU=NDIVA/8
      NV=NDIVB/8
      NW=NDIVC/8
      NU=(NU+1)*8
      NV=(NV+1)*8
      NW=(NW+1)*8
C----
C---- SEE IF THE MAP DOES NOT OVERFLOW THE DIMENSIONS OF GRID
C---- SEE THE REMARKS ON GVSBEL EARLIER IN THIS ROUTINE
C----
      MAXUVW=MAX(NU,MAX(NV,NW))
      IF (MAXUVW.GT.200) THEN
         CALL GVSBEL
         WRITE (6,1030) FLOAT(MAXUVW)/199.0
         GOTO 20
      END IF
C----
C---- INITIALIZE THE GRID
C---- THIS CAN BE REMOVED IF YOU SYSTEM AUTOMATICALLY SETS VARIABLES AT ZERO
C---- BEFORE STARTING EXECUTION.
C----
      DO 50 I1=1,NW
         DO 40 I2=1,NV
            DO 30 I3=1,NU
               GRID(I3,I2,I1)=0
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
C----
C---- GET UVWORG
C----
      UVWORG(1)=XMIN/ACELL
      UVWORG(2)=YMIN/BCELL
      UVWORG(3)=ZMIN/CCELL
C----
C---- IN CASE YOU WANT TO GENERATE A GENERAL MAP, YOU SHOULD INCLUDE THE
C---- FOLLOWING LINES (OR ANYTHING LIKE IT) AROUND HERE TO GET A HEADER
C---- ONTO THE MAP FROM WHICH THE ORIENTATION, ORIGIN, EXTENT, ETC. OF
C---- THE MAP CAN BE DEDUCED LATER.
C---- THE GRONINGEN STYLE MAPS ARE DEFINE AS FOLLOWS:
C---- THE MASTER FOURIER FILE (MFF):
C---- UNFORMATTED FILE STARTING WITH 8 INTRODUCTORY RECORDS FOLLOWED BY THE
C---- ELECTRON DENSITY IN Z-SECTIONS.
C---- 1) TITLE               (CHARATER*80)
C---- 2) A, B, C, ALPHA, BETA, GAMMA (ANGSTROMS AND DEGREES)
C---- 3) STEPU, STEPV, STEPW (STEPSIZE IN U, V, AND W DIRECTION IN ANGSTROMS)
C---- 4) IDIRU, IDIRV, IDIRW (ORIENTATION OF THE MAP; HERE ALWAYS 1, 2, 3)
C---- 5) NDIVA, NDIVB, NDIVC (FRACTIONAL STEP SIZE IN FOURIER CALCULATIONS)
C---- 6) NU, NV, NW          (NUMBER OF GRID POINTS ALONG THE 3 AXES)
C---- 7) UVWORG(3)           (FRACTIONAL COORDINATES OF THE FIRST POINT)
C---- 8) UVWMAT(3,3)         (MATRIX TO CONVERT GRIDPOINTS INTO FRACTIONAL
C----                         COORDINATES)
C----
C      WRITE (6,?) TITLE
C      WRITE (6,?) A,B,C,ALPHA,BETA,GAMMA
C      WRITE (6,?) STEPU,STEPV,STEPW
C      WRITE (6,?) 1,2,3
C      WRITE (6,?) NDIVA,NDIVB,NDIVC
C      WRITE (6,?) NU,NV,NW
C      WRITE (6,?) UVWORG
C      UVWMAT(1,1)=STEPA/ACELL
C      UVWMAT(2,2)=STEPB/BCELL
C      UVWMAT(3,3)=STEPC/CCELL
C      WRITE (6,?) UVWMAT
C----
C---- FILL THE GRID
C---- MANY PARAMETERS REMAIN CONSTANT THROUGHOUT THE CALCULATIONS. THESE ARE
C---- THEREFORE PRECALCULATED TO SAVE CPU TIME.
C----
      DX2=STEPA**2
      DY2=STEPB**2
      DZ2=STEPC**2
      DX=STEPA
      DY=STEPB
      DZ=STEPC
C----
C---- GVSNL6 IS A FORTRAN-77 COMPATIBLE SUBROUTINE TO WRITE EMPTY LINES AT
C---- UNIT 6 (THE TERMINAL)
C----
      CALL GVSNL6 (1)
      DO 60 IAT=1,NUMBER
C-------
C------- THIS COUNTER IS NICE OFCOURSE, BUT IT CAN EASILY BE REMOVED
C-------
         IF (MOD(IAT,25).EQ.0) WRITE (6,1040) IAT,NUMBER
         XM=X(IAT)
         YM=Y(IAT)
         ZM=Z(IAT)
         RAD2=VDWRAD(IAT)**2
         REX2=(VDWRAD(IAT)+PRBRAD)**2
         IX=NINT((X(IAT)-XMIN)/DX)
         IY=NINT((Y(IAT)-YMIN)/DY)
         IZ=NINT((Z(IAT)-ZMIN)/DZ)
         XG=XMIN+IX*DX
         YG=YMIN+IY*DY
         ZG=ZMIN+IZ*DZ
         XMDX2=2*XM*DX
         YMDY2=2*YM*DY
         ZMDZ2=2*ZM*DZ
C-------
C------- FILL THE GRIDPOINT CLOSEST TO THE PRESENT ATOM
C-------
         FU=(XM-XG)**2+(YM-YG)**2+(ZM-ZG)**2
         IF (FU.LE.RAD2) THEN
            GRID(IX+1,IY+1,IZ+1)=100
         ELSE
            ITEMP=GRID(IX+1,IY+1,IZ+1)
            GRID(IX+1,IY+1,IZ+1)=MAX(ITEMP,
     +       NINT(100*((REX2-FU)/(REX2-RAD2))))
         END IF
C-------
C------- NOW WLK IN ALL DIRECTIONS OVER THE GRID BUT INSIDE THE VAN DER WAALS
C------- RADIUS PLUS THE PROBE RADIUS AWAY FROM THE ATOM PRESENTLY INVESTIGATED
C-------
         CALL WLKXPS(FU,XG,YG,ZG,IX,IY,IZ,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         CALL WLKXNG(FU,XG,YG,ZG,IX,IY,IZ,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         CALL WLKYPS(FU,XG,YG,ZG,IX,IY,IZ,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         CALL WLKYNG(FU,XG,YG,ZG,IX,IY,IZ,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         CALL WLKZPS(FU,XG,YG,ZG,IX,IY,IZ,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         CALL WLKZNG(FU,XG,YG,ZG,IX,IY,IZ,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
   60 CONTINUE
C----
C---- IN CASE YOU WANT A GENERAL MAP YOU NEED SOME CODE LIKE:
C----
C     DO K=1,NW
C        WRITE (MAP,FORMAT) ((FLOAT(GRID(I,J,K)),I=1,NU),J=1,NV)
C        FORMAT ....
C     END DO
C     CALL GVSCLF (MAP)
C----
C---- FORMATS:
C---- REMOVE THE DOLLAR FROM FORMATS 1000 AND 1020 FOR FORTRAN-77 COMPATIBILITY
C---- FORMAT 1040 AND ITS ASSOCIATE WRITE CAN BE REMOVED WITHOUT HARM DONE.
C----
 1000 FORMAT ('$GIVE THE PROBE VDW RADIUS (1.4) ')
 1020 FORMAT ('$GIVE THE RESOLUTION (',F4.1,') ')
 1010 FORMAT('  ',F10.1,I10)
 1030 FORMAT (' ERROR. YOUR PROBLEM IS TO BIG FOR THE PROGRAM.',/
     + ' PLEASE REDUCE THE RESOLUTION AT LEAST BY A FACTOR OF',F7.4)
 1040 FORMAT ('+',I4,'/',I4,' ATOMS DONE')
 
      RETURN
      END
 
      SUBROUTINE WLKXNG (FU,XG,YG,ZG,IX,IY,IZ,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THIS SUBROUTINE WLKS OVER THE GRID, BUT ONLY IN NEGATIVE X DIRECTION ----
C---- HOWEVER, IT WILL STAY WITHIN THE LIMITS SET BY THE VAN DER WAALS      ----
C---- RADIUS PLUS THE PROBE RADIUS FOR THE ATOM PRESENTLY INVESTIGATED.     ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
C----
C---- REDUNDANT DECLARATIONS JUST FOR CLARITY:
C----
      REAL           X,Y,Z,VDWRAD,ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,
     +               STEPA,STEPB,STEPC,UVWORG,RESDEF,RESOLU,GVFMIN,
     +               GVFMAX,PROD,PRBRAD,BIGRAD,REX2,FU,RAD2,XR,YR,ZR,
     +               XM,YM,ZM,XG,YG,ZG,DX,DY,DZ,DX2,DY2,DZ2,CELL,
     +               XMDX2,YMDY2,ZMDZ2,XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,
     +               PROBE,FFX,F
      INTEGER        I,NUMBER,NDIVA,NDIVB,NDIVC,NU,NV,NW,ITEMP,IXMIN,
     +               IYMIN,IZMIN,IXMAX,IYMAX,IZMAX,MAXUVW,IOXYZ,IMXYZ,
     +               IGRID,IPLUS,INDEX,IUVW,MP,NX,NY,NZZ,IX,IY,IZ,
     +               GVFFFF,GVFFFC,J,K,II,JJ,KK,I1,I2,I3,ICT,NUMNAM,
     +               IAT,IIX,ICONTR
 
      PARAMETER      (MP=200)
      CHARACTER*30   FILE
      DIMENSION      IOXYZ(3),IUVW(3),IGRID(3),IMXYZ(3),UVWORG(3),
     +               CELL(6)
      BYTE GRID (200,200,200)
 
      FFX=FU+XMDX2+DX2-2*DX*XG
      XR=XG-DX
      IIX=IX-1
   10 CONTINUE
      IF (FFX.LT.REX2) THEN
         IF (FFX.LE.RAD2) THEN
            GRID(IIX+1,IY+1,IZ+1)=100
         ELSE
            ITEMP=GRID(IIX+1,IY+1,IZ+1)
            GRID(IIX+1,IY+1,IZ+1)=MAX(ITEMP,
     +       NINT(100*((REX2-FFX)/(REX2-RAD2))))
         END IF
         FFX=FFX+XMDX2+DX2-2*DX*XR
         XR=XR-DX
         IIX=IIX-1
         GOTO 10
      END IF
 
      RETURN
      END
      SUBROUTINE WLKXPS (FU,XG,YG,ZG,IX,IY,IZ,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THIS SUBROUTINE WLKS OVER THE GRID, BUT ONLY IN POSITIVE X DIRECTION ----
C---- HOWEVER, IT WILL STAY WITHIN THE LIMITS SET BY THE VAN DER WAALS      ----
C---- RADIUS PLUS THE PROBE RADIUS FOR THE ATOM PRESENTLY INVESTIGATED.     ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
C----
C---- REDUNDANT DECLARATIONS JUST FOR CLARITY:
C----
      REAL           X,Y,Z,VDWRAD,ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,
     +               STEPA,STEPB,STEPC,UVWORG,RESDEF,RESOLU,GVFMIN,
     +               GVFMAX,PROD,PRBRAD,BIGRAD,REX2,FU,RAD2,XR,YR,ZR,
     +               XM,YM,ZM,XG,YG,ZG,DX,DY,DZ,DX2,DY2,DZ2,CELL,
     +               XMDX2,YMDY2,ZMDZ2,XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,
     +               PROBE,FFX,F
      INTEGER        I,NUMBER,NDIVA,NDIVB,NDIVC,NU,NV,NW,ITEMP,IXMIN,
     +               IYMIN,IZMIN,IXMAX,IYMAX,IZMAX,MAXUVW,IOXYZ,IMXYZ,
     +               IGRID,IPLUS,INDEX,IUVW,MP,NX,NY,NZZ,IX,IY,IZ,
     +               GVFFFF,GVFFFC,J,K,II,JJ,KK,I1,I2,I3,ICT,NUMNAM,
     +               IAT,IIX,ICONTR
 
      PARAMETER      (MP=200)
      CHARACTER*30   FILE
      DIMENSION      IOXYZ(3),IUVW(3),IGRID(3),IMXYZ(3),UVWORG(3),
     +               CELL(6)
      BYTE GRID (200,200,200)
 
      FFX=FU-XMDX2+DX2+2*DX*XG
      XR=XG+DX
      IIX=IX+1
   10 CONTINUE
      IF (FFX.LT.REX2) THEN
         IF (FFX.LE.RAD2) THEN
            GRID(IIX+1,IY+1,IZ+1)=100
         ELSE
            ITEMP=GRID(IIX+1,IY+1,IZ+1)
            GRID(IIX+1,IY+1,IZ+1)=MAX(ITEMP,
     +       NINT(100*((REX2-FFX)/(REX2-RAD2))))
         END IF
         FFX=FFX-XMDX2+DX2+2*DX*XR
         XR=XR+DX
         IIX=IIX+1
         GOTO 10
      END IF
 
      RETURN
      END
 
      SUBROUTINE WLKYNG (FU,XG,YG,ZG,IX,IY,IZ,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THIS SUBROUTINE WLKS OVER THE GRID, BUT ONLY IN NEGATIVE Y DIRECTION ----
C---- HOWEVER, IT WILL STAY WITHIN THE LIMITS SET BY THE VAN DER WAALS      ----
C---- RADIUS PLUS THE PROBE RADIUS FOR THE ATOM PRESENTLY INVESTIGATED.     ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
C----
C---- REDUNDANT DECLARATIONS JUST FOR CLARITY:
C----
      REAL           X,Y,Z,VDWRAD,ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,
     +               STEPA,STEPB,STEPC,UVWORG,RESDEF,RESOLU,GVFMIN,
     +               GVFMAX,PROD,PRBRAD,BIGRAD,REX2,FU,RAD2,XR,YR,ZR,
     +               XM,YM,ZM,XG,YG,ZG,DX,DY,DZ,DX2,DY2,DZ2,CELL,
     +               XMDX2,YMDY2,ZMDZ2,XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,
     +               PROBE,FFX,F
      INTEGER        I,NUMBER,NDIVA,NDIVB,NDIVC,NU,NV,NW,ITEMP,IXMIN,
     +               IYMIN,IZMIN,IXMAX,IYMAX,IZMAX,MAXUVW,IOXYZ,IMXYZ,
     +               IGRID,IPLUS,INDEX,IUVW,MP,NX,NY,NZZ,IX,IY,IZ,
     +               GVFFFF,GVFFFC,J,K,II,JJ,KK,I1,I2,I3,ICT,NUMNAM,
     +               IAT,IIX,ICONTR
 
      PARAMETER      (MP=200)
      CHARACTER*30   FILE
      DIMENSION      IOXYZ(3),IUVW(3),IGRID(3),IMXYZ(3),UVWORG(3),
     +               CELL(6)
      BYTE GRID (200,200,200)
      SAVE F,YR,I
 
      F=FU+YMDY2+DY2-2*DY*YG
      YR=YG-DY
      I=IY-1
   10 CONTINUE
      IF (F.LT.REX2) THEN
         IF (F.LE.RAD2) THEN
            GRID(IX+1,I+1,IZ+1)=100
         ELSE
            ITEMP=GRID(IX+1,I+1,IZ+1)
            GRID(IX+1,I+1,IZ+1)=MAX(ITEMP,
     +       NINT(100*((REX2-F)/(REX2-RAD2))))
         END IF
C-------
C------- PER STEP IN Y DIRECTION ALL STEPS IN BOTH THE POSITIVE AND NEGATIVE
C------- X DIRECTION ARE TRIED. THIS WAS ORIGINALLY IMPLEMENTED AS A CALL
C------- TO A SUBROUTINE, BUT INLINE SOURCE CODE CAME OUT TO BE SIGNIFICANTLY
C------- FASTER DUE TO THE LARGE NUMBER OF PARAMETERS TO BE PASSED ON
C------- INCLUDE WLKXPOS
C-------
         FFX=F-XMDX2+DX2+2*DX*XG
         XR=XG+DX
         IIX=IX+1
   20    CONTINUE
         IF (FFX.LT.REX2) THEN
            IF (FFX.LE.RAD2) THEN
               GRID(IIX+1,I+1,IZ+1)=100
            ELSE
               ITEMP=GRID(IIX+1,I+1,IZ+1)
               GRID(IIX+1,I+1,IZ+1)=MAX(ITEMP,
     +          NINT(100*((REX2-FFX)/(REX2-RAD2))))
            END IF
            FFX=FFX-XMDX2+DX2+2*DX*XR
            XR=XR+DX
            IIX=IIX+1
            GOTO 20
         END IF
C-------
C------- PER STEP IN Y DIRECTION ALL STEPS IN BOTH THE POSITIVE AND NEGATIVE
C------- X DIRECTION ARE TRIED. THIS WAS ORIGINALLY IMPLEMENTED AS A CALL
C------- TO A SUBROUTINE, BUT INLINE SOURCE CODE CAME OUT TO BE SIGNIFICANTLY
C------- FASTER DUE TO THE LARGE NUMBER OF PARAMETERS TO BE PASSED ON
C------- INCLUDE WLKXNEG
C-------
         FFX=F+XMDX2+DX2-2*DX*XG
         XR=XG-DX
         IIX=IX-1
   30    CONTINUE
         IF (FFX.LT.REX2) THEN
            IF (FFX.LE.RAD2) THEN
               GRID(IIX+1,I+1,IZ+1)=100
            ELSE
               ITEMP=GRID(IIX+1,I+1,IZ+1)
               GRID(IIX+1,I+1,IZ+1)=MAX(ITEMP,
     +          NINT(100*((REX2-FFX)/(REX2-RAD2))))
            END IF
            FFX=FFX+XMDX2+DX2-2*DX*XR
            XR=XR-DX
            IIX=IIX-1
            GOTO 30
         END IF
         F=F+YMDY2+DY2-2*DY*YR
         YR=YR-DY
         I=I-1
         GOTO 10
      END IF
 
      RETURN
      END
      SUBROUTINE WLKYPS (FU,XG,YG,ZG,IX,IY,IZ,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THIS SUBROUTINE WLKS OVER THE GRID, BUT ONLY IN POSITIVE Y DIRECTION ----
C---- HOWEVER, IT WILL STAY WITHIN THE LIMITS SET BY THE VAN DER WAALS      ----
C---- RADIUS PLUS THE PROBE RADIUS FOR THE ATOM PRESENTLY INVESTIGATED.     ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
C----
C---- REDUNDANT DECLARATIONS JUST FOR CLARITY:
C----
      REAL           X,Y,Z,VDWRAD,ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,
     +               STEPA,STEPB,STEPC,UVWORG,RESDEF,RESOLU,GVFMIN,
     +               GVFMAX,PROD,PRBRAD,BIGRAD,REX2,FU,RAD2,XR,YR,ZR,
     +               XM,YM,ZM,XG,YG,ZG,DX,DY,DZ,DX2,DY2,DZ2,CELL,
     +               XMDX2,YMDY2,ZMDZ2,XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,
     +               PROBE,FFX,F
      INTEGER        I,NUMBER,NDIVA,NDIVB,NDIVC,NU,NV,NW,ITEMP,IXMIN,
     +               IYMIN,IZMIN,IXMAX,IYMAX,IZMAX,MAXUVW,IOXYZ,IMXYZ,
     +               IGRID,IPLUS,INDEX,IUVW,MP,NX,NY,NZZ,IX,IY,IZ,
     +               GVFFFF,GVFFFC,J,K,II,JJ,KK,I1,I2,I3,ICT,NUMNAM,
     +               IAT,IIX,ICONTR
 
      PARAMETER      (MP=200)
      CHARACTER*30   FILE
      DIMENSION      IOXYZ(3),IUVW(3),IGRID(3),IMXYZ(3),UVWORG(3),
     +               CELL(6)
      BYTE GRID (200,200,200)
      SAVE F,YR,I
 
      F=FU-YMDY2+DY2+2*DY*YG
      YR=YG+DY
      I=IY+1
   10 CONTINUE
      IF (F.LT.REX2) THEN
         IF (F.LE.RAD2) THEN
            GRID(IX+1,I+1,IZ+1)=100
         ELSE
            ITEMP=GRID(IX+1,I+1,IZ+1)
            GRID(IX+1,I+1,IZ+1)=MAX(ITEMP,
     +       NINT(100*((REX2-F)/(REX2-RAD2))))
         END IF
C-------
C------- PER STEP IN Y DIRECTION ALL STEPS IN BOTH THE POSITIVE AND NEGATIVE
C------- X DIRECTION ARE TRIED. THIS WAS ORIGINALLY IMPLEMENTED AS A CALL
C------- TO A SUBROUTINE, BUT INLINE SOURCE CODE CAME OUT TO BE SIGNIFICANTLY
C------- FASTER DUE TO THE LARGE NUMBER OF PARAMETERS TO BE PASSED ON
C------- INCLUDE WLKXPOS
C-------
         FFX=F-XMDX2+DX2+2*DX*XG
         XR=XG+DX
         IIX=IX+1
   20    CONTINUE
         IF (FFX.LT.REX2) THEN
            IF (FFX.LE.RAD2) THEN
               GRID(IIX+1,I+1,IZ+1)=100
            ELSE
               ITEMP=GRID(IIX+1,I+1,IZ+1)
               GRID(IIX+1,I+1,IZ+1)=MAX(ITEMP,
     +          NINT(100*((REX2-FFX)/(REX2-RAD2))))
            END IF
            FFX=FFX-XMDX2+DX2+2*DX*XR
            XR=XR+DX
            IIX=IIX+1
            GOTO 20
         END IF
C-------
C------- PER STEP IN Y DIRECTION ALL STEPS IN BOTH THE POSITIVE AND NEGATIVE
C------- X DIRECTION ARE TRIED. THIS WAS ORIGINALLY IMPLEMENTED AS A CALL
C------- TO A SUBROUTINE, BUT INLINE SOURCE CODE CAME OUT TO BE SIGNIFICANTLY
C------- FASTER DUE TO THE LARGE NUMBER OF PARAMETERS TO BE PASSED ON
C------- INCLUDE WLKXNEG
C-------
         FFX=F+XMDX2+DX2-2*DX*XG
         XR=XG-DX
         IIX=IX-1
   30    CONTINUE
         IF (FFX.LT.REX2) THEN
            IF (FFX.LE.RAD2) THEN
               GRID(IIX+1,I+1,IZ+1)=100
            ELSE
               ITEMP=GRID(IIX+1,I+1,IZ+1)
               GRID(IIX+1,I+1,IZ+1)=MAX(ITEMP,
     +          NINT(100*((REX2-FFX)/(REX2-RAD2))))
            END IF
            FFX=FFX+XMDX2+DX2-2*DX*XR
            XR=XR-DX
            IIX=IIX-1
            GOTO 30
         END IF
         F=F-YMDY2+DY2+2*DY*YR
         YR=YR+DY
         I=I+1
         GOTO 10
      END IF
 
      RETURN
      END
      SUBROUTINE WLKZNG (FU,XG,YG,ZG,IX,IY,IZ,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THIS SUBROUTINE WLKS OVER THE GRID, BUT ONLY IN NEGATIVE Z DIRECTION ----
C---- HOWEVER, IT WILL STAY WITHIN THE LIMITS SET BY THE VAN DER WAALS      ----
C---- RADIUS PLUS THE PROBE RADIUS FOR THE ATOM PRESENTLY INVESTIGATED.     ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
C----
C---- REDUNDANT DECLARATIONS JUST FOR CLARITY:
C----
      REAL           X,Y,Z,VDWRAD,ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,
     +               STEPA,STEPB,STEPC,UVWORG,RESDEF,RESOLU,GVFMIN,
     +               GVFMAX,PROD,PRBRAD,BIGRAD,REX2,FU,RAD2,XR,YR,ZR,
     +               XM,YM,ZM,XG,YG,ZG,DX,DY,DZ,DX2,DY2,DZ2,CELL,
     +               XMDX2,YMDY2,ZMDZ2,XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,
     +               PROBE,FFX,F
      INTEGER        I,NUMBER,NDIVA,NDIVB,NDIVC,NU,NV,NW,ITEMP,IXMIN,
     +               IYMIN,IZMIN,IXMAX,IYMAX,IZMAX,MAXUVW,IOXYZ,IMXYZ,
     +               IGRID,IPLUS,INDEX,IUVW,MP,NX,NY,NZZ,IX,IY,IZ,
     +               GVFFFF,GVFFFC,J,K,II,JJ,KK,I1,I2,I3,ICT,NUMNAM,
     +               IAT,IIX,ICONTR
 
      PARAMETER      (MP=200)
      CHARACTER*30   FILE
      DIMENSION      IOXYZ(3),IUVW(3),IGRID(3),IMXYZ(3),UVWORG(3),
     +               CELL(6)
      BYTE GRID (200,200,200)
      SAVE F,ZR,I
 
      F=FU+ZMDZ2+DZ2-2*DZ*ZG
      ZR=ZG-DZ
      I=IZ-1
   10 CONTINUE
      IF (F.LT.REX2) THEN
         IF (F.LE.RAD2) THEN
            GRID(IX+1,IY+1,I+1)=100
         ELSE
            ITEMP=GRID(IX+1,IY+1,I+1)
            GRID(IX+1,IY+1,I+1)=MAX(ITEMP,
     +       NINT(100*((REX2-F)/(REX2-RAD2))))
         END IF
C-------
C------- FOR EVERY STEP IN THE Z DIRECTION, TRY ALL STEPS IN THE POSITIVE AND
C------- NEGATIVE X AND Y DIRECTION.
C-------
         CALL WLKXPS (F,XG,YG,ZR,IX,IY,I,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         CALL WLKXNG (F,XG,YG,ZR,IX,IY,I,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         CALL WLKYPS (F,XG,YG,ZR,IX,IY,I,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         CALL WLKYNG (F,XG,YG,ZR,IX,IY,I,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         F=F+ZMDZ2+DZ2-2*DZ*ZR
         ZR=ZR-DZ
         I=I-1
         GOTO 10
      END IF
 
      RETURN
      END
      SUBROUTINE WLKZPS (FU,XG,YG,ZG,IX,IY,IZ,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THIS SUBROUTINE WLKS OVER THE GRID, BUT ONLY IN POSITIVE Y DIRECTION ----
C---- HOWEVER, IT WILL STAY WITHIN THE LIMITS SET BY THE VAN DER WAALS      ----
C---- RADIUS PLUS THE PROBE RADIUS FOR THE ATOM PRESENTLY INVESTIGATED.     ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
C----
C---- REDUNDANT DECLARATIONS JUST FOR CLARITY:
C----
      REAL           X,Y,Z,VDWRAD,ACELL,BCELL,CCELL,ALPHA,BETA,GAMMA,
     +               STEPA,STEPB,STEPC,UVWORG,RESDEF,RESOLU,GVFMIN,
     +               GVFMAX,PROD,PRBRAD,BIGRAD,REX2,FU,RAD2,XR,YR,ZR,
     +               XM,YM,ZM,XG,YG,ZG,DX,DY,DZ,DX2,DY2,DZ2,CELL,
     +               XMDX2,YMDY2,ZMDZ2,XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX,
     +               PROBE,FFX,F
      INTEGER        I,NUMBER,NDIVA,NDIVB,NDIVC,NU,NV,NW,ITEMP,IXMIN,
     +               IYMIN,IZMIN,IXMAX,IYMAX,IZMAX,MAXUVW,IOXYZ,IMXYZ,
     +               IGRID,IPLUS,INDEX,IUVW,MP,NX,NY,NZZ,IX,IY,IZ,
     +               GVFFFF,GVFFFC,J,K,II,JJ,KK,I1,I2,I3,ICT,NUMNAM,
     +               IAT,IIX,ICONTR
 
      PARAMETER      (MP=200)
      CHARACTER*30   FILE
      DIMENSION      IOXYZ(3),IUVW(3),IGRID(3),IMXYZ(3),UVWORG(3),
     +               CELL(6)
      BYTE GRID (200,200,200)
      SAVE F,ZR,I
 
      F=FU-ZMDZ2+DZ2+2*DZ*ZG
      ZR=ZG+DZ
      I=IZ+1
   10 CONTINUE
      IF (F.LT.REX2) THEN
         IF (F.LE.RAD2) THEN
            GRID(IX+1,IY+1,I+1)=100
         ELSE
            ITEMP=GRID(IX+1,IY+1,I+1)
            GRID(IX+1,IY+1,I+1)=MAX(ITEMP,
     +       NINT(100*((REX2-F)/(REX2-RAD2))))
         END IF
C-------
C------- FOR EVERY STEP IN THE Z DIRECTION, TRY ALL STEPS IN THE POSITIVE AND
C------- NEGATIVE X AND Y DIRECTION.
C-------
         CALL WLKXPS (F,XG,YG,ZR,IX,IY,I,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         CALL WLKXNG (F,XG,YG,ZR,IX,IY,I,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         CALL WLKYPS (F,XG,YG,ZR,IX,IY,I,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         CALL WLKYNG (F,XG,YG,ZR,IX,IY,I,GRID,XMDX2,YMDY2,ZMDZ2,
     +    DX,DY,DZ,DX2,DY2,DZ2,REX2,RAD2)
         F=F-ZMDZ2+DZ2+2*DZ*ZR
         ZR=ZR+DZ
         I=I+1
         GOTO 10
      END IF
 
      RETURN
      END
      SUBROUTINE GVSBEL
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- SUBROUTINE BELL RINGS THE BELL IN THE TERMINAL (CONTROL G) TO NOTIFY  ----
C---- THE USER OF AN IMPORTANT MESSAGE.                                     ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- W A R N I N G.  THIS SUBROUTINE IS VAX-VMS AND VT*** SPECIFIC         ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      BYTE      RING
      DATA      RING /"007/
 
      WRITE (6,1000) RING
 1000 FORMAT(1H+,A1)
 
      RETURN
      END
      REAL FUNCTION GVFMIN (ROW,NUM)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- DETERMINES THE VALUE OF THE SMALLEST NUMBER IN ROW. ROW IS A REAL     ----
C---- ARRAY AND IS DIMENSIONED ROW(NUM). NUM IS AN INTEGER.                 ----
C----                                                                       ----
C---- PLEASE BE AWARE THAT STRANGE THINGS CAN HAPPEN IN CASE YOUR ARRAYS    ----
C---- DO NOT HAVE THE PROPER LENGTH.                                        ----
C----                                                                       ----
C---- THIS ROUTINE IS WRITTEN IN FORTRAN 77.                                ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      REAL ROW(NUM)
 
      IF (NUM.GT.0) THEN
C-------
C------- IF NUM IS POSITIVE, PERFORM THE ACTION
C-------
         GVFMIN=ROW(1)
         DO 10 I=1,NUM
            GVFMIN=MIN(GVFMIN,ROW(I))
   10    CONTINUE
         RETURN
      ELSE IF (NUM.EQ.0) THEN
C-------
C------- IF NUM IS ZERO, SOMETHING COULD BE WRONG, SO WARN THE USER
C-------
         WRITE (6,1000)
         RETURN
      ELSE IF (NUM.LT.0) THEN
C-------
C------- IF NUM IS NEGATIVE, SOMETHING IS WRONG, SO WARN THE USER
C-------
         WRITE (6,1010) NUM
      END IF
 
 1000 FORMAT (' TLS> WARNING. GVFMIN CALLED WITH NUM = 0')
 1010 FORMAT (' TLS> ERROR. GVFMIN CALLED WITH NUM =',I6)
 
      RETURN
      END
      REAL FUNCTION GVFMAX (ROW,NUM)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- DETERMINES THE VALUE OF THE LARGEST NUMBER IN ROW. ROW IS A REAL      ----
C---- ARRAY AND IS DIMENSIONED ROW(NUM). NUM IS AN INTEGER.                 ----
C----                                                                       ----
C---- PLEASE BE AWARE THAT STRANGE THINGS CAN HAPPEN IN CASE YOUR ARRAYS    ----
C---- DO NOT HAVE THE PROPER LENGTH.                                        ----
C----                                                                       ----
C---- THIS ROUTINE IS WRITTEN IN FORTRAN 77.                                ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      REAL ROW(NUM)
 
      IF (NUM.GT.0) THEN
C-------
C------- IF NUM IS POSITIVE, PERFORM THE ACTION
C-------
         GVFMAX=ROW(1)
         DO 10 I=1,NUM
            GVFMAX=MAX(GVFMAX,ROW(I))
   10    CONTINUE
         RETURN
      ELSE IF (NUM.EQ.0) THEN
C-------
C------- IF NUM IS ZERO, SOMETHING COULD BE WRONG, SO WARN THE USER
C-------
         WRITE (6,1000)
         RETURN
      ELSE IF (NUM.LT.0) THEN
C-------
C------- IF NUM IS NEGATIVE, SOMETHING IS WRONG, SO WARN THE USER
C-------
         WRITE (6,1010) NUM
      END IF
 
 1000 FORMAT (' TLS> WARNING. GVFMAX CALLED WITH NUM = 0')
 1010 FORMAT (' TLS> ERROR. GVFMAX CALLED WITH NUM =',I6)
 
      RETURN
      END
      INTEGER FUNCTION GVFFFC (N,X)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- READS A CHARACTER*N FROM THE TERMINAL.                                ----
C---- THE FUNCTION RETURNS THE ACTUAL LENGTH OF THE CHARACTER STRING READ.  ----
C---- N IS AN INTEGER, BEING THE LENGTH OF THE CHARACTER TO BE READ.        ----
C---- X IS A CHARACTER STRING OF MAXIMUM LENGTH 80.                         ----
C---- WARNING: THIS FUNCTION IS VAX-VMS SPECIFIC.                           ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      CHARACTER*(*)  X
      INTEGER        GVFQQC
      BYTE           REC(80)
      COMMON /GVCIO1/LREC,NEXT,REC
C----
C---- READ THE INPUT LINE
C----
      CALL GVFQQR(5)
C----
C---- DECODE THE INPUT LINE
C----
      GVFFFC=GVFQQC(N,X)
C----
C---- PADD THE INPUT LINE WITH BLANKS
C----
      IF (GVFFFC.LT.N) THEN
         DO 100 I=GVFFFC+1,N
            X(I:I)=' '
  100    CONTINUE
      END IF
 
      RETURN
      END
      INTEGER FUNCTION GVFFFI(N,I)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- GVFFFI READS N INTEGERS FROM THE TERMINAL. THE FUNCTION RETURNS THE   ----
C---- ACTUAL NUMBER OF INTEGERS FOUND ON INPUT.                             ----
C---- N IS AN INTEGER.                                                      ----
C---- I IS AN ARRAY OF INTEGERS.                                            ----
C---- WARNING: THIS FUNCTION IS VAX-VMS SPECIFIC.                           ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      INTEGER    GVFQQI,I(N)
C----
C---- READ THE INPUT LINE FROM THE TERMINAL
C----
      CALL GVFQQR(5)
C----
C---- INTERPRET THE INPUT LINE AS N INTEGERS
C----
      GVFFFI = GVFQQI(N,I)
 
      RETURN
      END
      LOGICAL FUNCTION GVFYON (TEXT)
C///
C------------------------------------------------------------------------------
C----                                                                      ----
C---- PRINT THE INPUT QUESTION AND ACCEPT/TEST THE Y/N ANSWER              ----
C----                                                                      ----
C------------------------------------------------------------------------------
C----                                                                      ----
C---- YORN HAS TO BE DECLARED AS A LOGICAL IN THE CALLING SUBROUTINE       ----
C---- TEXT IS A CHARACTER STRING OF LENGTH 1-100 WHICH WILL BE PRINTED     ----
C----                                                                      ----
C---- WARNING, FORMAT 1000 HAS TO BE CHANGED TO MAKE THIS FORTRAN 77       ----
C----                                                                      ----
C------------------------------------------------------------------------------
C\\\
      CHARACTER*(*) TEXT
      CHARACTER*1   YORNIN,TXTROW(100)
      INTEGER       GVFLEN
C------------------------------------------------------------------------------
C----                                                                      ----
C---- PRINT THE TEXT                                                       ----
C----                                                                      ----
C------------------------------------------------------------------------------
      LENTXT=GVFLEN(TEXT)
      DO 5 I=1,LENTXT
         TXTROW(I)=TEXT(I:I)
    5 CONTINUE
      TXTROW(LENTXT+1)=' '
   10 WRITE (6,1000) (TXTROW(I),I=1,LENTXT+1)
C------------------------------------------------------------------------------
C----                                                                      ----
C---- ACCEPT/TEST THE USER INPUT                                           ----
C----                                                                      ----
C------------------------------------------------------------------------------
      READ (5,1010) YORNIN
      CALL GVSCUC (YORNIN)
      IF (YORNIN.NE.'N'.AND.YORNIN.NE.'Y') THEN
         CALL GVSBEL
         WRITE (6,1020)
         GOTO 10
      END IF
C----
C---- SET THE FUNCTION VALUE
C----
      GVFYON=YORNIN.EQ.'Y'
 
      RETURN
 
 1000 FORMAT ('$TLS> ',100A1)
 1010 FORMAT (A1)
 1020 FORMAT (' TLS> PLEASE GIVE Y OR N')
 
      END
      SUBROUTINE GVSOPF (TEXT,NUMUNT,TYPE)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- FILE 'TEXT' WILL BE OPENED AT UNIT NUMUNT.                            ----
C---- TEXT IS A CHARACTER STRING OF MAXIMUM LENGTH 80.                      ----
C---- NUMUNT IS A POSITIVE INTEGER.                                         ----
C---- TYPE IS VARIABLE LENGTH CHARACTER STRING.                             ----
C---- TYPE = 'OLD'     : OPEN AN EXISTING FILE, WARN UPON ERROR.            ----
C---- TYPE = 'NEW'     : CREATE A NEW FILE, WARN IF THE FILE ALREADY EXISTS ----
C---- TYPE = 'UNKNOWN' : OPEN AN EXISTING FILE, IF IT DOES NOT EXIST,       ----
C----                    THEN CREATE IT.                                    ----
C---- TYPE = 'WARN'    : OPEN A NEW FILE, IF IT ALREADY EXISTS, OPEN THE    ----
C----                    FILE WITH A HIGHER VERSION NUMBER (VAX!) AND GIVE  ----
C----                    A MESSAGE AT THE TERMINAL                          ----
C---- TYPE = 'NOWARN'  : OPEN A NEW FILE, IF IT ALREADY EXISTS, OPEN THE    ----
C----                    FILE WITH A HIGHER VERSION NUMBER (VAX!) AND GIVE  ----
C----                    NO MESSAGES AT THE TERMINAL                        ----
C---- IN THE ABOVE WARNING MEANS PRINT A MESSAGE, AND DO NOT PERFORM THE    ----
C---- ACTION.                                                               ----
C---- IF THE LAST CHARACTER OF TYPE IS AN * (EG 'OLD*', 'NEW*', 'UNKNOWN*', ----
C---- OR 'WARN*') AND NUMUNT IS STILL OPEN, IT IS CLOSED FIRST. OTHERWISE   ----
C---- IF NUMUNT IS OPEN UPON ENTERING THIS ROUTINE A ERROR MESSAGE IS GIVEN ----
C---- AND NO ACTION IS PERFORMED.                                           ----
C----                                                                       ----
C---- THIS ROUTINE IS WRITTEN IN FORTRAN-77.                                ----
C---- HOWEVER THIS PART OF FORTRAN-77 IS NOT IMPLEMENTED IDENTICAL ON ALL   ----
C---- COMPUTERS.                                                            ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      CHARACTER*(*) TEXT,TYPE
      LOGICAL       GVFFOP,ISIT
      INTEGER       GVFLEN
C----
C---- IF NUMUNT IS POSITIVE, PERFORM THE ACTION
C----
      IF (NUMUNT.GT.0) THEN
         LENTXT=GVFLEN(TEXT)
C-------
C------- DETERMINE WETHER TYPE ENDS AT AN *
C-------
         LENTYP=GVFLEN(TYPE)
C-------
C------- SEE WHAT TO DO IN CASE NUMUNT IS STILL OPEN
C-------
         IF (TYPE(LENTYP:LENTYP).EQ.'*') THEN
            IF (GVFFOP(NUMUNT)) CLOSE (UNIT=NUMUNT,STATUS='KEEP')
         ELSE
            IF (GVFFOP(NUMUNT)) THEN
               WRITE (6,1000) NUMUNT
               RETURN
            END IF
         END IF
C-------
C------- DEAL WITH THE CASE THAT TYPE = 'OLD'
C-------
         IF (TYPE(1:3).EQ.'OLD') THEN
            OPEN (UNIT=NUMUNT,FILE=TEXT,ACCESS='SEQUENTIAL',FORM=
     +       'FORMATTED',STATUS='OLD',READONLY,ERR=999)
C-------
C------- DEAL WITH THE CASE THAT TYPE = 'NEW'
C-------
         ELSE IF (TYPE(1:3).EQ.'NEW') THEN
C----------
C---------- FIRST TEST IF THIS FILE ALREADY EXISTS
C----------
            INQUIRE (FILE=TEXT,EXIST=ISIT)
            IF (ISIT) THEN
               WRITE (6,1010) (TEXT(I:I),I=1,LENTXT)
               RETURN
            END IF
            OPEN (UNIT=NUMUNT,FILE=TEXT,ACCESS='SEQUENTIAL',FORM=
     +       'FORMATTED',STATUS='NEW',ERR=998)
C-------
C------- DEAL WITH THE CASE THAT TYPE = 'WARN'
C-------
         ELSE IF (TYPE (1:4).EQ.'WARN') THEN
            INQUIRE (FILE=TEXT,EXIST=ISIT)
            IF (ISIT) THEN
               WRITE (6,1020) (TEXT(I:I),I=1,LENTXT)
            END IF
            OPEN (UNIT=NUMUNT,FILE=TEXT,ACCESS='SEQUENTIAL',FORM=
     +       'FORMATTED',STATUS='NEW',ERR=997)
C-------
C------- DEAL WITH THE CASE THAT TYPE = 'NOWARN'
C-------
         ELSE IF (TYPE (1:6).EQ.'NOWARN') THEN
            OPEN (UNIT=NUMUNT,FILE=TEXT,ACCESS='SEQUENTIAL',FORM=
     +       'FORMATTED',STATUS='NEW',ERR=997)
C-------
C------- DEAL WITH THE CASE THAT TYPE = 'UNKNOWN'
C-------
         ELSE IF (TYPE(1:7).EQ.'UNKNOWN') THEN
            OPEN (UNIT=NUMUNT,FILE=TEXT,ACCESS='SEQUENTIAL',FORM=
     +       'FORMATTED',STATUS='OLD',ERR=20)
            RETURN
C----------
C---------- THE FILE DOES NOT EXIST, SO CREATE IT
C----------
   20       CONTINUE
            CLOSE (UNIT=NUMUNT)
            OPEN (UNIT=NUMUNT,FILE=TEXT,ACCESS='SEQUENTIAL',FORM=
     +       'FORMATTED',STATUS='NEW',ERR=997)
C-------
C------- DEAL WITH INCORRECT TYPE PARAMETER
C-------
         ELSE
            WRITE (6,1030)
            WRITE (6,1090) (TEXT(I:I),I=1,LENTXT)
         END IF
         RETURN
      ELSE IF (NUMUNT.LE.0) THEN
C-------
C------- NEGATIVE OR ZERO UNIT NUMBERS ARE NOT ALLOWED
C-------
         WRITE (6,1080) NUMUNT
         WRITE (6,1090) (TEXT(I:I),I=1,LENTXT)
      END IF
C----
C---- ERRORS IN THE ABOVE OPEN STATEMENTS:
C----
  999 CONTINUE
      WRITE (6,1040)
      WRITE (6,1090) (TEXT(I:I),I=1,LENTXT)
      RETURN
 
  998 CONTINUE
      WRITE (6,1050)
      WRITE (6,1090) (TEXT(I:I),I=1,LENTXT)
      RETURN
 
  997 CONTINUE
      WRITE (6,1050)
      WRITE (6,1090) (TEXT(I:I),I=1,LENTXT)
      RETURN
 
  996 CONTINUE
      WRITE (6,1060)
      WRITE (6,1090) (TEXT(I:I),I=1,LENTXT)
      RETURN
 
 1000 FORMAT (' TLS> ERROR. UNIT NUMBER',I3,' IS ALREADY OPEN',
     + ' IN GVSOPF')
 1010 FORMAT (' TLS> ERROR. THERE ALREADY EXISTS A FILE CALLED:',80A1,
     + /' TLS> PLEASE USE TYPE IS `WARN` IF YOU WANT TO OPEN'
     + ' THIS FILE ANYWAY. (ERROR OCCURRED IN GVSOPF).')
 1020 FORMAT (' TLS> WARNING. THERE ALREADY EXISTS A FILE CALLED:',80A1,
     + /' TLS> A NEW FILE WITH THE SAME NAME, BUT A HIGHER EXTENSION',
     + ' WILL BE OPENED')
 1030 FORMAT (' TLS> ERROR. GVSOPF RECEIVED AN INCORRECT TYPE',
     + ' PARAMETER')
 1040 FORMAT (' TLS> ERROR. OPENING OLD FILE FAILED IN GVSOPF.',
     + ' (TYPE=`OLD`)')
 1050 FORMAT (' TLS> ERROR. OPENING NEW FILE FAILED IN GVSOPF.',
     + ' (TYPE=`NEW`)')
 1060 FORMAT (' TLS> ERROR. OPENING NEW FILE FAILED IN GVSOPF.',
     + ' (TYPE=`WARN`)')
 1070 FORMAT (' TLS> ERROR. OPENING FILE FAILED IN GVSOPF.',
     + ' (TYPE=`UNKNOWN`)')
 1080 FORMAT (' TLS> ERROR. GVSOPF CALLED WITH UNIT NUMBER=',I4)
 1090 FORMAT (' TLS> THE FILE NAME WAS:',80A1)
 
      END
      SUBROUTINE GVSSHL (TEXT)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- A LEFT SHIFT IS PERFORMED ON TEXT UNTIL ALL LEADING BLANCKS ARE GONE  ----
C---- TEXT SHOULD BE A CHARACTER STRING OF UNLIMITTED LENGTH                ----
C----                                                                       ----
C---- THIS ROUTINE IS WRITTEN IN FORTRAN 77.                                ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      CHARACTER*(*) TEXT
      INTEGER       GVFLEN
C----
C---- DETERMINE THE LENGTH OF THE TEXT
C----
      LENTXT=GVFLEN(TEXT)
      IF (LENTXT.GT.0) THEN
C-------
C------- IF THE TEXT HAS CHARACTERS, SEE IF IT STARTS WITH A BLANCK
C-------
         IF (TEXT(1:1).NE.' ') RETURN
C-------
C------- IF THERE ARE LEADING BLANCKS, LOOK FOR THE FIRST CHARACTER
C-------
         DO 30 I=2,LENTXT
            IF (TEXT(I:I).NE.' ') THEN
C-------------
C------------- ASAP THE FIRST NON BLANCK IS FOUND, SHIFT ALL CHARS TO THE LEFT
C-------------
               DO 10 J=1,LENTXT-I+1
                  TEXT(J:J)=TEXT(J+I-1:J+I-1)
   10          CONTINUE
C-------------
C------------- PADD THE RIGHT SIDE OF THE STRING WITH BLANCKS
C-------------
               DO 20 J=LENTXT-I+2,LENTXT
                  TEXT(J:J)=' '
   20          CONTINUE
               RETURN
            END IF
   30    CONTINUE
      ELSE
C-------
C------- GIVE A MESSAGE UPON EMPTY STRINGS
C-------
         WRITE (6,1000)
      END IF
 
 1000 FORMAT (' TLS> LENGTH OF INPUT TEXT IS ZERO IN GVSSHL')
 
      RETURN
      END
      LOGICAL FUNCTION GVFISF (FILNAM)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- BECOMES TRUE IF THE FILE FILNAM EXISTS.                               ----
C----                                                                       ----
C---- THIS FUNCTION IS WRITTEN IN FORTRAN 77.                               ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      CHARACTER*(*)  FILNAM
      INTEGER        GVFLEN
C----
C---- TEST FOR EMPTY FILE NAME
C----
      IF (GVFLEN(FILNAM).EQ.0) THEN
         ISFIL=.FALSE.
         RETURN
      END IF
C----
C---- SEE IF THE FILE EXISTS
C----
      INQUIRE (FILE=FILNAM,EXIST=GVFISF)
 
      RETURN
      END
      SUBROUTINE GVSCLF (NUMUNT)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- GVSCLF CLOSES THE OPEN AT UNIT=NUMUNT. CLOSING IS DONE WITH STATUS    ----
C---- IS 'KEEP'. A WARNING IS SEND OUT IF THE UNIT IS NOT OPEN.             ----
C----                                                                       ----
C---- THIS ROUTINE IS WRITTEN IN FORTRAN 77.                                ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      LOGICAL GVFFOP
C----
C---- IF NUMUNT IS POSITIVE, PERFORM THE ACTION
C----
      IF (NUMUNT.GT.0) THEN
C-------
C------- SEE WHAT TO DO IN CASE NUMUNT IS OPEN
C-------
         IF (GVFFOP(NUMUNT)) THEN
            CLOSE (UNIT=NUMUNT,STATUS='KEEP')
            RETURN
C-------
C------- SEE WHAT TO DO IN CASE NUMUNT IS ALREADY CLOSED
C-------
         ELSE
            WRITE (6,1000) NUMUNT
            RETURN
         END IF
      ELSE IF (NUMUNT.LE.0) THEN
C-------
C------- NEGATIVE OR ZERO UNIT NUMBERS ARE NOT ALLOWED
C-------
         WRITE (6,1010) NUMUNT
      END IF
 
 1000 FORMAT (' TLS> ERROR. UNIT NUMBER',I3,' IS ALREADY CLOSED',
     + ' IN GVSCLF')
 1010 FORMAT (' TLS> ERROR. GVSCLF CALLED WITH UNIT NUMBER=',I4)
 
      END
      SUBROUTINE GVSNL6 (NL)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- SUBROUTINE GVSNL6 SENDS NL EMPTY LINES TO OUTPUT UNIT 6               ----
C---- NL = AN INTEGER INDICATING THE NUMBER OF EMPTY LINES                  ----
C---- WARNING, IF UNIT 6 IS NOT AVAILABLE TO WRITE ON, FATALITIES OCCUR     ----
C----                                                                       ----
C---- THIS ROUTINE IS WRITTEN IN FORTRAN 77.                                ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
 
      IF (NL.GT.0) THEN
C-------
C------- IF NL IS POSITIVE, PERFORM THE ACTION
C-------
         DO 10 I=1,NL
            WRITE (6,1020)
   10    CONTINUE
         RETURN
      ELSE IF (NL.EQ.0) THEN
C-------
C------- IF NL IS ZERO, SOMETHING COULD BE WRONG, SO WARN THE USER
C-------
         WRITE (6,1000)
         RETURN
      ELSE IF (NL.LT.0) THEN
C-------
C------- IF NL IS NEGATIVE, SOMETHING IS WRONG, SO WARN THE USER
C-------
         WRITE (6,1010) NL
      END IF
 
 1000 FORMAT (' TLS> WARNING. GVSNL6 CALLED WITH NL = 0')
 1010 FORMAT (' TLS> ERROR. GVSNL6 CALLED WITH NL =',I6)
 1020 FORMAT(' ')
 
      RETURN
      END
      INTEGER FUNCTION GVFFFF(N,X)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- GVFFFF READS N REALS*4 FROM THE TERMINAL. THE FUNCTION RETURNS THE    ----
C---- ACTUAL NUMBER OF REALS FOUND ON INPUT.                                ----
C---- N IS AN INTEGER.                                                      ----
C---- X IS AN ARRAY OF REALS*4                                              ----
C---- WARNING: THIS FUNCTION IS VAX-VMS SPECIFIC.                           ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      REAL*4    X(N)
      INTEGER   GVFQQF
C----
C---- READ THE INPUT LINE
C----
      CALL GVFQQR(5)
C----
C---- INTERPRET THE INPUT LINE AS REALS
C----
      GVFFFF = GVFQQF(N,X)
 
      RETURN
      END
      LOGICAL FUNCTION GVFFOP (NUMUNT)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THIS FUNCTION BECOMES TRUE IF UNIT NUMUNT IS OPEN                     ----
C---- IF UNIT NUMUNT IS CLOSED, FALSE WILL BE RETURNED                      ----
C---- BEING OPENED DOES NOT MEAN THAT WRITE PERMISSION EXISTS               ----
C---- THEREFORE, THIS IS NOT A SUFFICIENT TEST TO SEE IF SOMETHING CAN BE   ----
C---- WRITTEN ONTO UNIT NUMUNT.                                             ----
C----                                                                       ----
C---- NUMUNT IS AN INTEGER                                                  ----
C---- GVFFOP SHOULD BE DECLARED LOGICAL                                     ----
C----                                                                       ----
C---- THIS ROUTINE IS WRITTEN IN FORTRAN 77.                                ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      LOGICAL ISOPEN
 
 
      IF (NUMUNT.GT.0) THEN
C-------
C------- IF NUMUNT IS POSITIVE, PERFORM THE ACTION
C-------
         INQUIRE (UNIT=NUMUNT,OPENED=ISOPEN)
         GVFFOP=ISOPEN
         RETURN
      ELSE IF (NUMUNT.LE.0) THEN
C-------
C------- NEGATIVE OR ZERO UNIT NUMBERS ARE NOT ALLOWED
C-------
         WRITE (6,1000) NUMUNT
         GVFFOP=.FALSE.
      END IF
 
 1000 FORMAT (' TLS> ERROR. GVFFOP CALLED WITH NUMUNT =',I6)
 
      RETURN
      END
      INTEGER FUNCTION GVFLEN (TEXT)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- RETURNS THE LENGTH OF TEXT                                            ----
C---- THE LENGTH IS DEFINED AS THE SEQUENCE NUMBER IN TEXT OF THE LAST      ----
C---- NON-BLANK CHARACTER. UNDEFINED CHARACTERS ARE ASSUMED TO BE BLANK.    ----
C----                                                                       ----
C---- TEXT IS A CHARACTER STRING OF ANY LENGTH                              ----
C---- DO NOT FORGET TO DECLARE GVFLEN AS AN INTEGER IN THE CALLING ROUTINE  ----
C----                                                                       ----
C---- THIS ROUTINE IS WRITTEN IN FORTRAN 77.                                ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      CHARACTER*(*) TEXT
 
      GVFLEN=0
      DO 10 I=LEN(TEXT),1,-1
         IF (TEXT(I:I).NE.' ') THEN
            GVFLEN=I
            RETURN
         END IF
   10 CONTINUE
 
      RETURN
      END
      INTEGER FUNCTION GVFQQC (N,X)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THIS ROUTINE INTERPRETS THE DATA LINE STORED IN ARRAY REC AS A        ----
C---- CHARACTER STRING OF LENGTH N. THE FUNCTION RETURNS THE NUMBER OF      ----
C---- CHARACTERS ACTUALLY READ IN.                                          ----
C---- NA IS AN INTEGER.                                                     ----
C---- AA IS A CHARACTER STRING OF MAXIMUM LENGTH 80                         ----
C---- WARNING: THIS FUNCTION IS VAX-VMS SPECIFIC.                           ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      CHARACTER*(*)   X
      BYTE            REC(80)
      COMMON /GVCIO1/ LREC, NEXT, REC
C----
C---- DETERMINE THE LENGTH OF THE CHARACTER STRING
C----
      GVFQQC = MIN0 (LREC-NEXT,N)
      IF (GVFQQC .LE. 0) GVFQQC = 0
      IF (GVFQQC.LE.0) RETURN
C----
C---- WRITE THE BYTES INTO THE CHARACTER STRING
C----
      DO 100  I = 1,GVFQQC
         WRITE (X(I:I),1000) REC(NEXT)
         NEXT = NEXT + 1
  100 CONTINUE
      RETURN
 
 1000 FORMAT (A1)
 
      END
      INTEGER FUNCTION GVFQQI (NI,II)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THIS ROUTINE INTERPRETS THE DATA LINE STORED IN ARRAY REC AS A SERIES ----
C---- OF INTEGERS. THE FUNCTION RETURNS THE NUMBER OF INTEGERS ACTUALLY     ----
C---- READ IN.                                                              ----
C---- NA IS AN INTEGER.                                                     ----
C---- II IS AN ARRAY OF INTEGERS.                                           ----
C---- WARNING: THIS FUNCTION IS VAX-VMS SPECIFIC.                           ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      COMMON    /GVCIO1/ LREC, NEXT, REC
      BYTE      REC(80)     ,NUM(10)     ,SPACE       ,TAB
      BYTE      COMMA       ,CHAR        ,MINUS
      LOGICAL   FOUND
      INTEGER   II(NI)
      INTEGER*4 IVALUE
      DATA      SPACE,TAB,COMMA,MINUS,NUM /1H ,"11,1H,,1H-,
     +          1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9/
 
      IA = 0
      LA = IABS(NI)
      IS = 1
  100 FOUND = .FALSE.
 
  110 IF (NEXT.GT.LREC) GOTO 180
      CHAR = REC(NEXT)
      NEXT = NEXT+1
 
      DO 120 I = 1,10
         IF (CHAR.EQ.NUM(I)) GOTO 140
  120 CONTINUE
      IF (CHAR.NE.MINUS) GOTO 130
      IS = -1
      IF (FOUND) GOTO 170
      GOTO 110
 
  130 IF (FOUND) GOTO 160
      IF (CHAR.NE.COMMA) GOTO 110
      IA = IA + 1
      IF (NI.GT.0) II(IA) = 0
      GOTO 170
 
  140 IF (FOUND) GOTO 150
      FOUND = .TRUE.
      IA = IA+1
      II(IA) = ISIGN(I-1,IS)
      GOTO 110
 
  150 IVALUE = 10*II(IA)+ISIGN(I-1,IS)
      IF (IVALUE.GT.32670 .OR. IVALUE.LT.-32670) THEN
         WRITE (6,1000)
         CALL GVSBEL
         RETURN
      END IF
      II(IA) = IVALUE
      GOTO 110
 
  160 IS = 1
  170 IF (CHAR.NE.SPACE .AND. CHAR.NE.COMMA) THEN
         CALL GVSBEL
         WRITE (6,1100)CHAR
      END IF
      IF (IA.LT.LA) GOTO 100
  180 GVFQQI = IA
      IF (IA.EQ.LA.OR.NI.LE.0) RETURN
      DO 190 I = IA+1,LA
         II(I) = 0
  190 CONTINUE
      RETURN
 
 1000 FORMAT (' TLS> Overflow of integer value ')
 1100 FORMAT (' TLS> Character is not a valid terminator (',A1,')',/)
 
      END
      INTEGER FUNCTION GVFQQR (LUN)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THE GVFQQ FUNCTIONS ARE USED TO READ TERMINAL OR FILE INPUT IN A FREE ----
C---- FORMAT.  GVFQQR READS A LINE FROM UNIT LUN AND STORES IT IN THE ARRAY ----
C---- REC. ROUTINES GVFQQA, GVFQQB, GVFQQC, GVFQQF AND GVFQQI THEN          ----
C---- INTERPRET THIS LINE AS ALPHANUMERIC, BYTE, CHARACTER, FLOATING POINT, ----
C---- OR INTEGER VALUES. THE FUNCTION RETURNS THE NUMBER OF CHARACTERS READ.----
C---- WARNING: THIS FUNCTION IS VAX-VMS SPECIFIC.                           ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      BYTE            REC(80), SPACE
      COMMON /GVCIO1/ LREC, NEXT, REC
      DATA            SPACE /1H /
C----
C---- READ THE INPUT LINE
C----
      READ  (LUN,1000,END=120,ERR=130) LREC,(REC(I),I=1,80)
C----
C---- DETERMINE THE LENGTH OF THE INPUT LINE
C----
      LREC = MIN0(LREC,79) + 1
C----
C---- CONVERT TO UPPER CASE
C----
      DO 100 I = 1,LREC
         IF (REC(I).GE. 97) REC(I) = REC(I)-32
  100 CONTINUE
C----
C---- PADD THE INNPUT LINE WITH BLANKS
C----
      DO 110 I = LREC,80
         REC(I) = SPACE
  110 CONTINUE
C----
C---- HERE UPON A NORMAL EXIT
C----
      NEXT = 1
      GVFQQR = LREC-1
      RETURN
C----
C---- HERE IF END OF INFORMATION WAS REACHED UPON READING
C----
  120 LREC = 0
      GVFQQR = -1
      RETURN
C----
C---- HERE UPON READING ERROR
C----
  130 WRITE (6,1100)
      LREC = 0
      GVFQQR = -1
      RETURN
 
 1000 FORMAT (Q,80A1)
 1100 FORMAT (' TLS> ERROR IN READING TERMINAL INPUT')
 
      END
      SUBROUTINE GVSCUC (TEXT)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- CONVERTS A LOWER CAST TEXT TO UPPER CAST                              ----
C---- NON CHARACTERS REMAIN UNALTERED                                       ----
C----                                                                       ----
C---- THIS ROUTINE IS BY NATURE OF ITS PURPOSE NOT WRITTEN IN FORTRAN 77.   ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      CHARACTER*(*) TEXT
      INTEGER       GVFLEN
      CHARACTER*1   UPPER(26)
      SAVE          UPPER
      DATA UPPER  /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     +             'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
C----
C---- INITIALIZE
C----
      LENTXT=GVFLEN(TEXT)
C----
C---- LOOP OVER THE LENGTH OF THE TEXT
C----
      IF (LENTXT.GT.0) THEN
         DO 20 L=1,LENTXT
C----------
C---------- FOR EVERY CHARACTER IN THE TEXT, LOOP OVER THE ALPHABET
C----------
            NUMLOW=INDEX('abcdefghijklmnopqrstuvwxyz',TEXT(L:L))
            IF (NUMLOW.NE.0) TEXT(L:L)=UPPER(NUMLOW)
   20    CONTINUE
      END IF
 
      RETURN
      END
      INTEGER FUNCTION GVFQQF(NF,FF)
C///
C-------------------------------------------------------------------------------
C----                                                                       ----
C---- THIS ROUTINE INTERPRETS THE DATA LINE STORED IN ARRAY REC AS A SERIES ----
C---- OF REALS*4. THE FUNCTION RETURNS THE NUMBER OF REALS ACTUALLY READ.   ----
C---- NA IS AN INTEGER.                                                     ----
C---- FF IS AN ARRAY OF REALS*4                                             ----
C---- WARNING: THIS FUNCTION IS VAX-VMS SPECIFIC.                           ----
C----                                                                       ----
C-------------------------------------------------------------------------------
C\\\
      COMMON    /GVCIO1/ LREC, NEXT, REC
      BYTE      REC(80)     ,NUM(10)     ,SPACE       ,TAB
      BYTE      COMMA       ,CHAR        ,MINUS       ,POINT
      BYTE      E
      LOGICAL   FOUND       ,FOUNDP      ,FOUNDE
      REAL*4    FF(NF)
      DATA      SPACE,TAB,COMMA,MINUS,POINT,E,NUM /1H ,"11,1H,,1H-,1H.,
     +          1HE,1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9/
C----
C---- INITIATE
C----
      IA = 0
      IS = 1
      LA = IABS(NF)
      IE = 0
  100 FOUND = .FALSE.
      FOUNDP = .FALSE.
      FOUNDE = .FALSE.
 
  110 IF (NEXT.GT.LREC) GOTO 240
 
      CHAR = REC(NEXT)
      NEXT = NEXT+1
 
      DO 120 I = 1,10
         IF (CHAR.EQ.NUM(I)) GOTO 180
  120 CONTINUE
      IF (CHAR.NE.MINUS) GOTO 130
      IS = -1
      IF (FOUND) GOTO 220
      GOTO 110
 
  130 IF (CHAR.NE.POINT) GOTO 170
      FOUNDP = .TRUE.
      GOTO 110
 
  160 FOUND = .FALSE.
      GOTO 110
 
  170 IF (FOUND) GOTO 210
      IF (CHAR.NE.COMMA) GOTO 110
      IA = IA+1
      IF (NF.GT.0) FF(IA)=0.0
      GOTO 230
 
  180 IF (FOUNDE) GOTO 200
      IF (FOUNDP) IE = IE-1
      IF (FOUND ) GOTO 190
      FOUND = .TRUE.
      IA = IA+1
      FF(IA) = ISIGN(I-1,IS)
      GOTO 110
 
  190 FF(IA) = 10.*FF(IA)+ISIGN(I-1,IS)
      GOTO 110
 
  200 IE = 10*IE + ISIGN(I-1,IS)
      FOUND = .TRUE.
      GOTO 110
 
  210 IS = 1
  220 IF (IE.EQ.0) GOTO 230
      FF(IA) = FF(IA)*10.**IE
      IE = 0
 
  230 IF (CHAR.NE.SPACE .AND. CHAR.NE.COMMA) THEN
         CALL GVSBEL
         WRITE (6,1000)CHAR
      END IF
      IF (IA.LT.LA) GOTO 100
 
  240 GVFQQF = IA
      IF (IA.EQ.LA.OR.NF.LE.0) RETURN
      DO 250 I = IA+1,LA
         FF(I) = 0.0
  250 CONTINUE
 
      RETURN
 
 1000 FORMAT(' TLS> CHARACTER IS NOT A VALID TERMINATOR (',A1,')',/)
 
      END

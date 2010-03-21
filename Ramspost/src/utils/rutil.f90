!############################# Change Log ##################################
! 4.3.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################
!
!     ******************************************************************
!
SUBROUTINE WINDUV(DD,FF,UU,VV)
DATA PI/3.14159/,PI180/.01745329/,V180PI/57.2957795/
UU=-FF*SIN(DD*PI180)
VV=-FF*COS(DD*PI180)
RETURN
ENTRY WINDDF(DD,FF,UU,VV)
U=UU
V=VV
FF=SQRT(U*U+V*V)
IF(ABS(U).LT.1.E-20)U=1.E-20
IF(ABS(V).LT.1.E-20)V=1.E-20
DD=ATAN2(-U,-V)*V180PI
IF(DD.LT.0.)DD=DD+360.
RETURN
END
!
!     **************************************************************
!
SUBROUTINE PRT2D(HORIZ,VERT,A,X,Y,IX,IY,I1,I2,J1,J2  &
                ,IFMT,TITLE,XX,YY,XLABL,YLABL)
CHARACTER*(*) IFMT,TITLE,XLABL,YLABL,HORIZ,VERT
CHARACTER*8 FMTX,FMTY
CHARACTER*133 B
CHARACTER*16 XLAB,YLAB,FMTB,FMTC
CHARACTER*80 FMTT,FMTT2
DIMENSION A(IX,IY),X(IX),Y(IY),XX(1),YY(1)
!
! for 15 columns      NCOL=105
NCOL=84
IFDW=7
MMAX=NCOL/IFDW
!                 See if the field is horizontally homogeneous (along
!                   the abscissa)
NDIFF=0
DO I=I1,I2
  XX(I-I1+1)=X(I)
  DO J=J1,J2
    YY(J-J1+1)=Y(J)
    A(I-I1+1,J-J1+1)=A(I,J)
    IF(A(I,J).NE.A(1,J-J1+1)) NDIFF=1
  ENDDO
ENDDO
XXMAX=MAX(ABS(X(I1)),ABS(X(I2)))
IF(XXMAX.LT.10.) FMTX='F7.3'
IF(XXMAX.GE.10.AND.XXMAX.LT.100.) FMTX='F7.2'
IF(XXMAX.GE.100.AND.XXMAX.LT.1000.) FMTX='F7.1'
IF(XXMAX.GE.1000.) FMTX='F7.0'
YYMAX=MAX(ABS(Y(J1)),ABS(Y(J2)))
IF(YYMAX.LT.10.) FMTY='F7.3'
IF(YYMAX.GE.10.AND.YYMAX.LT.100.) FMTY='F7.2'
IF(YYMAX.GE.100.AND.YYMAX.LT.1000.) FMTY='F7.1'
IF(YYMAX.GE.1000.) FMTY='F7.0'

!                  Set print window and number of pages accordingly
IF(NDIFF.EQ.0) THEN
  II1=1
  II2=1
  JJ1=1
  JJ2=J2-J1+1
  NPAGES=1
  NPPGE=MIN(MMAX,II2)
  IA=II1
  IB=II1
ELSE
  II1=1
  II2=I2-I1+1
  JJ1=1
  JJ2=J2-J1+1
  NPPGE=MIN(MMAX,II2)
  NPAGES=(II2-II1)/MMAX+1
  IA=II1
  IB=IA+NPPGE-1
ENDIF
!                      Set up formats for the coordinate line
IF(II1.EQ.II2)THEN
  FMTT='(/,1X,''CONSTANT ALONG ABSCISSA'',/)'
ELSE
  WRITE(FMTT,11) NPPGE,FMTX(1:4)
11      FORMAT('(/,1X,''PAGE : '',I3,//,2X,A1,''/'',A1,4X,',I3,A4,')')
  WRITE(FMTT2,12) NPPGE,FMTX(1:4)
12      FORMAT('(2X,A1,''/'',A1,4X,',I3,A4,')')
ENDIF
!
XLAB(1:8)=XLABL
YLAB(1:8)=YLABL
!                  Print top banner line
WRITE(6,220)
220 FORMAT(5(/),1X,'+',130('-'))
WRITE(6,21) TITLE,YLAB(1:8),XLAB(1:8)
21 FORMAT(' ! ',A68,'   Ordinate:',A8,'  Abscissa:',A8)
WRITE(6,221)
221 FORMAT(1X,'+',130('-'))
!
!                  Loop through number of pages
!
DO 1000 IPAGE=1,NPAGES
!
  IF(II1.EQ.II2)THEN
    WRITE(6,FMTT)
  ELSE
    WRITE(6,FMTT) IPAGE,VERT,HORIZ,(XX(II)+.00001,II=IA,IB)
  ENDIF
!
  B=' '
  WRITE(6,132)
  DO JJ=JJ1,JJ2
    JROW=JJ2-JJ+JJ1
    ICHST=1
    FMTC='(1X,'//FMTY//',A1)'
    WRITE(B(ICHST:ICHST+8),FMTC) YY(JROW)+.00001,'!'
    ICHST=ICHST+9
    DO II=IA,IB
      IF(A(II,JROW).EQ.0.)THEN
        B(ICHST:ICHST+IFDW-1)=' '
      ELSE
        FMTC='('//IFMT//')'
        WRITE(B(ICHST:ICHST+IFDW-1),FMTC) A(II,JROW)
      ENDIF
      ICHST=ICHST+IFDW
    ENDDO
    FMTC='(A1,'//FMTY//')'
    WRITE(B(ICHST:ICHST+7),FMTC) '!',YY(JROW)+.00001
    ICHEND=ICHST+7
!
    WRITE(FMTB,233) ICHEND
233     FORMAT('(A',I3,')')
    WRITE(6,FMTB) B(1:ICHEND)
!
  ENDDO
!
  WRITE(6,132)
132   FORMAT(8X,'+',105('-'),'+')
!
  IF(II1.NE.II2)THEN
    WRITE(6,FMTT2) VERT,HORIZ,(XX(II)+.00001,II=IA,IB)
  ENDIF
!
  IB=MIN(IB+NPPGE,II2)
  IA=IA+NPPGE
1000  CONTINUE
!
RETURN
END
!
!     *************************************************************
!
SUBROUTINE TRMGPH(A,L,TITLE,ITP,JTP)
DIMENSION ALAB(2)
CHARACTER*(*) TITLE
CHARACTER*8 ALAB,STR*1,BB*1,B*1
DIMENSION A(*),B(250),C(250)
DATA ALAB(1)/'   VALUE'/,ALAB(2) /'LOG10VAL'/
DATA STR/'-'/
IF(ITP.NE.1)THEN
  DO K=1,L
    IF(A(K).NE.0.)THEN
      C(K)=LOG10(ABS(A(K)))
    ELSE
      C(K)=-10
    ENDIF
  ENDDO
ELSE
  DO K=1,L
    C(K)=A(K)
  ENDDO
ENDIF
IMX=ISMAX(L,C,1)
IMN=ISMIN(L,C,1)
AMX=C(IMX)
AMN=C(IMN)
AINC=(AMX-AMN)*.1
PRINT 200,TITLE,C(IMX),IMX,C(IMN),IMN,ALAB(ITP),(K,K=5,L,5)
200   FORMAT(///,14X,A8,'   MAX IS',E17.10,' AT I =',I3,  &
'   MIN IS',E17.10,' AT I =',I3,  &
//,6X,A8,5X,20I5,//)
DO IL=1,11
  VAL=(11.-FLOAT(IL))*AINC+AMN
  DO K=1,L
    BB=' '
    IF(JTP.NE.1)THEN
      IF(C(K).GE.VAL) BB='*'
    ELSE
      IF(C(K).GE.VAL.AND.C(K).LT.VAL+AINC) BB='*'
    ENDIF
    B(K)=BB
  ENDDO
  PRINT 175,VAL,(B(K),K=1,L)
175     FORMAT(1X,E13.5,'  I  ',100A1)
ENDDO
PRINT 180,(STR,K=1,L)
180   FORMAT(19X,100A1)
RETURN
END

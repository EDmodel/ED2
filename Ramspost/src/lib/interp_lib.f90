!############################# Change Log ##################################
! 4.3.0.2
!
! 010328 JHC GDTOST2 ##
!            Added second horizontal interpolation scheme which is
!            a simple bilinear one to eliminate ringing errors. ##
! 001012 MJB HTINT ##
!            Declared passed args with length * instead of 1 to pass array
!            bounds compiler option for REVU. ##
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

SUBROUTINE TRNCL1(VCTRA,ZLEV,ZSURF,HTOP,VCTRB,VCTRC,VCTRD,NN)
DIMENSION VCTRA(1),VCTRB(1),VCTRC(1),VCTRD(1),ZLEV(1)

!     TRANSFORM TO Z COORDINATES FORM ZSTAR COORDINATES.

IF(ZSURF.EQ.0) RETURN
RTG=1.-ZSURF/HTOP
KABV=2
DO K=2,NN
VCTRC(K)=1./(RTG*(ZLEV(K)-ZLEV(K-1)))
VCTRD(K)=ZLEV(K)*RTG+ZSURF
ENDDO
VCTRD(1)=ZSURF
DO K=1,NN
KK1=KABV
DO KK=KK1,NN
KABV=KK
IF(VCTRD(KABV).GE.ZLEV(K)-1)GO TO 10
ENDDO
WRITE(6,5) ID,JD
5     FORMAT(' STOP IN TRNCL1',2I5)
STOP
10    CONTINUE
W1=MIN(1.,(VCTRD(KABV)-ZLEV(K))*VCTRC(KABV))
VCTRB(K)=VCTRA(KABV-1)*W1+VCTRA(KABV)*(1.-W1)
ENDDO
DO K=1,NN
VCTRA(K)=VCTRB(K)
ENDDO
RETURN
END

!     ******************************************************************

SUBROUTINE TRNCL2(VCTRA,ZLEV,ZSURF,HTOP,VCTRB,VCTRC,VCTRD,NN)
DIMENSION VCTRA(1),VCTRB(1),VCTRC(1),VCTRD(1),ZLEV(1)

!     TRANSFORM TO ZSTAR COORDINATES FROM Z COORDINATES

NNP=NN+1
IF(ZSURF.EQ.0) RETURN
RTG=1.-ZSURF/HTOP
KABV=2
DO K=2,NNP
VCTRC(K)=ZLEV(K)*RTG+ZSURF
VCTRD(K)=1./(ZLEV(K)-ZLEV(K-1))
ENDDO
VCTRC(1)=ZSURF
DO K=1,NN
KK1=KABV
DO KK=KK1,NNP
KABV=KK
IF(ZLEV(KABV)+.01.GE.VCTRC(K)) GO TO 10
ENDDO
WRITE(6,5) ID,JD
5     FORMAT(' STOP IN TRNCL2',2I5)
STOP
10    CONTINUE
W1=MIN((ZLEV(KABV)-VCTRC(K))*VCTRD(KABV),1.)
VCTRB(K)=W1*VCTRA(KABV-1)+(1.-W1)*VCTRA(KABV)
ENDDO
DO K=1,NN
VCTRA(K)=VCTRB(K)
ENDDO
RETURN
END

!     ******************************************************************

SUBROUTINE INTRP(A,B,ZVAL,Z,NZP)

!     LINEAR INTERPOLATION WITH Z

DIMENSION A(1),Z(1)
DO K=1,NZP
KABV=K
IF(ZVAL.LT.Z(K)) GO TO 10
ENDDO
10    CONTINUE
IF(KABV.EQ.1) KABV=2
B=(A(KABV)*(ZVAL-Z(KABV-1))+A(KABV-1)*(Z(KABV)-ZVAL))  &
/(Z(KABV)-Z(KABV-1))
RETURN
END

!     ******************************************************************

SUBROUTINE INTRRAP(A,B,PLNVAL,PLN,NZP)

!     LINEAR INTERPOLATION WITH LOG PRESSURE

DIMENSION PLN(1),A(1)
NZ=NZP-1
DO KK=1,NZ
K=NZP-KK
KABV=K+1
IF(PLNVAL.LT.PLN(K)) GO TO 10
ENDDO
10    CONTINUE
B=(A(KABV)*(PLNVAL-PLN(KABV-1))+A(KABV-1)*(PLN(KABV)-PLNVAL))  &
/(PLN(KABV)-PLN(KABV-1))
RETURN
END

!     ******************************************************************

SUBROUTINE BINOM(X1,X2,X3,X4,Y1,Y2,Y3,Y4,XXX,YYY)
COMMON/BIN/ITYPP,I0X,I1X,I2X,YOO
 YYY=1E30
 IF(X2.GT.1.E19.OR.X3.GT.1.E19.OR.  &
   Y2.GT.1.E19.OR.Y3.GT.1.E19)RETURN
WT1=(XXX-X3)/(X2-X3)
WT2=1.0-WT1
ISTEND=0
IF(Y4.LT.1.E19.AND.X4.LT.1.E19) GO TO 410
YZ22=WT1
YZ23=WT2
YZ24=0.0
ISTEND= 1
410   IF(Y1.LT.1.E19.AND.X1.LT.1.E19) GO TO 430
YZ11=0.0
YZ12=WT1
YZ13=WT2
IF(ISTEND.EQ.1)GO TO 480
GO TO 450
430   YZ11=(XXX-X2)*(XXX-X3)/((X1-X2)*(X1-X3))
YZ12=(XXX-X1)*(XXX-X3)/((X2-X1)*(X2-X3))
YZ13=(XXX-X1)*(XXX-X2)/((X3-X1)*(X3-X2))
IF(ISTEND.EQ.  1    ) GO TO 470
450   YZ22=(XXX-X3)*(XXX-X4)/((X2-X3)*(X2-X4))
YZ23=(XXX-X2)*(XXX-X4)/((X3-X2)*(X3-X4))
YZ24=(XXX-X2)*(XXX-X3)/((X4-X2)*(X4-X3))
470   YYY=WT1*(YZ11*Y1+YZ12*Y2+YZ13*Y3)+WT2*(YZ22*Y2+YZ23*Y3+YZ24*Y4)
 GO TO 490
480      YYY=WT1*Y2+WT2*Y3
490   YOO=YYY
RETURN
END

!     ******************************************************************

SUBROUTINE GDTOST(A,IX,IY,STAX,STAY,STAVAL)

!     SUBROUTINE TO RETURN STATIONS BACK-INTERPOLATED VALUES(STAVAL)
!     FROM UNIFORM GRID POINTS USING OVERLAPPING-QUADRATICS.
!     GRIDDED VALUES OF INPUT ARRAY A DIMENSIONED A(IX,IY),WHERE
!     IX=GRID POINTS IN X, IY = GRID POINTS IN Y .  STATION
!     LOCATION GIVEN IN TERMS OF GRID RELATIVE STATION X (STAX)
!     AND STATION COLUMN.
!     VALUES GREATER THAN 1.0E30 INDICATE MISSING DATA.

DIMENSION A(IX,IY),R(4),SCR(4)
IY1=INT(STAY)-1
IY2=IY1+3
IX1=INT(STAX)-1
IX2=IX1+3
STAVAL=1E30
FIYM2=FLOAT(IY1)-1
FIXM2=FLOAT(IX1)-1
II=0
DO 100 I=IX1,IX2
II=II+1
IF(I.GE.1.AND.I.LE.IX) GO TO 101
SCR(II)=1E30
GO TO 100
101   JJ=0
DO 111 J=IY1,IY2
JJ=JJ+1
IF(J.GE.1.AND.J.LE.IY) GO TO 112
R(JJ)=1E30
GO TO 111
112   R(JJ)=A(I,J)
111   CONTINUE
YY=STAY-FIYM2
CALL BINOM(1.,2.,3.,4.,R(1),R(2),R(3),R(4),YY,SCR(II))
100   CONTINUE
XX=STAX-FIXM2
CALL BINOM(1.,2.,3.,4.,SCR(1),SCR(2),SCR(3),SCR(4),XX,STAVAL)
RETURN
END

!     ******************************************************************

subroutine gdtost2(a,ix,iy,stax,stay,staval)

!     Subroutine to return stations back-interpolated values (staval)
!     from uniform grid points using bi-linear interpolation.
!     Gridded values of input array a dimensioned a(ix,iy), where
!     ix=grid points in x, iy = grid points in y.  Station
!     location given in terms of grid relative station x,y (stax,stay).

implicit none

! passed variables

integer :: ix, iy
real :: stax,stay,staval
real, dimension(ix,iy) :: a

! internal variables

integer :: i,j
real :: wtx1,wtx2,wty1,wty2


i = int(stax)
j = int(stay)
wtx2 = stax - float(i)
wty2 = stay - float(j)
wtx1 = 1. - wtx2
wty1 = 1. - wty2

staval = wtx1 * (wty1 * a(i  ,j  )+ wty2 * a(i  ,j+1)) &
       + wtx2 * (wty1 * a(i+1,j  )+ wty2 * a(i+1,j+1))
return
end subroutine gdtost2


!     ******************************************************************

SUBROUTINE HTINTCP(NZZ1,VCTRA,ELEVA,NZZ2,VCTRB,ELEVB,  &
 VT1C,VT2C,VT1D,VT2D,VT1E,VT2E)
DIMENSION VCTRA(1),VCTRB(1),ELEVA(1),ELEVB(1)
DIMENSION VT1C(1),VT1D(1),VT1E(1),VT2C(1),VT2D(1),VT2E(1)

L=1
DO 20 K=1,NZZ2
30 CONTINUE
IF(ELEVB(K).LT.ELEVA(1))GO TO 35
IF(ELEVB(K).GE.ELEVA(L).AND.ELEVB(K).LE.ELEVA(L+1))GO TO 35
L=L+1
IF(L.EQ.NZZ1)STOP 'htintcp'
GO TO 30
35 CONTINUE
WT=(ELEVB(K)-ELEVA(L))/(ELEVA(L+1)-ELEVA(L))
VCTRB(K)=VCTRA(L)+(VCTRA(L+1)-VCTRA(L))*WT
VT2C(K)=VT1C(L)+(VT1C(L+1)-VT1C(L))*WT
VT2D(K)=VT1D(L)+(VT1D(L+1)-VT1D(L))*WT
VT2E(K)=VT1E(L)+(VT1E(L+1)-VT1E(L))*WT
20 CONTINUE

RETURN
END

!     ******************************************************************

SUBROUTINE HTINT(NZZ1,VCTRA,ELEVA,NZZ2,VCTRB,ELEVB)
DIMENSION VCTRA(*),VCTRB(*),ELEVA(*),ELEVB(*)

L=1
DO 20 K=1,NZZ2
30 CONTINUE
IF(ELEVB(K).LT.ELEVA(1))GO TO 35
IF(ELEVB(K).GE.ELEVA(L).AND.ELEVB(K).LE.ELEVA(L+1))GO TO 35
IF(ELEVB(K).GT.ELEVA(NZZ1))GO TO 36
L=L+1
IF(L.EQ.NZZ1) then
  print *,'htint:nzz1',nzz1
  do kk=1,L
    print*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
  enddo
  stop 'htint'
endif
GO TO 30
35 CONTINUE
WT=(ELEVB(K)-ELEVA(L))/(ELEVA(L+1)-ELEVA(L))
VCTRB(K)=VCTRA(L)+(VCTRA(L+1)-VCTRA(L))*WT
GO TO 20
36 CONTINUE
WT=(ELEVB(K)-ELEVA(NZZ1))/(ELEVA(NZZ1-1)-ELEVA(NZZ1))
VCTRB(K)=VCTRA(NZZ1)+(VCTRA(NZZ1-1)-VCTRA(NZZ1))*WT
20 CONTINUE

RETURN
END

!     **************************************************************

SUBROUTINE HTINT2(NZZ1,VCTRA,ELEVA,NZZ2,VCTRB,ELEVB)
DIMENSION VCTRA(1),VCTRB(1),ELEVA(1),ELEVB(1)

!      htint for holding values of vctrb constant under eleva(1)

L=1
DO 20 K=1,NZZ2
30 CONTINUE
IF(ELEVB(K).LT.ELEVA(1))GO TO 34
IF(ELEVB(K).GE.ELEVA(L).AND.ELEVB(K).LE.ELEVA(L+1))GO TO 35
IF(ELEVB(K).GT.ELEVA(NZZ1))GO TO 36
L=L+1
IF(L.EQ.NZZ1)STOP 'htint2'
GO TO 30
34   CONTINUE
VCTRB(K)=VCTRA(1)
GO TO 20
35 CONTINUE
WT=(ELEVB(K)-ELEVA(L))/(ELEVA(L+1)-ELEVA(L))
VCTRB(K)=VCTRA(L)+(VCTRA(L+1)-VCTRA(L))*WT
GO TO 20
36 CONTINUE
WT=(ELEVB(K)-ELEVA(NZZ1))/(ELEVA(NZZ1-1)-ELEVA(NZZ1))
VCTRB(K)=VCTRA(NZZ1)+(VCTRA(NZZ1-1)-VCTRA(NZZ1))*WT
20 CONTINUE

RETURN
END

!     **************************************************************

SUBROUTINE AWTCMP(X,II,XX,IJ,ITYP,LLB,LRB,ICON,W,IORD)
DIMENSION W(6,6),X(1),XX(1),Q(6)

!     ITYP=0, NORMAL, LLB,LRB-LEFT,RIGHT BOUNDARIES
!          1, ISYM  , LLB-SYMMETRIC POINT,LRB-RIGHT BOUNDARY

DO L=1,6
   DO LL=1,6
      W(L,LL)=0.
   ENDDO
ENDDO

Q12=XX(IJ)

IF(ITYP.EQ.0.OR.(ITYP.EQ.1.AND.II.GE.LLB+2)) THEN
  IF(II-2.LT.LLB.OR.II+3.GT.LRB.OR.IORD.EQ.2) THEN
    W(1,3)=(X(II+1)-XX(IJ))/(X(II+1)-X(II))
    W(2,3)=.5/(X(II+1)-X(II))
    W(1,4)=-(X(II)-XX(IJ))/(X(II+1)-X(II))
    W(2,4)=-.5/(X(II+1)-X(II))
    RETURN
  ENDIF
  GOTO 3000
ENDIF

DO L=1,6
   Q(L)=X(II-3+L)
ENDDO

IF(ITYP.EQ.1.AND.II.LT.LLB+2) THEN
   Q(3)=X(II)
   Q(4)=X(II+1)
   Q(5)=X(II+2)
   IF(II.EQ.LLB) THEN
      Q(2)=2.*X(II)-X(II+1)
      Q(1)=2.*X(II)-X(II+2)
   ELSEIF(II.EQ.LLB+1) THEN
      Q(2)=X(II-1)
      Q(1)=2.*X(II)-X(II+1)
   ENDIF
ENDIF

3000 CONTINUE
GOTO (1000,2000) ICON+1
1000 CONTINUE

W(1,1)=  &
    ((Q(6)-Q12)*(Q(5)-Q12)*(Q(4)-Q12)*(Q(3)-Q12)*(Q(2)-Q12  &
 ))/((Q(6)-Q(1))*(Q(5)-Q(1))*(Q(4)-Q(1))*(Q(3)-Q(1))*(Q(2)  &
 -Q(1)))
W(2,1)=  &
    ((((Q(3)+Q(2)-2.*Q12)*Q(4)+(Q(2)-2.*Q12)*Q(3)-2.*Q(2)*  &
 Q12+3.*Q12**2)*Q(5)+((Q(2)-2.*Q12)*Q(3)-2.*Q(2)*Q12+3.*  &
 Q12**2)*Q(4)-(2.*Q(2)-3.*Q12)*Q(3)*Q12+3.*Q(2)*Q12**2-4.*  &
 Q12**3)*Q(6)+(((Q(2)-2.*Q12)*Q(3)-2.*Q(2)*Q12+3.*Q12**2)*  &
 Q(4)-(2.*Q(2)-3.*Q12)*Q(3)*Q12+3.*Q(2)*Q12**2-4.*Q12**3)*  &
 Q(5)-((2.*Q(2)-3.*Q12)*Q(3)-3.*Q(2)*Q12+4.*Q12**2)*Q(4)*  &
 Q12+(3.*Q(2)-4.*Q12)*Q(3)*Q12**2-4.*Q(2)*Q12**3+5.*Q12**4  &
 )/(2.*(Q(6)-Q(1))*(Q(5)-Q(1))*(Q(4)-Q(1))*(Q(3)-Q(1))*(Q(  &
 2)-Q(1)))
W(3,1)=  &
    (((Q(4)+Q(3)+Q(2)-3.*Q12)*Q(5)+(Q(3)+Q(2)-3.*Q12)*Q(4)  &
 +(Q(2)-3.*Q12)*Q(3)-3.*Q(2)*Q12+6.*Q12**2)*Q(6)+((Q(3)+Q(  &
 2)-3.*Q12)*Q(4)+(Q(2)-3.*Q12)*Q(3)-3.*Q(2)*Q12+6.*Q12**2)  &
 *Q(5)+((Q(2)-3.*Q12)*Q(3)-3.*Q(2)*Q12+6.*Q12**2)*Q(4)-3.*  &
 (Q(2)-2.*Q12)*Q(3)*Q12+6.*Q(2)*Q12**2-10.*Q12**3)/(3.*(Q(  &
 6)-Q(1))*(Q(5)-Q(1))*(Q(4)-Q(1))*(Q(3)-Q(1))*(Q(2)-Q(1)))
W(4,1)=  &
    ((Q(5)+Q(4)+Q(3)+Q(2)-4.*Q12)*Q(6)+(Q(4)+Q(3)+Q(2)-4.*  &
 Q12)*Q(5)+(Q(3)+Q(2)-4.*Q12)*Q(4)+(Q(2)-4.*Q12)*Q(3)-4.*Q  &
 (2)*Q12+10.*Q12**2)/(4.*(Q(6)-Q(1))*(Q(5)-Q(1))*(Q(4)-Q(1  &
 ))*(Q(3)-Q(1))*(Q(2)-Q(1)))
W(5,1)=  &
    (Q(6)+Q(5)+Q(4)+Q(3)+Q(2)-5.*Q12)/(5.*(Q(6)-Q(1))*(Q(5  &
 )-Q(1))*(Q(4)-Q(1))*(Q(3)-Q(1))*(Q(2)-Q(1)))
W(6,1)=  &
    1./(6.*(Q(6)-Q(1))*(Q(5)-Q(1))*(Q(4)-Q(1))*(Q(3)-Q(1))  &
 *(Q(2)-Q(1)))
W(1,2)=  &
    (-(Q(6)-Q12)*(Q(5)-Q12)*(Q(4)-Q12)*(Q(3)-Q12)*(Q(1)-  &
 Q12))/((Q(6)-Q(2))*(Q(5)-Q(2))*(Q(4)-Q(2))*(Q(3)-Q(2))*(Q  &
 (2)-Q(1)))
W(2,2)=  &
    (-((((Q(3)+Q(1)-2.*Q12)*Q(4)+(Q(1)-2.*Q12)*Q(3)-2.*Q(1  &
 )*Q12+3.*Q12**2)*Q(5)+((Q(1)-2.*Q12)*Q(3)-2.*Q(1)*Q12+3.*  &
 Q12**2)*Q(4)-(2.*Q(1)-3.*Q12)*Q(3)*Q12+3.*Q(1)*Q12**2-4.*  &
 Q12**3)*Q(6)+(((Q(1)-2.*Q12)*Q(3)-2.*Q(1)*Q12+3.*Q12**2)*  &
 Q(4)-(2.*Q(1)-3.*Q12)*Q(3)*Q12+3.*Q(1)*Q12**2-4.*Q12**3)*  &
 Q(5)-((2.*Q(1)-3.*Q12)*Q(3)-3.*Q(1)*Q12+4.*Q12**2)*Q(4)*  &
 Q12+(3.*Q(1)-4.*Q12)*Q(3)*Q12**2-4.*Q(1)*Q12**3+5.*Q12**4  &
 ))/(2.*(Q(6)-Q(2))*(Q(5)-Q(2))*(Q(4)-Q(2))*(Q(3)-Q(2))*(Q  &
 (2)-Q(1)))
W(3,2)=  &
    (-(((Q(4)+Q(3)+Q(1)-3.*Q12)*Q(5)+(Q(3)+Q(1)-3.*Q12)*Q(  &
 4)+(Q(1)-3.*Q12)*Q(3)-3.*Q(1)*Q12+6.*Q12**2)*Q(6)+((Q(3)+  &
 Q(1)-3.*Q12)*Q(4)+(Q(1)-3.*Q12)*Q(3)-3.*Q(1)*Q12+6.*Q12**  &
 2)*Q(5)+((Q(1)-3.*Q12)*Q(3)-3.*Q(1)*Q12+6.*Q12**2)*Q(4)-  &
 3.*(Q(1)-2.*Q12)*Q(3)*Q12+6.*Q(1)*Q12**2-10.*Q12**3))/(3.*  &
 (Q(6)-Q(2))*(Q(5)-Q(2))*(Q(4)-Q(2))*(Q(3)-Q(2))*(Q(2)-Q(1  &
 )))
W(4,2)=  &
    (-((Q(5)+Q(4)+Q(3)+Q(1)-4.*Q12)*Q(6)+(Q(4)+Q(3)+Q(1)-  &
 4.*Q12)*Q(5)+(Q(3)+Q(1)-4.*Q12)*Q(4)+(Q(1)-4.*Q12)*Q(3)-4.  &
 *Q(1)*Q12+10.*Q12**2))/(4.*(Q(6)-Q(2))*(Q(5)-Q(2))*(Q(4)-  &
 Q(2))*(Q(3)-Q(2))*(Q(2)-Q(1)))
W(5,2)=  &
    (-(Q(6)+Q(5)+Q(4)+Q(3)+Q(1)-5.*Q12))/(5.*(Q(6)-Q(2))*(  &
 Q(5)-Q(2))*(Q(4)-Q(2))*(Q(3)-Q(2))*(Q(2)-Q(1)))
W(6,2)=  &
    (-1.)/(6.*(Q(6)-Q(2))*(Q(5)-Q(2))*(Q(4)-Q(2))*(Q(3)-Q(  &
 2))*(Q(2)-Q(1)))
W(1,3)=  &
    ((Q(6)-Q12)*(Q(5)-Q12)*(Q(4)-Q12)*(Q(2)-Q12)*(Q(1)-Q12  &
 ))/((Q(6)-Q(3))*(Q(5)-Q(3))*(Q(4)-Q(3))*(Q(3)-Q(2))*(Q(3)  &
 -Q(1)))
W(2,3)=  &
    ((((Q(2)+Q(1)-2.*Q12)*Q(4)+(Q(1)-2.*Q12)*Q(2)-2.*Q(1)*  &
 Q12+3.*Q12**2)*Q(5)+((Q(1)-2.*Q12)*Q(2)-2.*Q(1)*Q12+3.*  &
 Q12**2)*Q(4)-(2.*Q(1)-3.*Q12)*Q(2)*Q12+3.*Q(1)*Q12**2-4.*  &
 Q12**3)*Q(6)+(((Q(1)-2.*Q12)*Q(2)-2.*Q(1)*Q12+3.*Q12**2)*  &
 Q(4)-(2.*Q(1)-3.*Q12)*Q(2)*Q12+3.*Q(1)*Q12**2-4.*Q12**3)*  &
 Q(5)-((2.*Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+4.*Q12**2)*Q(4)*  &
 Q12+(3.*Q(1)-4.*Q12)*Q(2)*Q12**2-4.*Q(1)*Q12**3+5.*Q12**4  &
 )/(2.*(Q(6)-Q(3))*(Q(5)-Q(3))*(Q(4)-Q(3))*(Q(3)-Q(2))*(Q(  &
 3)-Q(1)))
W(3,3)=  &
    (((Q(4)+Q(2)+Q(1)-3.*Q12)*Q(5)+(Q(2)+Q(1)-3.*Q12)*Q(4)  &
 +(Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+6.*Q12**2)*Q(6)+((Q(2)+Q(  &
 1)-3.*Q12)*Q(4)+(Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+6.*Q12**2)  &
 *Q(5)+((Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+6.*Q12**2)*Q(4)-3.*  &
 (Q(1)-2.*Q12)*Q(2)*Q12+6.*Q(1)*Q12**2-10.*Q12**3)/(3.*(Q(  &
 6)-Q(3))*(Q(5)-Q(3))*(Q(4)-Q(3))*(Q(3)-Q(2))*(Q(3)-Q(1)))
W(4,3)=  &
    ((Q(5)+Q(4)+Q(2)+Q(1)-4.*Q12)*Q(6)+(Q(4)+Q(2)+Q(1)-4.*  &
 Q12)*Q(5)+(Q(2)+Q(1)-4.*Q12)*Q(4)+(Q(1)-4.*Q12)*Q(2)-4.*Q  &
 (1)*Q12+10.*Q12**2)/(4.*(Q(6)-Q(3))*(Q(5)-Q(3))*(Q(4)-Q(3  &
 ))*(Q(3)-Q(2))*(Q(3)-Q(1)))
W(5,3)=  &
    (Q(6)+Q(5)+Q(4)+Q(2)+Q(1)-5.*Q12)/(5.*(Q(6)-Q(3))*(Q(5  &
 )-Q(3))*(Q(4)-Q(3))*(Q(3)-Q(2))*(Q(3)-Q(1)))
W(6,3)=  &
    1./(6.*(Q(6)-Q(3))*(Q(5)-Q(3))*(Q(4)-Q(3))*(Q(3)-Q(2))  &
 *(Q(3)-Q(1)))
W(1,4)=  &
    (-(Q(6)-Q12)*(Q(5)-Q12)*(Q(3)-Q12)*(Q(2)-Q12)*(Q(1)-  &
 Q12))/((Q(6)-Q(4))*(Q(5)-Q(4))*(Q(4)-Q(3))*(Q(4)-Q(2))*(Q  &
 (4)-Q(1)))
W(2,4)=  &
    (-((((Q(2)+Q(1)-2.*Q12)*Q(3)+(Q(1)-2.*Q12)*Q(2)-2.*Q(1  &
 )*Q12+3.*Q12**2)*Q(5)+((Q(1)-2.*Q12)*Q(2)-2.*Q(1)*Q12+3.*  &
 Q12**2)*Q(3)-(2.*Q(1)-3.*Q12)*Q(2)*Q12+3.*Q(1)*Q12**2-4.*  &
 Q12**3)*Q(6)+(((Q(1)-2.*Q12)*Q(2)-2.*Q(1)*Q12+3.*Q12**2)*  &
 Q(3)-(2.*Q(1)-3.*Q12)*Q(2)*Q12+3.*Q(1)*Q12**2-4.*Q12**3)*  &
 Q(5)-((2.*Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+4.*Q12**2)*Q(3)*  &
 Q12+(3.*Q(1)-4.*Q12)*Q(2)*Q12**2-4.*Q(1)*Q12**3+5.*Q12**4  &
 ))/(2.*(Q(6)-Q(4))*(Q(5)-Q(4))*(Q(4)-Q(3))*(Q(4)-Q(2))*(Q  &
 (4)-Q(1)))
W(3,4)=  &
    (-(((Q(3)+Q(2)+Q(1)-3.*Q12)*Q(5)+(Q(2)+Q(1)-3.*Q12)*Q(  &
 3)+(Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+6.*Q12**2)*Q(6)+((Q(2)+  &
 Q(1)-3.*Q12)*Q(3)+(Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+6.*Q12**  &
 2)*Q(5)+((Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+6.*Q12**2)*Q(3)-  &
 3.*(Q(1)-2.*Q12)*Q(2)*Q12+6.*Q(1)*Q12**2-10.*Q12**3))/(3.*  &
 (Q(6)-Q(4))*(Q(5)-Q(4))*(Q(4)-Q(3))*(Q(4)-Q(2))*(Q(4)-Q(1  &
 )))
W(4,4)=  &
    (-((Q(5)+Q(3)+Q(2)+Q(1)-4.*Q12)*Q(6)+(Q(3)+Q(2)+Q(1)-  &
 4.*Q12)*Q(5)+(Q(2)+Q(1)-4.*Q12)*Q(3)+(Q(1)-4.*Q12)*Q(2)-4.  &
 *Q(1)*Q12+10.*Q12**2))/(4.*(Q(6)-Q(4))*(Q(5)-Q(4))*(Q(4)-  &
 Q(3))*(Q(4)-Q(2))*(Q(4)-Q(1)))
W(5,4)=  &
    (-(Q(6)+Q(5)+Q(3)+Q(2)+Q(1)-5.*Q12))/(5.*(Q(6)-Q(4))*(  &
 Q(5)-Q(4))*(Q(4)-Q(3))*(Q(4)-Q(2))*(Q(4)-Q(1)))
W(6,4)=  &
    (-1.)/(6.*(Q(6)-Q(4))*(Q(5)-Q(4))*(Q(4)-Q(3))*(Q(4)-Q(  &
 2))*(Q(4)-Q(1)))
W(1,5)=  &
    ((Q(6)-Q12)*(Q(4)-Q12)*(Q(3)-Q12)*(Q(2)-Q12)*(Q(1)-Q12  &
 ))/((Q(6)-Q(5))*(Q(5)-Q(4))*(Q(5)-Q(3))*(Q(5)-Q(2))*(Q(5)  &
 -Q(1)))
W(2,5)=  &
    ((((Q(2)+Q(1)-2.*Q12)*Q(3)+(Q(1)-2.*Q12)*Q(2)-2.*Q(1)*  &
 Q12+3.*Q12**2)*Q(4)+((Q(1)-2.*Q12)*Q(2)-2.*Q(1)*Q12+3.*  &
 Q12**2)*Q(3)-(2.*Q(1)-3.*Q12)*Q(2)*Q12+3.*Q(1)*Q12**2-4.*  &
 Q12**3)*Q(6)+(((Q(1)-2.*Q12)*Q(2)-2.*Q(1)*Q12+3.*Q12**2)*  &
 Q(3)-(2.*Q(1)-3.*Q12)*Q(2)*Q12+3.*Q(1)*Q12**2-4.*Q12**3)*  &
 Q(4)-((2.*Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+4.*Q12**2)*Q(3)*  &
 Q12+(3.*Q(1)-4.*Q12)*Q(2)*Q12**2-4.*Q(1)*Q12**3+5.*Q12**4  &
 )/(2.*(Q(6)-Q(5))*(Q(5)-Q(4))*(Q(5)-Q(3))*(Q(5)-Q(2))*(Q(  &
 5)-Q(1)))
W(3,5)=  &
    (((Q(3)+Q(2)+Q(1)-3.*Q12)*Q(4)+(Q(2)+Q(1)-3.*Q12)*Q(3)  &
 +(Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+6.*Q12**2)*Q(6)+((Q(2)+Q(  &
 1)-3.*Q12)*Q(3)+(Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+6.*Q12**2)  &
 *Q(4)+((Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+6.*Q12**2)*Q(3)-3.*  &
 (Q(1)-2.*Q12)*Q(2)*Q12+6.*Q(1)*Q12**2-10.*Q12**3)/(3.*(Q(  &
 6)-Q(5))*(Q(5)-Q(4))*(Q(5)-Q(3))*(Q(5)-Q(2))*(Q(5)-Q(1)))
W(4,5)=  &
    ((Q(4)+Q(3)+Q(2)+Q(1)-4.*Q12)*Q(6)+(Q(3)+Q(2)+Q(1)-4.*  &
 Q12)*Q(4)+(Q(2)+Q(1)-4.*Q12)*Q(3)+(Q(1)-4.*Q12)*Q(2)-4.*Q  &
 (1)*Q12+10.*Q12**2)/(4.*(Q(6)-Q(5))*(Q(5)-Q(4))*(Q(5)-Q(3  &
 ))*(Q(5)-Q(2))*(Q(5)-Q(1)))
W(5,5)=  &
    (Q(6)+Q(4)+Q(3)+Q(2)+Q(1)-5.*Q12)/(5.*(Q(6)-Q(5))*(Q(5  &
 )-Q(4))*(Q(5)-Q(3))*(Q(5)-Q(2))*(Q(5)-Q(1)))
W(6,5)=  &
    1./(6.*(Q(6)-Q(5))*(Q(5)-Q(4))*(Q(5)-Q(3))*(Q(5)-Q(2))  &
 *(Q(5)-Q(1)))
W(1,6)=  &
    (-(Q(5)-Q12)*(Q(4)-Q12)*(Q(3)-Q12)*(Q(2)-Q12)*(Q(1)-  &
 Q12))/((Q(6)-Q(5))*(Q(6)-Q(4))*(Q(6)-Q(3))*(Q(6)-Q(2))*(Q  &
 (6)-Q(1)))
W(2,6)=  &
    (-((((Q(2)+Q(1)-2.*Q12)*Q(3)+(Q(1)-2.*Q12)*Q(2)-2.*Q(1  &
 )*Q12+3.*Q12**2)*Q(4)+((Q(1)-2.*Q12)*Q(2)-2.*Q(1)*Q12+3.*  &
 Q12**2)*Q(3)-(2.*Q(1)-3.*Q12)*Q(2)*Q12+3.*Q(1)*Q12**2-4.*  &
 Q12**3)*Q(5)+(((Q(1)-2.*Q12)*Q(2)-2.*Q(1)*Q12+3.*Q12**2)*  &
 Q(3)-(2.*Q(1)-3.*Q12)*Q(2)*Q12+3.*Q(1)*Q12**2-4.*Q12**3)*  &
 Q(4)-((2.*Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+4.*Q12**2)*Q(3)*  &
 Q12+(3.*Q(1)-4.*Q12)*Q(2)*Q12**2-4.*Q(1)*Q12**3+5.*Q12**4  &
 ))/(2.*(Q(6)-Q(5))*(Q(6)-Q(4))*(Q(6)-Q(3))*(Q(6)-Q(2))*(Q  &
 (6)-Q(1)))
W(3,6)=  &
    (-(((Q(3)+Q(2)+Q(1)-3.*Q12)*Q(4)+(Q(2)+Q(1)-3.*Q12)*Q(  &
 3)+(Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+6.*Q12**2)*Q(5)+((Q(2)+  &
 Q(1)-3.*Q12)*Q(3)+(Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+6.*Q12**  &
 2)*Q(4)+((Q(1)-3.*Q12)*Q(2)-3.*Q(1)*Q12+6.*Q12**2)*Q(3)-  &
 3.*(Q(1)-2.*Q12)*Q(2)*Q12+6.*Q(1)*Q12**2-10.*Q12**3))/(3.*  &
 (Q(6)-Q(5))*(Q(6)-Q(4))*(Q(6)-Q(3))*(Q(6)-Q(2))*(Q(6)-Q(1  &
 )))
W(4,6)=  &
    (-((Q(4)+Q(3)+Q(2)+Q(1)-4.*Q12)*Q(5)+(Q(3)+Q(2)+Q(1)-  &
 4.*Q12)*Q(4)+(Q(2)+Q(1)-4.*Q12)*Q(3)+(Q(1)-4.*Q12)*Q(2)-4.  &
 *Q(1)*Q12+10.*Q12**2))/(4.*(Q(6)-Q(5))*(Q(6)-Q(4))*(Q(6)-  &
 Q(3))*(Q(6)-Q(2))*(Q(6)-Q(1)))
W(5,6)=  &
    (-(Q(5)+Q(4)+Q(3)+Q(2)+Q(1)-5.*Q12))/(5.*(Q(6)-Q(5))*(  &
 Q(6)-Q(4))*(Q(6)-Q(3))*(Q(6)-Q(2))*(Q(6)-Q(1)))
W(6,6)=  &
    (-1.)/(6.*(Q(6)-Q(5))*(Q(6)-Q(4))*(Q(6)-Q(3))*(Q(6)-Q(  &
 2))*(Q(6)-Q(1)))

GOTO 240

2000 CONTINUE
W( 1, 1)= 1./60.
W( 2, 1)= 2./360.
W( 3, 1)=-1./48.
W( 4, 1)=-1./144.
W( 5, 1)= 1./240.
W( 6, 1)= 1./720.
W( 1, 2)=-8./60.
W( 2, 2)=-25./360.
W( 3, 2)= 7./48.
W( 4, 2)= 11./144.
W( 5, 2)=-3./240.
W( 6, 2)=-5./720.
W( 1, 3)=37./60.
W( 2, 3)=245./360.
W( 3, 3)=-6./48.
W( 4, 3)=-28/144.
W( 5, 3)= 2./240.
W( 6, 3)= 10./720.
W( 1, 4)=37./60.
W( 2, 4)=-245./360.
W( 3, 4)=-6./48.
W( 4, 4)= 28./144.
W( 5, 4)= 2./240.
W( 6, 4)=-10./720.
W( 1, 5)=-8./60.
W( 2, 5)= 25./360.
W( 3, 5)= 7./48.
W( 4, 5)=-11./144.
W( 5, 5)=-3./240.
W( 6, 5)= 5./720.
W( 1, 6)= 1./60.
W( 2, 6)=-2./360.
W( 3, 6)=-1./48.
W( 4, 6)= 1./144.
W( 5, 6)= 1./240.
W( 6, 6)=-1./720.
DO L=1,6
   DO LL=1,6
      W(LL,L)=W(LL,L)/(X(II+1)-X(II))**(LL-1)
   ENDDO
ENDDO

240 CONTINUE
IF(ITYP.EQ.1.AND.II.LT.LLB+2)THEN
IF(II.EQ.LLB)THEN
DO L=1,6
   W(L,4)=W(L,4)+W(L,2)
   W(L,5)=W(L,5)+W(L,1)
   W(L,2)=0.
   W(L,1)=0.
ENDDO
ELSEIF(II.EQ.LLB+1)THEN
DO L=1,6
   W(L,5)=W(L,5)+W(L,1)
   W(L,1)=0.
ENDDO
ENDIF
ENDIF

RETURN
END

!############################# Change Log ##################################
! 4.3.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

SUBROUTINE AZEROV(N1)
DIMENSION A1(*),A2(*),A3(*),A4(*),A5(*)
ENTRY AZERO(N1,A1)
   DO N=1,N1
      A1(N)=0.
   ENDDO
RETURN
ENTRY AZERO2(N1,A1,A2)
   DO N=1,N1
      A1(N)=0.
      A2(N)=0.
   ENDDO
RETURN
ENTRY AZERO3(N1,A1,A2,A3)
   DO N=1,N1
      A1(N)=0.
      A2(N)=0.
      A3(N)=0.
   ENDDO
RETURN
ENTRY AZERO4(N1,A1,A2,A3,A4)
   DO N=1,N1
      A1(N)=0.
      A2(N)=0.
      A3(N)=0.
      A4(N)=0.
   ENDDO
RETURN
ENTRY AZERO5(N1,A1,A2,A3,A4,A5)
   DO N=1,N1
      A1(N)=0.
      A2(N)=0.
      A3(N)=0.
      A4(N)=0.
      A5(N)=0.
   ENDDO
RETURN
END
!
SUBROUTINE AE1P1(NPTS,A,B,C)
DIMENSION A(NPTS),B(NPTS),C(NPTS)
DO I=1,NPTS
  A(I)=B(I)+C(I)
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE1M1(NPTS,A,B,C)
DIMENSION A(NPTS),B(NPTS),C(NPTS)
DO I=1,NPTS
  A(I)=B(I)-C(I)
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AEN1(NPTS,A,B)
DIMENSION A(NPTS),B(NPTS)
DO I=1,NPTS
  A(I)=-B(I)
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE1T1(NPTS,A,B,C)
DIMENSION A(NPTS),B(NPTS),C(NPTS)
DO I=1,NPTS
  A(I)=B(I)*C(I)
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE1(NPTS,A,B)
DIMENSION A(NPTS),B(NPTS)
DO I=1,NPTS
  A(I)=B(I)
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE1TN1(NPTS,A,B,C)
DIMENSION A(NPTS),B(NPTS),C(NPTS)
DO I=1,NPTS
  A(I)=-B(I)*C(I)
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE1T0(NPTS,A,B,C)
DIMENSION A(NPTS),B(NPTS)
DO I=1,NPTS
  A(I)=B(I)*C
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE1T0P1(NPTS,A,B,C,D)
DIMENSION A(NPTS),B(NPTS),D(NPTS)
DO I=1,NPTS
  A(I)=B(I)*C+D(I)
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE3(N1,N2,N3,I1,I2,J1,J2,K1,K2,A,B)
DIMENSION A(N1,N2,N3),B(N1,N2,N3)
DO J=J1,J2
  DO I=I1,I2
    DO K=K1,K2
      A(K,I,J)=B(K,I,J)
    ENDDO
  ENDDO
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE3P3(N1,N2,N3,I1,I2,J1,J2,K1,K2,A,B,C)
DIMENSION A(N1,N2,N3),B(N1,N2,N3),C(N1,N2,N3)
DO J=J1,J2
  DO I=I1,I2
    DO K=K1,K2
      A(K,I,J)=B(K,I,J)+C(K,I,J)
    ENDDO
  ENDDO
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE3M3(N1,N2,N3,I1,I2,J1,J2,K1,K2,A,B,C)
DIMENSION A(N1,N2,N3),B(N1,N2,N3),C(N1,N2,N3)
DO J=J1,J2
  DO I=I1,I2
    DO K=K1,K2
      A(K,I,J)=B(K,I,J)-C(K,I,J)
    ENDDO
  ENDDO
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE3T3(N1,N2,N3,I1,I2,J1,J2,K1,K2,A,B,C)
DIMENSION A(N1,N2,N3),B(N1,N2,N3),C(N1,N2,N3)
DO J=J1,J2
  DO I=I1,I2
    DO K=K1,K2
      A(K,I,J)=B(K,I,J)*C(K,I,J)
    ENDDO
  ENDDO
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE1P1P1(NPTS,A,B,C,F)
DIMENSION A(NPTS),B(NPTS),C(NPTS),F(NPTS)
DO I=1,NPTS
  A(I)=B(I)+C(I)+F(I)
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE1T1P1(NPTS,A,B,C,F)
DIMENSION A(NPTS),B(NPTS),C(NPTS),F(NPTS)
DO I=1,NPTS
  A(I)=B(I)*C(I)+F(I)
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE2(N2,N3,I1,I2,J1,J2,A,B)
DIMENSION A(N2,N3),B(N2,N3)
DO J=J1,J2
  DO I=I1,I2
    A(I,J)=B(I,J)
  ENDDO
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE3T3P3(N1,N2,N3,I1,I2,J1,J2,K1,K2,A,B,C,F)
DIMENSION A(N1,N2,N3),B(N1,N2,N3),C(N1,N2,N3),F(N1,N2,N3)
DO J=J1,J2
  DO I=I1,I2
    DO K=K1,K2
      A(K,I,J)=B(K,I,J)*C(K,I,J)+F(K,I,J)
    ENDDO
  ENDDO
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE3T0P3(N1,N2,N3,I1,I2,J1,J2,K1,K2,A,B,C,F)
DIMENSION A(N1,N2,N3),B(N1,N2,N3),F(N1,N2,N3)
DO J=J1,J2
  DO I=I1,I2
    DO K=K1,K2
      A(K,I,J)=B(K,I,J)*C+F(K,I,J)
    ENDDO
  ENDDO
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AEN3T0P3(N1,N2,N3,I1,I2,J1,J2,K1,K2,A,B,C,F)
DIMENSION A(N1,N2,N3),B(N1,N2,N3),F(N1,N2,N3)
DO J=J1,J2
  DO I=I1,I2
    DO K=K1,K2
      A(K,I,J)=-B(K,I,J)*C+F(K,I,J)
    ENDDO
  ENDDO
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE3M3D0(N1,N2,N3,I1,I2,J1,J2,K1,K2,A,B,C,F)
DIMENSION A(N1,N2,N3),B(N1,N2,N3),C(N1,N2,N3)
DO J=J1,J2
  DO I=I1,I2
    DO K=K1,K2
      A(K,I,J)=(B(K,I,J)-C(K,I,J))/F
    ENDDO
  ENDDO
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE A3E2(N1,N2,N3,I1,I2,J1,J2,K,A,B)
DIMENSION A(N1,N2,N3),B(N2,N3)
DO J=J1,J2
  DO I=I1,I2
    A(K,I,J)=B(I,J)
  ENDDO
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE A3E0(N1,N2,N3,I1,I2,J1,J2,K,A,B)
DIMENSION A(N1,N2,N3)
DO J=J1,J2
  DO I=I1,I2
    A(K,I,J)=B
  ENDDO
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE ALEBL(N1,N2,N3,KA,KB,A,B)
DIMENSION A(N1,N2,N3),B(N1,N2,N3)
DO J=1,N3
  DO I=1,N2
    A(KA,I,J)=B(KB,I,J)
  ENDDO
ENDDO
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE AE0(NPTS,A,B)
DIMENSION A(NPTS)
DO I=1,NPTS
  A(I)=B
ENDDO
RETURN
END
!
SUBROUTINE ADIVB(NNN,A,B,C)
DIMENSION A(1),B(1),C(1)
DO 1 NN=1,NNN
C(NN)=A(NN)/B(NN)
1 CONTINUE
RETURN
END
SUBROUTINE ATIMB(NNN,A,B,C)
DIMENSION A(1),B(1),C(1)
DO 1 NN=1,NNN
C(NN)=A(NN)*B(NN)
1 CONTINUE
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE TRID(VAR,CIM1,CI,CIP1,RHS,NPTS)
!
!     SOLVES A DIAGONALLY-DOMINANT TRIDIAGONAL MATRIX EQUATION BY
!     THE STANDARD QUICK METHOD
!
!     VAR   - VARIABLE BEING SOLVED FOR
!     CIM1  - VECTOR OF COEFFICIENTS AT THE I-1 POINT
!     CI    -   "     "       "       "  "  I     "
!     CIP1  -   "     "       "       "  "  I+1   "
!     RHS   -   "     "  THE RIGHT HAND SIDE OF THE EQUATION
!     NPTS  - NUMBER OF EQUATIONS IN THE MATRIX
!
!     WARNING: THE ARRAYS CIM1,CI, AND RHS ARE REUSED IN THIS ROUTINE.
!
DIMENSION VAR(1),CIM1(1),CI(1),CIP1(1),RHS(1)
!
CIP1(1)=CIP1(1)/CI(1)
DO 10 K=2,NPTS
CI(K)=CI(K)-CIM1(K)*CIP1(K-1)
CIP1(K)=CIP1(K)/CI(K)
10 CONTINUE
!
RHS(1)=RHS(1)/CI(1)
DO 20 K=2,NPTS
RHS(K)=(RHS(K)-CIM1(K)*RHS(K-1))/CI(K)
20 CONTINUE
!
VAR(NPTS)=RHS(NPTS)
DO 30 K=NPTS-1,1,-1
VAR(K)=RHS(K)-CIP1(K)*VAR(K+1)
30 CONTINUE
!
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE TRID2(VAR,CIM1,CI,CIP1,RHS,NPTS,SCR1,SCR2,SCR3)
!
!     SOLVES A DIAGONALLY-DOMINANT TRIDIAGONAL MATRIX EQUATION BY
!     THE STANDARD QUICK METHOD
!
!     VAR   - VARIABLE BEING SOLVED FOR
!     CIM1  - VECTOR OF COEFFICIENTS AT THE I-1 POINT
!     CI    -   "     "       "       "  "  I     "
!     CIP1  -   "     "       "       "  "  I+1   "
!     RHS   -   "     "  THE RIGHT HAND SIDE OF THE EQUATION
!     NPTS  - NUMBER OF EQUATIONS IN THE MATRIX
!     SCR1  - SCRATCH ARRAY AT LEAST NPTS LONG
!     SCR2  - SCRATCH ARRAY "    "    "    "
!     SCR3  - SCRATCH ARRAY "    "    "    "
!
DIMENSION VAR(1),CIM1(1),CI(1),CIP1(1),RHS(1)
DIMENSION SCR1(1),SCR2(1),SCR3(1)
!
SCR1(1)=CIP1(1)/CI(1)
SCR2(1)=CI(1)
DO 10 K=2,NPTS
SCR2(K)=CI(K)-CIM1(K)*SCR1(K-1)
SCR1(K)=CIP1(K)/SCR2(K)
10 CONTINUE
!
SCR3(1)=RHS(1)/SCR2(1)
DO 20 K=2,NPTS
SCR3(K)=(RHS(K)-CIM1(K)*SCR3(K-1))/SCR2(K)
20 CONTINUE
!
VAR(NPTS)=SCR3(NPTS)
DO 30 K=NPTS-1,1,-1
VAR(K)=SCR3(K)-SCR1(K)*VAR(K+1)
30 CONTINUE
!
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE UPDATE(N,A,FA,DT)
DIMENSION A(1),FA(1)
DO 10 NN=1,N
  A(NN)=A(NN)+FA(NN)*DT
10 CONTINUE
RETURN
END
!
!     ****************************************************************
!
SUBROUTINE ACCUM(NXYZ,ARR1,ARR2)
DIMENSION ARR1(1),ARR2(1)
DO N=1,NXYZ
  ARR1(N)=ARR1(N)+ARR2(N)
ENDDO
RETURN
END
!
!     ******************************************************************
!
!
SUBROUTINE ATOB(N,A,B)
DIMENSION A(1),B(1)
DO 100 I=1,N
B(I)=A(I)
100 CONTINUE
RETURN
END
!
!     ******************************************************************
!
SUBROUTINE ACNST(N,A,CNST)
DIMENSION A(1)
DO 10 NN=1,N
  A(NN)=CNST
10 CONTINUE
RETURN
END
!
!     ******************************************************************
!
FUNCTION VALUGP(N1,N2,N3,K,I,J,A)
DIMENSION A(N1,N2,N3)
  VALUGP=A(K,I,J)
RETURN
END
!
!     ******************************************************************
!
FUNCTION IVALUGP(N1,N2,N3,K,I,J,IA)
DIMENSION IA(N1,N2,N3)
  IVALUGP=IA(K,I,J)
RETURN
END

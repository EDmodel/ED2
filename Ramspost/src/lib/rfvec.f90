!############################# Change Log ##################################
! 4.3.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

!
!-----------------------------------------------------------------------
!        The following functions are the FORTRAN replacements for the
!          CRAY intrinsic vector functions that perform various tasks.
!          Since they are non-existent on
!          other machines, they need to be replaced by calls to that
!          machines's functions or simulated in standard FORTRAN.
!-----------------------------------------------------------------------
!
!       Return sum of vector.
!
FUNCTION SSUM(NN,VCTR,INC)
DIMENSION VCTR(1)
SUM=0.
NNN=NN*INC
DO 10 N=1,NNN,INC
  SUM=SUM+VCTR(N)
10 CONTINUE
SSUM=SUM
RETURN
END
! +------------------------------------------------------------------+
!
!       Return index of maximum of vector
!
FUNCTION ISMAX(NN,VCTR,INC)
DIMENSION VCTR(1)
ISM=0
SMAX=-1E10
DO 10 NNN=1,NN,INC
  IF(VCTR(NNN).GT.SMAX)THEN
    ISM=NNN
    SMAX=VCTR(NNN)
  ENDIF
10 CONTINUE
ISMAX=ISM
RETURN
END
! +------------------------------------------------------------------+
!
!       Return index of minimum of vector
!
FUNCTION ISMIN(NN,VCTR,INC)
DIMENSION VCTR(1)
ISM=0
SMIN=1E10
DO 10 NNN=1,NN,INC
  IF(VCTR(NNN).LT.SMIN)THEN
    ISM=NNN
    SMIN=VCTR(NNN)
  ENDIF
10 CONTINUE
ISMIN=ISM
RETURN
END
! +------------------------------------------------------------------+
!
!       Return VCT1 if VCT3 => 0., else VCT2.
!
FUNCTION CVMGP(VCT1,VCT2,VCT3)
IF(VCT3.GE.0)THEN
CVMGP=VCT1
ELSE
CVMGP=VCT2
ENDIF
RETURN
END
! +------------------------------------------------------------------+
!
!       Return VCT1 if VCT3 <= 0., else VCT2.
!
FUNCTION CVMGM(VCT1,VCT2,VCT3)
IF(VCT3.LT.0)THEN
CVMGM=VCT1
ELSE
CVMGM=VCT2
ENDIF
RETURN
END
! +------------------------------------------------------------------+
!
!       Return VCT1 if VCT3 = 0., else VCT2.
!
FUNCTION CVMGZ(VCT1,VCT2,VCT3)
IF(VCT3.EQ.0)THEN
CVMGZ=VCT1
ELSE
CVMGZ=VCT2
ENDIF
RETURN
END
! +------------------------------------------------------------------+
!
!       Return VCT1 if VCT3 NE 0., else VCT2.
!
FUNCTION CVMGN(VCT1,VCT2,VCT3)
IF(VCT3.NE.0)THEN
CVMGN=VCT1
ELSE
CVMGN=VCT2
ENDIF
RETURN
END

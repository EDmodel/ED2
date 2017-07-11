module shcu_vars_const

  use conv_coms, only : nkp      ! INTENT(IN)
!  INTEGER, PARAMETER :: nkp=100

  use grid_dims, only : maxgrds  ! INTENT(IN)

  ! Control Parameters. Some of the them in RAMSIN namelist
  integer :: nshcu, nnshcu(maxgrds)

!  COMMON/SHPARS /
  real    :: ENTF, ALHF, CAPE, EFIC, DCAPE, BCAPE, TCAPE, WSTAR
  integer :: KTOP

!  COMMON/SHPENV /
  real :: QVCON(NKP), AKVD(NKP), AKVDE(NKP), QVE(NKP)

!  COMMON/SHCON/
  real :: DSE(NKP), UHE(NKP), UHES(NKP), EVAPS(NKP), QVSE(NKP),  &
       RHE(NKP), GAMMA(NKP), DSC(NKP), WLC(NKP), UHC(NKP),       &
       DSCV(NKP), DSEV(NKP), DLDZBY2(NKP), DELZ(NKP),            &
       WC(NKP), VAPS(NKP), DTDT(NKP), DRDT(NKP), QVC(NKP),       &
       DQDT(NKP), DSC0(NKP), DSC0V(NKP), DSC0VM(NKP)

!  COMMON/ZI/
  real :: CL_CON(NKP), CL_PE(NKP)

!  COMMON/TMP/
  integer :: KZI

!  COMMON/SHCTES/
  real, parameter :: CPR=3.4965, CP=1004., P00=1E5, RCP=.286,  &
       ALVL3=2.5E6, ALIV=2.834E6, AKLV=2340.6, AKIV=2825.7, G=9.8, R=287.

end module shcu_vars_const

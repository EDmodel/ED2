! Module necessary to Grell Cumulus param.
! Variables used to define scratch array dimensions

module mem_grell_param

  !Memory for Grell's Cumulus scheme
  integer :: ngrids_cp, mgmxp, mgmyp, mgmzp
  integer, parameter ::               &
       maxiens = 2,  & !  Cloud spectral size
       maxens  = 3,                   & !3  ensemble one on cap_max
       maxens2 = 3,                   & !3  ensemble two on precip efficiency
       maxens_sh  = 3,                & !3  ensemble one on mbdt
       maxens2_sh = 1,                & !1  ensemble two on precip efficiency
       maxens3_sh = 10,               & ! 10 ensemble three done in cup_forcing_ens16
       maxens3 = 16                     !16 ensemble three done in cup_forcing_ens16
  integer :: ensdim,ensdim_sh                     !Ensemble dimension

  integer :: icoic                      ! Closure choice for deep
  integer :: icoic_sh                   ! Closure choice for shallow

  integer :: icbase  ! Choice of how to compute the PBL 
                     ! 1 - Maximum moist static energy
                     ! 2 - PBL top
  real, dimension(maxiens) :: depth_min ! Minimum depth that the cloud should have [m]
  real, dimension(maxiens) :: cap_maxs  ! Depth of inversion capping [mb]
  
  integer :: Flag_Grell = 0             ! = 0 Grell Arrays not allocated
  !     mgm*p not determined
  ! = 1 Grell Arrays not allocated
  !     mgm*p DETERMINED
  ! = 2 Grell Arrays ALLOCATED
  !     mgm*p DETERMINED

  character (len=2),dimension(maxiens) :: CLOSURE_TYPE  ! For new G.Grell Parameterization

contains

  subroutine define_memory(mmxp, mmyp, mmzp, ngrids, nnqparm, nnshcu)

    implicit none
    integer, dimension (*) ::   &
         mmxp,                  &  ! Number of points in X direction
         mmyp,                  &  ! Number of points in Y direction
         mmzp,                  &  ! Number of points in Z direction
         nnqparm,               &  ! Flag for cumulus parameterization
         nnshcu                    ! Flag for shallow cumulus parameterization

    ! indexed by number of grids
    ! The above integers data are passed by arguments and can be the amount
    ! of points in a node or the total points in a grid
    integer :: ngrids              ! Number of grids (nested)

    ! Local Variables
    integer :: i

    !Ensemble dimension
!srf - out 2004
    !ensdim=maxiens*maxens*maxens2*maxens3
     ensdim=1*maxens*maxens2*maxens3
     ! Allocate arrays based on options (if necessary)
     ![MLO - Shallow cumulus
     ensdim_sh=maxiens*maxens_sh*maxens2_sh*maxens3_sh
     !MLO
     !Memory for Grell's Cumulus Scheme
    mgmxp = 0
    mgmyp = 0
    mgmzp = 0
    ngrids_cp = 0
    do i=1, ngrids
       if (nnqparm(i) == 2 .or. nnshcu(i) == 2)  then
          mgmxp = max(mgmxp,mmxp(i))
          mgmyp = max(mgmyp,mmyp(i))
          mgmzp = max(mgmzp,mmzp(i))
          ngrids_cp = ngrids_cp + 1
       endif
    enddo

    Flag_Grell = 1  ! Seting the Flag

  end subroutine define_memory

end module mem_grell_param

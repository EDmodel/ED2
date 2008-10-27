! Module necessary to Grell Cumulus param.
! Variables used to define scratch array dimensions

module mem_grell_param

  use grell_coms, only: mgmzp
  !Memory for Grell's Cumulus scheme
  integer :: ngrids_cp, mgmxp, mgmyp
  integer, parameter ::                          &
       maxiens = 2,                   & !2 # of clouds with old grell (only 2 are allowed)
       maxens  = 3,                   & !3  ensemble one on cap_max
       maxens2 = 3,                   & !3  ensemble two on precip efficiency
       maxens_sh  = 3,                & !3  ensemble one on mbdt
       maxens2_sh = 1,                & !1  ensemble two on precip efficiency
       maxens3_sh = 10,               & ! 10 ensemble three done in cup_forcing_ens16
       maxens3 = 16                     !16 ensemble three done in cup_forcing_ens16
  integer :: ensdim,ensdim_sh                     !Ensemble dimension

  integer :: icoic                      ! Closure choice for deep
  integer :: icoic_sh                   ! Closure choice for shallow
  
  integer :: Flag_Grell = 0             ! = 0 Grell Arrays not allocated
  !     mgm*p not determined
  ! = 1 Grell Arrays not allocated
  !     mgm*p DETERMINED
  ! = 2 Grell Arrays ALLOCATED
  !     mgm*p DETERMINED


contains

  subroutine define_memory(mmxp, mmyp, mmzp, ngrids, nnqparm,ndeepest,nshallowest)

    implicit none
    integer, dimension (*) ::   &
         mmxp,                  &  ! Number of points in X direction
         mmyp,                  &  ! Number of points in Y direction
         mmzp,                  &  ! Number of points in Z direction
         nnqparm,               &  ! Flag for cumulus parameterization
         ndeepest,              &  ! Flag for deepest cumulus
         nshallowest               ! Flag for shallowest cumulus

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
       if (nnqparm(i) == 1 .and. (ndeepest(i) == 3 .or. nshallowest(i) == 3))  then
          mgmxp = max(mgmxp,mmxp(i))
          mgmyp = max(mgmyp,mmyp(i))
          mgmzp = max(mgmzp,mmzp(i))
          ngrids_cp = ngrids_cp + 1
       endif
    enddo

    Flag_Grell = 1  ! Seting the Flag

  end subroutine define_memory

end module mem_grell_param

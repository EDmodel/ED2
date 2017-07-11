module mem_tconv

  implicit none

  logical           :: trans_conv_alloc = .false.
  real, allocatable ::&
       se(:),         & ! environment scalar profile Z levels
       se_cup(:),     & ! environment scalar profile Z_cup levels
       sc_up(:),      & ! updraft   gas-phase  scalar profile
       sc_dn(:),      & ! DOwndraft gas-phase  scalar profile
       stcum1d(:),    & ! 1d convective tENDency
       dn01d(:),      & ! 1d air density
       sc_up_c(:),    & ! updraft   aqueous-phase scalar profile
       sc_dn_c(:),    & ! DOwndraft aqueous-phase scalar profile
       henry_coef(:), & ! Henry's constant for gases wet removal
       pw_up(:),      & ! updraft precitable gas/aer
       pw_dn(:)         ! DOwndraft precitable gas/aer

contains

  subroutine alloc_trans_conv(mgmzp)

    implicit none

    integer :: mgmzp

    if(trans_conv_alloc) then
       print *,'ERROR: trans_conv already allocated'
       print *,'Routine: trans_conv File: trans_conv.f90'
       print *,'Dir: .../shared/tools/brams20/src/rams/5.04/modules'
       stop
    end if

    allocate(se(mgmzp),se_cup(mgmzp),sc_up(mgmzp),sc_dn(mgmzp), &
         stcum1d(mgmzp),dn01d(mgmzp))

    allocate(sc_up_c(mgmzp),sc_dn_c(mgmzp),henry_coef(mgmzp), &
         pw_up(mgmzp),pw_dn(mgmzp) )

    trans_conv_alloc=.true.

  end subroutine alloc_trans_conv

  subroutine zero_tconv()

    implicit none

    se=.0
    se_cup=.0
    sc_up=.0
    sc_dn=.0
    stcum1d=.0
    dn01d=.0
    sc_up_c=.0
    sc_dn_c=.0
    henry_coef=.0
    pw_up=.0
    pw_dn=0.

  end subroutine zero_tconv

end module mem_tconv

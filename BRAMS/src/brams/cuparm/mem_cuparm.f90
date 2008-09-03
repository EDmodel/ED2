!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_cuparm

  use grid_dims

  type cuparm_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, dimension(:,:,:) :: &
          thsrc,rtsrc,thsrcf,rtsrcf,thsrcp,rtsrcp

     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer, dimension(:,:) :: &
          aconpr,conprr,conprrp,conprrf

  end type cuparm_vars

  type (cuparm_vars), allocatable :: cuparm_g(:), cuparmm_g(:)

  ! For CATT
  type (cuparm_vars), allocatable :: cuparm_g_sh(:), cuparmm_g_sh(:)

  integer, parameter :: maxcufiles=100, maxcugrids=10

  integer :: if_cuinv
  real(kind=8) :: tcu_beg, tcu_end
  real :: cu_til, cu_tel, tnudcu, wt_cu_grid(maxcugrids)
  character(len=128) :: cu_prefix

  character(len=128), dimension(maxcufiles) :: fnames_cu
  character(len=14) , dimension(maxcufiles) :: itotdate_cu
  real(kind=8), dimension(maxcufiles) :: cu_times

  integer :: ncufiles, ncufl
  real(kind=8) :: cutime1,cutime2

  integer, dimension(maxgrds) :: nnqparm
  real :: wcldbs,confrq

contains

  subroutine alloc_cuparm(cuparm,n1,n2,n3,ng)
![MLO - for shallow cumulus allocation to be correct
  use shcu_vars_const, only: nnshcu
!MLO]
    implicit none
    type (cuparm_vars) :: cuparm
    integer, intent(in) :: n1,n2,n3,ng

    ! Allocate arrays based on options (if necessary)

    if( nnqparm(ng)>= 0 .or. nnshcu(ng) >= 0 .or. if_cuinv == 1)  then
       allocate (cuparm%thsrc(n1,n2,n3))
       allocate (cuparm%rtsrc(n1,n2,n3))
       allocate (cuparm%aconpr(n2,n3))
       allocate (cuparm%conprr(n2,n3))
       if (if_cuinv == 1) then
          allocate (cuparm%thsrcp(n1,n2,n3))
          allocate (cuparm%rtsrcp(n1,n2,n3))
          allocate (cuparm%thsrcf(n1,n2,n3))
          allocate (cuparm%rtsrcf(n1,n2,n3))
          allocate (cuparm%conprrp(n2,n3))
          allocate (cuparm%conprrf(n2,n3))
       endif
    endif

    return
  end subroutine alloc_cuparm

  ! ----------------------------------------------------------------------

  subroutine nullify_cuparm(cuparm)

    implicit none
    type (cuparm_vars) :: cuparm

    if (associated(cuparm%thsrc))    nullify (cuparm%thsrc)
    if (associated(cuparm%rtsrc))    nullify (cuparm%rtsrc)
    if (associated(cuparm%thsrcp))    nullify (cuparm%thsrcp)
    if (associated(cuparm%rtsrcp))    nullify (cuparm%rtsrcp)
    if (associated(cuparm%thsrcf))    nullify (cuparm%thsrcf)
    if (associated(cuparm%rtsrcf))    nullify (cuparm%rtsrcf)
    if (associated(cuparm%aconpr))   nullify (cuparm%aconpr)
    if (associated(cuparm%conprr))   nullify (cuparm%conprr)
    if (associated(cuparm%conprrp))   nullify (cuparm%conprrp)
    if (associated(cuparm%conprrf))   nullify (cuparm%conprrf)

    return
  end subroutine nullify_cuparm

  subroutine dealloc_cuparm(cuparm)

    implicit none
    type (cuparm_vars) :: cuparm

    if (associated(cuparm%thsrc))    deallocate (cuparm%thsrc)
    if (associated(cuparm%rtsrc))    deallocate (cuparm%rtsrc)
    if (associated(cuparm%thsrcp))    deallocate (cuparm%thsrcp)
    if (associated(cuparm%rtsrcp))    deallocate (cuparm%rtsrcp)
    if (associated(cuparm%thsrcf))    deallocate (cuparm%thsrcf)
    if (associated(cuparm%rtsrcf))    deallocate (cuparm%rtsrcf)
    if (associated(cuparm%aconpr))   deallocate (cuparm%aconpr)
    if (associated(cuparm%conprr))   deallocate (cuparm%conprr)
    if (associated(cuparm%conprrp))   deallocate (cuparm%conprrp)
    if (associated(cuparm%conprrf))   deallocate (cuparm%conprrf)

    return
  end subroutine dealloc_cuparm

  ! ----------------------------------------------------------------------

  subroutine filltab_cuparm_sh(cuparm,cuparmm,imean,n1,n2,n3,ng)

    use var_tables

    implicit none
    type (cuparm_vars) :: cuparm,cuparmm
    integer, intent(in) :: imean,n1,n2,n3,ng
    integer :: npts

    ! Fill pointers to arrays into variable tables

    npts=n1*n2*n3

    if (associated(cuparm%thsrc))  &
         call vtables2 (cuparm%thsrc(1,1,1),cuparmm%thsrc(1,1,1)  &
         ,ng, npts, imean,  &
         'THSRC_SH :3:hist:anal:mpti:mpt3')
    if (associated(cuparm%rtsrc))  &
         call vtables2 (cuparm%rtsrc(1,1,1),cuparmm%rtsrc(1,1,1)  &
         ,ng, npts, imean,  &
         'RTSRC_SH :3:hist:anal:mpti:mpt3')
    if (associated(cuparm%thsrcp))  &
         call vtables2 (cuparm%thsrcp(1,1,1),cuparmm%thsrcp(1,1,1)  &
         ,ng, npts, imean,  &
         'THSRCP_SH :3:mpti:')
    if (associated(cuparm%rtsrcp))  &
         call vtables2 (cuparm%rtsrcp(1,1,1),cuparmm%rtsrcp(1,1,1)  &
         ,ng, npts, imean,  &
         'RTSRCP_SH :3:mpti:')
    if (associated(cuparm%thsrcf))  &
         call vtables2 (cuparm%thsrcf(1,1,1),cuparmm%thsrcf(1,1,1)  &
         ,ng, npts, imean,  &
         'THSRCF_SH :3:mpti:')
    if (associated(cuparm%rtsrcf))  &
         call vtables2 (cuparm%rtsrcf(1,1,1),cuparmm%rtsrcf(1,1,1)  &
         ,ng, npts, imean,  &
         'RTSRCF_SH :3:mpti:')

    npts=n2*n3
    if (associated(cuparm%aconpr))  &
         call vtables2 (cuparm%aconpr(1,1),cuparmm%aconpr(1,1)  &
         ,ng, npts, imean,  &
         'ACONPR_SH :2:hist:anal:mpti:mpt3')
    if (associated(cuparm%conprr))  &
         call vtables2 (cuparm%conprr(1,1),cuparmm%conprr(1,1)  &
         ,ng, npts, imean,  &
         'CONPRR_SH :2:hist:anal:mpt3')
    if (associated(cuparm%conprrp))  &
         call vtables2 (cuparm%conprrp(1,1),cuparmm%conprrp(1,1)  &
         ,ng, npts, imean,  &
         'CONPRRP_SH :2:mpti')
    if (associated(cuparm%conprrf))  &
         call vtables2 (cuparm%conprrf(1,1),cuparmm%conprrf(1,1)  &
         ,ng, npts, imean,  &
         'CONPRRF_SH :2:mpti')

    return
  end subroutine filltab_cuparm_sh

  ! ----------------------------------------------------------------------

  subroutine filltab_cuparm(cuparm,cuparmm,imean,n1,n2,n3,ng)

    use var_tables

    implicit none
    type (cuparm_vars) :: cuparm,cuparmm
    integer, intent(in) :: imean,n1,n2,n3,ng
    integer :: npts

    ! Fill pointers to arrays into variable tables

    npts=n1*n2*n3

    if (associated(cuparm%thsrc))  &
         call vtables2 (cuparm%thsrc(1,1,1),cuparmm%thsrc(1,1,1)  &
         ,ng, npts, imean,  &
         'THSRC :3:hist:anal:mpti:mpt3')
    if (associated(cuparm%rtsrc))  &
         call vtables2 (cuparm%rtsrc(1,1,1),cuparmm%rtsrc(1,1,1)  &
         ,ng, npts, imean,  &
         'RTSRC :3:hist:anal:mpti:mpt3')
    if (associated(cuparm%thsrcp))  &
         call vtables2 (cuparm%thsrcp(1,1,1),cuparmm%thsrcp(1,1,1)  &
         ,ng, npts, imean,  &
         'THSRCP :3:mpti:')
    if (associated(cuparm%rtsrcp))  &
         call vtables2 (cuparm%rtsrcp(1,1,1),cuparmm%rtsrcp(1,1,1)  &
         ,ng, npts, imean,  &
         'RTSRCP :3:mpti:')
    if (associated(cuparm%thsrcf))  &
         call vtables2 (cuparm%thsrcf(1,1,1),cuparmm%thsrcf(1,1,1)  &
         ,ng, npts, imean,  &
         'THSRCF :3:mpti:')
    if (associated(cuparm%rtsrcf))  &
         call vtables2 (cuparm%rtsrcf(1,1,1),cuparmm%rtsrcf(1,1,1)  &
         ,ng, npts, imean,  &
         'RTSRCF :3:mpti:')

    npts=n2*n3
    if (associated(cuparm%aconpr))  &
         call vtables2 (cuparm%aconpr(1,1),cuparmm%aconpr(1,1)  &
         ,ng, npts, imean,  &
         'ACONPR :2:hist:anal:mpti:mpt3')
    if (associated(cuparm%conprr))  &
         call vtables2 (cuparm%conprr(1,1),cuparmm%conprr(1,1)  &
         ,ng, npts, imean,  &
         'CONPRR :2:hist:anal:mpt3')
    if (associated(cuparm%conprrp))  &
         call vtables2 (cuparm%conprrp(1,1),cuparmm%conprrp(1,1)  &
         ,ng, npts, imean,  &
         'CONPRRP :2:mpti')
    if (associated(cuparm%conprrf))  &
         call vtables2 (cuparm%conprrf(1,1),cuparmm%conprrf(1,1)  &
         ,ng, npts, imean,  &
         'CONPRRF :2:mpti')

    return
  end subroutine filltab_cuparm

end module mem_cuparm

!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_varinit

  use grid_dims

  ! Memory for varfile, history, and condensate nudging

  type varinit_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, dimension(:,:,:) :: &
          varup,varvp,varpp,vartp,varrp,varop  &
          ,varuf,varvf,varpf,vartf,varrf,varof  &
          ,varwts &
          ,varrph,varrfh,varcph,varcfh

  end type varinit_vars

  type (varinit_vars), allocatable :: varinit_g(:), varinitm_g(:)


  integer, parameter :: maxnudfiles=44000  ! THIS SHOULD GET US TO ABOUT 5 YEARS

  integer :: nud_type, nnudfiles, nnudfl, nudlat

  character(len=str_len), dimension(maxnudfiles) :: fnames_nud
  character(len=14) , dimension(maxnudfiles) :: itotdate_nud
  real(kind=8), dimension(maxnudfiles) :: nud_times
  character(len=str_len) :: nud_hfile
  real :: tnudlat,tnudcent,tnudtop,znudtop
  real :: wt_nudge_uv,wt_nudge_th,wt_nudge_pi,wt_nudge_rt,wt_nudge_co2
  real :: wt_nudge_grid(maxgrds)
  real(kind=8) :: htime1, htime2

  integer :: igrid_match(maxgrds)

  !----------------------------------------------------------------------------
  integer :: nud_cond, ncondfiles, ncondfl
  character(len=str_len), dimension(maxnudfiles) :: fnames_cond
  character(len=14) , dimension(maxnudfiles) :: itotdate_cond
  real(kind=8), dimension(maxnudfiles) :: cond_times
  character(len=str_len) :: cond_hfile 
  real :: tcond_beg, tcond_end, wt_nudgec_grid(maxgrds),t_nudge_rc
  real(kind=8) :: condtime1, condtime2

  !----------------------------------------------------------------------------
  character(len=str_len), dimension(maxnudfiles) :: fnames_varf
  character(len=14) , dimension(maxnudfiles) :: itotdate_varf

  real(kind=8), dimension(maxnudfiles) :: varf_times

  character(len=str_len)       :: varfpfx
  ! Modif. by ALF

  real(kind=8)            :: vtime1,vtime2
  real                    :: vwait1,vwaittot
  integer                 :: nvarffiles,nvarffl

  character(len=14)       :: lastdate_iv

  !----------------------------------------------------------------------------

contains

  subroutine alloc_varinit(varinit,n1,n2,n3,ng,co2_on)

    implicit none
    type (varinit_vars) :: varinit
    integer, intent(in) :: n1,n2,n3,ng
    logical, intent(in) :: co2_on

    ! Allocate arrays based on options (if necessary)

    if( nud_type >= 1 ) then
       allocate (varinit%varup(n1,n2,n3))
       allocate (varinit%varvp(n1,n2,n3))
       allocate (varinit%varpp(n1,n2,n3))
       allocate (varinit%vartp(n1,n2,n3))
       allocate (varinit%varrp(n1,n2,n3))
       allocate (varinit%varuf(n1,n2,n3))
       allocate (varinit%varvf(n1,n2,n3))
       allocate (varinit%varpf(n1,n2,n3))
       allocate (varinit%vartf(n1,n2,n3))
       allocate (varinit%varrf(n1,n2,n3))                      
       allocate (varinit%varwts(n1,n2,n3))
       if (co2_on) then
          allocate (varinit%varop(n1,n2,n3))
          allocate (varinit%varof(n1,n2,n3))
       end if
    endif

    if (nud_cond == 1) then
       allocate (varinit%varcph(n1,n2,n3))
       allocate (varinit%varcfh(n1,n2,n3))                      
       allocate (varinit%varrph(n1,n2,n3))
       allocate (varinit%varrfh(n1,n2,n3))                      
    endif

    return
  end subroutine alloc_varinit


  subroutine nullify_varinit(varinit)

    implicit none
    type (varinit_vars) :: varinit


    nullify (varinit%varup)
    nullify (varinit%varvp)
    nullify (varinit%varpp)
    nullify (varinit%vartp)
    nullify (varinit%varrp)
    nullify (varinit%varop)
    nullify (varinit%varuf)
    nullify (varinit%varvf)
    nullify (varinit%varpf)
    nullify (varinit%vartf)
    nullify (varinit%varrf)
    nullify (varinit%varof)
    nullify (varinit%varwts)

    nullify (varinit%varcph)
    nullify (varinit%varcfh)
    nullify (varinit%varrph)
    nullify (varinit%varrfh)

    return
  end subroutine nullify_varinit

  subroutine dealloc_varinit(varinit)

    implicit none
    type (varinit_vars) :: varinit


    if (associated(varinit%varup))     deallocate (varinit%varup)
    if (associated(varinit%varvp))     deallocate (varinit%varvp)
    if (associated(varinit%varpp))     deallocate (varinit%varpp)
    if (associated(varinit%vartp))     deallocate (varinit%vartp)
    if (associated(varinit%varrp))     deallocate (varinit%varrp)
    if (associated(varinit%varop))     deallocate (varinit%varop)
    if (associated(varinit%varuf))     deallocate (varinit%varuf)
    if (associated(varinit%varvf))     deallocate (varinit%varvf)
    if (associated(varinit%varpf))     deallocate (varinit%varpf)
    if (associated(varinit%vartf))     deallocate (varinit%vartf)
    if (associated(varinit%varrf))     deallocate (varinit%varrf)
    if (associated(varinit%varof))     deallocate (varinit%varof)
    if (associated(varinit%varwts))    deallocate (varinit%varwts)

    if (associated(varinit%varcph))     deallocate (varinit%varcph)
    if (associated(varinit%varcfh))     deallocate (varinit%varcfh)
    if (associated(varinit%varrph))     deallocate (varinit%varrph)
    if (associated(varinit%varrfh))     deallocate (varinit%varrfh)

    return
  end subroutine dealloc_varinit


  subroutine filltab_varinit(varinit,varinitm,imean,n1,n2,n3,ng)

    use var_tables

    implicit none
    type (varinit_vars) :: varinit,varinitm
    integer, intent(in) :: imean,n1,n2,n3,ng
    integer :: npts
    real, pointer :: var,varm

    ! Fill pointers to arrays into variable tables

    npts=n1*n2*n3

    if (associated(varinit%varup))  &
         call vtables2 (varinit%varup,varinitm%varup  &
         ,ng, npts, imean,  &
         'VARUP :3:mpti')
    if (associated(varinit%varvp))  &
         call vtables2 (varinit%varvp,varinitm%varvp  &
         ,ng, npts, imean,  &
         'VARVP :3:mpti')
    if (associated(varinit%varpp))  &
         call vtables2 (varinit%varpp,varinitm%varpp  &
         ,ng, npts, imean,  &
         'VARPP :3:mpti')
    if (associated(varinit%vartp))  &
         call vtables2 (varinit%vartp,varinitm%vartp  &
         ,ng, npts, imean,  &
         'VARTP :3:mpti')
    if (associated(varinit%varrp))  &
         call vtables2 (varinit%varrp,varinitm%varrp  &
         ,ng, npts, imean,  &
         'VARRP :3:mpti')
    if (associated(varinit%varop))  &
         call vtables2 (varinit%varop,varinitm%varop  &
         ,ng, npts, imean,  &
         'VAROP :3:mpti')
    if (associated(varinit%varuf))  &
         call vtables2 (varinit%varuf,varinitm%varuf  &
         ,ng, npts, imean,  &
         'VARUF :3:mpti')
    if (associated(varinit%varvf))  &
         call vtables2 (varinit%varvf,varinitm%varvf  &
         ,ng, npts, imean,  &
         'VARVF :3:mpti')
    if (associated(varinit%varpf))  &
         call vtables2 (varinit%varpf,varinitm%varpf  &
         ,ng, npts, imean,  &
         'VARPF :3:mpti')
    if (associated(varinit%vartf))  &
         call vtables2 (varinit%vartf,varinitm%vartf  &
         ,ng, npts, imean,  &
         'VARTF :3:mpti')
    if (associated(varinit%varrf))  &
         call vtables2 (varinit%varrf,varinitm%varrf  &
         ,ng, npts, imean,  &
         'VARRF :3:mpti')
    if (associated(varinit%varof))  &
         call vtables2 (varinit%varof,varinitm%varof  &
         ,ng, npts, imean,  &
         'VAROF :3:mpti')
    if (associated(varinit%varwts))  &
         call vtables2 (varinit%varwts,varinitm%varwts  &
         ,ng, npts, imean,  &
         'VARWTS :3:mpti')

    if (nud_cond == 1) then               ! Inc. by ALF
       if (associated(varinit%varcph))  &
            call vtables2 (varinit%varcph,varinitm%varcph  &
            ,ng, npts, imean,  &
            'VARCPH :3:mpti')
       if (associated(varinit%varcfh))  &
            call vtables2 (varinit%varcfh,varinitm%varcfh  &
            ,ng, npts, imean,  &
            'VARCFH :3:mpti')
       if (associated(varinit%varrph))  &
            call vtables2 (varinit%varrph,varinitm%varrph  &
            ,ng, npts, imean,  &
            'VARRPH :3:mpti')
       if (associated(varinit%varrfh))  &
            call vtables2 (varinit%varrfh,varinitm%varrfh  &
            ,ng, npts, imean,  &
            'VARRFH :3:mpti')
    endif

    return
  end subroutine filltab_varinit

end module mem_varinit

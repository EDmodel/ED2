! Module necessary to GRELL param.

module mem_grell

  type grell_vars

     ! Variables to be dimensioned by (m2,m3)
     real, pointer, dimension(:,:) :: &
          UPMF,                       &
          DNMF,                       &
          XIACT_C,                    &
          XIACT_P,                    &
          XIERR,                      &
          XKDT,                       &
          XKTOP,                      &
          XJMIN,                      &
	  XK22,                       &		!For CATT
          XKBCON

     ! New variables for CATT
     ! Variables to be dimensioned by (m1,m2,m3)
     real, pointer, dimension(:,:,:) :: &
          lsfth, &
	  lsfrt

  end type grell_vars

  type (grell_vars), allocatable :: grell_g(:), grellm_g(:)
  !For CATT
  type (grell_vars), allocatable :: grell_g_sh(:), grellm_g_sh(:)

contains

  ! For CATT
  subroutine alloc_grell(grell, m1, m2, m3, ng)
    ! Modules necessary in V.5.04
    use mem_cuparm, only : nnqparm  ! INTENT(IN)
![MLO - For shallow cumulus without catt
    use shcu_vars_const, only : nnshcu  ! INTEN(IN)
!MLO]

    implicit none
    type (grell_vars) :: grell

    integer, intent(in) :: m1, m2, m3, ng


    !INCLUDE 'rcommons.h' ! Not necessary in V.5.04

    ! Allocate arrays based on options (if necessary)

    if( nnqparm(ng) == 2 .or. nnshcu(ng) == 2)  then
       allocate (grell%UPMF     (m2, m3))
       allocate (grell%DNMF     (m2, m3))
       allocate (grell%XIACT_C  (m2, m3))
       allocate (grell%XIACT_P  (m2, m3))
       allocate (grell%XIERR    (m2, m3))
       allocate (grell%XKDT     (m2, m3))
       allocate (grell%XKTOP    (m2, m3))
       allocate (grell%XKBCON   (m2, m3))
       allocate (grell%XJMIN    (m2, m3))
       !For CATT
       allocate (grell%XK22     (m2, m3))
       allocate (grell%lsfth(m1, m2, m3))
       allocate (grell%lsfrt(m1, m2, m3))
    endif

    return
  end subroutine alloc_grell

  subroutine nullify_grell(grell)

    implicit none
    type (grell_vars) :: grell

    if (associated(grell%UPMF))    nullify (grell%UPMF)
    if (associated(grell%DNMF))    nullify (grell%DNMF)
    if (associated(grell%XIACT_C)) nullify (grell%XIACT_C)
    if (associated(grell%XIACT_P)) nullify (grell%XIACT_P)
    if (associated(grell%XIERR))   nullify (grell%XIERR)
    if (associated(grell%XKDT))    nullify (grell%XKDT)
    if (associated(grell%XKTOP))   nullify (grell%XKTOP)
    if (associated(grell%XKBCON))  nullify (grell%XKBCON)
    if (associated(grell%XJMIN))   nullify (grell%XJMIN)
    !For CATT
    if (associated(grell%XK22))   nullify (grell%XK22)
    if (associated(grell%lsfth)) nullify (grell%lsfth)
    if (associated(grell%lsfrt)) nullify (grell%lsfrt)
    return
  end subroutine nullify_grell

  subroutine dealloc_grell(grell)

    implicit none
    type (grell_vars) :: grell

    if (associated(grell%UPMF))    deallocate (grell%UPMF)
    if (associated(grell%DNMF))    deallocate (grell%DNMF)
    if (associated(grell%XIACT_C)) deallocate (grell%XIACT_C)
    if (associated(grell%XIACT_P)) deallocate (grell%XIACT_P)
    if (associated(grell%XIERR))   deallocate (grell%XIERR)
    if (associated(grell%XKDT))    deallocate (grell%XKDT)
    if (associated(grell%XKTOP))   deallocate (grell%XKTOP)
    if (associated(grell%XKBCON))  deallocate (grell%XKBCON)
    if (associated(grell%XJMIN))   deallocate (grell%XJMIN)
    !For CATT
    if (associated(grell%XK22))   deallocate (grell%XK22)
    if (associated(grell%lsfth))   deallocate (grell%lsfth) !Lufla
    if (associated(grell%lsfrt))   deallocate (grell%lsfrt) !Lufla
    return
  end subroutine dealloc_grell

  ! For CATT
  subroutine filltab_grell(grell, grellm, imean, m1, m2, m3, ng)

    use var_tables

    !CATT
    use catt_start, only: CATT

    implicit none
    type (grell_vars) :: grell, grellm
    integer, intent(in) :: imean, m1,  m2, m3, ng
    integer :: npts

    ! Fill pointers to arrays into variable tables

    npts=m2*m3

    if (associated(grell%UPMF))  &
         call vtables2 (grell%UPMF(1,1),grellm%UPMF(1,1) &
         ,ng, npts, imean,  &
         'UPMF :2:hist:anal:mpti:mpt3')

    if (associated(grell%DNMF))  &
         call vtables2 (grell%DNMF(1,1),grellm%DNMF(1,1) &
         ,ng, npts, imean,  &
         'DNMF :2:hist:anal:mpti:mpt3')

    if (associated(grell%XIACT_C))  &
         call vtables2 (grell%XIACT_C(1,1),grellm%XIACT_C(1,1) &
         ,ng, npts, imean,  &
         'XIACT_C :2:hist:anal:mpti:mpt3')

    if (associated(grell%XIACT_P))  &
         call vtables2 (grell%XIACT_P(1,1),grellm%XIACT_P(1,1) &
         ,ng, npts, imean,  &
         'XIACT_P :2:hist:anal:mpti:mpt3')

    if (associated(grell%XIERR))  &
         call vtables2 (grell%XIERR(1,1),grellm%XIERR(1,1) &
         ,ng, npts, imean,  &
         'XIERR :2:hist:anal:mpti:mpt3')

    if (associated(grell%XKDT))  &
         call vtables2 (grell%XKDT(1,1),grellm%XKDT(1,1) &
         ,ng, npts, imean,  &
         'XKDT :2:hist:anal:mpti:mpt3')

    if (associated(grell%XKTOP))  &
         call vtables2 (grell%XKTOP(1,1),grellm%XKTOP(1,1) &
         ,ng, npts, imean,  &
         'XKTOP :2:hist:anal:mpti:mpt3')

    if (associated(grell%XKBCON))  &
         call vtables2 (grell%XKBCON(1,1),grellm%XKBCON(1,1) &
         ,ng, npts, imean,  &
         'XKBCON :2:hist:anal:mpti:mpt3')

    if (associated(grell%XJMIN))  &
         call vtables2 (grell%XJMIN(1,1),grellm%XJMIN(1,1) &
         ,ng, npts, imean,  &
         'XJMIN :2:hist:anal:mpti:mpt3')

    ! New variables for CATT
    if (CATT==1) then

       if (associated(grell%XK22))  &
            call vtables2 (grell%XK22(1,1),grellm%XK22(1,1) &
            ,ng, npts, imean,  &
            'XK22 :2:hist:anal:mpti:mpt3')

       ! 3D Arrays
       npts=m1*m2*m3

       if (associated(grell%lsfth))  &
            call vtables2 (grell%lsfth(1,1,1),grellm%lsfth(1,1,1) &
            ,ng, npts, imean,  &
            'lsfth :3:hist:anal:mpti:mpt3')

       if (associated(grell%lsfrt))  &
            call vtables2 (grell%lsfrt(1,1,1),grellm%lsfrt(1,1,1) &
            ,ng, npts, imean,  &
            'lsfrt :3:hist:anal:mpti:mpt3')

    else

       if (associated(grell%XK22))  &
            call vtables2 (grell%XK22(1,1),grellm%XK22(1,1) &
            ,ng, npts, imean,  &
            'XK22 :2:mpti:mpt3')

       ! 3D Arrays
       npts=m1*m2*m3

       if (associated(grell%lsfth))  &
            call vtables2 (grell%lsfth(1,1,1),grellm%lsfth(1,1,1) &
            ,ng, npts, imean,  &
            'lsfth :3:mpti:mpt3')

       if (associated(grell%lsfrt))  &
            call vtables2 (grell%lsfrt(1,1,1),grellm%lsfrt(1,1,1) &
            ,ng, npts, imean,  &
            'lsfrt :3:mpti:mpt3')

    endif

    return
  end subroutine filltab_grell


  ! For CATT
  subroutine filltab_grell_sh(grell, grellm, imean, m1, m2, m3, ng)

    use var_tables

    implicit none
    type (grell_vars) :: grell, grellm
    integer, intent(in) :: imean, m1,  m2, m3, ng
    integer :: npts

    ! Fill pointers to arrays into variable tables

    npts=m2*m3

    if (associated(grell%UPMF))  &
         call vtables2 (grell%UPMF(1,1),grellm%UPMF(1,1) &
         ,ng, npts, imean,  &
         'UPMFSH :2:hist:anal:mpti:mpt3')

    if (associated(grell%DNMF))  &
         call vtables2 (grell%DNMF(1,1),grellm%DNMF(1,1) &
         ,ng, npts, imean,  &
         'DNMFSH :2:hist:anal:mpti:mpt3')

    if (associated(grell%XIACT_C))  &
         call vtables2 (grell%XIACT_C(1,1),grellm%XIACT_C(1,1) &
         ,ng, npts, imean,  &
         'XIACT_CSH :2:hist:anal:mpti:mpt3')

    if (associated(grell%XIACT_P))  &
         call vtables2 (grell%XIACT_P(1,1),grellm%XIACT_P(1,1) &
         ,ng, npts, imean,  &
         'XIACT_PSH :2:hist:anal:mpti:mpt3')

    if (associated(grell%XIERR))  &
         call vtables2 (grell%XIERR(1,1),grellm%XIERR(1,1) &
         ,ng, npts, imean,  &
         'XIERRSH :2:hist:anal:mpti:mpt3')

    if (associated(grell%XKDT))  &
         call vtables2 (grell%XKDT(1,1),grellm%XKDT(1,1) &
         ,ng, npts, imean,  &
         'XKDTSH :2:hist:anal:mpti:mpt3')

    if (associated(grell%XKTOP))  &
         call vtables2 (grell%XKTOP(1,1),grellm%XKTOP(1,1) &
         ,ng, npts, imean,  &
         'XKTOPSH :2:hist:anal:mpti:mpt3')

    if (associated(grell%XKBCON))  &
         call vtables2 (grell%XKBCON(1,1),grellm%XKBCON(1,1) &
         ,ng, npts, imean,  &
         'XKBCONSH :2:hist:anal:mpti:mpt3')

    if (associated(grell%XJMIN))  &
         call vtables2 (grell%XJMIN(1,1),grellm%XJMIN(1,1) &
         ,ng, npts, imean,  &
         'XJMINSH :2:hist:anal:mpti:mpt3')

    ! New variables for CATT

    if (associated(grell%XK22))  &
         call vtables2 (grell%XK22(1,1),grellm%XK22(1,1) &
         ,ng, npts, imean,  &
         'XK22SH :2:hist:anal:mpti:mpt3')

    ! 3D Arrays
    npts=m1*m2*m3

    if (associated(grell%lsfth))  &
         call vtables2 (grell%lsfth(1,1,1),grellm%lsfth(1,1,1) &
         ,ng, npts, imean,  &
         'lsfthSH :3:hist:anal:mpti:mpt3')

    if (associated(grell%lsfrt))  &
         call vtables2 (grell%lsfrt(1,1,1),grellm%lsfrt(1,1,1) &
         ,ng, npts, imean,  &
         'lsfrtSH :3:hist:anal:mpti:mpt3')

    return
  end subroutine filltab_grell_sh

end module mem_grell

!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_turb

  use grid_dims

  type turb_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, dimension(:,:,:) :: &
          tkep,epsp,hkm,vkm,vkh,cdrag &
![MLO - For Nakanishi/Niino
         ,ltscale,sigw
!MLO]

     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer, dimension(:,:) :: &
          sflux_u,sflux_v,sflux_w,sflux_t,sflux_r,sflux_c &
![MLO - For Nakanishi/Niino
         ,lmo,pblhgt,akscal
     integer, pointer, dimension(:,:) :: kpbl
!MLO]

  end type turb_vars

  type (turb_vars), allocatable :: turb_g(:), turbm_g(:)

  integer :: if_urban_canopy,ihorgrad,ibruvais,ibotflx

  integer, dimension(maxgrds) :: idiffk

  real                     :: brunt ,rmax ,rmin
  real, dimension(maxgrds) :: zkhkm,xkhkm,csz,csx,akmin,akmax,hgtmin,hgtmax

contains
!==========================================================================================!
!==========================================================================================!
  subroutine alloc_turb(turb,n1,n2,n3,ng,co2_on)

    implicit none
    type (turb_vars) :: turb
    integer, intent(in) :: n1,n2,n3,ng
    logical, intent(in) :: co2_on

    ! Allocate arrays based on options (if necessary)
    
    ![MLO - Adding Nakanishi-Niino closure
    
    select case (idiffk(ng))
    case (1)
       allocate (turb%tkep(n1,n2,n3))
       allocate (turb%sigw(n1,n2,n3))
    case (4,5)
       allocate (turb%tkep(n1,n2,n3))
    case (6)
       allocate (turb%tkep(n1,n2,n3))
       allocate (turb%epsp(n1,n2,n3))
    case (7,8)
       allocate (turb%tkep(n1,n2,n3))
       allocate (turb%sigw(n1,n2,n3))
       allocate (turb%ltscale(n1,n2,n3))
       allocate (turb%pblhgt(n2,n3))
       allocate (turb%lmo(n2,n3))
    end select
    
    allocate (turb%kpbl(n2,n3))

    allocate (turb%hkm(n1,n2,n3))
    allocate (turb%vkm(n1,n2,n3))
    allocate (turb%vkh(n1,n2,n3))

    if(if_urban_canopy > 0) allocate (turb%cdrag(n1,n2,n3)) 

    allocate (turb%sflux_u(n2,n3))
    allocate (turb%sflux_v(n2,n3))
    allocate (turb%sflux_w(n2,n3))
    allocate (turb%sflux_t(n2,n3))
    allocate (turb%sflux_r(n2,n3))
    
    allocate (turb%akscal (n2,n3))

!    if (co2_on) allocate(turb%sflux_c(n2,n3))
    allocate(turb%sflux_c(n2,n3))

    return
  end subroutine alloc_turb
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine nullify_turb(turb)

    implicit none
    type (turb_vars) :: turb

    if (associated(turb%tkep    ))  nullify (turb%tkep    )
    if (associated(turb%epsp    ))  nullify (turb%epsp    )
    if (associated(turb%hkm     ))  nullify (turb%hkm     )
    if (associated(turb%vkm     ))  nullify (turb%vkm     )
    if (associated(turb%vkh     ))  nullify (turb%vkh     )
    if (associated(turb%cdrag   ))  nullify (turb%cdrag   )
    if (associated(turb%sflux_r ))  nullify (turb%sflux_r )
    if (associated(turb%sflux_u ))  nullify (turb%sflux_u )
    if (associated(turb%sflux_v ))  nullify (turb%sflux_v )
    if (associated(turb%sflux_w ))  nullify (turb%sflux_w )
    if (associated(turb%sflux_t ))  nullify (turb%sflux_t )
    if (associated(turb%sflux_c ))  nullify (turb%sflux_c )
    if (associated(turb%akscal  ))  nullify (turb%akscal  )
    if (associated(turb%ltscale ))  nullify (turb%ltscale )
    if (associated(turb%sigw    ))  nullify (turb%sigw    )
    if (associated(turb%pblhgt  ))  nullify (turb%pblhgt  )
    if (associated(turb%lmo     ))  nullify (turb%lmo     )
    if (associated(turb%kpbl    ))  nullify (turb%kpbl    )

    return
  end subroutine nullify_turb
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine dealloc_turb(turb)

    implicit none
    type (turb_vars) :: turb


    if (associated(turb%tkep    ))  deallocate (turb%tkep    )
    if (associated(turb%epsp    ))  deallocate (turb%epsp    )
    if (associated(turb%hkm     ))  deallocate (turb%hkm     )
    if (associated(turb%vkm     ))  deallocate (turb%vkm     )
    if (associated(turb%vkh     ))  deallocate (turb%vkh     )
    if (associated(turb%cdrag   ))  deallocate (turb%cdrag   )
    if (associated(turb%sflux_r ))  deallocate (turb%sflux_r )
    if (associated(turb%sflux_u ))  deallocate (turb%sflux_u )
    if (associated(turb%sflux_v ))  deallocate (turb%sflux_v )
    if (associated(turb%sflux_w ))  deallocate (turb%sflux_w )
    if (associated(turb%sflux_t ))  deallocate (turb%sflux_t )
    if (associated(turb%sflux_c ))  deallocate (turb%sflux_c )
    if (associated(turb%akscal  ))  deallocate (turb%akscal  )
    if (associated(turb%ltscale ))  deallocate (turb%ltscale )
    if (associated(turb%sigw    ))  deallocate (turb%sigw    )
    if (associated(turb%pblhgt  ))  deallocate (turb%pblhgt  )
    if (associated(turb%lmo     ))  deallocate (turb%lmo     )
    if (associated(turb%kpbl    ))  deallocate (turb%kpbl    )

    return
  end subroutine dealloc_turb
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine zero_turb(turb)
    use rconstants, only: tkmin
    implicit none
    type (turb_vars) :: turb


    if (associated(turb%tkep    ))  turb%tkep    = tkmin
    if (associated(turb%epsp    ))  turb%epsp    = tkmin
    if (associated(turb%hkm     ))  turb%hkm     = 0.
    if (associated(turb%vkm     ))  turb%vkm     = 0.
    if (associated(turb%vkh     ))  turb%vkh     = 0.
    if (associated(turb%cdrag   ))  turb%cdrag   = 0.
    if (associated(turb%sflux_u ))  turb%sflux_u = 0.
    if (associated(turb%sflux_v ))  turb%sflux_v = 0.
    if (associated(turb%sflux_w ))  turb%sflux_w = 0.
    if (associated(turb%sflux_t ))  turb%sflux_t = 0.
    if (associated(turb%sflux_r ))  turb%sflux_r = 0.
    if (associated(turb%sflux_c ))  turb%sflux_c = 0.
    if (associated(turb%akscal  ))  turb%akscal  = 0.
    if (associated(turb%ltscale ))  turb%ltscale = 0.
    if (associated(turb%sigw    ))  turb%sigw    = 0.
    if (associated(turb%pblhgt  ))  turb%pblhgt  = 0.
    if (associated(turb%lmo     ))  turb%lmo     = 0.
    if (associated(turb%kpbl    ))  turb%kpbl    = 0

    return
  end subroutine zero_turb
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine filltab_turb(turb,turbm,imean,n1,n2,n3,ng)

    use var_tables

    implicit none
    type (turb_vars) :: turb,turbm
    integer, intent(in) :: imean,n1,n2,n3,ng
    integer :: npts
    real, pointer :: var,varm

    ! Fill pointers to arrays into variable tables

    npts=n1*n2*n3

    if (associated(turb%tkep))  &
         call vtables2 (turb%tkep(1,1,1),turbm%tkep(1,1,1)  &
         ,ng, npts, imean,  &
         'TKEP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(turb%epsp))  &
         call vtables2 (turb%epsp(1,1,1),turbm%epsp(1,1,1)  &
         ,ng, npts, imean,  &
         'EPSP :3:hist:anal:mpti:mpt3:mpt1')

    if (associated(turb%hkm))  &
         call vtables2 (turb%hkm(1,1,1),turbm%hkm(1,1,1)  &
         ,ng, npts, imean,  &
         'HKM :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(turb%vkm))  &
         call vtables2 (turb%vkm(1,1,1),turbm%vkm(1,1,1)  &
         ,ng, npts, imean,  &
         'VKM :3:hist:mpti:mpt3:mpt1')
    if (associated(turb%vkh))  &
         call vtables2 (turb%vkh(1,1,1),turbm%vkh(1,1,1)  &
         ,ng, npts, imean,  &
         'VKH :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(turb%cdrag))  &
         call vtables2 (turb%cdrag(1,1,1),turbm%cdrag(1,1,1)  &
         ,ng, npts, imean,  &
         'CDRAG :3:hist:anal:mpti')
    if (associated(turb%ltscale)) &
       call vtables2 (turb%ltscale(1,1,1),turbm%ltscale(1,1,1) &
         ,ng, npts, imean, &
         'TL :3:hist:anal:mpti:mpt3')
    if (associated(turb%sigw)) &
       call vtables2 (turb%sigw(1,1,1),turbm%sigw(1,1,1) &
         ,ng, npts, imean, &
         'SIGW :3:hist:anal:mpti:mpt3')

    npts=n2*n3
    if (associated(turb%sflux_u))  &
         call vtables2 (turb%sflux_u(1,1),turbm%sflux_u(1,1)  &
         ,ng, npts, imean,  &
         'SFLUX_U :2:anal:mpt3:mpt1')
    if (associated(turb%sflux_v))  &
         call vtables2 (turb%sflux_v(1,1),turbm%sflux_v(1,1)  &
         ,ng, npts, imean,  &
         'SFLUX_V :2:anal:mpt3:mpt1')
    if (associated(turb%sflux_w))  &
         call vtables2 (turb%sflux_w(1,1),turbm%sflux_w(1,1)  &
         ,ng, npts, imean,  &
         'SFLUX_W :2:anal:mpt3')
    if (associated(turb%sflux_t))  &
         call vtables2 (turb%sflux_t(1,1),turbm%sflux_t(1,1)  &
         ,ng, npts, imean,  &
         'SFLUX_T :2:anal:mpt3')
    if (associated(turb%sflux_r))  &
         call vtables2 (turb%sflux_r(1,1),turbm%sflux_r(1,1)  &
         ,ng, npts, imean,  &
         'SFLUX_R :2:anal:mpt3')
    if (associated(turb%sflux_c))  &
         call vtables2 (turb%sflux_c(1,1),turbm%sflux_c(1,1)  &
         ,ng, npts, imean,  &
         'SFLUX_C :2:anal:mpt3')
    if (associated(turb%akscal))  &
         call vtables2 (turb%akscal(1,1),turbm%akscal(1,1)  &
         ,ng, npts, imean,  &
         'AKSCAL :2:hist:anal:mpti:mpt3')
    if (associated(turb%pblhgt)) &
       call vtables2 (turb%pblhgt(1,1),turbm%pblhgt(1,1) &
         ,ng, npts, imean, &
         'PBLHGT :2:hist:anal:mpti:mpt3')
    if (associated(turb%lmo)) &
       call vtables2 (turb%lmo(1,1),turbm%lmo(1,1) &
         ,ng, npts, imean, &
         'LMO    :2:hist:anal:mpti:mpt3')

    return
  end subroutine filltab_turb
!==========================================================================================!
!==========================================================================================!
end module mem_turb

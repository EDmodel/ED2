!============================= Change Log =================================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
module mem_turb

   use grid_dims, only : maxgrds ! ! intent(in)

   type turb_vars
      !----- Variables to be dimensioned by (nzp,nxp,nyp). --------------------------------!
      real   , dimension(:,:,:), pointer :: tkep
      real   , dimension(:,:,:), pointer :: epsp
      real   , dimension(:,:,:), pointer :: hkm
      real   , dimension(:,:,:), pointer :: vkm
      real   , dimension(:,:,:), pointer :: vkh
      real   , dimension(:,:,:), pointer :: cdrag
      real   , dimension(:,:,:), pointer :: ltscale
      real   , dimension(:,:,:), pointer :: sigw
      !----- Variables to be dimensioned by (nxp,nyp). ------------------------------------!
      real   , dimension(:,:)  , pointer :: sflux_u
      real   , dimension(:,:)  , pointer :: sflux_v
      real   , dimension(:,:)  , pointer :: sflux_w
      real   , dimension(:,:)  , pointer :: sflux_t
      real   , dimension(:,:)  , pointer :: sflux_r
      real   , dimension(:,:)  , pointer :: sflux_c
      real   , dimension(:,:)  , pointer :: lmo
      real   , dimension(:,:)  , pointer :: pblhgt
      real   , dimension(:,:)  , pointer :: akscal
      integer, dimension(:,:)  , pointer :: kpbl
      !------------------------------------------------------------------------------------!
   end type turb_vars

   !----- The main structures (instantaneous and averages). -------------------------------!
   type(turb_vars), dimension(:), allocatable :: turb_g
   type(turb_vars), dimension(:), allocatable :: turbm_g
   !---------------------------------------------------------------------------------------!


   !----- Variables to be filled by namelist variables. -----------------------------------!
   integer                                    :: if_urban_canopy
   integer                                    :: ihorgrad
   integer                                    :: ibruvais
   integer                                    :: ibotflx
   integer, dimension(maxgrds)                :: idiffk
   real   , dimension(maxgrds)                :: zkhkm
   real   , dimension(maxgrds)                :: xkhkm
   real   , dimension(maxgrds)                :: csz
   real   , dimension(maxgrds)                :: csx
   real   , dimension(maxgrds)                :: akmin
   real   , dimension(maxgrds)                :: akmax
   real   , dimension(maxgrds)                :: hgtmin
   real   , dimension(maxgrds)                :: hgtmax
   !---------------------------------------------------------------------------------------!


   !----- Other variables. ----------------------------------------------------------------!
   real                                       :: brunt
   real                                       :: rmax
   real                                       :: rmin
   !---------------------------------------------------------------------------------------!



   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine allocates the turbulence structure.                              !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_turb(turb,n1,n2,n3,ng)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (turb_vars), intent(inout) :: turb
      integer         , intent(in)    :: n1
      integer         , intent(in)    :: n2
      integer         , intent(in)    :: n3
      integer         , intent(in)    :: ng
      !------------------------------------------------------------------------------------!


      !----- Allocate arrays based on options (if necessary). -----------------------------!
      select case (idiffk(ng))
      case (1)
         allocate (turb%tkep    (n1,n2,n3))
         allocate (turb%sigw    (n1,n2,n3))
      case (4,5)
         allocate (turb%tkep    (n1,n2,n3))
      case (6)
         allocate (turb%tkep    (n1,n2,n3))
         allocate (turb%epsp    (n1,n2,n3))
      case (7,8)
         allocate (turb%tkep    (n1,n2,n3))
         allocate (turb%sigw    (n1,n2,n3))
         allocate (turb%ltscale (n1,n2,n3))
         allocate (turb%pblhgt     (n2,n3))
         allocate (turb%lmo        (n2,n3))
      end select

      if (if_urban_canopy > 0) then
         allocate (turb%cdrag   (n1,n2,n3))
      end if
      !------------------------------------------------------------------------------------!



      !----- Allocate the other arrays, which should be always allocated. -----------------!
      allocate (turb%kpbl          (n2,n3))
      allocate (turb%hkm        (n1,n2,n3))
      allocate (turb%vkm        (n1,n2,n3))
      allocate (turb%vkh        (n1,n2,n3))
      allocate (turb%sflux_u       (n2,n3))
      allocate (turb%sflux_v       (n2,n3))
      allocate (turb%sflux_w       (n2,n3))
      allocate (turb%sflux_t       (n2,n3))
      allocate (turb%sflux_r       (n2,n3))
      allocate (turb%sflux_c       (n2,n3))
      allocate (turb%akscal        (n2,n3))

      return
   end subroutine alloc_turb
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine nullifies all fields from the turbulence structure for a safe    !
   ! allocation.                                                                           !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_turb(turb)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (turb_vars), intent(inout) :: turb
      !------------------------------------------------------------------------------------!

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
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine de-allocates the turbulence structure.  This is used only        !
   ! dynamic allocation is on.                                                             !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_turb(turb)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (turb_vars), intent(inout) :: turb
      !------------------------------------------------------------------------------------!

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
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine flushes all arrays to zero (or to the minimum acceptable value), !
   ! so all variables are initialised when needed.                                         !
   !---------------------------------------------------------------------------------------!
   subroutine zero_turb(turb)
      use rconstants, only : tkmin ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (turb_vars), intent(inout) :: turb
      !------------------------------------------------------------------------------------!


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
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine fills the pointers to arrays into the variable table.            !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_turb(turb,turbm,imean,n1,n2,n3,ng)

      use var_tables

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (turb_vars), intent(in) :: turb
      type (turb_vars), intent(in) :: turbm
      integer         , intent(in) :: imean
      integer         , intent(in) :: n1
      integer         , intent(in) :: n2
      integer         , intent(in) :: n3
      integer         , intent(in) :: ng
      !----- Local variables. -------------------------------------------------------------!
      integer                      :: npts
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Attribute the pointers associated with the 3-D variables (n1,n2,n3)            !
      !------------------------------------------------------------------------------------!
      npts = n1 * n2 * n3

      if (associated(turb%tkep))                                                           &
         call vtables2(turb%tkep,turbm%tkep,ng,npts,imean                                  &
                      ,'TKEP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(turb%epsp))                                                           &
         call vtables2(turb%epsp,turbm%epsp,ng,npts,imean                                  &
           ,          ,'EPSP :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(turb%hkm))                                                            &
         call vtables2(turb%hkm,turbm%hkm,ng,npts,imean                                    &
                      ,'HKM :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(turb%vkm))                                                            &
         call vtables2(turb%vkm,turbm%vkm,ng,npts,imean                                    &
                      ,'VKM :3:hist:mpti:mpt3:mpt1')

      if (associated(turb%vkh))                                                            &
         call vtables2(turb%vkh,turbm%vkh,ng,npts,imean                                    &
                      ,'VKH :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(turb%cdrag))                                                          &
         call vtables2(turb%cdrag,turbm%cdrag,ng,npts,imean                                &
                      ,'CDRAG :3:hist:anal:mpti:mpt3')

      if (associated(turb%ltscale))                                                        &
         call vtables2(turb%ltscale,turbm%ltscale,ng,npts,imean                            &
                      ,'TL :3:hist:anal:mpti:mpt3')

      if (associated(turb%sigw))                                                           &
         call vtables2(turb%sigw,turbm%sigw,ng,npts,imean                                  &
                      ,'SIGW :3:hist:anal:mpti:mpt3')
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Attribute the pointers associated with the 2-D variables (n2,n3)               !
      !------------------------------------------------------------------------------------!
      npts = n2 * n3

      if (associated(turb%sflux_u))                                                        &
         call vtables2(turb%sflux_u,turbm%sflux_u,ng,npts,imean                            &
                      ,'SFLUX_U :2:anal:mpt3:mpt1')

      if (associated(turb%sflux_v))                                                        &
         call vtables2(turb%sflux_v,turbm%sflux_v,ng,npts,imean                            &
                      ,'SFLUX_V :2:anal:mpt3:mpt1')

      if (associated(turb%sflux_w))                                                        &
         call vtables2(turb%sflux_w,turbm%sflux_w,ng,npts,imean                            &
                      ,'SFLUX_W :2:anal:mpt3')

      if (associated(turb%sflux_t))                                                        &
         call vtables2(turb%sflux_t,turbm%sflux_t,ng,npts,imean                            &
                      ,'SFLUX_T :2:anal:mpt3')

      if (associated(turb%sflux_r))                                                        &
         call vtables2(turb%sflux_r,turbm%sflux_r,ng,npts,imean                            &
                      ,'SFLUX_R :2:anal:mpt3')

      if (associated(turb%sflux_c))                                                        &
         call vtables2(turb%sflux_c,turbm%sflux_c,ng,npts,imean                            &
                      ,'SFLUX_C :2:anal:mpt3')

      if (associated(turb%akscal))                                                         &
         call vtables2(turb%akscal,turbm%akscal,ng,npts,imean                              &
                      ,'AKSCAL :2:hist:anal:mpti:mpt3')

      if (associated(turb%pblhgt))                                                         &
         call vtables2(turb%pblhgt,turbm%pblhgt,ng,npts,imean                              &
                      ,'PBLHGT :2:hist:anal:mpti:mpt3')

      if (associated(turb%lmo))                                                            &
         call vtables2(turb%lmo,turbm%lmo,ng,npts,imean                                    &
                      ,'LMO    :2:hist:anal:mpti:mpt3')
      !------------------------------------------------------------------------------------!

      return
   end subroutine filltab_turb
   !=======================================================================================!
   !=======================================================================================!
end module mem_turb
!==========================================================================================!
!==========================================================================================!

!==========================================================================================!
!==========================================================================================!
!  Module mem_mass - This module contains variables associated with mass conservation      !
!                    check and mass flux related variables for Lagrangian Particle         !
!                    Dispersion Models (currently it is fitted for STILT).                 !
!------------------------------------------------------------------------------------------!
module mem_mass
   use grid_dims, only : maxgrds
   type mass_vars
      !------------------------------------------------------------------------------------!
      !   3-D variables (nzp,nxp,nyp)                                                      !
      !------------------------------------------------------------------------------------!
      !----- Full Exner function equation variables. --------------------------------------!
      real, pointer, dimension(:,:,:)   :: thvlast
      real, pointer, dimension(:,:,:)   :: lnthetav
      real, pointer, dimension(:,:,:)   :: lnthvadv
      real, pointer, dimension(:,:,:)   :: lnthvtend
      !----- Advective fluxes. ------------------------------------------------------------!
      real, pointer, dimension(:,:,:)   :: afxu
      real, pointer, dimension(:,:,:)   :: afxv
      real, pointer, dimension(:,:,:)   :: afxw
      !----- Averaged variables: Turbulence-related variables. ----------------------------!
      real, pointer, dimension(:,:,:)   :: ltscaleb
      real, pointer, dimension(:,:,:)   :: sigwb
      real, pointer, dimension(:,:,:)   :: tkepb
      !----- Averaged variables: Advective fluxes. ----------------------------------------!
      real, pointer, dimension(:,:,:)   :: afxub
      real, pointer, dimension(:,:,:)   :: afxvb
      real, pointer, dimension(:,:,:)   :: afxwb
      !------------------------------------------------------------------------------------!
      !   4-D variables (nzp,nxp,nyp,nclouds)                                              !
      !------------------------------------------------------------------------------------!
      !----- Convective fluxes. -----------------------------------------------------------!
      real, pointer, dimension(:,:,:,:) :: cfxup
      real, pointer, dimension(:,:,:,:) :: cfxdn
      real, pointer, dimension(:,:,:,:) :: dfxup
      real, pointer, dimension(:,:,:,:) :: efxup
      real, pointer, dimension(:,:,:,:) :: dfxdn
      real, pointer, dimension(:,:,:,:) :: efxdn
      !------------------------------------------------------------------------------------!
   end type

   type (mass_vars), allocatable :: mass_g(:), massm_g(:)
   integer                       :: iexev,imassflx
   real                          :: frqmassave

   !----- These variables control the time when the averages should be reset. -------------!
   real   , dimension(maxgrds)   :: etime_adve
   real   , dimension(maxgrds)   :: etime_turb

   contains
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine checks the most suitable time average. If the user is outputting      !
! lite analysis and the variables are in the lite analysis, then the mass flux averaging   !
! is done over the lite analysis. Otherwise it will be done at the full analysis           !
!------------------------------------------------------------------------------------------!
   subroutine define_frqmassave(ngrids,frqlite,frqanl,idiffk,maxlite,nlite_vars    &
                               ,lite_vars,mynum)
      implicit none
      integer, intent(in)                               :: ngrids     ! # of grids
      real   , intent(in)                               :: frqlite    ! Lite analysis freq.
      real   , intent(in)                               :: frqanl     ! Full analysis freq.
      integer, dimension(ngrids)           , intent(in) :: idiffk     ! Turbulence closure
      integer, intent(in)                               :: maxlite    ! Maximum # of lite
      integer, intent(in)                               :: nlite_vars ! Actual  # of lite
      character(len=32), dimension(maxlite), intent(in) :: lite_vars  ! Lite variables
      integer, intent(in)                               :: mynum      ! Node number

      integer                                           :: nlite      ! Maximum loop
      integer                                           :: l          ! Counter
      integer                                           :: ng         ! Counter
      integer                                           :: listed     ! Number of listed
      integer                                           :: benchmark  ! # of variables
      
      if (frqlite == 0.) then
         frqmassave=frqanl
      else
         !---------------------------------------------------------------------------------!
         !    frqlite will be used for average only if frqlite is set differently than     !
         ! zero and all needed variables are listed. This also depends on the configura-   !
         ! tion.                                                                           !
         !---------------------------------------------------------------------------------!
         if (any(idiffk == 7) .or. any(idiffk == 8)) then
            benchmark = 6
         elseif (any(idiffk == 1) .or. any(idiffk == 4) .or. &
                 any(idiffk == 5) .or. any(idiffk == 6)) then
            benchmark = 4
         else
            benchmark = 3
         end if
         nlite = min(maxlite,nlite_vars)
         listed = 0
         !---------------------------------------------------------------------------------!
         !   Loop through all variables.                                                   !
         !---------------------------------------------------------------------------------!
         do l=1,nlite
            select case (trim(lite_vars(l)))
            case ('AFXUB','AFXVB','AFXWB')
               listed=listed+1
            case ('TKEPB')
               if (any(idiffk == 1) .or. any(idiffk == 4) .or.  &
                   any(idiffk == 5) .or. any(idiffk == 6) .or.  &
                   any(idiffk == 7) .or. any(idiffk == 8)) then
                   listed=listed+1
               end if
            case ('TLB','SIGWB')
               if (any(idiffk == 7) .or. any(idiffk == 8)) then
                  listed = listed + 1
               end if
            end select
         end do
         
         if (listed == benchmark) then
            frqmassave = frqlite
         else
            frqmassave = frqanl
         end if
      end if
      
      !----- Initialise the variables which control the time averages. --------------------!
      etime_adve  (:) = 0.
      etime_turb  (:) = 0.
      !------------------------------------------------------------------------------------!

      return
   end subroutine define_frqmassave
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
   subroutine alloc_mass(mass,n1,n2,n3,nclouds,ng,nnqparmg,idiffkg)
      implicit none
      type (mass_vars) :: mass
      integer, intent(in) :: nnqparmg,idiffkg
      integer, intent(in) :: n1, n2, n3, nclouds, ng
 
   ! Allocate arrays based on options (if necessary)
      if (iexev == 2) then ! Full Exner function
         allocate (mass%thvlast(n1,n2,n3))
         allocate (mass%lnthetav(n1,n2,n3))
         allocate (mass%lnthvadv(n1,n2,n3))
         allocate (mass%lnthvtend(n1,n2,n3))
      end if
      if (imassflx == 1) then
         allocate (mass%afxu(n1,n2,n3))
         allocate (mass%afxv(n1,n2,n3))
         allocate (mass%afxw(n1,n2,n3))
         allocate (mass%afxub(n1,n2,n3))
         allocate (mass%afxvb(n1,n2,n3))
         allocate (mass%afxwb(n1,n2,n3))
         if (nnqparmg == 1) then
            allocate (mass%cfxup(n1,n2,n3,nclouds))
            allocate (mass%cfxdn(n1,n2,n3,nclouds))
            allocate (mass%dfxup(n1,n2,n3,nclouds))
            allocate (mass%efxup(n1,n2,n3,nclouds))
            allocate (mass%dfxdn(n1,n2,n3,nclouds))
            allocate (mass%efxdn(n1,n2,n3,nclouds))
         end if
         if (idiffkg /= 2 .and. idiffkg /= 3) then
            allocate (mass%tkepb(n1,n2,n3))
         end if
         if (idiffkg == 7 .or. idiffkg == 8) allocate (mass%ltscaleb(n1,n2,n3))
         if (idiffkg == 1 .or. idiffkg == 7 .or. idiffkg == 8) then
            allocate (mass%sigwb(n1,n2,n3))
         end if
      end if
      return
   end subroutine alloc_mass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
   subroutine nullify_mass(mass)
      implicit none
      type (mass_vars)   :: mass

      nullify (mass%thvlast    )
      nullify (mass%lnthetav   )
      nullify (mass%lnthvadv   ) 
      nullify (mass%lnthvtend  ) 
      nullify (mass%afxu       )
      nullify (mass%afxv       )
      nullify (mass%afxw       )
      nullify (mass%ltscaleb   )
      nullify (mass%sigwb      )
      nullify (mass%tkepb      )
      nullify (mass%afxub      )
      nullify (mass%afxvb      )
      nullify (mass%afxwb      )
      nullify (mass%cfxup      )
      nullify (mass%cfxdn      )
      nullify (mass%dfxup      )
      nullify (mass%efxup      )
      nullify (mass%dfxdn      )
      nullify (mass%efxdn      )
      return
   end subroutine nullify_mass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine flushes the average variables to zero whenever necessary. This is     !
! done at the output subroutine, right after the analysis (lite and full) are saved.       !
!------------------------------------------------------------------------------------------!
   subroutine zero_average_mass_adve(mass)
      implicit none
      type (mass_vars)   :: mass

      if (associated(mass%afxub   ))  mass%afxub   = 0.0
      if (associated(mass%afxvb   ))  mass%afxvb   = 0.0
      if (associated(mass%afxwb   ))  mass%afxwb   = 0.0

      return
   end subroutine zero_average_mass_adve
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine flushes the average variables to zero whenever necessary. This is     !
! done at the output subroutine, right after the analysis (lite and full) are saved.       !
!------------------------------------------------------------------------------------------!
   subroutine zero_average_mass_turb(mass)
      implicit none
      type (mass_vars)   :: mass

      if (associated(mass%ltscaleb))  mass%ltscaleb= 0.0
      if (associated(mass%sigwb   ))  mass%sigwb   = 0.0
      if (associated(mass%tkepb   ))  mass%tkepb   = 0.0

      return
   end subroutine zero_average_mass_turb
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
   subroutine zero_mass(mass)
      implicit none
      type (mass_vars)   :: mass

      if (associated(mass%thvlast    ))  mass%thvlast    = 0.0
      if (associated(mass%lnthetav   ))  mass%lnthetav   = 0.0
      if (associated(mass%lnthvadv   ))  mass%lnthvadv   = 0.0
      if (associated(mass%lnthvtend  ))  mass%lnthvtend  = 0.0
      if (associated(mass%afxu       ))  mass%afxu       = 0.0
      if (associated(mass%afxv       ))  mass%afxv       = 0.0
      if (associated(mass%afxw       ))  mass%afxw       = 0.0
      if (associated(mass%ltscaleb   ))  mass%ltscaleb   = 0.0
      if (associated(mass%sigwb      ))  mass%sigwb      = 0.0
      if (associated(mass%tkepb      ))  mass%tkepb      = 0.0
      if (associated(mass%afxub      ))  mass%afxub      = 0.0
      if (associated(mass%afxvb      ))  mass%afxvb      = 0.0
      if (associated(mass%afxwb      ))  mass%afxwb      = 0.0
      if (associated(mass%cfxup      ))  mass%cfxup      = 0.0
      if (associated(mass%cfxdn      ))  mass%cfxdn      = 0.0
      if (associated(mass%dfxup      ))  mass%dfxup      = 0.0
      if (associated(mass%efxup      ))  mass%efxup      = 0.0
      if (associated(mass%dfxdn      ))  mass%dfxdn      = 0.0
      if (associated(mass%efxdn      ))  mass%efxdn      = 0.0
      return
   end subroutine zero_mass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
   subroutine dealloc_mass(mass)
      implicit none
      type (mass_vars)   :: mass

      if (associated(mass%thvlast    ))  deallocate (mass%thvlast    ) 
      if (associated(mass%lnthetav   ))  deallocate (mass%lnthetav   ) 
      if (associated(mass%lnthvadv   ))  deallocate (mass%lnthvadv   )
      if (associated(mass%lnthvtend  ))  deallocate (mass%lnthvtend  )
      if (associated(mass%afxu       ))  deallocate (mass%afxu       )
      if (associated(mass%afxv       ))  deallocate (mass%afxv       )
      if (associated(mass%afxw       ))  deallocate (mass%afxw       )
      if (associated(mass%ltscaleb   ))  deallocate (mass%ltscaleb   )
      if (associated(mass%sigwb      ))  deallocate (mass%sigwb      )
      if (associated(mass%tkepb      ))  deallocate (mass%tkepb      )
      if (associated(mass%afxub      ))  deallocate (mass%afxub      )
      if (associated(mass%afxvb      ))  deallocate (mass%afxvb      )
      if (associated(mass%afxwb      ))  deallocate (mass%afxwb      )
      if (associated(mass%cfxup      ))  deallocate (mass%cfxup      )
      if (associated(mass%cfxdn      ))  deallocate (mass%cfxdn      )
      if (associated(mass%dfxup      ))  deallocate (mass%dfxup      )
      if (associated(mass%efxup      ))  deallocate (mass%efxup      )
      if (associated(mass%dfxdn      ))  deallocate (mass%dfxdn      )
      if (associated(mass%efxdn      ))  deallocate (mass%efxdn      )
      return
   end subroutine dealloc_mass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
   subroutine filltab_mass(mass,massm,imean,n1,n2,n3,nclds,ng)
 
      use var_tables
     
      implicit none
      type (mass_vars)   :: mass,massm
      integer             :: imean, n1, n2, n3, nclds, ng
      integer             :: npts
      real, pointer       :: var, varm

      npts=n1*n2*n3

      if (associated(mass%thvlast)) &
         call vtables2 (mass%thvlast,massm%thvlast &
                    ,ng, npts, imean, &
                    'THVLAST :3:hist:mpti:mpt3:mpt1')

      if (associated(mass%lnthvadv)) &
         call vtables2 (mass%lnthvadv,massm%lnthvadv &
                    ,ng, npts, imean, &
                    'LNTHVADV :3:mpti:mpt3:mpt1')

      if (associated(mass%lnthetav)) &
         call vtables2 (mass%lnthetav,massm%lnthetav &
                    ,ng, npts, imean, &
                    'LNTHETAV :3:mpti:mpt3:mpt1')
                    
      if (associated(mass%lnthvtend )) &
         call vtables2 (mass%lnthvtend,massm%lnthvtend &
                    ,ng, npts, imean, &
                    'LNTHVTEND :3:mpti:mpt3:mpt1')
     
      if (associated(mass%afxu)) &
         call vtables2 (mass%afxu,massm%afxu &
                    ,ng, npts, imean, &
                    'AFXU :3:hist:anal:mpti:mpt3')

      if (associated(mass%afxv)) &
         call vtables2 (mass%afxv,massm%afxv &
                    ,ng, npts, imean, &
                    'AFXV :3:hist:anal:mpti:mpt3')

      if (associated(mass%afxw)) &
         call vtables2 (mass%afxw,massm%afxw &
                    ,ng, npts, imean, &
                    'AFXW :3:hist:anal:mpti:mpt3')

      if (associated(mass%ltscaleb)) &
         call vtables2 (mass%ltscaleb,massm%ltscaleb &
                    ,ng, npts, imean, &
                    'TLB :3:hist:anal:mpti:mpt3')

      if (associated(mass%sigwb)) &
         call vtables2 (mass%sigwb,massm%sigwb &
                    ,ng, npts, imean, &
                    'SIGWB :3:hist:anal:mpti:mpt3')

      if (associated(mass%tkepb)) &
         call vtables2 (mass%tkepb,massm%tkepb &
                    ,ng, npts, imean, &
                    'TKEPB :3:hist:anal:mpti:mpt3')

      if (associated(mass%afxub)) &
         call vtables2 (mass%afxub,massm%afxub &
                    ,ng, npts, imean, &
                    'AFXUB :3:hist:anal:mpti:mpt3')

      if (associated(mass%afxvb)) &
         call vtables2 (mass%afxvb,massm%afxvb &
                    ,ng, npts, imean, &
                    'AFXVB :3:hist:anal:mpti:mpt3')

      if (associated(mass%afxwb)) &
         call vtables2 (mass%afxwb,massm%afxwb &
                    ,ng, npts, imean, &
                    'AFXWB :3:hist:anal:mpti:mpt3')

      !------------------------------------------------------------------------------------!
      ! 4-D variables                                                                      !
      !------------------------------------------------------------------------------------!
      npts=n1*n2*n3*nclds

      if (associated(mass%cfxup)) &
         call vtables2 (mass%cfxup,massm%cfxup &
                    ,ng, npts, imean, &
                    'CFXUP :8:hist:anal:mpti:mpt3')

      if (associated(mass%cfxdn)) &
         call vtables2 (mass%cfxdn,massm%cfxdn &
                    ,ng, npts, imean, &
                    'CFXDN :8:hist:anal:mpti:mpt3')

      if (associated(mass%dfxup)) &
         call vtables2 (mass%dfxup,massm%dfxup &
                    ,ng, npts, imean, &
                    'DFXUP :8:hist:anal:mpti:mpt3')

      if (associated(mass%efxup)) &
         call vtables2 (mass%efxup,massm%efxup &
                    ,ng, npts, imean, &
                    'EFXUP :8:hist:anal:mpti:mpt3')

      if (associated(mass%dfxdn)) &
         call vtables2 (mass%dfxdn,massm%dfxdn &
                    ,ng, npts, imean, &
                    'DFXDN :8:hist:anal:mpti:mpt3')

      if (associated(mass%efxdn)) &
         call vtables2 (mass%efxdn,massm%efxdn &
                    ,ng, npts, imean, &
                    'EFXDN :8:hist:anal:mpti:mpt3')

      return
   end subroutine filltab_mass
!==========================================================================================!
!==========================================================================================!
end module mem_mass

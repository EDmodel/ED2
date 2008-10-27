!==========================================================================================!
!==========================================================================================!
!  Module mem_mass - This module contains variables associated with mass conservation      !
!                    check and mass flux related variables for Lagrangian Particle         !
!                    Dispersion Models (currently it is fitted for STILT).                 !
!------------------------------------------------------------------------------------------!
module mem_mass
 
   type mass_vars
      real, pointer, dimension(:,:,:) :: & ! 3-D variables (nzp,nxp,nyp)
                              !Full Exner function equation variables
                               thvlast,thvadv,thetav,thvtend  &
                              !Advective fluxes
                              ,afxu, afxv, afxw               &
                              !Averaged variables:
                              !   Mixing layer related variables
                              ,ltscaleb, sigwb,tkepb          &
                              !   Advective fluxes
                              , afxub, afxvb, afxwb

     real, pointer, dimension(:,:,:,:) :: & 
                              !   Convective fluxes - (nzp,nxp,nyp,nclouds)
                              cfxup, cfxdn, dfxup, efxup &
                              ,dfxdn, efxdn

   end type

   type (mass_vars), allocatable :: mass_g(:), massm_g(:)
   integer                       :: iexev,imassflx
   real                          :: frqmassave 

   contains
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine checks the most suitable time average. If the user is outputting      !
! lite analysis and the variables are in the lite analysis, then the mass flux averaging   !
! is done over the lite analysis. Otherwise it will be done at the full analysis           !
!------------------------------------------------------------------------------------------!
   subroutine define_frqmassave(frqlite,frqanl,ngrids,idiffk,maxlite,nlite_vars,lite_vars)
      implicit none
      real   , intent(in)                               :: frqlite    ! Lite analysis freq.
      real   , intent(in)                               :: frqanl     ! Full analysis freq.
      integer, intent(in)                               :: ngrids     ! # of grids
      integer, dimension(ngrids)           , intent(in) :: idiffk     ! Turbulence closure
      integer, intent(in)                               :: maxlite    ! Maximum # of lite
      integer, intent(in)                               :: nlite_vars ! Actual  # of lite
      character(len=32), dimension(maxlite), intent(in) :: lite_vars  ! Lite variables

      integer                               :: nlite      ! Maximum loop
      integer                               :: l          ! Counter
      integer                               :: listed     ! Number of listed
      integer                               :: benchmark  ! Number of variables needed
      
      if (frqlite == 0.) then
         frqmassave=frqanl
         return
      else
         !---------------------------------------------------------------------------------!
         !    frqlite will be used for average only if frqlite is set differently than     !
         ! zero and all needed variables are listed. This also depends on the configura-   !
         ! tion.                                                                           !
         !---------------------------------------------------------------------------------!
         if (any(idiffk == 7)) then
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
                   any(idiffk == 7)) then
                   listed=listed+1
               end if
            case ('TLB','SIGWB')
               if (any(idiffk == 7)) then
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
         allocate (mass%thvadv(n1,n2,n3))
         allocate (mass%thetav(n1,n2,n3))
         allocate (mass%thvtend(n1,n2,n3))
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
         if (idiffkg == 7) allocate (mass%ltscaleb(n1,n2,n3))
         if (idiffkg == 1 .or. idiffkg == 7) allocate (mass%sigwb(n1,n2,n3))
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

      if (associated(mass%thvlast ))  nullify (mass%thvlast )
      if (associated(mass%thvadv  ))  nullify (mass%thvadv  )
      if (associated(mass%thetav  ))  nullify (mass%thetav  )
      if (associated(mass%thvtend ))  nullify (mass%thvtend )
      if (associated(mass%afxu    ))  nullify (mass%afxu    )
      if (associated(mass%afxv    ))  nullify (mass%afxv    )
      if (associated(mass%afxw    ))  nullify (mass%afxw    )
      if (associated(mass%ltscaleb))  nullify (mass%ltscaleb)
      if (associated(mass%sigwb   ))  nullify (mass%sigwb   )
      if (associated(mass%tkepb   ))  nullify (mass%tkepb   ) 
      if (associated(mass%afxub   ))  nullify (mass%afxub   )
      if (associated(mass%afxvb   ))  nullify (mass%afxvb   )
      if (associated(mass%afxwb   ))  nullify (mass%afxwb   )
      if (associated(mass%cfxup   ))  nullify (mass%cfxup   )
      if (associated(mass%cfxdn   ))  nullify (mass%cfxdn   )
      if (associated(mass%dfxup   ))  nullify (mass%dfxup   )
      if (associated(mass%efxup   ))  nullify (mass%efxup   )
      if (associated(mass%dfxdn   ))  nullify (mass%dfxdn   )
      if (associated(mass%efxdn   ))  nullify (mass%efxdn   )
      return
   end subroutine nullify_mass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine flushes the average variables to zero whenever necessary. This is     !
! done at the output subroutine, right after the analysis (lite and full) are saved.       !
!------------------------------------------------------------------------------------------!
   subroutine zero_average_mass(mass)
      implicit none
      type (mass_vars)   :: mass

      if (associated(mass%ltscaleb))  mass%ltscaleb= 0.0
      if (associated(mass%sigwb   ))  mass%sigwb   = 0.0
      if (associated(mass%tkepb   ))  mass%tkepb   = 0.0
      if (associated(mass%afxub   ))  mass%afxub   = 0.0
      if (associated(mass%afxvb   ))  mass%afxvb   = 0.0
      if (associated(mass%afxwb   ))  mass%afxwb   = 0.0

      return
   end subroutine zero_average_mass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
   subroutine zero_mass(mass)
      implicit none
      type (mass_vars)   :: mass

      if (associated(mass%thvlast ))  mass%thvlast = 0.0
      if (associated(mass%thvadv  ))  mass%thvadv  = 0.0
      if (associated(mass%thetav  ))  mass%thetav  = 0.0
      if (associated(mass%thvtend ))  mass%thvtend = 0.0
      if (associated(mass%afxu    ))  mass%afxu    = 0.0
      if (associated(mass%afxv    ))  mass%afxv    = 0.0
      if (associated(mass%afxw    ))  mass%afxw    = 0.0
      if (associated(mass%ltscaleb))  mass%ltscaleb= 0.0
      if (associated(mass%sigwb   ))  mass%sigwb   = 0.0
      if (associated(mass%tkepb   ))  mass%tkepb   = 0.0
      if (associated(mass%afxub   ))  mass%afxub   = 0.0
      if (associated(mass%afxvb   ))  mass%afxvb   = 0.0
      if (associated(mass%afxwb   ))  mass%afxwb   = 0.0
      if (associated(mass%cfxup   ))  mass%cfxup   = 0.0
      if (associated(mass%cfxdn   ))  mass%cfxdn   = 0.0
      if (associated(mass%dfxup   ))  mass%dfxup   = 0.0
      if (associated(mass%efxup   ))  mass%efxup   = 0.0
      if (associated(mass%dfxdn   ))  mass%dfxdn   = 0.0
      if (associated(mass%efxdn   ))  mass%efxdn   = 0.0
      return
   end subroutine zero_mass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
   subroutine dealloc_mass(mass)
      implicit none
      type (mass_vars)   :: mass

      if (associated(mass%thvlast ))  deallocate (mass%thvlast )
      if (associated(mass%thvadv  ))  deallocate (mass%thvadv  )
      if (associated(mass%thetav  ))  deallocate (mass%thetav  )
      if (associated(mass%thvtend ))  deallocate (mass%thvtend )
      if (associated(mass%afxu    ))  deallocate (mass%afxu    )
      if (associated(mass%afxv    ))  deallocate (mass%afxv    )
      if (associated(mass%afxw    ))  deallocate (mass%afxw    )
      if (associated(mass%ltscaleb))  deallocate (mass%ltscaleb)
      if (associated(mass%sigwb   ))  deallocate (mass%sigwb   )
      if (associated(mass%tkepb   ))  deallocate (mass%tkepb   )
      if (associated(mass%afxub   ))  deallocate (mass%afxub   )
      if (associated(mass%afxvb   ))  deallocate (mass%afxvb   )
      if (associated(mass%afxwb   ))  deallocate (mass%afxwb   )
      if (associated(mass%cfxup   ))  deallocate (mass%cfxup   )
      if (associated(mass%cfxdn   ))  deallocate (mass%cfxdn   )
      if (associated(mass%dfxup   ))  deallocate (mass%dfxup   )
      if (associated(mass%efxup   ))  deallocate (mass%efxup   )
      if (associated(mass%dfxdn   ))  deallocate (mass%dfxdn   )
      if (associated(mass%efxdn   ))  deallocate (mass%efxdn   )
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
         call vtables2 (mass%thvlast(1,1,1),massm%thvlast(1,1,1) &
                    ,ng, npts, imean, &
                    'THVLAST :3:hist:mpti:mpt3:mpt1')

      if (associated(mass%thvadv)) &
         call vtables2 (mass%thvadv(1,1,1),massm%thvadv(1,1,1) &
                    ,ng, npts, imean, &
                    'THVADV :3:mpti:mpt3:mpt1')

      if (associated(mass%thetav)) &
         call vtables2 (mass%thetav(1,1,1),massm%thetav(1,1,1) &
                    ,ng, npts, imean, &
                    'THETAV :3:mpti:mpt3:mpt1')
                    
      if (associated(mass%thvtend )) &
         call vtables2 (mass%thvtend(1,1,1),massm%thvtend(1,1,1) &
                    ,ng, npts, imean, &
                    'THVTEND :3:mpti:mpt3:mpt1')
     
      if (associated(mass%afxu)) &
         call vtables2 (mass%afxu(1,1,1),massm%afxu(1,1,1) &
                    ,ng, npts, imean, &
                    'AFXU :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%afxv)) &
         call vtables2 (mass%afxv(1,1,1),massm%afxv(1,1,1) &
                    ,ng, npts, imean, &
                    'AFXV :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%afxw)) &
         call vtables2 (mass%afxw(1,1,1),massm%afxw(1,1,1) &
                    ,ng, npts, imean, &
                    'AFXW :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%ltscaleb)) &
         call vtables2 (mass%ltscaleb(1,1,1),massm%ltscaleb(1,1,1) &
                    ,ng, npts, imean, &
                    'TLB :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%sigwb)) &
         call vtables2 (mass%sigwb(1,1,1),massm%sigwb(1,1,1) &
                    ,ng, npts, imean, &
                    'SIGWB :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%tkepb)) &
         call vtables2 (mass%tkepb(1,1,1),massm%tkepb(1,1,1) &
                    ,ng, npts, imean, &
                    'TKEPB :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%afxub)) &
         call vtables2 (mass%afxub(1,1,1),massm%afxub(1,1,1) &
                    ,ng, npts, imean, &
                    'AFXUB :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%afxvb)) &
         call vtables2 (mass%afxvb(1,1,1),massm%afxvb(1,1,1) &
                    ,ng, npts, imean, &
                    'AFXVB :3:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%afxwb)) &
         call vtables2 (mass%afxwb(1,1,1),massm%afxwb(1,1,1) &
                    ,ng, npts, imean, &
                    'AFXWB :3:hist:anal:mpti:mpt3:mpt1')

      !------------------------------------------------------------------------------------!
      ! 4-D variables                                                                      !
      !------------------------------------------------------------------------------------!
      npts=n1*n2*n3*nclds

      if (associated(mass%cfxup)) &
         call vtables2 (mass%cfxup(1,1,1,1),massm%cfxup(1,1,1,1) &
                    ,ng, npts, imean, &
                    'CFXUP :8:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%cfxdn)) &
         call vtables2 (mass%cfxdn(1,1,1,1),massm%cfxdn(1,1,1,1) &
                    ,ng, npts, imean, &
                    'CFXDN :8:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%dfxup)) &
         call vtables2 (mass%dfxup(1,1,1,1),massm%dfxup(1,1,1,1) &
                    ,ng, npts, imean, &
                    'DFXUP :8:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%efxup)) &
         call vtables2 (mass%efxup(1,1,1,1),massm%efxup(1,1,1,1) &
                    ,ng, npts, imean, &
                    'EFXUP :8:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%dfxdn)) &
         call vtables2 (mass%dfxdn(1,1,1,1),massm%dfxdn(1,1,1,1) &
                    ,ng, npts, imean, &
                    'DFXDN :8:hist:anal:mpti:mpt3:mpt1')

      if (associated(mass%efxdn)) &
         call vtables2 (mass%efxdn(1,1,1,1),massm%efxdn(1,1,1,1) &
                    ,ng, npts, imean, &
                    'EFXDN :8:hist:anal:mpti:mpt3:mpt1')

      return
   end subroutine filltab_mass
!==========================================================================================!
!==========================================================================================!
end module mem_mass

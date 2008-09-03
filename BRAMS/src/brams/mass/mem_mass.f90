!############################# Change Log ##################################
! 3.1.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003, 2004 - All Rights Reserved
!  Brazilian Regional Atmospheric Modeling System - BRAMS
!###########################################################################

module mem_mass
 
 use grid_dims
 use io_params, only: avgtim
 use mem_turb, only : idiffk
 
 type mass_vars
   real, pointer, dimension(:,:,:) :: &
                           !Full Exner function equation variables
                            thvlast,thvadv,thetav,thvtend  &
                           !Advective fluxes
                           ,afxu, afxv, afxw               &
                           !Averaged variables:
                           !   Mixing layer related variables
                           ,ltscaleb, sigwb,tkepb          &
                           !   Advective fluxes
                           , afxub, afxvb, afxwb           &
                           !   Convective fluxes
                           ,cfxup1, cfxdn1, dfxup1, efxup1 &
                           ,dfxdn1, efxdn1, cfxup2, dfxup2 &
                           ,efxup2 
 end type

 type (mass_vars), allocatable :: mass_g(:), massm_g(:)
 integer                        :: iexev,imassflx
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
         !----------------------------------------------------------------------------------!
         !    frqlite will be used for average only if frqlite is set differently than zero !
         ! and all needed variables are listed. This also depends on the configuration.     !
         !----------------------------------------------------------------------------------!
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

!------------------------------------------------------------------------------------------!
 subroutine alloc_mass(mass,n1,n2,n3,ng)
   implicit none
   type (mass_vars) :: mass
   integer, intent(in) :: n1, n2, n3, ng
 
! Allocate arrays based on options (if necessary)
   if (iexev == 2) then
     allocate (mass%thvlast(n1,n2,n3))
     allocate (mass%thvadv(n1,n2,n3))
     allocate (mass%thetav(n1,n2,n3))
     allocate (mass%thvtend(n1,n2,n3))
   end if
   if (imassflx == 1) then
     allocate (mass%afxu(n1,n2,n3))
     allocate (mass%afxv(n1,n2,n3))
     allocate (mass%afxw(n1,n2,n3))
     allocate (mass%cfxup1(n1,n2,n3))
     allocate (mass%cfxdn1(n1,n2,n3))
     allocate (mass%dfxup1(n1,n2,n3))
     allocate (mass%efxup1(n1,n2,n3))
     allocate (mass%dfxdn1(n1,n2,n3))
     allocate (mass%efxdn1(n1,n2,n3))
     allocate (mass%cfxup2(n1,n2,n3))
     allocate (mass%dfxup2(n1,n2,n3))
     allocate (mass%efxup2(n1,n2,n3))
     allocate (mass%afxub(n1,n2,n3))
     allocate (mass%afxvb(n1,n2,n3))
     allocate (mass%afxwb(n1,n2,n3))
     if (idiffk(ng) /= 2 .and. idiffk(ng) /= 3) then
       allocate (mass%tkepb(n1,n2,n3))
     end if
     if (idiffk(ng) == 7) then
       allocate (mass%ltscaleb(n1,n2,n3))
       allocate (mass%sigwb(n1,n2,n3))
     end if
   end if
   return
 end subroutine alloc_mass

!------------------------------------------------------------------------------------------!

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
   if (associated(mass%cfxup1  ))  nullify (mass%cfxup1  )
   if (associated(mass%cfxdn1  ))  nullify (mass%cfxdn1  )
   if (associated(mass%dfxup1  ))  nullify (mass%dfxup1  )
   if (associated(mass%efxup1  ))  nullify (mass%efxup1  )
   if (associated(mass%dfxdn1  ))  nullify (mass%dfxdn1  )
   if (associated(mass%efxdn1  ))  nullify (mass%efxdn1  )
   if (associated(mass%cfxup2  ))  nullify (mass%cfxup2  )
   if (associated(mass%dfxup2  ))  nullify (mass%dfxup2  )
   if (associated(mass%efxup2  ))  nullify (mass%efxup2  )
   return
 end subroutine nullify_mass

!------------------------------------------------------------------------------------------!

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
   if (associated(mass%cfxup1  ))  mass%cfxup1  = 0.0
   if (associated(mass%cfxdn1  ))  mass%cfxdn1  = 0.0
   if (associated(mass%dfxup1  ))  mass%dfxup1  = 0.0
   if (associated(mass%efxup1  ))  mass%efxup1  = 0.0
   if (associated(mass%dfxdn1  ))  mass%dfxdn1  = 0.0
   if (associated(mass%efxdn1  ))  mass%efxdn1  = 0.0
   if (associated(mass%cfxup2  ))  mass%cfxup2  = 0.0
   if (associated(mass%dfxup2  ))  mass%dfxup2  = 0.0
   if (associated(mass%efxup2  ))  mass%efxup2  = 0.0
   return
 end subroutine zero_mass

!------------------------------------------------------------------------------------------!

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
   if (associated(mass%cfxup1  ))  deallocate (mass%cfxup1  )
   if (associated(mass%cfxdn1  ))  deallocate (mass%cfxdn1  )
   if (associated(mass%dfxup1  ))  deallocate (mass%dfxup1  )
   if (associated(mass%efxup1  ))  deallocate (mass%efxup1  )
   if (associated(mass%dfxdn1  ))  deallocate (mass%dfxdn1  )
   if (associated(mass%efxdn1  ))  deallocate (mass%efxdn1  )
   if (associated(mass%cfxup2  ))  deallocate (mass%cfxup2  )
   if (associated(mass%dfxup2  ))  deallocate (mass%dfxup2  )
   if (associated(mass%efxup2  ))  deallocate (mass%efxup2  )
   return
 end subroutine dealloc_mass






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

!------------------------------------------------------------------------------------------!

 subroutine filltab_mass(mass,massm,imean,n1,n2,n3,ng)
 
   use var_tables
   
   implicit none
   type (mass_vars)   :: mass,massm
   integer             :: imean, n1, n2, n3, ng
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

   if (associated(mass%cfxup1)) &
      call vtables2 (mass%cfxup1(1,1,1),massm%cfxup1(1,1,1) &
                 ,ng, npts, imean, &
                 'CFXUP1 :3:hist:anal:mpti:mpt3:mpt1')

   if (associated(mass%cfxdn1)) &
      call vtables2 (mass%cfxdn1(1,1,1),massm%cfxdn1(1,1,1) &
                 ,ng, npts, imean, &
                 'CFXDN1 :3:hist:anal:mpti:mpt3:mpt1')

   if (associated(mass%dfxup1)) &
      call vtables2 (mass%dfxup1(1,1,1),massm%dfxup1(1,1,1) &
                 ,ng, npts, imean, &
                 'DFXUP1 :3:hist:anal:mpti:mpt3:mpt1')

   if (associated(mass%efxup1)) &
      call vtables2 (mass%efxup1(1,1,1),massm%efxup1(1,1,1) &
                 ,ng, npts, imean, &
                 'EFXUP1 :3:hist:anal:mpti:mpt3:mpt1')

   if (associated(mass%dfxdn1)) &
      call vtables2 (mass%dfxdn1(1,1,1),massm%dfxdn1(1,1,1) &
                 ,ng, npts, imean, &
                 'DFXDN1 :3:hist:anal:mpti:mpt3:mpt1')

   if (associated(mass%efxdn1)) &
      call vtables2 (mass%efxdn1(1,1,1),massm%efxdn1(1,1,1) &
                 ,ng, npts, imean, &
                 'EFXDN1 :3:hist:anal:mpti:mpt3:mpt1')

   if (associated(mass%cfxup2)) &
      call vtables2 (mass%cfxup2(1,1,1),massm%cfxup2(1,1,1) &
                 ,ng, npts, imean, &
                 'CFXUP2 :3:hist:anal:mpti:mpt3:mpt1')

   if (associated(mass%dfxup2)) &
      call vtables2 (mass%dfxup2(1,1,1),massm%dfxup2(1,1,1) &
                 ,ng, npts, imean, &
                 'DFXUP2 :3:hist:anal:mpti:mpt3:mpt1')

   if (associated(mass%efxup2)) &
      call vtables2 (mass%efxup2(1,1,1),massm%efxup2(1,1,1) &
                 ,ng, npts, imean, &
                 'EFXUP2 :3:hist:anal:mpti:mpt3:mpt1')
   return
 end subroutine filltab_mass

end module mem_mass

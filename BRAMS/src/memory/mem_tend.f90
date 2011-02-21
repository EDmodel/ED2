!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
module mem_tend


   !----- Tendency array structure. -------------------------------------------------------!
   type tend_vars
      real, dimension(:,:,:), pointer :: ut
      real, dimension(:,:,:), pointer :: vt
      real, dimension(:,:,:), pointer :: wt
      real, dimension(:,:,:), pointer :: pt
      real, dimension(:,:,:), pointer :: tht
      real, dimension(:,:,:), pointer :: rtt
      real, dimension(:,:,:), pointer :: co2t
      real, dimension(:,:,:), pointer :: rct
      real, dimension(:,:,:), pointer :: tket
      real, dimension(:,:,:), pointer :: epst
   end type

   !----- The tendency variable. ----------------------------------------------------------!
   type (tend_vars), dimension(:), allocatable :: tend_g

   !=======================================================================================!
   !=======================================================================================!



   contains


   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_tend(nmz,nmx,nmy,ngr,ntpts)

      use mem_basic     , only : basic_g   ! ! data type intent(in)
      use mem_turb      , only : turb_g    ! ! data type intent(in)
      use teb_spm_start , only : teb_spm   ! ! intent(in)
      use mem_gaspart   , only : gaspart_g ! ! data type intent(inout)
      use mem_emiss     , only : ichemi    & ! intent(in)
                               , isource   ! ! intent(in)
      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      integer, intent(in) :: nmz
      integer, intent(in) :: nmx
      integer, intent(in) :: nmy
      integer, intent(in) :: ngr
      integer, intent(in) :: ntpts
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Allocate arrays based on options (if necessary).  We are only checking grid 1, !
      ! since all grids must have same scalars defined.                                    !
      !------------------------------------------------------------------------------------!
      if (associated(basic_g(ngr)%up   ))   allocate (tend_g(ngr)%ut   (nmz,nmx,nmy))
      if (associated(basic_g(ngr)%vp   ))   allocate (tend_g(ngr)%vt   (nmz,nmx,nmy))
      if (associated(basic_g(ngr)%wp   ))   allocate (tend_g(ngr)%wt   (nmz,nmx,nmy))
      if (associated(basic_g(ngr)%pp   ))   allocate (tend_g(ngr)%pt   (nmz,nmx,nmy))
      if (associated(basic_g(ngr)%co2p ))   allocate (tend_g(ngr)%co2t (nmz,nmx,nmy))
      if (associated(basic_g(ngr)%thp  ))   allocate (tend_g(ngr)%tht  (nmz,nmx,nmy))
      if (associated(basic_g(ngr)%rtp  ))   allocate (tend_g(ngr)%rtt  (nmz,nmx,nmy))
      if (associated(turb_g(ngr)%tkep  ))   allocate (tend_g(ngr)%tket (nmz,nmx,nmy))
      if (associated(turb_g(ngr)%epsp  ))   allocate (tend_g(ngr)%epst (nmz,nmx,nmy))
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Define the TEB-related tendencies.  This could be updated to the grid-specific !
      ! domain.                                                                            !
      !------------------------------------------------------------------------------------!
      if (TEB_SPM==1 .and. isource == 1 .and. ngr == 1) then
         if (associated(gaspart_g(1)%pno) .and. (.not.associated(gaspart_g(1)%pnot)))      &
            allocate (gaspart_g(1)%pnot(ntpts))

         if (associated(gaspart_g(1)%pno2) .and. (.not.associated(gaspart_g(1)%pno2t)))    &
            allocate (gaspart_g(1)%pno2t(ntpts))

         if (associated(gaspart_g(1)%ppm25) .and. (.not.associated(gaspart_g(1)%ppm25t)))  &
            allocate (gaspart_g(1)%ppm25t(ntpts))

         if (associated(gaspart_g(1)%pco) .and. (.not.associated(gaspart_g(1)%pcot)))      &
            allocate (gaspart_g(1)%pcot(ntpts))

         if (associated(gaspart_g(1)%pso2) .and. (.not.associated(gaspart_g(1)%pso2t)))    &
            allocate (gaspart_g(1)%pso2t(ntpts))

         if (associated(gaspart_g(1)%pso4) .and. (.not.associated(gaspart_g(1)%pso4t)))    &
            allocate (gaspart_g(1)%pso4t(ntpts))

         if (associated(gaspart_g(1)%paer) .and. (.not.associated(gaspart_g(1)%paert)))    &
            allocate (gaspart_g(1)%paert(ntpts))

         if (associated(gaspart_g(1)%pvoc) .and. (.not.associated(gaspart_g(1)%pvoct)))    &
            allocate (gaspart_g(1)%pvoct(ntpts))

         if (ichemi==1) then
            if (associated(gaspart_g(1)%po3) .and. (.not.associated(gaspart_g(1)%po3t)))   &
               allocate (gaspart_g(1)%po3t(ntpts))

            if ( associated(gaspart_g(1)%prhco) .and.                                      &
                 (.not.associated(gaspart_g(1)%prhcot)))                                   &
               allocate (gaspart_g(1)%prhcot(ntpts))

            if (associated(gaspart_g(1)%pho2) .and. (.not.associated(gaspart_g(1)%pho2t))) &
               allocate (gaspart_g(1)%pho2t(ntpts))

            if (associated(gaspart_g(1)%po3p) .and. (.not.associated(gaspart_g(1)%po3pt))) &
               allocate (gaspart_g(1)%po3pt(ntpts))
            
            if (associated(gaspart_g(1)%po1d) .and. (.not.associated(gaspart_g(1)%po1dt))) &
               allocate (gaspart_g(1)%po1dt(ntpts))

            if (associated(gaspart_g(1)%pho) .and. (.not.associated(gaspart_g(1)%phot)))   &
               allocate (gaspart_g(1)%phot(ntpts))

            if (associated(gaspart_g(1)%proo) .and. (.not.associated(gaspart_g(1)%proot))) &
               allocate (gaspart_g(1)%proot(ntpts))

         end if

      elseif (teb_spm == 1 .and. isource == 1) then

         gaspart_g(ngr)%pnot   => gaspart_g(1)%pnot
         gaspart_g(ngr)%pno2t  => gaspart_g(1)%pno2t
         gaspart_g(ngr)%ppm25t => gaspart_g(1)%ppm25t
         gaspart_g(ngr)%pcot   => gaspart_g(1)%pcot
         gaspart_g(ngr)%pso2t  => gaspart_g(1)%pso2t
         gaspart_g(ngr)%pso4t  => gaspart_g(1)%pso4t
         gaspart_g(ngr)%paert  => gaspart_g(1)%paert
         gaspart_g(ngr)%pvoct  => gaspart_g(1)%pvoct

         if (ichemi==1) then
            gaspart_g(ngr)%po3t   => gaspart_g(1)%po3t
            gaspart_g(ngr)%prhcot => gaspart_g(1)%prhcot
            gaspart_g(ngr)%pho2t  => gaspart_g(1)%pho2t
            gaspart_g(ngr)%po3pt  => gaspart_g(1)%po3pt
            gaspart_g(ngr)%po1dt  => gaspart_g(1)%po1dt
            gaspart_g(ngr)%phot   => gaspart_g(1)%phot
            gaspart_g(ngr)%proot  => gaspart_g(1)%proot
         end if

      end if

      return
   end subroutine alloc_tend
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will reset the pointers of all elements of tend structure.        !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_tend(tend)
      use mem_scalar    , only : scalar_g  ! ! data type intent(inout)
      use teb_spm_start , only : teb_spm   ! ! intent(in)
      use mem_gaspart   , only : gaspart_g ! ! data type intent(inout)
      use mem_emiss     , only : ichemi    & ! intent(in)
                               , isource   ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(tend_vars), intent(inout) :: tend
      !------------------------------------------------------------------------------------!

      if (associated(tend%ut  ))  nullify (tend%ut  )
      if (associated(tend%vt  ))  nullify (tend%vt  )
      if (associated(tend%wt  ))  nullify (tend%wt  )
      if (associated(tend%pt  ))  nullify (tend%pt  )
      if (associated(tend%tht ))  nullify (tend%tht )
      if (associated(tend%rtt ))  nullify (tend%rtt )
      if (associated(tend%co2t))  nullify (tend%co2t)
      if (associated(tend%tket))  nullify (tend%tket)
      if (associated(tend%epst))  nullify (tend%epst)

      !----- TEB_SPM. ---------------------------------------------------------------------!
      if (teb_spm==1 .and. isource == 1) then
         if (associated(gaspart_g(1)%pnot  )) nullify (gaspart_g(1)%pnot  )
         if (associated(gaspart_g(1)%pno2t )) nullify (gaspart_g(1)%pno2t )
         if (associated(gaspart_g(1)%ppm25t)) nullify (gaspart_g(1)%ppm25t)
         if (associated(gaspart_g(1)%pcot  )) nullify (gaspart_g(1)%pcot  )
         if (associated(gaspart_g(1)%pso2t )) nullify (gaspart_g(1)%pso2t )
         if (associated(gaspart_g(1)%pso4t )) nullify (gaspart_g(1)%pso4t )
         if (associated(gaspart_g(1)%paert )) nullify (gaspart_g(1)%paert )
         if (associated(gaspart_g(1)%pvoct )) nullify (gaspart_g(1)%pvoct )

         if (ichemi==1) then
            if (associated(gaspart_g(1)%po3t  )) nullify (gaspart_g(1)%po3t  )
            if (associated(gaspart_g(1)%prhcot)) nullify (gaspart_g(1)%prhcot)
            if (associated(gaspart_g(1)%pho2t )) nullify (gaspart_g(1)%pho2t )
            if (associated(gaspart_g(1)%po3pt )) nullify (gaspart_g(1)%po3pt )
            if (associated(gaspart_g(1)%po1dt )) nullify (gaspart_g(1)%po1dt )
            if (associated(gaspart_g(1)%phot  )) nullify (gaspart_g(1)%phot  )
            if (associated(gaspart_g(1)%proot )) nullify (gaspart_g(1)%proot )
         end if
      end if

      return
   end subroutine nullify_tend
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will deallocate all elements of tend structure.                   !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_tend(tend)
      use mem_scalar    , only : scalar_g  ! ! data type intent(inout)
      use teb_spm_start , only : teb_spm   ! ! intent(in)
      use mem_gaspart   , only : gaspart_g ! ! data type intent(inout)
      use mem_emiss     , only : ichemi    & ! intent(in)
                               , isource   ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(tend_vars), intent(inout) :: tend
      !------------------------------------------------------------------------------------!

      if (associated(tend%ut  ))  nullify    (tend%ut  )
      if (associated(tend%vt  ))  deallocate (tend%vt  )
      if (associated(tend%wt  ))  deallocate (tend%wt  )
      if (associated(tend%pt  ))  deallocate (tend%pt  )
      if (associated(tend%tht ))  deallocate (tend%tht )
      if (associated(tend%rtt ))  deallocate (tend%rtt )
      if (associated(tend%co2t))  deallocate (tend%co2t)
      if (associated(tend%tket))  deallocate (tend%tket)
      if (associated(tend%epst))  deallocate (tend%epst)

      !----- TEB_SPM. ---------------------------------------------------------------------!
      if (teb_spm==1 .and. isource==1) then
         if (associated(gaspart_g(1)%pnot  )) deallocate (gaspart_g(1)%pnot  )
         if (associated(gaspart_g(1)%pno2t )) deallocate (gaspart_g(1)%pno2t )
         if (associated(gaspart_g(1)%ppm25t)) deallocate (gaspart_g(1)%ppm25t)
         if (associated(gaspart_g(1)%pcot  )) deallocate (gaspart_g(1)%pcot  )
         if (associated(gaspart_g(1)%pso2t )) deallocate (gaspart_g(1)%pso2t )
         if (associated(gaspart_g(1)%pso4t )) deallocate (gaspart_g(1)%pso4t )
         if (associated(gaspart_g(1)%paert )) deallocate (gaspart_g(1)%paert )
         if (associated(gaspart_g(1)%pvoct )) deallocate (gaspart_g(1)%pvoct )

         if (ichemi==1) then
            if (associated(gaspart_g(1)%po3t  )) deallocate (gaspart_g(1)%po3t  )
            if (associated(gaspart_g(1)%prhcot)) deallocate (gaspart_g(1)%prhcot)
            if (associated(gaspart_g(1)%pho2t )) deallocate (gaspart_g(1)%pho2t )
            if (associated(gaspart_g(1)%po3pt )) deallocate (gaspart_g(1)%po3pt )
            if (associated(gaspart_g(1)%po1dt )) deallocate (gaspart_g(1)%po1dt )
            if (associated(gaspart_g(1)%phot  )) deallocate (gaspart_g(1)%phot  )
            if (associated(gaspart_g(1)%proot )) deallocate (gaspart_g(1)%proot )
         end if
      end if

      return
   end subroutine dealloc_tend
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will fill pointers to arrays into scalar tables.                  !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_tend(nmz,nmx,nmy,naddsc,tend,basic,micro,turb,scalar,gaspart,ng)
      use mem_basic     , only : basic_vars   ! ! type
      use mem_micro     , only : micro_vars   ! ! type
      use mem_turb      , only : turb_vars    ! ! type
      use mem_scalar    , only : scalar_vars  ! ! type
      use teb_spm_start , only : TEB_SPM      ! ! intent(in)
      use mem_gaspart   , only : gaspart_vars ! ! type
      use mem_emiss     , only : ichemi       & ! intent(in)
                               , isource      ! ! intent(in)
      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      integer                               , intent(in)  :: nmz
      integer                               , intent(in)  :: nmx
      integer                               , intent(in)  :: nmy
      integer                               , intent(in)  :: naddsc
      type (tend_vars)                      , intent(in)  :: tend
      type (basic_vars)                     , intent(in)  :: basic
      type (micro_vars)                     , intent(in)  :: micro
      type (turb_vars)                      , intent(in)  :: turb
      type (scalar_vars) , dimension(naddsc), intent(in)  :: scalar
      type (gaspart_vars)                   , pointer     :: gaspart
      integer                               , intent(in)  :: ng
      !----- Local Variables. -------------------------------------------------------------!
      integer                                             :: nsc
      integer                                             :: npts
      character (len=7)                                   :: sname
      !------------------------------------------------------------------------------------!

      !------ Determine the number of points. ---------------------------------------------!
      npts = nmz * nmx * nmy

      if (associated(tend%tht)) then
         call vtables_scalar (npts,basic%thp,tend%tht,ng,'THP')
         call vtables_scalar_new (basic%thp,tend%tht,ng,'THP',npts)
      end if

      if (associated(tend%rtt)) then
         call vtables_scalar (npts,basic%rtp,tend%rtt,ng,'RTP')
         call vtables_scalar_new (basic%rtp,tend%rtt,ng,'RTP',npts)
      end if

      if (associated(tend%co2t)) then
         call vtables_scalar (npts,basic%co2p,tend%co2t,ng,'CO2P')
         call vtables_scalar_new (basic%co2p,tend%co2t,ng,'CO2P',npts)
      end if

      if (associated(tend%tket)) then
         call vtables_scalar (npts,turb%tkep,tend%tket,ng,'TKEP')
         call vtables_scalar_new (turb%tkep,tend%tket,ng,'TKEP',npts)
      end if

      if (associated(tend%epst)) then
         call vtables_scalar (npts,turb%epsp,tend%epst,ng,'EPSP')
         call vtables_scalar_new (turb%epsp,tend%epst,ng,'EPSP',npts)
      end if

      !----- TEB_SPM. ---------------------------------------------------------------------!
      if (TEB_SPM==1) then
         if (isource==1) then

            if (associated(gaspart%pnot)) then
               call vtables_scalar (npts,gaspart%pno,gaspart%pnot,ng,'PNO')
               call vtables_scalar_new (gaspart%pno,gaspart%pnot,ng,'PNO',npts)
            end if

            if (associated(gaspart%pno2t)) then
               call vtables_scalar (npts,gaspart%pno2,gaspart%pno2t,ng,'PNO2')
               call vtables_scalar_new (gaspart%pno2,gaspart%pno2t,ng,'PNO2',npts)
            end if

            if (associated(gaspart%ppm25t)) then
               call vtables_scalar (npts,gaspart%ppm25,gaspart%ppm25t,ng,'PM25')
               call vtables_scalar_new (gaspart%ppm25,gaspart%ppm25t,ng,'PM25',npts)
            end if

            if (associated(gaspart%pcot)) then
               call vtables_scalar (npts,gaspart%pco,gaspart%pcot,ng,'PCO')
               call vtables_scalar_new (gaspart%pco,gaspart%pcot,ng,'PCO',npts)
            end if
            
            if (associated(gaspart%pso2t)) then
               call vtables_scalar (npts,gaspart%pso2,gaspart%pso2t,ng,'PSO2')
               call vtables_scalar_new (gaspart%pso2,gaspart%pso2t,ng,'PSO2',npts)
            end if

            if (associated(gaspart%pso4t)) then
               call vtables_scalar (npts,gaspart%pso4,gaspart%pso4t,ng,'PSO4')
               call vtables_scalar_new (gaspart%pso4,gaspart%pso4t,ng,'PSO4',npts)
            end if

            if (associated(gaspart%paert)) then
               call vtables_scalar (npts,gaspart%paer,gaspart%paert,ng,'PAER')
               call vtables_scalar_new (gaspart%paer,gaspart%paert,ng,'PAER',npts)
            end if

            if (associated(gaspart%pvoct)) then
               call vtables_scalar (npts,gaspart%pvoc,gaspart%pvoct,ng,'PVOC')
               call vtables_scalar_new (gaspart%pvoc,gaspart%pvoct,ng,'PVOC',npts)
            end if

            if(ichemi==1) then
               if (associated(gaspart%po3t)) then
                  call vtables_scalar(npts,gaspart%po3,gaspart%po3t,ng,'PO3')
                  call vtables_scalar_new(gaspart%po3,gaspart%po3t,ng,'PO3',npts)
               end if

               if (associated(gaspart%prhcot)) then
                  call vtables_scalar(npts,gaspart%prhco,gaspart%prhcot,ng,'PRHCO')
                  call vtables_scalar_new(gaspart%prhco,gaspart%prhcot,ng,'PRHCO',npts)
               end if

               if (associated(gaspart%pho2t)) then
                  call vtables_scalar (npts,gaspart%pho2,gaspart%pho2t,ng,'PHO2')
                  call vtables_scalar_new (gaspart%pho2,gaspart%pho2t,ng,'PHO2',npts)
               end if

               if (associated(gaspart%po3pt)) then
                  call vtables_scalar (npts,gaspart%po3p,gaspart%po3pt,ng,'PO3P')
                  call vtables_scalar_new (gaspart%po3p,gaspart%po3pt,ng,'PO3P',npts)
               end if

               if (associated(gaspart%po1dt)) then
                  call vtables_scalar (npts,gaspart%po1d,gaspart%po1dt,ng,'PO1D')
                  call vtables_scalar_new (gaspart%po1d,gaspart%po1dt,ng,'PO1D',npts)
               end if

               if (associated(gaspart%phot)) then
                  call vtables_scalar (npts,gaspart%pho,gaspart%phot,ng,'PHO')
                  call vtables_scalar_new (gaspart%pho,gaspart%phot,ng,'PHO',npts)
               end if

               if (associated(gaspart%proot)) then
                  call vtables_scalar (npts,gaspart%proo,gaspart%proot,ng,'PROO')
                  call vtables_scalar_new (gaspart%proo,gaspart%proot,ng,'PROO',npts)
               end if
            end if
         end if
      end if

      do nsc=1,naddsc
         write(sname,'(a4,i3.3)') 'SCLP',nsc
         if (associated(scalar(nsc)%sclt)) then
            call vtables_scalar (npts,scalar(nsc)%sclp,scalar(nsc)%sclt,ng,sname)
            call vtables_scalar_new (scalar(nsc)%sclp,scalar(nsc)%sclt,ng,sname,npts)
         end if
      end do
   end subroutine filltab_tend
   !=======================================================================================!
   !=======================================================================================!
end module mem_tend
!==========================================================================================!
!==========================================================================================!

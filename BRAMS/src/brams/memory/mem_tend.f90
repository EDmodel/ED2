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
      real, dimension(:), pointer :: ut
      real, dimension(:), pointer :: vt
      real, dimension(:), pointer :: wt
      real, dimension(:), pointer :: pt
      real, dimension(:), pointer :: tht
      real, dimension(:), pointer :: rtt
      real, dimension(:), pointer :: co2t
      real, dimension(:), pointer :: rct
      real, dimension(:), pointer :: tket
      real, dimension(:), pointer :: epst
   end type

   !----- The tendency variable. ----------------------------------------------------------!
   type (tend_vars) :: tend

   !=======================================================================================!
   !=======================================================================================!



   contains


   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_tend(nmzp,nmxp,nmyp,ngrs,naddsc,proc_type)

      use mem_basic     , only : basic_g   ! ! data type intent(in)
      use mem_scalar    , only : scalar_g  ! ! data type intent(inout)
      use mem_turb      , only : turb_g    ! ! data type intent(in)
      use teb_spm_start , only : teb_spm   ! ! intent(in)
      use mem_gaspart   , only : gaspart_g ! ! data type intent(inout)
      use mem_emiss     , only : ichemi    & ! intent(in)
                               , isource   ! ! intent(in)
      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      integer                 , intent(in) :: ngrs
      integer                 , intent(in) :: proc_type
      integer                 , intent(in) :: naddsc
      integer, dimension(ngrs), intent(in) :: nmzp
      integer, dimension(ngrs), intent(in) :: nmxp
      integer, dimension(ngrs), intent(in) :: nmyp

      !----- Local Variables. -------------------------------------------------------------!
      integer                              :: ng
      integer                              :: ntpts
      integer                              :: nsc
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Find the maximum number of grid points needed for any grid.  Here we check    !
      ! which kind of node is this (master, slave, only node) to avoid allocating          !
      ! unecessary arrays.                                                                 !
      !------------------------------------------------------------------------------------!
      select case (proc_type)
      case (1) !----- Master node in a parallel run, no tendency is computed here... ------!
         ntpts = 1

      case default
         !----- Find the maximum number of points amongst all grids. ----------------------!
         ntpts = 0
         do ng=1,ngrs
            ntpts = max(nmxp(ng)*nmyp(ng)*nmzp(ng),ntpts)
         end do
      end select

      !------------------------------------------------------------------------------------!
      !     Allocate arrays based on options (if necessary).  We are only checking grid 1, !
      ! since all grids must have same scalars defined.                                    !
      !------------------------------------------------------------------------------------!
      if (associated(basic_g(1)%up   ))   allocate (tend%ut   (ntpts))
      if (associated(basic_g(1)%vp   ))   allocate (tend%vt   (ntpts))
      if (associated(basic_g(1)%wp   ))   allocate (tend%wt   (ntpts))
      if (associated(basic_g(1)%pp   ))   allocate (tend%pt   (ntpts))
      if (associated(basic_g(1)%co2p ))   allocate (tend%rtt  (ntpts))
      if (associated(basic_g(1)%thp  ))   allocate (tend%tht  (ntpts))
      if (associated(basic_g(1)%rtp  ))   allocate (tend%rtt  (ntpts))
      if (associated(turb_g(1)%tkep  ))   allocate (tend%tket (ntpts))
      if (associated(turb_g(1)%epsp  ))   allocate (tend%epst (ntpts))

      if (TEB_SPM==1) then
         if (isource==1) then
            if (associated(gaspart_g(1)%pno).and.                                          &
                 (.not.associated(gaspart_g(1)%pnot)))                                     &
                 allocate (gaspart_g(1)%pnot(ntpts))

            if (associated(gaspart_g(1)%pno2).and.                                         &
                 (.not.associated(gaspart_g(1)%pno2t)))                                    &
                 allocate (gaspart_g(1)%pno2t(ntpts))

            if (associated(gaspart_g(1)%ppm25).and.                                        &
                 (.not.associated(gaspart_g(1)%ppm25t)))                                   &
                 allocate (gaspart_g(1)%ppm25t(ntpts))
            
            if (associated(gaspart_g(1)%pco).and.                                          &
                 (.not.associated(gaspart_g(1)%pcot)))                                     &
                 allocate (gaspart_g(1)%pcot(ntpts))

            if (associated(gaspart_g(1)%pso2).and.                                         &
                 (.not.associated(gaspart_g(1)%pso2t)))                                    &
                 allocate (gaspart_g(1)%pso2t(ntpts))
            
            if (associated(gaspart_g(1)%pso4).and.                                         &
                 (.not.associated(gaspart_g(1)%pso4t)))                                    &
                 allocate (gaspart_g(1)%pso4t(ntpts))
            
            if (associated(gaspart_g(1)%paer).and.                                         &
                 (.not.associated(gaspart_g(1)%paert)))                                    &
                 allocate (gaspart_g(1)%paert(ntpts))
            
            if (associated(gaspart_g(1)%pvoc).and.                                         &
                 (.not.associated(gaspart_g(1)%pvoct)))                                    &
                 allocate (gaspart_g(1)%pvoct(ntpts))

            if (ichemi==1) then

               if (associated(gaspart_g(1)%po3).and.                                       &
                    (.not.associated(gaspart_g(1)%po3t)))                                  &
                    allocate (gaspart_g(1)%po3t(ntpts))
               
               if (associated(gaspart_g(1)%prhco).and.                                     &
                    (.not.associated(gaspart_g(1)%prhcot)))                                &
                    allocate (gaspart_g(1)%prhcot(ntpts))
               
               if (associated(gaspart_g(1)%pho2).and.                                      &
                    (.not.associated(gaspart_g(1)%pho2t)))                                 &
                    allocate (gaspart_g(1)%pho2t(ntpts))
               
               if (associated(gaspart_g(1)%po3p).and.                                      &
                    (.not.associated(gaspart_g(1)%po3pt)))                                 &
                    allocate (gaspart_g(1)%po3pt(ntpts))
               
               if (associated(gaspart_g(1)%po1d).and.                                      &
                    (.not.associated(gaspart_g(1)%po1dt)))                                 &
                    allocate (gaspart_g(1)%po1dt(ntpts))
               
               if (associated(gaspart_g(1)%pho).and.                                       &
                    (.not.associated(gaspart_g(1)%phot)))                                  &
                    allocate (gaspart_g(1)%phot(ntpts))
               
               if (associated(gaspart_g(1)%proo).and.                                      &
                    (.not.associated(gaspart_g(1)%proot)))                                 &
                    allocate (gaspart_g(1)%proot(ntpts))
            end if

            do ng=2,ngrs
               gaspart_g(ng)%pnot   => gaspart_g(1)%pnot
               gaspart_g(ng)%pno2t  => gaspart_g(1)%pno2t
               gaspart_g(ng)%ppm25t => gaspart_g(1)%ppm25t
               gaspart_g(ng)%pcot   => gaspart_g(1)%pcot
               gaspart_g(ng)%pso2t  => gaspart_g(1)%pso2t
               gaspart_g(ng)%pso4t  => gaspart_g(1)%pso4t
               gaspart_g(ng)%paert  => gaspart_g(1)%paert
               gaspart_g(ng)%pvoct  => gaspart_g(1)%pvoct
               if (ichemi==1) then
                  gaspart_g(ng)%po3t   => gaspart_g(1)%po3t
                  gaspart_g(ng)%prhcot => gaspart_g(1)%prhcot
                  gaspart_g(ng)%pho2t  => gaspart_g(1)%pho2t
                  gaspart_g(ng)%po3pt  => gaspart_g(1)%po3pt
                  gaspart_g(ng)%po1dt  => gaspart_g(1)%po1dt
                  gaspart_g(ng)%phot   => gaspart_g(1)%phot
                  gaspart_g(ng)%proot  => gaspart_g(1)%proot
               end if
            end do

         end if
         
      end if

      !------ Allocating the tendency arrays for additional scalars. ----------------------!
      do nsc=1,naddsc
         if (       associated(scalar_g(nsc,1)%sclp) .and.                                 &
             (.not. associated(scalar_g(nsc,1)%sclt))     )                                &
              allocate (scalar_g(nsc,1)%sclt(ntpts))
         do ng=2,ngrs
            scalar_g(nsc,ng)%sclt => scalar_g(nsc,1)%sclt
         end do
      end do

      return
   end subroutine alloc_tend
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will reset the pointers of all elements of tend structure.        !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_tend(naddsc)
      use mem_scalar    , only : scalar_g  ! ! data type intent(inout)
      use teb_spm_start , only : teb_spm   ! ! intent(in)
      use mem_gaspart   , only : gaspart_g ! ! data type intent(inout)
      use mem_emiss     , only : ichemi    & ! intent(in)
                               , isource   ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer, intent(in) :: naddsc
      !----- Local Variables. -------------------------------------------------------------!
      integer             :: nsc
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
      if (TEB_SPM==1) then
         if(isource==1)then
            if (associated(gaspart_g(1)%pnot  )) nullify (gaspart_g(1)%pnot  )
            if (associated(gaspart_g(1)%pno2t )) nullify (gaspart_g(1)%pno2t )
            if (associated(gaspart_g(1)%ppm25t)) nullify (gaspart_g(1)%ppm25t)
            if (associated(gaspart_g(1)%pcot  )) nullify (gaspart_g(1)%pcot  )
            if (associated(gaspart_g(1)%pso2t )) nullify (gaspart_g(1)%pso2t )
            if (associated(gaspart_g(1)%pso4t )) nullify (gaspart_g(1)%pso4t )
            if (associated(gaspart_g(1)%paert )) nullify (gaspart_g(1)%paert )
            if (associated(gaspart_g(1)%pvoct )) nullify (gaspart_g(1)%pvoct )

            if(ichemi==1)then
               if (associated(gaspart_g(1)%po3t  )) nullify (gaspart_g(1)%po3t  )
               if (associated(gaspart_g(1)%prhcot)) nullify (gaspart_g(1)%prhcot)
               if (associated(gaspart_g(1)%pho2t )) nullify (gaspart_g(1)%pho2t )
               if (associated(gaspart_g(1)%po3pt )) nullify (gaspart_g(1)%po3pt )
               if (associated(gaspart_g(1)%po1dt )) nullify (gaspart_g(1)%po1dt )
               if (associated(gaspart_g(1)%phot  )) nullify (gaspart_g(1)%phot  )
               if (associated(gaspart_g(1)%proot )) nullify (gaspart_g(1)%proot )
            end if
         end if
      end if
      !----- Additional scalars. ----------------------------------------------------------!
      do nsc=1,naddsc
         if (associated(scalar_g(nsc,1)%sclt)) nullify (scalar_g(nsc,1)%sclt)
      end do
           
      return
   end subroutine nullify_tend
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will deallocate all elements of tend structure.                   !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_tend(naddsc)
      use mem_scalar    , only : scalar_g  ! ! data type intent(inout)
      use teb_spm_start , only : teb_spm   ! ! intent(in)
      use mem_gaspart   , only : gaspart_g ! ! data type intent(inout)
      use mem_emiss     , only : ichemi    & ! intent(in)
                               , isource   ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer, intent(in) :: naddsc
      !----- Local Variables. -------------------------------------------------------------!
      integer             :: nsc
      !------------------------------------------------------------------------------------!

      if (associated(tend%ut  ))  nullify (tend%ut  )
      if (associated(tend%vt  ))  deallocate (tend%vt  )
      if (associated(tend%wt  ))  deallocate (tend%wt  )
      if (associated(tend%pt  ))  deallocate (tend%pt  )
      if (associated(tend%tht ))  deallocate (tend%tht )
      if (associated(tend%rtt ))  deallocate (tend%rtt )
      if (associated(tend%co2t))  deallocate (tend%co2t)
      if (associated(tend%tket))  deallocate (tend%tket)
      if (associated(tend%epst))  deallocate (tend%epst)

      !----- TEB_SPM. ---------------------------------------------------------------------!
      if (TEB_SPM==1) then
         if(isource==1)then
            if (associated(gaspart_g(1)%pnot  )) deallocate (gaspart_g(1)%pnot  )
            if (associated(gaspart_g(1)%pno2t )) deallocate (gaspart_g(1)%pno2t )
            if (associated(gaspart_g(1)%ppm25t)) deallocate (gaspart_g(1)%ppm25t)
            if (associated(gaspart_g(1)%pcot  )) deallocate (gaspart_g(1)%pcot  )
            if (associated(gaspart_g(1)%pso2t )) deallocate (gaspart_g(1)%pso2t )
            if (associated(gaspart_g(1)%pso4t )) deallocate (gaspart_g(1)%pso4t )
            if (associated(gaspart_g(1)%paert )) deallocate (gaspart_g(1)%paert )
            if (associated(gaspart_g(1)%pvoct )) deallocate (gaspart_g(1)%pvoct )

            if(ichemi==1)then
               if (associated(gaspart_g(1)%po3t  )) deallocate (gaspart_g(1)%po3t  )
               if (associated(gaspart_g(1)%prhcot)) deallocate (gaspart_g(1)%prhcot)
               if (associated(gaspart_g(1)%pho2t )) deallocate (gaspart_g(1)%pho2t )
               if (associated(gaspart_g(1)%po3pt )) deallocate (gaspart_g(1)%po3pt )
               if (associated(gaspart_g(1)%po1dt )) deallocate (gaspart_g(1)%po1dt )
               if (associated(gaspart_g(1)%phot  )) deallocate (gaspart_g(1)%phot  )
               if (associated(gaspart_g(1)%proot )) deallocate (gaspart_g(1)%proot )
            end if
         end if
      end if
      !----- Additional scalars. ----------------------------------------------------------!
      do nsc=1,naddsc
         if (associated(scalar_g(nsc,1)%sclt)) deallocate (scalar_g(nsc,1)%sclt)
      end do
           
      return
   end subroutine dealloc_tend
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will fill pointers to arrays into scalar tables.                  !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_tend(basic,micro,turb,scalar,gaspart,naddsc,ng)
      use mem_basic     , only : basic_vars   & ! type
      use mem_micro     , only : micro_vars   & ! type
      use mem_turb      , only : turb_vars    & ! type
      use mem_scalar    , only : scalar_vars  & ! type
      use teb_spm_start , only : TEB_SPM      & ! intent(in)
      use mem_gaspart   , only : gaspart_vars & ! type
      use mem_emiss     , only : ichemi       & ! intent(in)
                               , isource      ! ! intent(in)
      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      type (basic_vars)                , intent(in)  :: basic
      type (micro_vars)                , intent(in)  :: micro
      type (turb_vars)                 , intent(in)  :: turb
      type (scalar_vars) , dimension(:), intent(in)  :: scalar
      type (gaspart_vars)              , pointer     :: gaspart
      integer                          , intent(in)  :: naddsc
      integer                          , intent(in)  :: ng
      !----- Local Variables. -------------------------------------------------------------!
      integer                                        :: elements
      integer                                        :: nsc
      character (len=7)                              :: sname
      !------------------------------------------------------------------------------------!

      if (associated(tend%tht)) then
         call vtables_scalar (basic%thp(1,1,1),tend%tht(1),ng,'THP')
         elements = size(tend%tht)
         call vtables_scalar_new (basic%thp,tend%tht,ng,'THP',elements)
      end if

      if (associated(tend%rtt)) then
         call vtables_scalar (basic%rtp(1,1,1),tend%rtt(1),ng,'RTP')
         elements = size(tend%rtt)
         call vtables_scalar_new (basic%rtp,tend%rtt,ng,'RTP',elements)
      end if

      if (associated(tend%co2t)) then
         call vtables_scalar (basic%co2p(1,1,1),tend%co2t(1),ng,'CO2P')
         elements = size(tend%co2t)
         call vtables_scalar_new (basic%co2p,tend%co2t,ng,'CO2P',elements)
      end if

      if (associated(tend%tket)) then
         call vtables_scalar (turb%tkep(1,1,1),tend%tket(1),ng,'TKEP')
         elements = size(tend%tket)
         call vtables_scalar_new (turb%tkep,tend%tket,ng,'TKEP',elements)
      end if

      if (associated(tend%epst)) then
         call vtables_scalar (turb%epsp(1,1,1),tend%epst(1),ng,'EPSP')
         elements = size(tend%epst)
         call vtables_scalar_new (turb%epsp,tend%epst,ng,'EPSP',elements)
      end if

      !----- TEB_SPM. ---------------------------------------------------------------------!
      if (TEB_SPM==1) then
         if (isource==1) then

            if (associated(gaspart%pnot)) then
               call vtables_scalar (gaspart%pno(1,1,1),gaspart%pnot(1),ng,'PNO')
               elements = size(gaspart%pno)
               call vtables_scalar_new (gaspart%p,gaspart%pnot,ng,'PNO',elements)
            end if

            if (associated(gaspart%pno2t)) then
               call vtables_scalar (gaspart%pno2(1,1,1),gaspart%pno2t(1),ng,'PNO2')
               elements = size(gaspart%pno2)
               call vtables_scalar_new (gaspart%pno2,gaspart%pno2t,ng,'PNO2',elements)
            end if

            if (associated(gaspart%ppm25t)) then
               call vtables_scalar (gaspart%ppm25(1,1,1),gaspart%ppm25t(1),ng,'PM25')
               call vtables_scalar_new (gaspart%ppm25,gaspart%ppm25t,ng,'PM25',elements)
            end if

            if (associated(gaspart%pcot)) then
               call vtables_scalar (gaspart%pco(1,1,1),gaspart%pcot(1),ng,'PCO')
               call vtables_scalar_new (gaspart%pco,gaspart%pcot,ng,'PCO',elements)
            end if
            
            if (associated(gaspart%pso2t)) then
               call vtables_scalar (gaspart%pso2(1,1,1),gaspart%pso2t(1),ng,'PSO2')
               call vtables_scalar_new (gaspart%pso2,gaspart%pso2t,ng,'PSO2',elements)
            end if

            if (associated(gaspart%pso4t)) then
               call vtables_scalar (gaspart%pso4(1,1,1),gaspart%pso4t(1),ng,'PSO4')
               call vtables_scalar_new (gaspart%pso4,gaspart%pso4t,ng,'PSO4',elements)
            end if

            if (associated(gaspart%paert)) then
               call vtables_scalar (gaspart%paer(1,1,1),gaspart%paert(1),ng,'PAER')
               call vtables_scalar_new (gaspart%paer,gaspart%paert,ng,'PAER',elements)
            end if

            if (associated(gaspart%pvoct)) then
               call vtables_scalar (gaspart%pvoc(1,1,1),gaspart%pvoct(1),ng,'PVOC')
               call vtables_scalar_new (gaspart%pvoc,gaspart%pvoct,ng,'PVOC',elements)
            end if

            if(ichemi==1) then
               if (associated(gaspart%po3t)) then
                  call vtables_scalar(gaspart%po3(1,1,1),gaspart%po3t(1),ng,'PO3')
                  call vtables_scalar_new(gaspart%po3,gaspart%po3t,ng,'PO3',elements)
               end if

               if (associated(gaspart%prhcot)) then
                  call vtables_scalar(gaspart%prhco(1,1,1),gaspart%prhcot(1),ng,'PRHCO')
                  call vtables_scalar_new(gaspart%prhco,gaspart%prhcot,ng,'PRHCO',elements)
               end if

               if (associated(gaspart%pho2t)) then
                  call vtables_scalar (gaspart%pho2(1,1,1),gaspart%pho2t(1),ng,'PHO2')
                  call vtables_scalar_new (gaspart%pho2,gaspart%pho2t,ng,'PHO2',elements)
               end if

               if (associated(gaspart%po3pt)) then
                  call vtables_scalar (gaspart%po3p(1,1,1),gaspart%po3pt(1),ng,'PO3P')
                  call vtables_scalar_new (gaspart%po3p,gaspart%po3pt,ng,'PO3P',elements)
               end if

               if (associated(gaspart%po1dt)) then
                  call vtables_scalar (gaspart%po1d(1,1,1),gaspart%po1dt(1),ng,'PO1D')
                  call vtables_scalar_new (gaspart%po1d,gaspart%po1dt,ng,'PO1D',elements)
               end if

               if (associated(gaspart%phot)) then
                  call vtables_scalar (gaspart%pho(1,1,1),gaspart%phot(1),ng,'PHO')
                  call vtables_scalar_new (gaspart%pho,gaspart%phot,ng,'PHO',elements)
               end if

               if (associated(gaspart%proot)) then
                  call vtables_scalar (gaspart%proo(1,1,1),gaspart%proot(1),ng,'PROO')
                  call vtables_scalar_new (gaspart%proo,gaspart%proot,ng,'PROO',elements)
               end if
            end if
         end if
      end if

      do nsc=1,naddsc
         write(sname,'(a4,i3.3)') 'SCLP',nsc
         if (associated(scalar(nsc)%sclt)) then
            call vtables_scalar (scalar(nsc)%sclp(1,1,1),scalar(nsc)%sclt(1),ng,sname)
            call vtables_scalar_new (scalar(nsc)%sclp,scalar(nsc)%sclt,ng,sname,elements)
         end if
      end do
   end subroutine filltab_tend
   !=======================================================================================!
   !=======================================================================================!
end module mem_tend
!==========================================================================================!
!==========================================================================================!

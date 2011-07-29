!==========================================================================================!
!==========================================================================================!
!     This module contains the structures that holds the information that is exchanged     !
! between ED and BRAMS, plus any variable that may be used to set up the coupled runs that !
! should not be entirely part of ED, neither are used in any other part of BRAMS.          !
!------------------------------------------------------------------------------------------!
module mem_edcp

   !---------------------------------------------------------------------------------------!
   !    This variables are used to change the CO2 that ED receives and what is sent back   !
   ! to BRAMS.                                                                             !
   !---------------------------------------------------------------------------------------!
   real(kind=4)            :: co2_offset
   real(kind=4), parameter :: atm_co2_min = 70.
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    This variables contains the time for interpolation of the ED variables.            !
   !---------------------------------------------------------------------------------------!
   real(kind=8)            :: edtime1
   real(kind=8)            :: edtime2
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     This structure holds some fast scale flux variables.  They may not be updated     !
   ! every BRAMS time step, so an interpolation may be needed.  Therefore, we store two of !
   ! these structures, ed_fluxp_g, with the "past" information, and ed_fluxf_g, with the   !
   ! "future" information.  Only one structure is needed for water sites, as they are      !
   ! called every BRAMS time step.                                                         !
   !---------------------------------------------------------------------------------------!
   type ed_flux
      real(kind=4), dimension (:,:,:), pointer :: ustar
      real(kind=4), dimension (:,:,:), pointer :: tstar
      real(kind=4), dimension (:,:,:), pointer :: rstar
      real(kind=4), dimension (:,:,:), pointer :: cstar
      real(kind=4), dimension (:,:,:), pointer :: zeta
      real(kind=4), dimension (:,:,:), pointer :: ribulk
      real(kind=4), dimension (:,:,:), pointer :: sflux_u
      real(kind=4), dimension (:,:,:), pointer :: sflux_v
      real(kind=4), dimension (:,:,:), pointer :: sflux_t
      real(kind=4), dimension (:,:,:), pointer :: sflux_r
      real(kind=4), dimension (:,:,:), pointer :: sflux_c
      real(kind=4), dimension (:,:,:), pointer :: sflux_w
      real(kind=4), dimension (:,:,:), pointer :: albedt
      real(kind=4), dimension (:,:,:), pointer :: rlongup
      real(kind=4), dimension (:,:,:), pointer :: rk4step
   end type ed_flux

   type(ed_flux),pointer, dimension(:) :: ed_fluxp_g
   type(ed_flux),pointer, dimension(:) :: ed_fluxpm_g
   type(ed_flux),pointer, dimension(:) :: ed_fluxf_g
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     These variables are needed because the microphysics, cumulus, and ED have         !
   ! different update frequency amongst themselves.  We store the total accumulated        !
   ! precipitation at the previous ED call, so we can find the mean precipitation rate.    !
   !---------------------------------------------------------------------------------------!
   type ed_precip
      real(kind=4), dimension(:,:), pointer :: prev_aconpr
      real(kind=4), dimension(:,:), pointer :: prev_abulkpr
   end type ed_precip

   type (ed_precip), pointer, dimension(:) :: ed_precip_g
   type (ed_precip), pointer, dimension(:) :: ed_precipm_g ! Not that we really need, but
                                                           !    we add this to make sure
                                                           !    that the model will work.
   !---------------------------------------------------------------------------------------!



   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine allocates the flux structures.                                    !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_edflux(edflux,nxm,nym,nsite)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ed_flux), intent(inout) :: edflux
      integer      , intent(in)    :: nxm
      integer      , intent(in)    :: nym
      integer      , intent(in)    :: nsite
      !------------------------------------------------------------------------------------!


      !----- Nullify the structure first, for a safe allocation. --------------------------!
      call nullify_edflux(edflux)
      !------------------------------------------------------------------------------------!


      !----- Allocate the elements of the ed flux structure. ------------------------------!
      allocate(edflux%ustar     (nxm,nym,nsite)   )
      allocate(edflux%tstar     (nxm,nym,nsite)   )
      allocate(edflux%rstar     (nxm,nym,nsite)   )
      allocate(edflux%cstar     (nxm,nym,nsite)   )
      allocate(edflux%zeta      (nxm,nym,nsite)   )
      allocate(edflux%ribulk    (nxm,nym,nsite)   )
      allocate(edflux%sflux_u   (nxm,nym,nsite)   )
      allocate(edflux%sflux_v   (nxm,nym,nsite)   )
      allocate(edflux%sflux_r   (nxm,nym,nsite)   )
      allocate(edflux%sflux_c   (nxm,nym,nsite)   )
      allocate(edflux%sflux_t   (nxm,nym,nsite)   )
      allocate(edflux%sflux_w   (nxm,nym,nsite)   )
      allocate(edflux%albedt    (nxm,nym,nsite)   )
      allocate(edflux%rlongup   (nxm,nym,nsite)   )
      allocate(edflux%rk4step   (nxm,nym,nsite)   )

      return
   end subroutine alloc_edflux
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine nullifies the elements of the ed flux structure, for a safe      !
   ! allocation using pointers.                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_edflux(edflux)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ed_flux), intent(inout) :: edflux
      !------------------------------------------------------------------------------------!
     
      if (associated(edflux%ustar     ))  nullify(edflux%ustar     )
      if (associated(edflux%tstar     ))  nullify(edflux%tstar     )
      if (associated(edflux%rstar     ))  nullify(edflux%rstar     )
      if (associated(edflux%cstar     ))  nullify(edflux%cstar     )
      if (associated(edflux%zeta      ))  nullify(edflux%zeta      )
      if (associated(edflux%ribulk    ))  nullify(edflux%ribulk    )
      if (associated(edflux%sflux_u   ))  nullify(edflux%sflux_u   )
      if (associated(edflux%sflux_v   ))  nullify(edflux%sflux_v   )
      if (associated(edflux%sflux_r   ))  nullify(edflux%sflux_r   )
      if (associated(edflux%sflux_t   ))  nullify(edflux%sflux_t   )
      if (associated(edflux%sflux_c   ))  nullify(edflux%sflux_c   )
      if (associated(edflux%sflux_w   ))  nullify(edflux%sflux_w   )
      if (associated(edflux%albedt    ))  nullify(edflux%albedt    )
      if (associated(edflux%rlongup   ))  nullify(edflux%rlongup   )
      if (associated(edflux%rk4step   ))  nullify(edflux%rk4step   )

      return
   end subroutine nullify_edflux
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine assigns an initial value of zero for all the flux structures.     !
   !---------------------------------------------------------------------------------------!
   subroutine zero_edflux(edflux,ust0,dtref)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ed_flux), intent(inout) :: edflux
      real         , intent(in)    :: ust0
      real         , intent(in)    :: dtref
      !------------------------------------------------------------------------------------!

      if (associated(edflux%ustar     ))  edflux%ustar   = ust0
      if (associated(edflux%tstar     ))  edflux%tstar   = 0.0
      if (associated(edflux%rstar     ))  edflux%rstar   = 0.0
      if (associated(edflux%cstar     ))  edflux%cstar   = 0.0
      if (associated(edflux%zeta      ))  edflux%zeta    = 0.0
      if (associated(edflux%ribulk    ))  edflux%ribulk  = 0.0
      if (associated(edflux%sflux_u   ))  edflux%sflux_u = 0.0
      if (associated(edflux%sflux_v   ))  edflux%sflux_v = 0.0
      if (associated(edflux%sflux_r   ))  edflux%sflux_r = 0.0
      if (associated(edflux%sflux_t   ))  edflux%sflux_t = 0.0
      if (associated(edflux%sflux_c   ))  edflux%sflux_c = 0.0
      if (associated(edflux%sflux_w   ))  edflux%sflux_w = 0.0
      if (associated(edflux%albedt    ))  edflux%albedt  = 0.0
      if (associated(edflux%rlongup   ))  edflux%rlongup = 0.0
      if (associated(edflux%rk4step   ))  edflux%rk4step = dtref

      return
   end subroutine zero_edflux
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine frees the memory associated with the pointers, and should be      !
   ! called before the de-allocation of the actual structure.                              !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_edflux(edflux)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ed_flux), intent(inout) :: edflux
      !------------------------------------------------------------------------------------!
     
      if (associated(edflux%ustar     ))  deallocate(edflux%ustar     )
      if (associated(edflux%tstar     ))  deallocate(edflux%tstar     )
      if (associated(edflux%rstar     ))  deallocate(edflux%rstar     )
      if (associated(edflux%cstar     ))  deallocate(edflux%cstar     )
      if (associated(edflux%zeta      ))  deallocate(edflux%zeta      )
      if (associated(edflux%ribulk    ))  deallocate(edflux%ribulk    )
      if (associated(edflux%sflux_u   ))  deallocate(edflux%sflux_u   )
      if (associated(edflux%sflux_v   ))  deallocate(edflux%sflux_v   )
      if (associated(edflux%sflux_r   ))  deallocate(edflux%sflux_r   )
      if (associated(edflux%sflux_c   ))  deallocate(edflux%sflux_c   )
      if (associated(edflux%sflux_t   ))  deallocate(edflux%sflux_t   )
      if (associated(edflux%sflux_w   ))  deallocate(edflux%sflux_w   )
      if (associated(edflux%albedt    ))  deallocate(edflux%albedt    )
      if (associated(edflux%rlongup   ))  deallocate(edflux%rlongup   )
      if (associated(edflux%rk4step   ))  deallocate(edflux%rk4step   )

      return
   end subroutine dealloc_edflux
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the variable table with the previous precipitation values  !
   ! for ED.  This is because we need the information when we run with HISTORY.            !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_edflux(edflux,edfluxm,imean,nxm,nym,npatm,ng)
      use var_tables

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ed_flux), intent(inout) :: edflux
      type(ed_flux), intent(inout) :: edfluxm
      integer      , intent(in)    :: imean
      integer      , intent(in)    :: nxm
      integer      , intent(in)    :: nym
      integer      , intent(in)    :: npatm
      integer      , intent(in)    :: ng
      !----- Local variables. -------------------------------------------------------------!
      integer                      :: npts
      !------------------------------------------------------------------------------------!

      !----- Fill pointers to arrays into variable tables. --------------------------------!
      npts = nxm * nym * npatm

      if (associated(edflux%rk4step )) then
         call vtables2 (edflux%rk4step,edfluxm%rk4step,ng,npts,imean                       &
                       , 'RK4STEP :6:hist:mpti:mpt3')
      end if

      return
   end subroutine filltab_edflux
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine allocates the precipitation structures.                           !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_edprecip(edprec,nxm,nym)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ed_precip), intent(inout) :: edprec
      integer        , intent(in)    :: nxm
      integer        , intent(in)    :: nym
      !------------------------------------------------------------------------------------!


      !----- Nullify the structure first, for a safe allocation. --------------------------!
      call nullify_edprecip(edprec)
      !------------------------------------------------------------------------------------!


      !----- Allocate the elements of the ed precipitation structure. ---------------------!
      allocate(edprec%prev_aconpr (nxm,nym))
      allocate(edprec%prev_abulkpr(nxm,nym))

      return
   end subroutine alloc_edprecip
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine nullifies the elements of the ed precipitation structure, for a  !
   ! safe allocation using pointers.                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_edprecip(edprec)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ed_precip), intent(inout) :: edprec
      !------------------------------------------------------------------------------------!



      if (associated(edprec%prev_aconpr )) nullify(edprec%prev_aconpr )
      if (associated(edprec%prev_abulkpr)) nullify(edprec%prev_abulkpr)

      return
   end subroutine  nullify_edprecip
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine assigns an initial value of zero for the precipitation            !
   !  structure.                                                                           !
   !---------------------------------------------------------------------------------------!
   subroutine zero_edprecip(edprec)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ed_precip), intent(inout) :: edprec
      !------------------------------------------------------------------------------------!

      if (associated(edprec%prev_aconpr )) edprec%prev_aconpr  = 0.0
      if (associated(edprec%prev_abulkpr)) edprec%prev_abulkpr = 0.0

      return
   end subroutine  zero_edprecip
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine frees the memory associated with the pointers, and should be      !
   ! called before the de-allocation of the actual structure.                              !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_ed_precip(edprec)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ed_precip), intent(inout) :: edprec
      !------------------------------------------------------------------------------------!
      
      if (associated(edprec%prev_aconpr )) deallocate(edprec%prev_aconpr )
      if (associated(edprec%prev_abulkpr)) deallocate(edprec%prev_abulkpr)

      return
   end subroutine  dealloc_ed_precip
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This routine will fill the variable table with the previous precipitation values  !
   ! for ED.  This is because we need the information when we run with HISTORY.            !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_ed_precip(edprec,edprecm,imean,nxm,nym,ng)
      use var_tables

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(ed_precip), intent(inout) :: edprec
      type(ed_precip), intent(inout) :: edprecm
      integer        , intent(in)    :: imean
      integer        , intent(in)    :: nxm
      integer        , intent(in)    :: nym
      integer        , intent(in)    :: ng
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: npts
      !------------------------------------------------------------------------------------!

      !----- Fill pointers to arrays into variable tables. --------------------------------!
      npts = nxm * nym

      if (associated(edprec%prev_aconpr )) then
         call vtables2 (edprec%prev_aconpr,edprecm%prev_aconpr,ng,npts,imean               &
                       , 'PREV_ACONPR :2:hist:mpti:mpt3')
      end if
      if (associated(edprec%prev_abulkpr )) then
         call vtables2 (edprec%prev_abulkpr,edprecm%prev_abulkpr,ng,npts,imean             &
                       , 'PREV_ABULKPR :2:hist:mpti:mpt3')
      end if

      return
   end subroutine filltab_ed_precip
   !=======================================================================================!
   !=======================================================================================!
end module mem_edcp
!==========================================================================================!
!==========================================================================================!

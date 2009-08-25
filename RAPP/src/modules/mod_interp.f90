!==========================================================================================!
!==========================================================================================!
!    This module contains some interpolation parameters.  The interpolation method is bas- !
! ed on an objective analysis:                                                             !
!     Koch, S. E., M. desJardins, and P.J. Kocin, 1983: An interactive Barnes objective    !
! analysis for use with satellite and conventional data.                                   !
!------------------------------------------------------------------------------------------!
module mod_interp
   implicit none
   
   !---------------------------------------------------------------------------------------!
   !    Namelist variables                                                                 !
   !---------------------------------------------------------------------------------------!
   real         :: dtinc      ! Time step of integration of the radiation fluxes.
   real         :: gamma0     ! Parameter gamma.
   real         :: minweight  ! Minimum weight we will consider for the objective analysis
   !---------------------------------------------------------------------------------------!


   !----- Namelist based variables. -------------------------------------------------------!
   real         :: dtinc_fday ! Time step, in fraction of a day.
   real(kind=8) :: minweight8 ! Double precision version of minweight
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Parameters for precipitation downscaling.                                         !
   !---------------------------------------------------------------------------------------!
   !----- Number of elements in the local precipitation sampling CDF and inverse. ---------!
   integer, parameter :: nlocpcp  = 250
   real   , parameter :: nlocpcpi = 1. / real(nlocpcp)
   !----- Maximum precipitation rate. -----------------------------------------------------!
   real   , parameter :: max_local_precip = 1.e-2
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Indexing variables.                                                                !
   !---------------------------------------------------------------------------------------!
   integer :: mxgauss
   integer :: mygauss
   integer :: mxlola
   integer :: mylola
   integer :: xla
   integer :: xlz
   integer :: yla
   integer :: ylz
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Interpolation buffer structure.                                                    !
   !---------------------------------------------------------------------------------------!
   type interp_vars
      !----- Variables to be dimensioned by (NX_Gauss;NY_Gauss). --------------------------!
      real   , pointer, dimension(:,:)   :: deltan ! Grid resolution at this point.
      real   , pointer, dimension(:,:)   :: kappa0 ! Grid-dependent parameter kappa 0
      real   , pointer, dimension(:,:)   :: g0     ! The 1st pass of the intepolated field 
      !----- Variables to be dimensioned by (NX_Lon/Lat;NY_Lon/Lat). ----------------------!
      real   , pointer, dimension(:,:)   :: rm2    ! Square of distance between g and f
      real   , pointer, dimension(:,:)   :: weight !Weight of the objective analysis
      logical, pointer, dimension(:,:)   :: mask   ! Mask for skipping points far away.
      !----- Residual of 1st pass, it must be a (NX_Lon/Lat;NY_Lon/Lat;maxntimes). --------!
      real   , pointer, dimension(:,:,:) :: residu ! Residual after the first pass.
      !----- Variables for Gaussian-Lon/Lat grid matching. --------------------------------!
      integer, pointer, dimension(:,:) :: iwe      ! Corresponding i to the west
      integer, pointer, dimension(:,:) :: jno      ! Corresponding j to the north
      real   , pointer, dimension(:,:) :: dix      ! Distance from point to iwe
      real   , pointer, dimension(:,:) :: diy      ! Distance from point to jno
      !----- Variable for radiation interpolation. ----------------------------------------!
      real   , pointer, dimension(:)   :: zen_norm ! Normalised zenithal angle
      !------------------------------------------------------------------------------------!
   end type interp_vars

   type(interp_vars) :: interp_buffer
   !---------------------------------------------------------------------------------------!



   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_interp(interp,mtp)
      use mod_ioopts, only : radratio
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer          , intent(in)    :: mtp
      type(interp_vars), intent(inout) :: interp
      !------------------------------------------------------------------------------------!

      allocate(interp%deltan   (mxgauss , mygauss           ))
      allocate(interp%kappa0   (mxgauss , mygauss           ))
      allocate(interp%g0       (mxgauss , mygauss           ))

      allocate(interp%weight   (mxlola  , mylola            ))
      allocate(interp%rm2      (mxlola  , mylola            ))
      allocate(interp%mask     (mxlola  , mylola            ))

      allocate(interp%residu   (mxlola  , mylola  , mtp     ))
      
      allocate(interp%iwe      (mxlola  , mylola            ))
      allocate(interp%jno      (mxlola  , mylola            ))
      allocate(interp%dix      (mxlola  , mylola            ))
      allocate(interp%diy      (mxlola  , mylola            ))
      
      allocate(interp%zen_norm (                    radratio))

      return
   end subroutine alloc_interp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine nullify_interp(interp)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(interp_vars), intent(inout) :: interp  
      !------------------------------------------------------------------------------------!
      
      if (associated(interp%deltan  )) nullify(interp%deltan  )
      if (associated(interp%kappa0  )) nullify(interp%kappa0  )
      if (associated(interp%g0      )) nullify(interp%g0      )
      if (associated(interp%weight  )) nullify(interp%weight  )
      if (associated(interp%rm2     )) nullify(interp%rm2     )
      if (associated(interp%mask    )) nullify(interp%mask    )
      if (associated(interp%residu  )) nullify(interp%residu  )
      if (associated(interp%iwe     )) nullify(interp%iwe     )
      if (associated(interp%jno     )) nullify(interp%jno     )
      if (associated(interp%dix     )) nullify(interp%dix     )
      if (associated(interp%diy     )) nullify(interp%diy     )
      if (associated(interp%zen_norm)) nullify(interp%zen_norm)

      return
   end subroutine nullify_interp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine init_interp(interp,fullinit)
      use mod_ioopts , only : missflg_real & ! intent(in)
                            , missflg_int  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(interp_vars), intent(inout) :: interp  
      logical          , intent(in)    :: fullinit ! Flag for full initialisation.
      !------------------------------------------------------------------------------------!
      
      if (associated(interp%kappa0  )) interp%kappa0   = missflg_real
      if (associated(interp%g0      )) interp%g0       = missflg_real
      if (associated(interp%weight  )) interp%weight   = missflg_real
      if (associated(interp%rm2     )) interp%rm2      = missflg_real
      if (associated(interp%mask    )) interp%mask     = .false.
      if (associated(interp%residu  )) interp%residu   = missflg_real
      if (associated(interp%zen_norm)) interp%zen_norm = missflg_real

      !----- The grid resolution and matching information don't change during the run. ----!
      if (fullinit) then
         if (associated(interp%deltan )) interp%deltan = missflg_real

         if (associated(interp%iwe    )) interp%iwe    = missflg_int
         if (associated(interp%jno    )) interp%jno    = missflg_int
         if (associated(interp%dix    )) interp%dix    = missflg_real
         if (associated(interp%diy    )) interp%diy    = missflg_real
      end if


      return
   end subroutine init_interp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine dealloc_interp(interp,fulldeal)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(interp_vars), intent(inout) :: interp  
      logical          , intent(in)    :: fulldeal ! Flag for full initialisation.
      !------------------------------------------------------------------------------------!
      
      if (associated(interp%deltan  )) deallocate(interp%deltan  )
      if (associated(interp%kappa0  )) deallocate(interp%kappa0  )
      if (associated(interp%g0      )) deallocate(interp%g0      )
      if (associated(interp%weight  )) deallocate(interp%weight  )
      if (associated(interp%rm2     )) deallocate(interp%rm2     )
      if (associated(interp%mask    )) deallocate(interp%mask    )
      if (associated(interp%residu  )) deallocate(interp%residu  )
      if (associated(interp%iwe     )) deallocate(interp%iwe     )
      if (associated(interp%jno     )) deallocate(interp%jno     )
      if (associated(interp%dix     )) deallocate(interp%dix     )
      if (associated(interp%diy     )) deallocate(interp%diy     )
      if (associated(interp%zen_norm)) deallocate(interp%zen_norm)


      return
   end subroutine dealloc_interp
   !=======================================================================================!
   !=======================================================================================!
end module mod_interp
!==========================================================================================!
!==========================================================================================!


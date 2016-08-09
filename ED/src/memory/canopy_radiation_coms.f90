!==========================================================================================!
!==========================================================================================!
!    This module contains several parameters used in the canopy radiation solver.          !
!                                                                                          !
!  IMPORTANT: Do not initialize non-parameters in their modules - not all compilers will   !
!             actually initialize them.  Instead, assign them at init_can_rad_params sub-  !
!             routine (ed_params.f90).                                                     !
!------------------------------------------------------------------------------------------!
module canopy_radiation_coms

   use ed_max_dims , only : n_pft
   implicit none 


   !---------------------------------------------------------------------------------------!
   ! ICANRAD -- Specifies how canopy radiation is solved.  This variable sets both short-  !
   !            wave and longwave.                                                         !
   !            0.  Two-stream model (Medvigy 2006), with the possibility to apply         !
   !                finite crown area to direct shortwave radiation.                       !
   !            1.  Multiple-scattering model (Zhao and Qualls 2005,2006), with the        !
   !                possibility to apply finite crown area to all radiation fluxes.        !
   !---------------------------------------------------------------------------------------!
   integer :: icanrad
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     The following variables are temporary namelist variables used to control the      !
   ! radiation properties of leaves.                                                       !
   ! LTRANS_VIS   -- Leaf transmittance on visible.                                        !
   ! LTRANS_NIR   -- Leaf transmittance on near infrared.                                  !
   ! LREFLECT_VIS -- Leaf reflectance on visible.                                          !
   ! LREFLECT_NIR -- Leaf reflectance on near infrared.                                    !
   ! ORIENT_TREE  -- Leaf orientation parameter for tropical trees                         !
   ! ORIENT_GRASS -- Leaf orientation parameter for tropical grasses                       !
   ! CLUMP_TREE   -- Leaf clumping factor for tropical trees                               !
   ! CLUMP_GRASS  -- Leaf clumping factor for tropical grasses                             !
   !---------------------------------------------------------------------------------------!
   real :: ltrans_vis
   real :: ltrans_nir
   real :: lreflect_vis
   real :: lreflect_nir
   real :: orient_tree
   real :: orient_grass
   real :: clump_tree
   real :: clump_grass
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Fraction of diffuse solar radiation in the PAR band.  Used when you don't know    !
   ! the direct/diffuse breakdown. (parameters)                                            !
   !---------------------------------------------------------------------------------------!
   real :: fvis_beam_def
   real :: fvis_diff_def
   real :: fnir_beam_def
   real :: fnir_diff_def
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Structure with scratch variables for radiation (Thanks RGK!).                      !
   !---------------------------------------------------------------------------------------!
   type radscrtype
      integer         , pointer, dimension(:)   :: pft_array
      real(kind=8)    , pointer, dimension(:)   :: leaf_temp_array
      real(kind=8)    , pointer, dimension(:)   :: wood_temp_array
      real(kind=8)    , pointer, dimension(:)   :: lai_array
      real(kind=8)    , pointer, dimension(:)   :: wai_array 
      real(kind=8)    , pointer, dimension(:)   :: CA_array
      real(kind=8)    , pointer, dimension(:)   :: htop_array
      real(kind=8)    , pointer, dimension(:)   :: hbot_array
      real(kind=8)    , pointer, dimension(:)   :: par_level_beam
      real(kind=8)    , pointer, dimension(:)   :: par_level_diffd
      real(kind=8)    , pointer, dimension(:)   :: par_level_diffu
      real(kind=8)    , pointer, dimension(:)   :: light_level_array
      real(kind=8)    , pointer, dimension(:)   :: light_beam_level_array
      real(kind=8)    , pointer, dimension(:)   :: light_diff_level_array
      real            , pointer, dimension(:)   :: par_v_beam_array
      real            , pointer, dimension(:)   :: rshort_v_beam_array
      real            , pointer, dimension(:)   :: par_v_diffuse_array
      real            , pointer, dimension(:)   :: rshort_v_diffuse_array
      real            , pointer, dimension(:)   :: lw_v_array
      real            , pointer, dimension(:,:) :: radprof_array
   end type radscrtype
   type(radscrtype)   , pointer,dimension(:)    :: radscr(:)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     These are the normalised variables that will be used in the two-stream model.     !
   !---------------------------------------------------------------------------------------!
   real(kind=8) :: par_beam_norm
   real(kind=8) :: par_diff_norm
   real(kind=8) :: nir_beam_norm
   real(kind=8) :: nir_diff_norm
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Factors that define the orientation and clumping of leaves.                       !
   ! CLUMPING FACTOR - factor indicating the degree of clumpiness of leaves.               !
   ! ORIENT_FACTOR   - mean leaf orientation.                                              !
   !                     0 -- leaves are randomly oriented                                 !
   !                     1 -- all leaves are perfectly horizontal                          !
   !                    -1 -- all leaves are perfectly vertical.                           !
   ! PHI1            - The phi1 term from the CLM technical manual                         !
   ! PHI2            - The phi2 term from the CLM technical manual                         !
   ! MU_BAR          - average cosine of incidence angle for hemispheric (diffuse)         !
   !                   radiation (for both short wave and long wave)                       !
   !---------------------------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: clumping_factor
   real(kind=8), dimension(n_pft) :: orient_factor
   real(kind=8), dimension(n_pft) :: phi1
   real(kind=8), dimension(n_pft) :: phi2
   real(kind=8), dimension(n_pft) :: mu_bar
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Reflectance coefficients.                                                         !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_reflect_vis
   real(kind=8), dimension(n_pft) :: wood_reflect_vis
   !----- Near infrared. ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_reflect_nir
   real(kind=8), dimension(n_pft) :: wood_reflect_nir
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Transmittance coefficients.                                                       !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_trans_vis
   real(kind=8), dimension(n_pft) :: wood_trans_vis
   !----- Near infrared. ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_trans_nir
   real(kind=8), dimension(n_pft) :: wood_trans_nir
   !---------------------------------------------------------------------------------------!




   !----- Emissivity of the vegetation (TIR). ---------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_emiss_tir
   real(kind=8), dimension(n_pft) :: wood_emiss_tir
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Scattering coefficients.                                                          !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_scatter_vis
   real(kind=8), dimension(n_pft) :: wood_scatter_vis
   !----- Near infrared. ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_scatter_nir
   real(kind=8), dimension(n_pft) :: wood_scatter_nir
   !----- Thermal infrared. ---------------------------------------------------------------!
   !real(kind=8), dimension(n_pft) :: leaf_scatter_tir
   !real(kind=8), dimension(n_pft) :: wood_scatter_tir
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Fraction of diffuse radiation that is upscattered.                                !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_backscatter_vis
   real(kind=8), dimension(n_pft) :: wood_backscatter_vis
   !----- Near infrared. ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_backscatter_nir
   real(kind=8), dimension(n_pft) :: wood_backscatter_nir
   !----- Backscattering of thermal infrared. ---------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_backscatter_tir
   real(kind=8), dimension(n_pft) :: wood_backscatter_tir
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Snow pack properties.                                                              !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: snow_albedo_vis
   real(kind=4) :: snow_albedo_nir
   real(kind=4) :: snow_emiss_tir
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     The following variables control whether to call things that should be called      !
   ! when there is still some light.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4)    :: rshort_twilight_min
   real(kind=4)    :: cosz_min
   real(kind=8)    :: cosz_min8
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine allocates the scratch variables after all pointers have been     !
   ! nullified.                                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_radscratch(cradscr,maxcohort)
      
      use ed_max_dims          , only : n_radprof            ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(radscrtype), target     :: cradscr
      integer         , intent(in) :: maxcohort
      !------------------------------------------------------------------------------------!

      call nullify_radscratch(cradscr)

      allocate (cradscr%pft_array                 (          maxcohort))
      allocate (cradscr%leaf_temp_array           (          maxcohort))
      allocate (cradscr%wood_temp_array           (          maxcohort))
      allocate (cradscr%lai_array                 (          maxcohort))
      allocate (cradscr%wai_array                 (          maxcohort))
      allocate (cradscr%CA_array                  (          maxcohort))
      allocate (cradscr%htop_array                (          maxcohort))
      allocate (cradscr%hbot_array                (          maxcohort))
      allocate (cradscr%par_level_beam            (          maxcohort))
      allocate (cradscr%par_level_diffu           (          maxcohort))
      allocate (cradscr%par_level_diffd           (          maxcohort))
      allocate (cradscr%light_level_array         (          maxcohort))
      allocate (cradscr%light_beam_level_array    (          maxcohort))
      allocate (cradscr%light_diff_level_array    (          maxcohort))
      allocate (cradscr%par_v_beam_array          (          maxcohort))
      allocate (cradscr%rshort_v_beam_array       (          maxcohort))
      allocate (cradscr%par_v_diffuse_array       (          maxcohort))
      allocate (cradscr%rshort_v_diffuse_array    (          maxcohort))
      allocate (cradscr%lw_v_array                (          maxcohort))
      allocate (cradscr%radprof_array             (n_radprof,maxcohort))
      
      return
   end subroutine alloc_radscratch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine de-allocates the scratch variables.                              !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_radscratch(cradscr)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(radscrtype), target :: cradscr
      !------------------------------------------------------------------------------------!

      if (associated(cradscr%pft_array             )) deallocate(cradscr%pft_array        )
      if (associated(cradscr%leaf_temp_array       )) deallocate(cradscr%leaf_temp_array  )
      if (associated(cradscr%wood_temp_array       )) deallocate(cradscr%wood_temp_array  )
      if (associated(cradscr%lai_array             )) deallocate(cradscr%lai_array        )
      if (associated(cradscr%wai_array             )) deallocate(cradscr%wai_array        )
      if (associated(cradscr%CA_array              )) deallocate(cradscr%CA_array         )
      if (associated(cradscr%htop_array            )) deallocate(cradscr%htop_array       )
      if (associated(cradscr%hbot_array            )) deallocate(cradscr%hbot_array       )
      if (associated(cradscr%par_level_beam        )) deallocate(cradscr%par_level_beam   )
      if (associated(cradscr%par_level_diffu       )) deallocate(cradscr%par_level_diffu  )
      if (associated(cradscr%par_level_diffd       )) deallocate(cradscr%par_level_diffd  )
      if (associated(cradscr%light_level_array     )) deallocate(cradscr%light_level_array)
      if (associated(cradscr%light_beam_level_array))                                      &
                                                 deallocate(cradscr%light_beam_level_array)
      if (associated(cradscr%light_diff_level_array))                                      &
                                                 deallocate(cradscr%light_diff_level_array)
      if (associated(cradscr%par_v_beam_array      )) deallocate(cradscr%par_v_beam_array )
      if (associated(cradscr%rshort_v_beam_array   ))                                      &
                                                 deallocate(cradscr%rshort_v_beam_array   )
      if (associated(cradscr%par_v_diffuse_array   ))                                      &
                                                 deallocate(cradscr%par_v_diffuse_array   )
      if (associated(cradscr%rshort_v_diffuse_array))                                      &
                                                 deallocate(cradscr%rshort_v_diffuse_array)
      if (associated(cradscr%lw_v_array            )) deallocate(cradscr%lw_v_array       )
      if (associated(cradscr%radprof_array         )) deallocate(cradscr%radprof_array    )
      return
   end subroutine dealloc_radscratch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine nullifies the pointers, for a safe allocation.                   !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_radscratch(cradscr)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(radscrtype), target :: cradscr
      !------------------------------------------------------------------------------------!

      nullify(cradscr%pft_array                 )
      nullify(cradscr%leaf_temp_array           )
      nullify(cradscr%wood_temp_array           )
      nullify(cradscr%lai_array                 )
      nullify(cradscr%wai_array                 )
      nullify(cradscr%CA_array                  )
      nullify(cradscr%htop_array                )
      nullify(cradscr%hbot_array                )
      nullify(cradscr%par_level_beam            )
      nullify(cradscr%par_level_diffu           )
      nullify(cradscr%par_level_diffd           )
      nullify(cradscr%light_level_array         )
      nullify(cradscr%light_beam_level_array    )
      nullify(cradscr%light_diff_level_array    )
      nullify(cradscr%par_v_beam_array          )
      nullify(cradscr%rshort_v_beam_array       )
      nullify(cradscr%par_v_diffuse_array       )
      nullify(cradscr%rshort_v_diffuse_array    )
      nullify(cradscr%lw_v_array                )
      nullify(cradscr%radprof_array             )
      return
   end subroutine nullify_radscratch
   !=======================================================================================!
   !=======================================================================================!
end module canopy_radiation_coms
!==========================================================================================!
!==========================================================================================!

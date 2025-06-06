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
   ! ICANRAD -- Specifies how vertical canopy radiation is solved.  This variable sets     !
   !            both shortwave and longwave.                                               !
   !            0.  Two-stream model (Medvigy 2006), with the possibility to apply         !
   !                finite crown area to direct shortwave radiation.                       !
   !            1.  Multiple-scattering model (Zhao and Qualls 2005,2006), with the        !
   !                possibility to apply finite crown area to all radiation fluxes.        !
   !---------------------------------------------------------------------------------------!
   integer :: icanrad
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! IHRZRAD      -- Specifies how horizontal canopy radiation is solved.                  !
   !                 0.  Default ED-2.0: no horizontal patch shading.  All patches receive !
   !                     the same amount of light at the top.                              !
   !                 1.  A realized map of the plant community is built by randomly        !
   !                     assigning gaps associated with gaps (number of gaps proportional  !
   !                     to the patch area), and populating them with individuals,         !
   !                     respecting the cohort distribution in each patch.  The crown      !
   !                     closure index is calculated for the entire landscape and used     !
   !                     to change the amount of direct light reaching the top of the      !
   !                     canopy.  Patches are then split into 1-3 patches based on the     !
   !                     light condition (so expect simulations to be slower).  This       !
   !                     method is under development, suggestions on how to improve are    !
   !                     welcome.                                                          !
   !                 2.  Similar to option 1, except that height for trees with DBH >      !
   !                     DBH_crit are rescaled to calculate CCI.                           !
   !                 3.  Dummy horizontal canopy radiation.  This applies the same method  !
   !                     as 1 and 2 to split patches, but it does not change radiation     !
   !                     reaching the top of the canopy.  This is only useful to isolate   !
   !                     the effect of heterogeneous illumination from the patch count.    !
   !                 4.  Same as 0., but patch fusion takes into account the correction    !
   !                     for emergent trees.                                               !
   !---------------------------------------------------------------------------------------!
   integer :: ihrzrad
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
   ! real(kind=8), dimension(n_pft) :: leaf_scatter_tir
   ! real(kind=8), dimension(n_pft) :: wood_scatter_tir
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




   !---------------------------------------------------------------------------------------!
   !     The following variables control the method that allow light redistribution based  !
   ! on patch neighbourhood.  These are initialised in ed_xml_config.f90 or ed_params.f90. !
   !---------------------------------------------------------------------------------------!
   real(kind=4)    :: cci_radius   ! Maximum radius to calculate CCI               [     m]
   real(kind=4)    :: cci_pixres   ! Pixel resolution for TCH and CCI              [     m]
   real(kind=4)    :: cci_gapsize  ! Gap size                                      [     m]
   real(kind=4)    :: cci_gapmin   ! # of gaps associated with the smallest area   [   ---]
   integer         :: cci_nretn    ! "Return density" to generate the TCH map      [  1/m2]
   real(kind=4)    :: cci_hmax     ! Maximum height allowed in the CCI scheme      [     m]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These variables are derived from the properties above, they will be allocated     !
   ! during the initialisation step.                                                       !
   !---------------------------------------------------------------------------------------!
   !----- Total area of each single gap. --------------------------------------------------!
   real(kind=4)                              :: cci_gaparea
   !----- Number of grid points in x and y direction (pseudo-landscape). ------------------!
   integer                                   :: rls_nxy
   !----- Number of pixels in the pseudo-landscape. ---------------------------------------!
   integer                                   :: rls_npixel
   !----- Number of gaps in the pseudo-landscape. -----------------------------------------!
   integer                                   :: rls_ngap
   !----- 'raster' length along x/y axes. -------------------------------------------------!
   real(kind=4)                              :: rls_length
   !----- Number of pixels in each gap. ---------------------------------------------------!
   integer                                   :: gap_npixel
   !----- Total 'raster' landscape area. --------------------------------------------------!
   real(kind=4)                              :: rls_area
   !----- Use fixed thresholds to split patches by illumination classes? ------------------!
   logical                                   :: fixed_hrz_classes
   !----- Default thresholds in case fixed classes are to be used. ------------------------!
   real(kind=4)                              :: at_bright_def
   real(kind=4)                              :: at_dark_def
   !----- x of the 'raster' landscape. ----------------------------------------------------!
   real(kind=4), dimension(:,:), allocatable :: rls_x
   !----- y of the 'raster' landscape. ----------------------------------------------------!
   real(kind=4), dimension(:,:), allocatable :: rls_y
   !----- Top-of-canopy height. -----------------------------------------------------------!
   real(kind=4), dimension(:,:), allocatable :: rls_ztch
   !----- Crown closure index. ------------------------------------------------------------!
   real(kind=4), dimension(:,:), allocatable :: rls_cci
   !----- Absorption correction for incident beam radiation. ------------------------------!
   real(kind=4), dimension(:,:), allocatable :: rls_fbeam
   !----- Gap indices (zero is the default). ----------------------------------------------!
   integer     , dimension(:,:), allocatable :: rls_igp0
   integer     , dimension(:,:), allocatable :: rls_igp
   integer     , dimension(:,:), allocatable :: rls_ipa
   !----- Mask to decide which gaps can be used for any patch. ----------------------------!
   logical     , dimension(:,:), allocatable :: rls_mask
   !----- Gap origin. ---------------------------------------------------------------------!
   real(kind=4), dimension(:)  , allocatable :: gap_x0
   real(kind=4), dimension(:)  , allocatable :: gap_y0
   !----- Mean absorption correction for incident beam radiation. -------------------------!
   real(kind=4), dimension(:)  , allocatable :: gap_fbeam
   integer     , dimension(:)  , allocatable :: gap_nuse
   !----- Patch associated with the gap. --------------------------------------------------!
   integer     , dimension(:)  , allocatable :: gap_ipa
   !----- Auxiliary variable, patch index before shuffling, gap index after shuffling. ----!
   integer     , dimension(:)  , allocatable :: gap_idx
   !----- Mask to decide which gaps can be used for any patch. ----------------------------!
   logical     , dimension(:)  , allocatable :: gap_mask
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Hold these parameters as constants, the functional form may change soon.          !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: at0
   real(kind=4) :: at1
   real(kind=8) :: at08
   real(kind=8) :: at18
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

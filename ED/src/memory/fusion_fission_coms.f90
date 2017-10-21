!==========================================================================================!
!==========================================================================================!
! Module fusion_fission_coms.                                                              !
!     This modules contains some parameters to drive the patch and cohort merging and      !
! splitting, based on several criteria. These variables are tuners for such criteria.      !
!                                                                                          !
!     Unless the variables should be constants, you should never assign the initial value  !
! here. Instead, declare them here and add them in initialize_ff_coms to set up their      !
! default, initial values.                                                                 !
!==========================================================================================!
!==========================================================================================!
module fusion_fission_coms
   use ed_max_dims, only : str_len      ! ! intent(in)

   implicit none


   !---------------------------------------------------------------------------------------!
   !> IFUSION is a temporary flag for patch fusion. 
   !> 0 -- Original (ED-2.1) patch/cohort fusion routines.
   !> 1 -- Updated (ED-2.2) patch/cohort fusion routines.
   !> (MLO) The old scheme had some issues, mostly the fact that it used relative 
   !> differences in light levels, and had high rate of tolerance increase.  This would 
   !> eventually fuse patches with very different upper canopy.  Also, the patch fusion
   !> scheme did not check for number of remaining patches that were so small that would
   !> be terminated.  Both problems are much more likely to make a difference when initial
   !> conditions have a large number of patches (> 1000), which is quite common when using
   !> airborne lidar.  I am keeping the old scheme to avoid disrupting people's work, but
   !> eventually I would prefer deleting the routines. 
   !---------------------------------------------------------------------------------------!
   integer :: ifusion

   !----- Maximum number of iterations for patch fusion. ----------------------------------!
   integer :: niter_patfus

   !----- Minimum and maximum height (m) used in patch profiling --------------------------!
   real, dimension(:), allocatable :: hgt_class

   !----- Number of height bins in patch profiling ----------------------------------------!
   integer :: ff_nhgt
   !---------------------------------------------------------------------------------------!


   !---- Old patch fusion variables (slated to be deleted in the near future). ------------!
   real :: dark_cumlai_min
   real :: dark_cumlai_max
   real :: sunny_cumlai_min
   real :: sunny_cumlai_max
   real :: dark_cumlai_mult
   real :: sunny_cumlai_mult
   real :: light_toler_min
   real :: light_toler_max
   real :: light_toler_mult
   !---------------------------------------------------------------------------------------!


   !---- Old cohort fusion variables (slated to be deleted in the near future). -----------!
   real    :: fusetol           ! Cohort fusion tolerance on DBH (dimensionless) 
   real    :: fusetol_h         ! Cohort fusion tolerance on height (m) !
   real    :: lai_fuse_tol      ! Cohort fusion tolerance on LAI (m2 leaf/m2 ground)
   real    :: coh_tolerance_max ! Cohort maximum tolerance factor 
   logical :: fuse_relax        ! Flag to allow a less strict fusion test
   !---------------------------------------------------------------------------------------!


   !---- New patch fusion variables. ------------------------------------------------------!
   real :: pat_light_ext        ! Extinction coefficient for patch fusion.  This is more
                                !    like ED-1.0, but for simplicity we compare patch 
                                !    similarity using Beer's law. 
   real :: pat_light_tol_min    ! Minimum tolerance for patch light difference.
   real :: pat_light_tol_max    ! Maximum tolerance for patch light difference.
   real :: pat_light_tol_mult   ! Multiplier for the light tolerance.
   real :: pat_light_mxd_fac    ! Light tolerance for maximum deviation.
   real :: pat_diff_age_tol     ! Maximum age difference to be considered same age [yr].
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Minimum area to remain resolved.  This condition is normally met, except when    !
   ! initialising the simulation with massive amount of data (like airborne lidar data).   !
   !---------------------------------------------------------------------------------------!
   real    :: pat_min_area_remain
   !---------------------------------------------------------------------------------------!

   !---- Cohort fusion variables. ---------------------------------------------------------!
   integer :: niter_cohfus      ! Number of cohort fusion iterations.
   real    :: coh_size_tol_min  ! Minimum tolerance for relative size difference.
   real    :: coh_size_tol_max  ! Maximum tolerance for relative size difference.
   real    :: coh_size_tol_mult ! Multiplier for the relative size tolerance.
   !---------------------------------------------------------------------------------------!



   !----- Cohort splitting tolerance on LAI (m2 leaf/ m2 ground) --------------------------!
   real    :: lai_tol
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Correlation coefficient that is assumed between two patches and two cohorts when !
   ! they are fused.  This only affects the mean sum of squares.                           !
   !---------------------------------------------------------------------------------------!
   real    :: corr_patch
   real    :: corr_cohort
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Flag to decide whether we should print the full details of the patch fusion       !
   ! process.                                                                              !
   !---------------------------------------------------------------------------------------!
   logical                :: print_fuse_details
   character(len=str_len) :: fuse_prefix
   !---------------------------------------------------------------------------------------!

end Module fusion_fission_coms
!==========================================================================================!
!==========================================================================================!

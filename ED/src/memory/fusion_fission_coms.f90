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
   use ed_max_dims, only : str_len      & ! intent(in)
                         , n_dist_types ! ! intent(in)

   implicit none

   !----- Maximum number of iterations for patch fusion. ----------------------------------!
   integer :: niter_patfus

   !----- Minimum and maximum height (m) used in patch profiling --------------------------!
   real, dimension(:), allocatable :: hgt_class

   !----- Number of height bins in patch profiling ----------------------------------------!
   integer :: ff_nhgt

   !----- Cohort fusion tolerance on DBH (dimensionless) ----------------------------------!
   real    :: fusetol

   !----- Cohort fusion tolerance on height (m) -------------------------------------------!
   real    :: fusetol_h

   !----- Cohort fusion tolerance on LAI (m2 leaf/m2 ground) ------------------------------!
   real    :: lai_fuse_tol

   !----- Cohort splitting tolerance on LAI (m2 leaf/ m2 ground) --------------------------!
   real    :: lai_tol

   !----- Cohort maximum tolerance factor -------------------------------------------------!
   real    :: coh_tolerance_max

   !---- Patch fusion variables. ----------------------------------------------------------!
   real :: pat_light_ext        ! Extinction coefficient for patch fusion.  This is more
                                !    like ED-1.0, but for simplicity we compare patch 
                                !    similarity using Beer's law. 
   real :: pat_light_tol_min    ! Minimum tolerance for patch light difference
   real :: pat_light_tol_max    ! Maximum tolerance for patch light difference
   real :: pat_light_tol_mult   ! Multiplier for the light tolerance
   real :: pat_light_mxd_fac    ! Light tolerance for maximum deviation.
   real :: pat_diff_age_tol     ! Maximum age difference to be considered same age [yr]
   !---------------------------------------------------------------------------------------!

   !----- Flag to allow a less strict fusion test for short cohorts. ----------------------!
   logical :: fuse_relax

   !---------------------------------------------------------------------------------------!
   !      Correlation coefficient that is assumed between two patches and two cohorts when !
   ! they are fused.  This only affects the mean sum of squares.                           !
   !---------------------------------------------------------------------------------------!
   real    :: corr_patch
   real    :: corr_cohort
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Minimum age above which we disregard the disturbance type (land use) and assume  !
   ! old growth, thus allowing patch fusion to occur.                                      !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_dist_types) :: min_oldgrowth
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

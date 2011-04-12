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
   use ed_max_dims, only : str_len

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

   !----- Patch maximum tolerance factor --------------------------------------------------!
   real    :: dark_cumlai_min
   real    :: dark_cumlai_max
   real    :: sunny_cumlai_min
   real    :: sunny_cumlai_max
   real    :: dark_cumlai_mult
   real    :: sunny_cumlai_mult

   !----- Flag to allow a less strict fusion test for short cohorts. ----------------------!
   logical :: fuse_relax

   !----- Fractional tolerance for patch pft height comparisons (dimensionless) -----------!
   real    :: light_toler_min
   real    :: light_toler_max
   real    :: light_toler_mult

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

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
  
   implicit none

   !----- Minimum DBH class used in patch profiling ---------------------------------------!
   real :: min_dbh_class

   !----- Maximum DBH (cm) used in patch profiling ----------------------------------------!
   real :: maxdbh

   !----- Number of DBH bins in patch profiling -------------------------------------------!
   integer :: ff_ndbh

   !----- Minimum height class in patch profiling -----------------------------------------!
   real :: min_hgt_class

   !----- Cohort fusion tolerance on DBH (dimensionless) ----------------------------------!
   real :: fusetol

   !----- Cohort fusion tolerance on height (m) -------------------------------------------!
   real :: fusetol_h

   !----- Cohort fusion tolerance on LAI (m2 leaf/m2 ground) ------------------------------!
   real :: lai_fuse_tol

   !----- Cohort splitting tolerance on LAI (m2 leaf/ m2 ground) --------------------------!
   real :: lai_tol

   !----- Cohort maximum tolerance factor -------------------------------------------------!
   real :: coh_tolerance_max

   !----- Patch maximum tolerance factor --------------------------------------------------!
   real :: pat_tolerance_max

   !----- Flag to allow a less strict fusion test for short cohorts. ----------------------!
   logical :: fuse_relax

   !---------------------------------------------------------------------------------------!
   !    Minimum plant density for height bin to be used in height profile comparisons      !
   ! (plants/m2).                                                                          !
   !---------------------------------------------------------------------------------------!
   real :: ntol

   !----- Fractional tolerance for patch pft height comparisons (dimensionless) -----------!
   real :: profile_tol

   !----- Maximum patch age for fusion (years) --------------------------------------------!
   real :: max_patch_age

end Module fusion_fission_coms
!==========================================================================================!
!==========================================================================================!

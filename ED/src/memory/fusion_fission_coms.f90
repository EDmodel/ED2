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

   !----- Minimum DBH class used in patch profiling ---------------------------------------!
   real    :: min_dbh_class

   !----- Maximum DBH (cm) used in patch profiling ----------------------------------------!
   real    :: maxffdbh

   !----- Maximum height (m) used in patch profiling --------------------------------------!
   real    :: maxffhgt

   !----- Number of DBH bins in patch profiling -------------------------------------------!
   integer :: ff_ndbh

   !----- Number of height bins in patch profiling ----------------------------------------!
   integer :: ff_nhgt

   !----- Inverse of DBH bin class in patch profiling -------------------------------------!
   real    :: dffdbhi

   !----- Inverse of height bin class in patch profiling ----------------------------------!
   real    :: dffhgti

   !----- Minimum height class in patch profiling -----------------------------------------!
   real    :: min_hgt_class

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
   real    :: pat_tolerance_max

   !----- Flag to allow a less strict fusion test for short cohorts. ----------------------!
   logical :: fuse_relax

   !---------------------------------------------------------------------------------------!
   !    Minimum plant density for height bin to be used in height profile comparisons      !
   ! (plants/m2).                                                                          !
   !---------------------------------------------------------------------------------------!
   real    :: ntol

   !---------------------------------------------------------------------------------------!
   !    Minimum leaf area index for height bin to be used in height profile comparisons    !
   ! (plants/m2).                                                                          !
   !---------------------------------------------------------------------------------------!
   real    :: laimax_tol

   !----- Fractional tolerance for patch pft height comparisons (dimensionless) -----------!
   real    :: profile_tol

   !----- Maximum patch age for fusion (years) --------------------------------------------!
   real    :: max_patch_age
 
   !---------------------------------------------------------------------------------------!
   !     Maximum cumulative (potential) LAI.  Cohorts beneath the layer in which LAI       !
   ! reaches this value are ignored for patch fusion, because they are too dark anyway.    !
   !---------------------------------------------------------------------------------------!
   real    :: dark_cumlai
  
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

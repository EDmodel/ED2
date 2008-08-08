Module fusion_fission_coms

! See initialize_ff_coms for default, initial values

! Minimum biomass density (kgC/m2) required to form a new recruit.
real :: min_recruit_size

! Minimum DBH class used in patch profiling
real :: min_dbh_class

! Maximum DBH (cm) used in patch profiling
real :: maxdbh

! Number of DBH bins in patch profiling
integer, parameter :: ff_ndbh = 20

! Minimum height class in patch profiling
real :: min_hgt_class

! Cohort fusion tolerance on DBH (dimensionless)
real :: fusetol

! Cohort fusion tolerance on height (m)
real :: fusetol_h

! Cohort fusion tolerance on LAI (m2 leaf/m2 ground)
real :: lai_fuse_tol

! Cohort splitting tolerance on LAI (m2 leaf/ m2 ground)
real :: lai_tol

! Min plant density for height bin to be used in height profile comparisons (plants/m2)
real :: ntol

! Fractional tolerance for patch pft height comparisons (dimensionless)
real :: profile_tol

! Maximum patch age for fusion (years)
real :: max_patch_age

end Module fusion_fission_coms

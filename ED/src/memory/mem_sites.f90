Module mem_sites

  use ed_max_dims, only: max_soi, max_ed_regions

  implicit none

  integer :: n_soi
  real, dimension(max_soi) :: soi_lat
  real, dimension(max_soi) :: soi_lon
  integer :: n_ed_region
  real, dimension(max_ed_regions) :: ed_reg_latmin
  real, dimension(max_ed_regions) :: ed_reg_latmax
  real, dimension(max_ed_regions) :: ed_reg_lonmin
  real, dimension(max_ed_regions) :: ed_reg_lonmax
  real :: grid_res
  integer :: grid_type

  ! Restart resolution in case of ascii restart file
  real :: edres

  
  integer :: maxpatch   ! Benchmark maximum number of patches per site
  integer :: maxcohort  ! Benchmark maximum number of cohorts per patch

end Module mem_sites

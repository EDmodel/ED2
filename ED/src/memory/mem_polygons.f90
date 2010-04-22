!==========================================================================================!
!==========================================================================================!
!     This module contains some polygon-level global variables.                            !
!------------------------------------------------------------------------------------------!
module mem_polygons

   use ed_max_dims, only : max_poi         & ! intent(in)
                         , max_ed_regions  ! ! intent(in)

   implicit none

   !----- POI-type run variables. ---------------------------------------------------------!
   integer                         :: n_poi
   real, dimension(max_poi)        :: poi_lat
   real, dimension(max_poi)        :: poi_lon

   !----- Regional-type run variables. ----------------------------------------------------!
   integer                         :: n_ed_region
   real, dimension(max_ed_regions) :: ed_reg_latmin
   real, dimension(max_ed_regions) :: ed_reg_latmax
   real, dimension(max_ed_regions) :: ed_reg_lonmin
   real, dimension(max_ed_regions) :: ed_reg_lonmax

   !----- Regional-type grid structure variables. -----------------------------------------!
   real    :: grid_res
   integer :: grid_type

   !---------------------------------------------------------------------------------------!
   !    Restart resolution in case of ascii restart file (currently this is still used     !
   ! only by single-polygon, multiple site runs using restart, and could be phased out     !
   ! using the same method used in ed_read_ed10_ed20_history.f90.                          !
   !---------------------------------------------------------------------------------------!
   real :: edres

   !---------------------------------------------------------------------------------------!
   !    These variables are the sought maximum number of patches per site (MAXPATCH) and   !
   ! the maximum number of cohorts per patch (MAXCOHORT).  Notice that these numbers are   !
   ! not strict and the maximum number may actually exceed this when the model doesn't     !
   ! find any possibility of fusing patches or cohorts.                                    !
   !---------------------------------------------------------------------------------------!
   integer :: maxpatch
   integer :: maxcohort

end module mem_polygons
!==========================================================================================!
!==========================================================================================!

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
   !    Single polygon resolution in degrees.  This is only used to define the soil types  !
   ! for a multiple-site run.                                                              !
   !---------------------------------------------------------------------------------------!
   real, dimension(max_poi)        :: poi_res

   !---------------------------------------------------------------------------------------!
   !    Restart resolution in case of ascii restart file (currently this is still used     !
   ! only by single-polygon, multiple site runs using restart, and could be phased out     !
   ! using the same method used in ed_read_ed10_ed20_history.f90.                          !
   !---------------------------------------------------------------------------------------!
   real :: edres


   !---------------------------------------------------------------------------------------!
   !    This is the maximum allowed number of sites for a simulation.  This is a strict    !
   ! maximum, but the actual number of sites can be less in case there aren't enough       !
   ! different types of soil.                                                              !
   !---------------------------------------------------------------------------------------!
   integer :: maxsite

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

module grid_coms

  ! This module contains the dimensions involving ED and work structures. This is based
  ! on RAMS mem_grid.f90 module, with unused variables taken away. This does not contain
  ! a grid_g structure either, the only two variable we would need (glon/glat) were 
  ! merged with the work structure (ed_work_vars.f90)

  use ed_max_dims, only : nxpmax,nypmax,maxgrds
  

  real(kind=8)       :: time                ! Simulation time since the beginning
  real(kind=8)       :: timmax              ! Total time

  integer :: ngrids                         ! how many grids (by default in ED)
                                            ! this is n_soi+n_grids
  integer :: nzg                            ! Number of soil levels
  integer :: nzs                            ! Number of snow/surface water levels

  integer :: nxtnest(maxgrds)               ! next coarser grid number 
                                            ! (0 if grid is not nested)

  integer, target    :: nnxp(maxgrds)       ! grid cells at x direction
  integer            :: nstratx(maxgrds)    ! nest ratio for next coarser grid
  real               :: deltaxn(maxgrds)    ! delta x
  real               :: xtn(nxpmax,maxgrds) ! x coordinate of cell center on polar stereographic projection
  real               :: xmn(nxpmax,maxgrds) ! x coordinate of higher cell boundary on polar stereographic projection
  integer            :: ninest(maxgrds)     ! index on next coarser grid where this grid starts (lower southwest corner)
  !                                         ! ds to interpolate x direction from coarser to finner grids: 

  integer :: istp
  integer, target    :: nnyp(maxgrds)       ! grid cells at y direction
  integer            :: nstraty(maxgrds)    ! nest ratio for next coarser grid
  real               :: deltayn(maxgrds)    ! delta y
  real               :: ytn(nypmax,maxgrds) ! y coordinate of cell center on polar stereographic projection
  real               :: ymn(nypmax,maxgrds) ! y coordinate of higher cell boundary on polar stereographic projection
  integer            :: njnest(maxgrds)     ! index on next coarser grid where this grid starts (lower southwest corner)
  !                                         ! ds to interpolate y direction from coarser to finner grids: 

  real               :: polelat             ! Pole longitude read in the namelist
  real               :: polelon             ! Pole longitude read in the namelist
  real               :: platn(maxgrds)      ! pole latitude (degrees)
  real               :: plonn(maxgrds)      ! pole longitude (degrees)

  real               :: centlat(maxgrds)    ! grid center latitude (degrees)
  real               :: centlon(maxgrds)    ! grid center longitude (degrees)

  integer, parameter :: jdim=1              ! all horizontal grids are 1D (jdim=0) or 2D (jdim=1)

  integer            :: ngrid               ! Current grid
  integer            :: nxp                 ! Current grid number of points in x
  integer            :: nyp                 ! Current grid number of points in y
  real               :: deltax              ! Current grid resolution in x direction
  real               :: deltay              ! Current grid resolution in y direction
  real               :: xt(nxpmax)          ! Current grid xtn
  real               :: xm(nxpmax)          ! Current grid xmn
  real               :: yt(nypmax)          ! Current grid ytn
  real               :: ym(nypmax)          ! Current grid ymn


end module grid_coms

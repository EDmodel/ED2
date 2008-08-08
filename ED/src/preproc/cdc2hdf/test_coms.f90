Module test_coms

  implicit none

  type ncvar
     character(len=128) :: filename
     integer :: ncid
     integer :: varid
     real :: minlat
     real :: maxlat
     real :: latstep
     real :: minlon
     real :: maxlon
     integer :: nlat
     integer :: nlon
     real :: lonstep
     integer :: mintime
     integer :: maxtime
     integer :: N
     real :: offset
     real :: scale
     integer :: nodata
  end type ncvar

end Module test_coms

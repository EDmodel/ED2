program main

  use test_coms, only: ncvar

  implicit none

  character(len=256), parameter :: input_dir='input/'
  character(len=256) :: tmpvar_name
  type(ncvar) :: tmpvar

  tmpvar_name = trim(input_dir)//'uwnd.sig995.1948.nc'
  call new_ncvar(tmpvar_name, 'uwnd', tmpvar)

!  tmpvar = new_ncvar(tmpvar_name,"uwnd")

end program main

subroutine new_ncvar(filename, varname, var)
  use test_coms, only: ncvar
  use netcdf
  implicit none
  type(ncvar) :: var
  type(ncvar) :: v
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer :: nfstat

  ! open file
  nfstat = nf90_open(path=trim(filename), cmode=NF90_NOWRITE, ncid=v%ncid);
!  check_ncstat(ncstat, filename);


  return
end subroutine new_ncvar

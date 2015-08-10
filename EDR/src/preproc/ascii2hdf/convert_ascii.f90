! SUMMARY:
!    This program reads in an old-style ED2 met driver file for SOIs.  The
!       input data is in ASCII, the output data is in HDF5.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES YOU MAY NEED TO CHANGE:
!   first_year --- first year to process
!   last_year  --- last year to process
!   findir       --- directory of the input data
!   foutdir       --- directory of the output data
!   frqin ---  Frequency at which input (non-radiation) data was written 
!              out [seconds].
!   frqin_rad --- Frequency at which the input radiation data was written 
!                  out [seconds].
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
  implicit none

  ! USER CONTROL SECTION
  !------------------------------------------------------------------
  integer :: first_year=1997
  integer :: last_year=1997
  character(len=256) :: findir = 'input/'
  character(len=256) :: foutdir = 'output/'
  real :: frqin = 21600.0
  real :: frqin_rad = 900.0
  !------------------------------------------------------------------

  integer :: narg
  integer :: process_year
  character(len=256) :: tmp

  ! Optionally read arguements from command line (mcd)
  narg = IARGC ()
  if(narg .ge. 6) then
      CALL GETARG(1 , tmp)
      read(tmp,'(I4)') first_year
      CALL GETARG(2 , tmp)
      read(tmp,'(I4)') last_year
      CALL GETARG(3 , findir)
      CALL GETARG(4 , foutdir)
      CALL GETARG(5 , tmp)
      read(tmp,'(F20.3)') frqin
      CALL GETARG(6 , tmp)
      read(tmp,'(F20.3)') frqin_rad      
      print*,first_year,last_year,trim(findir),trim(foutdir),frqin,frqin_rad
  endif

  !process years
  do process_year = first_year,last_year
     call make_the_month(process_year, findir, foutdir, frqin, frqin_rad)
  enddo

end program main

!==============================================================
subroutine make_the_month(process_year, findir, foutdir, frqin, frqin_rad)

  use hdf5_utils

  implicit none
  
  real, intent(in) :: frqin
  real, intent(in) :: frqin_rad
  character(len=*), intent(in) :: findir
  character(len=*), intent(in) :: foutdir
  integer, intent(in) :: process_year
  integer :: times_per_day_nr
  integer :: times_per_day_r
  integer :: tsize_r
  integer :: tsize_nr
  integer :: m1
  integer :: m2
  integer :: process_month
  real, allocatable, dimension(:,:,:) :: hgt
  real, allocatable, dimension(:,:,:) :: uwnd
  real, allocatable, dimension(:,:,:) :: vwnd
  real, allocatable, dimension(:,:,:) :: air
  real, allocatable, dimension(:,:,:) :: shum
  real, allocatable, dimension(:,:,:) :: prss
  real, allocatable, dimension(:,:,:) :: prate
  real, allocatable, dimension(:,:,:) :: dlwrf
  real, allocatable, dimension(:,:,:) :: nbdsf
  real, allocatable, dimension(:,:,:) :: nddsf
  real, allocatable, dimension(:,:,:) :: vbdsf
  real, allocatable, dimension(:,:,:) :: vddsf
  real, allocatable, dimension(:,:,:) :: rshort
  real, allocatable, dimension(:,:,:) :: rshortd
  integer :: nlon
  integer :: nlat
  character(len=256) :: infile_nr
  character(len=256) :: infile_r
  character(len=256) :: outfile
  integer :: maxday
  integer :: looplat
  integer :: looplon
  integer :: itime
  real, dimension(11) :: var_in
  real, parameter :: p00=1.e5
  real, parameter :: cp = 1004.
  real, parameter :: rdry = 287.
  real, parameter :: cpi = 1. / cp
  real, parameter :: cpor = cp / rdry
  integer :: ndims_hdf
  integer, dimension(3) :: idims_hdf
  character(len=3), dimension(12) :: mnames = (/'JAN','FEB','MAR','APR',  &
       'MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/)

  !========================================================================

  times_per_day_nr = nint(86400.0/frqin)
  tsize_nr = 31 * times_per_day_nr
  times_per_day_r = nint(86400.0/frqin_rad)
  tsize_r = 31 * times_per_day_r

  m1 = 1
  m2 = 12

  nlon = 1
  nlat = 1

  ! Loop over months
  do process_month = m1,m2

     ! Allocate arrays for the first month
     if(process_month == m1)then
        allocate(hgt(tsize_nr,nlon,nlat))
        allocate(uwnd(tsize_nr,nlon,nlat))
        allocate(vwnd(tsize_nr,nlon,nlat))
        allocate(air(tsize_nr,nlon,nlat))
        allocate(shum(tsize_nr,nlon,nlat))
        allocate(prss(tsize_nr,nlon,nlat))
        allocate(prate(tsize_nr,nlon,nlat))
        allocate(dlwrf(tsize_nr,nlon,nlat))
        allocate(rshort(tsize_r,nlon,nlat))
        allocate(rshortd(tsize_r,nlon,nlat))
        allocate(nbdsf(tsize_r,nlon,nlat))
        allocate(nddsf(tsize_r,nlon,nlat))
        allocate(vbdsf(tsize_r,nlon,nlat))
        allocate(vddsf(tsize_r,nlon,nlat))
     endif     

     ! Open input data file
     print*,'trying time: ',process_month,process_year
     write(infile_nr,'(a,i2.2,i4.4,a)')trim(findir), process_month,  &
          process_year, '.dat'
     open(12,file=trim(infile_nr),form='formatted',status='old')
     write(infile_r,'(a,i2.2,i4.4,a)')trim(findir), process_month,  &
          process_year, '-rad.dat'
     open(13,file=trim(infile_r),form='formatted',status='old')

     ! Figure out how many days in this month
     maxday = 31
     if(process_month .eq.4     .or.    &
          process_month .eq.6     .or.  &
          process_month .eq.9     .or.  &
          process_month .eq.11)then
        maxday = 30
     elseif(process_month.eq.2)then
        if(0.25*process_year-int(0.25*process_year).lt.0.1)then
           maxday = 29
        else
           maxday = 28
        endif
     endif

     ! Read in data
     do looplat=1,nlat
        do looplon=1,nlon
           do itime = 1,times_per_day_nr * maxday
              read(12,*)var_in(1:11)
              hgt(itime,looplon,looplat)  = max(10.0,var_in(1)) ! Height
              uwnd(itime,looplon,looplat) = var_in(2)           ! Zonal wind
              vwnd(itime,looplon,looplat) = var_in(3)           ! Meridional wind

              ! MLO - Question, is the input variable really specific humidity or it is mixing ratio?
              shum(itime,looplon,looplat) = max(0.0,var_in(5))  ! Specific humidity

              !KIM - var_in(4): potential temperature
              !      var_in(6): Exner function
              !            air: temperature
              !           prss: Atmospheric pressure
              air(itime,looplon,looplat)   = var_in(4) * var_in(6) * cpi    ! Air temperature
              prss(itime,looplon,looplat)  = p00 * (var_in(6) * cpi)**cpor  ! Pressure
              prate(itime,looplon,looplat) = max(0.0,var_in(8)+var_in(9))   ! Modified to include pcpg in prate sum (mcd)
              dlwrf(itime,looplon,looplat) = max(0.0,var_in(11))            ! Downwelling longwave radiation flux
           enddo
           do itime = 1,times_per_day_r * maxday
              read(13,*)var_in(1:2)
              rshort(itime,looplon,looplat) = max(0.0,var_in(1))
              rshortd(itime,looplon,looplat) = max(0.0,var_in(2))
           enddo
        enddo
     enddo

     close(12)
     close(13)

     ! Open output file for precip and long wave
     write(outfile,'(a,i4.4,a)')trim(foutdir)//'ED_OL_', process_year,  &
          mnames(process_month)//'.h5'
     call shdf5_open_f(trim(outfile),'W',1)

     ! Write out the non-radiation data.
     ndims_hdf = 3
     idims_hdf(1) = times_per_day_nr * maxday
     idims_hdf(2) = nlon
     idims_hdf(3) = nlat

     call shdf5_orec_f(ndims_hdf, idims_hdf, 'hgt', rvara=hgt)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'tmp', rvara=air)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'pres', rvara=prss)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'sh', rvara=shum)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'ugrd', rvara=uwnd)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'vgrd', rvara=vwnd)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'prate', rvara=prate)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'dlwrf', rvara=dlwrf)
     
     ! Write out the radiation data
     idims_hdf(1) = times_per_day_r * maxday
     idims_hdf(2) = nlon
     idims_hdf(3) = nlat
     nbdsf = (rshort - rshortd) * 0.57
     nddsf = rshortd * 0.48
     vbdsf = (rshort - rshortd) * 0.43
     vddsf = rshortd * 0.52
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'nbdsf', rvara=nbdsf)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'nddsf', rvara=nddsf)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'vbdsf', rvara=vbdsf)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'vddsf', rvara=vddsf)

     ! Close the output file.
     call shdf5_close_f()

     if(process_month == m2)then

        deallocate(hgt)
        deallocate(uwnd)
        deallocate(vwnd)
        deallocate(air)
        deallocate(shum)
        deallocate(prss)
        deallocate(prate)
        deallocate(dlwrf)
        deallocate(rshort)
        deallocate(rshortd)
        deallocate(nbdsf)
        deallocate(nddsf)
        deallocate(vbdsf)
        deallocate(vddsf)

     endif

  enddo

  return

end subroutine make_the_month

! SUMMARY:
!   1)  This program reads in the following ASCII fields:
!      -- precipitation 
!      -- downward long wave radiation at the surface
!      -- NIR beam
!      -- NIR diffuse
!      -- PAR beam
!      -- PAR diffuse.
!   2)  The precipitation and long wave are then written out in HDF5 format.  
!   3)  The solar radiation fields are interpolated from longer (typically 
!   6-hourly resolution) to shorter time steps (typically about 15-60 minutes).
!   The interpolation conserves the total integral of radiation and 
!   accounting for variation in the solar zenith angle.  
!   4)  The interpolated solar radiation fields are then written out in HDF5.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES YOU MAY NEED TO CHANGE:
!   first_year --- first year to process
!   last_year  --- last year to process
!   fdir       --- directory of the input data (output files are placed in the
!                  same directory.
!   frqin ---  Frequency at which input data was written out [seconds].
!   frqout_rad --- Frequency at which you want the interpolated output 
!                  solar radiation data written [seconds].
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
  implicit none

  ! USER CONTROL SECTION
  !------------------------
  integer, parameter :: first_year=1948
  integer, parameter :: last_year=2007
  character(len=256), parameter :: fdir = 'amazon_met_driver/'
  real, parameter :: frqin = 21600.0 !  Frequency at which input data was 
                                     !  written out.
  real, parameter :: frqout_rad = 3600.0 !  Frequency at which you want the 
                                         !  output radiation written.
  !------------------------

  integer :: process_year

  do process_year = first_year,last_year
     call make_the_month(process_year, fdir, frqin, frqout_rad)
  enddo

end program main

!==============================================================
subroutine make_the_month(process_year, fdir, frqin, frqout_rad)

  use hdf5_utils

  implicit none
  
  real, intent(in) :: frqin
  real, intent(in) :: frqout_rad
  character(len=*), intent(in) :: fdir
  integer, intent(in) :: process_year
  integer :: m1
  integer :: m2
  integer :: process_month
  integer :: ndays
  character(len=256) :: infile
  character(len=256) :: radfile
  character(len=256) :: sfcfile
  integer :: idum1
  integer :: idum2
  integer :: nlat
  integer :: nlon
  real :: deltay
  real :: deltax
  real :: lat_min
  real :: lon_min
  integer :: maxday
  integer :: looplat
  integer :: looplon
  integer :: itime
  real :: lat
  real :: lon
  real :: time
  integer :: idate
  integer :: offset_from_utc
  integer :: date_utc
  integer :: id
  integer :: ib
  integer :: ih
  integer :: ntime
  integer, external :: julday
  real, external :: zen
  real, external :: dirfrac
  real :: time_utc
  character(len=3), dimension(12) :: mnames = (/'JAN','FEB','MAR','APR',  &
       'MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/)
  character(len=256) :: outfile
  integer :: ndims_hdf
  integer, dimension(3) :: idims_hdf
  integer :: tsize_input
  integer :: tsize_output
  integer :: times_per_day_in
  real :: time_inc_in
  integer :: tsize_out_per_in
  integer :: iout
  integer :: tindex
  real :: rad_norm

  real, allocatable, dimension(:) :: var_in
  real, allocatable, dimension(:,:,:) :: air_in
  real, allocatable, dimension(:,:,:) :: pres_in
  real, allocatable, dimension(:,:,:) :: shum_in
  real, allocatable, dimension(:,:,:) :: uwnd_in
  real, allocatable, dimension(:,:,:) :: hgt_in
  real, allocatable, dimension(:,:,:) :: vwnd_in
  real, allocatable, dimension(:,:,:) :: prate_in
  real, allocatable, dimension(:,:,:) :: dlwrf_in
  real, allocatable, dimension(:,:,:) :: nbdsf_in
  real, allocatable, dimension(:,:,:) :: nddsf_in
  real, allocatable, dimension(:,:,:) :: vbdsf_in
  real, allocatable, dimension(:,:,:) :: vddsf_in
  real, allocatable, dimension(:,:,:) :: nbdsf_output
  real, allocatable, dimension(:,:,:) :: nddsf_output
  real, allocatable, dimension(:,:,:) :: vbdsf_output
  real, allocatable, dimension(:,:,:) :: vddsf_output
  real, allocatable, dimension(:) :: zen_norm
  logical, external :: isleap
  !========================================================================

  times_per_day_in = nint(86400.0/frqin)
  allocate(var_in(11))
  tsize_out_per_in = nint(frqin / frqout_rad)
  allocate(zen_norm(tsize_out_per_in))


  m1 = 1
  m2 = 12

  ! Loop over months
  monthloop: do process_month = m1,m2
     


     ! Figure out how many days in this month
     select case (maxday)
     case (1,3,5,7,8,10,12)
        maxday = 31
     case (4,6,9,11)
        maxday = 30

     case(2)
        if(isleap(process_year))then
           maxday = 29
        else
           maxday = 28
        endif
     end select

     ! Maximum number of input times in a single month
     tsize_input = maxday * times_per_day_in
   
     ! Maximum number of output times in a single month
     tsize_output = maxday * nint(86400.0/frqout_rad)





     ! Open and read CDC header file
     infile = trim(fdir)//'SOUTHAM_HEADER'
     open(12,file=trim(infile),form='formatted',status='old')
     read(12,*)nlat,nlon
     read(12,*)deltay,deltax,lat_min,lon_min
     close(12)

     ! This is slower than the other previous style, but it should be safer. The variable will have
     ! exactly the same size it needs to have.
     allocate(air_in(tsize_input,nlon,nlat))
     allocate(pres_in(tsize_input,nlon,nlat))
     allocate(shum_in(tsize_input,nlon,nlat))
     allocate(uwnd_in(tsize_input,nlon,nlat))

     allocate(vwnd_in(tsize_input,nlon,nlat))
     allocate(prate_in(tsize_input,nlon,nlat))
     allocate(hgt_in(tsize_input,nlon,nlat))
     hgt_in = 10.0
     allocate(dlwrf_in(tsize_input,nlon,nlat))

     allocate(nbdsf_in(tsize_input,nlon,nlat))
     allocate(nddsf_in(tsize_input,nlon,nlat))
     allocate(vbdsf_in(tsize_input,nlon,nlat))
     allocate(vddsf_in(tsize_input,nlon,nlat))

     allocate(nbdsf_output(tsize_output,nlon,nlat))
     allocate(nddsf_output(tsize_output,nlon,nlat))
     allocate(vbdsf_output(tsize_output,nlon,nlat))
     allocate(vddsf_output(tsize_output,nlon,nlat))

     ! Open input data file (on the coarse time step)
     print*,'trying time: ',process_month,process_year
     write(infile,'(a,i4.4,a)')trim(fdir)//'SOUTHAM_', process_year,  &
          mnames(process_month)//'.dat'
     open(12,file=trim(infile),form='formatted',status='old')

     ! Read in data and convert units on the radiation
     do looplat=1,nlat
        do looplon=1,nlon
           do itime = 1,tsize_input
              read(12,*)var_in(1:11)
              prate_in(itime,looplon,looplat) = max(0.0,var_in(1))
              dlwrf_in(itime,looplon,looplat) = max(0.0,var_in(2))
              nbdsf_in(itime,looplon,looplat) = max(0.0,var_in(3))
              nddsf_in(itime,looplon,looplat) = max(0.0,var_in(4))
              vbdsf_in(itime,looplon,looplat) = max(0.0,var_in(5))
              vddsf_in(itime,looplon,looplat) = max(0.0,var_in(6))
              air_in(itime,looplon,looplat) = var_in(7)
              pres_in(itime,looplon,looplat) = var_in(8)
              shum_in(itime,looplon,looplat) = max(0.0,var_in(9))
              uwnd_in(itime,looplon,looplat) = var_in(10)
              vwnd_in(itime,looplon,looplat) = var_in(11)
           enddo
        enddo
     enddo

     close(12)

     ! Open output file for precip and long wave
     write(outfile,'(a,i4.4,a)')trim(fdir)//'SOUTHAM_OL1_', process_year,  &
          mnames(process_month)//'.h5'
     call shdf5_open_f(trim(outfile),'W',1)

     ! Write out the precip and longwave data.
     ndims_hdf = 3
     idims_hdf(1) = tsize_input
     idims_hdf(2) = nlon
     idims_hdf(3) = nlat

     call shdf5_orec_f(ndims_hdf, idims_hdf, 'hgt', rvara=hgt_in)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'tmp', rvara=air_in)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'pres', rvara=pres_in)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'sh', rvara=shum_in)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'ugrd', rvara=uwnd_in)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'vgrd', rvara=vwnd_in)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'prate', rvara=prate_in)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'dlwrf', rvara=dlwrf_in)
     
     ! Close the precip and long wave file.
     call shdf5_close_f()

     ! This code converts 6-hourly radiation to hourly radiation conserving
     ! total integral of radiation and accounting for variation in the 
     ! solar zenith angle.

     ! Loop over latitudes
     do looplat=1,nlat
        lat = lat_min + (nlat - looplat)*deltay

        ! Loop over longitudes
        do looplon=1,nlon
           lon = lon_min + (looplon-1)*deltax

           ! Loop over input data temporal bins
           time = 0.0  ! initial number of seconds elapsed in this month
           time_inc_in = 86400.0 / real(times_per_day_in) ! time inc [sec]
           do itime=1, tsize_input

              ! Compute the day of month and the number of seconds elapsed
              ! in this month.
              idate = int((itime-1)/times_per_day_in)+1

              ! Get normalization factors for this bin.

              call calculate_bin_normalization(process_year &
                   ,process_month &
                   ,idate,time,lat & 
                   ,lon, frqin, frqout_rad, tsize_out_per_in  &
                   ,rad_norm &
                   ,zen_norm)

              do iout = 1, tsize_out_per_in
                 tindex = (itime - 1) * tsize_out_per_in + iout
                 nbdsf_output(tindex,looplon,looplat) = rad_norm *   &
                      zen_norm(iout) * nbdsf_in(itime,looplon,looplat)
                 nddsf_output(tindex,looplon,looplat) = rad_norm *   &
                      zen_norm(iout) * nddsf_in(itime,looplon,looplat)
                 vbdsf_output(tindex,looplon,looplat) = rad_norm *   &
                      zen_norm(iout) * vbdsf_in(itime,looplon,looplat)
                 vddsf_output(tindex,looplon,looplat) = rad_norm *   &
                      zen_norm(iout) * vddsf_in(itime,looplon,looplat)
              enddo

              time = time + time_inc_in
           enddo
           
        enddo               ! lat, lon
     enddo

     ! Open output file for short wave
     write(outfile,'(a,i4.4,a)')trim(fdir)//'SOUTHAM_OL2_', process_year,  &
          mnames(process_month)//'.h5'
     call shdf5_open_f(trim(outfile),'W',1)

     ndims_hdf = 3
     idims_hdf(1) = tsize_output
     idims_hdf(2) = nlon
     idims_hdf(3) = nlat

     call shdf5_orec_f(ndims_hdf, idims_hdf, 'nbdsf', rvara=nbdsf_output)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'nddsf', rvara=nddsf_output)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'vbdsf', rvara=vbdsf_output)
     call shdf5_orec_f(ndims_hdf, idims_hdf, 'vddsf', rvara=vddsf_output)
     
     !print*,nbdsf_output(1:24,27,17)

     ! Close the short wave file
     call shdf5_close_f()

     ! Deallocating after every month
     deallocate(nbdsf_output)
     deallocate(nddsf_output)
     deallocate(vbdsf_output)
     deallocate(vddsf_output)

     deallocate(nbdsf_in)
     deallocate(nddsf_in)
     deallocate(vbdsf_in)
     deallocate(vddsf_in)

     deallocate(prate_in)
     deallocate(dlwrf_in)
     deallocate(uwnd_in)
     deallocate(vwnd_in)

     deallocate(air_in)
     deallocate(shum_in)
     deallocate(pres_in)
     deallocate(hgt_in)


  end do monthloop

  deallocate(var_in)
  deallocate(zen_norm)

  return

end subroutine make_the_month
!==================================================================
subroutine calculate_bin_normalization(year, month, date, time, lat,   &
     lon, frqin, frqout_rad, out_per_in, &
     radnorm, zennorm)

  implicit none

  integer, intent(in) :: month
  integer, intent(in) :: date
  integer, intent(in) :: year
  real, intent(in) :: lat
  real, intent(in) :: lon
  real, intent(in) :: frqin
  real, intent(in) :: frqout_rad
  real, intent(in) :: time
  integer, intent(in) :: out_per_in ! number of output steps in an input step
  
  real, parameter :: integration_inc=450.0  ! time step for integration
  real :: out_per_integ
  real :: in_per_integ
  integer :: n_integ_bins
  integer :: integ_per_out ! number of integration steps in an output step
  real :: normbin
  integer :: i
  integer :: iz
  integer :: ih
  real :: time_met
  real :: cosz
  real, external :: zen
  real :: norm_bin
  real :: max_cosz
  integer :: jmax_cosz
  real, intent(out) :: radnorm
  real, intent(out), dimension(out_per_in) :: zennorm

  n_integ_bins = nint(frqin / integration_inc)
  integ_per_out  = nint(frqout_rad / integration_inc)
  in_per_integ = integration_inc / frqin
  out_per_integ = integration_inc / frqout_rad
  
  norm_bin = 0.0
  zennorm(:) = 0.0
  max_cosz = -99.9

  time_met = time + 0.5 * integration_inc
  iz = 0 ! how many integration bins have we done in the current output bin?
  ih = 1 ! indicates which output bin we are in
  do i = 1, n_integ_bins ! loop over all integration bins in an input bin
     iz = iz + 1
     if(iz > integ_per_out)then
        iz = 1
        ih = ih + 1
     endif
     cosz = zen(month, date, year, time_met, lat, lon)
     if( cosz > max_cosz)then
        max_cosz = cosz
        jmax_cosz = ih
     endif
     cosz = max(0.0, cosz)
     zennorm(ih) = zennorm(ih) + cosz 
     norm_bin = norm_bin + cosz
     time_met = time_met + integration_inc
  enddo

  if(norm_bin > 0.0)then
     ! norm_bin is the average of cosz over the input bin.
     norm_bin = norm_bin * in_per_integ
     ! average of cosz over output bin.
     zennorm(:) = zennorm(:) * out_per_integ
     radnorm = 1.0 / norm_bin
  else
     radnorm = 1.0
     zennorm(:) = 0.0
     zennorm(jmax_cosz) = 1.0
  endif
  
  return
end subroutine calculate_bin_normalization

!-------------------------------------------------------------------
real function zen(month,date,year,time,lat,lon)
  implicit none
  integer month,date,year
  real time,lat,lon
  real cosz,pi180,declin,sdec,cdec,d0,d02,solfac,dayhr,radlat
  integer jday,julday
  real tempvar,tempvar2
  real cslcsd,snlsnd,gglon,dayhrr,hrangl,time_use
  
  jday = julday(month,date,year)
  
  time_use = time
  if(time.lt.0.0)then
     jday = jday - 1
     time_use = time + 86400.0
  endif
  if(time.gt.86400.0)then
     jday = jday + 1
     !         time_use = time - 86400.0
  endif
  
  pi180 = 3.1415/180.0
  declin = -23.5 * cos(6.283 / 365. * (jday + 9)) * pi180
  sdec = sin(declin)
  cdec = cos(declin)
  
  d0 = 6.2831853 * float(jday-1) / 365.0
  d02 = d0 * 2.0
  solfac = 1.000110 + 0.034221 * cos (d0) + 0.001280 * sin(d0) & 
       + 0.000719 * cos(d02) + 0.000077 * sin(d02)
  
  dayhr = time_use / 3600.0
  
  radlat = lat * pi180
  cslcsd = cos(radlat) * cdec
  snlsnd = sin(radlat) * sdec
  
  gglon = lon
  dayhrr = dayhr+gglon/15.0+24.0
  tempvar = (dayhr+gglon/15.0+24.0)/24.0 
  if(tempvar.gt.0.0)tempvar2 = 24.0*int(tempvar)
  if(tempvar.lt.0.0)tempvar2 = 24.0*(int(tempvar)+1)
  dayhrr = dayhrr - tempvar2
  
  hrangl = 15. * (dayhrr - 12.) * pi180
  
  cosz = snlsnd + cslcsd * cos(hrangl)
  
  zen=cosz

  return
end function zen

!----------------------------------------------------------------------

integer function julday(imonth,iday,iyear)
  
  integer imonth,iday,iyear
  
  integer julian_day
  integer           :: febdays
  logical, external :: isleap
  
  if (isleap(iyear)) then
     febdays = 29
  else
     febdays = 28
  end if
  
  
  julian_day= iday & 
       + min(1,max(0,imonth-1))*31 &  
       + min(1,max(0,imonth-2))*febdays &
       + min(1,max(0,imonth-3))*31   &
       + min(1,max(0,imonth-4))*30 &
       + min(1,max(0,imonth-5))*31   &
       + min(1,max(0,imonth-6))*30   &
       + min(1,max(0,imonth-7))*31   &
       + min(1,max(0,imonth-8))*31   &
       + min(1,max(0,imonth-9))*30   &
       + min(1,max(0,imonth-10))*31  &
       + min(1,max(0,imonth-11))*30  &
       + min(1,max(0,imonth-12))*31
      
  julday= julian_day
  return
end function julday

!-------------------------------------------------------------------
real function dirfrac(jday,time,rshort,lat)
  implicit none
  integer jday
  real dayangle,eccen,solartime,ket,rshort,time,lat,declination
  
  dayangle=2.0*3.14159*(jday-1.0)/365.2425
  declination=0.006918-0.399912*cos(dayangle)  &
       +0.070257*sin(dayangle)-0.006758*cos(2.0*dayangle) & 
       +0.000907*sin(2.0*dayangle)-0.002697*cos(3.0*dayangle) &  
       +0.00148*sin(3.0*dayangle)  
  eccen=1.00011+0.034221*cos(dayangle)+0.00128*sin(dayangle) &
       +0.000719*cos(2.0*dayangle)+0.000077*sin(2.0*dayangle)
  solartime=time*24.0-6.56-12.0
  ket=1367.0*eccen*(cos(declination)*cos(lat/360.0*2.0*3.14159) &
       *cos(15.0/360.0*2.0*3.14159*(solartime))  &
       +sin(declination)*sin((lat)/360.0*2.0*3.14159))
  dirfrac=rshort/ket

  if (dirfrac.gt.0.9) dirfrac=0.9
  if (dirfrac.lt.0.0) dirfrac=0.0
  return
end function dirfrac


!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
logical function isleap(year)
   !This function runs a check on whether the year is leap or not, based 
   ! on Gregorian calendar
   integer, intent(in) :: year
   isleap = (mod(year,400) == 0) .or.  &
            (mod(year,4) == 0 .and. mod(year,100) /= 0)

   return
end function isleap
!==========================================================================================!
!==========================================================================================!


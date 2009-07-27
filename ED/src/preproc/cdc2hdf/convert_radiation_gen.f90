!==========================================================================================!
!==========================================================================================!
! SUMMARY:                                                                                 !
!   1.  This program reads in the following ASCII fields:                                  !
!      - precipitation                                                                     !
!      - downward long wave radiation at the surface                                       !
!      - NIR beam                                                                          !
!      - NIR diffuse                                                                       !
!      - PAR beam                                                                          !
!      - PAR diffuse.                                                                      !
!      - Air temperature                                                                   !
!      - Pressure                                                                          !
!      - RELATIVE humidity                                                                 !
!      - Zonal wind                                                                        !
!      - Meridional wind                                                                   !
!   2.  The data are then written out in HDF5 format.                                      !
!   3.  The solar radiation fields are interpolated from longer (typically                 !
!        6-hourly resolution) to shorter time steps (typically about 15-60 minutes).       !
!        The interpolation conserves the total integral of radiation and                   !
!        accounting for variation in the solar zenith angle.                               !
!   4.  The interpolated solar radiation fields are then written out in HDF5.              !
!------------------------------------------------------------------------------------------!
! VARIABLES YOU MAY NEED TO CHANGE:                                                        !
!   first_year - first year to process                                                     !
!   last_year  - last year to process                                                      !
!   fdir       - directory of the input data (output files are placed in the               !
!                same directory.                                                           !
!   frqin      - Frequency at which input data was written out [seconds].                  !
!   frqout_rad - Frequency at which you want the interpolated output solar radiation data  !
!                written [seconds].                                                        !
!------------------------------------------------------------------------------------------!
program main
   implicit none

   !---------------------------------------------------------------------------------------!
   ! USER CONTROL SECTION                                                                  !
   !---------------------------------------------------------------------------------------!
   integer           , parameter :: first_year=1948
   integer           , parameter :: last_year=2008
   character(len=256), parameter :: fdir  = './ascii'
   character(len=256), parameter :: dpref = 'SOUTHAM'
   real              , parameter :: frqin = 21600.0     !  Frequency at which input data 
                                                        !      was written out.
   real              , parameter :: frqout_rad = 3600.0 !  Frequency at which you want the 
                                                        !      output radiation written.
   real              , parameter :: ref_hgt = 80.0      ! Reference height for Temperature,
                                                        !      humidity, and winds.
   !---------------------------------------------------------------------------------------!
   integer :: process_year

   do process_year = first_year,last_year
      call make_the_month(process_year,fdir,dpref,frqin,frqout_rad,ref_hgt)
   end do

end program main
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine make_the_month(process_year,fdir,dpref,frqin,frqout_rad,ref_hgt)

   use hdf5_utils
   use therm_lib   , only : eslif    ! ! function
   use consts_coms , only : ep       & ! intent(in)
                          , day_sec  ! ! intent(in)
   implicit none
   !----- Local constants. ----------------------------------------------------------------!
   character(len=3), dimension(12), parameter  :: mnames = (/'JAN','FEB','MAR','APR'       &
                                                            ,'MAY','JUN','JUL','AUG'       &
                                                            ,'SEP','OCT','NOV','DEC'/)
   integer                        , parameter  :: montha = 1
   integer                        , parameter  :: monthz = 12
   integer                        , parameter  :: ndh5   = 3
   !----- Arguments. ----------------------------------------------------------------------!
   real            , intent(in) :: frqin
   real            , intent(in) :: frqout_rad
   character(len=*), intent(in) :: fdir
   character(len=*), intent(in) :: dpref
   integer         , intent(in) :: process_year
   real            , intent(in) :: ref_hgt
   !----- Local variables. ----------------------------------------------------------------!
   character(len=256)                     :: infile
   character(len=256)                     :: outfile
   integer, dimension(ndh5)               :: idims_hdf
   integer                                :: idate
   integer                                :: ierr
   integer                                :: ilat
   integer                                :: ilon
   integer                                :: iout
   integer                                :: itime
   integer                                :: otime
   integer                                :: maxday
   integer                                :: ndays
   integer                                :: nlat
   integer                                :: nlon
   integer                                :: process_month
   integer                                :: times_per_day_in
   integer                                :: tsize_input
   integer                                :: tsize_output
   integer                                :: tsize_out_per_in
   real   , dimension(:,:,:), allocatable :: air_in
   real   , dimension(:,:,:), allocatable :: pres_in
   real   , dimension(:,:,:), allocatable :: shum_in
   real   , dimension(:,:,:), allocatable :: uwnd_in
   real   , dimension(:,:,:), allocatable :: hgt_in
   real   , dimension(:,:,:), allocatable :: vwnd_in
   real   , dimension(:,:,:), allocatable :: prate_in
   real   , dimension(:,:,:), allocatable :: dlwrf_in
   real   , dimension(:,:,:), allocatable :: nbdsf_in
   real   , dimension(:,:,:), allocatable :: nddsf_in
   real   , dimension(:,:,:), allocatable :: vbdsf_in
   real   , dimension(:,:,:), allocatable :: vddsf_in
   real   , dimension(:,:,:), allocatable :: nbdsf_output
   real   , dimension(:,:,:), allocatable :: nddsf_output
   real   , dimension(:,:,:), allocatable :: vbdsf_output
   real   , dimension(:,:,:), allocatable :: vddsf_output
   real   , dimension(:)    , allocatable :: zen_norm
   real                                   :: deltax
   real                                   :: deltay
   real                                   :: esat
   real                                   :: lat
   real                                   :: lon
   real                                   :: lat_min
   real                                   :: lon_min
   real                                   :: pvap
   real                                   :: rad_norm
   real                                   :: rhum
   real                                   :: time
   real                                   :: time_inc_in
   !----- External functions. -------------------------------------------------------------!
   integer, external :: julday
   logical, external :: isleap
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !  Find the number of times per day in the input and in the especial output
   !---------------------------------------------------------------------------------------!
   times_per_day_in = nint(day_sec/frqin)
   tsize_out_per_in = nint(frqin / frqout_rad)
   allocate(zen_norm(tsize_out_per_in))



   !----- Loop over months. ---------------------------------------------------------------!
   monthloop: do process_month = montha,monthz

      !----- Find the number of days the month has. ---------------------------------------!
      select case (maxday)
      case (1,3,5,7,8,10,12)
         maxday = 31
      case (4,6,9,11)
         maxday = 30
      case(2)
         if (isleap(process_year)) then
            maxday = 29
         else
            maxday = 28
         end if
      end select

      !----- Maximum number of input times in a single month. -----------------------------!
      tsize_input = maxday * times_per_day_in
    
      !----- Maximum number of output times in a single month. ----------------------------!
      tsize_output = maxday * nint(day_sec/frqout_rad)

      !----- Open and read CDC header file. -----------------------------------------------!
      infile = trim(fdir)//'/'//trim(dpref)//'_HEADER'
      open(12,file=trim(infile),form='formatted',status='old',iostat=ierr)

      if (ierr /= 0) then
         write (unit=*,fmt='(a)') '===== FATAL ERROR! ===================================='
         write (unit=*,fmt='(a)') ' -> Could not open header file '//trim(infile)//'...'
         write (unit=*,fmt='(a)') ' -> Make sure the file exists...'
         write (unit=*,fmt='(a)') '======================================================='
      end if

      read(unit=12,fmt=*) nlat,nlon
      read(unit=12,fmt=*) deltay,deltax,lat_min,lon_min
      close(12)

      !----- Allocate the variables with the correct time/longitude/latitude sizes. -------!
      allocate(air_in        (tsize_input,nlon,nlat))
      allocate(pres_in       (tsize_input,nlon,nlat))
      allocate(shum_in       (tsize_input,nlon,nlat))
      allocate(uwnd_in       (tsize_input,nlon,nlat))
      allocate(vwnd_in       (tsize_input,nlon,nlat))
      allocate(prate_in      (tsize_input,nlon,nlat))
      allocate(hgt_in        (tsize_input,nlon,nlat))
      allocate(dlwrf_in      (tsize_input,nlon,nlat))
      allocate(nbdsf_in      (tsize_input,nlon,nlat))
      allocate(nddsf_in      (tsize_input,nlon,nlat))
      allocate(vbdsf_in      (tsize_input,nlon,nlat))
      allocate(vddsf_in      (tsize_input,nlon,nlat))
      allocate(nbdsf_output (tsize_output,nlon,nlat))
      allocate(nddsf_output (tsize_output,nlon,nlat))
      allocate(vbdsf_output (tsize_output,nlon,nlat))
      allocate(vddsf_output (tsize_output,nlon,nlat))

      !---- Open input data file (on the coarse time step). -------------------------------!
      print*,'trying time: ',process_month,process_year
      write(infile,fmt='(a,i4.4,a)') trim(fdir)//'/'//trim(dpref)//'_',process_year        &
                                    ,mnames(process_month)//'.dat'
      open(unit=12,file=trim(infile),form='formatted',status='old')

      !----- Read in data and convert units on the radiation. -----------------------------!
      rlatloop: do ilat=1,nlat
         rlonloop: do ilon=1,nlon
            rtimeloop: do itime = 1,tsize_input
               read(unit=12,fmt=*)  prate_in(itime,ilon,ilat),dlwrf_in(itime,ilon,ilat)    &
                                   ,nbdsf_in(itime,ilon,ilat),nddsf_in(itime,ilon,ilat)    &
                                   ,vbdsf_in(itime,ilon,ilat),vddsf_in(itime,ilon,ilat)    &
                                   ,air_in(itime,ilon,ilat)  ,pres_in(itime,ilon,ilat)     &
                                   ,rhum                     ,uwnd_in(itime,ilon,ilat)     &
                                   ,vwnd_in(itime,ilon,ilat)

               !----- Height is constant, assigning the reference height. -----------------!
               hgt_in(itime,ilon,ilat)   = ref_hgt

               !----- Preventing negative values in some variables. -----------------------!
               prate_in(itime,ilon,ilat) = max(0.0,prate_in(itime,ilon,ilat))
               dlwrf_in(itime,ilon,ilat) = max(0.0,dlwrf_in(itime,ilon,ilat))
               nbdsf_in(itime,ilon,ilat) = max(0.0,nbdsf_in(itime,ilon,ilat))
               nddsf_in(itime,ilon,ilat) = max(0.0,nddsf_in(itime,ilon,ilat))
               vbdsf_in(itime,ilon,ilat) = max(0.0,vbdsf_in(itime,ilon,ilat))
               vddsf_in(itime,ilon,ilat) = max(0.0,vddsf_in(itime,ilon,ilat))
               rhum                      = min(1.0,max(0.0,0.01*rhum))

               !----- Converting relative humidity to specific humidity. ------------------!
               esat = eslif(air_in(itime,ilon,ilat))
               pvap = esat * rhum
               shum_in(itime,ilon,ilat) = ep * pvap                                        &
                                        / (pres_in(itime,ilon,ilat) - (1-ep) * pvap)
               
            end do rtimeloop
         end do rlonloop
      end do rlatloop

      close (unit=12,status='keep') 

      !----- Open output file for most variables. -----------------------------------------!
      write(outfile,'(a,i4.4,a)') trim(fdir)//'/'//trim(dpref)//'_OL1_', process_year      &
                                 ,mnames(process_month)//'.h5'
      call shdf5_open_f(trim(outfile),'W',1)

      !----- Write out the non-interpolated variables. ------------------------------------!
      idims_hdf(1) = tsize_input
      idims_hdf(2) = nlon
      idims_hdf(3) = nlat

      call shdf5_orec_f(ndh5, idims_hdf, 'hgt'  , rvara=hgt_in  )
      call shdf5_orec_f(ndh5, idims_hdf, 'tmp'  , rvara=air_in  )
      call shdf5_orec_f(ndh5, idims_hdf, 'pres' , rvara=pres_in )
      call shdf5_orec_f(ndh5, idims_hdf, 'sh'   , rvara=shum_in )
      call shdf5_orec_f(ndh5, idims_hdf, 'ugrd' , rvara=uwnd_in )
      call shdf5_orec_f(ndh5, idims_hdf, 'vgrd' , rvara=vwnd_in )
      call shdf5_orec_f(ndh5, idims_hdf, 'prate', rvara=prate_in)
      call shdf5_orec_f(ndh5, idims_hdf, 'dlwrf', rvara=dlwrf_in)
      
      !----- Close the output file that does not need time interpolation. -----------------!
      call shdf5_close_f()


      !------------------------------------------------------------------------------------!
      !      We now interpolate the input dataset to a more frequent time interval.  This  !
      ! interpolation conserves total integral of radiation and accounts for variation in  !
      ! the solar zenith angle.                                                            !
      !------------------------------------------------------------------------------------!
      ilatloop: do ilat=1,nlat
         !----- Latitude is from north to south. ------------------------------------------!
         lat = lat_min + (nlat - ilat) * deltay

         ilonloop: do ilon=1,nlon
            !----- Longitude is from west to east. ----------------------------------------!
            lon = lon_min + (ilon-1) * deltax

            !----- Initialise the elapsed time and time increment in this month [sec]. ----!
            time        = 0.0  
            time_inc_in = day_sec / real(times_per_day_in) ! time inc [sec]

            itimeloop: do itime=1, tsize_input

               !---------------------------------------------------------------------------!
               !     Compute the day of month and the number of seconds elapsed in this    !
               ! month.                                                                    !
               !---------------------------------------------------------------------------!
               idate = int((itime-1)/times_per_day_in)+1

               !----- Get normalization factors for this bin. -----------------------------!
               call calculate_bin_normalization(process_year,process_month,idate,time,lat  & 
                                               ,lon,frqin,frqout_rad,tsize_out_per_in      &
                                               ,rad_norm,zen_norm)

               otimeloop: do iout = 1, tsize_out_per_in
                  otime = (itime - 1) * tsize_out_per_in + iout
                  
                  nbdsf_output(otime,ilon,ilat) = rad_norm * zen_norm(iout)                &
                                                * nbdsf_in(itime,ilon,ilat)
                  nddsf_output(otime,ilon,ilat) = rad_norm * zen_norm(iout)                &
                                                * nddsf_in(itime,ilon,ilat)
                  vbdsf_output(otime,ilon,ilat) = rad_norm * zen_norm(iout)                &
                                                * vbdsf_in(itime,ilon,ilat)
                  vddsf_output(otime,ilon,ilat) = rad_norm * zen_norm(iout)                &
                                                * vddsf_in(itime,ilon,ilat)
               end do otimeloop

               time = time + time_inc_in

            end do itimeloop
         end do ilonloop
      end do ilatloop

      !----- Open output file for radiation fluxes. ---------------------------------------!
      write(outfile,'(a,i4.4,a)') trim(fdir)//'/'//trim(dpref)//'_OL2_', process_year      &
                                 ,mnames(process_month)//'.h5'
      call shdf5_open_f(trim(outfile),'W',1)

      idims_hdf(1) = tsize_output
      idims_hdf(2) = nlon
      idims_hdf(3) = nlat

      call shdf5_orec_f(ndh5, idims_hdf, 'nbdsf', rvara=nbdsf_output)
      call shdf5_orec_f(ndh5, idims_hdf, 'nddsf', rvara=nddsf_output)
      call shdf5_orec_f(ndh5, idims_hdf, 'vbdsf', rvara=vbdsf_output)
      call shdf5_orec_f(ndh5, idims_hdf, 'vddsf', rvara=vddsf_output)

      !----- Close the short wave file. ---------------------------------------------------!
      call shdf5_close_f()

      !----- Deallocating after every month. ----------------------------------------------!
      deallocate(nbdsf_output,nddsf_output,vbdsf_output,vddsf_output)
      deallocate(nbdsf_in,nddsf_in,vbdsf_in,vddsf_in,prate_in,dlwrf_in)
      deallocate(uwnd_in,vwnd_in,air_in,shum_in,pres_in,hgt_in)
   end do monthloop

   deallocate(zen_norm)

   return
end subroutine make_the_month
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine calculate_bin_normalization(year,month,date,time,lat,lon,frqin,frqout_rad       &
                                      ,out_per_in,radnorm,zennorm)
   implicit none

   !----- Local constants. ----------------------------------------------------------------!
   real                       , parameter   :: integration_inc = 450.0  ! integ. time step
   !----- Arguments. ----------------------------------------------------------------------!
   integer                    , intent(in)  :: month         ! Month
   integer                    , intent(in)  :: date          ! Day
   integer                    , intent(in)  :: year          ! Year
   integer                    , intent(in)  :: out_per_in    ! number of output steps in an
                                                             !    input step
   real                       , intent(in)  :: lat           ! Latitude
   real                       , intent(in)  :: lon           ! Longitude
   real                       , intent(in)  :: frqin         ! Input time interval
   real                       , intent(in)  :: frqout_rad    ! Output time interval for
                                                             !    the radiation variables
   real                       , intent(in)  :: time          ! Time
   real                       , intent(out) :: radnorm       ! Normalised radiation
   real, dimension(out_per_in), intent(out) :: zennorm       ! Normalised zenithal angle
   !----- Local variables. ----------------------------------------------------------------!
   real                                     :: out_per_integ ! Number of output steps per
                                                             !   integration step
   real                                     :: in_per_integ  ! Number of output steps per
                                                             !   integration step
   integer                                  :: n_integ_bins  ! Number of integration bins
   integer                                  :: integ_per_out ! number of integration steps 
                                                             !    in an output step
   integer                                  :: i
   integer                                  :: iz
   integer                                  :: ih
   real                                     :: time_met
   real                                     :: cosz
   real                                     :: norm_bin
   real                                     :: max_cosz
   integer                                  :: jmax_cosz
   !----- External functions. -------------------------------------------------------------!
   real                       , external    :: zen
   !---------------------------------------------------------------------------------------!

   !----- Finding some step countings. ----------------------------------------------------!
   n_integ_bins  = nint(frqin / integration_inc)
   integ_per_out = nint(frqout_rad / integration_inc)
   in_per_integ  = integration_inc / frqin
   out_per_integ = integration_inc / frqout_rad
  
   !----- Initialising some variables. ----------------------------------------------------!
   norm_bin   = 0.0
   zennorm(:) = 0.0
   max_cosz   = -huge(1.)

   time_met = time + 0.5 * integration_inc

   !----- How many integration bins have we done in the current output bin? ---------------!
   iz = 0 

   !----- Indicates which output bin we are in. -------------------------------------------!
   ih = 1 

   !----- Loop over all integration bins in an input bin. ---------------------------------!
   binloop: do i = 1, n_integ_bins
      iz = iz + 1
      if (iz > integ_per_out) then
         iz = 1
         ih = ih + 1
      end if

      !----- Finding the cosine of zenithal angle -----------------------------------------!
      cosz = zen(month,date,year,time_met,lat,lon)

      if (cosz > max_cosz) then
         max_cosz = cosz
         jmax_cosz = ih
      end if

      !------------------------------------------------------------------------------------!
      !     During the night time, cosz is actually negative, but for the radiation point  !
      ! of view negative cosz means that this time has zero contribution, so it is more    !
      ! accurate to make it zero rather than using the negative number.                    !
      !------------------------------------------------------------------------------------!
      cosz        = max(0.0, cosz)
      zennorm(ih) = zennorm(ih) + cosz 
      norm_bin    = norm_bin + cosz
      time_met    = time_met + integration_inc

   end do binloop

   !---------------------------------------------------------------------------------------!
   !     norm_bin is the average of cosz over the input bin.  If norm_bin is zero, this    !
   ! means we are analysing a nocturnal period of time, in which case we can't really find !
   ! the normalisation factors.                                                            !
   !---------------------------------------------------------------------------------------!
   if (norm_bin > 0.0) then
      norm_bin = norm_bin * in_per_integ
      ! average of cosz over output bin.
      zennorm(:) = zennorm(:) * out_per_integ
      radnorm = 1.0 / norm_bin
   else
      radnorm = 1.0
      zennorm(:) = 0.0
      zennorm(jmax_cosz) = 1.0
   end if
  
   return
end subroutine calculate_bin_normalization
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
real function zen(month,date,year,time,lat,lon)
   use consts_coms , only : pio1808  & ! intent(in)
                          , twopi8   & ! intent(in)
                          , day_sec  & ! intent(in)
                          , day_sec8 & ! intent(in)
                          , hr_sec8  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer     , intent(in) :: month
   integer     , intent(in) :: date
   integer     , intent(in) :: year
   real        , intent(in) :: time
   real        , intent(in) :: lat
   real        , intent(in) :: lon
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)             :: cosz
   real(kind=8)             :: declin
   real(kind=8)             :: sdec
   real(kind=8)             :: cdec
   real(kind=8)             :: dayhr
   real(kind=8)             :: radlat
   real(kind=8)             :: tempvar
   real(kind=8)             :: tempvar2
   real(kind=8)             :: cslcsd
   real(kind=8)             :: snlsnd
   real(kind=8)             :: gglon
   real(kind=8)             :: dayhrr
   real(kind=8)             :: hrangl
   real(kind=8)             :: time_use
   integer                  :: jday
   !----- External functions. -------------------------------------------------------------!
   integer     , external   :: julday
   real        , external   :: sngloff
   !---------------------------------------------------------------------------------------!
   
   !----- Finding the day of year ("Julian" day). -----------------------------------------!
   jday = julday(month,date,year)
   
   !----- Finding latitude in radians and putting longitude in double precision. ----------!
   radlat = dble(lat) * pio1808
   gglon  = dble(lon)

   !----- Since we add the time increment to time, it may happen to exceed one day... -----!
   if (time < 0.) then
      jday = jday - 1
      time_use = dble(time) + day_sec8
   elseif (time > day_sec) then
      jday = jday + 1
      time_use = dble(time) - day_sec8
   else 
      time_use = dble(time)
   end if

   !----- Finding the declination, its sine and cosine. -----------------------------------!
   declin = -2.35d1 * dcos(twopi8 / 3.65d2 * dble(jday + 9)) * pio1808
   sdec   = dsin(declin)
   cdec   = dcos(declin)

   !----- Making the time in hours. -------------------------------------------------------!
   dayhr  = time_use / hr_sec8

   cslcsd = dcos(radlat) * cdec
   snlsnd = dsin(radlat) * sdec
   
   dayhrr = dmod(dayhr+gglon/1.5d1+2.4d1,2.4d1)
   hrangl = 1.5d1 * (dayhrr - 1.2d1) * pio1808

   cosz   = snlsnd + cslcsd * cos(hrangl)
   zen    = sngloff(cosz,1.d-20)

   return
end function zen
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
integer function julday(month,day,year)
   implicit none
   !----- Local constants. ----------------------------------------------------------------!
   integer , dimension(12), parameter :: maxday = (/ 31, 28, 31, 30, 31, 30                &
                                                   , 31, 31, 30, 31, 30, 31 /)
   !----- Arguments. ----------------------------------------------------------------------!
   integer , intent(in) :: month
   integer , intent(in) :: day
   integer , intent(in) :: year
   !----- Local variables. ----------------------------------------------------------------!
   integer              :: im
   !----- External functions. -------------------------------------------------------------!
   logical , external   :: isleap
   !---------------------------------------------------------------------------------------!

   !----- Start by assigning the current day of the current month. ------------------------!
   julday = day

   !----- Then add the number of days from all months before this one. --------------------!
   do im=1,month-1
      julday = julday + maxday(im)
   end do

   !----- Finally, adding an extra day for leap years, only if it past February. ----------!
   if (isleap(year) .and. month > 2 ) then
      julday = julday + 1
   end if

   return
end function julday
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function runs a check on whether the year is leap or not, based on Gregorian    !
! calendar (97 leap years every 400 years).                                                !
!------------------------------------------------------------------------------------------!
logical function isleap(year)
   integer, intent(in) :: year
   isleap = (mod(year,400) == 0) .or. (mod(year,4) == 0 .and. mod(year,100) /= 0)
   return
end function isleap
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function converts the double precision variable into single, in a way to prevent !
! floating point exception when they are tiny.  In case the number is too small, less than !
! off, then the output value is flushed to 0.                                              !
!------------------------------------------------------------------------------------------!
real function sngloff(x,off)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=8), intent(in) :: x
   real(kind=8), intent(in) :: off
   !---------------------------------------------------------------------------------------!
   
   if (abs(x) < off) then
      sngloff = 0.
   else
      sngloff = sngl(x)
   end if

   return
end function sngloff
!==========================================================================================!
!==========================================================================================!

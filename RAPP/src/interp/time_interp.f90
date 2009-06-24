!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the time interpolation of solar fluxes.                 !
!------------------------------------------------------------------------------------------!
subroutine time_interp()
   use mod_model  , only : this_time     ! ! intent(inout)
   use mod_ioopts , only : radfrq        & ! intent(in)
                         , radratio      ! ! intent(in)
   use mod_interp , only : interp_buffer & ! intent(inout)
                         , mxgauss       & ! intent(in)
                         , mygauss       ! ! intent(in)
   use mod_grid   , only : grid_g        & ! intent(in)
                         , sstp          & ! intent(in)
                         , t_1st         & ! intent(in)
                         , tlast         ! ! intent(in)
   use mod_ncep   , only : ncep_g        ! ! intent(inout)
   use mod_time   , only : copy_time     ! ! subroutine
   use rconstants , only : day_sec       ! ! intent(in)
   implicit none 
   !----- Local variables. ----------------------------------------------------------------!
   integer :: t1
   integer :: t2
   integer :: ta
   integer :: tt
   integer :: tz
   integer :: xx
   integer :: yy
   integer :: it
   real    :: deltat
   real    :: fracday
   real    :: rad_norm
   !---------------------------------------------------------------------------------------!

   !----- This is the time advance in days. -----------------------------------------------!
   deltat = radfrq / day_sec

   write (unit=*,fmt='(a)') '     - Running the time interpolation...'

   !---------------------------------------------------------------------------------------!
   ! STEP 1.  Copy the data from the info table to the ncep structure.                     !
   !---------------------------------------------------------------------------------------!
   do t1=1,sstp(1)
      tt = t_1st(1) + t1 - 1
      call copy_time(this_time(tt,1),ncep_g(1)%when(t1))
   end do


   !---------------------------------------------------------------------------------------!
   ! STEP 2.  Fill the time array for the second grid, using the first grid as reference.  !
   !---------------------------------------------------------------------------------------!
   st1loop: do t2=1,sstp(2)
      t1   = 1 + (t2-1)/radratio
      it   = mod((t2-1),radratio)

      fracday = ncep_g(1)%when(t1)%fracday + it * deltat
      
      !------------------------------------------------------------------------------------!
      !    Finding the new day fraction and new day of year.  If we are dealing with       !
      ! December, this may become temporarily inconsistent, but it will be fixed in the    !
      ! fill_time_info subroutine.                                                         !
      !------------------------------------------------------------------------------------!
      ncep_g(2)%when(t2)%fracday = mod(fracday,1.)
      ncep_g(2)%when(t2)%doy     = ncep_g(1)%when(t1)%doy + int(fracday)
      ncep_g(2)%when(t2)%year    = ncep_g(1)%when(t1)%year

      call fill_time_info(ncep_g(2)%when(t2))
   end do st1loop

   !---------------------------------------------------------------------------------------!
   ! STEP 3.  Find the normalised factors for each bin of time, then perform interpol-     !
   !          ation.                                                                       !
   !---------------------------------------------------------------------------------------!
   st2loop: do t1=1,sstp(1)
      write (unit=*,fmt='(a,1x,2a)')                                                       &
                             '         [|] Interpolating:',ncep_g(1)%when(t1)%timestr,'...'

      yloop: do yy=1,mygauss
         xloop: do xx=1,mxgauss
            !----- Get normalisation factors for this time. -------------------------------!
            call calculate_normalised_factors(ncep_g(1)%when(t1),grid_g(1)%lon(xx,yy)      &
                                             ,grid_g(1)%lat(xx,yy),rad_norm)

            !----- Map the time equivalence, then fill the output arrays. -----------------!
            ta = (t1-1) * radratio + 1
            tz = ta + radratio - 1
            g2tloop: do tt=ta,tz
               t2 = tt - ta + 1
               ncep_g(2)%nbdsf(xx,yy,tt) = rad_norm * interp_buffer%zen_norm(t2)           &
                                         * ncep_g(1)%nbdsf(xx,yy,t1)
               ncep_g(2)%nddsf(xx,yy,tt) = rad_norm * interp_buffer%zen_norm(t2)           &
                                         * ncep_g(1)%nddsf(xx,yy,t1)
               ncep_g(2)%vbdsf(xx,yy,tt) = rad_norm * interp_buffer%zen_norm(t2)           &
                                         * ncep_g(1)%vbdsf(xx,yy,t1)
               ncep_g(2)%vddsf(xx,yy,tt) = rad_norm * interp_buffer%zen_norm(t2)           &
                                         * ncep_g(1)%vddsf(xx,yy,t1)
            end do g2tloop
         end do xloop
      end do yloop
   end do st2loop
   return
end subroutine time_interp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine calculate_normalised_factors(when,lon,lat,rad_norm)
   use mod_time    , only : time_stt      ! ! structure
   use mod_interp  , only : interp_buffer & ! intent(inout)
                          , dtinc_fday    ! ! intent(in)
   use mod_ioopts  , only : radratio      & ! intent(in)
                          , nsteps        & ! intent(in)
                          , nrads         & ! intent(in)
                          , nstepsi       & ! intent(in)
                          , nradsi        ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(time_stt) , intent(in)  :: when
   real           , intent(in)  :: lon
   real           , intent(in)  :: lat
   real           , intent(out) :: rad_norm
   !----- Local variables. ----------------------------------------------------------------!
   type(time_stt)               :: btime      ! Bin reference time
   real                         :: cosz       ! Cosine of zenithal angle
   real                         :: norm_bin   ! Normalisation factor
   real                         :: fracday    ! Fraction of one day
   integer                      :: tbmax_cosz ! Time in which the maximum cosz has happened
   integer                      :: ts         ! Time step counter
   integer                      :: tb         ! Output bin counter
   !----- External functions. -------------------------------------------------------------!
   real           , external    :: zen
   !---------------------------------------------------------------------------------------!

   !----- Initialising some variables. ----------------------------------------------------!
   norm_bin               = 0.0
   interp_buffer%zen_norm = 0.0

   
   !---------------------------------------------------------------------------------------!
   !    Finding the first bin day fraction and new day of year.  If we are dealing with    !
   ! December, this may become temporarily inconsistent, but it will be fixed in the       !
   ! fill_time_info subroutine.                                                            !
   !---------------------------------------------------------------------------------------!
   fracday       = when%fracday + 0.5 * dtinc_fday
   btime%fracday = mod(fracday,1.)
   btime%doy     = when%doy + int(fracday)
   btime%year    = when%year
   call fill_time_info(btime)


   timesteploop: do ts=1,nsteps
      tb = (ts-1)/nrads + 1

      !----- Finding the cosine of the zenithal angle. ------------------------------------!
      cosz = zen(btime%doy,btime%fracday,lon,lat)

      !------------------------------------------------------------------------------------!
      !     During the night time, cosz is actually negative, but for the radiation point  !
      ! of view negative cosz means that this time has zero contribution, so it is more    !
      ! accurate to make it zero rather than using the negative number.                    !
      !------------------------------------------------------------------------------------!
      cosz     = max(0.0, cosz)
      norm_bin = norm_bin + cosz
      interp_buffer%zen_norm(tb) = interp_buffer%zen_norm(tb) + cosz 

      !------------------------------------------------------------------------------------!
      !    Finding the next bin day fraction and new day of year.  If we are dealing with  !
      ! December, this may become temporarily inconsistent, but it will be fixed in the    !
      ! fill_time_info subroutine.                                                         !
      !------------------------------------------------------------------------------------!
      fracday       = btime%fracday + dtinc_fday
      btime%fracday = mod(fracday,1.)
      btime%doy     = btime%doy + int(fracday)
      call fill_time_info(btime)

   end do timesteploop
   !---------------------------------------------------------------------------------------!
   !     norm_bin is the average of cosz over the input bin.  If norm_bin is zero, this    !
   ! means we are analysing a nocturnal period of time, in which case we can't really find !
   ! the normalisation factors.                                                            !
   !---------------------------------------------------------------------------------------!
   if (norm_bin > 0.0) then
      norm_bin = norm_bin * nstepsi
      !----- Average of cosz over output bin. ---------------------------------------------!
      interp_buffer%zen_norm = interp_buffer%zen_norm * nradsi
      !----- Normalisation factor. --------------------------------------------------------!
      rad_norm = 1.0 / norm_bin
   else
      !----- Making up something, but the value will be zeroed. ---------------------------!
      tbmax_cosz = maxloc(interp_buffer%zen_norm(:),dim=1)
      rad_norm   = 1.0
      interp_buffer%zen_norm             = 0.0
      interp_buffer%zen_norm(tbmax_cosz) = 1.0
   end if

   return
end subroutine calculate_normalised_factors
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will fill the time structure with all information, given the day of   !
! year, fraction of day, and year are previously assigned.                                 !
!------------------------------------------------------------------------------------------!
subroutine fill_time_info(when)
   use mod_time  , only : time_stt ! ! intent(in)
   use rconstants, only : day_hr   & ! intent(in)
                        , day_min  & ! intent(in)
                        , day_sec  & ! intent(in)
                        , hr_min   & ! intent(in)
                        , hr_sec   & ! intent(in)
                        , min_sec  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(time_stt)    , intent(inout) :: when
   !----- Local variables. ----------------------------------------------------------------!
   integer                           :: days_year
   !----- External functions. -------------------------------------------------------------!
   character(len=15) , external      :: grads_dtstamp
   character(len=17) , external      :: rapp_dtstamp
   integer           , external      :: v5d_datestamp
   integer           , external      :: julday1000
   logical           , external      :: isleap
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    First we check if day of year didn't exceed the maximum.  In case it did, simplify !
   ! the doy/year combination to reduce to a normal calendar.                              !
   !---------------------------------------------------------------------------------------!
   !------ Check whether the year is leap or not and assign how many days this year has. --!
   if (isleap(when%year)) then
      days_year = 366
   else
      days_year = 365
   end if
   !----- Loop over years and simplify. ---------------------------------------------------!
   do while (when%doy > days_year)
      !----- Remove the excess of days and add one year. ----------------------------------!
      when%doy  = when%doy - days_year
      when%year = when%year + 1
      !----- Update the number of days of the new year. -----------------------------------!
      if (isleap(when%year)) then
         days_year = 366
      else
         days_year = 365
      end if
   end do
   !---------------------------------------------------------------------------------------!

   !----- Find the month and day. ---------------------------------------------------------!
   call doy_2_monday(when%doy,when%year,when%month,when%day,when%mmm)

   !---------------------------------------------------------------------------------------!
   !    Find the hours, minutes, and seconds, using the day fraction.                      !
   !---------------------------------------------------------------------------------------!
   when%hour = int(mod(when%fracday * day_hr ,day_hr ))
   when%minu = int(mod(when%fracday * day_min,hr_min ))
   when%seco = int(mod(when%fracday * day_sec,min_sec))
   
   !----- Finding the elapsed time. -------------------------------------------------------!
   when%elapsed = dble(julday1000(when%month,when%day,when%year)) + dble(when%fracday)
   
   
   !----- Finding some time stamps. -------------------------------------------------------!
   when%yyyyddd    = v5d_datestamp(when%year,when%doy)
   when%gradsstamp = grads_dtstamp(when%year,when%mmm,when%day,when%hour,when%minu)
   when%timestr    = rapp_dtstamp (when%year,when%month,when%day                           &
                                  ,when%hour,when%minu,when%seco)

   return
end subroutine fill_time_info
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function will compute the zenithal angle at a given time and location.          !
!------------------------------------------------------------------------------------------!
real function zen(doy,time,lon,lat)
   use rconstants , only : pio1808  & ! intent(in)
                         , twopi8   ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer     , intent(in) :: doy      ! Day of year
   real        , intent(in) :: time     ! Time of day (days)
   real        , intent(in) :: lon      ! Longitude
   real        , intent(in) :: lat      ! Latitude
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)             :: cosz     ! Cosine of zenithal angle
   real(kind=8)             :: declin   ! Declination
   real(kind=8)             :: sdec     ! Sine of declination
   real(kind=8)             :: cdec     ! Cosine of declination
   real(kind=8)             :: dayhr    ! Hour of day
   real(kind=8)             :: radlat   ! Latitude in radians
   real(kind=8)             :: cslcsd   ! Cosine of latitude times cosine of declination
   real(kind=8)             :: snlsnd   ! Sine of latitude times sine of declination
   real(kind=8)             :: gglon    ! Longitude in double precision
   real(kind=8)             :: dayhrr   ! Hour of day in radians
   real(kind=8)             :: hrangl   ! Hour angle
   !----- External functions. -------------------------------------------------------------!
   real        , external   :: sngloff  ! Safe double to single precision converter.
   !---------------------------------------------------------------------------------------!


   !----- Finding latitude in radians and putting longitude in double precision. ----------!
   radlat = dble(lat) * pio1808
   gglon  = dble(lon)
   !---------------------------------------------------------------------------------------!


   !----- Convert the fractional hour of day (a real number). -----------------------------!
   dayhr    = 2.4d1 * dble(time)
   !---------------------------------------------------------------------------------------!


   !----- Finding the declination, its sine and cosine. -----------------------------------!
   declin = -2.35d1 * dcos(twopi8 / 3.65d2 * dble(doy + 9)) * pio1808
   sdec   = dsin(declin)
   cdec   = dcos(declin)
   !---------------------------------------------------------------------------------------!


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

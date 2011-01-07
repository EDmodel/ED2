!==========================================================================================!
!==========================================================================================!
!    This subroutine calculates phenology factors for prescribed phenology schemes.        !
!------------------------------------------------------------------------------------------!
subroutine prescribed_leaf_state(lat,imonth,iyear,doy,green_leaf_factor,leaf_aging_factor  &
                                ,phen_pars)

   use phenology_coms , only : iphenys1         & ! intent(in)
                             , iphenysf         & ! intent(in)
                             , iphenyf1         & ! intent(in)
                             , iphenyff         & ! intent(in)
                             , prescribed_phen  ! ! intent(in)
   use ed_max_dims       , only : n_pft            ! ! intent(in)
   use pft_coms       , only : phenology        ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(prescribed_phen) , intent(in)  :: phen_pars
   integer               , intent(in)  :: iyear
   integer               , intent(in)  :: doy
   real                  , intent(in)  :: lat
   integer               , intent(in)  :: imonth
   real, dimension(n_pft), intent(out) :: green_leaf_factor
   real, dimension(n_pft), intent(out) :: leaf_aging_factor
   !----- Local variables -----------------------------------------------------------------!
   integer                             :: n_recycle_years
   integer                             :: my_year
   real                                :: elongf
   real                                :: delay
   real(kind=8)                        :: elonDen
   integer                             :: pft
   !----- Local constants -----------------------------------------------------------------!
   real                  , parameter   :: elonMin = 0.02
   !---------------------------------------------------------------------------------------!
  
   !---------------------------------------------------------------------------------------!
   !     This assumes dropping/flushing based on the day of year and hemisphere.           !
   ! + Northern Hemisphere: dropping between August 1 and December 31;                     !
   !                        flushing between January 1 and July 31.                        !
   ! + Southern Hemisphere: dropping between February 1 and July 31;                       !
   !                        flushing between August 1 and January 31.                      !
   !---------------------------------------------------------------------------------------!
   if( (lat >= 0.0 .and. imonth <= 7) .or.                                                 &
       (lat < 0.0  .and. (imonth > 7 .or. imonth == 1)) )then
      
      !----- Get the year. ----------------------------------------------------------------!
      n_recycle_years = iphenysf - iphenys1 + 1

      if (iyear > iphenysf) then
         my_year = mod(iyear-iphenys1,n_recycle_years) + 1
      elseif (iyear < iphenys1) then
         my_year = n_recycle_years - mod(iphenysf-iyear,n_recycle_years)
      else
         my_year = iyear - iphenys1 + 1
      end if

      !------------------------------------------------------------------------------------!
      !      Calculate the factors.  Precalc denominator and limit rate in order to        !
      ! increase numerical stability (MCD 10/23/08).                                       !
      !------------------------------------------------------------------------------------!
      elonDen = real((phen_pars%flush_a(my_year) * real(doy)),kind=8)                      &
              ** dble(max(phen_pars%flush_b(my_year),-100.))
      elonDen = 1.0d0 / (1.0d0 + elonDen)
      if(elonDen < 0.0001d0) then
         elongf = 0.0
      else
         elongf = sngl(elonDen)
      end if
      delay = elongf     
   else
      !------------------------------------------------------------------------------------!
      !      Leaves turning color.  Get the year.                                          !
      !------------------------------------------------------------------------------------!
      n_recycle_years = iphenyff - iphenyf1 + 1
      if (iyear > iphenyff) then
         my_year = mod(iyear-iphenyf1,n_recycle_years) + 1
      elseif (iyear < iphenyf1) then
         my_year = n_recycle_years - mod(iphenyff-iyear,n_recycle_years)
      else
         my_year = iyear - iphenyf1 + 1
      end if

      !----- Calculate the factors. -------------------------------------------------------!
      elongf = 1.0                                                                         &
             / (1.0 + (phen_pars%color_a(my_year) * real(doy))**phen_pars%color_b(my_year))
      delay = 1.0                                                                          &
            / (1.0 + (phen_pars%color_a(my_year) * real(doy) * 1.095)                      &
            **phen_pars%color_b(my_year))
   end if

   if(elongf < elonMin) elongf = 0.0

   !----- Load the values for each PFT. ---------------------------------------------------!
   do pft = 1, n_pft
      select case (phenology(pft))
      case (2)
         green_leaf_factor(pft) = elongf
         leaf_aging_factor(pft) = delay
      case default
         green_leaf_factor(pft) = 1.0
         leaf_aging_factor(pft) = 1.0
      end select
   end do
  
   return
end subroutine prescribed_leaf_state
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the number of chill and warming days in a month.            !
!  + Chill days  - number of days with average temperatures below 278.15 K;                !
!  + Degree days - sum of daily average temperatures above 278.15 K.                       !
!------------------------------------------------------------------------------------------!
subroutine update_thermal_sums(month, cpoly, isi, lat)
  
   use ed_state_vars ,only : polygontype & ! structure
                           , sitetype    ! ! structure
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype) , target     :: cpoly
   integer           , intent(in) :: isi
   integer           , intent(in) :: month
   real              , intent(in) :: lat
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)    , pointer    :: csite
   integer                        :: ipa
   !---------------------------------------------------------------------------------------!



   !----- Loop over patches. --------------------------------------------------------------!
   csite => cpoly%site(isi)

   do ipa = 1,csite%npatches

      !----- Minimum monthly temperature of the site. -------------------------------------!
      cpoly%min_monthly_temp(isi) = min(cpoly%min_monthly_temp(isi)                        &
                                       ,csite%avg_daily_temp(ipa))

      !----- Warm day, so check whether it is growing season and update the degree day... -!
      if (csite%avg_daily_temp(ipa) > 278.15) then
         !----- Update dgd only for growing season. ---------------------------------------!
         if ((lat >= 0.0 .and. month <= 8) .or.                                            &
             (lat <  0.0 .and. (month <= 2 .or. month >= 7))) then
            csite%sum_dgd(ipa) = csite%sum_dgd(ipa) + (csite%avg_daily_temp(ipa)-278.15)
         !----- Warm day during dropping season, set degree sum to zero... ----------------!
         else 
            csite%sum_dgd(ipa) = 0.0
         end if
      !---- Cold day, check whether it is dropping season and update chilling days... -----!
      elseif ((lat >= 0.0 .and. (month >= 11 .or. month <= 6)) .or.                        &
              (lat <  0.0 .and.  month >= 5)                 ) then 
         csite%sum_chd(ipa) = csite%sum_chd(ipa) + 1.0
      !---- Cold day, but not during dropping season, set chilling days to zero... --------!
      else 
         csite%sum_chd(ipa) = 0.0
      end if
   end do
   
   return
end subroutine update_thermal_sums
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine updates the turnover ratio and specific leaf area, taking into       !
! account the available radiation.                                                         !
!------------------------------------------------------------------------------------------!
subroutine update_turnover(cpoly, isi)
   use ed_state_vars  , only : polygontype        & ! structure
                             , sitetype           & ! structure
                             , patchtype          ! ! structure
   use pft_coms       , only : is_tropical        & ! intent(in)
                             , sla                & ! intent(in)
                             , leaf_turnover_rate ! ! intent(in)
   use phenology_coms , only : rad_turnover_int   & ! intent(in)
                             , rad_turnover_slope & ! intent(in)
                             , vm_tran            & ! intent(in)
                             , vm_slop            & ! intent(in)
                             , vm_amp             & ! intent(in)
                             , vm_min             ! ! intent(in)
   use consts_coms    , only : day_sec            ! ! intent(in)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype) , target     :: cpoly
   integer           , intent(in) :: isi
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)    , pointer    :: csite
   type(patchtype)   , pointer    :: cpatch
   integer                        :: ipa
   integer                        :: ico
   integer                        :: ipft
   real                           :: turnover0
   real                           :: llspan0
   real                           :: vm0
   !----- Local constants -----------------------------------------------------------------!
   real              , parameter  :: tfact10 = 0.1
   real              , parameter  :: tfact60 = 1./60
   real              , save       :: radcrit
   logical           , save       :: first_time=.true.
   !---------------------------------------------------------------------------------------!

   if (first_time) then
      first_time = .false.
      radcrit = - rad_turnover_int / rad_turnover_slope
   end if

   !----- Loop over patches. --------------------------------------------------------------!
   csite => cpoly%site(isi)
   patchloop: do ipa = 1,csite%npatches
     
      cpatch => csite%patch(ipa)
      cohortloop: do ico = 1,cpatch%ncohorts

         ipft = cpatch%pft(ico)
         
 !        write(unit=*,fmt='(a,1x,es12.5)') 'Rad_avg is       =', cpoly%rad_avg
         
         !----- Update turnover mulitplier. -----------------------------------------------!
         if (cpoly%rad_avg(isi) < radcrit) then
            turnover0 = 0.01
         else
            turnover0 = min(100.                                                           &
                           , max(0.01                                                      &
                                ,rad_turnover_int+rad_turnover_slope*cpoly%rad_avg(isi)))
         end if
         
   !               write(unit=*,fmt='(a,1x,es12.5)') 'New Turnover is       =', turnover0
                  
         cpatch%turnover_amp(ico) = (1.0 - tfact10) * cpatch%turnover_amp(ico)             &
                                  +        tfact10  * turnover0

         !----- Update leaf lifespan. -----------------------------------------------------!
         if (leaf_turnover_rate(ipft) > 0.) then
            llspan0       = 12.0 / (cpatch%turnover_amp(ico) * leaf_turnover_rate(ipft))
            if (llspan0 < 2.) then
                llspan0=2.
            elseif (llspan0 > 60.) then
                llspan0 = 60.
            end if
         else
            llspan0       = 9999.
         end if
         
   !               write(unit=*,fmt='(a,1x,es12.5)') 'llspan0 is       =', llspan0
         cpatch%llspan(ico) = (1.0 - tfact60) * cpatch%llspan(ico) + tfact60 * llspan0
   !               write(unit=*,fmt='(a,1x,es12.5)') 'llspan(ico) is       =', cpatch%llspan(ico)         

         !----- Update vm_bar. ------------------------------------------------------------!
         vm0               = vm_amp / (1.0 + (cpatch%llspan(ico)/vm_tran)**vm_slop) + vm_min
         cpatch%vm_bar(ico)= (1.0 - tfact60) * cpatch%vm_bar(ico) + tfact60 * vm0

         !----- Update the specific leaf area (SLA). --------------------------------------!
         if (is_tropical(ipft)) then
            cpatch%sla(ico) =  10.0                                                        &
                            ** ( 1.6923                                                    &
                               - 0.3305 *log10(12.0 / ( cpatch%turnover_amp(ico)           &
                                                      * leaf_turnover_rate(ipft)) ) )
         else
            cpatch%sla(ico) = sla(ipft)
         end if
      end do cohortloop

   end do patchloop

   return
end subroutine update_turnover
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function computes the length of daylight for a given latitude and day of year.  !
! The result is given in minutes.                                                          !
!------------------------------------------------------------------------------------------!
real function daylength(lat,doy)

   use consts_coms , only : pio180 & ! intent(in)
                          , twopi  ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real    , intent(in) :: lat
   integer , intent(in) :: doy
   !----- Local variables -----------------------------------------------------------------!
   real                 :: arg
   !---------------------------------------------------------------------------------------!

   arg = -tan(lat * pio180) * tan(-23.5 * pio180 * cos(twopi/365.0 * (float(doy)+9.0)))

   if (arg >= 1.0) then
      daylength = 0.0
   elseif (arg <= 1.0) then
      daylength = 1440.0
   else ! if (abs(arg) < 1.0) then
      daylength = 120.0 * acos(arg) / (15.0 * pio180)
   end if

   return
end function daylength
!==========================================================================================!
!==========================================================================================!

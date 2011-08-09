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
                             , prescribed_phen  & ! intent(in)
                             , elongf_min       ! ! intent(in)
   use ed_max_dims    , only : n_pft            ! ! intent(in)
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

   if(elongf < elongf_min) elongf = 0.0

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
         
         !write(unit=*,fmt='(a,1x,es12.5)') 'Rad_avg is       =', cpoly%rad_avg
         
         !----- Update turnover mulitplier. -----------------------------------------------!
         if (cpoly%rad_avg(isi) < radcrit) then
            turnover0 = 0.01
         else
            turnover0 = min(100.                                                           &
                           , max(0.01                                                      &
                                ,rad_turnover_int+rad_turnover_slope*cpoly%rad_avg(isi)))
         end if
         
         !         write(unit=*,fmt='(a,1x,es12.5)') 'New Turnover is       =', turnover0
                  
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
         
         !         write(unit=*,fmt='(a,1x,es12.5)') 'llspan0 is       =', llspan0
         cpatch%llspan(ico) = (1.0 - tfact60) * cpatch%llspan(ico) + tfact60 * llspan0
         !         write(unit=*,fmt='(a,1x,es12.5)') 'llspan(ico) is       =', cpatch%llspan(ico)         

         !----- Update vm_bar. ------------------------------------------------------------!
         vm0               = vm_amp / (1.0 + (cpatch%llspan(ico)/vm_tran)**vm_slop) + vm_min
         cpatch%vm_bar(ico)= (1.0 - tfact60) * cpatch%vm_bar(ico) + tfact60 * vm0

         !----- Update the specific leaf area (SLA). --------------------------------------!
         if (is_tropical(ipft) .and. ipft /= 17) then
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
!      This sub-routine will assign the initial potential available water and the          !
! phenology that has been assigned.  This sub-routine should be called whenever a new      !
! cohort is created, or at the initial run (except history).  The initial running average  !
! is simply the the instantaneous soil moisture variable.  For plants other than the       !
! drought-deciduous, the potential available water is found but it doesn't control the     !
! phenology, so we assign fully flushed leaves.                                            !
!------------------------------------------------------------------------------------------!
subroutine first_phenology(cgrid)
   use ed_state_vars , only : edtype           & ! structure
                            , polygontype      & ! structure
                            , sitetype         & ! structure
                            , patchtype        ! ! structure
   use ed_therm_lib  , only : calc_veg_hcap    ! ! function
   use allometry     , only : area_indices     ! ! subroutine
   use grid_coms     , only : nzg              ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)     , target     :: cgrid          ! Current grid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype), pointer    :: cpoly          ! Current polygon
   type(sitetype)   , pointer    :: csite          ! Current site
   type(patchtype)  , pointer    :: cpatch         ! Current patch
   integer                       :: ipy            ! Polygon counter
   integer                       :: isi            ! Site counter
   integer                       :: ipa            ! Patch counter
   integer                       :: ico            ! Cohort counter
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Loop over all cohorts in this grid.                                               !
   !---------------------------------------------------------------------------------------!
   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      siteloop: do isi=1,cpoly%nsites
         csite => cpoly%site(isi)
         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)
            cohortloop: do ico=1,cpatch%ncohorts

               !---------------------------------------------------------------------------!
               !    Find the initial guess for potential available water and elongation    !
               ! factor, then compute the equilibrium biomass of active tissues and        !
               ! storage.                                                                  !
               !---------------------------------------------------------------------------!
               call pheninit_balive_bstorage(nzg,csite,ipa,ico,cpoly%ntext_soil(:,isi)     &
                                            ,cpoly%green_leaf_factor(:,isi))
               !---------------------------------------------------------------------------!


               !----- Find LAI, WPA, WAI. -------------------------------------------------!
               call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico),cpatch%bdead(ico)    &
                                ,cpatch%balive(ico),cpatch%dbh(ico),cpatch%hite(ico)       &
                                ,cpatch%pft(ico),cpatch%sla(ico),cpatch%lai(ico)           &
                                ,cpatch%wpa(ico),cpatch%wai(ico),cpatch%crown_area(ico)    &
                                ,cpatch%bsapwood(ico)) 
               !---------------------------------------------------------------------------!


               !----- Find heat capacity and vegetation internal energy. ------------------!
               call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdead(ico),cpatch%bsapwood(ico) &
                                 ,cpatch%nplant(ico),cpatch%pft(ico)                       &
                                 ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
               cpatch%leaf_energy(ico) = cpatch%leaf_hcap(ico) * cpatch%leaf_temp(ico)
               cpatch%wood_energy(ico) = cpatch%wood_hcap(ico) * cpatch%wood_temp(ico)
               call is_resolvable(csite,ipa,ico,cpoly%green_leaf_factor(:,isi))
               !---------------------------------------------------------------------------!


            end do cohortloop
         end do patchloop
      end do siteloop
   end do polyloop

   return
end subroutine first_phenology
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine will assign the initial potential available water and the          !
! phenology that has been assigned, then find the biomass of active tissues and storage    !
! that is in equilibrium with the initial soil moisture.  This sub-routine should be       !
! called whenever a new cohort is planted or recruited, or at the initial run (except      !
! history).  The initial running average is simply the the instantaneous soil moisture     !
! variable.  For plants other than the drought-deciduous, the potential available water is !
! found but it doesn't control the phenology, so we assign the biomass that matches the    !
! fully flushed leaves.                                                                    !
!------------------------------------------------------------------------------------------!
subroutine pheninit_balive_bstorage(mzg,csite,ipa,ico,ntext_soil,green_leaf_factor)
   use ed_misc_coms  , only : ivegt_dynamics      ! ! intent(in)
   use ed_state_vars , only : sitetype            & ! structure
                            , patchtype           ! ! structure
   use soil_coms     , only : soil                & ! intent(in), look-up table
                            , slz                 & ! intent(in)
                            , dslz                ! ! intent(in)
   use phenology_coms, only : spot_phen           & ! intent(in)
                            , elongf_min          ! ! intent(in)
   use pft_coms      , only : phenology           & ! intent(in)
                            , q                   & ! intent(in)
                            , qsw                 ! ! intent(in)
   use ed_max_dims   , only : n_pft               ! ! intent(in)
   use allometry     , only : dbh2bl              ! ! function
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)                  , target     :: csite             ! Current site
   integer                         , intent(in) :: mzg               ! # of soil layers
   integer                         , intent(in) :: ipa               ! Current patch
   integer                         , intent(in) :: ico               ! Cohort counter
   integer       , dimension(mzg)  , intent(in) :: ntext_soil        ! Soil texture
   real          , dimension(n_pft), intent(in) :: green_leaf_factor ! Hardwood phenology
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)               , pointer    :: cpatch     ! Current patch
   integer                                    :: k          ! Layer counter
   integer                                    :: ipft       ! PFT type
   integer                                    :: nsoil      ! Alias for soil texture class
   real                                       :: salloc     ! balive:bleaf ratio
   real                                       :: salloci    ! bleaf:balive ratio
   real                                       :: bleaf_max  ! maximum bleaf
   real                                       :: balive_max ! balive if on-allometry
   real                                       :: psi_layer  ! Water potential of this layer
   real                                       :: psi_wilt   ! Wilting point potential
   real                                       :: psi_crit   ! Critical point potential
   !---------------------------------------------------------------------------------------!


   cpatch => csite%patch(ipa)

   ipft = cpatch%pft(ico)

   !---------------------------------------------------------------------------------------!
   !     Here we decide how to compute the mean available water fraction.                  !
   !---------------------------------------------------------------------------------------!
   if (spot_phen) then
      !----- Use soil potential to determine phenology. -----------------------------------!
      cpatch%paw_avg(ico) = 0.0
      do k=cpatch%krdepth(ico),mzg
         nsoil = ntext_soil(k)
         
         psi_layer = soil(nsoil)%slpots                                                    &
                   / (csite%soil_water(k,ipa) / soil(nsoil)%slmsts) ** soil(nsoil)%slbs
         psi_wilt  = soil(nsoil)%slpots                                                    &
                   / (soil(nsoil)%soilwp / soil(nsoil)%slmsts) ** soil(nsoil)%slbs
         psi_crit  = soil(nsoil)%slpots                                                    &
                   / (soil(nsoil)%soilld / soil(nsoil)%slmsts) ** soil(nsoil)%slbs

         cpatch%paw_avg(ico) = cpatch%paw_avg(ico)                                         &
                             + max(0.0, (psi_layer - psi_wilt))                            &
                             * dslz(k) / (psi_crit  - psi_wilt)
      end do
      cpatch%paw_avg(ico) = - cpatch%paw_avg(ico) / slz(cpatch%krdepth(ico))
   else 
      !----- Use soil moisture (mass) to determine phenology. -----------------------------!
      cpatch%paw_avg(ico) = 0.0
      do k = cpatch%krdepth(ico), mzg
         nsoil = ntext_soil(k)
         cpatch%paw_avg(ico) = cpatch%paw_avg(ico)                                         &
                             + max(0.0, (csite%soil_water(k,ipa) - soil(nsoil)%soilwp))    &
                                      * dslz(k) / (soil(nsoil)%soilld - soil(nsoil)%soilwp)
      end do
      cpatch%paw_avg(ico) = - cpatch%paw_avg(ico) / slz(cpatch%krdepth(ico))
   end if

   !---------------------------------------------------------------------------------------!
   !    We make the elongation factor 1.0 when we are not solving the vegetation dynamics, !
   ! otherwise we assign the normal values.                                                !
   !---------------------------------------------------------------------------------------!
   select case (ivegt_dynamics)
   case (0)
      cpatch%elongf(ico) = 1.0

   case default
      select case (phenology(ipft))
      case (1)
         if (cpatch%paw_avg(ico) < 1.0) then
            cpatch%elongf(ico) = 0.0
         else
            cpatch%elongf(ico) = 1.0
         end if
      case (4)
         cpatch%elongf(ico)  = max(0.0,min(1.0,cpatch%paw_avg(ico)))
      case default
         cpatch%elongf(ico)  = 1.0
      end select

   end select
   !---------------------------------------------------------------------------------------!

   !----- Set phenology status according to the elongation factor. ------------------------!
   if (cpatch%elongf(ico) >= 1.0) then
      cpatch%phenology_status(ico) = 0
   elseif (cpatch%elongf(ico) > elongf_min) then
      cpatch%phenology_status(ico) = -1
   else
      cpatch%phenology_status(ico) = 2
      cpatch%elongf(ico)           = 0.
   end if
   !---------------------------------------------------------------------------------------!



   !----- Compute the biomass of living tissues. ------------------------------------------!
   salloc               = 1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico)
   salloci              = 1.0 / salloc
   bleaf_max            = dbh2bl(cpatch%dbh(ico),cpatch%pft(ico))
   balive_max           = bleaf_max * salloc
   select case (cpatch%phenology_status(ico))
   case (2)
      cpatch%bleaf(ico)  = 0.
      cpatch%elongf(ico) = 0.
   case default
      cpatch%bleaf(ico) = bleaf_max * cpatch%elongf(ico) 
   end select
   cpatch%broot(ico)    = balive_max * q(ipft)   * salloci
   cpatch%bsapwood(ico) = balive_max * qsw(ipft) * cpatch%hite(ico) * salloci
   cpatch%balive(ico)   = cpatch%bleaf(ico) + cpatch%broot(ico) + cpatch%bsapwood(ico)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Here we account for part of the carbon that didn't go to the leaves.  At this      !
   ! point  we will be nice to the plants and leave all the carbon that didn't go to       !
   ! leaves in the storage.  This gives some extra chance for the plant whilst it          !
   ! conserves the total carbon.                                                           !
   !---------------------------------------------------------------------------------------!
   cpatch%bstorage(ico) = max(0.0, bleaf_max - cpatch%bleaf(ico))
   !---------------------------------------------------------------------------------------!

   return
end subroutine pheninit_balive_bstorage
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

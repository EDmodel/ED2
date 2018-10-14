!==========================================================================================!
!==========================================================================================!
!     This module contains sub-routines that are useful to determine phenology.            !
!------------------------------------------------------------------------------------------!
module phenology_aux

   implicit none
   contains


   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine calculates phenology factors for prescribed phenology schemes.     !
   !---------------------------------------------------------------------------------------!
   subroutine prescribed_leaf_state(lat,imonth,iyear,doy,green_leaf_factor                 &
                                   ,leaf_aging_factor,phen_pars)

      use phenology_coms , only : iphenys1         & ! intent(in)
                                , iphenysf         & ! intent(in)
                                , iphenyf1         & ! intent(in)
                                , iphenyff         & ! intent(in)
                                , prescribed_phen  & ! intent(in)
                                , elongf_min       ! ! intent(in)
      use ed_max_dims    , only : n_pft            ! ! intent(in)
      use pft_coms       , only : phenology        ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(prescribed_phen) , intent(in)  :: phen_pars
      integer               , intent(in)  :: iyear
      integer               , intent(in)  :: doy
      real                  , intent(in)  :: lat
      integer               , intent(in)  :: imonth
      real, dimension(n_pft), intent(out) :: green_leaf_factor
      real, dimension(n_pft), intent(out) :: leaf_aging_factor
      !----- Local variables --------------------------------------------------------------!
      integer                             :: n_recycle_years
      integer                             :: my_year
      real                                :: elongf
      real                                :: delay
      real(kind=8)                        :: elonDen
      integer                             :: pft
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     This assumes dropping/flushing based on the day of year and hemisphere.        !
      ! + Northern Hemisphere: dropping between August 1 and December 31;                  !
      !                        flushing between January 1 and July 31.                     !
      ! + Southern Hemisphere: dropping between February 1 and July 31;                    !
      !                        flushing between August 1 and January 31.                   !
      !------------------------------------------------------------------------------------!
      if( (lat >= 0.0 .and. imonth <= 7) .or.                                              &
          (lat < 0.0  .and. (imonth > 7 .or. imonth == 1)) )then

         !----- Get the year. -------------------------------------------------------------!
         n_recycle_years = iphenysf - iphenys1 + 1

         if (iyear > iphenysf) then
            my_year = mod(iyear-iphenys1,n_recycle_years) + 1
         elseif (iyear < iphenys1) then
            my_year = n_recycle_years - mod(iphenysf-iyear,n_recycle_years)
         else
            my_year = iyear - iphenys1 + 1
         end if

         !---------------------------------------------------------------------------------!
         !      Calculate the factors.  Precalc denominator and limit rate in order to     !
         ! increase numerical stability (MCD 10/23/08).                                    !
         !---------------------------------------------------------------------------------!
         elonDen = real((phen_pars%flush_a(my_year) * real(doy)),kind=8)                   &
                 ** dble(max(phen_pars%flush_b(my_year),-100.))
         elonDen = 1.0d0 / (1.0d0 + elonDen)
         if(elonDen < 0.0001d0) then
            elongf = 0.0
         else
            elongf = sngl(elonDen)
         end if
         delay = elongf
      else
         !---------------------------------------------------------------------------------!
         !      Leaves turning color.  Get the year.                                       !
         !---------------------------------------------------------------------------------!
         n_recycle_years = iphenyff - iphenyf1 + 1
         if (iyear > iphenyff) then
            my_year = mod(iyear-iphenyf1,n_recycle_years) + 1
         elseif (iyear < iphenyf1) then
            my_year = n_recycle_years - mod(iphenyff-iyear,n_recycle_years)
         else
            my_year = iyear - iphenyf1 + 1
         end if

         !----- Calculate the factors. ----------------------------------------------------!
         elongf =  1.0                                                                     &
                /  (1.0 + (phen_pars%color_a(my_year) * real(doy))                         &
                ** phen_pars%color_b(my_year))
         delay  =  1.0                                                                     &
                /  (1.0 + (phen_pars%color_a(my_year) * real(doy) * 1.095)                 &
                ** phen_pars%color_b(my_year))
      end if

      if(elongf < elongf_min) elongf = 0.0

      !----- Load the values for each PFT. ------------------------------------------------!
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
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes the number of chill and warming days in a month.         !
   !  + Chill days  - number of days with average temperatures below 278.15 K;             !
   !  + Degree days - sum of daily average temperatures above 278.15 K.                    !
   !---------------------------------------------------------------------------------------!
   subroutine update_thermal_sums(month, cpoly, isi, lat)

      use ed_state_vars ,only : polygontype & ! structure
                              , sitetype    ! ! structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype) , target     :: cpoly
      integer           , intent(in) :: isi
      integer           , intent(in) :: month
      real              , intent(in) :: lat
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)    , pointer    :: csite
      integer                        :: ipa
      !------------------------------------------------------------------------------------!



      !----- Loop over patches. -----------------------------------------------------------!
      csite => cpoly%site(isi)

      do ipa = 1,csite%npatches

         !----- Minimum monthly temperature of the site. ----------------------------------!
         cpoly%min_monthly_temp(isi) = min(cpoly%min_monthly_temp(isi)                     &
                                          ,csite%avg_daily_temp(ipa))

         !---------------------------------------------------------------------------------!
         !    Warm day, so check whether it is growing season and update the degree day... !
         !---------------------------------------------------------------------------------!
         if (csite%avg_daily_temp(ipa) > 278.15) then
            !----- Update dgd only for growing season. ------------------------------------!
            if ((lat >= 0.0 .and. month <= 8) .or.                                         &
                (lat <  0.0 .and. (month <= 2 .or. month >= 7))) then
               csite%sum_dgd(ipa) = csite%sum_dgd(ipa) + (csite%avg_daily_temp(ipa)-278.15)
            !----- Warm day during dropping season, set degree sum to zero... -------------!
            else
               csite%sum_dgd(ipa) = 0.0
            end if
         !---- Cold day, check whether it is dropping season and update chilling days... --!
         elseif ((lat >= 0.0 .and. (month >= 11 .or. month <= 6)) .or.                     &
                 (lat <  0.0 .and.  month >= 5)                 ) then
            csite%sum_chd(ipa) = csite%sum_chd(ipa) + 1.0
         !---- Cold day, but not during dropping season, set chilling days to zero... -----!
         else
            csite%sum_chd(ipa) = 0.0
         end if
      end do

      return
   end subroutine update_thermal_sums
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine updates the turnover ratio and specific leaf area, taking into    !
   ! account the available radiation.                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine update_turnover(cpoly, isi)
      use ed_state_vars  , only : polygontype        & ! structure
                                , sitetype           & ! structure
                                , patchtype          ! ! structure
      use pft_coms       , only : phenology          & ! intent(in)
                                , sla_s0             & ! intent(in)
                                , sla_s1             & ! intent(in)
                                , leaf_turnover_rate ! ! intent(in)
      use phenology_coms , only : radint             & ! intent(in)
                                , radslp             & ! intent(in)
                                , turnamp_wgt        & ! intent(in)
                                , turnamp_min        & ! intent(in)
                                , turnamp_max        & ! intent(in)
                                , radto_min          & ! intent(in)
                                , radto_max          & ! intent(in)
                                , llspan_wgt         & ! intent(in)
                                , llspan_min         & ! intent(in)
                                , llspan_max         & ! intent(in)
                                , llspan_inf         & ! intent(in)
                                , vm0_wgt            & ! intent(in)
                                , vm0_tran           & ! intent(in)
                                , vm0_slope          & ! intent(in)
                                , vm0_amp            & ! intent(in)
                                , vm0_min            ! ! intent(in)
      use consts_coms    , only : day_sec            ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype) , target     :: cpoly
      integer           , intent(in) :: isi
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)    , pointer    :: csite
      type(patchtype)   , pointer    :: cpatch
      integer                        :: ipa
      integer                        :: ico
      integer                        :: ipft
      real                           :: turnover_now
      real                           :: turnamp_now
      real                           :: llspan_now
      real                           :: vm0_now
      !------------------------------------------------------------------------------------!



      !----- Loop over patches. -----------------------------------------------------------!
      csite => cpoly%site(isi)
      patchloop: do ipa = 1,csite%npatches

         cpatch => csite%patch(ipa)
         cohortloop: do ico = 1,cpatch%ncohorts

            ipft = cpatch%pft(ico)


            !------------------------------------------------------------------------------!
            !     We must check whether the light phenology is to be applied for this PFT. !
            !------------------------------------------------------------------------------!
            select case (phenology(ipft))
            case (3)
               !---------------------------------------------------------------------------!
               !      Find the target turnover rate amplitude (turnamp_now).               !
               !---------------------------------------------------------------------------!
               if (cpoly%rad_avg(isi) <= radto_min) then
                  turnamp_now = turnamp_min
               elseif (cpoly%rad_avg(isi) >= radto_max) then
                  turnamp_now = turnamp_max
               else
                  turnamp_now = radint + radslp * cpoly%rad_avg(isi)
               end if
               !---------------------------------------------------------------------------!


               !------ The actual turnover amplitude is based on a running average. -------!
               cpatch%turnover_amp(ico) = (1.0 - turnamp_wgt) * cpatch%turnover_amp(ico)   &
                                        +        turnamp_wgt  * turnamp_now
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Update target leaf lifespan.                                          !
               !---------------------------------------------------------------------------!
               if (leaf_turnover_rate(ipft) > 0.) then
                  llspan_now = 12.0 / (cpatch%turnover_amp(ico) * leaf_turnover_rate(ipft))
                  !----- Make sure the life span is bounded. ------------------------------!
                  if ( llspan_now < llspan_min) then
                      llspan_now = llspan_min
                  elseif (llspan_now > llspan_max) then
                      llspan_now = llspan_max
                  end if
                  !------------------------------------------------------------------------!
               else
                  !---- Nothing lasts forever, so impose a maximum life span. -------------!
                  llspan_now = llspan_inf
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     The actual leaf lifespan is the weighted average.                     !
               !---------------------------------------------------------------------------!
               cpatch%llspan(ico) = (1.0 - llspan_wgt) * cpatch%llspan(ico)                &
                                  +        llspan_wgt  * llspan_now
               !---------------------------------------------------------------------------!



               !----- Update the running average of the photosythetic capacity (Vm0). -----!
               vm0_now = vm0_amp / (1.0 + (cpatch%llspan(ico) / vm0_tran)**vm0_slope)      &
                       + vm0_min
               cpatch%vm_bar(ico) = (1.0 - vm0_wgt) * cpatch%vm_bar(ico)                   &
                                  + vm0_wgt * vm0_now
               !---------------------------------------------------------------------------!



               !----- Update the specific leaf area (SLA). --------------------------------!
               turnover_now    = cpatch%turnover_amp(ico) * leaf_turnover_rate(ipft)
               cpatch%sla(ico) = sla_s0(ipft) * turnover_now ** sla_s1(ipft)
            case default
               !----- The default is to keep these variables the same. --------------------!
               continue
               !---------------------------------------------------------------------------!
            end select
         end do cohortloop

      end do patchloop

      return
   end subroutine update_turnover
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This sub-routine will assign the initial potential available water and the       !
   ! initial leaf phenology.  This sub-routine should be called at the beginning of the    !
   ! simulation, except when this is a full history restart.  The initial running average  !
   ! is simply the the instantaneous soil moisture variable.  For plants other than the    !
   ! drought-deciduous, the potential available water is found but it doesn't control the  !
   ! phenology, so we assign fully flushed leaves.                                         !
   !---------------------------------------------------------------------------------------!
   subroutine first_phenology(cgrid)
      use ed_state_vars , only : edtype           & ! structure
                               , polygontype      & ! structure
                               , sitetype         & ! structure
                               , patchtype        ! ! structure
      use ed_therm_lib  , only : calc_veg_hcap    ! ! function
      use ed_max_dims   , only : n_pft            ! ! intent(in)
      use allometry     , only : ed_balive        & ! function
                               , ed_biomass       & ! function
                               , area_indices     ! ! subroutine
      use grid_coms     , only : nzg              ! ! intent(in)
      use therm_lib     , only : cmtl2uext        ! ! function
      use stable_cohorts, only : is_resolvable    ! ! sub-routine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)                   , target      :: cgrid       ! Current grid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)              , pointer     :: cpoly      ! Current polygon
      type(sitetype)                 , pointer     :: csite      ! Current site
      type(patchtype)                , pointer     :: cpatch     ! Current patch
      integer                                      :: ipy        ! Polygon counter
      integer                                      :: isi        ! Site counter
      integer                                      :: ipa        ! Patch counter
      integer                                      :: ico        ! Cohort counter
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Loop over all sites in this grid.                                              !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         nullify(cpoly)
         cpoly => cgrid%polygon(ipy)
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)


            !------------------------------------------------------------------------------!
            !     Loop over all patches and cohorts.                                       !
            !------------------------------------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)
               cohortloop: do ico=1,cpatch%ncohorts

                  !------------------------------------------------------------------------!
                  !    Find the initial guess for potential available water and elongation !
                  ! factor, then compute the equilibrium biomass of active tissues and     !
                  ! storage.                                                               !
                  !------------------------------------------------------------------------!
                  call pheninit_balive_bstorage(nzg,cpatch%pft(ico),cpatch%krdepth(ico)    &
                                               ,cpatch%hite(ico),cpatch%dbh(ico)           &
                                               ,csite%soil_water(:,ipa)                    &
                                               ,cpoly%ntext_soil(:,isi)                    &
                                               ,cpatch%paw_avg(ico),cpatch%elongf(ico)     &
                                               ,cpatch%phenology_status(ico)               &
                                               ,cpatch%bleaf(ico),cpatch%broot(ico)        &
                                               ,cpatch%bsapwooda(ico)                      &
                                               ,cpatch%bsapwoodb(ico),cpatch%bbarka(ico)   &
                                               ,cpatch%bbarkb(ico),cpatch%bstorage(ico)    &
                                               ,cpatch%cb(:,ico),cpatch%cb_lightmax(:,ico) &
                                               ,cpatch%cb_moistmax(:,ico)                  &
                                               ,cpatch%cb_mlmax(:,ico),cpatch%cbr_bar(ico))
                  !------------------------------------------------------------------------!


                  !----- Find Balive. -----------------------------------------------------!
                  cpatch%balive(ico) = ed_balive(cpatch, ico)
                  !------------------------------------------------------------------------!


                  !----- Find LAI, WAI, and CAI. ------------------------------------------!
                  call area_indices(cpatch, ico)
                  !------------------------------------------------------------------------!


                  !----- Find heat capacity and vegetation internal energy. ---------------!
                  call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdeada(ico)                  &
                                    ,cpatch%bsapwooda(ico),cpatch%bbarka(ico)              &
                                    ,cpatch%nplant(ico),cpatch%pft(ico)                    &
                                    ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
                  cpatch%leaf_energy(ico) = cmtl2uext(cpatch%leaf_hcap (ico)               &
                                                     ,cpatch%leaf_water(ico)               &
                                                     ,cpatch%leaf_temp (ico)               &
                                                     ,cpatch%leaf_fliq (ico))
                  cpatch%wood_energy(ico) = cmtl2uext(cpatch%wood_hcap (ico)               &
                                                     ,cpatch%wood_water(ico)               &
                                                     ,cpatch%wood_temp (ico)               &
                                                     ,cpatch%wood_fliq (ico))
                  call is_resolvable(csite,ipa,ico)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Update biomass (other derived properties such as timber and bark   !
                  ! thickness are not immediately affected by changes in bleaf).           !
                  !------------------------------------------------------------------------!
                  cpatch%agb(ico)     = ed_biomass(cpatch,ico)
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine first_phenology
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This sub-routine will assign the initial potential available water and the       !
   ! phenology that has been assigned, then find the biomass of active tissues and storage !
   ! that is in equilibrium with the initial soil moisture.  This sub-routine should be    !
   ! called whenever a new cohort is planted or recruited, or at the initial run (except   !
   ! history).  The initial running average is simply the the instantaneous soil moisture  !
   ! variable.  For plants other than the drought-deciduous, the potential available water !
   ! is found but it doesn't control the phenology, so we assign the biomass that matches  !
   ! the fully flushed leaves.                                                             !
   !---------------------------------------------------------------------------------------!
   subroutine pheninit_balive_bstorage(mzg,ipft,kroot,height,dbh,soil_water,ntext_soil     &
                                      ,paw_avg,elongf,phenology_status                     &
                                      ,bleaf,broot,bsapwooda,bsapwoodb,bbarka,bbarkb       &
                                      ,bstorage,cb,cb_lightmax,cb_moistmax,cb_mlmax        &
                                      ,cbr_bar)
      use soil_coms      , only : soil                & ! intent(in), look-up table
                                , slz                 & ! intent(in)
                                , slzt                & ! intent(in)
                                , dslz                & ! intent(in)
                                , matric_potential    ! ! function
      use phenology_coms , only : spot_phen           & ! intent(in)
                                , elongf_min          ! ! intent(in)
      use pft_coms       , only : phenology           & ! intent(in)
                                , q                   & ! intent(in)
                                , qsw                 & ! intent(in)
                                , qbark               & ! intent(in)
                                , agf_bs              & ! intent(in)
                                , f_bstorage_init     ! ! intent(in)
      use ed_max_dims    , only : n_pft               ! ! intent(in)
      use consts_coms    , only : almost_zero         ! ! intent(in)
      use physiology_coms, only : iddmort_scheme      ! ! intent(in)
      use allometry      , only : size2bl             & ! function
                                , h2crownbh           ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer                  , intent(in)  :: mzg               ! # of soil layers
      integer                  , intent(in)  :: ipft              ! PFT type
      integer                  , intent(in)  :: kroot             ! Level of rooting depth
      real                     , intent(in)  :: height            ! Height
      real                     , intent(in)  :: dbh               ! DBH
      integer, dimension(mzg)  , intent(in)  :: ntext_soil        ! Soil texture
      real   , dimension(mzg)  , intent(in)  :: soil_water        ! Soil water
      real                     , intent(out) :: paw_avg           ! Pot. available water
      real                     , intent(out) :: elongf            ! Elongation factor
      integer                  , intent(out) :: phenology_status  ! phenology Flag
      real                     , intent(out) :: bleaf             ! Leaf biomass
      real                     , intent(out) :: broot             ! Root biomass
      real                     , intent(out) :: bsapwooda         ! AG Sapwood biomass 
      real                     , intent(out) :: bsapwoodb         ! BG Sapwood biomass 
      real                     , intent(out) :: bbarka            ! AG Bark biomass 
      real                     , intent(out) :: bbarkb            ! BG Bark biomass 
      real                     , intent(out) :: bstorage          ! Storage biomass
      real   , dimension(13)   , intent(out) :: cb                ! Carbon balance
      real   , dimension(13)   , intent(out) :: cb_lightmax       ! Light-plenty CB
      real   , dimension(13)   , intent(out) :: cb_moistmax       ! Moisture-plenty CB
      real   , dimension(13)   , intent(out) :: cb_mlmax          ! Moisture/light-plenty CB
      real                     , intent(out) :: cbr_bar           ! Relative CB
      !----- Local variables --------------------------------------------------------------!
      integer                                :: k                 ! Layer counter
      integer                                :: nsoil             ! Soil texture class
      real                                   :: bleaf_max         ! maximum bleaf
      real                                   :: balive_max        ! balive if on-allometry
      real                                   :: psi_layer         ! Water pot. of this layer
      real                                   :: mcheight          ! Mid-crown height
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Here we decide how to compute the mean available water fraction.               !
      !------------------------------------------------------------------------------------!
      if (spot_phen) then
         !----- Use soil potential to determine phenology. --------------------------------!
         paw_avg = 0.0
         do k=kroot,mzg
            nsoil     = ntext_soil(k)
            mcheight  = 0.5 * ( height + h2crownbh(height,ipft) )

            psi_layer = slzt(k) - mcheight + matric_potential(nsoil,soil_water(k))
            paw_avg   = paw_avg + max(0.0, (psi_layer - soil(nsoil)%slpotwp)) * dslz(k)    &
                                / (soil(nsoil)%slpotld  - soil(nsoil)%slpotwp)
         end do
         paw_avg = paw_avg / abs(slz(kroot))
      else
         !----- Use soil moisture (mass) to determine phenology. --------------------------!
         paw_avg = 0.0
         do k = kroot, mzg
            nsoil   = ntext_soil(k)
            paw_avg = paw_avg + max(0.0, (soil_water(k) - soil(nsoil)%soilwp)) * dslz(k)   &
                              / (soil(nsoil)%soilld - soil(nsoil)%soilwp)
         end do
         paw_avg = paw_avg / abs(slz(kroot))
      end if

      !------------------------------------------------------------------------------------!
      !    Find the elongation factor according to the phenology of this PFT.              !
      !------------------------------------------------------------------------------------!
      select case (phenology(ipft))
      case (1)
         if (paw_avg < 1.0) then
            elongf = 0.0
         else
            elongf = 1.0
         end if
      case (3,4)
         elongf = max(0.0,min(1.0,paw_avg))
      case default
         elongf = 1.0
      end select
      !------------------------------------------------------------------------------------!

      !----- Set phenology status according to the elongation factor. ---------------------!
      if (elongf >= 1.0) then
         phenology_status = 0
      elseif (elongf > elongf_min) then
         phenology_status = 1
      else
         phenology_status = -2
         elongf           = 0.
      end if
      !------------------------------------------------------------------------------------!



      !----- Compute the biomass of living tissues. ---------------------------------------!
      bleaf_max  = size2bl(dbh,height,ipft)
      balive_max = bleaf_max * ( 1.0 + q(ipft) + ( qsw(ipft) + qbark(ipft) ) * height )
      bleaf      = bleaf_max * elongf
      broot      = bleaf_max * q    (ipft)
      bsapwooda  = bleaf_max * qsw  (ipft) * height * agf_bs(ipft)
      bsapwoodb  = bleaf_max * qsw  (ipft) * height * (1.0 - agf_bs(ipft))
      bbarka     = bleaf_max * qbark(ipft) * height * agf_bs(ipft)
      bbarkb     = bleaf_max * qbark(ipft) * height * (1.0 - agf_bs(ipft))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Here we account for part of the carbon that didn't go to the leaves.  At this   !
      ! point  we will be nice to the plants and leave all the carbon that didn't go to    !
      ! leaves in the storage.  This gives some extra chance for the plant whilst it       !
      ! conserves the total carbon.                                                        !
      !------------------------------------------------------------------------------------!
      bstorage = max(f_bstorage_init(ipft),almost_zero) * balive_max                       &
               + max(0.0, bleaf_max - bleaf)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Initialise the carbon balance.  For initial conditions, we always assume       !
      ! storage biomass for the previous months so the scale is correct (carbon balance is !
      ! given in kgC/pl).  The current month carbon balance must be initialised            !
      ! consistently with the iddmort_scheme set by the user.                              !
      !------------------------------------------------------------------------------------!
      cb         (1:12)  = bstorage
      cb_lightmax(1:12)  = bstorage
      cb_moistmax(1:12)  = bstorage
      cb_mlmax   (1:12)  = bstorage
      select case (iddmort_scheme)
      case (0)
         !------ Storage is not accounted. ------------------------------------------------!
         cb         (13) = 0.0
         cb_lightmax(13) = 0.0
         cb_moistmax(13) = 0.0
         cb_mlmax   (13) = 0.0
         !---------------------------------------------------------------------------------!
      case (1)
         cb         (13) = bstorage
         cb_lightmax(13) = bstorage
         cb_moistmax(13) = bstorage
         cb_mlmax   (13) = bstorage
      end select
      cbr_bar            = 1.0
      !------------------------------------------------------------------------------------!

      return
   end subroutine pheninit_balive_bstorage
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the length of daylight for a given latitude and day of     !
   ! year. The result is given in minutes.                                                 !
   !---------------------------------------------------------------------------------------!
   real function daylength(lat,doy)

      use consts_coms , only : pio180 & ! intent(in)
                             , twopi  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real    , intent(in) :: lat
      integer , intent(in) :: doy
      !----- Local variables --------------------------------------------------------------!
      real                 :: arg
      !------------------------------------------------------------------------------------!

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
   !=======================================================================================!
   !=======================================================================================!
end module phenology_aux
!==========================================================================================!
!==========================================================================================!

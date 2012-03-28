!==========================================================================================!
!==========================================================================================!
!    This subroutine controls the changes in leaf biomass due to phenology.                !
!------------------------------------------------------------------------------------------!
subroutine phenology_driver(cgrid, doy, month, tfact)
   use ed_state_vars  , only : edtype                & ! structure
                             , polygontype           & ! structure
                             , sitetype              ! ! structure
   use phenology_coms , only : iphen_scheme          ! ! intent(in)
   use ed_misc_coms   , only : current_time          ! ! intent(in)
   use phenology_aux  , only : prescribed_leaf_state & ! subroutine
                             , update_thermal_sums   & ! subroutine
                             , update_turnover       ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)      , target      :: cgrid
   integer           , intent(in)  :: doy
   integer           , intent(in)  :: month
   real              , intent(in)  :: tfact
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer     :: cpoly
   type(sitetype)    , pointer     :: csite
   integer                         :: ipy
   integer                         :: isi
   integer                         :: ipa
   !---------------------------------------------------------------------------------------!

   do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         !---------------------------------------------------------------------------------!
         !     Get the patch-level average daily temperature, which is needed for mortal-  !
         ! ity, recruitment and some phenology schemes.                                    !
         !---------------------------------------------------------------------------------!
         do ipa = 1,csite%npatches
            csite%avg_daily_temp(ipa) = csite%avg_daily_temp(ipa) * tfact
         end do
         
         select case (iphen_scheme)
         case (-1,0,2)
            !------------------------------------------------------------------------------!
            !     Default predictive scheme (Botta et al.) or the modified drought         !
            ! deciduous phenology for broadleaf PFTs.                                      !
            !------------------------------------------------------------------------------!
            call update_thermal_sums(month, cpoly, isi, cgrid%lat(ipy))
            call update_phenology(doy,cpoly,isi,cgrid%lat(ipy))
            
         case (1)
            !----- Use prescribed phenology. ----------------------------------------------!
            call prescribed_leaf_state(cgrid%lat(ipy), current_time%month                  &
                                      ,current_time%year, doy                              &
                                      ,cpoly%green_leaf_factor(:,isi)                      &
                                      ,cpoly%leaf_aging_factor(:,isi),cpoly%phen_pars(isi)) 
            call update_phenology(doy,cpoly,isi,cgrid%lat(ipy))


         case (3)
            !----- KIM light-controlled predictive phenology scheme. ----------------------!
            call update_thermal_sums(month, cpoly, isi, cgrid%lat(ipy))
            call update_turnover(cpoly,isi)
            call update_phenology(doy,cpoly,isi,cgrid%lat(ipy))
         end select
      end do
   end do

   return
end subroutine phenology_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will compute only the parts of the phenology that don t change the    !
! ecosystem state.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine phenology_driver_eq_0(cgrid, doy, month, tfact)
   use ed_state_vars  , only : edtype                & ! structure
                             , polygontype           & ! structure
                             , sitetype              ! ! structure
   use phenology_coms , only : iphen_scheme          ! ! intent(in)
   use ed_misc_coms   , only : current_time          ! ! intent(in)
   use phenology_aux  , only : prescribed_leaf_state & ! subroutine
                             , update_thermal_sums   & ! subroutine
                             , update_turnover       ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)      , target      :: cgrid
   integer           , intent(in)  :: doy
   integer           , intent(in)  :: month
   real              , intent(in)  :: tfact
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer     :: cpoly
   type(sitetype)    , pointer     :: csite
   integer                         :: ipy
   integer                         :: isi
   integer                         :: ipa
   !---------------------------------------------------------------------------------------!

   do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         !---------------------------------------------------------------------------------!
         !     Get the patch-level average daily temperature, which is needed for mortal-  !
         ! ity, recruitment and some phenology schemes.                                    !
         !---------------------------------------------------------------------------------!
         do ipa = 1,csite%npatches
            csite%avg_daily_temp(ipa) = csite%avg_daily_temp(ipa) * tfact
         end do
         
         select case (iphen_scheme)
         case (-1,0,2)
            !------------------------------------------------------------------------------!
            !     Default predictive scheme (Botta et al.) or the modified drought         !
            ! deciduous phenology for broadleaf PFTs.                                      !
            !------------------------------------------------------------------------------!
            call update_thermal_sums(month, cpoly, isi, cgrid%lat(ipy))
            call update_phenology_eq_0(doy,cpoly,isi,cgrid%lat(ipy))
            
         case (1)
            !----- Use prescribed phenology. ----------------------------------------------!
            call prescribed_leaf_state(cgrid%lat(ipy), current_time%month                  &
                                      ,current_time%year, doy                              &
                                      ,cpoly%green_leaf_factor(:,isi)                      &
                                      ,cpoly%leaf_aging_factor(:,isi),cpoly%phen_pars(isi)) 
            call update_phenology_eq_0(doy,cpoly,isi,cgrid%lat(ipy))

         case (3)
            !----- KIM light-controlled predictive phenology scheme. ----------------------!
            call update_thermal_sums(month, cpoly, isi, cgrid%lat(ipy))
            call update_turnover(cpoly,isi)
            call update_phenology_eq_0(doy,cpoly,isi,cgrid%lat(ipy))
         end select
      end do
   end do

   return
end subroutine phenology_driver_eq_0
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine update_phenology(doy, cpoly, isi, lat)

   use ed_state_vars  , only : polygontype              & ! structure
                             , sitetype                 & ! structure
                             , patchtype                ! ! structure
   use grid_coms      , only : nzg                      ! ! intent(in)
   use pft_coms       , only : phenology                & ! intent(in)
                             , c2n_leaf                 & ! intent(in)
                             , q                        & ! intent(in)
                             , qsw                      & ! intent(in)
                             , l2n_stem                 & ! intent(in)
                             , c2n_stem                 & ! intent(in)
                             , c2n_storage              ! ! intent(in)
   use decomp_coms    , only : f_labile                 ! ! intent(in)
   use phenology_coms , only : retained_carbon_fraction & ! intent(in)
                             , iphen_scheme             & ! intent(in)
                             , elongf_min               ! ! intent(in)
   use consts_coms    , only : t3ple                    & ! intent(in)
                             , cice                     & ! intent(in)
                             , cliq                     & ! intent(in)
                             , alli                     ! ! intent(in)
   use ed_therm_lib   , only : calc_veg_hcap            & ! function
                             , update_veg_energy_cweh   ! ! subroutine
   use ed_max_dims    , only : n_pft                    ! ! intent(in)
   use ed_misc_coms   , only : current_time             ! ! intent(in)
   use allometry      , only : area_indices             & ! subroutine
                             , ed_biomass               & ! function
                             , size2bl                  & ! function
                             , dbh2bl                   ! ! function
   use phenology_aux  , only : daylength                ! ! function
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype)        , target     :: cpoly
   real                     , intent(in) :: lat
   integer                  , intent(in) :: doy    ! Day of year (1=Jan 1, 365/366=Dec 31)
   integer                  , intent(in) :: isi
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)           , pointer    :: csite
   type(patchtype)          , pointer    :: cpatch
   integer                               :: ipa
   integer                               :: ico
   integer                               :: isoil_lev
   integer                               :: kroot
   integer                               :: ipft
   logical                               :: leaf_out_cold
   logical                               :: drop_cold
   real, dimension(nzg)                  :: theta
   real                                  :: daylight
   real                                  :: delta_bleaf
   real                                  :: bl_max
   real                                  :: old_leaf_hcap
   real                                  :: old_wood_hcap
   real                                  :: salloci
   !----- Variables used only for debugging purposes. -------------------------------------!
   logical                  , parameter  :: printphen=.false.
   logical, dimension(n_pft), save       :: first_time=.true.
   !---------------------------------------------------------------------------------------!


   !----- Level to evaluate the soil temperature. -----------------------------------------!
   isoil_lev = nzg 

   !----- Calculate daylength for this gridcell. ------------------------------------------!
   daylight = daylength(lat, doy) 

   !----- Loop over patches. --------------------------------------------------------------!
   csite => cpoly%site(isi)

   patchloop: do ipa = 1,csite%npatches
      cpatch => csite%patch(ipa)

      !----- Re-initialize litter inputs. -------------------------------------------------!
      csite%fsc_in(ipa) = 0.0
      csite%fsn_in(ipa) = 0.0
      csite%ssc_in(ipa) = 0.0
      csite%ssl_in(ipa) = 0.0

      !----- Determine what phenology thresholds have been crossed. -----------------------!
      call phenology_thresholds(daylight,csite%soil_tempk(isoil_lev,ipa)                   &
                               ,csite%soil_water(:,ipa),cpoly%ntext_soil(:,isi)            &
                               ,csite%sum_chd(ipa),csite%sum_dgd(ipa),drop_cold            &
                               ,leaf_out_cold,theta,cpoly%lsl(isi))

      cohortloop: do ico = 1,cpatch%ncohorts
         ipft    = cpatch%pft(ico)
         kroot   = cpatch%krdepth(ico)
         salloci = 1. / (1.0 + qsw(ipft) * cpatch%hite(ico) + q(ipft))
         
         !----- Initially, we assume all leaves stay. -------------------------------------!
         cpatch%leaf_drop(ico) = 0.0
         
         !----- Initially, we assume elongation factor to be 1. ---------------------------!
         cpatch%elongf(ico)    = 1.0

         !----- Find cohort-specific thresholds. ------------------------------------------!
         select case (iphen_scheme)
         case (1)
            !----- Get cohort-specific thresholds for prescribed phenology. ---------------!
            call assign_prescribed_phen(cpoly%green_leaf_factor(ipft,isi)                  &
                                       ,cpoly%leaf_aging_factor(ipft,isi),cpatch%dbh(ico)  &
                                       ,cpatch%hite(ico),ipft,drop_cold,leaf_out_cold      &
                                       ,bl_max)
         case default
            !----- Drop_cold is computed in phenology_thresholds for Botta scheme. --------!
            if (drop_cold) bl_max = 0.0
         end select

         !---------------------------------------------------------------------------------!
         !---------------------------------------------------------------------------------!
         !     Here we decide what to do depending on the phenology habit. There are five  !
         ! different types:                                                                !
         ! 0. Evergreen         - neither cold nor drought makes these plants to drop      !
         !                        their leaves;                                            !
         ! 1. Drought deciduous - these plants will drop all leaves when drought           !
         !                        conditions happen. By drought conditions we mean a time  !
         !                        when the available water drops below a threshold;        !
         ! 2. Cold deciduous    - these plants will drop their leaves when cold conditions !
         !                        happen.                                                  !
         ! 3. Light phenology   - these plants will control their leaf structure with      !
         !                        the light history (last 10 days).  On top of that, they  !
         !                        are also drought-deciduous, similar to phenology 4;      !
         ! 4. Drought deciduous - similar to one, but the threshold is compared against    !
         !                        a 10-day running average rather than the instantaneous   !
         !                        value.                                                   !
         !---------------------------------------------------------------------------------!
         select case (phenology(ipft))
         case (0)
            !------------------------------------------------------------------------------!
            !    Evergreen, there is nothing to be done here, move on.                     !
            !------------------------------------------------------------------------------!
            continue

         case (1)
            !------------------------------------------------------------------------------!
            !     Drought deciduous.  Now we must check whether the plants still have      !
            ! enough water or it is too dry, or if there is some drought relief so leaves  !
            ! can start to grow again.                                                     ! 
            !------------------------------------------------------------------------------!
            if (theta(kroot) < 1.0) then

               !----- It is time to drop leaves. ------------------------------------------!
               if (cpatch%phenology_status(ico) < 2) then
                  !----- Update litter pools. ---------------------------------------------!
                  cpatch%leaf_drop(ico) = (1.0 - retained_carbon_fraction)                 &
                                        * cpatch%lai(ico) / cpatch%sla(ico)                &
                                        / cpatch%nplant(ico)

                  csite%fsc_in(ipa) = csite%fsc_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * f_labile(ipft)
                  csite%fsn_in(ipa) = csite%fsn_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * f_labile(ipft) / c2n_leaf(ipft)
                  csite%ssc_in(ipa) = csite%ssc_in(ipa) + cpatch%leaf_drop(ico)            &
                                    * (1.0 - f_labile(ipft))
                  csite%ssl_in(ipa) = csite%ssl_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * (1.0 - f_labile(ipft)) * l2n_stem / c2n_stem(ipft)

                  
                  !----- Update plant carbon pools. ---------------------------------------!
                  cpatch%balive(ico)   = cpatch%balive(ico) - cpatch%bleaf(ico)

                  cpatch%bstorage(ico) = cpatch%bstorage(ico) + cpatch%bleaf(ico)          &
                                       * retained_carbon_fraction

                  !------------------------------------------------------------------------!
                  !     Contribution due to the fact that c2n_leaf and c2n_storage may be  !
                  ! different.                                                             !
                  !------------------------------------------------------------------------!
                  csite%fsn_in(ipa)            = csite%fsn_in(ipa)                         &
                                               + cpatch%bleaf(ico) * cpatch%nplant(ico)    &
                                               * retained_carbon_fraction                  &
                                               * (1.0 / c2n_leaf(ipft) - 1.0 / c2n_storage)
                  cpatch%bleaf(ico)            = 0.0
                  cpatch%phenology_status(ico) = 2
                  cpatch%elongf(ico)           = 0.
                  cpatch%cb(13,ico)            = cpatch%cb(13,ico)                         &
                                               - cpatch%leaf_drop(ico)
                  cpatch%cb_max(13,ico)        = cpatch%cb_max(13,ico)                     &
                                               - cpatch%leaf_drop(ico)
               end if
               
            elseif(theta(kroot) > 1.0 .and. cpatch%phenology_status(ico) == 2) then
               
               !----- It is time to flush.  Change phenology_status will update -----------!
               !----- carbon pools in growth_balive       ---------------------------------!
               cpatch%phenology_status(ico) = 1

            end if  ! critical moisture

         case (2)
            !------------------------------------------------------------------------------!
            !    Cold deciduous.  Here we must check two possibilities:                    !
            !                                                                              !
            ! 1. It is cold, and the plants have leaves, so we flag them with              !
            !    phenology_status=0 (leaves not growing) and the plants will start losing  !
            !    their leaves;                                                             !
            ! 2. The plant has no leaves, but the temperature and light conditions are     !
            !    okay again, and leaves can start growing.                                 !
            !------------------------------------------------------------------------------!
            if (cpatch%phenology_status(ico) < 2 .and. drop_cold) then
            
               if (cpoly%green_leaf_factor(ipft,isi) < elongf_min) then
                   bl_max = 0.0
               end if
               
               delta_bleaf = cpatch%bleaf(ico) - bl_max

               if (delta_bleaf > 0.0) then
                  !------------------------------------------------------------------------!
                  !    Phenology_status = 0 means that the plant has leaves, but they are  !
                  ! not growing (not necessarily because the leaves are fully flushed).    !
                  !------------------------------------------------------------------------!
                  cpatch%phenology_status(ico) = -1
                  cpatch%leaf_drop(ico) = (1.0 - retained_carbon_fraction) * delta_bleaf
                  csite%fsc_in(ipa) = csite%fsc_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * f_labile(ipft)
                  csite%fsn_in(ipa) = csite%fsn_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * f_labile(ipft) / c2n_leaf(ipft)
                  csite%ssc_in(ipa) = csite%ssc_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * (1.0 - f_labile(ipft))
                  csite%ssl_in(ipa) = csite%ssl_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * (1.0 - f_labile(ipft)) * l2n_stem / c2n_stem(ipft)
                  
                  !----- Adjust plant carbon pools. ---------------------------------------!
                  cpatch%balive(ico)   = cpatch%balive(ico) - delta_bleaf
                  cpatch%bstorage(ico) = cpatch%bstorage(ico)                              &
                                       + retained_carbon_fraction * delta_bleaf

                  !------------------------------------------------------------------------!
                  !     Contribution due to the fact that c2n_leaf and c2n_storage may be  !
                  ! different.                                                             !
                  !------------------------------------------------------------------------!
                  csite%fsn_in(ipa)     = csite%fsn_in(ipa)                                &
                                        + delta_bleaf * cpatch%nplant(ico)                 &
                                        * retained_carbon_fraction                         &
                                        * (1.0 / c2n_leaf(ipft) - 1.0/c2n_storage)
                  cpatch%bleaf(ico)     = bl_max
                  cpatch%cb(13,ico)     = cpatch%cb(13,ico) - cpatch%leaf_drop(ico)
                  cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico) - cpatch%leaf_drop(ico)
               end if

               !----- Set status flag. ----------------------------------------------------!
               if (bl_max == 0.0) then
                  cpatch%phenology_status(ico) = 2
                  cpatch%elongf(ico) = 0. 
               end if
               
            elseif (cpatch%phenology_status(ico) == 2 .and. leaf_out_cold) then
               !---------------------------------------------------------------------------!
               !      Update the phenology status (1 means that leaves are growing),       !
               !---------------------------------------------------------------------------!
               cpatch%phenology_status(ico) = 1
               ! NML -> now all growth is done in growth_balive
               !cpatch%bleaf(ico) = cpoly%green_leaf_factor(ipft,isi) * cpatch%balive(ico)  &
               !                  * salloci
            end if


         case (3,4) 
            !------------------------------------------------------------------------------!
            !    Drought deciduous.  Here we must check two possibilities:                 !
            !                                                                              !
            ! 1. The soil has been dry recently, and the plants have leaves, so we flag    !
            !    them with phenology_status=0 (leaves not growing) and the plants will     !
            !    start losing their leaves;                                                !
            ! 2. The plant has no leaves, but the soil has started to come back to more    !
            !    moist conditions. Given this situation, leaves can start growing again.   !
            !------------------------------------------------------------------------------!

            cpatch%elongf(ico) = max(0.0, min (1.0, cpatch%paw_avg(ico)))
            bl_max    = cpatch%elongf(ico) * size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)

            !----- In case it is too dry, drop all the leaves... --------------------------!
            if (cpatch%elongf(ico) < elongf_min) then
               bl_max = 0.0
            end if
            
            delta_bleaf = cpatch%bleaf(ico) - bl_max

            if (delta_bleaf > 0.0 .and. cpatch%phenology_status(ico) < 2) then
               cpatch%phenology_status(ico) = -1 
               cpatch%leaf_drop(ico) = (1.0 - retained_carbon_fraction) * delta_bleaf
               csite%fsc_in(ipa) = csite%fsc_in(ipa)                                       &
                                 + cpatch%nplant(ico) * cpatch%leaf_drop(ico)              &
                                 * f_labile(ipft)
               csite%fsn_in(ipa) = csite%fsn_in(ipa)                                       &
                                 + cpatch%nplant(ico) * cpatch%leaf_drop(ico)              &
                                 * f_labile(ipft) / c2n_leaf(ipft)
               csite%ssc_in(ipa) = csite%ssc_in(ipa)                                       &
                                 + cpatch%nplant(ico) * cpatch%leaf_drop(ico)              &
                                 * (1.0-f_labile(ipft))
               csite%ssl_in(ipa) = csite%ssl_in(ipa)                                       &
                                 + cpatch%nplant(ico) * cpatch%leaf_drop(ico)              &
                                 * (1.0 - f_labile(ipft)) * l2n_stem / c2n_stem(ipft)

               !----- Adjust plant carbon pools. ------------------------------------------!
               cpatch%balive(ico)   = cpatch%balive(ico) - delta_bleaf
               cpatch%bstorage(ico) = cpatch%bstorage(ico) + retained_carbon_fraction      &
                                    * delta_bleaf
               !---------------------------------------------------------------------------!
               !     Contribution due to the fact that c2n_leaf and c2n_storage may be     !
               ! different.                                                                !
               !---------------------------------------------------------------------------!
               csite%fsn_in(ipa)     = csite%fsn_in(ipa) + delta_bleaf*cpatch%nplant(ico)  &
                                     * retained_carbon_fraction                            &
                                     * (1.0 / c2n_leaf(ipft) - 1.0/c2n_storage)
               
               cpatch%bleaf(ico)     = bl_max
               
               if (cpatch%bleaf(ico) == 0.0) then
                  !----- No leaves. -------------------------------------------------------!
                  cpatch%phenology_status(ico) = 2
                  cpatch%elongf(ico) = 0.
               end if
               
               cpatch%cb(13,ico)     = cpatch%cb(13,ico)     - cpatch%leaf_drop(ico)
               cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico) - cpatch%leaf_drop(ico)
            !------ Becoming slightly moister again, start flushing the leaves. -----------!

            elseif (cpatch%elongf(ico)            > elongf_min .and.                       &
                    cpatch%phenology_status(ico) /= 0           ) then  
               !----- Not in allometry but growing, allocate carbon in growth_balive. -----!
               cpatch%phenology_status(ico) = 1                
            end if


         end select
         !---------------------------------------------------------------------------------!

         !----- Update LAI, WAI, and CAI accordingly. -------------------------------------!
         call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico),cpatch%bdead(ico)          &
                          ,cpatch%balive(ico),cpatch%dbh(ico),cpatch%hite(ico)             &
                          ,cpatch%pft(ico),cpatch%sla(ico),cpatch%lai(ico)                 &
                          ,cpatch%wai(ico),cpatch%crown_area(ico),cpatch%bsapwooda(ico))

         !----- Update above-ground biomass. ----------------------------------------------!
         cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico),cpatch%bleaf(ico)                  &
                                     ,cpatch%bsapwooda(ico),cpatch%pft(ico)) 

         !---------------------------------------------------------------------------------!
         !    The leaf biomass of the cohort has changed, update the vegetation energy -   !
         ! using a constant temperature assumption.                                        !
         !---------------------------------------------------------------------------------!
         old_leaf_hcap       = cpatch%leaf_hcap(ico)
         old_wood_hcap       = cpatch%wood_hcap(ico)
         call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdead(ico),cpatch%bsapwooda(ico)       &
                           ,cpatch%nplant(ico),cpatch%pft(ico)                             &
                           ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
         call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap)
         call is_resolvable(csite,ipa,ico,cpoly%green_leaf_factor(:,isi))

         !----- Printing some debugging stuff if the code is set for it. ------------------!
         if (printphen) then
            ipft=cpatch%pft(ico)
            if (first_time(ipft)) then
               first_time(ipft) = .false.
               write (unit=40+ipft,fmt='(a10,5(1x,a12))')                                  &
                  &'      TIME','       PATCH','      COHORT','      NPLANT'               &
                  &            ,'   LEAF_DROP','       THETA'
            end if
         
            write (unit=40+ipft,fmt='(2(i2.2,a1),i4.4,2(1x,i12),4(1x,es12.5))')            &
                 current_time%month,'/',current_time%date,'/',current_time%year,ipa,ico    &
                ,cpatch%nplant(ico),cpatch%leaf_drop(ico),theta(kroot)
         end if
      end do cohortloop
   end do patchloop
   return
end subroutine update_phenology
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will make a simpler version of the leaf phenology:  it will change   !
! bleaf according to cold or drought, but it won't check whether there is enough carbon    !
! to change it, it will simply change it to the target bleaf.                              !
!------------------------------------------------------------------------------------------!
subroutine update_phenology_eq_0(doy, cpoly, isi, lat)

   use ed_state_vars  , only : polygontype              & ! structure
                             , sitetype                 & ! structure
                             , patchtype                ! ! structure
   use grid_coms      , only : nzg                      ! ! intent(in)
   use pft_coms       , only : phenology                & ! intent(in)
                             , q                        & ! intent(in)
                             , qsw                      ! ! intent(in)
   use decomp_coms    , only : f_labile                 ! ! intent(in)
   use phenology_coms , only : retained_carbon_fraction & ! intent(in)
                             , iphen_scheme             & ! intent(in)
                             , elongf_min               ! ! intent(in)
   use ed_therm_lib   , only : calc_veg_hcap            & ! function
                             , update_veg_energy_cweh   ! ! subroutine
   use ed_max_dims    , only : n_pft                    ! ! intent(in)
   use ed_misc_coms   , only : current_time             ! ! intent(in)
   use allometry      , only : area_indices             & ! subroutine
                             , ed_biomass               & ! function
                             , dbh2bl                   ! ! function
   use phenology_aux  , only : daylength                ! ! function

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype)        , target     :: cpoly
   real                     , intent(in) :: lat
   integer                  , intent(in) :: doy    ! Day of year (1=Jan 1, 365/366=Dec 31)
   integer                  , intent(in) :: isi
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)           , pointer    :: csite
   type(patchtype)          , pointer    :: cpatch
   integer                               :: ipa
   integer                               :: ico
   integer                               :: isoil_lev
   integer                               :: kroot
   integer                               :: ipft
   logical                               :: leaf_out_cold
   logical                               :: drop_cold
   real, dimension(nzg)                  :: theta
   real                                  :: daylight
   real                                  :: delta_bleaf
   real                                  :: bleaf_new
   real                                  :: old_leaf_hcap
   real                                  :: old_wood_hcap
   real                                  :: salloci
   !---------------------------------------------------------------------------------------!


   !----- Level to evaluate the soil temperature. -----------------------------------------!
   isoil_lev = nzg 

   !----- Calculate daylength for this gridcell. ------------------------------------------!
   daylight = daylength(lat, doy) 

   !----- Loop over patches. --------------------------------------------------------------!
   csite => cpoly%site(isi)

   patchloop: do ipa = 1,csite%npatches
      cpatch => csite%patch(ipa)

      !----- Determine what phenology thresholds have been crossed. -----------------------!
      call phenology_thresholds(daylight,csite%soil_tempk(isoil_lev,ipa)                   &
                               ,csite%soil_water(:,ipa),cpoly%ntext_soil(:,isi)            &
                               ,csite%sum_chd(ipa),csite%sum_dgd(ipa),drop_cold            &
                               ,leaf_out_cold,theta,cpoly%lsl(isi))

      cohortloop: do ico = 1,cpatch%ncohorts
         ipft    = cpatch%pft(ico)
         kroot   = cpatch%krdepth(ico)
         salloci = 1. / (1.0 + qsw(ipft) * cpatch%hite(ico) + q(ipft))
         
         !----- Initially, we assume all leaves stay. -------------------------------------!
         cpatch%leaf_drop(ico) = 0.0
         
         !----- Initially, we assume elongation factor to be 1. ---------------------------!
         cpatch%elongf(ico)    = 1.0

         !----- Find cohort-specific thresholds. ------------------------------------------!
         select case (iphen_scheme)
         case (1)
            !----- Get cohort-specific thresholds for prescribed phenology. ---------------!
            call assign_prescribed_phen(cpoly%green_leaf_factor(ipft,isi)                  &
                                       ,cpoly%leaf_aging_factor(ipft,isi),cpatch%dbh(ico)  &
                                       ,cpatch%hite(ico),ipft,drop_cold,leaf_out_cold      &
                                       ,bleaf_new)
         case default
            !----- Drop_cold is computed in phenology_thresholds for Botta scheme. --------!
            if (drop_cold) bleaf_new = 0.0
         end select

         !---------------------------------------------------------------------------------!
         !---------------------------------------------------------------------------------!
         !     Here we decide what to do depending on the phenology habit. There are five  !
         ! different types:                                                                !
         ! 0. Evergreen         - neither cold nor drought makes these plants to drop      !
         !                        their leaves;                                            !
         ! 1. Drought deciduous - these plants will drop all leaves when drought           !
         !                        conditions happen. By drought conditions we mean a time  !
         !                        when the available water drops below a threshold;        !
         ! 2. Cold deciduous    - these plants will drop their leaves when cold conditions !
         !                        happen.                                                  !
         ! 3. Light phenology   - these plants will control their leaf structure with      !
         !                        the light history (last 10 days).  On top of that, they  !
         !                        are also drought-deciduous, similar to phenology 4;      !
         ! 4. Drought deciduous - similar to one, but the threshold is compared against    !
         !                        a 10-day running average rather than the instantaneous   !
         !                        value.                                                   !
         !---------------------------------------------------------------------------------!
         select case (phenology(ipft))
         case (0)
            !------------------------------------------------------------------------------!
            !    Evergreen, there is nothing to be done here, move on.                     !
            !------------------------------------------------------------------------------!
            continue

         case (1)
            !------------------------------------------------------------------------------!
            !     Drought deciduous.  Now we must check whether the plants still have      !
            ! enough water or it is too dry, or if there is some drought relief so leaves  !
            ! can start to grow again.                                                     ! 
            !------------------------------------------------------------------------------!
            if (theta(kroot) < 1.0) then

               !----- It is time to drop leaves. ------------------------------------------!
               if (cpatch%phenology_status(ico) < 2) then
                  !----- Find the leaf drop. ----------------------------------------------!
                  cpatch%leaf_drop(ico) = (1.0 - retained_carbon_fraction)                 &
                                        * cpatch%lai(ico) / cpatch%sla(ico)                &
                                        / cpatch%nplant(ico)

                  !----- Update plant carbon pools. ---------------------------------------!
                  cpatch%balive(ico)   = cpatch%balive(ico) - cpatch%bleaf(ico)
                  cpatch%bleaf(ico)            = 0.0
                  cpatch%phenology_status(ico) = 2
                  cpatch%elongf(ico)           = 0.
                  cpatch%cb(13,ico)            = cpatch%cb(13,ico)                         &
                                               - cpatch%leaf_drop(ico)
                  cpatch%cb_max(13,ico)        = cpatch%cb_max(13,ico)                     &
                                               - cpatch%leaf_drop(ico)
               end if
               
            elseif(theta(kroot) > 1.0 .and. cpatch%phenology_status(ico) == 2) then
               
               !---------------------------------------------------------------------------! 
               !    It is time to flush.  Change phenology_status here, and make all       !
               ! leaves appear all of a sudden.                                            !
               !---------------------------------------------------------------------------! 
               cpatch%phenology_status(ico) = 1
               cpatch%bleaf(ico)            = dbh2bl(cpatch%dbh(ico),ipft)
               cpatch%balive(ico)           = cpatch%balive(ico) + cpatch%bleaf(ico)
               cpatch%elongf(ico)           = 1.0
            end if  ! critical moisture

         case (2)
            !------------------------------------------------------------------------------!
            !    Cold deciduous.  Here we must check two possibilities:                    !
            !                                                                              !
            ! 1. It is cold, and the plants have leaves, so we flag them with              !
            !    phenology_status=0 (leaves not growing) and the plants will start losing  !
            !    their leaves;                                                             !
            ! 2. The plant has no leaves, but the temperature and light conditions are     !
            !    okay again, and leaves can start growing.                                 !
            !------------------------------------------------------------------------------!
            if (cpatch%phenology_status(ico) < 2 .and. drop_cold) then
            
               if (cpoly%green_leaf_factor(ipft,isi) < elongf_min) then
                  bleaf_new = 0.0
               else
                  bleaf_new = cpoly%green_leaf_factor(ipft,isi)                            &
                            * dbh2bl(cpatch%dbh(ico),ipft)
               end if
               
               delta_bleaf = cpatch%bleaf(ico) - bleaf_new

               if (delta_bleaf > 0.0) then
                  !------------------------------------------------------------------------!
                  !    Phenology_status = 0 means that the plant has leaves, but they are  !
                  ! not growing (not necessarily because the leaves are fully flushed).    !
                  !------------------------------------------------------------------------!
                  cpatch%phenology_status(ico) = -1
                  cpatch%leaf_drop(ico) = (1.0 - retained_carbon_fraction) * delta_bleaf
                  !------------------------------------------------------------------------!



                  !----- Adjust plant carbon pools. ---------------------------------------!
                  cpatch%balive(ico) = cpatch%balive(ico) - delta_bleaf
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Contribution due to the fact that c2n_leaf and c2n_storage may be  !
                  ! different.                                                             !
                  !------------------------------------------------------------------------!
                  cpatch%bleaf(ico)     = bleaf_new
                  cpatch%cb(13,ico)     = cpatch%cb(13,ico) - cpatch%leaf_drop(ico)
                  cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico) - cpatch%leaf_drop(ico)
                  !------------------------------------------------------------------------!
               end if

               !----- Set status flag. ----------------------------------------------------!
               if (bleaf_new == 0.0) then
                  cpatch%phenology_status(ico) = 2
               end if
               !---------------------------------------------------------------------------!


            elseif (cpatch%phenology_status(ico) == 2 .and. leaf_out_cold) then
               !---------------------------------------------------------------------------!
               !      Update the phenology status (1 means that leaves are growing),       !
               !---------------------------------------------------------------------------!
               cpatch%phenology_status(ico) = 1
               cpatch%bleaf(ico)            = cpoly%green_leaf_factor(ipft,isi)            &
                                            * dbh2bl(cpatch%dbh(ico),ipft)
               cpatch%balive(ico)           = cpatch%balive(ico) + cpatch%bleaf(ico)
               !---------------------------------------------------------------------------!
            end if


         case (3,4) 
            !------------------------------------------------------------------------------!
            !    Drought deciduous.  Here we must check two possibilities:                 !
            !                                                                              !
            ! 1. The soil has been dry recently, and the plants have leaves, so we flag    !
            !    them with phenology_status=0 (leaves not growing) and the plants will     !
            !    start losing their leaves;                                                !
            ! 2. The plant has no leaves, but the soil has started to come back to more    !
            !    moist conditions. Given this situation, leaves can start growing again.   !
            !------------------------------------------------------------------------------!
            cpatch%elongf(ico) = max(0.0, min (1.0, cpatch%paw_avg(ico)))

            !------------------------------------------------------------------------------!
            !     In case it is too dry, drop all the leaves, otherwise scale the leaf     !
            ! biomass with elongation factor.                                              !
            !------------------------------------------------------------------------------!
            if (cpatch%elongf(ico) < elongf_min) then
               bleaf_new = 0.0
            else
               bleaf_new = cpatch%elongf(ico) * dbh2bl(cpatch%dbh(ico),ipft)
            end if
            
            delta_bleaf = cpatch%bleaf(ico) - bleaf_new

            if (delta_bleaf > 0.0 .and. cpatch%phenology_status(ico) < 2) then
               cpatch%phenology_status(ico) = -1 
               cpatch%leaf_drop(ico) = (1.0 - retained_carbon_fraction) * delta_bleaf

               !----- Adjust plant carbon pools. ------------------------------------------!
               cpatch%balive(ico)   = cpatch%balive(ico) - delta_bleaf

               !---------------------------------------------------------------------------!
               !     Contribution due to the fact that c2n_leaf and c2n_storage may be     !
               ! different.                                                                !
               !---------------------------------------------------------------------------!               
               cpatch%bleaf(ico)     = bleaf_new
               
               if (cpatch%bleaf(ico) == 0.0) then
                  !----- No leaves. -------------------------------------------------------!
                  cpatch%phenology_status(ico) = 2
                  cpatch%elongf(ico)           = 0.
               end if
               
               cpatch%cb(13,ico)     = cpatch%cb(13,ico)     - cpatch%leaf_drop(ico)
               cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico) - cpatch%leaf_drop(ico)
            !------ Becoming slightly moister again, start flushing the leaves. -----------!

            elseif (cpatch%elongf(ico)            > elongf_min .and.                       &
                    cpatch%phenology_status(ico) /= 0           ) then  
               !----- Not in allometry but growing, allocate carbon in growth_balive. -----!
               if (cpatch%elongf(ico) == 1.0) then
                  cpatch%phenology_status(ico) = 0
               else
                  cpatch%phenology_status(ico) = 1
               end if
               cpatch%bleaf(ico)            = bleaf_new
               cpatch%balive(ico)           = cpatch%balive(ico) - delta_bleaf
            end if


         end select
         !---------------------------------------------------------------------------------!
         !----- Update LAI, WAI, and CAI accordingly. -------------------------------------!
         call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico),cpatch%bdead(ico)          &
                          ,cpatch%balive(ico),cpatch%dbh(ico),cpatch%hite(ico)             &
                          ,cpatch%pft(ico),cpatch%sla(ico),cpatch%lai(ico)                 &
                          ,cpatch%wai(ico),cpatch%crown_area(ico),cpatch%bsapwooda(ico))

         !----- Update above-ground biomass. ----------------------------------------------!
         cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico),cpatch%bleaf(ico)                  &
                                     ,cpatch%bsapwooda(ico),cpatch%pft(ico)) 

         !---------------------------------------------------------------------------------!
         !    The leaf biomass of the cohort has changed, update the vegetation energy -   !
         ! using a constant temperature assumption.                                        !
         !---------------------------------------------------------------------------------!
         old_leaf_hcap       = cpatch%leaf_hcap(ico)
         old_wood_hcap       = cpatch%wood_hcap(ico)
         call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdead(ico),cpatch%bsapwooda(ico)      &
                           ,cpatch%nplant(ico),cpatch%pft(ico)                             &
                           ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
         call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap)
         call is_resolvable(csite,ipa,ico,cpoly%green_leaf_factor(:,isi))
      end do cohortloop
   end do patchloop
   return
end subroutine update_phenology_eq_0
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine establishes whether it´s time to drop leaves or start flushing them  !
! for cold deciduous or temperate drought deciduous.                                       !
! MLO. Shouldn´t we have a similar criterion for both tropical and temperate, based on a   !
!      long term dry condition?                                                            !
!------------------------------------------------------------------------------------------!
subroutine phenology_thresholds(daylight,soil_temp,soil_water,soil_class,sum_chd,sum_dgd   &
                               ,drop_cold,leaf_out_cold,theta,lsl)
   use grid_coms     , only : nzg          ! ! intent(in)
   use soil_coms     , only : soil         & ! intent(in)
                            , dslz         & ! intent(in)
                            , slz          ! ! intent(in)
   use phenology_coms, only : dl_tr        & ! intent(in)
                            , st_tr1       & ! intent(in)
                            , st_tr2       & ! intent(in)
                            , phen_a       & ! intent(in)
                            , phen_b       & ! intent(in)
                            , phen_c       & ! intent(in)
                            , iphen_scheme ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, dimension(nzg), intent(in)    :: soil_class    ! Soil class
   integer                , intent(in)    :: lsl           ! Lowest soil level
   real                   , intent(in)    :: daylight      ! Daytime Length
   real                   , intent(in)    :: soil_temp     ! 
   real   , dimension(nzg), intent(in)    :: soil_water    !
   real                   , intent(inout) :: sum_dgd       !
   real                   , intent(inout) :: sum_chd       !
   logical                , intent(out)   :: drop_cold     !
   logical                , intent(out)   :: leaf_out_cold !
   real   , dimension(nzg), intent(out)   :: theta         !
   !----- Local variables -----------------------------------------------------------------!
   real                                   :: gdd_threshold
   integer                                :: k1
   integer                                :: k2
   integer                                :: topsoil
   integer                                :: nsoil
   !---------------------------------------------------------------------------------------!

   !----- Initialize variables. -----------------------------------------------------------!
   drop_cold     = .false.
   leaf_out_cold = .false.
   theta(1:nzg)  = 0.0

   select case (iphen_scheme)
   !---------------------------------------------------------------------------------------!
   !     In case we use the prescribed phenology, we skip this and fill with the prescrib- !
   ! ed phenology.                                                                         !
   !---------------------------------------------------------------------------------------!
   case (1)
      drop_cold     = .false.
      leaf_out_cold = .false.

   !---------------------------------------------------------------------------------------!
   !     Otherwise, we are solving the phenology, determine whether or not it´s time to    !
   ! drop leaves or start flushing.                                                        !
   !---------------------------------------------------------------------------------------!
   case default
      !----- Too cold or too dark, time to shed leaves... ---------------------------------!
      drop_cold = (daylight <= dl_tr .and. soil_temp < st_tr1) .or.  soil_temp < st_tr2

      !----- Getting warmer again, time for leaves to come out again... -------------------!
      gdd_threshold = phen_a + phen_b * exp(phen_c * sum_chd)
      leaf_out_cold = sum_dgd >= gdd_threshold
   end select

   !---------------------------------------------------------------------------------------!
   !     Calculate average theta for drought deciduous PFTs.  The different  k1's          !
   ! represent different rooting depths.                                                   !
   !---------------------------------------------------------------------------------------!  
   theta(1:nzg) = 0.0
   topsoil      = soil_class(nzg)
   do k1 = lsl, nzg
      do k2 = k1,nzg
         nsoil     = soil_class(k2)
         theta(k1) = theta(k1)                                                             &
                   + (soil_water(k2)     - soil(nsoil)%soilwp) * dslz(k2)                  &
                   / (soil(nsoil)%soilld - soil(nsoil)%soilwp)
      end do
      theta(k1) = - theta(k1) / slz(k1)
   end do

   return
end subroutine phenology_thresholds
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine assigns the prescribed phenology.                                     !
!------------------------------------------------------------------------------------------!
subroutine assign_prescribed_phen(green_leaf_factor,leaf_aging_factor,dbh,height,pft       &
                                 ,drop_cold,leaf_out_cold,bl_max)
   use allometry     , only : dbh2bl
   use phenology_coms, only : elongf_min
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   logical, intent(out) :: drop_cold
   logical, intent(out) :: leaf_out_cold
   real   , intent(out) :: bl_max
   real   , intent(in)  :: green_leaf_factor
   real   , intent(in)  :: leaf_aging_factor
   real   , intent(in)  :: dbh
   real   , intent(in)  :: height
   integer, intent(in)  :: pft
   !---------------------------------------------------------------------------------------!

   drop_cold     = green_leaf_factor /= leaf_aging_factor
   leaf_out_cold = green_leaf_factor > elongf_min .and. (.not. drop_cold)
   bl_max        = green_leaf_factor * dbh2bl(dbh, pft)

   return
end subroutine assign_prescribed_phen
!==========================================================================================!
!==========================================================================================!

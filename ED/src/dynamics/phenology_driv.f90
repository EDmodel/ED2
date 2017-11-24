!==========================================================================================!
!==========================================================================================!
!     Module that controls changes in leaf biomass due to phenology.                       !
!------------------------------------------------------------------------------------------!
module phenology_driv
   contains

   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine controls the changes in leaf biomass due to phenology.             !
   !---------------------------------------------------------------------------------------!
   subroutine phenology_driver(cgrid, doy, month, tfact,veget_dyn_on)
      use ed_state_vars  , only : edtype                & ! structure
                                , polygontype           & ! structure
                                , sitetype              ! ! structure
      use phenology_coms , only : iphen_scheme          ! ! intent(in)
      use ed_misc_coms   , only : current_time          ! ! intent(in)
      use phenology_aux  , only : prescribed_leaf_state & ! subroutine
                                , update_thermal_sums   & ! subroutine
                                , update_turnover       ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)      , target      :: cgrid
      integer           , intent(in)  :: doy
      integer           , intent(in)  :: month
      real              , intent(in)  :: tfact
      logical           , intent(in)  :: veget_dyn_on
      !----- Local variables --------------------------------------------------------------!
      type(polygontype) , pointer     :: cpoly
      type(sitetype)    , pointer     :: csite
      integer                         :: ipy
      integer                         :: isi
      integer                         :: ipa
      !------------------------------------------------------------------------------------!

      do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !------------------------------------------------------------------------------!
            !     Get the patch-level average daily temperature, which is needed for       !
            ! mortality, recruitment and some phenology schemes.                           !
            !------------------------------------------------------------------------------!
            do ipa = 1,csite%npatches
               csite%avg_daily_temp(ipa) = csite%avg_daily_temp(ipa) * tfact
            end do
            
            select case (iphen_scheme)
            case (-1,0,2)
               !---------------------------------------------------------------------------!
               !     Default predictive scheme (Botta et al.) or the modified drought      !
               ! deciduous phenology for broadleaf PFTs.                                   !
               !---------------------------------------------------------------------------!
               call update_thermal_sums(month, cpoly, isi, cgrid%lat(ipy))
               call update_phenology(doy,cpoly,isi,cgrid%lat(ipy),veget_dyn_on)
               
            case (1)
               !----- Use prescribed phenology. -------------------------------------------!
               call prescribed_leaf_state(cgrid%lat(ipy), current_time%month               &
                                         ,current_time%year, doy                           &
                                         ,cpoly%green_leaf_factor(:,isi)                   &
                                         ,cpoly%leaf_aging_factor(:,isi)                   &
                                         ,cpoly%phen_pars(isi) ) 
               call update_phenology(doy,cpoly,isi,cgrid%lat(ipy),veget_dyn_on)


            case (3)
               !----- KIM light-controlled predictive phenology scheme. -------------------!
               call update_thermal_sums(month, cpoly, isi, cgrid%lat(ipy))
               call update_turnover(cpoly,isi)
               call update_phenology(doy,cpoly,isi,cgrid%lat(ipy),veget_dyn_on)
            end select
         end do
      end do

      return
   end subroutine phenology_driver
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine update_phenology(doy, cpoly, isi, lat,veget_dyn_on)
      use ed_state_vars  , only : polygontype              & ! structure
                                , sitetype                 & ! structure
                                , patchtype                ! ! structure
      use grid_coms      , only : nzg                      ! ! intent(in)
      use pft_coms       , only : phenology                & ! intent(in)
                                , c2n_leaf                 & ! intent(in)
                                , l2n_stem                 & ! intent(in)
                                , c2n_stem                 & ! intent(in)
                                , c2n_storage              & ! intent(in)
                                , f_labile_leaf            ! ! intent(in)
      use phenology_coms , only : retained_carbon_fraction & ! intent(in)
                                , iphen_scheme             & ! intent(in)
                                , elongf_min               & ! intent(in)
                                , elongf_flush             ! ! intent(in)
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
                                , size2bl                  ! ! function
      use phenology_aux  , only : daylength                ! ! function
      use stable_cohorts , only : is_resolvable            ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype)        , target     :: cpoly
      real                     , intent(in) :: lat
      integer                  , intent(in) :: doy  ! Day of year (1=Jan 1, 365/366=Dec 31)
      integer                  , intent(in) :: isi
      logical                  , intent(in) :: veget_dyn_on
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)           , pointer    :: csite
      type(patchtype)          , pointer    :: cpatch
      integer                               :: ipa
      integer                               :: ico
      integer                               :: isoil_lev
      integer                               :: kroot
      integer                               :: ipft
      logical                               :: leaf_out_cold
      logical                               :: drop_cold
      real                                  :: daylight
      real                                  :: delta_bleaf
      real                                  :: bl_max
      real                                  :: old_leaf_hcap
      real                                  :: old_wood_hcap
      real                                  :: elongf_try
      real                                  :: elongf_grow
      !----- Variables used only for debugging purposes. ----------------------------------!
      logical                  , parameter  :: printphen=.false.
      logical, dimension(n_pft), save       :: first_time=.true.
      !------------------------------------------------------------------------------------!


      !----- Level to evaluate the soil temperature. --------------------------------------!
      isoil_lev = nzg 

      !----- Calculate daylength for this gridcell. ---------------------------------------!
      daylight = daylength(lat, doy) 

      !----- Loop over patches. -----------------------------------------------------------!
      csite => cpoly%site(isi)

      patchloop: do ipa = 1,csite%npatches
         cpatch => csite%patch(ipa)

         !----- Re-initialize litter inputs. ----------------------------------------------!
         csite%fsc_in(ipa) = 0.0
         csite%fsn_in(ipa) = 0.0
         csite%ssc_in(ipa) = 0.0
         csite%ssl_in(ipa) = 0.0

         !----- Determine what phenology thresholds have been crossed. --------------------!
         call phenology_thresholds(daylight,csite%soil_tempk(isoil_lev,ipa)                &
                                  ,csite%sum_chd(ipa),csite%sum_dgd(ipa)                   &
                                  ,drop_cold,leaf_out_cold)

         cohortloop: do ico = 1,cpatch%ncohorts
            ipft    = cpatch%pft(ico)
            kroot   = cpatch%krdepth(ico)
            
            !----- Initially, we assume all leaves stay. ----------------------------------!
            cpatch%leaf_drop(ico) = 0.0

            !----- Find cohort-specific thresholds. ---------------------------------------!
            select case (iphen_scheme)
            case (1)
               !----- Get cohort-specific thresholds for prescribed phenology. ------------!
               call assign_prescribed_phen(cpoly%green_leaf_factor(ipft,isi)               &
                                          ,cpoly%leaf_aging_factor(ipft,isi)               &
                                          ,cpatch%dbh(ico),cpatch%hite(ico),ipft           &
                                          ,drop_cold,leaf_out_cold,bl_max)
            case default
               !----- Drop_cold is computed in phenology_thresholds for Botta scheme. -----!
               if (drop_cold) bl_max = 0.0
            end select

            !------------------------------------------------------------------------------!
            !------------------------------------------------------------------------------!
            !     Here we decide what to do depending on the phenology habit. There are    !
            ! five different types:                                                        !
            ! 0. Evergreen         - neither cold nor drought makes these plants to drop   !
            !                        their leaves;                                         !
            ! 1. Drought deciduous - these plants will drop all leaves when drought        !
            !                        conditions happen. By drought conditions we mean a    !
            !                        time when the available water drops below a threshold !
            ! 2. Cold deciduous    - these plants will drop their leaves when cold         !
            !                        conditions happen.                                    !
            ! 3. Light phenology   - these plants will control their leaf structure with   !
            !                        the light history (last 10 days).  They are also      !
            !                        drought-deciduous, similar to phenology 4;            !
            ! 4. Drought deciduous - similar to one, but the threshold is compared against !
            !                        a 10-day running average rather than the instant-     !
            !                        aneous value.                                         !
            !------------------------------------------------------------------------------!
            select case (phenology(ipft))
            case (0)
               !---------------------------------------------------------------------------!
               !    Evergreen, there is nothing to be done here except to assign maximum   !
               ! elongation factor.                                                        !
               !---------------------------------------------------------------------------!
               cpatch%elongf(ico)    = 1.0
               !---------------------------------------------------------------------------!

            case (1)
               !---------------------------------------------------------------------------!
               !     Drought deciduous.  Now we must check whether the plants still have   !
               ! enough water or it is too dry, or if there is some drought relief so      !
               ! leaves can start to grow again.                                           !
               !---------------------------------------------------------------------------!
               !----- This is the first guess for the new elongation factor. --------------!
               elongf_try  = max(0.0, min (1.0, cpatch%paw_avg(ico)))
               !---------------------------------------------------------------------------!



               if (elongf_try < 1.0 .and. cpatch%phenology_status(ico) /= -2) then
                  !----- It is time to drop leaves.  Drop all leaves. ---------------------!
                  cpatch%leaf_drop(ico) = (1.0 - retained_carbon_fraction)                 &
                                        * cpatch%bleaf(ico)
                  !------------------------------------------------------------------------!



                  !----- Update plant carbon pools. ---------------------------------------!
                  cpatch%balive  (ico) = cpatch%balive(ico) - cpatch%bleaf(ico)
                  cpatch%bleaf   (ico) = 0.0
                  cpatch%elongf  (ico) = 0.0
                  cpatch%phenology_status(ico) = -2
                  !------------------------------------------------------------------------!

                  !----- Update storage only if vegetation dynamics is on. ----------------!
                  if (veget_dyn_on) then
                     cpatch%bstorage(ico) = cpatch%bstorage(ico)                           &
                                          + cpatch%bleaf(ico) * retained_carbon_fraction
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Send the lost leaves to soil carbon and nitrogen pools.            !
                  !------------------------------------------------------------------------!
                  csite%fsc_in(ipa) = csite%fsc_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * f_labile_leaf(ipft)
                  csite%fsn_in(ipa) = csite%fsn_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * f_labile_leaf(ipft) / c2n_leaf(ipft)
                  csite%ssc_in(ipa) = csite%ssc_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * (1.0 - f_labile_leaf(ipft))
                  csite%ssl_in(ipa) = csite%ssl_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * (1.0 - f_labile_leaf(ipft))                          &
                                    * l2n_stem / c2n_stem(ipft)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Contribution due to the fact that c2n_leaf and c2n_storage may be  !
                  ! different.                                                             !
                  !------------------------------------------------------------------------!
                  csite%fsn_in(ipa) = csite%fsn_in(ipa)                                    &
                                    + cpatch%leaf_drop(ico) * cpatch%nplant(ico)           &
                                    * (1.0 / c2n_leaf(ipft) - 1.0 / c2n_storage)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Deduct the leaf drop from the carbon balance.                     !
                  !------------------------------------------------------------------------!
                  cpatch%cb          (13,ico)  = cpatch%cb          (13,ico)               &
                                               - cpatch%leaf_drop      (ico)
                  cpatch%cb_lightmax (13,ico)  = cpatch%cb_lightmax (13,ico)               &
                                               - cpatch%leaf_drop      (ico)
                  cpatch%cb_moistmax (13,ico)  = cpatch%cb_moistmax (13,ico)               &
                                               - cpatch%leaf_drop      (ico)
                  cpatch%cb_mlmax (13,ico)     = cpatch%cb_mlmax (13,ico)                  &
                                               - cpatch%leaf_drop (ico)
                  


                  !------------------------------------------------------------------------!

               elseif(elongf_try > 1.0 .and. cpatch%phenology_status(ico) == -2) then
                  !------------------------------------------------------------------------!
                  !      It is time to flush.  Change phenology_status will update carbon  !
                  ! pools in growth_balive.                                                !
                  !------------------------------------------------------------------------!
                  if (veget_dyn_on) then
                     cpatch%phenology_status(ico) = 1
                     cpatch%elongf          (ico) = 1.0
                  else
                     !---------------------------------------------------------------------!
                     !    When vegetation dynamics is off, instantaneously update tissues  !
                     !---------------------------------------------------------------------!
                     cpatch%bleaf(ico)  = size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
                     cpatch%balive(ico) = cpatch%balive(ico) + cpatch%bleaf(ico)
                     cpatch%phenology_status(ico) = 0
                     cpatch%elongf          (ico) = 1.0
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!
               end if  ! critical moisture

            case (2)
               !---------------------------------------------------------------------------!
               !    Cold deciduous.  Here we must check two possibilities:                 !
               !                                                                           !
               ! 1. It is cold, and the plants have leaves, so we flag them with           !
               !    phenology_status=0 (leaves not growing) and the plants will start      !
               !    losing their leaves;                                                   !
               ! 2. The plant has no leaves, but the temperature and light conditions are  !
               !    okay again, and leaves can start growing.                              !
               !---------------------------------------------------------------------------!
               if (cpatch%phenology_status(ico) /= -2 .and. drop_cold) then
                  if (cpoly%green_leaf_factor(ipft,isi) < elongf_min) then
                      bl_max = 0.0
                  end if
                  delta_bleaf = cpatch%bleaf(ico) - bl_max

                  if (delta_bleaf > 0.0) then
                     !---------------------------------------------------------------------!
                     !    Phenology_status = 0 means that the plant has leaves, but they   !
                     ! are not growing (not necessarily because the leaves are fully       !
                     ! flushed).                                                           !
                     !---------------------------------------------------------------------!
                     cpatch%phenology_status(ico) = -1
                     cpatch%leaf_drop(ico) = (1.0 - retained_carbon_fraction) * delta_bleaf
                     csite%fsc_in(ipa) = csite%fsc_in(ipa)                                 &
                                       + cpatch%nplant(ico) * cpatch%leaf_drop(ico)        &
                                       * f_labile_leaf(ipft)
                     csite%fsn_in(ipa) = csite%fsn_in(ipa)                                 &
                                       + cpatch%nplant(ico) * cpatch%leaf_drop(ico)        &
                                       * f_labile_leaf(ipft) / c2n_leaf(ipft)
                     csite%ssc_in(ipa) = csite%ssc_in(ipa)                                 &
                                       + cpatch%nplant(ico) * cpatch%leaf_drop(ico)        &
                                       * (1.0 - f_labile_leaf(ipft))
                     csite%ssl_in(ipa) = csite%ssl_in(ipa)                                 &
                                       + cpatch%nplant(ico) * cpatch%leaf_drop(ico)        &
                                       * (1.0 - f_labile_leaf(ipft))                       &
                                       * l2n_stem / c2n_stem(ipft)
                     
                     !----- Adjust plant carbon pools. ------------------------------------!
                     cpatch%balive(ico)   = cpatch%balive(ico) - delta_bleaf
                     if (veget_dyn_on) then
                        cpatch%bstorage(ico) = cpatch%bstorage(ico)                        &
                                             + retained_carbon_fraction * delta_bleaf
                     end if

                     !---------------------------------------------------------------------!
                     !     Contribution due to the fact that c2n_leaf and c2n_storage may  !
                     ! be different.                                                       !
                     !---------------------------------------------------------------------!
                     csite%fsn_in(ipa)     = csite%fsn_in(ipa)                             &
                                           + delta_bleaf * cpatch%nplant(ico)              &
                                           * retained_carbon_fraction                      &
                                           * (1.0 / c2n_leaf(ipft) - 1.0/c2n_storage)
                     cpatch%bleaf(ico)     = bl_max

                     !---------------------------------------------------------------------!
                     !      Deduct the leaf drop from the carbon balance.                  !
                     !---------------------------------------------------------------------!
                     cpatch%cb          (13,ico) = cpatch%cb         (13,ico)              &
                                                 - cpatch%leaf_drop     (ico)
                     cpatch%cb_lightmax (13,ico) = cpatch%cb_lightmax(13,ico)              &
                                                 - cpatch%leaf_drop     (ico)
                     cpatch%cb_moistmax (13,ico) = cpatch%cb_moistmax(13,ico)              &
                                                 - cpatch%leaf_drop     (ico)
                     cpatch%cb_mlmax (13,ico)    = cpatch%cb_mlmax (13,ico) &
                                                 - cpatch%leaf_drop (ico)


                     !---------------------------------------------------------------------!
                  end if

                  !----- Set status flag. -------------------------------------------------!
                  if (bl_max == 0.0) then
                     cpatch%phenology_status(ico) = -2
                     cpatch%elongf(ico) = 0.
                  else
                     cpatch%elongf(ico) = 1.0 ! It should become green_leaf_factor...
                  end if
                  
               elseif (.not. drop_cold .and. cpatch%phenology_status(ico) == -2            &
                       .and. leaf_out_cold) then
                  !------------------------------------------------------------------------!
                  !      Update the phenology status (1 means that leaves are growing),    !
                  !------------------------------------------------------------------------!
                  if (veget_dyn_on) then
                     cpatch%phenology_status(ico) = 1
                     cpatch%elongf          (ico) = 1.0 
                  else
                     cpatch%phenology_status(ico) = 0
                     cpatch%elongf          (ico) = 1.0 
                     cpatch%bleaf(ico)  = cpoly%green_leaf_factor(ipft,isi)                &
                                        * size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
                     cpatch%balive(ico) = cpatch%balive(ico) + cpatch%bleaf(ico)
                  end if
                  !------------------------------------------------------------------------!
               end if


            case (3,4) 
               !---------------------------------------------------------------------------!
               !    Drought deciduous.  Here we must check two possibilities:              !
               !                                                                           !
               ! 1. The soil has been dry recently, and the plants have leaves, so we flag !
               !    them with phenology_status=0 (leaves not growing) and the plants will  !
               !    start losing their leaves;                                             !
               ! 2. The plant has no leaves, but the soil has started to come back to more !
               !    moist conditions.  Given this situation, leaves can start growing      !
               !    again.                                                                 !
               !---------------------------------------------------------------------------!
               !----- This is the first guess for the new elongation factor. --------------!
               elongf_try  = max(0.0, min (1.0, cpatch%paw_avg(ico)))
               !----- If extremely dry, force the cohort to shed all leaves... ------------!
               if (elongf_try < elongf_min) elongf_try = 0.0
               !---------------------------------------------------------------------------!



               !----- Find the maximum allowed leaf biomass. ------------------------------!
               bl_max = elongf_try * size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Delta_bleaf is the difference between the current leaf biomass and    !
               ! the maximum permitted given the soil moisture conditions.  If delta_bleaf !
               ! is positive, it means that the plant has more leaves than it should.      !
               !---------------------------------------------------------------------------!
               delta_bleaf = cpatch%bleaf(ico) - bl_max
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Check whether drought is becoming more or less severe.                !
               !---------------------------------------------------------------------------!
               if (delta_bleaf > 0.0) then
                  !------------------------------------------------------------------------!
                  !    Drought conditions are becoming more severe, drop leaves.           !
                  !------------------------------------------------------------------------!
                  if (elongf_try >= elongf_min) then
                     cpatch%phenology_status(ico) = -1
                  else
                     cpatch%phenology_status(ico) = -2
                  end if
                  cpatch%leaf_drop (ico) = (1.0 - retained_carbon_fraction) * delta_bleaf
                  cpatch%elongf    (ico) = elongf_try
                  !----- Adjust plant carbon pools. ---------------------------------------!
                  cpatch%bleaf     (ico) = bl_max
                  cpatch%balive    (ico) = cpatch%balive(ico)   - delta_bleaf
                  if (veget_dyn_on) then
                     cpatch%bstorage  (ico) = cpatch%bstorage(ico)                         &
                                            + retained_carbon_fraction * delta_bleaf
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Send the lost leaves to soil carbon and nitrogen pools.            !
                  !------------------------------------------------------------------------!
                  csite%fsc_in(ipa) = csite%fsc_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * f_labile_leaf(ipft)
                  csite%fsn_in(ipa) = csite%fsn_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * f_labile_leaf(ipft) / c2n_leaf(ipft)
                  csite%ssc_in(ipa) = csite%ssc_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * ( 1.0 - f_labile_leaf(ipft) )
                  csite%ssl_in(ipa) = csite%ssl_in(ipa)                                    &
                                    + cpatch%nplant(ico) * cpatch%leaf_drop(ico)           &
                                    * ( 1.0 - f_labile_leaf(ipft) )                        &
                                    * l2n_stem / c2n_stem(ipft)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Contribution due to the fact that c2n_leaf and c2n_storage may be  !
                  ! different.                                                             !
                  !------------------------------------------------------------------------!
                  csite%fsn_in(ipa) = csite%fsn_in(ipa) + delta_bleaf*cpatch%nplant(ico)   &
                                    * retained_carbon_fraction                             &
                                    * (1.0 / c2n_leaf(ipft) - 1.0/c2n_storage)
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !      Deduct the leaf drop from the carbon balance.                     !
                  !------------------------------------------------------------------------!
                  cpatch%cb          (13,ico) = cpatch%cb          (13,ico)                &
                                              - cpatch%leaf_drop      (ico)
                  cpatch%cb_lightmax (13,ico) = cpatch%cb_lightmax (13,ico)                &
                                              - cpatch%leaf_drop      (ico)
                  cpatch%cb_moistmax (13,ico) = cpatch%cb_moistmax (13,ico)                &
                                              - cpatch%leaf_drop      (ico)
                  cpatch%cb_mlmax (13,ico)    = cpatch%cb_mlmax    (13,ico)                &
                                              - cpatch%leaf_drop      (ico)

                  !------------------------------------------------------------------------!
               elseif (cpatch%phenology_status(ico) /= 0) then
                  !------------------------------------------------------------------------!
                  !       Elongation factor could increase, but we first check whether it  !
                  ! is safe to do so based on the phenology status.                        !
                  !------------------------------------------------------------------------!
                  select case(cpatch%phenology_status(ico))
                  case (1)
                     !----- Leaves were already growing, keep growing. --------------------!
                     cpatch%elongf          (ico) = elongf_try
                     !---------------------------------------------------------------------!
                  case (-1,-2)
                     !---------------------------------------------------------------------!
                     !     Leaves were dropping or gone, we first check that conditions    !
                     ! are really improving before we turn on leaf production.             !
                     !---------------------------------------------------------------------!
                     elongf_grow = min(1.0,max(elongf_flush,cpatch%elongf(ico)+0.02))
                     if (elongf_try >= elongf_grow) then
                        cpatch%elongf          (ico) = elongf_try
                        cpatch%phenology_status(ico) = 1
                     end if
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Conditions are slightly more humid.  Let them grow.                !
                  !------------------------------------------------------------------------!
                  if (.not. veget_dyn_on) then
                     cpatch%bleaf           (ico) = bl_max
                     cpatch%balive          (ico) = cpatch%balive(ico) - delta_bleaf
                  end if
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!




            !----- Update LAI, WAI, and CAI accordingly. ----------------------------------!
            call area_indices(cpatch, ico)
            !------------------------------------------------------------------------------!




            !----- Update above-ground biomass. -------------------------------------------!
            cpatch%agb(ico) = ed_biomass(cpatch, ico)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    The leaf biomass of the cohort has changed, update the vegetation energy  !
            ! using a constant temperature assumption.                                     !
            !------------------------------------------------------------------------------!
            old_leaf_hcap       = cpatch%leaf_hcap(ico)
            old_wood_hcap       = cpatch%wood_hcap(ico)
            call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdead(ico),cpatch%bsapwooda(ico)   &
                              ,cpatch%bbark(ico),cpatch%nplant(ico),cpatch%pft(ico)        &
                              ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
            call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap)
            call is_resolvable(csite,ipa,ico)
            !------------------------------------------------------------------------------!

            !----- Printing some debugging stuff if the code is set for it. ---------------!
            if (printphen) then
               ipft=cpatch%pft(ico)
               if (first_time(ipft)) then
                  first_time(ipft) = .false.
                  write (unit=40+ipft,fmt='(a10,6(1x,a12))')                               &
                     &'      TIME','       PATCH','      COHORT','      NPLANT'            &
                     &            ,'   LEAF_DROP','     PAW_AVG','      ELONGF'
               end if
            
               write (unit=40+ipft,fmt='(2(i2.2,a1),i4.4,2(1x,i12),5(1x,es12.5))')         &
                    current_time%month,'/',current_time%date,'/',current_time%year,ipa,ico &
                   ,cpatch%nplant(ico),cpatch%leaf_drop(ico),cpatch%paw_avg(ico)           &
                   ,cpatch%elongf(ico)
            end if
         end do cohortloop
      end do patchloop
      return
   end subroutine update_phenology
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine establishes whether it's time to drop leaves or start flushing    !
   ! them for cold deciduous or temperate drought deciduous.                               !
   ! MLO. Shouldn't we have a similar criterion for both tropical and temperate, based on  !
   ! a long term dry condition?                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine phenology_thresholds(daylight,soil_temp,sum_chd,sum_dgd,drop_cold            &
                                  ,leaf_out_cold)
      use phenology_coms, only : dl_tr        & ! intent(in)
                               , st_tr1       & ! intent(in)
                               , st_tr2       & ! intent(in)
                               , phen_a       & ! intent(in)
                               , phen_b       & ! intent(in)
                               , phen_c       & ! intent(in)
                               , iphen_scheme ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real                   , intent(in)    :: daylight      ! Daytime Length
      real                   , intent(in)    :: soil_temp     ! 
      real                   , intent(inout) :: sum_dgd       !
      real                   , intent(inout) :: sum_chd       !
      logical                , intent(out)   :: drop_cold     !
      logical                , intent(out)   :: leaf_out_cold !
      !----- Local variables --------------------------------------------------------------!
      real                                   :: gdd_threshold
      !------------------------------------------------------------------------------------!

      !----- Initialize variables. --------------------------------------------------------!
      drop_cold     = .false.
      leaf_out_cold = .false.
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check which phenology scheme to use.                                           !
      !------------------------------------------------------------------------------------!
      select case (iphen_scheme)
      case (1)
         !---------------------------------------------------------------------------------!
         !     Prescribed phenology, skip this and fill with the prescribed phenology.     !
         !---------------------------------------------------------------------------------!
         drop_cold     = .false.
         leaf_out_cold = .false.
         !---------------------------------------------------------------------------------!

      case default
         !---------------------------------------------------------------------------------!
         !    Predicted phenology, check if this is the time to drop leaves or to start    !
         ! flushing.                                                                       !
         !---------------------------------------------------------------------------------!

         !----- Too cold or too dark, time to shed leaves... ------------------------------!
         drop_cold = (daylight <= dl_tr .and. soil_temp < st_tr1) .or.  soil_temp < st_tr2
         !---------------------------------------------------------------------------------!


         !----- Warm again, time to flush leaves. -----------------------------------------!
         gdd_threshold = phen_a + phen_b * exp(phen_c * sum_chd)
         leaf_out_cold = sum_dgd >= gdd_threshold
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine phenology_thresholds
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine assigns the prescribed phenology.                                  !
   !---------------------------------------------------------------------------------------!
   subroutine assign_prescribed_phen(green_leaf_factor,leaf_aging_factor,dbh,height,pft    &
                                    ,drop_cold,leaf_out_cold,bl_max)
      use allometry     , only : size2bl
      use phenology_coms, only : elongf_min
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      logical, intent(out) :: drop_cold
      logical, intent(out) :: leaf_out_cold
      real   , intent(out) :: bl_max
      real   , intent(in)  :: green_leaf_factor
      real   , intent(in)  :: leaf_aging_factor
      real   , intent(in)  :: dbh
      real   , intent(in)  :: height
      integer, intent(in)  :: pft
      !------------------------------------------------------------------------------------!

      drop_cold     = green_leaf_factor /= leaf_aging_factor
      leaf_out_cold = green_leaf_factor > elongf_min .and. (.not. drop_cold)
      bl_max        = green_leaf_factor * size2bl(dbh, height, pft)

      return
   end subroutine assign_prescribed_phen
   !=======================================================================================!
   !=======================================================================================!
end module phenology_driv
!==========================================================================================!
!==========================================================================================!

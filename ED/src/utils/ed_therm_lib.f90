!==========================================================================================!
!==========================================================================================!
module ed_therm_lib
  
   contains
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   ! This function calculates the total heat capacity in (J/K/m2) of the cohort            !
   ! leaf biomass.  This function is primarily used to calculate leaf temperatures         !
   ! based on leaf energy.                                                                 !
   !                                                                                       !
   ! Inputs:                                                                               !
   ! + BLEAF        - the leaf biomass of the cohort in kgC/plant.  NOTE: it will be       !
   !                  converted to kg/m2, by the product of leaf_carbon in units kg/plant  !
   !                  and number of plants per square meter, and the conversion of carbon  !
   !                  to total biomass.  The right hand side of the main equation accounts !
   !                  for the mass of insterstitial water and its ability to hold energy.  !
   ! + BDEAD        - the structural wood biomass of the cohort in kgC/plant.              !
   ! + BSAPWOODA    - the above ground sapwood biomass of the cohort, in kgC/plant.        !
   ! + NPLANTS      - the number of plants per m2.                                         !
   ! + PFT          - the plant functional type of the current cohort, which may serve     !
   !                  for defining different parameterizations of specific heat capacity   !
   !                                                                                       !
   ! Ouputs:                                                                               !
   ! + LEAF_HCAP    - the leaf heat capacity, in J/m2/K.                                   !
   ! + WOOD_HCAP    - the wood heat capacity, in J/m2/K.                                   !
   !                                                                                       !
   ! These methods follow the ways of Gu et al. 2007, with the only difference that for    !
   ! non-green biomass we dropped the temperature dependence and assumed T=T3ple, just to  !
   ! make it simpler.  See the module pft_coms.f90 for a description of the parameters,    !
   ! and see ed_params.f90 for the setting of these parameters.                            !
   !                                                                                       !
   ! Reference:                                                                            !
   !                                                                                       !
   ! Gu, L., T. Meyers, S. G. Pallardy, 2007: Influences of biomass heat and biochemical   !
   !      energy storages on the land surface fluxes and radiative temperature.            !
   !      J. Geophys. Res., v. 112, doi: 10.1029/2006JD007425.                             !
   !---------------------------------------------------------------------------------------!
   subroutine calc_veg_hcap(bleaf,bdead,bsapwooda,nplant,pft                           &
                           ,broot,dbh,leaf_rwc,wood_rwc                                &
                           ,leaf_hcap,wood_hcap)
      use consts_coms          , only : cliq                ! ! intent(in)
      use pft_coms             , only : c_grn_leaf_dry      & ! intent(in)
                                      , wat_dry_ratio_grn   & ! intent(in)
                                      , c_ngrn_biom_dry     & ! intent(in)
                                      , wat_dry_ratio_ngrn  & ! intent(in)
                                      , leaf_water_sat      & ! intent(in)
                                      , wood_water_sat      & ! intent(in)
                                      , wood_rwc_min        & ! intent(in)
                                      , delta_c             & ! intent(in)
                                      , agf_bs              & ! intent(in)
                                      , C2B                 & ! intent(in)
                                      , brf_wd              ! ! intent(in)

      use rk4_coms             , only : ibranch_thermo      ! ! intent(in)
      use physiology_coms      , only : plant_hydro_scheme  ! ! intent(in)
      use allometry            , only : dbh2sf              ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real    , intent(in)    :: bleaf         ! Biomass of leaves              [kgC/plant]
      real    , intent(in)    :: bdead         ! Biomass of structural wood     [kgC/plant]
      real    , intent(in)    :: bsapwooda     ! Biomass of above ground sapwood[kgC/plant]
      real    , intent(in)    :: nplant        ! Number of plants               [ plant/m2]
      integer , intent(in)    :: pft           ! Plant functional type          [     ----]
      real    , intent(in)    :: broot         ! Biomass of fine roots          [kgC/plant]
      real    , intent(in)    :: dbh           ! Diameter at Breast Height      [       cm]
      real    , intent(in)    :: leaf_rwc      ! Leaf relative water content    [      0-1]
      real    , intent(in)    :: wood_rwc      ! Wood relative water content    [      0-1]
      real    , intent(out)   :: leaf_hcap     ! Leaf heat capacity             [   J/m2/K]
      real    , intent(out)   :: wood_hcap     ! Wood heat capacity             [   J/m2/K]
      !----- Local variables --------------------------------------------------------------!
      real                    :: bwood         ! Wood biomass                   [kgC/plant]
      real                    :: spheat_leaf   ! Leaf specific heat             [   J/kg/K]
      real                    :: spheat_wood   ! Wood specific heat             [   J/kg/K]
      !------------------------------------------------------------------------------------!

      select case (plant_hydro_scheme)
      case (0)
          ! original ED-2.2
          ! No plant hydraulics
          !------------------------------------------------------------------------------------!
          !    Here we decide whether we compute the branch heat capacity or not.              !
          !------------------------------------------------------------------------------------!
          select case (ibranch_thermo)
          case (0)
             !----- Skip it, the user doesn't want to solve for branches. ---------------------!
             spheat_wood = 0.
             bwood       = 0.
          case default
             !----- Find branch/twig specific heat and biomass. -------------------------------!
             spheat_wood = (c_ngrn_biom_dry(pft) + wat_dry_ratio_ngrn(pft) * cliq)             &
                         / (1. + wat_dry_ratio_ngrn(pft)) + delta_c(pft)
             bwood       = brf_wd(pft) * (bsapwooda + bdead*agf_bs(pft))
          end select

          !----- Find the leaf specific heat. -------------------------------------------------!
          spheat_leaf = (c_grn_leaf_dry(pft) + wat_dry_ratio_grn(pft) * cliq)                  &
                      / (1. + wat_dry_ratio_grn(pft))

          !------------------------------------------------------------------------------------!
          !     The heat capacity is specific heat times the plant density times the leaf/wood !
          ! biomass.                                                                           !
          !------------------------------------------------------------------------------------!
          leaf_hcap = nplant * C2B * bleaf * spheat_leaf * (1. + wat_dry_ratio_grn (pft))
          wood_hcap = nplant * C2B * bwood * spheat_wood * (1. + wat_dry_ratio_ngrn(pft))
          !------------------------------------------------------------------------------------!
      case (1,2)
          ! with plant hydraulics
          !------------------------------------------------------------------------------------!
          !    Here we decide whether we compute the branch heat capacity or not.              !
          !------------------------------------------------------------------------------------!
          select case (ibranch_thermo)
          case (0)
             !----- Skip it, the user doesn't want to solve for branches. ---------------------!
             spheat_wood = 0.
             bwood       = 0.
          case default
             !----- Find branch/twig specific heat and biomass. -------------------------------!
             !Here we consider the relative water content of wood
             spheat_wood = (c_ngrn_biom_dry(pft) + wood_water_sat(pft) * wood_rwc * cliq)      &
                         / (1. + wood_water_sat(pft) * wood_rwc) + delta_c(pft) * wood_rwc
             bwood       = bdead * dbh2sf(dbh,pft) + broot
          end select

          !----- Find the leaf specific heat. -------------------------------------------------!
          spheat_leaf = (c_grn_leaf_dry(pft) + leaf_water_sat(pft) * leaf_rwc * cliq)       &
                      / (1. + leaf_water_sat(pft) * leaf_rwc)

          !------------------------------------------------------------------------------------!
          !     The heat capacity is specific heat times the plant density times the leaf/wood !
          ! biomass.                                                                           !
          !------------------------------------------------------------------------------------!
          leaf_hcap = nplant * C2B * bleaf * spheat_leaf * (1. + leaf_water_sat(pft) * leaf_rwc)
          wood_hcap = nplant * C2B * bwood * spheat_wood * (1. + wood_water_sat(pft) * wood_rwc)
          !------------------------------------------------------------------------------------!
      case (-1,-2)
          ! with plant hydraulics but fix the hcap to the value for saturated
          ! tissues
          !------------------------------------------------------------------------------------!
          !    Here we decide whether we compute the branch heat capacity or not.              !
          !------------------------------------------------------------------------------------!
          select case (ibranch_thermo)
          case (0)
             !----- Skip it, the user doesn't want to solve for branches. ---------------------!
             spheat_wood = 0.
             bwood       = 0.
          case default
             !----- Find branch/twig specific heat and biomass. -------------------------------!
             !Here we consider the relative water content of wood
             spheat_wood = (c_ngrn_biom_dry(pft) + wood_water_sat(pft) * cliq)      &
                         / (1. + wood_water_sat(pft)) + delta_c(pft)
             bwood       = bdead * dbh2sf(dbh,pft) + broot
          end select

          !----- Find the leaf specific heat. -------------------------------------------------!
          spheat_leaf = (c_grn_leaf_dry(pft) + leaf_water_sat(pft) * cliq)       &
                      / (1. + leaf_water_sat(pft))

          !------------------------------------------------------------------------------------!
          !     The heat capacity is specific heat times the plant density times the leaf/wood !
          ! biomass.                                                                           !
          !------------------------------------------------------------------------------------!
          leaf_hcap = nplant * C2B * bleaf * spheat_leaf * (1. + leaf_water_sat(pft))
          wood_hcap = nplant * C2B * bwood * spheat_wood * (1. + wood_water_sat(pft))
          !------------------------------------------------------------------------------------!

      end select

      return
   end subroutine calc_veg_hcap
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine updates the vegetation energy when the plant heat capacity has     !
   ! changed. This routine should be used only when leaf or structural biomass has         !
   ! changed, it should never be used in fast time steps.                                  !
   !                                                                                       !
   !     We look at leaf and wood separately, but the idea is the same.  When heat         !
   ! capacity is zero (i.e., no leaves or not solving branchwood thermodynamics), we       !
   ! cannot find the temperature using uextcm2tl because it is a singularity.  Notice that !
   ! this is different than skipping when cohorts are not resolvable...  If the cohort is  !
   ! not resolvable but still has some heat capacity, we should update internal energy     !
   ! using the traditional method, and NEVER force the heat capacity to be zero, otherwise !
   ! we violate the fact that heat capacity is a linear function of mass and this will     !
   ! cause problems during the fusion/splitting process.                                   !
   !                                                                                       !
   !    The "cweh" mean "consistent water&energy&hcap" assumption                          !
   !---------------------------------------------------------------------------------------!
   subroutine update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap)
      use ed_state_vars, only : sitetype   & ! structure
                              , patchtype  ! ! structure
      use therm_lib    , only : uextcm2tl  & ! subroutine
                              , cmtl2uext  ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype) , target     :: csite
      integer        , intent(in) :: ipa
      integer        , intent(in) :: ico
      real           , intent(in) :: old_leaf_hcap
      real           , intent(in) :: old_wood_hcap
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer    :: cpatch
      real                        :: new_temp
      real                        :: new_fliq
      integer                     :: kclosest
      integer                     :: k
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)


      !------------------------------------------------------------------------------------!
      !     Leaves.  Check whether heat capacity is zero or not.                           !
      !------------------------------------------------------------------------------------!
      if (cpatch%leaf_hcap(ico) == 0. ) then
         cpatch%leaf_energy(ico) = 0.
         cpatch%leaf_water(ico)  = 0.
         cpatch%leaf_fliq(ico)   = 0.
         if (cpatch%hite(ico) > csite%total_sfcw_depth(ipa)) then
            !----- Plant is exposed, set temperature to the canopy temperature. -----------!
            cpatch%leaf_temp(ico) = csite%can_temp(ipa)
         else
            !----- Find the snow layer that is the closest to where the leaves would be. --!
            do k = csite%nlev_sfcwater(ipa), 1, -1
               if (sum(csite%sfcwater_depth(1:k,ipa)) >= cpatch%hite(ico)) then
                  kclosest = k
               end if
            end do
            cpatch%leaf_temp(ico) = csite%sfcwater_tempk(kclosest,ipa)
         end if
         !---------------------------------------------------------------------------------!

      else
         !---------------------------------------------------------------------------------!
         !     Heat capacity is not zero.  Since we track leaf temperature and liquid      !
         ! fraction of water held by leaves, we can recalculate the internal energy by     !
         ! just switching the old heat capacity by the new one.                            !
         !---------------------------------------------------------------------------------!
         cpatch%leaf_energy(ico) = cmtl2uext(cpatch%leaf_hcap(ico),cpatch%leaf_water(ico)  &
                                            ,cpatch%leaf_temp(ico),cpatch%leaf_fliq(ico) )
         !---------------------------------------------------------------------------------!



         !----- This is a sanity check, it can be removed if it doesn't crash. ------------!
         call uextcm2tl(cpatch%leaf_energy(ico),cpatch%leaf_water(ico)                     &
                       ,cpatch%leaf_hcap(ico),new_temp,new_fliq)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     In case the temperature is different, give the user the bad news...         !
         !---------------------------------------------------------------------------------!
         if (abs(new_temp - cpatch%leaf_temp(ico)) > 0.1) then
            write(unit=*,fmt='(a)') '-----------------------------------------------------'
            write(unit=*,fmt='(a)') ' LEAF ENERGY CONSERVATION FAILED!:'
            write(unit=*,fmt='(a)') '-----------------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.5)') ' Old temperature:  ',cpatch%leaf_temp(ico)
            write(unit=*,fmt='(a,1x,es12.5)') ' New temperature:  ',new_temp
            write(unit=*,fmt='(a,1x,es12.5)') ' Old heat capacity:',old_leaf_hcap
            write(unit=*,fmt='(a,1x,es12.5)') ' New heat capacity:',cpatch%leaf_hcap(ico)
            write(unit=*,fmt='(a,1x,es12.5)') ' Leaf energy:      ',cpatch%leaf_energy(ico)
            write(unit=*,fmt='(a,1x,es12.5)') ' Leaf water:       ',cpatch%leaf_water(ico)
            write(unit=*,fmt='(a)') '-----------------------------------------------------'
            call fatal_error('Leaf energy is leaking!!!','update_veg_energy_cweh'          &
                            &,'ed_therm_lib.f90')
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Wood.  Check whether heat capacity is zero or not.                             !
      !------------------------------------------------------------------------------------!
      if (cpatch%wood_hcap(ico) == 0. ) then
         cpatch%wood_energy(ico) = 0.
         cpatch%wood_water(ico)  = 0.
         cpatch%wood_fliq(ico)   = 0.
         if (cpatch%hite(ico) > csite%total_sfcw_depth(ipa)) then
            !----- Plant is exposed, set temperature to the canopy temperature. -----------!
            cpatch%wood_temp(ico) = csite%can_temp(ipa)
         else
            !----- Find the snow layer that is the closest to where the leaves would be. --!
            do k = csite%nlev_sfcwater(ipa), 1, -1
               if (sum(csite%sfcwater_depth(1:k,ipa)) >= cpatch%hite(ico)) then
                  kclosest = k
               end if
            end do
            cpatch%wood_temp(ico) = csite%sfcwater_tempk(kclosest,ipa)
         end if
         !---------------------------------------------------------------------------------!

      else
         !---------------------------------------------------------------------------------!
         !     Heat capacity is not zero.  Since we track leaf temperature and liquid      !
         ! fraction of water held by leaves, we can recalculate the internal energy by     !
         ! just switching the old heat capacity by the new one.                            !
         !---------------------------------------------------------------------------------!
         cpatch%wood_energy(ico) = cmtl2uext(cpatch%wood_hcap(ico),cpatch%wood_water(ico)  &
                                            ,cpatch%wood_temp(ico),cpatch%wood_fliq (ico) )
         !---------------------------------------------------------------------------------!



         !----- This is a sanity check, it can be removed if it doesn't crash. ------------!
         call uextcm2tl(cpatch%wood_energy(ico),cpatch%wood_water(ico)                     &
                       ,cpatch%wood_hcap(ico),new_temp,new_fliq)
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     In case the temperature is different, give the user the bad news...         !
         !---------------------------------------------------------------------------------!
         if (abs(new_temp - cpatch%wood_temp(ico)) > 0.1) then
            write(unit=*,fmt='(a)') '-----------------------------------------------------'
            write(unit=*,fmt='(a)') ' WOOD ENERGY CONSERVATION FAILED!:'
            write(unit=*,fmt='(a)') '-----------------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.5)') ' Old temperature:  ',cpatch%wood_temp(ico)
            write(unit=*,fmt='(a,1x,es12.5)') ' New temperature:  ',new_temp
            write(unit=*,fmt='(a,1x,es12.5)') ' Old heat capacity:',old_wood_hcap
            write(unit=*,fmt='(a,1x,es12.5)') ' New heat capacity:',cpatch%wood_hcap(ico)
            write(unit=*,fmt='(a,1x,es12.5)') ' Wood energy:      ',cpatch%wood_energy(ico)
            write(unit=*,fmt='(a,1x,es12.5)') ' Wood water:       ',cpatch%wood_water(ico)
            write(unit=*,fmt='(a)') '-----------------------------------------------------'
            call fatal_error('Wood energy is leaking!!!','update_veg_energy_cweh'          &
                            &,'ed_therm_lib.f90')
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      return
   end subroutine update_veg_energy_cweh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This routine computes ground_shv, which is is the saturation mixing ratio of     !
   ! the top soil or snow surface and is used for dew formation and snow evaporation.      !
   ! References:                                                                           !
   !                                                                                       !
   ! P86  - Passerat de Silans, A., 1986: Transferts de masse et de chaleur dans un sol    !
   !        stratifi� soumis � une excitation amtosph�rique naturelle. Comparaison:        !
   !        Mod�les-exp�rience. Thesis, Institut National Polytechnique de Grenoble.       !
   !                                                                                       !
   ! NP89 - Noilhan, J., S. Planton, 1989: A simple parameterization of land surface       !
   !        processes for meteorological models. Mon. Wea. Rev., 117, 536-549.             !
   !                                                                                       !
   ! MN91 - Mahfouf, J. F., J. Noilhan, 1991: Comparative study of various formulations of !
   !        evaporation from bare soil using in situ data. J. Appl. Meteorol., 30,         !
   !        1354-1365.                                                                     !
   !                                                                                       !
   ! LP92 - Lee, T. J., R. A. Pielke, 1992: Estimating the soil surface specific humidity  !
   !        J. Appl. Meteorol., 31, 480-484.                                               !
   !                                                                                       !
   ! LP93 - Lee, T. J., R. A. Pielke, 1993: CORRIGENDUM, J. Appl. Meteorol., 32, 580.      !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine ed_grndvap(ksn,nsoil,topsoil_water,topsoil_temp,topsoil_fliq,sfcwater_temp   &
                        ,sfcwater_frac,can_prss,can_shv,ground_shv                         &
                        ,ground_ssh,ground_temp,ground_fliq,ggsoil)

      use canopy_air_coms, only : ied_grndvap       & ! intent(in)
                                , ggsoil0           & ! intent(in)
                                , kksoil            ! ! intent(in)
      use soil_coms      , only : soil              & ! intent(in)
                                , matric_potential  ! ! function
      use consts_coms    , only : pi1               & ! intent(in)
                                , wdns              & ! intent(in)
                                , gorh2o            & ! intent(in)
                                , lnexp_min         & ! intent(in)
                                , tiny_num          & ! intent(in)
                                , huge_num          ! ! intent(in)
      use therm_lib      , only : qslif             ! ! function
      use ed_max_dims    , only : n_pft             ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer     , intent(in)  :: ksn           ! # of surface water layers    [     ----]
      integer     , intent(in)  :: nsoil         ! Soil type                    [     ----]
      real(kind=4), intent(in)  :: topsoil_water ! Top soil water               [m�_h2o/m�]
      real(kind=4), intent(in)  :: topsoil_temp  ! Top soil temperature         [        K]
      real(kind=4), intent(in)  :: topsoil_fliq  ! Top soil liquid water frac.  [       --]
      real(kind=4), intent(in)  :: sfcwater_temp ! Snow/water temperature       [        K]
      real(kind=4), intent(in)  :: sfcwater_frac ! Snow/water liq. water frac.  [       --]
      real(kind=4), intent(in)  :: can_prss      ! canopy pressure              [       Pa]
      real(kind=4), intent(in)  :: can_shv       ! canopy vapour spec humidity  [kg_vap/kg]
      real(kind=4), intent(out) :: ground_shv    ! ground equilibrium spec hum  [kg_vap/kg]
      real(kind=4), intent(out) :: ground_ssh    ! sfc. saturation spec. hum.   [kg_vap/kg]
      real(kind=4), intent(out) :: ground_temp   ! Surface temperature          [        K]
      real(kind=4), intent(out) :: ground_fliq   ! Surface liquid water frac.   [       --]
      real(kind=4), intent(out) :: ggsoil        ! Soil conductance for evap.   [      m/s]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)              :: slpotvn       ! soil water potential         [        m]
      real(kind=4)              :: alpha         ! alpha term (Lee-Pielke,1992) [     ----]
      real(kind=4)              :: beta          ! beta term  (Lee-Pielke,1992) [     ----]
      real(kind=4)              :: lnalpha       ! ln(alpha)                    [     ----]
      real(kind=4)              :: smterm        ! soil moisture term           [     ----]
      real(kind=4)              :: topsoil_shv   ! ground equilibrium spec hum  [kg_vap/kg]
      real(kind=4)              :: topsoil_ssh   ! sfc. saturation spec. hum.   [kg_vap/kg]
      real(kind=4)              :: sfcwater_shv  ! ground equilibrium spec hum  [kg_vap/kg]
      real(kind=4)              :: sfcwater_ssh  ! sfc. saturation spec. hum.   [kg_vap/kg]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Topsoil_shv is the effective specific humidity of soil.  This value is a      !
      ! combination of the canopy air specific humidity, the saturation specific humidity  !
      ! at the soil temperature.  When the soil tends to dry air soil moisture,            !
      ! topsoil_shv tends to the canopy air space specific humidity, whereas topsoil_shv   !
      ! tends to the saturation value when the soil moisture is near or above field        !
      ! capacity.  These tendencies are determined by the alpha and beta parameters.       !
      !------------------------------------------------------------------------------------!
      !----- Compute the saturation specific humidity at top soil temperature. ------------!
      topsoil_ssh  = qslif(can_prss,topsoil_temp)
      !----- Determine alpha. -------------------------------------------------------------!
      slpotvn      = matric_potential(nsoil,topsoil_water)
      lnalpha      = gorh2o * slpotvn / topsoil_temp
      if (lnalpha > lnexp_min) then
         alpha   = exp(lnalpha)
      else
         alpha   = 0.0
      end if
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Determine Beta, following NP89.  However, because we want evaporation to be    !
      ! shut down when the soil approaches the dry air soil moisture, we offset both       !
      ! the soil moisture and field capacity to the soil moisture above dry air soil.      !
      ! This is necessary to avoid evaporation to be large just slightly above the dry air !
      ! soil, which would otherwise happen, especially for those soil types rich in clay.  !
      !------------------------------------------------------------------------------------!
      smterm = min(1.0, max(0.0, (topsoil_water       - soil(nsoil)%soilcp)                &
                               / (soil(nsoil)%sfldcap - soil(nsoil)%soilcp) ))
      beta   = 0.5 * (1.0 - cos (smterm * pi1))
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Decide which method to use to find the ground water vapour mixing ratio.       !
      !------------------------------------------------------------------------------------!
      select case (ied_grndvap)
      case (0)
         !----- LP92 method. --------------------------------------------------------------!
         topsoil_shv = max(can_shv, topsoil_ssh * alpha * beta + (1.0 - beta) * can_shv)
         ggsoil      = huge_num
         !---------------------------------------------------------------------------------!

      case (1)
         !----- MH91, test 1. -------------------------------------------------------------!
         topsoil_shv = max(can_shv, topsoil_ssh * beta)
         ggsoil      = huge_num
         !---------------------------------------------------------------------------------!

      case (2)
         !----- MH91, test 2. -------------------------------------------------------------!
         topsoil_shv = topsoil_ssh
         ggsoil      = ggsoil0 * exp(kksoil * smterm)
         !---------------------------------------------------------------------------------!

      case (3)
         !----- MH91, test 3. -------------------------------------------------------------!
         topsoil_shv = max(can_shv, topsoil_ssh * beta + (1.0 - beta) * can_shv)
         ggsoil      = huge_num
         !---------------------------------------------------------------------------------!

      case (4)
         !---------------------------------------------------------------------------------!
         !     MH91, test 4.                                                               !
         !---------------------------------------------------------------------------------!
         topsoil_shv = max(can_shv, topsoil_ssh * alpha)
         ggsoil      = ggsoil0 * exp(kksoil * smterm)
         !---------------------------------------------------------------------------------!

      case (5)
         !---------------------------------------------------------------------------------!
         !     Combination of NP89 and P86.                                                !
         !---------------------------------------------------------------------------------!
         topsoil_shv = max(can_shv, topsoil_ssh * beta)
         ggsoil      = ggsoil0 * exp(kksoil * smterm)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !----- Compute the saturation specific humidity at ground temperature. --------------!
      sfcwater_ssh = qslif(can_prss,sfcwater_temp)
      !----- The ground specific humidity in this case is just the saturation value. ------!
      sfcwater_shv = sfcwater_ssh
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      The properties are weighted averages with area being the weight.              !
      !------------------------------------------------------------------------------------!
      ground_temp = topsoil_temp
      ground_fliq = topsoil_fliq
      ground_ssh  = topsoil_ssh
      ground_shv  = topsoil_shv
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Ggsoil is corrected as the weighted averages of the resistances (assumed to  !
      ! be zero for temporary surface water).                                              !
      !------------------------------------------------------------------------------------!
      if (ggsoil /= huge_num .and. ksn > 0) then
         ggsoil      = min(huge_num, ggsoil / max(tiny_num,(1.0 - sfcwater_frac)))
      end if
      !------------------------------------------------------------------------------------!


      return
   end subroutine ed_grndvap
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This routine computes ground_shv, which is is the saturation mixing ratio of     !
   ! the top soil or snow surface and is used for dew formation and snow evaporation.      !
   ! References:                                                                           !
   !                                                                                       !
   ! P86  - Passerat de Silans, A., 1986: Transferts de masse et de chaleur dans un sol    !
   !        stratifi� soumis � une excitation amtosph�rique naturelle. Comparaison:        !
   !        Mod�les-exp�rience. Thesis, Institut National Polytechnique de Grenoble.       !
   !                                                                                       !
   ! NP89 - Noilhan, J., S. Planton, 1989: A simple parameterization of land surface       !
   !        processes for meteorological models. Mon. Wea. Rev., 117, 536-549.             !
   !                                                                                       !
   ! MN91 - Mahfouf, J. F., J. Noilhan, 1991: Comparative study of various formulations of !
   !        evaporation from bare soil using in situ data. J. Appl. Meteorol., 30,         !
   !        1354-1365.                                                                     !
   !                                                                                       !
   ! LP92 - Lee, T. J., R. A. Pielke, 1992: Estimating the soil surface specific humidity  !
   !        J. Appl. Meteorol., 31, 480-484.                                               !
   !                                                                                       !
   ! LP93 - Lee, T. J., R. A. Pielke, 1993: CORRIGENDUM, J. Appl. Meteorol., 32, 580.      !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine ed_grndvap8(ksn,topsoil_water,topsoil_temp,topsoil_fliq,sfcwater_temp        &
                         ,sfcwater_frac,can_prss,can_shv,ground_shv                        &
                         ,ground_ssh,ground_temp,ground_fliq,ggsoil)
      use canopy_air_coms, only : ied_grndvap       & ! intent(in)
                                , ggsoil08          & ! intent(in)
                                , kksoil8           ! ! intent(in)
      use soil_coms      , only : soil8             & ! intent(in)
                                , matric_potential8 ! ! function
      use consts_coms    , only : pi18              & ! intent(in)
                                , wdns8             & ! intent(in)
                                , gorh2o8           & ! intent(in)
                                , lnexp_min8        & ! intent(in)
                                , tiny_num8         & ! intent(in)
                                , huge_num8         ! ! intent(in)
      use therm_lib8     , only : qslif8            ! ! function
      use rk4_coms       , only : rk4site           ! ! intent(in)
      use grid_coms      , only : nzg               ! ! intent(in)
      use ed_max_dims    , only : n_pft             ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer     , intent(in)  :: ksn           ! # of surface water layers    [     ----]
      real(kind=8), intent(in)  :: topsoil_water ! Top soil water               [m�_h2o/m�]
      real(kind=8), intent(in)  :: topsoil_temp  ! Top soil temperature         [        K]
      real(kind=8), intent(in)  :: topsoil_fliq  ! Top soil liquid water frac.  [       --]
      real(kind=8), intent(in)  :: sfcwater_temp ! Snow/water temperature       [        K]
      real(kind=8), intent(in)  :: sfcwater_frac ! Snow/water liq. water frac.  [       --]
      real(kind=8), intent(in)  :: can_prss      ! canopy pressure              [       Pa]
      real(kind=8), intent(in)  :: can_shv       ! canopy vapour spec humidity  [kg_vap/kg]
      real(kind=8), intent(out) :: ground_shv    ! ground equilibrium spec hum  [kg_vap/kg]
      real(kind=8), intent(out) :: ground_ssh    ! sfc. saturation spec. hum.   [kg_vap/kg]
      real(kind=8), intent(out) :: ground_temp   ! Surface temperature          [        K]
      real(kind=8), intent(out) :: ground_fliq   ! Surface liquid water frac.   [       --]
      real(kind=8), intent(out) :: ggsoil        ! Soil conductance for evap.   [      m/s]
      !----- Local variables --------------------------------------------------------------!
      integer                   :: nsoil         ! Soil type                    [     ----]
      real(kind=8)              :: slpotvn       ! soil water potential         [        m]
      real(kind=8)              :: alpha         ! alpha term (Lee-Pielke,1992) [     ----]
      real(kind=8)              :: beta          ! beta term  (Lee-Pielke,1992) [     ----]
      real(kind=8)              :: lnalpha       ! ln(alpha)                    [     ----]
      real(kind=8)              :: smterm        ! soil moisture term           [     ----]
      real(kind=8)              :: topsoil_shv   ! ground equilibrium spec hum  [kg_vap/kg]
      real(kind=8)              :: topsoil_ssh   ! sfc. saturation spec. hum.   [kg_vap/kg]
      real(kind=8)              :: sfcwater_shv  ! ground equilibrium spec hum  [kg_vap/kg]
      real(kind=8)              :: sfcwater_ssh  ! sfc. saturation spec. hum.   [kg_vap/kg]
      !------------------------------------------------------------------------------------!


      !------ Soil type at the top layer. -------------------------------------------------!
      nsoil = rk4site%ntext_soil(nzg)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Topsoil_shv is the effective specific humidity of soil.  This value is a      !
      ! combination of the canopy air specific humidity, the saturation specific humidity  !
      ! at the soil temperature.  When the soil tends to dry air soil moisture,            !
      ! topsoil_shv tends to the canopy air space specific humidity, whereas topsoil_shv   !
      ! tends to the saturation value when the soil moisture is near or above field        !
      ! capacity.  These tendencies are determined by the alpha and beta parameters.       !
      !------------------------------------------------------------------------------------!
      !----- Compute the saturation specific humidity at top soil temperature. ------------!
      topsoil_ssh  = qslif8(can_prss,topsoil_temp)
      !----- Determine alpha. -------------------------------------------------------------!
      slpotvn      = matric_potential8(nsoil,topsoil_water)
      lnalpha      = gorh2o8 * slpotvn / topsoil_temp
      if (lnalpha > lnexp_min8) then
         alpha   = exp(lnalpha)
      else
         alpha   = 0.d0
      end if
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Determine Beta, following NP89.  However, because we want evaporation to be    !
      ! shut down when the soil approaches the dry air soil moisture, we offset both       !
      ! the soil moisture and field capacity to the soil moisture above dry air soil.      !
      ! This is necessary to avoid evaporation to be large just slightly above the dry air !
      ! soil, which would otherwise happen, especially for those soil types rich in clay.  !
      !------------------------------------------------------------------------------------!
      smterm = min(1.d0, max(0.d0, (topsoil_water        - soil8(nsoil)%soilcp)            &
                                 / (soil8(nsoil)%sfldcap - soil8(nsoil)%soilcp) ))
      beta   = 5.d-1 * (1.d0 - cos (smterm * pi18))
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Decide which method to use to find the ground water vapour mixing ratio.       !
      !------------------------------------------------------------------------------------!
      select case (ied_grndvap)
      case (0)
         !----- LP92 method. --------------------------------------------------------------!
         topsoil_shv = max(can_shv, topsoil_ssh * alpha * beta + (1.d0 - beta) * can_shv)
         ggsoil      = huge_num8
         !---------------------------------------------------------------------------------!

      case (1)
         !----- MH91, test 1. -------------------------------------------------------------!
         topsoil_shv = max(can_shv, topsoil_ssh * beta)
         ggsoil      = huge_num8
         !---------------------------------------------------------------------------------!

      case (2)
         !----- MH91, test 2. -------------------------------------------------------------!
         topsoil_shv = topsoil_ssh
         ggsoil      = ggsoil08 * exp(kksoil8 * smterm)
         !---------------------------------------------------------------------------------!

      case (3)
         !----- MH91, test 3. -------------------------------------------------------------!
         topsoil_shv = max(can_shv, topsoil_ssh * beta + (1.d0 - beta) * can_shv)
         ggsoil      = huge_num8
         !---------------------------------------------------------------------------------!

      case (4)
         !---------------------------------------------------------------------------------!
         !     MH91, test 4.                                                               !
         !---------------------------------------------------------------------------------!
         topsoil_shv = max(can_shv, topsoil_ssh * alpha)
         ggsoil      = ggsoil08 * exp(kksoil8 * smterm)
         !---------------------------------------------------------------------------------!

      case (5)
         !---------------------------------------------------------------------------------!
         !     Combination of NP89 and P86.                                                !
         !---------------------------------------------------------------------------------!
         topsoil_shv = max(can_shv, topsoil_ssh * beta)
         ggsoil      = ggsoil08 * exp(kksoil8 * smterm)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !----- Compute the saturation specific humidity at ground temperature. --------------!
      sfcwater_ssh = qslif8(can_prss,sfcwater_temp)
      !----- The ground specific humidity in this case is just the saturation value. ------!
      sfcwater_shv = sfcwater_ssh
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      The properties are weighted averages with area being the weight.              !
      !------------------------------------------------------------------------------------!
      ground_temp = topsoil_temp
      ground_fliq = topsoil_fliq
      ground_ssh  = topsoil_ssh
      ground_shv  = topsoil_shv
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Ggsoil is corrected as the weighted averages of the resistances (assumed to  !
      ! be zero for temporary surface water).                                              !
      !------------------------------------------------------------------------------------!
      if (ggsoil /= huge_num8 .and. ksn > 0) then
         ggsoil      = min(huge_num8, ggsoil / max(tiny_num8,(1.d0 - sfcwater_frac)))
      end if
      !------------------------------------------------------------------------------------!


      return
   end subroutine ed_grndvap8
   !=======================================================================================!
   !=======================================================================================!
 end module ed_therm_lib
!==========================================================================================!
!==========================================================================================!

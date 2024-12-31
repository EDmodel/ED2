!==========================================================================================!
!==========================================================================================!
!    MODULE RK4_DERIVS.F90  
!
!> \brief This module contains the derivatives of most short-term fluxes.
!> \details Three subroutines live in this module.  leaf_derivs is a wrapper; 
!!          leaftw_derivs solves mostly the ground derivatives (soils and temporary surface
!!          water), and canopy_derivs_two solves derivatives at the cohort scale and the 
!!          canopy air space at patch  level.  These routines are largely based on BRAMS
!!          LEAF-2/LEAF-3 model.
!> \author  Adapted from LEAF-2 by David Medvigy.  Further updates to make it compatible
!!          with LEAF-3 and to conserve energy and water by Ryan Knox and Marcos Longo.
!> \author  Converted to module to eliminate interfaces by Marcos Longo
!------------------------------------------------------------------------------------------!
module rk4_derivs
   contains
   !=======================================================================================!
   !=======================================================================================!
   ! Subroutine leaf_derivs                                                                !
   !                                                                                       !
   !     This subroutine finds the fast-scale derivatives at canopy, soil, and leaf sfc.   !
   ! This subroutine is based on LEAF-3, except that here only the derivative is computed, !
   ! whereas in LEAF-3 the actual step is done at once. This derivative will be used for   !
   ! the  Runge-Kutta integration step.                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine leaf_derivs(initp,dinitp,csite,ipa,ibuff,dt,is_hybrid)

      use rk4_coms               , only : rk4patchtype       & ! structure
                                        , rk4aux             & ! intent(out)
                                        , zero_rk4_aux       & ! sub-routine
                                        , zero_rk4_patch     & ! sub-routine
                                        , zero_rk4_cohort    ! ! sub-routine
      use ed_state_vars          , only : sitetype           & ! structure
                                        , polygontype        ! ! structure
      use grid_coms              , only : nzg                & ! intent(in)
                                        , nzs                ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target     :: initp     ! Structure with RK4 intermediate state
      type(rk4patchtype) , target     :: dinitp    ! Structure with RK4 derivatives
      type(sitetype)     , target     :: csite     ! This site (with previous values);
      integer            , intent(in) :: ipa       ! Patch ID
      integer            , intent(in) :: ibuff     ! The shared memory processor index
                                                   ! for the buffer space
      real(kind=8)       , intent(in) :: dt        ! Current time step
      logical            , intent(in) :: is_hybrid ! Flag to tell whether it is a hybrid
                                                   !    solver solution.
      !------------------------------------------------------------------------------------!


      !---- Flush all derivatives to zero. ------------------------------------------------!
      call zero_rk4_patch (dinitp)
      call zero_rk4_cohort(dinitp)
      !------------------------------------------------------------------------------------!


      !---- Flush auxiliary variables to zero. --------------------------------------------!
      call zero_rk4_aux(rk4aux(ibuff))
      !------------------------------------------------------------------------------------!


      !----- Find the derivatives. --------------------------------------------------------!
      call leaftw_derivs(nzg,nzs,initp,dinitp,csite,ipa,ibuff,dt,is_hybrid)
      !------------------------------------------------------------------------------------!

      return
   end subroutine leaf_derivs
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine leaftw_derivs(mzg,mzs,initp,dinitp,csite,ipa,ibuff,dt,is_hybrid)
      use ed_max_dims          , only : nzgmax                & ! intent(in)
                                      , nzsmax                ! ! intent(in)
      use consts_coms          , only : cliq8                 & ! intent(in)
                                      , cph2o8                & ! intent(in)
                                      , wdns8                 & ! intent(in)
                                      , wdnsi8                & ! intent(in)
                                      , lnexp_min8            & ! intent(in)
                                      , tiny_num8             ! ! intent(in)
      use soil_coms            , only : soil8                 & ! intent(in)
                                      , dslz8                 & ! intent(in)
                                      , dslzi8                & ! intent(in)
                                      , infiltration_method   & ! intent(in)
                                      , dslzti8               & ! intent(in)
                                      , slcons18              & ! intent(in)
                                      , slzt8                 & ! intent(in)
                                      , dslzt8                & ! intent(in)
                                      , ss                    & ! intent(in)
                                      , sin_sldrain8          & ! intent(in)
                                      , hydr_conduct8         ! ! function
      use rk4_coms             , only : checkbudget           & ! intent(in)
                                      , print_detailed        & ! intent(in)
                                      , rk4site               & ! intent(in)
                                      , rk4patchtype          & ! structure
                                      , rk4aux                ! ! intent(out)
      use ed_state_vars        , only : sitetype              & ! structure
                                      , patchtype             & ! structure
                                      , polygontype           ! ! structure
      use therm_lib8           , only : tl2uint8              ! ! functions
      use ed_misc_coms         , only : fast_diagnostics      ! ! intent(in)
      use physiology_coms      , only : h2o_plant_lim         & ! intent(in)
                                      , plant_hydro_scheme    ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype)  , target     :: initp     ! RK4 structure, intermediate step
      type(rk4patchtype)  , target     :: dinitp    ! RK4 structure, derivatives
      type(sitetype)      , target     :: csite     ! Current site (before integration)
      integer             , intent(in) :: ipa       ! Current patch number
      integer             , intent(in) :: ibuff     ! The shared memory processor index
                                                    ! for the buffer space
      integer             , intent(in) :: mzg       ! Number of ground layers
      integer             , intent(in) :: mzs       ! Max. number of TSW layers
      real(kind=8)        , intent(in) :: dt        ! Timestep
      logical             , intent(in) :: is_hybrid ! Hybrid solver?
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer    :: cpatch           ! Current patch
      integer                     :: ico              ! Cohort counter
      integer                     :: k                ! Level counter
      integer                     :: k1               ! Level counter
      integer                     :: k2               ! Level counter
      integer                     :: klsl             ! Alias for rk4site%lsl
      integer                     :: kben             ! Alias for layer beneath bottom
      integer                     :: ksn              ! # of temporary water/snow layers
      integer                     :: nsoil            ! Short for csite%soil_text(k,ipa)
      real(kind=8)                :: snden            ! Snow/water density
      real(kind=8)                :: hflxsc           ! TSW -> canopy heat flux
      real(kind=8)                :: wflxsc           ! TSW -> canopy water flux
      real(kind=8)                :: qwflxsc          ! TSW -> canopy latent heat flux
      real(kind=8)                :: hflxgc           ! Soil -> canopy heat flux
      real(kind=8)                :: wflxgc           ! Soil -> canopy water flux
      real(kind=8)                :: qwflxgc          ! Soil -> canopy latent heat flux
      real(kind=8)                :: dewgnd           ! Dew/frost flux to ground
      real(kind=8)                :: qdewgnd          ! Dew/frost heat flux to ground
      real(kind=8)                :: ddewgnd          ! Dew/frost density flux to ground
      real(kind=8)                :: wshed_tot        ! Water shedding flux
      real(kind=8)                :: qwshed_tot       ! Energy flux due to water shedding
      real(kind=8)                :: dwshed_tot       ! Depth flux due to water shedding
      real(kind=8)                :: throughfall_tot  ! Water shedding flux
      real(kind=8)                :: qthroughfall_tot ! Energy flux due to water shedding
      real(kind=8)                :: dthroughfall_tot ! Depth flux due to water shedding
      real(kind=8)                :: wilting_factor   ! Wilting factor
      real(kind=8)                :: ext_weight       ! Layer weight for transpiration
      real(kind=8)                :: wloss            ! Water loss due to transpiration
      real(kind=8)                :: wloss_tot        ! Total water loss amongst cohorts
      real(kind=8)                :: wvlmeloss_tot    ! Total water loss amongst cohorts
      real(kind=8)                :: qloss            ! Energy loss due to transpiration
      real(kind=8)                :: qloss_tot        ! Total energy loss amongst cohorts
      real(kind=8)                :: qvlmeloss_tot    ! Total energy loss amongst cohorts
      real(kind=8)                :: infilt           ! Surface infiltration rate
      real(kind=8)                :: qinfilt          ! Surface infiltration heat rate
      real(kind=8)                :: surface_water    ! Temp. variable. Available liquid
      real(kind=8)                :: avg_th_cond      ! Mean thermal conductivity
      real(kind=8)                :: avg_hydcond      ! Mean thermal conductivity
      real(kind=8)                :: avail_h2o_int_i  ! 1/available water (lyr k1)
      real(kind=8)                :: wloss_tot_k1     ! Total water loss (lyr k1)
      real(kind=8)                :: wloss_tot_k2     ! Total water loss (lyr k2)
      real(kind=8)                :: uint_water_k1    ! Intensive Internal Energy (lyr k1)
      real(kind=8)                :: uint_water_k2    ! Intensive Internal Energy (lyr k2)
      real(kind=8)                :: rshort_c_soil    ! Committed shortwave (top soil lyr)
      real(kind=8)                :: rshort_c_sfcw    ! Committed shortwave (top TSW lyr)
      real(kind=8)                :: rlong_c_soil     ! Committed longwave  (top soil lyr)
      !------------------------------------------------------------------------------------!


      !----- Set the pointer to the current patch. ----------------------------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!


      !----- Copy the # of surface water/snow layers and bottom layer to aliases ----------!
      ksn  = initp%nlev_sfcwater
      klsl = rk4site%lsl
      kben = klsl - 1
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Find the committed shortwave and longwave radiation.  Temporary surface water !
      ! (TSW) is dynamic and may disappear in between two thermodynamic steps (dtlsm).     !
      ! In contrast, radiation is calculated only once every dtlsm.  In this case, we must !
      ! redirect the radiation (committed radiation) to some other place, to ensure that   !
      ! energy is conserved.                                                               !
      !                                                                                    !
      ! Longwave radiation.  Net absorption is always applied to the TSW top layer.        !
      !                      Committed radiation is only applied when all TSW layers       !
      !                      disappear.                                                    !
      ! Shortwave radiation. Net absorption is applied to all layers.  In case some layers !
      !                      disappear, we add the net absorption of the extinct layers to !
      !                      the top TSW layer.  In case all TSW layers disappear, we      !
      !                      apply the committed radiation to the top soil layer.          !
      !------------------------------------------------------------------------------------!
      select case (initp%nlev_sfcwater)
      case (0)
         !------ No TSW layers left. Committed radiation goes to soil. --------------------!
         rlong_c_soil     = dble(csite%rlong_s(ipa))
         rshort_c_soil    = 0.d0
         do k = 1,mzs
            rshort_c_soil = rshort_c_soil + dble(csite%rshort_s(k,ipa))
         end do
         rshort_c_sfcw    = 0.d0
         !---------------------------------------------------------------------------------!
      case default
         !---------------------------------------------------------------------------------!
         !    There is still at least one TSW layer. Committed radiation goes to the top   !
         ! existing layer.                                                                 !
         !---------------------------------------------------------------------------------!
         rlong_c_soil     = 0.d0
         rshort_c_soil    = 0.d0
         rshort_c_sfcw    = 0.d0
         do k = ksn+1,mzs
            rshort_c_sfcw = rshort_c_sfcw + dble(csite%rshort_s(k,ipa))
         end do
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Compute the following variables:                                               !
      !                                                                                    !
      ! TH_COND_S           -- thermal conductivity - soil                        [ W/m/K] !
      ! TH_COND_P           -- thermal conductivity - pounding water or snow pack [ W/m/K] !
      ! HYDCOND             -- hydraulic conductivity                             [   m/s] !
      ! PSIPLUSZ            -- the total potential (water + gravitational)        [     m] !
      ! DRYSOIL             -- flag that tells whether this layer is totally dry  [   T|F] !
      ! SATSOIL             -- flag that tells whether this layer is saturated    [   T|F] !
      !------------------------------------------------------------------------------------!
      do k = klsl, mzg
         nsoil                      = rk4site%ntext_soil(k)
         rk4aux(ibuff)%th_cond_s(k) = ( soil8(nsoil)%thcond0                               &
                                      + soil8(nsoil)%thcond1 * initp%soil_water(k) )       &
                                    / ( soil8(nsoil)%thcond2                               &
                                      + soil8(nsoil)%thcond3 * initp%soil_water(k) )


         !----- Find the correction for (partially) frozen soil layers. -------------------!
         rk4aux(ibuff)%hydcond (k) = hydr_conduct8(k,nsoil,initp%soil_water(k)             &
                                                  ,initp%soil_fracliq(k))
         !---------------------------------------------------------------------------------!


         rk4aux(ibuff)%psiplusz(k) = slzt8(k) + initp%soil_mstpot(k)
         rk4aux(ibuff)%drysoil (k) = (initp%soil_water(k) - soil8(nsoil)%soilcp)           &
                                   * initp%soil_fracliq(k)                        <= 0.d0
         rk4aux(ibuff)%satsoil (k) = initp%soil_water(k) >= soil8(nsoil)%slmsts
      end do
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Find the thermal conductivity of the temporary surface water/snow.            !
      !------------------------------------------------------------------------------------!
      do k = 1, ksn
         if (initp%sfcwater_depth(k) > 0.d0 .and. initp%sfcwater_mass(k) > 0.d0) then
            snden = initp%sfcwater_mass(k) / initp%sfcwater_depth(k)
            rk4aux(ibuff)%th_cond_p(k) = ss(1) * exp(ss(2) * initp%sfcwater_tempk(k))      &
                                * (ss(3) + snden * (ss(4) + snden * (ss(5) + snden*ss(6))))
         else
            rk4aux(ibuff)%th_cond_p(k) = 0.d0
         end if
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Compute the following variables:                                               !
      !                                                                                    !
      ! AVAIL_H2O_LYR       -- the available water factor for this layer          [ kg/m2] !
      ! AVAIL_H2O_INT       -- the integral of AVAIL_H2O down to the layer        [ kg/m2] !
      !                                                                                    !
      ! Both AVAIL_H2O_LYR and AVAIL_H2O_INT depend on the water limitation method.        !
      !------------------------------------------------------------------------------------!
      select case (h2o_plant_lim)
      case (0,1)
         !---------------------------------------------------------------------------------!
         !     The available water factor is the fraction of water mass above wilting      !
         ! point, scaled by the liquid fraction.                                           !
         !---------------------------------------------------------------------------------!
         do k = mzg, klsl, -1
            nsoil                   = rk4site%ntext_soil(k)

            !----- Find the available water factor for this layer. ------------------------!
            rk4aux(ibuff)%avail_h2o_lyr(k) =                                               &
                                      max(0.d0,(initp%soil_water(k)-soil8(nsoil)%soilwp))  &
                                    * initp%soil_fracliq(k) * wdns8 * dslz8(k)
            if (rk4aux(ibuff)%avail_h2o_lyr(k) < tiny_num8) then
               rk4aux(ibuff)%avail_h2o_lyr(k) = 0.d0
            end if
            !------------------------------------------------------------------------------!

            !----- Add the factor from this layer to the integral. ------------------------!
            rk4aux(ibuff)%avail_h2o_int(k) = rk4aux(ibuff)%avail_h2o_int(k+1)              &
                                           + rk4aux(ibuff)%avail_h2o_lyr(k)
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!

      case default
         !---------------------------------------------------------------------------------!
         !     The available water factor is the soil moisture at field capacity minus     !
         ! wilting, scaled by the wilting factor, defined as a function of soil potential. !
         !---------------------------------------------------------------------------------!
         do k = mzg, klsl, -1
            nsoil                   = rk4site%ntext_soil(k)

            !----- Find the available water factor for this layer. ------------------------!
            wilting_factor          = (rk4aux(ibuff)%psiplusz(k) - soil8(nsoil)%slpotwp)   &
                                    / (soil8(nsoil)%slpotfc - soil8(nsoil)%slpotwp)
            rk4aux(ibuff)%avail_h2o_lyr(k) = min( 1.d0, max( 0.d0, wilting_factor ) )      &
                                           * initp%soil_fracliq(k)                         &
                                           * ( soil8(nsoil)%sfldcap-soil8(nsoil)%soilwp )  &
                                           * wdns8 * dslz8(k)
            if (rk4aux(ibuff)%avail_h2o_lyr(k) < tiny_num8) then
               rk4aux(ibuff)%avail_h2o_lyr(k) = 0.d0
            end if
            !------------------------------------------------------------------------------!


            !----- Add the factor from this layer to the integral. ------------------------!
            rk4aux(ibuff)%avail_h2o_int(k) = rk4aux(ibuff)%avail_h2o_int(k+1)              &
                                           + rk4aux(ibuff)%avail_h2o_lyr(k)
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Find the boundary condition for total potential beneath the bottom layer.       !
      !------------------------------------------------------------------------------------!
      nsoil = rk4site%ntext_soil(klsl)
      select case (rk4site%isoilbc)
      case (0)
         !---------------------------------------------------------------------------------!
         !    Bedrock.  Make the potential exactly the same as the bottom layer, and the   !
         ! flux will be zero.                                                              !
         !---------------------------------------------------------------------------------!
         initp%soil_water       (kben) = initp%soil_water   (klsl)
         initp%soil_mstpot      (kben) = initp%soil_mstpot  (klsl)
         initp%soil_fracliq     (kben) = initp%soil_fracliq (klsl)
         rk4aux(ibuff)%th_cond_s(kben) = rk4aux(ibuff)%th_cond_s   (klsl)
         rk4aux(ibuff)%hydcond  (kben) = rk4aux(ibuff)%hydcond     (klsl)
         rk4aux(ibuff)%psiplusz (kben) = rk4aux(ibuff)%psiplusz    (klsl)
         rk4aux(ibuff)%drysoil  (kben) = .true.
         rk4aux(ibuff)%satsoil  (kben) = .true.
         !---------------------------------------------------------------------------------!

      case (1)
         !---------------------------------------------------------------------------------!
         !     Free drainage.   Make the water potential at the layer beneath to be at the !
         ! same soil moisture as the bottom layer.                                         !
         !---------------------------------------------------------------------------------!
         initp%soil_water       (kben) = initp%soil_water   (klsl)
         initp%soil_mstpot      (kben) = initp%soil_mstpot  (klsl)
         initp%soil_fracliq     (kben) = initp%soil_fracliq (klsl)
         rk4aux(ibuff)%th_cond_s(kben) = rk4aux(ibuff)%th_cond_s   (klsl)
         rk4aux(ibuff)%hydcond  (kben) = rk4aux(ibuff)%hydcond     (klsl)
         rk4aux(ibuff)%psiplusz (kben) = slzt8(kben) + initp%soil_mstpot(kben)
         rk4aux(ibuff)%drysoil  (kben) = .false.
         rk4aux(ibuff)%satsoil  (kben) = .false.
         !---------------------------------------------------------------------------------!

      case (2)
         !---------------------------------------------------------------------------------!
         !     Lateral drainage (or reduced drainage).  Find the equivalent depth of the   !
         ! layer beneath as a function of the slope (sldrain), and assume the soil         !
         ! moisture and matric potential to be the same as the bottom layer.  Notice that  !
         ! when sldrain is zero this becomes the flat bedrock condition, and when sldrain  !
         ! is 90 degrees, then it becomes free drainage.                                   !
         !---------------------------------------------------------------------------------!
         initp%soil_water       (kben) = initp%soil_water   (klsl)
         initp%soil_mstpot      (kben) = initp%soil_mstpot  (klsl)
         initp%soil_fracliq     (kben) = initp%soil_fracliq (klsl)
         rk4aux(ibuff)%th_cond_s(kben) = rk4aux(ibuff)%th_cond_s   (klsl)
         rk4aux(ibuff)%hydcond  (kben) = rk4aux(ibuff)%hydcond     (klsl)
         rk4aux(ibuff)%psiplusz (kben) = slzt8(klsl) - dslzt8(klsl) * sin_sldrain8         &
                                       + initp%soil_mstpot(kben)
         rk4aux(ibuff)%drysoil  (kben) = .false.
         rk4aux(ibuff)%satsoil  (kben) = .false.
         !---------------------------------------------------------------------------------!

      case (3)
         !---------------------------------------------------------------------------------!
         !     Aquifer. Soil in the layer beneath is always at bubbling point.             !
         !---------------------------------------------------------------------------------!
         initp%soil_water       (kben) = soil8(nsoil)%slmsts
         initp%soil_mstpot      (kben) = soil8(nsoil)%slpots
         initp%soil_fracliq     (kben) = initp%soil_fracliq (klsl)
         rk4aux(ibuff)%th_cond_s(kben) = ( soil8(nsoil)%thcond0                            &
                                         + soil8(nsoil)%thcond1 * initp%soil_water(kben) ) &
                                       / ( soil8(nsoil)%thcond2                            &
                                         + soil8(nsoil)%thcond3 * initp%soil_water(kben) )
         rk4aux(ibuff)%hydcond  (kben) = slcons18(kben,nsoil)
         rk4aux(ibuff)%psiplusz (kben) = slzt8(kben) + initp%soil_mstpot(kben)
         rk4aux(ibuff)%drysoil  (kben) = .false.
         rk4aux(ibuff)%satsoil  (kben) = .false.
         !---------------------------------------------------------------------------------!

      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Get derivatives of vegetation and canopy air space variables, plus some       !
      ! fluxes that will be used for soil top boundary conditions and for transpiration.   !
      !------------------------------------------------------------------------------------!
      call canopy_derivs_two(mzg,initp,dinitp,csite,ipa,ibuff,hflxsc,wflxsc,qwflxsc,hflxgc &
                            ,wflxgc,qwflxgc,dewgnd,qdewgnd,ddewgnd,throughfall_tot         &
                            ,qthroughfall_tot,dthroughfall_tot,wshed_tot,qwshed_tot        &
                            ,dwshed_tot,dt,is_hybrid)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Calculate the sensible heat fluxes find soil and sfcwater internal sensible    !
      ! heat fluxes (h_flux_g and h_flux_s) [W/m2].  The convention is that fluxes are     !
      ! positive when they are upwards.                                                    !
      !     The average thermal conductivity is found using a log-linear interpolation     !
      ! between the mid-points of the consecutive layers.                                  !
      !------------------------------------------------------------------------------------!
      do k = klsl+1, mzg
         avg_th_cond                 =  rk4aux(ibuff)%th_cond_s(k-1)                       &
                                     *  ( rk4aux(ibuff)%th_cond_s(k)                       &
                                        / rk4aux(ibuff)%th_cond_s(k-1) )                   &
                                     ** ( dslz8(k-1) / ( dslz8(k-1) + dslz8(k) ) )
         rk4aux(ibuff)%h_flux_g(k)   = - avg_th_cond                                       &
                                       * (initp%soil_tempk(k) - initp%soil_tempk(k-1))     &
                                       * dslzti8(k)
         !------ Diagnostic sensible heat flux. -------------------------------------------!
         if (fast_diagnostics .or. print_detailed) then
            dinitp%avg_sensible_gg(k-1) = rk4aux(ibuff)%h_flux_g(k)
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If temporary water/snow layers exist, we compute them now.  Because the layers !
      ! may not cover the full area, we scale the fluxes by the fraction.                  !
      !------------------------------------------------------------------------------------!
      if (ksn >= 1) then
         !---------------------------------------------------------------------------------!
         !     The first layer is the interface between soil and TSW.  We account for the  !
         ! fluxes twice.                                                                   !
         !---------------------------------------------------------------------------------!
         avg_th_cond                   =  rk4aux(ibuff)%th_cond_s(mzg)                     &
                                       *  ( rk4aux(ibuff)%th_cond_p(1)                     &
                                          / rk4aux(ibuff)%th_cond_s(mzg) )                 &
                                       ** ( dslz8(mzg)                                     &
                                          / (initp%sfcwater_depth(1)+ dslz8(mzg)))
         rk4aux(ibuff)%h_flux_g(mzg+1) = - avg_th_cond                                     &
                                       * (initp%sfcwater_tempk(1) - initp%soil_tempk(mzg)) &
                                       / (5.d-1 * initp%sfcwater_depth(1) - slzt8(mzg) )
         rk4aux(ibuff)%h_flux_s(1)     = rk4aux(ibuff)%h_flux_g(mzg+1)
         do k = 2,ksn
            avg_th_cond                =  rk4aux(ibuff)%th_cond_p(k-1)                     &
                                       *  ( rk4aux(ibuff)%th_cond_p(k)                     &
                                          / rk4aux(ibuff)%th_cond_p(k-1))                  &
                                       ** ( initp%sfcwater_depth(k-1)                      &
                                          / ( initp%sfcwater_depth(k-1)                    &
                                            + initp%sfcwater_depth(k) ) )
            rk4aux(ibuff)%h_flux_s(k)  = - 2.d0 * avg_th_cond                              &
                                    * ( initp%sfcwater_tempk(k)-initp%sfcwater_tempk(k-1)) &
                                    / ( initp%sfcwater_depth(k)+initp%sfcwater_depth(k-1))
         end do
         !---------------------------------------------------------------------------------!

         ! This is just sensible heat flux loss from top surface layer
         ! The layer's energy budget (shown below) will include the other terms
         ! Convention is positive up
         rk4aux(ibuff)%h_flux_s        (ksn+1) = rk4aux(ibuff)%h_flux_s(ksn+1) + hflxsc

      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Add the irradiance and canopy fluxes.                                         !
      !------------------------------------------------------------------------------------!
      if (fast_diagnostics .or. print_detailed) then
         dinitp%avg_sensible_gg(mzg)   = hflxgc + qwflxgc - dble(csite%rlong_g(ipa))       &
                                       - dble(csite%rshort_g(ipa)) - rshort_c_soil         &
                                       - rlong_c_soil
      end if
      rk4aux(ibuff)%h_flux_g(mzg+1) = rk4aux(ibuff)%h_flux_g(mzg+1)                        &
                                    + hflxgc + qwflxgc - dble(csite%rlong_g(ipa))          &
                                    - dble(csite%rshort_g(ipa)) - rshort_c_soil            &
                                    - rlong_c_soil
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Update soil U values [J/m2/s] from sensible heat, upward water vapor (latent    !
      ! heat) and longwave fluxes. This excludes effects of dew/frost formation,           !
      ! precipitation, shedding, and percolation.                                          !
      !------------------------------------------------------------------------------------!
      do k = klsl,mzg
         dinitp%soil_energy(k) = dslzi8(k) * ( rk4aux(ibuff)%h_flux_g(k)                   &
                                             - rk4aux(ibuff)%h_flux_g(k+1) )
      end do
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! Update surface water U values [J/m2/s] from sensible heat and shortwave            !
      ! fluxes.  This excludes effects of dew/frost, latent heat flux, precipitation,      !
      ! shedding and percolation (and any mass fluxes). Right now, thermal radiation only  !
      ! affects the top layer.                                                             !
      !------------------------------------------------------------------------------------!
      do k = 1,ksn
        dinitp%sfcwater_energy(k) = rk4aux(ibuff)%h_flux_s(k)                              &
                                  - rk4aux(ibuff)%h_flux_s(k+1)                            &
                                  + dble(csite%rshort_s(k,ipa))
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Calculate the fluxes of water with their associated heat fluxes. Update top    !
      ! soil or snow moisture from evaporation only.                                       !
      !------------------------------------------------------------------------------------!
      if ( ksn > 0 ) then
         dinitp%sfcwater_mass  (ksn) =  dewgnd +  wshed_tot +  throughfall_tot -  wflxsc
         dinitp%sfcwater_energy(ksn) = dinitp%sfcwater_energy(ksn)                         &
                                     + dble(csite%rlong_s(ipa)) + rshort_c_sfcw            &
                                     + qdewgnd + qwshed_tot + qthroughfall_tot - qwflxsc
         dinitp%sfcwater_depth (ksn) = ddewgnd + dwshed_tot + dthroughfall_tot
      else
         dinitp%virtual_water        =  dewgnd +  wshed_tot +  throughfall_tot -  wflxsc
         dinitp%virtual_energy       = qdewgnd + qwshed_tot + qthroughfall_tot - qwflxsc
         dinitp%virtual_depth        = ddewgnd + dwshed_tot + dthroughfall_tot
      end if
      !------------------------------------------------------------------------------------!



      !------ Diagnostic variable for water flux, bypass the virtual/sfcw layers. ---------!
      if (fast_diagnostics .or. print_detailed) then
         dinitp%avg_smoist_gg(mzg) = rk4aux(ibuff)%w_flux_g(mzg+1)                         &
                                   + dewgnd +  wshed_tot +  throughfall_tot -  wflxsc      &
                                   - wflxgc
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Find amount of water transferred between soil layers (w_flux) [m] modulated by !
      ! the liquid water fraction.                                                         !
      !------------------------------------------------------------------------------------!
      rk4aux(ibuff)%w_flux_g(mzg+1) = wflxgc * wdnsi8 ! now in m/s
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     Alternate surface infiltration (MCD) based on surface conductivity not         !
      ! capacity.                                                                          !
      !------------------------------------------------------------------------------------!
      if (infiltration_method /= 0) then
         call fatal_error ('Running alt infiltation when we shouldn''t be'                 &
                          ,'leaftw_derivs','rk4_derivs.F90')

         if (initp%virtual_water /= 0.d0) then  !!process "virtural water" pool
            nsoil = rk4site%ntext_soil(mzg)
            select case (trim(soil8(nsoil)%method))
            case ('BDRK')
               continue
            case default
               infilt = - dslzi8(mzg) * 5.d-1                                              &
                        * hydr_conduct8(mzg,nsoil,initp%soil_water(mzg)                    &
                                       ,initp%soil_fracliq(mzg))                           &
                        * (rk4aux(ibuff)%psiplusz(mzg)-initp%virtual_water/2.d3)           & ! diff. in pot.
                        * 5.d-1 * (initp%soil_fracliq(mzg)+ initp%virtual_fracliq)         ! ! mean liquid fraction
               qinfilt = infilt * wdns8 * tl2uint8(initp%virtual_tempk,1.d0)
               !----- Adjust other rates accordingly --------------------------------------!
               rk4aux(ibuff)%w_flux_g (mzg+1) = rk4aux(ibuff)%w_flux_g(mzg+1)  + infilt
               rk4aux(ibuff)%qw_flux_g(mzg+1) = rk4aux(ibuff)%qw_flux_g(mzg+1) + qinfilt
               dinitp%virtual_water    = dinitp%virtual_water    - infilt*wdns8
               dinitp%virtual_energy   = dinitp%virtual_energy   - qinfilt
            end select
         end if  !! end virtual water pool
         if (initp%nlev_sfcwater >= 1) then !----- Process "snow" water pool --------------!
            surface_water = initp%sfcwater_mass(1)*initp%sfcwater_fracliq(1)*wdnsi8 !(m/m2)
            nsoil = rk4site%ntext_soil(mzg)
            select case (trim(soil8(nsoil)%method))
            case ('BDRK')
               continue
            case default
               !----- Calculate infiltration rate (m/s) -----------------------------------!
               infilt = - dslzi8(mzg) * 5.d-1                                              &
                        * hydr_conduct8(mzg,nsoil,initp%soil_water(mzg)                    &
                                       ,initp%soil_fracliq(mzg))                           &
                        * (rk4aux(ibuff)%psiplusz(mzg) - surface_water/2.d0)               & ! difference in potentials
                        * 5.d-1 * (initp%soil_fracliq(mzg) + initp%sfcwater_fracliq(1))    ! ! mean liquid fraction
               qinfilt = infilt * wdns8 * tl2uint8(initp%sfcwater_tempk(1),1.d0)
               !----- Adjust other rates accordingly --------------------------------------!
               rk4aux(ibuff)%w_flux_g(mzg+1)    = rk4aux(ibuff)%w_flux_g(mzg+1)  + infilt
               rk4aux(ibuff)%qw_flux_g(mzg+1)   = rk4aux(ibuff)%qw_flux_g(mzg+1) + qinfilt
               dinitp%sfcwater_mass(1)   = dinitp%sfcwater_mass(1)   - infilt*wdns8
               dinitp%sfcwater_energy(1) = dinitp%sfcwater_energy(1) - qinfilt
               dinitp%sfcwater_depth(1)  = dinitp%sfcwater_depth(1)  - infilt
            end select
         end if  ! End snow water pool
      end if  !! End alternate infiltration
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Here we solve for all layers including the bottommost one, which is now        !
      ! prepared according to the chosen boundary condition.  In this block, units are:    !
      ! W_Flux  = m/s                                                                      !
      ! QW_Flux = W/m2                                                                     !
      !------------------------------------------------------------------------------------!
      do k = klsl, mzg
         nsoil = rk4site%ntext_soil(k)
         select case (trim(soil8(nsoil)%method))
         case ('BDRK')
            rk4aux(ibuff)%w_flux_g(k) = 0.d0
         case default

            !----- Log-linear interpolation of hydraulic conductivity to layer interface. -!
            avg_hydcond =  rk4aux(ibuff)%hydcond(k-1)                                      &
                        *  ( rk4aux(ibuff)%hydcond(k) / rk4aux(ibuff)%hydcond(k-1) )       &
                        ** ( dslz8(k-1) / ( dslz8(k-1) + dslz8(k) ) )
            !------------------------------------------------------------------------------!



            !----- Find the potential flux. -----------------------------------------------!
            rk4aux(ibuff)%w_flux_g(k) = - avg_hydcond * ( rk4aux(ibuff)%psiplusz(k)        &
                                                        - rk4aux(ibuff)%psiplusz(k-1) )    &
                                      * dslzti8(k)
            !------------------------------------------------------------------------------!



            !----- Limit water transfers to prevent over-saturation and over-depletion. ---!
            if ( rk4aux(ibuff)%w_flux_g(k) >= 0. .and.                                     &
                (rk4aux(ibuff)%drysoil(k-1) .or. rk4aux(ibuff)%satsoil(k)) )    then
               rk4aux(ibuff)%w_flux_g(k) = 0.d0

            elseif( rk4aux(ibuff)%w_flux_g(k) < 0. .and.                                   &
                   (rk4aux(ibuff)%satsoil(k-1) .or. rk4aux(ibuff)%drysoil(k)) ) then
               rk4aux(ibuff)%w_flux_g(k) = 0.d0

            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Find the internal energy flux associated with the water flux.  This is     !
         ! sign-dependent because the temperature must be the source temperature.          !
         !---------------------------------------------------------------------------------!
         if (rk4aux(ibuff)%w_flux_g(k) > 0) then
            rk4aux(ibuff)%qw_flux_g(k) = rk4aux(ibuff)%w_flux_g(k) * wdns8                 &
                                       * tl2uint8(initp%soil_tempk(k-1),1.d0)
         else
            rk4aux(ibuff)%qw_flux_g(k) = rk4aux(ibuff)%w_flux_g(k) * wdns8                 &
                                       * tl2uint8(initp%soil_tempk(k)  ,1.d0)
         end if
         !---------------------------------------------------------------------------------!


         !----- Save the moisture flux in kg/m2/s. ----------------------------------------!
         if ( (fast_diagnostics .or. print_detailed) .and. (k /= 1) ) then
            dinitp%avg_smoist_gg(k-1) = rk4aux(ibuff)%w_flux_g(k) * wdns8 ! Diagnostic
         end if
         !---------------------------------------------------------------------------------!
      end do

      !------------------------------------------------------------------------------------!
      !     Find the drainage flux.  This is simply the flux at the bottom of the bottom-  !
      ! most layer.  Units are kg/m2/s for water flux and W/m2 for the associated internal !
      ! energy leak.  Notice that for the aquifer case we may actually have negative       !
      ! drainage, but that shouldn't affect the budget in any way (except that we are add- !
      ! ing water to the system).                                                          !
      !------------------------------------------------------------------------------------!
      if (fast_diagnostics .or. print_detailed) then
         dinitp%avg_drainage  = - rk4aux(ibuff)%w_flux_g (klsl) * wdns8
         dinitp%avg_qdrainage = - rk4aux(ibuff)%qw_flux_g(klsl)
      end if
      !----- Copy the variables to the budget arrays. -------------------------------------!
      if (checkbudget) then
         dinitp%wbudget_loss2drainage = - rk4aux(ibuff)%w_flux_g (klsl) * wdns8
         dinitp%ebudget_loss2drainage = - rk4aux(ibuff)%qw_flux_g(klsl)

         dinitp%wbudget_storage = dinitp%wbudget_storage - dinitp%wbudget_loss2drainage
         dinitp%ebudget_storage = dinitp%ebudget_storage - dinitp%ebudget_loss2drainage
      end if
      !------------------------------------------------------------------------------------!




      !----- Finally, update soil moisture and soil energy. -------------------------------!
      do k = klsl,mzg
         dinitp%soil_water(k)  = dinitp%soil_water(k)                                      &
                               + dslzi8(k) * (  rk4aux(ibuff)%w_flux_g(k)                  &
                                             -  rk4aux(ibuff)%w_flux_g(k+1) )
         dinitp%soil_energy(k) = dinitp%soil_energy(k)                                     &
                               + dslzi8(k) * ( rk4aux(ibuff)%qw_flux_g(k)                  &
                                             - rk4aux(ibuff)%qw_flux_g(k+1) )
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Check whether this simulation solves plant hydraulics.                        !
      !------------------------------------------------------------------------------------!
      select case (plant_hydro_scheme)
      case (0) ! No plant hydraulics
         !---- Update soil moisture and energy from transpiration/root uptake. ------------!
         if (rk4aux(ibuff)%any_resolvable) then
            !------------------------------------------------------------------------------!
            !    Loop over extracted water.                                                !
            !------------------------------------------------------------------------------!
            k1_transp_loop: do k1 = klsl, mzg
               !---------------------------------------------------------------------------!
               !     Transpiration happens only when there is some water left down to this !
               ! layer.                                                                    !
               !---------------------------------------------------------------------------!
               if (rk4aux(ibuff)%avail_h2o_int(k1) < tiny_num8) cycle k1_transp_loop
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !     Find inverse of available water (so we perform division only once).   !
               !---------------------------------------------------------------------------!
               avail_h2o_int_i = 1.d0 / rk4aux(ibuff)%avail_h2o_int(k1)
               !---------------------------------------------------------------------------!


               !-------- Integrate the total to be removed from this layer. ---------------!
               wloss_tot_k1 = 0.d0
               do ico=1,cpatch%ncohorts
                  wloss_tot_k1 = wloss_tot_k1 + rk4aux(ibuff)%extracted_water(ico,k1)
               end do
               !---------------------------------------------------------------------------!


               !----- Inner loop. ---------------------------------------------------------!
               k2_transp_loop: do k2=k1,mzg
                  !------------------------------------------------------------------------!
                  !     Skip calculation if no water is available in this layer.           !
                  !------------------------------------------------------------------------!
                  if ( trim(soil8(rk4site%ntext_soil(k2))%method) == 'BDRK'    .or.        &
                       rk4aux(ibuff)%avail_h2o_lyr(k2)            <  tiny_num8      ) then
                     cycle k2_transp_loop
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Find the contribution of layer k2 for the transpiration from        !
                  ! cohorts that reach layer k1.                                           !
                  !------------------------------------------------------------------------!
                  ext_weight = rk4aux(ibuff)%avail_h2o_lyr(k2) * avail_h2o_int_i
                  !------------------------------------------------------------------------!


                  !----- Save internal energy of water in layer k2. -----------------------!
                  uint_water_k2 = tl2uint8(initp%soil_tempk(k2),1.d0)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Find the loss of water from layer k2 due to cohorts that reach at   !
                  ! least layer k1.  Here we convert extracted water which is in kg/m2/s   !
                  ! (tloss) to m3/m3/s (wloss).  Also, find the internal energy loss       !
                  ! (qwloss) associated with the water loss.  Since plants can extract     !
                  ! liquid water only, the internal energy is assumed to be entirely in    !
                  ! liquid phase.  Because the actual conversion from liquid phase to      !
                  ! vapour happens at the leaf level, the internal energy must stay with   !
                  ! the leaves so energy is preserved.                                     !
                  !------------------------------------------------------------------------!
                  do ico=1,cpatch%ncohorts
                     !----- Find the soil water loss associated with this cohort. ---------!
                     wloss         = rk4aux(ibuff)%extracted_water(ico,k1) * ext_weight
                     qloss         = wloss * uint_water_k2
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !      Add the internal energy to the cohort.  This energy will be    !
                     ! eventually lost to the canopy air space because of transpiration,   !
                     ! but we will do it in two steps so we ensure energy is conserved.    !
                     !---------------------------------------------------------------------!
                     dinitp%leaf_energy(ico) = dinitp%leaf_energy(ico) + qloss
                     dinitp%veg_energy (ico) = dinitp%veg_energy (ico) + qloss
                     initp%hflx_lrsti  (ico) = initp%hflx_lrsti  (ico) + qloss
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!

                  !--------- Derive the total to be removed from tthis layer --------------!
                  wloss_tot     = wloss_tot_k1 * ext_weight
                  wloss_tot_k2  = wloss_tot    * dslzi8(k2)
                  wvlmeloss_tot = wloss_tot_k2 * wdnsi8
                  qvlmeloss_tot = wloss_tot_k2 * uint_water_k2
                  !------------------------------------------------------------------------!


                  !----- Update derivatives of water, energy, and transpiration. ----------!
                  dinitp%soil_water   (k2) = dinitp%soil_water   (k2) - wvlmeloss_tot
                  dinitp%soil_energy  (k2) = dinitp%soil_energy  (k2) - qvlmeloss_tot
                  if (fast_diagnostics .or. print_detailed) then
                     dinitp%avg_transloss(k2) = dinitp%avg_transloss(k2) - wloss_tot
                  end if
                  !------------------------------------------------------------------------!
               end do k2_transp_loop
               !---------------------------------------------------------------------------!
            end do k1_transp_loop
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      case default
         !---------------------------------------------------------------------------------!
         !     Track plant hydraulics.  In this case, we must always update soil water and !
         ! soil energy, even if the cohort is not resolvable.                              !
         !                                                                                 !
         !     Leaf and wood internal water are updated outside of the integrator, in      !
         ! rk4_driver.f90                                                                  !
         !---------------------------------------------------------------------------------!
         k1_transh_loop: do k1 = klsl, mzg    ! loop over soil layers
            wloss_tot      = 0.d0
            qloss_tot      = 0.d0
            uint_water_k1  = tl2uint8(initp%soil_tempk(k1),1.d0)


            !------------------------------------------------------------------------------!
            !  MLO -> XX.  I added this if to bypass contributions from this layer to      !
            !              transpiration when the soil is completely desiccated.  Do you   !
            !              foresee any problems?                                           !
            !------------------------------------------------------------------------------!
            if (rk4aux(ibuff)%drysoil(k1)) then
               cycle k1_transh_loop
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Integrate the total to be removed from this layer.  Add water and       !
            ! internal energy to each cohort, so water and energy are conserved.  The      !
            ! internal energy should go to wood.                                           !
            !                                                                              !
            ! MLO -> XX                                                                    !
            ! How does this work in case we are solving grasses?  Specifically, new        !
            ! grasses do not have any wood.                                                !
            ! XX -> MLO                                                                    !
            ! Update leaf_water_im2 and leaf_energy pool instead for new grasses? (or any  !
            ! condition with resolvable leaf and unresolvable wood) TODO                   !
            !------------------------------------------------------------------------------!
            do ico=1,cpatch%ncohorts
               !----- Find the soil water loss associated with this cohort. ---------------!
               wloss = dble(cpatch%wflux_gw_layer(k1,ico))                                 &
                     * dble(cpatch%nplant(ico))            ! kg/m2g/s
               if (wloss >= 0.d0) then
                  qloss = wloss * uint_water_k1            !  J/m2g/s
               else
                  qloss = wloss * tl2uint8(initp%wood_temp(ico),1.d0)
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Add the water and internal energy to the cohort.  This energy will   !
               ! be eventually lost to the canopy air space because of transpiration, but  !
               ! we will do it in two steps so we ensure energy is conserved.              !
               !---------------------------------------------------------------------------!
               dinitp%wood_water_im2(ico) = dinitp%wood_water_im2(ico) + wloss
               dinitp%veg_water_im2 (ico) = dinitp%veg_water_im2 (ico) + wloss
               dinitp%veg_energy    (ico) = dinitp%veg_energy    (ico) + qloss
               dinitp%wood_energy   (ico) = dinitp%wood_energy   (ico) + qloss
               initp%hflx_lrsti     (ico) = initp%hflx_lrsti     (ico) + qloss
               !---------------------------------------------------------------------------!


               !---- Count the total loss from the layer. ---------------------------------!
               wloss_tot = wloss_tot + wloss
               qloss_tot = qloss_tot + qloss
               !---------------------------------------------------------------------------!


               !----- Averaged fluxes. ----------------------------------------------------!
               if (fast_diagnostics .or. print_detailed) then
                  dinitp%avg_wflux_gw         (ico) = dinitp%avg_wflux_gw(ico) + wloss
                  dinitp%avg_wflux_gw_layer(k1,ico) = wloss
               end if
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!


            !---- Convert ground-wood water flow from kg/m2/s to m3/m3/s. -----------------!
            wvlmeloss_tot  = wloss_tot * wdnsi8 * dslzi8(k1)
            qvlmeloss_tot  = qloss_tot          * dslzi8(k1)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Update derivatives of water, energy, and transpiration.  Here notice     !
            ! that the avg_transloss actually represents total soil water loss due to      !
            ! plant uptake, not transpiration.                                             !
            !------------------------------------------------------------------------------!
            dinitp%soil_water   (k1) = dinitp%soil_water   (k1) - wvlmeloss_tot
            dinitp%soil_energy  (k1) = dinitp%soil_energy  (k1) - qvlmeloss_tot
            if (fast_diagnostics .or. print_detailed) then
               dinitp%avg_transloss(k1) = dinitp%avg_transloss(k1) - wloss_tot
            end if
            !------------------------------------------------------------------------------!
         end do k1_transh_loop
         !---------------------------------------------------------------------------------!

      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine leaftw_derivs
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine canopy_derivs_two(mzg,initp,dinitp,csite,ipa,ibuff,hflxsc,wflxsc,qwflxsc     &
                               ,hflxgc,wflxgc,qwflxgc,dewgndflx,qdewgndflx,ddewgndflx      &
                               ,throughfall_tot,qthroughfall_tot,dthroughfall_tot          &
                               ,wshed_tot,qwshed_tot,dwshed_tot,dt,is_hybrid)
      use rk4_coms              , only : rk4patchtype         & ! Structure
                                       , rk4site              & ! intent(in)
                                       , rk4aux               & ! intent(inout)
                                       , rk4eps               & ! intent(in)
                                       , effarea_heat         & ! intent(in)
                                       , effarea_evap         & ! intent(in)
                                       , effarea_transp       & ! intent(in)
                                       , tiny_offset          & ! intent(in)
                                       , rk4leaf_drywhc       & ! intent(in)
                                       , rk4leaf_maxwhc       & ! intent(in)
                                       , checkbudget          & ! intent(in)
                                       , print_detailed       & ! intent(in)
                                       , supersat_ok          & ! intent(in)
                                       , leaf_intercept       ! ! intent(in)
      use ed_state_vars         , only : sitetype             & ! Structure
                                       , patchtype            & ! Structure
                                       , polygontype          ! ! Structure
      use consts_coms           , only : twothirds8           & ! intent(in)
                                       , day_sec8             & ! intent(in)
                                       , grav8                & ! intent(in)
                                       , umol_2_kgC8          & ! intent(in)
                                       , pi18                 & ! intent(in)
                                       , halfpi8              & ! intent(in)
                                       , mmdry8               & ! intent(in)
                                       , mmdryi8              & ! intent(in)
                                       , wdns8                & ! intent(in)
                                       , wdnsi8               & ! intent(in)
                                       , fdnsi8               & ! intent(in)
                                       , t3ple8               & ! intent(in)
                                       , cpdry8               & ! intent(in)
                                       , cph2o8               & ! intent(in)
                                       , epi8                 & ! intent(in)
                                       , tsupercool_vap8      & ! intent(in)
                                       , huge_num8            ! ! intent(in)
      use soil_coms             , only : dslzi8               ! ! intent(in)
      use therm_lib8            , only : qslif8               & ! function
                                       , tq2enthalpy8         & ! function
                                       , tl2uint8             ! ! function
      use ed_misc_coms          , only : fast_diagnostics     ! ! intent(in)
      use canopy_struct_dynamics, only : vertical_vel_flux8   ! ! function
      use budget_utils          , only : compute_netrad8      ! ! function
      use physiology_coms       , only : plant_hydro_scheme   ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)     , target      :: csite            ! Current site
      type(rk4patchtype) , target      :: initp            ! RK4 structure, state vars
      type(rk4patchtype) , target      :: dinitp           ! RK4 structure, derivatives
      integer            , intent(in)  :: ipa              ! Current patch ID
      integer            , intent(in)  :: ibuff            ! Multithread ID
      integer            , intent(in)  :: mzg              ! Number of soil layers
      real(kind=8)       , intent(in)  :: dt               ! Timestep
      logical            , intent(in)  :: is_hybrid        ! Is this a hybrid call?
      real(kind=8)       , intent(out) :: hflxsc           ! Ground->canopy sens. heat flux
      real(kind=8)       , intent(out) :: wflxsc           ! Ground->canopy water flux
      real(kind=8)       , intent(out) :: qwflxsc          ! Ground->canopy latent heat flux
      real(kind=8)       , intent(out) :: hflxgc           ! Ground->canopy sens. heat flux
      real(kind=8)       , intent(out) :: wflxgc           ! Ground->canopy water flux
      real(kind=8)       , intent(out) :: qwflxgc          ! Ground->canopy latent heat flux
      real(kind=8)       , intent(out) :: dewgndflx        ! Dew/frost water flux
      real(kind=8)       , intent(out) :: qdewgndflx       ! Dew/frost heat flux
      real(kind=8)       , intent(out) :: ddewgndflx       ! Dew/frost density
      real(kind=8)       , intent(out) :: throughfall_tot  ! Throughfall rate
      real(kind=8)       , intent(out) :: qthroughfall_tot ! Throughfall energy flux
      real(kind=8)       , intent(out) :: dthroughfall_tot ! Throughfall depth flux
      real(kind=8)       , intent(out) :: wshed_tot        ! Water shed from leaves
      real(kind=8)       , intent(out) :: qwshed_tot       ! Internal energy of water shed
      real(kind=8)       , intent(out) :: dwshed_tot       ! Depth of water shed
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer     :: cpatch            ! Current patch
      logical                      :: is_dew_cp         ! Test whether to add dew to TSW
      logical                      :: is_dew_cs         ! Test whether to add dew to soil
      logical                      :: lwater_fine       ! Enough leaf water for transpir.
      logical                      :: wwater_fine       ! Enough wood water for transpir.
      logical                      :: transp_fine       ! Test whether to allow transpir.
      integer                      :: ico               ! Current cohort ID
      integer                      :: ksn               ! Number of TSW layers
      integer                      :: ipft              ! Shortcut for PFT type
      integer                      :: kroot             ! Level of the bottom of root is
      real(kind=8)                 :: lwater_im2_lwr    ! Minimum leaf water content
      real(kind=8)                 :: wwater_im2_lwr    ! Minimum wood water content
      real(kind=8)                 :: closedcan_frac    ! total fractional canopy coverage
      real(kind=8)                 :: transp            ! Cohort transpiration
      real(kind=8)                 :: wflux_wl          ! Wood-leaf water flow (sapflow)
      real(kind=8)                 :: cflxac            ! Atm->canopy carbon flux
      real(kind=8)                 :: wflxac            ! Atm->canopy water flux
      real(kind=8)                 :: hflxac            ! Atm->canopy sensible heat flux
      real(kind=8)                 :: eflxac            ! Atm->canopy Eq. Pot. temp flux
      real(kind=8)                 :: wflxlc_try        ! Intended flux leaf sfc -> canopy
      real(kind=8)                 :: wflxwc_try        ! Intended flux wood sfc -> canopy
      real(kind=8)                 :: shv_gradient      ! Term for psi_open/psi_closed
      real(kind=8)                 :: gleaf_open        ! Net leaf conductance (open)
      real(kind=8)                 :: gleaf_closed      ! Net leaf conductance (closed)
      real(kind=8)                 :: hflxlc            ! Leaf->canopy heat flux
      real(kind=8)                 :: hflxwc            ! Wood->canopy heat flux
      real(kind=8)                 :: sigmaw            !
      real(kind=8)                 :: wflxlc            ! Leaf sfc -> canopy water flux
      real(kind=8)                 :: wflxwc            ! Wood sfc -> canopy water flux
      real(kind=8)                 :: wshed             ! Water shed from leaves
      real(kind=8)                 :: qwshed            ! Internal energy of water shed
      real(kind=8)                 :: dwshed            ! Depth of water shed
      real(kind=8)                 :: throughfall       ! Extra throughfall from full coh.
      real(kind=8)                 :: qthroughfall      ! Its internal energy
      real(kind=8)                 :: dthroughfall      ! Its depth
      real(kind=8)                 :: taii              !
      real(kind=8)                 :: wood_shv          ! Sat. sp. hum. at wood sfc.
      real(kind=8)                 :: hflxlc_tot        ! Total leaf -> CAS heat flux
      real(kind=8)                 :: hflxwc_tot        ! Total wood -> CAS heat flux
      real(kind=8)                 :: transp_tot        ! Total transpiration (water)
      real(kind=8)                 :: qtransp_tot       ! Total transpiration (energy)
      real(kind=8)                 :: nee_tot           ! Total NEE (source > 0; sink < 0)
      real(kind=8)                 :: wflxlc_tot        ! Leaf -> CAS evaporation (water)
      real(kind=8)                 :: wflxwc_tot        ! Wood -> CAS evaporation (water)
      real(kind=8)                 :: qwflxlc_tot       ! Leaf -> CAS evaporation (energy)
      real(kind=8)                 :: qwflxwc_tot       ! Wood -> CAS evaporation (energy)
      real(kind=8)                 :: rhos_ustar        !
      real(kind=8)                 :: dmol_ustar        !
      real(kind=8)                 :: min_leaf_water    !
      real(kind=8)                 :: max_leaf_water    !
      real(kind=8)                 :: min_wood_water    !
      real(kind=8)                 :: max_wood_water    !
      real(kind=8)                 :: dew_now           !
      real(kind=8)                 :: sfcwater_ssh      ! Specific humidity at top layer
      real(kind=8)                 :: intercepted_max   ! Pot. interecepted rainfall
      real(kind=8)                 :: qintercepted_max  ! Int. energy of pot. intercept.
      real(kind=8)                 :: dintercepted_max  ! Depth of pot. interception
      real(kind=8)                 :: intercepted_tot   ! Actual intercepted rainfall
      real(kind=8)                 :: qintercepted_tot  ! Int. energy of act. intercept.
      real(kind=8)                 :: leaf_intercepted  ! Leaf interception
      real(kind=8)                 :: leaf_qintercepted ! Int. energy of leaf intercept.
      real(kind=8)                 :: wood_intercepted  ! Leaf interception
      real(kind=8)                 :: wood_qintercepted ! Int. energy of leaf intercept.
      real(kind=8)                 :: qwflxlc           ! Leaf -> CAS evaporation (energy)
      real(kind=8)                 :: qwflxwc           ! Wood -> CAS evaporation (energy)
      real(kind=8)                 :: qtransp           ! Transpiration (energy)
      real(kind=8)                 :: qwflux_wl         ! Wood-leaf sapflow (energy)
      real(kind=8)                 :: flux_area         ! Area between canopy and plant
      real(kind=8)                 :: a,b,c0            ! Temporary variables for solving
                                                        ! the CO2 ODE
      real(kind=8)                 :: max_dwdt          ! Used for capping leaf evap
      real(kind=8)                 :: op_rk4eps         ! Buffer for internal water 
      !----- Functions --------------------------------------------------------------------!
      real(kind=4), external           :: sngloff           ! Safe dble 2 single precision
      !------------------------------------------------------------------------------------!

      !----- First step, we assign the pointer for the current patch. ---------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Computing the fluxes from atmosphere to canopy.                                 !
      !------------------------------------------------------------------------------------!
      rhos_ustar = initp%can_rhos * initp%ustar        ! Aux. variable
      dmol_ustar = initp%can_dmol * initp%ustar        ! Aux. variable
      wflxac     = rhos_ustar     * initp%qstar        ! Water flux     [  kg/m2/s]
      eflxac     = rhos_ustar     * initp%estar        ! Enthalpy flux  [   J/m2/s]
      cflxac     = dmol_ustar     * initp%cstar        ! CO2            [umol/m2/s]
      !------ Sensible heat flux. ---------------------------------------------------------!
      hflxac     = eflxac + wflxac * cph2o8 * tsupercool_vap8
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Calculate fraction of closed canopy.                                           !
      !------------------------------------------------------------------------------------!
      closedcan_frac = max(0.d0,min(1.d0,1.d0-initp%opencan_frac))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Set all heat fluxes to zero, we only add them if they occur and are allowed.   !
      !------------------------------------------------------------------------------------!
      hflxsc         = 0.d0
      hflxgc         = 0.d0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Here we will determine the initial guess for throughfall precipitation and the !
      ! potential amount of interception.  Later we will adjust both values depending on   !
      ! whether some cohorts are already full or not, or if it is fine to exceed the       !
      ! maximum amount of water that a cohort can hold.                                    !
      !------------------------------------------------------------------------------------!
      if (rk4aux(ibuff)%any_resolvable) then
         taii = 0.d0
         cpatch => csite%patch(ipa)
         do ico = 1,cpatch%ncohorts
            taii = taii + initp%tai(ico)
         end do
         taii = 1.d0/taii

         !---------------------------------------------------------------------------------!
         !    If the canopy does not cover all of the ground, then it should not intercept !
         ! all of the water.                                                               !
         !---------------------------------------------------------------------------------!
         if (rk4site%pcpg > 0.d0) then
            if (leaf_intercept) then
               !----- Scale interception by canopy openess (MCD 01-12-09). ----------------!
               intercepted_max  = rk4site%pcpg  * closedcan_frac
               qintercepted_max = rk4site%qpcpg * closedcan_frac
               dintercepted_max = rk4site%dpcpg * closedcan_frac
            else
               !----- No interception (developer only). -----------------------------------!
               intercepted_max  = 0.d0
               qintercepted_max = 0.d0
               dintercepted_max = 0.d0
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    The first guess for through fall is the rainfall minus the maximum        !
            ! interception.  If some of the water can't be intercepted, we will add to the !
            ! through fall later.                                                          !
            !------------------------------------------------------------------------------!
            throughfall_tot  = rk4site%pcpg  - intercepted_max
            qthroughfall_tot = rk4site%qpcpg - qintercepted_max
            dthroughfall_tot = rk4site%dpcpg - dintercepted_max
            !------------------------------------------------------------------------------!
         else
            !----- No precipitation, nothing to be intercepted... -------------------------!
            intercepted_max  = 0.d0
            qintercepted_max = 0.d0
            dintercepted_max = 0.d0
            throughfall_tot  = 0.d0
            qthroughfall_tot = 0.d0
            dthroughfall_tot = 0.d0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!

      else
         !---------------------------------------------------------------------------------!
         !     If the TAI is very small or total patch vegetation heat capacity is too     !
         ! small, bypass vegetation computations.  Set throughfall precipitation heat and  !
         ! moisture to unintercepted values.                                               !
         !---------------------------------------------------------------------------------!
         intercepted_max  = 0.d0
         qintercepted_max = 0.d0
         dintercepted_max = 0.d0
         throughfall_tot  = rk4site%pcpg
         qthroughfall_tot = rk4site%qpcpg
         dthroughfall_tot = rk4site%dpcpg

         !---------------------------------------------------------------------------------!
         ! Note: If the condition of low TAI for the entire patch was met, then it does    !
         !       not matter what the individual cohorts are normalized by, because they    !
         !       are effectively zero. So make sure the inverse patch TAI is a nominal     !
         !       non-zero/non-infinite number. This will only be used when parsing out     !
         !       intercepted leaf water into shed water; in which case the intercepted     !
         !       water is zero anyway. So this is just to prevent FPEs.                    !
         !---------------------------------------------------------------------------------!
         taii = 0.d0
         do ico = 1,cpatch%ncohorts
            taii = taii + initp%tai(ico)
         end do
         if (taii == 0.d0) then
            taii = 1.d0
         else
            taii = 1.d0/taii
         end if
      end if
      !------------------------------------------------------------------------------------!
      !     Compute sensible heat and moisture fluxes between the ground and the canopy    !
      ! air space.  The ground may be either the top soil layer or the temporary surface   !
      ! water snow surface.  First we check whether the vapour fluxes to ground are        !
      ! positive or negative.  Later, three variables will define the flux between ground  !
      ! and canopy air space.  They are either positive or zero:                           !
      !                                                                                    !
      ! - wflxsc [kg/m2/s] is going to be the evaporation flux from pounding water or      !
      !                    snowpack to canopy.                                             !
      ! - wflxgc [kg/m2/s] is going to be the evaporation flux from top soil layer to      !
      !                    canopy.                                                         !
      ! - dewgnd [kg/m2/s] is going to be the dew/frost flux from canopy to ground (this   !
      !                    encompasses dew formation on both because dew always goes to    !
      !                    the virtual layer).                                             !
      !                                                                                    !
      !     Sensible heat is defined by two variables, hflxgc (soil) and hflxsc (TSW)      !
      ! [J/m2/s], which can be either positive or negative.                                !
      !------------------------------------------------------------------------------------!
      ksn = initp%nlev_sfcwater
      if (ksn > 0) then
         hflxsc       = initp%snowfac * initp%ggnet * initp%can_rhos * initp%can_cp        &
                      * (initp%sfcwater_tempk(ksn) - initp%can_temp)
         sfcwater_ssh = qslif8(initp%can_prss,initp%sfcwater_tempk(ksn))
         is_dew_cp    = sfcwater_ssh <= initp%can_shv
      else
         hflxsc       = 0.d0
         sfcwater_ssh = initp%can_shv
         is_dew_cp    = .false.
      end if
      hflxgc    = (1.d0 - initp%snowfac) * initp%ggnet * initp%can_rhos                    &
                * initp%can_cp * (initp%ground_temp - initp%can_temp)
      is_dew_cs = initp%ground_ssh <= initp%can_shv
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Set all water fluxes to zero, we only add them if they occur and are allowed.  !
      !------------------------------------------------------------------------------------!
      dewgndflx         = 0.d0
      qdewgndflx        = 0.d0
      ddewgndflx        = 0.d0
      wflxsc            = 0.d0
      qwflxsc           = 0.d0
      wflxgc            = 0.d0
      qwflxgc           = 0.d0
      initp%flag_wflxgc = 0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Here we will decide how to compute the evaporation and condensation fluxes     !
      ! between temporary surface water and canopy air space, based on the sign of water   !
      ! flux.                                                                              !
      !------------------------------------------------------------------------------------!
      if (is_dew_cp) then
         !----- Dew should be formed due to contact between top TSW and canopy air space. -!
         dew_now           = initp%snowfac * initp%ggnet * initp%can_rhos                  &
                           * ( initp%can_shv - sfcwater_ssh )
         dewgndflx         = dewgndflx  + dew_now
         qdewgndflx        = qdewgndflx                                                    &
                           + dew_now * tq2enthalpy8(initp%sfcwater_tempk(ksn),1.d0,.true.)
         ddewgndflx        = ddewgndflx                                                    &
                           + dew_now * (          initp%sfcwater_fracliq(ksn)   * wdnsi8   &
                                       + ( 1.d0 - initp%sfcwater_fracliq(ksn) ) * fdnsi8 )
         !----- Set flux flag. ------------------------------------------------------------!
         initp%flag_wflxgc = initp%flag_wflxgc + 1
         !---------------------------------------------------------------------------------!
      elseif (ksn > 0 .and. ( initp%can_rhv < 1.d0 .or. supersat_ok )) then
         !---------------------------------------------------------------------------------!
         !     Evaporation should occur, and canopy air space is not yet super-saturated   !
         ! (or the user is fine with super-saturation).                                    !
         !---------------------------------------------------------------------------------!
         wflxsc            = initp%snowfac * initp%ggnet * initp%can_rhos                  &
                           * ( sfcwater_ssh - initp%can_shv )
         qwflxsc           = wflxsc * tq2enthalpy8(initp%sfcwater_tempk(ksn),1.d0,.true.)
         !----- Set flux flag. ------------------------------------------------------------!
         initp%flag_wflxgc = initp%flag_wflxgc + 2
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Similarly, we will decide how to compute the evaporation and condensation      !
      ! fluxes between top soil layer and canopy air space, based on the sign of water     !
      ! flux.                                                                              !
      !------------------------------------------------------------------------------------!
      if (is_dew_cs) then
         !---- Dew should be formed due to contact between top soil and canopy air space. -!
         dew_now           = ( 1.d0 - initp%snowfac ) * initp%ggnet * initp%can_rhos       &
                           * ( initp%can_shv - initp%ground_ssh )
         dewgndflx         = dewgndflx  + dew_now
         qdewgndflx        = qdewgndflx                                                    &
                           + dew_now * tq2enthalpy8(initp%ground_temp,1.d0,.true.)
         ddewgndflx        = ddewgndflx                                                    &
                           + dew_now * (          initp%ground_fliq   * wdnsi8             &
                                       + ( 1.d0 - initp%ground_fliq ) * fdnsi8 )
         !----- Set flux flag. ------------------------------------------------------------!
         initp%flag_wflxgc = initp%flag_wflxgc + 4
         !---------------------------------------------------------------------------------!

      else if (rk4aux(ibuff)%drysoil(mzg)) then
         !---------------------------------------------------------------------------------!
         !    There should be evaporation, except that there is no water left to be        !
         ! extracted from the ground... Set both evaporation and condensation fluxes to    !
         ! zero.                                                                           !
         !---------------------------------------------------------------------------------!
         initp%flag_wflxgc = initp%flag_wflxgc + 7
         !---------------------------------------------------------------------------------!

      elseif (initp%can_rhv < 1.d0 .or. supersat_ok) then
         !---------------------------------------------------------------------------------!
         !     Evaporation should occur, and canopy air space is not yet super-saturated   !
         ! (or the user is fine with super-saturation).                                    !
         !---------------------------------------------------------------------------------!
         wflxgc            = ( 1.d0 - initp%snowfac ) * initp%ggnet * initp%can_rhos       &
                           * ( initp%ground_shv - initp%can_shv )                          &
                           * ( 1.d0 / (1.d0 + initp%ggnet / initp%ggsoil) )
         qwflxgc           = wflxgc * tq2enthalpy8(initp%ground_temp,1.d0,.true.)
         !----- Set flux flag. ------------------------------------------------------------!
         initp%flag_wflxgc = initp%flag_wflxgc + 10
         !---------------------------------------------------------------------------------!

      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      ! The implicit solver needs to know the mass flux from ground to canopy.             !
      !------------------------------------------------------------------------------------!
      initp%hflxsc  = hflxsc
      initp%wflxsc  = wflxsc
      initp%qwflxsc = qwflxsc
      initp%hflxgc  = hflxgc
      initp%wflxgc  = wflxgc
      initp%qwflxgc = qwflxgc
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop over the cohorts in the patch. Calculate energy fluxes with surrounding   !
      ! canopy air space, integrate cohort energy, calculate precipitation throughfall and !
      ! sum fluxes to the patch level. Initialize variables used to store sums over        !
      ! cohorts.                                                                           !
      !------------------------------------------------------------------------------------!
      hflxlc_tot       = 0.d0
      wflxlc_tot       = 0.d0
      qwflxlc_tot      = 0.d0
      transp_tot       = 0.d0
      qtransp_tot      = 0.d0
      wshed_tot        = 0.d0
      qwshed_tot       = 0.d0
      dwshed_tot       = 0.d0
      intercepted_tot  = 0.d0
      qintercepted_tot = 0.d0
      hflxwc_tot       = 0.d0
      wflxwc_tot       = 0.d0
      qwflxwc_tot      = 0.d0
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Initialise nee_tot with the patch-level variables (heterotrophic respiration    !
      ! commited growth and storage respiration).                                          !
      !------------------------------------------------------------------------------------!
      nee_tot    = initp%fgc_rh              + initp%fsc_rh                                &
                 + initp%stgc_rh             + initp%stsc_rh                               &
                 + initp%msc_rh              + initp%ssc_rh                                &
                 + initp%psc_rh                                                            &
                 + initp%commit_storage_resp + initp%commit_growth_resp
      !------------------------------------------------------------------------------------!

      cohortloop: do ico = 1,cpatch%ncohorts
         ipft = cpatch%pft(ico)

         !---------------------------------------------------------------------------------!
         !    Subtract the metabolic NPP of this cohort from patch-level NEE.              !
         !---------------------------------------------------------------------------------!
         nee_tot = nee_tot - ( initp%gpp(ico) - initp%leaf_resp(ico)                       &
                             - initp%root_resp(ico) - initp%stem_resp(ico))
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         ! LEAF BUDGET - Check whether the leaves of this cohort haven't been flagged as   !
         !               non-resolvable, i.e., it has leaves, belongs to a patch that is   !
         !               not too sparse, an it is not buried in snow.  We should compute   !
         !               energy and water at the cohort level only if the cohort is        !
         !               "safe".  Otherwise, we will set the leaf energy derivatives to    !
         !               zero, and convert any intercepted water into throughfall.  Later  !
         !               on, these "unsafe" cohorts will have their leaf energy set to     !
         !               equilibrium with the canopy air space (temperature).              !
         !---------------------------------------------------------------------------------!
         if (initp%leaf_resolvable(ico)) then

            !------ Define some shortcuts to indices --------------------------------------!
            ipft  = cpatch%pft(ico)
            kroot = cpatch%krdepth(ico)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Calculate interception by leaves by scaling the intercepted water by the !
            ! LAI of each cohort.  If this causes excess of water/ice over the leaf sur-   !
            ! face, no problem, the water will shed at adjust_veg_properties.              !
            !                                                                              !
            ! IMPORTANT, according to RGK, this block must come before the is_hybrid block !
            !            because hybrid predictive capping needs leaf_intercepted.         !
            !------------------------------------------------------------------------------!
            wshed             = 0.d0
            qwshed            = 0.d0
            dwshed            = 0.d0
            leaf_intercepted  = intercepted_max  * initp%lai(ico) * taii
            leaf_qintercepted = qintercepted_max * initp%lai(ico) * taii
            throughfall       = 0.d0
            qthroughfall      = 0.d0
            dthroughfall      = 0.d0
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Define the minimum leaf water to be considered, and the maximum amount   !
            ! possible.                                                                    !
            !------------------------------------------------------------------------------!
            min_leaf_water = rk4leaf_drywhc * initp%lai(ico)
            max_leaf_water = rk4leaf_maxwhc * initp%lai(ico)
            !------------------------------------------------------------------------------!



            !------ Calculate fraction of leaves covered with water. ----------------------!
            if (initp%leaf_water(ico) > min_leaf_water) then
               sigmaw = min(1.d0, (initp%leaf_water(ico)/max_leaf_water)**twothirds8)
            else
               sigmaw = 0.d0
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Here we must compute two different areas.  For transpiration, we want the !
            ! leaf area index only, because we assume transpiration to happen only through !
            ! leaves.  Evaporation of water/ice settled over the vegetation surface, or    !
            ! dew/frost formation must account the branches and stems as well.             !
            !------------------------------------------------------------------------------!
            !----- Evaporation/condensation "flux" ----------------------------------------!
            wflxlc_try = effarea_evap * initp%lai(ico) * initp%leaf_gbw(ico)               &
                       * (initp%lint_shv(ico) - initp%can_shv)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Computing the evapotranspiration or dew/frost deposition.                 !
            !------------------------------------------------------------------------------!
            if (wflxlc_try >= 0.d0) then
               !---------------------------------------------------------------------------!
               !    Probably evapotranspiration, as long as the canopy air is not          !
               ! saturated or the user doesn't mind that super-saturation occur.           !
               !---------------------------------------------------------------------------!
               if (supersat_ok .or. initp%can_rhv < 1.d0) then
                  !------------------------------------------------------------------------!
                  !     Evaporation, energy is scaled by liquid/ice partition (no phase    !
                  ! bias).  We scale by the relative area of leaves that is actually       !
                  ! covered with water.                                                    !
                  !------------------------------------------------------------------------!
                  wflxlc  = wflxlc_try * sigmaw
                  qwflxlc = wflxlc * tq2enthalpy8(initp%leaf_temp(ico),1.d0,.true.)



                  !------------------------------------------------------------------------!
                  !       This is called by the hybrid solver only.                        !
                  !------------------------------------------------------------------------!
                  if (is_hybrid) then

                     max_dwdt = initp%leaf_water(ico)/dt

                     !---------------------------------------------------------------------!
                     !     If we ever have shedding, force wshed to cap out at that        !
                     ! maximum leaf water. Assume this process happens before evaporation. !
                     !---------------------------------------------------------------------!
                     !! TURNING OFF SHEDDING FOR NOW

                     !! wshed  = max(0.d0,( (initp%leaf_water(ico) + leaf_intercepted*dt)  &
                     !!                   - max_leaf_water) / dt)
                     !! qwshed = wshed                                                     &
                     !!        * tl2uint8(initp%leaf_temp(ico),initp%leaf_fliq(ico))   
                     !! dwshed = wshed * ( initp%leaf_fliq(ico) * wdnsi8                   &
                     !!                    + (1.d0-initp%leaf_fliq(ico)) * fdnsi8)
                     !---------------------------------------------------------------------!


                     !----- Then constrain the amount that can be evaporated. -------------!
                     wflxlc = min(wflxlc,max_dwdt+leaf_intercepted-wshed)
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Transpiration, consider the one-sided leaf area rather than LAI,   !
                  ! when the PFT is amphistomatous, in which case it is double-sided.      !
                  ! Compute the water demand from both open closed and open stomata, but   !
                  ! first make sure that there is some water available for transpiration.  !
                  !     When running plant hydrodynamics, the soil water constraint should !
                  ! be Turned off. Instead, we use leaf water potential as water avail-    !
                  ! ability for transpiration.                                             !
                  !------------------------------------------------------------------------!
                  select case (plant_hydro_scheme)
                  case (0)
                     !----- Static plant hydraulics.  Check soil water availability. ------!
                     transp_fine = rk4aux(ibuff)%avail_h2o_int(kroot) >  0.d0
                     !---------------------------------------------------------------------!
                  case default
                     !----- Dynamic plant hydraulics. Check leaf water potential. ---------!
                     op_rk4eps      = 1.d0 + rk4eps
                     lwater_im2_lwr = op_rk4eps * rk4aux(ibuff)%rk4min_leaf_water_im2(ico)
                     wwater_im2_lwr = op_rk4eps * rk4aux(ibuff)%rk4min_wood_water_im2(ico)
                     lwater_fine    = initp%leaf_water_im2(ico) >= lwater_im2_lwr
                     wwater_fine    = initp%wood_water_im2(ico) >= wwater_im2_lwr
                     transp_fine    = lwater_fine .and.                                    &
                                      ( wwater_fine .or. ( .not. cpatch%is_small(ico) ) )
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Update leaf-level and patch-level transpiration rates.             !
                  !------------------------------------------------------------------------!
                  if (transp_fine) then
                     !---------------------------------------------------------------------!
                     !    Transpiration can occur, use stomatal conductance and leaf-level !
                     ! specific humidity deficit to obtain transpiration.                  !
                     !---------------------------------------------------------------------!
                     gleaf_open   = effarea_transp(ipft)                                   &
                                  * initp%leaf_gbw(ico) * initp%gsw_open(ico)              &
                                  / (initp%leaf_gbw(ico) + initp%gsw_open(ico) )
                     gleaf_closed = effarea_transp(ipft)                                   &
                                  * initp%leaf_gbw(ico) * initp%gsw_closed(ico)            &
                                  / ( initp%leaf_gbw(ico) + initp%gsw_closed(ico) )
                     shv_gradient = initp%lint_shv(ico) - initp%can_shv

                     dinitp%psi_open  (ico) = gleaf_open   * shv_gradient
                     dinitp%psi_closed(ico) = gleaf_closed * shv_gradient

                     transp = initp%lai(ico) * ( initp%fs_open(ico) * dinitp%psi_open(ico) &
                                               + (1.0d0 - initp%fs_open(ico))              &
                                               * dinitp%psi_closed(ico) )
                     !---------------------------------------------------------------------!
                  else
                     !---------------------------------------------------------------------!
                     !     Too dry for transpiration, shut down transpiration completely.  !
                     !---------------------------------------------------------------------!
                     dinitp%psi_open  (ico) = 0.d0
                     dinitp%psi_closed(ico) = 0.d0
                     transp                 = 0.d0
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Only liquid water is transpired, thus this is always the            !
                  ! condensation latent heat.                                              !
                  !------------------------------------------------------------------------!
                  qtransp = transp * tq2enthalpy8(initp%leaf_temp(ico),1.d0,.true.)
                  !------------------------------------------------------------------------!

               else
                  !----- Canopy is already saturated, no evapotranspiration is allowed. ---!
                  wflxlc                 = 0.d0
                  qwflxlc                = 0.d0
                  transp                 = 0.d0
                  qtransp                = 0.d0
                  dinitp%psi_open(ico)   = 0.d0
                  dinitp%psi_closed(ico) = 0.d0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!

            else
               !---------------------------------------------------------------------------!
               !     Dew/frost formation.                                                  !
               !---------------------------------------------------------------------------!
               wflxlc                 = wflxlc_try
               qwflxlc                = wflxlc                                             &
                                      * tq2enthalpy8(initp%leaf_temp(ico),1.d0,.true.)
               transp                 = 0.0d0
               qtransp                = 0.0d0
               dinitp%psi_open  (ico) = 0.0d0
               dinitp%psi_closed(ico) = 0.0d0
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!






            !----- We need to extract water from the soil equal to the transpiration. -----!
            rk4aux(ibuff)%extracted_water(ico,kroot) = transp                              &
                                                 + rk4aux(ibuff)%extracted_water(ico,kroot)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !   Calculate leaf-to-canopy sensible heat flux.  Always consider both sides   !
            ! of leaves.                                                                   !
            !------------------------------------------------------------------------------!
            flux_area = effarea_heat * initp%lai(ico)
            hflxlc    = flux_area    * initp%leaf_gbh(ico)                                 &
                      * (initp%leaf_temp(ico) - initp%can_temp)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Find the water energy balance for this cohort.                           !
            !------------------------------------------------------------------------------!
            dinitp%leaf_water(ico)  = leaf_intercepted     & ! Intercepted water
                                    - wflxlc               & ! Evaporation
                                    - wshed                ! ! Water shedding
            dinitp%leaf_energy(ico) = initp%rshort_l(ico)  & ! Absorbed SW radiation
                                    + initp%rlong_l(ico)   & ! Net thermal radiation
                                    - hflxlc               & ! Sensible heat flux
                                    - qwflxlc              & ! Evaporation
                                    - qwshed               & ! Water shedding
                                    - qtransp              & ! Transpiration
                                    + leaf_qintercepted    ! ! Intercepted water energy
            !------------------------------------------------------------------------------!


            initp%wflxlc(ico)     = wflxlc
            initp%wflxtr(ico)     = transp
            initp%hflx_lrsti(ico) = initp%rshort_l(ico) + initp%rlong_l(ico)               &
                                  - qwshed + leaf_qintercepted


            !------------------------------------------------------------------------------!
            !    If we are saving fast diagnostics, then we save the fluxes for this       !
            ! cohort.                                                                      !
            !------------------------------------------------------------------------------!
            if (fast_diagnostics .or. print_detailed) then
               dinitp%avg_sensible_lc      (ico)  = hflxlc
               dinitp%avg_vapor_lc         (ico)  = wflxlc
               dinitp%avg_transp           (ico)  = transp
               dinitp%avg_intercepted_al   (ico)  = leaf_intercepted
               dinitp%avg_wshed_lg         (ico)  = wshed
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    If the detailed output is tracked, then we save the fluxes for this       !
            ! cohort.                                                                      !
            !------------------------------------------------------------------------------!
            if (print_detailed) then
               dinitp%cfx_hflxlc      (ico)  = hflxlc
               dinitp%cfx_qwflxlc     (ico)  = qwflxlc
               dinitp%cfx_qwshed      (ico)  = qwshed
               dinitp%cfx_qtransp     (ico)  = qtransp
               dinitp%cfx_qintercepted(ico)  = leaf_qintercepted
            end if
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !    Add the contribution of this cohort to total heat and evapotranspiration. !
            !------------------------------------------------------------------------------!
            wflxlc_tot   = wflxlc_tot   + wflxlc
            qwflxlc_tot  = qwflxlc_tot  + qwflxlc
            hflxlc_tot   = hflxlc_tot   + hflxlc
            transp_tot   = transp_tot   + transp
            qtransp_tot  = qtransp_tot  + qtransp

            !------------------------------------------------------------------------------!
            !     Here we update the liquid/frozen water fluxes and their associated vari- !
            ! ables:                                                                       !
            !                                                                              !
            ! - wshed      : Water falling from vegetated canopy to soil surface;          !
            ! - intercepted: Precipitation that is intercepted by the vegetation;          !
            ! - throughfall: Precipitation that is never intercepted by the vegetation.    !
            !------------------------------------------------------------------------------!
            wshed_tot        = wshed_tot        + wshed
            qwshed_tot       = qwshed_tot       + qwshed
            dwshed_tot       = dwshed_tot       + dwshed
            intercepted_tot  = intercepted_tot  + leaf_intercepted
            qintercepted_tot = qintercepted_tot + leaf_qintercepted
            throughfall_tot  = throughfall_tot  + throughfall
            qthroughfall_tot = qthroughfall_tot + qthroughfall
            dthroughfall_tot = dthroughfall_tot + dthroughfall
         else
            !----- Leaf is not resolved. Set transpiration to zero. -----------------------!
            transp = 0.d0
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     If there is not enough leaf biomass to safely solve the leaf energy and  !
            ! water balances, set leaf fluxes and interception to zero.                    !
            !------------------------------------------------------------------------------!
            dinitp%leaf_energy(ico) = 0.d0
            dinitp%leaf_water (ico) = 0.d0
            dinitp%psi_open   (ico) = 0.d0
            dinitp%psi_closed (ico) = 0.d0
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    If we are saving fast diagnostics, then we save the fluxes for this       !
            ! cohort.                                                                      !
            !------------------------------------------------------------------------------!
            if (fast_diagnostics .or. print_detailed) then
               dinitp%avg_sensible_lc      (ico)  = 0.d0
               dinitp%avg_vapor_lc         (ico)  = 0.d0
               dinitp%avg_transp           (ico)  = 0.d0
               dinitp%avg_intercepted_al   (ico)  = 0.d0
               dinitp%avg_wshed_lg         (ico)  = 0.d0
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    If the detailed output is tracked, then we save the fluxes for this       !
            ! cohort.                                                                      !
            !------------------------------------------------------------------------------!
            if (print_detailed) then
               dinitp%cfx_hflxlc      (ico)  = 0.d0
               dinitp%cfx_qwflxlc     (ico)  = 0.d0
               dinitp%cfx_qwshed      (ico)  = 0.d0
               dinitp%cfx_qtransp     (ico)  = 0.d0
               dinitp%cfx_qintercepted(ico)  = 0.d0
            end if
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     Allow the complete bypass of precipitation if there are very few leaves. !
            ! Add this tiny amount to the throughfall.                                     !
            !------------------------------------------------------------------------------!
            throughfall_tot   = throughfall_tot  + intercepted_max  * initp%lai(ico) * taii
            qthroughfall_tot  = qthroughfall_tot + qintercepted_max * initp%lai(ico) * taii
            dthroughfall_tot  = dthroughfall_tot + dintercepted_max * initp%lai(ico) * taii
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!






         !---------------------------------------------------------------------------------!
         ! WOOD BUDGET - Check whether the wood of this cohort hasn't been flagged as non- !
         !               resolvable, i.e., the user wants wood thermodynamics, it belongs  !
         !               to a patch that is not too sparse, an it is not buried in snow.   !
         !               We should compute energy and water at the cohort level only if    !
         !               the cohort is "safe".  Otherwise, we will set the wood energy and !
         !               mass derivatives to zero, and convert any intercepted water into  !
         !               throughfall.  Later on, these "unsafe" cohorts will have their    !
         !               wood energy set to equilibrium with the canopy air space          !
         !               (temperature).                                                    !
         !---------------------------------------------------------------------------------!
         if (initp%wood_resolvable(ico)) then

            !------ Define some aliases to indices ----------------------------------------!
            ipft  = cpatch%pft(ico)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Calculate interception by wood by scaling the intercepted water by the   !
            ! WAI of each cohort.  If this causes excess of water/ice over the wood sur-   !
            ! face, no problem, the water will shed at adjust_veg_properties.              !
            !                                                                              !
            ! IMPORTANT, according to RGK, this block must come before the is_hybrid block !
            !            because hybrid predictive capping needs wood_intercepted.         !
            !------------------------------------------------------------------------------!
            wshed             = 0.d0
            qwshed            = 0.d0
            dwshed            = 0.d0
            wood_intercepted  = intercepted_max  * initp%wai(ico) * taii
            wood_qintercepted = qintercepted_max * initp%wai(ico) * taii
            throughfall       = 0.d0
            qthroughfall      = 0.d0
            dthroughfall      = 0.d0
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Define the minimum wood water to be considered, and the maximum amount   !
            ! possible.                                                                    !
            !------------------------------------------------------------------------------!
            min_wood_water = rk4leaf_drywhc * initp%wai(ico)
            max_wood_water = rk4leaf_maxwhc * initp%wai(ico)
            !------------------------------------------------------------------------------!

            !------ Calculate fraction of wood covered with water. ------------------------!
            if (initp%wood_water(ico) > min_wood_water) then
               sigmaw = min(1.d0, (initp%wood_water(ico)/max_wood_water)**twothirds8)
            else
               sigmaw = 0.d0
            end if

            !------------------------------------------------------------------------------!
            !    Here we must compute two different areas.  For transpiration, we want the !
            ! leaf area index only, because we assume transpiration to happen only through !
            ! leaves.  Evaporation of water/ice settled over the vegetation surface, or    !
            ! dew/frost formation must account the branches and stems as well.             !
            !------------------------------------------------------------------------------!
            !----- Find the wood specific humidity. ---------------------------------------!
            wood_shv   = qslif8(initp%can_prss,initp%wood_temp(ico))
            !----- Evaporation/condensation "flux" ----------------------------------------!
            wflxwc_try = initp%wai(ico) * initp%wood_gbw(ico) * (wood_shv - initp%can_shv)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Compute the evapotranspiration or dew/frost deposition.                   !
            !------------------------------------------------------------------------------!
            if (wflxwc_try >= 0.d0) then
               !---------------------------------------------------------------------------!
               !    Probably evapotranspiration, as long as the canopy air is not          !
               ! saturated or the user doesn't mind that super-saturation occur.           !
               !---------------------------------------------------------------------------!
               if (supersat_ok .or. initp%can_rhv < 1.d0) then
                  !------------------------------------------------------------------------!
                  !     Evaporation, energy is scaled by liquid/ice partition (no phase    !
                  ! bias).  We scale by the relative area of wood that is actually covered !
                  ! with water.                                                            !
                  !------------------------------------------------------------------------!
                  wflxwc  = wflxwc_try * sigmaw
                  qwflxwc = wflxwc * tq2enthalpy8(initp%wood_temp(ico),1.d0,.true.)
                  !------------------------------------------------------------------------!

               else
                  !----- Canopy is already saturated, no evapotranspiration is allowed. ---!
                  wflxwc                 = 0.d0
                  qwflxwc                = 0.d0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!

            else
               !---------------------------------------------------------------------------!
               !     Dew/frost formation. The deposition will conserve the liquid/ice      !
               ! partition (or use the default if there is no water).                      !
               !---------------------------------------------------------------------------!
               wflxwc                 = wflxwc_try
               qwflxwc                = wflxwc                                             &
                                      * tq2enthalpy8(initp%wood_temp(ico),1.d0,.true.)
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !       This is called by the hybrid solver only.                              !
            !------------------------------------------------------------------------------!
            if (is_hybrid) then

               max_dwdt = initp%wood_water(ico)/dt

               !---------------------------------------------------------------------------!
               !     If we ever have shedding, force wshed to cap out at that maximum leaf !
               ! water.  Assume this process happens before evaporation.                   !
               !---------------------------------------------------------------------------!
               !! TURNING OFF SHEDDING FOR NOW

               !! wshed  = max(0.d0,( (initp%wood_water(ico) + wood_intercepted*dt)        &
               !!                   - max_wood_water) / dt)
               !! qwshed = wshed                                                           &
               !!        * tl2uint8(initp%wood_temp(ico),initp%wood_fliq(ico))   
               !! dwshed = wshed * ( initp%wood_fliq(ico) * wdnsi8                         &
               !!                    + (1.d0-initp%wood_fliq(ico)) * fdnsi8)
               !---------------------------------------------------------------------------!


               !----- Then constrain the amount that can be evaporated. -------------------!
               wflxwc = min(wflxwc,max_dwdt+wood_intercepted-wshed)
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!





            !------------------------------------------------------------------------------!
            !   Calculate wood-to-canopy sensible heat flux.  Consider the full            !
            ! circumference of brances.                                                    !
            !------------------------------------------------------------------------------!
            flux_area = pi18 * initp%wai(ico)
            hflxwc    = flux_area * initp%wood_gbh(ico)                                    &
                      * (initp%wood_temp(ico) - initp%can_temp)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Find the water energy balance for this cohort.                           !
            !------------------------------------------------------------------------------!
            dinitp%wood_water(ico)  = wood_intercepted     & ! Intercepted water
                                    - wflxwc               & ! Evaporation
                                    - wshed                ! ! Water shedding
            dinitp%wood_energy(ico) = initp%rshort_w(ico)  & ! Absorbed SW radiation
                                    + initp%rlong_w(ico)   & ! Net thermal radiation
                                    - hflxwc               & ! Sensible heat flux
                                    - qwflxwc              & ! Evaporation
                                    - qwshed               & ! Water shedding
                                    + wood_qintercepted    ! ! Intercepted water energy
            !------------------------------------------------------------------------------!

            initp%wflxwc(ico) = wflxwc
            initp%hflx_wrsti(ico) = initp%rshort_w(ico) + initp%rlong_w(ico)               &
                                  - qwshed + wood_qintercepted



            !------------------------------------------------------------------------------!
            !    If we are saving fast diagnostics, then we save the fluxes for this       !
            ! cohort.                                                                      !
            !------------------------------------------------------------------------------!
            if (fast_diagnostics .or. print_detailed) then
               dinitp%avg_sensible_wc      (ico)  = hflxwc
               dinitp%avg_vapor_wc         (ico)  = wflxwc
               dinitp%avg_intercepted_aw   (ico)  = wood_intercepted
               dinitp%avg_wshed_wg         (ico)  = wshed
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    If the detailed output is tracked, then we save the fluxes for this       !
            ! cohort.                                                                      !
            !------------------------------------------------------------------------------!
            if (print_detailed) then
               dinitp%cfx_hflxwc      (ico)  = hflxwc
               dinitp%cfx_qwflxwc     (ico)  = qwflxwc
               dinitp%cfx_qwshed      (ico)  = dinitp%cfx_qwshed(ico) + qwshed
               dinitp%cfx_qintercepted(ico)  = dinitp%cfx_qintercepted(ico)                &
                                             + wood_qintercepted
            end if
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !    Add the contribution of this cohort to total heat and evapotranspiration. !
            !------------------------------------------------------------------------------!
            wflxwc_tot   = wflxwc_tot   + wflxwc
            qwflxwc_tot  = qwflxwc_tot  + qwflxwc
            hflxwc_tot   = hflxwc_tot   + hflxwc
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     Here we update the liquid/frozen water fluxes and their associated vari- !
            ! ables:                                                                       !
            !                                                                              !
            ! - wshed      : Water falling from vegetated canopy to soil surface;          !
            ! - intercepted: Precipitation that is intercepted by the vegetation;          !
            ! - throughfall: Precipitation that is never intercepted by the vegetation.    !
            !------------------------------------------------------------------------------!
            wshed_tot        = wshed_tot        + wshed
            qwshed_tot       = qwshed_tot       + qwshed
            dwshed_tot       = dwshed_tot       + dwshed
            intercepted_tot  = intercepted_tot  + wood_intercepted
            qintercepted_tot = qintercepted_tot + wood_qintercepted
            throughfall_tot  = throughfall_tot  + throughfall
            qthroughfall_tot = qthroughfall_tot + qthroughfall
            dthroughfall_tot = dthroughfall_tot + dthroughfall
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !     If there is not enough leaf biomass to safely solve the leaf energy and  !
            ! water balances, set leaf fluxes and interception to zero.                    !
            !------------------------------------------------------------------------------!
            dinitp%wood_energy(ico) = 0.d0
            dinitp%wood_water(ico)  = 0.d0
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    If we are saving fast diagnostics, then we save the fluxes for this       !
            ! cohort.                                                                      !
            !------------------------------------------------------------------------------!
            if (fast_diagnostics .or. print_detailed) then
               dinitp%avg_sensible_wc      (ico)  = 0.d0
               dinitp%avg_vapor_wc         (ico)  = 0.d0
               dinitp%avg_intercepted_aw   (ico)  = 0.d0
               dinitp%avg_wshed_wg         (ico)  = 0.d0
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    If the detailed output is tracked, then we save the fluxes for this       !
            ! cohort.                                                                      !
            !------------------------------------------------------------------------------!
            if (print_detailed) then
               dinitp%cfx_hflxwc      (ico)  = 0.d0
               dinitp%cfx_qwflxwc     (ico)  = 0.d0
            end if
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     Allow the complete bypass of precipitation if there are very few leaves. !
            ! Add this tiny amount to the throughfall.                                     !
            !------------------------------------------------------------------------------!
            throughfall_tot   = throughfall_tot  + intercepted_max  * initp%wai(ico) * taii
            qthroughfall_tot  = qthroughfall_tot + qintercepted_max * initp%wai(ico) * taii
            dthroughfall_tot  = dthroughfall_tot + dintercepted_max * initp%wai(ico) * taii
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Solve plant hydraulic derivatives. We always need to update the leaf/wood   !
         ! internal water when plant_hydro_scheme is non-zero.                             !
         !---------------------------------------------------------------------------------!
         select case (plant_hydro_scheme)
         case (0)
            !------------------------------------------------------------------------------!
            !    Plant hydraulics is disabled, water_im2 does not change, set the          !
            ! derivatives to zero.                                                         !
            !------------------------------------------------------------------------------!
            dinitp%leaf_water_im2(ico) = 0.d0
            dinitp%wood_water_im2(ico) = 0.d0
            dinitp%veg_water_im2 (ico) = 0.d0
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    If we are saving fast diagnostics, then we save the fluxes for this       !
            ! cohort.                                                                      !
            !------------------------------------------------------------------------------!
            if (fast_diagnostics .or. print_detailed) then
               dinitp%avg_wflux_wl (ico)  = 0.d0
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    If the detailed output is tracked, then we save the fluxes for this       !
            ! cohort.                                                                      !
            !------------------------------------------------------------------------------!
            if (print_detailed) then
               dinitp%cfx_qwflux_wl(ico)  = 0.d0
            end if
            !------------------------------------------------------------------------------!
         case default
            !------------------------------------------------------------------------------!
            !    Plant hydraulics is enabled, update internal water and internal energy.   !
            !------------------------------------------------------------------------------!

            !----- Update water (soil water uptake is accounted for in leaftw_derivs. -----!
            wflux_wl                   = dble(cpatch%wflux_wl(ico))                        &
                                       * dble(cpatch%nplant  (ico))
            dinitp%leaf_water_im2(ico) =   wflux_wl - transp 
            dinitp%wood_water_im2(ico) = - wflux_wl
            dinitp%veg_water_im2 (ico) = - transp
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Calculate the internal energy flow associated with wood-leaf flux        !
            ! consistent with the sign of the flux.                                        !
            !------------------------------------------------------------------------------!
            if (wflux_wl >= 0.d0) then
               qwflux_wl = wflux_wl * tl2uint8(initp%wood_temp(ico),1.d0)
            else
               qwflux_wl = wflux_wl * tl2uint8(initp%leaf_temp(ico),1.d0)
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Update energy.  Variable qwflux_wl is the internal energy flux associated !
            ! with sapflow.                                                                !
            !------------------------------------------------------------------------------!
            dinitp%leaf_energy(ico) = dinitp%leaf_energy(ico) + qwflux_wl
            dinitp%wood_energy(ico) = dinitp%wood_energy(ico) - qwflux_wl
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    If we are saving fast diagnostics, then we save the fluxes for this       !
            ! cohort.                                                                      !
            !------------------------------------------------------------------------------!
            if (fast_diagnostics .or. print_detailed) then
               dinitp%avg_wflux_wl (ico)  = wflux_wl
               if (print_detailed) then
                  dinitp%cfx_qwflux_wl(ico)  = qwflux_wl
               end if
            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


         !------ Find the combined leaf + wood derivative. --------------------------------!
         dinitp%veg_energy(ico) = dinitp%leaf_energy(ico) + dinitp%wood_energy(ico)
         dinitp%veg_water (ico) = dinitp%leaf_water (ico) + dinitp%wood_water (ico)
         !---------------------------------------------------------------------------------!
      end do cohortloop
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Update the intensive enthalpy, specific humidity, and CO2 mixing ratio.        !
      !------------------------------------------------------------------------------------!
      dinitp%can_enthalpy = ( hflxsc      + hflxgc      + hflxlc_tot                       &
                            + hflxwc_tot  + qwflxsc     + qwflxgc                          &
                            - qdewgndflx  + qwflxlc_tot + qwflxwc_tot                      &
                            + qtransp_tot + eflxac                    ) * initp%hcapcani
      dinitp%can_shv      = ( wflxsc      + wflxgc      - dewgndflx                        &
                            + wflxlc_tot  + wflxwc_tot  + transp_tot                       &
                            + wflxac                                  ) * initp%wcapcani
      dinitp%can_co2      = ( nee_tot     + cflxac                    ) * initp%ccapcani
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      if (.false.) then
      !if (is_hybrid) then

         a = ( nee_tot                                                                     &
             + initp%can_dmol*initp%ggbare*rk4site%atm_co2 )* initp%ccapcani

         b  = (initp%can_dmol*initp%ggbare) * initp%ccapcani
         c0 = initp%can_co2

         ! Calculate the effective derivative
         !      dinitp%can_co2 = ((a/b) + (c0-(a/b))*exp(-b*dt) - c0)/dt

         ! Calculate the effective cflxac term

         cflxac = (initp%can_dmol*initp%ggbare*initp%ccapcani)/dt                          &
                * (rk4site%atm_co2*dt - ((a/b)*dt - c0*exp(-b*dt)/b +                      &
                   c0/b + (a/b)*exp(-b*dt)/b - (a/b)/b  ))

         dinitp%can_co2 = ( nee_tot + cflxac) * initp%ccapcani

      end if
      !------------------------------------------------------------------------------------!

      initp%wflxac  = wflxac



      !------------------------------------------------------------------------------------!
      !     Water deficit.                                                                 !
      !------------------------------------------------------------------------------------!
      dinitp%water_deficit = - (wflxac + rk4site%pcpg)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Integrate diagnostic variables - These are not activated unless fast file-type !
      ! outputs are selected. This will speed up the integrator.                           !
      !------------------------------------------------------------------------------------!
      if (fast_diagnostics .or. print_detailed) then


         dinitp%avg_carbon_ac    = cflxac                       ! Carbon flx,  Atmo->Canopy
         dinitp%avg_carbon_st    = nee_tot + cflxac             ! Carbon storage flux
         dinitp%avg_sensible_ac  = hflxac                       ! Sens. heat,  Atmo->Canopy
         dinitp%avg_vapor_ac     = wflxac                       ! Lat.  heat,  Atmo->Canopy

         dinitp%avg_sensible_gc  = hflxsc + hflxgc              ! Sens. heat,  Grnd->Canopy
         dinitp%avg_vapor_gc     = wflxsc + wflxgc - dewgndflx  ! Lat.  heat,  Canopy->Grnd

         dinitp%avg_throughfall  = throughfall_tot              ! Throughfall,   Atmo->Grnd
         dinitp%avg_qthroughfall = qthroughfall_tot             ! Throughfall,   Atmo->Grnd
         !---------------------------------------------------------------------------------!

         !------ These are used to compute the averages of the star terms. ----------------!
         dinitp%avg_ustar = initp%ustar
         dinitp%avg_tstar = initp%tstar
         dinitp%avg_qstar = initp%qstar
         dinitp%avg_cstar = initp%cstar
         !---------------------------------------------------------------------------------!

      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update the budget variables.                                                   !
      !------------------------------------------------------------------------------------!
      if (checkbudget) then
         dinitp%co2budget_loss2atm = - cflxac
         dinitp%ebudget_loss2atm   = - eflxac
         dinitp%wbudget_loss2atm   = - wflxac
         dinitp%co2budget_storage  = dinitp%co2budget_storage + nee_tot + cflxac
         dinitp%ebudget_netrad     = compute_netrad8(csite,ipa)
         dinitp%ebudget_storage    = dinitp%ebudget_storage   + dinitp%ebudget_netrad      &
                                   + rk4site%qpcpg - dinitp%ebudget_loss2atm
         dinitp%wbudget_storage    = dinitp%wbudget_storage + rk4site%pcpg                 &
                                   - dinitp%wbudget_loss2atm
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     These variables below are virtual copies of the variables above, but are here  !
      ! for for use in the coupled model. They form the set of canopy-atmosphere fluxes    !
      ! that are used for turbulent closure. These variables are also zeroed and           !
      ! normalized every dtlsm timestep, the others are likely averaged over the analysis  !
      ! period.                                                                            !
      !------------------------------------------------------------------------------------!
      dinitp%upwp = -(initp%ustar*initp%ustar)
      dinitp%qpwp = -(initp%ustar*initp%qstar)
      dinitp%cpwp = -(initp%ustar*initp%cstar)
      dinitp%tpwp = -(initp%ustar*initp%tstar)
      dinitp%wpwp = vertical_vel_flux8(initp%zeta,initp%ustar)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If the single pond layer is too thin, we force equilibrium with top soil       !
      ! layer.  This is done in two steps: first, we transfer all the energy to the top    !
      ! soil layer.  Then, in adjust_sfcw_properties, we make them in thermal equilibrium. !
      ! We use the flag_sfcwater to transfer the energy rather than the actual mass be-    !
      ! cause if a layer is considered to thin at the beginning of the time step, we want  !
      ! it to keep the same status throughout the entire step.                             !
      !------------------------------------------------------------------------------------!
      select case (initp%flag_sfcwater)
      case (1)
         dinitp%soil_energy(mzg)   = dinitp%soil_energy(mzg)                               &
                                   + dinitp%sfcwater_energy(1) * dslzi8(mzg)
         dinitp%sfcwater_energy(1) = 0.d0
      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine canopy_derivs_two
   !=======================================================================================!
   !=======================================================================================!
end module rk4_derivs
!==========================================================================================!
!==========================================================================================!

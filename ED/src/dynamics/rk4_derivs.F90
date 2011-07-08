!==========================================================================================!
!==========================================================================================!
! Subroutine leaf_derivs                                                                   !
!                                                                                          !
!     This subroutine finds the fast-scale derivatives at canopy, soil, and leaf surface.  !
! This subroutine is based on LEAF-3, except that here only the derivative is computed,    !
! whereas in LEAF-3 the actual step is done at once. This derivative will be used for the  !
! Runge-Kutta integration step.                                                            !
!------------------------------------------------------------------------------------------!
subroutine leaf_derivs(initp,dinitp,csite,ipa)
  
   use rk4_coms               , only : rk4site            & ! intent(in)
                                     , rk4patchtype       ! ! structure
   use ed_state_vars          , only : sitetype           & ! structure
                                     , polygontype        ! ! structure
   use consts_coms            , only : cp8                & ! intent(in)
                                     , cpi8               ! ! intent(in)
   use grid_coms              , only : nzg                & ! intent(in)
                                     , nzs                ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: initp     ! Structure with RK4 intermediate state
   type(rk4patchtype) , target     :: dinitp    ! Structure with RK4 derivatives
   type(sitetype)     , target     :: csite     ! This site (with previous values);
   integer            , intent(in) :: ipa       ! Patch ID
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   ! Depending on the type of compilation, interfaces must be explicitly declared.         !
   !---------------------------------------------------------------------------------------!
#if USE_INTERF
   interface
      !------------------------------------------------------------------------------------!
      !    Subroutine that computes the canopy and leaf fluxes.                            ! 
      !------------------------------------------------------------------------------------!
      subroutine leaftw_derivs(mzg,mzs,initp,dinitp,csite,ipa)
         use rk4_coms      , only : rk4patchtype ! ! structure
         use ed_state_vars , only : sitetype,polygontype     ! ! structure
         implicit none
         !----- Arguments -----------------------------------------------------------------!
         type(rk4patchtype)  , target     :: initp  ! RK4 structure, intermediate step
         type(rk4patchtype)  , target     :: dinitp ! RK4 structure, derivatives
         type(sitetype)      , target     :: csite  ! Current site (before integration)
         integer             , intent(in) :: ipa    ! Current patch ID
         integer             , intent(in) :: mzg    ! Number of ground layers
         integer             , intent(in) :: mzs    ! Number of snow/ponding layers
      end subroutine leaftw_derivs
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!

   !----- Ensure that theta_Eiv and water storage derivatives are both zero. --------------!
   dinitp%ebudget_storage   = 0.d0
   dinitp%wbudget_storage   = 0.d0
   dinitp%co2budget_storage = 0.d0

   !----- Finding the derivatives. --------------------------------------------------------!
   call leaftw_derivs(nzg,nzs,initp,dinitp,csite,ipa)

   return
end subroutine leaf_derivs
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine leaftw_derivs(mzg,mzs,initp,dinitp,csite,ipa)
   use ed_max_dims          , only : nzgmax                & ! intent(in)
                                   , nzsmax                ! ! intent(in)
   use consts_coms          , only : alvl8                 & ! intent(in)
                                   , cliqvlme8             & ! intent(in)
                                   , tsupercool8           & ! intent(in)
                                   , wdns8                 & ! intent(in)
                                   , wdnsi8                & ! intent(in)
                                   , lnexp_min8            ! ! intent(in)
   use soil_coms            , only : soil8                 & ! intent(in)
                                   , slz8                  & ! intent(in)
                                   , dslz8                 & ! intent(in)
                                   , dslzi8                & ! intent(in)
                                   , infiltration_method   & ! intent(in)
                                   , dslzti8               & ! intent(in)
                                   , slcons18              & ! intent(in)
                                   , slzt8                 & ! intent(in)
                                   , ss                    & ! intent(in)
                                   , isoilbc               & ! intent(in)
                                   , freezecoef8           ! ! intent(in)
   use ed_misc_coms         , only : dtlsm                 & ! intent(in)
                                   , current_time          & ! intent(in)
                                   , fast_diagnostics      ! ! intent(in)
   use rk4_coms             , only : rk4eps                & ! intent(in)
                                   , rk4tiny_sfcw_mass     & ! intent(in)
                                   , checkbudget           & ! intent(in)
                                   , any_resolvable        & ! intent(in)
                                   , rk4site               & ! intent(in)
                                   , rk4patchtype          & ! structure
                                   , print_detailed        & ! intent(in)
                                   , rk4aux                & ! intent(out)
                                   , zero_rk4_aux          ! ! intent(in)
   use ed_state_vars        , only : sitetype              & ! structure
                                   , patchtype             & ! structure
                                   , polygontype
   use therm_lib8           , only : qtk8                  ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)  , target     :: initp            ! RK4 structure, intermediate step
   type(rk4patchtype)  , target     :: dinitp           ! RK4 structure, derivatives
   type(sitetype)      , target     :: csite            ! Current site (before integration)
   integer             , intent(in) :: ipa              ! Current patch number
   integer             , intent(in) :: mzg              ! Number of ground layers
   integer             , intent(in) :: mzs              ! Number of snow/ponding layers
   !----- Local variables -----------------------------------------------------------------!
   integer                          :: k                ! Level counter
   integer                          :: k1               ! Level counter
   integer                          :: k2               ! Level counter
   integer                          :: klsl             ! Alias for rk4site%lsl
   integer                          :: kben             ! Alias for layer beneath bottom
   integer                          :: ksn              ! # of temporary water/snow layers
   integer                          :: nsoil            ! Short for csite%soil_text(k,ipa)
   real(kind=8)                     :: wgpfrac          ! Fractional soil moisture
   real(kind=8)                     :: soilcond         ! Soil conductivity
   real(kind=8)                     :: snden            ! Snow/water density
   real(kind=8)                     :: hflxgc           ! Ground -> canopy heat flux
   real(kind=8)                     :: wflxgc           ! Ground -> canopy water flux
   real(kind=8)                     :: qwflxgc          ! Ground -> canopy latent heat flux
   real(kind=8)                     :: dewgnd           ! Dew/frost flux to ground
   real(kind=8)                     :: qdewgnd          ! Dew/frost heat flux to ground
   real(kind=8)                     :: ddewgnd          ! Dew/frost density flux to ground
   real(kind=8)                     :: wshed_tot        ! Water shedding flux
   real(kind=8)                     :: qwshed_tot       ! Energy flux due to water shedding
   real(kind=8)                     :: dwshed_tot       ! Depth flux due to water shedding
   real(kind=8)                     :: throughfall_tot  ! Water shedding flux
   real(kind=8)                     :: qthroughfall_tot ! Energy flux due to water shedding
   real(kind=8)                     :: dthroughfall_tot ! Depth flux due to water shedding
   real(kind=8)                     :: wgpmid           ! Soil in between layers
   real(kind=8)                     :: wloss            ! Water loss due to transpiration
   real(kind=8)                     :: qwloss           ! Energy loss due to transpiration
   real(kind=8)                     :: dqwt             ! Energy adjustment aux. variable
   real(kind=8)                     :: fracliq          ! Fraction of liquid water
   real(kind=8)                     :: tempk            ! Temperature
   real(kind=8)                     :: qwgoal           ! Goal energy for thin snow layers.
   real(kind=8)                     :: wprevious        ! Previous water content
   real(kind=8)                     :: infilt           ! Surface infiltration rate
   real(kind=8)                     :: qinfilt          ! Surface infiltration heat rate
   real(kind=8)                     :: snowdens         ! Snow density (kg/m2)
   real(kind=8)                     :: soilhcap         ! Soil heat capacity
   real(kind=8)                     :: int_sfcw_u       ! Intensive sfc. water internal en.
   real(kind=8)                     :: surface_water    ! Temp. variable. Available liquid 
                                                        !   water on the soil sfc (kg/m2)
   real(kind=8)                     :: avg_soil_fracliq ! Avg. fraction of liquid water
   real(kind=8)                     :: freezecor        ! Correction to conductivity for 
                                                        !    partially frozen soil.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! Depending on the type of compilation, interfaces must be explicitly declared.         !
   !---------------------------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine canopy_derivs_two(mzg,initp,dinitp,csite,ipa,hflxgc,wflxgc,qwflxgc        &
                                  ,dewgndflx,qdewgndflx,ddewgndflx,throughfall_tot         &
                                  ,qthroughfall_tot,dthroughfall_tot,wshed_tot,qwshed_tot  &
                                  ,dwshed_tot)
         use rk4_coms     , only: rk4patchtype  ! ! structure
         use ed_state_vars, only: sitetype      & ! structure
                                , patchtype     & ! structure
                                , polygontype   ! ! structure
         implicit none
         type (rk4patchtype), target      :: initp, dinitp
         type (sitetype)    , target      :: csite
         integer            , intent(in)  :: ipa
         integer            , intent(in)  :: mzg
         real(kind=8)       , intent(out) :: hflxgc
         real(kind=8)       , intent(out) :: wflxgc
         real(kind=8)       , intent(out) :: qwflxgc
         real(kind=8)       , intent(out) :: dewgndflx
         real(kind=8)       , intent(out) :: qdewgndflx
         real(kind=8)       , intent(out) :: ddewgndflx
         real(kind=8)       , intent(out) :: wshed_tot
         real(kind=8)       , intent(out) :: qwshed_tot
         real(kind=8)       , intent(out) :: dwshed_tot
         real(kind=8)       , intent(out) :: throughfall_tot
         real(kind=8)       , intent(out) :: qthroughfall_tot
         real(kind=8)       , intent(out) :: dthroughfall_tot
      end subroutine canopy_derivs_two
   end interface
#endif
   !---------------------------------------------------------------------------------------!

   !----- Copy the # of surface water/snow layers and bottom layer to aliases -------------!
   ksn  = initp%nlev_sfcwater
   klsl = rk4site%lsl
   kben = klsl - 1
  
   !---- Flush auxiliary variables to zero. -----------------------------------------------!
   call zero_rk4_aux()

   !----- Make sure derivatives are flushed to zero. --------------------------------------!
   dinitp%soil_energy(:)     = 0.0d0
   dinitp%soil_water(:)      = 0.0d0
   dinitp%sfcwater_depth(:)  = 0.0d0
   dinitp%sfcwater_energy(:) = 0.0d0
   dinitp%sfcwater_mass(:)   = 0.0d0
   dinitp%virtual_energy     = 0.0d0
   dinitp%virtual_water      = 0.0d0
   dinitp%virtual_depth      = 0.0d0
   dinitp%avg_transloss(:)   = 0.0d0
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Calculate water available to vegetation (in meters). SLZ is specified in RAMSIN.  !
   ! Each element of the array sets the value of the bottom of a corresponding soil layer. !
   ! Eg, SLZ = -2, -1, -0.5, -0.25.  There are four soil layers in this example; soil      !
   ! layer 1 goes from 2 meters below the surface to 1 meter below the surface.            !
   !---------------------------------------------------------------------------------------!
   nsoil = rk4site%ntext_soil(mzg)
   rk4aux%available_liquid_water(mzg) = dslz8(mzg)                                         &
                                      * max( 0.0d0                                         &
                                           , initp%soil_fracliq(mzg)                       &
                                           * (initp%soil_water(mzg)-soil8(nsoil)%soilwp) )
   do k = mzg - 1, klsl, -1
      nsoil = rk4site%ntext_soil(k)
      rk4aux%available_liquid_water(k) = rk4aux%available_liquid_water(k+1) + dslz8(k)     &
                                       * max( 0.0d0                                        &
                                            , initp%soil_fracliq(k)                        &
                                            * (initp%soil_water(k) - soil8(nsoil)%soilwp) )
   end do
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !     Compute the following variables:                                                  !
   ! - psiplusz      - the total potential (water + gravitational)   [     m]              !
   ! - hydcond       - hydraulic conductivity                        [   m/s]              !
   ! - soil_liq_wilt - the liquid water content above wilting point  [ m³/m³]              !
   ! - soil_liq_dry  - the liquid water content above dry air soil   [ m³/m³]              !
   !---------------------------------------------------------------------------------------!
   do k = klsl, mzg
      nsoil                   = rk4site%ntext_soil(k)
      wgpfrac                 = min(initp%soil_water(k)/soil8(nsoil)%slmsts, 1.d0)
      rk4aux%hydcond      (k) = slcons18(k,nsoil)                                          &
                              * wgpfrac ** (2.d0 * soil8(nsoil)%slbs + 3.d0)
      rk4aux%psiplusz     (k) = slzt8(k) + soil8(nsoil)%slpots                             &
                                         / wgpfrac ** soil8(nsoil)%slbs
      rk4aux%soil_liq_wilt(k) = max(0.d0, (initp%soil_water(k) - soil8(nsoil)%soilwp)      &
                                        * initp%soil_fracliq(k)                       )
      rk4aux%drysoil      (k) = (initp%soil_water(k) - soil8(nsoil)%soilcp)                &
                              * initp%soil_fracliq(k)                        <= 0.d0
      rk4aux%satsoil      (k) = initp%soil_water(k) >= soil8(nsoil)%slmsts
   end do
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !    Find the boundary condition for total potential beneath the bottom layer.          !
   !---------------------------------------------------------------------------------------!
   nsoil = rk4site%ntext_soil(klsl)
   select case (isoilbc)
   case (0)
      !------------------------------------------------------------------------------------!
      !    Bedrock.  Make the potential exactly the same as the bottom layer, and the flux !
      ! will be zero.                                                                      !
      !------------------------------------------------------------------------------------!
      initp%soil_water   (kben) = initp%soil_water   (klsl)
      initp%soil_fracliq (kben) = initp%soil_fracliq (klsl)
      rk4aux%hydcond     (kben) = rk4aux%hydcond     (klsl)
      rk4aux%psiplusz    (kben) = rk4aux%psiplusz    (klsl)
      rk4aux%drysoil     (kben) = .true.
      rk4aux%satsoil     (kben) = .true.
      !------------------------------------------------------------------------------------!

   case (1)
      !------------------------------------------------------------------------------------!
      !     Free drainage.   Make the water potential at the layer beneath to be at the    !
      ! same soil moisture as the bottom layer.                                            !
      !------------------------------------------------------------------------------------!
      initp%soil_water   (kben) = initp%soil_water   (klsl)
      initp%soil_fracliq (kben) = initp%soil_fracliq (klsl)

      wgpfrac                   = min(initp%soil_water(kben)/soil8(nsoil)%slmsts, 1.d0)
      rk4aux%hydcond     (kben) = rk4aux%hydcond     (klsl)
      rk4aux%psiplusz    (kben) = slzt8(kben) + soil8(nsoil)%slpots                        &
                                              / wgpfrac ** soil8(nsoil)%slbs
      rk4aux%drysoil     (kben) = .false.
      rk4aux%satsoil     (kben) = .false.
      !------------------------------------------------------------------------------------!

   case (2)
      !------------------------------------------------------------------------------------!
      !     Super drainage.   Make the water potential at the layer beneath to be at dry   !
      ! air soil, and find the corresponding hydraulic conductivity.                       !
      !------------------------------------------------------------------------------------!
      initp%soil_water   (kben) = soil8(nsoil)%soilcp
      initp%soil_fracliq (kben) = initp%soil_fracliq (klsl)

      wgpfrac                   = min(initp%soil_water(kben)/soil8(nsoil)%slmsts, 1.d0)
      rk4aux%hydcond     (kben) = slcons18(kben,nsoil)                                     &
                                * wgpfrac ** (2.d0 * soil8(nsoil)%slbs + 3.d0)
      rk4aux%psiplusz    (kben) = slzt8(kben) + soil8(nsoil)%slpots                        &
                                              / wgpfrac ** soil8(nsoil)%slbs
      rk4aux%drysoil     (kben) = .false.
      rk4aux%satsoil     (kben) = .false.
      !------------------------------------------------------------------------------------!

   case (3)
      !------------------------------------------------------------------------------------!
      !     Make the soil moisture in the layer beneath to be at field capacity.           !
      !------------------------------------------------------------------------------------!
      initp%soil_water   (kben) = soil8(nsoil)%slmsts
      initp%soil_fracliq (kben) = initp%soil_fracliq (klsl)

      wgpfrac                   = min(initp%soil_water(kben)/soil8(nsoil)%sfldcap, 1.d0)
      rk4aux%hydcond     (kben) = slcons18(kben,nsoil)                                     &
                                * wgpfrac ** (2.d0 * soil8(nsoil)%slbs + 3.d0)
      rk4aux%psiplusz    (kben) = slzt8(kben) + soil8(nsoil)%slpots                        &
                                              / wgpfrac ** soil8(nsoil)%slbs
      rk4aux%drysoil     (kben) = .false.
      rk4aux%satsoil     (kben) = .false.

   case (4)
      !------------------------------------------------------------------------------------!
      !     Aquifer.   Make the soil moisture in the layer beneath to be always saturated. !
      !------------------------------------------------------------------------------------!
      initp%soil_water   (kben) = soil8(nsoil)%slmsts
      initp%soil_fracliq (kben) = initp%soil_fracliq (klsl)

      wgpfrac                   = min(initp%soil_water(kben)/soil8(nsoil)%slmsts, 1.d0)
      rk4aux%hydcond     (kben) = slcons18(kben,nsoil)                                     &
                                * wgpfrac ** (2.d0 * soil8(nsoil)%slbs + 3.d0)
      rk4aux%psiplusz    (kben) = slzt8(kben) + soil8(nsoil)%slpots                        &
                                              / wgpfrac ** soil8(nsoil)%slbs
      rk4aux%drysoil     (kben) = .false.
      rk4aux%satsoil     (kben) = .false.
   end select
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !      Get derivatives of vegetation and canopy air space variables, plus some fluxes   !
   ! that will be used for soil top boundary conditions and for transpiration.             !
   !---------------------------------------------------------------------------------------!
   call canopy_derivs_two(mzg,initp,dinitp,csite,ipa,hflxgc,wflxgc,qwflxgc,dewgnd,qdewgnd  &
                         ,ddewgnd,throughfall_tot,qthroughfall_tot,dthroughfall_tot        &
                         ,wshed_tot,qwshed_tot,dwshed_tot)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Find the conductivity and resistance.  The if inside is to prevent division by   !
   ! zero, as bedrock porosity is zero.                                                    !
   !---------------------------------------------------------------------------------------!
   !------ Soil layers. -------------------------------------------------------------------!
   do k = klsl, mzg
      nsoil = rk4site%ntext_soil(k)
      if(nsoil /= 13)then
         wgpfrac  = min(initp%soil_water(k)/soil8(nsoil)%slmsts,1.d0)
         soilcond = soil8(nsoil)%soilcond0                                                 &
                  + wgpfrac * (soil8(nsoil)%soilcond1 + wgpfrac *  soil8(nsoil)%soilcond2)
      else
         soilcond=soil8(nsoil)%soilcond0
      end if
      rk4aux%rfactor(k) = dslz8(k) / soilcond
   end do
   !----- Snow/surface water --------------------------------------------------------------!
   do k = 1, ksn
      if(initp%sfcwater_depth(k) > 0.d0)then
         snden = initp%sfcwater_mass(k) / initp%sfcwater_depth(k)
         rk4aux%rfactor(k+mzg) = initp%sfcwater_depth(k)                                   &
                               / (ss(1) * exp(ss(2) * initp%sfcwater_tempk(k))             &
                               * (ss(3) + snden * (ss(4) + snden * (ss(5) + snden*ss(6)))))
      else
         rk4aux%rfactor(k+mzg) = 0.d0
      end if
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Calculate the sensible heat fluxes find soil and sfcwater internal sensible heat  !
   ! fluxes (hfluxgsc) [W/m2].  The convention is that hfluxgsc is positive when the flux  !
   ! is upwards.                                                                           !
   !---------------------------------------------------------------------------------------!
   do k = klsl+1, mzg
      rk4aux%hfluxgsc(k)          = - (initp%soil_tempk(k) - initp%soil_tempk(k-1))        &
                                    / ((rk4aux%rfactor(k) + rk4aux%rfactor(k-1)) * 5.d-1)
      !------ Diagnostic sensible heat flux. ----------------------------------------------!
      dinitp%avg_sensible_gg(k-1) = rk4aux%hfluxgsc(k)
   end do
   !----- If temporary water/snow layers exist, compute them now... -----------------------!
   if (ksn >= 1) then
      rk4aux%hfluxgsc(mzg+1)   = - (initp%sfcwater_tempk(1) - initp%soil_tempk(mzg))       &
                               / ((rk4aux%rfactor(mzg+1)   + rk4aux%rfactor(mzg)) * 5.d-1)
      do k = 2,ksn
         rk4aux%hfluxgsc(mzg+k) = - (initp%sfcwater_tempk(k) - initp%sfcwater_tempk(k-1))  &
                                  / ( (rk4aux%rfactor(mzg+k) + rk4aux%rfactor(mzg+k-1))    &
                                      * 5.d-1)
      end do
   end if
   !----- Heat flux (hfluxgsc) at soil or sfcwater top from longwave, sensible [W/m^2] ----!
   rk4aux%hfluxgsc(mzg+ksn+1)  = hflxgc + qwflxgc                                          &
                               - dble(csite%rlong_g(ipa)) - dble(csite%rlong_s(ipa))
   !----- Diagnostic sensible heat flux ---------------------------------------------------!
   dinitp%avg_sensible_gg(mzg) = rk4aux%hfluxgsc(mzg+ksn+1)
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !    Update soil U values [J/m³] from sensible heat, upward water vapor (latent heat)   !
   ! and longwave fluxes. This excludes effects of dew/frost formation, precipitation,     !
   ! shedding, and percolation.                                                            !
   !---------------------------------------------------------------------------------------!
   do k = klsl,mzg
      dinitp%soil_energy(k) = dslzi8(k) * (rk4aux%hfluxgsc(k) - rk4aux%hfluxgsc(k+1))
   end do
   !----- Update soil Q values [J/m³] from shortwave flux. --------------------------------!
   dinitp%soil_energy(mzg) = dinitp%soil_energy(mzg)                                       &
                           + dslzi8(mzg) * dble(csite%rshort_g(ipa))
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Update surface water U values [J/m²] from sensible heat, upward water vapor        !
   ! (latent heat), longwave, and shortwave fluxes.  This excludes effects of dew/frost    !
   ! formation, precipitation, shedding and percolation.                                   !
   !---------------------------------------------------------------------------------------!
   do k = 1,ksn
     dinitp%sfcwater_energy(k) = rk4aux%hfluxgsc(k+mzg) - rk4aux%hfluxgsc(k+1+mzg)         &
                               + dble(csite%rshort_s(k,ipa))
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Calculate the fluxes of water with their associated heat fluxes. Update top soil  !
   ! or snow moisture from evaporation only.                                               !
   !                                                                                       !
   !     New moisture, qw, and depth from dew/frost formation, precipitation, shedding,    !
   ! and percolation.  ksnnew is the layer that receives the new condensate that comes     !
   ! directly from the air above.  If there is no pre-existing snowcover, this is a        !
   ! temporary "snow" layer.                                                               !
   !---------------------------------------------------------------------------------------!
   rk4aux%w_flux(mzg+ksn+1)  = -  dewgnd -  wshed_tot -  throughfall_tot
   rk4aux%qw_flux(mzg+ksn+1) = - qdewgnd - qwshed_tot - qthroughfall_tot
   rk4aux%d_flux(ksn+1)      = - ddewgnd - dwshed_tot - dthroughfall_tot
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Transfer water downward through snow layers by percolation. Here we define:       !
   ! + fracliq: the fraction of liquid in the snowcover or surface water;                  !
   ! + wfree:   the quantity of that liquid in kg/m2 which is free (i.e. not attached to   !
   !            snowcover) and therefore available to soak into the layer below;           !
   ! + soilcap: the capacity of the top soil layer in kg/m2 to accept surface water.       !
   ! + flxliq:  the effective soaking flux in kg/m2/s over the timestep which is used to   !
   !            update the moisture in the top soil layer.                                 !
   !---------------------------------------------------------------------------------------!
   if (ksn > 0) then
      dinitp%sfcwater_mass(ksn)   = - rk4aux%w_flux(mzg+ksn+1) - wflxgc
      dinitp%sfcwater_energy(ksn) = dinitp%sfcwater_energy(ksn) - rk4aux%qw_flux(mzg+ksn+1)
      dinitp%sfcwater_depth(ksn)  = -rk4aux%d_flux(ksn+1)
   else if (rk4aux%w_flux(mzg+1) < 0.d0) then
      dinitp%virtual_energy       = - rk4aux%qw_flux(mzg+1)
      dinitp%virtual_water        = - rk4aux%w_flux(mzg+1)
      dinitp%virtual_depth        = - rk4aux%d_flux(1)
      rk4aux%qw_flux(mzg+1)       = 0.d0
      rk4aux%w_flux(mzg+1)        = wflxgc
   else
      rk4aux%w_flux(mzg+1)        = rk4aux%w_flux(mzg+1) + wflxgc
   end if
   !------ Diagnostic variable for water flux, bypass the virtual/sfcw layers. ------------!
   dinitp%avg_smoist_gg(mzg) = rk4aux%w_flux(mzg+ksn+1)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find amount of water transferred between soil layers (w_flux) [m] modulated by    !
   ! the liquid water fraction.                                                            !
   !---------------------------------------------------------------------------------------!
   rk4aux%w_flux(mzg+1) = rk4aux%w_flux(mzg+1) * wdnsi8 ! now in m/s
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !     Alternate surface infiltration (MCD) based on surface conductivity not capacity.  !
   !---------------------------------------------------------------------------------------!
   if (infiltration_method /= 0) then
      call fatal_error ('Running alt infiltation when we shouldn''t be'                    &
                       ,'leaftw_derivs','rk4_derivs.F90')

      if (initp%virtual_water /= 0.d0) then  !!process "virtural water" pool
         nsoil = rk4site%ntext_soil(mzg)
         if (nsoil /= 13) then
            infilt = -dslzi8(mzg)* 5.d-1 * slcons18(mzg,nsoil)                             &
                   * (initp%soil_water(mzg) / soil8(nsoil)%slmsts)                         &
                     **(2.d0 * soil8(nsoil)%slbs + 3.d0)                                   &
                   * (rk4aux%psiplusz(mzg)-initp%virtual_water/2.d3) & !diff. in pot.
                   * 5.d-1 * (initp%soil_fracliq(mzg)+ initp%virtual_fracliq) ! mean liquid fraction
            qinfilt = infilt * cliqvlme8 * (initp%virtual_tempk - tsupercool8)
            !----- Adjust other rates accordingly -----------------------------------------!
            rk4aux%w_flux(mzg+1)  = rk4aux%w_flux(mzg+1) + infilt
            rk4aux%qw_flux(mzg+1) = rk4aux%qw_flux(mzg+1)+ qinfilt
            dinitp%virtual_water  = dinitp%virtual_water  - infilt*wdns8
            dinitp%virtual_energy = dinitp%virtual_energy - qinfilt
         end if
      end if  !! end virtual water pool
      if (initp%nlev_sfcwater >= 1) then !----- Process "snow" water pool -----------------! 
         surface_water = initp%sfcwater_mass(1)*initp%sfcwater_fracliq(1)*wdnsi8 !(m/m2)
         nsoil = rk4site%ntext_soil(mzg)
         if (nsoil /= 13) then
            !----- Calculate infiltration rate (m/s) --------------------------------------!
            infilt = -dslzi8(mzg) * 5.d-1 * slcons18(mzg,nsoil)                            &
                   * (initp%soil_water(mzg) / soil8(nsoil)%slmsts)                         &
                     **(2.d0 * soil8(nsoil)%slbs + 3.d0)                                   &
                   * (rk4aux%psiplusz(mzg) - surface_water/2.d0) & !difference in potentials
                   * 5.d-1 * (initp%soil_fracliq(mzg) + initp%sfcwater_fracliq(1))
            qinfilt = infilt * cliqvlme8 * (initp%sfcwater_tempk(1) - tsupercool8)
            !----- Adjust other rates accordingly -----------------------------------------!
            rk4aux%w_flux(mzg+1)      = rk4aux%w_flux(mzg+1)      + infilt
            rk4aux%qw_flux(mzg+1)     = rk4aux%qw_flux(mzg+1)     + qinfilt 
            dinitp%sfcwater_mass(1)   = dinitp%sfcwater_mass(1)   - infilt*wdns8
            dinitp%sfcwater_energy(1) = dinitp%sfcwater_energy(1) - qinfilt
            dinitp%sfcwater_depth(1)  = dinitp%sfcwater_depth(1)  - infilt
         end if
      end if  ! End snow water pool
   end if  !! End alternate infiltration
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Here we solve for all layers including the bottommost one, which is now prepared  !
   ! according to the chosen boundary condition.  In this block, units are:                !
   ! W_Flux  = m/s                                                                         !
   ! QW_Flux = W/m2                                                                        !
   !---------------------------------------------------------------------------------------!
   do k = klsl, mzg
      nsoil = rk4site%ntext_soil(k)
      if (nsoil /= 13) then

         !----- Find the correction for (partially) frozen soil layers. -------------------!
         avg_soil_fracliq = 5.d-1 * (initp%soil_fracliq(k) + initp%soil_fracliq(k-1))
         freezecor        = 1.d1 ** ( max( lnexp_min8                                      &
                                         , - freezecoef8 * (1.d0 - avg_soil_fracliq)) )
         !---------------------------------------------------------------------------------!



         !----- Find the potential flux. --------------------------------------------------!
         rk4aux%w_flux(k) = - freezecor                                                    &
                            * 5.d-1 * (rk4aux%hydcond(k-1) + rk4aux%hydcond(k))            &
                            * (rk4aux%psiplusz(k) - rk4aux%psiplusz(k-1)) * dslzti8(k)
         !---------------------------------------------------------------------------------!



         !----- Limit water transfers to prevent over-saturation and over-depletion. ------!
         if ( rk4aux%w_flux(k) >= 0. .and.                                                 &
             (rk4aux%drysoil(k-1) .or. rk4aux%satsoil(k)) )    then
            rk4aux%w_flux(k) = 0.d0

         elseif( rk4aux%w_flux(k) < 0. .and.                                               &
                (rk4aux%satsoil(k-1) .or. rk4aux%drysoil(k)) ) then
            rk4aux%w_flux(k) = 0.d0

         end if
         !---------------------------------------------------------------------------------!

      else
         rk4aux%w_flux(k) = 0.d0
      end if
      !----- Only liquid water is allowed to flow, find qw_flux (W/m2) accordingly. -------!
      rk4aux%qw_flux(k)         = rk4aux%w_flux(k) * cliqvlme8                             &
                                * (initp%soil_tempk(k) - tsupercool8)
      !----- Save the moisture flux in kg/m2/s. -------------------------------------------!
      if (k /= 1) dinitp%avg_smoist_gg(k-1) = rk4aux%w_flux(k) * wdns8   ! Diagnostic
   end do

   !---------------------------------------------------------------------------------------!
   !     Find the drainage flux.  This is simply the flux at the bottom of the bottom-     !
   ! most layer.  Units are kg/m2/s for water flux and W/m2 for the associated internal    !
   ! energy leak.  Notice that for the aquifer case we may actually have negative drain-   !
   ! age, but that shouldn't affect the budget in any way (except that we are adding water !
   ! to the system).                                                                       !
   !---------------------------------------------------------------------------------------!
   dinitp%avg_drainage      = - rk4aux%w_flux(klsl) * wdns8
   dinitp%avg_drainage_heat = - rk4aux%qw_flux(klsl)
   !----- Copy the variables to the budget arrays. ----------------------------------------!
   if (checkbudget) then
      dinitp%wbudget_loss2drainage = dinitp%avg_drainage
      dinitp%ebudget_loss2drainage = dinitp%avg_drainage_heat

      dinitp%wbudget_storage = dinitp%wbudget_storage - dinitp%avg_drainage
      dinitp%ebudget_storage = dinitp%ebudget_storage - dinitp%avg_drainage_heat
   end if
   !---------------------------------------------------------------------------------------!




   !----- Finally, update soil moisture and soil energy. ----------------------------------!
   do k = klsl,mzg
      dinitp%soil_water(k)  = dinitp%soil_water(k)                                         &
                            + dslzi8(k) * (  rk4aux%w_flux(k)  -  rk4aux%w_flux(k+1) )
      dinitp%soil_energy(k) = dinitp%soil_energy(k)                                        &
                            + dslzi8(k) * ( rk4aux%qw_flux(k)  - rk4aux%qw_flux(k+1) )
   end do
   !---------------------------------------------------------------------------------------!




   !---- Update soil moisture and energy from transpiration/root uptake. ------------------!
   if (any_resolvable) then
      do k1 = klsl, mzg    ! loop over extracted water
         do k2=k1,mzg
            if (rk4site%ntext_soil(k2) /= 13) then
               if (rk4aux%available_liquid_water(k1) > 0.d0) then
                  wloss = wdnsi8 * rk4aux%extracted_water(k1)                              &
                        * rk4aux%soil_liq_wilt(k2) / rk4aux%available_liquid_water(k1)
                  dinitp%soil_water(k2) = dinitp%soil_water(k2) - dble(wloss)

                  !----- Energy: only liquid water is lost through transpiration. ---------!
                  qwloss = wloss * cliqvlme8 * (initp%soil_tempk(k2) - tsupercool8)
                  dinitp%soil_energy(k2)   = dinitp%soil_energy(k2)   - qwloss
                  dinitp%avg_transloss(k2) = dinitp%avg_transloss(k2) - wdns8*dble(wloss)  &
                             * dslz8(k2)
               end if
            end if
         end do
      end do
   end if
   !---------------------------------------------------------------------------------------!


   return
end subroutine leaftw_derivs
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine canopy_derivs_two(mzg,initp,dinitp,csite,ipa,hflxgc,wflxgc,qwflxgc,dewgndflx    &
                            ,qdewgndflx,ddewgndflx,throughfall_tot,qthroughfall_tot        &
                            ,dthroughfall_tot,wshed_tot,qwshed_tot,dwshed_tot)
   use rk4_coms              , only : rk4patchtype         & ! Structure
                                    , rk4site              & ! intent(in)
                                    , rk4aux               & ! intent(inout)
                                    , toocold              & ! intent(in)
                                    , toohot               & ! intent(in)
                                    , lai_to_cover         & ! intent(in)
                                    , effarea_heat         & ! intent(in)
                                    , effarea_evap         & ! intent(in)
                                    , effarea_transp       & ! intent(in)
                                    , zoveg                & ! intent(in)
                                    , zveg                 & ! intent(in)
                                    , wcapcan              & ! intent(in)
                                    , wcapcani             & ! intent(in)
                                    , hcapcani             & ! intent(in)
                                    , ccapcani             & ! intent(in)
                                    , any_resolvable       & ! intent(in)
                                    , tiny_offset          & ! intent(in)
                                    , rk4leaf_drywhc       & ! intent(in)
                                    , rk4leaf_maxwhc       & ! intent(in)
                                    , checkbudget          & ! intent(in)
                                    , print_detailed       & ! intent(in)
                                    , supersat_ok          & ! intent(in)
                                    , leaf_intercept       ! ! intent(in)
   use ed_state_vars         , only : sitetype             & ! Structure
                                    , patchtype            & ! Structure
                                    , polygontype
   use consts_coms           , only : alvl8                & ! intent(in)
                                    , cp8                  & ! intent(in)
                                    , cpi8                 & ! intent(in)
                                    , twothirds8           & ! intent(in)
                                    , day_sec8             & ! intent(in)
                                    , grav8                & ! intent(in)
                                    , alvi8                & ! intent(in)
                                    , alvl8                & ! intent(in)
                                    , alli8                & ! intent(in)
                                    , umol_2_kgC8          & ! intent(in)
                                    , pi18                 & ! intent(in)
                                    , halfpi8              & ! intent(in)
                                    , mmdry8               & ! intent(in)
                                    , mmdryi8              & ! intent(in)
                                    , wdns8                & ! intent(in)
                                    , wdnsi8               & ! intent(in)
                                    , fdnsi8               & ! intent(in)
                                    , t3ple8               & ! intent(in)
                                    , tsupercool8          & ! intent(in)
                                    , cice8                & ! intent(in)
                                    , cliq8                & ! intent(in)
                                    , epi8                 & ! intent(in)
                                    , huge_num8            ! ! intent(in)
   use soil_coms             , only : soil8                & ! intent(in)
                                    , dslzi8               & ! intent(in)
                                    , dewmax               ! ! intent(in)
   use therm_lib8            , only : rslif8               ! ! function
   use ed_misc_coms          , only : dtlsm                & ! intent(in)
                                    , fast_diagnostics     ! ! intent(in)
   use canopy_struct_dynamics, only : vertical_vel_flux8   ! ! function
   use pft_coms              , only : water_conductance    ! ! intent(in)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)     , target      :: csite             ! Current site
   type(rk4patchtype) , target      :: initp             ! RK4 structure, state vars
   type(rk4patchtype) , target      :: dinitp            ! RK4 structure, derivatives
   integer            , intent(in)  :: ipa               ! Current patch ID
   integer            , intent(in)  :: mzg               ! Current patch ID
   real(kind=8)       , intent(out) :: hflxgc            ! Ground->canopy sens. heat flux
   real(kind=8)       , intent(out) :: wflxgc            ! Ground->canopy water flux
   real(kind=8)       , intent(out) :: qwflxgc           ! Ground->canopy latent heat flux
   real(kind=8)       , intent(out) :: dewgndflx         ! Dew/frost water flux
   real(kind=8)       , intent(out) :: qdewgndflx        ! Dew/frost heat flux
   real(kind=8)       , intent(out) :: ddewgndflx        ! Dew/frost density
   real(kind=8)       , intent(out) :: throughfall_tot   ! Throughfall rate
   real(kind=8)       , intent(out) :: qthroughfall_tot  ! Throughfall energy flux
   real(kind=8)       , intent(out) :: dthroughfall_tot  ! Throughfall depth flux
   real(kind=8)       , intent(out) :: wshed_tot         ! Water shed from leaves
   real(kind=8)       , intent(out) :: qwshed_tot        ! Internal energy of water shed
   real(kind=8)       , intent(out) :: dwshed_tot        ! Depth of water shed
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)    , pointer     :: cpatch            ! Current patch
   integer                          :: ico               ! Current cohort ID
   integer                          :: k                 ! Soil layer counter
   integer                          :: ipft              ! Shortcut for PFT type
   integer                          :: kroot             ! Level of the bottom of root is
   real(kind=8)                     :: closedcan_frac    ! total fractional canopy coverage
   real(kind=8)                     :: transp            ! Cohort transpiration
   real(kind=8)                     :: cflxac            ! Atm->canopy carbon flux
   real(kind=8)                     :: wflxac            ! Atm->canopy water flux
   real(kind=8)                     :: hflxac            ! Atm->canopy sensible heat flux
   real(kind=8)                     :: eflxac            ! Atm->canopy Eq. Pot. temp flux
   real(kind=8)                     :: wflxlc_try        ! Intended flux leaf sfc -> canopy
   real(kind=8)                     :: wflxwc_try        ! Intended flux wood sfc -> canopy
   real(kind=8)                     :: c3lai             ! Term for psi_open/psi_closed
   real(kind=8)                     :: hflxlc            ! Leaf->canopy heat flux
   real(kind=8)                     :: hflxwc            ! Wood->canopy heat flux
   real(kind=8)                     :: rgnd              !
   real(kind=8)                     :: sigmaw            !
   real(kind=8)                     :: wflxlc            ! Leaf sfc -> canopy water flux
   real(kind=8)                     :: wflxwc            ! Wood sfc -> canopy water flux
   real(kind=8)                     :: cflxgc            !
   real(kind=8)                     :: wshed             ! Water shed from leaves
   real(kind=8)                     :: qwshed            ! Internal energy of water shed
   real(kind=8)                     :: dwshed            ! Depth of water shed
   real(kind=8)                     :: throughfall       ! Extra throughfall from full coh.
   real(kind=8)                     :: qthroughfall      ! Its internal energy
   real(kind=8)                     :: dthroughfall      ! Its depth
   real(kind=8)                     :: taii              !
   real(kind=8)                     :: wflx              !
   real(kind=8)                     :: wood_shv          ! Sat. sp. hum. at wood sfc.
   real(kind=8)                     :: hflxlc_tot        ! Total leaf -> CAS heat flux
   real(kind=8)                     :: hflxwc_tot        ! Total wood -> CAS heat flux
   real(kind=8)                     :: transp_tot        ! Total transpiration (water)
   real(kind=8)                     :: qtransp_tot       ! Total transpiration (energy)
   real(kind=8)                     :: cflxlc_tot        ! Total leaf -> CAS CO2 flux
   real(kind=8)                     :: wflxlc_tot        ! Leaf -> CAS evaporation (water)
   real(kind=8)                     :: wflxwc_tot        ! Wood -> CAS evaporation (water)
   real(kind=8)                     :: qwflxlc_tot       ! Leaf -> CAS evaporation (energy)
   real(kind=8)                     :: qwflxwc_tot       ! Wood -> CAS evaporation (energy)
   real(kind=8)                     :: rho_ustar         !
   real(kind=8)                     :: storage_decay     !
   real(kind=8)                     :: leaf_flux         !
   real(kind=8)                     :: min_leaf_water    !
   real(kind=8)                     :: max_leaf_water    !
   real(kind=8)                     :: min_wood_water    !
   real(kind=8)                     :: max_wood_water    !
   real(kind=8)                     :: maxfluxrate       !
   real(kind=8)                     :: intercepted_max   ! Pot. interecepted rainfall
   real(kind=8)                     :: qintercepted_max  ! Int. energy of pot. intercept.
   real(kind=8)                     :: dintercepted_max  ! Depth of pot. interception
   real(kind=8)                     :: intercepted_tot   ! Actual intercepted rainfall
   real(kind=8)                     :: qintercepted_tot  ! Int. energy of act. intercept.
   real(kind=8)                     :: leaf_intercepted  ! Leaf interception
   real(kind=8)                     :: leaf_qintercepted ! Int. energy of leaf intercept.
   real(kind=8)                     :: wood_intercepted  ! Leaf interception
   real(kind=8)                     :: wood_qintercepted ! Int. energy of leaf intercept.
   real(kind=8)                     :: qwflxlc           ! Leaf -> CAS evaporation (energy)
   real(kind=8)                     :: qwflxwc           ! Wood -> CAS evaporation (energy)
   real(kind=8)                     :: qtransp           ! Transpiration (energy)
   real(kind=8)                     :: water_demand      !
   real(kind=8)                     :: water_supply      !
   real(kind=8)                     :: flux_area         ! Area between canopy and plant
   !----- Functions -----------------------------------------------------------------------!
   real        , external           :: sngloff
   !---------------------------------------------------------------------------------------!



   !----- First step, we assign the pointer for the current patch. ------------------------!
   cpatch => csite%patch(ipa)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Computing the fluxes from atmosphere to canopy.                                    !
   !---------------------------------------------------------------------------------------!
   rho_ustar = initp%can_rhos * initp%ustar                          ! Aux. variable
   hflxac    = rho_ustar      * initp%tstar * initp%can_exner        ! Sensible Heat flux
   eflxac    = rho_ustar      * initp%estar * cp8 * initp%can_temp   ! Enthalpy flux
   wflxac    = rho_ustar      * initp%qstar                          ! Water flux
   cflxac    = rho_ustar      * initp%cstar * mmdryi8                ! CO2 flux [umol/m2/s]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Calculate fraction of closed canopy.                                              !
   !---------------------------------------------------------------------------------------!
   closedcan_frac = max(0.d0,min(1.d0,1.d0-initp%opencan_frac))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Here we will determine the initial guess for throughfall precipitation and the    !
   ! potential amount of interception.  Later we will adjust both values depending on      !
   ! whether some cohorts are already full or not, or if it is fine to exceed the maximum  !
   ! amount of water that a cohort can hold.                                               !
   !---------------------------------------------------------------------------------------!
   if (any_resolvable) then
      taii = 0.d0
      cpatch => csite%patch(ipa)
      do ico = 1,cpatch%ncohorts
         taii = taii + initp%tai(ico)
      end do
      taii = 1.d0/taii

      !------------------------------------------------------------------------------------!
      !    If the canopy does not cover all of the ground, then it should not intercept    !
      ! all of the water.                                                                  !
      !------------------------------------------------------------------------------------!
      if (rk4site%pcpg > 0.d0) then
         if (leaf_intercept) then
            !----- Scale interception by canopy openess (MCD 01-12-09). -------------------!
            intercepted_max  = rk4site%pcpg  * closedcan_frac
            qintercepted_max = rk4site%qpcpg * closedcan_frac
            dintercepted_max = rk4site%dpcpg * closedcan_frac
         else
            !----- No interception (developer only). --------------------------------------!
            intercepted_max  = 0.d0
            qintercepted_max = 0.d0
            dintercepted_max = 0.d0
         end if

         !---------------------------------------------------------------------------------!
         !    The first guess for through fall is the rainfall minus the maximum           !
         ! interception.  If some of the water can't be intercepted, we will add to the    !
         ! through fall later.                                                             !
         !---------------------------------------------------------------------------------!
         throughfall_tot  = rk4site%pcpg  - intercepted_max
         qthroughfall_tot = rk4site%qpcpg - qintercepted_max
         dthroughfall_tot = rk4site%dpcpg - dintercepted_max
      else
         !----- No precipitation, nothing to be intercepted... ----------------------------!
         intercepted_max  = 0.d0
         qintercepted_max = 0.d0
         dintercepted_max = 0.d0
         throughfall_tot  = 0.d0
         qthroughfall_tot = 0.d0
         dthroughfall_tot = 0.d0
      end if

   else
      !------------------------------------------------------------------------------------!
      !     If the TAI is very small or total patch vegetation heat capacity is too        !
      ! small, bypass vegetation computations.  Set throughfall precipitation heat and     !
      ! moisture to unintercepted values.                                                  !
      !------------------------------------------------------------------------------------!
      intercepted_max  = 0.d0
      qintercepted_max = 0.d0
      dintercepted_max = 0.d0
      throughfall_tot  = rk4site%pcpg
      qthroughfall_tot = rk4site%qpcpg
      dthroughfall_tot = rk4site%dpcpg

      !------------------------------------------------------------------------------------!
      ! Note: If the condition of low TAI for the entire patch was met, then it does not   !
      !       matter what the individual cohorts are normalized by, because they are       !
      !       effectively zero. So make sure the inverse patch TAI is a nominal non-zero/  !
      !       non-infinite number. This will only be used when parsing out intercepted     !
      !       leaf water into shed water; in which case the intercepted water is zero      !
      !       anyway. So this is just to prevent FPEs.                                     !
      !------------------------------------------------------------------------------------!
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
   !---------------------------------------------------------------------------------------!
   !     Compute sensible heat and moisture fluxes between the ground and the canopy air   !
   ! space.  The ground may be either the top soil layer or the temporary surface water    !
   ! snow surface.  wflx is the possible vapour flux, and can be either positive or        !
   ! negative.  After we compute it, we check its sign, because the value is good only     !
   ! when it is negative or when the surface is water/snow.  The ground calculation must   !
   ! be done differently.  Later, two variables will define the flux between ground and    !
   ! canopy air space:                                                                     !
   ! - wflxgc [kg/m2/s] is going to be the evaporation flux from ground to canopy.         !
   ! - dewgnd [kg/m2/s] is going to be the dew/frost flux from canopy to ground.           !
   !     Both are defined as positive quantities.  Sensible heat is defined by only one    !
   ! variable, hflxgc [J/m2/s], which can be either positive or negative.                  !
   !---------------------------------------------------------------------------------------!
   hflxgc = initp%ggnet * initp%can_rhos * cp8 * (initp%ground_temp - initp%can_temp)
   wflx   = initp%ggnet * initp%can_rhos       * (initp%ground_ssh  - initp%can_shv )
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !   Here we will decide how to compute the evaporation and condensation fluxes based on !
   ! the sign of wflx, the number of temporary surface water/snow layers, and whether the  !
   ! canopy air is saturation and whether the user is concerned about super-saturation.    !
   !---------------------------------------------------------------------------------------!
   if (wflx <= 0.d0) then
      !------------------------------------------------------------------------------------!
      !     Flux is negative, which means that dew/frost is going to form.  We must also   !
      ! decide whether dew or frost will form, and this is based on the current partition  !
      ! of the surface water phase.  The depth gain uses frost density rather than ice     !
      ! density based on MCD suggestion on 11/16/2009.                                     !
      !------------------------------------------------------------------------------------!
      dewgndflx  = - wflx
      qdewgndflx = dewgndflx * (alvi8 - initp%ground_fliq * alli8)
      ddewgndflx = dewgndflx                                                               &
                 * (initp%ground_fliq * wdnsi8 + (1.d0-initp%ground_fliq) * fdnsi8)
      !----- Set evaporation fluxes to zero. ----------------------------------------------!
      wflxgc     = 0.d0
      qwflxgc    = 0.d0

      !----- Set flux flag. ---------------------------------------------------------------!
      initp%flag_wflxgc = 1

   elseif (initp%can_rhv >= 1.d0 .and. (.not. supersat_ok)) then
      !------------------------------------------------------------------------------------!
      !     In principle evaporation could happen, but the canopy air space is saturated.  !
      ! The user doesn't want too much super-saturation, so we will not let any            !
      ! evaporation to occur.                                                              !
      !------------------------------------------------------------------------------------!
      dewgndflx  = 0.d0
      qdewgndflx = 0.d0
      ddewgndflx = 0.d0
      wflxgc     = 0.d0
      qwflxgc    = 0.d0

      !----- Set flux flag. ---------------------------------------------------------------!
      initp%flag_wflxgc = 2
   elseif (initp%nlev_sfcwater > 0) then
      !------------------------------------------------------------------------------------!
      !     Evaporation will happen, and there is a temporary surface water/snow layer     !
      ! above the soil, so wflx is still good, we simply use it, and make the dew/frost    !
      ! fluxes to be zero.                                                                 !
      !------------------------------------------------------------------------------------!
      dewgndflx  = 0.d0
      qdewgndflx = 0.d0
      ddewgndflx = 0.d0
      wflxgc     = wflx
      qwflxgc    = wflx * (alvi8 - initp%ground_fliq * alli8)

      !----- Set flux flag. ---------------------------------------------------------------!
      initp%flag_wflxgc = 3
   else if (rk4aux%drysoil(mzg)) then
      !------------------------------------------------------------------------------------!
      !   There should be evaporation, except that there is no water left to be extracted  !
      ! from the ground... Set both evaporation and condensation fluxes to zero.           !
      !------------------------------------------------------------------------------------!
      dewgndflx  = 0.d0
      qdewgndflx = 0.d0
      ddewgndflx = 0.d0
      wflxgc     = 0.d0
      qwflxgc    = 0.d0

      !----- Set flux flag. ---------------------------------------------------------------!
      initp%flag_wflxgc = 4
   else
      !------------------------------------------------------------------------------------!
      !    Evaporation will happen, and but water will come from the top most layer.  Wflx !
      ! cannot be used because the ground specific humidity is not the saturation specific !
      ! humidity at the soil temperature only, it depends on the canopy specific humidity  !
      ! itself and the soil moisture.                                                      !
      !------------------------------------------------------------------------------------!
      wflxgc     = initp%ggnet * initp%can_rhos * (initp%ground_shv - initp%can_shv)       &
                 * ( initp%ggsoil / (initp%ggnet + initp%ggsoil) )
      !----- Adjusting the flux accordingly to the surface fraction (no phase bias). ------!
      qwflxgc    = wflxgc * ( alvi8 - initp%ground_fliq * alli8)
      !----- Set condensation fluxes to zero. ---------------------------------------------!
      dewgndflx  = 0.d0
      qdewgndflx = 0.d0
      ddewgndflx = 0.d0

      !----- Set flux flag. ---------------------------------------------------------------!
      initp%flag_wflxgc = 5
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Loop over the cohorts in the patch. Calculate energy fluxes with surrounding      !
   ! canopy air space, integrate cohort energy, calculate precipitation throughfall and    !
   ! sum fluxes to the patch level. Initialize variables used to store sums over cohorts.  !
   !---------------------------------------------------------------------------------------!
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
   !---------------------------------------------------------------------------------------! 




   !---------------------------------------------------------------------------------------! 
   !    Heterotrophic respiration is a patch-level variable, so we initialise the total    !
   ! vegetation to canopy carbon flux with the total heterotrophic respiration due to the  !
   ! coarse wood debris, and remove that from the ground to canopy carbon flux to avoid    !
   ! double counting.                                                                      !
   !---------------------------------------------------------------------------------------! 
   cflxlc_tot       = initp%cwd_rh
   cflxgc           = initp%rh - initp%cwd_rh
   !---------------------------------------------------------------------------------------!
  
   cohortloop: do ico = 1,cpatch%ncohorts

      cflxgc = cflxgc + initp%root_resp(ico)

      !------------------------------------------------------------------------------------!
      !    Calculate 'decay' term of storage.                                              !
      !------------------------------------------------------------------------------------!
      storage_decay = initp%growth_resp(ico) + initp%storage_resp(ico)                     &
                    + initp%vleaf_resp(ico)
      cflxlc_tot    = cflxlc_tot + storage_decay
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! LEAF BUDGET - Check whether the leaves of this cohort haven't been flagged as non- !
      !               resolvable, i.e., it has leaves, belongs to a patch that is not too  !
      !               sparse, an it is not buried in snow.  We should compute energy and   !
      !               water at the cohort level only if the cohort is "safe".  Otherwise,  !
      !               we will set the leaf energy derivatives to zero, and convert any     !
      !               intercepted water into throughfall.  Later on, these "unsafe"        !
      !               cohorts will have their leaf energy set to equilibrium with the      !
      !               canopy air space (temperature).                                      !
      !------------------------------------------------------------------------------------!
      if (initp%leaf_resolvable(ico)) then

         !------ Defining some shortcuts to indices ---------------------------------------!
         ipft  = cpatch%pft(ico)
         kroot = cpatch%krdepth(ico)


         !------  Calculate leaf-level CO2 flux -------------------------------------------!
         leaf_flux = initp%gpp(ico) - initp%leaf_resp(ico)

         !------ Update CO2 flux from vegetation to canopy air space. ---------------------!
         cflxlc_tot = cflxlc_tot - leaf_flux

         !---------------------------------------------------------------------------------!
         !     Define the minimum leaf water to be considered, and the maximum amount      !
         ! possible.                                                                       !
         !---------------------------------------------------------------------------------!
         min_leaf_water = rk4leaf_drywhc * initp%lai(ico)
         max_leaf_water = rk4leaf_maxwhc * initp%lai(ico)

         !------ Calculate fraction of leaves covered with water. -------------------------!
         if (initp%leaf_water(ico) > min_leaf_water) then
            sigmaw = min(1.d0, (initp%leaf_water(ico)/max_leaf_water)**twothirds8)
         else
            sigmaw = 0.d0
         end if

         !---------------------------------------------------------------------------------!
         !    Here we must compute two different areas.  For transpiration, we want the    !
         ! leaf area index only, because we assume transpiration to happen only through    !
         ! leaves.  Evaporation of water/ice settled over the vegetation surface, or dew/  !
         ! frost formation must account the branches and stems as well.                    !
         !---------------------------------------------------------------------------------!
         !----- Evaporation/condensation "flux" -------------------------------------------!
         wflxlc_try = effarea_evap * initp%lai(ico) * initp%leaf_gbw(ico)                  &
                    * (initp%lint_shv(ico) - initp%can_shv)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Computing the evapotranspiration or dew/frost deposition.                    !
         !---------------------------------------------------------------------------------!
         if (wflxlc_try >= 0.d0) then  
            !------------------------------------------------------------------------------!
            !    Probably evapotranspiration, as long as the canopy air is not saturated   !
            ! or the user doesn't mind that super-saturation occur.                        !
            !------------------------------------------------------------------------------!
            if (supersat_ok .or. initp%can_rhv < 1.d0) then
               !---------------------------------------------------------------------------!
               !     Evaporation, energy is scaled by liquid/ice partition (no phase       !
               ! bias).  We scale by the relative area of leaves that is actually covered  !
               ! with water.                                                               !
               !---------------------------------------------------------------------------!
               wflxlc  = wflxlc_try * sigmaw
               qwflxlc = wflxlc * (alvi8 - initp%leaf_fliq(ico) * alli8)

               !---------------------------------------------------------------------------!
               !     Transpiration, consider the one-sided leaf area rather than LAI,      !
               ! when the PFT is amphistomatous, in which case it is double-sided.         !
               ! Compute the water demand from both open closed and open stomata, but      !
               ! first make sure that there is some water available for transpiration...   !
               !---------------------------------------------------------------------------!
               if (rk4aux%available_liquid_water(kroot) > 0.d0 ) then
                  c3lai = effarea_transp(ipft) * initp%lai(ico)                            &
                        * (initp%lint_shv(ico) - initp%can_shv) * initp%leaf_gbw(ico)

                  dinitp%psi_open(ico)   = c3lai * initp%gsw_open(ico)                     &
                                         / (initp%leaf_gbw(ico) + initp%gsw_open(ico))
                  dinitp%psi_closed(ico) = c3lai * initp%gsw_closed(ico)                   &
                                         / (initp%leaf_gbw(ico) + initp%gsw_closed(ico))

                  transp = initp%fs_open(ico) * dinitp%psi_open(ico)                       &
                         + (1.0d0 - initp%fs_open(ico)) * dinitp%psi_closed(ico)
               else
                  dinitp%psi_open(ico)   = 0.d0
                  dinitp%psi_closed(ico) = 0.d0
                  transp                 = 0.d0
               end if
               !---------------------------------------------------------------------------! 
               !    Only liquid water is transpired, thus this is always the condensation  !
               ! latent heat.                                                              !
               !---------------------------------------------------------------------------! 
               qtransp = transp * alvl8
               !---------------------------------------------------------------------------! 

            else
               !----- Canopy is already saturated, no evapotranspiration is allowed. ------!
               wflxlc                 = 0.d0
               qwflxlc                = 0.d0
               transp                 = 0.d0
               qtransp                = 0.d0
               dinitp%psi_open(ico)   = 0.d0
               dinitp%psi_closed(ico) = 0.d0
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!

         else
            !------------------------------------------------------------------------------!
            !     Dew/frost formation. The deposition will conserve the liquid/ice         !
            ! partition (or use the default if there is no water).                         !
            !------------------------------------------------------------------------------!
            wflxlc                 = wflxlc_try
            qwflxlc                = wflxlc * (alvi8 - initp%leaf_fliq(ico)*alli8)
            transp                 = 0.0d0
            qtransp                = 0.0d0
            dinitp%psi_open  (ico) = 0.0d0
            dinitp%psi_closed(ico) = 0.0d0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!



         !----- We need to extract water from the soil equal to the transpiration. --------!
         rk4aux%extracted_water(kroot) = rk4aux%extracted_water(kroot) + transp
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !   Calculate leaf-to-canopy sensible heat flux.  Always consider both sides of   !
         ! leaves.                                                                         !
         !---------------------------------------------------------------------------------!
         flux_area = effarea_heat * initp%lai(ico) ! + pi18 * initp%wai(ico)
         hflxlc    = flux_area    * initp%leaf_gbh(ico)                                    &
                   * (initp%leaf_temp(ico) - initp%can_temp)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Calculate interception by leaves by scaling the intercepted water by the    !
         ! LAI of each cohort.  If this causes excess of water/ice over the leaf surface,  !
         ! no problem, the water will shed at adjust_veg_properties.                       !
         !---------------------------------------------------------------------------------!
         wshed             = 0.d0
         qwshed            = 0.d0
         dwshed            = 0.d0
         leaf_intercepted  = intercepted_max  * initp%lai(ico) * taii
         leaf_qintercepted = qintercepted_max * initp%lai(ico) * taii
         throughfall       = 0.d0
         qthroughfall      = 0.d0
         dthroughfall      = 0.d0
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the water energy balance for this cohort.                              !
         !---------------------------------------------------------------------------------!
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
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    If the detailed output is tracked, then we save the fluxes for this cohort.  !
         !---------------------------------------------------------------------------------!
         if (print_detailed) then
            dinitp%cfx_hflxlc      (ico)  = hflxlc
            dinitp%cfx_qwflxlc     (ico)  = qwflxlc
            dinitp%cfx_qwshed      (ico)  = qwshed
            dinitp%cfx_qtransp     (ico)  = qtransp
            dinitp%cfx_qintercepted(ico)  = leaf_qintercepted
         end if
         !---------------------------------------------------------------------------------!

         !----- Add the contribution of this cohort to total heat and evapotranspiration. -!
         wflxlc_tot   = wflxlc_tot   + wflxlc
         qwflxlc_tot  = qwflxlc_tot  + qwflxlc
         hflxlc_tot   = hflxlc_tot   + hflxlc
         transp_tot   = transp_tot   + transp
         qtransp_tot  = qtransp_tot  + qtransp

         !---------------------------------------------------------------------------------!
         !     Here we update the liquid/frozen water fluxes and their associated vari-    !
         ! ables:                                                                          !
         !                                                                                 !
         ! - wshed      : Water falling from vegetated canopy to soil surface;             !
         ! - intercepted: Precipitation that is intercepted by the vegetation;             !
         ! - throughfall: Precipitation that is never intercepted by the vegetation.       !
         !---------------------------------------------------------------------------------!
         wshed_tot        = wshed_tot        + wshed 
         qwshed_tot       = qwshed_tot       + qwshed
         dwshed_tot       = dwshed_tot       + dwshed
         intercepted_tot  = intercepted_tot  + leaf_intercepted
         qintercepted_tot = qintercepted_tot + leaf_qintercepted
         throughfall_tot  = throughfall_tot  + throughfall
         qthroughfall_tot = qthroughfall_tot + qthroughfall
         dthroughfall_tot = dthroughfall_tot + dthroughfall
      else
         !---------------------------------------------------------------------------------! 
         !     If there is not enough leaf biomass to safely solve the leaf energy and     !
         ! water balances, set leaf fluxes and interception to zero.                       !
         !---------------------------------------------------------------------------------!
         dinitp%leaf_energy(ico) = 0.d0
         dinitp%leaf_water(ico)  = 0.d0
         dinitp%psi_open(ico)    = 0.d0
         dinitp%psi_closed(ico)  = 0.d0

         !---------------------------------------------------------------------------------!
         !    If the detailed output is tracked, then we save the fluxes for this cohort.  !
         !---------------------------------------------------------------------------------!
         if (print_detailed) then
            dinitp%cfx_hflxlc      (ico)  = 0.d0
            dinitp%cfx_qwflxlc     (ico)  = 0.d0
            dinitp%cfx_qwshed      (ico)  = 0.d0
            dinitp%cfx_qtransp     (ico)  = 0.d0
            dinitp%cfx_qintercepted(ico)  = 0.d0
         end if
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Allow the complete bypass of precipitation if there are very few leaves.    !
         ! Add this tiny amount to the throughfall.                                        !
         !---------------------------------------------------------------------------------!
         throughfall_tot   = throughfall_tot  + intercepted_max  * initp%lai(ico) * taii
         qthroughfall_tot  = qthroughfall_tot + qintercepted_max * initp%lai(ico) * taii
         dthroughfall_tot  = dthroughfall_tot + dintercepted_max * initp%lai(ico) * taii
      end if
      !------------------------------------------------------------------------------------!






      !------------------------------------------------------------------------------------!
      ! WOOD BUDGET - Check whether the wood of this cohort hasn't been flagged as non-    !
      !               resolvable, i.e., the user wants wood thermodynamics, it belongs to  !
      !               a patch that is not too sparse, an it is not buried in snow.  We     !
      !               should compute energy and water at the cohort level only if the      !
      !               cohort is "safe".  Otherwise, we will set the wood energy and mass   !
      !               derivatives to zero, and convert any intercepted water into          !
      !               throughfall.  Later on, these "unsafe" cohorts will have their wood  !
      !               energy set to equilibrium with the canopy air space (temperature).   !
      !------------------------------------------------------------------------------------!
      if (initp%wood_resolvable(ico)) then

         !------ Define some aliases to indices -------------------------------------------!
         ipft  = cpatch%pft(ico)

         !---------------------------------------------------------------------------------!
         !     Define the minimum wood water to be considered, and the maximum amount      !
         ! possible.                                                                       !
         !---------------------------------------------------------------------------------!
         min_wood_water = rk4leaf_drywhc * initp%wai(ico)
         max_wood_water = rk4leaf_maxwhc * initp%wai(ico)

         !------ Calculate fraction of wood covered with water. ---------------------------!
         if (initp%wood_water(ico) > min_wood_water) then
            sigmaw = min(1.d0, (initp%wood_water(ico)/max_wood_water)**twothirds8)
         else
            sigmaw = 0.d0
         end if

         !---------------------------------------------------------------------------------!
         !    Here we must compute two different areas.  For transpiration, we want the    !
         ! leaf area index only, because we assume transpiration to happen only through    !
         ! leaves.  Evaporation of water/ice settled over the vegetation surface, or dew/  !
         ! frost formation must account the branches and stems as well.                    !
         !---------------------------------------------------------------------------------!
         !----- Find the wood specific humidity. ------------------------------------------!
         wood_shv   = rslif8(initp%can_prss,initp%wood_temp(ico))
         wood_shv   = wood_shv / (1.d0 + wood_shv)
         !----- Evaporation/condensation "flux" -------------------------------------------!
         wflxwc_try = initp%wai(ico) * initp%wood_gbw(ico) * (wood_shv - initp%can_shv)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Compute the evapotranspiration or dew/frost deposition.                      !
         !---------------------------------------------------------------------------------!
         if (wflxwc_try >= 0.d0) then  
            !------------------------------------------------------------------------------!
            !    Probably evapotranspiration, as long as the canopy air is not saturated   !
            ! or the user doesn't mind that super-saturation occur.                        !
            !------------------------------------------------------------------------------!
            if (supersat_ok .or. initp%can_rhv < 1.d0) then
               !---------------------------------------------------------------------------!
               !     Evaporation, energy is scaled by liquid/ice partition (no phase       !
               ! bias).  We scale by the relative area of wood that is actually covered    !
               ! with water.                                                               !
               !---------------------------------------------------------------------------!
               wflxwc  = wflxwc_try * sigmaw
               qwflxwc = wflxwc * (alvi8 - initp%wood_fliq(ico) * alli8)
               !---------------------------------------------------------------------------!

            else
               !----- Canopy is already saturated, no evapotranspiration is allowed. ------!
               wflxwc                 = 0.d0
               qwflxwc                = 0.d0
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!

         else
            !------------------------------------------------------------------------------!
            !     Dew/frost formation. The deposition will conserve the liquid/ice         !
            ! partition (or use the default if there is no water).                         !
            !------------------------------------------------------------------------------!
            wflxwc                 = wflxwc_try
            qwflxwc                = wflxwc * (alvi8 - initp%wood_fliq(ico) * alli8)
            !------------------------------------------------------------------------------!

         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !   Calculate wood-to-canopy sensible heat flux.  Consider the full circumference !
         ! of brances.                                                                     !
         !---------------------------------------------------------------------------------!
         flux_area = pi18 * initp%wai(ico)
         hflxwc    = flux_area * initp%wood_gbh(ico)                                       &
                   * (initp%wood_temp(ico) - initp%can_temp)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Calculate interception by wood by scaling the intercepted water by the      !
         ! WAI of each cohort.  If this causes excess of water/ice over the wood surface,  !
         ! no problem, the water will shed at adjust_veg_properties.                       !
         !---------------------------------------------------------------------------------!
         wshed             = 0.d0
         qwshed            = 0.d0
         dwshed            = 0.d0
         wood_intercepted  = intercepted_max  * initp%wai(ico) * taii
         wood_qintercepted = qintercepted_max * initp%wai(ico) * taii
         throughfall       = 0.d0
         qthroughfall      = 0.d0
         dthroughfall      = 0.d0
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the water energy balance for this cohort.                              !
         !---------------------------------------------------------------------------------!
         dinitp%wood_water(ico)  = wood_intercepted     & ! Intercepted water 
                                 - wflxwc               & ! Evaporation
                                 - wshed                ! ! Water shedding
         dinitp%wood_energy(ico) = initp%rshort_w(ico)  & ! Absorbed SW radiation
                                 + initp%rlong_w(ico)   & ! Net thermal radiation
                                 - hflxwc               & ! Sensible heat flux
                                 - qwflxwc              & ! Evaporation
                                 - qwshed               & ! Water shedding
                                 + wood_qintercepted    ! ! Intercepted water energy
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    If the detailed output is tracked, then we save the fluxes for this cohort.  !
         !---------------------------------------------------------------------------------!
         if (print_detailed) then
            dinitp%cfx_hflxwc      (ico)  = hflxwc
            dinitp%cfx_qwflxwc     (ico)  = qwflxwc
            dinitp%cfx_qwshed      (ico)  = dinitp%cfx_qwshed(ico) + qwshed
            dinitp%cfx_qintercepted(ico)  = dinitp%cfx_qintercepted(ico)                   &
                                          + wood_qintercepted
         end if
         !---------------------------------------------------------------------------------!

         !----- Add the contribution of this cohort to total heat and evapotranspiration. -!
         wflxwc_tot   = wflxwc_tot   + wflxwc
         qwflxwc_tot  = qwflxwc_tot  + qwflxwc
         hflxwc_tot   = hflxwc_tot   + hflxwc

         !---------------------------------------------------------------------------------!
         !     Here we update the liquid/frozen water fluxes and their associated vari-    !
         ! ables:                                                                          !
         !                                                                                 !
         ! - wshed      : Water falling from vegetated canopy to soil surface;             !
         ! - intercepted: Precipitation that is intercepted by the vegetation;             !
         ! - throughfall: Precipitation that is never intercepted by the vegetation.       !
         !---------------------------------------------------------------------------------!
         wshed_tot        = wshed_tot        + wshed 
         qwshed_tot       = qwshed_tot       + qwshed
         dwshed_tot       = dwshed_tot       + dwshed
         intercepted_tot  = intercepted_tot  + wood_intercepted
         qintercepted_tot = qintercepted_tot + wood_qintercepted
         throughfall_tot  = throughfall_tot  + throughfall
         qthroughfall_tot = qthroughfall_tot + qthroughfall
         dthroughfall_tot = dthroughfall_tot + dthroughfall
      else
         !---------------------------------------------------------------------------------! 
         !     If there is not enough leaf biomass to safely solve the leaf energy and     !
         ! water balances, set leaf fluxes and interception to zero.                       !
         !---------------------------------------------------------------------------------!
         dinitp%wood_energy(ico) = 0.d0
         dinitp%wood_water(ico)  = 0.d0
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    If the detailed output is tracked, then we save the fluxes for this cohort.  !
         !---------------------------------------------------------------------------------!
         if (print_detailed) then
            dinitp%cfx_hflxwc      (ico)  = 0.d0
            dinitp%cfx_qwflxwc     (ico)  = 0.d0
         end if
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Allow the complete bypass of precipitation if there are very few leaves.    !
         ! Add this tiny amount to the throughfall.                                        !
         !---------------------------------------------------------------------------------!
         throughfall_tot   = throughfall_tot  + intercepted_max  * initp%wai(ico) * taii
         qthroughfall_tot  = qthroughfall_tot + qintercepted_max * initp%wai(ico) * taii
         dthroughfall_tot  = dthroughfall_tot + dintercepted_max * initp%wai(ico) * taii
      end if
      !------------------------------------------------------------------------------------!
   end do cohortloop
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !     Update the log of potential temperature (entropy), water vapour specific mass,    !
   ! and CO2 mixing ratio of the canopy air space.                                         !
   ! wcapcan: can_rhos * can_depth                  (water capacity)                       !
   ! hcapcan: can_rhos * can_depth * cp8 * can_temp (entropy capacity times temperature)   !
   ! ccapcan: can_rhos * can_depth * mmdryi         (carbon capacity)                      !
   !---------------------------------------------------------------------------------------!
   dinitp%can_lntheta = ( hflxgc + hflxlc_tot + hflxwc_tot + hflxac) * hcapcani
   dinitp%can_shv     = ( wflxgc - dewgndflx + wflxlc_tot                                  &
                        + wflxwc_tot + transp_tot +  wflxac)         * wcapcani
   dinitp%can_co2     = ( cflxgc + cflxlc_tot + cflxac)              * ccapcani
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Integrate diagnostic variables - These are not activated unless fast file-type    !
   ! outputs are selected. This will speed up the integrator.                              !
   !---------------------------------------------------------------------------------------!
   if (fast_diagnostics .or. print_detailed) then


      dinitp%avg_carbon_ac    = cflxac                      ! Carbon flx,  Atmo->Canopy
      dinitp%avg_sensible_ac  = hflxac                      ! Sens. heat,  Atmo->Canopy
      dinitp%avg_vapor_ac     = wflxac                      ! Lat.  heat,  Atmo->Canopy

      dinitp%avg_sensible_lc  = hflxlc_tot                  ! Sens. heat,  Leaf->Canopy
      dinitp%avg_vapor_lc     = wflxlc_tot                  ! Lat.  heat,  Leaf->Canopy

      dinitp%avg_sensible_wc  = hflxwc_tot                  ! Sens. heat,  Wood->Canopy
      dinitp%avg_vapor_wc     = wflxwc_tot                  ! Lat.  heat,  Wood->Canopy

      dinitp%avg_sensible_gc  = hflxgc                      ! Sens. heat,  Grnd->Canopy
      dinitp%avg_transp       = transp_tot                  ! Transpiration
      dinitp%avg_vapor_gc     = wflxgc                      ! Lat.  heat,  Grnd->Canopy
      dinitp%avg_dew_cg       = dewgndflx                   ! Lat.  heat,  Canopy->Grnd

      !----- Total evaporation to the canopy air space. -----------------------------------!
      dinitp%avg_evap         = wflxgc - dewgndflx + wflxlc_tot + wflxwc_tot


      dinitp%avg_wshed_vg     = wshed_tot                   ! Water shedding,Leaf->Grnd
      dinitp%avg_qwshed_vg    = qwshed_tot                  ! Water shedding,Leaf->Grnd
      dinitp%avg_intercepted  = intercepted_tot             ! Intercepted,   Atmo->Leaf
      dinitp%avg_qintercepted = qintercepted_tot            ! Intercepted,   Atmo->Lead
      dinitp%avg_throughfall  = throughfall_tot             ! Throughfall,   Atmo->Grnd
      dinitp%avg_qthroughfall = qthroughfall_tot            ! Throughfall,   Atmo->Grnd

   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Update the absorbed shortwave radiation and the net radiation.                   !
   !---------------------------------------------------------------------------------------!
   if (fast_diagnostics .or. checkbudget .or. print_detailed) then
      !----- Average shortwave radiation. -------------------------------------------------!
      dinitp%avg_rshort_gnd    = dble(csite%rshort_g(ipa))
      do k=1,initp%nlev_sfcwater
         dinitp%avg_rshort_gnd = dinitp%avg_rshort_gnd + dble(csite%rshort_s(k,ipa))
      end do
      !------------------------------------------------------------------------------------!

      !----- Average longwave radiation. --------------------------------------------------!
      dinitp%avg_rlong_gnd = dble(csite%rlong_g(ipa)) + dble(csite%rlong_s(ipa))
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Update the budget variables.                                                      !
   !---------------------------------------------------------------------------------------!
   if (checkbudget) then
      dinitp%co2budget_loss2atm = - cflxac
      dinitp%ebudget_loss2atm   = - eflxac
      dinitp%wbudget_loss2atm   = - wflxac
      dinitp%co2budget_storage  = dinitp%co2budget_storage + cflxgc + cflxlc_tot + cflxac
      dinitp%ebudget_storage    = dinitp%ebudget_storage   + eflxac                        &
                                + dinitp%avg_rshort_gnd    + dinitp%avg_rlong_gnd
      dinitp%wbudget_storage    = dinitp%wbudget_storage   + wflxac
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These variables below are virtual copies of the variables above, but are here for !
   ! for use in the coupled model. They form the set of canopy-atmosphere fluxes that are  !
   ! used for turbulent closure. These variables are also zeroed and normalized every      !
   ! dtlsm timestep, the others are likely averaged over the analysis period.              ! 
   !---------------------------------------------------------------------------------------!
   dinitp%upwp = -(initp%ustar*initp%ustar)
   dinitp%qpwp = -(initp%ustar*initp%qstar)
   dinitp%cpwp = -(initp%ustar*initp%cstar)
   dinitp%tpwp = -(initp%ustar*initp%tstar)
   dinitp%wpwp = vertical_vel_flux8(initp%zeta,initp%tstar,initp%ustar)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     If the single pond layer is too thin, we force equilibrium with top soil layer.   !
   ! This is done in two steps: first, we transfer all the energy to the top soil layer.   !
   ! Then, in adjust_sfcw_properties, we make them in thermal equilibrium.  We use the     !
   ! flag_sfcwater to transfer the energy rather than the actual mass because if a layer   !
   ! is considered to thin at the beginning of the time step, we want it to keep the same  !
   ! status throughout the entire step.                                                    !
   !---------------------------------------------------------------------------------------!
   select case (initp%flag_sfcwater)
   case (1)
      dinitp%soil_energy(mzg)   = dinitp%soil_energy(mzg)                                  &
                                + dinitp%sfcwater_energy(1) * dslzi8(mzg)
      dinitp%sfcwater_energy(1) = 0.d0
   end select
   !---------------------------------------------------------------------------------------!

   return
end subroutine canopy_derivs_two
!==========================================================================================!
!==========================================================================================!

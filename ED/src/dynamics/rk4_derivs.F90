!==========================================================================================!
!==========================================================================================!
! Subroutine leaf_derivs_ar                                                                !
!                                                                                          !
!     This subroutine finds the fast-scale derivatives at canopy, soil, and leaf surface.  !
! This subroutine is based on LEAF-3, except that here only the derivative is computed,    !
! whereas in LEAF-3 the actual step is done at once. This derivative will be used for the  !
! Runge-Kutta integration step.                                                            !
!------------------------------------------------------------------------------------------!
subroutine leaf_derivs_ar(initp,dinitp,csite,ipa,isi,ipy,rhos,prss,pcpg,qpcpg,dpcpg        &
                         ,atm_tmp,exner,geoht,vels,atm_shv,atm_co2,lsl)
  
   use ed_state_vars , only : sitetype     & ! structure
                            , rk4patchtype ! ! structure
   use consts_coms   , only : cp           ! ! intent(in)
   use grid_coms     , only : nzg          ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: initp   ! Structure with RK4 intermediate state
   type(rk4patchtype) , target     :: dinitp  ! Structure with RK4 derivatives
   type(sitetype)     , target     :: csite   ! This site (with previous values);
   integer            , intent(in) :: ipa     ! Patch ID
   integer            , intent(in) :: isi     ! Site ID
   integer            , intent(in) :: ipy     ! Polygon ID
   integer            , intent(in) :: lsl     ! Lowest soil level
   real               , intent(in) :: rhos    ! Air density
   real               , intent(in) :: prss    ! Air pressure
   real               , intent(in) :: pcpg    ! Precipitation rate
   real               , intent(in) :: qpcpg   ! Precipitation heat
   real               , intent(in) :: dpcpg   ! Precipitation heat
   real               , intent(in) :: atm_tmp ! Air temperature
   real               , intent(in) :: exner   ! Air Exner function
   real               , intent(in) :: geoht   ! Geopotential height
   real               , intent(in) :: vels    ! Wind Velocity
   real               , intent(in) :: atm_shv ! Air vapour mixing ratio
   real               , intent(in) :: atm_co2 ! Air CO2 mixing ratio
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   ! Depending on the type of compilation, interfaces must be explicitly declared.         !
   !---------------------------------------------------------------------------------------!
#if USE_INTERF
   interface
      !------------------------------------------------------------------------------------!
      !    Subroutine that computes the characteristic scales.                             !
      !------------------------------------------------------------------------------------!
      subroutine ed_stars(tha,rva,chia,thv,zpm,um,rough,ustar,rstar,tstar,cstar,can_shv    &
                         ,can_co2)
         implicit none
         real, intent(in) :: rough
         real, intent(out) :: ustar
         real, intent(out) :: rstar
         real, intent(out) :: tstar
         real, intent(out) :: cstar
         real, intent(in) :: can_shv
         real, intent(in) :: can_co2
         real ::tha,rva,chia,thv,zpm,um
      end subroutine ed_stars
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Subroutine that computes the canopy and leaf fluxes.                            ! 
      !------------------------------------------------------------------------------------!
      subroutine leaftw_derivs_ar(initp,dinitp,csite,ipa,isi,ipy,rhos,prss,pcpg,qpcpg      &
                                 ,dpcpg,atm_tmp,exner,geoht,lsl)
         use ed_state_vars , only : rk4patchtype & ! structure
                                  , sitetype     ! ! structure
         implicit none
         type(rk4patchtype) , target     :: initp  
         type(rk4patchtype) , target     :: dinitp 
         type(sitetype)     , target     :: csite
         integer            , intent(in) :: ipa,isi,ipy
         integer            , intent(in) :: lsl
         real               , intent(in) :: rhos
         real               , intent(in) :: prss
         real               , intent(in) :: pcpg
         real               , intent(in) :: qpcpg
         real               , intent(in) :: dpcpg
         real               , intent(in) :: atm_tmp
         real               , intent(in) :: exner
         real               , intent(in) :: geoht
      end subroutine leaftw_derivs_ar
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!

   !----- Resetting the energy budget. ----------------------------------------------------!
   dinitp%ebudget_latent = 0.0

   !----- Compute friction velocities. ----------------------------------------------------!
   call ed_stars(atm_tmp,atm_shv,atm_co2,initp%can_temp,geoht,vels,initp%rough,initp%ustar &
                ,initp%rstar,initp%tstar,initp%cstar,initp%can_shv,initp%can_co2)

   initp%tstar = initp%tstar * cp / exner

   call leaftw_derivs_ar(initp,dinitp,csite,ipa,isi,ipy,rhos,prss,pcpg,qpcpg,dpcpg,atm_tmp &
                        ,exner,geoht,lsl)

   !----- Nlev_sfcwater derivative... I doubt it's really used... -------------------------!
   dinitp%nlev_sfcwater = initp%nlev_sfcwater

   return
end subroutine leaf_derivs_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine leaftw_derivs_ar(initp,dinitp,csite,ipa,isi,ipy,rhos,prss,pcpg,qpcpg,dpcpg      &
                           ,atm_tmp,exner,geoht,lsl)
   use max_dims             , only : nzgmax               & ! intent(in)
                                   , nzsmax               ! ! intent(in)
   use consts_coms          , only : alvl                 & ! intent(in)
                                   , cliqvlme             & ! intent(in)
                                   , cpi                  & ! intent(in)
                                   , alvi                 & ! intent(in)
                                   , t3ple                & ! intent(in)
                                   , cliq                 & ! intent(in)
                                   , cice                 & ! intent(in)
                                   , tsupercool           & ! intent(in)
                                   , wdns                 & ! intent(in)
                                   , wdnsi                ! ! intent(in)
   use grid_coms            , only : nzg                  & ! intent(in)
                                   , nzs                  ! ! intent(in)
   use soil_coms            , only : soil                 & ! intent(in)
                                   , slz                  & ! intent(in)
                                   , dslz                 & ! intent(in)
                                   , dslzi                & ! intent(in)
                                   , water_stab_thresh    & ! intent(in)
                                   , infiltration_method  & ! intent(in)
                                   , dslzti               & ! intent(in)
                                   , slcons1              & ! intent(in)
                                   , slzt                 & ! intent(in)
                                   , min_sfcwater_mass    & ! intent(in)
                                   , ss                   & ! intent(in)
                                   , isoilbc              ! ! intent(in)
   use misc_coms            , only : dtlsm                & ! intent(in)
                                   , current_time         ! ! intent(in)
   use canopy_radiation_coms, only : lai_min              ! ! intent(in)
   use ed_state_vars        , only : sitetype             & ! structure
                                   , patchtype            & ! structure
                                   , rk4patchtype         ! ! structure
   use rk4_coms             , only : rk4eps               ! ! intent(in)
   use ed_therm_lib         , only : ed_grndvap           ! ! subroutine
   use therm_lib            , only : qtk                  & ! subroutine
                                   , qwtk                 & ! subroutine
                                   , qwtk8                ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)  , target     :: initp         ! RK4 structure, intermediate step
   type(rk4patchtype)  , target     :: dinitp        ! RK4 structure, derivatives
   type(sitetype)      , target     :: csite         ! Current site (before integration)
   integer                          :: ipa           ! Current patch ID
   integer                          :: isi           ! Current site ID
   integer                          :: ipy           ! Current polygon ID
   integer             , intent(in) :: lsl           ! Lowest soil level
   real                , intent(in) :: pcpg          ! Precipitation rate
   real                , intent(in) :: qpcpg         ! Precipitation heat
   real                , intent(in) :: dpcpg         ! Precipitation density
   real                , intent(in) :: exner         ! Exner function
   real                , intent(in) :: atm_tmp       ! Air temperature
   real                , intent(in) :: rhos          ! Air density
   real                , intent(in) :: prss          ! Pressure
   real                , intent(in) :: geoht         ! Geopotential height
   !----- Local variables -----------------------------------------------------------------!
   integer                          :: k, k1, k2     ! Level counters
   integer                          :: ksn           ! # of temporary water/snow layers
   integer                          :: nsoil         ! Short for csite%soil_text(k,ipa)
   real                             :: wgpfrac       ! Fractional soil moisture
   real                             :: soilcond      ! Soil conductivity
   real                             :: snden         ! Snow/water density
   real                             :: hflxgc        ! Ground -> canopy heat flux
   real                             :: wflxgc        ! Ground -> canopy water flux
   real                             :: qwflxgc       ! Ground -> canopy latent heat flux
   real                             :: dewgnd        ! Dew/frost flux to ground
   real                             :: qdewgnd       ! Dew/frost heat flux to ground
   real                             :: ddewgnd       ! Dew/frost density flux to ground
   real                             :: wshed         ! Water shedding flux
   real                             :: qwshed        ! Energy flux due to water shedding
   real                             :: dwshed        ! Density flux due to water shedding
   real, dimension(nzgmax+nzsmax)   :: rfactor       ! 
   real, dimension(nzgmax+nzsmax+1) :: hfluxgsc      ! Surface -> canopy heat flux
   real, dimension(nzg+nzs+1)       :: w_flux        ! Water flux (aux. variable)
   real, dimension(nzg+nzs+1)       :: qw_flux       ! Heat flux (aux. variable)
   real, dimension(nzs+1)           :: d_flux        ! Density flux
   real                             :: wgpmid        ! Soil in between layers
   real                             :: wloss         ! Water loss due to transpiration
   real                             :: qwloss        ! Energy loss due to transpiration
   real                             :: dqwt          ! Energy adjustment aux. variable
   real                             :: fracliq       ! Fraction of liquid water
   real                             :: tempk         ! Temperature
   real                             :: qwgoal        ! Goal energy for thin snow layers.
   real(kind=8)                     :: wprevious     ! Previous water content
   real                             :: infilt        ! Surface infiltration rate
   real                             :: qinfilt       ! Surface infiltration heat rate
   real                             :: snowdens      ! Snow density (kg/m2)
   real                             :: soilhcap      ! Soil heat capacity
   real                             :: surface_water ! Temp. variable. Availible liquid 
                                                     !   water on the soil surface (kg/m2)
   real                             :: freezeCor     ! Correction to conductivity for 
                                                     !    partially frozen soil.
   real                             :: surface_temp  ! Surface temperature.
   real                             :: surface_fliq  ! Liquid fraction at surface.
   real                             :: int_sfcw_u    ! Intensive sfc. water internal en.
   !----- Constants -----------------------------------------------------------------------!
   logical , parameter  :: debug = .false.  ! Debugging output flag (T/F)
   real    , parameter  :: freezeCoef = 7.0 ! Exponent in the frozen soil hydraulic 
                                            !    conductivity correction.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! Depending on the type of compilation, interfaces must be explicitly declared.         !
   !---------------------------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine canopy_derivs_two_ar(initp,dinitp,csite,ipa,isi,ipy,hflxgc,wflxgc,qwflxgc &
                                     ,dewgndflx,qdewgndflx,ddewgndflx,wshed_tot,qwshed_tot &
                                     ,dwshed_tot,rhos,prss,pcpg,qpcpg,dpcpg,exner,geoht    &
                                     ,atm_tmp,surface_temp,surface_fliq,lsl)
         use ed_state_vars, only: rk4patchtype  & ! structure
                                , sitetype      & ! structure
                                , patchtype     ! ! structure
         implicit none
         type (rk4patchtype) , target      :: initp, dinitp
         type (sitetype)     , target      :: csite
         integer             , intent(in)  :: ipa, isi, ipy, lsl
         real                , intent(in)  :: rhos, atm_tmp, prss, exner, geoht 
         real                , intent(in)  :: pcpg, qpcpg, dpcpg
         real                , intent(in)  :: surface_temp, surface_fliq
         real                , intent(out) :: hflxgc, wflxgc, qwflxgc
         real                , intent(out) :: dewgndflx, qdewgndflx, ddewgndflx
         real                , intent(out) :: wshed_tot, qwshed_tot, dwshed_tot
      end subroutine canopy_derivs_two_ar
   end interface
#endif
   !---------------------------------------------------------------------------------------!

   !----- Copying the # of surface water/snow layers to a shortcut ------------------------!
   ksn = initp%nlev_sfcwater
  
   !---- Initializing some vertical integration variables --------------------------------!
   w_flux  = 0.0
   qw_flux = 0.0
   d_flux  = 0.0

   !----- Initialize derivatives to zero -------------------------------------------------!
   dinitp%soil_energy(:)     = 0.0
   dinitp%soil_water(:)      = 0.0d+0
   dinitp%sfcwater_depth(:)  = 0.0
   dinitp%sfcwater_energy(:) = 0.0
   dinitp%sfcwater_mass(:)   = 0.0
   dinitp%virtual_heat       = 0.0
   dinitp%virtual_water      = 0.0
   dinitp%virtual_depth      = 0.0
   initp%extracted_water(:)  = 0.0
   dinitp%avg_smoist_gc(:)   = 0.0

   !----- Copying the soil texture flag to a shortcut -------------------------------------!
   nsoil = csite%ntext_soil(nzg,ipa)
   k = max(1,ksn)
   if (abs(initp%sfcwater_mass(k)) > min_sfcwater_mass) then
      int_sfcw_u = initp%sfcwater_energy(k)/initp%sfcwater_mass(k)
   else
      int_sfcw_u = 0.
   end if
   
   call ed_grndvap(ksn,nsoil,initp%soil_water(nzg),initp%soil_energy(nzg),int_sfcw_u,rhos  &
                  ,initp%can_shv,initp%ground_shv,initp%surface_ssh,surface_temp           &
                  ,surface_fliq)

   !---------------------------------------------------------------------------------------!
   !     Calculate water available to vegetation (in meters). SLZ is specified in RAMSIN.  !
   ! Each element of the array sets the value of the bottom of a corresponding soil layer. !
   ! Eg, SLZ = -2, -1, -0.5, -0.25.  There are four soil layers in this example; soil      !
   ! layer 1 goes from 2 meters below the surface to 1 meter below the surface.            !
   !---------------------------------------------------------------------------------------!
   initp%available_liquid_water(nzg) = dslz(nzg)                                           &
        * max(0.0,initp%soil_fracliq(nzg)*(sngl(initp%soil_water(nzg))-soil(nsoil)%soilcp))

   do k = nzg - 1, lsl, -1
      nsoil = csite%ntext_soil(k,ipa)
      initp%available_liquid_water(k) = initp%available_liquid_water(k+1) + dslz(k)        &
           * max(0.0,(sngl(initp%soil_water(k))-soil(nsoil)%soilcp)*initp%soil_fracliq(k))
   end do

   !---------------------------------------------------------------------------------------!
   !     Compute gravitational potential plus moisture potential, psi + z (psiplusz) [m],  !
   ! liquid water content (soil_liq) [m], and 99% the remaining water capacity (soilair99) !
   ! [m].                                                                                  !
   !---------------------------------------------------------------------------------------!
   do k = lsl, nzg
      nsoil = csite%ntext_soil(k,ipa)
      initp%psiplusz(k) = slzt(k) + soil(nsoil)%slpots                                     &
                        * (soil(nsoil)%slmsts / sngl(initp%soil_water(k)))                 & 
                        ** soil(nsoil)%slbs

      !------------------------------------------------------------------------------------!
      !    Soil liquid water must be converted to meters of liquid water per layer. This   !
      ! requires multiplication of volumetric water content, m3(water)/m3 must be          !
      ! multiplied by depth to get a depth of water.                                       !
      !------------------------------------------------------------------------------------!
      initp%soil_liq(k) = max(0.0                                                          &
                             ,( sngl(initp%soil_water(k)) - soil(nsoil)%soilcp)            &
                              * initp%soil_fracliq(k) )
      initp%soilair99(k) = (1.-rk4eps) * soil(nsoil)%slmsts - sngl(initp%soil_water(k))
      initp%soilair01(k) = (sngl(initp%soil_water(k)) - (1.+rk4eps) * soil(nsoil)%soilcp)  &
                         * initp%soil_fracliq(k)
   end do
   !---------------------------------------------------------------------------------------!

 
   !----- Get derivatives of canopy variables. --------------------------------------------!
   call canopy_derivs_two_ar(initp,dinitp,csite,ipa,isi,ipy,hflxgc,wflxgc,qwflxgc,dewgnd   &
                            ,qdewgnd,ddewgnd,wshed,qwshed,dwshed,rhos,prss,pcpg,qpcpg      &
                            ,dpcpg,exner,geoht,atm_tmp,surface_temp,surface_fliq,lsl)


   !---------------------------------------------------------------------------------------!
   !     Here we check whether it is bedrock or not. I think the reason for this check is  !
   ! to avoid problems with soil_water never being defined for bedrocks, otherwise simply  !
   ! assuming soilcond1 and soilcond2 would suffice.                                       !
   !---------------------------------------------------------------------------------------!
   do k = lsl, nzg
      nsoil = csite%ntext_soil(k,ipa)
      if(nsoil /= 13)then
         wgpfrac = min(sngl(initp%soil_water(k)) / soil(nsoil)%slmsts,1.0)
         soilcond = soil(nsoil)%soilcond0                                                  &
                  + wgpfrac * (soil(nsoil)%soilcond1 + wgpfrac *  soil(nsoil)%soilcond2)
      else
         soilcond=soil(nsoil)%soilcond0
      end if
      rfactor(k) = dslz(k) / soilcond
   end do


   !----- Snow/surface water --------------------------------------------------------------!
   do k = 1, ksn
      if(initp%sfcwater_depth(k) > 0.0)then
         snden = initp%sfcwater_mass(k) / initp%sfcwater_depth(k)
         rfactor(k+nzg) = initp%sfcwater_depth(k)                                          &
                        / (ss(1) * exp(ss(2) * initp%sfcwater_tempk(k))                    &
                        * (ss(3) + snden * (ss(4) + snden * (ss(5) + snden * ss(6)))))
      else
         rfactor(k+nzg) = 0.0
      end if
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Calculate the sensible heat fluxes find soil and sfcwater internal sensible heat  !
   ! fluxes (hfluxgsc) [W/m2].                                                             !
   !---------------------------------------------------------------------------------------!
   hfluxgsc(:) = 0.0
   do k = lsl+1, nzg
      hfluxgsc(k) = - (initp%soil_tempk(k) - initp%soil_tempk(k-1))                        &
                    / ((rfactor(k) + rfactor(k-1)) * .5)    
      dinitp%avg_sensible_gg(k-1) = hfluxgsc(k)  ! Diagnostic
   end do

   !----- If temporary water/snow layers exist, compute them now... -----------------------!
   if (ksn >= 1) then
      hfluxgsc(nzg+1) = - (initp%sfcwater_tempk(1) - initp%soil_tempk(nzg))                &
                        / ((rfactor(nzg+1)   + rfactor(nzg)) * .5)

      do k = 2,ksn
         hfluxgsc(nzg+k) = - (initp%sfcwater_tempk(k) - initp%sfcwater_tempk(k-1))         &
                         /   ((rfactor(nzg+k) + rfactor(nzg+k-1)) * .5)
      end do
   end if

   !----- Heat flux (hfluxgsc) at soil or sfcwater top from longwave, sensible [W/m^2] ----!
   hfluxgsc(nzg+ksn+1) = hflxgc + qwflxgc - csite%rlong_g(ipa) - csite%rlong_s(ipa)

   !----- Heat flux -----------------------------------------------------------------------!
   dinitp%avg_sensible_gg(nzg) = hfluxgsc(nzg+ksn+1) ! Diagnostic

   !---------------------------------------------------------------------------------------!
   !    Update soil U values [J/m³] from sensible heat, upward water vapor (latent heat)   !
   ! and longwave fluxes. This excludes effects of dew/frost formation, precipitation,     !
   ! shedding, and percolation.                                                            !
   !---------------------------------------------------------------------------------------!
   do k = lsl,nzg
      dinitp%soil_energy(k) = dslzi(k) * (hfluxgsc(k)- hfluxgsc(k+1))
   end do

   !----- Update soil Q values [J/m³] from shortwave flux. --------------------------------!
   dinitp%soil_energy(nzg) = dinitp%soil_energy(nzg) + dslzi(nzg) * csite%rshort_g(ipa)


   !---------------------------------------------------------------------------------------!
   !    Update surface water U values [J/m²] from sensible heat, upward water vapor        !
   ! (latent heat), longwave, and shortwave fluxes.  This excludes effects of dew/frost    !
   ! formation, precipitation, shedding and percolation.                                   !
   !---------------------------------------------------------------------------------------!
   do k = 1,ksn
     dinitp%sfcwater_energy(k) = hfluxgsc(k+nzg) - hfluxgsc(k+1+nzg)                       &
                               +  csite%rshort_s(k,ipa)
   end do

   !---------------------------------------------------------------------------------------!
   !     Calculate the fluxes of water with their associated heat fluxes. Update top soil  !
   ! or snow moisture from evaporation only.                                               !
   !                                                                                       !
   !     New moisture, qw, and depth from dew/frost formation, precipitation, shedding,    !
   ! and percolation.  ksnnew is the layer that receives the new condensate that comes     !
   ! directly from the air above.  If there is no pre-existing snowcover, this is a        !
   ! temporary "snow" layer.  This weird factor is essentially the one used by default in  !
   ! RAMS 4.3.0.  I'd rather use a different factor, huh?                                  !
   !---------------------------------------------------------------------------------------!
   w_flux(nzg+ksn+1)  = -  dewgnd -  wshed
   qw_flux(nzg+ksn+1) = - qdewgnd - qwshed
   d_flux(ksn+1)      =   ddewgnd + dwshed

   dinitp%avg_vapor_gc  = wflxgc   ! Diagnostic
   dinitp%avg_dew_cg    = dewgnd   ! Diagnostic
   dinitp%avg_qwshed_vg = qwshed   ! Diagnostic
   dinitp%avg_wshed_vg  = wshed    ! Diagnostic


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
      dinitp%sfcwater_mass(ksn)   = -w_flux(nzg+ksn+1) - wflxgc
      dinitp%sfcwater_energy(ksn) = dinitp%sfcwater_energy(ksn) - qw_flux(nzg+ksn+1)
      dinitp%sfcwater_depth(ksn) = -d_flux(ksn+1)
   else if (w_flux(nzg+1) < 0.0) then
      dinitp%virtual_heat  = -qw_flux(nzg+1)
      dinitp%virtual_water = -w_flux(nzg+1)
      dinitp%virtual_depth = -d_flux(1)
      qw_flux(nzg+1)       = 0.0
      w_flux(nzg+1)        = wflxgc
   else
      w_flux(nzg+1)        = w_flux(nzg+1) + wflxgc
   end if
   !---------------------------------------------------------------------------------------!

   dinitp%avg_smoist_gg(nzg) = w_flux(nzg+ksn+1)  ! Diagnostic



   !---------------------------------------------------------------------------------------!
   !     Find amount of water transferred between soil layers (w_flux) [m] modulated by    !
   ! the liquid water fraction.                                                            !
   !---------------------------------------------------------------------------------------!
   w_flux(nzg+1) = w_flux(nzg+1) * wdnsi ! now in m/s

   !---------------------------------------------------------------------------------------!
   !     Alternate surface infiltration (MCD) based on surface conductivity not capacity.  !
   !---------------------------------------------------------------------------------------!
   if (infiltration_method /= 0) then
      call fatal_error ('Running alt infiltation when we shouldn''t be',&
                       &'leaftw_derivs_ar','rk4_derivs.F90')

      if (initp%virtual_water /= 0.0) then  !!process "virtural water" pool
         nsoil = csite%ntext_soil(nzg,ipa)
         if (nsoil /= 13) then
            call qtk(initp%virtual_heat/initp%virtual_water,tempk,fracliq)
            infilt = -dslzi(nzg)* 0.5 * slcons1(nzg,nsoil)                                 &
                   * (sngl(initp%soil_water(nzg)) / soil(nsoil)%slmsts)                    &
                     **(2. * soil(nsoil)%slbs + 3.)                                        &
                   * (initp%psiplusz(nzg)-sngl(initp%virtual_water)/2000.0) & !diff. in pot.
                   * .5 * (initp%soil_fracliq(nzg)+ fracliq)         ! mean liquid fraction
            qinfilt = infilt * cliqvlme * (tempk - tsupercool)
            !----- Adjust other rates accordingly -----------------------------------------!
            w_flux(nzg+1)  = w_flux(nzg+1) + infilt
            qw_flux(nzg+1) = qw_flux(nzg+1)+ qinfilt
            dinitp%virtual_water = dinitp%virtual_water - infilt*wdns
            dinitp%virtual_heat  = dinitp%virtual_heat  - qinfilt
         end if
      end if  !! end virtual water pool
      if (initp%nlev_sfcwater >= 1) then !----- Process "snow" water pool -----------------! 
         call qtk(initp%sfcwater_energy(1)/initp%sfcwater_mass(1),tempk,fracliq)
         surface_water = initp%sfcwater_mass(1)*fracliq*wdnsi !(m/m2)
         nsoil = csite%ntext_soil(nzg,ipa)
         if (nsoil /= 13) then
            !----- Calculate infiltration rate (m/s) --------------------------------------!
            infilt = -dslzi(nzg) * 0.5 * slcons1(nzg,nsoil)                                &
                   * (sngl(initp%soil_water(nzg)) / soil(nsoil)%slmsts)                    &
                     **(2. * soil(nsoil)%slbs + 3.)                                        &
                   * (initp%psiplusz(nzg) - surface_water/2.0) &  !difference in potentials
                   * .5 * (initp%soil_fracliq(nzg)+ fracliq)    ! mean liquid fraction
            qinfilt = infilt * cliqvlme * (tempk - tsupercool)
            !----- Adjust other rates accordingly -----------------------------------------!
            w_flux(nzg+1)             = w_flux(nzg+1) + infilt
            qw_flux(nzg+1)            = qw_flux(nzg+1)+ qinfilt 
            dinitp%sfcwater_mass(1)   = dinitp%sfcwater_mass(1) - infilt*wdns
            dinitp%sfcwater_energy(1) = dinitp%sfcwater_energy(1) - qinfilt
            dinitp%sfcwater_depth(1)  = dinitp%sfcwater_depth(1) - infilt
         end if
      end if  ! End snow water pool
   end if  !! End alternate infiltration
   !---------------------------------------------------------------------------------------!


   !----- Keep qw_flux in W/m2. -----------------------------------------------------------!
   do k = lsl+1, nzg
      nsoil = csite%ntext_soil(k,ipa)
      if(nsoil /= 13 .and. csite%ntext_soil(k-1,ipa) /= 13)then

         wgpmid    = 0.5 * sngl(initp%soil_water(k) + initp%soil_water(k-1))
         freezeCor = 0.5 * (initp%soil_fracliq(k)+ initp%soil_fracliq(k-1))
         if(freezeCor < 1.0) freezeCor = 10.0**(-freezeCoef*(1.0-freezeCor))
         w_flux(k) = dslzti(k) * slcons1(k,nsoil)                                          &
                   * (wgpmid / soil(nsoil)%slmsts)**(2. * soil(nsoil)%slbs + 3.)           &
                   * (initp%psiplusz(k-1) - initp%psiplusz(k)) * freezeCor

        !----------------------------------------------------------------------------------!
        !      Limit water transfers to prevent over-saturation and over-depletion.        !
        !----------------------------------------------------------------------------------!
        if (w_flux(k) > 0.) then
           if (initp%soilair01(k-1) <= 0.0 .or. initp%soilair99(k) <= 0.0 ) w_flux(k) = 0.
        else
           if (initp%soilair99(k-1) <= 0.0 .or. initp%soilair01(k) <= 0.0 ) w_flux(k) = 0.
        end if
      end if
      !----- Only liquid water is allowed to flow, find qw_flux (W/m2) accordingly --------!
      qw_flux(k) = w_flux(k) * cliqvlme * (initp%soil_tempk(k) - tsupercool)
      dinitp%avg_smoist_gg(k-1) = w_flux(k)*wdns   ! Diagnostic
   end do

   !----- Boundary condition at the lowest soil level -------------------------------------!
   nsoil = csite%ntext_soil(lsl,ipa)
   if (nsoil /= 13 .and. isoilbc == 1) then
      !----- Free drainage ----------------------------------------------------------------!
      wgpmid      = sngl(initp%soil_water(lsl))
      freezeCor   = initp%soil_fracliq(lsl)
      if(freezeCor < 1.0) freezeCor = 10.0**(-freezeCoef*(1.0-freezeCor))
      w_flux(lsl) = dslzti(lsl) * slcons1(lsl,nsoil)                                       &
                  * (wgpmid / soil(nsoil)%slmsts)**(2. * soil(nsoil)%slbs + 3.)            &
                  * freezeCor
      !------------------------------------------------------------------------------------!
      !      Limit water transfers to prevent over-saturation and over-depletion.          !
      !------------------------------------------------------------------------------------!
      if (w_flux(lsl) > 0.) then
         if (initp%soilair99(lsl) <= 0.0) w_flux(lsl) = 0.
      else
         if (initp%soilair01(lsl) <= 0.0) w_flux(lsl) = 0.
      end if
      !----- Only liquid water is allowed to flow, find qw_flux (W/m2) accordingly --------!
      qw_flux(lsl) = w_flux(lsl) * cliqvlme * (initp%soil_tempk(lsl) - tsupercool)
   else
      !----- Bedrock, no flux accross it. -------------------------------------------------!
      w_flux(lsl)  = 0.
      qw_flux(lsl) = 0.
   end if

   !----- Finally, update soil moisture (impose minimum value of soilcp) and soil energy. -!
   do k = lsl,nzg
      dinitp%soil_water(k)  = dinitp%soil_water(k)                                         &
                            - dble(dslzi(k) * ( w_flux(k+1) -  w_flux(k)) )
      dinitp%soil_energy(k) =  dinitp%soil_energy(k)                                       &
                            - dslzi(k) * ( qw_flux(k+1) - qw_flux(k) )
   end do

   !---- Update soil moisture and energy from transpiration/root uptake. ------------------!
   if (csite%lai(ipa) > lai_min) then
      do k1 = lsl, nzg    ! loop over extracted water
         do k2=k1,nzg
            if (csite%ntext_soil(k2,ipa) /= 13) then
               if (initp%available_liquid_water(k1) > 0.0) then
                  
                  wloss = wdnsi * initp%extracted_water(k1)                                &
                        * initp%soil_liq(k2) / initp%available_liquid_water(k1)
                  
                  dinitp%soil_water(k2) = dinitp%soil_water(k2) - dble(wloss)
                  
                  !----- Energy: only liquid water is lost through transpiration. ---------!
                  qwloss = wloss * cliqvlme * (initp%soil_tempk(k2) - tsupercool)
                  dinitp%soil_energy(k2)   = dinitp%soil_energy(k2)   - qwloss
                  dinitp%avg_smoist_gc(k2) = dinitp%avg_smoist_gc(k2) - wdns*wloss
                  dinitp%ebudget_latent    = dinitp%ebudget_latent    + qwloss
               end if
            end if
         end do
      end do
   end if

   return
end subroutine leaftw_derivs_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine canopy_derivs_two_ar(initp,dinitp,csite,ipa,isi,ipy,hflxgc,wflxgc,qwflxgc       &
                               ,dewgndflx,qdewgndflx,ddewgndflx,wshed_tot,qwshed_tot       &
                               ,dwshed_tot,rhos,prss,pcpg,qpcpg,dpcpg,exner,geoht,atm_tmp  &
                               ,surface_temp,surface_fliq,lsl)
   use ed_state_vars         , only : rk4patchtype      & ! Structure
                                    , sitetype          & ! Structure
                                    , patchtype         ! ! Structure
   use consts_coms           , only : alvl              & ! intent(in)
                                    , cliq              & ! intent(in)
                                    , cice              & ! intent(in)
                                    , tsupercool        & ! intent(in)
                                    , cp                & ! intent(in)
                                    , cpi               & ! intent(in)
                                    , twothirds         & ! intent(in)
                                    , day_sec           & ! intent(in)
                                    , grav              & ! intent(in)
                                    , alvi              & ! intent(in)
                                    , alli              & ! intent(in)
                                    , umol_2_kgC        & ! intent(in)
                                    , mmdry             & ! intent(in)
                                    , mmdryi            & ! intent(in)
                                    , t3ple             & ! intent(in)
                                    , wdns              & ! intent(in)
                                    , wdnsi             & ! intent(in)
                                    , idns              ! ! intent(in)
   use grid_coms             , only : nzg               ! ! intent(in)
   use soil_coms             , only : soil              & ! intent(in)
                                    , dslzi             & ! intent(in)
                                    , water_stab_thresh & ! intent(in)
                                    , dewmax            ! ! intent(in)
   use canopy_radiation_coms , only : lai_min           ! ! intent(in)
   use canopy_air_coms       , only : min_veg_lwater    & ! intent(in)
                                    , max_veg_lwater    ! ! intent(in)
   use therm_lib             , only : qwtk              & ! subroutine
                                    , rslif             ! ! function
   use misc_coms             , only : dtlsm             ! ! intent(in)
   use ed_misc_coms          , only : fast_diagnostics  ! ! intent(in)
   use allometry             , only : dbh2ca            ! ! function
   use pft_coms              , only : water_conductance & ! intent(in)
                                    , q                 & ! intent(in)
                                    , qsw               ! ! intent(in)
   use rk4_coms              , only : debug             & ! intent(in)
                                    , toocold           & ! intent(in)
                                    , toohot            & ! intent(in)
                                    , lai_to_cover      & ! intent(in)
                                    , evap_area_one     & ! intent(in)
                                    , evap_area_two     ! ! intent(in)
                                       
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)     , target      :: csite          ! Current site
   type(rk4patchtype) , target      :: initp          ! RK4 structure, state vars
   type(rk4patchtype) , target      :: dinitp         ! RK4 structure, derivatives
   integer            , intent(in)  :: ipy            ! Current polygon ID
   integer            , intent(in)  :: isi            ! Current site ID
   integer            , intent(in)  :: ipa            ! Current patch ID
   integer            , intent(in)  :: lsl            ! Lowest soil level
   real               , intent(in)  :: rhos           ! Air density
   real               , intent(in)  :: atm_tmp        ! Air temperature
   real               , intent(in)  :: prss           ! Air pressure
   real               , intent(in)  :: pcpg           ! Precipitation rate
   real               , intent(in)  :: qpcpg          ! Precipitation heat
   real               , intent(in)  :: dpcpg          ! Precipitation density
   real               , intent(in)  :: exner          ! Air Exner function
   real               , intent(in)  :: geoht          ! Geopotential
   real               , intent(in)  :: surface_temp   ! Free surface temperature
   real               , intent(in)  :: surface_fliq   ! Free surface liquid fraction
   real               , intent(out) :: hflxgc         ! Ground->canopy sensible heat flux
   real               , intent(out) :: wflxgc         ! Ground->canopy water flux
   real               , intent(out) :: qwflxgc        ! Ground->canopy latent heat flux
   real               , intent(out) :: wshed_tot      ! Water shedding rate
   real               , intent(out) :: qwshed_tot     ! Water shedding energy flux
   real               , intent(out) :: dwshed_tot     ! Water shedding density flux
   real               , intent(out) :: dewgndflx      ! Dew/frost water flux
   real               , intent(out) :: qdewgndflx     ! Dew/frost heat flux
   real               , intent(out) :: ddewgndflx     ! Dew/frost density
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)    , pointer     :: cpatch           ! Current patch
   integer                          :: ico              ! Current cohort ID
   integer                          :: ipft             ! Shortcut for PFT type
   integer                          :: kroot            ! Level in which the root bottom is
   real                             :: can_frac         ! total fractional canopy coverage
   real                             :: transp           ! Cohort transpiration
   real                             :: cflxac           ! Air->canopy carbon flux
   real                             :: wflxac           ! Air->canopy water flux
   real                             :: hflxac           ! Air->canopy heat flux
   real                             :: c2               ! Coefficient (????)
   real                             :: c3               ! Coefficient (????)
   real                             :: hflxvc           ! Leaf->canopy heat flux
   real                             :: rasgnd           ! 
   real                             :: rbi              ! 
   real                             :: rd               !
   real                             :: sigmaw           !
   real                             :: wflxvc           !
   real                             :: wshed            !
   real                             :: qwshed           !
   real                             :: dwshed           !
   real                             :: zoveg,zveg       !
   real                             :: wcapcan          !
   real                             :: wcapcani         !
   real                             :: hcapcani         !
   real                             :: cflxgc           !
   real                             :: laii             !
   real                             :: wflx             !
   real                             :: qwflx            !
   real                             :: hflxvc_tot       !
   real                             :: transp_tot       !
   real                             :: cflxvc_tot       !
   real                             :: wflxvc_tot       !
   real                             :: qwflxvc_tot      !
   real                             :: rho_ustar        !
   real                             :: rdi              !
   real                             :: gpp_tot          !
   real                             :: storage_decay    !
   real                             :: leaf_flux        !
   real                             :: sat_shv          !
   real                             :: min_leaf_water   !
   real                             :: max_leaf_water   !
   real                             :: maxfluxrate      !
   real                             :: intercepted_tot  !
   real                             :: qintercepted_tot !
   real                             :: dintercepted_tot !
   real                             :: intercepted      !
   real                             :: qintercepted     !
   real                             :: qwflxvc          !
   real                             :: qtransp          !
   real                             :: water_demand     !
   real                             :: water_supply     !
   real                             :: veg_temp_sat     !
   !----- Functions -----------------------------------------------------------------------!
   real   , external                :: vertical_vel_flux  !
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Computing the fluxes from atmosphere to canopy.                                    !
   !---------------------------------------------------------------------------------------!
   rho_ustar = rhos      * initp%ustar         ! Aux. variable
   hflxac    = rho_ustar * initp%tstar * exner ! Sensible Heat flux
   wflxac    = rho_ustar * initp%rstar         ! Water flux
   cflxac    = rho_ustar * initp%cstar         ! CO2 flux
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Surface roughness parameters. Eventually I should account for snow factors here.  !
   !---------------------------------------------------------------------------------------!
   zoveg = csite%veg_rough(ipa)
   zveg = csite%veg_height(ipa)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Capacities of the canopy air space.                                               !
   !---------------------------------------------------------------------------------------!
   wcapcan = rhos * max(zveg,3.5)
   wcapcani = 1.0 / wcapcan
   hcapcani = cpi * wcapcani
   !---------------------------------------------------------------------------------------!

  
   !---------------------------------------------------------------------------------------!
   !     The following value of ground-canopy resistance for the nonvegetated (bare soil   !
   ! or water) surface is from John Garratt.  It is 5/ustar and replaces the one from old  !
   ! leaf.                                                                                 !
   !---------------------------------------------------------------------------------------!
   if(debug .and. abs(initp%ustar) < tiny(1.0)) print*,"USTAR = 0"
   rasgnd = 5. / initp%ustar

   !---------------------------------------------------------------------------------------!
   !Calculate fraction of open canopy                                                      !
   !---------------------------------------------------------------------------------------!
   cpatch => csite%patch(ipa)
   can_frac = 1.0
   do ico = 1,cpatch%ncohorts
      if(cpatch%lai(ico) > lai_min) then
         can_frac = can_frac                                                               &
                  * (1.0 - min(1.0                                                         &
                              ,cpatch%nplant(ico)*dbh2ca(cpatch%dbh(ico),cpatch%pft(ico))))
      end if
   end do
   can_frac = 1.0 - can_frac
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !   Obtaining the shedding water and t
   !---------------------------------------------------------------------------------------!
 
   if (csite%lai(ipa) > lai_min) then
      !------------------------------------------------------------------------------------!
      !     If vegetation is sufficiently abundant and not covered by snow, compute heat   !
      ! and moisture fluxes from vegetation to canopy, and flux resistance from soil or    !
      ! snow to canopy.                                                                    !
      !------------------------------------------------------------------------------------!
      c2 = max(0.,min(1., 0.509 * csite%lai(ipa)))
      rd = rasgnd * (1. - c2) + initp%rasveg * c2

      laii = 1.0/csite%lai(ipa)
      !------------------------------------------------------------------------------------!
      !    If the canopy does not cover all of the ground, then it should not intercept    !
      ! all of the water.                                                                  !
      !------------------------------------------------------------------------------------!
      if (pcpg>0.0) then
         !----- Scale interception by canopy openess (MCD 01-12-09). ----------------------!
         intercepted_tot  = pcpg  * can_frac
         qintercepted_tot = qpcpg * can_frac
         dintercepted_tot = dpcpg * can_frac
         !----- Energy and mass are extensive, this guarantees conservation ---------------!
         wshed_tot    = pcpg  - intercepted_tot
         qwshed_tot   = qpcpg - qintercepted_tot
         dwshed_tot   = dpcpg - dintercepted_tot
      else
         !----- No precipitation, nothing to be intercepted... ----------------------------!
         intercepted_tot  = 0.0
         qintercepted_tot = 0.0
         dintercepted_tot = 0.0
         wshed_tot        = 0.0
         qwshed_tot       = 0.0
         dwshed_tot       = 0.0
      end if

   else
      !------------------------------------------------------------------------------------!
      !     If the LAI is very small or if snow mostly covers the vegetation, bypass vege- !
      ! tation computations.  Set heat and moisture flux resistance rd  between the        !
      ! "canopy" and snow or soil surface to its bare soil value. Set shed precipitation   !
      ! heat and moisture to unintercepted values.                                         !
      !------------------------------------------------------------------------------------!
      rd               = rasgnd
      wshed_tot        = pcpg
      qwshed_tot       = qpcpg
      dwshed_tot       = dpcpg
      intercepted_tot  = 0.0
      qintercepted_tot = 0.0
      dintercepted_tot = 0.0

      !------------------------------------------------------------------------------------!
      ! Note: If the condition of low LAI for the entire patch was met, then it does not   !
      !       matter what the individual cohorts are normalized by, because they are       !
      !       effectively zero. So make sure the inverse patch LAI is a nominal non-zero/  !
      !       non-infinite number. This will only be used when parsing out intercepted     !
      !       leaf water into shed water; in which case the intercepted water is zero      !
      !       anyway. So this is just to prevent FPEs.                                     !
      !------------------------------------------------------------------------------------!
      laii = 1.0
   end if
 
   !----- Assign conductivity coefficient -------------------------------------------------!
   rdi = rhos / rd
  
   !---------------------------------------------------------------------------------------!
   !     Compute sensible heat and moisture fluxes between top soil or snow surface and    !
   ! canopy air.  wflxgc [kg/m2/s] is the upward vapor flux from soil or snow evaporation  !
   ! and dewgnd is the mass of dew that forms on the snow/soil surface this timestep; both !
   ! are defined as always positive or zero.                                               !
   !---------------------------------------------------------------------------------------!
   if (initp%nlev_sfcwater == 0) then
      hflxgc = cp * (initp%soil_tempk(nzg) - initp%can_temp) * rdi
   else
      hflxgc = cp * (initp%sfcwater_tempk(initp%nlev_sfcwater) - initp%can_temp) * rdi
   end if
  
   wflx  = (initp%surface_ssh - initp%can_shv) * rdi
   qwflx = wflx * (alvi - surface_fliq * alli)
   !---------------------------------------------------------------------------------------!


 
   !---------------------------------------------------------------------------------------!
   !     Calculate the dew flux.  Here the decision on whether dew or frost will form      !
   ! depends on the surface_temperature.                                                   !
   !---------------------------------------------------------------------------------------!
   dewgndflx  = max(0.0, -wflx)
   qdewgndflx = dewgndflx * (alvi - surface_fliq * alli)
   !---------------------------------------------------------------------------------------!
   !    I know this is a lame way to define frost density, however I couldn't find a good  !
   ! parametrisation for frost over leaves (just a bunch of engineering papers on frost    !
   ! formation over flat metal surfaces), so I decided to keep it simple and stupid. I     !
   ! will leave this as the first attempt, so if you know a better way to do it, feel free !
   ! to add it here.                                                                       !
   !---------------------------------------------------------------------------------------!
   ddewgndflx = dewgndflx / (surface_fliq * wdns + (1.-surface_fliq) * idns)


   !----- Temporary water/snow layers exist. ----------------------------------------------!
   if (initp%nlev_sfcwater > 0) then
      wflxgc  = max(0.,wflx)
      qwflxgc = max(0.,qwflx)
   !----- No surface water and not dry: evaporate from soil pores -------------------------!
   else if (initp%soilair01(nzg) > 0.0 ) then
      wflxgc = max( 0.0, (initp%ground_shv - initp%can_shv) * rdi)
      !----- Adjusting the flux accordingly to the surface fraction (no phase bias) -------!
      qwflxgc = wflxgc * ( alvi - surface_fliq * alli)
   !----- No surface water and really dry: don't evaporate at all -------------------------!
   else 
      wflxgc  = 0.0
      qwflxgc = 0.0
   endif
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Loop over the cohorts in the patch. Calculate energy fluxes with surrounding      !
   ! canopy air space, integrate cohort energy, calculate precipitation throughfall and    !
   ! sum fluxes to the patch level. Initialize variables used to store sums over cohorts.  !
   !---------------------------------------------------------------------------------------!
   hflxvc_tot  = 0.0               
   wflxvc_tot  = 0.0               
   qwflxvc_tot = 0.0               
   cflxvc_tot  = csite%cwd_rh(ipa) 
   transp_tot  = 0.0               
   cflxgc      = csite%rh(ipa) - csite%cwd_rh(ipa)
   gpp_tot     = 0.0
  
   cohortloop: do ico = 1,cpatch%ncohorts
      
      cflxgc = cflxgc + cpatch%root_respiration(ico)
      
      !------------------------------------------------------------------------------------!
      !    Calculate 'decay' term of storage (same for all) need to convert units from     !
      ! kgC/plant/day to umolC/m2/s.                                                       !
      !------------------------------------------------------------------------------------!
      storage_decay = ( cpatch%growth_respiration(ico) + cpatch%storage_respiration(ico)   &
                      + cpatch%vleaf_respiration(ico)) & 
                    * cpatch%nplant(ico) / (day_sec * umol_2_kgC)
      cflxvc_tot    = cflxvc_tot + storage_decay
      
      
      !------------------------------------------------------------------------------------!
      !     See if this cohort has leaves, if not set the leaf energy derivatives to zero, !
      ! and pass all throughfall to the ground. Later, those small cohorts will have their !
      ! leaf energy set to equilibrium with the canopy air space (temperature).            !
      !------------------------------------------------------------------------------------!
      if (cpatch%lai(ico) > lai_min .and. cpatch%hite(ico) > csite%total_snow_depth(ipa))  &
      then

         !------ Defining some shortcuts to indices ---------------------------------------!
         ipft  = cpatch%pft(ico)
         kroot = cpatch%krdepth(ico)


         !------  Calculate leaf-level CO2 flux -------------------------------------------!
         leaf_flux = cpatch%gpp(ico) - cpatch%leaf_respiration(ico)
         
         !------ Update CO2 flux from vegetation to canopy air space. ---------------------!
         cflxvc_tot = cflxvc_tot - leaf_flux

         !---------------------------------------------------------------------------------!
         !     Defining the minimum leaf water to be considered, and the maximum amount    !
         ! possible.                                                                       !
         !---------------------------------------------------------------------------------!
         min_leaf_water = min_veg_lwater*cpatch%lai(ico)
         max_leaf_water = max_veg_lwater*cpatch%lai(ico)
         
         !------ Calculate fraction of leaves covered with water. -------------------------!
         if(initp%veg_water(ico) > min_leaf_water)then
            sigmaw = min(1.0, (initp%veg_water(ico)/max_leaf_water)**twothirds )
         else
            sigmaw = 0.0
         end if


         !---------------------------------------------------------------------------------!
         !     Finding the saturation mixing ratio associated with leaf temperature.  The  !
         ! minimum is set to one to avoid FPE errors, the step will be rejected should     !
         ! this happen.                                                                    ! 
         !---------------------------------------------------------------------------------!
         veg_temp_sat = max(toocold,initp%veg_temp(ico))
         sat_shv=rslif(prss,veg_temp_sat)

         c3 = cpatch%lai(ico) * rhos * (sat_shv - initp%can_shv)
         rbi = 1.0 / cpatch%rb(ico)


         !---------------------------------------------------------------------------------!
         !    Computing the evapotranspiration or dew/frost deposition.                    !
         !---------------------------------------------------------------------------------!
         if (c3 >= 0.) then  
            !------------------------------------------------------------------------------!
            !    Evapotranspiration                                                        !
            !------------------------------------------------------------------------------!
            !------ Evaporation, energy is scaled by liquid/ice partition (no phase bias) -!
            wflxvc  = c3 * sigmaw * evap_area_one * rbi
            qwflxvc = wflxvc * (alvi - initp%veg_fliq(ico) * alli)

            cpatch%Psi_open(ico)   = c3 / (cpatch%rb(ico) + cpatch%rsw_open(ico)  )
            cpatch%Psi_closed(ico) = c3 / (cpatch%rb(ico) + cpatch%rsw_closed(ico))


            if(initp%available_liquid_water(kroot) > 0.0) then
               transp = cpatch%fs_open(ico) * cpatch%Psi_open(ico)                         &
                      + (1.0 - cpatch%fs_open(ico)) * cpatch%Psi_closed(ico)
           else
              transp = 0.0
            end if
            qtransp = transp * alvl
         else
            !------ Dew/frost formation ---------------------------------------------------!
            wflxvc                 = c3 * evap_area_one * rbi
            qwflxvc                = wflxvc * (alvi - initp%veg_fliq(ico)*alli)
            transp                 = 0.0
            qtransp                = 0.0
            cpatch%Psi_open(ico)   = 0.0
            cpatch%Psi_closed(ico) = 0.0
         end if

         !----- Diagnostic ----------------------------------------------------------------!
         dinitp%ebudget_latent = dinitp%ebudget_latent + qwflxvc + qtransp


         !----- We need to extract water from the soil equal to the transpiration. --------!
         initp%extracted_water(kroot) = initp%extracted_water(kroot) + transp


         !---------------------------------------------------------------------------------!
         !   Calculate vegetation-to-canopy sensible heat flux.                            !
         !---------------------------------------------------------------------------------!
         hflxvc = evap_area_two * cpatch%lai(ico) * cp * rhos * rbi                        &
                    * (initp%veg_temp(ico) - initp%can_temp)

         !---------------------------------------------------------------------------------!
         !     Calculate interception by leaves. Added RGK 11-2008, comments welcomed      !
         !                                                                                 !
         !  wflxvc accounts for evaporation and dew formation.  If the leaf has more water !
         ! than the carrying capacity, then it must flux all precipitation and dew. The    !
         ! leaf water may evaporate in every condition.                                    !
         !---------------------------------------------------------------------------------!
         if (initp%veg_water(ico) >= max_leaf_water) then
            !------------------------------------------------------------------------------!
            ! Case 1: Leaf has no space for rain. All rain/snow falls with the same        !
            !         density it fell. Dew and frost and old precipitation that were       !
            !         already there will likewise fall to bring it to the maximum amount   !
            !         of leaf water.                                                       !
            !------------------------------------------------------------------------------!
            wshed                 = intercepted_tot  * cpatch%lai(ico) * laii
            qwshed                = qintercepted_tot * cpatch%lai(ico) * laii
            dwshed                = dintercepted_tot * cpatch%lai(ico) * laii

            intercepted           = 0.
            qintercepted          = 0.
         else
            !------------------------------------------------------------------------------!
            ! Case 2: Leaf has space for rain. Rainfall and its internal energy accumulate !
            !         on the leaf.                                                         !
            !------------------------------------------------------------------------------!
            wshed                 = 0.0
            qwshed                = 0.0
            dwshed                = 0.0
            intercepted           = intercepted_tot  * cpatch%lai(ico) * laii
            qintercepted          = qintercepted_tot * cpatch%lai(ico) * laii
         end if
         
         dinitp%veg_water(ico) = - wflxvc + intercepted

         !---------------------------------------------------------------------------------!
         !     Find the energy balance for this cohort.                                    !
         !---------------------------------------------------------------------------------!
         dinitp%veg_energy(ico) = cpatch%rshort_v(ico) & ! Absorbed short wave radiation
                                + cpatch%rlong_v(ico)  & ! Net thermal radiation
                                - hflxvc               & ! Sensible heat flux
                                - qwflxvc              & ! Evaporative
                                - qtransp              & ! Transpiration.
                                + qintercepted         ! ! Intercepted water energy


         !----- Add the contribution of this cohort to total heat and evapotranspiration. -!
         wflxvc_tot  = wflxvc_tot  + wflxvc
         qwflxvc_tot = qwflxvc_tot + qwflxvc
         hflxvc_tot  = hflxvc_tot  + hflxvc
         transp_tot  = transp_tot  + transp

         !---------------------------------------------------------------------------------!
         ! wshed:  Water passing through vegetated canopy to soil surface                  !
         !         (enters virtual layer first), [kg/m2/s]                                 !
         !---------------------------------------------------------------------------------!
         wshed_tot  = wshed_tot  + wshed 
         qwshed_tot = qwshed_tot + qwshed
         dwshed_tot = dwshed_tot + dwshed

      else
         !---------------------------------------------------------------------------------! 
         !     If there are no leaves, leaf fluxes and interception don't exist.           !
         !---------------------------------------------------------------------------------!
         dinitp%veg_energy(ico) = 0.0
         dinitp%veg_water(ico)  = 0.0

         !---------------------------------------------------------------------------------!
         !     Allow the complete bypass of precipitation if there are no leaves.          !
         !                                                                                 !
         ! Added RGK 11-2008, comments welcomed:                                           !
         !          This will cause small deviations if the patch has LAI larger than      !
         !      lai_min, but only some of the cohorts do. If the ones that do not exceed   !
         !      the threshold of lai_min make up an appreciable fraction of the total LAI, !
         !      then we better account for it.                                             !
         !---------------------------------------------------------------------------------!
         wshed_tot  = wshed_tot  + intercepted_tot*cpatch%lai(ico)  * laii
         qwshed_tot = qwshed_tot + qintercepted_tot*cpatch%lai(ico) * laii
         dwshed_tot = dwshed_tot + dintercepted_tot*cpatch%lai(ico) * laii

      end if

   end do cohortloop


   !---------------------------------------------------------------------------------------!
   !     Update temperature and moisture of canopy.  hcapcan [J/m2/K] and wcapcan          !
   ! [kg_air/m2] are the heat and moisture capacities of the canopy.                       !
   !---------------------------------------------------------------------------------------!
   dinitp%can_temp = (hflxgc + hflxvc_tot + hflxac) * hcapcani
   dinitp%can_shv  = (wflxgc - dewgndflx + wflxvc_tot + transp_tot +  wflxac) * wcapcani


   !----- Update CO2 concentration in the canopy ------------------------------------------!
   dinitp%can_co2 = ( (cflxgc + cflxvc_tot)*mmdry + cflxac) * wcapcani


   !---------------------------------------------------------------------------------------!
   !     Integrate diagnostic variables - These are not activated unless fast file-type    !
   ! outputs are selected. This will speed up the integrator.                              !
   !---------------------------------------------------------------------------------------!
   if (fast_diagnostics) then

      dinitp%wbudget_loss2atm = - wflxac
      dinitp%ebudget_loss2atm = - hflxac
      dinitp%ebudget_latent   = dinitp%ebudget_latent -qdewgndflx + qwflxgc
      
      dinitp%co2budget_loss2atm = - cflxac * mmdryi
      dinitp%avg_gpp            = gpp_tot
      dinitp%avg_carbon_ac      = cflxac * mmdryi

      dinitp%avg_sensible_vc   = hflxvc_tot                     ! Sens. heat,  Leaf->Canopy
      dinitp%avg_sensible_2cas = hflxgc+hflxac+hflxvc_tot       ! Sens. heat,  All ->Canopy
      dinitp%avg_vapor_vc      = qwflxvc_tot                    ! Lat.  heat,  Leaf->Canopy
      dinitp%avg_sensible_gc   = hflxgc                         ! Sens. heat,  Gnd ->Canopy
      dinitp%avg_sensible_ac   = hflxac / exner                 ! Sens. heat,  Atmo->Canopy
      dinitp%avg_vapor_ac      = alvl*wflxac                    ! Lat.  heat,  Atmo->Canopy
      dinitp%avg_transp        = alvl*transp_tot                ! Transpiration
      dinitp%avg_evap          = qwflxgc-qdewgndflx+qwflxvc_tot ! Evaporation/Condensation
      dinitp%avg_sensible_tot  = (hflxgc + hflxvc_tot)          ! Total Sensible heat
      if (initp%nlev_sfcwater > 0) then
         dinitp%avg_netrad = csite%rlong_g(ipa) + csite%rlong_s(ipa) + csite%rshort_g(ipa) &
                           + sum(csite%rshort_s(1:initp%nlev_sfcwater,ipa),1)
      else
         dinitp%avg_netrad = csite%rlong_g(ipa) + csite%rlong_s(ipa) + csite%rshort_g(ipa)
      end if
   end if
  
   !---------------------------------------------------------------------------------------!
   !     These variables below are virtual copies of the variables above, but are here for !
   ! for use in the coupled model. They form the set of canopy-atmospher fluxes that are   !
   ! used for turbulent closure. These variables are also zeroed and normalized every      !
   ! dtlsm timestep, the others are likely averaged over the analysis period.              ! 
   !---------------------------------------------------------------------------------------!
   dinitp%upwp = -(initp%ustar**2)
   dinitp%rpwp = -(initp%ustar*initp%rstar)
   dinitp%tpwp = -(initp%ustar*initp%tstar)
   if(debug .and. abs(atm_tmp) < tiny(1.0)) print*,"atm_tmp = 0"
   dinitp%wpwp = vertical_vel_flux(grav * geoht * cpi * exner / atm_tmp &
                                  ,initp%tstar,initp%ustar)

   !----- If the single pond layer is too thin, force equilibrium with top soil layer -----!
   !if (initp%nlev_sfcwater == 1) then
   !   if (abs(initp%sfcwater_mass(1)) < water_stab_thresh) then
   !      dinitp%sfcwater_energy(nzg) = dinitp%sfcwater_energy(nzg)                         &
   !                                  + dinitp%sfcwater_energy(1) * dslzi(nzg)
   !      dinitp%sfcwater_energy(1) = 0.
   !   end if
   !end if

   return
end subroutine canopy_derivs_two_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the characteristic scales based on  Louis (1981) surface      !
! layer parameterization.                                                                  !
!------------------------------------------------------------------------------------------!
subroutine ed_stars(tha,rva,chia,thv,zpm,um,rough,ustar,rstar,tstar,cstar,can_shv,can_co2)

   use consts_coms     , only : grav   & ! intent(in)
                              , vonk   ! ! intent(in)
   use canopy_air_coms , only : ustmin & ! intent(in)
                              , ubmin  ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real, intent(in)  :: can_shv
   real, intent(in)  :: can_co2
   real, intent(in)  :: rough
   real, intent(out) :: ustar
   real, intent(out) :: rstar
   real, intent(out) :: tstar
   real, intent(out) :: cstar
   !----- Local variables -----------------------------------------------------------------!
   real              :: a2  ! Drag coefficient in neutral conditions, here same for h/m
   real              :: c1
   real              :: ri  ! The bulk richardson numer, eq. 3.45 in Garratt
   real              :: fm
   real              :: fh
   real              :: c2
   real              :: cm
   real              :: ch
   real              :: c3
   real              :: thv
   real              :: rva
   real              :: tha
   real              :: zpm
   real              :: um
   real              :: chia
   real              :: vels_pat
   !----- Constants -----------------------------------------------------------------------!
   real, parameter   :: b   = 5.0
   real, parameter   :: csm = 7.5
   real, parameter   :: csh = 5.0
   real, parameter   :: d   = 5.0
   !---------------------------------------------------------------------------------------!

   vels_pat = max(um,ubmin)
  
  
   a2 = (vonk / log(zpm / rough)) ** 2
   c1 = a2 * vels_pat

   ri = grav * zpm * (tha - thv)  / (.5 * (tha + thv) * vels_pat * vels_pat)
   if (tha - thv > 0.) then
      !----- Stable case ------------------------------------------------------------------!
      fm = 1. / (1. + (2. * b * ri / sqrt(1. + d * ri)))
      fh = 1. / (1. + (3. * b * ri * sqrt(1. + d * ri)))
   else
      !----- Unstable case ----------------------------------------------------------------!
      c2 = b * a2 * sqrt(zpm / rough * (abs(ri)))
      cm = csm * c2
      ch = csh * c2
      fm = (1. - 2. * b * ri / (1. + 2. * cm))
      fh = (1. - 3. * b * ri / (1. + 3. * ch))
   end if

   ustar = max(ustmin,sqrt(c1 * vels_pat * fm))
   c3 = c1 * fh / ustar
   rstar = c3 * (rva - can_shv)
   tstar = c3 * (tha - thv)
   cstar = c3 * (chia - can_co2)

   return
end subroutine ed_stars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
real function vertical_vel_flux(gzotheta,tstar,ustar)
  
   use consts_coms , only : vonk ! intent(in)
  
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real, intent(in)    :: ustar
   real, intent(in)    :: tstar
   real, intent(in)    :: gzotheta
   !----- Local variables -----------------------------------------------------------------!
   real :: zoverl
   real :: cx
   real :: psin
   !----- Constants -----------------------------------------------------------------------!
   real, parameter     :: wtol = 1.e-20
   !---------------------------------------------------------------------------------------!
  
  
   zoverl = gzotheta * vonk * tstar / (ustar * ustar)
  
   if (zoverl < 0.)then
      cx = zoverl * sqrt(sqrt(1. - 15. * zoverl))
   else
      cx = zoverl / (1.0 + 4.7 * zoverl)
   endif
  
   psin = sqrt((1.-2.86 * cx) / (1. + cx * (-5.39 + cx * 6.998 )))
   vertical_vel_flux = (0.27 * max(6.25 * (1. - cx) * psin,wtol) - 1.18 * cx * psin)       &
                     * ustar * ustar
  
   return
end function vertical_vel_flux
!==========================================================================================!
!==========================================================================================!

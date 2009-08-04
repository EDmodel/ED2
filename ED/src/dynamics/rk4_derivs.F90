!==========================================================================================!
!==========================================================================================!
! Subroutine leaf_derivs                                                                !
!                                                                                          !
!     This subroutine finds the fast-scale derivatives at canopy, soil, and leaf surface.  !
! This subroutine is based on LEAF-3, except that here only the derivative is computed,    !
! whereas in LEAF-3 the actual step is done at once. This derivative will be used for the  !
! Runge-Kutta integration step.                                                            !
!------------------------------------------------------------------------------------------!
subroutine leaf_derivs(initp,dinitp,csite,ipa,isi,ipy)
  
   use rk4_coms               , only : rk4met             & ! intent(in)
                                     , rk4patchtype       ! ! structure
   use ed_state_vars          , only : sitetype           ! ! structure
   use consts_coms            , only : cp8                & ! intent(in)
                                     , cpi8               ! ! intent(in)
   use grid_coms              , only : nzg                ! ! intent(in)
   use canopy_struct_dynamics , only : canopy_turbulence8 ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: initp     ! Structure with RK4 intermediate state
   type(rk4patchtype) , target     :: dinitp    ! Structure with RK4 derivatives
   type(sitetype)     , target     :: csite     ! This site (with previous values);
   integer            , intent(in) :: ipa       ! Patch ID
   integer            , intent(in) :: isi       ! Site ID
   integer            , intent(in) :: ipy       ! Polygon ID
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   ! Depending on the type of compilation, interfaces must be explicitly declared.         !
   !---------------------------------------------------------------------------------------!
#if USE_INTERF
   interface
      !------------------------------------------------------------------------------------!
      !    Subroutine that computes the canopy and leaf fluxes.                            ! 
      !------------------------------------------------------------------------------------!
      subroutine leaftw_derivs(initp,dinitp,csite,ipa,isi,ipy)
         use rk4_coms      , only : rk4patchtype ! ! structure
         use ed_state_vars , only : sitetype     ! ! structure
         implicit none
         type(rk4patchtype) , target     :: initp  
         type(rk4patchtype) , target     :: dinitp 
         type(sitetype)     , target     :: csite
         integer            , intent(in) :: ipa,isi,ipy
      end subroutine leaftw_derivs
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!

   !----- Resetting the energy budget. ----------------------------------------------------!
   dinitp%ebudget_latent = 0.0d0

   !----- Compute canopy turbulence properties. -------------------------------------------!
   call canopy_turbulence8(csite,initp,isi,ipa,.true.)

   !----- Finding the derivatives. --------------------------------------------------------!
   call leaftw_derivs(initp,dinitp,csite,ipa,isi,ipy)

   !----- Nlev_sfcwater derivative... I doubt it's really used... -------------------------!
   dinitp%nlev_sfcwater = initp%nlev_sfcwater

   return
end subroutine leaf_derivs
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine leaftw_derivs(initp,dinitp,csite,ipa,isi,ipy)
   use ed_max_dims             , only : nzgmax               & ! intent(in)
                                   , nzsmax               ! ! intent(in)
   use consts_coms          , only : alvl8                & ! intent(in)
                                   , cliqvlme8            & ! intent(in)
                                   , tsupercool8          & ! intent(in)
                                   , wdns8                & ! intent(in)
                                   , wdnsi8               ! ! intent(in)
   use grid_coms            , only : nzg                  & ! intent(in)
                                   , nzs                  ! ! intent(in)
   use soil_coms            , only : soil8                & ! intent(in)
                                   , slz8                 & ! intent(in)
                                   , dslz8                & ! intent(in)
                                   , dslzi8               & ! intent(in)
                                   , infiltration_method  & ! intent(in)
                                   , dslzti8              & ! intent(in)
                                   , slcons18             & ! intent(in)
                                   , slzt8                & ! intent(in)
                                   , ss                   & ! intent(in)
                                   , isoilbc              ! ! intent(in)
   use ed_misc_coms            , only : dtlsm                & ! intent(in)
                                   , current_time         ! ! intent(in)
   use rk4_coms             , only : rk4eps               & ! intent(in)
                                   , rk4min_sfcwater_mass & ! intent(in)
                                   , any_solvable         & ! intent(in)
                                   , rk4met               & ! intent(in)
                                   , rk4patchtype         ! ! structure
   use ed_state_vars        , only : sitetype             & ! structure
                                   , patchtype            ! ! structure
   use ed_therm_lib         , only : ed_grndvap8          ! ! subroutine
   use therm_lib            , only : qtk8                 & ! subroutine
                                   , qwtk8                ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)  , target     :: initp         ! RK4 structure, intermediate step
   type(rk4patchtype)  , target     :: dinitp        ! RK4 structure, derivatives
   type(sitetype)      , target     :: csite         ! Current site (before integration)
   integer             , intent(in) :: ipa           ! Current patch ID
   integer             , intent(in) :: isi           ! Current site ID
   integer             , intent(in) :: ipy           ! Current polygon ID
   !----- Local variables -----------------------------------------------------------------!
   integer                          :: k, k1, k2     ! Level counters
   integer                          :: ksn           ! # of temporary water/snow layers
   integer                          :: nsoil         ! Short for csite%soil_text(k,ipa)
   real(kind=8)                     :: wgpfrac       ! Fractional soil moisture
   real(kind=8)                     :: soilcond      ! Soil conductivity
   real(kind=8)                     :: snden         ! Snow/water density
   real(kind=8)                     :: hflxgc        ! Ground -> canopy heat flux
   real(kind=8)                     :: wflxgc        ! Ground -> canopy water flux
   real(kind=8)                     :: qwflxgc       ! Ground -> canopy latent heat flux
   real(kind=8)                     :: dewgnd        ! Dew/frost flux to ground
   real(kind=8)                     :: qdewgnd       ! Dew/frost heat flux to ground
   real(kind=8)                     :: ddewgnd       ! Dew/frost density flux to ground
   real(kind=8)                     :: wshed         ! Water shedding flux
   real(kind=8)                     :: qwshed        ! Energy flux due to water shedding
   real(kind=8)                     :: dwshed        ! Density flux due to water shedding
   real(kind=8)                     :: wgpmid        ! Soil in between layers
   real(kind=8)                     :: wloss         ! Water loss due to transpiration
   real(kind=8)                     :: qwloss        ! Energy loss due to transpiration
   real(kind=8)                     :: dqwt          ! Energy adjustment aux. variable
   real(kind=8)                     :: fracliq       ! Fraction of liquid water
   real(kind=8)                     :: tempk         ! Temperature
   real(kind=8)                     :: qwgoal        ! Goal energy for thin snow layers.
   real(kind=8)                     :: wprevious     ! Previous water content
   real(kind=8)                     :: infilt        ! Surface infiltration rate
   real(kind=8)                     :: qinfilt       ! Surface infiltration heat rate
   real(kind=8)                     :: snowdens      ! Snow density (kg/m2)
   real(kind=8)                     :: soilhcap      ! Soil heat capacity
   real(kind=8)                     :: int_sfcw_u    ! Intensive sfc. water internal en.
   real(kind=8)                     :: surface_water ! Temp. variable. Availible liquid 
                                                     !   water on the soil surface (kg/m2)
   real(kind=8)                     :: freezeCor     ! Correction to conductivity for 
                                                     !    partially frozen soil.
   real(kind=8), dimension(nzgmax+nzsmax)   :: rfactor       ! 
   real(kind=8), dimension(nzgmax+nzsmax+1) :: hfluxgsc      ! Surface -> canopy heat flux
   real(kind=8), dimension(nzg+nzs+1)       :: w_flux        ! Water flux (aux. variable)
   real(kind=8), dimension(nzg+nzs+1)       :: qw_flux       ! Heat flux (aux. variable)
   real(kind=8), dimension(nzs+1)           :: d_flux        ! Density flux
   !----- Constants -----------------------------------------------------------------------!
   logical     , parameter  :: debug = .false.   ! Debugging output flag (T/F)
   real(kind=8), parameter  :: freezeCoef = 7.d0 ! Exponent in the frozen soil hydraulic 
                                                 !    conductivity correction.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! Depending on the type of compilation, interfaces must be explicitly declared.         !
   !---------------------------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine canopy_derivs_two(initp,dinitp,csite,ipa,isi,ipy,hflxgc,wflxgc,qwflxgc &
                                     ,dewgndflx,qdewgndflx,ddewgndflx,wshed_tot,qwshed_tot &
                                     ,dwshed_tot)
         use rk4_coms     , only: rk4patchtype  ! ! structure
         use ed_state_vars, only: sitetype      & ! structure
                                , patchtype     ! ! structure
         implicit none
         type (rk4patchtype) , target      :: initp, dinitp
         type (sitetype)     , target      :: csite
         integer             , intent(in)  :: ipa, isi, ipy
         real(kind=8)        , intent(out) :: hflxgc, wflxgc, qwflxgc
         real(kind=8)        , intent(out) :: dewgndflx, qdewgndflx, ddewgndflx
         real(kind=8)        , intent(out) :: wshed_tot, qwshed_tot, dwshed_tot
      end subroutine canopy_derivs_two
   end interface
#endif
   !---------------------------------------------------------------------------------------!

   !----- Copying the # of surface water/snow layers to a shortcut ------------------------!
   ksn = initp%nlev_sfcwater
  
   !---- Initializing some vertical integration variables --------------------------------!
   w_flux  = 0.0d0
   qw_flux = 0.0d0
   d_flux  = 0.0d0

   !----- Initialize derivatives to zero -------------------------------------------------!
   dinitp%soil_energy(:)     = 0.0d0
   dinitp%soil_water(:)      = 0.0d0
   dinitp%sfcwater_depth(:)  = 0.0d0
   dinitp%sfcwater_energy(:) = 0.0d0
   dinitp%sfcwater_mass(:)   = 0.0d0
   dinitp%virtual_heat       = 0.0d0
   dinitp%virtual_water      = 0.0d0
   dinitp%virtual_depth      = 0.0d0
   initp%extracted_water(:)  = 0.0d0
   dinitp%avg_smoist_gc(:)   = 0.0d0

   !----- Copying the soil texture flag to a shortcut -------------------------------------!
   nsoil = csite%ntext_soil(nzg,ipa)
   k = max(1,ksn)
   if (abs(initp%sfcwater_mass(k)) > rk4min_sfcwater_mass) then
      int_sfcw_u = initp%sfcwater_energy(k)/initp%sfcwater_mass(k)
   else
      int_sfcw_u = 0.0d0
   end if
   
   call ed_grndvap8(ksn,nsoil,initp%soil_water(nzg),initp%soil_energy(nzg),int_sfcw_u      &
                   ,rk4met%rhos,initp%can_shv,initp%ground_shv,initp%surface_ssh           &
                   ,initp%surface_temp,initp%surface_fliq)

   !---------------------------------------------------------------------------------------!
   !     Calculate water available to vegetation (in meters). SLZ is specified in RAMSIN.  !
   ! Each element of the array sets the value of the bottom of a corresponding soil layer. !
   ! Eg, SLZ = -2, -1, -0.5, -0.25.  There are four soil layers in this example; soil      !
   ! layer 1 goes from 2 meters below the surface to 1 meter below the surface.            !
   !---------------------------------------------------------------------------------------!
   initp%available_liquid_water(nzg) = dslz8(nzg)                                          &
                                     * max(0.0d0                                           &
                                          ,initp%soil_fracliq(nzg)*(initp%soil_water(nzg)  &
                                          -soil8(nsoil)%soilcp))

   do k = nzg - 1, rk4met%lsl, -1
      nsoil = csite%ntext_soil(k,ipa)
      initp%available_liquid_water(k) = initp%available_liquid_water(k+1) + dslz8(k)       &
           *max(0.0d0,(initp%soil_water(k)-soil8(nsoil)%soilcp)*initp%soil_fracliq(k))
   end do

   !---------------------------------------------------------------------------------------!
   !     Compute gravitational potential plus moisture potential, psi + z (psiplusz) [m],  !
   ! liquid water content (soil_liq) [m], and 99% the remaining water capacity (soilair99) !
   ! [m].                                                                                  !
   !---------------------------------------------------------------------------------------!
   do k = rk4met%lsl, nzg
      nsoil = csite%ntext_soil(k,ipa)
      initp%psiplusz(k) = slzt8(k) + soil8(nsoil)%slpots                                   &
                        * (soil8(nsoil)%slmsts / initp%soil_water(k)) ** soil8(nsoil)%slbs

      !------------------------------------------------------------------------------------!
      !    Soil liquid water must be converted to meters of liquid water per layer. This   !
      ! requires multiplication of volumetric water content, m3(water)/m3 must be          !
      ! multiplied by depth to get a depth of water.                                       !
      !------------------------------------------------------------------------------------!
      initp%soil_liq(k) = max(0.0d0                                                        &
                             ,( initp%soil_water(k) - soil8(nsoil)%soilcp)                  &
                              * initp%soil_fracliq(k) )
      initp%soilair99(k) = (1.d0-rk4eps) * soil8(nsoil)%slmsts - initp%soil_water(k)
      initp%soilair01(k) = (initp%soil_water(k) - (1.+rk4eps) * soil8(nsoil)%soilcp)       &
                         * initp%soil_fracliq(k)
   end do
   !---------------------------------------------------------------------------------------!

 
   !----- Get derivatives of canopy variables. --------------------------------------------!
   call canopy_derivs_two(initp,dinitp,csite,ipa,isi,ipy,hflxgc,wflxgc,qwflxgc,dewgnd   &
                            ,qdewgnd,ddewgnd,wshed,qwshed,dwshed)

   !---------------------------------------------------------------------------------------!
   !     Here we check whether it is bedrock or not. I think the reason for this check is  !
   ! to avoid problems with soil_water never being defined for bedrocks, otherwise simply  !
   ! assuming soilcond1 and soilcond2 would suffice.                                       !
   !---------------------------------------------------------------------------------------!
   do k = rk4met%lsl, nzg
      nsoil = csite%ntext_soil(k,ipa)
      if(nsoil /= 13)then
         wgpfrac  = min(initp%soil_water(k) / soil8(nsoil)%slmsts,1.d0)
         soilcond = soil8(nsoil)%soilcond0                                                 &
                  + wgpfrac * (soil8(nsoil)%soilcond1 + wgpfrac *  soil8(nsoil)%soilcond2)
      else
         soilcond=soil8(nsoil)%soilcond0
      end if
      rfactor(k) = dslz8(k) / soilcond
   end do


   !----- Snow/surface water --------------------------------------------------------------!
   do k = 1, ksn
      if(initp%sfcwater_depth(k) > 0.d0)then
         snden = initp%sfcwater_mass(k) / initp%sfcwater_depth(k)
         rfactor(k+nzg) = initp%sfcwater_depth(k)                                          &
                        / (ss(1) * exp(ss(2) * initp%sfcwater_tempk(k))                    &
                        * (ss(3) + snden * (ss(4) + snden * (ss(5) + snden * ss(6)))))
      else
         rfactor(k+nzg) = 0.d0
      end if
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Calculate the sensible heat fluxes find soil and sfcwater internal sensible heat  !
   ! fluxes (hfluxgsc) [W/m2].                                                             !
   !---------------------------------------------------------------------------------------!
   hfluxgsc(:) = 0.d0
   do k = rk4met%lsl+1, nzg
      hfluxgsc(k) = - (initp%soil_tempk(k) - initp%soil_tempk(k-1))                        &
                    / ((rfactor(k) + rfactor(k-1)) * 5.d-1)    
      dinitp%avg_sensible_gg(k-1) = hfluxgsc(k)  ! Diagnostic
   end do

   !----- If temporary water/snow layers exist, compute them now... -----------------------!
   if (ksn >= 1) then
      hfluxgsc(nzg+1) = - (initp%sfcwater_tempk(1) - initp%soil_tempk(nzg))                &
                        / ((rfactor(nzg+1)   + rfactor(nzg)) * 5.d-1)

      do k = 2,ksn
         hfluxgsc(nzg+k) = - (initp%sfcwater_tempk(k) - initp%sfcwater_tempk(k-1))         &
                         /   ((rfactor(nzg+k) + rfactor(nzg+k-1)) * 5.d-1)
      end do
   end if

   !----- Heat flux (hfluxgsc) at soil or sfcwater top from longwave, sensible [W/m^2] ----!
   hfluxgsc(nzg+ksn+1) = hflxgc + qwflxgc                                                  &
                       - dble(csite%rlong_g(ipa)) - dble(csite%rlong_s(ipa))

   !----- Heat flux -----------------------------------------------------------------------!
   dinitp%avg_sensible_gg(nzg) = hfluxgsc(nzg+ksn+1) ! Diagnostic

   !---------------------------------------------------------------------------------------!
   !    Update soil U values [J/m³] from sensible heat, upward water vapor (latent heat)   !
   ! and longwave fluxes. This excludes effects of dew/frost formation, precipitation,     !
   ! shedding, and percolation.                                                            !
   !---------------------------------------------------------------------------------------!
   do k = rk4met%lsl,nzg
      dinitp%soil_energy(k) = dslzi8(k) * (hfluxgsc(k)- hfluxgsc(k+1))
   end do

   !----- Update soil Q values [J/m³] from shortwave flux. --------------------------------!
   dinitp%soil_energy(nzg) = dinitp%soil_energy(nzg)                                       &
                           + dslzi8(nzg) * dble(csite%rshort_g(ipa))


   !---------------------------------------------------------------------------------------!
   !    Update surface water U values [J/m²] from sensible heat, upward water vapor        !
   ! (latent heat), longwave, and shortwave fluxes.  This excludes effects of dew/frost    !
   ! formation, precipitation, shedding and percolation.                                   !
   !---------------------------------------------------------------------------------------!
   do k = 1,ksn
     dinitp%sfcwater_energy(k) = hfluxgsc(k+nzg) - hfluxgsc(k+1+nzg)                       &
                               + dble(csite%rshort_s(k,ipa))
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
      dinitp%sfcwater_depth(ksn)  = -d_flux(ksn+1)
   else if (w_flux(nzg+1) < 0.d0) then
      dinitp%virtual_heat  = -qw_flux(nzg+1)
      dinitp%virtual_water = -w_flux(nzg+1)
      dinitp%virtual_depth = -d_flux(1)
      qw_flux(nzg+1)       = 0.d0
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
   w_flux(nzg+1) = w_flux(nzg+1) * wdnsi8 ! now in m/s

   !---------------------------------------------------------------------------------------!
   !     Alternate surface infiltration (MCD) based on surface conductivity not capacity.  !
   !---------------------------------------------------------------------------------------!
   if (infiltration_method /= 0) then
      call fatal_error ('Running alt infiltation when we shouldn''t be',&
                       &'leaftw_derivs','rk4_derivs.F90')

      if (initp%virtual_water /= 0.d0) then  !!process "virtural water" pool
         nsoil = csite%ntext_soil(nzg,ipa)
         if (nsoil /= 13) then
            call qtk8(initp%virtual_heat/initp%virtual_water,tempk,fracliq)
            infilt = -dslzi8(nzg)* 5.d-1 * slcons18(nzg,nsoil)                             &
                   * (initp%soil_water(nzg) / soil8(nsoil)%slmsts)                         &
                     **(2.d0 * soil8(nsoil)%slbs + 3.d0)                                   &
                   * (initp%psiplusz(nzg)-initp%virtual_water/2.d3) & !diff. in pot.
                   * 5.d-1 * (initp%soil_fracliq(nzg)+ fracliq)       ! mean liquid fraction
            qinfilt = infilt * cliqvlme8 * (tempk - tsupercool8)
            !----- Adjust other rates accordingly -----------------------------------------!
            w_flux(nzg+1)  = w_flux(nzg+1) + infilt
            qw_flux(nzg+1) = qw_flux(nzg+1)+ qinfilt
            dinitp%virtual_water = dinitp%virtual_water - infilt*wdns8
            dinitp%virtual_heat  = dinitp%virtual_heat  - qinfilt
         end if
      end if  !! end virtual water pool
      if (initp%nlev_sfcwater >= 1) then !----- Process "snow" water pool -----------------! 
         call qtk8(initp%sfcwater_energy(1)/initp%sfcwater_mass(1),tempk,fracliq)
         surface_water = initp%sfcwater_mass(1)*fracliq*wdnsi8 !(m/m2)
         nsoil = csite%ntext_soil(nzg,ipa)
         if (nsoil /= 13) then
            !----- Calculate infiltration rate (m/s) --------------------------------------!
            infilt = -dslzi8(nzg) * 5.d-1 * slcons18(nzg,nsoil)                            &
                   * (initp%soil_water(nzg) / soil8(nsoil)%slmsts)                         &
                     **(2.d0 * soil8(nsoil)%slbs + 3.d0)                                   &
                   * (initp%psiplusz(nzg) - surface_water/2.d0) & !difference in potentials
                   * 5.d-1 * (initp%soil_fracliq(nzg)+ fracliq)   ! mean liquid fraction
            qinfilt = infilt * cliqvlme8 * (tempk - tsupercool8)
            !----- Adjust other rates accordingly -----------------------------------------!
            w_flux(nzg+1)             = w_flux(nzg+1)             + infilt
            qw_flux(nzg+1)            = qw_flux(nzg+1)            + qinfilt 
            dinitp%sfcwater_mass(1)   = dinitp%sfcwater_mass(1)   - infilt*wdns8
            dinitp%sfcwater_energy(1) = dinitp%sfcwater_energy(1) - qinfilt
            dinitp%sfcwater_depth(1)  = dinitp%sfcwater_depth(1)  - infilt
         end if
      end if  ! End snow water pool
   end if  !! End alternate infiltration
   !---------------------------------------------------------------------------------------!


   !----- Keep qw_flux in W/m2. -----------------------------------------------------------!
   do k = rk4met%lsl+1, nzg
      nsoil = csite%ntext_soil(k,ipa)
      if(nsoil /= 13 .and. csite%ntext_soil(k-1,ipa) /= 13)then

         wgpmid    = 5.d-1 * (initp%soil_water(k)   + initp%soil_water(k-1))
         freezeCor = 5.d-1 * (initp%soil_fracliq(k) + initp%soil_fracliq(k-1))
         if(freezeCor < 1.d0) freezeCor = 1.d1**(-freezeCoef*(1.d0-freezeCor))
         w_flux(k) = dslzti8(k) * slcons18(k,nsoil)                                        &
                   * (wgpmid / soil8(nsoil)%slmsts)**(2.d0 * soil8(nsoil)%slbs + 3.d0)     &
                   * (initp%psiplusz(k-1) - initp%psiplusz(k)) * freezeCor

        !----------------------------------------------------------------------------------!
        !      Limit water transfers to prevent over-saturation and over-depletion.        !
        !----------------------------------------------------------------------------------!
        if (w_flux(k) > 0.) then
           if (initp%soilair01(k-1) <= 0.d0 .or. initp%soilair99(k) <= 0.d0) w_flux(k)=0.d0
        else
           if (initp%soilair99(k-1) <= 0.d0 .or. initp%soilair01(k) <= 0.d0) w_flux(k)=0.d0
        end if
      end if
      !----- Only liquid water is allowed to flow, find qw_flux (W/m2) accordingly --------!
      qw_flux(k) = w_flux(k) * cliqvlme8 * (initp%soil_tempk(k) - tsupercool8)
      dinitp%avg_smoist_gg(k-1) = w_flux(k)*wdns8   ! Diagnostic
   end do

   !----- Boundary condition at the lowest soil level -------------------------------------!
   nsoil = csite%ntext_soil(rk4met%lsl,ipa)
   if (nsoil /= 13 .and. isoilbc == 1) then
      !----- Free drainage ----------------------------------------------------------------!
      wgpmid      = sngl(initp%soil_water(rk4met%lsl))
      freezeCor   = initp%soil_fracliq(rk4met%lsl)
      if(freezeCor < 1.d0) freezeCor = 1.d1**(-freezeCoef*(1.d0-freezeCor))
      w_flux(rk4met%lsl) = dslzti8(rk4met%lsl) * slcons18(rk4met%lsl,nsoil)                &
                         * (wgpmid/soil8(nsoil)%slmsts)**(2.d0 * soil8(nsoil)%slbs + 3.d0) &
                         * freezeCor

      !-----  Make it kg/s instead of m3. -------------------------------------------------!
      dinitp%avg_drainage = w_flux(rk4met%lsl) * wdns8
      
      !------------------------------------------------------------------------------------!
      !      Limit water transfers to prevent over-saturation and over-depletion.          !
      !------------------------------------------------------------------------------------!
      if (w_flux(rk4met%lsl) > 0.d0) then
         if (initp%soilair99(rk4met%lsl) <= 0.d0) w_flux(rk4met%lsl) = 0.d0
      else
         if (initp%soilair01(rk4met%lsl) <= 0.d0) w_flux(rk4met%lsl) = 0.d0
      end if
      !----- Only liquid water is allowed to flow, find qw_flux (W/m2) accordingly --------!
      qw_flux(rk4met%lsl) = w_flux(rk4met%lsl)                                             &
                          * cliqvlme8 * (initp%soil_tempk(rk4met%lsl) - tsupercool8)
   else
      !----- Bedrock, no flux accross it. -------------------------------------------------!
      w_flux(rk4met%lsl)  = 0.d0
      qw_flux(rk4met%lsl) = 0.d0
      dinitp%avg_drainage = 0.d0
   end if

   !----- Finally, update soil moisture (impose minimum value of soilcp) and soil energy. -!
   do k = rk4met%lsl,nzg
      dinitp%soil_water(k)  = dinitp%soil_water(k)                                         &
                            - dslzi8(k) * ( w_flux(k+1) -  w_flux(k)  )
      dinitp%soil_energy(k) =  dinitp%soil_energy(k)                                       &
                            - dslzi8(k) * ( qw_flux(k+1) - qw_flux(k) )
   end do

   !---- Update soil moisture and energy from transpiration/root uptake. ------------------!
   if (any_solvable) then
      do k1 = rk4met%lsl, nzg    ! loop over extracted water
         do k2=k1,nzg
            if (csite%ntext_soil(k2,ipa) /= 13) then
               if (initp%available_liquid_water(k1) > 0.d0) then
                  wloss = wdnsi8 * initp%extracted_water(k1)                               &
                        * initp%soil_liq(k2) / initp%available_liquid_water(k1)
                  dinitp%soil_water(k2) = dinitp%soil_water(k2) - dble(wloss)

                  !----- Energy: only liquid water is lost through transpiration. ---------!
                  qwloss = wloss * cliqvlme8 * (initp%soil_tempk(k2) - tsupercool8)
                  dinitp%soil_energy(k2)   = dinitp%soil_energy(k2)   - qwloss
                  dinitp%avg_smoist_gc(k2) = dinitp%avg_smoist_gc(k2) - wdns8*wloss
                  dinitp%ebudget_latent    = dinitp%ebudget_latent    + qwloss
               end if
            end if
         end do
      end do
   end if

   return
end subroutine leaftw_derivs
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine canopy_derivs_two(initp,dinitp,csite,ipa,isi,ipy,hflxgc,wflxgc,qwflxgc       &
                               ,dewgndflx,qdewgndflx,ddewgndflx,wshed_tot,qwshed_tot       &
                               ,dwshed_tot)
   use rk4_coms              , only : rk4patchtype         & ! Structure
                                    , rk4met               & ! intent(in)
                                    , debug                & ! intent(in)
                                    , toocold              & ! intent(in)
                                    , toohot               & ! intent(in)
                                    , lai_to_cover         & ! intent(in)
                                    , effarea_water        & ! intent(in)
                                    , effarea_heat         & ! intent(in)
                                    , zoveg                & ! intent(in)
                                    , zveg                 & ! intent(in)
                                    , wcapcan              & ! intent(in)
                                    , wcapcani             & ! intent(in)
                                    , ccapcani             & ! intent(in)
                                    , hcapcani             & ! intent(in)
                                    , any_solvable         & ! intent(in)
                                    , tiny_offset          & ! intent(in)
                                    , rk4water_stab_thresh & ! intent(in)
                                    , rk4dry_veg_lwater    & ! intent(in)
                                    , rk4fullveg_lwater    ! ! intent(in)
   use ed_state_vars         , only : sitetype             & ! Structure
                                    , patchtype            ! ! Structure
   use consts_coms           , only : alvl8                & ! intent(in)
                                    , cp8                  & ! intent(in)
                                    , cpi8                 & ! intent(in)
                                    , twothirds8           & ! intent(in)
                                    , day_sec8             & ! intent(in)
                                    , grav8                & ! intent(in)
                                    , alvi8                & ! intent(in)
                                    , alli8                & ! intent(in)
                                    , umol_2_kgC8          & ! intent(in)
                                    , pi18                 & ! intent(in)
                                    , mmdry8               & ! intent(in)
                                    , mmdryi8              & ! intent(in)
                                    , wdns8                & ! intent(in)
                                    , wdnsi8               & ! intent(in)
                                    , idns8                ! ! intent(in)
   use grid_coms             , only : nzg                  ! ! intent(in)
   use soil_coms             , only : soil8                & ! intent(in)
                                    , dslzi8               & ! intent(in)
                                    , dewmax               ! ! intent(in)
   use therm_lib             , only : qwtk                 & ! subroutine
                                    , rhovsil              ! ! function
   use ed_misc_coms             , only : dtlsm                ! ! intent(in)
   use ed_misc_coms          , only : fast_diagnostics     ! ! intent(in)
   use allometry             , only : dbh2ca               ! ! function
   use canopy_struct_dynamics, only : vertical_vel_flux8   ! ! function
   use pft_coms              , only : water_conductance    & ! intent(in)
                                    , q                    & ! intent(in)
                                    , qsw                  ! ! intent(in)
                                       
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)     , target      :: csite          ! Current site
   type(rk4patchtype) , target      :: initp          ! RK4 structure, state vars
   type(rk4patchtype) , target      :: dinitp         ! RK4 structure, derivatives
   integer            , intent(in)  :: ipy            ! Current polygon ID
   integer            , intent(in)  :: isi            ! Current site ID
   integer            , intent(in)  :: ipa            ! Current patch ID
   real(kind=8)       , intent(out) :: hflxgc         ! Ground->canopy sensible heat flux
   real(kind=8)       , intent(out) :: wflxgc         ! Ground->canopy water flux
   real(kind=8)       , intent(out) :: qwflxgc        ! Ground->canopy latent heat flux
   real(kind=8)       , intent(out) :: wshed_tot      ! Water shedding rate
   real(kind=8)       , intent(out) :: qwshed_tot     ! Water shedding energy flux
   real(kind=8)       , intent(out) :: dwshed_tot     ! Water shedding density flux
   real(kind=8)       , intent(out) :: dewgndflx      ! Dew/frost water flux
   real(kind=8)       , intent(out) :: qdewgndflx     ! Dew/frost heat flux
   real(kind=8)       , intent(out) :: ddewgndflx     ! Dew/frost density
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)    , pointer     :: cpatch           ! Current patch
   integer                          :: ico              ! Current cohort ID
   integer                          :: k                ! Soil layer counter
   integer                          :: ipft             ! Shortcut for PFT type
   integer                          :: kroot            ! Level in which the root bottom is
   real(kind=8)                     :: can_frac         ! total fractional canopy coverage
   real(kind=8)                     :: transp           ! Cohort transpiration
   real(kind=8)                     :: cflxac           ! Atm->canopy carbon flux
   real(kind=8)                     :: wflxac           ! Atm->canopy water flux
   real(kind=8)                     :: hflxac           ! Atm->canopy heat flux
   real(kind=8)                     :: c2               ! Coefficient (????)
   real                             :: c3lai            ! Coefficient (????)
   real(kind=8)                     :: c3tai            ! Coefficient (????)
   real(kind=8)                     :: hflxvc           ! Leaf->canopy heat flux
   real(kind=8)                     :: rasgnd           ! 
   real(kind=8)                     :: rbi              ! 
   real(kind=8)                     :: rd               !
   real(kind=8)                     :: sigmaw           !
   real(kind=8)                     :: wflxvc           !
   real(kind=8)                     :: wshed            !
   real(kind=8)                     :: qwshed           !
   real(kind=8)                     :: dwshed           !
   real(kind=8)                     :: cflxgc           !
   real(kind=8)                     :: taii             !
   real(kind=8)                     :: wflx             !
   real(kind=8)                     :: qwflx            !
   real(kind=8)                     :: hflxvc_tot       !
   real(kind=8)                     :: transp_tot       !
   real(kind=8)                     :: cflxvc_tot       !
   real(kind=8)                     :: wflxvc_tot       !
   real(kind=8)                     :: qwflxvc_tot      !
   real(kind=8)                     :: rho_ustar        !
   real(kind=8)                     :: rdi              !
   real(kind=8)                     :: storage_decay    !
   real(kind=8)                     :: leaf_flux        !
   real(kind=8)                     :: sat_shv          !
   real(kind=8)                     :: min_leaf_water   !
   real(kind=8)                     :: max_leaf_water   !
   real(kind=8)                     :: maxfluxrate      !
   real(kind=8)                     :: intercepted_tot  !
   real(kind=8)                     :: qintercepted_tot !
   real(kind=8)                     :: dintercepted_tot !
   real(kind=8)                     :: intercepted      !
   real(kind=8)                     :: qintercepted     !
   real(kind=8)                     :: qwflxvc          !
   real(kind=8)                     :: qtransp          !
   real(kind=8)                     :: water_demand     !
   real(kind=8)                     :: water_supply     !
   real(kind=8)                     :: gzotheta         !
   real(kind=8)                     :: flux_area        ! Area between canopy and plant
   real                             :: veg_temp_sat     !
   !----- Functions -----------------------------------------------------------------------!
   real        , external           :: sngloff
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Computing the fluxes from atmosphere to canopy.                                    !
   !---------------------------------------------------------------------------------------!
   rho_ustar = rk4met%rhos * initp%ustar                ! Aux. variable
   hflxac    = rho_ustar   * initp%tstar * rk4met%exner ! Sensible Heat flux
   wflxac    = rho_ustar   * initp%qstar                ! Water flux
   cflxac    = rho_ustar   * initp%cstar * mmdryi8      ! CO2 flux [umol/m2/s]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     The following value of ground-canopy resistance for the nonvegetated (bare soil   !
   ! or water) surface is from John Garratt.  It is 5/ustar and replaces the one from old  !
   ! leaf.                                                                                 !
   !---------------------------------------------------------------------------------------!
   if(debug .and. abs(initp%ustar) < tiny(1.d0)) print*,"USTAR = 0"
   rasgnd = 5.d0 / initp%ustar

   !---------------------------------------------------------------------------------------!
   !Calculate fraction of open canopy                                                      !
   !---------------------------------------------------------------------------------------!
   cpatch => csite%patch(ipa)
   can_frac = 1.d0
   do ico = 1,cpatch%ncohorts
      if(initp%solvable(ico)) then
         can_frac = can_frac                                                               &
                  * (1.d0 - min(1.d0                                                       &
                  ,dble(cpatch%nplant(ico))*dble(dbh2ca(cpatch%dbh(ico),cpatch%pft(ico)))))
      end if
   end do
   can_frac = 1.d0 - can_frac
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !   Obtaining the shedding water and t
   !---------------------------------------------------------------------------------------!
   if (any_solvable) then
      !------------------------------------------------------------------------------------!
      !     If vegetation is sufficiently abundant and not covered by snow, compute heat   !
      ! and moisture fluxes from vegetation to canopy, and flux resistance from soil or    !
      ! snow to canopy.                                                                    !
      !------------------------------------------------------------------------------------!
      c2 = max(0.d0,min(1.d0, 5.09d-1 * dble(csite%lai(ipa))))
      rd = rasgnd * (1.d0 - c2) + initp%rasveg * c2

      taii = 0.d0
      do ico = 1,cpatch%ncohorts
         taii = taii + initp%tai(ico)
      end do
      taii = 1.d0/taii

      !------------------------------------------------------------------------------------!
      !    If the canopy does not cover all of the ground, then it should not intercept    !
      ! all of the water.                                                                  !
      !------------------------------------------------------------------------------------!
      if (rk4met%pcpg > 0.d0) then
         !----- Scale interception by canopy openess (MCD 01-12-09). ----------------------!
         intercepted_tot  = rk4met%pcpg  * can_frac
         qintercepted_tot = rk4met%qpcpg * can_frac
         dintercepted_tot = rk4met%dpcpg * can_frac
         !----- Energy and mass are extensive, this guarantees conservation ---------------!
         wshed_tot    = rk4met%pcpg  - intercepted_tot
         qwshed_tot   = rk4met%qpcpg - qintercepted_tot
         dwshed_tot   = rk4met%dpcpg - dintercepted_tot
      else
         !----- No precipitation, nothing to be intercepted... ----------------------------!
         intercepted_tot  = 0.d0
         qintercepted_tot = 0.d0
         dintercepted_tot = 0.d0
         wshed_tot        = 0.d0
         qwshed_tot       = 0.d0
         dwshed_tot       = 0.d0
      end if

   else
      !------------------------------------------------------------------------------------!
      !     If the TAI is very small or total patch vegetation heat capacity is too        !
      ! small, bypass vegetation computations.  Set heat and moisture flux resistance rd   !
      ! between the "canopy" and snow or soil surface to its bare soil value. Set shed     !
      ! precipitation heat and moisture to unintercepted values.                           !
      !------------------------------------------------------------------------------------!
      rd               = rasgnd
      wshed_tot        = rk4met%pcpg
      qwshed_tot       = rk4met%qpcpg
      dwshed_tot       = rk4met%dpcpg
      intercepted_tot  = 0.d0
      qintercepted_tot = 0.d0
      dintercepted_tot = 0.d0

      !------------------------------------------------------------------------------------!
      ! Note: If the condition of low TAI for the entire patch was met, then it does not   !
      !       matter what the individual cohorts are normalized by, because they are       !
      !       effectively zero. So make sure the inverse patch TAI is a nominal non-zero/  !
      !       non-infinite number. This will only be used when parsing out intercepted     !
      !       leaf water into shed water; in which case the intercepted water is zero      !
      !       anyway. So this is just to prevent FPEs.                                     !
      !------------------------------------------------------------------------------------!
      taii = 1.d0
   end if
 
   !----- Assign conductivity coefficient -------------------------------------------------!
   rdi = rk4met%rhos / rd
  
   !---------------------------------------------------------------------------------------!
   !     Compute sensible heat and moisture fluxes between top soil or snow surface and    !
   ! canopy air.  wflxgc [kg/m2/s] is the upward vapor flux from soil or snow evaporation  !
   ! and dewgnd is the mass of dew that forms on the snow/soil surface this timestep; both !
   ! are defined as always positive or zero.                                               !
   !---------------------------------------------------------------------------------------!
   if (initp%nlev_sfcwater == 0) then
      hflxgc = cp8 * (initp%soil_tempk(nzg) - initp%can_temp) * rdi
   else
      hflxgc = cp8 * (initp%sfcwater_tempk(initp%nlev_sfcwater) - initp%can_temp) * rdi
   end if
  
   wflx  = (initp%surface_ssh - initp%can_shv) * rdi
   qwflx = wflx * (alvi8 - initp%surface_fliq * alli8)
   !---------------------------------------------------------------------------------------!


 
   !---------------------------------------------------------------------------------------!
   !     Calculate the dew flux.  Here the decision on whether dew or frost will form      !
   ! depends on the surface_temperature.                                                   !
   !---------------------------------------------------------------------------------------!
   dewgndflx  = max(0.d0, -wflx)
   qdewgndflx = dewgndflx * (alvi8 - initp%surface_fliq * alli8)
   !---------------------------------------------------------------------------------------!
   !    I know this is a lame way to define frost density, however I couldn't find a good  !
   ! parametrisation for frost over leaves (just a bunch of engineering papers on frost    !
   ! formation over flat metal surfaces), so I decided to keep it simple and stupid. I     !
   ! will leave this as the first attempt, so if you know a better way to do it, feel free !
   ! to add it here.                                                                       !
   !---------------------------------------------------------------------------------------!
   ddewgndflx = dewgndflx / (initp%surface_fliq*wdns8 + (1.d0-initp%surface_fliq)*idns8)


   !----- Temporary water/snow layers exist. ----------------------------------------------!
   if (initp%nlev_sfcwater > 0) then
      wflxgc  = max(0.d0,wflx)
      qwflxgc = max(0.d0,qwflx)
   !----- No surface water and not dry: evaporate from soil pores -------------------------!
   else if (initp%soilair01(nzg) > 0.d0 ) then
      wflxgc = max( 0.d0, (initp%ground_shv - initp%can_shv) * rdi)
      !----- Adjusting the flux accordingly to the surface fraction (no phase bias) -------!
      qwflxgc = wflxgc * ( alvi8 - initp%surface_fliq * alli8)
   !----- No surface water and really dry: don't evaporate at all -------------------------!
   else 
      wflxgc  = 0.d0
      qwflxgc = 0.d0
   endif
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Loop over the cohorts in the patch. Calculate energy fluxes with surrounding      !
   ! canopy air space, integrate cohort energy, calculate precipitation throughfall and    !
   ! sum fluxes to the patch level. Initialize variables used to store sums over cohorts.  !
   !---------------------------------------------------------------------------------------!
   hflxvc_tot  = 0.d0               
   wflxvc_tot  = 0.d0               
   qwflxvc_tot = 0.d0               
   cflxvc_tot  = csite%cwd_rh(ipa) 
   transp_tot  = 0.d0               
   cflxgc      = csite%rh(ipa) - csite%cwd_rh(ipa)
  
   cohortloop: do ico = 1,cpatch%ncohorts
      
      cflxgc = cflxgc + cpatch%root_respiration(ico)
      
      !------------------------------------------------------------------------------------!
      !    Calculate 'decay' term of storage (same for all) need to convert units from     !
      ! kgC/plant/day to umolC/m2/s.                                                       !
      !------------------------------------------------------------------------------------!
      storage_decay = ( dble(cpatch%growth_respiration(ico))                               &
                      + dble(cpatch%storage_respiration(ico))                              &
                      + dble(cpatch%vleaf_respiration(ico)))                               &
                    * dble(cpatch%nplant(ico)) / (day_sec8 * umol_2_kgC8)
      cflxvc_tot    = cflxvc_tot + storage_decay
      
      
      !------------------------------------------------------------------------------------!
      !     Check whether this this cohort hasn't been flagged as non-solvable, i.e., it   !
      ! has leaves, belongs to a patch that is not too sparse, an it is not buried in      !
      ! snow.  We should compute energy and water at the cohort level only if the cohort   !
      ! is "safe".  Otherwise, we will set the leaf energy derivatives to zero, and pass   !
      ! all throughfall to the ground.  Later, these "unsafe" cohorts will have their leaf !
      ! energy set to equilibrium with the canopy air space (temperature).                 !
      !------------------------------------------------------------------------------------!
      if (initp%solvable(ico)) then

         !------ Defining some shortcuts to indices ---------------------------------------!
         ipft  = cpatch%pft(ico)
         kroot = cpatch%krdepth(ico)


         !------  Calculate leaf-level CO2 flux -------------------------------------------!
         leaf_flux = dble(cpatch%gpp(ico)) - dble(cpatch%leaf_respiration(ico))

         !------ Update CO2 flux from vegetation to canopy air space. ---------------------!
         cflxvc_tot = cflxvc_tot - leaf_flux

         !---------------------------------------------------------------------------------!
         !     Defining the minimum leaf water to be considered, and the maximum amount    !
         ! possible. Here we use TAI because the intercepted area should be proportional   !
         ! to the projected area.                                                          !
         !---------------------------------------------------------------------------------!
         min_leaf_water = rk4dry_veg_lwater*initp%tai(ico)
         max_leaf_water = rk4fullveg_lwater*initp%tai(ico)
         
         !------ Calculate fraction of leaves covered with water. -------------------------!
         if(initp%veg_water(ico) > min_leaf_water)then
            sigmaw = min(1.d0, (initp%veg_water(ico)/max_leaf_water)**twothirds8)
         else
            sigmaw = 0.d0
         end if


         !---------------------------------------------------------------------------------!
         !     Finding the saturation specific humidity associated with leaf temperature.  !
         ! The minimum is set to one to avoid FPE errors, the step will be rejected should !
         ! this happen.                                                                    !
         !---------------------------------------------------------------------------------!
         veg_temp_sat = sngl(max(toocold,initp%veg_temp(ico)))
         sat_shv      = dble(rhovsil(veg_temp_sat)) / rk4met%rhos

         !---------------------------------------------------------------------------------!
         !    Here we must compute two different areas.  For transpiration, we want the    !
         ! leaf area index only, because we assume transpiration to happen only through    !
         ! leaves.  Evaporation of water/ice settled over the vegetation surface, or dew/  !
         ! frost formation must account the branches and stems as well.                    !
         !---------------------------------------------------------------------------------!
         !----- Transpiration "flux" ------------------------------------------------------!
         flux_area = initp%lai(ico)
         c3lai  = sngloff( flux_area * rk4met%rhos * (sat_shv - initp%can_shv),tiny_offset)
         !----- Evaporation/condensation "flux" -------------------------------------------!
         flux_area = effarea_water * initp%lai(ico) + pi18 * initp%wpa(ico)
         c3tai  = flux_area * rk4met%rhos * (sat_shv - initp%can_shv)
         rbi    = 1.d0 / initp%rb(ico)


         !---------------------------------------------------------------------------------!
         !    Computing the evapotranspiration or dew/frost deposition.                    !
         !---------------------------------------------------------------------------------!
         if (c3tai >= 0.d0) then  
            !------------------------------------------------------------------------------!
            !    Evapotranspiration                                                        !
            !------------------------------------------------------------------------------!
            !----- Evaporation, energy is scaled by liquid/ice partition (no phase bias). -!
            wflxvc  = c3tai * sigmaw * rbi
            qwflxvc = wflxvc * (alvi8 - initp%veg_fliq(ico) * alli8)
            !----- Transpiration, consider the leaf area rather than TAI. -----------------!
            if (initp%solvable(ico) .and. initp%available_liquid_water(kroot) > 0.d0 ) then
               cpatch%Psi_open(ico)   = c3lai / (cpatch%rb(ico) + cpatch%rsw_open(ico)  )
               cpatch%Psi_closed(ico) = c3lai / (cpatch%rb(ico) + cpatch%rsw_closed(ico))
               transp = dble(cpatch%fs_open(ico)) * dble(cpatch%Psi_open(ico))             &
                      + (1.0d0 - dble(cpatch%fs_open(ico))) * dble(cpatch%Psi_closed(ico))
           else
              cpatch%Psi_open(ico) = 0.
              cpatch%Psi_closed(ico) = 0.
              transp = 0.0d0
            end if
            qtransp = transp * alvl8
         else
            !------ Dew/frost formation ---------------------------------------------------!
            wflxvc                 = c3tai * rbi
            qwflxvc                = wflxvc  * (alvi8 - initp%veg_fliq(ico)*alli8)
            transp                 = 0.0d0
            qtransp                = 0.0d0
            cpatch%Psi_open(ico)   = 0.0d0
            cpatch%Psi_closed(ico) = 0.0d0
         end if

         !----- Diagnostic ----------------------------------------------------------------!
         dinitp%ebudget_latent = dinitp%ebudget_latent + qwflxvc + qtransp


         !----- We need to extract water from the soil equal to the transpiration. --------!
         initp%extracted_water(kroot) = initp%extracted_water(kroot) + transp


         !---------------------------------------------------------------------------------!
         !   Calculate vegetation-to-canopy sensible heat flux.  Consider the two sides of !
         ! leaves plus the actual projected branch area (not the effective), thus the pi   !
         ! factor (which to make it scalable with the cilinder.                            !
         !---------------------------------------------------------------------------------!
         flux_area = effarea_heat * initp%lai(ico) + pi18 * initp%wpa(ico)
         hflxvc    = flux_area * cp8 * rk4met%rhos * rbi                                   &
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
            wshed                 = intercepted_tot  * initp%tai(ico) * taii
            qwshed                = qintercepted_tot * initp%tai(ico) * taii
            dwshed                = dintercepted_tot * initp%tai(ico) * taii

            intercepted           = 0.d0
            qintercepted          = 0.d0
         else
            !------------------------------------------------------------------------------!
            ! Case 2: Leaf has space for rain. Rainfall and its internal energy accumulate !
            !         on the leaf.                                                         !
            !------------------------------------------------------------------------------!
            wshed                 = 0.d0
            qwshed                = 0.d0
            dwshed                = 0.d0
            intercepted           = intercepted_tot  * initp%tai(ico) * taii
            qintercepted          = qintercepted_tot * initp%tai(ico) * taii
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
         !     If there is not enough biomass to safely solve the vegetation energy        !
         ! balance, leaf fluxes and interception are set to zero.                          !
         !---------------------------------------------------------------------------------!
         dinitp%veg_energy(ico) = 0.d0
         dinitp%veg_water(ico)  = 0.d0

         !---------------------------------------------------------------------------------!
         !     Allow the complete bypass of precipitation if there are no leaves.          !
         !---------------------------------------------------------------------------------!
         wshed_tot  = wshed_tot  + intercepted_tot  * initp%tai(ico) * taii
         qwshed_tot = qwshed_tot + qintercepted_tot * initp%tai(ico) * taii
         dwshed_tot = dwshed_tot + dintercepted_tot * initp%tai(ico) * taii
      end if
   end do cohortloop


   !---------------------------------------------------------------------------------------!
   !     Update temperature and moisture of canopy.  hcapcan [J/m2/K] and wcapcan          !
   ! [kg_air/m2] are the heat and moisture capacities of the canopy.                       !
   !---------------------------------------------------------------------------------------!
   dinitp%can_temp = (hflxgc + hflxvc_tot + hflxac) * hcapcani
   dinitp%can_shv  = (wflxgc - dewgndflx + wflxvc_tot + transp_tot +  wflxac) * wcapcani


   !----- Update CO2 concentration in the canopy ------------------------------------------!
   dinitp%can_co2  = (cflxgc + cflxvc_tot + cflxac) * ccapcani


   !---------------------------------------------------------------------------------------!
   !     Integrate diagnostic variables - These are not activated unless fast file-type    !
   ! outputs are selected. This will speed up the integrator.                              !
   !---------------------------------------------------------------------------------------!
   if (fast_diagnostics) then

      dinitp%wbudget_loss2atm = - wflxac
      dinitp%ebudget_loss2atm = - hflxac
      dinitp%ebudget_latent   = dinitp%ebudget_latent -qdewgndflx + qwflxgc
      
      dinitp%co2budget_loss2atm = - cflxac
      dinitp%avg_carbon_ac      =   cflxac

      dinitp%avg_sensible_vc   = hflxvc_tot                     ! Sens. heat,  Leaf->Canopy
      dinitp%avg_sensible_2cas = hflxgc+hflxac+hflxvc_tot       ! Sens. heat,  All ->Canopy
      dinitp%avg_vapor_vc      = wflxvc_tot                     ! Lat.  heat,  Leaf->Canopy
      dinitp%avg_sensible_gc   = hflxgc                         ! Sens. heat,  Gnd ->Canopy
      dinitp%avg_sensible_ac   = hflxac / rk4met%exner          ! Sens. heat,  Atmo->Canopy
      dinitp%avg_vapor_ac      = wflxac                         ! Lat.  heat,  Atmo->Canopy
      dinitp%avg_transp        = transp_tot                     ! Transpiration
      dinitp%avg_evap          = qwflxgc-qdewgndflx+qwflxvc_tot ! Evaporation/Condensation
      dinitp%avg_sensible_tot  = (hflxgc + hflxvc_tot)          ! Total Sensible heat
      dinitp%avg_netrad = dble(csite%rlong_g(ipa)) + dble(csite%rlong_s(ipa))              &
                        + dble(csite%rshort_g(ipa))
      do k=1,initp%nlev_sfcwater
         dinitp%avg_netrad = dinitp%avg_netrad + dble(csite%rshort_s(k,ipa))
      end do
   end if
  
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
   if(debug .and. abs(rk4met%atm_tmp) < tiny(1.d0)) print*,"atm_tmp = 0"
   gzotheta = grav8 * rk4met%geoht * cpi8 * rk4met%exner / rk4met%atm_tmp
   dinitp%wpwp = vertical_vel_flux8(gzotheta,initp%tstar,initp%ustar)


   !---------------------------------------------------------------------------------------!
   !     If the single pond layer is too thin, force equilibrium with top soil layer.      !
   ! This is done in two steps: first, we don't transfer all the energy to the top soil    !
   ! layer.  Then, in redistribute_snow, we will make them in thermal equilibrium.      !
   !---------------------------------------------------------------------------------------!
   if (initp%nlev_sfcwater == 1) then
      if (abs(initp%sfcwater_mass(1)) < rk4water_stab_thresh) then
         dinitp%soil_energy(nzg)   = dinitp%soil_energy(nzg)                               &
                                   + dinitp%sfcwater_energy(1) * dslzi8(nzg)
         dinitp%sfcwater_energy(1) = 0.d0
      end if
   end if

   return
end subroutine canopy_derivs_two
!==========================================================================================!
!==========================================================================================!

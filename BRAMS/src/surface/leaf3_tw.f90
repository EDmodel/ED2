!==========================================================================================!
!==========================================================================================!
!      This sub-routine is just an interface for the main solver.  Some compilers like     !
! ifort 10 give segmentation violation if we call the subroutine directly from leaf3.f90.  !
!------------------------------------------------------------------------------------------!
subroutine leaf_prognostic(ifm,i,j,ip)
   use mem_leaf   , only : leaf_g     ! ! intent(inout)
   use mem_radiate, only : radiate_g  ! ! intent(inout)
   use mem_grid   , only : nzg        & ! intent(in)
                         , nzs        ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                , intent(in)    :: ifm
   integer                , intent(in)    :: i
   integer                , intent(in)    :: j
   integer                , intent(in)    :: ip
   !---------------------------------------------------------------------------------------!

   call leaftw(nzg,nzs                                                                     &
              ,leaf_g(ifm)%soil_water     (:,i,j,ip),leaf_g(ifm)%soil_energy    (:,i,j,ip) &
              ,leaf_g(ifm)%soil_text      (:,i,j,ip),leaf_g(ifm)%sfcwater_energy(:,i,j,ip) &
              ,leaf_g(ifm)%sfcwater_mass  (:,i,j,ip),leaf_g(ifm)%sfcwater_depth (:,i,j,ip) &
              ,leaf_g(ifm)%ustar            (i,j,ip),leaf_g(ifm)%tstar            (i,j,ip) &
              ,leaf_g(ifm)%rstar            (i,j,ip),leaf_g(ifm)%cstar            (i,j,ip) &
              ,leaf_g(ifm)%zeta             (i,j,ip),leaf_g(ifm)%ribulk           (i,j,ip) &
              ,leaf_g(ifm)%veg_albedo       (i,j,ip),leaf_g(ifm)%veg_fracarea     (i,j,ip) &
              ,leaf_g(ifm)%veg_lai          (i,j,ip),leaf_g(ifm)%veg_tai          (i,j,ip) &
              ,leaf_g(ifm)%veg_rough        (i,j,ip),leaf_g(ifm)%veg_height       (i,j,ip) &
              ,leaf_g(ifm)%veg_displace     (i,j,ip),leaf_g(ifm)%patch_area       (i,j,ip) &
              ,leaf_g(ifm)%patch_rough      (i,j,ip),leaf_g(ifm)%patch_wetind     (i,j,ip) &
              ,leaf_g(ifm)%leaf_class       (i,j,ip),leaf_g(ifm)%soil_rough       (i,j,ip) &
              ,leaf_g(ifm)%sfcwater_nlev    (i,j,ip),leaf_g(ifm)%stom_condct      (i,j,ip) &
              ,leaf_g(ifm)%ground_rsat      (i,j,ip),leaf_g(ifm)%ground_rvap      (i,j,ip) &
              ,leaf_g(ifm)%ground_temp      (i,j,ip),leaf_g(ifm)%ground_fliq      (i,j,ip) &
              ,leaf_g(ifm)%veg_water        (i,j,ip),leaf_g(ifm)%veg_hcap         (i,j,ip) &
              ,leaf_g(ifm)%veg_energy       (i,j,ip),leaf_g(ifm)%can_prss         (i,j,ip) &
              ,leaf_g(ifm)%can_theiv        (i,j,ip),leaf_g(ifm)%can_theta        (i,j,ip) &
              ,leaf_g(ifm)%can_rvap         (i,j,ip),leaf_g(ifm)%can_co2          (i,j,ip) &
              ,leaf_g(ifm)%sensible_gc      (i,j,ip),leaf_g(ifm)%sensible_vc      (i,j,ip) &
              ,leaf_g(ifm)%evap_gc          (i,j,ip),leaf_g(ifm)%evap_vc          (i,j,ip) &
              ,leaf_g(ifm)%transp           (i,j,ip),leaf_g(ifm)%gpp              (i,j,ip) &
              ,leaf_g(ifm)%plresp           (i,j,ip),leaf_g(ifm)%resphet          (i,j,ip) &
              ,leaf_g(ifm)%veg_ndvip        (i,j,ip),leaf_g(ifm)%veg_ndvic        (i,j,ip) &
              ,leaf_g(ifm)%veg_ndvif        (i,j,ip),radiate_g(ifm)%rshort        (i,j)    &
              ,radiate_g(ifm)%cosz          (i,j)   )

   return
end subroutine leaf_prognostic
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine leaftw(mzg,mzs,soil_water, soil_energy,soil_text,sfcwater_energy_int            &
                 ,sfcwater_mass,sfcwater_depth,ustar,tstar,rstar,cstar,zeta,ribulk         &
                 ,veg_albedo,veg_fracarea,veg_lai,veg_tai,veg_rough,veg_height             &
                 ,veg_displace,patch_area,patch_rough,patch_wetind,leaf_class,soil_rough   &
                 ,sfcwater_nlev,stom_condct,ground_rsat,ground_rvap,ground_temp            &
                 ,ground_fliq,veg_water,veg_hcap,veg_energy,can_prss,can_theiv,can_theta   &
                 ,can_rvap,can_co2,sensible_gc,sensible_vc,evap_gc,evap_vc,transp,gpp      &
                 ,plresp,resphet,veg_ndvip,veg_ndvic,veg_ndvif,rshort,cosz)

   use leaf_coms
   use mem_leaf
   use rconstants
   use mem_scratch
   use therm_lib   , only : qwtk        & ! subroutine
                          , qtk         & ! subroutine
                          , hpqz2temp   & ! function
                          , idealdenssh ! ! function
   use catt_start  , only : CATT        ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                , intent(in)    :: mzg
   integer                , intent(in)    :: mzs
   real   , dimension(mzg), intent(inout) :: soil_water
   real   , dimension(mzg), intent(inout) :: soil_energy
   real   , dimension(mzg), intent(inout) :: soil_text
   real   , dimension(mzs), intent(inout) :: sfcwater_energy_int
   real   , dimension(mzs), intent(inout) :: sfcwater_mass
   real   , dimension(mzs), intent(inout) :: sfcwater_depth
   real                   , intent(in)    :: ustar
   real                   , intent(in)    :: tstar
   real                   , intent(in)    :: rstar
   real                   , intent(in)    :: cstar
   real                   , intent(in)    :: zeta
   real                   , intent(in)    :: ribulk
   real                   , intent(in)    :: veg_albedo
   real                   , intent(in)    :: veg_fracarea
   real                   , intent(in)    :: veg_lai
   real                   , intent(in)    :: veg_tai
   real                   , intent(in)    :: veg_rough
   real                   , intent(in)    :: veg_height
   real                   , intent(in)    :: veg_displace
   real                   , intent(in)    :: patch_area
   real                   , intent(in)    :: patch_rough
   real                   , intent(in)    :: patch_wetind
   real                   , intent(in)    :: leaf_class
   real                   , intent(in)    :: soil_rough
   real                   , intent(inout) :: sfcwater_nlev
   real                   , intent(inout) :: stom_condct
   real                   , intent(inout) :: ground_rsat
   real                   , intent(inout) :: ground_rvap
   real                   , intent(inout) :: ground_temp
   real                   , intent(inout) :: ground_fliq
   real                   , intent(inout) :: veg_water
   real                   , intent(inout) :: veg_hcap
   real                   , intent(inout) :: veg_energy
   real                   , intent(inout) :: can_prss
   real                   , intent(inout) :: can_theiv
   real                   , intent(inout) :: can_theta
   real                   , intent(inout) :: can_rvap
   real                   , intent(inout) :: can_co2
   real                   , intent(inout) :: sensible_gc
   real                   , intent(inout) :: sensible_vc
   real                   , intent(inout) :: evap_gc
   real                   , intent(inout) :: evap_vc
   real                   , intent(inout) :: transp
   real                   , intent(inout) :: gpp
   real                   , intent(inout) :: plresp
   real                   , intent(inout) :: resphet
   real                   , intent(in)    :: veg_ndvip
   real                   , intent(in)    :: veg_ndvic
   real                   , intent(in)    :: veg_ndvif
   real                   , intent(in)    :: rshort
   real                   , intent(in)    :: cosz
   !----- Local variables. ----------------------------------------------------------------!
   integer                                :: k
   integer                                :: nveg
   integer                                :: ksn
   integer                                :: nsoil
   integer                                :: ksnnew
   integer                                :: newlayers
   integer                                :: nlayers
   integer                                :: kold
   integer                                :: ktrans
   integer                                :: kzs
   real                                   :: sfcw_temp
   real                                   :: sfcw_fliq
   real                                   :: stretch
   real                                   :: thik
   real                                   :: snden
   real                                   :: wfree
   real                                   :: soilhcap
   real                                   :: depthloss
   real                                   :: soilcap
   real                                   :: sndenmax
   real                                   :: sndenmin
   real                                   :: hold
   real                                   :: wtnew
   real                                   :: wtold
   real                                   :: wdiff
   real                                   :: wgpmid
   real                                   :: freezecor
   real                                   :: psiplusz_bl
   real                                   :: available_water
   real   , dimension(mzg)                :: available_layer
   real                                   :: extracted_water
   real                                   :: wloss
   real                                   :: qwloss
   real                                   :: soilcond
   real                                   :: wgpfrac
   !---------------------------------------------------------------------------------------!


   !----- Initialise some levels. ---------------------------------------------------------!
   do k = 1,mzg
      dslz   (k) = slz(k+1) - slz(k)
      dslzi  (k) = 1. / dslz(k)
      dslzidt(k) = dslzi(k) * dtll
      slzt   (k) = .5 * (slz(k) + slz(k+1))
   end do
   do k = 2,mzg
      dslzt   (k) = slzt(k) - slzt(k-1)
      dslzti  (k) = 1. / dslzt(k)
      dslztidt(k) = dslzti(k) * dtll
   end do

   !----- These must be defined for free drainage bc (RGK) --------------------------------!
   dslzt    (1) = 2.0*slz(1) - slzt(1)
   dslzti   (1) = 1./dslzt(1)
   dslztidt (1) = dslzti(1) * dtll


   nveg = nint(leaf_class)
   ksn = nint(sfcwater_nlev)



   !---------------------------------------------------------------------------------------!
   !    Remove water from all layers up to the root depth for transpiration.  Bottom       !
   ! k-index in root zone is kroot(nveg).  Limit soil moisture to values above soilwp.     !
   ! The available water profile is given in kg/m2.                                        !
   !---------------------------------------------------------------------------------------!
   ktrans = kroot(nveg)
   available_layer(:)   = 0.
   available_water      = 0.
   do k = mzg,ktrans,-1
      nsoil              = nint(soil_text(k))
      available_layer(k) = dslz(k) * wdns                                                  &
                         * max(0.0, soil_fracliq(k) * (soil_water(k)-soilwp(nsoil)))
      available_water    = available_water + available_layer(k)
   end do
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Compute gravitational potential plus moisture potential, psi + z (psiplusz) [m],  !
   ! liquid water content (soil_liq) [m], and 99% the remaining water capacity (soilair99) !
   ! [m].                                                                                  !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Compute gravitational potential plus moisture potential z + psi [m], liquid water !
   ! content [m], and half the remaining water capacity [m].                               !
   !---------------------------------------------------------------------------------------!
   do k = 1,mzg
      nsoil           = nint(soil_text(k))
      psiplusz(k)     = slzt(k) + slpots(nsoil)*(slmsts(nsoil)/soil_water(k))**slbs(nsoil)
      soil_liq(k)     = dslz(k) * max(0.,(soil_water(k) - soilwp(nsoil)) * soil_fracliq(k))
      soilair99(k)    = 0.99 * slmsts(nsoil) - soil_water(k)
      soilair01(k)    = (soil_water(k) - 1.01 * soilcp(nsoil)) * soil_fracliq(k)
      half_soilair(k) = (slmsts(nsoil) - soil_water(k)) * dslz(k) * .5
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Evaluate any exchanges of heat and moisture to or from vegetation, apply moist-  !
   ! ure and heat changes to vegetation, and evaluate the resistance parameter rgnd        !
   ! between canopy air and the top soil or snow surface.                                  !
   !---------------------------------------------------------------------------------------!
   call leaf_canopy(mzg,mzs,ksn,soil_energy,soil_water,soil_text,sfcwater_mass,ustar       &
                   ,tstar,rstar,cstar,zeta,ribulk,soil_rough,veg_rough,patch_rough         &
                   ,veg_height,veg_displace,veg_lai,veg_tai,veg_water,veg_hcap,veg_energy  &
                   ,leaf_class,veg_fracarea,stom_condct,can_prss,can_rvap,can_co2          &
                   ,sensible_gc,sensible_vc,evap_gc,evap_vc,transp,gpp,plresp,resphet      &
                   ,ground_rsat,ground_rvap,ground_temp,ground_fliq,available_water,rshort)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Compute soil and temporary surface water/snow heat resistance times layer depth   !
   ! (the rfactor variable).                                                               !
   !---------------------------------------------------------------------------------------!
   do k = 1,mzg
      nsoil     = nint(soil_text(k))
      wgpfrac   = min(1.0, soil_water(k) / slmsts(nsoil))
      soilcond  = soilcond0(nsoil)                                                         &
                + wgpfrac * (soilcond1(nsoil) + wgpfrac * soilcond2(nsoil))
      rfactor(k) = dslz(k) / soilcond
   end do
   do k = 1,ksn
      if (sfcwater_depth(k) > 0.0) then
         snden = sfcwater_mass(k) / sfcwater_depth(k)
         rfactor(k+mzg) = sfcwater_depth(k)                                                &
                        / ( ss(1) * exp(ss(2) * sfcwater_tempk(k))                         &
                          * ( ss(3) + snden * (ss(4) + snden * (ss(5) + snden * ss(6)) ) ))
      else
         rfactor(k+mzg) = 0.0
      end if
   end do
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !     Find the sensible heat fluxes find soil and surface water internal sensible heat  !
   ! fluxes (hfluxgsc) [W/m2].                                                             !
   !---------------------------------------------------------------------------------------!
   hfluxgsc(:) = 0.
   do k = 2,mzg
      hfluxgsc(k) = - (soil_tempk(k) - soil_tempk(k-1)) / ((rfactor(k) + rfactor(k-1)) * .5)
   end do
   !----- If temporary water/snow layers exist, compute them now... -----------------------!
   if (ksn >= 1) then
      hfluxgsc(mzg+1) = - (sfcwater_tempk(1) - soil_tempk(mzg))                            &
                        / ((rfactor(mzg+1)   + rfactor(mzg)) * 0.5)

      do k = 2,ksn
         hfluxgsc(mzg+k) = - (sfcwater_tempk(k) - sfcwater_tempk(k-1))                     &
                         /   ((rfactor(mzg+k) + rfactor(mzg+k-1)) * 0.5)
      end do
   end if
   !---------------------------------------------------------------------------------------!
   !     Heat flux at top soil layer or top temporary surface water/snow layer from long   !
   ! wave, sensible, and upward latent heat fluxes [W/m^2].                                !
   !---------------------------------------------------------------------------------------!
   hfluxgsc(ksn+1+mzg) = hflxgc + qwflxgc - rlonga_gs - rlongv_gs + rlonggs_v + rlonggs_a
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Update soil internal energy values [J/m3] and snow/surface water internal energy !
   ! values [J/m2] from sensible heat, upward water vapor (latent heat), longwave, and     !
   ! shortwave fluxes.  This excludes effects of dew/frost formation, precipitation, shed- !
   ! ding, and percolation.  Update top soil or snow moisture from evaporation only.       !
   !---------------------------------------------------------------------------------------!
   do k = 1,mzg
      soil_energy(k) = soil_energy(k) + dslzidt(k) * (hfluxgsc(k) - hfluxgsc(k+1))
   end do
   soil_energy(mzg)  = soil_energy(mzg) + dslzidt(mzg) * rshort_g

   do k = 1,ksn
      sfcwater_energy_ext(k) = sfcwater_energy_ext(k)                                      &
                             + dtll * (hfluxgsc(k+mzg) - hfluxgsc(k+1+mzg) + rshort_s(k))
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Reset all soil an temporar surface water fluxes.                                  !
   !---------------------------------------------------------------------------------------!
   w_flux (:) = 0.
   qw_flux(:) = 0.
   d_flux(:)  = 0.
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Add water from above to the top temporary surface water (snow) layer if there is  !
   ! one, otherwise dump it into the "virtual" layer.  It will be eliminated soon.         !
   !---------------------------------------------------------------------------------------!
   if (ksn > 0) then
      sfcwater_mass(ksn)       = sfcwater_mass(ksn)                                        &
                               + dewgnd_tot + wshed_tot + throughfall_tot - wflxgc * dtll
      sfcwater_energy_ext(ksn) = sfcwater_energy_ext(ksn)                                  &
                               + qdewgnd_tot + qwshed_tot + qthroughfall_tot
      sfcwater_depth(ksn)      = sfcwater_depth(ksn)                                       &
                               + ddewgnd_tot + dwshed_tot + dthroughfall_tot
      w_flux(mzg+1)            = 0.0
   else
      virtual_water  = virtual_water  +  dewgnd_tot +  wshed_tot +  throughfall_tot
      virtual_energy = virtual_energy + qdewgnd_tot + qwshed_tot + qthroughfall_tot
      virtual_depth  = virtual_depth  + ddewgnd_tot + dwshed_tot + dthroughfall_tot
      !----- w_flux should be in m. -------------------------------------------------------!
      w_flux(mzg+1)  = wflxgc * dtll * wdnsi
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find the water flux in each of the soil layers, and the corresponding internal    !
   ! energy flux.                                                                          !
   !---------------------------------------------------------------------------------------!
   do k = 2, mzg
      nsoil = nint(soil_text(k))

      wgpmid    = 0.5 * (soil_water(k)   + soil_water(k-1))
      freezecor = 0.5 * (soil_fracliq(k) + soil_fracliq(k-1))
      if (freezecor < 1.0) freezecor = 10.**(- freezecoef * (1.0 - freezecor))
      w_flux(k) = dslztidt(k) * slcons1(k,nsoil)                                           &
                * (wgpmid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)                      &
                * (psiplusz(k-1) - psiplusz(k)) * freezecor

      !------------------------------------------------------------------------------------!
      !     Limit water transfers to prevent over-saturation and over-depletion.  Compute  !
      ! q transfers between soil layers (qw_flux) [J/m2].                                  !
      !------------------------------------------------------------------------------------!
      if (w_flux(k) > 0.) then
         w_flux(k) = min(w_flux(k),soil_liq(k-1),half_soilair(k))
      else
         w_flux(k) = - min(-w_flux(k),soil_liq(k),half_soilair(k-1))
      endif
      qw_flux(k) = w_flux(k) * cliqvlme * (soil_tempk(k) - tsupercool)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Boundary condition at the lowest soil level.                                     !
   !---------------------------------------------------------------------------------------!
   nsoil = nint(soil_text(1))
   select case (isoilbc)
   case (0)
      !----- Bedrock, no flux across. -----------------------------------------------------!
      w_flux(1)  = 0.
   case (1)
      !------------------------------------------------------------------------------------!
      !     Free drainage, water movement of bottom soil layer is only under gravity, i.e. !
      ! the soil water content of boundary layer is equal to that of bottom soil layer.    !
      !------------------------------------------------------------------------------------!

      wgpmid    = soil_water(1)
      freezecor = soil_fracliq(1) 
      if (freezecor < 1.0) freezecor = 10.**(- freezecoef * (1.0 - freezecor))
      w_flux(1) = dslztidt(1) * slcons1(1,nsoil)                                           &
                * (wgpmid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)                      &
                * (slz(2) - slz(1)) * freezecor

      !------------------------------------------------------------------------------------!
      !     Limit water transfers to prevent over-saturation and over-depletion. Compute q !
      ! transfers between soil layers (qw_flux) [J/m2].                                     !
      !------------------------------------------------------------------------------------!
      if (w_flux(1) > 0.) then
         w_flux(1) = min(w_flux(1),half_soilair(1))
      else
         w_flux(1) = - min(-w_flux(1),soil_liq(1))
      end if
      !------------------------------------------------------------------------------------!
   case (2)
      !------------------------------------------------------------------------------------!
      !     Half drainage, water movement of bottom soil layer is only under gravity, i.e. !
      ! the soil water content of boundary layer is equal to that of bottom soil layer.    !
      !------------------------------------------------------------------------------------!
      wgpmid    = soil_water(1)
      freezecor = soil_fracliq(1) 
      if (freezecor < 1.0) freezecor = 10.**(- freezecoef * (1.0 - freezecor))
      w_flux(1) = dslztidt(1) * slcons1(1,nsoil)                                           &
                * (wgpmid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)                      &
                * (slz(2) - slz(1)) * freezecor * 0.5
      !------------------------------------------------------------------------------------!
      !     Limit water transfers to prevent over-saturation and over-depletion. Compute q !
      ! transfers between soil layers (qw_flux) [J/m2].                                     !
      !------------------------------------------------------------------------------------!
      if (w_flux(1) > 0.) then
         w_flux(1) = min(w_flux(1),half_soilair(1))
      else
         w_flux(1) = - min(-w_flux(1),soil_liq(1))
      end if
      !------------------------------------------------------------------------------------!
   case (3)
      !------------------------------------------------------------------------------------!
      !     Free drainage, water movement of bottom soil layer is under gravity and moist- !
      ! ure potential difference. The soil water content of boundary layer is equal to     !
      ! field capacity.                                                                    !
      !------------------------------------------------------------------------------------!
      wgpmid = soil_water(1)
      freezecor = soil_fracliq(1) 
      if (freezecor < 1.0) freezecor = 10.**(- freezecoef * (1.0 - freezecor))
      psiplusz_bl = slz(1)                                                                 &
                  + slpots(nsoil) * (slmsts(nsoil) / sfldcap(nsoil)) ** slbs(nsoil)

      if (psiplusz_bl <= psiplusz(1)) then
          w_flux(1) =  dslztidt(1) * slcons1(1,nsoil)                                      &
                    * (wgpmid/slmsts(nsoil)) ** (2.0 * slbs(nsoil) + 3.0)                  &
                    * (psiplusz(1) - psiplusz_bl) * freezecor
      else
         !----- Prevent bottom soil layer sucking water from the boundary layer. ----------!
         w_flux(1) = 0.0
      end if
   end select
   qw_flux(1) = w_flux(1) * cliqvlme * (soil_tempk(1) - tsupercool)
   !---------------------------------------------------------------------------------------!


   !----- Update soil moisture (impose minimum value of soilcp) and q value. --------------!
   do k = 1,mzg
      soil_water(k)  = soil_water(k)  - dslzi(k) * (w_flux(k+1)  - w_flux(k))
      soil_energy(k) = soil_energy(k) - dslzi(k) * (qw_flux(k+1) - qw_flux(k))
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Update soil moisture by extracting the water lost due to transpiration.           !
   !---------------------------------------------------------------------------------------!
   extracted_water = transp_tot * dtll
   if (extracted_water > 0. .and. available_water > 0.) then
      do k = ktrans,mzg
         wloss  = wdnsi * dslzi(k) * extracted_water * available_layer(k) / available_water
         qwloss = wloss * cliqvlme * (soil_tempk(k) - tsupercool)

         soil_water(k)  = soil_water(k)  - wloss
         soil_energy(k) = soil_energy(k) - qwloss
      end do
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Re-organise the snow layers, and shed the water from the virtual layer.           !
   !---------------------------------------------------------------------------------------!
   call leaf_adjust_sfcw(mzg,mzs,soil_energy,soil_water,soil_text,sfcwater_nlev            &
                        ,sfcwater_mass,sfcwater_depth)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find whether we should get rid of some water through runoff.                      !
   !---------------------------------------------------------------------------------------!
   ksn = nint(sfcwater_nlev)
   if (ksn > 0) then
      !------------------------------------------------------------------------------------!
      !     Get rid of some liquid water from the top layer through runoff.                !
      !------------------------------------------------------------------------------------!
      if (runoff_time > 0.0 .and. sfcwater_mass(ksn) > 0.0) then
         call qtk(sfcwater_energy_ext(ksn) / sfcwater_mass(ksn)                            &
                 ,sfcwater_tempk(ksn),sfcwater_fracliq(ksn))
         wloss  = max(0., min(1.0,dtll/runoff_time)                                        &
                        * (sfcwater_mass(ksn) - min_sfcwater_mass)                         &
                        * (sfcwater_fracliq(ksn) - 0.1) / 0.9)
         qwloss = wloss * cliq * (sfcwater_tempk(ksn) - tsupercool)
         sfcwater_energy_ext(ksn) = sfcwater_energy_ext(ksn) - qwloss
         sfcwater_mass(ksn)       = sfcwater_mass(ksn)       - wloss
         sfcwater_depth(ksn)      = sfcwater_depth(ksn)      - wloss * wdnsi
      end if
      !---------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaftw
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine initialises the soil and temporary surface water variables for a    !
! given patch, before the time step iteration loop.                                        !
!------------------------------------------------------------------------------------------!
subroutine leaf_soilsfcw_diag(ip,mzg,mzs,soil_energy,soil_water,soil_text,sfcwater_nlev    &
                             ,sfcwater_energy_int,sfcwater_mass,sfcwater_depth,initial)
   use leaf_coms , only : slcpd               & ! intent(in)
                        , soil_tempk          & ! intent(out)
                        , soil_fracliq        & ! intent(out)
                        , sfcwater_energy_ext & ! intent(inout)
                        , sfcwater_tempk      & ! intent(out)
                        , sfcwater_fracliq    & ! intent(out) 
                        , virtual_energy      & ! intent(out)
                        , virtual_water       & ! intent(out) 
                        , virtual_depth       ! ! intent(out)
   use therm_lib , only : qwtk                & ! sub-routine
                        , qtk                 ! ! sub-routine
   use rconstants, only : wdns                ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                , intent(in)    :: ip
   integer                , intent(in)    :: mzg
   integer                , intent(in)    :: mzs
   real   , dimension(mzg), intent(in)    :: soil_energy
   real   , dimension(mzg), intent(in)    :: soil_water
   real   , dimension(mzg), intent(in)    :: soil_text
   real                   , intent(in)    :: sfcwater_nlev
   real   , dimension(mzs), intent(inout) :: sfcwater_energy_int
   real   , dimension(mzs), intent(inout) :: sfcwater_mass
   real   , dimension(mzs), intent(in)    :: sfcwater_depth
   logical                , intent(in)    :: initial
   !----- Local variables. ----------------------------------------------------------------!
   integer                                :: k
   integer                                :: ksn
   integer                                :: nsoil
   !---------------------------------------------------------------------------------------!



   !----- If this is a water patch, leave the sub-routine. --------------------------------!
   if (ip == 1) return
   !---------------------------------------------------------------------------------------!


   !----- Find the soil temperature and liquid water fraction. ----------------------------!
   do k = 1,mzg
      nsoil = nint(soil_text(k))
      call qwtk(soil_energy(k),soil_water(k)*wdns,slcpd(nsoil)                             &
               ,soil_tempk(k),soil_fracliq(k))
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Find the temporary surface water temperature and liquid water fraction.  Also,     !
   ! find the extensive version of the internal energy and the fraction, and the fraction  !
   ! of exposed vegetation (by exposed I mean not buried in snow).                         !
   !---------------------------------------------------------------------------------------!
   ksn     = nint(sfcwater_nlev)
   if (initial) then
      !----- Initial call, find the extensive internal energy from the intensive. ---------!
      do k=1,ksn
         sfcwater_energy_ext(k) = sfcwater_energy_int(k) * sfcwater_mass(k)
         call qtk(sfcwater_energy_int(k),sfcwater_tempk(k),sfcwater_fracliq(k))
      end do
      !----- Fill the layers above with zeroes or dummy values. ---------------------------!
      if (ksn == 0) then
         sfcwater_energy_ext(1:mzs) = 0.
         sfcwater_energy_int(1:mzs) = 0.
         sfcwater_mass      (1:mzs) = 0.
         sfcwater_tempk     (1:mzs) = soil_tempk  (mzg)
         sfcwater_fracliq   (1:mzs) = soil_fracliq(mzg)
      else
         do k=ksn+1,mzs
            sfcwater_energy_ext(k) = 0.
            sfcwater_energy_int(k) = 0.
            sfcwater_mass      (k) = 0.
            sfcwater_tempk     (k) = sfcwater_tempk  (ksn)
            sfcwater_fracliq   (k) = sfcwater_fracliq(ksn)
         end do
      end if
   else
      !----- Convert extensive internal energy into intensive. ----------------------------!     
      do k=1,ksn
         sfcwater_energy_int(k) = sfcwater_energy_ext(k) / sfcwater_mass(k)
         call qtk(sfcwater_energy_int(k),sfcwater_tempk(k),sfcwater_fracliq(k))
      end do
      !----- Fill the layers above with zeroes or dummy values. ---------------------------!
      if (ksn == 0) then
         sfcwater_energy_ext(1:mzs) = 0.
         sfcwater_energy_int(1:mzs) = 0.
         sfcwater_mass      (1:mzs) = 0.
         sfcwater_tempk     (1:mzs) = soil_tempk  (mzg)
         sfcwater_fracliq   (1:mzs) = soil_fracliq(mzg)
      else
         do k=ksn+1,mzs
            sfcwater_energy_ext(k) = 0.
            sfcwater_energy_int(k) = 0.
            sfcwater_mass      (k) = 0.
            sfcwater_tempk     (k) = sfcwater_tempk  (ksn)
            sfcwater_fracliq   (k) = sfcwater_fracliq(ksn)
         end do
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !----- Reset the virtual layer. --------------------------------------------------------!
   virtual_energy = 0.
   virtual_water  = 0.
   virtual_depth  = 0.
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf_soilsfcw_diag
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine re-arranges the snow layers, and dump the water from the virtual    !
! layer to the place that can hold the water.                                              !
!------------------------------------------------------------------------------------------!
subroutine leaf_adjust_sfcw(mzg,mzs,soil_energy,soil_water,soil_text,sfcwater_nlev         &
                           ,sfcwater_mass,sfcwater_depth)
   use mem_leaf   , only : ipercol                  ! ! intent(in)
   use leaf_coms  , only : sfcwater_energy_ext      & ! intent(in)
                         , min_sfcwater_mass        & ! intent(in)
                         , min_sfcwater_depth       & ! intent(in)
                         , water_stab_thresh        & ! intent(in)
                         , virtual_water            & ! intent(inout)
                         , virtual_energy           & ! intent(inout)
                         , virtual_depth            & ! intent(inout)
                         , dslz                     & ! intent(in)
                         , dslzi                    & ! intent(in)
                         , slcpd                    & ! intent(in)
                         , slmsts                   & ! intent(in)
                         , thick                    & ! intent(in)
                         , thicknet                 ! ! intent(in)
   use mem_scratch, only : newsfcw_mass   => vctr14 & ! intent(out)
                         , newsfcw_energy => vctr16 & ! intent(out)
                         , newsfcw_depth  => vctr18 ! ! intent(out)
   use rconstants , only : wdns                     & ! intent(in)
                         , wdnsi                    & ! intent(in)
                         , cliq                     & ! intent(in)
                         , cice                     & ! intent(in)
                         , qliqt3                   & ! intent(in)
                         , tsupercool               ! ! intent(in)
   use therm_lib  , only : qtk                      & ! sub-routine
                         , qwtk                     ! ! sub-routine
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: mzg
   integer, intent(in) :: mzs
   real   , dimension(mzg), intent(inout) :: soil_energy
   real   , dimension(mzg), intent(inout) :: soil_water
   real   , dimension(mzg), intent(in)    :: soil_text
   real                   , intent(inout) :: sfcwater_nlev
   real   , dimension(mzs), intent(inout) :: sfcwater_mass
   real   , dimension(mzs), intent(inout) :: sfcwater_depth
   !----- Local variables. ----------------------------------------------------------------!
   integer                                :: kold
   integer                                :: newlayers
   integer                                :: nlayers
   integer                                :: ksn
   integer                                :: ksnnew
   integer                                :: k
   integer                                :: nsoil
   real                                   :: wtold
   real                                   :: wtnew
   real                                   :: wdiff
   real                                   :: sum_sfcw_mass
   real                                   :: sum_sfcw_energy
   real                                   :: sum_sfcw_depth
   real                                   :: energy_free
   real                                   :: wmass_free
   real                                   :: depth_free
   real                                   :: wmass_perc
   real                                   :: energy_perc
   real                                   :: depth_perc
   real                                   :: i_energy_try
   real                                   :: energy_try
   real                                   :: wmass_try
   real                                   :: depth_try
   real                                   :: temp_try
   real                                   :: fliq_try
   real                                   :: energy_tot
   real                                   :: wmass_tot
   real                                   :: hcapdry_tot
   real                                   :: wmass_room
   real                                   :: depthloss
   real                                   :: snden
   real                                   :: sndenmin
   real                                   :: sndenmax
   real                                   :: cr               ! snow water holding capacity
   real                                   :: gi               ! Partial density of ice
   !----- Constants -----------------------------------------------------------------------!
   real                   , parameter     :: crmin   = 0.03
   real                   , parameter     :: crmax   = 0.1
   real                   , parameter     :: ge      = 200.
   !---------------------------------------------------------------------------------------!


   !----- Copy the original number of temporary surface water layers to ksn. --------------!
   ksn       = nint(sfcwater_nlev)
   !---------------------------------------------------------------------------------------!



   !----- Copy the soil type at the topmost level to nsoil. -------------------------------!
   nsoil     = nint(soil_text(mzg))
   !---------------------------------------------------------------------------------------!
   


   !---------------------------------------------------------------------------------------!
   !      Determine the total amount of temporary surface water available as well as       !
   ! derived properties.                                                                   !
   !---------------------------------------------------------------------------------------!
   sum_sfcw_mass   = 0.0
   sum_sfcw_energy = 0.0
   sum_sfcw_depth  = 0.0
   do k=1,ksn
      sum_sfcw_mass   = sum_sfcw_mass   + sfcwater_mass      (k)
      sum_sfcw_energy = sum_sfcw_energy + sfcwater_energy_ext(k)
      sum_sfcw_depth  = sum_sfcw_depth  + sfcwater_depth     (k)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Check the total amount of water that has just fallen plus the amount that is al-  !
   ! ready sitting over the top soil layer.  We must do this as the first step because we  !
   ! may want to eliminate this water by binding it to the top soil layer in case there is !
   ! too little water.                                                                     !
   !---------------------------------------------------------------------------------------!
   if ((virtual_water + sum_sfcw_mass) < min_sfcwater_mass) then
      !------------------------------------------------------------------------------------!
      !     The mass of the potential new temporary surface water is within bounds but it  !
      ! is too small to be maintained.  We add both the virtual mass and the total surface !
      ! water and dump in the free water, but set ksnnew to zero so all the water is       !
      ! infiltrated in the top soil layer.                                                 !
      !------------------------------------------------------------------------------------!
      wmass_free             = virtual_water  + sum_sfcw_mass
      energy_free            = virtual_energy + sum_sfcw_energy
      depth_free             = virtual_depth  + sum_sfcw_depth
      !----- Reset both the temporary surface water and the virtual layer. ----------------!
      virtual_water          = 0.0
      virtual_energy         = 0.0
      virtual_depth          = 0.0
      sfcwater_mass(:)       = 0.0
      sfcwater_energy_ext(:) = 0.0
      sfcwater_depth(:)      = 0.0
      !----- Set ksnnew to zero to force all free water to go to the soil. ----------------!
      ksnnew                 = 0
   else
      !------------------------------------------------------------------------------------!
      !     The mass of the potential new temporary surface water could create at least    !
      ! one layer.  If there is already a temporary surface water or snow layer, the new   !
      ! amount is initially put there, otherwise, we attempt to create the first layer.    !
      !------------------------------------------------------------------------------------!
      wmass_free             = virtual_water
      energy_free            = virtual_energy
      depth_free             = virtual_depth
      virtual_water          = 0.0
      virtual_energy         = 0.0
      virtual_depth          = 0.0
      ksnnew                 = max(ksn,1)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Update the prognostic and diagnostic variables by adding the free standing water.  !
   ! Then we check the size of the temporary surface water layers, and update the          !
   ! temperature in a way that ensure the layer stability.  During this process, we        !
   ! update the total temporary surface water mass, energy, and depth, which will be used  !
   ! later in the sub-routine.                                                             !
   !---------------------------------------------------------------------------------------!
   sum_sfcw_mass   = 0.0
   sum_sfcw_energy = 0.0
   sum_sfcw_depth  = 0.0
   do k = ksnnew,1,-1
      !------------------------------------------------------------------------------------!
      !    Find the potential mass, energy, and depth of the temporary layer if all the    !
      ! free water became part of this layer.                                              !
      !------------------------------------------------------------------------------------!
      energy_try = sfcwater_energy_ext(k)   + energy_free
      wmass_try  = sfcwater_mass      (k)   + wmass_free
      depth_try  = sfcwater_depth     (k)   + depth_free
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    In case this is a single layer, and a very thin one, we may have a hard time    !
      ! achieving numerical stability.  We can treat this case like the leaf case, in      !
      ! the sense that the water sitting on the top of the surface is in thermal           !
      ! equilibrium with the surface.                                                      !
      !------------------------------------------------------------------------------------!
      if (ksnnew == 1 .and. wmass_try < water_stab_thresh) then
         !---------------------------------------------------------------------------------!
         !     Find the total internal energy of the combined pool (top soil layer plus    !
         ! the thin temporary surface water).  The units of soil properties are J/m3 for   !
         ! the internal energy, and m3/m3 for soil water, whilst the temporary surface     !
         ! water has units of J/m2 for internal energy and kg/m2 for mass.  We use the     !
         ! standard for the temporary surface water.                                       !
         !---------------------------------------------------------------------------------!
         energy_tot  = energy_try + soil_energy(mzg) * dslz(mzg)
         wmass_tot   = wmass_try  + soil_water(mzg)  * dslz(mzg) * wdns
         hcapdry_tot = slcpd(nsoil) * dslz(mzg)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Find the equilibrium temperature and liquid/ice partition.   Because we    !
         ! are assuming thermal equilibrium, the temperature and liquid fraction of the    !
         ! attempted layer is the same as the average temperature of the augmented pool.   !
         !---------------------------------------------------------------------------------!
         call qwtk(energy_tot,wmass_tot,hcapdry_tot,temp_try,fliq_try)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Re-compute the internal energy of the temporary layer, using the temperature !
         ! and fraction of liquid water distribution we have just found, keeping the mass  !
         ! constant.                                                                       !
         !---------------------------------------------------------------------------------!
         energy_try = wmass_try * (        fliq_try  * cliq * (temp_try - tsupercool)      &
                                  + (1.0 - fliq_try) * cice *  temp_try              )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Re-calculate the top soil internal energy, by removing the attempted surface !
         ! water energy from the total energy, and converting it back to J/m3.  The total  !
         ! amount of water does not need to be re-calculated at this time.                 !
         !---------------------------------------------------------------------------------!
         soil_energy(mzg)  = (energy_tot - energy_try) * dslzi(mzg)
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !      Layer is computationally stable, find temperature and liquid fraction of   !
         ! the attempted layer.                                                            !
         !---------------------------------------------------------------------------------!
         i_energy_try = energy_try / wmass_try
         call qtk(i_energy_try,temp_try,fliq_try)
        !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Determine a first guess for the amount of mass that can be lost from this      !
      ! layer through percolation (wmass_perc).                                            !
      !------------------------------------------------------------------------------------!
      select case (ipercol)
      case (0)
         !---------------------------------------------------------------------------------!
         !     Original method, from LEAF-3.  Shed liquid in excess of a 1:9               !
         ! liquid-to-ice ratio through percolation.                                        !
         !---------------------------------------------------------------------------------!
         wmass_perc  = max(0.0, wmass_try * (fliq_try - 0.1) / 0.9)
         !---------------------------------------------------------------------------------!
      case (1)
         !---------------------------------------------------------------------------------!
         !    Alternative "free" water calculation.                                        !
         !    Anderson (1976), NOAA Tech Report NWS 19.                                    !
         !---------------------------------------------------------------------------------!
         gi          = wmass_try/max(min_sfcwater_depth,depth_try) * (1.0 - fliq_try)
         Cr          = max(Crmin, Crmin + (Crmax - Crmin) * (ge - gi) / ge)
         wmass_perc  = max(0.0,wmass_try * (fliq_try - Cr / (1.0 + Cr)))
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Determinte whether the layer beneath the current one is another temporary      !
      ! surface water/snow layer, or the top soil layer.  In case it is the latter, we     !
      ! must check whether there is enough room for the percolate water to infiltrate      !
      ! (i.e., the soil will not become super-saturated), in which case we must reduce the !
      ! total amount of percolation.                                                       !
      !------------------------------------------------------------------------------------!
      if (k == 1) then
         !---------------------------------------------------------------------------------!
         !     Compute the available "room" for water at the top soil layer.  We must      !
         ! multiply by density and depth to make sure that the units match.                !
         !---------------------------------------------------------------------------------!
         wmass_room = max(0.0, slmsts(nsoil) - soil_water(mzg)) * wdns * dslz(mzg) 
         wmass_perc = min(wmass_perc,wmass_room)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Re-calculate the total water mass and energy of this temporary surface water.  !
      ! Here we must check whether the soil layer would be with too little mass, and if    !
      ! that is the case, we will eliminate the layer by forcing the tiny left-over to go  !
      ! to the layer beneath.                                                              !
      !------------------------------------------------------------------------------------!
      if (wmass_try - wmass_perc > min_sfcwater_mass) then
         !---------------------------------------------------------------------------------!
         !      Enough mass to keep this layer.                                            !
         !---------------------------------------------------------------------------------!
         !----- Compute the internal energy and depth associated with percolated water. ---!
         energy_perc = wmass_perc * cliq * (temp_try - tsupercool)
         depth_perc  = wmass_perc * wdnsi
         !----- Find the new water mass and energy for this layer. ------------------------!
         sfcwater_mass      (k) = wmass_try  - wmass_perc
         sfcwater_energy_ext(k) = energy_try - energy_perc
         !---------------------------------------------------------------------------------!
         !      Calculate density and depth of snow.  Start with the difference of depths, !
         ! but then we adjust it because the loss through percolation changes the ratio    !
         ! between ice and liquid in this layer                                            !
         !---------------------------------------------------------------------------------!
         sfcwater_depth     (k) = depth_try  - depth_perc
         snden                  = sfcwater_mass(k)                                         &
                                / max(min_sfcwater_depth,sfcwater_depth(k))
         sndenmax           = wdns
         sndenmin           = max(30., 200. * (wmass_free + wmass_perc) / sfcwater_mass(k))
         snden              = min(sndenmax, max(sndenmin,snden))
         sfcwater_depth (k) = sfcwater_mass(k) / snden
      else
         !---------------------------------------------------------------------------------!
         !      The layer would be too small, eliminate mass from this layer and send all  !
         ! mass to the layer beneath as percolated water.                                  !
         !---------------------------------------------------------------------------------!
         sfcwater_mass      (k) = 0.0
         sfcwater_energy_ext(k) = 0.0
         sfcwater_depth     (k) = 0.0
         wmass_perc             = wmass_try
         energy_perc            = energy_try
         depth_perc             = depth_try
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Integrate the total temporary surface water properties.                        !
      !------------------------------------------------------------------------------------!
      sum_sfcw_mass   = sum_sfcw_mass   + sfcwater_mass      (k)
      sum_sfcw_energy = sum_sfcw_energy + sfcwater_energy_ext(k)
      sum_sfcw_depth  = sum_sfcw_depth  + sfcwater_depth     (k)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The water available for the layer beneath is going to be the total percolated  !
      ! water.                                                                             !
      !------------------------------------------------------------------------------------!
      wmass_free  = wmass_perc
      energy_free = energy_perc
      depth_free  = depth_perc
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Add any remaining free water to the top soil layer.                               !
   !---------------------------------------------------------------------------------------!
   soil_water(mzg)  = soil_water(mzg)  + wmass_free  * dslzi(mzg) * wdnsi
   soil_energy(mzg) = soil_energy(mzg) + energy_free * dslzi(mzg)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Check the total amount of mass in the temporary surface water/snow, and adjust    !
   ! the number of layer accordingly.                                                      !
   !---------------------------------------------------------------------------------------!
   if (sum_sfcw_mass <= min_sfcwater_mass) then
      !----- Not enough water in the temporary surface water, eliminate all layers. -------!
      sfcwater_nlev = 0.
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      The total mass should be either zero or greater than rk4tiny_sfcw_mass,       !
      ! but, just in case, we add any remaining energy to the top soil layer.              !
      !------------------------------------------------------------------------------------!
      soil_water(mzg)  = soil_water(mzg)  + sum_sfcw_mass   * dslzi(mzg) * wdnsi
      soil_energy(mzg) = soil_energy(mzg) + sum_sfcw_energy * dslzi(mzg)
      !------------------------------------------------------------------------------------!

      !----- Loop all layers and re-set all extensive variables to zero. ------------------!
      do k = 1, mzs
         sfcwater_mass      (k) = 0.0
         sfcwater_energy_ext(k) = 0.0
         sfcwater_depth     (k) = 0.0
      end do
      !------------------------------------------------------------------------------------!
   elseif (sum_sfcw_mass < water_stab_thresh) then


      !----- Not much water in the temporary surface water, impose a single layer. --------!
      sfcwater_nlev = 1.
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    If the total amount of temporary surface water is not enough to make it stable, !
      ! we impose it to have a single layer with all the ponding/snow in there.            !
      !------------------------------------------------------------------------------------!
      sfcwater_mass       (1) = sum_sfcw_mass
      sfcwater_energy_ext (1) = sum_sfcw_energy
      sfcwater_depth      (1) = sum_sfcw_depth
      do k=2,mzs
         sfcwater_mass       (k) = 0.0
         sfcwater_energy_ext (k) = 0.0
         sfcwater_depth      (k) = 0.0
      end do
      !------------------------------------------------------------------------------------!


   else
      !---- Re-set the new layers buffer. -------------------------------------------------!
      newsfcw_mass  (:) = 0.
      newsfcw_energy(:) = 0.
      newsfcw_depth (:) = 0.


      !---- Check whether there is enough snow for a new layer. ---------------------------!
      nlayers   = ksnnew
      newlayers = 1
      do k = 1,ksnnew
         !---------------------------------------------------------------------------------!
         !     Check whether the layer as is meet the minimum requirements to stand as a   !
         ! new layer by itself.                                                            !
         !---------------------------------------------------------------------------------!
         if ( sfcwater_mass(k)                > min_sfcwater_mass        .and.             &
              water_stab_thresh * thicknet(k) <= sum_sfcw_mass           .and.             &
              sfcwater_energy_ext(k)          <  sfcwater_mass(k)*qliqt3       ) then
            newlayers = newlayers + 1
         end if
         !---------------------------------------------------------------------------------!
      end do

      !----- Newlayers is the new number of temporary surface water layers. ---------------!
      newlayers = min(newlayers, mzs, nlayers + 1)
      
      if (newlayers == 1) then
         newsfcw_mass  (1) = sum_sfcw_mass
         newsfcw_energy(1) = sum_sfcw_energy
         newsfcw_depth (1) = sum_sfcw_depth
      else
         kold  = 1
         wtnew = 1.0
         wtold = 1.0
         do k = 1,newlayers
            newsfcw_mass(k)   = sum_sfcw_mass * thick(k,newlayers)
            newsfcw_energy(k) = 0.0
            newsfcw_depth(k)  = 0.0
            !----- Find the properties of this new layer. ---------------------------------!
            find_layer: do

               !----- Difference between old and new snow ---------------------------------!
               wdiff = wtnew * newsfcw_mass(k) - wtold * sfcwater_mass(kold)  

               if (wdiff > 0.0) then
                  newsfcw_energy(k) = newsfcw_energy(k)                                    &
                                    + wtold * sfcwater_energy_ext(kold)
                  newsfcw_depth(k)  = newsfcw_depth(k)                                     &
                                    + wtold * sfcwater_depth(kold)
                  wtnew  = wtnew - wtold * sfcwater_mass(kold) / newsfcw_mass(k)
                  kold   = kold + 1
                  wtold  = 1.0
                  if (kold > nlayers) exit find_layer
               else
                  newsfcw_energy(k) = newsfcw_energy(k) + wtnew * newsfcw_mass(k)          &
                                    * sfcwater_energy_ext(kold)                            &
                                    / max(min_sfcwater_mass,sfcwater_mass(kold))
                  newsfcw_depth(k)  = newsfcw_depth(k)  + wtnew * newsfcw_mass(k)          &
                                    * sfcwater_depth(kold)                                 &
                                    / max(min_sfcwater_mass,sfcwater_mass(kold))
                  wtold = wtold - wtnew * newsfcw_mass(k)                                     &
                                / max(min_sfcwater_mass,sfcwater_mass(kold))
                  wtnew = 1.
                  exit find_layer
               end if
            end do find_layer
         end do
      end if

      !----- Update the water/snow layer prognostic properties. ---------------------------!
      sfcwater_nlev = real(newlayers)
      do k = 1,newlayers
         sfcwater_mass      (k) = newsfcw_mass(k)
         sfcwater_energy_ext(k) = newsfcw_energy(k)
         sfcwater_depth     (k) = newsfcw_depth(k)
      end do
      do k = newlayers + 1, mzs
         sfcwater_mass      (k) = 0.0
         sfcwater_energy_ext(k) = 0.0
         sfcwater_depth     (k) = 0.0
      end do
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf_adjust_sfcw
!==========================================================================================!
!==========================================================================================!

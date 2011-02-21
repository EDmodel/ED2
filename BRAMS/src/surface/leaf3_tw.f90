!==========================================================================================!
!==========================================================================================!
subroutine leaftw(mzg,mzs,soil_water, soil_energy,soil_text,sfcwater_mass                  &
                 ,sfcwater_depth,ustar,tstar,rstar,cstar,zeta,ribulk,veg_albedo            &
                 ,veg_fracarea,veg_lai,veg_tai,veg_rough,veg_height,patch_area,patch_rough &
                 ,patch_wetind,leaf_class,soil_rough,sfcwater_nlev,stom_condct,ground_rsat &
                 ,ground_rvap,ground_temp,ground_fliq,veg_water,veg_hcap,veg_energy        &
                 ,can_prss,can_theiv,can_theta,can_rvap,can_co2,sensible_gc,sensible_vc    &
                 ,evap_gc,evap_vc,transp,gpp,plresp,resphet,veg_ndvip,veg_ndvic,veg_ndvif  &
                 ,rshort,cosz,ip,i,j)

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
   integer                                :: ip
   integer                                :: k
   integer                                :: nveg
   integer                                :: ksn
   integer                                :: nsoil
   integer                                :: ksnnew
   integer                                :: newlayers
   integer                                :: nlayers
   integer                                :: kold
   integer                                :: ktrans
   integer                                :: nsl
   integer                                :: kzs
   integer                                :: i
   integer                                :: j
   real                                   :: sfcw_energy_int
   real                                   :: sfcw_temp
   real                                   :: sfcw_fliq
   real                                   :: stretch
   real                                   :: thik
   real                                   :: pf
   real                                   :: snden
   real                                   :: vegfracc
   real                                   :: qwfree
   real                                   :: wfree
   real                                   :: depthgain
   real                                   :: totsnow
   real                                   :: qw
   real                                   :: w
   real                                   :: qwt
   real                                   :: wt
   real                                   :: soilhcap
   real                                   :: fac
   real                                   :: wfreeb
   real                                   :: depthloss
   real                                   :: soilcap
   real                                   :: sndenmax
   real                                   :: sndenmin
   real                                   :: snowmin
   real                                   :: hold
   real                                   :: wtnew
   real                                   :: wtold
   real                                   :: wdiff
   real                                   :: watermid
   real                                   :: available_water
   real   , dimension(mzg)                :: available_layer
   real                                   :: extracted_water
   real                                   :: wloss
   real                                   :: qwloss
   real                                   :: soilcond
   real                                   :: waterfrac
   real                                   :: psiplusz_bl
   !----- Locally saved variables. --------------------------------------------------------!
   logical                , save          :: ncall=.true.
   real, dimension(20)    , save          :: thicknet
   real, dimension(20,20) , save          :: thick
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

   !----- Initialize snow thickness scaling array. ----------------------------------------!
   if (ncall) then
      ncall = .false.
      stretch = 2.0
      do kzs = 1,mzs
         thik = 1.0
         thicknet(kzs) = 0.0
         do k = 1,(kzs+1)/2
            thick(k,kzs) = thik
            thick(kzs+1-k,kzs) = thik
            thicknet(kzs) = thicknet(kzs) + 2. * thik
            thik = thik * stretch
         end do
         if ((kzs+1)/2 /= kzs/2) thicknet(kzs) = thicknet(kzs) - thik/stretch
         do k = 1,kzs
            thick(k,kzs) = thick(k,kzs) / thicknet(kzs)
         end do
      end do
   end if

   nveg = nint(leaf_class)
   ksn = nint(sfcwater_nlev)


   !---------------------------------------------------------------------------------------!
   !    Remove water from all layers up to the root depth for transpiration.  Bottom       !
   ! k-index in root zone is kroot(nveg).  Limit soil moisture to values above soilwp.     !
   ! The available water profile is given in kg/m2.                                        !
   !---------------------------------------------------------------------------------------!
   nsl    = nint(soil_text(mzg))
   ktrans = kroot(nveg)
   available_layer(:)   = 0.
   available_water      = 0.
   do k = mzg,ktrans,-1
      nsoil              = nint(soil_text(k))
      available_layer(k) = dslz(k) * wdns                                                  &
                         * max(0.0, fracliq(k) * (soil_water(k)-soilwp(nsoil)))
      available_water    = available_water + available_layer(k)
   end do


   !---------------------------------------------------------------------------------------!
   !      Evaluate any exchanges of heat and moisture to or from vegetation, apply moist-  !
   ! ure and heat changes to vegetation, and evaluate the resistance parameter rgnd        !
   ! between canopy air and the top soil or snow surface.                                  !
   !---------------------------------------------------------------------------------------!
   call leaf_canopy(mzg,mzs,ksn,soil_energy,soil_water,soil_text,sfcwater_mass,ustar       &
                   ,tstar,rstar,cstar,zeta,ribulk,soil_rough,veg_rough,patch_rough         &
                   ,veg_height,veg_lai,veg_tai,veg_water,veg_hcap,veg_energy,leaf_class    &
                   ,veg_fracarea,stom_condct,can_prss,can_theiv,can_theta,can_rvap,can_co2 &
                   ,sensible_gc,sensible_vc,evap_gc,evap_vc,transp,gpp,plresp,resphet      &
                   ,ground_rsat,ground_rvap,ground_temp,ground_fliq,available_water,rshort &
                   ,i,j,ip)

   !----- Compute soil and effective snow heat resistance times layer depth (rfactor). ----!
   do k = 1,mzg
      nsoil     = nint(soil_text(k))
      waterfrac = soil_water(k) / slmsts(nsoil)
      soilcond  = soilcond0(nsoil)                                                         &
                + waterfrac * (soilcond1(nsoil) + waterfrac * soilcond2(nsoil))
      rfactor(k) = dslz(k) / soilcond
   end do

   do k = 1,ksn
      snden = sfcwater_mass(k) / sfcwater_depth(k)
      rfactor(k+mzg) = sfcwater_depth(k) / (1.093e-3 * exp(.028 * tempk(k+mzg))            &
                    * (.030 + snden * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))
   end do

   !----- Find soil and snow internal sensible heat fluxes [W/m2]. ------------------------!
   hfluxgsc(1) = 0.
   do k = 2,mzg+ksn
      hfluxgsc(k) = - (tempk(k) - tempk(k-1)) / ((rfactor(k) + rfactor(k-1)) * .5)
   end do

   !---------------------------------------------------------------------------------------!
   !     Heat flux at soil or snow top from longwave, sensible, and upward latent heat     !
   ! fluxes [W/m^2].                                                                       !
   !---------------------------------------------------------------------------------------!
   hfluxgsc(ksn+1+mzg) = hflxgc + qwflxgc - rlonga_gs - rlongv_gs + rlonggs_v + rlonggs_a

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
      sfcw_energy_ext(k) = sfcw_energy_ext(k)                                              &
                         + dtll * (hfluxgsc(k+mzg) - hfluxgsc(k+1+mzg) + rshort_s(k)) 
   end do

   if (ksn == 0) then
      soil_water(mzg) = soil_water(mzg) - wdnsi * wflxgc * dslzidt(mzg)
   else
      sfcwater_mass(ksn) = max(0.,sfcwater_mass(ksn) - wflxgc * dtll)
   endif

   !---------------------------------------------------------------------------------------!
   !      New moisture, qw, and depth from dew/frost formation, precipitation, shedding,   !
   ! and percolation.  ksnnew is the layer that receives the new condensate that comes     !
   ! directly from the air above.  If there is no pre-existing snowcover, this is a        !
   ! temporary snow/surface water layer.                                                   !
   !---------------------------------------------------------------------------------------!
   if (pcpgl + wshed_tot + dewgnd_tot > min_sfcwater_mass) then
      ksnnew = max(ksn,1)
      vegfracc = 1. - veg_fracarea
      qwfree    = qpcpgl * vegfracc + (qdewgnd_tot + qwshed_tot ) * veg_fracarea
      wfree     =  pcpgl * vegfracc + (dewgnd_tot  + wshed_tot  ) * veg_fracarea
      depthgain = dpcpgl * vegfracc + (dewgnd_tot  + wshed_tot  ) * veg_fracarea * wdnsi
   else
      ksnnew    = ksn
      qwfree    = 0.
      wfree     = 0.
      depthgain = 0.
   end if

   if (ksnnew > 0) then

      !------------------------------------------------------------------------------------!
      !     Transfer water downward through snow layers by percolation. Fracliq is the     !
      ! fraction of liquid in the snowcover or surface water. Wfree is the quantity of     !
      ! that liquid in kg/m2 which is free (not attached to snowcover) and therefore       !
      ! available to soak into the layer below). Soilcap is the capacity of the top soil   !
      ! layer in kg/m2 to accept surface water.  Wfree in the lowest snow layer is limited !
      ! by this value.                                                                     !
      !------------------------------------------------------------------------------------!
      totsnow = 0.
      nsoil = nint(soil_text(mzg))

      do k = ksnnew,1,-1
         qw = sfcw_energy_ext(k) + qwfree
         w  = sfcwater_mass(k)   + wfree

         !---------------------------------------------------------------------------------!
         !     If this is the top layer and there is water, get rid of some water through  !
         ! runoff.                                                                         !
         !---------------------------------------------------------------------------------!
         if (runoff_time > 0.0 .and. w > 0. .and. k == ksnnew) then
            call qtk(qw/w,tempk(k+mzg),fracliq(k+mzg))
            wloss  = max(0.,min(1.0,dtll/runoff_time) * w * (fracliq(k+mzg) - 0.1) / 0.9)
            qwloss = wloss * cliq * (tempk(k+mzg) - tsupercool8)
            qw = qw - qwloss
            w  = w  - wloss
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     If (only) snow layer is too thin for computational stability, bring it to   !
         ! thermal equilibrium with the top soil layer by exchanging heat between them.    !
         !---------------------------------------------------------------------------------!
         if (ksnnew == 1 .and. sfcwater_mass(k) < 3.) then
            qwt      = qw + soil_energy(mzg) * dslz(mzg)
            wt       = w  + soil_water(mzg) * wdns * dslz(mzg)
            soilhcap = slcpd(nsoil) * dslz(mzg)
            call qwtk(qwt,wt,soilhcap,tempk(k+mzg),fracliq(k+mzg))
            qw = w * ( fracliq(k+mzg)*cliq*(tempk(k+mzg)-tsupercool)                       &
                     + (1.-fracliq(k+mzg))*cice*tempk(k+mzg)  )
            tempk(mzg)       = tempk(k+mzg)
            fracliq(mzg)     = fracliq(k+mzg)
            soil_energy(mzg) = (qwt - qw) * dslzi(mzg)
         else
            call qtk(qw/w,tempk(k+mzg),fracliq(k+mzg))
         end if

         !---------------------------------------------------------------------------------!
         !    Shed liquid in excess of a 1:9 liquid-to-ice ratio.  Limit this shed amount  !
         ! (wfreeb) in lowest snow layer to amount top soil layer can hold.                !
         !---------------------------------------------------------------------------------!
         wfreeb    = max (0., w * (fracliq(k+mzg) - .1) / 0.9)
         depthloss = wfreeb * wdnsi
         if (k == 1) then
            soilcap = wdns * max (0.,-slz(mzg) * (slmsts(nsoil) - soil_water(mzg)))
            wfreeb  = min(wfreeb, soilcap)
            
            qwfree  = wfreeb * cliq * (tempk(k+mzg) - tsupercool)
            soil_water(mzg) = soil_water(mzg) + wdnsi * wfreeb * dslzi(mzg)
            soil_energy(mzg) = soil_energy(mzg) + qwfree * dslzi(mzg)
         else
            qwfree = wfreeb * cliq * (tempk(k+mzg) - tsupercool)
         end if

         sfcwater_mass(k)  = w - wfreeb
         sfcwater_depth(k) = sfcwater_depth(k) + depthgain - depthloss

         if (sfcwater_mass(k) < min_sfcwater_mass) then
            sfcw_energy_ext(k) = 0.
            sfcwater_mass(k)   = 0.
            sfcwater_depth(k)  = 0.
         else
            totsnow = totsnow + sfcwater_mass(k)
            sfcw_energy_ext(k) = qw - qwfree
         end if

         !----- Temporary simple evolution of snow layer depth and density. ---------------!
         sfcwater_depth(k) = sfcwater_depth(k) * (1. - dtll / 1.e5)
         snden             = sfcwater_mass(k)  / max(min_sfcwater_depth,sfcwater_depth(k))
         sndenmax          =  wdns
         sndenmin          = max( 30., 200. * (wfree + wfreeb)                             &
                                     / max(min_sfcwater_mass,sfcwater_mass(k)))
         snden             = min (sndenmax, max (sndenmin, snden))
         sfcwater_depth(k) = sfcwater_mass(k) / snden
         wfree             = wfreeb
         depthgain         = depthloss
      end do

      !----- Re-distribute snow layers to maintain prescribed distribution of mass. -------!
      if (totsnow < min_sfcwater_mass) then
         sfcwater_nlev      = 0.
         sfcwater_mass(1)   = 0.
         sfcw_energy_ext(1) = 0.
         sfcwater_depth(1)  = 0.
      else
         nlayers   = ksnnew
         snowmin   = 3.0
         newlayers = 1
         do k = 2,mzs
            if (snowmin * thicknet(k) <= totsnow .and.                                     &
                sfcw_energy_ext(k) < qliqt3*sfcwater_mass(k)) then
               newlayers = newlayers + 1
            end if
         end do
         newlayers     = min (newlayers, mzs, nlayers+1)
         sfcwater_nlev = float(newlayers)

         kold  = 1
         wtnew = 1.
         wtold = 1.
         do k = 1,newlayers
            vctr14(k) = totsnow * thick(k,newlayers)
            vctr16(k) = 0.
            vctr18(k) = 0.
            inner_loop: do
               wdiff = wtnew * vctr14(k) - wtold * sfcwater_mass(kold)
               if (wdiff > 0.) then
                  vctr16(k) = vctr16(k) + wtold * sfcw_energy_ext(kold)
                  vctr18(k) = vctr18(k) + wtold * sfcwater_depth(kold)
                  wtnew     = wtnew - wtold * sfcwater_mass(kold) / vctr14(k)
                  kold      = kold + 1
                  wtold     = 1.
                  if (kold <= ksn) cycle inner_loop
               else
                  vctr16(k) = vctr16(k) + wtnew * vctr14(k) * sfcw_energy_ext(kold)        &
                            / max(min_sfcwater_mass,sfcwater_mass(kold))
                  vctr18(k) = vctr18(k) + wtnew * vctr14(k) * sfcwater_depth(kold)         &
                            / max(min_sfcwater_mass,sfcwater_mass(kold))
                  wtold     = wtold - wtnew * vctr14(k) / sfcwater_mass(kold)
                  wtnew     = 1.
               end if
               exit inner_loop
            end do inner_loop
         end do

         do k = 1,newlayers
            sfcwater_mass(k)   = vctr14(k)
            sfcw_energy_ext(k) = vctr16(k) 
            sfcwater_depth(k)  = vctr18(k)
         end do
      end if
   end if


   !---------------------------------------------------------------------------------------!
   !     Compute gravitational potential plus moisture potential z + psi [m], liquid water !
   ! content [m], and half the remaining water capacity [m].                               !
   !---------------------------------------------------------------------------------------!
   do k = 1,mzg
      nsoil = nint(soil_text(k))
      psiplusz(k) = slzt(k) + slpots(nsoil) * (slmsts(nsoil)/soil_water(k))**slbs(nsoil)
      soil_liq(k) = dslz(k) * max(0.,(soil_water(k) - soilcp(nsoil))*fracliq(k))
      half_soilair(k) = (slmsts(nsoil) - soil_water(k)) * dslz(k) * .5
   end do

   !---------------------------------------------------------------------------------------!
   !     Find amount of water transferred between soil layers (wflux) [m] modulated by the !
   ! liquid water fraction.                                                                !
   !---------------------------------------------------------------------------------------!
   wflux(mzg+1)  = 0.
   qwflux(mzg+1) = 0.

   !----- Boundary condition at the lowest soil level -------------------------------------!
   nsoil = nint(soil_text(1))
   select case (isoilbc)
   case (0)
      !----- Bedrock, no flux across. -----------------------------------------------------!
      wflux(1)  = 0.
      qwflux(1) = 0.
   case (1)
      !------------------------------------------------------------------------------------!
      !     Free drainage, water movement of bottom soil layer is only under gravity, i.e. !
      ! the soil water content of boundary layer is equal to that of bottom soil layer.    !
      !------------------------------------------------------------------------------------!
      watermid = soil_water(1)
      wflux(1) = dslztidt(1) * slcons1(1,nsoil)                                            &
               * (watermid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)                     &
               * (slz(2) - slz(1)) * fracliq(1)
      !------------------------------------------------------------------------------------!
      !     Limit water transfers to prevent over-saturation and over-depletion. Compute q !
      ! transfers between soil layers (qwflux) [J/m2].                                     !
      !------------------------------------------------------------------------------------!
      if (wflux(1) > 0.) then
         wflux(1) = min(wflux(1),half_soilair(1))
      else
         wflux(1) = - min(-wflux(1),soil_liq(1))
      end if
      qwflux(1) = wflux(1) * cliqvlme * (tempk(1) - tsupercool)
      !------------------------------------------------------------------------------------!
   case (2)
      !------------------------------------------------------------------------------------!
      !     Half drainage, water movement of bottom soil layer is only under gravity, i.e. !
      ! the soil water content of boundary layer is equal to that of bottom soil layer.    !
      !------------------------------------------------------------------------------------!
      watermid = soil_water(1)
      wflux(1) = dslztidt(1) * slcons1(1,nsoil)                                            &
               * (watermid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)                     &
               * (slz(2) - slz(1)) * fracliq(1) * 0.5
      !------------------------------------------------------------------------------------!
      !     Limit water transfers to prevent over-saturation and over-depletion. Compute q !
      ! transfers between soil layers (qwflux) [J/m2].                                     !
      !------------------------------------------------------------------------------------!
      if (wflux(1) > 0.) then
         wflux(1) = min(wflux(1),half_soilair(1))
      else
         wflux(1) = - min(-wflux(1),soil_liq(1))
      end if
      qwflux(1) = wflux(1) * cliqvlme * (tempk(1) - tsupercool)
      !------------------------------------------------------------------------------------!
   case (3)
      !------------------------------------------------------------------------------------!
      !     Free drainage, water movement of bottom soil layer is under gravity and moist- !
      ! ure potential difference. The soil water content of boundary layer is equal to     !
      ! field capacity.                                                                    !
      !------------------------------------------------------------------------------------!
      watermid = soil_water(1)
      psiplusz_bl = slz(1)                                                                 &
                  + slpots(nsoil) * (slmsts(nsoil) / sfldcap(nsoil)) ** slbs(nsoil)

      if (psiplusz_bl <= psiplusz(1)) then
          wflux(1) =  dslztidt(1) * slcons1(1,nsoil)                                       &
                   * (watermid/slmsts(nsoil)) ** (2.0 * slbs(nsoil) + 3.0)                 &
                   * (psiplusz(1) - psiplusz_bl) * fracliq(1)
      else
          !----- Prevent bottom soil layer sucking water from the boundary layer. ---------!
          wflux(1) = 0.0
      end if
   end select
      
   do k = 2,mzg
      nsoil    = nint(soil_text(k))
      watermid = 0.5 * (soil_water(k) + soil_water(k-1))
      wflux(k) = dslztidt(k) * slcons1(k,nsoil)                                            &
               * (watermid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)                     &
               * (psiplusz(k-1) - psiplusz(k)) * .5 * (fracliq(k) + fracliq(k-1))

      !------------------------------------------------------------------------------------!
      !     Limit water transfers to prevent over-saturation and over-depletion. Compute q !
      ! transfers between soil layers (qwflux) [J/m2].                                     !
      !------------------------------------------------------------------------------------!
      if (wflux(k) > 0.) then
         wflux(k) = min(wflux(k),soil_liq(k-1),half_soilair(k))
      else
         wflux(k) = - min(-wflux(k),soil_liq(k),half_soilair(k-1))
      endif
      qwflux(k) = wflux(k) * cliqvlme * (tempk(k) - tsupercool)

   end do

   !----- Update soil moisture (impose minimum value of soilcp) and q value. --------------!
   do k = 1,mzg
      nsoil = nint(soil_text(k))
      soil_water(k)  = max(soilcp(nsoil),soil_water(k) - dslzi(k) * (wflux(k+1)-wflux(k)))
      soil_energy(k) = soil_energy(k) - dslzi(k) * (qwflux(k+1)-qwflux(k))
   enddo

   !---------------------------------------------------------------------------------------!
   !     Update soil moisture by extracting the water lost due to transpiration.           !
   !---------------------------------------------------------------------------------------!
   extracted_water = transp_tot * dtll
   if (extracted_water > 0. .and. available_water > 0.) then
      do k = ktrans,mzg
         wloss  = wdnsi * dslzi(k) * extracted_water * available_layer(k) / available_water
         qwloss = wloss * cliqvlme * (tempk(k) - tsupercool)
         
         soil_water(k)  = soil_water(k)  - wloss
         soil_energy(k) = soil_energy(k) - qwloss
      end do
   end if


   !---------------------------------------------------------------------------------------!
   ! Compute ground vap mxrat for availability on next timestep; put into ground_rsat.     !
   !---------------------------------------------------------------------------------------!
   ksn = nint(sfcwater_nlev)
   if (ksn == 0) then
      sfcw_energy_int = 0.
   else
      if (sfcwater_mass(ksn) > min_sfcwater_mass) then
         sfcw_energy_int = sfcw_energy_ext(ksn) / sfcwater_mass(ksn)
      else 
         sfcw_energy_int = 0.
      end if
   end if
   call leaf_grndvap(soil_energy(mzg),soil_water(mzg),soil_text(mzg),sfcw_energy_int       &
                    ,sfcwater_nlev,can_rvap,can_prss,ground_rsat,ground_rvap,ground_temp   &
                    ,ground_fliq)

   return
end subroutine leaftw
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine initialises the extensive version of the surface water internal     !
! energy.  We will integrate the extensive variable rather than the intensive, since the   !
! latter has a strong dependence on sfcwater_mass and this increases errors.  The          !
! extensive variable is more accurate because the heat balance is done using extensive     !
! fluxes (W/m²), not (W/kg). After each iteration the intensive quantity is updated so     !
! routines like sfcrad will work fine.                                                     !
!------------------------------------------------------------------------------------------!
subroutine sfcw_int_2_ext(mzs,sfcwater_nlev,sfcwater_energy_int,sfcwater_mass)
   use leaf_coms, only : sfcw_energy_ext ! ! intent(out)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                , intent(in) :: mzs
   real                   , intent(in) :: sfcwater_nlev
   real   , dimension(mzs), intent(in) :: sfcwater_energy_int
   real   , dimension(mzs), intent(in) :: sfcwater_mass
   !----- Local variables. ----------------------------------------------------------------!
   integer                             :: ksn
   integer                             :: k
   !---------------------------------------------------------------------------------------!

   ksn = nint(sfcwater_nlev)

   do k=1,ksn
      sfcw_energy_ext(k) = sfcwater_energy_int(k) * sfcwater_mass(k)
   end do
   do k=ksn+1,mzs
      sfcw_energy_ext(k) = 0.
   end do

   return
end subroutine sfcw_int_2_ext
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine converts the updated internal energy surface water (in W/m²) to the !
! regular output unit (J/kg).                                                              !
!------------------------------------------------------------------------------------------!
subroutine sfcw_ext_2_int(mzs,sfcwater_nlev,sfcwater_energy_int,sfcwater_mass)
   use leaf_coms, only : sfcw_energy_ext ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                , intent(in)  :: mzs
   real                   , intent(in)  :: sfcwater_nlev
   real   , dimension(mzs), intent(out) :: sfcwater_energy_int
   real   , dimension(mzs), intent(in)  :: sfcwater_mass
   !----- Local variables. ----------------------------------------------------------------!
   integer                              :: ksn
   integer                              :: k
   !---------------------------------------------------------------------------------------!

   ksn = nint(sfcwater_nlev)

   do k=1,ksn
      sfcwater_energy_int(k) = sfcw_energy_ext(k) / sfcwater_mass(k)
   end do
   do k=ksn+1,mzs
      sfcwater_energy_int(k) = 0.
   end do

   return
end subroutine sfcw_ext_2_int
!==========================================================================================!
!==========================================================================================!

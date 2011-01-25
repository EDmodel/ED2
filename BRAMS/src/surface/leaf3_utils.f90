!==========================================================================================!
!==========================================================================================!
!                                                                                          !
!     This subroutine computes the characteristic scales, using one of the following surf- !
! ace layer parametrisation.                                                               !
!                                                                                          !
! 1. Based on L79;                                                                         !
! 2. Based on: OD95, but with some terms computed as in L79 and B71 to avoid singular-     !
!    ities.                                                                                !
! 3. Based on BH91, using an iterative method to find zeta, and using the modified         !
!    equation for stable layers.                                                           !
!                                                                                          !
! References:                                                                              !
! B71.  BUSINGER, J.A, et. al; Flux-Profile relationships in the atmospheric surface       !
!           layer. J. Atmos. Sci., 28, 181-189, 1971.                                      !
! L79.  LOUIS, J.F.; Parametric Model of vertical eddy fluxes in the atmosphere.           !
!           Boundary-Layer Meteor., 17, 187-202, 1979.                                     !
! BH91. BELJAARS, A.C.M.; HOLTSLAG, A.A.M.; Flux parameterization over land surfaces for   !
!           atmospheric models. J. Appl. Meteor., 30, 327-341, 1991.                       !
! OD95. ONCLEY, S.P.; DUDHIA, J.; Evaluation of surface fluxes from MM5 using observa-     !
!           tions.  Mon. Wea. Rev., 123, 3344-3357, 1995.                                  !
!------------------------------------------------------------------------------------------!
subroutine leaf_stars(theta_atm,theiv_atm,shv_atm,rvap_atm,co2_atm                         &
                     ,theta_can,theiv_can,shv_can,rvap_can,co2_can                         &
                     ,zref,uref,dtll,rough,ustar,tstar,estar,qstar,rstar,cstar             &
                     ,zeta,rib,r_aer)
   use mem_leaf  , only : istar      ! ! intent(in)
   use rconstants, only : grav       & ! intent(in)
                        , vonk       & ! intent(in)
                        , epim1      & ! intent(in)
                        , halfpi     ! ! intent(in)
   use leaf_coms , only : ustmin     & ! intent(in)
                        , bl79       & ! intent(in)
                        , csm        & ! intent(in)
                        , csh        & ! intent(in)
                        , dl79       & ! intent(in)
                        , ribmaxod95 & ! intent(in)
                        , ribmaxbh91 & ! intent(in)
                        , tprandtl   & ! intent(in)
                        , z0moz0h    & ! intent(in)
                        , z0hoz0m    & ! intent(in)
                        , psim       & ! function
                        , psih       & ! function
                        , zoobukhov  ! ! function
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real, intent(in)  :: theta_atm    ! Above canopy air pot. temperature        [        K]
   real, intent(in)  :: theiv_atm    ! Above canopy air eq. pot. temp           [        K]
   real, intent(in)  :: shv_atm      ! Above canopy vapour spec. hum.           [kg/kg_air]
   real, intent(in)  :: rvap_atm     ! Above canopy vapour mixing ratio         [kg/kg_air]
   real, intent(in)  :: co2_atm      ! Above canopy CO2 mixing ratio            [ µmol/mol]
   real, intent(in)  :: theta_can    ! Canopy air potential temperature         [        K]
   real, intent(in)  :: theiv_can    ! Canopy air equiv. potential temperature  [        K]
   real, intent(in)  :: shv_can      ! Canopy air vapour spec. humidity         [kg/kg_air]
   real, intent(in)  :: rvap_can     ! Canopy air vapour mixing ratio           [kg/kg_air]
   real, intent(in)  :: co2_can      ! Canopy air CO2 mixing ratio              [ µmol/mol]
   real, intent(in)  :: zref         ! Height at reference point                [        m]
   real, intent(in)  :: uref         ! Wind speed at reference height           [      m/s]
   real, intent(in)  :: dtll         ! Time step                                [        m]
   real, intent(in)  :: rough        ! z0, the roughness                        [        m]
   real, intent(out) :: ustar        ! U*, friction velocity                    [      m/s]
   real, intent(out) :: tstar        ! Temperature friction scale               [        K]
   real, intent(out) :: estar        ! Theta_Eiv friction scale                 [        K]
   real, intent(out) :: qstar        ! Specific humidity friction scale         [kg/kg_air]
   real, intent(out) :: rstar        ! Vapour mixing ratio friction scale       [kg/kg_air]
   real, intent(out) :: cstar        ! CO2 mixing ratio                         [ µmol/mol]
   real, intent(out) :: zeta         ! z/(Obukhov length).                      [      ---]
   real, intent(out) :: rib          ! Bulk richardson number.                  [      ---]
   real, intent(out) :: r_aer        ! Aerodynamic resistance                   [      s/m]
   !----- Local variables, common to both models. -----------------------------------------!
   logical           :: stable       ! Stable state
   real              :: zoz0m        ! zref/rough(momentum)
   real              :: lnzoz0m      ! ln[zref/rough(momentum)]
   real              :: zoz0h        ! zref/rough(heat)
   real              :: lnzoz0h      ! ln[zref/rough(heat)]
   real              :: c3           ! coefficient to find the other stars
   real              :: delz         !
   real              :: d_vel        !
   real              :: vel_new      !
   !----- Local variables, used by L79. ---------------------------------------------------!
   real              :: a2           ! Drag coefficient in neutral conditions
   real              :: c1           ! a2 * vels
   real              :: fh           ! Stability parameter for heat
   real              :: fm           ! Stability parameter for momentum
   real              :: c2           ! Part of the c coefficient common to momentum & heat.
   real              :: cm           ! c coefficient times |Rib|^1/2 for momentum.
   real              :: ch           ! c coefficient times |Rib|^1/2 for heat.
   real              :: ee           ! (z/z0)^1/3 -1. for eqn. 20 w/o assuming z/z0 >> 1.
   !----- Local variables, used by OD95 and/or BH91. --------------------------------------!
   real              :: zeta0m       ! roughness(momentum)/(Obukhov length).
   real              :: zeta0h       ! roughness(heat)/(Obukhov length).
   real              :: ribold       ! Bulk richardson number.
   !----- Aux. environment conditions. ----------------------------------------------------!
   real              :: thetav_atm   ! Atmos. virtual potential temperature     [        K]
   real              :: thetav_can   ! Canopy air virtual pot. temperature      [        K]
   !----- External functions. -------------------------------------------------------------!
   real, external    :: cbrt         ! Cubic root
   !---------------------------------------------------------------------------------------!

   !----- Finding the variables common to both methods. -----------------------------------!
   thetav_atm = theta_atm * (1. + epim1 * shv_atm)
   thetav_can = theta_can * (1. + epim1 * shv_can)
   zoz0m      = zref/rough
   lnzoz0m    = log(zoz0m)
   zoz0h      = z0moz0h * zoz0m
   lnzoz0h    = log(zoz0h)
   rib        = 2.0 * grav * zref * (thetav_atm-thetav_can)                                &
              / ( (thetav_atm+thetav_can) * uref * uref)
   stable     = thetav_atm >= thetav_can

   !---------------------------------------------------------------------------------------!
   !     Here we find u* and the coefficient to find the other stars based on the chosen   !
   ! surface model.                                                                        !
   !---------------------------------------------------------------------------------------!
   select case (istar)
   case (1)
      !------------------------------------------------------------------------------------!
      !     Here we will use L79 model, the BRAMS default.                                 !
      !------------------------------------------------------------------------------------!

      !----- Compute the a-square factor and the coefficient to find theta*. --------------!
      a2   = (vonk / lnzoz0m) ** 2.
      c1   = a2 * uref

      if (stable) then
         !---------------------------------------------------------------------------------!
         !     Stable case.                                                                !
         !---------------------------------------------------------------------------------!
         fm = 1.0 / (1.0 + (2.0 * bl79 * rib / sqrt(1.0 + dl79 * rib)))
         fh = 1.0 / (1.0 + (3.0 * bl79 * rib * sqrt(1.0 + dl79 * rib)))

      else
         !---------------------------------------------------------------------------------!
         !     Unstable case.  The only difference from the original method is that we no  !
         ! longer assume z >> z0, so the "c" coefficient uses the full z/z0 term.          !
         !---------------------------------------------------------------------------------!
         ee = cbrt(zoz0m) - 1.
         c2 = bl79 * a2 * ee * sqrt(ee * abs(rib))
         cm = csm * c2
         ch = csh * c2
         fm = (1.0 - 2.0 * bl79 * rib / (1.0 + 2.0 * cm))
         fh = (1.0 - 3.0 * bl79 * rib / (1.0 + 3.0 * ch))
      end if
      r_aer = 1. / (a2 * uref * fh)

      !----- Finding ustar, making sure it is not too small. ------------------------------!
      ustar = max(ustmin,sqrt(c1 * uref * fm))

      !----- Finding the coefficient to scale the other stars. ----------------------------!
      c3 = c1 * fh / ustar
      !------------------------------------------------------------------------------------!


      !----- Compute zeta from u* and T* --------------------------------------------------!
      zeta = grav * vonk * c3 * (theta_atm - theta_can) / (theta_atm * ustar * ustar)


   case (2,4)
      !------------------------------------------------------------------------------------!
      ! 2. Here we use the model proposed by OD95, the standard for MM5, but with some     !
      !    terms that were computed in B71 (namely, the "0" terms). which prevent sin-     !
      !    gularities.  Since we use OD95 to estimate zeta, which avoids the computation   !
      !    of the Obukhov length L , we can't compute zeta0 by its definition(z0/L). How-  !
      !    ever we know zeta, so zeta0 can be written as z0/z * zeta.                      !
      ! 4. We use the model proposed by BH91, but we find zeta using the approximation     !
      !    given by OD95.                                                                  !
      !------------------------------------------------------------------------------------!

      !----- Make sure that the bulk Richardson number is not above ribmax. ---------------!
      rib = min(rib,ribmaxod95)
      
      !----- We now compute the stability correction functions. ---------------------------!
      if (stable) then
         !----- Stable case. --------------------------------------------------------------!
         zeta  = rib * lnzoz0m / (1.1 - 5.0 * rib)
      else
         !----- Unstable case. ------------------------------------------------------------!
         zeta = rib * lnzoz0m
      end if
      zeta0m = rough * zeta / zref

      !----- Finding the aerodynamic resistance similarly to L79. -------------------------!
      r_aer = tprandtl * (lnzoz0m - psih(zeta,stable) + psih(zeta0m,stable))               &
                       * (lnzoz0m - psim(zeta,stable) + psim(zeta0m,stable))               &
                       / (vonk * vonk * uref)

      !----- Finding ustar, making sure it is not too small. ------------------------------!
      ustar = max (ustmin, vonk * uref                                                     &
                         / (lnzoz0m - psim(zeta,stable) + psim(zeta0m,stable)))

      !----- Finding the coefficient to scale the other stars. ----------------------------!
      c3    = vonk / (tprandtl * (lnzoz0m - psih(zeta,stable) + psih(zeta0m,stable)))

      !------------------------------------------------------------------------------------!

   case (3)
      !------------------------------------------------------------------------------------!
      !      Here we use the model proposed by BH91, which is almost the same as the OD95  !
      ! method, with the two following (important) differences.                            !
      ! 1. Zeta (z/L) is actually found using the iterative method.                        !
      ! 2. Stable functions are computed in a more generic way.  BH91 claim that the       !
      !    oft-used approximation (-beta*zeta) can cause poor ventilation of the stable    !
      !    layer, leading to decoupling between the atmosphere and the canopy air space    !
      !    and excessive cooling.                                                          !
      ! 3. Here we distinguish the fluxes between roughness for momentum and for heat, as  !
      !    BH91 did.                                                                       !
      !------------------------------------------------------------------------------------!

      !----- Make sure that the bulk Richardson number is not above ribmax. ---------------!
      ribold = rib
      rib    = min(rib,ribmaxbh91)


      !----- We now compute the stability correction functions. ---------------------------!
      zeta   = zoobukhov(rib,zref,rough,zoz0m,lnzoz0m,zoz0h,lnzoz0h,stable)
      zeta0m = rough * zeta / zref
      zeta0h = z0hoz0m * zeta0m

      !----- Finding the aerodynamic resistance similarly to L79. -------------------------!
      r_aer = tprandtl * (lnzoz0h - psih(zeta,stable) + psih(zeta0h,stable))               &
                       * (lnzoz0m - psim(zeta,stable) + psim(zeta0m,stable))               &
                       / (vonk * vonk * uref)

      !----- Finding ustar, making sure it is not too small. ------------------------------!
      ustar = max (ustmin, vonk * uref                                                     &
                         / (lnzoz0m - psim(zeta,stable) + psim(zeta0m,stable)))


      !----- Finding the coefficient to scale the other stars. ----------------------------!
      c3    = vonk / (tprandtl * (lnzoz0h - psih(zeta,stable) + psih(zeta0h,stable)))
      !------------------------------------------------------------------------------------!

   end select

   !----- Finding all stars. --------------------------------------------------------------!
   tstar = c3 *    (theta_atm - theta_can)
   estar = c3 * log(theiv_atm / theiv_can)
   qstar = c3 *    (shv_atm   - shv_can  )
   rstar = c3 *    (rvap_atm  - rvap_can )
   cstar = c3 *    (co2_atm   - co2_can  )

   return

   if (abs(tstar) < 1.e-7) tstar = 0.
   if (abs(estar) < 1.e-7) estar = 0.
   if (abs(qstar) < 1.e-7) qstar = 0.
   if (abs(rstar) < 1.e-7) rstar = 0.
   if (abs(cstar) < 1.e-7) cstar = 0.

   !---------------------------------------------------------------------------------------!
   !     Limit ustar so that the flux cannot take more than 1/2 velocity in a timestep.    !
   !---------------------------------------------------------------------------------------!
   delz  = 2. * zref
   d_vel = - ustar * ustar * dtll / delz
   vel_new = uref + d_vel
   if (vel_new < .5 * uref) then
      d_vel = .5 * uref
      ustar = sqrt(d_vel * delz / dtll)
   end if
   !---------------------------------------------------------------------------------------!
   return
end subroutine leaf_stars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine computes the turbulent fluxes of momentum, heat and moisture from the    !
! surface layer using the  Manton-Cotton algebraic surface layer equations.                !
!------------------------------------------------------------------------------------------!
subroutine sfclmcv(ustar,tstar,rstar,cstar,zeta,vels,vels_pat,ups,vps,patch_area           &
                  ,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r,sflux_c,g_urban)
   use rconstants
   use teb_spm_start , only : teb_spm ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real , intent(in)    :: ustar,tstar,rstar,cstar,zeta
   real , intent(in)    :: vels,vels_pat,ups,vps,patch_area
   real , intent(inout) :: sflux_u,sflux_v,sflux_w,sflux_t,sflux_r,sflux_c
   real , intent(inout) :: g_urban
   !----- Local variables. ----------------------------------------------------------------!
   real                 :: cosine1,sine1,vtscr,cx,psin
   !----- Local constants. ----------------------------------------------------------------!
   real , parameter     :: wtol = 1.e-20
   !---------------------------------------------------------------------------------------!

   cosine1 = ups / vels_pat
   sine1   = vps / vels_pat

   vtscr = ustar * patch_area

   !----- Check whether TEB is used. If it is, skip the lines unless g_urban is zero. -----!
   if (teb_spm /= 1 .or. nint(g_urban) == 0) then
      sflux_u = sflux_u - ustar * vtscr * cosine1
      sflux_v = sflux_v - ustar * vtscr * sine1
      sflux_t = sflux_t - tstar * vtscr
      sflux_r = sflux_r - rstar * vtscr
   end if

   !----- TEB currently doesn't save CO2, so compute sflux_c outside the if statement. ----!
   sflux_c = sflux_c - cstar * vtscr

   !----- Define cx based on the layer stability. -----------------------------------------!
   if (zeta < 0.)then
      cx = zeta * sqrt(sqrt(1. - 15. * zeta))
   else
      cx = zeta / (1.0 + 4.7 * zeta)
   end if

   psin    = sqrt((1.-2.86 * cx) / (1. + cx * (-5.39 + cx * 6.998 )))
   sflux_w = sflux_w + (0.27 * max(6.25 * (1. - cx) * psin,wtol) - 1.18 * cx * psin)       &
                     * ustar * vtscr

   return
end subroutine sfclmcv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the ground properties.  By ground we mean the top soil      !
! layer if no temporary surface water/snow exists, or the top temporary surface water/snow !
! layer if it exists.                                                                      !
!------------------------------------------------------------------------------------------!
subroutine leaf_grndvap(topsoil_energy,topsoil_water,topsoil_text,sfcw_energy_int          &
                       ,sfcwater_nlev,can_rvap,can_prss,ground_rsat,ground_rvap            &
                       ,ground_temp,ground_fliq)

   use leaf_coms  , only : slcpd       & ! intent(in)
                         , slpots      & ! intent(in)
                         , slmsts      & ! intent(in)
                         , soilcp      & ! intent(in)
                         , slbs        & ! intent(in)
                         , sfldcap     ! ! intent(in)
   use rconstants , only : gorh2o      & ! intent(in)
                         , pi1         & ! intent(in)
                         , wdns        & ! intent(in)
                         , lnexp_min   ! ! intent(in)
   use therm_lib  , only : rslif       & ! function
                         , qwtk        & ! function
                         , qtk         ! ! function
   use mem_leaf   , only : betapower   ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real, intent(in)  :: topsoil_energy   ! Top soil internal energy             [     J/m³]
   real, intent(in)  :: topsoil_water    ! Top soil water content               [m³_h2o/m³]
   real, intent(in)  :: topsoil_text     ! Top soil texture class               [      ---]
   real, intent(in)  :: sfcw_energy_int  ! Soil internal energy                 [     J/kg]
   real, intent(in)  :: sfcwater_nlev    ! # active levels of surface water     [      ---]
   real, intent(in)  :: can_rvap         ! Canopy vapour mixing ratio           [kg_vap/kg]
   real, intent(in)  :: can_prss         ! Canopy pressure                      [       Pa]
   real, intent(out) :: ground_rsat      ! Surface (saturation) mixing ratio    [kg_vap/kg]
   real, intent(out) :: ground_rvap      ! Ground equilibrium mixing ratio      [kg_vap/kg]
   real, intent(out) :: ground_temp      ! Surface temperature                  [        K]
   real, intent(out) :: ground_fliq      ! Fraction of sfc H2O in liquid phase  [      ---]
   !----- Local variables. ----------------------------------------------------------------!
   integer           :: ksn              ! # active levels of surface water
   integer           :: nsoil            ! Soil texture class                   [      ---]
   real              :: slpotvn          ! soil water potential                 [        m]
   real              :: alpha            ! "alpha" term in Lee and Pielke (1993)
   real              :: beta             ! "beta" term in Lee and Pielke (1993)
   real              :: lnalpha          ! ln(alpha)
   real              :: smterm           ! soil moisture term                   [     ----]
   !---------------------------------------------------------------------------------------!


   !----- Set the number of temporary surface water (or snow) layers. ---------------------!
   ksn = nint(sfcwater_nlev)


   !---------------------------------------------------------------------------------------!
   !    Ground_rsat is the saturation mixing ratio of the top soil/snow surface and is     !
   ! used for dew formation and snow evaporation.  
   !---------------------------------------------------------------------------------------!
   select case (ksn)
   case (0)
      !------------------------------------------------------------------------------------!
      !      Without snowcover or water ponding, ground_shv is the effective specific      !
      ! humidity of soil and is used for soil evaporation.  This value is a combination of !
      ! the canopy air specific humidity, the saturation specific humidity at the soil     !
      ! temperature.  When the soil tends to dry air soil moisture, ground_shv tends to    !
      ! the canopy air space specific humidity, whereas it tends to the saturation value   !
      ! when the soil moisture is near or above field capacity.  These tendencies will be  !
      ! determined by the alpha and beta parameters.                                       !
      !------------------------------------------------------------------------------------!
      nsoil = nint(topsoil_text)
      call qwtk(topsoil_energy,topsoil_water*wdns,slcpd(nsoil),ground_temp,ground_fliq)
      !----- Compute the saturation mixing ratio at ground temperature. -------------------!
      ground_rsat = rslif(can_prss,ground_temp)
      !----- Determine alpha. -------------------------------------------------------------!
      slpotvn  = slpots(nsoil) * (slmsts(nsoil) / topsoil_water) ** slbs(nsoil)
      lnalpha  = gorh2o * slpotvn / ground_temp
      if (lnalpha > lnexp_min) then
         alpha = exp(lnalpha)
      else
         alpha = 0.0
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Determine Beta, following NP89.  However, because we want evaporation to be    !
      ! shut down when the soil approaches the dry air soil moisture, we offset both the   !
      ! soil moisture and field capacity to the soil moisture above dry air soil.  This is !
      ! necessary to avoid evaporation to be large just slightly above the dry air soil,   !
      ! which was happening especially for those clay-rich soil types.                     !
      !------------------------------------------------------------------------------------!
      smterm     = (topsoil_water - soilcp(nsoil)) / (sfldcap(nsoil) - soilcp(nsoil))
      beta       = (.5 * (1. - cos (min(1.,smterm) * pi1))) ** betapower
      !----- Use the expression from LP92 to determine the specific humidity. -------------!
      ground_rvap = ground_rsat * alpha * beta + (1. - beta) * can_rvap
      !------------------------------------------------------------------------------------!

   case default
      !------------------------------------------------------------------------------------!
      !    If a temporary layer exists, we use the top layer as the surface.  Since this   !
      ! is "pure" water or snow, we let it evaporate freely.  We can understand  this as   !
      ! the limit of alpha and beta tending to one.                                        !
      !------------------------------------------------------------------------------------!
      call qtk(sfcw_energy_int,ground_temp,ground_fliq)
      !----- Compute the saturation specific humidity at ground temperature. --------------!
      ground_rsat = rslif(can_prss,ground_temp)
      !----- The ground specific humidity in this case is just the saturation value. ------!
      ground_rvap = ground_rsat
      !------------------------------------------------------------------------------------!
   end select

   return
end subroutine leaf_grndvap
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will aplly the boundary condition to all leaf variables at the        !
! absolute domain edges.                                                                   !
!------------------------------------------------------------------------------------------!
subroutine leaf_bcond(m2,m3,mzg,mzs,npat,jdim,soil_water,sfcwater_mass,soil_energy         &
                     ,sfcwater_energy,soil_text,sfcwater_depth,ustar,tstar,rstar,cstar     &
                     ,zeta,ribulk,veg_albedo,veg_fracarea,veg_lai,veg_tai,veg_rough        &
                     ,veg_height,patch_area,patch_rough,patch_wetind,leaf_class,soil_rough &
                     ,sfcwater_nlev,stom_resist,ground_rsat,ground_rvap,ground_temp        &
                     ,ground_fliq,veg_water,veg_hcap,veg_energy,can_prss,can_theiv         &
                     ,can_theta,can_rvap,can_co2,sensible_gc,sensible_vc,evap_gc,evap_vc   &
                     ,transp,gpp,plresp,resphet,veg_ndvip,veg_ndvic,veg_ndvif)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                        , intent(in)    :: m2,m3,mzg,mzs,npat,jdim
   real, dimension(mzg,m2,m3,npat), intent(inout) :: soil_water,soil_energy,soil_text
   real, dimension(mzs,m2,m3,npat), intent(inout) :: sfcwater_mass,sfcwater_energy
   real, dimension(mzs,m2,m3,npat), intent(inout) :: sfcwater_depth
   real, dimension(m2,m3,npat)    , intent(inout) :: ustar,tstar,rstar,cstar,zeta,ribulk
   real, dimension(m2,m3,npat)    , intent(inout) :: veg_albedo,veg_fracarea
   real, dimension(m2,m3,npat)    , intent(inout) :: veg_lai,veg_tai,veg_rough,veg_height
   real, dimension(m2,m3,npat)    , intent(inout) :: patch_area,patch_rough,patch_wetind
   real, dimension(m2,m3,npat)    , intent(inout) :: leaf_class,soil_rough,sfcwater_nlev
   real, dimension(m2,m3,npat)    , intent(inout) :: stom_resist,ground_rsat,ground_rvap
   real, dimension(m2,m3,npat)    , intent(inout) :: ground_temp,ground_fliq
   real, dimension(m2,m3,npat)    , intent(inout) :: veg_water,veg_energy,veg_hcap
   real, dimension(m2,m3,npat)    , intent(inout) :: can_prss,can_theiv,can_theta
   real, dimension(m2,m3,npat)    , intent(inout) :: can_rvap,can_co2
   real, dimension(m2,m3,npat)    , intent(inout) :: sensible_gc,sensible_vc
   real, dimension(m2,m3,npat)    , intent(inout) :: evap_gc,evap_vc,transp
   real, dimension(m2,m3,npat)    , intent(inout) :: gpp,plresp,resphet
   real, dimension(m2,m3,npat)    , intent(inout) :: veg_ndvip,veg_ndvic,veg_ndvif
   !----- Local variables. ----------------------------------------------------------------!
   integer                                        :: i,j,k,ipat
   !---------------------------------------------------------------------------------------!
   do ipat = 1,npat
      do j = 1,m3

         ustar          (1,j,ipat) = ustar            (2,j,ipat)
         tstar          (1,j,ipat) = tstar            (2,j,ipat)
         rstar          (1,j,ipat) = rstar            (2,j,ipat)
         cstar          (1,j,ipat) = cstar            (2,j,ipat)
         zeta           (1,j,ipat) = zeta             (2,j,ipat)
         ribulk         (1,j,ipat) = ribulk           (2,j,ipat)
         veg_fracarea   (1,j,ipat) = veg_fracarea     (2,j,ipat)
         veg_lai        (1,j,ipat) = veg_lai          (2,j,ipat)
         veg_tai        (1,j,ipat) = veg_tai          (2,j,ipat)
         veg_rough      (1,j,ipat) = veg_rough        (2,j,ipat)
         veg_height     (1,j,ipat) = veg_height       (2,j,ipat)
         patch_area     (1,j,ipat) = patch_area       (2,j,ipat)
         patch_rough    (1,j,ipat) = patch_rough      (2,j,ipat)
         patch_wetind   (1,j,ipat) = patch_wetind     (2,j,ipat)
         leaf_class     (1,j,ipat) = leaf_class       (2,j,ipat)
         soil_rough     (1,j,ipat) = soil_rough       (2,j,ipat)
         sfcwater_nlev  (1,j,ipat) = sfcwater_nlev    (2,j,ipat)
         stom_resist    (1,j,ipat) = stom_resist      (2,j,ipat)
         ground_rsat    (1,j,ipat) = ground_rsat      (2,j,ipat)
         ground_rvap    (1,j,ipat) = ground_rvap      (2,j,ipat)
         ground_temp    (1,j,ipat) = ground_temp      (2,j,ipat)
         ground_fliq    (1,j,ipat) = ground_fliq      (2,j,ipat)
         veg_water      (1,j,ipat) = veg_water        (2,j,ipat)
         veg_hcap       (1,j,ipat) = veg_hcap         (2,j,ipat)
         veg_energy     (1,j,ipat) = veg_energy       (2,j,ipat)
         can_prss       (1,j,ipat) = can_prss         (2,j,ipat)
         can_theiv      (1,j,ipat) = can_theiv        (2,j,ipat)
         can_theta      (1,j,ipat) = can_theta        (2,j,ipat)
         can_rvap       (1,j,ipat) = can_rvap         (2,j,ipat)
         can_co2        (1,j,ipat) = can_co2          (2,j,ipat)
         sensible_gc    (1,j,ipat) = sensible_gc      (2,j,ipat)
         sensible_vc    (1,j,ipat) = sensible_vc      (2,j,ipat)
         evap_gc        (1,j,ipat) = evap_gc          (2,j,ipat)
         evap_vc        (1,j,ipat) = evap_vc          (2,j,ipat)
         transp         (1,j,ipat) = transp           (2,j,ipat)
         gpp            (1,j,ipat) = gpp              (2,j,ipat)
         plresp         (1,j,ipat) = plresp           (2,j,ipat)
         resphet        (1,j,ipat) = resphet          (2,j,ipat)
         veg_ndvip      (1,j,ipat) = veg_ndvip        (2,j,ipat)
         veg_ndvic      (1,j,ipat) = veg_ndvic        (2,j,ipat)
         veg_ndvif      (1,j,ipat) = veg_ndvif        (2,j,ipat)

         ustar         (m2,j,ipat) = ustar         (m2-1,j,ipat)
         tstar         (m2,j,ipat) = tstar         (m2-1,j,ipat)
         rstar         (m2,j,ipat) = rstar         (m2-1,j,ipat)
         cstar         (m2,j,ipat) = cstar         (m2-1,j,ipat)
         zeta          (m2,j,ipat) = zeta          (m2-1,j,ipat)
         ribulk        (m2,j,ipat) = ribulk        (m2-1,j,ipat)
         veg_albedo    (m2,j,ipat) = veg_albedo    (m2-1,j,ipat)
         veg_fracarea  (m2,j,ipat) = veg_fracarea  (m2-1,j,ipat)
         veg_lai       (m2,j,ipat) = veg_lai       (m2-1,j,ipat)
         veg_tai       (m2,j,ipat) = veg_tai       (m2-1,j,ipat)
         veg_rough     (m2,j,ipat) = veg_rough     (m2-1,j,ipat)
         veg_height    (m2,j,ipat) = veg_height    (m2-1,j,ipat)
         patch_area    (m2,j,ipat) = patch_area    (m2-1,j,ipat)
         patch_rough   (m2,j,ipat) = patch_rough   (m2-1,j,ipat)
         patch_wetind  (m2,j,ipat) = patch_wetind  (m2-1,j,ipat)
         leaf_class    (m2,j,ipat) = leaf_class    (m2-1,j,ipat)
         soil_rough    (m2,j,ipat) = soil_rough    (m2-1,j,ipat)
         sfcwater_nlev (m2,j,ipat) = sfcwater_nlev (m2-1,j,ipat)
         stom_resist   (m2,j,ipat) = stom_resist   (m2-1,j,ipat)
         ground_rsat   (m2,j,ipat) = ground_rsat   (m2-1,j,ipat)
         ground_rvap   (m2,j,ipat) = ground_rvap   (m2-1,j,ipat)
         ground_temp   (m2,j,ipat) = ground_temp   (m2-1,j,ipat)
         ground_fliq   (m2,j,ipat) = ground_fliq   (m2-1,j,ipat)
         veg_water     (m2,j,ipat) = veg_water     (m2-1,j,ipat)
         veg_hcap      (m2,j,ipat) = veg_hcap      (m2-1,j,ipat)
         veg_energy    (m2,j,ipat) = veg_energy    (m2-1,j,ipat)
         can_prss      (m2,j,ipat) = can_prss      (m2-1,j,ipat)
         can_theiv     (m2,j,ipat) = can_theiv     (m2-1,j,ipat)
         can_theta     (m2,j,ipat) = can_theta     (m2-1,j,ipat)
         can_rvap      (m2,j,ipat) = can_rvap      (m2-1,j,ipat)
         can_co2       (m2,j,ipat) = can_co2       (m2-1,j,ipat)
         sensible_gc   (m2,j,ipat) = sensible_gc   (m2-1,j,ipat)
         sensible_vc   (m2,j,ipat) = sensible_vc   (m2-1,j,ipat)
         evap_gc       (m2,j,ipat) = evap_gc       (m2-1,j,ipat)
         evap_vc       (m2,j,ipat) = evap_vc       (m2-1,j,ipat)
         transp        (m2,j,ipat) = transp        (m2-1,j,ipat)
         gpp           (m2,j,ipat) = gpp           (m2-1,j,ipat)
         plresp        (m2,j,ipat) = plresp        (m2-1,j,ipat)
         resphet       (m2,j,ipat) = resphet       (m2-1,j,ipat)
         veg_ndvip     (m2,j,ipat) = veg_ndvip     (m2-1,j,ipat)
         veg_ndvic     (m2,j,ipat) = veg_ndvic     (m2-1,j,ipat)
         veg_ndvif     (m2,j,ipat) = veg_ndvif     (m2-1,j,ipat)

         do k = 1,mzg
            soil_water       (k,1,j,ipat) = soil_water         (k,2,j,ipat)
            soil_energy      (k,1,j,ipat) = soil_energy        (k,2,j,ipat)
            soil_text        (k,1,j,ipat) = soil_text          (k,2,j,ipat)

            soil_water      (k,m2,j,ipat) = soil_water      (k,m2-1,j,ipat)
            soil_energy     (k,m2,j,ipat) = soil_energy     (k,m2-1,j,ipat)
            soil_text       (k,m2,j,ipat) = soil_text       (k,m2-1,j,ipat)
         end do

         do k = 1,mzs
            sfcwater_mass    (k,1,j,ipat) = sfcwater_mass      (k,2,j,ipat)
            sfcwater_energy  (k,1,j,ipat) = sfcwater_energy    (k,2,j,ipat)
            sfcwater_depth   (k,1,j,ipat) = sfcwater_depth     (k,2,j,ipat)

            sfcwater_mass   (k,m2,j,ipat) = sfcwater_mass   (k,m2-1,j,ipat)
            sfcwater_energy (k,m2,j,ipat) = sfcwater_energy (k,m2-1,j,ipat)
            sfcwater_depth  (k,m2,j,ipat) = sfcwater_depth  (k,m2-1,j,ipat)
         end do
      end do


      if (jdim == 1) then
         do i = 1,m2
            ustar          (i,1,ipat) = ustar            (i,2,ipat)
            tstar          (i,1,ipat) = tstar            (i,2,ipat)
            rstar          (i,1,ipat) = rstar            (i,2,ipat)
            cstar          (i,1,ipat) = cstar            (i,2,ipat)
            zeta           (i,1,ipat) = zeta             (i,2,ipat)
            ribulk         (i,1,ipat) = ribulk           (i,2,ipat)
            veg_albedo     (i,1,ipat) = veg_albedo       (i,2,ipat)
            veg_fracarea   (i,1,ipat) = veg_fracarea     (i,2,ipat)
            veg_lai        (i,1,ipat) = veg_lai          (i,2,ipat)
            veg_tai        (i,1,ipat) = veg_tai          (i,2,ipat)
            veg_rough      (i,1,ipat) = veg_rough        (i,2,ipat)
            veg_height     (i,1,ipat) = veg_height       (i,2,ipat)
            patch_area     (i,1,ipat) = patch_area       (i,2,ipat)
            patch_rough    (i,1,ipat) = patch_rough      (i,2,ipat)
            patch_wetind   (i,1,ipat) = patch_wetind     (i,2,ipat)
            leaf_class     (i,1,ipat) = leaf_class       (i,2,ipat)
            soil_rough     (i,1,ipat) = soil_rough       (i,2,ipat)
            sfcwater_nlev  (i,1,ipat) = sfcwater_nlev    (i,2,ipat)
            stom_resist    (i,1,ipat) = stom_resist      (i,2,ipat)
            ground_rsat    (i,1,ipat) = ground_rsat      (i,2,ipat)
            ground_rvap    (i,1,ipat) = ground_rvap      (i,2,ipat)
            ground_temp    (i,1,ipat) = ground_temp      (i,2,ipat)
            ground_fliq    (i,1,ipat) = ground_fliq      (i,2,ipat)
            veg_water      (i,1,ipat) = veg_water        (i,2,ipat)
            veg_hcap       (i,1,ipat) = veg_hcap         (i,2,ipat)
            veg_energy     (i,1,ipat) = veg_energy       (i,2,ipat)
            can_prss       (i,1,ipat) = can_prss         (i,2,ipat)
            can_theiv      (i,1,ipat) = can_theiv        (i,2,ipat)
            can_theta      (i,1,ipat) = can_theta        (i,2,ipat)
            can_rvap       (i,1,ipat) = can_rvap         (i,2,ipat)
            can_co2        (i,1,ipat) = can_co2          (i,2,ipat)
            sensible_gc    (i,1,ipat) = sensible_gc      (i,2,ipat)
            sensible_vc    (i,1,ipat) = sensible_vc      (i,2,ipat)
            evap_gc        (i,1,ipat) = evap_gc          (i,2,ipat)
            evap_vc        (i,1,ipat) = evap_vc          (i,2,ipat)
            transp         (i,1,ipat) = transp           (i,2,ipat)
            gpp            (i,1,ipat) = gpp              (i,2,ipat)
            plresp         (i,1,ipat) = plresp           (i,2,ipat)
            resphet        (i,1,ipat) = resphet          (i,2,ipat)
            veg_ndvip      (i,1,ipat) = veg_ndvip        (i,2,ipat)
            veg_ndvic      (i,1,ipat) = veg_ndvic        (i,2,ipat)
            veg_ndvif      (i,1,ipat) = veg_ndvif        (i,2,ipat)

            ustar         (i,m3,ipat) = ustar         (i,m3-1,ipat)
            tstar         (i,m3,ipat) = tstar         (i,m3-1,ipat)
            rstar         (i,m3,ipat) = rstar         (i,m3-1,ipat)
            cstar         (i,m3,ipat) = cstar         (i,m3-1,ipat)
            zeta          (i,m3,ipat) = zeta          (i,m3-1,ipat)
            ribulk        (i,m3,ipat) = ribulk        (i,m3-1,ipat)
            veg_albedo    (i,m3,ipat) = veg_albedo    (i,m3-1,ipat)
            veg_fracarea  (i,m3,ipat) = veg_fracarea  (i,m3-1,ipat)
            veg_lai       (i,m3,ipat) = veg_lai       (i,m3-1,ipat)
            veg_tai       (i,m3,ipat) = veg_tai       (i,m3-1,ipat)
            veg_rough     (i,m3,ipat) = veg_rough     (i,m3-1,ipat)
            veg_height    (i,m3,ipat) = veg_height    (i,m3-1,ipat)
            patch_area    (i,m3,ipat) = patch_area    (i,m3-1,ipat)
            patch_rough   (i,m3,ipat) = patch_rough   (i,m3-1,ipat)
            patch_wetind  (i,m3,ipat) = patch_wetind  (i,m3-1,ipat)
            leaf_class    (i,m3,ipat) = leaf_class    (i,m3-1,ipat)
            soil_rough    (i,m3,ipat) = soil_rough    (i,m3-1,ipat)
            sfcwater_nlev (i,m3,ipat) = sfcwater_nlev (i,m3-1,ipat)
            stom_resist   (i,m3,ipat) = stom_resist   (i,m3-1,ipat)
            ground_rsat   (i,m3,ipat) = ground_rsat   (i,m3-1,ipat)
            ground_rvap   (i,m3,ipat) = ground_rvap   (i,m3-1,ipat)
            ground_temp   (i,m3,ipat) = ground_temp   (i,m3-1,ipat)
            ground_fliq   (i,m3,ipat) = ground_fliq   (i,m3-1,ipat)
            veg_water     (i,m3,ipat) = veg_water     (i,m3-1,ipat)
            veg_hcap      (i,m3,ipat) = veg_hcap      (i,m3-1,ipat)
            veg_energy    (i,m3,ipat) = veg_energy    (i,m3-1,ipat)
            can_prss      (i,m3,ipat) = can_prss      (i,m3-1,ipat)
            can_theiv     (i,m3,ipat) = can_theiv     (i,m3-1,ipat)
            can_theta     (i,m3,ipat) = can_theta     (i,m3-1,ipat)
            can_rvap      (i,m3,ipat) = can_rvap      (i,m3-1,ipat)
            can_co2       (i,m3,ipat) = can_co2       (i,m3-1,ipat)
            sensible_gc   (i,m3,ipat) = sensible_gc   (i,m3-1,ipat)
            sensible_vc   (i,m3,ipat) = sensible_vc   (i,m3-1,ipat)
            evap_gc       (i,m3,ipat) = evap_gc       (i,m3-1,ipat)
            evap_vc       (i,m3,ipat) = evap_vc       (i,m3-1,ipat)
            transp        (i,m3,ipat) = can_co2       (i,m3-1,ipat)
            gpp           (i,m3,ipat) = gpp           (i,m3-1,ipat)
            plresp        (i,m3,ipat) = plresp        (i,m3-1,ipat)
            resphet       (i,m3,ipat) = resphet       (i,m3-1,ipat)
            veg_ndvip     (i,m3,ipat) = veg_ndvip     (i,m3-1,ipat)
            veg_ndvic     (i,m3,ipat) = veg_ndvic     (i,m3-1,ipat)
            veg_ndvif     (i,m3,ipat) = veg_ndvif     (i,m3-1,ipat)

            do k = 1,mzg
               soil_water       (k,i,1,ipat) = soil_water         (k,i,2,ipat)
               soil_energy      (k,i,1,ipat) = soil_energy        (k,i,2,ipat)
               soil_text        (k,i,1,ipat) = soil_text          (k,i,2,ipat)

               soil_water      (k,i,m3,ipat) = soil_water      (k,i,m3-1,ipat)
               soil_energy     (k,i,m3,ipat) = soil_energy     (k,i,m3-1,ipat)
               soil_text       (k,i,m3,ipat) = soil_text       (k,i,m3-1,ipat)
            end do

            do k = 1,mzs
               sfcwater_mass    (k,i,1,ipat) = sfcwater_mass      (k,i,2,ipat)
               sfcwater_energy  (k,i,1,ipat) = sfcwater_energy    (k,i,2,ipat)
               sfcwater_depth   (k,i,1,ipat) = sfcwater_depth     (k,i,2,ipat)

               sfcwater_mass   (k,i,m3,ipat) = sfcwater_mass   (k,i,m3-1,ipat)
               sfcwater_energy (k,i,m3,ipat) = sfcwater_energy (k,i,m3-1,ipat)
               sfcwater_depth  (k,i,m3,ipat) = sfcwater_depth  (k,i,m3-1,ipat)
            end do
         end do
      end if

   end do
   return
end subroutine leaf_bcond
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine sfc_fields(m1,m2,m3,ia,iz,ja,jz,jd,thp,theta,rv,rtp,co2p,up,vp,dn0,pp,pi0,rtgt  &
                     ,zt,zm,thils2,ths2,rvs2,rtps2,co2s2,ups2,vps2,pis2,dens2,zts2)

   use leaf_coms
   use rconstants

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: m1
   integer                  , intent(in)    :: m2
   integer                  , intent(in)    :: m3
   integer                  , intent(in)    :: ia
   integer                  , intent(in)    :: iz
   integer                  , intent(in)    :: ja
   integer                  , intent(in)    :: jz
   integer                  , intent(in)    :: jd
   real, dimension(m1,m2,m3), intent(in)    :: thp
   real, dimension(m1,m2,m3), intent(in)    :: theta
   real, dimension(m1,m2,m3), intent(in)    :: rv
   real, dimension(m1,m2,m3), intent(in)    :: rtp
   real, dimension(m1,m2,m3), intent(in)    :: co2p
   real, dimension(m1,m2,m3), intent(in)    :: up
   real, dimension(m1,m2,m3), intent(in)    :: vp
   real, dimension(m1,m2,m3), intent(in)    :: dn0
   real, dimension(m1,m2,m3), intent(in)    :: pp
   real, dimension(m1,m2,m3), intent(in)    :: pi0
   real, dimension(m1)      , intent(in)    :: zt
   real, dimension(m1)      , intent(in)    :: zm
   real, dimension(m2,m3)   , intent(in)    :: rtgt
   real, dimension(m2,m3)   , intent(inout) :: thils2
   real, dimension(m2,m3)   , intent(inout) :: ths2
   real, dimension(m2,m3)   , intent(inout) :: rvs2
   real, dimension(m2,m3)   , intent(inout) :: rtps2
   real, dimension(m2,m3)   , intent(inout) :: co2s2
   real, dimension(m2,m3)   , intent(inout) :: ups2
   real, dimension(m2,m3)   , intent(inout) :: vps2
   real, dimension(m2,m3)   , intent(inout) :: pis2
   real, dimension(m2,m3)   , intent(inout) :: dens2
   real, dimension(m2,m3)   , intent(inout) :: zts2
   !----- Local variables. ----------------------------------------------------------------!
   integer                                  :: i
   integer                                  :: j
   !---------------------------------------------------------------------------------------!



   !----- Compute surface atmospheric conditions. -----------------------------------------!
   do j = ja,jz
      do i = ia,iz
         thils2(i,j) = thp  (2,i,j)
         ths2  (i,j) = theta(2,i,j)
         rvs2  (i,j) = rv   (2,i,j)
         rtps2 (i,j) = rtp  (2,i,j)
         co2s2 (i,j) = co2p (2,i,j)
         ups2  (i,j) = (up(2,i-1,j) + up(2,i,j))  * .5
         vps2  (i,j) = (vp(2,i,j-jd) + vp(2,i,j)) * .5
         zts2  (i,j) = (zt(2)-zm(1)) * rtgt(i,j)
         pis2  (i,j) = 0.5 * (pp(1,i,j) + pi0(1,i,j) + pp(2,i,j) + pi0(2,i,j))
         dens2 (i,j) = (dn0(1,i,j) + dn0(2,i,j))  * .5
      end do
   end do

   return
end subroutine sfc_fields
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine sfc_fields_adap(m1,m2,m3,ia,iz,ja,jz,jd,flpu,flpv,flpw,topma,aru,arv,thp,theta  &
                          ,rv,rtp,co2p,up,vp,dn0,pp,pi0,zt,zm,dzt,thils2,ths2,rvs2,rtps2   &
                          ,co2s2,ups2,vps2,pis2,dens2,zts2)

   use leaf_coms
   use rconstants

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: m1
   integer                  , intent(in)    :: m2
   integer                  , intent(in)    :: m3
   integer                  , intent(in)    :: ia
   integer                  , intent(in)    :: iz
   integer                  , intent(in)    :: ja
   integer                  , intent(in)    :: jz
   integer                  , intent(in)    :: jd
   real, dimension(m1,m2,m3), intent(in)    :: aru
   real, dimension(m1,m2,m3), intent(in)    :: arv
   real, dimension(m1,m2,m3), intent(in)    :: thp
   real, dimension(m1,m2,m3), intent(in)    :: theta
   real, dimension(m1,m2,m3), intent(in)    :: rv
   real, dimension(m1,m2,m3), intent(in)    :: rtp
   real, dimension(m1,m2,m3), intent(in)    :: co2p
   real, dimension(m1,m2,m3), intent(in)    :: up
   real, dimension(m1,m2,m3), intent(in)    :: vp
   real, dimension(m1,m2,m3), intent(in)    :: dn0
   real, dimension(m1,m2,m3), intent(in)    :: pp
   real, dimension(m1,m2,m3), intent(in)    :: pi0
   real, dimension(m1)      , intent(in)    :: zt
   real, dimension(m1)      , intent(in)    :: zm
   real, dimension(m1)      , intent(in)    :: dzt
   real, dimension(m2,m3)   , intent(in)    :: flpu
   real, dimension(m2,m3)   , intent(in)    :: flpv
   real, dimension(m2,m3)   , intent(in)    :: flpw
   real, dimension(m2,m3)   , intent(in)    :: topma
   real, dimension(m2,m3)   , intent(inout) :: thils2
   real, dimension(m2,m3)   , intent(inout) :: ths2
   real, dimension(m2,m3)   , intent(inout) :: rvs2
   real, dimension(m2,m3)   , intent(inout) :: rtps2
   real, dimension(m2,m3)   , intent(inout) :: co2s2
   real, dimension(m2,m3)   , intent(inout) :: ups2
   real, dimension(m2,m3)   , intent(inout) :: vps2
   real, dimension(m2,m3)   , intent(inout) :: pis2
   real, dimension(m2,m3)   , intent(inout) :: dens2
   real, dimension(m2,m3)   , intent(inout) :: zts2
   !----- Local variables. ----------------------------------------------------------------!
   integer                                  :: i
   integer                                  :: j
   integer                                  :: k1
   integer                                  :: k2
   integer                                  :: k3
   real                                     :: topma_t
   real                                     :: wtw
   real                                     :: wtu1
   real                                     :: wtu2
   real                                     :: wtv1
   real                                     :: wtv2
   !---------------------------------------------------------------------------------------!


   !----- Compute surface atmospheric conditions. -----------------------------------------!

   do j = ja,jz
      do i = ia,iz
         k2 = nint(flpw(i,j))
         k1 = k2 - 1
         k3 = k2 + 1

         topma_t = .25 * (topma(i,j) + topma(i-1,j) + topma(i,j-jd) + topma(i-1,j-jd))

         !----- Weights for lowest predicted points, relative to points above them. -------!
         wtw  = (zm(k2) - topma_t) * dzt(k2)
         wtu1 = aru(nint(flpu(i-1,j)),i-1,j)   / aru(nint(flpu(i-1,j))+1,i-1,j)
         wtu2 = aru(nint(flpu(i,j)),i,j)       / aru(nint(flpu(i,j))+1,i,j)
         wtv1 = arv(nint(flpv(i,j-jd)),i,j-jd) / arv(nint(flpv(i,j-jd))+1,i,j-jd)
         wtv2 = arv(nint(flpv(i,j)),i,j)       / arv(nint(flpv(i,j))+1,i,j)

         !----- Interpolate the values to with height. ------------------------------------!
         thils2(i,j) =  wtw * thp  (k2,i,j) + (1. - wtw)  * thp  (k3,i,j)
         ths2  (i,j) =  wtw * theta(k2,i,j) + (1. - wtw)  * theta(k3,i,j)
         rvs2  (i,j) =  wtw * rv   (k2,i,j) + (1. - wtw)  * rv   (k3,i,j)
         rtps2 (i,j) =  wtw * rtp  (k2,i,j) + (1. - wtw)  * rtp  (k3,i,j)

         co2s2(i,j) =  wtw * co2p(k2,i,j)  + (1. - wtw)  * co2p(k3,i,j)

         ups2(i,j) = (wtu1        * up(nint(flpu(i-1,j)),i-1,j)                            &
                   +  (1. - wtu1) * up(nint(flpu(i-1,j))+1,i-1,j)                          &
                   +  wtu2        * up(nint(flpu(i,j)),i,j)                                &
                   +  (1. - wtu2) * up(nint(flpu(i,j))+1,i,j)) * .5

         vps2(i,j) = (wtv1        * vp(nint(flpv(i,j-jd)),i,j-jd)                          &
                   +  (1. - wtv1) * vp(nint(flpv(i,j-jd))+1,i,j-jd)                        &
                   +  wtv2        * vp(nint(flpv(i,j)),i,j)                                &
                   +  (1. - wtv2) * vp(nint(flpv(i,j))+1,i,j)) * .5

         zts2(i,j) = (wtw * (zt(k2) - zm(k1)) + (1. - wtw) * (zt(k3) - zm(k2)))

         if (wtw >= .5) then
            pis2(i,j)  = ( (wtw - .5) * (pp(k1,i,j) + pi0(k1,i,j))                         &
                         + (1.5 - wtw) * (pp(k2,i,j) + pi0(k2,i,j)))
            dens2(i,j) = (wtw - .5)  * dn0(k1,i,j) + (1.5 - wtw) * dn0(k2,i,j)
         else
            pis2(i,j)  = ( (wtw + .5) * (pp(k2,i,j) + pi0(k2,i,j))                         &
                         + (.5 - wtw) * (pp(k3,i,j) + pi0(k3,i,j)))
            dens2(i,j) = (wtw + .5) * dn0(k2,i,j) + (.5 - wtw) * dn0(k3,i,j)
         end if

      end do
   end do

   return
end subroutine sfc_fields_adap
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine sfc_pcp(nqparm,i,j,cuparm,micro)
   use mem_micro , only : micro_vars  ! ! structure
   use mem_cuparm, only : cuparm_vars & ! structure
                        , nclouds     ! ! intent(in)
   use leaf_coms , only : pcpgl       & ! intent(in)
                        , qpcpgl      & ! intent(in)
                        , dpcpgl      & ! intent(in)
                        , dtll        & ! intent(in)
                        , dtll_factor & ! intent(in)
                        , atm_temp    ! ! intent(in)
   use therm_lib , only : bulk_on     ! ! intent(in)
   use rconstants, only : cice        & ! intent(in)
                        , cliq        & ! intent(in)
                        , tsupercool  & ! intent(in)
                        , t3ple       & ! intent(in)
                        , wdnsi       ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer           , intent(in) :: nqparm
   integer           , intent(in) :: i
   integer           , intent(in) :: j
   type (cuparm_vars), intent(in) :: cuparm
   type (micro_vars) , intent(in) :: micro
   !----- Local variables. ----------------------------------------------------------------!
   integer                        :: icld
   real                           :: pcpgcum
   real                           :: fice
   real                           :: sndeni
   !---------------------------------------------------------------------------------------!

   !----- Initialise the precipitation variables. -----------------------------------------!
   pcpgl  = 0.
   qpcpgl = 0.
   dpcpgl = 0.

   !----- Add cumulus parametrisation precipitation if it is being used. ------------------!
   if (nqparm > 0) then

      !----- First find the total precipitation. ------------------------------------------!
      pcpgcum = 0.
      do icld=1,nclouds
         pcpgcum = pcpgcum + cuparm%conprr(i,j,icld) * dtll
      end do
      pcpgl = pcpgl + pcpgcum
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !  Precipitation "depth". Snow fraction and density derived from                     !
      !  Jin et al 1999 Hydrol Process. 13:2467-2482 Table 2                               !
      !  [[modified 11/16/09 by MCD]]                                                      !
      !------------------------------------------------------------------------------------!
      if (atm_temp > (t3ple + 2.5)) then
         !----- Rain only. ----------------------------------------------------------------!
         fice    = 0.0
         sndeni  = 1. / 189.0

      elseif (atm_temp <= (t3ple + 2.5) .and. atm_temp  > (t3ple + 2.0) ) then
         !---------------------------------------------------------------------------------!
         !     60% snow, 40% rain. (N.B. May not be appropriate for sub-tropical           !
         ! regions where the level of the melting layer is higher...).                     !
         !---------------------------------------------------------------------------------!
         fice    = 0.6
         sndeni  = 1. / 189.0

      elseif (atm_temp <= (t3ple + 2.0) .and. atm_temp > t3ple ) then
         !---------------------------------------------------------------------------------!
         !     Increasing the fraction of snow. (N.B. May not be appropriate for           !
         ! sub-tropical regions where the level of the melting layer is higher...).        !
         !---------------------------------------------------------------------------------!
         fice   = min(1.0, 1. + (54.62 - 0.2*atm_temp))
         sndeni = 1. / (50.0+1.7*(atm_temp-258.15)**1.5 )

      elseif (atm_temp <= t3ple .and. atm_temp > (t3ple - 15.0)) then
         !----- Below freezing point, snow only. ------------------------------------------!
         fice   = 1.0
         sndeni = 1. / (50.0+1.7*(atm_temp-258.15)**1.5 )

      else ! if (atm_temp < (t3ple - 15.0)) then
         !----- Below freezing point, snow only. ------------------------------------------!
         fice   = 1.0
         sndeni = 1. / 50.
      end if
      dpcpgl  = dpcpgl + pcpgcum * ((1.0-fice) * wdnsi + fice * sndeni)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Set internal energy.  This will be the precipitation times the specific        !
      ! internal energy of water (above or at triple point) multiplied by the liquid       !
      ! fraction plus the specific internal energy of ice (below or at the triple point)   !
      ! multiplied by the ice fraction.                                                    !
      !------------------------------------------------------------------------------------!
      qpcpgl = qpcpgl                                                                      &
             + pcpgcum  * ( (1.0-fice) * cliq * ( max(t3ple,atm_temp) - tsupercool)        &
                          +      fice  * cice *   min(atm_temp,t3ple)               )
      !------------------------------------------------------------------------------------!
   end if

   !----- Add microphysics precipitation if the bulk microphysics is being used. ----------!
   if (bulk_on) then
      pcpgl  = pcpgl  + dtll_factor * micro%pcpg(i,j)
      qpcpgl = qpcpgl + dtll_factor * micro%qpcpg(i,j)
      dpcpgl = dpcpgl + dtll_factor * micro%dpcpg(i,j)

   end if

   return
end subroutine sfc_pcp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the NDVI-related variables such as LAI, TAI, roughness,     !
! and albedo.                                                                              !
!------------------------------------------------------------------------------------------!
subroutine vegndvi(ifm,patch_area,leaf_class,veg_fracarea,veg_lai,veg_tai,veg_rough        &
                  ,veg_height,veg_albedo,veg_ndvip,veg_ndvic,veg_ndvif)

   use leaf_coms
   use rconstants
   use io_params
   use mem_grid
   use catt_start, only: catt  ! INTENT(IN)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                         , intent(in)    :: ifm
   real                            , intent(in)    :: patch_area
   real                            , intent(in)    :: leaf_class
   real                            , intent(in)    :: veg_height
   real                            , intent(in)    :: veg_ndvip
   real                            , intent(in)    :: veg_ndvif
   real                            , intent(out)   :: veg_lai
   real                            , intent(out)   :: veg_tai
   real                            , intent(out)   :: veg_fracarea
   real                            , intent(out)   :: veg_rough
   real                            , intent(out)   :: veg_albedo
   real                            , intent(inout) :: veg_ndvic
   !----- Local variables. ----------------------------------------------------------------!
   integer                                       :: nveg
   real                                          :: timefac_ndvi
   real                                          :: sr
   real                                          :: fpar
   real                                          :: dead_lai
   real                                          :: green_frac
   !----- Local constants. ----------------------------------------------------------------!
   real                            , parameter   :: sr_min=1.081
   real                            , parameter   :: fpar_min=.001
   real                            , parameter   :: fpar_max=.950
   real                            , parameter   :: fpcon=-.3338082
   real                            , parameter   :: ccc=-2.9657
   real                            , parameter   :: bz=.91
   real                            , parameter   :: hz=.0075
   real                            , parameter   :: extinc_veg=.5
   !----- Locally saved variables. --------------------------------------------------------!
   logical                         , save        :: nvcall = .true.
   real, dimension(nvtyp+nvtyp_teb), save        :: dfpardsr
   !---------------------------------------------------------------------------------------!


   !----- Initialise the dfpardsr array, which will be used to compute LAI. ---------------!
   if (nvcall) then
      nvcall = .false.
      do nveg = 1,(nvtyp+nvtyp_teb)
         dfpardsr(nveg) = (fpar_max - fpar_min) / (sr_max(nveg) - sr_min)
      end do
   end if

   !----- Alias for vegetation class. -----------------------------------------------------!
   nveg = nint(leaf_class)

   !---------------------------------------------------------------------------------------!
   !    We only compute LAI and related variables for those vegetation types that can hold !
   ! some actual vegetation.                                                               !
   !---------------------------------------------------------------------------------------!
   if (tai_max(nveg) < .1) then
      veg_lai      = 0.
      veg_tai      = 0.
      veg_rough    = 0.
      veg_albedo   = 0.
      veg_fracarea = 0.

   else
      !----- Compute the interpolation factor. --------------------------------------------!
      if (iupdndvi == 0) then
         timefac_ndvi = 0.
      else
         timefac_ndvi = sngl((time - ndvitime1(ifm)) / (ndvitime2(ifm) - ndvitime1(ifm)))
      end if

      !  Time-interpolate ndvi to get current value veg_ndvic(i,j) for this patch
      !  Limit ndvi to prevent values > .99 to prevent division by zero.

      

      if (iuselai == 1) then
         veg_ndvic = max(0.0,veg_ndvip + (veg_ndvif - veg_ndvip) * timefac_ndvi)
      else
         !---------------------------------------------------------------------------------!
         !     Time-interpolate ndvi to get current value veg_ndvic(i,j) for this patch.   !
         ! Limit NDVI to prevent values > .99 to prevent division by zero.                 !
         !---------------------------------------------------------------------------------!
         veg_ndvic = max(0.99,veg_ndvip + (veg_ndvif - veg_ndvip) * timefac_ndvi)


         !---------------------------------------------------------------------------------!
         ! Based on SRF's suggestion - MODIS sometimes has weird values of NDVI for        !
         !                             evergreen forests, so we impose a lower threshold   !
         !                             for this vegetation type.                           !
         ! (Maybe we should try using something better than NDVI, perhaps EVI?)            !
         !---------------------------------------------------------------------------------!
         if (nveg == 7) veg_ndvic = max(0.7,veg_ndvic)
         
      end if


      if (iuselai == 1) then

         veg_lai = veg_ndvic

      else

         !----- Compute "simple ratio" and limit between sr_min and sr_max(nveg). ---------!
         sr = min(sr_max(nveg), max(sr_min, (1. + veg_ndvic) / (1. - veg_ndvic) ) )


         !----- Compute fpar. -------------------------------------------------------------!
         fpar = fpar_min + (sr - sr_min) * dfpardsr(nveg)

         !---------------------------------------------------------------------------------!
         !      Compute green leaf area index (veg_lai), dead leaf area index (dead_lai),  !
         ! total area index (tai), and green fraction.                                     !
         !---------------------------------------------------------------------------------!
         veg_lai    = glai_max(nveg) * (       veg_clump(nveg)  * fpar / fpar_max          &
                                       + (1. - veg_clump(nveg)) * alog(1. - fpar) * fpcon )

      end if
         
         
      dead_lai   = (glai_max(nveg) - veg_lai) * dead_frac(nveg)
      veg_tai    = veg_lai + sai(nveg) + dead_lai
      green_frac = veg_lai / veg_tai

      !----- Compute vegetation roughness height, albedo, and fractional area. ------------!
      veg_rough    = veg_height * (1. - bz * exp(-hz * veg_tai))
      veg_albedo   = albv_green(nveg) * green_frac + albv_brown(nveg) * (1. - green_frac)
      veg_fracarea = veg_frac(nveg) * (1. - exp(-extinc_veg * veg_tai))
   end if

   return
end subroutine vegndvi
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This routine is called by the radiation parameterization and by LEAF.  It computes   !
! net surface albedo plus radiative exchange between the atmosphere, vegetation, and the   !
! snow/ground given previously computed downward longwave and shortwave fluxes from the    !
! atmosphere.  Also computed are functions of snowcover that are required for the above    !
! radiation calculations as well as other calculations in LEAF.                            !
!     The shortwave parameterizations are only valid if the cosine of the zenith angle is  !
! greater than .03 .  Water albedo from Atwater and Bell (1981) alg, als, and alv are the  !
! albedos of the ground, snow, and vegetation (als needs a better formula based on age of  !
! the surface snow).  absg and vctr32 are the actual fractions of shortwave incident on    !
! snow plus ground that get absorbed by the ground and each snow layer, respectively.      !
! They currently use the variable fractrans, which is the fraction of light transmitted    !
! through each layer based on mass per square meter.  algs is the resultant albedo from    !
! snow plus ground.                                                                        !
!------------------------------------------------------------------------------------------!
subroutine sfcrad(mzg,mzs,ip,soil_energy,soil_water,soil_text,sfcwater_energy              &
                 ,sfcwater_mass,sfcwater_depth,patch_area,leaf_class,veg_height            &
                 ,veg_fracarea,veg_albedo,sfcwater_nlev,veg_energy,veg_water,veg_hcap      &
                 ,can_prss,can_theiv,can_theta,can_rvap,rshort,rlong,albedt,rlongup,cosz   &
                 ,g_urban, etown, albtown, tstown)
   use mem_leaf
   use leaf_coms
   use rconstants
   use mem_scratch
   use node_mod     , only : mynum         ! ! intent(in)
   use therm_lib    , only : qwtk          & ! subroutine
                           , qtk           & ! subroutine
                           , idealdenssh   ! ! function
   use catt_start   , only : CATT          ! ! intent(in)
   use teb_spm_start, only : TEB_SPM       ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                , intent(in)    :: mzg,mzs,ip
   real   , dimension(mzg), intent(in)    :: soil_energy,soil_water,soil_text
   real   , dimension(mzs), intent(in)    :: sfcwater_energy,sfcwater_depth ,sfcwater_mass
   real                   , intent(in)    :: patch_area,leaf_class,veg_height,veg_fracarea
   real                   , intent(in)    :: veg_albedo,sfcwater_nlev
   real                   , intent(in)    :: veg_energy,veg_water,veg_hcap
   real                   , intent(in)    :: can_prss,can_theiv,can_theta,can_rvap
   real                   , intent(in)    :: rshort,rlong,cosz
   real                   , intent(in)    :: g_urban, etown, albtown, tstown
   real                   , intent(inout) :: albedt,rlongup
   !----- Local variables. ----------------------------------------------------------------!
   integer                                :: k,m,nsoil,nveg,ksn
   real                                   :: alb,vfc,fcpct,alg,rad,als,fractrans
   real                                   :: absg,algs,emv,emgs,gslong,vlong,alv
   !----- Local constants. ----------------------------------------------------------------!
   character(len=9)      , parameter   :: fmti='(a,1x,i6)'
   character(len=13)     , parameter   :: fmtf='(a,1x,es12.5)'
   character(len=3)      , parameter   :: fmtc='(a)'
   character(len=9)      , parameter   :: fmtl='(a,1x,l1)'
   character(len=9)      , parameter   :: fmth='(7(a,1x))'
   character(len=23)     , parameter   :: fmts='(2(i5,1x),5(es12.5,1x))'
   character(len=9)      , parameter   :: fmte='(5(a,1x))'
   character(len=23)     , parameter   :: fmtw='(1(i5,1x),4(es12.5,1x))'
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     First we update the canopy air properties.                                        !
   !---------------------------------------------------------------------------------------!
   can_exner    = cp  * (p00i * can_prss) ** rocp
   can_lntheta  = log(can_theta)
   can_temp     = cpi * can_theta * can_exner
   can_shv      = can_rvap / (can_rvap + 1.)
   can_rhos     = idealdenssh(can_prss,can_temp,can_shv)

   if (ip == 1) then
      !----- Compute the albedo and upward longwave for water patches. --------------------!
      if (cosz > .03) then
         alb = min(max(-.0139 + .0467 * tan(acos(cosz)),.03),.999)
         albedt = albedt + patch_area * alb
      end if
      call qtk(soil_energy(mzg),tempk(mzg),fracliq(mzg))
      rlongup = rlongup + patch_area * stefan * tempk(mzg) ** 4

   elseif (isfcl == 0) then
      !------ Not running a land surface model, use prescribed value of can_temp. ---------!
      albedt  = albedt  + patch_area * albedo
      rlongup = rlongup + patch_area * stefan * can_temp ** 4
   else
      !------ Running an actual land surface model... -------------------------------------!


      !------ Diagnose vegetation temperature and surface water liquid fraction. ----------!
      call qwtk(veg_energy,veg_water,veg_hcap,veg_temp,veg_fliq)

      !------ Diagnose soil temperature and liquid fraction. ------------------------------!
      do k = 1,mzg
         nsoil = nint(soil_text(k))
         call qwtk(soil_energy(k),soil_water(k)*wdns,slcpd(nsoil),tempk(k),fracliq(k))
      end do

      !------ Diagnose snow temperature and the influence of snow covering veg. -----------!
      nveg = nint(leaf_class)
      ksn  = nint(sfcwater_nlev)
      snowfac = 0.
      do k = 1,ksn
         snowfac = snowfac + sfcwater_depth(k)
         call qtk(sfcwater_energy(k),tempk(k+mzg),fracliq(k+mzg))
      end do
      snowfac = min(.99, snowfac / max(.001,veg_height))

      !------ Defining the exposed area. --------------------------------------------------!
      vf = veg_fracarea * (1. - snowfac)
      vfc = 1. - vf

      !------------------------------------------------------------------------------------!
      !     Shortwave radiation calculations.                                              !
      !------------------------------------------------------------------------------------!
      nsoil=nint(soil_text(mzg))
      fcpct = soil_water(mzg) / slmsts(nsoil)
      if (fcpct > .5) then
         alg = .14
      else
         alg = .31 - .34 * fcpct
      end if
      alv = veg_albedo

      rad = 1.
      if (ksn > 0) then
         !------ als = .14 (the wet soil value) for all-liquid. ---------------------------!
         als = .5 - .36 * fracliq(ksn+mzg)
         rad = 1. - als
      end if
      do k = ksn,1,-1
         fractrans = exp(-20. * sfcwater_depth(k))
         vctr32(k) = rad * (1. - fractrans)
         rad = rad * fractrans
      end do
      absg = (1. - alg) * rad
      algs = 1. - absg
      do k = ksn,1,-1
         algs = algs - vctr32(k)
         rshort_s(k) = rshort * vfc * vctr32(k)
      end do
      rshort_g = rshort * vfc * absg
      rshort_v = rshort * vf * (1. - alv + vfc * algs)
      alb      = vf * alv + vfc * vfc * algs

      !----- Adding urban contribution if running TEB. ------------------------------------!
      if (teb_spm==1) then
         if (nint(g_urban) == 0) then
            albedt = albedt + patch_area * alb
         else
            albedt = albedt + patch_area * albtown
         endif
      else
         albedt = albedt + patch_area * alb
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Longwave radiation calculations.                                               !
      !------------------------------------------------------------------------------------!
      emv  = emisv(nveg)
      emgs = emisg(nsoil)
      if (ksn > 0) emgs = 1.0
      gslong = emgs * stefan * tempk(ksn+mzg) ** 4
      vlong  = emv * stefan * veg_temp ** 4

      rlonga_v  = rlong  * vf * (emv + vfc * (1. - emgs))
      rlonga_gs = rlong  * vfc * emgs
      rlongv_gs = vlong  * vf * emgs
      rlongv_a  = vlong  * vf * (2. - emgs - vf + emgs * vf)
      rlonggs_v = gslong * vf * emv
      rlonggs_a = gslong * vfc
      rlonga_a  = rlong  * (vf * (1. - emv) + vfc * vfc * (1. - emgs))
      
      !----- Sanity check. ----------------------------------------------------------------!
      if (rlonga_v /= rlonga_v) then
         write (unit=*,fmt=fmtc) '------------ Longwave radiation is screwed. ------------'
         write (unit=*,fmt=fmti) ' - PATCH        = ',ip
         write (unit=*,fmt=fmti) ' - LEAF_CLASS   = ',nveg
         write (unit=*,fmt=fmti) ' - KSN          = ',ksn
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - RLONGA_V     = ',rlonga_v
         write (unit=*,fmt=fmtf) ' - RLONG        = ',rlong
         write (unit=*,fmt=fmtf) ' - EMV          = ',emv
         write (unit=*,fmt=fmtf) ' - VF           = ',vf
         write (unit=*,fmt=fmtf) ' - VFC          = ',vfc
         write (unit=*,fmt=fmtf) ' - EMGS         = ',emgs
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - CAN_THETA    = ',can_theta
         write (unit=*,fmt=fmtf) ' - CAN_RVAP     = ',can_rvap
         write (unit=*,fmt=fmtf) ' - CAN_PRSS     = ',can_prss
         write (unit=*,fmt=fmtf) ' - CAN_THEIV    = ',can_theiv
         write (unit=*,fmt=fmtf) ' - CAN_TEMP     = ',can_temp
         write (unit=*,fmt=fmtf) ' - CAN_RHOS     = ',can_rhos
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - VEG_ENERGY   = ',veg_energy
         write (unit=*,fmt=fmtf) ' - VEG_WATER    = ',veg_water
         write (unit=*,fmt=fmtf) ' - VEG_HCAP     = ',veg_hcap
         write (unit=*,fmt=fmtf) ' - VEG_TEMP     = ',veg_temp
         write (unit=*,fmt=fmtf) ' - VEG_FLIQ     = ',veg_fliq
         write (unit=*,fmt=fmtc) ' '
         if (ip == 1) then
            write (unit=*,fmt=fmtf) ' - SST_ENERGY   = ',soil_energy(mzg)
            write (unit=*,fmt=fmtf) ' - SST_TEMP     = ',tempk(mzg)
            write (unit=*,fmt=fmtf) ' - SST_FLIQ     = ',fracliq(mzg)
         write (unit=*,fmt=fmtc) ' '
         else
            write (unit=*,fmt=fmth) '    K','NSOIL',' SOIL_ENERGY','  SOIL_WATER'         &
                                                   ,'   SOIL_HCAP','   SOIL_TEMP'         &
                                                   ,'   SOIL_FLIQ'
            do k=1,mzg
               nsoil = nint(soil_text(k))
               write(unit=*,fmt=fmts) k,nsoil,soil_energy(k),soil_water(k),slcpd(nsoil)   &
                                             ,tempk(k),fracliq(k)
            end do
            write (unit=*,fmt=fmtc) ' '
            if (ksn > 0) then
               write (unit=*,fmt=fmte) '    K',' SFCW_ENERGY','   SFCW_MASS'              &
                                              ,'   SFCW_TEMP','   SFCW_FLIQ'
               do k=1,mzg
                  nsoil = nint(soil_text(k))
                  write(unit=*,fmt=fmtw) k,sfcwater_energy(k),sfcwater_mass(k)            &
                                          ,tempk(k+mzg),fracliq(k+mzg)
               end do
               write (unit=*,fmt=fmtc) ' '
            end if
         end if
         write (unit=*,fmt=fmtc) '-----------------------------------------------------'
         
      end if

      !----- Adding urban contribution if running TEB. ------------------------------------!
      if (teb_spm==1) then
         if (nint(g_urban) == 0) then
            rlongup = rlongup + patch_area * (rlongv_a + rlonggs_a + rlonga_a)
         else
            rlongup = rlongup + patch_area * etown * stefan * tstown**4
         endif
      else
         rlongup = rlongup + patch_area * (rlongv_a + rlonggs_a + rlonga_a)
      endif

      !------------------------------------------------------------------------------------!
      !      In case rlong is not computed, zero out all longwave fluxes other than        !
      ! rlongup.  [On the first timestep, radiative fluxes may not be available until      !
      ! microphysics is called, and zeroing these fluxes avoids the imbalance of having    !
      ! upward longwave without downward longwave in LEAF-3.  Also, this allows LEAF-3 to  !
      ! run without radiation for all timesteps, if desired for testing.].                 !
      !------------------------------------------------------------------------------------!
      if (rlong < .1) then
         rlonga_v  = 0.
         rlonga_gs = 0.
         rlongv_gs = 0.
         rlongv_a  = 0.
         rlonggs_v = 0.
         rlonggs_a = 0.
         rlonga_a  = 0.
      end if
      !------------------------------------------------------------------------------------!
   end if

   return
end subroutine sfcrad
!==========================================================================================!
!==========================================================================================!

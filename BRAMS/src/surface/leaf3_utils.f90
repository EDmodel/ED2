!==========================================================================================!
!==========================================================================================!
!                                                                                          !
!     This subroutine computes the characteristic scales, using one of the following surf- !
! ace layer parametrisation.                                                               !
!                                                                                          !
!                                                                                          !
! 1. Based on L79;                                                                         !
! 2. Based on: OD95, but with some terms computed as in L79 and B71 to avoid singular-     !
!    ities (now using the iterative method to find zeta).                                  !
! 3. Based on BH91, using an iterative method to find zeta, and using the modified         !
!    equation for stable layers.                                                           !
! 4. Based on CLM04, with special functions for very stable and very stable case, even     !
!    though we use a different functional form for very unstable case for momentum.        !
!    This is ensure that phi_m decreases monotonically as zeta becomes more negative.      !
!    We use a power law of order of -1/6 instead.                                          !
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
! CLM04. OLESON, K. W., et al.; Technical description of the community land model (CLM)    !
!           NCAR Technical Note NCAR/TN-461+STR, Boulder, CO, May 2004.                    !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine leaf3_stars(theta_atm,enthalpy_atm,shv_atm,rvap_atm,co2_atm                     &
                      ,theta_can,enthalpy_can,shv_can,rvap_can,co2_can                     &
                      ,zref,dheight,uref,dtl3,rough,ustar,tstar,estar,qstar,rstar,cstar    &
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
                        , ribmax     & ! intent(in)
                        , tprandtl   & ! intent(in)
                        , z0moz0h    & ! intent(in)
                        , z0hoz0m    & ! intent(in)
                        , ggbare     & ! intent(out)
                        , psim       & ! function
                        , psih       & ! function
                        , zoobukhov  ! ! function
   use node_mod  , only : mynum      ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real, intent(in)  :: theta_atm    ! Above canopy air pot. temperature        [        K]
   real, intent(in)  :: enthalpy_atm ! Above canopy air specific enthalpy       [     J/kg]
   real, intent(in)  :: shv_atm      ! Above canopy vapour spec. hum.           [kg/kg_air]
   real, intent(in)  :: rvap_atm     ! Above canopy vapour mixing ratio         [kg/kg_air]
   real, intent(in)  :: co2_atm      ! Above canopy CO2 mixing ratio            [ µmol/mol]
   real, intent(in)  :: theta_can    ! Canopy air potential temperature         [        K]
   real, intent(in)  :: enthalpy_can ! Canopy air specific enthalpy             [     J/kg]
   real, intent(in)  :: shv_can      ! Canopy air vapour spec. humidity         [kg/kg_air]
   real, intent(in)  :: rvap_can     ! Canopy air vapour mixing ratio           [kg/kg_air]
   real, intent(in)  :: co2_can      ! Canopy air CO2 mixing ratio              [ µmol/mol]
   real, intent(in)  :: zref         ! Height at reference point                [        m]
   real, intent(in)  :: dheight      ! Displacement height                      [        m]
   real, intent(in)  :: uref         ! Wind speed at reference height           [      m/s]
   real, intent(in)  :: dtl3         ! Time step                                [        m]
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
   real              :: uuse         ! Wind for too stable cases (Rib > Ribmax)
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
   real              :: utotal       ! Total wind (actual + convective)
   real              :: uconv        ! Convective velocity
   real              :: uconv_prev   ! Previous convective velocity
   real              :: change       ! Difference in convective velocity
   integer           :: icnt         ! Iteration counter
   !----- Aux. environment conditions. ----------------------------------------------------!
   real              :: thetav_atm   ! Atmos. virtual potential temperature     [        K]
   real              :: thetav_can   ! Canopy air virtual pot. temperature      [        K]
   !----- External functions. -------------------------------------------------------------!
   real, external    :: cbrt          ! Cubic root
   real, external    :: leaf3_sflux_w ! Surface flux in the w direction
   !---------------------------------------------------------------------------------------!


   !----- Find the variables common to both methods. --------------------------------------!
   thetav_atm = theta_atm * (1. + epim1 * shv_atm)
   thetav_can = theta_can * (1. + epim1 * shv_can)
   zoz0m      = (zref-dheight)/rough
   lnzoz0m    = log(zoz0m)
   zoz0h      = z0moz0h * zoz0m
   lnzoz0h    = log(zoz0h)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the Bulk Richardson number.  For stable cases, and for L79 in both cases,    !
   ! this will be the definitive RiB, whilst this is the first guess, which will be        !
   ! corrected by the convective velocity in the other unstable cases.                     !
   !---------------------------------------------------------------------------------------!
   rib        = 2.0 * grav * (zref-dheight-rough) * (thetav_atm-thetav_can)                &
              / ( (thetav_atm+thetav_can) * uref * uref)
   stable     = thetav_atm >= thetav_can
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Correct the bulk Richardson number in case it's too stable.  We also define a      !
   ! stable case correction to bring down the stars other than ustar, so the flux doesn't  !
   ! increase for stabler cases (it remains constant).                                     !
   !---------------------------------------------------------------------------------------!
   if (rib > ribmax) then
      uuse = sqrt(rib/ribmax) * uref
      rib  = ribmax
   else
      uuse = uref
   end if
   !---------------------------------------------------------------------------------------!




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
      c1   = a2 * uuse

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
      r_aer = 1. / (a2 * uuse * fh)

      !----- Finding ustar, making sure it is not too small. ------------------------------!
      ustar = max(ustmin,sqrt(c1 * uuse * fm))

      !----- Finding the coefficient to scale the other stars. ----------------------------!
      c3 = c1 * fh / ustar
      !------------------------------------------------------------------------------------!


      !----- Compute zeta from u* and T* --------------------------------------------------!
      zeta = grav * vonk * c3 * (theta_atm - theta_can) / (theta_atm * ustar * ustar)
      !------------------------------------------------------------------------------------!


   case default
      !------------------------------------------------------------------------------------!
      ! 2. Here we use the model proposed by OD95, the standard for MM5, but with some     !
      !    terms that were computed in B71 (namely, the "0" terms), which prevent sin-     !
      !    gularities.                                                                     !
      !    However we know zeta, so zeta0 can be written as z0/z * zeta.                   !
      ! 3. Here we use the model proposed by BH91, which is almost the same as the OD95    !
      !    method, except that the stable functions are computed in a more generic way.    !
      !    BH91 claim that the oft-used approximation (-beta*zeta) can cause poor          !
      !    ventilation of the stable layer, leading to decoupling between the atmo-        !
      !    sphere and the canopy air space and excessive cooling                           !
      ! 4. Here we use a similar approach as in CLM04, excepth that the momentum flux      !
      !    gradient function for the unstable case for momentum is switched by a power     !
      !    of -1/6 (kind of the square of the heat one).  This is to guarantee that        !
      !    the psi function doesn't have local maxima/minima.                              !
      !------------------------------------------------------------------------------------!
      !----- Initialise uconv. ------------------------------------------------------------!
      uconv      = 0.0
      !----- Check if we need to go through the iterative process. ------------------------!
      if (stable) then
         !----- We now compute the stability correction functions. ------------------------!
         zeta   = zoobukhov(rib,zref-dheight,rough,zoz0m,lnzoz0m,zoz0h,lnzoz0h,stable)
         zeta0m = rough * zeta / (zref-dheight)
         zeta0h = z0hoz0m * zeta0m

         !----- Find the aerodynamic resistance similarly to L79. -------------------------!
         r_aer = tprandtl * (lnzoz0h - psih(zeta,stable) + psih(zeta0h,stable))            &
                          * (lnzoz0m - psim(zeta,stable) + psim(zeta0m,stable))            &
                          / (vonk * vonk * uuse)

         !----- Find ustar, making sure it is not too small. ------------------------------!
         ustar = max (ustmin, vonk * uuse                                                  &
                            / (lnzoz0m - psim(zeta,stable) + psim(zeta0m,stable)))


         !----- Find the coefficient to scale the other stars. ----------------------------!
         c3    = vonk / (tprandtl * (lnzoz0h - psih(zeta,stable) + psih(zeta0h,stable)))
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !    Unstable case.  Here we run a few iterations to make sure we correct the     !
         ! bulk Richardson number.  This is really a simple correction, so we don't need   !
         ! uconv to be totally in equilibrium.                                             !
         !---------------------------------------------------------------------------------!
         unstable: do icnt=1,6
            !----- Update total winds. ----------------------------------------------------!
            uconv_prev = uconv
            utotal     = sqrt(uuse*uuse + uconv_prev * uconv_prev)
            !------------------------------------------------------------------------------!


            !----- Update the Bulk Richardson number. -------------------------------------!
            rib        = 2.0 * grav * (zref-dheight-rough) * (thetav_atm-thetav_can)       &
                       / ( (thetav_atm+thetav_can) * utotal)
            !------------------------------------------------------------------------------!


            !----- We now compute the stability correction functions. ---------------------!
            zeta   = zoobukhov(rib,zref-dheight,rough,zoz0m,lnzoz0m,zoz0h,lnzoz0h,stable)
            !------------------------------------------------------------------------------!


            !----- Find the coefficient to scale the other stars. -------------------------!
            zeta0m = rough * zeta / (zref-dheight)
            zeta0h = z0hoz0m * zeta0m
            !------------------------------------------------------------------------------!


            !----- Find ustar, making sure it is not too small. ---------------------------!
            ustar = max (ustmin, vonk * uuse                                               &
                               / (lnzoz0m - psim(zeta,stable) + psim(zeta0m,stable)))
            !------------------------------------------------------------------------------!


            !----- Find the coefficient to scale the other stars. -------------------------!
            c3    = vonk                                                                   &
                  / (tprandtl * (lnzoz0h - psih(zeta,stable) + psih(zeta0h,stable)))
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Use potential virtual temperature here because convection is related to  !
            ! buoyancy.                                                                    !
            !------------------------------------------------------------------------------!
            tstar  = c3 * (thetav_atm - thetav_can )
            !------------------------------------------------------------------------------!


            !----- Estimate the convective velocity. --------------------------------------!
            uconv = leaf3_sflux_w(zeta,tstar,ustar) / ustar
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     We are only after a rough estimate of this velocity, so if the differ-   !
            ! ence is less than the RK4 tolerance, then it's enough.                       !
            !------------------------------------------------------------------------------!
            change = 2.0 * abs(uconv-uconv_prev) / (abs(uconv) + abs(uconv_prev))
            if (change < 0.01) exit unstable
            !------------------------------------------------------------------------------!
         end do unstable
      end if
   end select

   !----- Finding all stars. --------------------------------------------------------------!
   tstar = c3 * (theta_atm    - theta_can   )
   estar = c3 * (enthalpy_atm - enthalpy_can)
   qstar = c3 * (shv_atm      - shv_can     )
   rstar = c3 * (rvap_atm     - rvap_can    )
   cstar = c3 * (co2_atm      - co2_can     )

   if (abs(tstar) < 1.e-7) tstar = 0.
   if (abs(estar) < 1.e-7) estar = 0.
   if (abs(qstar) < 1.e-7) qstar = 0.
   if (abs(rstar) < 1.e-7) rstar = 0.
   if (abs(cstar) < 1.e-7) cstar = 0.

   !---------------------------------------------------------------------------------------!
   !    Compute the ground conductance.  This equation is similar to the original, except  !
   ! that we don't assume the ratio between the gradient and the characteristic scale to   !
   ! be 0.2; instead we use the actual ratio that is computed here.                        !
   !---------------------------------------------------------------------------------------!
   ggbare = c3 * ustar
   !---------------------------------------------------------------------------------------!


   return
end subroutine leaf3_stars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine computes the turbulent fluxes of momentum, heat and moisture from the    !
! surface layer using the  Manton-Cotton algebraic surface layer equations.                !
!------------------------------------------------------------------------------------------!
subroutine leaf3_sfclmcv(ustar,tstar,rstar,cstar,zeta,vels_pat,ups,vps,patch_area          &
                        ,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r,sflux_c)
   use rconstants
   use leaf_coms     , only : g_urban ! ! intent(in)
   use teb_spm_start , only : teb_spm ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real , intent(in)    :: ustar
   real , intent(in)    :: tstar
   real , intent(in)    :: rstar
   real , intent(in)    :: cstar
   real , intent(in)    :: zeta
   real , intent(in)    :: vels_pat
   real , intent(in)    :: ups
   real , intent(in)    :: vps
   real , intent(in)    :: patch_area
   real , intent(inout) :: sflux_u
   real , intent(inout) :: sflux_v
   real , intent(inout) :: sflux_w
   real , intent(inout) :: sflux_t
   real , intent(inout) :: sflux_r
   real , intent(inout) :: sflux_c
   !----- Local variables. ----------------------------------------------------------------!
   real                 :: cosine1
   real                 :: sine1
   real                 :: vtscr
   !----- Local constants. ----------------------------------------------------------------!
   real , parameter     :: wtol = 1.e-20
   !----- External functions. -------------------------------------------------------------!
   real , external      :: leaf3_sflux_w
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

   !----- Define vertical flux. -----------------------------------------------------------!
   sflux_w = sflux_w + leaf3_sflux_w(zeta,tstar,ustar) * patch_area

   return
end subroutine leaf3_sfclmcv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    Vertical flux, as in:                                                                 !
!                                                                                          !
!   Manton, M. J., Cotton, W. R., 1977: Parameterization of the atmospheric surface        !
!      layer.  J. Atm. Sci., 34, 331-334.                                                  !
!------------------------------------------------------------------------------------------!
real function leaf3_sflux_w(zeta,tstar,ustar)
   use rconstants , only : vonk ! intent(in)
  
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real, intent(in)    :: zeta
   real, intent(in)    :: ustar
   real, intent(in)    :: tstar
   !----- Local variables -----------------------------------------------------------------!
   real                :: cx
   real                :: psin
   !----- Constants -----------------------------------------------------------------------!
   real, parameter     :: wtol = 1.e-20
   !---------------------------------------------------------------------------------------!

   if (zeta < 0.0)then
      cx = zeta * sqrt(sqrt(1.0 - 15.0 * zeta))
   else
      cx = zeta / (1.0 + 4.7 * zeta)
   endif
  
   psin = sqrt((1.0-2.86 * cx) / (1.0 + cx * (-5.390 + cx * 6.9980 )))
   leaf3_sflux_w = ( 0.27 * max(6.25 * (1.0 - cx) * psin,wtol)                             &
                       - 1.180 * cx * psin) * ustar * ustar
  
   return
end function leaf3_sflux_w
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the ground properties.  By ground we mean the top soil      !
! layer if no temporary surface water/snow exists, or the top temporary surface water/snow !
! layer if it exists.                                                                      !
!                                                                                          !
! References:                                                                              !
!                                                                                          !
! Passerat de Silans, A., 1986: Transferts de masse et de chaleur dans un sol stratifié    !
!     soumis à une excitation amtosphérique naturelle. Comparaison: Modèles-expérience.    !
!     Thesis, Institut National Polytechnique de Grenoble. (P86)                           !
!                                                                                          !
! Noilhan, J., S. Planton, 1989: Simple parameterization of land surface processes for     !
!     meteorological models, Mon. Wea. Rev., 117, 536-549. (NP89)                          !
!                                                                                          !
! Mahfouf, J. F., J. Noilhan, 1991: Comparative study of various formulations of           !
!     evaporation from bare soil using in situ data. J. Appl. Meteorol., 30, 1354-1365.    !
!     (MN91)                                                                               !
!                                                                                          !
! Lee, T. J., R. A. Pielke, 1992: Estimating the soil surface specific humidity. J. Appl.  !
!     Meteorol., 31, 480-484. (LP92)                                                       !
!                                                                                          !
! Lee, T. J., R. A. Pielke, 1993: Corrigendum. J. Appl. Meteorol., 32, 580-580. (LP93)     !
!------------------------------------------------------------------------------------------!
subroutine leaf3_grndvap(topsoil_energy,topsoil_water,topsoil_text,sfcwater_energy_int     &
                        ,sfcwater_nlev,can_rvap,can_prss,ground_rsat,ground_rvap           &
                        ,ground_temp,ground_fliq)

   use leaf_coms  , only : slcpd                  & ! intent(in)
                         , slpots                 & ! intent(in)
                         , slmsts                 & ! intent(in)
                         , soilcp                 & ! intent(in)
                         , slbs                   & ! intent(in)
                         , sfldcap                & ! intent(in)
                         , snowfac                & ! intent(in)
                         , ggsoil                 & ! intent(in)
                         , ggsoil0                & ! intent(in)
                         , kksoil                 & ! intent(in)
                         , can_shv                & ! intent(inout)
                         , leaf3_matric_potential ! ! intent(inout)
   use rconstants , only : gorh2o                 & ! intent(in)
                         , pi1                    & ! intent(in)
                         , wdns                   & ! intent(in)
                         , lnexp_min              & ! intent(in)
                         , huge_num               & ! intent(in)
                         , tiny_num               & ! intent(in)
                         , toodry                 ! ! intent(in)
   use therm_lib  , only : qslif                  & ! function
                         , uextcm2tl              & ! function
                         , uint2tl                ! ! function
   use mem_leaf   , only : igrndvap               ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real, intent(in)  :: topsoil_energy      ! Top soil internal energy          [     J/m³]
   real, intent(in)  :: topsoil_water       ! Top soil water content            [m³_h2o/m³]
   real, intent(in)  :: topsoil_text        ! Top soil texture class            [      ---]
   real, intent(in)  :: sfcwater_energy_int ! Soil internal energy              [     J/kg]
   real, intent(in)  :: sfcwater_nlev       ! # active levels of surface water  [      ---]
   real, intent(in)  :: can_rvap            ! Canopy vapour mixing ratio        [kg_vap/kg]
   real, intent(in)  :: can_prss            ! Canopy pressure                   [       Pa]
   real, intent(out) :: ground_rsat         ! Surface (saturation) mixing ratio [kg_vap/kg]
   real, intent(out) :: ground_rvap         ! Ground equilibrium mixing ratio   [kg_vap/kg]
   real, intent(out) :: ground_temp         ! Surface temperature               [        K]
   real, intent(out) :: ground_fliq         ! Frac. of sfc H2O in liquid phase  [      ---]
   !----- Local variables. ----------------------------------------------------------------!
   integer           :: ksn                 ! # active levels of surface water
   integer           :: nsoil               ! Soil texture class                [      ---]
   real(kind=4)      :: slpotvn             ! soil water potential              [        m]
   real(kind=4)      :: alpha               ! "alpha" term in LP92
   real(kind=4)      :: beta                ! "beta" term in LP92
   real(kind=4)      :: lnalpha             ! ln(alpha)
   real(kind=4)      :: smterm              ! soil moisture term                [     ----]
   real(kind=4)      :: ground_shv          ! ground equilibrium spec hum       [kg_vap/kg]
   real(kind=4)      :: ground_ssh          ! sfc. saturation spec. hum.        [kg_vap/kg]
   real(kind=4)      :: sfcwater_temp       ! Surface temperature               [        K]
   real(kind=4)      :: sfcwater_shv        ! ground equilibrium spec hum       [kg_vap/kg]
   real(kind=4)      :: sfcwater_ssh        ! sfc. saturation spec. hum.        [kg_vap/kg]
   real(kind=4)      :: sfcwater_fliq       ! Frac. of sfc H2O in liquid phase  [      ---]
   !---------------------------------------------------------------------------------------!


   !----- Set the number of temporary surface water (or snow) layers. ---------------------!
   ksn = nint(sfcwater_nlev)
   !---------------------------------------------------------------------------------------!


   !----- Find top soil temperature and canopy air space specific humidity. ---------------!
   nsoil = nint(topsoil_text)
   call uextcm2tl(topsoil_energy,topsoil_water*wdns,slcpd(nsoil),ground_temp,ground_fliq)
   can_shv = can_rvap / (1.0 + can_rvap)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Topsoil_shv is the effective specific humidity of soil.  This value is a         !
   ! combination of the canopy air specific humidity, the saturation specific humidity at  !
   ! the soil temperature.  When the soil tends to dry air soil moisture, topsoil_shv      !
   ! tends to the canopy air space specific humidity, whereas topsoil_shv tends to the     !
   ! saturation value when the soil moisture is near or above field capacity.  These       !
   ! tendencies are determined by the alpha and beta parameters.                           !
   !---------------------------------------------------------------------------------------!
   !----- Compute the saturation specific humidity at top soil temperature. ---------------!
   ground_ssh  = qslif(can_prss,ground_temp)
   !----- Determine alpha. ----------------------------------------------------------------!
   slpotvn      = leaf3_matric_potential(nsoil,topsoil_water)
   lnalpha      = gorh2o * slpotvn / ground_temp
   if (lnalpha > lnexp_min) then
      alpha   = exp(lnalpha)
   else
      alpha   = 0.0
   end if
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !     Determine Beta, following NP89.  However, because we want evaporation to be shut  !
   ! down when the soil approaches the dry air soil moisture, we offset both the soil      !
   ! moisture and field capacity to the soil moisture above dry air soil.  This is         !
   ! necessary to avoid evaporation to be large just slightly above the dry air soil,      !
   ! which would otherwise happen, especially for those soil types rich in clay.           !
   !---------------------------------------------------------------------------------------!
   smterm = min(1.0, max(0.0, (topsoil_water  - soilcp(nsoil))                             &
                            / (sfldcap(nsoil) - soilcp(nsoil)) ))
   beta   = 0.5 * (1.0 - cos (smterm * pi1))
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !     Decide which method to use to find the ground water vapour mixing ratio.          !
   !---------------------------------------------------------------------------------------!
   select case (igrndvap)
   case (0)
      !----- LP92 method. -----------------------------------------------------------------!
      ground_shv = max(can_shv, ground_ssh * alpha * beta + (1.0 - beta) * can_shv)
      ggsoil     = huge_num
      !------------------------------------------------------------------------------------!

   case (1)
      !----- MH91, test 1. ----------------------------------------------------------------!
      ground_shv = max(can_shv, ground_ssh * beta)
      ggsoil     = huge_num
      !------------------------------------------------------------------------------------!

   case (2)
      !----- MH91, test 2. ----------------------------------------------------------------!
      ground_shv = ground_ssh
      ggsoil     = ggsoil0 * exp(kksoil * smterm)
      !------------------------------------------------------------------------------------!

   case (3)
      !----- MH91, test 3. ----------------------------------------------------------------!
      ground_shv = max(can_shv, ground_ssh * beta + (1.0 - beta) * can_shv)
      ggsoil     = huge_num
      !------------------------------------------------------------------------------------!

   case (4)
      !------------------------------------------------------------------------------------!
      !     MH91, test 4.                                                                  !
      !------------------------------------------------------------------------------------!
      ground_shv = max(can_shv, ground_ssh * alpha)
      ggsoil     = ggsoil0 * exp(kksoil * smterm)
      !------------------------------------------------------------------------------------!

   case (5)
      !------------------------------------------------------------------------------------!
      !     Combination of NP89 and P86.                                                   !
      !------------------------------------------------------------------------------------!
      ground_shv = max(can_shv, ground_ssh * beta)
      ggsoil     = ggsoil0 * exp(kksoil * smterm)
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!


   !----- Ground mixing ratio and saturation mixing ratio. --------------------------------!
   ground_rvap = max(toodry, ground_shv / (1.0 - ground_shv))
   ground_rsat = max(toodry, ground_ssh / (1.0 - ground_ssh))
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !   Update conductance in case there is snowpack.                                       !
   !---------------------------------------------------------------------------------------!
   if (ggsoil /= huge_num .and. ksn > 0) then
      call uint2tl(sfcwater_energy_int,sfcwater_temp,sfcwater_fliq)
      !----- Compute the saturation specific humidity at ground temperature. --------------!
      sfcwater_ssh = qslif(can_prss,sfcwater_temp)
      !----- The ground specific humidity in this case is just the saturation value. ------!
      sfcwater_shv = sfcwater_ssh
      !------------------------------------------------------------------------------------!

      ggsoil      = min(huge_num, ggsoil / max(tiny_num,(1.0 - snowfac)))
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf3_grndvap
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
!     This sub-routine finds the precipitation rate.  Unlike original LEAF-3, now we store !
! the rates, not the total, so we can easily scale for each nested time step.              !
!------------------------------------------------------------------------------------------!
subroutine sfc_pcp(m2,m3,mcld,ia,iz,ja,jz,dtlt,theta2,exner2,conprr,bulkpcpg,bulkqpcpg     &
                  ,bulkdpcpg,leafpcpg,leafqpcpg,leafdpcpg)
   use rconstants, only : t3ple        & ! intent(in)
                        , t00          & ! intent(in)
                        , hr_sec       & ! intent(in)
                        , wdnsi        ! ! intent(in)
   use node_mod  , only : mynum        ! ! intent(in)
   use grid_dims , only : str_len      ! ! intent(in)
   use therm_lib , only : tl2uint      & ! function
                        , extheta2temp ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                       , intent(in)    :: m2
   integer                       , intent(in)    :: m3
   integer                       , intent(in)    :: mcld
   integer                       , intent(in)    :: ia
   integer                       , intent(in)    :: iz
   integer                       , intent(in)    :: ja
   integer                       , intent(in)    :: jz
   real                          , intent(in)    :: dtlt
   real   , dimension(m2,m3)     , intent(in)    :: theta2
   real   , dimension(m2,m3)     , intent(in)    :: exner2
   real   , dimension(m2,m3,mcld), intent(in)    :: conprr
   real   , dimension(m2,m3)     , intent(in)    :: bulkpcpg
   real   , dimension(m2,m3)     , intent(in)    :: bulkqpcpg
   real   , dimension(m2,m3)     , intent(in)    :: bulkdpcpg
   real   , dimension(m2,m3)     , intent(inout) :: leafpcpg
   real   , dimension(m2,m3)     , intent(inout) :: leafqpcpg
   real   , dimension(m2,m3)     , intent(inout) :: leafdpcpg
   !----- Local variables. ----------------------------------------------------------------!
   character(len=str_len)                        :: rainfile
   integer                                       :: i
   integer                                       :: j
   integer                                       :: icld
   real                                          :: rain_temp
   real                                          :: fice
   real                                          :: sndeni
   real                                          :: cumpcpg
   real                                          :: cumqpcpg
   real                                          :: cumdpcpg
   !----- Local constants. ----------------------------------------------------------------!
   logical                       , parameter     :: print_rain = .false.
   character(len=10)             , parameter     :: fmth       = '(16(a,1x))' 
   character(len=25)             , parameter     :: fmtb       = '(2(i5,1x),14(es12.5,1x))' 
   !----- Locally saved variables. --------------------------------------------------------!
   logical                       , save          :: first_time = .true. 
   !---------------------------------------------------------------------------------------!


   !----- Re-create the file in case this is the first time the routine is called. --------!
   if (print_rain .and. first_time) then
      first_time = .false.
      write (rainfile,fmt='(a,i3.3,a)') 'rainleaf-',mynum,'.txt'

      open  (unit=63,file=trim(rainfile),status='replace',action='write')
      write (unit=63,fmt=fmth) '    I','    J','        DTLT'                              &
                               ,'       THETA','       EXNER','       RTEMP'               &
                               ,'      CONPRR','    CUM_PCPG','   CUM_QPCPG'               &
                               ,'   CUM_DPCPG','   BULK_PCPG','  BULK_QPCPG'               &
                               ,'  BULK_DPCPG','   LEAF_PCPG','  LEAF_QPCPG'               &
                               ,'  LEAF_DPCPG'
      close (unit=63,status='keep')
   end if
   !---------------------------------------------------------------------------------------!



   !----- Initialise the precipitation variables. -----------------------------------------!
   leafpcpg (:,:) = 0.
   leafqpcpg(:,:) = 0.
   leafdpcpg(:,:) = 0.


   !----- Loop over the horizontal domain. ------------------------------------------------!
   jloop: do j=ja,jz
      iloop: do i=ia,iz

         !----- Estimate the precipitation temperature. -----------------------------------!
         rain_temp     = extheta2temp(exner2(i,j),theta2(i,j))


         !----- Integrate precipitation rate. ---------------------------------------------!
         cumpcpg = 0.
         cloop: do icld =1,mcld
            cumpcpg = cumpcpg + conprr(i,j,icld)
         end do cloop



         !---------------------------------------------------------------------------------!
         !  Precipitation "depth". Snow fraction and density derived from                  !
         !  Jin et al 1999 Hydrol Process. 13:2467-2482 Table 2                            !
         !  [[modified 11/16/09 by MCD]]                                                   !
         !---------------------------------------------------------------------------------!
         if (rain_temp > (t3ple + 2.5)) then
            !----- Rain only. -------------------------------------------------------------!
            fice    = 0.0
            sndeni  = 1. / 189.0

         elseif (rain_temp <= (t3ple + 2.5) .and. rain_temp  > (t3ple + 2.0) ) then
            !------------------------------------------------------------------------------!
            !     60% snow, 40% rain. (N.B. May not be appropriate for sub-tropical        !
            ! regions where the level of the melting layer is higher...).                  !
            !------------------------------------------------------------------------------!
            fice    = 0.6
            sndeni  = 1. / 189.0

         elseif (rain_temp <= (t3ple + 2.0) .and. rain_temp > t3ple ) then
            !------------------------------------------------------------------------------!
            !     Increasing the fraction of snow. (N.B. May not be appropriate for        !
            ! sub-tropical regions where the level of the melting layer is higher...).     !
            !------------------------------------------------------------------------------!
            fice   = min(1.0, 1. + (54.62 - 0.2*rain_temp))
            sndeni = 1. / (50.0+1.7*(rain_temp-258.15)**1.5 )

         elseif (rain_temp <= t3ple .and. rain_temp > (t3ple - 15.0)) then
            !----- Below freezing point, snow only. ---------------------------------------!
            fice   = 1.0
            sndeni = 1. / (50.0+1.7*(rain_temp-258.15)**1.5 )

         else ! if (rain_temp < (t3ple - 15.0)) then
            !----- Below freezing point, snow only. ---------------------------------------!
            fice   = 1.0
            sndeni = 1. / 50.
         end if
         cumdpcpg  = cumpcpg * ((1.0-fice) * wdnsi + fice * sndeni)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Set internal energy.  This will be the precipitation times the specific     !
         ! internal energy of water (above or at triple point) multiplied by the liquid    !
         ! fraction plus the specific internal energy of ice (below or at the triple       !
         ! point) multiplied by the ice fraction.                                          !
         !---------------------------------------------------------------------------------!
         cumqpcpg = cumpcpg  * ( (1.0-fice) * tl2uint(max(t3ple,rain_temp),1.0)            &
                               +      fice  * tl2uint(min(rain_temp,t3ple),0.0) )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Leaf total precipitation, and associated internal energy and depth, is the  !
         ! total amount integrated over one leaf time step.                                !
         !---------------------------------------------------------------------------------!
         leafpcpg(i,j)  = cumpcpg   + bulkpcpg (i,j) / dtlt
         leafqpcpg(i,j) = cumqpcpg  + bulkqpcpg(i,j) / dtlt
         leafdpcpg(i,j) = cumdpcpg  + bulkdpcpg(i,j) / dtlt
         !---------------------------------------------------------------------------------!

         if (leafpcpg(i,j) > 0.0 .and. print_rain) then
            write (rainfile,fmt='(a,i3.3,a)') 'rainleaf-',mynum,'.txt'
            open  (unit=63,file=trim(rainfile),status='old',action='write'                 &
                  ,position='append')
            write (unit=63,fmt=fmtb)  i,j,dtlt,theta2(i,j),exner2(i,j)                     &
                                    , rain_temp-t00,conprr(i,j,1)*hr_sec                   &
                                    , cumpcpg      ,cumqpcpg      ,cumdpcpg*1000.          &
                                    , bulkpcpg(i,j),bulkqpcpg(i,j),bulkdpcpg(i,j)*1000.    &
                                    , leafpcpg(i,j),leafqpcpg(i,j),leafdpcpg(i,j)*1000.
            close (unit=63,status='keep')
         end if
      end do iloop
   end do jloop

   return
end subroutine sfc_pcp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine updates several vegetation properties such as LAI, heat capacity,    !
! roughness, displacement height, and internal energy.                                     !
!------------------------------------------------------------------------------------------!
subroutine veg_misc_update(ifm,patch_area,leaf_class,veg_fracarea,veg_lai,veg_tai          &
                          ,veg_rough,veg_height,veg_displace,veg_albedo,veg_ndvip          &
                          ,veg_ndvic,veg_ndvif,veg_agb,veg_energy,veg_water,veg_hcap       &
                          ,psibar_10d)

   use leaf_coms , only : nvtyp          & ! intent(in)
                        , nvtyp_teb      & ! intent(in)
                        , veg_temp       & ! intent(inout)
                        , veg_fliq       & ! intent(inout)
                        , timefac_ndvi   & ! intent(in)
                        , tai_min        & ! intent(in)
                        , sr_max         & ! intent(in)
                        , tai_max        & ! intent(in)
                        , veg_clump      & ! intent(in)
                        , glai_max       & ! intent(in)
                        , dead_frac      & ! intent(in)
                        , phenology      & ! intent(in)
                        , sai            & ! intent(in)
                        , vh2dh          & ! intent(in)
                        , albv_green     & ! intent(in)
                        , albv_brown     & ! intent(in)
                        , veg_frac       & ! intent(in)
                        , hcapveg_ref    & ! intent(in)
                        , hcapveg_hmin   & ! intent(in)
                        , sla_0          & ! intent(in)
                        , sla_m          & ! intent(in)
                        , agb_am14_a     & ! intent(in)
                        , agb_am14_b     & ! intent(in)
                        , gu_spheat_leaf & ! intent(in)
                        , gu_spheat_wood ! ! intent(in)
   use therm_lib , only : uextcm2tl      & ! function
                        , cmtl2uext      ! ! function
   use io_params , only : ndvitime1      & ! intent(in)
                        , ndvitime2      & ! intent(in)
                        , iupdndvi       & ! intent(in)
                        , ndviflg        & ! intent(in)
                        , iuselai        ! ! intent(in)
   use mem_leaf  , only : isfcl          ! ! intent(in)
   use mem_grid  , only : time           ! ! intent(in)
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
   real                            , intent(out)   :: veg_displace
   real                            , intent(out)   :: veg_rough
   real                            , intent(out)   :: veg_albedo
   real                            , intent(inout) :: veg_ndvic
   real                            , intent(out)   :: veg_agb
   real                            , intent(inout) :: veg_energy
   real                            , intent(in)    :: veg_water
   real                            , intent(inout) :: veg_hcap
   real                            , intent(in)    :: psibar_10d
   !----- Local variables. ----------------------------------------------------------------!
   integer                                         :: nveg
   real                                            :: sr
   real                                            :: fpar
   real                                            :: dead_lai
   real                                            :: green_frac
   real                                            :: bleaf
   real                                            :: btwig
   real                                            :: spheat_leaf
   real                                            :: spheat_wood
   real                                            :: vhgt_dum
   !----- Local constants. ----------------------------------------------------------------!
   real                            , parameter     :: sr_min     =  1.081
   real                            , parameter     :: fpar_min   =  0.001
   real                            , parameter     :: fpar_max   =  0.950
   real                            , parameter     :: fpcon      = -0.3338082
   real                            , parameter     :: ccc        = -2.9657
   real                            , parameter     :: bz         =  0.91
   real                            , parameter     :: hz         =  0.0075
   real                            , parameter     :: extinc_veg = 0.75
   real                            , parameter     :: fbranch    = 0.20
   !----- Locally saved variables. --------------------------------------------------------!
   logical                         , save          :: nvcall     = .true.
   real, dimension(nvtyp+nvtyp_teb), save          :: dfpardsr
   logical                         , save          :: first_time = .true.
   !---------------------------------------------------------------------------------------!


   !----- Initialise the dfpardsr array, which will be used to compute LAI. ---------------!
   if (nvcall) then
      nvcall = .false.
      do nveg = 1,(nvtyp+nvtyp_teb)
         dfpardsr(nveg) = (fpar_max - fpar_min) / (sr_max(nveg) - sr_min)
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- Find the time interpolation factor for updating NDVI or LAI. --------------------!
   if (iupdndvi == 0) then
      timefac_ndvi = 0.
   else if ( ndvitime1(ifm) == ndvitime2(ifm) ) then
      write(unit=*,fmt="(a)"          ) "-------------------------------------------------"
      write(unit=*,fmt="(a,1x,i14)"   ) " IFM      = ",ifm
      write(unit=*,fmt="(a,1x,es14.7)") " NDVITIME1 = ",ndvitime1(ifm)
      write(unit=*,fmt="(a,1x,es14.7)") " NDVITIME2 = ",ndvitime2(ifm)
      write(unit=*,fmt="(a)")           "-------------------------------------------------"
      call abort_run("NDVI/LAI times must be different!","vegndvi","leaf3_utils.f90")
   else
      timefac_ndvi = sngl((time - ndvitime1(ifm)) / (ndvitime2(ifm) - ndvitime1(ifm)))
   end if
   !---------------------------------------------------------------------------------------!



   !----- Alias for vegetation class. -----------------------------------------------------!
   nveg = nint(leaf_class)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    We only compute LAI and related variables for those vegetation types that can hold !
   ! some actual vegetation.                                                               !
   !---------------------------------------------------------------------------------------!
   if (tai_max(nveg) < tai_min) then
      veg_lai      = 0.
      veg_tai      = 0.
      veg_rough    = 0.
      veg_albedo   = 0.
      veg_fracarea = 0.
   else
      
      !------------------------------------------------------------------------------------!
      !  Time-interpolate NDVI to get current value veg_ndvic(i,j) for this patch.  Limit  !
      ! NDVI to prevent values > .99 to prevent division by zero.                          !
      !------------------------------------------------------------------------------------!

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
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Here we decide which way we are going to compute LAI.  If NDVIFLG is 0 or 2,   !
      ! we use LEAF-3 phenology, otherwise we use NDVI to prescribe LAI.                   !
      !------------------------------------------------------------------------------------!
      select case (ndviflg(ifm))
      case (1)
         !------ We've read information from files, use it. -------------------------------!

         if (iuselai == 1) then
            !----- Input data were LAI, copy it. ------------------------------------------!
            veg_lai = veg_ndvic

         else

            !----- Compute "simple ratio" and limit between sr_min and sr_max(nveg). ------!
            sr = min(sr_max(nveg), max(sr_min, (1. + veg_ndvic) / (1. - veg_ndvic) ) )


            !----- Compute fpar. ----------------------------------------------------------!
            fpar = fpar_min + (sr - sr_min) * dfpardsr(nveg)

            !------------------------------------------------------------------------------!
            !      Compute green leaf area index (veg_lai), dead leaf area index           !
            ! (dead_lai), total area index (tai), and green fraction.                      !
            !------------------------------------------------------------------------------!
            veg_lai    = glai_max(nveg) * (       veg_clump(nveg)  * fpar / fpar_max       &
                                          + (1. - veg_clump(nveg)) * alog(1. - fpar)       &
                                          * fpcon )
            !------------------------------------------------------------------------------!
         end if

      case default
         !---------------------------------------------------------------------------------!
         !    We haven't read anything, use the table value, scaled by phenology if needed !
         ! by this PFT.                                                                    !
         !---------------------------------------------------------------------------------!
         select case (phenology(nveg))
         case (4)
            !----- Use drought/cold phenology. --------------------------------------------!
            veg_lai = glai_max(nveg) * max(0.02,min(1.0,psibar_10d))
         case default
            !----- Evergreen. -------------------------------------------------------------!
            veg_lai = glai_max(nveg)
            !------------------------------------------------------------------------------!
         end select
      end select
      dead_lai   = (glai_max(nveg) - veg_lai) * dead_frac(nveg)
      veg_tai    = veg_lai + sai(nveg) + dead_lai
      green_frac = veg_lai / veg_tai
      !------------------------------------------------------------------------------------!


      !----- Compute vegetation roughness height, albedo, and fractional area. ------------!
      veg_rough    = veg_height * (1. - bz * exp(-hz * veg_tai))
      veg_displace = veg_height * vh2dh
      veg_albedo   = albv_green(nveg) * green_frac + albv_brown(nveg) * (1. - green_frac)
      veg_fracarea = veg_frac(nveg) * (1. - exp(-extinc_veg * veg_tai))
   end if
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !    Now we update the AGB, energy, and heat capacity.                                  !
   !---------------------------------------------------------------------------------------!
   veg_agb      = agb_am14_a * veg_height ** agb_am14_b
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Save temperature and liquid water fraction, before we update heat capacity.        !
   !---------------------------------------------------------------------------------------!
   call uextcm2tl(veg_energy,veg_water,veg_hcap,veg_temp,veg_fliq)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Calculate the heat capacity.                                                      !
   !---------------------------------------------------------------------------------------!
   if (isfcl == 4 .and. veg_tai >= tai_min .and. sla_0(nveg) > 0.0) then
      !------------------------------------------------------------------------------------!
      !     Heat capacity is based on ED for leaves, using CLM-4 SLA profile.  For the     !
      ! woody component, we use the TCH->AGB equation from Asner and Mascaro (2014).       !
      !------------------------------------------------------------------------------------!
      if (sla_m(nveg) == 0.) then
         bleaf = veg_lai / sla_0(nveg)
      else
         bleaf = 2.0 * log(1.0 + sla_m(nveg) * veg_lai / sla_0(nveg) ) / sla_m(nveg)
      end if
      btwig    = 2.0 * fbranch * veg_agb
      veg_hcap = bleaf * gu_spheat_leaf + btwig * gu_spheat_wood
      !------------------------------------------------------------------------------------!

      !if (first_time) then
      !   first_time = .false.
      !   write(unit=61,fmt='(a5,7(1x,a12))')                                               &
      !          ' NVEG','     VEG_LAI','     VEG_TAI','     VEG_AGB','       BLEAF'        &
      !                 ,'       BTWIG','    VEG_HCAP','    VHGT_DUM'
      !end if
      !vhgt_dum = hcapveg_ref * max(veg_height,hcapveg_hmin)
      !write(unit=61,fmt='(i5,7(1x,es12.5))')                                               &
      !   nveg,veg_lai,veg_tai,veg_agb,bleaf,btwig,veg_hcap,vhgt_dum
   else
      !----- Use the old style. -----------------------------------------------------------!
      veg_hcap = hcapveg_ref * max(veg_height,hcapveg_hmin)
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Since vegetation heat capacity may have changed, we must update energy.           !
   !---------------------------------------------------------------------------------------!
   if (veg_hcap > 0.0) then
      veg_energy = cmtl2uext(veg_hcap,veg_water,veg_temp,veg_fliq)
   else
      veg_energy = 0.0
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine veg_misc_update
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function determines the wind at a given height, given that the stars are al-     !
! ready known, as well as the Richardson number and the zetas to find the wind at the top  !
! of the canopy.  The result is the TAI weighted average wind profile, following the       !
! exponential decay proposed by Leuning (1995).                                            !
!------------------------------------------------------------------------------------------!
real(kind=4) function leaf3_reduced_wind(ustar,zeta,rib,zref,dheight,height,rough,vfarea   &
                                        ,stai)
   use rconstants     , only : vonk      & ! intent(in)
                             , lnexp_min ! ! intent(in)
   use leaf_coms      , only : bl79      & ! intent(in)
                             , csm       & ! intent(in)
                             , csh       & ! intent(in)
                             , dl79      & ! intent(in)
                             , ugbmin    & ! intent(in)
                             , tai_min   & ! intent(in)
                             , psim      ! ! function
   use mem_leaf       , only : istar     ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4), intent(in) :: ustar        ! Friction velocity                   [    m/s]
   real(kind=4), intent(in) :: zeta         ! Normalised height                   [    ---]
   real(kind=4), intent(in) :: rib          ! Bulk Richardson number              [    ---]
   real(kind=4), intent(in) :: zref         ! Reference height                    [      m]
   real(kind=4), intent(in) :: dheight      ! Displacement height                 [      m]
   real(kind=4), intent(in) :: height       ! Height to determine the red. wind   [      m]
   real(kind=4), intent(in) :: rough        ! Roughness scale                     [      m]
   real(kind=4), intent(in) :: vfarea       ! Vegetation fraction                 [  m2/m2]
   real(kind=4), intent(in) :: stai         ! Exposed TAI                         [  m2/m2]
   !----- Local variables. ----------------------------------------------------------------!
   logical                  :: stable       ! Canopy air space is stable          [    T|F]
   real(kind=4)             :: zetah        ! Zeta for h=height                   [    ---]
   real(kind=4)             :: zeta0        ! Zeta for h=rough                    [    ---]
   real(kind=4)             :: hoz0         ! ((h-d0)/z0)                         [    ---]
   real(kind=4)             :: lnhoz0       ! ln ((h-d0)/z0)                      [    ---]
   real(kind=4)             :: a2           ! Drag coeff. in neutral conditions
   real(kind=4)             :: fm           ! Stability parameter for momentum
   real(kind=4)             :: c2           ! Part of the c coefficient.
   real(kind=4)             :: cm           ! c coefficient times |Rib|^1/2
   real(kind=4)             :: ee           ! (z/z0)^1/3 -1. for eqn. 20 (L79)
   real(kind=4)             :: veg_wind_top ! Wind at the top of the canopy       [    m/s]
   real(kind=4)             :: lnexp_now    ! Exponential term                    [    ---]
   !----- External functions. -------------------------------------------------------------!
   real(kind=4), external   :: cbrt      ! Cubic root
   !---------------------------------------------------------------------------------------!



   !----- Define whether the layer is stable or not. --------------------------------------!
   stable    = rib >= 0.
   !---------------------------------------------------------------------------------------!



   !----- Find the log for the log-height interpolation of wind. --------------------------!
   hoz0      = (height-dheight)/rough
   lnhoz0    = log(hoz0)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The wind at a given height is found by using the same definition of wind speed at !
   ! a given height.                                                                       !
   !---------------------------------------------------------------------------------------!
   select case (istar)
   case (1) !---- Louis (1979) method. ----------------------------------------------------!

      !----- Compute the a-square factor and the coefficient to find theta*. --------------!
      a2   = vonk * vonk / (lnhoz0 * lnhoz0)

      if (stable) then
         !----- Stable case ---------------------------------------------------------------!
         fm = 1.0 / (1.0 + (2.0 * bl79 * rib / sqrt(1.0 + dl79 * rib)))

      else
         !---------------------------------------------------------------------------------!
         !     Unstable case.  The only difference from the original method is that we no  !
         ! longer assume z >> z0, so the "c" coefficient uses the full z/z0 term.          !
         !---------------------------------------------------------------------------------!
         ee = cbrt(hoz0) - 1.
         c2 = bl79 * a2 * ee * sqrt(ee * abs(rib))
         cm = csm * c2
         fm = (1.0 - 2.0 * bl79 * rib / (1.0 + 2.0 * cm))
      end if
      
      !----- Find the wind. ---------------------------------------------------------------!
      veg_wind_top = (ustar/vonk) * (lnhoz0/sqrt(fm))

   case default  !----- Other methods. ----------------------------------------------------!

      !----- Determine zeta for the sought height and for roughness height. ---------------!
      zetah = zeta * (height-dheight) / (zref-dheight)
      zeta0 = zeta * rough            / (zref-dheight)
      !------------------------------------------------------------------------------------!

      veg_wind_top = (ustar/vonk) * (lnhoz0 - psim(zetah,stable) + psim(zeta0,stable))

   end select
   !---------------------------------------------------------------------------------------!



   !----- Impose the wind to be more than the minimum. ------------------------------------!
   if (stai > tai_min .and. vfarea > 0.1) then
      lnexp_now = max(lnexp_min,-stai/vfarea)
      leaf3_reduced_wind = max( veg_wind_top                                               &
                              * ( 1. - vfarea + vfarea * vfarea * exp(lnexp_now) / stai )  &
                              , ugbmin )
   else
      leaf3_reduced_wind = max(veg_wind_top, ugbmin)
   end if
   !---------------------------------------------------------------------------------------!


   return
end function leaf3_reduced_wind
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine computes the aerodynamic conductance between leaf and canopy air    !
! space for both heat and water vapour, based on:                                          !
!                                                                                          !
! L95 - Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf           !
!       nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to    !
!       canopies.  Plant, Cell and Environ., 18, 1183-1200.                                !
! M08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,        !
!       3rd. edition, Academic Press, Amsterdam, 418pp.  (Mostly Chapter 10).              !
!                                                                                          !
! Notice that the units are somewhat different from L95.                                   !
! - gbh is in J/(K m2 s), and                                                              !
! - gbw is in kg_H2O/m2/s.                                                                 !
!------------------------------------------------------------------------------------------!
subroutine leaf3_aerodynamic_conductances(leaf_class,veg_wind,veg_temp,can_temp,can_shv    &
                                         ,can_rhos,can_cp)
   use leaf_coms , only : leaf_width & ! intent(in)
                        , aflat_turb & ! intent(in)
                        , aflat_lami & ! intent(in)
                        , nflat_turb & ! intent(in)
                        , nflat_lami & ! intent(in)
                        , bflat_turb & ! intent(in)
                        , bflat_lami & ! intent(in)
                        , mflat_turb & ! intent(in)
                        , mflat_lami & ! intent(in)
                        , gbh_2_gbw  & ! intent(in)
                        , gbh        & ! intent(in)
                        , gbw        ! ! intent(in)
   use rconstants, only : t00        & ! intent(in)
                        , grav       & ! intent(in)
                        , kin_visc0  & ! intent(in)
                        , dkin_visc  & ! intent(in)
                        , th_diff0   & ! intent(in)
                        , dth_diff   ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4)   , intent(in)  :: leaf_class      ! Vegetation class           [      ---]
   real(kind=4)   , intent(in)  :: veg_wind        ! Wind at cohort height      [      m/s]
   real(kind=4)   , intent(in)  :: veg_temp        ! Leaf temperature           [        K]
   real(kind=4)   , intent(in)  :: can_temp        ! Canopy air temperature     [        K]
   real(kind=4)   , intent(in)  :: can_shv         ! Canopy air spec. hum.      [    kg/kg]
   real(kind=4)   , intent(in)  :: can_rhos        ! Canopy air density         [    kg/m³]
   real(kind=4)   , intent(in)  :: can_cp          ! Canopy air spec. heat      [   J/kg/K]
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: iveg            ! Vegetation class           [      ---]
   real(kind=4)                 :: lwidth          ! Leaf width                 [        m]
   real(kind=4)                 :: kin_visc        ! Kinematic viscosity        [     m²/s]
   real(kind=4)                 :: th_diff         ! Kinematic viscosity        [     m²/s]
   real(kind=4)                 :: th_expan        ! Thermal expansion          [      1/K]
   real(kind=4)                 :: gr_coeff        ! grav*th_expan/kin_visc²    [   1/K/m³]
   real(kind=4)                 :: grashof         ! Grashof number             [      ---]
   real(kind=4)                 :: reynolds        ! Reynolds number            [      ---]
   real(kind=4)                 :: nusselt_lami    ! Nusselt number (laminar)   [      ---]
   real(kind=4)                 :: nusselt_turb    ! Nusselt number (turb.)     [      ---]
   real(kind=4)                 :: nusselt         ! Nusselt number             [      ---]
   real(kind=4)                 :: forced_gbh_mos  ! Forced convection cond.    [      m/s]
   real(kind=4)                 :: free_gbh_mos    ! Free convection cond.      [      m/s]
   real(kind=4)                 :: gbh_mos         ! Total convection cond.     [      m/s]
   !---------------------------------------------------------------------------------------!


   !----- Vegetation class. ---------------------------------------------------------------!
   iveg   = nint(leaf_class)
   !---------------------------------------------------------------------------------------!


   !----- Save the leaf width of this PFT. ------------------------------------------------!
   lwidth = leaf_width(iveg)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Compute kinematic viscosity, thermal diffusivity, and expansion coefficient as    !
   ! functions of temperature.  Here we use the canopy air space temperature because       !
   ! this is the representative temperature of the fluid.                                  !
   !                                                                                       !
   !     Kinematic viscosity and thermal diffusivity are determined from MU08, see         !
   ! discussion on page 32.  Thermal expansion is assumed to be of an ideal gas (1/T),     !
   ! like in Dufour and van Mieghem (1975), for example.                                   !
   !---------------------------------------------------------------------------------------!
   th_expan = 1.0 / can_temp
   !----- kin_visc and th_diff are assumed linear functions of temperature. ---------------!
   kin_visc = kin_visc0 * ( 1.0 + dkin_visc * ( can_temp - t00 ) )
   th_diff  = th_diff0  * ( 1.0 + dth_diff  * ( can_temp - t00 ) )
   !---------------------------------------------------------------------------------------!
   !    Grashof coefficient (a*g/nu²) in MU08's equation 10.8.                             !
   !---------------------------------------------------------------------------------------!
   gr_coeff = th_expan * grav  / ( kin_visc * kin_visc )
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the conductance, in m/s, associated with forced convection.                  !
   !---------------------------------------------------------------------------------------!
   !----- 1. Compute the Reynolds number. -------------------------------------------------!
   reynolds        = veg_wind * lwidth / th_diff
   !----- 2. Compute the Nusselt number for both the laminar and turbulent case. ----------!
   nusselt_lami    = aflat_lami * reynolds ** nflat_lami
   nusselt_turb    = aflat_turb * reynolds ** nflat_turb
   !----- 4. The right Nusselt number is the largest. -------------------------------------!
   nusselt         = max(nusselt_lami,nusselt_turb)
   !----- 5. The conductance is given by MU08 - equation 10.4 -----------------------------!
   forced_gbh_mos  = th_diff * nusselt / lwidth
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the conductance, in m/s,  associated with free convection.                   !
   !---------------------------------------------------------------------------------------!
   !----- 1. Find the Grashof number. -----------------------------------------------------!
   grashof         = gr_coeff  * abs(veg_temp - can_temp) * lwidth * lwidth * lwidth
   !----- 2. Compute the Nusselt number for both the laminar and turbulent case. ----------!
   nusselt_lami    = bflat_lami * grashof ** mflat_lami
   nusselt_turb    = bflat_turb * grashof ** mflat_turb
   !----- 4. The right Nusselt number is the largest. -------------------------------------!
   nusselt         = max(nusselt_lami,nusselt_turb)
   !----- 5. The conductance is given by MU08 - equation 10.4 -----------------------------!
   free_gbh_mos    = th_diff * nusselt / lwidth
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     The heat conductance for the thermodynamic budget is the sum of conductances,     !
   ! because we assume both forms of convection happen parallelly.  The conversion from    !
   ! heat to water conductance (in m/s) can be found in L95, page 1198, after equation E5. !
   ! For the ED purposes, the output variables are converted to the units of entropy and   !
   ! water fluxes [J/K/m²/s and kg/m²/s, respectively].                                    !
   !---------------------------------------------------------------------------------------!
   gbh_mos = free_gbh_mos + forced_gbh_mos
   gbh     =             gbh_mos * can_rhos * can_cp
   gbw     = gbh_2_gbw * gbh_mos * can_rhos
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf3_aerodynamic_conductances
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine copies some atmospheric fields from the 2-D arrays to the common   !
! module variable.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine leaf3_atmo1d(m2,m3,i,j,thp,theta,rv,rtp,co2p,up,vp,pitot,dens,height            &
                       ,pcpg,qpcpg,dpcpg)
   use leaf_coms , only : ubmin        & ! intent(in)
                        , atm_up       & ! intent(out)
                        , atm_vp       & ! intent(out)
                        , atm_thil     & ! intent(out)
                        , atm_theta    & ! intent(out)
                        , atm_rvap     & ! intent(out)
                        , atm_rtot     & ! intent(out)
                        , atm_shv      & ! intent(out)
                        , geoht        & ! intent(out)
                        , atm_exner    & ! intent(out)
                        , atm_co2      & ! intent(out)
                        , atm_prss     & ! intent(out)
                        , atm_rhos     & ! intent(out)
                        , atm_vels     & ! intent(out)
                        , atm_temp     & ! intent(out)
                        , atm_theiv    & ! intent(out)
                        , atm_vpdef    & ! intent(out)
                        , pcpgl        & ! intent(out)
                        , qpcpgl       & ! intent(out)
                        , dpcpgl       ! ! intent(out)
   use rconstants, only : srtwo        ! ! intent(in)
   use therm_lib , only : thetaeiv     & ! function
                        , vpdefil      & ! function
                        , exner2press  & ! function
                        , extheta2temp ! ! function

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in) :: m2
   integer                  , intent(in) :: m3
   integer                  , intent(in) :: i
   integer                  , intent(in) :: j
   real   , dimension(m2,m3), intent(in) :: thp
   real   , dimension(m2,m3), intent(in) :: theta
   real   , dimension(m2,m3), intent(in) :: rv
   real   , dimension(m2,m3), intent(in) :: rtp
   real   , dimension(m2,m3), intent(in) :: co2p
   real   , dimension(m2,m3), intent(in) :: up
   real   , dimension(m2,m3), intent(in) :: vp
   real   , dimension(m2,m3), intent(in) :: pitot
   real   , dimension(m2,m3), intent(in) :: dens
   real   , dimension(m2,m3), intent(in) :: height
   real   , dimension(m2,m3), intent(in) :: pcpg
   real   , dimension(m2,m3), intent(in) :: qpcpg
   real   , dimension(m2,m3), intent(in) :: dpcpg
   !----- Local variables. ----------------------------------------------------------------!
   real                                  :: wfact
   real                                  :: vels_1st
   real                                  :: tmp
   !---------------------------------------------------------------------------------------!


   !----- Find the wind speed and make sure it is above a minimum. ------------------------!
   vels_1st     = sqrt(up(i,j)*up(i,j) + vp(i,j)*vp(i,j))
   if (vels_1st == 0.0) then
      atm_vels  = ubmin
      atm_up    = 0.5 * srtwo * ubmin
      atm_vp    = atm_up
   else
      wfact     = max(1.0, ubmin / vels_1st)
      atm_vels  = vels_1st * wfact
      atm_up    = up(i,j)  * wfact
      atm_vp    = vp(i,j)  * wfact
   end if
   !---------------------------------------------------------------------------------------!


   !----- Copy the other values stored at the 2-D arrays to the leaf common. --------------!
   atm_thil     = thp(i,j)
   atm_theta    = theta(i,j)
   atm_rvap     = rv(i,j)
   atm_rtot     = rtp(i,j)
   atm_shv      = atm_rvap / (1. + atm_rvap)
   geoht        = height(i,j)
   atm_exner    = pitot(i,j)
   atm_co2      = co2p(i,j)
   atm_prss     = exner2press(atm_exner)
   atm_temp     = extheta2temp(atm_exner,atm_theta)
   atm_theiv    = thetaeiv(atm_thil,atm_prss,atm_temp,atm_rvap,atm_rtot)
   atm_vpdef    = vpdefil(atm_prss,atm_temp,atm_shv,.true.)
   pcpgl        = pcpg(i,j)
   qpcpgl       = qpcpg(i,j)
   dpcpgl       = dpcpg(i,j)
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf3_atmo1d
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine assigns various canopy air space variables for the case in which   !
! leaf is not solved.                                                                      !
!------------------------------------------------------------------------------------------!
subroutine leaf0(m2,m3,mpat,i,j,can_theta,can_rvap,can_co2,can_prss,can_theiv,can_vpdef    &
                ,patch_area)
   use mem_leaf  , only : dthcon         & ! intent(in)
                        , drtcon         & ! intent(in)
                        , pctlcon        ! ! intent(in)
   use leaf_coms , only : atm_theta      & ! intent(in)
                        , atm_rvap       & ! intent(in)
                        , atm_co2        & ! intent(in)
                        , atm_prss       & ! intent(in)
                        , atm_shv        & ! intent(in)
                        , geoht          & ! intent(in)
                        , min_patch_area & ! intent(in)
                        , can_depth_min  & ! intent(in)
                        , can_shv        & ! intent(out)
                        , can_rsat       & ! intent(out)
                        , can_rhv        & ! intent(out)
                        , can_exner      & ! intent(out)
                        , can_temp       & ! intent(out)
                        , can_enthalpy   & ! intent(out)
                        , can_rhos       & ! intent(out)
                        , can_cp         ! ! intent(out)
   use rconstants, only : ep             & ! intent(in)
                        , cpdry          & ! intent(in)
                        , cph2o          ! ! intent(in)
   use therm_lib , only : thetaeiv       & ! function
                        , vpdefil        & ! function
                        , rslif          & ! function
                        , reducedpress   & ! function
                        , idealdenssh    & ! function
                        , press2exner    & ! function
                        , extheta2temp   & ! function
                        , tq2enthalpy    ! ! function

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                       , intent(in)    :: m2
   integer                       , intent(in)    :: m3
   integer                       , intent(in)    :: mpat
   integer                       , intent(in)    :: i
   integer                       , intent(in)    :: j
   real   , dimension(m2,m3,mpat), intent(inout) :: can_theta
   real   , dimension(m2,m3,mpat), intent(inout) :: can_rvap
   real   , dimension(m2,m3,mpat), intent(inout) :: can_co2
   real   , dimension(m2,m3,mpat), intent(inout) :: can_prss
   real   , dimension(m2,m3,mpat), intent(inout) :: can_theiv
   real   , dimension(m2,m3,mpat), intent(inout) :: can_vpdef
   real   , dimension(m2,m3,mpat), intent(inout) :: patch_area
   !---------------------------------------------------------------------------------------!



   !----- Fill the canopy properties. -----------------------------------------------------!
   can_theta(i,j,2)  = atm_theta - dthcon
   can_rvap(i,j,2)   = atm_rvap  - drtcon
   can_co2(i,j,2)    = atm_co2

   can_shv           = can_rvap(i,j,2) / (1. + can_rvap(i,j,2))
   
   can_prss(i,j,2)   = reducedpress(atm_prss,atm_theta,atm_shv,geoht,can_theta(i,j,2)      &
                                   ,can_shv,can_depth_min)

   can_exner         = press2exner(can_prss(i,j,2))
   can_temp          = extheta2temp(can_exner,can_theta(i,j,2))

   can_rsat          = rslif(can_prss(i,j,2),can_temp)

   can_rhv           = can_rvap(i,j,2) * (ep + can_rsat)                                   &
                     / ( can_rsat * (ep + can_rvap(i,j,2)))

   can_theiv(i,j,2)  = thetaeiv(can_theta(i,j,2),can_prss(i,j,2),can_temp,can_rvap(i,j,2)  &
                               ,can_rvap(i,j,2))

   can_vpdef(i,j,2)  = vpdefil(can_prss(i,j,2),can_temp,can_shv,.true.)

   can_enthalpy      = tq2enthalpy(can_temp,can_shv,.true.)
   can_rhos          = idealdenssh(can_prss(i,j,2),can_temp,can_shv)
   can_cp            = (1.0 - can_shv) * cpdry + can_shv * cph2o
   !---------------------------------------------------------------------------------------!



   !----- Impose area to be bounded. ------------------------------------------------------!
   patch_area(i,j,1) = min(1.0,max(min_patch_area,1.0-pctlcon))
   patch_area(i,j,2) = 1.0 - patch_area(i,j,1)
   !---------------------------------------------------------------------------------------!


   return
end subroutine leaf0
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine computes the roughness for a given patch.                           !
!------------------------------------------------------------------------------------------!
subroutine leaf3_roughness(ip,veg_fracarea,patch_area,ustar,topzo,veg_rough,soil_rough     &
                          ,patch_rough)
   use mem_leaf , only : isfcl          ! ! intent(in)
   use leaf_coms, only : min_patch_area & ! intent(in)
                       , z0fac_water    & ! intent(in)
                       , min_waterrough & ! intent(in)
                       , snowfac        & ! intent(in)
                       , snowrough      ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)  :: ip
   real   , intent(in)  :: veg_fracarea
   real   , intent(in)  :: patch_area
   real   , intent(in)  :: ustar
   real   , intent(in)  :: topzo
   real   , intent(in)  :: veg_rough
   real   , intent(in)  :: soil_rough
   real   , intent(out) :: patch_rough
   !----- Local variables. ----------------------------------------------------------------!
   real                 :: summer_rough
   !---------------------------------------------------------------------------------------!


   select case (ip)
   case (1)
      !------------------------------------------------------------------------------------!
      !    For water surfaces (patch 1), compute roughness length based on previous ustar. !
      !------------------------------------------------------------------------------------!
      patch_rough = max(z0fac_water * ustar ** 2, min_waterrough)
      !------------------------------------------------------------------------------------!
   case default
      select case (isfcl)
      case (0)
         !----- This is just to dump something in the roughness, not really used. ---------!
         patch_rough  = snowrough
         !---------------------------------------------------------------------------------!
      case default
         !----- Possibly land, and with sufficient area. -------------------------------------!
         summer_rough = max( topzo                                                         &
                           , veg_rough * veg_fracarea + soil_rough * (1.0 - veg_fracarea) )
         patch_rough  = summer_rough * (1. - snowfac) + snowrough * snowfac
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!
   end select

   return
end subroutine leaf3_roughness
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine normalises the accumulated fluxes and albedo seen by atmosphere     !
! over one full BRAMS timestep (dtlt).                                                     !
!------------------------------------------------------------------------------------------!
subroutine normal_accfluxes(m2,m3,mpat,ia,iz,ja,jz,atm_rhos,patch_area,sflux_u,sflux_v     &
                           ,sflux_w,sflux_t,sflux_r,sflux_c,albedt,rlongup)
   use leaf_coms  , only : dtl3_factor    & ! intent(in)
                         , min_patch_area ! ! intent(in)
   use mem_radiate, only : iswrtyp        & ! intent(in)
                         , ilwrtyp        ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                       , intent(in)    :: m2
   integer                       , intent(in)    :: m3
   integer                       , intent(in)    :: mpat
   integer                       , intent(in)    :: ia
   integer                       , intent(in)    :: iz
   integer                       , intent(in)    :: ja
   integer                       , intent(in)    :: jz
   real   , dimension(m2,m3)     , intent(in)    :: atm_rhos
   real   , dimension(m2,m3,mpat), intent(in)    :: patch_area
   real   , dimension(m2,m3)     , intent(inout) :: sflux_u
   real   , dimension(m2,m3)     , intent(inout) :: sflux_v
   real   , dimension(m2,m3)     , intent(inout) :: sflux_w
   real   , dimension(m2,m3)     , intent(inout) :: sflux_t
   real   , dimension(m2,m3)     , intent(inout) :: sflux_r
   real   , dimension(m2,m3)     , intent(inout) :: sflux_c
   real   , dimension(m2,m3)     , intent(inout) :: albedt
   real   , dimension(m2,m3)     , intent(inout) :: rlongup
   !----- Local variables. ----------------------------------------------------------------!
   integer                                       :: i
   integer                                       :: j
   integer                                       :: p
   real                                          :: rho_dtlt
   real                                          :: solarea
   real                                          :: solarea_i
   logical                                       :: rad_on
   !---------------------------------------------------------------------------------------!


   !----- Save the check in a logical variable to speed up the loops. ---------------------!
   rad_on = iswrtyp > 0 .or. ilwrtyp > 0
   !---------------------------------------------------------------------------------------!


   !----- Horizontal loops. ---------------------------------------------------------------!
   latloop: do j=ja,jz
      lonloop: do i=ia,iz

         solarea = patch_area(i,j,1)
         do p=2,mpat
            if (patch_area(i,j,p) >= min_patch_area) solarea = solarea + patch_area(i,j,p)
         end do
         solarea_i = 1.0 / solarea

         rho_dtlt = atm_rhos(i,j) * dtl3_factor

         sflux_u(i,j) = sflux_u(i,j) * rho_dtlt * solarea_i
         sflux_v(i,j) = sflux_v(i,j) * rho_dtlt * solarea_i
         sflux_w(i,j) = sflux_w(i,j) * rho_dtlt * solarea_i
         sflux_t(i,j) = sflux_t(i,j) * rho_dtlt * solarea_i
         sflux_r(i,j) = sflux_r(i,j) * rho_dtlt * solarea_i
         sflux_c(i,j) = sflux_c(i,j) * rho_dtlt * solarea_i

          if (rad_on) then
             albedt (i,j) = albedt (i,j) * dtl3_factor * solarea_i
             rlongup(i,j) = rlongup(i,j) * dtl3_factor * solarea_i
          end if

      end do lonloop
   end do latloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine normal_accfluxes
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine decides whether this patch should be solved or not.                 !
!------------------------------------------------------------------------------------------!
subroutine leaf3_solve_veg(ip,mzs,leaf_class,soil_rough,patch_area,veg_fracarea,veg_tai    &
                          ,sfcwater_nlev,sfcwater_mass,sfcwater_depth,initial)
   use leaf_coms , only : min_patch_area    & ! intent(in)
                        , min_sfcwater_mass & ! intent(in)
                        , tai_max           & ! intent(in)
                        , tai_min           & ! intent(in)
                        , snowfac_max       & ! intent(in)
                        , ny07_eq04_a       & ! intent(in)
                        , ny07_eq04_m       & ! intent(in)
                        , snowfac           & ! intent(inout)
                        , resolvable        ! ! intent(inout)
   use rconstants, only : wdns              & ! intent(in)
                        , fsdns             & ! intent(in)
                        , fsdnsi            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                , intent(in) :: ip
   integer                , intent(in) :: mzs
   real                   , intent(in) :: leaf_class
   real                   , intent(in) :: soil_rough
   real                   , intent(in) :: patch_area
   real                   , intent(in) :: veg_fracarea
   real                   , intent(in) :: veg_tai
   real                   , intent(in) :: sfcwater_nlev
   real   , dimension(mzs), intent(in) :: sfcwater_mass
   real   , dimension(mzs), intent(in) :: sfcwater_depth
   logical                , intent(in) :: initial
   !----- Local variables. ----------------------------------------------------------------!
   integer                             :: nveg
   integer                             :: k
   integer                             :: ksn
   real                                :: total_sfcw_depth
   real                                :: total_sfcw_mass
   real                                :: bulk_sfcw_dens
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the fraction of the canopy covered in snow.  I could not find any            !
   ! reference for the original method (commented out), so I implemented the method used   !
   ! in CLM-4, which is based on:                                                          !
   !                                                                                       !
   ! Niu, G.-Y., and Z.-L. Yang (2007), An observation-based formulation of snow cover     !
   !    fraction and its evaluation over large North American river basins,                !
   !    J. Geophys. Res., 112, D21101, doi:10.1029/2007JD008674                            !
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the area covered by snow.  This should be done every time.                   !
   !---------------------------------------------------------------------------------------!
   ksn     = nint(sfcwater_nlev)
   total_sfcw_depth = 0.0
   total_sfcw_mass  = 0.0
   do k=1,ksn
      total_sfcw_depth = total_sfcw_depth + sfcwater_depth(k)
      total_sfcw_mass  = total_sfcw_mass  + sfcwater_mass (k)
   end do
   if (total_sfcw_mass > min_sfcwater_mass) then
      bulk_sfcw_dens = max( fsdns, min( wdns, total_sfcw_mass / total_sfcw_depth))
      snowfac        = max( 0.0, min( 0.99                                                 &
                              , tanh( total_sfcw_depth                                     &
                                    / ( ny07_eq04_a * max(0.001,soil_rough)                &
                                      * (bulk_sfcw_dens * fsdnsi) ** ny07_eq04_m ) ) ) )
   else
      snowfac = 0.0
   end if
   ! snowfac = min(.99, total_sfcw_depth / max(.001,veg_height))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Now we must check whether this is the initial call or just an update.             !
   !---------------------------------------------------------------------------------------!
   if (initial) then
      !---- First call.  Decide whether the vegetation can be solved. ---------------------!
      if (ip == 1) then
         resolvable = .false.
      else
         nveg = nint(leaf_class)
         resolvable = tai_max(nveg)            >= tai_min .and.                            &
                      veg_tai * (1. - snowfac) >= tai_min
      end if
   elseif (resolvable) then
      !------------------------------------------------------------------------------------!
      !     Call in the middle of the step.  This can go only from resolvable to non-      !
      ! resolvable, in case snow or water has buried or drowned the plants.  The single    !
      ! direction is to avoid calling vegetation properties that cannot start to be        !
      ! computed in the middle of the step.                                                !
      !------------------------------------------------------------------------------------!
      nveg       = nint(leaf_class)
      resolvable = tai_max(nveg)            >= tai_min .and.                               &
                   veg_tai * (1. - snowfac) >= tai_min
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf3_solve_veg
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     This sub-routine updates the 10-day running average of the relative soil potential   !
! for phenology in the next time step.                                                     !
!------------------------------------------------------------------------------------------!
subroutine update_psibar(m2,m3,mzg,npat,ia,iz,ja,jz,dtime,soil_energy,soil_water,soil_text &
                        ,leaf_class,psibar_10d)
   use therm_lib , only : uextcm2tl ! ! subroutine
   use mem_leaf  , only : slz       & ! intent(in)
                        , dtleaf    ! ! intent(in)
   use leaf_coms , only : slzt      & ! intent(in)
                        , slmsts    & ! intent(in)
                        , slpots    & ! intent(in)
                        , slbs      & ! intent(in)
                        , slcpd     & ! intent(in)
                        , kroot     & ! intent(in)
                        , psild     & ! intent(in)
                        , psiwp     ! ! intent(in)
   use rconstants, only : wdns      & ! intent(in)
                        , day_sec   ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                           , intent(in)    :: m2
   integer                           , intent(in)    :: m3
   integer                           , intent(in)    :: mzg
   integer                           , intent(in)    :: npat
   integer                           , intent(in)    :: ia
   integer                           , intent(in)    :: iz
   integer                           , intent(in)    :: ja
   integer                           , intent(in)    :: jz
   real                              , intent(in)    :: dtime
   real   , dimension(mzg,m2,m3,npat), intent(in)    :: soil_energy
   real   , dimension(mzg,m2,m3,npat), intent(in)    :: soil_water
   real   , dimension(mzg,m2,m3,npat), intent(in)    :: soil_text
   real   , dimension    (m2,m3,npat), intent(in)    :: leaf_class
   real   , dimension    (m2,m3,npat), intent(inout) :: psibar_10d
   !----- Local variables. ----------------------------------------------------------------!
   integer                                           :: i
   integer                                           :: j
   integer                                           :: k
   integer                                           :: ip
   integer                                           :: nsoil
   integer                                           :: nveg
   real                                              :: available_water
   real                                              :: psi_layer
   real                                              :: wgpfrac
   real                                              :: soil_temp
   real                                              :: soil_fliq
   real                                              :: weight
   !---------------------------------------------------------------------------------------!


   !----- Find the weight. ----------------------------------------------------------------!
   weight = min(dtime,dtleaf) / (10. * day_sec)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Loop over all points, but skip water patches.                                      !
   !---------------------------------------------------------------------------------------!
   yloop: do j=ja,jz
      xloop: do i=ia,iz
         patchloop: do ip=2,npat
            nveg = nint(leaf_class(i,j,ip))
            available_water = 0.0
            do k = kroot(nveg),mzg
               nsoil = nint(soil_text(k,i,j,ip))

               if (nsoil /= 13) then
                  !----- Find the liquid fraction, which will scale available water. ------!
                  call uextcm2tl(soil_energy(k,i,j,ip),soil_water(k,i,j,ip)*wdns           &
                                ,slcpd(nsoil),soil_temp,soil_fliq)
                  !------------------------------------------------------------------------!



                  !----- Find the water potential of this layer. --------------------------!
                  wgpfrac         = min(soil_water(k,i,j,ip)/slmsts(nsoil), 1.0)
                  psi_layer       = slzt(k) + slpots(nsoil) / wgpfrac ** slbs(nsoil)
                  !------------------------------------------------------------------------!



                  !----- Add the contribution of this layer, based on the potential. ------!
                  available_water = available_water                                        &
                                  + max(0., (psi_layer    - psiwp(nsoil))                  &
                                          / (psild(nsoil) - psiwp(nsoil)) )                &
                                  * soil_fliq * (slz(k+1)-slz(k))
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!



            !----- Normalise the available water. -----------------------------------------!
            available_water      = available_water / abs(slz(kroot(nveg)))
            !------------------------------------------------------------------------------!



            !----- Move the moving average. -----------------------------------------------!
            psibar_10d(i,j,ip) = available_water    * weight                               &
                               + psibar_10d(i,j,ip) * (1.0 - weight)
            !------------------------------------------------------------------------------!

         end do patchloop
      end do xloop
   end do yloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine update_psibar
!==========================================================================================!
!==========================================================================================!

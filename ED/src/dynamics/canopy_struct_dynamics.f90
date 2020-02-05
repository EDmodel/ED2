!==========================================================================================!
!==========================================================================================!
!    This module contains some subroutines used to derive turbulence-related properties of !
! the canopy air space.                                                                    !
!------------------------------------------------------------------------------------------!
module canopy_struct_dynamics

   contains
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     The following subroutine initialises the canopy turbulence variables.  Most of    !
   ! these variables can be initialised with zeroes or dummy values, without the need to   !
   ! run a formal calculation.  This should simplify the code.                             !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_turbulence_init(cpoly,isi,ipa)
      use ed_state_vars    , only : polygontype          & ! structure
                                  , sitetype             & ! structure
                                  , patchtype            ! ! structure
      use met_driver_coms  , only : met_driv_state       ! ! structure
      use canopy_air_coms  , only : isfclyrm             & ! intent(in)
                                  , ugbmin               & ! intent(in)
                                  , ustmin               & ! intent(in)
                                  , vh2vr                & ! intent(in)
                                  , vh2dh                & ! intent(in)
                                  , exar                 & ! intent(in)
                                  , tprandtl             & ! intent(in)
                                  , zd98_a               & ! intent(in)
                                  , zd98_b               & ! intent(in)
                                  , zd98_emax            ! ! intent(in)
      use consts_coms      , only : vonk                 & ! intent(in)
                                  , cpdry                & ! intent(in)
                                  , cph2o                & ! intent(in)
                                  , kin_visc0            & ! intent(in)
                                  , dkin_visc            & ! intent(in)
                                  , t00                  ! ! intent(in)
      use soil_coms        , only : snow_rough           & ! intent(in)
                                  , soil_rough           ! ! intent(in)
      use rk4_coms         , only : tiny_offset          ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype)   , target      :: cpoly         ! Current polygon
      integer             , intent(in)  :: isi           ! Site index
      integer             , intent(in)  :: ipa           ! Patch index
      !----- Pointers ---------------------------------------------------------------------!
      type(sitetype)      , pointer     :: csite         ! Current site
      type(patchtype)     , pointer     :: cpatch        ! Current patch
      type(met_driv_state), pointer     :: cmet          ! Current meteorology
      !----- Local variables --------------------------------------------------------------!
      integer      :: ico           !> Cohort index                             [      ---]
      integer      :: ipft          !> PFT alias                                [      ---]
      real         :: can_cp        !> Canopy air space specific heat           [   J/kg/K]
      real         :: uh            !> Wind speed at the canopy top (z=h)       [      m/s]
      real         :: factv         !> Wind-dependent term for old rasveg       [      ---]
      real         :: aux           !> Aux. variable                            [      ---]
      real         :: extinct_half  !> Wind extinction coefficient at half lyr  [      ---]
      real         :: extinct_full  !> Full Wind extinction coefficient         [      ---]
      real         :: snowfac_can   !> fraction of canopy covered in snow       [      ---]
      real         :: kin_visc      !> Kinematic viscosity                      [     m2/s]
      real         :: z0hoz0m       !> Heat:momentum roughness ratio            [      ---]
      real         :: z0moz0h       !> Momentum:heat roughness ratio            [      ---]
      real         :: zstar         !> Ref. elevation above displacement height [        m]
      real         :: zoz0m         !> zstar / roughness (momentum)             [      ---]
      real         :: lnzoz0m       !> log of zoz0m                             [      ---]
      real         :: zoz0h         !> zstar / roughness (heat)                 [      ---]
      real         :: lnzoz0h       !> log of zoz0h                             [      ---]
      real         :: c3            !> coefficient for CAS conductance          [      ---]
      real(kind=8) :: can_temp8     !> CAS temperature (double precision)       [        K]
      real(kind=8) :: can_rhos8     !> CAS air density (double precision)       [    kg/m3]
      real(kind=8) :: can_cp8       !> CAS specific capacity (double precision) [   J/kg/K]
      real(kind=8) :: veg_wind8     !> Wind speed at cohort height (dbl. prec.) [      m/s]
      real(kind=8) :: leaf_temp8    !> Leaf temperature (double precision)      [        K]
      real(kind=8) :: wood_temp8    !> Wood temperature (double precision)      [        K]
      real(kind=8) :: veg_gbh8      !> Heat  conductance                        [ J/K/m2/s]
      real(kind=8) :: veg_gbw8      !> Water conductance                        [  kg/m2/s]
      real(kind=8) :: veg_reynolds8 !> Reynolds number                          [      ---]
      real(kind=8) :: veg_grashof8  !> Grashop number                           [      ---]
      real(kind=8) :: veg_nussfree8 !> Nusselt number (free convection)         [      ---]
      real(kind=8) :: veg_nussforc8 !> Nusselt number (forced convection)       [      ---]
      !----- External functions. ----------------------------------------------------------!
      real(kind=4), external :: sngloff !> Safe dble->sngl precision converter
      !------------------------------------------------------------------------------------!


      !----- Assign some pointers. --------------------------------------------------------!
      csite  => cpoly%site(isi)
      cmet   => cpoly%met(isi)
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the specific heat at constant pressure for this canopy air space.         !
      !------------------------------------------------------------------------------------!
      can_cp     = (1.0 - csite%can_shv(ipa)) * cpdry + csite%can_shv(ipa) * cph2o
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the fraction of the canopy covered in snow (original snowfac function)    !
      ! I think canopy roughness may need to be re-thought, but this was necessary         !
      ! for turbulence & CO2 mixing to not occasionally fail sanity checks in young patches!
      ! (CR)
      ! Calculation of surface roughness should use the planar snow fraction, as this      !
      ! should be entailing an area weighted average.  If there is instability, I argue we !
      ! should address how the horizontal planar fraction is calculated, not use the       !
      ! vertical as a surrogate. (RGK) I Will leave as is for now. But this needs to be    !
      ! fixed.                                                                             !
      !------------------------------------------------------------------------------------!
      snowfac_can     = min(0.99,csite%total_sfcw_depth(ipa)/csite%veg_height(ipa))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find roughness based on the combination of bare ground, snowpack and           !
      ! vegetation (the latter two may or may not exist).                                  !
      !------------------------------------------------------------------------------------!
      if (cpatch%ncohorts == 0) then

         !----- Calculate the surface roughness inside the canopy. ------------------------!
         csite%rough       (ipa) = soil_rough * (1.0 - snowfac_can)                        &
                                 + snow_rough * snowfac_can
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Find the roughness as the average between the bare ground and vegetated     !
         ! ground, and apply the snow cover to further scale it.  The weighting factors    !
         ! are the fraction of open canopy and the fraction of the canopy buried in snow.  !
         !---------------------------------------------------------------------------------!
         csite%rough(ipa) = snow_rough * snowfac_can                                       &
                          + ( soil_rough           * csite%opencan_frac(ipa)               &
                            + csite%veg_rough(ipa) * (1.0 - csite%opencan_frac(ipa)) )     &
                          * (1.0 - snowfac_can)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------ Initialise displacement height as a simple function of net roughness. -------!
      csite%veg_displace(ipa) = vh2dh * csite%rough(ipa) / vh2vr
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the characteristic scales (a.k.a. stars).  Assume equilibrium conditions  !
      ! and set all but ustar to zero.  Ustar cannot be zero, assume the minimum           !
      ! acceptable value.                                                                  !
      !------------------------------------------------------------------------------------!
      csite%ustar (ipa) = ustmin
      csite%qstar (ipa) = 0.
      csite%tstar (ipa) = 0.
      csite%cstar (ipa) = 0.
      csite%zeta  (ipa) = 0.
      csite%ribulk(ipa) = 0.
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the ratio between momentum and heat roughness scales.                      !
      !------------------------------------------------------------------------------------!
      select case (isfclyrm)
      case (4)
         !---------------------------------------------------------------------------------!
         !      Use previous u* to estimate the ratio between roughness scales for         !
         ! momentum (z0m) and roughness scale for enthalpy/vapour/CO2 (z0h), using the     !
         ! parametrisation from ZD98.                                                      ! 
         !---------------------------------------------------------------------------------!
         kin_visc = kin_visc0 * ( 1.0 + dkin_visc * ( csite%can_temp(ipa) - t00 ) )
         z0hoz0m  = exp(- min( zd98_emax                                                   &
                             , zd98_a * (csite%ustar(ipa)*csite%rough(ipa)/kin_visc)       &
                                      ** zd98_b                                      ) )
         z0moz0h  = 1.0 / z0hoz0m
         !---------------------------------------------------------------------------------!
      case default
         !---------------------------------------------------------------------------------!
         !      Assume the same roughness scheme.  Mind that this could be off by multiple !
         ! orders of magnitude.                                                            ! 
         !---------------------------------------------------------------------------------!
         z0hoz0m  = 1.0
         z0moz0h  = 1.0 / z0hoz0m
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !----- Find bare-ground conductance assuming neutral conditions. --------------------!
      zstar             = cmet%geoht - csite%veg_displace(ipa)
      zoz0m             = (cmet%geoht-csite%veg_displace(ipa)) / csite%rough(ipa)
      lnzoz0m           = log(zoz0m)
      zoz0h             = z0moz0h * zoz0m
      lnzoz0h           = log(zoz0h)
      c3                = vonk / (tprandtl * log(zoz0h))
      csite%ggbare(ipa) = c3 * csite%ustar(ipa)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use a simplified initalisation for conductances.                               !
      !------------------------------------------------------------------------------------!
      if (cpatch%ncohorts == 0 .or. snowfac_can >= .9 .or. csite%opencan_frac(ipa) > .999) &
      then
         csite%ggveg(ipa) = 0.0
         csite%ggnet(ipa) = csite%ggbare(ipa)
      else
         factv            = log((cmet%geoht - csite%veg_displace(ipa)) / csite%rough(ipa)) &
                          / (vonk * vonk * cmet%vels)
         aux              = exp(exar * (1. - (csite%veg_displace(ipa) + csite%rough(ipa))  &
                                           / csite%veg_height(ipa)) )
         csite%ggveg(ipa) = (exar * (csite%veg_height(ipa)-csite%veg_displace(ipa)))       &
                          / (factv * csite%veg_height(ipa) * (exp(exar) - aux))
         csite%ggnet(ipa) = csite%ggbare(ipa) * csite%ggveg(ipa)                           &
                          / ( csite%ggveg(ipa) + (1. - csite%opencan_frac(ipa))            &
                                               * csite%ggbare(ipa))
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Convert patch-level variables to double precision (they will be used to        !
      ! compute cohort conductances).                                                      !
      !------------------------------------------------------------------------------------!
      can_temp8 = dble(csite%can_temp(ipa))
      can_rhos8 = dble(csite%can_rhos(ipa))
      can_cp8   = dble(can_cp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Assign a simplified initial wind profile, based on the exponetial decay with    !
      ! LAI, as in Leuning (1995), with a modification to account for crown clumping.      !
      !------------------------------------------------------------------------------------!
      uh = csite%ustar(ipa) * lnzoz0m / vonk
      do ico=1,cpatch%ncohorts
         ipft = cpatch%pft(ico)

         !----- Find the extinction coefficients. -----------------------------------------!
         extinct_half = cpatch%crown_area(ico)                                             &
                      * exp(- 0.25 * cpatch%lai(ico) / cpatch%crown_area(ico))             &
                      + (1. - cpatch%crown_area(ico))
         extinct_full = cpatch%crown_area(ico)                                             &
                      * exp(- 0.50 * cpatch%lai(ico) / cpatch%crown_area(ico))             &
                      + (1. - cpatch%crown_area(ico))
         !---------------------------------------------------------------------------------!


         !----- Assume that wind is at the middle of the thin crown. ----------------------!
         cpatch%veg_wind(ico) = max(ugbmin, uh * extinct_half)
         uh                   = uh * extinct_full
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Convert a few cohort variables to double precision (needed for              !
         ! conductances).                                                                  !
         !---------------------------------------------------------------------------------!
         veg_wind8  = dble(cpatch%veg_wind (ico))
         leaf_temp8 = dble(cpatch%leaf_temp(ico))
         wood_temp8 = dble(cpatch%wood_temp(ico))
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Find the aerodynamic conductances for heat and water at the leaf and branch- !
         ! wood boundary layers.                                                           !
         !---------------------------------------------------------------------------------!
         call leaf_aerodynamic_conductances8(ipft,veg_wind8,leaf_temp8,can_temp8           &
                                            ,can_rhos8,can_cp8,veg_gbh8,veg_gbw8           &
                                            ,veg_reynolds8,veg_grashof8,veg_nussfree8      &
                                            ,veg_nussforc8 )
         cpatch%leaf_gbh(ico) = sngloff(veg_gbh8,tiny_offset)
         cpatch%leaf_gbw(ico) = sngloff(veg_gbw8,tiny_offset)
         aux = 1. / cpatch%leaf_gbw(ico)
         call wood_aerodynamic_conductances8(ipft,veg_wind8,wood_temp8,can_temp8           &
                                           ,can_rhos8,can_cp8,veg_gbh8,veg_gbw8            &
                                           ,veg_reynolds8,veg_grashof8,veg_nussfree8       &
                                           ,veg_nussforc8)
         cpatch%wood_gbh(ico) = sngloff(veg_gbh8,tiny_offset)
         cpatch%wood_gbw(ico) = sngloff(veg_gbw8,tiny_offset)
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      return
   end subroutine canopy_turbulence_init
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This is the double precision version of the above sub-routine.                    !
   !                                                                                       !
   !     The following subroutine calculates the turbulent transfer of scalar quantities   !
   ! using Monin Obukhov Similarity Theory under a few different flavors.  The basic       !
   ! assumptions are that the ground surface and canopy space is modeled as a rough wall   !
   ! boundary above which exists a dynamic sublayer, above which exists the free atmo-     !
   ! sphere.  We assume a logarithmic decay of wind speed over the dynamic sublayer, which !
   ! is due to the assumption of constant shear (ustar) in the dynamic sublayer.           !
   !                                                                                       !
   ! 0.  We apply the wind profile based on the similarity theory for the tallest cohort,  !
   !     with further reduction due to the presence of obstacles (trees), similar to       !
   !     Leuning et al. (1995) but considering a finite crown area.  The aerodynamic       !
   !     resistance from ground to canopy air is found in a similar way as in LEAF-3.      !
   !     Displacement height and roughness are linear functions of the vegetation height,  !
   !     where vegetation height is the weighted average of cohort heights, the weights    !
   !     being the basal area.                                                             !
   !     References:                                                                       !
   !                                                                                       !
   !     Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf          !
   !       nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to !
   !       canopies.  Plant, Cell and Environ., 18, 1183-1200.                             !
   !                                                                                       !
   ! 1. This is similar to option 0, where displacement height and rougness are linear     !
   !    functions of mean canopy height, where mean canopy height is defined by the        !
   !    dead-biomass weighted average of canopy height.  Wind speed and diffusivity decay  !
   !    exponentially in the canopy on parameter alpha.  HOWEVER, if the windspeed at      !
   !    reference height is less than the canopy height (WHICH IS REALLY INNAPROPRIATE TO  !
   !    BEGIN WITH, BUT OH WELL), then revert to a calculation of ustar that uses no       !
   !    displacement height.  After ustar is calculated, then use displacement height in   !
   !    the calculation of the eddy length scale at the canopy top to find diffusivity at  !
   !    canopy top. THIS IS THE METHOD USED IN LEAF3 AND ED2.0.                            !
   !                                                                                       !
   ! 2. Follow the methods of Massman 1997.  The canopy top is taken as the height of the  !
   !    tallest cohort in the patch. The displacement height is determined by a center of  !
   !    pressure method. The wind speed inside the canopy decays exponentially on          !
   !    cumulative foliar drag elements, where the equation is scaled such that it solves  !
   !    for the balance of integrated canopy shear with dissipation.  To close the         !
   !    equations M97 assumes that the surface drag (U(h)/Ustar) is a function of the      !
   !    cumulative foliar drag surfaces.  With this, we can explicitly solve for roughness !
   !    length. Assumptions are made about the profile of canopy drag as a function of     !
   !    cohort attributes, please see the technical manual by RGK2009 for more info.       !
   !    Please also refer to the original paper by Massman for information regarding the   !
   !    basic principles of the closure scheme.                                            !
   !                                                                                       !
   !    Massman, W.J., 1997: An analytical one-Dimensional model of momentum transfer by   !
   !        vegetation of arbitrary structure. Boundary-Layer Meteorol., 83, 407-421.      !
   !    Wohlfahrt, G., and A. Cernusca, 2002: Momentum transfer by a mountain meadow       !
   !        canopy: a simulation analysis based on Massman's (1997) model.  Boundary-Layer !
   !        Meteorol., 103, 391-407.                                                       !
   !                                                                                       !
   ! 3.  This is related to option 2, but using a second-order clousure.                   !
   !                                                                                       !
   !    Massman, W. J., and J. C. Weil, 1999: An analytical one-dimension second-order     !
   !        closure model turbulence statistics and the Lagrangian time scale within and   !
   !        above plant canopies of arbitrary structure.  Boundary-Layer Meteorol., 91,    !
   !        81-107.                                                                        !
   !                                                                                       !
   ! 4.  This option is almost the same as option 1, except that the ground conductance is !
   !     found following the CLM technical note, equations (5.98-5.100).                   !
   !                                                                                       !
   !     Oleson, K. W., et al.; Technical description of the community land model (CLM)    !
   !        NCAR Technical Note NCAR/TN-461+STR, Boulder, CO, May 2004.                    !
   !                                                                                       !
   !     Ultimately, this routine solves for the resistance of water vapor and sensible    !
   ! heat from the soil surface to canopy air space and from the leaf surfaces to canopy   !
   ! air space.                                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_turbulence8(csite,initp,ipa,ibuff)
      use ed_state_vars    , only : polygontype          & ! structure
                                  , sitetype             & ! structure
                                  , patchtype            ! ! structure
      use rk4_coms         , only : rk4patchtype         & ! structure
                                  , rk4site              & ! intent(in)
                                  , tiny_offset          & ! intent(in)
                                  , ibranch_thermo       ! ! intent(in)
      use grid_coms        , only : nzg                  ! ! intent(in)
      use canopy_air_coms  , only : icanturb             & ! intent(in), can. turb. scheme
                                  , f_bndlyr_init        & ! intent(in)
                                  , ugbmin8              & ! intent(in)
                                  , vh2vr8               & ! intent(in)
                                  , vh2dh8               & ! intent(in)
                                  , veg_height_min8      & ! intent(in)
                                  , exar8                & ! intent(in)
                                  , cdrag08              & ! intent(in)
                                  , cdrag18              & ! intent(in)
                                  , cdrag28              & ! intent(in)
                                  , cdrag38              & ! intent(in)
                                  , c1_m978              & ! intent(in)
                                  , c2_m978              & ! intent(in)
                                  , c3_m978              & ! intent(in)
                                  , kvwake8              & ! intent(in)
                                  , alpha_m97_8          & ! intent(in)
                                  , alpha_mw99_8         & ! intent(in)
                                  , infunc_8             & ! intent(in)
                                  , gamma_mw99_8         & ! intent(in)
                                  , nu_mw99_8            & ! intent(in)
                                  , cs_dense08           & ! intent(in)
                                  , gamma_clm48          & ! intent(in)
                                  , tprandtl8            & ! intent(in)
                                  , dpsimdzeta8          ! ! function
      use canopy_layer_coms, only : ncanlyr              & ! intent(in)
                                  , dzcan8               & ! intent(in)
                                  , zztop0i8             & ! intent(in)
                                  , ehgti8               & ! intent(in)
                                  , zztop8               & ! intent(in)
                                  , zzbot8               & ! intent(in)
                                  , zzmid8               & ! intent(in)
                                  , canstr               &
                                  , zero_canopy_layer    ! ! subroutine
      use consts_coms      , only : cpdry                & ! intent(in)
                                  , vonk8                & ! intent(in)
                                  , grav8                & ! intent(in)
                                  , t008                 & ! intent(in)
                                  , epim18               & ! intent(in)
                                  , sqrt2o28             & ! intent(in)
                                  , srthree8             & ! intent(in)
                                  , onethird8            & ! intent(in)
                                  , twothirds8           & ! intent(in)
                                  , th_diff08            & ! intent(in)
                                  , dth_diff8            & ! intent(in)
                                  , lnexp_min8           & ! intent(in)
                                  , lnexp_max8           ! ! intent(in)
      use soil_coms        , only : snow_rough8          & ! intent(in)
                                  , soil_rough8          ! ! intent(in)
      use pft_coms         , only : is_grass             ! ! intent(in)
      use allometry        , only : h2crownbh            & ! function
                                  , size2bl              ! ! function
      use ed_misc_coms     , only : igrass               ! ! intent(in)
      use phenology_coms   , only : elongf_min           ! ! intent(in)
      use physiology_coms  , only : gbh_2_gbw            ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)     , target     :: csite         ! Current site
      type(rk4patchtype) , target     :: initp         ! Current integrator
      integer            , intent(in) :: ipa           ! Patch loop
      integer            , intent(in) :: ibuff         ! Multithread ID
      !----- Pointers ---------------------------------------------------------------------!
      type(patchtype)    , pointer    :: cpatch        ! Current patch
      !----- Local variables --------------------------------------------------------------!
      integer        :: ico          ! Cohort loop
      integer        :: ipft         ! PFT alias
      integer        :: k            ! Elevation index
      integer        :: ksn          ! Number of temporary sfc. water layers    [      ---]
      integer        :: zcan         ! Index of canopy top elevation
      integer        :: kafull       ! First layer fully occupied by crown
      integer        :: kzfull       ! Last layer fully occupied by crown
      integer        :: kapartial    ! First layer partially occupied by crown
      integer        :: kzpartial    ! Last layer partially occupied by crown
      integer        :: zels         ! Index of roughness height
      integer        :: nalpha_fails ! Counter for number of failed attemps     [      ---]
      integer        :: kk           ! Counter                                  [      ---]
      logical        :: stable       ! Stable canopy air space
      logical        :: acomp        ! Flag to check for convergence            [      T|F]
      real(kind=8)   :: rasveg       ! Resistance of vegetated ground           [      s/m]
      real(kind=8)   :: atm_thetav   ! Free atmosphere virtual potential temp.  [        K]
      real(kind=8)   :: can_thetav   ! Free atmosphere virtual potential temp.  [        K]
      real(kind=8)   :: ldga_bk      ! Cumulative zeta function                 [      ---]
      real(kind=8)   :: lyrhalf      ! Half the contrib. of this layer to zeta  [      1/m]
      real(kind=8)   :: sigmakh      ! Kh coefficient at z=h                    [        m]
      real(kind=8)   :: K_top        ! Diffusivity at canopy top z=h            [     m2/s]
      real(kind=8)   :: kdiff        ! Diffusivity                              [     m2/s]
      real(kind=8)   :: surf_rough   ! Roughness length of the bare ground
                                     !     at canopy bottom                     [        m]
      real(kind=8)   :: uh           ! Wind speed at the canopy top (z=h)       [      m/s]
      real(kind=8)   :: factv        ! Wind-dependent term for old rasveg
      real(kind=8)   :: aux          ! Aux. variable
      real(kind=8)   :: ustarouh     ! The ratio of ustar over u(h)             [      ---]
      real(kind=8)   :: nn           ! In-canopy wind attenuation scal. param.  [      ---]
      real(kind=8)   :: waiuse       ! Wood area index                          [    m2/m2]
      real(kind=8)   :: htopcrown    ! height at the top of the crown           [        m]
      real(kind=8)   :: hmidcrown    ! Height at the middle of the crown        [        m]
      real(kind=8)   :: hbotcrown    ! Height at the bottom of the crown        [        m]
      real(kind=8)   :: htop         ! Height of the topmost layer              [        m]
      real(kind=8)   :: zetatop      ! Dimensionless height at the topmost lyr. [      ---]
      real(kind=8)   :: d0ohgt       ! d0/height                                [      ---]
      real(kind=8)   :: z0ohgt       ! z0/height                                [      ---]
      real(kind=8)   :: ladcohort    ! Leaf Area Density of this cohort         [    m2/m3]
      real(kind=8)   :: extinct_half ! Wind extinction coefficient at half lyr  [      ---]
      real(kind=8)   :: extinct_full ! Full Wind extinction coefficient         [      ---]
      real(kind=8)   :: elenscale    ! Eddy lenght scale                        [        m]
      real(kind=8)   :: alpha_eq10   ! Alpha (may be tweaked for convergence)   [      ---]
      real(kind=8)   :: lam          ! Mixed term from MW99                     [      ---]
      real(kind=8)   :: b1_mw99      ! B1 term from MW99                        [      ---]
      real(kind=8)   :: nddfun       ! Normalised drag density function         [      ---]
      real(kind=8)   :: ure          ! Wind speed averaged of sfc mixing length [      m/s]
      real(kind=8)   :: sigstar      ! ustar norm. canopy velocity variance     [      ---]
      real(kind=8)   :: sigstar3     ! Cubic of sigstar                         [      ---]
      real(kind=8)   :: sigcomm      ! Common term for variance calculation     [      m/s]
      real(kind=8)   :: sigma_uou2   ! Square of (sigma u / u)                  [      ---]
      real(kind=8)   :: sigma_vou2   ! Square of (sigma v / u)                  [      ---]
      real(kind=8)   :: sigma_wou2   ! Square of (sigma w / u)                  [      ---]
      real(kind=8)   :: turbi        ! Mean turbulent intensity                 [      m/s]
      real(kind=8)   :: can_reynolds ! Reynolds number of the sfc. mixing layer [      ---]
      real(kind=8)   :: ground_temp  ! Ground temperature                       [      ---]
      real(kind=8)   :: stab_clm4    ! Stability parameter (CLM4, eq. 5.104)    [      ---]
      logical        :: dry_grasses  ! Flag to check whether LAI+WAI is zero    [      ---]
      real(kind=8)   :: tai_drygrass ! TAI for when a grass-only patch is dry   [    m2/m2]
      real(kind=8)   :: c3_lad       ! c3 * lad for estimating drag coefficient [      ---]
      real(kind=8)   :: c3_cumldrag  ! c3 * cumulative drag                     [      ---]
      real(kind=8)   :: snowfac_can  ! percent vertical canopy covered in snow
      !------ External procedures ---------------------------------------------------------!
      real(kind=8), external :: cbrt8    ! Cubic root that works for negative numbers
      real(kind=4), external :: sngloff  ! Safe double -> simple precision.
      !------------------------------------------------------------------------------------!

      !----- Assign some pointers. --------------------------------------------------------!
      cpatch=>csite%patch(ipa)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the fraction of the canopy covered in snow (original snowfac function)    !
      !     I think canopy roughness may need to be re-thought, but this was necessary     !
      ! for turbulence & CO2 mixing to not occasionally fail sanity checks in young patches!
      !------------------------------------------------------------------------------------!
      snowfac_can     = min(9.9d-1,initp%total_sfcw_depth/initp%veg_height)
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     Find the virtual potential temperatures and decide whether the canopy air is   !
      ! stable or not.                                                                     !
      !------------------------------------------------------------------------------------!
      atm_thetav = rk4site%atm_theta * (1.d0 + epim18 * rk4site%atm_shv)
      can_thetav = initp%can_theta * (1.d0 + epim18 * initp%can_shv  )
      stable     = atm_thetav >= can_thetav

      !------------------------------------------------------------------------------------!
      !     If there is no vegetation in this patch, then we apply turbulence to bare      !
      ! soil, no d0 and exit.                                                              !
      !------------------------------------------------------------------------------------!
      if (cpatch%ncohorts == 0) then

         !----- Calculate the surface roughness and displacement height. ------------------!
         initp%rough        = soil_rough8 *(1.d0 - snowfac_can)                            &
                            + snow_rough8 * snowfac_can
         initp%veg_displace = vh2dh8 * initp%rough / vh2vr8
         !---------------------------------------------------------------------------------!



         !----- Find the characteristic scales (a.k.a. stars). ----------------------------!
         call ed_stars8(rk4site%atm_theta,initp%atm_enthalpy,rk4site%atm_shv               &
                       ,rk4site%atm_co2,initp%can_theta ,initp%can_enthalpy,initp%can_shv  &
                       ,initp%can_co2,initp%can_temp,rk4site%geoht,initp%veg_displace      &
                       ,rk4site%atm_ustar,initp%vels,initp%rough,initp%ustar,initp%tstar   &
                       ,initp%estar,initp%qstar,initp%cstar,initp%zeta,initp%ribulk        &
                       ,initp%ggbare)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      This patch is empty, so we can't define a conductance for vegetated        !
         ! grounds.  Assign it to be zero, and set the net conductance to be the bare      !
         ! ground.                                                                         !
         !---------------------------------------------------------------------------------!
         initp%ggveg = 0.d0
         initp%ggnet = initp%ggbare
         !---------------------------------------------------------------------------------!

         return
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Reset scratch variables in canopy_layer_coms.                                  !
      !------------------------------------------------------------------------------------!
      call zero_canopy_layer('canopy_turbulence8',canstr(ibuff))
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     In case we do have cohorts, choose which method we use to compute the          !
      ! resistance.                                                                        !
      !------------------------------------------------------------------------------------!
      select case (icanturb)
      !------------------------------------------------------------------------------------!
      ! 0.  LEAF-3 Case: This approach is very similar to older implementations of ED-2,   !
      !     and it is very similar to LEAF-3, and to option 1, except that option 1        !
      !     computes the vegetated ground conductance and wind profile differently.        !
      ! 4.  CLM-4 Case: This is the same as option 0., except that the ground conductance  !
      !     is found using equations (5.99-5.104) from CLM technical note (except the bare !
      !     ground, which is still found using the similarity theory).                     !
      !------------------------------------------------------------------------------------!
      case (0,4)

         !---------------------------------------------------------------------------------!
         !     The roughness is found by combining two weighted averages.  The first one   !
         ! checks the fraction of the patch that has closed canopy, and averages between   !
         ! soil and vegetation roughness.  The other is the fraction of the vegetation     !
         ! that is covered in snow.                                                        !
         !---------------------------------------------------------------------------------!
         initp%rough = snow_rough8 * snowfac_can                                           &
                     + ( soil_rough8     * initp%opencan_frac                              &
                       + initp%veg_rough * (1.d0 - initp%opencan_frac) )                   &
                     * (1.d0 - snowfac_can)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Get ustar for the ABL, assume it is a dynamic shear layer that generates a !
         ! logarithmic profile of velocity.                                                !
         !---------------------------------------------------------------------------------!
         call ed_stars8(rk4site%atm_theta,initp%atm_enthalpy,rk4site%atm_shv               &
                       ,rk4site%atm_co2,initp%can_theta ,initp%can_enthalpy,initp%can_shv  &
                       ,initp%can_co2,initp%can_temp,rk4site%geoht,initp%veg_displace      &
                       ,rk4site%atm_ustar,initp%vels,initp%rough,initp%ustar,initp%tstar   &
                       ,initp%estar,initp%qstar,initp%cstar,initp%zeta,initp%ribulk        &
                       ,initp%ggbare)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the aerodynamic resistance due to vegetated ground.                    !
         !---------------------------------------------------------------------------------!
         if (snowfac_can < 9.d-1) then
            select case (icanturb)
            case (0)
               !----- LEAF-3 method. ------------------------------------------------------!
               factv        = log((rk4site%geoht-initp%veg_displace) / initp%rough)        &
                            / (vonk8 * vonk8 * initp%vels)
               aux          = exp(exar8 * (1.d0 - (initp%veg_displace + initp%rough)       &
                                                / initp%veg_height) )
               initp%ggveg  = (exar8 * (initp%veg_height - initp%veg_displace))            &
                            / (factv * initp%veg_displace * (exp(exar8) - aux))
            case (4)
               !---------------------------------------------------------------------------!
               !     Similar to CLM-4, although we use the ground conductance for bare     !
               ! ground from the similarity theory and the fraction of open canopy to      !
               ! determine the weights.  The conductance for dense canopy is simply        !
               ! Cs_dense * u*, but Cs_dense depends on the stability.                     !
               !---------------------------------------------------------------------------!
               !----- Find which ground temperature to use. -------------------------------!
               ksn = csite%nlev_sfcwater(ipa)
               if (ksn == 0) then
                  ground_temp = csite%soil_tempk(nzg,ipa)
               else
                  ground_temp = csite%sfcwater_tempk(ksn,ipa)
               end if
               !----- Stability factor (0 when it is unstable). ---------------------------!
               stab_clm4 = max( 0.d0, min(1.d1, 2.d0 * grav8 * initp%veg_height            &
                                              * ( initp%can_temp - ground_temp )           &
                                              / ( ( initp%can_temp + ground_temp)          &
                                              * initp%ustar * initp%ustar ) ) )
               !----- The "veg" conductance is equivalent to CLM4 dense canopy. -----------!
               initp%ggveg = cs_dense08 * initp%ustar / (1.d0 + gamma_clm48 * stab_clm4)
            end select
         else
            initp%ggveg = 0.d0
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the wind profile.  This is done by scaling the wind at the top of the  !
         ! canopy in a method based on Leuning et al. (1995), but with some important      !
         ! differences, depending on the crown model chosen.                               !
         !---------------------------------------------------------------------------------!
         !----- Find the wind at the top of the canopy. -----------------------------------!
         uh = reduced_wind8(initp%ustar,initp%zeta,initp%ribulk,rk4site%geoht              &
                           ,initp%veg_displace,dble(cpatch%hite(1)),initp%rough)

         !---------------------------------------------------------------------------------!
         !     In this version we still base ourselves on the Leuning et al. (1995) model, !
         ! assuming each cohort to be on top of each other (i.e. no ties).  We also assume !
         ! extinction to be limited to the finite crown area, so the bottom can get a bit  !
         ! windier.                                                                        !
         !---------------------------------------------------------------------------------!
         do ico=1,cpatch%ncohorts
            !----- Find the extinction coefficients. --------------------------------------!
            extinct_half = initp%crown_area(ico)                                           &
                         * exp(- 2.5d-1 * initp%lai(ico) / initp%crown_area(ico))          &
                         + (1.d0 - initp%crown_area(ico))
            extinct_full = initp%crown_area(ico)                                           &
                         * exp(- 5.d-1 * initp%lai(ico) / initp%crown_area(ico))           &
                         + (1.d0 - initp%crown_area(ico))
            !----- Assume that wind is at the middle of the thin crown. -------------------!
            initp%veg_wind(ico) = max(ugbmin8, uh * extinct_half)
            uh                  = uh * extinct_full
         end do
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !   Find the aerodynamic conductances.                                            !
         !---------------------------------------------------------------------------------!
         do ico=1,cpatch%ncohorts

            ipft  = cpatch%pft(ico)


            !------ Leaf boundary layer conductance. --------------------------------------!
            if (initp%leaf_resolvable(ico)) then
               !---------------------------------------------------------------------------!
               !    Find the aerodynamic conductances for heat and water at the leaf       !
               ! boundary layer.                                                           !
               !---------------------------------------------------------------------------!
               call leaf_aerodynamic_conductances8(ipft,initp%veg_wind(ico)                &
                                                  ,initp%leaf_temp(ico),initp%can_temp     &
                                                  ,initp%can_rhos                          &
                                                  ,initp%can_cp                            &
                                                  ,initp%leaf_gbh(ico),initp%leaf_gbw(ico) &
                                                  ,initp%leaf_reynolds(ico)                &
                                                  ,initp%leaf_grashof(ico)                 &
                                                  ,initp%leaf_nussfree(ico)                &
                                                  ,initp%leaf_nussforc(ico) )
               cpatch%leaf_gbh(ico) = sngloff(initp%leaf_gbh(ico),tiny_offset)
               cpatch%leaf_gbw(ico) = sngloff(initp%leaf_gbw(ico),tiny_offset)
               !---------------------------------------------------------------------------!
            else
               !----- Set numbers to zero. ------------------------------------------------!
               initp%leaf_reynolds(ico) = 0.d0
               initp%leaf_grashof (ico) = 0.d0
               initp%leaf_nussfree(ico) = 0.d0
               initp%leaf_nussforc(ico) = 0.d0
               !---------------------------------------------------------------------------!


               !----- Leaf conductances cannot be zero.  Set to non-zero defaults. --------!
               cpatch%leaf_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
               cpatch%leaf_gbh(ico) = cpatch%leaf_gbw(ico) / gbh_2_gbw * cpdry
               initp%leaf_gbw (ico) = dble(cpatch%leaf_gbw(ico))
               initp%leaf_gbh (ico) = dble(cpatch%leaf_gbh(ico))
               !---------------------------------------------------------------------------!
            end if
            !------ Wood boundary layer conductance. --------------------------------------!
            if (initp%wood_resolvable(ico)) then
               !---------------------------------------------------------------------------!
               !    Find the aerodynamic conductances for heat and water at the leaf       !
               ! boundary layer.                                                           !
               !---------------------------------------------------------------------------!
               call wood_aerodynamic_conductances8(ipft                                    &
                                                  ,initp%veg_wind(ico)                     &
                                                  ,initp%wood_temp(ico),initp%can_temp     &
                                                  ,initp%can_rhos                          &
                                                  ,initp%can_cp                            &
                                                  ,initp%wood_gbh(ico),initp%wood_gbw(ico) &
                                                  ,initp%wood_reynolds(ico)                &
                                                  ,initp%wood_grashof(ico)                 &
                                                  ,initp%wood_nussfree(ico)                &
                                                  ,initp%wood_nussforc(ico) )
               cpatch%wood_gbh(ico) = sngloff(initp%wood_gbh(ico),tiny_offset)
               cpatch%wood_gbw(ico) = sngloff(initp%wood_gbw(ico),tiny_offset)
               !---------------------------------------------------------------------------!
            else
               !----- Set numbers to zero. ------------------------------------------------!
               initp%wood_reynolds(ico) = 0.d0
               initp%wood_grashof (ico) = 0.d0
               initp%wood_nussfree(ico) = 0.d0
               initp%wood_nussforc(ico) = 0.d0
               !---------------------------------------------------------------------------!


               !----- Leaf conductances cannot be zero.  Set to non-zero defaults. --------!
               cpatch%wood_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
               cpatch%wood_gbh(ico) = cpatch%wood_gbw(ico) / gbh_2_gbw * cpdry
               initp%wood_gbw (ico) = dble(cpatch%wood_gbw(ico))
               initp%wood_gbh (ico) = dble(cpatch%wood_gbh(ico))
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do

         !---------------------------------------------------------------------------------!
         !     Find the net ground conductance.  The net conductance is derived from the   !
         ! net resistance, which is, in turn, the weighted average of the resistances in   !
         ! bare and vegetated grounds.                                                     !
         !---------------------------------------------------------------------------------!
         if (initp%opencan_frac > 9.99d-1 .or. snowfac_can >= 9.d-1) then
            initp%ggnet = initp%ggbare
         else
            initp%ggnet = initp%ggbare * initp%ggveg                                       &
                        / ( initp%ggveg + (1.d0 - initp%opencan_frac) * initp%ggbare )
         end if
         !---------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      ! ED-2.1 Case: There uses a slightly different approach for wind profile and ground  !
      !              conductance.                                                          !
      !------------------------------------------------------------------------------------!
      case (1)
         !---------------------------------------------------------------------------------!
         !     The roughness is found by combining two weighted averages.  The first one   !
         ! checks the fraction of the patch that has closed canopy, and averages between   !
         ! soil and vegetation roughness.  The other is the fraction of the vegetation     !
         ! that is covered in snow.                                                        !
         !---------------------------------------------------------------------------------!
         initp%rough = snow_rough8 * snowfac_can                                           &
                     + ( soil_rough8     * initp%opencan_frac                              &
                       + initp%veg_rough * (1.d0 - initp%opencan_frac) )                   &
                     * (1.d0 - snowfac_can)
         !---------------------------------------------------------------------------------!

         !----- Calculate the soil surface roughness inside the canopy. -------------------!
         surf_rough = soil_rough8 * (1.d0 - snowfac_can) + snow_rough8 * snowfac_can
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Get ustar for the ABL, assume it is a dynamic shear layer that generates a !
         ! logarithmic profile of velocity.                                                !
         !---------------------------------------------------------------------------------!
         call ed_stars8(rk4site%atm_theta,initp%atm_enthalpy,rk4site%atm_shv               &
                       ,rk4site%atm_co2,initp%can_theta ,initp%can_enthalpy,initp%can_shv  &
                       ,initp%can_co2,initp%can_temp,rk4site%geoht,initp%veg_displace      &
                       ,rk4site%atm_ustar,initp%vels,initp%rough,initp%ustar,initp%tstar   &
                       ,initp%estar,initp%qstar,initp%cstar,initp%zeta,initp%ribulk        &
                       ,initp%ggbare)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     The surface conductance in the sub-canopy layer is the integration of       !
         ! K, from the rough soil surface, to the reference point in the canopy            !
         ! where the state variable is integrated (canopy top ~ h).                        !
         !---------------------------------------------------------------------------------!
         K_top = vonk8 * initp%ustar * (initp%veg_height - initp%veg_displace)
         initp%ggveg = (exar8 * K_top)                                                     &
                     / ( exp(exar8) * initp%veg_height                                     &
                       * (exp(-exar8 * surf_rough/initp%veg_height) - exp(-exar8)))
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     This part of the code finds the geometry of the canopy air space, the       !
         ! structure of the vegetation and its attenuation effects and the heat and water  !
         ! capacities.                                                                     !
         !---------------------------------------------------------------------------------!
         !----- Top of canopy wind speed. -------------------------------------------------!
         uh   = reduced_wind8(initp%ustar,initp%zeta,initp%ribulk,rk4site%geoht            &
                             ,initp%veg_displace,dble(cpatch%hite(1)),initp%rough)
         htop = dble(cpatch%hite(1))
         do ico=1,cpatch%ncohorts

            !----- Alias for PFT type. ----------------------------------------------------!
            ipft  = cpatch%pft(ico)
            !------------------------------------------------------------------------------!



            !----- Estimate the height center of the crown. -------------------------------!
            htopcrown = dble(cpatch%hite(ico))
            hbotcrown = dble(h2crownbh(cpatch%hite(ico),ipft))
            hmidcrown = 5.d-1 * (htopcrown + hbotcrown)
            !------------------------------------------------------------------------------!



            !----- Calculate the wind speed at height z. ----------------------------------!
            initp%veg_wind(ico) = max(ugbmin8                                              &
                                     , uh * exp(-exar8 * (1.d0 - hmidcrown/ htop)))

            !------ Leaf boundary layer conductance. --------------------------------------!
            if (initp%leaf_resolvable(ico)) then
               !---------------------------------------------------------------------------!
               !    Find the aerodynamic conductances for heat and water at the leaf       !
               ! boundary layer.                                                           !
               !---------------------------------------------------------------------------!
               call leaf_aerodynamic_conductances8(ipft,initp%veg_wind(ico)                &
                                                  ,initp%leaf_temp(ico),initp%can_temp     &
                                                  ,initp%can_rhos                          &
                                                  ,initp%can_cp                            &
                                                  ,initp%leaf_gbh(ico),initp%leaf_gbw(ico) &
                                                  ,initp%leaf_reynolds(ico)                &
                                                  ,initp%leaf_grashof(ico)                 &
                                                  ,initp%leaf_nussfree(ico)                &
                                                  ,initp%leaf_nussforc(ico) )
               cpatch%leaf_gbh(ico) = sngloff(initp%leaf_gbh(ico),tiny_offset)
               cpatch%leaf_gbw(ico) = sngloff(initp%leaf_gbw(ico),tiny_offset)
               !---------------------------------------------------------------------------!
            else
               !----- Set numbers to zero. ------------------------------------------------!
               initp%leaf_reynolds(ico) = 0.d0
               initp%leaf_grashof (ico) = 0.d0
               initp%leaf_nussfree(ico) = 0.d0
               initp%leaf_nussforc(ico) = 0.d0
               !---------------------------------------------------------------------------!


               !----- Leaf conductances cannot be zero.  Set to non-zero defaults. --------!
               cpatch%leaf_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
               cpatch%leaf_gbh(ico) = cpatch%leaf_gbw(ico) / gbh_2_gbw * cpdry
               initp%leaf_gbw (ico) = dble(cpatch%leaf_gbw(ico))
               initp%leaf_gbh (ico) = dble(cpatch%leaf_gbh(ico))
               !---------------------------------------------------------------------------!
            end if
            !------ Wood boundary layer conductance. --------------------------------------!
            if (initp%wood_resolvable(ico)) then
               !---------------------------------------------------------------------------!
               !    Find the aerodynamic conductances for heat and water at the leaf       !
               ! boundary layer.                                                           !
               !---------------------------------------------------------------------------!
               call wood_aerodynamic_conductances8(ipft                                    &
                                                  ,initp%veg_wind(ico)                     &
                                                  ,initp%wood_temp(ico),initp%can_temp     &
                                                  ,initp%can_rhos                          &
                                                  ,initp%can_cp                            &
                                                  ,initp%wood_gbh(ico),initp%wood_gbw(ico) &
                                                  ,initp%wood_reynolds(ico)                &
                                                  ,initp%wood_grashof(ico)                 &
                                                  ,initp%wood_nussfree(ico)                &
                                                  ,initp%wood_nussforc(ico) )
               cpatch%wood_gbh(ico) = sngloff(initp%wood_gbh(ico),tiny_offset)
               cpatch%wood_gbw(ico) = sngloff(initp%wood_gbw(ico),tiny_offset)
               !---------------------------------------------------------------------------!
            else
               !----- Set numbers to zero. ------------------------------------------------!
               initp%wood_reynolds(ico) = 0.d0
               initp%wood_grashof (ico) = 0.d0
               initp%wood_nussfree(ico) = 0.d0
               initp%wood_nussforc(ico) = 0.d0
               !---------------------------------------------------------------------------!


               !----- Leaf conductances cannot be zero.  Set to non-zero defaults. --------!
               cpatch%wood_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
               cpatch%wood_gbh(ico) = cpatch%wood_gbw(ico) / gbh_2_gbw * cpdry
               initp%wood_gbw (ico) = dble(cpatch%wood_gbw(ico))
               initp%wood_gbh (ico) = dble(cpatch%wood_gbh(ico))
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the net ground conductance.  The net conductance is derived from the   !
         ! net resistance, which is, in turn, the weighted average of the resistances in   !
         ! bare and vegetated grounds.                                                     !
         !---------------------------------------------------------------------------------!
         if (initp%opencan_frac > 9.99d-1 .or. snowfac_can >= 9.d-1) then
            initp%ggnet = initp%ggbare
         else
            initp%ggnet = initp%ggbare * initp%ggveg                                       &
                        / (initp%ggveg + (1.d0 - initp%opencan_frac) * initp%ggbare )
         end if
         !---------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Use the methods of Massman (1997) or Massman and Weil (1999).  Using option 2 !
      ! will make the within-canopy drag and sheltering factors to be constant (default    !
      ! Massman 1997), whereas option 3 will allow the values to vary and use the second   !
      ! order closure as in Massman and Weil (1999).                                       !
      !------------------------------------------------------------------------------------!
      case (2:3)
         !---------------------------------------------------------------------------------!
         !    Find the top layer and the top height.                                       !
         !---------------------------------------------------------------------------------!
         zcan = min(ncanlyr,ceiling((dble(cpatch%hite(1)) * zztop0i8)**ehgti8))
         htop = zztop8(zcan)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Loop through cohorts, and integrate the leaf area density.  Notice that we  !
         ! will be solving all cohorts here, because even those with no leaves have        !
         ! effects on the canopy air turbulence.  We only need to decide whether we have   !
         ! branches or not, and we do it outside the loop to make it more efficient.       !
         !---------------------------------------------------------------------------------!
         !----- Reset the leaf area density array. ----------------------------------------!
         canstr(ibuff)%lad8(:) = 0.d0
         select case (ibranch_thermo)
         case (0)
            !------------------------------------------------------------------------------!
            !     Although we are not solving branch thermodynamics, we have some drag     !
            ! effect due to the branches.  Assume branch area index to be 15% of on-       !
            ! -allometry LAI.  Even for grasses we may want to keep some residual due to   !
            ! dead material that remains around.                                           !
            !------------------------------------------------------------------------------!
            do ico=1,cpatch%ncohorts
               ipft = cpatch%pft(ico)


               !------ Estimate WAI. ------------------------------------------------------!
               if (is_grass(ipft) .and. igrass==1) then
                   !--use actual leaf mass for grass
                   waiuse = 1.d-1 * initp%nplant(ico) * dble(cpatch%sla(ico))              &
                          * dble(cpatch%bleaf(ico))
               else
                   !--use dbh for trees
                   waiuse = 1.d-1 * initp%nplant(ico) * dble(cpatch%sla(ico))              &
                          * dble(size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft))
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the heights, and compute the LAD of this cohort.                 !
               !---------------------------------------------------------------------------!
               htopcrown = dble(cpatch%hite(ico))
               hbotcrown = dble(h2crownbh(cpatch%hite(ico),cpatch%pft(ico)))
               ladcohort = (initp%lai(ico) + waiuse) / (htopcrown - hbotcrown)
               kapartial = min(ncanlyr,floor  ((hbotcrown * zztop0i8)**ehgti8) + 1)
               kafull    = min(ncanlyr,ceiling((hbotcrown * zztop0i8)**ehgti8) + 1)
               kzpartial = min(ncanlyr,ceiling((htopcrown * zztop0i8)**ehgti8))
               kzfull    = min(ncanlyr,floor  ((htopcrown * zztop0i8)**ehgti8))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Add the LAD for the full layers.                                      !
               !---------------------------------------------------------------------------!
               do k = kafull,kzfull
                  canstr(ibuff)%lad8(k) = canstr(ibuff)%lad8(k) + ladcohort
               end do
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Add the LAD for the partial layers.  The only special case is when   !
               ! they are both the same layer, which must be done separately.              !
               !---------------------------------------------------------------------------!
               if (kapartial == kzpartial) then
                  canstr(ibuff)%lad8(kapartial) = canstr(ibuff)%lad8(kapartial)            &
                                  + ladcohort * (htopcrown         - hbotcrown        )    &
                                              / (zztop8(kapartial) - zzbot8(kapartial))
               else
                  canstr(ibuff)%lad8(kapartial) = canstr(ibuff)%lad8(kapartial)            &
                                  + ladcohort * (zztop8(kapartial) - hbotcrown        )    &
                                              / (zztop8(kapartial) - zzbot8(kapartial))
                  canstr(ibuff)%lad8(kzpartial) = canstr(ibuff)%lad8(kzpartial)            &
                                  + ladcohort * (htopcrown         - zzbot8(kzpartial))    &
                                              / (zztop8(kzpartial) - zzbot8(kzpartial))
               end if
               !---------------------------------------------------------------------------!
            end do

         case default

            !------------------------------------------------------------------------------!
            !     Branches are turned on.  In arid places, there is a chance that all      !
            ! cohorts will be grasses with phenology status is set to 2.  This creates a   !
            ! singularity, so we must check whether this is the case.                      !
            !------------------------------------------------------------------------------!
            dry_grasses = sum(cpatch%lai(:)+cpatch%wai(:)) == 0.0
            !------------------------------------------------------------------------------!


            !----- Use the default wood area index. ---------------------------------------!
            do ico=1,cpatch%ncohorts
               ipft = cpatch%pft(ico)

               !---------------------------------------------------------------------------!
               !     Find the heights, and compute the LAD of this cohort.                 !
               !---------------------------------------------------------------------------!
               htopcrown = dble(cpatch%hite(ico))
               hbotcrown = dble(h2crownbh(cpatch%hite(ico),ipft))
               if (dry_grasses) then
                  !------------------------------------------------------------------------!
                  !     Dry grasses only.  Create a pseudo TAI so it won't be a            !
                  ! singularity.                                                           !
                  !------------------------------------------------------------------------!
                  tai_drygrass = dble( elongf_min                                          &
                                     * size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft))
                  ladcohort    = tai_drygrass / (htopcrown - hbotcrown)
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !     At least one plant has branches or leaves, use the real stuff      !
                  ! instead.                                                               !
                  !------------------------------------------------------------------------!
                  ladcohort = (initp%lai(ico) + initp%wai(ico)) / (htopcrown - hbotcrown)
                  !------------------------------------------------------------------------!
               end if
               kapartial = min(ncanlyr,floor  ((hbotcrown * zztop0i8)**ehgti8) + 1)
               kafull    = min(ncanlyr,ceiling((hbotcrown * zztop0i8)**ehgti8) + 1)
               kzpartial = min(ncanlyr,ceiling((htopcrown * zztop0i8)**ehgti8))
               kzfull    = min(ncanlyr,floor  ((htopcrown * zztop0i8)**ehgti8))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Add the LAD for the full layers.                                      !
               !---------------------------------------------------------------------------!
               do k = kafull,kzfull
                  canstr(ibuff)%lad8(k) = canstr(ibuff)%lad8(k) + ladcohort
               end do
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Add the LAD for the partial layers.  The only special case is when   !
               ! they are both the same layer, which must be done separately.              !
               !---------------------------------------------------------------------------!
               if (kapartial == kzpartial) then
                  canstr(ibuff)%lad8(kapartial) = canstr(ibuff)%lad8(kapartial)            &
                                  + ladcohort * (htopcrown         - hbotcrown        )    &
                                              / (zztop8(kapartial) - zzbot8(kapartial))
               else
                  canstr(ibuff)%lad8(kapartial) = canstr(ibuff)%lad8(kapartial)            &
                                  + ladcohort * (zztop8(kapartial) - hbotcrown        )    &
                                              / (zztop8(kapartial) - zzbot8(kapartial))
                  canstr(ibuff)%lad8(kzpartial) = canstr(ibuff)%lad8(kzpartial)            &
                                  + ladcohort * (htopcrown         - zzbot8(kzpartial))    &
                                              / (zztop8(kzpartial) - zzbot8(kzpartial))
               end if
               !---------------------------------------------------------------------------!
            end do
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     In this loop we find the drag coefficient, the sheltering factor, and the   !
         ! cumulative leaf drag area.  This must be done in a separate loop because we     !
         ! need the full integral of the leaf area density before we determine these       !
         ! variables.                                                                      !
         !---------------------------------------------------------------------------------!
         select case (icanturb)
         case (2)
            !----- Drag when there are no plants. -----------------------------------------!
            canstr(ibuff)%cdrag8  (:)  = cdrag18 + 5.d-1 * cdrag28
            ldga_bk      = 0.d0
            do k = 1,zcan
               !---------------------------------------------------------------------------!
               !     Add the contribution of this layer to Massman's zeta (which we call   !
               ! cumldrag here to not confuse with the other zeta from the similarity      !
               ! theory.  We integrate in three steps so we save the value in the middle   !
               ! of the layer.                                                             !
               !                                                                           !
               !     We use a re-fit of the original equation by Wohlfahrt and Cernusca    !
               ! (2002) because it is simpler and because the coefficients listed in the   !
               ! paper do not yield to the results shown in their fig. 3 (probably a typo) !
               ! Their cdeff is cdrag / pshelter, here we fix pshelter = 1 and dump the    !
               ! ratio to cdrag.                                                           !
               !---------------------------------------------------------------------------!
               c3_lad       = &
                     max(lnexp_min8,min(lnexp_max8,cdrag38 * canstr(ibuff)%lad8(k)))
               canstr(ibuff)%cdrag8   (k) = cdrag18 + cdrag28 / (1.d0 + exp(c3_lad))
               canstr(ibuff)%pshelter8(k) = 1.d0
               lyrhalf      = 5.d-1 * canstr(ibuff)%lad8(k) * canstr(ibuff)%cdrag8(k) /    &
                     canstr(ibuff)%pshelter8(k) * dzcan8(k)
               canstr(ibuff)%cumldrag8(k) = ldga_bk + lyrhalf
               ldga_bk      = ldga_bk + 2.d0 * lyrhalf
               !---------------------------------------------------------------------------!
            end do
         case (3)
            !----- Constant drag. ---------------------------------------------------------!
            canstr(ibuff)%cdrag8   (:) = cdrag08
            ldga_bk      = 0.d0
            do k = 1,zcan
               !---------------------------------------------------------------------------!
               !     Add the contribution of this layer to Massman's zeta (which we call   !
               ! cumldrag here to not confuse with the other zeta from the similarity      !
               ! theory).  We integrate in three steps so we save the value in the middle  !
               ! of the layer.                                                             !
               !     Notice that pshelter is the inverse of the parametrisation shown in   !
               ! Massman (1997) (a typo according to personal communication between Ryan   !
               ! and Massman), since the shelter factor should be always >= 1.             !
               ! Alpha_m97 is no longer 0.4 * h, but 5 so it becomes a constant.           !
               !---------------------------------------------------------------------------!
               canstr(ibuff)%pshelter8(k) = 1.d0 + alpha_m97_8 * canstr(ibuff)%lad8(k)
               lyrhalf      = 5.d-1 * canstr(ibuff)%lad8(k) * canstr(ibuff)%cdrag8(k) /    &
                     canstr(ibuff)%pshelter8(k) * dzcan8(k)
               canstr(ibuff)%cumldrag8(k) = ldga_bk + lyrhalf
               ldga_bk      = ldga_bk + 2.d0 * lyrhalf
               !---------------------------------------------------------------------------!
            end do
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Find the ratio between u* and u at the top cohort, using Massman's equation  !
         ! (6).                                                                            !
         !---------------------------------------------------------------------------------!
         c3_cumldrag = min(lnexp_max8                                                      &
                          ,max(lnexp_min8,c3_m978 * canstr(ibuff)%cumldrag8(zcan)))
         ustarouh    = (c1_m978 - c2_m978 * exp(- c3_cumldrag))
         !---------------------------------------------------------------------------------!



         !----- NN is Massman's n, the coefficient of attenuation. ------------------------!
         nn = 5.d-1 * canstr(ibuff)%cumldrag8(zcan) / (ustarouh * ustarouh)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the ratio between zero-plane displacement height and height, and the   !
         ! ratio between roughness and height.  The actual values will be defined based    !
         ! on the patch-level averaged height, not the height of the tallest cohort,       !
         ! because it may be an extremely sparse cohort that has little impact on the      !
         ! drag.                                                                           !
         !---------------------------------------------------------------------------------!
         d0ohgt = 1.d0
         do k=1,zcan
            d0ohgt = d0ohgt - dzcan8(k) / htop                                             &
                            * exp(-2.d0 * nn *(1.d0 - canstr(ibuff)%cumldrag8(k) /         &
                            canstr(ibuff)%cumldrag8(zcan)))
         end do
         z0ohgt = (1.d0 - d0ohgt) * min(1.d0,exp(- vonk8 / ustarouh + infunc_8))
         !---------------------------------------------------------------------------------!




         !----- Find the actual displacement height and roughness. ------------------------!
         initp%veg_displace = max( vh2dh8 * veg_height_min8, d0ohgt * initp%veg_height)
         initp%rough        = max( vh2vr8 * veg_height_min8, z0ohgt * initp%veg_height)
         !---------------------------------------------------------------------------------!


         !----- Calculate ustar, tstar, qstar, and cstar. ---------------------------------!
         call ed_stars8(rk4site%atm_theta,initp%atm_enthalpy,rk4site%atm_shv               &
                       ,rk4site%atm_co2,initp%can_theta ,initp%can_enthalpy,initp%can_shv  &
                       ,initp%can_co2,initp%can_temp,rk4site%geoht,initp%veg_displace      &
                       ,rk4site%atm_ustar,initp%vels,initp%rough,initp%ustar,initp%tstar   &
                       ,initp%estar,initp%qstar,initp%cstar,initp%zeta,initp%ribulk        &
                       ,initp%ggbare)
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Find the wind profile.                                                      !
         !---------------------------------------------------------------------------------!
         !----- Top of canopy wind speed. -------------------------------------------------!
         uh = reduced_wind8(initp%ustar,initp%zeta,initp%ribulk,rk4site%geoht              &
                           ,initp%veg_displace,htop,initp%rough)
         !----- Get the wind profile. -----------------------------------------------------!
         do k=1,zcan
            !----- Normalised drag density fraction and wind for this layer. --------------!
            nddfun      = 1.d0 - canstr(ibuff)%cumldrag8(k) / canstr(ibuff)%cumldrag8(zcan)
            canstr(ibuff)%windlyr8(k) = max(ugbmin8, uh * exp(- nn * nddfun))
         end do
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Calculate the leaf level aerodynamic resistance.                            !
         !---------------------------------------------------------------------------------!
         do ico=1,cpatch%ncohorts
            ipft = cpatch%pft(ico)

            !----- Find the crown relevant heights. ---------------------------------------!
            htopcrown = dble(cpatch%hite(ico))
            hbotcrown = dble(h2crownbh(cpatch%hite(ico),cpatch%pft(ico)))
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find the layer indices for the bottom and top of the crown.              !
            !------------------------------------------------------------------------------!
            kapartial = min(ncanlyr,floor  ((hbotcrown * zztop0i8)**ehgti8) + 1)
            kafull    = min(ncanlyr,ceiling((hbotcrown * zztop0i8)**ehgti8) + 1)
            kzpartial = min(ncanlyr,ceiling((htopcrown * zztop0i8)**ehgti8))
            kzfull    = min(ncanlyr,floor  ((htopcrown * zztop0i8)**ehgti8))
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Add the LAD for the full layers.                                         !
            !------------------------------------------------------------------------------!
            if ( kapartial == kzpartial ) then
               !----- Cohort crown is in a single layer, copy the layer wind speed. -------!
               initp%veg_wind(ico) = canstr(ibuff)%windlyr8(kapartial)
            else
               !---------------------------------------------------------------------------!
               !      Cohort spans through multiple layers.  Use the average, weighted by  !
               ! the thickness of the layer.                                               !
               !---------------------------------------------------------------------------!
               !----- Partial layers (bottom and top). ------------------------------------!
               initp%veg_wind(ico) = &
                     canstr(ibuff)%windlyr8(kapartial) * (zztop8(kapartial) - hbotcrown)   &
                   + canstr(ibuff)%windlyr8(kzpartial) * (htopcrown - zzbot8(kzpartial))
               do k = kafull,kzfull
                  initp%veg_wind(ico) = initp%veg_wind(ico) +                              &
                        canstr(ibuff)%windlyr8(k) * dzcan8(k)
               end do
               !----- Divide by the total crown length to obtain the average wind. --------!
               initp%veg_wind(ico) = initp%veg_wind(ico) / (htopcrown - hbotcrown)
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !------ Leaf boundary layer conductance. --------------------------------------!
            if (initp%leaf_resolvable(ico)) then
               !---------------------------------------------------------------------------!
               !    Find the aerodynamic conductances for heat and water at the leaf       !
               ! boundary layer.                                                           !
               !---------------------------------------------------------------------------!
               call leaf_aerodynamic_conductances8(ipft,initp%veg_wind(ico)                &
                                                  ,initp%leaf_temp(ico),initp%can_temp     &
                                                  ,initp%can_rhos                          &
                                                  ,initp%can_cp                            &
                                                  ,initp%leaf_gbh(ico),initp%leaf_gbw(ico) &
                                                  ,initp%leaf_reynolds(ico)                &
                                                  ,initp%leaf_grashof(ico)                 &
                                                  ,initp%leaf_nussfree(ico)                &
                                                  ,initp%leaf_nussforc(ico) )
               cpatch%leaf_gbh(ico) = sngloff(initp%leaf_gbh(ico),tiny_offset)
               cpatch%leaf_gbw(ico) = sngloff(initp%leaf_gbw(ico),tiny_offset)
               !---------------------------------------------------------------------------!
            else
               !----- Set numbers to zero. ------------------------------------------------!
               initp%leaf_reynolds(ico) = 0.d0
               initp%leaf_grashof (ico) = 0.d0
               initp%leaf_nussfree(ico) = 0.d0
               initp%leaf_nussforc(ico) = 0.d0
               !---------------------------------------------------------------------------!


               !----- Leaf conductances cannot be zero.  Set to non-zero defaults. --------!
               cpatch%leaf_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
               cpatch%leaf_gbh(ico) = cpatch%leaf_gbw(ico) / gbh_2_gbw * cpdry
               initp%leaf_gbw (ico) = dble(cpatch%leaf_gbw(ico))
               initp%leaf_gbh (ico) = dble(cpatch%leaf_gbh(ico))
               !---------------------------------------------------------------------------!
            end if
            !------ Wood boundary layer conductance. --------------------------------------!
            if (initp%wood_resolvable(ico)) then
               !---------------------------------------------------------------------------!
               !    Find the aerodynamic conductances for heat and water at the leaf       !
               ! boundary layer.                                                           !
               !---------------------------------------------------------------------------!
               call wood_aerodynamic_conductances8(ipft                                    &
                                                  ,initp%veg_wind(ico)                     &
                                                  ,initp%wood_temp(ico),initp%can_temp     &
                                                  ,initp%can_rhos                          &
                                                  ,initp%can_cp                            &
                                                  ,initp%wood_gbh(ico),initp%wood_gbw(ico) &
                                                  ,initp%wood_reynolds(ico)                &
                                                  ,initp%wood_grashof(ico)                 &
                                                  ,initp%wood_nussfree(ico)                &
                                                  ,initp%wood_nussforc(ico) )
               cpatch%wood_gbh(ico) = sngloff(initp%wood_gbh(ico),tiny_offset)
               cpatch%wood_gbw(ico) = sngloff(initp%wood_gbw(ico),tiny_offset)
               !---------------------------------------------------------------------------!
            else
               !----- Set numbers to zero. ------------------------------------------------!
               initp%wood_reynolds(ico) = 0.d0
               initp%wood_grashof (ico) = 0.d0
               initp%wood_nussfree(ico) = 0.d0
               initp%wood_nussforc(ico) = 0.d0
               !---------------------------------------------------------------------------!


               !----- Leaf conductances cannot be zero.  Set to non-zero defaults. --------!
               cpatch%wood_gbw(ico) = f_bndlyr_init * cpatch%leaf_gsw(ico)
               cpatch%wood_gbh(ico) = cpatch%wood_gbw(ico) / gbh_2_gbw * cpdry
               initp%wood_gbw (ico) = dble(cpatch%wood_gbw(ico))
               initp%wood_gbh (ico) = dble(cpatch%wood_gbh(ico))
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !     Decide how to solve the aerodynamic resistance based on the user's option.  !
         !---------------------------------------------------------------------------------!
         select case (icanturb)
         case (2)
            !------------------------------------------------------------------------------!
            !     Find the bulk resistance by integrating the inverse of the diffusivity   !
            ! for each layer (e.g. Sellers et al. (1986)).  We assumed Km = sigma u, which !
            ! is the preferred approach by Sellers et al. (1986), and divide by the        !
            ! Prandtl number to convert Km to Kh.                                          !
            !------------------------------------------------------------------------------!
            rasveg  = 0.d0
            zetatop = initp%zeta * (htop          - initp%veg_displace)                    &
                                 / (rk4site%geoht - initp%veg_displace)
            sigmakh = vonk8 * initp%ustar * htop * (1.d0 - d0ohgt)                         &
                    / ( tprandtl8 * uh * ( 1.d0 - zetatop * dpsimdzeta8(zetatop,stable) ) )
            do k=1,zcan
               !---------------------------------------------------------------------------!
               !    Find the normalised drag density fraction and wind for this layer.     !
               !---------------------------------------------------------------------------!
               Kdiff      = sigmakh * canstr(ibuff)%windlyr8(k) + kvwake8
               rasveg     = rasveg + dzcan8(k) / Kdiff
            end do
            initp%ggveg = 1.d0 / rasveg
            !------------------------------------------------------------------------------!
         case (3)
            !------------------------------------------------------------------------------!
            !     Use MW99 second order analytical solution to sigma_w and turb intensity. !
            ! To calculate the Reynolds number, we will use the mean velocity over the     !
            ! depth of the mixing length scale (which is estimated as the roughness).      !
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     Eddy length scale at the soil surface (estimates other than roughness    !
            ! are up for debate, but be warned this length scale is rather stable.         !
            !------------------------------------------------------------------------------!
            elenscale = initp%rough
            zels      = min(ncanlyr,ceiling((elenscale * zztop0i8)**ehgti8))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Fill in with the wind using the similarity theory in case zels is       !
            ! greater than zcan.                                                           !
            !------------------------------------------------------------------------------!
            do k=zcan+1,zels
               canstr(ibuff)%windlyr8(k) =  reduced_wind8(initp%ustar,initp%zeta           &
                                                         ,initp%ribulk,rk4site%geoht       &
                                                         ,initp%veg_displace,zzmid8(k)     &
                                                         ,initp%rough)
            end do
            !------------------------------------------------------------------------------!


            !----- Initialise alpha, failure counters, and the logical flag. --------------!
            alpha_eq10   = alpha_mw99_8
            nalpha_fails = 0
            acomp        = .true.
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !   Go through this loop until we find a safe solution or if we give up and    !
            ! fall back to the first-order.                                                !
            !------------------------------------------------------------------------------!
            afail: do

               lam = srthree8 * nu_mw99_8(1) / alpha_eq10

               !----- Arbitrary coefficient in analytical solution. -----------------------!
               b1_mw99 = - (9.d0 * ustarouh)                                               &
                       / ( 2.d0 * alpha_eq10 * nu_mw99_8(1)                                &
                         * (2.25d0 - lam * lam * (initp%ustar/uh) ** 4))

               ure   = 0.d0
               turbi = 0.d0
               do k=1,zels
                  !------------------------------------------------------------------------!
                  !    Find the normalised drag density fraction and wind for this layer.  !
                  !------------------------------------------------------------------------!
                  nddfun      = 1.d0 - canstr(ibuff)%cumldrag8(k) /                        &
                        canstr(ibuff)%cumldrag8(zcan)

                  !------------------------------------------------------------------------!
                  !    Integrate the wind speed.  It will be normalised outside the loop.  !
                  !------------------------------------------------------------------------!
                  ure        = ure + canstr(ibuff)%windlyr8(k) * dzcan8(k)

                  !----- Sigstar, as in equation 10 of MW99. ------------------------------!
                  sigstar3   = &
                        nu_mw99_8(3) * exp( - lam * canstr(ibuff)%cumldrag8(zcan) * nddfun)&
                        + b1_mw99 * ( exp( - 3.d0 * nn * nddfun)                      &
                        - exp( - lam * canstr(ibuff)%cumldrag8(zcan) * nddfun))
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      It is possible (but highly unlikely) that sigstar3 could be       !
                  ! negative.  Even though cubic roots of negative number are fine, it     !
                  ! wouldn't make sense to have sigma_u/sigma_v/sigma_w to be negative.    !
                  !------------------------------------------------------------------------!
                  if (sigstar3 < 0.d0) then
                     nalpha_fails = nalpha_fails + 1
                     alpha_eq10   = alpha_eq10 - 0.001
                     if (nalpha_fails == 10) then
                        !------------------------------------------------------------------!
                        !     Find the bulk resistance by integrating the inverse of the   !
                        ! diffusivity for each layer (e.g. Sellers et al. (1986)).  We     !
                        ! assumed Km = sigma u, which is the preferred approach by Sellers !
                        ! et al. (1986), and divide by the Prandtl number to convert Km to !
                        ! Kh.                                                              !
                        !------------------------------------------------------------------!
                        rasveg  = 0.d0
                        zetatop = initp%zeta * (htop          - initp%veg_displace)        &
                                             / (rk4site%geoht - initp%veg_displace)
                        sigmakh = vonk8 * initp%ustar * htop * (1.d0 - d0ohgt)             &
                                / ( tprandtl8 * uh                                         &
                                  * ( 1.d0 - zetatop * dpsimdzeta8(zetatop,stable) ) )
                        do kk=1,zcan
                           Kdiff       = sigmakh * canstr(ibuff)%windlyr8(kk) + kvwake8
                           rasveg      = rasveg + dzcan8(kk) / Kdiff
                        end do
                        initp%ggveg    = 1.d0 / rasveg
                        exit afail
                     end if

                     cycle afail
                     !---------------------------------------------------------------------!
                  else
                     sigstar  = cbrt8(sigstar3)
                     sigcomm  = initp%ustar * sigstar * nu_mw99_8(1)
                     sigma_uou2 = (sigcomm * gamma_mw99_8(1) /                             &
                           canstr(ibuff)%windlyr8(k)) ** 2
                     sigma_vou2 = (sigcomm * gamma_mw99_8(2) /                             &
                           canstr(ibuff)%windlyr8(k)) ** 2
                     sigma_wou2 = (sigcomm * gamma_mw99_8(3) /                             &
                           canstr(ibuff)%windlyr8(k)) ** 2
                     turbi    = turbi                                                      &
                              + sqrt(onethird8 * (sigma_uou2 + sigma_vou2 + sigma_wou2))   &
                              * dzcan8(k)
                     exit afail
                  end if
               end do
            end do afail
            ure   = ure   / zztop8(zels)
            turbi = turbi / zztop8(zels)

            !------ Normalise the Reynolds number by diffusivity. -------------------------!
            can_reynolds = ure / ( th_diff08 * (1.d0 + dth_diff8 * (initp%can_temp-t008)))


            !------------------------------------------------------------------------------!
            !     Find the aerodynamic conductance based on MW99.                          !
            !------------------------------------------------------------------------------!
            initp%ggveg = sqrt(elenscale)                                                  &
                        * (1.d0 + 2.d0 * turbi) * tprandtl8 ** (-twothirds8)               &
                        / sqrt(can_reynolds)                                               &
                        * uh * sqrt(ure / uh)
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !     Find the net ground conductance.  The net conductance is derived from the   !
         ! net resistance, which is, in turn, the weighted average of the resistances in   !
         ! bare and vegetated grounds.                                                     !
         !---------------------------------------------------------------------------------!
         if (initp%opencan_frac > 9.99d-1 .or. snowfac_can >= 9.d-1) then
            initp%ggnet = initp%ggbare
         else
            initp%ggnet = initp%ggbare * initp%ggveg                                       &
                        / (initp%ggveg + (1.d0 - initp%opencan_frac) * initp%ggbare )
         end if
         !---------------------------------------------------------------------------------!

      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine canopy_turbulence8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Double precision version of the above-described subroutine.  This subroutine       !
   ! computes the characteristic scales based on surface layer parameterization.  Assume   !
   ! that stability is calculated on the potential gradient from the surface to the        !
   ! atmosphere at reference height.  Apply the instability parameters to calculate the    !
   ! friction velocity.  Use the instability parameters for momentum and scalars, that     !
   ! were calculated over the distance from surface to reference height to determine the   !
   ! heat, moisture and carbon flux rates at the canopy to atmosphere at reference height. !
   !    Two models are available, and the user can choose between them by setting the      !
   ! variable ISFCLYRM in ED2IN (if running the coupled model, this is done in ISTAR).     !
   !                                                                                       !
   ! 0. Same as option 4 (see below), but u* is prescribed.  Not intended for long-term    !
   !    simulations (where vegetation dynamics must affect u*), this should be used for    !
   !    testing only.                                                                      !
   ! 1. Based on L79;                                                                      !
   ! 2. Based on: OD95, but with some terms computed as in L79 and B71 to avoid singular-  !
   !    ities (now using the iterative method to find zeta).                               !
   ! 3. Based on BH91, using an iterative method to find zeta, and using the modified      !
   !    equation for stable layers.                                                        !
   ! 4. Based on CLM04, with special functions for very stable and very stable case, even  !
   !    though we use a different functional form for very unstable case for momentum.     !
   !    This is ensure that phi_m decreases monotonically as zeta becomes more negative.   !
   !    We use a power law of order of -1/6 instead.                                       !
   !                                                                                       !
   ! References:                                                                           !
   ! B71.  BUSINGER, J.A, et. al; Flux-Profile relationships in the atmospheric surface    !
   !           layer. J. Atmos. Sci., 28, 181-189, 1971.                                   !
   ! L79.  LOUIS, J.F.; Parametric Model of vertical eddy fluxes in the atmosphere.        !
   !           Boundary-Layer Meteor., 17, 187-202, 1979.                                  !
   ! BH91. BELJAARS, A.C.M.; HOLTSLAG, A.A.M.; Flux parameterization over land surfaces    !
   !           for atmospheric models. J. Appl. Meteor., 30, 327-341, 1991.                !
   ! OD95. ONCLEY, S.P.; DUDHIA, J.; Evaluation of surface fluxes from MM5 using observa-  !
   !           tions.  Mon. Wea. Rev., 123, 3344-3357, 1995.                               !
   ! CLM04. OLESON, K. W., et al.; Technical description of the community land model (CLM) !
   !           NCAR Technical Note NCAR/TN-461+STR, Boulder, CO, May 2004.                 !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine ed_stars8(theta_atm,enthalpy_atm,shv_atm,co2_atm                             &
                       ,theta_can,enthalpy_can,shv_can,co2_can,temp_can                    &
                       ,zref,dheight,atm_ustar,uref,rough,ustar,tstar,estar,qstar,cstar    &
                       ,zeta,rib,ggbare)
      use consts_coms     , only : grav8            & ! intent(in)
                                 , vonk8            & ! intent(in)
                                 , epim18           & ! intent(in)
                                 , halfpi8          & ! intent(in)
                                 , t008             & ! intent(in)
                                 , kin_visc08       & ! intent(in)
                                 , dkin_visc8       ! ! intent(in)
      use canopy_air_coms , only : isfclyrm         & ! intent(in)
                                 , ustmin8          & ! intent(in)
                                 , bl798            & ! intent(in)
                                 , csm8             & ! intent(in)
                                 , csh8             & ! intent(in)
                                 , dl798            & ! intent(in)
                                 , ribmax8          & ! intent(in)
                                 , tprandtl8        & ! intent(in)
                                 , zd98_a8          & ! intent(in)
                                 , zd98_b8          & ! intent(in)
                                 , zd98_emax8       & ! intent(in)
                                 , psim8            & ! function
                                 , psih8            & ! function
                                 , zoobukhov8       & ! function
                                 , zoobukhov_ustar8 ! ! function
      use rk4_coms        , only : rk4eps           ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)    :: theta_atm    ! Free air pot. temp.         [        K]
      real(kind=8), intent(in)    :: enthalpy_atm ! Free air spec. enthalpy     [ J/kg_air]
      real(kind=8), intent(in)    :: shv_atm      ! Free air vapour spec. hum.  [kg/kg_air]
      real(kind=8), intent(in)    :: co2_atm      ! Free air CO2 mixing ratio   [ umol/mol]
      real(kind=8), intent(in)    :: theta_can    ! CAS pot. temperature        [        K]
      real(kind=8), intent(in)    :: enthalpy_can ! CAS specific enthalpy       [ J/kg_air]
      real(kind=8), intent(in)    :: shv_can      ! CAS vapour spec. hum.       [kg/kg_air]
      real(kind=8), intent(in)    :: co2_can      ! CAS CO2 mixing ratio        [ umol/mol]
      real(kind=8), intent(in)    :: temp_can     ! CAS temperature             [        K]
      real(kind=8), intent(in)    :: zref         ! Height at reference point   [        m]
      real(kind=8), intent(in)    :: dheight      ! Zero-plane displace. hgt.   [        m]
      real(kind=8), intent(in)    :: atm_ustar    ! Prescribed u*               [      m/s]
      real(kind=8), intent(in)    :: uref         ! Wind speed at ref. hgt.     [      m/s]
      real(kind=8), intent(in)    :: rough        ! Roughness                   [        m]
      real(kind=8), intent(inout) :: ustar        ! U*, friction velocity       [      m/s]
      real(kind=8), intent(out)   :: qstar        ! Spec. humidity turb. scale  [kg/kg_air]
      real(kind=8), intent(out)   :: tstar        ! Temp. turbulence scale      [        K]
      real(kind=8), intent(out)   :: estar        ! Spec. enthalpy turb. scale  [ J/kg_air]
      real(kind=8), intent(out)   :: cstar        ! CO2 mix. ratio turb. scale  [ umol/mol]
      real(kind=8), intent(out)   :: zeta         ! z/(Obukhov length).         [    -----]
      real(kind=8), intent(out)   :: rib          ! Bulk richardson number.     [    -----]
      real(kind=8), intent(out)   :: ggbare       ! Ground conductance          [      m/s]
      !----- Local variables --------------------------------------------------------------!
      logical                     :: stable       ! Stable state                [      T|F]
      real(kind=8)                :: kin_visc     ! Kinematic viscosity         [     m2/s]
      real(kind=8)                :: z0hoz0m      ! z0h:z0m ratio               [    -----]
      real(kind=8)                :: z0moz0h      ! z0m:z0h ratio               [    -----]
      real(kind=8)                :: zoz0m        ! zref/rough(momentum)        [    -----]
      real(kind=8)                :: lnzoz0m      ! ln[zref/rough(momentum)]    [    -----]
      real(kind=8)                :: zoz0h        ! zref/rough(heat)            [    -----]
      real(kind=8)                :: lnzoz0h      ! ln[zref/rough(heat)]        [    -----]
      real(kind=8)                :: c3           ! aux. coefficient            [    -----]
      real(kind=8)                :: uuse         ! Wind when (Rib > Ribmax)    [      m/s]
      real(kind=8)                :: kuoustar     ! k * u / u*
      !----- Local variables, used by L79. ------------------------------------------------!
      real(kind=8)                :: a2           ! Drag coefficient in neutral conditions
      real(kind=8)                :: c1           ! a2 * vels
      real(kind=8)                :: fm           ! Stability parameter for momentum
      real(kind=8)                :: fh           ! Stability parameter for heat
      real(kind=8)                :: c2           ! Part of the c coefficient common
                                                  !      to momentum & heat.
      real(kind=8)                :: cm           ! c times |Rib|^1/2 for momentum.
      real(kind=8)                :: ch           ! c times |Rib|^1/2 for heat.
      real(kind=8)                :: ee           ! (z/z0)^1/3 -1. for eqn. 20
      !----- Local variables, used by others. ---------------------------------------------!
      real(kind=8)                :: zstar        ! Reference height above displacement
      real(kind=8)                :: zeta0m       ! roughness(momentum)/(Obukhov length).
      real(kind=8)                :: zeta0h       ! roughness(heat)/(Obukhov length).
      real(kind=8)                :: utotal       ! Total wind (actual + convective)
      real(kind=8)                :: uconv        ! Convective velocity
      real(kind=8)                :: uconv_prev   ! Previous convective velocity
      real(kind=8)                :: change       ! Difference in convective velocity
      integer                     :: icnt         ! Iteration counter
      !----- Aux. environment conditions. -------------------------------------------------!
      real(kind=8)                :: thetav_atm   ! Atmos. virtual pot. temp.     [      K]
      real(kind=8)                :: thetav_can   ! Canopy air virtual pot. temp. [      K]
      !----- External functions. ----------------------------------------------------------!
      real(kind=8), external      :: cbrt8        ! Cube root
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the ratio between momentum and heat roughness scales.                      !
      !------------------------------------------------------------------------------------!
      select case (isfclyrm)
      case (4)
         !---------------------------------------------------------------------------------!
         !      Use previous u* to estimate the ratio between roughness scales for         !
         ! momentum (z0m) and roughness scale for enthalpy/vapour/CO2 (z0h), using the     !
         ! parametrisation from ZD98.                                                      ! 
         !---------------------------------------------------------------------------------!
         kin_visc = kin_visc08 * ( 1.d0 + dkin_visc8 * ( temp_can - t008 ) )
         z0hoz0m  = exp(- min( zd98_emax8, zd98_a8 * (ustar*rough/kin_visc)**zd98_b8 ) )
         z0moz0h  = 1.d0 / z0hoz0m
         !---------------------------------------------------------------------------------!
      case default
         !---------------------------------------------------------------------------------!
         !      Assume the same roughness scheme.  Mind that this could be off by multiple !
         ! orders of magnitude.                                                            ! 
         !---------------------------------------------------------------------------------!
         z0hoz0m  = 1.d0
         z0moz0h  = 1.d0 / z0hoz0m
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !----- Find the variables common to both methods. -----------------------------------!
      thetav_atm = theta_atm * (1.d0 + epim18 * shv_atm)
      thetav_can = theta_can * (1.d0 + epim18 * shv_can)
      zstar      = zref - dheight
      zoz0m      = zstar / rough
      lnzoz0m    = log(zoz0m)
      zoz0h      = z0moz0h * zoz0m
      lnzoz0h    = log(zoz0h)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the Bulk Richardson number.  For stable cases, and for L79 in both cases, !
      ! this will be the definitive RiB, whilst this is the first guess, which will be     !
      ! corrected by the convective velocity in the other unstable cases.                  !
      !------------------------------------------------------------------------------------!
      rib        = 2.d0 * grav8 * (zstar-rough) * (thetav_atm-thetav_can)                  &
                 / ( (thetav_atm+thetav_can) * uref * uref)
      stable     = thetav_atm >= thetav_can
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Correct the bulk Richardson number in case it's too stable.  We also define a   !
      ! stable case correction to bring down the stars other than ustar, so the flux       !
      ! doesn't increase for stabler cases (it remains constant).                          !
      !------------------------------------------------------------------------------------!
      if (rib > ribmax8) then
         uuse = sqrt(rib / ribmax8) * uref
         rib  = ribmax8
      else
         uuse = uref
      end if
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Here we find u* and the coefficient to find the other stars based on the       !
      ! chosen surface model.                                                              !
      !------------------------------------------------------------------------------------!
      select case (isfclyrm)
      case (0)
         !----- Use prescribed u* instead of calculating it. ------------------------------!
         ustar    = max(ustmin8,atm_ustar)
         kuoustar = vonk8 * uref / ustar
         !---------------------------------------------------------------------------------!


         !----- Find the dimensionless height. --------------------------------------------!
         zeta   = zoobukhov_ustar8(rib,zstar,rough,zoz0h,lnzoz0h,kuoustar,stable)
         !---------------------------------------------------------------------------------!


         !----- Find the dimensionless roughness. -----------------------------------------!
         zeta0m = rough * zeta / zstar
         zeta0h = z0hoz0m * zeta0m
         !---------------------------------------------------------------------------------!


         !----- Find the coefficient to scale the other stars. ----------------------------!
         c3    = vonk8                                                                     &
               / (tprandtl8 * (lnzoz0h - psih8(zeta,stable) + psih8(zeta0h,stable)))
         !---------------------------------------------------------------------------------!

      case (1)

         !----- Compute the a-square factor and the coefficient to find theta*. -----------!
         a2   = vonk8 * vonk8 / (lnzoz0m * lnzoz0m)
         c1   = a2 * uuse

         if (stable) then
            !----- Stable case ------------------------------------------------------------!

            fm = 1.d0 / (1.d0 + (2.d0 * bl798 * rib / sqrt(1.d0 + dl798 * rib)))
            fh = 1.d0 / (1.d0 + (3.d0 * bl798 * rib * sqrt(1.d0 + dl798 * rib)))

         else

            !------------------------------------------------------------------------------!
            !     Unstable case.  The only difference from the original method is that we  !
            ! no longer assume z >> z0, so the "c" coefficient uses the full z/z0 term.    !
            !------------------------------------------------------------------------------!
            ee = cbrt8(zoz0m) - 1.d0
            c2 = bl798 * a2 * ee * sqrt(ee * abs(rib))
            cm = csm8 * c2
            ch = csh8 * c2
            fm = (1.d0 - 2.d0 * bl798 * rib / (1.d0 + 2.d0 * cm))
            fh = (1.d0 - 3.d0 * bl798 * rib / (1.d0 + 3.d0 * ch))
         end if

         !----- Finding ustar, making sure it is not too small. ---------------------------!
         ustar = max(ustmin8,sqrt(c1 * uuse * fm))
         !----- Finding the coefficient to scale the other stars. -------------------------!
         c3 = c1 * fh / ustar
         !---------------------------------------------------------------------------------!

         !----- Compute zeta from u* and T* -----------------------------------------------!
         zeta = grav8 * vonk8 * c3 * (thetav_atm - thetav_can)                             &
              / (thetav_atm * ustar * ustar)


      case default
         !---------------------------------------------------------------------------------!
         ! 2. Here we use the model proposed by OD95, the standard for MM5, but with some  !
         !    terms that were computed in B71 (namely, the "0" terms), which prevent sin-  !
         !    gularities.                                                                  !
         !    However we know zeta, so zeta0 can be written as z0/z * zeta.                !
         ! 3. Here we use the model proposed by BH91, which is almost the same as the OD95 !
         !    method, except that the stable functions are computed in a more generic way. !
         !    BH91 claim that the oft-used approximation (-beta*zeta) can cause poor       !
         !    ventilation of the stable layer, leading to decoupling between the atmo-     !
         !    sphere and the canopy air space and excessive cooling                        !
         ! 4. Here we use a similar approach as in CLM04, excepth that the momentum flux   !
         !    gradient function for the unstable case for momentum is switched by a power  !
         !    of -1/6 (kind of the square of the heat one).  This is to guarantee that     !
         !    the psi function doesn't have local maxima/minima.                           !
         !---------------------------------------------------------------------------------!
         !----- Initialise uconv. ---------------------------------------------------------!
         uconv      = 0.0
         !----- Check if we need to go through the iterative process. ---------------------!
         if (stable) then
            !------------------------------------------------------------------------------!
            !     Stable case, we don't need to find convective velocity, so we don't need !
            ! the iterative case.                                                          !
            !------------------------------------------------------------------------------!
            !----- Find the dimensionless height. -----------------------------------------!
            zeta   = zoobukhov8(rib,zstar,rough,z0moz0h,zoz0m,lnzoz0m,zoz0h,lnzoz0h,stable)
            !------------------------------------------------------------------------------!


            !----- Find the dimensionless roughness. --------------------------------------!
            zeta0m = rough * zeta / zstar
            zeta0h = z0hoz0m * zeta0m
            !------------------------------------------------------------------------------!



            !----- Find ustar, make sure it is not too small. -----------------------------!
            ustar = max (ustmin8, vonk8 * uuse                                             &
                                / (lnzoz0m - psim8(zeta,stable) + psim8(zeta0m,stable)))
            !------------------------------------------------------------------------------!



            !----- Find the coefficient to scale the other stars. -------------------------!
            c3    = vonk8                                                                  &
                  / (tprandtl8 * (lnzoz0h - psih8(zeta,stable) + psih8(zeta0h,stable)))
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !    Unstable case.  Here we run a few iterations to make sure we correct the  !
            ! bulk Richardson number.  This is really a simple correction, so we don't     !
            ! need uconv to be totally in equilibrium.                                     !
            !------------------------------------------------------------------------------!
            unstable: do icnt=1,6
               !----- Update total winds. -------------------------------------------------!
               uconv_prev = uconv
               utotal     = sqrt(uuse*uuse + uconv_prev * uconv_prev)
               !---------------------------------------------------------------------------!

               !----- Update the Bulk Richardson number. ----------------------------------!
               rib        = 2.d0 * grav8 * (zstar-rough) * (thetav_atm-thetav_can)         &
                          / ( (thetav_atm+thetav_can) * utotal * utotal)
               !---------------------------------------------------------------------------!


               !----- Find the dimensionless height. --------------------------------------!
               zeta   = zoobukhov8(rib,zstar,rough,z0moz0h,zoz0m,lnzoz0m,zoz0h,lnzoz0h     &
                                  ,stable)
               !---------------------------------------------------------------------------!


               !----- Find the dimensionless roughness. -----------------------------------!
               zeta0m = rough * zeta / zstar
               zeta0h = z0hoz0m * zeta0m
               !---------------------------------------------------------------------------!


               !----- Find ustar, making sure it is not too small. ------------------------!
               ustar = max (ustmin8, vonk8 * uuse                                          &
                                   / (lnzoz0m - psim8(zeta,stable) + psim8(zeta0m,stable)))
               !---------------------------------------------------------------------------!


               !----- Find the coefficient to scale the other stars. ----------------------!
               c3    = vonk8                                                               &
                     / (tprandtl8 * (lnzoz0h - psih8(zeta,stable) + psih8(zeta0h,stable)))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Use potential virtual temperature here because convection is related  !
               ! to buoyancy.                                                              !
               !---------------------------------------------------------------------------!
               tstar  = c3 * (thetav_atm - thetav_can )
               !---------------------------------------------------------------------------!


               !----- Estimate the convective velocity. -----------------------------------!
               uconv = vertical_vel_flux8(zeta,ustar) / ustar
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     We are only after a rough estimate of this velocity, so if the        !
               ! difference is less than the RK4 tolerance, then it's enough.              !
               !---------------------------------------------------------------------------!
               change = 2.d0 * abs(uconv-uconv_prev) / (abs(uconv) + abs(uconv_prev))
               if (change < rk4eps) exit unstable
               !---------------------------------------------------------------------------!
            end do unstable
         end if
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !----- Compute the other scales. ----------------------------------------------------!
      qstar = c3 * (shv_atm      - shv_can      )
      tstar = c3 * (theta_atm    - theta_can    )
      estar = c3 * (enthalpy_atm - enthalpy_can )
      cstar = c3 * (co2_atm      - co2_can      )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Compute the bare ground conductance.  This equation is similar to the original, !
      ! except that we don't assume the ratio between the gradient and the characteristic  !
      ! scale to be 0.2; instead we use the actual ratio that is computed here.            !
      !------------------------------------------------------------------------------------!
      ggbare = c3 * ustar
      !------------------------------------------------------------------------------------!

      return
   end subroutine ed_stars8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function determines the wind at a given height, given that the stars are al-  !
   ! ready known, as well as the Richardson number and the zetas.                          !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function reduced_wind8(ustar,zeta,rib,zref,dheight,height,rough)
      use consts_coms    , only : vonk8     ! ! intent(in)
      use canopy_air_coms, only : isfclyrm  & ! intent(in)
                                , bl798     & ! intent(in)
                                , csm8      & ! intent(in)
                                , dl798     & ! intent(in)
                                , ugbmin8   & ! intent(in)
                                , psim8     ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: ustar     ! Friction velocity                   [    m/s]
      real(kind=8), intent(in) :: zeta      ! Normalised height                   [    ---]
      real(kind=8), intent(in) :: rib       ! Bulk Richardson number              [    ---]
      real(kind=8), intent(in) :: zref      ! Reference height                    [      m]
      real(kind=8), intent(in) :: dheight   ! Displacement height                 [      m]
      real(kind=8), intent(in) :: height    ! Height to determine the red. wind   [      m]
      real(kind=8), intent(in) :: rough     ! Roughness scale                     [      m]
      !----- Local variables. -------------------------------------------------------------!
      logical                  :: stable    ! Canopy air space is stable          [    T|F]
      real(kind=8)             :: zetah     ! Zeta for h=height                   [    ---]
      real(kind=8)             :: zeta0     ! Zeta for h=rough                    [    ---]
      real(kind=8)             :: hoz0      ! ((h-d0)/z0)                         [    ---]
      real(kind=8)             :: lnhoz0    ! ln ((h-d0)/z0)                      [    ---]
      real(kind=8)             :: a2        ! Drag coeff. in neutral conditions
      real(kind=8)             :: fm        ! Stability parameter for momentum
      real(kind=8)             :: c2        ! Part of the c coefficient.
      real(kind=8)             :: cm        ! c coefficient times |Rib|^1/2
      real(kind=8)             :: ee        ! (z/z0)^1/3 -1. for eqn. 20 (L79)
      !----- External functions. ----------------------------------------------------------!
      real(kind=8), external   :: cbrt8     ! Cubic root
      !------------------------------------------------------------------------------------!



      !----- Define whether the layer is stable or not. -----------------------------------!
      stable    = rib >= 0.d0
      !------------------------------------------------------------------------------------!



      !----- Find the log for the log-height interpolation of wind. -----------------------!
      hoz0      = (height-dheight)/rough
      lnhoz0    = log(hoz0)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The wind at a given height is found by using the same definition of wind speed !
      ! at a given height.                                                                 !
      !------------------------------------------------------------------------------------!
      select case (isfclyrm)
      case (1) !---- Louis (1979) method. -------------------------------------------------!

         !----- Compute the a-square factor and the coefficient to find theta*. -----------!
         a2   = vonk8 * vonk8 / (lnhoz0 * lnhoz0)

         if (stable) then
            !----- Stable case ------------------------------------------------------------!
            fm = 1.d0 / (1.d0 + (2.d0 * bl798 * rib / sqrt(1.d0 + dl798 * rib)))

         else
            !------------------------------------------------------------------------------!
            !     Unstable case.  The only difference from the original method is that we  !
            ! no longer assume z >> z0, so the "c" coefficient uses the full z/z0 term.    !
            !------------------------------------------------------------------------------!
            ee = cbrt8(hoz0) - 1.d0
            c2 = bl798 * a2 * ee * sqrt(ee * abs(rib))
            cm = csm8 * c2
            fm = (1.d0 - 2.d0 * bl798 * rib / (1.d0 + 2.d0 * cm))
         end if

         !----- Find the wind. ------------------------------------------------------------!
         reduced_wind8 = (ustar/vonk8) * (lnhoz0/sqrt(fm))

      case default  !----- Other methods. -------------------------------------------------!

         !----- Determine zeta for the sought height and for roughness height. ------------!
         zetah = zeta * (height-dheight) / (zref-dheight)
         zeta0 = zeta * rough            / (zref-dheight)
         !---------------------------------------------------------------------------------!

         reduced_wind8 = (ustar/vonk8)                                                     &
                       * (lnhoz0 - psim8(zetah,stable) + psim8(zeta0,stable))

      end select
      !------------------------------------------------------------------------------------!



      !----- Impose the wind to be more than the minimum. ---------------------------------!
      reduced_wind8 = max(reduced_wind8,ugbmin8)
      !------------------------------------------------------------------------------------!

      return
   end function reduced_wind8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Vertical flux, as in:                                                              !
   !                                                                                       !
   !   Manton, M. J., Cotton, W. R., 1977: Parameterization of the atmospheric surface     !
   !      layer.  J. Atm. Sci., 34, 331-334.                                               !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function vertical_vel_flux8(zeta,ustar)
      use consts_coms , only : vonk8 ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)    :: zeta
      real(kind=8), intent(in)    :: ustar
      !----- Local variables --------------------------------------------------------------!
      real(kind=8) :: cx
      real(kind=8) :: psin
      !----- Constants --------------------------------------------------------------------!
      real(kind=8), parameter     :: wtol = 1.d-20
      !------------------------------------------------------------------------------------!


      if (zeta < 0.d0)then
         cx = zeta * sqrt(sqrt(1.d0 - 1.5d1 * zeta))
      else
         cx = zeta / (1.0d0 + 4.7d0 * zeta)
      endif

      psin = sqrt((1.d0-2.86d0 * cx) / (1.d0 + cx * (-5.39d0 + cx * 6.998d0 )))
      vertical_vel_flux8 = ( 2.7d-1 * max(6.25d0 * (1.d0 - cx) * psin,wtol)                &
                           - 1.18d0 * cx * psin) * ustar * ustar

      return
   end function vertical_vel_flux8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Calculate some canopy air space properties, such as the total mass and depth, and !
   ! also the total capacities (carbon, water, and heat).                                  !
   !     Canopy air space depth must change because we are imposing that pressure and      !
   ! density are constants. Since mass and temperature can freely change, in order to      !
   ! guarantee that the ideal gas law and the first law of thermodynamic apply, the volume !
   ! (or the canopy depth) must be allowed to change over time, so work can be done by the !
   ! canopy or into the canopy.                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine can_whccap8(can_rhos,can_dmol,can_depth,wcapcan,hcapcan,ccapcan              &
                        ,wcapcani,hcapcani,ccapcani)
      use consts_coms, only : mmdry8  & ! intent(in)
                            , mmdryi8 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)  :: can_rhos  ! Canopy air mass density            [  kg/m3]
      real(kind=8), intent(in)  :: can_dmol  ! Canopy dry-air molar density       [ mol/m3]
      real(kind=8), intent(in)  :: can_depth ! Depth of canopy air space          [      m]
      real(kind=8), intent(out) :: wcapcan   ! Water capacity - canopy air space  [  kg/m2]
      real(kind=8), intent(out) :: hcapcan   ! Enthalpy capacity - CAS            [  kg/m2]
      real(kind=8), intent(out) :: ccapcan   ! CO2 capacity - CAS                 [ mol/m2]
      real(kind=8), intent(out) :: wcapcani  ! Inverse of water capacity          [  m2/kg]
      real(kind=8), intent(out) :: hcapcani  ! Inverse of enthalpy capcity        [  m2/kg]
      real(kind=8), intent(out) :: ccapcani  ! Inverse of CO2 capacity            [ m2/mol]
      !------------------------------------------------------------------------------------!

      !----- Find the water capacity and its inverse. -------------------------------------!
      wcapcan  = can_rhos * can_depth
      wcapcani = 1.d0 / wcapcan
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !       Because we track specific enthalpy [J/kg], the value is the same as water    !
      ! capacity.                                                                          !
      !------------------------------------------------------------------------------------!
      hcapcan  = wcapcan
      hcapcani = wcapcani
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The CO2 capacity must be in mol/m2 rather than kg/m2, since CO2 variable is in !
      ! umol/mol.                                                                          !
      !------------------------------------------------------------------------------------!
      ccapcan  = can_dmol * can_depth
      ccapcani = 1.d0 / ccapcan
      !------------------------------------------------------------------------------------!

      return
   end subroutine can_whccap8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine computes the aerodynamic conductance between leaf and canopy     !
   ! air space for both heat and water vapour, based on:                                   !
   !                                                                                       !
   ! L95 - Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf        !
   !       nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to !
   !       canopies.  Plant, Cell and Environ., 18, 1183-1200.                             !
   ! M08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,     !
   !       3rd. edition, Academic Press, Amsterdam, 418pp.  (Mostly Chapter 10).           !
   !                                                                                       !
   ! Notice that the units are somewhat different from L95.                                !
   ! - gbh is in J/(K m2 s), and                                                           !
   ! - gbw is in kg_H2O/m2/s.                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine leaf_aerodynamic_conductances8(ipft,veg_wind,leaf_temp,can_temp              &
                                            ,can_rhos,can_cp,leaf_gbh,leaf_gbw,reynolds    &
                                            ,grashof,nusselt_free,nusselt_forced)
      use pft_coms       , only : leaf_width    ! ! intent(in)
      use canopy_air_coms, only : aflat_lami8   & ! intent(in)
                                , nflat_lami8   & ! intent(in)
                                , aflat_turb8   & ! intent(in)
                                , nflat_turb8   & ! intent(in)
                                , bflat_lami8   & ! intent(in)
                                , mflat_lami8   & ! intent(in)
                                , bflat_turb8   & ! intent(in)
                                , mflat_turb8   & ! intent(in)
                                , gbhmos_min8   ! ! intent(in)
      use consts_coms    , only : t008          & ! intent(in)
                                , grav8         & ! intent(in)
                                , kin_visc08    & ! intent(in)
                                , dkin_visc8    & ! intent(in)
                                , th_diff08     & ! intent(in)
                                , dth_diff8     ! ! intent(in)
      use physiology_coms, only : gbh_2_gbw8    ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                      :: ipft            ! Plant functional type   [      ---]
      real(kind=8)   , intent(in)  :: veg_wind        ! Wind at cohort height   [      m/s]
      real(kind=8)   , intent(in)  :: leaf_temp       ! Leaf temperature        [        K]
      real(kind=8)   , intent(in)  :: can_temp        ! Canopy air temperature  [        K]
      real(kind=8)   , intent(in)  :: can_rhos        ! Canopy air density      [    kg/m2]
      real(kind=8)   , intent(in)  :: can_cp          ! Canopy air spec. heat   [   J/kg/K]
      real(kind=8)   , intent(out) :: leaf_gbh        ! Heat  conductance       [ J/K/m2/s]
      real(kind=8)   , intent(out) :: leaf_gbw        ! Water conductance       [  kg/m2/s]
      real(kind=8)   , intent(out) :: grashof         ! Grashof number          [      ---]
      real(kind=8)   , intent(out) :: reynolds        ! Reynolds number         [      ---]
      real(kind=8)   , intent(out) :: nusselt_free    ! Nusselt number (free)   [      ---]
      real(kind=8)   , intent(out) :: nusselt_forced  ! Nusselt number (forced) [      ---]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                 :: lwidth          ! Leaf width              [        m]
      real(kind=8)                 :: kin_visc        ! Kinematic viscosity     [     m2/s]
      real(kind=8)                 :: th_diff         ! Kinematic viscosity     [     m2/s]
      real(kind=8)                 :: th_expan        ! Thermal expansion       [      1/K]
      real(kind=8)                 :: gr_coeff        ! grav*th_expan/kin_visc2 [   1/K/m2]
      real(kind=8)                 :: nusselt_lami    ! Nusselt number (laminar)[      ---]
      real(kind=8)                 :: nusselt_turb    ! Nusselt number (turb.)  [      ---]
      real(kind=8)                 :: forced_gbh_mos  ! Forced convection cond. [      m/s]
      real(kind=8)                 :: free_gbh_mos    ! Free convection cond.   [      m/s]
      real(kind=8)                 :: gbh_mos         ! Total convection cond.  [      m/s]
      !------------------------------------------------------------------------------------!


      !----- Save the leaf width of this PFT. ---------------------------------------------!
      lwidth = dble(leaf_width(ipft))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Compute kinematic viscosity, thermal diffusivity, and expansion coefficient as !
      ! functions of temperature.  Here we use the canopy air space temperature because    !
      ! this is the representative temperature of the fluid.                               !
      !                                                                                    !
      !     Kinematic viscosity and thermal diffusivity are determined from MU08, see      !
      ! discussion on page 32.  Thermal expansion is assumed to be of an ideal gas (1/T),  !
      ! like in Dufour and van Mieghem (1975), for example.                                !
      !------------------------------------------------------------------------------------!
      th_expan = 1.d0 / can_temp
      !----- kin_visc and th_diff are assumed linear functions of temperature. ------------!
      kin_visc = kin_visc08 * ( 1.d0 + dkin_visc8 * ( can_temp - t008 ) )
      th_diff  = th_diff08  * ( 1.d0 + dth_diff8  * ( can_temp - t008 ) )
      !------------------------------------------------------------------------------------!
      !    Grashof coefficient (a*g/nu2) in MU08's equation 10.8.                          !
      !------------------------------------------------------------------------------------!
      gr_coeff = th_expan * grav8  / ( kin_visc * kin_visc )
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the conductance, in m/s, associated with forced convection.               !
      !------------------------------------------------------------------------------------!
      !----- 1. Compute the Reynolds number. ----------------------------------------------!
      reynolds        = veg_wind * lwidth / th_diff
      !----- 2. Compute the Nusselt number for both the laminar and turbulent case. -------!
      nusselt_lami    = aflat_lami8 * reynolds ** nflat_lami8
      nusselt_turb    = aflat_turb8 * reynolds ** nflat_turb8
      !----- 3. The right Nusselt number is the largest of the both. ----------------------!
      nusselt_forced  = max(nusselt_lami,nusselt_turb)
      !----- 4. The conductance is given by MU08 - equation 10.4 --------------------------!
      forced_gbh_mos  = th_diff * nusselt_forced / lwidth
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the conductance, in m/s,  associated with free convection.                !
      !------------------------------------------------------------------------------------!
      !----- 1. Find the Grashof number. --------------------------------------------------!
      grashof         = gr_coeff * abs(leaf_temp - can_temp) * lwidth * lwidth * lwidth
      !----- 2. Compute the Nusselt number for both the laminar and turbulent case. -------!
      nusselt_lami    = bflat_lami8 * grashof ** mflat_lami8
      nusselt_turb    = bflat_turb8 * grashof ** mflat_turb8
      !----- 3. The right Nusselt number is the largest of the both. ----------------------!
      nusselt_free    = max(nusselt_lami,nusselt_turb)
      !----- 4. The conductance is given by MU08 - equation 10.4 --------------------------!
      free_gbh_mos    = th_diff * nusselt_free / lwidth
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The heat conductance for the thermodynamic budget is the sum of conductances,  !
      ! because we assume both forms of convection happen parallelly.  The conversion from !
      ! heat to water conductance (in m/s) can be found in L95, page 1198, after equation  !
      ! E5.  For the ED purposes, the output variables are converted to the units of       !
      ! entropy and water fluxes [J/K/m2/s and kg/m2/s, respectively].                     !
      !------------------------------------------------------------------------------------!
      gbh_mos  = max(gbhmos_min8, free_gbh_mos + forced_gbh_mos)
      leaf_gbh =              gbh_mos * can_rhos * can_cp
      leaf_gbw = gbh_2_gbw8 * gbh_mos * can_rhos
      !------------------------------------------------------------------------------------!

      return
   end subroutine leaf_aerodynamic_conductances8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine computes the aerodynamic conductance between wood and canopy     !
   ! air space for both heat and water vapour, based on:                                   !
   !                                                                                       !
   ! L95 - Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf        !
   !       nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to !
   !       canopies.  Plant, Cell and Environ., 18, 1183-1200.                             !
   ! M08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,     !
   !       3rd. edition, Academic Press, Amsterdam, 418pp.  (Mostly Chapter 10).           !
   !                                                                                       !
   ! Notice that the units are somewhat different from L95.                                !
   ! - gbh is in J/(K m2 s), and                                                           !
   ! - gbw is in kg_H2O/m2/s.                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine wood_aerodynamic_conductances8(ipft,veg_wind,wood_temp,can_temp,can_rhos     &
                                            ,can_cp,wood_gbh,wood_gbw,reynolds,grashof     &
                                            ,nusselt_free,nusselt_forced)
      use canopy_air_coms, only : ocyli_lami8   & ! intent(in)
                                , acyli_lami8   & ! intent(in)
                                , ncyli_lami8   & ! intent(in)
                                , acyli_turb8   & ! intent(in)
                                , ocyli_turb8   & ! intent(in)
                                , ncyli_turb8   & ! intent(in)
                                , bcyli_lami8   & ! intent(in)
                                , mcyli_lami8   & ! intent(in)
                                , bcyli_turb8   & ! intent(in)
                                , mcyli_turb8   & ! intent(in)
                                , gbhmos_min8   ! ! intent(in)
      use consts_coms    , only : t008          & ! intent(in)
                                , grav8         & ! intent(in)
                                , kin_visc08    & ! intent(in)
                                , dkin_visc8    & ! intent(in)
                                , th_diff08     & ! intent(in)
                                , dth_diff8     ! ! intent(in)
      use physiology_coms, only : gbh_2_gbw8    ! ! intent(in)
      use pft_coms       , only : branch_diam   ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer        , intent(in)  :: ipft            ! Current PFT             [       --]
      real(kind=8)   , intent(in)  :: veg_wind        ! Wind at cohort height   [      m/s]
      real(kind=8)   , intent(in)  :: wood_temp       ! Wood temperature        [        K]
      real(kind=8)   , intent(in)  :: can_temp        ! Canopy air temperature  [        K]
      real(kind=8)   , intent(in)  :: can_rhos        ! Canopy air density      [    kg/m2]
      real(kind=8)   , intent(in)  :: can_cp          ! Canopy air spec. heat   [   J/kg/K]
      real(kind=8)   , intent(out) :: wood_gbh        ! Heat  conductance       [ J/K/m2/s]
      real(kind=8)   , intent(out) :: wood_gbw        ! Water conductance       [  kg/m2/s]
      real(kind=8)   , intent(out) :: grashof         ! Grashof number          [      ---]
      real(kind=8)   , intent(out) :: reynolds        ! Reynolds number         [      ---]
      real(kind=8)   , intent(out) :: nusselt_free    ! Nusselt number (free)   [      ---]
      real(kind=8)   , intent(out) :: nusselt_forced  ! Nusselt number (forced) [      ---]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                 :: kin_visc        ! Kinematic viscosity     [     m2/s]
      real(kind=8)                 :: th_diff         ! Kinematic viscosity     [     m2/s]
      real(kind=8)                 :: th_expan        ! Thermal expansion       [      1/K]
      real(kind=8)                 :: gr_coeff        ! grav*th_expan/kin_visc2 [   1/K/m2]
      real(kind=8)                 :: nusselt_lami    ! Nusselt number (laminar)[      ---]
      real(kind=8)                 :: nusselt_turb    ! Nusselt number (turb.)  [      ---]
      real(kind=8)                 :: forced_gbh_mos  ! Forced convection cond. [      m/s]
      real(kind=8)                 :: free_gbh_mos    ! Free convection cond.   [      m/s]
      real(kind=8)                 :: gbh_mos         ! Total convection cond.  [      m/s]
      real(kind=8)                 :: w_diam          ! Typical branch diameter [        m]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Characteristic branch diameter.  We currently assume a fixed size for         !
      ! simplicity, but we could make it dynamic by applying a branch size distribution    !
      ! model.                                                                             !
      !------------------------------------------------------------------------------------!
      w_diam = dble(branch_diam(ipft))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Compute kinematic viscosity, thermal diffusivity, and expansion coefficient as !
      ! functions of temperature.  Here we use the canopy air space temperature because    !
      ! this is the representative temperature of the fluid.                               !
      !                                                                                    !
      !     Kinematic viscosity and thermal diffusivity are determined from MU08, see      !
      ! discussion on page 32.  Thermal expansion is assumed to be of an ideal gas (1/T),  !
      ! like in Dufour and van Mieghem (1975), for example.                                !
      !------------------------------------------------------------------------------------!
      th_expan = 1.d0 / can_temp
      !----- kin_visc and th_diff are assumed linear functions of temperature. ------------!
      kin_visc = kin_visc08 * ( 1.d0 + dkin_visc8 * ( can_temp - t008 ) )
      th_diff  = th_diff08  * ( 1.d0 + dth_diff8  * ( can_temp - t008 ) )
      !------------------------------------------------------------------------------------!
      !    Grashof coefficient (a*g/nu2) in MU08's equation 10.8.                          !
      !------------------------------------------------------------------------------------!
      gr_coeff = th_expan * grav8  / ( kin_visc * kin_visc )
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the conductance, in m/s, associated with forced convection.               !
      !------------------------------------------------------------------------------------!
      !----- 1. Compute the Reynolds number. ----------------------------------------------!
      reynolds        = veg_wind * w_diam / th_diff
      !----- 2. Compute the Nusselt number for both the laminar and turbulent case. -------!
      nusselt_lami    = ocyli_lami8 + acyli_lami8 * reynolds ** ncyli_lami8
      nusselt_turb    = ocyli_turb8 + acyli_turb8 * reynolds ** ncyli_turb8
      !----- 3. The right Nusselt number is the largest of the both. ----------------------!
      nusselt_forced  = max(nusselt_lami,nusselt_turb)
      !----- 5. The conductance is given by MU08 - equation 10.4 --------------------------!
      forced_gbh_mos  = th_diff * nusselt_forced / w_diam
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the conductance, in m/s,  associated with free convection.                !
      !------------------------------------------------------------------------------------!
      !----- 1. Find the Grashof number. --------------------------------------------------!
      grashof         = gr_coeff * abs(wood_temp - can_temp) * w_diam * w_diam * w_diam
      !----- 2. Compute the Nusselt number for both the laminar and turbulent case. -------!
      nusselt_lami    = bcyli_lami8 * grashof ** mcyli_lami8
      nusselt_turb    = bcyli_turb8 * grashof ** mcyli_turb8
      !----- 3. The right Nusselt number is the largest of the both. ----------------------!
      nusselt_free    = max(nusselt_lami,nusselt_turb)
      !----- 5. The conductance is given by MU08 - equation 10.4 --------------------------!
      free_gbh_mos    = th_diff * nusselt_free / w_diam
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The heat conductance for the thermodynamic budget is the sum of conductances,  !
      ! because we assume both forms of convection happen parallelly.  The conversion from !
      ! heat to water conductance (in m/s) can be found in L95, page 1198, after equation  !
      ! E5.  For the ED purposes, the output variables are converted to the units of       !
      ! entropy and water fluxes [J/K/m2/s and kg/m2/s, respectively].                     !
      !------------------------------------------------------------------------------------!
      gbh_mos  = max(gbhmos_min8, free_gbh_mos + forced_gbh_mos)
      wood_gbh =              gbh_mos * can_rhos * can_cp
      wood_gbw = gbh_2_gbw8 * gbh_mos * can_rhos
      !------------------------------------------------------------------------------------!

      return
   end subroutine wood_aerodynamic_conductances8
   !=======================================================================================!
   !=======================================================================================!
end module canopy_struct_dynamics
!==========================================================================================!
!==========================================================================================!

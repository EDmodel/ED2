!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This module contains some variables used in LEAF-3.                                   !
!------------------------------------------------------------------------------------------!
module leaf_coms

   use grid_dims
   use rconstants, only: grav      & ! intent(in)
                       , vonk      & ! intent(in)
                       , twothirds ! ! intent(in)

   integer :: niter_leaf   ! ! number of leaf timesteps

   real    :: dtll         & ! leaf timestep
            , dtll_factor  & ! leaf timestep factor (leaf timestep / model timestep)
            , dtllowcc     & ! leaf timestep   / (can_depth * can_rhos)
            , dtlloccc     & ! mmdry * leaf timestep   / (can_depth * can_rhos)

            , atm_up       & ! U velocity at top of surface layer                [     m/s]
            , atm_vp       & ! V velocity at top of surface layer                [     m/s]
            , atm_theta    & ! potential temperature at top of surface layer     [       K]
            , atm_temp     & ! temperature at top of surface layer               [       K]
            , atm_rvap     & ! vapor mixing ratio at top of surface layer        [   kg/kg]
            , atm_shv      & ! specific humidity at top of surface layer         [   kg/kg]
            , atm_co2      & ! CO2 mixing ratio at top of surface layer          [µmol/mol]
            , atm_enthalpy & ! atmospheric enthalpy                              [    J/kg]
            , geoht        & ! height at top of surface layer                    [       m]
            , atm_exner    & ! "Exner" function at surface (Exner/cp)            [     ---]
            , atm_prss     & ! pressure at surface                               [      Pa]
            , vels         & ! wind speed at top of surface layer                [       m]
            , vels_pat     & ! vels with patch-dep. ubmin imposed as minimum     [     m/s]
            , pcpgl        & ! precip mass in leaf timestep                      [   kg/m²]
            , qpcpgl       & ! precip energy in leaf timestep                    [    J/m²]
            , dpcpgl       & ! precip depth in leaf timestep                     [       m]
            , pcpgc        & ! precip mass in canopy timestep                    [   kg/m²]
            , qpcpgc       & ! precip energy in canopy timestep                  [    J/m²]
            
            , snowfac      & ! fraction of vegetation height covered by sfcwater [     ---]
            , vf           & ! product of veg_fracarea and (1-snowfac)           [     ---]
            , can_temp     & ! canopy air temperature                            [       K]
            , can_rhos     & ! canopy air density                                [   kg/m³]
            , can_shv      & ! canopy air specific humidity                      [   kg/kg]
            , can_enthalpy & ! canopy air enthalpy                               [    J/kg]
            , can_depth    & ! canopy depth                                      [       m]
            , veg_temp     & ! vegetation temperature                            [       K]
            , veg_fliq     & ! liquid fraction of vegetation surface water       [      --]
            , estar        & ! enthalpy characteristic friction scale            [    J/kg]
            , qstar        & ! specific humidity characteristic friction scale   [   kg/kg]
            , timefac_sst  & ! time interpolation factor for SST                 [     ---]
            
            , rb           & ! vegetation aerodynamic resistance
            , rd           & ! canopy to ground aerodynamic resistance
            , rdi          & ! canopy to ground aerodynamic conductance
            , rho_ustar    & ! canopy density time friction velocity
            , rshort_g     & ! net SW radiation absorbed by grnd
            , rshort_v     & ! net SW radiation absorbed by veg
            , rshort_a     & ! net SW radiation reflected to atm by veg plus grnd
            
            , rlonga_v     & ! net atm LW radiation absorbed by veg
            , rlonga_gs    & ! net atm LW radiation absorbed by grnd OR snow
            , rlongv_gs    & ! net veg LW radiation absorbed by grnd OR snow
            , rlongv_a     & ! net veg LW radiation to atm
            , rlonggs_v    & ! net grnd OR snow LW radiation absorbed by veg
            , rlonggs_a    & ! net grnd OR snow LW radiation to atm
            , rlonga_a     & ! net atm LW radiation refl. to atm by veg plus grnd OR snow

            , hflxgc       & ! sensible heat from ground to canopy              [   J/m²/s]
            , wflxgc       & ! water vapor from ground to canopy                [  kg/m²/s]
            , qwflxgc      & ! latent heat from ground to canopy                [   J/m²/s]
            , cflxgc       & ! carbon from ground to canopy                     [µmol/m²/s]
            , eflxac       & ! enthalpy flux from atmosphere to canopy          [   J/m²/s]
            , wflxac       & ! water flux from atmosphere to canopy             [  kg/m²/s]
            , cflxac       & ! carbon flux from atmosphere to canopy            [µmol/m²/s]
            , hflxvc       & ! sensible heat from vegetation to canopy          [   J/m²/s]
            , wflxvc       & ! water vapor from veg. to canopy (evaporation)    [  kg/m²/s]
            , qwflxvc      & ! latent heat from veg. to canopy (evaporation)    [   J/m²/s]
            , cflxvc       & ! carbon from vegetation to canopy                 [µmol/m²/s]
            , cflxcv       & ! carbon from canopy to vegetation                 [µmol/m²/s]
            , transp_loc   & ! water flux due to transpiration                  [  kg/m²/s]
            , qtransp_loc  & ! latent heat flux due to transpiration            [  kg/m²/s]
            , transp_tot   & ! water flux due to transpiration in 1 leaf ts.    [  kg/m²/s]
            , wshed_tot    & ! water shed from vegetation to ground             [  kg/m²/s]
            , qwshed_tot   & ! energy from shed water                           [   J/m²/s]
            , dewgnd_tot   & ! dew formation on ground                          [  kg/m²/s]
            , qdewgnd_tot  ! ! energy from dew formation on ground              [   J/m²/s]

   real, save, allocatable, dimension(:) ::  &
              dslz            & ! soil layer thickness at T point
            , dslzi           & ! (1. / soil layer thickness at T point)
            , dslzidt         & ! (dtll / soil layer thickness at T point)
            , slzt            & ! soil depth at T point
            , dslzt           & ! soil layer thickness at M point
            , dslzti          & ! (1. / soil layer thickness at M point)
            , dslztidt        & ! (dtll / soil layer thickness at M point)

            , rshort_s        & ! net SW radiation absorbed by snow layers (mzs)
            , sfcw_energy_ext & ! Extensive sfc. water internal energy (mzs)
            , tempk           & ! diagnosed temp (K) of soil and sfcwater
            , fracliq         & ! diagnosed liquid fraction (soil_water & sfcwater_mass)
      
            , hfluxgsc        & ! sensible heat flux between soil, sfcwater, canopy
            , psiplusz        & ! soil water potential plus geopotential [m]
            , half_soilair    & ! half of available airspace in soil [m]
            , rfactor         & ! soil, sfcwater thermal resistance    
            , wflux           & ! soil water flux [m]
            , soil_liq        & ! soil liquid water content [m]
            , qwflux          ! ! soil energy flux from water flux [J/m2] 
   !---------------------------------------------------------------------------------------!



   !----- Some dimensions -----------------------------------------------------------------!
   integer, parameter :: nstyp     = 12 ! # of soil types
   integer, parameter :: nvtyp     = 20 ! # of land use types
   integer, parameter :: nvtyp_teb =  1 ! # of TEB extra land use types (21 - Very urban).
   !---------------------------------------------------------------------------------------!



   !----- Soil properties -----------------------------------------------------------------!
   real, dimension(nstyp)           :: slden,slcpd,slbs,slcond,sfldcap,slcons,slmsts,slpots
   real, dimension(nstyp)           :: ssand,sclay,sorgan,sporo,soilwp,soilcp,slfc,emisg
   real, dimension(nstyp)           :: slcons00,slcons0,fhydraul
   real, dimension(nstyp)           :: soilcond0,soilcond1,soilcond2
   real, dimension(nzgmax,nstyp)    :: slcons1
   !---------------------------------------------------------------------------------------!



   !----- Plant properties ----------------------------------------------------------------!
   real   , dimension(nvtyp+nvtyp_teb)        :: albv_green,albv_brown,emisv,sr_max,tai_max
   real   , dimension(nvtyp+nvtyp_teb)        :: sai,veg_clump,veg_frac,veg_ht,glai_max
   real   , dimension(nvtyp+nvtyp_teb)        :: rcmin,dead_frac
   integer, dimension(nvtyp+nvtyp_teb)        :: kroot
   real   , dimension(nzgmax,nvtyp+nvtyp_teb) :: root
   !---------------------------------------------------------------------------------------!


   !----- Other variables -----------------------------------------------------------------!
   real                             :: cmin,corg,cwat,cair,cka,ckw
   !---------------------------------------------------------------------------------------!


   !----- Roughness -----------------------------------------------------------------------!
   real, parameter :: z0fac_water    = .016 / grav ! Coefficient before ustar²
   real, parameter :: min_waterrough = .0001       ! Minimum water roughness height
   real, parameter :: waterrough     = .0001       ! Water roughness height
   real, parameter :: snowrough      = .001        ! Snow roughness height
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Vegetation heat capacity.  A constant, but it could be scaled by LAI, height,     !
   ! etc...                                                                                !
   !---------------------------------------------------------------------------------------!
   real, parameter :: hcapveg_ref  = 3.e3      ! [J/m2/K]
   real, parameter :: hcapveg_hmin = 1.5
   !---------------------------------------------------------------------------------------!
   !     Minimum canopy air space depth.                                                   !
   !---------------------------------------------------------------------------------------!
   real, parameter :: can_depth_min = 5.0 ! [J/m2/K]

   !----- Some constants to ensure the model good behaviour -------------------------------!
   real, parameter :: min_sfcwater_mass  = 1.e-6
   real, parameter :: min_sfcwater_depth = 1.e-9
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Speed-related minimum values we will consider.                                    !
   !---------------------------------------------------------------------------------------!
   real, parameter :: ubmin    = 0.25  ! Minimum velocity                        [     m/s]
   real, parameter :: ustmin   = 0.025 ! Minimum ustar                           [     m/s]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Constants for surface layer models.                                              !
   !---------------------------------------------------------------------------------------!
   !----- Louis (1979) model. -------------------------------------------------------------!
   real, parameter :: bl79       = 5.0    ! b prime parameter
   real, parameter :: csm        = 7.5    ! C* for momentum (eqn. 20, not co2 char. scale)
   real, parameter :: csh        = 5.0    ! C* for heat (eqn.20, not co2 char. scale)
   real, parameter :: dl79       = 5.0    ! ???
   !----- Oncley and Dudhia (1995) model. -------------------------------------------------!
   real, parameter :: bbeta      = 5.0    ! Beta used by Businger et al. (1971)
   real, parameter :: ribmaxod95 = 0.20   ! Maximum bulk Richardson number
   !----- Beljaars and Holtslag (1991) model. ---------------------------------------------!
   real, parameter :: abh91       = -1.00         ! -a from equation  (28) and (32)
   real, parameter :: bbh91       = -twothirds    ! -b from equation  (28) and (32)
   real, parameter :: cbh91       =  5.0          !  c from equations (28) and (32)
   real, parameter :: dbh91       =  0.35         !  d from equations (28) and (32)
   real, parameter :: ebh91       = -twothirds    ! the 2/3 factor in equation (32)
   real, parameter :: fbh91       =  1.50         ! exponent in equation (32)
   real, parameter :: cod         = cbh91/dbh91   ! c/d
   real, parameter :: bcod        = bbh91 * cod   ! b*c/d
   real, parameter :: fm1         = fbh91 - 1.0   ! f-1
   real, parameter :: ate         = abh91 * ebh91 ! a * e
   real, parameter :: atetf       = ate   * fbh91 ! a * e * f
   real, parameter :: z0moz0h     = 1.0           ! z0(M)/z0(h)
   real, parameter :: z0hoz0m     = 1. / z0moz0h  ! z0(M)/z0(h)
   real, parameter :: ribmaxbh91  = 6.00          ! Maximum bulk Richardson number
   !----- Used by OD95 and BH91. ----------------------------------------------------------!
   real, parameter :: gamm       = 14.0   ! Gamma used by Businger et al. (1971) - momentum.
   real, parameter :: gamh       = 14.0   ! Gamma used by Businger et al. (1971) - heat.
   real, parameter :: tprandtl   = 1.00   ! Turbulent Prandtl number.
   real, parameter :: vkopr      = vonk/tprandtl ! von Karman / turbulent Prandtl
   !---------------------------------------------------------------------------------------!



   !----- Parameters that used to be in LEAF-3 (leaftw) -----------------------------------!


   !---------------------------------------------------------------------------------------!
   real, parameter                   :: exar=2.5
   real, parameter                   :: covr=2.16


   !---------------------------------------------------------------------------------------!
   !     Note: c1=261.5*sqrt((1.-exp(-2.*exar))/(2.*exar))                                 !
   !     from Lee's dissertation, Eq. 3.36.  The factor of 261.5 is                        !
   !     100 * ln((h-d)/zo) / vonk   where d = .63 * h and zo = .13 * h.                   !
   !     The factor of 100 is 1/L in Eq. 3.37.  Thus, c1 * ustar is the                    !
   !     total expression inside the radical in Eq. 3.37.                                  !
   !---------------------------------------------------------------------------------------!
   real, parameter                   :: c1=116.6
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   real, parameter                   :: brad =   196.0,  srad =   0.047
   real, parameter                   :: btlo =   281.5,  stlo =   0.26
   real, parameter                   :: bthi =   310.1,  sthi =  -0.124
   real, parameter                   :: bvpd =  4850.0,  svpd =  -0.0051
   real, parameter                   :: bswp = -1.07e6,  sswp =   7.42e-6
   !---------------------------------------------------------------------------------------!

   contains
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine allocates some specific scratch arrays used by LEAF-3.              !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_leafcol(nzg,nzs)

      implicit none
      integer, intent(in) :: nzg,nzs

      !----- Allocate leaf column arrays. -------------------------------------------------!
      allocate (dslz            (nzg)      )
      allocate (dslzi           (nzg)      )
      allocate (dslzidt         (nzg)      )
      allocate (slzt            (nzg)      )
      allocate (dslzt           (nzg)      )
      allocate (dslzti          (nzg)      )
      allocate (dslztidt        (nzg)      )

      allocate (rshort_s        (nzs)      )
      allocate (sfcw_energy_ext (nzs)      )
      allocate (tempk           (nzg+nzs)  )
      allocate (fracliq         (nzg+nzs)  )

      allocate (hfluxgsc        (nzg+nzs+1))
      allocate (psiplusz        (nzg)      )
      allocate (half_soilair    (nzg)      )
      allocate (rfactor         (nzg+nzs)  )
      allocate (wflux           (nzg+1)    )
      allocate (qwflux          (nzg+1)    )
      allocate (soil_liq        (nzg)      )
               
      return
   end subroutine alloc_leafcol   
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the stability  correction function for momentum.            !
   !---------------------------------------------------------------------------------------!
   real function psim(zeta,stable)
      use rconstants, only : halfpi
      use mem_leaf  , only : istar
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: zeta   ! z/L, z is the height, and L the Obukhov length [ ---]
      logical, intent(in) :: stable ! Flag... This surface layer is stable           [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: xx
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (istar)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            psim = - bbeta * zeta 
         case (3,4) !----- Beljaars and Holtslag (1991). ----------------------------------!
            psim = abh91 * zeta                                                            &
                 + bbh91 * (zeta - cod) * exp(max(-38.,-dbh91 * zeta))                     &
                 + bcod
         end select
      else
         !----- Unstable case, both papers use the same expression. -----------------------!
         xx   = sqrt(sqrt(1.0 - gamm * zeta))
         psim = log(0.125 * (1.0+xx) * (1.0+xx) * (1.0 + xx*xx)) - 2.0*atan(xx) + halfpi
      end if
      return
   end function psim
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the stability  correction function for heat (and vapour,    !
   ! and carbon dioxide too.)                                                              !
   !---------------------------------------------------------------------------------------!
   real function psih(zeta,stable)
      use mem_leaf  , only : istar
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: zeta   ! z/L, z is the height, and L the Obukhov length [ ---]
      logical, intent(in) :: stable ! Flag... This surface layer is stable           [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: yy
      !------------------------------------------------------------------------------------!

      if (stable) then
         select case (istar)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            psih = - bbeta * zeta 
         case (3,4) !----- Beljaars and Holtslag (1991). ----------------------------------!
            psih = 1.0 - (1.0 + ate * zeta)**fbh91                                         &
                 + bbh91 * (zeta - cod) * exp(max(-38.,-dbh91 * zeta)) + bcod
         end select
      else
         !----- Unstable case, both papers use the same expression. -----------------------!
         yy   = sqrt(1.0 - gamh * zeta)
         psih = log(0.25 * (1.0+yy) * (1.0+yy))
      end if
   end function psih
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the derivative of the stability correction function for     !
   ! momentum with respect to zeta.                                                        !
   !---------------------------------------------------------------------------------------!
   real function dpsimdzeta(zeta,stable)
      use mem_leaf  , only : istar
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: zeta   ! z/L, z is the height, and L the Obukhov length [ ---]
      logical, intent(in) :: stable ! Flag... This surface layer is stable           [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: xx
      !------------------------------------------------------------------------------------!

      if (stable) then
         select case (istar)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            dpsimdzeta = - bbeta 
         case (3,4) !----- Beljaars and Holtslag (1991). ----------------------------------!
            dpsimdzeta = abh91 + bbh91 * (1.0 - dbh91 * zeta + cbh91)                      &
                               * exp(max(-38.,-dbh91 * zeta))
         end select
      else
         !----- Unstable case, both papers use the same expression. -----------------------!
         xx         = sqrt(sqrt(1.0 - gamm * zeta))
         dpsimdzeta = - gamm / (xx * (1.0+xx) * (1.0 + xx*xx)) 
      end if

      return
   end function dpsimdzeta
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the derivative of the stability correction function for     !
   ! heat/moisture/CO2 with respect to zeta.                                               !
   !---------------------------------------------------------------------------------------!
   real function dpsihdzeta(zeta,stable)
      use mem_leaf  , only : istar
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: zeta   ! z/L, z is the height, and L the Obukhov length [ ---]
      logical, intent(in) :: stable ! Flag... This surface layer is stable           [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: yy
      !------------------------------------------------------------------------------------!

      if (stable) then
         select case (istar)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            dpsihdzeta = - bbeta
         case (3,4) !----- Beljaars and Holtslag (1991). ----------------------------------!
            dpsihdzeta = - atetf * (1.0 + ate * zeta)**fm1                                 &
                         + bbh91 * (1.0 - dbh91 * zeta + cbh91)                            &
                         * exp(max(-38.,-dbh91 * zeta))
         end select
      else
         !----- Unstable case, both papers use the same expression. -----------------------!
         yy   = sqrt(1.0 - gamh * zeta)
         dpsihdzeta = -gamh / (yy * (1.0 + yy))
      end if

      return
   end function dpsihdzeta
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the value of zeta for a given Richardson number, reference    !
   ! height and the roughness scale.  This is solved by using the definition of Obukhov    !
   ! length scale as stated in Louis (1979) equation (10), modified to define z/L rather   !
   ! than L.  The solution is found  iteratively since it's not a simple function to       !
   ! invert.  It tries to use Newton's method, which should take care of most cases.  In   !
   ! the unlikely case in which Newton's method fails, switch back to modified Regula      !
   ! Falsi method (Illinois).                                                              !
   !---------------------------------------------------------------------------------------!
   real function zoobukhov(rib,zref,rough,zoz0m,lnzoz0m,zoz0h,lnzoz0h,stable)
      use therm_lib, only : toler  & ! intent(in)
                          , maxfpo & ! intent(in)
                          , maxit  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: rib       ! Bulk Richardson number                    [   ---]
      real   , intent(in) :: zref      ! Reference height                          [     m]
      real   , intent(in) :: rough     ! Roughness length scale                    [     m]
      real   , intent(in) :: zoz0m     ! zref/roughness(momentum)                  [   ---]
      real   , intent(in) :: lnzoz0m   ! ln[zref/roughness(momentum)]              [   ---]
      real   , intent(in) :: zoz0h     ! zref/roughness(heat)                      [   ---]
      real   , intent(in) :: lnzoz0h   ! ln[zref/roughness(heat)]                  [   ---]
      logical, intent(in) :: stable    ! Flag... This surface layer is stable      [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: fm        ! lnzoz0m - psim(zeta) + psim(zeta0m)       [   ---]
      real                :: fh        ! lnzoz0h - psih(zeta) + psih(zeta0h)       [   ---]
      real                :: dfmdzeta  ! d(fm)/d(zeta)                             [   ---]
      real                :: dfhdzeta  ! d(fh)/d(zeta)                             [   ---]
      real                :: z0moz     ! Roughness(momentum) / Reference height    [   ---]
      real                :: zeta0m    ! Roughness(momentum) / Obukhov length      [   ---]
      real                :: z0hoz     ! Roughness(heat) / Reference height        [   ---]
      real                :: zeta0h    ! Roughness(heat) / Obukhov length          [   ---]
      real                :: zetaa     ! Smallest guess (or previous guess)        [   ---]
      real                :: zetaz     ! Largest guess (or new guess in Newton's)  [   ---]
      real                :: deriv     ! Function Derivative                       [   ---]
      real                :: fun       ! Function for which we seek a root.        [   ---]
      real                :: funa      ! Smallest guess function.                  [   ---]
      real                :: funz      ! Largest guess function.                   [   ---]
      real                :: delta     ! Aux. var --- 2nd guess for bisection      [   ---]
      real                :: zetamin   ! Minimum zeta for stable case.             [   ---]
      real                :: zetamax   ! Maximum zeta for unstable case.           [   ---]
      real                :: zetasmall ! Zeta dangerously close to zero            [   ---]
      integer             :: itb       ! Iteration counters                        [   ---]
      integer             :: itn       ! Iteration counters                        [   ---]
      integer             :: itp       ! Iteration counters                        [   ---]
      logical             :: converged ! Flag... The method converged!             [   T|F]
      logical             :: zside     ! Flag... I'm on the z-side.                [   T|F]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     First thing: if the bulk Richardson number is zero or almost zero, then we     !
      ! rather just assign z/L to be the one given by Oncley and Dudhia (1995).  This      !
      ! saves time and also avoids the risk of having zeta with the opposite sign.         !
      !------------------------------------------------------------------------------------!
      zetasmall = vkopr * rib * min(lnzoz0m,lnzoz0h)
      if (rib <= 0. .and. zetasmall > - z0moz0h * toler) then
         zoobukhov = vkopr * rib * lnzoz0m
         return
      elseif (rib > 0. .and. zetasmall < z0moz0h * toler) then
         zoobukhov = zetasmall / (1.1 - 5.0 * rib)
         return
      else
         zetamin    =  toler
         zetamax    = -toler
      end if

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(60a1)') ('-',itn=1,60)
      !write(unit=89,fmt='(5(a,1x,f11.4,1x),a,l1)')                                         &
      !   'Input values: Rib =',rib,'zref=',zref,'rough=',rough,'zoz0=',zoz0                &
      !           ,'lnzoz0=',lnzoz0,'stable=',stable
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!


      !----- Defining some values that won't change during the iterative method. ----------!
      z0moz = 1. / zoz0m
      z0hoz = 1. / zoz0h

      !------------------------------------------------------------------------------------!
      !     First guess, using Oncley and Dudhia (1995) approximation for unstable case.   !
      ! We won't use the stable case to avoid FPE or zeta with opposite sign when          !
      ! Ri > 0.20.                                                                         !
      !------------------------------------------------------------------------------------!
      zetaa = vkopr * rib * lnzoz0m

      !----- Finding the function and its derivative. -------------------------------------!
      zeta0m   = zetaa * z0moz
      zeta0h   = zetaa * z0hoz
      fm       = lnzoz0m - psim(zetaa,stable) + psim(zeta0m,stable)
      fh       = lnzoz0h - psih(zetaa,stable) + psih(zeta0h,stable)
      dfmdzeta = z0moz * dpsimdzeta(zeta0m,stable) - dpsimdzeta(zetaa,stable)
      dfhdzeta = z0hoz * dpsihdzeta(zeta0h,stable) - dpsihdzeta(zetaa,stable)
      funa     = vkopr * rib * fm * fm / fh - zetaa
      deriv    = vkopr * rib * (2. * fm * dfmdzeta * fh - fm * fm * dfhdzeta)              &
               / (fh * fh) - 1.

      !----- Copying just in case it fails at the first iteration. ------------------------!
      zetaz = zetaa
      fun   = funa

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
      !   '1STGSS: itn=',0,'bisection=',.false.,'zetaz=',zetaz,'fun=',fun,'fm=',fm          &
      !  ,'fh=',fh,'dfmdzeta=',dfmdzeta,'dfhdzeta=',dfhdzeta,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      !----- Enter Newton's method loop. --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1, maxfpo/6
         !---------------------------------------------------------------------------------!
         !     Newton's method converges fast when it's on the right track, but there are  !
         ! cases in which it becomes ill-behaved.  Two situations are known to cause       !
         ! trouble:                                                                        !
         ! 1.  If the derivative is tiny, the next guess can be too far from the actual    !
         !     answer;                                                                     !
         ! 2.  For this specific problem, when zeta is too close to zero.  In this case    !
         !     the derivative will tend to infinity at this point and Newton's method is   !
         !     not going to perform well and can potentially enter in a weird behaviour or !
         !     lead to the wrong answer.  In any case, so we rather go with bisection.     !
         !---------------------------------------------------------------------------------!
         if (abs(deriv) < toler) then
            exit newloop
         elseif(stable .and. (zetaz - fun/deriv < zetamin)) then
            exit newloop
         elseif((.not. stable) .and. (zetaz - fun/deriv > zetamax)) then
            exit newloop
         end if

         !----- Copying the previous guess ------------------------------------------------!
         zetaa = zetaz
         funa  = fun
         !----- New guess, its function and derivative evaluation -------------------------!
         zetaz = zetaa - fun/deriv

         zeta0m   = zetaz * z0moz
         zeta0h   = zetaz * z0hoz
         fm       = lnzoz0m - psim(zetaz,stable) + psim(zeta0m,stable)
         fh       = lnzoz0h - psih(zetaz,stable) + psih(zeta0h,stable)
         dfmdzeta = z0moz * dpsimdzeta(zeta0m,stable) - dpsimdzeta(zetaz,stable)
         dfhdzeta = z0hoz * dpsihdzeta(zeta0h,stable) - dpsihdzeta(zetaz,stable)
         fun      = vkopr * rib * fm * fm / fh - zetaz
         deriv    = vkopr * rib * (2. * fm * dfmdzeta * fh - fm * fm * dfhdzeta)           &
                  / (fh * fh) - 1.

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
         !   'NEWTON: itn=',itn,'bisection=',.false.,'zetaz=',zetaz,'fun=',fun,'fm=',fm        &
         !  ,'fh=',fh,'dfmdzeta=',dfmdzeta,'dfhdzeta=',dfhdzeta,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         converged = abs(zetaz-zetaa) < toler * abs(zetaz)

         if (converged) then
            zoobukhov = 0.5 * (zetaa+zetaz)
            return
         elseif (fun == 0.0) then !---- Converged by luck. --------------------------------!
            zoobukhov = zetaz
            return
         end if
      end do newloop

      !------------------------------------------------------------------------------------!
      !     If we reached this point then it's because Newton's method failed or it has    !
      ! become too dangerous.  We use the Regula Falsi (Illinois) method, which is just a  !
      ! fancier bisection.  For this we need two guesses, and the guesses must have        !
      ! opposite signs.                                                                    !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.0) then
         funz  = fun
         zside = .true. 
      else
         if (abs(fun-funa) < 100. * toler * abs(zetaa)) then
            if (stable) then
               delta = max(0.5 * abs(zetaa-zetamin),100. * toler * abs(zetaa))
            else
               delta = max(0.5 * abs(zetaa-zetamax),100. * toler * abs(zetaa))
            end if
         else
            if (stable) then
               delta = max(abs(funa * (zetaz-zetaa)/(fun-funa))                            &
                          ,100. * toler * abs(zetaa)                                       &
                          ,0.5 * abs(zetaa-zetamin))
            else
               delta = max(abs(funa * (zetaz-zetaa)/(fun-funa))                            &
                          ,100. * toler * abs(zetaa)                                       &
                          ,0.5 * abs(zetaa-zetamax))
            end if
         end if
         if (stable) then
            zetaz = max(zetamin,zetaa + delta)
         else
            zetaz = min(zetamax,zetaa + delta)
         end if
         zside = .false.
         zgssloop: do itp=1,maxfpo
            if (stable) then
               zetaz    = max(zetamin,zetaa + real((-1)**itp * (itp+3)/2) * delta)
            else
               zetaz    = min(zetamax,zetaa + real((-1)**itp * (itp+3)/2) * delta)
            end if
            zeta0m   = zetaz * z0moz
            zeta0h   = zetaz * z0hoz
            fm       = lnzoz0m - psim(zetaz,stable) + psim(zeta0m,stable)
            fh       = lnzoz0h - psih(zetaz,stable) + psih(zeta0h,stable)
            funz     = vkopr * rib * fm * fm / fh - zetaz
            zside    = funa * funz < 0.0
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                 &
            !   '2NDGSS: itp=',itp,'zside=',zside,'zetaa=',zetaa,'zetaz=',zetaz             &
            !  ,'funa=',funa,'funz=',funz,'fm=',fm,'fh=',fh,'delta=',delta
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') '=================================================='
            write (unit=*,fmt='(a)') '    No second guess for you...'
            write (unit=*,fmt='(a)') '=================================================='
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zref   =',zref   ,'rough  =',rough
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'lnzoz0m=',lnzoz0m,'lnzoz0h=',lnzoz0h
            write (unit=*,fmt='(1(a,1x,es14.7,1x))') 'rib    =',rib
            write (unit=*,fmt='(1(a,1x,l1,1x))')     'stable =',stable
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'fun    =',fun    ,'delta  =',delta
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zetaa  =',zetaa  ,'funa   =',funa
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zetaz  =',zetaz  ,'funz   =',funz
            call abort_run('Failed finding the second guess for regula falsi'              &
                            ,'zoobukhov','leaf_coms.f90')
         end if
      end if

      !----- Now we are ready to start the regula falsi method. ---------------------------!
      bisloop: do itb=itn,maxfpo
         zoobukhov = (funz*zetaa-funa*zetaz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(zoobukhov-zetaa) < toler * abs(zoobukhov)
         if (converged) exit bisloop

         !------ Finding the new function -------------------------------------------------!
         zeta0m   = zoobukhov * z0moz
         zeta0h   = zoobukhov * z0hoz
         fm       = lnzoz0m - psim(zoobukhov,stable) + psim(zeta0m,stable)
         fh       = lnzoz0h - psih(zoobukhov,stable) + psih(zeta0h,stable)
         fun      = vkopr * rib * fm * fm / fh - zoobukhov

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
         !   'REGULA: itn=',itb,'bisection=',.true.,'zetaa=',zetaa,'zetaz=',zetaz,'fun=',fun   &
         !  ,'funa=',funa,'funz=',funz,'fm=',fm,'fh=',fh
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

         !------ Defining my new interval based on the intermediate value theorem. --------!
         if (fun*funa < 0. ) then
            zetaz = zoobukhov
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa = funa * 0.5
            !----- We just updated zside, setting zside to true. --------------------------!
            zside = .true.
         else
            zetaa = zoobukhov
            funa  = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz = funz * 0.5
            !----- We just updated aside, setting aside to true. --------------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' Zeta finding didn''t converge!!!'
         write (unit=*,fmt='(a,1x,i5,1x,a)') ' I gave up, after',maxfpo,'iterations...'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Input values.'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rib             [   ---] =',rib
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zref            [     m] =',zref
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rough           [     m] =',rough
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoz0m           [   ---] =',zoz0m
         write (unit=*,fmt='(a,1x,f12.4)' ) 'lnzoz0m         [   ---] =',lnzoz0m
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoz0h           [   ---] =',zoz0h
         write (unit=*,fmt='(a,1x,f12.4)' ) 'lnzoz0h         [   ---] =',lnzoz0h
         write (unit=*,fmt='(a,1x,l1)'    ) 'stable          [   T|F] =',stable
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome (downdraft values).'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetaa           [   ---] =',zetaa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetaz           [   ---] =',zetaz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fun             [   ---] =',fun
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fm              [   ---] =',fm
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fh              [   ---] =',fh
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funa            [   ---] =',funa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funz            [   ---] =',funz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'deriv           [   ---] =',deriv
         write (unit=*,fmt='(a,1x,es12.4)') 'toler           [   ---] =',toler
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [   ---] ='                   &
                                                            ,abs(zetaz-zetaa)/abs(zetaz)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoobukhov       [   ---] =',zoobukhov
         write (unit=*,fmt='(a)') '-------------------------------------------------------'

         call abort_run('Zeta didn''t converge, giving up!!!'                              &
                         ,'zoobukhov','leaf_coms.f90')
      end if

      return
   end function zoobukhov
   !=======================================================================================!
   !=======================================================================================!
end module leaf_coms
!==========================================================================================!
!==========================================================================================!




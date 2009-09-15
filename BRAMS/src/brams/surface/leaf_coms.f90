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
   use rconstants, only: grav

   integer :: niter_leaf   & ! number of leaf timesteps in model long timestep
            , niter_can    ! ! number of canopy timesteps in leaf timestep

   real    :: dtll         & ! leaf timestep
            , dtll_factor  & ! leaf timestep factor (leaf timestep / model timestep)
            , dtlc         & ! canopy timestep
            , dtlc_factor  & ! canopy timestep factor (canopy timestep / leaf timestep)
            , dtllowcc     & ! leaf timestep   / (can_depth * can_rhos)
            , dtlcowcc     & ! canopy timestep / (can_depth * can_rhos)
            , dtlloccc     & ! mmdry * leaf timestep   / (can_depth * can_rhos)
            , dtlcoccc     & ! mmdry * canopy timestep / (can_depth * can_rhos)

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
            , gzotheta     & ! (g*z/theta) at top of surface layer (for Ri#)     [  J/kg/K]
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
            , veg_temp     & ! vegetation temperature                            [       K]
            , veg_fliq     & ! liquid fraction of vegetation surface water       [      --]
            , ground_temp  & ! ground temperature                                [       K]
            , ground_fliq  & ! ground water liquid fraction                      [     ---]
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
   real, dimension(nstyp)           :: ssand,sclay,sorgan,sporo,soilcp,slfc,emisg,slcons00
   real, dimension(nstyp)           :: slcons0,fhydraul
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

   !----- Canopy depth --------------------------------------------------------------------!
   real, parameter :: can_depth = 20.0 ! [metres]

   !---------------------------------------------------------------------------------------!
   !     Vegetation heat capacity.  A constant, but it could be scaled by LAI, height,     !
   ! etc...                                                                                !
   !---------------------------------------------------------------------------------------!
   real, parameter :: hcapveg = 3.e4  ! [J/m2/K]

   !----- Some constants to ensure the model good behaviour -------------------------------!
   real, parameter :: min_sfcwater_mass  = 1.e-6
   real, parameter :: min_sfcwater_depth = 1.e-9
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Speed-related minimum values we will consider.                                    !
   !---------------------------------------------------------------------------------------!
   real, parameter :: ubmin    = 0.25 ! Minimum velocity                         [     m/s]
   real, parameter :: ustmin   = 0.1  ! Minimum ustar                            [     m/s]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Constants for surface layer models.                                              !
   !---------------------------------------------------------------------------------------!
   !----- Louis (1979) model. -------------------------------------------------------------!
   real, parameter   :: bl79     = 5.0    ! b prime parameter
   real, parameter   :: csm      = 7.5    ! C* for momentum (eqn. 20, not co2 char. scale)
   real, parameter   :: csh      = 5.0    ! C* for heat (eqn.20, not co2 char. scale)
   real, parameter   :: dl79     = 5.0    ! ???
   !----- Oncley and Dudhia (1995) model. -------------------------------------------------!
   real, parameter   :: bbeta    = 5.0    ! Beta used by Businger et al. (1971)
   real, parameter   :: gamm     = 15.0   ! Gamma used by Businger et al. (1971) - momentum.
   real, parameter   :: gamh     = 9.0    ! Gamma used by Businger et al. (1971) - heat.
   real, parameter   :: rri      = 1./.74 ! 1/R value, used by Businger and Louis.
   real, parameter   :: ribmax   = 0.20   ! Maximum bulk Richardson number
   real, parameter   :: tprandtl = 1.00   ! Turbulent Prandtl number.
   !---------------------------------------------------------------------------------------!



   !----- Parameters that used to be in LEAF-3 (leaftw) -----------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Parameters for new soil heat conductivity (8/17/00):  Move to sfcdata later.       !
   !---------------------------------------------------------------------------------------!
   real, dimension(nstyp), parameter :: soilcond0 = (/ 0.30, 0.30, 0.29, 0.27, 0.28, 0.28  &
                                                     , 0.26, 0.27, 0.27, 0.25, 0.25, 0.06 /)
   real, dimension(nstyp), parameter :: soilcond1 = (/ 4.80, 4.66, 4.27, 3.47, 3.63, 3.78  &
                                                     , 2.73, 3.23, 3.32, 2.58, 2.40, 0.46 /)
   real, dimension(nstyp), parameter :: soilcond2 = (/-2.70,-2.60,-2.31,-1.74,-1.85,-1.96  &
                                                     ,-1.20,-1.56,-1.63,-1.09,-0.96, 0.00 /)
   !---------------------------------------------------------------------------------------!


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
end module leaf_coms
!==========================================================================================!
!==========================================================================================!




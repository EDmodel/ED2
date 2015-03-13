!==========================================================================================!
!==========================================================================================!
!    This module contains several parameters used by the Runge-Kutta integrator scheme.    !
! This module is intended to host the numerical method variables only, variables related   !
! to physical, hydrological or ecological properties should be stored somewhere else.      !
!                                                                                          !
!   MLO. I attempted to describe some variables, not sure about all of them, if the        !
!        definition is wrong, please correct it. Thanks!                                   !
!------------------------------------------------------------------------------------------!
module rk4_coms
   use ed_max_dims, only : n_pft   & ! intent(in)
                         , nzgmax  & ! intent(in)
                         , str_len ! ! intent(in)

   implicit none

   !=======================================================================================!
   !=======================================================================================!
   !     Structure used for the integration using the RK4 method.  Real variables are      !
   ! stored in double precision, so we can solve pools whose order of magnitude are very   !
   ! different.                                                                            !
   !---------------------------------------------------------------------------------------!
   type rk4patchtype

      !----- Canopy air variables. -------------------------------------------------------!
      real(kind=8)                        :: can_enthalpy ! Canopy sp. enthalpy  [    J/kg]
      real(kind=8)                        :: can_theta    ! Pot. Temperature     [       K]
      real(kind=8)                        :: can_temp     ! Temperature          [       K]
      real(kind=8)                        :: can_shv      ! Specific humidity    [   kg/kg]
      real(kind=8)                        :: can_ssh      ! Sat. spec. humidity  [   kg/kg]
      real(kind=8)                        :: can_rhv      ! Relative humidity    [     ---]
      real(kind=8)                        :: can_co2      ! CO_2                 [�mol/mol]
      real(kind=8)                        :: can_depth    ! Canopy depth         [       m]
      real(kind=8)                        :: can_rhos     ! Canopy air density   [   kg/m�]
      real(kind=8)                        :: can_prss     ! Pressure             [      Pa]
      real(kind=8)                        :: can_exner    ! Exner function       [  J/kg/K]
      real(kind=8)                        :: can_cp       ! Specific heat        [  J/kg/K]
      !-----------------------------------------------------------------------------------!

      !------ Above Canopy Air Variables -------------------------------------------------!
      real(kind=8)                        :: vels         ! wind speed [m/s]
      real(kind=8)                        :: atm_enthalpy ! Specific enthalpy of the air
                                                          ! above the canopy, dependent
                                                          ! on canopy geometry, so not a
                                                          ! site variable        [    J/kg]

      !----- Vegetation properties. -------------------------------------------------------!
      real(kind=8)                        :: veg_height   ! Vegetation height    [       m]
      real(kind=8)                        :: veg_displace ! 0-plane displacement [       m]
      real(kind=8)                        :: veg_rough    ! Vegetation roughness [       m]
      !------------------------------------------------------------------------------------!


      !----- Fraction of open canopy. -----------------------------------------------------!
      real(kind=8)                        :: opencan_frac   ! Frac. of open canopy
      !------------------------------------------------------------------------------------!



      !----- Total depth of temporary surface water / snow. -------------------------------!
      real(kind=8)                        :: total_sfcw_depth
      real(kind=8)                        :: total_sfcw_mass
      !------------------------------------------------------------------------------------!



      !----- Fraction of open canopy buried in snow. --------------------------------------!
      real(kind=8)                        :: snowfac
      !------------------------------------------------------------------------------------!



      !----- Ground -> Canopy flux type. --------------------------------------------------!
      real(kind=8)                        :: ggbare       ! Cond. of bare ground [     m/s]
      real(kind=8)                        :: ggveg        ! Cond. of veg. ground [     m/s]
      real(kind=8)                        :: ggnet        ! Net ground  conduct. [     m/s]
      real(kind=8)                        :: ggsoil       ! Soil evap. conduct.  [     m/s]
      integer                             :: flag_wflxgc  ! Flag for water flux.
      !------------------------------------------------------------------------------------!



      !----- Soil variables. --------------------------------------------------------------!
      real(kind=8), dimension(:), pointer :: soil_energy  ! Internal energy       [   J/m�]
      real(kind=8), dimension(:), pointer :: soil_mstpot  ! Soil matric potential [      m]
      real(kind=8), dimension(:), pointer :: soil_tempk   ! Specific humidity     [      K]
      real(kind=8), dimension(:), pointer :: soil_fracliq ! Liquid fraction       [   ----]
      real(kind=8), dimension(:), pointer :: soil_water   ! Water content         [  m�/m�]
      !------------------------------------------------------------------------------------!



      !----- Temporary surface water variables. -------------------------------------------!
      integer                             :: nlev_sfcwater    ! # of layers       [   ----]
      integer                             :: flag_sfcwater    ! Status flag       [  -----]
      real(kind=8), dimension(:), pointer :: sfcwater_depth   ! Depth             [      m]
      real(kind=8), dimension(:), pointer :: sfcwater_mass    ! Mass              [  kg/m�]
      real(kind=8), dimension(:), pointer :: sfcwater_energy  ! Internal energy   [   J/m�]
      real(kind=8), dimension(:), pointer :: sfcwater_tempk   ! Temperature       [      K]
      real(kind=8), dimension(:), pointer :: sfcwater_fracliq ! Liquid fraction   [   ----]
      !------------------------------------------------------------------------------------!



      !----- Virtual layer variables. -----------------------------------------------------!
      real(kind=8)                        :: virtual_water    ! Mass              [  kg/m�]
      real(kind=8)                        :: virtual_energy   ! Internal energy   [  kg/m�]
      real(kind=8)                        :: virtual_depth    ! Depth             [      m]
      real(kind=8)                        :: virtual_tempk    ! Temperature       [      K]
      real(kind=8)                        :: virtual_fracliq  ! Liquid fraction   [      K]
      !------------------------------------------------------------------------------------!



      !----- Surface variables. -----------------------------------------------------------!
      real(kind=8)                        :: ground_shv   ! Ground sp. humidity   [  kg/kg]
      real(kind=8)                        :: ground_ssh   ! Ground sat. humidity  [  kg/kg]
      real(kind=8)                        :: ground_temp  ! Ground temperature    [      K]
      real(kind=8)                        :: ground_fliq  ! Ground liquid frac.   [    ---]
      real(kind=8)                        :: rough        ! Roughness             [      m]
      !------------------------------------------------------------------------------------!



      !----- Characteristic scale. --------------------------------------------------------!
      real(kind=8)                        :: ustar  ! Momentum                    [    m/s]
      real(kind=8)                        :: cstar  ! Carbon mixing ratio         [�mol/m�]
      real(kind=8)                        :: tstar  ! Temperature                 [      K]
      real(kind=8)                        :: qstar  ! Water vapour spec. humidity [  kg/kg]
      real(kind=8)                        :: estar  ! Eq. potential temperature   [      K]
      real(kind=8)                        :: zeta   ! z / Obukhov length          [    ---]
      real(kind=8)                        :: ribulk ! Bulk Richardson number      [    ---]
      !------------------------------------------------------------------------------------!



      !----- Vertical fluxes. -------------------------------------------------------------!
      real(kind=8)                        :: upwp 
      real(kind=8)                        :: qpwp
      real(kind=8)                        :: cpwp
      real(kind=8)                        :: tpwp
      real(kind=8)                        :: wpwp
      !------------------------------------------------------------------------------------!



      !----- Resistance variables. --------------------------------------------------------!
      real(kind=8)                        :: rasveg
      real(kind=8)                        :: root_res_fac
      !------------------------------------------------------------------------------------!



      !----- Heterotrophic respiration.[�mol/m�/s] ----------------------------------------!
      real(kind=8)                        :: cwd_rh
      real(kind=8)                        :: rh
      !------------------------------------------------------------------------------------!



      !----- Leaf (cohort-level) variables. -----------------------------------------------!
      real(kind=8), pointer, dimension(:) :: leaf_energy     ! Internal energy  [     J/m�]
      real(kind=8), pointer, dimension(:) :: leaf_water      ! Sfc. water mass  [    kg/m�]
      real(kind=8), pointer, dimension(:) :: leaf_temp       ! Temperature      [        K]
      real(kind=8), pointer, dimension(:) :: leaf_fliq       ! Liquid fraction  [      ---]
      real(kind=8), pointer, dimension(:) :: leaf_hcap       ! Heat capacity    [   J/m�/K]
      real(kind=8), pointer, dimension(:) :: leaf_reynolds   ! Reynolds number  [      ---]
      real(kind=8), pointer, dimension(:) :: leaf_grashof    ! Grashof number   [      ---]
      real(kind=8), pointer, dimension(:) :: leaf_nussfree   ! Nusselt # (free) [      ---]
      real(kind=8), pointer, dimension(:) :: leaf_nussforc   ! Nusselt # (forc.)[      ---]
      real(kind=8), pointer, dimension(:) :: lint_shv        ! Interc. sp. hum. [    kg/kg]
      logical     , pointer, dimension(:) :: leaf_resolvable ! resolve leaves?  [      T|F]
      real(kind=8), pointer, dimension(:) :: leaf_gbh        ! Bnd.lyr. condct. [ J/K/m�/s]
      real(kind=8), pointer, dimension(:) :: leaf_gbw        ! Bnd.lyr. condct. [  kg/m�/s]
      real(kind=8), pointer, dimension(:) :: gsw_open        ! Sto. cond. (op.) [ J/K/m�/s]
      real(kind=8), pointer, dimension(:) :: gsw_closed      ! Sto. cond. (cl.) [  kg/m�/s]
      real(kind=8), pointer, dimension(:) :: rshort_l        ! Absorbed SWRad.  [   J/m�/s]
      real(kind=8), pointer, dimension(:) :: rlong_l         ! Absorbed LWRad.  [   J/m�/s]
      !------------------------------------------------------------------------------------!


      !----- Wood (cohort-level) variables. -----------------------------------------------!
      real(kind=8), pointer, dimension(:) :: wood_energy     ! Internal energy  [     J/m�]
      real(kind=8), pointer, dimension(:) :: wood_water      ! Sfc. water mass  [    kg/m�]
      real(kind=8), pointer, dimension(:) :: wood_temp       ! Temperature      [        K]
      real(kind=8), pointer, dimension(:) :: wood_fliq       ! Liquid fraction  [      ---]
      real(kind=8), pointer, dimension(:) :: wood_hcap       ! Heat capacity    [   J/m�/K]
      real(kind=8), pointer, dimension(:) :: wood_reynolds   ! Reynolds number  [      ---]
      real(kind=8), pointer, dimension(:) :: wood_grashof    ! Grashof number   [      ---]
      real(kind=8), pointer, dimension(:) :: wood_nussfree   ! Nusselt # (free) [      ---]
      real(kind=8), pointer, dimension(:) :: wood_nussforc   ! Nusselt # (forc.)[      ---]
      logical     , pointer, dimension(:) :: wood_resolvable ! resolve wood?    [      T|F]
      real(kind=8), pointer, dimension(:) :: wood_gbh        ! Bnd.lyr. condct. [ J/K/m�/s]
      real(kind=8), pointer, dimension(:) :: wood_gbw        ! Bnd.lyr. condct. [  kg/m�/s]
      real(kind=8), pointer, dimension(:) :: rshort_w        ! Absorbed SWRad.  [   J/m�/s]
      real(kind=8), pointer, dimension(:) :: rlong_w         ! Absorbed LWRad.  [   J/m�/s]
      !------------------------------------------------------------------------------------!



      !----- Leaf+branchwood (cohort-level) variables. ------------------------------------!
      real(kind=8), pointer, dimension(:) :: veg_energy      ! Internal energy  [     J/m�]
      real(kind=8), pointer, dimension(:) :: veg_water       ! Sfc. water mass  [    kg/m�]
      real(kind=8), pointer, dimension(:) :: veg_hcap        ! Heat capacity    [   J/m�/K]
      logical     , pointer, dimension(:) :: veg_resolvable  ! resolve leaves?  [      T|F]
      !------------------------------------------------------------------------------------!



      !----- General cohort-level properties. ---------------------------------------------!
      real(kind=8), pointer, dimension(:) :: nplant       ! Plant density       [ plant/m�]
      real(kind=8), pointer, dimension(:) :: veg_wind      ! Cohort-level wind  [      m/s]
      real(kind=8), pointer, dimension(:) :: lai          ! Leaf area index     [    m�/m�]
      real(kind=8), pointer, dimension(:) :: wai          ! Wood area index     [    m�/m�]
      real(kind=8), pointer, dimension(:) :: tai          ! Tree area index     [    m�/m�]
      real(kind=8), pointer, dimension(:) :: crown_area   ! Crown area          [    m�/m�]
      real(kind=8), pointer, dimension(:) :: elongf       ! Elongation factor   [     ----]
      real(kind=8), pointer, dimension(:) :: psi_open     ! Water demand (op.)  [kg/m�lf/s]
      real(kind=8), pointer, dimension(:) :: psi_closed   ! Water demand (clos.)[kg/m�lf/s]
      real(kind=8), pointer, dimension(:) :: fs_open      ! Frac. of op. stom.  [      ---]
      real(kind=8), pointer, dimension(:) :: gpp          ! Gross primary prod. [�mol/m�/s]
      real(kind=8), pointer, dimension(:) :: leaf_resp    ! Leaf respiration    [�mol/m�/s]
      real(kind=8), pointer, dimension(:) :: root_resp    ! Root respiration    [�mol/m�/s]
      real(kind=8), pointer, dimension(:) :: growth_resp  ! Growth respiration  [�mol/m�/s]
      real(kind=8), pointer, dimension(:) :: storage_resp ! Storage respiration [�mol/m�/s]
      real(kind=8), pointer, dimension(:) :: vleaf_resp   ! Virtual leaf resp.  [�mol/m�/s]


      !------ Variables used for hybrid stepping -----------------------------------------!
      real(kind=8)                        :: hflxgc
      real(kind=8)                        :: hflxsc
      real(kind=8)                        :: wflxgc
      real(kind=8)                        :: wflxsc
      real(kind=8)                        :: wflxac
      real(kind=8)                        :: qwflxgc
      real(kind=8)                        :: qwflxsc
      real(kind=8), pointer, dimension(:) :: wflxtr
      real(kind=8), pointer, dimension(:) :: wflxlc
      real(kind=8), pointer, dimension(:) :: wflxwc
      real(kind=8), pointer, dimension(:) :: hflx_lrsti ! heat gained from rnet,shed,tr,int
      real(kind=8), pointer, dimension(:) :: hflx_wrsti !
      !-----------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Fast time flux diagnostic variables.  These variables may be turned off under  !
      ! different conditions.                                                              !
      !------------------------------------------------------------------------------------!
      real(kind=8) :: avg_rshort_gnd    ! Total absorbed SW radiation  (ground)
      real(kind=8) :: avg_par_gnd       ! Total absorbed PAR radiation (ground)
      real(kind=8) :: avg_rlong_gnd     ! Net absorbed LW radiation    (ground)
      !----- Water fluxes -----------------------------------------------------------------!
      real(kind=8) :: avg_vapor_gc      ! Ground     -> canopy air:  evaporation flux
      real(kind=8) :: avg_vapor_ac      ! Free atm.  -> canopy air:  vapour flux
      !----- Mass and energy input due to throughfall precipitation. ----------------------!
      real(kind=8) :: avg_throughfall   ! Mass flux
      real(kind=8) :: avg_qthroughfall  ! Energy flux
      !----- Sensible heat flux -----------------------------------------------------------!
      real(kind=8) :: avg_sensible_gc   ! Ground    -> canopy air
      real(kind=8) :: avg_sensible_ac   ! Free atm. -> canopy air
      real(kind=8) :: avg_heatstor_veg  ! Heat storage in vegetation
      !----- Carbon flux ------------------------------------------------------------------!
      real(kind=8) :: avg_carbon_ac     ! Free atm. -> canopy air
      real(kind=8) :: avg_carbon_st     ! Canopy air storage flux
      !----- Characteristic gradient scales (stars). --------------------------------------!
      real(kind=8) :: avg_ustar         ! Friction velocity
      real(kind=8) :: avg_tstar         ! Temperature
      real(kind=8) :: avg_qstar         ! Specific humidity
      real(kind=8) :: avg_cstar         ! CO2 mixing ratio
      !----- Soil fluxes ------------------------------------------------------------------!
      real(kind=8),pointer,dimension(:) :: avg_smoist_gg     ! Moisture flux between layers
      real(kind=8),pointer,dimension(:) :: avg_transloss     ! Transpired soil moisture sink
      real(kind=8),pointer,dimension(:) :: avg_sensible_gg   ! Soil heat flux between layers
      real(kind=8)                      :: avg_drainage      ! Drainage at the bottom.
      real(kind=8)                      :: avg_qdrainage     ! Drainage at the bottom.
      !----- Cohort-level fluxes. ---------------------------------------------------------!
      real(kind=8),pointer,dimension(:) :: avg_sensible_lc    ! Sensible heat (Leaf-CAS)
      real(kind=8),pointer,dimension(:) :: avg_sensible_wc    ! Sensible heat (Wood-CAS)
      real(kind=8),pointer,dimension(:) :: avg_vapor_lc       ! Water flux (Leaf sfc - CAS)
      real(kind=8),pointer,dimension(:) :: avg_vapor_wc       ! Water flux (Wood sfc - CAS)
      real(kind=8),pointer,dimension(:) :: avg_transp         ! Transpiration (Leaf - CAS)
      real(kind=8),pointer,dimension(:) :: avg_intercepted_al ! Leaf interception
      real(kind=8),pointer,dimension(:) :: avg_intercepted_aw ! Wood interception
      real(kind=8),pointer,dimension(:) :: avg_wshed_lg       ! Leaf shedding
      real(kind=8),pointer,dimension(:) :: avg_wshed_wg       ! Wood shedding
      !----- Water deficit. ---------------------------------------------------------------!
      real(kind=8) :: water_deficit     ! Step water deficit
      !------------------------------------------------------------------------------------!
      !     Fast time flux variables for each time step.  These variables will be defined  !
      ! only when the user is debugging.                                                   !
      !------------------------------------------------------------------------------------!
      real(kind=8) :: flx_rshort_gnd    ! Absorbed SW radiation
      real(kind=8) :: flx_par_gnd       ! Absorbed PAR
      real(kind=8) :: flx_rlong_gnd     ! Absorbed LW radiation
      !----- Water fluxes -----------------------------------------------------------------!
      real(kind=8) :: flx_vapor_lc      ! Leaf       -> canopy air:  evap./cond. flux
      real(kind=8) :: flx_vapor_wc      ! Wood       -> canopy air:  evap./cond. flux
      real(kind=8) :: flx_vapor_gc      ! Ground     -> canopy air:  evaporation flux
      real(kind=8) :: flx_vapor_ac      ! Free atm.  -> canopy air:  vapour flux
      real(kind=8) :: flx_transp        ! Transpiration
      !----- Mass and energy fluxes due to water shedding. --------------------------------!
      real(kind=8) :: flx_wshed_vg      ! Mass flux
      real(kind=8) :: flx_qwshed_vg     ! Energy flux
      !----- Mass and energy input due to intercepted water. ------------------------------!
      real(kind=8) :: flx_intercepted   ! Mass flux
      real(kind=8) :: flx_qintercepted  ! Energy flux
      !----- Mass and energy input due to throughfall precipitation. ----------------------!
      real(kind=8) :: flx_throughfall   ! Mass flux
      real(kind=8) :: flx_qthroughfall  ! Energy flux
      !----- Sensible heat flux -----------------------------------------------------------!
      real(kind=8) :: flx_sensible_lc   ! Leaf      -> canopy air
      real(kind=8) :: flx_sensible_wc   ! Wood      -> canopy air
      real(kind=8) :: flx_sensible_gc   ! Ground    -> canopy air
      real(kind=8) :: flx_sensible_ac   ! Free atm. -> canopy air
      real(kind=8) :: flx_heatstor_veg  ! Heat storage in vegetation
      !----- Carbon flux ------------------------------------------------------------------!
      real(kind=8) :: flx_carbon_ac     ! Free atm. -> canopy air
      real(kind=8) :: flx_carbon_st     ! Canopy CO2 storage flux
      !----- Soil fluxes ------------------------------------------------------------------!
      real(kind=8),pointer,dimension(:) :: flx_smoist_gg     ! Moisture flux between layers
      real(kind=8),pointer,dimension(:) :: flx_transloss     ! Transpired soil moisture sink
      real(kind=8),pointer,dimension(:) :: flx_sensible_gg   ! Soil heat flux between layers
      real(kind=8)                      :: flx_drainage      ! Drainage at the bottom.
      real(kind=8)                      :: flx_qdrainage     ! Drainage at the bottom.
      !----- Cohort-level fluxes. ---------------------------------------------------------!
      real(kind=8),pointer,dimension(:) :: cfx_hflxlc        ! Sensible heat
      real(kind=8),pointer,dimension(:) :: cfx_hflxwc        ! Sensible heat
      real(kind=8),pointer,dimension(:) :: cfx_qwflxlc       ! Latent heat - Evaporation
      real(kind=8),pointer,dimension(:) :: cfx_qwflxwc       ! Latent heat - Evaporation
      real(kind=8),pointer,dimension(:) :: cfx_qwshed        ! Int. en. of shed water
      real(kind=8),pointer,dimension(:) :: cfx_qtransp       ! Latent heat - Transpiration
      real(kind=8),pointer,dimension(:) :: cfx_qintercepted  ! Int. en. of intercept. H2O
      !----- Full budget variables --------------------------------------------------------!
      real(kind=8) :: co2budget_storage
      real(kind=8) :: co2budget_loss2atm
      real(kind=8) :: ebudget_storage
      real(kind=8) :: ebudget_netrad
      real(kind=8) :: ebudget_loss2atm
      real(kind=8) :: ebudget_loss2drainage
      real(kind=8) :: ebudget_loss2runoff
      real(kind=8) :: wbudget_storage
      real(kind=8) :: wbudget_loss2atm
      real(kind=8) :: wbudget_loss2drainage
      real(kind=8) :: wbudget_loss2runoff
   end type rk4patchtype
   !---------------------------------------------------------------------------------------!

   
   type bdf2patchtype
      real(kind=8)                        :: can_temp
      real(kind=8),pointer,dimension(:)   :: leaf_temp
      real(kind=8),pointer,dimension(:)   :: wood_temp
   end type bdf2patchtype


   !---------------------------------------------------------------------------------------!
   !    Structure with atmospheric and some other site-level data that is often used       !
   ! inside the Runge-Kutta integrator and do not change over time.                        !
   !---------------------------------------------------------------------------------------!
   type rk4sitetype
      integer                         :: lsl
      integer     , dimension(nzgmax) :: ntext_soil
      real(kind=8), dimension(n_pft)  :: green_leaf_factor
      real(kind=8)                    :: atm_rhos
!!      real(kind=8)                    :: vels
      real(kind=8)                    :: atm_ustar
!!      real(kind=8)                    :: atm_enthalpy  ! this is actually
                                                         ! affected by canopy geometry so
                                                         ! is really patch level
!!      real(kind=8)                    :: atm_tmp_zcan  ! Same, patch level, not used either
      real(kind=8)                    :: atm_tmp
      real(kind=8)                    :: atm_theta
      real(kind=8)                    :: atm_theiv
      real(kind=8)                    :: atm_vpdef
      real(kind=8)                    :: atm_shv
      real(kind=8)                    :: atm_rhv
      real(kind=8)                    :: atm_co2
      real(kind=8)                    :: zoff
      real(kind=8)                    :: atm_exner
      real(kind=8)                    :: pcpg
      real(kind=8)                    :: qpcpg
      real(kind=8)                    :: dpcpg
      real(kind=8)                    :: atm_prss
      real(kind=8)                    :: geoht
      real(kind=8)                    :: rshort
      real(kind=8)                    :: rlong
      real(kind=8)                    :: par_beam
      real(kind=8)                    :: par_diffuse
      real(kind=8)                    :: nir_beam
      real(kind=8)                    :: nir_diffuse
      real(kind=8)                    :: lon
      real(kind=8)                    :: lat
      real(kind=8)                    :: cosz
   end type rk4sitetype
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Structure with auxiliary gridded variables that will be used by the derivatives.   !
   !---------------------------------------------------------------------------------------!
   type rk4auxtype
      !----- Total potential [m]. ---------------------------------------------------------!
      real(kind=8), dimension(:)  , pointer :: psiplusz
      !----- Hydraulic conductivity [m/s]. ------------------------------------------------!
      real(kind=8), dimension(:)  , pointer :: hydcond
      !----- Available water factor at this layer [n/d]. ----------------------------------!
      real(kind=8), dimension(:)  , pointer :: avail_h2o_lyr
      !----- Integral of available water factor from top to this layer [n/d]. -------------!
      real(kind=8), dimension(:)  , pointer :: avail_h2o_int
      !----- Extracted water by transpiration [kg/m�]. ------------------------------------!
      real(kind=8), dimension(:,:), pointer :: extracted_water
      !----- Thermal conductivity [W/m/K]. ------------------------------------------------!
      real(kind=8), dimension(:)  , pointer :: th_cond_s
      real(kind=8), dimension(:)  , pointer :: th_cond_p
      !----- Sensible heat flux at staggered layer (k = k-1/2) [W/m�]. --------------------!
      real(kind=8), dimension(:)  , pointer :: h_flux_g
      real(kind=8), dimension(:)  , pointer :: h_flux_s
      !----- Water flux at staggered layers (k = k-1/2) [kg/m�/s]. ------------------------!
      real(kind=8), dimension(:)  , pointer :: w_flux_g
      !----- Latent heat flux at staggered layers (k=k-1/2) [W/m�]. -----------------------!
      real(kind=8), dimension(:)  , pointer :: qw_flux_g
      !----- Tests to check if the soil is too dry or too wet. ----------------------------!
      logical     , dimension(:)  , pointer :: drysoil
      logical     , dimension(:)  , pointer :: satsoil

      !----- Flag to tell whether there is at least one "resolvable" cohort in this patch -!
      logical :: any_resolvable
      !------------------------------------------------------------------------------------!
      
      !-------------------------------------------------------------------------------------!
      !      Canopy air space capacities.  These variables are used to convert the intensive!
      ! version of canopy air space prognostic variables (specific enthalpy, water vapour   !
      ! specific humidity and CO2 mixing ratio) into extensive variables.                   ! 
      !-------------------------------------------------------------------------------------!
      real(kind=8) :: wcapcan  ! Water capacity                             [  kg_air/m�gnd]
      real(kind=8) :: hcapcan  ! Enthalpy capacity                          [  kg_air/m�gnd]
      real(kind=8) :: ccapcan  ! CO2 capacity                               [ mol_air/m�gnd]
      real(kind=8) :: wcapcani ! Inverse of water capacity                  [  m�gnd/kg_air]
      real(kind=8) :: hcapcani ! Inverse of enthalpy capacity               [  m�gnd/kg_air]
      real(kind=8) :: ccapcani ! Inverse of CO2 capacity                    [ m�gnd/mol_air]


      
      !----- Minimums and Maximum Bounds on integration that change with each step ------- !
      real(kind=8)                          :: rk4min_can_prss
      real(kind=8)                          :: rk4max_can_prss
      real(kind=8)                          :: rk4min_can_theta
      real(kind=8)                          :: rk4max_can_theta
      real(kind=8)                          :: rk4min_can_enthalpy
      real(kind=8)                          :: rk4max_can_enthalpy

      real(kind=8), dimension(:),  pointer :: rk4min_soil_water
      real(kind=8), dimension(:),  pointer :: rk4max_soil_water

   end type rk4auxtype
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Structure with all the necessary buffers.                                             !
   !---------------------------------------------------------------------------------------!
   type integration_vars
      type(rk4patchtype), pointer :: initp   ! The current state
      type(rk4patchtype), pointer :: dinitp  ! The derivatives
      type(rk4patchtype), pointer :: yscal   ! The scale for prognostic variables
      type(rk4patchtype), pointer :: y       ! 
      type(rk4patchtype), pointer :: dydx    ! 
      type(rk4patchtype), pointer :: yerr    ! The error for the current guess
      type(rk4patchtype), pointer :: ytemp   ! Temporary
      type(rk4patchtype), pointer :: ak2     ! 
      type(rk4patchtype), pointer :: ak3     !
      type(rk4patchtype), pointer :: ak4     ! 
      type(rk4patchtype), pointer :: ak5     ! 
      type(rk4patchtype), pointer :: ak6     ! 
      type(rk4patchtype), pointer :: ak7     ! 
      type(bdf2patchtype), pointer :: yprev  ! Previous state
   end type integration_vars
   !---------------------------------------------------------------------------------------!


   !----- This is the actual integration buffer structure. --------------------------------!

   type(integration_vars), dimension(:), pointer :: integration_buff(:)
   type(rk4sitetype)                             :: rk4site
   type(rk4auxtype),dimension(:), pointer        :: rk4aux(:)

   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !    The following variable will be loaded from the user's namelist.                    !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   ! IBRANCH_THERMO -- This determines whether branches should be included in the          !
   !                   0.  No branches in energy/radiation (ED-2.1 default);               !
   !                   1.  Branches are accounted in the energy and radiation.  Branchwood !
   !                       and leaf are treated separately in the canopy radiation scheme, !
   !                       but solved as a single pool in the biophysics integration.      !
   !                   2.  Similar to 1, but branches are treated as separate pools in the !
   !                       biophysics (thus doubling the number of prognostic variables).  !
   !---------------------------------------------------------------------------------------!
   integer                   :: ibranch_thermo 
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! RK4_TOLERANCE -- This is the relative tolerance for Runge-Kutta or Heun's             !
   !                  integration.  Larger numbers will make runs go faster, at the        !
   !                  expense of being less accurate.  Currently the valid range is        !
   !                  between 1.e-7 and 1.e-1, but recommended values are between 1.e-4    !
   !                  and 1.e-2.                                                           !
   !---------------------------------------------------------------------------------------!
   real                      :: rk4_tolerance
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! IPERCOL -- This controls percolation and infiltration.                                !
   !            0.  Default method.  Assumes soil conductivity constant and for the        !
   !                temporary surface water, it sheds liquid in excess of a 1:9 liquid-    !
   !                -to-ice ratio through percolation.  Temporary surface water exists     !
   !                only if the top soil layer is at saturation.                           !
   !            1.  Constant soil conductivity, and it uses the percolation model as in    !
   !                Anderson (1976) NOAA technical report NWS 19.  Temporary surface       !
   !                water may exist after a heavy rain event, even if the soil doesn't     !
   !                saturate.  Recommended value.                                          !
   !            2.  Soil conductivity decreases with depth even for constant soil moisture !
   !                , otherwise it is the same as 1.                                       !
   !---------------------------------------------------------------------------------------!
   integer                   :: ipercol
   !---------------------------------------------------------------------------------------!

   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    The following variables are parameters that do not depend on the specific run.     !
   !---------------------------------------------------------------------------------------!

   !----- Small number, to avoid singularities. -------------------------------------------!
   real(kind=8), parameter :: tiny_offset = 1.d-20
   !----- Huge number, to bypass errmax checks. -------------------------------------------!
   real(kind=8), parameter :: huge_offset = 1.d30
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     These are all RK4 integrator factors.  These should never be changed, so they are !
   ! declared as parameters and will not be initialised at init_rk4_coms (ed_params.f90)   !
   !---------------------------------------------------------------------------------------!
   real(kind=8), parameter :: rk4_a2  = 2.d-1
   real(kind=8), parameter :: rk4_a3  = 3.d-1
   real(kind=8), parameter :: rk4_a4  = 6.d-1
   real(kind=8), parameter :: rk4_a5  = 1.d0
   real(kind=8), parameter :: rk4_a6  = 8.75d-1
   real(kind=8), parameter :: rk4_b21 = 2.d-1
   real(kind=8), parameter :: rk4_b31 = 3.d0/4.d1
   real(kind=8), parameter :: rk4_b32 = 9.d0/4.d1
   real(kind=8), parameter :: rk4_b41 = 3.d-1
   real(kind=8), parameter :: rk4_b42 = -9.d-1
   real(kind=8), parameter :: rk4_b43 = 1.2d0  
   real(kind=8), parameter :: rk4_b51 = -1.1d1/5.4d1
   real(kind=8), parameter :: rk4_b52 = 2.5d0
   real(kind=8), parameter :: rk4_b53 = -7.d1/2.7d1
   real(kind=8), parameter :: rk4_b54 = 3.5d1/2.7d1
   real(kind=8), parameter :: rk4_b61 = 1.631d3/5.52960d4
   real(kind=8), parameter :: rk4_b62 = 1.750d2/5.120d2
   real(kind=8), parameter :: rk4_b63 = 5.750d2/1.38240d4
   real(kind=8), parameter :: rk4_b64 = 4.42750d4/1.105920d5
   real(kind=8), parameter :: rk4_b65 = 2.530d2/4.0960d3
   real(kind=8), parameter :: rk4_c1  = 3.70d1/3.780d2
   real(kind=8), parameter :: rk4_c3  = 2.500d2/6.210d2
   real(kind=8), parameter :: rk4_c4  = 1.250d2/5.940d2
   real(kind=8), parameter :: rk4_c6  = 5.120d2/1.7710d3
   real(kind=8), parameter :: rk4_dc5 = -2.770d2/1.43360d4
   real(kind=8), parameter :: rk4_dc1 = rk4_c1-2.8250d3/2.76480d4
   real(kind=8), parameter :: rk4_dc3 = rk4_c3-1.85750d4/4.83840d4
   real(kind=8), parameter :: rk4_dc4 = rk4_c4-1.35250d4/5.52960d4
   real(kind=8), parameter :: rk4_dc6 = rk4_c6-2.5d-1
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     These are all Heun integrator factors.  These should never be changed, so they    !
   ! are declared as parameters and will not be initialised at init_rk4_coms               !
   ! (ed_params.f90)                                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=8), parameter :: heun_a2  = 1.d0
   real(kind=8), parameter :: heun_b21 = 1.d0
   real(kind=8), parameter :: heun_c1  = 5.d-1
   real(kind=8), parameter :: heun_c2  = 5.d-1
   real(kind=8), parameter :: heun_dc1 = heun_c1 - 1.d0
   real(kind=8), parameter :: heun_dc2 = heun_c2 - 0.d0
   !---------------------------------------------------------------------------------------!
   !     Maybe these could go to ed_params.f90 initialization.  Leaving them here for the  !
   ! time being.                                                                           !
   !---------------------------------------------------------------------------------------!
   real(kind=8), parameter :: safety =  9.d-1
   real(kind=8), parameter :: pgrow =  -2.d-1
   real(kind=8), parameter :: pshrnk = -2.5d-1
   real(kind=8), parameter :: errcon =  1.89d-4
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     The following variables are assigned at init_rk4_params (ed_params.f90).          !
   !---------------------------------------------------------------------------------------!

   !----- The following variables control the Runge-Kutta performance. --------------------!
   integer      :: maxstp         ! Maximum number of intermediate steps.
   real(kind=8) :: rk4eps         ! Desired relative accuracy
   real(kind=8) :: rk4epsi        ! Inverse of rk4eps
   real(kind=8) :: rk4eps2        ! rk4eps * rk4eps
   real(kind=8) :: hmin           ! Minimum step size
   logical      :: print_diags    ! Flag to print the diagnostic check.
   logical      :: checkbudget    ! Flag to decide whether we will check whether the 
                                  !    budgets close every time step (and stop the run if 
                                  !    they don't) or if we will skip this part.
   logical      :: record_err     ! Flag to keep track of which variable is causing the 
                                  !    most errors in the integrator.
   logical      :: print_detailed ! Flag to determine whether to print the patch 
                                  !     thermodynamic state every time step, including the
                                  !     intermediate ones.  This will create one file for
                                  !     each patch, so it is not recommended to be used 
                                  !     accross months, as the patches change.
   logical      :: print_thbnd    ! Flag to keep track of which variable is causing the
                                  !     most errors in the integrator.
   logical      :: print_budget   ! Flag that tells whether the budget is to be printed
                                  !     to a file.
   !---------------------------------------------------------------------------------------!



   !----- Constants used in rk4_derivs ----------------------------------------------------!
   logical      :: supersat_ok    ! It is fine for evaporation and transpiration to
                                  !    occur even if this causes the canopy air to 
                                  !    be super-saturated                          [   T|F]
                                  !    (N.B. Super-saturation can occur even if 
                                  !     supersat_ok is .false., but in this case
                                  !     only mixing with free atmosphere can cause
                                  !     the super-saturation).
   logical      :: leaf_intercept ! This flag is to turn on and on the leaf interception.  
                                  !    Except for developer tests, this variable should be 
                                  !    always true.  
   logical      :: debug          ! Verbose output for debug                       [   T|F]
   real(kind=8) :: toocold        ! Minimum temperature for saturation spec. hum.  [     K]
   real(kind=8) :: toohot         ! Maximum temperature for saturation spec. hum.  [     K]
   real(kind=8) :: lai_to_cover   ! Canopies with LAI less than this number are assumed to
                                  !     be open, ie, some fraction of the rain-drops can 
                                  !     reach the soil/litter layer unimpeded.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    The following variables are the bounds of what we could consider barely reasonable !
   ! during the integration process. If values go beyond these limits, reject the step     !
   ! right away.  These numbers are filled in init_rk4_params (ed_params.f90).  Most of    !
   ! them are in the same units we are used to, but note that the last ones have different !
   ! units to make it scalable.                                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=8) :: rk4min_can_temp      ! Minimum canopy    temperature         [        K]
   real(kind=8) :: rk4max_can_temp      ! Maximum canopy    temperature         [        K]
   real(kind=8) :: rk4min_can_shv       ! Minimum canopy    specific humidity   [kg/kg_air]
   real(kind=8) :: rk4max_can_shv       ! Maximum canopy    specific humidity   [kg/kg_air]
   real(kind=8) :: rk4min_can_rhv       ! Minimum canopy    relative humidity   [      ---]
   real(kind=8) :: rk4max_can_rhv       ! Maximum canopy    relative humidity   [      ---]
   real(kind=8) :: rk4min_can_co2       ! Minimum canopy    CO2 mixing ratio    [ �mol/mol]
   real(kind=8) :: rk4max_can_co2       ! Maximum canopy    CO2 mixing ratio    [ �mol/mol]
   real(kind=8) :: rk4min_soil_temp     ! Minimum soil      temperature         [        K]
   real(kind=8) :: rk4max_soil_temp     ! Maximum soil      temperature         [        K]
   real(kind=8) :: rk4min_veg_temp      ! Minimum leaf      temperature         [        K]
   real(kind=8) :: rk4max_veg_temp      ! Maximum leaf      temperature         [        K]
   real(kind=8) :: rk4min_sfcw_temp     ! Minimum snow/pond temperature         [        K]
   real(kind=8) :: rk4max_sfcw_temp     ! Maximum snow/pond temperature         [        K]
   !----- The following variables have units different from the actual value. -------------!
   real(kind=8) :: rk4min_veg_lwater    ! Maximum leaf      surface water mass  [kg/m�leaf]
   real(kind=8) :: rk4min_sfcw_moist    ! Maximum snow/pond water mass          [m�/m�soil]
   real(kind=8) :: rk4min_virt_moist    ! Minimum virtual pool mass             [m�/m�soil]
   !----- The following variables will be defined in sfcdata_ed (ed_init.f90). ------------!
   real(kind=8) :: rk4min_sfcw_mass     ! Minimum snow/pond    mass             [    kg/m�]
   real(kind=8) :: rk4min_virt_water    ! Minimum virtual pool mass             [    kg/m�]

   !----- MOVED THE FOLLOWING BLOCK TO RK4AUX BECAUSE WITH SMP WE NEED MULTIPLES
   !----- The following variables will be defined every time step. ------------------------!
!   real(kind=8) :: rk4min_can_theta     ! Minimum canopy    potential temp.     [        K]
!   real(kind=8) :: rk4max_can_theta     ! Maximum canopy    potential temp.     [        K]
!   real(kind=8) :: rk4min_can_enthalpy  ! Minimum canopy    enthalpy            [     J/m2]
!   real(kind=8) :: rk4max_can_enthalpy  ! Maximum canopy    enthalpy            [     J/m2]
!   real(kind=8) :: rk4min_can_prss      ! Minimum canopy    pressure            [       Pa]
!   real(kind=8) :: rk4max_can_prss      ! Maximum canopy    pressure            [       Pa]
   !---------------------------------------------------------------------------------------!


   !----- The following variables are double precision version of some bounds. ------------!
   real(kind=8) :: rk4tiny_sfcw_mass     ! Min. non-negligible snow/pond mass   [    kg/m�]
   real(kind=8) :: rk4water_stab_thresh  ! Min. mass for a layer to be stable   [    kg/m�]
   real(kind=8) :: rk4snowmin            ! Min. snow mass required for new lyr  [    kg/m�]
   real(kind=8) :: rk4leaf_drywhc        ! Min. non-negligible leaf water mass  [kg/m�leaf]
   real(kind=8) :: rk4leaf_maxwhc        ! Max. leaf water mass possible        [kg/m�leaf]
   real(kind=8) :: rk4tiny_sfcw_depth    ! Minimum snow/pond depth              [        m]
   !---------------------------------------------------------------------------------------!

   !----- MOVED THE FOLLOWING BLOCK TO RK4AUX BECAUSE WITH SMP WE NEED MULTIPLES
  !---------------------------------------------------------------------------------------!
  !     The following variables are going to be allocated at the first time step, and     !
   ! they will be re-calculated every time the integrator is called.                       !
   ! rk4min_soil_water                  ! Minimum soil moisture                 [    m�/m�]!
   ! rk4max_soil_water                  ! Maximum soil moisture                 [    m�/m�]!
   !---------------------------------------------------------------------------------------!
!   real(kind=8), dimension(:), allocatable :: rk4min_soil_water
!   real(kind=8), dimension(:), allocatable :: rk4max_soil_water
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     These variables are assigned during the integration process. They are not set at  !
   ! ed_params.f90, and their values may change during the simulation.                     !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Integration limits for time. Since 'derivs' do not explicitly depend on time it   !
   ! doesn't really matter what this is as long as tend-tbeg makes sense.                  !
   !---------------------------------------------------------------------------------------!
   real(kind=8)   :: tbeg
   real(kind=8)   :: tend
   real(kind=8)   :: dtrk4
   real(kind=8)   :: dtrk4i
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     These variables are assigned in ed_params.f90.  Heat area should be 2.0 for all   !
   ! PFTs (two sides of the leaves exchange heat), and the evaporation area should be 1.0  !
   ! for all PFTs (only one side of the leaf is usually covered by water).  The transpir-  !
   ! ation area should be 1.0 for hypostomatous leaves, and 2.0 for symmetrical (pines)    !
   ! and amphistomatous (araucarias) leaves.  Sometimes heat and evaporation are multi-    !
   ! plied  by 1.2 and 2.2 to account for branches and twigs.  This is not recommended,    !
   ! though, because branches and twigs do not contribute to heat storage when             !
   ! ibranch_thermo is set to zero, and they are otherwise accounted through the wood area !
   ! index.                                                                                !
   !---------------------------------------------------------------------------------------!
   real(kind=8)                   :: effarea_heat   ! Heat area: related to 2*LAI
   real(kind=8)                   :: effarea_evap   ! Evaporation area: related to LAI
   real(kind=8), dimension(n_pft) :: effarea_transp ! Evaporation area: related to LAI
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Flag to determine whether the patch is too sparsely populated to be computed at   !
   ! the cohort level.  The decision is made based on the difference in order of magnitude !
   ! between the patch "natural" leaf heat capacity and the minimum heat capacity for the  !
   ! Runge-Kutta solver (checked at copy_patch_init).                                      !
   !---------------------------------------------------------------------------------------!
!!!   logical :: toosparse  ! NEVER USED RGK
   !---------------------------------------------------------------------------------------!

   !----- Flag to tell whether there is at least one "resolvable" cohort in this patch ----!
!!!   logical :: any_resolvable   ! MOVED TO RK4AUX FOR MULTI-THREADING
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Canopy air space capacities.  These variables are used to convert the intensive  !
   ! version of canopy air space prognostic variables (specific enthalpy, water vapour     !
   ! specific humidity and CO2 mixing ratio) into extensive variables.                     ! 
   !---------------------------------------------------------------------------------------!
   !! MOVED TO RK4AUX FOR MULTI-THREADING
!   real(kind=8) :: wcapcan  ! Water capacity                               [  kg_air/m�gnd]
!   real(kind=8) :: hcapcan  ! Enthalpy capacity                            [  kg_air/m�gnd]
!   real(kind=8) :: ccapcan  ! CO2 capacity                                 [ mol_air/m�gnd]
!   real(kind=8) :: wcapcani ! Inverse of water capacity                    [  m�gnd/kg_air]
!   real(kind=8) :: hcapcani ! Inverse of enthalpy capacity                 [  m�gnd/kg_air]
!   real(kind=8) :: ccapcani ! Inverse of CO2 capacity                      [ m�gnd/mol_air]
   !---------------------------------------------------------------------------------------!

   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      Integrator error statistics.                                                     !
   !---------------------------------------------------------------------------------------!
   !----- Number of variables other than soil and surface that will be analysed. ----------!
   integer                          , parameter   :: nerrfix = 21
   !---------------------------------------------------------------------------------------!


   !----- Total number of variables that will be analysed. --------------------------------!
   integer                                        :: nerr
   !---------------------------------------------------------------------------------------!


   !----- Default file name to be given to the error files. -------------------------------!
   character(len=str_len)                         :: errmax_fout
   character(len=str_len)                         :: sanity_fout
   character(len=str_len)                         :: thbnds_fout
   character(len=str_len)                         :: detail_pref
   character(len=str_len)                         :: budget_pref
   !---------------------------------------------------------------------------------------!


   !----- The error counter and label. ----------------------------------------------------!
   integer(kind=8)  , dimension(:,:), allocatable :: integ_err
   character(len=13), dimension(:)  , allocatable :: integ_lab
   !---------------------------------------------------------------------------------------!


   !----- Offset needed for each of the following variables. ------------------------------!
   integer                                        :: osow ! Soil water.
   integer                                        :: osoe ! Soil energy.
   integer                                        :: oswe ! Surface water energy.
   integer                                        :: oswm ! Surface water mass.
   !---------------------------------------------------------------------------------------!
   !=======================================================================================!
   !=======================================================================================!


   contains


     ! ==========================================================
     ! The next three subroutines are for allocating
     ! integration memory for the previous time-step's
     ! leaf, wood and canopy temperature.  This is needed
     ! only in the BDF2 implicit solver method.
     ! =========================================================

     subroutine allocate_bdf2_patch(y,maxcohort)
       
       implicit none
       type(bdf2patchtype), target :: y
       integer :: maxcohort
       
       allocate(y%leaf_temp(maxcohort))
       allocate(y%wood_temp(maxcohort))
       
       return
     end subroutine allocate_bdf2_patch


     subroutine deallocate_bdf2_patch(y)

       implicit none
       type(bdf2patchtype),target :: y
       
       deallocate(y%leaf_temp)
       deallocate(y%wood_temp)
       return
     end subroutine deallocate_bdf2_patch

     
     subroutine nullify_bdf2_patch(y)
       
       implicit none
       type(bdf2patchtype), target :: y
       
       nullify(y%leaf_temp)
       nullify(y%wood_temp)
       return
     end subroutine nullify_bdf2_patch



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will perform the temporary patch allocation for the RK4           !
   ! integration.                                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine allocate_rk4_patch(y)
      use grid_coms     , only : nzg          & ! intent(in)
                               , nzs          ! ! intent(in)
      implicit none
      !----- Argument ---------------------------------------------------------------------!
      type(rk4patchtype), target :: y
      !------------------------------------------------------------------------------------!

      call nullify_rk4_patch(y)

      allocate(y%soil_energy            (0:nzg))
      allocate(y%soil_mstpot            (0:nzg))
      allocate(y%soil_water             (0:nzg))
      allocate(y%soil_fracliq           (0:nzg))
      allocate(y%soil_tempk             (0:nzg))

      allocate(y%sfcwater_energy          (nzs))
      allocate(y%sfcwater_mass            (nzs))
      allocate(y%sfcwater_depth           (nzs))
      allocate(y%sfcwater_fracliq         (nzs))
      allocate(y%sfcwater_tempk           (nzs))

      !------------------------------------------------------------------------------------!
      !     Diagnostics - for now we will always allocate the diagnostics, even if they    !
      !                   aren't used.                                                     !
      !------------------------------------------------------------------------------------!
      allocate(y%avg_sensible_gg          (nzg))
      allocate(y%avg_smoist_gg            (nzg))
      allocate(y%avg_transloss            (nzg))
      allocate(y%flx_sensible_gg          (nzg))
      allocate(y%flx_smoist_gg            (nzg))
      allocate(y%flx_transloss            (nzg))

      call zero_rk4_patch(y)

      return
   end subroutine allocate_rk4_patch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will nullify all pointers, to make a safe allocation.             !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_rk4_patch(y)
      implicit none
      !----- Argument ---------------------------------------------------------------------!
      type(rk4patchtype), target :: y
      !------------------------------------------------------------------------------------!

      nullify(y%soil_energy )
      nullify(y%soil_mstpot )
      nullify(y%soil_water  )
      nullify(y%soil_fracliq)
      nullify(y%soil_tempk  )

      nullify(y%sfcwater_energy )
      nullify(y%sfcwater_mass   )
      nullify(y%sfcwater_depth  )
      nullify(y%sfcwater_fracliq)
      nullify(y%sfcwater_tempk  )

      !------------------------------------------------------------------------------------!
      !     Diagnostics - for now we will always allocate the diagnostics, even if they    !
      !                   aren't used.                                                     !
      !------------------------------------------------------------------------------------!
      nullify(y%avg_smoist_gg  )
      nullify(y%avg_transloss  )
      nullify(y%avg_sensible_gg)
      nullify(y%flx_smoist_gg  )
      nullify(y%flx_transloss  )
      nullify(y%flx_sensible_gg)

      return
   end subroutine nullify_rk4_patch
   !=======================================================================================!
   !=======================================================================================!





  subroutine zero_bdf2_patch(y)

     implicit none
     type(bdf2patchtype),target :: y
     
     y%can_temp  = 0.d0
     y%leaf_temp = 0.d0
     y%wood_temp = 0.d0

     return
   end subroutine zero_bdf2_patch
     


   !=======================================================================================!
   !=======================================================================================!
   !    Forcing all variables to be zero.                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine zero_rk4_patch(y)
      use grid_coms     , only : nzg          & ! intent(in)
                               , nzs          ! ! intent(in)
      implicit none
      !----- Argument ---------------------------------------------------------------------!
      type(rk4patchtype), target :: y
      !------------------------------------------------------------------------------------!

      y%co2budget_storage              = 0.d0
      y%co2budget_loss2atm             = 0.d0
      y%ebudget_storage                = 0.d0
      y%ebudget_netrad                 = 0.d0
      y%ebudget_loss2atm               = 0.d0
      y%ebudget_loss2drainage          = 0.d0
      y%ebudget_loss2runoff            = 0.d0
      y%wbudget_storage                = 0.d0
      y%wbudget_loss2atm               = 0.d0
      y%wbudget_loss2drainage          = 0.d0
      y%wbudget_loss2runoff            = 0.d0
     
      y%can_temp                       = 0.d0
      y%can_shv                        = 0.d0
      y%can_ssh                        = 0.d0
      y%can_rhv                        = 0.d0
      y%can_co2                        = 0.d0
      y%can_theta                      = 0.d0
      y%can_enthalpy                   = 0.d0
      y%can_depth                      = 0.d0
      y%can_rhos                       = 0.d0
      y%can_prss                       = 0.d0
      y%can_exner                      = 0.d0
      y%can_cp                         = 0.d0
      y%veg_height                     = 0.d0

!      y%vels                           = 0.d0
!      y%atm_enthalpy                   = 0.d0

      y%veg_displace                   = 0.d0
      y%veg_rough                      = 0.d0
      y%opencan_frac                   = 0.d0
      y%total_sfcw_depth               = 0.d0
      y%total_sfcw_mass                = 0.d0
      y%snowfac                        = 0.d0
      y%ggbare                         = 0.d0
      y%ggveg                          = 0.d0
      y%ggnet                          = 0.d0
      y%ggsoil                         = 0.d0
      y%flag_wflxgc                    = -1
      y%virtual_water                  = 0.d0
      y%virtual_energy                 = 0.d0
      y%virtual_depth                  = 0.d0
      y%virtual_tempk                  = 0.d0
      y%virtual_fracliq                = 0.d0
     
      y%ground_shv                     = 0.d0
      y%ground_ssh                     = 0.d0
      y%ground_temp                    = 0.d0
      y%ground_fliq                    = 0.d0
      y%nlev_sfcwater                  = 0
      y%flag_sfcwater                  = 0

      y%rough                          = 0.d0

      y%ustar                          = 0.d0
      y%cstar                          = 0.d0
      y%tstar                          = 0.d0
      y%qstar                          = 0.d0
      y%estar                          = 0.d0

      y%zeta                           = 0.d0
      y%ribulk                         = 0.d0

      y%rasveg                         = 0.d0
      y%root_res_fac                   = 0.d0
      y%cwd_rh                         = 0.d0
      y%rh                             = 0.d0

      y%upwp                           = 0.d0
      y%wpwp                           = 0.d0
      y%tpwp                           = 0.d0
      y%qpwp                           = 0.d0
      y%cpwp                           = 0.d0
      
      y%rasveg                         = 0.d0

      y%avg_ustar                      = 0.d0
      y%avg_tstar                      = 0.d0
      y%avg_qstar                      = 0.d0
      y%avg_cstar                      = 0.d0

      y%avg_carbon_ac                  = 0.d0
      y%avg_carbon_st                  = 0.d0
      y%avg_vapor_gc                   = 0.d0
      y%avg_throughfall                = 0.d0
      y%avg_vapor_ac                   = 0.d0
      y%avg_rshort_gnd                 = 0.d0
      y%avg_par_gnd                    = 0.d0
      y%avg_rlong_gnd                  = 0.d0
      y%avg_qthroughfall               = 0.d0
      y%avg_sensible_gc                = 0.d0
      y%avg_sensible_ac                = 0.d0
      y%avg_heatstor_veg               = 0.d0

      y%avg_drainage                   = 0.d0
      y%avg_qdrainage                  = 0.d0

      y%water_deficit                  = 0.d0

      y%flx_carbon_ac                  = 0.d0
      y%flx_carbon_st                  = 0.d0
      y%flx_vapor_lc                   = 0.d0
      y%flx_vapor_wc                   = 0.d0
      y%flx_vapor_gc                   = 0.d0
      y%flx_wshed_vg                   = 0.d0
      y%flx_intercepted                = 0.d0
      y%flx_throughfall                = 0.d0
      y%flx_vapor_ac                   = 0.d0
      y%flx_transp                     = 0.d0
      y%flx_rshort_gnd                 = 0.d0
      y%flx_par_gnd                    = 0.d0
      y%flx_rlong_gnd                  = 0.d0
      y%flx_sensible_lc                = 0.d0
      y%flx_sensible_wc                = 0.d0
      y%flx_qwshed_vg                  = 0.d0
      y%flx_qintercepted               = 0.d0
      y%flx_qthroughfall               = 0.d0
      y%flx_sensible_gc                = 0.d0
      y%flx_sensible_ac                = 0.d0
      y%flx_heatstor_veg               = 0.d0

      y%hflxgc                         = 0.d0
      y%hflxsc                         = 0.d0
      y%wflxgc                         = 0.d0
      y%wflxsc                         = 0.d0
      y%wflxac                         = 0.d0
      y%qwflxgc                        = 0.d0
      y%qwflxsc                        = 0.d0

      y%flx_drainage                   = 0.d0
      y%flx_qdrainage                  = 0.d0

      !----- The following variables are pointers, check whether they are linked or not. --!
      if(associated(y%soil_energy           ))   y%soil_energy     (:)            = 0.d0
      if(associated(y%soil_mstpot           ))   y%soil_mstpot     (:)            = 0.d0
      if(associated(y%soil_tempk            ))   y%soil_tempk      (:)            = 0.d0
      if(associated(y%soil_fracliq          ))   y%soil_fracliq    (:)            = 0.d0
      if(associated(y%soil_water            ))   y%soil_water      (:)            = 0.d0
     
      if(associated(y%sfcwater_depth        ))   y%sfcwater_depth  (:)            = 0.d0
      if(associated(y%sfcwater_mass         ))   y%sfcwater_mass   (:)            = 0.d0
      if(associated(y%sfcwater_energy       ))   y%sfcwater_energy (:)            = 0.d0
      if(associated(y%sfcwater_tempk        ))   y%sfcwater_tempk  (:)            = 0.d0
      if(associated(y%sfcwater_fracliq      ))   y%sfcwater_fracliq(:)            = 0.d0

      if(associated(y%avg_smoist_gg         ))   y%avg_smoist_gg   (:)            = 0.d0
      if(associated(y%avg_transloss         ))   y%avg_transloss   (:)            = 0.d0
      if(associated(y%avg_sensible_gg       ))   y%avg_sensible_gg (:)            = 0.d0

      if(associated(y%flx_smoist_gg         ))   y%flx_smoist_gg   (:)            = 0.d0
      if(associated(y%flx_transloss         ))   y%flx_transloss   (:)            = 0.d0
      if(associated(y%flx_sensible_gg       ))   y%flx_sensible_gg (:)            = 0.d0

      return
   end subroutine zero_rk4_patch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will perform the temporary patch deallocation.                     !
   !---------------------------------------------------------------------------------------!
   subroutine deallocate_rk4_patch(y)
      implicit none
      !----- Argument ---------------------------------------------------------------------!
      type(rk4patchtype), target :: y
      !------------------------------------------------------------------------------------!

      if (associated(y%soil_energy ))            deallocate(y%soil_energy )
      if (associated(y%soil_mstpot ))            deallocate(y%soil_mstpot )
      if (associated(y%soil_water  ))            deallocate(y%soil_water  )
      if (associated(y%soil_fracliq))            deallocate(y%soil_fracliq)
      if (associated(y%soil_tempk  ))            deallocate(y%soil_tempk  )

      if (associated(y%sfcwater_energy ))        deallocate(y%sfcwater_energy )
      if (associated(y%sfcwater_mass   ))        deallocate(y%sfcwater_mass   )
      if (associated(y%sfcwater_depth  ))        deallocate(y%sfcwater_depth  )
      if (associated(y%sfcwater_fracliq))        deallocate(y%sfcwater_fracliq)
      if (associated(y%sfcwater_tempk  ))        deallocate(y%sfcwater_tempk  )

      ! Diagnostics
      if (associated(y%avg_smoist_gg  ))         deallocate(y%avg_smoist_gg  )
      if (associated(y%avg_transloss  ))         deallocate(y%avg_transloss  )
      if (associated(y%avg_sensible_gg))         deallocate(y%avg_sensible_gg)
      if (associated(y%flx_smoist_gg  ))         deallocate(y%flx_smoist_gg  )
      if (associated(y%flx_transloss  ))         deallocate(y%flx_transloss  )
      if (associated(y%flx_sensible_gg))         deallocate(y%flx_sensible_gg)

      return
   end subroutine deallocate_rk4_patch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will allocate the cohorts of the temporary patch.                  !
   !---------------------------------------------------------------------------------------!

   subroutine allocate_rk4_coh(maxcohort,y)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target     :: y
      integer            , intent(in) :: maxcohort
      !------------------------------------------------------------------------------------!
      
      call nullify_rk4_cohort(y)


      allocate(y%leaf_energy      (maxcohort))
      allocate(y%leaf_water       (maxcohort))
      allocate(y%leaf_temp        (maxcohort))
      allocate(y%leaf_fliq        (maxcohort))
      allocate(y%leaf_hcap        (maxcohort))
      allocate(y%leaf_reynolds    (maxcohort))
      allocate(y%leaf_grashof     (maxcohort))
      allocate(y%leaf_nussfree    (maxcohort))
      allocate(y%leaf_nussforc    (maxcohort))
      allocate(y%lint_shv         (maxcohort))
      allocate(y%leaf_resolvable  (maxcohort))
      allocate(y%leaf_gbh         (maxcohort))
      allocate(y%leaf_gbw         (maxcohort))
      allocate(y%gsw_open         (maxcohort))
      allocate(y%gsw_closed       (maxcohort))
      allocate(y%rshort_l         (maxcohort))
      allocate(y%rlong_l          (maxcohort))
      allocate(y%wood_energy      (maxcohort))
      allocate(y%wood_water       (maxcohort))
      allocate(y%wood_temp        (maxcohort))
      allocate(y%wood_fliq        (maxcohort))
      allocate(y%wood_hcap        (maxcohort))
      allocate(y%wood_reynolds    (maxcohort))
      allocate(y%wood_grashof     (maxcohort))
      allocate(y%wood_nussfree    (maxcohort))
      allocate(y%wood_nussforc    (maxcohort))
      allocate(y%wood_resolvable  (maxcohort))
      allocate(y%wood_gbh         (maxcohort))
      allocate(y%wood_gbw         (maxcohort))
      allocate(y%rshort_w         (maxcohort))
      allocate(y%rlong_w          (maxcohort))
      allocate(y%veg_energy       (maxcohort))
      allocate(y%veg_water        (maxcohort))
      allocate(y%veg_hcap         (maxcohort))
      allocate(y%veg_resolvable   (maxcohort))
      allocate(y%nplant           (maxcohort))
      allocate(y%veg_wind         (maxcohort))
      allocate(y%lai              (maxcohort))
      allocate(y%wai              (maxcohort))
      allocate(y%tai              (maxcohort))
      allocate(y%crown_area       (maxcohort))
      allocate(y%elongf           (maxcohort))
      allocate(y%psi_open         (maxcohort))
      allocate(y%psi_closed       (maxcohort))
      allocate(y%fs_open          (maxcohort))
      allocate(y%gpp              (maxcohort))
      allocate(y%leaf_resp        (maxcohort))
      allocate(y%root_resp        (maxcohort))
      allocate(y%growth_resp      (maxcohort))
      allocate(y%storage_resp     (maxcohort))
      allocate(y%vleaf_resp       (maxcohort))

      allocate(y%wflxlc           (maxcohort))
      allocate(y%wflxwc           (maxcohort))
      allocate(y%wflxtr           (maxcohort))
      allocate(y%hflx_wrsti       (maxcohort))
      allocate(y%hflx_lrsti       (maxcohort))

      allocate(y%cfx_hflxlc       (maxcohort))
      allocate(y%cfx_hflxwc       (maxcohort))
      allocate(y%cfx_qwflxlc      (maxcohort))
      allocate(y%cfx_qwflxwc      (maxcohort))
      allocate(y%cfx_qwshed       (maxcohort))
      allocate(y%cfx_qtransp      (maxcohort))
      allocate(y%cfx_qintercepted (maxcohort))

      allocate(y%avg_sensible_lc   (maxcohort))
      allocate(y%avg_sensible_wc   (maxcohort))
      allocate(y%avg_vapor_lc      (maxcohort))
      allocate(y%avg_vapor_wc      (maxcohort))
      allocate(y%avg_transp        (maxcohort))
      allocate(y%avg_intercepted_al(maxcohort))
      allocate(y%avg_intercepted_aw(maxcohort))
      allocate(y%avg_wshed_lg      (maxcohort))
      allocate(y%avg_wshed_wg      (maxcohort))


      call zero_rk4_cohort(y)

      return
   end subroutine allocate_rk4_coh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will nullify the cohort pointers for a safe allocation.            !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_rk4_cohort(y)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype), target :: y
      !------------------------------------------------------------------------------------!

      nullify(y%leaf_energy      )
      nullify(y%leaf_water       )
      nullify(y%leaf_temp        )
      nullify(y%leaf_fliq        )
      nullify(y%leaf_hcap        )
      nullify(y%leaf_reynolds    )
      nullify(y%leaf_grashof     )
      nullify(y%leaf_nussfree    )
      nullify(y%leaf_nussforc    )
      nullify(y%lint_shv         )
      nullify(y%leaf_resolvable  )
      nullify(y%leaf_gbh         )
      nullify(y%leaf_gbw         )
      nullify(y%gsw_open         )
      nullify(y%gsw_closed       )
      nullify(y%rshort_l         )
      nullify(y%rlong_l          )
      nullify(y%wood_energy      )
      nullify(y%wood_water       )
      nullify(y%wood_temp        )
      nullify(y%wood_fliq        )
      nullify(y%wood_hcap        )
      nullify(y%wood_reynolds    )
      nullify(y%wood_grashof     )
      nullify(y%wood_nussfree    )
      nullify(y%wood_nussforc    )
      nullify(y%wood_resolvable  )
      nullify(y%wood_gbh         )
      nullify(y%wood_gbw         )
      nullify(y%rshort_w         )
      nullify(y%rlong_w          )
      nullify(y%veg_energy       )
      nullify(y%veg_water        )
      nullify(y%veg_hcap         )
      nullify(y%veg_resolvable   )
      nullify(y%nplant           )
      nullify(y%veg_wind         )
      nullify(y%lai              )
      nullify(y%wai              )
      nullify(y%tai              )
      nullify(y%crown_area       )
      nullify(y%elongf           )
      nullify(y%psi_open         )
      nullify(y%psi_closed       )
      nullify(y%fs_open          )
      nullify(y%gpp              )
      nullify(y%leaf_resp        )
      nullify(y%root_resp        )
      nullify(y%growth_resp      )
      nullify(y%storage_resp     )
      nullify(y%vleaf_resp       )

      nullify(y%wflxlc           )
      nullify(y%wflxwc           )
      nullify(y%wflxtr           )
      nullify(y%hflx_lrsti       )
      nullify(y%hflx_wrsti       )

      nullify(y%cfx_hflxlc       )
      nullify(y%cfx_hflxwc       )
      nullify(y%cfx_qwflxlc      )
      nullify(y%cfx_qwflxwc      )
      nullify(y%cfx_qwshed       )
      nullify(y%cfx_qtransp      )
      nullify(y%cfx_qintercepted )



      nullify(y%avg_sensible_lc   )
      nullify(y%avg_sensible_wc   )
      nullify(y%avg_vapor_lc      )
      nullify(y%avg_vapor_wc      )
      nullify(y%avg_transp        )
      nullify(y%avg_intercepted_al)
      nullify(y%avg_intercepted_aw)
      nullify(y%avg_wshed_lg      )
      nullify(y%avg_wshed_wg      )

      return
   end subroutine nullify_rk4_cohort
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will initialize the cohort variables with zeroes.                  !
   !---------------------------------------------------------------------------------------!
   subroutine zero_rk4_cohort(y)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype), target :: y
      !------------------------------------------------------------------------------------!

      if (associated(y%leaf_energy      )) y%leaf_energy      = 0.d0
      if (associated(y%leaf_water       )) y%leaf_water       = 0.d0
      if (associated(y%leaf_temp        )) y%leaf_temp        = 0.d0
      if (associated(y%leaf_fliq        )) y%leaf_fliq        = 0.d0
      if (associated(y%leaf_hcap        )) y%leaf_hcap        = 0.d0
      if (associated(y%leaf_reynolds    )) y%leaf_reynolds    = 0.d0
      if (associated(y%leaf_grashof     )) y%leaf_grashof     = 0.d0
      if (associated(y%leaf_nussfree    )) y%leaf_nussfree    = 0.d0
      if (associated(y%leaf_nussforc    )) y%leaf_nussforc    = 0.d0
      if (associated(y%lint_shv         )) y%lint_shv         = 0.d0
      if (associated(y%leaf_resolvable  )) y%leaf_resolvable  = .false.
      if (associated(y%leaf_gbh         )) y%leaf_gbh         = 0.d0
      if (associated(y%leaf_gbw         )) y%leaf_gbw         = 0.d0
      if (associated(y%gsw_open         )) y%gsw_open         = 0.d0
      if (associated(y%gsw_closed       )) y%gsw_closed       = 0.d0
      if (associated(y%rshort_l         )) y%rshort_l         = 0.d0
      if (associated(y%rlong_l          )) y%rlong_l          = 0.d0
      if (associated(y%wood_energy      )) y%wood_energy      = 0.d0
      if (associated(y%wood_water       )) y%wood_water       = 0.d0
      if (associated(y%wood_temp        )) y%wood_temp        = 0.d0
      if (associated(y%wood_fliq        )) y%wood_fliq        = 0.d0
      if (associated(y%wood_hcap        )) y%wood_hcap        = 0.d0
      if (associated(y%wood_reynolds    )) y%wood_reynolds    = 0.d0
      if (associated(y%wood_grashof     )) y%wood_grashof     = 0.d0
      if (associated(y%wood_nussfree    )) y%wood_nussfree    = 0.d0
      if (associated(y%wood_nussforc    )) y%wood_nussforc    = 0.d0
      if (associated(y%wood_resolvable  )) y%wood_resolvable  = .false.
      if (associated(y%wood_gbh         )) y%wood_gbh         = 0.d0
      if (associated(y%wood_gbw         )) y%wood_gbw         = 0.d0
      if (associated(y%rshort_w         )) y%rshort_w         = 0.d0
      if (associated(y%rlong_w          )) y%rlong_w          = 0.d0
      if (associated(y%veg_energy       )) y%veg_energy       = 0.d0
      if (associated(y%veg_water        )) y%veg_water        = 0.d0
      if (associated(y%veg_hcap         )) y%veg_hcap         = 0.d0
      if (associated(y%veg_resolvable   )) y%veg_resolvable   = .false.
      if (associated(y%nplant           )) y%nplant           = 0.d0
      if (associated(y%veg_wind         )) y%veg_wind         = 0.d0
      if (associated(y%lai              )) y%lai              = 0.d0
      if (associated(y%wai              )) y%wai              = 0.d0
      if (associated(y%tai              )) y%tai              = 0.d0
      if (associated(y%crown_area       )) y%crown_area       = 0.d0
      if (associated(y%elongf           )) y%elongf           = 0.d0
      if (associated(y%psi_open         )) y%psi_open         = 0.d0
      if (associated(y%psi_closed       )) y%psi_closed       = 0.d0
      if (associated(y%fs_open          )) y%fs_open          = 0.d0
      if (associated(y%gpp              )) y%gpp              = 0.d0
      if (associated(y%leaf_resp        )) y%leaf_resp        = 0.d0
      if (associated(y%root_resp        )) y%root_resp        = 0.d0
      if (associated(y%growth_resp      )) y%growth_resp      = 0.d0
      if (associated(y%storage_resp     )) y%storage_resp     = 0.d0
      if (associated(y%vleaf_resp       )) y%vleaf_resp       = 0.d0

      if (associated(y%wflxlc           )) y%wflxlc           = 0.d0
      if (associated(y%wflxwc           )) y%wflxwc           = 0.d0
      if (associated(y%wflxtr           )) y%wflxtr           = 0.d0
      if (associated(y%hflx_wrsti       )) y%hflx_wrsti       = 0.d0      
      if (associated(y%hflx_lrsti       )) y%hflx_lrsti       = 0.d0

      if (associated(y%cfx_hflxlc       )) y%cfx_hflxlc       = 0.d0
      if (associated(y%cfx_hflxwc       )) y%cfx_hflxwc       = 0.d0
      if (associated(y%cfx_qwflxlc      )) y%cfx_qwflxlc      = 0.d0
      if (associated(y%cfx_qwflxwc      )) y%cfx_qwflxwc      = 0.d0
      if (associated(y%cfx_qwshed       )) y%cfx_qwshed       = 0.d0
      if (associated(y%cfx_qtransp      )) y%cfx_qtransp      = 0.d0
      if (associated(y%cfx_qintercepted )) y%cfx_qintercepted = 0.d0


      if(associated(y%avg_sensible_lc   )) y%avg_sensible_lc    = 0.d0
      if(associated(y%avg_sensible_wc   )) y%avg_sensible_wc    = 0.d0
      if(associated(y%avg_vapor_lc      )) y%avg_vapor_lc       = 0.d0
      if(associated(y%avg_vapor_wc      )) y%avg_vapor_wc       = 0.d0
      if(associated(y%avg_transp        )) y%avg_transp         = 0.d0
      if(associated(y%avg_intercepted_al)) y%avg_intercepted_al = 0.d0
      if(associated(y%avg_intercepted_aw)) y%avg_intercepted_aw = 0.d0
      if(associated(y%avg_wshed_lg      )) y%avg_wshed_lg       = 0.d0
      if(associated(y%avg_wshed_wg      )) y%avg_wshed_wg       = 0.d0

      return
   end subroutine zero_rk4_cohort
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will deallocate the cohorts of the temporary patch.                !
   !---------------------------------------------------------------------------------------!
   subroutine deallocate_rk4_coh(y)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype), target :: y
      !------------------------------------------------------------------------------------!

      if (associated(y%leaf_energy      )) deallocate(y%leaf_energy       )
      if (associated(y%leaf_water       )) deallocate(y%leaf_water        )
      if (associated(y%leaf_temp        )) deallocate(y%leaf_temp         )
      if (associated(y%leaf_fliq        )) deallocate(y%leaf_fliq         )
      if (associated(y%leaf_hcap        )) deallocate(y%leaf_hcap         )
      if (associated(y%leaf_reynolds    )) deallocate(y%leaf_reynolds     )
      if (associated(y%leaf_grashof     )) deallocate(y%leaf_grashof      )
      if (associated(y%leaf_nussfree    )) deallocate(y%leaf_nussfree     )
      if (associated(y%leaf_nussforc    )) deallocate(y%leaf_nussforc     )
      if (associated(y%lint_shv         )) deallocate(y%lint_shv          )
      if (associated(y%leaf_resolvable  )) deallocate(y%leaf_resolvable   )
      if (associated(y%leaf_gbh         )) deallocate(y%leaf_gbh          )
      if (associated(y%leaf_gbw         )) deallocate(y%leaf_gbw          )
      if (associated(y%gsw_open         )) deallocate(y%gsw_open          )
      if (associated(y%gsw_closed       )) deallocate(y%gsw_closed        )
      if (associated(y%rshort_l         )) deallocate(y%rshort_l          )
      if (associated(y%rlong_l          )) deallocate(y%rlong_l           )
      if (associated(y%wood_energy      )) deallocate(y%wood_energy       )
      if (associated(y%wood_water       )) deallocate(y%wood_water        )
      if (associated(y%wood_temp        )) deallocate(y%wood_temp         )
      if (associated(y%wood_fliq        )) deallocate(y%wood_fliq         )
      if (associated(y%wood_hcap        )) deallocate(y%wood_hcap         )
      if (associated(y%wood_reynolds    )) deallocate(y%wood_reynolds     )
      if (associated(y%wood_grashof     )) deallocate(y%wood_grashof      )
      if (associated(y%wood_nussfree    )) deallocate(y%wood_nussfree     )
      if (associated(y%wood_nussforc    )) deallocate(y%wood_nussforc     )
      if (associated(y%wood_resolvable  )) deallocate(y%wood_resolvable   )
      if (associated(y%wood_gbh         )) deallocate(y%wood_gbh          )
      if (associated(y%wood_gbw         )) deallocate(y%wood_gbw          )
      if (associated(y%rshort_w         )) deallocate(y%rshort_w          )
      if (associated(y%rlong_w          )) deallocate(y%rlong_w           )
      if (associated(y%veg_energy       )) deallocate(y%veg_energy        )
      if (associated(y%veg_water        )) deallocate(y%veg_water         )
      if (associated(y%veg_hcap         )) deallocate(y%veg_hcap          )
      if (associated(y%veg_resolvable   )) deallocate(y%veg_resolvable    )
      if (associated(y%nplant           )) deallocate(y%nplant            )
      if (associated(y%veg_wind         )) deallocate(y%veg_wind          )
      if (associated(y%lai              )) deallocate(y%lai               )
      if (associated(y%wai              )) deallocate(y%wai               )
      if (associated(y%tai              )) deallocate(y%tai               )
      if (associated(y%crown_area       )) deallocate(y%crown_area        )
      if (associated(y%elongf           )) deallocate(y%elongf            )
      if (associated(y%psi_open         )) deallocate(y%psi_open          )
      if (associated(y%psi_closed       )) deallocate(y%psi_closed        )
      if (associated(y%fs_open          )) deallocate(y%fs_open           )
      if (associated(y%gpp              )) deallocate(y%gpp               )
      if (associated(y%leaf_resp        )) deallocate(y%leaf_resp         )
      if (associated(y%root_resp        )) deallocate(y%root_resp         )
      if (associated(y%growth_resp      )) deallocate(y%growth_resp       )
      if (associated(y%storage_resp     )) deallocate(y%storage_resp      )
      if (associated(y%vleaf_resp       )) deallocate(y%vleaf_resp        )

      if (associated(y%wflxlc           )) deallocate(y%wflxlc            )
      if (associated(y%wflxwc           )) deallocate(y%wflxwc            )
      if (associated(y%wflxtr           )) deallocate(y%wflxtr            )
      if (associated(y%hflx_lrsti       )) deallocate(y%hflx_lrsti        )
      if (associated(y%hflx_wrsti       )) deallocate(y%hflx_wrsti        )


      if (associated(y%cfx_hflxlc       )) deallocate(y%cfx_hflxlc        )
      if (associated(y%cfx_hflxwc       )) deallocate(y%cfx_hflxwc        )
      if (associated(y%cfx_qwflxlc      )) deallocate(y%cfx_qwflxlc       )
      if (associated(y%cfx_qwflxwc      )) deallocate(y%cfx_qwflxwc       )
      if (associated(y%cfx_qwshed       )) deallocate(y%cfx_qwshed        )
      if (associated(y%cfx_qtransp      )) deallocate(y%cfx_qtransp       )
      if (associated(y%cfx_qintercepted )) deallocate(y%cfx_qintercepted  )



      if(associated(y%avg_sensible_lc   )) deallocate(y%avg_sensible_lc   )
      if(associated(y%avg_sensible_wc   )) deallocate(y%avg_sensible_wc   )
      if(associated(y%avg_vapor_lc      )) deallocate(y%avg_vapor_lc      )
      if(associated(y%avg_vapor_wc      )) deallocate(y%avg_vapor_wc      )
      if(associated(y%avg_transp        )) deallocate(y%avg_transp        )
      if(associated(y%avg_intercepted_al)) deallocate(y%avg_intercepted_al)
      if(associated(y%avg_intercepted_aw)) deallocate(y%avg_intercepted_aw)
      if(associated(y%avg_wshed_lg      )) deallocate(y%avg_wshed_lg      )
      if(associated(y%avg_wshed_wg      )) deallocate(y%avg_wshed_wg      )
      return
   end subroutine deallocate_rk4_coh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Re-set the instantaneous flux variables to zero.                                   !
   !---------------------------------------------------------------------------------------!
   subroutine reset_rk4_fluxes(y)
      use grid_coms     , only : nzg          & ! intent(in)
                               , nzs          ! ! intent(in)
      implicit none
      !----- Argument ---------------------------------------------------------------------!
      type(rk4patchtype), target :: y
      !------------------------------------------------------------------------------------!

      y%flx_carbon_ac                  = 0.d0
      y%flx_vapor_lc                   = 0.d0
      y%flx_vapor_wc                   = 0.d0
      y%flx_vapor_gc                   = 0.d0
      y%flx_wshed_vg                   = 0.d0
      y%flx_intercepted                = 0.d0
      y%flx_throughfall                = 0.d0
      y%flx_vapor_ac                   = 0.d0
      y%flx_transp                     = 0.d0
      y%flx_rshort_gnd                 = 0.d0
      y%flx_par_gnd                    = 0.d0
      y%flx_rlong_gnd                  = 0.d0
      y%flx_sensible_lc                = 0.d0
      y%flx_sensible_wc                = 0.d0
      y%flx_qwshed_vg                  = 0.d0
      y%flx_qintercepted               = 0.d0
      y%flx_qthroughfall               = 0.d0
      y%flx_sensible_gc                = 0.d0
      y%flx_sensible_ac                = 0.d0
      y%flx_heatstor_veg               = 0.d0

      y%flx_drainage                   = 0.d0
      y%flx_qdrainage                  = 0.d0

      !----- Reset soil fluxes only when they are allocated. ------------------------------!
      if(associated(y%flx_smoist_gg   )) y%flx_smoist_gg    (:) = 0.d0
      if(associated(y%flx_transloss   )) y%flx_transloss    (:) = 0.d0
      if(associated(y%flx_sensible_gg )) y%flx_sensible_gg  (:) = 0.d0
      !----- Reset cohort-level energy fluxes when they are allocated. --------------------!
      if(associated(y%cfx_hflxlc      )) y%cfx_hflxlc       (:) = 0.d0
      if(associated(y%cfx_hflxwc      )) y%cfx_hflxwc       (:) = 0.d0
      if(associated(y%cfx_qwflxlc     )) y%cfx_qwflxlc      (:) = 0.d0
      if(associated(y%cfx_qwflxwc     )) y%cfx_qwflxwc      (:) = 0.d0
      if(associated(y%cfx_qwshed      )) y%cfx_qwshed       (:) = 0.d0
      if(associated(y%cfx_qtransp     )) y%cfx_qtransp      (:) = 0.d0
      if(associated(y%cfx_qintercepted)) y%cfx_qintercepted (:) = 0.d0

      return
   end subroutine reset_rk4_fluxes
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This sub-routine normalises the fluxes by dividing it by the actual time step.     !
   !---------------------------------------------------------------------------------------!
   subroutine norm_rk4_fluxes(y,hdid)
      use grid_coms     , only : nzg          & ! intent(in)
                               , nzs          ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(rk4patchtype), target        :: y
      real(kind=8)      , intent(in)    :: hdid
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                      :: hdidi
      !------------------------------------------------------------------------------------!

      !----- Find the inverse of the time step. -------------------------------------------!
      hdidi = 1.d0 / hdid

      y%flx_carbon_ac     = y%flx_carbon_ac     * hdidi
      y%flx_vapor_lc      = y%flx_vapor_lc      * hdidi
      y%flx_vapor_wc      = y%flx_vapor_wc      * hdidi
      y%flx_vapor_gc      = y%flx_vapor_gc      * hdidi
      y%flx_wshed_vg      = y%flx_wshed_vg      * hdidi
      y%flx_intercepted   = y%flx_intercepted   * hdidi
      y%flx_throughfall   = y%flx_throughfall   * hdidi
      y%flx_vapor_ac      = y%flx_vapor_ac      * hdidi
      y%flx_transp        = y%flx_transp        * hdidi
      y%flx_rshort_gnd    = y%flx_rshort_gnd    * hdidi
      y%flx_par_gnd       = y%flx_par_gnd       * hdidi
      y%flx_rlong_gnd     = y%flx_rlong_gnd     * hdidi
      y%flx_sensible_lc   = y%flx_sensible_lc   * hdidi
      y%flx_sensible_wc   = y%flx_sensible_wc   * hdidi
      y%flx_qwshed_vg     = y%flx_qwshed_vg     * hdidi
      y%flx_qintercepted  = y%flx_qintercepted  * hdidi
      y%flx_qthroughfall  = y%flx_qthroughfall  * hdidi
      y%flx_sensible_gc   = y%flx_sensible_gc   * hdidi
      y%flx_sensible_ac   = y%flx_sensible_ac   * hdidi
      y%flx_heatstor_veg  = y%flx_heatstor_veg  * hdidi

      y%flx_drainage      = y%flx_drainage      * hdidi
      y%flx_qdrainage     = y%flx_qdrainage     * hdidi

      if(associated(y%flx_smoist_gg   )) y%flx_smoist_gg    (:) = y%flx_smoist_gg   (:)    &
                                                                * hdidi
      if(associated(y%flx_transloss   )) y%flx_transloss    (:) = y%flx_transloss   (:)    &
                                                                * hdidi
      if(associated(y%flx_sensible_gg )) y%flx_sensible_gg  (:) = y%flx_sensible_gg (:)    &
                                                                * hdidi
      !----- Cohort-level energy fluxes. --------------------------------------------------!
      if(associated(y%cfx_hflxlc      )) y%cfx_hflxlc       (:) = y%cfx_hflxlc      (:)    &
                                                                * hdidi
      if(associated(y%cfx_hflxwc      )) y%cfx_hflxwc       (:) = y%cfx_hflxwc      (:)    &
                                                                * hdidi
      if(associated(y%cfx_qwflxlc     )) y%cfx_qwflxlc      (:) = y%cfx_qwflxlc     (:)    &
                                                                * hdidi
      if(associated(y%cfx_qwflxwc     )) y%cfx_qwflxwc      (:) = y%cfx_qwflxwc     (:)    &
                                                                * hdidi
      if(associated(y%cfx_qwshed      )) y%cfx_qwshed       (:) = y%cfx_qwshed      (:)    &
                                                                * hdidi
      if(associated(y%cfx_qtransp     )) y%cfx_qtransp      (:) = y%cfx_qtransp     (:)    &
                                                                * hdidi
      if(associated(y%cfx_qintercepted)) y%cfx_qintercepted (:) = y%cfx_qintercepted(:)    &
                                                                * hdidi

      return
   end subroutine norm_rk4_fluxes
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will allocate the auxiliary variables.                             !
   !---------------------------------------------------------------------------------------!
   subroutine allocate_rk4_aux(crk4aux,mzg,mzs,mcoh)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer            , intent(in) :: mzg
      integer            , intent(in) :: mzs
      integer            , intent(in) :: mcoh
      !------------------------------------------------------------------------------------!
      type(rk4auxtype),target :: crk4aux


      !------ Nullify the structure for a safe allocation. --------------------------------!
      call nullify_rk4_aux(crk4aux)
      !------------------------------------------------------------------------------------!

      allocate(crk4aux%psiplusz                (          0:mzg) )
      allocate(crk4aux%hydcond                 (          0:mzg) )
      allocate(crk4aux%th_cond_s               (          0:mzg) )
      allocate(crk4aux%th_cond_p               (          0:mzs) )
      allocate(crk4aux%drysoil                 (          0:mzg) )
      allocate(crk4aux%satsoil                 (          0:mzg) )
      allocate(crk4aux%avail_h2o_lyr           (          mzg+1) )
      allocate(crk4aux%avail_h2o_int           (          mzg+1) )
      allocate(crk4aux%h_flux_g                (          mzg+1) )
      allocate(crk4aux%h_flux_s                (          mzs+1) )
      allocate(crk4aux%w_flux_g                (          mzg+1) )
      allocate(crk4aux%qw_flux_g               (          mzg+1) )
      allocate(crk4aux%extracted_water         (mcoh,     mzg+1) )

      allocate(crk4aux%rk4min_soil_water       (            mzg) )
      allocate(crk4aux%rk4max_soil_water       (            mzg) )


      !------ Flush the variables within this structure to zero. --------------------------!
      call zero_rk4_aux(crk4aux)
      !------------------------------------------------------------------------------------!
      return
   end subroutine allocate_rk4_aux
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will nullify the auxiliary pointers for a safe allocation.         !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_rk4_aux(crk4aux)
     
     implicit none
     type(rk4auxtype),target :: crk4aux
     
     nullify(crk4aux%psiplusz                )
     nullify(crk4aux%hydcond                 )
     nullify(crk4aux%th_cond_s               )
     nullify(crk4aux%th_cond_p               )
     nullify(crk4aux%drysoil                 )
     nullify(crk4aux%satsoil                 )
     nullify(crk4aux%avail_h2o_lyr           )
     nullify(crk4aux%avail_h2o_int           )
     nullify(crk4aux%extracted_water         )
     nullify(crk4aux%h_flux_g                )
     nullify(crk4aux%h_flux_s                )
     nullify(crk4aux%w_flux_g                )
     nullify(crk4aux%qw_flux_g               )
     nullify(crk4aux%rk4min_soil_water       )
     nullify(crk4aux%rk4max_soil_water       )
     
     return
   end subroutine nullify_rk4_aux
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will flush all auxiliary variables to zero.                        !
   !---------------------------------------------------------------------------------------!
   subroutine zero_rk4_aux(crk4aux)
     implicit none
     type(rk4auxtype),target :: crk4aux
     
     if (associated(crk4aux%psiplusz        )) crk4aux%psiplusz        (:  ) = 0.d0
     if (associated(crk4aux%hydcond         )) crk4aux%hydcond         (:  ) = 0.d0
     if (associated(crk4aux%th_cond_s       )) crk4aux%th_cond_s       (:  ) = 0.d0
     if (associated(crk4aux%th_cond_p       )) crk4aux%th_cond_p       (:  ) = 0.d0
     if (associated(crk4aux%drysoil         )) crk4aux%drysoil         (:  ) = .false.
     if (associated(crk4aux%satsoil         )) crk4aux%satsoil         (:  ) = .false.
     if (associated(crk4aux%avail_h2o_lyr   )) crk4aux%avail_h2o_lyr   (:  ) = 0.d0
     if (associated(crk4aux%avail_h2o_int   )) crk4aux%avail_h2o_int   (:  ) = 0.d0
     if (associated(crk4aux%h_flux_g        )) crk4aux%h_flux_g        (:  ) = 0.d0
     if (associated(crk4aux%h_flux_s        )) crk4aux%h_flux_s        (:  ) = 0.d0
     if (associated(crk4aux%w_flux_g        )) crk4aux%w_flux_g        (:  ) = 0.d0
     if (associated(crk4aux%qw_flux_g       )) crk4aux%qw_flux_g       (:  ) = 0.d0
     if (associated(crk4aux%extracted_water )) crk4aux%extracted_water (:,:) = 0.d0
     
!     crk4aux%rk4min_can_prss          = 0.d0
!     crk4aux%rk4max_can_prss          = 0.d0
!     crk4aux%rk4min_can_theta        = 0.d0
!     crk4aux%rk4max_can_theta        = 0.d0
!     crk4aux%rk4min_can_enthalpy     = 0.d0
!     crk4aux%rk4max_can_enthalpy     = 0.d0

!     if (associated(crk4aux%rk4min_soil_water))crk4aux%rk4min_soil_water(:)  = 0.d0
!     if (associated(crk4aux%rk4min_soil_water))crk4aux%rk4max_soil_water(:)  = 0.d0

     return
   end subroutine zero_rk4_aux
   !=======================================================================================!
   !=======================================================================================!
   





   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will deallocate the auxiliary variables.                           !
   !---------------------------------------------------------------------------------------!
   subroutine deallocate_rk4_aux(crk4aux)
     implicit none
     type(rk4auxtype),target :: crk4aux
     
     if (associated(crk4aux%psiplusz        )) deallocate(crk4aux%psiplusz        )
     if (associated(crk4aux%hydcond         )) deallocate(crk4aux%hydcond         )
     if (associated(crk4aux%th_cond_s       )) deallocate(crk4aux%th_cond_s       )
     if (associated(crk4aux%th_cond_p       )) deallocate(crk4aux%th_cond_p       )
     if (associated(crk4aux%drysoil         )) deallocate(crk4aux%drysoil         )
     if (associated(crk4aux%satsoil         )) deallocate(crk4aux%satsoil         )
     if (associated(crk4aux%avail_h2o_lyr   )) deallocate(crk4aux%avail_h2o_lyr   )
     if (associated(crk4aux%avail_h2o_int   )) deallocate(crk4aux%avail_h2o_int   )
     if (associated(crk4aux%h_flux_g        )) deallocate(crk4aux%h_flux_g        )
     if (associated(crk4aux%h_flux_s        )) deallocate(crk4aux%h_flux_s        )
     if (associated(crk4aux%w_flux_g        )) deallocate(crk4aux%w_flux_g        )
     if (associated(crk4aux%qw_flux_g       )) deallocate(crk4aux%qw_flux_g       )
     if (associated(crk4aux%extracted_water )) deallocate(crk4aux%extracted_water )
     
     if (associated(crk4aux%rk4min_soil_water)) deallocate(crk4aux%rk4min_soil_water)
     if (associated(crk4aux%rk4max_soil_water)) deallocate(crk4aux%rk4max_soil_water)

     return
   end subroutine deallocate_rk4_aux
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function allocates the integration error counter.                            !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_integ_err()
      use grid_coms, only : nzg & ! intent(in)
                          , nzs ! ! intent(in)
      implicit none
      !------------------------------------------------------------------------------------!
      !     Find the total number of variables to be analysed, and find the offset for     !
      ! soil and surface water (snow) properties.                                          !
      !------------------------------------------------------------------------------------!
      nerr = nerrfix + 2 * nzg + 2 * nzs
      osow = nerrfix
      osoe = osow + nzg
      oswe = osoe + nzg
      oswm = oswe + nzs
      !------------------------------------------------------------------------------------!

      if (.not. allocated(integ_err)) allocate(integ_err(nerr,2))
      if (.not. allocated(integ_lab)) allocate(integ_lab(nerr))

      return

   end subroutine alloc_integ_err
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine flushes the integrator error counter to zero.                    !
   !---------------------------------------------------------------------------------------!
   subroutine reset_integ_err()
      implicit none
      if (allocated(integ_err)) integ_err(:,:) = 0_8
      return
   end subroutine reset_integ_err
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function creates the labels for the error analysis.                          !
   !---------------------------------------------------------------------------------------!
   subroutine assign_err_label()
      use grid_coms, only : nzg & ! intent(in)
                          , nzs ! ! intent(in)
      implicit none
      !----- Local constants. -------------------------------------------------------------!
      character(len=13), dimension(nerrfix), parameter :: err_lab_fix = (/                 &
           'CAN_ENTHALPY ','CAN_THETA    ','CAN_SHV      ','CAN_TEMP     ','CAN_PRSS     ' &
          ,'CAN_CO2      ','LEAF_WATER   ','LEAF_ENERGY  ','WOOD_WATER   ','WOOD_ENERGY  ' &
          ,'VIRT_HEAT    ','VIRT_WATER   ','CO2B_STORAGE ','CO2B_LOSS2ATM','EB_NETRAD    ' &
          ,'EB_LOSS2ATM  ','WATB_LOSS2ATM','ENB_LOSS2DRA ','WATB_LOSS2DRA','ENB_STORAGE  ' &
          ,'WATB_STORAGE '/)
      !----- Local variables. -------------------------------------------------------------!
      integer                                          :: n
      character(len=13)                                :: err_lab_loc
      !------------------------------------------------------------------------------------!

      !----- First we fill the labels for the size-independent variables. -----------------!
      do n=1,nerrfix
         integ_lab(n) = adjustr(err_lab_fix(n))
      end do

      !----- Soil water. ------------------------------------------------------------------!
      do n=1,nzg
         write(err_lab_loc,fmt='(a,i2.2)') 'SOIL_WATER_',n
         integ_lab(osow+n) = adjustr(err_lab_loc)
      end do

      !----- Soil energy. -----------------------------------------------------------------!
      do n=1,nzg
         write(err_lab_loc,fmt='(a,i2.2)') 'SOIL_ENER_',n
         integ_lab(osoe+n) = adjustr(err_lab_loc)
      end do

      !----- Surface water energy. --------------------------------------------------------!
      do n=1,nzs
         write(err_lab_loc,fmt='(a,i2.2)') 'SFCW_ENER_',n
         integ_lab(oswe+n) = adjustr(err_lab_loc)
      end do

      !----- Surface water energy. --------------------------------------------------------!
      do n=1,nzs
         write(err_lab_loc,fmt='(a,i2.2)') 'SFCW_MASS_',n
         integ_lab(oswm+n) = adjustr(err_lab_loc)
      end do

      return
   end subroutine assign_err_label
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine find_derived_thbounds(can_rhos,can_theta,can_temp,                           &
                                    can_shv,can_prss,can_depth)
      use grid_coms   , only : nzg           ! ! intent(in)
      use consts_coms , only : rdry8         & ! intent(in)
                             , epim18        & ! intent(in)
                             , ep8           & ! intent(in)
                             , mmdryi8       & ! intent(in)
                             , hr_sec        & ! intent(in)
                             , min_sec       ! ! intent(in)

      use therm_lib8  , only : press2exner8  & ! function
                             , extemp2theta8 & ! function
                             , tq2enthalpy8  & ! function
                             , thetaeiv8     & ! function
                             , thetaeivs8    & ! function
                             , idealdenssh8  & ! function
                             , reducedpress8 ! ! function
      use soil_coms   , only : soil8         ! ! intent(in)
      use ed_misc_coms, only : current_time  ! ! intent(in)
      !$ use omp_lib
      
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8)                , intent(in) :: can_rhos
      real(kind=8)                , intent(in) :: can_theta
      real(kind=8)                , intent(in) :: can_temp
      real(kind=8)                , intent(in) :: can_shv
      real(kind=8)                , intent(in) :: can_prss
      real(kind=8)                , intent(in) :: can_depth
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                             :: can_prss_try
      real(kind=8)                             :: can_exner_try
      real(kind=8)                             :: can_theta_try
      real(kind=8)                             :: can_enthalpy_try
      integer                                  :: k
      integer                                  :: hour
      integer                                  :: minute
      integer                                  :: second
      integer                                  :: nsoil
      integer                                  :: ibuff
      !----- Local parameters and locally saved variables. --------------------------------!
      logical                     , save       :: firsttime    = .true.
      !------------------------------------------------------------------------------------!

      ibuff=1
      !$ ibuff = OMP_get_thread_num()+1

      !------------------------------------------------------------------------------------!
      !     When we call this sub-routine for the first time, we must allocate the sanity  !
      ! check bounds for soil moisture.  Also, if the debugger of the sanity check is      !
      ! going to be used, we write the file header.                                        !
      !------------------------------------------------------------------------------------!
!      if (firsttime) then
!
!         allocate(rk4min_soil_water(nzg))
!         allocate(rk4max_soil_water(nzg))
!
!         if (print_thbnd) then
!            open (unit=39,file=trim(thbnds_fout),status='replace',action='write')
!            write(unit=39,fmt='(16(a,1x))')  '        YEAR','       MONTH','         DAY'  &
!                                            ,'        HOUR','        MINU','        SECO'  &
!                                            ,'    MIN_TEMP','    MAX_TEMP','     MIN_SHV'  &
!                                            ,'     MAX_SHV','   MIN_THETA','   MAX_THETA'  &
!                                            ,'    MIN_PRSS','    MAX_PRSS','MIN_ENTHALPY'  &
!                                            ,'MAX_ENTHALPY'
!            close(unit=39,status='keep')
!         end if
!         firsttime = .false.
!      end if

      !------------------------------------------------------------------------------------!
      !     Find the bounds for pressure.  To avoid the pressure range to be too relaxed,  !
      ! switch one of the dependent variable a time, and use the current values for the    !
      ! other.  In addition, we force pressure to be bounded between the reduced pressure  !
      ! in case the reference height was different by the order of 10%.                    !
      !------------------------------------------------------------------------------------!
      !----- 1. Initial value, the most extreme one. --------------------------------------!
      rk4aux(ibuff)%rk4min_can_prss = &
            reducedpress8(rk4site%atm_prss,rk4site%atm_theta,rk4site%atm_shv &
            ,9.d-1*rk4site%geoht,can_theta,can_shv,can_depth)
      rk4aux(ibuff)%rk4max_can_prss = &
            reducedpress8(rk4site%atm_prss,rk4site%atm_theta,rk4site%atm_shv &
            ,1.1d0*rk4site%geoht,can_theta,can_shv,can_depth)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the bounds for potential temperature.  To avoid the pressure range to be  !
      ! too relaxed, switch one of the dependent variable a time, and use the current      !
      ! values for the other.                                                              !
      !------------------------------------------------------------------------------------!
      !----- 1. Initial value, the most extreme one. --------------------------------------!
      rk4aux(ibuff)%rk4min_can_theta   =  huge(1.d0)
      rk4aux(ibuff)%rk4max_can_theta   = -huge(1.d0)
      !----- 2. Minimum temperature. ------------------------------------------------------!
      can_exner_try      = press2exner8(can_prss)
      can_theta_try      = extemp2theta8(can_exner_try,rk4min_can_temp)
      rk4aux(ibuff)%rk4min_can_theta   = min(rk4aux(ibuff)%rk4min_can_theta,can_theta_try)
      rk4aux(ibuff)%rk4max_can_theta   = max(rk4aux(ibuff)%rk4max_can_theta,can_theta_try)
      !----- 3. Maximum temperature. ------------------------------------------------------!
      can_exner_try      = press2exner8(can_prss)
      can_theta_try      = extemp2theta8(can_exner_try,rk4max_can_temp)
      rk4aux(ibuff)%rk4min_can_theta   = min(rk4aux(ibuff)%rk4min_can_theta,can_theta_try)
      rk4aux(ibuff)%rk4max_can_theta   = max(rk4aux(ibuff)%rk4max_can_theta,can_theta_try)
      !----- 4. Minimum pressure. ---------------------------------------------------------!
      can_exner_try      = press2exner8(rk4aux(ibuff)%rk4min_can_prss)
      can_theta_try      = extemp2theta8(can_exner_try,can_temp)
      rk4aux(ibuff)%rk4min_can_theta   = min(rk4aux(ibuff)%rk4min_can_theta,can_theta_try)
      rk4aux(ibuff)%rk4max_can_theta   = max(rk4aux(ibuff)%rk4max_can_theta,can_theta_try)
      !----- 5. Maximum pressure. ---------------------------------------------------------!
      can_exner_try      = press2exner8(rk4aux(ibuff)%rk4max_can_prss)
      can_theta_try      = extemp2theta8(can_exner_try,can_temp)
      rk4aux(ibuff)%rk4min_can_theta   = min(rk4aux(ibuff)%rk4min_can_theta,can_theta_try)
      rk4aux(ibuff)%rk4max_can_theta   = max(rk4aux(ibuff)%rk4max_can_theta,can_theta_try)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Minimum and maximum enthalpy.                                                 !
      !------------------------------------------------------------------------------------!
      !----- 1. Initial value, the most extreme one. --------------------------------------!
      rk4aux(ibuff)%rk4min_can_enthalpy = huge(1.d0)
      rk4aux(ibuff)%rk4max_can_enthalpy = - huge(1.d0)
      !----- 2. Minimum temperature. ------------------------------------------------------!
      can_enthalpy_try    = tq2enthalpy8(rk4min_can_temp,can_shv,.true.)
      rk4aux(ibuff)%rk4min_can_enthalpy = min(rk4aux(ibuff)%rk4min_can_enthalpy,can_enthalpy_try)
      rk4aux(ibuff)%rk4max_can_enthalpy = max(rk4aux(ibuff)%rk4max_can_enthalpy,can_enthalpy_try)
      !----- 3. Maximum temperature. ------------------------------------------------------!
      can_enthalpy_try    = tq2enthalpy8(rk4max_can_temp,can_shv,.true.)
      rk4aux(ibuff)%rk4min_can_enthalpy = min(rk4aux(ibuff)%rk4min_can_enthalpy,can_enthalpy_try)
      rk4aux(ibuff)%rk4max_can_enthalpy = max(rk4aux(ibuff)%rk4max_can_enthalpy,can_enthalpy_try)
      !----- 4. Minimum specific humidity. ------------------------------------------------!
      can_enthalpy_try    = tq2enthalpy8(can_temp,rk4min_can_shv,.true.)
      rk4aux(ibuff)%rk4min_can_enthalpy = min(rk4aux(ibuff)%rk4min_can_enthalpy,can_enthalpy_try)
      rk4aux(ibuff)%rk4max_can_enthalpy = max(rk4aux(ibuff)%rk4max_can_enthalpy,can_enthalpy_try)
      !----- 5. Maximum specific humidity. ------------------------------------------------!
      can_enthalpy_try    = tq2enthalpy8(can_temp,rk4max_can_shv,.true.)
      rk4aux(ibuff)%rk4min_can_enthalpy = min(rk4aux(ibuff)%rk4min_can_enthalpy,can_enthalpy_try)
      rk4aux(ibuff)%rk4max_can_enthalpy = max(rk4aux(ibuff)%rk4max_can_enthalpy,can_enthalpy_try)
      !------------------------------------------------------------------------------------!

      if (print_thbnd) then
         hour   = floor(nint(current_time%time) / hr_sec)
         minute = floor((nint(current_time%time) - hour * hr_sec) / min_sec)
         second = mod(nint(current_time%time) - hour * hr_sec - minute * min_sec,min_sec)

         open (unit=39,file=trim(thbnds_fout),status='old',action='write'                  &
                      ,position='append')
         write (unit=39,fmt='(6(i12,1x),10(es12.5,1x))')                                   &
                                 current_time%year, current_time%month,  current_time%date &
                              ,               hour,             minute,             second &
                              ,    rk4min_can_temp,    rk4max_can_temp,     rk4min_can_shv &
                              ,     rk4max_can_shv,   rk4aux(ibuff)%rk4min_can_theta       &
                              , rk4aux(ibuff)%rk4max_can_theta, rk4aux(ibuff)%rk4min_can_prss &
                              , rk4aux(ibuff)%rk4max_can_prss,rk4aux(ibuff)%rk4min_can_enthalpy &
                              , rk4aux(ibuff)%rk4max_can_enthalpy
         close (unit=39,status='keep')
      end if

      do k = rk4site%lsl, nzg
         nsoil = rk4site%ntext_soil(k)
         rk4aux(ibuff)%rk4min_soil_water(k) = soil8(nsoil)%soilcp * (1.d0 - rk4eps)
         rk4aux(ibuff)%rk4max_soil_water(k) = soil8(nsoil)%slmsts * (1.d0 + rk4eps)
      end do

      return
   end subroutine find_derived_thbounds
   !=======================================================================================!
   !=======================================================================================!
end module rk4_coms
!==========================================================================================!
!==========================================================================================!

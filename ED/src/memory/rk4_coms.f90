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
                         , str_len ! ! intent(in)

   implicit none

   !=======================================================================================!
   !=======================================================================================!
   !     Structure used for the integration using the RK4 method.  Real variables are      !
   ! stored in double precision, so we can solve pools whose order of magnitude are very   !
   ! different.                                                                            !
   !---------------------------------------------------------------------------------------!
   type rk4patchtype

      !----- Canopy air variables. --------------------------------------------------------!
      real(kind=8)                        :: can_theiv    ! Eq. pot. temperature [       K]
      real(kind=8)                        :: can_lntheta  ! Log (theta)          [     ---]
      real(kind=8)                        :: can_theta    ! Pot. Temperature     [       K]
      real(kind=8)                        :: can_temp     ! Temperature          [       K]
      real(kind=8)                        :: can_shv      ! Specific humidity    [   kg/kg]
      real(kind=8)                        :: can_ssh      ! Sat. spec. humidity  [   kg/kg]
      real(kind=8)                        :: can_rvap     ! Vapour mixing ratio  [   kg/kg]
      real(kind=8)                        :: can_rhv      ! Relative humidity    [     ---]
      real(kind=8)                        :: can_co2      ! CO_2                 [µmol/mol]
      real(kind=8)                        :: can_depth    ! Canopy depth         [       m]
      real(kind=8)                        :: can_rhos     ! Canopy air density   [   kg/m³]
      real(kind=8)                        :: can_prss     ! Pressure             [      Pa]
      real(kind=8)                        :: can_exner    ! Exner function       [  J/kg/K]
      !------------------------------------------------------------------------------------!



      !----- Fraction of open canopy. -----------------------------------------------------!
      real(kind=8)                        :: opencan_frac ! Frac. of open canopy [     ---]
      !------------------------------------------------------------------------------------!



      !----- Ground -> Canopy flux type. --------------------------------------------------!
      real(kind=8)                        :: ggbare       ! Cond. of bare ground [     m/s]
      real(kind=8)                        :: ggveg        ! Cond. of veg. ground [     m/s]
      real(kind=8)                        :: ggnet        ! Net ground  conduct. [     m/s]
      integer                             :: flag_wflxgc  ! Flag for water flux.
      !------------------------------------------------------------------------------------!



      !----- Soil variables. --------------------------------------------------------------!
      real(kind=8), dimension(:), pointer :: soil_energy  ! Internal energy       [   J/m³]
      real(kind=8), dimension(:), pointer :: soil_tempk   ! Specific humidity     [      K]
      real(kind=8), dimension(:), pointer :: soil_fracliq ! Liquid fraction       [   ----]
      real(kind=8), dimension(:), pointer :: soil_water   ! Water content         [  m³/m³]
      real(kind=8), dimension(:), pointer :: soil_restz   ! Resistance term       [    m/s]
      real(kind=8), dimension(:), pointer :: psiplusz     ! Water potential       [      m]
      real(kind=8), dimension(:), pointer :: soilair99    ! 99% of field capacity [  m³/m³]
      real(kind=8), dimension(:), pointer :: soilair01    !  1% above wilting pt. [  m³/m³]
      real(kind=8), dimension(:), pointer :: soil_liq     ! Liquid fraction       [  kg/m²]
      real(kind=8), dimension(:), pointer :: available_liquid_water !             [  kg/m²]
      real(kind=8), dimension(:), pointer :: extracted_water        !             [  kg/m²]
      !------------------------------------------------------------------------------------!



      !----- Temporary surface water variables. -------------------------------------------!
      integer                             :: nlev_sfcwater    ! # of layers       [   ----]
      integer                             :: flag_sfcwater    ! Status flag       [  -----]
      real(kind=8), dimension(:), pointer :: sfcwater_depth   ! Depth             [      m]
      real(kind=8), dimension(:), pointer :: sfcwater_mass    ! Mass              [  kg/m²]
      real(kind=8), dimension(:), pointer :: sfcwater_energy  ! Internal energy   [   J/m²]
      real(kind=8), dimension(:), pointer :: sfcwater_tempk   ! Temperature       [      K]
      real(kind=8), dimension(:), pointer :: sfcwater_fracliq ! Liquid fraction   [   ----]
      real(kind=8), dimension(:), pointer :: sfcwater_restz   ! Resistance term   [   ----]
      !------------------------------------------------------------------------------------!



      !----- Virtual layer variables. -----------------------------------------------------!
      real(kind=8)                        :: virtual_water    ! Mass              [  kg/m²]
      real(kind=8)                        :: virtual_energy   ! Internal energy   [  kg/m²]
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
      real(kind=8)                        :: cstar  ! Carbon mixing ratio         [µmol/m³]
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



      !----- Heterotrophic respiration.[µmol/m²/s] ----------------------------------------!
      real(kind=8)                        :: cwd_rh
      real(kind=8)                        :: rh
      !------------------------------------------------------------------------------------!



      !----- Leaf (cohort-level) variables. -----------------------------------------------!
      real(kind=8), pointer, dimension(:) :: veg_energy   ! Internal energy     [     J/m²]
      real(kind=8), pointer, dimension(:) :: veg_water    ! Surface water mass  [    kg/m²]
      real(kind=8), pointer, dimension(:) :: veg_temp     ! Temperature         [        K]
      real(kind=8), pointer, dimension(:) :: veg_fliq     ! Liquid fraction     [      ---]
      real(kind=8), pointer, dimension(:) :: hcapveg      ! Heat capacity       [   J/m²/K]
      real(kind=8), pointer, dimension(:) :: veg_wind     ! Wind felt by cohort [      m/s]
      real(kind=8), pointer, dimension(:) :: veg_reynolds ! Reynolds number     [      ---]
      real(kind=8), pointer, dimension(:) :: veg_grashof  ! Grashof number      [      ---]
      real(kind=8), pointer, dimension(:) :: veg_nussfree ! Nusselt # (free)    [      ---]
      real(kind=8), pointer, dimension(:) :: veg_nussforc ! Nusselt # (forced)  [      ---]
      real(kind=8), pointer, dimension(:) :: lint_shv     ! Intercell. sp. hum. [    kg/kg]
      real(kind=8), pointer, dimension(:) :: nplant       ! Plant density       [ plant/m²]
      real(kind=8), pointer, dimension(:) :: lai          ! Leaf area index     [    m²/m²]
      real(kind=8), pointer, dimension(:) :: wai          ! Wood area index     [    m²/m²]
      real(kind=8), pointer, dimension(:) :: wpa          ! Wood projected area [    m²/m²]
      real(kind=8), pointer, dimension(:) :: tai          ! Tree area index     [    m²/m²]
      real(kind=8), pointer, dimension(:) :: crown_area   ! Crown area          [    m²/m²]
      real(kind=8), pointer, dimension(:) :: gbh          ! Leaf b.lyr. condct. [ J/K/m²/s]
      real(kind=8), pointer, dimension(:) :: gbw          ! Leaf b.lyr. condct. [  kg/m²/s]
      real(kind=8), pointer, dimension(:) :: gsw_open     ! Sto. condct. (op.)  [ J/K/m²/s]
      real(kind=8), pointer, dimension(:) :: gsw_closed   ! Sto. condct. (cl.)  [  kg/m²/s]
      real(kind=8), pointer, dimension(:) :: psi_open     ! Water demand (op.)  [  kg/m²/s]
      real(kind=8), pointer, dimension(:) :: psi_closed   ! Water demand (clos.)[  kg/m²/s]
      real(kind=8), pointer, dimension(:) :: fs_open      ! Frac. of op. stom.  [      ---]
      logical     , pointer, dimension(:) :: resolvable   ! resolve this cohort [      T|F]
      real(kind=8), pointer, dimension(:) :: gpp          ! Gross primary prod. [µmol/m²/s]
      real(kind=8), pointer, dimension(:) :: leaf_resp    ! Leaf respiration    [µmol/m²/s]
      real(kind=8), pointer, dimension(:) :: root_resp    ! Root respiration    [µmol/m²/s]
      real(kind=8), pointer, dimension(:) :: growth_resp  ! Growth respiration  [µmol/m²/s]
      real(kind=8), pointer, dimension(:) :: storage_resp ! Storage respiration [µmol/m²/s]
      real(kind=8), pointer, dimension(:) :: vleaf_resp   ! Virtual leaf resp.  [µmol/m²/s]
      real(kind=8), pointer, dimension(:) :: rshort_v     ! Net absorbed SWRad. [   J/m²/s]
      real(kind=8), pointer, dimension(:) :: rlong_v      ! Net absorbed LWRad. [   J/m²/s]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Fast time flux diagnostic variables.  These variables may be turned off under  !
      ! different conditions.                                                              !
      !------------------------------------------------------------------------------------!
      real(kind=8) :: avg_netrad        ! Net radiation
      !----- Water fluxes -----------------------------------------------------------------!
      real(kind=8) :: avg_vapor_vc ! Leaf       -> canopy air:  evap./cond. flux
      real(kind=8) :: avg_dew_cg   ! Canopy     -> ground    :  condensation flux
      real(kind=8) :: avg_vapor_gc ! Ground     -> canopy air:  evaporation flux
      real(kind=8) :: avg_vapor_ac ! Free atm.  -> canopy air:  vapour flux
      real(kind=8) :: avg_transp   ! Transpiration
      real(kind=8) :: avg_evap     ! Evaporation
      !----- Mass and energy fluxes due to water shedding. --------------------------------!
      real(kind=8) :: avg_wshed_vg      ! Mass flux
      real(kind=8) :: avg_qwshed_vg     ! Energy flux
      !----- Mass and energy input due to intercepted water. ------------------------------!
      real(kind=8) :: avg_intercepted   ! Mass flux
      real(kind=8) :: avg_qintercepted  ! Energy flux
      !----- Mass and energy input due to throughfall precipitation. ----------------------!
      real(kind=8) :: avg_throughfall   ! Mass flux
      real(kind=8) :: avg_qthroughfall  ! Energy flux
      !----- Sensible heat flux -----------------------------------------------------------!
      real(kind=8) :: avg_sensible_vc   ! Leaf      -> canopy air
      real(kind=8) :: avg_sensible_gc   ! Ground    -> canopy air
      real(kind=8) :: avg_sensible_ac   ! Free atm. -> canopy air
      real(kind=8) :: avg_heatstor_veg  ! Heat storage in vegetation
      !----- Carbon flux ------------------------------------------------------------------!
      real(kind=8) :: avg_carbon_ac     ! Free atm. -> canopy air
      !----- Soil fluxes ------------------------------------------------------------------!
      real(kind=8),pointer,dimension(:) :: avg_smoist_gg     ! Moisture flux between layers
      real(kind=8),pointer,dimension(:) :: avg_smoist_gc     ! Transpired soil moisture sink
      real(kind=8),pointer,dimension(:) :: avg_sensible_gg   ! Soil heat flux between layers
      real(kind=8)                      :: avg_drainage      ! Drainage at the bottom.
      real(kind=8)                      :: avg_drainage_heat ! Drainage at the bottom.
      !------------------------------------------------------------------------------------!
      !     Fast time flux variables for each time step.  These variables will be defined  !
      ! only when the user is debugging.                                                   !
      !------------------------------------------------------------------------------------!
      real(kind=8) :: flx_netrad        ! Net radiation
      !----- Water fluxes -----------------------------------------------------------------!
      real(kind=8) :: flx_vapor_vc      ! Leaf       -> canopy air:  evap./cond. flux
      real(kind=8) :: flx_dew_cg        ! Canopy     -> ground    :  condensation flux
      real(kind=8) :: flx_vapor_gc      ! Ground     -> canopy air:  evaporation flux
      real(kind=8) :: flx_vapor_ac      ! Free atm.  -> canopy air:  vapour flux
      real(kind=8) :: flx_transp        ! Transpiration
      real(kind=8) :: flx_evap          ! Evaporation
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
      real(kind=8) :: flx_sensible_vc   ! Leaf      -> canopy air
      real(kind=8) :: flx_sensible_gc   ! Ground    -> canopy air
      real(kind=8) :: flx_sensible_ac   ! Free atm. -> canopy air
      real(kind=8) :: flx_heatstor_veg  ! Heat storage in vegetation
      !----- Carbon flux ------------------------------------------------------------------!
      real(kind=8) :: flx_carbon_ac     ! Free atm. -> canopy air
      !----- Soil fluxes ------------------------------------------------------------------!
      real(kind=8),pointer,dimension(:) :: flx_smoist_gg     ! Moisture flux between layers
      real(kind=8),pointer,dimension(:) :: flx_smoist_gc     ! Transpired soil moisture sink
      real(kind=8),pointer,dimension(:) :: flx_sensible_gg   ! Soil heat flux between layers
      real(kind=8)                      :: flx_drainage      ! Drainage at the bottom.
      real(kind=8)                      :: flx_drainage_heat ! Drainage at the bottom.
      !----- Cohort-level fluxes. ---------------------------------------------------------!
      real(kind=8),pointer,dimension(:) :: cfx_hflxvc        ! Sensible heat
      real(kind=8),pointer,dimension(:) :: cfx_qwflxvc       ! Latent heat - Evaporation
      real(kind=8),pointer,dimension(:) :: cfx_qwshed        ! Int. en. of shed water
      real(kind=8),pointer,dimension(:) :: cfx_qtransp       ! Latent heat - Transpiration
      real(kind=8),pointer,dimension(:) :: cfx_qintercepted  ! Int. en. of intercept. H2O
      !----- Full budget variables --------------------------------------------------------!
      real(kind=8) :: co2budget_storage
      real(kind=8) :: co2budget_loss2atm
      real(kind=8) :: ebudget_storage
      real(kind=8) :: ebudget_loss2atm
      real(kind=8) :: ebudget_loss2drainage
      real(kind=8) :: ebudget_loss2runoff
      real(kind=8) :: wbudget_storage
      real(kind=8) :: wbudget_loss2atm
      real(kind=8) :: wbudget_loss2drainage
      real(kind=8) :: wbudget_loss2runoff
   end type rk4patchtype
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Structure with atmospheric and some other site-level data that is often used       !
   ! inside the Runge-Kutta integrator and do not change over time.                        !
   !---------------------------------------------------------------------------------------!
   type rk4sitetype
      integer                        :: lsl
      real(kind=8), dimension(n_pft) :: green_leaf_factor
      real(kind=8)                   :: atm_rhos
      real(kind=8)                   :: vels
      real(kind=8)                   :: atm_tmp
      real(kind=8)                   :: atm_theta
      real(kind=8)                   :: atm_theiv
      real(kind=8)                   :: atm_lntheta
      real(kind=8)                   :: atm_shv
      real(kind=8)                   :: atm_rvap
      real(kind=8)                   :: atm_rhv
      real(kind=8)                   :: atm_co2
      real(kind=8)                   :: zoff
      real(kind=8)                   :: atm_exner
      real(kind=8)                   :: pcpg
      real(kind=8)                   :: qpcpg
      real(kind=8)                   :: dpcpg
      real(kind=8)                   :: atm_prss
      real(kind=8)                   :: geoht
      real(kind=8)                   :: rshort
      real(kind=8)                   :: rlong
      real(kind=8)                   :: lon
      real(kind=8)                   :: lat
   end type rk4sitetype
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
   end type integration_vars
   !---------------------------------------------------------------------------------------!


   !----- This is the actual integration buffer structure. --------------------------------!
   type(integration_vars) :: integration_buff
   type(rk4sitetype)      :: rk4site
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    The following variable will be loaded from the user's namelist.                    !
   !---------------------------------------------------------------------------------------!
   integer                   :: ibranch_thermo ! This flag tells whether we consider the 
                                               !    branch actively affecting heat capacity
                                               !    and radiation interception. 
                                               ! 0 - no (default);
                                               ! 1 - yes (under development/test);

   real                      :: rk4_tolerance  ! The RK4 tolerance (or epsilon).  While 
                                               !    rk4eps is the actual variable used in
                                               !    Runge-Kutta (a double precision vari-
                                               !    able), rk4tol is the one given at the 
                                               !    namelist (a single precision variable).

   integer                   :: ipercol        ! This flag controls which percolation 
                                               !    scheme we should use
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
   !---------------------------------------------------------------------------------------!



   !----- Constants used in rk4_derivs ----------------------------------------------------!
   logical      :: supersat_ok    ! It is fine for evaporation and transpiration to
                                  !    occur even if this causes the canopy air to 
                                  !    be super-saturated                          [   T|F]
                                  !    (N.B. Super-saturation can occur even if 
                                  !     supersat_ok is .false., but in this case
                                  !     only mixing with free atmosphere can cause
                                  !     the super-saturation).
   logical      :: force_idealgas ! The integrator will adjust pressure every time 
                                  !    step, including the internal ones, to make 
                                  !    sure the ideal gas is respected.  If set to
                                  !    false, it will keep pressure constant 
                                  !    within on DTLSM time step, and not bother
                                  !    forcing the canopy air space to respect the
                                  !    ideal gas                                   [   T|F]
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
   !    These two parameter will scale the cohort heat capacity inside the RK4 integrator. !
   ! These are used only when the biophysics is set to be the same as in ED-2.1, and       !
   ! should not be used in standard simulations.                                           !
   !---------------------------------------------------------------------------------------!
   real(kind=8) :: hcapveg_ref         ! Reference value                          [ J/m³/K]
   real(kind=8) :: min_height          ! Minimum vegetation height                [      m]
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
   real(kind=8) :: rk4min_can_rvap      ! Minimum canopy    mixing ratio        [kg/kg_air]
   real(kind=8) :: rk4max_can_rvap      ! Maximum canopy    mixing ratio        [kg/kg_air]
   real(kind=8) :: rk4min_can_rhv       ! Minimum canopy    relative humidity   [      ---]
   real(kind=8) :: rk4max_can_rhv       ! Maximum canopy    relative humidity   [      ---]
   real(kind=8) :: rk4min_can_co2       ! Minimum canopy    CO2 mixing ratio    [ µmol/mol]
   real(kind=8) :: rk4max_can_co2       ! Maximum canopy    CO2 mixing ratio    [ µmol/mol]
   real(kind=8) :: rk4min_soil_temp     ! Minimum soil      temperature         [        K]
   real(kind=8) :: rk4max_soil_temp     ! Maximum soil      temperature         [        K]
   real(kind=8) :: rk4min_veg_temp      ! Minimum leaf      temperature         [        K]
   real(kind=8) :: rk4max_veg_temp      ! Maximum leaf      temperature         [        K]
   real(kind=8) :: rk4min_sfcw_temp     ! Minimum snow/pond temperature         [        K]
   real(kind=8) :: rk4max_sfcw_temp     ! Maximum snow/pond temperature         [        K]
   !----- The following variables have units different from the actual value. -------------!
   real(kind=8) :: rk4min_veg_lwater    ! Maximum leaf      surface water mass  [kg/m²leaf]
   real(kind=8) :: rk4min_sfcw_moist    ! Maximum snow/pond water mass          [m³/m³soil]
   real(kind=8) :: rk4min_virt_moist    ! Minimum virtual pool mass             [m³/m³soil]
   !----- The following variables will be defined in sfcdata_ed (ed_init.f90). ------------!
   real(kind=8) :: rk4min_sfcw_mass     ! Minimum snow/pond    mass             [    kg/m²]
   real(kind=8) :: rk4min_virt_water    ! Minimum virtual pool mass             [    kg/m²]
   !----- The following variables will be defined every time step. ------------------------!
   real(kind=8) :: rk4min_can_theta     ! Minimum canopy    potential temp.     [        K]
   real(kind=8) :: rk4max_can_theta     ! Maximum canopy    potential temp.     [        K]
   real(kind=8) :: rk4min_can_lntheta   ! Minimum canopy    log of theta        [      ---]
   real(kind=8) :: rk4max_can_lntheta   ! Maximum canopy    log of theta        [      ---]
   real(kind=8) :: rk4min_can_theiv     ! Minimum canopy    eq. pot. temp.      [        K]
   real(kind=8) :: rk4max_can_theiv     ! Maximum canopy    eq. pot. temp.      [        K]
   real(kind=8) :: rk4min_can_prss      ! Minimum canopy    pressure            [       Pa]
   real(kind=8) :: rk4max_can_prss      ! Maximum canopy    pressure            [       Pa]
   !---------------------------------------------------------------------------------------!


   !----- The following variables are double precision version of some bounds. ------------!
   real(kind=8) :: rk4tiny_sfcw_mass     ! Min. non-negligible snow/pond mass   [    kg/m²]
   real(kind=8) :: rk4water_stab_thresh  ! Min. mass for a layer to be stable   [    kg/m²]
   real(kind=8) :: rk4snowmin            ! Min. snow mass required for new lyr  [    kg/m²]
   real(kind=8) :: rk4dry_veg_lwater     ! Min. non-negligible leaf water mass  [kg/m²leaf]
   real(kind=8) :: rk4fullveg_lwater     ! Max. leaf water mass possible        [kg/m²leaf]
   real(kind=8) :: rk4tiny_sfcw_depth    ! Minimum snow/pond depth              [        m]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     The following variables are going to be allocated at the first time step, and     !
   ! they will be re-calculated every time the integrator is called.                       !
   ! rk4min_soil_water                  ! Minimum soil moisture                 [    m³/m³]!
   ! rk4max_soil_water                  ! Maximum soil moisture                 [    m³/m³]!
   !---------------------------------------------------------------------------------------!
   real(kind=8), dimension(:), allocatable :: rk4min_soil_water
   real(kind=8), dimension(:), allocatable :: rk4max_soil_water
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
   ! Runge-Kutta solver (checked at copy_patch_init).                                   !
   !---------------------------------------------------------------------------------------!
   logical :: toosparse
   
   !----- Flag to tell whether there is at least one "resolvable" cohort in this patch ----!
   logical :: any_resolvable

   !----- Canopy water and heat capacity variables. ---------------------------------------!
   real(kind=8)    :: zoveg
   real(kind=8)    :: zveg
   real(kind=8)    :: wcapcan
   real(kind=8)    :: wcapcani
   real(kind=8)    :: hcapcani
   real(kind=8)    :: ccapcani
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      Integrator error statistics.                                                     !
   !---------------------------------------------------------------------------------------!
   !----- Number of variables other than soil and surface that will be analysed. ----------!
   integer                          , parameter   :: nerrfix = 18

   !----- Total number of variables that will be analysed. --------------------------------!
   integer                                        :: nerr

   !----- Default file name to be given to the error files. -------------------------------!
   character(len=str_len)                         :: errmax_fout
   character(len=str_len)                         :: sanity_fout
   character(len=str_len)                         :: thbnds_fout
   character(len=str_len)                         :: detail_pref

   !----- The error counter and label. ----------------------------------------------------!
   integer(kind=8)  , dimension(:,:), allocatable :: integ_err
   character(len=13), dimension(:)  , allocatable :: integ_lab

   !----- Offset needed for each of the following variables. ------------------------------!
   integer                                        :: osow ! Soil water.
   integer                                        :: osoe ! Soil energy.
   integer                                        :: oswe ! Surface water energy.
   integer                                        :: oswm ! Surface water mass.
   !=======================================================================================!
   !=======================================================================================!


   contains



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

      allocate(y%soil_energy(nzg))
      allocate(y%soil_water(nzg))
      allocate(y%soil_fracliq(nzg))
      allocate(y%soil_tempk(nzg))
      allocate(y%available_liquid_water(nzg))
      allocate(y%extracted_water(nzg))
      allocate(y%psiplusz(nzg))
      allocate(y%soilair99(nzg))
      allocate(y%soilair01(nzg))
      allocate(y%soil_liq(nzg))
      allocate(y%soil_restz(nzg))

      allocate(y%sfcwater_energy(nzs))
      allocate(y%sfcwater_mass(nzs))
      allocate(y%sfcwater_depth(nzs))
      allocate(y%sfcwater_fracliq(nzs))
      allocate(y%sfcwater_tempk(nzs))
      allocate(y%sfcwater_restz(nzs))

      !------------------------------------------------------------------------------------!
      !     Diagnostics - for now we will always allocate the diagnostics, even if they    !
      !                   aren't used.                                                     !
      !------------------------------------------------------------------------------------!
      allocate(y%avg_sensible_gg(nzg))
      allocate(y%avg_smoist_gg(nzg))
      allocate(y%avg_smoist_gc(nzg))
      allocate(y%flx_sensible_gg(nzg))
      allocate(y%flx_smoist_gg(nzg))
      allocate(y%flx_smoist_gc(nzg))

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

      nullify(y%soil_energy)
      nullify(y%soil_water)
      nullify(y%soil_fracliq)
      nullify(y%soil_tempk)
      nullify(y%available_liquid_water)
      nullify(y%extracted_water)
      nullify(y%psiplusz)
      nullify(y%soilair99)
      nullify(y%soilair01)
      nullify(y%soil_liq)
      nullify(y%soil_restz)

      nullify(y%sfcwater_energy)
      nullify(y%sfcwater_mass)
      nullify(y%sfcwater_depth)
      nullify(y%sfcwater_fracliq)
      nullify(y%sfcwater_tempk)
      nullify(y%sfcwater_restz)

      !------------------------------------------------------------------------------------!
      !     Diagnostics - for now we will always allocate the diagnostics, even if they    !
      !                   aren't used.                                                     !
      !------------------------------------------------------------------------------------!
      nullify(y%avg_smoist_gg)
      nullify(y%avg_smoist_gc)
      nullify(y%avg_sensible_gg)
      nullify(y%flx_smoist_gg)
      nullify(y%flx_smoist_gc)
      nullify(y%flx_sensible_gg)

      return
   end subroutine nullify_rk4_patch
   !=======================================================================================!
   !=======================================================================================!






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
      y%ebudget_loss2atm               = 0.d0
      y%ebudget_loss2drainage          = 0.d0
      y%ebudget_loss2runoff            = 0.d0
      y%wbudget_storage                = 0.d0
      y%wbudget_loss2atm               = 0.d0
      y%wbudget_loss2drainage          = 0.d0
      y%wbudget_loss2runoff            = 0.d0
     
      y%can_temp                       = 0.d0
      y%can_rvap                       = 0.d0
      y%can_shv                        = 0.d0
      y%can_ssh                        = 0.d0
      y%can_rhv                        = 0.d0
      y%can_co2                        = 0.d0
      y%can_theta                      = 0.d0
      y%can_theiv                      = 0.d0
      y%can_lntheta                    = 0.d0
      y%can_depth                      = 0.d0
      y%can_rhos                       = 0.d0
      y%can_prss                       = 0.d0
      y%can_exner                      = 0.d0
      y%opencan_frac                   = 0.d0
      y%ggbare                         = 0.d0
      y%ggveg                          = 0.d0
      y%ggnet                          = 0.d0
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
     

      y%avg_carbon_ac                  = 0.d0
      y%avg_vapor_vc                   = 0.d0
      y%avg_dew_cg                     = 0.d0
      y%avg_vapor_gc                   = 0.d0
      y%avg_wshed_vg                   = 0.d0
      y%avg_intercepted                = 0.d0
      y%avg_throughfall                = 0.d0
      y%avg_vapor_ac                   = 0.d0
      y%avg_transp                     = 0.d0
      y%avg_evap                       = 0.d0
      y%avg_netrad                     = 0.d0
      y%avg_sensible_vc                = 0.d0
      y%avg_qwshed_vg                  = 0.d0
      y%avg_qintercepted               = 0.d0
      y%avg_qthroughfall               = 0.d0
      y%avg_sensible_gc                = 0.d0
      y%avg_sensible_ac                = 0.d0
      y%avg_heatstor_veg               = 0.d0

      y%avg_drainage                   = 0.d0
      y%avg_drainage_heat              = 0.d0

      y%flx_carbon_ac                  = 0.d0
      y%flx_vapor_vc                   = 0.d0
      y%flx_dew_cg                     = 0.d0
      y%flx_vapor_gc                   = 0.d0
      y%flx_wshed_vg                   = 0.d0
      y%flx_intercepted                = 0.d0
      y%flx_throughfall                = 0.d0
      y%flx_vapor_ac                   = 0.d0
      y%flx_transp                     = 0.d0
      y%flx_evap                       = 0.d0
      y%flx_netrad                     = 0.d0
      y%flx_sensible_vc                = 0.d0
      y%flx_qwshed_vg                  = 0.d0
      y%flx_qintercepted               = 0.d0
      y%flx_qthroughfall               = 0.d0
      y%flx_sensible_gc                = 0.d0
      y%flx_sensible_ac                = 0.d0
      y%flx_heatstor_veg               = 0.d0

      y%flx_drainage                   = 0.d0
      y%flx_drainage_heat              = 0.d0

      !----- The following variables are pointers, check whether they are linked or not. --!
      if(associated(y%soil_energy           ))   y%soil_energy(:)                 = 0.d0
      if(associated(y%soil_tempk            ))   y%soil_tempk(:)                  = 0.d0
      if(associated(y%soil_fracliq          ))   y%soil_fracliq(:)                = 0.d0
      if(associated(y%soil_water            ))   y%soil_water(:)                  = 0.d0
      if(associated(y%available_liquid_water))   y%available_liquid_water(:)      = 0.d0
      if(associated(y%extracted_water       ))   y%extracted_water(:)             = 0.d0
      if(associated(y%psiplusz              ))   y%psiplusz(:)                    = 0.d0
      if(associated(y%soilair99             ))   y%soilair99(:)                   = 0.d0
      if(associated(y%soilair01             ))   y%soilair01(:)                   = 0.d0
      if(associated(y%soil_liq              ))   y%soil_liq(:)                    = 0.d0
      if(associated(y%soil_restz            ))   y%soil_restz(:)                  = 0.d0
     
      if(associated(y%sfcwater_depth        ))   y%sfcwater_depth(:)              = 0.d0
      if(associated(y%sfcwater_mass         ))   y%sfcwater_mass(:)               = 0.d0
      if(associated(y%sfcwater_energy       ))   y%sfcwater_energy(:)             = 0.d0
      if(associated(y%sfcwater_tempk        ))   y%sfcwater_tempk(:)              = 0.d0
      if(associated(y%sfcwater_fracliq      ))   y%sfcwater_fracliq(:)            = 0.d0
      if(associated(y%sfcwater_restz        ))   y%sfcwater_restz(:)              = 0.d0

      if(associated(y%avg_smoist_gg         ))   y%avg_smoist_gg(:)               = 0.d0
      if(associated(y%avg_smoist_gc         ))   y%avg_smoist_gc(:)               = 0.d0
      if(associated(y%avg_sensible_gg       ))   y%avg_sensible_gg(:)             = 0.d0

      if(associated(y%flx_smoist_gg         ))   y%flx_smoist_gg(:)               = 0.d0
      if(associated(y%flx_smoist_gc         ))   y%flx_smoist_gc(:)               = 0.d0
      if(associated(y%flx_sensible_gg       ))   y%flx_sensible_gg(:)             = 0.d0

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

      if (associated(y%soil_energy))             deallocate(y%soil_energy)
      if (associated(y%soil_water))              deallocate(y%soil_water)
      if (associated(y%soil_fracliq))            deallocate(y%soil_fracliq)
      if (associated(y%soil_tempk))              deallocate(y%soil_tempk)
      if (associated(y%available_liquid_water))  deallocate(y%available_liquid_water)
      if (associated(y%extracted_water))         deallocate(y%extracted_water)
      if (associated(y%psiplusz))                deallocate(y%psiplusz)
      if (associated(y%soilair99))               deallocate(y%soilair99)
      if (associated(y%soilair01))               deallocate(y%soilair01)
      if (associated(y%soil_liq))                deallocate(y%soil_liq)
      if (associated(y%soil_restz))              deallocate(y%soil_restz)

      if (associated(y%sfcwater_energy))         deallocate(y%sfcwater_energy)
      if (associated(y%sfcwater_mass))           deallocate(y%sfcwater_mass)
      if (associated(y%sfcwater_depth))          deallocate(y%sfcwater_depth)
      if (associated(y%sfcwater_fracliq))        deallocate(y%sfcwater_fracliq)
      if (associated(y%sfcwater_tempk))          deallocate(y%sfcwater_tempk)
      if (associated(y%sfcwater_restz))          deallocate(y%sfcwater_restz)

      ! Diagnostics
      if (associated(y%avg_smoist_gg))           deallocate(y%avg_smoist_gg)
      if (associated(y%avg_smoist_gc))           deallocate(y%avg_smoist_gc)
      if (associated(y%avg_sensible_gg))         deallocate(y%avg_sensible_gg)
      if (associated(y%flx_smoist_gg))           deallocate(y%flx_smoist_gg)
      if (associated(y%flx_smoist_gc))           deallocate(y%flx_smoist_gc)
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

      allocate(y%veg_energy      (maxcohort))
      allocate(y%veg_water       (maxcohort))
      allocate(y%veg_temp        (maxcohort))
      allocate(y%veg_fliq        (maxcohort))
      allocate(y%hcapveg         (maxcohort))
      allocate(y%veg_wind        (maxcohort))
      allocate(y%veg_reynolds    (maxcohort))
      allocate(y%veg_grashof     (maxcohort))
      allocate(y%veg_nussfree    (maxcohort))
      allocate(y%veg_nussforc    (maxcohort))
      allocate(y%lint_shv        (maxcohort))
      allocate(y%nplant          (maxcohort))
      allocate(y%lai             (maxcohort))
      allocate(y%wai             (maxcohort))
      allocate(y%wpa             (maxcohort))
      allocate(y%tai             (maxcohort))
      allocate(y%crown_area      (maxcohort))
      allocate(y%gbh             (maxcohort))
      allocate(y%gbw             (maxcohort))
      allocate(y%gsw_open        (maxcohort))
      allocate(y%gsw_closed      (maxcohort))
      allocate(y%psi_open        (maxcohort))
      allocate(y%psi_closed      (maxcohort))
      allocate(y%fs_open         (maxcohort))
      allocate(y%resolvable      (maxcohort))
      allocate(y%gpp             (maxcohort))
      allocate(y%leaf_resp       (maxcohort))
      allocate(y%root_resp       (maxcohort))
      allocate(y%growth_resp     (maxcohort))
      allocate(y%storage_resp    (maxcohort))
      allocate(y%vleaf_resp      (maxcohort))
      allocate(y%rshort_v        (maxcohort))
      allocate(y%rlong_v         (maxcohort))
      allocate(y%cfx_hflxvc      (maxcohort))  
      allocate(y%cfx_qwflxvc     (maxcohort))  
      allocate(y%cfx_qwshed      (maxcohort))  
      allocate(y%cfx_qtransp     (maxcohort))  
      allocate(y%cfx_qintercepted(maxcohort))  

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
          
      nullify(y%veg_energy      )
      nullify(y%veg_water       )
      nullify(y%veg_temp        )
      nullify(y%veg_fliq        )
      nullify(y%hcapveg         )
      nullify(y%veg_wind        )
      nullify(y%veg_reynolds    )
      nullify(y%veg_grashof     )
      nullify(y%veg_nussfree    )
      nullify(y%veg_nussforc    )
      nullify(y%lint_shv        )
      nullify(y%nplant          )
      nullify(y%lai             )
      nullify(y%wai             )
      nullify(y%wpa             )
      nullify(y%tai             )
      nullify(y%crown_area      )
      nullify(y%gbh             )
      nullify(y%gbw             )
      nullify(y%gsw_open        )
      nullify(y%gsw_closed      )
      nullify(y%psi_open        )
      nullify(y%psi_closed      )
      nullify(y%fs_open         )
      nullify(y%resolvable      )
      nullify(y%gpp             )
      nullify(y%leaf_resp       )
      nullify(y%root_resp       )
      nullify(y%growth_resp     )
      nullify(y%storage_resp    )
      nullify(y%vleaf_resp      )
      nullify(y%rshort_v        )
      nullify(y%rlong_v         )
      nullify(y%cfx_hflxvc      )
      nullify(y%cfx_qwflxvc     )
      nullify(y%cfx_qwshed      )
      nullify(y%cfx_qtransp     )
      nullify(y%cfx_qintercepted)

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

      if(associated(y%veg_energy      )) y%veg_energy       = 0.d0
      if(associated(y%veg_water       )) y%veg_water        = 0.d0
      if(associated(y%veg_temp        )) y%veg_temp         = 0.d0
      if(associated(y%veg_fliq        )) y%veg_fliq         = 0.d0
      if(associated(y%hcapveg         )) y%hcapveg          = 0.d0
      if(associated(y%veg_wind        )) y%veg_wind         = 0.d0
      if(associated(y%veg_reynolds    )) y%veg_reynolds     = 0.d0
      if(associated(y%veg_grashof     )) y%veg_grashof      = 0.d0
      if(associated(y%veg_nussfree    )) y%veg_nussfree     = 0.d0
      if(associated(y%veg_nussforc    )) y%veg_nussforc     = 0.d0
      if(associated(y%lint_shv        )) y%lint_shv         = 0.d0
      if(associated(y%nplant          )) y%nplant           = 0.d0
      if(associated(y%lai             )) y%lai              = 0.d0
      if(associated(y%wai             )) y%wai              = 0.d0
      if(associated(y%wpa             )) y%wpa              = 0.d0
      if(associated(y%tai             )) y%tai              = 0.d0
      if(associated(y%crown_area      )) y%crown_area       = 0.d0
      if(associated(y%gbh             )) y%gbh              = 0.d0
      if(associated(y%gbw             )) y%gbw              = 0.d0
      if(associated(y%gsw_open        )) y%gsw_open         = 0.d0
      if(associated(y%gsw_closed      )) y%gsw_closed       = 0.d0
      if(associated(y%psi_open        )) y%psi_open         = 0.d0
      if(associated(y%psi_closed      )) y%psi_closed       = 0.d0
      if(associated(y%fs_open         )) y%fs_open          = 0.d0
      if(associated(y%resolvable      )) y%resolvable       = .false.
      if(associated(y%gpp             )) y%gpp              = 0.d0
      if(associated(y%leaf_resp       )) y%leaf_resp        = 0.d0
      if(associated(y%root_resp       )) y%root_resp        = 0.d0
      if(associated(y%growth_resp     )) y%growth_resp      = 0.d0
      if(associated(y%storage_resp    )) y%storage_resp     = 0.d0
      if(associated(y%vleaf_resp      )) y%vleaf_resp       = 0.d0
      if(associated(y%rshort_v        )) y%rshort_v         = 0.d0
      if(associated(y%rlong_v         )) y%rlong_v          = 0.d0
      if(associated(y%cfx_hflxvc      )) y%cfx_hflxvc       = 0.d0
      if(associated(y%cfx_qwflxvc     )) y%cfx_qwflxvc      = 0.d0
      if(associated(y%cfx_qwshed      )) y%cfx_qwshed       = 0.d0
      if(associated(y%cfx_qtransp     )) y%cfx_qtransp      = 0.d0
      if(associated(y%cfx_qintercepted)) y%cfx_qintercepted = 0.d0

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

      if(associated(y%veg_energy      )) deallocate(y%veg_energy      )
      if(associated(y%veg_water       )) deallocate(y%veg_water       )
      if(associated(y%veg_temp        )) deallocate(y%veg_temp        )
      if(associated(y%veg_fliq        )) deallocate(y%veg_fliq        )
      if(associated(y%hcapveg         )) deallocate(y%hcapveg         )
      if(associated(y%veg_wind        )) deallocate(y%veg_wind        )
      if(associated(y%veg_reynolds    )) deallocate(y%veg_reynolds    )
      if(associated(y%veg_grashof     )) deallocate(y%veg_grashof     )
      if(associated(y%veg_nussfree    )) deallocate(y%veg_nussfree    )
      if(associated(y%veg_nussforc    )) deallocate(y%veg_nussforc    )
      if(associated(y%lint_shv        )) deallocate(y%lint_shv        )
      if(associated(y%nplant          )) deallocate(y%nplant          )
      if(associated(y%lai             )) deallocate(y%lai             )
      if(associated(y%wai             )) deallocate(y%wai             )
      if(associated(y%wpa             )) deallocate(y%wpa             )
      if(associated(y%tai             )) deallocate(y%tai             )
      if(associated(y%crown_area      )) deallocate(y%crown_area      )
      if(associated(y%gbh             )) deallocate(y%gbh             )
      if(associated(y%gbw             )) deallocate(y%gbw             )
      if(associated(y%gsw_open        )) deallocate(y%gsw_open        )
      if(associated(y%gsw_closed      )) deallocate(y%gsw_closed      )
      if(associated(y%psi_open        )) deallocate(y%psi_open        )
      if(associated(y%psi_closed      )) deallocate(y%psi_closed      )
      if(associated(y%fs_open         )) deallocate(y%fs_open         )
      if(associated(y%resolvable      )) deallocate(y%resolvable      )
      if(associated(y%gpp             )) deallocate(y%gpp             )
      if(associated(y%leaf_resp       )) deallocate(y%leaf_resp       )
      if(associated(y%root_resp       )) deallocate(y%root_resp       )
      if(associated(y%growth_resp     )) deallocate(y%growth_resp     )
      if(associated(y%storage_resp    )) deallocate(y%storage_resp    )
      if(associated(y%vleaf_resp      )) deallocate(y%vleaf_resp      )
      if(associated(y%rshort_v        )) deallocate(y%rshort_v        )
      if(associated(y%rlong_v         )) deallocate(y%rlong_v         )
      if(associated(y%cfx_hflxvc      )) deallocate(y%cfx_hflxvc      )
      if(associated(y%cfx_qwflxvc     )) deallocate(y%cfx_qwflxvc     )
      if(associated(y%cfx_qwshed      )) deallocate(y%cfx_qwshed      )
      if(associated(y%cfx_qtransp     )) deallocate(y%cfx_qtransp     )
      if(associated(y%cfx_qintercepted)) deallocate(y%cfx_qintercepted)

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
      y%flx_vapor_vc                   = 0.d0
      y%flx_dew_cg                     = 0.d0
      y%flx_vapor_gc                   = 0.d0
      y%flx_wshed_vg                   = 0.d0
      y%flx_intercepted                = 0.d0
      y%flx_throughfall                = 0.d0
      y%flx_vapor_ac                   = 0.d0
      y%flx_transp                     = 0.d0
      y%flx_evap                       = 0.d0
      y%flx_netrad                     = 0.d0
      y%flx_sensible_vc                = 0.d0
      y%flx_qwshed_vg                  = 0.d0
      y%flx_qintercepted               = 0.d0
      y%flx_qthroughfall               = 0.d0
      y%flx_sensible_gc                = 0.d0
      y%flx_sensible_ac                = 0.d0
      y%flx_heatstor_veg               = 0.d0

      y%flx_drainage                   = 0.d0
      y%flx_drainage_heat              = 0.d0

      !----- Reset soil fluxes only when they are allocated. ------------------------------!
      if(associated(y%flx_smoist_gg   )) y%flx_smoist_gg    (:) = 0.d0
      if(associated(y%flx_smoist_gc   )) y%flx_smoist_gc    (:) = 0.d0
      if(associated(y%flx_sensible_gg )) y%flx_sensible_gg  (:) = 0.d0
      !----- Reset cohort-level energy fluxes when they are allocated. --------------------!
      if(associated(y%cfx_hflxvc      )) y%cfx_hflxvc       (:) = 0.d0
      if(associated(y%cfx_qwflxvc     )) y%cfx_qwflxvc      (:) = 0.d0
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
      y%flx_vapor_vc      = y%flx_vapor_vc      * hdidi
      y%flx_dew_cg        = y%flx_dew_cg        * hdidi
      y%flx_vapor_gc      = y%flx_vapor_gc      * hdidi
      y%flx_wshed_vg      = y%flx_wshed_vg      * hdidi
      y%flx_intercepted   = y%flx_intercepted   * hdidi
      y%flx_throughfall   = y%flx_throughfall   * hdidi
      y%flx_vapor_ac      = y%flx_vapor_ac      * hdidi
      y%flx_transp        = y%flx_transp        * hdidi
      y%flx_evap          = y%flx_evap          * hdidi
      y%flx_netrad        = y%flx_netrad        * hdidi
      y%flx_sensible_vc   = y%flx_sensible_vc   * hdidi
      y%flx_qwshed_vg     = y%flx_qwshed_vg     * hdidi
      y%flx_qintercepted  = y%flx_qintercepted  * hdidi
      y%flx_qthroughfall  = y%flx_qthroughfall  * hdidi
      y%flx_sensible_gc   = y%flx_sensible_gc   * hdidi
      y%flx_sensible_ac   = y%flx_sensible_ac   * hdidi
      y%flx_heatstor_veg  = y%flx_heatstor_veg  * hdidi

      y%flx_drainage      = y%flx_drainage      * hdidi
      y%flx_drainage_heat = y%flx_drainage_heat * hdidi

      if(associated(y%flx_smoist_gg   )) y%flx_smoist_gg    (:) = y%flx_smoist_gg   (:)    &
                                                                * hdidi
      if(associated(y%flx_smoist_gc   )) y%flx_smoist_gc    (:) = y%flx_smoist_gc   (:)    &
                                                                * hdidi
      if(associated(y%flx_sensible_gg )) y%flx_sensible_gg  (:) = y%flx_sensible_gg (:)    &
                                                                * hdidi
      !----- Cohort-level energy fluxes. --------------------------------------------------!
      if(associated(y%cfx_hflxvc      )) y%cfx_hflxvc       (:) = y%cfx_hflxvc      (:)    &
                                                                * hdidi
      if(associated(y%cfx_qwflxvc     )) y%cfx_qwflxvc      (:) = y%cfx_qwflxvc     (:)    &
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
           'CAN_THEIV    ','CAN_THETA    ','CAN_SHV      ','CAN_TEMP     ','CAN_PRSS     ' &
          ,'CAN_CO2      ','VEG_WATER    ','VEG_ENERGY   ','VIRT_HEAT    ','VIRT_WATER   ' &
          ,'CO2B_STORAGE ','CO2B_LOSS2ATM','EB_LOSS2ATM  ','WATB_LOSS2ATM','ENB_LOSS2DRA ' &
          ,'WATB_LOSS2DRA','ENB_STORAGE  ','WATB_STORAGE '/)
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
   subroutine find_derived_thbounds(lsl,mzg,can_rhos,can_theta,can_temp,can_shv,can_rvap   &
                                   ,can_prss,can_depth,nsoil)
      use grid_coms   , only : nzg           ! ! intent(in)
      use consts_coms , only : p008          & ! intent(in)
                             , rocp8         & ! intent(in)
                             , cp8           & ! intent(in)
                             , rdry8         & ! intent(in)
                             , epim18        & ! intent(in)
                             , ep8           & ! intent(in)
                             , mmdryi8       & ! intent(in)
                             , day_sec       & ! intent(in)
                             , hr_sec        & ! intent(in)
                             , min_sec       ! ! intent(in)
      use therm_lib8  , only : thetaeiv8     & ! function
                             , thetaeivs8    & ! function
                             , idealdenssh8  & ! function
                             , reducedpress8 & ! function
                             , eslif8        & ! function
                             , rslif8        ! ! function
      use soil_coms   , only : soil8         ! ! intent(in)
      use ed_misc_coms, only : current_time  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                     , intent(in) :: lsl
      integer                     , intent(in) :: mzg
      real(kind=8)                , intent(in) :: can_rhos
      real(kind=8)                , intent(in) :: can_theta
      real(kind=8)                , intent(in) :: can_temp
      real(kind=8)                , intent(in) :: can_shv
      real(kind=8)                , intent(in) :: can_rvap
      real(kind=8)                , intent(in) :: can_prss
      real(kind=8)                , intent(in) :: can_depth
      integer     , dimension(mzg), intent(in) :: nsoil
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                             :: can_prss_try
      real(kind=8)                             :: can_theta_try
      real(kind=8)                             :: can_theiv_try
      integer                                  :: k
      integer                                  :: hour
      integer                                  :: minute
      integer                                  :: second
      !----- Local parameters and locally saved variables. --------------------------------!
      logical                     , save       :: firsttime    = .true.
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     When we call this sub-routine for the first time, we must allocate the sanity  !
      ! check bounds for soil moisture.  Also, if the debugger of the sanity check is      !
      ! going to be used, we write the file header.                                        !
      !------------------------------------------------------------------------------------!
      if (firsttime) then

         allocate(rk4min_soil_water(mzg))
         allocate(rk4max_soil_water(mzg))

         if (print_thbnd) then
            open (unit=39,file=trim(thbnds_fout),status='replace',action='write')
            write(unit=39,fmt='(16(a,1x))')  '        YEAR','       MONTH','         DAY'  &
                                            ,'        HOUR','        MINU','        SECO'  &
                                            ,'    MIN_TEMP','    MAX_TEMP','     MIN_SHV'  &
                                            ,'     MAX_SHV','   MIN_THETA','   MAX_THETA'  &
                                            ,'   MIN_THEIV','   MAX_THEIV','    MIN_PRSS'  &
                                            ,'    MAX_PRSS'
            close(unit=39,status='keep')
         end if
         firsttime = .false.
      end if



      !------------------------------------------------------------------------------------!
      !     Find the bounds for pressure.  To avoid the pressure range to be too relaxed,  !
      ! switch one of the dependent variable a time, and use the current values for the    !
      ! other.  In addition, we force pressure to be bounded between the reduced pressure  !
      ! in case the reference height was different by the order of 10%.                    !
      !------------------------------------------------------------------------------------!
      !----- 1. Initial value, the most extreme one. --------------------------------------!
      rk4min_can_prss = reducedpress8(rk4site%atm_prss,rk4site%atm_theta,rk4site%atm_shv   &
                                     ,5.d-1*rk4site%geoht,can_theta,can_shv,can_depth)
      rk4max_can_prss = reducedpress8(rk4site%atm_prss,rk4site%atm_theta,rk4site%atm_shv   &
                                     ,1.1d0*rk4site%geoht,can_theta,can_shv,can_depth)
      !----- 2. Minimum temperature. ------------------------------------------------------!
      can_prss_try    = can_rhos * rdry8 * rk4min_can_temp * (1.d0 + epim18 * can_shv)
      rk4min_can_prss = min(rk4min_can_prss,can_prss_try)
      rk4max_can_prss = max(rk4max_can_prss,can_prss_try)
      !----- 3. Maximum temperature. ------------------------------------------------------!
      can_prss_try    = can_rhos * rdry8 * rk4max_can_temp * (1.d0 + epim18 * can_shv)
      rk4min_can_prss = min(rk4min_can_prss,can_prss_try)
      rk4max_can_prss = max(rk4max_can_prss,can_prss_try)
      !----- 4. Minimum specific humidity. ------------------------------------------------!
      can_prss_try    = can_rhos * rdry8 * can_temp * (1.d0 + epim18 * rk4min_can_shv)
      rk4min_can_prss = min(rk4min_can_prss,can_prss_try)
      rk4max_can_prss = max(rk4max_can_prss,can_prss_try)
      !----- 5. Maximum specific humidity. ------------------------------------------------!
      can_prss_try    = can_rhos * rdry8 * can_temp * (1.d0 + epim18 * rk4max_can_shv)
      rk4min_can_prss = min(rk4min_can_prss,can_prss_try)
      rk4max_can_prss = max(rk4max_can_prss,can_prss_try)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the bounds for potential temperature.  To avoid the pressure range to be  !
      ! too relaxed, switch one of the dependent variable a time, and use the current      !
      ! values for the other.                                                              !
      !------------------------------------------------------------------------------------!
      !----- 1. Initial value, the most extreme one. --------------------------------------!
      rk4min_can_theta   =  huge(1.d0)
      rk4max_can_theta   = -huge(1.d0)
      !----- 2. Minimum temperature. ------------------------------------------------------!
      can_theta_try      = rk4min_can_temp * (p008 / can_prss) ** rocp8
      rk4min_can_theta   = min(rk4min_can_theta,can_theta_try)
      rk4max_can_theta   = max(rk4max_can_theta,can_theta_try)
      !----- 3. Maximum temperature. ------------------------------------------------------!
      can_theta_try      = rk4max_can_temp * (p008 / can_prss) ** rocp8
      rk4min_can_theta   = min(rk4min_can_theta,can_theta_try)
      rk4max_can_theta   = max(rk4max_can_theta,can_theta_try)
      !----- 4. Minimum pressure. ---------------------------------------------------------!
      can_theta_try      = can_temp * (p008 / rk4min_can_prss) ** rocp8
      rk4min_can_theta   = min(rk4min_can_theta,can_theta_try)
      rk4max_can_theta   = max(rk4max_can_theta,can_theta_try)
      !----- 5. Maximum pressure. ---------------------------------------------------------!
      can_theta_try      = can_temp * (p008 / rk4max_can_prss) ** rocp8
      rk4min_can_theta   = min(rk4min_can_theta,can_theta_try)
      rk4max_can_theta   = max(rk4max_can_theta,can_theta_try)
      !----- 6. Find the logarithms. ------------------------------------------------------!
      rk4min_can_lntheta = log(rk4min_can_theta)
      rk4max_can_lntheta = log(rk4max_can_theta)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Minimum and maximum ice-vapour equivalent potential temperature.              !
      !------------------------------------------------------------------------------------!
      !----- 1. Initial value, the most extreme one. --------------------------------------!
      rk4min_can_theiv = rk4min_can_theta
      rk4max_can_theiv = -huge(1.d0)
      !----- 2. Maximum temperature. ------------------------------------------------------!
      can_theta_try    = rk4max_can_temp * (p008 / can_prss) ** rocp8
      can_theiv_try    = thetaeivs8(can_theta_try,rk4max_can_temp,can_rvap,0.d0,0.d0)
      rk4max_can_theiv = max(rk4max_can_theiv,can_theiv_try)
      !----- 3. Minimum pressure. ---------------------------------------------------------!
      can_theta_try    = can_temp * (p008 / rk4min_can_prss) ** rocp8
      can_theiv_try    = thetaeivs8(can_theta_try,can_temp,can_rvap,0.d0,0.d0)
      rk4max_can_theiv = max(rk4max_can_theiv,can_theiv_try)
      !----- 4. Maximum vapour mixing ratio. ----------------------------------------------!
      can_theta_try    = can_temp * (p008 / can_prss) ** rocp8
      can_theiv_try    = thetaeivs8(can_theta_try,can_temp,rk4max_can_rvap,0.d0,0.d0)
      rk4max_can_theiv = max(rk4max_can_theiv,can_theiv_try)
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
                              ,     rk4min_can_shv,   rk4min_can_theta,   rk4max_can_theta &
                              ,   rk4min_can_theiv,   rk4max_can_theiv,    rk4min_can_prss &
                              ,    rk4max_can_prss
         close (unit=39,status='keep')
      end if


      do k = lsl, mzg
         rk4min_soil_water(k) = soil8(nsoil(k))%soilcp * (1.d0 - rk4eps)
         rk4max_soil_water(k) = soil8(nsoil(k))%slmsts * (1.d0 + rk4eps)
      end do

      return
   end subroutine find_derived_thbounds
   !=======================================================================================!
   !=======================================================================================!
end module rk4_coms
!==========================================================================================!
!==========================================================================================!

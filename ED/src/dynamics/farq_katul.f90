!==========================================================================================!
!==========================================================================================!
! MODULE: PLANT_HYDRO
!> \brief Solve Farquhar Photosynthesis model together with Katul optimization-based
!> stomata conductance model.
!> \details The optimization framework is to maximize A - stoma_lambda * E every
!> model time step. Currently the optimization is acheived using regula falsi. The module only
!includes light-limited and RuBP-limited conditions\n
!> The references are:\n
!>       Xu et al. (2016) Diversity in plant hydraulic traits explains seasonal and
!> inter-annual variations in vegetation dynamics in seasonally dry tropical
!> forests.  New Phytologist\n
!>\n
!>       Manzoni, S., G. Vico, et al. (2011). "Optimizing stomatal conductance for maximum
!> carbon gain under water stress: a meta-analysis across plant functional types
!> and climates." Functional Ecology 25(3): 456-467.\n
!>\n
!>       Vico, G., S. Manzoni, et al. (2013). "A perspective on optimal leaf
!> stomatal conductance under CO2 and light co-limitations." Agricultural
!> and Forest Meteorology 182--183(0): 191-199.\n
!>\n
!>       Katul, G., S. Manzoni, et al. (2010). "A stomatal optimization theory to describe
!> the effects of atmospheric CO2 on leaf photosynthesis and transpiration."
!> Annals of Botany 105(3): 431-442.
!> \author Xiangtao Xu, 18 MAY 2020
!==========================================================================================!
!==========================================================================================!
module farq_katul

   contains
   !=======================================================================================!
   !=======================================================================================!
   ! SUBROUTINE KATUL_LPHYS       
   !> \brief   Main driver to calculate Farquhar-Katul photosynthesis system.
   !> Alternative to lphysio_full in farq_leuning.
   !> \details The realized Vcmax is reduced under low leaf water potential while leaf dark
   !> respiration keeps the same.
   !> \author Xiangtao Xu, 15 Feb. 2018
   !---------------------------------------------------------------------------------------!

   subroutine katul_lphys(ib,can_prss,can_rhos,can_shv,can_co2,ipft,leaf_par,leaf_temp     &
                         ,lint_shv,green_leaf_factor,leaf_aging_factor,vm_bar,rd_bar       &
                         ,leaf_gbw,leaf_psi,dmax_leaf_psi,A_open,A_closed,A_light,A_rubp,A_co2           &
                         ,gsw_open,gsw_closed,lsfc_shv_open,lsfc_shv_closed                &
                         ,lsfc_co2_open,lsfc_co2_closed,lint_co2_open,lint_co2_closed      &
                         ,leaf_resp,vmout,jmout,tpmout,jactout,comppout,limit_flag)
      
    use rk4_coms       , only : tiny_offset              & ! intent(in)
                              , effarea_transp           ! ! intent(in)
    use c34constants   , only : thispft                  & ! intent(out)
                              , met                      & ! intent(out)
                              , aparms                   ! ! intent(out)
    use farq_leuning   , only : comp_photo_tempfun       ! function
    use pft_coms       , only : photosyn_pathway         & ! intent(in)
                              , Vm0                      & ! intent(in)
                              , vm_low_temp              & ! intent(in)
                              , vm_high_temp             & ! intent(in)
                              , vm_hor                   & ! intent(in)
                              , vm_q10                   & ! intent(in)
                              , vm_decay_elow            & ! intent(in)
                              , vm_decay_ehigh           & ! intent(in)
                              , Jm0                      & ! intent(in)
                              , jm_low_temp              & ! intent(in)
                              , jm_high_temp             & ! intent(in)
                              , jm_hor                   & ! intent(in)
                              , jm_q10                   & ! intent(in)
                              , jm_decay_elow            & ! intent(in)
                              , jm_decay_ehigh           & ! intent(in)
                              , TPm0                     & ! intent(in)
                              , rd_low_temp              & ! intent(in)
                              , rd_high_temp             & ! intent(in)
                              , rd_hor                   & ! intent(in)
                              , rd_q10                   & ! intent(in)
                              , rd_decay_elow            & ! intent(in)
                              , rd_decay_ehigh           & ! intent(in)
                              , cuticular_cond           & ! intent(in)
                              , quantum_efficiency       & ! intent(in)
                              , curvpar_electron         & ! intent(in)
                              , qyield_psII              & ! intent(in)
                              , leaf_psi_tlp             & ! intent(in)
                              , stoma_lambda             & ! intent(in)
                              , stoma_beta               ! ! intent(in)
   use consts_coms    ,  only : mmh2oi8                  & ! intent(in)
                              , mmh2o8                   & ! intent(in)
                              , mmdryi8                  & ! intent(in)
                              , mmdry8                   & ! intent(in)
                              , ep8                      & ! intent(in)
                              , epi8                     & ! intent(in)
                              , t008                     & ! intent(in)
                              , umol_2_mol8              & ! intent(in)
                              , mol_2_umol8              & ! intent(in)
                              , Watts_2_Ein8             ! ! intent(in)
    use physiology_coms, only : gbw_2_gbc8               & ! intent(in)
                              , gsc_2_gsw8               & ! intent(in)
                              , h2o_plant_lim            & ! intent(in)
                              , o2_ref8                  ! ! intent(in)

    implicit none

      !------ Arguments. ------------------------------------------------------------------!
      integer     , intent(in)    :: ib                ! Multithread buffer
      real(kind=4), intent(in)    :: can_prss          ! Canopy air pressure    [       Pa]
      real(kind=4), intent(in)    :: can_rhos          ! Canopy air density     [    kg/m2]
      real(kind=4), intent(in)    :: can_shv           ! Canopy air sp. hum.    [    kg/kg]
      real(kind=4), intent(in)    :: can_co2           ! Canopy air CO2         [ umol/mol]
      integer     , intent(in)    :: ipft              ! Plant functional type  [      ---]
      real(kind=4), intent(in)    :: leaf_par          ! Absorbed PAR           [     W/m2]
      real(kind=4), intent(in)    :: leaf_temp         ! Leaf temperature       [        K]
      real(kind=4), intent(in)    :: lint_shv          ! Leaf interc. sp. hum.  [    kg/kg]
      real(kind=4), intent(in)    :: green_leaf_factor ! Frac. of on-allom. gr. [      ---]
      real(kind=4), intent(in)    :: leaf_aging_factor ! Ageing parameter       [      ---]
      real(kind=4), intent(in)    :: vm_bar            ! Average Vm function    [umol/m2/s]
      real(kind=4), intent(in)    :: rd_bar            ! Average Rd function    [umol/m2/s]
      real(kind=4), intent(in)    :: leaf_gbw          ! B.lyr. cnd. of H2O     [  kg/m2/s]
      real(kind=4), intent(in)    :: leaf_psi          ! leaf water potential   [        m]
      real(kind=4), intent(in)    :: dmax_leaf_psi     ! Daily maximum leaf water potential   [        m]
      real(kind=4), intent(out)   :: A_open            ! Photosyn. rate (op.)   [umol/m2/s]
      real(kind=4), intent(out)   :: A_closed          ! Photosyn. rate (cl.)   [umol/m2/s]
      real(kind=4), intent(out)   :: A_light           ! Photosyn. rate (light) [umol/m2/s]
      real(kind=4), intent(out)   :: A_rubp            ! Photosyn. rate (RuBP)  [umol/m2/s]
      real(kind=4), intent(out)   :: A_co2             ! Photosyn. rate (CO2)   [umol/m2/s]
      real(kind=4), intent(out)   :: gsw_open          ! St. cnd. of H2O  (op.) [  kg/m2/s]
      real(kind=4), intent(out)   :: gsw_closed        ! St. cnd. of H2O  (cl.) [  kg/m2/s]
      real(kind=4), intent(out)   :: lsfc_shv_open     ! Leaf sfc. sp.hum.(op.) [    kg/kg]
      real(kind=4), intent(out)   :: lsfc_shv_closed   ! Leaf sfc. sp.hum.(cl.) [    kg/kg]
      real(kind=4), intent(out)   :: lsfc_co2_open     ! Leaf sfc. CO2    (op.) [ umol/mol]
      real(kind=4), intent(out)   :: lsfc_co2_closed   ! Leaf sfc. CO2    (cl.) [ umol/mol]
      real(kind=4), intent(out)   :: lint_co2_open     ! Intercell. CO2   (op.) [ umol/mol]
      real(kind=4), intent(out)   :: lint_co2_closed   ! Intercell. CO2   (cl.) [ umol/mol]
      real(kind=4), intent(out)   :: leaf_resp         ! Leaf respiration rate  [umol/m2/s]
      real(kind=4), intent(out)   :: vmout             ! Max. Rubisco capacity  [umol/m2/s]
      real(kind=4), intent(out)   :: jmout             ! Max. electron transp.  [umol/m2/s]
      real(kind=4), intent(out)   :: tpmout            ! Max. triose phoshphate [umol/m2/s]
      real(kind=4), intent(out)   :: jactout           ! Act. electron transp.  [umol/m2/s]
      real(kind=4), intent(out)   :: comppout          ! GPP compensation point [ umol/mol]
      integer     , intent(out)   :: limit_flag        ! Photosyn. limit. flag  [      ---]

      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                :: f_plastic8        ! Plasticity factor
      real(kind=8)                :: lambda8           ! marginal water use efficiency
      real(kind=4)                :: water_stress_factor !factor to represent water stress on photosystems
      real(kind=8)                :: opt_ci            ! ci of the optial solution [mol CO2/mol air]
      real(kind=8)                :: opt_fc            ! optimal CO2 assimilation rate [molCO2/m2/s]
      real(kind=8)                :: opt_gsc           ! optimal stomatal conductance  [mol/m2/s]
      real(kind=8)                :: opt_ci_closed     ! ci under closed stomata   [mol CO2/mol air]
      real(kind=8)                :: opt_fc_closed     ! fc under closed stomata       [molCO2/m2/s]
      real(kind=8)                :: opt_gsc_closed    ! gsc under closed stomata      [mol/m2/s]
      real(kind=8)                :: opt_fc_rubp       ! RuBisco limited assimilation  [molCO2/m2/s]
      real(kind=8)                :: opt_fc_light      ! Light   limited assimilation  [molCO2/m2/s]
      real(kind=8)                :: opt_fc_3rd        ! TPU/CO2 limited assimilation  [molCO2/m2/s]
      !----- External function. -----------------------------------------------------------!
      real(kind=4)    , external  :: sngloff           ! Safe double -> single precision
      !------------------------------------------------------------------------------------!
      
      !----- Initialise limit_flag to night time value. -----------------------------------!
      limit_flag = 0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Copy the meteorological forcing to the "met" structure.  Notice that some      !
      ! variables go through unit conversions, and all variables are converted to double   !
      ! precision.                                                                         !
      !------------------------------------------------------------------------------------!
      !----- 1. Variables that remain with the same units. --------------------------------!
      met(ib)%leaf_temp    = dble(leaf_temp)
      met(ib)%can_rhos     = dble(can_rhos )
      met(ib)%can_prss     = dble(can_prss )
      met(ib)%can_o2       = o2_ref8
      !----- 2. Convert specific humidity to mol/mol. -------------------------------------!
      met(ib)%can_shv      = epi8 * dble(can_shv) 
      !----- 3. Convert CO2 to mol/mol. ---------------------------------------------------!
      met(ib)%can_co2      = dble(can_co2) * umol_2_mol8
      !----- 4. Convert W/m2 to mol/m2/s. -------------------------------------------------!
      met(ib)%par          = dble(leaf_par) * Watts_2_Ein8
      !------------------------------------------------------------------------------------!
      !  5. Intercellular specific humidity, which is assumed to be at saturation          !
      !     given the leaf temperature.  We convert it to mol/mol.                         !
      !------------------------------------------------------------------------------------!
      met(ib)%lint_shv    = epi8 * dble(lint_shv)
      !------------------------------------------------------------------------------------!
      !  6. Find the conductivities for water and carbon.  The input for water is in       !
      !     kg/m2/s, and here we convert to mol/m2/s.  The convertion coefficient from     !
      !     water to carbon dioxide comes from M09's equation B14.  Here we multiply by    !
      !     the effective area for transpiration, depending on whether the leaves of this  !
      !     plant functional type are hypo-stomatous, symmetrical, or amphistomatous.      !
      !------------------------------------------------------------------------------------!
      met(ib)%blyr_cond_h2o = dble(leaf_gbw)  * mmdryi8 * effarea_transp(ipft)
      met(ib)%blyr_cond_co2 = gbw_2_gbc8 * met(ib)%blyr_cond_h2o
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Load physiological parameters that are PFT-dependent to the thispft structure. !
      ! Convert all variables to mol and Kelvin, when needed.                              !
      !------------------------------------------------------------------------------------!
      thispft(ib)%photo_pathway    = photosyn_pathway(ipft)
      thispft(ib)%b                = dble(cuticular_cond(ipft)) * umol_2_mol8
      thispft(ib)%vm_low_temp      = dble(vm_low_temp(ipft))  + t008
      thispft(ib)%vm_high_temp     = dble(vm_high_temp(ipft)) + t008
      thispft(ib)%vm_hor           = dble(vm_hor(ipft))
      thispft(ib)%vm_q10           = dble(vm_q10(ipft))
      thispft(ib)%vm_decay_elow    = dble(vm_decay_elow(ipft))
      thispft(ib)%vm_decay_ehigh   = dble(vm_decay_ehigh(ipft))
      thispft(ib)%jm_low_temp      = dble(jm_low_temp(ipft))  + t008
      thispft(ib)%jm_high_temp     = dble(jm_high_temp(ipft)) + t008
      thispft(ib)%jm_hor           = dble(jm_hor(ipft))
      thispft(ib)%jm_q10           = dble(jm_q10(ipft))
      thispft(ib)%jm_decay_elow    = dble(jm_decay_elow(ipft))
      thispft(ib)%jm_decay_ehigh   = dble(jm_decay_ehigh(ipft))
      thispft(ib)%rd_low_temp      = dble(rd_low_temp(ipft))  + t008
      thispft(ib)%rd_high_temp     = dble(rd_high_temp(ipft)) + t008
      thispft(ib)%rd_hor           = dble(rd_hor(ipft))
      thispft(ib)%rd_q10           = dble(rd_q10(ipft))
      thispft(ib)%rd_decay_elow    = dble(rd_decay_elow(ipft))
      thispft(ib)%rd_decay_ehigh   = dble(rd_decay_ehigh(ipft))
      thispft(ib)%alpha0           = dble(quantum_efficiency(ipft))
      thispft(ib)%curvpar          = dble(curvpar_electron(ipft))
      thispft(ib)%phi_psII         = dble(qyield_psII     (ipft))
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Set Vm0 and find terms that typically depend upon Vm0 (Rd0, Jm0, TPm0).        !
      ! This may be the default parameters, but in case trait plasticity is enabled, they  !
      ! must be down-regulated.  Convert the resulting terms to mol/m2/s.                  !
      !     If no plasticity is applied, vm_bar and rd_bar will be qual to vm0 and rd0     !
      ! and f_plastic8 will be 1. Therefore it also works for this scenario.               !
      !     Jm0, and TPm0 are all scaled with the vm_bar:Vm0 ratio.                        !
      !------------------------------------------------------------------------------------!
      thispft(ib)%vm0  = dble(vm_bar) * umol_2_mol8
      thispft(ib)%rd0  = dble(rd_bar) * umol_2_mol8
      f_plastic8       = dble(vm_bar) / dble(vm0(ipft))
      thispft(ib)%jm0  = f_plastic8 * dble(jm0 (ipft)) * umol_2_mol8
      thispft(ib)%TPm0 = f_plastic8 * dble(TPm0(ipft)) * umol_2_mol8
      !------------------------------------------------------------------------------------!



      lambda8 = dble(stoma_lambda(ipft) * can_co2 / 400.)
      !stomata marginal water use efficiency
      !the co2 correcting factor is based on Manzoni et al. paper

      !update photosynthetic parameters with water stress
      select case (h2o_plant_lim)
      case (0,1,2,3)
          ! use fsw to account for water stress in photosyn_driv
          water_stress_factor = 1.
      case (4)
          ! leaf water potential will influence stomata optimization 
          ! at two different scales
          ! (1) [instantaneous] RUBP regeneration will be limited due 
          !     to low activity of ATP synthesas while the amount of 
          !     Rubisco does not change (Tezara et al. 2001 Science)
          ! (2) [daily or longer] stomatal sensitivity to water loss 
          !     would increase (Manzoni et al. 2011 Func. Ecol.)
            
          ! Jm0 will decrease by ~10% when leaf_psi
          ! is equal to leaf_psi_tlp [a test value, no real data to 
          ! parameterize it. However, leaf functions usually start to
          ! deteriorate after turgor is lost.
          water_stress_factor = max(1e-6,               &
                                    min(1.0,            &
                                        1. / (1. +      &
                                0.1 * (leaf_psi / leaf_psi_tlp(ipft)) ** 6.0)))
          lambda8 = lambda8 * dble(exp(stoma_beta(ipft) * dmax_leaf_psi))
      end select
         
      !thispft(ib)%vm0  = thispft(ib)%vm0 * dble(water_stress_factor)
      thispft(ib)%jm0  = thispft(ib)%jm0 * dble(water_stress_factor)
      !thispft(ib)%TPm0  = thispft(ib)%TPm0 * dble(water_stress_factor)
      thispft(ib)%alpha0  = thispft(ib)%alpha0 * dble(water_stress_factor)

      !------------------------------------------------------------------------------------!
      !     Compute the photosynthesis and leaf respiration parameters that depend on      !
      ! temperature, according to the parameters defined in the "thispft" structure and    !
      ! the functional form chosen by the user.  The variables that are defined there are: !
      ! - alpha     - the quantum yield, which may be a function of temperature.           !
      ! - Vm        - the maximum capacity of Rubisco to perform the carboxylase function. !
      ! - Jm        - the maximum electron transport rate.                                 !
      ! - J         - the actual electron transport rate.                                  !
      ! - TPm       - the maximum triose phosphate utilisation rate.                       !
      ! - leaf_resp - the leaf respiration.                                                !
      ! - compp     - the CO2 compensation point for gross photosynthesis (Gamma*)         !
      ! - kco2      - Michaelis-Mentel coefficient for CO2                                 !
      ! - ko2       - Michaelis-Mentel coefficient for O2.                                 !
      !------------------------------------------------------------------------------------!
      call comp_photo_tempfun(ib,leaf_aging_factor,green_leaf_factor)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! Perform Optimization
      !------------------------------------------------------------------------------------!
      call optimization_solver8(ib,ipft,lambda8,                                           &
                                opt_ci,opt_fc,opt_gsc,                                     &
                                opt_ci_closed,opt_fc_closed,opt_gsc_closed,                &
                                opt_fc_rubp,opt_fc_light,opt_fc_3rd,                       &
                                limit_flag)
      !------------------------------------------------------------------------------------!
      

      !------------------------------------------------------------------------------------!
      !    Copy the solution to the standard output variables.  Here we convert the values !
      ! back to the standard ED units                                                      !
      !------------------------------------------------------------------------------------!
      !----- Carbon demand, convert them to [umol/m2/s]. ----------------------------------!
      A_closed       = sngloff(opt_fc_closed              * mol_2_umol8 , tiny_offset)
      A_open         = sngloff(opt_fc                     * mol_2_umol8 , tiny_offset)
      A_light        = sngloff(opt_fc_light               * mol_2_umol8 , tiny_offset)
      A_rubp         = sngloff(opt_fc_rubp                * mol_2_umol8 , tiny_offset)
      A_co2          = sngloff(opt_fc_3rd                 * mol_2_umol8 , tiny_offset)
      !----- Stomatal resistance, convert the conductances to [kg/m2/s]. ------------------!
      gsw_closed     = sngloff(opt_gsc_closed * gsc_2_gsw8 * mmdry8 / effarea_transp(ipft)  &
                              , tiny_offset)
      gsw_open       = sngloff(opt_gsc * gsc_2_gsw8 * mmdry8 / effarea_transp(ipft)  &
                              , tiny_offset)
      !----- Leaf surface specific humidity, convert them to [kg/kg]. ---------------------!
      lsfc_shv_closed = sngloff((opt_gsc_closed * gsc_2_gsw8 * met(ib)%lint_shv    &
                                +met(ib)%blyr_cond_h2o * met(ib)%can_shv)   &
                                / (opt_gsc_closed * gsc_2_gsw8 + met(ib)%blyr_cond_h2o) * ep8         , tiny_offset)
      lsfc_shv_open = sngloff((opt_gsc * gsc_2_gsw8 * met(ib)%lint_shv    &
                                +met(ib)%blyr_cond_h2o * met(ib)%can_shv)   &
                                / (opt_gsc * gsc_2_gsw8 + met(ib)%blyr_cond_h2o) * ep8         , tiny_offset)
      !----- Leaf surface CO2 concentration, convert them to [umol/mol]. ------------------!
      lsfc_co2_closed = sngloff((met(ib)%can_co2 - opt_fc_closed / met(ib)%blyr_cond_co2)     * mol_2_umol8 , tiny_offset)
      lsfc_co2_open   = sngloff((met(ib)%can_co2 - opt_fc / met(ib)%blyr_cond_co2)     * mol_2_umol8 , tiny_offset)
      !----- Intercellular carbon dioxide concentration, convert them to [umol/mol]. ------!
      lint_co2_closed = sngloff(opt_ci_closed     * mol_2_umol8 , tiny_offset)
      lint_co2_open   = sngloff(opt_ci       * mol_2_umol8 , tiny_offset)
      !----- Leaf respiration [umol/m2/s]. ------------------------------------------------!
      leaf_resp       = sngloff(aparms(ib)%leaf_resp      * mol_2_umol8 , tiny_offset)
      !----- Maximum Rubisco capacity to perform the carboxylase function [umol/m2/s]. ----!
      vmout           = sngloff(aparms(ib)%vm             * mol_2_umol8 , tiny_offset)
      !----- Maximum electron transport rate [umol/m2/s]. ---------------------------------!
      jmout           = sngloff(aparms(ib)%jm             * mol_2_umol8 , tiny_offset)
      !----- Maximum triose phosphate utilisation rate [umol/m2/s]. -----------------------!
      tpmout          = sngloff(aparms(ib)%tpm            * mol_2_umol8 , tiny_offset)
      !----- Actual electron transport rate [umol/m2/s]. ----------------------------------!
      jactout         = sngloff(aparms(ib)%jact           * mol_2_umol8 , tiny_offset)
      !----- Gross photosynthesis compensation point, convert it to [umol/mol]. -----------!
      comppout        = sngloff(aparms(ib)%compp          * mol_2_umol8 , tiny_offset)
      !------------------------------------------------------------------------------------!
      return
   end subroutine katul_lphys
   !=======================================================================================!
   !=======================================================================================!



   !=======================================================================================!
   !=======================================================================================!
   ! SUBROUTINE OPTIMIZATION_SOLVER8
   !> \brief   Solver for the stomatal optimization problem
   !> \details Use Regular Falsi to search for the optimum
   !> \author Xiangtao Xu, 19 MAY 2018
   !---------------------------------------------------------------------------------------!
    subroutine optimization_solver8(ib,ipft,lambda,                                        &
                                    opt_ci,opt_fc,opt_gsc,                                 &
                                    opt_ci_closed,opt_fc_closed, opt_gsc_closed,           &
                                    opt_fc_rubp,opt_fc_light,opt_fc_3rd,                   &
                                    limit_flag)
    use c34constants,    only : met                      & ! intent(in)
                              , aparms                   & ! intent(in) 
                              , thispft                  ! ! intent(in)
    use farq_leuning,    only : find_twilight_min        ! ! function
    use physiology_coms, only : gsw_2_gsc8               ! ! intent(in)
    use consts_coms,     only : tiny_num8                ! ! intent(in)
    use ed_misc_coms   , only : current_time             ! ! intent(in)
    implicit none
        !------ Arguments. ------------------------------------------------------------------!
        integer     , intent(in)    :: ib           !! Multithread ID
        integer     , intent(in)    :: ipft         !! PFT for debugging purpose only
        real(kind=8), intent(in)    :: lambda       !! Marginal water use efficiency i.e. the langrangian multiplier in the optimization problem
        real(kind=8), intent(out)   :: opt_ci       !! Intercellular CO2 under optimized gsc [molCO2/molAir]
        real(kind=8), intent(out)   :: opt_fc       !! CO2 assimilation rate under opt gsc   [molCO2/m2/s]
        real(kind=8), intent(out)   :: opt_gsc      !! The optimal gsc                       [mol/m2/s]
        real(kind=8), intent(out)   :: opt_ci_closed!! Intercellular CO2 under closed gsc [molCO2/molAir]
        real(kind=8), intent(out)   :: opt_gsc_closed !! gsc under closed stomata [mol/m2/s]
        real(kind=8), intent(out)   :: opt_fc_closed!! Assimilation rate under closed stomata[molCO2/m2/s]
        real(kind=8), intent(out)   :: opt_fc_rubp  !! Rubisco limited assimilation rate     [molCO2/m2/s]
        real(kind=8), intent(out)   :: opt_fc_light !! light limited assimilation rate       [molCO2/m2/s]
        real(kind=8), intent(out)   :: opt_fc_3rd   !! TPU/CO2 limited assimilation rate     [molCO2/m2/s]
        integer,      intent(out)   :: limit_flag   !! Flag for photosynthesislimitation

        !------ Local Variables  ------------------------------------------------------------!
        real(kind=8), parameter     :: gsc_max = 1.5  ! maximum gsc allowed
        real(kind=8)                :: gsc_min        ! minimum gsc allowed
        real(kind=8)                :: par_twilight_min ! Minimum daytime radiation [mol/m2/s]
        logical                     :: is_resolvable
        real(kind=8)                :: rfx_lower      ! lower boundary X of regula falsi
        real(kind=8)                :: rfx_upper      ! upper boundary X of regula falsi
        real(kind=8)                :: rfy_lower      ! lower boundary Y of regula falsi
        real(kind=8)                :: rfy_upper      ! upper boundary Y of regula falsi
        real(kind=8)                :: rfx_new        ! new boundary X of regula falsi
        real(kind=8)                :: rfy_new        ! new boundary Y of regula falsi
        integer                     :: rf_side        ! side of current x relative to the root in regula falsi
        integer                     :: iter           ! iteration index
        integer, parameter          :: iter_max = 600 ! Maximum number of iteration
        real(kind=8), parameter     :: dg_min = 1.d-4 ! Tolerance of difference
        real(kind=8)                :: test_gsc       
        real(kind=8)                :: test_fc
        real(kind=8)                :: test_fe
        real(kind=8)                :: test_ci
        real(kind=8)                :: test_dcidg
        real(kind=8)                :: test_dfcdg
        real(kind=8)                :: test_dfedg
        real(kind=8)                :: test_fc_light
        real(kind=8)                :: test_fc_rubp
        real(kind=8)                :: test_fc_3rd
        real(kind=8)                :: opt_ci_light
        real(kind=8)                :: opt_ci_rubp
        real(kind=8)                :: opt_ci_3rd
        real(kind=8)                :: opt_gsc_light
        real(kind=8)                :: opt_gsc_rubp
        real(kind=8)                :: opt_gsc_3rd
        character(len=256)          :: limit_case

        ! for debugging purposes
        integer                     :: k
        logical, parameter          :: debug_flag = .false.
        !------------------------------------------------------------------------------------!
        gsc_min = thispft(ib)%b * gsw_2_gsc8 
  
  
     
        ! There are extreme cases when blyr_cond_co2 is too small or light is too low
        ! In this case, we do not solve optimization and choose to close stomata
        par_twilight_min = find_twilight_min(ib)
        is_resolvable = (met(ib)%blyr_cond_co2 > tiny_num8) .and. (met(ib)%par >= par_twilight_min)
        if (is_resolvable) then
            ! If resolvable use Regular Falsi to solve the optimization problem
            ! the purpose is to find a root for dfcdg - lambda * dfedg = 0

            ! loop over different photosyn_pathways

            !=======================!
            ! First RUBP limitation
            !======================!
            limit_case='RUBP'
            ! initial range of gsc is cuticular_gsc and gsc_max
            rfx_lower = gsc_min / 2.d0 ! a very small value, half of cuticular conductance
            rfx_upper = gsc_max * 20.d0 ! a very large value
  
            ! calculate the y values for rfx_lower
            call photosynthesis_stomata_solver8(ib,rfx_lower,limit_case,                          &
                                                test_ci,test_fc,test_fe,                          &
                                                test_dcidg,test_dfcdg,test_dfedg)
            rfy_lower = test_dfcdg - lambda * test_dfedg 
  
  
            ! do it again for rfx_upper
            call photosynthesis_stomata_solver8(ib,rfx_upper,limit_case,                          &
                                                test_ci,test_fc,test_fe,                          &
                                                test_dcidg,test_dfcdg,test_dfedg)
            rfy_upper = test_dfcdg - lambda * test_dfedg 
  
  
            ! Start iteration
            iter = 0
            ! check whether the y values have the same sign
            if (rfy_lower * rfy_upper >= 0.d0) then
                ! In this case, there is no root within the given range
                ! if rfy_lower is positive, we take the value of rfx_upper
                ! else we take the value of rfx_lower
                if  (rfy_lower > 0.) then
                    test_gsc = rfx_upper
                else
                    test_gsc = rfx_lower
                endif
            else
                ! There is at least one root
                ! Run regula falsi
                rf_side = 0 !
                do iter = 1, iter_max
                    ! exit condition
                    if (abs(rfx_lower - rfx_upper) .le. dg_min) then
                        exit
                    endif
  
                    ! update rfx and rfy with Illinois Method
                    rfx_new = (rfx_lower * rfy_upper - rfx_upper * rfy_lower) / (rfy_upper - rfy_lower)
                    call photosynthesis_stomata_solver8(ib,rfx_new,limit_case,                            &
                                                        test_ci,test_fc,test_fe,                          &
                                                        test_dcidg,test_dfcdg,test_dfedg)

                    rfy_new = test_dfcdg - lambda * test_dfedg
  
                    if (rfy_new * rfy_lower > 0.d0) then
                        ! the new point has the same sign as the lower
                        ! update the lower
                        rfx_lower = rfx_new
                        rfy_lower = rfy_new
  
                        ! Illinois Method, improve efficiency
                        if (rf_side == -1) then
                            rfy_upper = rfy_upper / 2.d0
                        endif
  
                        rf_side = -1
                    elseif (rfy_new * rfy_upper > 0.d0) then
                        ! the new point has the same sign as the upper
                        ! update the lower
                        rfx_upper = rfx_new
                        rfy_upper = rfy_new
  
                        ! Illinois Method, improve efficiency
                        if (rf_side == 1) then
                            rfy_lower = rfy_lower / 2.d0
                        endif
  
                        rf_side = 1
  
                    else
                        ! numerically they are the same
                        exit
                    endif
                enddo
  
                test_gsc = (rfx_lower + rfx_upper) / 2.d0
            endif
  
            ! final gsc should be bounded within gsc_min and gsc_max
            test_gsc = max(gsc_min,min(test_gsc,gsc_max))
  
            ! calculate the realized fc, ci
            call photosynthesis_stomata_solver8(ib,test_gsc,limit_case,                           &
                                                test_ci,test_fc,test_fe,                          &
                                                test_dcidg,test_dfcdg,test_dfedg)

            ! save the optimal gsc and carbon fluxes
            opt_fc_rubp = test_fc
            opt_gsc_rubp = test_gsc
            opt_ci_rubp = test_ci
  
            !=======================!
            ! Second LIGHT limitation
            !======================!
            limit_case='LIGHT'
            ! initial range of gsc is cuticular_gsc and gsc_max
            rfx_lower = gsc_min / 2.d0 ! a very small value, half of cuticular conductance
            rfx_upper = gsc_max * 20.d0 ! a very large value
  
            ! calculate the y values for rfx_lower
            call photosynthesis_stomata_solver8(ib,rfx_lower,limit_case,                          &
                                                test_ci,test_fc,test_fe,                          &
                                                test_dcidg,test_dfcdg,test_dfedg)
            rfy_lower = test_dfcdg - lambda * test_dfedg 
  
  
            ! do it again for rfx_upper
            call photosynthesis_stomata_solver8(ib,rfx_upper,limit_case,                          &
                                                test_ci,test_fc,test_fe,                          &
                                                test_dcidg,test_dfcdg,test_dfedg)
            rfy_upper = test_dfcdg - lambda * test_dfedg 
  
  
            ! Start iteration
            iter = 0
            ! check whether the y values have the same sign
            if (rfy_lower * rfy_upper >= 0.d0) then
                ! In this case, there is no root within the given range
                ! if rfy_lower is positive, we take the value of rfx_upper
                ! else we take the value of rfx_lower
                if  (rfy_lower > 0.) then
                    test_gsc = rfx_upper
                else
                    test_gsc = rfx_lower
                endif
            else
                ! There is at least one root
                ! Run regula falsi
                rf_side = 0 !
                do iter = 1, iter_max
                    ! exit condition
                    if (abs(rfx_lower - rfx_upper) .le. dg_min) then
                        exit
                    endif
  
                    ! update rfx and rfy with Illinois Method
                    rfx_new = (rfx_lower * rfy_upper - rfx_upper * rfy_lower) / (rfy_upper - rfy_lower)
                    call photosynthesis_stomata_solver8(ib,rfx_new,limit_case,                            &
                                                        test_ci,test_fc,test_fe,                          &
                                                        test_dcidg,test_dfcdg,test_dfedg)

                    rfy_new = test_dfcdg - lambda * test_dfedg
  
                    if (rfy_new * rfy_lower > 0.d0) then
                        ! the new point has the same sign as the lower
                        ! update the lower
                        rfx_lower = rfx_new
                        rfy_lower = rfy_new
  
                        ! Illinois Method, improve efficiency
                        if (rf_side == -1) then
                            rfy_upper = rfy_upper / 2.d0
                        endif
  
                        rf_side = -1
                    elseif (rfy_new * rfy_upper > 0.d0) then
                        ! the new point has the same sign as the upper
                        ! update the lower
                        rfx_upper = rfx_new
                        rfy_upper = rfy_new
  
                        ! Illinois Method, improve efficiency
                        if (rf_side == 1) then
                            rfy_lower = rfy_lower / 2.d0
                        endif
  
                        rf_side = 1
  
                    else
                        ! numerically they are the same
                        exit
                    endif
                enddo
  
                test_gsc = (rfx_lower + rfx_upper) / 2.d0
            endif
  
            ! final gsc should be bounded within gsc_min and gsc_max
            test_gsc = max(gsc_min,min(test_gsc,gsc_max))
  
            ! calculate the realized fc, ci
            call photosynthesis_stomata_solver8(ib,test_gsc,limit_case,                           &
                                                test_ci,test_fc,test_fe,                          &
                                                test_dcidg,test_dfcdg,test_dfedg)

            ! save the optimal gsc and carbon fluxes
            opt_fc_light = test_fc
            opt_gsc_light = test_gsc
            opt_ci_light = test_ci
 
            !=======================!
            ! third TPU or CO2 limitation
            !======================!
            select case (thispft(ib)%photo_pathway)
            case (3)
                limit_case='TPU'
            case (4)
                limit_case='CO2'
            end select
            ! initial range of gsc is cuticular_gsc and gsc_max
            rfx_lower = gsc_min / 2.d0 ! a very small value, half of cuticular conductance
            rfx_upper = gsc_max * 20.d0 ! a very large value
  
            ! calculate the y values for rfx_lower
            call photosynthesis_stomata_solver8(ib,rfx_lower,limit_case,                          &
                                                test_ci,test_fc,test_fe,                          &
                                                test_dcidg,test_dfcdg,test_dfedg)
            rfy_lower = test_dfcdg - lambda * test_dfedg 
  
  
            ! do it again for rfx_upper
            call photosynthesis_stomata_solver8(ib,rfx_upper,limit_case,                          &
                                                test_ci,test_fc,test_fe,                          &
                                                test_dcidg,test_dfcdg,test_dfedg)
            rfy_upper = test_dfcdg - lambda * test_dfedg 
  
  
            ! Start iteration
            iter = 0
            ! check whether the y values have the same sign
            if (rfy_lower * rfy_upper >= 0.d0) then
                ! In this case, there is no root within the given range
                ! if rfy_lower is positive, we take the value of rfx_upper
                ! else we take the value of rfx_lower
                if  (rfy_lower > 0.) then
                    test_gsc = rfx_upper
                else
                    test_gsc = rfx_lower
                endif
            else
                ! There is at least one root
                ! Run regula falsi
                rf_side = 0 !
                do iter = 1, iter_max
                    ! exit condition
                    if (abs(rfx_lower - rfx_upper) .le. dg_min) then
                        exit
                    endif
  
                    ! update rfx and rfy with Illinois Method
                    rfx_new = (rfx_lower * rfy_upper - rfx_upper * rfy_lower) / (rfy_upper - rfy_lower)
                    call photosynthesis_stomata_solver8(ib,rfx_new,limit_case,                            &
                                                        test_ci,test_fc,test_fe,                          &
                                                        test_dcidg,test_dfcdg,test_dfedg)

                    rfy_new = test_dfcdg - lambda * test_dfedg
  
                    if (rfy_new * rfy_lower > 0.d0) then
                        ! the new point has the same sign as the lower
                        ! update the lower
                        rfx_lower = rfx_new
                        rfy_lower = rfy_new
  
                        ! Illinois Method, improve efficiency
                        if (rf_side == -1) then
                            rfy_upper = rfy_upper / 2.d0
                        endif
  
                        rf_side = -1
                    elseif (rfy_new * rfy_upper > 0.d0) then
                        ! the new point has the same sign as the upper
                        ! update the lower
                        rfx_upper = rfx_new
                        rfy_upper = rfy_new
  
                        ! Illinois Method, improve efficiency
                        if (rf_side == 1) then
                            rfy_lower = rfy_lower / 2.d0
                        endif
  
                        rf_side = 1
  
                    else
                        ! numerically they are the same
                        exit
                    endif
                enddo
  
                test_gsc = (rfx_lower + rfx_upper) / 2.d0
            endif
  
            ! final gsc should be bounded within gsc_min and gsc_max
            test_gsc = max(gsc_min,min(test_gsc,gsc_max))
  
            ! calculate the realized fc, ci
            call photosynthesis_stomata_solver8(ib,test_gsc,limit_case,                           &
                                                test_ci,test_fc,test_fe,                          &
                                                test_dcidg,test_dfcdg,test_dfedg)

            ! save the optimal gsc and carbon fluxes
            opt_fc_3rd = test_fc
            opt_gsc_3rd = test_gsc
            opt_ci_3rd = test_ci
 
            !------------------------------------------------------------------------------------!
            ! Compare fc and set limit flag
            !------------------------------------------------------------------------------------!
            if ((opt_fc_light <= opt_fc_rubp) .and. &
                (opt_fc_light <= opt_fc_3rd)) then
                ! light-limited
                opt_fc = opt_fc_light
                opt_ci = opt_ci_light
                opt_gsc = opt_gsc_light
                limit_flag = 1
            elseif (opt_fc_rubp <= opt_fc_3rd) then
                ! RubisCO limited
                opt_fc = opt_fc_rubp
                opt_ci = opt_ci_rubp
                opt_gsc = opt_gsc_rubp
                limit_flag = 2
            else
                ! Triose Phosphate Utilisation or CO2 limited
                opt_fc = opt_fc_3rd
                opt_ci = opt_ci_3rd
                opt_gsc = opt_gsc_3rd
                limit_flag = 3
            endif
            !------------------------------------------------------------------------------------!

       
            opt_fc_closed = -aparms(ib)%leaf_resp
            opt_gsc_closed = gsc_min
            opt_ci_closed = met(ib)%can_co2 - opt_fc_closed /                   &
                            (opt_gsc_closed * met(ib)%blyr_cond_co2) *          &
                            (opt_gsc_closed + met(ib)%blyr_cond_co2)
        else  ! not resolvable
            limit_flag   = 0 ! night-time limitation
            opt_gsc      = gsc_min ! cuticular conductance
            opt_fc       = -aparms(ib)%leaf_resp
            opt_ci       = met(ib)%can_co2 + aparms(ib)%leaf_resp / opt_gsc
            opt_fc_rubp  = -aparms(ib)%leaf_resp
            opt_fc_light = -aparms(ib)%leaf_resp
            opt_fc_3rd   = -aparms(ib)%leaf_resp
            opt_fc_closed = -aparms(ib)%leaf_resp
            opt_gsc_closed = gsc_min
            opt_ci_closed = met(ib)%can_co2 - opt_fc_closed /                   &
                            (opt_gsc_closed * met(ib)%blyr_cond_co2) *          &
                            (opt_gsc_closed + met(ib)%blyr_cond_co2)
        endif
  

        
        if (debug_flag) then
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'Katul Stomatal Scheme Quality Check:'
            write (unit=*,fmt='(a,1x,i9)')   ' + HOUR:                ',current_time%hour
            write (unit=*,fmt='(a,1x,i9)')   ' + PFT:                ',ipft
            write (unit=*,fmt='(a,1x,es12.4)')   ' + Vm:                  ',aparms(ib)%vm
            write (unit=*,fmt='(a,1x,es12.4)')   ' + Jm:                  ',aparms(ib)%jm
            write (unit=*,fmt='(a,1x,es12.4)')   ' + Tpm:                 ',aparms(ib)%tpm
            write (unit=*,fmt='(a,1x,es12.4)')   ' + PAR:                 ',met(ib)%par
            write (unit=*,fmt='(a,1x,es12.4)')   ' + Leaf_temp:           ',met(ib)%leaf_temp
            write (unit=*,fmt='(a,1x,es12.4)')   ' + fc_rubp:             ',opt_fc_rubp
            write (unit=*,fmt='(a,1x,es12.4)')   ' + fc_light:            ',opt_fc_light
            write (unit=*,fmt='(a,1x,es12.4)')   ' + fc_3rd:              ',opt_fc_3rd
            write (unit=*,fmt='(a,1x,es12.4)')   ' + gsc_open:            ',opt_gsc    
            write (unit=*,fmt='(a,1x,es12.4)')   ' + ci:                  ',opt_ci     
            write (unit=*,fmt='(a,1x,es12.4)')   ' + lambda:              ',lambda
            write (unit=*,fmt='(a,1x,i9)')       ' + LIMIT_FLAG:          ',limit_flag
            write (unit=*,fmt='(a,1x,es12.4)')   ' + rfy_lower:       ',rfy_lower
            write (unit=*,fmt='(a,1x,es12.4)')   ' + rfy_upper:       ',rfy_upper
        endif

        return

    end subroutine optimization_solver8
   !=======================================================================================!
   !=======================================================================================!
      
   !=======================================================================================!
   !=======================================================================================!
   ! SUBROUTINE PHOTOSYNTHESIS_STOMATA_SOLVER8
   !> \brief   Solver for photosynthesis and its derivatives wrt. gsc
   !> \details Both C3 and C4 have three cases of limitation:\n
   !> C3: Light, RuBisCO, TPU\n
   !> C4: Light, RuBisCO, CO2\n
   !> \author Xiangtao Xu, 19 MAY 2018
   !---------------------------------------------------------------------------------------!
    subroutine photosynthesis_stomata_solver8(ib,gsc,limit_case,                            &
                                              ci,fc,fe,dcidg,dfcdg,dfedg)
    use c34constants,    only : met                      & ! intent(in)
                              , aparms                   & ! intent(in) 
                              , thispft                  ! ! intent(in)
    use physiology_coms, only : gsc_2_gsw8               & ! intent(in)
                              , klowco28                 & ! intent(in)
                              , iphysiol                 ! ! intent(in)
    use consts_coms,     only : tiny_num8                ! ! intent(in)
    implicit none
        !------ Arguments. ------------------------------------------------------------------!
        integer     , intent(in)        :: ib           !! Multithread ID
        real(kind=8), intent(in)        :: gsc          !! input stomatal conductance for CO2,   [mol/m2/s]
        character(len=*), intent(in)    :: limit_case   !! flag telling which case we are solving
        real(kind=8), intent(out)   :: ci           !! Intercellular CO2 under gsc           [molCO2/molAir]
        real(kind=8), intent(out)   :: fc           !! Realized assimilation rate under gsc  [molCO2/m2/s]
        real(kind=8), intent(out)   :: fe           !! transpiration                         [molH2O/m2/s]
        real(kind=8), intent(out)   :: dcidg        !! derivative of ci wrt. gsc
        real(kind=8), intent(out)   :: dfcdg        !! derivative of fc wrt. gsc
        real(kind=8), intent(out)   :: dfedg        !! derivative of fe wrt. gsc

        !------ Local Variables  ------------------------------------------------------------!
        real(kind=8)                :: gsbc         !! CO2 conductance from canopy air space to leaf (stomata + boundary layer)
        real(kind=8)                :: gsbw         !! H2O conductance from canopy air space to leaf (stomata + boundary layer)
        real(kind=8)                :: k1,k2        !! Variable used in photosynthesis equation
        real(kind=8)                :: a,b,c        !! Coefficients of the quadratic equation to solve ci
        real(kind=8)                :: rad          !! sqrt(b2-4ac)
        real(kind=8)                :: dbdg,dcdg    !! derivatives of b,c wrt. gsc
        real(kind=8)                :: ci_rubp      !! ci for rubp-limited scenario
        real(kind=8)                :: dcidg_rubp   !! derivative of ci wrt. gsc for rubp-limited scenario
        real(kind=8)                :: dfcdg_rubp   !! derivative of fc wrt. gsc for rubp-limited scenario
        real(kind=8)                :: ci_light     !! ci for light-limited scenario
        real(kind=8)                :: dcidg_light  !! derivative of ci wrt. gsc for light-limited scenario
        real(kind=8)                :: dfcdg_light  !! derivative of fc wrt. gsc for light-limited scenario
        real(kind=8)                :: ci_3rd       !! ci for TPU/CO2-limited scenario
        real(kind=8)                :: dcidg_3rd    !! derivative of ci wrt. gsc for TPU/CO2-limited scenario
        real(kind=8)                :: dfcdg_3rd    !! derivative of fc wrt. gsc for TPU/CO2-limited scenario

        !------------------------------------------------------------------------------------!


        !------------------------------------------------------------------------------------!
        ! First calculate fc, ci, and their derivatives
        ! A = k1 * (ci - compp) / (ci + k2) - leaf_resp
        ! ci ** 2 + b * ci + c = 0.
        !------------------------------------------------------------------------------------!
        gsbc = gsc * met(ib)%blyr_cond_co2 / (gsc + met(ib)%blyr_cond_co2) 
        gsbw = gsc * gsc_2_gsw8 * met(ib)%blyr_cond_h2o / (gsc * gsc_2_gsw8 + met(ib)%blyr_cond_h2o) 
        select case (thispft(ib)%photo_pathway)
        case (3)
            ! C3 plant
            select case (trim(limit_case))
            case ('RUBP')
            !-------------!
            ! RUBP limited
            !-------------!
            k1 = aparms(ib)%vm
            k2 = aparms(ib)%kco2 * (1.d0 + met(ib)%can_o2 / aparms(ib)%ko2)

            a = 1.d0
            b = (k2 - met(ib)%can_co2 + (k1 - aparms(ib)%leaf_resp) / gsbc)
            c = (-k1 * aparms(ib)%compp - k2 * aparms(ib)%leaf_resp) / gsbc - k2 * met(ib)%can_co2

            ! solve the quadratic equation
            rad = sqrt(b ** 2 - 4.d0 * a * c)
            ci = - (b - rad) / (2.d0 * a)
            fc = gsbc * (met(ib)%can_co2 - ci)

            ! calculate derivatives
            dbdg = -(k1 - aparms(ib)%leaf_resp) / gsc ** 2 ! Note that gbc cancelled out
            dcdg = (k1 * aparms(ib)%compp + k2 * aparms(ib)%leaf_resp) / gsc ** 2 !Note that gbc cancelled out

            if (abs(rad) < tiny_num8 ) then
                ! rad is effectively zero
                dcidg = -5.d-1 * dbdg
            else
                dcidg = -5.d-1 * dbdg + (5.d-1 * b * dbdg - dcdg) / rad
            endif
            dfcdg = (gsbc / gsc) ** 2 * (met(ib)%can_co2 - ci) &
                    + gsbc * (-1.d0 * dcidg)

            case ('LIGHT')
            !-------------!
            ! Light limited
            !-------------!
            k1 = aparms(ib)%jact / 4.d0
            k2 = 2.d0 * aparms(ib)%compp

            a = 1.d0
            b = (k2 - met(ib)%can_co2 + (k1 - aparms(ib)%leaf_resp) / gsbc)
            c = (-k1 * aparms(ib)%compp - k2 * aparms(ib)%leaf_resp) / gsbc - k2 * met(ib)%can_co2

            rad = sqrt(b ** 2 - 4.d0 * a * c)
            ci = - (b - rad) / (2.d0 * a)
            fc = gsbc * (met(ib)%can_co2 - ci)

            ! calculate derivatives
            dbdg = -(k1 - aparms(ib)%leaf_resp) / gsc ** 2 ! Note that gbc cancelled out
            dcdg = (k1 * aparms(ib)%compp + k2 * aparms(ib)%leaf_resp) / gsc ** 2 !Note that gbc cancelled out

            if (abs(rad) < tiny_num8 ) then
                ! rad is effectively zero
                dcidg = -5.d-1 * dbdg
            else
                dcidg = -5.d-1 * dbdg + (5.d-1 * b * dbdg - dcdg) / rad
            endif
            dfcdg = (gsbc / gsc) ** 2 * (met(ib)%can_co2 - ci) &
                    + gsbc * (-1.d0 * dcidg)


            case ('TPU')
            !-------------!
            ! TPU limited
            !-------------!
            ! only account for TPU limiation when iphysio is 1 or 3
            ! otherwise set fc to be a huge value so that TPU is always unlimited
                select case (iphysiol)
                case (0,2)
                    fc = huge(1.)
                    ci = tiny(1.)
                    dcidg = tiny(1.)
                    dfcdg = 0.
                case (1,3)
                    ! This is a linear case
                    fc = 3.d0 * aparms(ib)%tpm - aparms(ib)%leaf_resp
                    ci = met(ib)%can_co2 - fc / gsbc 

                    ! calculate derivatives
                    dcidg = fc / gsc ** 2 ! note gbc cancelledout
                    dfcdg = 0. ! constant fc
                end select

            end select
        case (4)
            ! C4 plant
            !-------------!
            select case (trim(limit_case))
            case ('RUBP')
            ! RUBP limited, linear case
            !-------------!
            fc = aparms(ib)%vm - aparms(ib)%leaf_resp
            ci = met(ib)%can_co2 - fc / gsbc 

            ! calculate derivatives
            dcidg = fc / gsc ** 2 ! note gbc cancelledout
            dfcdg = 0. ! constant fc
            
            case ('LIGHT')
            !-------------!
            ! Light limited, linear case
            !-------------!
            fc = aparms(ib)%jact / 4.d0 - aparms(ib)%leaf_resp
            ci = met(ib)%can_co2 - fc / gsbc

            ! calculate derivatives
            dcidg = fc / gsc ** 2 ! note gbc cancelledout
            dfcdg = 0. ! constant fc

            case ('CO2')
            !-------------!
            ! CO2 limited, linear case
            !-------------!
            ci = (aparms(ib)%leaf_resp + gsbc * met(ib)%can_co2) / (klowco28 * aparms(ib)%vm + gsbc)
            fc = gsbc * (met(ib)%can_co2 - ci)

            ! calculate derivatives
            dcidg = (klowco28 * aparms(ib)%vm - aparms(ib)%leaf_resp) &
                    / (klowco28 * aparms(ib)%vm + gsbc) ** 2            &
                    * (gsbc / gsc) ** 2
            dfcdg = (gsbc / gsc) ** 2 * (met(ib)%can_co2 - ci) &
                    + gsbc * (-1.d0 * dcidg)
            end select

        end select
        !------------------------------------------------------------------------------------!

        !------------------------------------------------------------------------------------!
        ! Second, Calculate fe and dfedg
        !------------------------------------------------------------------------------------!
        fe = gsbw * (met(ib)%lint_shv - met(ib)%can_shv)
        dfedg = (met(ib)%lint_shv - met(ib)%can_shv) &
              * (gsbw / (gsc * gsc_2_gsw8))**2 * gsc_2_gsw8  
        !------------------------------------------------------------------------------------!




        return
    end subroutine photosynthesis_stomata_solver8 
   !=======================================================================================!
   !=======================================================================================!


  end module farq_katul
  !=======================================================================================!
  !=======================================================================================!

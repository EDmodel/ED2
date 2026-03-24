!==========================================================================================!
!==========================================================================================!
!     This module contains the Photosynthesis model.  The references are:                  !
!                                                                                          !
! - M09 - Medvigy, D.M., S. C. Wofsy, J. W. Munger, D. Y. Hollinger, P. R. Moorcroft,      !
!         2009: Mechanistic scaling of ecosystem function and dynamics in space and time:  !
!         Ecosystem Demography model version 2.  J. Geophys. Res., 114, G01002,            !
!         doi:10.1029/2008JG000812.                                                        !
! - M06 - Medvigy, D.M., 2006: The State of the Regional Carbon Cycle: results from a      !
!         constrained coupled ecosystem-atmosphere model, 2006.  Ph.D. dissertation,       !
!         Harvard University, Cambridge, MA, 322pp.                                        !
! - M01 - Moorcroft, P. R., G. C. Hurtt, S. W. Pacala, 2001: A method for scaling          !
!         vegetation dynamics: the ecosystem demography model, Ecological Monographs, 71,  !
!         557-586.                                                                         !
! - F96 - Foley, J. A., I. Colin Prentice, N. Ramankutty, S. Levis, D. Pollard, S. Sitch,  !
!         A. Haxeltime, 1996: An integrated biosphere model of land surface processes,     !
!         terrestrial carbon balance, and vegetation dynamics. Glob. Biogeochem. Cycles,   !
!         10, 603-602.                                                                     !
! - L95 - Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf         !
!         nitrogen, photosynthesis, conductance, and transpiration: scaling from leaves to !
!         canopies. Plant, Cell, and Environ., 18, 1183-1200.                              !
! - F80 - Farquhar, G. D., S. von Caemmerer, J. A. Berry, 1980: A biochemical model of     !
!         photosynthetic  CO2 assimilation in leaves of C3 species. Planta, 149, 78-90.    !
! - C91 - Collatz, G. J., J. T. Ball, C. Grivet, J. A. Berry, 1991: Physiology and         !
!         environmental regulation of stomatal conductance, photosynthesis and transpir-   !
!         ation: a model that includes a laminar boundary layer, Agric. and Forest         !
!         Meteorol., 54, 107-136.                                                          !
! - C92 - Collatz, G. J., M. Ribas-Carbo, J. A. Berry, 1992: Coupled photosynthesis-       !
!         stomatal conductance model for leaves of C4 plants.  Aust. J. Plant Physiol.,    !
!         19, 519-538.                                                                     !
! - E78 - Ehleringer, J. R., 1978: Implications of quantum yield differences on the        !
!         distributions of C3 and C4 grasses.  Oecologia, 31, 255-267.                     !
!                                                                                          !
!------------------------------------------------------------------------------------------!
module farq_leuning
   use therm_lib8, only : toler8

   !---------------------------------------------------------------------------------------!
   !     This is a flag used in various sub-routines and functions and denote that we      !
   ! should ignore the result.                                                             !
   !---------------------------------------------------------------------------------------!
   real(kind=8), parameter :: discard = huge(1.d0)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    This is the tolerance for iterative methods.  Some of the functions are quite flat !
   ! so it is a good idea to use a somewhat more strict tolerance than the ones used in    !
   ! therm_lib8.                                                                           !
   !---------------------------------------------------------------------------------------!
   ! real(kind=8), parameter :: tolerfl8 = 1.d-10
   real(kind=8), parameter :: tolerfl8 = toler8
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Since this is more strict than therm_lib, the maximum number of attempts for the   !
   ! Regula Falsi method should be also increased.                                         !
   !---------------------------------------------------------------------------------------!
   integer     , parameter :: maxfpofl = 320
   !---------------------------------------------------------------------------------------!



   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !      This is the main driver for the photosynthesis model.                            !
   !---------------------------------------------------------------------------------------!
   subroutine lphysiol_full(ib,can_prss,can_rhos,can_shv,can_co2,ipft,leaf_par,leaf_temp   &
                           ,lint_shv,green_leaf_factor,leaf_aging_factor,vm_bar,rd_bar     &
                           ,leaf_gbw,leaf_D0,A_open,A_closed,A_light,A_rubp,A_co2,gsw_open &
                           ,gsw_closed,lsfc_shv_open,lsfc_shv_closed,lsfc_co2_open         &
                           ,lsfc_co2_closed,lint_co2_open,lint_co2_closed,leaf_resp,vmout  &
                           ,jmout,tpmout,jactout,comppout,limit_flag)
      use rk4_coms       , only : tiny_offset              & ! intent(in)
                                , effarea_transp           ! ! intent(in)
      use c34constants   , only : thispft                  & ! intent(out)
                                , met                      & ! intent(out)
                                , aparms                   & ! intent(out)
                                , stclosed                 & ! intent(inout)
                                , stopen                   & ! intent(inout)
                                , lightlim                 & ! intent(inout)
                                , rubpsatlim               & ! intent(inout)
                                , thirdlim                 ! ! intent(inout)
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
                                , stomatal_slope           & ! intent(in)
                                , quantum_efficiency       & ! intent(in)
                                , curvpar_electron         & ! intent(in)
                                , qyield_psII              ! ! intent(in)
      use physiology_coms, only : gbw_2_gbc8               & ! intent(in)
                                , o2_ref8                  ! ! intent(in)
      use consts_coms    , only : mmh2oi8                  & ! intent(in)
                                , mmh2o8                   & ! intent(in)
                                , mmdryi8                  & ! intent(in)
                                , mmdry8                   & ! intent(in)
                                , ep8                      & ! intent(in)
                                , epi8                     & ! intent(in)
                                , t008                     & ! intent(in)
                                , umol_2_mol8              & ! intent(in)
                                , mol_2_umol8              & ! intent(in)
                                , Watts_2_Ein8             ! ! intent(in)
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
      real(kind=4), intent(in)    :: leaf_D0           ! VPD factor for closure [  mol/mol]
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
      thispft(ib)%D0               = dble(leaf_D0)
      thispft(ib)%b                = dble(cuticular_cond(ipft)) * umol_2_mol8
      thispft(ib)%m                = dble(stomatal_slope(ipft))
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
      !     Call the main solver.                                                          !
      !------------------------------------------------------------------------------------!
      call photosynthesis_exact_solver(ib,limit_flag)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Copy the solution to the standard output variables.  Here we convert the values !
      ! back to the standard ED units                                                      !
      !------------------------------------------------------------------------------------!
      !----- Carbon demand, convert them to [umol/m2/s]. ----------------------------------!
      A_closed       = sngloff(stclosed(ib)%co2_demand    * mol_2_umol8 , tiny_offset)
      A_open         = sngloff(stopen(ib)%co2_demand      * mol_2_umol8 , tiny_offset)
      A_light        = sngloff(lightlim(ib)%co2_demand    * mol_2_umol8 , tiny_offset)
      A_rubp         = sngloff(rubpsatlim(ib)%co2_demand  * mol_2_umol8 , tiny_offset)
      A_co2          = sngloff(thirdlim(ib)%co2_demand    * mol_2_umol8 , tiny_offset)
      !----- Stomatal resistance, convert the conductances to [kg/m2/s]. ------------------!
      gsw_closed     = sngloff(stclosed(ib)%stom_cond_h2o * mmdry8 / effarea_transp(ipft)  &
                              , tiny_offset)
      gsw_open       = sngloff(stopen(ib)%stom_cond_h2o   * mmdry8 / effarea_transp(ipft)  &
                              , tiny_offset)
      !----- Leaf surface specific humidity, convert them to [kg/kg]. ---------------------!
      lsfc_shv_closed = sngloff(stclosed(ib)%lsfc_shv     * ep8         , tiny_offset)
      lsfc_shv_open   = sngloff(stopen(ib)%lsfc_shv       * ep8         , tiny_offset)
      !----- Leaf surface CO2 concentration, convert them to [umol/mol]. ------------------!
      lsfc_co2_closed = sngloff(stclosed(ib)%lsfc_co2     * mol_2_umol8 , tiny_offset)
      lsfc_co2_open   = sngloff(stopen(ib)%lsfc_co2       * mol_2_umol8 , tiny_offset)
      !----- Intercellular carbon dioxide concentration, convert them to [umol/mol]. ------!
      lint_co2_closed = sngloff(stclosed(ib)%lint_co2     * mol_2_umol8 , tiny_offset)
      lint_co2_open   = sngloff(stopen(ib)%lint_co2       * mol_2_umol8 , tiny_offset)
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
   end subroutine lphysiol_full
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes the maximum capacity of Rubisco to perform the           !
   ! carboxylase function at a certain temperature, given the PFT and phenology            !
   ! properties, and, in case of light controlled phenology, the leaf life span and the    !
   ! average reference value of the Vm function (Vm0).  The output variables are stored in !
   ! the aparms structure, as they are going to be used in other sub-routines. Both Vm and !
   ! the reference value Vm0 have units of umol/m2/s                                       !
   !     Compute the photosynthesis and leaf respiration parameters that depend on temper- !
   ! ature, according to the parameters defined in the "thispft" structure and the         !
   ! functional form chosen by the user.  The variables that are defined there are:        !
   ! - alpha     - the quantum yield, which may be a function of temperature.              !
   ! - Vm        - the maximum capacity of Rubisco to perform the carboxylase function.    !
   ! - Jm        - the maximum electron transport rate.                                    !
   ! - J         - the actual electron transport rate.                                     !
   ! - TPm       - the maximum triose phosphate utilisation rate.                          !
   ! - leaf_resp - the leaf respiration.                                                   !
   ! - compp     - the CO2 compensation point for gross photosynthesis (Gamma*)            !
   ! - kco2      - Michaelis-Mentel coefficient for CO2                                    !
   ! - ko2       - Michaelis-Mentel coefficient for O2.                                    !
   !---------------------------------------------------------------------------------------!
   subroutine comp_photo_tempfun(ib,leaf_aging_factor,green_leaf_factor)
      use physiology_coms, only : iphysiol              & ! intent(in)
                                , quantum_efficiency_T  & ! intent(in)
                                , qyield08              & ! intent(in)
                                , qyield18              & ! intent(in)
                                , qyield28              & ! intent(in)
                                , ehleringer_alpha0c8   & ! intent(in)
                                , compp_refval8         & ! intent(out)
                                , compp_hor8            & ! intent(out)
                                , compp_q108            & ! intent(out)
                                , kco2_refval8          & ! intent(out)
                                , kco2_hor8             & ! intent(out)
                                , kco2_q108             & ! intent(out)
                                , ko2_refval8           & ! intent(out)
                                , ko2_hor8              & ! intent(out)
                                , ko2_q108              ! ! intent(out)
      use c34constants   , only : met                   & ! intent(in) 
                                , aparms                & ! intent(in)
                                , thispft               ! ! intent(in)
      use consts_coms    , only : lnexp_min8            & ! intent(in)
                                , lnexp_max8            & ! intent(in)
                                , rmol8                 & ! intent(in)
                                , t008                  & ! intent(in)
                                , mol_2_umol8           ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      integer     , intent(in) :: ib                ! Multithread ID 
      real(kind=4), intent(in) :: leaf_aging_factor ! Ageing factor             [      ---]
      real(kind=4), intent(in) :: green_leaf_factor ! Greeness (prescr. phen.)  [      ---]
      !------ Local variables. ------------------------------------------------------------!
      real(kind=8)             :: vm_nocorr !> Vm  with no correction           [ mol/m2/s]
      real(kind=8)             :: jm_nocorr !> Jm  with no correction           [ mol/m2/s]
      real(kind=8)             :: rd_nocorr !> Rd  with no correction           [ mol/m2/s]
      real(kind=8)             :: lnexplow  !> Low temperature exponent         [      ---]
      real(kind=8)             :: lnexphigh !> High temperature exponent        [      ---]
      real(kind=8)             :: tlow_fun  !> Low temperature bound            [      ---]
      real(kind=8)             :: thigh_fun !> High temperature bound           [      ---]
      real(kind=8)             :: greeness  !> Leaf "Greeness"                  [   0 to 1]
      real(kind=8)             :: ipsII     !> Light used by photosystem II     [ mol/m2/s]
      real(kind=8)             :: ipsII_max !> critical value for ipsII         [      ---]
      real(kind=8)             :: aterm     !> 'a' term for quadratic equation  [      ---]
      real(kind=8)             :: bterm     !> 'b' term for quadratic equation  [      ---]
      real(kind=8)             :: cterm     !> 'c' term for quadratic equation  [      ---]
      real(kind=8)             :: jroot1    !> root for quadratic equation      [      ---]
      real(kind=8)             :: jroot2    !> root for quadratic equation      [      ---]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     If this plant functional type has cold phenology associated with it, we must   !
      ! determine the correction term for Vm.  Otherwise, we set the greeness to 1.        !
      !------------------------------------------------------------------------------------!
      if (leaf_aging_factor > 0.01 .and. green_leaf_factor > 0.0001) then
         greeness = dble(leaf_aging_factor) / dble(green_leaf_factor)
      else
         greeness = 1.d0
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Decide whether to make alpha dependent on temperature or not.  This is only     !
      ! done for C3 plants, C4 plants always have constant alpha in our model.             !
      !------------------------------------------------------------------------------------!
      select case (thispft(ib)%photo_pathway)
      case (3)
         select case(quantum_efficiency_T)
         case (0)
            !----- The user wants constant alpha for C3 plants. ---------------------------!
            aparms(ib)%alpha = thispft(ib)%alpha0
            !------------------------------------------------------------------------------!

         case (1)
            !------------------------------------------------------------------------------!
            !    The user wants alpha as a function of temperature.  The equation below is !
            ! equation 2 from E78, modified to take in temperature in Kelvin.              !
            !------------------------------------------------------------------------------!
            if (met(ib)%leaf_temp > t008) then
               aparms(ib)%alpha = qyield08 + met(ib)%leaf_temp                             &
                                           * (qyield18   + qyield28 * met(ib)%leaf_temp)
            else
               !---------------------------------------------------------------------------!
               !     We don't apply the equation when the temperature is below freezing.   !
               ! If this is the case, we use the value at 0 degC instead.                  !
               !---------------------------------------------------------------------------!
               aparms(ib)%alpha = ehleringer_alpha0c8
            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!

      case (4)
         aparms(ib)%alpha = thispft(ib)%alpha0

      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Here we decide which temperature-dependence to use, the "Arrhenius" form as in !
      ! F96, or the power law as in C91.  Solve all the temperature functions first, then  !
      ! we apply the temperature correction term.                                          !
      !------------------------------------------------------------------------------------!
      select case (iphysiol)
      case (0,1)
         !---------------------------------------------------------------------------------!
         !     Use  the "Arrhenius" form, like in F96 and M01.                             !
         !---------------------------------------------------------------------------------!


         !----- Find Vm, with no correction. ----------------------------------------------!
         vm_nocorr = greeness * arrhenius(met(ib)%leaf_temp,thispft(ib)%vm0                &
                                         ,thispft(ib)%vm_hor)
         !---------------------------------------------------------------------------------!


         !----- Find Jm, with no correction. ----------------------------------------------!
         jm_nocorr = greeness * arrhenius(met(ib)%leaf_temp,thispft(ib)%jm0                &
                                         ,thispft(ib)%jm_hor)
         !---------------------------------------------------------------------------------!


         !----- Find Rd, with no correction. ----------------------------------------------!
         rd_nocorr = greeness * arrhenius(met(ib)%leaf_temp,thispft(ib)%rd0                &
                                         ,thispft(ib)%rd_hor)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    If this is a C3 plant, find the compensation point, and the Michaelis-Mentel !
         ! constants for CO2 and O2.  Otherwise, assign zeroes for compensation point and  !
         ! the Michaelis-Mentel constant for CO2.  The oxygen one should have no impact,   !
         ! but we always assign it to avoid divisions by zero.                             !
         !---------------------------------------------------------------------------------!
         select case (thispft(ib)%photo_pathway)
         case (3)
            aparms(ib)%compp = arrhenius(met(ib)%leaf_temp,compp_refval8,compp_hor8)
            aparms(ib)%kco2  = arrhenius(met(ib)%leaf_temp,kco2_refval8 ,kco2_hor8 )
         case (4)
            aparms(ib)%compp = 0.d0
            aparms(ib)%kco2  = 0.d0
         end select
         aparms(ib)%ko2   = arrhenius(met(ib)%leaf_temp,ko2_refval8,ko2_hor8)
         !---------------------------------------------------------------------------------!




      case (2,3)
         !---------------------------------------------------------------------------------!
         !     Use  the "Collatz" (Q10) form, like in C91 and C92.                         !
         !---------------------------------------------------------------------------------!


         !----- Find Vm, with no correction. ----------------------------------------------!
         vm_nocorr = greeness * collatz(met(ib)%leaf_temp,thispft(ib)%vm0                  &
                                       ,thispft(ib)%vm_q10)
         !---------------------------------------------------------------------------------!


         !----- Find Jm, with no correction. ----------------------------------------------!
         jm_nocorr = greeness * collatz(met(ib)%leaf_temp,thispft(ib)%jm0                  &
                                       ,thispft(ib)%jm_q10)
         !---------------------------------------------------------------------------------!


         !----- Find Rd, with no correction. ----------------------------------------------!
         rd_nocorr = greeness * collatz(met(ib)%leaf_temp,thispft(ib)%rd0                  &
                                       ,thispft(ib)%rd_q10)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    If this is a C3 plant, find the compensation point, and the Michaelis-Mentel !
         ! constants for CO2 and O2.  Otherwise, assign zeroes for compensation point and  !
         ! the Michaelis-Mentel constant for CO2.  The oxygen one should have no impact,   !
         ! but we always assign it to avoid divisions by zero.                             !
         !---------------------------------------------------------------------------------!
         select case (thispft(ib)%photo_pathway)
         case (3)
            aparms(ib)%compp = collatz(met(ib)%leaf_temp,compp_refval8,compp_q108)
            aparms(ib)%kco2  = collatz(met(ib)%leaf_temp,kco2_refval8 ,kco2_q108 )
         case (4)
            aparms(ib)%compp = 0.d0
            aparms(ib)%kco2  = 0.d0
         end select
         aparms(ib)%ko2   = collatz(met(ib)%leaf_temp,ko2_refval8,ko2_q108)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Apply further temperature correction to shut down Vm at extreme temperatures.   !
      ! In order to avoid floating point exceptions, we check whether the temperature will !
      ! make the exponential too large or too small.                                       !
      !------------------------------------------------------------------------------------!
      !----- Low temperature. -------------------------------------------------------------!
      lnexplow  = thispft(ib)%vm_decay_elow  * (thispft(ib)%vm_low_temp -met(ib)%leaf_temp)
      lnexplow  = max(lnexp_min8,min(lnexp_max8,lnexplow ))
      tlow_fun  = 1.d0 +  exp(lnexplow )
      !----- High temperature. ------------------------------------------------------------!
      lnexphigh = thispft(ib)%vm_decay_ehigh * (met(ib)%leaf_temp-thispft(ib)%vm_high_temp)
      lnexphigh = max(lnexp_min8,min(lnexp_max8,lnexphigh))
      thigh_fun = 1.d0 + exp(lnexphigh)
      !------ Correct Vm. -----------------------------------------------------------------!
      aparms(ib)%vm = vm_nocorr / (tlow_fun * thigh_fun)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Apply further temperature correction to shut down Jm at extreme temperatures.   !
      ! In order to avoid floating point exceptions, we check whether the temperature will !
      ! make the exponential too large or too small.                                       !
      !------------------------------------------------------------------------------------!
      !----- Low temperature. -------------------------------------------------------------!
      lnexplow  = thispft(ib)%jm_decay_elow  * (thispft(ib)%jm_low_temp -met(ib)%leaf_temp)
      lnexplow  = max(lnexp_min8,min(lnexp_max8,lnexplow ))
      tlow_fun  = 1.d0 +  exp(lnexplow )
      !----- High temperature. ------------------------------------------------------------!
      lnexphigh = thispft(ib)%jm_decay_ehigh * (met(ib)%leaf_temp-thispft(ib)%jm_high_temp)
      lnexphigh = max(lnexp_min8,min(lnexp_max8,lnexphigh))
      thigh_fun = 1.d0 + exp(lnexphigh)
      !------ Correct jm. -----------------------------------------------------------------!
      aparms(ib)%jm = jm_nocorr / (tlow_fun * thigh_fun)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Apply further temperature correction to shut down Rd at extreme temperatures.   !
      ! In order to avoid floating point exceptions, we check whether the temperature will !
      ! make the exponential too large or too small.                                       !
      !------------------------------------------------------------------------------------!
      !----- Low temperature. -------------------------------------------------------------!
      lnexplow  = thispft(ib)%rd_decay_elow  * (thispft(ib)%rd_low_temp -met(ib)%leaf_temp)
      lnexplow  = max(lnexp_min8,min(lnexp_max8,lnexplow ))
      tlow_fun  = 1.d0 +  exp(lnexplow )
      !----- High temperature. ------------------------------------------------------------!
      lnexphigh = thispft(ib)%rd_decay_ehigh * (met(ib)%leaf_temp-thispft(ib)%rd_high_temp)
      lnexphigh = max(lnexp_min8,min(lnexp_max8,lnexphigh))
      thigh_fun = 1.d0 + exp(lnexphigh)
      !------ Correct leaf respiration. ---------------------------------------------------!
      aparms(ib)%leaf_resp = rd_nocorr / (tlow_fun * thigh_fun)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the maximum triose phosphate utilisation rate based on Vm and reference   !
      ! values.                                                                            !
      !------------------------------------------------------------------------------------!
      aparms(ib)%tpm = aparms(ib)%vm * thispft(ib)%TPm0 / thispft(ib)%vm0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the actual electron transport electron rate.  This is solved only when    !
      ! iphysiol is 1 or 3.  For the old ED-method (cases 0 and 2), we assume that         !
      ! J = 4 * alpha * PAR, which makes the light-limited case the same as the original.  !
      !------------------------------------------------------------------------------------!
      select case (iphysiol)
      case (0,2)
         !----- Quantum efficiency method. ------------------------------------------------!
         aparms(ib)%jact = 4.d0 * aparms(ib)%alpha * met(ib)%par
         !---------------------------------------------------------------------------------!
      case (1,3)
         !---------------------------------------------------------------------------------!
         !      Check whether there is a risk of negative discriminant.  This occurs only  !
         ! when the curvature parameter is greater than 1.0.                               !
         !---------------------------------------------------------------------------------!
         if (thispft(ib)%curvpar > 1.d0) then
            ipsII_max = aparms(ib)%jm                                                      &
                      * ( 2.d0 * thispft(ib)%curvpar - 1.d0                                &
                        - sqrt(thispft(ib)%curvpar * (thispft(ib)%curvpar-1.d0) ) )
         else
            ipsII_max = discard
         end if
         !---------------------------------------------------------------------------------!


         !----- Find absorbed light used by photosystem II. -------------------------------!
         ipsII = min(ipsII_max,0.5 * thispft(ib)%phi_psII * met(ib)%par)
         !---------------------------------------------------------------------------------!


         !----- Find the quadratic terms. -------------------------------------------------!
         aterm = thispft(ib)%curvpar
         bterm = - (ipsII + aparms(ib)%jm)
         cterm = ipsII * aparms(ib)%jm
         !---------------------------------------------------------------------------------!



         !----- Solve the quadratic equation. ---------------------------------------------!
         call solve_quadratic8(aterm,bterm,cterm,discard,jroot1,jroot2)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Make sure at least one root is valid.                                       !
         !---------------------------------------------------------------------------------!
         if ( ( jroot1 /= discard ) .or. ( jroot2 /= discard ) ) then
            aparms(ib)%jact = min( jroot1, jroot2 )
         else
            !----- Broadcast the bad news. ------------------------------------------------!
            write(unit=*,fmt='(a)'          ) ''
            write(unit=*,fmt='(a)'          ) ''
            write(unit=*,fmt='(a)'          ) '==========================================='
            write(unit=*,fmt='(a)'          ) '==========================================='
            write(unit=*,fmt='(a)'          ) '   Discriminant is negative.'
            write(unit=*,fmt='(a)'          ) '-------------------------------------------'
            write(unit=*,fmt='(a,1x,f12.5)' ) 'PAR       = ', met(ib)%par * mol_2_umol8
            write(unit=*,fmt='(a,1x,f12.5)' ) 'QYIELD    = ', thispft(ib)%phi_psII
            write(unit=*,fmt='(a,1x,f12.5)' ) 'CURVATURE = ', thispft(ib)%curvpar
            write(unit=*,fmt='(a,1x,f12.5)' ) 'JMAX      = ', aparms(ib)%jm   * mol_2_umol8
            write(unit=*,fmt='(a,1x,f12.5)' ) 'JM0       = ', thispft(ib)%jm0 * mol_2_umol8
            write(unit=*,fmt='(a,1x,f12.5)' ) 'LEAF_TEMP = ', met(ib)%leaf_temp + t008
            write(unit=*,fmt='(a,1x,es12.5)') 'A         = ', aterm
            write(unit=*,fmt='(a,1x,es12.5)') 'B         = ', bterm
            write(unit=*,fmt='(a,1x,es12.5)') 'C         = ', cterm
            write(unit=*,fmt='(a,1x,es12.5)') 'DISCR     = ', bterm*bterm-4.d0*aterm*cterm
            write(unit=*,fmt='(a)'          ) '==========================================='
            write(unit=*,fmt='(a)'          ) '==========================================='
            write(unit=*,fmt='(a)'          ) ''
            write(unit=*,fmt='(a)'          ) ''
            call fatal_error(' Negative discriminant when seeking the J rate.'             &
                            ,'comp_photo_tempfun','farq_leuning.f90')
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine comp_photo_tempfun
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine is the main driver for the C3 photosynthesis.                     !
   !---------------------------------------------------------------------------------------!
   subroutine photosynthesis_exact_solver(ib,limit_flag)
      use c34constants   , only : met              & ! intent(in)
                                , thispft          & ! intent(in)
                                , stopen           & ! intent(inout)
                                , stclosed         & ! intent(inout)
                                , rubpsatlim       & ! intent(inout)
                                , thirdlim         & ! intent(inout)
                                , lightlim         & ! intent(inout)
                                , copy_solution    ! ! intent(in)
      use physiology_coms, only : iphysiol         ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      integer     , intent(in)  :: ib               ! Multithread ID
      integer     , intent(out) :: limit_flag       ! Flag with limiting case    [     ---]
      !------ Local variables. ------------------------------------------------------------!
      real(kind=8)              :: par_twilight_min ! Minimum daytime radiation  [mol/m2/s]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Initialise the limitation flag with -1.  This will mean that none of the      !
      ! solutions worked and it is daytime.  We will replace the flag depending on the     !
      ! success of solutions.                                                              !
      !------------------------------------------------------------------------------------!
      limit_flag = -1
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Initialise the parameters to compute the carbon demand for the case where the !
      ! stomata are closed.  We call the solver for the specific case in which the stomata !
      ! are closed, because in this case neither the carbon demand nor the stomatal water  !
      ! conductance depend on the intercellular CO2 mixing ratio.                          !
      !------------------------------------------------------------------------------------!
      call set_co2_demand_params(ib,'CLOSED')
      call solve_aofixed_case(ib,stclosed(ib))
      if (.not. stclosed(ib)%success) then
         call fatal_error ('Solution failed for closed case'                               &
                          ,'photosynthesis_exact_solver','farq_leuning.f90')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    There is enough light to be considered at least dawn or dusk, so we go with     !
      ! photosynthesis.  The closed stomata is already solved, so now we must solve the    !
      ! open stomata case only.                                                            !
      !    During the day, conditions may lead to light, rubisco, or CO2 limitation on the !
      ! photosynthesis.  Following M01, we solve all cases and determine which ones gives  !
      ! the lowest carbon demand.  The actual solution is going to be simply the one with  !
      ! the lowest carbon_demand.                                                          !
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    First we check whether it is at least dawn or dusk.  In case it is not, no      !
      ! photosynthesis should happen, so we copy the closed case stomata values to the     !
      ! open case.  Limit_flag becomes 0, which is the flag for night time limitation.     !
      !------------------------------------------------------------------------------------!
      par_twilight_min = find_twilight_min(ib)
      !------------------------------------------------------------------------------------!



      !----- Check the amount of light before running the light-limited case. -------------!
      if (met(ib)%par < par_twilight_min) then
         limit_flag = 0
         call copy_solution(stclosed(ib),lightlim(ib))
      else
         !----- Light above minimum, run the RubP-regeneration (aka light) limited case. --!
         call set_co2_demand_params(ib,'LIGHT')
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Day time: choose the appropriate solver depending on the photosynthetic     !
         ! pathway.                                                                        !
         !---------------------------------------------------------------------------------!
         select case(thispft(ib)%photo_pathway)
         case (3)
            call solve_iterative_case(ib,lightlim(ib))
         case (4)
            call solve_aofixed_case(ib,lightlim(ib))
         end select
         !---------------------------------------------------------------------------------!
         !     In case success was returned as "false", this means that the light-limited  !
         ! case didn't converge (there was no root).  If this is the case we don't let     !
         ! this solution to be selected (but do not give up because this case may fail     !
         ! when light is excessive).                                                       !
         !---------------------------------------------------------------------------------!
         if (lightlim(ib)%success) then
            limit_flag              = 99
         else
            lightlim(ib)%lint_co2   = -discard
            lightlim(ib)%co2_demand = + discard
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    2. The RubP-saturated (aka Rubisco-limited) case.                               !
      !------------------------------------------------------------------------------------!
      !----- Update the CO2 demand function parameters for Rubisco limitation. ------------!
      call set_co2_demand_params(ib,'RUBP_SAT')
      !----- Choose the appropriate solver depending on the kind of photosynthesis. -------!
      select case(thispft(ib)%photo_pathway)
      case (3)
         call solve_iterative_case(ib,rubpsatlim(ib))
      case (4)
         call solve_aofixed_case(ib,rubpsatlim(ib))
      end select
      !------------------------------------------------------------------------------------!
      !     In case success was returned as "false", this means that the RubP-             !
      ! saturated case didn't converge (there was no root).  If this is the case we don't  !
      ! let this solution to be selected.  Furthermore, if this is C3, it means that       !
      ! the solution is invalid, in which case we close the stomata.                       !
      !------------------------------------------------------------------------------------!
      if (rubpsatlim(ib)%success) then
         limit_flag    = 99
      else
         select case (thispft(ib)%photo_pathway)
         case (3)
            limit_flag = -2
         case (4)
            rubpsatlim(ib)%lint_co2   = -discard
            rubpsatlim(ib)%co2_demand = +discard
         end select
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !  3. The third case.  This is the triose phosphate utilisation limitation for C3    !
      !     plants, and the low CO2 concentration case for C4 plants.                      !
      !------------------------------------------------------------------------------------!
      !----- Choose the appropriate solver depending on the kind of photosynthesis. -------!
      select case(thispft(ib)%photo_pathway)
      case (3)
         !---------------------------------------------------------------------------------!
         !     The triose phosphate utilisation limitation was not originally solved in    !
         ! ED-1/ED-2.0, so we skip it for back-compatibility.                              !
         !---------------------------------------------------------------------------------!
         select case (iphysiol)
         case (0,2)
            !----- Old-style solution, ignore triose phosphate case. ----------------------!
            call copy_solution(rubpsatlim(ib),thirdlim(ib))
            !------------------------------------------------------------------------------!
         case (1,3)
            !----- Set parameters for TPU limitation and solve it. ------------------------!
            call set_co2_demand_params(ib,'TPU')
            call solve_aofixed_case(ib,thirdlim(ib))
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      case (4)
         !---------------------------------------------------------------------------------!
         !    C4, both carbon_demand and stomatal conductance of water depend on the       !
         ! inter-cellular CO2 concentration.  We must find all these three variables       !
         ! simultaneously, using an iterative method.                                      !
         !---------------------------------------------------------------------------------!
         call set_co2_demand_params(ib,'CO2')
         call solve_iterative_case(ib,thirdlim(ib))
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!
      !     In case success was returned as "false", this means that the CO2-limited case  !
      ! didn't converge (there was no root).  If this is the case we don't let             !
      ! this solution to be selected.  Furthermore, if this is C4, it means that the       !
      !  solution is invalid, so we close the stomata.                                     !
      !------------------------------------------------------------------------------------!
      if (thirdlim(ib)%success) then
         limit_flag = 99
      else
         select case (thispft(ib)%photo_pathway)
         case (3)
            thirdlim(ib)%lint_co2   = -discard
            thirdlim(ib)%co2_demand = +discard
         case (4)
            limit_flag = -3
         end select
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     The actual solution will be the one with the lowest carbon demand.  In case    !
      ! it is night time or we do not have good solutions, we close all stomata.           !
      !------------------------------------------------------------------------------------!
      select case (limit_flag)
      case (99)
         !---------------------------------------------------------------------------------!
         !      At least one plausible solution was found, pick the most limiting one.     !
         !---------------------------------------------------------------------------------!
         if (lightlim(ib)%co2_demand <= rubpsatlim(ib)%co2_demand .and.                    &
             lightlim(ib)%co2_demand <= thirdlim(ib)%co2_demand              )  then
            !----- RubP-regeneration (aka light) is limiting photosynthesis. --------------!
            call copy_solution(lightlim(ib),stopen(ib))
            limit_flag = 1
            !------------------------------------------------------------------------------!
         elseif (rubpsatlim(ib)%co2_demand <= thirdlim(ib)%co2_demand        ) then
            !----- RubP is saturated. -----------------------------------------------------!
            call copy_solution(rubpsatlim(ib),stopen(ib))
            limit_flag = 2
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !     Triose Phosphate Utilisation(C3) or CO2 (C4) is limiting photosynthesis. !
            !------------------------------------------------------------------------------!
            call copy_solution(thirdlim(ib),stopen(ib))
            limit_flag = 3
         end if
         !---------------------------------------------------------------------------------!
      case default
         !---------------------------------------------------------------------------------!
         !     No suitable solution found, because it is either night time or the iter-    !
         ! ative cases (RubP-sat for C3; CO2 for C4) failed.  Close stomata.               !
         !---------------------------------------------------------------------------------!
         call copy_solution(stclosed(ib),stopen  (ib))
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    To avoid FPE problems, we copy the closed case to unsuccessful solutions.       !
      !------------------------------------------------------------------------------------!
      if (.not. lightlim(ib)%success  ) call copy_solution(stclosed(ib),lightlim  (ib))
      if (.not. rubpsatlim(ib)%success) call copy_solution(stclosed(ib),rubpsatlim(ib))
      if (.not. thirdlim(ib)%success  ) call copy_solution(stclosed(ib),thirdlim  (ib))
      !------------------------------------------------------------------------------------!

      return
   end subroutine photosynthesis_exact_solver
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine solves the model for the case where the carbon demand doesn't     !
   ! depend on the internal carbon.  This is simpler than the iterative case because we    !
   ! can solve through a quadratic equation for the stomatal conductance for water vapour. !
   !---------------------------------------------------------------------------------------!
   subroutine solve_aofixed_case(ib,answer)
      use c34constants   , only : solution_vars     & ! structure
                                , met               & ! intent(in)
                                , thispft           & ! intent(in)
                                , aparms            ! ! intent(in)
      use physiology_coms, only : gsw_2_gsc8        & ! intent(in)
                                , c34smax_gsw8      & ! intent(in)
                                , c34smin_lint_co28 & ! intent(in)
                                , c34smax_lint_co28 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                          :: ib
      type(solution_vars), intent(out) :: answer
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                     :: qterm1
      real(kind=8)                     :: qterm2
      real(kind=8)                     :: qterm3
      real(kind=8)                     :: aquad
      real(kind=8)                     :: bquad
      real(kind=8)                     :: cquad
      real(kind=8)                     :: gswroot1
      real(kind=8)                     :: gswroot2
      real(kind=8)                     :: ciroot1
      real(kind=8)                     :: ciroot2
      logical                          :: bounded1
      logical                          :: bounded2
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !   1. Initialise the success flag as false, and switch it to true only when (and    !
      !      if) a valid and plausible solution is found.                                  !
      !------------------------------------------------------------------------------------!
      answer%success = .false.
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   2. Since the carbon demand doesn't depend on the intercellular CO2, compute it   !
      !      using the first guess.                                                        !
      !------------------------------------------------------------------------------------!
      answer%co2_demand = calc_co2_demand(ib,met(ib)%can_co2)
      !------------------------------------------------------------------------------------!



      !----- 3. Compute the leaf surface CO2. ---------------------------------------------!
      answer%lsfc_co2  = met(ib)%can_co2 - answer%co2_demand / met(ib)%blyr_cond_co2
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !       In case the leaf surface CO2 is unrealistic (i.e. negative), don't bother    !
      ! looking for a solution.                                                            !
      !------------------------------------------------------------------------------------!
      if (answer%lsfc_co2 < c34smin_lint_co28) return
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   3. Here we check the sign of the carbon demand.                                  !
      !------------------------------------------------------------------------------------!
      if (answer%co2_demand <= 0.d0) then
         !---------------------------------------------------------------------------------!
         !     If carbon demand is zero or negative, this means that light is below the    !
         ! light compensation point, so all stomata should remain closed.                  !
         !---------------------------------------------------------------------------------!
         answer%stom_cond_h2o = thispft(ib)%b
         answer%stom_cond_co2 = gsw_2_gsc8 * answer%stom_cond_h2o
         answer%lint_co2      = answer%lsfc_co2 - answer%co2_demand / answer%stom_cond_co2
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Carbon demand is positive, look for a solution.                             !
         !---------------------------------------------------------------------------------!
         !----- Find auxiliary coefficients to compute the quadratic terms. ---------------!
         qterm1 = (met(ib)%can_co2 - aparms(ib)%compp) * met(ib)%blyr_cond_co2             &
                - answer%co2_demand
         qterm2 = (thispft(ib)%d0 + met(ib)%lint_shv - met(ib)%can_shv)                    &
                * met(ib)%blyr_cond_h2o
         qterm3 = thispft(ib)%m * answer%co2_demand * thispft(ib)%d0                       &
                * met(ib)%blyr_cond_co2
         !----- Find the coefficients for the quadratic equation. -------------------------!
         aquad = qterm1 * thispft(ib)%d0
         bquad = qterm1 * qterm2 - aquad * thispft(ib)%b - qterm3
         cquad = - qterm1 * qterm2 * thispft(ib)%b - qterm3 * met(ib)%blyr_cond_h2o
         !----- Solve the quadratic equation for gsw. -------------------------------------!
         call solve_quadratic8(aquad,bquad,cquad,discard,gswroot1,gswroot2)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Check both solutions, and decide which one makes sense.  Once the right     !
         ! one is determined, compute the stomatal resistance for CO2 and the inter-       !
         ! cellular CO2 concentration.  In case both make solutions make sense (unlikely), !
         ! we decide the root based on the intercellular CO2.                              !
         !---------------------------------------------------------------------------------!
         bounded1 = gswroot1 >= thispft(ib)%b .and. gswroot1 <= c34smax_gsw8
         bounded2 = gswroot2 >= thispft(ib)%b .and. gswroot2 <= c34smax_gsw8
         if (bounded1 .and. bounded2) then
            !----- Both solutions are valid, warn the user as this should never happen. ---!
            ciroot1 = answer%lsfc_co2 - answer%co2_demand / (gsw_2_gsc8 * gswroot1)
            ciroot2 = answer%lsfc_co2 - answer%co2_demand / (gsw_2_gsc8 * gswroot2)

            bounded1 = ciroot1 >= c34smin_lint_co28 .and.                                  &
                       ciroot1 <= min(c34smax_lint_co28,met(ib)%can_co2)
            bounded2 = ciroot2 >= c34smin_lint_co28 .and.                                  &
                       ciroot2 <= min(c34smax_lint_co28,met(ib)%can_co2)
             
            if (bounded1 .and. bounded2) then
               !----- Both intercellular CO2 work, pick the highest and warn the user. ----!
               if (ciroot1 >= ciroot2) then
                  answer%stom_cond_h2o = gswroot1
                  answer%stom_cond_co2 = gsw_2_gsc8 * answer%stom_cond_h2o
                  answer%lint_co2      = ciroot1
               else
                  answer%stom_cond_h2o = gswroot2
                  answer%stom_cond_co2 = gsw_2_gsc8 * answer%stom_cond_h2o
                  answer%lint_co2      = ciroot2
               end if
            elseif (bounded1) then
               answer%stom_cond_h2o = gswroot1
               answer%stom_cond_co2 = gsw_2_gsc8 * answer%stom_cond_h2o
               answer%lint_co2      = ciroot1
            elseif (bounded2) then
               answer%stom_cond_h2o = gswroot2
               answer%stom_cond_co2 = gsw_2_gsc8 * answer%stom_cond_h2o
               answer%lint_co2      = ciroot2
            else
               return
            end if
         elseif (bounded1) then
            !----- First root is the only one that makes sense. ---------------------------!
            answer%stom_cond_h2o = gswroot1
            answer%stom_cond_co2 = gsw_2_gsc8 * answer%stom_cond_h2o
            answer%lint_co2      = answer%lsfc_co2                                         &
                                 - answer%co2_demand / answer%stom_cond_co2

         elseif (bounded2) then
            !----- Second root is the only one that makes sense. --------------------------!
            answer%stom_cond_h2o = gswroot2
            answer%stom_cond_co2 = gsw_2_gsc8 * answer%stom_cond_h2o
            answer%lint_co2      = answer%lsfc_co2                                         &
                                 - answer%co2_demand / answer%stom_cond_co2
         else
            !----- None of the solutions are bounded.  This solution failed. --------------!
            return
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   8. Lastly, find the surface water specific humidity.                             !
      !------------------------------------------------------------------------------------!
      answer%lsfc_shv = ( answer%stom_cond_h2o * met(ib)%lint_shv                          &
                        + met(ib)%blyr_cond_h2o    * met(ib)%can_shv  )                    &
                      / ( met(ib)%blyr_cond_h2o + answer%stom_cond_h2o)
      !------------------------------------------------------------------------------------!


      !----- Successful solution, update the logical flag before quitting. ----------------!
      answer%success = .true.
      !------------------------------------------------------------------------------------!

      return
   end subroutine solve_aofixed_case
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine will solve the case in which both the carbon demand and the      !
   ! stomatal conductance of water are functions of the intercellular CO2.  This function  !
   ! can be used for the simpler cases too, but it's not advisable because this is based   !
   ! on iterative methods, which makes the solution slower.                                !
   !     The iterative method is designed to use Newton's method as the default, and this  !
   ! should take care of most cases.  In case Newton's method fails, it will fall back to  !
   ! the modified Regula Falsi method (Illinois) and look for guesses with opposite sign.  !
   ! If the method fails finding the pair, it means that there is no viable solution       !
   ! within this range, so the method quits and return the error message.                  !
   !---------------------------------------------------------------------------------------!
   subroutine solve_iterative_case(ib,answer)
      use c34constants   , only : solution_vars     & ! structure
                                , met               ! ! intent(in)
      use physiology_coms, only : gsw_2_gsc8        ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                          :: ib        ! Multithread ID
      type(solution_vars), intent(out) :: answer    ! The strutcure with the answer
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                     :: ci        ! Intercellular CO2         [  mol/mol]
      real(kind=8)                     :: cia       ! Smallest/previous guess   [  mol/mol]
      real(kind=8)                     :: ciz       ! Largest/new guess         [  mol/mol]
      real(kind=8)                     :: deriv     ! Function derivative       [  mol/mol]
      real(kind=8)                     :: fun       ! Function evaluation       [mol2/mol2]
      real(kind=8)                     :: funa      ! Smallest  guess function  [mol2/mol2]
      real(kind=8)                     :: funz      ! Largest   guess function  [mol2/mol2]
      real(kind=8)                     :: delta     ! Aux. var to find 2nd guess[         ]
      real(kind=8)                     :: cimin     ! Minimum intercell. CO2    [  mol/mol]
      real(kind=8)                     :: cimax     ! Maximum intercell. CO2    [  mol/mol]
      integer                          :: itn       ! Iteration counter         [      ---]
      integer                          :: itb       ! Iteration counter         [      ---]
      logical                          :: converged ! A solution was found     [      T|F]
      logical                          :: zside     ! Check for 1-sided appr.   [      ---]
      logical                          :: hitmin    ! 2nd guess tried minimum   [      ---]
      logical                          :: hitmax    ! 2nd guess tried maximum   [      ---]
      logical                          :: bounded   ! Guess range is bounded.   [      T|F]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Initialise the convergence and success flags.  Set both to false, and change  !
      ! them only when (and if) a real, plausible answer is found.  In case such           !
      ! solution doesn't exist, we report back the failure so this case is not selected.   !
      !------------------------------------------------------------------------------------!
      converged      = .false.
      answer%success = .false.
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Determine the minimum and maximum intercellular CO2 that will still produce  a  !
      ! positive and bounded conductance, which should be above the cuticular conductance  !
      ! (otherwise there is no reason to keep stomata opened.                              !
      !------------------------------------------------------------------------------------!
      call find_lint_co2_bounds(ib,cimin,cimax,bounded)
      if (bounded) then
         !---------------------------------------------------------------------------------!
         !     We have found bounds, find put the first guess in the middle and find the   !
         ! function evaluation and derivative of this guess.                               !
         !---------------------------------------------------------------------------------!
         cia = sqrt(cimin*cimax)
         call iter_solver_step(ib,.true.,cia,funa,deriv)
      else
         !----- No reasonable bound was found, we don't even bother solving this case. ----!
         return
      end if
      !------------------------------------------------------------------------------------!



      !----- Copy to the new guess just in case it fails at the first iteration -----------!
      ciz = cia
      fun = funa
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Also check whether we have been lucky and have hit the jackpot (answer is     !
      ! ciz).                                                                              !
      !------------------------------------------------------------------------------------!
      if (fun == 0.d0) then 
         ci        = ciz
         converged = .true.
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Enter Newton's method loop in case we haven't found the answer.                !
      !------------------------------------------------------------------------------------!
      if (.not. converged) then
         newloop: do itn = 1,nint(maxfpofl/6.)
            !------------------------------------------------------------------------------!
            !    In case the derivative is bad, we give up on Newton's and go with Regula  !
            ! Falsi.                                                                       !
            !------------------------------------------------------------------------------!
            if (abs(deriv) < tolerfl8) exit newloop
            !------------------------------------------------------------------------------!


            !----- Copy the previous guess. -----------------------------------------------!
            cia   = ciz
            funa  = fun
            !------------------------------------------------------------------------------!


            !----- Update guess. ----------------------------------------------------------!
            ciz      = cia - fun/deriv
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Check whether the method converged.                                      !
            !------------------------------------------------------------------------------!
            converged = 2.d0 * abs(cia-ciz) < tolerfl8 * (abs(cia)+abs(ciz))
            !------------------------------------------------------------------------------!

 
            !------------------------------------------------------------------------------!
            !    At this point we test the current status of the guess.                    !
            !------------------------------------------------------------------------------!
            if (ciz < cimin .or. ciz > cimax) then
               !----- This guess went off-bounds, we give up on Newton's. -----------------!
               converged = .false.
               exit newloop
               !---------------------------------------------------------------------------!
            elseif (converged) then
               !----- Converged, find the root as the mid-point. --------------------------!
               ci  = 5.d-1 * (cia+ciz)
               exit newloop
               !---------------------------------------------------------------------------!
            else
               !----- Not there yet, update the function evaluation and the derivative. ---!
               call iter_solver_step(ib,.true.,ciz,fun,deriv)
               if (fun == 0.d0) then 
                  !----- We have actually hit the jackpot, the answer is ciz. -------------!
                  ci        = ciz
                  converged = .true.
                  exit newloop
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do newloop
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      if (.not. converged) then 

         !---------------------------------------------------------------------------------!
         !     If we have reached this point it means that Newton's method failed.  We     !
         ! switch to the Regula Falsi instead.  The first step is to find two solutions of !
         ! the opposite side.                                                              !
         !---------------------------------------------------------------------------------!
         if (ciz < cimin .or. ciz > cimax) then
            !----- The guess is outside the range, discard it and start over. -------------!
            cia      = sqrt(cimin*cimax)
            call iter_solver_step(ib,.false.,cia,funa,deriv)
            zside    = .false.
         else if (funa * fun < 0.d0) then
            funz     = fun
            zside    = .true.
         else
            zside    = .false.
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     If we still don't have a good second guess, look for one.                   !
         !---------------------------------------------------------------------------------!
         if (.not. zside) then
            !----- Find the extrapolation term to try to hit the other side. --------------!
            delta = 1.d-2 * min(abs(cimax-cia),abs(cimin-cia),abs(cimax-cimin))
            !------------------------------------------------------------------------------!



            !----- First attempt. ---------------------------------------------------------!
            ciz   = cia + delta
            zside = .false.
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Second guess seeker loop.                                                !
            !------------------------------------------------------------------------------!
            hitmin = .false.
            hitmax = .false.
            itb    = 0
            zgssloop: do
               itb = itb + 1
               if (hitmin .and. hitmax) then
                  !------------------------------------------------------------------------!
                  !     We searched through the entire range of ci, and we couldn't find   !
                  ! any pair of roots of the opposite sign, it's likely that there is no   !
                  ! solution, so we give up.                                               !
                  !------------------------------------------------------------------------!
                  return
               end if
               !---------------------------------------------------------------------------!

               ciz = cia + dble((-1)**itb * (itb+3)/2) * delta
               if (ciz < cimin) then
                   !-----------------------------------------------------------------------!
                   !    We have hit the minimum.  Force it to be the minimum, and make the !
                   ! hitmin flag true.                                                     !
                   !-----------------------------------------------------------------------!
                   ciz    = cimin
                   hitmin = .true.
               elseif (ciz > cimax) then
                   ciz    = cimax
                   hitmax = .true.
               end if
               !---------------------------------------------------------------------------!
               !     Compute the function evaluate and check signs.                        !
               !---------------------------------------------------------------------------!
               call iter_solver_step(ib,.false.,ciz,funz,deriv)
               zside = funa*funz < 0.d0
               if (zside) exit zgssloop
            end do zgssloop
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     If we have reached this point, it means that there is a root that hasn't    !
         ! been found yet, but at least we know that it is between cia and ciz.  Use the   !
         ! modified Regula Falsi (Illinois) method to find it.                             !
         !---------------------------------------------------------------------------------!
         fpoloop: do itb=itn,maxfpofl
            ci = (funz * cia - funa * ciz) / (funz-funa)

            !------------------------------------------------------------------------------!
            !     Now that we updated the guess, check whether they are really close.      !
            ! In case they are, this means that it converged, we can use this as our root. !
            !------------------------------------------------------------------------------!
            converged = 2.d0 * abs(ci - cia) < tolerfl8 * (abs(cia)+abs(ciz))
            if (converged) then
               exit fpoloop
            end if
            !------------------------------------------------------------------------------!


            !----- Find the new function evaluation. --------------------------------------!
            call iter_solver_step(ib,.false.,ci,fun,deriv)
            !------------------------------------------------------------------------------!


            !------ Define the new interval based on the intermediate value theorem. ------!
            if (fun == 0.d0) then
               !----- We have actually hit the jackpot, the answer is ciz. ----------------!
               converged = .true.
               exit fpoloop
               !---------------------------------------------------------------------------!
            elseif (fun*funa < 0.d0 ) then
               ciz = ci
               funz  = fun
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa=funa * 5.d-1
               !----- We just updated zside, set zside to true. ---------------------------!
               zside = .true.
               !---------------------------------------------------------------------------!
            else
               cia = ci
               funa   = fun
               !----- If we are updating aside again, modify aside (Illinois method) ------!
               if (.not. zside) funz=funz * 5.d-1
               !----- We just updated aside, set zside to false. --------------------------!
               zside = .false.
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do fpoloop
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     In case it converged, we can compute the answer.                               !
      !------------------------------------------------------------------------------------!
      if (converged) then
         !----- 1. Intercellular CO2, we just utilise the answer we've just got. ----------!
         answer%lint_co2      = ci
         !----- 3. Compute the CO2 demand. ------------------------------------------------!
         answer%co2_demand    = calc_co2_demand(ib,answer%lint_co2)
         !----- 4. Compute the actual stomatal conductance of water and CO2. --------------!
         answer%stom_cond_h2o = calc_stom_cond_h2o(ib,answer%lint_co2,answer%co2_demand)
         answer%stom_cond_co2 = gsw_2_gsc8 * answer%stom_cond_h2o
         !----- 5. Compute the leaf surface CO2. ------------------------------------------!
         answer%lsfc_co2  = met(ib)%can_co2 - answer%co2_demand / met(ib)%blyr_cond_co2
         !----- 6. Compute the leaf surface vapour mixing ratio. --------------------------!
         answer%lsfc_shv  = ( answer%stom_cond_h2o * met(ib)%lint_shv                      &
                            + met(ib)%blyr_cond_h2o    * met(ib)%can_shv  )                &
                          / ( met(ib)%blyr_cond_h2o + answer%stom_cond_h2o)
         !----- 7. Update the success flag. -----------------------------------------------!
         answer%success   = .true.
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine solve_iterative_case
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will adjust the terms used to compute the CO2 demand function     !
   ! (A_open / A_closed) for the light-, and rubisco-limited cases, as well as the closed  !
   ! stomata case (when most terms are not going to be used in any case).                  !
   !     The CO2 demand function A is defined in the most generic format at M09's          !
   ! equations B1 (open) and B2 (closed).  For the actual solver, M09's equation B2 is     !
   ! substituted by the more generic format, M06's equation B2.                            !
   !---------------------------------------------------------------------------------------!
   subroutine set_co2_demand_params(ib,whichlim)
      use c34constants   , only : aparms       & ! intent(inout)
                                , thispft      & ! intent(in)
                                , met          ! ! intent(in)
      use physiology_coms, only : klowco28     ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                      :: ib        ! Multithread ID
      character(len=*), intent(in) :: whichlim  ! A flag telling which case we are
                                                !   about to solve
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Define the parameters based on this call.                                       !
      !------------------------------------------------------------------------------------!
      select case(trim(whichlim))
      case ('CLOSED')
         !----- Closed stomata case, or night time.  These are the same for C3 and C4. ----!
         aparms(ib)%rho   = 0.d0
         aparms(ib)%sigma = 0.d0
         aparms(ib)%xi    = 0.d0
         aparms(ib)%tau   = 1.d0
         aparms(ib)%nu    = - aparms(ib)%leaf_resp
      case default
         !---------------------------------------------------------------------------------!
         !     Open stomata case, so now we distinguish between C3 and C4 as their         !
         ! functional forms are different.                                                 !
         !---------------------------------------------------------------------------------!
         select case (thispft(ib)%photo_pathway)
         case (3)
            !------------------------------------------------------------------------------!
            !     C3 case.  Decide which case to solve.                                    !
            !------------------------------------------------------------------------------!
            select case (trim(whichlim))
            case ('LIGHT')
               !---- RubP regeneration (aka Light) limitation case. -----------------------!
               aparms(ib)%rho   =  aparms(ib)%jact
               aparms(ib)%sigma = -aparms(ib)%jact * aparms(ib)%compp
               aparms(ib)%xi    = 4.d0
               aparms(ib)%tau   = 8.d0 * aparms(ib)%compp
               aparms(ib)%nu    = -aparms(ib)%leaf_resp
               !---------------------------------------------------------------------------!
            case ('RUBP_SAT')
               !----- RubP-saturated rate of photosynthesis case. -------------------------!
               aparms(ib)%rho   =  aparms(ib)%vm
               aparms(ib)%sigma = -aparms(ib)%vm * aparms(ib)%compp
               aparms(ib)%xi    = 1.d0
               aparms(ib)%tau   = aparms(ib)%kco2 * (1.d0 + met(ib)%can_o2 / aparms(ib)%ko2)
               aparms(ib)%nu    = -aparms(ib)%leaf_resp
               !---------------------------------------------------------------------------!
            case ('TPU')
               !----- Triose phosphate utilisation limitation case. -----------------------!
               aparms(ib)%rho   = 0.d0
               aparms(ib)%sigma = 3.d0 * aparms(ib)%tpm
               aparms(ib)%xi    = 0.d0
               aparms(ib)%tau   = 1.d0
               aparms(ib)%nu    = -aparms(ib)%leaf_resp
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!
         case (4)
            !------------------------------------------------------------------------------!
            !     C4 case.  Decide which case to solve.                                    !
            !------------------------------------------------------------------------------!
            select case(trim(whichlim))
            case ('LIGHT')
               !---- RubP regeneration (aka Light) limitation case. -----------------------!
               aparms(ib)%rho   = 0.d0
               aparms(ib)%sigma = aparms(ib)%jact
               aparms(ib)%xi    = 0.d0
               aparms(ib)%tau   = 4.d0
               aparms(ib)%nu    = - aparms(ib)%leaf_resp
               !---------------------------------------------------------------------------!

            case ('RUBP_SAT')
               !----- RubP-saturated rate of photosynthesis case. -------------------------!
               aparms(ib)%rho   = 0.d0
               aparms(ib)%sigma = aparms(ib)%vm
               aparms(ib)%xi    = 0.d0
               aparms(ib)%tau   = 1.d0
               aparms(ib)%nu    = -aparms(ib)%leaf_resp
               !---------------------------------------------------------------------------!

            case ('CO2')
               !----- CO2-limited for low CO2 concentration case. -------------------------!
               aparms(ib)%rho   = klowco28 * aparms(ib)%vm
               aparms(ib)%sigma = 0.d0
               aparms(ib)%xi    = 0.d0
               aparms(ib)%tau   = 1.d0
               aparms(ib)%nu    = -aparms(ib)%leaf_resp
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!
      return
   end subroutine set_co2_demand_params
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine finds the updated function evaluation and derivative for the      !
   ! iterative step for Newton's (or Regula Falsi) method.  The "function" here is a       !
   ! combination of the definition of the water stomatal conductance for open stomata      !
   ! (F96's equation 13), after substituting Ds by a combination of M09's equations B13    !
   ! and B16 so water demand is eliminated(Psi_open), and incorporating F96's equation 14  !
   ! to eliminate the surface carbon.  This function has the property of having one root   !
   ! that corresponds to the intercellular CO2 concentration.                              !
   !  The logical variable is used to decide whether to compute the derivatives or not.    !
   ! They are necessary only when it is a Newton's method call.                            !
   !---------------------------------------------------------------------------------------!
   subroutine iter_solver_step(ib,newton,lint_co2,fun,deriv)
      use c34constants, only : thispft & ! intent(in)
                             , aparms  & ! intent(in)
                             , met     ! ! intent(in)
      
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                   :: ib                  ! Multithread ID
      logical     , intent(in)  :: newton              ! Newton's method step  [       T|F]
      real(kind=8), intent(in)  :: lint_co2            ! Intercell. CO2 conc.  [   mol/mol]
      real(kind=8), intent(out) :: fun                 ! Function evaluation   [       ---]
      real(kind=8), intent(out) :: deriv               ! Derivative.           [   mol/mol]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)              :: stom_cond_h2o       ! Water conductance     [  mol/m2/s]
      real(kind=8)              :: co2_demand          ! CO2 demand            [  mol/m2/s]
      real(kind=8)              :: co2_demand_prime    ! Deriv. of CO2 demand  [1/mol/m2/s]
      real(kind=8)              :: stom_cond_h2o_prime ! Deriv. of H2O cond.   [1/mol/m2/s]
      real(kind=8)              :: efun1               ! 1st term
      real(kind=8)              :: efun2               ! 2nd term
      real(kind=8)              :: efun3               ! 3rd term
      real(kind=8)              :: eprime1             ! Derivative of the 1st term
      real(kind=8)              :: eprime2             ! Derivative of the 2nd term
      real(kind=8)              :: eprime3             ! Derivative of the 3rd term
      !------------------------------------------------------------------------------------!

      !----- Find the CO2 demand. ---------------------------------------------------------!
      co2_demand       = calc_co2_demand(ib,lint_co2)
      !----- Find the stomatal conductance of water. --------------------------------------!
      stom_cond_h2o    = calc_stom_cond_h2o(ib,lint_co2,co2_demand)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the function components, then the function evaluation.                    !
      !------------------------------------------------------------------------------------!
      efun1 = (stom_cond_h2o - thispft(ib)%b) / (thispft(ib)%m * co2_demand)
      efun2 = (met(ib)%can_co2 - aparms(ib)%compp - co2_demand/met(ib)%blyr_cond_co2)
      efun3 = 1.d0 + ( met(ib)%blyr_cond_h2o * (met(ib)%lint_shv - met(ib)%can_shv)        &
                     / (thispft(ib)%d0 * (met(ib)%blyr_cond_h2o + stom_cond_h2o)))
      fun   = efun1 * efun2 * efun3 - 1.d0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case this is a Newton's step, we must also compute the derivatives.  Other- !
      ! wise, we assign any number, but it should not be used.                             !
      !------------------------------------------------------------------------------------!
      if (newton) then
         !----- CO2 demand. ---------------------------------------------------------------!
         co2_demand_prime    = calc_co2_demand_prime(ib,lint_co2,co2_demand)
         !----- stomatal conductance of water. --------------------------------------------!
         stom_cond_h2o_prime = calc_stom_cond_h2o_prime(ib,lint_co2,stom_cond_h2o          &
                                                       ,co2_demand,co2_demand_prime)
         !----- Function components. ------------------------------------------------------!
         eprime1 = ( stom_cond_h2o_prime * co2_demand                                      &
                   - co2_demand_prime * (stom_cond_h2o - thispft(ib)%b) )                  &
                 / (thispft(ib)%m * co2_demand * co2_demand)

         eprime2 = - co2_demand_prime / met(ib)%blyr_cond_co2

         eprime3 = - (efun3 - 1.d0) *  stom_cond_h2o_prime                                 &
                 / ( met(ib)%blyr_cond_h2o + stom_cond_h2o )

         deriv   = eprime1 * efun2   * efun3                                               &
                 + efun1   * eprime2 * efun3                                               &
                 + efun1   * efun2   * eprime3
      else
         deriv = 0.d0
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine iter_solver_step
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the CO2 demand given the intercellular CO2 concentration.  !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function calc_co2_demand(ib,lint_co2)
      use c34constants, only : aparms ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: ib       ! Multithread ID
      real(kind=8), intent(in) :: lint_co2 ! Intercellular CO2 concentration    [  mol/mol]
      !------------------------------------------------------------------------------------!

      calc_co2_demand = (aparms(ib)%rho * lint_co2 + aparms(ib)%sigma)                     &
                      / (aparms(ib)%xi  * lint_co2 + aparms(ib)%tau  ) + aparms(ib)%nu
      return
   end function calc_co2_demand
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the derivative of the CO2 regarding the intercellular CO2  !
   ! concentration.                                                                        !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function calc_co2_demand_prime(ib,lint_co2,co2_demand)
      use c34constants, only : aparms ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: ib         ! Multithread ID
      real(kind=8), intent(in) :: lint_co2   ! Intercellular CO2 concentration  [  mol/mol]
      real(kind=8), intent(in) :: co2_demand ! CO2 demand                       [ mol/m2/s]
      !------------------------------------------------------------------------------------!
      
      calc_co2_demand_prime = ( aparms(ib)%rho                                             &
                              - aparms(ib)%xi * (co2_demand - aparms(ib)%nu))              &
                            / ( aparms(ib)%xi  * lint_co2 + aparms(ib)%tau  )
      return
   end function calc_co2_demand_prime
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the stomatal conductance of water given the intercellular  !
   ! CO2 concentration and the CO2 demand.                                                 !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function calc_stom_cond_h2o(ib,lint_co2,co2_demand)
      use c34constants   , only : met        ! ! intent(in)
      use physiology_coms, only : gsw_2_gsc8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: ib         ! Multithread ID
      real(kind=8), intent(in) :: lint_co2   ! Intercellular CO2 concentration  [  mol/mol]
      real(kind=8), intent(in) :: co2_demand ! CO2 demand                       [ mol/m2/s]
      !------------------------------------------------------------------------------------!

      calc_stom_cond_h2o = met(ib)%blyr_cond_co2 * co2_demand                              &
                         / (gsw_2_gsc8 * ( (met(ib)%can_co2 - lint_co2)                    &
                                         * met(ib)%blyr_cond_co2  - co2_demand ) )

      return
   end function calc_stom_cond_h2o
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the stomatal conductance of water derivative given the     !
   ! intercellular CO2 concentration, the water stomatal conductance, and the CO2 demand   !
   ! and its derivative.                                                                   !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function calc_stom_cond_h2o_prime(ib,lint_co2,stom_cond_h2o,co2_demand     &
                                                 ,co2_demand_prime)
      use c34constants   , only : met        ! ! intent(in)
      use physiology_coms, only : gsw_2_gsc8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: ib               ! Multithread ID
      real(kind=8), intent(in) :: lint_co2         ! Intercell. CO2 conc.      [   mol/mol]
      real(kind=8), intent(in) :: stom_cond_h2o    ! Water                     [  mol/m2/s]
      real(kind=8), intent(in) :: co2_demand       ! CO2 demand                [  mol/m2/s]
      real(kind=8), intent(in) :: co2_demand_prime ! Derivative of CO2 demand  [1/mol/m2/s]
      !------------------------------------------------------------------------------------!

      calc_stom_cond_h2o_prime = stom_cond_h2o                                             &
                               * ( co2_demand_prime / co2_demand                           &
                                 + (met(ib)%blyr_cond_co2 + co2_demand_prime)              &
                                 / (gsw_2_gsc8 * ( (met(ib)%can_co2 - lint_co2)            &
                                                 * met(ib)%blyr_cond_co2 - co2_demand) ))

      return
   end function calc_stom_cond_h2o_prime
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine computes the minimum and maximum intercellular carbon dioxide    !
   ! concentration that we should seek the solution.  Both the function from which we seek !
   ! the root and the stomatal conductance function have singularities, so we don't want   !
   ! the intercellular carbon dioxide to cross these singularities because the root-find-  !
   ! ing method assumes continuity.                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine find_lint_co2_bounds(ib,cimin,cimax,bounded)
      use c34constants   , only : thispft           & ! intent(in)
                                , aparms            & ! intent(in)
                                , met               ! ! intent(in)
      use physiology_coms, only : gbw_2_gbc8        & ! intent(in)
                                , c34smin_lint_co28 ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                    :: ib      ! Multithread ID
      real(kind=8), intent(out)  :: cimin   ! Minimum intercellular CO2         [  mol/mol]
      real(kind=8), intent(out)  :: cimax   ! Maximum intercellular CO2         [  mol/mol]
      logical     , intent(out)  :: bounded ! This problem is bounded           [      T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)               :: aquad   ! Quadratic coefficient             [      ---]
      real(kind=8)               :: bquad   ! Linear coefficient                [  mol/mol]
      real(kind=8)               :: cquad   ! Intercept                         [mol2/mol2]
      real(kind=8)               :: ciAo    ! Ci at singularity where Aopen=0   [  mol/mol]
      real(kind=8)               :: cigsw   ! Ci at sing. where gbc(ca-ci)=Ao   [  mol/mol]
      real(kind=8)               :: ciQ     ! Ci at singularity where qi-qs=-D0 [  mol/mol]
      real(kind=8)               :: xtmp    ! variable for ciQ equation
      real(kind=8)               :: ytmp    ! variable for ciQ
      real(kind=8)               :: ztmp    ! variable for ciQ
      real(kind=8)               :: wtmp    ! variable for ciQ
      real(kind=8)               :: ciroot1 ! 1st root for the quadratic eqn.   [  mol/mol]
      real(kind=8)               :: ciroot2 ! 2nd root for the quadratic eqn.   [  mol/mol]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! First case: This check will find when Aopen goes to 0., which causes a singularity !
      ! in the function of which we are looking for a root.                                !
      !------------------------------------------------------------------------------------!
      if (aparms(ib)%rho == 0.d0 .and. aparms(ib)%xi == 0.d0) then
         ciAo = c34smin_lint_co28
      else
         ciAo = - (aparms(ib)%tau * aparms(ib)%nu + aparms(ib)%sigma)                      &
              / (aparms(ib)%xi  * aparms(ib)%nu + aparms(ib)%rho  )
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! Second case: This check finds which ci causes the terms [gbc (ca -ci) - Aopen] to  !
      ! be 0.   This will cause a singularity in the gsw function.                         !
      !------------------------------------------------------------------------------------!
      !----- 1. Define the coefficients for the quadratic equation. -----------------------!
      aquad = met(ib)%blyr_cond_co2 * aparms(ib)%xi
      bquad = aparms(ib)%xi * (aparms(ib)%nu - met(ib)%blyr_cond_co2 * met(ib)%can_co2 )   &
            + met(ib)%blyr_cond_co2 * aparms(ib)%tau + aparms(ib)%rho
      cquad = aparms(ib)%tau * (aparms(ib)%nu - met(ib)%blyr_cond_co2 * met(ib)%can_co2 )  &
            + aparms(ib)%sigma
      !----- 2. Solve the quadratic equation. ---------------------------------------------!
      call solve_quadratic8(aquad,bquad,cquad,-discard,ciroot1,ciroot2)
      !----- 3. Save the largest value (or flag them as invalid when both are invalid). ---!
      cigsw=max(ciroot1, ciroot2)
      if (cigsw == -discard) cigsw = discard
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      ! Third case: This will find the intercellular CO2 mixing ratio that causes the      !
      ! qi-qs term to be equal to -D0, which also creates a singularity in the stomatal    !
      ! conductance.                                                                       !
      !------------------------------------------------------------------------------------!
      !----- 1. Find some auxiliary variables. --------------------------------------------!
      xtmp = met(ib)%blyr_cond_h2o                                                         &
           * ( met(ib)%can_shv - met(ib)%lint_shv - thispft(ib)%d0 ) / thispft(ib)%d0
      ytmp = met(ib)%blyr_cond_co2 + xtmp * gbw_2_gbc8
      ztmp = xtmp * gbw_2_gbc8 * met(ib)%blyr_cond_co2
      wtmp = ztmp * met(ib)%can_co2 - ytmp * aparms(ib)%nu

      !----- 2. Define the coefficients for the quadratic equation. -----------------------!
      aquad = ztmp * aparms(ib)%xi
      bquad = ytmp * aparms(ib)%rho - aparms(ib)%xi * wtmp + ztmp * aparms(ib)%tau
      cquad = - aparms(ib)%tau * wtmp + ytmp * aparms(ib)%sigma

      !----- 3. Solve the quadratic equation. ---------------------------------------------!
      call solve_quadratic8(aquad,bquad,cquad,-discard,ciroot1,ciroot2)

      !------------------------------------------------------------------------------------!
      !      4. Save the largest of the values.  In case both were discarded, we switch it !
      ! to the positive discard so this will never be chosen.                              !
      !------------------------------------------------------------------------------------!
      ciQ=max(ciroot1, ciroot2)
      if (ciQ == -discard) ciQ = discard
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Make sure that there is no singularity within Cimin and Cimax.  From previous  !
      ! tests, we know that cimin is the one associated with the case in which Aopen goes  !
      ! to zero, and the maximum is the minimum between gsw case, q case, and the canopy   !
      ! air CO2.                                                                           !
      !------------------------------------------------------------------------------------!
      ciroot1 = max(c34smin_lint_co28,ciAo)
      ciroot2 = min(ciQ, cigsw,met(ib)%can_co2)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The actual bounds are slightly squeezed so the edges will not be at the the    !
      ! singularities.                                                                     !
      !------------------------------------------------------------------------------------!
      cimin = ciroot1 + 1.d-5 * (ciroot2 - ciroot1)
      cimax = ciroot2 - 1.d-5 * (ciroot2 - ciroot1)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    In the last part we make sure that the two guesses are positively defined and   !
      ! make sense.  If both bounds are non-positive, or if the minimum ci is greater than !
      ! the maximum cimax, the bounds make no sense, and the solution is not bounded.      !
      !------------------------------------------------------------------------------------!
      if (cimin > cimax) then
         cimin = c34smin_lint_co28
         cimax = c34smin_lint_co28
         bounded = .false.
      elseif (cimin <= 0.d0 .and. cimax <= c34smin_lint_co28) then
         cimin = c34smin_lint_co28
         cimax = c34smin_lint_co28
         bounded = .false.
      elseif (cimin <= 0.d0) then
         cimin = c34smin_lint_co28
         bounded = .true.
      else
         bounded = .true.
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine find_lint_co2_bounds
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This is the Arrhenius function, written as a function of a reference temper-     !
   ! ature, given by tphysref8.  The output variable will have the same unit as the        !
   ! reference value.                                                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function arrhenius(temp,refval,hor)
      use physiology_coms, only : tphysrefi8 ! ! intent(in)  
      use consts_coms    , only : lnexp_min8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp   ! Temperature                              [    K]
      real(kind=8), intent(in) :: refval ! Reference value at T=Tphysref            [  any]
      real(kind=8), intent(in) :: hor    ! "Activation energy" / Gas constant       [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: lnexp  ! Term that will go to the exponential     [  ---]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the term that goes to the exponential term, and check its size.  This is  !
      ! to avoid floating point exceptions due to overflow or underflow.                   !
      !------------------------------------------------------------------------------------!
      lnexp = hor * (tphysrefi8 - 1.d0/temp)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If the exponential factor is tiny, make it zero, otherwise compute the actual  !
      ! function.                                                                          !
      !------------------------------------------------------------------------------------!
      if (lnexp < lnexp_min8) then
         arrhenius = 0.d0
      else
         arrhenius = refval * exp(lnexp)
      end if
      !------------------------------------------------------------------------------------!

      return
   end function arrhenius
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This is the temperature-dependence for the physiology parameters, represented as !
   ! in function, written as a in Collatz et al. (1991).                                   !
   ! The output variable will have the same unit as the reference value at tphysref.       !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function collatz(temp,refval,q10)
      use physiology_coms, only : tphysref8  & ! intent(in)
                                , fcoll8     ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp      ! Temperature                           [    K]
      real(kind=8), intent(in) :: refval    ! Pre-factor coefficient                [  any]
      real(kind=8), intent(in) :: q10       ! Exponential coefficient               [    K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     If the exponential factor is tiny, make it zero, otherwise compute the actual  !
      ! function.                                                                          !
      !------------------------------------------------------------------------------------!
      collatz = refval * q10 ** (fcoll8 * (temp - tphysref8))
      !------------------------------------------------------------------------------------!

      return
   end function collatz
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Find the minimum amount of radiation for which we will consider daytime.  This    !
   ! cannot be a constant because it depends on the PFT and the environmental conditions.  !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function find_twilight_min(ib)
      use c34constants   , only : aparms              & ! intent(in)
                                , thispft             & ! intent(in)
                                , met                 ! ! intent(in)
      use physiology_coms, only : iphysiol            & ! intent(in)
                                , par_lightcompp_max8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: ib !> Multithread ID
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: j0 !> Critical electron transport [mol/m2/s]
      !------------------------------------------------------------------------------------!


      !----- Electron transport that makes Ao = 0. ----------------------------------------!
      j0 = 4.d0 * aparms(ib)%leaf_resp * (met(ib)%can_co2 + 2.d0 * aparms(ib)%compp)       &
                                       / (met(ib)%can_co2 -        aparms(ib)%compp)
      !------------------------------------------------------------------------------------!


      select case (iphysiol)
      case (0,2)
         !----- Linear equation relating PAR and J (Foley/Collatz style). -----------------!
         find_twilight_min = 2.5d-1 * j0 / aparms(ib)%alpha
         !---------------------------------------------------------------------------------!
      case (1,3)
         !---------------------------------------------------------------------------------!
         !      J is defined from quadratic equation (van Caemmerer/CLM 4.5 style).        !
         !---------------------------------------------------------------------------------!

         !----- Compare J0 with Jmax, to avoid singularities. -----------------------------!
         if (j0 < aparms(ib)%jm) then
            !------------------------------------------------------------------------------!
            !    Safe case, use the inverse of the quadratic equation to find the light    !
            ! compensation point.                                                          !
            !------------------------------------------------------------------------------!
            find_twilight_min = 2.d0 * j0 * (aparms(ib)%jm - thispft(ib)%curvpar * j0)     &
                              / ( thispft(ib)%phi_psII * (aparms(ib)%jm - j0) )
            !------------------------------------------------------------------------------!
         else
            !----- Set to a large number so it cannot be selected. ------------------------!
            find_twilight_min = max(par_lightcompp_max8,1.1d0*met(ib)%par)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end select

      return
   end function find_twilight_min
   !=======================================================================================!
   !=======================================================================================!
end module farq_leuning
!==========================================================================================!
!==========================================================================================!

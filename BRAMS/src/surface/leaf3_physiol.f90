!==========================================================================================!
!==========================================================================================!
!    This module contains some structures and subroutines used by the model to solve the   !
! photosynthesis.  The references for this model are given by:                             !
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
module leaf3_physiol

   implicit none

   !---------------------------------------------------------------------------------------!
   !     This is a flag used in various sub-routines and functions and denote that we      !
   ! should ignore the result.                                                             !
   !---------------------------------------------------------------------------------------!
   real(kind=8), parameter :: discard = huge(1.d0)
   !---------------------------------------------------------------------------------------!



   !----- Small number, to avoid singularities. -------------------------------------------!
   real(kind=8), parameter :: tiny_offphys = 1.d-20
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    This is the tolerance for iterative methods.  Some of the functions are quite flat !
   ! so it is a good idea to use a somewhat more strict tolerance than the ones used in    !
   ! therm_lib8.                                                                           !
   !---------------------------------------------------------------------------------------!
   real(kind=8), parameter :: tolerfl8 = 1.d-10
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Since this is more strict than therm_lib, the maximum number of attempts for the   !
   ! Regula Falsi method should be also increased.                                         !
   !---------------------------------------------------------------------------------------!
   integer     , parameter :: maxfpofl = 320
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     This structure contains some PFT-dependent physiological and empirical constants  !
   ! for the photosynthesis model.  The units inside the photosynthesis model are not the  !
   ! same as the input data.  Everything that is in µmol becomes mol here, and everything  !
   ! that is in °C becomes Kelvin here.                                                    !
   !---------------------------------------------------------------------------------------!
   type leaf3_fl_parm
      integer      :: photo_pathway ! The photosynthetic pathway                 [     ---]
      real(kind=8) :: D0            ! Threshold to close stomata when it's dry   [ mol/mol]
      real(kind=8) :: m             ! Empirical const. for conductance calc.     [     ---]
      real(kind=8) :: b             ! Cuticular conductance                      [mol/m²/s]
      real(kind=8) :: vm0           ! Ref. value for maximum rubisco activity    [mol/m²/s]
      real(kind=8) :: vm_low_temp   ! Low temp. threshold for rapid decline      [       K]
      real(kind=8) :: vm_high_temp  ! High temp. threshold for rapid decline     [       K]
      real(kind=8) :: vm_q10        ! Vm base for C91 physiology                 [      --]
      real(kind=8) :: vm_decay_e    ! Decay factor for Vm correction as in F96   [     1/K]
      real(kind=8) :: rd0           ! Term used for respiration                  [mol/m²/s]
      real(kind=8) :: rd_low_temp   ! Low temp. threshold for rapid decline      [       K]
      real(kind=8) :: rd_high_temp  ! High temp. threshold for rapid decline     [       K]
      real(kind=8) :: rd_q10        ! Rd base for C91 physiology                 [      --]
      real(kind=8) :: rd_decay_e    ! Decay factor for Rd correction as in F96   [     1/K]
   end type leaf3_fl_parm
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This structure contains the canopy air (plus some leaf) conditions.  The short    !
   ! name given in David Medvigy's thesis is given in parenthesis, and the units used in   !
   ! the solver are given in brackets.                                                     !
   !---------------------------------------------------------------------------------------!
   type leaf3_met_inp
      real(kind=8) :: can_o2        ! Canopy air oxygen concentration   ([O2] ) [  mol/mol]
      real(kind=8) :: can_shv       ! Canopy air specific humidity      (ea   ) [  mol/mol]
      real(kind=8) :: can_co2       ! Canopy air CO2 mixing ratio       (Ca   ) [  mol/mol]
      real(kind=8) :: can_rhos      ! Canopy air density                (dens ) [    kg/m³]
      real(kind=8) :: can_prss      ! Canopy air pressure               (press) [       Pa]
      real(kind=8) :: leaf_temp     ! Leaf temperature                  (T    ) [        K]
      real(kind=8) :: lint_shv      ! Leaf intercell. vap. mix. ratio   (eL   ) [  mol/mol]
      real(kind=8) :: par           ! Photosyntethically active rad.    (L    ) [ mol/m²/s]
      real(kind=8) :: blyr_cond_co2 ! Leaf bnd. layer conduct. for CO2  (gbc  ) [ mol/m²/s]
      real(kind=8) :: blyr_cond_h2o ! Leaf bnd. layer conduct. for H2O  (gbw  ) [ mol/m²/s]
   end type leaf3_met_inp
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     This structure contains the variables used to define the optimal leaf-level       !
   ! carbon dioxide demand of photosynthesis.  The names follow the original notation      !
   ! found in M06 equations 2.62, 2.63, 2.67, and B.2.                                     !
   !---------------------------------------------------------------------------------------!
   type leaf3_co2_dem_vars
      real(kind=8) :: alpha     ! Quantum yield                                 [      ---]
      real(kind=8) :: vm        ! Rubisco capacity to perform carboxylase       [ mol/m²/s]
      real(kind=8) :: compp     ! Compens. pt. for gross photos.    (Gamma*)    [  mol/mol]
      real(kind=8) :: kco2      ! Michaelis-Mentel coefficient for CO2          [  mol/mol]
      real(kind=8) :: ko2       ! Michaelis-Mentel coefficient for O2           [  mol/mol]
      real(kind=8) :: leaf_resp ! Leaf respiration                              [ mol/m²/s]
      real(kind=8) :: sigma     ! sigma from equation B.2                       [ mol/m²/s]
      real(kind=8) :: rho       ! rho from equation B.2                         [ mol/m²/s]
      real(kind=8) :: xi        ! coefficient for Ci in the denominator         [  mol/mol]
      real(kind=8) :: tau       ! tau from equation B.2                         [  mol/mol]
      real(kind=8) :: nu        ! The negative of leaf respiration              [ mol/m²/s]
   end type leaf3_co2_dem_vars
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     This structure contains the variables that will be found for both open and closed !
   ! stomata cases, and in case of the open, for the Rubisco- and light-limited cases,     !
   ! plus the CO2 case for C4 photosynthesis.                                              !
   !---------------------------------------------------------------------------------------!
   type leaf3_fl_sol
      real(kind=8) :: lsfc_shv      ! Leaf surface specific humidity     (es )  [  mol/mol]
      real(kind=8) :: lsfc_co2      ! Leaf surface CO2 mixing ratio      (Cs )  [  mol/mol]
      real(kind=8) :: lint_co2      ! Leaf internal CO2 mixing ratio     (Ci )  [  mol/mol]
      real(kind=8) :: co2_demand    ! Carbon demand                      (A  )  [ mol/m²/s]
      real(kind=8) :: stom_cond_h2o ! Stomatal water conductance         (gsw)  [ mol/m²/s]
      real(kind=8) :: stom_cond_co2 ! Stomatal CO2 conductance           (gsc)  [ mol/m²/s]
   end type leaf3_fl_sol
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    These are the structures used throughout the solution.  They are kept here so they !
   ! can be accessed by any subroutine.                                                    !
   !---------------------------------------------------------------------------------------!
   type(leaf3_fl_parm     ) :: thispft
   type(leaf3_met_inp     ) :: met
   type(leaf3_co2_dem_vars) :: aparms
   type(leaf3_fl_sol      ) :: stopen
   type(leaf3_fl_sol      ) :: stclosed
   type(leaf3_fl_sol      ) :: rubiscolim
   type(leaf3_fl_sol      ) :: co2lim
   type(leaf3_fl_sol      ) :: lightlim
   !---------------------------------------------------------------------------------------!



   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !      This is the main driver for the photosynthesis model.                            !
   !---------------------------------------------------------------------------------------!
   subroutine leaf3_farquhar_leuning(can_prss,can_rhos,can_shv,can_co2,leaf_class,leaf_par &
                                    ,leaf_temp,lint_shv,vm0,rd0,a_open,a_closed,gsw_open   &
                                    ,gsw_closed,leaf_resp)
      use leaf_coms      , only : gbw                      & ! intent(in)
                                , pathway                  & ! intent(in)
                                , leaf_supply              & ! intent(in)
                                , dr_gamma                 & ! intent(in)
                                , gsw_0                    & ! intent(in)
                                , gsw_m                    & ! intent(in)
                                , quantum_yield            & ! intent(in)
                                , vm0_qten                 & ! intent(in)
                                , vm0_dec                  & ! intent(in)
                                , rd0_qten                 & ! intent(in)
                                , rd0_dec                  & ! intent(in)
                                , phys_low_temp            & ! intent(in)
                                , phys_high_temp           & ! intent(in)
                                , gbw                      & ! intent(in)
                                , stom_side                & ! intent(in)
                                , c34smin_lint_co28        & ! intent(in)
                                , c34smax_lint_co28        & ! intent(in)
                                , gbh_2_gbw8               & ! intent(in)
                                , gbw_2_gbc8               & ! intent(in)
                                , o2_ref8                  & ! intent(in)
                                , compp_refval8            & ! intent(in)
                                , compp_q108               & ! intent(in)
                                , kco2_refval8             & ! intent(in)
                                , kco2_q108                & ! intent(in)
                                , ko2_refval8              & ! intent(in)
                                , ko2_q108                 & ! intent(in)
                                , stom_side                ! ! intent(in)
      use rconstants     , only : mmh2oi8                  & ! intent(in)
                                , mmh2o8                   & ! intent(in)
                                , mmdryi8                  & ! intent(in)
                                , mmdry8                   & ! intent(in)
                                , ep8                      & ! intent(in)
                                , epi8                     & ! intent(in)
                                , t008                     & ! intent(in)
                                , umol_2_mol8              & ! intent(in)
                                , mol_2_umol8              & ! intent(in)
                                , Watts_2_Ein8             & ! intent(in)
                                , lnexp_min8               & ! intent(in)
                                , lnexp_max8               & ! intent(in)
                                , rmol8                    ! ! intent(in)

      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real(kind=4), intent(in)    :: can_prss       ! Canopy air pressure       [       Pa]
      real(kind=4), intent(in)    :: can_rhos       ! Canopy air density        [    kg/m³]
      real(kind=4), intent(in)    :: can_shv        ! Canopy air sp. hum.       [    kg/kg]
      real(kind=4), intent(in)    :: can_co2        ! Canopy air CO2            [ µmol/mol]
      real(kind=4), intent(in)    :: leaf_class     ! Plant functional type     [      ---]
      real(kind=4), intent(in)    :: leaf_par       ! Absorbed PAR              [     W/m²]
      real(kind=4), intent(in)    :: leaf_temp      ! Leaf temperature          [        K]
      real(kind=4), intent(in)    :: lint_shv       ! Leaf interc. sp. hum.     [    kg/kg]
      real(kind=4), intent(in)    :: vm0            ! Reference Vm              [µmol/m²/s]
      real(kind=4), intent(in)    :: rd0            ! Reference Rd              [µmol/m²/s]
      real(kind=4), intent(out)   :: A_open         ! Photosyn. rate (op.)      [µmol/m²/s]
      real(kind=4), intent(out)   :: A_closed       ! Photosyn. rate (cl.)      [µmol/m²/s]
      real(kind=4), intent(out)   :: gsw_open       ! St. cnd. of H2O  (op.)    [  kg/m²/s]
      real(kind=4), intent(out)   :: gsw_closed     ! St. cnd. of H2O  (cl.)    [  kg/m²/s]
      real(kind=4), intent(out)   :: leaf_resp      ! Leaf respiration rate     [µmol/m²/s]
      !------ Local variables. ------------------------------------------------------------!
      integer                     :: nveg           ! Plant functional type     [      ---]
      real(kind=8)                :: stom_side8     ! # of sides with stomata   [      1|2]
      real(kind=8)                :: vm_nocorr      ! Vm  with no correction    [ mol/m²/s]
      real(kind=8)                :: rd_nocorr      ! Rd  with no correction    [ mol/m²/s]
      real(kind=8)                :: lnexplow       ! Low temperature exponent  [      ---]
      real(kind=8)                :: lnexphigh      ! High temperature exponent [      ---]
      real(kind=8)                :: tlow_fun       ! Low temperature bound     [      ---]
      real(kind=8)                :: thigh_fun      ! High temperature bound    [      ---]
      !----- External function. -----------------------------------------------------------!
      real(kind=4)    , external  :: sngloff           ! Safe double -> single precision
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Aliases to some useful variables.                                              !
      !------------------------------------------------------------------------------------!
      nveg       = nint(leaf_class      )
      stom_side8 = dble(stom_side (nveg))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Copy the meteorological forcing to the "met" structure.  Notice that some      !
      ! variables go through unit conversions, and all variables are converted to double   !
      ! precision.                                                                         !
      !------------------------------------------------------------------------------------!
      !----- 1. Variables that remain with the same units. --------------------------------!
      met%leaf_temp     = dble(leaf_temp)
      met%can_rhos      = dble(can_rhos )
      met%can_prss      = dble(can_prss )
      met%can_o2        = o2_ref8
      !----- 2. Convert specific humidity to mol/mol. -------------------------------------!
      met%can_shv       = epi8 * dble(can_shv) 
      !----- 3. Convert CO2 to mol/mol. ---------------------------------------------------!
      met%can_co2       = dble(can_co2) * umol_2_mol8
      !----- 4. Convert W/m2 to mol/m2/s. -------------------------------------------------!
      met%par           = dble(leaf_par) * Watts_2_Ein8
      !------------------------------------------------------------------------------------!
      !  5. Intercellular specific humidity, which is assumed to be at saturation          !
      !     given the leaf temperature.  We convert it to mol/mol.                         !
      !------------------------------------------------------------------------------------!
      met%lint_shv      = epi8 * dble(lint_shv)
      !------------------------------------------------------------------------------------!
      !  6. Find the conductivities for water and carbon.  The input for water is in       !
      !     kg/m²/s, and here we convert to mol/m²/s.  The convertion coefficient from     !
      !     water to carbon dioxide comes from M09's equation B14.  Here we multiply by    !
      !     the effective area for transpiration, depending on whether the leaves of this  !
      !     plant functional type are hypo-stomatous, symmetrical, or amphistomatous.      !
      !------------------------------------------------------------------------------------!
      met%blyr_cond_h2o = dble(gbw)  * mmdryi8 * stom_side8
      met%blyr_cond_co2 = gbw_2_gbc8 * met%blyr_cond_h2o
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Load physiological parameters that are PFT-dependent to the thispft structure. !
      ! Convert all variables to mol and Kelvin, when needed.                              !
      !------------------------------------------------------------------------------------!
      thispft%photo_pathway = pathway(nveg)
      thispft%D0            = dble(leaf_supply   (nveg))
      thispft%b             = dble(gsw_0         (nveg)) * umol_2_mol8
      thispft%m             = dble(gsw_m         (nveg))
      thispft%vm0           = dble(vm0)                  * umol_2_mol8
      thispft%vm_low_temp   = dble(phys_low_temp (nveg)) + t008
      thispft%vm_high_temp  = dble(phys_high_temp(nveg)) + t008
      thispft%vm_q10        = dble(vm0_qten      (nveg))
      thispft%vm_decay_e    = dble(vm0_dec       (nveg))
      thispft%rd0           = dble(rd0)                  * umol_2_mol8
      thispft%rd_low_temp   = dble(phys_low_temp (nveg)) + t008
      thispft%rd_high_temp  = dble(phys_high_temp(nveg)) + t008
      thispft%rd_q10        = dble(rd0_qten      (nveg))
      thispft%rd_decay_e    = dble(rd0_dec       (nveg))
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Compute the photosynthesis and leaf respiration parameters that depend on      !
      ! temperature, according to the parameters defined in the "thispft" structure.  The  !
      ! variables that are defined there are:                                              !
      ! - alpha     - the quantum yield.                                                   !
      ! - Vm        - the maximum capacity of Rubisco to perform the carboxylase function. !
      ! - leaf_resp - the leaf respiration.                                                !
      ! - compp     - the CO2 compensation point for gross photosynthesis (Gamma*)         !
      ! - kco2      - Michaelis-Mentel coefficient for CO2                                 !
      ! - ko2       - Michaelis-Mentel coefficient for O2.                                 !
      !------------------------------------------------------------------------------------!

         !----- Quantum yield.  Use the reference value. ----------------------------------!
         aparms%alpha = dble(quantum_yield (nveg))
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find Vm, the photosynthetic capacity.                                       !
         !---------------------------------------------------------------------------------!
         !----- Find Vm using the Collatz equation, with no correction. -------------------!
         vm_nocorr = leaf3_qten_fun(met%leaf_temp,thispft%vm0,thispft%vm_q10)
         !----- Low temperature. ----------------------------------------------------------!
         lnexplow  = thispft%vm_decay_e * (thispft%vm_low_temp  - met%leaf_temp)
         lnexplow  = max(lnexp_min8,min(lnexp_max8,lnexplow))
         tlow_fun  = 1.d0 +  exp(lnexplow)
         !----- High temperature. ---------------------------------------------------------!
         lnexphigh = thispft%vm_decay_e * (met%leaf_temp - thispft%vm_high_temp)
         lnexphigh = max(lnexp_min8,min(lnexp_max8,lnexphigh))
         thigh_fun = 1.d0 + exp(lnexphigh)
         !------ Correct Vm. --------------------------------------------------------------!
         aparms%vm = vm_nocorr / (tlow_fun * thigh_fun)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find leaf respiration.                                                      !
         !---------------------------------------------------------------------------------!
         !----- Find Rd using the Collatz equation, with no correction. -------------------!
         rd_nocorr = leaf3_qten_fun(met%leaf_temp,thispft%rd0,thispft%rd_q10)
         !---------------------------------------------------------------------------------!
         !----- Low temperature. ----------------------------------------------------------!
         lnexplow  = thispft%rd_decay_e * (thispft%rd_low_temp  - met%leaf_temp)
         lnexplow  = max(lnexp_min8,min(lnexp_max8,lnexplow))
         tlow_fun  = 1.d0 +  exp(lnexplow)
         !----- High temperature. ---------------------------------------------------------!
         lnexphigh = thispft%rd_decay_e * (met%leaf_temp - thispft%rd_high_temp)
         lnexphigh = max(lnexp_min8,min(lnexp_max8,lnexphigh))
         thigh_fun = 1.d0 + exp(lnexphigh)
         !------ Correct Rd. --------------------------------------------------------------!
         aparms%leaf_resp = rd_nocorr / (tlow_fun * thigh_fun)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    If this is a C3 plant, find the compensation point, and the Michaelis-Mentel !
         ! constants for CO2 and O2.  Otherwise, assign zeroes for compensation point and  !
         ! the Michaelis-Mentel constant for CO2.  The oxygen one should have no impact,   !
         ! but we always assign it to avoid divisions by zero.                             !
         !---------------------------------------------------------------------------------!
         select case (thispft%photo_pathway)
         case (3)
            aparms%compp = leaf3_qten_fun(met%leaf_temp,compp_refval8,compp_q108)
            aparms%kco2  = leaf3_qten_fun(met%leaf_temp,kco2_refval8 ,kco2_q108 )
         case (4)
            aparms%compp = 0.d0
            aparms%kco2  = 0.d0
         end select
         aparms%ko2   = leaf3_qten_fun(met%leaf_temp,ko2_refval8,ko2_q108)
         !---------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Call the main solver.                                                          !
      !------------------------------------------------------------------------------------!
      call leaf3_photo_exact_solver()
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Copy the solution to the standard output variables.  Here we convert the values !
      ! back to the standard LEAF-3 units                                                  !
      !------------------------------------------------------------------------------------!
      !----- Carbon demand, convert them to [µmol/m²/s]. ----------------------------------!
      A_closed    = sngloff(stclosed%co2_demand    * mol_2_umol8        , tiny_offphys)
      A_open      = sngloff(stopen%co2_demand      * mol_2_umol8        , tiny_offphys)
      !----- Stomatal resistance, convert the conductances to [kg/m²/s]. ------------------!
      gsw_closed  = sngloff(stclosed%stom_cond_h2o * mmdry8 / stom_side8, tiny_offphys)
      gsw_open    = sngloff(stopen%stom_cond_h2o   * mmdry8 / stom_side8, tiny_offphys)
      !----- Leaf respiration [µmol/m²/s]. ------------------------------------------------!
      leaf_resp   = sngloff(aparms%leaf_resp       * mol_2_umol8        , tiny_offphys)
      !------------------------------------------------------------------------------------!


      return

   end subroutine leaf3_farquhar_leuning
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine is the main driver for the C3 photosynthesis.                     !
   !---------------------------------------------------------------------------------------!
   subroutine leaf3_photo_exact_solver()
      use leaf_coms, only : c34smax_gsw8     ! ! intent(in)
      implicit none
      !------ Local variables. ------------------------------------------------------------!
      logical                   :: success          ! The solver succeeded.      [     T|F]
      real(kind=8)              :: par_twilight_min ! Minimum daytime radiation  [mol/m2/s]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Initialise the parameters to compute the carbon demand for the case where the !
      ! stomata are closed.  We call the solver for the specific case in which the stomata !
      ! are closed, because in this case neither the carbon demand nor the stomatal water  !
      ! conductance depend on the intercellular CO2 mixing ratio.                          !
      !------------------------------------------------------------------------------------!
      call leaf3_demand_parm_co2('CLOSED')
      call leaf3_solve_aofixed(stclosed,success)
      if (.not. success) then
         call abort_run ('Solution failed for closed case'                                 &
                        ,'leaf3_photo_exact_solver','leaf3_physiol.f90')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    First we check whether it is at least dawn or dusk.  In case it is not, no      !
      ! photosynthesis should happen, so we copy the closed case stomata values to the     !
      ! open case.  Limit_flag becomes 0, which is the flag for night time limitation.     !
      !------------------------------------------------------------------------------------!
      par_twilight_min = leaf3_twilight_min()
      if (met%par < par_twilight_min) then
         call leaf3_copy_solution(stclosed,stopen)
         return
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
      !   1. The light-limited (aka PAR) case.                                             !
      !------------------------------------------------------------------------------------!
      !----- Update the CO2 demand function parameters for light limitation. --------------!
      call leaf3_demand_parm_co2('LIGHT')
      !----- Choose the appropriate solver depending on the kind of photosynthesis. -------!
      select case(thispft%photo_pathway)
      case (3)
         call leaf3_solve_iterative(lightlim,success)
      case (4)
         call leaf3_solve_aofixed(lightlim,success)
      end select
      !------------------------------------------------------------------------------------!
      !     In case success was returned as "false", this means that the light-limited     !
      ! case didn't converge (there was no root).  If this is the case we give up and      !
      ! close all stomata, as there was no viable state for stomata to remain opened.      !
      !------------------------------------------------------------------------------------!
      if (.not. success) then
         call leaf3_copy_solution(stclosed,stopen)
         return
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    2. The Rubisco-limited (aka Vm) case.                                           !
      !------------------------------------------------------------------------------------!
      !----- Update the CO2 demand function parameters for Rubisco limitation. ------------!
      call leaf3_demand_parm_co2('RUBISCO')
      !----- Choose the appropriate solver depending on the kind of photosynthesis. -------!
      select case(thispft%photo_pathway)
      case (3)
         call leaf3_solve_iterative(rubiscolim,success)
      case (4)
         call leaf3_solve_aofixed(rubiscolim,success)
      end select
      !------------------------------------------------------------------------------------!
      !     In case success was returned as "false", this means that the Rubisco-limited   !
      ! case didn't converge (there was no root).  If this is the case we give up and      !
      ! close all stomata, as there was no viable state for stomata to remain opened.      !
      !------------------------------------------------------------------------------------!
      if (.not. success) then
         call leaf3_copy_solution(stclosed,stopen)
         return
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      ! 3. The low CO2 concentration case (aka the CO2 case).                              !
      !------------------------------------------------------------------------------------!
      !----- Update the CO2 demand function parameters for CO2 limitation. ----------------!
      call leaf3_demand_parm_co2('CO2')
      !----- Choose the appropriate solver depending on the kind of photosynthesis. -------!
      select case(thispft%photo_pathway)
      case (3)
         !---------------------------------------------------------------------------------!
         !    C3, there is no CO2 limitation in this formulation.  Copy the closed stomata !
         ! case, but assign a large number for carbon_demand, so this will never be        !
         ! chosen.                                                                         !
         !---------------------------------------------------------------------------------!
         call leaf3_copy_solution(stclosed,co2lim)
         co2lim%co2_demand = discard
         success           = .true.

      case (4)
         !---------------------------------------------------------------------------------!
         !    C4, both carbon_demand and stomatal conductance of water depend on the       !
         ! inter-cellular CO2 concentration.  We must find all these three variables       !
         ! simultaneously, using an iterative method.                                      !
         !---------------------------------------------------------------------------------!
         call leaf3_solve_iterative(co2lim,success)
      end select
      !------------------------------------------------------------------------------------!
      !     In case success was returned as "false", this means that the CO2-limited case  !
      ! didn't converge (there was no root).  If this is the case we give up and close all !
      ! stomata, as there was no viable state for stomata to remain opened.                !
      !------------------------------------------------------------------------------------!
      if (.not. success) then
         call leaf3_copy_solution(stclosed,stopen)
         co2lim%co2_demand = discard
         return
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If we have reached this point, it means that we found solutions for all the    !
      ! cases.  The actual solution will be the one with the lowest carbon demand.         !
      !------------------------------------------------------------------------------------!
      if (lightlim%co2_demand <= rubiscolim%co2_demand .and.                               &
          lightlim%co2_demand <= co2lim%co2_demand              )  then
         !----- Light is the strongest limitation. ----------------------------------------!
         call leaf3_copy_solution(lightlim,stopen)
      elseif (rubiscolim%co2_demand <  lightlim%co2_demand .and.                           &
              rubiscolim%co2_demand <= co2lim%co2_demand        ) then
         !----- Rubisco is the strongest limitation. --------------------------------------!
         call leaf3_copy_solution(rubiscolim,stopen)
      else
         !----- CO2 is the strongest limitation. ------------------------------------------!
         call leaf3_copy_solution(co2lim,stopen)
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine leaf3_photo_exact_solver
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine solves the model for the case where the carbon demand doesn't     !
   ! depend on the internal carbon.  This is simpler than the iterative case because we    !
   ! can solve through a quadratic equation for the stomatal conductance for water vapour. !
   !---------------------------------------------------------------------------------------!
   subroutine leaf3_solve_aofixed(answer,success)
      use leaf_coms, only : gsw_2_gsc8        & ! intent(in)
                          , c34smax_gsw8      & ! intent(in)
                          , c34smin_lint_co28 & ! intent(in)
                          , c34smax_lint_co28 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(leaf3_fl_sol) , intent(out) :: answer
      logical            , intent(out) :: success
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                     :: qterm1
      real(kind=8)                     :: qterm2
      real(kind=8)                     :: qterm3
      real(kind=8)                     :: aquad
      real(kind=8)                     :: bquad
      real(kind=8)                     :: cquad
      real(kind=8)                     :: discr
      real(kind=8)                     :: restot
      real(kind=8)                     :: gswroot1
      real(kind=8)                     :: gswroot2
      real(kind=8)                     :: ciroot1
      real(kind=8)                     :: ciroot2
      logical                          :: bounded1
      logical                          :: bounded2
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   1. Initialise the success flag as true.  In case we have trouble solving this    !
      !      case, we switch the flag to false before we quit.                             !
      !------------------------------------------------------------------------------------!
      success = .true.
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   2. Since the carbon demand doesn't depend on the intercellular CO2, compute it   !
      !      using the first guess.                                                        !
      !------------------------------------------------------------------------------------!
      answer%co2_demand = leaf3_co2_demand(met%can_co2)
      !------------------------------------------------------------------------------------!



      !----- 3. Compute the leaf surface CO2. ---------------------------------------------!
      answer%lsfc_co2  = met%can_co2 - answer%co2_demand / met%blyr_cond_co2
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   3. Here we check the sign of the carbon demand.                                  !
      !------------------------------------------------------------------------------------!
      if (answer%co2_demand <= 0.d0) then
         !---------------------------------------------------------------------------------!
         !     If carbon demand is zero or negative, this means that light is below the    !
         ! light compensation point, so all stomata should remain closed.                  !
         !---------------------------------------------------------------------------------!
         answer%stom_cond_h2o = thispft%b
         answer%stom_cond_co2 = gsw_2_gsc8 * answer%stom_cond_h2o
         answer%lint_co2      = answer%lsfc_co2 - answer%co2_demand / answer%stom_cond_co2
      else
         !---------------------------------------------------------------------------------!
         !     Carbon demand is positive, look for a solution.                             !
         !---------------------------------------------------------------------------------!
         !----- Find auxiliary coefficients to compute the quadratic terms. ---------------!
         qterm1 = (met%can_co2 - aparms%compp) * met%blyr_cond_co2 - answer%co2_demand
         qterm2 = (thispft%d0 + met%lint_shv - met%can_shv) * met%blyr_cond_h2o
         qterm3 = thispft%m * answer%co2_demand * thispft%d0 * met%blyr_cond_co2
         !----- Find the coefficients for the quadratic equation. -------------------------!
         aquad = qterm1 * thispft%d0
         bquad = qterm1 * qterm2 - aquad * thispft%b - qterm3
         cquad = - qterm1 * qterm2 * thispft%b - qterm3 * met%blyr_cond_h2o
         !----- Solve the quadratic equation for gsw. -------------------------------------!
         if (aquad == 0.d0) then
            !----- Not really a quadratic equation. ---------------------------------------!
            gswroot1 = -cquad / bquad
            gswroot2 = discard
         else
            !----- A quadratic equation, find the discriminant. ---------------------------!
            discr = bquad * bquad - 4.d0 * aquad * cquad
            !----- Decide what to do based on the discriminant. ---------------------------!
            if (discr == 0.d0) then
               !----- Double root. --------------------------------------------------------!
               gswroot1 = - bquad / (2.d0 * aquad)
               gswroot2 = discard
            elseif (discr > 0.d0) then
               !----- Two distinct roots. -------------------------------------------------!
               gswroot1 = (- bquad - sqrt(discr)) / (2.d0 * aquad)
               gswroot2 = (- bquad + sqrt(discr)) / (2.d0 * aquad)
            else
               !----- None of the roots are real, this solution failed. -------------------!
               success = .false.
               return
            end if
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Check both solutions, and decide which one makes sense.  Once the right     !
         ! one is determined, compute the stomatal resistance for CO2 and the inter-       !
         ! cellular CO2 concentration.  In case both make solutions make sense (unlikely), !
         ! we decide the root based on the intercellular CO2.                              !
         !---------------------------------------------------------------------------------!
         bounded1 = gswroot1 >= thispft%b .and. gswroot1 <= c34smax_gsw8
         bounded2 = gswroot2 >= thispft%b .and. gswroot2 <= c34smax_gsw8
         if (bounded1 .and. bounded2) then
            !----- Both solutions are valid, warn the user as this should never happen. ---!
            ciroot1 = answer%lsfc_co2 - answer%co2_demand / (gsw_2_gsc8 * gswroot1)
            ciroot2 = answer%lsfc_co2 - answer%co2_demand / (gsw_2_gsc8 * gswroot2)

            bounded1 = ciroot1 >= c34smin_lint_co28 .and.                                  &
                       ciroot1 <= min(c34smax_lint_co28,met%can_co2)
            bounded2 = ciroot2 >= c34smin_lint_co28 .and.                                  &
                       ciroot2 <= min(c34smax_lint_co28,met%can_co2)
             
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
               success = .false.
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
            success = .false.
            return
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   8. Lastly, find the surface water specific humidity.                             !
      !------------------------------------------------------------------------------------!
      answer%lsfc_shv = ( answer%stom_cond_h2o * met%lint_shv                              &
                        + met%blyr_cond_h2o    * met%can_shv  )                            &
                      / ( met%blyr_cond_h2o + answer%stom_cond_h2o)
      !------------------------------------------------------------------------------------!

      return
   end subroutine leaf3_solve_aofixed
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
   subroutine leaf3_solve_iterative(answer,converged)
      use leaf_coms, only : gsw_2_gsc8        ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(leaf3_fl_sol) , intent(out) :: answer    ! The strutcure with the answer
      logical            , intent(out) :: converged ! A solution was found     [      T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                     :: ci        ! Intercellular CO2         [  mol/mol]
      real(kind=8)                     :: cia       ! Smallest/previous guess   [  mol/mol]
      real(kind=8)                     :: ciz       ! Largest/new guess         [  mol/mol]
      real(kind=8)                     :: deriv     ! Function derivative       [  mol/mol]
      real(kind=8)                     :: fun       ! Function evaluation       [mol²/mol²]
      real(kind=8)                     :: funa      ! Smallest  guess function  [mol²/mol²]
      real(kind=8)                     :: funz      ! Largest   guess function  [mol²/mol²]
      real(kind=8)                     :: delta     ! Aux. var to find 2nd guess[         ]
      real(kind=8)                     :: cimin     ! Minimum intercell. CO2    [  mol/mol]
      real(kind=8)                     :: cimax     ! Maximum intercell. CO2    [  mol/mol]
      integer                          :: itn       ! Iteration counter         [      ---]
      integer                          :: itb       ! Iteration counter         [      ---]
      logical                          :: zside     ! Check for 1-sided appr.   [      ---]
      logical                          :: hitmin    ! 2nd guess tried minimum   [      ---]
      logical                          :: hitmax    ! 2nd guess tried maximum   [      ---]
      logical                          :: bounded   ! Guess range is bounded.   [      T|F]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Initialise the convergence flag.  Here we start realistic, ops, I mean,       !
      ! pessimistic, and assume that we are failing.  This will be switched to true only   !
      ! if a real answer is found, so in case we quit the sub-routine due to impossible    !
      ! solution, the result will be failure.                                              !
      !------------------------------------------------------------------------------------!
      converged = .false.
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Determine the minimum and maximum intercellular CO2 that will still produce  a  !
      ! positive and bounded conductance, which should be above the cuticular conductance  !
      ! (otherwise there is no reason to keep stomata opened.                              !
      !------------------------------------------------------------------------------------!
      call leaf3_lint_co2_bounds(cimin,cimax,bounded)
      if (bounded) then
         !---------------------------------------------------------------------------------!
         !     We have found bounds, find put the first guess in the middle and find the   !
         ! function evaluation and derivative of this guess.                               !
         !---------------------------------------------------------------------------------!
         cia = sqrt(cimin*cimax)
         call leaf3_iter_step(.true.,cia,funa,deriv)
      else
         !----- No reasonable bound was found, we don't even bother solving this case. ----!
         return
      end if
      !------------------------------------------------------------------------------------!



      !----- Copy to the new guess just in case it fails at the first iteration -----------!
      ciz = cia
      fun = funa
      !------------------------------------------------------------------------------------!

      if (fun == 0.d0) then 
         !----- We have actually hit the jackpot, the answer is ciz. ----------------------!
         ci        = ciz
         converged = .true.
      end if


      !------------------------------------------------------------------------------------!
      !     Enter Newton's method loop in case we haven't found the answer.                !
      !------------------------------------------------------------------------------------!
      if (.not. converged) then
         newloop: do itn = 1,maxfpofl/6
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
               exit newloop
            elseif (converged) then
               !----- Converged, find the root as the mid-point. --------------------------!
               ci  = 5.d-1 * (cia+ciz)
               exit newloop
            else
               !----- Not there yet, update the function evaluation and the derivative. ---!
               call leaf3_iter_step(.true.,ciz,fun,deriv)
               if (fun == 0.d0) then 
                  !----- We have actually hit the jackpot, the answer is ciz. -------------!
                  ci        = ciz
                  converged = .true.
                  exit newloop
               end if
            end if
            !------------------------------------------------------------------------------!
         end do newloop
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
            call leaf3_iter_step(.false.,cia,funa,deriv)
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
               call leaf3_iter_step(.false.,ciz,funz,deriv)
               zside = funa*funz < 0.d0
               if (zside) exit zgssloop
            end do zgssloop
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     If we have reached this points, it means that there is a root that hasn't   !
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
            call leaf3_iter_step(.false.,ci,fun,deriv)


            !------ Define the new interval based on the intermediate value theorem. ------!
            if (fun == 0.d0) then
               !----- We have actually hit the jackpot, the answer is ciz. ----------------!
               converged = .true.
               exit fpoloop
            elseif (fun*funa < 0.d0 ) then
               ciz = ci
               funz  = fun
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa=funa * 5.d-1
               !----- We just updated zside, set zside to true. ---------------------------!
               zside = .true.
            else
               cia = ci
               funa   = fun
               !----- If we are updating aside again, modify aside (Illinois method) ------!
               if (.not. zside) funz=funz * 5.d-1
               !----- We just updated aside, set zside to false. --------------------------!
               zside = .false.
            end if
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
         answer%co2_demand    = leaf3_co2_demand(answer%lint_co2)
         !----- 4. Compute the actual stomatal conductance of water and CO2. --------------!
         answer%stom_cond_h2o = leaf3_stom_cond_h2o(answer%lint_co2,answer%co2_demand)
         answer%stom_cond_co2 = gsw_2_gsc8 * answer%stom_cond_h2o
         !----- 5. Compute the leaf surface CO2. ------------------------------------------!
         answer%lsfc_co2      = met%can_co2 - answer%co2_demand / met%blyr_cond_co2
         !----- 6. Compute the leaf surface vapour mixing ratio. --------------------------!
         answer%lsfc_shv      = ( answer%stom_cond_h2o * met%lint_shv                      &
                                + met%blyr_cond_h2o    * met%can_shv  )                    &
                              / ( met%blyr_cond_h2o + answer%stom_cond_h2o)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine leaf3_solve_iterative
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
   subroutine leaf3_demand_parm_co2(whichlim)
      use leaf_coms, only : klowco28     ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=*), intent(in) :: whichlim      ! A flag telling which case we are
                                                    !   about to solve
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Define the parameters based on this call.                                       !
      !------------------------------------------------------------------------------------!
      select case(trim(whichlim))
      case ('CLOSED')
         !----- Closed stomata case, or night time.  These are the same for C3 and C4. ----!
         aparms%rho   = 0.d0
         aparms%sigma = 0.d0
         aparms%xi    = 0.d0
         aparms%tau   = 1.d0
         aparms%nu    = - aparms%leaf_resp
      case default
         !---------------------------------------------------------------------------------!
         !     Open stomata case, so now we distinguish between C3 and C4 as their         !
         ! functional forms are different.                                                 !
         !---------------------------------------------------------------------------------!
         select case (thispft%photo_pathway)
         case (3)
            !------------------------------------------------------------------------------!
            !     C3 case.  Decide whether this is the light- or Rubisco-limited case.     !
            !------------------------------------------------------------------------------!
            select case (trim(whichlim))
            case ('LIGHT')
               !---- Light-limited case. --------------------------------------------------!
               aparms%rho   =  aparms%alpha * met%par
               aparms%sigma = -aparms%alpha * met%par * aparms%compp
               aparms%xi    = 1.d0
               aparms%tau   = 2.d0 * aparms%compp
               aparms%nu    = -aparms%leaf_resp

            case ('RUBISCO')
               !----- Rubisco-limited rate of photosynthesis case. ------------------------!
               aparms%rho   =  aparms%vm
               aparms%sigma = -aparms%vm * aparms%compp
               aparms%xi    = 1.d0
               aparms%tau   = aparms%kco2 * (1.d0 + met%can_o2 / aparms%ko2)
               aparms%nu    = -aparms%leaf_resp
            !------------------------------------------------------------------------------!

            case ('CO2')
               !----- CO2-limited for low CO2 concentration case. -------------------------!
               aparms%rho   = 0.d0
               aparms%sigma = 5.d-1 * aparms%vm
               aparms%xi    = 0.d0
               aparms%tau   = 1.d0
               aparms%nu    = -aparms%leaf_resp

            end select
         case (4)
            !------------------------------------------------------------------------------!
            !     C4 case.  There are three possibilities, the light-limited, the Rubisco- !
            ! limited, and the CO2-limited cases.                                          !
            !------------------------------------------------------------------------------!
            select case(trim(whichlim))
            case ('LIGHT')
               !----- Light-limited case. -------------------------------------------------!
               aparms%rho   = 0.d0
               aparms%sigma = aparms%alpha * met%par
               aparms%xi    = 0.d0
               aparms%tau   = 1.d0
               aparms%nu    = - aparms%leaf_resp

            case ('RUBISCO')
               !----- Rubisco-limited rate of photosynthesis case. ------------------------!
               aparms%rho   = 0.d0
               aparms%sigma = aparms%vm
               aparms%xi    = 0.d0
               aparms%tau   = 1.d0
               aparms%nu    = -aparms%leaf_resp

            case ('CO2')
               !----- CO2-limited for low CO2 concentration case. -------------------------!
               aparms%rho   = klowco28 * aparms%vm
               aparms%sigma = 0.d0
               aparms%xi    = 0.d0
               aparms%tau   = 1.d0
               aparms%nu    = -aparms%leaf_resp

            end select
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!
      return
   end subroutine leaf3_demand_parm_co2
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
   subroutine leaf3_iter_step(newton,lint_co2,fun,deriv)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      logical     , intent(in)  :: newton              ! Newton's method step  [       T|F]
      real(kind=8), intent(in)  :: lint_co2            ! Intercell. CO2 conc.  [   mol/mol]
      real(kind=8), intent(out) :: fun                 ! Function evaluation   [       ---]
      real(kind=8), intent(out) :: deriv               ! Derivative.           [   mol/mol]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)              :: stom_cond_h2o       ! Water conductance     [  mol/m²/s]
      real(kind=8)              :: co2_demand          ! CO2 demand            [  mol/m²/s]
      real(kind=8)              :: co2_demand_prime    ! Deriv. of CO2 demand  [1/mol/m²/s]
      real(kind=8)              :: stom_cond_h2o_prime ! Deriv. of H2O cond.   [1/mol/m²/s]
      real(kind=8)              :: efun1               ! 1st term
      real(kind=8)              :: efun2               ! 2nd term
      real(kind=8)              :: efun3               ! 3rd term
      real(kind=8)              :: eprime1             ! Derivative of the 1st term
      real(kind=8)              :: eprime2             ! Derivative of the 2nd term
      real(kind=8)              :: eprime3             ! Derivative of the 3rd term
      !------------------------------------------------------------------------------------!


      !----- Find the CO2 demand. ---------------------------------------------------------!
      co2_demand       = leaf3_co2_demand(lint_co2)
      !----- Find the stomatal conductance of water. --------------------------------------!
      stom_cond_h2o    = leaf3_stom_cond_h2o(lint_co2,co2_demand)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the function components, then the function evaluation.                    !
      !------------------------------------------------------------------------------------!
      efun1 = (stom_cond_h2o - thispft%b) / (thispft%m * co2_demand)
      efun2 = (met%can_co2 - aparms%compp - co2_demand/met%blyr_cond_co2)
      efun3 = 1.d0 + ( met%blyr_cond_h2o * (met%lint_shv - met%can_shv)                    &
                     / (thispft%d0 * (met%blyr_cond_h2o + stom_cond_h2o)))
      fun   = efun1 * efun2 * efun3 - 1.d0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case this is a Newton's step, we must also compute the derivatives.  Other- !
      ! wise, we assign any number, but it should not be used.                             !
      !------------------------------------------------------------------------------------!
      if (newton) then
         !----- CO2 demand. ---------------------------------------------------------------!
         co2_demand_prime    = leaf3_co2_demand_prime(lint_co2,co2_demand)
         !----- stomatal conductance of water. --------------------------------------------!
         stom_cond_h2o_prime = leaf3_stom_cond_h2o_prime(lint_co2,stom_cond_h2o,co2_demand &
                                                        ,co2_demand_prime)
         !----- Function components. ------------------------------------------------------!
         eprime1 = ( stom_cond_h2o_prime * co2_demand                                      &
                   - co2_demand_prime * (stom_cond_h2o - thispft%b) )                      &
                 / (thispft%m * co2_demand * co2_demand)

         eprime2 = - co2_demand_prime / met%blyr_cond_co2

         eprime3 = - (efun3 - 1.d0) *  stom_cond_h2o_prime                                 &
                 / ( met%blyr_cond_h2o + stom_cond_h2o )

         deriv   = eprime1 * efun2   * efun3                                               &
                 + efun1   * eprime2 * efun3                                               &
                 + efun1   * efun2   * eprime3
      else
         deriv = 0.d0
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine leaf3_iter_step
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the CO2 demand given the intercellular CO2 concentration.  !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function leaf3_co2_demand(lint_co2)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: lint_co2 ! Intercellular CO2 concentration    [  mol/mol]
      !------------------------------------------------------------------------------------!
      
      leaf3_co2_demand = (aparms%rho * lint_co2 + aparms%sigma)                            &
                       / (aparms%xi  * lint_co2 + aparms%tau  ) + aparms%nu
      return
   end function leaf3_co2_demand
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the derivative of the CO2 regarding the intercellular CO2  !
   ! concentration.                                                                        !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function leaf3_co2_demand_prime(lint_co2,co2_demand)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: lint_co2   ! Intercellular CO2 concentration  [  mol/mol]
      real(kind=8), intent(in) :: co2_demand ! CO2 demand                       [ mol/m²/s]
      !------------------------------------------------------------------------------------!
      
      leaf3_co2_demand_prime = (aparms%rho - aparms%xi * (co2_demand - aparms%nu))         &
                             / (aparms%xi  * lint_co2 + aparms%tau  )
      return
   end function leaf3_co2_demand_prime
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the stomatal conductance of water given the intercellular  !
   ! CO2 concentration and the CO2 demand.                                                 !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function leaf3_stom_cond_h2o(lint_co2,co2_demand)
      use leaf_coms, only : gsw_2_gsc8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: lint_co2   ! Intercellular CO2 concentration  [  mol/mol]
      real(kind=8), intent(in) :: co2_demand ! CO2 demand                       [ mol/m²/s]
      !------------------------------------------------------------------------------------!

      leaf3_stom_cond_h2o = met%blyr_cond_co2 * co2_demand                                 &
                          / (gsw_2_gsc8 * ( (met%can_co2 - lint_co2) * met%blyr_cond_co2   &
                                          - co2_demand ) )

      return
   end function leaf3_stom_cond_h2o
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the stomatal conductance of water derivative given the     !
   ! intercellular CO2 concentration, the water stomatal conductance, and the CO2 demand   !
   ! and its derivative.                                                                   !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function leaf3_stom_cond_h2o_prime(lint_co2,stom_cond_h2o,co2_demand       &
                                                  ,co2_demand_prime)
      use leaf_coms, only : gsw_2_gsc8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: lint_co2         ! Intercell. CO2 conc.      [   mol/mol]
      real(kind=8), intent(in) :: stom_cond_h2o    ! Water                     [  mol/m²/s]
      real(kind=8), intent(in) :: co2_demand       ! CO2 demand                [  mol/m²/s]
      real(kind=8), intent(in) :: co2_demand_prime ! Derivative of CO2 demand  [1/mol/m²/s]
      !------------------------------------------------------------------------------------!

      leaf3_stom_cond_h2o_prime = stom_cond_h2o                                            &
                                * ( co2_demand_prime / co2_demand                           &
                                 + (met%blyr_cond_co2 + co2_demand_prime)                  &
                                 / (gsw_2_gsc8 * ( (met%can_co2 - lint_co2)                &
                                                 * met%blyr_cond_co2 - co2_demand) ))

      return
   end function leaf3_stom_cond_h2o_prime
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
   subroutine leaf3_lint_co2_bounds(cimin,cimax,bounded)
      use leaf_coms, only : gsw_2_gsc8        & ! intent(in)
                          , gbw_2_gbc8        & ! intent(in)
                          , c34smin_lint_co28 & ! intent(in)
                          , c34smax_lint_co28 & ! intent(in)
                          , c34smax_gsw8      ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(out)  :: cimin   ! Minimum intercellular CO2         [  mol/mol]
      real(kind=8), intent(out)  :: cimax   ! Maximum intercellular CO2         [  mol/mol]
      logical     , intent(out)  :: bounded ! This problem is bounded           [      T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)               :: gsw     ! The stom. conductance for water   [ mol/m²/s]
      real(kind=8)               :: restot  ! Total resistance (bnd.lyr.+stom.) [ m² s/mol]
      real(kind=8)               :: aquad   ! Quadratic coefficient             [      ---]
      real(kind=8)               :: bquad   ! Linear coefficient                [  mol/mol]
      real(kind=8)               :: cquad   ! Intercept                         [mol²/mol²]
      real(kind=8)               :: ciAo    ! Ci at singularity where Aopen=0   [  mol/mol]
      real(kind=8)               :: cigsw   ! Ci at sing. where gbc(ca-ci)=Ao   [  mol/mol]
      real(kind=8)               :: ciQ     ! Ci at singularity where qi-qs=-D0 [  mol/mol]
      real(kind=8)               :: xtmp    ! variable for ciQ equation
      real(kind=8)               :: ytmp    ! variable for ciQ
      real(kind=8)               :: ztmp    ! variable for ciQ
      real(kind=8)               :: wtmp    ! variable for ciQ
      real(kind=8)               :: discr   ! The discriminant of the quad. eq. [mol²/mol²]
      real(kind=8)               :: ciroot1 ! 1st root for the quadratic eqn.   [  mol/mol]
      real(kind=8)               :: ciroot2 ! 2nd root for the quadratic eqn.   [  mol/mol]
      real(kind=8), dimension(2) :: cibnds  ! Good root for the low gsw case    [  mol/mol]
      logical                    :: ok1     ! 1st root is okay.                 [      T|F]
      logical                    :: ok2     ! 2nd root is okay.                 [      T|F]
      integer                    :: ibnd    ! Loop for low and high conductance.
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! First case: This check will find when Aopen goes to 0., which causes a singularity !
      ! in the function of which we are looking for a root.                                !
      !------------------------------------------------------------------------------------!
      ciAo = - (aparms%tau * aparms%nu + aparms%sigma)                                     &
           / (aparms%xi  * aparms%nu + aparms%rho  )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! Second case: This check finds which ci causes the terms [gbc (ca -ci) - Aopen] to  !
      ! be 0.   This will cause a singularity in the gsw function.                         !
      !------------------------------------------------------------------------------------!
      !----- 1. Define the coefficients for the quadratic equation. -----------------------!
      aquad = met%blyr_cond_co2 * aparms%xi
      bquad = aparms%xi * (aparms%nu - met%blyr_cond_co2 * met%can_co2 )                   &
            + met%blyr_cond_co2 * aparms%tau + aparms%rho
      cquad = aparms%tau * (aparms%nu - met%blyr_cond_co2 * met%can_co2 ) + aparms%sigma
      !----- 2. Decide whether this is a true quadratic case or not. ----------------------!
      if (aquad /= 0.d0) then
         !---------------------------------------------------------------------------------!
         !     This is a true quadratic case, the first step is to find the discriminant.  !
         !---------------------------------------------------------------------------------!
         discr = bquad * bquad - 4.d0 * aquad * cquad
         if (discr == 0.d0) then
            !------------------------------------------------------------------------------!
            !      Discriminant is zero, both roots are the same.  We save only one, and   !
            ! make the other negative, which will make the guess discarded.                !
            !------------------------------------------------------------------------------!
            ciroot1      = - bquad / (2.d0 * aquad)
            ciroot2      = - discard
         elseif (discr > 0.d0) then
            ciroot1 = (- bquad + sqrt(discr)) / (2.d0 * aquad)
            ciroot2 = (- bquad - sqrt(discr)) / (2.d0 * aquad)
         else
            !----- Discriminant is negative.  Impossible to solve. ------------------------!
            ciroot1      = - discard
            ciroot2      = - discard
         end if
      else
         !---------------------------------------------------------------------------------!
         !    This is a linear case, the xi term is zero.  There is only one number that   !
         ! works for this case.                                                            !
         !---------------------------------------------------------------------------------!
         ciroot1      = - cquad / bquad
         ciroot2      = - discard
         !----- Not used, just for the debugging process. ---------------------------------!
         discr        = bquad * bquad
      end if
      !------------------------------------------------------------------------------------!
      !     Save the largest of the values.  In case both were discarded, we switch it to  !
      ! the positive discard so this will never be chosen.                                 !
      !------------------------------------------------------------------------------------!
      cigsw=max(ciroot1, ciroot2)
      if (cigsw == -discard) cigsw = discard
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      ! Third case: This will find the intercellular CO2 mixing ratio that causes the      !
      ! qi-qs term to be equal to -D0, which also creates a singularity in the stomatal    !
      ! conductance.                                                                       !
      !------------------------------------------------------------------------------------!
      !----- 1. Find some auxiliary variables. --------------------------------------------!
      xtmp = met%blyr_cond_h2o * ( met%can_shv - met%lint_shv - thispft%d0)/thispft%d0
      ytmp = met%blyr_cond_co2 + xtmp * gbw_2_gbc8
      ztmp = xtmp * gbw_2_gbc8 * met%blyr_cond_co2
      wtmp = ztmp * met%can_co2 - ytmp * aparms%nu

      !----- 2. Define the coefficients for the quadratic equation. -----------------------!
      aquad = ztmp * aparms%xi
      bquad = ytmp * aparms%rho - aparms%xi * wtmp + ztmp * aparms%tau
      cquad = - aparms%tau * wtmp + ytmp * aparms%sigma

      !----- 3. Decide whether this is a true quadratic case or not. ----------------------!
      if (aquad /= 0.d0) then
         !---------------------------------------------------------------------------------!
         !     This is a true quadratic case, the first step is to find the discriminant.  !
         !---------------------------------------------------------------------------------!
         discr = bquad * bquad - 4.d0 * aquad * cquad
         if (discr == 0.d0) then
            !------------------------------------------------------------------------------!
            !      Discriminant is zero, both roots are the same.  We save only one, and   !
            ! make the other negative, which will make the guess discarded.                !
            !------------------------------------------------------------------------------!
            ciroot1 = - bquad / (2.d0 * aquad)
            ciroot2 = -discard
         elseif (discr > 0.d0) then
            ciroot1 = (- bquad + sqrt(discr)) / (2.d0 * aquad)
            ciroot2 = (- bquad - sqrt(discr)) / (2.d0 * aquad)
         else
            !----- Discriminant is negative.  Impossible to solve. ------------------------!
            ciroot1      = -discard
            ciroot2      = -discard
         end if
      else
         !---------------------------------------------------------------------------------!
         !    This is a linear case, the xi term is zero.  There is only one number        !
         ! that works for this case.                                                       !
         !---------------------------------------------------------------------------------!
         ciroot1 = - cquad / bquad
         ciroot2 = -discard
         !----- Not used, just for the debugging process. ---------------------------------!
         discr = bquad * bquad
      end if
      !------------------------------------------------------------------------------------!
      !     Save the largest of the values.  In case both were discarded, we switch it to  !
      ! the positive discard so this will never be chosen.                                 !
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
      ciroot2 = min(ciQ, cigsw,met%can_co2)
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

      return
   end subroutine leaf3_lint_co2_bounds
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This is the temperature-dependence for the physiology parameters, represented as !
   ! in function, written as a in Collatz et al. (1991).                                   !
   ! The output variable will have the same unit as the reference value at tphysref.       !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function leaf3_qten_fun(temp,refval,q10)
      use leaf_coms, only : temp0_q108  & ! intent(in)
                          , slope_q108  ! ! intent(in)
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
      leaf3_qten_fun = refval * q10 ** (slope_q108 * (temp - temp0_q108))
      !------------------------------------------------------------------------------------!

      return
   end function leaf3_qten_fun
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Find the minimum amount of radiation for which we will consider daytime.  This    !
   ! cannot be a constant because it depends on the PFT and the environmental conditions.  !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function leaf3_twilight_min()
      use leaf_coms, only : gsw_2_gsc8
      implicit none
      
      leaf3_twilight_min = ( aparms%leaf_resp  * (met%can_co2 + 2.d0 * aparms%compp) )     &
                         / ( aparms%alpha      * (met%can_co2 -        aparms%compp) )
                        

      return
   end function leaf3_twilight_min
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine copies a solution structure from one variable to the other.       !
   !---------------------------------------------------------------------------------------!
   subroutine leaf3_copy_solution(source_sol,target_sol)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(leaf3_fl_sol), intent(in)  :: source_sol ! The solution to be copied.
      type(leaf3_fl_sol), intent(out) :: target_sol ! The solution to be filled.
      !------------------------------------------------------------------------------------!


      !----- Copy all elements from one structure to the other.
      target_sol%lsfc_shv      = source_sol%lsfc_shv
      target_sol%lsfc_co2      = source_sol%lsfc_co2
      target_sol%lint_co2      = source_sol%lint_co2
      target_sol%co2_demand    = source_sol%co2_demand
      target_sol%stom_cond_h2o = source_sol%stom_cond_h2o
      target_sol%stom_cond_co2 = source_sol%stom_cond_co2
      !------------------------------------------------------------------------------------!


      return
   end subroutine leaf3_copy_solution
   !=======================================================================================!
   !=======================================================================================!
end module leaf3_physiol
!==========================================================================================!
!==========================================================================================!

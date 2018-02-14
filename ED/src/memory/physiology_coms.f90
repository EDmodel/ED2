!==========================================================================================!
!==========================================================================================!
!     This module contains several parameters used in the Farquar Leuning photosynthesis   !
! solver.  The following references are used as a departing point:                         !
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
module physiology_coms
   use ed_max_dims, only : str_len ! ! intent(in)

   implicit none


   !=======================================================================================!
   !=======================================================================================!
   !     Variables that are defined by the user in the namelist.                           !
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! IPHYSIOL --  This variable will determine the functional form that will control how   !
   !              the various parameters will vary with temperature, and how the CO2       !
   !              compensation point for gross photosynthesis (Gamma*) will be found.      !
   !              Options are:                                                             !
   !                                                                                       !
   ! 0 -- Original ED-2.1, we use the "Arrhenius" function as in Foley et al. (1996) and   !
   !      Moorcroft et al. (2001).  Gamma* is found using the parameters for tau as in     !
   !      Foley et al. (1996).                                                             !
   ! 1 -- Modified ED-2.1.  In this case Gamma* is found using the Michaelis-Mentel        !
   !      coefficients for CO2 and O2, as in Farquhar et al. (1980) and in CLM.            !
   ! 2 -- Collatz et al. (1991).  We use the power (Q10) equations, with Collatz et al.    !
   !      parameters for compensation point, and the Michaelis-Mentel coefficients.  The   !
   !      correction for high and low temperatures are the same as in Moorcroft et al.     !
   !      (2001).                                                                          !
   ! 3 -- Same as 2, except that we find Gamma* as in Farquhar et al. (1980) and in CLM.   !
   !---------------------------------------------------------------------------------------!
   integer                :: iphysiol 
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! N_PLANT_LIM -- This controls whether plant photosynthesis can be limited by nitrogen. !
   !                0.  No limitation                                                      !
   !                1.  ED-2.1 nitrogen limitation model.                                  !
   !---------------------------------------------------------------------------------------!
   integer                :: n_plant_lim
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! H2O_PLANT_LIM -- this determines whether plant photosynthesis can be limited by       !
   !                  soil moisture, the FSW, defined as FSW = Supply / (Demand + Supply). !
   !                                                                                       !
   !                  Demand is always the transpiration rates in case soil moisture is    !
   !                  not limiting (the psi_0 term times LAI).  The supply is determined   !
   !                  by Kw * nplant * Broot * Available_Water, and the definition of      !
   !                  available water changes depending on H2O_PLANT_LIM:                  !
   !                  0.  Force FSW = 1 (effectively available water is infinity).         !
   !                  1.  Available water is the total soil water above wilting point,     !
   !                      integrated across all layers within the rooting zone.            !
   !                  2.  Available water is the soil water at field capacity minus        !
   !                      wilting point, scaled by the so-called wilting factor:           !
   !                      (psi(k) - (H - z(k)) - psi_wp) / (psi_fc - psi_wp)               !
   !                      where psi is the matric potentital at layer k, z is the layer    !
   !                      depth, H it the crown height and psi_fc and psi_wp are the       !
   !                      matric potentials at wilting point and field capacity.           !
   !---------------------------------------------------------------------------------------!
   integer               :: h2o_plant_lim 
   !---------------------------------------------------------------------------------------!



   integer :: plant_hydro_scheme
   !---------------------------------------------------------------------------------------!
   !< PLANT_HYDRO_SCHEME -- Whether to track plant hydrodynamics.\n
   !<                       0. No hydraulics (leaf and wood are always saturated)\n
   !<                       1. Track plant hydrodynamics using parameters from
   !<  Christofferson et al. 2016 GMD\n
   !<                       2. Track plant hydrodynamics using parameters from
   !<  Xu et al. 2016 New Phytologist\n
   !<                      -1. Same as 1 but leaf/wood heat capacity does not
   !< change with internal water content.
   !<                      -2. Same as 2 but leaf/wood heat capacity does not
   !< change with internal water content.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !  IDDMORT_SCHEME -- This flag determines whether storage should be accounted in the    !
   !                    carbon balance.                                                    !
   !                    0 -- Carbon balance is done in terms of fluxes only.  This is the  !
   !                         default in ED-2.1                                             !
   !                    1 -- Carbon balance is offset by the storage pool.  Plants will be !
   !                         in negative carbon balance only when they run out of storage  !
   !                         and are still losing more carbon than gaining.                !
   !  CBR_SCHEME     -- This flag determines which carbon stress scheme is used:           !
   !                    0 -- Single stress. CBR = cb/cb_mlmax                              !
   !                         cb_mlmax is the carbon balance in full sun and no moisture    !
   !                         limitation                                                    !
   !                    1 -- Co-limitation from light and moisture                         !
   !                         CBR_LIGHT = cb/cb_lightmax and CBR_MOIST = cb/cb_moistmax     !
   !                         CBR_LIGHT and CBR_MOIST are then weighted according to        !
   !                         DDMORT_CONST (below)                                          !
   !                    2 -- Leibig Style, i.e. limitation from either light or moisture   !
   !                         depending on which is lower at a given point in time          !
   !                         CBR = cb/max(cb_lightmax, cb_moistmax)                        !
   !                                                                                       !
   !  DDMORT_CONST   -- CBR_Scheme = 1 only                                                !
   !                    This constant (k) determines the relative contribution of light    !
   !                    and soil moisture to the density-dependent mortality rate.  Values !
   !                    range from 0 (soil moisture only) to 1 (light only).               !
   !                                                                                       !
   !                                    mort1                                              !
   !                  mu_DD = -------------------------                                    !
   !                            1 + exp [ mort2 * cr ]                                     !
   !                                                                                       !
   !                                 CB                      CB                            !
   !                  cr    = k ------------- + (1 - k) -------------                      !
   !                             CB_lightmax             CB_watermax                       !
   !---------------------------------------------------------------------------------------!
   integer              :: iddmort_scheme
   real                 :: ddmort_const
   integer              :: cbr_scheme
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     The following variables are factors that control photosynthesis and respiration.  !
   ! Notice that some of them are relative values whereas others are absolute.             !
   !                                                                                       !
   ! VMFACT_C3     -- Factor multiplying the default Vm0 for C3 plants (1.0 = default).    !
   ! VMFACT_C4     -- Factor multiplying the default Vm0 for C4 plants (1.0 = default).    !
   ! MPHOTO_TRC3   -- Stomatal slope (M) for tropical C3 plants                            !
   ! MPHOTO_TEC3   -- Stomatal slope (M) for conifers and temperate C3 plants              !
   ! MPHOTO_C4     -- Stomatal slope (M) for C4 plants.                                    !
   ! BPHOTO_BLC3   -- cuticular conductance for broadleaf C3 plants  [umol/m2/s]           !
   ! BPHOTO_NLC3   -- cuticular conductance for needleleaf C3 plants [umol/m2/s]           !
   ! BPHOTO_C4     -- cuticular conductance for C4 plants            [umol/m2/s]           !
   ! KW_GRASS      -- Water conductance for trees, in m2/yr/kgC_root.  This is used only   !
   !                  when H2O_PLANT_LIM is not 0.                                         !
   ! KW_TREE       -- Water conductance for grasses, in m2/yr/kgC_root.  This is used only !
   !                  when H2O_PLANT_LIM is not 0.                                         !
   ! GAMMA_C3      -- The dark respiration factor (gamma) for C3 plants.  Subtropical      !
   !                  conifers will be scaled by GAMMA_C3 * 0.028 / 0.02                   !
   ! GAMMA_C4      -- The dark respiration factor (gamma) for C4 plants.                   !
   ! D0_GRASS      -- The transpiration control in gsw (D0) for ALL grasses.               !
   ! D0_TREE       -- The transpiration control in gsw (D0) for ALL trees.                 !
   ! ALPHA_C3      -- Quantum yield of ALL C3 plants.  This is only applied when           !
   !                  QUANTUM_EFFICIENCY_T = 0.                                            !
   ! ALPHA_C4      -- Quantum yield of C4 plants.  This is always applied.                 !
   ! KLOWCO2IN     -- The coefficient that controls the PEP carboxylase limited rate of    !
   !                  carboxylation for C4 plants.                                         !
   ! RRFFACT       -- Factor multiplying the root respiration factor for ALL PFTs.         !
   !                  (1.0 = default).                                                     !
   ! GROWTHRESP    -- The actual growth respiration factor (C3/C4 tropical PFTs only).     !
   !                  (1.0 = default).                                                     !
   ! LWIDTH_GRASS  -- Leaf width for grasses, in metres.  This controls the leaf boundary  !
   !                  layer conductance (gbh and gbw).                                     !
   ! LWIDTH_BLTREE -- Leaf width for trees, in metres.  This controls the leaf boundary    !
   !                  layer conductance (gbh and gbw).  This is applied to broadleaf trees !
   !                  only.                                                                !
   ! LWIDTH_NLTREE -- Leaf width for trees, in metres.  This controls the leaf boundary    !
   !                  layer conductance (gbh and gbw).  This is applied to conifer trees   !
   !                  only.                                                                !
   ! Q10_C3        -- Q10 factor for C3 plants (used only if IPHYSIOL is set to 2 or 3).   !
   ! Q10_C4        -- Q10 factor for C4 plants (used only if IPHYSIOL is set to 2 or 3).   !
   !---------------------------------------------------------------------------------------!
   real(kind=4)               :: vmfact_c3
   real(kind=4)               :: vmfact_c4
   real(kind=4)               :: mphoto_trc3
   real(kind=4)               :: mphoto_tec3
   real(kind=4)               :: mphoto_c4
   real(kind=4)               :: bphoto_blc3
   real(kind=4)               :: bphoto_nlc3
   real(kind=4)               :: bphoto_c4
   real(kind=4)               :: kw_grass
   real(kind=4)               :: kw_tree
   real(kind=4)               :: gamma_c3
   real(kind=4)               :: gamma_c4
   real(kind=4)               :: d0_grass
   real(kind=4)               :: d0_tree
   real(kind=4)               :: alpha_c3
   real(kind=4)               :: alpha_c4
   real(kind=4)               :: klowco2in
   real(kind=4)               :: rrffact
   real(kind=4)               :: growthresp
   real(kind=4)               :: lwidth_grass
   real(kind=4)               :: lwidth_bltree
   real(kind=4)               :: lwidth_nltree
   real(kind=4)               :: q10_c3
   real(kind=4)               :: q10_c4
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     This flag determines whether to use temperature dep quantum efficiency            !
   !---------------------------------------------------------------------------------------!
   integer                    :: quantum_efficiency_T
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Parameters for the the physiology solver.  These number are assigned at sub-      !
   ! routine init_physiology_params (ed_params.f90).  The minimum conductance is not       !
   ! defined because it must be the cuticular conductance.                                 !
   !---------------------------------------------------------------------------------------!
   !----- Bounds for the new C3 solver. ---------------------------------------------------!
   real(kind=4) :: c34smin_lint_co2 ! Minimum carbon dioxide concentration      [  mol/mol]
   real(kind=4) :: c34smax_lint_co2 ! Maximum carbon dioxide concentration      [  mol/mol]
   real(kind=4) :: c34smax_gsw      ! Maximum stom. conductance for water vap.  [ mol/m²/s]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This is the minimum threshold for the photosynthetically active radiation, in     !
   ! µmol/m²/s to consider non-night time conditions (day time or twilight).               !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: par_twilight_min ! Minimum non-nocturnal PAR.                [ mol/m²/s]
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Parameters that control debugging output.  These are also assigned in the sub-    !
   ! routine init_physiology_params (ed_params.f90).                                       !
   !---------------------------------------------------------------------------------------!
   !----- I should print detailed debug information. --------------------------------------!
   logical                :: print_photo_debug
   !----- File name prefix for the detailed information in case of debugging. -------------!
   character(len=str_len) :: photo_prefix
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     These are constants obtained in Leuning et al. (1995). to convert different       !
   ! conductivities.                                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: gbh_2_gbw ! heat  to water  - leaf boundary layer
   real(kind=4) :: gbw_2_gbc ! water to carbon - leaf boundary layer
   real(kind=4) :: gsw_2_gsc ! water to carbon - stomata
   real(kind=4) :: gsc_2_gsw ! carbon to water - stomata
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This parameter is from F80, and is the ratio between the turnover number for the  !
   ! oxylase function and the turnover number for the carboxylase function.                !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: kookc
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Coefficients for the various Arrhenius and Collatz function calculations of the   !
   ! terms that control the photosynthesis rates and do not depend on the PFT as long as   !
   ! it is C3.  Such terms are compensation point, and the Michaelis-Mentel coefficients   !
   ! for CO2 and O2.                                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: tphysref      ! Reference temperature
   real(kind=4) :: tphysrefi     ! 1./tarrh
   real(kind=4) :: fcoll         ! Factor for the Collatz temperature dependence

   !----- Compensation point neglecting the respiration (Gamma*). -------------------------!
   real(kind=4) :: compp_refval  ! Reference value                                [mol/mol]
   real(kind=4) :: compp_hor     ! "Activation energy" when using Arrhenius       [      K]
   real(kind=4) :: compp_q10     ! Base value when using Collatz                  [    ---]

   !------ These terms are used to find the Michaelis-Menten coefficient for CO2. ---------!
   real(kind=4) :: kco2_refval   ! Reference value                                [mol/mol]
   real(kind=4) :: kco2_hor      ! "Activation energy" when using Arrhenius       [      K]
   real(kind=4) :: kco2_q10      ! Base value when using Collatz                  [    ---]
   !------ These terms are used to find the Michaelis-Menten coefficient for O2. ----------!
   real(kind=4) :: ko2_refval    ! Reference value                                [mol/mol]
   real(kind=4) :: ko2_hor       ! "Activation energy" when using Arrhenius       [      K]
   real(kind=4) :: ko2_q10       ! Base value when using Collatz                  [    ---]
   !---------------------------------------------------------------------------------------!
   !    The following parameter is the k coefficient in F96 that is used to determine the  !
   ! CO2-limited photosynthesis for C4 grasses.                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: klowco2       ! coefficient for low CO2                        [  1/mol]
   !---------------------------------------------------------------------------------------!


   !----- The following parameter is the prescribed O2 concentration. ---------------------!
   real(kind=4) :: o2_ref        ! Oxygen concentration                           [mol/mol]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The following parameters are used to define the quantum yield for C3 plants as a  !
   ! function of temperature.  The reference for this response curve is:                   ! 
   !                                                                                       !
   ! Ehleringer, J. R., 1978: Implications of quantum yield differences on the             !
   !    distributions of C3 and C4 grasses.  Oecologia, 31, 255-267.                       !
   !                                                                                       !
   ! ehleringer_alpha0c is alpha at 0 degrees Celsius, with we use for temperatures below  !
   ! 0 Celsius.                                                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: qyield0
   real(kind=4) :: qyield1
   real(kind=4) :: qyield2
   real(kind=4) :: ehleringer_alpha0c
   !---------------------------------------------------------------------------------------!


   !----- The double precision version of the previous variables. -------------------------!
   real(kind=8) :: c34smin_lint_co28
   real(kind=8) :: c34smax_lint_co28
   real(kind=8) :: c34smax_gsw8
   real(kind=8) :: gbh_2_gbw8
   real(kind=8) :: gbw_2_gbc8
   real(kind=8) :: gsw_2_gsc8
   real(kind=8) :: gsc_2_gsw8
   real(kind=8) :: kookc8
   real(kind=8) :: tphysref8
   real(kind=8) :: tphysrefi8
   real(kind=8) :: fcoll8
   real(kind=8) :: compp_refval8
   real(kind=8) :: compp_hor8
   real(kind=8) :: compp_q108
   real(kind=8) :: kco2_refval8
   real(kind=8) :: kco2_hor8
   real(kind=8) :: kco2_q108
   real(kind=8) :: ko2_refval8
   real(kind=8) :: ko2_hor8
   real(kind=8) :: ko2_q108
   real(kind=8) :: klowco28
   real(kind=8) :: par_twilight_min8
   real(kind=8) :: o2_ref8
   real(kind=8) :: qyield08
   real(kind=8) :: qyield18
   real(kind=8) :: qyield28
   real(kind=8) :: ehleringer_alpha0c8
   !---------------------------------------------------------------------------------------!

end module physiology_coms
!==========================================================================================!
!==========================================================================================!

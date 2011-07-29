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


   !---------------------------------------------------------------------------------------!
   !     Variables that are defined by the user in the namelist.                           !
   !---------------------------------------------------------------------------------------!
   !----- This flag controls which physiology scheme we should use. -----------------------!
   integer                :: iphysiol !  0 -- Original ED-2.1, with temperature functions
                                      !       described as in F96
                                      !  1 -- Original ED-2.1, but with compensation point
                                      !       as a function of the Michaelis-Mentel 
                                      !       constants for CO2 and O2 (like in F80)
                                      !  2 -- Temperature response is based on C91/C92,
                                      !       with temperature decay for Vm0 and leaf 
                                      !       respiration for both C3 and C4 as in C91
                                      !  3 -- Same as 2, except that the compensation point
                                      !       were found as functions of the Michaelis-
                                      !       Mentel constants for CO2 and O2 (like in F80)
   !----- This flag controls whether to use the exact or small perturbation solution. -----!
   integer                :: istoma_scheme
   !----- This flag controls whether the plants should be limited by nitrogen. ------------!
   integer                :: n_plant_lim
   !---------------------------------------------------------------------------------------!
   !     This parameter will decide whether the fraction of open stomata should be         !
   ! calculated through the original method or the new empirical relation.                 !
   !---------------------------------------------------------------------------------------!
   integer               :: h2o_plant_lim ! 0 -- No plant limitation
                                          ! 1 -- Original ED-2.1 model
                                          ! 2 -- CLM-based function.
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     These parameter will decide scaling factors for stomatal slope, Vm0, Kw, and the  !
   ! dark respiration constant.                                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=4)               :: vmfact
   real(kind=4)               :: mfact
   real(kind=4)               :: kfact
   real(kind=4)               :: gamfact
   real(kind=4)               :: d0fact
   real(kind=4)               :: alphafact
   real(kind=4)               :: lwfact
   real(kind=4)               :: thioff
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

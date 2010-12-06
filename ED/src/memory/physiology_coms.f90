!==========================================================================================!
!==========================================================================================!
!     This module contains several parameters used in the Farquar Leuning photosynthesis   !
! solver.  The following references are used as a departing point:                         !
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
!------------------------------------------------------------------------------------------!
module physiology_coms
   use ed_max_dims, only : str_len ! ! intent(in)

   implicit none

   !----- Max roots is used for the old bracketing method. --------------------------------!
   integer, parameter     :: maxroots  = 5
   !----- xdim is the dimension of the vector we are solving. -----------------------------!
   integer, parameter     :: xdim =2
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Variables that are defined by the user in the namelist.                           !
   !---------------------------------------------------------------------------------------!
   !----- This flag controls whether to use the exact or small perturbation solution. -----!
   integer                :: istoma_scheme
   !----- This flag controls whether the plants should be limited by nitrogen. ------------!
   integer                :: n_plant_lim
   !---------------------------------------------------------------------------------------!
   !     This parameter will decide whether the fraction of open stomata should be         !
   ! calculated through the original method or the new empirical relation.                 !
   !---------------------------------------------------------------------------------------!
   integer               :: h2o_plant_lim
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
   !     Coefficients for the various Arrhenius function calculations and the maximum      !
   ! capacity of Rubisco to perform the carboxylase function.                              !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: tarrh         ! Arrhenius' reference temperature
   real(kind=4) :: tarrhi        ! 1./tarrh
   real(kind=4) :: compp_refkin  ! Reference ratio of "kinetic parameters describing the
                                 !    partioning of enzyme activity to carboxylase or
                                 !    oxylase function" (F96)  
   real(kind=4) :: compp_ecoeff  ! Reference exponential coefficient for compensation point
   real(kind=4) :: vm_ecoeff     ! Reference exp. coefficient for carboxylase function.
   real(kind=4) :: vm_tempcoeff  ! Coefficient to control the temperature dependence
   !------ These terms are used to find the Michaelis-Menten coefficient for CO2. ---------!
   real(kind=4) :: kco2_prefac   ! Reference CO2 concentration                    [mol/mol]
   real(kind=4) :: kco2_ecoeff   ! Reference exponential coefficient              [      K]
   !------ These terms are used to find the Michaelis-Menten coefficient for O2. ----------!
   real(kind=4) :: ko2_prefac    ! Reference O2 concentration                     [mol/mol]
   real(kind=4) :: ko2_ecoeff    ! Reference exponential coefficient              [      K]
   !---------------------------------------------------------------------------------------!
   !    The following parameter is the k coefficient in F96 that is used to determine the  !
   ! CO2-limited photosynthesis for C4 grasses.                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: klowco2       ! coefficient for low CO2                        [  1/mol]
   !---------------------------------------------------------------------------------------!


   !----- The following parameter is the prescribed O2 concentration. ---------------------!
   real(kind=4) :: o2_ref        ! Oxygen concentration                           [mol/mol]


   !----- The double precision version of the previous variables. -------------------------!
   real(kind=8) :: c34smin_lint_co28
   real(kind=8) :: c34smax_lint_co28
   real(kind=8) :: c34smax_gsw8
   real(kind=8) :: gbh_2_gbw8 
   real(kind=8) :: gbw_2_gbc8 
   real(kind=8) :: gsw_2_gsc8 
   real(kind=8) :: gsc_2_gsw8 
   real(kind=8) :: tarrh8
   real(kind=8) :: tarrhi8
   real(kind=8) :: compp_refkin8
   real(kind=8) :: compp_ecoeff8
   real(kind=8) :: vm_ecoeff8
   real(kind=8) :: vm_tempcoeff8
   real(kind=8) :: kco2_prefac8
   real(kind=8) :: kco2_ecoeff8
   real(kind=8) :: ko2_prefac8
   real(kind=8) :: ko2_ecoeff8
   real(kind=8) :: klowco28
   real(kind=8) :: par_twilight_min8
   real(kind=8) :: o2_ref8
   !---------------------------------------------------------------------------------------!

end module physiology_coms
!==========================================================================================!
!==========================================================================================!

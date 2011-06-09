!==========================================================================================!
!==========================================================================================!
!    This module contains some structures used by the model to solve the photosynthesis.   !
! The references for this model are given by:                                              !
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
!------------------------------------------------------------------------------------------!
module c34constants
   implicit none

   !---------------------------------------------------------------------------------------!
   !     This structure contains some PFT-dependent physiological and empirical constants  !
   ! for the photosynthesis model.  The units inside the photosynthesis model are not the  !
   ! same as the input data.  Everything that is in µmol becomes mol here, and everything  !
   ! that is in °C becomes Kelvin here.                                                    !
   !---------------------------------------------------------------------------------------!
   type farq_consts
      integer      :: photo_pathway ! The photosynthetic pathway                 [     ---]
      real(kind=8) :: D0            ! Threshold to close stomata when it's dry   [ mol/mol]
      real(kind=8) :: m             ! Empirical const. for conductance calc.     [     ---]
      real(kind=8) :: b             ! Cuticular conductance                      [mol/m²/s]
      real(kind=8) :: alpha0        ! Quantum efficiency                         [     ---]
      real(kind=8) :: vm0           ! Ref. value for maximum rubisco activity    [mol/m²/s]
      real(kind=8) :: vm_low_temp   ! Low temp. threshold for rapid decline      [       K]
      real(kind=8) :: vm_high_temp  ! High temp. threshold for rapid decline     [       K]
      real(kind=8) :: vm_hor        ! "Activation energy" for Vm                 [       K]
      real(kind=8) :: vm_q10        ! Vm base for C91 physiology                 [      --]
      real(kind=8) :: vm_decay_a    ! Decay factor for Vm correction as in C91   [   J/mol]
      real(kind=8) :: vm_decay_b    ! Decay factor for Vm correction as in C91   [ J/mol/K]
      real(kind=8) :: vm_decay_e    ! Decay factor for Vm correction as in F96   [     1/K]
      real(kind=8) :: rd0           ! Term used for respiration                  [mol/m²/s]
      real(kind=8) :: rd_low_temp   ! Low temp. threshold for rapid decline      [       K]
      real(kind=8) :: rd_high_temp  ! High temp. threshold for rapid decline     [       K]
      real(kind=8) :: rd_hor        ! "Activation energy" for Vm                 [       K]
      real(kind=8) :: rd_q10        ! Rd base for C91 physiology                 [      --]
      real(kind=8) :: rd_decay_e    ! Decay factor for Rd correction as in F96   [     1/K]
   end type farq_consts
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This structure contains the canopy air (plus some leaf) conditions.  The short    !
   ! name given in David Medvigy's thesis is given in parenthesis, and the units used in   !
   ! the solver are given in brackets.                                                     !
   !---------------------------------------------------------------------------------------!
   type metinp_vars
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
   end type metinp_vars
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     This structure contains the variables used to define the optimal leaf-level       !
   ! carbon dioxide demand of photosynthesis.  The names follow the original notation      !
   ! found in M06 equations 2.62, 2.63, 2.67, and B.2.                                     !
   !---------------------------------------------------------------------------------------!
   type carb_demand_vars
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
   end type carb_demand_vars
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     This structure contains the variables that will be found for both open and closed !
   ! stomata cases, and in case of the open, for the Rubisco- and light-limited cases,     !
   ! plus the CO2 case for C4 photosynthesis.                                              !
   !---------------------------------------------------------------------------------------!
   type solution_vars
      real(kind=8) :: lsfc_shv      ! Leaf surface specific humidity     (es )  [  mol/mol]
      real(kind=8) :: lsfc_co2      ! Leaf surface CO2 mixing ratio      (Cs )  [  mol/mol]
      real(kind=8) :: lint_co2      ! Leaf internal CO2 mixing ratio     (Ci )  [  mol/mol]
      real(kind=8) :: co2_demand    ! Carbon demand                      (A  )  [ mol/m²/s]
      real(kind=8) :: stom_cond_h2o ! Stomatal water conductance         (gsw)  [ mol/m²/s]
      real(kind=8) :: stom_cond_co2 ! Stomatal CO2 conductance           (gsc)  [ mol/m²/s]
   end type solution_vars
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    These are the structures used throughout the solution.  They are kept here so they !
   ! can be accessed by any subroutine.                                                    !
   !---------------------------------------------------------------------------------------!
   type(farq_consts  )    :: thispft
   type(metinp_vars  )    :: met
   type(carb_demand_vars) :: aparms
   type(solution_vars)    :: stopen
   type(solution_vars)    :: stclosed
   type(solution_vars)    :: rubiscolim
   type(solution_vars)    :: co2lim
   type(solution_vars)    :: lightlim
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   This structure contains the variables that may be copied to the model standard      !
   ! output.                                                                               !
   !---------------------------------------------------------------------------------------!
   type stoma_data
      integer         :: recalc=1   !THIS SHOULD BE INIT IN ED_PARAMS
      real(kind=4)    :: T_L
      real(kind=4)    :: e_A
      real(kind=4)    :: PAR
      real(kind=4)    :: rb_factor
      real(kind=4)    :: prss
      real(kind=4)    :: phenology_factor
      real(kind=4)    :: gsw_open
      integer         :: ilimit
      
      real(kind=4)    :: T_L_residual
      real(kind=4)    :: e_a_residual
      real(kind=4)    :: par_residual
      real(kind=4)    :: rb_residual
      real(kind=4)    :: prss_residual
      real(kind=4)    :: leaf_residual
      real(kind=4)    :: gsw_residual
   end type stoma_data
   !---------------------------------------------------------------------------------------!


   !------ The number of stomatal attributes. ---------------------------------------------!
   integer, parameter :: n_stoma_atts = 16
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine copies a solution structure from one variable to the other.       !
   !---------------------------------------------------------------------------------------!
   subroutine copy_solution(source_sol,target_sol)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(solution_vars), intent(in)  :: source_sol ! The solution to be copied.
      type(solution_vars), intent(out) :: target_sol ! The solution to be filled.
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
   end subroutine copy_solution
   !=======================================================================================!
   !=======================================================================================!
end module c34constants
!==========================================================================================!
!==========================================================================================!

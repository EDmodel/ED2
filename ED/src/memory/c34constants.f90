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
! - L95 - Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf         !
!         nitrogen, photosynthesis, conductance, and transpiration: scaling from leaves to !
!         canopies. Plant, Cell, and Environ., 18, 1183-1200.                              !
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
      real(kind=8) :: alpha         ! Quantum efficiency                         [     ---]
      real(kind=8) :: gamma         ! Proportionality constant for leaf resp.    [     ---]
      real(kind=8) :: vm0_photo     ! Ref. value for maximum rubisco activity    [mol/m²/s]
      real(kind=8) :: vm0_resp      ! Term used for respiration                  [mol/m²/s]
      real(kind=8) :: vm_low_temp   ! Low temp. threshold for rapid decline      [       K]
      real(kind=8) :: vm_high_temp  ! High temp. threshold for rapid decline     [       K]
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
      real(kind=8) :: compp         ! Compens. pt. for gross photos.    (GAMMA) [  mol/mol]
   end type metinp_vars
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     This structure contains the variables used to define the optimal leaf-level       !
   ! carbon dioxide demand of photosynthesis.  The names follow the original notation      !
   ! found in M06 equations 2.62, 2.63, 2.67, and B.2.                                     !
   !---------------------------------------------------------------------------------------!
   type carb_demand_vars
      real(kind=8) :: sigma     ! sigma from equation B.2
      real(kind=8) :: rho       ! rho from equation B.2
      real(kind=8) :: vm        ! Max. capacity of Rubisco to perform the carboxylase func.
      real(kind=8) :: leaf_resp ! Leaf respiration
      real(kind=8) :: xi        ! coefficient for Ci in the denominator (always 1 in M09).
      real(kind=8) :: tau       ! tau from equation B.2
      real(kind=8) :: nu        ! The negative of leaf respiration
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

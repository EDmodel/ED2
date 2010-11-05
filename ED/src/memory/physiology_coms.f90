!==========================================================================================!
!==========================================================================================!
!     This module contains several variables used in the Farquar Leuning photosynthesis    !
! solver.                                                                                  !
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


   !---------------------------------------------------------------------------------------!
   !     This flag controls whether the internal carbon and conductivity should be solved  !
   ! using the new method (2-dimensional Newton, solving both equations simultaneously),   !
   ! or the original method (Brent method with bracketing on gsw).  This is assigned at    !
   ! sub-routine init_physiology_params (ed_params.f90).                                   !
   !---------------------------------------------------------------------------------------!
   logical :: new_c3_solver
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters for the 2-D Newton's method with line searching).  These number are    !
   ! assigned at sub-routine init_physiology_params (ed_params.f90).                       !
   !---------------------------------------------------------------------------------------!
   !----- Bounds for the new C3 solver. ---------------------------------------------------!
   real(kind=4)  :: c34smin_ci   ! Minimum carbon dioxide concentration          [ mol/mol]
   real(kind=4)  :: c34smax_ci   ! Maximum carbon dioxide concentration          [ mol/mol]
   real(kind=4)  :: c34smin_gsw  ! Minimum conductivity                          [     m/s]
   real(kind=4)  :: c34smax_gsw  ! Maximum conductivity                          [     m/s]
   !----- Relative scale for the nudging the guess when the solver gets stuck. ------------!
   real(kind=4)  :: nudgescal  
   !---------------------------------------------------------------------------------------!
   !      The parameter used as a fraction of the average rate of decrease of fn2 in the   !
   ! line search method.  Press' suggested value is 1.e-4.                                 !
   !---------------------------------------------------------------------------------------!
   real(kind=4)  :: alfls
   !----- Maximum normalised value of a step for a line search. ---------------------------!
   real(kind=4)  :: normstmax
   !----- A small number that is almost the floating point accuracy, but somewhat larger. -!
   real(kind=4)  :: tinynum
   !---------------------------------------------------------------------------------------!
   !      A large number that can still be squared without causing overflow.  Usually      !
   ! the reciprocal of tinynum.                                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=4)  :: hugenum
   !----- Maximum number of iterations before giving up. ----------------------------------!
   integer       :: maxmdne
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

end module physiology_coms
!==========================================================================================!
!==========================================================================================!

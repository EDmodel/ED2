!============================= Change Log =================================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This module contains structures that will set up the global variable table in BRAMS. !
!------------------------------------------------------------------------------------------!
module var_tables

   !----- Maximum number of variables of all types (3d + 2d + leaf + CARMA + Grell). ------!
   integer, parameter :: maxvars=3200
   !---------------------------------------------------------------------------------------!



   !----- Define data type for main variable table. ---------------------------------------!
   type var_tables_r
      
      real, dimension(:), pointer :: var_p     ! Pointer to instantaneous variable
      real, dimension(:), pointer :: var_m     ! Pointer to mean variable
      integer                     :: npts      ! Number of points
      integer                     :: idim_type ! Variable type
      integer                     :: ihist     ! Variable is written to history
      integer                     :: ianal     ! Variable is written to analysis
      integer                     :: imean     ! Variable is written to mean analysis
      integer                     :: ilite     ! Variable is written to light analysis
      integer                     :: impti     ! Variable is sent out during initialisation
      integer                     :: impt1     ! Lateral bnd. cond. at regular time step
      integer                     :: impt2     ! Lateral bnd. cond. at acoustic time step
      integer                     :: impt3     ! Variable is to be sent back at analysis
      integer                     :: iadvt     ! Advection time step (T variables)
      integer                     :: iadvu     ! Advection time step (U variables)
      integer                     :: iadvv     ! Advection time step (V variables)
      integer                     :: iadvw     ! Advection time step (W variables)
      integer                     :: irecycle  ! Recycle variable 
      character(len=16)           :: name      ! Variable name
   end type var_tables_r
   !---------------------------------------------------------------------------------------!


   !----- Main variable table allocated to (maxvars,maxgrds). -----------------------------!
   type(var_tables_r), dimension(:,:), allocatable :: vtab_r
   !---------------------------------------------------------------------------------------!


   !----- "nvgrids" is "ngrids", for convenience. -----------------------------------------!
   integer :: nvgrids
   !---------------------------------------------------------------------------------------!


   !----- Number of variables for each grid, allocated to "ngrids". -----------------------!
   integer, dimension(:), allocatable :: num_var
   !---------------------------------------------------------------------------------------!



   !----- Define data type for scalar variable table. -------------------------------------!
   type scalar_table
      real, dimension(:), pointer :: var_p
      real, dimension(:), pointer :: var_t
      character (len=16) :: name
      real, dimension(:), pointer :: a_var_p
      real, dimension(:), pointer :: a_var_t
   end type scalar_table
   !---------------------------------------------------------------------------------------!



   !----- Scalar variable table allocated to (maxsclr,maxgrds). ---------------------------!
   type(scalar_table), dimension(:,:), allocatable :: scalar_tab
   !---------------------------------------------------------------------------------------!


   !----- Number of scalars for each grid, allocated to "ngrids". -------------------------!
   integer, dimension(:), allocatable :: num_scalar
   !---------------------------------------------------------------------------------------!


end module var_tables
!==========================================================================================!
!==========================================================================================!


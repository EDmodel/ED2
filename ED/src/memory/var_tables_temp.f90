!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module var_tables

    !    Maximum number of variables of all types (3d + 2d + leaf)
integer, parameter :: maxvars=1050

    !    Define data type for main variable table

type var_tables_r
   
   real, pointer :: var_p,var_m
   integer :: idim_type
   integer,dimension(6) :: dimensions
   integer :: ndims
   integer :: ihist,ianal,imean,ilite,impti,impt1,impt2,impt3,irecycle
   character (len=16) :: name
   ![ED2-MLO
   integer :: imont,idail
   !ED2-MLO]

   
end type var_tables_r

    !    Main variable table allocated to (maxvars,maxgrds)
type(var_tables_r), allocatable :: vtab_r(:,:)

    !    "nvgrids" is "ngrids", for convenience
integer :: nvgrids

    !    number of variables for each grid, allocated to "ngrids"
integer, allocatable :: num_var(:)


    !    Define data type for scalar variable table

type scalar_table
   
   real, pointer :: var_p,var_t
   character (len=16) :: name
   ! ALF
   real, pointer :: a_var_p(:), a_var_t(:)

end type

    !    Scalar variable table allocated to (maxsclr,maxgrds)
type(scalar_table), allocatable :: scalar_tab(:,:)


    !    number of scalars for each grid, allocated to "ngrids"
integer, allocatable :: num_scalar(:)


End Module var_tables


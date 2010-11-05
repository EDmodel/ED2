!==========================================================================================!
!==========================================================================================!
!    This module contains some structures used by the model to solve the photosynthesis.   !
!------------------------------------------------------------------------------------------!
module c34constants
   implicit none

   type farqdata
      real :: D0
      real :: alpha
      real :: gamma
      real :: m
      real :: b
   end type farqdata

   type metdat
      real :: ea
      real :: ca
      real :: rn
      real :: tl
      real :: par
      real :: gbc
      real :: gbw
      real :: ta
      real :: el
      real :: compp
      real :: eta
      real :: gbci
   end type metdat

   type glim
      real :: sigma
      real :: rho
      real :: vm
      real :: tau
      real :: nu
      real :: k1
      real :: k2
   end type glim

   type solution
      real, dimension(2,2) :: es
      real, dimension(2,2) :: ci
      real, dimension(2,2) :: cs
      real, dimension(2,2) :: a
      real, dimension(2,2) :: gsw
      real                 :: eps
      real                 :: gsw2_1st
      real                 :: ci2_1st
      integer              :: ninterval
   end type solution

   type stoma_data
      integer :: recalc=1   !THIS SHOULD BE INIT IN ED_PARAMS
      real    :: T_L
      real    :: e_A
      real    :: PAR
      real    :: rb_factor
      real    :: prss
      real    :: phenology_factor
      real    :: gsw_open
      integer :: ilimit
      
      real    :: T_L_residual
      real    :: e_a_residual
      real    :: par_residual
      real    :: rb_residual
      real    :: prss_residual
      real    :: leaf_residual
      real    :: gsw_residual
   end type stoma_data

   !------ The number of stomatal attributes. ---------------------------------------------!
   integer, parameter :: n_stoma_atts = 16
   !---------------------------------------------------------------------------------------!
end module c34constants
!==========================================================================================!
!==========================================================================================!

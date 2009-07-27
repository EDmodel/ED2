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
     real :: ea,ca,rn,tl,par,gbc,gbw,ta,el,compp,eta,gbci
  end type metdat

  type glim
     real :: sigma,rho,vm,tau,nu,k1,k2
  end type glim

  type solution
     real, dimension(2,2) :: es,ci,cs,a,gsw
     real :: eps
     integer :: ninterval
  end type solution

  Type stoma_data

     integer :: recalc=1   !THIS SHOULD BE INIT IN ED_PARAMS
     real :: T_L
     real :: e_A
     real :: PAR
     real :: rb_factor
     real :: prss
     real :: phenology_factor
     real :: gsw_open
     integer :: ilimit
     
     real :: T_L_residual
     real :: e_a_residual
     real :: par_residual
     real :: rb_residual
     real :: prss_residual
     real :: leaf_residual
     real :: gsw_residual
     
  End Type stoma_data

  ! the number of stomatal attributes
  integer,parameter :: n_stoma_atts = 16



end module c34constants

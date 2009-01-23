!==========================================================================================!
!==========================================================================================!
!    This module contains several parameters used by the Runge-Kutta integrator scheme.    !
! This module is intended to host the numerical method variables only, variables related   !
! to physical, hydrological or ecological properties should be stored somewhere else.      !
!                                                                                          !
!   MLO. I attempted to describe some variables, not sure about all of them, if the        !
!        definition is wrong, please correct it. Thanks!                                   !
!------------------------------------------------------------------------------------------!
module rk4_coms
   implicit none



   !---------------------------------------------------------------------------------------!
   !    The following variables are parameters that do not depend on the specific run.     !
   !---------------------------------------------------------------------------------------!
   
   !----- Maximum number of intermediate steps --------------------------------------------!
   integer, parameter :: maxstp = 100000000

   !----- Small number, to avoid singularities --------------------------------------------!
   real   , parameter :: tiny_offset = 1.0e-20

   !----- rk4eps is the desired accuracy (former eps), and rk4epsi is its reciprocal ------!
   real   , parameter :: rk4eps  = 0.01
   real   , parameter :: rk4epsi = 1./rk4eps
   
   !----- hmin is the minimum step size ---------------------------------------------------!
   real   , parameter :: hmin = 1.e-9
   !---------------------------------------------------------------------------------------!

   real   , parameter :: safety = 0.9
   real   , parameter :: pgrow = -0.2
   real   , parameter :: pshrnk = -0.25
   real   , parameter :: errcon = 1.89e-4


   !---------------------------------------------------------------------------------------!
   !     Integration limits for time. Since 'derivs' do not explicitly depend on time it   !
   ! doesn't really matter what this is as long as tend-tbeg makes sense.                  !
   !---------------------------------------------------------------------------------------!
   real   :: tbeg
   real   :: tend
   real   :: dtrk4
   real   :: dtrk4i

   !---------------------------------------------------------------------------------------!
   
end module rk4_coms

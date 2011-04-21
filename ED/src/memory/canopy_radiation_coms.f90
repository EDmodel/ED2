!==========================================================================================!
!==========================================================================================!
!    This module contains several parameters used in the canopy radiation solver.          !
!                                                                                          !
!  IMPORTANT: Do not initialize non-parameters in their modules - not all compilers will   !
!             actually initialize them.  Instead, assign them at init_can_rad_params sub-  !
!             routine (ed_params.f90).                                                     !
!------------------------------------------------------------------------------------------!
module canopy_radiation_coms

   use ed_max_dims , only : n_pft
   implicit none 


   !---------------------------------------------------------------------------------------!
   !     Leaf angle distribution parameter (dimensionless).  Let mu' be the cosine of leaf !
   ! angle and G(mu') be the distribution of mu'.  Then, mubar = (integral from 0 to 1)    !
   ! (d mu'   mu' / G(mu')).  See, for example, Dickinson 1983.                            !
   !---------------------------------------------------------------------------------------!
   real(kind=8) :: mubar
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Fraction of solar radiation in the PAR band.  Used when you don't know the        !
   ! direct/diffuse breakdown.                                                             !
   !---------------------------------------------------------------------------------------!
   real :: visible_fraction
   !---------------------------------------------------------------------------------------!


   !----- Fraction of direct solar radiation in the PAR band. -----------------------------!
   real :: visible_fraction_dir
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Fraction of diffuse solar radiation in the PAR band.  Used when you don't know    !
   ! the direct/diffuse breakdown.                                                         !
   !---------------------------------------------------------------------------------------!
   real :: visible_fraction_dif
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     The following three values are based on Dickinson (1983).                         !
   !---------------------------------------------------------------------------------------!
   real :: leaf_reflect_nir  !----- Fraction of scattered NIR that is reflected. ----------!
   real :: leaf_trans_nir    !----- Fraction of scattered NIR that is transmitted. --------!
   real :: leaf_scatter_nir  !----- Sum of reflected plus scattered NIR. ------------------!

   !---------------------------------------------------------------------------------------!
   !     The following two values are from Baldocchi et al.                                !
   !---------------------------------------------------------------------------------------!
   !----- Fraction of scattered PAR that is reflected (temperate trees). ------------------!
   real, parameter :: leaf_reflect_vis_temperate = 0.11
   !----- Fraction of scattered PAR that is transmitted (temperate trees). ----------------!
   real, parameter :: leaf_trans_vis_temperate = 0.16  
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     The following two values are from Poorter et al.                                  !
   !---------------------------------------------------------------------------------------!
   !----- Fraction of scattered PAR that is reflected (tropical trees). -------------------!
   real, parameter :: leaf_reflect_vis_tropics = 0.062  
   !----- Fraction of scattered PAR that is transmitted (tropical trees). -----------------!
   real, parameter :: leaf_trans_vis_tropics = 0.028  
   !---------------------------------------------------------------------------------------!

   !----- Backscatter parameter for diffuse NIR. ------------------------------------------!
   real :: diffuse_backscatter_nir 

   !----- Sum of reflected plus scattered PAR. --------------------------------------------!
   real, dimension(n_pft) :: leaf_scatter_vis

   !----- Backscatter parameter for diffuse PAR. ------------------------------------------!
   real, dimension(n_pft) :: diffuse_backscatter_vis

   !----- Emissivity of the vegetation. ---------------------------------------------------!
   real(kind=8), dimension(n_pft) :: emis_v
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Crown model flag (set at the namelist):                                           !
   !     0. Original;                                                                      !
   !     1. Finite-crown mixing model.                                                     !
   !---------------------------------------------------------------------------------------!
   integer :: crown_mod
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The following variables control whether to call things that should be called      !
   ! when there is still some light.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4)    :: rshort_twilight_min
   real(kind=4)    :: cosz_min
   real(kind=8)    :: cosz_min8
   !---------------------------------------------------------------------------------------!
end module canopy_radiation_coms
!==========================================================================================!
!==========================================================================================!

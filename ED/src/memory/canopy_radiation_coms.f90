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

   !----- Converts PAR radiation from  watts to Einsteins (units are Ein/watt of PAR). ----!
   real, parameter :: Watts2Ein = 4.6e-6
   real, parameter :: Ein2Watts = 1./Watts2Ein
   !---------------------------------------------------------------------------------------!


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
   !    The following values are used to decide whether the radiation should be called or  !
   ! not.  For a simple test, I am currently using TAI rather than LAI, but I'm aware that !
   ! branches and twigs probably should be treated differently.                            !
   !---------------------------------------------------------------------------------------!
   real :: lai_min       !----- Minimum LAI used in the radiative transfer scheme. --------!
   real :: tai_min       !----- Minimum TAI used in the radiative transfer scheme. --------!
   !---------------------------------------------------------------------------------------!
   ! Fraction of bleaf_max below which we skip photosynthesis, radiation and heat balance. !
   !---------------------------------------------------------------------------------------!
   real :: blfac_min
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     The following variables are minimum values that are considered acceptable.
   !---------------------------------------------------------------------------------------!
   real :: rlong_min     !----- Minimum allowed downward longwave radiation [W/m²] --------!
   real :: veg_temp_min  !----- Minimum allowed vegetation temperature      [   K] --------!
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Crown model flag (set at the namelist):                                           !
   !     0. Original;                                                                      !
   !     1. Finite-crown mixing model.                                                     !
   !---------------------------------------------------------------------------------------!
   integer :: crown_mod
   !---------------------------------------------------------------------------------------!

end module canopy_radiation_coms
!==========================================================================================!
!==========================================================================================!

Module canopy_radiation_coms

  ! DO NOT INITIALIZE NON-PARAMETERS IN THEIR MODULES - NOT ALL COMPILERS WILL ACTUALLY INITIALIZE THEM

  ! See "initialize_canopy_radiation_params" for initial values

use max_dims, only: n_pft

implicit none 
real :: mubar            ! Leaf angle distribution parameter (dimensionless).
                         !Let mu' be the cosine of leaf angle and G(mu') be the 
                         !distribution of mu'.  Then, mubar = (integral from 0 to 1) 
                         !(d mu'   mu' / G(mu')).  See, for example, Dickinson 1983.

real, parameter :: Watts2Ein = 4.6e-6 ! Converts PAR radiation from  watts 
                                      !to Einsteins (units are Ein/watt of PAR).

real :: visible_fraction        ! fraction of solar radiation in the PAR band.  
                                !Used when you don't know the direct/diffuse breakdown.

real :: visible_fraction_dir    ! fraction of direct solar radiation in the PAR band.

real :: visible_fraction_dif     ! fraction of diffuse solar radiation in the PAR band. 
                                 ! Used when you don't know the direct/diffuse breakdown.

! The following two values are based on Dickinson (1983)
real :: leaf_reflect_nir  !  fraction of scattered NIR that is reflected.  

real :: leaf_trans_nir    !  fraction of scattered NIR that is transmitted.

real :: leaf_scatter_nir ! sum of reflected plus scattered NIR

! The following two values are from Baldocchi et al.
real, parameter :: leaf_reflect_vis_temperate = 0.11  !  fraction of scattered 
                                                      ! PAR that is reflected.  
                                                  ! Value obtained for temperate trees.

real, parameter :: leaf_trans_vis_temperate = 0.16  
!  fraction of scattered PAR that is transmitted.  Value obtained for temperate trees.

! The following two values are from Poorter et al.
real, parameter :: leaf_reflect_vis_tropics = 0.062  
!  fraction of scattered PAR that is reflected.  Value obtained for tropical trees.

real, parameter :: leaf_trans_vis_tropics = 0.028  
!  fraction of scattered PAR that is transmitted.  Value obtained for tropical trees.

real :: diffuse_backscatter_nir !  Backscatter parameter for diffuse NIR

real :: lai_min       ! Minimum LAI used in the radiative transfer scheme.

real, dimension(n_pft) :: leaf_scatter_vis  ! sum of reflected plus scattered PAR

real, dimension(n_pft) :: diffuse_backscatter_vis !  Backscatter parameter for diffuse PAR

real(kind=8), dimension(n_pft) :: emis_v  ! Emissivity of the vegetation

real :: rlong_min     !! minimum allowed downward longwave radiation (w/m2)
real :: veg_temp_min  !! minimum allowed vegetation temperature

End Module canopy_radiation_coms

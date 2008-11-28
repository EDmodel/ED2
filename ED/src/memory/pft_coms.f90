module pft_coms

use max_dims, only: n_pft
     !================================================================
     ! Plant functional type:
     ! 1 - C4 grass
     ! 2 - tropical early successional tree
     ! 3 - tropical mid successional tree
     ! 4 - tropical late successional tree
     ! 5 - C3 grass
     ! 6 - northern pines (temperate)
     ! 7 - southern pines (temperate)
     ! 8 - late successional conifers
     ! 9 - early successional cold-deciduous hardwood
     ! 10 - mid successional cold-deciduous hardwood
     ! 11 - late successional cold-deciduous hardwood
     ! 12 - c3 pasture
     ! 13 - c3 crop (e.g.,wheat, rice, soybean) 
     ! 14 - c4 pasture
     ! 15 - c4 crop (e.g.,corn/maize)


! DO NOT INITIALIZE NON-PARAMETERS IN THEIR MODULES - NOT ALL COMPILERS WILL ACTUALLY INITIALIZE THEM

!--------------------
! GENERAL
!--------------------
integer, dimension(n_pft) :: include_these_pft 
integer, dimension(n_pft) :: include_pft ! Set to 1 if you want to include this PFT; 0 otherwise.

! This is the list of grass PFTs that may be included in agricultural patches. Only PFTs included here and 
! at the include_these_pft will be used for agricultural patches.
integer, dimension(n_pft) :: grass_pft

! This flag specifies what non-agricutural PFTs (i.e., grass)  can grow on agriculture patches.  Set 
! to 1 if you want to include this PFT on agriculture patches
integer, dimension(n_pft) :: include_pft_ag

! This flag determines what to do at the PFT initialization. This option is ignored for 
! near-bare ground simulations.
! 0. Stop if a undesired PFT is at the restart file;
! 1. Include the PFT in the list of usable pfts (recompute include_these_pft and include_pft)
! 2. Ignore the cohort and keep going.
integer :: pft_1st_check

!--------------------
! PHOTOSYNTHESIS AND STOMATAL CONDUCTANCE
!--------------------
real, dimension(n_pft) :: D0 ! (mol H2O/mol air).   Stomata begin to rapidly close once the difference between intercellular and boundary layer H2O mixing ratios exceed this value.

real, dimension(n_pft) :: Vm_low_temp ! Temperature (C) below which leaf metabolic activity begins to rapidly decline.

real, dimension(n_pft) :: Vm0  ! maximum photosynthetic capacity at a reference temperature. (umol/m2/s)

real, dimension(n_pft) :: stomatal_slope !  Slope of the Ball/Berry stomatal conductance-photosynthesis relationship.

real, dimension(n_pft) :: cuticular_cond ! Intercept of the Ball/Berry stomatal conductance relationship.  (umol/m2/s)

real, dimension(n_pft) :: quantum_efficiency  ! efficiency of using PAR to fix CO2 (mol CO2 / Einstein)

integer, dimension(n_pft) :: photosyn_pathway ! specifies photosynthetic pathway.  3 corresponds to C3, 4 corresponds to C4.

!--------------------
! RESPIRATION AND TURNOVER
!--------------------
real, dimension(n_pft) :: growth_resp_factor  !  Determines level of growth respiration.  Starting with accumulated photosynthesis (P), leaf (Rl) and root respiration (Rr) are first subtracted.  Then, the growth respiration = (growth_resp_factor) * (P - Rl - Rr).

real, dimension(n_pft) :: leaf_turnover_rate ! This is 1/(leaf life span).  Units are 1/year.

real, dimension(n_pft) :: root_turnover_rate ! This is 1/(fine root life span).  Units are 1/year.

real, dimension(n_pft) :: dark_respiration_factor ! Sets the rate of dark (i.e., leaf) respiration.  (dimensionless; it is relative to Vm0.)

real, dimension(n_pft) :: storage_turnover_rate ! Turnover rate of plant storage pools (1/year).

real, dimension(n_pft) :: root_respiration_factor ! (umol CO2)/(kg fine roots)/second

!--------------------
! MORTALITY AND SURVIVORSHIP
!--------------------
real, dimension(n_pft) :: mort1 ! Sets the time scale at which plants out of carbon balance suffer mortality (1/year)

real, dimension(n_pft) :: mort2 ! Determines how poor the carbon balance needs to be before plants suffer large mortality rates.

real, dimension(n_pft) :: mort3 ! Density-independent mortality rate (1/years)

real :: frost_mort    ! Determines how rapidly trees die if it is too cold for them (1/year)

real, dimension(n_pft) :: seedling_mortality ! Fraction of seedlings that suffer mortality without becoming a recruit. 

real, dimension(n_pft) :: treefall_s_gtht ! Survivorship fraction for trees with heights greater than treefall_hite_threshold (see disturbance_coms.f90)

real, dimension(n_pft) :: treefall_s_ltht ! Survivorship fraction for trees with heights less than treefall_hite_threshold (see disturbance_coms.f90)

real, dimension(n_pft) :: plant_min_temp ! Below this temperature, mortality rapidly increases.

!--------------------
! NITROGEN AND WATER REQUIREMENTS  -- See "initialize_pft_nitro_params"
!--------------------

real :: c2n_slow               ! Carbon to Nitrogen ratio, slow pool.
real :: c2n_structural         ! Carbon to Nitrogen ratio, structural pool.
real :: c2n_storage            ! Carbon to Nitrogen ratio, storage pool.
real :: c2n_stem               ! Carbon to Nitrogen ratio, structural stem.
real :: l2n_stem               ! Carbon to Nitrogen ratio, structural stem.
real, dimension(n_pft) :: c2n_leaf ! Leaf carbon to nitrogen ratio
real, dimension(n_pft) :: c2n_recruit ! Recruit carbon to nitrogen ratio
real :: C2B                    !  Carbon-to-biomass ratio of plant tissues.
real :: agf_bs                 ! fraction of structural stem that is assumed to be above ground.
real :: plant_N_supply_scale   ! Supply coefficient for plant nitrogen uptake (m2/(kgC fine root)/day)
real, dimension(n_pft) :: water_conductance  ! Supply coefficient for plant water uptake:  (kg H2O) (m2 ground) / (m3 H2O) / (kgC root) / (seconds)

!--------------------
! ALLOCATION AND ALLOMETRY
!--------------------
real, dimension(n_pft) :: rho  ! wood density.  Used only for tropical PFTs (kg/m3).
real, dimension(n_pft) :: SLA ! specific leaf area (m2 leaf / kg C)
real, dimension(n_pft) :: q ! Ratio of (kg fine roots) / (kg leaves)
real, dimension(n_pft) :: qsw ! Ratio of (kg sapwood) / (kg leaves)
real, dimension(n_pft) :: hgt_min ! minimum height of an individual (m)
real, dimension(n_pft) :: b1Ht  !  DBH-height allometry intercept (m).  Temperate PFTs only.
real, dimension(n_pft) :: b2Ht  !  DBH-height allometry slope (1/cm).  Temperate PFTs only.
real, dimension(n_pft) :: b1Bs  !  DBH-stem allometry intercept (kg stem biomass / plant * cm^{-b2Bs}).  Temperate PFTs only.
real, dimension(n_pft) :: b2Bs  !  DBH-stem allometry slope (dimensionless).  Temperate PFTs only.
real, dimension(n_pft) :: b1Bl  !  DBH-leaf allometry intercept (kg leaf biomass / plant * cm^{-b2Bl}).  Temperate PFTs only.
real, dimension(n_pft) :: b2Bl  !  DBH-leaf allometry slope (dimensionless).  Temperate PFTs only.
real, dimension(n_pft) :: max_dbh !  Maximum DBH attainable by this PFT (cm)

!--------------------
! LEAF HABIT AND PHYSICAL PROPERTIES
!--------------------
integer, dimension(n_pft) :: phenology ! Indicates leaf habit.  0 - evergreen coniferous; 1 - drought deciduous; 2 - cold deciduous.

real, dimension(n_pft) :: clumping_factor ! 0-1 factor indicating degree of clumpiness of leaves and shoots.

real, dimension(n_pft) :: leaf_width  ! leaf width used to compute the aerodynamic resistance (m).
! NEW PARAMETERS 11-26-08

real, dimension(n_pft) :: c_grn_leaf_dry ! Specific heat capacity of dry leaf biomass (J/kg/K)

real, dimension(n_pft) :: c_ngrn_biom_dry! Specific heat capacity of dry non-green biomass(J/kg/K)

real, dimension(n_pft) :: wat_dry_ratio_grn  ! Ratio of water to dry mass in green leaves

real, dimension(n_pft) :: wat_dry_ratio_ngrn ! Ratio of water to dry mass in non-green biomass


real :: hcap_stem_fraction               ! Fraction of structural biomass, that
                                         ! that is included in the calculation of vegetation
                                         ! heat capacity, used for dianosing leaf temerpatures
                                         ! from vegetation energy. Ie %veg_energy(ico) and
                                         ! %veg_temp(ico)


!--------------------
! REPRODUCTION AND RECRUITMENT
!--------------------
real, dimension(n_pft) :: r_fract ! fraction of (positive) carbon balance devoted to reproduction.

real, dimension(n_pft) :: seed_rain ! external input of seeds (kgC/m2/year)

real, dimension(n_pft) :: nonlocal_dispersal !  Fraction of seed dispersal that is gridcell-wide.

real, dimension(n_pft) :: repro_min_h ! Minimum height plants need to attain before allocating to reproduction.



end module pft_coms

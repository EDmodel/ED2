!==========================================================================================!
!==========================================================================================!
!   This module contains a list of plant-functional type dependent properties.             !
!                                                                                          !
! IMPORTANT: DO NOT INITIALIZE PARAMETERS IN THEIR MODULES - NOT ALL COMPILERS WILL        !
!            ACTUALLY INITIALIZE THEM.  See "init_pft_*_coms" (ed_params.f90) to check     !
!            the default values.                                                           !
!==========================================================================================!
!==========================================================================================!
module pft_coms

   use ed_max_dims, only: n_pft

   !---------------------------------------------------------------------------------------!
   ! PFT | Name                            | Grass | Liana | Tropical | Savannah | Conifer !
   !-----+---------------------------------+-------+-------+----------+----------+---------!
   !   1 | C4 grass                        |   yes |    no |      yes |       no |      no !
   !   2 | Early tropical                  |    no |    no |      yes |       no |      no !
   !   3 | Mid tropical                    |    no |    no |      yes |       no |      no !
   !   4 | Late tropical                   |    no |    no |      yes |       no |      no !
   !   5 | Temperate C3 grass              |   yes |    no |       no |       no |      no !
   !   6 | Northern pines                  |    no |    no |       no |       no |     yes !
   !   7 | Southern pines                  |    no |    no |       no |       no |     yes !
   !   8 | Late conifers                   |    no |    no |       no |       no |     yes !
   !   9 | Early temperate deciduous       |    no |    no |       no |       no |      no !
   !  10 | Mid temperate deciduous         |    no |    no |       no |       no |      no !
   !  11 | Late temperate deciduous        |    no |    no |       no |       no |      no !
   !  12 | Early tropical savannah         |    no |    no |      yes |       no |      no !
   !  13 | Mid tropical savannah           |    no |    no |      yes |       no |      no !
   !  14 | Late tropical savannah          |    no |    no |      yes |       no |      no !
   !  15 | Araucaria                       |    no |    no |      yes |       no |     yes !
   !  16 | Tropical C3 grass               |   yes |    no |      yes |       no |      no !
   !  17 | Liana                           |    no |   yes |      yes |       no |     yes !
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!
   !     Variables that are provided by the user's namelist.  They control PFT habits such !
   ! as which PFT should be used for agriculture, which one goes for forest plantation.    !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     This variable is provided by the user through namelist, and contains the list of  !
   ! PFTs he or she wants to use.                                                          !
   !---------------------------------------------------------------------------------------!
   integer, dimension(n_pft) :: include_these_pft

   !---------------------------------------------------------------------------------------!
   !     This flag determines what to do at the PFT initialization.  This option is        !
   ! ignored for near-bare ground simulations:                                             !
   !  0. Stop if a undesired PFT is at the restart file;                                   !
   !  1. Include the PFT in the list of usable pfts (recompute include_these_pft and       !
   !     include_pft);                                                                     !
   !  2. Ignore the cohort and keep going.                                                 !
   !---------------------------------------------------------------------------------------!
   integer :: pft_1st_check

   !---------------------------------------------------------------------------------------!
   !     These are the flags that indicate which PFTs should be used for pastures, agri-   !
   ! culture and plantation stockS. They are currently a single PFT, but they should       !
   ! become vectors eventually (and thus multiple PFTs can be used...).                    !
   !---------------------------------------------------------------------------------------!
   integer :: pasture_stock 
   integer :: agri_stock 
   integer :: plantation_stock
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     The following parameters aren't PFT-dependent, but they are used to determine     !
   ! PFT-dependent properties.                                                             !
   !---------------------------------------------------------------------------------------!
   !----- Carbon-to-biomass ratio of plant tissues. ---------------------------------------!
   real :: C2B
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    The following variables are flags that control which PFTs mare used in general,    !
   ! for agriculture and grasses.                                                          !
   !---------------------------------------------------------------------------------------!
   
   
   !---------------------------------------------------------------------------------------!
   !    This is the list of PFTs that are included.  0 means off, 1 means on.              !
   !---------------------------------------------------------------------------------------!
   logical, dimension(n_pft) :: include_pft
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    This flag specifies which PFTs are allowed to grow on pasture patches.             !
   ! Zero means the PFT is forbidden, 1 means that the PFT is allowed.                     !
   !---------------------------------------------------------------------------------------!
   logical, dimension(n_pft) :: include_pft_pt
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    This flag specifies which PFTs are allowed to grow on agriculture patches.         !
   ! Zero means the PFT is forbidden, 1 means that the PFT is allowed.                     !
   !---------------------------------------------------------------------------------------!
   logical, dimension(n_pft) :: include_pft_ag
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    This flag specifies which PFTs are allowed to grow on forest plantation patches.   !
   ! Zero means the PFT is forbidden, 1 means that the PFT is allowed.                     !
   !---------------------------------------------------------------------------------------!
   logical, dimension(n_pft) :: include_pft_fp
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    The following logical flags will describe the main life forms (grasses, lianas,    !
   ! conifers) as well as the likely landscape where they are to be found                  !
   ! (tropical/temperate, savannah/forest).  Note that ED2 does not restrict location of   !
   ! any PFT: if you really want, you could include Eastern hemlocks in an Amazonian run,  !
   ! or a Cecropia in the tundra.                                                          !
   !---------------------------------------------------------------------------------------!
   logical, dimension(n_pft)    :: is_tropical
   logical, dimension(n_pft)    :: is_savannah
   logical, dimension(n_pft)    :: is_conifer
   logical, dimension(n_pft)    :: is_grass
   logical, dimension(n_pft)    :: is_liana
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Photosynthesis and stomatal conductance properties.                               !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !   Stomata begin to rapidly close once the difference between intercellular and        !
   ! boundary layer H2O mixing ratios exceed this value. [mol_H2O/mol_air].                !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: D0
   !---------------------------------------------------------------------------------------!

   !----- Temperature [C] below which leaf metabolic activity begins to rapidly decline. --!
   real, dimension(n_pft) :: Vm_low_temp 

   !----- Temperature [C] above which leaf metabolic activity begins to rapidly decline. --!
   real, dimension(n_pft) :: Vm_high_temp 

   !----- Decay factors for the exponential correction. -----------------------------------!
   real, dimension(n_pft) :: Vm_decay_elow
   real, dimension(n_pft) :: Vm_decay_ehigh

   !----- Maximum carboxylation rate at the reference temperature [umol/m2/s]. ------------!
   real, dimension(n_pft) :: Vm0 
   !----- Parameters used by the model that predicts Vm0 based on SLA. --------------------!
   real, dimension(n_pft) :: Vm0_v0
   real, dimension(n_pft) :: Vm0_v1
   !----- Parameters used by the new stomatal conductance model (Farquhar-Katul). ---------!
   real, dimension(n_pft) :: Vcmax25      ! Vm0 at 25degC
   !----- Parameters currently used only by trait plasticity. -----------------------------!
   real, dimension(n_pft) :: kplastic_vm0 ! Expansion factor for Vm0 (extinction if < 0).
   !----- Exponent for Vm in the Arrhenius equation [K]. ----------------------------------!
   real, dimension(n_pft) :: Vm_hor 

   !----- Base (Q10 term) for Vm in Collatz equation. -------------------------------------!
   real, dimension(n_pft) :: Vm_q10

   !----- Temperature [C] below which leaf metabolic activity begins to rapidly decline. --!
   real, dimension(n_pft) :: Jm_low_temp 

   !----- Temperature [C] above which leaf metabolic activity begins to rapidly decline. --!
   real, dimension(n_pft) :: Jm_high_temp 

   !----- Decay factors for the exponential correction. -----------------------------------!
   real, dimension(n_pft) :: Jm_decay_elow
   real, dimension(n_pft) :: Jm_decay_ehigh

   !----- Maximum electron transport rate at the reference temperature [umol/m2/s]. -------!
   real, dimension(n_pft) :: Jm0 

   !----- Parameters used by the new stomatal conductance model (Farquhar-Katul). ---------!
   real, dimension(n_pft) :: Jmax25      ! Jm0 at 25degC

   !----- Exponent for Jm in the Arrhenius equation [K]. ----------------------------------!
   real, dimension(n_pft) :: Jm_hor 

   !----- Base (Q10 term) for Jm in Collatz equation. -------------------------------------!
   real, dimension(n_pft) :: Jm_q10

   !----- Maximum triose phosphate utilisation rate at the ref. temperature [umol/m2/s]. --!
   real, dimension(n_pft) :: TPm0

   !----- Temperature [C] below which leaf metabolic activity begins to rapidly decline. --!
   real, dimension(n_pft) :: Rd_low_temp 

   !----- Temperature [C] above which leaf metabolic activity begins to rapidly decline. --!
   real, dimension(n_pft) :: Rd_high_temp 

   !----- Decay factors for the exponential correction. -----------------------------------!
   real, dimension(n_pft) :: Rd_decay_elow
   real, dimension(n_pft) :: Rd_decay_ehigh

   !----- Maximum respiration factor at the reference temperature [umol/m2/s]. ------------!
   real, dimension(n_pft) :: Rd0

   !----- Exponent for Rd in the Arrhenius equation [K]. ----------------------------------!
   real, dimension(n_pft) :: Rd_hor 

   !----- Base (Q10 term) for respiration in Collatz equation. ----------------------------!
   real, dimension(n_pft) :: Rd_q10

   !----- Slope of the Ball/Berry stomatal conductance-photosynthesis relationship. -------!
   real, dimension(n_pft) :: stomatal_slope

   !----- Intercept of the Ball/Berry stomatal conductance relationship [umol/m2/s]. ------!
   real, dimension(n_pft) :: cuticular_cond

   !----- Efficiency of using PAR to fix CO2 [ ----]. -------------------------------------!
   real, dimension(n_pft) :: quantum_efficiency

   !----- Curvature parameter for determining the electron transport rate (aka J) [ ---]. -!
   real, dimension(n_pft) :: curvpar_electron

   !----- Quantum yield of photosystem II [ ----]. ----------------------------------------!
   real, dimension(n_pft) :: qyield_psII

   !----- Specifies photosynthetic pathway.  3 corresponds to C3, 4 corresponds to C4. ----!
   integer, dimension(n_pft) :: photosyn_pathway
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Respiration and turnover properties.                                              !
   !---------------------------------------------------------------------------------------!
  
   !---------------------------------------------------------------------------------------!
   !   This variable determines level of growth respiration.  Starting with accumulated    !
   ! photosynthesis (P), leaf (Rl) and root respiration (Rr) are first subtracted.  Then,  !
   ! the growth respiration = (growth_resp_factor) * (P - Rl - Rr).                        !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: growth_resp_factor

   !----- This is the inverse of leaf life span [1/year]. ---------------------------------!
   real, dimension(n_pft) :: leaf_turnover_rate
   !----- Parameters that control the leaf turnover rate (TRAIT_PLASTICITY < 0). ----------!
   real, dimension(n_pft) :: eplastic_vm0  ! Expansion/extinction exponents for turnover
   real, dimension(n_pft) :: eplastic_sla  ! (extinction if < 1)
   !----- Parameters that control the leaf turnover rate (TRAIT_PLASTICITY > 0). ----------!
   real, dimension(n_pft) :: kplastic_LL   ! Expansion factor for leaf longevity
   !----- This is the inverse of fine root life span [1/year]. ----------------------------!
   real, dimension(n_pft) :: root_turnover_rate

   !----- This is the inverse of bark life span [1/year]. ---------------------------------!
   real, dimension(n_pft) :: bark_turnover_rate

   !---------------------------------------------------------------------------------------!
   !    This variable sets the rate of dark (i.e., leaf) respiration.  It is dimensionless !
   ! because it is relative to Vm0.  This is only needed in case Rd0 is not provided       !
   ! through xml.                                                                          !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: dark_respiration_factor
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    This variable sets the rate of electron transport.  It is dimensionless            !
   ! because it is relative to Vm0.  This is only needed in case Jm0 is not provided       !
   ! through xml.                                                                          !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: electron_transport_factor
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    This variable sets the rate of triose phosphate utilisation.  It is dimensionless  !
   ! because it is relative to Vm0.  This is only needed in case Jm0 is not provided       !
   ! through xml.                                                                          !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: triose_phosphate_factor
   !---------------------------------------------------------------------------------------!


   !----- Turnover rate of plant storage pools [1/year]. ----------------------------------!
   real, dimension(n_pft) :: storage_turnover_rate 

   !---------------------------------------------------------------------------------------!
   !    This variable sets the contribution of roots to respiration at the reference       !
   ! temperature of 15C.  Its units is umol_CO2/kg_fine_roots/s.                           !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: root_respiration_factor 

   !----- Temperature [C] below which root metabolic activity begins to rapidly decline. --!
   real, dimension(n_pft) :: rrf_low_temp 

   !----- Temperature [C] above which root metabolic activity begins to rapidly decline. --!
   real, dimension(n_pft) :: rrf_high_temp 

   !----- Decay factors for the exponential correction. -----------------------------------!
   real, dimension(n_pft) :: rrf_decay_elow
   real, dimension(n_pft) :: rrf_decay_ehigh

   !----- Exponent for Rr in the Arrhenius equation [K]. ----------------------------------!
   real, dimension(n_pft) :: rrf_hor 

   !----- Base (Q10 term) for respiration in Collatz equation. ----------------------------!
   real, dimension(n_pft) :: rrf_q10

   !----- Temperature [C] below which storage respiration begins to rapidly decline. ------!
   real, dimension(n_pft) :: strf_low_temp 

   !----- Temperature [C] above which storage respiration begins to rapidly decline. ------!
   real, dimension(n_pft) :: strf_high_temp 

   !----- Decay factor for the exponential correction. ------------------------------------!
   real, dimension(n_pft) :: strf_decay_elow
   real, dimension(n_pft) :: strf_decay_ehigh

   !----- Exponent for Rr in the Arrhenius equation [K]. ----------------------------------!
   real, dimension(n_pft) :: strf_hor 

   !----- Base (Q10 term) for respiration in Collatz equation. ----------------------------!
   real, dimension(n_pft) :: strf_q10
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   Mortality and survivorship parameters.                                              !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     This variable shifts the inflection point for the relative carbon balance curve.  !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: mort0

   !---------------------------------------------------------------------------------------!
   !     This variable controls the time scale at which plants out of carbon balance       !
   ! suffer mortality [1/years].                                                           !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: mort1

   !---------------------------------------------------------------------------------------!
   !     This variable determines how poor the carbon balance needs to be before plants    !
   ! suffer large mortality rates.                                                         !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: mort2

   !---------------------------------------------------------------------------------------!
   !     This variable controls the density-independent mortality rate due to ageing       !
   ! [1/years].                                                                            !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: mort3 

   !---------------------------------------------------------------------------------------!
   ! This is the way to initialize mort3 through hard parameterization (ie. non dependent) !
   ! the atlernative is to let mort3 for tropical pfts be controlled by wood density (rho),!
   ! a slope parameter (m3_slope) and a scale parameter (m3_scale).                        !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: mort3_pft_init

   real :: m3_scale
   real :: m3_slope

   !---------------------------------------------------------------------------------------!
   !     This variable sets up the relative carbon balance when plants are experiencing    !
   ! severe stress (i.e., when the maximum carbon balance is negative due to severe light  !
   ! or water stress).                                                                     !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: cbr_severe_stress

   !---------------------------------------------------------------------------------------!
   !     This variable determines how rapidly trees die if it is too cold for them         !
   ! [1/years].                                                                            !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: frost_mort  

   !----- Fraction of seedlings that suffer mortality without becoming a recruit. ---------!
   real, dimension(n_pft) :: seedling_mortality

   !---------------------------------------------------------------------------------------!
   !     Survivorship fraction for trees with heights greater than treefall_hite_threshold !
   ! (see disturbance_coms.f90).                                                           !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: treefall_s_gtht

   !---------------------------------------------------------------------------------------!
   !     Survivorship fraction for trees with heights less than treefall_hite_threshold    !
   ! (see disturbance_coms.f90).                                                           !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: treefall_s_ltht



   !---------------------------------------------------------------------------------------!
   !    Temporary parameters to predict fire survivorship from bark thickness.  In the     !
   ! future survivorship should also depend on fire intensity.                             !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: fire_s_min
   real, dimension(n_pft) :: fire_s_max
   real, dimension(n_pft) :: fire_s_inter
   real, dimension(n_pft) :: fire_s_slope
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Survivorship fraction for plants near felled trees.                               !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: felling_s_ltharv
   real, dimension(n_pft) :: felling_s_gtharv

   !---------------------------------------------------------------------------------------!
   !     Survivorship fractions for plants at skid trails and roads.  Both smaller and     !
   ! larger than the minimum harvesting size.                                              !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: skid_s_ltharv
   real, dimension(n_pft) :: skid_s_gtharv
   !---------------------------------------------------------------------------------------!

   !----- Below this temperature, mortality rapidly increases. ----------------------------!
   real, dimension(n_pft) :: plant_min_temp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   ! Nitrogen and water requirements  -- see "initialize_pft_nitro_params".                !
   !---------------------------------------------------------------------------------------!
   !----- Carbon to Nitrogen ratio, storage pool. -----------------------------------------!
   real :: c2n_storage
   !----- Carbon to Nitrogen ratio, structural stem. --------------------------------------!
   real, dimension(n_pft) :: c2n_stem
   !----- Carbon to Nitrogen ratio, structural stem. --------------------------------------!
   real :: l2n_stem
   !----- Leaf carbon to nitrogen ratio. --------------------------------------------------!
   real, dimension(n_pft) :: c2n_leaf
   !----- Recruit carbon to nitrogen ratio. -----------------------------------------------!
   real, dimension(n_pft) :: c2n_recruit 
   !----- Fraction of structural stem that is assumed to be above ground. -----------------!
   real, dimension(n_pft) :: agf_bs
   !----- Fraction of above-ground wood biomass that is in the branches and twigs. --------!
   real, dimension(n_pft) :: brf_wd
   !----- Supply coefficient for plant nitrogen uptake [m2/kgC_fine_root/day].  -----------!
   real :: plant_N_supply_scale
   !---------------------------------------------------------------------------------------!
   !    Supply coefficient for plant water uptake [m2_ground/kgC_root/sec].                !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: water_conductance  
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   ! Plant hydrodynamics -- see "initialize_pft_hydro_params".                             !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: leaf_water_cap
   !< Leaf hydaulic capacitance [kg H2O/kg biomass/m ]. This variable is assumed as
   !< constants for now

   real, dimension(n_pft) :: wood_water_cap
   !< Wood hydaulic capacitance [kg H2O/kg biomass/m ]. This variable is assumed as
   !< constants for now

   real, dimension(n_pft) :: leaf_water_sat
   !< Leaf water content at saturation (&Psi;=0, rwc=1.) [kg H2O/kg biomass]
   
   real, dimension(n_pft) :: wood_water_sat
   !< Wood water content at saturation (&Psi;=0, rwc=1.) [kg H2O/kg biomass]
   
   real, dimension(n_pft) :: bark_water_sat
   !< Wood water content at saturation (&Psi;=0, rwc=1.) [kg H2O/kg biomass]
   
   real, dimension(n_pft) :: leaf_rwc_min  
   !< Leaf minimum relative water content or leaf residual fraction [-]
   
   real, dimension(n_pft) :: leaf_psi_min
   !< Leaf minimum water potential based on leaf_rwc_min [m]

   real, dimension(n_pft) :: wood_rwc_min  
   !< Sapwood minimum relative water content or Sapwood residual fraction [-]

   real, dimension(n_pft) :: wood_psi_min
   !< Sapwood minimum water potential based on leaf_rwc_min [m]

   real, dimension(n_pft) :: leaf_psi_tlp
   !< Leaf water potential at turgor loss point [m]

   real, dimension(n_pft) :: wood_psi_tlp
   !< Sapwood water potential at turgor loss point [m]

   real, dimension(n_pft) :: leaf_psi_osmotic
   !< Leaf osmotic water potential at saturation [m]

   real, dimension(n_pft) :: wood_psi_osmotic
   !< Sapwood osmotic water potential at saturation [m]

   real, dimension(n_pft) :: leaf_elastic_mod
   !< Leaf bulk elastic modulus [MPa]                    

   real, dimension(n_pft) :: wood_elastic_mod
   !< Sapwood bulk elastic modulus [MPa]                   

   real, dimension(n_pft) :: wood_Kmax     
   !< Maximum hydraulic conductivity of the stem [kg H2O / m / s]       
   
   real, dimension(n_pft) :: wood_Kexp     
   !< Exponent for the hydraulic vulnerability curve of stem conductivity under
   !< the Weibull function 1/(1+(psi/psi50) ** Kexp_stem) [-]

   real, dimension(n_pft) :: wood_psi50
   !< Water potential at which 50% of stem conductivity is lost [m]     

   real, dimension(n_pft) :: vessel_curl_factor
   !< Ratio of actual vessel length (water conducting length) to tree height [-]

   real, dimension(n_pft) :: stoma_lambda
   !< Marginal water use efficiency under well-watered conditions
   !< in the optimization based stomatal model [umol/mol/kPa]
   !< (Katul et al. 2010 Annals of Botany, Manoni et al. 2011 Functional Ecology)

   real, dimension(n_pft) :: stoma_beta
   !< Sensitivity of stoma_lambda to leaf water potential [m-1]

   real, dimension(n_pft) :: stoma_psi_b  
   !< Water potential scaler to modify stomatal conductance under water stress from 
   !< Powell et al. 2017  [ m]
   real, dimension(n_pft) :: stoma_psi_c  
   !< Exponent to modify stomatal conductance under water stress from 
   !< Powell et al. 2017  [ unitless]

   ! Parameters for new drought phenology
   integer, dimension(n_pft) :: high_psi_threshold
   !< Threshold of consecutive wet days to grow new leaves   [# of days]

   integer, dimension(n_pft) :: low_psi_threshold
   !< Threshold of consecutive dry days to grow new leaves   [# of days]

   real, dimension(n_pft) :: leaf_shed_rate
   !< Rate of leaf shedding if low_psi_threshold is crossed  [unitless]

   real, dimension(n_pft) :: leaf_grow_rate
   !< Rate of leaf growing if high_psi_threshold is crossed  [unitless]

   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   ! Allocation and allometry.                                                             !
   !---------------------------------------------------------------------------------------!
   !----- Wood density.  Used only for tropical PFTs and grasses [ g/cm3]. ----------------!
   real   , dimension(n_pft)    :: rho
   !----- Specific Leaf Area (m2leaf/kg_C]. -----------------------------------------------!
   real   , dimension(n_pft)    :: SLA
   !----- Parameters used by the model that predicts SLA based on leaf life span. ---------!
   real   , dimension(n_pft)    :: sla_s0
   real   , dimension(n_pft)    :: sla_s1
   !----- Parameters that control the SLA plasticity. -------------------------------------!
   real   , dimension(n_pft)    :: kplastic_SLA ! Expansion factor for SLA
   real   , dimension(n_pft)    :: LMA_slope    ! Slope for LMA:Height relationship
   !----- Parameter that limits plasticity (Vm, SLA, LTOR). -------------------------------!
   real   , dimension(n_pft)    :: laimax_plastic
   !----- The initialization parameters for SLA:  SLA = sla_pft_init for non-trop PFTs
   real   , dimension(n_pft)    :: sla_pft_init
   !----- Mass ratio between fine root and leaves [kg_fine_roots]/[kg_leaves]. ------------!
   real   , dimension(n_pft)    :: q
   !----- Mass ratio between sapwood and leaves [kg_sapwood]/[kg_leaves]/[m]. -------------!
   real   , dimension(n_pft)    :: qsw
   real   , dimension(n_pft)    :: sapwood_ratio ! AREA ratio
   !----- Mass ratio between bark and leaves [kg_bark]/[kg_leaves]/[m]. -------------------!
   real   , dimension(n_pft)    :: qbark
   !----- Density ratio between bark and wood [(g cm-3)_bark/(g cm-3)_wood]. --------------!
   real   , dimension(n_pft)    :: qrhob
   !----- Specific Root Area (m2root area/kg_C]. ------------------------------------------!
   real   , dimension(n_pft)    :: SRA
   !----- Root vertical profile parameter. Fraction of root biomass below max root depth --!
   real   , dimension(n_pft)    :: root_beta
   !---------------------------------------------------------------------------------------!
   !     DBH-height allometry intercept (m).  Notice that this variable has different      !
   ! meaning between temperate and tropical PFTs.                                          !
   !---------------------------------------------------------------------------------------!
   real   , dimension(n_pft)    :: b1Ht
   !---------------------------------------------------------------------------------------!
   !     DBH-height allometry slope (1/cm).  Notice that this variable has different       !
   ! meaning between temperate and tropical PFTs.                                          !
   !---------------------------------------------------------------------------------------!
   real   , dimension(n_pft)    :: b2Ht
   !---------------------------------------------------------------------------------------!
   !     Reference height for DBH -> height allometry.  They may not be used for tropical  !
   ! PFTs but when they are, they have a different meaning from the temperate allometry.   !
   !---------------------------------------------------------------------------------------!
   real   , dimension(n_pft)    :: hgt_ref
   !----- DBH-stem allometry intercept.  All PFTs. ----------------------------------------!
   real   , dimension(n_pft)    :: b1Bs_small
   !----- DBH-stem allometry slope (dimensionless).  All PFTs. ----------------------------!
   real   , dimension(n_pft)    :: b2Bs_small
   !----- DBH-stem allometry intercept for large DBH cohorts. -----------------------------!
   real   , dimension(n_pft)    :: b1Bs_large
   !----- DBH-stem allometry slope for large DBH cohorts. ---------------------------------!
   real   , dimension(n_pft)    :: b2Bs_large
   !----- stem-DBH allometry intercept.  All PFTs. ----------------------------------------!
   real   , dimension(n_pft)    :: d1DBH_small
   !----- stem-DBH allometry slope (dimensionless).  All PFTs. ----------------------------!
   real   , dimension(n_pft)    :: d2DBH_small
   !----- stem-DBH allometry intercept for large DBH cohorts. -----------------------------!
   real   , dimension(n_pft)    :: d1DBH_large
   !----- stem-DBH allometry slope for large DBH cohorts. ---------------------------------!
   real   , dimension(n_pft)    :: d2DBH_large
   !----- DBH-leaf allometry intercept. All PFTs. -----------------------------------------!
   real   , dimension(n_pft)    :: b1Bl
   !----- DBH-leaf allometry slope. All PFTs. ---------------------------------------------!
   real   , dimension(n_pft)    :: b2Bl
   !----- Leaf-DBH allometry intercept. All PFTs. -----------------------------------------!
   real   , dimension(n_pft)    :: l1DBH
   !----- Leaf-DBH allometry slope. All PFTs. ---------------------------------------------!
   real   , dimension(n_pft)    :: l2DBH
   !----- DBH-crown allometry intercept.  All PFTs. ---------------------------------------!
   real   , dimension(n_pft)    :: b1Ca
   !----- DBH-crown allometry slope.  All PFTs. -------------------------------------------!
   real   , dimension(n_pft)    :: b2Ca
   !----- DBH-WAI allometry intercept. All PFTs. ------------------------------------------!
   real   , dimension(n_pft)    :: b1WAI
   !----- DBH-WAI allometry slope. All PFTs. ----------------------------------------------!
   real   , dimension(n_pft)    :: b2WAI
   !----- DBH-bark thickness slope.  All PFTs. --------------------------------------------!
   real   , dimension(n_pft)    :: b1Xb
   !----- DBH-sapwood thickness slope.  All PFTs. -----------------------------------------!
   real   , dimension(n_pft)    :: b1Xs

   real   , dimension(n_pft)    :: b1SA
   !< DBH-sapwood area allometry intercept
   real   , dimension(n_pft)    :: b2SA
   !< DBH-sapwood area allometry slope
   !----- Minimum DBH attainable by this PFT. ---------------------------------------------!
   real   , dimension(n_pft)    :: min_dbh
   !----- Critical DBH for height/bdead, point in which plants stop growing vertically. ---!
   real   , dimension(n_pft)    :: dbh_crit
   !----- Prescribed DBH for the big leaf model, that allows a reasonable LAI/biomass. ----!
   real   , dimension(n_pft)    :: dbh_bigleaf
   !----- Minimum Bdead attainable by this PFT. -------------------------------------------!
   real   , dimension(n_pft)    :: min_bdead
   !----- Critical Bdead, point in which plants stop growing vertically. ------------------!
   real   , dimension(n_pft)    :: bdead_crit
   !----- Critical Bleaf, maximum allocation to leaves. -----------------------------------!
   real   , dimension(n_pft)    :: bleaf_crit
   !----- Critical balive, maximum allocation to living tissues. --------------------------!
   real   , dimension(n_pft)    :: balive_crit
   !----- Critical balive+bdead ("Everything but storage"). -------------------------------!
   real   , dimension(n_pft)    :: bevery_crit
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   ! Leaf habit and physical properties.                                                   !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   ! Phenology indicates the leaf habit regarding phenology:                               !
   ! 0. Evergreen coniferous;                                                              !
   ! 1. Drought deciduous;                                                                 !
   ! 2. Cold deciduous;                                                                    !
   ! 3. Light controlled;                                                                  !
   ! 4. Drought deciduous - based on 10day average.                                        !
   !---------------------------------------------------------------------------------------!
   integer, dimension(n_pft) :: phenology 


   !----- Leaf width [m], which is used to compute the leaf boundary layer conductance. ---!
   real, dimension(n_pft) :: leaf_width
   !---------------------------------------------------------------------------------------!


   !----- This is the characteristic branch diameter [m]. ---------------------------------!
   real, dimension(n_pft) :: branch_diam
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Parameters to find the crown length, which will be used to find the height of the !
   ! bottom of the crown.                                                                  !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: b1Cl
   real, dimension(n_pft) :: b2Cl

   !---------------------------------------------------------------------------------------!
   !     Parameters to find the volume and the root depth.                                 !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: b1Vol
   real, dimension(n_pft) :: b2Vol
   real, dimension(n_pft) :: b1Rd
   real, dimension(n_pft) :: b2Rd

   !---------------------------------------------------------------------------------------!
   !     Parameters to find the effective functional root depth, following B18.  This is   !
   ! used only by tropical trees and only when variable use_efrd_trtree is set to .true.   !
   !                                                                                       !
   ! Reference:                                                                            !
   !                                                                                       !
   ! Brum M, Vadeboncoeur MA, Ivanov V, Asbjornsen H, Saleska S, Alves LF, Penha D,        !
   !    Dias JD, Aragao LEOC, Barros F, Bittencourt P, Pereira L, Oliveira RS, 2018.       !
   !    Hydrological niche segregation defines forest structure and drought                !
   !    tolerance strategies in a seasonal Amazonian forest. J. Ecol., in press.           !
   !    doi:10.1111/1365-2745.13022 (B18).                                                 !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: d18O_ref
   real, dimension(n_pft) :: b1d18O
   real, dimension(n_pft) :: b2d18O
   real, dimension(n_pft) :: b1Efrd
   real, dimension(n_pft) :: b2Efrd
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   Parameters used for vegetation heat capacity calculation.                           !
   !---------------------------------------------------------------------------------------!
   !----- Specific heat capacity of dry leaf biomass [J/kg/K]. ----------------------------!
   real, dimension(n_pft) :: c_grn_leaf_dry
   !----- Specific heat capacity of dry non-green wood biomass [J/kg/K]. ------------------!
   real, dimension(n_pft) :: c_ngrn_wood_dry
   !----- Specific heat capacity of dry non-green bark biomass [J/kg/K]. ------------------!
   real, dimension(n_pft) :: c_ngrn_bark_dry
   !-----  Correction-term for energy storage in wood-water bond. -------------------------!
   real, dimension(n_pft) :: delta_c_wood
   real, dimension(n_pft) :: delta_c_bark
   !----- Net specific heat capacity of leaves, sapwood, heartwood, and bark. -------------!
   real, dimension(n_pft) :: cleaf
   real, dimension(n_pft) :: csapw
   real, dimension(n_pft) :: cdead
   real, dimension(n_pft) :: cbark
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Reproduction and recruitment.                                                     !
   !---------------------------------------------------------------------------------------!
   !----- Initial plant density in a near-bare-ground run [plant/m2]. ---------------------!
   real   , dimension(n_pft)    :: init_density
   !----- Initial maximum LAI in a near-bare-ground run [m2/m2] - Big leaf only. ----------!
   real   , dimension(n_pft)    :: init_laimax
   !----- Minimum height of an individual [m]. --------------------------------------------!
   real   , dimension(n_pft)    :: hgt_min
   !----- Maximum height of an individual [m]. --------------------------------------------!
   real   , dimension(n_pft)    :: hgt_max
   !----- Minimum biomass density [kgC/m2] required to form a new recruit. ----------------!
   real   , dimension(n_pft) :: min_recruit_size
   !----- Amount of biomass [kgC] in one tree, used for 'big-leaf' ED. --------------------!
   real   , dimension(n_pft) :: one_plant_c
   !---------------------------------------------------------------------------------------!
   !    Fraction of (positive) carbon balance devoted to storage (unwise to set this to    !
   ! anything other than zero unless storage turnover rate is adjusted accordingly).       !
   !---------------------------------------------------------------------------------------!
   real   , dimension(n_pft) :: st_fract
   !----- Reproduction allocation function (true = bang; false = asymptote). --------------!
   logical, dimension(n_pft) :: r_bang
   !----- Fraction of (positive) carbon balance devoted to reproduction (bang). -----------!
   real   , dimension(n_pft) :: r_fract
   !----- Curvature term for asymptote. ---------------------------------------------------!
   real   , dimension(n_pft) :: r_cv50
   !----- External input of seeds [kgC/m2/year]. ------------------------------------------!
   real   , dimension(n_pft) :: seed_rain
   !----- Fraction of seed dispersal that is gridcell-wide. -------------------------------!
   real   , dimension(n_pft) :: nonlocal_dispersal !  
   !----- Minimum height plants need to attain before allocating to reproduction. ---------!
   real   , dimension(n_pft) :: repro_min_h
   !----- Minimum DBH plants need to attain before allocating to reproduction. ------------!
   real   , dimension(n_pft) :: repro_min_dbh
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     Variables to control carbon stocks during initialisation.                         !
   !---------------------------------------------------------------------------------------!
   !------ Storage biomass, relative to on-allometry living biomass. ----------------------!
   real, dimension(n_pft) :: f_bstorage_init
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     The following variables control the cohort existence/termination.                 !
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !    Minimum size (measured as biomass of living and structural tissues) allowed in a   !
   ! cohort.  Cohorts with less biomass than this are going to be terminated.              !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: min_cohort_size
   !---------------------------------------------------------------------------------------!
   !    The following variable is the absolute minimum cohort population that a cohort can !
   ! have.  This should be used only to avoid nplant=0, but IMPORTANT: this will lead to a !
   ! ridiculously small cohort almost guaranteed to be extinct and SHOULD BE USED ONLY IF  !
   ! THE AIM IS TO ELIMINATE THE COHORT.                                                   !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: negligible_nplant
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     The following varible will be used to "turn off" the biophysics for extremely low !
   ! biomass cohorts, that otherwise could shrink the time step in case they were solved.  !
   ! We use heat capacity rather than leaf/wood area index as the threshold because the    !
   ! heat capacity what will control the time step.  Also, in case of leaves, the bio-     !
   ! physics "turn off" will kill the cohorts, because they won't be able to do photo-     !
   ! synthesis.                                                                            !
   !=======================================================================================!
   !=======================================================================================!
   real, dimension(n_pft) :: veg_hcap_min
   !=======================================================================================!
   !=======================================================================================!


   !---------------------------------------------------------------------------------------!
   !     Labile fraction of different tissues.  This is the fraction of biomass that goes  !
   ! to the fast soil pools as opposed to the structural soil carbon.                      !
   ! MLO.  Migrated the parameters from decomp_coms to here because these are              !
   ! PFT-dependent parameters.  Also split the fractions for "leaf" and "stem".            !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_pft) :: f_labile_leaf !  Leaf and fire root
   real, dimension(n_pft) :: f_labile_stem !  Sapwood, bark, heartwood
   !---------------------------------------------------------------------------------------!






   !=======================================================================================!
   !=======================================================================================!
   !     Liana-specific parameters.                                                        !
   !=======================================================================================!
   !=======================================================================================!
   real :: h_edge          !< maximum height advantage for lianas
   real :: liana_dbh_crit  !< liana specific critical dbh
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !     Look-up table used to convert total woody biomass to DBH.                         !
   !---------------------------------------------------------------------------------------!
   integer                                   :: nbt_lut
   real(kind=4), dimension(:,:), allocatable :: dbh_lut
   real(kind=4), dimension(:,:), allocatable :: bleaf_lut
   real(kind=4), dimension(:,:), allocatable :: bdead_lut
   real(kind=4), dimension(:,:), allocatable :: balive_lut
   real(kind=4), dimension(:,:), allocatable :: bevery_lut
   logical     , dimension(:)  , allocatable :: le_mask_lut
   logical     , dimension(:)  , allocatable :: ge_mask_lut
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     The following variable is an identifier for the PFT, for some of the debugging    !
   ! output.                                                                               !
   !---------------------------------------------------------------------------------------!
   character(len=16), dimension(n_pft) :: pft_name16
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Structure containing information needed for recruits.                              !
   !---------------------------------------------------------------------------------------!
   type recruittype
      integer             :: pft
      integer             :: krdepth
      integer             :: phenology_status
      real                :: leaf_temp
      real                :: wood_temp
      real                :: leaf_temp_pv
      real                :: wood_temp_pv
      real                :: leaf_vpdef
      real                :: hite
      real                :: dbh
      real                :: bdeada
      real                :: bdeadb
      real                :: bleaf
      real                :: broot
      real                :: bsapwooda
      real                :: bsapwoodb
      real                :: bbarka
      real                :: bbarkb
      real                :: balive
      real                :: paw_avg
      real                :: elongf
      real                :: bstorage
      real                :: nplant
      real, dimension(13) :: cb
      real, dimension(13) :: cb_lightmax
      real, dimension(13) :: cb_moistmax
      real, dimension(13) :: cb_mlmax
      real                :: cbr_bar
   end type recruittype
   !=======================================================================================!
   !=======================================================================================!
   
   
   contains



   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine simply resets a recruittype structure.                             !
   !---------------------------------------------------------------------------------------!
   subroutine zero_recruit(maxp,recruit)
      implicit none
      !----- Argument. --------------------------------------------------------------------!
      integer                           , intent(in)  :: maxp
      type(recruittype), dimension(maxp), intent(out) :: recruit
      !----- Local variable. --------------------------------------------------------------!
      integer                                         :: p
      integer                                         :: imon
      !------------------------------------------------------------------------------------!

      do p=1,maxp
         recruit(p)%pft                  = 0
         recruit(p)%krdepth              = 0
         recruit(p)%phenology_status     = 0
         recruit(p)%leaf_temp            = 0.
         recruit(p)%wood_temp            = 0.
         recruit(p)%leaf_temp_pv         = 0.
         recruit(p)%wood_temp_pv         = 0.
         recruit(p)%leaf_vpdef           = 0.
         recruit(p)%hite                 = 0.
         recruit(p)%dbh                  = 0.
         recruit(p)%bdeada               = 0.
         recruit(p)%bdeadb               = 0.
         recruit(p)%bleaf                = 0.
         recruit(p)%broot                = 0.
         recruit(p)%bsapwooda            = 0.
         recruit(p)%bsapwoodb            = 0.
         recruit(p)%bbarka               = 0.
         recruit(p)%bbarkb               = 0.
         recruit(p)%balive               = 0.
         recruit(p)%paw_avg              = 0.
         recruit(p)%elongf               = 0.
         recruit(p)%bstorage             = 0.
         recruit(p)%nplant               = 0.
         do imon=1,13
            recruit(p)%cb         (imon) = 0.
            recruit(p)%cb_lightmax(imon) = 0.
            recruit(p)%cb_moistmax(imon) = 0.
            recruit(p)%cb_mlmax   (imon) = 0.
         end do
         recruit(p)%cbr_bar              = 0.
      end do

      return
   end subroutine zero_recruit
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine simply copies a recruittype structure.                             !
   !---------------------------------------------------------------------------------------!
   subroutine copy_recruit(recsource,rectarget)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(recruittype), intent(in)  :: recsource
      type(recruittype), intent(out) :: rectarget
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: imon
      !------------------------------------------------------------------------------------!

      rectarget%pft              = recsource%pft
      rectarget%krdepth          = recsource%krdepth
      rectarget%phenology_status = recsource%phenology_status
      rectarget%leaf_temp        = recsource%leaf_temp
      rectarget%wood_temp        = recsource%wood_temp
      rectarget%leaf_temp_pv     = recsource%leaf_temp_pv
      rectarget%wood_temp_pv     = recsource%wood_temp_pv
      rectarget%leaf_vpdef       = recsource%leaf_vpdef
      rectarget%hite             = recsource%hite
      rectarget%dbh              = recsource%dbh
      rectarget%bdeada           = recsource%bdeada
      rectarget%bdeadb           = recsource%bdeadb
      rectarget%bleaf            = recsource%bleaf
      rectarget%broot            = recsource%broot
      rectarget%bsapwooda        = recsource%bsapwooda
      rectarget%bsapwoodb        = recsource%bsapwoodb
      rectarget%bbarka           = recsource%bbarka
      rectarget%bbarkb           = recsource%bbarkb
      rectarget%balive           = recsource%balive
      rectarget%paw_avg          = recsource%paw_avg
      rectarget%elongf           = recsource%elongf
      rectarget%bstorage         = recsource%bstorage
      rectarget%nplant           = recsource%nplant

      do imon=1,13
         rectarget%cb         (imon) = recsource%cb         (imon)
         rectarget%cb_lightmax(imon) = recsource%cb_lightmax(imon)
         rectarget%cb_moistmax(imon) = recsource%cb_moistmax(imon)
         rectarget%cb_mlmax   (imon) = recsource%cb_mlmax   (imon)
      end do

      rectarget%cbr_bar          = recsource%cbr_bar

      return
   end subroutine copy_recruit
   !=======================================================================================!
   !=======================================================================================!
end module pft_coms
!==========================================================================================!
!==========================================================================================!

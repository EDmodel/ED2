!==========================================================================================!
!==========================================================================================!
!     This is the main loader of ecosystem parameters.  Since some compilers do not under- !
! stand the assignment in the modules when the variable is not a constant (parameter),     !
! this is the safest way to guarantee it will read something (not to mention that makes    !
! compilation much faster when you want to test the sensitivity of one number).            !
!------------------------------------------------------------------------------------------!
subroutine load_ed_ecosystem_params()

   use ed_max_dims , only : n_pft               ! ! intent(in)
   use pft_coms    , only : include_these_pft   & ! intent(in)
                          , pasture_stock       & ! intent(in)
                          , agri_stock          & ! intent(in)
                          , plantation_stock    & ! intent(in)
                          , pft_name16          & ! intent(out)
                          , is_tropical         & ! intent(out)
                          , is_savannah         & ! intent(out)
                          , is_conifer          & ! intent(out)
                          , is_grass            & ! intent(out)
                          , is_liana            & ! intent(out)
                          , include_pft         & ! intent(out)
                          , include_pft_pt      & ! intent(out)
                          , include_pft_ag      & ! intent(out)
                          , include_pft_fp      ! ! intent(out)
   use disturb_coms, only : ianth_disturb       ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer :: p
   !---------------------------------------------------------------------------------------!



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
   !   9 | Early hardwood                  |    no |    no |       no |       no |      no !
   !  10 | Mid hardwood                    |    no |    no |       no |       no |      no !
   !  11 | Late hardwood                   |    no |    no |       no |       no |      no !
   !  12 | Early savannah                  |    no |    no |      yes |       no |      no !
   !  13 | Mid savannah                    |    no |    no |      yes |       no |      no !
   !  14 | Late savannah                   |    no |    no |      yes |       no |      no !
   !  15 | Araucaria                       |    no |    no |      yes |       no |     yes !
   !  16 | Tropical C3 grass               |   yes |    no |      yes |       no |      no !
   !  17 | Liana                           |    no |   yes |      yes |       no |     yes !
   !---------------------------------------------------------------------------------------!


   !----- Name the PFTs (no spaces, please). ----------------------------------------------!
   pft_name16( 1) = 'C4_grass        '
   pft_name16( 2) = 'Early_tropical  '
   pft_name16( 3) = 'Mid_tropical    '
   pft_name16( 4) = 'Late_tropical   '
   pft_name16( 5) = 'C3_grass        '
   pft_name16( 6) = 'North_pine      '
   pft_name16( 7) = 'South_pine      '
   pft_name16( 8) = 'Late_conifer    '
   pft_name16( 9) = 'Early_hardwood  '
   pft_name16(10) = 'Mid_hardwood    '
   pft_name16(11) = 'Late_hardwood   '
   pft_name16(12) = 'Early_savannah  '
   pft_name16(13) = 'Mid_savannah    '
   pft_name16(14) = 'Late_savannah   '
   pft_name16(15) = 'Araucaria       '
   pft_name16(16) = 'Subtrop_C3_grass'
   pft_name16(17) = 'Liana           '
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    This flag should be used to define whether the plant is tropical/subtropical or    !
   ! not.                                                                                  !
   !---------------------------------------------------------------------------------------!
   is_tropical(1:4)   = .true.
   is_tropical(5:11)  = .false.
   is_tropical(12:14) = .true.
   is_tropical(15)    = .true.
   is_tropical(16)    = .true.
   is_tropical(17)    = .true.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------! 
   !    This flag should be used to define whether the plant is a savannah tree or not.    !
   ! Currently this is used only to define fire adaptation, i.e. bark thickness            !
   ! parameters.  Grasses are not set as savannah because they have no bark                !
   !---------------------------------------------------------------------------------------! 
   is_savannah(1:11)  = .false.
   is_savannah(12:14) = .true.
   is_savannah(15:17) = .false.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !           This flag is used to define whether the plant is a liana or not             !
   !---------------------------------------------------------------------------------------!
   is_liana(1:16)  = .false.
   is_liana(17)    = .true.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------! 
   !    This flag should be used to define whether the plant is conifer or flowering.      !
   !---------------------------------------------------------------------------------------! 
   is_conifer(1:5)   = .false.
   is_conifer(6:8)   = .true.
   is_conifer(9:14)  = .false.
   is_conifer(15)    = .true.
   is_conifer(16:17) = .false.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------! 
   !    This flag should be used to define whether the plant is tree or grass.             !
   !---------------------------------------------------------------------------------------! 
   is_grass(1)     = .true.
   is_grass(2:4)   = .false.
   is_grass(5)     = .true.
   is_grass(6:15)  = .false.
   is_grass(16)    = .true.
   is_grass(17)    = .false.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Include_pft: flag specifying to whether you want to include a plant functional     !
   !                 type (T) or whether you want it excluded (F) from the simulation.     !
   !---------------------------------------------------------------------------------------!
   include_pft    = .false.
   do p=1,n_pft
      if (include_these_pft(p) >  0 .and. include_these_pft(p) <= n_pft) then
         include_pft(include_these_pft(p)) = .true.
      end if
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Only the PFTs listed in pasture_stock are allowed in pasture patches.  For the    !
   ! time being this means a single PFT, but it could change in the future.                !
   !---------------------------------------------------------------------------------------!
   include_pft_pt                = .false.
   include_pft_pt(pasture_stock) = is_grass(pasture_stock) .and. include_pft(pasture_stock)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Only the PFTs listed in agri_stock are allowed in agriculture patches.  For the   !
   ! time being this means a single PFT, but it could change in the future.                !
   !---------------------------------------------------------------------------------------!
   include_pft_ag             = .false.
   include_pft_ag(agri_stock) = is_grass(agri_stock) .and. include_pft(agri_stock)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Only the PFTs listed in plantation_stock are allowed in agriculture patches.  For !
   ! the time being this means a single PFT, but it could change in the future.            !
   !---------------------------------------------------------------------------------------!
   include_pft_fp                   = .false.
   include_pft_fp(plantation_stock) = ( .not. is_grass(plantation_stock) ) .and.           &
                                      include_pft(plantation_stock)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Warn the user in case the PFT choice for pasture, agriculture or forest          !
   ! plantation was inconsistent.                                                          !
   !---------------------------------------------------------------------------------------!
   if (count(include_pft_pt) == 0 .and. ianth_disturb /= 0) then
      call warning ('PFT defined in pasture_stock is not included in include_these_pft,'// &
                    ' your cattle will starve and your ranch will not be profitable...'    &
                   ,'load_ecosystem_params','ed_params.f90')
   end if
   if (count(include_pft_ag) == 0 .and. ianth_disturb /= 0) then
      call warning ('PFT defined in agri_stock is not included in include_these_pft,'//    &
         ' your croplands will be barren and not very profitable...'            &
         ,'load_ecosystem_params','ed_params.f90')
   end if
   if (count(include_pft_fp) == 0 .and. ianth_disturb /= 0) then
      call warning ('PFT defined in plantation_stock is not listed in include_these_pft,'//&
         ' your forest plantation will be barren and not very profitable ...'   &
         ,'load_ecosystem_params','ed_params.f90')
   end if
   !---------------------------------------------------------------------------------------!




   !----- Load several parameters ---------------------------------------------------------!
   call init_decomp_params()
   call init_ff_coms()
   call init_disturb_params()
   call init_physiology_params()
   call init_met_params()
   call init_lapse_params()
   call init_hydro_coms()
   call init_soil_coms()
   call init_phen_coms()
   call init_ed_misc_coms()
   call init_hrzshade_params()
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Assign many PFT-dependent parameters.  Here the order may matter, so think twice  !
   ! before changing the order.                                                            !
   !---------------------------------------------------------------------------------------!
   !----- Allometry and some plant traits. ------------------------------------------------!
   call init_pft_alloc_params()
   !----- Photosynthesis and leaf respiration. --------------------------------------------!
   call init_pft_photo_params()
   !----- Root and heterotrophic respiration. ---------------------------------------------!
   call init_pft_resp_params()
   !----- Mortality. ----------------------------------------------------------------------!
   call init_pft_mort_params()
   !----- Nitrogen. -----------------------------------------------------------------------!
   call init_pft_nitro_params()
   !----- Hydraulics ----------------------------------------------------------------------!
   call init_pft_hydro_params()
   !----- Specific heat -------------------------------------------------------------------!
   call init_pft_spheat_params()
   !----- Phenology leaf properties. ------------------------------------------------------!
   call init_pft_phen_params()
   !----- Reproduction. -------------------------------------------------------------------!
   call init_pft_repro_params()
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Assign canopy properties.  This must be done after defining the PFT stuff.        !
   ! Again, check the routines and think twice before changing the order.                  !
   !---------------------------------------------------------------------------------------!
   !----- Canopy turbulence and aerodynamic resistance. -----------------------------------!
   call init_can_air_params()
   !----- Canopy radiation properties. ----------------------------------------------------!
   call init_can_rad_params()
   !----- Canopy splitting into height layers. --------------------------------------------!
   call init_can_lyr_params()
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     This should be always the last one, since it depends on variables assigned in the !
   ! previous init_????_params.                                                            !
   !---------------------------------------------------------------------------------------!
   call init_dt_thermo_params()
   !---------------------------------------------------------------------------------------!

   return
end subroutine load_ed_ecosystem_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine initialises parameters associated with necromass and soil C (and     !
! soil N) dynamics.                                                                        !
!------------------------------------------------------------------------------------------!
subroutine init_decomp_params()
   use soil_coms   , only : slz                         ! ! intent(in)
   use grid_coms   , only : nzg                         ! ! intent(in)
   use consts_coms , only : yr_day                      & ! intent(in)
                          , t00                         ! ! intent(in)
   use decomp_coms , only : decomp_scheme               & ! intent(in)
                          , resp_opt_water              & ! intent(out)
                          , resp_water_below_opt        & ! intent(out)
                          , resp_water_above_opt        & ! intent(out)
                          , resp_temperature_increase   & ! intent(out)
                          , N_immobil_supply_scale      & ! intent(out)
                          , e_lignin                    & ! intent(out)
                          , r_fsc                       & ! intent(out)
                          , r_stsc_l                    & ! intent(out)
                          , r_stsc_o                    & ! intent(out)
                          , r_msc_int                   & ! intent(out)
                          , r_msc_slp                   & ! intent(out)
                          , r_ssc                       & ! intent(out)
                          , r_psc                       & ! intent(out)
                          , decay_rate_stsc             & ! intent(out)
                          , decay_rate_fsc              & ! intent(out)
                          , decay_rate_msc              & ! intent(out)
                          , decay_rate_ssc              & ! intent(out)
                          , decay_rate_psc              & ! intent(out)
                          , fx_msc_psc_int              & ! intent(out)
                          , fx_msc_psc_slp              & ! intent(out)
                          , fx_ssc_psc_int              & ! intent(out)
                          , fx_ssc_psc_slp              & ! intent(out)
                          , rh_lloyd_1                  & ! intent(out)
                          , rh_lloyd_2                  & ! intent(out)
                          , rh_lloyd_3                  & ! intent(out)
                          , rh_decay_low                & ! intent(out)
                          , rh_decay_high               & ! intent(out)
                          , rh_low_temp                 & ! intent(out)
                          , rh_high_temp                & ! intent(out)
                          , rh_decay_dry                & ! intent(out)
                          , rh_decay_wet                & ! intent(out)
                          , rh_dry_smoist               & ! intent(out)
                          , rh_wet_smoist               & ! intent(out)
                          , rh_active_depth             & ! intent(out)
                          , rh_moyano12_a0              & ! intent(out)
                          , rh_moyano12_a1              & ! intent(out)
                          , rh_moyano12_a2              & ! intent(out)
                          , rh0                         & ! intent(out)
                          , rh_q10                      & ! intent(out)
                          , rh_p_smoist                 & ! intent(out)
                          , rh_p_oxygen                 & ! intent(out)
                          , agf_fsc                     & ! intent(out)
                          , agf_stsc                    & ! intent(out)
                          , f0_msc                      & ! intent(out)
                          , f0_psc                      & ! intent(out)
                          , nbg_nlim_fsc                & ! intent(out)
                          , nbg_nlim_stsc               & ! intent(out)
                          , nbg_nlim_ssc                ! ! intent(out)
   implicit none
   !---------------------------------------------------------------------------------------!

   
   resp_opt_water            = 0.8938
   resp_water_below_opt      = 5.0786
   resp_water_above_opt      = 4.5139
   resp_temperature_increase = 0.0757
   N_immobil_supply_scale    = 40.0 / yr_day


   !---------------------------------------------------------------------------------------!
   !  CENTURY parameter that controls the effect of lignin on structural decomposition.    !
   !---------------------------------------------------------------------------------------!
   e_lignin = 3.0
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Fraction of decay that is lost as CO2.                                            !
   !---------------------------------------------------------------------------------------!
   select case (decomp_scheme)
   case (5)
      r_fsc           = 0.55
      r_stsc_l        = 0.20
      r_stsc_o        = 0.50
      r_msc_int       = 0.60
      r_msc_slp       = 0.17
      r_ssc           = 0.55
      r_psc           = 0.55
   case default
      r_fsc           = 1.0
      r_stsc_l        = 0.3
      r_stsc_o        = 0.3
      r_msc_int       = 0.0
      r_msc_slp       = 0.0
      r_ssc           = 1.0
      r_psc           = 1.0
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Decay rates for necromass and soil carbon pools.                                     !
   !---------------------------------------------------------------------------------------!
   select case (decomp_scheme)
   case (0,1)
      !------------------------------------------------------------------------------------!
      ! MLO.  After talking to Paul, it seems the decay rate for the slow carbon pool is   !
      !       artificially high for when nitrogen limitation is turned on.  If it is       !
      !       turned off, however, then the slow carbon will disappear very quickly.  I    !
      !       don't want to mess other people's results, so I will change the rate only    !
      !       when DECOMP_SCHEME is 2 (see below).  I think this should be applied to all  !
      !       schemes, but I will let the users of these schemes to decide.                !
      !------------------------------------------------------------------------------------!
      decay_rate_fsc  =  11.0 / yr_day    ! former K2
      decay_rate_stsc =   4.5 / yr_day    ! former K1
      decay_rate_msc  =   0.0 / yr_day    ! Inexistent.  Microbes are not solved here.
      decay_rate_ssc  = 100.2 / yr_day    ! former K3
      decay_rate_psc  =   0.0 / yr_day    ! Inexistent.  Passive is not solved here.
      !------------------------------------------------------------------------------------!
   case (2)
      !----- Update rate for slow (so slow is slower than fast...). -----------------------!
      decay_rate_fsc  =  11.0 / yr_day    ! former K2
      decay_rate_stsc =   4.5 / yr_day    ! former K1
      decay_rate_msc  =   0.0 / yr_day    ! Inexistent.  Microbes are not solved here.
      decay_rate_ssc  =   0.2 / yr_day    ! former K3
      decay_rate_psc  =   0.0 / yr_day    ! Inexistent.  Passive is not solved here.
      !------------------------------------------------------------------------------------!
   case (5)
      !------------------------------------------------------------------------------------!
      !    CENTURY model, closer to the implementation by B98.  Note that the original ED2 !
      ! has fewer carbon pools, and the names do not exactly match B98 to be consistent    !
      ! with the original ED-2.0 implementation. The table below has the current           !
      ! translation.                                                                       !
      !                                                                                    !
      !      |----------------------+---------------------------+--------------------|     !
      !      |  B98                 |    DECOMP_SCHEME = 0      | DECOMP_SCHEME = 5  |     !
      !      |----------------------+---------------------------+--------------------|     !
      !      | Metabolic Litter     | Fast                      | Fast               |     !
      !      | Structural Litter    | Structural                | Structural         |     !
      !      | Microbes (fast SOM)  | ------ (lumped with Fast) | Microbial          |     !
      !      | Slow SOM             | Slow                      | Slow (Humified)    |     !
      !      | Passive SOM          | ------ (lumped with Slow) | Passive            |     !
      !      |----------------------+---------------------------+--------------------|     !
      !                                                                                    !
      ! Note that ED-1.0 did have Passive SOM but it was literally passive (i.e. no flux   !
      ! in or out).  Also note that we do not split the microbial (fast SOM) between       !
      ! above- and below-ground, but we do distinguish the metabolic and structural litter !
      ! because of fires (which only burn the above-ground fraction).  The current         !
      ! parameters are similar to CENTURY and are related to the half-life time shown in   !
      ! B98's Fig. 3.  Note, however, the typo in their definition in their caption, it    !
      ! should be ln(2), not log2(e).  The actual numbers were obtained from M96.  Because !
      ! most studies report soil carbon for the top 30 cm (as opposed to 20cm), we         !
      ! decreased the default decay rates for microbial, humified and passive pools by     !
      ! 15%, following the suggestion by M96, except for the passive pool, because the     !
      ! typical rate seems to be higher than M96 (e.g. see B98's Fig. 3).  For the passive !
      ! pool, we assumed a half-life of about 70 years.                                    !
      !                                                                                    !
      ! References:                                                                        !
      !                                                                                    !
      ! Bolker BM, Pacala SW, and Parton WJ. 1998. Linear analysis of soil decomposition:  !
      !    insights from the CENTURY model. Ecol. Appl., 8(2):425-439.                     !
      !    doi:10.1890/1051- 0761(1998)008[0425:LAOSDI]2.0.CO;2 (B98).                     !
      !                                                                                    !
      ! Metherell AK, Harding LA, Cole CV, Parton WJ. 1996. CENTURY Soil Organic Matter    !
      !    model environment.  Agroecosystem Version 4.0. Great Plains System Research     !
      !    Unit. Techinical Report No 4. USDA-ARS, Fort Collins CO.                        !
      !    https://www2.nrel.colostate.edu/projects/century/MANUAL/html_manual/man96.html. !
      !    Accessed on 9 Aug 2018 (M96).                                                   !
      !                                                                                    !
      !------------------------------------------------------------------------------------!
      decay_rate_fsc  =  12.0  / yr_day ! 16.7 / yr_day
      decay_rate_stsc =  1.5   / yr_day
      decay_rate_msc  =  6.0   / yr_day
      decay_rate_ssc  =  0.15  / yr_day
      decay_rate_psc  =  0.012 / yr_day 
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      When soil decays, a fraction is lost as CO2 (heterotrophic respiration).  The    !
   ! remaining decayed carbon is partitioned between pools.  Parameters define the         !
   ! fraction that goes to passive carbon.  Other fractions are internally defined to      !
   ! ensure that all fluxes are conserved.                                                 !
   !---------------------------------------------------------------------------------------!
   fx_msc_psc_int = 0.003
   fx_msc_psc_slp = 0.032
   fx_ssc_psc_int = 0.003
   fx_ssc_psc_slp = 0.009
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Parameters used for the Lloyd and Taylor (1994) temperature dependence.          !
   !---------------------------------------------------------------------------------------!
   rh_lloyd_1 = 308.56
   rh_lloyd_2 = 1./56.02
   rh_lloyd_3 = 227.15
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The following variables are used when DECOMP_SCHEME is 3 or 4. (based on M12).    !
   !                                                                                       !
   !  Moyano FE, Vasilyeva N, Bouckaert L, Cook F, Craine J, Curiel Yuste J, Don A,        !
   !     Epron D, Formanek P, Franzluebbers A et al. 2012. The moisture response of soil   !
   !     heterotrophic respiration: interaction with soil properties. Biogeosciences, 9:   !
   !     1173-1182. doi:10.5194/bg-9-1173-2012 (M12).                                      !
   !---------------------------------------------------------------------------------------!
   rh_moyano12_a0 = -0.3195897
   rh_moyano12_a1 =  4.0893
   rh_moyano12_a2 = -3.1681
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Parameters used for the ED-1.0/CENTURY based functions of temperature and soil   !
   ! moisture (DECOMP_SCHEME =2).                                                          !
   !---------------------------------------------------------------------------------------!
   rh_decay_low   = 0.24
   rh_decay_high  = 0.60
   rh_low_temp    = 18.0 + t00
   rh_high_temp   = 45.0 + t00
   rh_decay_dry   = 12.0 ! 18.0
   rh_decay_wet   = 36.0 ! 36.0
   rh_dry_smoist  = 0.48 ! 0.36
   rh_wet_smoist  = 0.98 ! 0.96
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Parameters used for the new temperature scaling of , based on CLM (K13), but     !
   ! with a power factor for water limitation, and a simplified approach for anoxic        !
   ! environment (used when DECOMP_SCHEME=5).                                              !
   !                                                                                       !
   !  Koven CD, Riley WJ, Subin ZM, Tang JY, Torn MS, Collins WD, Bonan GB, Lawrence DM,   !
   !     Swenson SC. 2013. The effect of vertically resolved soil biogeochemistry and      !
   !     alternate soil C and N models on C dynamics of CLM4. Biogeosciences, 10:          !
   !     7109-7131. doi:10.5194/bg-10-7109-2013 (K13).                                     !
   !---------------------------------------------------------------------------------------!
   rh0          = 0.700 ! 0.701 ! 0.425
   rh_q10       = 1.500 ! 1.500 ! 1.893
   rh_p_smoist  = 1.600 ! 0.836 ! 0.606
   rh_p_oxygen  = 0.450 ! 0.404 ! 0.164
   !---------------------------------------------------------------------------------------!




   !----- Determine the top layer to consider for heterotrophic respiration. --------------!
   select case (decomp_scheme)
   case (2)
      !---- Use 20cm for back-compatibility with ED-2.2. ----------------------------------!
      rh_active_depth = -0.20
      !------------------------------------------------------------------------------------!
   case (5)
      !------------------------------------------------------------------------------------!
      !     To be consistent with most measurements, we set the active soil layer to the   !
      ! top 30 cm.  This is somewhat deeper than the implicit depth in the CENTURY model   !
      ! (20 cm), but the decay rates have been adjusted to account for deeper soil.  Most  !
      ! data in the tropics are provided for the top 30 cm, hence our choice.              !
      !------------------------------------------------------------------------------------!
      rh_active_depth = -0.30
      !------------------------------------------------------------------------------------!
   case default
      !----- Use the top soil layer (back-compatibility with ED-1). -----------------------!
      rh_active_depth = slz(nzg)
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     These variables are used for initialisation of above- and below-ground necromass  !
   ! pools.                                                                                !
   !---------------------------------------------------------------------------------------!
   agf_fsc            = 0.5
   agf_stsc           = 0.7
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     These variables define the partition of SOM into microbial, passive, and          !
   ! humified.  These fractions are used only during initialisation (NBG or PSS/CSS files) !
   ! and only when DECOMP_SCHEME is 5.  In this case, the fraction of soil C that is       !
   ! humified ("slow") will be 1.0 - f_msc - f_psc.   The default fractions are based on   !
   ! the CENTURY manual, but can be changed through XML.                                   !
   !                                                                                       !
   ! https://www2.nrel.colostate.edu/projects/century/MANUAL/html_manual/man96.html        !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   f0_msc   = 0.03
   f0_psc   = 0.35
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters to initialise soil carbon pools for near bare ground simulations when  !
   ! N limitation is turned on.  These used to be hardcoded in ed_nbg_init.f90, and hard-  !
   ! coded parameters are deprecated.                                                      !
   !---------------------------------------------------------------------------------------!
   nbg_nlim_fsc  = 0.2
   nbg_nlim_stsc = 10.0 ! Isn't this too high for NBG?
   nbg_nlim_ssc  = 0.01
   !---------------------------------------------------------------------------------------!

   return

end subroutine init_decomp_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine assigns the fusion and splitting parameters.                         !
!------------------------------------------------------------------------------------------!
subroutine init_ff_coms
   use fusion_fission_coms, only : ifusion                   & ! intent(in)
                                 , ff_nhgt                   & ! intent(out)
                                 , niter_patfus              & ! intent(out)
                                 , fusetol                   & ! intent(out), old_fusion
                                 , fusetol_h                 & ! intent(out), old_fusion
                                 , lai_fuse_tol              & ! intent(out), old_fusion
                                 , coh_tolerance_max         & ! intent(out), old_fusion
                                 , dark_cumlai_min           & ! intent(out), old_fusion
                                 , dark_cumlai_max           & ! intent(out), old_fusion
                                 , dark_cumlai_mult          & ! intent(out), old_fusion
                                 , sunny_cumlai_min          & ! intent(out), old_fusion
                                 , sunny_cumlai_max          & ! intent(out), old_fusion
                                 , sunny_cumlai_mult         & ! intent(out), old_fusion
                                 , light_toler_min           & ! intent(out), old_fusion
                                 , light_toler_max           & ! intent(out), old_fusion
                                 , light_toler_mult          & ! intent(out), old_fusion
                                 , fuse_relax                & ! intent(out), old_fusion
                                 , lai_tol                   & ! intent(out)
                                 , pat_light_ext             & ! intent(out)
                                 , pat_light_tol_min         & ! intent(out)
                                 , pat_light_tol_max         & ! intent(out)
                                 , pat_light_tol_mult        & ! intent(out)
                                 , pat_light_mxd_fac         & ! intent(out)
                                 , pat_diff_age_tol          & ! intent(out)
                                 , pat_min_area_remain       & ! intent(out)
                                 , niter_cohfus              & ! intent(out)
                                 , coh_size_tol_min          & ! intent(out)
                                 , coh_size_tol_max          & ! intent(out)
                                 , coh_size_tol_mult         & ! intent(out)
                                 , corr_patch                & ! intent(out)
                                 , corr_cohort               & ! intent(out)
                                 , print_fuse_details        & ! intent(out)
                                 , fuse_prefix               & ! intent(out)
                                 , pat_laimax_fine           ! ! intent(out)
   use consts_coms        , only : onethird                  & ! intent(out)
                                 , twothirds                 & ! intent(in)
                                 , onesixth                  ! ! intent(in)
   use ed_max_dims        , only : n_dist_types              ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   real              :: exp_cohfus
   real              :: exp_patfus
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Parameters that control number of iterations.  More iterations mean slower         !
   ! increase in tolerance, which normally allows more fusion to occur at relatively       !
   ! stricter tolerance, although it may increase computational burden.                    !
   !                                                                                       !
   ! niter_patfus       -- number of patch fusion iterations                               !
   ! exp_patfus         -- exponential factor, used to determine the incremental           !
   !                       multiplication factor.                                          !
   ! niter_cohfus       -- number of cohort fusion iterations                              !
   ! exp_cohfus         -- exponential factor, used to determine the incremental           !
   !                       multiplication factor.                                          !
   !---------------------------------------------------------------------------------------!
   niter_cohfus = 100
   exp_cohfus   = 1. / (niter_cohfus - 1.0) 
   niter_patfus = 100
   exp_patfus   = 1. / (niter_patfus-1.0) 
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Old patch fusion variables (slated to be deleted in the near future).             !
   !---------------------------------------------------------------------------------------!
   dark_cumlai_min    = 5.5      ! Minimum cumulative LAI to be ignored (under storey)
   dark_cumlai_max    = 8.0      ! Maximum cumulative LAI to be ignored (under storey)
   sunny_cumlai_min   = 0.1      ! Minimum cumulative LAI to be ignored (top canopy)
   sunny_cumlai_max   = 0.3      ! Maximum cumulative LAI to be ignored (top canopy)
   light_toler_min    = 0.01     ! Minimum cumulative LAI to be ignored (under storey)
   light_toler_max    = onethird ! Maximum cumulative LAI to be ignored (under storey)
   !----- Multiplication factors. ---------------------------------------------------------!
   sunny_cumlai_mult  = (sunny_cumlai_max/sunny_cumlai_min)**exp_patfus
   dark_cumlai_mult   = (dark_cumlai_min /dark_cumlai_max )**exp_patfus
   light_toler_mult   = (light_toler_max /light_toler_min )**exp_patfus
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Old cohort fusion variables (slated to be deleted in the near future).            !
   !---------------------------------------------------------------------------------------!
   fusetol           =  0.4    ! Cohort fusion tolerance on DBH (dimensionless) 
   fusetol_h         =  0.5    ! Cohort fusion tolerance on height (m) !
   lai_fuse_tol      =  0.8    ! Cohort fusion tolerance on LAI (m2 leaf/m2 ground)
   coh_tolerance_max = 10.0    ! Cohort maximum tolerance factor 
   fuse_relax        = .false. ! Flag to allow a less strict fusion test
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Maximum LAI that a cohort is allowed to have.  This cap ensures that self-thin-  !
   ! ning mechanisms work in the model.  This parameter is used in fuse_cohorts and        !
   ! split_cohorts.                                                                        !
   !---------------------------------------------------------------------------------------!
   lai_tol            = 1.0
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Number of height classes (classes will be initialised in                         !
   ! init_derived_params_after_xml.                                                        !
   !---------------------------------------------------------------------------------------!
   select case (ifusion)
   case (1)
      ff_nhgt = 19
   case default
      ff_nhgt = 8
   end select
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Patch fusion layers were relocated to init_derived_params_after_xml, so the       !
   ! default numbers are based on the height of the tallest possible PFT.                  !
   !                                                                                       !
   ! pat_light_ext      -- extinction coefficient for light level.  Typically this should  !
   !                       be 0.5.                                                         !
   ! pat_light_tol_min  -- Minimum tolerance for average profile                           !
   ! pat_light_tol_max  -- Maximum tolerance for average profile                           !
   ! pat_light_tol_mult -- Factor that increments tolerance (derived from previous vari-   !
   !                       ables).                                                         !
   ! pat_light_mxd_fac  -- Maximum deviation from average tolerance (e.g. 1.25 means that  !
   !                       the maximum difference in light levels can be 25% greater than  !
   !                       tolerance for average maximum.                                  !
   !---------------------------------------------------------------------------------------!
   pat_light_ext      = 0.5
   pat_light_tol_min  = twothirds * 0.01
   pat_light_tol_max  = 0.10
   pat_light_tol_mult = (pat_light_tol_max/pat_light_tol_min)**exp_patfus
   pat_light_mxd_fac  = 1.50
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Cohort fusion variables.                                                          !
   !                                                                                       !
   ! coh_size_tol_min  -- Minimum tolerance for relative difference in size.               !
   ! coh_size_tol_max  -- Maximum tolerance for relative difference in size.               !
   ! coh_size_tol_mult -- Factor that increments tolerance (derived from previous vari-    !
   !                       ables).                                                         !
   !---------------------------------------------------------------------------------------!
   coh_size_tol_min  = twothirds * 0.01
   coh_size_tol_max  = twothirds * 0.10
   coh_size_tol_mult = (coh_size_tol_max/coh_size_tol_min)**exp_cohfus
   !---------------------------------------------------------------------------------------!



   !----- Maximum age difference allowed for two patches being considered same age [yr]. --!
   pat_diff_age_tol   = 0.999 / 12.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Minimum area to remain resolved.  This condition is normally met, except when    !
   ! initialising the simulation with massive amount of data (like airborne lidar data).   !
   !---------------------------------------------------------------------------------------!
   pat_min_area_remain = 0.90
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Coefficient of correlation assumed between two patches and cohorts that are      !
   ! about to be fused.                                                                    !
   !---------------------------------------------------------------------------------------!
   corr_patch  = 1.0
   corr_cohort = 1.0
   !---------------------------------------------------------------------------------------!


   !----- The following flag switches detailed debugging on. ------------------------------!
   print_fuse_details = .false.
   fuse_prefix        = 'patch_fusion_'
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Maximum patch-level LAI to be considered realistic during initialisation.  In     !
   ! very few cases, the airborne lidar initialisation algorithm predicts unreasonable     !
   ! total LAI (often in places with barely any return above the minimum height            !
   ! considered).  Because airborne lidar has tens of thousands of patches, it is hard to  !
   ! spot individual patches that didn't work, so we terminate these patches.              !
   !---------------------------------------------------------------------------------------!
   pat_laimax_fine = 20.
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_ff_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_disturb_params

   use disturb_coms , only : treefall_disturbance_rate & ! intent(in)
                           , include_fire              & ! intent(in)
                           , treefall_hite_threshold   & ! intent(out)
                           , does_hite_limit_tfpatch   & ! intent(out)
                           , forestry_on               & ! intent(out)
                           , agriculture_on            & ! intent(out)
                           , plantation_year           & ! intent(out)
                           , plantation_rotation       & ! intent(out)
                           , min_harvest_biomass       & ! intent(out)
                           , mature_harvest_age        & ! intent(out)
                           , min_oldgrowth             & ! intent(out)
                           , fire_dryness_threshold    & ! intent(out)
                           , fire_smoist_depth         & ! intent(out)
                           , fuel_height_max           & ! intent(out)
                           , f_combusted_fast_c        & ! intent(out)
                           , f_combusted_struct_c      & ! intent(out)
                           , f_combusted_fast_n        & ! intent(out)
                           , f_combusted_struct_n      & ! intent(out)
                           , k_fire_first              & ! intent(out)
                           , min_plantation_frac       & ! intent(out)
                           , max_plantation_dist       ! ! intent(out)
   use consts_coms  , only : erad                      & ! intent(in)
                           , pio180                    & ! intent(in)
                           , tiny_num                  & ! intent(in)
                           , huge_num                  ! ! intent(in)
   use soil_coms    , only : slz                       ! ! intent(in)
   use grid_coms    , only : nzg                       ! ! intent(in)
   implicit none

   !----- Only trees above this height create a gap when they fall. -----------------------!
   treefall_hite_threshold = 10.0

   !----- Flag to decide whether or not to limit disturbance to patches with tall trees. --!
   does_hite_limit_tfpatch = .true.

   !----- Set to 1 if to do forest harvesting. --------------------------------------------!
   forestry_on = 0

   !----- Set to 1 if to do agriculture. --------------------------------------------------!
   agriculture_on = 0

   !----- Earliest year at which plantations occur. ---------------------------------------!
   plantation_year = 1960

   !----- Number of years that a plantation requires to reach maturity. -------------------!
   plantation_rotation = 25.0

   !----- Minimum site biomass, below which harvest is skipped. ---------------------------!
   min_harvest_biomass = 0.001

   !----- Years that a non-plantation patch requires to reach maturity. -------------------!
   mature_harvest_age = 50.0

   !---------------------------------------------------------------------------------------!
   !     If include_fire is 1, then fire may occur if total (ground + underground) water   !
   ! converted to meters falls below this threshold.                                       !
   !---------------------------------------------------------------------------------------!
   fire_dryness_threshold = 0.2
   !---------------------------------------------------------------------------------------!

   !----- Maximum depth that will be considered in the average soil -----------------------!
   fire_smoist_depth     = -0.50
   !---------------------------------------------------------------------------------------!

   !----- Cut-off for fuel counting (used only when include_fire is 3). -------------------!
   fuel_height_max       = 2.0
   !---------------------------------------------------------------------------------------!

   !----- Fraction of biomass and necromass that are combusted and lost to air. -----------!
   select case (include_fire)
   case (3)
      f_combusted_fast_c   = 0.8
      f_combusted_struct_c = 0.5
      f_combusted_fast_n   = 0.72
      f_combusted_struct_n = 0.45
   case default
      f_combusted_fast_c   = 0.0
      f_combusted_struct_c = 0.0
      f_combusted_fast_n   = 0.0
      f_combusted_struct_n = 0.0
   end select
   !---------------------------------------------------------------------------------------!



   !----- Determine the top layer to consider for fires in case include_fire is 2 or 3. ---!
   kfireloop: do k_fire_first=nzg-1,1,-1
      if (slz(k_fire_first) < fire_smoist_depth) exit kfireloop
   end do kfireloop
   k_fire_first = k_fire_first + 1
   !---------------------------------------------------------------------------------------!



   !----- Minimum plantation fraction to consider the site a plantation. ------------------!
   min_plantation_frac = 0.125

   !---------------------------------------------------------------------------------------!
   !    Maximum distance to the current polygon that we still consider the file grid point !
   ! to be representative of the polygon.  The value below is 1.25 degree at the Equator.  !
   !---------------------------------------------------------------------------------------!
   max_plantation_dist = 1.25 * erad * pio180
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the minimum age above which we disregard the disturbance type because the    !
   ! patch can be considered old growth.                                                   !
   !---------------------------------------------------------------------------------------!
   !----- Non-cultivated patches: use the mean age for tree fall disturbances. ------------!
   if (abs(treefall_disturbance_rate) > tiny_num) then
      min_oldgrowth(:) = 1. / abs(treefall_disturbance_rate)
   else
      min_oldgrowth(:) = huge_num
   end if
   !----- Cultivated lands should never be considered old-growth. -------------------------!
   min_oldgrowth(1) = huge_num
   min_oldgrowth(2) = huge_num
   min_oldgrowth(8) = huge_num
   !---------------------------------------------------------------------------------------!

   return

end subroutine init_disturb_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine initialises various PFT-independent variables for the leaf physiology !
! model.  Some useful references for this sub-routine:                                     !
!                                                                                          !
! - M09 - Medvigy, D.M., S. C. Wofsy, J. W. Munger, D. Y. Hollinger, P. R. Moorcroft,      !
!         2009: Mechanistic scaling of ecosystem function and dynamics in space and time:  !
!         Ecosystem Demography model version 2.  J. Geophys. Res., 114, G01002,            !
!         doi:10.1029/2008JG000812.                                                        !
! - M06 - Medvigy, D.M., 2006: The State of the Regional Carbon Cycle: results from a      !
!         constrained coupled ecosystem-atmosphere model, 2006.  Ph.D. dissertation,       !
!         Harvard University, Cambridge, MA, 322pp.                                        !
! - M01 - Moorcroft, P. R., G. C. Hurtt, S. W. Pacala, 2001: A method for scaling          !
!         vegetation dynamics: the ecosystem demography model, Ecological Monographs, 71,  !
!         557-586.                                                                         !
! - F96 - Foley, J. A., I. Colin Prentice, N. Ramankutty, S. Levis, D. Pollard, S. Sitch,  !
!         A. Haxeltime, 1996: An integrated biosphere model of land surface processes,     !
!         terrestrial carbon balance, and vegetation dynamics. Glob. Biogeochem. Cycles,   !
!         10, 603-602.                                                                     !
! - L95 - Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf         !
!         nitrogen, photosynthesis, conductance, and transpiration: scaling from leaves to !
!         canopies. Plant, Cell, and Environ., 18, 1183-1200.                              !
! - F80 - Farquhar, G. D., S. von Caemmerer, J. A. Berry, 1980: A biochemical model of     !
!         photosynthetic  CO2 assimilation in leaves of C3 species. Planta, 149, 78-90.    !
! - C91 - Collatz, G. J., J. T. Ball, C. Grivet, J. A. Berry, 1991: Physiology and         !
!         environmental regulation of stomatal conductance, photosynthesis and transpir-   !
!         ation: a model that includes a laminar boundary layer, Agric. and Forest         !
!         Meteorol., 54, 107-136.                                                          !
! - C92 - Collatz, G. J., M. Ribas-Carbo, J. A. Berry, 1992: Coupled photosynthesis-       !
!         stomatal conductance model for leaves of C4 plants.  Aust. J. Plant Physiol.,    !
!         19, 519-538.                                                                     !
! - E78 - Ehleringer, J. R., 1978: Implications of quantum yield differences on the        !
!         distributions of C3 and C4 grasses.  Oecologia, 31, 255-267.                     !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine init_physiology_params()
   use detailed_coms  , only : idetailed           ! ! intent(in)
   use ed_misc_coms, only    : ffilout             ! ! intent(in)
   use physiology_coms, only : iphysiol            & ! intent(in)
                             , klowco2in           & ! intent(in)
                             , c34smin_lint_co2    & ! intent(out)
                             , c34smax_lint_co2    & ! intent(out)
                             , c34smax_gsw         & ! intent(out)
                             , gbh_2_gbw           & ! intent(out)
                             , gbw_2_gbc           & ! intent(out)
                             , gsw_2_gsc           & ! intent(out)
                             , gsc_2_gsw           & ! intent(out)
                             , kookc               & ! intent(out)
                             , tphysref            & ! intent(out)
                             , tphysrefi           & ! intent(out)
                             , fcoll               & ! intent(out)
                             , compp_refval        & ! intent(out)
                             , compp_hor           & ! intent(out)
                             , compp_q10           & ! intent(out)
                             , kco2_refval         & ! intent(out)
                             , kco2_hor            & ! intent(out)
                             , kco2_q10            & ! intent(out)
                             , ko2_refval          & ! intent(out)
                             , ko2_hor             & ! intent(out)
                             , ko2_q10             & ! intent(out)
                             , klowco2             & ! intent(out)
                             , o2_ref              & ! intent(out)
                             , par_lightcompp_max  & ! intent(out)
                             , qyield0             & ! intent(out)
                             , qyield1             & ! intent(out)
                             , qyield2             & ! intent(out)
                             , ehleringer_alpha0c  & ! intent(out)
                             , c34smin_lint_co28   & ! intent(out)
                             , c34smax_lint_co28   & ! intent(out)
                             , c34smax_gsw8        & ! intent(out)
                             , gbh_2_gbw8          & ! intent(out)
                             , gbw_2_gbc8          & ! intent(out)
                             , gsw_2_gsc8          & ! intent(out)
                             , gsc_2_gsw8          & ! intent(out)
                             , kookc8              & ! intent(out)
                             , tphysref8           & ! intent(out)
                             , tphysrefi8          & ! intent(out)
                             , fcoll8              & ! intent(out)
                             , compp_refval8       & ! intent(out)
                             , compp_hor8          & ! intent(out)
                             , compp_q108          & ! intent(out)
                             , kco2_refval8        & ! intent(out)
                             , kco2_hor8           & ! intent(out)
                             , kco2_q108           & ! intent(out)
                             , ko2_refval8         & ! intent(out)
                             , ko2_hor8            & ! intent(out)
                             , ko2_q108            & ! intent(out)
                             , klowco28            & ! intent(out)
                             , par_lightcompp_max8 & ! intent(out)
                             , o2_ref8             & ! intent(out)
                             , qyield08            & ! intent(out)
                             , qyield18            & ! intent(out)
                             , qyield28            & ! intent(out)
                             , ehleringer_alpha0c8 & ! intent(out)
                             , print_photo_debug   & ! intent(out)
                             , photo_prefix        ! ! intent(out)
   use consts_coms    , only : umol_2_mol          & ! intent(in)
                             , t00                 & ! intent(in)
                             , rmol                & ! intent(in)
                             , mmdoc               & ! intent(in)
                             , mmcod               & ! intent(in)
                             , mmo2                & ! intent(in)
                             , mmdryi              & ! intent(in)
                             , Watts_2_Ein         & ! intent(in)
                             , prefsea             & ! intent(in)
                             , solar               ! ! intent(in)
   implicit none
   !------ Local constants. ---------------------------------------------------------------!
   real, parameter :: ehl0           =  8.10e-2 ! Intercept                    (E78, eqn 2)
   real, parameter :: ehl1           = -5.30e-5 ! Linear coefficient           (E78, eqn 2)
   real, parameter :: ehl2           = -1.90e-5 ! Quadratic coefficient        (E78, eqn 2)
   real, parameter :: tau_refval_f96 =  4500.   ! Reference tau                (F96)
   real, parameter :: tau_refval_c91 =  2600.   ! Reference tau                (C91)
   real, parameter :: tau_hor        = -5000.   ! "Activation energy" for tau  (F96)
   real, parameter :: tau_q10        =  0.57    ! "Q10" term for tau           (C91)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Bounds for internal carbon and water stomatal conductance.                        !
   !---------------------------------------------------------------------------------------!
   c34smin_lint_co2 = 0.5   * umol_2_mol ! Minimum carbon dioxide concentration [  mol/mol]
   c34smax_lint_co2 = 1200. * umol_2_mol ! Maximum carbon dioxide concentration [  mol/mol]
   c34smax_gsw      = 1.e+2              ! Max. stomatal conductance (water)    [ mol/m2/s]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Many parameters in the model are temperature-dependent and we must provide a      !
   ! reference value and shape parameters.  Here we define this reference temperature.     !
   ! IMPORTANT: Some schemes use 15 C (F96-based in particular), whilst others use 25 C    !
   ! (C91-based in particular).  Be sure to make the necessary conversions of these        !
   ! references otherwise you may have surprises.                                          !
   !---------------------------------------------------------------------------------------!
   tphysref     = 15.0+t00
   tphysrefi    = 1./tphysref
   fcoll        = 0.10
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The following parameter is the concentration of oxygen, assumed constant through- !
   ! out the integration.  The value is in mol/mol.                                        !
   !---------------------------------------------------------------------------------------!
   o2_ref       = 0.209
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This parameter is from F80, and is the ratio between the turnover number for the  !
   ! oxygenase function and the turnover number for the carboxylase function.  For methods !
   ! 1 and 3, we use the ration from von Caemmerer (2000).                                 !
   !---------------------------------------------------------------------------------------!
   select case (iphysiol)
   case (1,3)
      kookc = 0.25
   case default
      kookc = 0.21
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Define the Michaelis-Mentel constants for CO2 and O2, and the CO2 compensation   !
   ! point without accounting leaf respiration (Gamma*).  Here we define the reference     !
   ! value and the compensation point parameters differently depending on the              !
   ! photosynthesis method.                                                                !
   !---------------------------------------------------------------------------------------!
   !----- Common parameters. --------------------------------------------------------------!
   select case (iphysiol)
   case (0)
      !------------------------------------------------------------------------------------!
      !    Original ED1, taken from IBIS (F96).                                            !
      !------------------------------------------------------------------------------------!

      !------ Parameters for CO2 Michaelis-Mentel constant. -------------------------------!
      kco2_hor    = 6000.             ! "Activation energy"                      [       K]
      kco2_q10    = 2.1               ! Q10 term                                 [     ---]
      kco2_refval = 150. * umol_2_mol ! Reference at 15 degC.                    [ mol/mol]
      !------------------------------------------------------------------------------------!


      !------ Parameters for O2 Michaelis-Mentel constant. --------------------------------!
      ko2_hor     = 1400.             ! "Activation energy"                      [       K]
      ko2_q10     = 1.2               ! Q10 term                                 [     ---]
      ko2_refval  = 0.250             ! Reference at 15 degC.                    [ mol/mol]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   Use the default CO2 compensation point values from F96/M01, based on tau.        !
      !   Gamma* = [O2] / (2. * tau)                                                       !
      !------------------------------------------------------------------------------------!
      compp_hor     = -tau_hor                       ! "Activation energy"       [       K]
      compp_q10     = 1. / tau_q10                   ! Q10 term                  [     ---]
      compp_refval  = o2_ref / (2. * tau_refval_f96) ! Reference at 15 degC.     [ mol/mol]
      !------------------------------------------------------------------------------------!

   case (2)
      !------------------------------------------------------------------------------------!
      !    Use default CO2 and O2 reference values from C91.  Here we must convert the     !
      ! reference values from Pa to mol/mol, and the reference must be set to 15 C.        !
      !------------------------------------------------------------------------------------!

      !------ Parameters for CO2 Michaelis-Mentel constant. -------------------------------!
      kco2_hor    = 6000.                            ! "Activation energy"       [       K]
      kco2_q10    = 2.1                              ! Q10 term                  [     ---]
      kco2_refval = 30. * mmcod / prefsea / kco2_q10 ! Reference at 15 degC.     [ mol/mol]
      !------------------------------------------------------------------------------------!


      !------ Parameters for O2 Michaelis-Mentel constant. --------------------------------!
      ko2_hor     = 1400.                                    ! "Act. energy"     [       K]
      ko2_q10     = 1.2                                      ! Q10 term          [     ---]
      ko2_refval  = 3.e4 * mmo2 * mmdryi / prefsea / ko2_q10 ! Ref. at 15 degC.  [ mol/mol]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   Use the default CO2 compensation point values from C91, based on tau.            !
      !   Gamma* = [O2] / (2. * tau)                                                       !
      !------------------------------------------------------------------------------------!
      compp_hor     = - tau_hor                                ! "Activation en. [       K]
      compp_q10     = 1. / tau_q10                             ! "Q10" term      [     ---]
      compp_refval  = o2_ref * tau_q10 / (2. * tau_refval_c91) ! Reference value [ mol/mol]
      !------------------------------------------------------------------------------------!

   case (1,3)
      !------------------------------------------------------------------------------------!
      !    Use default CO2 and O2 reference values from von Caemmerer (2000) and define    !
      ! compp assuming that Vomax = kookc * Vcmax, as she did.  Here we must convert the   !
      ! reference values from Pa to mol/mol, and the reference must be set to 15 C.        !
      !------------------------------------------------------------------------------------!

      !------ Parameters for CO2 Michaelis-Mentel constant. -------------------------------!
      kco2_hor    = 59360. * tphysref / (rmol * (t00+25.)) ! "Activation energy" [       K]
      kco2_q10    = 2.24                             ! Q10 term                  [     ---]
      kco2_refval = 40.4 / prefsea / kco2_q10        ! Reference at 15 degC.     [ mol/mol]
      !------------------------------------------------------------------------------------!


      !------ Parameters for O2 Michaelis-Mentel constant. --------------------------------!
      ko2_hor     = 59360. * tphysref / (rmol * (t00+25.)) ! "Activation energy" [       K]
      ko2_q10     = 1.63                                   ! Q10 term            [     ---]
      ko2_refval  = 24800. / prefsea / ko2_q10             ! Ref. at 15 degC.    [ mol/mol]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   Find the compensation point that is consistent with KCO2 and KO2, assuming that  !
      ! Vomax = 0.25 Vcmax (von Caemmerer 2000).                                           !
      !------------------------------------------------------------------------------------!
      compp_hor    = kco2_hor - ko2_hor                        ! "Activation E"  [       K]
      compp_q10    = kco2_q10 / ko2_q10                        ! "Q10" term      [     ---]
      compp_refval = kookc * o2_ref * kco2_refval                                          &
                   / (2. * ko2_refval)                         ! Reference value [ mol/mol]
      !------------------------------------------------------------------------------------!

   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    The following parameter is the k coefficient in Foley et al. (1996) that is used   !
   ! to determine the CO2-limited photosynthesis for C4 grasses.  Notice that Foley et al. !
   ! (1996) applied Vm0 so the slope is a function of temperature like in Collatz (1992).  !
   !---------------------------------------------------------------------------------------!
   klowco2      = klowco2in 
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These are constants obtained in Leuning et al. (1995) and Collatz et al. (1991)   !
   ! to convert different conductivities.                                                  !
   !---------------------------------------------------------------------------------------!
   gbh_2_gbw  = 1.075           ! heat  to water  - leaf boundary layer
   gbw_2_gbc  = 1.0 / 1.4       ! water to carbon - leaf boundary layer
   gsw_2_gsc  = 1.0 / 1.6       ! water to carbon - stomata
   gsc_2_gsw  = 1./gsw_2_gsc    ! carbon to water - stomata
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This is used to impose high light compensation point for when the electron        !
   ! transport rate (J0) that causes net assimilation rate to be zero is greater than      !
   ! Jmax.  This typically occurs at extremely high temperatures, in which case the        !
   ! stomata should remain closed.                                                         !
   !---------------------------------------------------------------------------------------!
   par_lightcompp_max = solar * 0.4 * Watts_2_Ein
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This is the quantum yield response curve, using equation 2 from E78:              !
   !     The coefficients here are different because his function used temperature in      !
   ! Celsius and here we give the coefficients for temperature in Kelvin.  We don't use    !
   ! this equation for temperatures below 0C, instead we find the equivalent at 0 Celsius  !
   ! and assume it is constant below this temperature.                                     !
   !---------------------------------------------------------------------------------------!
   qyield0            = ehl0 - ehl1 * t00 + ehl2 * t00 * t00
   qyield1            = ehl1 - 2. * ehl2 * t00
   qyield2            = ehl2
   ehleringer_alpha0c = ehl0
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Find the double precision version of the variables above.                          !
   !---------------------------------------------------------------------------------------!
   c34smin_lint_co28   = dble(c34smin_lint_co2  )
   c34smax_lint_co28   = dble(c34smax_lint_co2  )
   c34smax_gsw8        = dble(c34smax_gsw       )
   gbh_2_gbw8          = dble(gbh_2_gbw         )
   gbw_2_gbc8          = dble(gbw_2_gbc         )
   gsw_2_gsc8          = dble(gsw_2_gsc         )
   gsc_2_gsw8          = dble(gsc_2_gsw         )
   kookc8              = dble(kookc             )
   tphysref8           = dble(tphysref          )
   tphysrefi8          = dble(tphysrefi         )
   fcoll8              = dble(fcoll             )
   compp_refval8       = dble(compp_refval      )
   compp_hor8          = dble(compp_hor         )
   compp_q108          = dble(compp_q10         )
   kco2_refval8        = dble(kco2_refval       )
   kco2_hor8           = dble(kco2_hor          )
   kco2_q108           = dble(kco2_q10          )
   ko2_refval8         = dble(ko2_refval        )
   ko2_hor8            = dble(ko2_hor           )
   ko2_q108            = dble(ko2_q10           )
   klowco28            = dble(klowco2           )
   o2_ref8             = dble(o2_ref            )
   par_lightcompp_max8 = dble(par_lightcompp_max)
   qyield08            = dble(qyield0           )
   qyield18            = dble(qyield1           )
   qyield28            = dble(qyield2           )
   ehleringer_alpha0c8 = dble(ehleringer_alpha0c)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters that control debugging output.                                         !
   !---------------------------------------------------------------------------------------!
   !----- I should print detailed debug information. --------------------------------------!
   print_photo_debug = btest(idetailed,1)
   !----- File name prefix for the detailed information in case of debugging. -------------!
   photo_prefix      = trim(ffilout)//'_photo_state_'
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_physiology_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine defines the minimum and maximum acceptable values in the meteoro-    !
! logical forcing.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine init_met_params()
   use ed_misc_coms   , only : dtlsm           ! ! intent(in)
   use met_driver_coms, only : met_land_min    & ! intent(out)
                             , rshort_min      & ! intent(out)
                             , rshort_max      & ! intent(out)
                             , rlong_min       & ! intent(out)
                             , rlong_max       & ! intent(out)
                             , dt_radinterp    & ! intent(out)
                             , atm_tmp_min     & ! intent(out)
                             , atm_tmp_max     & ! intent(out)
                             , atm_shv_min     & ! intent(out)
                             , atm_shv_max     & ! intent(out)
                             , atm_rhv_min     & ! intent(out)
                             , atm_rhv_max     & ! intent(out)
                             , atm_co2_min     & ! intent(out)
                             , atm_co2_max     & ! intent(out)
                             , prss_min        & ! intent(out)
                             , prss_max        & ! intent(out)
                             , pcpg_min        & ! intent(out)
                             , pcpg_max        & ! intent(out)
                             , vels_min        & ! intent(out)
                             , vels_max        & ! intent(out)
                             , geoht_min       & ! intent(out)
                             , geoht_max       & ! intent(out)
                             , print_radinterp & ! intent(out)
                             , vbdsf_file      & ! intent(out)
                             , vddsf_file      & ! intent(out)
                             , nbdsf_file      & ! intent(out)
                             , nddsf_file      ! ! intent(out)

   !----- Minimum land fraction for a met driver point to be considered land. -------------!
   met_land_min = 0.5

   !----- Minimum and maximum acceptable shortwave radiation [W/m2]. ----------------------!
   rshort_min  = 0.
   rshort_max  = 1500.
   !----- Minimum and maximum acceptable longwave radiation [W/m2]. -----------------------!
   rlong_min   = 40.
   rlong_max   = 850.
   !----- Minimum and maximum acceptable air temperature    [   K]. -----------------------!
   atm_tmp_min = 184.     ! Lowest temperature ever measured, in Vostok Basin, Antarctica
   atm_tmp_max = 350.     ! About 78degC, or > 20degC higher than Death Valley's record
   !----- Minimum and maximum acceptable air specific humidity [kg_H2O/kg_air]. -----------!
   atm_shv_min = 1.0e-6   ! That corresponds to a relative humidity of 0.1% at 1000hPa
   atm_shv_max = 8.0e-2   ! That corresponds to a dew point of 41degC at 1000hPa.
   !----- Minimum and maximum acceptable CO2 mixing ratio [umol/mol]. ---------------------!
   atm_co2_min =   10.    ! Very low limit to allow Tonzi TS to run.  Plants may find it 
                          ! hard to thrive with such low CO2.
   atm_co2_max = 2000.    ! This should allow any RCP simulation to run with no problem
   !----- Minimum and maximum acceptable pressure [Pa]. -----------------------------------!
   prss_min =  45000.     ! It may crash if you run a simulation in Mt. Everest.
   prss_max = 110000.     ! It may crash if you run a simulation under water.
   !----- Minimum and maximum acceptable precipitation rates [kg/m2/s]. -------------------!
   pcpg_min     = 0.0     ! No negative precipitation is allowed
   pcpg_max     = 0.1111  ! This is a precipitation rate of 400mm/hr.
   !----- Minimum and maximum acceptable wind speed [m/s]. --------------------------------!
   vels_min     =  0.0    ! No negative wind is acceptable.
   vels_max     = 85.0    ! Maximum sustained winds recorded during Typhoon Tip (1970).
   !----- Minimum and maximum reference heights [m]. --------------------------------------!
   geoht_min    =   1.0   ! This should be above-canopy measurement, but 1.0 is okay for
                          !     grasslands...
   geoht_max    = 350.0   ! This should be not that much above the canopy, but tall towers
                          !     do exist...
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Minimum and maximum acceptable relative humidity (fraction).  This is not going   !
   ! cause the simulation to crash, instead it will just impose these numbers to the       !
   ! meteorological forcing.                                                               !
   !---------------------------------------------------------------------------------------!
   atm_rhv_min = 5.e-3 ! 0.5%
   atm_rhv_max = 1.0   ! 100.0%.  Although canopy air space can experience super-
                       !    saturation, we don't allow the air above to be super-saturated.
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Time step used to perform the daytime average of the secant of the zenith angle.  !
   !---------------------------------------------------------------------------------------!
   dt_radinterp = dtlsm    ! Value in seconds.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !   These variables control the detailed interpolation output (for debugging only).     !
   !---------------------------------------------------------------------------------------!
   print_radinterp = .false.
   vbdsf_file      = 'visible_beam.txt'
   vddsf_file      = 'visible_diff.txt'
   nbdsf_file      = 'near_infrared_beam.txt'
   nddsf_file      = 'near_infrared_diff.txt'
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_met_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine defines defaults for lapse rates.                                    !
!------------------------------------------------------------------------------------------!
subroutine init_lapse_params()

   use met_driver_coms, only : lapse             & ! intent(out)
                             , atm_tmp_intercept & ! intent(out)
                             , atm_tmp_slope     & ! intent(out)
                             , prec_intercept    & ! intent(out)
                             , prec_slope        & ! intent(out)
                             , humid_scenario    ! ! intent(out)

   lapse%geoht        = 0.0
   lapse%atm_ustar    = 0.0
   lapse%vels         = 0.0
   lapse%atm_tmp      = 0.0
   lapse%atm_theta    = 0.0
   lapse%atm_theiv    = 0.0
   lapse%atm_vpdef    = 0.0
   lapse%atm_shv      = 0.0
   lapse%prss         = 0.0
   lapse%pcpg         = 0.0
   lapse%atm_co2      = 0.0
   lapse%rlong        = 0.0
   lapse%nir_beam     = 0.0
   lapse%nir_diffuse  = 0.0
   lapse%par_beam     = 0.0
   lapse%par_diffuse  = 0.0
   lapse%pptnorm      = 0.0

   atm_tmp_intercept = 0.0
   atm_tmp_slope     = 1.0
   prec_intercept    = 0.0
   prec_slope        = 1.0
   humid_scenario    = 0

   return
end subroutine init_lapse_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_hydro_coms

   use ed_misc_coms  , only : ied_init_mode       ! ! intent(in)
   use hydrology_coms, only : useTOPMODEL         & ! intent(out)
                            , useRUNOFF           & ! intent(out)
                            , HydroOutputPeriod   & ! intent(out)
                            , MoistRateTuning     & ! intent(out)
                            , MoistSatThresh      & ! intent(out)
                            , Moist_dWT           & ! intent(out)
                            , FracLiqRunoff       & ! intent(out)
                            , GrassLAIMax         & ! intent(out)
                            , inverse_runoff_time ! ! intent(out)

   implicit none

   select case (ied_init_mode)
   case (3)
      ! Signifies a restart from an ED2 history file
      useTOPMODEL = 1
      useRUNOFF   = 0
   case default
      useTOPMODEL = 0
      useRUNOFF = 0
   end select

   HydroOutputPeriod = 96 !! multiples of dtlsm

   MoistRateTuning = 1.0

   MoistSatThresh = 0.95

   Moist_dWT = 2.0

   FracLiqRunoff = 0.5

   GrassLAImax = 4.0

   inverse_runoff_time = 0.1

   return
end subroutine init_hydro_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    Subroutine that initialises most of the soil parameters.                              !
!                                                                                          !
! MLO: This sub-routine formerly initiliased the soil and soil8 tables, but this creates   !
!      a problem for some HISTORY runs (especially those with multiple sites that modify   !
!      the default parameters).  The soil table is now initialised in a separate sub-      !
!      routine (ed_gen_soil_table), after loading HISTORY variables, or right before or    !
!      right after reading the initial conditions.                                         !
!------------------------------------------------------------------------------------------!
subroutine init_soil_coms
   use soil_coms      , only : isoilflg              & ! intent(in)
                             , nslcon                & ! intent(in)
                             , soil_hydro_scheme     & ! intent(in)
                             , slxclay               & ! intent(in)
                             , slxsand               & ! intent(in)
                             , slsoc                 & ! intent(in)
                             , slph                  & ! intent(in)
                             , slcec                 & ! intent(in)
                             , sldbd                 & ! intent(in)
                             , slxkey_ref            & ! intent(out)
                             , slhydro_ref           & ! intent(out)
                             , slxclay_ref           & ! intent(out)
                             , slxsilt_ref           & ! intent(out)
                             , slxsand_ref           & ! intent(out)
                             , slsoc_ref             & ! intent(out)
                             , slph_ref              & ! intent(out)
                             , slcec_ref             & ! intent(out)
                             , sldbd_ref             & ! intent(out)
                             , soilcol               & ! intent(in)
                             , soilcol_class         & ! type
                             , water_stab_thresh     & ! intent(out)
                             , snowmin               & ! intent(out)
                             , dewmax                & ! intent(out)
                             , soil_rough            & ! intent(out)
                             , snow_rough            & ! intent(out)
                             , ny07_eq04_a           & ! intent(out)
                             , ny07_eq04_m           & ! intent(out)
                             , tiny_sfcwater_mass    & ! intent(out)
                             , infiltration_method   & ! intent(out)
                             , soil_rough8           & ! intent(out)
                             , snow_rough8           & ! intent(out)
                             , ny07_eq04_a8          & ! intent(out)
                             , ny07_eq04_m8          & ! intent(out)
                             , freezecoef            & ! intent(out)
                             , freezecoef8           & ! intent(out)
                             , sldrain               & ! intent(out)
                             , sldrain8              & ! intent(out)
                             , sin_sldrain           & ! intent(out)
                             , sin_sldrain8          & ! intent(out)
                             , hydcond_min           & ! intent(out)
                             , hydcond_min8          ! ! intent(out)
   use grid_coms      , only : ngrids                ! ! intent(in)
   use consts_coms    , only : wdns                  & ! intent(in)
                             , day_sec               & ! intent(in)
                             , pio180                & ! intent(in)
                             , pio1808               ! ! intent(in)

   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   logical      :: update_slx                        ! Update texture fractions?  [    T|F]
   !----- Local constants. ----------------------------------------------------------------!
   real(kind=4), parameter :: residual_K  =  1.e-5   ! min. hydr. cond. (RS02)  [mm/day]
   !---------------------------------------------------------------------------------------!


   !----- Initialise some standard variables. ---------------------------------------------!
   water_stab_thresh   = 10.0           ! Minimum mass to be considered stable    [   kg/m2]
   snowmin             = 5.0            ! Minimum mass needed to create new layer [   kg/m2]
   dewmax              = 3.0e-5         ! Maximum dew flux rate (deprecated)      [ kg/m2/s]
   soil_rough          = 0.01           ! Soil roughness height                   [       m]
   snow_rough          = 0.0024         ! Snowcover roughness height              [       m]
   tiny_sfcwater_mass  = 1.0e-3         ! Minimum mass in temporary layers        [   kg/m2]
   infiltration_method = 0              ! Infiltration method, used in rk4_derivs [     0|1]
   freezecoef          = 7.0 * log(10.) ! Coeff. for infiltration of frozen water [     ---]
   hydcond_min         = residual_K / wdns /day_sec
                                        ! Residual hydraulic conductivity         [     m/s]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters for fraction covered with snow, which is based on:                     !
   !                                                                                       !
   ! Niu, G.-Y., and Z.-L. Yang (2007), An observation-based formulation of snow cover     !
   !    fraction and its evaluation over large North American river basins,                !
   !    J. Geophys. Res., 112, D21101, doi:10.1029/2007JD008674                            !
   !                                                                                       !
   !    These are the parameters in equation 4.  Fresh snow density is defined at          !
   ! consts_coms.f90                                                                       !
   !---------------------------------------------------------------------------------------!
   ny07_eq04_a = 2.5   ! the coefficient next to soil roughness.
   ny07_eq04_m = 1.0   ! m
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Removed the hardcoded initialisation of the entire structure.  Instead, we set the !
   ! texture for every class, then use the equations to populate the structure.            !
   !---------------------------------------------------------------------------------------!
   slxkey_ref ( :) = (/'  Sa',' LSa',' SaL',' SiL','   L','SaCL','SiCL','  CL'             &
                      ,' SaC',' SiC','   C','Peat','BdRk','  Si','  CC',' CSa',' CSi' /)
   slxsand_ref( :) = (/ 0.920, 0.825, 0.660, 0.200, 0.410, 0.590, 0.100, 0.320             &
                      , 0.520, 0.060, 0.200, 0.200, 0.333, 0.075, 0.100, 0.375, 0.125 /)
   slxclay_ref( :) = (/ 0.030, 0.060, 0.110, 0.160, 0.170, 0.270, 0.340, 0.340             &
                      , 0.420, 0.470, 0.600, 0.200, 0.333, 0.050, 0.800, 0.525, 0.525 /)
   slxsilt_ref( :) = 1. - slxsand_ref(:) - slxclay_ref(:)
   slhydro_ref( :) = soil_hydro_scheme
   slhydro_ref(12) = 12
   slhydro_ref(13) = 13
   slsoc_ref  ( :) = slsoc
   slph_ref   ( :) = slph
   slcec_ref  ( :) = slcec
   sldbd_ref  ( :) = sldbd
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Check whether or not to overwrite soil texture.                                   !
   !---------------------------------------------------------------------------------------!
   update_slx = any(isoilflg(1:ngrids) == 2) .and. slxclay > 0. .and. slxsand > 0. .and.   &
                (slxclay + slxsand) <= 1.
   if (update_slx) then
      slxsand_ref(nslcon) = slxsand
      slxclay_ref(nslcon) = slxclay
      slxsilt_ref(nslcon) = 1. - slxsand - slxclay
      slxkey_ref (nslcon) = 'User'
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Fill in the albedo information regarding the soil colour classes.                 !
   !---------------------------------------------------------------------------------------!
   !                    |    Dry soil   |   Saturated   | Emis. |                          !
   !   Soil class       |---------------+---------------|-------|                          !
   !                    |   VIS |   NIR |   VIS |   NIR |   TIR |                          !
   !---------------------------------------------------------------------------------------!
   soilcol = (/                                                     & !
       soilcol_class   (    0.36,   0.61,   0.25,   0.50,   0.98 )  & ! 01 - Brightest
      ,soilcol_class   (    0.34,   0.57,   0.23,   0.46,   0.98 )  & ! 02
      ,soilcol_class   (    0.32,   0.53,   0.21,   0.42,   0.98 )  & ! 03
      ,soilcol_class   (    0.31,   0.51,   0.20,   0.40,   0.98 )  & ! 04
      ,soilcol_class   (    0.30,   0.49,   0.19,   0.38,   0.98 )  & ! 05
      ,soilcol_class   (    0.29,   0.48,   0.18,   0.36,   0.98 )  & ! 06
      ,soilcol_class   (    0.28,   0.45,   0.17,   0.34,   0.98 )  & ! 07
      ,soilcol_class   (    0.27,   0.43,   0.16,   0.32,   0.98 )  & ! 08
      ,soilcol_class   (    0.26,   0.41,   0.15,   0.30,   0.98 )  & ! 09
      ,soilcol_class   (    0.25,   0.39,   0.14,   0.28,   0.98 )  & ! 10
      ,soilcol_class   (    0.24,   0.37,   0.13,   0.26,   0.98 )  & ! 11
      ,soilcol_class   (    0.23,   0.35,   0.12,   0.24,   0.98 )  & ! 12
      ,soilcol_class   (    0.22,   0.33,   0.11,   0.22,   0.98 )  & ! 13
      ,soilcol_class   (    0.20,   0.31,   0.10,   0.20,   0.98 )  & ! 14
      ,soilcol_class   (    0.18,   0.29,   0.09,   0.18,   0.98 )  & ! 15
      ,soilcol_class   (    0.16,   0.27,   0.08,   0.16,   0.98 )  & ! 16
      ,soilcol_class   (    0.14,   0.25,   0.07,   0.14,   0.98 )  & ! 17
      ,soilcol_class   (    0.12,   0.23,   0.06,   0.12,   0.98 )  & ! 18
      ,soilcol_class   (    0.10,   0.21,   0.05,   0.10,   0.98 )  & ! 19
      ,soilcol_class   (    0.08,   0.16,   0.04,   0.08,   0.98 )  & ! 20 - Darkest
      ,soilcol_class   (    0.00,   0.00,   0.00,   0.00,   0.98 )  & ! 21 - ED-2.1, unused
      /)
   !---------------------------------------------------------------------------------------!



   !----- Double precision of additional scalar variables. --------------------------------!
   soil_rough8  = dble(soil_rough )
   snow_rough8  = dble(snow_rough )
   ny07_eq04_a8 = dble(ny07_eq04_a)
   ny07_eq04_m8 = dble(ny07_eq04_m)
   freezecoef8  = dble(freezecoef )
   hydcond_min8 = dble(hydcond_min)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the double precision version of the drainage slope, and find and save the    !
   ! sine of it.                                                                           !
   !---------------------------------------------------------------------------------------!
   sldrain8     = dble(sldrain)
   sin_sldrain  = sin(sldrain  * pio180 )
   sin_sldrain8 = sin(sldrain8 * pio1808)
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_soil_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_phen_coms
   use consts_coms   , only : erad                     & ! intent(in)
                            , pio180                   ! ! intent(in)
   use phenology_coms, only : thetacrit                & ! intent(in)
                            , retained_carbon_fraction & ! intent(out)
                            , root_phen_factor         & ! intent(out)
                            , f_psi_xdry               & ! intent(out)
                            , elongf_min               & ! intent(out)
                            , elongf_flush             & ! intent(out)
                            , spot_phen                & ! intent(out)
                            , dl_tr                    & ! intent(out)
                            , st_tr1                   & ! intent(out)
                            , st_tr2                   & ! intent(out)
                            , phen_a                   & ! intent(out)
                            , phen_b                   & ! intent(out)
                            , phen_c                   & ! intent(out)
                            , max_phenology_dist       & ! intent(out)
                            , radavg_window            & ! intent(out)
                            , turnamp_window           & ! intent(out)
                            , turnamp_min              & ! intent(out)
                            , turnamp_max              & ! intent(out)
                            , llspan_window            & ! intent(out)
                            , llspan_min               & ! intent(out)
                            , llspan_max               & ! intent(out)
                            , llspan_inf               & ! intent(out)
                            , vm0_window               & ! intent(out)
                            , vm0_tran                 & ! intent(out)
                            , vm0_slope                & ! intent(out)
                            , vm0_amp                  & ! intent(out)
                            , vm0_min                  & ! intent(out)
                            , sla_window               ! ! intent(out)
   use ed_misc_coms  , only : economics_scheme         ! ! intent(out)
   implicit none

   !---------------------------------------------------------------------------------------!
   !     Before plants drop their leaves, they retain this fraction of their leaf carbon   !
   ! and nitrogen and put it into storage.                                                 !
   !---------------------------------------------------------------------------------------!
   retained_carbon_fraction = 0.5  ! XX-> Meta-analysis from Vergutz et al. 2012 reports ~0.25
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Factor that controls the fine-root "elongation factor" relative to leaf           !
   ! elongation factor.  This is currently applied only for PFTs with phenology = 5.       !
   !                                                                                       !
   ! e_root = (e_leaf + root_phen_factor - 1) / root_phen_factor.                          !
   !                                                                                       !
   ! root_phen_factor > 1.  Fine roots will senesce more slowly than leaf shedding.        !
   ! root_phen_factor = 1.  Fine root elongation factor will be the same as for leaves.    !
   ! root_phen_factor < 1.  Fine roots will senesce more rapidly than leaf shedding.       !
   ! root_phen_factor = 0.  Special flag to disable fine-root phenology.                   !
   ! root_phen_factor < 0.  Non-sensical, currently assume the same as 0.                  !
   !---------------------------------------------------------------------------------------!
   root_phen_factor = 2.0
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Threshold for shedding all leaves when leaf water potential is very low. .       !
   !---------------------------------------------------------------------------------------!
   f_psi_xdry               = 0.95
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Minimum elongation factor before plants give up completely and shed all remain-  !
   ! ing leaves.                                                                           !
   !---------------------------------------------------------------------------------------!
   elongf_min               = 0.05
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Minimum elongation factor that allows plants to start flushing out new leaves if !
   ! they are drought deciduous and have been losing leaves.                               !
   !---------------------------------------------------------------------------------------!
   elongf_flush             = 0.25
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Flag that checks whether to Use soil potential rather than soil moisture to drive !
   ! phenology.                                                                            !
   !---------------------------------------------------------------------------------------!
   spot_phen                = thetacrit < 0.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Leaf offset parameters are from:                                                 !
   !      White et al. 1997, Global Biogeochemical Cycles 11(2) 217-234                    !
   !---------------------------------------------------------------------------------------!
   dl_tr                    = 655.0
   st_tr1                   = 284.3
   st_tr2                   = 275.15
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Phenology parameters for cold deciduous trees:                                   !
   !      Botta et al. 2000, Global Change Biology, 6, 709--725                            !
   !---------------------------------------------------------------------------------------!
   phen_a                   = -68.0
   phen_b                   = 638.0
   phen_c                   = -0.01
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     This variable is the maximum distance between the coordinates of a prescribed     !
   ! phenology file and the actual polygon that we will still consider close enough to be  !
   ! representative.  If the user wants to run with prescribed phenology and the closest   !
   ! file is farther away from the polygon than the number below, the simulation will      !
   ! stop.                                                                                 !
   !---------------------------------------------------------------------------------------!
   max_phenology_dist       = 1.25 * erad * pio180
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Variables controlling the light phenology.                                       !
   !---------------------------------------------------------------------------------------!
   select case (economics_scheme)
   case (1)
      !----- Radiation window for running average [days] ----------------------------------!
      radavg_window  = 14.
      !----- Turnover window for running average [days] -----------------------------------!
      turnamp_window = radavg_window
      !----- Minimum instantaneous turnover rate amplitude [n/d]. -------------------------!
      turnamp_min    = 0.01
      !----- Maximum instantaneous turnover rate amplitude [n/d]. -------------------------!
      turnamp_max    = 100.
      !----- Lifespan window for running average [days]. ----------------------------------!
      llspan_window   = 14.
      !----- Minimum instantaneous life span [months]. ------------------------------------!
      llspan_min      = 2.0
      !----- Maximum instantaneous life span [months]. ------------------------------------!
      llspan_max      = 60.
      !----- Instantaneous life span in case the turnover rate is 0. ----------------------!
      llspan_inf      = 9999.
      !----- Vm0 window for running average [days]. ---------------------------------------!
      vm0_window      = llspan_window
      !----- SLA window for running average [days]. ---------------------------------------!
      sla_window      = llspan_window
      !------------------------------------------------------------------------------------!
   case default
      !----- Radiation window for running average [days] ----------------------------------!
      radavg_window  = 10.
      !----- Turnover window for running average [days] -----------------------------------!
      turnamp_window = radavg_window
      !----- Minimum instantaneous turnover rate amplitude [n/d]. -------------------------!
      turnamp_min    = 0.01
      !----- Maximum instantaneous turnover rate amplitude [n/d]. -------------------------!
      turnamp_max    = 100.
      !----- Lifespan window for running average [days]. ----------------------------------!
      llspan_window   = 60.
      !----- Minimum instantaneous life span [months]. ------------------------------------!
      llspan_min      = 2.0
      !----- Maximum instantaneous life span [months]. ------------------------------------!
      llspan_max      = 60.
      !----- Instantaneous life span in case the turnover rate is 0. ----------------------!
      llspan_inf      = 9999.
      !----- Vm0 window for running average [days]. ---------------------------------------!
      vm0_window      = llspan_window
      !----- SLA window for running average [days]. ---------------------------------------!
      sla_window      = llspan_window
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Parameters to define the instantaneous Vm0 as a function of leaf life span.      !
   ! These parameters are only used when economics_scheme is zero.                         !
   !---------------------------------------------------------------------------------------!
   vm0_tran        = 1.98   ! 8.5
   vm0_slope       = 6.53   ! 7.0
   vm0_amp         = 57.2   ! 42.0
   vm0_min         = 7.31   ! 18.0
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_phen_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine assigns values for some variables that are in ed_misc_coms, which    !
! wouldn't fit in any of the other categories.                                             !
!------------------------------------------------------------------------------------------!
subroutine init_ed_misc_coms
   use ed_max_dims  , only : n_pft                & ! intent(in)
                           , n_dbh                & ! intent(in)
                           , n_age                ! ! intent(in)
   use consts_coms  , only : erad                 & ! intent(in)
                           , pio180               ! ! intent(in)
   use ed_misc_coms , only : burnin               & ! intent(out)
                           , restart_target_year  & ! intent(out)
                           , use_target_year      & ! intent(out)
                           , maxage               & ! intent(out)
                           , dagei                & ! intent(out)
                           , maxdbh               & ! intent(out)
                           , ddbhi                & ! intent(out)
                           , vary_elev            & ! intent(out)
                           , vary_hyd             & ! intent(out)
                           , vary_rad             & ! intent(out)
                           , max_thsums_dist      & ! intent(out)
                           , max_poihist_dist     & ! intent(out)
                           , max_poi99_dist       & ! intent(out)
                           , suppress_h5_warnings & ! intent(out)
                           , use_efrd_trtree      ! ! intent(out)
   implicit none


   !----- Flags that allow components of subgrid heterogeneity to be turned on/off --------!
   vary_elev = 1
   vary_rad  = 1
   vary_hyd  = 1
   !---------------------------------------------------------------------------------------!


   !----- Number of years to ignore demography when starting a run. -----------------------!
   burnin = 0
   !---------------------------------------------------------------------------------------!


   !----- Year to read when parsing pss/css with multiple years. --------------------------!
   restart_target_year = 2000
   !---------------------------------------------------------------------------------------!


   !----- Flag specifying whether to search for a target year in pss/css. -----------------!
   use_target_year = 0
   !---------------------------------------------------------------------------------------!



   !----- Maximum age [yr] to split into classes. -----------------------------------------!
   maxage = 200.
   !---------------------------------------------------------------------------------------!



   !----- Maximum DBH [cm] to be split into classes. --------------------------------------!
   maxdbh = 100.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The inverse of bin classes will depend on max??? and n_???, leaving one class for !
   ! when the number exceeds the maximum.                                                  !
   !---------------------------------------------------------------------------------------!
   dagei = real(n_age-1) / maxage
   ddbhi = real(n_dbh-1) / maxdbh
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Maximum distance to the current polygon that we still consider the file grid point !
   ! to be representative of the polygon.  The value below is 1.25 degree at the Equator.  !
   !---------------------------------------------------------------------------------------!
   max_thsums_dist    = 1.25 * erad * pio180
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Alternative method for mixing 1 grid and POI's.  Only use the grid if their is   !
   ! NOT an POI  within a user specified resolution.  Remember, this assumes there is only !
   ! 1 gridded file, and it is the first file when ied_init_mode is set to 99  (Developer  !
   ! use only).                                                                            !
   !---------------------------------------------------------------------------------------!
   max_poi99_dist     = 5.0 * erad * pio180
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      This variable is used for the history start initialisation.  This sets the       !
   ! maximum acceptable distance between the expected polygon and the polygon found in the !
   ! history file.  Units: m.                                                              !
   !---------------------------------------------------------------------------------------!
   max_poihist_dist   = 250.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      If you don't want to read a million warnings about certain initialization        !
   ! variables not being available in the history input file, set this to .true. .  It's   !
   ! better for new users to see what is missing though.                                   !
   !---------------------------------------------------------------------------------------!
   suppress_h5_warnings = .true.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Temporary flag: decide whether to use the effective functional rooting depth     !
   ! based on delta 18O (true) following Brum et al. (2018); or the direct size-rooting    !
   ! depth function (false).  Note that both approaches ultimately depend on size.         !
   !---------------------------------------------------------------------------------------!
   use_efrd_trtree = .false.
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_ed_misc_coms
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will assign some parameters used by the horizontal shading scheme.    !
!------------------------------------------------------------------------------------------!
subroutine init_hrzshade_params()

   use canopy_radiation_coms , only : cci_radius                  & ! intent(out)
                                    , cci_pixres                  & ! intent(out)
                                    , cci_gapsize                 & ! intent(out)
                                    , cci_gapmin                  & ! intent(out)
                                    , cci_nretn                   & ! intent(out)
                                    , cci_hmax                    & ! intent(out)
                                    , fixed_hrz_classes           & ! intent(out)
                                    , at_bright_def               & ! intent(out)
                                    , at_dark_def                 & ! intent(out)
                                    , at0                         & ! intent(out)
                                    , at1                         & ! intent(out)
                                    , at08                        & ! intent(out)
                                    , at18                        ! ! intent(out)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=4) :: cci_bright_def
   real(kind=4) :: cci_dark_def
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     The following variables control the method that allow light redistribution based  !
   ! on patch neighbourhood                                                                !
   !---------------------------------------------------------------------------------------!
   cci_radius   = 10.0 ! Maximum radius to calculate CCI                           [     m]
   cci_pixres   =  1.0 ! Pixel resolution for TCH and CCI                          [     m]
   cci_gapsize  = 14.0 ! Gap size                                                  [     m]
   cci_gapmin   = 30.0 ! # of gaps associated with the smallest area               [   ---]
   cci_nretn    = 30   ! "Return density" to generate the TCH map                  [  1/m2]
   cci_hmax     = 70.0 ! Maximum height allowed in the CCI scheme                  [     m]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Coefficients that control CCI vs. light-correction curve.                         !
   !---------------------------------------------------------------------------------------!
   at0       =  3.012569
   at1       = -0.0044086
   at08      = dble(at0)
   at18      = dble(at1)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Parameters that control the gap splitting regarding light conditions.             !
   !---------------------------------------------------------------------------------------!
   fixed_hrz_classes = .true. ! If true, use thresholds below instead of quantiles  [  T|F]
   cci_bright_def    = 75.    ! Default threshold below which gaps are bright       [  ---]
   cci_dark_def      = 150.   ! Default threshold above which gaps are dark         [  ---]
   !---------------------------------------------------------------------------------------!


   !----- Derive the default fbeam thresholds for bright and dark. ------------------------!
   at_bright_def = exp(at0+at1*cci_bright_def)
   at_dark_def   = exp(at0+at1*cci_dark_def  )
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_hrzshade_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_alloc_params()

   use pft_coms       , only : is_tropical           & ! intent(in)
                             , is_savannah           & ! intent(in)
                             , is_conifer            & ! intent(in)
                             , is_liana              & ! intent(in)
                             , is_grass              & ! intent(in)
                             , rho                   & ! intent(out)
                             , SLA                   & ! intent(out)
                             , leaf_turnover_rate    & ! intent(out)
                             , q                     & ! intent(out)
                             , qsw                   & ! intent(out)
                             , qbark                 & ! intent(out)
                             , qrhob                 & ! intent(out)
                             , SRA                   & ! intent(out)
                             , root_beta             & ! intent(out)
                             , init_density          & ! intent(out)
                             , init_laimax           & ! intent(out)
                             , agf_bs                & ! intent(out)
                             , brf_wd                & ! intent(out)
                             , hgt_min               & ! intent(out)
                             , hgt_ref               & ! intent(out)
                             , hgt_max               & ! intent(out)
                             , min_dbh               & ! intent(out)
                             , dbh_crit              & ! intent(out)
                             , dbh_bigleaf           & ! intent(out)
                             , min_bdead             & ! intent(out)
                             , bdead_crit            & ! intent(out)
                             , b1Ht                  & ! intent(out)
                             , b2Ht                  & ! intent(out)
                             , b1Bs_small            & ! intent(out)
                             , b2Bs_small            & ! intent(out)
                             , b1Bs_large            & ! intent(out)
                             , b2Bs_large            & ! intent(out)
                             , b1Ca                  & ! intent(out)
                             , b2Ca                  & ! intent(out)
                             , b1Cl                  & ! intent(out)
                             , b2Cl                  & ! intent(out)
                             , b1Rd                  & ! intent(out)
                             , b2Rd                  & ! intent(out)
                             , d18O_ref              & ! intent(out)
                             , b1d18O                & ! intent(out)
                             , b2d18O                & ! intent(out)
                             , b1Efrd                & ! intent(out)
                             , b2Efrd                & ! intent(out)
                             , b1Vol                 & ! intent(out)
                             , b2Vol                 & ! intent(out)
                             , b1Bl                  & ! intent(out)
                             , b2Bl                  & ! intent(out)
                             , b1WAI                 & ! intent(out)
                             , b2WAI                 & ! intent(out)
                             , b1SA                  & ! intent(out)
                             , b2SA                  & ! intent(out)
                             , b1Xs                  & ! intent(out)
                             , b1Xb                  & ! intent(out)
                             , C2B                   & ! intent(out)
                             , sla_s0                & ! intent(out)
                             , sla_s1                & ! intent(out)
                             , kplastic_SLA          & ! intent(out)
                             , eplastic_vm0          & ! intent(out)
                             , eplastic_sla          & ! intent(out)
                             , kplastic_LL           & ! intent(out)
                             , laimax_plastic        & ! intent(out)
                             , LMA_slope             & ! intent(out)
                             , sapwood_ratio         & ! intent(out)
                             , f_bstorage_init       & ! intent(out)
                             , leaf_width            & ! intent(out)
                             , branch_diam           & ! intent(out)
                             , h_edge                & ! intent(out)
                             , liana_dbh_crit        & ! intent(out)
                             , nbt_lut               ! ! intent(out)
   use allometry      , only : h2dbh                 & ! function
                             , dbh2h                 & ! function
                             , size2bd               & ! function
                             , size2bl               ! ! function
   use consts_coms    , only : onethird              & ! intent(in)
                             , onesixth              & ! intent(in)
                             , twothirds             & ! intent(in)
                             , huge_num              & ! intent(in)
                             , pi1                   ! ! intent(in)
   use ed_max_dims    , only : n_pft                 & ! intent(in)
                             , str_len               & ! intent(in)
                             , undef_real            ! ! intent(in)
   use ed_misc_coms   , only : iallom                & ! intent(in)
                             , economics_scheme      & ! intent(in)
                             , ibigleaf              & ! intent(in)
                             , ivegt_dynamics        ! ! intent(in)
   use canopy_air_coms, only : lwidth_grass          & ! intent(in)
                             , lwidth_bltree         & ! intent(in)
                             , lwidth_nltree         ! ! intent(in)

   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer                   :: ipft
   real   , dimension(n_pft) :: abas
   real                      :: aux
   real                      :: init_bleaf
   real                      :: eta_f16
   real                      :: eta_c_f16
   real                      :: asal_bar
   real                      :: hgt_max_trop
   real, dimension(2)        :: c14f15_bs_xx
   !----- Constants shared by both bdead and bleaf (tropical PFTs) ------------------------!
   real                  , parameter :: a1          =  -1.981
   real                  , parameter :: b1          =   1.047
   real                  , parameter :: dcrit       = 100.0
   !----- Constants used by bdead only (tropical PFTs) ------------------------------------!
   real                  , parameter :: c1d         =   0.572
   real                  , parameter :: d1d         =   0.931
   real                  , parameter :: a2d         =  -1.086
   real                  , parameter :: b2d         =   0.876
   real                  , parameter :: c2d         =   0.604
   real                  , parameter :: d2d         =   0.871
   !----- Constants used by bleaf only (tropical PFTs) ------------------------------------!
   real                  , parameter :: c1l         =  -0.584
   real                  , parameter :: d1l         =   0.550
   real                  , parameter :: a2l         =  -4.111
   real                  , parameter :: b2l         =   0.605
   real                  , parameter :: c2l         =   0.848
   real                  , parameter :: d2l         =   0.438
   !----- Constants shared by both bdead and bleaf (tropical PFTs) ------------------------!
   !---------------------------------------------------------------------------------------!
   !     MLO.   These are the new parameters obtained by adjusting a curve that is similar !
   !            to the modified Chave's equation to include wood density effect on the     !
   !            DBH->AGB allometry as described by:                                        !
   !                                                                                       !
   !            Baker, T. R., and co-authors, 2004: Variation in wood density determines   !
   !               spatial patterns in Amazonian forest biomass.  Glob. Change Biol., 10,  !
   !               545-562.                                                                !
   !                                                                                       !
   !     The "a" parameters were obtaining by splitting balive and bdead at the same ratio !
   ! as the original ED-2.1 allometry, and optimising a function of the form               !
   ! B? = (rho / a3) * exp [a1 + a2 * ln(DBH)]                                             !
   !     The "z" parameters were obtaining by using the original balive and computing      !
   ! bdead as the difference between the total biomass and the original balive.            !
   !---------------------------------------------------------------------------------------!
   real, dimension(3)    , parameter :: ndead_small = (/-1.2639530, 2.4323610, 1.8018010 /)
   real, dimension(3)    , parameter :: ndead_large = (/-0.8346805, 2.4255736, 2.6822805 /)
   real, dimension(3)    , parameter :: nleaf       = (/ 0.0192512, 0.9749494, 2.5858509 /)
   real, dimension(2)    , parameter :: ncrown_area = (/ 0.1184295, 1.0521197            /)
   !---------------------------------------------------------------------------------------!
   !   Coefficients for leaf and structural biomass (iallom = 3 or 5).  For adult          !
   ! individuals, we use the pantropical allometric equation from C14 that estimates AGB   !
   ! and the leaf biomass from an allometric equation derived from F15 data (tropical      !
   ! forest, wild flowering trees only), and the size- and site-dependent stratified       !
   ! sampling and aggregation (J17).  Total individual leaf area was fitted, so to get     !
   ! biomass we must divide by SLA.  The C2B term is added here but is removed when the    !
   ! coefficients are set.                                                                 !
   !                                                                                       !
   !  References:                                                                          !
   !                                                                                       !
   !   Chave, J, Rejou-Mechain M, Burquez A, Chidumayo E, Colgan MS, Delitti WB, Duque A,  !
   !      Eid T, Fearnside PM, Goodman RC et al. 2014. Improved allometric models to       !
   !      estimate the aboveground biomass of tropical trees. Glob. Change Biol., 20(10),  !
   !      3177-3190. doi:10.1111/gcb.12629 (C14).                                          !
   !                                                                                       !
   !   Falster DS, Duursma RA, Ishihara MI, Barneche DR, FitzJohn RG, Vahammar A, Aiba M,  !
   !      Ando M, Anten N, Aspinwall MJ. 2015. BAAD: a biomass and allometry database for  !
   !      woody plants. Ecology, 96 (5):1445-1445. doi:10.1890/14-1889.1 (F15).            !
   !                                                                                       !
   !   Jucker T, Caspersen J, Chave J, Antin C, Barbier N, Bongers F, Dalponte M,          !
   !      van Ewijk KY, Forrester DI, Haeni M et al. 2017. Allometric equations for        !
   !      integrating remote sensing imagery into forest monitoring programmes.            !
   !      Glob. Change Biol., 23(1):177-190. doi:10.1111/gcb.13388 (J17).                  !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   real, dimension(2)    , parameter :: c14f15_bl_xx  = (/ 0.23384770,0.6410495 /)
   real, dimension(3)    , parameter :: c14f15_la_wd  = (/-0.5874,0.5679,0.5476 /)
   real, dimension(3)    , parameter :: c14f15_ht_xx  = (/0.5709,-0.1007,0.6734 /)
   real, dimension(2)    , parameter :: c14f15_bs_tf  = (/ 0.06080334,1.0044785 /)
   real, dimension(2)    , parameter :: c14f15_bs_sv  = (/ 0.05602791,1.0093501 /)
   real, dimension(2)    , parameter :: c14f15_bs_gr  = (/ 1.0e-5, 1.0 /) * c14f15_bl_xx
   real                  , parameter :: SLA_ref       = 17.419
   real                  , parameter :: rho_ref       = 0.610
   !---------------------------------------------------------------------------------------!


   !----- Carbon-to-biomass ratio of plant tissues. ---------------------------------------!
   C2B    = 2.0
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Liana-specific parameters, move here so they are properly initialised.            !
   ! MLO - Manfredo, is there a reason two define two dbh_crit for lianas? Couldn't this   !
   !       number simply replace dbh_crit for lianas?                                      !
   !---------------------------------------------------------------------------------------!
   h_edge = 0.5          !< maximum height advantage for lianas
   liana_dbh_crit = 26.0 !< liana specific critical dbh
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Wood density.                                                                     !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      if (is_grass(ipft)) then ! Grasses. Dummy value for some WAI estimates.
         rho(ipft) = 0.08
      elseif (is_liana(ipft)) then ! BCI traits
         rho(ipft) = 0.46
      elseif (is_tropical(ipft) .and. is_conifer(ipft)) then ! Sub-tropical conifers
         rho(ipft) = 0.52 ! From TRY
      elseif ((.not. is_tropical(ipft)) .and. (.not. is_grass(ipft))) then
         !---------------------------------------------------------------------------------!
         !   Mid-latitude PFTs (not used in the model, defined for completeness by XX).    !
         !---------------------------------------------------------------------------------!
         select case (ipft)
         case ( 6) ! Northern pines, use white pine - ponderosa pine
            rho(ipft) = 0.45
         case ( 7) ! Southern pines, use loblolly pine
            rho(ipft) = 0.57
         case ( 8) ! Late-successional conifers, use black spruce
            rho(ipft) = 0.45
         case ( 9) ! Early-successional hardwood, use aspen
            rho(ipft) = 0.42
         case (10) ! Mid-succesional hardwood, use red oak
            rho(ipft) = 0.74
         case (11) ! Late-successional hardwood, use sugar maple
            rho(ipft) = 0.70
         end select
         !---------------------------------------------------------------------------------!
      else
         !----- Tropical broadleaf trees.  These must be defined individually. ------------!
         select case (economics_scheme)
         case (1)
            !------------------------------------------------------------------------------!
            !     Test: use TRY+GLOPNET data base and cluster analysis to define PFTs.     !
            !------------------------------------------------------------------------------!
            select case (ipft)
            case ( 1, 5,16) ! Grasses
               rho(ipft) = 0.080
            case (    2,12) ! Early-successional tropical/savannah.
               rho(ipft) = 0.450
            case (    3,13) ! Mid-successional tropical/savannah.
               rho(ipft) = 0.615
            case (    4,14) ! Late-successional tropical/savannah.
               rho(ipft) = 0.790
            case default ! Just in case some PFT was forgotten, use global average
               rho(ipft) = rho_ref
            end select
            !------------------------------------------------------------------------------!
         case default
            select case (ipft)
            case ( 1, 5,16) ! Grasses
               rho(ipft) = 0.20
            case (    2,12)  ! Early-successional tropical/savannah.
               rho(ipft) = 0.53 ! 0.40
            case (    3,13)  ! Mid-successional tropical/savannah.
               rho(ipft) = 0.71 ! 0.60
            case (    4,14)  ! Late-successional tropical/savannah.
               rho(ipft) = 0.90 ! 0.87
            case default ! Just in case some PFT was forgotten, use global average
               rho(ipft) = rho_ref
            end select
         end select
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Leaf turnover rate.  We initialise it here despite not being allocation because   !
   ! this parameter may be needed to define SLA.                                           !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      if (is_liana(ipft)) then ! Lianas
         leaf_turnover_rate(ipft) = 1.27
      elseif (is_tropical(ipft) .and. is_conifer(ipft)) then ! Sub-tropical conifers
         leaf_turnover_rate(ipft) = 0.04160842 ! From TRY
      elseif (is_conifer(ipft)) then ! Temperate conifers
         leaf_turnover_rate(ipft) = onethird
      elseif (is_grass(ipft) .and. (.not. is_tropical(ipft))) then ! Temperate grasses
         leaf_turnover_rate(ipft) = 2.0
      elseif (.not. is_tropical(ipft)) then ! Hardwoods. Phenology drives turnover.
         leaf_turnover_rate(ipft) = 0.0
      elseif (economics_scheme == 1) then
         !---------------------------------------------------------------------------------!
         !      Trait trade-off and cluster analysis uses SLA because of the abundance of  !
         ! SLA data in the TRY+GLOPNET data bases.  Here we provide the leaf turnover rate !
         ! based on SLA + standard major axis, so SLA will be the same as the original.    !
         ! This also ensures SLA can be predicted from dynamic leaf turnover rate.         !
         !                                                                                 !
         ! We must set case by case.                                                       !
         !---------------------------------------------------------------------------------!
         select case (ipft)
         case (1,16)     ! C4 grass
            leaf_turnover_rate(ipft) = 2.0
         case (2,12)     ! Early-successional tropical/savannah.
            leaf_turnover_rate(ipft) = 1.55581030
         case (3,13)     ! Mid-successional tropical/savannah.
            leaf_turnover_rate(ipft) = 0.80621772
         case (4,14)     ! Late-successional tropical/savannah.
            leaf_turnover_rate(ipft) = 0.40146228
         case default ! Just in case
            leaf_turnover_rate(ipft) = 1.3141913
         end select
         !---------------------------------------------------------------------------------!
      else
         !---- Tropical plants, we must set case by case. ---------------------------------!
         select case (ipft)
         case (1,16)  ! Grasses
            leaf_turnover_rate(ipft) = 2.00
         case (2,12)  ! Early- successional
            leaf_turnover_rate(ipft) = 1.00
         case (3,13)  ! Mid-successional
            leaf_turnover_rate(ipft) = 0.50
         case (4,14)  ! Late-successional
            leaf_turnover_rate(ipft) = onethird
         case default ! Just in case
            leaf_turnover_rate(ipft) = 0.50
         end select
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Set specific leaf area (SLA, m2leaf/kgC).  The curve relating SLA and leaf        !
   ! turnover rate for ECONOMICS_SCHEME = 1 came a combination of multiple data sets (W04, !
   ! C09, K11, B17, N17).  The model fitting for ECONOMICS_SCHEME = 0 was developed by     !
   ! K12.                                                                                  !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   ! Bahar, N. H. A. et al. Leaf-level photosynthetic capacity in lowland Amazonian and    !
   !    high-elevation Andean tropical moist forests of Peru. New Phytol.,                 !
   !    214(3):1002-1018, May 2017. doi:10.1111/nph.14079 (B17).                           !
   !                                                                                       !
   ! Chave, J., D. Coomes, S. Jansen, S. L. Lewis, N. G. Swenson, and A. E. Zanne. Towards !
   !    a worldwide wood economics spectrum. Ecol. Lett., 12(4):351-366, Apr 2009.         !
   !    doi:10.1111/j.1461-0248.2009.01285.x (C09).                                        !
   !                                                                                       !
   ! Kattge, J., S. Diaz, S. Lavorel, et al., TRY -- a global database of plant traits.    !
   !    Glob. Change Biol., 17 (9): 2905-2935, Sep 2011.                                   !
   !    doi:10.1111/j.1365-2486.2011.02451.x (K11).                                        !
   !                                                                                       !
   ! Kim, Y., R. G. Knox, M. Longo, D. Medvigy, L. R. Hutyra, E. H. Pyle, S. C. Wofsy,     !
   !    R. L. Bras, and P. R. Moorcroft. Seasonal carbon dynamics and water fluxes in an   !
   !    Amazon rainforest. Glob. Change Biol., 18 (4):1322-1334, Apr 2012.                 !
   !    doi:10.1111/j.1365-2486.2011.02629.x (K12).                                        !
   !                                                                                       !
   ! Norby, R. J., L. Gu, I. C. Haworth, A. M. Jensen, B. L. Turner, A. P. Walker,         !
   !    J. M. Warren, D. J. Weston, C. Xu, and K. Winter. Informing models through         !
   !    empirical relationships between foliar phosphorus, nitrogen and photosynthesis     !
   !    across diverse woody species in tropical forests of Panama.                        !
   !    New Phytol., 215 (4):1425-1437, Sep 2017. doi:10.1111/nph.14319 (N17).             !
   !                                                                                       !
   ! Wright, I. J., P. B. Reich, M. Westoby, et al., The worldwide leaf economics          !
   !    spectrum. Nature, 428(6985):821-827, Apr 2004. doi:10.1038/nature02403 (W04).      !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      if (is_tropical(ipft) .and. is_conifer(ipft)) then
         !----- Sub-tropical conifers, use median from TRY. -------------------------------!
         sla_s0(ipft) = 6.324555
         sla_s1(ipft) = 0.0
         !---------------------------------------------------------------------------------!
      elseif (is_tropical(ipft)) then
         !----- Tropical trees, check economics_scheme. -----------------------------------!
         select case (economics_scheme)
         case (1)
            !----- Standard Major Axis derived from TRY/GLOPNET/RAINFOR/NGEE-Tropics. -----!
            if (is_grass(ipft)) then ! Grasses, use trait data base
               sla_s0(ipft) = 15.159000
               sla_s1(ipft) = 0.8637294
            else ! Broadleaf trees, use trait data base
               sla_s0(ipft) = 21.667770
               sla_s1(ipft) = 0.4314902
            end if
            !------------------------------------------------------------------------------!
         case default
            !----- Original ED-2.1 scheme. ------------------------------------------------!
            if (is_grass(ipft)) then
               sla_s0(ipft) = 22.7
               sla_s1(ipft) = 0.0
            else ! Broadleaf trees, use trait data base
               sla_s0(ipft) = exp(log(.1*C2B) + 2.4 * log(10.) - 0.46 * log(12.))
               sla_s1(ipft) = 0.46
            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !    Temperate trees, each PFT must be initialised separately.                    !
         !---------------------------------------------------------------------------------!
         select case (ipft)
         case (5)     ! Temperate C3 grass.
            sla_s0(ipft) = 22.0
         case (6)     ! Northern pines. 
            sla_s0(ipft) = 6.0
         case (7)     ! Southern pines.
            sla_s0(ipft) = 9.0
         case (8)     ! Late conifers. 
            sla_s0(ipft) = 10.0
         case (9)     ! Early hardwood.
            sla_s0(ipft) = 30.0
         case (10)    ! Mid hardwood. 
            sla_s0(ipft) = 24.2
         case (11)    ! Late hardwood. 
            sla_s0(ipft) = 60.0 ! Possible unit conversion issue? In ED-1 it used to be 30.
         case default ! Just in case. 
            sla_s0(ipft) = SLA_ref
         end select
         !---------------------------------------------------------------------------------!

         !----- Ensure SLA = sla_s0. ------------------------------------------------------!
         sla_s1(ipft) = 0.0
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!


   !----- Apply turnover:SLA relationship (but check for zero turnover to avoid FPE). -----!
   SLA(:) = merge( sla_s0(:)                                                               &
                 , sla_s0(:) * leaf_turnover_rate(:)**sla_s1(:)                            &
                 , leaf_turnover_rate(:)*sla_s1(:) == 0.        )
   !---------------------------------------------------------------------------------------!




   !----- Fraction of structural stem that is assumed to be above ground. -----------------!
   agf_bs(:) = 0.7
   !---------------------------------------------------------------------------------------!


   !----- Ratio between fine roots and leaves [kg_fine_roots/kg_leaves] -------------------!
   q(:) = merge(1.0,merge(0.3463,1.1274,is_conifer(:)),is_tropical(:) .or. is_grass(:))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   KIM: ED1/ED2 codes and Moorcroft et al. had the incorrect ratio.                    !
   !   MLO: The ratio is corrected only for tropical PFTs using iallom = 3 or 5.  To       !
   !        extend this fix to other PFTs, one must refit parameters for other tissues     !
   !        (e.g. bdead), so the total AGB is consistent with the original allometric      !
   !        equation for AGB.                                                              !
   !                                                                                       !
   !        For the PFTs that were updated, we combine the pipe model with the data from   !
   !        CA08 and shape parameter from F16 to derive the ratio.                         !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   ! Calvo-Alvarado, J. C., N. G. McDowell, and R. H. Waring. Allometric relationships     !
   !    predicting foliar biomass and leaf area:sapwood area ratio from tree height in     !
   !    five Costa Rican rain forest species. Tree Physiol., 28 (11):1601-1608, Sep 2008.  !
   !    doi:10.1093/treephys/28.11.1601. (CA08)                                            !
   !                                                                                       !
   ! Falster, D. S., R. G. FitzJohn, A. Brannstrom, U. Dieckmann, and M. Westoby.  plant:  !
   !    A package for modelling forest trait ecology and evolution.  Methods Ecol. Evol.,  !
   !    7(2):136-146, Feb 2016. doi:10.1111/2041-210X.12525. (F16)                         !
   !                                                                                       !
   ! McDowell, N., H. Barnard, B. Bond, T. Hinckley, R. Hubbard, H. Ishii, B. Kostner,     !
   !    F. Magnani, J. Marshall, F. Meinzer, N. Phillips, M. Ryan, and D. Whitehead. The   !
   !    relationship between tree height and leaf area: sapwood area ratio. Oecologia,     !
   !    132(1):12-20, Jun 2002. doi:10.1007/s00442-002-0904-x. (MD02)                      !
   !                                                                                       !
   ! Rosell, J. A., S. Gleason, R. Mendez-Alonzo, Y. Chang, and M. Westoby. Bark           !
   !    functional ecology: evidence for tradeoffs, functional coordination, and environ-  !
   !    ment producing bark diversity. New Phytol., 201(2): 486-497, Jan 2014.             !
   !    doi:10.1111/nph.12541. (R14)                                                       !
   !                                                                                       !
   ! Yokozawa M., and T. Hara. Foliage profile, size structure and stem diameter-plant     !
   !    height relationship in crowded plant populations. Ann. Bot.-London, 76(3):271-285, !
   !    Sep 1995. doi:10.1006/anbo.1995.1096. (YH95)                                       !
   !                                                                                       !
   !      Leaf-to-sapwood area ratio (Al:As, or As/Al) for conifers was estimated from the !
   ! average of all conifers listed in Table 1 of MD02, weighted by number of individuals. !
   ! Broadleaf is currently tropical-only, and was obtained from the average of all points !
   ! from Fig. 1 of CA08 (obtained from extracting data from the figure itself).           !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   select case (iallom)
   case (3,5)
      do ipft=1,n_pft
         if (is_liana(ipft)) then
            !------------------------------------------------------------------------------!
            !      Lianas.  For the time being I did not change to preserve total biomass. !
            ! Mind that this number may be off by one or two orders of magnitude.          !
            !------------------------------------------------------------------------------!
            sapwood_ratio(ipft) = 3900.0
            !------------------------------------------------------------------------------!
         elseif (is_tropical(ipft) .and. is_grass(ipft)) then
            !------------------------------------------------------------------------------!
            !     Tropical grasses.  Only a small fraction should be sapwood.              !
            !------------------------------------------------------------------------------!
            sapwood_ratio(ipft) = SLA(ipft) / 1.0e-5
            !------------------------------------------------------------------------------!
         elseif (is_tropical(ipft)) then
            !------------------------------------------------------------------------------!
            !     Define the shape parameter eta.  Eta describes the vertical distribution !
            ! of leaf area within the crown.  According to F16 and the original reference  !
            ! (YH95), eta=1 is typical of conifers, whereas eta=12 is closer to the        !
            ! profiles observed in angiosperms.  Araucarias somewhat intermediate, so we   !
            ! set eta=5.                                                                   !
            !                                                                              !
            !     asal_bar is the leaf-to-sapwood area ratio (Al:As, or As/Al) for         !
            ! conifers was estimated from the average of all conifers listed in Table 1 of !
            ! MD02, weighted by number of individuals.  Broadleaf is currently tropical-   !
            ! only, and was obtained from the average of all points from Fig. 1 of CA08    !
            ! (obtained from extracting data from the figure itself)                       !
            !------------------------------------------------------------------------------!
            if (is_conifer(ipft)) then ! Sub-tropical needleleaf
               eta_f16  = 5.0
               asal_bar = 3.709184e-05
            else ! Broadleaf plant (trees/grasses).
               eta_f16  = 12.0
               asal_bar = 7.400139e-05
            end if
            !------------------------------------------------------------------------------!


            !------ Eta_c_f16 is calculated following F16. --------------------------------!
            eta_c_f16 = 1.0 - 2.0 / (1.0+eta_f16) + 1.0 / (1 + 2.0 * eta_f16)
            !------------------------------------------------------------------------------!


            !------ Sapwood ratio. --------------------------------------------------------!
            sapwood_ratio(ipft) = 1.0 / ( eta_c_f16 * asal_bar * rho(ipft) * 1000. / C2B )
            !------------------------------------------------------------------------------!

         else
            !------------------------------------------------------------------------------!
            !      Default ED-1 ratio.  Mind that this number may be off by one or two     !
            ! orders of magnitude.                                                         !
            !------------------------------------------------------------------------------!
            sapwood_ratio(ipft) = 3900.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!

   case default
      !------------------------------------------------------------------------------------!
      !      Default ED-1 ratio.  Mind that this number may be off by one or two orders of !
      ! magnitude...                                                                       !
      !------------------------------------------------------------------------------------!
      sapwood_ratio(:) = 3900.0
      !------------------------------------------------------------------------------------!

   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    SLA-dependent ratio between sapwood and leaves [kg_sapwood/kg_leaves]              !
   !---------------------------------------------------------------------------------------!
   qsw(:)    = SLA(:) / sapwood_ratio(:)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Specific Root Area, based on:                                                      !
   !                                                                                       !
   ! Metcalfe, D. B., P. Meir, L. E. O. C. Aragao, A. C. L. da Costa, A. P. Braga,         !
   !    P. H. L. Goncalves, J. d. A. Silva Junior, S. S. de Almeida, L. A. Dawson,         !
   !    Y. Malhi, and M. Williams (2008), The effects of water availability on root growth !
   !    and morphology in an Amazon rainforest, Plant Soil, 311(1-2), 189-199,             !
   !    doi:10.1007/s11104-008-9670-9.                                                     !
   !---------------------------------------------------------------------------------------!
   SRA(:)   = 24. * 2. ! m2/kgC --> this is from Amazon
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Root vertical distribution shape factor, based on:                                 !
   !                                                                                       !
   ! Jackson, R. B., J. Canadell, J. R. Ehleringer, H. A. Mooney, O. E. Sala, and          !
   !    E. D. Schulze (1996), A global analysis of root distributions for terrestrial      !
   !    biomes, Oecologia, 108(3), 389-411, doi:10.1007/BF00333714.                        !
   !                                                                                       !
   !    The reference only reports for ecosystem level root profile but shows that root    !
   ! cumulative biomass generally follows an exponential distribution.  Here, we are       !
   ! applying this pattern to cohort level, ignoring the potential root niche separation   !
   ! that results in smaller fraction of roots in top layers relative to deeper layers/    !
   !                                                                                       !
   !    It is assumed that the root has an exponential distribution, with only ROOT_BETA   !
   ! fraction of roots below maximum rooting depth. Increasing the ROOT_BETA will increase !
   ! the root profile fraction from deeper layers.                                         !
   !                                                                                       !
   !    The root fraction (Y) above depth D cm for a cohort with max rooting depth as      !
   !  D_max (cm) can be calculated as:                                                     !
   !                                                                                       !
   !  Y = ( 1. - (root_beta) ** (D / D_max) ) / (1 - root_beta)                            !
   !                                                                                       !
   !                                                                                       !
   !  MLO (2020-10-27): I added the denominator (1 - root_beta) to ensure that Y at        !
   !                    D=D_max is always 1, regardless of the value of root_beta, as      !
   !                    long as root_beta < 1.                                             !
   !                                                                                       !
   !    Suggested values range from 0.0001 to 0.1.                                         !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   root_beta(:)   =   0.1
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !   qrhob is the ratio between bark density and wood density.  For tropical broadleaf   !
   ! trees, we use fitted curves using SMA on the data available in P14 (Tab. S1).  For    !
   ! other PFTs, we use the average ratio based on R14, which is similar to the average    !
   ! of all data points of P14, but it has more types of ecosystems.                       !
   !                                                                                       !
   ! Poorter, L., A. McNeil, V.-H. Hurtado, H. H. T. Prins, and F. E. Putz. Bark traits    !
   !    and life-history strategies of tropical dry- and moist forest trees.               !
   !    Funct. Ecol., 28(1):232-242, Feb 2014. doi:10.1111/1365-2435.12158 (P14)           !
   !                                                                                       !
   ! Rosell, J. A., S. Gleason, R. Mendez-Alonzo, Y. Chang, and M. Westoby. Bark           !
   !    functional ecology: evidence for tradeoffs, functional coordination, and environ-  !
   !    ment producing bark diversity. New Phytol., 201(2): 486-497, Jan 2014.             !
   !    doi:10.1111/nph.12541. (R14)                                                       !
   !---------------------------------------------------------------------------------------!
   qrhob(:) = merge( 0.49 / 0.61                                                           &
                   , exp(0.6966550 - 1.602123 * rho(:))                                    &
                   , is_grass(:) .or. is_conifer(:) .or. (.not. is_tropical(:)) )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Set bark thickness and carbon allocation to bark.  This is currently done only    !
   ! for tropical trees when IALLOM = 3 or 5, because all biomass pools must be corrected  !
   ! to ensure that total aboveground biomass is consistent with the allometric equations. !
   ! This may and should be changed in the future.                                         !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   ! Meinzer, F. C., G. Goldstein, and J. L. Andrade. Regulation of water flux through     !
   !    tropical forest canopy trees: Do universal rules apply? Tree Physiol.,             !
   !    21(1):19-26, Jan 2001. doi:10.1093/treephys/21.1.19. (M01)                         !
   !                                                                                       !
   ! de Mattos, P. P., A. T. dos Santos, H. Rivera, Y. M. M. de Oliveira, M. A. D. Rosot,  !
   !    and M. C. Garrastazu. Growth of Araucaria angustifolia in the Embrapa/Epagri       !
   !    forest reserve, Cacador, SC, Brazil. Pesq. Flor. Bras., 55(2):107-114, Jul 2007.   !
   !    URL http://pfb.cnpf.embrapa.br/pfb/index.php/pfb/ article/view/124. In Portuguese. !
   !    (M07)                                                                              !
   !                                                                                       !
   ! Lawes, M. J. , J. J. Midgley, and P. J. Clarke. Costs and benefits of relative bark   !
   !    thickness in relation to fire damage: a savanna/forest contrast. J. Ecol.,         !
   !    101(2):517-524, Dec 2013. doi:10.1111/1365-2745.12035. (L13)                       !
   !                                                                                       !
   ! Falster, D. S., R. G. FitzJohn, A. Brannstrom, U. Dieckmann, and M. Westoby.          !
   !    plant: A package for modelling forest trait ecology and evolution.  Methods Ecol.  !
   !    Evol., 7(2):136-146, Feb 2016. doi:10.1111/2041-210X.12525. (F16)                  !
   !                                                                                       !
   ! b1Xs - slope of the dbh to sapwood thickness curve.  Currently this is used only to   !
   !        define biomass allocation to bark consistent with bark thickness.              !
   !                                                                                       !
   ! b1Xb - slope of the dbh to bark thickness curve.                                      !
   ! qbark - ratio between leaf biomass and bark biomass per unit height.                  !
   !---------------------------------------------------------------------------------------!
   select case (iallom)
   case (3,5)
      !------ New allometry, use estimate based on M01. -----------------------------------!
      b1Xs(:) = 0.315769481
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Define bark thickness parameter based on life form.                            !
      !------------------------------------------------------------------------------------!
      do ipft=1,n_pft
         if (.not. is_tropical(ipft)) then ! Non-tropical trees (bark is not set)
            b1Xb(ipft) = 0.0
         elseif (is_grass(ipft)) then ! Grasses (bark is not set)
            b1Xb(ipft) = 0.0
         elseif (is_liana(ipft)) then ! Lianas (bark is not set)
            b1Xb(ipft) = 0.0
         elseif (is_conifer (ipft)) then ! Araucarias, use slope from M07, Table 1
            b1Xb(ipft) = 0.03936468
         elseif (is_savannah(ipft)) then ! Avg. Slope of savannah trees (L13)
            b1Xb(ipft) = 0.128
         else ! Avg. Slope of forest trees (L13)
            b1Xb(ipft) = 0.019
         end if
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Define the area ratio between sapwood and bark, assuming that sapwood and bark !
      ! are concentric rings.                                                              !
      !------------------------------------------------------------------------------------!
      abas(:)  = b1Xb(:) * (1.0 - b1Xb(:)) / ( b1Xs(:) * (1. + b1Xs(:) - 2. * b1Xb(:)) )
      qbark(:) = qrhob(:) * abas(:) * qsw(:)
      !------------------------------------------------------------------------------------!

   case default
      !------ Old allometry, exclude bark. ------------------------------------------------!
      b1Xs (:) = 0.315769481
      b1Xb (:) = 0.0
      abas (:) = 0.0
      qbark(:) = 0.0
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   DBH/height allometry parameters.                                                    !
   !                                                                                       !
   !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!    !
   !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!    !
   !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!    !
   !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!    !
   !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!    !
   !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!    !
   !                                                                                       !
   !   b1Ht, b2Ht, and hgt_ref are parameters that have different meaning for tropical and !
   ! temperate PFTs, and the meaning for tropical PFTs depends on IALLOM.                  !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      if (is_liana(ipft)) then
         !---- Liana allometry: 10cm lianas are 35m tall ----------------------------------!
         hgt_ref(ipft) = 61.7
         b1Ht   (ipft) = 0.1136442
         b2Ht   (ipft) = 0.8675
         !---------------------------------------------------------------------------------!
      elseif (is_tropical(ipft)) then
         !---------------------------------------------------------------------------------!
         !    Tropical trees have different parameters and functional forms depending on   !
         ! the allometry.                                                                  !
         !---------------------------------------------------------------------------------!
         select case (iallom)
         case (0:1)
            !----- Regular log-log fit, b1 is the intercept and b2 is the slope. ----------!
            b1Ht   (ipft) = 0.37 * log(10.0)
            b2Ht   (ipft) = 0.64
            !----- hgt_ref is not used. ---------------------------------------------------!
            hgt_ref(ipft) = 0.0
            !------------------------------------------------------------------------------!
         case (2)
            !------------------------------------------------------------------------------!
            ! Weibull function --- H = Hinf * (1-exp(-b1*DBH^b2)) --- proposed by:         !
            !                                                                              !
            ! Poorter, L., L. Bongers, and F. Bongers. Architecture of 54 moist-forest     !
            !    tree species: traits, trade-offs, and functional groups. Ecology,         !
            !    87(5):1289-1301, May 2006.                                                !
            !    doi:10.1890/0012- 9658(2006)87[1289:AOMTST]2.0.CO;2.                      !
            !------------------------------------------------------------------------------!
            !----- b1Ht is their "a" and b2Ht is their "b". -------------------------------!
            b1Ht   (ipft) = 0.0352
            b2Ht   (ipft) = 0.694
            !----- hgt_ref is their "Hmax". -----------------------------------------------!
            hgt_ref(ipft) = 61.7
            !------------------------------------------------------------------------------!
        case (3,5)
            !------------------------------------------------------------------------------!
            !     Allometric equation based on the fitted curve using the Sustainable      !
            ! Landscapes data set (L16) and the size- and site-dependent stratified        !
            ! sampling and aggregation (J17), as described in (L20).  This relationship is !
            ! fitted using Standardised Major Axis (SMA) so the same parameters can be     !
            ! used for y=f(x) and x=f(y).  This is particularly useful when initialising   !
            ! the model with airborne lidar data (L20).  Because it would be extremely     !
            ! cumbersome to derive a SMA-based regression based on Weibull function, we    !
            ! use a log-linear relationship.  The maximum height is based on the 99%       !
            ! quantile of all trees measured by the SL team.                               !
            !                                                                              !
            ! References:                                                                  !
            !                                                                              !
            ! Jucker T, Caspersen J, Chave J, Antin C, Barbier N, Bongers F, Dalponte M,   !
            !    van Ewijk KY, Forrester DI, Haeni M et al. 2017. Allometric equations for !
            !    integrating remote sensing imagery into forest monitoring programmes.     !
            !    Glob. Change Biol., 23(1):177-190. doi:10.1111/gcb.13388 (J17).           !
            !                                                                              !
            ! Longo M, Keller M, dos-Santos MN, Leitold V, Pinage ER, Baccini A,           !
            !    Saatchi S, Nogueira EM, Batistella M , Morton DC. 2016. Aboveground       !
            !    biomass variability across intact and degraded forests in the Brazilian   !
            !    Amazon.  Global Biogeochem. Cycles, 30(11):1639-1660.                     !
            !    doi:10.1002/2016GB005465 (L16).                                           !
            !                                                                              !
            ! Longo M, Saatchi SS, Keller M, Bowman KW, Ferraz A, Moorcroft PR, Morton D,  !
            !    Bonal D, Brando P, Burban B et al. 2020. Impacts of degradation on water, !
            !    energy, and carbon cycling of the Amazon tropical forests. J. Geophys.    !
            !    Res.-Biogeosci., 125: e2020JG005677. doi:10.1029/2020JG005677 (L20).      !
            !------------------------------------------------------------------------------!
            b1Ht   (ipft) = 1.139963
            b2Ht   (ipft) = 0.564899
            !----- hgt_ref is not used. ---------------------------------------------------!
            hgt_ref(ipft) = 0.0
            !------------------------------------------------------------------------------!
         case (4)
            !------------------------------------------------------------------------------!
            !  Allometric equation based on Chave et al. 2014 and Falster et al. 2015      !
            !------------------------------------------------------------------------------!
            b1Ht   (ipft) = c14f15_ht_xx(1) + c14f15_ht_xx(2) * log(rho(ipft))
            b2Ht   (ipft) = c14f15_ht_xx(3)
            !----- hgt_ref is not used. ---------------------------------------------------!
            hgt_ref(ipft) = 0.0
         end select
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Temperate PFTs, each one has a different model fitting.  Reference:         !
         !                                                                                 !
         ! Albani, M., D. Medvigy, G. C. Hurtt, and P. R. Moorcroft. The contributions of  !
         !    land-use change, CO2 fertilization, and climate variability to the eastern   !
         !    US carbon sink. Glob. Change Biol., 12(12):2370-2390, Dec 2006.              !
         !    doi:10.1111/j.1365-2486.2006.01254.x.                                        !
         !---------------------------------------------------------------------------------!
         select case (ipft)
         case (5)     ! Temperate C3 grass. 
            b1Ht   (ipft) =  0.4778
            b2Ht   (ipft) = -0.75
            hgt_ref(ipft) =  0.0
         case (6,7)   ! Northern and Southern Pines. 
            b1Ht   (ipft) = 27.14
            b2Ht   (ipft) = -0.03884
            hgt_ref(ipft) = 1.3
         case (8)     ! Late conifers. 
            b1Ht   (ipft) = 22.79
            b2Ht   (ipft) = -0.04445 
            hgt_ref(ipft) = 1.3
         case (9)     ! Early hardwood. 
            b1Ht   (ipft) = 22.6799
            b2Ht   (ipft) = -0.06534
            hgt_ref(ipft) = 1.3
         case (10)    ! Mid hardwood. 
            b1Ht   (ipft) = 25.18
            b2Ht   (ipft) = -0.04964
            hgt_ref(ipft) = 1.3
         case (11)    ! Late hardwood. 
            b1Ht   (ipft) = 23.3874
            b2Ht   (ipft) = -0.05404
            hgt_ref(ipft) = 1.3
         case default ! Forgotten PFT.
            b1Ht   (ipft) = 25.18
            b2Ht   (ipft) = -0.04964
            hgt_ref(ipft) = 1.3
         end select
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Minimum and maximum height allowed for each cohort.                                !
   !---------------------------------------------------------------------------------------!
   select case (iallom)
   case (3,4,5)
      !------------------------------------------------------------------------------------!
      !    This value corresponds to the 99% quantile of all trees measured by the         !
      ! Sustainable Landscapes.                                                            !
      !------------------------------------------------------------------------------------!
      hgt_max_trop = 46.0
      !------------------------------------------------------------------------------------!
   case default
      hgt_max_trop = 35.0
   end select
   hgt_min(:) = merge( merge(0.50        ,0.50         ,is_grass(:))                       &
                     , merge(0.15        ,hgt_ref+0.2  ,is_grass(:))                       &
                     , is_tropical(:)                                )
   hgt_max(:) = merge( merge(1.50        ,hgt_max_trop ,is_grass(:))                       &
                     , merge(0.95*b1Ht(:),0.999*b1Ht(:),is_grass(:))                       &
                     , is_tropical(:)                                )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   MIN_DBH     -- minimum DBH allowed for the PFT.                                     !
   !   DBH_CRIT    -- minimum DBH that brings the PFT to its tallest possible height.      !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      min_dbh    (ipft) = h2dbh(hgt_min(ipft),ipft)
      dbh_crit   (ipft) = h2dbh(hgt_max(ipft),ipft)
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    This is the typical DBH that all big leaf plants will have.  Because the big-leaf  !
   ! ED doesn't really solve individuals, the typical DBH should be one that makes a good  !
   ! ratio between LAI and biomass.  This is a tuning parameter and right now the initial  !
   ! guess is about 1/3 of the critical DBH for trees.                                     !
   !---------------------------------------------------------------------------------------!
   select case (iallom)
   case (0,1)
      !----- Critical DBH for all PFTs. ---------------------------------------------------!
      dbh_bigleaf(:) = dbh_crit(:)
      !------------------------------------------------------------------------------------!
   case default
      do ipft=1,n_pft
         if (is_grass(ipft)) then
            !----- Grasses: critical DBH. -------------------------------------------------!
            dbh_bigleaf(ipft) = dbh_crit(ipft)
            !------------------------------------------------------------------------------!
         else
            !----- Trees: 1/3 of the critical DBH. ----------------------------------------!
            dbh_bigleaf(ipft) = dbh_crit(ipft) * onethird
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     DBH-crown allometry.                                                              !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      if (is_liana(ipft)) then
         !----- Lianas. -------------------------------------------------------------------!
         b1Ca(ipft) = exp(ncrown_area(1))
         b2Ca(ipft) = 1.26254364
         !---------------------------------------------------------------------------------!
      elseif (is_tropical(ipft)) then
         !----- Tropical PFT.  Decide parameters based on iallom. -------------------------!
         select case (iallom)
         case (0,1)
            !------------------------------------------------------------------------------!
            !    Coefficients from Poorter et al. (2006), after being transformed so the   !
            ! allometry is a function of DBH.                                              !
            !                                                                              !
            ! Poorter, L., L. Bongers, and F. Bongers. Architecture of 54 moist-forest     !
            !    tree species: traits, trade-offs, and functional groups. Ecology,         !
            !    87(5):1289-1301, May 2006.                                                !
            !    doi:10.1890/0012- 9658(2006)87[1289:AOMTST]2.0.CO;2.                      !
            !------------------------------------------------------------------------------!
            b1Ca(ipft) = exp(-1.853) * exp(b1Ht(ipft)) ** 1.888
            b2Ca(ipft) = b2Ht(ipft) * 1.888
            !------------------------------------------------------------------------------!
         case (2)
            !------------------------------------------------------------------------------!
            !    Coefficients also based on Poorter et al. (2006), but using refitted      !
            ! coefficients to use DBH as the predictive variable.                          !
            !------------------------------------------------------------------------------!
            b1Ca(ipft) = exp(ncrown_area(1))
            b2Ca(ipft) = ncrown_area(2)
            !------------------------------------------------------------------------------!
         case (3,4,5)
            !------------------------------------------------------------------------------!
            !     Allometry using the Sustainable Landscapes data.                         !
            !------------------------------------------------------------------------------!
            !                                                                              !
            ! Longo M, Saatchi SS, Keller M, Bowman KW, Ferraz A, Moorcroft PR, Morton D,  !
            !    Bonal D, Brando P, Burban B et al. 2020. Impacts of degradation on water, !
            !    energy, and carbon cycling of the Amazon tropical forests. J. Geophys.    !
            !    Res.-Biogeosci., 125: e2020JG005677. doi:10.1029/2020JG005677 (L20).      !
            !                                                                              !
            !    Equation was derived from forest inventory measurements carried out at    !
            ! multiple locations in the Brazilian Amazon, and fitted using a               !
            ! heteroscedastic least squares approach.                                      !
            !                                                                              !
            !     The functional form uses both DBH and Height.                            !
            !                                                                              !
            !                       CA = b1Ca * (DBH^2 * Hgt)^b2Ca                         !
            !                       m2            cm      m                                !
            !                                                                              !
            ! Total number of trees: 17072                                                 !
            ! b1Ca    = 0.370 (95% CI: [0.346;0.398])                                      !
            ! b2Ca    = 0.464 (95% CI: [0.457;0.472])                                      !
            ! R2      = 0.521                                                              !
            ! RMSE    = 29.78                                                              !
            !------------------------------------------------------------------------------!
            b1Ca(ipft) = 0.370
            b2Ca(ipft) = 0.464
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Temperate PFT, use Dietze and Clark (2008).                                 !
         !---------------------------------------------------------------------------------!
         b1Ca(ipft) = 2.490154
         b2Ca(ipft) = 0.8068806
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These are used to compute the crown length, which will then be used to find the   !
   ! height of the bottom of the crown.  This allometry is based on:                       !
   !                                                                                       !
   ! Poorter L., L. Bongers, F. Bongers, 2006: Architecture of 54 moist-forest tree        !
   !     species: traits, trade-offs, and functional groups. Ecology, 87, 1289-1301.       !
   !                                                                                       !
   !    For iallom = 3 or 5, we use the allometric equation based on the Sustainable       !
   ! Landscapes data set.                                                                  !
   !                                                                                       !
   ! Longo M, Saatchi SS, Keller M, Bowman KW, Ferraz A, Moorcroft PR, Morton D, Bonal D,  !
   !    Brando P, Burban B et al. 2020. Impacts of degradation on water, energy, and       !
   !    carbon cycling of the Amazon tropical forests. J. Geophys. Res.-Biogeosci., 125:   !
   !    e2020JG005677. doi:10.1029/2020JG005677 (L20).                                     !
   !                                                                                       !
   !    Equation was derived from forest inventory measurements carried out at multiple    !
   ! locations in the Brazilian Amazon, and fitted using a heteroscedastic least           !
   ! squares approach.                                                                     !
   !                                                                                       !
   ! Total number of trees: 16064                                                          !
   ! b1Cl    = 0.298 (95% CI: [0.288;0.306])                                               !
   ! b2Cl    = 1.032 (95% CI: [1.022;1.044])                                               !
   ! R2      = 0.673                                                                       !
   ! RMSE    = 2.29                                                                        !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      !----- Check life form. -------------------------------------------------------------!
      if (is_grass(ipft)) then
         b1Cl(ipft) = 0.99
         b2Cl(ipft) = 1.00
      elseif (is_tropical(ipft)) then
         !----- Tropical PFTs: check allometry settings. ----------------------------------!
         select case (iallom)
         case (3,4,5)
            b1Cl(ipft) = 0.29754
            b2Cl(ipft) = 1.0324
         case default
            b1Cl(ipft) = 0.3106775
            b2Cl(ipft) = 1.098
         end select
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Temperate PFTs: right now we are using the same allometry as the old        !
         ! tropical (Poorter et al. 2006).  We probably need to switch these coefficients  !
         ! by appropriate reference.                                                       !
         !---------------------------------------------------------------------------------!
         b1Cl(ipft) = 0.3106775
         b2Cl(ipft) = 1.098
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters for DBH -> Bleaf allometry.                                            !
   !                                                                                       !
   !   IALLOM = 0,1,2  --  Bleaf = b1Bl * DBH^b2Bl                                         !
   !   IALLOM = 3,4,5  --  leaf_A= b1Bl * (DBH*DBH*Height)^b2Bl                            !
   !                       b1Bl is a function of wood density (IALLOM=4 only).             !
   !                       For IALLOM=3,4,5, leaf biomass will depend on SLA.              !
   !                                                                                       !
   !   The coefficients and thresholds depend on the PFT and allometric equations.  In     !
   ! addition to the coefficients, we define the dbh point that defines adult cohorts as   !
   ! opposed to seedlings, and the associated leaf biomass.                                !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      if (is_liana(ipft)) then
          !----- Liana leaf Biomass (Putz, 1983) ------------------------------------------!
          b1Bl (ipft) = 0.0856
          b2Bl (ipft) = 2.0
          !--------------------------------------------------------------------------------!
      elseif (is_tropical(ipft)) then
         select case(iallom)
         case (0,1)
            !------------------------------------------------------------------------------!
            !      ED-2.0 allometry, based on:                                             !
            !                                                                              !
            !   Saldarriaga, J. G., D. C. West, M. L. Tharp, and C. Uhl.  Long-term        !
            !      chronosequence of forest succession in the upper Rio Negro of Colombia  !
            !      and Venezuela.  J. Ecol., 76, 4, 938-958, 1988.                         !
            !------------------------------------------------------------------------------!
            b1Bl (ipft) = exp(a1 + c1l * b1Ht(ipft) + d1l * log(rho(ipft)))
            aux         = ( (a2l - a1) + b1Ht(ipft) * (c2l - c1l) + log(rho(ipft))   &
                        * (d2l - d1l) ) * (1.0/log(dcrit))
            b2Bl (ipft) = C2B * b2l + c2l * b2Ht(ipft) + aux
            !------------------------------------------------------------------------------!
         case (2)
            !------------------------------------------------------------------------------!
            !     ED-2.1 allometry, based on:                                              !
            !                                                                              !
            !   Calvo-Alvarado, J. C., N. G. McDowell, and R. H. Waring.  Tree Physiol.,   !
            !      28, 11, 1601-1608, 2008.                                                !
            !                                                                              !
            !   Cole, T. G., J. J. Ewel.  Allometric equations for four valuable tropical  !
            !      tree species.  Forest Ecol. Manag., 229, 1--3, 351-360, 2006.           !
            !------------------------------------------------------------------------------!
            b1Bl (ipft) = C2B * exp(nleaf(1)) * rho(ipft) / nleaf(3)
            b2Bl (ipft) = nleaf(2)
            !------------------------------------------------------------------------------!
        case (3,5) 
            !------------------------------------------------------------------------------!
            !    Allometry based on the BAAD data based (F15) and described in (L20).  We  !
            ! only used leaves from wild tropical, flowering trees, and applied a          !
            ! stratified sample by DBH class and location and cross-validation, following  !
            ! (J17).                                                                       !
            !                                                                              !
            ! References:                                                                  !
            !                                                                              !
            ! Falster DS, Duursma RA, Ishihara MI, Barneche DR, FitzJohn RG, Vahammar A,   !
            !    Aiba M, Ando M, Anten N, Aspinwall MJ. 2015. BAAD: a biomass and          !
            !    allometry database for woody plants. Ecology, 96 (5):1445-1445.           !
            !    doi:10.1890/14-1889.1 (F15).                                              !
            !                                                                              !
            ! Jucker T, Caspersen J, Chave J, Antin C, Barbier N, Bongers F, Dalponte M,   !
            !    van Ewijk KY, Forrester DI, Haeni M et al. 2017. Allometric equations for !
            !    integrating remote sensing imagery into forest monitoring programmes.     !
            !    Glob. Change Biol., 23(1):177-190. doi:10.1111/gcb.13388 (J17).           !
            !                                                                              !
            ! Longo M, Saatchi SS, Keller M, Bowman KW, Ferraz A, Moorcroft PR, Morton D,  !
            !    Bonal D, Brando P, Burban B et al. 2020. Impacts of degradation on water, !
            !    energy, and carbon cycling of the Amazon tropical forests. J. Geophys.    !
            !    Res.-Biogeosci., 125: e2020JG005677. doi:10.1029/2020JG005677 (L20).      !
            !------------------------------------------------------------------------------!
            b1Bl(ipft) = c14f15_bl_xx(1)
            b2Bl(ipft) = c14f15_bl_xx(2)
            !------------------------------------------------------------------------------!
        case (4)
            !------------------------------------------------------------------------------!
            !    Allometry based on the BAAD data based (F15).  We only used leaves from   !
            ! wild tropical, note that b1Bl has the unit of m2 leaf under this scenario    !
            ! and will be converted to leaf carbon using SLA in size2bl.                   !
            !                                                                              !
            ! Falster DS, Duursma RA, Ishihara MI, Barneche DR, FitzJohn RG, Vahammar A,   !
            !    Aiba M, Ando M, Anten N, Aspinwall MJ. 2015. BAAD: a biomass and          !
            !    allometry database for woody plants. Ecology, 96 (5):1445-1445.           !
            !    doi:10.1890/14-1889.1 (F15).                                              !
            !------------------------------------------------------------------------------!
            b1Bl(ipft) = exp( c14f15_la_wd(1) + c14f15_la_wd(2) * log(rho(ipft)))
            b2Bl(ipft) = c14f15_la_wd(3)
            !------------------------------------------------------------------------------!

         end select
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !      Temperate PFTs.  Each class has specific optimised parameters, and there   !
         ! is no distinction between small and large cohorts.                              !
         !---------------------------------------------------------------------------------!
         select case (ipft)
         case (5)   ! Temperate C3 grass. 
            b1Bl(ipft) = 0.08
            b2Bl(ipft) = 1.0
         case (6,7) ! Northern and Southern Pines. 
            b1Bl(ipft) = 0.024
            b2Bl(ipft) = 1.899
         case (8)   ! Late conifers. 
            b1Bl(ipft) = 0.0454
            b2Bl(ipft) = 1.6829
         case (9)   ! Early hardwood. 
            b1Bl(ipft) = 0.0129
            b2Bl(ipft) = 1.7477
         case (10) ! Mid hardwood. 
            b1Bl(ipft) = 0.048
            b2Bl(ipft) = 1.455
         case (11) ! Late hardwood. 
            b1Bl(ipft) = 0.017
            b2Bl(ipft) = 1.731
         case default ! Just in case
            b1Bl(ipft) = 0.046
            b2Bl(ipft) = 1.930
         end select
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters for DBH -> Bdead allometry.                                            !
   !                                                                                       !
   !   IALLOM = 0, 1, 2                                                                    !
   !                                                                                       !
   !           { b1Bs_small * DBH^b2Bs_small  , if dbh < dbh_crit                          !
   !   Bdead = {                                                                           !
   !           { b1Bs_large * DBH^b2Bl_large  , if dbh > dbh_crit                          !
   !                                                                                       !
   !   IALLOM = 3, 4, 5                                                                    !
   !                                                                                       !
   !   Bdead = b1Bs_small * (DBH^2 * Height) ^ b2Bs_small                                  !
   !                                                                                       !
   !   The coefficients and thresholds depend on the PFT and allometric equations.         !
   !---------------------------------------------------------------------------------------!
   !------- Fill in the tropical PFTs, which are functions of wood density. ---------------!
   do ipft=1,n_pft
      if (is_liana(ipft)) then
         !---------------------------------------------------------------------------------!
         ! Liana allometry                                                                 !
         !                                                                                 !
         !    This is an experimental value for dead biomass. The Schnitzer article        !
         ! provide an allometric equation to relate AGB to DBH. Since here we are dealing  !
         ! with DB instead of AGB I also subtracted the leaf biomass from AGB. Bl dbh      !
         ! allometry is the same used in size2bl from Putz. I have dropped the intercept   !
         ! though. This is because  otherwise one would have DB != when plant has DBH = 0. !
         ! At this point one would have                                                    !
         ! dbh2bd = (exp(-1.484 + 2.657 * log(dbh)) - 0.0856 * dbh * dbh) / C2B            !
         ! Moreover since Schnitzer gives AGB we have to divide by agf_bs so that          !
         ! when we calculate the above ground fraction multiplying by agf_bs we get the    !
         ! correct AGB. I fitted this equation dbh2bd with a form dbh2bd = a*(dbh)**b      !
         ! because otherwise the formula is not invertible (we need bd2dbh). Fit converges !
         ! nicely with                                                                     !
         ! a= 0.13745 and b=2.69373 .                                                      !
         !                                                                                 !
         ! MLO -- qsw is not zero for lianas, so you may be overestimating total AGB in    !
         !        ED-2, because bdeada = agb - bleaf - bsapwooda - bbarka.                 !
         !---------------------------------------------------------------------------------!
         b1Bs_small(ipft)  = 0.2749
         b1Bs_large(ipft)  = b1Bs_small(ipft)
         b2Bs_small(ipft)  = 2.69373
         b2Bs_large(ipft)  = b2Bs_small(ipft)
         !---------------------------------------------------------------------------------!
      elseif (is_tropical(ipft)) then
         select case (iallom)
         case (0)
            !---- ED-2.1 allometry. -------------------------------------------------------!
            b1Bs_small(ipft)  = exp(a1 + c1d * b1Ht(ipft) + d1d * log(rho(ipft)))
            b1Bs_large(ipft)  = exp(a1 + c1d * log(hgt_max(ipft)) + d1d * log(rho(ipft)))

            aux               = ( (a2d - a1) + b1Ht(ipft) * (c2d - c1d) + log(rho(ipft))   &
                                * (d2d - d1d)) * (1.0/log(dcrit))
            b2Bs_small(ipft)  = C2B * b2d + c2d * b2Ht(ipft) + aux

            aux               = ( (a2d - a1) + log(hgt_max(ipft)) * (c2d - c1d)            &
                                + log(rho(ipft)) * (d2d - d1d)) * (1.0/log(dcrit))
            b2Bs_large(ipft)  = C2B * b2d + aux
            !------------------------------------------------------------------------------!

         case (1,2)
            !------------------------------------------------------------------------------!
            !     Based on modified Chave et al. (2001) allometry.                         !
            !                                                                              !
            ! Chave, J., B. Riera, M.-A. Dubois. Estimation of biomass in a neotropical    !
            !    forest of French Guiana: spatial and temporal variability.                !
            !    J. Trop. Ecol., 17(1):79-96, Jan 2001. doi:10.1017/S0266467401001055.     !
            !------------------------------------------------------------------------------!
            b1Bs_small (ipft) = C2B * exp(ndead_small(1)) * rho(ipft) / ndead_small(3)
            b2Bs_small (ipft) = ndead_small(2)
            b1Bs_large (ipft) = C2B * exp(ndead_large(1)) * rho(ipft) / ndead_large(3)
            b2Bs_large (ipft) = ndead_large(2)
            !------------------------------------------------------------------------------!
         case (3,4,5)
            !------------------------------------------------------------------------------!
            ! Trees:   set parameters based on Chave et al. (2014).                        !
            ! Grasses: set numbers to small values, too keep bdead at a minimum but still  !
            !          growing as grasses grow (so dbh can be inferred from bdead).        !
            !                                                                              !
            ! Chave, J.,M. Rejou-Mechain, A. Burquez, et al. Improved allometric models    !
            !    to estimate the aboveground biomass of tropical trees. Glob. Change       !
            !    Biol., 20(10):3177-3190, Oct 2014. doi:10.1111/gcb.12629.                 !
            !------------------------------------------------------------------------------!
            if (is_grass(ipft)) then
               c14f15_bs_xx = c14f15_bs_gr
            elseif (is_savannah(ipft)) then
               c14f15_bs_xx = c14f15_bs_sv
            else
               c14f15_bs_xx = c14f15_bs_tf
            end if
            b1Bs_small(ipft) = c14f15_bs_xx(1) * rho(ipft) ** c14f15_bs_xx(2)
            b2Bs_small(ipft) = c14f15_bs_xx(2)
            b1Bs_large(ipft) = b1Bs_small(ipft)
            b2Bs_large(ipft) = b2Bs_small(ipft)
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !      Temperate PFTs.  Each class has specific optimised parameters, and there   !
         ! is no distinction between small and large cohorts.                              !
         !---------------------------------------------------------------------------------!
         select case (ipft)
         case (5)   ! Temperate C3 grass. 
            b1Bs_small(ipft) = 1.0e-5
            b2Bs_small(ipft) = 1.0
         case (6,7) ! Northern and Southern Pines. 
            b1Bs_small(ipft) = 0.147
            b2Bs_small(ipft) = 2.238
         case (8)   ! Late conifers. 
            b1Bs_small(ipft) = 0.1617
            b2Bs_small(ipft) = 2.1536
         case (9)   ! Early hardwood. 
            b1Bs_small(ipft) = 0.02648
            b2Bs_small(ipft) = 2.95954
         case (10)  ! Mid hardwood. 
            b1Bs_small(ipft) = 0.1617
            b2Bs_small(ipft) = 2.4572
         case (11)  ! Late hardwood. 
            b1Bs_small(ipft) = 0.235
            b2Bs_small(ipft) = 2.2518
         end select
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Duplicate coefficients.                                                      !
         !---------------------------------------------------------------------------------!
         b1Bs_large (ipft) = b1Bs_small(ipft)
         b2Bs_large (ipft) = b2Bs_small(ipft)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     In case we run big leaf model with IALLOM set to 0 or 1, we must change some of   !
   ! the allometric parameters.                                                            !
   !---------------------------------------------------------------------------------------!
   if (ibigleaf == 1 .and. (iallom == 0 .or. iallom == 1)) then
      b1Bl      ( 1) = 0.04538826
      b1Bl      ( 2) = 0.07322115
      b1Bl      ( 3) = 0.07583497
      b1Bl      ( 4) = 0.08915847
      b1Bl      (12) = b1Bl(2)
      b1Bl      (13) = b1Bl(3)
      b1Bl      (14) = b1Bl(4)
      b1Bl      (15) = 0.07322115
      b1Bl      (16) = 0.04538826
      b1Bl      (17) = 0.07322115
      
      b2Bl      ( 1) = 1.316338
      b2Bl      ( 2) = 1.509083
      b2Bl      ( 3) = 1.646576
      b2Bl      ( 4) = 1.663773
      b2Bl      (12) = b2Bl(2)
      b2Bl      (13) = b2Bl(3)
      b2Bl      (14) = b2Bl(4)
      b2Bl      (15) = 1.509083
      b2Bl      (16) = 1.316338
      b2Bl      (17) = 1.509083
      
      b1Bs_small( 1) = 0.05291854
      b1Bs_small( 2) = 0.15940854
      b1Bs_small( 3) = 0.21445616
      b1Bs_small( 4) = 0.26890751
      b1Bs_small(12) = b1Bs_small(2)
      b1Bs_small(13) = b1Bs_small(3)
      b1Bs_small(14) = b1Bs_small(4)
      b1Bs_small(15) = 0.15940854
      b1Bs_small(16) = 0.05291854
      b1Bs_small(17) = 0.15940854
      
      b2Bs_small( 1) = 3.706955
      b2Bs_small( 2) = 2.342587
      b2Bs_small( 3) = 2.370640
      b2Bs_small( 4) = 2.254336
      b2Bs_small(12) = b2Bs_small(2)
      b2Bs_small(13) = b2Bs_small(3)
      b2Bs_small(14) = b2Bs_small(4)
      b2Bs_small(15) = 2.342587
      b2Bs_small(16) = 3.706955
      b2Bs_small(17) = 2.342587
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Fill in variables that are derived from bdead allometry.                          !
   !---------------------------------------------------------------------------------------!
   do ipft = 1, n_pft
      !------------------------------------------------------------------------------------!
      ! -- MIN_BDEAD is the minimum structural biomass possible.  This is used in the      !
      !    initialisation only, to prevent cohorts to be less than the minimum size due to !
      !    change in allometry.                                                            !
      ! -- BDEAD_CRIT corresponds to BDEAD when DBH is exactly at DBH_CRIT.  This is       !
      !    used to determine which b1Bs/b2Bs pair to use.                                  !
      !------------------------------------------------------------------------------------!
      min_bdead (ipft) = size2bd(min_dbh (ipft),hgt_min(ipft),ipft)
      bdead_crit(ipft) = size2bd(dbh_crit(ipft),hgt_max(ipft),ipft)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    WAI parameters, the choice depends on IALLOM.                                      !
   !---------------------------------------------------------------------------------------!
   select case (iallom)
   case (3,4,5)
      !------------------------------------------------------------------------------------!
      !    WAI is defined as a fraction of (potential) LAI.   This is just a refit of      !
      ! allometry 2 but using DBH*DBH*Height as predictor for consistency.                 !
      !------------------------------------------------------------------------------------!
      b1WAI(:) = merge( 0.0                                                                &
                      , merge( merge(0.01148449,0.00378399,is_conifer(:))                  &
                             , merge(0.0553*0.5,0.0192*0.5,is_conifer(:))                  &
                             , is_tropical(:) .and. (.not. is_liana(:) )  )                &
                      , is_grass(:)                                       )
      b2WAI(:) = merge( 1.0                                                                &
                      , merge( merge(0.77075160,0.81667933,is_conifer(:))                  &
                             , merge(    1.9769,    2.0947,is_conifer(:))                  &
                             , is_tropical(:) .and. (.not. is_liana(:) )  )                &
                      , is_grass(:)                                       )
      !------------------------------------------------------------------------------------!
   case default
      !------------------------------------------------------------------------------------!
      !    Use the equation by:                                                            !
      !                                                                                    !
      ! Hormann, G., S. Irrgan, H. Jochheim, M. Lukes, H. Meesenburg, J. Muller,           !
      !    B. Scheler, J. Scherzer, G. Schuler, B. Schultze, B. Strohbach, F. Suckow,      !
      !    M. Wegehenkel, and G. Wessolek.   Wasserhaushalt von waldokosystemen:           !
      !    methodenleitfaden zur bestimmung der wasserhaushaltskomponenten auf level       !
      !    II-flachen. Technical note, Bundesministerium fur Verbraucherschutz, Ernahrung  !
      !    und Landwirtschaft (BMVEL), Bonn, Germany, 2003.                                !
      !    URL http://www.wasklim.de/download/Methodenband.pdf.                            !
      !------------------------------------------------------------------------------------!
      b1WAI(:) = merge(0.0 ,merge(0.0553*0.5,0.0192*0.5,is_conifer(:)),is_grass(:))
      b2WAI(:) = merge(1.0 ,merge(    1.9769,    2.0947,is_conifer(:)),is_grass(:))
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    brf_wd is the fraction of above-ground wood that is assumed to be in branches and  !
   ! twigs.  Here we must be careful to make sure that the fraction is 0 in case WAI is    !
   ! going to be zero (e.g. grasses).                                                      !
   !---------------------------------------------------------------------------------------!
   brf_wd(:) = merge(0.00,0.16,is_grass(:))
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Sapwood area allometry for plants                                                  !
   !    - Grasses: 100% of sapwood area over basal area but sapwood area (not used)        !
   !    - Temperate angiosperms: 30% of sapwood (oaks)                                     !
   !    - Tropical angiosperms:  (M01).                                                    !
   !    - Conifers: 30% of sapwood. (loblolly pine)                                        !
   !    Users are welcome to update those numbuers based on literatures or measurements    !
   !                                                                                       !
   ! Meinzer, F. C., G. Goldstein, and J. L. Andrade (2001), Regulation of water flux      !
   !    through tropical forest canopy trees: Do universal rules apply?, Tree Physiol.,    !
   !    21(1), 19-26, doi:10.1093/treephys/21.1.19 (M01).                                  !
   !---------------------------------------------------------------------------------------!
   select case (iallom)
   case (4)
       ! Data from BAAD and Christoffersen et al. 2016 GMD
       b1SA(:) = merge( 1.0                                                                &
                      , merge( 0.30, merge(0.6572,0.30,is_tropical(:)), is_conifer(:) )     &
                      , is_grass(:)                                                   )
       b2SA(:) = merge(1.8530,2.0,is_tropical(:) .and. (.not. is_grass(:)))  
   case default
       b1SA(:) = merge( 1.0                                                                &
                      , merge( 0.30, merge(1.582,0.30,is_tropical(:)), is_conifer(:) )     &
                      , is_grass(:)                                                   )
       b2SA(:) = merge(1.764,2.0,is_tropical(:) .and. (.not. is_grass(:)))  
   end select
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !     Commercial volume of trees (stem/bole volume, in m3).  The current equation is a  !
   ! re-fit from Nogueira et al. (2008) so a single set of parameters can be used.  Their  !
   ! equation is for tropical trees only, so temperate and boreal forests may need a       !
   ! different set of parameters.  Grasses are assumed to have no commercial volume.       !
   !                                                                                       !
   ! Nogueira, E. M., et al. Estimates of forest biomass in the Brazilian Amazon: new      !
   !    allometric equations and adjustments to biomass from wood-volume inventories.      !
   !    Forest Ecol. Manag., 256(11), 1853-1867, Nov. 2008,                                !
   !    doi:10.1016/j.foreco.2008.07.022.                                                  !
   !---------------------------------------------------------------------------------------!
   b1Vol(:) = merge(0.0,3.528e-5,is_grass(:))
   b2Vol(:) = merge(1.0,0.976   ,is_grass(:))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     DBH-Root depth allometry.  Check which allometry to use.  The original volume     !
   ! equation in ED-2.1 didn't have any references (at least I never found it), so the     !
   ! rooting depth equation now incorporates the original "volume".  The current volume    !
   ! allometry is based on tropical trees, and it estimates the commercial volume (stem    !
   ! volume).                                                                              !
   !---------------------------------------------------------------------------------------!
   select case (iallom)
   case (0)
      !------------------------------------------------------------------------------------!
      !      Grasses are assumed to have constant rooting depth.  Tree rooting depth is a  !
      ! function of dbh and height.                                                        !
      !------------------------------------------------------------------------------------!
      b1Rd(:) = merge( -0.700                                                              &
                     , - exp(0.545*log(10.)) * (0.65 * pi1 * 0.11 * 0.11)**0.277           &
                     , is_grass(:) )
      b2Rd(:) = merge( 0.000,0.277,is_grass(:))
      !------------------------------------------------------------------------------------!
   case (1,2)
      !------------------------------------------------------------------------------------!
      !     Simple fit, based on that the soil  moisture depletion during the dry season   !
      ! in Tapajos is most noticeable down to 4-6 m (C13, p. 187, Fig. 3).  These          !
      ! coefficients make rooting depth to be 0.5m for recruits and 5m for 35-m tall       !
      ! trees.                                                                             !
      !                                                                                    !
      ! References:                                                                        !
      !                                                                                    !
      ! Christoffersen BO. 2013. The ecohydrological mechanisms of resilience and          !
      !    vulnerability of Amazonian tropical forests to water stress. Ph.d.              !
      !    dissertation, University of Arizona, Tucson, AZ, USA.                           !
      !    URL http://hdl.handle.net/10150/293566 (C13).                                   !
      !------------------------------------------------------------------------------------!
      b1Rd(:)  = -1.1140580
      b2Rd(:)  =  0.4223014
      !------------------------------------------------------------------------------------!
    case (3)
      !------------------------------------------------------------------------------------!
      !    Test allometry, similar to 2, but based on D*D*H.  The curve loosely fits B18   !
      ! for large trees, and Xiangtao's fit based on unpublished data from excavation in   !
      ! Panama (mostly small trees).  This is not a full, formal optimisation because the  !
      ! data from Panama were not available and the data from B18 is indirect, so there is !
      ! plenty of room for improvement.  The curve is exactly the same as iallom=2 for     !
      ! non-tropical trees.                                                                !
      !                                                                                    !
      ! References:                                                                        !
      !                                                                                    !
      ! Brum M, Vadeboncoeur MA, Ivanov V, Asbjornsen H, Saleska S, Alves LF,              !
      !    Penha D, Dias JD, Aragao LEOC, Barros F, Bittencourt P, Pereira L,              !
      !    Oliveira RS, 2018.  Hydrological niche segregation defines forest               !
      !    structure and drought tolerance strategies in a seasonal Amazonian              !
      !    forest. J. Ecol., in press. doi:10.1111/1365-2745.13022 (B18).                  !
      !------------------------------------------------------------------------------------!
      b1Rd(:) = merge( merge(-2.572,-0.947,is_grass(:))                                    &
                     , -1.1140580                                                          &
                     , is_tropical(:) .and. (.not. is_liana(:)) )
      b2Rd(:) = merge( merge(0.246,0.148,is_grass(:))                                      &
                     , +0.4223014                                                          &
                     , is_tropical(:) .and. (.not. is_liana(:)) )
      !------------------------------------------------------------------------------------!
   case (4,5)
      !------------------------------------------------------------------------------------!
      !    Test allometry based on excavation data in Panama based on  H.                  !
      !    Multiply it by 2 so that a 40 m tree can get access to water below 5m depth     !
      !------------------------------------------------------------------------------------!
      b1Rd(:) = -0.609 * 2.
      b2Rd(:) = 0.580

   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Alternative rooting depth parameter for tropical trees.  Use the allometric model  !
   ! to obtain the Effective Functional Rooting Depth based on B18.  We made a slight      !
   ! modification in their equation relating delta 18O and depth:                          !
   !                                                                                       !
   !    depth = exp(a + b * d18O^2)                                                        !
   !                                                                                       !
   ! because it fits the data better than the original equation without the square, and it !
   ! avoids extremely shallow soils for small trees.  We also use a heteroscedastic least  !
   ! squares, using the algorithm developed by L16.                                        !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   ! Brum M, Vadeboncoeur MA, Ivanov V, Asbjornsen H, Saleska S, Alves LF, Penha D,        !
   !    Dias JD, Aragao LEOC, Barros F, Bittencourt P, Pereira L, Oliveira RS, 2018.       !
   !    Hydrological niche segregation defines forest structure and drought tolerance      !
   !    strategies in a seasonal Amazonian forest. J. Ecol., in press.                     !
   !    doi:10.1111/1365-2745.13022 (B18).                                                 !
   !---------------------------------------------------------------------------------------!
   d18O_ref(:) = -5.356
   b1d18O  (:) = 0.0516
   b2d18O  (:) = 1.0
   b1Efrd  (:) = -2.436
   b2Efrd  (:) = 0.1822
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Initial storage pool, relative to on-allometry living biomass, to be given to the  !
   ! PFTs when the model is run using INITIAL conditions.                                  !
   !---------------------------------------------------------------------------------------!
   select case (ivegt_dynamics)
   case (0)
      f_bstorage_init(:) = onesixth
   case default
      f_bstorage_init(:) = 1.00
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Initial density of plants, for near-bare-ground simulations [# of individuals/m2]  !
   !---------------------------------------------------------------------------------------!
   select case (ibigleaf)
   case (0)
      !----- Size and age structure. ------------------------------------------------------!
      select case (iallom)
      case (0,1)
         init_density(:) = merge(1.0,0.1,is_grass(:))
      case default
         init_density(:) = 0.1
      end select
      !------------------------------------------------------------------------------------!

      !----- Define a non-sense number. ---------------------------------------------------!
      init_laimax(:)   = huge_num
      !------------------------------------------------------------------------------------!

   case(1)
      !----- Big leaf. 1st we set the maximum initial LAI for each PFT. -------------------!
      init_laimax(:)   = 0.1
      do ipft=1,n_pft
         init_bleaf         = size2bl(dbh_bigleaf(ipft),hgt_max(ipft),SLA(ipft),ipft)
         init_density(ipft) = init_laimax(ipft) / (init_bleaf * SLA(ipft))
      end do
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!



   !----- Leaf width [m].  This controls the boundary layer conductance. ------------------!
   leaf_width(:) = merge( lwidth_grass                                                     &
                        , merge(lwidth_nltree,lwidth_bltree,is_conifer(:))                 &
                        , is_grass(:) )
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Characteristic branch diameter [m].  This controls the boundary layer             !
   ! conductance.  Currently we assume the same size as leaves, although we could use      !
   ! branch size distribution functions to estimate it dynamically.                        !
   !---------------------------------------------------------------------------------------!
   branch_diam(:) = leaf_width(:)
   !---------------------------------------------------------------------------------------!


   !----- Define the number of bins for the look-up tables. -------------------------------!
   nbt_lut = 10000
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Parameters for the trait-plasticity model.  These control SLA, for the controls   !
   ! on Vm0, check init_pft_photo_params.                                                  !
   !---------------------------------------------------------------------------------------!
   !     This controls the expansion/extinction factor for SLA when                        !
   ! TRAIT_PLASTICITY_SCHEME > 0.  The default depends upon the reference SLA and is       !
   ! initialised in init_derived_params_after_xml.  This relationship is empirical, so we  !
   !  allow it to be initialised through XML too.                                          !
   !---------------------------------------------------------------------------------------!
   kplastic_SLA (:) = undef_real
   !---------------------------------------------------------------------------------------!
   !     This controls the expansion/reduction exponent for leaf turnover rate when using  !
   ! trait plasticity when TRAIT_PLASTICITY_SCHEME is negative.  Currently we start from   !
   ! Eq. 1 of X17 (originally from K91), by substituting b with their SMA fit, and fitting !
   ! a curve relating Aa with Vcmax25_m, which was obtained by simulating ED-2 for         !
   ! multiple average diurnal cycles at tower sites in South America, using multiple       !
   ! Vcmax25_m values found in trait data bases (TRY/GLOPNET/RAINFOR/NGEE-Tropics).        !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   ! Kikuzawa K. 1991. A cost-benefit analysis of leaf habit and leaf longevity of trees   !
   !    and their geographical pattern. Am. Nat. 138(5): 1250-1263. doi:10.1086/285281     !
   !    (K91).                                                                             !
   !                                                                                       !
   ! Xu X, Medvigy D, Wright SJ, Kitajima K, Wu J, Albert LP, Martins GA, Saleska SR,      !
   !    Pacala SW. 2017. Variations of leaf longevity in tropical moist forests predicted  !
   !    by a trait-driven carbon optimality model. Ecol. Lett. 20(9): 1097-1106.           !
   !    doi:10.1111/ele.12804 (X17).                                                       !
   !---------------------------------------------------------------------------------------!
   eplastic_vm0  (:) = 0.5 * (-1.36 - 0.8399)
   eplastic_sla  (:) = 0.5 * (-1.36 - 1.0000)
   !---------------------------------------------------------------------------------------!
   !     This controls the expansion/extinction factor for leaf longevity when             !
   ! TRAIT_PLASTICITY_SCHEME > 0.  The default depends upon the reference leaf turnover    !
   ! rate so it is initialised in init_derived_params_after_xml.  This relationship is     !
   ! empirical, so we allow it to be initialised through XML too.                          !
   !---------------------------------------------------------------------------------------!
   kplastic_LL   (:) = undef_real
   !----- Maximum LAI to consider for plasticity (TRAIT_PLASTICITY_SCHEME is positive). ---!
   laimax_plastic(:) = 6.0
   !----- Linearised slope for when TRAIT_PLASTICITY_SCHEME is negative. ------------------!
   LMA_slope     (:) = 0.015  ! linearized slope (Trait_plasticity_scheme < 0)
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_pft_alloc_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_photo_params()
   use ed_max_dims    , only : n_pft                     & ! intent(in)
                             , undef_real                ! ! intent(in)
   use ed_misc_coms   , only : ibigleaf                  & ! intent(in)
                             , economics_scheme          ! ! intent(in)
   use pft_coms       , only : is_tropical               & ! intent(in)
                             , is_conifer                & ! intent(in)
                             , is_grass                  & ! intent(in)
                             , is_liana                  & ! intent(in)
                             , SLA                       & ! intent(in)
                             , C2B                       & ! intent(in)
                             , D0                        & ! intent(out)
                             , Vm0                       & ! intent(out)
                             , Vm0_v0                    & ! intent(out)
                             , Vm0_v1                    & ! intent(out)
                             , Vm_low_temp               & ! intent(out)
                             , Vm_high_temp              & ! intent(out)
                             , Vm_decay_elow             & ! intent(out)
                             , Vm_decay_ehigh            & ! intent(out)
                             , Vm_hor                    & ! intent(out)
                             , Vm_q10                    & ! intent(out)
                             , kplastic_vm0              & ! intent(out)
                             , Jm0                       & ! intent(out)
                             , Jm_low_temp               & ! intent(out)
                             , Jm_high_temp              & ! intent(out)
                             , Jm_decay_elow             & ! intent(out)
                             , Jm_decay_ehigh            & ! intent(out)
                             , Jm_hor                    & ! intent(out)
                             , Jm_q10                    & ! intent(out)
                             , TPm0                      & ! intent(out)
                             , Rd0                       & ! intent(out)
                             , kplastic_rd0              & ! intent(out)
                             , Rd_low_temp               & ! intent(out)
                             , Rd_high_temp              & ! intent(out)
                             , Rd_decay_elow             & ! intent(out)
                             , Rd_decay_ehigh            & ! intent(out)
                             , Rd_hor                    & ! intent(out)
                             , Rd_q10                    & ! intent(out)
                             , stomatal_slope            & ! intent(out)
                             , cuticular_cond            & ! intent(out)
                             , quantum_efficiency        & ! intent(out)
                             , qyield_psII               & ! intent(out)
                             , curvpar_electron          & ! intent(out)
                             , photosyn_pathway          & ! intent(out)
                             , dark_respiration_factor   & ! intent(out)
                             , electron_transport_factor & ! intent(out)
                             , triose_phosphate_factor   & ! intent(out)
                             , water_conductance         ! ! intent(out)
   use consts_coms    , only : t00                       & ! intent(in)
                             , rmol                      & ! intent(in)
                             , twothirds                 & ! intent(in)
                             , umol_2_mol                & ! intent(in)
                             , yr_sec                    ! ! intent(in)
   use physiology_coms, only : iphysiol                  & ! intent(in)
                             , vmfact_c3                 & ! intent(in)
                             , vmfact_c4                 & ! intent(in)
                             , mphoto_trc3               & ! intent(in)
                             , mphoto_tec3               & ! intent(in)
                             , mphoto_c4                 & ! intent(in)
                             , bphoto_blc3               & ! intent(in)
                             , bphoto_nlc3               & ! intent(in)
                             , bphoto_c4                 & ! intent(in)
                             , gamma_c3                  & ! intent(in)
                             , gamma_c4                  & ! intent(in)
                             , d0_grass                  & ! intent(in)
                             , d0_tree                   & ! intent(in)
                             , alpha_c3                  & ! intent(in)
                             , alpha_c4                  & ! intent(in)
                             , kw_grass                  & ! intent(in)
                             , kw_tree                   & ! intent(in)
                             , q10_c3                    & ! intent(in)
                             , q10_c4                    & ! intent(in)
                             , tphysref                  ! ! intent(in)
   implicit none
   !---------------------------------------------------------------------------------------!


   !----- Local variables. ----------------------------------------------------------------!
   real(kind=4)                   :: ssfact       !> Correction factor for Vcmax
   real(kind=4)                   :: vmexpo_pft   !> Correction factor dor Rdmax
   real(kind=4)                   :: gamma_pft    !> Rdmax / Vcmax
   real(kind=4)                   :: a_pft        !> Parameter for gamma estimation
   real(kind=4)                   :: b_pft        !> Parameter for gamma estimation
   integer                        :: ipft         !> PFT index
   !----- Local parameters, based on Atkin et al. (2015), Table S3 (Rdark,m). -------------!
   real(kind=4), parameter        :: a_c3grss = -1.962            !> 5d, C3H
   real(kind=4), parameter        :: b_c3grss =  1.247            !> 5d, C3H
   real(kind=4), parameter        :: a_c4grss =  0.30103000-1.962 !> Make it twice C3H
   real(kind=4), parameter        :: b_c4grss =  1.247            !> Make it twice C3H
   real(kind=4), parameter        :: a_bltrop = -1.533            !> 5f, TWQ >= 25C
   real(kind=4), parameter        :: b_bltrop =  1.022            !> 5f, TWQ >= 25C
   real(kind=4), parameter        :: a_bltemp = -0.862            !> 5f, 15C <= TWQ < 25C
   real(kind=4), parameter        :: b_bltemp =  0.753            !> 5f, 15C <= TWQ < 25C
   real(kind=4), parameter        :: a_needle = -0.366            !> 5d, NlT
   real(kind=4), parameter        :: b_needle =  0.494            !> 5d, NlT
   !---------------------------------------------------------------------------------------!


   !----- Photosynthetic pathway (3 is C3; 4 is C4). --------------------------------------!
   do ipft=1,n_pft
      if (is_grass(ipft)) then
         select case(ipft)
         case (1) ! C4 grass
            photosyn_pathway(ipft) = 4
         case default ! C3 grasses
            photosyn_pathway(ipft) = 3
         end select
      else ! Trees, C3 photosynthesis
         photosyn_pathway(ipft) = 3
      end if
   end do
   !---------------------------------------------------------------------------------------!


   !----- Critical VPD, for stomata closure due to dry air. -------------------------------!
   D0(:) = merge(d0_grass,d0_tree,is_grass(:))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Temperature thresholds for photosynthetic capacity.                               !
   !---------------------------------------------------------------------------------------!
   select case (iphysiol)
   case (0,2)
      Vm_low_temp (:) = merge(4.7137,10.0,is_conifer(:) .or. (.not. is_tropical(:)))
      Vm_high_temp(:) = 45.0
   case (1,3)
      Vm_low_temp (:) = merge( 4.7137                                                      &
                             , merge(15.0,10.0,photosyn_pathway(:) == 4)                   &
                             , is_conifer(:) .or. (.not. is_tropical(:)) )
      Vm_high_temp(:) = merge(42.5,45.0,photosyn_pathway(:) == 4)
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !>    Vm_decay_elow and Vm_decay_ehigh are the correction terms for high and low 
   !> temperatures when running the original ED-2.1 correction as in Moorcroft et al.
   !> (2001).
   !---------------------------------------------------------------------------------------!
   select case (iphysiol)
   case (0,2)
      Vm_decay_elow (:) = 0.4
      Vm_decay_ehigh(:) = 0.4
   case (1,3)
      Vm_decay_elow (:) = 0.3
      Vm_decay_ehigh(:) = 0.3
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Vm_hor is the Arrhenius "activation energy" divided by the universal gas         !
   ! constant.  Vm_q10 is the base for the Collatz approach.                               !
   !---------------------------------------------------------------------------------------!
   select case (iphysiol)
   case (0,2)
      !------------------------------------------------------------------------------------!
      !  Default parameters (M01/M09/L19).                                                 !
      !                                                                                    !
      ! Longo M, Knox RG, Medvigy DM, Levine NM, Dietze MC, Kim Y, Swann ALS, Zhang K,     !
      !    Rollinson CR, Bras RL et al. 2019. The biophysics, ecology, and biogeochemistry !
      !    of functionally diverse, vertically and horizontally heterogeneous ecosystems:  !
      !    the Ecosystem Demography model, version 2.2 -- part 1: Model description.       !
      !    Geosci. Model Dev., 12: 4309-4346. doi:10.5194/gmd-12-4309-2019 (L19).          !
      !                                                                                    !
      ! Medvigy DM, Wofsy SC, Munger JW, Hollinger DY , Moorcroft PR. 2009. Mechanistic    !
      !    scaling of ecosystem function and dynamics in space and time: Ecosystem         !
      !    demography model version 2. J. Geophys. Res.-Biogeosci., 114: G01002.           !
      !    doi:10.1029/2008JG000812 (M09).                                                 !
      !                                                                                    !
      ! Moorcroft PR, Hurtt GC , Pacala SW. 2001. A method for scaling vegetation          !
      !    dynamics: The Ecosystem Demography model (ED). Ecol. Monogr., 71: 557-586.      !
      !    doi:10.1890/0012- 9615(2001)071[0557:AMFSVD]2.0.CO;2 (M01).                     !
      !------------------------------------------------------------------------------------!
      vm_hor(:) = 3000.
      vm_q10(:) = merge(q10_c4,q10_c3,photosyn_pathway(:) == 4)
      !------------------------------------------------------------------------------------!
   case (1,3)
      !------------------------------------------------------------------------------------!
      !  Use values from vC00.                                                             !
      !                                                                                    !
      ! von Caemmerer S. 2000. Biochemical models of leaf photosynthesis. No. 2 in         !
      !    Techniques in Plant Sciences. CSIRO Publishing, Collingwood, VIC, Australia.    !
      !    doi:10.1006/anbo.2000.1296 (vC00).                                              !
      !------------------------------------------------------------------------------------!
      vm_hor(:) = 58520. * tphysref / (rmol * (t00+25.))
      vm_q10(:) = 2.21
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Vm0 is the maximum photosynthesis capacity in umol/m2/s.  Notice that depending   !
   ! on the size structure (SAS or Big Leaf), there is an addition factor multiplied.      !
   !---------------------------------------------------------------------------------------!
   !----- Find the additional factor to multiply Vm0. -------------------------------------!
   select case (ibigleaf)
   case (0) ! SAS, use only the modification from the namelist. 
      ssfact = 1.0
   case (1)
      ssfact = 3.0
   end select
   !---- Define Vm0 (case by case). -------------------------------------------------------!
   do ipft=1,n_pft
      if (is_tropical(ipft) .and. (.not. is_conifer(ipft)) .and. (.not. is_liana(ipft)) )  &
      then
         select case (economics_scheme)
         case (1)
            !------------------------------------------------------------------------------!
            !     New tropical parameters based on multiple data sets (K11,W04, B17, and   !
            ! N17) and using a standard major axis on log-transformed SLA and              !
            ! photosynthesis traits.  The area-based photosynthesis traits were poorly     !
            ! correlated with SLA, but the mass-based showed stronger correlation.         !
            ! Because ED-2 needs the area-based parameters, converted the relationship     !
            ! to area based by dividing the mass-based Vcmax with SLA (tooking care of     !
            ! the necessary unit conversions).                                             !
            !                                                                              !
            ! References                                                                   !
            !                                                                              !
            !                                                                              !
            ! Bahar, N. H. A. et al. Leaf-level photosynthetic capacity in lowland         !
            !    Amazonian and high-elevation Andean tropical moist forests of Peru. New   !
            !    Phytol., 214(3):1002-1018, May 2017. doi:10.1111/nph.14079 (B17).         !
            !                                                                              !
            ! Kattge, J., S. Diaz, S. Lavorel, et al., TRY -- a global database of plant   !
            !    traits. Glob. Change Biol., 17 (9): 2905-2935, Sep 2011.                  !
            !    doi:10.1111/j.1365-2486.2011.02451.x (K11).                               !
            !                                                                              !
            ! Norby, R. J., L. Gu, I. C. Haworth, A. M. Jensen, B. L. Turner,              !
            !    A. P. Walker, J. M. Warren, D. J. Weston, C. Xu, and K. Winter. Informing !
            !    models through empirical relationships between foliar phosphorus,         !
            !    nitrogen and photosynthesis across diverse woody species in tropical      !
            !    forests of Panama.  New Phytol., 215 (4):1425-1437, Sep 2017.             !
            !    doi:10.1111/nph.14319 (N17).                                              !
            !                                                                              !
            ! Rowland, L. et al. Scaling leaf respiration with nitrogen and phosphorus in  !
            !    tropical forests across two continents. New Phytol., 214(3):1064-1077,    !
            !    May 2017. doi:10.1111/nph.13992 (R17).                                    !
            !                                                                              !
            ! Wright, I. J., P. B. Reich, M. Westoby, et al., The worldwide leaf           !
            !    economics spectrum. Nature, 428(6985):821-827, Apr 2004.                  !
            !    doi:10.1038/nature02403 (W04).                                            !
            !------------------------------------------------------------------------------!
            select case (ipft)
            case (1)     ! C4 grass. Use CLM defaults
               Vm0_v0(ipft) = 23.348416
               Vm0_v1(ipft) = 0.0
            case (16)    ! Subtropical C3 grass.  Use CLM defaults
               Vm0_v0(ipft) = 35.384615
               Vm0_v1(ipft) = 0.0
            case default ! Tropical trees.  Divide by Q10 to bring it to 15C.
               Vm0_v0(ipft) = 9.22326 / vm_q10(ipft)
               Vm0_v1(ipft) = 0.50175
            end select
            !------------------------------------------------------------------------------!
         case default
            !------ Original parameters. --------------------------------------------------!
            select case (ipft)
            case (1)     ! C4 grass. 
               Vm0_v0(ipft) = 12.50
            case (2,12)  ! Early-successional tropical tree
               Vm0_v0(ipft) = 18.75
            case (3,13)  ! Mid-successional tropical tree
               Vm0_v0(ipft) = 12.50
            case (4,14)  ! Late-successional tropical tree
               Vm0_v0(ipft) = 6.25
            case (16)
               Vm0_v0(ipft) = 18.75
            case default !  Just in case. 
               Vm0_v0(ipft) = 15.625
            end select
            Vm0_v1(ipft) = 0.0
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      else
         !------ Temperate PFTs. ----------------------------------------------------------!
         select case (ipft)
         case (5)     ! C3 grass. 
            Vm0_v0(ipft) = 18.300000
         case (6,7)   ! Pines (N/S). 
            Vm0_v0(ipft) = 11.350000
         case (8)     ! Late conifers. 
            Vm0_v0(ipft) = 4.540000
         case (9)     ! Early hardwood. 
            Vm0_v0(ipft) = 20.387075
         case (10)    ! Mid hardwood. 
            Vm0_v0(ipft) = 17.454687
         case (11)    ! Late hardwood.
            Vm0_v0(ipft) = 6.981875
         case (15)    ! Araucaria
            Vm0_v0(ipft) = 10.
         case (17)    ! Liana
            Vm0_v0(ipft) = 9.097
         case default !  Just in case. 
            Vm0_v0(ipft) = 15.625
         end select
         Vm0_v1(ipft) = 0.0
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Correct Vm0 based on SLA dependence, input settings, and the big leaf correction  !
   ! factor.  This will work even for plants whose Vm0 is not an explicit function of SLA  !
   ! (i.e. when vm0_v1 is 0), as long as SLA itself is not 0 (and it shouldn't be).        !
   !---------------------------------------------------------------------------------------!
   Vm0(:) = vm0_v0(:) * SLA(:) ** vm0_v1(:)                                                &
          * ssfact * merge(vmfact_c4,vmfact_c3,photosyn_pathway(:) == 4)
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !       Initialise Rd, Jm, and TPm terms with undefined values.  In case they are not   !
   ! provided in XML, we assign defaults, gamma-based values in subroutine                 !
   ! init_derived_params_after_xml.                                                        !
   !---------------------------------------------------------------------------------------!
   select case (iphysiol)
   case (0,2)
      !----- The default is to use the same curves for Vm, Jm (2) and Rd. -----------------!
      Jm_low_temp   (:) = undef_real
      Jm_high_temp  (:) = undef_real
      Jm_decay_elow (:) = undef_real
      Jm_decay_ehigh(:) = undef_real
      Jm_hor        (:) = undef_real
      Jm_q10        (:) = undef_real
      Rd_low_temp   (:) = undef_real
      Rd_high_temp  (:) = undef_real
      Rd_decay_elow (:) = undef_real
      Rd_decay_ehigh(:) = undef_real
      Rd_hor        (:) = undef_real
      Rd_q10        (:) = undef_real
      !------------------------------------------------------------------------------------!
   case (1,3)
      !------------------------------------------------------------------------------------!
      ! Jm: use Q10/hor from von Caemmerer (2000).  Fill other terms after xml             !
      ! initialisation.                                                                    !
      !------------------------------------------------------------------------------------!
      Jm_low_temp   (:) = undef_real
      Jm_high_temp  (:) = undef_real
      Jm_decay_elow (:) = undef_real
      Jm_decay_ehigh(:) = undef_real
      Jm_hor        (:) = 37000. * tphysref / (rmol * (t00+25.))
      Jm_q10        (:) = 1.65
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! Rd: use Q10/hor from von Caemmerer (2000).  Fill other terms after xml             !
      ! initialisation.                                                                    !
      !------------------------------------------------------------------------------------!
      Rd_low_temp   (:) = undef_real
      Rd_high_temp  (:) = undef_real
      Rd_decay_elow (:) = undef_real
      Rd_decay_ehigh(:) = undef_real
      Rd_hor        (:) = 66400. * tphysref / (rmol * (t00+25.))
      Rd_q10        (:) = Vm_q10(:)
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!




   !----- The reference values are defined by default as fractions of Vm0. ----------------!
   Rd0           (:) = undef_real
   Jm0           (:) = undef_real
   TPm0          (:) = undef_real
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Respiration terms.  Here we must check whether this will run Foley-based or        !
   ! Collatz-based photosynthesis, because the respiration/Vm ratio is constant in the     !
   ! former but not necessarily in the latter.                                             !
   !                                                                                       !
   !    Dark_respiration_factor is the lower-case gamma in Moorcroft et al. (2001).        !
   !  In case gamma is set to zero, we use the ratio derived from Atkin et al. (2015).     !
   !                                                                                       !
   ! Atkin, O. K., et al. (2015). Global variability in leaf respiration in relation to    !
   !    climate, plant functional types and leaf traits. New Phytol., 206(2), 614-636,     !
   !    10.1111/nph.13253.                                                                 !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      !----- Set Rd parameters for this PFT. ----------------------------------------------!
      if (is_grass(ipft) .and. (photosyn_pathway(ipft) == 3)) then
         a_pft     = a_c3grss
         b_pft     = b_c3grss
         gamma_pft = gamma_c3
      elseif (is_grass(ipft)) then
         a_pft     = a_c4grss
         b_pft     = b_c4grss
         gamma_pft = gamma_c4
      elseif (is_conifer(ipft)) then
         a_pft     = a_needle
         b_pft     = b_needle
         gamma_pft = gamma_c3
      elseif (is_tropical(ipft)) then
         a_pft     = a_bltrop
         b_pft     = b_bltrop
         gamma_pft = gamma_c3
      else
         a_pft     = a_bltemp
         b_pft     = b_bltemp
         gamma_pft = gamma_c3
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Decide whether to use Atkin's definition of Rd:Vm ratio or the original one.   !
      !------------------------------------------------------------------------------------!
      if (gamma_pft == 0.) then
         !---------------------------------------------------------------------------------!
         !     The ratio depends on the method: Collatz may have different q10 factors,    !
         ! and in ED the reference temperature is 15 degC, not 25 degC.  The q10 numbers   !
         ! make sure that the comparison is carried out at 25 degC, but the reference is   !
         ! defined at 15 degC.                                                             !
         !---------------------------------------------------------------------------------!
         select case (iphysiol)
         case (0,1)
            vmexpo_pft = SLA(ipft) * Vm0(ipft) / C2B
         case (2,3)
            vmexpo_pft = SLA(ipft) * vm_q10(ipft) * Vm0(ipft) / C2B
         end select
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Define dark_respiration_factor (Rd:Vm ratio) from Atkin et al. (2015).     !
         !---------------------------------------------------------------------------------!
         dark_respiration_factor(ipft) = 10.0**a_pft * vmexpo_pft ** (b_pft-1.)
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !      gamma_pft is the dark_respiration_factor (Rd:Vm ratio).                    !
         !---------------------------------------------------------------------------------!
         dark_respiration_factor(ipft) = gamma_pft
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Define the Jm0:Vm0 and TPm0:Vm0 ratios.  The default values are based on         !
   ! von Caemmerer (2000), and they are ignored in case Jm0 and TPm0 are provided through  !
   ! XML.  Jm0:Vm0 is set to 2 when iphysiol is 0 or 2 (old ED style, it assumes           !
   ! J = Jm = 2Vm).  Also, the TPm0:Vm0 ratio is ignored when iphysiol is 0 or 2 (old ED   !
   ! style, triose phosphate utilisation limitation is ignored).                           !
   !                                                                                       !
   ! Jmax/Vcmax ratio and Tpmax/Vcmax ratios were obtained from data available in B17 and  !
   ! N17 (tropical forests only, mind you).                                                !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   ! Bahar, N. H. A. et al. Leaf-level photosynthetic capacity in lowland Amazonian and    !
   !    high-elevation Andean tropical moist forests of Peru. New Phytol.,                 !
   !    214(3):1002-1018, May 2017. doi:10.1111/nph.14079 (B17).                           !
   !                                                                                       !
   ! Norby, R. J., L. Gu, I. C. Haworth, A. M. Jensen, B. L. Turner, A. P. Walker,         !
   !    J. M. Warren, D. J. Weston, C. Xu, and K. Winter. Informing models through         !
   !    empirical relationships between foliar phosphorus, nitrogen and photosynthesis     !
   !    across diverse woody species in tropical forests of Panama.  New Phytol.,          !
   !    215 (4):1425-1437, Sep 2017. doi:10.1111/nph.14319 (N17).                          !
   !---------------------------------------------------------------------------------------!
   select case (iphysiol)
   case (0,2)
      electron_transport_factor(:) = 2.0
      triose_phosphate_factor  (:) = 0.0
   case (1,3)
      electron_transport_factor(:) = 1.767
      triose_phosphate_factor  (:) = merge( 0.0, 0.110,photosyn_pathway(:) == 4)
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    The slope factor of the stomatal conductance (the "m" term in L95).                !
   !---------------------------------------------------------------------------------------!
   stomatal_slope(:) = merge( mphoto_trc3, mphoto_tec3                                     &
                            , (is_tropical(:) .or. is_grass(:)) .and. (.not. is_conifer(:)))
   stomatal_slope(:) = merge( mphoto_c4,stomatal_slope(:),photosyn_pathway(:) == 4)
   !---------------------------------------------------------------------------------------!



   !----- Define the residual stomatal conductance (aka the b term, given in umol/m2/s). --!
   cuticular_cond(:) = merge( bphoto_c4                                                    &
                            , merge(bphoto_nlc3,bphoto_blc3,is_conifer(:))                 &
                            , photosyn_pathway(:) == 4 )
   !---------------------------------------------------------------------------------------!


   !------ Set quantum yield fraction. ----------------------------------------------------!
   quantum_efficiency(:) = merge(alpha_c4,alpha_c3,photosyn_pathway(:) == 4)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters to determine the rate of electron transport.  qyield_psII is the       !
   ! quantum yield of the photosystem II, and curvpar_electron is the curvature parameter  !
   ! for the quadratic equation.  Both parameters are used only when iphysiol is 1 or 3.   !
   !---------------------------------------------------------------------------------------!
   qyield_psII     (:) = merge(0.625,0.850,photosyn_pathway(:) == 4)
   curvpar_electron(:) = 0.70
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The KW parameter. Medvigy et al. (2009) and Moorcroft et al. (2001), and the      !
   ! namelist, give the number in m2/yr/kg_C_root.  Here we must convert it to             !
   !  m2/s/kg_C_root.                                                                      !
   !---------------------------------------------------------------------------------------!
   water_conductance(:) = merge(kw_grass,kw_tree,is_grass(:)) / yr_sec
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Parameter for the trait-plasticity model for Vm0 and Rd0 (for SLA, check          !
   ! init_pft_alloc_params).  This controls the extinction factor for Vm0.  The default    !
   ! depends upon Vcmax25 (Lloyd et al. 2010, Biogeosciences) and is initialised in        !
   ! init_derived_params_after_xml.  However, this relationship is empirical, so we also   !
   ! allow it to be initialised through XML.                                               !
   !---------------------------------------------------------------------------------------!
   kplastic_vm0(:) = undef_real
   kplastic_rd0(:) = undef_real
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_pft_photo_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_resp_params()
   use ed_max_dims    , only : n_pft                     & ! intent(in)
                             , undef_real                ! ! intent(in)
   use physiology_coms, only : iphysiol                  & ! intent(in)
                             , rrffact                   & ! intent(in)
                             , growthresp                ! ! intent(in)
   use pft_coms       , only : is_tropical               & ! intent(in)
                             , is_grass                  & ! intent(in)
                             , is_conifer                & ! intent(in)
                             , leaf_turnover_rate        & ! intent(in)
                             , growth_resp_factor        & ! intent(out)
                             , root_turnover_rate        & ! intent(out)
                             , bark_turnover_rate        & ! intent(out)
                             , storage_turnover_rate     & ! intent(out)
                             , root_respiration_factor   & ! intent(out)
                             , rrf_low_temp              & ! intent(out)
                             , rrf_high_temp             & ! intent(out)
                             , rrf_decay_elow            & ! intent(out)
                             , rrf_decay_ehigh           & ! intent(out)
                             , rrf_hor                   & ! intent(out)
                             , rrf_q10                   & ! intent(out)
                             , stem_respiration_factor   & ! intent(out)
                             , stem_resp_size_scaler     & ! intent(out)
                             , srf_low_temp              & ! intent(out)
                             , srf_high_temp             & ! intent(out)
                             , srf_decay_elow            & ! intent(out)
                             , srf_decay_ehigh           & ! intent(out)
                             , srf_hor                   & ! intent(out)
                             , srf_q10                   & ! intent(out)
                             , f_labile_leaf             & ! intent(out)
                             , f_labile_stem             ! ! intent(out)
   use ed_misc_coms   , only : iallom                    & ! intent(out)
                             , economics_scheme          ! ! intent(in)
   use consts_coms    , only : onesixth                  & ! intent(in)
                             , onethird                  & ! intent(in)
                             , t00                       ! ! intent(in)
   implicit none
   !------ Local variables. ---------------------------------------------------------------!
   integer                 :: ipft
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     GPP:Growth respiration ratio.                                                     !
   !---------------------------------------------------------------------------------------!
   growth_resp_factor(:) = merge(0.4503    ,merge(onethird,0.0,is_grass(:)),is_conifer (:))
   growth_resp_factor(:) = merge(growthresp,          growth_resp_factor(:),is_tropical(:))
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Root turnover rate. Currently values are the same as leaf turnover rate, except   !
   ! for temperate trees whose values come from M09 optimization.  For tropical trees, the !
   ! default (ED-1/ED-2.0 style) is to assume that leaf and fine root turnover rate are    !
   ! the same.  In case ECONOMICS_SCHEME=1, the ratio is defined based on very limited,    !
   ! plot-level aggregated data from M11's review (assuming that NPP = maintenance costs). !
   !                                                                                       !
   ! Malhi, Y., C. Doughty, and D. Galbraith. The allocation of ecosystem net primary      !
   !     productivity in tropical forests. Philos. Trans. R. Soc. B-Biol. Sci.,            !
   !     366(1582):3225-3245, Nov 2011. doi:10.1098/rstb.2011.0062.                        !
   !                                                                                       !
   ! Medvigy, D. M., S. C. Wofsy, J. W. Munger, D. Y. Hollinger, and P. R. Moorcroft.      !
   !    Mechanistic scaling of ecosystem function and dynamics in space and time:          !
   !    Ecosystem demography model version 2. J. Geophys. Res.-Biogeosci., 114(G1):G01002, !
   !    Jan 2009. doi:10.1029/2008JG000812.                                                !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      !------------------------------------------------------------------------------------!
      !      Temperate trees have optimised turnover rates.  Fill each PFT separately.     !
      !------------------------------------------------------------------------------------!
      select case (ipft)
      case (6)     ! Northern pines. 
         root_turnover_rate(ipft) = 3.927218 ! 0.333
      case (7)     ! Southern pines. 
         root_turnover_rate(ipft) = 4.117847 ! 0.333
      case (8)     ! Late conifers. 
         root_turnover_rate(ipft) = 3.800132 ! 0.333
      case (9)     ! Early hardwoods. 
         root_turnover_rate(ipft) = 5.772506
      case (10)    ! Mid hardwoods. 
         root_turnover_rate(ipft) = 5.083700
      case (11)    ! Late hardwoods. 
         root_turnover_rate(ipft) = 5.070992
      case default ! Grasses, tropical trees, and lianas
         select case (economics_scheme)
         case (1)
            root_turnover_rate(ipft) = 0.9 * leaf_turnover_rate(ipft)
         case default
            root_turnover_rate(ipft) = leaf_turnover_rate(ipft)
         end select
      end select
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Bark turnover rate.   A first-guess estimate is the ratio between bark increment   !
   ! and bark thickness.  Increment rates came from CD17.  Derived values were of the      !
   ! order of 0.40/yr, and for the time being we use this number for all PFTs.  Actual     !
   ! data on bark turnover could be helpful.                                               !
   !                                                                                       !
   ! Charles-Dominique, T., G. F. Midgley, and W. J. Bond. Fire frequency filters species  !
   !    by bark traits in a savanna-forest mosaic. J. Veg. Sci., 28(4):728-735, Jul 2017.  !
   !    doi:10.1111/jvs.12528.                                                             !
   !---------------------------------------------------------------------------------------!
   bark_turnover_rate(:) = 0.40
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Storage turnover rate.  Number for tropical trees is based on tuning for the time  !
   ! being.                                                                                !
   !---------------------------------------------------------------------------------------!
    storage_turnover_rate(:) = merge( merge(onethird,onesixth,is_grass(:))                 &
                                    , merge(0.0,0.6243,is_grass(:) .or. is_conifer(:))     &
                                    , is_tropical(:) )
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Labile fraction of leaf-like tissues (leaves and fine roots) and wood-like tissues !
   ! (sapwood, heartwood, bark).  These names are rather confusing but I adopted the same  !
   ! names already in use in c2n factors.                                                  !
   !---------------------------------------------------------------------------------------!
   select case (iallom)
   case (2,3,4,5)
      !------------------------------------------------------------------------------------!
      !   For tropical leaves/fine roots, assume the metabolic/structural ratio obtained   !
      ! by B17.  For grasses and temperate plants, use B17 equation and R96 values for     !
      !lignin and C:N ratio. For now, assume sapwood, bark, and heartwood are 100%         !
      ! structural (even though some studies use mixed structural and metabolic, see S09). !
      !                                                                                    !
      ! References:                                                                        !
      !                                                                                    !
      ! Brechet L, Le Dantec V, Ponton S , Goret JY, Sayer E, Bonal D, Freycon V, Roy J,   !
      !    Epron D. 2017. Short- and long-term influence of litter quality and quantity on !
      !    simulated heterotrophic soil respiration in a lowland tropical forest.          !
      !    Ecosystems, 20(6):1190-1204. doi:10.1007/s10021-016-0104-x (B17).               !
      !                                                                                    !
      ! Bolker BM, Pacala SW, and Parton WJ. 1998. Linear analysis of soil decomposition:  !
      !    insights from the CENTURY model. Ecol. Appl., 8(2):425-439.                     !
      !    doi:10.1890/1051- 0761(1998)008[0425:LAOSDI]2.0.CO;2 (B98).                     !
      !                                                                                    !
      ! Randerson JT, Thompson MV, Malmstrom CM, Field CB, and Fung IY. 1996. Substrate    !
      !    limitations for heterotrophs: Implications for models that estimate the         !
      !    seasonal cycle of atmospheric CO2. Global Biogeochem. Cycles, 10(4):585-602.    !
      !    doi:10.1029/96GB01981 (R96).                                                    !
      !                                                                                    !
      ! Shevliakova E, Pacala SW, Malyshev S, Hurtt GC, Milly PCD, Caspersen JP,           !
      !    Sentman LT, Fisk JP, Wirth C, and Crevoisier C. 2009.  Carbon cycling under 300 !
      !    years of land use change: Importance of the secondary vegetation sink. Global   !
      !    Biogeochem. Cycles, 23(2):GB2022. doi:10.1029/2007GB003176 (S09).               !
      !------------------------------------------------------------------------------------!
      f_labile_leaf(:) = merge( 0.78                                                       &
                              , merge( 0.32                                                &
                                     , merge(0.75,0.50,is_conifer(:))                      &
                                     , is_tropical(:)                 )                    &
                              , is_grass(:)                             )
      f_labile_stem(:) = 0.0
      !------------------------------------------------------------------------------------!
   case default
      f_labile_leaf(:) = merge( 0.79                                                       &
                              , merge(1.00,0.79,is_tropical(:) .or. is_grass(:))           &
                              , is_conifer(:)                                    )
      f_labile_stem(:) = 0.0
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    This variable sets the contribution of roots to respiration at the reference       !
   ! temperature of 15C.  Its units is umol_CO2/kg_fine_roots/s.                           !
   !---------------------------------------------------------------------------------------!
   select case (iphysiol)
   case (0,1)
      !----- Arrhenius function. ----------------------------------------------------------!
      root_respiration_factor(:) = 0.528 * rrffact
      !------------------------------------------------------------------------------------!
   case (2,3)
      !----- Q10 function. ----------------------------------------------------------------!
      root_respiration_factor(:) = 0.2455212 * rrffact
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Initialise root respiration terms with undefined, and attribute default values    !
   ! only in case they are not defined through XML (check init_derived_params_after_xml).  !
   !---------------------------------------------------------------------------------------!
   rrf_low_temp   (:) = undef_real
   rrf_high_temp  (:) = undef_real
   rrf_decay_elow (:) = undef_real
   rrf_decay_ehigh(:) = undef_real
   rrf_hor        (:) = undef_real
   rrf_q10        (:) = undef_real
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Initialise stem respiration terms with undefined, and attribute default values    !
   ! only in case they are not defined through XML (check init_derived_params_after_xml).  !
   !---------------------------------------------------------------------------------------!
   stem_respiration_factor(:) = 10. ** (-0.672 - 0.19) / 2.2 ! umol_CO2/m2_stem_area/s from Chambers et al. 2004
   stem_resp_size_scaler(:)   = 0.0041  ! cm_DBH -1, value based on troipcal trees
   srf_low_temp   (:) = undef_real
   srf_high_temp  (:) = undef_real
   srf_decay_elow (:) = undef_real
   srf_decay_ehigh(:) = undef_real
   srf_hor        (:) = undef_real
   srf_q10        (:) = undef_real
   !---------------------------------------------------------------------------------------!


   return
end subroutine init_pft_resp_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine assigns some PFT-dependent parameters that control mortality rates.  !
!------------------------------------------------------------------------------------------!
subroutine init_pft_mort_params()
   use ed_max_dims ,    only : n_pft                      ! ! intent(in)
   use pft_coms    ,    only : is_grass                   & ! intent(in)
                             , is_tropical                & ! intent(in)
                             , is_conifer                 & ! intent(in)
                             , is_liana                   & ! intent(in)
                             , mort0                      & ! intent(out)
                             , mort1                      & ! intent(out)
                             , mort2                      & ! intent(out)
                             , mort3                      & ! intent(out)
                             , hydro_mort0                & ! intent(out)
                             , hydro_mort1                & ! intent(out)
                             , cbr_severe_stress          & ! intent(out)
                             , rho                        & ! intent(out)
                             , seedling_mortality         & ! intent(out)
                             , treefall_s_gtht            & ! intent(out)
                             , treefall_s_ltht            & ! intent(out)
                             , fire_s_min                 & ! intent(out)
                             , fire_s_max                 & ! intent(out)
                             , fire_s_inter               & ! intent(out)
                             , fire_s_slope               & ! intent(out)
                             , plant_min_temp             & ! intent(out)
                             , frost_mort                 ! ! intent(out)
   use consts_coms ,    only : t00                        & ! intent(in)
                             , lnexp_max                  & ! intent(in)
                             , onethird                   & ! intent(in)
                             , twothirds                  ! ! intent(in)
   use ed_misc_coms,    only : ibigleaf                   & ! intent(in)
                             , economics_scheme           ! ! intent(in)
   use disturb_coms,    only : include_fire               & ! intent(in)
                             , time2canopy                & ! intent(in)
                             , treefall_disturbance_rate  ! ! intent(in)
   use physiology_coms, only : carbon_mortality_scheme    & ! intent(in)
                             , hydraulic_mortality_scheme ! ! intent(in)
   implicit none

   !----- Local variables. ----------------------------------------------------------------!
   real                   :: aquad
   real                   :: bquad
   real                   :: cquad
   real                   :: discr
   real                   :: lambda_ref
   real                   :: lambda_eff
   real                   :: leff_neg
   real                   :: leff_pos
   real, dimension(n_pft) :: rho_use
   integer                :: ipft
   !----- Local constants for C18 mortality (see below). ----------------------------------!
   real, parameter        :: rho_c18      =  0.600     ! C18 ref. wood density
   real, parameter        :: alpha0_c18   =  0.0498774 ! scale for alpha
   real, parameter        :: alpha1_c18   = -0.5598224 ! exponent for alpha
   real, parameter        :: beta0_c18    = 23.6258866 ! scale for beta
   real, parameter        :: beta1_c18    =  0.3672854 ! exponent for beta
   real, parameter        :: gamma0_c18   =  0.0110928 ! scale for gamma
   real, parameter        :: gamma1_c18   = -2.2347380 ! exponent for gamma
   real, parameter        :: delta_c18    =  0.8998337 ! average correction term
   real, parameter        :: epsilon_ed2  =  1.20      ! quick scaling factor to transform 
                                                       !    CBAREL into growth
   real, parameter        :: extinct_mort = 10.0       ! Maximum mortality rate for C18
                                                       !    This should be high to cause
                                                       !    near-extinction but not so 
                                                       !    high that it would be 
                                                       !    difficult to display in output
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Frost mortality is assumed the same for all PFTs (hardiness will be given by     !
   ! plant_min_temp.                                                                       !
   !---------------------------------------------------------------------------------------!
   frost_mort(:) = 3.0
   do ipft=1,n_pft
      if (is_tropical(ipft) .and. (.not. is_conifer(ipft))) then
         !----- Tropical trees, grasses, and lianas. --------------------------------------!
         plant_min_temp(ipft) = t00 + 2.5
         !---------------------------------------------------------------------------------!
      else
         !----- Sub-tropical and temperate plants.  PFT-specific values. ------------------!
         select case (ipft)
         case (5,6,9)
            plant_min_temp(ipft) = t00 - 80.0
         case (7)
            plant_min_temp(ipft) = t00 - 10.0
         case (8)
            plant_min_temp(ipft) = t00 - 60.0
         case (10,11)
            plant_min_temp(ipft) = t00 - 20.0
         case default
            plant_min_temp(ipft) = t00 - 15.0
         end select
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     The following variables control the density-dependent mortality rates.  The       !
   ! functional form depends on the trait data base, be sure to familiarise with the       !
   ! differences before changing these parameters.  Also, keep in mind that both           !
   ! approaches have caveats.                                                              !
   !                                                                                       !
   !   CARBON_MORTALITY_SCHEME = 0                                                         !
   !  ----------------------                                                               !
   !     This follows ED-1 original formulation (M01):                                     !
   !                                                                                       !
   !                     DD = mort1 / (1 + exp(mort2 * (CBAREL - mort0))                   !
   !                                                                                       !
   !     Parameters are not based on literature, but tuning.  They allow mortality to      !
   ! quickly increase to maximum (set by mort1) when carbon balance is negative.  The      !
   ! term mort0 (typically negative) allows to offset the curve so mortality is not        !
   ! too high when carbon balance is still positive but low.                               !
   !                                                                                       !
   !   CARBON_MORTALITY_SCHEME = 1                                                         !
   !  ----------------------                                                               !
   !     This follows the trait-dependent exponential model by C18:                        !
   !                                                                                       !
   !                     DD = mort1 * exp( - mort2 * (CBAREL-mort0) )                      !
   !                                                                                       !
   !     Note that in the original C18 model, growth rates are used instead of CB.  We     !
   ! translate the growth rates into CBAREL using a simple conversion of 1.20.  This       !
   ! factor came from a quick look on one simulation for Paracou, so mind that it is       !
   ! highly speculative.  Using CBAREL makes catastrophic mortality in case carbon balance !
   ! is negative.  Hazard rates are unbounded in this case, which may need to be           !
   ! revisited.                                                                            !
   !                                                                                       !
   !   CARBON_MORTALITY_SCHEME = 2                                                         !
   !  ----------------------                                                               !
   !     This follows the trait-dependent exponential model by C18:                        !
   !                                                                                       !
   !                     DD = mort1 * exp( - mort2 * annual grwoth rate)                   !
   !                                                                                       !
   !     This option uses the raw growth-based relationship reported in C18.of CB. Note    !
   ! that the maximum mortality is bounded. However, if cohorts have no alive biomass      !
   ! they will be flagged as inviable and removed in fuse_fiss_utils.                      !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   !                                                                                       !
   !  Moorcroft, PR, Hurtt GC, and Pacala SW. A method for scaling vegetation dynamics:    !
   !     The Ecosystem Demography model (ED). Ecol. Monogr., 71(4):557-586, Nov 2001.      !
   !     doi:10.1890/0012- 9615(2001)071[0557:AMFSVD]2.0.CO;2 (M01).                       !
   !                                                                                       !
   !  Camac, JS, Condit R, FitzJohn RG, McCalman L, Steinberg D, Westoby M, Wright SJ,     !
   !     Falster DS. Partitioning mortality into growth-dependent and growth-independent   !
   !     hazards across 203 tropical tree species, Proc. Natl. Acad. Sci. U. S. A., 115    !
   !     (49), 12459-12464, Dec 2018.  doi:10.1073/pnas.1721040115 (C18).                  !
   !---------------------------------------------------------------------------------------!
   select case (carbon_mortality_scheme)
   case (2)
      rho_use(:) = merge(rho_c18                                                           &
                        ,rho(:)                                                            &
                        ,is_grass(:) .or. is_liana(:) .or. (.not. is_tropical(:)) )
      mort0(:) = 0.0
      mort1(:) = delta_c18   * alpha0_c18 * (rho_use(:)/rho_c18) ** alpha1_c18
      mort2(:) = beta0_c18  * (rho_use(:)/rho_c18) ** beta1_c18
   case (1)
      rho_use(:) = merge(rho_c18                                                           &
                        ,rho(:)                                                            &
                        ,is_grass(:) .or. is_liana(:) .or. (.not. is_tropical(:)) )
      mort0(:) = 0.0
      mort1(:) = delta_c18   * alpha0_c18 * (rho_use(:)/rho_c18) ** alpha1_c18
      mort2(:) = epsilon_ed2 * beta0_c18  * (rho_use(:)/rho_c18) ** beta1_c18
   case default
      mort0(:) = merge(-0.30,  0.0,is_tropical(:))
      mort1(:) = merge(  1.0,  1.0,is_tropical(:))
      mort2(:) = merge( 20.0, 20.0,is_tropical(:))
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Variable mort3 controls the density-independent mortality rate due to ageing.     !
   ! This value is a constant in "hazard rates" (i.e. m in dn/dt = - m * n).               !
   !                                                                                       !
   !     For tropical trees, parameter choice will also depend on economics_scheme.        !
   !                                                                                       !
   ! With lianas the idea is to give it the normal mort3 mortality that is basically rho   !
   ! dependent. Then in the next phase I will try to add a limit to the maximum liana load !
   ! that a tree can born. When the number of lianas hosted by a given tree exceeds the    !
   ! avreage values given in O. Phillips (2005) I will kill the lianas(...and the trees?..)!
   ! That will increase mortality as a consequence. I should later check if the resulting  !
   ! mortality is in line with what stated in the aforementioned article (lianas have +-6% !
   ! turnover rate, 3 times that of trees...).                                             !
   ! MLO to MdPB: Liana load mortality should _NOT_ be part of mort3, because you are      !
   !              describing a density-dependent mortality.  We may need a mort4 parameter !
   !              for that.                                                                !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      if (is_liana(ipft)) then ! Lianas, taken from O. Phillips (2005). 
         mort3(ipft) = 0.06311576
      elseif (is_tropical(ipft) .and. is_conifer(ipft)) then
         mort3(ipft) = 0.00111 ! Based on TRY (but likely guesstimated).
      elseif (is_tropical(ipft)) then
          select case (carbon_mortality_scheme)
          case (1,2)
             !----- Updated parameters. ---------------------------------------------------!
             if (is_grass(ipft)) then
                !--------------------------------------------------------------------------!
                !    Median plant lifespan of the grasses and herbs included in the TRY    !
                ! database (K11, see reference below).                                     !
                !--------------------------------------------------------------------------!
                mort3(ipft) = 0.124
                !--------------------------------------------------------------------------!
             else
                !--------------------------------------------------------------------------!
                !      The slope was defined based on the 50-ha inventory data from Barro  !
                ! Colorado (C12), using the growth-independent mortality function from     !
                ! C18, which they call "baseline mortality".                               !
                !                                                                          !
                ! Note: we may be double counting treefall mortality rates, but mort3      !
                ! would become negative for mid- and late-successional PFTs when using a   !
                ! typical 0.014 yr-1 (which may be too large anyway).                      !
                !--------------------------------------------------------------------------!
                mort3(ipft) = delta_c18   * gamma0_c18 * (rho(ipft)/rho_c18) ** gamma1_c18
                !--------------------------------------------------------------------------!
             end if
             !-----------------------------------------------------------------------------!
          case default
             !----- Default parameters (ED-1.0). ------------------------------------------!
             if (is_grass(ipft)) then ! Tropical grasses.  Same as temperate grasses
                mort3(ipft) = 0.066
             else                     ! Tropical trees, use ED-1 approach.
                mort3(ipft) = 0.15 * (1. - rho(ipft) / 0.90)
             end if
             !-----------------------------------------------------------------------------!
          end select
      else
         !----- Temperate trees, use optimised values from Medvigy et al. (2009). ---------!
         select case (ipft)
         case (5)     ! Temperate C3 grasses.
            mort3(ipft) = 0.066
         case (6)     ! Northern pines. 
            mort3(ipft) = 0.0033928
         case (7)     ! Southern pines. 
            mort3(ipft) = 0.0043
         case (8)     ! Late-successional conifers.
            mort3(ipft) = 0.0023568
         case (9)     ! Early-successional hardwoods. 
            mort3(ipft) = 0.006144
         case (10)    ! Mid-successional hardwoods.
            mort3(ipft) = 0.003808
         case (11)    ! Late-successional hardwoods.
            mort3(ipft) = 0.00428 ! Late-successional greater than than mid-successional?
         case default ! Forgotten PFT
            mort3(ipft) = 0.005
         end select
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      This is the default mortality for when the maximum carbon balance is negative.   !
   !  Default in ED-1.0 and ED-2.0 was to assume zero, an alternative is to assume maximum !
   !  mortality.                                                                           !
   !---------------------------------------------------------------------------------------!
   select case (carbon_mortality_scheme)
   case (2)
      cbr_severe_stress(:) = 0.0 ! not used in this case
   case (1)
      cbr_severe_stress(:) = mort0(:) - 1.0 / mort2(:) * log(extinct_mort/mort1(:))
   case default
      cbr_severe_stress(:) = log(epsilon(1.0)) / mort2(:)
   end select
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     The following variables control the hydraulic failure mortality rates. It is      !
   ! highly experimental and assumes a log-linear relationship between mortality rates and !
   ! percentage loss of conductance (PLC).                                                 !
   ! Note that hydraulic failure mortality for grasses are set to be always zero because   !
   ! they are treated with no stems in the model                                           !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   select case (hydraulic_mortality_scheme)
   case (1)
       ! Assume cohort will die in half month when PLC is 1. (100%, all conductance is lost)
       ! and die in 1. year when PLC is 0.6 (Estimated from Adams et al. 2017 and 
       ! Hammond et al. 2019)
       hydro_mort0(:) = merge(0.0,24.0,is_grass(:))
       do ipft=1,n_pft
           if (is_grass(ipft)) then
               hydro_mort1(ipft) = 1.0
           else
               hydro_mort1(ipft) = - log(hydro_mort0(ipft)) / log(0.6)
           endif
       enddo
   case default
       ! ED 2.2 default, no hydraulic failure mortality
       hydro_mort0(:) = 0.0
       hydro_mort1(:) = 1.0
   end select
       
 
   !---------------------------------------------------------------------------------------!
   !     Here we check whether we need to re-calculate the treefall disturbance rate so it !
   ! is consistent with the time to reach the canopy.                                      !
   !---------------------------------------------------------------------------------------!
   if (treefall_disturbance_rate == 0.) then
      !------ No disturbance rate, set time to reach canopy to infinity. ------------------!
      time2canopy = huge(1.)
      lambda_ref  = 0.
      lambda_eff  = 0.
      !------------------------------------------------------------------------------------!

   else
      lambda_ref = treefall_disturbance_rate

      if (time2canopy > 0.) then
         !---------------------------------------------------------------------------------!
         !     We are not going to knock down trees as soon as the patch is created;       !
         ! instead, we will wait until the patch age is older than time2canopy.  We want,  !
         ! however, to make the mean patch age to be 1/treefall_disturbance_rate.  The     !
         ! equation below can be retrieved by integrating the steady-state probability     !
         ! distribution function.  The equation is quadratic and the discriminant will     !
         ! never be zero and the treefall_disturbance_rate will be always positive because !
         ! the values of time2canopy and treefall_disturbance_rate have already been       !
         ! tested in ed_opspec.F90.                                                        !
         !---------------------------------------------------------------------------------!
         aquad    = time2canopy * time2canopy * lambda_ref  - 2. * time2canopy
         bquad    = 2. * time2canopy * lambda_ref - 2.
         cquad    = 2. * lambda_ref
         !------ Find the discriminant. ---------------------------------------------------!
         discr    = bquad * bquad - 4. * aquad * cquad
         leff_neg = - 0.5 * (bquad - sqrt(discr)) / aquad
         leff_pos = - 0.5 * (bquad + sqrt(discr)) / aquad
         !---------------------------------------------------------------------------------!
         !      Use the maximum value, but don't let the value to be too large otherwise   !
         ! the negative exponential will cause underflow.                                  !
         !---------------------------------------------------------------------------------!
         lambda_eff = min(lnexp_max,max(leff_neg,leff_pos))
         !---------------------------------------------------------------------------------!
      else
         lambda_eff = lambda_ref
      end if

   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Print out the summary.                                                             !
   !---------------------------------------------------------------------------------------!
   write (unit=*,fmt='(a)')           '----------------------------------------'
   write (unit=*,fmt='(a)')           '  Treefall disturbance parameters:'
   write (unit=*,fmt='(a,1x,es12.5)') '  - LAMBDA_REF  =',lambda_ref
   write (unit=*,fmt='(a,1x,es12.5)') '  - LAMBDA_EFF  =',lambda_eff
   write (unit=*,fmt='(a,1x,es12.5)') '  - TIME2CANOPY =',time2canopy
   write (unit=*,fmt='(a)')           '----------------------------------------'
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !     Seedling mortality must be redefined for big leaf runs: this is necessary because !
   ! big leaf plants don't grow in diameter, but in "population".                          !
   !---------------------------------------------------------------------------------------!
   select case (ibigleaf)
   case (0)
      seedling_mortality(:) = 0.95
   case (1)
      seedling_mortality(:) = merge(0.9500,onethird,is_grass(:))
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Treefall survivorship fraction.                                                  !
   !---------------------------------------------------------------------------------------!
   !----- Trees taller than treefall_hite_threshold (liana survivorship: Putz 1983). ------!
   treefall_s_gtht(:) = merge(0.80,0.00,is_liana(:))
   !----- Trees shorter than treefall_hite_threshold. -------------------------------------!
   select case (economics_scheme)
   case (1)
      !------------------------------------------------------------------------------------!
      !    Higher survivorship.  Tree survivorship is based on measurements at Paracou.    !
      ! Estimate is the average ratio between secondary and primary tree fall mortality.   !
      !                                                                                    !
      ! Ferry, B. et al.  Higher treefall rates on slopes and waterlogged soils result     !
      !    in lower stand biomass and productivity in a tropical rain forest. J. Veg.      !
      !    Sci., 98(1), 106-116, 2010. doi:10.1111/j.1365-2745.2009.01604.x                !
      !                                                                                    !
      !    Grasses are assigned high survivorship as they are unlikely to be crushed to    !
      ! death.                                                                             !
      !------------------------------------------------------------------------------------!
      treefall_s_ltht(:) = merge(0.80,0.30,is_grass(:) .or. is_liana(:))
      !------------------------------------------------------------------------------------!
   case default
      !------------------------------------------------------------------------------------!
      !    Original parameters.  MLO: Is it reasonable that small lianas have lower        !
      ! survivorship than tall lianas?.                                                    !
      !------------------------------------------------------------------------------------!
      treefall_s_ltht(:) = merge(0.25,0.10,is_grass(:))
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Fire survivorship fraction.  These variables will be replaced by bark thickness  !                                                    !
   !---------------------------------------------------------------------------------------!
   select case (include_fire)
   case (3)
      fire_s_min  (:) = merge(0.2,0.2,is_grass(:))
      fire_s_max  (:) = merge(0.2,0.9,is_grass(:))
      fire_s_inter(:) =  1.5
      fire_s_slope(:) = -1.0
   case default
      fire_s_min  (:) =  0.0
      fire_s_max  (:) =  0.0
      fire_s_inter(:) =  1.5
      fire_s_slope(:) = -1.0
   end select
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_pft_mort_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_pft_nitro_params()
   use pft_coms    , only : Vm0                  & ! intent(in)
                          , SLA                  & ! intent(in)
                          , is_conifer           & ! intent(in)
                          , is_tropical          & ! intent(in)
                          , is_liana             & ! intent(in)
                          , c2n_leaf             & ! intent(out)
                          , c2n_storage          & ! intent(out)
                          , c2n_stem             & ! intent(out)
                          , l2n_stem             & ! intent(out)
                          , plant_N_supply_scale ! ! intent(out)
   use decomp_coms , only : c2n_slow             & ! intent(out)
                          , c2n_fast_0           & ! intent(out)
                          , c2n_structural       ! ! intent(out)
   use ed_max_dims , only : n_pft                ! ! intent(in)
   use ed_misc_coms, only : economics_scheme     ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer :: ipft
   real    :: vm0_ref
   !---------------------------------------------------------------------------------------!
   pftloop: do ipft=1,n_pft
      !------ Temperate trees. Use non-optimised Vm0 values. ------------------------------!
      ! MLO - There must be something uninitialised.  This works when I use if, but it 
      !       fails with -O3 when I use select case (which normally should be safer).
      if (ipft == 6 .or. ipft == 7) then ! Northern and Southern Pines. 
         vm0_ref = 15.625
      elseif (ipft == 8) then     ! Late conifers. 
         vm0_ref = 6.25
      elseif (ipft == 9) then     ! Early hardwood. 
         vm0_ref = 18.25
      elseif (ipft == 10) then   ! Mid hardwood. 
         vm0_ref = 15.625
      elseif (ipft == 11) then   ! Late hardwood. 
         vm0_ref = 6.25
      else ! Use actual Vm0 in case we missed some PFT. 
         vm0_ref = Vm0(ipft)
      end if
      !------------------------------------------------------------------------------------!


      !----- Ratio. -----------------------------------------------------------------------!
      select case (economics_scheme)
      case (1)
         if (is_liana(ipft) .or. (.not. is_tropical(ipft))) then
            !----- Use ED-1 default. ------------------------------------------------------!
            c2n_leaf(ipft) = 1000.0 / ( (0.11289 + 0.12947 *   vm0_ref) * SLA(ipft) )
            !------------------------------------------------------------------------------!
         elseif (is_conifer(ipft)) then
            !----- Araucaria, use values from TRY. ----------------------------------------!
            c2n_leaf(ipft) = 86.29189
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !    SMA fit obtained using the available data from W04, B17, and N17, with    !
            ! SLA complemented with data from the TRY data base (K11).  The model fit for  !
            ! grasses/herbs/shrubs were almost the same as for flowering trees, therefore  !
            ! a single fit for all was used.                                               !
            !                                                                              !
            ! References:                                                                  !
            !                                                                              !
            ! Bahar, N. H. A. et al. Leaf-level photosynthetic capacity in lowland         !
            !    Amazonian and high-elevation Andean tropical moist forests of Peru. New   !
            !    Phytol., 214(3):1002-1018, May 2017. doi:10.1111/nph.14079 (B17).         !
            !                                                                              !
            ! Kattge, J., S. Diaz, S. Lavorel, et al., TRY -- a global database of plant   !
            !    traits. Glob. Change Biol., 17 (9): 2905-2935, Sep 2011.                  !
            !    doi:10.1111/j.1365-2486.2011.02451.x (K11).                               !
            !                                                                              !
            ! Norby, R. J., L. Gu, I. C. Haworth, A. M. Jensen, B. L. Turner,              !
            !    A. P. Walker, J. M. Warren, D. J. Weston, C. Xu, and K. Winter. Informing !
            !    models through empirical relationships between foliar phosphorus,         !
            !    nitrogen and photosynthesis across diverse woody species in tropical      !
            !    forests of Panama.  New Phytol., 215 (4):1425-1437, Sep 2017.             !
            !    doi:10.1111/nph.14319 (N17).                                              !
            !                                                                              !
            ! Wright, I. J., P. B. Reich, M. Westoby, et al., The worldwide leaf           !
            !    economics spectrum. Nature, 428(6985):821-827, Apr 2004.                  !
            !    doi:10.1038/nature02403 (W04).                                            !
            !------------------------------------------------------------------------------!
            c2n_leaf(ipft) = 337.959 / SLA(ipft) ** 0.834527
            !------------------------------------------------------------------------------!
         end if
      case default
         c2n_leaf(ipft) = 1000.0 / ( (0.11289 + 0.12947 *   vm0_ref) * SLA(ipft) )
      end select
      !------------------------------------------------------------------------------------!
   end do pftloop
   !---------------------------------------------------------------------------------------!

   c2n_stem(:)          = 150.0 ! Carbon to Nitrogen ratio, structural stem.
   c2n_storage          = 150.0 ! Carbon to Nitrogen ratio, storage pool.
   l2n_stem             = 150.0 ! Carbon to Lignin ratio, structural stem.
   plant_N_supply_scale = 0.5 


   c2n_fast_0           = 30.0  ! Carbon to Nitrogen ratio, slow pool.
   c2n_slow             = 10.0  ! Carbon to Nitrogen ratio, slow pool.
   c2n_structural       = 150.0 ! Carbon to Nitrogen ratio, structural pool.

   return
end subroutine init_pft_nitro_params
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
!  SUBROUTINE: INIT_PFT_HYDRO_PARAMS   
!> \brief This subroutine initializes plant hydraulic parameters.
!> \warning Some plant hydraulic paramters are depedent on SLA and LMA. So this
!> subroutine should always be after init_pft_alloc_params
!> \author Xiangtao Xu, Jan. 27 2018
!==========================================================================================!
subroutine init_pft_hydro_params()

   use consts_coms    , only : grav                 & ! intent(in)
                             , wdns                 ! ! intent(in)
   use ed_max_dims    , only : n_pft                ! ! intent(in)
   use physiology_coms, only : plant_hydro_scheme   ! ! intent(in)
   use plant_hydro    , only : rwc2psi              & ! subroutine
                             , psi2rwc              ! ! subroutine
   use ed_misc_coms   , only : economics_scheme     ! ! intent(in)
   use pft_coms       , only : is_grass             & ! intent(in)
                             , is_conifer           & ! intent(in)
                             , is_liana             & ! intent(in)
                             , is_tropical          & ! intent(in)
                             , photosyn_pathway     & ! intent(in)
                             , SLA                  & ! intent(in)
                             , rho                  & ! intent(in)
                             , C2B                  & ! intent(in)
                             , vessel_curl_factor   & ! intent(out)
                             , leaf_water_cap       & ! intent(out)
                             , wood_water_cap       & ! intent(out)
                             , leaf_water_sat       & ! intent(out)
                             , wood_water_sat       & ! intent(out)
                             , bark_water_sat       & ! intent(out)
                             , leaf_rwc_min         & ! intent(out)
                             , wood_rwc_min         & ! intent(out)
                             , leaf_psi_min         & ! intent(out)
                             , wood_psi_min         & ! intent(out)
                             , leaf_psi_osmotic     & ! intent(out)
                             , wood_psi_osmotic     & ! intent(out)
                             , leaf_elastic_mod     & ! intent(out)
                             , wood_elastic_mod     & ! intent(out)
                             , leaf_psi_tlp         & ! intent(out)
                             , wood_psi_tlp         & ! intent(out)
                             , wood_Kmax            & ! intent(out)
                             , wood_Kexp            & ! intent(out)
                             , wood_psi50           & ! intent(out)
                             , stoma_lambda         & ! intent(out)
                             , stoma_beta           & ! intent(out)
                             , stoma_psi_b          & ! intent(out)
                             , stoma_psi_c          ! ! intent(out)

   implicit none
   !------ Local variables. ---------------------------------------------------------------!
   integer                   :: ipft
   real   , dimension(n_pft) :: rwc_tlp_wood ! RWC at turgor loss point for sapwood
   !real   , dimension(n_pft) :: leaf_density ! density of leaf tissue [kg/m3]
   real   , dimension(n_pft) :: LMA          ! leaf mass per area     [ g/m2]
   !real   , dimension(n_pft) :: Amax_25      ! estimated max. photosynthetic rates at 25C
   real   , dimension(n_pft) :: rho_bnd      ! Bounded wood density, to avoid FPE.
   logical, dimension(n_pft) :: is_troptree  ! Flag to select only tropical trees.
   !------ Local parameters. --------------------------------------------------------------!
   real   , parameter        :: MPa2m   = wdns / grav
   real   , parameter        :: rho_min = 0.35 ! Minimum wood density
   real   , parameter        :: rho_max = 0.95 ! Maximum wood density
   real   , parameter        :: f_cap   = 0.07 ! Fraction of capillary water.
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   ! References for this sub-routine.:                                                     !
   !                                                                                       !
   ! Bartlett MK, Scoffoni C , Sack L. 2012. The determinants of leaf turgor loss point    !
   !    and predic- tion of drought tolerance of species and biomes: a global meta-        !
   !    analysis. Ecol. Lett., 15: 393-405. doi:10.1111/j.1461-0248.2012.01751.x (B12).    !
   !                                                                                       !
   ! Christoffersen BO, Gloor M, Fauset S, Fyllas NM, Galbraith DR, Baker TR, Kruijt B,    !
   !    Rowland L, Fisher RA, Binks OJ et al. 2016. Linking hydraulic traits to tropical   !
   !    forest function in a size- structured and trait-driven model (TFS v.1-Hydro).      !
   !    Geosci. Model Dev., 9: 4227-4255. doi:10.5194/gmd-9- 4227-2016 (C16).              !
   !                                                                                       !
   ! Forest Products Laboratory. Wood handbook - wood as an engineering material. General  !
   !    Technical Report FPL-GTR-190, U.S. Department of Agriculture, Madison, WI, 2010.   !
   !    doi:10.2737/FPL-GTR-190 (FPL10)                                                    !
   !                                                                                       !
   ! Gu, L, T. Meyers, S. G. Pallardy, P. J. Hanson, B. Yang, M. Heuer, K. P. Hosman,     !
   !    Q. Liu, J. S. Riggs, D. Sluss, and S. D. Wullschleger. Influences of biomass heat  !
   !    and biochemical energy storages on the land surface fluxes and radiative temper-   !
   !    ature. J. Geophys. Res., 112(D2):D02107, Jan 2007. doi:10.1029/2006JD007425 (G07)  !
   !                                                                                       !
   ! Kursar TA, Engelbrecht BMJ, Burke A, Tyree MT, EI Omari B , Giraldo JP. 2009.         !
   !    Tolerance to low leaf water status of tropical tree seedlings is related to        !
   !    drought performance and distribution. Funct. Ecol., 23: 93-102.                    !
   !    doi:10.1111/j.1365-2435.2008.01483.x (K09).                                        !
   !                                                                                       !
   ! Lin, Y.-S., B. E. Medlyn, R. A. Duursma, I. C. Prentice, H. Wang, S. Baig, D. Eamus,  !
   ! V. R. deDios, P. Mitchell, D. S. Ellsworth, M. O. de Beeck, G. Wallin, J. Uddling,    !
   ! L. Tarvainen, M.-L. Linderson, L. A. Cernusak, J. B. Nippert, T. W. Ocheltree,        !
   ! D. T. Tissue, N. K. Martin-StPaul, A. Rogers, J. M. Warren, P. De Angelis, K. Hikosaka!
   !, Q. Han, Y. Onoda, T. E. Gimeno, C. V. M. Barton, J. Bennie, D. Bonal, A. Bosc, M. Low!
   !, C. Macinins-Ng, A. Rey, L. Rowland, S. A. Setterfield, S. Tausz-Posch,               !
   ! J. Zaragoza-Castells, M. S. J. Broadmeadow, J. E. Drake, M. Freeman, O. Ghannoum,     !
   ! L. B. Hutley, J. W. Kelly, K. Kikuzawa, P. Kolari, K. Koyama, J.-M. Limousin, P. Meir,!
   ! A. C. Lola da Costa, T. N. Mikkelsen, N. Salinas, W. Sun, and L. Wingate. 2015.       ! 
   ! Optimal stomatal behaviour around the world. Nature Clim. Change 5:459-464. (L15)     !
   !                                                                                       !
   ! Manzoni S, Vico G, Katul G, Fay PA, Polley W, Palmroth S , Porporato A. 2011.         !
   !    Optimizing stomatal conductance for maximum carbon gain under water stress: a      !
   !    meta-analysis across plant functional types and climates. Funct. Ecol., 25:        !
   !    456-467. doi:10.1111/j.1365-2435.2010.01822.x (M11).                               !
   !                                                                                       !
   ! Niinemets U. 2001. Global-scale climatic controls of leaf dry mass per area, density, !
   !    and thickness in trees and shrubs. Ecology, 82: 453-469.                           !
   !    doi:10.1890/0012-9658(2001)082[0453:GSCCOL]2.0.CO;2 (N01).                         !
   !                                                                                       !
   ! Powell TL, Wheeler JK, de Oliveira AAR, da Costa ACL, Saleska SR, Meir P , Moorcroft  !
   !    PR. 2017. Differences in xylem and leaf hydraulic traits explain differences in    !
   !    drought tolerance among mature Amazon rainforest trees. Glob. Change Biol., 23:    !
   !    4280-4293. doi:10.1111/gcb.13731 (P17).                                            !
   !                                                                                       !
   ! Sack L, Cowan PD, Jaikumar N , Holbrook NM. 2003. The 'hydrology'  of leaves:         !
   !    coordination of structure and function in temperate woody species. Plant Cell      !
   !    Environ., 26: 1343-1356. doi:10.1046/j.0016-8025.2003.01058.x (S03).               !
   !                                                                                       !
   ! Scholz FG, Phillips NG, Bucci SJ, Meinzer FC , Goldstein G. 2011. Hydraulic           !
   !    capacitance: Bio- physics and functional significance of internal water sources in !
   !    relation to tree size. In: Size- and Age- Related Changes in Tree Structure and    !
   !    Function (eds. Meinzer FC., Lachenbruch B. & Dawson TE.). Springer Netherlands,    !
   !    Dordrecht, pp. 341-361. doi:10.1007/978-94-007-1242-3 13 (S11).                    !
   !                                                                                       !!
   ! Steppe K, De Pauw DJW, Lemeur R , Vanrolleghem PA. 2006. A mathematical model linking !
   !    tree sap flow dynamics to daily stem diameter fluctuations and radial stem growth. !
   !    Tree Physiol., 26: 257-273. doi:10.1093/treephys/26.3.257 (S06).                   !
   !                                                                                       !
   ! Xu X, Medvigy D, Powers JS, Becknell JM , Guan K. 2016. Diversity in plant hydraulic  !
   !    traits explains seasonal and inter-annual variations of vegetation dynamics in     !
   !    seasonally dry tropical forests. New Phytol., 212: 80-95. doi:10.1111/nph.14009    !
   !    (X16).                                                                             !
   !---------------------------------------------------------------------------------------!


   !----- Bounded wood density. -----------------------------------------------------------!
   rho_bnd(:) = max(rho_min,min(rho_max,rho(:)))
   !---------------------------------------------------------------------------------------!


   !----- Flag for selecting tropical trees. ----------------------------------------------!
   is_troptree(:) = is_tropical(:)        .and. (.not. is_grass(:)) .and.                  &
                    (.not. is_conifer(:)) .and. (.not. is_liana(:))
   !---------------------------------------------------------------------------------------!


   !----- This is kind of arbitrary, total hydraulic path is 50% more than tree height. ---!
   vessel_curl_factor(:)  = 1.5 
   !---------------------------------------------------------------------------------------!




   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !---------------------------------------------------------------------------------------!
   !     Plant hydraulic traits are from C16.  The data are from tropical species. Need to !
   ! update them for temperate PFTs                                                        !
   !                                                                                       !
   !     Don't panic on the number of parameters here.  Most of them are used to calculate !
   ! a few traits that control plant hydraulic calculations.   These traits are listed in   !
   ! the next section.                                                                     !
   !---------------------------------------------------------------------------------------!
   LMA(:) = 1.e3 * C2B / SLA(:)
   !---------------------------------------------------------------------------------------!

   !----- Leaf osmotic potential at saturation [m]. ---------------------------------------!
   leaf_psi_osmotic(:) = MPa2m * (-0.04 - 1.51 * rho_bnd(:) - 0.0067 * LMA(:))
   !---------------------------------------------------------------------------------------!

   !----- Leaf bulk elastic modulus [MPa]. ------------------------------------------------!
   leaf_elastic_mod(:) = 2.5 + 37.5 / (1. + exp(-8. * rho_bnd(:) + 5.7))
   !---------------------------------------------------------------------------------------!

   !----- Leaf minimum relative water content, or residual fraction. ----------------------!
   leaf_rwc_min(:) = 0.01 * leaf_elastic_mod(:) + 0.17
   !---------------------------------------------------------------------------------------!

   !----- Leaf turgor loss point [m]. -----------------------------------------------------!
   leaf_psi_tlp(:) = leaf_elastic_mod(:)                                                   &
                      * (leaf_psi_osmotic(:) / MPa2m)                                      &
                      / (leaf_elastic_mod(:) + leaf_psi_osmotic(:) / MPa2m)                &
                      * MPa2m
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Leaf water content at saturation (rwc = 1.) [kg H2O / kg biomass].                 !
   !                                                                                       !
   !    First calculate leaf_density from Fig. 4 in N01.  The  calculated                  !
   ! leaf_water_sat is comparable to measurements from K09.                                !
   !---------------------------------------------------------------------------------------!
   !leaf_density(:)   = max(0.1 * 1.e3, (leaf_elastic_mod(:) - 2.03) / 25.4 * 1.e3)
   !leaf_water_sat(:) = (-2.32e4 / LMA(:) + 782.)                                           &
   !                     * (1. / (-0.21 * log(1.e4 / LMA(:)) + 1.43) - 1.)                  &
   !                     / leaf_density(:)
   ! TODO: remove above after test
   ! The previous two equations will yield large values for leaf_water_sat per leaf area (>0.2
   ! kg/m2 while the observed is ~ 0.12), probably due to large uncertainties in both equations

   ! Now we use data from Powers and Tiffin 2009 to calculate leaf_water_sat from wood density
   ! This will generate much reasonable values for leaf_water_sat. We add a correcting factor of 1.3
   ! so that the average leaf water content from Powers & Tiffin is the same as K09
   leaf_water_sat(:) = 1.3 * 2.57 * exp(-0.94 * rho_bnd(:)) ! R2 = 0.24
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Sapwood traits.                                                                    !
   !---------------------------------------------------------------------------------------!
   !----- Sapwood osmotic potential at saturation [m]. ------------------------------------!
   wood_psi_osmotic(:) = (0.52 - 4.16 * rho_bnd(:)) * MPa2m
   !---------------------------------------------------------------------------------------!

   !----- Sapwood bulk elastic modulus [MPa]. ---------------------------------------------!
   wood_elastic_mod(:) = sqrt(1.02 * exp(8.5 * rho_bnd(:)) - 2.89)
   !---------------------------------------------------------------------------------------!

   !----- Sapwood minimum relative water content, or residual fraction. -------------------!
   rwc_tlp_wood(:) = 1.00 - (1.00 - 0.75 * rho_bnd(:)) / (2.74 + 2.01 * rho_bnd(:))
   wood_rwc_min(:) = wood_elastic_mod(:) * (1.00 - rwc_tlp_wood(:) - f_cap)             &
                        / (wood_psi_osmotic(:) / MPa2M * (1.00 - f_cap) ) + 1.00 - f_cap 

   ! note that there is a typo in the original equation, the last epsilon_x before fcap is extra
   !---------------------------------------------------------------------------------------!

   !----- Wood water content at saturation (rwc = 1.) [kg H2O / kg biomass]. --------------!
   wood_water_sat(:) = (1. - rho_bnd(:) / 1.53) * wdns / (rho_bnd(:) * 1.e3)
   !---------------------------------------------------------------------------------------!

   !----- Bark water content at saturation (rwc = 1.) [kg H2O / kg biomass]. --------------!
   bark_water_sat(:) = wood_water_sat(:)
   !---------------------------------------------------------------------------------------!

   !----- Wood turgor loss point. ---------------------------------------------------------!
   wood_psi_tlp(:) = wood_elastic_mod(:) * (wood_psi_osmotic(:) / MPa2m)                   &
                   / (wood_elastic_mod(:) + wood_psi_osmotic(:) / MPa2m)                   &
                   * MPa2m
   !---------------------------------------------------------------------------------------!
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!






   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !     Key traits that drive plant hydrodynamics.  Again based on C16 (see above).       !
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !    Calculate capacitance based on linearizing the water retention curve from B12.     !
   ! Assume water retention curve is linear from 0. to 4 * turgor loss point               !
   ! [kg H2O / kg biomass / m].                                                            !
   !---------------------------------------------------------------------------------------!
   leaf_water_cap(:) = (1. - leaf_rwc_min(:))                                              &
                     * (1. - leaf_psi_osmotic(:) / (4. * leaf_psi_tlp(:)) )                &
                     * leaf_water_sat(:) / (4. * abs(leaf_psi_tlp(:)) )

   wood_water_cap(:) = (1. - wood_rwc_min(:))                                              &
                     * (1. - wood_psi_osmotic(:) / (4. * wood_psi_tlp(:)) )                &
                     * wood_water_sat(:) / (4. * abs(wood_psi_tlp(:)) )

   !---------------------------------------------------------------------------------------!



   !----- Wood P50 [m]. -------------------------------------------------------------------!
   wood_psi50(:) = (-1.09 - (3.57 * rho_bnd(:) ** 1.73)) * MPa2m 
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   ! This is only an estimate. 2.4 is Q10, converting to Vcmax_25. The 4.1 factor is a     !
   ! conversion factor from Vcmax to Amax at ~25degC.                                      !
   !---------------------------------------------------------------------------------------!
   !Amax_25(:) = Vm0(:) * 2.4 / 4.1 ! umol/m2/s
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Sapwood maximum conductivity [kg/m2/s].  This is estimated from Figure S2.2       !
   ! of C16.                                                                               !
   !---------------------------------------------------------------------------------------!
   !wood_Kmax(:)  = exp(2.11 - 20.05 * rho_bnd(:) / Amax_25(:)) / MPa2m 
   ! TODO: remove above after test
   wood_Kmax(:) = exp( 2.32 - 2.27 * rho_bnd(:)                                            &
                     - 0.48 * log(-wood_psi50(:) / MPa2m)                                  &
                     + 0.5 * 0.89) / MPa2m
   ! from analysis of the Gleason et al. and Xu et al. data.
   ! Again this makes more sense....

   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     C16 only reports the slope of PLC at psi50.  We determine the slope               !
   ! (- wood_Kexp / (4 * wood_psi50)) by calculating the derivatives of the PLC            !
   ! function. We then back-calculate wood_Kexp from the slope here.                       !
   !---------------------------------------------------------------------------------------!
   wood_Kexp(:)  = 0.544 * 4. * (-wood_psi50(:) / MPa2m) ** (-0.17)
   !---------------------------------------------------------------------------------------!
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!


   !---------------------------------------------------------------------------------------!
   !      Overwrite some parameters if PLANT_HYDRO_SCHEME is 2, using the meta-analysis    !
   ! from X16.  This is tested for tropical tree PFTs only, so we only update the values   !
   ! for these trees.                                                                      !
   !---------------------------------------------------------------------------------------!
   select case (plant_hydro_scheme)
   case (2)
      !------------------------------------------------------------------------------------!
      !      Capacitance is estimated from S11.  Since the capacitance is treated as a     !
      ! constant in the model, it should be way smaller than lab/field measured            !
      ! capacitance (S03; S06).  Thus, here Cap_leaf is multiplied by 1/2 and Cap_stem     !
      ! is multipled by 1/3.                                                               !
      !------------------------------------------------------------------------------------!
      leaf_water_cap(:) = merge( 3.e-3 * SLA(:) / C2B / MPa2m / 2.                         &
                               , leaf_water_cap(:)                                         &
                               , is_troptree(:) )
      wood_water_cap(:) = merge( min(400.,max(50., -700. * (rho_bnd(:) - 0.3) + 400.))     &
                               / (rho_bnd(:) * 1.e3) / MPa2m / 3.                          &
                               , wood_water_cap(:)                                         &
                               , is_troptree(:) )
      !------------------------------------------------------------------------------------!




      !----- Copied from default values (G07). --------------------------------------------!
      leaf_water_sat(:) = merge(1.85,leaf_water_sat(:),is_troptree(:))
      wood_water_sat(:) = merge(0.70,wood_water_sat(:),is_troptree(:))
      bark_water_sat(:) = merge(0.70,bark_water_sat(:),is_troptree(:))
      !------------------------------------------------------------------------------------!



      !----- Set some rwc_min so that psi_min makes sense. --------------------------------!
      leaf_rwc_min(:) = merge(0.5 ,leaf_rwc_min(:),is_troptree(:))
      wood_rwc_min(:) = merge(0.05,wood_rwc_min(:),is_troptree(:))
      !------------------------------------------------------------------------------------!



      !----- Additional parameters. -------------------------------------------------------!
      leaf_psi_tlp(:) = merge( ( - 4.59 + 0.62 * log(SLA(:))                               &
                               - 1.15 * log(rho_bnd(:)) ) * MPa2m                          &
                             , leaf_psi_tlp(:)                                             &
                             , is_troptree(:) )
      wood_psi50  (:) = merge( (-3. * rho_bnd(:) - 0.599) * MPa2m                          &
                             , wood_psi50(:)                                               &
                             , is_troptree(:) )
      wood_Kmax   (:) = merge( exp(-2.455 * rho_bnd(:) + 2.348 + 0.5 * 0.6186) / MPa2m     &
                             , wood_Kmax   (:)                                             &
                             , is_troptree(:) )
      wood_Kexp   (:) = merge(4.,wood_Kexp(:),is_troptree(:))
      !------------------------------------------------------------------------------------!
   case (0)
      !------------------------------------------------------------------------------------!
      !     Use default values for back-compatibility.                                     !
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! Leaf Water:oven-dry mass ratio.                                                    !
      !                                                                                    !
      ! Tropical  -- Average well-watered values of wd from K09.                           !
      ! Temperate -- check with MCD.  Original value was 1.5, from G07 but it has been     !
      !              replaced by 2.5.  Perhaps to account for the C:B ratio, but if this   !
      !              is the case, it is accounting for it twice.                           !
      !------------------------------------------------------------------------------------!
      leaf_water_sat (:) = merge(1.85,2.50,is_tropical(:))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! Wood and bark water:oven-dry mass ratio.                                           !
      !                                                                                    !
      ! Tropical  -- Numbers were obtained from P14 (Table S1).  P14's WWC and BWC were    !
      !              converted to water:oven-dry ratio (XWDR = XWC / (1 - XWC)).  Their    !
      !              data suggest that water content depends on wood density but were      !
      !              indistinguishable between moist and dry forests.  Values for tropical !
      !              trees follow the fitted curve using SMA.  This is used only when      !
      !              ECONOMICS_SCHEME=1.                                                   !
      !                                                                                    !
      ! Temperate -- Use values from FPL10.                                                !
      !------------------------------------------------------------------------------------!
      select case (economics_scheme)
      case (1)
         do ipft=1,n_pft
            if (is_grass(ipft) .or. is_conifer(ipft) .or. (.not. is_tropical(ipft))) then
               wood_water_sat(ipft) = 0.7
               bark_water_sat(ipft) = 0.7
            else
               wood_water_sat(ipft) = exp(1.5018230 - 3.137476 * rho_bnd(ipft))
               bark_water_sat(ipft) = exp(1.9892840 - 3.174365 * rho_bnd(ipft))
            end if
         end do
      case default
         wood_water_sat(:) = 0.7
         bark_water_sat(:) = 0.7
      end select
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Equivalent leaf/wood minimum water potential calculated from leaf_rwc with the    !
   ! assumption of constant capacitance [m].                                               !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      !----- Make minimum water potential consistent with relative water content. ---------!
      call rwc2psi(leaf_rwc_min(ipft),wood_rwc_min(ipft),ipft,leaf_psi_min(ipft)        &
                  ,wood_psi_min(ipft))
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!


   !----- Parameters related with stomatal conductance, estimated from L15 and M11. ---------------!

   ! stoma_lambda has a unit of mol CO2 / mol H2O
   ! 2.85e-2 / (3.77) ** 2 will generate most reasonal seasonality in GPP and ET over tropical
   ! forests. Therefore, we treat 2.85e-2 is equivalent to a g1 of 1 in L11
   ! lambda ~ 1/sqrt(g1)
   ! therefore, we convert PFT average g1 values in L11 to stoma_lambda
   stoma_lambda(:) = merge( merge( 1.62                     &  ! C4 grass
                                 , 4.5                      &  ! C3 grass
                                 , photosyn_pathway(:) == 4)&
                           ,merge( 3.77                     &  ! tropical trees
                                 , merge( 2.35              & ! conifer trees
                                        , 4.64              & ! temperate/boreal deciduous
                                        ,is_conifer(:))     &
                                 ,is_tropical(:))           &
                           ,is_grass(:))

   stoma_lambda = 2.85e-2 / stoma_lambda ** 2


   ! stoma_beta is based on Table 2 in M11
   stoma_beta(:) = merge( -1.1                              & ! Grass
                          ,merge( -0.13                     & ! tropical trees
                                , merge( -1.06              & ! conifer trees
                                       , -1.34              & ! temperate/boreal deciduous
                                       ,is_conifer(:))     &
                                ,is_tropical(:))           &
                          ,is_grass(:)) / MPa2m
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Use parameters from P17.                                                           !
   !    For tropical intolerant, use -1.67 * MPa2m, 3.                                     !
   !    For tropical tolerant, use -2.83 * MPa2m, 3.5                                      !
   !---------------------------------------------------------------------------------------!
   stoma_psi_b(:) = leaf_psi_tlp(:)   ! default
   stoma_psi_c(:) = 3.
   !---------------------------------------------------------------------------------------!


   return
end subroutine init_pft_hydro_params
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!  SUBROUTINE: INIT_PFT_SPHEAT_PARAMS   
!> \brief This subroutine initialises specific heat parameters for plants.
!> \details This subroutine initialises specific heat of leaves, wood, and bark, and it
!>          must be called after init_pft_hydro_params.
!> \author Marcos Longo, 10 September 2019
!
! References for this subroutine:
!
! Forest Products Laboratory. Wood handbook - wood as an engineering material. General
!    Technical Report FPL-GTR-190, U.S. Department of Agriculture, Madison, WI, 2010.
!    doi:10.2737/FPL-GTR-190 (FPL10)
!
! Gu, L., T. Meyers, S. G. Pallardy, P. J. Hanson, B. Yang, M. Heuer, K. P. Hosman,
!    Q. Liu, J. S. Riggs, D. Sluss, and S. D. Wullschleger. Influences of biomass heat
!    and biochemical energy storages on the land surface fluxes and radiative temper-
!    ature. J. Geophys. Res., 112(D2):D02107, Jan 2007. doi:10.1029/2006JD007425 (G07)
!
! Jones, H. G. Plants and Microclimate: A quantitative approach to environmental plant
!    physiology. Cambridge Univ. Press, Cambridge, UK, 3rd edition, Jan 2014.
!    doi:10.1017/CBO9780511845727 (J14)
!
! Kursar, T. A., B. M. J. Engelbrecht, A. Burke, M. T. Tyree, B. EI Omari,
!    and J. P. Giraldo. Tolerance to low leaf water status of tropical tree seedlings
!    is related to drought performance and distribution. Funct. Ecol., 23(1):93-102,
!    Feb 2009. doi:10.1111/j.1365-2435.2008.01483.x (K09)
!
! Poorter, L., A. McNeil, V.-H. Hurtado, H. H. T. Prins, and F. E. Putz. Bark traits
!    and life-history strategies of tropical dry- and moist forest trees.
!    Funct. Ecol., 28(1):232-242, Feb 2014. doi:10.1111/1365-2435.12158 (P14)
!------------------------------------------------------------------------------------------!
subroutine init_pft_spheat_params()
   use consts_coms    , only : t00                  ! ! intent(in)
   use pft_coms       , only : wood_water_sat       & ! intent(in)
                             , bark_water_sat       & ! intent(in)
                             , c_grn_leaf_dry       & ! intent(out)
                             , c_ngrn_wood_dry      & ! intent(out)
                             , c_ngrn_bark_dry      & ! intent(out)
                             , delta_c_wood         & ! intent(out)
                             , delta_c_bark         ! ! intent(out)

   implicit none
   !---------------------------------------------------------------------------------------!
   !      Local parameters used to define specific heat of wood and bark, following the    !
   ! tref   -- Reference temperature (Kelvin) for specific heat properties.                !
   ! wdr_fs -- Water:Dry ratio at fiber saturation.  This is the maximum moisture that     !
   !           affects specific heat due to water-wood bond according to FPL10.  We        !
   !           currently assume 0.30, the suggested value.                                 !
   !                                                                                       !
   ! Reference:                                                                            !
   !                                                                                       !
   ! Forest Products Laboratory. Wood handbook - wood as an engineering material. General  !
   !    Technical Report FPL-GTR-190, U.S. Department of Agriculture, Madison, WI, 2010.   !
   !    doi:10.2737/FPL-GTR-190 (FPL10)                                                    !
   !---------------------------------------------------------------------------------------!
   real, parameter :: tref   = t00 + 15.
   real, parameter :: wdr_fs = 0.30
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Specific heat of dry materials (J kg-1 K-1).                                       !
   !---------------------------------------------------------------------------------------!
   !----- Leaves, value from G07. ---------------------------------------------------------!
   c_grn_leaf_dry(:) = 3218.0
   !----- Wood and bark, values from FPL10. -----------------------------------------------!
   c_ngrn_wood_dry(:) = 103.1 + 3.867 * tref
   c_ngrn_bark_dry(:) = 103.1 + 3.867 * tref
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    The oven-dry wood heat capacity and the specific heat correction for water-wood    !
   ! bonding come both from FPL10, previous version cited by G07.  Following FPL10, the    !
   ! maximum moisture to affect the water-wood bonding is the fiber saturation, above      !
   ! which additional water is considered free water.                                      !
   !---------------------------------------------------------------------------------------!
   where (wood_water_sat(:) > wdr_fs)
      delta_c_wood(:) = 1.e5 * wdr_fs                                                      &
                      * ( - 0.06191 + 2.36e-4 * tref - 1.33e-2 * wdr_fs )
   elsewhere
      delta_c_wood(:) = 1.e5 * wood_water_sat(:)                                           &
                      * ( - 0.06191 + 2.36e-4 * tref - 1.33e-2 * wood_water_sat(:) )
   end where
   where (bark_water_sat(:) > wdr_fs)
      delta_c_bark(:) = 1.e5 * wdr_fs                                                      &
                      * ( - 0.06191 + 2.36e-4 * tref - 1.33e-2 * wdr_fs )
   elsewhere
      delta_c_bark(:) = 1.e5 * bark_water_sat(:)                                           &
                      * ( - 0.06191 + 2.36e-4 * tref - 1.33e-2 * bark_water_sat(:) )
   end where
   !---------------------------------------------------------------------------------------!


   return
end subroutine init_pft_spheat_params
!==========================================================================================!
!==========================================================================================!








!==========================================================================================!
!==========================================================================================!
!   This subroutine sets up some PFT and leaf dependent properties.                        !
!------------------------------------------------------------------------------------------!
subroutine init_pft_phen_params()
   use phenology_coms , only : iphen_scheme         ! ! intent(in)
   use ed_misc_coms   , only : igrass               ! ! intent(in)
   use pft_coms       , only : phenology            & ! intent(out)
                             , is_grass             & ! intent(in)
                             , is_conifer           & ! intent(in)
                             , is_tropical          & ! intent(in)
                             , high_psi_threshold   & ! intent(out)
                             , low_psi_threshold    & ! intent(out)
                             , storage_reflush_times& ! intent(out)
                             , leaf_shed_rate       & ! intent(out)
                             , leaf_grow_rate       ! ! intent(out)
   use ed_max_dims    , only : n_pft                ! ! intent(in)
   implicit none

   !----- Local variables. ----------------------------------------------------------------!
   integer            :: ipft
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Tree phenology is the same for both cases, but in the new grass allometry they    !
   ! must be evergreens.                                                                   !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      if (is_conifer(ipft)) then  ! Conifers. Currently they are always evergreen
         phenology(ipft) = 0
      elseif (is_grass(ipft) .and. igrass == 1) then ! New grasses must be evergreen
         phenology(ipft) = 0
      elseif (.not. (is_tropical(ipft) .or. is_grass(ipft)) ) then ! Cold deciduous
         phenology(ipft) = 2
      else
         !----- Lianas, Tropical broadleaf trees, and old-scheme (aka bonsai) grasses. ----!
         select case (iphen_scheme)
         case (-1) ! Assume that they are all evergreen
            phenology(ipft) = 0
         case (0,1) ! Old drought-deciduous scheme 
            phenology(ipft) = 1
         case (2)   ! New drought-deciduous scheme
            phenology(ipft) = 4
         case (3)
            !------------------------------------------------------------------------------!
            !     Environmental Light phenology scheme.                                    !
            !                                                                              !
            ! Kim, Y., R. G. Knox, M. Longo, D. Medvigy, L. R. Hutyra, E. H. Pyle,         !
            !    S. C. Wofsy, R. L. Bras, and P. R. Moorcroft. Seasonal carbon dynamics    !
            !    and water fluxes in an Amazon rainforest. Glob. Change Biol.,             !
            !    18(4):1322-1334, Apr 2012. doi:10.1111/j.1365-2486.2011.02629.x.          !
            !------------------------------------------------------------------------------!
            phenology(ipft) = 3
            !------------------------------------------------------------------------------!
         case (4)
            !------------------------------------------------------------------------------!
            !   Xiangtao Xu's plant-hydraulics-driven drought phenology.                   !
            !------------------------------------------------------------------------------!
            phenology(ipft) = 5
            !------------------------------------------------------------------------------!
         case (5)
            !------------------------------------------------------------------------------!
            !   Combined light and plant-hydraulics driven drought phenology.              !
            !------------------------------------------------------------------------------!
            phenology(ipft) = 6
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Storage Reserves to reflush leaves and fine roots                                 !
   ! set it to be 2. if run with hydraulics-driven drought phenology                       !
   ! set it to be 0 otherwise (default)                                                    !
   !---------------------------------------------------------------------------------------!
   storage_reflush_times(:) = merge(2.0,0.0,phenology(:) == 5)

   !---------------------------------------------------------------------------------------!
   !     Phenology-related parameters for phenology(ipft) = 5.                             !
   !---------------------------------------------------------------------------------------!
   high_psi_threshold(:) = 10
   low_psi_threshold (:) = 10
   leaf_shed_rate    (:) = 1./20.
   leaf_grow_rate    (:) = 1./20.
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_pft_phen_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine sets some reproduction-related parameters.                            !
!------------------------------------------------------------------------------------------!
subroutine init_pft_repro_params()

   use pft_coms      , only : hgt_min            & ! intent(in)
                            , hgt_max            & ! intent(in)
                            , is_tropical        & ! intent(in)
                            , is_conifer         & ! intent(in)
                            , is_liana           & ! intent(in)
                            , is_grass           & ! intent(in)
                            , seed_rain          & ! intent(out)
                            , r_bang             & ! intent(out)
                            , r_fract            & ! intent(out)
                            , r_cv50             & ! intent(out)
                            , st_fract           & ! intent(out)
                            , nonlocal_dispersal & ! intent(out)
                            , repro_min_h        ! ! intent(out)
   use ed_max_dims   , only : n_pft              & ! intent(in)
                            , undef_real         ! ! intent(in)
   use ed_misc_coms  , only : economics_scheme   ! ! intent(in)
   use consts_coms   , only : onesixth           & ! intent(in)
                            , onethird           ! ! intent(in)
   use phenology_coms, only : repro_scheme       ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer                 :: ipft
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Allocation to storage (amount that is save in each month, not going to           !
   ! reproduction or growth).  Currently only tropical trees with ECONOMICS_SCHEME=1       !
   ! maintain storage pools.  The fraction is tuned to make the pool somewhat closer to    !
   ! the typical storage/biomass ratio (e.g. MV16).                                        !
   !                                                                                       !
   ! Martinez-Vilalta, J, A. Sala, D. Asensio, L. Galiano, G. Hoch, S. Palacio,            !
   !    F. I. Piper, and F. Lloret. Dynamics of non-structural carbohydrates in            !
   !    terrestrial plants: a global synthesis. Ecol. Monogr., 86(4):495-516, Nov 2016.    !
   !    doi:10.1002/ecm.1231.                                                              !
   !---------------------------------------------------------------------------------------!
   select case (economics_scheme)
   case (1)
      st_fract(:) = merge(0.0,onesixth                                                     &
                         ,is_grass(:) .or. is_liana(:) .or. (.not. is_tropical(:)))
   case default
      st_fract(:) = 0.0
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Reproduction parameters.                                                         !
   !                                                                                       !
   ! repro_min_h - Minimum height for a PFT to be considered mature for reproduction.      !
   ! r_bang      - Logical variable that decides whether reproduction should follow a      !
   !               "partial/big bang" (same notation as W15) or asymptote.                 !
   ! r_fract     - This is only used when r_bang is .true..  Fraction of carbon balance    !
   !               (bstorage) that is allocated to reproduction given that the plant       !
   !               height is greater than repro_min_h.                                     !
   !                                                                                       !
   !               Details on the default values: for trees, we take repro_min_h as the    !
   !               average value reported in W05 (mind that the study is for tropical      !
   !               trees).  The grass strategy depends on economics_scheme.  The original  !
   !               scheme assumes that they always allocate 30% to growth ("partial        !
   !               bang"), whereas when ECONOMICS_SCHEME=1 the plants will allocate 100%   !
   !               to reproduction once they reach the maximum height (W15's "big bang").  !
   !               The fraction allocated to reproduction for tropical trees was updated   !
   !               to match F10's number, or the default ED-1.0 numbers (0.30, following   !
   !               M01).                                                                   !
   ! r_cv50      - This is only used when r_bang is .false..  For plants with asymptote    !
   !               reproduction (r_bang = .false.), this is the dimensionless curvature    !
   !               parameter.  The reproduction allocation converges to "big bang" when    !
   !               r_cv50 approaches 0, and it shuts down reproduction when r_cv50         !
   !               approaches infinity.                                                    !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   !                                                                                       !
   !  Moorcroft, P. R., G. C. Hurtt, and S. W. Pacala. A method for scaling vegetation     !
   !     dynamics: The Ecosystem Demography model (ED). Ecol. Monogr., 71(4):557-586,      !
   !     Nov 2001. doi:10.1890/0012- 9615(2001)071[0557:AMFSVD]2.0.CO;2.                   !
   !                                                                                       !
   !  Wright, S. J., M. A. Jaramillo, J. Pavon, R. Condit, S. P. Hubbell, and              !
   !     R. B. Foster. Reproductive size thresholds in tropical trees: variation among     !
   !     individuals, species and forests. J. Trop. Ecol., 21(3):307-315, May 2005.        !
   !     doi:10.1017/S0266467405002294. (W05).                                             !
   !                                                                                       !
   !  Fisher, R., N. McDowell, D. Purves, P. Moorcroft, S. Sitch, P. Cox, C. Huntingford,  !
   !    P. Meir, and F. Ian Woodward. Assessing uncertainties in a second-generation       !
   !    dynamic vegetation model caused by ecological scale limitations. New Phytol.,      !
   !    187(3):666--681, Aug 2010. doi:10.1111/j.1469-8137.2010.03340.x.  (F10)            !
   !                                                                                       !
   !  Wenk, E. H., and D. S. Falster. Quantifying and understanding reproductive           !
   !    allocation schedules in plants. Ecol. Evol., 5(23):5521--5538, Nov 2015.           !
   !    doi:10.1002/ece3.1802. (W15)                                                       !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   select case (economics_scheme)
   case (1)
      !----- New parameters. --------------------------------------------------------------!
      repro_min_h(:) = merge(hgt_max(:),                           18.0,is_grass(:))
      r_fract    (:) = merge(       1.0,merge(0.37,0.30,is_tropical(:)),is_grass(:))
      !------------------------------------------------------------------------------------!
   case default
      !----- Original parameters. ---------------------------------------------------------!
      repro_min_h(:) = merge(hgt_min(:),18.0,is_grass(:))
      r_fract    (:) = 0.30
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!
   !      Asymptote parameters.  This is a temporary sensitivity test, these options will  !
   ! be combined soon into a single option (repro_scheme=3).                               !
   !---------------------------------------------------------------------------------------!
   select case (repro_scheme)
   case (3)
      r_cv50(:) = onethird
      r_bang(:) = is_grass(:) .or. is_liana(:) .or. (.not. is_tropical(:))
   case default
      r_cv50(:) = onethird
      r_bang(:) = .true.
   end select
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Fraction of seeds randomly dispersed.                                             !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      if (is_tropical(ipft) .and. is_conifer(ipft)) then
         !----- Sub-tropical needleleaf: assume the same values as pines. -----------------!
         nonlocal_dispersal(ipft) = 0.766
         !---------------------------------------------------------------------------------!
      elseif (is_tropical(ipft) .or. is_grass(ipft)) then
         !----- Tropical trees or grasses. Assume 100% random dispersal. ------------------!
         nonlocal_dispersal(ipft) = 1.00
         !---------------------------------------------------------------------------------!
      else
         !----- Temperate broadleaf trees. ------------------------------------------------!
         select case (ipft)
         case (6:7)   ! Pines. 
            nonlocal_dispersal(ipft) = 0.766
         case (8)     ! Late conifers. 
            nonlocal_dispersal(ipft) = 0.001
         case (9)     ! Early hardwood. 
            nonlocal_dispersal(ipft) = 1.000
         case (10)    ! Mid hardwood. 
            nonlocal_dispersal(ipft) = 0.325
         case (11)    ! Late hardwood. 
            nonlocal_dispersal(ipft) = 0.074
         case default ! This shouldn't happen. 
            nonlocal_dispersal(ipft) = 1.000
         end select
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Seed rain: temporarily set this parameter to undefined.  In case it is not        !
   ! initialised by XML, it will be initialised in init_derived_params_after_xml.          !
   !---------------------------------------------------------------------------------------!
   seed_rain(:) = undef_real
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_pft_repro_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will assign some canopy air related parameters.                       !
!------------------------------------------------------------------------------------------!
subroutine init_can_air_params()
   use consts_coms    , only : onethird              & ! intent(in)
                             , twothirds             & ! intent(in)
                             , onesixth              & ! intent(in)
                             , vonk                  ! ! intent(in)
   use pft_coms       , only : hgt_min               ! ! intent(in)
   use rk4_coms       , only : tiny_offset           ! ! intent(in)
   use canopy_air_coms, only : psim8                 & ! function
                             , psih8                 & ! function
                             , ugbmin                & ! intent(in)
                             , ubmin                 & ! intent(in)
                             , ustmin                & ! intent(in)
                             , gamm                  & ! intent(in)
                             , gamh                  & ! intent(in)
                             , tprandtl              & ! intent(in)
                             , vh2vr                 & ! intent(out)
                             , vh2dh                 & ! intent(out)
                             , ribmax                & ! intent(out)
                             , f_bndlyr_init         & ! intent(out)
                             , leaf_drywhc           & ! intent(out)
                             , leaf_maxwhc           & ! intent(out)
                             , gbhmos_min            & ! intent(out)
                             , gbhmos_min8           & ! intent(out)
                             , veg_height_min        & ! intent(out)
                             , veg_height_min8       & ! intent(out)
                             , minimum_canopy_depth  & ! intent(out)
                             , minimum_canopy_depth8 & ! intent(out)
                             , exar                  & ! intent(out)
                             , covr                  & ! intent(out)
                             , exar8                 & ! intent(out)
                             , ez                    & ! intent(out)
                             , ustmin8               & ! intent(out)
                             , ugbmin8               & ! intent(out)
                             , ubmin8                & ! intent(out)
                             , ez8                   & ! intent(out)
                             , vh2vr8                & ! intent(out)
                             , vh2dh8                & ! intent(out)
                             , cdrag0                & ! intent(out)
                             , cdrag1                & ! intent(out)
                             , cdrag2                & ! intent(out)
                             , cdrag3                & ! intent(out)
                             , pm0                   & ! intent(out)
                             , c1_m97                & ! intent(out)
                             , c2_m97                & ! intent(out)
                             , c3_m97                & ! intent(out)
                             , kvwake                & ! intent(out)
                             , alpha_m97             & ! intent(out)
                             , alpha_mw99            & ! intent(out)
                             , gamma_mw99            & ! intent(out)
                             , nu_mw99               & ! intent(out)
                             , infunc                & ! intent(out)
                             , cs_dense0             & ! intent(out)
                             , gamma_clm4            & ! intent(out)
                             , cdrag08               & ! intent(out)
                             , cdrag18               & ! intent(out)
                             , cdrag28               & ! intent(out)
                             , cdrag38               & ! intent(out)
                             , pm08                  & ! intent(out)
                             , c1_m978               & ! intent(out)
                             , c2_m978               & ! intent(out)
                             , c3_m978               & ! intent(out)
                             , kvwake8               & ! intent(out)
                             , alpha_m97_8           & ! intent(out)
                             , alpha_mw99_8          & ! intent(out)
                             , gamma_mw99_8          & ! intent(out)
                             , nu_mw99_8             & ! intent(out)
                             , infunc_8              & ! intent(out)
                             , cs_dense08            & ! intent(out)
                             , bl79                  & ! intent(out)
                             , csm                   & ! intent(out)
                             , csh                   & ! intent(out)
                             , dl79                  & ! intent(out)
                             , beta_s                & ! intent(out)
                             , abh91                 & ! intent(out)
                             , bbh91                 & ! intent(out)
                             , cbh91                 & ! intent(out)
                             , dbh91                 & ! intent(out)
                             , ebh91                 & ! intent(out)
                             , fbh91                 & ! intent(out)
                             , cod                   & ! intent(out)
                             , bcod                  & ! intent(out)
                             , fm1                   & ! intent(out)
                             , ate                   & ! intent(out)
                             , atetf                 & ! intent(out)
                             , beta_vs               & ! intent(out)
                             , chim                  & ! intent(out)
                             , chih                  & ! intent(out)
                             , zetac_um              & ! intent(out)
                             , zetac_uh              & ! intent(out)
                             , zetac_sm              & ! intent(out)
                             , zetac_sh              & ! intent(out)
                             , zetac_umi             & ! intent(out)
                             , zetac_uhi             & ! intent(out)
                             , zetac_smi             & ! intent(out)
                             , zetac_shi             & ! intent(out)
                             , zetac_umi16           & ! intent(out)
                             , zetac_uhi13           & ! intent(out)
                             , psimc_um              & ! intent(out)
                             , psihc_uh              & ! intent(out)
                             , zd98_a                & ! intent(out)
                             , zd98_b                & ! intent(out)
                             , zd98_emax             & ! intent(out)
                             , bl798                 & ! intent(out)
                             , csm8                  & ! intent(out)
                             , csh8                  & ! intent(out)
                             , dl798                 & ! intent(out)
                             , beta_s8               & ! intent(out)
                             , gamm8                 & ! intent(out)
                             , gamh8                 & ! intent(out)
                             , ribmax8               & ! intent(out)
                             , tprandtl8             & ! intent(out)
                             , abh918                & ! intent(out)
                             , bbh918                & ! intent(out)
                             , cbh918                & ! intent(out)
                             , dbh918                & ! intent(out)
                             , ebh918                & ! intent(out)
                             , fbh918                & ! intent(out)
                             , cod8                  & ! intent(out)
                             , bcod8                 & ! intent(out)
                             , fm18                  & ! intent(out)
                             , ate8                  & ! intent(out)
                             , atetf8                & ! intent(out)
                             , beta_vs8              & ! intent(out)
                             , chim8                 & ! intent(out)
                             , chih8                 & ! intent(out)
                             , zetac_um8             & ! intent(out)
                             , zetac_uh8             & ! intent(out)
                             , zetac_sm8             & ! intent(out)
                             , zetac_sh8             & ! intent(out)
                             , zetac_umi8            & ! intent(out)
                             , zetac_uhi8            & ! intent(out)
                             , zetac_smi8            & ! intent(out)
                             , zetac_shi8            & ! intent(out)
                             , zetac_umi168          & ! intent(out)
                             , zetac_uhi138          & ! intent(out)
                             , psimc_um8             & ! intent(out)
                             , psihc_uh8             & ! intent(out)
                             , zd98_a8               & ! intent(out)
                             , zd98_b8               & ! intent(out)
                             , zd98_emax8            & ! intent(out)
                             , aflat_turb            & ! intent(out)
                             , aflat_lami            & ! intent(out)
                             , bflat_turb            & ! intent(out)
                             , bflat_lami            & ! intent(out)
                             , nflat_turb            & ! intent(out)
                             , nflat_lami            & ! intent(out)
                             , mflat_turb            & ! intent(out)
                             , mflat_lami            & ! intent(out)
                             , ocyli_turb            & ! intent(out)
                             , ocyli_lami            & ! intent(out)
                             , acyli_turb            & ! intent(out)
                             , acyli_lami            & ! intent(out)
                             , bcyli_turb            & ! intent(out)
                             , bcyli_lami            & ! intent(out)
                             , ncyli_turb            & ! intent(out)
                             , ncyli_lami            & ! intent(out)
                             , mcyli_turb            & ! intent(out)
                             , mcyli_lami            & ! intent(out)
                             , aflat_turb8           & ! intent(out)
                             , aflat_lami8           & ! intent(out)
                             , bflat_turb8           & ! intent(out)
                             , bflat_lami8           & ! intent(out)
                             , nflat_turb8           & ! intent(out)
                             , nflat_lami8           & ! intent(out)
                             , mflat_turb8           & ! intent(out)
                             , mflat_lami8           & ! intent(out)
                             , ocyli_turb8           & ! intent(out)
                             , ocyli_lami8           & ! intent(out)
                             , acyli_turb8           & ! intent(out)
                             , acyli_lami8           & ! intent(out)
                             , bcyli_turb8           & ! intent(out)
                             , bcyli_lami8           & ! intent(out)
                             , ncyli_turb8           & ! intent(out)
                             , ncyli_lami8           & ! intent(out)
                             , mcyli_turb8           & ! intent(out)
                             , mcyli_lami8           & ! intent(out)
                             , ggsoil0               & ! intent(out)
                             , kksoil                & ! intent(out)
                             , ggsoil08              & ! intent(out)
                             , kksoil8               ! ! intent(out)
   implicit none
   !----- External functions. -------------------------------------------------------------!
   real   , external :: cbrt
   real   , external :: sngloff
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Minimum leaf water content to be considered.  Values smaller than this will be     !
   ! flushed to zero.  This value is in kg/[m2 tree], so it will be scaled by (LAI+WAI)    !
   ! where needed be.                                                                      !
   !---------------------------------------------------------------------------------------!
   leaf_drywhc = 5.e-4 * leaf_maxwhc
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Variables to define the vegetation aerodynamic conductance.  They are currently  !
   ! not PFT dependent.                                                                    !
   !---------------------------------------------------------------------------------------!
   gbhmos_min  = 1.e-9
   gbhmos_min8 = dble(gbhmos_min)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! veg_height_min       - This is the minimum vegetation height allowed [m].  Vegetation !
   !                        height is used to calculate drag coefficients and patch        !
   !                        roughness.                                                     !
   ! minimum_canopy_depth - This is the minimum canopy depth allowed [m].  Canopy depth    !
   !                        is used to calculate the heat and moisture storage capacity in !
   !                        the canopy air space.                                          !
   !---------------------------------------------------------------------------------------!
   veg_height_min        = minval(hgt_min)
   minimum_canopy_depth  = 5.0  ! alternative: minval(hgt_min)
   !---------------------------------------------------------------------------------------!

   !----- This is the dimensionless exponential wind atenuation factor. -------------------!
   exar  = 2.5
   !---------------------------------------------------------------------------------------!

   !----- This is the scaling factor of tree area index (not sure if it is used...) -------!
   covr = 2.16
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Parameters for surface layer models.                                             !
   !---------------------------------------------------------------------------------------!
   !----- Vegetation roughness:vegetation height ratio. -----------------------------------!
   vh2vr    = 0.13
   !----- Displacement height:vegetation height ratio. ------------------------------------!
   vh2dh    = 0.63
   !---------------------------------------------------------------------------------------!




   !----- Louis (1979) model. -------------------------------------------------------------!
   bl79        = 5.0    ! b prime parameter
   csm         = 7.5    ! C* for momentum (eqn. 20, not co2 char. scale)
   csh         = 5.0    ! C* for heat (eqn.20, not co2 char. scale)
   dl79        = 5.0    ! ???
   !----- Oncley and Dudhia (1995) model. -------------------------------------------------!
   beta_s       = 5.0          ! Beta
   !----- Beljaars and Holtslag (1991) model. ---------------------------------------------!
   abh91       = -1.00         ! -a from equation  (28) and (32)
   bbh91       = -twothirds    ! -b from equation  (28) and (32)
   cbh91       =  5.0          !  c from equations (28) and (32)
   dbh91       =  0.35         !  d from equations (28) and (32)
   ebh91       = -twothirds    ! - factor multiplying a*zeta in equation (32)
   fbh91       =  1.50         ! exponent in equation (32)
   cod         = cbh91/dbh91   ! c/d
   bcod        = bbh91 * cod   ! b*c/d
   fm1         = fbh91 - 1.0   ! f-1
   ate         = abh91 * ebh91 ! a * e
   atetf       = ate   * fbh91 ! a * e * f
   !----- Similar to CLM (2004), but with different phi_m for very unstable case. ---------!
   zetac_um    = -1.5
   zetac_uh    = -0.5
   zetac_sm    =  1.0
   zetac_sh    =  zetac_sm
   !----- Define chim and chih so the functions are continuous. ---------------------------!
   chim        = (-zetac_um) ** onesixth / sqrt(sqrt(1.0 - gamm * zetac_um))
   chih        = cbrt(-zetac_uh) / sqrt(1.0 - gamh * zetac_uh)
   beta_vs     = 1.0 - (1.0 - beta_s) * zetac_sm
   !----- Define derived values to speed up the code a little. ----------------------------!
   zetac_umi   = 1.0 / zetac_um
   zetac_uhi   = 1.0 / zetac_uh
   zetac_smi   = 1.0 / zetac_sm
   zetac_shi   = 1.0 / zetac_sh
   zetac_umi16 = 1.0 / (- zetac_um) ** onesixth
   zetac_uhi13 = 1.0 / cbrt(-zetac_uh)


   !----- Parameters for the z0m:z0h ratio, following Zeng and Dickinson (1998). ----------!
   zd98_a     = 0.13 * tprandtl
   zd98_b     = 0.45
   zd98_emax  = 10.
   zd98_a8    = dble(zd98_a   )
   zd98_b8    = dble(zd98_b   )
   zd98_emax8 = dble(zd98_emax)
   !---------------------------------------------------------------------------------------!


   
   !----- Legacy variable, we can probably remove it. -------------------------------------!
   ez  = 0.172
   !---------------------------------------------------------------------------------------!







   !---------------------------------------------------------------------------------------!
   !      Parameters for the aerodynamic resistance between the leaf (flat surface) and    !
   ! wood (kind of cylinder surface), and the canopy air space.  These are the A, B, n,    !
   ! and m parameters that define the Nusselt number for forced and free convection, at    !
   ! equations 10.7 and 10.9.  The parameters are found at the appendix A, table A.5(a)    !
   ! and A.5(b).                                                                           !
   !                                                                                       !
   ! M08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,     !
   !       3rd. edition, Academic Press, Amsterdam, 418pp.  (Mostly Chapter 10).           !
   !                                                                                       !
   ! The coefficient B for flat plates under turbulent flow was changed to 0.19 so the     !
   !     transition from laminar to turbulent regime will happen at Gr ~ 100,000, the      !
   !     number suggested by M08.                                                          !
   !---------------------------------------------------------------------------------------!
   aflat_lami = 0.600    ! A (forced convection), laminar   flow
   nflat_lami = 0.500    ! n (forced convection), laminar   flow
   aflat_turb = 0.032    ! A (forced convection), turbulent flow
   nflat_turb = 0.800    ! n (forced convection), turbulent flow
   bflat_lami = 0.500    ! B (free   convection), laminar   flow
   mflat_lami = 0.250    ! m (free   convection), laminar   flow
   bflat_turb = 0.190    ! B (free   convection), turbulent flow
   mflat_turb = onethird ! m (free   convection), turbulent flow
   ocyli_lami = 0.320    ! intercept (forced convection), laminar   flow
   acyli_lami = 0.510    ! A (forced convection), laminar   flow
   ncyli_lami = 0.520    ! n (forced convection), laminar   flow
   ocyli_turb = 0.000    ! intercept (forced convection), turbulent flow
   acyli_turb = 0.240    ! A (forced convection), turbulent flow
   ncyli_turb = 0.600    ! n (forced convection), turbulent flow
   bcyli_lami = 0.480    ! B (free   convection), laminar   flow
   mcyli_lami = 0.250    ! m (free   convection), laminar   flow
   bcyli_turb = 0.090    ! B (free   convection), turbulent flow
   mcyli_turb = onethird ! m (free   convection), turbulent flow
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Define the variables that are going to be used by Massman (1997) and Massman and  !
   ! Weil (1999).  Full reference:                                                         !
   !                                                                                       !
   ! Massman, W. J., 1997: An analytical one-dimensional model of momentum transfer by     !
   !    vegetation of arbitrary structure.  Boundary-Layer Meteorol., 83, 407-421.         !
   !                                                                                       !
   ! Massman, W. J., and J. C. Weil, 1999: An analytical one-dimension second-order clos-  !
   !    ure model turbulence statistics and the Lagrangian time scale within and above     !
   !    plant canopies of arbitrary structure.  Boundary-Layer Meteorol., 91, 81-107.      !
   !                                                                                       !
   ! Wohlfahrt, G., and A. Cernusca, 2002: Momentum transfer by a mountain meadow canopy:  !
   !    a simulation analysis based on Massman's (1997) model.  Boundary-Layer Meteorol.,  !
   !    103, 391-407.
   !---------------------------------------------------------------------------------------!
   !----- Fluid drag coefficient for turbulent flow in leaves. ----------------------------!
   cdrag0    = 0.2
   !----- Values from re-fit of the data used by WC02. ------------------------------------!
   cdrag1    = 0.086
   cdrag2    = 1.192
   cdrag3    = 0.480
   !----- Sheltering factor of fluid drag on canopies. ------------------------------------!
   pm0       = 1.0
   !----- Surface drag parameters (Massman 1997). -----------------------------------------!
   c1_m97    = 0.320
   c2_m97    = 0.264
   c3_m97    = 15.1
   !----- Eddy diffusivity due to Von Karman Wakes in gravity flows. ----------------------!
   kvwake    = 0.001
   !---------------------------------------------------------------------------------------!
   !     Alpha factors to produce the profile of sheltering factor and within canopy drag, !
   ! as suggested by Massman (1997) and Massman and Weil (1999).                           !
   !---------------------------------------------------------------------------------------!
   alpha_m97  = 5.00
   alpha_mw99 = 0.03
   !---------------------------------------------------------------------------------------!
   !      Parameter to represent the roughness sublayer effect.  According to Massman,     !
   ! assuming this to be zero means that the sublayer effects will be ignored.  Otherwise  !
   ! Raupach (1994) tried values up to 0.316.                                              !
   !---------------------------------------------------------------------------------------!
   infunc    = 0.193
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters for CLM, at equation 5.103 of CLM-4 techical note.                     !
   !     Oleson, K. W., et al.; Technical description of version 4.0 of the community land !
   !        model (CLM) NCAR Technical Note NCAR/TN-478+STR, Boulder, CO, April 2010.      !
   !---------------------------------------------------------------------------------------!
   cs_dense0  = 0.004
   gamma_clm4 = 0.5
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !  Gamma and nu are the parameters that close equation 10 in Massman and Weil (1999).   !
   !---------------------------------------------------------------------------------------!
   gamma_mw99 = (/2.4, 1.9, 1.25/)
   nu_mw99(1) = 1.0 / sqrt(sum(gamma_mw99(:)*gamma_mw99(:)))
   nu_mw99(2) = ( 1.0 - 3. * (nu_mw99(1)*gamma_mw99(3)) * (nu_mw99(1)*gamma_mw99(3)) )     &
              / ( 6.0 * nu_mw99(1) * nu_mw99(1) * nu_mw99(1) )
   nu_mw99(3) = 1.0 / ( nu_mw99(1) * nu_mw99(1) * nu_mw99(1) )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   Soil conductance terms, from:                                                       !
   !                                                                                       !
   ! Passerat de Silans, A., 1986: Transferts de masse et de chaleur dans un sol stratifie !
   !     soumis a une excitation amtospherique naturelle. Comparaison: Modeles-experience. !
   !     Thesis, Institut National Polytechnique de Grenoble. (P86)                        !
   !                                                                                       !
   ! retrieved from:                                                                       !
   ! Mahfouf, J. F., J. Noilhan, 1991: Comparative study of various formulations of        !
   !     evaporation from bare soil using in situ data. J. Appl. Meteorol., 30, 1354-1365. !
   !     (MN91)                                                                            !
   !                                                                                       !
   !     Please notice that the values are inverted because we compute conductance, not    !
   ! resistance.                                                                           !
   !---------------------------------------------------------------------------------------!
   ggsoil0 = 1. / 38113.
   kksoil  = 13.515
   !---------------------------------------------------------------------------------------!


   !-----  Multiplication factor for initial leaf/wood boundary layer conductance. --------!
   f_bndlyr_init = 10.0
   !---------------------------------------------------------------------------------------!



   !----- Set the double precision variables. ---------------------------------------------!
   veg_height_min8       = dble(veg_height_min      )
   minimum_canopy_depth8 = dble(minimum_canopy_depth)
   exar8                 = dble(exar                )
   ubmin8                = dble(ubmin               )
   ugbmin8               = dble(ugbmin              )
   ustmin8               = dble(ustmin              )
   ez8                   = dble(ez                  )
   vh2vr8                = dble(vh2vr               )
   vh2dh8                = dble(vh2dh               )
   bl798                 = dble(bl79                )
   csm8                  = dble(csm                 )
   csh8                  = dble(csh                 )
   dl798                 = dble(dl79                )
   beta_s8               = dble(beta_s              )
   gamm8                 = dble(gamm                )
   gamh8                 = dble(gamh                )
   ribmax8               = dble(ribmax              )
   tprandtl8             = dble(tprandtl            )
   abh918                = dble(abh91               )
   bbh918                = dble(bbh91               )
   cbh918                = dble(cbh91               )
   dbh918                = dble(dbh91               )
   ebh918                = dble(ebh91               )
   fbh918                = dble(fbh91               )
   cod8                  = dble(cod                 )
   bcod8                 = dble(bcod                )
   fm18                  = dble(fm1                 )
   ate8                  = dble(ate                 )
   atetf8                = dble(atetf               )
   aflat_lami8           = dble(aflat_lami          )
   nflat_lami8           = dble(nflat_lami          )
   aflat_turb8           = dble(aflat_turb          )
   nflat_turb8           = dble(nflat_turb          )
   bflat_lami8           = dble(bflat_lami          )
   mflat_lami8           = dble(mflat_lami          )
   bflat_turb8           = dble(bflat_turb          )
   mflat_turb8           = dble(mflat_turb          )
   ocyli_lami8           = dble(ocyli_lami          )
   acyli_lami8           = dble(acyli_lami          )
   ncyli_lami8           = dble(ncyli_lami          )
   ocyli_turb8           = dble(ocyli_turb          )
   acyli_turb8           = dble(acyli_turb          )
   ncyli_turb8           = dble(ncyli_turb          )
   bcyli_lami8           = dble(bcyli_lami          )
   mcyli_lami8           = dble(mcyli_lami          )
   bcyli_turb8           = dble(bcyli_turb          )
   mcyli_turb8           = dble(mcyli_turb          )
   cdrag08               = dble(cdrag0              )
   cdrag18               = dble(cdrag1              )
   cdrag28               = dble(cdrag2              )
   cdrag38               = dble(cdrag3              )
   pm08                  = dble(pm0                 )
   c1_m978               = dble(c1_m97              )
   c2_m978               = dble(c2_m97              )
   c3_m978               = dble(c3_m97              )
   kvwake8               = dble(kvwake              )
   alpha_m97_8           = dble(alpha_m97           )
   alpha_mw99_8          = dble(alpha_mw99          )
   gamma_mw99_8          = dble(gamma_mw99          )
   nu_mw99_8             = dble(nu_mw99             )
   infunc_8              = dble(infunc              )
   cs_dense08            = dble(cs_dense0           )
   ggsoil08              = dble(ggsoil0             )
   kksoil8               = dble(kksoil              )
   zetac_um8             = dble(zetac_um            )
   zetac_uh8             = dble(zetac_uh            )
   zetac_sm8             = dble(zetac_sm            )
   zetac_sh8             = dble(zetac_sh            )
   chim8                 = dble(chim                )
   chih8                 = dble(chih                )
   beta_vs8              = dble(beta_vs             )
   zetac_umi8            = dble(zetac_umi           )
   zetac_uhi8            = dble(zetac_uhi           )
   zetac_smi8            = dble(zetac_smi           )
   zetac_shi8            = dble(zetac_shi           )
   zetac_umi168          = dble(zetac_umi16         )
   zetac_uhi138          = dble(zetac_uhi13         )
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Initialise these values with dummies, it will be updated after we define the      !
   ! functions.                                                                            !
   !---------------------------------------------------------------------------------------!
   psimc_um8  = psim8(zetac_um8,.false.)
   psimc_um   = sngloff(psimc_um8,tiny_offset)
   psihc_uh8  = psih8(zetac_uh8,.false.)
   psihc_uh   = sngloff(psihc_uh8,tiny_offset)
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_can_air_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will assign some radiation related parameters.                        !
!------------------------------------------------------------------------------------------!
subroutine init_can_rad_params()

   use canopy_radiation_coms , only : ltrans_vis                  & ! intent(in)
                                    , ltrans_nir                  & ! intent(in)
                                    , lreflect_vis                & ! intent(in)
                                    , lreflect_nir                & ! intent(in)
                                    , orient_tree                 & ! intent(in)
                                    , orient_grass                & ! intent(in)
                                    , clump_tree                  & ! intent(in)
                                    , clump_grass                 & ! intent(in)
                                    , leaf_reflect_nir            & ! intent(out)
                                    , leaf_trans_nir              & ! intent(out)
                                    , leaf_reflect_vis            & ! intent(out)
                                    , leaf_trans_vis              & ! intent(out)
                                    , leaf_emiss_tir              & ! intent(out)
                                    , clumping_factor             & ! intent(out)
                                    , orient_factor               & ! intent(out)
                                    , wood_reflect_nir            & ! intent(out)
                                    , wood_trans_nir              & ! intent(out)
                                    , wood_reflect_vis            & ! intent(out)
                                    , wood_trans_vis              & ! intent(out)
                                    , wood_emiss_tir              & ! intent(out)
                                    , fvis_beam_def               & ! intent(out)
                                    , fvis_diff_def               & ! intent(out)
                                    , fnir_beam_def               & ! intent(out)
                                    , fnir_diff_def               & ! intent(out)
                                    , snow_albedo_vis             & ! intent(out)
                                    , snow_albedo_nir             & ! intent(out)
                                    , snow_emiss_tir              & ! intent(out)
                                    , rshort_twilight_min         & ! intent(out)
                                    , cosz_min                    & ! intent(out)
                                    , cosz_min8                   ! ! intent(out)
   use pft_coms              , only : is_grass                    & ! intent(in)
                                    , is_tropical                 & ! intent(in)
                                    , is_conifer                  ! ! intent(in)
   use consts_coms           , only : pio180                      & ! intent(in)
                                    , twothirds8                  ! ! intent(in)
   use ed_max_dims           , only : n_pft                       ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer :: ipft
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      The following parameters are used to split the shortwave radiation into visible  !
   ! and near-infrared radiation.                                                          !
   !---------------------------------------------------------------------------------------!
   fvis_beam_def = 0.43
   fnir_beam_def = 1.0 - fvis_beam_def
   fvis_diff_def = 0.52
   fnir_diff_def = 1.0 - fvis_diff_def
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Clumping factor.  This factor indicates the degree of clumpiness of leaves.       !a
   !  0 -- black hole                                                                      !
   !  1 -- homogeneous, no clumping.                                                       !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      if (is_conifer(ipft)) then ! Conifers (subtropical and temperate). 
         clumping_factor(ipft) = 7.350d-1
         !---------------------------------------------------------------------------------!
      elseif (is_tropical(ipft) .and. is_grass(ipft)) then ! Tropical grasses. 
         clumping_factor(ipft) = dble(clump_grass)
      elseif (is_tropical(ipft)) then ! Lianas and tropical trees. 
         clumping_factor(ipft) = dble(clump_tree)
      else ! Temperate broadleaf (trees pr grasses). 
         clumping_factor(ipft) = 8.400d-1
      end if
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Orientation factor.  The numbers come from CLM, and the original value from      !
   ! ED-2.1 used to 0.  This works in the following way:                                   !
   !  0 -- leaves are randomly oriented                                                    !
   !  1 -- all leaves are perfectly horizontal                                             !
   ! -1 -- all leaves are perfectly vertical.                                              !
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      if (.not. is_tropical(ipft)) then ! Temperate PFTs (grasses and trees). 
         orient_factor(ipft) = 0.d0
      elseif (is_conifer(ipft)) then ! Araucaria, (CLM value for evergreen needleleaf).
         orient_factor(ipft) = 1.0d-2
      elseif (is_grass(ipft)) then ! Tropical grasses. 
         orient_factor(ipft) = dble(orient_grass)
      else ! Lianas and tropical broadleaf trees. 
         orient_factor(ipft) = dble(orient_tree)
      end if
   end do
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !      Emissivity on Thermal infra-red (TIR).                                           !
   !---------------------------------------------------------------------------------------!
   leaf_emiss_tir(:) = merge(9.60d-1,merge(9.70d-1,9.50d-1,is_conifer(:)),is_grass(:))
   wood_emiss_tir(:) = merge(9.60d-1,9.00d-1,is_grass(:))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Leaf reflectance.                                                                !
   !      Values for temperate PFTs were left as they were.  Tropical and sub-tropical     !
   ! PFTs use the parameters from CLM.  I checked the values against some published and    !
   ! they seem similar at a first glance, at least closer than the original values, which  !
   ! looked like the visible ignoring the green band.                                      !
   !                                                                                       !
   ! Tropical / Subtropical values for visible came from:                                  !
   ! - Poorter, L., S. F. Oberbauer, D. B. Clark, 1995: Leaf optical properties along a    !
   !      vertical gradient in a tropical rainforest in Costa Rica. American J. of Botany, !
   !      82, 1257-1263.                                                                   !
   ! Tropical values for NIR were estimated from:                                          !
   ! - Roberts, D. A., B. W. Nelson, J. B. Adams, F. Palmer, 1998: Spectral changes with   !
   !      leaf aging in Amazon caatinga. Trees, 12, 315-325.                               !
   !---------------------------------------------------------------------------------------!
   leaf_reflect_vis(:) = merge( merge(9.00d-2,dble(lreflect_vis),is_conifer(:))            &
                              , 1.10d-1                                                    &
                              , is_tropical(:) )
   leaf_reflect_nir(:) = merge( dble(lreflect_nir)                                         &
                              , 5.77d-1                                                    &
                              , is_tropical(:) .and. (.not. is_conifer(:)) )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Wood reflectance, using values based on:                                         !
   !                                                                                       !
   ! Asner, G., 1998: Biophysical and biochemical sources of variability in canopy         !
   !     reflectance. Remote Sensing of Environment, 64, 234-253.                          !
   !                                                                                       !
   ! Commented values are from CLM, but they were quite high.                              !
   !---------------------------------------------------------------------------------------!
   wood_reflect_vis(:) = merge(1.60d-1,1.10d-1,is_grass(:))
   wood_reflect_nir(:) = 2.50d-1
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Leaf transmittance.                                                              !
   !      Values for temperate PFTs were left as they were.  Tropical and sub-tropical     !
   ! PFTs use the parameters from CLM.  I checked the values against some published and    !
   ! they seem similar at a first glance, at least closer than the original values, which  !
   ! looked like the visible ignoring the green band.                                      !
   !                                                                                       !
   ! Tropical / Subtropical values for visible came from:                                  !
   ! - Poorter, L., S. F. Oberbauer, D. B. Clark, 1995: Leaf optical properties along a    !
   !      vertical gradient in a tropical rainforest in Costa Rica. American J. of Botany, !
   !      82, 1257-1263.                                                                   !
   ! Tropical values for NIR were estimated from:                                          !
   ! - Roberts, D. A., B. W. Nelson, J. B. Adams, F. Palmer, 1998: Spectral changes with   !
   !      leaf aging in Amazon caatinga. Trees, 12, 315-325.                               !
   !---------------------------------------------------------------------------------------!
   leaf_trans_vis(:) = merge( merge(5.00d-2,dble(ltrans_vis),is_conifer(:))                &
                            , 1.60d-1                                                      &
                            , is_tropical(:) )
   leaf_trans_nir(:) = merge( dble(ltrans_nir)                                             &
                            , 2.48d-1                                                      &
                            , is_tropical(:) .and. (.not. is_conifer(:)) )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Wood transmittance, using the parameters from CLM.                               !
   !---------------------------------------------------------------------------------------!
   wood_trans_vis(:) = merge(2.80d-2,1.00d-3,is_grass(:))
   wood_trans_nir(:) = merge(2.48d-1,1.00d-3,is_grass(:))
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !     Optical properties for snow.  Values are a first guess, and a more thorough snow  !
   ! model that takes snow age and snow melt into account (like in CLM-4 or ECHAM-5) are   !
   ! very welcome.                                                                         !
   !                                                                                       !
   !  References for current snow values:                                                  !
   !  Roesch, A., et al., 2002: Comparison of spectral surface albedos and their           !
   !      impact on the general circulation model simulated surface climate.  J.           !
   !      Geophys. Res.-Atmosph., 107(D14), 4221, 10.1029/2001JD000809.                    !
   !      Average between minimum and maximum snow albedo on land, af = 0. and af=1.       !
   !                                                                                       !
   !  Oleson, K.W., et al., 2010: Technical description of version 4.0 of the              !
   !      Community Land Model (CLM). NCAR Technical Note NCAR/TN-478+STR.                 !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   snow_albedo_vis = 0.518
   snow_albedo_nir = 0.435
   snow_emiss_tir  = 0.970
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     These variables are the thresholds for things that should be computed during the  !
   ! day time hours only.                                                                  !
   !---------------------------------------------------------------------------------------!
   rshort_twilight_min = 0.5
   cosz_min            = cos(89.*pio180) !cos(89.5*pio180)
   cosz_min8           = dble(cosz_min)
   !---------------------------------------------------------------------------------------!


   return
end subroutine init_can_rad_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine initialises some of the canopy layer variables.  These are used     !
! by sub-routines that need to calculate the canopy properties by layers rather than by    !
! cohorts (or when both must be considered).                                               !
!------------------------------------------------------------------------------------------!
subroutine init_can_lyr_params()
   use canopy_layer_coms, only : tai_lyr_max                 & ! intent(out)
                               , ncanlyr                     & ! intent(out)
                               , ncanlyrp1                   & ! intent(out)
                               , ncanlyrt2                   & ! intent(out)
                               , zztop0                      & ! intent(out)
                               , zztop08                     & ! intent(out)
                               , zztop0i                     & ! intent(out)
                               , zztop0i8                    & ! intent(out)
                               , ehgt                        & ! intent(out)
                               , ehgt8                       & ! intent(out)
                               , ehgti                       & ! intent(out)
                               , ehgti8                      & ! intent(out)
                               , dzcan                       & ! intent(out)
                               , dzcan8                      & ! intent(out)
                               , zztop                       & ! intent(out)
                               , zzmid                       & ! intent(out)
                               , zzbot                       & ! intent(out)
                               , zztop8                      & ! intent(out)
                               , zzmid8                      & ! intent(out)
                               , zzbot8                      & ! intent(out)
                               , alloc_canopy_layer          ! ! subroutine
   use pft_coms         , only : hgt_min                     & ! intent(in)
                               , hgt_max                     ! ! intent(in)
   use consts_coms      , only : onethird                    & ! intent(in)
                               , onethird8                   ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer    :: ilyr
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Set the maximum tai that each layer is allowed to have.                          !
   !---------------------------------------------------------------------------------------!
   tai_lyr_max = 1.0
   !---------------------------------------------------------------------------------------!



   !----- Find the layer thickness and the number of layers needed. -----------------------!
   ncanlyr   = 100
   ncanlyrp1 = ncanlyr + 1
   ncanlyrt2 = ncanlyr * 2
   zztop0    = onethird  * minval(hgt_min)
   zztop08   = onethird8 * dble(minval(hgt_min))
   zztop0i   = 1.   / zztop0
   zztop0i8  = 1.d0 / zztop08
   ehgt      = log(maxval(hgt_max)/zztop0)        / log(real(ncanlyr))
   ehgt8     = log(dble(maxval(hgt_max))/zztop08) / log(dble(ncanlyr))
   ehgti     = 1./ ehgt
   ehgti8    = 1.d0 / ehgt8
   !---------------------------------------------------------------------------------------!



   !----- Allocate the variables. ---------------------------------------------------------!
   call alloc_canopy_layer()
   !---------------------------------------------------------------------------------------!



   !----- Define the layer heights. -------------------------------------------------------!
   do ilyr =1,ncanlyr
      zztop (ilyr) = zztop0 * real(ilyr  ) ** ehgt
      zzbot (ilyr) = zztop0 * real(ilyr-1) ** ehgt
      dzcan (ilyr) = zztop(ilyr) - zzbot(ilyr)
      zzmid (ilyr) = 0.5 * (zzbot(ilyr) + zztop(ilyr))
      zztop8(ilyr) = zztop0  * dble(ilyr  ) ** ehgt8
      zzbot8(ilyr) = zztop08 * dble(ilyr-1) ** ehgt8
      dzcan8(ilyr) = zztop8(ilyr) - zzbot8(ilyr)
      zzmid8(ilyr) = 5.d-1 * (zzbot8(ilyr) + zztop8(ilyr))
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_can_lyr_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine assigns various parameters for the thermodynamics solver (Euler,      !
! Heun, Runge-Kutta or Hybrid).  It uses many values previously assigned in other          !
! parameter initialisation, so this should be the last one called.                         !
!------------------------------------------------------------------------------------------!
subroutine init_dt_thermo_params()
   use soil_coms      , only : water_stab_thresh      & ! intent(in)
                             , snowmin                & ! intent(in)
                             , tiny_sfcwater_mass     ! ! intent(in)
   use canopy_air_coms, only : leaf_drywhc            & ! intent(in)
                             , leaf_maxwhc            ! ! intent(in)
   use ed_misc_coms   , only : dtlsm                  & ! intent(in)
                             , ffilout                & ! intent(in)
                             , nsub_euler             & ! intent(in)
                             , integration_scheme     & ! intent(in)
                             , dteuler                ! ! intent(out)
   use consts_coms    , only : wdnsi8                 & ! intent(in)
                             , r_tol_trunc            ! ! intent(in)
   use detailed_coms  , only : idetailed              ! ! intent(in)
   use pft_coms       , only : is_conifer             ! ! intent(in)
   use budget_utils   , only : tol_subday_budget      & ! intent(in)
                             , tol_carbon_budget      ! ! intent(in)
   use rk4_coms       , only : rk4_tolerance          & ! intent(in)
                             , maxstp                 & ! intent(out)
                             , rk4eps                 & ! intent(out)
                             , rk4eps2                & ! intent(out)
                             , rk4epsi                & ! intent(out)
                             , hmin                   & ! intent(out)
                             , print_diags            & ! intent(out)
                             , checkbudget            & ! intent(out)
                             , debug                  & ! intent(out)
                             , toocold                & ! intent(out)
                             , toohot                 & ! intent(out)
                             , lai_to_cover           & ! intent(out)
                             , rk4min_veg_temp        & ! intent(out)
                             , rk4water_stab_thresh   & ! intent(out)
                             , rk4tiny_sfcw_mass      & ! intent(out)
                             , rk4tiny_sfcw_depth     & ! intent(out)
                             , rk4leaf_drywhc         & ! intent(out)
                             , rk4leaf_maxwhc         & ! intent(out)
                             , rk4snowmin             & ! intent(out)
                             , rk4min_can_temp        & ! intent(out)
                             , rk4max_can_temp        & ! intent(out)
                             , rk4min_can_shv         & ! intent(out)
                             , rk4max_can_shv         & ! intent(out)
                             , rk4max_can_rhv         & ! intent(out)
                             , rk4min_can_co2         & ! intent(out)
                             , rk4max_can_co2         & ! intent(out)
                             , rk4min_soil_temp       & ! intent(out)
                             , rk4max_soil_temp       & ! intent(out)
                             , rk4min_veg_temp        & ! intent(out)
                             , rk4max_veg_temp        & ! intent(out)
                             , rk4min_veg_lwater      & ! intent(out)
                             , rk4min_sfcw_temp       & ! intent(out)
                             , rk4max_sfcw_temp       & ! intent(out)
                             , rk4min_sfcw_moist      & ! intent(out)
                             , rk4min_virt_moist      & ! intent(out)
                             , effarea_heat           & ! intent(out)
                             , effarea_evap           & ! intent(out)
                             , effarea_transp         & ! intent(out)
                             , leaf_intercept         & ! intent(out)
                             , supersat_ok            & ! intent(out)
                             , record_err             & ! intent(out)
                             , print_detailed         & ! intent(out)
                             , print_budget           & ! intent(out)
                             , print_thbnd            & ! intent(out)
                             , errmax_fout            & ! intent(out)
                             , sanity_fout            & ! intent(out)
                             , thbnds_fout            & ! intent(out)
                             , detail_pref            & ! intent(out)
                             , budget_pref            ! ! intent(out)
   use ed_misc_coms   , only : fast_diagnostics       ! ! intent(inout)
   implicit none

   !---------------------------------------------------------------------------------------!
   !     Define the maximum time step for forward Euler solver.  Because forward Euler is  !
   ! lower order, the time step should be typically less than dtlsm.                       !
   !---------------------------------------------------------------------------------------!
   dteuler = dtlsm / real(nsub_euler)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Copy some variables to the Runge-Kutta counterpart (double precision).            !
   !---------------------------------------------------------------------------------------!
   rk4water_stab_thresh  = dble(water_stab_thresh )
   rk4tiny_sfcw_mass     = dble(tiny_sfcwater_mass)
   rk4leaf_drywhc        = dble(leaf_drywhc       )
   rk4leaf_maxwhc        = dble(leaf_maxwhc       )
   rk4snowmin            = dble(snowmin           )
   rk4tiny_sfcw_depth    = rk4tiny_sfcw_mass  * wdnsi8
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    The following variables control the performance of the Runge-Kutta and Euler       !
   ! integration schemes. Think twice before changing them...                              !
   !---------------------------------------------------------------------------------------!
   maxstp      = 100000000           ! Maximum number of intermediate steps.
   rk4eps      = dble(rk4_tolerance) ! The desired accuracy.
   rk4epsi     = 1.d0/rk4eps         ! The inverse of desired accuracy.
   rk4eps2     = rk4eps**2           ! square of the accuracy
   hmin        = 1.d-7               ! The minimum step size.
   print_diags = .false.             ! Flag to print the diagnostic check.
   checkbudget = .true.              ! Flag to check CO2, water, and energy budgets every
                                     !     time step and stop the run in case any of these
                                     !     budgets don't close.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Tolerances.  Following Stefan Olin's suggestion on the ED-2.2 model description  !
   ! paper, we use a stricter tolerance, by default the truncation tolerance (about 1e-5). !
   ! For carbon, we use 10 times the the values for the energy and water because the       !
   ! solver uses single-precision.  In case we use hybdrid, we relax tolerance for the     !
   ! time being.  We should identify the causes of leakage in that scheme in the future.   !
   !                                                                                       !
   ! Update: There are a few cases in which the strict tolerance is not working.  Most of  !
   ! them seem to be associated with the long-term vegetation dynamics, but there are a    !
   ! few cases that the thermodynamics is also causing crashes (I think still related to   !
   ! temporary surface water).  This still needs to be addressed.  For the time being, I   !
   ! am relaxing the tolerance so people can run their simulations.                        !
   !---------------------------------------------------------------------------------------!
   select case (integration_scheme)
   case (3)
      tol_subday_budget = 100. * r_tol_trunc
      tol_carbon_budget = 100. * r_tol_trunc
   case default
      tol_subday_budget = 10   * r_tol_trunc
      tol_carbon_budget = 100. * r_tol_trunc
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Miscellaneous constants used in rk4_derivs.                                       !
   !---------------------------------------------------------------------------------------!
   debug         = .false.  ! Verbose output for debug                             [   T|F]
   toocold       = 1.5315d2 ! Minimum temperature for saturation specific hum.     [     K]
   toohot        = 3.6315d2 ! Maximum temperature for saturation specific hum.     [     K]
   lai_to_cover  = 1.5d0    ! Canopies with LAI less than this number  are assumed to be
                            !     open, ie, some fraction of the rain-drops can reach
                            !    the soil/litter layer unimpeded.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Variables used to keep track on the error.  We use the idetailed flag to          !
   ! determine whether to create the output value or not.                                  !
   !---------------------------------------------------------------------------------------!
   !------ Detailed budget (every DTLSM). -------------------------------------------------!
   print_budget   = btest(idetailed,0)
   if (print_budget) checkbudget = .true.
   !------ Detailed output from the integrator (every HDID). ------------------------------!
   print_detailed = btest(idetailed,2)
   !------ Thermodynamic boundaries for sanity check (every HDID). ------------------------!
   print_thbnd    = btest(idetailed,3)
   !------ Daily error statistics (count how often a variable shrunk the time step). ------!
   record_err     = btest(idetailed,4)
   !---------------------------------------------------------------------------------------!
   errmax_fout    = 'error_max_count.txt'    ! File with the maximum error count
   sanity_fout    = 'sanity_check_count.txt' ! File with the sanity check count
   thbnds_fout    = 'thermo_bounds.txt'      ! File with the thermodynamic boundaries.
   detail_pref    = 'thermo_state_'          ! Prefix for the detailed thermodynamic file
   budget_pref    = 'budget_state_'          ! File with the thermodynamic boundaries.
   !---------------------------------------------------------------------------------------!



   !----- Append the same prefix used for analysis files. ---------------------------------!
   detail_pref = trim(ffilout)//'_'//trim(detail_pref)
   budget_pref = trim(ffilout)//'_'//trim(budget_pref)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Assigning some default values for the bounds at the sanity check.  Units are      !
   ! usually the standard, but a few of them are defined differently so they can be scaled !
   ! depending on the cohort and soil grid definitions.                                    !
   !---------------------------------------------------------------------------------------!
   rk4min_can_temp   =  1.7400d2  ! Minimum canopy    temperature               [        K]
   rk4max_can_temp   =  3.6000d2  ! Maximum canopy    temperature               [        K]
   rk4min_can_shv    =  1.0000d-8 ! Minimum canopy    specific humidity         [kg/kg_air]
   rk4max_can_shv    =  8.0000d-2 ! Maximum canopy    specific humidity         [kg/kg_air]
   rk4max_can_rhv    =  1.1000d0  ! Maximum canopy    relative humidity (**)    [      ---]
   rk4min_can_co2    =  1.0000d0  ! Minimum canopy    CO2 mixing ratio          [ umol/mol]
   rk4max_can_co2    =  5.0000d4  ! Maximum canopy    CO2 mixing ratio          [ umol/mol]
   rk4min_soil_temp  =  1.7400d2  ! Minimum soil      temperature               [        K]
   rk4max_soil_temp  =  3.6000d2  ! Maximum soil      temperature               [        K]
   rk4min_veg_temp   =  1.7400d2  ! Minimum leaf      temperature               [        K]
   rk4max_veg_temp   =  3.6000d2  ! Maximum leaf      temperature               [        K]
   rk4min_sfcw_temp  =  1.7400d2  ! Minimum snow/pond temperature               [        K]
   rk4max_sfcw_temp  =  3.6000d2  ! Maximum snow/pond temperature               [        K]
   !.......................................................................................!
   ! (**) Please, don't be too strict here.  The model currently doesn't have radiation    !
   !      fog, so supersaturation may happen.  This is a problem we may want to address in !
   !      the future, though...                                                            !
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Minimum water mass at the leaf surface.  This is given in kg/m2leaf rather than   !
   ! kg/m2ground, so we scale it with LAI.                                                 !
   !---------------------------------------------------------------------------------------!
   rk4min_veg_lwater = -rk4leaf_drywhc            ! Minimum leaf water mass     [kg/m2leaf]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    The minimum mass of surface water and virtual layer are given in m3/m3 rather than !
   ! kg/m2.  This is because there will be exchange between the top soil layer and the     !
   ! layers above in case the mass goes below the minimum.  Since this would make the im-  !
   ! pact of such exchange dependent on the soil depth, we assign the scale a function of  !
   ! the top layer thickness.                                                              !
   !---------------------------------------------------------------------------------------!
   rk4min_sfcw_moist = -5.0000d-4                  ! Minimum water mass allowed.
   rk4min_virt_moist = -5.0000d-4                  ! Minimum water allowed at virtual pool.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     These variables are assigned in ed_params.f90.  Heat area should be 2.0 for all   !
   !                                                                                       !
   ! PFTs (two sides of the leaves exchange heat), and the evaporation area should be 1.0  !
   ! for all PFTs (only one side of the leaf is usually covered by water).  Sometimes heat !
   ! and evaporation are multiplied  by 1.2 and 2.2 to account for branches and twigs.     !
   !  This is not recommended, though, because branches and twigs do not contribute to     !
   ! heat storage when ibranch_thermo is set to zero, and they are otherwise accounted     !
   ! through the wood area index.                                                          !
   !---------------------------------------------------------------------------------------!
   effarea_heat   = 2.d0 ! Heat area: related to 2*LAI
   effarea_evap   = 1.d0 ! Evaporation area: related to LAI
   !---------------------------------------------------------------------------------------!
   !     Transpiration.  Possible values are.                                              !
   !                                                                                       !
   ! 1.d0 - Hypostomatous.                                                                 !
   ! 2.d0 - Symmetrical or Amphistomatous. (currently this means conifers)                 !
   !---------------------------------------------------------------------------------------!
   effarea_transp(:) = merge(2.d0,1.d0,is_conifer(:))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    This flag is used to control evaporation and transpiration when the air is         !
   ! saturated or super-saturated.  If supersat_ok is TRUE, then evaporation and           !
   ! transpiration will continue to happen even if the air is super-saturated at the       !
   ! canopy air temperature (but not at the soil or vegetation temperature).  Otherwise,   !
   ! evaporation and transpiration will be interrupted until the air becomes sub-saturated !
   ! again.  The air can still become super-saturated because mixing with the free atmo-   !
   ! sphere will not stop, but that is unlikely.                                           !
   !---------------------------------------------------------------------------------------!
   supersat_ok = .false.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This flag is to turn on and on the leaf interception.  Except for developer       !
   ! tests, this variable should be always true.                                           !
   !---------------------------------------------------------------------------------------!
   leaf_intercept = .true.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Update fast_diagnostics in case checkbudget is set to true.                      !
   !---------------------------------------------------------------------------------------!
   fast_diagnostics = fast_diagnostics .or. checkbudget
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_dt_thermo_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This sub-routine seeks XML files. In case it exists, it reads values and and over-     !
! writes parameters that are found in the XML file.                                        !
!------------------------------------------------------------------------------------------!
subroutine overwrite_with_xml_config(thisnode)
   use ed_max_dims   , only : n_pft
   use ed_misc_coms  , only : iedcnfgf
   use hydrology_coms, only : useTOPMODEL
   use soil_coms     , only : isoilbc
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: thisnode
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: max_pft_xml
   logical             :: iamhere
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Check whether the user provided XML and whether the file exists.                   !
   !---------------------------------------------------------------------------------------!

   if (iedcnfgf /= '') then
   
      !----- Test whether the file exists. ------------------------------------------------!
      inquire (file=iedcnfgf,exist=iamhere)
      !------------------------------------------------------------------------------------!


      if (iamhere) then
         !----- First, determine number of pft's defined in xml file. ---------------------!
         call count_pft_xml_config(trim(iedcnfgf),max_pft_xml)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Make sure XML does not try to read PFT indices that are greater than n_pft.  !
         !---------------------------------------------------------------------------------!
         if (max_pft_xml > n_pft) then
            write(unit=*,fmt='(a)'   ) '=================================================='
            write(unit=*,fmt='(a)'   ) '=================================================='
            write(unit=*,fmt='(a)'   ) ''
            write(unit=*,fmt='(a,a)' ) ' - XML file    -- ',trim(iedcnfgf)
            write(unit=*,fmt='(a,i6)') ' - N_PFT (ED2) -- ',n_pft
            write(unit=*,fmt='(a,i6)') ' - N_PFT (XML) -- ',max_pft_xml
            write(unit=*,fmt='(a)'   ) ''
            write(unit=*,fmt='(a)'   ) '    Number of PFTs required by the XML configu-'
            write(unit=*,fmt='(a)'   ) ' ration file exceeds the memory available.'
            write(unit=*,fmt='(a)'   ) ' Please change n_pft in module ed_max_dims and'
            write(unit=*,fmt='(a)'   ) ' recompile ED-2.2.'
            write(unit=*,fmt='(a)'   ) '=================================================='
            call fatal_error('Too many PFTs','overwrite_with_xml_config','ed_params.f90')
         end if
         !---------------------------------------------------------------------------------!


         !----- Update parameter defaults from XML. ---------------------------------------!
         call read_ed_xml_config(trim(iedcnfgf))
         !---------------------------------------------------------------------------------!


         !----- Recalculate any derived values based on xml. ------------------------------!
         if(useTOPMODEL == 1) then
            isoilbc = 0
            print*,"TOPMODEL enabled, setting ISOILBC to 0"
         end if
         !---------------------------------------------------------------------------------!


         !----- Write out a copy of the settings. -----------------------------------------!
         !call write_ed_xml_config() always call write_ed_xml_config in ed_driver after
         !init_derived_params_after_xml()
         !---------------------------------------------------------------------------------!
      else if (thisnode == 1) then
         !----- XML file not found, warn the user. ----------------------------------------!
         write(unit=*,fmt='(a)'   ) '=================================================='
         write(unit=*,fmt='(a)'   ) '=================================================='
         write(unit=*,fmt='(a)'   ) '   WARNING! WARNING! WARNING! WARNING! WARNING!   '
         write(unit=*,fmt='(a)'   ) '--------------------------------------------------'
         write(unit=*,fmt='(a)'   ) ''
         write(unit=*,fmt='(a)'   ) '    XML file wasn''t found. Using default'
         write(unit=*,fmt='(a)'   ) ' parameters in ED2.                       '
         write(unit=*,fmt='(a)'   ) ' (You provided '//trim(iedcnfgf)//').'
         write(unit=*,fmt='(a)'   ) ' '
         write(unit=*,fmt='(a)'   ) '=================================================='
         write(unit=*,fmt='(a)'   ) ' '
         write(unit=*,fmt='(a)'   ) ' '
         !call write_ed_xml_config()
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end if  !! end XML
   !---------------------------------------------------------------------------------------!

   ! call write_ed_xml_config()

   return
end subroutine overwrite_with_xml_config
!==========================================================================================!
!==========================================================================================!









!==========================================================================================!
!==========================================================================================!
!    This subroutine initialises parameters that depend on other input parameters.  This   !
! is an incipient sub-routine, many parameters should be moved from previous routines to   !
! here.  Contributions are welcome and very needed.                                        !
!                                                                                          !
!                                                                                          !
! Examples of parameters to include here:                                                  !
! ----------------------------------------                                                 !
!                                                                                          !
! - Scattering coefficients, which depend on reflectivity and transmissivity.              !
! - Specific heat of wood, which depends on water:dry ratio and water-wood bonding.        !
! - Negligible_nplant, that must be linked to the minimum PFT size to avoid surprises.     !
!                                                                                          !
!                                                                                          !
! Examples of parameters that SHOULD NOT be included here:                                 !
! ---------------------------------------------------------                                !
!                                                                                          !
! - Leaf turnover as a function of SLA (Joe's empirical relationship based on forest A     !
!   may be different from Jane's, which was based on forest B).                            !
!                                                                                          !
!                                                                                          !
! Examples of parameters that are initialised here only if not initialised in xml:         !
! ---------------------------------------------------------------------------------        !
!                                                                                          !
! - seed_rain, which must be greater than the minimum PFT size otherwise it will never     !
!   occur.                                                                                 !
! - Rd0, which may depend on dark_respiration_factor (aka gamma), or could be defined      !
!   directly.  Eventually dark_respiration_factor should be phased out and Rd0 should be   !
!   the only variable, but we keep the proportionality factor for legacy.                  !
!                                                                                          !
!                                                                                          !
! IMPORTANT.  Variables flagged as "intent(in)" must appear in XML, whereas variables      !
!    flagged with "intent(out)" must NOT appear in the xml.  Variables flagged as          !
!    "intent(inout)" may appear in the xml but must be initialised with dummy values in    !
!    one of the subroutines above.  Recommended values: undef_real and related variables   !
!    from ed_max_dims, so it is easy to track.                                             !
!------------------------------------------------------------------------------------------!
subroutine init_derived_params_after_xml()
   use detailed_coms        , only : idetailed                 ! ! intent(in)
   use soil_coms            , only : slz                       ! ! intent(in)
   use grid_coms            , only : nzg                       ! ! intent(in)
   use ed_misc_coms         , only : ibigleaf                  & ! intent(in)
                                   , iallom                    & ! intent(in)
                                   , ivegt_dynamics            & ! intent(in)
                                   , radfrq                    & ! intent(in)
                                   , economics_scheme          & ! intent(in)
                                   , lianas_included           ! ! intent(out)
   use ed_max_dims          , only : n_pft                     & ! intent(in)
                                   , str_len                   & ! intent(in)
                                   , undef_real                ! ! intent(in)
   use consts_coms          , only : onesixth                  & ! intent(in)
                                   , twothirds                 & ! intent(in)
                                   , solar                     & ! intent(in)
                                   , t008                      & ! intent(in)
                                   , cliq                      & ! intent(in)
                                   , day_sec                   & ! intent(in)
                                   , yr_sec                    & ! intent(in)
                                   , almost_zero               & ! intent(in)
                                   , lnexp_min8                & ! intent(in)
                                   , lnexp_max8                & ! intent(in)
                                   , tiny_num8                 ! ! intent(in)
   use physiology_coms      , only : iphysiol                  & ! intent(in)
                                   , trait_plasticity_scheme   ! ! intent(in)
   use pft_coms             , only : include_pft               & ! intent(in)
                                   , is_tropical               & ! intent(in)
                                   , is_grass                  & ! intent(in)
                                   , is_savannah               & ! intent(in)
                                   , is_conifer                & ! intent(in)
                                   , is_liana                  & ! intent(in)
                                   , nbt_lut                   & ! intent(in)
                                   , init_density              & ! intent(in)
                                   , init_laimax               & ! intent(in)
                                   , c2n_leaf                  & ! intent(in)
                                   , c2n_stem                  & ! intent(in)
                                   , c2n_storage               & ! intent(in)
                                   , b1Ht                      & ! intent(in)
                                   , b2Ht                      & ! intent(in)
                                   , hgt_min                   & ! intent(in)
                                   , hgt_ref                   & ! intent(in)
                                   , b1Bl                      & ! intent(in)
                                   , b2Bl                      & ! intent(in)
                                   , b1Bs_small                & ! intent(in)
                                   , b2Bs_small                & ! intent(in)
                                   , b1Bs_large                & ! intent(in)
                                   , b2Bs_large                & ! intent(in)
                                   , b1Ca                      & ! intent(in)
                                   , b2Ca                      & ! intent(in)
                                   , b1WAI                     & ! intent(in)
                                   , b2WAI                     & ! intent(in)
                                   , b1SA                      & ! intent(in)
                                   , b2SA                      & ! intent(in)
                                   , b1Xs                      & ! intent(in)
                                   , b1Xb                      & ! intent(in)
                                   , b1Rd                      & ! intent(in)
                                   , b2Rd                      & ! intent(in)
                                   , d18O_ref                  & ! intent(in)
                                   , b1d18O                    & ! intent(in)
                                   , b2d18O                    & ! intent(in)
                                   , b1Efrd                    & ! intent(in)
                                   , b2Efrd                    & ! intent(in)
                                   , min_dbh                   & ! intent(in)
                                   , dbh_bigleaf               & ! intent(in)
                                   , agf_bs                    & ! intent(in)
                                   , q                         & ! intent(in)
                                   , qsw                       & ! intent(in)
                                   , qrhob                     & ! intent(in)
                                   , qbark                     & ! intent(in)
                                   , rho                       & ! intent(in)
                                   , sla                       & ! intent(in)
                                   , sra                       & ! intent(in)
                                   , root_beta                 & ! intent(in)
                                   , pft_name16                & ! intent(in)
                                   , dbh_bigleaf               & ! intent(in)
                                   , f_bstorage_init           & ! intent(in)
                                   , c_grn_leaf_dry            & ! intent(in)
                                   , c_ngrn_wood_dry           & ! intent(in)
                                   , c_ngrn_bark_dry           & ! intent(in)
                                   , leaf_water_sat            & ! intent(in)
                                   , wood_water_sat            & ! intent(in)
                                   , bark_water_sat            & ! intent(in)
                                   , delta_c_wood              & ! intent(in)
                                   , delta_c_bark              & ! intent(in)
                                   , D0                        & ! intent(in)
                                   , Vm_low_temp               & ! intent(in)
                                   , Vm_high_temp              & ! intent(in)
                                   , Vm_decay_elow             & ! intent(in)
                                   , Vm_decay_ehigh            & ! intent(in)
                                   , Vm0                       & ! intent(in)
                                   , Vm_hor                    & ! intent(in)
                                   , Vm_q10                    & ! intent(in)
                                   , stomatal_slope            & ! intent(in)
                                   , leaf_width                & ! intent(in)
                                   , cuticular_cond            & ! intent(in)
                                   , quantum_efficiency        & ! intent(in)
                                   , curvpar_electron          & ! intent(in)
                                   , qyield_psII               & ! intent(in)
                                   , photosyn_pathway          & ! intent(in)
                                   , dark_respiration_factor   & ! intent(in)
                                   , root_respiration_factor   & ! intent(in)
                                   , electron_transport_factor & ! intent(in)
                                   , triose_phosphate_factor   & ! intent(in)
                                   , water_conductance         & ! intent(in)
                                   , growth_resp_factor        & ! intent(in)
                                   , leaf_turnover_rate        & ! intent(in)
                                   , root_turnover_rate        & ! intent(in)
                                   , storage_turnover_rate     & ! intent(in)
                                   , bark_turnover_rate        & ! intent(in)
                                   , eplastic_vm0              & ! intent(in)
                                   , eplastic_sla              & ! intent(in)
                                   , kplastic_LL               & ! intent(in)
                                   , mort0                     & ! intent(in)
                                   , mort1                     & ! intent(in)
                                   , mort2                     & ! intent(in)
                                   , mort3                     & ! intent(in)
                                   , seedling_mortality        & ! intent(in)
                                   , treefall_s_ltht           & ! intent(in)
                                   , fire_s_min                & ! intent(in)
                                   , fire_s_max                & ! intent(in)
                                   , fire_s_inter              & ! intent(in)
                                   , fire_s_slope              & ! intent(in)
                                   , st_fract                  & ! intent(in)
                                   , r_fract                   & ! intent(in)
                                   , r_bang                    & ! intent(in)
                                   , r_cv50                    & ! intent(in)
                                   , nonlocal_dispersal        & ! intent(in)
                                   , seed_rain                 & ! intent(in)
                                   , f_labile_leaf             & ! intent(in)
                                   , f_labile_stem             & ! intent(in)
                                   , high_psi_threshold        & ! intent(in)
                                   , low_psi_threshold         & ! intent(in)
                                   , leaf_shed_rate            & ! intent(in)
                                   , leaf_grow_rate            & ! intent(in)
                                   , vessel_curl_factor        & ! intent(in)
                                   , leaf_water_cap            & ! intent(in)
                                   , wood_water_cap            & ! intent(in)
                                   , leaf_rwc_min              & ! intent(in)
                                   , wood_rwc_min              & ! intent(in)
                                   , leaf_psi_min              & ! intent(in)
                                   , wood_psi_min              & ! intent(in)
                                   , leaf_psi_osmotic          & ! intent(in)
                                   , wood_psi_osmotic          & ! intent(in)
                                   , leaf_elastic_mod          & ! intent(in)
                                   , wood_elastic_mod          & ! intent(in)
                                   , leaf_psi_tlp              & ! intent(in)
                                   , wood_psi_tlp              & ! intent(in)
                                   , wood_Kmax                 & ! intent(in)
                                   , wood_Kexp                 & ! intent(in)
                                   , wood_psi50                & ! intent(in)
                                   , stoma_lambda              & ! intent(in)
                                   , stoma_beta                & ! intent(in)
                                   , stoma_psi_b               & ! intent(in)
                                   , stoma_psi_c               & ! intent(in)
                                   , C2B                       & ! intent(in)
                                   , hgt_max                   & ! intent(inout)
                                   , repro_min_h               & ! intent(inout)
                                   , dbh_crit                  & ! intent(inout)
                                   , bleaf_crit                & ! intent(inout)
                                   , bdead_crit                & ! intent(inout)
                                   , bevery_crit               & ! intent(inout)
                                   , seed_rain                 & ! intent(inout)
                                   , kplastic_SLA              & ! intent(inout)
                                   , kplastic_vm0              & ! intent(inout)
                                   , Jm_low_temp               & ! intent(inout)
                                   , Jm_high_temp              & ! intent(inout)
                                   , Jm_decay_elow             & ! intent(inout)
                                   , Jm_decay_ehigh            & ! intent(inout)
                                   , Jm0                       & ! intent(inout)
                                   , Jm_hor                    & ! intent(inout)
                                   , Jm_q10                    & ! intent(inout)
                                   , Jm0                       & ! intent(inout)
                                   , TPm0                      & ! intent(inout)
                                   , Rd_low_temp               & ! intent(inout)
                                   , Rd_high_temp              & ! intent(inout)
                                   , Rd_decay_elow             & ! intent(inout)
                                   , Rd_decay_ehigh            & ! intent(inout)
                                   , Rd0                       & ! intent(inout)
                                   , kplastic_rd0              & ! intent(inout)
                                   , Rd_hor                    & ! intent(inout)
                                   , Rd_q10                    & ! intent(inout)
                                   , rrf_low_temp              & ! intent(inout)
                                   , rrf_high_temp             & ! intent(inout)
                                   , rrf_decay_elow            & ! intent(inout)
                                   , rrf_decay_ehigh           & ! intent(inout)
                                   , rrf_hor                   & ! intent(inout)
                                   , rrf_q10                   & ! intent(inout)
                                   , srf_low_temp              & ! intent(inout)
                                   , srf_high_temp             & ! intent(inout)
                                   , srf_decay_elow            & ! intent(inout)
                                   , srf_decay_ehigh           & ! intent(inout)
                                   , srf_hor                   & ! intent(inout)
                                   , srf_q10                   & ! intent(inout) 
                                   , bleaf_crit                & ! intent(inout)
                                   , ddh_allom                 & ! intent(out)
                                   , d1DBH_small               & ! intent(out)
                                   , d2DBH_small               & ! intent(out)
                                   , d1DBH_large               & ! intent(out)
                                   , d2DBH_large               & ! intent(out)
                                   , l1DBH                     & ! intent(out)
                                   , l2DBH                     & ! intent(out)
                                   , balive_crit               & ! intent(out)
                                   , one_plant_c               & ! intent(out)
                                   , min_recruit_size          & ! intent(out)
                                   , min_cohort_size           & ! intent(out)
                                   , negligible_nplant         & ! intent(out)
                                   , repro_min_dbh             & ! intent(out)
                                   , c2n_recruit               & ! intent(out)
                                   , veg_hcap_min              & ! intent(out)
                                   , cleaf                     & ! intent(out)
                                   , csapw                     & ! intent(out)
                                   , cdead                     & ! intent(out)
                                   , cbark                     & ! intent(out)
                                   , dbh_lut                   & ! intent(out)
                                   , bleaf_lut                 & ! intent(out)
                                   , balive_lut                & ! intent(out)
                                   , bdead_lut                 & ! intent(out)
                                   , bevery_lut                & ! intent(out)
                                   , le_mask_lut               & ! intent(out)
                                   , ge_mask_lut               & ! intent(out)
                                   , Vcmax25                   & ! intent(out)
                                   , Jmax25                    & ! intent(out)
                                   , small_rwc_min             & ! intent(out)
                                   , small_psi_min             ! ! intent(out)
   use fusion_fission_coms  , only : ifusion                   & ! intent(in)
                                   , ff_nhgt                   & ! intent(in)
                                   , hgt_class                 ! ! intent(out)
   use allometry            , only : h2dbh                     & ! function
                                   , dbh2h                     & ! function
                                   , size2bl                   & ! function
                                   , size2bd                   & ! function
                                   , size2prd                  ! ! function
   use ed_therm_lib         , only : calc_veg_hcap             ! ! function
   use canopy_radiation_coms, only : ihrzrad                   & ! intent(in)
                                   , cci_hmax                  & ! intent(in)
                                   , leaf_trans_vis            & ! intent(in)
                                   , leaf_reflect_vis          & ! intent(in)
                                   , wood_trans_vis            & ! intent(in)
                                   , wood_reflect_vis          & ! intent(in)
                                   , leaf_trans_nir            & ! intent(in)
                                   , leaf_reflect_nir          & ! intent(in)
                                   , wood_trans_nir            & ! intent(in)
                                   , wood_reflect_nir          & ! intent(in)
                                   , leaf_emiss_tir            & ! intent(in)
                                   , wood_emiss_tir            & ! intent(in)
                                   , orient_factor             & ! intent(in)
                                   , leaf_scatter_vis          & ! intent(out)
                                   , leaf_backscatter_vis      & ! intent(out)
                                   , wood_scatter_vis          & ! intent(out)
                                   , wood_backscatter_vis      & ! intent(out)
                                   , leaf_scatter_nir          & ! intent(out)
                                   , leaf_backscatter_nir      & ! intent(out)
                                   , wood_scatter_nir          & ! intent(out)
                                   , wood_backscatter_nir      & ! intent(out)
                                   , leaf_backscatter_tir      & ! intent(out)
                                   , wood_backscatter_tir      & ! intent(out)
                                   , phi1                      & ! intent(out)
                                   , phi2                      & ! intent(out)
                                   , mu_bar                    ! ! intent(out)
   use rk4_coms             , only : effarea_heat              & ! intent(in)
                                   , effarea_evap              & ! intent(in)
                                   , effarea_transp            ! ! intent(in)
   use decomp_coms          , only : f0_msc                    & ! intent(in)
                                   , f0_psc                    & ! intent(in)
                                   , rh_active_depth           & ! intent(in)
                                   , rh0                       & ! intent(in)
                                   , rh_q10                    & ! intent(in)
                                   , f0_ssc                    & ! intent(out)
                                   , k_rh_active               & ! intent(out)
                                   , rh08                      & ! intent(out)
                                   , rh_q108                   ! ! intent(out)
   use phenology_coms       , only : iphen_scheme              & ! intent(in)
                                   , repro_scheme              & ! intent(in)
                                   , radint                    & ! intent(in)
                                   , radslp                    & ! intent(in)
                                   , radavg_window             & ! intent(in)
                                   , turnamp_window            & ! intent(in)
                                   , turnamp_min               & ! intent(in)
                                   , turnamp_max               & ! intent(in)
                                   , llspan_window             & ! intent(in)
                                   , vm0_window                & ! intent(in)
                                   , sla_window                & ! intent(in)
                                   , radavg_wgt                & ! intent(out)
                                   , turnamp_wgt               & ! intent(out)
                                   , llspan_wgt                & ! intent(out)
                                   , radto_min                 & ! intent(out)
                                   , radto_max                 & ! intent(out)
                                   , vm0_wgt                   & ! intent(out)
                                   , sla_wgt                   ! ! intent(out)
   use farq_leuning         , only : arrhenius                 & ! function
                                   , collatz                   ! ! function
   use plant_hydro          , only : psi2rwc                   & ! function
                                   , rwc2psi                   ! ! function
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   character(len=2)                  :: char_pathway
   integer                           :: ipft
   integer                           :: ihgt
   integer                           :: ilut
   logical                           :: print_zero_table
   real(kind=8)                      :: exp_dbh8
   real(kind=8)                      :: dbh_mult8
   real(kind=8)                      :: dbh_now8
   real                              :: dbh_now
   real                              :: dbh
   real                              :: huge_dbh
   real                              :: huge_height
   real                              :: height_now
   real                              :: bleaf_now
   real                              :: bdead_now
   real                              :: bleaf_min
   real                              :: broot_min
   real                              :: bsapwood_min
   real                              :: bbark_min
   real                              :: balive_min
   real                              :: bdead_min
   real                              :: bstorage_min
   real                              :: bfast_min
   real                              :: bstruct_min
   real                              :: rdepth_min
   real                              :: bleaf_max
   real                              :: broot_max
   real                              :: bsapwood_max
   real                              :: bbark_max
   real                              :: balive_max
   real                              :: bdead_max
   real                              :: bstorage_max
   real                              :: rdepth_max
   real                              :: bleaf_bl
   real                              :: broot_bl
   real                              :: bsapwood_bl
   real                              :: bbark_bl
   real                              :: balive_bl
   real                              :: bdead_bl
   real                              :: bstorage_bl
   real                              :: leaf_hcap_min
   real                              :: wood_hcap_min
   real                              :: lai_min
   real                              :: nplant_res_min
   real                              :: max_hgt_max
   real                              :: leaf_psi_swap
   real                              :: wood_psi_swap
   real                              :: leaf_rwc_small
   real                              :: wood_rwc_small
   real                              :: Rdark25
   real(kind=8)                      :: temp25C8
   real(kind=8)                      :: Vcmax258
   real(kind=8)                      :: Jmax258
   real(kind=8)                      :: Rdark258
   real(kind=8)                      :: refval8
   real(kind=8)                      :: hor8
   real(kind=8)                      :: q108
   real(kind=8)                      :: decay_elow8
   real(kind=8)                      :: decay_ehigh8
   real(kind=8)                      :: low_temp8
   real(kind=8)                      :: high_temp8
   real(kind=8)                      :: lnexplow8
   real(kind=8)                      :: tlow_fun8
   real(kind=8)                      :: lnexphigh8
   real(kind=8)                      :: thigh_fun8
   !----- Local constants. ----------------------------------------------------------------!
   real(kind=4)          , parameter :: kplastic_ref_lai = 4.0 ! used for trait_plasticity == 3
   character(len=str_len), parameter :: zero_table_fn = 'pft_sizes.txt'
   character(len=str_len), parameter :: photo_file    = 'photo_param.txt'
   character(len=str_len), parameter :: allom_file    = 'allom_param.txt'
   character(len=str_len), parameter :: strat_file    = 'strategy_param.txt'
   !----- External functions. -------------------------------------------------------------!
   real                  , external  :: sngloff
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Decide whether to write the table with the sizes.                                 !
   !---------------------------------------------------------------------------------------!
   print_zero_table = btest(idetailed,5)
   !---------------------------------------------------------------------------------------!


   !----- 25 degC, to find Vcmax25. -------------------------------------------------------!
   temp25C8 = 2.5d1 + t008
   !---------------------------------------------------------------------------------------!



   !------ Make sure the soil carbon fractions add up to one. -----------------------------!
   if (f0_msc < 0. .or. f0_psc < 0.0 .or. (f0_msc + f0_psc) > 1.0) then
      write (unit=*,fmt='(a)')          '-------------------------------------------------'
      write (unit=*,fmt='(a)')          ' F0_MSC and F0_PSC must be fractions (0-1)'
      write (unit=*,fmt='(a)')          '    and their sum cannot exceed 1.0'
      write (unit=*,fmt='(a)')          ''
      write (unit=*,fmt='(a)')          ' Current values: '
      write (unit=*,fmt='(a,1x,f12.5)') ' F0_MSC = ',f0_msc
      write (unit=*,fmt='(a,1x,f12.5)') ' F0_PSC = ',f0_psc
      write (unit=*,fmt='(a)')          '-------------------------------------------------'
      write (unit=*,fmt='(a)')          ''
      call fatal_error('Invalid partition of soil C pools. Fix ed_params.f90 or XML.'      &
                      ,'init_derived_params_after_xml','ed_params.f90')
   end if
   f0_ssc = 1.0 - f0_msc - f0_psc
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Determine the bottomost layer to consider for environmental regulation of hetero- !
   ! trophic respiration.                                                                  !
   !---------------------------------------------------------------------------------------!
   k_rh_loop: do k_rh_active=nzg-1,1,-1
      if (slz(k_rh_active) < rh_active_depth) exit k_rh_loop
   end do k_rh_loop
   k_rh_active = k_rh_active + 1
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Double precision of heterotrophic respiration parameters.                         !
   !---------------------------------------------------------------------------------------!
   rh08    = dble(rh0)
   rh_q108 = dble(rh_q10)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Set derived parameters for light-controlled phenology (Kim et al., 2012).        !
   !---------------------------------------------------------------------------------------!
   !----- Radiation weight, the inverse of the window, corrected by radfrq. ---------------!
   radavg_wgt     = radfrq / ( radavg_window * day_sec )
   !----- Turnover weight, the inverse of the window. -------------------------------------!
   turnamp_wgt    = 1. / turnamp_window
   !----- Lifespan weight, the inverse of the window. -------------------------------------!
   llspan_wgt      = 1. / llspan_window
   !----- Vm0 weight, the inverse of the window. ------------------------------------------!
   vm0_wgt         = 1. / vm0_window
   !----- SLA weight, the inverse of the window. ------------------------------------------!
   sla_wgt         = 1. / sla_window
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Minimum and maximum radiation should be defined according to the                 !
   ! economics_scheme, as the model formulation is different.  But set dummy values for    !
   ! both in case this simulation is not using light-controlled phenology.                 !
   !---------------------------------------------------------------------------------------!
   select case (iphen_scheme)
   case (3,5)
      !------------------------------------------------------------------------------------!
      !    Light phenology is enabled.                                                     !
      !------------------------------------------------------------------------------------!
      select case (economics_scheme)
      case (1)
         !---------------------------------------------------------------------------------!
         !      Turnover amplitude is a log-linear function of radiation, based on litter  !
         ! fall data directly related to radiation.                                        !
         !---------------------------------------------------------------------------------!
         radto_min = (turnamp_min / radint) ** (1./radslp)
         radto_max = (turnamp_max / radint) ** (1./radslp)
         !---------------------------------------------------------------------------------!
      case default
         !---------------------------------------------------------------------------------!
         !      Turnover amplitude is a linear function of radiation, like the original    !
         ! approach.                                                                       !
         !---------------------------------------------------------------------------------!
         radto_min       = (turnamp_min - radint) / radslp
         radto_max       = (turnamp_max - radint) / radslp
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!
   case default
      !------------------------------------------------------------------------------------!
      !    Light phenology is disabled.  Set dummy values for minimum and maximum          !
      ! radiation so the turnover amplitude is never calculated.                           !
      !------------------------------------------------------------------------------------!
      radto_min = solar
      radto_max = solar
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find net specific heat for leaves, sapwood, heartwood, and bark.  Here we modify  !
   ! the calculation from earlier versions of ED-2.2.  We now separate oven-dry biomass    !
   ! and internal water for leaves and sapwood, as the latter is dynamic.  For now heart-  !
   ! wood and bark are assumed to have constant internal water content, so we incorporate  !
   ! the specific heat of internal water to the bulk specific heat.                        !
   !     In addition, we ignore the wood-water bonding heat capacity for sapwood.  This is !
   !  a relative small effect (2-5% of heat capacity), and it only changes when water      !
   !  content is very low.  Keeping this term may add some non-linearities and would       !
   ! require dynamic updates --- not a big deal, and we may add this effect in future      !
   ! versions.                                                                             !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   ! Forest Products Laboratory. 2010. Wood handbook -- wood as an engineering material.   !
   !    General Technical Report FPL-GTR-190, U.S. Department of Agriculture, Madison, WI. !
   !    doi:10.2737/FPL-GTR-190 (F10).                                                     !
   !                                                                                       !
   ! Gu L, Meyers T, Pallardy SG, Hanson PJ, Yang B, Heuer M, Hosman KP, Liu Q, Riggs JS,  !
   !    Sluss D et al. 2007. Influences of biomass heat and biochemical energy storages on !
   !    the land surface fluxes and radiative temperature. J. Geophys. Res., 112: D02107.  !
   !    doi:10.1029/2006JD007425 (G07).                                                    !
   !                                                                                       !
   ! Xu X, Medvigy D, Powers JS, Becknell JM , Guan K. 2016. Diversity in plant hydraulic  !
   !    traits explains seasonal and inter-annual variations of vegetation dynamics in     !
   !    seasonally dry tropical forests. New Phytol., 212: 80-95. doi:10.1111/nph.14009    !
   !    (X16).                                                                             !
   !---------------------------------------------------------------------------------------!
   cleaf(:) = c_grn_leaf_dry (:)
   csapw(:) = c_ngrn_wood_dry(:)
   cdead(:) = (c_ngrn_wood_dry(:) + wood_water_sat(:) * cliq) / (1. + wood_water_sat(:))   &
            + delta_c_wood(:)
   cbark(:) = (c_ngrn_bark_dry(:) + bark_water_sat(:) * cliq) / (1. + bark_water_sat(:))   &
            + delta_c_bark(:)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Hgt_max of temperate trees cannot exceed b1Ht, and cannot exceed hgt_ref for     !
   ! tropical trees (IALLOM = 2).                                                          !
   !---------------------------------------------------------------------------------------!
   select case (iallom)
   case (2)
      where(is_tropical(:) .and. (hgt_max(:) >= 0.99 * hgt_ref(:)))
         hgt_max(:) = 0.99 * hgt_ref(:)
      end where
   end select
   where ( (.not. is_tropical(:)) .and. hgt_max(:) >= 0.99 * b1Ht(:))
       hgt_max(:) = 0.99 * b1Ht(:)
   end where
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Make sure that repro_min_h is bounded between hgt_min and hgt_max.  This will     !
   ! avoid floating point exceptions, or surprises when the user only partially sets       !
   ! allometry parameters with XML.                                                        !
   !---------------------------------------------------------------------------------------!
   repro_min_h(:) = merge(repro_min_h(:),hgt_min(:),repro_min_h(:) >= hgt_min(:))
   repro_min_h(:) = merge( repro_min_h(:)                                                  &
                         , merge(hgt_max(:),0.5*(hgt_min(:)+hgt_max(:)),is_grass(:))       &
                         , repro_min_h(:) <= hgt_max(:)                              )
   !---------------------------------------------------------------------------------------!




   !------ Find corresponding DBH. --------------------------------------------------------!
   do ipft=1,n_pft
      repro_min_dbh(ipft) = h2dbh(repro_min_h(ipft),ipft)
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    In case the user does not want reproduction (or in case vegetation dynamics is set !
   ! to zero, set seedling mortality to one, so nothing will become recruit.               !
   !---------------------------------------------------------------------------------------!
   if ( repro_scheme == 0 .or. ivegt_dynamics == 0 ) then
      seedling_mortality(:) = 1.0
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The minimum recruitment size and the recruit carbon to nitrogen ratio.  Both      !
   ! parameters actually depend on which PFT we are solving, since grasses always have     !
   ! significantly less biomass.                                                           !
   !---------------------------------------------------------------------------------------!
   if (print_zero_table) then
      open  (unit=61,file=trim(zero_table_fn),status='replace',action='write')
      write (unit=61,fmt='(40(a,1x))')                '  PFT',        'NAME            '   &
                                              ,'     HGT_MIN','         DBH'               &
                                              ,'  RDEPTH_MIN','   BLEAF_MIN'               &
                                              ,'   BROOT_MIN','BSAPWOOD_MIN'               &
                                              ,'   BBARK_MIN','  BALIVE_MIN'               &
                                              ,'   BDEAD_MIN','BSTORAGE_MIN'               &
                                              ,'    BLEAF_BL','    BROOT_BL'               &
                                              ,' BSAPWOOD_BL','    BBARK_BL'               &
                                              ,'   BALIVE_BL','    BDEAD_BL'               &
                                              ,' BSTORAGE_BL','   BLEAF_MAX'               &
                                              ,'   BROOT_MAX','BSAPWOOD_MAX'               &
                                              ,'   BBARK_MAX','  BALIVE_MAX'               &
                                              ,'   BDEAD_MAX','BSTORAGE_MAX'               &
                                              ,'   INIT_DENS','   SEED_RAIN'               &
                                              ,'MIN_REC_SIZE','MIN_COH_SIZE'               &
                                              ,'         SLA','VEG_HCAP_MIN'               &
                                              ,'     LAI_MIN',' REPRO_MIN_H'               &
                                              ,'     HGT_MAX','REPR_MIN_DBH'               &
                                              ,'    DBH_CRIT', ' RDEPTH_MAX'               &
                                              ,' ONE_PLANT_C',' NEGL_NPLANT'

   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Derive additional parameters.                                                      !
   !---------------------------------------------------------------------------------------!
   do ipft = 1,n_pft
      !----- Re-define the "critical" parameters as the numbers may have changed. ---------!
      dbh_crit    (ipft) = h2dbh(hgt_max(ipft),ipft)
      bleaf_crit  (ipft) = size2bl(dbh_crit(ipft),hgt_max(ipft),SLA(ipft),ipft)
      bdead_crit  (ipft) = size2bd(dbh_crit(ipft),hgt_max(ipft),ipft)
      balive_crit (ipft) = bleaf_crit (ipft)                                               &
                         * (1. + q(ipft) + (qsw(ipft)+qbark(ipft)) * hgt_max(ipft) )
      bevery_crit (ipft) = balive_crit(ipft) + bdead_crit(ipft)
      !------------------------------------------------------------------------------------!


      !----- Set allometric formula. ------------------------------------------------------!
      select case (iallom)
      case (3,4,5)
         ddh_allom(ipft) = is_tropical(ipft) .and. (.not. is_liana(ipft))
      case default
         ddh_allom(ipft) = .false.
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the stem-DBH and leaf-DBH parameters.  These are based on the inverse of  !
      ! the size2bd and size2bl functions, and to be consistent, they cannot be            !
      ! initialised through XML.                                                           !
      !------------------------------------------------------------------------------------!
      if (ddh_allom(ipft)) then
         !---------------------------------------------------------------------------------!
         !    Incorporate both heartwood and height allometric equations to derive DBH.    !
         !---------------------------------------------------------------------------------!
         d2DBH_small(ipft) = 1.  / ( ( 2. + b2Ht(ipft) ) * b2Bs_small(ipft) )
         d1DBH_small(ipft) = ( C2B                                                         &
                             / (b1Bs_small(ipft) * exp(b1Ht(ipft) * b2Bs_small(ipft)) ) )  &
                             ** d2DBH_small(ipft)
         d2DBH_large(ipft) = 1.  / ( ( 2. + b2Ht(ipft) ) * b2Bs_large(ipft) )
         d1DBH_large(ipft) = ( C2B                                                         &
                             / (b1Bs_large(ipft) * exp(b1Ht(ipft) * b2Bs_large(ipft)) ) )  &
                             ** d2DBH_large(ipft)
         !---------------------------------------------------------------------------------!


         !------ Inverse of the leaf biomass function. ------------------------------------!
         l2DBH(ipft) = 1.  / ( ( 2. + b2Ht(ipft) ) * b2Bl(ipft) )
         l1DBH(ipft) = ( 1. / (b1Bl(ipft) * exp(b1Ht(ipft) * b2Bl(ipft)) ) ) ** l2DBH(ipft)
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Bdead is not a direct function of height, simple inversion is sufficient.   !
         !---------------------------------------------------------------------------------!
         d2DBH_small(ipft)  = 1.0  / b2Bs_small(ipft)
         d1DBH_small(ipft)  = (C2B / b1Bs_small(ipft)) ** d2DBH_small(ipft)
         d2DBH_large(ipft)  = 1.0  / b2Bs_large(ipft)
         d1DBH_large(ipft)  = (C2B / b1Bs_large(ipft)) ** d2DBH_large(ipft)
         !---------------------------------------------------------------------------------!


         !------ Inverse of the leaf biomass function. ------------------------------------!
         l2DBH(ipft)  =   1.0 / b2Bl(ipft)
         l1DBH(ipft)  = ( C2B / b1Bl(ipft) ) ** l2DBH(ipft)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!





      !----- Find the DBH and carbon pools associated with a newly formed recruit. --------!
      dbh          = h2dbh(hgt_min(ipft),ipft)
      rdepth_min   = size2prd(hgt_min(ipft),dbh,ipft)
      bleaf_min    = size2bl(dbh,hgt_min(ipft),SLA(ipft),ipft) 
      broot_min    = bleaf_min * q(ipft)
      bsapwood_min = bleaf_min * qsw(ipft)   * hgt_min(ipft)
      bbark_min    = bleaf_min * qbark(ipft) * hgt_min(ipft)
      balive_min   = bleaf_min + broot_min + bsapwood_min + bbark_min
      bdead_min    = size2bd(dbh,hgt_min(ipft),ipft)
      bstorage_min = max(almost_zero,f_bstorage_init(ipft)) * balive_min
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   Find the maximum bleaf and bdead supported.  This is to find the negligible      !
      ! nplant so we ensure that the cohort is always terminated if its mortality rate is  !
      ! very high.                                                                         !
      !------------------------------------------------------------------------------------!
      huge_dbh     = 3. * dbh_crit(ipft)
      huge_height  = dbh2h(ipft, dbh_crit(ipft))
      rdepth_max   = size2prd(hgt_max(ipft),dbh_crit(ipft),ipft)
      bleaf_max    = size2bl(huge_dbh,huge_height,SLA(ipft),ipft)
      broot_max    = bleaf_max * q(ipft)
      bsapwood_max = bleaf_max * qsw(ipft)   * huge_height
      bbark_max    = bleaf_max * qbark(ipft) * huge_height
      balive_max   = bleaf_max + broot_max + bsapwood_max + bbark_max
      bdead_max    = size2bd(huge_dbh,huge_height,ipft)
      bstorage_max = max(almost_zero,f_bstorage_init(ipft)) * balive_max
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Biomass of one individual plant at recruitment.                                 !
      !------------------------------------------------------------------------------------!
      bleaf_bl     = size2bl(dbh_bigleaf(ipft),hgt_min(ipft),SLA(ipft),ipft) 
      broot_bl     = bleaf_bl * q(ipft)
      bsapwood_bl  = bleaf_bl * qsw(ipft)   * hgt_max(ipft)
      bbark_bl     = bleaf_bl * qbark(ipft) * hgt_max(ipft)
      balive_bl    = bleaf_bl + broot_bl + bsapwood_bl + bbark_bl
      bdead_bl     = size2bd(dbh_bigleaf(ipft),hgt_min(ipft),ipft)
      bstorage_bl  = max(almost_zero,f_bstorage_init(ipft)) * balive_bl
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the carbon value for one plant.  If SAS approximation, we define it as    !
      ! the carbon of a seedling, and for the big leaf case we assume the typical big leaf !
      ! plant size.                                                                        !
      !------------------------------------------------------------------------------------!
      select case (ibigleaf)
      case (0)
         one_plant_c(ipft) = bdead_min + balive_min + bstorage_min
      case (1)
         one_plant_c(ipft) = bdead_bl  + balive_bl  + bstorage_bl
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Minimum sizes are allometry-dependent.  IALLOM = 3 defines the sizes based on  !
      ! heat capacity, because this variable considerably affects the model speed.  Other  !
      ! allometry sets define the minimum sizes as before, for back-compability.           !
      !------------------------------------------------------------------------------------!
      select case (iallom)
      case (3,4,5)
         !---------------------------------------------------------------------------------!
         !     New method, each PFT has a minimum resolvable density. The fraction ensures !
         ! that plants start as resolvable.                                                !
         !---------------------------------------------------------------------------------!
         nplant_res_min     = 0.1 * init_density(ipft)
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     The following variable is the minimum heat capacity of either the leaf, or  !
         ! the branches, or the combined pool that is solved by the biophysics.  Value is  !
         ! in J/m2/K.  Because leaves are the pools that can determine the fate of the     !
         ! tree, and all PFTs have leaves (but not branches), we only consider the leaf    !
         ! heat capacity only for the minimum value.                                       !
         !---------------------------------------------------------------------------------!
         call calc_veg_hcap(bleaf_min,agf_bs(ipft)*bdead_min,agf_bs(ipft)*bsapwood_min     &
                           ,agf_bs(ipft)*bbark_min,nplant_res_min,ipft,leaf_hcap_min       &
                           ,wood_hcap_min)
         veg_hcap_min(ipft) = leaf_hcap_min
         lai_min            = nplant_res_min * bleaf_min * sla(ipft)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    The definition of the minimum recruitment size is the minimum amount of      !
         ! biomass in kgC/m2 is available for new recruits.  It does not make much sense   !
         ! to throw in new recruits that cannot be resolved (as they would be initialised  !
         ! and terminated), so we define the initial size to be greater than the minimum   !
         ! resolvable nplant.                                                              !
         !---------------------------------------------------------------------------------!
         min_recruit_size(ipft) = init_density(ipft) * one_plant_c(ipft)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Minimum size (measured as biomass of living and structural tissues) allowed  !
         ! in a cohort.  Cohorts with less biomass than this are going to be terminated.   !
         !---------------------------------------------------------------------------------!
         min_cohort_size(ipft)  = 0.5 * nplant_res_min * one_plant_c(ipft)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Seed_rain is the density of seedling that will be added from somewhere else. !
         ! By default, this variable is initialised as a function of the cohort's minimum  !
         ! size.  Setting the size to be 1/10 of the minimum size to allow the model to    !
         ! try establishing cohorts in different months.                                   !
         !---------------------------------------------------------------------------------! 
         if (seed_rain(ipft) == undef_real) then
            seed_rain(ipft)  = 0.1 * min_recruit_size(ipft) / one_plant_c(ipft)
         end if
         !---------------------------------------------------------------------------------! 
      case default
         !----- Old method, nplant_res_min is a fraction of the smallest initial density. -!
         nplant_res_min     = 0.1 * minval(init_density)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    The definition of the minimum recruitment size is the minimum amount of      !
         ! biomass in kgC/m2 is available for new recruits.  For the time being we use the !
         ! near-bare ground state value as the minimum recruitment size, but this may      !
         ! change depending on how well it goes.                                           !
         !---------------------------------------------------------------------------------!
         min_recruit_size(ipft) = nplant_res_min * one_plant_c(ipft)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Minimum size (measured as biomass of living and structural tissues) allowed  !
         ! in a cohort.  Cohorts with less biomass than this are going to be terminated.   !
         !---------------------------------------------------------------------------------!
         min_cohort_size(ipft)  = 0.1 * min_recruit_size(ipft)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     The following variable is the minimum heat capacity of either the leaf, or  !
         ! the branches, or the combined pool that is solved by the biophysics.  Value is  !
         ! in J/m2/K.  Because leaves are the pools that can determine the fate of the     !
         ! tree, and all PFTs have leaves (but not branches), we only consider the leaf    !
         ! heat capacity only for the minimum value.                                       !
         !---------------------------------------------------------------------------------!
         call calc_veg_hcap(bleaf_min,agf_bs(ipft)*bdead_min,agf_bs(ipft)*bsapwood_min     &
                           ,agf_bs(ipft)*bbark_min,init_density(ipft),ipft,leaf_hcap_min   &
                           ,wood_hcap_min)
         veg_hcap_min(ipft) = onesixth * leaf_hcap_min
         lai_min            = onesixth * init_density(ipft) * bleaf_min * sla(ipft)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Seed_rain is the density of seedling that will be added from somewhere else. !
         ! By default, this variable is initialised as a function of the cohort's minimum  !
         ! size to ensure it allows for reintroduction.  In case this has been initialised !
         ! through xml, then don't change the values.                                      !
         !---------------------------------------------------------------------------------! 
         if (seed_rain(ipft) == undef_real) seed_rain(ipft)  = 0.1 * init_density(ipft)
         !---------------------------------------------------------------------------------! 
      end select
      !------------------------------------------------------------------------------------! 


      !------------------------------------------------------------------------------------! 
      !    The following variable is the absolute minimum cohort population that a cohort  !
      ! can have.  This should be used only to avoid nplant=0, but IMPORTANT: this will    !
      ! lead to a ridiculously small cohort almost guaranteed to be extinct and SHOULD BE  !
      ! USED ONLY IF THE AIM IS TO ELIMINATE THE COHORT.                                   !
      !------------------------------------------------------------------------------------! 
      negligible_nplant(ipft) = onesixth * min_cohort_size(ipft)                           &
                              / (bdead_max + balive_max + bstorage_max)
      !------------------------------------------------------------------------------------! 


      !----- Find the recruit carbon to nitrogen ratio. -----------------------------------!
      bfast_min   = f_labile_leaf(ipft) * ( bleaf_min + broot_min)                         &
                  + f_labile_stem(ipft) * ( bsapwood_min + bbark_min + bdead_min)
      bstruct_min = (1.0 - f_labile_leaf(ipft)) * ( bleaf_min + broot_min )                &
                  + (1.0 - f_labile_stem(ipft)) * ( bsapwood_min + bbark_min + bdead_min )
      c2n_recruit(ipft) = (balive_min + bdead_min + bstorage_min)                          &
                        / ( bfast_min    / c2n_leaf(ipft)                                  &
                          + bstruct_min  / c2n_stem(ipft)                                  &
                          + bstorage_min / c2n_storage    )
      !------------------------------------------------------------------------------------! 


      !------------------------------------------------------------------------------------!
      !     Add PFT parameters to the reference table.                                     !
      !------------------------------------------------------------------------------------! 
      if (print_zero_table) then
         write (unit=61,fmt='(i5,1x,a16,1x,9(f12.8,1x),28(f12.5,1x),1(es12.5,1x))')        &
                                                     ipft,pft_name16(ipft),hgt_min(ipft)   &
                                                    ,dbh,rdepth_min,bleaf_min,broot_min    &
                                                    ,bsapwood_min,bbark_min,balive_min     &
                                                    ,bdead_min,bstorage_min,bleaf_bl       &
                                                    ,broot_bl,bsapwood_bl,bbark_bl         &
                                                    ,balive_bl,bdead_bl,bstorage_bl        &
                                                    ,bleaf_max,broot_max,bsapwood_max      &
                                                    ,bbark_max,balive_max,bdead_max        &
                                                    ,bstorage_max,init_density(ipft)       &
                                                    ,seed_rain(ipft)                       &
                                                    ,min_recruit_size(ipft)                &
                                                    ,min_cohort_size(ipft)                 &
                                                    ,sla(ipft),veg_hcap_min(ipft)          &
                                                    ,lai_min,repro_min_h(ipft)             &
                                                    ,hgt_max(ipft),repro_min_dbh(ipft)     &
                                                    ,dbh_crit(ipft),rdepth_max             &
                                                    ,one_plant_c(ipft)                     &
                                                    ,negligible_nplant(ipft)
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Neatly close the parameter table.                                                  !
   !---------------------------------------------------------------------------------------!
   if (print_zero_table) then
      close (unit=61,status='keep')
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Define the patch fusion layers based on the maximum height amongst PFTs.           !
   !---------------------------------------------------------------------------------------!
   select case (ifusion)
   case (1)
      !----- ED-2.2 fusion scheme.  -------------------------------------------------------!
      select case (ihrzrad)
      case (2,4)
         max_hgt_max   = max(cci_hmax,maxval(hgt_max))
      case default
         max_hgt_max   = maxval(hgt_max)
      end select
      allocate (hgt_class(ff_nhgt))
      do ihgt=1,ff_nhgt
         hgt_class(ihgt) = real(ihgt-1) * 0.055 * max_hgt_max
      end do
      !------------------------------------------------------------------------------------!
   case default
      !------------------------------------------------------------------------------------!
      !      ED-2.1 fusion scheme. Likely problems: it inherently assumes that the tallest !
      ! PFT is at 35m, and it may be too relaxed at the upper canopy.                      !
      !------------------------------------------------------------------------------------!
      allocate (hgt_class(ff_nhgt))
      hgt_class( 1) =  0.0
      hgt_class( 2) =  2.5
      hgt_class( 3) =  7.5
      hgt_class( 4) = 12.5
      hgt_class( 5) = 17.5
      hgt_class( 6) = 22.5
      hgt_class( 7) = 27.5
      hgt_class( 8) = 32.5
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Scattering coefficients.  Contrary to ED-2.1, these values are based on the       !
   ! description by by Sellers (1985) and the CLM technical manual, which includes the     !
   ! leaf orientation factor in the backscattering.  This DOES NOT reduce to ED-2.1 case   !
   ! when the leaf orientation is random.                                                  !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Forward scattering.                                                               !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   leaf_scatter_vis(:) = leaf_reflect_vis(:) + leaf_trans_vis(:)
   wood_scatter_vis(:) = wood_reflect_vis(:) + wood_trans_vis(:)
   !----- Near infrared (NIR). ------------------------------------------------------------!
   leaf_scatter_nir(:) = leaf_reflect_nir(:) + leaf_trans_nir(:)
   wood_scatter_nir(:) = wood_reflect_nir(:) + wood_trans_nir(:)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Back-scattering coefficients following CLM.                                      !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   leaf_backscatter_vis(:) = ( leaf_scatter_vis(:)                                         &
                           + 2.5d-1 * ( leaf_reflect_vis(:) - leaf_trans_vis(:)   )        &
                           * ( 1.d0 + orient_factor(:)) ** 2 )                             &
                           / ( 2.d0 * leaf_scatter_vis )
   wood_backscatter_vis(:) = ( wood_scatter_vis(:)                                         &
                          + 2.5d-1 * ( wood_reflect_vis(:) - wood_trans_vis(:)   )         &
                          * ( 1.d0 + orient_factor(:)) ** 2 )                              &
                          / ( 2.d0 * wood_scatter_vis )
   !----- Near infrared (NIR). ------------------------------------------------------------!
   leaf_backscatter_nir(:) = ( leaf_scatter_nir(:)                                         &
                           + 2.5d-1 * ( leaf_reflect_nir(:) - leaf_trans_nir(:)   )        &
                           * ( 1.d0 + orient_factor(:)) ** 2 )                             &
                           / ( 2.d0 * leaf_scatter_nir )
   wood_backscatter_nir(:) = ( wood_scatter_nir(:)                                         &
                          + 2.5d-1 * ( wood_reflect_nir(:) - wood_trans_nir(:)   )         &
                          * ( 1.d0 + orient_factor(:)) ** 2 )                              &
                          / ( 2.d0 * wood_scatter_nir )
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Thermal scattering coefficients.  Contrary to ED-2.1, these values are based on   !     !
   ! the description by S85 and O13 and the CLM technical manual, which includes the leaf  !
   ! orientation factor in the backscattering.  It also assumes that transmittance is      !
   ! zero, similarly to ZQ06.  This DOES NOT reduce to ED-2.1 case when the leaf           !
   ! orientation is random.  The source of the original parameters is not clear, though.   !
   !---------------------------------------------------------------------------------------!
   leaf_backscatter_tir(:) = 5.d-1 + 1.25d-1 * (1 + orient_factor(:)) ** 2
   wood_backscatter_tir(:) = 5.d-1 + 1.25d-1 * (1 + orient_factor(:)) ** 2
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Light extinction coefficients.   These are found following CLM technical manual,  !
   ! and the values fall back to ED-2.0 defaults when orient_factor is zero.               !
   !---------------------------------------------------------------------------------------!
   phi1(:) = 5.d-1 - orient_factor(:) * ( 6.33d-1 + 3.3d-1 * orient_factor(:) )
   phi2(:) = 8.77d-1 * (1.d0 - 2.d0 * phi1(:))
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Find the average inverse diffuse optical depth per unit leaf and stem area.       !
   ! We follow CLM technical manual, equation 3.4 only when the orientation factor is      !
   ! non-zero.   Otherwise, we make it 1.d0, which is the limit of that equation when      !
   ! phi2 approaches zero.                                                                 !
   !---------------------------------------------------------------------------------------!
   where (orient_factor(:) == 0.d0)
      mu_bar(:) = 1.d0
   elsewhere
      mu_bar(:) = ( 1.d0 - phi1(:) * log(1.d0 + phi2(:) / phi1(:)) / phi2(:) ) / phi2(:)
   end where
   !---------------------------------------------------------------------------------------!







   !---------------------------------------------------------------------------------------!
   !       The electron transport and respiration terms may or may not be the same as the  !
   ! the Vm terms.  In case they are not provided through XML and have not been previously !
   ! assigned, copy over the Vm-equivalent parameters.                                     !
   !---------------------------------------------------------------------------------------!
   where (Jm_low_temp (:) == undef_real)
      Jm_low_temp (:) = Vm_low_temp (:)
   end where
   where (Jm_high_temp(:) == undef_real)
      Jm_high_temp(:) = Vm_high_temp(:)
   end where
   where (Jm_decay_elow(:) == undef_real)
      Jm_decay_elow(:) = Vm_decay_elow(:)
   end where
   where (Jm_decay_ehigh(:) == undef_real)
      Jm_decay_ehigh(:) = Vm_decay_ehigh(:)
   end where
   where (Jm_hor      (:) == undef_real)
      Jm_hor      (:) = Vm_hor      (:)
   end where
   where (Jm_q10      (:) == undef_real)
      Jm_q10      (:) = 0.8 * Vm_q10      (:)
      ! Jm usually shows reduced temperature sensitivity
      ! than Vcmax (e.g. Slot et al. 2017 PCE)
   end where

   where (Rd_low_temp (:) == undef_real)
      Rd_low_temp (:) = Vm_low_temp (:)
   end where
   where (Rd_high_temp(:) == undef_real)
      Rd_high_temp(:) = Vm_high_temp(:)
   end where
   where (Rd_decay_elow(:) == undef_real)
      Rd_decay_elow(:) = Vm_decay_elow(:)
   end where
   where (Rd_decay_ehigh(:) == undef_real)
      Rd_decay_ehigh(:) = Vm_decay_ehigh(:)
   end where
   where (Rd_hor      (:) == undef_real)
      Rd_hor      (:) = Vm_hor      (:)
   end where
   where (Rd_q10      (:) == undef_real)
      Rd_q10      (:) = Vm_q10      (:)
   end where
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Reference values for electron transport, triose phosphate utilisation, and dark    !
   ! respiration.  We only fill those that have not been initialised trhough XML.  We must !
   ! also check the photosynthesis scheme.  Unless IPHYSIOL is set to 0 (original ED-2.0), !
   ! we take differences in Q10 into account to correct for the fact that these ratios are !
   ! defined for T=25C whereas ED2 uses T=15C as reference.  Although IPHYSIOL also uses   !
   ! the Arrhenius formulation, we use Q10 for correction because it is simpler, feel free !
   ! to correct this, or simply populate Q10 with the sought correction factor as Q10 will !
   ! not be used for anything else.  Note that we currently impose that Tpm:Vm ratio is    !
   ! constant, and thus we do not have tpm_q10.  This could change in the future in case   !
   ! evidence for different Q10 factor appears in the literature.                          !
   !---------------------------------------------------------------------------------------!
   select case (iphysiol)
   case (0)
      where (Jm0(:) == undef_real)
         Jm0(:) = electron_transport_factor(:) * Vm0(:)
      end where
      where (Rd0(:) == undef_real)
         Rd0(:) = dark_respiration_factor(:)   * Vm0(:)
      end where
      where (TPm0(:) == undef_real)
         TPm0(:) = triose_phosphate_factor(:)  * Vm0(:)
      end where
   case default
      where (Jm0(:) == undef_real)
         Jm0(:) = electron_transport_factor(:) * Vm0(:) * vm_q10(:) / jm_q10(:)
      end where
      where (Rd0(:) == undef_real)
         Rd0(:) = dark_respiration_factor(:)   * Vm0(:) * vm_q10(:) / rd_q10(:)
      end where
      where (TPm0(:) == undef_real)
         TPm0(:) = triose_phosphate_factor(:)  * Vm0(:)
      end where
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     In case root respiration parameters have not been initialised through XML, we     !
   ! follow Moorcroft et al. (2001) and assign the same numbers used for leaf respiration. !
   !---------------------------------------------------------------------------------------!
   !----- Temperature [degC] below which root metabolic activity rapidly declines. --------!
   where (rrf_low_temp(:)  == undef_real)
      rrf_low_temp(:)  = rd_low_temp(:)
   end where
   !----- Temperature [degC] above which root metabolic activity rapidly declines. --------!
   where (rrf_high_temp(:) == undef_real)
      rrf_high_temp(:) = rd_high_temp(:)
   end where
   !----- Decay factor for the exponential correction. ------------------------------------!
   where (rrf_decay_elow(:)   == undef_real)
      rrf_decay_elow(:)   = rd_decay_elow(:)
   end where
   !----- Decay factor for the exponential correction. ------------------------------------!
   where (rrf_decay_ehigh(:)   == undef_real)
      rrf_decay_ehigh(:)  = rd_decay_ehigh(:)
   end where
   !----- Exponent for Rr in the Arrhenius equation [K]. ----------------------------------!
   where (rrf_hor(:)       == undef_real)
      rrf_hor(:)       = rd_hor(:)
   end where
   !----- Base (Q10 term) for respiration in Collatz equation . ---------------------------!
   where (rrf_q10(:)       == undef_real)
      rrf_q10(:)       = rd_q10(:)
   end where
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     In case stem respiration parameters have not been initialised through XML, we     !
   ! assign the same numbers used for leaf respiration.                                    !
   !---------------------------------------------------------------------------------------!
   !----- Temperature [degC] below which stem metabolic activity rapidly declines. --------!
   where (srf_low_temp(:)  == undef_real)
      srf_low_temp(:)  = rd_low_temp(:)
   end where
   !----- Temperature [degC] above which stem metabolic activity rapidly declines. --------!
   where (srf_high_temp(:) == undef_real)
      srf_high_temp(:) = rd_high_temp(:)
   end where
   !----- Decay factor for the exponential correction. ------------------------------------!
   where (srf_decay_elow(:)   == undef_real)
      srf_decay_elow(:)   = rd_decay_elow(:)
   end where
   !----- Decay factor for the exponential correction. ------------------------------------!
   where (srf_decay_ehigh(:)   == undef_real)
      srf_decay_ehigh(:)  = rd_decay_ehigh(:)
   end where
   !----- Exponent for Rr in the Arrhenius equation [K]. ----------------------------------!
   where (srf_hor(:)       == undef_real)
      srf_hor(:)       = rd_hor(:)
   end where
   !----- Base (Q10 term) for respiration in Collatz equation . ---------------------------!
   where (srf_q10(:)       == undef_real)
      srf_q10(:)       = rd_q10(:)
   end where
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find carboxylation rate, dark respiration, and electron transport rate at 25 degC.!
   !---------------------------------------------------------------------------------------!
   do ipft=1,n_pft
      !------ Convert values to double precision to leverage the existing functions. ------!
      refval8      = dble(Vm0           (ipft))
      hor8         = dble(vm_hor        (ipft))
      q108         = dble(vm_q10        (ipft))
      decay_elow8  = dble(vm_decay_elow (ipft))
      decay_ehigh8 = dble(vm_decay_ehigh(ipft))
      low_temp8    = dble(vm_low_temp   (ipft)) + t008
      high_temp8   = dble(vm_high_temp  (ipft)) + t008
      !------------------------------------------------------------------------------------!



      !----- Find the physiology-scheme dependent, uncorrected Vcmax at 25degC.  ----------!
      select case (iphysiol)
      case (0,1)
         Vcmax258 = arrhenius(temp25C8,refval8,hor8)
      case (2,3)
         Vcmax258 = collatz(temp25C8,refval8,q108)
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the temperature correction factors to shut down Vm at extreme temper-      !
      ! atures.  This should have minimal effect at 25degC, though.  It is unlikely, but   !
      ! to avoid floating point exceptions, we also check whether the temperature will     !
      ! make the exponential too large or too small.                                       !
      !------------------------------------------------------------------------------------!
      !----- Low temperature. -------------------------------------------------------------!
      lnexplow8  = decay_elow8  * (low_temp8 - temp25C8)
      lnexplow8  = max(lnexp_min8,min(lnexp_max8,lnexplow8))
      tlow_fun8  = 1.d0 +  exp(lnexplow8)
      !----- High temperature. ------------------------------------------------------------!
      lnexphigh8 = decay_ehigh8 * (temp25C8 - high_temp8)
      lnexphigh8 = max(lnexp_min8,min(lnexp_max8,lnexphigh8))
      thigh_fun8 = 1.d0 + exp(lnexphigh8)
      !----- Correct Vcmax. ---------------------------------------------------------------!
      Vcmax258   = Vcmax258 / (tlow_fun8 * thigh_fun8)
      !------------------------------------------------------------------------------------!


      !------ Convert to single precision. ------------------------------------------------!
      Vcmax25(ipft) = sngloff(Vcmax258,tiny_num8)
      !------------------------------------------------------------------------------------!

      !------ Convert values to double precision to leverage the existing functions. ------!
      refval8      = dble(Rd0           (ipft))
      hor8         = dble(Rd_hor        (ipft))
      q108         = dble(Rd_q10        (ipft))
      decay_elow8  = dble(Rd_decay_elow (ipft))
      decay_ehigh8 = dble(Rd_decay_ehigh(ipft))
      low_temp8    = dble(Rd_low_temp   (ipft)) + t008
      high_temp8   = dble(Rd_high_temp  (ipft)) + t008
      !------------------------------------------------------------------------------------!



      !----- Find the physiology-scheme dependent, uncorrected Vcmax at 25degC.  ----------!
      select case (iphysiol)
      case (0,1)
         Rdark258 = arrhenius(temp25C8,refval8,hor8)
      case (2,3)
         Rdark258 = collatz(temp25C8,refval8,q108)
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the temperature correction factors to shut down Rd at extreme temper-      !
      ! atures.  This should have minimal effect at 25degC, though.  It is unlikely, but   !
      ! to avoid floating point exceptions, we also check whether the temperature will     !
      ! make the exponential too large or too small.                                       !
      !------------------------------------------------------------------------------------!
      !----- Low temperature. -------------------------------------------------------------!
      lnexplow8  = decay_elow8  * (low_temp8 - temp25C8)
      lnexplow8  = max(lnexp_min8,min(lnexp_max8,lnexplow8))
      tlow_fun8  = 1.d0 +  exp(lnexplow8)
      !----- High temperature. ------------------------------------------------------------!
      lnexphigh8 = decay_ehigh8 * (temp25C8 - high_temp8)
      lnexphigh8 = max(lnexp_min8,min(lnexp_max8,lnexphigh8))
      thigh_fun8 = 1.d0 + exp(lnexphigh8)
      !----- Correct Vcmax. ---------------------------------------------------------------!
      Rdark258   = Rdark258 / (tlow_fun8 * thigh_fun8)
      !------------------------------------------------------------------------------------!


      !------ Convert to single precision. ------------------------------------------------!
      Rdark25    = sngloff(Rdark258,tiny_num8)
      !------------------------------------------------------------------------------------!


      !------ Convert values to double precision to leverage the existing functions. ------!
      refval8      = dble(Jm0           (ipft))
      hor8         = dble(jm_hor        (ipft))
      q108         = dble(jm_q10        (ipft))
      decay_elow8  = dble(jm_decay_elow (ipft))
      decay_ehigh8 = dble(jm_decay_ehigh(ipft))
      low_temp8    = dble(jm_low_temp   (ipft)) + t008
      high_temp8   = dble(jm_high_temp  (ipft)) + t008
      !------------------------------------------------------------------------------------!



      !----- Find the physiology-scheme dependent, uncorrected Jmax at 25degC.  -----------!
      select case (iphysiol)
      case (0,1)
         Jmax258 = arrhenius(temp25C8,refval8,hor8)
      case (2,3)
         Jmax258 = collatz(temp25C8,refval8,q108)
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the temperature correction factors to shut down Vm at extreme temper-      !
      ! atures.  This should have minimal effect at 25degC, though.  It is unlikely, but   !
      ! to avoid floating point exceptions, we also check whether the temperature will     !
      ! make the exponential too large or too small.                                       !
      !------------------------------------------------------------------------------------!
      !----- Low temperature. -------------------------------------------------------------!
      lnexplow8  = decay_elow8  * (low_temp8 - temp25C8)
      lnexplow8  = max(lnexp_min8,min(lnexp_max8,lnexplow8))
      tlow_fun8  = 1.d0 +  exp(lnexplow8)
      !----- High temperature. ------------------------------------------------------------!
      lnexphigh8 = decay_ehigh8 * (temp25C8 - high_temp8)
      lnexphigh8 = max(lnexp_min8,min(lnexp_max8,lnexphigh8))
      thigh_fun8 = 1.d0 + exp(lnexphigh8)
      !----- Correct Vcmax. ---------------------------------------------------------------!
      Jmax258    = Jmax258 / (tlow_fun8 * thigh_fun8)
      !------------------------------------------------------------------------------------!


      !------ Convert to single precision. ------------------------------------------------!
      Jmax25(ipft) = sngloff(Jmax258,tiny_num8)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    In case kplastic_vm0 was not initialised through XML, use the function from     !
      ! L10.  This is not the same fit, instead it was found using robust standardised     !
      ! major axis to reduce leverage.                                                     !
      !                                                                                    !
      ! Reference:                                                                         !
      !                                                                                    !
      ! Lloyd J, Patino S, Paiva RQ, Nardoto GB, Quesada CA, Santos AJB, Baker TR,         !
      !    Brand WA, Hilke I, Gielmann H, Raessler M, Luizao FJ, Martinelli LA,            !
      !    Mercado LM. 2009. Optimisation of photosynthetic carbon gain and within-canopy  !
      !    gradients of associated foliar traits for Amazon forest trees. Biogeosciences   !
      !    7(6): 1833-1859, doi:10.5194/bg-7-1833-2010 (L10).                              !
      !------------------------------------------------------------------------------------!
      if (kplastic_vm0(ipft) == undef_real) then
         !----- Use the "low" variables as placeholders for double precision. -------------!
         lnexplow8 = -2.788d0 + 1.439d-2 * Vcmax258
         lnexplow8 = max(lnexp_min8,min(lnexp_max8,lnexplow8))
         !---- Set the value to negative so Vm0 decreases in the understorey. -------------!
         tlow_fun8 = - 1.d0 * exp(lnexplow8)
         !---------------------------------------------------------------------------------!

         !----- Set kplastic for Vm0. -----------------------------------------------------!
         kplastic_vm0(ipft) = sngloff(tlow_fun8,tiny_num8)
         !---------------------------------------------------------------------------------!
         
         !---------------------------------------------------------------------------------!
         ! rewrite plasticity for tropical trees based on BCI data if 
         ! trait_plasticity_scheme is 3
         !---------------------------------------------------------------------------------!
         select case (trait_plasticity_scheme)
         case (3)
            if (is_tropical(ipft) .and. (.not. is_grass(ipft))) then
                kplastic_vm0(ipft) = - 1.0 * (0.811 * log(Vcmax25(ipft)) - 2.22)           &
                                   / kplastic_ref_lai
            endif
         end select
         !---------------------------------------------------------------------------------!

      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    In case kplastic_rd0 was not initialised through XML, use the function from     !
      ! L10.  This is not the same fit, instead it was found using robust standardised     !
      ! major axis to reduce leverage.                                                     !
      !                                                                                    !
      ! Reference: Unpublished trait data from BCI                                         !
      !                                                                                    !
      !------------------------------------------------------------------------------------!
      if (kplastic_rd0(ipft) == undef_real) then


         !----- Default: assume the same as the Vcmax decay. ------------------------------!
         kplastic_rd0(ipft) = kplastic_vm0(ipft)
         !---------------------------------------------------------------------------------!
 

         !---------------------------------------------------------------------------------!
         !     Rewrite plasticity for tropical trees based on BCI data in case             !
         ! trait_plasticity_scheme is 3.                                                   !
         !---------------------------------------------------------------------------------!
         select case (trait_plasticity_scheme)
         case (3)
            !----- Make sure this is applied to tropical trees only. ----------------------!
            if (is_tropical(ipft) .and. (.not. is_grass(ipft))) then
                kplastic_rd0(ipft) = - 1.0 * (0.559 * log(Rdark25) + 0.82)                 &
                                   / kplastic_ref_lai
            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    In case kplastic_SLA was not initialised through XML, use an empirical          !
      ! function.  This function is obtained using the data from KN2016, assuming Beer's   !
      ! law to find the LAI that would cause the drop of light conditions from 40 mol/m2/d !
      ! to 3 mol/m2/d (their standard light conditions), and assuming that changes in SLA  !
      ! would follow a exponential decay/expansion, akin to Vm0.  Only tropical trees were !
      ! used for the standardised major axis regression, and we used a linear-log function !
      ! because kplastic_SLA may switch sign: in fact KN2016 observations had negative     !
      ! values for very high-canopy-level SLA values.                                      !
      !                                                                                    !
      ! Reference:                                                                         !
      !                                                                                    !
      ! Keenan TF, Niinemets U. 2016. Global leaf trait estimates biased due to plasticity !
      !    in the shade. Nat. Plants, 3:16201, Dec 2016. doi:10.1038/nplants.2016.201      !
      !    (KN16).                                                                         !
      !------------------------------------------------------------------------------------!
      if (kplastic_SLA(ipft) == undef_real) then
         kplastic_SLA(ipft) = 0.462 - 0.1239 * log(SLA(ipft))
            
         !---------------------------------------------------------------------------------!
         ! rewrite plasticity for tropical trees based on BCI data if 
         ! trait_plasticity_scheme is 3
         !---------------------------------------------------------------------------------!
         select case (trait_plasticity_scheme)
         case (3)
            if (is_tropical(ipft) .and. (.not. is_grass(ipft))) then
                kplastic_SLA(ipft) = (0.214 * log(1. / SLA(ipft) * 2000.) - 0.088)         &
                                   / kplastic_ref_lai
            endif
         end select
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    In case kplastic_LL was not initialised through XML, use an empirical           !
      ! function.  This function is obtained using digitsed data from RK2016 (Figure 5a),  !
      ! assuming Beer's law to find the LAI that would cause the drop of light conditions  !
      ! to 0.8% of the full illumination  (their reference light levels), and assuming     !
      ! that changes in leaf longevity (LL) would follow a exponential decay/expansion,    !
      ! akin to Vm0.  We used a linear-log function because kplastic_LL may switch sign.   !
      ! It did not happen in RK2016 but some individuals had near-zero decay, and none of  !
      ! their species had extremely high SLA like KN2016.                                  !
      !                                                                                    !
      ! Reference:                                                                         !
      !                                                                                    !
      ! Russo S, Kitajima K. 2016. The ecophysiology of leaf lifespan in tropical forests: !
      !    Adaptive and plastic responses to environmental heterogeneity. In: Goldstein G, !
      !    Santiago, LS, eds. Tropical Tree Physiology: Adaptations and Responses in a     !
      !    Changing Environment, 357-383. Springer International Publishing, Cham,         !
      !    Switzerland. doi:10.1007/978-3-319-27422- 5 17.                                 !
      !------------------------------------------------------------------------------------!
      if (kplastic_LL(ipft) == undef_real .and. leaf_turnover_rate(ipft) > 0.) then
         kplastic_LL(ipft) = 0.2126 - 0.062 * log(12./leaf_turnover_rate(ipft))
         
         !---------------------------------------------------------------------------------!
         ! rewrite plasticity for tropical trees based on BCI data if 
         ! trait_plasticity_scheme is 3
         !---------------------------------------------------------------------------------!
         select case (trait_plasticity_scheme)
         case (3)
            if (is_tropical(ipft) .and. (.not. is_grass(ipft))) then
                kplastic_LL(ipft) = -1.0 * (0.504 * log(Vcmax25(ipft) * SLA(ipft) / 2000.) &
                                           - 0.401) / kplastic_ref_lai
            endif
         end select

      elseif (kplastic_LL(ipft) == undef_real) then
         kplastic_LL(ipft) = 0.0
      end if

      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Make safe small-tree thresholds.                                               !
      !------------------------------------------------------------------------------------!
      !----- Find the "swapped" minimum psi (the swapped order below is not a bug). -------!
      call rwc2psi(wood_rwc_min(ipft),leaf_rwc_min(ipft),ipft,leaf_psi_swap,wood_psi_swap)
      small_psi_min(ipft) = max( leaf_psi_min(ipft),wood_psi_min(ipft)                     &
                               , leaf_psi_swap     ,wood_psi_swap      )
      call psi2rwc(small_psi_min(ipft),small_psi_min(ipft),ipft                            &
                  ,leaf_rwc_small,wood_rwc_small)
      small_rwc_min(ipft) = min(leaf_rwc_small,wood_rwc_small)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!




   !----- Build the look-up table for btotal. ---------------------------------------------!
   if (nbt_lut < 2) then
      write(unit=*,fmt='(a,1x,i5)') ' NBT_LUT = ',nbt_lut
      call fatal_error(' NBT_LUT must be at least 2.','init_derived_params_after_xml'      &
                      ,'ed_params.f90')
   else
      allocate(dbh_lut    (nbt_lut,n_pft))
      allocate(bleaf_lut  (nbt_lut,n_pft))
      allocate(balive_lut (nbt_lut,n_pft))
      allocate(bdead_lut  (nbt_lut,n_pft))
      allocate(bevery_lut (nbt_lut,n_pft))
      allocate(le_mask_lut(nbt_lut      ))     ! Aux variable used by b?2dbh
      allocate(ge_mask_lut(nbt_lut      ))     ! Aux variable used by b?2dbh
      exp_dbh8 = 1.d0 / (dble(nbt_lut) - 1.d0)
      do ipft=1,n_pft
         dbh_mult8 = (dble(dbh_crit(ipft))/dble(min_dbh(ipft))) **exp_dbh8
         do ilut = 1, nbt_lut
            dbh_now8               = dble(min_dbh(ipft)) * dbh_mult8 ** (ilut-1)
            dbh_now                = sngloff(dbh_now8,tiny_num8)
            height_now             = dbh2h(ipft, dbh_now)
            bleaf_now              = size2bl(dbh_now,height_now,SLA(ipft),ipft)
            ! NOTE:
            ! When IALLOM == 4, bleaf becomes plastic to light conditions
            ! bleaf_now would therefore become inaccurate (usually biased high)
            ! This would underestimate the fraction used for structural growth at monthly scale
            ! But should be ok at annual scale.
            bdead_now              = size2bd(dbh_now,height_now,ipft)
            dbh_lut    (ilut,ipft) = dbh_now
            bleaf_lut  (ilut,ipft) = bleaf_now
            balive_lut (ilut,ipft) = bleaf_now                                             &
                                   * (1. + q(ipft) + (qsw(ipft) + qbark(ipft))*height_now)
            bdead_lut  (ilut,ipft) = bdead_now
            bevery_lut (ilut,ipft) = balive_lut(ilut,ipft) + bdead_lut(ilut,ipft)
         end do
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Check whether lianas were included in this simulation.                            !
   !---------------------------------------------------------------------------------------!
   lianas_included = .false.
   do ipft=1,n_pft
      lianas_included = lianas_included .or. (include_pft(ipft) .and. is_liana(ipft))
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Decide which variables to write in the output based on iphysiol.                   !
   !---------------------------------------------------------------------------------------!
   if (print_zero_table) then
      open (unit=18,file=trim(photo_file),status='replace',action='write')
      select case (iphysiol)
      case (0,1)
         !---------------------------------------------------------------------------------!
         !     Arrhenius-based model, print Arrhenius reference and skip Q10.              !
         !---------------------------------------------------------------------------------!
         write(unit=18,fmt='(44(1x,a))') '         PFT','    TROPICAL','       GRASS'      &
                                        ,'     PATHWAY','         SLA','          D0'      &
                                        ,'     VCMAX25','         VM0','      JMAX25'      &
                                        ,'         JM0','        TPM0','         RD0'      &
                                        ,'     RD0:VM0','     JM0:VM0','    TPM0:VM0'      &
                                        ,'    VM_TCOLD','     VM_THOT','  VM_EXP_LOW'      &
                                        ,' VM_EXP_HIGH','      VM_HOR','    JM_TCOLD'      &
                                        ,'     JM_THOT','  JM_EXP_LOW',' JM_EXP_HIGH'      &
                                        ,'      VM_HOR','    RD_TCOLD','     RD_THOT'      &
                                        ,'  RD_EXP_LOW',' RD_EXP_HIGH','      RD_HOR'      &
                                        ,'    ST_SLOPE','         GS0','Q_EFFICIENCY'      &
                                        ,' CV_ELECTRON',' QYIELD_PSII','      KWROOT'      &
                                        ,'  LEAF_WIDTH','       RRFF0','KPLASTIC_VM0'      &
                                        ,'KPLASTIC_RD0','KPLASTIC_SLA',' KPLASTIC_LL'      &
                                        ,'EPLASTIC_VM0','EPLASTIC_SLA'
                                        
         do ipft=1,n_pft
            write(char_pathway,fmt='(a,i1)') 'C',photosyn_pathway(ipft)

            write (unit=18,fmt='(8x,i5,2(12x,l1),11x,a2,40(1x,f12.6))')                    &
                           ipft,is_tropical(ipft),is_grass(ipft),char_pathway,SLA(ipft)    &
                          ,D0(ipft),Vcmax25(ipft),Vm0(ipft),Jmax25(ipft),Jm0(ipft)         &
                          ,TPm0(ipft),Rd0(ipft),dark_respiration_factor(ipft)              &
                          ,electron_transport_factor(ipft),triose_phosphate_factor(ipft)   &
                          ,Vm_low_temp(ipft),Vm_high_temp(ipft),Vm_decay_elow(ipft)        &
                          ,Vm_decay_ehigh(ipft),vm_hor(ipft),Jm_low_temp(ipft)             &
                          ,Jm_high_temp(ipft),Jm_decay_elow(ipft),Jm_decay_ehigh(ipft)     &
                          ,jm_hor(ipft),Rd_low_temp(ipft),Rd_high_temp(ipft)               &
                          ,Rd_decay_elow(ipft),Rd_decay_ehigh(ipft),Rd_hor(ipft)           &
                          ,stomatal_slope(ipft),cuticular_cond(ipft)                       &
                          ,quantum_efficiency(ipft),curvpar_electron(ipft)                 &
                          ,qyield_psII(ipft),water_conductance(ipft)*yr_sec                &
                          ,leaf_width(ipft),root_respiration_factor(ipft)                  &
                          ,kplastic_vm0(ipft),kplastic_rd0(ipft)                           &
                          ,kplastic_sla(ipft),kplastic_LL(ipft)                            &
                          ,eplastic_vm0(ipft),eplastic_sla(ipft)
         end do
         !---------------------------------------------------------------------------------!
      case (2,3)
         !---------------------------------------------------------------------------------!
         !     Collatz-based model, print Q10 instead of Arrhenius reference.              !
         !---------------------------------------------------------------------------------!
         write(unit=18,fmt='(44(1x,a))') '         PFT','    TROPICAL','       GRASS'      &
                                        ,'     PATHWAY','         SLA','          D0'      &
                                        ,'     VCMAX25','         VM0','      JMAX25'      &
                                        ,'         JM0','        TPM0','         RD0'      &
                                        ,'     RD0:VM0','     JM0:VM0','    TPM0:VM0'      &
                                        ,'    VM_TCOLD','     VM_THOT','  VM_EXP_LOW'      &
                                        ,' VM_EXP_HIGH','      VM_Q10','    JM_TCOLD'      &
                                        ,'     JM_THOT','  JM_EXP_LOW',' JM_EXP_HIGH'      &
                                        ,'      JM_Q10','    RD_TCOLD','     RD_THOT'      &
                                        ,'  RD_EXP_LOW',' RD_EXP_HIGH','      RD_Q10'      &
                                        ,'    ST_SLOPE','         GS0','Q_EFFICIENCY'      &
                                        ,' CV_ELECTRON',' QYIELD_PSII','      KWROOT'      &
                                        ,'  LEAF_WIDTH','       RRFF0','KPLASTIC_VM0'      &
                                        ,'KPLASTIC_RD0','KPLASTIC_SLA',' KPLASTIC_LL'      &
                                        ,'EPLASTIC_VM0','EPLASTIC_SLA'
         do ipft=1,n_pft
            write(char_pathway,fmt='(a,i1)') 'C',photosyn_pathway(ipft)

            write (unit=18,fmt='(8x,i5,2(12x,l1),11x,a2,39(1x,f12.6))')                    &
                           ipft,is_tropical(ipft),is_grass(ipft),char_pathway,SLA(ipft)    &
                          ,D0(ipft),Vcmax25(ipft),Vm0(ipft),Jmax25(ipft),Jm0(ipft)         &
                          ,TPm0(ipft),Rd0(ipft),dark_respiration_factor(ipft)              &
                          ,electron_transport_factor(ipft),triose_phosphate_factor(ipft)   &
                          ,Vm_low_temp(ipft),Vm_high_temp(ipft),Vm_decay_elow(ipft)        &
                          ,Vm_decay_ehigh(ipft),vm_q10(ipft),Jm_low_temp(ipft)             &
                          ,Jm_high_temp(ipft),Jm_decay_elow(ipft),Jm_decay_ehigh(ipft)     &
                          ,jm_q10(ipft),Rd_low_temp(ipft),Rd_high_temp(ipft)               &
                          ,Rd_decay_elow(ipft),Rd_decay_ehigh(ipft),Rd_q10(ipft)           &
                          ,stomatal_slope(ipft),cuticular_cond(ipft)                       &
                          ,quantum_efficiency(ipft),curvpar_electron(ipft)                 &
                          ,qyield_psII(ipft),water_conductance(ipft)*yr_sec                &
                          ,leaf_width(ipft),root_respiration_factor(ipft)                  &
                          ,kplastic_vm0(ipft),kplastic_rd0(ipft)                           &
                          ,kplastic_sla(ipft),kplastic_LL(ipft)                            &
                          ,eplastic_vm0(ipft),eplastic_sla(ipft)
         end do
         !---------------------------------------------------------------------------------!
      end select
      close(unit=18,status='keep')
   end if
   !---------------------------------------------------------------------------------------!


   !----- Print allometric coefficients. --------------------------------------------------!
   if (print_zero_table) then
      open (unit=18,file=trim(allom_file),status='replace',action='write')
      write(unit=18,fmt='(55(1x,a))') '          PFT','     TROPICAL','        GRASS'      &
                                     ,'      CONIFER','     SAVANNAH','        LIANA'      &
                                     ,'    DDH_ALLOM','          RHO','         B1HT'      &
                                     ,'         B2HT','      HGT_REF','         B1BL'      &
                                     ,'         B2BL','   B1BS_SMALL','   B2BS_SMALL'      &
                                     ,'   B1BS_LARGE','   B2BS_LARGE','  D1DBH_SMALL'      &
                                     ,'  D2DBH_SMALL','  D1DBH_LARGE','  D2DBH_LARGE'      &
                                     ,'        L1DBH','        L2DBH','         B1CA'      &
                                     ,'         B2CA','        B1WAI','        B2WAI'      &
                                     ,'         B1SA','         B2SA','         B1RD'      &
                                     ,'         B2RD','         B1XS','         B1XB'      &
                                     ,'      HGT_MIN','      HGT_MAX','      MIN_DBH'      &
                                     ,'     DBH_CRIT','  DBH_BIGLEAF','   BDEAD_CRIT'      &
                                     ,'   BLEAF_CRIT','  BALIVE_CRIT','  BEVERY_CRIT'      &
                                     ,'    INIT_DENS','          SLA',' F_BSTOR_INIT'      &
                                     ,'            Q','          QSW','        QBARK'      &
                                     ,'        QRHOB','     d18O_REF','      B1_D18O'      &
                                     ,'      B2_D18O','      B1_EFRD','      B2_EFRD'      &
                                     ,'  INIT_LAIMAX'


      do ipft=1,n_pft
         write (unit=18,fmt='(9x,i5,6(13x,l1),47(1x,f13.6),1(1x,es13.6))')                 &
                        ipft,is_tropical(ipft),is_grass(ipft),is_conifer(ipft)             &
                       ,is_savannah(ipft),is_liana(ipft),ddh_allom(ipft),rho(ipft)         &
                       ,b1Ht(ipft),b2Ht(ipft),hgt_ref(ipft),b1Bl(ipft),b2Bl(ipft)          &
                       ,b1Bs_small(ipft),b2Bs_small(ipft),b1Bs_large(ipft)                 &
                       ,b2Bs_large(ipft),d1DBH_small(ipft),d2DBH_small(ipft)               &
                       ,d1DBH_large(ipft),d2DBH_large(ipft),l1DBH(ipft),l2DBH(ipft)        &
                       ,b1Ca(ipft),b2Ca(ipft),b1WAI(ipft),b2WAI(ipft),b1SA(ipft)           &
                       ,b2SA(ipft),b1Rd(ipft),b2Rd(ipft),b1Xs(ipft),b1Xb(ipft)             &
                       ,hgt_min(ipft),hgt_max(ipft),min_dbh(ipft),dbh_crit(ipft)           &
                       ,dbh_bigleaf(ipft),bdead_crit(ipft),bleaf_crit(ipft)                &
                       ,balive_crit(ipft),bevery_crit(ipft),init_density(ipft),sla(ipft)   &
                       ,f_bstorage_init(ipft),q(ipft),qsw(ipft),qbark(ipft),qrhob(ipft)    &
                       ,d18O_ref(ipft),b1d18O(ipft),b2d18O(ipft),b1Efrd(ipft),b2Efrd(ipft) &
                       ,init_laimax(ipft)
      end do
      close(unit=18,status='keep')
   end if


   !----- Print trait coefficients. -------------------------------------------------------!
   if (print_zero_table) then
      open (unit=19,file=trim(strat_file),status='replace',action='write')
      write(unit=19,fmt='( 98(1x,a))') '          PFT','     TROPICAL','        GRASS'     &
                                      ,'      CONIFER','     SAVANNAH','        LIANA'     &
                                      ,'       R_BANG','          RHO','          SLA'     &
                                      ,'          SRA','    ROOT_BETA','          VM0'     &
                                      ,'   F_DARKRESP','  F_GROW_RESP','     LEAF_TOR'     &
                                      ,'     ROOT_TOR','     BARK_TOR','  STORAGE_TOR'     &
                                      ,' FLABILE_LEAF',' FLABILE_STEM','        MORT0'     &
                                      ,'        MORT1','        MORT2','        MORT3'     &
                                      ,'    SEED_MORT',' TFALL_S_GTHT'                     &
                                      ,'   FIRE_S_MIN','   FIRE_S_MAX',' FIRE_S_INTER'     &
                                      ,' FIRE_S_SLOPE','     ST_FRACT','      R_FRACT'     &
                                      ,'       R_CV50','  NONLOC_DISP','    SEED_RAIN'     &
                                      ,'     EFF_HEAT','     EFF_EVAP','   EFF_TRANSP'     &
                                      ,'   LTRANS_VIS',' LREFLECT_VIS','   WTRANS_VIS'     &
                                      ,' WREFLECT_VIS','   LTRANS_NIR',' LREFLECT_NIR'     &
                                      ,'   WTRANS_NIR',' WREFLECT_NIR','   LEMISS_TIR'     &
                                      ,'   WEMISS_TIR','  ORIENT_FACT','    LSCAT_VIS'     &
                                      ,'   LBSCAT_VIS','    WSCAT_VIS','   WBSCAT_VIS'     &
                                      ,'    LSCAT_NIR','   LBSCAT_NIR','    WSCAT_NIR'     &
                                      ,'   WBSCAT_NIR','   LBSCAT_TIR','   WBSCAT_TIR'     &
                                      ,'         PHI1','         PHI2','       MU_BAR'     &
                                      ,'        CLEAF','        CSAPW','        CDEAD'     &
                                      ,'        CBARK','     C2N_LEAF','     C2N_STEM'     &
                                      ,'  C2N_STORAGE','  C2N_RECRUIT',' LEAF_SHED_RT'     &
                                      ,' LEAF_GROW_RT','  VESSEL_CURL',' LEAF_H2O_CAP'     &
                                      ,' WOOD_H2O_CAP',' LEAF_H2O_SAT',' WOOD_H2O_SAT'     &
                                      ,' LEAF_RWC_MIN',' WOOD_RWC_MIN','SMALL_RWC_MIN'     &
                                      ,' LEAF_PSI_MIN',' WOOD_PSI_MIN','SMALL_PSI_MIN'     &
                                      ,' LEAF_PSI_OSM',' WOOD_PSI_OSM',' LEAF_ELA_MOD'     &
                                      ,' WOOD_ELA_MOD',' LEAF_PSI_TLP',' WOOD_PSI_TLP'     &
                                      ,'    WOOD_KMAX','    WOOD_KEXP','   WOOD_PSI50'     &
                                      ,' STOMA_LAMBDA','   STOMA_BETA','  STOMA_PSI_B'     &
                                      ,'  STOMA_PSI_C',' HIGH_PSI_THR','  LOW_PSI_THR'

      do ipft=1,n_pft
         write (unit=19,fmt='(9x,i5,6(13x,l1),89(1x,f13.6),2(1x,i13))')                    &
                        ipft,is_tropical(ipft),is_grass(ipft),is_conifer(ipft)             &
                       ,is_savannah(ipft),is_liana(ipft),r_bang(ipft),rho(ipft),SLA(ipft)  &
                       ,SRA(ipft),root_beta(ipft),Vm0(ipft),dark_respiration_factor(ipft)  &
                       ,growth_resp_factor(ipft),leaf_turnover_rate(ipft)                  &
                       ,root_turnover_rate(ipft),bark_turnover_rate(ipft)                  &
                       ,storage_turnover_rate(ipft),f_labile_leaf(ipft)                    &
                       ,f_labile_stem(ipft),mort0(ipft),mort1(ipft),mort2(ipft)            &
                       ,mort3(ipft),seedling_mortality(ipft),treefall_s_ltht(ipft)         &
                       ,fire_s_min(ipft),fire_s_max(ipft)                                  &
                       ,fire_s_inter(ipft),fire_s_slope(ipft),st_fract(ipft),r_fract(ipft) &
                       ,r_cv50(ipft),nonlocal_dispersal(ipft),seed_rain(ipft)              &
                       ,effarea_heat,effarea_evap,effarea_transp(ipft)                     &
                       ,leaf_trans_vis(ipft),leaf_reflect_vis(ipft),wood_trans_vis(ipft)   &
                       ,wood_reflect_vis(ipft),leaf_trans_nir(ipft),leaf_reflect_nir(ipft) &
                       ,wood_trans_nir(ipft),wood_reflect_nir(ipft),leaf_emiss_tir(ipft)   &
                       ,wood_emiss_tir(ipft),orient_factor(ipft)                           &
                       ,leaf_scatter_vis(ipft),leaf_backscatter_vis(ipft)                  &
                       ,wood_scatter_vis(ipft),wood_backscatter_vis(ipft)                  &
                       ,leaf_scatter_nir(ipft),leaf_backscatter_nir(ipft)                  &
                       ,wood_scatter_nir(ipft),wood_backscatter_nir(ipft)                  &
                       ,leaf_backscatter_tir(ipft),wood_backscatter_tir(ipft)              &
                       ,phi1(ipft),phi2(ipft),mu_bar(ipft),cleaf(ipft),csapw(ipft)         &
                       ,cdead(ipft),cbark(ipft),c2n_leaf(ipft),c2n_stem(ipft),c2n_storage  &
                       ,c2n_recruit(ipft),leaf_shed_rate(ipft),leaf_grow_rate(ipft)        &
                       ,vessel_curl_factor(ipft),leaf_water_cap(ipft),wood_water_cap(ipft) &
                       ,leaf_water_sat(ipft),wood_water_sat(ipft),leaf_rwc_min(ipft)       &
                       ,wood_rwc_min(ipft),small_rwc_min(ipft),leaf_psi_min(ipft)          &
                       ,wood_psi_min(ipft),small_psi_min(ipft),leaf_psi_osmotic(ipft)      &
                       ,wood_psi_osmotic(ipft),leaf_elastic_mod(ipft)                      &
                       ,wood_elastic_mod(ipft),leaf_psi_tlp(ipft),wood_psi_tlp(ipft)       &
                       ,wood_Kmax(ipft),wood_Kexp(ipft),wood_psi50(ipft)                   &
                       ,stoma_lambda(ipft),stoma_beta(ipft),stoma_psi_b(ipft)              &
                       ,stoma_psi_c(ipft),high_psi_threshold(ipft),low_psi_threshold(ipft)
      end do
      close(unit=19,status='keep')
   end if



   return
end subroutine init_derived_params_after_xml
!==========================================================================================!
!==========================================================================================!

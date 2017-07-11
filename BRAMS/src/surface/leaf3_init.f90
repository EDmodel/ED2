!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!
subroutine sfcdata()
   use mem_grid
   use mem_leaf
   use leaf_coms

   implicit none
   !----- Local structures. ---------------------------------------------------------------!
   type soil_class
      !------------------------------------------------------------------------------------!
      !  Soil Characteristics (see Clapp & Hornberger, 1978; McCumber & Pielke, 1981;      !
      !                            Pielke, 1984; Tremback & Kessler, 1985)                 !
      !------------------------------------------------------------------------------------!
      real(kind=4) :: slpots     ! Soil moisture potential at saturation         [       m]
      real(kind=4) :: slmsts     ! Soil moisture at saturation                   [   m3/m3]
      real(kind=4) :: slbs       ! B exponent                                    [     n/d]
      real(kind=4) :: slcpd      ! Specific heat of dry soil                     [  J/m3/K]
      real(kind=4) :: soilcp     ! Dry soil capacity (at -3.1MPa)                [   m3/m3]
      real(kind=4) :: soilwp     ! Wilting point capacity (at -1.5MPa)           [   m3/m3]
      real(kind=4) :: slcons     ! hydraulic conductivity at saturation          [     m/s]
      real(kind=4) :: slcons0    ! Surface value for slcons                      [     m/s]
      real(kind=4) :: thcond0    ! Thermal conductance, term 3                   [        ]
      real(kind=4) :: thcond1    ! Thermal conductance, term 3                   [        ]
      real(kind=4) :: thcond2    ! Thermal conductance, term 3                   [        ]
      real(kind=4) :: thcond3    ! Thermal conductance, term 3                   [        ]
      real(kind=4) :: sfldcap    ! Soil field capacity                           [   m3/m3]
      real(kind=4) :: xsand      ! Percentage of sand                            [     ---]
      real(kind=4) :: xclay      ! Percentage of clay                            [     ---]
      real(kind=4) :: xsilt      ! Percentage of silt                            [     ---]
      real(kind=4) :: xrobulk    ! Bulk density                                  [     ---]
      real(kind=4) :: slden      ! "Dry" soil density (porosity)                 [   kg/m3]
      real(kind=4) :: soilld     ! Soil threshold for leaf drop                  [   m3/m3]
      real(kind=4) :: soilfr     ! Soil threshold for fire (unused)              [   m3/m3]
      real(kind=4) :: psiwp      ! Soil potential at wilting point               [       m]
      real(kind=4) :: psifc      ! Soil potential at field capacity              [       m]
      real(kind=4) :: psild      ! Soil potential threshold for leaf drop        [       m]
      real(kind=4) :: psifr      ! Soil potential threshold for fire (unused)    [       m]
   end type soil_class
   !.......................................................................................!
   type soilcol_class
      !------------------------------------------------------------------------------------!
      !     Define soil colour structure.                                                  !
      !------------------------------------------------------------------------------------!
      real(kind=4) :: alb_vis_dry ! Dry soil albedo       - Visible              [     ---]
      real(kind=4) :: alb_nir_dry ! Dry soil albedo       - Near Infrared        [     ---]
      real(kind=4) :: alb_vis_wet ! Saturated soil albedo - Visible              [     ---]
      real(kind=4) :: alb_nir_wet ! Saturated soil albedo - Near Infrared        [     ---]
   end type soilcol_class
   !.......................................................................................!
   type vegt_class
      !----- LEAF-3 biophysical parameters by land use class number. ----------------------!
      real(kind=4) :: albv_green
      real(kind=4) :: albv_brown
      real(kind=4) :: emisv
      real(kind=4) :: sr_max
      real(kind=4) :: tai_max
      real(kind=4) :: sai
      real(kind=4) :: veg_clump
      real(kind=4) :: veg_frac
      real(kind=4) :: veg_ht
      real(kind=4) :: rootdep
      real(kind=4) :: dead_frac
      real(kind=4) :: gsw_max
      real(kind=4) :: leaf_width
      real(kind=4) :: stom_side
      real(kind=4) :: phenology
      real(kind=4) :: clumping_factor
      real(kind=4) :: orient_factor
      real(kind=4) :: leaf_reflect_vis
      real(kind=4) :: leaf_reflect_nir
      real(kind=4) :: wood_reflect_vis
      real(kind=4) :: wood_reflect_nir
      real(kind=4) :: leaf_trans_vis
      real(kind=4) :: leaf_trans_nir
      real(kind=4) :: wood_trans_vis
      real(kind=4) :: wood_trans_nir
      real(kind=4) :: leaf_emiss_tir
      real(kind=4) :: wood_emiss_tir
      real(kind=4) :: pathway
      real(kind=4) :: leaf_supply
      real(kind=4) :: dr_gamma
      real(kind=4) :: gsw_0
      real(kind=4) :: gsw_m
      real(kind=4) :: quantum_yield
      real(kind=4) :: cn_l
      real(kind=4) :: f_lnr
      real(kind=4) :: fun_nitro
      real(kind=4) :: sla_0
      real(kind=4) :: sla_m
      real(kind=4) :: vm0_qten
      real(kind=4) :: vm0_dec
      real(kind=4) :: rd0_qten
      real(kind=4) :: rd0_dec
      real(kind=4) :: rr0_0
      real(kind=4) :: rr0_qten
      real(kind=4) :: rr0_dec
      real(kind=4) :: gr_factor
      real(kind=4) :: phys_low_temp
      real(kind=4) :: phys_high_temp
   end type vegt_class
   !----- Local variables. ----------------------------------------------------------------!
   type(soil_class)   , dimension(nstyp)            :: soilparms
   type(vegt_class)   , dimension(nvtyp+nvtyp_teb)  :: bioparms
   type(soilcol_class), dimension(nscol)            :: scolparms
   real               , dimension(nstyp)            :: xrobulk
   integer                                          :: k
   integer                                          :: nnn
   real                                             :: romin
   real                                             :: roorg
   real                                             :: slfcap
   real                                             :: refdepth
   real                                             :: tmin
   real                                             :: ratio
   real                                             :: xmin
   real                                             :: slz0
   real                                             :: ezg
   !------ Local constants. ---------------------------------------------------------------!
   real(kind=4)       , parameter                   :: f_nr      = 7.16
   real(kind=4)       , parameter                   :: alpha_r25 = 60.0
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Maximum TAI is based on average LAI by class after:                               !
   !                                                                                       !
   ! Asner, G. P., J. M. O. Scurlock, J. A. Hicke, 2003: Global synthesis of leaf area     !
   !    index observations: implications for ecological and remote sensing studies.        !
   !    Glob. Ecol. Biogeogr., 12, 191--205, doi: 10.1046/j.1466-822X.2003.00026.x         !
   !                                                                                       !
   ! SAI:TAI ratio came from the following publications:                                   !
   !                                                                                       !
   ! Breda, N. J. J., 2003: Ground-based measurements of leaf area index: a review of      !
   !    methods, instruments, and current controversies.  J. Exp. Bot. 54(392),            !
   !    2403--2417, doi: 10.1093/jxb/erg263.                                               !
   !                                                                                       !
   ! Gower, S. T., C. J. Kucharik, J. M. Norman, 1999: Direct and indirect estimation of   !
   !    leaf area index, f_APAR, and net primary production of terrestrial ecosystems.     !
   !    Remote Sens. Environ., 70(1), 29--51, doi: 10.1016/S0034-4257(99)00056-5           !
   !                                                                                       !
   ! Hirano, T., R. Hirata, Y. Fujinuma, N. Saigusa, S. Yamamoto, Y. Harazono, M. Takada,  !
   !    K. Inukai, G. Inoue, 2003: CO2 and water vapour exchange of a larch fores in       !
   !    northern Japan.  Tellus B, 55(2), 244--257, doi: 10.1034/j.1600-0889.2003.00063.x  !
   !                                                                                       !
   ! Olivas, P. C., S. F. Oberbauer, D. B. Clark, D. A. Clark, M. G. Ryan, J. J. O'Brien,  !
   !    H. Ordonez, 2013: Comparison of direct and indirect methods for assessing leaf     !
   !    area index across a tropical rain forest landscape.  Agric. For. Meteorol., 177,   !
   !    110--116, doi: 10.1016/j.agrformet.2013.04.010                                     !
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !    Declare the biophysics parameters.                                                 !
   ! ( 1st line)   albv_green    albv_brown          emisv          sr_max        tai_max  !
   ! ( 2nd line)          sai     veg_clump       veg_frac          veg_ht        rootdep  !
   ! ( 3rd line)    dead_frac       gsw_max     leaf_width       stom_side      phenology  !
   ! ( 4th line)     clumping        orient   leaf_ref_vis    leaf_ref_nir   wood_ref_vis  !
   ! ( 5th line) wood_ref_nir  leaf_tra_vis   leaf_tra_nir    wood_tra_vis   wood_tra_nir  !
   ! ( 6th line) leaf_emi_tir  wood_emi_tir        pathway     leaf_supply       dr_gamma  !
   ! ( 7th line)        gsw_0         gsw_m  quantum_yield            cn_l          f_lnr  !
   ! ( 8th line)    fun_nitro         sla_0          sla_m        vm0_qten        vm0_dec  !
   ! ( 9th line)     rd0_qten       rd0_dec          rr0_0        rr0_qten        rr0_dec  !
   ! (10th line)    gr_factor    phys_low_t    phys_high_t                                 !
   !---------------------------------------------------------------------------------------!
   bioparms = (/                                                                           &
   !-----  0. Ocean.  Not really used, and this is class 0.  Let's skip it. ---------------!
   !    vegt_class(      0.00,       0.00,            0.000,            0.0,           0.0 &
   !              ,       0.0,        0.0,            0.000,            0.0,           0.0 &
   !              ,       0.0,        0.0,            0.000,            0.0,           0.0 &
   !              ,       0.0,        0.0,            0.000,            0.0,           0.0 &
   !              ,       0.0,        0.0,            0.000,            0.0,           0.0 &
   !              ,       0.0,        0.0,            0.000,            0.0,           0.0 &
   !              ,       0.0,        0.0,            0.000,            0.0,           0.0 &
   !              ,       0.0,        0.0,            0.000,            0.0,           0.0 &
   !              ,       0.0,        0.0,            0.000,            0.0,           0.0 &
   !              ,       0.0,        0.0,            0.000)                               &
   !-----  1. Lakes, rivers, streams. -----------------------------------------------------!
       vegt_class(      0.00,       0.00,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000)                                &
   !-----  2. Ice cap/glacier. ------------------------------------------------------------!
      ,vegt_class(      0.00,       0.00,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000)                                &
   !-----  3. Desert, bare soil. ----------------------------------------------------------!
      ,vegt_class(      0.00,       0.00,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000,            0.0,           0.0  &
                 ,       0.0,        0.0,            0.000)                                &
   !-----  4. Evergreen needleleaf tree. --------------------------------------------------!
      ,vegt_class(      0.14,       0.24,            0.970,            5.4,           6.6  &
                 ,       1.1,        1.0,            0.800,           20.0,           1.5  &
                 ,       0.0,     0.0020,            0.050,            2.0,           0.0  &
                 ,      0.63,       0.01,            0.070,           0.35,          0.16  &
                 ,      0.39,       0.05,            0.100,          0.001,         0.001  &
                 ,      0.97,       0.90,            3.000,          0.016,        0.0145  &
                 ,     1000.,        7.2,            0.060,           35.0,         0.050  &
                 ,      0.72,       10.0,            1.250,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      4.714,           45.000)                                &
   !-----  5. Deciduous needleleaf tree. --------------------------------------------------!
      ,vegt_class(      0.14,       0.24,            0.950,            5.4,           5.9  &
                 ,       1.4,        1.0,            0.800,           22.0,           1.5  &
                 ,       0.0,     0.0020,            0.050,            2.0,           4.0  &
                 ,      0.78,       0.01,            0.070,           0.35,          0.16  &
                 ,      0.39,       0.05,            0.100,          0.001,         0.001  &
                 ,      0.97,       0.90,            3.000,          0.016,        0.0145  &
                 ,     1000.,        7.2,            0.060,           25.0,         0.080  &
                 ,      0.79,        24.,            3.000,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      4.714,           45.000)                                &
   !-----  6. Deciduous broadleaf tree. ---------------------------------------------------!
      ,vegt_class(      0.20,       0.24,            0.950,            6.2,           5.6  &
                 ,       1.1,        0.0,            0.800,           22.0,           1.5  &
                 ,       0.0,     0.0020,            0.100,            1.0,           4.0  &
                 ,     0.705,       0.25,            0.100,           0.45,          0.16  &
                 ,      0.39,       0.05,            0.250,          0.001,         0.001  &
                 ,      0.97,       0.90,            3.000,          0.016,        0.0145  &
                 ,    10000.,        7.2,            0.080,           25.0,         0.090  &
                 ,      0.66,        30.,            4.000,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      4.714,           45.000)                                &
   !-----  7. "Evergreen" broadleaf tree (but drought deciduous). -------------------------!
      ,vegt_class(      0.12,       0.18,            0.950,            4.1,           5.4  &
                 ,       0.6,        0.0,            0.900,           32.0,           2.5  &
                 ,       0.0,     0.0035,            0.200,            1.0,           4.0  &
                 ,     0.640,       0.10,            0.100,           0.46,          0.11  &
                 ,      0.25,       0.05,            0.230,          0.001,         0.001  &
                 ,      0.95,       0.90,            3.000,          0.016,        0.0145  &
                 ,    10000.,        9.0,            0.080,           30.0,         0.060  &
                 ,      0.83,       12.0,            1.500,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      8.000,           45.000)                                &
   !-----  8. Short grass. ----------------------------------------------------------------!
      ,vegt_class(      0.13,       0.30,            0.960,            5.1,           1.8  &
                 ,       0.1,        0.0,            0.750,            0.3,           0.7  &
                 ,       0.7,     0.0100,            0.100,            1.0,           4.0  &
                 ,      0.78,      -0.30,            0.100,           0.46,          0.16  &
                 ,      0.25,       0.05,            0.230,          0.028,         0.248  &
                 ,      0.96,       0.96,            4.000,          0.016,        0.0350  &
                 ,    10000.,        5.2,            0.055,           25.0,         0.090  &
                 ,      0.64,       30.0,            0.000,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      8.000,           45.000)                                &
   !-----  9. Tall grass. -----------------------------------------------------------------!
      ,vegt_class(      0.24,       0.43,            0.960,            5.1,           2.6  &
                 ,       0.1,        0.0,            0.800,            1.2,           1.0  &
                 ,       0.7,     0.0100,            0.100,            1.0,           4.0  &
                 ,      0.78,      -0.30,            0.100,           0.46,          0.16  &
                 ,      0.25,       0.05,            0.230,          0.028,         0.248  &
                 ,      0.96,       0.96,            4.000,          0.016,        0.0350  &
                 ,    10000.,        5.2,            0.055,           25.0,         0.090  &
                 ,      0.64,       30.0,            0.000,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      8.000,           45.000)                                &
   !----- 10. Semi-desert. ----------------------------------------------------------------!
      ,vegt_class(      0.24,       0.24,            0.960,            5.1,           1.7  &
                 ,       0.4,        1.0,            0.200,            0.7,           1.0  &
                 ,       0.0,     0.0020,            0.030,            1.0,           4.0  &
                 ,      0.78,        0.0,            0.100,           0.46,          0.16  &
                 ,      0.25,       0.05,            0.230,          0.028,         0.248  &
                 ,      0.96,       0.96,            4.000,          0.016,        0.0350  &
                 ,     3160.,        5.2,            0.055,           25.0,         0.090  &
                 ,      0.64,       30.0,            0.000,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      8.000,           45.000)                                &
   !----- 11. Tundra. ---------------------------------------------------------------------!
      ,vegt_class(      0.20,       0.24,            0.950,            5.1,           1.9  &
                 ,       0.2,        1.0,            0.600,            0.2,           1.0  &
                 ,       0.0,     0.0200,            0.030,            1.0,           4.0  &
                 ,      0.75,      -0.30,            0.110,           0.35,          0.31  &
                 ,      0.53,       0.05,            0.340,          0.120,         0.250  &
                 ,      0.96,       0.96,            3.000,          0.016,        0.0145  &
                 ,    10000.,        9.0,            0.060,           25.0,         0.090  &
                 ,      0.68,       30.0,            0.000,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      4.714,           45.000)                                &
   !----- 12. Evergreen shrub. ------------------------------------------------------------!
      ,vegt_class(      0.14,       0.24,            0.970,            5.1,           2.8  &
                 ,       0.7,        1.0,            0.700,            1.0,           1.0  &
                 ,       0.0,     0.0020,            0.050,            1.0,           0.0  &
                 ,      0.77,       0.01,            0.100,           0.46,          0.11  &
                 ,      0.25,       0.05,            0.230,          0.001,         0.001  &
                 ,      0.95,       0.95,            3.000,          0.016,        0.0145  &
                 ,    10000.,       9.00,            0.080,           30.0,          0.06  &
                 ,      0.62,       12.0,            0.000,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      8.000,           45.000)                                &
   !----- 13. Deciduous shrub. ------------------------------------------------------------!
      ,vegt_class(      0.20,       0.28,            0.970,            5.1,           2.8  &
                 ,       0.7,        1.0,            0.700,            1.0,           1.0  &
                 ,       0.0,     0.0020,            0.050,            1.0,           4.0  &
                 ,      0.74,       0.25,            0.100,           0.45,          0.16  &
                 ,      0.39,       0.05,            0.250,          0.001,         0.001  &
                 ,      0.97,       0.90,            3.000,          0.016,        0.0145  &
                 ,    10000.,       9.00,            0.080,           25.0,          0.09  &
                 ,      0.60,       30.0,            0.000,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      4.714,           45.000)                                &
   !----- 14. Mixed woodland. -------------------------------------------------------------!
      ,vegt_class(      0.16,       0.24,            0.960,            6.2,           6.1  &
                 ,       1.1,        0.5,            0.800,           22.0,           1.5  &
                 ,       0.0,     0.0020,            0.080,            1.0,           4.0  &
                 ,      0.70,       0.05,            0.100,           0.46,          0.11  &
                 ,      0.25,       0.05,            0.230,          0.001,         0.001  &
                 ,      0.95,       0.90,            3.000,          0.016,        0.0145  &
                 ,    10000.,       8.10,            0.080,           30.0,          0.07  &
                 ,      0.68,       20.0,            2.625,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      4.714,           45.000)                                &
   !----- 15. Crop/mixed farming, C3 grassland. -------------------------------------------!
      ,vegt_class(      0.22,       0.40,            0.950,            5.1,           3.8  &
                 ,       0.2,        0.0,            0.850,            1.0,           1.0  &
                 ,       0.0,     0.0100,            0.100,            1.0,           4.0  &
                 ,      0.78,      -0.30,            0.100,           0.46,          0.16  &
                 ,      0.25,       0.05,            0.230,          0.028,         0.248  &
                 ,      0.96,       0.96,            3.000,          0.016,        0.0145  &
                 ,    10000.,       9.00,            0.060,           25.0,          0.09  &
                 ,      0.61,       30.0,            0.000,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      4.714,           45.000)                                &
   !----- 16. Irrigated crop. -------------------------------------------------------------!
      ,vegt_class(      0.18,       0.40,            0.950,            5.1,           3.8  &
                 ,       0.2,        0.0,            0.800,            1.1,           1.0  &
                 ,       0.0,     0.0020,            0.100,            1.0,           0.0  &
                 ,      0.78,      -0.30,            0.100,           0.46,          0.16  &
                 ,      0.25,       0.05,            0.230,          0.028,         0.248  &
                 ,      0.96,       0.96,            3.000,          0.016,        0.0145  &
                 ,    10000.,       9.00,            0.060,           25.0,          0.09  &
                 ,      0.61,       30.0,            0.000,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      8.000,           45.000)                                &
   !----- 17. Bog or marsh. ---------------------------------------------------------------!
      ,vegt_class(      0.12,       0.43,            0.980,            5.1,           3.8  &
                 ,       0.2,        0.0,            0.800,            1.6,           1.0  &
                 ,       0.0,     0.0020,            0.200,            1.0,           0.0  &
                 ,      0.80,      -0.30,            0.100,           0.46,          0.11  &
                 ,      0.25,       0.05,            0.230,          0.001,         0.001  &
                 ,      0.95,       0.90,            3.000,          0.016,        0.0145  &
                 ,    10000.,        9.0,            0.080,           30.0,         0.060  &
                 ,      0.83,       12.0,            1.500,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      8.000,           45.000)                                &
   !----- 18. Wooded grassland. -----------------------------------------------------------!
      ,vegt_class(      0.13,       0.30,            0.960,            5.1,           4.0  &
                 ,       0.5,        0.0,            0.800,            7.0,           1.5  &
                 ,       0.0,     0.0100,            0.080,            1.0,           4.0  &
                 ,      0.70,      0.000,            0.100,           0.46,          0.11  &
                 ,      0.25,       0.05,            0.230,          0.001,         0.001  &
                 ,      0.95,       0.90,            3.000,          0.016,        0.0145  &
                 ,    10000.,        9.0,            0.080,           27.5,         0.075  &
                 ,      0.72,       21.0,            0.750,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      8.000,           45.000)                                &
   !----- 19. Urban and built up. ---------------------------------------------------------!
      ,vegt_class(      0.20,       0.36,            0.900,            5.1,           2.0  &
                 ,       0.3,        0.0,            0.500,            6.0,           0.8  &
                 ,       0.0,     0.0020,            0.050,            1.0,           0.0  &
                 ,      0.70,       0.00,            0.100,           0.46,          0.11  &
                 ,      0.25,       0.05,            0.230,          0.001,         0.001  &
                 ,      0.95,       0.90,            3.000,          0.016,        0.0145  &
                 ,    10000.,        9.0,            0.080,           30.0,         0.060  &
                 ,      0.83,       12.0,            1.500,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      8.000,           45.000)                                &
   !----- 20. Wetland evergreen broadleaf tree. -------------------------------------------!
      ,vegt_class(      0.17,       0.24,            0.950,            4.1,           7.1  &
                 ,       0.8,        0.0,            0.900,           32.0,           1.5  &
                 ,       0.0,     0.0020,            0.200,            1.0,           0.0  &
                 ,      0.67,       0.10,            0.100,           0.46,          0.11  &
                 ,      0.25,       0.05,            0.230,          0.001,         0.001  &
                 ,      0.95,       0.90,            3.000,          0.016,        0.0145  &
                 ,    10000.,        9.0,            0.080,           30.0,         0.060  &
                 ,      0.83,       12.0,            1.500,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      8.000,           45.000)                                &
   !----- 21. Very urban. -----------------------------------------------------------------!
      ,vegt_class(      0.16,       0.24,            0.960,            5.1,           1.2  &
                 ,       0.2,        1.0,            0.200,           20.0,           1.5  &
                 ,       0.0,     0.0020,            0.050,            1.0,           0.0  &
                 ,      0.70,       0.00,            0.100,           0.46,          0.11  &
                 ,      0.25,       0.05,            0.230,          0.001,         0.001  &
                 ,      0.95,       0.90,            3.000,          0.016,        0.0145  &
                 ,    10000.,        9.0,            0.080,           30.0,         0.060  &
                 ,      0.83,       12.0,            1.500,            2.4,           0.4  &
                 ,       2.4,        0.4,            0.672,            2.4,           0.4  &
                 ,     0.333,      8.000,           45.000)                                &
   /)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Assign the vector variables from the look-up table.  These variables will be used !
   ! by the biophysics solver in LEAF-3.                                                   !
   !---------------------------------------------------------------------------------------!
   do nnn = 1,(nvtyp+nvtyp_teb)
      !----- Variables obtained directly from the look up table above. --------------------!
      albv_green       (nnn) = bioparms(nnn)%albv_green
      albv_brown       (nnn) = bioparms(nnn)%albv_brown
      emisv            (nnn) = bioparms(nnn)%emisv
      sr_max           (nnn) = bioparms(nnn)%sr_max
      tai_max          (nnn) = bioparms(nnn)%tai_max
      sai              (nnn) = bioparms(nnn)%sai
      veg_clump        (nnn) = bioparms(nnn)%veg_clump
      veg_frac         (nnn) = bioparms(nnn)%veg_frac
      veg_ht           (nnn) = bioparms(nnn)%veg_ht
      dead_frac        (nnn) = bioparms(nnn)%dead_frac
      gsw_max          (nnn) = bioparms(nnn)%gsw_max
      leaf_width       (nnn) = bioparms(nnn)%leaf_width
      stom_side        (nnn) = bioparms(nnn)%stom_side
      phenology        (nnn) = nint(bioparms(nnn)%phenology)
      clumping_factor  (nnn) = dble(bioparms(nnn)%clumping_factor)
      orient_factor    (nnn) = dble(bioparms(nnn)%orient_factor)
      leaf_reflect_vis (nnn) = dble(bioparms(nnn)%leaf_reflect_vis)
      leaf_reflect_nir (nnn) = dble(bioparms(nnn)%leaf_reflect_nir)
      wood_reflect_vis (nnn) = dble(bioparms(nnn)%wood_reflect_vis)
      wood_reflect_nir (nnn) = dble(bioparms(nnn)%wood_reflect_nir)
      leaf_trans_vis   (nnn) = dble(bioparms(nnn)%leaf_trans_vis)
      leaf_trans_nir   (nnn) = dble(bioparms(nnn)%leaf_trans_nir)
      wood_trans_vis   (nnn) = dble(bioparms(nnn)%wood_trans_vis)
      wood_trans_nir   (nnn) = dble(bioparms(nnn)%wood_trans_nir)
      leaf_emiss_tir   (nnn) = dble(bioparms(nnn)%leaf_emiss_tir)
      wood_emiss_tir   (nnn) = dble(bioparms(nnn)%wood_emiss_tir)
      pathway          (nnn) = nint(bioparms(nnn)%pathway)
      leaf_supply      (nnn) = bioparms(nnn)%leaf_supply
      dr_gamma         (nnn) = bioparms(nnn)%dr_gamma
      gsw_0            (nnn) = bioparms(nnn)%gsw_0
      gsw_m            (nnn) = bioparms(nnn)%gsw_m
      quantum_yield    (nnn) = bioparms(nnn)%quantum_yield
      cn_l             (nnn) = bioparms(nnn)%cn_l
      f_lnr            (nnn) = bioparms(nnn)%f_lnr
      fun_nitro        (nnn) = bioparms(nnn)%fun_nitro
      sla_0            (nnn) = bioparms(nnn)%sla_0
      sla_m            (nnn) = bioparms(nnn)%sla_m
      vm0_qten         (nnn) = bioparms(nnn)%vm0_qten
      vm0_dec          (nnn) = bioparms(nnn)%vm0_dec
      rd0_qten         (nnn) = bioparms(nnn)%rd0_qten
      rd0_dec          (nnn) = bioparms(nnn)%rd0_dec
      rr0_0            (nnn) = bioparms(nnn)%rr0_0
      rr0_qten         (nnn) = bioparms(nnn)%rr0_qten
      rr0_dec          (nnn) = bioparms(nnn)%rr0_dec
      gr_factor        (nnn) = bioparms(nnn)%gr_factor
      phys_low_temp    (nnn) = bioparms(nnn)%phys_low_temp
      phys_high_temp   (nnn) = bioparms(nnn)%phys_high_temp
      !------------------------------------------------------------------------------------!



      !----- Maximum leaf area index. -----------------------------------------------------!
      glai_max(nnn)   = tai_max(nnn) - sai(nnn)
      !------------------------------------------------------------------------------------!

      !----- Not used... ------------------------------------------------------------------!
      root(1,nnn)  = 0.
      !------------------------------------------------------------------------------------!

      !----- Find the bottom layer that this patch can access water. ----------------------!
      kroot(nnn)   = nzg
      do k = nzg-1,1,-1
         if (slz(k+1) > -bioparms(nnn)%rootdep) kroot(nnn) = k
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Find the Vm0 coefficient, which is constant and will be divided by SLA to     !
      ! find the actual Vm0.                                                               !
      !------------------------------------------------------------------------------------!
      if (cn_l(nnn) == 0.0) then
         vm0_0(nnn) = 0.0
      else
         vm0_0(nnn) = f_lnr(nnn) * f_nr * alpha_r25 / cn_l(nnn)
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Forward scattering.                                                            !
      !------------------------------------------------------------------------------------!
      !----- Visible (PAR). ---------------------------------------------------------------!
      leaf_scatter_vis(nnn) = leaf_reflect_vis(nnn) + leaf_trans_vis(nnn)
      wood_scatter_vis(nnn) = wood_reflect_vis(nnn) + wood_trans_vis(nnn)
      !----- Near infrared (NIR). ---------------------------------------------------------!
      leaf_scatter_nir(nnn) = leaf_reflect_nir(nnn) + leaf_trans_nir(nnn)
      wood_scatter_nir(nnn) = wood_reflect_nir(nnn) + wood_trans_nir(nnn)
      !----- Thermal infrared (TIR). ------------------------------------------------------!
      leaf_scatter_tir(nnn) = 1.d0 - leaf_emiss_tir(nnn)
      wood_scatter_tir(nnn) = 1.d0 - wood_emiss_tir(nnn)
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !      Back-scattering coefficients following CLM.                                   !
      !------------------------------------------------------------------------------------!
      !----- Visible (PAR). ---------------------------------------------------------------!
      if (leaf_scatter_vis(nnn) == 0.d0) then
         leaf_backscatter_vis(nnn) = 5.d-1
      else
         leaf_backscatter_vis(nnn) = ( leaf_scatter_vis(nnn)                               &
                                      + 2.5d-1                                             &
                                      * ( leaf_reflect_vis(nnn) - leaf_trans_vis(nnn)   )  &
                                      * ( 1.d0 + orient_factor(nnn)) ** 2 )                &
                                    / ( 2.d0 * leaf_scatter_vis(nnn) )
      end if
      if (wood_scatter_vis(nnn) == 0.d0) then
         wood_backscatter_vis(nnn) = 5.d-1
      else
         wood_backscatter_vis(nnn) = ( wood_scatter_vis(nnn)                               &
                                      + 2.5d-1                                             &
                                      * ( wood_reflect_vis(nnn) - wood_trans_vis(nnn)   )  &
                                      * ( 1.d0 + orient_factor(nnn)) ** 2 )                &
                                    / ( 2.d0 * wood_scatter_vis(nnn) )
      end if
      !----- Near infrared (NIR). ---------------------------------------------------------!
      if (leaf_scatter_nir(nnn) == 0.d0) then
         leaf_backscatter_nir(nnn) = 5.d-1
      else
         leaf_backscatter_nir(nnn) = ( leaf_scatter_nir(nnn)                               &
                                      + 2.5d-1                                             &
                                      * ( leaf_reflect_nir(nnn) - leaf_trans_nir(nnn)   )  &
                                      * ( 1.d0 + orient_factor(nnn)) ** 2 )                &
                                    / ( 2.d0 * leaf_scatter_nir(nnn) )
      end if
      if (wood_scatter_nir(nnn) == 0.d0) then
         wood_backscatter_nir(nnn) = 5.d-1
      else
         wood_backscatter_nir(nnn) = ( wood_scatter_nir(nnn)                               &
                                      + 2.5d-1                                             &
                                      * ( wood_reflect_nir(nnn) - wood_trans_nir(nnn)   )  &
                                      * ( 1.d0 + orient_factor(nnn)) ** 2 )                &
                                    / ( 2.d0 * wood_scatter_nir(nnn) )
      end if
      !------------------------------------------------------------------------------------!
      !      Thermal infra-red (TIR): Here we use the same expression from CLM manual,     !
      ! further assuming that the transmittance is zero like Zhao and Qualls (2006) did,   !
      ! the backscattering coefficient becomes a function of the leaf orientation only.    !
      ! We don't have different orientation factor for wood (we could have), so we assume  !
      ! them to be the same as leaves.                                                     !
      !------------------------------------------------------------------------------------!
      leaf_backscatter_tir(nnn) = 5.d-1 + 1.25d-1 * (1 + orient_factor(nnn)) ** 2
      wood_backscatter_tir(nnn) = 5.d-1 + 1.25d-1 * (1 + orient_factor(nnn)) ** 2
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Light extinction coefficients.   These are found following CLM technical       !
      ! manual, and the values fall back to ED-2.0 defaults when orient_factor is zero.    !
      !---------------------------------------------------------------------------------------!
      phi1(nnn) = 5.00d-1 - orient_factor(nnn) * ( 6.33d-1 + 3.3d-1 * orient_factor(nnn) )
      phi2(nnn) = 8.77d-1 * (1.d0 - 2.d0 * phi1(nnn))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the average inverse diffuse optical depth per unit leaf and stem area.    !
      ! We follow CLM technical manual, equation 3.4 only when the orientation factor is   !
      ! non-zero.   Otherwise, we make it 1.d0, which is the limit of that equation when   !
      ! phi2 approaches zero.                                                              !
      !------------------------------------------------------------------------------------!
      if (orient_factor(nnn) == 0.d0) then
         mu_bar(nnn) = 1.d0
      else
         mu_bar(nnn) = ( 1.d0 - phi1(nnn) * log(1.d0 + phi2(nnn)/phi1(nnn)) / phi2(nnn) )  &
                     / phi2(nnn)
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   ! Define the soil parameters.                                                           !
   ! (1st line)          slpots        slmsts          slbs      slcpd        soilcp       !
   ! (2nd line)          soilwp        slcons       slcons0    thcond0       thcond1       !
   ! (3rd line)         thcond2       thcond3       sfldcap      xsand         xclay       !
   ! (4th line)           xsilt       xrobulk      slden        soilld        soilfr       !
   ! (5th line)         slpotwp       slpotfc    slpotld       slpotfr                     !
   !---------------------------------------------------------------------------------------!
   soilparms = (/                                                                          &
      !----- 1. Sand. ---------------------------------------------------------------------!
       soil_class( -0.049831046,     0.373250,     3.295000,  1342809.,  0.026183447       &
                 ,  0.032636854,  2.446421e-5,  0.000500000, 0.9546011,    0.5333047       &
                 ,    0.6626306,   -0.4678112,  0.132130936,     0.920,        0.030       &
                 ,        0.050,        1200.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,        0.000,     0.000              )      &
      !----- 2. Loamy sand. ---------------------------------------------------------------!
      ,soil_class( -0.067406224,     0.385630,     3.794500,  1326165.,  0.041560499       &
                 ,  0.050323046,  1.776770e-5,  0.000600000, 0.9279457,    0.5333047       &
                 ,    0.6860126,   -0.4678112,  0.155181959,     0.825,        0.060       &
                 ,        0.115,        1250.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,     0.000,        0.000              )      &
      !----- 3. Sandy loam. ---------------------------------------------------------------!
      ,soil_class( -0.114261521,     0.407210,     4.629000,  1295982.,  0.073495043       &
                 ,  0.085973722,  1.022660e-5,  0.000769000, 0.8826064,    0.5333047       &
                 ,    0.7257838,   -0.4678112,  0.194037750,     0.660,        0.110       &
                 ,        0.230,        1300.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,     0.000,        0.000              )      &
      !----- 4. Silt loam. ----------------------------------------------------------------!
      ,soil_class( -0.566500112,     0.470680,     5.552000,  1191975.,  0.150665475       &
                 ,  0.171711257,  2.501101e-6,  0.000010600, 0.7666418,    0.5333047       &
                 ,    0.8275072,   -0.4678112,  0.273082063,     0.200,        0.160       &
                 ,        0.640,        1400.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,     0.000,        0.000              )      &
      !----- 5. Loam. ---------------------------------------------------------------------!
      ,soil_class( -0.260075834,     0.440490,     5.646000,  1245546.,  0.125192234       &
                 ,  0.142369513,  4.532431e-6,  0.002200000, 0.8168244,    0.5333047       &
                 ,    0.7834874,   -0.4678112,  0.246915025,     0.410,        0.170       &
                 ,        0.420,        1350.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,     0.000,        0.000              )      &
      !----- 6. Sandy clay loam. ----------------------------------------------------------!
      ,soil_class( -0.116869181,     0.411230,     7.162000,  1304598.,  0.136417267       &
                 ,  0.150969505,  6.593731e-6,  0.001500000, 0.8544779,    0.5333047       &
                 ,    0.7504579,   -0.4678112,  0.249629687,     0.590,        0.270       &
                 ,        0.140,        1350.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,        0.000,     0.000              )      &
      !----- 7. Silty clay loam. ----------------------------------------------------------!
      ,soil_class( -0.627769194,     0.478220,     8.408000,  1193778.,  0.228171947       &
                 ,  0.248747504,  1.435262e-6,  0.000107000, 0.7330059,    0.5333047       &
                 ,    0.8570124,   -0.4678112,  0.333825332,     0.100,        0.340       &
                 ,        0.560,        1500.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,     0.000,        0.000              )      &
      !----- 8. Clayey loam. --------------------------------------------------------------!
      ,soil_class( -0.281968114,     0.446980,     8.342000,  1249582.,  0.192624431       &
                 ,  0.210137962,  2.717260e-6,  0.002200000, 0.7847168,    0.5333047       &
                 ,    0.8116520,   -0.4678112,  0.301335491,     0.320,        0.340       &
                 ,        0.340,        1450.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,        0.000,     0.000              )      &
      !----- 9. Sandy clay. ---------------------------------------------------------------!
      ,soil_class( -0.121283019,     0.415620,     9.538000,  1311396.,  0.182198910       &
                 ,  0.196607427,  4.314507e-6,  0.000002167, 0.8273339,    0.5333047       &
                 ,    0.7742686,   -0.4678112,  0.286363001,     0.520,        0.420       &
                 ,        0.060,        1450.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,        0.000,     0.000              )      &
      !----- 10. Silty clay. --------------------------------------------------------------!
      ,soil_class( -0.601312179,     0.479090,    10.461000,  1203168.,  0.263228486       &
                 ,  0.282143846,  1.055191e-6,  0.000001033, 0.7164724,    0.5333047       &
                 ,    0.8715154,   -0.4678112,  0.360319788,     0.060,        0.470       &
                 ,        0.470,        1650.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,     0.000,        0.000              )      &
      !----- 11. Clay. --------------------------------------------------------------------!
      ,soil_class( -0.299226464,     0.454400,    12.460000,  1259466.,  0.259868987       &
                 ,  0.275459057,  1.307770e-6,  0.000001283, 0.7406805,    0.5333047       &
                 ,    0.8502802,   -0.4678112,  0.353255209,     0.200,        0.600       &
                 ,        0.200,        1700.,     1600.,        0.000,        0.000       &
                 ,        0.000,        0.000,     0.000,        0.000              )      &
      !----- 12. Peat. --------------------------------------------------------------------!
      ,soil_class( -0.534564359,     0.469200,     6.180000,   874000.,  0.167047523       &
                 ,  0.187868805,  2.357930e-6,  0.000008000, 0.7644011,    0.5333047       &
                 ,    0.8294728,   -0.4678112,  0.285709966,    0.2000,       0.2000       &
                 ,       0.6000,         500.,         300.,     0.000,        0.000       &
                 ,        0.000,        0.000,        0.000,     0.000              )      &
      !----- 13. Bedrock. -----------------------------------------------------------------!
      ,soil_class(    0.0000000,     0.000000,     0.000000,  2130000.,  0.000000000       &
                 ,  0.000000000,  0.000000e+0,  0.000000000, 1.3917897,    0.5333047       &
                 ,    0.2791318,   -0.4678112,  0.000000001,        0.0000        &
                 ,       0.0000,       0.0000,           0.,        0.,        0.000       &
                 ,        0.000,        0.000,        0.000,     0.000,        0.000)      &
      !----- 14. Silt. --------------------------------------------------------------------!
      ,soil_class( -1.047128548,     0.492500,     3.862500,  1143842.,  0.112299080       &
                 ,  0.135518820,  2.046592e-6,  0.000010600, 0.7425839,    0.5333047       &
                 ,    0.8486106,   -0.4678112,  0.245247642,     0.075,        0.050       &
                 ,        0.875,        1400.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,        0.000,     0.000              )      &
      !----- 15. Heavy clay. --------------------------------------------------------------!
      ,soil_class( -0.322106879,     0.461200,    15.630000,  1264547.,  0.296806035       &
                 ,  0.310916364,  7.286705e-7,  0.000001283, 0.7057374,    0.5333047       &
                 ,    0.8809321,   -0.4678112,  0.382110712,     0.100,        0.800       &
                 ,        0.100,        1700.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,        0.000,     0.000              )      &
      !----- 16. Clayey sand. -------------------------------------------------------------!
      ,soil_class( -0.176502150,     0.432325,    11.230000,  1292163.,  0.221886929       &
                 ,  0.236704039,  2.426785e-6,  0.000001283, 0.7859325,    0.5333047       &
                 ,    0.8105855,   -0.4678112,  0.320146708,     0.375,        0.525       &
                 ,        0.100,        1700.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,        0.000,     0.000              )      &
      !----- 17. Clayey silt. -------------------------------------------------------------!
      ,soil_class( -0.438278332,     0.467825,    11.305000,  1228490.,  0.261376708       &
                 ,  0.278711303,  1.174982e-6,  0.000001283, 0.7281197,    0.5333047       &
                 ,    0.8612985,   -0.4678112,  0.357014719,     0.125,        0.525       &
                 ,        0.350,        1700.,        1600.,     0.000,        0.000       &
                 ,        0.000,        0.000,        0.000,     0.000              )      &
   /)
   !---------------------------------------------------------------------------------------!


   !------ Thermal conductivity in J/msK. -------------------------------------------------!
   cka   = 0.418684 * 0.0615
   ckw   = 0.418684 * 1.45
   romin = 2655.0
   roorg = 1300.0


   !------ Set the top soil depth to be zero. ---------------------------------------------!
   slz(nzg+1) = 0.
   slfcap     = -10. / 3.

   !------ The reference depth varies depending on the surface model. ---------------------!
   select case (isfcl)
   case (2)
      refdepth   = -2.0
   case default
      refdepth   = -0.5
   end select

   do k = 1,nzg
      slzt   (k) = .5 * (slz(k) + slz(k+1))
   end do

   !----- Find the exponential increase factor. -------------------------------------------!
   ezg    = log(slz(1)/slz(nzg)) / log(real(nzg))
   slz0   = slz(1) * (real(nzg+1)/real(nzg))**ezg

   slzt(0) = .5 * (slz0 + slz(1))



   !----- Determine the top layer to consider for heterotrophic respiration. --------------!
   k_rh_loop: do k_hetresp=nzg-1,1,-1
     if (slz(k_hetresp) < depth_hetresp) exit k_rh_loop
   end do k_rh_loop
   k_hetresp = k_hetresp + 1
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Populate the parameters for all soil types.                                       !
   !---------------------------------------------------------------------------------------!
   do nnn = 1,nstyp
      slcons0(nnn) = soilparms(nnn)%slcons0
      if (nnn /= 13) then
         fhydraul(nnn) = log (soilparms(nnn)%slcons / soilparms(nnn)%slcons0) / refdepth
      else
         fhydraul(nnn) = 0.
      end if


      select case (isfcl)
      case (2)
         !----- TOPMODEL form - large at surface, and exponential decrease with depth. ----!
         do k = 1,nzg
            slcons1(k,nnn) = soilparms(nnn)%slcons0 * exp(slz(k) * fhydraul(nnn))
         end do
         slcons1(0,nnn) = soilparms(nnn)%slcons0 * exp(slz0 * fhydraul(nnn))
      case default
         select case (ipercol)
         case (0,1)
            !----- ORIGINAL form - const with depth. --------------------------------------!
            do k=0,nzg
               slcons1(k,nnn) = soilparms(nnn)%slcons
            end do

         case (2)
            !------------------------------------------------------------------------------!
            !    TOPMODEL form, similar to CLM.  Here we use the same definition of slcons !
            ! from Cosby et al. (1984) because it has a stronger spread and it accounts    !
            ! for sand and clay contents.                                                  !
            !------------------------------------------------------------------------------!
            do k=0,nzg
               slcons1(k,nnn) = soilparms(nnn)%slcons * exp ( - slzt(k) / refdepth)
            end do
         end select
      end select
      !------ Copy the other parameters to the vectors. -----------------------------------!
      slpots   (nnn) = soilparms(nnn)%slpots
      slmsts   (nnn) = soilparms(nnn)%slmsts
      slbs     (nnn) = soilparms(nnn)%slbs
      slcons   (nnn) = soilparms(nnn)%slcons
      slcons00 (nnn) = soilparms(nnn)%slcons0
      slcpd    (nnn) = soilparms(nnn)%slcpd
      slden    (nnn) = soilparms(nnn)%slden
      sfldcap  (nnn) = soilparms(nnn)%sfldcap
      soilcp   (nnn) = soilparms(nnn)%soilcp
      soilwp   (nnn) = soilparms(nnn)%soilwp
      thcond0  (nnn) = soilparms(nnn)%thcond0
      thcond1  (nnn) = soilparms(nnn)%thcond1
      thcond2  (nnn) = soilparms(nnn)%thcond2
      thcond3  (nnn) = soilparms(nnn)%thcond3
      xsand    (nnn) = soilparms(nnn)%xsand
      xclay    (nnn) = soilparms(nnn)%xclay
      xsilt    (nnn) = soilparms(nnn)%xsilt
      xrobulk  (nnn) = soilparms(nnn)%xrobulk
      !------------------------------------------------------------------------------------!

      !----- Find potential at the wilting point. -----------------------------------------!
      if (nnn /= 13) then
         psild    (nnn) = -76.47872
         psiwp    (nnn) = soilparms(nnn)%slpots                                            &
                           / (soilparms(nnn)%soilwp / soilparms(nnn)%slmsts)               &
                           ** soilparms(nnn)%slbs
         psifc    (nnn) = soilparms(nnn)%slpots                                            &
                           / (soilparms(nnn)%sfldcap / soilparms(nnn)%slmsts)              &
                           ** soilparms(nnn)%slbs
      else
         psild    (nnn) = 0.0
         psiwp    (nnn) = 0.0
         psifc    (nnn) = 0.0
      end if
      !------------------------------------------------------------------------------------!

      !----- Define the emmisivity and slfc. ----------------------------------------------!
      emisg(nnn) = .98
      if (nnn /= 13) then
         slfc(nnn) = slmsts(nnn) * (slfcap / slpots(nnn)) ** (-1. / slbs(nnn))
      else
         slfc(nnn) = 0.0
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Fill in the albedo information regarding the soil colour classes.                 !
   !---------------------------------------------------------------------------------------!
   !                    |    Dry soil   |   Saturated   |                                  !
   !   Soil class       |---------------+---------------|                                  !
   !                    |   VIS |   NIR |   VIS |   NIR |                                  !
   !---------------------------------------------------------------------------------------!
   scolparms = (/                                            & !
       soilcol_class   (    0.36,   0.61,   0.25,   0.50  )  & ! 01 - Brightest
      ,soilcol_class   (    0.34,   0.57,   0.23,   0.46  )  & ! 02
      ,soilcol_class   (    0.32,   0.53,   0.21,   0.42  )  & ! 03
      ,soilcol_class   (    0.31,   0.51,   0.20,   0.40  )  & ! 04
      ,soilcol_class   (    0.30,   0.49,   0.19,   0.38  )  & ! 05
      ,soilcol_class   (    0.29,   0.48,   0.18,   0.36  )  & ! 06
      ,soilcol_class   (    0.28,   0.45,   0.17,   0.34  )  & ! 07
      ,soilcol_class   (    0.27,   0.43,   0.16,   0.32  )  & ! 08
      ,soilcol_class   (    0.26,   0.41,   0.15,   0.30  )  & ! 09
      ,soilcol_class   (    0.25,   0.39,   0.14,   0.28  )  & ! 10
      ,soilcol_class   (    0.24,   0.37,   0.13,   0.26  )  & ! 11
      ,soilcol_class   (    0.23,   0.35,   0.12,   0.24  )  & ! 12
      ,soilcol_class   (    0.22,   0.33,   0.11,   0.22  )  & ! 13
      ,soilcol_class   (    0.20,   0.31,   0.10,   0.20  )  & ! 14
      ,soilcol_class   (    0.18,   0.29,   0.09,   0.18  )  & ! 15
      ,soilcol_class   (    0.16,   0.27,   0.08,   0.16  )  & ! 16
      ,soilcol_class   (    0.14,   0.25,   0.07,   0.14  )  & ! 17
      ,soilcol_class   (    0.12,   0.23,   0.06,   0.12  )  & ! 18
      ,soilcol_class   (    0.10,   0.21,   0.05,   0.10  )  & ! 19
      ,soilcol_class   (    0.08,   0.16,   0.04,   0.08  )  & ! 20 - Darkest
      ,soilcol_class   (    0.00,   0.00,   0.00,   0.00  )  & ! 21 - ED-2.1, unused
   /)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Fill in the vectors.                                                               !
   !---------------------------------------------------------------------------------------!
   do nnn = 1,nscol
      alb_vis_dry(nnn) = scolparms(nnn)%alb_vis_dry
      alb_nir_dry(nnn) = scolparms(nnn)%alb_nir_dry
      alb_vis_wet(nnn) = scolparms(nnn)%alb_vis_wet
      alb_nir_wet(nnn) = scolparms(nnn)%alb_nir_wet
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine sfcdata
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!      Fill the snow_mass and snow_depth arrays with a default value of 0.  This default   !
! is used when snowcover information is not read from varfiles.                            !
!------------------------------------------------------------------------------------------!
subroutine snowinit(n2,n3,snow_mass,snow_depth)
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer                  , intent(in)    :: n2
   integer                  , intent(in)    :: n3
   real   , dimension(n2,n3), intent(inout) :: snow_mass
   real   , dimension(n2,n3), intent(inout) :: snow_depth
   !------ Local variables. ---------------------------------------------------------------!
   integer                                  :: i
   integer                                  :: j
   !---------------------------------------------------------------------------------------!

   do j = 1,n3
      do i = 1,n2
         snow_mass(i,j) = 0.
         snow_depth(i,j) = 0.
      end do
   end do

   return
end subroutine snowinit
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!      This subroutine fills the arrays PATCH_AREA, leaf_class, and SOIL_TEXT horizontally !
! homogeneously with default values that are defined in the RAMSIN namelist file.  These   !
! fields comprise the land/sea surface data types that are normally available on standard  !
! RAMS datasets.  The default values assigned here may be overridden by:                   !
! (1) interpolation from coarser grids,                                                    !
! (2) specifying new hardcoded values in subroutine sfcinit_user in the file ruser.f90,    !
! (3) reading data from the standard RAMS datasets.                                        !
!------------------------------------------------------------------------------------------!
subroutine sfcinit_file(n2,n3,mzg,npat,ifm,patch_area,leaf_class,soil_color,soil_text)
   use mem_leaf
   use rconstants
   use leaf_coms , only : min_patch_area

   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer                           , intent(in)    :: n2
   integer                           , intent(in)    :: n3
   integer                           , intent(in)    :: mzg
   integer                           , intent(in)    :: npat
   integer                           , intent(in)    :: ifm
   real   , dimension(n2,n3,npat)    , intent(inout) :: soil_color
   real   , dimension(mzg,n2,n3,npat), intent(inout) :: soil_text
   real   , dimension(n2,n3,npat)    , intent(inout) :: patch_area
   real   , dimension(n2,n3,npat)    , intent(inout) :: leaf_class
   !------ Local variables. ---------------------------------------------------------------!
   integer                                           :: i
   integer                                           :: j
   integer                                           :: k
   integer                                           :: ipat
   !---------------------------------------------------------------------------------------!



   !----- Re-define the percentage of land cover, so it's always bounded. -----------------!
   pctlcon = max(0.,min(pctlcon, 1. - min_patch_area))


   !---------------------------------------------------------------------------------------!
   !     Create two patches, one with water, the other with land.                          !
   !---------------------------------------------------------------------------------------!
   do j = 1,n3
      do i = 1,n2
         !----- Water patch. --------------------------------------------------------------!
         patch_area(i,j,1) = 1. - pctlcon
         leaf_class(i,j,1) = 0.
         do k = 1,mzg
            soil_text(k,i,j,1) = 0.               ! patch 1
         end do
         soil_color(i,j,1) = 0.

         !----- Land patch. ---------------------------------------------------------------!
         patch_area(i,j,2) = pctlcon
         leaf_class(i,j,2) = float(nvgcon)
         do k = 1,mzg
            soil_text(k,i,j,2) = float(nslcon)
         end do
         soil_color(i,j,2) = isoilcol
      end do
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Leave any other patch that may exist empty.                                       !
   !---------------------------------------------------------------------------------------!
   do ipat = 3,npat
      do j = 1,n3
         do i = 1,n2
            patch_area(i,j,ipat) = 0.
            leaf_class(i,j,ipat) = leaf_class(i,j,2)

            do k = 1,mzg
               soil_text(k,i,j,ipat) = soil_text(k,i,j,2)
            end do

            soil_color(i,j,ipat) = soil_color(i,j,2)

         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!
   return
end subroutine sfcinit_file
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine fills the primary LEAF3 arrays for which standard RAMS data files   !
! do not exist with initial default values.  Many of the initial values are horizontally   !
! homogeneous, although some depend on atmospheric conditions.  The default values         !
! assigned here may be overridden by:                                                      !
! 1. Specification from coarser grids;                                                     !
! 2. Specifying new values in subroutine sfcinit_nofile_user in the file ruser.f90.        !
!------------------------------------------------------------------------------------------!
subroutine sfcinit_nofile(n1,n2,n3,mzg,mzs,npat,ifm,theta,pi0,pp,rv,co2p,seatp,seatf       &
                         ,soil_water,soil_energy,psibar_10d,soil_color,soil_text           &
                         ,sfcwater_mass,sfcwater_energy,sfcwater_depth,ustar,tstar,rstar   &
                         ,cstar,zeta,ribulk,veg_fracarea,veg_agb,veg_lai,veg_tai,veg_rough &
                         ,veg_height,veg_displace,veg_albedo,patch_area,patch_rough        &
                         ,patch_wetind,leaf_class,soil_rough,sfcwater_nlev,stom_condct     &
                         ,ground_rsat,ground_rvap,ground_temp,ground_fliq,veg_water        &
                         ,veg_hcap,veg_energy,can_prss,can_theiv,can_vpdef,can_theta       &
                         ,can_rvap,can_co2,hflxac_out,wflxac_out,qwflxac_out,eflxac_out    &
                         ,cflxac_out,hflxgc_out,wflxgc_out,qwflxgc_out,hflxvc_out          &
                         ,wflxvc_out,qwflxvc_out,transp_out,qtransp_out,intercepted_out    &
                         ,qintercepted_out,wshed_out,qwshed_out,throughfall_out            &
                         ,qthroughfall_out,runoff_out,qrunoff_out,drainage_out             &
                         ,qdrainage_out,gpp_out,plresp_out,resphet_out,growresp            &
                         ,veg_ndvip,veg_ndvic,veg_ndvif,snow_mass,snow_depth,rshort_gnd    &
                         ,rlong_gnd,cosz,rlongup,albedt,rvv,prsv,piv,vt2da,vt2db           &
                         ,glat,glon,zot,flpw,rtgt)
   use mem_grid
   use mem_leaf
   use leaf_coms
   use io_params
   use rconstants
   use therm_lib  , only : reducedpress & ! function
                         , thetaeiv     & ! function
                         , vpdefil      & ! function
                         , press2exner  & ! function
                         , exner2press  & ! function
                         , extheta2temp & ! function
                         , cmtl2uext    & ! function
                         , tl2uint      & ! function
                         , tq2enthalpy  ! ! function

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                        , intent(in)    :: n1
   integer                        , intent(in)    :: n2
   integer                        , intent(in)    :: n3
   integer                        , intent(in)    :: mzg
   integer                        , intent(in)    :: mzs
   integer                        , intent(in)    :: npat
   integer                        , intent(in)    :: ifm
   real, dimension(n1,n2,n3)      , intent(in)    :: theta
   real, dimension(n1,n2,n3)      , intent(in)    :: pi0
   real, dimension(n1,n2,n3)      , intent(in)    :: pp
   real, dimension(n1,n2,n3)      , intent(in)    :: rv
   real, dimension(n1,n2,n3)      , intent(in)    :: co2p
   real, dimension(   n2,n3)      , intent(in)    :: seatp
   real, dimension(   n2,n3)      , intent(in)    :: seatf
   real, dimension(   n2,n3)      , intent(in)    :: snow_mass
   real, dimension(   n2,n3)      , intent(in)    :: snow_depth
   real, dimension(   n2,n3)      , intent(in)    :: flpw
   real, dimension(   n2,n3)      , intent(in)    :: rtgt
   real, dimension(   n2,n3)      , intent(in)    :: cosz
   real, dimension(mzg,n2,n3,npat), intent(inout) :: soil_water
   real, dimension(mzg,n2,n3,npat), intent(inout) :: soil_energy
   real, dimension(    n2,n3,npat), intent(inout) :: psibar_10d
   real, dimension(    n2,n3,npat), intent(inout) :: soil_color
   real, dimension(mzg,n2,n3,npat), intent(inout) :: soil_text
   real, dimension(mzs,n2,n3,npat), intent(inout) :: sfcwater_mass
   real, dimension(mzs,n2,n3,npat), intent(inout) :: sfcwater_energy
   real, dimension(mzs,n2,n3,npat), intent(inout) :: sfcwater_depth
   real, dimension(    n2,n3,npat), intent(inout) :: ustar
   real, dimension(    n2,n3,npat), intent(inout) :: tstar
   real, dimension(    n2,n3,npat), intent(inout) :: rstar
   real, dimension(    n2,n3,npat), intent(inout) :: cstar
   real, dimension(    n2,n3,npat), intent(inout) :: zeta
   real, dimension(    n2,n3,npat), intent(inout) :: ribulk
   real, dimension(    n2,n3,npat), intent(inout) :: veg_fracarea
   real, dimension(    n2,n3,npat), intent(inout) :: veg_agb
   real, dimension(    n2,n3,npat), intent(inout) :: veg_lai
   real, dimension(    n2,n3,npat), intent(inout) :: veg_tai
   real, dimension(    n2,n3,npat), intent(inout) :: veg_rough
   real, dimension(    n2,n3,npat), intent(inout) :: veg_height
   real, dimension(    n2,n3,npat), intent(inout) :: veg_displace
   real, dimension(    n2,n3,npat), intent(inout) :: veg_albedo
   real, dimension(    n2,n3,npat), intent(inout) :: patch_area
   real, dimension(    n2,n3,npat), intent(inout) :: patch_rough
   real, dimension(    n2,n3,npat), intent(inout) :: patch_wetind
   real, dimension(    n2,n3,npat), intent(inout) :: leaf_class
   real, dimension(    n2,n3,npat), intent(inout) :: soil_rough
   real, dimension(    n2,n3,npat), intent(inout) :: sfcwater_nlev
   real, dimension(    n2,n3,npat), intent(inout) :: stom_condct
   real, dimension(    n2,n3,npat), intent(inout) :: ground_rsat
   real, dimension(    n2,n3,npat), intent(inout) :: ground_rvap
   real, dimension(    n2,n3,npat), intent(inout) :: ground_temp
   real, dimension(    n2,n3,npat), intent(inout) :: ground_fliq
   real, dimension(    n2,n3,npat), intent(inout) :: veg_water
   real, dimension(    n2,n3,npat), intent(inout) :: veg_energy
   real, dimension(    n2,n3,npat), intent(inout) :: veg_hcap
   real, dimension(    n2,n3,npat), intent(inout) :: can_prss
   real, dimension(    n2,n3,npat), intent(inout) :: can_theiv
   real, dimension(    n2,n3,npat), intent(inout) :: can_vpdef
   real, dimension(    n2,n3,npat), intent(inout) :: can_theta
   real, dimension(    n2,n3,npat), intent(inout) :: can_rvap
   real, dimension(    n2,n3,npat), intent(inout) :: can_co2
   real, dimension(    n2,n3,npat), intent(inout) :: hflxac_out
   real, dimension(    n2,n3,npat), intent(inout) :: wflxac_out
   real, dimension(    n2,n3,npat), intent(inout) :: qwflxac_out
   real, dimension(    n2,n3,npat), intent(inout) :: eflxac_out
   real, dimension(    n2,n3,npat), intent(inout) :: cflxac_out
   real, dimension(    n2,n3,npat), intent(inout) :: hflxgc_out
   real, dimension(    n2,n3,npat), intent(inout) :: wflxgc_out
   real, dimension(    n2,n3,npat), intent(inout) :: qwflxgc_out
   real, dimension(    n2,n3,npat), intent(inout) :: hflxvc_out
   real, dimension(    n2,n3,npat), intent(inout) :: wflxvc_out
   real, dimension(    n2,n3,npat), intent(inout) :: qwflxvc_out
   real, dimension(    n2,n3,npat), intent(inout) :: transp_out
   real, dimension(    n2,n3,npat), intent(inout) :: qtransp_out
   real, dimension(    n2,n3,npat), intent(inout) :: intercepted_out
   real, dimension(    n2,n3,npat), intent(inout) :: qintercepted_out
   real, dimension(    n2,n3,npat), intent(inout) :: wshed_out
   real, dimension(    n2,n3,npat), intent(inout) :: qwshed_out
   real, dimension(    n2,n3,npat), intent(inout) :: throughfall_out
   real, dimension(    n2,n3,npat), intent(inout) :: qthroughfall_out
   real, dimension(    n2,n3,npat), intent(inout) :: runoff_out
   real, dimension(    n2,n3,npat), intent(inout) :: qrunoff_out
   real, dimension(    n2,n3,npat), intent(inout) :: drainage_out
   real, dimension(    n2,n3,npat), intent(inout) :: qdrainage_out
   real, dimension(    n2,n3,npat), intent(inout) :: gpp_out
   real, dimension(    n2,n3,npat), intent(inout) :: plresp_out
   real, dimension(    n2,n3,npat), intent(inout) :: resphet_out
   real, dimension(    n2,n3,npat), intent(inout) :: growresp
   real, dimension(    n2,n3,npat), intent(inout) :: veg_ndvip
   real, dimension(    n2,n3,npat), intent(inout) :: veg_ndvic
   real, dimension(    n2,n3,npat), intent(inout) :: veg_ndvif
   real, dimension(    n2,n3,npat), intent(inout) :: rshort_gnd
   real, dimension(    n2,n3,npat), intent(inout) :: rlong_gnd
   real, dimension(    n2,n3)     , intent(out)   :: rlongup
   real, dimension(    n2,n3)     , intent(out)   :: albedt
   real, dimension(    n2,n3)     , intent(inout) :: rvv
   real, dimension(    n2,n3)     , intent(inout) :: prsv
   real, dimension(    n2,n3)     , intent(inout) :: piv
   real, dimension(    n2,n3)     , intent(inout) :: vt2da
   real, dimension(    n2,n3)     , intent(inout) :: vt2db
   real, dimension(    n2,n3)     , intent(in)    :: glat
   real, dimension(    n2,n3)     , intent(in)    :: glon
   real, dimension(    n2,n3)     , intent(inout) :: zot
   !----- Local variables. ----------------------------------------------------------------!
   integer                                        :: k2
   integer                                        :: i
   integer                                        :: j
   integer                                        :: k
   integer                                        :: ipat
   integer                                        :: nveg
   integer                                        :: nsoil
   real                                           :: soil_temp
   real                                           :: soil_fliq
   real                                           :: available_water
   real                                           :: psi_layer
   !----- External functions. -------------------------------------------------------------!
   real                           , external      :: soil_idx2water
   !---------------------------------------------------------------------------------------!



   !----- Set up some scratch variables. --------------------------------------------------!
   g_urban    = 0.
   emis_town  = 0.
   alb_town   = 0.
   ts_town    = 0.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Time interpolation factor for updating SST.  Use calculation of continuous model !
   ! time.                                                                                 !
   !---------------------------------------------------------------------------------------!
   if (iupdsst == 0) then
      timefac_sst = 0.
   else
      timefac_sst = (time - ssttime1(ifm)) / (ssttime2(ifm) - ssttime1(ifm))
   end if
   if (iupdndvi == 0) then
      timefac_ndvi = 0.
   else
      timefac_ndvi = (time - ndvitime1(ifm)) / (ndvitime2(ifm) - ndvitime1(ifm))
   end if

   albedt(:,:)  = 0.0
   rlongup(:,:) = 0.0

   jloop: do j = 1,n3
      iloop: do i = 1,n2
         k2=nint(flpw(i,j))
         piv(i,j)  = 0.5 * (pi0(k2-1,i,j) + pi0(k2,i,j) + pp(k2-1,i,j) + pp(k2,i,j))
         prsv(i,j) = exner2press(piv(i,j))
         geoht     = (zt(k2)-zm(k2-1)) * rtgt(i,j)

         atm_theta = theta(k2,i,j)
         atm_prss  = prsv(i,j)
         atm_shv   = rv(k2,i,j) / (rv(k2,i,j) + 1.)

         patch_rough(i,j,1) = waterrough

         !---------------------------------------------------------------------------------!
         !     Canopy properties.  Copy conserved variables from lowest atmospheric grid,  !
         ! and compute pressure and temperature.                                           !
         !---------------------------------------------------------------------------------!
         can_prss(i,j,1)    = reducedpress(atm_prss,atm_theta,atm_shv,geoht                &
                                          ,atm_theta,atm_shv,can_depth)
         can_exner          = press2exner(can_prss(i,j,1))
         can_theta(i,j,1)   = theta(k2,i,j)
         can_rvap(i,j,1)    = rv(k2,i,j)
         can_shv            = can_rvap(i,j,1) / (can_rvap(i,j,1) + 1.)
         can_co2(i,j,1)     = co2p(k2,i,j)
         can_temp           = extheta2temp(can_exner,atm_theta)
         can_theiv(i,j,1)   = thetaeiv(can_theta(i,j,1),can_prss(i,j,1),can_temp           &
                                      ,can_rvap(i,j,1),can_rvap(i,j,1))
         can_vpdef(i,j,1)   = vpdefil (can_prss(i,j,1),can_temp,can_shv,.true.)
         can_enthalpy       = tq2enthalpy(can_temp,can_shv,.true.)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the temperature immediately above the canopy air space and the         !
         ! specific enthalpy associated with it.                                           !
         !---------------------------------------------------------------------------------!
         atm_temp_zcan      = extheta2temp(can_exner,atm_theta)
         atm_enthalpy       = tq2enthalpy(atm_temp_zcan,atm_shv,.true.)
         !---------------------------------------------------------------------------------!

         !----- Water patch, so we set vegetation properties to zero. ---------------------!
         veg_energy(i,j,1)  = 0.0
         veg_water (i,j,1)  = 0.0
         veg_hcap  (i,j,1)  = 0.0
         
         !----- Unless we are solving water lilies or phytoplankton, these should be 0. ---!
         veg_lai   (i,j,1)  = 0.0
         veg_tai   (i,j,1)  = 0.0
         veg_agb   (i,j,1)  = 0.0

         !----- Soil properties. Except for top layer energy, everything is set to zero. --!
         soil_water (:,i,j,1) = 1.
         psibar_10d   (i,j,1) = 1.
         do k=1,mzg
            soil_energy(k,i,j,1) = tl2uint( seatp(i,j)+(seatf(i,j)-seatp(i,j))*timefac_sst &
                                          , 1.0)
         end do
         
         !----- Fluxes.  Initially they should be all zero. -------------------------------!
         hflxac_out      (i,j,1) = 0.0
         wflxac_out      (i,j,1) = 0.0
         qwflxac_out     (i,j,1) = 0.0
         eflxac_out      (i,j,1) = 0.0
         cflxac_out      (i,j,1) = 0.0
         hflxgc_out      (i,j,1) = 0.0
         wflxgc_out      (i,j,1) = 0.0
         qwflxgc_out     (i,j,1) = 0.0
         hflxvc_out      (i,j,1) = 0.0
         wflxvc_out      (i,j,1) = 0.0
         qwflxvc_out     (i,j,1) = 0.0
         transp_out      (i,j,1) = 0.0
         qtransp_out     (i,j,1) = 0.0
         intercepted_out (i,j,1) = 0.0
         qintercepted_out(i,j,1) = 0.0
         wshed_out       (i,j,1) = 0.0
         qwshed_out      (i,j,1) = 0.0
         throughfall_out (i,j,1) = 0.0
         qthroughfall_out(i,j,1) = 0.0
         runoff_out      (i,j,1) = 0.0
         qrunoff_out     (i,j,1) = 0.0
         drainage_out    (i,j,1) = 0.0
         qdrainage_out   (i,j,1) = 0.0
         gpp_out         (i,j,1) = 0.0
         plresp_out      (i,j,1) = 0.0
         resphet_out     (i,j,1) = 0.0
         growresp        (i,j,1) = 0.0

         !----- Above-ground biomass.  This should be always 0 for water patches. ---------!
         veg_agb  (i,j,1) = 0.0

         !---------------------------------------------------------------------------------!
         !     We now loop over the land patches.  Some properties such as canopy          !
         ! properties are simply copied from the water patch.                              !
         !---------------------------------------------------------------------------------!
         patchloop1: do ipat = 2,npat

            nveg = nint(leaf_class(i,j,ipat))

            soil_rough  (i,j,ipat)  = zrough
            patch_rough (i,j,ipat)  = max(zrough,zot(i,j))
            veg_rough   (i,j,ipat)  = vh2vr * veg_ht(nveg)

            veg_height  (i,j,ipat)  = veg_ht(nveg)
            veg_displace(i,j,ipat)  = vh2dh * veg_ht(nveg)
            veg_albedo  (i,j,ipat)  = albv_green(nveg)
            stom_condct (i,j,ipat)  = 1.e-6

            !------------------------------------------------------------------------------!
            !     We cannot allow the vegetation height to be above the first level.  If   !
            ! that happens, we stop the run and ask the user to coarsen delta-z.           !
            !------------------------------------------------------------------------------!
            can_depth = max(veg_height(i,j,ipat),can_depth_min)
            if (can_depth + 6. > geoht) then
               write (unit=*,fmt='(a)') '================================================='
               write (unit=*,fmt='(a)') '   DELTA-Z is too fine, the first level is'
               write (unit=*,fmt='(a)') ' beneath or too close to the top of the canopy...'
               write (unit=*,fmt='(a)') ' Try coarsening it and run again.'
               write (unit=*,fmt='(a)') ' '
               write (unit=*,fmt='(a,1x,i5)'    ) 'I          =',i
               write (unit=*,fmt='(a,1x,i5)'    ) 'J          =',j
               write (unit=*,fmt='(a,1x,i5)'    ) 'P          =',ipat
               write (unit=*,fmt='(a,1x,i5)'    ) 'NVEG       =',nveg
               write (unit=*,fmt='(a,1x,es12.5)') 'GLON       =',glon(i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'GLAT       =',glat(i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'GEOHT      =',geoht
               write (unit=*,fmt='(a,1x,es12.5)') 'VEG_HEIGHT =',veg_height(i,j,ipat)
               write (unit=*,fmt='(a,1x,es12.5)') 'CAN_DEPTH  =',can_depth
               write (unit=*,fmt='(a)') '================================================='
               call abort_run('Delta-z is too thin!','sfcinit_nofile','leaf3_init.f90')
            end if

            !----- Canopy air properties.  Initially we just assume same properties. ------!
            can_prss  (i,j,ipat) = can_prss (i,j,1)
            can_theta (i,j,ipat) = can_theta(i,j,1)
            can_theiv (i,j,ipat) = can_theiv(i,j,1)
            can_vpdef (i,j,ipat) = can_vpdef(i,j,1)
            can_rvap  (i,j,ipat) = can_rvap (i,j,1)
            can_co2   (i,j,ipat) = can_co2  (i,j,1)

            !----- Fluxes.  Initially they should be all zero. ----------------------------!
            hflxac_out      (i,j,ipat) = 0.0
            wflxac_out      (i,j,ipat) = 0.0
            qwflxac_out     (i,j,ipat) = 0.0
            eflxac_out      (i,j,ipat) = 0.0
            cflxac_out      (i,j,ipat) = 0.0
            hflxgc_out      (i,j,ipat) = 0.0
            wflxgc_out      (i,j,ipat) = 0.0
            qwflxgc_out     (i,j,ipat) = 0.0
            hflxvc_out      (i,j,ipat) = 0.0
            wflxvc_out      (i,j,ipat) = 0.0
            qwflxvc_out     (i,j,ipat) = 0.0
            transp_out      (i,j,ipat) = 0.0
            qtransp_out     (i,j,ipat) = 0.0
            intercepted_out (i,j,ipat) = 0.0
            qintercepted_out(i,j,ipat) = 0.0
            wshed_out       (i,j,ipat) = 0.0
            qwshed_out      (i,j,ipat) = 0.0
            throughfall_out (i,j,ipat) = 0.0
            qthroughfall_out(i,j,ipat) = 0.0
            runoff_out      (i,j,ipat) = 0.0
            qrunoff_out     (i,j,ipat) = 0.0
            drainage_out    (i,j,ipat) = 0.0
            qdrainage_out   (i,j,ipat) = 0.0
            gpp_out         (i,j,ipat) = 0.0
            plresp_out      (i,j,ipat) = 0.0
            resphet_out     (i,j,ipat) = 0.0
            growresp        (i,j,ipat) = 0.0


            available_water = 0.0
            do k = 1,mzg

               nsoil = nint(soil_text(k,i,j,ipat))

               !---------------------------------------------------------------------------!
               !     For persistent wetlands (bogs, marshes, fens, swamps) and irrigated   !
               ! crops, initialize with saturated soil.  Currently, this corresponds to    !
               ! leaf classes 16, 17, and 20.  Otherwise, use the user-defined input       !
               ! profile.                                                                  !
               !---------------------------------------------------------------------------!
               select case (nint(leaf_class(i,j,ipat)))
               case (16,17,20)
                  soil_water(k,i,j,ipat) = slmsts(nsoil)
               case default
                  soil_water(k,i,j,ipat) = soil_idx2water(slmstr(k),nsoil)
               end select

               !---------------------------------------------------------------------------!
               !     By default, initialize soil internal energy at a temperature equal to !
               ! can_temp + stgoff(k).  If the temperature is initially below triple       !
               ! point, we assume all soil water to be frozen, otherwise we assume all     !
               ! water to be liquid.  If soil temperature is exactly at the triple point,  !
               ! assume that the fraction is 50%.                                          !
               !---------------------------------------------------------------------------!
               soil_temp = can_temp + stgoff(k)
               if (soil_temp == t3ple) then
                  soil_fliq = 0.5
               elseif (soil_temp > t3ple) then
                  soil_fliq = 1.0
               else
                  soil_fliq = 0.0
               end if

               soil_energy(k,i,j,ipat) = cmtl2uext( slcpd(nsoil)                           &
                                                  , soil_water(k,i,j,ipat)                 &
                                                  , soil_temp,soil_fliq )

               !------ Integrate the relative potential. ----------------------------------!
               if (k >= kroot(nveg) .and. nsoil /= 13) then
                  psi_layer       = leaf3_matric_potential(nsoil,soil_water(k,i,j,ipat))
                  available_water = available_water                                        &
                                  + max(0., (psi_layer-psiwp(nsoil))                       &
                                          / (psild(nsoil) - psiwp(nsoil)) )                &
                                  * (slz(k+1) - slz(k))
               end if
               !---------------------------------------------------------------------------!
            end do

            !----- Normalise the available water. -----------------------------------------!
            available_water      = available_water / abs(slz(kroot(nveg)))
            psibar_10d(i,j,ipat) = available_water
            !------------------------------------------------------------------------------!

            !------ Surface water, if any, will initially occupy just the first level. ----!
            do k = 1,mzs
               sfcwater_mass(k,i,j,ipat) = 0.
               sfcwater_energy(k,i,j,ipat) = 0.
               sfcwater_depth(k,i,j,ipat) = 0.
            end do

            !------------------------------------------------------------------------------!
            !    For persistent wetlands (bogs, marshes, fens, swamps), initialise with 10 !
            ! cm water depth and temperature not colder than the triple point.  Currently, !
            ! this corresponds to leaf classes 17 and  20.  Glaciers (leaf class 2) will   !
            ! be initialised with a thick layer of compact ice (currently 6 m), not warmer !
            ! than the triple point, so we make sure there will enough ice there.  The ice !
            ! will be eventually split into several layers.                                !
            !------------------------------------------------------------------------------!
            select case (nint(leaf_class(i,j,ipat)))
            case (2)
               sfcwater_depth (1,i,j,ipat) = 6.
               sfcwater_mass  (1,i,j,ipat) = idns * sfcwater_depth(1,i,j,ipat)
               sfcwater_energy(1,i,j,ipat) = tl2uint(min(t3ple,can_temp),0.0)
            case (17,20)
               sfcwater_depth (1,i,j,ipat) = .1
               sfcwater_mass  (1,i,j,ipat) = wdns * sfcwater_depth(1,i,j,ipat)
               sfcwater_energy(1,i,j,ipat) = tl2uint(max(t3ple,can_temp),1.0)
            end select

            !------------------------------------------------------------------------------!
            !    If there is initial snow information, add the information into the first  !
            ! layer.                                                                       !
            !------------------------------------------------------------------------------!
            if (snow_mass(i,j) > 0.) then
               sfcwater_energy(1,i,j,ipat) = ( sfcwater_energy(1,i,j,ipat)                 &
                                             * sfcwater_mass  (1,i,j,ipat)                 &
                                             + min(t3ple,can_temp)*cice*snow_mass(i,j) )   &
                                           / (sfcwater_mass(1,i,j,ipat) + snow_mass(i,j))
               sfcwater_mass  (1,i,j,ipat) = sfcwater_mass (1,i,j,ipat) + snow_mass(i,j)
               !---------------------------------------------------------------------------!
               !    Add a depth of 5x the equivalent of liquid depth so the extra snow is  !
               ! fresh and fluffy.                                                         !
               !---------------------------------------------------------------------------!
               sfcwater_depth(1,i,j,ipat)  = sfcwater_depth(1,i,j,ipat)                    &
                                           + snow_mass(i,j) * 5. * wdnsi
            end if

            !----- We must have either nothing or a single layer at this point. -----------!
            if (sfcwater_mass(1,i,j,ipat) > 0.) then
               sfcwater_nlev(i,j,ipat) = 1.
            else
               sfcwater_nlev(i,j,ipat) = 0.
            end if

            !------------------------------------------------------------------------------!
            !     Compute some vegetation properties, such as LAI, heat capacity, and      !
            ! update energy whilst preserving the temperature and liquid fraction.  We     !
            ! assign a dummy heat capacity because the routine will transform the energy   !
            ! back to temperature.                                                         !
            !------------------------------------------------------------------------------!
            veg_temp = can_temp
            if (veg_temp > t00) then
               veg_fliq = 1.0
            elseif (veg_temp < t00) then
               veg_fliq = 0.0
            else
               veg_fliq = 0.5
            end if
            veg_hcap  (i,j,ipat) = hcapveg_ref * max(veg_height(i,j,ipat),hcapveg_hmin)
            veg_water (i,j,ipat) = 0.
            veg_energy(i,j,ipat) = cmtl2uext(veg_hcap(i,j,ipat),veg_water(i,j,ipat)        &
                                            ,veg_temp,veg_fliq)
            call veg_misc_update( ifm                                                      &
                                , patch_area      (i,j,ipat) , leaf_class       (i,j,ipat) &
                                , veg_fracarea    (i,j,ipat) , veg_lai          (i,j,ipat) &
                                , veg_tai         (i,j,ipat) , veg_rough        (i,j,ipat) &
                                , veg_height      (i,j,ipat) , veg_displace     (i,j,ipat) &
                                , veg_albedo      (i,j,ipat) , veg_ndvip        (i,j,ipat) &
                                , veg_ndvic       (i,j,ipat) , veg_ndvif        (i,j,ipat) &
                                , veg_agb         (i,j,ipat) , veg_energy       (i,j,ipat) &
                                , veg_water       (i,j,ipat) , veg_hcap         (i,j,ipat) &
                                , psibar_10d      (i,j,ipat) )
            !------------------------------------------------------------------------------!




            call leaf3_grndvap( soil_energy(mzg,i,j,ipat) , soil_water     (mzg,i,j,ipat)  &
                              , soil_text  (mzg,i,j,ipat) , sfcwater_energy(  1,i,j,ipat)  &
                              , sfcwater_nlev  (i,j,ipat) , can_rvap           (i,j,ipat)  &
                              , can_prss       (i,j,ipat) , ground_rsat        (i,j,ipat)  &
                              , ground_rvap    (i,j,ipat) , ground_temp        (i,j,ipat)  &
                              , ground_fliq    (i,j,ipat) )
         end do patchloop1

         !---------------------------------------------------------------------------------!
         !     In this second loop we initialise the radiation, but before we must define  !
         ! some other variables.                                                           !
         !---------------------------------------------------------------------------------!
         patchloop2: do ipat = 1,npatch
            !------------------------------------------------------------------------------!
            !     Initialise various canopy air space and temporary surface water vari-    !
            ! ables.                                                                       !
            !------------------------------------------------------------------------------!
            if (ipat == 1) then
               call leaf3_ocean_diag(ifm,mzg,seatp (i,j),seatf (i,j),soil_energy(:,i,j,1))
            end if
            call leaf3_soilsfcw_diag(ipat,mzg,mzs,soil_energy       (:,i,j,ipat)           &
                                                 ,soil_water        (:,i,j,ipat)           &
                                                 ,soil_text         (:,i,j,ipat)           &
                                                 ,sfcwater_nlev     (  i,j,ipat)           &
                                                 ,sfcwater_energy   (:,i,j,ipat)           &
                                                 ,sfcwater_mass     (:,i,j,ipat)           &
                                                 ,sfcwater_depth    (:,i,j,ipat)           &
                                                 ,.true.                                   )

            call leaf3_solve_veg(ipat,mzs,leaf_class                  (i,j,ipat)           &
                                         ,soil_rough                (  i,j,ipat)           &
                                         ,patch_area                  (i,j,ipat)           &
                                         ,veg_fracarea                (i,j,ipat)           &
                                         ,veg_tai                   (  i,j,ipat)           &
                                         ,sfcwater_nlev               (i,j,ipat)           &
                                         ,sfcwater_mass             (:,i,j,ipat)           &
                                         ,sfcwater_depth            (:,i,j,ipat)           &
                                         ,.true.                                           )

            call leaf3_can_diag(ipat,can_theta        (i,j,ipat)                           &
                                    ,can_theiv        (i,j,ipat)                           &
                                    ,can_vpdef        (i,j,ipat)                           &
                                    ,can_rvap         (i,j,ipat)                           &
                                    ,leaf_class       (i,j,ipat)                           &
                                    ,can_prss         (i,j,ipat)                           &
                                    ,.true.                                                )

            call leaf3_veg_diag(veg_energy(i,j,ipat),veg_water(i,j,ipat),veg_hcap(i,j,ipat))

            call leaf3_sfcrad( mzg, mzs, ipat                                              &
                             , soil_water    (:,i,j,ipat) , soil_color    (  i,j,ipat)     &
                             , soil_text     (:,i,j,ipat) , sfcwater_depth(:,i,j,ipat)     &
                             , patch_area    (  i,j,ipat) , veg_fracarea  (  i,j,ipat)     &
                             , leaf_class    (  i,j,ipat) , veg_lai       (  i,j,ipat)     &
                             , veg_tai       (  i,j,ipat) , veg_albedo    (  i,j,ipat)     &
                             , sfcwater_nlev (  i,j,ipat) , 0.                             &
                             , 0.                         , 0.                             &
                             , 0.                         , 0.                             &
                             , 0.                         , 0.                             &
                             , cosz          (  i,j     ) , albedt        (  i,j     )     &
                             , rlongup       (  i,j     ) , rshort_gnd    (  i,j,ipat)     &
                             , rlong_gnd     (  i,j,ipat) , .true.                         )
         end do patchloop2
      end do iloop
   end do jloop

   return
end subroutine sfcinit_nofile
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine maps the input datp classes to a smaller set datq which represents  !
! the full set of LEAF-2 or LEAF-3 classes for which LSP values are defined.               !
!------------------------------------------------------------------------------------------!
subroutine leaf_datp_datq(datp,datq)
   use teb_spm_start, only : teb_spm ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real   , intent(in)       :: datp      ! Input data, the Olson Global Ecosystems dataset
   integer, intent(out)      :: datq      ! Output data, the LEAF-3 classification. 
   !----- Local variable. -----------------------------------------------------------------!
   integer, dimension(0:100) :: oge2leaf3 ! Conversion table
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Olson Global Ecosystems dataset OGE_2 (96 classes) mapped to LEAF-3 classes (see  !
   ! leaf3_document). 97 and 98 are not used, 99 is Goode Homolosine empty data, and 100   !
   ! is missing data. They are all mapped to ocean.                                        !
   !---------------------------------------------------------------------------------------!
   oge2leaf3 = (/                                      0, &  !   0 -   0
                  19,  8,  4,  5,  6,  7,  9,  3, 11, 16, &  !   1 -  10
                  10,  2, 20,  0,  0, 12, 13, 14, 18,  4, &  !  11 -  20
                   4,  4, 14, 14,  6,  6,  4,  7,  7, 15, &  !  21 -  30
                  15,  6,  7,  7, 15, 16, 16, 16, 16,  8, &  !  31 -  40
                   8,  8, 18, 17, 17, 12, 12,  7, 10,  3, &  !  41 -  50
                  10, 10, 11, 14, 18, 18, 18, 18, 13,  6, &  !  51 -  60
                   5,  4, 11, 12,  0,  0,  0,  0,  3,  2, &  !  61 -  70
                   3, 20,  0, 17, 17, 17,  4, 14,  7,  3, &  !  71 -  80
                   3,  3,  3,  3,  3,  3,  8, 12,  7,  6, &  !  81 -  90
                  18, 15, 15, 15, 19,  0,  0,  0,  0,  0  /) !  91 - 100
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     If TEB_SPM, we switch OGE classes 95 and 50 to 21 (very urban).  Otherwise, leave !
   ! the standard.                                                                         !
   !---------------------------------------------------------------------------------------!
   if (teb_spm == 1) then
      oge2leaf3(95) = 21
      oge2leaf3(50) = 21
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert the data.  In the unlikely case that the dataset lies outside the range,  !
   ! re-assign it to ocean.                                                                !
   !---------------------------------------------------------------------------------------!
   select case (nint(datp))
   case (0:95)
      datq = oge2leaf3(nint(datp))
   case default
      datq = 0
   end select
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf_datp_datq
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine maps the input datp soil classes to a smaller set datsoil which     !
! represents the full set of LEAF-2 classes for which soil parameter values are defined.   !
! This sub-routine is now deprecated as we converted the soil types in the input files to  !
! have the same classification.                                                            !
!      This is currently deprecated, because the input files should have indices already   !
! converted to the LEAF-3 classes.                                                         !
!------------------------------------------------------------------------------------------!
subroutine leaf_datp_datsoil(datp,datsoil)


   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real   , intent(in)       :: datp      ! Input data, the FAO dataset
   integer, intent(out)      :: datsoil   ! Output data, the LEAF-3 classification. 
   !----- Local variable. -----------------------------------------------------------------!
   integer, dimension(0:133) :: fao2leaf3 ! Conversion table
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      (Bob 9/14/2000) This table maps FAO soil units numbered 0 to 132, plus our own   !
   ! missing value designated 133, to the USDA soil textural classes.  FAO classes [0]     !
   ! (ocean), [1, 7, 27, 42, 66, 69, 77, 78, 84, 98, 113] (soils not used in original FAO  !
   ! dataset), [132] (water), and [133] (our missing value) are all mapped to a default    !
   ! class of sandy clay loam in case they happen to correspond to a land surface area in  !
   ! the landuse dataset that RAMS uses to define land area.  We wrote missing value class !
   ! 133 to the RAMS FAO files  whenever a negative value (which is out of range of        !
   ! defined values) was found in the original FAO dataset, which occurred in about 2.6%   !
   ! of the pixels.  For the remaining FAO classes, a cross reference table to Zobler soil !
   ! texture classes that was provided, plus our own cross referencing table from Zobler   !
   ! to USDA classes listed below, provides the mapping from FAO to USDA.  In this         !
   ! mapping, we use only organic USDA classes and omit nonorganic classes since most      !
   ! soils do contain organic matter, and organic content information is not provided in   !
   ! the Zobler classes.                                                                   !
   !                                                                                       !
   !          |--------------------------------+--------------------------------|          !
   !          |          Zobler Class          |           USDA Class           |          !
   !          |--------------------------------+--------------------------------|          !
   !          | 1  Coarse                      |  2  Loamy sand                 |          !
   !          | 2  Medium                      |  4  Silt loam                  |          !
   !          | 3  Fine                        |  8  Clay loam                  |          !
   !          | 4  Coarse-medium               |  3  Sandy loam                 |          !
   !          | 5  Coarse-fine                 |  6  Sandy clay loam            |          !
   !          | 6  Medium-fine                 |  7  Silty clay loam            |          !
   !          | 7  Coarse-medium-fine          |  6  Sandy clay loam            |          !
   !          | 8  Organic matter              |  5  Loam                       |          !
   !          |--------------------------------+--------------------------------|          !
   !          |                                | ... not used:                  |          !
   !          |                                |  1  Sand                       |          !
   !          |                                |  9  Sandy clay                 |          !
   !          |                                | 10  Silty clay                 |          !
   !          |                                | 11  Clay                       |          !
   !          |                                | 12  Peat                       |          !
   !          |--------------------------------+--------------------------------|          !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!

   fao2leaf3 = (/                                      6  &  !   0 -   0
                ,  6,  4,  4,  7,  7,  8,  6,  4,  4,  4  &  !   1 -  10
                ,  7,  4,  4,  4,  8,  4,  8,  4,  4,  8  &  !  11 -  20
                ,  4,  2,  4,  4,  4,  4,  6,  8,  8,  8  &  !  21 -  30
                ,  4,  8,  8,  2,  6,  4,  7,  4,  4,  3  &  !  31 -  40
                ,  4,  6,  7,  4,  4,  4,  4,  4,  4,  4  &  !  41 -  50
                ,  4,  4,  4,  4,  4,  4,  2,  4,  4,  2  &  !  51 -  60
                ,  4,  3,  4,  2,  7,  6,  4,  4,  6,  8  &  !  61 -  70
                ,  8,  7,  2,  5,  4,  5,  6,  6,  4,  2  &  !  71 -  80
                ,  2,  2,  4,  6,  2,  2,  2,  2,  2,  4  &  !  81 -  90
                ,  2,  2,  2,  4,  2,  4,  3,  6,  2,  7  &  !  91 - 100
                ,  4,  4,  4,  8,  8,  8,  3,  7,  4,  4  &  ! 101 - 110
                ,  4,  3,  6,  4,  2,  4,  4,  4,  2,  2  &  ! 111 - 120
                ,  2,  4,  6,  4,  4,  7,  7,  6,  3,  2  &  ! 121 - 130
                ,  2,  6,  6                              /) ! 131 - 133

   datsoil = fao2leaf3(nint(datp))

   return
end subroutine leaf_datp_datsoil
!==========================================================================================!
!==========================================================================================!

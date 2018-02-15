!==========================================================================================!
!==========================================================================================!
! MODULE: PLANT_HYDRO
!> \brief Solve Farquhar Photosynthesis model together with Katul optimization-based
!> stomata conductance model.
!> \details The optimization framework is to maximize A - stoma_lambda * E every
!> model time step. Currently the optimization is acheived using Newton's method
!> with initial values as stom. cond. from last timestep.\n
!> The references are:\n
!>       Xu et al. (2016) Diversity in plant hydraulic traits explains seasonal and
!> inter-annual variations in vegetation dynamics in seasonally dry tropical
!> forests.  New Phytologist\n
!>\n
!>       Manzoni, S., G. Vico, et al. (2011). "Optimizing stomatal conductance for maximum
!> carbon gain under water stress: a meta-analysis across plant functional types
!> and climates." Functional Ecology 25(3): 456-467.\n
!>\n
!>       Vico, G., S. Manzoni, et al. (2013). "A perspective on optimal leaf
!> stomatal conductance under CO2 and light co-limitations." Agricultural
!> and Forest Meteorology 182–183(0): 191-199.\n
!>\n
!>       Katul, G., S. Manzoni, et al. (2010). "A stomatal optimization theory to describe
!> the effects of atmospheric CO2 on leaf photosynthesis and transpiration."
!> Annals of Botany 105(3): 431-442.
!> \author Xiangtao Xu, 15 Feb. 2018
!==========================================================================================!
!==========================================================================================!
Module farq_katul

Contains
   !=======================================================================================!
   !=======================================================================================!
   ! SUBROUTINE KATUL_LPHYS       
   !> \brief   Main driver to calculate Farquhar-Katul photosynthesis system.
   !> Alternative to lphysio_full in farq_leuning.
   !> \details The realized Vcmax is reduced under low leaf water potential while leaf dark
   !> respiration keeps the same.
   !> \author Xiangtao Xu, 15 Feb. 2018
   !---------------------------------------------------------------------------------------!

  subroutine katul_lphys(can_prss,can_shv,can_co2,ipft,leaf_par,leaf_temp                 &
                        ,lint_shv,green_leaf_factor,leaf_aging_factor,llspan,vm0in        &
                        ,leaf_gbw,leaf_psi,last_gV,last_gJ,A_open,A_closed,gsw_open       &
                        ,gsw_closed,lsfc_shv_open,lsfc_shv_closed,lsfc_co2_open           &
                        ,lsfc_co2_closed,lint_co2_open,lint_co2_closed,leaf_resp          &
                        ,vmout,comppout,limit_flag)
      
    use rk4_coms       , only : tiny_offset              & ! intent(in)
                              , effarea_transp           ! ! intent(in)
    use pft_coms       , only : photosyn_pathway         & ! intent(in)
                              , phenology                & ! intent(in)
                              , Vm0                      & ! intent(in)
                              , Rd0                      & ! intent(in)
                              , vm_hor                   & ! intent(in)
                              , vm_low_temp              & ! intent(in)
                              , vm_high_temp             & ! intent(in)
                              , vm_decay_e               & ! intent(in)
                              , vm_q10                   & ! intent(in)
                              , rd_hor                   & ! intent(in)
                              , rd_low_temp              & ! intent(in)
                              , rd_high_temp             & ! intent(in)
                              , rd_decay_e               & ! intent(in)
                              , rd_q10                   & ! intent(in)
                              , cuticular_cond           & ! intent(in)
                              , dark_respiration_factor  & ! intent(in)
                              , quantum_efficiency       & ! intent(in)
                              , leaf_psi_tlp             & ! intent(in)
                              , stoma_lambda             & ! intent(in)
                              , stoma_beta               & ! intent(in)
                              , leaf_psi_min             ! ! intent(in)
    use consts_coms    , only : t00                      & ! intent(in)
                              , umol_2_mol               & ! intent(in)
                              , mmdry1000                & ! intent(in)
                              , mmh2o                    & ! intent(in)
                              , mmdry                    & ! intent(in)
                              , mmdryi                   & ! intent(in)
                              , epi                      & ! intent(in)
                              , ep                       ! ! intent(in)
    use therm_lib      , only : eslf                     ! ! function
    use ed_misc_coms   , only : current_time             ! ! intent(in)
    use physiology_coms, only : gbw_2_gbc                & ! intent(in)
                              , gsw_2_gsc                & ! intent(in)
                              , iphysiol                 & ! intent(in)
                              , h2o_plant_lim            & ! intent(in)
                              , o2_ref                   & ! intent(in)
                              , kco2_refval              & ! intent(in)
                              , kco2_hor                 & ! intent(in)
                              , kco2_q10                 & ! intent(in)
                              , ko2_refval               & ! intent(in)
                              , ko2_hor                  & ! intent(in)
                              , ko2_q10                  & ! intent(in)
                              , compp_refval             & ! intent(in)
                              , compp_q10                & ! intent(in)
                              , compp_hor               
    use phenology_coms , only : vm0_tran                 & ! intent(in)
                              , vm0_slope                & ! intent(in)
                              , vm0_amp                  & ! intent(in)
                              , vm0_min                  ! ! intent(in)

    implicit none

      !------ Arguments. ------------------------------------------------------------------!
      real(kind=4), intent(in)    :: can_prss          ! Canopy air pressure    [       Pa]
      real(kind=4), intent(in)    :: can_shv           ! Canopy air sp. hum.    [    kg/kg]
      real(kind=4), intent(in)    :: can_co2           ! Canopy air CO2         [ µmol/mol]
      integer     , intent(in)    :: ipft              ! Plant functional type  [      ---]
      real(kind=4), intent(in)    :: leaf_par          ! Absorbed PAR           [     W/m²]
      real(kind=4), intent(in)    :: leaf_temp         ! Leaf temperature       [        K]
      real(kind=4), intent(in)    :: lint_shv          ! Leaf interc. sp. hum.  [    kg/kg]
      real(kind=4), intent(in)    :: green_leaf_factor ! Frac. of on-allom. gr. [      ---]
      real(kind=4), intent(in)    :: leaf_aging_factor ! Ageing parameter       [      ---]
      real(kind=4), intent(in)    :: llspan            ! Leaf life span         [     mnth]
      real(kind=4), intent(in)    :: vm0in             ! Input Vm0              [µmol/m²/s]
      real(kind=4), intent(in)    :: leaf_gbw          ! B.lyr. cnd. of H2O     [  kg/m²/s]
      real(kind=4), intent(in)    :: leaf_psi          ! leaf water potential   [        m]
      real(kind=4), intent(inout) :: last_gV           ! gs for last timestep   [  kg/m²/s]
      real(kind=4), intent(inout) :: last_gJ           ! gs for last timestep   [  kg/m²/s]
      real(kind=4), intent(out)   :: A_open            ! Photosyn. rate (op.)   [µmol/m²/s]
      real(kind=4), intent(out)   :: A_closed          ! Photosyn. rate (cl.)   [µmol/m²/s]
      real(kind=4), intent(out)   :: gsw_open          ! St. cnd. of H2O  (op.) [  kg/m²/s]
      real(kind=4), intent(out)   :: gsw_closed        ! St. cnd. of H2O  (cl.) [  kg/m²/s]
      real(kind=4), intent(out)   :: lsfc_shv_open     ! Leaf sfc. sp.hum.(op.) [    kg/kg] 
      real(kind=4), intent(out)   :: lsfc_shv_closed   ! Leaf sfc. sp.hum.(cl.) [    kg/kg]
      real(kind=4), intent(out)   :: lsfc_co2_open     ! Leaf sfc. CO2    (op.) [ µmol/mol]
      real(kind=4), intent(out)   :: lsfc_co2_closed   ! Leaf sfc. CO2    (cl.) [ µmol/mol]
      real(kind=4), intent(out)   :: lint_co2_open     ! Intercell. CO2   (op.) [ µmol/mol]
      real(kind=4), intent(out)   :: lint_co2_closed   ! Intercell. CO2   (cl.) [ µmol/mol]
      real(kind=4), intent(out)   :: leaf_resp         ! Leaf respiration rate  [µmol/m²/s]
      real(kind=4), intent(out)   :: vmout             ! Max. Rubisco capacity  [µmol/m²/s]
      real(kind=4), intent(out)   :: comppout          ! GPP compensation point [ µmol/mol]
      integer     , intent(out)   :: limit_flag        ! Photosyn. limit. flag  [      ---]
      !----- External function. -----------------------------------------------------------!
      real(kind=4)    , external      :: sngloff     ! Safe double -> single precision
      !----- Local Variables    -----------------------------------------------------------!
      real(kind=4)                :: blyr_cond_h2o      ! leaf boundary layer h2o conductance
      real(kind=4)                :: blyr_cond_co2      ! leaf boundary layer co2 conductance
      real(kind=4)                :: stom_cond_h2o      ! leaf stomata h2o conductance
      real(kind=4)                :: stom_cond_co2      ! leaf stomata co2 conductance
      real(kind=4)                :: leaf_o2            ! O2 concentration  mmol/mol
      real(kind=4)                :: par                ! PAR          micromol/m2/s
      real(kind=4)                :: leaf_temp_degC     ! leaf_temperature in degree C
      real(kind=4)                :: can_vpr_prss       ! canopy vapor pressure in kPa
      real(kind=4)                :: Vcmax              ! current Vcmax  umol/m2/s
      real(kind=4)                :: Vcmax25            ! current Vcmax at 25 degC umol/m2/s
      real(kind=4)                :: Vcmax15            ! current Vcmax at 15 degC umol/m2/s
      real(kind=4)                :: Jmax               ! current Jmax  umol/m2/s
      real(kind=4)                :: Jmax25             ! current Jmax at 25 degC umol/m2/s
      real(kind=4)                :: Jmax15             ! current Jmax at 15 degC umol/m2/s
      real(kind=4)                :: Jrate              ! current Jrate umol/m2/s
      real(kind=4)                :: Rdark              ! current dark respiration rate umol/m2/s
      real(kind=4)                :: cuticular_gsc      ! current cuticular_conductance for CO2 umol/m2/s
      real(kind=4)                :: lambda             ! current lambda factor    numo/mol/kPa
      real(kind=4)                :: down_factor        ! photosynthetic down-regulation factor
      real(kind=4)                :: aero_resistance    ! aerodynamic resistance
      real(kind=4)                :: accepted_gsc       ! gsc solved from optimization scheme umol/m2/s
      real(kind=4)                :: accepted_fc        ! CO2 flux solved from optimization scheme umol/m2/s
      real(kind=4)                :: accepted_ci        ! Internal CO2 concentration ppm
      real(kind=4)                :: leaf_vpr_prss      ! leaf internal vapor pressure  in kPa
      real(kind=4)                :: cp                 ! constants in photosynthesis  umol/mol
      real(kind=4)                :: kc                 ! constants in photosynthesis  umol/mol
      real(kind=4)                :: ko                 ! constants in photosynthesis  umol/mol
      real(kind=4)                :: dfcdg              ! variables used in optimization scheme
      real(kind=4)                :: dfedg              ! variables used in optimization scheme
      real(kind=4)                :: d2fcdg2            ! variables used in optimization scheme
      real(kind=4)                :: d2fedg2            ! variables used in optimization scheme
      real(kind=4)                :: delta_g            ! variables used in optimization scheme
      real(kind=4)                :: a1gk,a2gk          ! variables used in optimization scheme
      real(kind=4)                :: k1ci,k2ci,k3ci,k4ci! variables used in optimization scheme
      real(kind=4)                :: test_gsc           ! variables used in optimization scheme
      real(kind=4)                :: test_gV            ! variables used in optimization scheme
      real(kind=4)                :: test_fcV           ! variables used in optimization scheme
      real(kind=4)                :: test_ciV           ! variables used in optimization scheme
      real(kind=4)                :: test_gJ            ! variables used in optimization scheme
      real(kind=4)                :: test_fcJ           ! variables used in optimization scheme
      real(kind=4)                :: test_ciJ           ! variables used in optimization scheme
      integer                     :: iter               ! variables used in optimization scheme
      real(kind=4)                :: testfc             ! variables used in optimization scheme
      real(kind=4)                :: testci             ! variables used in optimization scheme
      real(kind=4)                :: greeness           ! Leaf "Greeness"           [   0 to 1]
      real           ,parameter   :: Jmax_vmhor_coef = 5./7.  ! fraction of Jmax  vmhor to Vcmax vmhor estimated from Kattge et al. 2007
      integer                     :: k
      logical                     :: is_resolvable
      logical, parameter          :: quality_check = .false.
    

      !-------------------    Define some constants....
      leaf_o2           = o2_ref * 1000.                             ! convert to mmol/mol
      par               = leaf_par * 4.6                             ! convert to micromol/m2/s
      leaf_temp_degC    = leaf_temp - t00                            ! convert to degC
      leaf_vpr_prss     = eslf(leaf_temp) * 0.001                    ! in kPa
      can_vpr_prss      = can_shv * can_prss * 0.001                 ! in kPa
  
    
      !------------------------------------------------------------------------------------!
      ! correcting for light-phenology
      !------------------------------------------------------------------------------------!
      select case(phenology(ipft))
      case (3)
         !------ Light-controlled phenology. ----------------------------------------------!
         Vcmax15    = vm0_amp / (1.0 + (llspan/vm0_tran)**vm0_slope) + vm0_min
      case default
         !------ Other phenologies, no distinction on Vm0. --------------------------------!
         Vcmax15    = vm0in
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! Calculate temperature dependence and other derived variables
      !------------------------------------------------------------------------------------!
      select case (iphysiol) 
      case (0,1)
          ! We go for the Arrhenius form as in farq_leuning module
          Vcmax25 = Vcmax15 &
                  * mod_arrhenius(298.15,                 &
                                  vm_hor(ipft),           &
                                  vm_low_temp(ipft),      &
                                  vm_high_temp(ipft),     &
                                  vm_decay_e(ipft),       &
                                  .true.)

          Jmax25 = Vcmax25 * 1.97 ! Leuning 1997, Wullschlegger et al. 1991
                                  ! Other values like 1.67 is also reported
                                  ! based on the temperature dependence used

          Jmax15 = Jmax25 &
                 /  mod_arrhenius(298.15,                 &
                        vm_hor(ipft) * Jmax_vmhor_coef,   & 
                        vm_low_temp(ipft),                &
                        vm_high_temp(ipft),               &
                        vm_decay_e(ipft),                 &
                        .true.)
    
          ! calculate the Vcmax Jmax and Rd at current T
          Vcmax = Vcmax15                              &
                * mod_arrhenius(leaf_temp,             &
                       vm_hor(ipft),                   &
                       vm_low_temp(ipft),              &
                       vm_high_temp(ipft),             &
                       vm_decay_e(ipft),               &
                       .true.)

          Jmax = Jmax15                                &
               * mod_arrhenius(leaf_temp,              &
                     vm_hor(ipft) * Jmax_vmhor_coef,   &
                     vm_low_temp(ipft),                &
                     vm_high_temp(ipft),               &
                     vm_decay_e(ipft),                 &
                     .true.)

          Rdark = Vcmax15                              &
                * dark_respiration_factor(ipft)        &
                * mod_arrhenius(leaf_temp,             &
                      rd_hor(ipft),                    &        
                      rd_low_temp(ipft),               &
                      rd_high_temp(ipft),              &
                      rd_decay_e(ipft),                &
                      .true.)

          cp = compp_refval * mod_arrhenius(leaf_temp,compp_hor,0.,0.,0.,.false.) / umol_2_mol
          kc = kco2_refval * mod_arrhenius(leaf_temp,kco2_hor,0.,0.,0.,.false.)   / umol_2_mol
          ko = kco2_refval * mod_arrhenius(leaf_temp,ko2_hor,0.,0.,0.,.false.)    / umol_2_mol

      case (2,3)
          ! Use Q10 from Collatz et al. 1991 
          Vcmax25 = Vcmax15                             &
                  * mod_collatz(298.15,                 &
                                vm_q10(ipft),           &
                                vm_low_temp(ipft),      &
                                vm_high_temp(ipft),     &
                                vm_decay_e(ipft),       &
                                .true.)

          Jmax25 = Vcmax25 * 1.97 ! Leuning 1997, Wullschlegger et al. 1991
                                  ! Other values like 1.67 is also reported
                                  ! based on the temperature dependence used

          Jmax15 = Jmax25 &
                 /  mod_collatz(298.15,                 &
                        vm_q10(ipft),                   & ! assume Jmax has the same q10 as Vcmax
                        vm_low_temp(ipft),              &
                        vm_high_temp(ipft),             &
                        vm_decay_e(ipft),               &
                        .true.)
    
          ! calculate the Vcmax Jmax and Rd at current T
          Vcmax = Vcmax15                              &
                * mod_collatz(leaf_temp,               &
                       vm_q10(ipft),                   &
                       vm_low_temp(ipft),              &
                       vm_high_temp(ipft),             &
                       vm_decay_e(ipft),               &
                       .true.)

          Jmax = Jmax15                                &
               * mod_collatz(leaf_temp,                &
                     vm_q10(ipft),                     &
                     vm_low_temp(ipft),                &
                     vm_high_temp(ipft),               &
                     vm_decay_e(ipft),                 &
                     .true.)

          Rdark = Vcmax15                              &
                * dark_respiration_factor(ipft)        &
                * mod_collatz(leaf_temp,               &
                      rd_q10(ipft),                    &        
                      rd_low_temp(ipft),               &
                      rd_high_temp(ipft),              &
                      rd_decay_e(ipft),                &
                      .true.)

          cp = compp_refval * mod_collatz(leaf_temp,compp_q10,0.,0.,0.,.false.) / umol_2_mol
          kc = kco2_refval * mod_collatz(leaf_temp,kco2_q10,0.,0.,0.,.false.)   / umol_2_mol
          ko = kco2_refval * mod_collatz(leaf_temp,ko2_q10,0.,0.,0.,.false.)    / umol_2_mol

      case (4)
          ! Vcmax, Jmax Temperature dependence according to Harley et al. 1991
          ! CLM parameters for Rd Bonan et al. 2011
          ! This was used in Xu et al. 2016 New Phyt.

          Vcmax25 = Vcmax15                               &
                  / harley_arrhenius(15.+273.15,298.15,   &
                                     58.52,          & ! Hv
                                     0.710,          & ! Sv
                                     220.0)            ! Hd
          Jmax25 = Vcmax25 * 1.97

          ! calculate the Vcmax Jmax and Rd at current T
          Vcmax = Vcmax25                                 &
                * harley_arrhenius(leaf_temp,298.15,        &
                                58.52,                  & ! Hv
                                0.710,                  & ! Sv
                                220.0)                    ! Hd


          ! Jmax
          Jmax = Jmax25                                   &
               * harley_arrhenius(leaf_temp,298.15,&
                                37.10,                  & ! Hv
                                0.710,                  & ! Sv
                                220.)                     ! Hd

          ! Rd
          Rdark = Vcmax25 * dark_respiration_factor(ipft)   &
                * harley_arrhenius(leaf_temp,298.15,      &
                                64.5,                       & ! Hv
                                0.71,                       & ! Sv
                                220.)                         ! Hd
  
          cp = exp(19.02-37830./(8.314* leaf_temp ))
          kc = exp(38.05-79430./(8.314* leaf_temp ))
          ko = exp(20.30-36380./(8.314* leaf_temp ))

      end select
    
      !------------------------------------------------------------------------------------!

      ! calculate greeness
      if (leaf_aging_factor > 0.01 .and. green_leaf_factor > 0.0001) then
         greeness = leaf_aging_factor / green_leaf_factor
      else
         greeness = 1.0
      end if

      !------------------------------------------------------------------------------------!
      ! correcting for water stress impact on realized Vcmax
      !------------------------------------------------------------------------------------!
      select case (h2o_plant_lim)
      case (0,1,2,3)
          ! use fsw to account for water stress outside of this module
          down_factor = 1.
          lambda =  stoma_lambda(ipft) * can_co2 / 400.
      case (4)
          ! down scale Vcmax, Jmax, lambda using leaf_psi
          ! parameters are kind of arbitrary from Xu et al. 2016 New Phyt.
      	  down_factor   = max(1e-6,min(1.0,1. / (1. + (leaf_psi / leaf_psi_tlp(ipft)) ** 6.0)))
          lambda =  stoma_lambda(ipft) * can_co2 / 400. * exp(stoma_beta(ipft) * leaf_psi)
      end select

      Jmax      = Jmax * down_factor * greeness
      Vcmax     = Vcmax * down_factor * greeness


      !              max(0.0,min(1.0, (leaf_psi - leaf_psi_min(ipft)) &
      !              / (TLP(ipft) - leaf_psi_min(ipft))))
      ! when leaf water potential is too low, there can't be any water loss...
      ! mimicking the cease of soil evaporation near sh

       
      !------------------------------------------------------------------------------------!

      
      !------------------------------------------------------------------------------------!
      ! Solve the optimization
      !------------------------------------------------------------------------------------!
      cuticular_gsc = cuticular_cond(ipft) * 1.0e-6 ! convert to umol/m2/s

      ! initialize limit_flag as 0
      limit_flag      = 0
      ! Solve the quadratic function for light-limited photosynthesis
      if (photosyn_pathway(ipft) == 3) then
          Jrate = ((Jmax+0.385*par) -   &
                  sqrt((Jmax+0.385*par)**2-4.*0.7*0.385*Jmax*par))/1.4
      elseif (photosyn_pathway(ipft) == 4) then
          ! this is still under test
          Jrate = quantum_efficiency(ipft) * par  
      endif

      ! calcualte aerodynamic resistance
      if (leaf_gbw > 0.) then
          aero_resistance = mmdry / (leaf_gbw / gbw_2_gbc)
      else
          aero_resistance = 1e10
      endif

      is_resolvable = (Jmax /= 0.) .and. (Vcmax /= 0.) .and. &
                      (aero_resistance < 1e8) .and.(cuticular_gsc > 1e-8)

      if (is_resolvable) then
        ! 1. Rubisco-limited photosynthesis
        a1gk = Vcmax
        a2gk = kc * (1. + leaf_o2 / ko)
        k1ci = a1gk / can_co2 - Rdark / can_co2
        k2ci = a1gk * aero_resistance / can_co2 - Rdark * aero_resistance / can_co2 - 1. + a2gk / can_co2
        k3ci = -a1gk*cp/can_co2/can_co2-Rdark*a2gk/can_co2/can_co2
        k4ci = -a1gk*cp/can_co2*aero_resistance/can_co2-Rdark*a2gk/can_co2*aero_resistance/can_co2-a2gk/can_co2
       

        if(leaf_vpr_prss > can_vpr_prss) then
            ! Use newton's method to find the zero point of 
            ! start with gsw from last time
            if ((.not. isnan(last_gV)) .and. last_gV > 1e-10) then
               test_gsc = last_gV
            else
               test_gsc = cuticular_gsc
            endif
           
            do iter = 1, 500
                ! calculate dfcdg - dfedg
                call fluxsolver(test_gsc, aero_resistance, can_co2, k1ci, k2ci, k3ci, k4ci, testci, testfc)
                call deriv_fc(test_gsc, aero_resistance, k1ci, k2ci, k3ci, k4ci, can_co2, testci, dfcdg)
            
                dfedg = lambda / gsw_2_gsc * (leaf_vpr_prss - can_vpr_prss) / &
                       (1. + gbw_2_gbc/gsw_2_gsc * test_gsc * aero_resistance) ** 2

                call deriv_dfcdg(test_gsc, aero_resistance, k1ci, k2ci, k3ci, k4ci,can_co2,testci,dfcdg,d2fcdg2)
                d2fedg2 = - 2. * gbw_2_gbc/gsw_2_gsc * aero_resistance / &
                          (gbw_2_gbc/gsw_2_gsc * test_gsc * aero_resistance + 1) * dfedg

                ! calculate the derivative of dfcdg - dfedg
                if (d2fcdg2 - d2fedg2 == 0) then
                    delta_g = 0.
                else
                    delta_g = - (dfcdg - dfedg) / (d2fcdg2 - d2fedg2)
                endif

                ! control exit

                if (abs(dfcdg - dfedg) < 1e-4 .or.                              &  ! converge
                    (test_gsc < cuticular_gsc .and. dfcdg-dfedg < 0.) .or.      &  ! close stomatal
                    (test_gsc < cuticular_gsc .and. isnan(dfcdg-dfedg)) .or.    &  ! close stomatal
                    (test_gsc + delta_g <= 0.)                 .or.             &  ! unrealistic values
                    (isnan(delta_g) < 0.)                      .or.             &  ! unrealistic values
                    (delta_g == 0.)                            .or.             &  ! trapped
                    (test_gsc > 1.0 .and. dfcdg-dfedg > 0.)                     &  ! fully open stomatal
                   ) then
                    exit
                endif

                test_gsc = test_gsc + delta_g

            enddo

            ! check the case that a negative or no optimal value is found
            if (test_gsc < cuticular_gsc .or. isnan(test_gsc)) then 
                test_gsc = cuticular_gsc
            endif
                
            ! Check the case that a large gsc value is found
            if (test_gsc > 1.0) then
               if (par > 0. .and. dfcdg-dfedg > 0.) then ! light
                    test_gsc = 1.0
               else ! dark
                    test_gsc = cuticular_gsc 
               endif
            endif

            call fluxsolver(test_gsc, aero_resistance, can_co2, k1ci, k2ci, k3ci, k4ci, testci, testfc)

            last_gV = test_gsc
            
            test_gV = test_gsc
            test_fcV = testfc
            test_ciV = testci
          
        else
          ! canopy air is saturated...
            test_ciV = can_co2
            test_gV = 1000. / aero_resistance
            call fluxsolver(test_gV, aero_resistance, can_co2, k1ci, k2ci, k3ci, k4ci, testci,test_fcV)
        endif


        ! 2. Repeat calculation for light-limited case
        Jrate = Jrate * 0.25
        a1gk = Jrate
        a2gk = 2. * cp
        k1ci = a1gk / can_co2 - Rdark / can_co2
        k2ci = a1gk * aero_resistance / can_co2 - Rdark * aero_resistance / can_co2 - 1. + a2gk / can_co2
        k3ci = -a1gk*cp/can_co2/can_co2-Rdark*a2gk/can_co2/can_co2
        k4ci = -a1gk*cp/can_co2*aero_resistance/can_co2-Rdark*a2gk/can_co2*aero_resistance/can_co2-a2gk/can_co2
    
        if(leaf_vpr_prss > can_vpr_prss)then
        ! Use newton's method to find the zero point of 
        ! start with gsw from last time
            if ((.not. isnan(last_gJ)) .and. last_gJ > 1e-10) then
               test_gsc = last_gJ
            else
               test_gsc = cuticular_gsc 
            endif

            do iter = 1, 500
                ! calculate dfcdg - dfedg
                call fluxsolver(test_gsc, aero_resistance, can_co2, k1ci, k2ci, k3ci, k4ci, testci, testfc)
                call deriv_fc(test_gsc, aero_resistance, k1ci, k2ci, k3ci, k4ci, can_co2, testci, dfcdg)
            
                dfedg = lambda / gsw_2_gsc * (leaf_vpr_prss - can_vpr_prss) / &
                        (1. + gbw_2_gbc/gsw_2_gsc * test_gsc * aero_resistance) ** 2
                

                call deriv_dfcdg(test_gsc, aero_resistance, k1ci, k2ci, k3ci, k4ci,can_co2,testci,dfcdg,d2fcdg2)
                d2fedg2 = -2. * gbw_2_gbc/gsw_2_gsc * aero_resistance / &
                          (gbw_2_gbc/gsw_2_gsc * test_gsc * aero_resistance + 1) * dfedg
    
                ! calculate the derivative of dfcdg - dfedg
                if (d2fcdg2 - d2fedg2 == 0.) then
                    delta_g = 0.
                else
                    delta_g = - (dfcdg - dfedg) / (d2fcdg2 - d2fedg2)
                endif

                ! control exit
                if (abs(dfcdg - dfedg) < 1e-4 .or.                              &   ! converge
                    (test_gsc < cuticular_gsc .and. dfcdg-dfedg < 0.) .or.      &   ! close stomatal
                    (test_gsc < cuticular_gsc .and. isnan(dfcdg-dfedg)) .or.    &   ! close stomatal
                    (test_gsc + delta_g <= 0.)                 .or.             &   ! unrealistic values
                    (isnan(delta_g) < 0.)                   .or.                &   ! unrealistic values
                    (delta_g == 0.)                 .or.                        &   ! trapped
                    (test_gsc > 1.0 .and. dfcdg-dfedg > 0.)                     &   ! fullyopen stomatal
                   ) then
                            exit
                endif

                test_gsc = test_gsc + delta_g

            enddo

            if (test_gsc < cuticular_gsc .or. isnan(test_gsc)) then
                test_gsc = cuticular_gsc
            endif
            
            if (test_gsc > 1.0) then 
                if (par > 0. .and. dfcdg-dfedg>0.) then ! light
                    test_gsc = 1.0
                else ! dark
                    test_gsc = cuticular_gsc 
                endif
            endif

            call fluxsolver(test_gsc, aero_resistance, can_co2, k1ci, k2ci, k3ci, k4ci, testci, testfc)
            last_gJ = test_gsc
            test_gJ = test_gsc
            test_fcJ = testfc
            test_ciJ = testci

        else
            test_ciJ = can_co2
            test_gJ = 1000. / aero_resistance
            call fluxsolver(test_gJ, aero_resistance, can_co2, k1ci, k2ci, k3ci, k4ci, testci,test_fcJ)
        endif


        if(test_fcV < test_fcJ)then
            accepted_fc  = test_fcV
            accepted_gsc = test_gV
            accepted_ci  = test_ciV
            limit_flag = 2 ! limited by RuBisCo
        else
            accepted_fc   = test_fcJ
            accepted_gsc  = test_gJ
            accepted_ci   = test_ciJ
            limit_flag = 1 ! limited by light
        endif

    else  ! not resolvable
        accepted_gsc     = cuticular_gsc
        accepted_fc      = -Rdark
        accepted_ci      = can_co2

    endif

    if (quality_check) then
       write (unit=*,fmt='(80a)')         ('=',k=1,80)
       write (unit=*,fmt='(a)')           'Katul Stomatal Scheme Quality Check:'
       write (unit=*,fmt='(a,1x,i9)')   ' + HOUR:                ',current_time%hour
       write (unit=*,fmt='(a,1x,i9)')   ' + PFT:                 ',ipft
       write (unit=*,fmt='(a,1x,l9)')   ' + RESOLVABLE:          ',is_resolvable
       write (unit=*,fmt='(a,1x,es12.4)')   ' + PSI_LEAF:            ',leaf_psi
       write (unit=*,fmt='(a,1x,es12.4)')   ' + PAR:                 ',par
       write (unit=*,fmt='(a,1x,es12.4)')   ' + Vcmax25:             ',Vcmax25
       write (unit=*,fmt='(a,1x,es12.4)')   ' + Jmax25:              ',Jmax25
       write (unit=*,fmt='(a,1x,es12.4)')   ' + test_gV:             ',test_gV
       write (unit=*,fmt='(a,1x,es12.4)')   ' + test_gJ:             ',test_gJ
       write (unit=*,fmt='(a,1x,es12.4)')   ' + test_fcV:            ',test_fcV
       write (unit=*,fmt='(a,1x,es12.4)')   ' + test_fcJ:            ',test_fcJ
       write (unit=*,fmt='(a,1x,es12.4)')   ' + aero_resistance      ',aero_resistance
       write (unit=*,fmt='(a,1x,es12.4)')   ' + cuticular_gsc        ',cuticular_gsc
    endif
!-------------------- copy the solution to output-----------------------!
    A_closed        = -Rdark                ! umol/m2/s
    A_open          = accepted_fc           ! umol/m2/s
    leaf_resp       = Rdark

    gsw_closed      = cuticular_gsc / gsw_2_gsc  * mmdry / effarea_transp(ipft)  ! convert to kg/m2/s
    gsw_open        = accepted_gsc / gsw_2_gsc  * mmdry / effarea_transp(ipft)  ! convert to kg/m2/s

    !------------------- these variables are not tracked....
    blyr_cond_h2o   = leaf_gbw * mmdryi * effarea_transp(ipft)
    stom_cond_h2o   = cuticular_gsc / gsw_2_gsc
    lsfc_shv_closed = ( stom_cond_h2o * lint_shv + blyr_cond_h2o * can_shv)                 &
                    / ( stom_cond_h2o + blyr_cond_h2o)  
    stom_cond_h2o   = accepted_gsc / gsw_2_gsc
    lsfc_shv_open   = ( stom_cond_h2o * lint_shv + blyr_cond_h2o * can_shv)                 &
                    / ( stom_cond_h2o + blyr_cond_h2o)  

    blyr_cond_co2   = blyr_cond_h2o * gbw_2_gbc
    lsfc_co2_open   = can_co2 - A_open / blyr_cond_co2
    lsfc_co2_closed = can_co2 - A_closed / blyr_cond_co2

    stom_cond_co2   = cuticular_gsc
    lint_co2_closed = lsfc_co2_closed - A_closed / stom_cond_co2
    stom_cond_co2   = accepted_gsc
    lint_co2_open   = lsfc_co2_open   - A_open   / stom_cond_co2

    vmout           = Vcmax
    comppout        = cp * umol_2_mol

    return

  end subroutine katul_lphys

  !==========================================================

  
  !=======================================================================================!
  !=======================================================================================!
  ! SUBROUTINE FLUXSOLVER
  !> \brief Solve carbon flux (fc)
  !---------------------------------------------------------------------------------------!
  subroutine fluxsolver(g, ra, ca, k1, k2, k3, k4, ci,fc)
    implicit none
    
    real, intent(in) :: g, k1,k2, k3,k4, ra, ca
    real, intent(out) :: ci, fc
    real :: cip, cim, rad
    
    rad = sqrt((k1/g+k2)**2 - 4. * (k3/g + k4))
    
    cip = ca * (-(k1/g+k2) + rad)/2.
    cim = ca * (-(k1/g+k2) - rad)/2.
    
    ci = cip
    
    fc = (ca - ci) / (1./g + ra)
    
    return
  end subroutine fluxsolver
  !=======================================================================================!
  
  !=======================================================================================!
  !=======================================================================================!
  ! SUBROUTINE DERIV_FC
  !> \brief Calculate the derivation of carbon flux with respect to g
  !---------------------------------------------------------------------------------------!
  subroutine deriv_fc(g, ra, k1, k2, k3, k4, ca, ci, dfcdg)
    implicit none
    
    real, intent(in) :: g, ra, k1, k2, k3, k4, ca, ci
    real, intent(out) :: dfcdg
    real :: dcidg, myroot
    
    myroot = sqrt((k1/g+k2)**2-4.*(k3/g+k4))
    
    dcidg = ca * (0.5*k1/g**2 +   &
         0.25/myroot*(-2.*k1**2/g**3-2.*k1*k2/g**2+4.*k3/g**2))
    
    dfcdg = ((1./g+ra)*(-dcidg) - (ca-ci)*(-1./g**2))/(1./g+ra)**2
    
    return
  end subroutine deriv_fc
  !======================================================

  !=======================================================================================!
  !=======================================================================================!
  ! SUBROUTINE DERIV_DFCDG
  !> \brief Calculate the derivation of dfcdg with respect to g
  !---------------------------------------------------------------------------------------!
  subroutine deriv_dfcdg(g, ra, k1, k2, k3, k4,ca,ci,dfcdg,d2fcdg2)
    implicit none
    real, intent(in) :: g,ra,k1,k2,k3,k4,ca,ci,dfcdg
    real, intent(out) :: d2fcdg2
    real :: t1, t2

    t1 = 4. * k3 / g - 2. * k1 * (k2 + k1/g) / g
    t2 = sqrt((k1/g + k2) ** 2 - 4. * (k3/g + k4))

    d2fcdg2 = 1. / ((ra + 1./g) * g ** 2) * ( &
                ca / 2. * (t1 ** 2 / (4. * t2 ** 3) - &
                (k1 ** 2 / g ** 2 - t1) / t2 + 2. * k1 * g) + &
                2. * dfcdg - 2. * (ca - ci) / (g * (ra + 1./g)))

    return
  end subroutine deriv_dfcdg
  !======================================================

  !=======================================================================================!
  !=======================================================================================!
  ! FUNCTION HARLEY_ARRHENIUS
  !> \brief   Arrhenius equation for photosynthetic temperature dependence based
  !> on Harley et al. 1991
  !> \details This function does not consider low temperature cut-off... 
  !---------------------------------------------------------------------------------------!
  real(kind=4) function harley_arrhenius(Tleaf,Tref,Hv,Sv,Hd)
  implicit none
      real, intent(in) :: Tleaf  ! K
      real, intent(in) :: Tref   ! K
      real, intent(in) :: Hv     ! kJ/mol
      real, intent(in) :: Sv     ! kJ/mol/K
      real, intent(in) :: Hd     ! kJ/mol

      real,  parameter :: R = 8.314e-3 !kJ/mol/K

      harley_arrhenius = exp(Hv/(R * Tref)    &
             * (1 - Tref/Tleaf))    &
             / (1 + exp((Sv * Tleaf - Hd) / (R * Tleaf)))
      return

  end function harley_arrhenius
  !=======================================================================================!
  !=======================================================================================!

  !=======================================================================================!
  !=======================================================================================!
  ! FUNCTION MOD_ARRHENIUS     
  !> \brief   Arrhenius equation for photosynthetic temperature dependence
  !> with high and low temperature modification
  !---------------------------------------------------------------------------------------!
  real(kind=4) function mod_arrhenius(Tleaf,hor,Tlow,Thigh,decay_e,is_decay) 
      use physiology_coms, only : tphysrefi ! ! intent(in)
      use consts_coms,     only : lnexp_min & ! intent(in)
                                , lnexp_max & ! intent(in)
                                , t00       ! ! intent(in)
      implicit none
      !-------------------- Arguments. -----------------------------------!
      real(kind=4), intent(in) :: Tleaf   ! leaf temperature [K]
      real(kind=4), intent(in) :: hor     ! activation energy / gas constant [K]
      real(kind=4), intent(in) :: Tlow    ! low tempeature threshold [degC]
      real(kind=4), intent(in) :: Thigh   ! high tempeature threshold [degC]
      real(kind=4), intent(in) :: decay_e ! Decay rate under low/high temperature
      logical, intent(in)      :: is_decay! whether to include decay

      !-------------------- Local vars -----------------------------------!
      real(kind=4)  :: lnexp      ! term that go to the exponential
      real(kind=4)  :: lnexplow   ! term that go to the exponential
      real(kind=4)  :: lnexphigh  ! term that go to the exponential

      !------------------------------------------------------------------------------------!
      !     Find the term that goes to the exponential term, and check its size.  This is  !
      ! to avoid floating point exceptions due to overflow or underflow.                   !
      !------------------------------------------------------------------------------------!
      lnexp = hor * (tphysrefi - 1.0/Tleaf)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     If the exponential factor is tiny, make it zero, otherwise compute the actual  !
      ! function.                                                                          !
      !------------------------------------------------------------------------------------!
      if (lnexp < lnexp_min) then
         mod_arrhenius = 0.
      else
         mod_arrhenius = exp(lnexp)
      end if
      !------------------------------------------------------------------------------------!

      if (is_decay) then
          !---------------------------------------------------------------------------------!
          !    Compute the functions that will control the Vm function for low and high     !
          ! temperature.  In order to avoid floating point exceptions, we check whether the !
          ! temperature will make the exponential too large or too small.                   !
          !---------------------------------------------------------------------------------!
          !----- Low temperature. ----------------------------------------------------------!
          lnexplow  = decay_e * (Tlow  - (Tleaf - t00))
          lnexplow  = max(lnexp_min,min(lnexp_max,lnexplow))
          !----- High temperature. ---------------------------------------------------------!
          lnexphigh = decay_e * ((Tleaf-t00) - Thigh)
          lnexphigh = max(lnexp_min,min(lnexp_max,lnexphigh))
          !---------------------------------------------------------------------------------!

          mod_arrhenius = mod_arrhenius / ( (1. + exp(lnexplow)) * (1. + exp(lnexphigh)))
      endif

      return

  end function mod_arrhenius
  !=======================================================================================!
  !=======================================================================================!


  !=======================================================================================!
  !=======================================================================================!
  ! FUNCTION MOD_COLLATZ
  !> \brief Photosynthetic temperature dependence based on Collatz et al. 1991,
  !> using Q10 with high and low temperature modification
  !---------------------------------------------------------------------------------------!
  real(kind=4) function mod_collatz(temp,q10,Tlow,Thigh,decay_e,is_decay)
     use physiology_coms, only : tphysref  & ! intent(in)
                               , fcoll     ! ! intent(in)
     use consts_coms,     only : lnexp_min & ! intent(in)
                               , lnexp_max & ! intent(in)
                               , t00       ! ! intent(in)
     implicit none
     !----- Arguments. -------------------------------------------------------------------!
     real(kind=4), intent(in) :: temp      ! Temperature                           [    K]
     real(kind=4), intent(in) :: q10       ! Exponential coefficient               [    K]
     real(kind=4), intent(in) :: Tlow    ! low tempeature threshold [degC]
     real(kind=4), intent(in) :: Thigh   ! high tempeature threshold [degC]
     real(kind=4), intent(in) :: decay_e ! Decay rate under low/high temperature
     logical, intent(in)      :: is_decay! whether to include decay
     !-------------------- Local vars -----------------------------------!
     real(kind=4)  :: lnexphigh  ! term that go to the exponential
     real(kind=4)  :: lnexplow   ! term that go to the exponential

     !------------------------------------------------------------------------------------!
     !     If the exponential factor is tiny, make it zero, otherwise compute the actual  !
     ! function.                                                                          !
     !------------------------------------------------------------------------------------!
     mod_collatz = q10 ** (fcoll * (temp - tphysref))
     !------------------------------------------------------------------------------------!
     if (is_decay) then
         !---------------------------------------------------------------------------------!
         !    Compute the functions that will control the Vm function for low and high     !
         ! temperature.  In order to avoid floating point exceptions, we check whether the !
         ! temperature will make the exponential too large or too small.                   !
         !---------------------------------------------------------------------------------!
         !----- Low temperature. ----------------------------------------------------------!
         lnexplow  = decay_e * (Tlow  - (temp - t00))
         lnexplow  = max(lnexp_min,min(lnexp_max,lnexplow))
         !----- High temperature. ---------------------------------------------------------!
         lnexphigh = decay_e * ((temp-t00) - Thigh)
         lnexphigh = max(lnexp_min,min(lnexp_max,lnexphigh))
         !---------------------------------------------------------------------------------!

         mod_collatz = mod_collatz / ( (1. + exp(lnexplow)) * (1. + exp(lnexphigh)))
     endif

     return
  end function mod_collatz
  !=======================================================================================!
  !=======================================================================================!

end Module farq_katul

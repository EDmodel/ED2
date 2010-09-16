!==========================================================================================!
!    Module mem_scratch_grell - This module contains all Grell's scratch variables that do !
! not depend on ensemble dimensions. This module is shared between shallow and deep        !
! convection.                                                                              !
!==========================================================================================!
module mem_scratch_grell

   implicit none

   !------ Scalars, mostly grid definitions -----------------------------------------------!
   integer ::  mkx         & ! Number of cloud points this grid has;
              ,lpw         & ! Lower point in the vertical for this grid;
              ,kgoff       & ! Offset between BRAMS and Grell's grids
              ,kpbl        & ! PBL when running Nakanishi/Niino
              ,ktpse       ! ! maximum height allowed for cloud top
   real    ::  tscal_kf    ! ! Kain-Fritsch (1990) time scale, for ensemble.
   !---------------------------------------------------------------------------------------!



   !------ 1D dependence (mgmzp), grid-related stuff --------------------------------------!
   real, dimension(:), allocatable :: &
            z                         & ! Height above surface at model levels          [m]
           ,z_cup                     & ! Height above surface at cloud levels          [m]
           ,dzd_cld                   & ! Layer thickness for downdraft calculation     [m]
           ,dzu_cld                   ! ! Layer thickness for updraft calculation       [m]
   !---------------------------------------------------------------------------------------!


   !------ Scalars, surface variables -----------------------------------------------------!
   real    ::  co2sur    & ! Surface: CO2 mixing ratio                             [   ppm]
              ,exnersur  & ! Surface: Exner function                               [J/kg/K]
              ,psur      & ! Surface: pressure                                     [    Pa]
              ,qtotsur   & ! Surface: mixing ratio                                 [ kg/kg]
              ,qvapsur   & ! Surface: mixing ratio                                 [ kg/kg]
              ,qliqsur   & ! Surface: mixing ratio                                 [ kg/kg]
              ,qicesur   & ! Surface: mixing ratio                                 [ kg/kg]
              ,tsur      & ! Surface: temperature                                  [     K]
              ,theivsur  & ! Surface: ice-vapour equivalent potential temperature  [     K]
              ,thilsur   ! ! Surface: ice-liquid potential temperature             [     K]
   !---------------------------------------------------------------------------------------!


   !------ Scalars, variables at current time step ----------------------------------------!
   real    ::  mconv     ! ! Column integrated mass flux convergence             [ kg/m²/s]
   !---------------------------------------------------------------------------------------!

   !------ 1D dependence (mgmzp), variables with all forcings but convection --------------!
     real, dimension(:), allocatable ::          &
            dco2dt         & ! Temporary CO2 mixing ratio tendency               [   ppm/s]
           ,dqtotdt        & ! Temporary total mixing ratio tendency             [ kg/kg/s]
           ,dthildt        & ! Temporary temperature tendency                    [     K/s]
           ,dtkedt         & ! Temporary TKE tendency                            [  J/kg/s]
           ,co2            & ! CO2 Mixing ratio                                  [     ppm]
           ,exner          & ! Exner function                                    [  J/kg/K]
           ,omeg           & ! Omega - Lagrangian pressure tendency              [    Pa/s]
           ,p              & ! Pressure                                          [      Pa]
           ,rho            & ! Air density                                       [   kg/m³]
           ,qtot           & ! Total mixing ratio                                [   kg/kg]
           ,qvap           & ! Water vapour mixing ratio                         [   kg/kg]
           ,qliq           & ! Liquid water mixing ratio                         [   kg/kg]
           ,qice           & ! Ice mixing ratio                                  [   kg/kg]
           ,t              & ! Temperature                                       [       K]
           ,thil           & ! Ice-liquid potential temperature                  [       K]
           ,theiv          & ! Ice-vapour equivalent potential temperature       [       K]
           ,tke            & ! Turbulent kinetic energy                          [    J/kg]
           ,sigw           & ! Vertical velocity standard deviation              [     m/s]
           ,uwind          & ! Zonal wind speed                                  [     m/s]
           ,vwind          & ! Meridional wind speed                             [     m/s]
           ,wwind          ! ! Vertical velocity                                 [     m/s]
   !---------------------------------------------------------------------------------------!



   !------ 1D dependence (mgmzp), variables at current time step --------------------------!
   real, dimension(:), allocatable :: &
            co20      & ! CO2 mixing ratio                                       [     ppm]
           ,exner0    & ! Exner function                                         [  J/kg/K]
           ,p0        & ! Pressure with forcing                                  [     hPa]
           ,qtot0     & ! Total mixing ratio                                     [   kg/kg]
           ,qvap0     & ! Water vapour mixing ratio                              [   kg/kg]
           ,qliq0     & ! Liquid water mixing ratio                              [   kg/kg]
           ,qice0     & ! Ice mixing ratio                                       [   kg/kg]
           ,rho0      & ! Air density                                            [   kg/m³]
           ,t0        & ! Temperature                                            [       K]
           ,thil0     & ! Ice-liquid potential temperature                       [       K]
           ,theiv0    & ! Ice-vapour equivalent potential temperature            [       K]
           ,tke0      ! ! Turbulent Kinetic Energy                               [    J/kg]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Environment variables at Grell's staggered grid.                                  !
   !---------------------------------------------------------------------------------------!
   !----- Variables based on values before any forcing is applied -------------------------!
   real, dimension(:), allocatable ::  &
            co20_cup      & ! CO2 mixing ratio                                   [     ppm]
           ,exner0_cup    & ! Exner function                                     [    J/kg]
           ,p0_cup        & ! Pressure                                           [      Pa]
           ,qtot0_cup     & ! Total mixing ratio                                 [   kg/kg]
           ,qvap0_cup     & ! Water vapour mixing ratio                          [   kg/kg]
           ,qliq0_cup     & ! Liquid water                                       [   kg/kg]
           ,qice0_cup     & ! Mixing ratio                                       [   kg/kg]
           ,qsat0_cup     & ! Saturation mixing ratio                            [   kg/kg]
           ,rho0_cup      & ! Density                                            [   kg/m³]
           ,t0_cup        & ! Temperature                                        [       K]
           ,thil0_cup     & ! Ice-liquid potential temperature                   [       K]
           ,theiv0_cup    & ! Ice-vapour potential temperature                   [       K]
           ,theivs0_cup   ! ! Sat. Ice-vapour equiv. pot. temp.                  [       K]
   !----- Variables with all time tendencies but this convection. -------------------------!
   real, dimension(:), allocatable :: &
            co2_cup       & ! CO2 mixing ratio                                   [     ppm]
           ,exner_cup     & ! Exner function                                     [    J/kg]
           ,p_cup         & ! Pressure                                           [      Pa]
           ,qtot_cup      & ! Total mixing ratio                                 [   kg/kg]
           ,qvap_cup      & ! Water vapour mixing ratio                          [   kg/kg]
           ,qliq_cup      & ! Liquid water                                       [   kg/kg]
           ,qice_cup      & ! Mixing ratio                                       [   kg/kg]
           ,qsat_cup      & ! Saturation mixing ratio                            [   kg/kg]
           ,rho_cup       & ! Density                                            [   kg/m³]
           ,t_cup         & ! Temperature                                        [       K]
           ,thil_cup      & ! Ice-liquid potential temperature                   [       K]
           ,theiv_cup     & ! Ice-vapour equiv. pot. temp.                       [       K]
           ,theivs_cup    ! ! Sat. Ice-vapour equiv. pot. temp.                  [       K]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Downdraft-related variables                                                        !
   !---------------------------------------------------------------------------------------!
   !----- Based on current values (no forcing) --------------------------------------------!
   real, dimension(:), allocatable :: &
            co20d_cld     & ! CO2 water mixing ratio                             [     ppm]
           ,dby0d         & ! Buoyancy acceleration                              [    m/s²]
           ,pw0d_cld      & ! Condensation                                       [   kg/kg]
           ,qtot0d_cld    & ! Total water mixing ratio                           [   kg/kg]
           ,qvap0d_cld    & ! Vapour mixing ratio                                [   kg/kg]
           ,qliq0d_cld    & ! Liquid water mixing ratio                          [   kg/kg]
           ,qice0d_cld    & ! Ice mixing ratio                                   [   kg/kg]
           ,qsat0d_cld    & ! Sat. mixing ratio                                  [   kg/kg]
           ,rho0d_cld     & ! Density                                            [   kg/m³]
           ,t0d_cld       & ! Temperature                                        [       K]
           ,theiv0d_cld   & ! Ice-vapour equiv. pot. temp.                       [       K]
           ,thil0d_cld    ! ! Ice-liquid pot. temperature                        [       K]
   !----- Based on values with forcing ----------------------------------------------------!
   real, dimension(:), allocatable :: &
            cdd           & ! Norm. dndraft detrainment rate                     [     ---]
           ,co2d_cld      & ! CO2 water mixing ratio                             [     ppm]
           ,dbyd          & ! Buoyancy acceleration                              [    m/s²]
           ,etad_cld      & ! normalised downdraft mass flux                     [     ---]
           ,mentrd_rate   & ! Norm. dndraft entrainment rate                     [     ---]
           ,pwd_cld       & ! Condensation                                       [   kg/kg]
           ,qtotd_cld     & ! Total water mixing ratio                           [   kg/kg]
           ,qvapd_cld     & ! Vapour mixing ratio                                [   kg/kg]
           ,qliqd_cld     & ! Liquid water mixing ratio                          [   kg/kg]
           ,qiced_cld     & ! Ice mixing ratio                                   [   kg/kg]
           ,qsatd_cld     & ! Sat. mixing ratio                                  [   kg/kg]
           ,rhod_cld      & ! Density                                            [   kg/m³]
           ,td_cld        & ! Temperature                                        [       K]
           ,theivd_cld    & ! Ice-vapour equiv. pot. temp.                       [       K]
           ,thild_cld     ! ! Ice-liquid pot. temperature                        [       K]
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Updraft-related variables                                                          !
   !---------------------------------------------------------------------------------------!
   !----- Based on current values (no forcing) --------------------------------------------!
   real, dimension(:), allocatable :: &
            co20u_cld     & ! CO2 water mixing ratio                             [     ppm]
           ,dby0u         & ! Buoyancy acceleration                              [    m/s²]
           ,pw0u_cld      & ! Condensation                                       [   kg/kg]
           ,qtot0u_cld    & ! Total water mixing ratio                           [   kg/kg]
           ,qvap0u_cld    & ! Vapour mixing ratio                                [   kg/kg]
           ,qliq0u_cld    & ! Liquid water mixing ratio                          [   kg/kg]
           ,qice0u_cld    & ! Ice mixing ratio                                   [   kg/kg]
           ,qsat0u_cld    & ! Sat. mixing ratio                                  [   kg/kg]
           ,rho0u_cld     & ! Density                                            [   kg/m³]
           ,t0u_cld       & ! Temperature                                        [       K]
           ,thil0u_cld    & ! Ice-liquid potential temperature                   [       K]
           ,theiv0u_cld   ! ! Ice-vapour equiv. pot. temp.                       [       K]
   !----- Based on values with forcing ----------------------------------------------------!
   real, dimension(:), allocatable :: &
            cdu           & ! Norm. updraft detrainment rate                     [     ---]
           ,co2u_cld      & ! CO2 water mixing ratio                             [     ppm]
           ,dbyu          & ! Buoyancy acceleration                              [    m/s²]
           ,etau_cld      & ! Normalised updraft mass flux                       [     ---]
           ,mentru_rate   & ! Norm. updraft entrainment rate                     [     ---]
           ,pwu_cld       & ! Condensation                                       [   kg/kg]
           ,qtotu_cld     & ! Total water mixing ratio                           [   kg/kg]
           ,qvapu_cld     & ! Vapour mixing ratio                                [   kg/kg]
           ,qliqu_cld     & ! Liquid water mixing ratio                          [   kg/kg]
           ,qiceu_cld     & ! Ice mixing ratio                                   [   kg/kg]
           ,qsatu_cld     & ! Sat. mixing ratio                                  [   kg/kg]
           ,rhou_cld      & ! Density                                            [   kg/m³]
           ,tu_cld        & ! Temperature                                        [       K]
           ,thilu_cld     & ! Ice-liquid potential temperature                   [       K]
           ,theivu_cld    ! ! Ice-vapour equiv. pot. temp.                       [       K]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Variables modified by the arbitrary mass flux of the subensemble                   !
   !---------------------------------------------------------------------------------------!
   !----- Model levels --------------------------------------------------------------------!
   real, dimension(:), allocatable :: &
            x_co2        & ! CO2 mixing ratio                                     [    ppm]
           ,x_qtot       & ! Total mixing ratio                                   [  kg/kg]
           ,x_qvap       & ! Vapour mixing ratio                                  [  kg/kg]
           ,x_qliq       & ! Liquid water mixing ratio                            [  kg/kg]
           ,x_qice       & ! Ice mixing ratio                                     [  kg/kg]
           ,x_t          & ! Temperature                                          [      K]
           ,x_theiv      & ! Ice-vapour equiv. pot. temp.                         [      K]
           ,x_thil       ! ! Ice-liquid potential temperature.                    [      K]
   !----- Cloud levels --------------------------------------------------------------------!
   real, dimension(:), allocatable :: &
            x_co2_cup    & ! CO2 mixing ratio                                     [    ppm]
           ,x_exner_cup  & ! Exner function                                       [   J/kg]
           ,x_p_cup      & ! Pressure                                             [     Pa]
           ,x_qtot_cup   & ! Total mixing ratio                                   [  kg/kg]
           ,x_qvap_cup   & ! Water vapour mixing ratio                            [  kg/kg]
           ,x_qliq_cup   & ! Liquid water                                         [  kg/kg]
           ,x_qice_cup   & ! Mixing ratio                                         [  kg/kg]
           ,x_qsat_cup   & ! Saturation mixing ratio                              [  kg/kg]
           ,x_rho_cup    & ! Density                                              [  kg/m³]
           ,x_t_cup      & ! Temperature                                          [      K]
           ,x_thil_cup   & ! Ice-liquid potential temperature                     [      K]
           ,x_theiv_cup  & ! Ice-vapour equiv. pot. temp.                         [      K]
           ,x_theivs_cup ! ! Sat. Ice-vapour equiv. pot. temp.                    [      K]
   !----- Downdraft variables -------------------------------------------------------------!
   real, dimension(:), allocatable :: &
            x_co2d_cld   & ! CO2 mixing ratio                                     [    ppm]
           ,x_dbyd       & ! Buoyancy acceleration                                [   m/s²]
           ,x_pwd_cld    & ! Condensation                                         [  kg/kg]
           ,x_qtotd_cld  & ! Total water mixing ratio                             [  kg/kg]
           ,x_qvapd_cld  & ! Vapour mixing ratio                                  [  kg/kg]
           ,x_qliqd_cld  & ! Liquid water mixing ratio                            [  kg/kg]
           ,x_qiced_cld  & ! Ice mixing ratio                                     [  kg/kg]
           ,x_qsatd_cld  & ! Sat. mixing ratio                                    [  kg/kg]
           ,x_rhod_cld   & ! Density                                              [  kg/m³]
           ,x_td_cld     & ! Temperature                                          [      K]
           ,x_theivd_cld & ! Ice-vapour equiv. pot. temp.                         [      K]
           ,x_thild_cld  ! ! Ice-liquid potential temperature                     [      K]
   !----- Updraft variables ---------------------------------------------------------------!
   real, dimension(:), allocatable :: &
            x_co2u_cld   & ! CO2 mixing ratio                                     [    ppm]
           ,x_dbyu       & ! Buoyancy acceleration                                [   m/s²]
           ,x_pwu_cld    & ! Condensation                                         [  kg/kg]
           ,x_qtotu_cld  & ! Total water mixing ratio                             [  kg/kg]
           ,x_qvapu_cld  & ! Vapour mixing ratio                                  [  kg/kg]
           ,x_qliqu_cld  & ! Liquid water mixing ratio                            [  kg/kg]
           ,x_qiceu_cld  & ! Ice mixing ratio                                     [  kg/kg]
           ,x_qsatu_cld  & ! Sat. mixing ratio                                    [  kg/kg]
           ,x_rhou_cld   & ! Density                                              [  kg/m³]
           ,x_tu_cld     & ! Temperature                                          [      K]
           ,x_thilu_cld  & ! Ice-liquid potential temperature                     [      K]
           ,x_theivu_cld ! ! Ice-vapour equiv. pot. temp.                         [      K]
   !---------------------------------------------------------------------------------------!



   !----- Scalars. ------------------------------------------------------------------------!
   logical :: comp_dn   ! Flag for downdraft thermodynamics.                     [     ---]
   integer :: ierr      ! Flag for convection error                              [     ---]
   integer :: klod      ! Level of origin of downdrafts                          [     ---]
   integer :: klou      ! Level of origin of updrafts                            [     ---]
   integer :: klcl      ! Lifting condensation level                             [     ---]
   integer :: klfc      ! Level of free convection                               [     ---]
   integer :: kdet      ! Top of downdraft detrainemnt layer                     [     ---]
   integer :: kstabi    ! cloud stable layer base                                [     ---]
   integer :: kstabm    ! cloud stable layer top                                 [     ---]
   integer :: klnb      ! cloud level of neutral buoyancy                        [     ---]
   integer :: ktop      ! cloud top                                              [     ---]
   real    :: aa0u      ! ! Updraft work function                                [    J/kg]
   real    :: aa0d      ! ! Updraft work function                                [    J/kg]
   real    :: pwav0     ! ! Integrated condensation                              [   kg/kg]
   real    :: pwav      ! ! Integrated condensation with forcing                 [   kg/kg]
   real    :: pwev0     ! ! Integrated evaporation                               [   kg/kg]
   real    :: pwev      ! ! Integrated evaporation  with forcing                 [   kg/kg]
   real    :: x_aau     ! ! Updraft work function                                [    J/kg]
   real    :: x_aad     ! ! Updraft work function                                [    J/kg]
   real    :: x_pwav    ! ! Integrated condensation                              [   kg/kg]
   real    :: x_pwev    ! ! Integrated evaporation                               [   kg/kg]
   real    :: wbuoymin0 ! ! Minimum buoyant velocity                             [     m/s]
   real    :: wbuoymin  ! Minimum buoyancy velocity with forcing                 [     m/s]
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_scratch_grell(mgmzp)

      implicit none
      integer, intent(in)                    :: mgmzp

      allocate (z             (mgmzp))
      allocate (z_cup         (mgmzp))
      allocate (dzd_cld       (mgmzp))
      allocate (dzu_cld       (mgmzp))

      allocate (dco2dt        (mgmzp))
      allocate (dqtotdt       (mgmzp))
      allocate (dthildt       (mgmzp))
      allocate (dtkedt        (mgmzp))
      allocate (co2           (mgmzp))
      allocate (exner         (mgmzp))
      allocate (omeg          (mgmzp))
      allocate (p             (mgmzp))
      allocate (rho           (mgmzp))
      allocate (qtot          (mgmzp))
      allocate (qvap          (mgmzp))
      allocate (qliq          (mgmzp))
      allocate (qice          (mgmzp))
      allocate (t             (mgmzp))
      allocate (thil          (mgmzp))
      allocate (theiv         (mgmzp))
      allocate (tke           (mgmzp))
      allocate (sigw          (mgmzp))
      allocate (uwind         (mgmzp))
      allocate (vwind         (mgmzp))
      allocate (wwind         (mgmzp))

      allocate (co20          (mgmzp))
      allocate (exner0        (mgmzp))
      allocate (p0            (mgmzp))
      allocate (qtot0         (mgmzp))
      allocate (qvap0         (mgmzp))
      allocate (qliq0         (mgmzp))
      allocate (qice0         (mgmzp))
      allocate (rho0          (mgmzp))
      allocate (t0            (mgmzp))
      allocate (thil0         (mgmzp))
      allocate (theiv0        (mgmzp))
      allocate (tke0          (mgmzp))

      allocate (co20_cup      (mgmzp))
      allocate (exner0_cup    (mgmzp))
      allocate (p0_cup        (mgmzp))
      allocate (qtot0_cup     (mgmzp))
      allocate (qvap0_cup     (mgmzp))
      allocate (qliq0_cup     (mgmzp))
      allocate (qice0_cup     (mgmzp))
      allocate (qsat0_cup     (mgmzp))
      allocate (rho0_cup      (mgmzp))
      allocate (t0_cup        (mgmzp))
      allocate (thil0_cup     (mgmzp))
      allocate (theiv0_cup    (mgmzp))
      allocate (theivs0_cup   (mgmzp))

      allocate (co2_cup       (mgmzp))
      allocate (exner_cup     (mgmzp))
      allocate (p_cup         (mgmzp))
      allocate (qtot_cup      (mgmzp))
      allocate (qvap_cup      (mgmzp))
      allocate (qliq_cup      (mgmzp))
      allocate (qice_cup      (mgmzp))
      allocate (qsat_cup      (mgmzp))
      allocate (rho_cup       (mgmzp))
      allocate (t_cup         (mgmzp))
      allocate (thil_cup      (mgmzp))
      allocate (theiv_cup     (mgmzp))
      allocate (theivs_cup    (mgmzp))

      allocate (co20d_cld     (mgmzp))
      allocate (dby0d         (mgmzp))
      allocate (pw0d_cld      (mgmzp))
      allocate (qtot0d_cld    (mgmzp))
      allocate (qvap0d_cld    (mgmzp))
      allocate (qliq0d_cld    (mgmzp))
      allocate (qice0d_cld    (mgmzp))
      allocate (qsat0d_cld    (mgmzp))
      allocate (rho0d_cld     (mgmzp))
      allocate (t0d_cld       (mgmzp))
      allocate (theiv0d_cld   (mgmzp))
      allocate (thil0d_cld    (mgmzp))

      allocate (cdd           (mgmzp))
      allocate (co2d_cld      (mgmzp))
      allocate (dbyd          (mgmzp))
      allocate (etad_cld      (mgmzp))
      allocate (mentrd_rate   (mgmzp))
      allocate (pwd_cld       (mgmzp))
      allocate (qtotd_cld     (mgmzp))
      allocate (qvapd_cld     (mgmzp))
      allocate (qliqd_cld     (mgmzp))
      allocate (qiced_cld     (mgmzp))
      allocate (qsatd_cld     (mgmzp))
      allocate (rhod_cld      (mgmzp))
      allocate (td_cld        (mgmzp))
      allocate (theivd_cld    (mgmzp))
      allocate (thild_cld     (mgmzp))

      allocate (co20u_cld     (mgmzp))
      allocate (dby0u         (mgmzp))
      allocate (pw0u_cld      (mgmzp))
      allocate (qtot0u_cld    (mgmzp))
      allocate (qvap0u_cld    (mgmzp))
      allocate (qliq0u_cld    (mgmzp))
      allocate (qice0u_cld    (mgmzp))
      allocate (qsat0u_cld    (mgmzp))
      allocate (rho0u_cld     (mgmzp))
      allocate (t0u_cld       (mgmzp))
      allocate (thil0u_cld    (mgmzp))
      allocate (theiv0u_cld   (mgmzp))

      allocate (cdu           (mgmzp))
      allocate (co2u_cld      (mgmzp))
      allocate (dbyu          (mgmzp))
      allocate (etau_cld      (mgmzp))
      allocate (mentru_rate   (mgmzp))
      allocate (pwu_cld       (mgmzp))
      allocate (qtotu_cld     (mgmzp))
      allocate (qvapu_cld     (mgmzp))
      allocate (qliqu_cld     (mgmzp))
      allocate (qiceu_cld     (mgmzp))
      allocate (qsatu_cld     (mgmzp))
      allocate (rhou_cld      (mgmzp))
      allocate (tu_cld        (mgmzp))
      allocate (thilu_cld     (mgmzp))
      allocate (theivu_cld    (mgmzp))

      allocate (x_co2         (mgmzp))
      allocate (x_qtot        (mgmzp))
      allocate (x_qvap        (mgmzp))
      allocate (x_qliq        (mgmzp))
      allocate (x_qice        (mgmzp))
      allocate (x_t           (mgmzp))
      allocate (x_theiv       (mgmzp))
      allocate (x_thil        (mgmzp))

      allocate (x_co2_cup     (mgmzp))
      allocate (x_exner_cup   (mgmzp))
      allocate (x_p_cup       (mgmzp))
      allocate (x_qtot_cup    (mgmzp))
      allocate (x_qvap_cup    (mgmzp))
      allocate (x_qliq_cup    (mgmzp))
      allocate (x_qice_cup    (mgmzp))
      allocate (x_qsat_cup    (mgmzp))
      allocate (x_rho_cup     (mgmzp))
      allocate (x_t_cup       (mgmzp))
      allocate (x_thil_cup    (mgmzp))
      allocate (x_theiv_cup   (mgmzp))
      allocate (x_theivs_cup  (mgmzp))

      allocate (x_co2d_cld    (mgmzp))
      allocate (x_dbyd        (mgmzp))
      allocate (x_pwd_cld     (mgmzp))
      allocate (x_qtotd_cld   (mgmzp))
      allocate (x_qvapd_cld   (mgmzp))
      allocate (x_qliqd_cld   (mgmzp))
      allocate (x_qiced_cld   (mgmzp))
      allocate (x_qsatd_cld   (mgmzp))
      allocate (x_rhod_cld    (mgmzp))
      allocate (x_td_cld      (mgmzp))
      allocate (x_theivd_cld  (mgmzp))
      allocate (x_thild_cld   (mgmzp))

      allocate (x_co2u_cld    (mgmzp))
      allocate (x_dbyu        (mgmzp))
      allocate (x_pwu_cld     (mgmzp))
      allocate (x_qtotu_cld   (mgmzp))
      allocate (x_qvapu_cld   (mgmzp))
      allocate (x_qliqu_cld   (mgmzp))
      allocate (x_qiceu_cld   (mgmzp))
      allocate (x_qsatu_cld   (mgmzp))
      allocate (x_rhou_cld    (mgmzp))
      allocate (x_tu_cld      (mgmzp))
      allocate (x_thilu_cld   (mgmzp))
      allocate (x_theivu_cld  (mgmzp))

      return
   end subroutine alloc_scratch_grell
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine dealloc_scratch_grell()

      implicit none
     
      if(allocated(z             )) deallocate (z             )
      if(allocated(z_cup         )) deallocate (z_cup         )
      if(allocated(dzd_cld       )) deallocate (dzd_cld       )
      if(allocated(dzu_cld       )) deallocate (dzu_cld       )

      if(allocated(dco2dt        )) deallocate (dco2dt        )
      if(allocated(dqtotdt       )) deallocate (dqtotdt       )
      if(allocated(dthildt       )) deallocate (dthildt       )
      if(allocated(dtkedt        )) deallocate (dtkedt        )
      if(allocated(co2           )) deallocate (co2           )
      if(allocated(exner         )) deallocate (exner         )
      if(allocated(omeg          )) deallocate (omeg          )
      if(allocated(p             )) deallocate (p             )
      if(allocated(qtot          )) deallocate (qtot          )
      if(allocated(qvap          )) deallocate (qvap          )
      if(allocated(qliq          )) deallocate (qliq          )
      if(allocated(qice          )) deallocate (qice          )
      if(allocated(rho           )) deallocate (rho           )
      if(allocated(t             )) deallocate (t             )
      if(allocated(thil          )) deallocate (thil          )
      if(allocated(theiv         )) deallocate (theiv         )
      if(allocated(tke           )) deallocate (tke           )
      if(allocated(sigw          )) deallocate (sigw          )
      if(allocated(uwind         )) deallocate (uwind         )
      if(allocated(vwind         )) deallocate (vwind         )
      if(allocated(wwind         )) deallocate (wwind         )

      if(allocated(co20          )) deallocate (co20          )
      if(allocated(exner0        )) deallocate (exner0        )
      if(allocated(p0            )) deallocate (p0            )
      if(allocated(qtot0         )) deallocate (qtot0         )
      if(allocated(qvap0         )) deallocate (qvap0         )
      if(allocated(qliq0         )) deallocate (qliq0         )
      if(allocated(qice0         )) deallocate (qice0         )
      if(allocated(rho0          )) deallocate (rho0          )
      if(allocated(t0            )) deallocate (t0            )
      if(allocated(thil0         )) deallocate (thil0         )
      if(allocated(theiv0        )) deallocate (theiv0        )
      if(allocated(tke0          )) deallocate (tke0          )

      if(allocated(co20_cup      )) deallocate (co20_cup      )
      if(allocated(exner0_cup    )) deallocate (exner0_cup    )
      if(allocated(p0_cup        )) deallocate (p0_cup        )
      if(allocated(qtot0_cup     )) deallocate (qtot0_cup     )
      if(allocated(qvap0_cup     )) deallocate (qvap0_cup     )
      if(allocated(qliq0_cup     )) deallocate (qliq0_cup     )
      if(allocated(qice0_cup     )) deallocate (qice0_cup     )
      if(allocated(qsat0_cup     )) deallocate (qsat0_cup     )
      if(allocated(rho0_cup      )) deallocate (rho0_cup      )
      if(allocated(t0_cup        )) deallocate (t0_cup        )
      if(allocated(thil0_cup     )) deallocate (thil0_cup     )
      if(allocated(theiv0_cup    )) deallocate (theiv0_cup    )
      if(allocated(theivs0_cup   )) deallocate (theivs0_cup   )

      if(allocated(co2_cup       )) deallocate (co2_cup       )
      if(allocated(exner_cup     )) deallocate (exner_cup     )
      if(allocated(p_cup         )) deallocate (p_cup         )
      if(allocated(qtot_cup      )) deallocate (qtot_cup      )
      if(allocated(qvap_cup      )) deallocate (qvap_cup      )
      if(allocated(qliq_cup      )) deallocate (qliq_cup      )
      if(allocated(qice_cup      )) deallocate (qice_cup      )
      if(allocated(qsat_cup      )) deallocate (qsat_cup      )
      if(allocated(rho_cup       )) deallocate (rho_cup       )
      if(allocated(t_cup         )) deallocate (t_cup         )
      if(allocated(thil_cup      )) deallocate (thil_cup      )
      if(allocated(theiv_cup     )) deallocate (theiv_cup     )
      if(allocated(theivs_cup    )) deallocate (theivs_cup    )

      if(allocated(co20d_cld     )) deallocate (co20d_cld     )
      if(allocated(dby0d         )) deallocate (dby0d         )
      if(allocated(pw0d_cld      )) deallocate (pw0d_cld      )
      if(allocated(qtot0d_cld    )) deallocate (qtot0d_cld    )
      if(allocated(qvap0d_cld    )) deallocate (qvap0d_cld    )
      if(allocated(qliq0d_cld    )) deallocate (qliq0d_cld    )
      if(allocated(qice0d_cld    )) deallocate (qice0d_cld    )
      if(allocated(qsat0d_cld    )) deallocate (qsat0d_cld    )
      if(allocated(rho0d_cld     )) deallocate (rho0d_cld     )
      if(allocated(t0d_cld       )) deallocate (t0d_cld       )
      if(allocated(theiv0d_cld   )) deallocate (theiv0d_cld   )
      if(allocated(thil0d_cld    )) deallocate (thil0d_cld    )

      if(allocated(cdd           )) deallocate (cdd           )
      if(allocated(co2d_cld      )) deallocate (co2d_cld      )
      if(allocated(dbyd          )) deallocate (dbyd          )
      if(allocated(etad_cld      )) deallocate (etad_cld      )
      if(allocated(mentrd_rate   )) deallocate (mentrd_rate   )
      if(allocated(pwd_cld       )) deallocate (pwd_cld       )
      if(allocated(qtotd_cld     )) deallocate (qtotd_cld     )
      if(allocated(qvapd_cld     )) deallocate (qvapd_cld     )
      if(allocated(qliqd_cld     )) deallocate (qliqd_cld     )
      if(allocated(qiced_cld     )) deallocate (qiced_cld     )
      if(allocated(qsatd_cld     )) deallocate (qsatd_cld     )
      if(allocated(rhod_cld      )) deallocate (rhod_cld      )
      if(allocated(td_cld        )) deallocate (td_cld        )
      if(allocated(theivd_cld    )) deallocate (theivd_cld    )
      if(allocated(thild_cld     )) deallocate (thild_cld     )

      if(allocated(co20u_cld     )) deallocate (co20u_cld     )
      if(allocated(dby0u         )) deallocate (dby0u         )
      if(allocated(pw0u_cld      )) deallocate (pw0u_cld      )
      if(allocated(qtot0u_cld    )) deallocate (qtot0u_cld    )
      if(allocated(qvap0u_cld    )) deallocate (qvap0u_cld    )
      if(allocated(qliq0u_cld    )) deallocate (qliq0u_cld    )
      if(allocated(qice0u_cld    )) deallocate (qice0u_cld    )
      if(allocated(qsat0u_cld    )) deallocate (qsat0u_cld    )
      if(allocated(rho0u_cld     )) deallocate (rho0u_cld     )
      if(allocated(t0u_cld       )) deallocate (t0u_cld       )
      if(allocated(thil0u_cld    )) deallocate (thil0u_cld    )
      if(allocated(theiv0u_cld   )) deallocate (theiv0u_cld   )

      if(allocated(cdu           )) deallocate (cdu           )
      if(allocated(co2u_cld      )) deallocate (co2u_cld      )
      if(allocated(dbyu          )) deallocate (dbyu          )
      if(allocated(etau_cld      )) deallocate (etau_cld      )
      if(allocated(mentru_rate   )) deallocate (mentru_rate   )
      if(allocated(pwu_cld       )) deallocate (pwu_cld       )
      if(allocated(qtotu_cld     )) deallocate (qtotu_cld     )
      if(allocated(qvapu_cld     )) deallocate (qvapu_cld     )
      if(allocated(qliqu_cld     )) deallocate (qliqu_cld     )
      if(allocated(qiceu_cld     )) deallocate (qiceu_cld     )
      if(allocated(qsatu_cld     )) deallocate (qsatu_cld     )
      if(allocated(rhou_cld      )) deallocate (rhou_cld      )
      if(allocated(tu_cld        )) deallocate (tu_cld        )
      if(allocated(thilu_cld     )) deallocate (thilu_cld     )
      if(allocated(theivu_cld    )) deallocate (theivu_cld    )

      if(allocated(x_co2         )) deallocate (x_co2         )
      if(allocated(x_qtot        )) deallocate (x_qtot        )
      if(allocated(x_qvap        )) deallocate (x_qvap        )
      if(allocated(x_qliq        )) deallocate (x_qliq        )
      if(allocated(x_qice        )) deallocate (x_qice        )
      if(allocated(x_t           )) deallocate (x_t           )
      if(allocated(x_theiv       )) deallocate (x_theiv       )
      if(allocated(x_thil        )) deallocate (x_thil        )

      if(allocated(x_co2_cup     )) deallocate (x_co2_cup     )
      if(allocated(x_exner_cup   )) deallocate (x_exner_cup   )
      if(allocated(x_p_cup       )) deallocate (x_p_cup       )
      if(allocated(x_qtot_cup    )) deallocate (x_qtot_cup    )
      if(allocated(x_qvap_cup    )) deallocate (x_qvap_cup    )
      if(allocated(x_qliq_cup    )) deallocate (x_qliq_cup    )
      if(allocated(x_qice_cup    )) deallocate (x_qice_cup    )
      if(allocated(x_qsat_cup    )) deallocate (x_qsat_cup    )
      if(allocated(x_rho_cup     )) deallocate (x_rho_cup     )
      if(allocated(x_t_cup       )) deallocate (x_t_cup       )
      if(allocated(x_thil_cup    )) deallocate (x_thil_cup    )
      if(allocated(x_theiv_cup   )) deallocate (x_theiv_cup   )
      if(allocated(x_theivs_cup  )) deallocate (x_theivs_cup  )

      if(allocated(x_co2d_cld    )) deallocate (x_co2d_cld    )
      if(allocated(x_dbyd        )) deallocate (x_dbyd        )
      if(allocated(x_pwd_cld     )) deallocate (x_pwd_cld     )
      if(allocated(x_qtotd_cld   )) deallocate (x_qtotd_cld   )
      if(allocated(x_qvapd_cld   )) deallocate (x_qvapd_cld   )
      if(allocated(x_qliqd_cld   )) deallocate (x_qliqd_cld   )
      if(allocated(x_qiced_cld   )) deallocate (x_qiced_cld   )
      if(allocated(x_qsatd_cld   )) deallocate (x_qsatd_cld   )
      if(allocated(x_rhod_cld    )) deallocate (x_rhod_cld    )
      if(allocated(x_td_cld      )) deallocate (x_td_cld      )
      if(allocated(x_theivd_cld  )) deallocate (x_theivd_cld  )
      if(allocated(x_thild_cld   )) deallocate (x_thild_cld   )

      if(allocated(x_co2u_cld    )) deallocate (x_co2u_cld    )
      if(allocated(x_dbyu        )) deallocate (x_dbyu        )
      if(allocated(x_pwu_cld     )) deallocate (x_pwu_cld     )
      if(allocated(x_qtotu_cld   )) deallocate (x_qtotu_cld   )
      if(allocated(x_qvapu_cld   )) deallocate (x_qvapu_cld   )
      if(allocated(x_qliqu_cld   )) deallocate (x_qliqu_cld   )
      if(allocated(x_qiceu_cld   )) deallocate (x_qiceu_cld   )
      if(allocated(x_qsatu_cld   )) deallocate (x_qsatu_cld   )
      if(allocated(x_rhou_cld    )) deallocate (x_rhou_cld    )
      if(allocated(x_tu_cld      )) deallocate (x_tu_cld      )
      if(allocated(x_thilu_cld   )) deallocate (x_thilu_cld   )
      if(allocated(x_theivu_cld  )) deallocate (x_theivu_cld  )

      return
   end subroutine dealloc_scratch_grell
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   recursive subroutine zero_scratch_grell(fullz)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer, intent(in) :: fullz ! Indicates where I should stop resetting. 
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Here we decide what to reset based on the given flag.                          !
      !------------------------------------------------------------------------------------!
      select case (fullz)
      case (0)
         !---------------------------------------------------------------------------------!
         !     Here we reset all grid and large-scale dependent variables.                 !
         !---------------------------------------------------------------------------------!
         if(allocated(z             )) z             = 0.
         if(allocated(z_cup         )) z_cup         = 0.
         if(allocated(dzd_cld       )) dzd_cld       = 0.
         if(allocated(dzu_cld       )) dzu_cld       = 0.

         if(allocated(dco2dt        )) dco2dt        = 0.
         if(allocated(dqtotdt       )) dqtotdt       = 0.
         if(allocated(dthildt       )) dthildt       = 0.
         if(allocated(dtkedt        )) dtkedt        = 0.

         if(allocated(co2           )) co2           = 0.
         if(allocated(exner         )) exner         = 0.
         if(allocated(omeg          )) omeg          = 0.
         if(allocated(p             )) p             = 0.
         if(allocated(qtot          )) qtot          = 0.
         if(allocated(qvap          )) qvap          = 0.
         if(allocated(qliq          )) qliq          = 0.
         if(allocated(qice          )) qice          = 0.
         if(allocated(rho           )) rho           = 0.
         if(allocated(t             )) t             = 0.
         if(allocated(thil          )) thil          = 0.
         if(allocated(theiv         )) theiv         = 0.
         if(allocated(tke           )) tke           = 0.
         if(allocated(sigw          )) sigw          = 0.
         if(allocated(wwind         )) wwind         = 0.

         if(allocated(co20          )) co20          = 0.
         if(allocated(exner0        )) exner0        = 0.
         if(allocated(p0            )) p0            = 0.
         if(allocated(qtot0         )) qtot0         = 0.
         if(allocated(qvap0         )) qvap0         = 0.
         if(allocated(qliq0         )) qliq0         = 0.
         if(allocated(qice0         )) qice0         = 0.
         if(allocated(rho0          )) rho0          = 0.
         if(allocated(t0            )) t0            = 0.
         if(allocated(thil0         )) thil0         = 0.
         if(allocated(theiv0        )) theiv0        = 0.
         if(allocated(tke0          )) tke0          = 0.

         !---------------------------------------------------------------------------------!
         ! Flushing scalars, we don't need to check for allocation here...                 !
         !---------------------------------------------------------------------------------!
         !----- Integer variables ---------------------------------------------------------!
         mkx               = 0
         lpw               = 0
         kgoff             = 0
         kpbl              = 0
         ktpse             = 0
         !----- Real variables. -----------------------------------------------------------!
         tscal_kf          = 0.
         mconv             = 0.
         co2sur            = 0.
         exnersur          = 0.
         psur              = 0.
         qtotsur           = 0.
         qvapsur           = 0.
         qliqsur           = 0.
         qicesur           = 0.
         tsur              = 0.
         theivsur          = 0.
         thilsur           = 0.
         !---------------------------------------------------------------------------------!

      case(1)
         !---------------------------------------------------------------------------------!
         !      Although these variables aren't really static control, they must be reset  !
         ! here because they are recalculated inside the cloud loop.
         !---------------------------------------------------------------------------------!
         if(allocated(uwind         )) uwind         = 0.
         if(allocated(vwind         )) vwind         = 0.

         !---------------------------------------------------------------------------------!
         !      Here we reset all static control variables.  It is also reset at the       !
         ! dynamic control because we "borrow" some variables to be used there as aliases  !
         ! to those long ensemble variable names.                                          !
         !---------------------------------------------------------------------------------!
         if(allocated(co20_cup      )) co20_cup      = 0.
         if(allocated(exner0_cup    )) exner0_cup    = 0.
         if(allocated(p0_cup        )) p0_cup        = 0.
         if(allocated(qtot0_cup     )) qtot0_cup     = 0.
         if(allocated(qvap0_cup     )) qvap0_cup     = 0.
         if(allocated(qliq0_cup     )) qliq0_cup     = 0.
         if(allocated(qice0_cup     )) qice0_cup     = 0.
         if(allocated(qsat0_cup     )) qsat0_cup     = 0.
         if(allocated(rho0_cup      )) rho0_cup      = 0.
         if(allocated(t0_cup        )) t0_cup        = 0.
         if(allocated(thil0_cup     )) thil0_cup     = 0.
         if(allocated(theiv0_cup    )) theiv0_cup    = 0.
         if(allocated(theivs0_cup   )) theivs0_cup   = 0.

         if(allocated(co2_cup       )) co2_cup       = 0.
         if(allocated(exner_cup     )) exner_cup     = 0.
         if(allocated(p_cup         )) p_cup         = 0.
         if(allocated(qtot_cup      )) qtot_cup      = 0.
         if(allocated(qvap_cup      )) qvap_cup      = 0.
         if(allocated(qliq_cup      )) qliq_cup      = 0.
         if(allocated(qice_cup      )) qice_cup      = 0.
         if(allocated(qsat_cup      )) qsat_cup      = 0.
         if(allocated(rho_cup       )) rho_cup       = 0.
         if(allocated(t_cup         )) t_cup         = 0.
         if(allocated(thil_cup      )) thil_cup      = 0.
         if(allocated(theiv_cup     )) theiv_cup     = 0.
         if(allocated(theivs_cup    )) theivs_cup    = 0.

         if(allocated(co20d_cld     )) co20d_cld     = 0.
         if(allocated(dby0d         )) dby0d         = 0.
         if(allocated(pw0d_cld      )) pw0d_cld      = 0.
         if(allocated(qtot0d_cld    )) qtot0d_cld    = 0.
         if(allocated(qvap0d_cld    )) qvap0d_cld    = 0.
         if(allocated(qliq0d_cld    )) qliq0d_cld    = 0.
         if(allocated(qice0d_cld    )) qice0d_cld    = 0.
         if(allocated(qsat0d_cld    )) qsat0d_cld    = 0.
         if(allocated(rho0d_cld     )) rho0d_cld     = 0.
         if(allocated(t0d_cld       )) t0d_cld       = 0.
         if(allocated(theiv0d_cld   )) theiv0d_cld   = 0.
         if(allocated(thil0d_cld    )) thil0d_cld    = 0.

         if(allocated(cdd           )) cdd           = 0.
         if(allocated(co2d_cld      )) co2d_cld      = 0.
         if(allocated(dbyd          )) dbyd          = 0.
         if(allocated(etad_cld      )) etad_cld      = 0.
         if(allocated(mentrd_rate   )) mentrd_rate   = 0.
         if(allocated(pwd_cld       )) pwd_cld       = 0.
         if(allocated(qtotd_cld     )) qtotd_cld     = 0.
         if(allocated(qvapd_cld     )) qvapd_cld     = 0.
         if(allocated(qliqd_cld     )) qliqd_cld     = 0.
         if(allocated(qiced_cld     )) qiced_cld     = 0.
         if(allocated(qsatd_cld     )) qsatd_cld     = 0.
         if(allocated(rhod_cld      )) rhod_cld      = 0.
         if(allocated(td_cld        )) td_cld        = 0.
         if(allocated(theivd_cld    )) theivd_cld    = 0.
         if(allocated(thild_cld     )) thild_cld     = 0.

         if(allocated(co20u_cld     )) co20u_cld     = 0.
         if(allocated(dby0u         )) dby0u         = 0.
         if(allocated(pw0u_cld      )) pw0u_cld      = 0.
         if(allocated(qtot0u_cld    )) qtot0u_cld    = 0.
         if(allocated(qvap0u_cld    )) qvap0u_cld    = 0.
         if(allocated(qliq0u_cld    )) qliq0u_cld    = 0.
         if(allocated(qice0u_cld    )) qice0u_cld    = 0.
         if(allocated(qsat0u_cld    )) qsat0u_cld    = 0.
         if(allocated(rho0u_cld     )) rho0u_cld     = 0.
         if(allocated(t0u_cld       )) t0u_cld       = 0.
         if(allocated(thil0u_cld    )) thil0u_cld    = 0.
         if(allocated(theiv0u_cld   )) theiv0u_cld   = 0.

         if(allocated(cdu           )) cdu           = 0.
         if(allocated(co2u_cld      )) co2u_cld      = 0.
         if(allocated(dbyu          )) dbyu          = 0.
         if(allocated(etau_cld      )) etau_cld      = 0.
         if(allocated(mentru_rate   )) mentru_rate   = 0.
         if(allocated(pwu_cld       )) pwu_cld       = 0.
         if(allocated(qtotu_cld     )) qtotu_cld     = 0.
         if(allocated(qvapu_cld     )) qvapu_cld     = 0.
         if(allocated(qliqu_cld     )) qliqu_cld     = 0.
         if(allocated(qiceu_cld     )) qiceu_cld     = 0.
         if(allocated(qsatu_cld     )) qsatu_cld     = 0.
         if(allocated(rhou_cld      )) rhou_cld      = 0.
         if(allocated(tu_cld        )) tu_cld        = 0.
         if(allocated(thilu_cld     )) thilu_cld     = 0.
         if(allocated(theivu_cld    )) theivu_cld    = 0.

         !---------------------------------------------------------------------------------!
         ! Flushing scalars, we don't need to check for allocation here...                 !
         !---------------------------------------------------------------------------------!
         !----- Integer variables ---------------------------------------------------------!
         comp_dn           = .true.
         ierr              = 0
         klod              = 0
         klou              = 0
         klcl              = 0
         klfc              = 0
         kdet              = 0
         kstabi            = 0
         kstabm            = 0
         klnb              = 0
         ktop              = 0
         !----- Real variables. -----------------------------------------------------------!
         aa0u              = 0.
         aa0d              = 0.
         pwav0             = 0.
         pwav              = 0.
         pwev0             = 0.
         pwev              = 0.
         wbuoymin0         = 0.
         wbuoymin          = 0.
         !---------------------------------------------------------------------------------!

      case (2) 
         !---------------------------------------------------------------------------------!
         !     Here we reset only the variables that are modified by an arbitrary cloud    !
         ! flux of a given type of cloud.                                                  !
         !---------------------------------------------------------------------------------!
         if(allocated(x_co2         )) x_co2         = 0.
         if(allocated(x_qtot        )) x_qtot        = 0.
         if(allocated(x_qvap        )) x_qvap        = 0.
         if(allocated(x_qliq        )) x_qliq        = 0.
         if(allocated(x_qice        )) x_qice        = 0.
         if(allocated(x_t           )) x_t           = 0.
         if(allocated(x_theiv       )) x_theiv       = 0.
         if(allocated(x_thil        )) x_thil        = 0.

         if(allocated(x_co2_cup     )) x_co2_cup     = 0.
         if(allocated(x_exner_cup   )) x_exner_cup   = 0.
         if(allocated(x_p_cup       )) x_p_cup       = 0.
         if(allocated(x_qtot_cup    )) x_qtot_cup    = 0.
         if(allocated(x_qvap_cup    )) x_qvap_cup    = 0.
         if(allocated(x_qliq_cup    )) x_qliq_cup    = 0.
         if(allocated(x_qice_cup    )) x_qice_cup    = 0.
         if(allocated(x_qsat_cup    )) x_qsat_cup    = 0.
         if(allocated(x_rho_cup     )) x_rho_cup     = 0.
         if(allocated(x_t_cup       )) x_t_cup       = 0.
         if(allocated(x_thil_cup    )) x_thil_cup    = 0.
         if(allocated(x_theiv_cup   )) x_theiv_cup   = 0.
         if(allocated(x_theivs_cup  )) x_theivs_cup  = 0.

         if(allocated(x_co2d_cld    )) x_co2d_cld    = 0.
         if(allocated(x_dbyd        )) x_dbyd        = 0.
         if(allocated(x_pwd_cld     )) x_pwd_cld     = 0.
         if(allocated(x_qtotd_cld   )) x_qtotd_cld   = 0.
         if(allocated(x_qvapd_cld   )) x_qvapd_cld   = 0.
         if(allocated(x_qliqd_cld   )) x_qliqd_cld   = 0.
         if(allocated(x_qiced_cld   )) x_qiced_cld   = 0.
         if(allocated(x_qsatd_cld   )) x_qsatd_cld   = 0.
         if(allocated(x_rhod_cld    )) x_rhod_cld    = 0.
         if(allocated(x_td_cld      )) x_td_cld      = 0.
         if(allocated(x_theivd_cld  )) x_theivd_cld  = 0.
         if(allocated(x_thild_cld   )) x_thild_cld   = 0.

         if(allocated(x_co2u_cld    )) x_co2u_cld    = 0.
         if(allocated(x_dbyu        )) x_dbyu        = 0.
         if(allocated(x_pwu_cld     )) x_pwu_cld     = 0.
         if(allocated(x_qtotu_cld   )) x_qtotu_cld   = 0.
         if(allocated(x_qvapu_cld   )) x_qvapu_cld   = 0.
         if(allocated(x_qliqu_cld   )) x_qliqu_cld   = 0.
         if(allocated(x_qiceu_cld   )) x_qiceu_cld   = 0.
         if(allocated(x_qsatu_cld   )) x_qsatu_cld   = 0.
         if(allocated(x_rhou_cld    )) x_rhou_cld    = 0.
         if(allocated(x_tu_cld      )) x_tu_cld      = 0.
         if(allocated(x_thilu_cld   )) x_thilu_cld   = 0.
         if(allocated(x_theivu_cld  )) x_theivu_cld  = 0.
         !---------------------------------------------------------------------------------!
         ! Flushing scalars, we don't need to check for allocation here...                 !
         !---------------------------------------------------------------------------------!
         !----- Real variables. -----------------------------------------------------------!
         x_aad             = 0.
         x_aau             = 0.
         x_pwav            = 0.
         x_pwev            = 0.
         !---------------------------------------------------------------------------------!
      
      case (3) 
         !---------------------------------------------------------------------------------!
         !     Everything is zeroed, call the subroutine itself three times.               !
         !---------------------------------------------------------------------------------!
         call zero_scratch_grell(0)
         call zero_scratch_grell(1)
         call zero_scratch_grell(2)

      end select

      return

   end subroutine zero_scratch_grell
   !=======================================================================================!
   !=======================================================================================!
end module mem_scratch_grell
!==========================================================================================!
!==========================================================================================!

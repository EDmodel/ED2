!==========================================================================================!
! grell_cupar_main.f90                                                                     !
!                                                                                          !
!    This file contains the main Grell's subroutine. This subroutine is the same for both  !
! shallow and deep convection, the only difference is that shallow convection bypasses     !
! everything related to downdrafts.                                                        !
!    Whenever applicable, the follow convention is used:                                   !
!       [+] "?_cup"  : These are variables defined at the Grell level (staggered).         !
!       [+] "?d_cld" : These are cloud variables associated with downdrafts.               !
!       [+] "?u_cld" : These are cloud variables associated with updrafts.                 !
!       [+] "x_?"    : Variables modified by a fudging factor, for ensemble statistics     !
!       [+] "?0?"    : Whenever a variable like the ones above contains a 0, it means that !
!                      these variables are based on current values. If this is a BRAMS     !
!                      variable, it means the value before applying the tendency,          !
!                      otherwise it is a derived value based on the current large scale    !
!                      These are used only if Grell (1993) or Kain-Fritsch (1990)          !
!                      parameterizations are used.                                         !
!                                                                                          !
!------------------------------------------------------------------------------------------!
!==========================================================================================!
!==========================================================================================!
subroutine grell_cupar_main(closure_type,comp_down,comp_noforc_cldwork,comp_modif_thermo   &
                           ,checkmass,iupmethod,maxens_cap,maxens_dyn,maxens_eff           &
                           ,maxens_lsf,mgmzp,cap_maxs,cap_max_increment,wnorm_max          &
                           ,wnorm_increment,dtime,depth_min,depth_max,edtmax,edtmin        &
                           ,inv_ensdim,masstol,max_heat,pmass_left,radius,relheight_down   &
                           ,zkbmax,zcutdown,z_detr,edt_eff,dellatheiv_eff,dellathil_eff    &
                           ,dellaqtot_eff,pw_eff,dnmf_ens,upmf_ens,aad,aau,edt,dnmf,upmf   &
                           ,mynum)
               
   use mem_scratch_grell, only : &
      !----- Variables defined at the initialization process ------------------------------!
       dens_curr    & ! intent(in)  - Current density to strenghten updraft       [kg/m²/s]
      ,dzd_cld      & ! intent(in)  - Delta-height for downdraft calculations     [      m]
      ,dzu_cld      & ! intent(in)  - Delta-height for updraft calculations       [      m]
      ,exner0       & ! intent(in)  - Current Exner function                      [ J/kg/K]
      ,exner        & ! intent(in)  - Forced Exner funtion                        [ J/kg/K]
      ,exnersur     & ! intent(in)  - Surface Exner function                      [ J/kg/K]
      ,kpbl         & ! intent(in)  - Level of PBL top (Nakanishi/Niino only)     [    ---]
      ,mkx          & ! intent(in)  - Number of vertical levels                   [    ---]
      ,mconv        & ! intent(in)  - Integrated moisture convergence             [kg/m²/s]
      ,omeg         & ! intent(in)  - Vertical velocity in pressure, dp/dt:       [   Pa/s]
      ,p0           & ! intent(in)  - Current pressure                            [     Pa]
      ,p            & ! intent(in)  - Pressure with forcing                       [     Pa]
      ,psur         & ! intent(in)  - Surface pressure                            [     Pa]
      ,prev_dnmf    & ! intent(in)  - Previous call downdraft mass flux           [kg/m²/s]
      ,qtot0        & ! intent(in)  - Current total mixing ratio                  [  kg/kg]
      ,qtot         & ! intent(in)  - Total mixing ratio with forcing.            [  kg/kg]
      ,qtotsur      & ! intent(in)  - surface total mixing ratio                  [  kg/kg]
      ,qvap0        & ! intent(in)  - Current water vapour mixing ratio           [  kg/kg]
      ,qvap         & ! intent(in)  - Water vapour mixing ratio with forcing.     [  kg/kg]
      ,qvapsur      & ! intent(in)  - surface water vapour mixing ratio           [  kg/kg]
      ,qliq0        & ! intent(in)  - Current liquid water mixing ratio           [  kg/kg]
      ,qliq         & ! intent(in)  - Liquid water mixing ratio with forcing.     [  kg/kg]
      ,qliqsur      & ! intent(in)  - surface liquid water mixing ratio           [  kg/kg]
      ,qice0        & ! intent(in)  - Current ice mixing ratio                    [  kg/kg]
      ,qice         & ! intent(in)  - Ice mixing ratio with forcing.              [  kg/kg]
      ,qicesur      & ! intent(in)  - surface ice mixing ratio                    [  kg/kg]
      ,rho          & ! intent(in)  - Air density                                 [  kg/m³]
      ,theiv0       & ! intent(in)  - Current ice-vapour equiv. pot. temp.        [      K]
      ,theiv        & ! intent(in)  - Ice-vapour equiv. pot. temp. with forcing   [      K]
      ,theivsur     & ! intent(in)  - Surface ice-vapour equiv. pot. temp.        [      K]
      ,thil0        & ! intent(in)  - Current ice-liquid potential temp.          [      K]
      ,thil         & ! intent(in)  - Ice-liquid potential temp. with forcing     [      K]
      ,thilsur      & ! intent(in)  - Surface ice-liquid potential temp.          [      K]
      ,t0           & ! intent(in)  - Current ice-liquid potential temp.          [      K]
      ,t            & ! intent(in)  - Temperature with forcing.                   [      K]
      ,tsur         & ! intent(in)  - Surface temperature                         [      K]
      ,tke0         & ! intent(in)  - Turbulent Kinetic Energy                    [   J/kg]
      ,tke          & ! intent(in)  - Turbulent Kinetic Energy with forcing       [   J/kg]
      ,tscal_kf     & ! intent(in)  - Time scale for Kain-Fritsch (1990)          [      s]
      ,sigw         & ! intent(in)  - Vertical velocity standard deviation        [    m/s]
      ,upconv       & ! intent(in)  - Flag: is convection happening upstream?     [    ---]
      ,uwind        & ! intent(in)  - Zonal wind                                  [    m/s]
      ,vwind        & ! intent(in)  - Meridional wind                             [    m/s]
      ,wwind        & ! intent(in)  - Mean vertical velocity                      [    m/s]
      ,z            & ! intent(in)  - Height                                      [      m]
      ,z_cup        & ! intent(in)  - Height at cloud levels                      [      m]
      !------ Variables defined here ------------------------------------------------------!
      ,cdd          & ! intent(out) - Normalized downdraft detr. function         [    1/m]
      ,cdu          & ! intent(out) - Normalized updraft detr. function           [    1/m]
      ,dbyd         & ! intent(out) - Downdraft buoyancy acceleration             [   m/s²]
      ,dbyu         & ! intent(out) - Updraft buoyancy acceleration               [   m/s²]
      ,etad_cld     & ! intent(out) - Normalized downdraft mass flux              [   ----]
      ,etau_cld     & ! intent(out) - Normalized updraft Mass flux                [   ----]
      ,ierr         & ! intent(out) - Flag for convection error                   [   ----]
      ,jmin         & ! intent(out) - Downdraft originating level                 [   ----]
      ,k22          & ! intent(out) - Updraft originating level                   [   ----]
      ,kbcon        & ! intent(out) - Level of free convection                    [   ----]
      ,kdet         & ! intent(out) - Top of downdraft detrainment layer          [   ----]
      ,kstabi       & ! intent(out) - Stable layer bottom                         [   ----]
      ,kstabm       & ! intent(out) - Stable layer top                            [   ----]
      ,ktop         & ! intent(out) - Cloud top                                   [   ----]
      ,mentrd_rate  & ! intent(out) - Normalized downdraft entrainment rate       [    1/m]
      ,mentru_rate  & ! intent(out) - Normalized updraft entrainment rate         [    1/m]
      ,outthil      & ! intent(out) - Ice-liquid pot. temp. tendency              [    K/s]
      ,outqtot      & ! intent(out) - Total H20 mixing ratio tendency             [kg/kg/s]
      ,precip       & ! intent(out) - Precipitation rate                          [kg/m²/s]
      ,rhod_cld     & ! intent(out) - Downdraft density                           [  kg/m³]
      ,rhou_cld     & ! intent(out) - Updraft density                             [  kg/m³]
      ,qliqd_cld    & ! intent(out) - Downdraft liquid water mixing ratio         [  kg/kg]
      ,qliqu_cld    & ! intent(out) - Updraft liquid water mixing ratio           [  kg/kg]
      ,qiced_cld    & ! intent(out) - Downdraft ice mixing ratio                  [  kg/kg]
      ,qiceu_cld    & ! intent(out) - Updraft ice mixing ratio                    [  kg/kg]
      ,wbuoymin     ! ! intent(out) - Minimum buoyant vertical velocity           [    m/s]

   use rconstants, only: cp,cpi,alvl,alvi,ep,epi,rgas,g,ttripoli,t3ple,day_sec,toodry

   implicit none
   
   !---------------------------------------------------------------------------------------!
   ! List of arguments                                                                     !
   !---------------------------------------------------------------------------------------!
   !----- Character, containing the closure type(s) for dynamic control -------------------!
   character(len=2), intent(in) :: closure_type! Short name to define the method.
   !----- Logical flags, to bypass uncessary steps ----------------------------------------!
   logical , intent(in) :: comp_down           ! I will compute downdrafts.           [T/F]
   logical , intent(in) :: comp_noforc_cldwork ! I will compute no forced cloud work  [T/F]
   logical , intent(in) :: comp_modif_thermo   ! I will compute the x_* variables     [T/F]
   logical , intent(in) :: checkmass           ! I will check mass balance            [T/F]
   !----- Integer, global dimensions ------------------------------------------------------!
   integer , intent(in) :: iupmethod        ! Method to find the updraft originatin level
   integer , intent(in) :: maxens_cap       ! Ensemble size on static control (cap_maxs)
   integer , intent(in) :: maxens_dyn       ! Ensemble size on dynamic control 
   integer , intent(in) :: maxens_eff       ! Ensemble size on precipitation efficiency
   integer , intent(in) :: maxens_lsf       ! Ensemble size on large scale perturbations
   integer , intent(in) :: mgmzp            ! Vertical grid size
   integer , intent(in) :: mynum            ! Node ID, for debugging purposes only

   !----- Real, parameters to define the cloud --------------------------------------------!
   real    , intent(in) :: cap_maxs         ! Maximum depth of capping inversion     [ hPa]
   real    , intent(in) :: cap_max_increment! Extra cap_maxs due to upstream conv.   [ hPa]
   real    , intent(in) :: dtime            ! Time step                              [   s]
   real    , intent(in) :: depth_min        ! Minimum cloud depth to qualify it      [   m]
   real    , intent(in) :: depth_max        ! Maximum cloud depth to qualify it      [   m]
   real    , intent(in) :: edtmax           ! Maximum epsilon (dnmf/upmf)            [ ---]
   real    , intent(in) :: edtmin           ! Minimum epsilon (dnmf/upmf)            [ ---]
   real    , intent(in) :: inv_ensdim       ! Inverse of ensemble dimension size     [ ---]
   real    , intent(in) :: masstol          ! Maximum mass leak allowed to happen    [ ---]
   real    , intent(in) :: max_heat         ! Maximum heating scale                  [K/dy]
   real    , intent(in) :: pmass_left       ! Fraction of mass left at the ground    [ ---]
   real    , intent(in) :: radius           ! Radius, for entrainment rate.          [   m]
   real    , intent(in) :: relheight_down   ! Relative height for downdraft origin   [ ---]
   real    , intent(in) :: wnorm_max        ! Maximum normalised w to be considered  [ ---]
   real    , intent(in) :: wnorm_increment  ! Extra wnorm for ensemble               [ ---]
   real    , intent(in) :: zkbmax           ! Top height for updrafts to originate   [   m]
   real    , intent(in) :: zcutdown         ! Top height for downdrafts to originate [   m]
   real    , intent(in) :: z_detr           ! Top height for downdraft detrainment   [   m]

   !----- Ensemble variable, fraction between reference downdraft and updraft -------------!
   real, dimension(maxens_eff,maxens_cap)      , intent(out) :: & 
                           edt_eff          ! Grell's epsilon                     [    ---]
   !----- Ensemble variables, changes of the thermodynamic properties per unit of mass. ---!
   real, dimension(mgmzp,maxens_eff,maxens_cap), intent(inout) :: &
                           dellatheiv_eff & ! Ice-vapour equiv. pot. temperature  [      K]
                          ,dellathil_eff  & ! Ice-liquid potential temperature    [      K]
                          ,dellaqtot_eff  & ! Total mixing ratio                  [  kg/kg]
                          ,pw_eff         ! ! Fall-out water                      [  kg/kg]
   !----- Ensemble variables, reference mass fluxes, in [kg/m²/s] -------------------------!
   real, dimension(maxens_dyn,maxens_lsf,maxens_eff,maxens_cap), intent(out) :: & 
                           dnmf_ens       & ! Reference downdraft flux            [kg/m²/s]
                          ,upmf_ens       ! ! Reference updraft flux              [kg/m²/s]
   !----- Output variables, the reference mass fluxes. Others will be found afterwards. ---!
   real                  , intent(out)   :: aad  ! Downdraft work function        [   J/kg]
   real                  , intent(out)   :: aau  ! Updraft work function          [   J/kg]
   real                  , intent(out)   :: edt  ! dnmf/upmf                      [    ---]
   real                  , intent(out)   :: dnmf ! Reference downward mass flux   [kg/m²/s]
   real                  , intent(out)   :: upmf ! Reference upward mass flux     [kg/m²/s]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Local variables.                                                                   !
   !---------------------------------------------------------------------------------------!
   !----- Counters and level holders ------------------------------------------------------!
   integer                   :: icap          ! Static control (cap_max) counter;
   integer                   :: iedt          ! Precipitation efficiency counter;
   integer                   :: imbp          ! Response variable fudge counter;
   integer                   :: k             ! Generic vertical level counter;
   integer                   :: kbmin         ! Auxiliary variable;
   integer                   :: kbmax         ! Top level that updrafts can generate;
   integer                   :: kzdown        ! Actual maximum origin for downdrafts;
   integer                   :: capoffset     ! Offset for perturbation of cap_max
   integer                   :: mboffset      ! Offset for the perturbation of forcing;

   !---------------------------------------------------------------------------------------!
   !     Environment variables at Grell's staggered grid.                                  !
   !---------------------------------------------------------------------------------------!
   !----- Variables based on values before any forcing is applied -------------------------!
   real   , dimension(mgmzp) :: exner0_cup    ! Exner function                    [   J/kg]
   real   , dimension(mgmzp) :: p0_cup        ! Pressure                          [     Pa]
   real   , dimension(mgmzp) :: qtot0_cup     ! Total mixing ratio                [  kg/kg]
   real   , dimension(mgmzp) :: qvap0_cup     ! Water vapour mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: qliq0_cup     ! Liquid water                      [  kg/kg]
   real   , dimension(mgmzp) :: qice0_cup     ! Mixing ratio                      [  kg/kg]
   real   , dimension(mgmzp) :: qsat0_cup     ! Saturation mixing ratio           [  kg/kg]
   real   , dimension(mgmzp) :: rho0_cup      ! Density                           [  kg/m³]
   real   , dimension(mgmzp) :: t0_cup        ! Temperature                       [      K]
   real   , dimension(mgmzp) :: thil0_cup     ! Ice-liquid potential temperature  [      K]
   real   , dimension(mgmzp) :: theiv0_cup    ! Ice-vapour potential temperature  [      K]
   real   , dimension(mgmzp) :: theivs0_cup   ! Sat. Ice-vapour equiv. pot. temp. [      K]
   !----- Variables with all time tendencies but this convection. -------------------------!
   real   , dimension(mgmzp) :: exner_cup     ! Exner function                    [   J/kg]
   real   , dimension(mgmzp) :: p_cup         ! Pressure                          [     Pa]
   real   , dimension(mgmzp) :: qtot_cup      ! Total mixing ratio                [  kg/kg]
   real   , dimension(mgmzp) :: qvap_cup      ! Water vapour mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: qliq_cup      ! Liquid water                      [  kg/kg]
   real   , dimension(mgmzp) :: qice_cup      ! Mixing ratio                      [  kg/kg]
   real   , dimension(mgmzp) :: qsat_cup      ! Saturation mixing ratio           [  kg/kg]
   real   , dimension(mgmzp) :: rho_cup       ! Density                           [  kg/m³]
   real   , dimension(mgmzp) :: t_cup         ! Temperature                       [      K]
   real   , dimension(mgmzp) :: thil_cup      ! Ice-liquid potential temperature  [      K]
   real   , dimension(mgmzp) :: theiv_cup     ! Ice-vapour equiv. pot. temp.      [      K]
   real   , dimension(mgmzp) :: theivs_cup    ! Sat. Ice-vapour equiv. pot. temp. [      K]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Downdraft-related variables                                                        !
   !---------------------------------------------------------------------------------------!
   !----- Based on current values (no forcing) --------------------------------------------!
   real   , dimension(mgmzp) :: dby0d         ! Buoyancy acceleration             [   m/s²]
   real   , dimension(mgmzp) :: pw0d_cld      ! Condensation                      [  kg/kg]
   real   , dimension(mgmzp) :: qtot0d_cld    ! Total water mixing ratio          [  kg/kg]
   real   , dimension(mgmzp) :: qvap0d_cld    ! Vapour mixing ratio               [  kg/kg]
   real   , dimension(mgmzp) :: qliq0d_cld    ! Liquid water mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: qice0d_cld    ! Ice mixing ratio                  [  kg/kg]
   real   , dimension(mgmzp) :: qsat0d_cld    ! Sat. mixing ratio                 [  kg/kg]
   real   , dimension(mgmzp) :: rho0d_cld     ! Density                           [  kg/m³]
   real   , dimension(mgmzp) :: t0d_cld       ! Temperature                       [      K]
   real   , dimension(mgmzp) :: theiv0d_cld   ! Ice-vapour equiv. pot. temp.      [      K]
   real   , dimension(mgmzp) :: thil0d_cld    ! Ice-liquid pot. temperature       [      K]
   real                      :: aa0d          ! Updraft work function             [   J/kg]
   real                      :: pwev0         ! Integrated evaporation            [  kg/kg]
   !----- Based on values with forcing ----------------------------------------------------!
   real   , dimension(mgmzp) :: pwd_cld       ! Condensation                      [  kg/kg]
   real   , dimension(mgmzp) :: qtotd_cld     ! Total water mixing ratio          [  kg/kg]
   real   , dimension(mgmzp) :: qvapd_cld     ! Vapour mixing ratio               [  kg/kg]
   real   , dimension(mgmzp) :: qsatd_cld     ! Sat. mixing ratio                 [  kg/kg]
   real   , dimension(mgmzp) :: td_cld        ! Temperature                       [      K]
   real   , dimension(mgmzp) :: theivd_cld    ! Ice-vapour equiv. pot. temp.      [      K]
   real   , dimension(mgmzp) :: thild_cld     ! Ice-liquid pot. temperature       [      K]
   real                      :: pwev          ! Integrated evaporation            [  kg/kg]
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Updraft-related variables                                                          !
   !---------------------------------------------------------------------------------------!
   !----- Based on current values (no forcing) --------------------------------------------!
   real   , dimension(mgmzp) :: dby0u         ! Buoyancy acceleration             [   m/s²]
   real   , dimension(mgmzp) :: pw0u_cld      ! Condensation                      [  kg/kg]
   real   , dimension(mgmzp) :: qtot0u_cld    ! Total water mixing ratio          [  kg/kg]
   real   , dimension(mgmzp) :: qvap0u_cld    ! Vapour mixing ratio               [  kg/kg]
   real   , dimension(mgmzp) :: qliq0u_cld    ! Liquid water mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: qice0u_cld    ! Ice mixing ratio                  [  kg/kg]
   real   , dimension(mgmzp) :: qsat0u_cld    ! Sat. mixing ratio                 [  kg/kg]
   real   , dimension(mgmzp) :: rho0u_cld     ! Density                           [  kg/m³]
   real   , dimension(mgmzp) :: t0u_cld       ! Temperature                       [      K]
   real   , dimension(mgmzp) :: thil0u_cld    ! Ice-liquid potential temperature  [      K]
   real   , dimension(mgmzp) :: theiv0u_cld   ! Ice-vapour equiv. pot. temp.      [      K]
   real                      :: aa0u          ! Updraft work function             [   J/kg]
   real                      :: pwav0         ! Integrated condensation           [  kg/kg]
   real                      :: wbuoymin0     ! Minimum buoyant velocity          [    m/s]
   !----- Based on values with forcing ----------------------------------------------------!
   real   , dimension(mgmzp) :: pwu_cld       ! Condensation                      [  kg/kg]
   real   , dimension(mgmzp) :: qtotu_cld     ! Total water mixing ratio          [  kg/kg]
   real   , dimension(mgmzp) :: qvapu_cld     ! Vapour mixing ratio               [  kg/kg]
   real   , dimension(mgmzp) :: qsatu_cld     ! Sat. mixing ratio                 [  kg/kg]
   real   , dimension(mgmzp) :: tu_cld        ! Temperature                       [      K]
   real   , dimension(mgmzp) :: thilu_cld     ! Ice-liquid potential temperature  [      K]
   real   , dimension(mgmzp) :: theivu_cld    ! Ice-vapour equiv. pot. temp.      [      K]
   real                      :: pwav          ! Integrated condensation           [  kg/kg]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Combined updraft/downdraft                                                         !
   !---------------------------------------------------------------------------------------!
   real                      :: aatot0        ! Current total cloud work function [   J/kg]
   real                      :: aatot         ! Total cloud work func. (forced)   [   J/kg]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Flag for convection error when using the "0"  variables.                           !
   !---------------------------------------------------------------------------------------!
   integer :: ierr0     ! Flag for convection error
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Variables modified by the arbitrary mass flux                                      !
   !---------------------------------------------------------------------------------------!
   !----- Model levels --------------------------------------------------------------------!
   real   , dimension(mgmzp) :: x_qtot        ! Total mixing ratio                [  kg/kg]
   real   , dimension(mgmzp) :: x_qvap        ! Vapour mixing ratio               [  kg/kg]
   real   , dimension(mgmzp) :: x_qliq        ! Liquid water mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: x_qice        ! Ice mixing ratio                  [  kg/kg]
   real   , dimension(mgmzp) :: x_t           ! Temperature                       [      K]
   real   , dimension(mgmzp) :: x_theiv       ! Ice-vapour equiv. pot. temp.      [      K]
   real   , dimension(mgmzp) :: x_thil        ! Ice-liquid potential temperature. [      K]
   !----- Cloud levels --------------------------------------------------------------------!
   real   , dimension(mgmzp) :: x_exner_cup   ! Exner function                    [   J/kg]
   real   , dimension(mgmzp) :: x_p_cup       ! Pressure                          [     Pa]
   real   , dimension(mgmzp) :: x_qtot_cup    ! Total mixing ratio                [  kg/kg]
   real   , dimension(mgmzp) :: x_qvap_cup    ! Water vapour mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: x_qliq_cup    ! Liquid water                      [  kg/kg]
   real   , dimension(mgmzp) :: x_qice_cup    ! Mixing ratio                      [  kg/kg]
   real   , dimension(mgmzp) :: x_qsat_cup    ! Saturation mixing ratio           [  kg/kg]
   real   , dimension(mgmzp) :: x_rho_cup     ! Density                           [  kg/m³]
   real   , dimension(mgmzp) :: x_t_cup       ! Temperature                       [      K]
   real   , dimension(mgmzp) :: x_thil_cup    ! Ice-liquid potential temperature  [      K]
   real   , dimension(mgmzp) :: x_theiv_cup   ! Ice-vapour equiv. pot. temp.      [      K]
   real   , dimension(mgmzp) :: x_theivs_cup  ! Sat. Ice-vapour equiv. pot. temp. [      K]
   !----- Downdraft variables -------------------------------------------------------------!
   real   , dimension(mgmzp) :: x_dbyd        ! Buoyancy acceleration             [   m/s²]
   real   , dimension(mgmzp) :: x_pwd_cld     ! Condensation                      [  kg/kg]
   real   , dimension(mgmzp) :: x_qtotd_cld   ! Total water mixing ratio          [  kg/kg]
   real   , dimension(mgmzp) :: x_qvapd_cld   ! Vapour mixing ratio               [  kg/kg]
   real   , dimension(mgmzp) :: x_qliqd_cld   ! Liquid water mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: x_qiced_cld   ! Ice mixing ratio                  [  kg/kg]
   real   , dimension(mgmzp) :: x_qsatd_cld   ! Sat. mixing ratio                 [  kg/kg]
   real   , dimension(mgmzp) :: x_rhod_cld    ! Density                           [  kg/m³]
   real   , dimension(mgmzp) :: x_td_cld      ! Temperature                       [      K]
   real   , dimension(mgmzp) :: x_theivd_cld  ! Ice-vapour equiv. pot. temp.      [      K]
   real   , dimension(mgmzp) :: x_thild_cld   ! Ice-liquid potential temperature  [      K]
   real                      :: x_aad         ! Updraft work function             [   J/kg]
   real                      :: x_pwev        ! Integrated evaporation            [  kg/kg]
   real                      :: x_wbuoymin    ! Minimum buoyant velocity          [    m/s]
   !----- Updraft variables ---------------------------------------------------------------!
   real   , dimension(mgmzp) :: x_dbyu        ! Buoyancy acceleration             [   m/s²]
   real   , dimension(mgmzp) :: x_pwu_cld     ! Condensation                      [  kg/kg]
   real   , dimension(mgmzp) :: x_qtotu_cld   ! Total water mixing ratio          [  kg/kg]
   real   , dimension(mgmzp) :: x_qvapu_cld   ! Vapour mixing ratio               [  kg/kg]
   real   , dimension(mgmzp) :: x_qliqu_cld   ! Liquid water mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: x_qiceu_cld   ! Ice mixing ratio                  [  kg/kg]
   real   , dimension(mgmzp) :: x_qsatu_cld   ! Sat. mixing ratio                 [  kg/kg]
   real   , dimension(mgmzp) :: x_rhou_cld    ! Density                           [  kg/m³]
   real   , dimension(mgmzp) :: x_tu_cld      ! Temperature                       [      K]
   real   , dimension(mgmzp) :: x_thilu_cld   ! Ice-liquid potential temperature  [      K]
   real   , dimension(mgmzp) :: x_theivu_cld  ! Ice-vapour equiv. pot. temp.      [      K]
   real                      :: x_aau         ! Updraft work function             [   J/kg]
   real                      :: x_pwav        ! Integrated condensation           [  kg/kg]
   !----- Combined variables --------------------------------------------------------------!
   real                      :: x_aatot       ! Combined cloud work function      [   J/kg]
   !----- Other parameters ----------------------------------------------------------------!
   integer                   :: x_ierr        ! Flag for convection error
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Miscellaneous parameters                                                           !
   !---------------------------------------------------------------------------------------!
   real :: cap_max    ! Actual "Depth" of inversion capping                       [     Pa]
   real :: mbprime    ! Current arbitrary mass flux fudging factor for ensemble   [kg/m²/s]
   real :: one_b      ! (1-b), as defined in Krishnamurti et al. (1983)           [    ---]
   real :: wbuoy_max  ! Maximum acceptable buoyant velocity                       [    m/s]
   real :: zktop      ! Actual maximum height that downdrafts can originate.      [      m]
   real :: pmass_kept ! 1.-pmass_left                                             [    ---]
   !----- Scratch array -------------------------------------------------------------------!
   real, dimension(mgmzp) :: scrvar    ! Scratch variable
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !  Setting up the mass that is kept by downdraft                                        !
   !---------------------------------------------------------------------------------------!
   pmass_kept = 1.-pmass_left



   !---------------------------------------------------------------------------------------!
   !    Initialise reference downdraft and updraft reference.                              !
   !---------------------------------------------------------------------------------------!
   dnmf = 0.
   upmf = 0.


   !---------------------------------------------------------------------------------------!
   !    This is the beginning of the outermost loop, on static control. Currently it is    !
   ! set for cap_maxs sensitivity but any other parameter could be used. Also in case you  !
   ! feel like including more terms, it is just a matter of including more loops.          !
   !                                                                                       !
   !    IMPORTANT: This loop is done from the end to the first element because for all     !
   !               ensemble tests, 1 is the reference, so to avoid storing larger          !
   !               dimension variables, I simply start from the end and when I am done     !
   !               the variables will have the values I need.                              !
   !---------------------------------------------------------------------------------------!
   stacloop: do icap=maxens_cap,1,-1
   
      capoffset= (-1)**icap * icap / 2 !----- This gives 0,1,-1,2,-2 for icap=1,2,3... ----!
      if (cap_maxs > 0.) then
         !---------------------------------------------------------------------------------!
         ! A1. Define cap_max for this member if the user opted by cap_maxs (i.e.,         !
         !     cap_max > 0). The reference given by the user will be used as               !
         !     reference, but the actual cap_max will be increasedd in case upstream       !
         !     convection has happened and the user asked to consider upstream convection. !
         !     Also, cap_max will be perturbed for members other than the first.           !
         !---------------------------------------------------------------------------------!
         if (upconv) then
            cap_max =  max(cap_max_increment,cap_maxs + (capoffset+1)*cap_max_increment)
         else
            cap_max =  max(cap_max_increment,cap_maxs + capoffset * cap_max_increment)
         end if
         wbuoy_max = 0.
      else
      !------------------------------------------------------------------------------------!
      ! A2. If the user opted by the probability, then I use capoffset to perturb the      !
      !     threshold in which I should accept convection. threwnorm is the maximum norm-  !
      !     alized value of vertical wind that will give a cumulative probability density  !
      !     function (CPDF) greater than what the user asked, assuming that the vertical   !
      !     velocity with average equals to wp and standard deviation equals to sigw (from !
      !     Mellor-Yamada closures).                                                       ! 
      !------------------------------------------------------------------------------------!
         cap_max = cap_maxs
         if (upconv) then
            wbuoy_max = wnorm_max + (capoffset+1) * wnorm_increment
         else
            wbuoy_max = wnorm_max + capoffset * wnorm_increment
         end if
      end if


      !------------------------------------------------------------------------------------!
      ! B. This flag states whether convection happened or failed. If ierr=0, then         !
      !    convection  happened, otherwise something kept it to develop. Whenever the      !
      !    parameterization finds a condition to impeach convection, a different number is !
      !    assigned and the static control is cycled. The error table is the following:    !
      !                                                                                    !
      !  1. Too early in the simulation for convection to be called.                       !
      !  2. The level updrafts origin is way too high;                                     !
      !  3. Inversion capping may be too thick so the LFC is probably out of reach;        !
      !  4. Couldn't find a layer with negative buoyancy for downdrafts;                   !
      !  5. Couldn't find a suitable cloud top, it would be above the domain top;          !
      !  6. This cloud would be too thin to fall in this spectral type;                    !
      !  7. This cloud would be too thick to fall in this spectral type;                   !
      !  8. Forced downdraft layer would have positive buoyancy, which doesn't make sense; !
      !  9. Downdraft would require more water than what is available to stay saturated;   !
      ! 10. Cloud work function associated with updraft is zero.                           !
      ! 11. Reference upward mass flux is zero.                                            !
      ! 12. Downdrafts would happen below the updrafts origin.                             !
      !------------------------------------------------------------------------------------!
      ierr   = 0

      !------------------------------------------------------------------------------------!
      ! C. Initialize Entrainment and Detrainment variables.                               !
      !------------------------------------------------------------------------------------!
      mentrd_rate(1:mkx) = .2/radius          ! Entrainment rate associated with dndrafts
      mentru_rate(1:mkx) = mentrd_rate(1:mkx) ! Entrainment rate associated with updrafts
      cdd(1:mkx) = 0.                         ! Detrainment function associated w/ dndrafts
      cdu(1:mkx) = 0.1*mentru_rate(1:mkx)     ! Detrainment function associated w/ updrafts

      !------------------------------------------------------------------------------------!
      ! D. Initialize mass fluxes, since they may never be computed in case convection     !
      !    fails.                                                                          !
      !------------------------------------------------------------------------------------!
      etad_cld(1:mkx) = 0.
      etau_cld(1:mkx) = 0.

      !------------------------------------------------------------------------------------!
      ! E. Initialize Cloud work functions associated with updrafts and downdrafts.        !
      !------------------------------------------------------------------------------------!
      aau    = 0.
      aad    = 0.
      aatot  = 0.
      edt    = 0.
   
      !------------------------------------------------------------------------------------!
      ! F. Initialize all level-related cloud variables. They are all output variables and !
      !    some may never be assigned when convection fails. Except for kstabm, I am set-  !
      !    ting all to a non-sense value to make the point that convection did not happen. !
      !------------------------------------------------------------------------------------!
      jmin   = 0
      k22    = 0
      kbcon  = 0
      kdet   = 0
      ktop   = 0
      kstabi = 0
      kstabm = mkx -2

      !------------------------------------------------------------------------------------!
      ! G. Calculate all thermodynamic properties at the cloud level. The cloud levels are !
      !    staggered in relation to BRAMS model.                                           !
      !------------------------------------------------------------------------------------!
      call grell_thermo_cldlev(mkx,mgmzp,z_cup,exner,thil,t,qtot,qliq,qice,exnersur        &
                              ,thilsur,tsur,qtotsur,qliqsur,qicesur,exner_cup,p_cup,t_cup  &
                              ,thil_cup,qtot_cup,qvap_cup,qliq_cup,qice_cup,qsat_cup       &
                              ,rho_cup,theiv_cup,theivs_cup)

      !------------------------------------------------------------------------------------!
      ! H. Initialize updraft temperature and drafts liquid mixing ratio, in case          !
      !    convection fails. These variables are used in other modules, such as CATT and   !
      !    Harrington.                                                                     !
      !------------------------------------------------------------------------------------!
      qliqd_cld(1:mkx) = 0.
      qliqu_cld(1:mkx) = 0.
      qiced_cld(1:mkx) = 0.
      qiceu_cld(1:mkx) = 0.



      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! I. Finding some updraft properties, namely the levels in which updrafts originate, !
      !    the level of free convection and the cloud top, and the updraft-related energy, !
      !    mass, and moisture properties.                                                  !
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! 1. Find the topmost level in which updrafts can originate based on this grid.      !
      !------------------------------------------------------------------------------------!
      kbmaxloop: do kbmax=1,mkx
         if (z_cup(kbmax) > zkbmax) exit kbmaxloop
      end do kbmaxloop

      !------------------------------------------------------------------------------------!
      ! 2. Determine a first guess for k22, the level of origin of updrafts. Check if that !
      !    didn't prevent convection to happen.                                            !
      !------------------------------------------------------------------------------------!
      call grell_updraft_origin(mkx,mgmzp,iupmethod,kpbl,kbmax,z,tke,qice,qliq,theiv_cup   &
                               ,ierr,k22)


      if (ierr /= 0) cycle stacloop
      !------------------------------------------------------------------------------------!
      ! 3. Finding the level of free convection. Two important points here:                !
      !    a. This call may end up preventing convection, so I must check after the call   !
      !    b. This subroutine may also affect the updraft originating level.               !
      !------------------------------------------------------------------------------------!
      call grell_find_cloud_lfc(mkx,mgmzp,kbmax,cap_max,wbuoy_max,wwind,sigw,exner_cup     &
                               ,p_cup,theiv_cup,thil_cup,t_cup,qtot_cup,qvap_cup,qliq_cup  &
                               ,qice_cup,qsat_cup,rho_cup,dzd_cld,mentru_rate,theivu_cld   &
                               ,thilu_cld,tu_cld,qtotu_cld,qvapu_cld,qliqu_cld,qiceu_cld   &
                               ,qsatu_cld,rhou_cld,dbyu,k22,ierr,kbcon,wbuoymin)
      if (ierr /= 0) cycle stacloop

      !------------------------------------------------------------------------------------!
      ! 4. Finding the minimum saturation thetae_iv. This will be the bottom of the stable !
      !    layer.                                                                          !
      !------------------------------------------------------------------------------------!
      kstabi=(kbcon - 1) + minloc(theivs_cup(kbcon:kstabm),dim=1)

      !------------------------------------------------------------------------------------!
      ! 5. Increasing the detrainment in stable layers provided that there is such layer.  !
      !    this rate increases linearly until a maximum value, currently set to 10 times   !
      !    the entrainment rate.                                                           !
      !------------------------------------------------------------------------------------!
      do k=kstabi,kstabm-1
         cdu(k) = min(10.*mentru_rate(k) , cdu(k-1) + 1.5 * mentru_rate(k))
      end do

      !------------------------------------------------------------------------------------!
      ! 6. Calculate the incloud ice-vapour equivalent potential temperature               !
      !------------------------------------------------------------------------------------!
      call grell_theiv_updraft(mkx,mgmzp,k22,kbcon,cdu,mentru_rate,theiv,theiv_cup,dzu_cld &
                              ,theivu_cld)
                           
      !------------------------------------------------------------------------------------!
      ! 7. Finding the cloud top. Since this may keep convection to happen, I check        !
      !    whether I should move on or break here.                                         !
      !------------------------------------------------------------------------------------!
      call grell_find_cloud_top(mkx,mgmzp,kbcon,cdu,mentru_rate,theiv_cup,theivs_cup       &
                               ,dzu_cld,theivu_cld,dbyu,ierr,ktop)
      if (ierr /= 0) cycle stacloop

      !------------------------------------------------------------------------------------!
      ! 8. Also check whether this cloud qualifies to be in this spectrum size. It need    !
      !    to be thicker than the minimum value provided by the user.                      !
      !------------------------------------------------------------------------------------!
      !if (z_cup(ktop)-z_cup(kbcon) < depth_min) then
      !   ierr = 6
      !   cycle stacloop
      !elseif (z_cup(ktop)-z_cup(kbcon) > depth_max) then
      !   ierr = 7
      !   cycle stacloop
      !end if
      if (z_cup(ktop)-z_cup(k22) < depth_min) then
         ierr = 6
         cycle stacloop
      elseif (z_cup(ktop)-z_cup(k22) > depth_max) then
         ierr = 7
         cycle stacloop
      end if

      !------------------------------------------------------------------------------------!
      ! 9. Finding the normalized mass fluxes associated with updrafts. Since we are using !
      !    height-based vertical coordinate, there is no need to find the forced           !
      !    normalized mass flux, they'll be the same, so just copy it afterwards.          !
      !------------------------------------------------------------------------------------!
      call grell_nms_updraft(mkx,mgmzp,k22,kbcon,ktop,mentru_rate,cdu,dzu_cld,etau_cld)

      !------------------------------------------------------------------------------------!
      !10. Finding the moisture profiles associated with updrafts.                         !
      !------------------------------------------------------------------------------------!
      call grell_most_thermo_updraft(comp_down,mkx,mgmzp,kbcon,ktop,cdu,mentru_rate,qtot   &
                                    ,p_cup,exner_cup,thil_cup,t_cup,qtot_cup,qvap_cup      &
                                    ,qliq_cup,qice_cup,qsat_cup,rho_cup,theivu_cld         &
                                    ,etau_cld,dzu_cld,thilu_cld,tu_cld,qtotu_cld,qvapu_cld &
                                    ,qliqu_cld,qiceu_cld,qsatu_cld,rhou_cld,dbyu,pwu_cld   &
                                    ,pwav)

      !------------------------------------------------------------------------------------!
      !11. Finding the cloud work function associated with updrafts. If this cloud doesn't !
      !    produce cloud work, break the run, we don't simulate lazy clouds in this model. !
      !------------------------------------------------------------------------------------!
      call grell_cldwork_updraft(mkx,mgmzp,kbcon,ktop,dbyu,dzu_cld,etau_cld,aau)
      if (aau == 0.) then
         ierr = 10
         cycle stacloop
      end if

      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!



      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! J. Finding the downdraft counterpart of the above properties, namely where the     !
      !    downdrafts detrain all its mass, where they originate, their mass, energy and   !
      !    moisture properties. This should be done only when it is a cloud that has       !
      !    downdrafts. !                                                                   !
      !------------------------------------------------------------------------------------!
      if (comp_down) then
         !---------------------------------------------------------------------------------!
         ! 1. The level in which downdrafts should detrain (almost) all its mass.          !
         !---------------------------------------------------------------------------------!
         kdetloop: do kdet=1,mkx
            if (z_cup(kdet) > z_detr) exit kdetloop
         end do kdetloop

         !---------------------------------------------------------------------------------!
         ! 2. The level in which downdrafts will originate in this cloud. Since this may   !
         !    lead to an impossible cloud, check for success.                              !
         !---------------------------------------------------------------------------------!
         call grell_find_downdraft_origin(mkx,mgmzp,k22,ktop,relheight_down,zcutdown,z_cup &
                                         ,theivs_cup,dzd_cld,ierr,kdet,jmin)
         if (ierr /= 0) cycle stacloop

         !---------------------------------------------------------------------------------!
         ! 3. Now that we know kdet, change detrainment rate below this level.             !
         !---------------------------------------------------------------------------------!
         do k=kdet,1,-1
            cdd(k) = mentrd_rate(k)+(1. - (pmass_kept*z_cup(k) + pmass_left*z_cup(kdet))   &
                   / (pmass_kept*z_cup(k+1) + pmass_left*z_cup(kdet)) )/dzd_cld(k)
         end do


         !---------------------------------------------------------------------------------!
         ! 4. Normalized mass fluxes.                                                      !
         !---------------------------------------------------------------------------------!
         call grell_nms_downdraft(mkx,mgmzp,kdet,jmin,mentrd_rate,cdd,z_cup,dzd_cld        &
                                 ,etad_cld)

         !---------------------------------------------------------------------------------!
         ! 5. Moist static energy. This may prevent convection to happen, check for        !
         !    success.                                                                     !
         !---------------------------------------------------------------------------------!
         call grell_theiv_downdraft(mkx,mgmzp,jmin,cdd,mentrd_rate,theiv,theiv_cup         &
                                   ,theivs_cup,dzd_cld,ierr,theivd_cld)
         if (ierr /= 0) cycle stacloop

         !---------------------------------------------------------------------------------!
         ! 6. Moisture properties of downdraft                                             !
         !---------------------------------------------------------------------------------!
         call grell_most_thermo_downdraft(mkx,mgmzp,jmin,qtot,mentrd_rate,cdd,p_cup        &
                                         ,exner_cup,thil_cup,t_cup,qtot_cup,qvap_cup       &
                                         ,qliq_cup,qice_cup,qsat_cup,rho_cup,pwav          &
                                         ,theivd_cld,etad_cld,dzd_cld,thild_cld,td_cld     &
                                         ,qtotd_cld,qvapd_cld,qliqd_cld,qiced_cld          &
                                         ,qsatd_cld,rhod_cld,dbyd,pwd_cld,pwev,ierr)
         if (ierr /= 0) cycle stacloop
         !---------------------------------------------------------------------------------!
         ! 7. Compute the downdraft strength in terms of windshear. Remembering that edt   !
         !    is the fraction between downdraft and updraft.                               !
         !---------------------------------------------------------------------------------!
         call grell_efficiency_ensemble(mkx,mgmzp,maxens_eff,k22,kbcon,ktop,edtmin,edtmax  &
                                       ,pwav,pwev,z_cup,uwind,vwind,dzd_cld,edt_eff(1,icap))

         !---------------------------------------------------------------------------------!
         ! 8. Checking for water availability and evaporation consistency: we assume that  !
         !    downdraft is always saturated, and it gets the moisture from the rain. How-  !
         !    ever, it cannot require more rain than what is available, so if that would   !
         !    be the case, convection will be skipped this time. This is probably more     !
         !    strict than what we would need. Grell (1993) for example, would still allow  !
         !    the cloud to exist, but without downdraft, which may be too soft, since      !
         !    precipitation would increase in an ambient that is probably too dry to       !
         !    produce precipitation (my guess is that a mid-term, a non-precipitating      !
         !    cloud would be probably the most correct way to deal with it).               !
         !---------------------------------------------------------------------------------!
         if (pwav + edt_eff(1,icap) * pwev < 0. ) then
            ierr = 9
            cycle stacloop
         end if

         !---------------------------------------------------------------------------------!
         ! 9. Computing cloud work function associated with downdrafts                     !
         !---------------------------------------------------------------------------------!
         call grell_cldwork_downdraft(mkx,mgmzp,jmin,dbyd,dzd_cld,etad_cld,aad)

      else
         !---------------------------------------------------------------------------------!
         !    Now we assume that if a cloud doesn't have dowdrafts, it doesn't precipitate !
         ! either. In the future we may revisit this assumption. For now, if the cloud is  !
         ! shallow, then we set up some of the above parameters for non existent dowdraft  !
         ! (i.e. no mass flux, downdraft thermodynamic variables coincident to cloud       !
         ! levels etc.). kdet and jmin were already initialized to zero.                   !
         !---------------------------------------------------------------------------------!
         etad_cld     = 0.
         dbyd         = 0.
         thild_cld    = thil_cup
         theivd_cld   = theiv_cup
         qtotd_cld    = qtot_cup
         qliqd_cld    = qliq_cup
         qiced_cld    = qice_cup
         pwd_cld      = 0.
         pwev         = 0.
         edt_eff      = 0.
         aad          = 0.
      end if
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!


      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! K. Finding the non-forced cloud work, which will require computing most            !
      !    subroutines again. This will go and compute even if the tests would prevent     !
      !    the cloud to happen. This part will be skipped if the user is asking for moist- !
      !    ure convergence only.                                                           !
      !------------------------------------------------------------------------------------!
      if (comp_noforc_cldwork) then
         !---------------------------------------------------------------------------------!
         ! i.     Initialising error flag                                                  !
         !---------------------------------------------------------------------------------!
         ierr0   = 0
      
         !---------------------------------------------------------------------------------!
         ! ii.    Initialise detrainment and cloud work variables.                         !
         !---------------------------------------------------------------------------------!
         aa0d   = 0.
         aa0u   = 0.
         aatot0 = 0.

         !---------------------------------------------------------------------------------!
         ! iii.   Calculate all thermodynamic properties at the cloud level and initialize !
         !        draft thermodynamic properties                                           !
         !---------------------------------------------------------------------------------!
         call grell_thermo_cldlev(mkx,mgmzp,z_cup,exner0,thil0,t0,qtot0,qliq0,qice0        &
                                 ,exnersur,thilsur,tsur,qtotsur,qliqsur,qicesur,exner0_cup &
                                 ,p0_cup,t0_cup,thil0_cup,qtot0_cup,qvap0_cup,qliq0_cup    &
                                 ,qice0_cup,qsat0_cup,rho0_cup,theiv0_cup,theivs0_cup)

         !---------------------------------------------------------------------------------!
         ! iv.    Calculate the thermodynamic properties below the level of free convec-   !
         !        tion.                                                                    !
         !---------------------------------------------------------------------------------!
         call grell_buoy_below_lfc(mkx,mgmzp,k22,kbcon,exner0_cup,p0_cup,theiv0_cup        &
                                  ,thil0_cup,t0_cup,qtot0_cup,qvap0_cup,qliq0_cup          &
                                  ,qice0_cup,qsat0_cup,rho0_cup,theiv0u_cld,thil0u_cld     &
                                  ,t0u_cld,qtot0u_cld,qvap0u_cld,qliq0u_cld,qice0u_cld     &
                                  ,qsat0u_cld,rho0u_cld,dby0u)

         !---------------------------------------------------------------------------------!
         ! v.     Calculate the incloud ice-vapour equivalent potential temperature        !
         !---------------------------------------------------------------------------------!
         call grell_theiv_updraft(mkx,mgmzp,k22,kbcon,cdu,mentru_rate,theiv0,theiv0_cup    &
                                 ,dzu_cld,theiv0u_cld)

         !---------------------------------------------------------------------------------!
         ! vi.    Finding the moisture profiles associated with updrafts.                  !
         !---------------------------------------------------------------------------------!
         call grell_most_thermo_updraft(comp_down,mkx,mgmzp,kbcon,ktop,cdu,mentru_rate     &
                                       ,qtot0,p0_cup,exner0_cup,thil0_cup,t0_cup,qtot0_cup &
                                       ,qvap0_cup,qliq0_cup,qice0_cup,qsat0_cup,rho0_cup   &
                                       ,theiv0u_cld,etau_cld,dzu_cld,thil0u_cld,t0u_cld    &
                                       ,qtot0u_cld,qvap0u_cld,qliq0u_cld,qice0u_cld        &
                                       ,qsat0u_cld,rho0u_cld,dby0u,pw0u_cld,pwav0)
         !---------------------------------------------------------------------------------!
         ! vii.   Finding the cloud work function                                          !
         !---------------------------------------------------------------------------------!
         call grell_cldwork_updraft(mkx,mgmzp,kbcon,ktop,dby0u,dzu_cld,etau_cld,aa0u)
         
         !----- Downdraft properties ------------------------------------------------------!
         if (comp_down) then

            !------------------------------------------------------------------------------!
            ! viii.  Moist static energy.                                                  !
            !------------------------------------------------------------------------------!
            call grell_theiv_downdraft(mkx,mgmzp,jmin,cdd,mentrd_rate,theiv0,theiv0_cup    &
                                      ,theivs0_cup,dzd_cld,ierr0,theiv0d_cld)
            !------------------------------------------------------------------------------!
            ! ix.    Moisture properties of downdraft                                      !
            !------------------------------------------------------------------------------!
            call grell_most_thermo_downdraft(mkx,mgmzp,jmin,qtot0,mentrd_rate,cdd,p0_cup   &
                                            ,exner0_cup,thil0_cup,t0_cup,qtot0_cup         &
                                            ,qvap0_cup,qliq0_cup,qice0_cup,qsat0_cup       &
                                            ,rho0_cup,pwav0,theiv0d_cld,etad_cld,dzd_cld   &
                                            ,thil0d_cld,t0d_cld,qtot0d_cld,qvap0d_cld      &
                                            ,qliq0d_cld,qice0d_cld,qsat0d_cld,rho0d_cld    &
                                            ,dby0d,pw0d_cld,pwev0,ierr0)

            !------------------------------------------------------------------------------!
            ! x.     Downdraft cloud work function.                                        !
            !------------------------------------------------------------------------------!
            call grell_cldwork_downdraft(mkx,mgmzp,jmin,dby0d,dzd_cld,etad_cld,aa0d)

         end if

      end if
      !------------------------------------------------------------------------------------!
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!



      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! L. Now we perform the big loop across the precipitation efficiency ensembles.      !
      !------------------------------------------------------------------------------------!
      effloop: do iedt=1,maxens_eff
         
         !---------------------------------------------------------------------------------!
         ! 1. Copying the current epsilon to a shortcut. This is an output variable, how-  !
         !    ever, it will be used as scratch here. At the feedback time the value will   !
         !    be overwritten by the actual output value.                                   !
         !---------------------------------------------------------------------------------!
         edt  = edt_eff(iedt,icap)

         !---------------------------------------------------------------------------------!
         ! 2. Compute the total cloud work function for this edt member                    !
         !---------------------------------------------------------------------------------!
         if (comp_noforc_cldwork) aatot0 = aa0u + edt * aa0d
         aatot  = aau + edt * aad


         !---------------------------------------------------------------------------------!
         ! 3. Compute the changes in thetae_iv,theta_il, and total mixing ratio at the     !
         !    bottom, considering downdraft effects. If this cloud doesn't have down-      !
         !    drafts, it should be zero so skip it.                                        !
         !---------------------------------------------------------------------------------!
         if (comp_down) then
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,theiv,p_cup          &
                                         ,theiv_cup,mentrd_rate,cdd,dzd_cld,etad_cld       &
                                         ,theivd_cld,dellatheiv_eff(1,iedt,icap))
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,qtot,p_cup,qtot_cup  &
                                         ,mentrd_rate,cdd,dzd_cld,etad_cld,qtotd_cld       &
                                         ,dellaqtot_eff(1,iedt,icap))
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,thil,p_cup,thil_cup  &
                                         ,mentrd_rate,cdd,dzd_cld,etad_cld,thild_cld       &
                                         ,dellathil_eff(1,iedt,icap))
         end if
         
         !---------------------------------------------------------------------------------!
         ! 4. Compute the changes in thetae_iv and total mixing ratio at other levels. In  !
         !    case downdrafts are absent, this will work fine too, since the mass fluxes   !
         !    will be all zero.                                                            !
         !---------------------------------------------------------------------------------!
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,k22,kbcon,jmin,ktop   &
                                   ,theiv,p_cup,theiv_cup,mentrd_rate,mentru_rate,cdd,cdu  &
                                   ,dzd_cld,etad_cld,etau_cld,theivd_cld,theivu_cld        &
                                   ,dellatheiv_eff(1:mgmzp,iedt,icap))
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,k22,kbcon,jmin,ktop   &
                                   ,qtot,p_cup,qtot_cup,mentrd_rate,mentru_rate,cdd,cdu    &
                                   ,dzd_cld,etad_cld,etau_cld,qtotd_cld,qtotu_cld          &
                                   ,dellaqtot_eff(1:mgmzp,iedt,icap))
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,k22,kbcon,jmin,ktop   &
                                   ,thil,p_cup,thil_cup,mentrd_rate,mentru_rate,cdd,cdu    &
                                   ,dzd_cld,etad_cld,etau_cld,thild_cld,thilu_cld          &
                                   ,dellathil_eff(1:mgmzp,iedt,icap))


         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write(unit=16,fmt='(a)') '-------------------------------------------------------'
         !write(unit=16,fmt='(2(a,1x,i5,1x))') '   Densities: iedt=',iedt,'icap=',icap
         !write(unit=16,fmt='(4(a,1x,f10.4,1x))') '     k22      =',   z(k22)               &
         !                                            ,'kbcon    =',   z(kbcon)             &
         !                                            ,'jmin     =',   z(jmin)              &
         !                                            ,'ktop     =',   z(ktop)
         !write(unit=16,fmt='(a5,2(3x,5(1x,a11)))')                                         &
         !   'Level','   RHO0_CUP','  RHO0D_CLD','  RHO0U_CLD','      DBY0D','      DBY0U'  &
         !          ,'    RHO_CUP','   RHOD_CLD','   RHOU_CLD','       DBYD','       DBYU'
         !do k=1,mkx
         !   write(unit=16,fmt='(i5,2(3x,5(1x,f11.5)))')                                    &
         !       k,rho0_cup(k),rho0d_cld(k),rho0u_cld(k),dby0d(k),dby0u(k)                  &
         !        ,rho_cup (k),rhod_cld (k),rhou_cld (k),dbyd (k),dbyu (k)
         !end do
         !write(unit=16,fmt='(a)') '-------------------------------------------------------'
         !write(unit=16,fmt='(a)') ' '
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         !---------------------------------------------------------------------------------!
         ! 5. Computing the rainfall flux, which will be simply all the fall-out water     !
         !    that wasn't consumed by the downdrafts.                                      !
         !---------------------------------------------------------------------------------!
         do k=1,ktop
            !------ Precipitation ---------------------------------------------------------!
            pw_eff(k,iedt,icap) =  pwu_cld(k) + edt * pwd_cld(k)
         end do

         !---------------------------------------------------------------------------------!
         ! 6. This is another big loop, across the large scale force members. We modify    !
         !    the thermodynamic profile with a known mass flux (small one) and check how   !
         !    much change in the cloud work we get. Then we can use this scale to assess   !
         !    what mass flux the convection should have.                                   !
         !---------------------------------------------------------------------------------!
         mbprimeloop: do imbp=maxens_lsf,1,-1

            !------------------------------------------------------------------------------!
            ! 6a. Initialize some ensemble-related scratch variables. mbprime is an        !
            !     arbitrary extra mass flux that is applied to the forcing, to get the     !
            !     upward mass flux for all cloud work related parametrisations. This is    !
            !     also a variable that we perturb to get different members for the         !
            !     ensemble.                                                                !
            !------------------------------------------------------------------------------!
            mboffset = (-1)**imbp * imbp /2 !----- This gives 0, 1, -1, 2, -2, ... --------!
            mbprime=(4.+real(mboffset))*1.e-3

            !------------------------------------------------------------------------------!
            ! 6b. Defining the Kuo (1974) moistening factor b. Following Krishnamurti      !
            !     et al. (1983), this value is around 0.3, with significant oscillations   !
            !     between 0.0 and 0.8. Here we will play it closer to the average,         !
            !     oscillating it                                                           !
            !     between 0.2 and 0.4, equivalent to oscillate (1-b) between 0.6 and 0.8.  !
            !------------------------------------------------------------------------------!
            one_b = 0.7+0.2*real(mboffset)

            !------------------------------------------------------------------------------!
            ! 6c. Applying the fudging to some thermodynamic variables                     !
            !------------------------------------------------------------------------------!

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            ! write(unit=30+mynum,fmt='(a)')                '-------------------------------------------------------------------------------------------------'
            ! write(unit=30+mynum,fmt='(3(a,1x,f10.4,1x))') '§§§§ mbprime  =',  mbprime,'edt     =',        edt,'tscal_kf =',  tscal_kf
            ! write(unit=30+mynum,fmt='(3(a,1x,f10.4,1x))') '     dens_curr=',dens_curr,'cap_max =',    cap_max,'dtime    =',     dtime
            ! write(unit=30+mynum,fmt='(4(a,1x,i10,1x))')   '     k22      =',      k22,'kbcon   =',      kbcon,'jmin     =',      jmin,'ktop    =',      ktop
            ! write(unit=30+mynum,fmt='(3(a,1x,f10.4,1x))') '     aa0d     =',     aa0d,'aa0u    =',       aa0u,'aatot0   =',    aatot0
            ! write(unit=30+mynum,fmt='(3(a,1x,f10.4,1x))') '     aad      =',      aad,'aau     =',        aau,'aatot    =',     aatot
            ! write(unit=30+mynum,fmt='(a)'               ) ' '
            ! write(unit=30+mynum,fmt='(a5,1x,6(a12,1x))') '    K','       THEIV','        QTOT','        THIL','    DELTHEIV','     DELQTOT','     DELTHIL'
            ! do k=1,mkx-1
            !    write(unit=30+mynum,fmt='(i5,1x,6(es12.5,1x))') k,qtot(k),qtot_cup(k),qtotd_cld(k),
            ! end do
            ! write(unit=30+mynum,fmt='(a)')                '-------------------------------------------------------------------------------------------------'
            ! write(unit=30+mynum,fmt='(a)')                ' '
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            do k=1,mkx-1
               x_theiv(k) = theiv(k)          + mbprime*dtime * dellatheiv_eff(k,iedt,icap)
               x_qtot(k)  = max(toodry,qtot(k)+ mbprime*dtime * dellaqtot_eff(k,iedt,icap))
               x_thil(k)  = thil(k)           + mbprime*dtime * dellathil_eff(k,iedt,icap)
            end do
            !----- Boundary condition at mkx, just use environment ------------------------!
            x_theiv (mkx) = theiv (mkx)
            x_qtot  (mkx) = qtot  (mkx)
            x_thil  (mkx) = thil  (mkx)

            !------------------------------------------------------------------------------!
            ! 6d. Initialise some variables.                                               !
            !------------------------------------------------------------------------------!
            x_ierr   = 0
            !----- Cloud work -------------------------------------------------------------!
            x_aad    = 0.
            x_aau    = 0.
            !----- Ensemble variables -----------------------------------------------------!
            x_aatot  = 0.
        
            !------------------------------------------------------------------------------!
            ! OBS: Contrary to before, I won't return or cycle in case convection fails    !
            !      here.                                                                   !
            !------------------------------------------------------------------------------!
            modif_comp_if: if (comp_modif_thermo) then
               !---------------------------------------------------------------------------!
               ! 6e. Compute the modified structure, finding the consistent set and then   !
               !     interpolating them to the cloud levels.                               !
               !---------------------------------------------------------------------------!
               do k=1,mkx
                  x_t(k)    = t   (k)
                  call thil2tqall(x_thil(k),exner(k),p(k),x_qtot(k),x_qliq(k),x_qice(k)    &
                                 ,x_t(k),x_qvap(k),scrvar(k))
               end do
               call grell_thermo_cldlev(mkx,mgmzp,z_cup,exner,x_thil,x_t,x_qtot,x_qliq     &
                                       ,x_qice,exnersur,thilsur,tsur,qtotsur,qliqsur       &
                                       ,qicesur,x_exner_cup,x_p_cup,x_t_cup,x_thil_cup     &
                                       ,x_qtot_cup,x_qvap_cup,x_qliq_cup,x_qice_cup        &
                                       ,x_qsat_cup,x_rho_cup,x_theiv_cup,x_theivs_cup      )

               !---------------------------------------------------------------------------!
               ! 6f. Finding the updraft thermodynamics between the updraft origin and the !
               !     level of free convection.                                             !
               !---------------------------------------------------------------------------!
               call grell_buoy_below_lfc(mkx,mgmzp,k22,kbcon,x_exner_cup,x_p_cup           &
                                        ,x_theiv_cup,x_thil_cup,x_t_cup,x_qtot_cup         &
                                        ,x_qvap_cup,x_qliq_cup,x_qice_cup,x_qsat_cup       &
                                        ,x_rho_cup,x_theivu_cld,x_thilu_cld,x_tu_cld       &
                                        ,x_qtotu_cld,x_qvapu_cld,x_qliqu_cld,x_qiceu_cld   &
                                        ,x_qsatu_cld,x_rhou_cld,x_dbyu)

               !---------------------------------------------------------------------------!
               ! 6g. Compute the updraft profiles ice-vapour equivalent potential temper-  !
               !     ature.                                                                !
               !---------------------------------------------------------------------------!
               call grell_theiv_updraft(mkx,mgmzp,k22,kbcon,cdu,mentru_rate,x_theiv        &
                                       ,x_theiv_cup,dzu_cld,x_theivu_cld)

               !---------------------------------------------------------------------------!
               ! 6h. Getting the updraft moisture profile                                  !
               !---------------------------------------------------------------------------!
               call grell_most_thermo_updraft(comp_down,mkx,mgmzp,kbcon,ktop,cdu           &
                                             ,mentru_rate,x_qtot,x_p_cup,x_exner_cup       &
                                             ,x_thil_cup,x_t_cup,x_qtot_cup,x_qvap_cup     &
                                             ,x_qliq_cup,x_qice_cup,x_qsat_cup,x_rho_cup   &
                                             ,x_theivu_cld,etau_cld,dzu_cld,x_thilu_cld    &
                                             ,x_tu_cld,x_qtotu_cld,x_qvapu_cld,x_qliqu_cld &
                                             ,x_qiceu_cld,x_qsatu_cld,x_rhou_cld,x_dbyu    &
                                             ,x_pwu_cld,x_pwav)

               !---------------------------------------------------------------------------!
               ! 6i. Recalculating the updraft cloud work                                  !
               !---------------------------------------------------------------------------!
               call grell_cldwork_updraft(mkx,mgmzp,kbcon,ktop,x_dbyu,dzu_cld,etau_cld     &
                                         ,x_aau)

               modif_down: if (comp_down) then

                  !------------------------------------------------------------------------!
                  ! 6j. Finding moist static energy                                        !
                  !------------------------------------------------------------------------!
                  call grell_theiv_downdraft(mkx,mgmzp,jmin,cdd,mentrd_rate,x_theiv        &
                                            ,x_theiv_cup,x_theivs_cup,dzd_cld,x_ierr       &
                                            ,x_theivd_cld)

                  !------------------------------------------------------------------------!
                  ! 6k. Moisture properties                                                !
                  !------------------------------------------------------------------------!
                  call grell_most_thermo_downdraft(mkx,mgmzp,jmin,x_qtot,mentrd_rate,cdd   &
                                                  ,x_p_cup,x_exner_cup,x_thil_cup,x_t_cup  &
                                                  ,x_qtot_cup,x_qvap_cup,x_qliq_cup        &
                                                  ,x_qice_cup,x_qsat_cup,x_rho_cup,x_pwav  &
                                                  ,x_theivd_cld,etad_cld,dzd_cld           &
                                                  ,x_thild_cld,x_td_cld,x_qtotd_cld        &
                                                  ,x_qvapd_cld,x_qliqd_cld,x_qiced_cld     &
                                                  ,x_qsatd_cld,x_rhod_cld,x_dbyd,x_pwd_cld &
                                                  ,x_pwev,x_ierr)

                  !------------------------------------------------------------------------!
                  ! 6l. Computing cloud work function associated with downdrafts           !
                  !------------------------------------------------------------------------!
                  call grell_cldwork_downdraft(mkx,mgmzp,jmin,x_dbyd,dzd_cld,etad_cld,x_aad)
               end if modif_down
            end if modif_comp_if
            !------------------------------------------------------------------------------!
            ! 6m. If the fudged field produced a cloud, we compute the total cloud work    !
            !     function for this particular test. If not, or if the user is using       !
            !     moisture convergence closure only, this will be a silly 0 + 0.           !
            !------------------------------------------------------------------------------!
            x_aatot   = x_aau + edt * x_aad

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write(unit=60+mynum,fmt='(a)')                '-------------------------------------------------------------------------------------------------'
            !write(unit=60+mynum,fmt='(3(a,1x,f10.4,1x))') '§§§§ mbprime  =',  mbprime,'edt     =',        edt,'tscal_kf =',  tscal_kf
            !write(unit=60+mynum,fmt='(3(a,1x,f10.4,1x))') '     dens_curr=',dens_curr,'cap_max =',    cap_max,'dtime    =',     dtime
            !write(unit=60+mynum,fmt='(4(a,1x,f10.4,1x))') '     k22      =',   z(k22),'kbcon   =',   z(kbcon),'jmin     =',   z(jmin),'ktop    =',   z(ktop)
            !write(unit=60+mynum,fmt='(3(a,1x,f10.4,1x))') '     aa0d     =',     aa0d,'aa0u    =',       aa0u,'aatot0   =',    aatot0
            !write(unit=60+mynum,fmt='(3(a,1x,f10.4,1x))') '     aad      =',      aad,'aau     =',        aau,'aatot    =',     aatot
            !write(unit=60+mynum,fmt='(3(a,1x,f10.4,1x))') '     x_aad    =',    x_aad,'x_aau   =',      x_aau,'x_aatot  =',   x_aatot
            !write(unit=60+mynum,fmt='(a)'               ) ' '
            !write(unit=60+mynum,fmt='(a)'               ) '    Dynamic control:'
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!



            !------------------------------------------------------------------------------!
            ! 6n. Running the dynamic control ensemble, using the different combinations   !
            !     of cloud efficiency and fudged variables. The reference upward mass flux !
            !     will then have maxens_eff × maxens_lsf × maxens_dyn different values.    !
            !------------------------------------------------------------------------------!
            call grell_dyncontrol_ensemble(closure_type,comp_down,mgmzp,maxens_dyn,dtime   &
                                          ,tscal_kf,dens_curr,prev_dnmf,mconv,k22,kbcon    &
                                          ,ktop,omeg,x_p_cup,edt,mbprime,one_b,aatot0      &
                                          ,aatot,x_aatot,pwav,pwev                         &
                                          ,dnmf_ens(1,imbp,iedt,icap)                      &
                                          ,upmf_ens(1,imbp,iedt,icap))
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write(unit=60+mynum,fmt='(2(a,1x,f10.4,1x))') '     dnmf   =',dnmf_ens(1,imbp,iedt,icap),'upmf   =',upmf_ens(1,imbp,iedt,icap)
            !write(unit=60+mynum,fmt='(a)')                '-------------------------------------------------------------------------------------------------'
            !write(unit=60+mynum,fmt='(a)'               ) ' '
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         end do mbprimeloop
      end do effloop
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!
   end do stacloop

   !---------------------------------------------------------------------------------------!
   ! 7. Now that all ensemble members are computed, we will proceed with the feedback      !
   !    to the large scale. If ierr is activated, then I don't need to compute this step.  !
   !---------------------------------------------------------------------------------------!
   if (ierr == 0) then 
      call grell_feedback(comp_down,mgmzp,maxens_cap,maxens_eff,maxens_lsf,maxens_dyn      &
                         ,inv_ensdim,max_heat,ktop,edt_eff,dellathil_eff,dellaqtot_eff     &
                         ,pw_eff,dnmf_ens,upmf_ens,ierr,upmf,dnmf,edt,outthil,outqtot      &
                         ,precip)
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine grell_cupar_main
!==========================================================================================!
!==========================================================================================!

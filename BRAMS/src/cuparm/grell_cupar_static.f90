!==========================================================================================!
!==========================================================================================!
!    This subroutine will compute all the variables that depend on the static control, and !
! create the arrays of ensemble members that depend on precipitation efficiency and the    !
! cap_maxs variation.  This subroutine works for any cloud, and the only difference is     !
! that if this is a shallow cloud, then it will bypass everything related to downdrafts.   !
! Whenever applicable, the follow convention is used:                                      !
!       [+] "?_cup"  : These are variables defined at the Grell level (staggered).         !
!       [+] "?d_cld" : These are cloud variables associated with downdrafts.               !
!       [+] "?u_cld" : These are cloud variables associated with updrafts.                 !
!       [+] "x_?"    : Variables modified by an arbitrary mass flux, for ensemble          !
!                      statistics and interaction between clouds.                          !
!       [+] "?0?"    : Whenever a variable like the ones above contains a 0, it means that !
!                      these variables are based on current values. If this is a BRAMS     !
!                      variable, it means the value before applying the tendency,          !
!                      otherwise it is a derived value based on the current large scale    !
!                      These are used only if Grell (1993) or Kain-Fritsch (1990)          !
!                      parameterizations are used.                                         !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine grell_cupar_static(comp_noforc_cldwork,checkmass,iupmethod,maxens_cap           &
                             ,maxens_eff,mgmzp,cap_maxs,cap_maxs_increment,wnorm_max       &
                             ,wnorm_increment,depth_min,depth_max,edtmax,edtmin,masstol    &
                             ,pmass_left,radius,relheight_down,zkbmax,zcutdown,z_detr      &
                             ,cld2prec,prec_cld,aad,aau,dellatheiv_eff,dellathil_eff       &
                             ,dellaqtot_eff,dellaco2_eff,pw_eff,edt_eff,aatot0_eff         &
                             ,aatot_eff,ierr_cap,comp_down_cap,klod_cap,klou_cap,klcl_cap  &
                             ,klfc_cap,kdet_cap,kstabi_cap,kstabm_cap,klnb_cap,ktop_cap    &
                             ,pwav_cap,pwev_cap,wbuoymin_cap,cdd_cap,cdu_cap               &
                             ,mentrd_rate_cap,mentru_rate_cap,dbyd_cap,dbyu_cap            &
                             ,etad_cld_cap,etau_cld_cap,rhod_cld_cap,rhou_cld_cap          &
                             ,qliqd_cld_cap,qliqu_cld_cap,qiced_cld_cap,qiceu_cld_cap,i,j  &
                             ,icld,mynum)
   use grid_dims        , only : &
       str_len      ! ! intent(in) - Typical string length. 
   use mem_ensemble     , only : &
       ensemble_vars             ! ! structure - The ensemble scratch structure. ----------!

   use mem_scratch_grell, only : &
      !----- Variables defined at the initialization process ------------------------------!
       co20         & ! intent(in)  - Current CO2 mixing ratio                    [    ppm]
      ,co2          & ! intent(in)  - CO2 mixing ratio with forcing.              [    ppm]
      ,co2sur       & ! intent(in)  - surface CO2 mixing ratio                    [    ppm]
      ,dzd_cld      & ! intent(in)  - Delta-height for downdraft calculations     [      m]
      ,dzu_cld      & ! intent(in)  - Delta-height for updraft calculations       [      m]
      ,exner0       & ! intent(in)  - Current Exner function                      [ J/kg/K]
      ,exner        & ! intent(in)  - Forced Exner funtion                        [ J/kg/K]
      ,exnersur     & ! intent(in)  - Surface Exner function                      [ J/kg/K]
      ,kpbl         & ! intent(in)  - Level of PBL top (Nakanishi/Niino only)     [    ---]
      ,ktpse        & ! intent(in)  - Maximum level allowed for cloud top         [    ---]
      ,mkx          & ! intent(in)  - Number of vertical levels                   [    ---]
      ,p0           & ! intent(in)  - Current pressure                            [     Pa]
      ,p            & ! intent(in)  - Pressure with forcing                       [     Pa]
      ,psur         & ! intent(in)  - Surface pressure                            [     Pa]
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
      ,rho          & ! intent(in)  - Air density                                 [  kg/m設
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
      ,sigw         & ! intent(in)  - Vertical velocity standard deviation        [    m/s]
      ,uwind        & ! intent(in)  - Zonal wind                                  [    m/s]
      ,vwind        & ! intent(in)  - Meridional wind                             [    m/s]
      ,wwind        & ! intent(in)  - Mean vertical velocity                      [    m/s]
      ,z            & ! intent(in)  - Height                                      [      m]
      ,z_cup        & ! intent(in)  - Height at cloud levels                      [      m]
      ,aa0d         & ! intent(out) - Dndraft Updraft work function               [   J/kg]
      ,aa0u         & ! intent(out) - Updraft work function                       [   J/kg]
      ,cdd          & ! intent(out) - Norm. dndraft detrainment rate              [    ---]
      ,cdu          & ! intent(out) - Norm. updraft detrainment rate              [    ---]
      ,co20_cup     & ! intent(out) - CO2 mixing ratio                            [    ppm]
      ,co2_cup      & ! intent(out) - CO2 mixing ratio with forcing               [    ppm]
      ,co2d_cld     & ! intent(out) - Dndraft CO2 mixing ratio w/ forcing         [    ppm]
      ,co20d_cld    & ! intent(out) - Dndraft CO2 mixing ratio                    [    ppm]
      ,co20u_cld    & ! intent(out) - Updraft CO2 water mixing ratio              [    ppm]
      ,co2u_cld     & ! intent(out) - Updraft CO2 water mixing ratio w/ forcing   [    ppm]
      ,comp_dn      & ! intent(out) - Downdraft thermodynamics is computed.       [    ---]
      ,dby0d        & ! intent(out) - Dndraft Buoyancy acceleration               [   m/s淫
      ,dbyd         & ! intent(out) - Dndraft Buoyancy acceleration w/ forcing    [   m/s淫
      ,dby0u        & ! intent(out) - Updraft Buoyancy acceleration               [   m/s淫
      ,dbyu         & ! intent(out) - Updraft Buoyancy acceleration w/ forcing    [   m/s淫
      ,etad_cld     & ! intent(out) - Normalised downdraft mass flux w/ forcing   [    ---]
      ,etau_cld     & ! intent(out) - Normalised updraft mass flux                [    ---]
      ,exner0_cup   & ! intent(out) - Exner function                              [   J/kg]
      ,exner_cup    & ! intent(out) - Exner function with forcing                 [   J/kg]
      ,ierr         & ! intent(out) - Flag for convection error                   [    ---]
      ,klod         & ! intent(out) - Level of origin of downdrafts               [    ---]
      ,klou         & ! intent(out) - Level of origin of updrafts                 [    ---]
      ,klcl         & ! intent(out) - Lifting condensation level                  [    ---]
      ,klfc         & ! intent(out) - Level of free convection                    [    ---]
      ,kdet         & ! intent(out) - Top of downdraft detrainemnt layer          [    ---]
      ,kstabi       & ! intent(out) - cloud stable layer base                     [    ---]
      ,kstabm       & ! intent(out) - cloud stable layer top                      [    ---]
      ,klnb         & ! intent(out) - level of neutral buoyancy                   [    ---]
      ,ktop         & ! intent(out) - cloud top                                   [    ---]
      ,mentrd_rate  & ! intent(out) - Norm. dndraft entrainment rate w/ forcing   [    ---]
      ,mentru_rate  & ! intent(out) - Norm. updraft entrainment rate              [    ---]
      ,p0_cup       & ! intent(out) - Pressure                                    [     Pa]
      ,p_cup        & ! intent(out) - Pressure with forcing                       [     Pa]
      ,pw0d_cld     & ! intent(out) - Dndraft Condensation                        [  kg/kg]
      ,pwd_cld      & ! intent(out) - Dndraft Condensation w/ forcing             [  kg/kg]
      ,pw0u_cld     & ! intent(out) - Updraft Condensation                        [  kg/kg]
      ,pwu_cld      & ! intent(out) - Updraft Condensation w/ forcing             [  kg/kg]
      ,pwav0        & ! intent(out) - Updraft Integrated condensation             [  kg/kg]
      ,pwav         & ! intent(out) - Dndraft Integrated condensation w/ forcing  [  kg/kg]
      ,pwev0        & ! intent(out) - Dndraft Integrated evaporation              [  kg/kg]
      ,pwev         & ! intent(out) - Updraft Integrated evaporation  w/ forcing  [  kg/kg]
      ,qtot0_cup    & ! intent(out) - Total mixing ratio                          [  kg/kg]
      ,qtot_cup     & ! intent(out) - Total mixing ratio with forcing             [  kg/kg]
      ,qtot0d_cld   & ! intent(out) - Dndraft Total water mixing ratio            [  kg/kg]
      ,qtotd_cld    & ! intent(out) - Dndraft Total water mixing ratio w/ forcing [  kg/kg]
      ,qtot0u_cld   & ! intent(out) - Updraft Total water mixing ratio            [  kg/kg]
      ,qtotu_cld    & ! intent(out) - Updraft Total water mixing ratio w/ forcing [  kg/kg]
      ,qvap0_cup    & ! intent(out) - Water vapour mixing ratio                   [  kg/kg]
      ,qvap_cup     & ! intent(out) - Water vapour mixing ratio with forcing      [  kg/kg]
      ,qvap0d_cld   & ! intent(out) - Dndraft Vapour mixing ratio                 [  kg/kg]
      ,qvapd_cld    & ! intent(out) - Dndraft Vapour mixing ratio w/ forcing      [  kg/kg]
      ,qvap0u_cld   & ! intent(out) - Updraft Vapour mixing ratio                 [  kg/kg]
      ,qvapu_cld    & ! intent(out) - Updraft Vapour mixing ratio w/ forcing      [  kg/kg]
      ,qliq0_cup    & ! intent(out) - Liquid water                                [  kg/kg]
      ,qliq_cup     & ! intent(out) - Liquid water with forcing                   [  kg/kg]
      ,qliq0d_cld   & ! intent(out) - Dndraft Liquid water mixing ratio           [  kg/kg]
      ,qliqd_cld    & ! intent(out) - Dndraft Liquid water mix. ratio w/ forcing  [  kg/kg]
      ,qliqu_cld    & ! intent(out) - Updraft Liquid water mix. ratio w/ forcing  [  kg/kg]
      ,qliq0u_cld   & ! intent(out) - Updraft Liquid water mixing ratio           [  kg/kg]
      ,qice0_cup    & ! intent(out) - Mixing ratio                                [  kg/kg]
      ,qice_cup     & ! intent(out) - Mixing ratio with forcing                   [  kg/kg]
      ,qice0d_cld   & ! intent(out) - Dndraft Ice mixing ratio                    [  kg/kg]
      ,qiced_cld    & ! intent(out) - Dndraft Ice mixing ratio w/ forcing         [  kg/kg]
      ,qice0u_cld   & ! intent(out) - Updraft Ice mixing ratio                    [  kg/kg]
      ,qiceu_cld    & ! intent(out) - Updraft Ice mixing ratio w/ forcing         [  kg/kg]
      ,qsat0_cup    & ! intent(out) - Saturation mixing ratio                     [  kg/kg]
      ,qsat_cup     & ! intent(out) - Saturation mixing ratio with forcing        [  kg/kg]
      ,qsat0d_cld   & ! intent(out) - Dndraft Sat. mixing ratio                   [  kg/kg]
      ,qsatd_cld    & ! intent(out) - Dndraft Sat. mixing ratio w/ forcing        [  kg/kg]
      ,qsat0u_cld   & ! intent(out) - Updraft Sat. mixing ratio                   [  kg/kg]
      ,qsatu_cld    & ! intent(out) - Updraft Sat. mixing ratio  w/ forcing       [  kg/kg]
      ,rho0_cup     & ! intent(out) - Density                                     [  kg/m設
      ,rho_cup      & ! intent(out) - Density with forcing                        [  kg/m設
      ,rho0d_cld    & ! intent(out) - Dndraft Density                             [  kg/m設
      ,rhod_cld     & ! intent(out) - Dndraft Density w/ forcing                  [  kg/m設
      ,rho0u_cld    & ! intent(out) - Updraft Density                             [  kg/m設
      ,rhou_cld     & ! intent(out) - Updraft Density w/ forcing                  [  kg/m設
      ,t0_cup       & ! intent(out) - Temperature                                 [      K]
      ,t_cup        & ! intent(out) - Temperature with forcing                    [      K]
      ,t0d_cld      & ! intent(out) - Dndraft Temperature                         [      K]
      ,td_cld       & ! intent(out) - Dndraft Temperature w/ forcing              [      K]
      ,t0u_cld      & ! intent(out) - Updraft Temperature                         [      K]
      ,tu_cld       & ! intent(out) - Updraft Temperature w/ forcing              [      K]
      ,thil0_cup    & ! intent(out) - Ice-liquid potential temperature            [      K]
      ,thil_cup     & ! intent(out) - Ice-liquid pot. temperature w. forcing      [      K]
      ,thil0d_cld   & ! intent(out) - Dndraft Ice-liquid pot. temperature         [      K]
      ,thild_cld    & ! intent(out) - Dndraft Theta-IL with forcing               [      K]
      ,thil0u_cld   & ! intent(out) - Updraft Ice-liquid potential temperature    [      K]
      ,thilu_cld    & ! intent(out) - Updraft THETA-il with forcing               [      K]
      ,theiv0_cup   & ! intent(out) - Ice-vapour potential temperature            [      K]
      ,theiv_cup    & ! intent(out) - Ice-vapour equiv. pot. temp. with forcing   [      K]
      ,theiv0d_cld  & ! intent(out) - Dndraft Ice-vapour equiv. pot. temp.        [      K]
      ,theivd_cld   & ! intent(out) - Dndraft Theta-Eiv with forcing              [      K]
      ,theiv0u_cld  & ! intent(out) - Updraft Ice-vapour equiv. pot. temp.        [      K]
      ,theivu_cld   & ! intent(out) - Updraft THETA-Eiv  with forcing             [      K]
      ,theivs0_cup  & ! intent(out) - Sat. Ice-vapour equiv. pot. temp.           [      K]
      ,theivs_cup   & ! intent(out) - Sat. Ice-vapour equiv. pot. temp. w/ forcing[      K]
      ,wbuoymin0    & ! intent(out) - Updraft Minimum buoyant velocity            [    m/s]
      ,wbuoymin     ! ! intent(out) - Minimum buoyancy velocity                   [    ---]
   use rconstants, only : toodry & ! intent(in)
                        , t00    ! ! intent(in)
   use therm_lib , only : toler
   implicit none
   
   !---------------------------------------------------------------------------------------!
   ! List of arguments                                                                     !
   !---------------------------------------------------------------------------------------!
   !----- Logical flags, to bypass uncessary steps ----------------------------------------!
   logical , intent(in) :: comp_noforc_cldwork ! I will compute no forced cloud work  [T/F]
   logical , intent(in) :: checkmass           ! I will check mass balance            [T/F]
   !----- Integer, global dimensions ------------------------------------------------------!
   integer , intent(in) :: iupmethod         ! Method to find the updraft originatin level
   integer , intent(in) :: maxens_cap        ! Ensemble size on static control (cap_maxs)
   integer , intent(in) :: maxens_eff        ! Ensemble size on precipitation efficiency
   integer , intent(in) :: mgmzp             ! Vertical grid size
   integer , intent(in) :: i                 ! Node X position
   integer , intent(in) :: j                 ! Node Y position
   integer , intent(in) :: icld              ! Cloud type
   integer , intent(in) :: mynum             ! Node ID, for debugging purposes only

   !----- Real, parameters to define the cloud --------------------------------------------!
   real    , intent(in) :: cap_maxs          ! Maximum depth of capping inversion    [ hPa]
   real    , intent(in) :: cap_maxs_increment! Extra cap_maxs due to upstream conv.  [ hPa]
   real    , intent(in) :: depth_min         ! Minimum cloud depth to qualify it     [   m]
   real    , intent(in) :: depth_max         ! Maximum cloud depth to qualify it     [   m]
   real    , intent(in) :: edtmax            ! Maximum epsilon (dnmf/upmf)           [ ---]
   real    , intent(in) :: edtmin            ! Minimum epsilon (dnmf/upmf)           [ ---]
   real    , intent(in) :: masstol           ! Maximum mass leak allowed to happen   [ ---]
   real    , intent(in) :: pmass_left        ! Fraction of mass left at the ground   [ ---]
   real    , intent(in) :: radius            ! Radius, for entrainment rate.         [   m]
   real    , intent(in) :: relheight_down    ! Relative height for downdraft origin  [ ---]
   real    , intent(in) :: wnorm_max         ! Maximum normalised w to be considered [ ---]
   real    , intent(in) :: wnorm_increment   ! Extra wnorm for ensemble              [ ---]
   real    , intent(in) :: zkbmax            ! Top height for updrafts to originate  [   m]
   real    , intent(in) :: zcutdown          ! Top height for dndrafts to originate  [   m]
   real    , intent(in) :: z_detr            ! Top height for downdraft detrainment  [   m]
   real    , intent(in) :: cld2prec          ! Fraction of cloud that becomes prec.  [   m]

   !----- Logical, parameters to define the cloud -----------------------------------------!
   logical , intent(in) :: prec_cld          ! Precipitating cloud.                  [ T|F]

   !----- The ensemble scratch structure variables. ---------------------------------------!
   real   , dimension(mgmzp,maxens_eff,maxens_cap), intent(inout) ::                       &
            dellatheiv_eff  & ! Change in ice-liquid potential temperature  [       K/kg/s]
           ,dellathil_eff   & ! Change in ice-liquid potential temperature  [       K/kg/s]
           ,dellaqtot_eff   & ! Change in total mixing ratio                [   kg/kg/kg/s]
           ,dellaco2_eff    & ! Change in CO2 mixing ratio                  [痠ol/mol/kg/s]
           ,pw_eff          ! ! Water that doesn't evaporate (aka rain).    [   kg/kg/kg/s]
   real   , dimension(maxens_eff,maxens_cap), intent(inout) ::                             &
            edt_eff         & ! Precipitation efficiency for each member    [          ---]
           ,aatot0_eff      & ! Current total cloud work function           [         J/kg]
           ,aatot_eff       ! ! Total cloud work func. (forced)             [         J/kg]
   logical, dimension(maxens_cap), intent(inout) ::                                        &
            comp_down_cap   ! ! I will compute downdrafts.           [T/F]
   integer, dimension(maxens_cap), intent(inout) ::                                        &
            ierr_cap        & ! Convection failure flag.                    [          ---]
           ,klod_cap        & ! Level in which downdrafts originate         [          ---]
           ,klou_cap        & ! Level in which updrafts originate           [          ---]
           ,klcl_cap        & ! Lifting condensation level                  [          ---]
           ,klfc_cap        & ! Level of free convection                    [          ---]
           ,kdet_cap        & ! Top of downdraft detrainemnt layer          [          ---]
           ,kstabi_cap      & ! cloud stable layer base                     [          ---]
           ,kstabm_cap      & ! cloud stable layer top                      [          ---]
           ,klnb_cap        & ! Level of neutral buoyancy                   [          ---]
           ,ktop_cap        ! ! cloud top                                   [          ---]
   real   , dimension(maxens_cap), intent(inout) ::                                        &
            pwav_cap        & ! Integrated condensation                     [        kg/kg]
           ,pwev_cap        & ! Integrated evaporation                      [        kg/kg]
           ,wbuoymin_cap    ! ! Minimum buoyant velocity                    [          m/s]
   real   , dimension(mgmzp,maxens_cap), intent(inout) ::                                  &
            cdd_cap         & ! normalised downdraft detrainment rate       [          ---]
           ,cdu_cap         & ! normalised updraft detrainment rate         [          ---]
           ,mentrd_rate_cap & ! normalised downdraft entrainment rate       [          ---]
           ,mentru_rate_cap & ! normalised updraft entrainment rate         [          ---]
           ,dbyd_cap        & ! Buoyancy associated with downdrafts         [         m/s淫
           ,dbyu_cap        & ! Buoyancy associated with updrafts           [         m/s淫
           ,etad_cld_cap    & ! normalised downdraft mass flux              [          ---]
           ,etau_cld_cap    & ! normalised updraft mass flux                [          ---]
           ,rhod_cld_cap    & ! Downdraft density                           [        kg/m設
           ,rhou_cld_cap    & ! Updraft density                             [        kg/m設
           ,qliqd_cld_cap   & ! Liquid water mixing ratio at downdraft      [        kg/kg]
           ,qliqu_cld_cap   & ! Liquid water mixing ratio at updraft        [        kg/kg]
           ,qiced_cld_cap   & ! Ice mixing ratio at downdraft               [        kg/kg]
           ,qiceu_cld_cap   ! ! Ice mixing ratio at updraft                 [        kg/kg]
   !---------------------------------------------------------------------------------------!

   !----- Output variables, the cloud work functions. Others will be found afterwards. ----!
   real                  , intent(inout) :: aad  ! Downdraft work function        [   J/kg]
   real                  , intent(inout) :: aau  ! Updraft work function          [   J/kg]
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
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !    Flag for convection error when using the "0"  variables.                           !
   !---------------------------------------------------------------------------------------!
   integer :: ierr0     ! Flag for convection error
   integer :: klnb0     ! klnb, but it is just a dummy variable.
   integer :: ktop0     ! ktop, but it is just a dummy variable.
   integer :: iun       ! Unit number, for debugging only
   integer :: l         ! Counter for drawing lines
   logical :: comp_dn0  ! Flag for downdraft
   !---------------------------------------------------------------------------------------!
   !    Miscellaneous parameters                                                           !
   !---------------------------------------------------------------------------------------!
   real :: cap_max    ! Actual "Depth" of inversion capping                       [     Pa]
   real :: wbuoy_max  ! Maximum acceptable buoyant velocity                       [    m/s]
   real :: pmass_kept ! 1.-pmass_left                                             [    ---]
   real :: edt        ! dnmf/upmf                                                 [    ---]
   !----- Scratch array -------------------------------------------------------------------!
   real, dimension(mgmzp) :: scrvar    ! Scratch variable
   !----- String for sanity check. --------------------------------------------------------!
   character(len=str_len) :: which
   !----- Parameter to print debug stuff. -------------------------------------------------!
   logical          , parameter :: print_debug = .false.
   character(len=11)            :: levname
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !  Setting up the mass that is kept by downdraft                                        !
   !---------------------------------------------------------------------------------------!
   pmass_kept = 1.-pmass_left

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
         !     reference, but cap_max will be perturbed for members other than the first.  !
         !---------------------------------------------------------------------------------!
         cap_max   =  max(cap_maxs_increment,cap_maxs + capoffset * cap_maxs_increment)
         wbuoy_max = 0.
      else
         !---------------------------------------------------------------------------------!
         ! A2. If the user opted by the probability, then I use capoffset to perturb the   !
         !     threshold in which I should accept convection. wnorm_max is the maximum     !
         !     normalised value of vertical wind that will give a cumulative probability   !
         !     density function (CPDF) greater than what the user asked, assuming that the !
         !     vertical velocity with average equals to wp and standard deviation equals   !
         !     to sigw (from Mellor-Yamada closures).                                      ! 
         !---------------------------------------------------------------------------------!
         cap_max   = cap_maxs
         wbuoy_max = wnorm_max + capoffset * wnorm_increment
      end if


      !------------------------------------------------------------------------------------!
      ! B. Initialise ierr to zero, thus allowing convection to happen.  But also setting  !
      !    wbuoymin to an unrealistically large value, so it will not make sense in case   !
      !    convection fails.                                                               !
      !------------------------------------------------------------------------------------!
      ierr     = 0
      wbuoymin = 1.e20

      !------------------------------------------------------------------------------------!
      ! C. Initialise entrainment and detrainment variables.                               !
      !------------------------------------------------------------------------------------!
      mentrd_rate = .2/radius             ! Entrain. rate  associated w/ dndrafts
      mentru_rate = mentrd_rate           ! Entrain. rate  associated w/ updrafts
      cdd         = 0.                    ! Detrain. fctn. associated w/ dndrafts
      !----- Detrainment fraction associated with updrafts. -------------------------------!
      if (prec_cld) then
         cdu         = 0.1*mentru_rate
      else
         !---------------------------------------------------------------------------------!
         !     If convection is shallow, then we follow Grell et al. (1994) assumption on  !
         ! strong lateral mixing.  This avoids the already large lateral mixing to over-   !
         ! shoot and create unreasonable values of THETAE_IV and RTOT.                     !
         !---------------------------------------------------------------------------------!
         cdu         = mentru_rate
      end if

      !------------------------------------------------------------------------------------!
      ! D. Initialize mass fluxes and buoyancy, since they may never be computed in case   !
      !    convection fails.                                                               !
      !------------------------------------------------------------------------------------!
      etad_cld = 0.
      etau_cld = 0.
      dbyd     = 0.
      dbyu     = 0.

      !------------------------------------------------------------------------------------!
      ! E. Initialize Cloud work functions associated with updrafts and downdrafts.        !
      !------------------------------------------------------------------------------------!
      aau       = 0.
      aad       = 0.
      edt       = 0.
   
      !------------------------------------------------------------------------------------!
      ! F. Initialize all level-related cloud variables. They are all output variables and !
      !    some may never be assigned when convection fails. Except for kstabm, I am set-  !
      !    ting all to a non-sense value to make the point that convection did not happen. !
      !------------------------------------------------------------------------------------!
      klod    = 0
      klou    = 1
      klcl    = 0
      klfc    = 0
      kdet    = 0
      klnb    = 0
      kstabi  = 0
      kstabm  = mkx -2
      comp_dn = prec_cld

      !------------------------------------------------------------------------------------!
      ! G. Calculate all thermodynamic properties at the cloud level. The cloud levels are !
      !    staggered in relation to BRAMS model.  We must check whether the result is      !
      !    reasonable.                                                                     !
      !------------------------------------------------------------------------------------!
      call grell_thermo_cldlev(mkx,mgmzp,z_cup,exner,thil,t,qtot,qliq,qice,co2,exnersur    &
                              ,thilsur,tsur,qtotsur,qliqsur,qicesur,co2sur,exner_cup,p_cup &
                              ,t_cup,thil_cup,qtot_cup,qvap_cup,qliq_cup,qice_cup,qsat_cup &
                              ,co2_cup,rho_cup,theiv_cup,theivs_cup)
      write(which,fmt='(2(a,i4.4))') 'extrap_environment_icap=',icap,'_icld=',icld
      call grell_sanity_check(mkx,mgmzp,z_cup,p_cup,exner_cup,theiv_cup,thil_cup,t_cup     &
                             ,qtot_cup,qvap_cup,qliq_cup,qice_cup,co2_cup,rho_cup          &
                             ,which)

      !------------------------------------------------------------------------------------!
      ! H. Initialise drafts liquid mixing ratio and density.  These will be the           !
      !    convective values in case fails.                                                !
      !------------------------------------------------------------------------------------!
      qliqd_cld = qliq_cup
      qliqu_cld = qliq_cup
      qiced_cld = qice_cup
      qiceu_cld = qice_cup
      rhod_cld  = rho_cup
      rhou_cld  = rho_cup


      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! I. Find some updraft properties, namely the levels in which updrafts originate,    !
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
      ! 2. Determine a first guess for klou, the level of origin of updrafts. Check if that !
      !    didn't prevent convection to happen.                                            !
      !------------------------------------------------------------------------------------!
      call grell_updraft_origin(mkx,mgmzp,iupmethod,kpbl,kbmax,z,wwind,sigw,tke,qice,qliq  &
                               ,theiv_cup,ierr,klou)
      !------------------------------------------------------------------------------------!


      if (ierr /= 0) then
         ierr_cap(icap) = ierr
         cycle stacloop
      end if

      !------------------------------------------------------------------------------------!
      ! 3. Find the level of free convection.  Two important points here:                  !
      !    a. This call may end up preventing convection, so I must check after the call   !
      !    b. This subroutine may also affect the updraft originating level.               !
      !------------------------------------------------------------------------------------!
      call grell_find_cloud_lfc(mkx,mgmzp,kbmax,cap_max,wbuoy_max,wwind,sigw,z_cup         &
                               ,exner_cup,p_cup,theiv_cup,thil_cup,t_cup,qtot_cup,qvap_cup &
                               ,qliq_cup,qice_cup,qsat_cup,co2_cup,rho_cup,dzd_cld         &
                               ,mentru_rate,theivu_cld,thilu_cld,tu_cld,qtotu_cld          &
                               ,qvapu_cld,qliqu_cld,qiceu_cld,qsatu_cld,co2u_cld,rhou_cld  &
                               ,dbyu,klou,ierr,klcl,klfc,wbuoymin)

      if (ierr /= 0) then
         ierr_cap(icap) = ierr
         cycle stacloop
      end if

      !------------------------------------------------------------------------------------!
      ! 4. Find the minimum saturation thetae_iv. This will be the bottom of the stable    !
      !    layer.                                                                          !
      !------------------------------------------------------------------------------------!
      kstabi=(klfc - 1) + minloc(theivs_cup(klfc:kstabm),dim=1)

      !------------------------------------------------------------------------------------!
      ! 5. Increase the detrainment in stable layers provided that there is such layer.    !
      !    this rate increases linearly until a maximum value, currently set to 10 times   !
      !    the entrainment rate.  If the cloud is sufficiently small, we further simplify  !
      !    and assume detrainment to be equal to entrainment (otherwise the detrainment    !
      !    can become unrealistically large.                                               !
      !------------------------------------------------------------------------------------!
      if (prec_cld) then
         do k=kstabi,kstabm-1
            cdu(k) = min(10.*mentru_rate(k) , cdu(k-1) + 1.5 * mentru_rate(k))
         end do
      end if

      !------------------------------------------------------------------------------------!
      ! 6. Calculate the incloud ice-vapour equivalent potential temperature               !
      !------------------------------------------------------------------------------------!
      call grell_theiv_updraft(mkx,mgmzp,klou,klfc,cdu,mentru_rate,theiv,theiv_cup,dzu_cld &
                              ,theivu_cld)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! 7. Find the normalized mass fluxes associated with updrafts.  Since we are using   !
      !    height-based vertical coordinate, there is no need to find the forced           !
      !    normalized mass flux, they'll be the same, so just copy it afterwards.          !
      !------------------------------------------------------------------------------------!
      call grell_nms_updraft(mkx,mgmzp,klou,klfc,ktpse,mentru_rate,cdu,dzu_cld,etau_cld)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! 8. Find the moisture profiles associated with updrafts, and check whether the      !
      !    profile makes sense or not.                                                     !
      !------------------------------------------------------------------------------------!
      write(which,fmt='(2(a,i4.4))') 'extrap_updraft_icap=',icap,'_icld=',icld
      call grell_most_thermo_updraft(prec_cld,.true.,mkx,mgmzp,klfc,ktpse,cld2prec,cdu     &
                                    ,mentru_rate,qtot,co2,z_cup,p_cup,exner_cup            &
                                    ,theiv_cup,thil_cup,t_cup,qtot_cup,qvap_cup,qliq_cup   &
                                    ,qice_cup,qsat_cup,co2_cup,rho_cup,theivu_cld,etau_cld &
                                    ,dzu_cld,thilu_cld,tu_cld,qtotu_cld,qvapu_cld          &
                                    ,qliqu_cld,qiceu_cld,qsatu_cld,co2u_cld,rhou_cld,dbyu  &
                                    ,pwu_cld,pwav,klnb,ktop,ierr,which)
      call grell_sanity_check(mkx,mgmzp,z_cup,p_cup,exner_cup,theivu_cld,thilu_cld,tu_cld  &
                             ,qtotu_cld,qvapu_cld,qliqu_cld,qiceu_cld,co2u_cld,rhou_cld    &
                             ,which)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 9. Check whether we found a cloud top. Since this may keep convection to           !
      !    happen, I check whether I should move on or break here.  Also check whether     !
      !    this cloud qualifies to be in this spectrum size. It need to be thicker than    !
      !    the minimum value provided by the user, and there must be some condensation     !
      !    (we don't solve dry convection, that should be done by the turbulence scheme).  !
      !------------------------------------------------------------------------------------!
      if (ierr /= 0) then
         ierr_cap(icap) = ierr
         cycle stacloop
      elseif (z_cup(ktop)-z_cup(klou) < depth_min) then
         ierr_cap(icap) = 6
         cycle stacloop
      elseif (z_cup(ktop)-z_cup(klou) > depth_max) then
         ierr_cap(icap) = 7
         cycle stacloop
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !10. Find the cloud work function associated with updrafts. If this cloud doesn't    !
      !    produce cloud work, break the run, we don't simulate lazy clouds in this model. !
      !------------------------------------------------------------------------------------!
      call grell_cldwork_updraft(mkx,mgmzp,klou,ktop,dbyu,dzu_cld,etau_cld,aau)
      if (aau == 0.) then
         ierr_cap(icap) = 10
         cycle stacloop
      end if
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!




      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! J. Find the downdraft counterpart of the above properties, namely where the        !
      !    downdrafts detrain all its mass, where they originate, their mass, energy and   !
      !    moisture properties. This should be done only when it is a cloud that has       !
      !    downdrafts.  In case we cannot find a suitable level in which downdrafts        !
      !    originate, we allow clouds to happen, but with no downdrafts.  Once downdrafts  !
      !    are assumed to exist, they are integral part of clouds, so if the cloud doesn't !
      !    behave well, the entire cloud is forbidden.  Downdrafts will never be allowed   !
      !    in non-precipitating clouds.                                                    !
      !------------------------------------------------------------------------------------!
      if (comp_dn) then
         !---------------------------------------------------------------------------------!
         ! 1. The level in which downdrafts should detrain (almost) all its mass.          !
         !---------------------------------------------------------------------------------!
         kdetloop: do kdet = 1,mkx
            if (z_cup(kdet) > z_detr) exit kdetloop
         end do kdetloop

         !---------------------------------------------------------------------------------!
         ! 2. The level in which downdrafts will originate in this cloud.                  !
         !---------------------------------------------------------------------------------!
         call grell_find_downdraft_origin(mkx,mgmzp,klou,ktop,relheight_down,zcutdown      &
                                         ,z_cup,theivs_cup,dzd_cld,comp_dn,kdet,klod)
      end if

      !------------------------------------------------------------------------------------!
      !     In case we didn't find a good downdraft level, comp_dn has become .false..     !
      ! Check again.                                                                       !
      !------------------------------------------------------------------------------------!
      if (comp_dn) then
         !---------------------------------------------------------------------------------!
         ! 3. Now that we know kdet, change detrainment rate below this level.             !
         !------------------------------------------------------------------------------------!
         do k=kdet,1,-1
            cdd(k) = mentrd_rate(k)                                                        &
                   + (1. - (pmass_kept*z_cup(k) + pmass_left*z_cup(kdet))                  &
                     / (pmass_kept*z_cup(k+1) + pmass_left*z_cup(kdet)) )/dzd_cld(k)
         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         ! 4. Normalised mass fluxes.                                                      !
         !---------------------------------------------------------------------------------!
         call grell_nms_downdraft(mkx,mgmzp,kdet,klod,mentrd_rate,cdd,z_cup,dzd_cld        &
                                 ,etad_cld)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         ! 5. Moist static energy.                                                         !
         !---------------------------------------------------------------------------------!
         call grell_theiv_downdraft(mkx,mgmzp,klod,cdd,mentrd_rate,theiv,theiv_cup         &
                                   ,theivs_cup,dzd_cld,theivd_cld)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         ! 6. Moisture properties of downdraft, and its sanity check.  Besides the thermo- !
         !    dynamics, we must check whether buoyancy makes sense, in case it doesn't we  !
         !    will assign an error flag to this cloud and don't let it happen.             !
         !---------------------------------------------------------------------------------!
         write(which,fmt='(2(a,i4.4))') 'extrap_downdraft_icap=',icap,'_icld=',icld
         call grell_most_thermo_downdraft(mkx,mgmzp,klod,qtot,co2,mentrd_rate,cdd,z_cup    &
                                         ,p_cup,exner_cup,thil_cup,t_cup,qtot_cup,qvap_cup &
                                         ,qliq_cup,qice_cup,qsat_cup,co2_cup,rho_cup,pwav  &
                                         ,theivd_cld,etad_cld,dzd_cld,thild_cld,td_cld     &
                                         ,qtotd_cld,qvapd_cld,qliqd_cld,qiced_cld          &
                                         ,qsatd_cld,co2d_cld,rhod_cld,dbyd,pwd_cld,pwev    &
                                         ,ierr,which)
         call grell_sanity_check(mkx,mgmzp,z_cup,p_cup,exner_cup,theivd_cld,thild_cld      &
                                ,td_cld,qtotd_cld,qvapd_cld,qliqd_cld,qiced_cld,co2d_cld   &
                                ,rhod_cld,which)
         if (ierr /= 0) then
            ierr_cap(icap) = ierr
            cycle stacloop
         end if
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         ! 7. Compute the downdraft strength in terms of windshear.  Remembering that edt  !
         !    is also the fraction between downdraft and updraft.                          !
         !---------------------------------------------------------------------------------!
         call grell_efficiency_ensemble(mkx,mgmzp,maxens_eff,klou,klfc,klnb,edtmin,edtmax  &
                                       ,pwav,pwev,z_cup,uwind,vwind,dzd_cld                &
                                       ,edt_eff(:,icap),icld,icap)
         !------------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         ! 8. Check for water availability and evaporation consistency: we assume that     !
         !    downdraft is always saturated, and it gets the moisture from the rain.  How- !
         !    ever, it cannot require more rain than what is available, so if that would   !
         !    be the case, we don't allow this cloud to exist.                             !
         !------------------------------------------------------------------------------------!
         ddcheckloop: do iedt = 1, maxens_eff
            if (pwav + edt_eff(iedt,icap) * pwev < 0. ) then
               ierr_cap(icap) = 9
               cycle stacloop
            end if
         end do ddcheckloop
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         ! 9. Compute cloud work function associated with downdrafts.                      !
         !---------------------------------------------------------------------------------!
         call grell_cldwork_downdraft(mkx,mgmzp,klod,dbyd,dzd_cld,etad_cld,aad)
         !---------------------------------------------------------------------------------!

      else
         !---------------------------------------------------------------------------------!
         !    If the cloud doesn't have dowdrafts, we set some key variables to innocuous  !
         ! values, just to avoid them to be undefined.  kdet and klod were already         !
         ! initialised to zero.                                                            !
         !---------------------------------------------------------------------------------!
         etad_cld           = 0.
         dbyd               = 0.
         thild_cld          = thil_cup
         theivd_cld         = theiv_cup
         qtotd_cld          = qtot_cup
         qliqd_cld          = qliq_cup
         qiced_cld          = qice_cup
         co2d_cld           = co2_cup
         pwd_cld            = 0.
         pwev               = 0.
         aad                = 0.
         do iedt=1,maxens_eff
            edt_eff(iedt,icap) = 0.
         end do
      end if
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!


      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! K. Find the non-forced cloud work, which will require computing most subroutines   !
      !    again.  This will go and compute even if the tests would prevent the cloud to   !
      !    happen. This part will be skipped if the user is asking for moisture            !
      !    convergence only.                                                               !
      !------------------------------------------------------------------------------------!
      if (comp_noforc_cldwork) then
         !---------------------------------------------------------------------------------!
         ! i.     Initialise error flag.                                                   !
         !---------------------------------------------------------------------------------!
         ierr0    = 0
         comp_dn0 = comp_dn
      
         !---------------------------------------------------------------------------------!
         ! ii.    Initialise detrainment and cloud work variables.                         !
         !---------------------------------------------------------------------------------!
         aa0d      = 0.
         aa0u      = 0.

         !---------------------------------------------------------------------------------!
         ! iii.   Calculate all thermodynamic properties at the cloud level and initialize !
         !        draft thermodynamic properties (sanity check too...).                    !
         !---------------------------------------------------------------------------------!
         call grell_thermo_cldlev(mkx,mgmzp,z_cup,exner0,thil0,t0,qtot0,qliq0,qice0,co20   &
                                 ,exnersur,thilsur,tsur,qtotsur,qliqsur,qicesur,co2sur     &
                                 ,exner0_cup,p0_cup,t0_cup,thil0_cup,qtot0_cup,qvap0_cup   &
                                 ,qliq0_cup,qice0_cup,qsat0_cup,co20_cup,rho0_cup          &
                                 ,theiv0_cup,theivs0_cup)
         write(which,fmt='(2(a,i4.4))') 'zero_environment_icap=',icap,'_icld=',icld
         call grell_sanity_check(mkx,mgmzp,z_cup,p0_cup,exner0_cup,theiv0_cup,thil0_cup    &
                                ,t0_cup,qtot0_cup,qvap0_cup,qliq0_cup,qice0_cup,co20_cup   &
                                ,rho0_cup,which)

         !---------------------------------------------------------------------------------!
         ! iv.    Calculate the thermodynamic properties below the level of free convec-   !
         !        tion.                                                                    !
         !---------------------------------------------------------------------------------!
         write(which,fmt='(2(a,i4.4))') 'buoy_below_lfc_icap=',icap,'_icld=',icld
         call grell_buoy_below_lfc(mkx,mgmzp,klou,klfc,z_cup,exner0_cup,p0_cup,theiv0_cup  &
                                  ,thil0_cup,t0_cup,qtot0_cup,qvap0_cup,qliq0_cup          &
                                  ,qice0_cup,qsat0_cup,co20_cup,rho0_cup,theiv0u_cld       &
                                  ,thil0u_cld,t0u_cld,qtot0u_cld,qvap0u_cld,qliq0u_cld     &
                                  ,qice0u_cld,qsat0u_cld,co20u_cld,rho0u_cld,dby0u,which)

         !---------------------------------------------------------------------------------!
         ! v.     Calculate the incloud ice-vapour equivalent potential temperature        !
         !---------------------------------------------------------------------------------!
         call grell_theiv_updraft(mkx,mgmzp,klou,klfc,cdu,mentru_rate,theiv0,theiv0_cup    &
                                 ,dzu_cld,theiv0u_cld)

         !---------------------------------------------------------------------------------!
         ! vi.    Find the moisture profiles associated with updrafts, and check whether   !
         !        they make sense or not.                                                  !
         !---------------------------------------------------------------------------------!
         write(which,fmt='(2(a,i4.4))') 'zero_updraft_icap=',icap
         call grell_most_thermo_updraft(comp_down_cap(icap),.false.,mkx,mgmzp,klfc,ktop    &
                                       ,cld2prec,cdu,mentru_rate,qtot0,co20,z_cup,p0_cup   &
                                       ,exner0_cup,theiv0_cup,thil0_cup,t0_cup,qtot0_cup   &
                                       ,qvap0_cup,qliq0_cup,qice0_cup,qsat0_cup,co20_cup   &
                                       ,rho0_cup,theiv0u_cld,etau_cld,dzu_cld,thil0u_cld   &
                                       ,t0u_cld,qtot0u_cld,qvap0u_cld,qliq0u_cld           &
                                       ,qice0u_cld,qsat0u_cld,co20u_cld,rho0u_cld,dby0u    &
                                       ,pw0u_cld,pwav0,klnb0,ktop0,ierr0,which)
         call grell_sanity_check(mkx,mgmzp,z_cup,p0_cup,exner0_cup,theiv0u_cld,thil0u_cld  &
                                ,t0u_cld,qtot0u_cld,qvap0u_cld,qliq0u_cld,qice0u_cld       &
                                ,co20u_cld,rho0u_cld,which)

         !---------------------------------------------------------------------------------!
         ! vii.   Find the cloud work function.                                            !
         !---------------------------------------------------------------------------------!
         call grell_cldwork_updraft(mkx,mgmzp,klou,ktop,dby0u,dzu_cld,etau_cld,aa0u)
         
         !----- Downdraft properties ------------------------------------------------------!
         if (comp_dn) then

            !------------------------------------------------------------------------------!
            ! viii.  Moist static energy.                                                  !
            !------------------------------------------------------------------------------!
            call grell_theiv_downdraft(mkx,mgmzp,klod,cdd,mentrd_rate,theiv0,theiv0_cup    &
                                      ,theivs0_cup,dzd_cld,theiv0d_cld)
            !------------------------------------------------------------------------------!
            ! ix.    Moisture properties of downdraft                                      !
            !------------------------------------------------------------------------------!
            write(which,fmt='(a,i4.4)') 'zero_downdraft_icap=',icap
            call grell_most_thermo_downdraft(mkx,mgmzp,klod,qtot0,co20,mentrd_rate,cdd     &
                                            ,z_cup,p0_cup,exner0_cup,thil0_cup,t0_cup      &
                                            ,qtot0_cup,qvap0_cup,qliq0_cup,qice0_cup       &
                                            ,qsat0_cup,co20_cup,rho0_cup,pwav0,theiv0d_cld &
                                            ,etad_cld,dzd_cld,thil0d_cld,t0d_cld           &
                                            ,qtot0d_cld,qvap0d_cld,qliq0d_cld,qice0d_cld   &
                                            ,qsat0d_cld,co20d_cld,rho0d_cld,dby0d,pw0d_cld &
                                            ,pwev0,ierr0,which)
            call grell_sanity_check(mkx,mgmzp,z_cup,p0_cup,exner0_cup,theiv0d_cld          &
                                   ,thil0d_cld,t0d_cld,qtot0d_cld,qvap0d_cld,qliq0d_cld    &
                                   ,qice0d_cld,co20d_cld,rho0d_cld,which)

            !------------------------------------------------------------------------------!
            ! x.     Downdraft cloud work function.                                        !
            !------------------------------------------------------------------------------!
            call grell_cldwork_downdraft(mkx,mgmzp,klod,dby0d,dzd_cld,etad_cld,aa0d)
         end if

      end if
      !------------------------------------------------------------------------------------!
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!



      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! L. Now we perform the loop across the precipitation efficiency ensembles to find   !
      !    the normalised tendency of conserved variables due to convection for each       !
      !    ensemble member.                                                                !
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
         if (comp_noforc_cldwork) aatot0_eff(iedt,icap) = aa0u + edt * aa0d
         aatot_eff(iedt,icap)  = aau + edt * aad


         !---------------------------------------------------------------------------------!
         ! 3. Compute the changes in thetae_iv,theta_il, and total mixing ratio at the     !
         !    bottom, considering downdraft effects. If this cloud doesn't have down-      !
         !    drafts, it should be zero so skip it.                                        !
         !---------------------------------------------------------------------------------!
         if (comp_dn) then
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,theiv,p_cup          &
                                         ,theiv_cup,mentrd_rate,cdd,dzd_cld,etad_cld       &
                                         ,theivd_cld,dellatheiv_eff(1,iedt,icap))
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,qtot,p_cup,qtot_cup  &
                                         ,mentrd_rate,cdd,dzd_cld,etad_cld,qtotd_cld       &
                                         ,dellaqtot_eff(1,iedt,icap))
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,thil,p_cup,thil_cup  &
                                         ,mentrd_rate,cdd,dzd_cld,etad_cld,thild_cld       &
                                         ,dellathil_eff(1,iedt,icap))
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,co2,p_cup,co2_cup    &
                                         ,mentrd_rate,cdd,dzd_cld,etad_cld,co2d_cld        &
                                         ,dellaco2_eff(1,iedt,icap))
         end if
         
         !---------------------------------------------------------------------------------!
         ! 4. Compute the changes in thetae_iv and total mixing ratio at other levels. In  !
         !    case downdrafts are absent, this will work fine too, since the mass fluxes   !
         !    will be all zero.                                                            !
         !---------------------------------------------------------------------------------!
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,klou,klfc,klod,ktop   &
                                   ,theiv,p_cup,theiv_cup,mentrd_rate,mentru_rate,cdd      &
                                   ,cdu,dzd_cld,etad_cld,etau_cld,theivd_cld,theivu_cld    &
                                   ,dellatheiv_eff(:,iedt,icap))
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,klou,klfc,klod,ktop   &
                                   ,qtot,p_cup,qtot_cup,mentrd_rate,mentru_rate,cdd,cdu    &
                                   ,dzd_cld,etad_cld,etau_cld,qtotd_cld,qtotu_cld          &
                                   ,dellaqtot_eff(:,iedt,icap))
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,klou,klfc,klod,ktop   &
                                   ,thil,p_cup,thil_cup,mentrd_rate,mentru_rate,cdd,cdu    &
                                   ,dzd_cld,etad_cld,etau_cld,thild_cld,thilu_cld          &
                                   ,dellathil_eff(:,iedt,icap))
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,klou,klfc,klod,ktop   &
                                   ,co2,p_cup,co2_cup,mentrd_rate,mentru_rate,cdd,cdu      &
                                   ,dzd_cld,etad_cld,etau_cld,co2d_cld,co2u_cld            &
                                   ,dellaco2_eff(:,iedt,icap))

         !---------------------------------------------------------------------------------!
         ! 5. Computing the rainfall flux, which will be simply all the fall-out water     !
         !    that wasn't consumed by the downdrafts.                                      !
         !---------------------------------------------------------------------------------!
         do k=1,ktop
            !------ Precipitation ---------------------------------------------------------!
            pw_eff(k,iedt,icap) =  pwu_cld(k) + edt * pwd_cld(k)
         end do


         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         if (print_debug) then
            iun=icld+70
            write(unit=iun,fmt='(504a)')             ('=',l=1,504)
            write(unit=iun,fmt='(504a)')             ('=',l=1,504)
            write(unit=iun,fmt='((a,1x,i5,1x))'    ) ' I      =',i
            write(unit=iun,fmt='((a,1x,i5,1x))'    ) ' J      =',j
            write(unit=iun,fmt='((a,1x,i5,1x))'    ) ' IEDT   =',iedt
            write(unit=iun,fmt='((a,1x,i5,1x))'    ) ' ICAP   =',icap
            write(unit=iun,fmt='(a,1x,f10.4,1x,i5)') ' KLOU   =',z(klou),klou
            write(unit=iun,fmt='(a,1x,f10.4,1x,i5)') ' KLCL   =',z(klcl),klcl
            write(unit=iun,fmt='(a,1x,f10.4,1x,i5)') ' KLFC   =',z(klfc),klfc
            write(unit=iun,fmt='(a,1x,f10.4,1x,i5)') ' KDET   =',z(kdet),kdet
            write(unit=iun,fmt='(a,1x,f10.4,1x,i5)') ' KLOD   =',z(klod),klod
            write(unit=iun,fmt='(a,1x,f10.4,1x,i5)') ' KLNB   =',z(klnb),klnb
            write(unit=iun,fmt='(a,1x,f10.4,1x,i5)') ' KTOP   =',z(ktop),ktop
            write(unit=iun,fmt='(a,1x,f10.4,1x,i5)') ' PWAV   =',pwav
            write(unit=iun,fmt='(a,1x,f10.4,1x,i5)') ' PWEV   =',pwev
            write(unit=iun,fmt='(a,1x,f10.4,1x,i5)') ' EDT    =',edt
            write(unit=iun,fmt='(504a)')          ('-',l=1,504)
            write(unit=iun,fmt='(a,3(1x,a))'    ) ' Type'    ,'    AA_DOWN','      AA_UP'  &
                                                             ,'     AA_TOT'
            write(unit=iun,fmt='(a,3(1x,f11.4))') ' Zero',aa0d,aa0u,aatot0_eff(iedt,icap)
            write(unit=iun,fmt='(a,3(1x,f11.4))') ' One ',aad ,aau ,aatot_eff (iedt,icap)
            write(unit=iun,fmt='(504a)')          ('-',l=1,504)
            write(unit=iun,fmt='(42(1x,a))')    '      LEVEL','      Z_CUP','      P_CUP'  &
                                 ,'       THIL','   THIL_CUP','  THILD_CLD','  THILU_CLD'  &
                                 ,'      THEIV','  THEIV_CUP',' THEIVD_CLD',' THEIVU_CLD'  &
                                 ,'          T','      T_CUP','     TD_CLD','     TU_CLD'  &
                                 ,'       QTOT','   QTOT_CUP','  QTOTD_CLD','  QTOTU_CLD'  &
                                 ,'        CO2','    CO2_CUP','   CO2D_CLD','   CO2U_CLD'  &
                                 ,'       QVAP','   QVAP_CUP','  QVAPD_CLD','  QVAPU_CLD'  &
                                 ,'       QCON','   QCON_CUP','  QCOND_CLD','  QCONU_CLD'  &
                                 ,'        RHO','    RHO_CUP','   RHOD_CLD','   RHOU_CLD'  &
                                 ,'    PWD_CLD','    PWU_CLD','     PW_TOT'                &
                                 ,' DELLA_THIL','DELLA_THEIV',' DELLA_QTOT','  DELLA_CO2'
            write(unit=iun,fmt='(504a)')          ('-',l=1,504)
            do k=mkx,1,-1
               if (k == ktop) then
                  levname = '        TOP'
               elseif (k == klnb) then
                  levname = '        LNB'
               elseif (k == klod) then
                  levname = '        LOD'
               elseif (k == kdet) then
                  levname = '        DET'
               elseif (k == klfc) then
                  levname = '        LFC'
               elseif (k == klou) then
                  levname = '        LOU'
               elseif (k == klcl) then
                  levname = '        LCL'
               else
                  levname = '           '
               end if
            
               write(unit=iun,fmt='(1x,a,41(1x,f11.4))') levname,z_cup(k),0.01*p_cup(k)    &
                    ,      thil(k),      thil_cup(k),      thild_cld(k),      thilu_cld(k) &
                    ,     theiv(k),     theiv_cup(k),     theivd_cld(k),     theivu_cld(k) &
                    ,     t(k)-t00,     t_cup(k)-t00,     td_cld(k)-t00,     tu_cld(k)-t00 &
                    ,1000.*qtot(k),1000.*qtot_cup(k),1000.*qtotd_cld(k),1000.*qtotu_cld(k) &
                    ,       co2(k),       co2_cup(k),       co2d_cld(k),       co2u_cld(k) &
                    ,1000.*qvap(k),1000.*qvap_cup(k),1000.*qvapd_cld(k),1000.*qvapu_cld(k) &
                    ,1000.*(qliq(k)+qice(k)),1000.*(qliq_cup(k)+qice_cup(k))               &
                    ,1000.*(qliqd_cld(k)+qiced_cld(k)),1000.*(qliqu_cld(k)+qiceu_cld(k))   &
                    ,       rho(k),       rho_cup(k),       rhod_cld(k),       rhou_cld(k) &
                    ,1000.*pwd_cld(k),1000.*pwu_cld(k),1000.*pw_eff(k,iedt,icap)           &
                    ,dellathil_eff(k,iedt,icap)     ,dellatheiv_eff(k,iedt,icap)           &
                    ,1000*dellaqtot_eff(k,iedt,icap),dellaco2_eff(k,iedt,icap)
            end do
            write(unit=iun,fmt='(504a)')          ('-',l=1,504)
            write(unit=iun,fmt='(42(1x,a))')    '      LEVEL','      Z_CUP','      P_CUP'  &
                                 ,'       THIL','   THIL_CUP','  THILD_CLD','  THILU_CLD'  &
                                 ,'      THEIV','  THEIV_CUP',' THEIVD_CLD',' THEIVU_CLD'  &
                                 ,'          T','      T_CUP','     TD_CLD','     TU_CLD'  &
                                 ,'       QTOT','   QTOT_CUP','  QTOTD_CLD','  QTOTU_CLD'  &
                                 ,'        CO2','    CO2_CUP','   CO2D_CLD','   CO2U_CLD'  &
                                 ,'       QVAP','   QVAP_CUP','  QVAPD_CLD','  QVAPU_CLD'  &
                                 ,'       QCON','   QCON_CUP','  QCOND_CLD','  QCONU_CLD'  &
                                 ,'        RHO','    RHO_CUP','   RHOD_CLD','   RHOU_CLD'  &
                                 ,'    PWD_CLD','    PWU_CLD','     PW_TOT'                &
                                 ,' DELLA_THIL','DELLA_THEIV',' DELLA_QTOT','  DELLA_CO2'
            write(unit=iun,fmt='(504a)')          ('=',l=1,504)
            write(unit=iun,fmt='(504a)')          ('=',l=1,504)
            write(unit=iun,fmt='(a)') ' '
            write(unit=iun,fmt='(a)') ' '
            write(unit=iun,fmt='(a)') ' '
            write(unit=iun,fmt='(a)') ' '
            write(unit=iun,fmt='(a)') ' '
            write(unit=iun,fmt='(a)') ' '
         end if
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      end do effloop
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!

      !------------------------------------------------------------------------------------!
      ! M.  Before we cycle to the next static control, we copy the variables that depend  !
      !     on the static control.  In case this static control didn't happen and we don't !
      !     reach this point, no problem, the values will keep the original value (zero).  !
      !------------------------------------------------------------------------------------!
      ierr_cap     (icap)     = 0
      comp_down_cap(icap)     = comp_dn
      klod_cap(icap)          = klod
      klou_cap(icap)          = klou
      klcl_cap(icap)          = klcl
      klfc_cap(icap)          = klfc
      kdet_cap(icap)          = kdet
      kstabi_cap(icap)        = kstabi
      kstabm_cap(icap)        = kstabm
      klnb_cap(icap)          = klnb
      ktop_cap(icap)          = ktop

      pwav_cap(icap)          = pwav
      pwev_cap(icap)          = pwev
      wbuoymin_cap(icap)      = wbuoymin
      
      do k=1,mkx
         cdd_cap(k,icap)         = cdd(k)
         cdu_cap(k,icap)         = cdu(k)
         mentrd_rate_cap(k,icap) = mentrd_rate(k)
         mentru_rate_cap(k,icap) = mentru_rate(k)
         dbyd_cap(k,icap)        = dbyd(k)
         dbyu_cap(k,icap)        = dbyu(k)
         etad_cld_cap(k,icap)    = etad_cld(k)
         etau_cld_cap(k,icap)    = etau_cld(k)
         rhod_cld_cap(k,icap)    = rhod_cld(k)
         rhou_cld_cap(k,icap)    = rhou_cld(k)
         qliqd_cld_cap(k,icap)   = qliqd_cld(k)
         qliqu_cld_cap(k,icap)   = qliqu_cld(k)
         qiced_cld_cap(k,icap)   = qiced_cld(k)
         qiceu_cld_cap(k,icap)   = qiceu_cld(k)
      end do

   end do stacloop

   return
end subroutine grell_cupar_static
!==========================================================================================!
!==========================================================================================!
